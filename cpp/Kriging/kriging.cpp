#include "kriging.h"
#include <math.h>
//#include "..\HOPS\hops.cpp"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "src/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include "src/interpolation.h"
#include <nlopt.h>


using namespace Eigen;



//gaussian model, not used i think
static double gaussModel(double h, double r, double c0, double b) {
	double y;
	double a = r / 2.0;
	y = b + c0 * (1 - exp(-(h * h) / (a * a)));
	return y;
}

//from adjoint method, calculates qoi (scattered field)
std::complex<double> Kriging::get_QoI2(Domain& dom_forward, Domain& dom_adjoint) {
	std::complex<double> qoi = 0.0;
	for (int i = 0; i < dom_forward.cAlpha.size(); ++i) {
		qoi += dom_forward.cAlpha[i] * std::conj(dom_adjoint.scatter1.cGr[1][i]);
	}
	return qoi;
}

// complex double divided by double
static inline std::complex<double> operator/(std::complex<double> c1, double d1) {
	std::complex<double> output = { c1.real() / d1, c1.imag() / d1 };
	return output;
}


//currently using for real part ie material parameters
//find_lowest takes an array and a value and returns the lowest index greater than the value
static int find_lowest(std::vector<std::complex<double>> array, std::complex<double> value) {
	double min = 1e15;
	int index;
	for (int i = 0; i < array.size(); i++) {
		if (value.real() - array[i].real() < 0) {
			//min = abs(array[i].real() - value);
			index = i;
			break;
		}
	}
	return index;

}

//as opposed to find_lowest, find_index returns the closest index value in the array
static int find_index(std::vector<std::complex<double>> array, std::complex<double> value) {
	double min = 1e15;
	int index;
	for (int i = 0; i < array.size(); i++) {
		if (abs(value.real() - array[i].real()) < min) {
			min = abs(array[i].real() - value.real());
			index = i;
		}
	}
	return index;

}

//complex comparitor for if (complex double a < complex double b)
//currently only returns real part since this is the uncertain parameter
static bool complexComparitor(std::complex<double> a, std::complex<double> b) {
	return real(a) <= real(b);
}

//insert a value into a complex array
static void insert_complex(std::vector<std::complex<double>>& array, int index, std::complex<double> value) {
	std::vector<std::complex<double>> tempArray(array.size() + 1);

	int count = 0;
	for (int i = 0; i < tempArray.size(); ++i) {
		if (i == index) {
			tempArray[i] = value;
		}
		else {
			tempArray[i] = array[count];
			count++;
		}
	}
	//array.clear();
	array.resize(tempArray.size());
	array = tempArray;

	//for (int i = 0; i < array.size(); ++i) {
	//	array[i] = tempArray[i];
	//	std::cout << array[i] << std::endl;
	//}

}

//L2 norm, used for material_list split (only needed for HOPS, not Kriging)
static double norm(std::vector < std::complex<double>> v1) {
	std::complex<double> norm = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		norm += v1[i] * std::conj(v1[i]);
	}
	double normd = abs(norm);
	return sqrt(normd);

}

//complex subtraction operator
static inline std::vector<std::complex<double>> operator-(std::vector < std::complex<double>>& v1, std::vector < std::complex<double>>& v2) {
	std::vector<std::complex<double>> vout;
	for (int i = 0; i < v1.size(); ++i) {
		vout.push_back(v1[i] - v2[i]);
	}
	return vout;
}


//matrix is filled according to pairs of the reference points.  M(i,j) = {X(i), X(j), distance between X(i) X(j), ect}
static std::vector<std::vector<std::complex<double>>> buildDMatrix2(std::vector<std::complex<double>> x, std::vector<std::complex<double>> yI, std::vector<std::complex<double>> yJ) {
	std::vector<std::vector<std::complex<double>>> dMat(x.size() * x.size());
	int index = 0;
	for (int i = 0; i < x.size(); ++i) {
		double RCSi = sqrt(yI[i].real() * yI[i].real() + yI[i].imag() * yI[i].imag()) / 4 / 3.14159;
		for (int j = 0; j < x.size(); ++j) {
			double RCSj = sqrt(yJ[j].real() * yJ[j].real() + yJ[j].imag() * yJ[j].imag()) / 4 / 3.14159;

			//std::cout << "index: " << index << " x.size(): " << x.size() << std::endl;
			dMat[index].push_back(double(i));
			dMat[index].push_back(double(j));
			dMat[index].push_back(abs(x[i].real() - x[j].real()));

			//3rd term is the variogram term
			//will square the term when building the variogram instead of here
			//dMat[index].push_back(abs(RCSi - RCSj));// *(RCSi - RCSj));
			dMat[index].push_back(abs(yI[i] - yJ[j]));

			dMat[index].push_back(yI[i]);
			dMat[index].push_back(yI[j]);
			dMat[index].push_back(yI[i].real() - yJ[j].real());
			dMat[index].push_back(yI[i].imag() - yJ[j].imag());

			index++;
		}
	}
	return dMat;
}

//derivative calculation using finite differences
static std::vector<double> firstDerivative(std::vector<double> function, std::vector<double> distance) {
	std::vector<double> derivative;
	for (int i = 0; i < function.size(); ++i) {

		if (i <= 2) {
			derivative.push_back((function[i + 1] - function[i])  / (distance[i + 1] - distance[i]));
		}
		else if (i == function.size() - 1) {
			derivative.push_back((function[i - 1] - function[i]) / (distance[i - 1] - distance[i]));
		}
		else {
			derivative.push_back((function[i + 1] - function[i - 1]) / (distance[i + 1] - distance[i - 1]));
		}
	}
	return derivative;
}


//builds the empirical variogram
static std::vector<double> buildVariogram(std::vector<std::vector<std::complex<double>>> dMat, std::vector<std::complex<double>> x, std::vector<double>& edges) {

	std::vector<double> binCounts(edges.size()-1);
	std::vector<double> binIdx(x.size()*x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	std::vector<double> variogram(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					 //variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					variogram[int(binIdx[j])] += real(dMat[j][3] * dMat[j][3]);
				}
			}
		}
	}
	for (int i = 0; i < variogram.size(); ++i) {

		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		variogram[i] = variogram[i] / (2.0 * binCounts[i]);
	}
	return variogram;
	
}

//builds a complex empirical variogram
static std::vector<std::complex<double>> buildVariogramComplex(std::vector<std::vector<std::complex<double>>> dMat, std::vector<std::complex<double>> x, std::vector<double>& edges) {

	std::vector<double> binCounts(edges.size() - 1);
	std::vector<double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	std::vector<std::complex<double>> variogram(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					variogram[int(binIdx[j])] += { abs(dMat[j][6] * dMat[j][6]), abs(dMat[j][7] * dMat[j][7]) };
				}
			}
		}
	}
	for (int i = 0; i < variogram.size(); ++i) {

		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		variogram[i] = variogram[i] / (2.0 * binCounts[i]);
	}
	return variogram;

}

//builds empirical covariance function
static std::vector<double> buildCoVariogram(std::vector<std::vector<std::complex<double>>> dMat, std::vector<std::complex<double>> x, std::vector<double>& edges) {

	std::complex<double> sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += dMat[i][4];
	}
	std::complex<double> avg = sum / x.size();

	std::vector<double> binCounts(edges.size() - 1);
	std::vector<double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	//calculate average
	std::vector<double> average(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					average[i] += real(dMat[j][4] + dMat[j][5]);
				}
			}
		}
	}
	for (int i = 0; i < average.size(); ++i) {
		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		average[i] = average[i] / (2.0 * binCounts[i]);
	}
	
	//calculate covariogram
	std::vector<double> covariogram(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					//covariogram[i] += abs((dMat[j][4] - average[i]) * (dMat[j][5] - average[i]));
					covariogram[int(binIdx[j])] += abs((dMat[j][4] - avg) * (dMat[j][5] - avg));
				}
			}
		}
	}
	for (int i = 0; i < covariogram.size(); ++i) {

		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		covariogram[i] = covariogram[i] / (binCounts[i]);
		//covariogram[i] = 0.0;
	}
	return covariogram;

}

//builds a complex empirical covariance function
static std::vector<std::complex<double>> buildCoVariogramComplex(std::vector<std::vector<std::complex<double>>> dMat, std::vector<std::complex<double>> x, std::vector<double>& edges) {

	std::complex<double> sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += dMat[i][4];
	}
	std::complex<double> avg = sum / x.size();

	std::vector<double> binCounts(edges.size() - 1);
	std::vector<double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	//calculate average
	std::vector<double> averageReal(edges.size() - 1);
	std::vector<double> averageImag(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					averageReal[int(binIdx[j])] += (real(dMat[j][4]) + real(dMat[j][5]));
					averageImag[int(binIdx[j])] += (imag(dMat[j][4]) + imag(dMat[j][5]));
				}
			}
		}
	}

	for (int i = 0; i < averageReal.size(); ++i) {
		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		averageReal[i] = averageReal[i] / (2.0 * binCounts[i]);
		averageImag[i] = averageImag[i] / (2.0 * binCounts[i]);
	}

	std::vector<std::complex<double>> average(averageReal.size());
	for (int i = 0; i < averageReal.size(); ++i) {
		average[i] = { averageReal[i], averageImag[i] };
	}

	//calculate covariogram
	std::vector<std::complex<double>> covariogram(edges.size() - 1);
	std::vector<double> covariogramReal(edges.size() - 1);
	std::vector<double> covariogramImag(edges.size() - 1);

	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					//covariogram[int(binIdx[j])] += abs((dMat[j][4] - average[i]) * (dMat[j][5] - average[i]));
					//covariogram[int(binIdx[j])] += ((dMat[j][4] - average[i]) * (dMat[j][5] - average[i]));
					covariogramReal[int(binIdx[j])] += (dMat[j][4].real() - averageReal[i]) * (dMat[j][5].real() - averageReal[i]);
					covariogramImag[int(binIdx[j])] += (dMat[j][4].imag() - averageImag[i]) * (dMat[j][5].imag() - averageImag[i]);

				}
			}
		}
	}
	for (int i = 0; i < covariogram.size(); ++i) {
		covariogram[i] = { covariogramReal[i], covariogramImag[i] };
		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		covariogram[i] = covariogram[i] / (binCounts[i]);
	}
	return covariogram;

}

//correlation function attempt
static std::vector<double> buildCorrelation(std::vector<std::vector<std::complex<double>>> dMat, std::vector<std::complex<double>> x, std::vector<double>& edges, std::vector<std::complex<double>> referenceVals, std::vector<double> d) {

	std::complex<double> sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += dMat[i][4];
	}
	std::complex<double> avg = sum / x.size();

	std::vector<double> binCounts(edges.size() - 1);
	std::vector<double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	//calculate average
	std::vector<double> average(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					average[int(binIdx[j])] += real(dMat[j][4] + dMat[j][5]);
				}
			}
		}
	}
	for (int i = 0; i < average.size(); ++i) {
		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		average[i] = average[i] / (2.0 * binCounts[i]);
	}

	//calculate correlation
	std::vector<double> correlation(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					//covariogram[int(binIdx[j])] += abs((dMat[j][4] - average[i]) * (dMat[j][5] - average[i]));
					correlation[int(binIdx[j])] += abs((dMat[j][4] - avg) * (dMat[j][5] - avg));
				}
			}
		}
	}
	for (int i = 0; i < correlation.size(); ++i) {

		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		correlation[i] = correlation[i] / (binCounts[i]);
	}
	return correlation;

}

//old variance calculation
static double calcVariance2(std::vector<double> variogramFit, std::vector<double> variogramGradFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	double variance = 0.0;
	const int size = xSample.size();
	MatrixXd varMat(size, size);
	MatrixXd varMatGrad(size, size);

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			varMat(i, j) = variogramFit[index];
			varMatGrad(i, j) = variogramGradFit[index];
		}
	}
	MatrixXd weightVec(size, 1);
	MatrixXd weightVecGrad(size, 1);
	MatrixXd covVec(size, 1);
	for (int i = 0; i < size; ++i) {
		double xDist = abs(xSample[i].real() - xVal);
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		covVec(i, 0) = variogramFit[index];
		weightVec(i, 0) = weights[i];
	}

	std::cout << "weights size: " << weights.size() << "size: " << size << std::endl;

	for (int i = size; i < size * 2.0; ++i) {
		weightVecGrad(i - size, 0) = weights[i];
	}


	//std::cout << "weightVec\n" << weightVec << std::endl;
	//std::cout << "varMat\n" << varMat << std::endl;
	//std::cout << "covVec\n" << covVec << std::endl;

	//MatrixXd varianceVec = weightVec.transpose() * varMat * weightVec + weightVecGrad.transpose() * varMatGrad * weightVecGrad - covVec.transpose() * weightVec - weightVec.transpose() * covVec;
	MatrixXd varianceVec = weightVec.transpose() * varMat * weightVec - covVec.transpose() * weightVec - weightVec.transpose() * covVec;
	//MatrixXd varianceVec = 


	//std::cout << "old variance: " << varianceVec;
	variance = abs(varianceVec(0, 0));

	return variance;
}

//old variance calculation
static double calcVariance(std::vector<double> variogramFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	double variance = 0.0;
	const int size = xSample.size();
	MatrixXd varMat(size, size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			varMat(i, j) = variogramFit[index];
		}
	}
	MatrixXd weightVec(size, 1);
	MatrixXd covVec(size, 1);
	for (int i = 0; i < size; ++i) {
		double xDist = abs(xSample[i].real() - xVal);
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		covVec(i, 0) = variogramFit[index];
		weightVec(i, 0) = weights[i];
	}

	//std::cout << "weightVec\n" << weightVec << std::endl;
	//std::cout << "varMat\n" << varMat << std::endl;
	//std::cout << "covVec\n" << covVec << std::endl;

	MatrixXd varianceVec = weightVec.transpose() * varMat * weightVec - covVec.transpose() * weightVec - weightVec.transpose() * covVec;
	//std::cout << "old variance: " << varianceVec;
	variance = abs(varianceVec(0, 0));

	return variance;
}

//Kriging Matrix Solution
static void kriging(std::vector<double> variogramFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	const int size = xSample.size();
	//Matrix <double, Dynamic, Dynamic> krigMat;
	MatrixXd krigMat(size+1, size+1);
	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;
	for (int i = 0; i < size + 1; ++i) {
		for (int j = 0; j < size + 1; ++j) {


			if (i == size) {
				krigMat(i, j) = 1.0;
			}
			else if (j == size) {
				krigMat(i, j) = 1.0;
			}
			else {
				double dist = abs(xSample[i].real() - xSample[j].real());
				auto it = std::lower_bound(d.begin(), d.end(), dist);
				int index = it - d.begin();
				krigMat(i, j) = variogramFit[index];
			}
			//std::cout << " mat(i, j): " << krigMat(i, j) << std::endl;
		}
	}
	krigMat(size, size) =  0.0;
	
	//std::cout << "krigMat: " << krigMat << std::endl;

	MatrixXd krigVec(size + 1, 1);

	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		krigVec(i, 0) = variogramFit[index];
	}

	krigVec(size, 0) = 1.0;
	//std::cout << " krigVec: " << krigVec << std::endl;


	MatrixXd weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	for (int i = 0; i < size; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i,0));
	}
	//std::cout << weightMat << std::endl;
	//std::cout << std::endl;

}

double fact(int n) {
	double factorial;
	for (int i = 1; i <= n; ++i) {
		factorial *= double(i);
	}
	return factorial;
}

//attempt at Kriging matrix solution but with Bezier basis functions
static void bezierKriging(std::vector<double> variogramFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();

	//if (M > 18) { M = 18; }

	//Matrix <double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd fMat(size, M);
	MatrixXd zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
		}
	}

	int n = M;
	// F Matrix


	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMat(i, j) = 0.0;
		}
	}

	MatrixXd upper(size, size + M);
	MatrixXd lower(M, size + M);
	MatrixXd krigMat(size + M, size + M);

	upper << sMat, fMat;
	lower << fMat.transpose(), zeroMat;

	krigMat << upper, lower;

	//std::cout << "krigMat: " << krigMat << std::endl;

	MatrixXd cVec(size, 1);
	MatrixXd fVec(M, 1);
	MatrixXd krigVec(size + M, 1);

	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal), i - M/2);
	}

	krigVec << cVec, fVec;

	MatrixXd weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	//std::cout << krigMat << std::endl;
	//std::cout << krigMat.inverse() << std::endl;
	//std::cout << weightMat << std::endl;

	for (int i = 0; i < size; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0));
	}
	//std::cout << weightMat << std::endl;
	//std::cout << std::endl;

}

//Kriging Matrix solution using Fourier Basis functions
static void fourierKriging(std::vector<double> variogramFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();
	//int M = 25;

	if (M > 4) { M = 4; }

	//Matrix <double, Dynamic, Dynamic> krigMat;
	Matrix <std::complex<double>, Dynamic, Dynamic> sMat(size, size);
	Matrix <std::complex<double>, Dynamic, Dynamic> fMat(size, M);
	Matrix <std::complex<double>, Dynamic, Dynamic> zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
		}
	}


	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			std::complex<double> num(cos(2 * 3.14159 * j * (xSample[i].real() - average) / average), sin(2 * 3.14159 * j * (xSample[i].real() - average) / average));
			fMat(i, j) = num;
		}
	}

	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMat(i, j) = 0.0;
		}
	}

	Matrix <std::complex<double>, Dynamic, Dynamic> upper(size, size + M);
	Matrix <std::complex<double>, Dynamic, Dynamic> lower(M, size + M);
	Matrix <std::complex<double>, Dynamic, Dynamic> krigMat(size + M, size + M);

	upper << sMat, fMat;
	lower << fMat.transpose(), zeroMat;

	krigMat << upper, lower;

	//std::cout << "krigMat: " << krigMat << std::endl;

	Matrix <std::complex<double>, Dynamic, Dynamic> cVec(size, 1);
	Matrix <std::complex<double>, Dynamic, Dynamic> fVec(M, 1);
	Matrix <std::complex<double>, Dynamic, Dynamic> krigVec(size + M, 1);

	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
	}

	for (int i = 0; i < M; ++i) {

		std::complex<double> num(cos(2 * 3.14159 * i * (xVal - average) / average), sin(2 * 3.14159 * i * (xVal - average) / average));
		fVec(i, 0) = num;

	}

	krigVec << cVec, fVec;

	//std::cout << "krig matrix:\n";
	//std::cout << krigMat << std::endl;
	//std::cout << "krig vec:\n";
	//std::cout << krigVec << std::endl;

	Matrix <std::complex<double>, Dynamic, Dynamic> weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	for (int i = 0; i < size; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0).real());
	}
	//std::cout << "weight vector: " << weightMat << std::endl;
	//std::cout << std::endl;

}

//Taylor Kriging Matrix Solution
static void taylorKriging(std::vector<double> variogramFit, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();

	//if (M > 3) { M = 18; }

	//Matrix <double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd fMat(size, M);
	MatrixXd zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
		}
	}


	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			fMat(i, j) = std::pow((xSample[i].real() - average), j);

		}
	}

	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMat(i, j) = 0.0;
		}
	}

	MatrixXd upper(size, size + M);
	MatrixXd lower(M, size + M);
	MatrixXd krigMat(size + M, size + M);

	upper << sMat, fMat;
	lower << fMat.transpose(), zeroMat;

	krigMat << upper, lower;

	MatrixXd cVec(size, 1);
	MatrixXd fVec(M, 1);
	MatrixXd krigVec(size + M, 1);

	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal - average), i);
	}

	krigVec << cVec, fVec;

	MatrixXd weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	for (int i = 0; i < size; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0));
	}
	//std::cout << weightMat << std::endl;
	//std::cout << std::endl;

}

//Gradient Enhanced Taylor Kriging Matrix Solution
static void taylorKrigingHigherOrder(std::vector<double> variogramFit, std::vector<double> variogramFitFirstD, std::vector<double> variogramFitSecondD, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal, double theta) {
	//double conditionNum = 1.0e8;

	//order of the basis functions is M - 1; differs from literature	
	//M = 1;
	int M = xSample.size() * 2.0;
	//M = 1;
	//if (M > 3)
		//M = 3;

	const int size = xSample.size();
	//MatrixXd upper(size, size + size + M);
	//MatrixXd mid(size, size + size + M);
	//MatrixXd lower(M, size + size + M);
	//MatrixXd krigMat(size + size + M, size + size + M);
	//MatrixXd sMat(size, size);
	//MatrixXd sMatFirstDerivative(size, size);
	//MatrixXd sMatFirstDerivative2(size, size);
	//MatrixXd sMatSecondDerivative(size, size);
	//MatrixXd fMatDerivative(size, M);
	//MatrixXd fMat(size, M);
	//MatrixXd zeroMatM(M, M);
	//MatrixXd zeroMatN(size, size);
	//MatrixXd cVec(size, 1);
	//MatrixXd cVecDerivative(size, 1);
	//MatrixXd fVec(M, 1);
	//MatrixXd krigVec(size + size + M, 1);

	//while (conditionNum > 1e6) {
		//size is N in calculations

		MatrixXd upper(size, size + size + M);
		MatrixXd mid(size, size + size + M);
		MatrixXd lower(M, size + size + M);
		MatrixXd krigMat(size + size + M, size + size + M);
		MatrixXd sMat(size, size);
		MatrixXd sMatFirstDerivative(size, size);
		MatrixXd sMatFirstDerivativeNeg(size, size);
		MatrixXd sMatSecondDerivative(size, size);
		MatrixXd fMatDerivative(size, M);
		MatrixXd fMat(size, M);
		MatrixXd zeroMatM(M, M);
		MatrixXd zeroMatN(size, size);
		MatrixXd cVec(size, 1);
		MatrixXd cVecDerivative(size, 1);
		MatrixXd fVec(M, 1);
		MatrixXd krigVec(size + size + M, 1);

		//calculate sample average
		double sum = 0.0;
		for (int i = 0; i < xSample.size(); ++i) {
			sum += xSample[i].real();
		}
		double average = sum / xSample.size();

		// S matrix
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				double dist = abs(xSample[i].real() - xSample[j].real());
				auto it = std::lower_bound(d.begin(), d.end(), dist);
				int index = it - d.begin();
				sMat(i, j) = variogramFit[index];


				sMatFirstDerivative(i, j) = 2.0 * theta * (xSample[i].real() - xSample[j].real()) * variogramFit[index];
				sMatFirstDerivativeNeg(i, j) = -2.0 * theta * (xSample[i].real() - xSample[j].real()) * variogramFit[index];

				//if (xSample[i].real() - xSample[j].real() > 0) {
				//	sMatFirstDerivative(i, j) = variogramFitFirstD[index];
				//	sMatFirstDerivative2(j, i) = variogramFitFirstD[index];
				//}
				//else {
				//	sMatFirstDerivative(i, j) = variogramFitFirstD[index];
				//	sMatFirstDerivative2(j, i) = variogramFitFirstD[index];
				//}


				//sMatFirstDerivative(i, j) = 0.0;
				//sMatSecondDerivative(i, j) = -variogramFitSecondD[index];
				sMatSecondDerivative(i, j) = -2.0 * theta * variogramFit[index] * (2.0 * theta * xSample[i].real() * xSample[i].real() - 4 * theta * xSample[i].real() * xSample[j].real() + 2 * theta * xSample[j].real() * xSample[j].real() + 1);
			}
		}

		// F Matrix
		double factor;
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < M; ++j) {

				fMat(i, j) = std::pow((xSample[i].real() - average), j);
				fMatDerivative(i, j) = j * std::pow((xSample[i].real() - average), j - 1);

				//if (j == 0)
				//	fMatDerivative(i, j) = 0.0;
				//else
			}
		}

		//MxM zero mat
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < M; ++j) {
				zeroMatM(i, j) = 0.0;
			}
		}

		//NxN zero Mat
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				zeroMatN(i, j) = 0.0;
			}
		}

		upper << sMat, sMatFirstDerivative, fMat;

		mid << -sMatFirstDerivative, sMatSecondDerivative, fMatDerivative;

		lower << fMat.transpose(), fMatDerivative.transpose(), zeroMatM;

		krigMat << upper, mid, lower;

		JacobiSVD<MatrixXd> svd(krigMat);
		double conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
		//std::cout << "condition number: \n" << conditionNum << std::endl;

		//if (conditionNum > 1e6) {
		//	M -= 1;
		//	continue;
		//}
		std::cout << "final M: " << M << std::endl;

		//C vector
		for (int i = 0; i < size; ++i) {
			//std::cout << "i: " << i;
			double xDist = abs(xVal - xSample[i].real());
			auto it = std::lower_bound(d.begin(), d.end(), xDist);
			int index = it - d.begin();
			cVec(i, 0) = variogramFit[index];
			
			cVecDerivative(i, 0) = -2.0 * theta * (xSample[i].real() - xVal) * variogramFit[index];

			//if (xVal > xSample[i].real()) {
			//	cVecDerivative(i, 0) = -variogramFitFirstD[index];
			//}
			//else {
			//	cVecDerivative(i, 0) = variogramFitFirstD[index];
			//}
			
			//cVecDerivative(i, 0) = 0.0;
		}

		for (int i = 0; i < M; ++i) {

			fVec(i, 0) = std::pow((xVal - average), i);
		}

		//std::cout << "krigVec: \n" << krigVec << std::endl;
		krigVec << cVec, cVecDerivative, fVec;


		MatrixXd weightMat = krigMat.inverse() * krigVec;
		//std::cout << "weightMat: " << weightMat << std::endl;

		for (int i = 0; i < size * 2; ++i) {
			//std::cout << weightMat(i, 0) << " ";
			weights.push_back(weightMat(i, 0));
		}

		//for (int i = 0; i < size * 2; ++i) {
		//	//std::cout << weightMat(i, 0) << " ";
		//	std::cout << " heres the real values" << weights[i] << std::endl;
		//}
		//std::cout << xVal << std::endl;

		//std::cout << "weightMat: \n" << weightMat << std::endl;
		//std::cout << std::endl;
	//}



}

//Gradient Enhanced Taylor Kriging Matrix solution with Complex weights/covariance function
static void taylorKrigingHigherOrderComplex(std::vector<std::complex<double>> variogramFit, std::vector<std::complex<double>> variogramFitFirstD, std::vector<std::complex<double>> variogramFitSecondD, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<std::complex<double>>& weights, double xVal, std::vector<double> variogramFitGrad) {

	//size is N in calculations
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();
	//M = 1;
	//if (M > 3) { M = 18; }
	//M = 0;
	//Matrix <double, Dynamic, Dynamic> krigMat;
	MatrixXcd sMat(size, size);
	MatrixXcd sMatFirstDerivative(size, size);
	MatrixXcd sMatFirstDerivative2(size, size);
	MatrixXcd sMatSecondDerivative(size, size);
	MatrixXcd fMatDerivative(size, M);

	MatrixXcd fMat(size, M);
	MatrixXcd zeroMatM(M, M);
	MatrixXcd zeroMatN(size, size);
	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	double average = sum / xSample.size();


	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
			sMatFirstDerivative(i, j) = variogramFitFirstD[index];
			sMatFirstDerivative2(j, i) = variogramFitFirstD[index];

			//sMatFirstDerivative(i, j) = 0.0;
			sMatSecondDerivative(i, j) = variogramFitSecondD[index];
		}
	}

	//sMat = sMat / sMat(0, 0);
	//sMatFirstDerivative = sMatFirstDerivative / sMat(0, 0);
	//sMatFirstDerivative2 = sMatFirstDerivative2 / sMat(0, 0);
	//sMatSecondDerivative = sMatSecondDerivative / sMat(0, 0);

	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {

			fMat(i, j) = std::pow((xSample[i].real() - average), j);

			if (j == 0)
				fMatDerivative(i, j) = 0.0;
			else
				fMatDerivative(i, j) = j * std::pow((xSample[i].real() - average), j - 1);

		}
	}

	//MxM zero mat
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMatM(i, j) = 0.0;
		}
	}

	//NxN zero Mat
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			zeroMatN(i, j) = 0.0;
		}
	}
	//need to fix starting here
	MatrixXcd upper(size, size + size + M);
	MatrixXcd mid(size, size + size + M);
	MatrixXcd lower(M, size + size + M);

	MatrixXcd krigMat(size + size + M, size + size + M);


	upper << sMat, sMatFirstDerivative2, fMat;

	mid << sMatFirstDerivative, sMatSecondDerivative, fMatDerivative;

	lower << fMat.transpose(), fMatDerivative.transpose(), zeroMatM;

	krigMat << upper, mid, lower;

	JacobiSVD<MatrixXcd> svd(krigMat);
	double conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//std::cout << "condition number: \n" << conditionNum << std::endl;
	MatrixXd identity = MatrixXd::Identity(size + size + M, size + size + M);
	if (conditionNum > 7e4) {
		krigMat = krigMat + identity * conditionNum / 2e7;
	}
	JacobiSVD<MatrixXcd> svd2(krigMat);
	conditionNum = svd2.singularValues()(0) / svd2.singularValues()(svd2.singularValues().size() - 1);
	//std::cout << "New condition number: \n" << conditionNum << std::endl;


	MatrixXcd cVec(size, 1);
	MatrixXcd cVecDerivative(size, 1);
	MatrixXcd fVec(M, 1);
	MatrixXcd krigVec(size + size + M, 1);

	//C vector
	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
		cVecDerivative(i, 0) = variogramFitFirstD[index];
		//cVecDerivative(i, 0) = 0.0;
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal - average), i);
	}

	//std::cout << "krigVec: \n" << krigVec << std::endl;
	krigVec << cVec, cVecDerivative, fVec;

	MatrixXcd weightMat = krigMat.inverse() * krigVec;
	std::cout << "weightMat: " << weightMat << std::endl;
	for (int i = 0; i < size * 2; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0));
	}

}

//GETK solution using alternate Kriging solution (seperable mean interpolation and stochastic kriging update terms)
std::complex<double> taylorKrigingHigherOrder3(std::vector<double> variogramFit, std::vector<double> variogramFitFirstD, std::vector<double> variogramFitSecondD, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal, std::vector<double> variogramFitGrad, std::vector<std::complex<double>> yVals, std::vector<std::complex<double>> yValsGrad, double theta) {


	//size is N in calculations
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size()*2-1;
	//M = 1;
		//if (M > 3) { M = 18; }
	//M = 0;
	//Matrix <double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd sMatFirstDerivative(size, size);
	MatrixXd sMatFirstDerivativeNeg(size, size);
	MatrixXd sMatSecondDerivative(size, size);
	MatrixXd fMatDerivative(size, M);

	MatrixXd fMat(size, M);
	MatrixXd zeroMatM(M, M);
	MatrixXd zeroMatN(size, size);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];


			sMatFirstDerivative(i, j) = 2.0 * theta * (xSample[i].real() - xSample[j].real()) * variogramFit[index];
			sMatFirstDerivativeNeg(i, j) = -2.0 * theta * (xSample[i].real() - xSample[j].real()) * variogramFit[index];

			//if (xSample[i].real() - xSample[j].real() > 0) {
			//	sMatFirstDerivative(i, j) = variogramFitFirstD[index];
			//	sMatFirstDerivative2(j, i) = variogramFitFirstD[index];
			//}
			//else {
			//	sMatFirstDerivative(i, j) = variogramFitFirstD[index];
			//	sMatFirstDerivative2(j, i) = variogramFitFirstD[index];
			//}


			//sMatFirstDerivative(i, j) = 0.0;
			//sMatSecondDerivative(i, j) = -variogramFitSecondD[index];
			sMatSecondDerivative(i, j) = -2.0 * theta * variogramFit[index] * (2.0 * theta * xSample[i].real() * xSample[i].real() - 4.0 * theta * xSample[i].real() * xSample[j].real() + 2.0 * theta * xSample[j].real() * xSample[j].real() + 1);
		}
	}

	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {

			fMat(i, j) = std::pow((xSample[i].real() - average), j);
			fMatDerivative(i, j) = j * std::pow((xSample[i].real() - average), j - 1);
			//if (j == 0)
				//fMatDerivative(i, j) = 0.0;
			//else
				

			//if (j == 0) {
			//	fMat(i, j) = 1.0;
			//}
			//else if (j == 1) {
			//	fMat(i, j) = xSample[i].real() - average;
			//}
			//else if (j == 2) {
			//	fMat(i, j) = (xSample[i].real() - average) * (xSample[i].real() - average);
			//}
		}
	}

	//MxM zero mat
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMatM(i, j) = 0.0;
		}
	}

	//NxN zero Mat
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			zeroMatN(i, j) = 0.0;
		}
	}

	//need to fix starting here
	MatrixXd upper(size, size + size + M);
	MatrixXd mid(size, size + size + M);
	MatrixXd lower(M, size + size + M);

	MatrixXd krigMat(size + size + M, size + size + M);


	upper << sMat, -sMatFirstDerivative, fMat;

	mid << sMatFirstDerivative, sMatSecondDerivative, fMatDerivative;

	lower << fMat.transpose(), fMatDerivative.transpose(), zeroMatM;

	krigMat << upper, mid, lower;



	//new method//////////////////////////////////////
	MatrixXd cMatUpper(size, size + size);
	MatrixXd cMatLower(size, size + size);
	MatrixXd cMat(size + size, size + size);
	cMatUpper << sMat, -sMatFirstDerivative;
	cMatLower << sMatFirstDerivative, sMatSecondDerivative;
	cMat << cMatUpper, cMatLower;
	//std::cout << "cMat: \n" << cMat << std::endl;

	MatrixXd bMat(size+size, M);
	bMat << fMat, fMatDerivative;
	//////////////////////////////////////////////////


	MatrixXd cVec(size, 1);
	MatrixXd cVecDerivative(size, 1);
	MatrixXd fVec(M, 1);
	MatrixXd krigVec(size + size + M, 1);
	MatrixXd yVec(size + size, 1);

	//C vector
	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];

		cVecDerivative(i, 0) = 2.0 * theta * (xSample[i].real() - xVal) * variogramFit[index];

		yVec(i, 0) = abs(yVals[i]);
		yVec(i + size, 0) = abs(yValsGrad[i]);


		//if (xVal > xSample[i].real()) {
		//	cVecDerivative(i, 0) = -variogramFitFirstD[index];
		//}
		//else {
		//	cVecDerivative(i, 0) = variogramFitFirstD[index];
		//}

		//cVecDerivative(i, 0) = 0.0;
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal - average), i);
	}



	//More new method//////////////////////////
	krigVec << cVec, cVecDerivative, fVec;

	MatrixXd cVecFull(size + size, 1);
	MatrixXd bVec(M, 1);

	cVecFull << cVec, cVecDerivative;
	bVec << fVec;


	VectorXd beta(M, 1);
	beta = (bMat.transpose() * cMat.inverse() * bMat).inverse() * bMat.transpose() * cMat.inverse() * yVec;

	//std::cout << "beta: \n" << beta << std::endl;

	//std::cout << "cMat: \n" << cMat << std::endl;
	//std::cout << "bMat^T*cMat^-1: \n" << bMat.transpose()*cMat.inverse() << std::endl;

	MatrixXd y(1, 1);
	double betaConstant = beta(0, 0);

	VectorXd ones = VectorXd::Ones(size);
	VectorXd zeros = VectorXd::Zero(size);
	VectorXd onesZeros(size + size);
	onesZeros << ones, zeros;



	//average term
	y = beta.transpose() * bVec;
	//full kriging
	y = beta.transpose() * bVec + cVecFull.transpose() * cMat.inverse() * (yVec - bMat * beta);
	//full kriging
	// 
	// 
	std::cout << "kriging update term:\n" << cVecFull.transpose() * cMat.inverse() * (yVec - bMat * beta) << std::endl;
	//y = (bVec - cVec.transpose() * cMat.inverse() * bMat) * beta + cVec.transpose() * yVec;
	//std::cout << "update term: \n" << cVecFull.transpose() * cMat.inverse() * (yVec - s) << std::endl;
	//y(0,0) = 0.0;
	

	//////////////////////////////////////



	MatrixXd weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	for (int i = 0; i < size * 2; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0));
	}

	//for (int i = 0; i < size * 2; ++i) {
	//	//std::cout << weightMat(i, 0) << " ";
	//	std::cout << " heres the real values" << weights[i] << std::endl;
	//}
	//std::cout << xVal << std::endl;
	//std::cout << weightMat << std::endl;
	//std::cout << std::endl;

	return y(0, 0);

}

//sum of weights*reference points for final output of kriging
static std::complex<double> krigSum(std::vector<double> weights, std::vector<std::complex<double>> referenceVals) {
	std::complex<double> sum = 0.0;
	for (int i = 0; i < weights.size(); ++i) {
		sum += weights[i] * referenceVals[i];
	}
	return sum;
}

//sum of weights*ref + gradWeights*gradRef
static std::complex<double> krigSumHigherOrder(std::vector<double> weights, std::vector<std::complex<double>> referenceVals, std::vector<std::complex<double>> gradient, double xVal, double xAvg) {
	std::complex<double> sum = 0.0;
	std::complex<double> sum2 = 0.0;


	//std::cout << "started krigSumHigher, referenceVals.size()  "  << referenceVals.size() << " gradient.size(): " << gradient.size() << "  weights.size(): " << weights.size() << std::endl;

	for (int i = 0; i < weights.size(); ++i) {

		if (i < referenceVals.size()) {
			sum += referenceVals[i] * weights[i];
		}
		else
			sum += gradient[i-referenceVals.size()] * weights[i];
		//sum += referenceVals[i] * weights[i] + gradient[i] * weights[i + referenceVals.size()];// *(xVal - xAvg);

		//std::cout << "referenceVals[i]: " << referenceVals[i] << "  weights[i]: " << weights[i] << "  gradient[i]: " << gradient[i] << "  weights[i + size]: " << weights[i + referenceVals.size()] << std::endl;

	}


	return sum;
}

//GETK Kriging weights sum for complex weights
static std::complex<double> krigSumHigherOrderComplex(std::vector<std::complex<double>> weights, std::vector<std::complex<double>> referenceVals, std::vector<std::complex<double>> gradient, double xVal, double xAvg) {
	std::complex<double> sum = 0.0;
	std::complex<double> sum2 = 0.0;


	//std::cout << "started krigSumHigher, referenceVals.size()  "  << referenceVals.size() << " gradient.size(): " << gradient.size() << "  weights.size(): " << weights.size() << std::endl;

	for (int i = 0; i < weights.size(); ++i) {

		if (i < referenceVals.size()) {
			sum += referenceVals[i] * weights[i];
		}
		else
			sum += gradient[i - referenceVals.size()] * weights[i];
		//sum += referenceVals[i] * weights[i] + gradient[i] * weights[i + referenceVals.size()];// *(xVal - xAvg);

		//std::cout << "referenceVals[i]: " << referenceVals[i] << "  weights[i]: " << weights[i] << "  gradient[i]: " << gradient[i] << "  weights[i + size]: " << weights[i + referenceVals.size()] << std::endl;

	}


	return sum;
}

//trim nan values from empirical variogram/covariance function
static void trimVariogram(std::vector<double>& variogram, std::vector<double>& edges) {
	int idx = 0;
	std::vector<double> variogramNew;
	std::vector<double> edgesNew;

	//std::cout << "edges size: " << edges.size() << " variogram size: " << variogram.size() << std::endl;

	for (int i = 0; i < edges.size() - 1; i++) {
		//std::cout << "edges[i]: " << edges[i];
		if (!isnan(double(variogram[i]))) {
			variogramNew.push_back(variogram[i]);
			edgesNew.push_back(edges[i+1]);
		}
	}
	//std::cout << std::endl;

	for (int i = 0; i < edgesNew.size(); ++i) {
		//std::cout << "edgesNew: " << edgesNew[i] << " variogram new: " << variogramNew[i];
	}
	//std::cout << std::endl;

	variogram = variogramNew;
	edges = edgesNew;

	//variogram[0] = 0.0;

}

//trim nan values from complex empirical covariance/variogram
static void trimVariogramComplex(std::vector<std::complex<double>>& variogram, std::vector<double>& edges) {
	int idx = 0;
	std::vector<std::complex<double>> variogramNew;
	std::vector<double> edgesNew;

	//std::cout << "edges size: " << edges.size() << " variogram size: " << variogram.size() << std::endl;

	for (int i = 0; i < edges.size() - 1; i++) {
		//std::cout << "edges[i]: " << edges[i];
		if (!isnan(double(variogram[i].real())) && !isnan(double(variogram[i].real()))) {
			variogramNew.push_back(variogram[i]);
			edgesNew.push_back(edges[i + 1]);
		}
	}
	//std::cout << std::endl;

	for (int i = 0; i < edgesNew.size(); ++i) {
		//std::cout << "edgesNew: " << edgesNew[i] << " variogram new: " << variogramNew[i];
	}
	//std::cout << std::endl;

	variogram = variogramNew;
	edges = edgesNew;

	//variogram[0] = 0.0;

}

//gaussian fit for empirical variogram
void function_fit_gaussian(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(-x[0] * x[0] / (c[1] * c[1])));
}

//gaussian fit for empirical covariance
void function_fit_gaussianCo(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (std::exp(-x[0] * x[0] / (c[1] * c[1])));
}

//exponential fit for empirical variogram
void function_fit_exp(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(-x[0] / (c[1])));
}

//exponential fit for empirical covariance
void function_fit_expCo(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (std::exp(-x[0] / (c[1])));
}

//stable fit for empirical variogram
void function_fit_stable(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(std::pow(-x[0] / c[1], int(c[3]))));
}

//linear fit for empirical variogram
void function_fit_linear(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[0] * x[0] + c[1];
}

//cubic fit for empirical variogram
void function_fit_cubicNew(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = 1.0 / 3.0 * std::pow((1.0 - x[0] * c[0]), 6) * (35.0 * x[0] * x[0] * c[0] * c[0] + 180.0 * x[0] * c[0] + 3);
}

//spherical fit for empirical variogram
void function_fit_spherical(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	if (x[0] < c[1])
		func = c[2] + c[0] * (1.5 * x[0] / c[1] - 0.5 * x[0] * x[0] * x[0] / c[1] / c[1] / c[1]);
	else
		func = c[2] + c[0];
}


//fit one of the models for variogram
std::vector<double> fitVariogramModel(std::vector<double> variogram, std::vector<double> edges, std::vector<double> d, double size) {
	alglib::real_2d_array x;
	alglib::real_1d_array y;
	y.setlength(variogram.size());
	x.setlength(variogram.size(), 1);

	for (int i = 0; i < variogram.size(); ++i) {
		y[i] = variogram[i];
		x[i][0] = edges[i];
	}

	double epsx = 0.000001;
	//double epsx = 0;
	alglib::ae_int_t maxits = 0;
	alglib::ae_int_t info;
	alglib::lsfitstate state;
	alglib::lsfitreport rep;
	double diffstep = 0.0001;
	std::vector<double> fitVariogram(d.size());


		//////////////////////////linear fit////////////////////////////////////////////
		//alglib::real_1d_array c;
		//c.setlength(2);
		//c[0] = *std::max_element(variogram.begin(), variogram.end());
		//c[1] = edges[edges.size() - 1] / 2.0;

		//alglib::real_1d_array bndl;
		//alglib::real_1d_array bndu;
		//bndl.setlength(2);
		//bndu.setlength(2);
		//bndl[0] = 0.0;
		//bndl[1] = 0.0;
		//bndu[0] = c[0];
		//bndu[1] = c[1] * 2.0;

		//alglib::lsfitcreatef(x, y, c, diffstep, state);
		//alglib::lsfitsetcond(state, epsx, maxits);
		//alglib::lsfitsetbc(state, bndl, bndu);
		//alglib::lsfitfit(state, function_fit_linear);
		//alglib::lsfitresults(state, info, c, rep);
		////std::cout << "c: " << c[0] << "  a: " << c[1] << std::endl;


		//for (int i = 0; i < d.size(); ++i) {
		//	fitVariogram[i] = c[0] * d[i] + c[1];
		//}


		///////////////////////////squared exp fit//////////////////////////////////////////
	
		//alglib::real_1d_array c;
		//c.setlength(3);

		////sill parameter
		//c[0] = *std::max_element(variogram.begin(), variogram.end());
		////range parameter
		//c[1] = edges[edges.size() - 1] / 2.0;
		////y offset
		//c[2] = 0.0;

		//alglib::real_1d_array bndl;
		//alglib::real_1d_array bndu;
		//bndl.setlength(3);
		//bndu.setlength(3);
		//bndl[0] = 0.0;
		//bndl[1] = 0.0;
		//bndl[2] = 0.0;
		//bndu[0] = c[0];
		//bndu[1] = c[1] * 2.0;
		//bndu[2] = bndu[1] * 0.99;

		//alglib::lsfitcreatef(x, y, c, diffstep, state);
		//alglib::lsfitsetcond(state, epsx, maxits);
		//alglib::lsfitsetbc(state, bndl, bndu);
		//alglib::lsfitfit(state, function_fit_gaussianCo);
		//alglib::lsfitresults(state, info, c, rep);
		////std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		//for (int i = 0; i < d.size(); ++i) {
		//	fitVariogram[i] = c[2] + c[0] * (std::exp(-d[0] * d[0] / (c[1] * c[1])));

		//	if (i != 0 && isnan(fitVariogram[i]))
		//		fitVariogram[i] = fitVariogram[i - 1];
		//}
		//if (isnan(fitVariogram[0]))
		//	fitVariogram[0] = fitVariogram[1];

		/////////////////////////////gaussian fit//////////////////////////////////////////
	
		alglib::real_1d_array c;
		c.setlength(3);

		//sill parameter
		c[0] = *std::max_element(variogram.begin(), variogram.end());
		//range parameter
		c[1] = edges[edges.size() - 1] / 2.0;
		//y offset
		c[2] = 0.0;

		alglib::real_1d_array bndl;
		alglib::real_1d_array bndu;
		bndl.setlength(3);
		bndu.setlength(3);
		bndl[0] = 0.0;
		bndl[1] = 0.0;
		bndl[2] = 0.0;
		bndu[0] = c[0];
		bndu[1] = c[1] * 2.0;
		bndu[2] = bndu[1] * 0.99;

		alglib::lsfitcreatef(x, y, c, diffstep, state);
		alglib::lsfitsetcond(state, epsx, maxits);
		alglib::lsfitsetbc(state, bndl, bndu);
		alglib::lsfitfit(state, function_fit_gaussian);
		alglib::lsfitresults(state, info, c, rep);
		//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		for (int i = 0; i < d.size(); ++i) {
			fitVariogram[i] = c[2] + c[0] * (1 - std::exp(-d[i] * d[i] / c[1] / c[1]));

			if (i != 0 && isnan(fitVariogram[i]))
				fitVariogram[i] = fitVariogram[i - 1];
		}
		if (isnan(fitVariogram[0]))
			fitVariogram[0] = fitVariogram[1];

		
		/////////////////////////////	//Stable fit///////////////////////////////////////////////////
		//alglib::real_1d_array c;
		//c.setlength(4);

		////sill parameter
		//c[0] = *std::max_element(variogram.begin(), variogram.end());
		////range parameter
		//c[1] = edges[edges.size() - 1] / 3.0;
		////y offset
		//c[2] = 0.0;
		//c[3] = 2.0;

		//alglib::real_1d_array bndl;
		//alglib::real_1d_array bndu;
		//bndl.setlength(4);
		//bndu.setlength(4);
		//bndl[0] = 0.0;
		//bndl[1] = 0.0;
		//bndl[2] = 0.0;
		//bndl[3] = 0.0;
		//bndu[0] = c[0];
		//bndu[1] = c[1] * 3.0;
		//bndu[2] = bndu[1] * 0.99;
		//bndu[3] = 10.0;

		//alglib::lsfitcreatef(x, y, c, diffstep, state);
		//alglib::lsfitsetcond(state, epsx, maxits);
		//alglib::lsfitsetbc(state, bndl, bndu);
		//alglib::lsfitfit(state, function_fit_stable);
		//alglib::lsfitresults(state, info, c, rep);
		//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		//for (int i = 0; i < d.size(); ++i) {
		//	fitVariogram[i] = c[2] + c[0] * (1.0 - std::exp(std::pow(-d[i] / c[1], int(c[3]))));

		//}

		////////////////////exponential fit////////////////////////////////////////
		//alglib::real_1d_array c;
		//c.setlength(3);

		////sill parameter
		//c[0] = *std::max_element(variogram.begin(), variogram.end());
		////range parameter
		//c[1] = edges[edges.size() - 1] / 3.0;
		////y offset
		//c[2] = 0.0;

		//alglib::real_1d_array bndl;
		//alglib::real_1d_array bndu;
		//bndl.setlength(3);
		//bndu.setlength(3);
		//bndl[0] = 0.0;
		//bndl[1] = 0.0;
		//bndl[2] = 0.0;
		//bndu[0] = c[0];
		//bndu[1] = c[1] * 3.0;
		//bndu[2] = bndu[1] * 0.99;

		//alglib::lsfitcreatef(x, y, c, diffstep, state);
		//alglib::lsfitsetcond(state, epsx, maxits);
		//alglib::lsfitsetbc(state, bndl, bndu);
		//alglib::lsfitfit(state, function_fit_exp);
		//alglib::lsfitresults(state, info, c, rep);
		//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		//for (int i = 0; i < d.size(); ++i) {
		//	fitVariogram[i] = c[2] + c[0] * (1 - std::exp(-d[i] / c[1]));

		//	if (i != 0 && isnan(fitVariogram[i]))
		//		fitVariogram[i] = fitVariogram[i - 1];
		//}
		//if (isnan(fitVariogram[0]))
		//	fitVariogram[0] = fitVariogram[1];

		//////////////////////////////////spherical fit///////////////////////////////////////
		/*alglib::real_1d_array c;
		 c.setlength(3);
		c[0] = *std::max_element(variogram.begin(), variogram.end());
		c[1] = edges[edges.size() - 1];
		c[2] = 0.0;
 		alglib::real_1d_array bndl;
		alglib::real_1d_array bndu;
		bndl.setlength(3);
		bndu.setlength(3);

		bndl[0] = *std::min_element(variogram.begin(), variogram.end());
		bndl[1] = edges[0];
		bndl[2] = 0.0;
		bndu[0] = c[0];
		bndu[1] = c[1];
		bndu[2] = bndu[1]*0.99;

		alglib::lsfitcreatef(x, y, c, diffstep, state);
		alglib::lsfitsetcond(state, epsx, maxits);
		alglib::lsfitsetbc(state, bndl, bndu);
		alglib::lsfitfit(state, function_fit_spherical);
		alglib::lsfitresults(state, info, c, rep);
		std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		for (int i = 0; i < d.size(); ++i) {
			if (d[i] < c[1])
				fitVariogram[i] = c[2] + c[0] * (1.5 * d[i] / c[1] - 0.5 * d[i] * d[i] * d[i] / c[1] / c[1] / c[1]);
			else
				fitVariogram[i] = c[2] + c[0];

		}*/

		//-////////////////////////////////////////////////////////////////////////////////////////////
	double max = *std::max_element(fitVariogram.begin(), fitVariogram.end());
	for (int i = 0; i < fitVariogram.size(); ++i) {
		fitVariogram[i] = fitVariogram[i] / max;
	}

	return fitVariogram;

}

//fit one of the models for covariance function
std::vector<double> fitCoVariogramModel(std::vector<double> variogram, std::vector<double> edges, std::vector<double> d, double size) {
	alglib::real_2d_array x;
	alglib::real_1d_array y;
	y.setlength(variogram.size());
	x.setlength(variogram.size(), 1);

	for (int i = 0; i < variogram.size(); ++i) {
		y[i] = variogram[i];
		x[i][0] = edges[i];
	}

	std::cout << "start of fit method\n";

	double epsx = 0.00001;
	//double epsx = 0;
	alglib::ae_int_t maxits = 0;
	alglib::ae_int_t info;
	alglib::lsfitstate state;
	alglib::lsfitreport rep;
	double diffstep = 0.0001;
	std::vector<double> fitVariogram(d.size());


	//////////////////////////linear fit////////////////////////////////////////////
	
	//alglib::real_1d_array c = "[0.25, 0.25]";
	//alglib::lsfitcreatef(x, y, c, diffstep, state);
	//alglib::lsfitsetcond(state, epsx, maxits);
	//alglib::lsfitfit(state, function_fit_linear);
	//alglib::lsfitresults(state, info, c, rep);
	//std::cout << "c: " << c[0] << "  a: " << c[1] << std::endl;


	//for (int i = 0; i < d.size(); ++i) {
	//	fitVariogram[i] = c[0] * d[i] + c[1];
	//}
	
	/////////////////////////gaussian fit//////////////////////////////////////////
	
	alglib::real_1d_array c;
	c.setlength(3);

	//sill parameter
	c[0] = *std::max_element(variogram.begin(), variogram.end());
	//range parameter
	c[1] = edges[edges.size() - 1] / 2.0;
	//y offset
	c[2] = 0.0;

	alglib::real_1d_array bndl;
	alglib::real_1d_array bndu;
	bndl.setlength(3);
	bndu.setlength(3);
	bndl[0] = 0.0;
	bndl[1] = 0.0;
	bndl[2] = 0.0;
	bndu[0] = c[0];
	bndu[1] = c[1] * 0.7;
	bndu[2] = bndu[1] * 0.99;


	alglib::lsfitcreatef(x, y, c, diffstep, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetbc(state, bndl, bndu);

	alglib::lsfitfit(state, function_fit_gaussianCo);
	alglib::lsfitresults(state, info, c, rep);

	std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	for (int i = 0; i < d.size(); ++i) {
		fitVariogram[i] = c[2] + c[0] * (std::exp(-d[i] * d[i] / c[1] / c[1]));

		if (i != 0 && isnan(fitVariogram[i]))
			fitVariogram[i] = fitVariogram[i - 1];
	}
	if (isnan(fitVariogram[0]))
		fitVariogram[0] = fitVariogram[1];
	
	////////////////////////Cubic fit////////////////////////////////////////////////
	//alglib::real_1d_array c;
	//c.setlength(1);

	//alglib::lsfitcreatef(x, y, c, diffstep, state);
	//alglib::lsfitsetcond(state, epsx, maxits);
	////alglib::lsfitsetbc(state, bndl, bndu);
	//alglib::lsfitfit(state, function_fit_cubicNew);
	//alglib::lsfitresults(state, info, c, rep);
	//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	//for (int i = 0; i < d.size(); ++i) {
	//	if (d[i] < 1.0 / c[0])
	//		fitVariogram[i] = 1.0 / 3.0 * std::pow((1.0 - d[i]*c[0]), 6)*(35.0 * d[i] * d[i] * c[0] * c[0] + 180.0 * d[i] * c[0] + 3);
	//	else
	//		fitVariogram[i] = 0.0;

	//}


	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////	//Stable fit///////////////////////////////////////////////////
	//alglib::real_1d_array c;
	//c.setlength(4);

	////sill parameter
	//c[0] = *std::max_element(variogram.begin(), variogram.end());
	////range parameter
	//c[1] = edges[edges.size() - 1] / 3.0;
	////y offset
	//c[2] = 0.0;
	//c[3] = 2.0;

	//alglib::real_1d_array bndl;
	//alglib::real_1d_array bndu;
	//bndl.setlength(4);
	//bndu.setlength(4);
	//bndl[0] = 0.0;
	//bndl[1] = 0.0;
	//bndl[2] = 0.0;
	//bndl[3] = 0.0;
	//bndu[0] = c[0];
	//bndu[1] = c[1] * 3.0;
	//bndu[2] = bndu[1] * 0.99;
	//bndu[3] = 10.0;

	//alglib::lsfitcreatef(x, y, c, diffstep, state);
	//alglib::lsfitsetcond(state, epsx, maxits);
	//alglib::lsfitsetbc(state, bndl, bndu);
	//alglib::lsfitfit(state, function_fit_stable);
	//alglib::lsfitresults(state, info, c, rep);
	//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	//for (int i = 0; i < d.size(); ++i) {
	//	fitVariogram[i] = c[2] + c[0] * (1.0 - std::exp(std::pow(-d[i] / c[1], int(c[3]))));

	//}

	////////////////////exponential fit////////////////////////////////////////
	//alglib::real_1d_array c;
	//c.setlength(3);

	////sill parameter
	//c[0] = *std::max_element(variogram.begin(), variogram.end());
	////range parameter
	//c[1] = edges[edges.size() - 1] / 3.0;
	////y offset
	//c[2] = 0.0;

	//alglib::real_1d_array bndl;
	//alglib::real_1d_array bndu;
	//bndl.setlength(3);
	//bndu.setlength(3);
	//bndl[0] = 0.0;
	//bndl[1] = 0.0;
	//bndl[2] = 0.0;
	//bndu[0] = c[0];
	//bndu[1] = c[1] * 3.0;
	//bndu[2] = bndu[1] * 0.99;

	//alglib::lsfitcreatef(x, y, c, diffstep, state);
	//alglib::lsfitsetcond(state, epsx, maxits);
	////alglib::lsfitsetbc(state, bndl, bndu);
	//alglib::lsfitfit(state, function_fit_expCo);
	//alglib::lsfitresults(state, info, c, rep);
	//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	//for (int i = 0; i < d.size(); ++i) {
	//	fitVariogram[i] = c[2] + c[0] * (std::exp(-d[i] / c[1]));

	//	if (i != 0 && isnan(fitVariogram[i]))
	//		fitVariogram[i] = fitVariogram[i - 1];
	//}
	//if (isnan(fitVariogram[0]))
	//	fitVariogram[0] = fitVariogram[1];

	//////////////////////////////////spherical fit///////////////////////////////////////
	/*alglib::real_1d_array c;
	 c.setlength(3);
	c[0] = *std::max_element(variogram.begin(), variogram.end());
	c[1] = edges[edges.size() - 1];
	c[2] = 0.0;
	alglib::real_1d_array bndl;
	alglib::real_1d_array bndu;
	bndl.setlength(3);
	bndu.setlength(3);

	bndl[0] = *std::min_element(variogram.begin(), variogram.end());
	bndl[1] = edges[0];
	bndl[2] = 0.0;
	bndu[0] = c[0];
	bndu[1] = c[1];
	bndu[2] = bndu[1]*0.99;

	alglib::lsfitcreatef(x, y, c, diffstep, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetbc(state, bndl, bndu);
	alglib::lsfitfit(state, function_fit_spherical);
	alglib::lsfitresults(state, info, c, rep);
	std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	for (int i = 0; i < d.size(); ++i) {
		if (d[i] < c[1])
			fitVariogram[i] = c[2] + c[0] * (1.5 * d[i] / c[1] - 0.5 * d[i] * d[i] * d[i] / c[1] / c[1] / c[1]);
		else
			fitVariogram[i] = c[2] + c[0];

	}*/
	//////////////////////////////////////////////////////////////////////////////

	double max = *std::max_element(fitVariogram.begin(), fitVariogram.end());
	for (int i = 0; i < fitVariogram.size(); ++i) {
		//fitVariogram[i] = fitVariogram[i] / max;
	}

	return fitVariogram;

}

//data struct for MLE
struct data
{
	const MatrixXd* X;
	const VectorXd* y;
	const MatrixXd* B;
	const VectorXd* yGrad;
};

//MLE covariance function
double covariance(double d, double theta, double sigma2) {
	//make sure to change all 3
	//gaussian fit
	return sigma2 * exp(-d * d * theta);
	//exp fit
	//return sigma2 * exp(-d * theta);
}

//calculate a vector
//a vector contains the coefficients of the regression for each corresponding basis function
//when a vector is multiplied by basis function vector this will give the regression (but no Kriging update)
void aCalc(MatrixXd X, VectorXd y, MatrixXd B, double theta, VectorXd &a) {

	//calculate correlation matrix R
	int n = X.rows();
	MatrixXd R = MatrixXd::Zero(n, n);
	VectorXd ones = VectorXd::Ones(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double dist = (X.row(i) - X.row(j)).norm();
			R(i, j) = exp(-theta * dist);
		}
	}
	//this is actually mu(x)
	//lower order version
	a = (ones.transpose() * R.inverse() * ones).inverse()* ones.transpose()* R.inverse()* y;
	//higher order version
	//a = (B.transpose() * R.inverse() * B.transpose()).inverse() * B.transpose() * R.inverse() * y;

}

//basis function matrix B
void basisCalc(MatrixXd X, VectorXd y, MatrixXd& B, VectorXd yGrad) {
	int n = y.rows();
	int m = n;
	MatrixXd R = MatrixXd::Zero(n, m);
	//calculate sample average
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += X(i,0);
	}
	double average = sum / n;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {

			B(i, j) = std::pow((X(i, 0) - average), j);
		}
	}

}

//best variance calculation (using MLE hyperparameters)
void var_calc(MatrixXd X, VectorXd y, double& sigma2, double theta, MatrixXd B, VectorXd yGrad) {

	//calculate correlation matrix R
	int n = X.rows();
	MatrixXd R = MatrixXd::Zero(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double dist = (X.row(i) - X.row(j)).norm();
			//gaussian fit
			R(i, j) = exp(-theta * dist * dist);
			//exp fit
			//R(i, j) = exp(-theta * dist);
		}
	}

	JacobiSVD<MatrixXcd> svd(R);
	double conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

	//while (conditionNum > 1e3) {
	//	R += MatrixXd::Identity(n, n) * 1e-5;
	//	JacobiSVD<MatrixXcd> svd(R);
	//	conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//}
	
	///////////no basis functions///////////
	VectorXd ones = VectorXd::Ones(n);
	double mu = (ones.transpose() * R.inverse() * ones).inverse() * ones.transpose() * R.inverse() * y;
	double y_Rinv_y = (y - ones*mu).transpose() * R.inverse() * (y - ones*mu);

	//////////basis function version://////////////////////
	//VectorXd a = (B.transpose() * R.inverse() * B).inverse() * B.transpose() * R.inverse() * y;
	//double y_Rinv_y = (y - a).transpose() * R.inverse() * (y - a);
	////////////////////////////////////////////////////////



	sigma2 = y_Rinv_y / n;
	//sigma2 = (y - B * a.transpose()).transpose();
	//sigma2 = (y - B * a).transpose() * R.inverse() * (y - B * a);
}

//objective function to be minimized.  Comes from loglikelihood min terms
double obj_func(unsigned n, const double* x, double* grad, void* data_ptr)
{
	data* d = reinterpret_cast<data*>(data_ptr);
	const MatrixXd& X = *(d->X);
	const VectorXd& y = *(d->y);
	const MatrixXd& B = *(d->B);
	double theta = x[0];
	int m = X.rows();
	MatrixXd R = MatrixXd::Zero(m, m);

	// Calculate the correlation matrix using a Gaussian correlation function
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double dist = (X.row(i) - X.row(j)).norm();
			//gauss fit
			R(i, j) = exp(-theta * dist * dist);
			//exp fit
			//R(i, j) = exp(-theta * dist);
		}
	}

	JacobiSVD<MatrixXcd> svd(R);
	double conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//std::cout << "condition number for R: \n" << conditionNum << std::endl;

	//while (conditionNum > 1e3) {
	//	R += MatrixXd::Identity(m, m) * 1e-5;
	//	JacobiSVD<MatrixXcd> svd(R);
	//	conditionNum = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	//	std::cout << "NEW condition number for R: \n" << conditionNum << std::endl;

	//}
	//// Add a small diagonal term for numerical stability
	////R += 1e-8 * MatrixXd::Identity(m, m);

	////conditional matrix conditioning
	////if (conditionNum > 7e4) {
	////	R += conditionNum / 1e5 * MatrixXd::Identity(m,m);
	////}

	////different form of conditional matrix conditioning
	//R += conditionNum / 1e8 * MatrixXd::Identity(m, m);
	//JacobiSVD<MatrixXcd> svd2(R);
	//conditionNum = svd2.singularValues()(0) / svd2.singularValues()(svd2.singularValues().size() - 1);
	//std::cout << "New condition number for R: \n" << conditionNum << std::endl;

	// Calculate the negative log likelihood


	//double log_det = 2.0 * (R.diagonal().array().log().sum());

	double log_det = log(R.determinant());

	///////////no basis functions///////////
	//VectorXd ones = VectorXd::Ones(m);
	//double mu = (ones.transpose() * R.inverse() * ones).inverse() * ones.transpose() * R.inverse() * y;
	//double y_Rinv_y = (y - ones*mu).transpose() * R.inverse() * (y - ones*mu);
	
	//////////basis function version://////////////////////
	VectorXd a = (B.transpose() * R.inverse() * B).inverse() * B.transpose() * R.inverse() * y;
	double y_Rinv_y = (y - a).transpose() * R.inverse() * (y - a);
	////////////////////////////////////////////////////////

	//variance
	double sigma2 = y_Rinv_y / m;

	//maximization function
	double log_likelihood = 0.5 * (log_det + log(sigma2)*m);// +m * log(2 * 3.14159));
	return log_likelihood;
}

//MLE minimization for hyperparameter terms, using NLOpt c++ package
void MLE_estimation(const MatrixXd& X, const VectorXd& y, double& theta, MatrixXd& B, VectorXd yGrad)
{
	int n = X.cols();
	data d = { &X, &y, &B, &yGrad };

	// Create the optimization problem
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_SBPLX, 1);
	//opt.set_min_objective(obj_func, &d);

	// Set the bounds and initial guess for the parameters
	double lb[1] = { 1e-5 };
	double ub[1] = { 1e5 };
	double x0[1] = { 1.0 };

	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, obj_func, &d);
	nlopt_set_xtol_rel(opt, 1e-6);

	// Optimize the parameters
	double min_f;
	if (nlopt_optimize(opt, x0, &min_f) < 0)
	{
		std::cout << "NLOPT optimization failed!\n" ;
	}
	else
	{
		theta = x0[0];
	}
}

//calculates qoi and qoi gradients
std::complex<double> Kriging::sensitivity_to_epsr(std::string & file_name, std::vector<std::complex<double>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::complex<double> referenceFreq, std::complex<double>& gradient)
{
	//get orig QoI
	//get entries of Bint
	//get part of Gint
	bool check_q = false;
	bool plot = false;
	bool sens = false;
	bool error = false;
	bool refine = false;
	bool basis_error = false;
	bool useAdjoint = false;
	bool higher_order = false; 
	bool check_results = false;
	//mesh inputs
	std::string mesh_name = file_name;
	//get the forward solution
	auto t1 = std::chrono::high_resolution_clock::now();
	Domain dom_forward(mesh_name);

	run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, false, referenceFreq.real());
	std::cout << "Forward solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count() << std::endl;
	//get the adjoint solution
	Domain dom_adjoint(mesh_name);
	run::_standard(dom_adjoint, mesh_name, plot, true, higher_order, check_results, false, referenceFreq.real());
	//get the reference qoi
	double k0 = dom_forward.scatter1.K0[1];
	std::complex<double> reference_qoi = get_QoI(dom_forward, dom_adjoint);
	//std::complex<double> reference_qoi_compare = dom_forward.scatter1.
	std::cout << "reference qoi: " << reference_qoi << " reference mat: " << referenceFreq.real() << std::endl;
	
	std::complex<double> eps_diff = 0.0;
	std::complex<double> reference_material = 0.0;
	double reference_k0;
	std::complex<double> new_qoi;
	for (auto e = dom_forward.elements.begin(); e != dom_forward.elements.end(); ++e) {
		if (e->materials.region == 1) {
			reference_material = e->materials.epsr_list[0][0]; //only works for homogeneous scatterers at the moment!
			break;
		}
	}

	//get the remainder of the estimate
	std::complex<double> dQoI = 0.0;

	std::complex<double> dQoI2 = 0.0, dQoI1 = 0.0, dQoI2_ver1 = 0.0, dQoI2_ver2 = 0.0;
	//loop through the stored Bint pieces
	for (auto e = dom_forward.elements.begin(); e != dom_forward.elements.end(); ++e) {

		////////////////////Material Perturbation//////////////////////////////////////
		if (e->materials.region == 1) {
			// ------------------------------------------------------------------------------
			for (int i_unkn = e->unknownsStart; i_unkn <= e->unknownsEnd; ++i_unkn) {
				int icon = abs(dom_forward.vectorD[i_unkn]);
				dQoI += e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
			//	dQoI1 += eps_diff * e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
				for (int j_unkn = e->unknownsStart; j_unkn <= e->unknownsEnd; ++j_unkn) {
					int jcon = abs(dom_adjoint.vectorD[j_unkn]);
					int i0 = i_unkn - e->unknownsStart + 1;
					int j0 = j_unkn - e->unknownsStart + 1;
					//dQoI += eps_diff * k0*k0*e->cBrPAK(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]);
					//dQoI2_ver2 += eps_diff * k0*k0*e->cBrPAK(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]); //this version is much worse, about double the error
				//	dQoI += eps_diff * e->cBr_HOPS(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]);
					//dQoI2_ver1 += eps_diff * e->cBr_HOPS(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]); //good version
					dQoI += e->cBr_HOPS(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]);
				}
			}
			/////////////////Material Perturbation///////////////////////////////////////////////////
		}
		//--------------------------------------------------------------------------------------
	}

	std::cout << "dQoI: " << dQoI << std::endl;
	gradient = dQoI;


	double kReal = 2 * 3.14159 / (2.99E8 / referenceFreq.real());
	std::complex<double> referenceK = (kReal, 0.0);

	//now have the base QoI modifier, need to scale by the change in perturbed material parameter
	//for (int i = 0; i < epsr_list.size(); ++i) {
	//	/////////////////////////////////////////Frequency Perturbation/////////////////////////////////////////////
	//	//double k0_list = epsr_list[i].real() * 2.0 * 3.14159 / 3e8;
	//	//eps_diff = -1.0*k0_list + k0;
	//	//eps_diff = (k0_list / k0 - 1.0) * k0;
	//	//---------------------------------------------------------------------------------------------------------

	//	///////////////////////////////////////Material Perturbation//////////////////////////////////////////////
	//	eps_diff = epsr_list[i] - referenceFreq;
	//	///-------------------------------------------------------------------------------------------------------


	//	updated_qoi.push_back(reference_qoi + eps_diff * dQoI);// +eps_diff * eps_diff * dQoI2 / 2.0);


	//	//std::cout << reference_qoi + eps_diff * dQoI << std::endl;
	//	//std::cout << "qoi: " << reference_qoi + eps_diff * dQoI << std::endl;
	//	
	//}

	return reference_qoi;
}

std::complex<double> Kriging::get_QoI(Domain & dom_forward, Domain & dom_adjoint)
{
	std::complex<double> qoi = 0.0;
	for (int i = 0; i < dom_forward.scatter1.cGr[1].size(); ++i) {
		qoi += dom_forward.scatter1.cGr[1][i] * std::conj(dom_adjoint.cAlpha[i]);
	}
	return qoi;
}

//monte carlo method for use in generating comparison plot
void Kriging::monte_carlo_instance(std::string & file_name)
{
	//compute the forward solution
	//then compute just the RHS for the adjoint solution (to get the QoI term to compute the monte carlo)
	bool check_q = false;
	bool plot = false;
	bool sens = false;
	bool error = false;
	bool refine = false;
	bool basis_error = false;
	bool higher_order = false;
	bool check_results = false;
	//mesh inputs
	std::string mesh_name = file_name;
	//get the forward solution
	Domain dom_forward(mesh_name);
	run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, false, 0.0);
	//get RHS
	Domain dom_adjoint_rhs(mesh_name);
	run::_RHS_compute(dom_adjoint_rhs, mesh_name, plot, true, higher_order, check_results);
	//get QoI
	std::complex<double> qoi = Kriging::get_QoI2(dom_forward, dom_adjoint_rhs);
}

void Kriging::multi_HOPS_epsr(std::string& file_name)
{


	//load in the list of random materials being tested
	std::vector<std::complex<double>> material_list;

	//old material range
	//std::ifstream materials_in("../ioFiles/input/materials_list.txt");

	//new material range with monte carlo gradient
	std::ifstream materials_in("../ioFiles/input/materials_list4.txt");
	std::string line;
	//std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		//material_list.push_back(std::complex<double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
		material_list.push_back(std::complex<double>(std::stod(data[0]), std::stod(data[1])));
	}
	materials_in.close();

	//for sweeping ref inputs, otherwise use next segment for mat or freq perturbation with set references---------------------------------
	//std::vector<double> referencesReal;
	//std::string refString = "../reference_files_mat/input_ref_list.txt";
	//std::ifstream refs_in(refString);
	//std::string line2;
	////std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	//while (std::getline(refs_in, line2)) {
	//	auto data = functions::split(line2, ' ');
	//	//material_list.push_back(std::complex<double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
	//	referencesReal.push_back(std::stod(data[0]));
	//}
	//refs_in.close();
	//std::vector<std::complex<double>> referencesFull;
	//for (int i = 0; i < referencesReal.size(); i++) {
	//	referencesFull.push_back(std::complex<double>(referencesReal[i], -2.0));
	//}
	//--------------------------------------------------------------------------------------

	//sort the material input list
	sort(material_list.begin(), material_list.end(), complexComparitor);

	//acceptable error
	double delta = 0.0008;

	//convergence condition
	double stdLim = 0.0001;

	//iteration limit
	int itLim = 26;

	//material_list.txt
	//double matStd = 1.0;
	//double matMean = 4.5;

	//material_list2.txt
	double matStd = 1.5;
	double matMean = 6.0;

	std::vector<double> stdDevRCSList;
	std::vector<std::complex<double>> stdDevList;


	///////////////////////Material Perturbation/////////////////////////////////////////////////////////////////////

	//new material range (materials_list3 I think)
	//std::vector<std::complex<double>> referencesFull = { {3.000000, -2.000000},{3.100000, -2.000000},{3.200000, -2.000000},{3.300000, -2.000000},{3.400000, -2.000000},{3.500000, -2.000000},{3.600000, -2.000000},{3.700000, -2.000000},{3.800000, -2.000000},{3.900000, -2.000000},{4.000000, -2.000000},{4.100000, -2.000000},{4.200000, -2.000000},{4.300000, -2.000000},{4.400000, -2.000000},{4.500000, -2.000000},{4.600000, -2.000000},{4.700000, -2.000000},{4.800000, -2.000000},{4.900000, -2.000000},{5.000000, -2.000000},{5.100000, -2.000000},{5.200000, -2.000000},{5.300000, -2.000000},{5.400000, -2.000000},{5.500000, -2.000000},{5.600000, -2.000000},{5.700000, -2.000000},{5.800000, -2.000000},{5.900000, -2.000000},{6.000000, -2.000000},{6.100000, -2.000000},{6.200000, -2.000000},{6.300000, -2.000000},{6.400000, -2.000000},{6.500000, -2.000000},{6.600000, -2.000000},{6.700000, -2.000000},{6.800000, -2.000000},{6.900000, -2.000000},{7.000000, -2.000000},{7.100000, -2.000000},{7.200000, -2.000000},{7.300000, -2.000000},{7.400000, -2.000000},{7.500000, -2.000000},{7.600000, -2.000000},{7.700000, -2.000000},{7.800000, -2.000000},{7.900000, -2.000000},{8.000000, -2.000000},{8.100000, -2.000000},{8.200000, -2.000000},{8.300000, -2.000000},{8.400000, -2.000000},{8.500000, -2.000000},{8.600000, -2.000000},{8.700000, -2.000000},{8.800000, -2.000000},{8.900000, -2.000000},{9.000000, -2.000000},{9.100000, -2.000000},{9.200000, -2.000000},{9.300000, -2.000000},{9.400000, -2.000000},{9.500000, -2.000000},{9.600000, -2.000000},{9.700000, -2.000000},{9.800000, -2.000000},{9.900000, -2.000000},{10.000000, -2.000000} };
	//----------------------------------------------------------------------------------------------------------------

	//material range 4/5
	std::vector<std::complex<double>> referencesFull = { {1.500000, -2.000000}, { 1.600000, -2.000000 }, { 1.700000, -2.000000 }, { 1.800000, -2.000000 }, { 1.900000, -2.000000 }, { 2.000000, -2.000000 }, { 2.100000, -2.000000 }, { 2.200000, -2.000000 }, { 2.300000, -2.000000 }, { 2.400000, -2.000000 }, { 2.500000, -2.000000 }, { 2.600000, -2.000000 }, { 2.700000, -2.000000 }, { 2.800000, -2.000000 }, { 2.900000, -2.000000 }, { 3.000000, -2.000000 }, { 3.100000, -2.000000 }, { 3.200000, -2.000000 }, { 3.300000, -2.000000 }, { 3.400000, -2.000000 }, { 3.500000, -2.000000 }, { 3.600000, -2.000000 }, { 3.700000, -2.000000 }, { 3.800000, -2.000000 }, { 3.900000, -2.000000 }, { 4.000000, -2.000000 }, { 4.100000, -2.000000 }, { 4.200000, -2.000000 }, { 4.300000, -2.000000 }, { 4.400000, -2.000000 }, { 4.500000, -2.000000 }, { 4.600000, -2.000000 }, { 4.700000, -2.000000 }, { 4.800000, -2.000000 }, { 4.900000, -2.000000 }, { 5.000000, -2.000000 }, { 5.100000, -2.000000 }, { 5.200000, -2.000000 }, { 5.300000, -2.000000 }, { 5.400000, -2.000000 }, { 5.500000, -2.000000 }, { 5.600000, -2.000000 }, { 5.700000, -2.000000 }, { 5.800000, -2.000000 }, { 5.900000, -2.000000 }, { 6.000000, -2.000000 }, { 6.100000, -2.000000 }, { 6.200000, -2.000000 }, { 6.300000, -2.000000 }, { 6.400000, -2.000000 }, { 6.500000, -2.000000 }, { 6.600000, -2.000000 }, { 6.700000, -2.000000 }, { 6.800000, -2.000000 }, { 6.900000, -2.000000 }, { 7.000000, -2.000000 }, { 7.100000, -2.000000 }, { 7.200000, -2.000000 }, { 7.300000, -2.000000 }, { 7.400000, -2.000000 }, { 7.500000, -2.000000 }, { 7.600000, -2.000000 }, { 7.700000, -2.000000 }, { 7.800000, -2.000000 }, { 7.900000, -2.000000 }, { 8.000000, -2.000000 }, { 8.100000, -2.000000 }, { 8.200000, -2.000000 }, { 8.300000, -2.000000 }, { 8.400000, -2.000000 }, { 8.500000, -2.000000 }, { 8.600000, -2.000000 }, { 8.700000, -2.000000 }, { 8.800000, -2.000000 }, { 8.900000, -2.000000 }, { 9.000000, -2.000000 }, { 9.100000, -2.000000 }, { 9.200000, -2.000000 }, { 9.300000, -2.000000 }, { 9.400000, -2.000000 }, { 9.500000, -2.000000 }, { 9.600000, -2.000000 }, { 9.700000, -2.000000 }, { 9.800000, -2.000000 }, { 9.900000, -2.000000 }, { 10.000000, -2.000000 }, { 10.100000, -2.000000 }, { 10.200000, -2.000000 }, { 10.300000, -2.000000 }, { 10.400000, -2.000000 }, { 10.500000, -2.000000 }, { 10.600000, -2.000000 }, { 10.700000, -2.000000 }, { 10.800000, -2.000000 }, { 10.900000, -2.000000 }, { 11.000000, -2.000000 } };
	//------------------------------------------------------------------------------------------------

	//size of referencesFull
	double size = referencesFull.size();

	//index of first/middle/last referencesFull to add
	int loRefIndex = find_index(referencesFull, material_list[0]);
	if (loRefIndex < 0) { loRefIndex = 0; }

	int hiRefIndex = find_lowest(referencesFull, material_list[material_list.size() - 1]);
	if (hiRefIndex > size - 1) { hiRefIndex = size - 1; }

	loRefIndex -= 1;

	int midRefIndex = (loRefIndex + hiRefIndex) / 2;

	//initial set of the references to be passed into HOPS (first, middle, last) requires odd number of references in the full vector
	//std::vector<std::complex<double>> references = { referencesFull[0], referencesFull[(size - 1) / 2], referencesFull[size - 1] };
	//references = { referencesFull[5], referencesFull[10], referencesFull[15], referencesFull[20], referencesFull[25], referencesFull[30],referencesFull[35], referencesFull[40], referencesFull[45], referencesFull[50], referencesFull[55], referencesFull[60], referencesFull[65], referencesFull[70] };
	//references = { referencesFull[5], referencesFull[9], referencesFull[13], referencesFull[17], referencesFull[21], referencesFull[25], referencesFull[29], referencesFull[33], referencesFull[37], referencesFull[41], referencesFull[45], referencesFull[49], referencesFull[53], referencesFull[57], referencesFull[61], referencesFull[65], referencesFull[70] };
	//references = referencesFull;

	//the current references being used for Kriging
	std::vector<std::complex<double>> references;
	std::vector<int> refIndex;
	//references = { referencesFull[loRefIndex], referencesFull[midRefIndex], referencesFull[hiRefIndex] };
	//refIndex = { loRefIndex, midRefIndex, hiRefIndex };

	////easy sweep condition
	refIndex.clear();
	references.clear();
	//list [43, 32, 24, 19, 16, 14, 12, 11, 9, ]
	for (int i = 0; i < 95; i += 48) {
		refIndex.push_back(i);
		references.push_back(referencesFull[i]);
	}
	if (refIndex[refIndex.size() - 1] != 95) {
		refIndex.push_back(95);
		references.push_back(referencesFull[95]);
	}
	////-------
	//references.push_back(referencesFull[25]);
	//references.push_back(referencesFull[70]);
	//refIndex.push_back(25);
	//refIndex.push_back(70);


	//qoi_list and HOPS_splitting were from HOPS UQ
	std::vector<std::vector<std::complex<double>>> qoi_list(references.size());
	std::vector<std::vector<std::complex<double>>> HOPS_splitting(references.size());

	//refIndex is the index of the ref_i file so that it can be called correctly from the reference_files folder
	//std::vector<int> refIndex = { 0, int(size - 1) / 2, int(size - 1) };
	//refIndex = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70 };
	//refIndex = { 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 70 };


	//toggle to end the while loop
	bool toggle = true;
	bool whileCond = true;
	//some declarations
	std::string in_file_name;

	int iterations = 0;


	std::vector<std::vector<double>> RCS(references.size());
	std::complex<double> sum = 0.0;
	double sumRCS = 0.0;

	std::complex<double> mean;
	double meanRCS;

	std::complex<double> stdDev = 0.0;
	double stdDevRCS = 0.0;
	//std::vector<double> pdfScaleFactor(references.size());

	//initiate vectors
	std::vector<std::vector<double>> weights(material_list.size());
	std::vector<std::vector<double>> weightsGrad(material_list.size());
	std::vector<std::vector<double>> weightsHigher(material_list.size());
	std::vector<std::vector<std::complex<double>>>(material_list.size());
	std::vector<std::vector<double>> weightsHigherReal(material_list.size());
	std::vector<std::vector<double>> weightsHigherImag(material_list.size());
	std::vector<std::vector<std::complex<double>>> weightsHigherComplex(material_list.size());

	std::vector<std::complex<double>> reconstruction(material_list.size());
	std::vector<std::complex<double>> reconstructionReal(material_list.size());
	std::vector<std::complex<double>> reconstructionImag(material_list.size());
	std::vector<std::complex<double>> reconstructionHigher(material_list.size());
	std::vector<std::complex<double>> reconstructionUpdate(material_list.size());
	std::vector<std::complex<double>> reconstructionGrad(material_list.size());
	//finding weights for inserted reference value
	//current insertion is garbage code
	int bestIndex;
	double bestAvg = 0.0;
	double bestVar = 0.0;
	int insertIndex;

	//std::complex<double> bestReference;
	//std::complex<double> bestReferenceVal;
	std::vector<std::complex<double>> referencesTempBest;
	std::vector<std::complex<double>> referenceValsTempBest;
	std::vector<std::vector<double>> weightsTemp(material_list.size());
	std::vector<double> referencesReal(references.size());
	std::vector<std::complex<double>> referencesTemp(references.size());
	std::vector<std::complex<double>> referenceVals(references.size());
	std::vector<std::complex<double>> referenceValsReal(references.size());
	std::vector<std::complex<double>> referenceValsImag(references.size());

	std::vector < std::complex<double>> gradient(references.size());

	bool maxVar = true;
	bool avgVar = true;
	double varAvg;


	//while loop continues to adaptively refine until given convergence condition is met
	while (whileCond) {
		std::cout << "_____________________________________________________________Pass in the for loop_______________________________________________\n";
		if (toggle == false) {
			whileCond = false;
		}

		//need to get rid of qoi_list (not neccessary for kriging), but I need to first get rid of its use in sensitivity_to_epsr method.
		
		//reset vectors for adaptive refinement
		qoi_list.clear();
		qoi_list.resize(references.size());
		weights.clear();
		weights.resize(material_list.size());
		weightsGrad.clear();
		weightsGrad.resize(material_list.size());
		weightsHigher.clear();
		weightsHigher.resize(material_list.size());
		reconstruction.clear();
		reconstruction.resize(material_list.size());
		//referenceVals.clear();
		//referenceVals.resize(references.size());

		//new set of references (should include same plus one new)
		for (int i = 0; i < references.size(); ++i) {
			std::cout << "references[i]: " << references[i] << " ref Index[i]: " << refIndex[i] << std::endl;
		}

		std::cout << "ref size: " << references.size() << std::endl;

		//iterating through references now and get the qoi and gradients
		for (int i = 0; i < references.size(); ++i) {
			//if (HOPS_splitting[i].size() == 0) continue;
			in_file_name = file_name;



			//setting so that new passes in the while loop only run solve for new points---------------------------------------------
			//only works if starting with 3 points
			// adaptive refinement method
			//if (references.size() == 3) {			
			//	std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange4\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//	std::system(command.c_str());
			//	std::cout << "Testing for reference number: " << refIndex[i] << std::endl;

			//	referenceVals[i] = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i], gradient[i]);
			//	referenceValsReal[i] = referenceVals[i].real();
			//	referenceValsImag[i] = referenceVals[i].imag();
			//}
			//else {
			//	std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange4\\ref_" + std::to_string(refIndex[insertIndex]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//	std::system(command.c_str());
			//	std::cout << "Testing for reference number: " << refIndex[insertIndex] << std::endl;

			//	std::complex<double> val;

			//	//insert a zero before calling sensitivity so that gradient[i] is in the correct place
			//	insert_complex(gradient, insertIndex, 0.0);

			//	val = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[insertIndex], qoi_list[insertIndex], references[insertIndex], gradient[insertIndex]);
			//	insert_complex(referenceVals, insertIndex, val);
			//	insert_complex(referenceValsReal, insertIndex, val.real());
			//	insert_complex(referenceValsImag, insertIndex, val.imag());

			//	
			//	break;
			//}
			//--------------------------------------------------------------------------------------------------------------
			

			//this is for even split--------------------------------------------------------------------------
			//will solve all references again, not just new point
			std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange4\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			std::system(command.c_str());
			std::cout << "Testing for reference number: " << refIndex[i] << std::endl;

			referenceVals[i] = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i], gradient[i]);
			//-----------------------------------------------------------------------


		}

		//adding references point based on error (old method)
		double loError = 0.0;
		double hiError = 0.0;

		std::vector<double> loErrorVec;
		std::vector<double> hiErrorVec;

		//subdivide highest error
		int hiErrorIndex = 0;
		double hiErrorCheck = 0.0;


		//val is the endpoint of HOPS_splitting
		int val = HOPS_splitting.size() - 1;
		std::vector<std::complex<double>> referencesNew = references;
		std::vector<int> refIndexNew = refIndex;
		
		//calculate variogram (empirical method)-----------------------------------------

		//main params for variogram (this and fit method)
		double nBins = 10;
		double maxDist = abs(references[references.size() - 1].real() - references[0].real());// / 2.0;
		//maxDist = 10.0;
		//std::cout << "check 3------------------------\n";
		//double binTol = maxDist / nBins;
		std::vector<double> edges;
		std::vector<double> edgesGrad;
		std::vector<double> edgesComplex;
		std::vector<double> edgesGradComplex;
		std::vector<double> edgesValGrad;
		std::vector<double> edgesGradVal;

		for (int i = 0; i <= nBins; i++) {
			edges.push_back(i * maxDist / nBins);
			edgesComplex.push_back(i * maxDist / nBins);
			edgesGradComplex.push_back(i * maxDist / nBins);
			edgesValGrad.push_back(i * maxDist / nBins);
			edgesGrad.push_back(i * maxDist / nBins);
			edgesGradVal.push_back(i * maxDist / nBins);

		}

		//linspace for the x values in the variogramFit
		int fitSize = 15000;
		std::vector<double> d;
		for (int i = 0; i <= fitSize; ++i) {
			d.push_back(i * maxDist / fitSize);
		}

		//build matrix of distance values
		std::vector<std::vector<std::complex<double>>> dMat = buildDMatrix2(references, referenceVals, referenceVals);
		std::vector<std::vector<std::complex<double>>> dMatGrad = buildDMatrix2(references, gradient, gradient);
		std::vector<std::vector<std::complex<double>>> dMatGradVal = buildDMatrix2(references, gradient, referenceVals); // these are approx 0
		std::vector<std::vector<std::complex<double>>> dMatValGrad = buildDMatrix2(references, referenceVals, gradient); //

		//build variogram matrix and solve for empirical variogram
		std::vector<double> variogram = buildCoVariogram(dMat, references, edges);
		std::vector<double> variogramGrad = buildCoVariogram(dMatGrad, references, edges);
		std::vector<double> variogramGradVal = buildCoVariogram(dMatGradVal, references, edges);
		std::vector<double> variogramValGrad = buildCoVariogram(dMatValGrad, references, edges);

		std::vector<std::complex<double>> variogramComplex = buildCoVariogramComplex(dMat, references, edges);
		//std::vector<std::complex<double>> variogramGradComplex = buildVariogramComplex(dMatGrad, references, edges);

		//trim any nan values
		trimVariogram(variogram, edges);
		trimVariogram(variogramGrad, edgesGrad);
		trimVariogram(variogramGradVal, edgesGradVal);
		trimVariogram(variogramValGrad, edgesValGrad);

		trimVariogramComplex(variogramComplex, edgesComplex);
		//trimVariogramComplex(variogramGradComplex, edgesGradComplex);

		//fit model to empirical variogram
		std::vector<double> variogramFit = fitCoVariogramModel(variogram, edges, d, references.size());
		std::vector<double> variogramFitGrad = fitCoVariogramModel(variogramGrad, edgesGrad, d, references.size());
		std::vector<double> variogramFitGradVal = fitCoVariogramModel(variogramGradVal, edgesGradVal, d, references.size());
		std::vector<double> variogramFitValGrad = fitCoVariogramModel(variogramValGrad, edgesValGrad, d, references.size());

		std::vector<double> variogramComplexReal(variogramComplex.size());
		std::vector<double> variogramComplexImag(variogramComplex.size());
		//std::vector<double> variogramGradComplexReal(variogramGradComplex.size());
		//std::vector<double> variogramGradComplexImag(variogramGradComplex.size());

		for (int i = 0; i < variogramComplex.size(); ++i) {
			std::cout << variogramComplex[i].real() << "  " << variogramComplex[i].imag() << "\n";
			variogramComplexReal[i] = variogramComplex[i].real();
			variogramComplexImag[i] = variogramComplex[i].imag();
		}

		std::vector<double> variogramFitReal = fitCoVariogramModel(variogramComplexReal, edgesComplex, d, references.size());
		std::vector<double> variogramFitImag = fitCoVariogramModel(variogramComplexImag, edgesComplex, d, references.size());
		//std::vector<double> variogramFitGradReal = fitVariogramModel(variogramGradComplexReal, edges, d, references.size());
		//std::vector<double> variogramFitGradImag = fitVariogramModel(variogramGradComplexImag, edges, d, references.size());
		// end of empirical variogram methods ----------------------------------------------------------------------------------



		// working MLE ----------------------------------------------------------------------------------------------------------------

		//theta is the hyperparameter
		double theta;
	

		//input data, 1 since dimensionality is currently 1
		MatrixXd X(referenceVals.size(), 1);
		VectorXd y(referenceVals.size());
		VectorXd yGrad(referenceVals.size());
		MatrixXd R(referenceVals.size(), referenceVals.size());
		VectorXd absRefVals(referenceVals.size());

		for (int k = 0; k < referenceVals.size(); k++) {
			X(k, 0) = references[k].real();
			y(k) = abs(referenceVals[k]);
			//y(k + referenceVals.size()) = abs(gradient[k]);
			yGrad(k) = abs(gradient[k]);
			absRefVals(k) = abs(referenceVals[k]);
		}

		int n = referenceVals.size();
		int m = n;

		//calculate basis function matrix
		MatrixXd B = MatrixXd::Zero(n, m);
		basisCalc(X, y, B, yGrad);

		std::cout << "basis calc works\n";

		MLE_estimation(X, y, theta, B, yGrad);

		std::cout << "MLE successful theta:" << theta << "\n";


		std::cout << "aCalc successful\n";

		//calculate variance estimate
		double sigma2;
		var_calc(X, y, sigma2, theta, B, yGrad);

		//currently not needed when done this way
		////calculate basis function coefficient a
		////size(a) = m
		//VectorXd a(referenceVals.size());
		//aCalc(X, y, B, theta, a);
		///////////////////////////////////


		//theta = theta / 4.0;
		std::cout << "values for variance: " << sigma2 << "  value for theta: " << theta << std::endl;

		for (int k = 0; k < d.size(); k++) {
			variogramFit[k] = covariance(d[k], theta, sigma2);
		}

		// end of MLE -----------------------------------------------------------------------------

		//calculate variogram derivatives
		std::vector<double> variogramFitFirstD = firstDerivative(variogramFit, d);
		std::vector<double> variogramFitSecondD = firstDerivative(variogramFitFirstD, d);



		//------------------------------- output empirical variogram and fitted variogram for checking----------------------------
		
		std::string outvar = "../iofiles/output/paper/variogram/fit" + std::to_string(iterations + 3) + ".txt";
		std::ofstream var_out(outvar);
		for (int i = 0; i < d.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			var_out << d[i] << "\t" << variogramFit[i] << std::endl;
		}
		var_out.close();

		std::string outvarfit = "../iofiles/output/paper/variogram/fitFirstD" + std::to_string(iterations + 3) + ".txt";
		std::ofstream var_fit_out(outvarfit);
		for (int i = 0; i < d.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			var_fit_out << d[i] << "\t" << variogramFitFirstD[i] << std::endl;
		}
		var_fit_out.close();

		std::string outvarfit2 = "../iofiles/output/paper/variogram/fitSecondD" + std::to_string(iterations + 3) + ".txt";
		std::ofstream var_fit2_out(outvarfit2);
		for (int i = 0; i < d.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			var_fit2_out << d[i] << "\t" << variogramFitSecondD[i] << std::endl;
		}
		var_fit2_out.close();

		///----------------------------------------------------------------------------------------

		varAvg = 0.0;
		//---------------------  Kriging System Solution --------------------------------------------
		double xAvg;
		double xAvgReal;
		double xAvgImag;

		std::cout << "made it to reconstruction\n";


		for (int i = 0; i < material_list.size(); ++i) {

			//currently being used to fit gradient values////////////
			taylorKriging(variogramFitGrad, d, references, weightsGrad[i], material_list[i].real());
			reconstructionGrad[i] = krigSum(weightsGrad[i], gradient);
			////////////////////////////////////////////////
			
			/////////////GETK///////////////
			//taylorKrigingHigherOrder(variogramFit, variogramFitFirstD, variogramFitSecondD, d, references, weightsHigher[i], material_list[i].real(), theta);
			//reconstruction[i] = krigSumHigherOrder(weightsHigher[i], referenceVals, gradient, material_list[i].real(), xAvg);
			reconstruction[i] = taylorKrigingHigherOrder3(variogramFit, variogramFitFirstD, variogramFitSecondD, d, references, weightsHigher[i], material_list[i].real(), variogramFitGrad, referenceVals, gradient, theta);
			/////////////////////////////////////////////
			
			///////////TK//////////////////////
			//taylorKriging(variogramFit, d, references, weights[i], material_list[i].real());
			//reconstruction[i] = krigSum(weights[i], referenceVals);
			////////////////////////////////			
		}

		bestVar = 0;

		std::cout << "made it out of reconstruction\n";


		double matMean = 4.5;
		double matStd = 1.0;
		//------------------------------------------------------------------------------------------------------------------
		
		//find max variance for next iteration--------------------------------------------------------------------------------------------------
		for (int i = 0; i < referencesFull.size(); ++i) {
			double variance = 0.0;
			double avg = 0;
			double dist = 0;

			//skip values of referencesFull that are already contained in references
			if (std::find(references.begin(), references.end(), referencesFull[i]) != references.end()) {
				std::cout << "references[" << i << "] found \n";
				continue;
			}

			////skip values of referencesFull outside the range
			// not currently using RefIndex
			//if (i < loRefIndex || i >  hiRefIndex)
			//	continue;

			//resize temporary vectors
			weightsTemp.clear();
			referencesTemp.clear();
			weightsTemp.resize(material_list.size());
			referencesTemp.resize(references.size());
			referencesTemp = references;
			//referencesReal.resize(references.size());

			//find the index of insertion point
			int index = find_lowest(references, referencesFull[i]);

			//insert referenceFull value into referencesTemp
			insert_complex(referencesTemp, index, referencesFull[i]);

			//get index of the added value
			int weightIndex = find_index(referencesTemp, referencesFull[i]);


			std::vector<std::complex<double>> referenceValsTemp(referencesTemp.size());

			//REDO MLE WITH NEW VECTOR

			//theta is the hyperparameter
			double theta;


			//input data, 1 since dimensionality is currently 1
			MatrixXd X(referenceVals.size(), 1);

			//yData (QoI)
			VectorXd y(referenceVals.size());

			//gradient data
			VectorXd yGrad(referenceVals.size());

			//correlation matrix
			MatrixXd R(referenceVals.size(), referenceVals.size());

			//absolute values of reference data
			VectorXd absRefVals(referenceVals.size());

			for (int k = 0; k < referenceVals.size(); k++) {
				X(k, 0) = references[k].real();
				y(k) = abs(referenceVals[k]);
				//y(k + referenceVals.size()) = abs(gradient[k]);
				yGrad(k) = abs(gradient[k]);
				absRefVals(k) = abs(referenceVals[k]);
			}

			int n = referenceVals.size();
			int m = n;

			//calculate basis function matrix
			MatrixXd B = MatrixXd::Zero(n, m);
			basisCalc(X, y, B, yGrad);

			std::cout << "basis calc works\n";

			MLE_estimation(X, y, theta, B, yGrad);

			std::cout << "MLE successful theta:" << theta << "\n";


			std::cout << "aCalc successful\n";

			//calculate variance estimate
			double sigma2;
			var_calc(X, y, sigma2, theta, B, yGrad);

			if (variance > bestVar){

				bestIndex = i;
				//bestAvg = avg;
				bestVar = variance;
				//bestReference = referencesFull[i];
				//std::cout << "best value---------------------------------------------------: " << referencesFull[i] << std::endl;
				//bestReferenceVal = referenceValsTemp[weightIndex];
				
				referencesTempBest.clear();
				referencesTempBest.resize(referencesTemp.size());
				referenceValsTempBest.clear();
				referenceValsTempBest.resize(referenceValsTemp.size());
				referencesTempBest = referencesTemp;

				//std::cout << "referencesValsTempBest[0] Assigned: " << referenceValsTempBest[0] << std::endl;

				referenceValsTempBest = referenceValsTemp;
			}
			
			//------------------------------------------------
		}

	// now just need to add the cooresponding index into references and recompute the reconstruction
		//insert referencesFull[bestIndex] into references.
		//recompute
		//auto it = references.begin();
		// 
		//set new reference values;

		//average value of sample points

		//variance = variance / integral;

		//std::cout << "referencesValsTempBest[0] outside loop: " << referenceValsTempBest[0] << std::endl;

		double z0 = 0.0;
		for (int i = 0; i < referenceValsTempBest.size(); ++i) {
			z0 += referenceValsTempBest[i].real();
		}

		std::cout << "z0 before division: " << z0 << "  divide factor: " << referenceValsTempBest.size() << std::endl;

		z0 /= referenceValsTempBest.size();



		std::cout << "bestVar: " << bestVar << "   varAvg: " << varAvg << "   z0: " << z0 << std::endl;
		std::cout << " max variance condition: " <<  1 / z0 * bestVar << std::endl;
		std::cout << "average variance condition: "  << (1 / z0 * sqrt(varAvg / (referencesFull.size() - references.size()))) << std::endl;
		//stopping criteria
		/*if ((1 / z0 * sqrt(varAvg / (referencesFull.size()-references.size()))) < 0.00015 && (1 / z0 * bestVar) < 0.00025) {*/
		if ((1 / z0 * sqrt(varAvg / (referencesFull.size() - references.size()))) < 0.00011 && (1 / z0 * bestVar) < 0.00021) {

			//toggle = false;
			std::cout << "convergence condition met\n";
			//break;
		}

		auto refIt = std::lower_bound(refIndex.begin(), refIndex.end(), bestIndex);
		insertIndex = refIt - refIndex.begin();
		refIndex.insert(refIt, bestIndex);

		references.clear();
		references.resize(referencesTempBest.size());

		for (int i = 0; i < referencesTempBest.size(); ++i) {
			std::cout << "referencesTempBest[i]: " << referencesTempBest[i] << "  refIndx[i]: " << refIndex[i] << std::endl;
		}



		references = referencesTempBest;

		for (int i = 0; i < references.size(); ++i) {
			std::cout << "references[i]: " << references[i] << std::endl;
		}
	

		////////////done to subdivide lowest error////////////////////
		//auto it = references.begin();
		//auto it2 = refIndex.begin();
		//int refFullIndex = (refIndex[hiErrorIndex + 1] - refIndex[hiErrorIndex]) / 2 + refIndex[hiErrorIndex];
		//references.insert(it + hiErrorIndex + 1, referencesFull[refFullIndex]);
		//refIndex.insert(it2 + hiErrorIndex + 1, refFullIndex);


		//iteration limit//
		iterations++;
		//if (iterations >= itLim) {
		//	std::cout << "iteration condition met\n";
		//	break;
		//}
		///-------------------------------------------------------------

		//AR sweep output
		std::string outF = "../ioFiles/output/paper/mat4/GETK_gauss_MLE_nullBasis/qoi_kriging" + std::to_string(iterations) + ".txt";
		std::ofstream qoi_dist_out(outF);
		for (int i = 0; i < material_list.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			if (i == 0)
				std::cout << std::setprecision(12) << "precision check" << reconstruction[i].real() << std::endl;
			qoi_dist_out << std::setprecision(12) << material_list[i].real() << " " << material_list[i].imag() << " " << reconstruction[i].real() << " " << reconstruction[i].imag() << " " << reconstructionGrad[i].real() << " " << reconstructionGrad[i].imag() << std::endl;
		}
		qoi_dist_out.close();
		//qoi_dist_out.close();

	}

	std::cout << "out of the while loop=============================\n";



	////output one time
	////output results to file
	//std::ofstream qoi_dist_out("../ioFiles/output/kriging/qoi_kriging_14.txt");
	//for (int i = 0; i < material_list.size(); ++i) {
	//	if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
	//		continue;
	//	}
	//	qoi_dist_out << material_list[i].real() << " " << material_list[i].imag() << " " << reconstruction[i].real() << " " << reconstruction[i].imag() << std::endl;
	//}
	//qoi_dist_out.close();
}

double Kriging::normalFunction(double x, double mean, double std) {
	double arg = (x - mean) / std;
	double exponential = exp(-0.5 * arg * arg);
	double integral = 1.0 / (std * sqrt(2.0 * 3.14159)) * exponential;
	return integral;
}

double Kriging::trap_integral(double intervalMid, double mean, double std, double step) {
	double count = 50.0;
	double intervalBeg = intervalMid - count * step / 2.0;
	double intervalEnd = intervalMid + count * step / 2.0;
	double integral = 0.5 * (Kriging::normalFunction(intervalBeg, mean, std) + Kriging::normalFunction(intervalEnd, mean, std));
	for (int i = 0; i < (int)count; ++i) {
		integral += Kriging::normalFunction(intervalBeg + step * i, mean, std);
	}
	integral *= step;
	return integral;

}
/*
void Kriging::multi_HOPS_multi_epsr(std::string & file_name) //mutli HOPS for variation of multiple material values
{
	//load in the list of random materials being tested
	std::vector<std::vector<std::complex<double>>> material_list;
	//std::ifstream materials_in("../monte_carlo_files_highdim/materials_list.txt");
	std::ifstream materials_in("../ioFiles/input/frequencies_list1.txt");
	std::string line;
	

	//get the monte carlo points
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		std::vector<std::complex<double>> material_list_temp;
		for (int mati = 0; mati < data.size(); mati+=1) {
			// b/c the data is irrelevant for air and PML layer, set those to zero
			if (mati < 64)
				material_list_temp.push_back(std::complex<double>(std::stod(data[mati]),0.0));
			else
				material_list_temp.push_back(0.0);
		}
		material_list.push_back(material_list_temp);
	}
	materials_in.close();
	
	////generate references in c++
	//double sample_start = 2.56 - .21; //3 std
	//double sample_end = 2.56 + .21;//3 std
	//int num_refs = 1;
	//std::vector<std::vector<std::complex<double>>> references(num_refs);
	//double increment = (sample_end - sample_start) / (double(num_refs-1.0));
	//for (int refs = 0; refs < num_refs; ++refs) {
	//	double val = sample_start + refs*increment;
	//	if (num_refs == 1)
	//		val = 2.56;
	//	for (int numel = 0; numel < 256; ++numel) {
	//		if (numel < 64)
	//			references[refs].push_back(val);
	//		else
	//			references[refs].push_back(0.0);
	//	}
	//}
	
	//generate using FILE IO
	std::vector<std::vector<std::complex<double>>> references;
	std::ifstream references_in("../monte_carlo_files_highdim/references_highdim.txt");
	while (std::getline(references_in, line)) {
		auto data = functions::split(line, ' ');
		std::vector<std::complex<double>> material_list_temp;
		for (int mati = 0; mati < data.size(); mati += 1) {
			// b/c the data is irrelevant for air and PML layer, set those to zero
			if (mati < 64)
				material_list_temp.push_back(std::complex<double>(std::stod(data[mati]), 0.0));
			else
				material_list_temp.push_back(0.0);
		}
		if (data.size() < 256) {
			for (int mati = 64; mati < 256; ++mati) {
				material_list_temp.push_back(0.0);
			}
		}
		references.push_back(material_list_temp);
	}
	references_in.close();

	//cycle through and give the appropriate ones to the number of hops references

	
	
	std::vector<std::vector<std::vector<std::complex<double>>>> HOPS_splitting(references.size());
	
	int mat_counter = 0;
	

	for (int j = 0; j < material_list.size(); ++j) {
		//find the closest reference point
		double min_dist = 1.0e15;
		int mindex = 0;
		for (int i = 0; i < references.size(); ++i) {
			
				double dist = norm(references[i] - material_list[j]);
				if (dist < min_dist) { 
					mindex = i; 
					min_dist = dist;
				}
			
		}
		
			HOPS_splitting[mindex].push_back(material_list[j]);
			mat_counter++;
			
		
	}
	std::vector<std::vector<std::complex<double>>> qoi_list(references.size());
	//iterating through references now and get the qois
	for (int i = 0; i < references.size(); ++i) {
		if (HOPS_splitting[i].size() == 0) continue;
		//if (i != 2) continue; //only testing the mean reference right now
		//load the correct input file into the references folder
		std::string command = "..\\file_input_converter ..\\reference_files_2im\\ref_0.in ..\\exampleFiles\\" + file_name + "\\";
		std::system(command.c_str());
		std::cout << "Testing for reference number: " << i << std::endl;
		Kriging::sensitivity_to_multi_epsr(file_name, HOPS_splitting[i], qoi_list[i], references[i]);
	}
	//output results to file
	//std::ofstream qoi_dist_out("../exampleFiles/debug_data/" + file_name + "/results/qoi_HOPS_dist.txt");

	//test
	std::ofstream qoi_dist_out("../ioFiles/output/qoi_HOPS_dist1.txt");
	qoi_dist_out << std::setprecision(16);
	std::cout << qoi_list.size() << std::endl;
	for (int i = 0; i < qoi_list.size(); i++) {
		for (int j = 0; j < qoi_list[i].size(); ++j) {
			//std::cout << "i: " << i << " j: " << j << std::endl;
			qoi_dist_out << qoi_list[i][j].real() << " " << qoi_list[i][j].imag() << std::endl;
		}
	}
	qoi_dist_out.close();
}
*/

/*
void Kriging::sensitivity_to_multi_epsr(std::string & file_name, std::vector<std::vector<std::complex<double>>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::vector<std::complex<double>>& reference_values)
{
	//get orig QoI
	//get entries of Bint
	//get part of Gint
	bool check_q = false;
	bool plot = false;
	bool sens = false;
	bool error = false;
	bool refine = false;
	bool basis_error = false;
	bool useAdjoint = false;
	bool higher_order = false;
	bool check_results = false;
	//mesh inputs
	std::string mesh_name = file_name;
	//get the forward solution
	auto t1 = std::chrono::high_resolution_clock::now();
	Domain dom_forward(mesh_name);
	run::_standard_custom_materials(dom_forward, mesh_name, plot, false, higher_order, check_results, false, reference_values);
	std::cout << "Forward solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count() << std::endl;
	//get the adjoint solution
	Domain dom_adjoint(mesh_name);
	run::_standard_custom_materials(dom_adjoint, mesh_name, plot, true, higher_order, check_results, false, reference_values);
	//get the reference qoi
	std::complex<double> reference_qoi = get_QoI(dom_forward, dom_adjoint);
	std::cout << "Reference qoi: " << reference_qoi << std::endl;

	double k0 = dom_forward.scatter1.K0[1];
	//get the remainder of the estimate
	std::vector<std::complex<double>> dQoI(epsr_list.size(), 0.0);
	
	
	//loop through the stored Bint pieces
	for (auto e = dom_forward.elements.begin(); e != dom_forward.elements.end(); ++e) {

		if (e->materials.region == 1) {
			//if in scatterer
			for (int i_unkn = e->unknownsStart; i_unkn <= e->unknownsEnd; ++i_unkn) {
				int icon = abs(dom_forward.vectorD[i_unkn]);
				for (int epsr_val = 0; epsr_val < epsr_list.size(); ++epsr_val) {
					std::complex<double> eps_diff = epsr_list[epsr_val][e->index-1] / e->materials.epsr_list[0][0] - 1.0;
					dQoI[epsr_val] += eps_diff*e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
				}
				//	dQoI1 += eps_diff * e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
				for (int j_unkn = e->unknownsStart; j_unkn <= e->unknownsEnd; ++j_unkn) {
					int jcon = abs(dom_adjoint.vectorD[j_unkn]);
					int i0 = i_unkn - e->unknownsStart + 1;
					int j0 = j_unkn - e->unknownsStart + 1;
					
					for (int epsr_val = 0; epsr_val < epsr_list.size(); ++epsr_val) {
						std::complex<double> eps_diff = epsr_list[epsr_val][e->index - 1] / e->materials.epsr_list[0][0] - 1.0;
						dQoI[epsr_val] += eps_diff * e->cBr_HOPS(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]);
					}
				}
			}
		}
	}
	
	for (int i = 0; i < epsr_list.size(); ++i) {
		//eps_diff = epsr_list[i] / reference_material - 1.0;
		updated_qoi.push_back(reference_qoi + dQoI[i]);
		std::cout << reference_qoi + dQoI[i] << std::endl;
		
	}
	
}
*/

