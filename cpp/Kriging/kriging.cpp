#include "kriging.h"
#include <math.h>
//#include "..\HOPS\hops.cpp"
#include <Eigen/Dense>
#include "src/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include "src/interpolation.h"


using namespace Eigen;

static long double gaussModel(long double h, long double r, long double c0, long double b) {
	long double y;
	long double a = r / 2.0;
	y = b + c0 * (1 - exp(-(h * h) / (a * a)));
	return y;
}


std::complex<long double> Kriging::get_QoI2(Domain& dom_forward, Domain& dom_adjoint) {
	std::complex<long double> qoi = 0.0;
	for (int i = 0; i < dom_forward.cAlpha.size(); ++i) {
		qoi += dom_forward.cAlpha[i] * std::conj(dom_adjoint.scatter1.cGr[1][i]);
	}
	return qoi;
}

static inline std::complex<long double> operator/(std::complex<long double> c1, long double d1) {
	std::complex<long double> output = { c1.real() / d1, c1.imag() / d1 };
	return output;
}


//currently using for real part ie material parameters
static int find_lowest(std::vector<std::complex<long double>> array, std::complex<long double> value) {
	long double min = 1e15;
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

static int find_index(std::vector<std::complex<long double>> array, std::complex<long double> value) {
	long double min = 1e15;
	int index;
	for (int i = 0; i < array.size(); i++) {
		if (abs(value.real() - array[i].real()) < min) {
			min = abs(array[i].real() - value.real());
			index = i;
		}
	}
	return index;

}

static bool complexComparitor(std::complex<long double> a, std::complex<long double> b) {
	return real(a) <= real(b);
}

static void insert_complex(std::vector<std::complex<long double>>& array, int index, std::complex<long double> value) {
	std::vector<std::complex<long double>> tempArray(array.size() + 1);

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

static long double norm(std::vector < std::complex<long double>> v1) {
	std::complex<long double> norm = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		norm += v1[i] * std::conj(v1[i]);
	}
	long double normd = abs(norm);
	return sqrt(normd);

}

//std::vector<std::complex<long double>> vectorInsert(int position, std::complex<long double> value, std::vector<std::complex<long double>> vector) {
//	std::vector<std::complex<long double>> newVector(vector.size() + 1);
//
//	for (int i = 0; i < newVector.size(); ++i) {
//		while (i < position) {
//			newVector[i] = vector[i];
//		}
//		if (i == position) {
//			newVector[i] = value;
//		}
//		while (i > position) {
//			newVector[i] = vector[i - 1];
//		}
//	}
//
//	return newVector;
//
//}

static inline std::vector<std::complex<long double>> operator-(std::vector < std::complex<long double>>& v1, std::vector < std::complex<long double>>& v2) {
	std::vector<std::complex<long double>> vout;
	for (int i = 0; i < v1.size(); ++i) {
		vout.push_back(v1[i] - v2[i]);
	}
	return vout;
}

static std::vector<std::vector<std::complex<long double>>> buildDMatrix(std::vector<std::complex<long double>> x, std::vector<std::complex<long double>> y) {
	std::vector<std::vector<std::complex<long double>>> dMat(x.size() * x.size());
	int index = 0;
	for (int i = 0; i < x.size(); ++i) {
		long double RCSi = sqrt(y[i].real() * y[i].real() + y[i].imag() * y[i].imag()) / 4 / 3.14159;
		for (int j = 0; j < x.size(); ++j) {
			long double RCSj = sqrt(y[j].real() * y[j].real() + y[j].imag() * y[j].imag()) / 4 / 3.14159;
			
			//std::cout << "index: " << index << " x.size(): " << x.size() << std::endl;
			dMat[index].push_back(double(i));
			dMat[index].push_back(double(j));
			dMat[index].push_back(abs(x[i].real() - x[j].real()));

			//will square the term when building the variogram instead of here
			dMat[index].push_back(abs(RCSi - RCSj));// *(RCSi - RCSj));
			dMat[index].push_back(RCSi);
			dMat[index].push_back(RCSj);
			//dMat[index].push_back(abs(y[i] - y[j]));
			index++;
		}
	}
	return dMat;
}

static std::vector<std::vector<std::complex<long double>>> buildDMatrix2(std::vector<std::complex<long double>> x, std::vector<std::complex<long double>> yI, std::vector<std::complex<long double>> yJ) {
	std::vector<std::vector<std::complex<long double>>> dMat(x.size() * x.size());
	int index = 0;
	for (int i = 0; i < x.size(); ++i) {
		long double RCSi = sqrt(yI[i].real() * yI[i].real() + yI[i].imag() * yI[i].imag()) / 4 / 3.14159;
		for (int j = 0; j < x.size(); ++j) {
			long double RCSj = sqrt(yJ[j].real() * yJ[j].real() + yJ[j].imag() * yJ[j].imag()) / 4 / 3.14159;

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

static std::vector<long double> firstDerivative(std::vector<long double> function, std::vector<long double> distance) {
	std::vector<long double> derivative;
	for (int i = 0; i < function.size(); ++i) {

		if (i == 0) {
			derivative.push_back(function[i + 1] - function[i] / (distance[i + 1] - distance[i]));
		}
		else {
			derivative.push_back(function[i] - function[i - 1] / (distance[i] - distance[i - 1]));
		}

		//if (i > 2) {
		//	derivative.push_back(function[i + 1] - function[i] / (distance[i + 1] - distance[i]));
		//}
		//else if (i > function.size() - 3) {
		//	derivative.push_back(function[i] - function[i - 1] / (distance[i] - distance[i - 1]));
		//}
		//else
		//	derivative.push_back((function[i - 2] - 8.0 * function[i - 1] + 8.0 * function[i + 1] - function[i + 2]) / 12 * (distance[i] - distance[i-1]));
	}
	return derivative;
}

static std::vector<long double> buildVariogram(std::vector<std::vector<std::complex<long double>>> dMat, std::vector<std::complex<long double>> x, std::vector<long double>& edges) {

	std::vector<long double> binCounts(edges.size()-1);
	std::vector<long double> binIdx(x.size()*x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	std::vector<long double> variogram(edges.size() - 1);
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

static std::vector<std::complex<long double>> buildVariogramComplex(std::vector<std::vector<std::complex<long double>>> dMat, std::vector<std::complex<long double>> x, std::vector<long double>& edges) {

	std::vector<long double> binCounts(edges.size() - 1);
	std::vector<long double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	std::vector<std::complex<long double>> variogram(edges.size() - 1);
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

static std::vector<long double> buildCoVariogram(std::vector<std::vector<std::complex<long double>>> dMat, std::vector<std::complex<long double>> x, std::vector<long double>& edges) {

	std::complex<long double> sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += dMat[i][4];
	}
	std::complex<long double> avg = sum / x.size();

	std::vector<long double> binCounts(edges.size() - 1);
	std::vector<long double> binIdx(x.size() * x.size());
	for (int i = 0; i < edges.size(); ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (dMat[j][2].real() < edges[i + 1] && dMat[j][2].real() >= edges[i]) {
				binCounts[i] += 1.0;
				binIdx[j] = i;
			}
		}
	}

	//calculate average
	std::vector<long double> average(edges.size() - 1);
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

	//calculate covariogram
	std::vector<long double> covariogram(edges.size() - 1);
	for (int i = 0; i < edges.size() - 1.0; ++i) {
		for (int j = 0; j < binIdx.size(); ++j) {
			if (binCounts[i] != 0) {
				if (binIdx[j] == i) {
					//variogram[int(binIdx[j])] += std::sqrt(dMat[j][3]).real();
					//covariogram[int(binIdx[j])] += abs((dMat[j][4] - average[i]) * (dMat[j][5] - average[i]));
					covariogram[int(binIdx[j])] += abs((dMat[j][4] - avg) * (dMat[j][5] - avg));
				}
			}
		}
	}
	for (int i = 0; i < covariogram.size(); ++i) {

		//variogram[i] = std::pow(variogram[i] / binCounts[i], 4) / 2.0 / (0.457 + 0.494 / binCounts[i] + 0.045 / binCounts[i] / binCounts[i]);
		covariogram[i] = covariogram[i] / (binCounts[i]);
	}
	return covariogram;

}

static long double calcVariance2(std::vector<long double> variogramFit, std::vector<long double> variogramGradFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	long double variance = 0.0;
	const int size = xSample.size();
	MatrixXd varMat(size, size);
	MatrixXd varMatGrad(size, size);

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
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
		long double xDist = abs(xSample[i].real() - xVal);
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

	MatrixXd varianceVec = weightVec.transpose() * varMat * weightVec + weightVecGrad.transpose() * varMatGrad * weightVecGrad - covVec.transpose() * weightVec - weightVec.transpose() * covVec;
	//MatrixXd varianceVec = weightVec.transpose() * varMat * weightVec - covVec.transpose() * weightVec - weightVec.transpose() * covVec;

	//std::cout << "old variance: " << varianceVec;
	variance = abs(varianceVec(0, 0));

	//scale relative
	//variance = variance / xVal;
	//std::cout << "variance:\n" << varianceVec << std::endl;


	//variance = 0.0;
	////from textbook source:
	//for (int i = 0; i < size; ++i) {
	//	for (int j = 0; j < size; ++j) {

	//		long double dist = abs(xSample[i].real() - xSample[j].real());
	//		auto it = std::lower_bound(d.begin(), d.end(), dist);
	//		int index = it - d.begin();
	//		
	//		variance -= weights[i] * weights[j] * variogramFit[index];
	//	}
	//}

	//for (int i = 0; i < size; ++i) {

	//	long double dist = abs(xSample[i].real() - xVal);
	//	auto it = std::lower_bound(d.begin(), d.end(), dist);
	//	int index = it - d.begin();

	//	variance += 2.0 * weights[i] * variogramFit[index];
	//}

	////std::cout << "\tnew variance: " << variance << std::endl;


	return variance;
}

static long double calcVariance(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	long double variance = 0.0;
	const int size = xSample.size();
	MatrixXd varMat(size, size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			varMat(i, j) = variogramFit[index];
		}
	}
	MatrixXd weightVec(size, 1);
	MatrixXd covVec(size, 1);
	for (int i = 0; i < size; ++i) {
		long double xDist = abs(xSample[i].real() - xVal);
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

	//scale relative
	//variance = variance / xVal;
	//std::cout << "variance:\n" << varianceVec << std::endl;


	//variance = 0.0;
	////from textbook source:
	//for (int i = 0; i < size; ++i) {
	//	for (int j = 0; j < size; ++j) {

	//		long double dist = abs(xSample[i].real() - xSample[j].real());
	//		auto it = std::lower_bound(d.begin(), d.end(), dist);
	//		int index = it - d.begin();
	//		
	//		variance -= weights[i] * weights[j] * variogramFit[index];
	//	}
	//}

	//for (int i = 0; i < size; ++i) {

	//	long double dist = abs(xSample[i].real() - xVal);
	//	auto it = std::lower_bound(d.begin(), d.end(), dist);
	//	int index = it - d.begin();

	//	variance += 2.0 * weights[i] * variogramFit[index];
	//}

	////std::cout << "\tnew variance: " << variance << std::endl;


	return variance;
}

static void kriging(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	const int size = xSample.size();
	//Matrix <long double, Dynamic, Dynamic> krigMat;
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
				long double dist = abs(xSample[i].real() - xSample[j].real());
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
		long double xDist = abs(xVal - xSample[i].real());
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

long double fact(int n) {
	long double factorial;
	for (int i = 1; i <= n; ++i) {
		factorial *= double(i);
	}
	return factorial;
}

static void bezierKriging(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();

	//if (M > 18) { M = 18; }

	//Matrix <long double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd fMat(size, M);
	MatrixXd zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	long double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	long double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
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
		long double xDist = abs(xVal - xSample[i].real());
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

static void fourierKriging(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();
	//int M = 25;

	if (M > 4) { M = 4; }

	//Matrix <long double, Dynamic, Dynamic> krigMat;
	Matrix <std::complex<long double>, Dynamic, Dynamic> sMat(size, size);
	Matrix <std::complex<long double>, Dynamic, Dynamic> fMat(size, M);
	Matrix <std::complex<long double>, Dynamic, Dynamic> zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	long double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	long double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
		}
	}


	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			//fMat(i, j) = std::pow((xSample[i].real() - average), j);

			//if (j == 0)
			//	fMat(i, j) = 1.0;
			//else if (j % 2 == 0) {
			//	fMat(i, j) = std::cos(2 * 3.14159 * j * xSample[i].real() / 12.0);
			//	//std::cout << "cos function: " << std::cos(2 * 3.14159 * j * xSample[i].real()/12.0) << "  j: " << j << "  xSample[i].real(): " << xSample[i].real() << std::endl;
			//}
			//else {
			//	fMat(i, j) = std::sin(2 * 3.14159 * j * xSample[i].real() / 12.0);
			//	//std::cout << "sin function: " << std::sin(2 * 3.14159 * j * xSample[i].real()/12.0) << "  j: " << j << "  xSample[i].real(): " << xSample[i].real() << std::endl;
			//}


			std::complex<long double> num(cos(2 * 3.14159 * j * (xSample[i].real() - average) / average), sin(2 * 3.14159 * j * (xSample[i].real() - average) / average));
			fMat(i, j) = num;

			//fMat(i, j) = std::exp(-std::pow((xSample[i].real() - average), 2)*j);



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

	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			zeroMat(i, j) = 0.0;
		}
	}

	Matrix <std::complex<long double>, Dynamic, Dynamic> upper(size, size + M);
	Matrix <std::complex<long double>, Dynamic, Dynamic> lower(M, size + M);
	Matrix <std::complex<long double>, Dynamic, Dynamic> krigMat(size + M, size + M);

	upper << sMat, fMat;
	lower << fMat.transpose(), zeroMat;

	krigMat << upper, lower;

	//std::cout << "krigMat: " << krigMat << std::endl;

	Matrix <std::complex<long double>, Dynamic, Dynamic> cVec(size, 1);
	Matrix <std::complex<long double>, Dynamic, Dynamic> fVec(M, 1);
	Matrix <std::complex<long double>, Dynamic, Dynamic> krigVec(size + M, 1);

	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		long double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
	}

	for (int i = 0; i < M; ++i) {

		std::complex<long double> num(cos(2 * 3.14159 * i * (xVal - average) / average), sin(2 * 3.14159 * i * (xVal - average) / average));
		fVec(i, 0) = num;
		//fVec(i, 0) = std::exp(-std::pow(xVal - average, 2)*i);
		//if (i == 0)
		//	fVec(i, 0) = 1.0;
		//else if (i % 2 == 0)
		//	fVec(i, 0) = std::cos(2 * 3.14159 * i * xVal / 12.0);
		//else
		//	fVec(i, 0) = std::sin(2 * 3.14159 * i * xVal / 12.0);


	}

	krigVec << cVec, fVec;

	//std::cout << "krig matrix:\n";
	//std::cout << krigMat << std::endl;
	//std::cout << "krig vec:\n";
	//std::cout << krigVec << std::endl;

	Matrix <std::complex<long double>, Dynamic, Dynamic> weightMat = krigMat.inverse() * krigVec;
	//std::cout << "weightMat: " << weightMat << std::endl;

	for (int i = 0; i < size; ++i) {
		//std::cout << weightMat(i, 0) << " ";
		weights.push_back(weightMat(i, 0).real());
	}
	//std::cout << "weight vector: " << weightMat << std::endl;
	//std::cout << std::endl;

}

static void taylorKriging(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();

	//if (M > 3) { M = 18; }

	//Matrix <long double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd fMat(size, M);
	MatrixXd zeroMat(M, M);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	long double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	long double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
		}
	}


	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			fMat(i, j) = std::pow((xSample[i].real() - average), j);

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
		long double xDist = abs(xVal - xSample[i].real());
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

static void taylorKrigingHigherOrder(std::vector<long double> variogramFit, std::vector<long double> variogramFitFirstD, std::vector<long double> variogramFitSecondD, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal, std::vector<long double> variogramFitGrad) {
	
	//size is N in calculations
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();

	//if (M > 3) { M = 18; }

	//Matrix <long double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd sMatFirstDerivative(size, size);
	MatrixXd sMatSecondDerivative(size, size);
	MatrixXd fMatDerivative(size, M);
	MatrixXd fMat(size, M);
	MatrixXd zeroMatM(M, M);
	MatrixXd zeroMatN(size, size);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	long double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	long double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
			sMatFirstDerivative(i, j) = variogramFitFirstD[index];
			sMatSecondDerivative(i, j) = variogramFitSecondD[index];
		}
	}

	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			fMat(i, j) = std::pow((xSample[i].real() - average), j);


			if (j == 0)
				fMatDerivative(i, j) = 0.0;
			else
				fMatDerivative(i, j) = j * std::pow((xSample[i].real() - average), j - 1);

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

	upper << sMat, sMatFirstDerivative, fMat;

	mid << sMatFirstDerivative, sMatSecondDerivative, fMatDerivative;

	lower << fMat.transpose(), fMatDerivative.transpose(), zeroMatM;

	krigMat << upper, mid, lower;

	//std::cout << "krigMat: " << krigMat << std::endl;

	MatrixXd cVec(size, 1);
	MatrixXd cVecDerivative(size, 1);
	MatrixXd fVec(M, 1);
	MatrixXd krigVec(size + size + M, 1);

	//C vector
	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		long double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
		cVecDerivative(i, 0) = variogramFitFirstD[index];
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal - average), i);
	}

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

	//std::cout << weightMat << std::endl;
	//std::cout << std::endl;

}

static long double taylorKrigingHigherOrder2(std::vector<long double> variogramFit, std::vector<long double> variogramFitValGrad, std::vector<long double> variogramFitGradVal, std::vector<long double> variogramFitGrad, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal) {

	//size is N in calculations
	const int size = xSample.size();

	//order of the basis functions is M - 1; differs from literature
	int M = xSample.size();
	//trying for "ordinary kriging"
	//int M = 1;

	//if (M > 3) { M = 18; }

	//Matrix <long double, Dynamic, Dynamic> krigMat;
	MatrixXd sMat(size, size);
	MatrixXd sMatFirstDerivative(size, size);
	MatrixXd sMatSecondDerivative(size, size);
	MatrixXd fMatDerivative(size, M);
	MatrixXd fMat(size, M);
	MatrixXd zeroMatM(M, M);
	MatrixXd zeroMatN(size, size);

	//std::cout << "size: " << size << " rows: " << krigMat.rows() << " cols: " << krigMat.cols() << std::endl;

	//calculate sample average
	long double sum = 0.0;
	for (int i = 0; i < xSample.size(); ++i) {
		sum += xSample[i].real();
	}
	long double average = sum / xSample.size();

	// S matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			long double dist = abs(xSample[i].real() - xSample[j].real());
			auto it = std::lower_bound(d.begin(), d.end(), dist);
			int index = it - d.begin();
			sMat(i, j) = variogramFit[index];
			//sMatFirstDerivative(i, j) = variogramFitValGrad[index];
			sMatFirstDerivative(i, j) = 0.0;
			sMatSecondDerivative(i, j) = variogramFitGrad[index];
		}
	}

	// F Matrix
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < M; ++j) {


			fMat(i, j) = std::pow((xSample[i].real() - average), j);
			//fMat(i, j) = 1.0;

			fMatDerivative(i, j) = j*std::pow((xSample[i].real() - average), j-1);
			//if (j == 0)
			//	fMatDerivative(i, j) = 0.0;
			//fMatDerivative(i, j) = 0.0;

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

	upper << sMat, sMatFirstDerivative, fMat;

	mid << sMatFirstDerivative, sMatSecondDerivative, fMatDerivative;

	lower << fMat.transpose(), fMatDerivative.transpose(), zeroMatM;

	krigMat << upper, mid, lower;

	//std::cout << "krigMat: " << krigMat << std::endl;

	MatrixXd cVec(size, 1);
	MatrixXd cVecDerivative(size, 1);
	MatrixXd fVec(M, 1);
	MatrixXd krigVec(size + size + M, 1);

	//C vector
	for (int i = 0; i < size; ++i) {
		//std::cout << "i: " << i;
		long double xDist = abs(xVal - xSample[i].real());
		auto it = std::lower_bound(d.begin(), d.end(), xDist);
		int index = it - d.begin();
		cVec(i, 0) = variogramFit[index];
		//cVecDerivative(i, 0) = variogramFitGradVal[index];
		cVecDerivative(i, 0) = 0.0;
	}

	for (int i = 0; i < M; ++i) {
		fVec(i, 0) = std::pow((xVal - average), i);
		//fVec(i, 0) = 1.0;
	}

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
	//std::cout << "matrix: \n";
	//std::cout << krigMat << std::endl;
	//std::cout << "vector: \n";
	//std::cout << krigVec << std::endl;
	//std::cout << "weights: \n";
	//std::cout << weightMat << std::endl;
	////std::cout << std::endl;
	return average;
}


static std::complex<long double> krigSum(std::vector<long double> weights, std::vector<std::complex<long double>> referenceVals) {
	std::complex<long double> sum = 0.0;
	for (int i = 0; i < weights.size(); ++i) {
		sum += weights[i] * referenceVals[i];
	}
	return sum;
}

static std::complex<long double> krigSumHigherOrder(std::vector<long double> weights, std::vector<std::complex<long double>> referenceVals, std::vector<std::complex<long double>> gradient, long double xVal, long double xAvg) {
	std::complex<long double> sum = 0.0;
	std::complex<long double> sum2 = 0.0;


	//std::cout << "started krigSumHigher, referenceVals.size()  "  << referenceVals.size() << " gradient.size(): " << gradient.size() << "  weights.size(): " << weights.size() << std::endl;

	for (int i = 0; i < referenceVals.size(); ++i) {

		//if (i < referenceVals.size()) {
		//	sum += referenceVals[i] * weights[i];
		//}
		//else
		//	sum += gradient[i-referenceVals.size()] * weights[i];
		sum += referenceVals[i] * weights[i] + gradient[i] * weights[i + referenceVals.size()];// *(xVal - xAvg);

		//std::cout << "referenceVals[i]: " << referenceVals[i] << "  weights[i]: " << weights[i] << "  gradient[i]: " << gradient[i] << "  weights[i + size]: " << weights[i + referenceVals.size()] << std::endl;

	}


	return sum;
}

static void trimVariogram(std::vector<long double>& variogram, std::vector<long double>& edges) {
	int idx = 0;
	std::vector<long double> variogramNew;
	std::vector<long double> edgesNew;

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

static void trimVariogramComplex(std::vector<std::complex<long double>>& variogram, std::vector<long double>& edges) {
	int idx = 0;
	std::vector<std::complex<long double>> variogramNew;
	std::vector<long double> edgesNew;

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

void function_fit_gaussian(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(-x[0] * x[0] / (c[1] * c[1])));
}

void function_fit_gaussianCo(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (std::exp(-x[0] * x[0] / (c[1] * c[1])));
}

void function_fit_exp(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(-x[0] / (c[1])));
}

void function_fit_expCo(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (std::exp(-x[0] / (c[1])));
}

void function_fit_stable(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[2] + c[0] * (1.0 - std::exp(std::pow(-x[0] / c[1], int(c[3]))));
}

void function_fit_linear(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	func = c[0] * x[0] + c[1];
}

void function_fit_spherical(const alglib::real_1d_array& c, const alglib::real_1d_array& x, double& func, void* ptr) {
	if (x[0] < c[1])
		func = c[2] + c[0] * (1.5 * x[0] / c[1] - 0.5 * x[0] * x[0] * x[0] / c[1] / c[1] / c[1]);
	else
		func = c[2] + c[0];
}

std::vector<long double> fitVariogram(std::vector<long double> variogram, std::vector<long double> edges, std::vector<long double> d) {
	long double xSum = 0.0;
	long double ySum = 0.0;
	long double x2Sum = 0.0;
	long double xySum = 0.0;
	for (int i = 0; i < variogram.size(); ++i) {
		xSum += edges[i];
		ySum += variogram[i];
		x2Sum += x2Sum + (edges[i] * edges[i]);
		xySum += xySum + (edges[i] * variogram[i]);
	}
	long double slope = (edges.size() * xySum - xSum * ySum) / (edges.size() * x2Sum - xSum * xSum);
	long double intercept = (x2Sum * ySum - xSum * xySum) / (x2Sum * edges.size() - xSum * xSum);
	std::vector<long double> fitVariogram(d.size());
	
	//std::cout << "m: " << slope << std::endl;

	for (int i = 0; i < d.size(); ++i) {
		fitVariogram[i] = slope * d[i] + intercept;
	}
	return fitVariogram;
}


std::vector<long double> fitVariogramModel(std::vector<long double> variogram, std::vector<long double> edges, std::vector<long double> d, long double size) {
	alglib::real_2d_array x;
	alglib::real_1d_array y;
	y.setlength(variogram.size());
	x.setlength(variogram.size(), 1);

	for (int i = 0; i < variogram.size(); ++i) {
		y[i] = variogram[i];
		x[i][0] = edges[i];
	}

	long double epsx = 0.000001;
	//long double epsx = 0;
	alglib::ae_int_t maxits = 0;
	alglib::ae_int_t info;
	alglib::lsfitstate state;
	alglib::lsfitreport rep;
	long double diffstep = 0.0001;
	std::vector<long double> fitVariogram(d.size());


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

		///////////////////////////gaussian fit//////////////////////////////////////////
	
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
		//alglib::real_1d_array c;
		// c.setlength(3);
		//c[0] = *std::max_element(variogram.begin(), variogram.end());
		//c[1] = edges[edges.size() - 1];
		//c[2] = 0.0;
 	//	alglib::real_1d_array bndl;
		//alglib::real_1d_array bndu;
		//bndl.setlength(3);
		//bndu.setlength(3);

		//bndl[0] = *std::min_element(variogram.begin(), variogram.end());
		//bndl[1] = edges[0];
		//bndl[2] = 0.0;
		//bndu[0] = c[0];
		//bndu[1] = c[1];
		//bndu[2] = bndu[1]*0.99;

		//alglib::lsfitcreatef(x, y, c, diffstep, state);
		//alglib::lsfitsetcond(state, epsx, maxits);
		//alglib::lsfitsetbc(state, bndl, bndu);
		//alglib::lsfitfit(state, function_fit_spherical);
		//alglib::lsfitresults(state, info, c, rep);
		//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

		//for (int i = 0; i < d.size(); ++i) {
		//	if (d[i] < c[1])
		//		fitVariogram[i] = c[2] + c[0] * (1.5 * d[i] / c[1] - 0.5 * d[i] * d[i] * d[i] / c[1] / c[1] / c[1]);
		//	else
		//		fitVariogram[i] = c[2] + c[0];

		//}
	long double max = *std::max_element(fitVariogram.begin(), fitVariogram.end());
	for (int i = 0; i < fitVariogram.size(); ++i) {
		//fitVariogram[i] = fitVariogram[i] / max;
	}

	return fitVariogram;

}

std::vector<long double> fitCoVariogramModel(std::vector<long double> variogram, std::vector<long double> edges, std::vector<long double> d, long double size) {
	alglib::real_2d_array x;
	alglib::real_1d_array y;
	y.setlength(variogram.size());
	x.setlength(variogram.size(), 1);

	for (int i = 0; i < variogram.size(); ++i) {
		y[i] = variogram[i];
		x[i][0] = edges[i];
	}

	std::cout << "start of fit method\n";

	long double epsx = 0.000001;
	//long double epsx = 0;
	alglib::ae_int_t maxits = 0;
	alglib::ae_int_t info;
	alglib::lsfitstate state;
	alglib::lsfitreport rep;
	long double diffstep = 0.0001;
	std::vector<long double> fitVariogram(d.size());


	//////////////////////////linear fit////////////////////////////////////////////
	
	alglib::real_1d_array c = "[0.25, 0.25]";
	alglib::lsfitcreatef(x, y, c, diffstep, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitfit(state, function_fit_linear);
	alglib::lsfitresults(state, info, c, rep);
	std::cout << "c: " << c[0] << "  a: " << c[1] << std::endl;


	for (int i = 0; i < d.size(); ++i) {
		fitVariogram[i] = c[0] * d[i] + c[1];
	}
	
	/////////////////////////gaussian fit//////////////////////////////////////////
	/*
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
	//alglib::lsfitsetbc(state, bndl, bndu);

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
	*/

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
	//alglib::real_1d_array c;
	// c.setlength(3);
	//c[0] = *std::max_element(variogram.begin(), variogram.end());
	//c[1] = edges[edges.size() - 1];
	//c[2] = 0.0;
//	alglib::real_1d_array bndl;
	//alglib::real_1d_array bndu;
	//bndl.setlength(3);
	//bndu.setlength(3);

	//bndl[0] = *std::min_element(variogram.begin(), variogram.end());
	//bndl[1] = edges[0];
	//bndl[2] = 0.0;
	//bndu[0] = c[0];
	//bndu[1] = c[1];
	//bndu[2] = bndu[1]*0.99;

	//alglib::lsfitcreatef(x, y, c, diffstep, state);
	//alglib::lsfitsetcond(state, epsx, maxits);
	//alglib::lsfitsetbc(state, bndl, bndu);
	//alglib::lsfitfit(state, function_fit_spherical);
	//alglib::lsfitresults(state, info, c, rep);
	//std::cout << "c: " << c[0] << "  a: " << c[1] << " b: " << c[2] << std::endl;

	//for (int i = 0; i < d.size(); ++i) {
	//	if (d[i] < c[1])
	//		fitVariogram[i] = c[2] + c[0] * (1.5 * d[i] / c[1] - 0.5 * d[i] * d[i] * d[i] / c[1] / c[1] / c[1]);
	//	else
	//		fitVariogram[i] = c[2] + c[0];

	//}

	long double max = *std::max_element(fitVariogram.begin(), fitVariogram.end());
	for (int i = 0; i < fitVariogram.size(); ++i) {
		fitVariogram[i] = fitVariogram[i] / max;
	}

	return fitVariogram;

}


std::complex<long double> Kriging::sensitivity_to_epsr(std::string & file_name, std::vector<std::complex<long double>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::complex<long double> referenceFreq, std::complex<long double>& gradient)
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
	long double k0 = dom_forward.scatter1.K0[1];
	std::complex<long double> reference_qoi = get_QoI(dom_forward, dom_adjoint);
	//std::complex<long double> reference_qoi_compare = dom_forward.scatter1.
	std::cout << "reference qoi: " << reference_qoi << " reference mat: " << referenceFreq.real() << std::endl;
	
	std::complex<long double> eps_diff = 0.0;
	std::complex<long double> reference_material = 0.0;
	long double reference_k0;
	std::complex<long double> new_qoi;
	for (auto e = dom_forward.elements.begin(); e != dom_forward.elements.end(); ++e) {
		if (e->materials.region == 1) {
			reference_material = e->materials.epsr_list[0][0]; //only works for homogeneous scatterers at the moment!
			break;
		}
	}

	//get the remainder of the estimate
	std::complex<long double> dQoI = 0.0;

	std::complex<long double> dQoI2 = 0.0, dQoI1 = 0.0, dQoI2_ver1 = 0.0, dQoI2_ver2 = 0.0;
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
					//dQoI2_ver2 += eps_diff * k0*k0*e->cBrPAK(i0, j0)*dom_forward.cAlpha[icon] * std::conj(dom_adjoint.cAlpha[jcon]); //this version is much worse, about long double the error
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


	long double kReal = 2 * 3.14159 / (2.99E8 / referenceFreq.real());
	std::complex<long double> referenceK = (kReal, 0.0);

	//now have the base QoI modifier, need to scale by the change in perturbed material parameter
	//for (int i = 0; i < epsr_list.size(); ++i) {
	//	/////////////////////////////////////////Frequency Perturbation/////////////////////////////////////////////
	//	//long double k0_list = epsr_list[i].real() * 2.0 * 3.14159 / 3e8;
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

std::complex<long double> Kriging::get_QoI(Domain & dom_forward, Domain & dom_adjoint)
{
	std::complex<long double> qoi = 0.0;
	for (int i = 0; i < dom_forward.scatter1.cGr[1].size(); ++i) {
		qoi += dom_forward.scatter1.cGr[1][i] * std::conj(dom_adjoint.cAlpha[i]);
	}
	return qoi;
}
void Kriging::monte_carlo_instance(std::string & file_name)
{
	//compute the forward solution
	//then compute just the RHS for the adjoint solution (to get the QoI term to compare w/ HOPS)
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
	std::complex<long double> qoi = Kriging::get_QoI2(dom_forward, dom_adjoint_rhs);
}

void Kriging::multi_HOPS_epsr(std::string& file_name)
{


	//load in the list of random materials being tested
	std::vector<std::complex<long double>> material_list;

	//old material range
	//std::ifstream materials_in("../ioFiles/input/materials_list.txt");

	//new material range with monte carlo gradient
	std::ifstream materials_in("../ioFiles/input/materials_list3.txt");
	std::string line;
	//std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		//material_list.push_back(std::complex<long double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
		material_list.push_back(std::complex<long double>(std::stod(data[0]), std::stod(data[1])));
	}
	materials_in.close();

	//for sweeping ref inputs, otherwise use next segment for mat or freq perturbation with set references---------------------------------
	//std::vector<long double> referencesReal;
	//std::string refString = "../reference_files_mat/input_ref_list.txt";
	//std::ifstream refs_in(refString);
	//std::string line2;
	////std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	//while (std::getline(refs_in, line2)) {
	//	auto data = functions::split(line2, ' ');
	//	//material_list.push_back(std::complex<long double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
	//	referencesReal.push_back(std::stod(data[0]));
	//}
	//refs_in.close();
	//std::vector<std::complex<long double>> referencesFull;
	//for (int i = 0; i < referencesReal.size(); i++) {
	//	referencesFull.push_back(std::complex<long double>(referencesReal[i], -2.0));
	//}

	sort(material_list.begin(), material_list.end(), complexComparitor);

	//acceptable error
	long double delta = 0.0008;

	//convergence condition
	long double stdLim = 0.0001;

	//iteration limit
	int itLim = 26;

	//material_list.txt
	//long double matStd = 1.0;
	//long double matMean = 4.5;

	//material_list2.txt
	long double matStd = 1.0;
	long double matMean = 4.5;

	std::vector<long double> stdDevRCSList;
	std::vector<std::complex<long double>> stdDevList;


	///////////////////////Material Perturbation/////////////////////////////////////////////////////////////////////
	//std::vector<std::complex<long double>> references = { {4.0, -2.0} };
	//std::vector<std::complex<long double>> referencesFull = { {1.000000, -2.000000},{1.666667, -2.000000},{2.333333, -2.000000},{3.000000, -2.000000},{3.666667, -2.000000},{4.333333, -2.000000},{5.000000, -2.000000},{5.666667, -2.000000},{6.333333, -2.000000},{7.000000, -2.000000} };
	//std::vector<std::complex<long double>> referencesFull = { {1.000000, -2.000000},{1.300000, -2.000000},{1.600000, -2.000000},{1.900000, -2.000000},{2.200000, -2.000000},{2.500000, -2.000000},{2.800000, -2.000000},{3.100000, -2.000000},{3.400000, -2.000000},{3.700000, -2.000000},{4.000000, -2.000000},{4.300000, -2.000000},{4.600000, -2.000000},{4.900000, -2.000000},{5.200000, -2.000000},{5.500000, -2.000000},{5.800000, -2.000000},{6.100000, -2.000000},{6.400000, -2.000000},{6.700000, -2.000000},{7.000000, -2.000000} };
	//std::vector<std::complex<long double>> referencesFull = { {2.000000, -2.000000},{2.250000, -2.000000},{2.500000, -2.000000},{2.750000, -2.000000},{3.000000, -2.000000},{3.250000, -2.000000},{3.500000, -2.000000},{3.750000, -2.000000},{4.000000, -2.000000},{4.250000, -2.000000},{4.500000, -2.000000},{4.750000, -2.000000},{5.000000, -2.000000},{5.250000, -2.000000},{5.500000, -2.000000},{5.750000, -2.000000},{6.000000, -2.000000},{6.250000, -2.000000},{6.500000, -2.000000},{6.750000, -2.000000},{7.000000, -2.000000} };
	
	
	//old frequency and material range
	//std::vector<std::complex<long double>> referencesFull = { {1.500000, -2.000000}, { 1.600000, -2.000000 }, { 1.700000, -2.000000 }, { 1.800000, -2.000000 }, { 1.900000, -2.000000 }, { 2.000000, -2.000000 }, { 2.100000, -2.000000 }, { 2.200000, -2.000000 }, { 2.300000, -2.000000 }, { 2.400000, -2.000000 }, { 2.500000, -2.000000 }, { 2.600000, -2.000000 }, { 2.700000, -2.000000 }, { 2.800000, -2.000000 }, { 2.900000, -2.000000 }, { 3.000000, -2.000000 }, { 3.100000, -2.000000 }, { 3.200000, -2.000000 }, { 3.300000, -2.000000 }, { 3.400000, -2.000000 }, { 3.500000, -2.000000 }, { 3.600000, -2.000000 }, { 3.700000, -2.000000 }, { 3.800000, -2.000000 }, { 3.900000, -2.000000 }, { 4.000000, -2.000000 }, { 4.100000, -2.000000 }, { 4.200000, -2.000000 }, { 4.300000, -2.000000 }, { 4.400000, -2.000000 }, { 4.500000, -2.000000 }, { 4.600000, -2.000000 }, { 4.700000, -2.000000 }, { 4.800000, -2.000000 }, { 4.900000, -2.000000 }, { 5.000000, -2.000000 }, { 5.100000, -2.000000 }, { 5.200000, -2.000000 }, { 5.300000, -2.000000 }, { 5.400000, -2.000000 }, { 5.500000, -2.000000 }, { 5.600000, -2.000000 }, { 5.700000, -2.000000 }, { 5.800000, -2.000000 }, { 5.900000, -2.000000 }, { 6.000000, -2.000000 }, { 6.100000, -2.000000 }, { 6.200000, -2.000000 }, { 6.300000, -2.000000 }, { 6.400000, -2.000000 }, { 6.500000, -2.000000 }, { 6.600000, -2.000000 }, { 6.700000, -2.000000 }, { 6.800000, -2.000000 }, { 6.900000, -2.000000 }, { 7.000000, -2.000000 }, { 7.100000, -2.000000 }, { 7.200000, -2.000000 }, { 7.300000, -2.000000 }, { 7.400000, -2.000000 }, { 7.500000, -2.000000 }, { 7.600000, -2.000000 }, { 7.700000, -2.000000 }, { 7.800000, -2.000000 }, { 7.900000, -2.000000 }, { 8.000000, -2.000000 } };

	//new material range
	std::vector<std::complex<long double>> referencesFull = { {3.000000, -2.000000},{3.100000, -2.000000},{3.200000, -2.000000},{3.300000, -2.000000},{3.400000, -2.000000},{3.500000, -2.000000},{3.600000, -2.000000},{3.700000, -2.000000},{3.800000, -2.000000},{3.900000, -2.000000},{4.000000, -2.000000},{4.100000, -2.000000},{4.200000, -2.000000},{4.300000, -2.000000},{4.400000, -2.000000},{4.500000, -2.000000},{4.600000, -2.000000},{4.700000, -2.000000},{4.800000, -2.000000},{4.900000, -2.000000},{5.000000, -2.000000},{5.100000, -2.000000},{5.200000, -2.000000},{5.300000, -2.000000},{5.400000, -2.000000},{5.500000, -2.000000},{5.600000, -2.000000},{5.700000, -2.000000},{5.800000, -2.000000},{5.900000, -2.000000},{6.000000, -2.000000},{6.100000, -2.000000},{6.200000, -2.000000},{6.300000, -2.000000},{6.400000, -2.000000},{6.500000, -2.000000},{6.600000, -2.000000},{6.700000, -2.000000},{6.800000, -2.000000},{6.900000, -2.000000},{7.000000, -2.000000},{7.100000, -2.000000},{7.200000, -2.000000},{7.300000, -2.000000},{7.400000, -2.000000},{7.500000, -2.000000},{7.600000, -2.000000},{7.700000, -2.000000},{7.800000, -2.000000},{7.900000, -2.000000},{8.000000, -2.000000},{8.100000, -2.000000},{8.200000, -2.000000},{8.300000, -2.000000},{8.400000, -2.000000},{8.500000, -2.000000},{8.600000, -2.000000},{8.700000, -2.000000},{8.800000, -2.000000},{8.900000, -2.000000},{9.000000, -2.000000},{9.100000, -2.000000},{9.200000, -2.000000},{9.300000, -2.000000},{9.400000, -2.000000},{9.500000, -2.000000},{9.600000, -2.000000},{9.700000, -2.000000},{9.800000, -2.000000},{9.900000, -2.000000},{10.000000, -2.000000} };
	//----------------------------------------------------------------------------------------------------------------

	////////////////Frequency Perturbation//////////////////////////////////////////////////////////////////
	//std::vector<std::complex<long double>> references = { 10000000.000000, 12142857.142857, 14285714.285714, 16428571.428571, 18571428.571429, 20714285.714286, 22857142.857143, 25000000.000000, 27142857.142857, 29285714.285714, 31428571.428571, 33571428.571429, 35714285.714286, 37857142.857143, 40000000.000000};
	//std::vector<std::complex<long double>> references = { 10000000.000000, 17500000.000000, 25000000.000000, 32500000.000000, 40000000.000000 };

	//std::vector<std::complex<long double>> references = { 10000000.000000, 12727272.727273, 15454545.454545, 18181818.181818, 20909090.909091, 23636363.636364, 26363636.363636, 29090909.090909, 31818181.818182, 34545454.545455, 37272727.272727, 40000000.000000 };
	//std::vector<std::complex<long double>> references = { 10000000.000000, 11500000.000000, 13000000.000000, 14500000.000000, 16000000.000000, 17500000.000000, 19000000.000000, 20500000.000000, 22000000.000000, 23500000.000000, 25000000.000000, 26500000.000000, 28000000.000000, 29500000.000000, 31000000.000000, 32500000.000000, 34000000.000000, 35500000.000000, 37000000.000000, 38500000.000000, 40000000.000000};
	//std::vector<std::complex<long double>> references = { 15000000,15416666.6666667,15833333.3333333,16250000,16666666.6666667,17083333.3333333,17500000,17916666.6666667,18333333.3333333,18750000,19166666.6666667,19583333.3333333,20000000,20416666.6666667,20833333.3333333,21250000,21666666.6666667,22083333.3333333,22500000,22916666.6666667,23333333.3333333,23750000,24166666.6666667,24583333.3333333,25000000,25416666.6666667,25833333.3333333,26250000,26666666.6666667,27083333.3333333,27500000,27916666.6666667,28333333.3333333,28750000,29166666.6666667,29583333.3333333,30000000,30416666.6666667,30833333.3333333,31250000,31666666.6666667,32083333.3333333,32500000,32916666.6666667,33333333.3333333,33750000,34166666.6666667,34583333.3333333,35000000,35416666.6666667,35833333.3333333,36250000,36666666.6666667,37083333.3333333,37500000,37916666.6666667,38333333.3333333,38750000,39166666.6666667,39583333.3333333,40000000 };
	////-----------------------------------------------------------------------------------------------------

	long double size = referencesFull.size();

	int loRefIndex = find_index(referencesFull, material_list[0]);
	if (loRefIndex < 0) { loRefIndex = 0; }

	int hiRefIndex = find_lowest(referencesFull, material_list[material_list.size() - 1]);
	if (hiRefIndex > size - 1) { hiRefIndex = size - 1; }

	int midRefIndex = (loRefIndex + hiRefIndex) / 2;

	//initial set of the references to be passed into HOPS (first, middle, last) requires odd number of references in the full vector
	//std::vector<std::complex<long double>> references = { referencesFull[0], referencesFull[(size - 1) / 2], referencesFull[size - 1] };
	std::vector<std::complex<long double>> references = { referencesFull[loRefIndex], referencesFull[midRefIndex], referencesFull[hiRefIndex] };
	//references = { referencesFull[5], referencesFull[10], referencesFull[15], referencesFull[20], referencesFull[25], referencesFull[30],referencesFull[35], referencesFull[40], referencesFull[45], referencesFull[50], referencesFull[55], referencesFull[60], referencesFull[65], referencesFull[70] };
	//references = { referencesFull[5], referencesFull[9], referencesFull[13], referencesFull[17], referencesFull[21], referencesFull[25], referencesFull[29], referencesFull[33], referencesFull[37], referencesFull[41], referencesFull[45], referencesFull[49], referencesFull[53], referencesFull[57], referencesFull[61], referencesFull[65], referencesFull[70] };
	//references = referencesFull;


	std::vector<std::vector<std::complex<long double>>> qoi_list(references.size());
	std::vector<std::vector<std::complex<long double>>> HOPS_splitting(references.size());

	//refIndex is the index of the ref_i file so that it can be called correctly from the reference_files folder
	std::vector<int> refIndex = { loRefIndex, midRefIndex, hiRefIndex };
	//std::vector<int> refIndex = { 0, int(size - 1) / 2, int(size - 1) };
	//refIndex = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70 };
	//refIndex = { 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 70 };



	//for (int i = 0; i < 71; i++) {
	//	refIndex.push_back(i);
	//}

	////easy sweep condition
	//refIndex.clear();
	//references.clear();
	//for (int i = 0; i < 71; i += 2) {
	//	refIndex.push_back(i);
	//	references.push_back(referencesFull[i]);
	//}
	//if (refIndex[refIndex.size() - 1] != 70) {
	//	refIndex.push_back(70);
	//	references.push_back(70);
	//}
	////-------


	//toggle to end the while loop
	bool toggle = true;
	bool whileCond = true;
	//some declarations
	std::string in_file_name;

	int iterations = 0;

	std::vector<std::vector<long double>> RCS(references.size());
	std::complex<long double> sum = 0.0;
	long double sumRCS = 0.0;

	std::complex<long double> mean;
	long double meanRCS;

	std::complex<long double> stdDev = 0.0;
	long double stdDevRCS = 0.0;
	//std::vector<long double> pdfScaleFactor(references.size());

	std::vector<std::vector<long double>> weights(material_list.size());
	std::vector<std::vector<long double>> weightsGrad(material_list.size());
	std::vector<std::vector<long double>> weightsHigher(material_list.size());

	std::vector<std::complex<long double>> reconstruction(material_list.size());
	std::vector<std::complex<long double>> reconstructionHigher(material_list.size());
	std::vector<std::complex<long double>> reconstructionUpdate(material_list.size());
	std::vector<std::complex<long double>> reconstructionGrad(material_list.size());
	//finding weights for inserted reference value
	//current insertion is garbage code
	int bestIndex;
	long double bestAvg = 0.0;
	long double bestVar = 0.0;
	int insertIndex;

	//std::complex<long double> bestReference;
	//std::complex<long double> bestReferenceVal;
	std::vector<std::complex<long double>> referencesTempBest;
	std::vector<std::complex<long double>> referenceValsTempBest;
	std::vector<std::vector<long double>> weightsTemp(material_list.size());
	std::vector<long double> referencesReal(references.size());
	std::vector<std::complex<long double>> referencesTemp(references.size());
	std::vector<std::complex<long double>> referenceVals(references.size());
	std::vector < std::complex<long double>> gradient(references.size());

	bool maxVar = true;
	bool avgVar = true;
	long double varAvg;



	while (whileCond) {
		std::cout << "_____________________________________________________________Pass in the for loop_______________________________________________\n";
		if (toggle == false) {
			whileCond = false;
		}
		//need to get rid of qoi_list (not neccessary for kriging), but I need to first get rid of its use in sensitivity_to_epsr method.
		
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

		for (int i = 0; i < references.size(); ++i) {
			std::cout << "references[i]: " << references[i] << " ref Index[i]: " << refIndex[i] << std::endl;

		}

		std::cout << "ref size: " << references.size() << std::endl;

		//iterating through references now and get the qois
		for (int i = 0; i < references.size(); ++i) {
			//if (HOPS_splitting[i].size() == 0) continue;
			in_file_name = file_name;

			//////////////////////////////Frequency Perturbation///////////////////////////////////////////////////////////////////////////////////////////////
			//std::string command = "..\\file_input_converter ..\\reference_files_freq_low_HOPS\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//std::string command = "..\\file_input_converter ..\\reference_files_freq\\ref_3.in ..\\exampleFiles\\" + in_file_name + "\\";
			//---------------------------------------------------------------------------------------------------------------------------------------------------

			/////////////////////////////Material Perturbation///////////////////////////////////////////////////////////////////////////////
			//std::string command = "..\\file_input_converter ..\\reference_files_mat\\ref_0.in ..\\exampleFiles\\" + in_file_name + "\\";
			
			//--------------------------------------------------------------------------------------------------------------------------------------------------


			//setting so that new passes in the while loop only run solve for new points
			//only works if starting with 3 points
			if (references.size() == 3) {			
				std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange3\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
				std::system(command.c_str());
				std::cout << "Testing for reference number: " << refIndex[i] << std::endl;

				referenceVals[i] = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i], gradient[i]);
			}
			else {
				std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange3\\ref_" + std::to_string(refIndex[insertIndex]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
				std::system(command.c_str());
				std::cout << "Testing for reference number: " << refIndex[insertIndex] << std::endl;

				std::complex<long double> val;

				//insert a zero before calling sensitivity so that gradient[i] is in the correct place
				insert_complex(gradient, insertIndex, 0.0);

				val = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[insertIndex], qoi_list[insertIndex], references[insertIndex], gradient[insertIndex]);
				insert_complex(referenceVals, insertIndex, val);
				
				break;
			}
		
			////this is for even split
			//std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange3\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//std::system(command.c_str());
			//std::cout << "Testing for reference number: " << refIndex[i] << std::endl;

			//referenceVals[i] = Kriging::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i], gradient[i]);


		}

		for (int i = 0; i < references.size(); ++i) {
			std::cout << "references: " << references[i].real() << "  referenceVals[i]: " << referenceVals[i].real() << "  gradient[i]: " << gradient[i] << std::endl;
		}
		
		//adding references point based on error

		long double loError = 0.0;
		long double hiError = 0.0;

		std::vector<long double> loErrorVec;
		std::vector<long double> hiErrorVec;

		//subdivide highest error
		int hiErrorIndex = 0;
		long double hiErrorCheck = 0.0;


		//val is the endpoint of HOPS_splitting
		int val = HOPS_splitting.size() - 1;
		std::vector<std::complex<long double>> referencesNew = references;
		std::vector<int> refIndexNew = refIndex;
		
		//calculate variogram_______________________________________________________________

		//main params for variogram (this and fit method)
		long double nBins = 10;
		long double maxDist = abs(references[references.size() - 1].real() - references[0].real());// / 2.0;

		//std::cout << "check 3------------------------\n";
		//long double binTol = maxDist / nBins;
		std::vector<long double> edges;
		std::vector<long double> edgesGrad;
		std::vector<long double> edgesValGrad;
		std::vector<long double> edgesGradVal;


		for (int i = 0; i <= nBins; i++) {
			edges.push_back(i * maxDist / nBins);
			edgesValGrad.push_back(i * maxDist / nBins);
			edgesGrad.push_back(i * maxDist / nBins);
			edgesGradVal.push_back(i * maxDist / nBins);

		}


		//linspace for the x values in the variogramFit
		int fitSize = 15000;
		std::vector<long double> d;
		for (int i = 0; i <= fitSize; ++i) {
			d.push_back(i * maxDist / fitSize);
		}

		//build matrix of distance values
		std::vector<std::vector<std::complex<long double>>> dMat = buildDMatrix2(references, referenceVals, referenceVals);
		std::vector<std::vector<std::complex<long double>>> dMatGrad = buildDMatrix2(references, gradient, gradient);
		std::vector<std::vector<std::complex<long double>>> dMatGradVal = buildDMatrix2(references, gradient, referenceVals); // these are approx 0
		std::vector<std::vector<std::complex<long double>>> dMatValGrad = buildDMatrix2(references, referenceVals, gradient); //


		std::cout << "buildDMatrix2 works\n";


		//maybe make variogram for real and imag gradient parts

		//build variogram matrix and solve for empirical variogram
		std::vector<long double> variogram = buildVariogram(dMat, references, edges);
		std::vector<long double> variogramGrad = buildVariogram(dMatGrad, references, edges);
		std::vector<long double> variogramGradVal = buildVariogram(dMatGradVal, references, edges);
		std::vector<long double> variogramValGrad = buildVariogram(dMatValGrad, references, edges);

		std::vector<std::complex<long double>> variogramComplex = buildVariogramComplex(dMat, references, edges);
		std::vector<std::complex<long double>> variogramGradComplex = buildVariogramComplex(dMatGrad, references, edges);

		std::cout << "buildVariogram works\n";

		//trim any nan values
		trimVariogram(variogram, edges);
		trimVariogram(variogramGrad, edgesGrad);
		trimVariogram(variogramGradVal, edgesGradVal);
		trimVariogram(variogramValGrad, edgesValGrad);

		trimVariogramComplex(variogramComplex, edges);
		trimVariogramComplex(variogramGradComplex, edgesGrad);


		//fit model to empirical variogram
		std::vector<long double> variogramFit = fitVariogramModel(variogram, edges, d, references.size());
		std::vector<long double> variogramFitGrad = fitVariogramModel(variogramGrad, edgesGrad, d, references.size());
		std::vector<long double> variogramFitGradVal = fitVariogramModel(variogramGradVal, edgesGradVal, d, references.size());
		std::vector<long double> variogramFitValGrad = fitVariogramModel(variogramValGrad, edgesValGrad, d, references.size());

		std::vector<long double> variogramComplexReal(variogramComplex.size());
		std::vector<long double> variogramComplexImag(variogramComplex.size());
		std::vector<long double> variogramGradComplexReal(variogramGradComplex.size());
		std::vector<long double> variogramGradComplexImag(variogramGradComplex.size());
		for (int i = 0; i < variogram.size(); ++i) {
			variogramComplexReal[i] = variogramComplex[i].real();
			variogramComplexImag[i] = variogramComplex[i].imag();
			variogramGradComplexReal[i] = variogramGradComplex[i].real();
			variogramGradComplexImag[i] = variogramGradComplex[i].imag();
		}

		std::vector<long double> variogramFitReal = fitVariogramModel(variogramComplexReal, edges, d, references.size());
		std::vector<long double> variogramFitImag = fitVariogramModel(variogramComplexImag, edges, d, references.size());
		std::vector<long double> variogramFitGradReal = fitVariogramModel(variogramGradComplexReal, edges, d, references.size());
		std::vector<long double> variogramFitGradImag = fitVariogramModel(variogramGradComplexImag, edges, d, references.size());


		//std::vector<long double> variogramFit = fitVariogram(variogram, edges, d, references.size());
		//std::vector<long double> variogramFitGrad = fitVariogram(variogramGrad, edgesGrad, d, references.size());
		//std::vector<long double> variogramFitGradVal = fitVariogram(variogramGradVal, edgesGradVal, d, references.size());
		//std::vector<long double> variogramFitValGrad = fitVariogram(variogramValGrad, edgesValGrad, d, references.size());

		//std::vector<long double> variogramFitGradVal = variogramFitValGrad;

		//calculate variogram derivatives
		std::vector<long double> variogramFitFirstD = firstDerivative(variogramFit, d);
		std::vector<long double> variogramFitSecondD = firstDerivative(variogramFitFirstD, d);


		//------------------------------- output empirical variogram and fitted variogram for checking----------------------------
		
		std::string outvar = "../iofiles/output/taylorkriging/sweepGauss/covariogram" + std::to_string(iterations + 3) + ".txt";
		std::ofstream var_out(outvar);
		for (int i = 0; i < edges.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			var_out << edges[i] << "\t" << variogramValGrad[i] << std::endl;
		}
		var_out.close();

		std::string outvarfit = "../iofiles/output/taylorkriging/sweepGauss/covariogramFitExp" + std::to_string(iterations + 3) + ".txt";
		std::ofstream var_fit_out(outvarfit);
		for (int i = 0; i < d.size(); ++i) {
			//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
			//	continue;
			//}
			var_fit_out << d[i] << "\t" << variogramFitValGrad[i] << std::endl;
		}
		var_fit_out.close();

		///----------------------------------------------------------------------------------------

		varAvg = 0.0;
		//---------------------  Kriging System Solution --------------------------------------------

		long double xAvg;
		for (int i = 0; i < material_list.size(); ++i) {

			taylorKriging(variogramFitGrad, d, references, weightsGrad[i], material_list[i].real());
			

			//taylorKrigingHigherOrder(std::vector<long double> variogramFit, std::vector<long double> variogramFitFirstD, std::vector<long double> variogramFitSecondD, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal, std::vector<long double> variogramFitGrad)
			//taylorKriging(variogramFit, d, references, weights[i], material_list[i].real());
			//taylorKrigingHigherOrder(variogramFit, variogramFitFirstD, variogramFitSecondD, d, references, weightsHigher[i], material_list[i].real(), variogramFitGrad);
			xAvg = taylorKrigingHigherOrder2(variogramFit, variogramFitValGrad, variogramFitGradVal, variogramFitGrad, d, references, weightsHigher[i], material_list[i].real());
			

			//reconstruction[i] = krigSum(weights[i], referenceVals);
			reconstruction[i] = krigSumHigherOrder(weightsHigher[i], referenceVals, gradient, material_list[i].real(), xAvg);
			
			reconstructionGrad[i] = krigSum(weightsGrad[i], gradient);
		}

		std::cout << "out of kriging system loop\n";

		std::cout << "referencesFull.size(): " << referencesFull.size() << std::endl;

		bestVar = 0;

		long double matMean = 4.5;
		long double matStd = 1.0;
		//------------------------------------------------------------------------------------------------------------------

		//doesnt work
		////------------------------Iterative update to Cov function------------------------------------------------------
		//bool covUpdate = true;

		//while (covUpdate) {
		//	//build matrix of distance values
		//	std::vector<std::vector<std::complex<long double>>> dMatUpdate = buildDMatrix2(material_list, reconstruction, reconstruction);
		//	std::vector<std::vector<std::complex<long double>>> dMatGradUpdate = buildDMatrix2(material_list, reconstructionGrad, reconstructionGrad);


		//	//build variogram matrix and solve for empirical variogram
		//	std::vector<long double> variogramUpdate = buildVariogram(dMatUpdate, material_list, edges);
		//	std::vector<long double> variogramGradUpdate = buildVariogram(dMatGradUpdate, material_list, edges);

		//	//trim any nan values
		//	trimVariogram(variogramUpdate, edges);
		//	trimVariogram(variogramGradUpdate, edgesGrad);


		//	//fit model to empirical variogram
		//	std::vector<long double> variogramFitUpdate = fitVariogramModel(variogramUpdate, edges, d, material_list.size());
		//	std::vector<long double> variogramFitGradUpdate = fitVariogramModel(variogramGradUpdate, edgesGrad, d, material_list.size());


		//	for (int i = 0; i < material_list.size(); ++i) {
		//		//taylorKrigingHigherOrder(std::vector<long double> variogramFit, std::vector<long double> variogramFitFirstD, std::vector<long double> variogramFitSecondD, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<long double>& weights, long double xVal, std::vector<long double> variogramFitGrad)
		//		//taylorKriging(variogramFit, d, references, weights[i], material_list[i].real());
		//		//taylorKrigingHigherOrder(variogramFit, variogramFitFirstD, variogramFitSecondD, d, references, weightsHigher[i], material_list[i].real(), variogramFitGrad);
		//		xAvg = taylorKrigingHigherOrder2(variogramFitUpdate, variogramFitValGrad, variogramFitGradVal, variogramFitGradUpdate, d, references, weightsHigher[i], material_list[i].real());


		//		//reconstruction[i] = krigSum(weights[i], referenceVals);
		//		reconstructionUpdate[i] = krigSumHigherOrder(weightsHigher[i], referenceVals, gradient, material_list[i].real(), xAvg);

		//	}

		//	std::string outF = "../ioFiles/output/taylorKriging/sweepHigherUpdate3/qoi_kriging" + std::to_string(iterations + 2) + ".txt";
		//	std::ofstream qoi_dist_out(outF);
		//	for (int i = 0; i < material_list.size(); ++i) {
		//		//if (material_list[i].real() < 2.0 || material_list[i].real() > 7.0) {
		//		//	continue;
		//		//}
		//		if (i == 0)
		//			std::cout << std::setprecision(12) << "precision check" << reconstruction[i].real() << std::endl;
		//		qoi_dist_out << std::setprecision(12) << material_list[i].real() << " " << material_list[i].imag() << " " << reconstructionUpdate[i].real() << " " << reconstructionUpdate[i].imag() << " " << reconstructionGrad[i].real() << " " << reconstructionGrad[i].imag() << std::endl;
		//	}
		//	qoi_dist_out.close();
		//	//qoi_dist_out.close();
		//	break;


		//}
		////--------------------------------------------END Iterative Update-----------------------------------------------





		//find max variance for next iteration--------------------------------------------------------------------------------------------------
		for (int i = 0; i < referencesFull.size(); ++i) {
			long double variance = 0.0;
			long double avg = 0;
			long double dist = 0;

			//skip values of referencesFull that are already contained in references
			if (std::find(references.begin(), references.end(), referencesFull[i]) != references.end()) {
				std::cout << "references[" << i << "] found \n";
				continue;
			}

			//skip values of referencesFull outside the range
			if (i < loRefIndex || i >  hiRefIndex)
				continue;

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

			//for (int j = 0; j < referencesTemp.size(); ++j) {
			//	std::cout << "referencesTemp[i]: " << referencesTemp[j] << std::endl;
			//}

			//std::cout << "weight index: " << weightIndex << std::endl;

			//std::cout << "referencesFull[i](added ref): " << referencesFull[i] << "referencesTemp[weightIndex](should be same): " << referencesTemp[weightIndex] << std::endl;
			//long double variance = 0;
			//find new references variance:
			//for (int i = 0; i < referencesTemp.size(); ++i) {
				//variance += (referencesTemp[i].real() - 4.5) * (referencesTemp[i].real() - 4.5);
			//}
			//variance /= referencesTemp.size();

			
			//---------variogram---------------

			std::vector<long double> edges;
			for (int j = 0; j <= nBins; j++) {
				edges.push_back(j * maxDist / nBins);
			}

			//linspace for the x values in the variogramFit
			int fitSize = 15000;
			std::vector<long double> d;
			for (int j = 0; j <= fitSize; ++j) {
				d.push_back(j * maxDist / fitSize);
			}

			//build referenceVals from previous reconstruction
			//not sure if this works.  Should find the closest material parameter in material list to referencesTemp.
			//it should then fill referenceValsTemp with the cooresponding index from the kriging reconstruction

			//build matrix of distance values
			std::vector<std::vector<std::complex<long double>>> dMat = buildDMatrix2(references, referenceVals, referenceVals);

			//build variogram matrix and solve for empirical variogram
			std::vector<long double> variogram = buildVariogram(dMat, references, edges);

			//trim any nan values
			trimVariogram(variogram, edges);

			//fit model to empirical variogram
			std::vector<long double> variogramFit = fitVariogramModel(variogram, edges, d, references.size());


			std::vector<std::complex<long double>> referenceValsTemp(referencesTemp.size());
			for (int j = 0; j < referenceValsTemp.size(); j++) {
				//maybe just set the added reference value to reconstruction (others already calculated explicitly)
				int index = find_index(material_list, referencesTemp[j].real());
				referenceValsTemp[j] = reconstruction[index];
				//std::cout << "referencesTemp[j](value looking to match): " << referencesTemp[j] << " material_list[index](index should match): " << material_list[index] << std::endl;
			}

			//find the weights with referencesFull[i] as the new test point x0
			//taylorKriging(variogramFit, d, references, weightsTemp[i], referencesFull[i].real());
			xAvg = taylorKrigingHigherOrder2(variogramFit, variogramFitValGrad, variogramFitGradVal, variogramFitGrad, d, references, weightsTemp[i], referencesFull[i].real());


			//calculate the variance of the added referencesFull point
			variance = calcVariance2(variogramFit, variogramFitGrad, d, references, weightsTemp[i], referencesFull[i].real());
			varAvg += variance * variance;

			//std::cout << "variance: " << variance << std::endl;

			long double integral = Kriging::trap_integral(referencesFull[i].real(), matMean, matStd, material_list[1].real() - material_list[0].real());

			//std::cout << "   integral: " << integral << std::endl;

			//weighting function
			//variance *= integral;


			////_____________________________OLD METHOD for calculating variance____________________________________________________
			//std::vector<std::vector<long double>> dMat = buildDMatrix(referencesTemp, referenceValsTemp);

			//std::vector<long double> variogram = buildVariogram(dMat, referencesTemp, edges);
			//trimVariogram(variogram, edges);
			//std::vector<long double> variogramFit = fitVariogram(variogram, edges, d);
			//

			//for (int j = 0; j < material_list.size(); ++j) {
			//	kriging(variogramFit, d, referencesTemp, weightsTemp[j], material_list[j].real());
			//	
			//	//relative distance between material_list[j] entry and referencesFull[i].
			//	//the greater the distance, the less each weight. ie weight corresponds to 1/d.  multiplying by (relative) distance should fix?
			//	//to compensate, 
			//	//dist = abs(abs(material_list[j].real() - referencesFull[i].real()) / abs(referencesFull[i].real() + material_list[j].real()));
			//	variance += calcVariance(variogramFit, d, referencesTemp, weightsTemp[j], material_list[j].real());

			//	//avg += abs(weightsTemp[j][weightIndex] * dist);
			//	//std::cout << avg << std::endl;
			// 
			//	//std::cout << "weightsTemp[j][weightIndex]: " << weightsTemp[j][weightIndex] << std::endl;
			//}
		

			//std::cout << "var: " << variance << " material_list.size(): " << material_list.size() << std::endl;
			//variance /= material_list.size();
			//avg = avg / material_list.size();
			////____________________________________________________________________________________

			std::cout << "referencesFull[i]: " << referencesFull[i] << "  variance: " << variance << std::endl;

			//std::cout << "avg: " << avg << std::endl;

			//if (avg >= bestAvg) {
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

		long double z0 = 0.0;
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
		std::string outF = "../ioFiles/output/taylorKriging/sweepGaussHigher3/qoi_kriging" + std::to_string(iterations + 2) + ".txt";
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

long double Kriging::normalFunction(long double x, long double mean, long double std) {
	long double arg = (x - mean) / std;
	long double exponential = exp(-0.5 * arg * arg);
	long double integral = 1.0 / (std * sqrt(2.0 * 3.14159)) * exponential;
	return integral;
}

long double Kriging::trap_integral(long double intervalMid, long double mean, long double std, long double step) {
	long double count = 50.0;
	long double intervalBeg = intervalMid - count * step / 2.0;
	long double intervalEnd = intervalMid + count * step / 2.0;
	long double integral = 0.5 * (Kriging::normalFunction(intervalBeg, mean, std) + Kriging::normalFunction(intervalEnd, mean, std));
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
	std::vector<std::vector<std::complex<long double>>> material_list;
	//std::ifstream materials_in("../monte_carlo_files_highdim/materials_list.txt");
	std::ifstream materials_in("../ioFiles/input/frequencies_list1.txt");
	std::string line;
	

	//get the monte carlo points
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		std::vector<std::complex<long double>> material_list_temp;
		for (int mati = 0; mati < data.size(); mati+=1) {
			// b/c the data is irrelevant for air and PML layer, set those to zero
			if (mati < 64)
				material_list_temp.push_back(std::complex<long double>(std::stod(data[mati]),0.0));
			else
				material_list_temp.push_back(0.0);
		}
		material_list.push_back(material_list_temp);
	}
	materials_in.close();
	
	////generate references in c++
	//long double sample_start = 2.56 - .21; //3 std
	//long double sample_end = 2.56 + .21;//3 std
	//int num_refs = 1;
	//std::vector<std::vector<std::complex<long double>>> references(num_refs);
	//long double increment = (sample_end - sample_start) / (long double(num_refs-1.0));
	//for (int refs = 0; refs < num_refs; ++refs) {
	//	long double val = sample_start + refs*increment;
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
	std::vector<std::vector<std::complex<long double>>> references;
	std::ifstream references_in("../monte_carlo_files_highdim/references_highdim.txt");
	while (std::getline(references_in, line)) {
		auto data = functions::split(line, ' ');
		std::vector<std::complex<long double>> material_list_temp;
		for (int mati = 0; mati < data.size(); mati += 1) {
			// b/c the data is irrelevant for air and PML layer, set those to zero
			if (mati < 64)
				material_list_temp.push_back(std::complex<long double>(std::stod(data[mati]), 0.0));
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

	
	
	std::vector<std::vector<std::vector<std::complex<long double>>>> HOPS_splitting(references.size());
	
	int mat_counter = 0;
	

	for (int j = 0; j < material_list.size(); ++j) {
		//find the closest reference point
		long double min_dist = 1.0e15;
		int mindex = 0;
		for (int i = 0; i < references.size(); ++i) {
			
				long double dist = norm(references[i] - material_list[j]);
				if (dist < min_dist) { 
					mindex = i; 
					min_dist = dist;
				}
			
		}
		
			HOPS_splitting[mindex].push_back(material_list[j]);
			mat_counter++;
			
		
	}
	std::vector<std::vector<std::complex<long double>>> qoi_list(references.size());
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
void Kriging::sensitivity_to_multi_epsr(std::string & file_name, std::vector<std::vector<std::complex<long double>>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::vector<std::complex<long double>>& reference_values)
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
	std::complex<long double> reference_qoi = get_QoI(dom_forward, dom_adjoint);
	std::cout << "Reference qoi: " << reference_qoi << std::endl;

	long double k0 = dom_forward.scatter1.K0[1];
	//get the remainder of the estimate
	std::vector<std::complex<long double>> dQoI(epsr_list.size(), 0.0);
	
	
	//loop through the stored Bint pieces
	for (auto e = dom_forward.elements.begin(); e != dom_forward.elements.end(); ++e) {

		if (e->materials.region == 1) {
			//if in scatterer
			for (int i_unkn = e->unknownsStart; i_unkn <= e->unknownsEnd; ++i_unkn) {
				int icon = abs(dom_forward.vectorD[i_unkn]);
				for (int epsr_val = 0; epsr_val < epsr_list.size(); ++epsr_val) {
					std::complex<long double> eps_diff = epsr_list[epsr_val][e->index-1] / e->materials.epsr_list[0][0] - 1.0;
					dQoI[epsr_val] += eps_diff*e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
				}
				//	dQoI1 += eps_diff * e->cGr_eps_el[icon] * std::conj(dom_adjoint.cAlpha[icon]);
				for (int j_unkn = e->unknownsStart; j_unkn <= e->unknownsEnd; ++j_unkn) {
					int jcon = abs(dom_adjoint.vectorD[j_unkn]);
					int i0 = i_unkn - e->unknownsStart + 1;
					int j0 = j_unkn - e->unknownsStart + 1;
					
					for (int epsr_val = 0; epsr_val < epsr_list.size(); ++epsr_val) {
						std::complex<long double> eps_diff = epsr_list[epsr_val][e->index - 1] / e->materials.epsr_list[0][0] - 1.0;
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