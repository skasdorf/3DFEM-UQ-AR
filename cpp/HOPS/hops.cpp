#include "hops.h"
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;
#include "src/stdafx.h"
//#include <stdlib.h>
//#include <stdio.h>
#include "src/interpolation.h"
#include "src/ap.h"

static bool complexComparitor(std::complex<long double> a, std::complex<long double> b) {
	return real(a) <= real(b);
}

std::complex<long double> get_QoI2(Domain& dom_forward, Domain& dom_adjoint) {
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

//std::vector<std::vector<long double>> buildDMatrix(std::vector<std::complex<long double>> x, std::vector<std::complex<long double>> y) {
//	std::vector<std::vector<long double>> dMat(x.size());
//	int index = 0;
//	for (int i = 0; i < x.size(); ++i) {
//		long double RCSi = sqrt(y[i].real() * y[i].real() + y[i].imag() * y[i].imag()) / 4 / 3.14159;
//		for (int j = 0; j < x.size(); ++j) {
//			long double RCSj = sqrt(y[j].real() * y[j].real() + y[j].imag() * y[j].imag()) / 4 / 3.14159;
//			dMat[index].push_back(i);
//			dMat[index].push_back(j);
//			dMat[index].push_back(abs(x[i].real() - x[j].real()));
//			dMat[index].push_back((RCSi - RCSj) * (RCSi - RCSj));
//			index++;
//		}
//	}
//	return dMat;
//}
//
//std::vector<long double> buildVariogram(std::vector<std::vector<long double>> dMat, std::vector<std::complex<long double>> x, std::vector<long double>& edges) {
//
//	std::vector<long double> binCounts(edges.size()-1);
//	std::vector<long double> binIdx(x.size()*x.size());
//	for (int i = 0; i < edges.size(); ++i) {
//		for (int j = 0; j < binIdx.size(); ++j) {
//			if (dMat[j][2] < edges[i + 1] && dMat[j][2] >= edges[i]) {
//				binCounts[i] += 1.0;
//				binIdx[j] = i;
//			}
//		}
//	}
//
//	std::vector<long double> variogram(edges.size());
//	for (int i = 0; i < edges.size() - 1.0; ++i) {
//		for (int j = 0; j < binIdx.size(); ++j) {
//			if (binCounts[i] != 0) {
//				if (binIdx[j] == i) {
//					variogram[int(binIdx[j])] += dMat[j][3];
//				}
//			}
//		}
//	}
//	for (int i = 0; i < variogram.size(); ++i) {
//		variogram[i] = variogram[i] / (2.0 * binCounts[i]);
//	}
//	return variogram;
//	
//}
//
//std::vector<long double> kriging(std::vector<long double> variogramFit, std::vector<long double> d, std::vector<std::complex<long double>> xSample, std::vector<std::complex<long double>> ySample) {
//	const int size = xSample.size();
//	//Matrix <long double, Dynamic, Dynamic> krigMat;
//	MatrixXd krigMat(size, size);
//	for (int i = 0; i < d.size(); ++i) {
//		for (int j = 0; j < d.size(); ++j) {
//			if (i == d.size() - 1) {
//				krigMat.setConstant(i, j, 1.0);
//			}
//			if (j == d.size() - 1) {
//				krigMat.setConstant(i, j, 1.0);
//			}
//			long double dist = abs(xSample[i] - xSample[j]);
//			auto it = std::lower_bound(d.begin(), d.end(), dist);
//			int index = it - d.begin();
//			krigMat.setConstant(i, j, variogramFit[index]);
//		}
//	}
//	krigMat.setConstant(d.size(), d.size(), 0.0);
//
//
//
//
//}
//
//void trimVariogram(std::vector<long double>& variogram, std::vector<long double>& edges) {
//	int idx = 0;
//	std::vector<long double> variogramNew;
//	std::vector<long double> edgesNew;
//	for (int i = 0; i < variogram.size(); i++) {
//		if (!isnan(variogram[i])) {
//			variogramNew[idx] = variogram[i];
//			edgesNew[idx] = edges[i];
//		}
//	}
//	variogram = variogramNew;
//	edges = edgesNew;
//
//}
//
//std::vector<long double> fitVariogram(std::vector<long double> variogram, std::vector<long double> edges, std::vector<long double> d) {
//	long double xSum = 0.0;
//	long double ySum = 0.0;
//	long double x2Sum = 0.0;
//	long double xySum = 0.0;
//	for (int i = 0; i < variogram.size(); ++i) {
//		xSum += edges[i];
//		ySum += variogram[i];
//		x2Sum += x2Sum + (edges[i] * edges[i]);
//		xySum += xySum + (edges[i] * variogram[i]);
//	}
//	long double slope = (edges.size() * xySum - xSum * ySum) / (edges.size() * x2Sum - xSum * xSum);
//	std::vector<long double> fitVariogram(d.size());
//
//
//	for (int i = 0; i < d.size(); ++i) {
//		fitVariogram[i] = slope * d[i];
//	}
//	return fitVariogram;
//}

std::complex<long double> HOPS::sensitivity_to_epsr(std::string & file_name, std::vector<std::complex<long double>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::complex<long double> referenceFreq, std::complex<long double>& gradient)
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

	long double kReal = 2 * 3.14159 / (2.99E8 / referenceFreq.real());
	std::complex<long double> referenceK = (kReal, 0.0);
	//now have the base QoI modifier, need to scale by the change in perturbed material parameter
	for (int i = 0; i < epsr_list.size(); ++i) {
		/////////////////////////////////////////Frequency Perturbation/////////////////////////////////////////////
		//long double k0_list = epsr_list[i].real() * 2.0 * 3.14159 / 3e8;
		//eps_diff = -1.0*k0_list + k0;
		//eps_diff = (k0_list / k0 - 1.0) * k0;
		//---------------------------------------------------------------------------------------------------------

		///////////////////////////////////////Material Perturbation//////////////////////////////////////////////
		eps_diff = epsr_list[i] - referenceFreq;
		///-------------------------------------------------------------------------------------------------------

		std::complex<long double> dQoIMag = (real(dQoI) * real(dQoI) + imag(dQoI) * imag(dQoI));
		std::complex<long double> eps_diffMag = (real(eps_diff) * real(eps_diff) + imag(eps_diff) * imag(eps_diff));

		std::complex<long double> eps_diff_scaled = (real(eps_diff) * real(eps_diff) + imag(eps_diff) * imag(eps_diff));
		std::complex<long double> reference_qoi_scaled = (real(reference_qoi) * real(reference_qoi) + imag(reference_qoi) * imag(reference_qoi));

		updated_qoi.push_back(reference_qoi + eps_diff * dQoI / eps_diffMag);// +eps_diff * eps_diff * dQoI2 / 2.0);

		std::cout << "gradient magnitude: " << dQoIMag << std::endl;
		//std::cout << reference_qoi + eps_diff * dQoI << std::endl;
		//std::cout << "qoi: " << reference_qoi + eps_diff * dQoI << std::endl;
		
	}

	std::cout << "gradient: " << dQoI << "   referenceFreq: " << referenceFreq << "   referenceQoi: " << reference_qoi << "\n";

	gradient = dQoI;
	return reference_qoi;
}

std::complex<long double> HOPS::get_QoI(Domain & dom_forward, Domain & dom_adjoint)
{
	std::complex<long double> qoi = 0.0;
	for (int i = 0; i < dom_forward.scatter1.cGr[1].size(); ++i) {
		qoi += dom_forward.scatter1.cGr[1][i] * std::conj(dom_adjoint.cAlpha[i]);
	}
	return qoi;
}
void HOPS::monte_carlo_instance(std::string & file_name)
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



	//this works for just output value///
	//mesh inputs
	//std::string mesh_name = file_name;
	////get the forward solution
	//Domain dom_forward(mesh_name);
	//run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, false, 0.0);
	////get RHS
	//Domain dom_adjoint_rhs(mesh_name);
	//run::_RHS_compute(dom_adjoint_rhs, mesh_name, plot, true, higher_order, check_results);
	////get QoI
	//std::complex<long double> qoi = get_QoI2(dom_forward, dom_adjoint_rhs);

	//instead for value and gradient we need to call sensitivity to epsr
	std::string mesh_name = file_name;



}


//dual here just means QoI solution and QoI gradient
//probably should come up with a better name
//dual typically means some sort of adjoint transformation problem
void HOPS::monte_carlo_dual(std::string& file_name)
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



	//this works for just output value///
	// 
	//mesh inputs
	//std::string mesh_name = file_name;
	////get the forward solution
	//Domain dom_forward(mesh_name);
	//run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, false, 0.0);
	////get RHS
	//Domain dom_adjoint_rhs(mesh_name);
	//run::_RHS_compute(dom_adjoint_rhs, mesh_name, plot, true, higher_order, check_results);
	////get QoI
	//std::complex<long double> qoi = get_QoI2(dom_forward, dom_adjoint_rhs);
	// 
	//instead for value and gradient we need to call sensitivity to epsr



	//materials list
	std::vector<std::complex<long double>> material_list;
	std::ifstream materials_in("../ioFiles/input/materials_list5.txt");
	std::string line;
	//std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		//material_list.push_back(std::complex<long double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
		material_list.push_back(std::complex<long double>(std::stod(data[0]), std::stod(data[1])));
		//std::cout << std::stod(data[0]) << std::endl;
	}
	materials_in.close();

	//sort(material_list.begin(), material_list.end(), complexComparitor);

	std::vector<std::vector<std::complex<long double>>> HOPS_splitting(material_list.size());
	std::vector<std::vector<std::complex<long double>>> qoi_list(material_list.size());
	//sort(material_list.begin(), material_list.end(), complexComparitor);

	std::vector<std::complex<long double>> gradient(material_list.size());
	std::vector<std::complex<long double>> referenceVals(material_list.size());
	std::vector<std::complex<long double>> references(material_list.size());

	//for (int i = 0; i < material_list.size(); ++i) {
	//	HOPS_splitting[i] = material_list[i];
	//}


	//std::cout << "made it here 1 ...........................\n";

	//iterating through references now and get the qois
	for (int i = 0; i < 1000; ++i) {
		//if (HOPS_splitting[i].size() == 0) continue;
		std::string in_file_name = file_name;

		//////////////////////////////Frequency Perturbation///////////////////////////////////////////////////////////////////////////////////////////////
		//std::string command = "..\\file_input_converter ..\\reference_files_freq_low_HOPS\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
		//std::string command = "..\\file_input_converter ..\\reference_files_freq\\ref_3.in ..\\exampleFiles\\" + in_file_name + "\\";
		//---------------------------------------------------------------------------------------------------------------------------------------------------

		//std::cout << "made it here 2 ...........................\n";

		/////////////////////////////Material Perturbation///////////////////////////////////////////////////////////////////////////////
		//std::string command = "..\\file_input_converter ..\\reference_files_mat\\ref_0.in ..\\exampleFiles\\" + in_file_name + "\\";
		//std::string command = "..\\file_input_converter ..\\reference_files_mat_newCarlo4\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
		std::string command = "..\\file_input_converter ..\\..\\..\\FEM-Matlab\\Meshing_of_Spheres_dl\\Meshing_of_Spheres\\monte_carlo_files_mat\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + in_file_name + "\\";

		//--------------------------------------------------------------------------------------------------------------------------------------------------

		std::cout << "material_list[i]: " << material_list[i] << std::endl;

		std::system(command.c_str());
		std::cout << "Testing for reference number: " << i << std::endl;
		referenceVals[i] = HOPS::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], material_list[i], gradient[i]);

		//std::cout << "made it here 4 ...........................\n";

		//output one time
		//output results to file
		std::ofstream qoi_dist_out("../ioFiles/output/MonteCarlo/MC5.txt", std::ios::out | std::ios::app);
		qoi_dist_out << std::setprecision(12) << material_list[i].real() << " " << material_list[i].imag() << " " << referenceVals[i].real() << " " << referenceVals[i].imag() << " " << gradient[i].real() << " " << gradient[i].imag() << std::endl;


		qoi_dist_out.close();

	}






	//output one time
//output results to file
//std::ofstream qoi_dist_out("../ioFiles/output/sweep/weight/qoi_HOPS_materialAR0.txt");
//int index = 0;
//for (int i = 0; i < qoi_list.size(); i++) {
//	//need to remove the max/min test values for each of the batches
//	for (int j = 0; j < qoi_list[i].size()-2; ++j) {
//		if (index > HOPS_splitting[index].size()) index++;
//		//std::cout << "i: " << i << " j: " << j << std::endl;
//		qoi_dist_out << HOPS_splitting[i][j].real() << " " << HOPS_splitting[i][j].imag() << " " << qoi_list[i][j].real() << " " << qoi_list[i][j].imag() << std::endl;
//	}
//}
//qoi_dist_out.close();


}




//different method////
void HOPS::multi_HOPS_epsr(std::string & file_name)
{


	//load in the list of random materials being tested
	std::vector<std::complex<long double>> material_list;
	std::ifstream materials_in("../ioFiles/input/materials_list4.txt");
	std::string line;
	//std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, ' ');
		//material_list.push_back(std::complex<long double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
		material_list.push_back(std::complex<long double>(std::stod(data[0]), std::stod(data[1])));
	}
	materials_in.close();

	sort(material_list.begin(), material_list.end(), complexComparitor);

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

	std::cout << "material_list size: " << material_list.size() << std::endl;

	//acceptable error
	long double delta = 0.0008;

	//convergence condition
	long double stdLim = 0.0001;

	//iteration limit
	int itLim = 20;

	//material_list.txt
	//long double matStd = 1.0;
	//long double matMean = 4.5;

	//material_list2.txt
	long double matStd = 1.0;
	long double matMean = 4.5;



	std::vector<long double> stdDevRCSList;
	std::vector<std::complex<long double>> stdDevList;


	///////////////////////Material Perturbation/////////////////////////////////////////////////////////////////////

	//old frequency and material range
	//std::vector<std::complex<long double>> referencesFull = { {1.500000, -2.000000}, { 1.600000, -2.000000 }, { 1.700000, -2.000000 }, { 1.800000, -2.000000 }, { 1.900000, -2.000000 }, { 2.000000, -2.000000 }, { 2.100000, -2.000000 }, { 2.200000, -2.000000 }, { 2.300000, -2.000000 }, { 2.400000, -2.000000 }, { 2.500000, -2.000000 }, { 2.600000, -2.000000 }, { 2.700000, -2.000000 }, { 2.800000, -2.000000 }, { 2.900000, -2.000000 }, { 3.000000, -2.000000 }, { 3.100000, -2.000000 }, { 3.200000, -2.000000 }, { 3.300000, -2.000000 }, { 3.400000, -2.000000 }, { 3.500000, -2.000000 }, { 3.600000, -2.000000 }, { 3.700000, -2.000000 }, { 3.800000, -2.000000 }, { 3.900000, -2.000000 }, { 4.000000, -2.000000 }, { 4.100000, -2.000000 }, { 4.200000, -2.000000 }, { 4.300000, -2.000000 }, { 4.400000, -2.000000 }, { 4.500000, -2.000000 }, { 4.600000, -2.000000 }, { 4.700000, -2.000000 }, { 4.800000, -2.000000 }, { 4.900000, -2.000000 }, { 5.000000, -2.000000 }, { 5.100000, -2.000000 }, { 5.200000, -2.000000 }, { 5.300000, -2.000000 }, { 5.400000, -2.000000 }, { 5.500000, -2.000000 }, { 5.600000, -2.000000 }, { 5.700000, -2.000000 }, { 5.800000, -2.000000 }, { 5.900000, -2.000000 }, { 6.000000, -2.000000 }, { 6.100000, -2.000000 }, { 6.200000, -2.000000 }, { 6.300000, -2.000000 }, { 6.400000, -2.000000 }, { 6.500000, -2.000000 }, { 6.600000, -2.000000 }, { 6.700000, -2.000000 }, { 6.800000, -2.000000 }, { 6.900000, -2.000000 }, { 7.000000, -2.000000 }, { 7.100000, -2.000000 }, { 7.200000, -2.000000 }, { 7.300000, -2.000000 }, { 7.400000, -2.000000 }, { 7.500000, -2.000000 }, { 7.600000, -2.000000 }, { 7.700000, -2.000000 }, { 7.800000, -2.000000 }, { 7.900000, -2.000000 }, { 8.000000, -2.000000 } };

	//new material range
	//std::vector<std::complex<long double>> referencesFull = { {1.000000, -2.000000},{1.100000, -2.000000},{1.200000, -2.000000},{1.300000, -2.000000},{1.400000, -2.000000},{1.500000, -2.000000},{1.600000, -2.000000},{1.700000, -2.000000},{1.800000, -2.000000},{1.900000, -2.000000},{2.000000, -2.000000},{2.100000, -2.000000},{2.200000, -2.000000},{2.300000, -2.000000},{2.400000, -2.000000},{2.500000, -2.000000},{2.600000, -2.000000},{2.700000, -2.000000},{2.800000, -2.000000},{2.900000, -2.000000},{3.000000, -2.000000},{3.100000, -2.000000},{3.200000, -2.000000},{3.300000, -2.000000},{3.400000, -2.000000},{3.500000, -2.000000},{3.600000, -2.000000},{3.700000, -2.000000},{3.800000, -2.000000},{3.900000, -2.000000},{4.000000, -2.000000},{4.100000, -2.000000},{4.200000, -2.000000},{4.300000, -2.000000},{4.400000, -2.000000},{4.500000, -2.000000},{4.600000, -2.000000},{4.700000, -2.000000},{4.800000, -2.000000},{4.900000, -2.000000},{5.000000, -2.000000},{5.100000, -2.000000},{5.200000, -2.000000},{5.300000, -2.000000},{5.400000, -2.000000},{5.500000, -2.000000},{5.600000, -2.000000},{5.700000, -2.000000},{5.800000, -2.000000},{5.900000, -2.000000},{6.000000, -2.000000},{6.100000, -2.000000},{6.200000, -2.000000},{6.300000, -2.000000},{6.400000, -2.000000},{6.500000, -2.000000},{6.600000, -2.000000},{6.700000, -2.000000},{6.800000, -2.000000},{6.900000, -2.000000},{7.000000, -2.000000},{7.100000, -2.000000},{7.200000, -2.000000},{7.300000, -2.000000},{7.400000, -2.000000},{7.500000, -2.000000},{7.600000, -2.000000},{7.700000, -2.000000},{7.800000, -2.000000},{7.900000, -2.000000},{8.000000, -2.000000} };
	//std::vector<std::complex<long double>> referencesFull = { {3.000000, -2.000000},{3.100000, -2.000000},{3.200000, -2.000000},{3.300000, -2.000000},{3.400000, -2.000000},{3.500000, -2.000000},{3.600000, -2.000000},{3.700000, -2.000000},{3.800000, -2.000000},{3.900000, -2.000000},{4.000000, -2.000000},{4.100000, -2.000000},{4.200000, -2.000000},{4.300000, -2.000000},{4.400000, -2.000000},{4.500000, -2.000000},{4.600000, -2.000000},{4.700000, -2.000000},{4.800000, -2.000000},{4.900000, -2.000000},{5.000000, -2.000000},{5.100000, -2.000000},{5.200000, -2.000000},{5.300000, -2.000000},{5.400000, -2.000000},{5.500000, -2.000000},{5.600000, -2.000000},{5.700000, -2.000000},{5.800000, -2.000000},{5.900000, -2.000000},{6.000000, -2.000000},{6.100000, -2.000000},{6.200000, -2.000000},{6.300000, -2.000000},{6.400000, -2.000000},{6.500000, -2.000000},{6.600000, -2.000000},{6.700000, -2.000000},{6.800000, -2.000000},{6.900000, -2.000000},{7.000000, -2.000000},{7.100000, -2.000000},{7.200000, -2.000000},{7.300000, -2.000000},{7.400000, -2.000000},{7.500000, -2.000000},{7.600000, -2.000000},{7.700000, -2.000000},{7.800000, -2.000000},{7.900000, -2.000000},{8.000000, -2.000000},{8.100000, -2.000000},{8.200000, -2.000000},{8.300000, -2.000000},{8.400000, -2.000000},{8.500000, -2.000000},{8.600000, -2.000000},{8.700000, -2.000000},{8.800000, -2.000000},{8.900000, -2.000000},{9.000000, -2.000000},{9.100000, -2.000000},{9.200000, -2.000000},{9.300000, -2.000000},{9.400000, -2.000000},{9.500000, -2.000000},{9.600000, -2.000000},{9.700000, -2.000000},{9.800000, -2.000000},{9.900000, -2.000000},{10.000000, -2.000000} };
	std::vector<std::complex<long double>> referencesFull = { {1.500000, -2.000000}, { 1.600000, -2.000000 }, { 1.700000, -2.000000 }, { 1.800000, -2.000000 }, { 1.900000, -2.000000 }, { 2.000000, -2.000000 }, { 2.100000, -2.000000 }, { 2.200000, -2.000000 }, { 2.300000, -2.000000 }, { 2.400000, -2.000000 }, { 2.500000, -2.000000 }, { 2.600000, -2.000000 }, { 2.700000, -2.000000 }, { 2.800000, -2.000000 }, { 2.900000, -2.000000 }, { 3.000000, -2.000000 }, { 3.100000, -2.000000 }, { 3.200000, -2.000000 }, { 3.300000, -2.000000 }, { 3.400000, -2.000000 }, { 3.500000, -2.000000 }, { 3.600000, -2.000000 }, { 3.700000, -2.000000 }, { 3.800000, -2.000000 }, { 3.900000, -2.000000 }, { 4.000000, -2.000000 }, { 4.100000, -2.000000 }, { 4.200000, -2.000000 }, { 4.300000, -2.000000 }, { 4.400000, -2.000000 }, { 4.500000, -2.000000 }, { 4.600000, -2.000000 }, { 4.700000, -2.000000 }, { 4.800000, -2.000000 }, { 4.900000, -2.000000 }, { 5.000000, -2.000000 }, { 5.100000, -2.000000 }, { 5.200000, -2.000000 }, { 5.300000, -2.000000 }, { 5.400000, -2.000000 }, { 5.500000, -2.000000 }, { 5.600000, -2.000000 }, { 5.700000, -2.000000 }, { 5.800000, -2.000000 }, { 5.900000, -2.000000 }, { 6.000000, -2.000000 }, { 6.100000, -2.000000 }, { 6.200000, -2.000000 }, { 6.300000, -2.000000 }, { 6.400000, -2.000000 }, { 6.500000, -2.000000 }, { 6.600000, -2.000000 }, { 6.700000, -2.000000 }, { 6.800000, -2.000000 }, { 6.900000, -2.000000 }, { 7.000000, -2.000000 }, { 7.100000, -2.000000 }, { 7.200000, -2.000000 }, { 7.300000, -2.000000 }, { 7.400000, -2.000000 }, { 7.500000, -2.000000 }, { 7.600000, -2.000000 }, { 7.700000, -2.000000 }, { 7.800000, -2.000000 }, { 7.900000, -2.000000 }, { 8.000000, -2.000000 }, { 8.100000, -2.000000 }, { 8.200000, -2.000000 }, { 8.300000, -2.000000 }, { 8.400000, -2.000000 }, { 8.500000, -2.000000 }, { 8.600000, -2.000000 }, { 8.700000, -2.000000 }, { 8.800000, -2.000000 }, { 8.900000, -2.000000 }, { 9.000000, -2.000000 }, { 9.100000, -2.000000 }, { 9.200000, -2.000000 }, { 9.300000, -2.000000 }, { 9.400000, -2.000000 }, { 9.500000, -2.000000 }, { 9.600000, -2.000000 }, { 9.700000, -2.000000 }, { 9.800000, -2.000000 }, { 9.900000, -2.000000 }, { 10.000000, -2.000000 }, { 10.100000, -2.000000 }, { 10.200000, -2.000000 }, { 10.300000, -2.000000 }, { 10.400000, -2.000000 }, { 10.500000, -2.000000 }, { 10.600000, -2.000000 }, { 10.700000, -2.000000 }, { 10.800000, -2.000000 }, { 10.900000, -2.000000 }, { 11.000000, -2.000000 } };

	
	//----------------------------------------------------------------------------------------------------------------


	////////////////Frequency Perturbation//////////////////////////////////////////////////////////////////
	//std::vector<std::complex<long double>> references = { 10000000.000000, 12142857.142857, 14285714.285714, 16428571.428571, 18571428.571429, 20714285.714286, 22857142.857143, 25000000.000000, 27142857.142857, 29285714.285714, 31428571.428571, 33571428.571429, 35714285.714286, 37857142.857143, 40000000.000000};
	//std::vector<std::complex<long double>> references = { 10000000.000000, 17500000.000000, 25000000.000000, 32500000.000000, 40000000.000000 };

	//std::vector<std::complex<long double>> references = { 10000000.000000, 12727272.727273, 15454545.454545, 18181818.181818, 20909090.909091, 23636363.636364, 26363636.363636, 29090909.090909, 31818181.818182, 34545454.545455, 37272727.272727, 40000000.000000 };
	//std::vector<std::complex<long double>> references = { 10000000.000000, 11500000.000000, 13000000.000000, 14500000.000000, 16000000.000000, 17500000.000000, 19000000.000000, 20500000.000000, 22000000.000000, 23500000.000000, 25000000.000000, 26500000.000000, 28000000.000000, 29500000.000000, 31000000.000000, 32500000.000000, 34000000.000000, 35500000.000000, 37000000.000000, 38500000.000000, 40000000.000000};
	//std::vector<std::complex<long double>> references = { 15000000,15416666.6666667,15833333.3333333,16250000,16666666.6666667,17083333.3333333,17500000,17916666.6666667,18333333.3333333,18750000,19166666.6666667,19583333.3333333,20000000,20416666.6666667,20833333.3333333,21250000,21666666.6666667,22083333.3333333,22500000,22916666.6666667,23333333.3333333,23750000,24166666.6666667,24583333.3333333,25000000,25416666.6666667,25833333.3333333,26250000,26666666.6666667,27083333.3333333,27500000,27916666.6666667,28333333.3333333,28750000,29166666.6666667,29583333.3333333,30000000,30416666.6666667,30833333.3333333,31250000,31666666.6666667,32083333.3333333,32500000,32916666.6666667,33333333.3333333,33750000,34166666.6666667,34583333.3333333,35000000,35416666.6666667,35833333.3333333,36250000,36666666.6666667,37083333.3333333,37500000,37916666.6666667,38333333.3333333,38750000,39166666.6666667,39583333.3333333,40000000 };
	////-----------------------------------------------------------------------------------------------------

	long double size = referencesFull.size();

	//initial set of the references to be passed into HOPS (first, middle, last) requires odd number of references in the full vector
	std::vector<std::complex<long double>> references;
	std::vector<int> refIndex;
	//references = { referencesFull[0], referencesFull[5], referencesFull[10], referencesFull[15], referencesFull[20], referencesFull[25], referencesFull[30], referencesFull[35], referencesFull[40], referencesFull[45], referencesFull[50], referencesFull[55], referencesFull[60], referencesFull[65], referencesFull[70]};
	//refIndex = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70};
	references = referencesFull;

	////easy sweep condition
	refIndex.clear();
	references.clear();
	for (int i = 0; i < 95; i += 25) {
		refIndex.push_back(i);
		references.push_back(referencesFull[i]);
	}
	if (refIndex[refIndex.size() - 1] != 95) {
		refIndex.push_back(95);
		references.push_back(referencesFull[95]);
	}
	////-------


	//references = { referencesFull[5], referencesFull[10], referencesFull[15], referencesFull[20], referencesFull[25], referencesFull[30],referencesFull[35], referencesFull[40], referencesFull[45], referencesFull[50], referencesFull[55], referencesFull[60], referencesFull[65], referencesFull[70] };
	//references = { referencesFull[5], referencesFull[9], referencesFull[13], referencesFull[17], referencesFull[21], referencesFull[25], referencesFull[29], referencesFull[33], referencesFull[37], referencesFull[41], referencesFull[45], referencesFull[49], referencesFull[53], referencesFull[57], referencesFull[61], referencesFull[65], referencesFull[70] };

	//std::vector<std::complex<long double>> references = { referencesFull[0], referencesFull[(size - 1) / 2], referencesFull[size - 1] };
	//references = { referencesFull[0], referencesFull[1], referencesFull[2], referencesFull[4], referencesFull[5], referencesFull[6], referencesFull[8], referencesFull[9], referencesFull[10], referencesFull[12], referencesFull[13], referencesFull[14], referencesFull[16], referencesFull[17], referencesFull[18], referencesFull[20], referencesFull[21], referencesFull[22], referencesFull[24], referencesFull[25], referencesFull[26], referencesFull[28], referencesFull[29], referencesFull[30], referencesFull[32], referencesFull[33], referencesFull[34], referencesFull[36], referencesFull[37], referencesFull[38], referencesFull[40], referencesFull[41], referencesFull[42], referencesFull[44], referencesFull[45], referencesFull[46], referencesFull[48], referencesFull[49], referencesFull[50], referencesFull[52], referencesFull[53], referencesFull[54], referencesFull[56], referencesFull[57], referencesFull[58], referencesFull[60], referencesFull[61], referencesFull[62], referencesFull[64]};
	
	//references = referencesFull;
	//refIndex = { 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 70 };


	for (int i = 0; i < refIndex.size(); ++i) {
		std::cout << "i: " << i << "  references[i]: " << references[i] << "  refIndex[i]: " << refIndex[i] << std::endl;
	}

	////easy sweep condition
	//refIndex.clear();
	//references.clear();
	//for (int i = 0; i < 71; i += 1) {
	//	refIndex.push_back(i);
	//	references.push_back(referencesFull[i]);
	//}
	//if (refIndex[refIndex.size() - 1] != 70) {
	//	refIndex.push_back(70);
	//	references.push_back(70);
	//}
	////-------


	//gradient for use in spline interpolation
	std::vector <std::complex<long double>> gradient(references.size());
	std::vector <std::complex<long double>> gradient2(references.size());
	//references = { referencesFull[0], referencesFull[22], referencesFull[44], referencesFull[65] };
	
	//std::vector<std::complex<long double>> references = referencesFull;
	std::vector<std::vector<std::complex<long double>>> qoi_list(references.size());
	std::vector<std::vector<std::complex<long double>>> HOPS_splitting(references.size());

	//refIndex is the index of the ref_i file so that it can be called correctly from the reference_files folder
	//std::vector<int> refIndex = { 0, int(size - 1) / 2, int(size - 1) };
	//refIndex = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22, 24, 25, 26, 28, 29, 30, 32, 33, 34, 36, 37, 38, 40, 41, 42, 44, 45, 46, 48, 49, 50, 52, 53, 54, 56, 57, 58, 60, 61, 62, 64};




	for (int i = 0; i < refIndex.size(); ++i) {
		std::cout << "NEW i: " << i << "  references[i]: " << references[i] << "  refIndex[i]: " << refIndex[i] << std::endl;
	}

	//toggle to end the while loop
	bool toggle = true;
	
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

	std::cout << "MAterial list size: " << material_list.size() << std::endl;

	while (toggle) {
		std::cout << "_____________________________________________________________Pass in the for loop_______________________________________________\n";

		qoi_list.clear();
		HOPS_splitting.clear();
		RCS.clear();

		qoi_list.resize(references.size());
		HOPS_splitting.resize(references.size());
		RCS.resize(references.size());


		for (int i = 0; i < references.size(); ++i) {
			std::cout << "references[i]: " << references[i] << "  refIndex: " << refIndex[i] << std::endl;
		}


		//diff on the low and high side of the reference point (subdivisions mean that it isn't always the same)
		std::vector<long double> diffLo(references.size());
		std::vector<long double> diffHi(references.size());
		for (int i = 0; i < diffHi.size() - 1; ++i) {
			diffHi[i] = abs(references[i + 1].real() - references[i].real()) / 2.0;
			if (i == 0)
				diffLo[i] = diffHi[i];
			else
				diffLo[i] = abs(references[i - 1].real() - references[i].real()) / 2.0;

		}

		//this sets the diff for the last element diffLo[end] = diffLo[end-1] diffHi[end] shouldn't be needed provided the final reference is higher than any monte carlo point
		diffLo[diffLo.size() - 1] = diffHi[diffLo.size() - 2];
		diffHi[diffHi.size() - 1] = diffHi[diffLo.size() - 2];
		//long double diff = (references[1].real() - references[0].real()) / 2.0;

		int mat_counter = 0;

		for (int j = 0; j < material_list.size(); ++j) {

			//hops splitting is N segments CENTERED around N reference points
			long double min_dist = 1.0e15;
			int min_index = -1;

			for (int i = 0; i < references.size(); ++i) {
				//check if random var falls within upper lim
				long double dist = norm(references[i] - material_list[j]);
				if (dist < min_dist)
				{
					min_index = i;
					min_dist = dist;
				}
			}
			assert(min_index >= 0 && "Suitable reference point not found!");
			HOPS_splitting[min_index].push_back(material_list[j]);
			++mat_counter;

		}

		////hops splitting is N-1 segments BETWEEN N reference points
		//int rIdx = 0;
		////checking for monte carlo security
		//for (int j = 0; j < material_list.size(); ++j) {

		//	if (material_list[j].real() <= references[rIdx + 1].real()) {

		//		if (mat_counter == 993) {
		//			std::cout << material_list[j] << std::endl;
		//		}

		//		HOPS_splitting[rIdx].push_back(material_list[j]);
		//		mat_counter++;
		//	}
		//	else {
		//		rIdx++;
		//		j--;
		//	}
		//}

		//old method for adding to HOPS_splitting
		//for (int j = 0; j < material_list.size(); ++j) {
		//	//if (material_list[j].real() < 15.5e6) { continue; }
		//	for (int i = 0; i < references.size(); ++i) {
		//		//check if random var falls within upper lim
		//		if ((material_list[j].real() - references[i].real()) >= 0.0) {
		//			if ((material_list[j].real() - references[i].real()) <= diffHi[i]) {
		//				HOPS_splitting[i].push_back(material_list[j]);
		//				mat_counter++;
		//				break;
		//			}
		//		}
		//		else {
		//			//check if random var falls within lower lim
		//			if ((material_list[j].real() - references[i].real()) >= -diffLo[i]) {
		//				HOPS_splitting[i].push_back(material_list[j]);
		//				mat_counter++;

		//				//if (material_list[j].real() > max[i]) {
		//				//	max[i] = material_list[j].real();
		//				//	maxIndex[i] = HOPS_splitting[i].size() - 1;
		//				//}
		//				//if (material_list[j].real() < min[i]) {
		//				//	min[i] = material_list[j].real();
		//				//	min[i] = HOPS_splitting[i].size() - 1;
		//				//}

		//				break;
		//			}
		//		}

		//	}

		//}

		//these are the test values to be compared generated within HOPS.  HOPS_splitting[size-2] = min, HOPS_splitting[size-1] = max
		// ie first push_back is min value, second push_back is max value, within a given segment

		for (int i = 0; i < references.size(); ++i) {

			//HOPS_splitting[i].push_back(references[i - 1]);
			//HOPS_splitting[i].push_back(references[i + 1]);


			if (i == 0) {
				//new, not sure, used to be zero
				HOPS_splitting[i].push_back(references[i]);
				HOPS_splitting[i].push_back((references[i] + references[i + 1]) / 2.0);
			}
			else if (i == references.size()-1) {

				HOPS_splitting[i].push_back((references[i] + references[i - 1]) / 2.0);
				//new, not sure
				HOPS_splitting[i].push_back(references[i]);
			}
			else {
				HOPS_splitting[i].push_back((references[i] + references[i - 1]) / 2.0);
				HOPS_splitting[i].push_back((references[i] + references[i + 1]) / 2.0);

			}

			std::cout << "i: " << i << "  hops_splitting[i][end-1](max): " << HOPS_splitting[i][HOPS_splitting[i].size() - 1] << "  hops_splitting[i][end-2](min): " << HOPS_splitting[i][HOPS_splitting[i].size() - 2] << std::endl;

			//original method
			//if (i == 0)
			//	HOPS_splitting[i].push_back(references[i]);
			//else
			//	HOPS_splitting[i].push_back(references[i - 1]);

			//if (i == references.size()-1)
			//	HOPS_splitting[i].push_back(references[i]);
			//else
			//	HOPS_splitting[i].push_back(references[i + 1]);


		}


		std::cout << "mat_counter: " << mat_counter << std::endl;

		//for (int i = 0; i < references.size(); i++) {
		//	std::cout << "new references: " << references[i] << std::endl;
		//}

		std::vector<std::complex<long double>> referenceVals(references.size());

		//iterating through references now and get the qois
		for (int i = 0; i < references.size(); ++i) {
			if (HOPS_splitting[i].size() == 0) continue;
			in_file_name = file_name;

			//////////////////////////////Frequency Perturbation///////////////////////////////////////////////////////////////////////////////////////////////
			//std::string command = "..\\file_input_converter ..\\reference_files_freq_low_HOPS\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//std::string command = "..\\file_input_converter ..\\reference_files_freq\\ref_3.in ..\\exampleFiles\\" + in_file_name + "\\";
			//---------------------------------------------------------------------------------------------------------------------------------------------------

			/////////////////////////////Material Perturbation///////////////////////////////////////////////////////////////////////////////
			//std::string command = "..\\file_input_converter ..\\reference_files_mat\\ref_0.in ..\\exampleFiles\\" + in_file_name + "\\";
			std::string command = "..\\file_input_converter ..\\reference_files_mat_newRange4\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//--------------------------------------------------------------------------------------------------------------------------------------------------

			std::system(command.c_str());
			std::cout << "Testing for reference number: " << refIndex[i] << std::endl;
			referenceVals[i] = HOPS::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i], gradient[i]);

		}

		//end condition based on convergence
		long double size = 0.0;
		for (int i = 0; i < qoi_list.size(); i++) {
			for (int j = 0; j < qoi_list[i].size() - 2; j++) {
				sum = sum + qoi_list[i][j];
				RCS[i].push_back(sqrt(qoi_list[i][j].real() * qoi_list[i][j].real() + qoi_list[i][j].imag() * qoi_list[i][j].imag()) / (4 * 3.14159));
				sumRCS = sumRCS + RCS[i][j];
				size = size + 1.0;
			}
		}

		mean = sum / size;
		meanRCS = sumRCS / size;


		for (int i = 0; i < qoi_list.size(); i++) {
			for (int j = 0; j < qoi_list[i].size() - 2; ++j) {
				stdDev = stdDev + (qoi_list[i][j] - mean) * (qoi_list[i][j] - mean);
				stdDevRCS = stdDevRCS + (RCS[i][j] - meanRCS) * (RCS[i][j] - meanRCS);

			}
		}

		stdDev = sqrt(stdDev / size);
		stdDevRCS = sqrt(stdDevRCS / size);


		stdDevList.push_back(stdDev);
		stdDevRCSList.push_back(stdDevRCS);

		long double devChange = abs(stdDevList[stdDevList.size() - 1] - stdDevList[stdDevList.size() - 2]) / abs(stdDevList[stdDevList.size() - 1] + stdDevList[stdDevList.size() - 2]) / 2.0;
		long double devChangeRCS = abs(stdDevRCSList[stdDevRCSList.size() - 1] - stdDevRCSList[stdDevRCSList.size() - 2]) / abs(stdDevRCSList[stdDevRCSList.size() - 1] + stdDevRCSList[stdDevRCSList.size() - 2]) / 2.0;

		std::cout << "Mean: " << stdDevRCS << " std dev: " << meanRCS << std::endl;

		if (devChangeRCS < stdLim) {
			std::cout << "Convergence condition met" << std::endl;
			break;
		}


		//adding references point based on error

		long double loError = 0.0;
		long double hiError = 0.0;

		std::vector<long double> loErrorVec;
		std::vector<long double> hiErrorVec;

		//subdivide highest error
		int hiErrorIndex = 0;
		long double hiErrorCheck = 0.0;

		toggle = false;

		//val is the endpoint of HOPS_splitting
		int val = HOPS_splitting.size() - 1;
		std::vector<std::complex<long double>> referencesNew = references;
		std::vector<int> refIndexNew = refIndex;
		//int index = 0;
		// 


		
	////calculate variogram
	//long double nBins = 10;

	//long double maxDist = abs(references[references.size() - 1].real() - references[0].real());
	////long double binTol = maxDist / nBins;
	//std::vector<long double> edges;
	//for (int i = 0; i <= nBins; i++) {
	//	edges.push_back(i * maxDist / nBins);
	//}

	////linspace for the x values in the variogramFit
	//int fitSize = 100;
	//std::vector<long double> d;
	//for (int i = 0; i <= fitSize; ++i) {
	//	d.push_back(i * maxDist / fitSize);
	//}


	//std::vector<std::vector<long double>> dMat = buildDMatrix(references, referenceVals);
	//std::vector<long double> variogram = buildVariogram(dMat, references, edges);
	//trimVariogram(variogram, edges);
	//std::vector<long double> variogramFit = fitVariogram(variogram, edges, d);


	//	for (int i = 0; i < referencesFull.size(); ++i) {
	//		std::vector<std::complex<long double>> referencesTemp = references;
	//		auto it = std::upper_bound(referencesTemp.cbegin(), referencesTemp.cend(), referencesFull[i]);
	//		referencesTemp.insert(it, referencesFull[i]);


	//	}


		//iteration limit//
		iterations++;
		if (iterations >= itLim) {
			std::cout << "iteration condition met\n";
			break;
		}
		///-------------------------------------------------------------

		std::cout << "made it here\n";

		/////////////// Cubic Spline//////////////////////////////////////////

		alglib::real_1d_array x;
		alglib::real_1d_array yReal;
		alglib::real_1d_array yImag;
		x.setlength(references.size());
		yReal.setlength(references.size());
		yImag.setlength(references.size());


		for (int i = 0; i < references.size(); ++i) {
			x[i] = references[i].real();
			yReal[i] = referenceVals[i].real();
			yImag[i] = referenceVals[i].imag();

		}

		alglib::spline1dinterpolant sReal;
		alglib::spline1dinterpolant sImag;


		alglib::ae_int_t bound = 1;

		for (int i = 0; i < references.size(); ++i) {
			std::cout << references[i].real() << " " << yReal[i] << " \n";
		}

		long double re;
		long double ima;
		for (int i = 0; i < qoi_list.size() - 1; i++) {
			std::cout << "made it here2\n";
			alglib::spline1dbuildcubic(x, yReal, references.size(), bound, gradient[i].real(), bound, gradient[i+1].real(), sReal);
			alglib::spline1dbuildcubic(x, yImag, references.size(), bound, gradient[i].imag(), bound, gradient[i+1].imag(), sImag);

			//need to remove the max/min test values for each of the batches
			for (int j = 0; j < qoi_list[i].size() - 2; ++j) {
				
				re = alglib::spline1dcalc(sReal, HOPS_splitting[i][j].real());
				//std::cout << "made it here3\n";
				ima = alglib::spline1dcalc(sImag, HOPS_splitting[i][j].real());
				//std::cout << "made it here4\n";

				//qoi_list[i][j] = { re, ima };
				//std::cout << "made it here5\n";

				//std::cout << "i: " << i << " j: " << j << std::endl;
			}
		}

		//for (int i = 0; i < references.size(); i++) {
		//	if (i == references.size() - 1) {
		//		gradient2[i] = (gradient[i] - gradient[i - 1]) / references[i].real();
		//	}
		//	else {
		//		gradient2[i] = (gradient[i + 1] - gradient[i]) / references[i].real();
		//	}
		//}
		



		//////---------------------------------------------------------------


		//AR sweep output
		std::string outF = "../ioFiles/output/sweep/firstOrder4/qoi_HOPS_materialAR" + std::to_string(iterations + 2) + ".txt";
		//std::string outF = "../ioFiles/output/sweep/splineGrad/qoi_HOPS_materialAR15.txt";
		std::ofstream qoi_dist_out(outF);
		int index = 0;

		//qoi_list.size() - 1 for spline
		for (int i = 0; i < qoi_list.size(); i++) {
			//need to remove the max/min test values for each of the batches


			for (int j = 0; j < qoi_list[i].size() - 2; ++j) {
				if (index > HOPS_splitting[index].size()) index++;

				//std::cout << "i: " << i << " j: " << j << std::endl;
				qoi_dist_out << std::setprecision(12) << HOPS_splitting[i][j].real() << " " << HOPS_splitting[i][j].imag() << " " << qoi_list[i][j].real() << " " << qoi_list[i][j].imag() << std::endl;
			}
		}

		qoi_dist_out.close();
		//just quits (no iteration)
		break;
	}


	std::cout << "out of the while loop=============================\n";


	//output one time
	//output results to file
	//std::ofstream qoi_dist_out("../ioFiles/output/sweep/weight/qoi_HOPS_materialAR0.txt");
	//int index = 0;
	//for (int i = 0; i < qoi_list.size(); i++) {
	//	//need to remove the max/min test values for each of the batches
	//	for (int j = 0; j < qoi_list[i].size()-2; ++j) {
	//		if (index > HOPS_splitting[index].size()) index++;
	//		//std::cout << "i: " << i << " j: " << j << std::endl;
	//		qoi_dist_out << HOPS_splitting[i][j].real() << " " << HOPS_splitting[i][j].imag() << " " << qoi_list[i][j].real() << " " << qoi_list[i][j].imag() << std::endl;
	//	}
	//}
	//qoi_dist_out.close();
}

std::vector<std::complex<long double>> spline(std::vector<std::complex<long double>> gradient2, std::vector<std::complex<long double>> references, std::vector<std::complex<long double>> referenceVals) {

}

long double HOPS::normalFunction(long double x, long double mean, long double std) {
	long double arg = (x - mean) / std;
	long double exponential = exp(-0.5 * arg * arg);
	long double integral = 1.0 / (std * sqrt(2.0 * 3.14159)) * exponential;
	return integral;
}

long double HOPS::trap_integral(long double intervalBeg, long double intervalEnd, long double mean, long double std) {
	long double count = 250.0;
	long double step = (intervalEnd - intervalBeg) / count;
	long double integral = 0.5 * (HOPS::normalFunction(intervalBeg, mean, std) + HOPS::normalFunction(intervalEnd, mean, std));
	for (int i = 0; i < (int)count; ++i) {
		integral += HOPS::normalFunction(intervalBeg + step * i, mean, std);
	}
	integral *= step;
	return integral;

}
/*
void HOPS::multi_HOPS_multi_epsr(std::string & file_name) //mutli HOPS for variation of multiple material values
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
		HOPS::sensitivity_to_multi_epsr(file_name, HOPS_splitting[i], qoi_list[i], references[i]);
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

void HOPS::sensitivity_to_multi_epsr(std::string & file_name, std::vector<std::vector<std::complex<long double>>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::vector<std::complex<long double>>& reference_values)
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