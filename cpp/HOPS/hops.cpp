#include "hops.h"
std::complex<double> get_QoI2(Domain& dom_forward, Domain& dom_adjoint) {
	std::complex<double> qoi = 0.0;
	for (int i = 0; i < dom_forward.cAlpha.size(); ++i) {
		qoi += dom_forward.cAlpha[i] * std::conj(dom_adjoint.scatter1.cGr[1][i]);
	}
	return qoi;
}

inline std::complex<double> operator/(std::complex<double> c1, double d1) {
	std::complex<double> output = { c1.real() / d1, c1.imag() / d1 };
	return output;
}

double norm(std::vector < std::complex<double>> v1) {
	std::complex<double> norm = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		norm += v1[i] * std::conj(v1[i]);
	}
	double normd = abs(norm);
	return sqrt(normd);

}

std::vector<std::complex<double>> vectorInsert(int position, std::complex<double> value, std::vector<std::complex<double>> vector) {
	std::vector<std::complex<double>> newVector(vector.size() + 1);

	for (int i = 0; i < newVector.size(); ++i) {
		while (i < position) {
			newVector[i] = vector[i];
		}
		if (i == position) {
			newVector[i] = value;
		}
		while (i > position) {
			newVector[i] = vector[i - 1];
		}
	}

	return newVector;

}

inline std::vector<std::complex<double>> operator-(std::vector < std::complex<double>>& v1, std::vector < std::complex<double>>& v2) {
	std::vector<std::complex<double>> vout;
	for (int i = 0; i < v1.size(); ++i) {
		vout.push_back(v1[i] - v2[i]);
	}
	return vout;
}
std::complex<double> HOPS::sensitivity_to_epsr(std::string & file_name, std::vector<std::complex<double>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::complex<double> referenceFreq)
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

	std::cout << "reference qoi: " << reference_qoi << " scattered field: " << std::endl;
	
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
	std::cout << "current dQoI: " << dQoI << std::endl;
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
	double kReal = 2 * 3.14159 / (2.99E8 / referenceFreq.real());
	std::complex<double> referenceK = (kReal, 0.0);
	//now have the base QoI modifier, need to scale by the change in perturbed material parameter
	for (int i = 0; i < epsr_list.size(); ++i) {
		/////////////////////////////////////////Frequency Perturbation/////////////////////////////////////////////
		//double k0_list = epsr_list[i].real() * 2.0 * 3.14159 / 3e8;
		//eps_diff = -1.0*k0_list + k0;
		//eps_diff = (k0_list / k0 - 1.0) * k0;
		//---------------------------------------------------------------------------------------------------------

		///////////////////////////////////////Material Perturbation//////////////////////////////////////////////
		eps_diff = epsr_list[i] - referenceFreq;
		///-------------------------------------------------------------------------------------------------------

		updated_qoi.push_back(reference_qoi + eps_diff * dQoI);// +eps_diff * eps_diff * dQoI2 / 2.0);
		//std::cout << reference_qoi + eps_diff * dQoI << std::endl;
		//std::cout << "qoi: " << reference_qoi + eps_diff * dQoI << std::endl;
		
	}
	return reference_qoi;
}

std::complex<double> HOPS::get_QoI(Domain & dom_forward, Domain & dom_adjoint)
{
	std::complex<double> qoi = 0.0;
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
	//mesh inputs
	std::string mesh_name = file_name;
	//get the forward solution
	Domain dom_forward(mesh_name);
	run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, false, 0.0);
	//get RHS
	Domain dom_adjoint_rhs(mesh_name);
	run::_RHS_compute(dom_adjoint_rhs, mesh_name, plot, true, higher_order, check_results);
	//get QoI
	std::complex<double> qoi = get_QoI2(dom_forward, dom_adjoint_rhs);
}

void HOPS::multi_HOPS_epsr(std::string & file_name)
{


	//load in the list of random materials being tested
	std::vector<std::complex<double>> material_list;
	std::ifstream materials_in("../ioFiles/input/materials_list.txt");
	std::string line;
	//std::cout << "Current rounding material values to match the shitty output from MATLAB!" << std::endl;
	while (std::getline(materials_in, line)) {
		auto data = functions::split(line, '\t');
		//material_list.push_back(std::complex<double>(std::stod(data[0]), roundf(std::stod(data[1])*100.0)/100.0));
		material_list.push_back(std::complex<double>(std::stod(data[0]), std::stod(data[1])));
	}
	materials_in.close();

	//acceptable error
	double delta = 0.0008;

	//convergence condition
	double stdLim = 0.1;

	//iteration limit
	int itLim = 7;

	std::vector<double> stdDevRCSList;
	std::vector<std::complex<double>> stdDevList;


	///////////////////////Material Perturbation/////////////////////////////////////////////////////////////////////
	//std::vector<std::complex<double>> references = { {4.0, -2.0} };
	std::vector<std::complex<double>> referencesFull = { {1.000000, -2.000000},{1.666667, -2.000000},{2.333333, -2.000000},{3.000000, -2.000000},{3.666667, -2.000000},{4.333333, -2.000000},{5.000000, -2.000000},{5.666667, -2.000000},{6.333333, -2.000000},{7.000000, -2.000000} };
	//----------------------------------------------------------------------------------------------------------------

	////////////////Frequency Perturbation//////////////////////////////////////////////////////////////////
	//std::vector<std::complex<double>> references = { 10000000.000000, 12142857.142857, 14285714.285714, 16428571.428571, 18571428.571429, 20714285.714286, 22857142.857143, 25000000.000000, 27142857.142857, 29285714.285714, 31428571.428571, 33571428.571429, 35714285.714286, 37857142.857143, 40000000.000000};
	//std::vector<std::complex<double>> references = { 10000000.000000, 17500000.000000, 25000000.000000, 32500000.000000, 40000000.000000 };

	//std::vector<std::complex<double>> references = { 10000000.000000, 12727272.727273, 15454545.454545, 18181818.181818, 20909090.909091, 23636363.636364, 26363636.363636, 29090909.090909, 31818181.818182, 34545454.545455, 37272727.272727, 40000000.000000 };
	//std::vector<std::complex<double>> references = { 10000000.000000, 11500000.000000, 13000000.000000, 14500000.000000, 16000000.000000, 17500000.000000, 19000000.000000, 20500000.000000, 22000000.000000, 23500000.000000, 25000000.000000, 26500000.000000, 28000000.000000, 29500000.000000, 31000000.000000, 32500000.000000, 34000000.000000, 35500000.000000, 37000000.000000, 38500000.000000, 40000000.000000};
	//std::vector<std::complex<double>> references = { 15000000,15416666.6666667,15833333.3333333,16250000,16666666.6666667,17083333.3333333,17500000,17916666.6666667,18333333.3333333,18750000,19166666.6666667,19583333.3333333,20000000,20416666.6666667,20833333.3333333,21250000,21666666.6666667,22083333.3333333,22500000,22916666.6666667,23333333.3333333,23750000,24166666.6666667,24583333.3333333,25000000,25416666.6666667,25833333.3333333,26250000,26666666.6666667,27083333.3333333,27500000,27916666.6666667,28333333.3333333,28750000,29166666.6666667,29583333.3333333,30000000,30416666.6666667,30833333.3333333,31250000,31666666.6666667,32083333.3333333,32500000,32916666.6666667,33333333.3333333,33750000,34166666.6666667,34583333.3333333,35000000,35416666.6666667,35833333.3333333,36250000,36666666.6666667,37083333.3333333,37500000,37916666.6666667,38333333.3333333,38750000,39166666.6666667,39583333.3333333,40000000 };
	////-----------------------------------------------------------------------------------------------------

	double size = referencesFull.size();

	//initial set of the references to be passed into HOPS (first, middle, last) requires odd number of references in the full vector
	std::vector<std::complex<double>> references = { referencesFull[0], referencesFull[(size - 1) / 2], referencesFull[size - 1] };
	//std::vector<std::complex<double>> references = referencesFull;
	std::vector<std::vector<std::complex<double>>> qoi_list(references.size());
	std::vector<std::vector<std::complex<double>>> HOPS_splitting(references.size());

	//refIndex is the index of the ref_i file so that it can be called correctly from the reference_files folder
	std::vector<int> refIndex = { 0, int(size - 1) / 2, int(size - 1) };
	//std::vector<int> refIndex = {  };


	//toggle to end the while loop
	bool toggle = true;
	
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
		std::vector<double> diffLo(references.size());
		std::vector<double> diffHi(references.size());
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
		//double diff = (references[1].real() - references[0].real()) / 2.0;

		int mat_counter = 0;


		for (int j = 0; j < material_list.size(); ++j) {
			//if (material_list[j].real() < 15.5e6) { continue; }
			for (int i = 0; i < references.size(); ++i) {
				//check if random var falls within upper lim
				if ((material_list[j].real() - references[i].real()) >= 0.0) {
					if ((material_list[j].real() - references[i].real()) <= diffHi[i]) {
						HOPS_splitting[i].push_back(material_list[j]);
						mat_counter++;
						break;
					}
				}
				else {
					//check if random var falls within lower lim
					if ((material_list[j].real() - references[i].real()) >= -diffLo[i]) {
						HOPS_splitting[i].push_back(material_list[j]);
						mat_counter++;

						//if (material_list[j].real() > max[i]) {
						//	max[i] = material_list[j].real();
						//	maxIndex[i] = HOPS_splitting[i].size() - 1;
						//}
						//if (material_list[j].real() < min[i]) {
						//	min[i] = material_list[j].real();
						//	min[i] = HOPS_splitting[i].size() - 1;
						//}

						break;
					}
				}

			}

		}

		//these are the test values to be compared generated within HOPS.  HOPS_splitting[size-2] = min, HOPS_splitting[size-1] = max
		// ie first push_back is min value, second push_back is max value, within a given segment

		for (int i = 0; i < references.size(); ++i) {

			//HOPS_splitting[i].push_back(references[i - 1]);
			//HOPS_splitting[i].push_back(references[i + 1]);


			if (i == 0) {
				HOPS_splitting[i].push_back(0.0);
				HOPS_splitting[i].push_back((references[i] + references[i + 1]) / 2.0);
			}
			else if (i == references.size() - 1) {
				HOPS_splitting[i].push_back((references[i] + references[i - 1]) / 2.0);
				HOPS_splitting[i].push_back(0.0);
			}
			else {
				HOPS_splitting[i].push_back((references[i] + references[i - 1]) / 2.0);
				HOPS_splitting[i].push_back((references[i] + references[i + 1]) / 2.0);

			}
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

		for (int i = 0; i < references.size(); i++) {
			std::cout << "i: " << i << " (min/max)  HOPS_Splitting[i][size-2]: " << HOPS_splitting[i][HOPS_splitting[i].size() - 2] << "  HOPS_Splitting[i][size-1]: " << HOPS_splitting[i][HOPS_splitting[i].size() - 1] << std::endl;
		}

		for (int i = 0; i < references.size(); i++) {
			for (int j = HOPS_splitting[i].size() - 2; j < HOPS_splitting[i].size(); j++) {
				std::cout << "New check: HOPS_Splitting[i][j]" << i << " " << j << " " << HOPS_splitting[i][j] << std::endl;
			}
		}


		std::cout << "mat_counter: " << mat_counter << std::endl;

		//for (int i = 0; i < references.size(); i++) {
		//	std::cout << "new references: " << references[i] << std::endl;
		//}

		std::vector<std::complex<double>> referenceVals(references.size());

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
			std::string command = "..\\file_input_converter ..\\reference_files_mat\\ref_" + std::to_string(refIndex[i]) + ".in ..\\exampleFiles\\" + in_file_name + "\\";
			//--------------------------------------------------------------------------------------------------------------------------------------------------

			std::system(command.c_str());
			std::cout << "Testing for reference number: " << refIndex[i] << std::endl;
			referenceVals[i] = HOPS::sensitivity_to_epsr(in_file_name, HOPS_splitting[i], qoi_list[i], references[i]);

		}

		//end condition based on convergence
		double size = 0.0;
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

		double devChange = abs(stdDevList[stdDevList.size() - 1] - stdDevList[stdDevList.size() - 2]) / abs(stdDevList[stdDevList.size() - 1] + stdDevList[stdDevList.size() - 2]) / 2.0;
		double devChangeRCS = abs(stdDevRCSList[stdDevRCSList.size() - 1] - stdDevRCSList[stdDevRCSList.size() - 2]) / abs(stdDevRCSList[stdDevRCSList.size() - 1] + stdDevRCSList[stdDevRCSList.size() - 2]) / 2.0;

		std::cout << "Mean: " << stdDevRCS << " std dev: " << meanRCS << std::endl;

		if (devChangeRCS < stdLim) {
			std::cout << "Convergence condition met" << std::endl;
			break;
		}


		//adding references point based on error

		double loError = 0.0;
		double hiError = 0.0;

		//subdivide highest error
		int hiErrorIndex = 0;
		double hiErrorCheck = 0.0;

		toggle = false;

		//val is the endpoint of HOPS_splitting
		int val = HOPS_splitting.size() - 1;
		std::vector<std::complex<double>> referencesNew = references;
		std::vector<int> refIndexNew = refIndex;
		int index = 0;
		std::cout << "++++++++++++++++++++++val: " << val << std::endl;
		//std::cout << "references.size(): " << references.size() << std::endl;


		for (int i = 0; i < references.size(); i++) {
			for (int j = HOPS_splitting[i].size() - 2; j < HOPS_splitting[i].size(); j++) {
				std::cout << "New check___2: HOPS_Splitting[i][j]" << i << " " << j << " " << HOPS_splitting[i][j] << std::endl;
			}
		}


		for (int i = 0; i < references.size() - 1; ++i) {


			//check for min or max end value of HOPS_splitting
			int jMin = HOPS_splitting[i+1].size() - 2;
			int jMax = HOPS_splitting[i].size() - 1;

			//stupid add iterators, fuck you
			auto it = references.begin();
			auto it2 = refIndex.begin();
			//if (i == references.size() - 1) { break; }

			/*if (i == 0)
				loError = 0.0;
			else*/
			//loError = std::abs(referenceVals[i-1] - qoi_list[i][references.size() - 2]) / std::abs(referenceVals[i-1]);


			// 
			//shouldn't need this with only doing hiError at the moment
			//if (i == index - 1) {
			//	hiError = 0.0;
			//}
			//else

			//qoi_list[i][val] should be max, qoi_list[i][val-1] is min

			//[i][val] is hi side on current point, [i+1][val-1] is lo side on the next point


			//  max - min(i+1) / 
			hiError = std::abs(qoi_list[i][jMax] - qoi_list[i + 1][jMin]) / (std::abs(qoi_list[i][jMax] + qoi_list[i + 1][jMin]) / 2.0);

			std::cout << "hiError: " << hiError << "  i: " << i << " j: " << jMin << "  refIndex[i]: " << refIndex[i] << " splitting[i][j+1] " << HOPS_splitting[i][jMax] << "  splitting[i+1][j] " << HOPS_splitting[i + 1][jMin] << std::endl;



			//this is done to subdivide highest error
			if (hiError > hiErrorCheck) {


				//check to see if the value is already contained in the vector
				int refFullIndex = (refIndex[i + 1] - refIndex[i]) / 2 + refIndex[i];
				if (std::find(refIndex.begin(), refIndex.end(), refFullIndex) != refIndex.end()) {
					std::cout << "found a duplicate\n";
					break;
				}

				hiErrorIndex = i;
				hiErrorCheck = hiError;
			}
			//old method

			////////////////////////////////this method subdivides everywhere///////////////////////////////////

			//hiError = std::abs(referenceVals[i + 1] - qoi_list[i][val]) / std::abs(referenceVals[i + 1]);
			//			
			//if (hiError > delta){
			//	toggle = true;
			//	
			//	int refFullIndex = (refIndex[i + index + 1] - refIndex[i + index]) / 2 + refIndex[i + index];
			//	references.insert(it + i + index + 1, referencesFull[refFullIndex]);
			//	
			//    refIndex.insert(it2 + i + index + 1, refFullIndex);
			//	//account for the extra element in references
			//	index++;
			//	std::cout << "i: " << i << "    hi error: " << hiError << "     references[i+1]: " << references[i+1] << "    references[i]: " << references[i] << std::endl;
			//}
			///------------------------------------------------------------------------------------------------

			//don't need currently
			//if (loError < delta) {
			//	toggle = true;
			//	auto it = references.begin() + i;
			//	references.insert(it, (references[i] - references[i - 1]) / 2.0);
			//	refIndex.insert(it, (refIndex[i] - refIndex[i - 1]) / 2);
			//}
	}

		for (int i = 0; i < references.size(); i++) {
			for (int j = HOPS_splitting[i].size() - 2; j < HOPS_splitting[i].size(); j++) {
				std::cout << "New check___3____________________________: HOPS_Splitting[i][j]" << i << " " << j << " " << HOPS_splitting[i][j] << std::endl;
			}
		}

		////////////done to subdivide lowest error////////////////////
		auto it = references.begin();
		auto it2 = refIndex.begin();
		int refFullIndex = (refIndex[hiErrorIndex + 1] - refIndex[hiErrorIndex]) / 2 + refIndex[hiErrorIndex];
		references.insert(it + hiErrorIndex + 1, referencesFull[refFullIndex]);
		refIndex.insert(it2 + hiErrorIndex + 1, refFullIndex);
		toggle = true;

		iterations++;
		if (iterations >= itLim) {
			std::cout << "iteration condition met\n";
			break;
		}
		///-------------------------------------------------------------

	}

	std::cout << "out of the while loop=============================\n";

	//output results to file
	std::ofstream qoi_dist_out("../ioFiles/output/sweep/qoi_HOPS_materialAR6.txt");
	int index = 0;
	for (int i = 0; i < qoi_list.size(); i++) {
		//need to remove the max/min test values for each of the batches
		for (int j = 0; j < qoi_list[i].size()-2; ++j) {
			if (index > HOPS_splitting[index].size()) index++;
			//std::cout << "i: " << i << " j: " << j << std::endl;
			qoi_dist_out << HOPS_splitting[i][j].real() << " " << HOPS_splitting[i][j].imag() << " " << qoi_list[i][j].real() << " " << qoi_list[i][j].imag() << std::endl;
		}
	}
	qoi_dist_out.close();
}

void HOPS::multi_HOPS_multi_epsr(std::string & file_name) //mutli HOPS for variation of multiple material values
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
		/*std::string command = "C:\\Users\\jjh_3\\git\\FEMCUDA\\cpp\\file_input_converter C:\\Users\\jjh_3\\git\\FEMCUDA\\cpp\\reference_files\\ref_" + std::to_string(i) + ".in C:\\Users\\jjh_3\\git\\FEMCUDA\\cpp\\examplesFiles\\" + file_name + "\\";
		std::system(command.c_str());*/
		/*std::string command = "..\\file_input_converter ..\\reference_files_2im\\ref_" + std::to_string(i) + ".in ..\\exampleFiles\\" + file_name + "\\";*/
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

void HOPS::sensitivity_to_multi_epsr(std::string & file_name, std::vector<std::vector<std::complex<double>>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::vector<std::complex<double>>& reference_values)
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