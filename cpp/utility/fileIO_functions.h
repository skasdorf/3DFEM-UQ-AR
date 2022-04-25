//#include <vector>
//#include <complex>
//#include <fstream>
//#include "functions.h"
//namespace filefunctions {
//	void read_previous_solve(std::string mesh_name, std::vector<std::complex<double>>& cAlphaForLower, std::vector<std::complex<double>>& cAlphaAdjoint, std::vector<std::complex<double>>& cGrhigher){
//		std::ifstream file;
//		std::string line;
//		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha.txt"); //lower order cAlpha computation
//																					  //getline(file, line);
//		std::vector<std::string> line_input;
//		while (getline(file, line)) {
//			line_input = functions::split(line, ' ');
//			cAlphaForLower.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
//
//		}
//		file.close();
//		//
//		file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt"); //higher order for cAlphaAdjoint
//		
//		while (getline(file, line)) {
//			line = line.substr(1, line.size() - 2);
//			line_input = functions::split(line, ',');
//			cAlphaAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
//
//		}
//		file.close();
//		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt"); //higher order for cGr from forwars solve
//	
//		while (getline(file, line)) {
//			line = line.substr(1, line.size() - 2);
//			line_input = functions::split(line, ',');
//			cGrhigher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1]))); //only works for one freq atm
//
//		}
//		file.close();
//
//	}
//}
