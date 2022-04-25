//input file converter: takes FEMIN.DAT and creates a input files with the information we require
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "fileFacet.h"
#include "materials_input.h"
std::vector<std::string> split(const std::string &s, char delim) {
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> tokens;
	while (std::getline(ss, item, delim)) {
		if (item == "") continue;
		tokens.push_back(item);
	}
	return tokens;
}
int main(int argc, char** argv) {
	if (argc <= 1) std::cout << "Using defaults!" << std::endl;
	std::cout << "Num args: " << argc << std::endl;
	for (int i = 0; i < argc; ++i) {
		std::cout << argv[i] << std::endl;
	}
	//agr0 = file location of executable
	//arg1 = input file name
	//arg2 = output FOLDER location
	std::string line;
	std::string line2;
	std::string line3;

	std::ifstream file_in;
	std::ofstream file_out;

	std::ifstream file_in_10; //for order_u, etc
	std::ifstream file_in_6; //for element type: pml, etc.
	std::ifstream file_in_8; // for bc
	//file_in.open("../FEM_Mesh.in");
	//file_in_10.open("../FEM_Mesh.in");
	//file_in_6.open("../FEM_Mesh.in");

	/*file_in.open("../FEM_Mesh.in");
	file_in_10.open("../FEM_Mesh.in");
	file_in_6.open("../FEM_Mesh.in");
	file_in_8.open("../FEM_Mesh.in");*/

	//std::string input_file_name = "../FEM_Mesh_cube.in";
	//std::string input_file_name = "../FEM_Mesh_2nd_1element_split.in";

	std::string input_file_name = "../sphere_d.in";
	std::string location = "../";
	if (argc > 1) {
		input_file_name = argv[1];
		location = argv[2];
	}
	file_in.open(input_file_name);
	file_in_10.open(input_file_name);
	file_in_6.open(input_file_name);
	file_in_8.open(input_file_name);
	//int j = 0;
	while (getline(file_in_6, line) && line.find("6. 3D ELEMENT") == std::string::npos) {
		//j++;
		//skip line
	}
	while (getline(file_in_10, line) && line.find("10. 3D NUMERICS") == std::string::npos) {
		//skip line
	}
	while (getline(file_in_8, line) && line.find("8. 3D BOUNDARY") == std::string::npos) {
		//skips lines
	}
	for (int i = 0; i < 23; ++i) {
		getline(file_in_8, line);
	}
	
	for (int i = 0; i < 45; ++i) {
		getline(file_in, line);
	}
	std::replace(line.begin(), line.end(), '\t', ' ');
	std::vector<std::string> basic_data = split(line, ' ');
	int number_of_elements = std::stoi(basic_data[2]);
	for (int i = 0; i < 11; ++i) {
		getline(file_in, line);
	}
	for (int i = 0; i < 35; ++i) {
		getline(file_in_6, line);
	}
	for (int i = 0; i < 17; ++i) {
		getline(file_in_10, line);
	}
	//copy the coordinate lines and send to a file name geometry.txt
	file_out.open(location + "geometry.txt");
	file_out << "!Geometry.txt" << std::endl;
	while (getline(file_in, line) && std::string::npos != line.find_first_of("0123456789")) {
		std::replace(line.begin(), line.end(), '\t', ' ');
		file_out << line << std::endl;
		//std::vector<std::string> test = split(line, ' ');
	}
	file_out.close();
	file_out.open(location + "elements.txt");
	file_out << "!Elements.txt" << std::endl;
	//skip 21 lines
	for (int i = 0; i < 20; ++i) {
		getline(file_in, line);
	}
	//get vertices and nodes that makeup the elements
	facetContainer facc;
	while (getline(file_in, line) && std::string::npos != line.find_first_of("0123456789")
		&& getline(file_in_6, line2) && std::string::npos != line2.find_first_of("0123456789")
		&& getline(file_in_10, line3) && std::string::npos != line3.find_first_of("0123456789")) {
		std::replace(line.begin(), line.end(), '\t', ' ');
		std::replace(line2.begin(), line2.end(), '\t', ' ');
		std::replace(line3.begin(), line3.end(), '\t', ' ');

		std::vector<std::string> sample = split(line, ' ');
		std::vector<std::string> sample2 = split(line2, ' ');
		std::vector<std::string> sample3 = split(line3, ' ');
		//first entry is the index
//**********************************************************************************************
		int this_element = std::stoi(sample[0]);
		int this_pml_type = std::stoi(sample2[3]);
		file_out << "element: " << sample[0] << " " << sample2[3] << " " << sample[1] << std::endl;
		file_out << "expansion: " << sample3[1] << " " << sample3[2] << " " << sample3[3] << std::endl;
		file_out << "quadrature: " << sample3[4] << " " << sample3[5] << " " << sample3[6] << std::endl;

		//**********************************************************************************************
		file_out << "vertices: ";
		std::vector<unsigned int> vertices_temp;
		//for the u,v,w designation portion (which is wrong), for order 2 the corners are located at 0,2,6,8,18,20,24,26
		/*for (int i = 2; i < 10; ++i) {
			vertices_temp.push_back(std::stoi(sample[i]));
			file_out << sample[i] << " ";
		}*/
		if (sample[1] == "2") {
			file_out << sample[2] << " " << sample[4] << " " << sample[8] << " " << sample[10] << " "
				<< sample[20] << " " << sample[22] << " " << sample[26] << " " << sample[28];

			vertices_temp.push_back(std::stoi(sample[2]));
			vertices_temp.push_back(std::stoi(sample[4]));
			vertices_temp.push_back(std::stoi(sample[8]));
			vertices_temp.push_back(std::stoi(sample[10]));
			vertices_temp.push_back(std::stoi(sample[20]));
			vertices_temp.push_back(std::stoi(sample[22]));
			vertices_temp.push_back(std::stoi(sample[26]));
			vertices_temp.push_back(std::stoi(sample[28]));
			//get nodes, goes from rest of the entries on the line, starting at 11
			file_out << std::endl << "nodes:";
			for (unsigned int i = 3; i < int(sample.size()); ++i) {
				if (i != 2 && i != 4 && i != 8 && i != 10 && i != 20 && i != 22 && i != 26 && i != 28) {
					file_out << " " << sample[i];
				}
			}
			file_out << std::endl;
			file_out << "all:";
			for (int i = 2; i < int(sample.size()); ++i) {
				file_out << " " << sample[i];
			}
			file_out << std::endl;
		}
		else if (sample[1] == "1"){
			file_out << sample[2] << " " << sample[3] << " " << sample[4] << " " << sample[5] << " "
				<< sample[6] << " " << sample[7] << " " << sample[8] << " " << sample[9];

			vertices_temp.push_back(std::stoi(sample[2]));
			vertices_temp.push_back(std::stoi(sample[3]));
			vertices_temp.push_back(std::stoi(sample[4]));
			vertices_temp.push_back(std::stoi(sample[5]));
			vertices_temp.push_back(std::stoi(sample[6]));
			vertices_temp.push_back(std::stoi(sample[7]));
			vertices_temp.push_back(std::stoi(sample[8]));
			vertices_temp.push_back(std::stoi(sample[9]));
			//get nodes, goes from rest of the entries on the line, starting at 11
			file_out << std::endl << "nodes:";
			for (unsigned int i = 0; i < int(sample.size()); ++i) {
				if (i != 2 && i != 4 && i != 8 && i != 10 && i != 20 && i != 22 && i != 26 && i != 28) {
					file_out << " " << sample[i];
				}
			}
			file_out << std::endl;
			file_out << "all:";
			for (int i = 2; i < int(sample.size()); ++i) {
				file_out << " " << sample[i];
			}
			file_out << std::endl;
		}
		else {
			int M = std::stoi(sample[1]);
			for (int iM = 0; iM < 8; ++iM) {
				int offset = 0;
				if (iM == 1) offset = M;
				else if (iM == 2) offset = M + M * M;
				else if (iM == 3) offset = 2*M +  M * M;
				else if (iM == 4) offset = M * (M + 1)*(M + 1);
				else if (iM == 5) offset = M * (M + 1)*(M + 1) + M;
				else if (iM == 6) offset = M * (M + 1)*(M + 1) + M + M*M;
				else if (iM == 7) offset = M * (M + 1)*(M + 1) + 2*M + M*M;
				file_out << sample[offset + 2] << " ";
				vertices_temp.push_back(std::stoi(sample[offset + 2]));
			}
			file_out << std::endl << "nodes:";
			int offset1 = M + 2;
			int offset2 = M + M * M + 2;
			int offset3 = 2*M + M * M + 2;
			int offset4= M * (M + 1)*(M + 1) + 2;
			int offset5 = M * (M + 1)*(M + 1) + M + 2;
			int offset6 = M * (M + 1)*(M + 1) + M + M * M + 2;
			int offset7 = M * (M + 1)*(M + 1) + 2*M + M * M + 2;
			for (unsigned int i = 3; i < int(sample.size()); ++i) {
				if (i != offset1 && i != offset2 && i != offset3 && i != offset4 && i != offset5 && i != offset6 && i != offset7) {
					file_out << " " << sample[i];
				}
			}
			file_out << std::endl;
			file_out << "all:";
			for (int i = 2; i < int(sample.size()); ++i) {
				file_out << " " << sample[i];
			}
			file_out << std::endl;
		}
		//********************************************	
			//add facets to elements using vertices
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		
	
		fileFacet f1 = fileFacet(vertices_temp[0], vertices_temp[2], vertices_temp[4], vertices_temp[6], 1);
		fileFacet f2 = fileFacet(vertices_temp[1], vertices_temp[3], vertices_temp[5], vertices_temp[7], 2);
		fileFacet f3 = fileFacet(vertices_temp[0], vertices_temp[1], vertices_temp[4], vertices_temp[5], 3);
		fileFacet f4 = fileFacet(vertices_temp[2], vertices_temp[3], vertices_temp[6], vertices_temp[7], 4);
		fileFacet f5 = fileFacet(vertices_temp[0], vertices_temp[1], vertices_temp[2], vertices_temp[3], 5);
		fileFacet f6 = fileFacet(vertices_temp[4], vertices_temp[5], vertices_temp[6], vertices_temp[7], 6);
		//push the facet to the container, checks if unique and returns the location then add the other params to it
		file_out << "facets: ";
		int j = facc.addFacet(f1);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << " ";
		j = facc.addFacet(f2);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << " ";
		j = facc.addFacet(f3);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << " ";
		j = facc.addFacet(f4);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << " ";
		j = facc.addFacet(f5);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << " ";
		j = facc.addFacet(f6);
		facc.setParams(j, this_element, this_pml_type);
		file_out << j + 1 << std::endl;
	}
		//********************************************	
		file_out.close();
		//*****************************************
		//output file for facets
		//count number of facets to match the fortran code (just count the ones that are shared between two elements
		
		file_out.open(location + "facets.txt");
		file_out << "!Facets.txt" << std::endl;
		std::vector<fileFacet> facet_list = facc.facets;
		int number_of_facets = 0;
		/*for (auto i = facet_list.begin(); i != facet_list.end(); ++i) {
			auto element_list = i->element_indices;
			if (element_list.size() == 2) {
				++number_of_facets;
			}
		}*/
		for (auto i = facet_list.begin(); i != facet_list.end(); ++i) {
			file_out << "facet: " << (i - facet_list.begin()+1);
			auto pml_list = i->pml_of_elements;
			int type_print;
			if (pml_list.size() < 2) {
				if (pml_list[0] == 2) type_print = 2; //pec void
				else if (pml_list[0] == 3) type_print = 3; //FEM-MoM boundary
				else {
					type_print = 2;
				}
			}
			else {
				if (pml_list[0] == 0 && pml_list[1] == 0) type_print = 0; //domain-domain
				else if (pml_list[0] == 0 && pml_list[1] == 1 || pml_list[1] == 0 && pml_list[0]==1 ) type_print = 1;
				else type_print = -1;
			}
			file_out << " " << type_print << std::endl;
			file_out << "elements: ";
			auto element_list = i->element_indices;
			if (element_list.size() < 2) {
				file_out << element_list[0] << " " << -1 << std::endl;
			}
			else {
				file_out << element_list[0] << " " << element_list[1] << std::endl;
			}
			file_out << "vertices: " << i->vertex_indices[0] << " " << i->vertex_indices[1] << " " << i->vertex_indices[2] << " " << i->vertex_indices[3] << std::endl;
			//if designations are zero and zero then is an exterior facet
			file_out << "designations: " << (i->designation_1).first << "," << (i->designation_1).second << " " 
				<< (i->designation_2).first << "," << (i->designation_2).second << std::endl;

			if (i->local_facet_indices.size() > 1)
			file_out << "local_facet_indices: " << i->local_facet_indices[0] << " " << i->local_facet_indices[1] << std::endl;
			else file_out << "local_facet_indices: " << i->local_facet_indices[0] << " " << "-1" << std::endl;
		}

		//*****************************************
		file_out.close();



		//for materials
		file_out.open(location + "materials.txt");
		file_out << "!Materials.txt" << std::endl;
		for (int i = 0; i < 46; ++i) {
			getline(file_in, line);
		}
		std::vector<Material> materials(number_of_elements+1);
		std::vector<std::string> material_lines_eps;
		std::vector <std::string> material_lines_mu;
		//material lines has 27 entries, one for each sample, where the first entry corresponds to the first node/vertex
		//append each pair of columns to the respective sample point entries
		bool is_matrix = false;
		int num_recorded = 1;
		int elem_index=-1;
		int hcode=-1;
		int icode=-1;
		int previous_element = 0;
		while (getline(file_in, line) && std::string::npos != line.find_first_of("0123456789")) {
			std::replace(line.begin(), line.end(), '\t', ' ');
			std::vector<std::string> sample = split(line, ' ');
			int current_element = std::stoi(sample[0]);
			if (current_element > previous_element) {
				previous_element = current_element;
				//adding just the basics
				materials[current_element].hcode = std::stoi(sample[1]);
				materials[current_element].icode = std::stoi(sample[2]);
				materials[current_element].pmlcode = std::stoi(sample[3]);
				materials[current_element].epsr = std::complex<double>(std::stod(sample[4]), std::stod(sample[5]));
				materials[current_element].mur = std::complex<double>(std::stod(sample[6]), std::stod(sample[7]));
				materials[current_element].Kuvw = -1;
				materials[current_element].sym = -1;

			}
			else {
				//use the first line to tell how many lines to get next
				material_lines_eps.clear();
				material_lines_mu.clear();
				materials[current_element].Kuvw = std::stoi(sample[1]);
				materials[current_element].sym = std::stoi(sample[2]);
				int size = (materials[current_element].Kuvw + 1)*(materials[current_element].Kuvw + 1)*(materials[current_element].Kuvw + 1);
				material_lines_eps.resize(size);
				material_lines_mu.resize(size);
				//Kuvw tells number of entries per line
				//sym tells number of lines (0 => 18 lines total (9 for eps), 1 => 12 lines total (6 for eps), 2 => 6 lines total (3 for eps)
				int num_lines = 9 - materials[current_element].sym * 3;
				for (int i = 0; i < num_lines; ++i) {
					getline(file_in, line);
					std::replace(line.begin(), line.end(), '\t', ' ');
					std::vector<std::string> sample = split(line, ' ');
					int z = 0;
					for (int j = 0; j < size; ++j) {
						material_lines_eps[j] += " " + sample[z] + "," + sample[z + 1];
						z += 2;
					}
				}

				for (int i = 0; i < num_lines; ++i) {
					getline(file_in, line);
					std::replace(line.begin(), line.end(), '\t', ' ');
					std::vector<std::string> sample = split(line, ' ');
					int z = 0;
					for (int j = 0; j < size; ++j) {
						material_lines_mu[j] += " " + sample[z] + "," + sample[z + 1];
						z += 2;
					}
				}
				materials[current_element].eps_list = material_lines_eps;
				materials[current_element].mu_list = material_lines_mu;
			}
			
			

		
		}
		//printing now
		
		for (auto i = ++materials.begin(); i < materials.end(); ++i) {
			file_out << (i - materials.begin()) << " " << i->hcode << " " << i->icode << " " << i->pmlcode << std::endl;
			file_out << i->Kuvw << " " << i->sym << std::endl;
			if (i->Kuvw > -1) { 
				file_out << "eps:" << std::endl;
					for (auto j = (i->eps_list.begin()); j != i->eps_list.end(); ++j) {
						file_out << *j << std::endl;
					}
					file_out << "mu:" << std::endl;
					for (auto j = (i->mu_list.begin()); j != i->mu_list.end(); ++j) {
						file_out << *j << std::endl;
					}
					//file_out << std::endl;
			}
			else {
				file_out << "eps: " << std::endl<<" " <<i->epsr.real() << "," << i->epsr.imag() << std::endl;
				file_out << "mu: " << std::endl << " "<<i->mur.real() << "," << i->mur.imag() << std::endl;
			}
		}
		file_out.close();
		file_out.open(location + "boundary.txt");
		file_out << "!boundary.txt";
		while (getline(file_in_8, line) && std::string::npos != line.find_first_of("0123456789")) {
			file_out << std::endl;
			std::replace(line.begin(), line.end(), '\t', ' ');
			std::vector <std::string> input_line = split(line, ' ');
			//file_out << line << std::endl;
			int v1_bc = std::stoi(input_line[1]);
			int v2_bc = std::stoi(input_line[2]);
			int v3_bc = std::stoi(input_line[3]);
			int v4_bc = std::stoi(input_line[4]);
			int vin[4] = {v1_bc, v2_bc, v3_bc, v4_bc};
			std::sort(std::begin(vin), std::end(vin));
			for (auto fac_it = facc.facets.begin(); fac_it != facc.facets.end(); ++fac_it) {
				int v1 = fac_it->vertex_indices[0];
				int v2 = fac_it->vertex_indices[1];
				int v3 = fac_it->vertex_indices[2];
				int v4 = fac_it->vertex_indices[3];
				int vin2[4] = {v1,v2,v3,v4};
				std::sort(std::begin(vin2), std::end(vin2));

				if (vin[0] == vin2[0] && vin[1] == vin2[1] && vin[2] == vin2[2] && vin[3] == vin2[3]) {
					file_out << input_line[0] << " " << fac_it - facc.facets.begin() << " " << input_line[5];
				}
			}
		}
		file_out.close();
		return 1;
	}
