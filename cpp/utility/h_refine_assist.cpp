#include "h_refine_assist.h"

unsigned int facetContainer::addFacet(fileFacet f1) {
	int v1;
	int v2;
	fileFacet f1_temp = f1;
	std::sort((f1_temp.vertex_indices).begin(), (f1_temp.vertex_indices).end());
	for (auto i = facets.begin(); i != facets.end(); ++i) {

		fileFacet f_orig_temp = *i;
		if (f_orig_temp.element_indices.size() > 2) continue;
		std::sort((f_orig_temp.vertex_indices).begin(), (f_orig_temp.vertex_indices).end());

		if (f_orig_temp.vertex_indices[0] == f1_temp.vertex_indices[0]
			&& f_orig_temp.vertex_indices[1] == f1_temp.vertex_indices[1]
			&& f_orig_temp.vertex_indices[2] == f1_temp.vertex_indices[2]
			&& f_orig_temp.vertex_indices[3] == f1_temp.vertex_indices[3]) {
			//make the ept table entry
			int face_number = f1.local_facet_indices[0];
			//check the permutations from that face
			Tables table;
			std::vector<std::vector<int>> fac_table = table.face_table[face_number - 1];
			std::vector<std::vector<int>> fac_table_i = table.face_table[i->local_facet_indices[0] - 1];
			std::vector<std::vector<int>> fac_table_rel = table.face_table_rel_ind;
			bool not_found = true;
			for (auto permu = fac_table_rel.begin(); permu != fac_table_rel.end(); ++permu) {
				if (f1.vertex_indices[(*permu)[0]] == i->vertex_indices[0]) {
					if (f1.vertex_indices[(*permu)[1]] == i->vertex_indices[1]) {
						//v1 = (*permu)[4];
						//v2 = (*permu)[5];
						v1 = fac_table[int(permu - fac_table_rel.begin())][4];
						v2 = fac_table[int(permu - fac_table_rel.begin())][5];
						not_found = false;
						break;
					}
					else {
						//	v1 = (*(++permu))[4];
						//	v2 = (*permu)[5];
						v1 = fac_table[int(permu - fac_table_rel.begin()) + 1][4];
						v2 = fac_table[int(permu - fac_table_rel.begin()) + 1][5];
						not_found = false;
						break;
					}
				}
				else ++permu;

			}
			
			i->local_facet_indices.push_back(f1.local_facet_indices[0]);
			i->designation_1 = std::make_pair(fac_table_i[0][4], fac_table_i[0][5]);
			i->designation_2 = std::make_pair(v1, v2);
			return (i - facets.begin());
		}
	}
	//never found a match
	facets.push_back(f1);
	return (facets.size() - 1);
}

void facetContainer::setParams(unsigned int indx, unsigned int element_idx, int pml_type) {
	facets[indx].element_indices.push_back(element_idx);
	facets[indx].pml_of_elements.push_back(pml_type);
}

void h_refine_assist::write_refine(std::string mesh_name) {
	//new geom file
	std::string line;
	std::string line_basic;
	std::ifstream basic_in;
	std::ofstream basic_out;
	basic_in.open("../Basic_refined.dat");
	std::getline(basic_in, line);
	getline(basic_in, line_basic);
	std::vector<std::string> basic = functions::split(line_basic, ' ');
	int number_of_elements = std::stoi(basic[2]);
	std::ifstream geom_in;
	std::ofstream geom_out;
	geom_in.open("../Geometry_refined.dat");
	getline(geom_in, line);
	geom_out.open("../exampleFiles/" + mesh_name + "/geometry_refined.txt");
	geom_out << "!Geom file" << std::endl;
	int counter = 1;
	while (getline(geom_in, line)) {
		geom_out << counter << " " << line << std::endl;
		++counter;
	}
	geom_in.close();
	geom_out.close();

	//new elements file
	std::ifstream elem_in;
	std::ofstream elem_out;
	elem_in.open("../Element_refined.dat");
	getline(elem_in, line);
	elem_out.open("../exampleFiles/" + mesh_name + "/elements_refined.txt");
	elem_out << "!Elem file" << std::endl;
	facetContainer facc;
	std::vector<int> vertices_temp;
	int elem_order = 0;
	int this_element = 0;
	int this_pml_type = 0;
	//std::vector<int> 
	while (getline(elem_in, line)) {
		std::vector<std::string> sample = functions::split(line, ' ');
		if (line[0] == 'I') {
			vertices_temp = std::vector<int>(8);
			elem_out << "element: " << sample[1] << " " << sample[2] << " " << sample[3] << std::endl;
			elem_order = std::stoi(sample[3]);
			this_element = std::stoi(sample[1]);
			this_pml_type = std::stoi(sample[2]);
		}
		else if (line[0] == 'B') {
			elem_out << "expansion: " << sample[1] << " " << sample[2] << " " << sample[3] << std::endl;
		}
		else if (line[0] == 'Q') {
			elem_out << "quadrature: " << sample[1] << " " << sample[2] << " " << sample[3] << std::endl;
		}
		else if (line[0] == 'N') {
			elem_out << "vertices: ";
			if (elem_order == 2) {
				elem_out << sample[1] << " " << sample[3] << " " << sample[7] << " " << sample[9] << " "
					<< sample[19] << " " << sample[21] << " " << sample[25] << " " << sample[27] << std::endl;

				vertices_temp[0] = std::stoi(sample[1]); 
				vertices_temp[1] = std::stoi(sample[3]); 
				vertices_temp[2] = std::stoi(sample[7]); 
				vertices_temp[3] = std::stoi(sample[9]);
				vertices_temp[4] = std::stoi(sample[19]);
				vertices_temp[5] = std::stoi(sample[21]);
				vertices_temp[6] = std::stoi(sample[25]);
				vertices_temp[7] = std::stoi(sample[27]);

				elem_out << "nodes:";
				for (unsigned int i = 1; i < int(sample.size()); ++i) {
					if (i != 1 && i != 3 && i != 7 && i != 9 && i != 19 && i != 21 && i != 25 && i != 27) {
						elem_out << " " << sample[i];
					}
				}

 			}
			else if (elem_order == 1) {
				elem_out << sample[1] << " " << sample[2] << " " << sample[3] << " " << sample[4] << " "
					<< sample[5] << " " << sample[6] << " " << sample[7] << " " << sample[8] << std::endl;

				vertices_temp[0] = std::stoi(sample[1]);
				vertices_temp[1] = std::stoi(sample[2]);
				vertices_temp[2] = std::stoi(sample[3]);
				vertices_temp[3] = std::stoi(sample[4]);
				vertices_temp[4] = std::stoi(sample[5]);
				vertices_temp[5] = std::stoi(sample[6]);
				vertices_temp[6] = std::stoi(sample[7]);
				vertices_temp[7] = std::stoi(sample[8]);

				elem_out << "nodes: " << vertices_temp[0];
			}
			else {
				std::cout << "GLOBAL WARNING: unsupported elem order!!!!!!!!";
			}
			elem_out << std::endl << "all:";
			for (int i = 1; i < sample.size(); ++i) {
				elem_out << " " << sample[i];
			}
			elem_out << std::endl;

			//now build facets
			fileFacet f1 = fileFacet(vertices_temp[0], vertices_temp[2], vertices_temp[4], vertices_temp[6], 1);
			fileFacet f2 = fileFacet(vertices_temp[1], vertices_temp[3], vertices_temp[5], vertices_temp[7], 2);
			fileFacet f3 = fileFacet(vertices_temp[0], vertices_temp[1], vertices_temp[4], vertices_temp[5], 3);
			fileFacet f4 = fileFacet(vertices_temp[2], vertices_temp[3], vertices_temp[6], vertices_temp[7], 4);
			fileFacet f5 = fileFacet(vertices_temp[0], vertices_temp[1], vertices_temp[2], vertices_temp[3], 5);
			fileFacet f6 = fileFacet(vertices_temp[4], vertices_temp[5], vertices_temp[6], vertices_temp[7], 6);

			//push the facet to the container, checks if unique and returns the location then add the other params to it
			elem_out << "facets: ";
			int j = facc.addFacet(f1);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << " ";
			j = facc.addFacet(f2);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << " ";
			j = facc.addFacet(f3);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << " ";
			j = facc.addFacet(f4);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << " ";
			j = facc.addFacet(f5);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << " ";
			j = facc.addFacet(f6);
			facc.setParams(j, this_element, this_pml_type);
			elem_out << j + 1 << std::endl;
		}
		
	}
	elem_out.close();

	std::ofstream fac_out;
	fac_out.open("../exampleFiles/" + mesh_name + "/facets_refined.txt");
	fac_out << "!facet file" << std::endl;
	std::vector<fileFacet> facet_list = facc.facets;
	int number_of_facets = 0;
	/*for (auto i = facet_list.begin(); i != facet_list.end(); ++i) {
	auto element_list = i->element_indices;
	if (element_list.size() == 2) {
	++number_of_facets;
	}
	}*/
	for (auto i = facet_list.begin(); i != facet_list.end(); ++i) {
		fac_out << "facet: " << (i - facet_list.begin() + 1);
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
			else if (pml_list[0] == 0 && pml_list[1] == 1 || pml_list[1] == 0 && pml_list[0] == 1) type_print = 1;
			else type_print = -1;
		}
		fac_out << " " << type_print << std::endl;
		fac_out << "elements: ";
		auto element_list = i->element_indices;
		if (element_list.size() < 2) {
			fac_out << element_list[0] << " " << -1 << std::endl;
		}
		else {
			fac_out << element_list[0] << " " << element_list[1] << std::endl;
		}
		fac_out << "vertices: " << i->vertex_indices[0] << " " << i->vertex_indices[1] << " " << i->vertex_indices[2] << " " << i->vertex_indices[3] << std::endl;
		//if designations are zero and zero then is an exterior facet
		fac_out << "designations: " << (i->designation_1).first << "," << (i->designation_1).second << " "
			<< (i->designation_2).first << "," << (i->designation_2).second << std::endl;

		if (i->local_facet_indices.size() > 1)
			fac_out << "local_facet_indices: " << i->local_facet_indices[0] << " " << i->local_facet_indices[1] << std::endl;
		else fac_out << "local_facet_indices: " << i->local_facet_indices[0] << " " << "-1" << std::endl;
	}
	//materials doesnt need to be redone
	std::ifstream boundary_in("../Boundary_refined.dat");
	std::getline(boundary_in, line);
	std::ofstream boundary_out;
	boundary_out.open("../exampleFiles/" + mesh_name + "/boundary_refined.txt");
	boundary_out << "!boundary.txt";
	while (getline(boundary_in, line)) {
		boundary_out << std::endl;
		std::replace(line.begin(), line.end(), '\t', ' ');
		std::vector <std::string> input_line = functions::split(line, ' ');
		//file_out << line << std::endl;
		int v1_bc = std::stoi(input_line[1]);
		int v2_bc = std::stoi(input_line[2]);
		int v3_bc = std::stoi(input_line[3]);
		int v4_bc = std::stoi(input_line[4]);
		int vin[4] = { v1_bc, v2_bc, v3_bc, v4_bc };
		std::sort(std::begin(vin), std::end(vin));
		for (auto fac_it = facc.facets.begin(); fac_it != facc.facets.end(); ++fac_it) {
			int v1 = fac_it->vertex_indices[0];
			int v2 = fac_it->vertex_indices[1];
			int v3 = fac_it->vertex_indices[2];
			int v4 = fac_it->vertex_indices[3];
			int vin2[4] = { v1,v2,v3,v4 };
			std::sort(std::begin(vin2), std::end(vin2));

			if (vin[0] == vin2[0] && vin[1] == vin2[1] && vin[2] == vin2[2] && vin[3] == vin2[3]) {
				boundary_out << input_line[0] << " " << fac_it - facc.facets.begin() << " " << input_line[5];
			}
		}
	}
	boundary_out.close();
	//make controls file

	std::ofstream control_file;
	control_file.open("../exampleFiles/" + mesh_name + "/control.txt");
	control_file << "Controls file" << std::endl;
	control_file << "Adjoint: false" << std::endl << "Readprevious: false" << std::endl << "elements file name: elements_refined.txt" << std::endl << "facets file name: facets_refined.txt" << std::endl << "geometry file name: geometry_refined.txt" << std::endl
		<< "materials file name: materials_refined.txt" << std::endl << "conditions file name: boundary_refined.txt"
		<< std::endl << "monostatic scattering file name: scattering.txt" << std::endl << "input_unknowns_global file name: unknonwsAdjoint.dat"
		<< std::endl << "basic_structure: " << line_basic << std::endl;

		//done with file stuff
}