#pragma once
#include "../structure/Domain.h"
#include "../utility/Plotter.h"
#include <iostream>
#include <assert.h>
inline double vector_length(Point& v1, Point& v2) {
	double x = (v1.x - v2.x);
	double y = (v1.y - v2.y);
	double z = (v1.z - v2.z);
	return sqrt(x*x + y*y + z*z);
}

class fileIOCon {
public:
	Domain* dom;
	bool higher_adjoint;
	int extra;
	fileIOCon(Domain* dom, bool higher_adjoint, int extra) {
		this->dom = dom;
		this->higher_adjoint = higher_adjoint;
		this->extra = extra;
	}
	/*void setDom(Domain& dom) {
		this->dom = dom;
	}*/
	/*void getPlottingInfo(std::string& mesh_name, PlotController plot) {

	}*/
	void output_connectivity(std::string& output_file) {
		//loop through every element saving elements it is connected to, and the average of the 8 vertices
		std::ofstream file;
		file.open(output_file);
		for (auto e = dom->elements.begin(); e != dom->elements.end(); ++e) {
			//loop through facets
			for (auto f = e->facet_indices.begin(); f != e->facet_indices.end(); ++f) {
				Facet f1 = dom->facets[*f - 1];
				f1.element_indices.first == e->index ? file << f1.element_indices.second << " " : file << f1.element_indices.first << " ";
			}
			//get average of the vertices
			double x(0), y(0), z(0);
			for (auto p = e->vertex_indices.begin(); p != e->vertex_indices.end(); ++p) {
				Point p1 = dom->nodes[*p - 1];
				x += p1.x;
				y += p1.y;
				z += p1.z;
			}
			file << x / 8.0 << " " << y / 8.0 << " " << z / 8.0 << std::endl;
		}
	}

	void file_input(std::string& mesh_name_in) {
		std::string mesh_name = mesh_name_in;
		const std::string path = "../exampleFiles/" + mesh_name_in + "/";
		// StructureControl sc = StructureControl();
		dom->sc.readControlFile(mesh_name);
		std::string a = path + dom->sc.fileNames["geometry file name"];
		std::string b = path + dom->sc.fileNames["elements file name"];
		std::string c = path + dom->sc.fileNames["facets file name"];
		std::string d = path + dom->sc.fileNames["materials file name"];
		std::string e = path + dom->sc.fileNames["conditions file name"];
		std::string f = path + dom->sc.fileNames["input_unknowns_global file name"];

		readGeometry(a);
		readElements(b);
		readFacet(c);
		readMaterial(d);
		readBoundary(e);
		if (dom->sc.readPrevious) readUknownsGlobal(e);
		readMonoStaticScattering(path + dom->sc.fileNames["monostatic scattering file name"]);
		dom->dimct = int(dom->elements.size());
	}
	//void Domain::readBasicData(std::string filename) {}
	void readUknownsGlobal(std::string& filename) {
		std::string line;
		std::ifstream file;
		file.open(filename);
		getline(file, line);
		std::vector<std::string> line_input = functions::split(line, ' ');
		dom->matDimConIn = std::stoi(line_input[0]);
		//dom->cAlphaIn.resize(dom->matDimCon);
		while (getline(file, line)) {
			line_input = functions::split(line, ' ');
			dom->cAlphaIn.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
		}
	}
	void readElements(std::string& inputFile) {

		std::string line;
		std::ifstream file;

		file.open(inputFile);

		getline(file, line);
		int elem_tracker = 0;
		while (getline(file, line)) {
			std::vector<std::string> line_input = functions::split(line, ' ');
			unsigned int index_temp = stoi(line_input[1]);
			int type = stoi(line_input[2]);
			unsigned int geometric_order = std::stoi(line_input[3]);
			getline(file, line);
			line_input = functions::split(line, ' ');
			std::vector<int> expansion_temp;
			if (higher_adjoint) {

			/*	std::cout << "Making all elements one order higher! See line 88 fileiOcon" << std::endl;
				expansion_temp.push_back(std::stoi(line_input[1])+1);
				expansion_temp.push_back(std::stoi(line_input[2])+1);
				expansion_temp.push_back(std::stoi(line_input[3])+1);
*/
				//std::cout << "Making all elements an order higher! See line 88 fileiOcon" << std::endl;
				expansion_temp.push_back(std::stoi(line_input[1])+1);
				expansion_temp.push_back(std::stoi(line_input[2])+1);
				expansion_temp.push_back(std::stoi(line_input[3])+1);

			}
			else {
				//std::cout << "Making elements lower! See line 101 fileiOcon" << std::endl;
				expansion_temp.push_back(std::stoi(line_input[1])); 
				expansion_temp.push_back(std::stoi(line_input[2]));
				expansion_temp.push_back(std::stoi(line_input[3]));
			}
			/*if (elem_tracker > -1) {
				expansion_temp[0] += 1;
				expansion_temp[1] += 1;
				expansion_temp[2] += 1;
			}*/
			++elem_tracker;
			getline(file, line);
			line_input = functions::split(line, ' ');
			std::vector<int> quadrature_temp;
			quadrature_temp.push_back(std::stoi(line_input[1]));
			quadrature_temp.push_back(std::stoi(line_input[2]));
			quadrature_temp.push_back(std::stoi(line_input[3]));
			getline(file, line);
			std::vector<int> vertices_temp;
			line_input = functions::split(line, ' ');
			for (auto i = ++line_input.begin(); i != line_input.end(); ++i) {
				vertices_temp.push_back(std::stoi(*i));
			}
			getline(file, line);
			line_input = functions::split(line, ' ');
			std::vector<int> nodes_temp;
			for (auto i = ++line_input.begin(); i != line_input.end(); ++i) {
				nodes_temp.push_back(std::stoi(*i));
			}
			//getline(file, line);
			//line_input = functions::split(line, ' ');

			getline(file, line);
			line_input = functions::split(line, ' ');
			std::vector<int> all_temp;
			for (auto i = ++line_input.begin(); i != line_input.end(); ++i) {
				all_temp.push_back(std::stoi(*i));
			}
			getline(file, line);
			line_input = functions::split(line, ' ');
			std::vector<int> facets_temp;
			for (auto i = ++line_input.begin(); i != line_input.end(); ++i) {
				facets_temp.push_back(std::stoi(*i));
			}
			/*if (type == 1) {
				expansion_temp[0] += 1;
				expansion_temp[1] += 1;
				expansion_temp[2] += 1;
			}*/

			dom->elements.push_back(Element(index_temp, type, geometric_order, vertices_temp, nodes_temp, facets_temp, quadrature_temp, expansion_temp, all_temp));
			//if (type == 1) (dom->elements.end() - 1)->refine = true;
			//compute the element diameter
			(*--(dom->elements.end())).h = vector_length(dom->nodes[vertices_temp[vertices_temp.size()-1]-1], dom->nodes[vertices_temp[0]-1]);
		}
		
		/*double maxh = dom->elements[0].h;
		for (int i = 1; i < dom->elements.size(); ++i) {
			if (dom->elements[i].h > maxh) maxh = dom->elements[i].h;
		}*/
		file.close();
	}
	void readMaterial(std::string& filename) {
		std::string line;
		std::ifstream file;
		std::vector< std::vector<std::complex<double>>> epsr_storage;
		std::vector< std::vector<std::complex<double>>> mur_storage;
		std::vector< std::complex<double>> epsr_temp;
		std::vector<std::complex<double>> mur_temp;
		//read materials
		file.open(filename);
		//getline(file, line); //discard first line
		bool first_scan = true;
		getline(file, line);
		int index = -1;
		while (getline(file, line)) {
			++index;
			/*if (index == 446) {
				int pause = 0;
			}*/
			epsr_storage.clear();
			mur_storage.clear();
			std::vector<std::string> element_line = functions::split(line, ' ');
			int element_index = std::stoi(element_line[0]);
			int hcode = std::stoi(element_line[1]);
			int icode = std::stoi(element_line[2]);
			int pmlcode = std::stoi(element_line[3]);
			getline(file, line);
			element_line = functions::split(line, ' ');
			int Kuvw = std::stoi(element_line[0]);
			int sym = std::stoi(element_line[1]);
			int size = (Kuvw + 1)*(Kuvw + 1)*(Kuvw + 1);
			getline(file, line);
			//for (int i = 0; i < size; ++i) {
			int i = 0;
			do {
				epsr_temp.clear();
				//epsr_storage.clear();
				getline(file, line);
				element_line = functions::split(line, ' ');
				for (auto j = element_line.begin(); j != element_line.end(); ++j) {
					std::vector<std::string> compl_num = functions::split(*j, ',');
					epsr_temp.push_back(std::complex<double>(std::stod(compl_num[0]), std::stod(compl_num[1])));
				}
				epsr_storage.push_back(epsr_temp);
				++i;
			} while (i < size);
			getline(file, line);
			i = 0;
			do {
				//for (int i = 0; i < size; ++i) {

				mur_temp.clear();

				//mur_storage.clear();
				getline(file, line);
				element_line = functions::split(line, ' ');
				for (auto j = element_line.begin(); j != element_line.end(); ++j) {
					std::vector<std::string> compl_num = functions::split(*j, ',');
					mur_temp.push_back(std::complex<double>(std::stod(compl_num[0]), std::stod(compl_num[1])));
				}
				mur_storage.push_back(mur_temp);
				++i;
			} while (i < size);
			int kuvwA = 0;
			int region;
			if (icode == 0) kuvwA = Kuvw;
			if ((pmlcode != 1) && (hcode == 1) && (icode == 1) && (abs(mur_storage[0][0] - 1.0) < .0000001) && (abs(epsr_storage[0][0] - 1.0) < .0000001)) region = 0;
			else if (pmlcode == 1) region = 2;
			else region = 1;
			//	dom->elements[element_index - 1].materials.push_back(Material(hcode, icode, pmlcode, Kuvw, sym, epsr_storage, mur_storage, kuvwA));
			/*std::cout << "Warning: materials are modified from original (changed region 0)" << std::endl;
			if (region == 1) {
				epsr_storage = { {2.22} };
			}*/
			dom->elements[element_index - 1].materials.hcode = hcode;
			dom->elements[element_index - 1].materials.icode = icode;
			dom->elements[element_index - 1].materials.pmlcode = pmlcode;
			dom->elements[element_index - 1].materials.Kuvw = Kuvw;
			dom->elements[element_index - 1].materials.sym = sym;
			dom->elements[element_index - 1].materials.epsr_list = epsr_storage;
			dom->elements[element_index - 1].materials.mur_list = mur_storage;
			dom->elements[element_index - 1].materials.KuvwA = kuvwA;
			dom->elements[element_index - 1].materials.region = region;


		}
	}
	void readGeometry(std::string& filename) {
		//Line to store read data and input file stream for reading file
		std::string line;
		std::ifstream file;

		//read geometry with the format: idx,x,y,z
		file.open(filename);
		//read first line and discard
		std::getline(file, line);
		while (std::getline(file, line)) {
			std::vector<std::string> coord = functions::split(line, ' ');
			//coord[0] just the index, doesn't matter since they'll be in order anyway
			
			(dom->nodes).push_back(Point(std::stod(coord[1]), std::stod(coord[2]), std::stod(coord[3])));
		}
		file.close();
	}

	void readFacet(std::string& filename) {
		std::string line;
		std::ifstream file;

		//read geometry with the format: idx,x,y,z
		file.open(filename);
		//read first line and discard
		std::getline(file, line);
		Facet newFacet = Facet();
		while (std::getline(file, line)) {
			std::string attribute = functions::split(line, ':')[0];
			std::vector<std::string> tokens = functions::split(functions::split(line, ':')[1], ' ');
			if (attribute[0] == 'f') {
				newFacet.index = std::stoi(tokens[0]);
				newFacet.type = std::stoi(tokens[1]);
			}
			else if (attribute[0] == 'e') {
				newFacet.element_indices = std::make_pair(std::stoi(tokens[0]), std::stoi(tokens[1]));
			}
			else if (attribute[0] == 'd') {
				// newFacet.designations = std::make_pair(std::stoi(tokens[0]), std::stoi(tokens[1]));
				//  facets.push_back(newFacet);
				//   newFacet = Facet();
				std::vector<std::string> desig1 = functions::split(tokens[0], ',');
				std::vector<std::string> desig2 = functions::split(tokens[1], ',');
				newFacet.designations1 = std::make_pair(std::stoi(desig1[0]), std::stoi(desig1[1]));
				newFacet.designations2 = std::make_pair(std::stoi(desig2[0]), std::stoi(desig2[1]));
			}
			else if (attribute[0] == 'v') {
				newFacet.vertices[0] = std::stoi(tokens[0]);
				newFacet.vertices[1] = std::stoi(tokens[1]);
				newFacet.vertices[2] = std::stoi(tokens[2]);
				newFacet.vertices[3] = std::stoi(tokens[3]);
			}
			else {
				newFacet.local_facet_indices = std::make_pair(std::stoi(tokens[0]), std::stoi(tokens[1]));
				dom->facets.push_back(newFacet);
				newFacet = Facet();
			}


		}
		file.close();
	}
	void readBoundary(std::string filename) {
		/*          -1 = PEC,
		0 = ABC(1st),
		1 = SYM,
		5 = PML Boundary(face between PML element and normal Element)
		*/
		std::string line;
		std::ifstream file;
		//read geometry with the format: idx,x,y,z
		file.open(filename);
		//read first line and discard
		std::getline(file, line);
		while (std::getline(file, line)) {
			std::vector < std::string> input_line = functions::split(line, ' ');
			//first is the index of the boundary condiditon
			//second is the index of the facet (starting at zero!)
			//third is the type of boundary condition
			dom->facets[std::stoi(input_line[1])].boundary_condition = std::stoi(input_line[2]);
		}
	}

	void readMonoStaticScattering(std::string filename) {

		std::string line;
		std::ifstream file;
		file.open(filename);
		std::getline(file, line);

		std::getline(file, line);
		std::vector<std::string> input_line = functions::split(line, ' ');
		dom->scatter1.thStart = std::stod(input_line[0]);
		dom->scatter1.thStop = std::stod(input_line[1]);
		dom->scatter1.nTh = std::stoi(input_line[2]);
		dom->scatter1.nPh = std::stoi(input_line[5]);
		dom->scatter1.phStart = std::stod(input_line[3]);
		dom->scatter1.phStop = std::stod(input_line[4]);
		std::getline(file, line);
		input_line = functions::split(line, ' ');
		//dom->scatter1.waveNum = std::stoi(input_line[0]);
		//dom->scatter1.theta = std::stod(input_line[1]);
		//dom->scatter1.phi = std::stod(input_line[2]);
		dom->scatter1.waveTheta = { 0, std::complex<double>(std::stod(input_line[3]), std::stod(input_line[4])) };
		dom->scatter1.wavePhi = { 0, std::complex<double>(std::stod(input_line[5]), std::stod(input_line[6])) };

	}
};