#include "fileIO2.h"

void fileIO::read_geom(Domain & dom, std::string geom_file)
{
	//Line to store read data and input file stream for reading file
	std::string line;
	std::ifstream file(geom_file);
	//read first line and discard
	std::getline(file, line);
	while (std::getline(file, line)) {
		std::vector<std::string> coord = functions::split(line, ' ');
		//coord[0] just the index, doesn't matter since they'll be in order anyway
		(dom.nodes).push_back(Point(std::stod(coord[0]), std::stod(coord[1]), std::stod(coord[2])));
	}
	file.close();
}

void fileIO::read_elements(Domain & dom, std::string element_file)
{

}

void fileIO::read_material(Domain & dom, std::string material_file)
{
	std::string line;
	std::ifstream file(material_file);
	std::vector< std::vector<std::complex<double>>> epsr_storage;
	std::vector< std::vector<std::complex<double>>> mur_storage;
	std::vector< std::complex<double>> epsr_temp;
	std::vector<std::complex<double>> mur_temp;
	//getline(file, line); //discard first line
	bool first_scan = true;
	getline(file, line);
	while (getline(file, line)) {
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
		dom.elements[element_index - 1].materials.hcode = hcode;
		dom.elements[element_index - 1].materials.icode = icode;
		dom.elements[element_index - 1].materials.pmlcode = pmlcode;
		dom.elements[element_index - 1].materials.Kuvw = Kuvw;
		dom.elements[element_index - 1].materials.sym = sym;
		dom.elements[element_index - 1].materials.epsr_list = epsr_storage;
		dom.elements[element_index - 1].materials.mur_list = mur_storage;
		dom.elements[element_index - 1].materials.KuvwA = kuvwA;
		dom.elements[element_index - 1].materials.region = region;


	}
}
