
#include "StructureControl.h"

StructureControl::StructureControl(){

}
//TODO: add number of boundary conditions (nBCs) to file
void StructureControl::readControlFile(std::string mesh_file_in) {

    std::string line;
    std::ifstream file;

    //read geometry with the format: idx,x,y,z
    file.open("../exampleFiles/"+ mesh_file_in+ "/control.txt");
    //read first line and discard
    std::getline(file,line);
    while(std::getline(file,line)){
        std::string attribute = functions::split(line, ':')[0];
        std::vector<std::string> tokens = functions::split(functions::split(line, ':')[1],' ');
        if(attribute[0] == 'A') {
            useAdjoint = false;
            if(tokens[0][0] == 't') {
                useAdjoint = true;
            }
        }
        else if(attribute[0] == 'R') {
            readPrevious = false;
            if(tokens[0][0] == 't') {
                readPrevious = true;
            }
        }
		else if (attribute[0] == 'b') {
			mode_of_operation = std::stoi(tokens[0]);
			nodes = std::stoi(tokens[1]);
			num_elem = std::stoi(tokens[2]);
			nbc = std::stoi(tokens[3]);
			ngen = std::stoi(tokens[4]);
			numberOfWaves = std::stoi(tokens[5]);
			nports = std::stoi(tokens[6]);
			fstart = std::stoi(tokens[7]);
			fstop = std::stoi(tokens[8]);
			nfr = std::stoi(tokens[9]);
		}
        else {
            fileNames[attribute] = tokens[0];
        }
       

    }
    file.close();

}
void StructureControl::calculate_tetafi() {

}