#ifndef CPP_STRUCTURECONTROL_H
#define CPP_STRUCTURECONTROL_H

#include "../utility/functions.h"
#include <map>

class StructureControl{
public:
    bool useAdjoint;
    bool readPrevious;
    std::map<std::string, std::string> fileNames;
	int mode_of_operation;
	int nodes;
	int num_elem;
	int nbc;
	int ngen;
	int numberOfWaves;
	int nports;
	int fstart;
	int fstop;
	int nfr;

    StructureControl();
    void readControlFile(std::string mesh_file_in);
	void calculate_tetafi();
};


#endif //CPP_STRUCTURECONTROL_H
