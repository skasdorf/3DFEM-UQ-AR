#ifndef FACET_H
#define FACET_H

#include <vector>

class Facet {
public:
	int index; //index of the facet
	int type; // pml or pmc boundary etc
	std::pair<int, int> element_indices;
	std::vector<int> vertices = std::vector<int>(4);
	std::pair<int, int> designations1;
	std::pair<int, int> designations2;
	std::pair<int, int> local_facet_indices;
	//for boundary_condition:
/* -1 = PEC,
	0 = default boundary,

	5 = PML Boundary(face between PML element and normal Element)
*/
	int boundary_condition = 0;
	
	Facet();
	std::vector<std::vector<int>> MAKEEPT();
};

#endif