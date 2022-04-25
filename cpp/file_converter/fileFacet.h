
//file facet class
//designation type: needs to be looked at more
//We need to keep track of both designations, plus keep track of all the possible permutations of each face
//As it is currently, local_facet_indices lists the local facet of the index (1-6) to denote what sort of face, then from this I believe we can later
//exgtract the desired designations from the facet based on the required permutations
/*Facet1 = N1,2,3,4
Facet2 = N1,2,5,6
Facet3 = N1,3,5,7
Facet4 = N2,4,6,8
Facet5 = N3,4,7,8
Facet6 = N5,6,7,8

Where N1,...,N8 come from the ordering in the input file (N1 = -u, -v, -w; N2 = +u, -v, -w; ...)

in short, describes if the facet is a +-u, +-v, +-w facet
*/
#include "../structure/Tables.h"
#include <assert.h>
class fileFacet {
public:
	std::vector<unsigned int> vertex_indices;
	std::vector<unsigned int> element_indices;
	std::vector<int> pml_of_elements;
	std::vector<int> local_facet_indices;
	std::pair<int, int> designation_1;
	std::pair<int, int> designation_2;
	fileFacet(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, int local_facet_index) {
		vertex_indices.push_back(v1);
		vertex_indices.push_back(v2);
		vertex_indices.push_back(v3);
		vertex_indices.push_back(v4);
		local_facet_indices.push_back(local_facet_index);
	}
	fileFacet(std::vector<unsigned int> vertex_indices) {
		this->vertex_indices = vertex_indices;
	}
	fileFacet(std::vector<unsigned int> vertex_indices, std::vector<unsigned int> element_indices, std::vector<int> pml_of_elements, std::vector<int> local_facet_indices) {
		this->vertex_indices = vertex_indices;
		this->element_indices = element_indices;
		this->pml_of_elements = pml_of_elements;
		this->local_facet_indices = local_facet_indices;

	}
};

class facetContainer {
public:
	std::vector<fileFacet> facets;
	unsigned int addFacet(fileFacet f1) {
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
				std::vector<std::vector<int>> fac_table = table.face_table[face_number-1];
				std::vector<std::vector<int>> fac_table_i = table.face_table[i->local_facet_indices[0]-1];
				std::vector<std::vector<int>> fac_table_rel = table.face_table_rel_ind;
				bool not_found = true;
				for (auto permu = fac_table_rel.begin(); permu != fac_table_rel.end(); ++permu) {
					if (f1.vertex_indices[(*permu)[0]] == i->vertex_indices[0]){
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
							v1 = fac_table[int(permu - fac_table_rel.begin())+1][4];
							v2 = fac_table[int(permu - fac_table_rel.begin())+1][5];
							not_found = false;
							break;
						}
					}
					else ++permu;

				}
				assert(!not_found && "Topological error! No matching permutation!");
				i->local_facet_indices.push_back(f1.local_facet_indices[0]);
				i->designation_1 = std::make_pair(fac_table_i[0][4], fac_table_i[0][5]);
				i->designation_2 = std::make_pair(v1,v2);
				return (i - facets.begin());
			}
		}
		//never found a match
		facets.push_back(f1);
		return (facets.size() - 1);
	}
	void setParams(unsigned int indx, unsigned int element_idx, int pml_type) {
		facets[indx].element_indices.push_back(element_idx);
		facets[indx].pml_of_elements.push_back(pml_type);
		//facets[indx].local_facet_indices.push_back(local_facet_index);
	}
	std::vector<fileFacet> getFacets() {
		return facets;
	}
};