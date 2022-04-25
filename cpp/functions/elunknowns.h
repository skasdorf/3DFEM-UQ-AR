//calculate my elements
#ifndef ELUNKNOWNS_H
#define ELUNKNOWNS_H
#include "../structure/Domain.h"
#include "../utility/constants.h"
namespace elunknowns {
	int calculate_my_elements(Domain& doma) {
		int elUnknownDis = 0;
		for (auto e = doma.elements.begin(); e != doma.elements.end(); ++e) {
			elUnknownDis += e->disconnected_dimension*e->disconnected_dimension;
		}
		return elUnknownDis;
	}
};
#endif ELUNKNOWNS_H