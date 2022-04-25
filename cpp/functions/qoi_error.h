#pragma once
//qoi_error file
#include "../structure/Domain.h"

#include "../structure/Scatter.h"
namespace qoi {
	void qoi_error(Domain& dom_h, Domain& dom_l, std::string mesh_name);
}