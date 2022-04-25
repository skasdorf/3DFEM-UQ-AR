#pragma once
#include "../structure/Domain.h"
#include "../file_converter/fileFacet.h"

namespace fileIO {
	void read_geom(Domain& dom, std::string geom_file);
	void read_basic(Domain& dom, std::string basic_file);
	void read_elements(Domain& dom, std::string element_file);
	void read_material(Domain& dom, std::string material_file);
	void read_boundary(Domain& dom, std::string geometry);
}