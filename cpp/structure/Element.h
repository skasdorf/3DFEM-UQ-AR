//element class

/*
example element:
element: 187 0 3
expansion: 4 3 4
vertices: 6 14 908 266 17 76 162 91
nodes: 56 24 10 789 182 956 273 727 68 926 126 734 89 91 537 823 472 296 58 19 62 563 968 152 727 846 92
Quadrature: 14 13 14
facets: 10 984 234 183 8 387
*/


#ifndef ELEMENT_H
#define ELEMENT_H

#include "../utility/constants.h"
#include "Material.h"
#include "adjacent.h"

using dcomplex = std::complex<double>;

class Element {
public:
	bool refine = false;
	unsigned int index;
	int type; //0 - regular, 1 - PML
	//int region;
	int disconnected_dimension;
	int geom_order;
	int unknownsStart;
	int unknownsEnd;
	matrix2d<double> rs; //starts at 1
	int nRs;
	double h; //the diameter of the element
	adjacent adjacent_list;
	//std::vector<Error_struct> error_structure;
	std::complex<double> element_error;
	std::vector<int> vertex_indices;
	std::vector<int> node_indices;
	std::vector<int> all_indices;
	matrix3d<double> jacobian; //starts at 0
	std::vector<int> expansion; //nu, nv, nw
	std::vector<int> facet_indices;
	std::vector<int> quadrature; //nGL_u, nGL_v, nGL_w
	std::vector<int> abcList = {0,0,0,0,0,0};
	std::vector<int> scatter_boundary_list = { 0,0,0,0,0,0 };
	std::vector<int> PML_boundary_list = { 0,0,0,0,0,0 };
	std::vector<std::complex<double>> cGr_eps_el;
	Material materials; //0 is air, 1 is scatterer, 2 is PML
	matrix2d<double> fuPowersLagr;
	matrix2d<double> fvPowersLagr;
	matrix2d<double> fwPowersLagr;
	matrix4d<double> rMatrix;
	matrix2d<double> uPowers; 
	matrix2d<double> vPowers;
	matrix2d<double> wPowers;
	matrix2d<double> fuPowers;
	matrix2d<double> fvPowers;
	matrix2d<double> fwPowers;
	matrix2d<double> fpuPowers;
	matrix2d<double> fpvPowers;
	matrix2d<double> fpwPowers;
	matrix4d<double> auMatrix;
	matrix4d<double> avMatrix;
	matrix4d<double> awMatrix;
	matrix2d<dcomplex> cArPAK;
	matrix2d<dcomplex> cBrPAK;
	matrix2d<dcomplex> cSrPAK;
	matrix2d<dcomplex> cBr_HOPS;
	matrix4d<std::complex<double>> EpsRelInt;
	matrix4d<std::complex<double>> MuRelIntInv;
	matrix4d<std::complex<double>> MuRelInt;
	std::vector<double> wglu, wglv, wglw; //weights


	//Element();
	Element(unsigned int index, int type, int geom_order, std::vector<int> vertex_indices, std::vector<int> node_indices, 
		std::vector<int> facet_indices, std::vector<int> quadrature, std::vector<int> expansion, std::vector<int> all_indices) {
		this->index = index;
		this->type = type;
		this->geom_order = geom_order;
		this->vertex_indices = vertex_indices;
		this->node_indices = node_indices;
		this->facet_indices = facet_indices;
		this->quadrature = quadrature;
		this->expansion = expansion;
		this->all_indices = all_indices;
	}
	Element() {
		this->index = 0;
	}
};
#endif