/*#define PI 3.1415926535897932385
 *#define eps0 8.85419e-12
 *
 * Domain
 *
 * Stores the structural information of the geometry in use.
 * The information is provided via txt files supplied by the control.txt.txt file
 *
 * Takes the following input:
 * 1)
 *
*/

#ifndef DOMAIN_H
#define DOMAIN_H
#include <fstream>
#include <unordered_map>
#include "../eigen-eigen-5a0156e40feb/Eigen/Core.h"
#include "../eigen-eigen-5a0156e40feb/Eigen/Sparse.h"
#include "Tables.h"
#include "Element.h"
#include <iomanip>
#include "Facet.h"
#include "Material.h"
#include "Point.h"
#include "../utility/functions.h"
#include "../utility/fileIO_functions.h"
#include "../utility/memory_functions.h"
#include "../control/StructureControl.h"
#include "Scatter.h"
#include "../utility/sparseSolver.h"
#include "../functions/unitVectorsM.h"
#include "../integrals/Integral_g.h"
#include "../integrals/Integral_c.h"
#include "../integrals/Integral_d.h"
//#include "../integrals/findIntegrals.h"
//#include "../functions/basis.h"
#include "../utility/constants.h"
//#include "../functions/sparsePacking.h"
#include "../structure/FrequencySweep.h"
#include "../functions/sensitivity.h"
#include "../structure/error_struct.h"

#include "../utility/qoi_helper.h"
#include "../utility/additional_basis.h"
#include "../integrals/BasisEval.h"


class Domain {
public:
	Domain() {

	}
	//std::vector<Element> elements;
	std::vector<std::complex<double>> cBrSparse, cGr_k0_half;
	std::vector<int> iRow, jCol;
	bool error = false;
	bool check_results;
	bool higher_order;
	std::string mesh_name;
	std::vector<Facet> facets;
	std::vector<Point> vertices;
	std::vector<Point> nodes;
	//std::vector<Material> materials;
	std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> i0_table;
	std::vector<Element> elements;
	Tables vertfac_table = Tables();
	int dimct;
	int matDimDis;
	int matDimCon;
	std::vector<std::vector<adBasis>> p_ref_basis;
	//std::vector<std::vector<Error_struct>> error_structure;
	std::vector<std::vector<std::vector<int>>> nTotal; //NOTE: Goes from 1 to N, NOT 0 to N-1!
	std::vector<std::vector<int>> ntCum; //NOTE: Goes from 1 to N, NOT 0 to N-1!
	StructureControl sc;
	std::vector<int> vectorD;
	std::vector<std::vector<int>> eiUVWijk; //NOTE: Goes from 1 to N, NOT 0 to N-1!
	std::vector<std::vector<int>> unknown_description;
	int matDimConIn;
	std::vector<std::complex<double>> cAlphaIn;
	std::vector<std::complex<double>> cAlphaForward;
	std::vector<std::complex<double>> cAlphaAdjoint;
	std::vector<std::complex<double>> cAlphaAdjointLower;
	Scatter scatter1;
	std::vector<int> vectorDUnique;
	std::vector<std::complex<double>> cAlpha;
	std::vector<std::complex<double>> cGr_eps_vec;
	std::vector<std::complex<double>> cBr_eps_vec;
	Domain(std::string mesh_name_in);
	void makeNT(); 
	void definers_defrmnls();
	int DADDR(int element_index, int face_index, int coord, int a, int b); 
	void connect_elements();
	//void Domain::make_bcct(int number_of_elements, int nBCS);
	std::vector<int> make_bcept(Facet f);
	void bc_coeff_setting(std::vector<int> bcept, Facet f);
	void set_bc_elements();
	std::vector<int> prenumunknowns();
	void fill_eiuvwijk();
	void fill_eiuvwijk_refine(int order_dif);
	void make_unknown_description();
	void vectorDSignPW(int addr1, int addr2, int signPW, int b, std::vector<int> &vectorD);
	void startMatching(int numberOfElements, int e1, int e2, std::vector<std::vector<int>>& ept, std::vector<int> & vectorD);
	bool check_q();
	void findPowersLagrange(const int& kuvw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw, matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr);
	void clear_all_submatrices(std::vector<Element>& elements);
	

	double trilagrangeoduvw(const int & m, const int& n, const int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr);
	void find_eps_mu_hermitian(const int& size1, const int& nglu, const int& nglv, const int& nglw, matrix4d < std::complex<double>>& MuRelIntInv);
	void find_eps_mu_matrix(const int& kuvw, const int& nglu, const int& nglv, const int& nglw, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr, const int& size1, const std::vector<std::vector<std::complex<double>>>& muRel, matrix4d<std::complex<double>>& muRelInt);
	/*void find_integrals(Element element, int iuvwh, int iuvw, int ih, int jh, int kh, int i, int j, int k,
		std::complex<double>& cSolA, std::complex<double>& cSolB, std::complex<double>& cSolS,
		std::complex<double>& cSolStotal, int size1, matrix4d<std::complex<double>>& MuRelInt,
		matrix4d<std::complex<double>>& MuRelIntInv, matrix4d<std::complex<double>>& EpsRelInt);*/
	void findelements_abs_sparse();
	void _RHS_ONLY();
	void Domain::SPARSE_PACKING(std::vector<int>& iRow,
		std::vector<int>&jCol, int& nNZ, int noAbcBCs,
		std::vector<dcomplex>& cArSparse,
		std::vector<dcomplex>& cBrSparse,
		std::vector<dcomplex>& cSrSparse);
	void Domain::spherical_sensitivity(const int& myElementCount, std::vector<Element>& elements, const int& matDimDis,
		const std::vector<int>& vectorD, const std::vector<std::vector<int>>& eiUVWijk, const int& matDimCon, const std::vector<std::complex<double>>& cAlphaForward,
		std::vector<std::complex<double>>& cAlphaAdjoint, std::complex<double>& cSensitivity, const int& numberOfFrequencies,
		std::vector<double>& k0List, const int& indicatorBasisType, Scatter* scatter1, const std::vector<Facet>& facets);

	
	void Domain::element_error_f(int u_order, int v_order, int w_order, Domain& dom2);
	void Domain::element_error_basis(int u_order, int v_order, int w_order, Domain& dom2);
	/*void Domain::findPowers(const int& nu, const int& nv, const int& nw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
		matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
		matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers);*/
	void Domain::_HOPS_EPS(int u_order, int v_order, int w_order, Domain& dom2);
	void Domain::element_error_improved(Domain& dom2, std::vector<std::complex<double>>& cAlphaFor, std::vector<std::complex<double>>& cAlphaAdj);
	void Domain::basis_error_improved(Domain& dom2, std::vector<std::complex<double>>& cAlphaFor, std::vector<std::complex<double>>& cAlphaAdj);
};


#endif

