#pragma once
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include "../utility/functions.h"
#include <Eigen\Core>
#include "../structure/Domain.h"
#include "../functions/unitVectorsM.h"
#include "../utility/additional_basis.h"
#include "../structure/Element.h"

namespace refinement {
	std::vector<int> pick_elements_to_refine(bool p_refine, std::string qoi_errors_file, const int& n, std::vector<int>& coarsers, const int& p); // true = use p refine, false = h refine
																						
	std::complex<double> add_qoi_error(std::string qoi_error, std::complex<double> qoi);
	std::vector<int> pick_elements_to_refine(bool p_refine, std::string qoi_errors_file, const int& n, std::vector<int>& coarsers, const int& p); // true = use p refine, false = h refine
	std::vector<int> magnitude_refine(bool p_refine, std::string qoi_errors_file, const int& n);
	void p_refine(std::vector<int> refine_elements, Domain& dom, bool neigh, int const& increment);
	void track_refine(std::vector<int> refine_elements, Domain& dom);
	void h_refine2(std::vector<int> refine_elements, Domain& dom, double factor);
	void basis_add(Domain& dom, std::vector<int> refine_elements, int increaser);
	void adjacent_check(Domain& dom);
	void basis_check(Domain& dom);
	void h_refineFin(std::vector<int> refine_elements, Domain& dom, double factor);
	void h_refine_1st_order(std::vector<int> refine_elements, Domain& dom, double factor);
	void ga_refine(Domain& dom, std::vector<int>& refine_elements);
	void set_refine(Domain & dom, std::vector<int>& refine_elements, int higher);
	std::vector<int> magnitude_refine_prop(std::vector<std::complex<double>>& qoi_errors, const double& tol);
	std::vector<int> msg_refine_prop(std::vector<std::complex<double>>& qoi_errors, const double& frac);
}

struct indexer {
	/*std::vector<int> e0 = { 1,55,28,4,58,31,7,61,34,10,64,37,3,67,40,16,70,70,43,19,73,46,22,76,49,25,79,52 };
	std::vector<int> e1 = { 30,57,3,33,60,6,36,63,9,39,66,12,42,63,15,45,72,18,48,75,21,51,78,24,54,81,27 };
	std::vector<int> e2 = { 1,2,3,55,56,57,28,29,30,10,11,12,64,65,66,37,38,39,19,20,21,73,74,75,46,47,48 };
	std::vector<int> e3 = { 34,35,363,61,62,63,7,8,9,43,44,45,70,71,72,16,17,18,52,53,54,79,80,81,25,26,27 };
	std::vector<int> e4 = { 1,2,3,4,5,66,7,8,9,55,56,57,58,59,60,61,62,63,28,29,30,31,32,33,34,35,36 };
	std::vector<int> e5 = { 46,47,48,49,50,51,52,53,54,73,74,75,76,77,78,79,80,81,19,20,21,22,23,24,25,26,27 };
	std::vector<int> e_center = { 28,29,30,31,32,33,34,35,36,37,38,39,40,14,42,43,44,45,46,47,48,49,50,51,52,53,54 };*/
	/*std::vector<std::vector<int>> indices = { { 1,55,28,4,58,31,7,61,34,10,64,37,13,67,40,16,70,43,19,73,46,22,76,49,25,79,52 },
	{ 30,57,3,33,60,6,36,63,9,39,66,12,42,63,15,45,72,18,48,75,21,51,78,24,54,81,27 } ,
	{ 1,2,3,55,56,57,28,29,30,10,11,12,64,65,66,37,38,39,19,20,21,73,74,75,46,47,48 },
	{ 34,35,36,61,62,63,7,8,9,43,44,45,70,71,72,16,17,18,52,53,54,79,80,81,25,26,27 },
	{ 1,2,3,4,5,6,7,8,9,55,56,57,58,59,60,61,62,63,28,29,30,31,32,33,34,35,36 },
	{ 46,47,48,49,50,51,52,53,54,73,74,75,76,77,78,79,80,81,19,20,21,22,23,24,25,26,27 },
	{ 28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54 } };*/
	/*std::vector<std::vector<int>> indices = { {1,54,28,4,57,31,7,60,34,10,63,37,13,66,40,16,68,42,19,71,45,22,74,48,25,77,51},
	{30,56,3,33,59,6,36,62,9,39,65,12,41,67,15,44,70,18,47,73,21,50,76,24,53,79,27},
	{1,2,3,54,55,56,28,29,30,10,11,12,63,64,65,37,38,39,19,20,21,71,72,73,45,46,47},
	{34,35,36,60,61,62,7,8,9,42,43,44,68,69,70,16,17,18,51,52,53,77,78,79,25,26,27},
	{1,2,3,4,5,6,7,8,9,54,55,56,57,58,59,60,61,62,28,29,30,31,32,33,34,35,36},
	{45,46,47,48,49,50,51,52,53,71,72,73,74,75,76,77,78,79,19,20,21,22,23,24,25,26,27},
	{28,29,30,31,32,33,34,35,36,37,38,39,40,14,41,42,43,44,45,46,47,48,49,50,51,52,53} };*/

	std::vector<std::vector<int>> indices = { { 1,54,28,4,57,31,7,60,34,10,63,37,13,66,40,16,68,42,19,71,45,22,74,48,25,77,51},{9,62,36,6,59,33,3,56,30,18,70,44,15,67,41,12,65,39,27,79,53,24,76,50,21,73,47 },{ 1,2,3,54,55,56,28,29,30,10,11,12,63,64,65,37,38,39,19,20,21,71,72,73,45,46,47 },{ 9,8,7,62,21,60,36,35,34,18,17,16,70,69,68,44,43,42,27,26,25,79,78,77,53,52,51 },{ 1,2,3,4,5,6,7,8,9,54,55,56,57,58,59,60,61,62,28,29,30,31,32,33,34,35,36 },{ 27,26,25,24,23,22,21,20,19,79,78,77,76,75,74,73,72,71,53,52,51,50,49,48,47,46,45 },{ 28,29,30,31,32,33,34,35,36,37,38,39,40,14,41,42,43,44,45,46,47,48,49,50,51,52,53 } };
	std::vector < std::vector<int>> first_indices = { {1,9,3,11,5,13,7,15}, {10,2,12,4,14,6,16,8}, {1,2,9,10,5,6,13,14}, {11,12,3,4,15,16,7,8}, {1,2,3,4,9,10,11,12}, {13,14,15,16,5,6,7,8}, {9,10,11,12,13,14,15,16} };

};