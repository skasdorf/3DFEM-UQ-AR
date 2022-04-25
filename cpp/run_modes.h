#pragma once
#include <iostream>
#include "structure/Domain.h"
#include "utility/sparseSolver.h"
#include "utility/fileIOCon.h"
#include "utility/Plotter.h"
#include "structure\element_error.h"
#include "utility\refinement.h"
#include "utility\h_refine_assist.h"
#include <unordered_map>
#include "utility\refinement_tests.h"
#include <chrono>
#include <thread>
#include "postprocessing\postprocessing.h"
#include <random>

namespace run {
	void _standard(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, double freq);
	void _RHS_compute(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results);
	void convert_and_run(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, int ref_dex);
	void _element_error(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double tol);
	void _standard_set(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, std::vector<int>& orders);
	void _element_error_AF(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double tot_error, double frac);
	void _element_error_MSG(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double tot_error, double frac);
	void _standard_custom_materials(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, std::vector<std::complex<double>>& epsr_values);
}