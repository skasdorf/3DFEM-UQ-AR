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
#include "HOPS\hops.h"
#include "Kriging\kriging.h"
#include "Bayes\bayes.h"
void refine_main(std::vector<int>& refiners, std::complex<double>& qoi, const int& n, int& numu, char ref_type) {
	/////////////////////////
	//Run options:
	bool check_q = false;
	bool plot = false;
	bool sens = false;
	bool error = false;
	bool refine = true;
	bool basis_error = false;
	bool useAdjoint = false;
	bool higher_order = false; //right now higher_order off reduces order by 1 (3->2), on tests with order = 3
	bool check_results = false;
	//////////////////////////////////



	int extra = 0; //increase in the overall expansion order for every element
	if (higher_order) extra = 1;
	/////////////////////////////////
	std::complex<double> cSensitivity;
	std::string mesh_name = "standard";
	Domain dom(mesh_name);
	Domain dom2(mesh_name);
	//////////////////

	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;
	if (refine) {
		dom2.sc.useAdjoint = true;
		dom2.higher_order = higher_order;
	}
	else {
		dom2.sc.useAdjoint = useAdjoint;
		dom2.higher_order = false; //for non refinement, needs lower order solve
	}
	/////////////////
	if (dom.sc.useAdjoint && !check_q) dom.mesh_name = mesh_name + "/adjoint";
	if (dom2.sc.useAdjoint && !check_q) dom2.mesh_name = mesh_name + "/adjoint";
	fileIOCon fileIO(&dom, higher_order, extra);
	fileIOCon fileIO2(&dom2, dom2.higher_order, extra);
	dom.check_results = check_results;
	dom2.check_results = check_results;



	fileIO.file_input(mesh_name);
	fileIO2.file_input(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;

	dom2.sc.useAdjoint = true;
	dom2.higher_order = higher_order;

	refinement::p_refine(refiners, dom, false,1);
	refinement::p_refine(refiners, dom2, false,1);
	//refinement::adjacent_check(dom);
	//refinement::adjacent_check(dom2);

	//dom.file_input("standard");

	//fix reflections with new basis
	//refinement::basis_check(dom);
	//refinement::basis_check(dom2);

	//
	//refinement::basis_add(dom, refiners, 1);
	//refinement::basis_add(dom2, refiners, 1);
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	//dom.fill_eiuvwijk_refine(1);
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	if (ref_type == 'h') {
		refinement::h_refine2(refiners, dom, 0.3);
	}
	//solve forward and adjoint problems
	dom.findelements_abs_sparse();
	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha_.txt", true);
		plotter.plotField(mesh_name, false);
	}
	std::vector<double> rcs_list = postproc::post(dom);
	std::ofstream outFile("rcs_list.txt");
	for (const auto &rcs : rcs_list) outFile << rcs << "\n";
	//setup dom2 for refinement solve
	dom2.makeNT();
	dom2.connect_elements();
	dom2.set_bc_elements();
	dom2.fill_eiuvwijk();
	//dom2.fill_eiuvwijk_refine(1);
	std::vector<int> vectorDUnique2 = dom2.prenumunknowns();
	dom2.make_unknown_description();
	dom2.scatter1.waveSetup(dom2.matDimCon, dom2.sc.nfr, dom2.sc.fstart, dom2.sc.fstop, dom2.sc.numberOfWaves);
	dom2.scatter1.calcThetaPhi();
	dom2.scatter1.calcFrequencies();
	dom2.definers_defrmnls();
	dom2.findelements_abs_sparse();
	//get QoI
	if (higher_order == false) qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
	else qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
	numu = dom.cAlpha.size();

}
void _run_refine_RCS(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results,
	std::vector<int>& refiners, std::vector<double>& rcs_list, const int& n, int& numu, std::vector<int>& coarsers) {
	int extra = 0;
	if (higher_order) extra = 1;
	Domain dom(mesh_name);

	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;

	/////////////////
	if (dom.sc.useAdjoint) dom.mesh_name = mesh_name + "/adjoint";

	fileIOCon fileIO(&dom, higher_order, extra);
	dom.check_results = check_results;

	std::cout << "Running analysis now: " << std::endl;
	std::cout << "Adjoint: " << useAdjoint << std::endl << "Higher order basis: " << higher_order << std::endl;

	dom.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
	fileIO.file_input(mesh_name);
	

	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;

	refinement::p_refine(refiners, dom, false,1);
	if (coarsers.size() > 0) {
		std::vector<int> coarsers1(coarsers.begin(), coarsers.begin() + 5);
		std::vector<int> coarsers2(coarsers.begin() + 5, coarsers.end());
		refinement::p_refine(coarsers1, dom, false, -2);
		refinement::p_refine(coarsers2, dom, false, -1);
	}
	
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	std::string connectivity = "../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt";
	fileIO.output_connectivity(connectivity);
	dom.findelements_abs_sparse();
	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");
	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", true);
		plotter.plotField(mesh_name, false);
	}

	rcs_list = postproc::post(dom);
	numu = dom.cAlpha.size();
	//std::ofstream outFile("rcs_list.txt");
	//for (const auto &rcs : rcs_list) outFile << rcs << "\n";
}
void _run_standard(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results) {
	int extra = 0;
	if (higher_order) extra = 1;
	Domain dom(mesh_name);
	//Domain dom2(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;
	//dom2.sc.useAdjoint = useAdjoint;
	//dom2.higher_order = false;
	/////////////////
	if (dom.sc.useAdjoint) dom.mesh_name = mesh_name + "/adjoint";
	//if (dom2.sc.useAdjoint) dom2.mesh_name = mesh_name + "/adjoint";
	fileIOCon fileIO(&dom, higher_order, extra);
	//fileIOCon fileIO2(&dom2, dom2.higher_order, extra);
	dom.check_results = check_results;
	//dom2.check_results = check_results;
	std::cout << "Running analysis now: " << std::endl;
	std::cout << "Adjoint: " << useAdjoint << std::endl << "Higher order basis: " << higher_order << std::endl;

	dom.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
	fileIO.file_input(mesh_name);
	//fileIO2.file_input(mesh_name);

	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;
	//dom2.sc.useAdjoint = useAdjoint;
	//dom2.higher_order = false;
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	dom.findelements_abs_sparse();
	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");
	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", true);
		plotter.plotField(mesh_name, false);
	}
	std::vector<double> rcs_list = postproc::post(dom);
	std::ofstream outFile("rcs_list.txt");
	for (const auto &rcs : rcs_list) outFile << rcs << "\n";
	std::cout << "Num unknowns: " << dom.cAlpha.size() << std::endl;
}
void _run_element_error(std::string mesh_name) {
	bool check_q = false;
	bool plot = false;
	bool sens = false;
	bool error = true;
	bool refine = false;
	bool basis_error = false;
	bool useAdjoint = false;
	bool higher_order = true; //right now higher_order off reduces order by 1 (3->2), on tests with order = 3
	bool check_results = false;
	//////////////////////////////////



	int extra = 0;
	if (higher_order) extra = 1;
	/////////////////////////////////
	std::complex<double> cSensitivity;
	//std::string mesh_name = "standard";
	Domain dom(mesh_name);
	Domain dom2(mesh_name);
	//////////////////

	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;
	if (refine) {
		dom2.sc.useAdjoint = true;
		dom2.higher_order = higher_order;
	}
	else {
		dom2.sc.useAdjoint = useAdjoint;
		dom2.higher_order = false;
	}
	/////////////////
	if (dom.sc.useAdjoint && !check_q) dom.mesh_name = mesh_name + "/adjoint";
	if (dom2.sc.useAdjoint && !check_q) dom2.mesh_name = mesh_name + "/adjoint";
	fileIOCon fileIO(&dom, higher_order, extra);
	fileIOCon fileIO2(&dom2, dom2.higher_order, extra);
	dom.check_results = check_results;
	dom2.check_results = check_results;
	std::cout << "Running analysis now: " << std::endl;
	std::cout << "Adjoint: " << useAdjoint << std::endl << "Higher order basis: " << higher_order << std::endl
		<< "Check_q: " << check_q << std::endl << "Calc sens: " << sens << std::endl << "Calc QoI error: " << error << std::endl;
	dom.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;

	//fileIO.setDom(dom);
	fileIO.file_input(mesh_name);
	fileIO2.file_input(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.error = error;
	if (refine) {
		dom2.sc.useAdjoint = true;
		dom2.higher_order = higher_order;
	}
	else {
		dom2.sc.useAdjoint = useAdjoint;
		dom2.higher_order = false;
	}
	//dom.file_input("standard");
	if (!check_q || sens || error || refine) {
		dom.makeNT();
		dom.connect_elements();
		dom.set_bc_elements();
		dom.fill_eiuvwijk();
		std::vector<int> vectorDUnique = dom.prenumunknowns();
		dom.make_unknown_description();
		dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
		dom.scatter1.calcThetaPhi();
		dom.scatter1.calcFrequencies();
		dom.definers_defrmnls();
		if (!plot) {
			dom2.makeNT();
			dom2.connect_elements();
			dom2.set_bc_elements();
			dom2.fill_eiuvwijk();
			std::vector<int> vectorDUnique2 = dom2.prenumunknowns();
			dom2.make_unknown_description();
			dom2.scatter1.waveSetup(dom2.matDimCon, dom2.sc.nfr, dom2.sc.fstart, dom2.sc.fstop, dom2.sc.numberOfWaves);
			dom2.scatter1.calcThetaPhi();
			dom2.scatter1.calcFrequencies();
			dom2.definers_defrmnls();
		}
	}
	std::complex<double> qoi;

	if (check_q) {
		/*qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
		std::cout << "QoI (lower order): " << qoi << std::endl;*/

		std::complex<double> qoi_2 = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
		std::cout << "QoI (higher order): " << qoi_2 << std::endl;
		//	
		//	std::cout << "QoI (lower) + QoI error per element: " << qoi_error::add_qoi_error("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt", qoi) << std::endl;

	}
	else if (error) dom.element_error_f(2, 2, 2, dom2); //input the order of each element in forward solve, (WIP)
	else if (basis_error) dom.element_error_basis(2, 2, 2, dom2);
	else if (!plot) dom.findelements_abs_sparse();

	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");

	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", true);
		plotter.plotField(mesh_name, false);
	}
}
void _run_check_q(std::string mesh_name, int qoi_control) {
	std::complex<double> qoi;
	if (qoi_control == 0) {
		qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
		std::cout << "QoI (lower order): " << qoi << std::endl;
	}
	if (qoi_control == 1) {
		std::complex<double> qoi_2 = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
		std::cout << "QoI (higher order): " << qoi_2 << std::endl;
	}
	/*if (qoi_control == 2) {
		std::cout << "QoI (lower) + QoI error per element: " << qoi_error::add_qoi_error("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt", qoi) << std::endl;
	}*/

}
\



void _eval_qoi(std::string mesh_name, bool higher_order, std::complex<double>& qoi) {
	if (higher_order == false) qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
	else qoi = qoi_error::check_q("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
}
void _run_h_refine_setup(std::string mesh_name, std::vector<int>& refine_elements,  bool useAdjoint, bool higher_order) {
	int extra = 0;
	if (higher_order) extra = 1;
	Domain dom(mesh_name);
	Domain dom2(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;
	/////////////////
	if (dom.sc.useAdjoint) dom.mesh_name = mesh_name + "/adjoint";
	fileIOCon fileIO(&dom, higher_order, extra);

	fileIO.file_input(mesh_name);
	
	dom.check_results = false;
	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	//tests::test(refine_elements, dom, 0.2);
	
	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");
	if (dom.elements[0].geom_order == 2) {
		tests::split_elements(dom, refine_elements, 0.3);
		//refinement::h_refineFin(refine_elements, dom, 0.2);
	}
	else if (dom.elements[0].geom_order == 1) {
		//refinement::h_refine_1st_order(refine_elements, dom, 0.2);
	}
	h_refine_assist href;
	href.write_refine(mesh_name + "_refine");
}
void _run_h_refine(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, std::vector<double>& rcs_list) {
	int extra = 0;
	if (higher_order) extra = 1;
	Domain dom(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;
	/////////////////
	if (dom.sc.useAdjoint) dom.mesh_name = mesh_name + "/adjoint";
	fileIOCon fileIO(&dom, higher_order, extra);
	dom.check_results = check_results;
	std::cout << "Running analysis now: " << std::endl;
	std::cout << "Adjoint: " << useAdjoint << std::endl << "Higher order basis: " << higher_order << std::endl;
	dom.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
	fileIO.file_input(mesh_name);
	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	dom.findelements_abs_sparse();
	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");
	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", true);
		plotter.plotField(mesh_name, false);
	}

	/*std::vector<double> rcs_list = postproc::post(dom);
	std::ofstream outFile("rcs_list.txt");
	for (const auto &rcs : rcs_list) outFile << rcs << "\n";*/

	rcs_list = postproc::post(dom);

}

void _run_refine_GA(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results,
	std::vector<int>& refiners, std::vector<double>& rcs_list, int& numu) {
	int extra = 0;
	if (higher_order) extra = 1;
	Domain dom(mesh_name);

	dom.sc.useAdjoint = useAdjoint;
	dom.higher_order = higher_order;

	/////////////////
	if (dom.sc.useAdjoint) dom.mesh_name = mesh_name + "/adjoint";

	fileIOCon fileIO(&dom, higher_order, extra);
	dom.check_results = check_results;

	std::cout << "Running analysis now: " << std::endl;
	std::cout << "Adjoint: " << useAdjoint << std::endl << "Higher order basis: " << higher_order << std::endl;

	dom.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
	fileIO.file_input(mesh_name);


	dom.sc.useAdjoint = useAdjoint;
	dom.error = false;

	/*refinement::p_refine(refiners, dom, false, 1);
	if (coarsers.size() > 0) {
		std::vector<int> coarsers1(coarsers.begin(), coarsers.begin() + 5);
		std::vector<int> coarsers2(coarsers.begin() + 5, coarsers.end());
		refinement::p_refine(coarsers1, dom, false, -2);
		refinement::p_refine(coarsers2, dom, false, -1);
	}*/
	refinement::ga_refine(dom, refiners);
	dom.makeNT();
	dom.connect_elements();
	dom.set_bc_elements();
	dom.fill_eiuvwijk();
	std::vector<int> vectorDUnique = dom.prenumunknowns();
	dom.make_unknown_description();
	dom.scatter1.waveSetup(dom.matDimCon, dom.sc.nfr, dom.sc.fstart, dom.sc.fstop, dom.sc.numberOfWaves);
	dom.scatter1.calcThetaPhi();
	dom.scatter1.calcFrequencies();
	dom.definers_defrmnls();
	std::string connectivity = "../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt";
	fileIO.output_connectivity(connectivity);
	dom.findelements_abs_sparse();
	//fileIO.output_connectivity("../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt");
	if (plot) {
		PlotController plotter(&dom);
		plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", true);
		plotter.plotField(mesh_name, false);
	}

	rcs_list = postproc::post(dom);
	numu = dom.cAlpha.size();
	//std::ofstream outFile("rcs_list.txt");
	//for (const auto &rcs : rcs_list) outFile << rcs << "\n";
}


////main for multi_HOPS testing...
//int main() {
//	//solver inputs
//	bool check_q = false;
//	bool plot = false;
//	bool sens = false;
//	bool error = false;
//	bool refine = false;
//	bool basis_error = false;
//	bool useAdjoint = false;
//	bool higher_order = false; //right now higher_order off reduces order by 1 (3->2), on tests with order = 3
//	bool check_results = false;
//	//mesh inputs
//	std::string mesh_name = "references";
//	std::vector<std::complex<double>> epsr_list = { std::complex<double>(6.0, -1.84) };
//	std::vector<std::complex<double>> updated_qoi;
//	int i = 0;
//	HOPS::multi_HOPS_epsr(mesh_name);
//
//	return 0;
//}


int main() {
		//solver inputs
		bool check_q = false;
		bool plot = false;
		bool sens = false;
		bool error = false;
		bool refine = false;
		bool basis_error = false;
		bool useAdjoint = false;
		bool higher_order = false; //right now higher_order off reduces order by 1 (3->2), on tests with order = 3
		bool check_results = false;
		//mesh inputs
		std::string mesh_name = "references";
		std::vector<std::complex<double>> epsr_list = { std::complex<double>(6.0, -1.84) };
		std::vector<std::complex<double>> updated_qoi;
		int i = 0;
		//Domain dom_forward(mesh_name);
		//run::_standard(dom_forward, mesh_name, plot, false, higher_order, check_results, true);
		//HOPS::multi_HOPS_multi_epsr(mesh_name);
		
		//HOPS::multi_HOPS_epsr(mesh_name);
		Kriging::multi_HOPS_epsr(mesh_name);
		//HOPS::monte_carlo_dual(mesh_name);
		//Bayes::posterior(mesh_name);


		//run::_element_error(mesh_name, plot, useAdjoint, higher_order, check_results, .2);
		//run::_element_error_AF(mesh_name, plot, useAdjoint, higher_order, check_results, .2, 2.0/3.0);
		return 0;
}



