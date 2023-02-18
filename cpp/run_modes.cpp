#include "run_modes.h"
double max(std::vector<std::complex<double>>& vec) {
	double maximum = abs(vec[0]);
	for (int i = 1; i < vec.size(); ++i) {
		if (abs(vec[i]) > maximum) {
			maximum = abs(vec[i]);
		}
	}
	return maximum;
}
int minint(std::vector<int>& vec) {
	int min = vec[0];
	for (int i = 1; i < vec.size(); ++i) {
		if (vec[i] < min) min = vec[i];
	}
	return min;
}
int maxint(std::vector<int>& vec) {
	int max = vec[0];
	for (int i = 1; i < vec.size(); ++i) {
		if (vec[i] > max) max = vec[i];
	}
	return max;
}
void run::_standard(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, double freq)
{
	int extra = 0;
	if (higher_order) extra = 1;
	//dom = Domain(mesh_name);
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

	/////////////////////////////////Frequency Perturbation///////////////////////////////////////////////////
	//control frequency from input list instead of control file
	//dom.sc.fstart = freq;
	//dom.sc.fstop = freq;
	//--------------------------------------------------------------------------------------------------------

	/////////////////////////////////Material Perturbation///////////////////////////////////////////////////
	dom.sc.fstart = 70e6;
	dom.sc.fstop = 70e6;
	//--------------------------------------------------------------------------------------------------------

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
	if (!useAdjoint && get_RCS) {
		std::vector<double> rcs_list = postproc::post(dom);
		std::ofstream outFile("../exampleFiles/debug_data/" + mesh_name + "/results/rcs_list.txt");
		for (const auto& rcs : rcs_list) outFile << rcs << "\n";
		for (const auto& rcs : rcs_list) std::cout << rcs << "\n";
	}
	std::cout << "Num unknowns: " << dom.cAlpha.size() << std::endl;

}
void run::_standard_set(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, std::vector<int>& orders)
{
	int extra = 0;
	if (higher_order) extra = 1;
	//dom = Domain(mesh_name);
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
	//update the element orders
	int order_boost = 0;
	if (higher_order) order_boost = 1;
	refinement::set_refine(dom, orders, order_boost);
	//std::cout << "Adjusting NGL to be optimal!" << std::endl;
	//int maxorder = maxint(orders);
	//if (higher_order) maxorder += order_boost;
	//int nglset = ceil(double(maxorder + 1.0) / 2.0)+1;
	////if (nglset < 3) nglset = 3;
	//for (auto e = dom.elements.begin(); e != dom.elements.end(); ++e) {
	//	e->quadrature[0] = nglset;
	//	e->quadrature[1] = nglset;
	//	e->quadrature[2] = nglset;
	//	if (e->type == 1) {
	//		e->quadrature[0] = 9;
	//		e->quadrature[1] = 9;
	//		e->quadrature[2] = 9;
	//	}
	//}
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
	if (!useAdjoint && get_RCS) {
		std::vector<double> rcs_list = postproc::post(dom);
		std::ofstream outFile("../exampleFiles/debug_data/" + mesh_name + "/results/rcs_list.txt");
		for (const auto& rcs : rcs_list) outFile << rcs << "\n";
	}
	std::cout << "Num unknowns: " << dom.cAlpha.size() << std::endl;

}
void run::_RHS_compute(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results)
{
	int extra = 0;
	if (higher_order) extra = 1;
	//dom = Domain(mesh_name);
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
	dom._RHS_ONLY();

}

void run::convert_and_run(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, int ref_dex)
{
	int extra = 0;
	if (higher_order) extra = 1;
	//dom = Domain(mesh_name);
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

	//get reference files
	std::string command = "..\\file_input_converter ..\\reference_files_2im\\ref_" + std::to_string(ref_dex) + ".in ..\\exampleFiles\\" + mesh_name + "\\";
	//std::string command = "..\\file_input_converter ..\\reference_files\\ref_" + std::to_string(ref_dex) + ".in ..\\exampleFiles\\" + mesh_name + "\\";
	std::system(command.c_str());
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
	if (!useAdjoint && get_RCS) {
		std::vector<double> rcs_list = postproc::post(dom);
		std::ofstream outFile("../exampleFiles/debug_data/" + mesh_name + "/results/rcs_list.txt");
		for (const auto& rcs : rcs_list) outFile << rcs << "\n";
	}
	std::cout << "Num unknowns: " << dom.cAlpha.size() << std::endl;
}

void run::_element_error(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double toli) {
	//////////////////////////////////
	//vector of current (low-order) element expansion
	std::ofstream ref_file("refinement_iterations.txt");
	std::ofstream orders_file("refinement_orders.txt");
	std::ofstream outFile("refinement_rcs.txt");
	ref_file << std::setprecision(16);
	outFile << std::setprecision(16);
	std::vector<double> tols = { 0.05 };
	for (auto runs = 0; runs < tols.size(); ++runs) {
		double tol = tols[runs];
		ref_file << tol << std::endl;
		orders_file << tol << std::endl;
		std::vector<int> orders(256, 1);
		double errortol = 0.0;
		int iteration = 0;
		do {
			std::cout << "***********Iteration: " << iteration++ << std::endl;
			//get the forward solution
			Domain dom_forward_low(mesh_name);
			run::_standard_set(dom_forward_low, mesh_name, plot, false, false, check_results, false, orders);
			//dom_forward_low.clear_all_submatrices(dom_forward_low.elements);
			/*auto rcs_list = postproc::post(dom_forward_low);
			for (const auto &rcs : rcs_list) outFile << rcs << " ";
			outFile << std::endl;*/
			//print out cAlpha
			/*std::ofstream out_debug("debug_calpha.txt");
			for (int cI = 0; cI < dom_forward_low.cAlpha.size(); ++cI) {
				out_debug << dom_forward_low.cAlpha[cI] << std::endl;
			}
			out_debug.close();*/
			//find higher order adjoint solution
			Domain dom_adjoint_high(mesh_name);
			run::_standard_set(dom_adjoint_high, mesh_name, plot, true, true, check_results, false, orders);

			//find higher order RHS for forward problem
			Domain dom_high(mesh_name);
			//run::_RHS_compute(dom_forward_high, mesh_name, plot, false, true, check_results);
			auto adj_high_soln = dom_adjoint_high.cAlpha;
			dom_adjoint_high = Domain();//clear adjoint memory
			bool error = true;
			int extra = 0;
			if (higher_order) extra = 1;
			/////////////////////////////////
			std::complex<double> cSensitivity;

			dom_high.sc.useAdjoint = false;
			dom_high.higher_order = true;

			fileIOCon fileIO(&dom_high, higher_order, extra);
			dom_high.check_results = check_results;
			std::cout << "Running analysis now: " << std::endl;
			dom_high.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
			fileIO.file_input(mesh_name);
			dom_high.error = error;
			int order_boost = 1;
			refinement::set_refine(dom_high, orders, order_boost);
			////adjust NGL
			//std::cout << "Adjusting NGL to be optimal!" << std::endl;
			//int maxorder = maxint(orders);
			//if (dom_high.higher_order) maxorder += order_boost;
			//int nglset = ceil(double(maxorder + 1.0) / 2.0);
			//for (auto e = dom_high.elements.begin(); e != dom_high.elements.end(); ++e) {
			//	e->quadrature[0] = nglset;
			//	e->quadrature[1] = nglset;
			//	e->quadrature[2] = nglset;
			//	if (e->type == 1) {
			//		e->quadrature[0] = 9;
			//		e->quadrature[1] = 9;
			//		e->quadrature[2] = 9;
			//	}
			//}
			dom_high.makeNT();
			dom_high.connect_elements();
			dom_high.set_bc_elements();
			dom_high.fill_eiuvwijk();
			std::vector<int> vectorDUnique = dom_high.prenumunknowns();
			dom_high.make_unknown_description();
			dom_high.scatter1.waveSetup(dom_high.matDimCon, dom_high.sc.nfr, dom_high.sc.fstart, dom_high.sc.fstop, dom_high.sc.numberOfWaves);
			dom_high.scatter1.calcThetaPhi();
			dom_high.scatter1.calcFrequencies();
			dom_high.definers_defrmnls();

			dom_high.element_error_improved(dom_forward_low, dom_forward_low.cAlpha, adj_high_soln); //input the order of each element in forward solve, (WIP)
			std::vector<std::complex<double>> errors(dom_high.elements.size());
			for (int i = 0; i < dom_high.elements.size(); ++i) {
				errors[i] = dom_high.elements[i].element_error;
			}
			errortol = max(errors);
			if (errortol < tol) {
				std::cout << "Max element error: " << std::setprecision(16) << errortol << std::endl;
				std::complex<double> toterrorf = 0.0;
				for (auto el = dom_high.elements.begin(); el != dom_high.elements.end(); ++el) {
					toterrorf += el->element_error;
				}
				std::cout << "Tot error: " << abs(toterrorf) << std::endl;
				std::cout << "error below tol!" << std::endl;

				ref_file << errortol << " " << abs(toterrorf) << std::endl;
				for (int o = 0; o < orders.size(); ++o) {
					orders_file << orders[o] << " ";
				}
				orders_file << dom_forward_low.cAlpha.size() << std::endl;
				/*	orders_file.close();
					ref_file.close();*/
				std::vector < std::complex<double>> errorterm = { 0.0,0.0,toterrorf };
				auto rcs_list = postproc::postplus(dom_forward_low, errorterm);
				for (const auto& rcs : rcs_list) outFile << rcs << " ";
				outFile << std::endl;
				continue;
			}
			//figure out what the errors need to be set to..
			bool is_same = true;
			for (int i = 0; i < dom_high.elements.size(); ++i) {


				double h = dom_high.elements[i].h;

				double hi_scal = pow(tol / abs(dom_high.elements[i].element_error), 1.0 / orders[i]);
				std::cout << "Error: " << abs(dom_high.elements[i].element_error) << " h: " << h << " hscal: " << hi_scal << std::endl;


				//double p = orders[i] + (log(tol) - log(abs(dom_high.elements[i].element_error))) / abs(log(h));
				//if (tol - abs(dom_high.elements[i].element_error) < 0) p *= -1.0;
				//double p = log(tol / abs(dom_high.elements[i].element_error)) / log(hi_scal);
				//std::cout << "P double val: " << p << std::endl;
				//int pi = abs(round(p));

				/*std::cout << "Pi val: " << pi << std::endl;
				std::cout << "Hi/h val: " << hi_scal << std::endl;
				std::cout << "orig h val: " << h << std::endl;*/
				//if (pi == 0) pi = 1;
				//if (pi > 5) pi = 5;
				//if (orders[i] != abs(pi)) is_same = false;
				//double p = log(hi_scal) / log(tol / abs(dom_high.elements[i].element_error));
				double p = abs(log(tol / abs(dom_high.elements[i].element_error))) * h;

				int pi = ceil(p);
				std::cout << p << std::endl;
				if (hi_scal < 1.0) {
					orders[i] += pi;
					is_same = false;
					if (orders[i] > 6) {
						std::cout << "Exceeded expected!" << std::endl;
					}
				}
				/*if (hi_scal < 1.0) {
					++orders[i];
					is_same = false;
				}
				if (hi_scal < 0.5) ++orders[i];*/


				//if (hi_scal < 0.1) ++orders[i];
				//orders[i] = abs(pi);
				//errors[i] = dom_high.elements[i].element_error;

			}
			errortol = max(errors);
			std::cout << std::setprecision(16) << "Max element Error: " << errortol << std::endl;
			std::complex<double> toterror = 0.0;
			for (auto el = dom_high.elements.begin(); el != dom_high.elements.end(); ++el) {
				toterror += el->element_error;
			}
			std::cout << "Tot error: " << abs(toterror) << std::endl;
			std::cout << "Min order: " << minint(orders) << std::endl;
			std::cout << "Max order: " << maxint(orders) << std::endl;
			ref_file << errortol << " " << abs(toterror) << std::endl;
			/*if (minint(orders) == 4) {
				std::cout << "Reached max expansion order!" << std::endl;
				break;
			}*/
			for (int o = 0; o < orders.size(); ++o) {
				orders_file << orders[o] << " ";
			}
			orders_file << dom_forward_low.cAlpha.size() << std::endl;
			if (is_same) {
				//repeated orders
				std::cout << "Repeated order!" << std::endl;
				break;
			}
		} while (errortol > tol);
	}
	ref_file.close();
	orders_file.close();
	outFile.close();

}

void run::_element_error_AF(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double tot_error, double frac) {
	//////////////////////////////////
	//vector of current (low-order) element expansion
	std::ofstream ref_file("refinement_iterations.txt");
	std::ofstream orders_file("refinement_orders.txt");
	std::ofstream outFile("refinement_rcs.txt");
	ref_file << std::setprecision(16);
	outFile << std::setprecision(16);

	std::vector<double> tols = { 10.0 };
	for (auto runs = 0; runs < tols.size(); ++runs) {
		double tol = tols[runs];
		ref_file << tol << std::endl;
		orders_file << tol << std::endl;
		std::vector<int> orders(256, 1);
		//orders = {5,4,4,5,5,4,3,4,5,4,3,4,5,4,4,5,5,4,4,5,5,5,5,5,5,5,5,5,5,4,4,5,5,4,4,5,5,5,5,5,5,5,5,5,5,4,4,5,5,4,4,5,5,4,3,4,5,4,3,4,5,4,4,5,3,4,3,3,3,3,3,3,3,3,3,3,3,4,4,3,5,3,3,4,4,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,3,3,5,3,3,3,3,3,3,3,3,3,2,2,3,3,3,2,3,3,3,2,3,3,3,2,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,2,2,3,3,3,3,2,3,3,2,2,3,3,2,2,2,2,2,3,1,1,2,2,2,2,2,3,3,2,1,3,3,2,1,2,2,2,1};
		double errortol = 0.0;
		int iteration = 0;
		do {
			std::cout << "***********Iteration: " << iteration++ << std::endl;
			//get the forward solution
			Domain dom_forward_low(mesh_name);
			run::_standard_set(dom_forward_low, mesh_name, plot, false, false, check_results, false, orders);
			auto rcs_list = postproc::post(dom_forward_low);
			for (const auto& rcs : rcs_list) outFile << rcs << " ";
			outFile << std::endl;
			//print out cAlpha
			/*std::ofstream out_debug("debug_calpha.txt");
			for (int cI = 0; cI < dom_forward_low.cAlpha.size(); ++cI) {
			out_debug << dom_forward_low.cAlpha[cI] << std::endl;
			}
			out_debug.close();*/
			//find higher order adjoint solution
			Domain dom_adjoint_high(mesh_name);
			run::_standard_set(dom_adjoint_high, mesh_name, plot, true, true, check_results, false, orders);
			//find higher order RHS for forward problem
			Domain dom_high(mesh_name);
			//run::_RHS_compute(dom_forward_high, mesh_name, plot, false, true, check_results);
			auto adj_high_soln = dom_adjoint_high.cAlpha;
			dom_adjoint_high = Domain();//clear adjoint memory
			bool error = true;
			int extra = 0;
			if (higher_order) extra = 1;
			/////////////////////////////////
			std::complex<double> cSensitivity;

			dom_high.sc.useAdjoint = false;
			dom_high.higher_order = true;

			fileIOCon fileIO(&dom_high, higher_order, extra);
			dom_high.check_results = check_results;
			std::cout << "Running analysis now: " << std::endl;
			dom_high.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
			fileIO.file_input(mesh_name);
			dom_high.error = error;
			int order_boost = 1;
			refinement::set_refine(dom_high, orders, order_boost);

			dom_high.makeNT();
			dom_high.connect_elements();
			dom_high.set_bc_elements();
			dom_high.fill_eiuvwijk();
			std::vector<int> vectorDUnique = dom_high.prenumunknowns();
			dom_high.make_unknown_description();
			dom_high.scatter1.waveSetup(dom_high.matDimCon, dom_high.sc.nfr, dom_high.sc.fstart, dom_high.sc.fstop, dom_high.sc.numberOfWaves);
			dom_high.scatter1.calcThetaPhi();
			dom_high.scatter1.calcFrequencies();
			dom_high.definers_defrmnls();
			std::string connectivity = "../exampleFiles/debug_data/" + mesh_name + "/results/connectivity.txt";
			fileIO.output_connectivity(connectivity);
			dom_high.basis_error_improved(dom_forward_low, dom_forward_low.cAlpha, adj_high_soln);
			dom_high.element_error_improved(dom_forward_low, dom_forward_low.cAlpha, adj_high_soln); //input the order of each element in forward solve, (WIP)
			PlotController plotter(&dom_high);
			plotter.plotterIO("../exampleFiles/" + mesh_name + "/plotting.txt", "../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt", false);
			plotter.plotField(mesh_name, false);
			std::vector<std::complex<double>> errors(dom_high.elements.size());
			std::complex<double> abserror_sum = 0.0, error_sum;
			for (int i = 0; i < dom_high.elements.size(); ++i) {
				errors[i] = dom_high.elements[i].element_error;
				abserror_sum += abs(errors[i]);
				error_sum += errors[i];
			}
			errortol = abs(error_sum);
			double max_elem_error = max(errors);
			if (errortol < tol) {

				std::cout << "Tot error: " << abs(error_sum) << std::endl;
				std::cout << "error below tol!" << std::endl;

				ref_file << max_elem_error << " " << abs(error_sum) << std::endl;
				for (int o = 0; o < orders.size(); ++o) {
					orders_file << orders[o] << " ";
				}
				orders_file << dom_forward_low.cAlpha.size() << std::endl;
				/*	orders_file.close();
				ref_file.close();*/
				std::vector < std::complex<double>> errorterm = { 0.0,0.0,error_sum };
				auto rcs_list = postproc::postplus(dom_forward_low, errorterm);
				for (const auto& rcs : rcs_list) outFile << rcs << " ";
				outFile << std::endl;
				continue;
			}
			//figure out what the errors need to be set to..
			auto refiners = refinement::magnitude_refine_prop(errors, tol);
			for (int refi = 0; refi < refiners.size(); ++refi) {
				++orders[refiners[refi]];
			}
			std::cout << std::setprecision(16) << "Max element Error: " << max_elem_error << std::endl;
			std::complex<double> toterror = 0.0;
			for (auto el = dom_high.elements.begin(); el != dom_high.elements.end(); ++el) {
				toterror += el->element_error;
			}
			std::cout << "Tot error: " << abs(toterror) << std::endl;
			std::cout << "Min order: " << minint(orders) << std::endl;
			std::cout << "Max order: " << maxint(orders) << std::endl;
			ref_file << max_elem_error << " " << abs(toterror) << std::endl;
			/*if (minint(orders) == 4) {
			std::cout << "Reached max expansion order!" << std::endl;
			break;
			}*/
			for (int o = 0; o < orders.size(); ++o) {
				orders_file << orders[o] << " ";
			}
			orders_file << dom_forward_low.cAlpha.size() << std::endl;

		} while (errortol > tol);
	}
	ref_file.close();
	orders_file.close();
	outFile.close();

}

void run::_element_error_MSG(std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, double tot_error, double frac) {
	//////////////////////////////////
	//vector of current (low-order) element expansion
	std::ofstream ref_file("refinement_iterations.txt");
	std::ofstream orders_file("refinement_orders.txt");
	std::ofstream outFile("refinement_rcs.txt");
	ref_file << std::setprecision(16);
	outFile << std::setprecision(16);
	std::vector<int> orders(256, 1);
	std::vector<double> tols = { 10.0, 5.0, 1.0, 0.5 };
	for (auto runs = 0; runs < tols.size(); ++runs) {
		double tol = tols[runs];
		ref_file << tol << std::endl;
		orders_file << tol << std::endl;

		double errortol = 0.0;
		int iteration = 0;
		do {
			std::cout << "***********Iteration: " << iteration++ << std::endl;
			//get the forward solution
			Domain dom_forward_low(mesh_name);
			run::_standard_set(dom_forward_low, mesh_name, plot, false, false, check_results, false, orders);
			/*auto rcs_list = postproc::post(dom_forward_low);
			for (const auto &rcs : rcs_list) outFile << rcs << " ";
			outFile << std::endl;*/
			//print out cAlpha
			/*std::ofstream out_debug("debug_calpha.txt");
			for (int cI = 0; cI < dom_forward_low.cAlpha.size(); ++cI) {
			out_debug << dom_forward_low.cAlpha[cI] << std::endl;
			}
			out_debug.close();*/
			//find higher order adjoint solution
			Domain dom_adjoint_high(mesh_name);
			run::_standard_set(dom_adjoint_high, mesh_name, plot, true, true, check_results, false, orders);
			//find higher order RHS for forward problem
			Domain dom_high(mesh_name);
			//run::_RHS_compute(dom_forward_high, mesh_name, plot, false, true, check_results);
			auto adj_high_soln = dom_adjoint_high.cAlpha;
			dom_adjoint_high = Domain();//clear adjoint memory
			bool error = true;
			int extra = 0;
			if (higher_order) extra = 1;
			/////////////////////////////////
			std::complex<double> cSensitivity;

			dom_high.sc.useAdjoint = false;
			dom_high.higher_order = true;

			fileIOCon fileIO(&dom_high, higher_order, extra);
			dom_high.check_results = check_results;
			std::cout << "Running analysis now: " << std::endl;
			dom_high.check_results ? std::cout << "Checking results w/ FORTRAN!" << std::endl : std::cout << "Not checking results w/ FORTRAN!" << std::endl;
			fileIO.file_input(mesh_name);
			dom_high.error = error;
			int order_boost = 1;
			refinement::set_refine(dom_high, orders, order_boost);

			dom_high.makeNT();
			dom_high.connect_elements();
			dom_high.set_bc_elements();
			dom_high.fill_eiuvwijk();
			std::vector<int> vectorDUnique = dom_high.prenumunknowns();
			dom_high.make_unknown_description();
			dom_high.scatter1.waveSetup(dom_high.matDimCon, dom_high.sc.nfr, dom_high.sc.fstart, dom_high.sc.fstop, dom_high.sc.numberOfWaves);
			dom_high.scatter1.calcThetaPhi();
			dom_high.scatter1.calcFrequencies();
			dom_high.definers_defrmnls();

			dom_high.element_error_improved(dom_forward_low, dom_forward_low.cAlpha, adj_high_soln); //input the order of each element in forward solve, (WIP)
			std::vector<std::complex<double>> errors(dom_high.elements.size());
			std::complex<double> abserror_sum = 0.0, error_sum;
			for (int i = 0; i < dom_high.elements.size(); ++i) {
				errors[i] = dom_high.elements[i].element_error;
				abserror_sum += abs(errors[i]);
				error_sum += errors[i];
			}
			errortol = abs(error_sum);
			if (errortol < tol) {

				std::cout << "Tot error: " << abs(error_sum) << std::endl;
				std::cout << "error below tol!" << std::endl;

				ref_file << errortol << " " << abs(error_sum) << std::endl;
				for (int o = 0; o < orders.size(); ++o) {
					orders_file << orders[o] << " ";
				}
				orders_file << dom_forward_low.cAlpha.size() << std::endl;
				/*	orders_file.close();
				ref_file.close();*/
				std::vector < std::complex<double>> errorterm = { 0.0,0.0,error_sum };
				auto rcs_list = postproc::postplus(dom_forward_low, errorterm);
				for (const auto& rcs : rcs_list) outFile << rcs << " ";
				outFile << std::endl;
				continue;
			}
			//figure out what the errors need to be set to..
			auto refiners = refinement::magnitude_refine_prop(errors, frac);
			for (int refi = 0; refi < refiners.size(); ++refi) {
				++orders[refiners[refi]];
			}
			std::cout << std::setprecision(16) << "Max element Error: " << errortol << std::endl;
			std::complex<double> toterror = 0.0;
			for (auto el = dom_high.elements.begin(); el != dom_high.elements.end(); ++el) {
				toterror += el->element_error;
			}
			std::cout << "Tot error: " << abs(toterror) << std::endl;
			std::cout << "Min order: " << minint(orders) << std::endl;
			std::cout << "Max order: " << maxint(orders) << std::endl;
			ref_file << errortol << " " << abs(toterror) << std::endl;
			/*if (minint(orders) == 4) {
			std::cout << "Reached max expansion order!" << std::endl;
			break;
			}*/
			for (int o = 0; o < orders.size(); ++o) {
				orders_file << orders[o] << " ";
			}
			orders_file << dom_forward_low.cAlpha.size() << std::endl;

		} while (errortol > tol);
	}
	ref_file.close();
	orders_file.close();
	outFile.close();

}
void run::_standard_custom_materials(Domain& dom, std::string mesh_name, bool plot, bool useAdjoint, bool higher_order, bool check_results, bool get_RCS, std::vector<std::complex<double>>& epsr_values)
{
	int extra = 0;
	if (higher_order) extra = 1;
	//dom = Domain(mesh_name);
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
	//
	for (int el = 0; el < dom.elements.size(); ++el) {
		if (dom.elements[el].type == 0) {
			if (dom.elements[el].materials.region == 1) {
				dom.elements[el].materials.epsr_list[0][0] = epsr_values[el];
			}
		}
	}
	//
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
	if (!useAdjoint && get_RCS) {
		std::vector<double> rcs_list = postproc::post(dom);
		std::ofstream outFile("../exampleFiles/debug_data/" + mesh_name + "/results/rcs_list.txt");
		for (const auto& rcs : rcs_list) outFile << rcs << "\n";
	}
	std::cout << "Num unknowns: " << dom.cAlpha.size() << std::endl;

}