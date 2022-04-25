#include "refinement.h"

//Point operator*(const Point& p1, const double scalar) {
//	return Point(p1.x*scalar, p1.y*scalar, p1.z*scalar);
//}
//Point operator/(Point& p1, const double scalar) {
//	return Point(p1.x / scalar, p1.y / scalar, p1.z / scalar);
//}
//Point operator+(Point& p1, const Point p2) {
//	return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
//}
struct entry {
public:
	double val;
	int row;
	int col;
	entry(double vali, int rowi, int coli) {
		val = vali;
		row = rowi;
		col = coli;
	}
};


std::complex<double> refinement::add_qoi_error(std::string qoi_error_file, std::complex<double> qoi) {
	std::ifstream file;
	std::string line;
	std::complex<double> qoi_error;
	file.open(qoi_error_file);
	while (getline(file, line)) {
		auto line_input = functions::split(line, ' ');
		qoi_error += std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1]));
	}
	file.close();
	return qoi_error + qoi;
}

bool operator<(const entry& a, const entry&b) {
	if (a.val < b.val) return true;
	else return false;
}
std::vector<int> refinement::pick_elements_to_refine(bool p_refine, std::string qoi_errors_file, const int& n, std::vector<int>& coarsers, const int& p) {
	if (n == 0) return {};
	std::ifstream file;
	file.open(qoi_errors_file);
	std::string line;
	std::vector < std::complex<double>> qoi_errors;
	while (getline(file, line)) {
		auto line_input = functions::split(line, ' ');
		qoi_errors.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	//Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> k= Eigen::Map<Eigen::Matrix<std::complex<double>>, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
	Eigen::VectorXcd k = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
	std::vector<entry> cVec_vec;
	//Eigen::Matrix<entry, Eigen::Dynamic, 1> cVec;
	auto vert = k.replicate(1, qoi_errors.size());
	auto horz = vert.transpose();
	auto C = (vert + horz).cwiseAbs();

	for (int i = 0; i < C.rows(); ++i) {
		for (int j = 0; j < C.cols(); ++j) {
			if (i >= j) continue;

			cVec_vec.push_back(entry(C(i, j), i, j));
		}
	}
	std::sort(cVec_vec.rbegin(), cVec_vec.rend());
	std::vector<int> refine_elements;
	int num_pairs = 0;
	int pair_max = n;
	//refiners
	for (auto it = cVec_vec.begin(); it != cVec_vec.end(); ++it) {
		if (pair_max == 0) break;
		int e1 = it->col;
		int e2 = it->row;
		if (!(std::find(refine_elements.begin(), refine_elements.end(), e1) != refine_elements.end())
			&& !(std::find(refine_elements.begin(), refine_elements.end(), e2) != refine_elements.end())) {
			refine_elements.push_back(e1);
			refine_elements.push_back(e2);
			++num_pairs;
		}
		if (num_pairs >= pair_max) break;
	}
	//coarsers
	std::sort(cVec_vec.begin(), cVec_vec.end());
	int num_coarsers = 0;
	int coarse_pairs = p;
	for (auto it = cVec_vec.begin(); it != cVec_vec.end(); ++it) {
		if (coarse_pairs == 0) break;
		int e1 = it->col;
		int e2 = it->row;
		if (!(std::find(coarsers.begin(), coarsers.end(), e1) != coarsers.end())
			&& !(std::find(coarsers.begin(), coarsers.end(), e2) != coarsers.end())) {
			coarsers.push_back(e1);
			coarsers.push_back(e2);
			++num_coarsers;
		}
		if (num_coarsers >= coarse_pairs) break;
	}
	return refine_elements;
}

std::vector<int> refinement::magnitude_refine(bool p_refine, std::string qoi_errors_file, const int& n) {
	if (n == 0) return{};
	std::ifstream file;
	file.open(qoi_errors_file);
	std::string line;
	std::vector < std::complex<double>> qoi_errors;
	while (getline(file, line)) {
		auto line_input = functions::split(line, ' ');
		qoi_errors.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	std::vector<entry> cVec;
	int j = 0;
	for (int i = 0; i < qoi_errors.size(); ++i) {
		cVec.push_back(entry(abs(qoi_errors[i]), i, j));
	}
	std::sort(cVec.rbegin(), cVec.rend());
	std::vector<int> refine_elements;
	//int i = 0; //counter
	for (int i = 0; i < n; ++i){

		refine_elements.push_back(cVec[i].row);
	}
	return refine_elements;
}

std::vector<int> refinement::magnitude_refine_prop(std::vector<std::complex<double>>& qoi_errors, const double& tol) {

	std::vector<entry> cVec;
	double tot = 0.0;
	int j = 0;
	for (int i = 0; i < qoi_errors.size(); ++i) {
		cVec.push_back(entry(abs(qoi_errors[i]), i, j));
		tot += abs(qoi_errors[i]);
	}
	double frac = 1.0 - tol / tot;
	std::sort(cVec.rbegin(), cVec.rend());
	std::vector<int> refine_elements;
	//int i = 0; //counter
	double sum = 0.0;
	int i = 0;
	do {
		refine_elements.push_back(cVec[i].row);
		sum += abs(cVec[i++].val);
	} while (sum < frac*tot);
	return refine_elements;
}
std::vector<int> refinement::msg_refine_prop(std::vector<std::complex<double>>& qoi_errors, const double& frac) {
	
	//Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> k= Eigen::Map<Eigen::Matrix<std::complex<double>>, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
	Eigen::VectorXcd k = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
	std::vector<entry> cVec_vec;
	//Eigen::Matrix<entry, Eigen::Dynamic, 1> cVec;
	auto vert = k.replicate(1, qoi_errors.size());
	auto horz = vert.transpose();
	auto C = (vert + horz).cwiseAbs();

	for (int i = 0; i < C.rows(); ++i) {
		for (int j = 0; j < C.cols(); ++j) {
			if (i >= j) continue;

			cVec_vec.push_back(entry(C(i, j), i, j));
		}
	}
	std::sort(cVec_vec.rbegin(), cVec_vec.rend());
	std::vector<int> refine_elements;
	
	//refiners
	for (auto it = cVec_vec.begin(); it != cVec_vec.end(); ++it) {
		
		int e1 = it->col;
		int e2 = it->row;
		if (!(std::find(refine_elements.begin(), refine_elements.end(), e1) != refine_elements.end())
			&& !(std::find(refine_elements.begin(), refine_elements.end(), e2) != refine_elements.end())) {
			refine_elements.push_back(e1);
			refine_elements.push_back(e2);
		}
		
	}
	//coarsers
	
	return refine_elements;
}
void refinement::p_refine(std::vector<int> refine_elements, Domain& dom, bool neigh, const int& increment)
{
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		bool refine = true;
		/*if (dom.elements[*it].materials.pmlcode == 1) { 
			refine = false;
			*it = -1;
			continue;
		}*/
		/*for (auto f_it = dom.elements[*it].facet_indices.begin(); f_it != dom.elements[*it].facet_indices.end(); ++f_it) {
			if (dom.facets[*f_it - 1].type == -1) { 
				refine = false; 
				*it = -1;
				break;
			}
		}*/
		if (refine) {
			if (dom.elements[*it].refine == true) continue;
			else {
				dom.elements[*it].refine = true;
				dom.elements[*it].expansion[0] += increment;
				dom.elements[*it].expansion[1] += increment;
				dom.elements[*it].expansion[2] += increment;
			}
			if (neigh) {
				for (auto f_it = dom.elements[*it].facet_indices.begin(); f_it != dom.elements[*it].facet_indices.end(); ++f_it) {
					int neigh_elem = *it == dom.facets[*f_it - 1].element_indices.first - 1 ? dom.facets[*f_it - 1].element_indices.second - 1 : dom.facets[*f_it - 1].element_indices.first - 1;
					if (neigh_elem < 0) continue;
					if (dom.elements[neigh_elem].refine == true) continue;
					else {
						dom.elements[neigh_elem].refine = true;
						dom.elements[neigh_elem].expansion[0] += increment;
						dom.elements[neigh_elem].expansion[1] += increment;
						dom.elements[neigh_elem].expansion[2] += increment;
					}
				}
			}
		}
	}
}

void refinement::track_refine(std::vector<int> refine_elements, Domain& dom)
{
	for (int i = 0; i < dom.elements.size(); ++i) {
		Element e1 = dom.elements[i];
		int refine_sum;
		for (auto fac = e1.facet_indices.begin(); fac != e1.facet_indices.end(); ++fac) {
			Facet f1 = dom.facets[*fac - 1];
			Element e2 = dom.elements[f1.element_indices.second - 1];
			if (f1.element_indices.second == e1.index) e2 = dom.elements[f1.element_indices.first - 1];

		}
	}
}
Point operator+(const Point& p1, const Point& p2) {
	return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}
Point operator-(const Point& p1, const Point& p2) {
	return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}
Point operator*(const double scal, const Point& p1) {
	return Point(scal*p1.x, scal*p1.y, scal*p1.z);
}
Point find_halfway(Point p1, Point p2) {
	//returns the point inbetween 2 points
	Point p_ret;
	p_ret = p1 + (1.0 / 2.0)*(p2 - p1);
	return p_ret;
}
void findPowersLagrange(const int& kuvw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
	matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr, const int& ending) {
	//right now we have kuvw -- for now the order is the same in each direction -- this could be generalized later

	double x, xj, xm, f;
	for (int i = 0; i < ending; ++i) {

		x = xglu[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0*j - kuvw) / kuvw;
				if (j != m) {
					f = f*(x - xj) / (xm - xj);
				}
			}//for j
			fuPowersLagr(i, m) = f;
		}//for m



		x = xglv[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0*j - kuvw) / kuvw;
				if (j != m) {
					f = f*(x - xj) / (xm - xj);
				}
			}//for j
			fvPowersLagr(i, m) = f;
		}//for m



		x = xglw[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0*j - kuvw) / kuvw;
				if (j != m) {
					f = f*(x - xj) / (xm - xj);
				}
			}//for j
			fwPowersLagr(i, m) = f;
		}//for m
	}


}//void findPowersLagrange
double trilagrangeoduvw(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr) {
	return fuPowersLagr(m, i)*fvPowersLagr(n, j)*fwPowersLagr(l, k);
}
void find_eps_mu_matrix(const int& kuvw, const int& num_coords, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr,
	const int& size1, const std::vector<std::vector<std::complex<double>>>& muRel, matrix2d<std::complex<double>>& muRelInt) {
	//std::cout << "entered eps_mu_matrix" << std::endl;
	std::complex<double> cMu;
	for (auto s = 1; s <= size1; ++s) { //matrix entries
		for (int m = 0; m < num_coords; ++m) {
			//for (int n = 0; n <= nglv + 1; ++n) {
			//	for (int l = 0; l <= nglw + 1; ++l) {
			cMu = std::complex<double>(0, 0);
			for (int mat_loc = 0; mat_loc < muRel.size(); ++mat_loc) { //nodes
				int kuvwp = kuvw + 1;
				double f = trilagrangeoduvw(m, m, m, mat_loc % (kuvwp), int(floor(mat_loc / kuvwp)) % kuvwp,
					int(floor(mat_loc / (kuvwp*kuvwp))), fuPowersLagr, fvPowersLagr, fwPowersLagr);

				muRelInt(s, m) += muRel[mat_loc][s - 1] * f;
			}
			//	muRelInt(s, m, n, l) = cMu;
			//}
			//	}

		}

	}
	//std::cout << "exited eps_mu_match" << std::endl;
}
void refinement::h_refine_1st_order(std::vector<int> refine_elements, Domain& dom, double factor) {
	std::vector<Point> points;
	points = { Point(-1.0, -1.0, -1.0), 
		Point(1.0, -1.0, -1.0),
		Point(-1.0, 1.0, -1.0),
		Point(1.0, 1.0, -1.0) };
	for (int i = 0; i < 4; ++i) {
		points.push_back(points[i] + Point(0.0, 0.0, 2.0));
	}
	int size = points.size();
	for (int i = 0; i < size; ++i) {
			points.push_back(points[i] * factor);
	}
	indexer index_list;
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		Element e1 = dom.elements[*it - 1];
		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!

		if (e1.materials.hcode != 1 || e1.materials.icode != 1) continue;

		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!
		std::vector<Element> new_elements(6);

		int kuvw;
		int size1;
		if (e1.materials.icode == 0) {
			kuvw = e1.materials.KuvwA;
			size1 = 9 - 3 * e1.materials.sym;
		}
		else {
			kuvw = e1.materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> x_par, y_par, z_par;
		for (int par_it = 8; par_it < points.size(); ++par_it) {
			x_par.push_back(points[par_it].x);
			y_par.push_back(points[par_it].y);
			z_par.push_back(points[par_it].z);
		}
		matrix2d < std::complex<double>> epsrNew = matrix2d<std::complex<double>>(size1, x_par.size());
		matrix2d < std::complex<double>> murNew = matrix2d<std::complex<double>>(size1, x_par.size());
		if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
			epsrNew(1, 1) = e1.materials.epsr_list[0][0];
			murNew(1, 1) = e1.materials.mur_list[0][0];
		}
		else {
			//material params at one face already known, get material params at other locations
			matrix2d<double> fuPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fvPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fwPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);

			findPowersLagrange(kuvw, x_par, y_par, z_par, fuPowersLagr, fvPowersLagr, fwPowersLagr, x_par.size());
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.epsr_list, epsrNew);
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.mur_list, murNew);
		}
		//add the points to the original set of points
		std::vector<double> r(4);
		//std::vector<Point> real_domain;
		std::vector<int> real_node_indices(points.size());
		for (auto p = 0; p < points.size(); ++p) {
			if (p < 8) real_node_indices[p] = e1.all_indices[p] - 1;
			//else if (p == 40) real_node_indices[p] = e1.all_indices[13]-1;
			else {
				int temp_index = dom.nodes.size();
				r = unitVectorsM::findR(e1.geom_order, points[p].x, points[p].y, points[p].z, e1.rs, e1.nRs);
				dom.nodes.push_back(Point(r[1], r[2], r[3]));
				real_node_indices[p] = temp_index;
			}
		}
		for (int e_new = 0; e_new < 7; ++e_new) {
			Element e_test;
			std::vector<int> order_pt_indices;
			std::vector<std::vector<std::complex<double>>> new_epsr;
			std::vector<std::vector<std::complex<double>>> new_mur;
			bool do_extra_mats = false;
			if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
				new_epsr = e1.materials.epsr_list;
				new_mur = e1.materials.mur_list;
			}
			else {
				do_extra_mats = true;
				new_epsr = std::vector<std::vector<std::complex<double>>>(e1.materials.epsr_list.size(), std::vector<std::complex<double>>(e1.materials.epsr_list[0].size()));
				new_mur = new_epsr;
			}
			std::vector<int> indices = index_list.first_indices[e_new];
			for (int i = 0; i < indices.size(); ++i) {
				order_pt_indices.push_back(real_node_indices[indices[i] - 1] + 1);
			}


			//now have vertices matched up properly
			//turn materials matrix into the proper form

			e_test.all_indices = order_pt_indices;
			e_test.expansion = e1.expansion;
			e_test.quadrature = e1.quadrature;
			e_test.materials.mur_list = new_mur;
			e_test.materials.epsr_list = new_epsr;
			e_test.index = dom.elements.size() + 1;
			e_test.materials.hcode = e1.materials.hcode;
			e_test.materials.icode = e1.materials.icode;
			e_test.materials.pmlcode = e1.materials.pmlcode;
			e_test.materials.Kuvw = e1.materials.Kuvw;
			e_test.materials.sym = e1.materials.sym;
			e_test.geom_order = e1.geom_order;
			if (e_new < 6) dom.elements.push_back(e_test);
			else {
				e_test.index = *it;
				dom.elements[*it - 1] = e_test;
			}
		}
	}
	//print out all necessary information to then redo the mesh
	std::ofstream file_elem, file_geom, file_bc, file_basic, file_material;
	file_geom.open("../Geometry_refined.dat");
	file_geom << "!Refined geometry!" << std::endl;
	file_elem.open("../Element_refined.dat");
	file_elem << "!Refined elements!" << std::endl;
	file_bc.open("../Boundary_refined.dat");
	file_bc << "!Refined boundary!" << std::endl;
	file_basic.open("../Basic_refined.dat");
	file_basic << "!Basic info file" << std::endl;
	file_material.open("../exampleFiles/" + dom.mesh_name + "_refine/materials_refined.txt");
	file_material << "!Material file" << std::endl;


	for (auto dom_points = dom.nodes.begin(); dom_points != dom.nodes.end(); ++dom_points) {
		file_geom << dom_points->x << " " << dom_points->y << " " << dom_points->z << std::endl;
	}
	file_geom.close();
	for (auto elem_p = dom.elements.begin(); elem_p != dom.elements.end(); ++elem_p) {
		file_elem << "Index: " << elem_p->index << " " << elem_p->materials.pmlcode << " " << elem_p->geom_order;
		file_elem << std::endl << "Basis: " << elem_p->expansion[0] << " " << elem_p->expansion[1] << " " << elem_p->expansion[2];
		file_elem << std::endl << "Quadrature: " << elem_p->quadrature[0] << " " << elem_p->quadrature[1] << " " << elem_p->quadrature[2] << std::endl;
		file_elem << "Nodes: ";
		for (auto elm_nodes = elem_p->all_indices.begin(); elm_nodes != elem_p->all_indices.end(); ++elm_nodes) {
			file_elem << *elm_nodes << " ";
		}


		file_elem << std::endl;
		file_material << elem_p->index << " " << elem_p->materials.hcode << " " << elem_p->materials.icode << " " << elem_p->materials.pmlcode << std::endl;
		file_material << elem_p->materials.Kuvw << " " << elem_p->materials.sym << std::endl;
		file_material << "eps:";// << std::endl;

		for (auto mat_ent = elem_p->materials.epsr_list.begin(); mat_ent != elem_p->materials.epsr_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl << "mu:";
		for (auto mat_ent = elem_p->materials.mur_list.begin(); mat_ent != elem_p->materials.mur_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl;
	}
	file_elem.close();
	file_material.close();
	//save boundary condition info
	int bc_index = 1;
	for (auto fac_bc = dom.facets.begin(); fac_bc != dom.facets.end(); ++fac_bc) {
		if (fac_bc->boundary_condition == 0) continue;
		file_bc << bc_index++ << " " << fac_bc->vertices[0] << " " << fac_bc->vertices[1] << " " << fac_bc->vertices[2] << " " << fac_bc->vertices[3] << " " << fac_bc->boundary_condition << std::endl;
	}
	file_bc.close();

	//MODE - NODES - NOEL - NBC - NGEN - NWAVES - NPORTS - FSTART - FSTOP - NFR
	file_basic << dom.sc.mode_of_operation << " " << dom.nodes.size() << " " << dom.elements.size() << " " << dom.sc.nbc << " " << dom.sc.ngen << " " << dom.sc.numberOfWaves << " " << dom.sc.nports << " " << dom.sc.fstart << " " << dom.sc.fstop << " " << dom.sc.nfr << std::endl;
	file_basic.close();
}

void refinement::ga_refine(Domain & dom, std::vector<int>& refine_elements)
{
	for (int e = 0; e < dom.elements.size(); ++e) {
		int order = refine_elements[e];
		dom.elements[e].expansion[0] = order;
		dom.elements[e].expansion[1] = order;
		dom.elements[e].expansion[2] = order;
	}
	std::cout << "Set elements using genetic algorithm!" << std::endl;
}

void refinement::set_refine(Domain & dom, std::vector<int>& refine_elements, int higher)
{
	for (int e = 0; e < dom.elements.size(); ++e) {
		int order = refine_elements[e] + higher;
		dom.elements[e].expansion[0] = order;
		dom.elements[e].expansion[1] = order;
		dom.elements[e].expansion[2] = order;
	}
	std::cout << "Set elements using genetic algorithm!" << std::endl;
}

void refinement::h_refineFin(std::vector<int> refine_elements, Domain& dom, double factor) {
	std::vector<Point> points;
	points = {Point(-1.0, -1.0, -1.0), Point(0.0, -1.0, -1.0), 
		Point(1.0, -1.0, -1.0), Point(-1.0, 0.0, -1.0),Point(0.0, 0.0, -1.0),
	Point(1.0, 0.0, -1.0) ,Point(-1.0, 1.0, -1.0) ,Point(0.0, 1.0, -1.0), 
		Point(1.0, 1.0, -1.0) };
	for (int i = 0; i < 9; ++i) {
		points.push_back(points[i] + Point(0.0, 0.0, 1.0));
	}
	for (int i = 0; i < 9; ++i) {
		points.push_back(points[i] + Point(0.0, 0.0, 2.0));
	}
	int size = points.size();
	for (int i = 0; i < size; ++i) {
		if (i != 13)
		points.push_back(points[i] * factor);
	}
	int dif = 27;
	for (int i = 0; i < 27; ++i) {
		if (i != 13)
		points.push_back(find_halfway(points[i], points[i + dif]));
		else dif -= 1;
	}
	indexer index_list;
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		Element e1 = dom.elements[*it - 1];
		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!

		if (e1.materials.hcode != 1 || e1.materials.icode != 1) continue;

		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!
		std::vector<Element> new_elements(6);

		int kuvw;
		int size1;
		if (e1.materials.icode == 0) {
			kuvw = e1.materials.KuvwA;
			size1 = 9 - 3 * e1.materials.sym;
		}
		else {
			kuvw = e1.materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> x_par, y_par, z_par;
		for (int par_it = 27; par_it < points.size(); ++par_it) {
			x_par.push_back(points[par_it].x);
			y_par.push_back(points[par_it].y);
			z_par.push_back(points[par_it].z);
		}
		matrix2d < std::complex<double>> epsrNew = matrix2d<std::complex<double>>(size1, x_par.size());
		matrix2d < std::complex<double>> murNew = matrix2d<std::complex<double>>(size1, x_par.size());
		if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
			epsrNew(1, 1) = e1.materials.epsr_list[0][0];
			murNew(1, 1) = e1.materials.mur_list[0][0];
		}
		else {
			//material params at one face already known, get material params at other locations
			matrix2d<double> fuPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fvPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fwPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);

			findPowersLagrange(kuvw, x_par, y_par, z_par, fuPowersLagr, fvPowersLagr, fwPowersLagr, x_par.size());
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.epsr_list, epsrNew);
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.mur_list, murNew);
		}
		//add the points to the original set of points
		std::vector<double> r(4);
		//std::vector<Point> real_domain;
		std::vector<int> real_node_indices(points.size());
		for (auto p = 0; p < points.size(); ++p) {
			if (p < 27) real_node_indices[p] = e1.all_indices[p] - 1;
			//else if (p == 40) real_node_indices[p] = e1.all_indices[13]-1;
			else {
				int temp_index = dom.nodes.size();
				r = unitVectorsM::findR(e1.geom_order, points[p].x, points[p].y, points[p].z, e1.rs, e1.nRs);
				dom.nodes.push_back(Point(r[1], r[2], r[3]));
				real_node_indices[p] = temp_index;
			}
		}
		//dom.nodes.insert(std::end(dom.nodes), std::begin(real_domain), std::end(real_domain));
		for (int e_new = 0; e_new < 7; ++e_new) {
			Element e_test;
			std::vector<int> order_pt_indices;
			std::vector<std::vector<std::complex<double>>> new_epsr;
			std::vector<std::vector<std::complex<double>>> new_mur;
			bool do_extra_mats = false;
			if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
				new_epsr = e1.materials.epsr_list;
				new_mur = e1.materials.mur_list;
			}
			else {
				do_extra_mats = true;
				new_epsr = std::vector<std::vector<std::complex<double>>>(e1.materials.epsr_list.size(), std::vector<std::complex<double>>(e1.materials.epsr_list[0].size()));
				new_mur = new_epsr;
			}
			std::vector<int> indices = index_list.indices[e_new];
			for (int i = 0; i < indices.size(); ++i) {
				order_pt_indices.push_back(real_node_indices[indices[i]-1] + 1);
			}
			
		
			//now have vertices matched up properly
			//turn materials matrix into the proper form

			e_test.all_indices = order_pt_indices;
			e_test.expansion = e1.expansion;
			e_test.quadrature = e1.quadrature;
			e_test.materials.mur_list = new_mur;
			e_test.materials.epsr_list = new_epsr;
			e_test.index = dom.elements.size() + 1;
			e_test.materials.hcode = e1.materials.hcode;
			e_test.materials.icode = e1.materials.icode;
			e_test.materials.pmlcode = e1.materials.pmlcode;
			e_test.materials.Kuvw = e1.materials.Kuvw;
			e_test.materials.sym = e1.materials.sym;
			e_test.geom_order = e1.geom_order;
			if (e_new < 6) dom.elements.push_back(e_test);
			else { 
				e_test.index = *it;
				dom.elements[*it - 1] = e_test; 
			}
		}
	}
	//print out all necessary information to then redo the mesh
	std::ofstream file_elem, file_geom, file_bc, file_basic, file_material;
	file_geom.open("../Geometry_refined.dat");
	file_geom << "!Refined geometry!" << std::endl;
	file_elem.open("../Element_refined.dat");
	file_elem << "!Refined elements!" << std::endl;
	file_bc.open("../Boundary_refined.dat");
	file_bc << "!Refined boundary!" << std::endl;
	file_basic.open("../Basic_refined.dat");
	file_basic << "!Basic info file" << std::endl;
	file_material.open("../exampleFiles/" + dom.mesh_name + "_refine/materials_refined.txt");
	file_material << "!Material file" << std::endl;


	for (auto dom_points = dom.nodes.begin(); dom_points != dom.nodes.end(); ++dom_points) {
		file_geom << dom_points->x << " " << dom_points->y << " " << dom_points->z << std::endl;
	}
	file_geom.close();
	for (auto elem_p = dom.elements.begin(); elem_p != dom.elements.end(); ++elem_p) {
		file_elem << "Index: " << elem_p->index << " " << elem_p->materials.pmlcode << " " << elem_p->geom_order;
		file_elem << std::endl << "Basis: " << elem_p->expansion[0] << " " << elem_p->expansion[1] << " " << elem_p->expansion[2];
		file_elem << std::endl << "Quadrature: " << elem_p->quadrature[0] << " " << elem_p->quadrature[1] << " " << elem_p->quadrature[2] << std::endl;
		file_elem << "Nodes: ";
		
		for (auto elm_nodes = elem_p->all_indices.begin(); elm_nodes != elem_p->all_indices.end(); ++elm_nodes) {

			file_elem << *elm_nodes << " ";
		}


		file_elem << std::endl;
		file_material << elem_p->index << " " << elem_p->materials.hcode << " " << elem_p->materials.icode << " " << elem_p->materials.pmlcode << std::endl;
		file_material << elem_p->materials.Kuvw << " " << elem_p->materials.sym << std::endl;
		file_material << "eps:";// << std::endl;

		for (auto mat_ent = elem_p->materials.epsr_list.begin(); mat_ent != elem_p->materials.epsr_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl << "mu:";
		for (auto mat_ent = elem_p->materials.mur_list.begin(); mat_ent != elem_p->materials.mur_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl;
	}
	file_elem.close();
	file_material.close();
	//save boundary condition info
	int bc_index = 1;
	for (auto fac_bc = dom.facets.begin(); fac_bc != dom.facets.end(); ++fac_bc) {
		if (fac_bc->boundary_condition == 0) continue;
		file_bc << bc_index++ << " " << fac_bc->vertices[0] << " " << fac_bc->vertices[1] << " " << fac_bc->vertices[2] << " " << fac_bc->vertices[3] << " " << fac_bc->boundary_condition << std::endl;
	}
	file_bc.close();

	//MODE - NODES - NOEL - NBC - NGEN - NWAVES - NPORTS - FSTART - FSTOP - NFR
	file_basic << dom.sc.mode_of_operation << " " << dom.nodes.size() << " " << dom.elements.size() << " " << dom.sc.nbc << " " << dom.sc.ngen << " " << dom.sc.numberOfWaves << " " << dom.sc.nports << " " << dom.sc.fstart << " " << dom.sc.fstop << " " << dom.sc.nfr << std::endl;
	file_basic.close();
}
void refinement::h_refine2(std::vector<int> refine_elements, Domain& dom, double factor) {
	//std::vector<Element> new_elements;
	std::vector<Point> p_central;
	p_central.resize(27);
	//build central cube CORNER VERTICES
	p_central[0] = (Point(-1.0*factor, -1.0*factor, -1.0*factor));
	p_central[2] = (Point(1.0*factor, -1.0*factor, -1.0*factor));
	p_central[6] = (Point(-1.0*factor, 1.0*factor, -1.0*factor));
	p_central[8] = (Point(1.0*factor, 1.0*factor, -1.0*factor));
	p_central[18] = (Point(-1.0*factor, -1.0*factor, 1.0*factor));
	p_central[20] = (Point(1.0*factor, -1.0*factor, 1.0*factor));
	p_central[24] = (Point(-1.0*factor, 1.0*factor, 1.0*factor));
	p_central[26] = (Point(1.0*factor, 1.0*factor, 1.0*factor));
	//build just the vertices on faces
	p_central[1] = find_halfway(p_central[0], p_central[2]);
	p_central[3] = find_halfway(p_central[0], p_central[6]);
	p_central[4] = find_halfway(p_central[0], p_central[8]);
	p_central[5] = find_halfway(p_central[2], p_central[8]);
	p_central[7] = find_halfway(p_central[6], p_central[8]);
	p_central[9] = find_halfway(p_central[0], p_central[18]);
	p_central[10] = find_halfway(p_central[0], p_central[20]);
	p_central[11] = find_halfway(p_central[2], p_central[20]);
	p_central[12] = find_halfway(p_central[0], p_central[24]);
	p_central[14] = find_halfway(p_central[2], p_central[26]);
	p_central[13] = find_halfway(p_central[12], p_central[14]);
	p_central[15] = find_halfway(p_central[6], p_central[24]);
	p_central[16] = find_halfway(p_central[6], p_central[26]);
	p_central[17] = find_halfway(p_central[8], p_central[26]);
	p_central[19] = find_halfway(p_central[18], p_central[20]);
	p_central[21] = find_halfway(p_central[18], p_central[24]);
	p_central[22] = find_halfway(p_central[18], p_central[26]);
	p_central[23] = find_halfway(p_central[20], p_central[26]);
	p_central[25] = find_halfway(p_central[24], p_central[26]);
	
	//need to make the intermediate elements for the 6 additional elems
	std::vector<Point> interm_points(27);
	for (int lines = 0; lines < p_central.size(); ++lines) {
		if (lines == 13) continue; //this node is inside
		interm_points[lines] = (find_halfway(p_central[lines], (1.0/factor)*p_central[lines]));
	}
	std::vector<std::vector<int>> elem_pts(6); //stored per face index
	for (int i = 0; i < 6; ++i) {
		//std::vector<int> elem_cent_pts, elem_interm_pts, elem_interp_pts;

		 //account for pts on the central element, the interp points from the orig eleme and the new interm pts

		//elem zero gets 0,3,6,9,12,15,18,21,24
		//elem one gets 2,5,8,11,14,17,20,23,26
		/*
		elem0 and 1 gets lowest + 3 each time for 9 vertices from central

		elem 2 and 3 gets lowest, lowest + 1, lowest+2, then increment by 9 and repeat twice

		elem 4 and 5 gets lowest up to lowest + 8
		*/
		switch (i) {
		case 0:
		case 1:
		{
			int start1 = i * 2;
			for (int k = 0; k < 9; ++k) {
				elem_pts[i].push_back(start1 + 3 * k);
			}
			break;
		}
		case 2:
		case 3:
		{
			int start2 = (i - 2) * 6;
			int adjust = 0;
			for (int k = 0; k < 3; ++k) {
				adjust = k * 9;
				elem_pts[i].push_back(start2 + adjust);
				elem_pts[i].push_back(start2 + adjust + 1);
				elem_pts[i].push_back(start2 + adjust + 2);
			}
			break;
		}
		case 4:
		case 5:
		{
			int start3 = (i - 4) * 18;
			for (int k = 0; k < 9; ++k) {
				elem_pts[i].push_back(start3 + k);
			}
			break;
		}
		}
	}
	//iterate through elements now
	//std::vector<Element> overall_new_elements;
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		//iterate through elements to refine
		//in an element, take a smaller sized cube in the parametric domain in its center
		//map the vertices of these coordinates to the real domain
		//connect the vertices of each face to each face, forming new elements
		//preserve the interp points!
		Element e1 = dom.elements[*it - 1];
		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!

		if (e1.materials.hcode != 1 || e1.materials.icode != 1) continue;

		//temp condition b/c material isn't debugged yet!!!!!!!!!!!!!!
		std::vector<Element> new_elements(6);

		int kuvw;
		int size1;
		if (e1.materials.icode == 0) {
			kuvw = e1.materials.KuvwA;
			size1 = 9 - 3 * e1.materials.sym;
		}
		else {
			kuvw = e1.materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> x_par, y_par, z_par;
		for (int par_it = 0; par_it < p_central.size(); ++par_it) {
			x_par.push_back(p_central[par_it].x);
			y_par.push_back(p_central[par_it].y);
			z_par.push_back(p_central[par_it].z);
		}
		for (int par_it = 0; par_it < interm_points.size(); ++par_it) {
			x_par.push_back(interm_points[par_it].x);
			y_par.push_back(interm_points[par_it].y);
			z_par.push_back(interm_points[par_it].z);
		}
		matrix2d < std::complex<double>> epsrNew = matrix2d<std::complex<double>>(size1, x_par.size());
		matrix2d < std::complex<double>> murNew = matrix2d<std::complex<double>>(size1, x_par.size());
		if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
			epsrNew(1,1) = e1.materials.epsr_list[0][0];
			murNew(1,1) = e1.materials.mur_list[0][0];
		}
		else {
			//material params at one face already known, get material params at other locations
			matrix2d<double> fuPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fvPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			matrix2d<double> fwPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
			
			findPowersLagrange(kuvw, x_par, y_par, z_par, fuPowersLagr, fvPowersLagr, fwPowersLagr, x_par.size());
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.epsr_list, epsrNew);
			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.mur_list, murNew);
		}
		//add the points to the original set of points
		std::vector<double> r(4);
		std::vector<Point> real_domain;
		for (auto p = p_central.begin(); p != p_central.end(); ++p) {
			if (p - p_central.begin() == 13) continue;
			
			r = unitVectorsM::findR(e1.geom_order, p->x, p->y, p->z, e1.rs, e1.nRs);
			real_domain.push_back(Point(r[1], r[2], r[3]));
		}
		for (auto p = interm_points.begin(); p != interm_points.end(); ++p) {
			if (p - interm_points.begin() == 13) continue;
			r = unitVectorsM::findR(e1.geom_order, p->x, p->y, p->z, e1.rs, e1.nRs);
			real_domain.push_back(Point(r[1], r[2], r[3]));
		}
		int cent_start = dom.nodes.size()+1;
		int interm_start = dom.nodes.size() + 27;
		dom.nodes.insert(std::end(dom.nodes), std::begin(real_domain), std::end(real_domain));
		//add the very central point
		auto p = p_central[13];
		r = unitVectorsM::findR(e1.geom_order, p.x, p.y, p.z, e1.rs, e1.nRs);
		dom.nodes.push_back(Point(r[1], r[2], r[3]));
		for (int e_new = 0; e_new < 6; ++e_new) {
			Element e_test;
			std::vector<int> order_pt_indices;
			std::vector<std::vector<std::complex<double>>> new_epsr;
			std::vector<std::vector<std::complex<double>>> new_mur;
			bool do_extra_mats = false;
			if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
				new_epsr = e1.materials.epsr_list;
				new_mur = e1.materials.mur_list;
			}
			else {
				do_extra_mats = true;
				new_epsr = std::vector<std::vector<std::complex<double>>>(e1.materials.epsr_list.size(), std::vector<std::complex<double>>(e1.materials.epsr_list[0].size()));
				new_mur = new_epsr;
			}
			

			if (e_new % 2 == 0) {

				//is even
				//order starts at orig
				int loc = 0;
				switch (e_new) {
				case 0:
				{
					//orig, interm, central
					//all x3

					for (int order = 0; order < 9; ++order) {
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order]]);
						order_pt_indices.push_back(elem_pts[e_new][order] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order] + cent_start);
						if (do_extra_mats) {
							std::vector<std::complex<double>> temp_vec = (e1.materials.epsr_list[elem_pts[e_new][order]]);
							//temp_vec.push_back;
							new_epsr[order * 3] = temp_vec;
							temp_vec = e1.materials.mur_list[elem_pts[e_new][order]];
							new_mur[order * 3] = temp_vec;
						}
					}
					break;
				}
				case 2:
				{
					//orig orig orig
					//interm interm interm
					//central central central
					//all x 3
					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 1]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 2]]);
						//
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
						//
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
					}
					break;
				}
				case 4:
				{
					//orig orig orig x 3
					//interm interm interm x 3
					//central central central x 3
					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 1]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 2]]);
					}
					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
					}
					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
					}


					break;
				}
					
				}
			}
			else {
				//is odd
				//order starts at central
				switch (e_new) {
				case 1:
				{
					for (int order = 0; order < 9; ++order) {

						order_pt_indices.push_back(elem_pts[e_new][order] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order] + interm_start);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order]]);
					}
					break;
				}
				case 3:
				{
					for (int order = 0; order < 3; ++order) {

						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
						//

						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
						//
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 1]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 2]]);


					}
					break;
				}
				case 5:
				{
					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
					}

					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
					}

					for (int order = 0; order < 3; ++order) {
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 1]]);
						order_pt_indices.push_back(e1.all_indices[elem_pts[e_new][order * 3 + 2]]);
					}

					break;
				}
				}
			}
			//now have vertices matched up properly
			//turn materials matrix into the proper form
			
			e_test.all_indices = order_pt_indices;
			e_test.expansion = e1.expansion;
			e_test.quadrature = e1.quadrature;
			e_test.materials.mur_list = new_mur;
			e_test.materials.epsr_list = new_epsr;
			e_test.index = dom.elements.size()+1;
			e_test.materials.hcode = e1.materials.hcode;
			e_test.materials.icode = e1.materials.icode;
			e_test.materials.pmlcode = e1.materials.pmlcode;
			e_test.materials.Kuvw = e1.materials.Kuvw;
			e_test.materials.sym = e1.materials.sym;
			e_test.geom_order = e1.geom_order;
			dom.elements.push_back(e_test);
		}

		//do the inner element last
		Element inner_element;
		std::vector<int> inner_indices;
		for (int order = 0; order < 27; ++order) {
			if (order != 13)
				inner_indices.push_back(cent_start + order);
			else inner_indices.push_back(dom.nodes.size() + 1);
		}
		inner_element.all_indices = inner_indices;
		inner_element.expansion = e1.expansion;
		inner_element.quadrature = e1.quadrature;
		////////////////////////////////////////////////fix this portion!!!!!!!!!!!!!! need to adjust not constant mat params
		//////////////////////////////////////////////////////////////////
		inner_element.materials.mur_list = e1.materials.mur_list;
		inner_element.materials.epsr_list = e1.materials.epsr_list;
		///////////////////////////////////////////////////////////////////
		inner_element.index = e1.index;
		inner_element.materials.hcode = e1.materials.hcode;
		inner_element.materials.icode = e1.materials.icode;
		inner_element.materials.pmlcode = e1.materials.pmlcode;
		inner_element.materials.Kuvw = e1.materials.Kuvw;
		inner_element.materials.sym = e1.materials.sym;
		inner_element.geom_order = e1.geom_order;
		//kill spawning element
		dom.elements[*it - 1] = inner_element;
	}
	//print out all necessary information to then redo the mesh
	std::ofstream file_elem, file_geom, file_bc, file_basic, file_material;
	file_geom.open("../Geometry_refined.dat");
	file_geom << "!Refined geometry!" << std::endl;
	file_elem.open("../Element_refined.dat");
	file_elem << "!Refined elements!" << std::endl;
	file_bc.open("../Boundary_refined.dat");
	file_bc << "!Refined boundary!" << std::endl;
	file_basic.open("../Basic_refined.dat");
	file_basic << "!Basic info file" << std::endl;
	file_material.open("../exampleFiles/" + dom.mesh_name + "_refine/materials_refined.txt");
	file_material << "!Material file" << std::endl;
	

	for (auto dom_points = dom.nodes.begin(); dom_points != dom.nodes.end(); ++dom_points) {
		file_geom << dom_points->x << " " << dom_points->y << " " << dom_points->z << std::endl;
	}
	file_geom.close();
	for (auto elem_p = dom.elements.begin(); elem_p != dom.elements.end(); ++elem_p) {
		file_elem << "Index: " << elem_p->index << " " << elem_p->materials.pmlcode << " " << elem_p->geom_order;
		file_elem << std::endl << "Basis: " << elem_p->expansion[0] << " " << elem_p->expansion[1] << " " << elem_p->expansion[2];
		file_elem << std::endl << "Quadrature: " << elem_p->quadrature[0] << " " << elem_p->quadrature[1] << " " << elem_p->quadrature[2] << std::endl;
		file_elem << "Nodes: ";
		
		for (auto elm_nodes = elem_p->all_indices.begin(); elm_nodes != elem_p->all_indices.end(); ++elm_nodes) { 
			file_elem << *elm_nodes << " ";
		}
		
		
		file_elem << std::endl;
		file_material << elem_p->index << " " << elem_p->materials.hcode << " " << elem_p->materials.icode << " " << elem_p->materials.pmlcode << std::endl;
			file_material << elem_p->materials.Kuvw << " " << elem_p->materials.sym << std::endl;
		file_material << "eps:";// << std::endl;
		
		for (auto mat_ent = elem_p->materials.epsr_list.begin(); mat_ent != elem_p->materials.epsr_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl << "mu:";
		for (auto mat_ent = elem_p->materials.mur_list.begin(); mat_ent != elem_p->materials.mur_list.end(); ++mat_ent) {
			file_material << std::endl;
			//iterates through nodes
			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
				//gets values
				file_material << " " << val_mat->real() << "," << val_mat->imag();
			}
		}
		file_material << std::endl;
	}
	file_elem.close();
	file_material.close();
	//save boundary condition info
	int bc_index = 1;
	for (auto fac_bc = dom.facets.begin(); fac_bc != dom.facets.end(); ++fac_bc) {
		if (fac_bc->boundary_condition == 0) continue;
		file_bc << bc_index++ << " " << fac_bc->vertices[0] << " " << fac_bc->vertices[1] << " " << fac_bc->vertices[2] << " " << fac_bc->vertices[3] << " " << fac_bc->boundary_condition << std::endl;
	}
	file_bc.close();

	//MODE - NODES - NOEL - NBC - NGEN - NWAVES - NPORTS - FSTART - FSTOP - NFR
	file_basic << dom.sc.mode_of_operation << " " << dom.nodes.size() << " " << dom.elements.size() << " " << dom.sc.nbc << " " << dom.sc.ngen << " " << dom.sc.numberOfWaves << " " << dom.sc.nports << " " << dom.sc.fstart << " " << dom.sc.fstop << " " << dom.sc.nfr << std::endl;
	file_basic.close();
}




void make_all_indices(std::vector<int>& all_indices, int& start, Domain& dom) {
	//note: 0, 2, 6, 8, 18, 20, 24, 26 are always included
	for (auto it = all_indices.begin(); it != all_indices.end(); ++it) {
		if (*it == 0) {
			//need to find a new value
			int index = it - all_indices.begin();
			int loc, loc6, loc8;
			switch (index) {
			case 1:
				break;
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[2]]) / 2.0);
				*it = loc;
			case 3:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[6]]) / 2.0);
				*it = loc;
				break;
			case 4:
				if (all_indices[5] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[3]] + dom.nodes[all_indices[5]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[8]]) / 2.0);
					all_indices[5] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[3]] + dom.nodes[all_indices[5]]) / 2.0);
					*it = loc;
				}
				break;
			case 5:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[8]]) / 2.0);
				*it = loc;
				break;

			case 7:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[8]]) / 2.0);
				*it = loc;
				break;

			case 9:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[18]]) / 2.0);
				*it = loc;
				break;
			case 10:
				if (all_indices[11] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[11]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[20]]) / 2.0);
					all_indices[11] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[11]]) / 2.0);
					*it = loc;
				}
				break;
			case 11:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[20]]) / 2.0);
				*it = loc;
				break;
			case 12:
				if (all_indices[15] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[15]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[24]]) / 2.0);
					all_indices[15] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[15]]) / 2.0);
					*it = loc;
				}
				break;
			case 13:
				if (all_indices[14] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
					*it = loc;
				}
				else {
					if (all_indices[17] != 0) {
						loc6 = dom.nodes.size();
						dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
						all_indices[14] = loc6;
						loc = dom.nodes.size();
						dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
						*it = loc;
					}
					else {
						loc8 = dom.nodes.size();
						dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
						all_indices[17] = loc8;
						loc6 = dom.nodes.size();
						dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
						all_indices[14] = loc6;
						loc = dom.nodes.size();
						dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
						*it = loc;
					}
				}
				break;
			case 14:
				if (all_indices[17] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
					all_indices[17] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
					*it = loc;
				}
				break;
			case 15:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[24]]) / 2.0);
				*it = loc;
				break;
			case 16:
				if (all_indices[17] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[15]] + dom.nodes[all_indices[17]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
					all_indices[17] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[5]] + dom.nodes[all_indices[17]]) / 2.0);
					*it = loc;
				}
				break;
			case 17:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
				*it = loc;
				break;

			case 19:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[18]] + dom.nodes[all_indices[20]]) / 2.0);
				*it = loc;
				break;

			case 21:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[18]] + dom.nodes[all_indices[24]]) / 2.0);
				*it = loc;
				break;
			case 22:
				if (all_indices[25] != 0) {
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[19]] + dom.nodes[all_indices[25]]) / 2.0);
					*it = loc;
				}
				else {
					loc6 = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[24]] + dom.nodes[all_indices[26]]) / 2.0);
					all_indices[25] = loc6;
					loc = dom.nodes.size();
					dom.nodes.push_back((dom.nodes[all_indices[19]] + dom.nodes[all_indices[25]]) / 2.0);
					*it = loc;
				}
				break;
			case 23:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[20]] + dom.nodes[all_indices[26]]) / 2.0);
				*it = loc;
				break;
			case 25:
				loc = dom.nodes.size();
				dom.nodes.push_back((dom.nodes[all_indices[24]] + dom.nodes[all_indices[26]]) / 2.0);
				*it = loc;
				break;



			}
		}
		/*for (int i = 0; i < 3; ++i) {
		Point p1, p2, p3, p4, p5, p6, p7, p8, p9;
		switch (i) {
		case 0:
		p1 = nodes[vertex_indices[0]]];
		p3 = nodes[vertex_indices[1]];
		p7 = nodes[vertex_indices[2]];
		p9 = nodes[vertex_indices[3]];
		break;
		case 1:
		p1 = (nodes[vertex_indices[0]] + nodes[vertex_indices[4]]) / 2.0;
		p3 = (nodes[vertex_indices[1]] + nodes[vertex_indices[5]]) / 2.0;
		p7 = (nodes[vertex_indices[2]] + nodes[vertex_indices[6]]) / 2.0;
		p9 = (nodes[vertex_indices[3]] + nodes[vertex_indices[7]]) / 2.0;
		break;
		case 2:
		p1 = nodes[vertex_indices[4]];
		p3 = nodes[vertex_indices[5]];
		p7 = nodes[vertex_indices[6]];
		p9 = nodes[vertex_indices[7]];
		break;
		}





		}*/
	}
}


void refinement::basis_check(Domain& dom) {
	//returns data structure for extra unknowns???
	std::vector<std::vector<adBasis>> new_basis(dom.elements.size());
	for (auto e = dom.elements.begin(); e != dom.elements.end(); ++e) {
		
		for (auto fac = e->facet_indices.begin(); fac != e->facet_indices.end(); ++fac) {
			int local_index_1;
			int local_index_2;
			int el_1 = dom.facets[*fac-1].element_indices.first;
			int el_2 = dom.facets[*fac-1].element_indices.second;
			if (e->index != el_1) {
			el_2 = el_1; 
			el_1 = e->index;
			local_index_1 = dom.facets[*fac-1].local_facet_indices.second;
			local_index_2 = dom.facets[*fac-1].local_facet_indices.first;
			
			}
			else {
				local_index_1 = dom.facets[*fac-1].local_facet_indices.first;
				local_index_2 = dom.facets[*fac-1].local_facet_indices.second;
			}
			if (local_index_1 < 1 || local_index_2 < 1) continue;
			int order_diff = dom.elements[el_2-1].expansion[ceil(local_index_2/2.0) - 1] - e->expansion[ceil(local_index_1 / 2.0) - 1];
			if (order_diff > 0) {
				//need to record that we have additional unknowns to add
				int dir1;
				int dir2;
				int pol_1 = local_index_1 % 2;
				int pol_2 = local_index_2 % 2;
				if (local_index_1 == 1 || local_index_1 == 2) { //-+u face if minus face,need one that is non zero on minus
					//need to add v,w basis fcns
				
					for (int ord = 0; ord < order_diff; ++ord) {
						//new_basis[e->index - 1].push_back(adBasis(e->index, 2, e->expansion[0] + ord, (pol_1 + 1) % 2, 0));
						//new_basis[e->index - 1].push_back(adBasis(e->index, 3, e->expansion[0] + ord, 0, (pol_1 + 1) % 2));
						for (int inner = 0; inner <= order_diff+e->expansion[0]; ++inner ){
							new_basis[e->index - 1].push_back(adBasis(e->index, 2, (pol_1 + 1) % 2, e->expansion[0] + ord, inner));
							new_basis[e->index - 1].push_back(adBasis(e->index, 2, (pol_1 + 1) % 2, inner, e->expansion[0] + ord));
							new_basis[e->index - 1].push_back(adBasis(e->index, 3, (pol_1 + 1) % 2, inner, e->expansion[0] + ord));
							new_basis[e->index - 1].push_back(adBasis(e->index, 3, (pol_1 + 1) % 2, e->expansion[0] + ord, inner));
						}
						
						
						//add other basis functions
					}
				}
				else if (local_index_1 == 3 || local_index_1 == 4) { //-+v face
					//need to add u,w basis fcns
					for (int ord = 0; ord < order_diff; ++ord) {
						//new_basis[e->index - 1].push_back(adBasis(e->index, 1,  (pol_1 + 1) % 2, e->expansion[1] + ord, 0));
						//new_basis[e->index - 1].push_back(adBasis(e->index, 3,  0, e->expansion[1] + ord, (pol_1 + 1) % 2));
						for (int inner = 0; inner <= order_diff + e->expansion[0]; ++inner) {
							new_basis[e->index - 1].push_back(adBasis(e->index, 1, e->expansion[1] + ord, (pol_1 + 1) % 2, inner));
							new_basis[e->index - 1].push_back(adBasis(e->index, 3, inner, (pol_1 + 1) % 2, e->expansion[1] + ord));
							new_basis[e->index - 1].push_back(adBasis(e->index, 1, inner, (pol_1 + 1) % 2, e->expansion[1] + ord));
							new_basis[e->index - 1].push_back(adBasis(e->index, 3, e->expansion[1] + ord, (pol_1 + 1) % 2, inner));
						}
					}
				}
				else { //-+w face
					//need to add u,v basis fcns
					for (int ord = 0; ord <= order_diff; ++ord) {
						//new_basis[e->index - 1].push_back(adBasis(e->index, 1, (pol_1 + 1) % 2, 0, e->expansion[1] + ord));
						//new_basis[e->index - 1].push_back(adBasis(e->index, 2, 0, (pol_1 + 1) % 2, e->expansion[1] + ord));
						for (int inner = 0; inner <= order_diff + e->expansion[0]; ++inner) {
							new_basis[e->index - 1].push_back(adBasis(e->index, 1, e->expansion[2] + ord, inner, (pol_1 + 1) % 2));
							new_basis[e->index - 1].push_back(adBasis(e->index, 2, inner, e->expansion[2] + ord, (pol_1 + 1) % 2));
							new_basis[e->index - 1].push_back(adBasis(e->index, 1, inner, e->expansion[2] + ord, (pol_1 + 1) % 2));
							new_basis[e->index - 1].push_back(adBasis(e->index, 2, e->expansion[2] + ord, inner, (pol_1 + 1) % 2));
						}
					}
				}
			}
		}
	}
	dom.p_ref_basis = new_basis;
}

void refinement::basis_add(Domain& dom, std::vector<int> refine_elements, int increaser ) {
	//returns data structure for extra unknowns???
	std::vector<std::vector<adBasis>> new_basis(dom.elements.size());
	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
		dom.elements[*it].refine = true;
		//add the candy functions
		//for u
			for (int ord = 0; ord < dom.elements[*it].expansion[0] + increaser; ++ord) {
				//new_basis[e->index - 1].push_back(adBasis(e->index, 2, e->expansion[0] + ord, (pol_1 + 1) % 2, 0));
				//new_basis[e->index - 1].push_back(adBasis(e->index, 3, e->expansion[0] + ord, 0, (pol_1 + 1) % 2));
				for (int vb = 0; vb <= dom.elements[*it].expansion[1] + increaser; ++vb) {
					for (int wb = 0; wb <= dom.elements[*it].expansion[2] + increaser; ++wb) {
						if (ord < dom.elements[*it].expansion[0]) { //if v was already included, only take higher order
							if (vb <= dom.elements[*it].expansion[1] && wb <= dom.elements[*it].expansion[2]) continue;
							else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 1, ord, vb, wb));
						}
						else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 1, ord, vb, wb));
					}
				}
		}
			//for v
			for (int ord = 0; ord < dom.elements[*it].expansion[1] + increaser; ++ord) {
			
				for (int ub = 0; ub <= dom.elements[*it].expansion[0] + increaser; ++ub) {
					for (int wb = 0; wb <= dom.elements[*it].expansion[2] + increaser; ++wb) {
						if (ord < dom.elements[*it].expansion[1]) { //if v was already included, only take higher order
							if (ub <= dom.elements[*it].expansion[0] && wb <= dom.elements[*it].expansion[2]) continue;
							else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 2, ub, ord, wb));
							
						}
						else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 2, ub, ord, wb));
					}
				}
			}
			//for w
			for (int ord = 0; ord < dom.elements[*it].expansion[2] + increaser; ++ord) {
				
				for (int ub = 0; ub <= dom.elements[*it].expansion[0] + increaser; ++ub) {
					for (int vb = 0; vb <= dom.elements[*it].expansion[1] + increaser; ++vb) {
						if (ord < dom.elements[*it].expansion[2]) { //if w was already included, only take higher order
							if (ub <= dom.elements[*it].expansion[0] && vb <= dom.elements[*it].expansion[1]) continue;
							else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 3, ub, vb, ord));
						}
						else new_basis[dom.elements[*it].index - 1].push_back(adBasis(dom.elements[*it].index, 3, ub, vb, ord));
					}
				}
			}
	}
	dom.p_ref_basis = new_basis;
}

void refinement::adjacent_check(Domain& dom) {
	//returns data structure for extra unknowns???
	for (auto e = dom.elements.begin(); e != dom.elements.end(); ++e) {
		if (e->refine == false) continue;
		for (auto fac = e->facet_indices.begin(); fac != e->facet_indices.end(); ++fac) {
			int local_index_1;
			int local_index_2;
			int el_1 = dom.facets[*fac - 1].element_indices.first;
			int el_2 = dom.facets[*fac - 1].element_indices.second;
			
			if (e->index != el_1) {
				el_2 = el_1;
				el_1 = e->index;
				local_index_1 = dom.facets[*fac - 1].local_facet_indices.second;
				local_index_2 = dom.facets[*fac - 1].local_facet_indices.first;

			}
			else {
				local_index_1 = dom.facets[*fac - 1].local_facet_indices.first;
				local_index_2 = dom.facets[*fac - 1].local_facet_indices.second;
			}
			if (el_2 < 0) continue;
			if (dom.elements[el_2 - 1].refine == false) continue;
			if (local_index_1 < 1 || local_index_2 < 1) continue;
			switch (local_index_1) {
			case 1:
				e->adjacent_list.u_dir += 1;
				break;
			case 2:
				e->adjacent_list.u_dir += 2;
				break;
			case 3:
				e->adjacent_list.v_dir += 1;
				break;
			case 4:
				e->adjacent_list.v_dir += 2;
				break;
			case 5:
				e->adjacent_list.w_dir += 1;
				break;
			case 6:
				e->adjacent_list.w_dir += 2;
				break;
			}
		}
	}
	
}


