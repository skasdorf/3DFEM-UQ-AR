#include "element_error.h"

//Point operator*(const Point& p1, const double scalar) {
//	return Point(p1.x*scalar, p1.y*scalar, p1.z*scalar);
//}
//Point operator/(Point& p1, const double scalar) {
//	return Point(p1.x / scalar, p1.y / scalar, p1.z / scalar);
//}
//Point operator+(Point& p1, const Point p2) {
//	return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
//}
//struct entry {
//public:
//	double val;
//	int row;
//	int col;
//	entry(double vali, int rowi, int coli) {
//		val = vali;
//		row = rowi;
//		col = coli;
//	}
//};


std::complex<double> qoi_error::check_q(std::string cAlphaAdjointFile, std::string cGrFile) {
	std::ifstream file;
	std::string line;
	std::vector<std::complex<double>> cGrForward, cAlphaAdjoint;

	file.open(cAlphaAdjointFile);
	while (getline(file, line)) {
		auto line_input = functions::split(line, ' ');
		cAlphaAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

	}
	file.close();
	file.open(cGrFile);
	while (getline(file, line)) {
		auto line_input = functions::split(line, ' ');
		cGrForward.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	file.close();
	std::complex<double> qoi_error = 0;
	for (int i = 0; i < cGrForward.size(); ++i) {
		qoi_error += std::conj(cAlphaAdjoint[i])*cGrForward[i];
	}
	return qoi_error;
}

//std::complex<double> qoi_error::add_qoi_error(std::string qoi_error_file, std::complex<double> qoi){
//	std::ifstream file;
//	std::string line;
//	std::complex<double> qoi_error;
//	file.open(qoi_error_file);
//	while (getline(file, line)) {
//		auto line_input = functions::split(line, ' ');
//		qoi_error += std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1]));
//	}
//	file.close();
//	return qoi_error + qoi;
//}

//bool operator<(const entry& a, const entry&b) {
//	if (a.val < b.val) return true;
//	else return false;
//}
//std::vector<int> qoi_error::pick_elements_to_refine(bool p_refine, std::string qoi_errors_file, const int& n){
//	std::ifstream file;
//	file.open(qoi_errors_file);
//	std::string line;
//	std::vector < std::complex<double>> qoi_errors;
//	while (getline(file, line)) {
//		auto line_input = functions::split(line, ' ');
//		qoi_errors.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
//	}
//	//Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> k= Eigen::Map<Eigen::Matrix<std::complex<double>>, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
//	Eigen::VectorXcd k = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(qoi_errors.data(), qoi_errors.size());
//	std::vector<entry> cVec_vec;
//	//Eigen::Matrix<entry, Eigen::Dynamic, 1> cVec;
//	auto vert = k.replicate(1, qoi_errors.size());
//	auto horz = vert.transpose();
//	auto C = (vert + horz).cwiseAbs();
//	
//	for (int i = 0; i < C.rows(); ++i) {
//		for (int j = 0; j < C.cols(); ++j) {
//			if (i >= j) continue;
//			
//			cVec_vec.push_back(entry(C(i, j), i, j));
//		}
//	}
//	std::sort(cVec_vec.rbegin(), cVec_vec.rend());
//	std::vector<int> refine_elements;
//	int num_pairs = 0;
//	int pair_max = n;
//	for (auto it = cVec_vec.begin(); it != cVec_vec.end(); ++it) {
//		
//		int e1 = it->col;
//		int e2 = it->row;
//		if (!(std::find(refine_elements.begin(), refine_elements.end(), e1) != refine_elements.end())
//			&& !(std::find(refine_elements.begin(), refine_elements.end(), e2) != refine_elements.end())) {
//			refine_elements.push_back(e1);
//			refine_elements.push_back(e2);
//			++num_pairs;
//		}
//		if (num_pairs >= pair_max) break;
//	}
//	return refine_elements;
//}

//void qoi_error::p_refine(std::vector<int> refine_elements, Domain& dom)
//{
//	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
//		dom.elements[*it - 1].refine = true;
//		dom.elements[*it - 1].expansion[0] = 4; //note: rn assumes starting order is 3, only sets to 4 for repeatability
//		dom.elements[*it - 1].expansion[1] = 4;
//		dom.elements[*it - 1].expansion[2] = 4;
//	}
//}

//void qoi_error::track_refine(std::vector<int> refine_elements, Domain& dom)
//{
//	for (int i = 0; i < dom.elements.size(); ++i) {
//		Element e1 = dom.elements[i];
//		int refine_sum;
//		for (auto fac = e1.facet_indices.begin(); fac != e1.facet_indices.end(); ++fac) {
//			Facet f1 = dom.facets[*fac-1];
//			Element e2 = dom.elements[f1.element_indices.second - 1];
//			if (f1.element_indices.second == e1.index) e2 = dom.elements[f1.element_indices.first - 1];
//			
//		}
//	}
//}

//Point operator+(const Point& p1, const Point& p2) {
//	return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
//}
//Point operator-(const Point& p1, const Point& p2) {
//	return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
//}
//Point operator*(const double scal, const Point& p1) {
//	return Point(scal*p1.x, scal*p1.y, scal*p1.z);
//}
//Point find_halfway(Point p1, Point p2) {
//	//returns the point inbetween 2 points
//	Point p_ret;
//	p_ret = p1 + (1.0 / 2.0)*(p2 - p1);
//	return p_ret;
//}
//void findPowersLagrange(const int& kuvw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
//	matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr, const int& ending) {
//	//right now we have kuvw -- for now the order is the same in each direction -- this could be generalized later
//
//	double x, xj, xm, f;
//	for (int i = 0; i < ending; ++i) {
//
//		x = xglu[i]; //Gauss-Legendre Integration point
//		for (int m = 0; m <= kuvw; m++) {
//			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//			f = 1.0;
//			for (int j = 0; j <= kuvw; j++) {
//				xj = (2.0*j - kuvw) / kuvw;
//				if (j != m) {
//					f = f*(x - xj) / (xm - xj);
//				}
//			}//for j
//			fuPowersLagr(i, m) = f;
//		}//for m
//
//
//
//		x = xglv[i]; //Gauss-Legendre Integration point
//		for (int m = 0; m <= kuvw; m++) {
//			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//			f = 1.0;
//			for (int j = 0; j <= kuvw; j++) {
//				xj = (2.0*j - kuvw) / kuvw;
//				if (j != m) {
//					f = f*(x - xj) / (xm - xj);
//				}
//			}//for j
//			fvPowersLagr(i, m) = f;
//		}//for m
//
//
//
//		x = xglw[i]; //Gauss-Legendre Integration point
//		for (int m = 0; m <= kuvw; m++) {
//			xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//			f = 1.0;
//			for (int j = 0; j <= kuvw; j++) {
//				xj = (2.0*j - kuvw) / kuvw;
//				if (j != m) {
//					f = f*(x - xj) / (xm - xj);
//				}
//			}//for j
//			fwPowersLagr(i, m) = f;
//		}//for m
//	}
//	
//
//}//void findPowersLagrange
//double trilagrangeoduvw(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr) {
//	return fuPowersLagr(m, i)*fvPowersLagr(n, j)*fwPowersLagr(l, k);
//}
//void find_eps_mu_matrix(const int& kuvw, const int& num_coords, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr,
//	const int& size1, const std::vector<std::vector<std::complex<double>>>& muRel, matrix2d<std::complex<double>>& muRelInt) {
//	//std::cout << "entered eps_mu_matrix" << std::endl;
//	std::complex<double> cMu;
//	for (auto s = 1; s <= size1; ++s) { //matrix entries
//		for (int m = 0; m < num_coords; ++m) {
//			//for (int n = 0; n <= nglv + 1; ++n) {
//			//	for (int l = 0; l <= nglw + 1; ++l) {
//					cMu = std::complex<double>(0, 0);
//					for (int mat_loc = 0; mat_loc < muRel.size(); ++mat_loc) { //nodes
//						int kuvwp = kuvw + 1;
//						double f = trilagrangeoduvw(m, m, m, mat_loc % (kuvwp), int(floor(mat_loc / kuvwp)) % kuvwp,
//							int(floor(mat_loc / (kuvwp*kuvwp))), fuPowersLagr, fvPowersLagr, fwPowersLagr);
//
//						muRelInt(s, m) += muRel[mat_loc][s - 1] * f;
//					}
//					//	muRelInt(s, m, n, l) = cMu;
//				//}
//		//	}
//
//		}
//
//	}
//	//std::cout << "exited eps_mu_match" << std::endl;
//}
////void copy_mats(std::vector<std::complex<double>>& eps, std::vector<std::complex<double>>& mur) {
////
////}
////void h_refine2(std::vector<int> refine_elements, Domain& dom, double factor) {
////	std::vector<Element> new_elements;
////	std::vector<Point> p_central;
////	p_central.resize(27);
////	//build central cube CORNER VERTICES
////	p_central[0] = (Point(-1.0*factor, -1.0*factor, -1.0*factor));
////	p_central[2] = (Point(1.0*factor, -1.0*factor, -1.0*factor));
////	p_central[6] = (Point(-1.0*factor, 1.0*factor, -1.0*factor));
////	p_central[8] = (Point(1.0*factor, 1.0*factor, -1.0*factor));
////	p_central[18] = (Point(-1.0*factor, -1.0*factor, 1.0*factor));
////	p_central[20] = (Point(1.0*factor, -1.0*factor, 1.0*factor));
////	p_central[24] = (Point(-1.0*factor, 1.0*factor, 1.0*factor));
////	p_central[26] = (Point(1.0*factor, 1.0*factor, 1.0*factor));
////	//build just the vertices on faces
////	p_central[1] = find_halfway(p_central[0], p_central[2]);
////	p_central[3] = find_halfway(p_central[0], p_central[6]);
////	p_central[4] = find_halfway(p_central[0], p_central[8]);
////	p_central[5] = find_halfway(p_central[2], p_central[8]);
////	p_central[7] = find_halfway(p_central[6], p_central[8]);
////	p_central[9] = find_halfway(p_central[0], p_central[18]);
////	p_central[10] = find_halfway(p_central[0], p_central[20]);
////	p_central[11] = find_halfway(p_central[2], p_central[20]);
////	p_central[12] = find_halfway(p_central[0], p_central[24]);
////	p_central[14] = find_halfway(p_central[2], p_central[26]);
////	p_central[15] = find_halfway(p_central[6], p_central[24]);
////	p_central[16] = find_halfway(p_central[6], p_central[26]);
////	p_central[17] = find_halfway(p_central[8], p_central[26]);
////	p_central[19] = find_halfway(p_central[18], p_central[20]);
////	p_central[21] = find_halfway(p_central[18], p_central[24]);
////	p_central[22] = find_halfway(p_central[18], p_central[26]);
////	p_central[23] = find_halfway(p_central[20], p_central[26]);
////	
////	p_central[25] = find_halfway(p_central[24], p_central[26]);
////	//need to make the intermediate elements for the 6 additional elems
////	std::vector<Point> interm_points(27);
////	for (int lines = 0; lines < p_central.size(); ++lines) {
////		if (lines == 13) continue; //this node is inside
////		interm_points[lines] = (find_halfway(p_central[lines], (1.0/factor)*p_central[lines]));
////	}
////	std::vector<std::vector<int>> elem_pts(6); //stored per face index
////	for (int i = 0; i < 6; ++i) {
////		//std::vector<int> elem_cent_pts, elem_interm_pts, elem_interp_pts;
////
////		 //account for pts on the central element, the interp points from the orig eleme and the new interm pts
////
////		//elem zero gets 0,3,6,9,12,15,18,21,24
////		//elem one gets 2,5,8,11,14,17,20,23,26
////		/*
////		elem0 and 1 gets lowest + 3 each time for 9 vertices from central
////
////		elem 2 and 3 gets lowest, lowest + 1, lowest+2, then increment by 9 and repeat twice
////
////		elem 4 and 5 gets lowest up to lowest + 8
////		*/
////		switch (i) {
////		case 0:
////		case 1:
////			int start = i * 2;
////			for (int k = 0; k < 9; ++k) {
////				elem_pts[i].push_back(start + 3 * k);
////			}
////			break;
////		case 2:
////		case 3:
////			int start = (i - 2) * 6;
////			int adjust = 0;
////			for (int k = 0; k < 3; ++k) {
////				adjust = k * 9;
////				elem_pts[i].push_back(start + adjust);
////				elem_pts[i].push_back(start + adjust + 1);
////				elem_pts[i].push_back(start + adjust + 2);
////			}
////			break;
////		case 4:
////		case 5:
////			int start = (i - 4) * 18;
////			for (int k = 0; k < 9; ++k) {
////				elem_pts[i].push_back(start + k);
////			}
////			break;
////		}
////	}
////	//iterate through elements now
////	std::vector<Element> overall_new_elements;
////	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
////		//iterate through elements to refine
////		//in an element, take a smaller sized cube in the parametric domain in its center
////		//map the vertices of these coordinates to the real domain
////		//connect the vertices of each face to each face, forming new elements
////		//preserve the interp points!
////		Element e1 = dom.elements[*it - 1];
////		std::vector<Element> new_elements[6];
////
////		int kuvw;
////		int size1;
////		if (e1.materials.icode == 0) {
////			kuvw = e1.materials.KuvwA;
////			size1 = 9 - 3 * e1.materials.sym;
////		}
////		else {
////			kuvw = e1.materials.Kuvw;
////			size1 = 1;
////		}
////		std::vector<double> x_par, y_par, z_par;
////		for (int par_it = 0; par_it < p_central.size(); ++par_it) {
////			x_par.push_back(p_central[par_it].x);
////			y_par.push_back(p_central[par_it].y);
////			z_par.push_back(p_central[par_it].z);
////		}
////		for (int par_it = 0; par_it < interm_points.size(); ++par_it) {
////			x_par.push_back(interm_points[par_it].x);
////			y_par.push_back(interm_points[par_it].y);
////			z_par.push_back(interm_points[par_it].z);
////		}
////		matrix2d < std::complex<double>> epsrNew = matrix2d<std::complex<double>>(size1, x_par.size());
////		matrix2d < std::complex<double>> murNew = matrix2d<std::complex<double>>(size1, x_par.size());
////		if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
////			epsrNew(1,1) = e1.materials.epsr_list[0][0];
////			murNew(1,1) = e1.materials.mur_list[0][0];
////		}
////		else {
////			//material params at one face already known, get material params at other locations
////			matrix2d<double> fuPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
////			matrix2d<double> fvPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
////			matrix2d<double> fwPowersLagr = matrix2d<double>(x_par.size() + 1, kuvw);
////			
////			findPowersLagrange(kuvw, x_par, y_par, z_par, fuPowersLagr, fvPowersLagr, fwPowersLagr, x_par.size());
////			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.epsr_list, epsrNew);
////			find_eps_mu_matrix(kuvw, x_par.size(), fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e1.materials.mur_list, murNew);
////		}
////		//add the points to the original set of points
////		std::vector<double> r(4);
////		std::vector<Point> real_domain;
////		for (auto p = p_central.begin(); p != p_central.end(); ++p) {
////			if (p - p_central.begin() == 13) continue;
////			r = unitVectorsM::findR(e1.geom_order, p->x, p->y, p->z, e1.rs, e1.nRs);
////			real_domain.push_back(Point(r[1], r[2], r[3]));
////		}
////		for (auto p = interm_points.begin(); p != interm_points.end(); ++p) {
////			r = unitVectorsM::findR(e1.geom_order, p->x, p->y, p->z, e1.rs, e1.nRs);
////			real_domain.push_back(Point(r[1], r[2], r[3]));
////		}
////		int cent_start = dom.nodes.size();
////		int interm_start = dom.nodes.size() + 25;
////		dom.nodes.insert(std::end(dom.nodes), std::begin(real_domain), std::end(real_domain));
////		for (int e_new = 0; e_new < 6; ++e_new) {
////			Element e_test;
////			std::vector<int> order_pt_indices;
////			std::vector<std::vector<std::complex<double>>> new_epsr;
////			std::vector<std::vector<std::complex<double>>> new_mur;
////			bool do_extra_mats = false;
////			if (e1.materials.hcode == 1 && e1.materials.icode == 1) {
////				new_epsr = e1.materials.epsr_list;
////				new_mur = e1.materials.mur_list;
////			}
////			else {
////				do_extra_mats = true;
////				new_epsr = std::vector<std::vector<std::complex<double>>>(e1.materials.epsr_list.size(), std::vector<std::complex<double>>(e1.materials.epsr_list[0].size()));
////				new_mur = new_epsr;
////			}
////			if (e_new % 2 == 0) {
////
////				//is even
////				//order starts at orig
////				int loc = 0;
////				switch (e_new) {
////				case 0:
////					//orig, interm, central
////					//all x3
////					
////					for (int order = 0; order < 9; ++order) {
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order]]);
////						order_pt_indices.push_back(elem_pts[e_new][order] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order]+cent_start);
////						if (do_extra_mats) {
////							std::vector<std::complex<double>> temp_vec = (e1.materials.epsr_list[elem_pts[e_new][order]]);
////							//temp_vec.push_back;
////							new_epsr[order * 3] = temp_vec;
////							temp_vec = e1.materials.mur_list[elem_pts[e_new][order]];
////							new_mur[order * 3] = temp_vec;
////					}
////					break;
////				case 2:
////					//orig orig orig
////					//interm interm interm
////					//central central central
////					//all x 3
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order*3]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order*3 + 1]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order*3 + 2]]);
////						//
////						order_pt_indices.push_back(elem_pts[e_new][order*3] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order*3 + 1] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order*3 + 2] + interm_start);
////						//
////						order_pt_indices.push_back(elem_pts[e_new][order*3] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order*3 + 1] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order*3 + 2] + cent_start);
////					}
////					break;
////				case 4:
////					//orig orig orig x 3
////					//interm interm interm x 3
////					//central central central x 3
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order*3]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 1]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 2]]);
////					}
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
////					}
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
////					}
////
////
////					break;
////				}
////			}
////			else {
////				//is odd
////				//order starts at central
////				switch (e_new) {
////				case 1:
////					for (int order = 0; order < 9; ++order) {
////						
////						order_pt_indices.push_back(elem_pts[e_new][order] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order] + interm_start);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order]]);
////					}
////					break;
////				case 3:
////					for (int order = 0; order < 3; ++order) {
////						
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
////						//
////						
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
////						//
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 1]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 2]]);
////						
////						
////					}
////					break;
////				case 5:
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + cent_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + cent_start);
////					}
////					
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(elem_pts[e_new][order * 3] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 1] + interm_start);
////						order_pt_indices.push_back(elem_pts[e_new][order * 3 + 2] + interm_start);
////					}
////
////					for (int order = 0; order < 3; ++order) {
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 1]]);
////						order_pt_indices.push_back(e1.node_indices[elem_pts[e_new][order * 3 + 2]]);
////					}
////					
////					break;
////				}
////			}
////			//now have vertices matched up properly
////			//turn materials matrix into the proper form
////			
////
////		}
////
////		//do the inner element last
////	}
////	//print out all necessary information to then redo the mesh
////	std::ofstream file_elem, file_geom, file_bc;
////	file_geom.open("../Geometry_refined.dat");
////	file_geom << "!Refined geometry!" << std::endl;
////	file_elem.open("../Element_refined.dat");
////	file_elem << "!Refined elements!" << std::endl;
////	file_bc.open("../Boundary_refined.dat");
////	file_bc << "!Refined boundary!" << std::endl;
////	for (auto dom_points = dom.nodes.begin(); dom_points != dom.nodes.end(); ++dom_points) {
////		file_geom << dom_points->x << " " << dom_points->y << " " << dom_points->z << std::endl;
////	}
////	file_geom.close();
////	for (auto elem_p = dom.elements.begin(); elem_p != dom.elements.end(); ++elem_p) {
////		file_elem << "Index: " << elem_p->index << std::endl;
////		file_elem << "Nodes: ";
////		for (auto elm_nodes = elem_p->all_indices.begin(); elm_nodes != elem_p->all_indices.end(); ++elm_nodes) {
////			file_elem << *elm_nodes << " ";
////		}
////		file_elem << std::endl << "Basis orders: " << elem_p->expansion[0] << " " << elem_p->expansion[1] << elem_p->expansion[2] << std::endl;
////		file_elem << std::endl << "Quadrature orders: " << elem_p->quadrature[0] << " " << elem_p->quadrature[1] << " " << elem_p->quadrature[2] << std::endl;
////		file_elem << std::endl << "Materials: " << elem_p->materials.hcode << " " << elem_p->materials.icode
////			<< " " << elem_p->materials.pmlcode << " " << elem_p->materials.Kuvw << elem_p->materials.sym << std::endl;
////		file_elem << "Eprs: ";// << std::endl;
////		for (auto mat_ent = elem_p->materials.epsr_list.begin(); mat_ent != elem_p->materials.epsr_list.end(); ++mat_ent) {
////			file_elem << std::endl;
////			//iterates through nodes
////			for (auto val_mat = mat_ent->begin(); val_mat != mat_ent->end(); ++val_mat) {
////				//gets values
////				file_elem << val_mat->real() << "," << val_mat->imag() << " ";
////			}
////		}
////	}
////	file_elem.close();
////	//save boundary condition info
////	for (auto fac_bc = dom.facets.begin(); fac_bc != dom.facets.end(); ++fac_bc) {
////		if (fac_bc->boundary_condition == 0) continue;
////		file_bc << fac_bc->vertices[0] << " " << fac_bc->vertices[1] << " " << fac_bc->vertices[2] << " " << fac_bc->vertices[3] << " " << fac_bc->boundary_condition << std::endl;
////	}
////	
////}


//void qoi_error::h_refine(std::vector<int> refine_elements, Domain & dom, double factor)
//{
//	//double factor = 1.0 / 3.0;
//	std::vector<Element> new_elements;
//	std::vector<Point> p_points;
//	/*Point p1(-1.0*factor, -1.0*factor, -1.0*factor);
//	Point p2(1.0*factor, -1.0*factor, -1.0*factor);
//	Point p3(-1.0*factor, 1.0*factor, -1.0*factor);
//	Point p4(1.0*factor, 1.0*factor, -1.0*factor);
//
//	Point p5(-1.0*factor, -1.0*factor, 1.0*factor);
//	Point p6(1.0*factor, -1.0*factor, 1.0*factor);
//	Point p7(-1.0*factor, 1.0*factor, 1.0*factor);
//	Point p8(1.0*factor, 1.0*factor, 1.0*factor);*/
//
//	p_points.push_back(Point(-1.0*factor, -1.0*factor, -1.0*factor));
//	p_points.push_back(Point(1.0*factor, -1.0*factor, -1.0*factor));
//	p_points.push_back(Point(-1.0*factor, 1.0*factor, -1.0*factor));
//	p_points.push_back(Point(1.0*factor, 1.0*factor, -1.0*factor));
//
//	p_points.push_back(Point(-1.0*factor, -1.0*factor, 1.0*factor));
//	p_points.push_back(Point(1.0*factor, -1.0*factor, 1.0*factor));
//	p_points.push_back(Point(-1.0*factor, 1.0*factor, 1.0*factor));
//	p_points.push_back(Point(1.0*factor, 1.0*factor, 1.0*factor));
//	std::vector<Point> real_domain;
////	facetContainer facc;
//	for (auto it = refine_elements.begin(); it != refine_elements.end(); ++it) {
//		//iterate through elements to refine
//		//in an element, take a smaller sized cube in the parametric domain in its center
//		//map the vertices of these coordinates to the real domain
//		//connect the vertices of each face to each face, forming new elements
//		//preserve the interp points!
//		Element e1 = dom.elements[*it - 1];
//		std::vector<double> r(4);
//		for (auto p = p_points.begin(); p != p_points.end(); ++p) {
//			r = unitVectorsM::findR(e1.geom_order, p->x, p->y, p->z, e1.rs, e1.nRs);
//			real_domain.push_back(Point(r[1], r[2], r[3]));
//		}
//		//partion old interp points
//		//add the new vertices
//		int starting = dom.nodes.size();
//		dom.nodes.insert(dom.nodes.end(), real_domain.begin(), real_domain.end());
//	
//		int num_elements = dom.elements.size();
//		std::vector<Element> new_elements;
//		for (int el_i = 0; el_i < 6; ++el_i) {
//			Element e_in;
//			e_in.index = (++num_elements);
//			e_in.type = e1.type;
//			e_in.geom_order = e1.geom_order;
//			e_in.quadrature = e1.quadrature;
//			e_in.expansion = e1.expansion;
//			e_in.all_indices.resize((e_in.geom_order + 1)*(e_in.geom_order + 1)*(e_in.geom_order + 1));
//			if (el_i == 0) {
//				//load the copied face
//				e_in.all_indices[0] = e1.all_indices[0];
//				e_in.all_indices[3] = e1.all_indices[3];
//				e_in.all_indices[6] = e1.all_indices[6];
//				e_in.all_indices[9] = e1.all_indices[9];
//				e_in.all_indices[12] = e1.all_indices[12];
//				e_in.all_indices[15] = e1.all_indices[15];
//				e_in.all_indices[18] = e1.all_indices[18];
//				e_in.all_indices[21] = e1.all_indices[21];
//				e_in.all_indices[24] = e1.all_indices[24];
//				//load 4 other vertices
//				e_in.all_indices[2] = starting;
//				e_in.all_indices[8] = starting+2;
//				e_in.all_indices[20] = starting+4;
//				e_in.all_indices[26] = starting+6;
//			}
//			else if (el_i == 1) {
//				e_in.all_indices[2] = e1.all_indices[2];
//				e_in.all_indices[5] = e1.all_indices[5];
//				e_in.all_indices[8] = e1.all_indices[8];
//				e_in.all_indices[11] = e1.all_indices[11];
//				e_in.all_indices[14] = e1.all_indices[14];
//				e_in.all_indices[17] = e1.all_indices[17];
//				e_in.all_indices[20] = e1.all_indices[20];
//				e_in.all_indices[23] = e1.all_indices[23];
//				e_in.all_indices[26] = e1.all_indices[26];
//				//load 4 other vertices
//				e_in.all_indices[0] = starting+1;
//				e_in.all_indices[6] = starting + 3;
//				e_in.all_indices[18] = starting + 5;
//				e_in.all_indices[24] = starting + 7;
//			}
//			else if (el_i == 2) {
//				e_in.all_indices[2] = e1.all_indices[2];
//				e_in.all_indices[5] = e1.all_indices[5];
//				e_in.all_indices[8] = e1.all_indices[8];
//				e_in.all_indices[11] = e1.all_indices[11];
//				e_in.all_indices[14] = e1.all_indices[14];
//				e_in.all_indices[17] = e1.all_indices[17];
//				e_in.all_indices[20] = e1.all_indices[20];
//				e_in.all_indices[23] = e1.all_indices[23];
//				e_in.all_indices[26] = e1.all_indices[26];
//				//load 4 other vertices
//				e_in.all_indices[0] = starting + 1;
//				e_in.all_indices[6] = starting + 3;
//				e_in.all_indices[18] = starting + 5;
//				e_in.all_indices[24] = starting + 7;
//			}
//			else if (el_i == 3) {
//				e_in.all_indices[2] = e1.all_indices[2];
//				e_in.all_indices[5] = e1.all_indices[5];
//				e_in.all_indices[8] = e1.all_indices[8];
//				e_in.all_indices[11] = e1.all_indices[11];
//				e_in.all_indices[14] = e1.all_indices[14];
//				e_in.all_indices[17] = e1.all_indices[17];
//				e_in.all_indices[20] = e1.all_indices[20];
//				e_in.all_indices[23] = e1.all_indices[23];
//				e_in.all_indices[26] = e1.all_indices[26];
//				//load 4 other vertices
//				e_in.all_indices[0] = starting + 1;
//				e_in.all_indices[6] = starting + 3;
//				e_in.all_indices[18] = starting + 5;
//				e_in.all_indices[24] = starting + 7;
//			}
//			else if (el_i == 4) {
//				e_in.all_indices[2] = e1.all_indices[2];
//				e_in.all_indices[5] = e1.all_indices[5];
//				e_in.all_indices[8] = e1.all_indices[8];
//				e_in.all_indices[11] = e1.all_indices[11];
//				e_in.all_indices[14] = e1.all_indices[14];
//				e_in.all_indices[17] = e1.all_indices[17];
//				e_in.all_indices[20] = e1.all_indices[20];
//				e_in.all_indices[23] = e1.all_indices[23];
//				e_in.all_indices[26] = e1.all_indices[26];
//				//load 4 other vertices
//				e_in.all_indices[0] = starting + 1;
//				e_in.all_indices[6] = starting + 3;
//				e_in.all_indices[18] = starting + 5;
//				e_in.all_indices[24] = starting + 7;
//			}
//			else if (el_i == 5) {
//				e_in.all_indices[2] = e1.all_indices[2];
//				e_in.all_indices[5] = e1.all_indices[5];
//				e_in.all_indices[8] = e1.all_indices[8];
//				e_in.all_indices[11] = e1.all_indices[11];
//				e_in.all_indices[14] = e1.all_indices[14];
//				e_in.all_indices[17] = e1.all_indices[17];
//				e_in.all_indices[20] = e1.all_indices[20];
//				e_in.all_indices[23] = e1.all_indices[23];
//				e_in.all_indices[26] = e1.all_indices[26];
//				//load 4 other vertices
//				e_in.all_indices[0] = starting + 1;
//				e_in.all_indices[6] = starting + 3;
//				e_in.all_indices[18] = starting + 5;
//				e_in.all_indices[24] = starting + 7;
//			}
//
//			//set up known vertices;
//			
//		}
//		//output to file the necessary info
//		//e2.index = (++num_elements);
//		//e3.index = (++num_elements);
//		//e4.index = (++num_elements);
//		//e5.index = (++num_elements);
//		//e6.index = (++num_elements);
//		//e7.index = (++num_elements);
//
//		//e2.type = e1.type;
//		//e3.type = e1.type;
//		//e4.type = e1.type;
//		//e5.type = e1.type;
//		//e6.type = e1.type;
//		//e7.type = e1.type;
//
//		//e2.geom_order = e1.geom_order;
//		//e3.geom_order = e1.geom_order;
//		//e4.geom_order = e1.geom_order;
//		//e5.geom_order = e1.geom_order;
//		//e6.geom_order = e1.geom_order;
//		//e7.geom_order = e1.geom_order;
//
//		//e2.quadrature = e1.quadrature;
//		//e2.expansion = e1.expansion;
//		//e3.quadrature = e1.quadrature;
//		//e3.expansion = e1.expansion;
//		//e4.quadrature = e1.quadrature;
//		//e4.expansion = e1.expansion;
//		//e5.quadrature = e1.quadrature;
//		//e5.expansion = e1.expansion;
//		//e6.quadrature = e1.quadrature;
//		//e6.expansion = e1.expansion;
//		//e7.quadrature = e1.quadrature;
//		//e7.expansion = e1.expansion;
//		//
//		////need vertex indices, node indices, facet indices, all indices
//		//e2.vertex_indices = {e1.vertex_indices[0], starting, e1.vertex_indices[2], starting + 2, e1.vertex_indices[4], starting + 4, e1.vertex_indices[6], starting + 6};
//		//e3.vertex_indices = {  starting+1,e1.vertex_indices[1],  starting + 3,e1.vertex_indices[3],  starting + 5, e1.vertex_indices[5], starting + 7, e1.vertex_indices[7] };
//		//e4.vertex_indices = { e1.vertex_indices[0],  e1.vertex_indices[1],starting, starting + 1, e1.vertex_indices[4],  e1.vertex_indices[5],starting + 4, starting + 5 };
//		//e5.vertex_indices = {  starting + 2,  starting + 3, e1.vertex_indices[2],e1.vertex_indices[3],  starting + 6, starting + 7, e1.vertex_indices[6],e1.vertex_indices[7]  };
//		//e6.vertex_indices = { e1.vertex_indices[0],e1.vertex_indices[1],e1.vertex_indices[2],e1.vertex_indices[3], starting,  starting + 1,  starting + 2,  starting + 3 };
//		//e7.vertex_indices = { starting + 4, starting + 5, starting + 6, starting + 7, e1.vertex_indices[4],  e1.vertex_indices[5],  e1.vertex_indices[6],  e1.vertex_indices[7]};
//
//		//e2.all_indices = { e1.vertex_indices[0], starting, e1.vertex_indices[2], starting + 2, e1.vertex_indices[4], starting + 4, e1.vertex_indices[6], starting + 6 };
//
//	//	Element(index_temp, type, geometric_order, vertices_temp, nodes_temp, facets_temp, quadrature_temp, expansion_temp, all_temp));
//		//fileFacet f1 = fileFacet(starting + 0, starting+2, starting+4, starting+6, 1);
//		//fileFacet f2 = fileFacet(starting+1, starting+3, starting+5, starting+7, 2);
//		//fileFacet f3 = fileFacet(starting+0, starting+1, starting+4, starting+5, 3);
//		//fileFacet f4 = fileFacet(starting+2, starting+3, starting+6, starting+7, 4);
//		//fileFacet f5 = fileFacet(starting+0, starting+1, starting+2, starting+3, 5);
//		//fileFacet f6 = fileFacet(starting+4, starting+5, starting+6, starting+7, 6);
//
//		
//
//		/*Facet f1;
//		f1.vertices = { starting, starting + 2, starting + 4, starting + 6 };
//		f1.local_facet_indices.first = 1;
//		f1.index = dom.facets.size();
//		f1.type = 0;
//		f1.element_indices.first = e1.index;
//		f1.element_indices.second = dom.elements.size() + f1.local_facet_indices.first;*/
//
//
//
//
//	}
//}
//
//void make_all_indices(std::vector<int>& all_indices, int& start, Domain& dom) {
//	//note: 0, 2, 6, 8, 18, 20, 24, 26 are always included
//	for (auto it = all_indices.begin(); it != all_indices.end(); ++it) {
//		if (*it == 0) {
//			//need to find a new value
//			int index = it - all_indices.begin();
//			int loc, loc6, loc8;
//			switch (index) {
//			case 1:
//				break;
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[2]]) / 2.0);
//				*it = loc;
//			case 3:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[6]]) / 2.0);
//				*it = loc;
//				break;
//			case 4:
//				if (all_indices[5] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[3]] + dom.nodes[all_indices[5]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[8]]) / 2.0);
//					all_indices[5] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[3]] + dom.nodes[all_indices[5]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 5:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[8]]) / 2.0);
//				*it = loc;
//				break;
//
//			case 7:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[8]]) / 2.0);
//				*it = loc;
//				break;
//
//			case 9:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[0]] + dom.nodes[all_indices[18]]) / 2.0);
//				*it = loc;
//				break;
//			case 10:
//				if (all_indices[11] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[11]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[20]]) / 2.0);
//					all_indices[11] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[11]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 11:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[2]] + dom.nodes[all_indices[20]]) / 2.0);
//				*it = loc;
//				break;
//			case 12:
//				if (all_indices[15] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[15]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[24]]) / 2.0);
//					all_indices[15] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[9]] + dom.nodes[all_indices[15]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 13:
//				if (all_indices[14] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					if (all_indices[17] != 0) {
//						loc6 = dom.nodes.size();
//						dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
//						all_indices[14] = loc6;
//						loc = dom.nodes.size();
//						dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
//						*it = loc;
//					}
//					else {
//						loc8 = dom.nodes.size();
//						dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
//						all_indices[17] = loc8;
//						loc6 = dom.nodes.size();
//						dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
//						all_indices[14] = loc6;
//						loc = dom.nodes.size();
//						dom.nodes.push_back((dom.nodes[all_indices[12]] + dom.nodes[all_indices[14]]) / 2.0);
//						*it = loc;
//					}
//				}
//				break;
//			case 14:
//				if (all_indices[17] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
//					all_indices[17] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[11]] + dom.nodes[all_indices[17]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 15:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[6]] + dom.nodes[all_indices[24]]) / 2.0);
//				*it = loc;
//				break;
//			case 16:
//				if (all_indices[17] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[15]] + dom.nodes[all_indices[17]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
//					all_indices[17] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[5]] + dom.nodes[all_indices[17]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 17:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[8]] + dom.nodes[all_indices[26]]) / 2.0);
//				*it = loc;
//				break;
//
//			case 19:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[18]] + dom.nodes[all_indices[20]]) / 2.0);
//				*it = loc;
//				break;
//
//			case 21:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[18]] + dom.nodes[all_indices[24]]) / 2.0);
//				*it = loc;
//				break;
//			case 22:
//				if (all_indices[25] != 0) {
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[19]] + dom.nodes[all_indices[25]]) / 2.0);
//					*it = loc;
//				}
//				else {
//					loc6 = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[24]] + dom.nodes[all_indices[26]]) / 2.0);
//					all_indices[25] = loc6;
//					loc = dom.nodes.size();
//					dom.nodes.push_back((dom.nodes[all_indices[19]] + dom.nodes[all_indices[25]]) / 2.0);
//					*it = loc;
//				}
//				break;
//			case 23:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[20]] + dom.nodes[all_indices[26]]) / 2.0);
//				*it = loc;
//				break;
//			case 25:
//				loc = dom.nodes.size();
//				dom.nodes.push_back((dom.nodes[all_indices[24]] + dom.nodes[all_indices[26]]) / 2.0);
//				*it = loc;
//				break;
//
//
//
//			}
//		}
//		/*for (int i = 0; i < 3; ++i) {
//			Point p1, p2, p3, p4, p5, p6, p7, p8, p9;
//			switch (i) {
//			case 0:
//				p1 = nodes[vertex_indices[0]]];
//				p3 = nodes[vertex_indices[1]];
//				p7 = nodes[vertex_indices[2]];
//				p9 = nodes[vertex_indices[3]];
//				break;
//			case 1:
//				p1 = (nodes[vertex_indices[0]] + nodes[vertex_indices[4]]) / 2.0;
//				p3 = (nodes[vertex_indices[1]] + nodes[vertex_indices[5]]) / 2.0;
//				p7 = (nodes[vertex_indices[2]] + nodes[vertex_indices[6]]) / 2.0;
//				p9 = (nodes[vertex_indices[3]] + nodes[vertex_indices[7]]) / 2.0;
//				break;
//			case 2:
//				p1 = nodes[vertex_indices[4]];
//				p3 = nodes[vertex_indices[5]];
//				p7 = nodes[vertex_indices[6]];
//				p9 = nodes[vertex_indices[7]];
//				break;
//			}
//
//
//
//
//
//			}*/
//	}
//}
//



