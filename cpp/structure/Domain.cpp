#include "Domain.h"
#include <chrono>
#include <thread>
//#include <iostream>
void el_clear(std::vector<Element>& elements) {
	for (auto e = elements.begin(); e != elements.end(); ++e) {
		std::vector<Element>().swap(elements);
	}
}
void Domain::clear_all_submatrices(std::vector<Element>& elements) {
	for (auto e = elements.begin(); e != elements.end(); ++e) {
		std::vector<std::vector<dcomplex>>().swap(e->cArPAK.matrix);
		std::vector<std::vector<dcomplex>>().swap(e->cBrPAK.matrix);
		std::vector<std::vector<dcomplex>>().swap(e->cBr_HOPS.matrix);

		std::vector<std::vector<double>>().swap(e->fuPowersLagr.matrix);
		std::vector<std::vector<double>>().swap(e->fvPowersLagr.matrix);
		std::vector<std::vector<double>>().swap(e->fwPowersLagr.matrix);

		std::vector<std::vector<double>>().swap(e->uPowers.matrix);
		std::vector<std::vector<double>>().swap(e->vPowers.matrix);
		std::vector<std::vector<double>>().swap(e->wPowers.matrix);

		std::vector<std::vector<double>>().swap(e->fuPowers.matrix);
		std::vector<std::vector<double>>().swap(e->fvPowers.matrix);
		std::vector<std::vector<double>>().swap(e->fwPowers.matrix);

		std::vector<std::vector<double>>().swap(e->fpuPowers.matrix);
		std::vector<std::vector<double>>().swap(e->fpvPowers.matrix);
		std::vector<std::vector<double>>().swap(e->fpwPowers.matrix);

		/*std::vector<std::vector<std::vector<std::vector<double>>>>().swap(e->auMatrix.matrix);
		std::vector<std::vector<std::vector<std::vector<double>>>>().swap(e->avMatrix.matrix);
		std::vector<std::vector<std::vector<std::vector<double>>>>().swap(e->awMatrix.matrix);*/

	}
}
void read_previous_solve(std::string mesh_name, std::vector<std::complex<double>>& cAlphaForLower, std::vector<std::complex<double>>& cAlphaAdjoint, std::vector<std::complex<double>>& cGrhigher) {
	std::ifstream file;
	std::string line;
	file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha.txt"); //lower order cAlpha computation
																				  //getline(file, line);
	std::vector<std::string> line_input;
	while (getline(file, line)) {
		line_input = functions::split(line, ' ');
		cAlphaForLower.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	file.close();
	//
	file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt"); //higher order for cAlphaAdjoint

	while (getline(file, line)) {
		line_input = functions::split(line, ' ');
		cAlphaAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	file.close();
	file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt"); //higher order for cGr from forwars solve

	while (getline(file, line)) {
		line_input = functions::split(line, ' ');
		cGrhigher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1]))); //only works for one freq atm

	}
	file.close();

}

Domain::Domain(std::string mesh_name_in) {
	mesh_name = mesh_name_in;
}
void Domain::makeNT() {
	//bool isAdjoint = false;
	bool isAdjoint = this->sc.useAdjoint;
	//if (isAdjoint) mesh_name = mesh_name + "/adjoint";
	//returns the MatDimDis

	std::vector<std::vector<std::vector<int>>> nTotal(4, std::vector<std::vector<int>>(this->elements.size() + 1, std::vector<int>(5))); //starts at 1
	std::vector<std::vector<int>> ntCum(this->elements.size() + 1, std::vector<int>(4)); //starts at 1

	//int matDimDis = 0;

	int sum = 0;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		int current_element = e->index;
		int sume = 0;
		std::vector<int> order_temp;
		for (int iuvw = 1; iuvw <= 3; ++iuvw) {

			nTotal[iuvw][current_element][1] = e->expansion[0] + 1;
			nTotal[iuvw][current_element][2] = e->expansion[1] + 1;
			nTotal[iuvw][current_element][3] = e->expansion[2] + 1;
			nTotal[iuvw][current_element][iuvw] += -1;

			nTotal[iuvw][current_element][4] = nTotal[iuvw][current_element][1] * nTotal[iuvw][current_element][2] * nTotal[iuvw][current_element][3];
			ntCum[current_element][iuvw] = sum;
			sume += nTotal[iuvw][current_element][4];
			sum += nTotal[iuvw][current_element][4];
		}


		e->disconnected_dimension = sume;
	}
	this->ntCum = ntCum;
	this->nTotal = nTotal;
	this->matDimDis = sum;
	//debugging checks
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/nTotal_mat.txt", nTotal, 1)) {
			std::cout << "ERROR! nTotal does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! nTotal matches FORTRAN!" << std::endl;
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/ntCum_mat.txt", ntCum, 1)) {
			std::cout << "ERROR! ntCum does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! ntCum matches FORTRAN!" << std::endl;
	}

	////matDimDis = sum; //might need to be saved externally, otherwise remove
	////return sum;
}

void Domain::definers_defrmnls() {
	/////testing
	//matrix2d<double> testMat(27, 3);
	//std::vector<double> insert = { 0, 1, 2, 3 };
	//testMat(1) = insert;
	/////testing
	int imax;
	double um, uk, kr, imen, koef;
	int mm, nn, ll, m, n, k, l, broj, jmax, imnl, jmin;
	double ep = 0.0000001;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		int tempOrders = (e->geom_order + 1) * (e->geom_order + 1) * (e->geom_order + 1);
		e->nRs = tempOrders;
		e->rs = matrix2d<double>(tempOrders, 3);
		imax = tempOrders; // note: this assumes that ku, kv, kw are equal to each other 
		for (int icnt = 1; icnt <= imax; ++icnt) {
			k = e->geom_order + 1;
			mm = (icnt - 1) % k;
			ll = (icnt - 1) / k;
			nn = (ll % k);
			k *= k;
			ll = (icnt - 1) / k;
			kr = double(e->geom_order);
			imen = 1.0;
			if (mm == 0) {
				um = -1.0;
				for (k = 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else if (mm == e->geom_order) {
				um = 1;
				for (k = 0; k < e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else {
				um = (2 * mm - kr) / kr;
				for (k = 0; k < mm; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
				for (k = mm + 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			kr = double(e->geom_order);
			if (nn == 0) {
				um = -1.0;
				for (k = 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else if (nn == e->geom_order) {
				um = 1.0;
				for (k = 0; k <= e->geom_order - 1; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else {
				um = (2 * nn - kr) / kr;
				for (k = 0; k < nn; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
				for (k = nn + 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			kr = 1 * e->geom_order;
			if (ll == 0) {
				um = -1.0;
				for (k = 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else if (ll == e->geom_order) {
				um = 1;
				for (k = 0; k < e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			else {
				um = (2 * ll - kr) / kr;
				for (k = 0; k <= ll - 1; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
				for (k = ll + 1; k <= e->geom_order; ++k) {
					uk = (2 * k - kr) / kr;
					imen = imen * (um - uk);
				}
			}
			jmax = pow(2, e->geom_order) - 1;
			for (int jcnt = 0; jcnt <= jmax; ++jcnt) {
				m = 0;
				broj = jcnt;
				koef = 1.0;
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != mm) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++m;
						broj = broj / 2;
					}
				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != nn) {
						uk = (kr - 2 * k) / kr;
						koef = koef * uk;
						//if (abs(koef - 0.0) < ep) koef = 0;
					}
				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != ll) {
						uk = (kr - 2 * k) / kr;
						koef = koef * uk;
						//if (abs(koef - 0.0) < ep) koef = 0;
					}
				}
				imnl = m + 1;

				//Point ptest = (this->nodes)[e->all_indices[icnt]];
				//std::vector<double> temp = this->nodes[e->all_indices[icnt-1]].toVector()*(koef / imen);
				(e->rs)(imnl) = (e->rs)(imnl) + this->nodes[e->all_indices[icnt - 1] - 1].toVector() * (koef / imen);


			}
			jmin = pow(2, e->geom_order);
			jmax = jmin * jmin - 1;
			for (int jcnt = jmin; jcnt <= jmax; ++jcnt) {
				m = 0;
				n = 0;
				broj = jcnt;
				koef = 1.0;
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != mm) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++m;
						broj = broj / 2;
					}
				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != nn) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++n;
						broj = broj / 2;
					}
				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != ll) {
						uk = (kr - 2 * k) / kr;
						koef = koef * uk;
						//if (abs(koef - 0.0) < ep) koef = 0;
					}
				}
				imnl = (e->geom_order + 1) * n + m + 1;


				(e->rs)(imnl) = (e->rs)(imnl) + this->nodes[e->all_indices[icnt - 1] - 1].toVector() * (koef / imen);


			}
			jmin = (pow(2, e->geom_order));
			jmin = jmin * jmin;
			jmax = jmin * pow(2, e->geom_order) - 1;
			for (int jcnt = jmin; jcnt <= jmax; ++jcnt) {
				if (jcnt == 21) {
					int stop = 0;
				}
				m = 0;
				n = 0;
				l = 0;
				broj = jcnt;
				koef = 1;
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != mm) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//	if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++m;
						broj = broj / 2;
					}

				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != nn) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//	if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++n;
						broj = broj / 2;
					}
				}
				kr = e->geom_order;
				for (k = 0; k <= e->geom_order; ++k) {
					if (k != ll) {
						if (broj % 2 == 0) {
							uk = (kr - 2 * k) / kr;
							koef = koef * uk;
							//if (abs(koef - 0.0) < ep) koef = 0;
						}
						else ++l;
						broj = broj / 2;
					}
				}
				imnl = (e->geom_order + 1);
				imnl = imnl * imnl * l + imnl * n + m + 1;

				(e->rs)(imnl) = (e->rs)(imnl) + this->nodes[e->all_indices[icnt - 1] - 1].toVector() * (koef / imen);

			}
		}
	}
	////remove this part below once done coding 
	////********************************
	//std::vector<double> rs_mat;
	//for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
	//	for (int i = 1; i <= e->nRs; ++i) {
	//		for (int j = 1; j < 4; ++j) {
	//			rs_mat.push_back(e->rs(i, j));
	//		}
	//	}
	//}
	//if (check_results) {
	//	if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/rs_mat.txt", rs_mat, 0)) {
	//		std::cout << "ERROR! rs does not match FORTRAN!" << std::endl;
	//	}
	//	else std::cout << "Success! rs matches FORTRAN!" << std::endl;
	//}
	////*********************************
}
int Domain::DADDR(int element_index, int face_index, int coord, int a, int b) {
	//element_index -= 1;
	int addr;
	switch (face_index) {
	case 1:
		switch (coord) {
		case 2:
			addr = this->ntCum[element_index][2] + a * nTotal[2][element_index][3] + b + 1;
			break;
		case 3:
			addr = this->ntCum[element_index][3] + a * nTotal[3][element_index][3] + b + 1;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
			break;
		}
		break;
	case 2:
		switch (coord) {
		case 2:
			addr = ntCum[element_index][2] + nTotal[2][element_index][2] * nTotal[2][element_index][3] + a * nTotal[2][element_index][3] + b + 1;
			break;
		case 3:
			addr = ntCum[element_index][3] + nTotal[3][element_index][2] * nTotal[3][element_index][3] + a * nTotal[3][element_index][3] + b + 1;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
			break;
		}
		break;
	case 3:
		switch (coord) {
		case 1:
			addr = ntCum[element_index][1] + a * nTotal[1][element_index][2] * nTotal[1][element_index][3] + b + 1;
			break;
		case 3:
			addr = ntCum[element_index][3] + a * nTotal[3][element_index][2] * nTotal[3][element_index][3] + b + 1;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
		}
		break;
	case 4:
		switch (coord) {
		case 1:
			addr = ntCum[element_index][1] + a * nTotal[1][element_index][2] * nTotal[1][element_index][3] + nTotal[1][element_index][3] + b + 1;
			break;
		case 3:
			addr = ntCum[element_index][3] + a * nTotal[3][element_index][2] * nTotal[3][element_index][3] + nTotal[3][element_index][3] + b + 1;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
			break;
		}
		break;
	case 5:
		switch (coord) {
		case 1:
			addr = ntCum[element_index][1] + a * nTotal[1][element_index][2] * nTotal[1][element_index][3] + b * nTotal[1][element_index][3] + 1;
			break;
		case 2:
			addr = ntCum[element_index][2] + a * nTotal[2][element_index][2] * nTotal[2][element_index][3] + b * nTotal[2][element_index][3] + 1;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
			break;
		}
		break;
	case 6:
		switch (coord) {
		case 1:
			addr = ntCum[element_index][1] + a * nTotal[1][element_index][2] * nTotal[1][element_index][3] + b * nTotal[1][element_index][3] + 2;
			break;
		case 2:
			addr = ntCum[element_index][2] + a * nTotal[2][element_index][2] * nTotal[2][element_index][3] + b * nTotal[2][element_index][3] + 2;
			break;
		default:
			std::cout << "FAILURE" << std::endl;
			break;
		}

	}
	if (addr == 494) {
		int pause = 0;
	}
	return addr;
}
void Domain::connect_elements() {
	int number_of_elements = this->elements.size();
	std::vector<int> vectorD = std::vector<int>(matDimDis + 1); //adding extra space since it is being indexed starting at 1 instead
	for (int i = 1; i <= matDimDis; ++i) {
		(vectorD)[i] = i;
	}
	//for every shared face, compute the ept and then compute the coef matching
	int facet_tracker = 0;
	for (auto fac_it = this->facets.begin(); fac_it != this->facets.end(); ++fac_it) {
		if (fac_it->element_indices.second <= 0) continue;
		std::vector<std::vector<int>> eqParTab = (*fac_it).MAKEEPT();

		//		//uncomment line below when coeffmatching is finished//		  //
		startMatching(number_of_elements, fac_it->element_indices.first, fac_it->element_indices.second, eqParTab, vectorD);
		++facet_tracker;
		for (int i = 1; i <= matDimDis; ++i) {
			if (vectorD[i] != i) {
				int pause = 0;
			}
		}
	}
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/vectorD_mat.txt", vectorD, 1)) {
			std::cout << "ERROR! vectorD does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! vectorD matches FORTRAN!" << std::endl;
	}
	this->vectorD = vectorD;
}

std::vector<int> Domain::make_bcept(Facet f) {
	std::vector<int> bcept(3);
	//note that from BCCT (implicitly calculated) vertices are always in order of increasing local index
	std::vector<std::vector<int>> fac_table = Tables().face_table[f.local_facet_indices.first - 1];
	//std::vector<std::vector<int>> fac_table_rel = Tables().face_table_rel_ind;
	bcept[0] = std::abs(fac_table[0][4]);
	bcept[1] = std::abs(fac_table[0][5]);
	bcept[2] = f.local_facet_indices.first;
	return bcept;
}
void Domain::bc_coeff_setting(std::vector<int> bcept, Facet f) {
	int indirectAddress;
	std::vector<int> ja = { 0,2,2,1,1,1,1 };
	std::vector<int> jb = { 0, 3, 3, 3, 3, 2, 2 };
	//int local_facet_index = bcept[2];
	int local_facet_index = f.local_facet_indices.first;
	int ja1 = ja[local_facet_index];
	int jb1 = jb[local_facet_index];
	for (int i = 0; i < 2; ++i) {
		int c = bcept[i];
		int amax = nTotal[c][f.element_indices.first][ja1] - 1;
		int bmax = nTotal[c][f.element_indices.first][jb1] - 1;
		for (int a = 0; a <= amax; ++a) {
			for (int b = 0; b <= bmax; ++b) {
				int addr = DADDR(f.element_indices.first, f.local_facet_indices.first, c, a, b);
				if (addr == 205) {
					int stop = 0;
				}
				indirectAddress = std::abs(this->vectorD[addr]); //vector has FORTRAN indexing
				if (indirectAddress != 0) this->vectorD[indirectAddress] = 0;
			}
		}
	}

}
void Domain::set_bc_elements() {
	//for every fac with boundary condition, call makebcept
	//select the case depending on whether -1, 5,  etc.
	for (auto fac_it = this->facets.begin(); fac_it != this->facets.end(); ++fac_it) {
		std::vector<int> bcept = make_bcept(*fac_it);
		if (fac_it->boundary_condition == 0) continue; //default case: domain to domain
		else if (fac_it->boundary_condition == -1) {

			bc_coeff_setting(bcept, *fac_it);
		}
		else if (fac_it->boundary_condition == 5) {
			if (this->elements[fac_it->element_indices.first - 1].type != 1) this->elements[fac_it->element_indices.first - 1].PML_boundary_list[bcept[2] - 1] = 1;
			if (fac_it->element_indices.second > 0)
				if (this->elements[fac_it->element_indices.second - 1].type != 1) this->elements[fac_it->element_indices.second - 1].PML_boundary_list[bcept[2] - 1] = 1;
		}
		else EXIT_FAILURE;
	}
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/vectorD_mat_after_BC.txt", vectorD, 1)) {
			std::cout << "ERROR! vectorD does not match FORTRAN vectorD after BC!" << std::endl;
		}
		else std::cout << "Success! vectorD matches FORTRAN after BC!" << std::endl;
	}
}
std::vector<int> Domain::prenumunknowns() {
	//make a copy of vectorD
	int matDimCon = 0;
	//int negativeOne = -1;
	int sign;
	int vectorDi;
	vectorDUnique = std::vector<int>(matDimDis + 1); //commented because it is undefined (thanks blake)
	auto vectorD_original = vectorD;
	for (int i = 1; i <= matDimDis; ++i) {
		int absVectorD = std::abs(vectorD[i]);
		if (absVectorD != 0) {
			if (absVectorD == i) {
				++matDimCon;
				vectorD[i] = matDimCon;
			}
			else {
				vectorDi = vectorD[i];
				sign = 1;
				if (vectorD[absVectorD] < 0) {
					sign = -1 * sign;
				}
				else if (vectorD[absVectorD] == 0) {
					vectorD[i] = 0;
				}
				if (vectorDi < 0) sign = -1 * sign;
				vectorD[i] = sign * vectorD[absVectorD];
			}
		}
	}
	int matDimCount = 1;
	for (int i = 1; i <= matDimDis; ++i) {
		if (vectorD[i] != 0) {
			if (vectorD[i] == matDimCount) ++matDimCount;
			else {
				vectorDi = vectorD_original[i];
				vectorDUnique[i] = std::abs(vectorDi);
			}
		}
	}
	this->matDimCon = matDimCon;
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/vectorD_mat_prenum.txt", vectorD, 1)) {
			std::cout << "ERROR! vectorD does not match FORTRAN! after prenum" << std::endl;
		}
		else std::cout << "Success! vectorD matches FORTRAN after prenum!" << std::endl;

		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/vectorDunique_mat.txt", vectorDUnique, 1)) {
			std::cout << "ERROR! vectorDUnique does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! VectorDUnique matches FORTRAN!" << std::endl;
	}
	return vectorDUnique;

}

void Domain::fill_eiuvwijk() {
	//NOTE: eiuvwijk goes from 1 to N, NOT 0 to N-1
	std::vector<std::vector<int>> eiUVWijk(matDimDis + 1, std::vector<int>(6));
	int iDis = 0;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		(*e).unknownsStart = iDis + 1;
		for (int iuvw = 1; iuvw <= 3; ++iuvw) {
			switch (iuvw) {
			case 1:
				for (int i = 0; i < e->expansion[0]; ++i) { //nu
					for (int j = 0; j <= e->expansion[1]; ++j) { //nv
						for (int k = 0; k <= e->expansion[2]; ++k) { //nw
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;

						}
					}
				}
				break;
			case 2:
				for (int i = 0; i <= e->expansion[0]; ++i) {
					for (int j = 0; j < e->expansion[1]; ++j) {
						for (int k = 0; k <= e->expansion[2]; ++k) {
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;

						}
					}
				}
				break;
			case 3:
				for (int i = 0; i <= e->expansion[0]; ++i) {
					for (int j = 0; j <= e->expansion[1]; ++j) {
						for (int k = 0; k < e->expansion[2]; ++k) {
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;

						}
					}
				}
				break;
			}
		}


		e->unknownsEnd = iDis;
	}
	this->eiUVWijk = eiUVWijk;
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/eiuvwijk_mat.txt", eiUVWijk, 1)) {
			std::cout << "ERROR! eiUVWijk does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! eiUVWijk matches FORTRAN!" << std::endl;
	}
}
void Domain::fill_eiuvwijk_refine(int order_dif) {
	//NOTE: eiuvwijk goes from 1 to N, NOT 0 to N-1
	std::vector<std::vector<int>> eiUVWijk(matDimDis + 1, std::vector<int>(6));
	int iDis = 0;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		(*e).unknownsStart = iDis + 1;
		for (int iuvw = 1; iuvw <= 3; ++iuvw) {
			switch (iuvw) {
			case 1:
				for (int i = 0; i < e->expansion[0]; ++i) { //nu
					for (int j = 0; j <= e->expansion[1]; ++j) { //nv
						for (int k = 0; k <= e->expansion[2]; ++k) { //nw
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;
							if (e->refine) {
								int vectorD_storage = this->vectorD[iDis];
								if (i >= e->expansion[0] - order_dif) {
									if (j < 2 || k < 2) this->vectorD[iDis] = 0;

								}
								if (j > e->expansion[1] - order_dif) {
									if (k < 2) this->vectorD[iDis] = 0;

								}
								if (k > e->expansion[2] - order_dif) {
									if (j < 2) this->vectorD[iDis] = 0;
								}
								//adding the tangential components 
								//check if v or w connection exists
								//if (i > 1) continue;
								//if (e->adjacent_list.v_dir > 0) {
								//	if (k < 2) continue;
								//	if (e->adjacent_list.v_dir < 2) { //add -v face fcn
								//		if (j != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.v_dir < 3) {//add +v face fcn
								//		if (j != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}if (e->adjacent_list.w_dir > 0) {
								//	if (j < 2) continue;
								//	if (e->adjacent_list.w_dir < 2) {//add -w face fcn
								//		if (k != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.w_dir < 3) { //add +w face fcn
								//		if (k != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}
							}
						}
					}
				}
				break;
			case 2:
				for (int i = 0; i <= e->expansion[0]; ++i) {
					for (int j = 0; j < e->expansion[1]; ++j) {
						for (int k = 0; k <= e->expansion[2]; ++k) {
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;
							if (e->refine) {
								int vectorD_storage = this->vectorD[iDis];
								if (j >= e->expansion[0] - order_dif) {
									if (i < 2 || k < 2) this->vectorD[iDis] = 0;
									continue;
								}
								if (i > e->expansion[0] - order_dif) {
									if (k < 2) this->vectorD[iDis] = 0;
									continue;
								}
								if (k > e->expansion[2] - order_dif) {
									if (i < 2) this->vectorD[iDis] = 0;
								}
								//adding the tangential components 
								//check if u or w connection exists
								//if (j > 1) continue;
								//if (e->adjacent_list.u_dir > 0) {
								//	if (k < 2) continue;
								//	if (e->adjacent_list.u_dir < 2) { //add -v face fcn
								//		if (i != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.u_dir < 3) {//add +v face fcn
								//		if (i != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}if (e->adjacent_list.w_dir > 0) {
								//	if (i < 2) continue;
								//	if (e->adjacent_list.w_dir < 2) {//add -w face fcn
								//		if (k != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.w_dir < 3) { //add +w face fcn
								//		if (k != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}
							}
						}
					}
				}
				break;
			case 3:
				for (int i = 0; i <= e->expansion[0]; ++i) {
					for (int j = 0; j <= e->expansion[1]; ++j) {
						for (int k = 0; k < e->expansion[2]; ++k) {
							++iDis;
							eiUVWijk[iDis][1] = e->index;
							eiUVWijk[iDis][2] = iuvw;
							eiUVWijk[iDis][3] = i;
							eiUVWijk[iDis][4] = j;
							eiUVWijk[iDis][5] = k;
							if (e->refine) {
								int vectorD_storage = this->vectorD[iDis];
								if (k >= e->expansion[2] - order_dif) {
									if (i < 2 || j < 2) this->vectorD[iDis] = 0;
									continue;
								}
								if (j > e->expansion[1] - order_dif) {
									if (i < 2) this->vectorD[iDis] = 0;
									continue;
								}
								if (i > e->expansion[0] - order_dif) {
									if (j < 2) this->vectorD[iDis] = 0;
								}
								//adding the tangential components 
								//check if u or v connection exists
								//if (k > 1) continue;
								//if (e->adjacent_list.u_dir > 0) {
								//	if (j < 2) continue;
								//	if (e->adjacent_list.u_dir < 2) { //add -v face fcn
								//		if (i != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.u_dir < 3) {//add +v face fcn
								//		if (i != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}if (e->adjacent_list.v_dir > 0) {
								//	if (i < 2) continue;
								//	if (e->adjacent_list.v_dir < 2) {//add -w face fcn
								//		if (j != 0) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}if (e->adjacent_list.v_dir < 3) { //add +w face fcn
								//		if (j != 1) continue;
								//		this->vectorD[iDis] = vectorD_storage;
								//	}
								//}
							}
						}
					}
				}
				break;
			}
		}
		e->unknownsEnd = iDis;
	}
	this->eiUVWijk = eiUVWijk;
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/eiuvwijk_mat.txt", eiUVWijk, 1)) {
			std::cout << "ERROR! eiUVWijk does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! eiUVWijk matches FORTRAN!" << std::endl;
	}
}

void Domain::make_unknown_description() {
	std::vector<std::vector<int>> unkn_desc(this->matDimCon + 1, std::vector<int>(2));
	int el;
	int j;
	int index;
	for (int i = 1; i <= matDimDis; ++i) {
		if (vectorD[i] != 0) {
			el = eiUVWijk[i][1];
			j = std::abs(vectorD[i]);
			unkn_desc[j][1] = unkn_desc[j][1] + 1;
			//index = 2 * unkn_desc[j][1];
			//unkn_desc[j][index] = el;
			//unkn_desc[j][index + 1] = i;
			unkn_desc[j].push_back(el);
			unkn_desc[j].push_back(i);
		}
	}
	this->unknown_description = unkn_desc;
	if (check_results) {
		if (!functions::_debug_matrix("../exampleFiles/debug_data/" + mesh_name + "/unknowns_desc_mat.txt", unknown_description, 1)) {
			std::cout << "ERROR! unknown_description does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! unknown_description matches FORTRAN!" << std::endl;
	}
}

void Domain::vectorDSignPW(int addr1, int addr2, int signPW, int b, std::vector<int>& vectorD) {

	if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
		functions::MIMAISWAP(&addr1, &addr2);
	}

	if (signPW == -1 && b % 2 == 0) {
		vectorD[addr2] = -vectorD[addr1];
	}
	else {
		vectorD[addr2] = vectorD[addr1];
	}
}

void Domain::startMatching(int numberOfElements, int e1, int e2, std::vector<std::vector<int>>& ept, std::vector<int>& vectorD) {
	//e1 -= 1; //nTotal etc is index starting at 0
	//e2 -= 1;

	/* Declare Variables Section */
	if (e2 == 4 || e1 == 4) {
		int stop = 0;
	}
	std::vector<int> ja = { 0,2,2,1,1,1,1 };
	std::vector<int> jb = { 0,3,3,3,3,2,2 };



	int a1, a2, b1, b2, x1, x2, y1, y2, addr1, addr2;
	int alo, ahi, blo, bhi, f1, c1, f2, c2, c2mod;
	int znak1, znak2, ja1, jb1, ja2, jb2;

	/* End Declare Variables Section */

	/* We are looping through the first two columns of ept*/
	for (int ic = 1; ic <= 2; ++ic) {
		int tempzero = 0;
		int tempone = 1;

		int signRT = 1;
		int signPW = 1;
		int pw = ept[4][ic]; //power values of the ic column; possible values = {1,2}
		int rt = ept[5][ic]; //rooftop values of the ic column; possible values = {1,2}
		int sw = ept[3][ic]; //swap values of the ic column; possible values = {0,1}


		if (ept[2][rt] < 0) {
			signRT = -1;
		}
		else {
			signRT = 1;
		}

		if (ept[2][pw] < 0) {
			signPW = -1;
		}
		else {
			signPW = 1;
		}

		f1 = ept[1][3];
		c1 = ept[1][ic];
		f2 = ept[2][3];
		c2 = ept[2][ic];

		c2mod = std::abs(c2);

		ja1 = ja[f1];
		jb1 = jb[f1];

		if (sw == 1) {
			ja2 = jb[f2];
			jb2 = ja[f2];
		}
		else {
			ja2 = ja[f2];
			jb2 = jb[f2];
		}

		a1 = nTotal[c1][e1][ja1] - 1;
		a2 = nTotal[c2mod][e2][ja2] - 1;
		b1 = nTotal[c1][e1][jb1] - 1;
		b2 = nTotal[c2mod][e2][jb2] - 1;



		if (a2 < a1) {
			alo = a2;
			ahi = a1;
		}
		else {
			alo = a1;
			ahi = a2;
		}



		if (b2 < b1) {
			blo = b2;
			bhi = b1;
		}
		else {
			blo = b1;
			bhi = b2;
		}


		switch (rt) {
		case 1:
			for (int b = 0; b <= blo; ++b) {
				switch (sw) {
				case 1:
					switch (signRT) {
					case 1:
						x1 = b;
						x2 = b;
						y1 = 0;
						y2 = 1;
						break;
					case -1:
						x1 = b;
						x2 = b;
						y1 = 1;
						y2 = 0;
						break;
					}
					break;
				case 0:
					switch (signRT) {
					case 1:
						x1 = 0;
						x2 = 1;
						y1 = b;
						y2 = b;
						break;
					case -1:
						x1 = 1;
						x2 = 0;
						y1 = b;
						y2 = b;
						break;
					}
					break;
				}

				addr2 = DADDR(e2, f2, c2mod, x1, y1);
				addr1 = DADDR(e1, f1, c1, tempzero, b);
				vectorDSignPW(addr1, addr2, signPW, b, vectorD);
				addr2 = DADDR(e2, f2, c2mod, x2, y2);
				addr1 = DADDR(e1, f1, c1, tempone, b);
				vectorDSignPW(addr1, addr2, signPW, b, vectorD);
			}

			for (int b = blo + 1; b <= bhi; ++b) {
				if (b2 < b1) {
					addr1 = DADDR(e1, f1, c1, tempzero, b);
					vectorD[addr1] = 0;
					addr1 = DADDR(e1, f1, c1, tempone, b);
					vectorD[addr1] = 0;
				}
				else {
					switch (sw) {
					case 1:
						x1 = b;
						x2 = b;
						y1 = 0;
						y2 = 1;
						break;
					case 0:
						x1 = 0;
						x2 = 1;
						y1 = b;
						y2 = b;
						break;
					}
					addr2 = DADDR(e2, f2, c2mod, x1, y1);
					vectorD[addr2] = 0;
					addr2 = DADDR(e2, f2, c2mod, x2, y2);
					vectorD[addr2] = 0;
				}
			}

			for (int a = 2; a <= alo; ++a) {
				for (int b = 0; b <= blo; ++b) {
					switch (sw) {
					case 1:
						addr2 = DADDR(e2, f2, c2mod, b, a);
						addr1 = DADDR(e1, f1, c1, a, b);

						if (signPW == -1 && b % 2 == 0) {
							znak1 = -1;
						}
						else {
							znak1 = 1;
						}
						if (signRT == -1 && a % 2 != 0) {
							znak2 = -1;
						}
						else {
							znak2 = 1;
						}

						if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
							functions::MIMAISWAP(&addr1, &addr2);
						}
						vectorD[addr2] = znak1 * znak2 * vectorD[addr1];
						break;
					case 0:
						addr2 = DADDR(e2, f2, c2mod, a, b);
						addr1 = DADDR(e1, f1, c1, a, b);

						if (signPW == -1 && b % 2 == 0) {
							znak1 = -1;
						}
						else {
							znak1 = 1;
						}
						if (signRT == -1 && a % 2 != 0) {
							znak2 = -1;
						}
						else {
							znak2 = 1;
						}
						if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
							functions::MIMAISWAP(&addr1, &addr2);
						}
						vectorD[addr2] = znak1 * znak2 * vectorD[addr1];
						break;
					}
				}
				for (int b = blo + 1; b <= bhi; ++b) {
					if (b2 < b1) {
						addr1 = DADDR(e1, f1, c1, a, b);
						vectorD[addr1] = 0;
					}
					else {
						switch (sw) {
						case 1:
							addr2 = DADDR(e2, f2, c2mod, b, a);
							vectorD[addr2] = 0;
							break;
						case 0:
							addr2 = DADDR(e2, f2, c2mod, a, b);
							vectorD[addr2] = 0;
							break;
						}
					}
				}
			}

			for (int a = alo + 1; a <= ahi; ++a) {
				if (a2 < a1) {
					for (int b = 0; b <= b1; ++b) {
						addr1 = DADDR(e1, f1, c1, a, b);
						vectorD[addr1] = 0;
					}
				}
				else {
					for (int b = 0; b <= b2; ++b) {
						switch (sw) {
						case 1:
							addr2 = DADDR(e2, f2, c2mod, b, a);
							vectorD[addr2] = 0;
							break;
						case 0:
							addr2 = DADDR(e2, f2, c2mod, a, b);
							vectorD[addr2] = 0;
						}
					}
				}
			}
			break; //case 1

		case 2:
			for (int a = 0; a <= alo; ++a) {
				switch (sw) {
				case 1:
					switch (signRT) {
					case 1:
						x1 = 0;
						x2 = 1;
						y1 = a;
						y2 = a;
						break;
					case -1:
						x1 = 1;
						x2 = 0;
						y1 = a;
						y2 = a;
						break;
					}
					break;
				case 0:
					switch (signRT) {
					case 1:
						x1 = a;
						x2 = a;
						y1 = 0;
						y2 = 1;
						break;
					case -1:
						x1 = a;
						x2 = a;
						y1 = 1;
						y2 = 0;
						break;
					}
					break;
				}


				addr2 = DADDR(e2, f2, c2mod, x1, y1);
				addr1 = DADDR(e1, f1, c1, a, tempzero);

				if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
					functions::MIMAISWAP(&addr1, &addr2);
				}

				if (signPW == -1 && (a % 2) == 0) {
					vectorD[addr2] = -vectorD[addr1];
				}
				else {
					vectorD[addr2] = vectorD[addr1];
				}

				addr2 = DADDR(e2, f2, c2mod, x2, y2);
				addr1 = DADDR(e1, f1, c1, a, tempone);

				if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
					functions::MIMAISWAP(&addr1, &addr2);
				}

				if (signPW == -1 && a % 2 == 0) {
					vectorD[addr2] = -vectorD[addr1];
				}
				else {
					vectorD[addr2] = vectorD[addr1];
				}
			}


			for (int a = alo + 1; a <= ahi; ++a) {
				if (a2 < a1) {
					addr1 = DADDR(e1, f1, c1, a, tempzero);
					vectorD[addr1] = 0;
					addr1 = DADDR(e1, f1, c1, a, tempone);
					vectorD[addr1] = 0;
				}
				else {
					switch (sw) {
					case 1:
						x1 = 0;
						x2 = 1;
						y1 = a;
						y2 = a;
						break;
					case 0:
						x1 = a;
						x2 = a;
						y1 = 0;
						y2 = 1;
						break;
					}
					addr2 = DADDR(e2, f2, c2mod, x1, y1);
					vectorD[addr2] = 0;
					addr2 = DADDR(e2, f2, c2mod, x2, y2);
					vectorD[addr2] = 0;
				}
			}

			for (int b = 2; b <= blo; ++b) {
				for (int a = 0; a <= alo; ++a) {
					switch (sw) {
					case 1:
						addr2 = DADDR(e2, f2, c2mod, b, a);
						addr1 = DADDR(e1, f1, c1, a, b);
						if (signPW == -1 && a % 2 == 0) {
							znak1 = -1;
						}
						else {
							znak1 = 1;
						}
						if (signRT == -1 && b % 2 != 0) {
							znak2 = -1;
						}
						else {
							znak2 = 1;
						}

						if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
							functions::MIMAISWAP(&addr1, &addr2);
						}

						vectorD[addr2] = znak1 * znak2 * vectorD[addr1];
						break;
					case 0:
						addr2 = DADDR(e2, f2, c2mod, a, b);
						addr1 = DADDR(e1, f1, c1, a, b);
						if (signPW == -1 && a % 2 == 0) {
							znak1 = -1;
						}
						else {
							znak1 = 1;
						}
						if (signRT == -1 && b % 2 != 0) {
							znak2 = -1;
						}
						else {
							znak2 = 1;
						}


						if (std::abs(vectorD[addr2]) < std::abs(vectorD[addr1])) {
							functions::MIMAISWAP(&addr1, &addr2);
						}


						vectorD[addr2] = znak1 * znak2 * vectorD[addr1];
						break;
					}
				}

				for (int a = alo + 1; a <= ahi; ++a) {
					if (a2 < a1) {
						addr1 = DADDR(e1, f1, c1, a, b);
						vectorD[addr1] = 0;
					}
					else {
						switch (sw) {
						case 1:
							addr2 = DADDR(e2, f2, c2mod, b, a);
							vectorD[addr2] = 0;
							break;
						case 0:
							addr2 = DADDR(e2, f2, c2mod, a, b);
							vectorD[addr2] = 0;
							break;
						}
					}
				}
			}

			for (int b = blo + 1; b <= bhi; ++b) {
				if (b2 < b1) {
					for (int a = 0; a <= a1; ++a) {
						addr1 = DADDR(e1, f1, c1, a, b);
						vectorD[addr1] = 0;
					}
				}
				else {
					for (int a = 0; a <= a2; a++) {
						switch (sw) {
						case 1:
							addr2 = DADDR(e2, f2, c2mod, b, a);
							vectorD[addr2] = 0;
							break;
						case 0:
							addr2 = DADDR(e2, f2, c2mod, a, b);
							vectorD[addr2] = 0;
							break;
						}
					}
				}
			}
			break; //case 2


		default:
			std::cout << "Error in matching coefficients!\n";

			break;
		}
	}
}


//powers are getting moved to BasisEval
//skipping legendre powers for now since it is not used in the example we have

//void findPowersLegendre(int nu, int nv, int nw, int nglu, int nglv, int glw, std::vector<double> xglu, std::vector<double> xglv, std::vector<double> xglw,
//	matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
//	matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
//	matrix2d<double> legendre; //starts at 0
//	legendre(0, 0) = 1; legendre(1, 1) = 1;
//	legendre(2, 0) = -.5; legendre(2, 2) = 1.5;
//	legendre(3, 1) = -1.5;
//}

void Domain::findPowersLagrange(const int& kuvw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
	matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr) {
	//right now we have kuvw -- for now the order is the same in each direction -- this could be generalized later

	double x, xj, xm, f;

	for (int i = 0; i <= nglu + 1; i++) { // for u
		x = xglu[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0 * m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0 * j - kuvw) / kuvw;
				if (j != m) {
					f = f * (x - xj) / (xm - xj);
				}
			}//for j
			fuPowersLagr(i, m) = f;
		}//for m
	}//for i

	for (int i = 0; i <= nglv + 1; i++) { // for v
		x = xglv[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0 * m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0 * j - kuvw) / kuvw;
				if (j != m) {
					f = f * (x - xj) / (xm - xj);
				}
			}//for j
			fvPowersLagr(i, m) = f;
		}//for m
	}//for i

	for (int i = 0; i <= nglw + 1; i++) { // for w
		x = xglw[i]; //Gauss-Legendre Integration point
		for (int m = 0; m <= kuvw; m++) {
			xm = (2.0 * m - kuvw) / kuvw; //Lagrange Interpolation point
			f = 1.0;
			for (int j = 0; j <= kuvw; j++) {
				xj = (2.0 * j - kuvw) / kuvw;
				if (j != m) {
					f = f * (x - xj) / (xm - xj);
				}
			}//for j
			fwPowersLagr(i, m) = f;
		}//for m
	}//for i

}//void findPowersLagrange

double Domain::trilagrangeoduvw(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr) {
	return fuPowersLagr(m, i) * fvPowersLagr(n, j) * fwPowersLagr(l, k);
}


void Domain::find_eps_mu_matrix(const int& kuvw, const int& nglu, const int& nglv, const int& nglw, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr,
	const int& size1, const std::vector<std::vector<std::complex<double>>>& muRel, matrix4d<std::complex<double>>& muRelInt) {
	//std::cout << "entered eps_mu_matrix" << std::endl;
	std::complex<double> cMu;
	for (auto s = 1; s <= size1; ++s) { //matrix entries
		for (int m = 0; m <= nglu + 1; ++m) {
			for (int n = 0; n <= nglv + 1; ++n) {
				for (int l = 0; l <= nglw + 1; ++l) {
					cMu = std::complex<double>(0, 0);
					for (int mat_loc = 0; mat_loc < muRel.size(); ++mat_loc) { //nodes
						int kuvwp = kuvw + 1;
						double f = trilagrangeoduvw(m, n, l, mat_loc % (kuvwp), int(floor(mat_loc / kuvwp)) % kuvwp,
							int(floor(mat_loc / (kuvwp * kuvwp))), fuPowersLagr, fvPowersLagr, fwPowersLagr);

						muRelInt(s, m, n, l) += muRel[mat_loc][s - 1] * f;
					}
					//	muRelInt(s, m, n, l) = cMu;
				}
			}

		}

	}
	//std::cout << "exited eps_mu_match" << std::endl;
}

void Domain::find_eps_mu_hermitian(const int& size1, const int& nglu, const int& nglv, const int& nglw, matrix4d < std::complex<double>>& MuRelIntInv) {
	for (int u = 0; u <= nglu + 1; ++u) {
		for (int v = 0; v <= nglv + 1; ++v) {
			for (int w = 0; w <= nglw + 1; ++w) {
				if (size1 == 1) MuRelIntInv(1, u, v, w).imag(-MuRelIntInv(1, u, v, w).imag());
				else if (size1 == 3) {
					MuRelIntInv(1, u, v, w).imag(-MuRelIntInv(1, u, v, w).imag());
					MuRelIntInv(2, u, v, w).imag(-MuRelIntInv(2, u, v, w).imag());
					MuRelIntInv(3, u, v, w).imag(-MuRelIntInv(3, u, v, w).imag());
				}
				else if (size1 == 6) {
					MuRelIntInv(1, u, v, w).imag(-MuRelIntInv(1, u, v, w).imag());
					MuRelIntInv(2, u, v, w).imag(-MuRelIntInv(2, u, v, w).imag());
					MuRelIntInv(3, u, v, w).imag(-MuRelIntInv(3, u, v, w).imag());
					MuRelIntInv(4, u, v, w).imag(-MuRelIntInv(4, u, v, w).imag());
					MuRelIntInv(5, u, v, w).imag(-MuRelIntInv(5, u, v, w).imag());
					MuRelIntInv(6, u, v, w).imag(-MuRelIntInv(6, u, v, w).imag());
				}
				else {
					double t1_im = -MuRelIntInv(4, u, v, w).imag();
					double t2_im = -MuRelIntInv(7, u, v, w).imag();
					double t3_im = -MuRelIntInv(8, u, v, w).imag();
					MuRelIntInv(1, u, v, w).imag(-MuRelIntInv(1, u, v, w).imag());
					MuRelIntInv(4, u, v, w).imag(-MuRelIntInv(2, u, v, w).imag());
					MuRelIntInv(5, u, v, w).imag(-MuRelIntInv(5, u, v, w).imag());
					MuRelIntInv(7, u, v, w).imag(-MuRelIntInv(3, u, v, w).imag());
					MuRelIntInv(8, u, v, w).imag(-MuRelIntInv(6, u, v, w).imag());
					MuRelIntInv(9, u, v, w).imag(-MuRelIntInv(9, u, v, w).imag());
					MuRelIntInv(2, u, v, w).imag(t1_im);
					MuRelIntInv(3, u, v, w).imag(t2_im);
					MuRelIntInv(6, u, v, w).imag(t3_im);

				}
			}
		}
	}
}

void Domain::SPARSE_PACKING(std::vector<int>& iRow,
	std::vector<int>& jCol, int& nNZ, int noAbcBCs,
	std::vector<dcomplex>& cArSparse,
	std::vector<dcomplex>& cBrSparse,
	std::vector<dcomplex>& cSrSparse) {
	cArSparse.push_back(0);
	cBrSparse.push_back(0);
	cSrSparse.push_back(0);
	jCol.push_back(0);
	iRow.push_back(0);
	//noAbcBCs should be 0

	//std::vector<Element>* elements = &(dom.elements);
	int myElementCount = this->elements.size();
	std::vector<int> IsContained = std::vector<int>(myElementCount + 1);



	int e1, e2, i0, j0, i1, j1, lim1, D1, D2, C1, C2;
	int index1, index2, k0, k1, k2, k3;
	int FoundFlag;



	for (auto elem = ++this->elements.begin(); elem != this->elements.end(); ++elem) {


		lim1 = elem->unknownsStart - 1;
		for (int i = elem->unknownsStart; i <= elem->unknownsEnd; ++i) {
			i0 = i - lim1;
			D1 = this->vectorDUnique[i];
			//C1 = The index of where the basis unknown goes globally
			C1 = std::abs(this->vectorD[i]);
			//Skip if this unknown is zero
			if (D1 == 0) {
				continue;
			}

			for (int j = elem->unknownsStart; j <= elem->unknownsEnd; ++j) {
				//How far from the beginning element we are for the testing function
				j0 = j - lim1;
				D2 = this->vectorDUnique[j];
				//C2 = the index of where the testing unknown goes globally
				C2 = std::abs(this->vectorD[j]);
				if (D2 != 0) {
					e2 = 0;

					//Index 1 and Index 2 tell us if the unknowns appear in more than one element
					index1 = this->unknown_description[C1][1];
					index2 = this->unknown_description[C2][1];

					//If both unknowns (basis and testing) appear in the global system at some point
					if ((index1 >= 1) && (index2 >= 1)) {

						//Whether or not we have found this entry on a different element
						FoundFlag = 0;
						//For all elements in which either unknown appears
						for (int k0 = 1; k0 <= index1; ++k0) {
							for (int k1 = 1; k1 <= index2; ++k1) {
								//If both unknowns appear in a previous element, which is owned by this process,
								//  those entries should be copied over to the previous element
								//std::cout << "C1: " << C1 << " C2: " << C2 << " 2*k0: " << 2*k0 << " 2*k2: " << 2*k2 << std::endl;
								if ((this->unknown_description[C1][2 * k0] == this->unknown_description[C2][2 * k1]) &&
									(this->unknown_description[C1][2 * k0] < elem->index) && (FoundFlag == 0)) //&&
									/*(IsContained[unknown_description(C1, 2 * k0)] == 1))*/ {
									//The element it should be copied to
									e2 = this->unknown_description[C1][2 * k0];
									//The global disconnected unknown on that element (basis)
									k2 = k0;
									//The global disconnected unknown on that element (testing)
									k3 = k1;
									FoundFlag = 1;
								}
							}
						}
					}


					//If we found somewhere else to put those entries, copy them over to the previous element
					if (e2 != 0) {

						i1 = this->unknown_description[C1][2 * k2 + 1] - this->elements[e2 - 1].unknownsStart + 1;
						j1 = this->unknown_description[C2][2 * k3 + 1] - this->elements[e2 - 1].unknownsStart + 1;
						//std::cout << "e2: " << e2 << std::endl;
						//std::cout << "i1: " << i1 << std::endl;
						//std::cout << "j1: " << j1 << std::endl;
						//this->elements[e2-1].error_structure.push_back(Error_struct(e2, e1, i1, j1, this->elements[e2 - 1].cArPAK(i1, j1), this->elements[e2 - 1].cBrPAK(i1, j1), 0));
						this->elements[e2 - 1].cArPAK(i1, j1) =
							this->elements[e2 - 1].cArPAK(i1, j1) + elem->cArPAK(i0, j0);
						this->elements[e2 - 1].cBrPAK(i1, j1) =
							this->elements[e2 - 1].cBrPAK(i1, j1) + elem->cBrPAK(i0, j0);
						if (noAbcBCs > 0) {
							this->elements[e2 - 1].cSrPAK(i1, j1) =
								this->elements[e2 - 1].cSrPAK(i1, j1) + elem->cSrPAK(i0, j0);
						}


						elem->cArPAK(i0, j0) = 0.0;
						elem->cBrPAK(i0, j0) = 0.0;
						if (noAbcBCs > 0) {
							elem->cSrPAK(i0, j0) = { 0, 0 };
						}
					}
				}
			}
		}
	}



	for (auto elem = this->elements.begin(); elem != this->elements.end(); ++elem) {
		//int error_tracker = 0; //keeps track of index for error_struct

		lim1 = elem->unknownsStart - 1;

		for (int i = elem->unknownsStart; i <= elem->unknownsEnd; ++i) {
			i0 = i - lim1;
			D1 = abs(this->vectorD[i]);
			if (D1 == 0) { continue; }

			for (int j = i; j <= elem->unknownsEnd; ++j) {
				j0 = j - lim1;
				D2 = abs(this->vectorD[j]);
				if (D2 != 0) {

					dcomplex czero = { 0,0 };
					if ((elem->cArPAK(i0, j0) != czero) ||
						(elem->cBrPAK(i0, j0) != czero)) {

						if (D1 <= D2) {

							nNZ = nNZ + 1;
							//cArSparse[nNZ] = elem->cArPAK(i0, j0);
							//cBrSparse[nNZ] = elem->cBrPAK(i0, j0);
							//elem->error_structure[elem->index - 1].nNz = nNZ;
						/*	for (auto err_it = elem->error_structure.begin(); err_it != elem->error_structure.end(); ++err_it) {
								if (i0 == err_it->i0 && j0 == err_it->j0) {
										err_it->nNz = nNZ;
										break;
								}
							}*/
							cArSparse.push_back(elem->cArPAK(i0, j0));
							cBrSparse.push_back(elem->cBrPAK(i0, j0));
							if (noAbcBCs > 0) { //cSrSparse[nNZ] = elem->cSrPAK(i0, j0);
								cSrSparse.push_back(elem->cSrPAK(i0, j0));
							}

							//iRow[nNZ] = D1;
							//jCol[nNZ] = D2;
							iRow.push_back(D1);
							jCol.push_back(D2);

							if (D1 != D2) {
								nNZ = nNZ + 1;
								//cArSparse[nNZ] = elem->cArPAK(j0, i0);
								//cBrSparse[nNZ] = elem->cBrPAK(j0, i0);
								cArSparse.push_back(elem->cArPAK(j0, i0));
								cBrSparse.push_back(elem->cBrPAK(j0, i0));
								if (noAbcBCs > 0) { //cSrSparse[nNZ] = elem->cSrPAK(j0, i0); 
									cSrSparse.push_back(elem->cSrPAK(j0, i0));
								}

								//iRow[nNZ] = D2;
							//	jCol[nNZ] = D1;
								iRow.push_back(D2);
								jCol.push_back(D1);
							}
						}
						else {

							nNZ = nNZ + 1;
							//	cArSparse[nNZ] = elem->cArPAK(i0, j0);
							//	cBrSparse[nNZ] = elem->cBrPAK(i0, j0);
								/*for (auto err_it = elem->error_structure.begin(); err_it != elem->error_structure.end(); ++err_it) {
									if (i0 == err_it->i0 && j0 == err_it->j0) {
										err_it->nNz = nNZ;
										break;
									}
								}*/
							cArSparse.push_back(elem->cArPAK(i0, j0));
							cBrSparse.push_back(elem->cBrPAK(i0, j0));
							if (noAbcBCs > 0) { //cSrSparse[nNZ] = elem->cSrPAK(i0, j0); 
								cSrSparse.push_back(elem->cSrPAK(i0, j0));
							}

							//iRow[nNZ] = D2;
							//jCol[nNZ] = D1;

							iRow.push_back(D2);
							jCol.push_back(D1);

							if (D1 != D2) {
								nNZ = nNZ + 1;
								//cArSparse[nNZ] = elem->cArPAK(j0, i0);
								//cBrSparse[nNZ] = elem->cBrPAK(j0, i0);
								cArSparse.push_back(elem->cArPAK(j0, i0));
								cBrSparse.push_back(elem->cBrPAK(j0, i0));
								if (noAbcBCs > 0) //cSrSparse[nNZ] = elem->cSrPAK(j0, i0);
									cSrSparse.push_back(elem->cSrPAK(j0, i0));

								//iRow[nNZ] = D1;
								//jCol[nNZ] = D2;
								iRow.push_back(D1);
								jCol.push_back(D2);
							}
						}
					}
				}
			}
		}
	}




}





void Domain::findelements_abs_sparse() {
	std::cout << "Entered matrix filling!" << std::endl;
	int dummy = 0;
	int nNz = 0;
	//bool isAdjoint = false;
	bool isAdjoint = this->sc.useAdjoint;
	//if (isAdjoint) mesh_name = mesh_name + "/adjoint";
	int tracker = 0;
	double k0 = this->scatter1.K0[1];
	int pause = 0;
	//auto t1 = std::chrono::high_resolution_clock::now();
	//i0_table.resize(this->elements.size());

	//temp to test basiseval
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	//end temp to test basis eval

	//for (auto e = elements.begin(); e != elements.end(); ++e) {
#pragma omp parallel for num_threads(6)
	for (int e_index = 0; e_index < elements.size(); ++e_index) {
		//	std::cout << e_index << std::endl;
		auto e = &elements[e_index];
		e->cGr_eps_el = std::vector<std::complex<double>>(this->scatter1.cGr[1].size(), 0.0);
		//(i0_table[e->index - 1]) = std::vector<std::vector<std::vector<std::vector<int>>>>(4, std::vector<std::vector<std::vector<int>>>(e->expansion[0]+1, std::vector<std::vector<int>>(e->expansion[1]+1, std::vector<int>(e->expansion[2]+1))));
		//auto t1 = std::chrono::high_resolution_clock::now();
		Integral_c c_integrator(&(*e));
		Integral_d d_integrator(&(*e));
		Integral_g gInts(&(*e));

		e->cArPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		e->cBrPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		e->cSrPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		e->cBr_HOPS = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);


		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1; //kolund basis

		e->uPowers = eval.get_u_samples(nglu, nu)[0];
		e->vPowers = eval.get_v_samples(nglv, nv)[0];
		e->wPowers = eval.get_w_samples(nglw, nw)[0];
		e->fuPowers = eval.get_u_samples(nglu, nu)[1];
		e->fvPowers = eval.get_v_samples(nglv, nv)[1];
		e->fwPowers = eval.get_w_samples(nglw, nw)[1];
		e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
		e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
		e->fpwPowers = eval.get_w_samples(nglw, nw)[2];

		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		//if indicator = 1, then kolun, if indicator = 2, then legendre
		//if (indicatorBasisType == 1)
			//kolund basis functions
//			finduPowers(nu, nglu, xglu, e->uPowers, e->fuPowers, e->fpuPowers);
//			findvPowers(nv, nglv, xglv, e->vPowers, e->fvPowers, e->fpvPowers);
//			findwPowers(nw, nglw, xglw, e->wPowers, e->fwPowers, e->fpwPowers);


		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}


			}
		}
		int iCon;
		//#pragma omp parallel for
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {

			iCon = abs(vectorD[kMat]);
			//if (iCon == 0) continue;
			int eh = eiUVWijk[kMat][1];
			int iuvwh = eiUVWijk[kMat][2];
			int ih = eiUVWijk[kMat][3];
			int jh = eiUVWijk[kMat][4];
			int kh = eiUVWijk[kMat][5];
			if (iuvwh == 1 && ih == 3 && jh == 0 && kh == 1) {
				int pause = 0;
			}
			//if (iuvwh == 1 && ih == 0 && jh == 4 && kh == 4) {
			//	int pause = 0;
			//}
			if (iCon == 0) continue;
			int i0, j0;
			if (isAdjoint) {
				for (int face = 1; face <= 6; ++face) {
					if (((this->elements[eh - 1]).abcList[face - 1] != 1) &&
						this->elements[eh - 1].PML_boundary_list[face - 1] != 1) {
						continue;
					}
					if (((this->elements[eh - 1]).PML_boundary_list[face - 1] == 1) && this->elements[eh - 1].materials.region == 2) {
						continue;
					}
					gInts.set_e(&(this->elements[eh - 1]));
					gInts.findGWave(&scatter1, vectorD, isAdjoint, iuvwh, ih, jh, kh, size1, kMat, iCon, face);
				}
			}
			else {
				int face = 1;
				if (e->materials.region == 1) {
					//should be for element eh
					gInts.set_e(&(this->elements[eh - 1]));
					gInts.findGWave(&scatter1, vectorD, isAdjoint, iuvwh, ih, jh, kh, size1,//. . .
						kMat, iCon, face);

				}
			}

			for (int lMat = kMat; lMat <= e->unknownsEnd; ++lMat) {

				int jCon = abs(vectorD[lMat]);
				if (jCon == 0) continue;
				int el = eiUVWijk[lMat][1];
				int iuvw = eiUVWijk[lMat][2];
				int i = eiUVWijk[lMat][3];
				int j = eiUVWijk[lMat][4];
				int k = eiUVWijk[lMat][5];
				if (eh != el) {
					exit(EXIT_FAILURE);
				}
				std::complex<double> cSolStotal(0, 0), cSolA, cSolB, cSolS;


				//c_integrator.findC(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolA, size1);
				//std::cout << std::setprecision(15) << cSolA << std::endl;

				c_integrator.findC(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolA, size1, eval);


				d_integrator.findD(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolB, size1);
				//std::cout << std::setprecision(15) << cSolB << std::endl;
				tracker++;
				i0 = kMat - e->unknownsStart + 1;
				j0 = lMat - e->unknownsStart + 1;
				//perturbed k0
				//std::cout << "k0 factor: " << k0 << "cEps factor: " << cEps << std::endl;


				//////////////////////////////Frequency Perturbation///////////////////////////////////////////////
				//if (iHomCode == 1 && aCode == 1)
				//	e->cBr_HOPS(i0, j0) = 2.0*k0*cSolB*cEps;
				//else
				//	// If iHomCode != 1 || aCode != 1, then the cEps term is already included in the integration
				//	// of cSolB (as we have tensor/matrix material parameters) -- See below in the filling of cBrPAK
				//	e->cBr_HOPS(i0, j0) = 2.0 * k0 * cSolB;
				////--------------------------------------------------------------------------------------------------

				////////////////////////////Material Perturbation////////////////////////////////////////////////
				e->cBr_HOPS(i0, j0) = k0 * k0 * cSolB;
				////-------------------------------------------------------------------------------------------------

				// Due to symmetry, copy the value over
				e->cBr_HOPS(j0, i0) = e->cBr_HOPS(i0, j0);


				if (iHomCode == 1 && aCode == 1) {
					cSolA = cSolA / cMu;
					cSolB = cEps * cSolB;
				}
				if (std::signbit(double(vectorD[lMat])) == std::signbit(double(vectorD[kMat]))) {
					//if (vectorD[lMat]/vectorD[kMat] > 0){
					e->cArPAK(i0, j0) = cSolA;
					e->cBrPAK(i0, j0) = cSolB;
					//	e->cBr_HOPS(i0, j0) = cSolB;
					if (cSolStotal != std::complex<double>(0, 0)) {
						e->cSrPAK(i0, j0) = cSolStotal;
					}
				}
				else {
					e->cArPAK(i0, j0) = -cSolA;
					e->cBrPAK(i0, j0) = -cSolB;
					//e->cBr_HOPS(i0, j0) = -cSolB;
					if (cSolStotal != std::complex<double>(0, 0)) {
						e->cSrPAK(i0, j0) = -cSolStotal;
					}
				}
				if (i0 != j0) {
					if (std::signbit(double(vectorD[lMat])) == std::signbit(double(vectorD[kMat]))) { //for some reason, multiplying them and checking > 0 fails?????
					//	if (vectorD[lMat]/vectorD[kMat] > 0){
						e->cArPAK(j0, i0) = cSolA;
						e->cBrPAK(j0, i0) = cSolB;
						//e->cBr_HOPS(j0, i0) = cSolB;
						if (cSolStotal != std::complex<double>(0, 0)) {
							e->cSrPAK(j0, i0) = cSolStotal;
						}
					}
					else {
						e->cArPAK(j0, i0) = -cSolA;
						e->cBrPAK(j0, i0) = -cSolB;
						//	e->cBr_HOPS(j0, i0) = -cSolB;
						if (cSolStotal != std::complex<double>(0, 0)) {
							e->cSrPAK(j0, i0) = -cSolStotal;
						}
					}
				}
			}
			//i0_table[e->index - 1][iuvwh][ih][jh][kh] = kMat;

		}

		/*std::cout << "Element num: " << e->index << " took: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count()
			<< " milliseconds\n";*/
	}
	//output i0_table


	//check that cBrPak, cArPak are that same
	if (check_results) {
		std::string line;
		std::ifstream file;
		std::ifstream file2;
		file.open("../exampleFiles/debug_data/" + mesh_name + "/cArPAK_mat.txt");
		file2.open("../exampleFiles/debug_data/" + mesh_name + "/cBrPAK_mat.txt");
		bool messed_up = false;

		file.close(); file2.close();
		file.open("../exampleFiles/debug_data/" + mesh_name + "/cGr_mat.txt");
		for (int waveNo = 1; waveNo <= this->scatter1.numberOfWaves; ++waveNo) {
			for (int fcount = 1; fcount <= this->scatter1.nF; ++fcount) {
				for (int i0_it = 1; i0_it <= this->matDimCon; ++i0_it) {
					std::getline(file, line);
					std::vector<std::string> input_line = functions::split(line, ' ');
					if (abs(this->scatter1.cGr[fcount][i0_it].real() - std::stod(input_line[0])) > 1e-8) {
						std::cout << "Error for cGr at " << fcount << " and i0_it = " << i0_it << std::endl;
						messed_up = true;
						//std::cin.get();
					}
					if (abs(this->scatter1.cGr[fcount][i0_it].imag() - std::stod(input_line[1])) > 1e-8) {
						std::cout << "Error for cGr at " << fcount << " and i0_it = " << i0_it << std::endl;
						messed_up = true;
						//std::cin.get();
					}
				}
			}
			messed_up ? std::cout << "cGr is messed up!" << std::endl : std::cout << "cGr is okay!" << std::endl;
		}
		file.close();
	}
	//std::cout << "If no error was reported, it was successful!" << std::endl;
	//sparsepacking function
	//std::vector<int> iRow, jCol;
	std::vector<std::complex<double>> cArSparse, cSrSparse, cFEMrSparse;

	SPARSE_PACKING(iRow, jCol, nNz, 0, cArSparse, this->cBrSparse, cSrSparse);
	//now to do freq sweep
	double eps = 1e-7;
	//el_clear(this->elements); //clears elements from memory
	FrequencySweep fsweep;
	for (int fcount = 1; fcount <= this->scatter1.nF; ++fcount) {
		fsweep.getCMatrix(this->scatter1, cFEMrSparse, cArSparse, this->cBrSparse, nNz, fcount);
		sparseSolver solver;
		SpMat sparseMatrix(this->matDimCon, this->matDimCon);
		this->cAlpha.resize(this->matDimCon + 1);
		auto solver_start = std::chrono::high_resolution_clock::now();
		solver.doEverything(cFEMrSparse, iRow, jCol, nNz, this->matDimCon, this->scatter1.cGr[fcount], sparseMatrix, this->cAlpha);
		std::cout << "solving time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - solver_start).count() << std::endl;

		//compare with fortran here
		if (check_results) {
			bool wrong_calpha = false;
			std::ifstream file;
			std::string line;
			file.open("../exampleFiles/debug_data/" + mesh_name + "/cAlpha_mat.txt");
			for (int linecount = 1; linecount < this->cAlpha.size(); ++linecount) {
				getline(file, line);
				std::vector<std::string> input_line = functions::split(line, ' ');
				//first entry is a space
				if (abs(this->cAlpha[linecount].real() - std::stod(input_line[0])) > eps || abs(this->cAlpha[linecount].imag() - std::stod(input_line[1])) > eps) {
					//	std::cout << "value at: " << linecount << " in cAlpha (starting at 1) is wrong!" << std::endl;
					wrong_calpha = true;
				}
			}
			wrong_calpha ? std::cout << "Wrong cAlpha results! Does not match FORTRAN for frequency count: " << fcount << "!" << std::endl : std::cout << "cAlpha results for frequency count: " << fcount << " are same!" << std::endl;
			file.close();
		}
		std::cout << "Not writing to file, see line 2155 domain.cpp!" << std::endl;
		bool write_to_file = false;
		if (write_to_file) {
			std::ofstream file_out;
			std::cout << "Writing cAlpha to file..." << std::endl;
			if (!higher_order) {
				file_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha.txt");
				for (int i = 0; i < this->cAlpha.size(); ++i) {
					file_out << std::setprecision(15) << this->cAlpha[i].real() << " " << this->cAlpha[i].imag() << std::endl;
				}
			}
			else {
				file_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha_higher.txt");
				for (int i = 0; i < this->cAlpha.size(); ++i) {
					file_out << std::setprecision(15) << this->cAlpha[i].real() << " " << this->cAlpha[i].imag() << std::endl;
				}
			}
			file_out.close();
			std::cout << "Finished writing cAlpha to file." << std::endl;
			std::cout << "Writing cGr to file..." << std::endl;
			if (higher_order) {
				file_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
				for (int i = 0; i < this->scatter1.cGr[fcount].size(); ++i) {
					file_out << std::setprecision(15) << this->scatter1.cGr[fcount][i].real() << " " << this->scatter1.cGr[fcount][i].imag() << std::endl;
				}
				file_out.close();
				std::cout << "Finished writing cGr (higher order) to file." << std::endl;
			}
			else {
				file_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
				for (int i = 0; i < this->scatter1.cGr[fcount].size(); ++i) {
					file_out << std::setprecision(15) << this->scatter1.cGr[fcount][i].real() << " " << this->scatter1.cGr[fcount][i].imag() << std::endl;
				}
				file_out.close();
				std::cout << "Finished writing cGr to file." << std::endl;
			}
		}
		//std::endl;
	}
	//std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count()
	std::cout << "Finished testing for mesh: " << mesh_name << std::endl;
}
void Domain::_RHS_ONLY() {
	std::cout << "Entered matrix filling!" << std::endl;
	int dummy = 0;
	int nNz = 0;
	//bool isAdjoint = false;
	bool isAdjoint = this->sc.useAdjoint;
	//if (isAdjoint) mesh_name = mesh_name + "/adjoint";
	int tracker = 0;

	int pause = 0;
	//auto t1 = std::chrono::high_resolution_clock::now();
	//i0_table.resize(this->elements.size());

	//temp to test basiseval
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	//end temp to test basis eval

	//for (auto e = elements.begin(); e != elements.end(); ++e) {
#pragma omp parallel for num_threads(6)
	for (int e_index = 0; e_index < elements.size(); ++e_index) {
		auto e = &elements[e_index];
		//(i0_table[e->index - 1]) = std::vector<std::vector<std::vector<std::vector<int>>>>(4, std::vector<std::vector<std::vector<int>>>(e->expansion[0]+1, std::vector<std::vector<int>>(e->expansion[1]+1, std::vector<int>(e->expansion[2]+1))));
		//auto t1 = std::chrono::high_resolution_clock::now();
		Integral_c c_integrator(&(*e));
		Integral_d d_integrator(&(*e));
		Integral_g gInts(&(*e));

		e->cArPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		e->cBrPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);
		e->cSrPAK = matrix2d<std::complex<double>>(e->disconnected_dimension + 1, e->disconnected_dimension);

		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);


		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		e->uPowers = eval.get_u_samples(nglu, nu)[0];
		e->vPowers = eval.get_v_samples(nglv, nv)[0];
		e->wPowers = eval.get_w_samples(nglw, nw)[0];
		e->fuPowers = eval.get_u_samples(nglu, nu)[1];
		e->fvPowers = eval.get_v_samples(nglv, nv)[1];
		e->fwPowers = eval.get_w_samples(nglw, nw)[1];
		e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
		e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
		e->fpwPowers = eval.get_w_samples(nglw, nw)[2];
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);
		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);

		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}


			}
		}
		int iCon;
		//#pragma omp parallel for
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {

			iCon = abs(vectorD[kMat]);
			//if (iCon == 0) continue;
			int eh = eiUVWijk[kMat][1];
			int iuvwh = eiUVWijk[kMat][2];
			int ih = eiUVWijk[kMat][3];
			int jh = eiUVWijk[kMat][4];
			int kh = eiUVWijk[kMat][5];
			if (iuvwh == 1 && ih == 3 && jh == 0 && kh == 1) {
				int pause = 0;
			}
			//if (iuvwh == 1 && ih == 0 && jh == 4 && kh == 4) {
			//	int pause = 0;
			//}
			if (iCon == 0) continue;
			int i0, j0;
			if (isAdjoint) {
				for (int face = 1; face <= 6; ++face) {
					if (((this->elements[eh - 1]).abcList[face - 1] != 1) &&
						this->elements[eh - 1].PML_boundary_list[face - 1] != 1) {
						continue;
					}
					if (((this->elements[eh - 1]).PML_boundary_list[face - 1] == 1) && this->elements[eh - 1].materials.region == 2) {
						continue;
					}
					gInts.set_e(&(this->elements[eh - 1]));
					gInts.findGWave(&scatter1, vectorD, isAdjoint, iuvwh, ih, jh, kh, size1, kMat, iCon, face);
				}
			}
			else {
				int face = 1;
				if (e->materials.region == 1) {
					//should be for element eh
					gInts.set_e(&(this->elements[eh - 1]));
					gInts.findGWave(&scatter1, vectorD, isAdjoint, iuvwh, ih, jh, kh, size1,//. . .
						kMat, iCon, face);
				}
			}
		}

	}
}
bool Domain::check_q() {
	std::cout << "Checking Q (loading results from file)..." << std::endl;
	//bool check = false; 	
	std::ifstream file;
	std::string line;
	//load in the cAlphaForward
	//std::vector<std::complex<double>> cAlphaForward, cAlphaAdjoint, cGrForward, cGrAdjoint;
	std::vector<std::complex<double>> cGrForward, cGrAdjoint;

	file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha.txt");
	getline(file, line);
	std::vector<std::string> line_input;
	while (getline(file, line)) {
		line = line.substr(1, line.size() - 2);
		line_input = functions::split(line, ' ');
		cAlphaForward.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

	}
	file.close();
	file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha.txt");
	getline(file, line);
	while (getline(file, line)) {
		line = line.substr(1, line.size() - 2);
		line_input = functions::split(line, ' ');
		cAlphaAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

	}
	file.close();
	file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr.txt");
	getline(file, line);
	while (getline(file, line)) {
		line = line.substr(1, line.size() - 2);
		line_input = functions::split(line, ',');
		cGrForward.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	file.close();
	file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cGr.txt");
	getline(file, line);
	while (getline(file, line)) {
		line = line.substr(1, line.size() - 2);
		line_input = functions::split(line, ',');
		cGrAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
	}
	file.close();
	//do the summations for forward and adjoint
	std::complex<double> sumForward(0, 0), sumAdjoint(0, 0);
	for (int i = 0; i < cAlphaForward.size(); ++i) {
		sumForward += cAlphaForward[i] * std::conj(cGrAdjoint[i]);
	}
	for (int i = 0; i < cAlphaAdjoint.size(); ++i) {
		sumAdjoint += std::conj(cAlphaAdjoint[i]) * cGrForward[i];
	}
	//note: takes conjugate of the adjoint sum due to def of inner product
	/*if ((abs(sumAdjoint.real() - sumForward.real()) < 1e-5) && (abs(sumForward.imag() + sumAdjoint.imag()) < 1e-5)) {
		std::cout << "Qoi result: " << std::setprecision(15) << sumForward << std::endl;*/
	std::cout << "Qoi result: " << std::setprecision(15) << sumAdjoint << "and " << sumForward << std::endl;
	if (error) {


		//load qoi error results
		std::complex<double> total_qoi_error(0, 0);
		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt");
		while (getline(file, line)) {
			auto input_line = functions::split(line, ' ');
			sumAdjoint.real(sumAdjoint.real() + std::stod(input_line[0]));
			sumAdjoint.imag(sumAdjoint.imag() + std::stod(input_line[1]));
			total_qoi_error.real(total_qoi_error.real() + (std::stod(input_line[0])));
			total_qoi_error.imag(total_qoi_error.imag() + (std::stod(input_line[1])));
		}
		std::cout << "Qoi error total: " << std::setprecision(15) << total_qoi_error << std::endl;
		std::cout << "Qoi with Qoi Error added: " << std::setprecision(15) << sumAdjoint << std::endl;
		file.close();

		std::vector<std::complex<double>> cGrForward_higher, cGrAdjoint_higher, cAlphaForward_higher, cAlphaAdjoint_higher;

		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cAlpha_higher.txt");
		getline(file, line);
		std::vector<std::string> line_input;
		while (getline(file, line)) {
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ' ');
			cAlphaForward_higher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

		}
		file.close();
		file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cAlpha_higher.txt");
		getline(file, line);
		while (getline(file, line)) {
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ',');
			cAlphaAdjoint_higher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

		}
		file.close();
		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt");
		getline(file, line);
		while (getline(file, line)) {
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ',');
			cGrForward_higher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
		}
		file.close();
		file.open("../exampleFiles/debug_data/" + mesh_name + "/adjoint/results/cGr_higher_order.txt");
		getline(file, line);
		while (getline(file, line)) {
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ',');
			cGrAdjoint_higher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));
		}
		file.close();
		//do the summations for forward and adjoint
		std::complex<double> sumForward_higher(0, 0), sumAdjoint_higher(0, 0);
		for (int i = 0; i < cAlphaForward_higher.size(); ++i) {
			sumForward_higher += cAlphaForward_higher[i] * std::conj(cGrAdjoint_higher[i]);
		}
		for (int i = 0; i < cAlphaAdjoint_higher.size(); ++i) {
			sumAdjoint_higher += std::conj(cAlphaAdjoint_higher[i]) * cGrForward_higher[i];
		}
		//note: takes conjugate of the adjoint sum due to def of inner product
		//if ((abs(sumAdjoint_higher.real() - sumForward_higher.real()) < 1e-5) && (abs(sumForward_higher.imag() + sumAdjoint_higher.imag()) < 1e-5)) {
		std::cout << "Qoi result (higher): " << std::setprecision(15) << sumAdjoint_higher << std::endl;
		//}
		return true;
	}
	else return false;
}

void findPowersSens(const int& nu, const int& nv, const int& nw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
	matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
	matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
	//upowers goes from (0 to nglu+1, 0 to nu)
	//vpower (0 to nglv + 1, 0 to nv)
	//wPowers (0 to nglw + 1, 0 to nw)
	//fuPowers, fpuPowers same as uPowers
	//fvPowers, fpvPowers same as vPowers
	double x;
	for (int i = 0; i <= nglu + 1; ++i) {
		x = xglu[i];
		uPowers(i, 0) = 1.0;
		uPowers(i, 1) = x;
		fuPowers(i, 0) = 1.0 - x;
		fuPowers(i, 1) = 1.0 + x;
		fpuPowers(i, 0) = -1.0;
		fpuPowers(i, 1) = 1.0;
		for (int j = 2; j <= nu; ++j) {
			uPowers(i, j) = uPowers(i, j - 1) * x;
			if (j % 2 == 0) {
				fuPowers(i, j) = uPowers(i, j) - 1.0;
				fpuPowers(i, j) = j * uPowers(i, j - 1);
			}
			else {
				fuPowers(i, j) = uPowers(i, j) - x;
				fpuPowers(i, j) = j * uPowers(i, j - 1) - 1.0;
			}
		}
	}
	for (int i = 0; i <= nglv + 1; ++i) {
		x = xglv[i];
		vPowers(i, 0) = 1.0;
		vPowers(i, 1) = x;
		fvPowers(i, 0) = 1.0 - x;
		fvPowers(i, 1) = 1.0 + x;
		fpvPowers(i, 0) = -1.0;
		fpvPowers(i, 1) = 1.0;
		for (int j = 2; j <= nv; ++j) {
			vPowers(i, j) = vPowers(i, j - 1) * x;
			if (j % 2 == 0) {
				fvPowers(i, j) = vPowers(i, j) - 1.0;
				fpvPowers(i, j) = j * vPowers(i, j - 1);
			}
			else {
				fvPowers(i, j) = vPowers(i, j) - x;
				fpvPowers(i, j) = j * vPowers(i, j - 1) - 1.0;
			}
		}
	}
	for (int i = 0; i <= nglw + 1; ++i) {
		x = xglw[i];
		wPowers(i, 0) = 1.0;
		wPowers(i, 1) = x;
		fwPowers(i, 0) = 1.0 - x;
		fwPowers(i, 1) = 1.0 + x;
		fpwPowers(i, 0) = -1.0;
		fpwPowers(i, 1) = 1.0;
		for (int j = 2; j <= nw; ++j) {
			wPowers(i, j) = wPowers(i, j - 1) * x;
			if (j % 2 == 0) {
				fwPowers(i, j) = wPowers(i, j) - 1.0;
				fpwPowers(i, j) = j * wPowers(i, j - 1);
			}
			else {
				fwPowers(i, j) = wPowers(i, j) - x;
				fpwPowers(i, j) = j * wPowers(i, j - 1) - 1.0;
			}
		}
	}

}

void Domain::spherical_sensitivity(const int& myElementCount, std::vector<Element>& elements, const int& matDimDis,
	const std::vector<int>& vectorD, const std::vector<std::vector<int>>& eiUVWijk, const int& matDimCon, const std::vector<std::complex<double>>& cAlphaForward,
	std::vector<std::complex<double>>& cAlphaAdjoint, std::complex<double>& cSensitivity, const int& numberOfFrequencies,
	std::vector<double>& k0List, const int& indicatorBasisType, Scatter* scatter1, const std::vector<Facet>& facets) {
	//note: cAlphaForward == forward solve results (potentially read in from file)
	//note: cAlphaAdjoint == adjoint solve results read in from file (potentially just computed) 
	std::complex<double> cj(0, 1);
	cSensitivity = 0;
	int kuvw = 0;

	int size1 = 1;
	int iHomCode = 1;
	int aCode = 1;
	int integralType = 1;
	//currently only works with one single wave
	int waveNumber = 1;
	int freqNo = 1;
	double val1 = -scatter1->K0[freqNo];
	auto cEteta = scatter1->waveTheta[waveNumber];
	auto cEfi = scatter1->wavePhi[waveNumber];
	auto iFiAr = scatter1->iPhiAr[waveNumber];
	auto iThetaAr = scatter1->iThetaAr[waveNumber];
	auto nAr = scatter1->nAr[waveNumber];
	std::vector<double> iR(4);
	std::vector<std::complex<double>> asec2(4), eVector(4), hVector(4), asec3(4);
	asec2[1] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][1] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][1]);
	asec2[2] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][2] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][2]);
	asec2[3] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][3] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][3]);
	eVector[1] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][1] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][1];
	eVector[2] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][2] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][2];
	eVector[3] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][3] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][3];
	products::cross(scatter1->nAr[waveNumber], eVector, hVector);
	asec3[1] = -cj * hVector[1];
	asec3[2] = -cj * hVector[2];
	asec3[3] = -cj * hVector[3];
	iR[1] = -scatter1->nAr[waveNumber][1];
	iR[2] = -scatter1->nAr[waveNumber][2];
	iR[3] = -scatter1->nAr[waveNumber][3];
	std::complex<double> cEpsilon;
	std::complex<double> cMu;
	std::complex<double> cSolB, cSolGTotal;
	for (auto e = elements.begin(); e != elements.end(); ++e) {
		//while (e->index != 433) { ++e;}
		if (e->materials.hcode != 1 || e->materials.icode != 1) continue; //not sure about this line
		//if (aCode == 0) {
		//	kuvw = e->materials.KuvwA;
		//	//size1 = 9 - 3 * e->materials.sym;
		//}
		//else {
		kuvw = e->materials.Kuvw;
		//size1 = 1;
	//}
		Integral_g gInts(&(*e));

		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		//get the faces of the element that are connected to another element

		for (auto face = e->facet_indices.begin(); face != e->facet_indices.end(); ++face) {
			//std::cout << "Elem: " << e->index << " Face: " << *face << std::endl;
			int connectedTo; //the index of the other element
			int connectedToType;
			Facet this_facet = facets[*face - 1];
			int iFace;

			/*this_facet.element_indices.first == e->index ? connectedTo = this_facet.element_indices.second
				: connectedTo = this_facet.element_indices.first;*/
			if (this_facet.element_indices.first == e->index) {
				connectedTo = this_facet.element_indices.second;
				iFace = this_facet.local_facet_indices.first;
			}
			else {
				connectedTo = this_facet.element_indices.first;
				iFace = this_facet.local_facet_indices.second;
			}

			//if not connected directly to another elements
			if (connectedTo != -1) { //might need to -1, instead of 0
				connectedToType = elements[connectedTo - 1].materials.region;
				if (connectedToType == 0 && e->materials.region == 1) { //might need to change 0 and 1 to be different
					continue;
				}
			}
			else {
				continue;
			}
			int nglu = e->quadrature[0];
			int nglv = e->quadrature[1];
			int nglw = e->quadrature[2];
			if (iFace == 1 || iFace == 2) nglu = 1;
			else if (iFace == 3 || iFace == 4) nglv = 1;
			else if (iFace == 5 || iFace == 6) nglw = 1;
			//call new G-L integration params
			std::vector<double> xglu, xglv, xglw, wglu, wglv, wglw;
			functions::gaussk(nglu, xglu, e->wglu);
			functions::gaussk(nglv, xglv, e->wglv);
			functions::gaussk(nglw, xglw, e->wglw);
			e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
			e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
			e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
			e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
			e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
			e->uPowers = matrix2d<double>(nglu + 1, nu);
			e->vPowers = matrix2d<double>(nglv + 1, nv);
			e->wPowers = matrix2d<double>(nglw + 1, nw);
			e->fuPowers = matrix2d<double>(nglu + 1, nu);
			e->fvPowers = matrix2d<double>(nglv + 1, nv);
			e->fwPowers = matrix2d<double>(nglw + 1, nw);
			e->fpuPowers = matrix2d<double>(nglu + 1, nu);
			e->fpvPowers = matrix2d<double>(nglv + 1, nv);
			e->fpwPowers = matrix2d<double>(nglw + 1, nw);
			e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
			e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
			e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);
			e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
			e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
			e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
			//modify the xglu,v,w vectors
			switch (iFace) {
			case 1:
				xglu[1] = -1;
				break;
			case 2:
				xglu[1] = 1;
				break;
			case 3:
				xglv[1] = -1;
				break;
			case 4:
				xglv[1] = 1;
				break;
			case 5:
				xglw[1] = -1;
				break;
			case 6:
				xglw[1] = 1;
				break;
			}
			unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
				e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix,
				e->awMatrix, e->jacobian, e->nRs, integralType);
			switch (indicatorBasisType) { //current set to only do kolundzija
			case 1:
				findPowersSens(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers,
					e->vPowers, e->wPowers, e->fuPowers, e->fvPowers,
					e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);
				break;
			case 2:
				//do legendre powers
				break;
			}
			if (iHomCode == 1 && aCode == 1) {
				cMu = std::complex<double>(1.0, 0.0);
				if (e->materials.hcode == 1 && e->materials.icode == 1) {
					cEpsilon = e->materials.epsr_list[0][0];
				}
				else {
					cEpsilon = 1;
					e->materials.hcode = 1;
					e->materials.icode = 1;
				}
				//e->MuRelInt(1, 1, 1, 1) = std::complex<double>(1.0, 0.0);
				e->MuRelInt(1, 1, 1, 1) = 1;
				e->EpsRelInt(1, 1, 1, 1) = cEpsilon;
			}
			for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
				int iCon = abs(vectorD[kMat]);
				if (iCon == 0) continue;
				int eh = eiUVWijk[kMat][1];
				int iuvwh = eiUVWijk[kMat][2];
				int ih = eiUVWijk[kMat][3];
				int jh = eiUVWijk[kMat][4];
				int kh = eiUVWijk[kMat][5];
				gInts.set_e(&elements[eh - 1]);
				cSolGTotal = 0;
				//gInts.findGInc_EpsOnly(scatter1, vectorD, iuvwh, ih, jh, kh, size1, k0List[freqNo], freqNo, cSolGTotal, asec2, hVector, asec3, val1);
				gInts.findGInc_EpsOnly(scatter1, vectorD, iuvwh, ih, jh, kh, size1, waveNumber, freqNo, cSolGTotal, asec2, hVector, asec3, val1);
				cSensitivity += cSolGTotal * std::conj(cAlphaAdjoint[iCon - 1]);
				for (int lMat = kMat; lMat <= e->unknownsEnd; ++lMat) {
					int jCon = abs(vectorD[lMat]);
					if (jCon == 0) {
						continue;
					}
					int el = eiUVWijk[lMat][1];
					int iuvw = eiUVWijk[lMat][2];
					int i = eiUVWijk[lMat][3];
					int j = eiUVWijk[lMat][4];
					int k = eiUVWijk[lMat][5];
					cSolB = 0;
					// ------>fix			//findD
					Integral_d d_integrator(&elements[el - 1]);
					d_integrator.findD(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolB, size1);
					cSensitivity += cSolB * cAlphaForward[jCon - 1] * std::conj(cAlphaAdjoint[iCon - 1]) * (cEpsilon - std::complex<double>(1.0, 0.0));

				}
			}

		}
	}
	std::cout << std::endl << "Sensitivity: " << std::setprecision(15) << cSensitivity << std::endl;
}



//fix find powers - done.
void Domain::element_error_f(int u_order, int v_order, int w_order, Domain& dom2) {
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	std::ifstream file;
	std::string line;
	std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
	read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
	double k0 = this->scatter1.K0[1];
	std::complex<double> g_sum, e_sum;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//if (e->index != 82) continue;
		///////////////////////////////////////////
#pragma region r1

		std::cout << "For element " << e->index << std::endl;
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		bool isAdjoint = false;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);


		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1; //kolund basis

		e->uPowers = matrix2d<double>(nglu + 1, nu);
		e->vPowers = matrix2d<double>(nglv + 1, nv);
		e->wPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowers = matrix2d<double>(nglu + 1, nu);
		e->fvPowers = matrix2d<double>(nglv + 1, nv);
		e->fwPowers = matrix2d<double>(nglw + 1, nw);
		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		//if indicator = 1, then kolun, if indicator = 2, then legendre
		if (indicatorBasisType == 1) {
			e->uPowers = eval.get_u_samples(nglu, nu)[0];
			e->vPowers = eval.get_v_samples(nglv, nv)[0];
			e->wPowers = eval.get_w_samples(nglw, nw)[0];
			e->fuPowers = eval.get_u_samples(nglu, nu)[1];
			e->fvPowers = eval.get_v_samples(nglv, nv)[1];
			e->fwPowers = eval.get_w_samples(nglw, nw)[1];
			e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
			e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
			e->fpwPowers = eval.get_w_samples(nglw, nw)[2];
		}
		//kolund basis functions
		//findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);


		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}


			}
		}
#pragma endregion ini_stuff 
		//////////////////////////////////////////	
		Integral_g IntG(&(*e));
		Integral_c Intc(&(*e));
		Integral_d Intd(&(*e));
		//iterate through higher order unknowns
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			int eh = this->eiUVWijk[kMat][1];
			int iuvwh = this->eiUVWijk[kMat][2];
			int ih = this->eiUVWijk[kMat][3];
			int jh = this->eiUVWijk[kMat][4];
			int kh = this->eiUVWijk[kMat][5];
			int D1 = abs(this->vectorD[kMat]);
			//if non-pml: compute G integral part
			if (qoi_error::is_higher(u_order, v_order, w_order, iuvwh, ih, jh, kh)) {
				if (e->materials.region == 1) {
					int face = 1;
					std::complex<double> cGr_prev = this->scatter1.cGr[1][D1];
					IntG.findGWave(&scatter1, vectorD, false, iuvwh, ih, jh, kh, size1, kMat, D1, face);
					e->element_error += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					g_sum += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				}
				//iterate through lower order unknowns
				for (int lMat = dom2.elements[e->index - 1].unknownsStart; lMat <= dom2.elements[e->index - 1].unknownsEnd; ++lMat) {
					int el = dom2.eiUVWijk[lMat][1];
					int iuvwl = dom2.eiUVWijk[lMat][2];
					int il = dom2.eiUVWijk[lMat][3];
					int jl = dom2.eiUVWijk[lMat][4];
					int kl = dom2.eiUVWijk[lMat][5];
					int D2 = abs(dom2.vectorD[lMat]);

					std::complex<double> cSolA, cSolB;
					//Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1);
					//Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1, eval);
					Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					if (std::signbit(double(dom2.vectorD[lMat])) != std::signbit(double(vectorD[kMat]))) {
						cSolA = -cSolA;
						cSolB = -cSolB;
					}
					if (iHomCode == 1 && aCode == 1) {
						cSolA = cSolA / cMu;
						cSolB = cEps * cSolB;
					}
					e->element_error -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e->element_error += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					e_sum -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e_sum += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

				}
			}
		}
	}
	std::cout << "g sum: " << g_sum << std::endl;
	std::complex<double> sumAdjoint;
	for (int i = 0; i < cAlphaAdj.size(); ++i) {
		sumAdjoint += std::conj(cAlphaAdj[i]) * this->scatter1.cGr[1][i];
	}
	std::cout << "e sum: " << e_sum << std::endl;

	std::cout << "total error: " << e_sum + g_sum << std::endl;
	//print out qoi errors per element
	std::cout << "Saving QoI Error per element..." << std::endl;
	std::ofstream file_qoi_out;
	file_qoi_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt");
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//one element per line: real imag
		file_qoi_out << e->element_error.real() << " " << e->element_error.imag() << std::endl;
	}
	std::cout << "Done saving QoI error." << std::endl;
}
//fix find powers in element error
void Domain::element_error_basis(int u_order, int v_order, int w_order, Domain& dom2) {
	std::ifstream file;
	std::string line;
	std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
	read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
	double k0 = this->scatter1.K0[1];
	std::complex<double> g_sum, e_sum;
	//std::vector<std::unordered_map<std::string, std::complex<double>>> basis_maps(this->elements.size());
	std::vector<std::complex<double>> basis_error_coeffs(this->matDimDis + 1);
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//if (e->index != 648) continue;
		///////////////////////////////////////////
#pragma region r1

		std::cout << "For element " << e->index << std::endl;
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		bool isAdjoint = false;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);


		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1; //kolund basis

		e->uPowers = matrix2d<double>(nglu + 1, nu);
		e->vPowers = matrix2d<double>(nglv + 1, nv);
		e->wPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowers = matrix2d<double>(nglu + 1, nu);
		e->fvPowers = matrix2d<double>(nglv + 1, nv);
		e->fwPowers = matrix2d<double>(nglw + 1, nw);
		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		//if indicator = 1, then kolun, if indicator = 2, then legendre
		//if (indicatorBasisType == 1)
			//kolund basis functions
		//	findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);
//
//
//
		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}


			}
		}
#pragma endregion ini_stuff 
		//////////////////////////////////////////	
		Integral_g IntG(&(*e));
		Integral_c Intc(&(*e));
		Integral_d Intd(&(*e));
		//iterate through higher order unknowns
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			int eh = this->eiUVWijk[kMat][1];
			int iuvwh = this->eiUVWijk[kMat][2];
			int ih = this->eiUVWijk[kMat][3];
			int jh = this->eiUVWijk[kMat][4];
			int kh = this->eiUVWijk[kMat][5];
			int D1 = abs(this->vectorD[kMat]);
			//if non-pml: compute G integral part
			if (qoi_error::is_higher(u_order, v_order, w_order, iuvwh, ih, jh, kh)) {
				////convert ih, jh, kh to hash
				//std::string hash;
				//std::ostringstream convert;
				//convert << ih << "," << jh << "," << kh;
				//hash = convert.str();
				//auto map_it = basis_maps[e->index - 1].find(hash);
				//if (map_it == basis_maps[e->index - 1].end()) basis_maps[e->index - 1].insert({ hash, std::complex<double>(0, 0) });
				if (e->materials.region == 1) {
					int face = 1;
					std::complex<double> cGr_prev = this->scatter1.cGr[1][D1];
					IntG.findGWave(&scatter1, vectorD, false, iuvwh, ih, jh, kh, size1, kMat, D1, face);
					e->element_error += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					g_sum += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					basis_error_coeffs[kMat] += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					//basis_maps[e->index-1][hash] += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				}
				//iterate through lower order unknowns
				for (int lMat = dom2.elements[e->index - 1].unknownsStart; lMat <= dom2.elements[e->index - 1].unknownsEnd; ++lMat) {
					int el = dom2.eiUVWijk[lMat][1];
					int iuvwl = dom2.eiUVWijk[lMat][2];
					int il = dom2.eiUVWijk[lMat][3];
					int jl = dom2.eiUVWijk[lMat][4];
					int kl = dom2.eiUVWijk[lMat][5];
					int D2 = abs(dom2.vectorD[lMat]);

					std::complex<double> cSolA, cSolB;
					/*Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1);
					Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);*/
					if (std::signbit(double(dom2.vectorD[lMat])) != std::signbit(double(vectorD[kMat]))) {
						cSolA = -cSolA;
						cSolB = -cSolB;
					}
					if (iHomCode == 1 && aCode == 1) {
						cSolA = cSolA / cMu;
						cSolB = cEps * cSolB;
					}
					e->element_error -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e->element_error += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					e_sum -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e_sum += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					basis_error_coeffs[kMat] -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					basis_error_coeffs[kMat] += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					//basis_maps[e->index - 1][hash] -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					//basis_maps[e->index - 1][hash] += k0 * k0* cSolB*cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);



				}
			}
		}
	}
	/*std::cout << "g sum: " << g_sum << std::endl;
	std::complex<double> sumAdjoint;
	for (int i = 0; i < cAlphaAdj.size(); ++i) {
		sumAdjoint += std::conj(cAlphaAdj[i]) * this->scatter1.cGr[1][i];
	}
	std::cout << "e sum: " << e_sum << std::endl;
	std::cout << "total error: " << e_sum + g_sum << std::endl;*/
	//print out qoi errors per element
	std::cout << "Saving QoI Error per basis..." << std::endl;
	std::ofstream file_qoi_out;
	file_qoi_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt");
	for (auto coef = basis_error_coeffs.begin(); coef != basis_error_coeffs.end(); ++coef) {
		//one element per line: real imag
		file_qoi_out << (*coef).real() << " " << (*coef).imag() << std::endl;
	}
	std::cout << "Done saving QoI error per basis." << std::endl;
}


void Domain::_HOPS_EPS(int u_order, int v_order, int w_order, Domain& dom2) {
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	std::ifstream file;
	std::string line;
	std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
	read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
	double k0 = this->scatter1.K0[1];
	std::complex<double> g_sum, e_sum;
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//if (e->index != 82) continue;
		///////////////////////////////////////////

		std::cout << "For element " << e->index << std::endl;
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		bool isAdjoint = false;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);


		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1; //kolund basis

		e->uPowers = matrix2d<double>(nglu + 1, nu);
		e->vPowers = matrix2d<double>(nglv + 1, nv);
		e->wPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowers = matrix2d<double>(nglu + 1, nu);
		e->fvPowers = matrix2d<double>(nglv + 1, nv);
		e->fwPowers = matrix2d<double>(nglw + 1, nw);
		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		//if indicator = 1, then kolun, if indicator = 2, then legendre
		if (indicatorBasisType == 1) {
			e->uPowers = eval.get_u_samples(nglu, nu)[0];
			e->vPowers = eval.get_v_samples(nglv, nv)[0];
			e->wPowers = eval.get_w_samples(nglw, nw)[0];
			e->fuPowers = eval.get_u_samples(nglu, nu)[1];
			e->fvPowers = eval.get_v_samples(nglv, nv)[1];
			e->fwPowers = eval.get_w_samples(nglw, nw)[1];
			e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
			e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
			e->fpwPowers = eval.get_w_samples(nglw, nw)[2];
		}
		//kolund basis functions
		//findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);


		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}


			}
		}
		//////////////////////////////////////////	
		Integral_g IntG(&(*e));
		Integral_c Intc(&(*e));
		Integral_d Intd(&(*e));
		//iterate through higher order unknowns
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			if (e->materials.region == 1) {
				int eh = this->eiUVWijk[kMat][1];
				int iuvwh = this->eiUVWijk[kMat][2];
				int ih = this->eiUVWijk[kMat][3];
				int jh = this->eiUVWijk[kMat][4];
				int kh = this->eiUVWijk[kMat][5];
				int D1 = abs(this->vectorD[kMat]);
				//if non-pml: compute G integral part
				//if (qoi_error::is_higher(u_order, v_order, w_order, iuvwh, ih, jh, kh)) {
					//if (e->materials.region == 1) {
				int face = 1;
				std::complex<double> cGr_prev = this->scatter1.cGr[1][D1];
				IntG.findGWave(&scatter1, vectorD, false, iuvwh, ih, jh, kh, size1, kMat, D1, face);
				e->element_error += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				g_sum += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				//	}
					//iterate through lower order unknowns
				for (int lMat = dom2.elements[e->index - 1].unknownsStart; lMat <= dom2.elements[e->index - 1].unknownsEnd; ++lMat) {
					int el = dom2.eiUVWijk[lMat][1];
					int iuvwl = dom2.eiUVWijk[lMat][2];
					int il = dom2.eiUVWijk[lMat][3];
					int jl = dom2.eiUVWijk[lMat][4];
					int kl = dom2.eiUVWijk[lMat][5];
					int D2 = abs(dom2.vectorD[lMat]);

					std::complex<double> cSolA, cSolB;
					//Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1);
					//Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1, eval);
					Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					if (std::signbit(double(dom2.vectorD[lMat])) != std::signbit(double(vectorD[kMat]))) {
						cSolA = -cSolA;
						cSolB = -cSolB;
					}
					/*if (iHomCode == 1 && aCode == 1) {
						cSolA = cSolA / cMu;
						cSolB = cEps * cSolB;
					}*/
					e->element_error -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e->element_error += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					e_sum -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e_sum += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

				}
				//}
			}
		}
	}
	std::cout << "g sum: " << g_sum << std::endl;
	std::complex<double> sumAdjoint;
	for (int i = 0; i < cAlphaAdj.size(); ++i) {
		sumAdjoint += std::conj(cAlphaAdj[i]) * this->scatter1.cGr[1][i];
	}
	std::cout << "e sum: " << e_sum << std::endl;

	std::cout << "total error: " << e_sum + g_sum << std::endl;
	//print out qoi errors per element
	std::cout << "Saving QoI Error per element..." << std::endl;
	std::ofstream file_qoi_out;
	file_qoi_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt");
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//one element per line: real imag
		file_qoi_out << e->element_error.real() << " " << e->element_error.imag() << std::endl;
	}
	std::cout << "Done saving QoI error." << std::endl;
}



void Domain::element_error_improved(Domain& dom2, std::vector<std::complex<double>>& cAlphaFor, std::vector<std::complex<double>>& cAlphaAdj) {
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	std::ifstream file;
	std::string line;
	//std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
	//read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
	double k0 = this->scatter1.K0[1];
	std::complex<double> g_sum, e_sum;
#pragma omp parallel for num_threads(8)
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		e->cGr_eps_el = std::vector<std::complex<double>>(this->scatter1.cGr[1].size(), 0.0);
		int u_order = dom2.elements[e->index - 1].expansion[0];
		int v_order = dom2.elements[e->index - 1].expansion[1];
		int w_order = dom2.elements[e->index - 1].expansion[2];
		//std::cout << "For element " << e->index << std::endl;
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		bool isAdjoint = false;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1;
		e->uPowers = matrix2d<double>(nglu + 1, nu);
		e->vPowers = matrix2d<double>(nglv + 1, nv);
		e->wPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowers = matrix2d<double>(nglu + 1, nu);
		e->fvPowers = matrix2d<double>(nglv + 1, nv);
		e->fwPowers = matrix2d<double>(nglw + 1, nw);
		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		if (indicatorBasisType == 1) {
			e->uPowers = eval.get_u_samples(nglu, nu)[0];
			e->vPowers = eval.get_v_samples(nglv, nv)[0];
			e->wPowers = eval.get_w_samples(nglw, nw)[0];
			e->fuPowers = eval.get_u_samples(nglu, nu)[1];
			e->fvPowers = eval.get_v_samples(nglv, nv)[1];
			e->fwPowers = eval.get_w_samples(nglw, nw)[1];
			e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
			e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
			e->fpwPowers = eval.get_w_samples(nglw, nw)[2];
		}
		//kolund basis functions
		//findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);


		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}
			}
		}
		//////////////////////////////////////////	
		Integral_g IntG(&(*e));
		Integral_c Intc(&(*e));
		Integral_d Intd(&(*e));
		//iterate through higher order unknowns
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			int eh = this->eiUVWijk[kMat][1];
			int iuvwh = this->eiUVWijk[kMat][2];
			int ih = this->eiUVWijk[kMat][3];
			int jh = this->eiUVWijk[kMat][4];
			int kh = this->eiUVWijk[kMat][5];
			int D1 = abs(this->vectorD[kMat]);
			//if non-pml: compute G integral part
			if (qoi_error::is_higher(u_order, v_order, w_order, iuvwh, ih, jh, kh)) {
				if (e->materials.region == 1) {
					int face = 1;
					std::complex<double> cGr_prev = this->scatter1.cGr[1][D1];
					IntG.findGWave(&scatter1, vectorD, false, iuvwh, ih, jh, kh, size1, kMat, D1, face);
					e->element_error += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					//g_sum += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				}
				//iterate through lower order unknowns
				for (int lMat = dom2.elements[e->index - 1].unknownsStart; lMat <= dom2.elements[e->index - 1].unknownsEnd; ++lMat) {
					int el = dom2.eiUVWijk[lMat][1];
					int iuvwl = dom2.eiUVWijk[lMat][2];
					int il = dom2.eiUVWijk[lMat][3];
					int jl = dom2.eiUVWijk[lMat][4];
					int kl = dom2.eiUVWijk[lMat][5];
					int D2 = abs(dom2.vectorD[lMat]);

					std::complex<double> cSolA, cSolB;
					Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1, eval);
					Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					if (std::signbit(double(dom2.vectorD[lMat])) != std::signbit(double(vectorD[kMat]))) {
						cSolA = -cSolA;
						cSolB = -cSolB;
					}
					if (iHomCode == 1 && aCode == 1) {
						cSolA = cSolA / cMu;
						cSolB = cEps * cSolB;
					}
					e->element_error -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					e->element_error += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					//e_sum -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					//e_sum += k0 * k0* cSolB*cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

				}
			}
		}
	}
	//std::cout << "g sum: " << g_sum << std::endl;
	std::complex<double> sumAdjoint;
	for (int i = 0; i < cAlphaAdj.size(); ++i) {
		sumAdjoint += std::conj(cAlphaAdj[i]) * this->scatter1.cGr[1][i];
	}
	//std::cout << "e sum: " << e_sum << std::endl;

	//std::cout << "total error: " << e_sum + g_sum << std::endl;
	//print out qoi errors per element
	std::cout << "Saving QoI Error per element..." << std::endl;
	std::ofstream file_qoi_out;
	file_qoi_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_error.txt");
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		//one element per line: real imag
		file_qoi_out << e->element_error.real() << " " << e->element_error.imag() << std::endl;
	}
	std::cout << "Done saving QoI error." << std::endl;
}

void Domain::basis_error_improved(Domain& dom2, std::vector<std::complex<double>>& cAlphaFor, std::vector<std::complex<double>>& cAlphaAdj) {
	BasisEval eval;
	if (eval.basisType == 1) { //legendre basis
		eval.setup_legendre();
	}
	std::ifstream file;
	std::string line;
	//std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
	//read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
	double k0 = this->scatter1.K0[1];
	std::complex<double> g_sum, e_sum;
	std::vector<std::complex<double>> basis_error_coeffs(this->matDimDis + 1);
#pragma omp parallel for num_threads(8)
	for (auto e = this->elements.begin(); e != this->elements.end(); ++e) {
		e->cGr_eps_el = std::vector<std::complex<double>>(this->scatter1.cGr[1].size(), 0.0);
		int u_order = dom2.elements[e->index - 1].expansion[0];
		int v_order = dom2.elements[e->index - 1].expansion[1];
		int w_order = dom2.elements[e->index - 1].expansion[2];
		//std::cout << "For element " << e->index << std::endl;
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];
		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		int iHomCode = e->materials.hcode;
		int aCode = e->materials.icode;
		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		std::complex<double> cEps;
		std::complex<double> cMu;
		bool isAdjoint = false;
		int kuvw;
		int size1;
		if (aCode == 0) {
			kuvw = e->materials.KuvwA;
			size1 = 9 - 3 * e->materials.sym;
		}
		else {
			kuvw = e->materials.Kuvw;
			size1 = 1;
		}
		std::vector<double> xglu, xglv, xglw; //coords

		functions::gaussk(nglu, xglu, e->wglu);
		functions::gaussk(nglv, xglv, e->wglv);
		functions::gaussk(nglw, xglw, e->wglw);
		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);
		int indicatorBasisType = 1;
		e->uPowers = matrix2d<double>(nglu + 1, nu);
		e->vPowers = matrix2d<double>(nglv + 1, nv);
		e->wPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowers = matrix2d<double>(nglu + 1, nu);
		e->fvPowers = matrix2d<double>(nglv + 1, nv);
		e->fwPowers = matrix2d<double>(nglw + 1, nw);
		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);

		if (indicatorBasisType == 1) {
			e->uPowers = eval.get_u_samples(nglu, nu)[0];
			e->vPowers = eval.get_v_samples(nglv, nv)[0];
			e->wPowers = eval.get_w_samples(nglw, nw)[0];
			e->fuPowers = eval.get_u_samples(nglu, nu)[1];
			e->fvPowers = eval.get_v_samples(nglv, nv)[1];
			e->fwPowers = eval.get_w_samples(nglw, nw)[1];
			e->fpuPowers = eval.get_u_samples(nglu, nu)[2];
			e->fpvPowers = eval.get_v_samples(nglv, nv)[2];
			e->fpwPowers = eval.get_w_samples(nglw, nw)[2];
		}
		//kolund basis functions
		//findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);


		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);


		if ((iHomCode == 1) && (aCode == 1)) {
			cEps = e->materials.epsr_list[0][0]; //homogen
			cMu = e->materials.mur_list[0][0];//homogen
			e->MuRelInt(1, 1, 1, 1) = cMu;
			e->EpsRelInt(1, 1, 1, 1) = cEps;
			if (isAdjoint == true) {

				cMu.imag(-cMu.imag());
				cEps.imag(-cEps.imag());
				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
			}
		}
		else {
			if (kuvw != -1) {
				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);

				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);

				if (aCode == 0) {
					if (iHomCode == 0) {
						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
						}
					}
					else {
						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
						if (isAdjoint) {
							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
						}
					}
				}
				else {
					iHomCode = 1; //homogeneous, anisotropic
					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
					if (isAdjoint) {
						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
					}
				}
			}
		}
		//////////////////////////////////////////	
		Integral_g IntG(&(*e));
		Integral_c Intc(&(*e));
		Integral_d Intd(&(*e));
		//iterate through higher order unknowns
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			int eh = this->eiUVWijk[kMat][1];
			int iuvwh = this->eiUVWijk[kMat][2];
			int ih = this->eiUVWijk[kMat][3];
			int jh = this->eiUVWijk[kMat][4];
			int kh = this->eiUVWijk[kMat][5];
			int D1 = abs(this->vectorD[kMat]);
			//if non-pml: compute G integral part
			if (qoi_error::is_higher(u_order, v_order, w_order, iuvwh, ih, jh, kh)) {
				if (e->materials.region == 1) {
					int face = 1;
					std::complex<double> cGr_prev = this->scatter1.cGr[1][D1];
					IntG.findGWave(&scatter1, vectorD, false, iuvwh, ih, jh, kh, size1, kMat, D1, face);
					basis_error_coeffs[kMat] += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					//e->element_error += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
					//g_sum += (this->scatter1.cGr[1][D1] - cGr_prev) * std::conj(cAlphaAdj[D1]);
				}
				//iterate through lower order unknowns
				for (int lMat = dom2.elements[e->index - 1].unknownsStart; lMat <= dom2.elements[e->index - 1].unknownsEnd; ++lMat) {
					int el = dom2.eiUVWijk[lMat][1];
					int iuvwl = dom2.eiUVWijk[lMat][2];
					int il = dom2.eiUVWijk[lMat][3];
					int jl = dom2.eiUVWijk[lMat][4];
					int kl = dom2.eiUVWijk[lMat][5];
					int D2 = abs(dom2.vectorD[lMat]);

					std::complex<double> cSolA, cSolB;
					Intc.findC(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolA, size1, eval);
					Intd.findD(iuvwh, iuvwl, ih, jh, kh, il, jl, kl, cSolB, size1);
					if (std::signbit(double(dom2.vectorD[lMat])) != std::signbit(double(vectorD[kMat]))) {
						cSolA = -cSolA;
						cSolB = -cSolB;
					}
					if (iHomCode == 1 && aCode == 1) {
						cSolA = cSolA / cMu;
						cSolB = cEps * cSolB;
					}
					/*	e->element_error -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
						e->element_error += k0 * k0* cSolB*cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);*/

						//e_sum -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
						//e_sum += k0 * k0* cSolB*cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);

					basis_error_coeffs[kMat] -= cSolA * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
					basis_error_coeffs[kMat] += k0 * k0 * cSolB * cAlphaFor[D2] * std::conj(cAlphaAdj[D1]);
				}
			}
		}
	}
	std::cout << "Saving QoI Error per basis..." << std::endl;
	std::ofstream file_qoi_out;
	file_qoi_out.open("../exampleFiles/debug_data/" + mesh_name + "/results/qoi_basis_coeff.txt");
	for (auto coef = basis_error_coeffs.begin(); coef != basis_error_coeffs.end(); ++coef) {
		//one element per line: real imag
		file_qoi_out << (*coef).real() << " " << (*coef).imag() << std::endl;
	}
	std::cout << "Done saving QoI error per basis." << std::endl;
}