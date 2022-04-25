
#include <iostream>
#include "functions.h"

namespace functions {

	std::vector<std::string> split(const std::string &s, char delim) {
		std::stringstream ss(s);
		std::string item;
		std::vector<std::string> tokens;
		while (std::getline(ss, item, delim)) {
		//	if (item.find_first_of(" ") != std::string::npos) continue;
			if (item.empty()) continue; //kills spaces
			tokens.push_back(item);
		}
		return tokens;
	}

	void MIMAISWAP(int * a, int * b) {
		int temp;
		temp = *a;
		*a = *b;
		*b = temp;
	}
	bool _debug_matrix(std::string filename, std::vector<int>& matrix, int start) {
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			std::getline(file, line);
			std::vector<std::string> input_line = functions::split(line, ' ');
			if (matrix[i] != std::stoi(input_line[0])) {
				return false;
			}
		}
		return true;
	}
	bool _debug_matrix2(std::string filename, std::vector<int>& matrix, int start) {
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			std::getline(file, line);
			std::vector<std::string> input_line = functions::split(line, ' ');
			if (matrix[i] != std::stoi(input_line[0])) {
				return false;
			}
		}
		return true;
	}
	bool _debug_matrix(std::string filename, std::vector<std::vector<int>>& matrix, int start) {
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			for (int j = start; j < matrix[start].size(); ++j) {
				std::getline(file, line);
				std::vector<std::string> input_line = functions::split(line, ' ');
				if (matrix[i][j] != std::stoi(input_line[0])) { 
					int stop = 0;
					return false; }
			}
		}
		return true;
	}
	bool _debug_matrix(std::string filename, std::vector<std::vector<std::vector<int>>>& matrix, int start) {
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			for (int j = start; j < matrix[start].size(); ++j) {
				for (int k = start; k < matrix[start][start].size(); ++k) {
					std::getline(file, line);
					std::vector<std::string> input_line = functions::split(line, ' ');
					if (matrix[i][j][k] != std::stoi(input_line[0])) { 
						return false;
					}
				}
			}
		}
		return true;
	}
	bool _debug_matrix(std::string filename, std::vector<double>& matrix, int start) {
		double eps = 0.000001;
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			std::getline(file, line);
			std::vector<std::string> input_line = functions::split(line, ' ');
			if (abs(matrix[i] - std::stod(input_line[0])) > eps) {
				int stop = 0;
				return false;
			}
		}
		return true;
	}
	bool _debug_matrix(std::string filename, std::vector<std::vector<double>>& matrix, int start) {
		double eps = 1e-7;
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			for (int j = start; j < matrix[start].size(); ++j) {
				std::getline(file, line);
				std::vector<std::string> input_line = functions::split(line, ' ');
				if (abs(matrix[i][j] - std::stod(input_line[0])) > eps) {
					int stop = 0;
					return false;
				}
			}
		}
		return true;
	}
	bool _debug_matrix(std::string filename, std::vector<std::vector<std::vector<double>>>& matrix, int start) {
		double eps = 1e-7;
		std::string line;
		std::ifstream file;
		file.open(filename);
		for (int i = start; i < matrix.size(); ++i) {
			for (int j = start; j < matrix[start].size(); ++j) {
				for (int k = start; k < matrix[start][start].size(); ++k) {
					std::getline(file, line);
					std::vector<std::string> input_line = functions::split(line, ' ');
					if (abs(matrix[i][j][k] - std::stod(input_line[0])) > eps) { 
						return false; }
				}
			}
		}
		return true;
	}
	void gaussk(int n, std::vector<double> &x, std::vector<double>& w) {
		//x starts at 0 goes to n+1
		//w starts at 1 goes to n
		//x contains the coordinates (in one direction (u,v, or w))
		//w contains the weights for those pts
		x.clear();
		x.resize(n + 2);
		w.clear();
		w.resize(n + 1);
		int it, k, m, nd;
		double fx, h, p, pk, pkm1, pkp1, u, v;
		double t, x0, t1, den, d1, dpn, d2pn, d3pn, d4pn, dp;
	
		if (n > 100) n = 100;
		m = (n + 1) / 2;
		double e1 = n*(n + 1);
		for (int i = 1; i <= m; ++i) {
			t = (4 * i - 1)*6.28318530717958624 / (8 * n + 4);
			x0 = (1.0 - (1.0 - 1.0 / n) / (8.0 * n*n));
			x0=x0*std::cos(t);
			pkm1 = 1.0;
			pk = x0;
			for (int k = 2; k <= n; ++k) {
				t1 = x0*pk;
				pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1;
				pkm1 = pk;
				pk = pkp1;
			}
			den = 1.0 - x0*x0;
			d1 = n*(pkm1 - x0*pk);
			dpn = d1 / den;
			d2pn = (2.0 * x0*dpn - e1*pk) / den;
			d3pn = (4.0 * x0*d2pn + (2.0 - e1)*dpn) / den;
			d4pn = (6.0 * x0*d3pn + (6.0 - e1)*d2pn) / den;
			u = pk / dpn;
			v = d2pn / dpn;
			h = -u*(1 + .5*u*(v + u*(v*v - u*d3pn / (3.0 * dpn))));
			p = pk + h*(dpn + .5*h*(d2pn + h / 3.0 * (d3pn + .25*h*d4pn)));
			dp = dpn + h*(d2pn + .5*h*(d3pn + h*d4pn / 3.0));
			h = h - p / dp;
			x[i] = x0 + h;
			fx = d1 - h*e1*(pk + .5*h*(dpn + h / 3 * (d2pn + .25*h*(d3pn + .2*h*d4pn))));
			w[i] = 2.0 * (1.0 - x[i] * x[i]) / (fx*fx);
		}
		if (2*m > n) x[m] = 0;
		nd = (n + 1) / 2;
		it = n;
		for (int i = 1; i <= nd; ++i) {
			x[it] = x[i];
			x[i] = -x[it];
			w[it] = w[i];
			w[i] = w[it];
			it = it - 1;
		}
		x[0] = -1.0;
		x[n + 1] = 1.0;
	}
	void find_mur_inv(const int& size1, const matrix4d<std::complex<double>>& MuRel, const int& nglu, const int& nglv, const int& nglw, matrix4d<std::complex<double>>& MuRelInv) {
		std::complex<double> det, a, b, c, d, e, f, g, h, i;
		std::complex<double> A1, B1, C1, D1, E1, F1, G1, H1, K1;

		for (int k = 1; k <= nglu; ++k) {
			for (int l = 1; l <= nglv; ++l) {
				for (int m = 1; m <= nglw; ++m) {
					if (size1 == 3) {//only diag matrix
						if ((MuRelInv(1, k, l, m) == std::complex<double>(0, 0))) exit(EXIT_FAILURE); //non invertible
						else {
							MuRelInv(1, k, l, m) = double(1) / MuRel(1, k, l, m);
							MuRelInv(2, k, l, m) = double(1) / MuRel(2, k, l, m);
							MuRelInv(3, k, l, m) = double(1) / MuRel(3, k, l, m);
						}
					}
					else {
						a = MuRel(1, k, l, m);
						b = MuRel(2, k, l, m);
						c = MuRel(3, k, l, m);
						i = MuRel(size1, k, l, m);
						if (size1 == 6) {
							d = MuRel(2, k, l, m);
							e = MuRel(4, k, l, m);
							f = MuRel(5, k, l, m);
							g = MuRel(3, k, l, m);
							h = MuRel(5, k, l, m);
						}
						else if (size1 == 9) {
							d = MuRel(4, k, l, m);
							e = MuRel(5, k, l, m);
							f = MuRel(6, k, l, m);
							g = MuRel(7, k, l, m);
							h = MuRel(8, k, l, m);
						}
						det = a*(e*i - f*h) + b*(f*g - i*d) + c*(d*h - e*g);
						if (det == std::complex<double>(0, 0)) exit(EXIT_FAILURE); //non invertible matrix
						A1 = (e*i - f*h);
						B1 = (f*g - d*i);
						C1 = (d*h - e*g);
						E1 = (a*i - c*g);
						F1 = (g*b - a*h);
						K1 = (a*e - b*d);
						if (size1 == 9) {
							D1 = (c*h - b*i);
							G1 = (b*f - c*e);
							H1 = (c*d - a*f);
						}
						MuRelInv(1, k, l, m) = A1 / det;
						MuRelInv(2, k, l, m) = B1 / det;
						MuRelInv(3, k, l, m) = C1 / det;
						MuRelInv(size1, k, l, m) = K1 / det;
						if (size1 == 6) {
							MuRelInv(4, k, l, m) = E1 / det;
							MuRelInv(5, k, l, m) = F1 / det;
						}
						else if (size1 == 9) {
							MuRelInv(4, k, l, m) = D1 / det;
							MuRelInv(5, k, l, m) = E1 / det;
							MuRelInv(6, k, l, m) = F1 / det;
							MuRelInv(7, k, l, m) = G1 / det;
							MuRelInv(8, k, l, m) = H1 / det;
						}
					}
				}
			}
		}
	}
	double trilagrangeoduvw(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, const matrix2d<double>& fwPowersLagr, const int& nglu, const int& nglv, const int& nglw) {
		return fuPowersLagr(m, i)*fvPowersLagr(n, j)*fwPowersLagr(l, k);
	}
}

	void functions::read_previous_solve(std::string mesh_name, std::vector<std::complex<double>>& cAlphaForLower, std::vector<std::complex<double>>& cAlphaAdjoint, std::vector<std::complex<double>>& cGrhigher){
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
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ',');
			cAlphaAdjoint.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1])));

		}
		file.close();
		file.open("../exampleFiles/debug_data/" + mesh_name + "/results/cGr_higher_order.txt"); //higher order for cGr from forwars solve
	
		while (getline(file, line)) {
			line = line.substr(1, line.size() - 2);
			line_input = functions::split(line, ',');
			cGrhigher.push_back(std::complex<double>(std::stod(line_input[0]), std::stod(line_input[1]))); //only works for one freq atm

		}
		file.close();

	}