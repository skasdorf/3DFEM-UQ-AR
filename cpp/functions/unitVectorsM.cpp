#include "unitVectorsM.h"

namespace unitVectorsM {

	std::vector<double> findRA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		//all vectors start at 1
		std::vector<double> r(4);
		
		int ijk = 0;
		double unai, vnaj, wnak;
		for (int k = 0; k <= elementOrders[3]; ++k) {
			if (k > 0) wnak *= w;
			else wnak = 1.0;
			for (int j = 0; j <= elementOrders[2]; ++j) {
				if (j > 0) vnaj *= v;
				else vnaj = 1.0;
				for (int i = 0; i <= elementOrders[1]; ++i) {
					if (i > 0) unai *= u;
					else unai = 1.0;
					++ijk;
					for (int col = 1; col <= 3; ++col) {
						r[col] += rs(ijk, col) * unai*vnaj*wnak;
					}
				}
				
				
			}
		}
		return r;
	}

	std::vector<double> findR(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {

		std::vector<int> elementOrdersA = { 0, elementOrders, elementOrders,elementOrders };
		return unitVectorsM::findRA(elementOrdersA, u, v, w, rs, nRs);
	}

	std::vector<double> findAUA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		double unai, vnaj, wnak;
		std::vector<double> au(4);
		int ijk = 0;
		for (int k = 0; k <= elementOrders[3]; ++k) {
			if (k > 0) wnak *= w;
			else wnak = 1;
			for (int j = 0; j <= elementOrders[2]; ++j) {
				if (j > 0) vnaj *= v;
				else vnaj = 1;
				for (int i = 0; i <= elementOrders[1] - 1; ++i) {
					if (i > 0) unai *= u;
					else {
						unai = 1;
						++ijk;
					}
					++ijk;
					for (int col = 1; col <= 3; ++col) {
						au[col] += (i + 1)*rs(ijk, col) * unai*vnaj*wnak;
					}
				}
			}
		}
		return au;
	}

	std::vector<double> findAU(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		std::vector<int> elementOrdersA = { 0, elementOrders, elementOrders, elementOrders };
		return findAUA(elementOrdersA, u, v, w, rs, nRs);
	}

	std::vector<double> findAVA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		double unai, vnaj, wnak;
		int ijk = 0;
		std::vector<double> av(4);
		for (int k = 0; k <= elementOrders[3]; ++k) {
			if (k > 0) wnak *= w;
			else wnak = 1;
			for (int j = 0; j <= elementOrders[2] - 1; ++j) {
				if (j > 0) vnaj *= v;
				else {
					vnaj = 1;
					ijk += elementOrders[1] + 1;
				}
				for (int i = 0; i <= elementOrders[1]; ++i) {
					if (i > 0) unai *= u;
					else unai = 1;
					++ijk;
					for (int col = 1; col <= 3; ++col) {
						av[col] += (j + 1)*rs(ijk,col) * unai*vnaj*wnak;
					}
				}
			}
		}
		return av;
	}

	std::vector<double> findAV(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		std::vector<int> elementOrdersA = { 0, elementOrders, elementOrders, elementOrders };
		return findAVA(elementOrdersA, u, v, w, rs, nRs);
	}
	std::vector<double> findAWA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		std::vector<double> aw(4);
		double unai, vnaj, wnak;
		int ijk = (elementOrders[1] + 1)*(elementOrders[2] + 1);
		for (int k = 0; k <= elementOrders[3] - 1; ++k) {
			if (k > 0) wnak *= w;
			else wnak = 1;
			for (int j = 0; j <= elementOrders[2]; ++j) {
				if (j > 0) vnaj *= v;
				else vnaj = 1;
				for (int i = 0; i <= elementOrders[1]; ++i) {
					if (i > 0) unai *= u;
					else unai = 1;
					++ijk;
					for (int col = 1; col <= 3; ++col) {
						aw[col] += (k + 1)*rs(ijk,col) * unai*vnaj*wnak;
					}
				}
			}
		}
		return aw;
	}
	std::vector<double> findAW(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs) {
		std::vector<int> elementOrdersA = { 0, elementOrders, elementOrders, elementOrders };
		return findAWA(elementOrdersA, u, v, w, rs, nRs);
	}
	void unitaryvectors(const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, 
		const std::vector<double>& xglw, const int& elementOrders, const matrix2d<double>& rs, matrix4d<double>& rMatrix, 
		matrix4d<double>& auMatrix, matrix4d<double>& avMatrix, matrix4d<double>& awMatrix, matrix3d<double>& jacobianMatrix, const int& nRs, const int& integralType) {
		// Matrix(0:nglu+1,0:nglv+1,0:nglw+1,3) for 4d matrices (i.e. last dim is length of 4 since it starts index at 1, others start at 0
		//jacobianMatrix(0:nglu+1,0:nglv+1,0:nglw+1) all dimensions start at 0 for jacobian
		std::vector<double> r;
		std::vector<double> au;
		std::vector<double> av;
		std::vector<double> aw;
		int face;
		double jacobian;
		for (int m = 0; m <= nglu + 1; ++m) {
			for (int n = 0; n <= nglv + 1; ++n) {
				for (int l = 0; l <= nglw + 1; ++l) {
					r = unitVectorsM::findR(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
					rMatrix(m, n, l) = r;
				}
			}
		}

		for (int m = 0; m <= nglu + 1; ++m) {
			for (int n = 0; n <= nglv + 1; ++n) {
				for (int l = 0; l <= nglw + 1; ++l) {
					au = findAU(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
					auMatrix(m, n, l, 1) = au[1];
					auMatrix(m, n, l, 2) = au[2];
					auMatrix(m, n, l, 3) = au[3];

				}
			}
		}
		
		for (int m = 0; m <= nglu + 1; ++m) {
			for (int n = 0; n <= nglv + 1; ++n) {
				for (int l = 0; l <= nglw + 1; ++l) {
					av = findAV(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
					avMatrix(m, n, l, 1) = av[1];
					avMatrix(m, n, l, 2) = av[2];
					avMatrix(m, n, l, 3) = av[3];
				}
			}
		}
		for (int m = 0; m <= nglu + 1; ++m) {
			for (int n = 0; n <= nglv + 1; ++n) {
				for (int l = 0; l <= nglw + 1; ++l) {
					aw = findAW(elementOrders, xglu[m], xglv[n], xglw[l], rs, nRs);
					awMatrix(m, n, l, 1) = aw[1];
					awMatrix(m, n, l, 2) = aw[2];
					awMatrix(m, n, l, 3) = aw[3];
				}
			}
		}
		switch (integralType) {
		case 1: //surface integral
			if (nglu == 1) face = 1;
			else if (nglv == 1) face = 2;
			else if (nglw == 1) face = 3;
			else exit(EXIT_FAILURE);
			for (int m = 0; m <= nglu + 1; ++m) {
				for (int n = 0; n <= nglv + 1; ++n) {
					for (int l = 0; l <= nglw + 1; ++l) {
						for (int coord = 1; coord <= 3; ++coord) {
							au[coord] = auMatrix(m, n, l, coord);
							av[coord] = avMatrix(m, n, l, coord);
							aw[coord] = awMatrix(m, n, l, coord);
						}
						//jacobian = find_jacobian_surface(au, av, aw, face);
						//jacobianMatrix(m, n, l) = jacobian;
						find_jacobian_surface(au, av, aw, face, jacobianMatrix(m,n,l));
					//	jacobianMatrix(m, n, l) = jacobian;

					}
				}
			}
			break;
		case 2:
			for (int m = 0; m <= nglu + 1; ++m) {
				for (int n = 0; n <= nglv + 1; ++n) {
					for (int l = 0; l <= nglw + 1; ++l) {
						for (int coord = 1; coord <= 3; ++coord) {
							au[coord] = auMatrix(m, n, l, coord);
							av[coord] = avMatrix(m, n, l, coord);
							aw[coord] = awMatrix(m, n, l, coord);
						}
						//jacobian = find_jacobian(au, av, aw);
						//jacobianMatrix(m, n, l) = jacobian;
						find_jacobian(au, av, aw, jacobianMatrix(m, n, l));
					}
				}
			}
		}

	}

	void find_jacobian(const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw, double& jacobian) {
		std::vector<double> vec(4);
		products::cross(au, av, vec);
		products::dot(vec, aw, jacobian); //returns jacobian
	}
	void find_jacobian_surface(const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw, const int& face, double& result) {
		std::vector<double> vec1(4);
		//std::vector<double> vec2;
		//double jacobianTemp;

		switch (face) {
		case 1:
			//vec1 = products::cross(av, aw);
			products::cross(av, aw, vec1);
			break;
		case 2:
			//vec1 = products::cross(au, aw);
			products::cross(au, aw, vec1);
			break;
		case 3:
		//	vec1 = products::cross(au, av);
			products::cross(au, av, vec1);
			break;
		}
		//vec2 = vec1;
		products::dot(vec1, vec1, result);
		result = std::sqrt(result);
	}
};