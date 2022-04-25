//Not being used ATM, b/c of weird include issues!


//
//#ifndef BASIS_H
//#define BASIS_H
//#include "../utility/constants.h"
//
//
//namespace basis {
//	void findPowers(int nu, int nv, int nw, int nglu, int nglv, int nglw, std::vector<double> xglu, std::vector<double> xglv, std::vector<double> xglw,
//		matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
//		matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
//		//upowers goes from (0 to nglu+1, 0 to nu)
//		//vpower (0 to nglv + 1, 0 to nv)
//		//wPowers (0 to nglw + 1, 0 to nw)
//		//fuPowers, fpuPowers same as uPowers
//		//fvPowers, fpvPowers same as vPowers
//		double x;
//		for (int i = 0; i <= nglu + 1; ++i) {
//			x = xglu[i];
//			uPowers(i, 0) = 1;
//			uPowers(i, 1) = x;
//			fuPowers(i, 0) = 1 - x;
//			fuPowers(i, 1) = 1 + x;
//			fpuPowers(i, 0) = -1;
//			fpuPowers(i, 1) = 1;
//			for (int j = 2; j <= nu; ++j) {
//				uPowers(i, j) = uPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fuPowers(i, j) = uPowers(i, j - 1)*x;
//					fpuPowers(i, j) = j*uPowers(i, j - 1);
//				}
//				else {
//					fuPowers(i, j) = uPowers(i, j) - x;
//					fpuPowers(i, j) = j*  uPowers(i, j - 1) - 1;
//				}
//			}
//		}
//		for (int i = 0; i <= nglv + 1; ++i) {
//			x = xglv[i];
//			vPowers(i, 0) = 1;
//			vPowers(i, 1) = x;
//			fvPowers(i, 0) = 1 - x;
//			fvPowers(i, 1) = 1 + x;
//			fpvPowers(i, 0) = -1;
//			for (int j = 2; j <= nv; ++j) {
//				vPowers(i, j) = vPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fvPowers(i, j) = vPowers(i, j) - 1;
//					fpvPowers(i, j) = j* vPowers(i, j - 1);
//				}
//				else {
//					fvPowers(i, j) = vPowers(i, j) - x;
//					fpvPowers(i, j) = j*vPowers(i, j - 1) - 1;
//				}
//			}
//		}
//		for (int i = 0; i <= nglw + 1; ++i) {
//			x = xglw[i];
//			wPowers(i, 0) = 1;
//			wPowers(i, 1) = x;
//			fwPowers(i, 0) = 1 - x;
//			fwPowers(i, 1) = 1 + x;
//			fpwPowers(i, 0) = -1;
//			fpwPowers(i, 1) = 1;
//			for (int j = 2; j <= nw; ++j) {
//				wPowers(i, j) = wPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fwPowers(i, j) = wPowers(i, j) - 1;
//					fpwPowers(i, j) = j*wPowers(i, j - 1);
//				}
//				else {
//					fwPowers(i, j) = wPowers(i, j) - x;
//					fpwPowers(i, j) = j*wPowers(i, j - 1) - 1;
//				}
//			}
//		}
//
//	}
//
//	//skipping legendre powers for now since it is not used in the example we have
//
//	//void findPowersLegendre(int nu, int nv, int nw, int nglu, int nglv, int glw, std::vector<double> xglu, std::vector<double> xglv, std::vector<double> xglw,
//	//	matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
//	//	matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
//	//	matrix2d<double> legendre; //starts at 0
//	//	legendre(0, 0) = 1; legendre(1, 1) = 1;
//	//	legendre(2, 0) = -.5; legendre(2, 2) = 1.5;
//	//	legendre(3, 1) = -1.5;
//	//}
//
//	void findPowersLagrange(int kuvw, int nglu, int nglv, int nglw, std::vector<double> xglu, std::vector<double> xglv, std::vector<double> xglw,
//		matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr,
//		matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers) {
//		//right now we have kuvw -- for now the order is the same in each direction -- this could be generalized later
//
//		double x, xj, xm, f;
//
//		for (int i = 0; i <= nglu + 1; i++) { // for u
//			x = xglu[i]; //Gauss-Legendre Integration point
//			for (int m = 0; m <= kuvw; m++) {
//				xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//				f = 1.0;
//				for (int j = 0; j <= kuvw; j++) {
//					xj = (2.0*j - kuvw) / kuvw;
//					if (j != m) {
//						f = f*(x - xj) / (xm - xj);
//					}
//				}//for j
//				fuPowersLagr(i, m) = f;
//			}//for m
//		}//for i
//
//		for (int i = 0; i <= nglv + 1; i++) { // for v
//			x = xglv[i]; //Gauss-Legendre Integration point
//			for (int m = 0; m <= kuvw; m++) {
//				xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//				f = 1.0;
//				for (int j = 0; j <= kuvw; j++) {
//					xj = (2.0*j - kuvw) / kuvw;
//					if (j != m) {
//						f = f*(x - xj) / (xm - xj);
//					}
//				}//for j
//				fvPowersLagr(i, m) = f;
//			}//for m
//		}//for i
//
//		for (int i = 0; i <= nglw + 1; i++) { // for w
//			x = xglw[i]; //Gauss-Legendre Integration point
//			for (int m = 0; m <= kuvw; m++) {
//				xm = (2.0*m - kuvw) / kuvw; //Lagrange Interpolation point
//				f = 1.0;
//				for (int j = 0; j <= kuvw; j++) {
//					xj = (2.0*j - kuvw) / kuvw;
//					if (j != m) {
//						f = f*(x - xj) / (xm - xj);
//					}
//				}//for j
//				fwPowersLagr(i, m) = f;
//			}//for m
//		}//for i
//
//	}//void findPowersLagrange
//
//	double trilagrangeoduvw(int m, int n, int l, int i, int j, int k, matrix2d<double> fuPowersLagr, matrix2d<double> fvPowersLagr, matrix2d<double> fwPowersLagr) {
//		return fuPowersLagr(m, i)*fvPowersLagr(n, j)*fwPowersLagr(l, k);
//	}
//
//
//	void find_eps_mu_matrix(int kuvw, int nglu, int nglv, int nglw, matrix2d<double>& fuPowersLagr, matrix2d<double>& fvPowersLagr, matrix2d<double>& fwPowersLagr,
//		int size1, std::vector<std::vector<std::complex<double>>>& muRel, matrix4d<std::complex<double>> muRelInt) {
//		std::complex<double> cMu;
//		for (auto s = 0; s < size1; ++s) { //matrix entries
//			for (int m = 0; m <= nglu + 1; ++m) {
//				for (int n = 0; n <= nglv + 1; ++n) {
//					for (int l = 0; l <= nglw + 1; ++l) {
//						cMu = std::complex<double>(0, 0);
//						for (int mat_loc = 0; mat_loc < muRel.size(); ++mat_loc) { //nodes
//							double f = trilagrangeoduvw(m, n, l, mat_loc % (kuvw + 1), int(floor(mat_loc / (kuvw + 1))) % (kuvw + 1),
//								int(floor(mat_loc / ((kuvw + 1)*(kuvw + 1)))), fuPowersLagr, fvPowersLagr, fwPowersLagr);
//
//							cMu += muRel[mat_loc][s] * f;
//						}
//						muRelInt(s, m, n, l) = cMu;
//					}
//				}
//
//			}
//
//		}
//	}
//
//};
//#endif //BASIS_H
