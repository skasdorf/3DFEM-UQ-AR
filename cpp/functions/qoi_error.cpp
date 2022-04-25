//#include "qoi_error.h"
//
//void qoi::qoi_error(Domain & dom_h, Domain & dom_l, std::string mesh_name)
//{
//	//read in cAlpha (lower), cAlpha (adjoint), and cGr (forward - higher)
//	
//	std::vector < std::complex<double>> cGr_higher, cAlphaFor, cAlphaAdj, cAlphaFormod;
//	functions::read_previous_solve(mesh_name, cAlphaFor, cAlphaAdj, cGr_higher);
//	double k0 = dom_h.scatter1.K0[1] * dom_h.scatter1.K0[1];
//
//	for (auto e = dom_h.elements.begin(); e != dom_h.elements.end(); ++e) {
//		std::cout << "For element " << e->index << std::endl;
//		int nu = e->expansion[0];
//		int nv = e->expansion[1];
//		int nw = e->expansion[2];
//		int nglu = e->quadrature[0];
//		int nglv = e->quadrature[1];
//		int nglw = e->quadrature[2];
//		int iHomCode = e->materials.hcode;
//		int aCode = e->materials.icode;
//		e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
//		std::complex<double> cEps;
//		std::complex<double> cMu;
//		bool isAdjoint = false;
//		int kuvw;
//		int size1;
//		if (aCode == 0) {
//			kuvw = e->materials.KuvwA;
//			size1 = 9 - 3 * e->materials.sym;
//		}
//		else {
//			kuvw = e->materials.Kuvw;
//			size1 = 1;
//		}
//		std::vector<double> xglu, xglv, xglw; //coords
//
//		functions::gaussk(nglu, xglu, e->wglu);
//		functions::gaussk(nglv, xglv, e->wglv);
//		functions::gaussk(nglw, xglw, e->wglw);
//		e->rMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
//		e->auMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
//		e->avMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
//		e->awMatrix = matrix4d<double>(nglu + 1, nglv + 1, nglw + 1, 3);
//
//
//		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
//			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
//			e->jacobian, e->nRs, 2);
//		int indicatorBasisType = 1; //kolund basis
//
//		e->uPowers = matrix2d<double>(nglu + 1, nu);
//		e->vPowers = matrix2d<double>(nglv + 1, nv);
//		e->wPowers = matrix2d<double>(nglw + 1, nw);
//		e->fuPowers = matrix2d<double>(nglu + 1, nu);
//		e->fvPowers = matrix2d<double>(nglv + 1, nv);
//		e->fwPowers = matrix2d<double>(nglw + 1, nw);
//		e->fpuPowers = matrix2d<double>(nglu + 1, nu);
//		e->fpvPowers = matrix2d<double>(nglv + 1, nv);
//		e->fpwPowers = matrix2d<double>(nglw + 1, nw);
//		e->fuPowersLagr = matrix2d<double>(nglu + 1, kuvw);
//		e->fvPowersLagr = matrix2d<double>(nglv + 1, kuvw);
//		e->fwPowersLagr = matrix2d<double>(nglw + 1, kuvw);
//
//		//if indicator = 1, then kolun, if indicator = 2, then legendre
//		if (indicatorBasisType == 1)
//			//kolund basis functions
//			dom_h.findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, e->uPowers, e->vPowers, e->wPowers, e->fuPowers, e->fvPowers, e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);
//
//
//		e->MuRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
//		e->EpsRelInt = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
//		e->MuRelIntInv = matrix4d<std::complex<double>>(size1, nglu + 1, nglv + 1, nglw + 1);
//
//
//		if ((iHomCode == 1) && (aCode == 1)) {
//			cEps = e->materials.epsr_list[0][0]; //homogen
//			cMu = e->materials.mur_list[0][0];//homogen
//			e->MuRelInt(1, 1, 1, 1) = cMu;
//			e->EpsRelInt(1, 1, 1, 1) = cEps;
//			if (isAdjoint == true) {
//
//				cMu.imag(-cMu.imag());
//				cEps.imag(-cEps.imag());
//				e->EpsRelInt(1, 1, 1, 1).imag(e->EpsRelInt(1, 1, 1, 1).imag());
//				e->MuRelInt(1, 1, 1, 1).imag(e->MuRelInt(1, 1, 1, 1).imag());
//			}
//		}
//		else {
//			if (kuvw != -1) {
//				findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr);
//
//				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.epsr_list, e->EpsRelInt);
//				find_eps_mu_matrix(kuvw, nglu, nglv, nglw, e->fuPowersLagr, e->fvPowersLagr, e->fwPowersLagr, size1, e->materials.mur_list, e->MuRelInt);
//
//				if (aCode == 0) {
//					if (iHomCode == 0) {
//						functions::find_mur_inv(size1, e->MuRelInt, nglu, nglv, nglw, e->MuRelIntInv);
//						if (isAdjoint) {
//							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->MuRelIntInv);
//							find_eps_mu_hermitian(size1, nglu, nglv, nglw, e->EpsRelInt);
//						}
//					}
//					else {
//						functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
//						if (isAdjoint) {
//							find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
//							find_eps_mu_hermitian(size1, 1, 1, 1, e->EpsRelInt);
//						}
//					}
//				}
//				else {
//					iHomCode = 1; //homogeneous, anisotropic
//					e->EpsRelInt(1, 1, 1, 1) = e->materials.epsr_list[0][0];
//					e->MuRelInt(1, 1, 1, 1) = e->materials.mur_list[0][0];
//					functions::find_mur_inv(size1, e->MuRelInt, 1, 1, 1, e->MuRelIntInv);
//					if (isAdjoint) {
//						find_eps_mu_hermitian(size1, 1, 1, 1, e->MuRelIntInv);
//						e->EpsRelInt(1, 1, 1, 1).imag(-e->EpsRelInt(1, 1, 1, 1).imag());
//						e->MuRelInt(1, 1, 1, 1).imag(-e->MuRelInt(1, 1, 1, 1).imag());
//					}
//				}
//
//
//			}
//		}
//
//		Integral_g IntG(&(*e));
//		Integral_c Intc(&(*e));
//		Integral_d Intd(&(*e));
//	}
//}
