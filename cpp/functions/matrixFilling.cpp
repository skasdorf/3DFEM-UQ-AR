//#include "matrixFilling.h"
//
//
//
//namespace matrixFilling {
//	void findelements_abs_sparse(Domain& dom) {
//		bool isAdjoint = false; //get from domain class
//		matrix2d<double> fuPowersLagr;
//		matrix2d<double> fvPowersLagr;
//		matrix2d<double> fwPowersLagr;
//		for (auto e = elements.begin(); e != elements.end(); ++e) {
//			int nu = e->expansion[0];
//			int nv = e->expansion[1];
//			int nw = e->expansion[2];
//			int nglu = e->quadrature[0];
//			int nglv = e->quadrature[1];
//			int nglw = e->quadrature[2];
//			int iHomCode = e->materials.hcode;
//			int aCode = e->materials.icode;
//			std::complex<double> cEps;
//			std::complex<double> cMu;
//			int kuvw;
//			int size1;
//			if (aCode == 0) {
//				kuvw = e->kuvwA;
//				size1 = 9 - 3 * e->materials.sym;
//			}
//			else {
//				kuvw = e->kuvw;
//				size1 = 1;
//			}
//			std::vector<double> xglu, xglv, xglw; //coords
//			std::vector<double> wglu, wglv, wglw; //weights
//			functions::gaussk(nglu, xglu, wglu);
//			functions::gaussk(nglv, xglv, wglv);
//			functions::gaussk(nglw, xglw, wglw);
//			matrix4d<double> rMatrix;
//			matrix4d<double> auMatrix;
//			matrix4d<double> avMatrix;
//			matrix4d<double> awMatrix;
//			//matrix3d<double>& jacobianMatrix
//			matrix2d<double> rs;
//			unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
//				e->geom_order, e->rs,rMatrix, auMatrix, avMatrix, awMatrix, 
//				e->jacobian, e->nRs, 2);
//			int indicatorBasisType = 1; //kolund basis
//			matrix2d<double> uPowers, vPowers, wPowers, fuPowers, fvPowers, fwPowers, fpuPowers, fpvPowers, fpwPowers;
//			 //if indicator = 1, then kolun, if indicator = 2, then legendre
//			if (indicatorBasisType == 1)
//				//kolund basis functions
////				basis::findPowers(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw, uPowers, vPowers, wPowers, fuPowers, fvPowers, fwPowers, fpuPowers, fpvPowers, fpwPowers);
//			//else if (indicatorBasisType == 2)
//				//not currently implemented
//				//legendre
//			//	basis::findPowersLegendre(nu, nv, nw, nglu, nglv, xglu, xglv, xglw, uPowers)
//			std::complex<double> EpsRelInt = cEps;
//			std::complex<double> MuRelInt = cMu;
//			if ((iHomCode == 1) && (aCode == 1)) {
//				cEps = e->materials.epsr_list[0][0]; //homogen
//				cMu = e->materials.mur_list[0][0];//homogen
//				
//				if (isAdjoint == true) {
//					cMu = std::conj(cMu);
//					cEps = std::conj(cEps);
////					EpsRelInt = std::conj(EpsRelInt);
//					MuRelInt = std::conj(MuRelInt);
//				}
//			}
//			else {
//				if (kuvw != 0) { //make sure that we acatully want this to be the case, since currenty k=-1 is the default case
//					//basis::findPowersLagrange(kuvw, nglu, nglv, nglw, xglu, xglv, xglw, fuPowersLagr, fvPowersLagr, fwPowersLagr);
//					//basis::find_eps_mu_matrix(kuvw, nglu, nglv, nglw, fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e->materials.epsilonRel, EpsRelInt);
//					//basis::find_eps_mu_matrix(kuvw, nglu, nglv, nglw, fuPowersLagr, fvPowersLagr, fwPowersLagr, size1, e->materials.muRel, MuRelInt);
//					if (aCode == 0) {
//						if (iHomCode == 0) {
//		//					functions::find_mur_inv(size1, MuRelInt, nglu, nglv, nglw, MuRelIntInv);
//							if (isAdjoint) {
//								//do nothing for now
//								int nothing = 0;
//							}
//						}
//						else {
//			//				functions::find_mur_inv(size1, MuRelInt, iOne, iOne, iOne, MuRelIntInv);
//							if (isAdjoint) {
//								//do nothing for now
//								int nothing = 0;
//							}
//						}
//					}
//					iHomCode = 1; //homogen, anisotro
//				//	EpsRelInt = e->materials.epsilonRel;
//				//	MuRelInt = e->materials.muRel;
//				//	functions::find_mur_inv(size1, MuRelInt, 1, 1, 1, MuRelIntInv);
//					if (isAdjoint) {
//						//do nothing for now
//						int nothing = 0;
//					}
//				}
//			}
//			int iCon;
//			for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
//				iCon = abs(vectorD[kMat]);
//				if (iCon == 0) continue;
//				int eh = eiUVWijk[kMat][1];
//				int iuvwh = eiUVWijk[kMat][2];
//				int ih = eiUVWijk[kMat][3];
//				int jh = eiUVWijk[kMat][4];
//				int kh = eiUVWijk[kMat][5];
//				if (isAdjoint) {
//					//do nothing for now
//					int nothing = 0;
//				}
//				else {
//					int face = 1;
//					if (e->region == 1) {
//						Integral_g gInts;
//					//	gInts.findGWave(dom, eh, iuvwh, ih, jh, kh, scatter1.nWaves, scatter1.nF, );
//					}
//				}
//				for (int lMat = kMat; lMat <= e->unknownsEnd; ++lMat) {
//					int jCon = abs(vectorD[lMat]);
//					if (jCon == 0) continue;
//					int e = eiUVWijk[lMat][1];
//					int iuvw = eiUVWijk[lMat][2];
//					int i = eiUVWijk[lMat][3];
//					int j = eiUVWijk[lMat][4];
//					int k = eiUVWijk[lMat][5];
//					if (eh != e) {
//						exit(EXIT_FAILURE);
//					}
//					//cSolStotal = std::complex<double>(0, 0);
//
//				}
//			}
//		}
//	}
//};