//#include "../structure/Element.h"
//#include "../structure/Scatter.h"
//#include "../structure/Facet.h"
//#include "../functions/unitVectorsM.h"
//#include "../integrals/Integral_g.h"
//
//namespace sensitivity {
//
//	void findPowersSens(const int& nu, const int& nv, const int& nw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
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
//			uPowers(i, 0) = 1.0;
//			uPowers(i, 1) = x;
//			fuPowers(i, 0) = 1.0 - x;
//			fuPowers(i, 1) = 1.0 + x;
//			fpuPowers(i, 0) = -1.0;
//			fpuPowers(i, 1) = 1.0;
//			for (int j = 2; j <= nu; ++j) {
//				uPowers(i, j) = uPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fuPowers(i, j) = uPowers(i, j) - 1.0;
//					fpuPowers(i, j) = j*uPowers(i, j - 1);
//				}
//				else {
//					fuPowers(i, j) = uPowers(i, j) - x;
//					fpuPowers(i, j) = j*  uPowers(i, j - 1) - 1.0;
//				}
//			}
//		}
//		for (int i = 0; i <= nglv + 1; ++i) {
//			x = xglv[i];
//			vPowers(i, 0) = 1.0;
//			vPowers(i, 1) = x;
//			fvPowers(i, 0) = 1.0 - x;
//			fvPowers(i, 1) = 1.0 + x;
//			fpvPowers(i, 0) = -1.0;
//			fpvPowers(i, 1) = 1.0;
//			for (int j = 2; j <= nv; ++j) {
//				vPowers(i, j) = vPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fvPowers(i, j) = vPowers(i, j) - 1.0;
//					fpvPowers(i, j) = j* vPowers(i, j - 1);
//				}
//				else {
//					fvPowers(i, j) = vPowers(i, j) - x;
//					fpvPowers(i, j) = j*vPowers(i, j - 1) - 1.0;
//				}
//			}
//		}
//		for (int i = 0; i <= nglw + 1; ++i) {
//			x = xglw[i];
//			wPowers(i, 0) = 1.0;
//			wPowers(i, 1) = x;
//			fwPowers(i, 0) = 1.0 - x;
//			fwPowers(i, 1) = 1.0 + x;
//			fpwPowers(i, 0) = -1.0;
//			fpwPowers(i, 1) = 1.0;
//			for (int j = 2; j <= nw; ++j) {
//				wPowers(i, j) = wPowers(i, j - 1)*x;
//				if (j % 2 == 0) {
//					fwPowers(i, j) = wPowers(i, j) - 1.0;
//					fpwPowers(i, j) = j*wPowers(i, j - 1);
//				}
//				else {
//					fwPowers(i, j) = wPowers(i, j) - x;
//					fpwPowers(i, j) = j*wPowers(i, j - 1) - 1.0;
//				}
//			}
//		}
//
//	}
//	
//	void spherical_sensitivity(const int& myElementCount, std::vector<Element>& elements, const int& matDimDis,
//		const std::vector<int>& vectorD, const std::vector<std::vector<int>>& eiUVWijk, const int* matDimCon, const std::vector<std::complex<double>>& cAlphaForward,
//		std::vector<std::complex<double>>& cAlphaAdjoint, std::complex<double>& cSensitivity, const int& numberOfFrequencies,  
//		std::vector<double>& k0List, const int& indicatorBasisType, Scatter* scatter1, const std::vector<Facet>& facets) {
//		//note: cAlphaForward == forward solve results (potentially read in from file)
//		//note: cAlphaAdjoint == adjoint solve results read in from file (potentially just computed) 
//		std::complex<double> cj(0, 1);
//		cSensitivity = 0;
//		int kuvw = 0;
//		int size1 = 1;
//		int iHomCode = 1;
//		int aCode = 1;
//		int integralType = 1;
//		//currently only works with one single wave
//		int waveNumber = 1;
//		int freqNo = 1;
//		double val1 = -scatter1->K0[freqNo];
//		auto cEteta = scatter1->waveTheta[waveNumber];
//		auto cEfi = scatter1->wavePhi[waveNumber];
//		auto iFiAr = scatter1->iPhiAr[waveNumber];
//		auto iThetaAr = scatter1->iThetaAr[waveNumber];
//		auto nAr = scatter1->nAr[waveNumber];
//		std::vector<double> iR(4);
//		std::vector<std::complex<double>> asec2(4), eVector(4), hVector(4), asec3(4);
//		asec2[1] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][1] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][1]);
//		asec2[2] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][2] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][2]);
//		asec2[3] = (scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][3] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][3]);
//		eVector[1] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][1] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][1];
//		eVector[2] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][2] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][2];
//		eVector[3] = scatter1->wavePhi[waveNumber] * scatter1->iPhiAr[waveNumber][3] + scatter1->waveTheta[waveNumber] * scatter1->iThetaAr[waveNumber][3];
//		products::cross(scatter1->nAr[waveNumber], eVector, hVector);
//		asec3[1] = -cj*hVector[1];
//		asec3[2] = -cj*hVector[2];
//		asec3[3] = -cj*hVector[3];
//		iR[1] = -scatter1->nAr[waveNumber][1];
//		iR[2] = -scatter1->nAr[waveNumber][2];
//		iR[3] = -scatter1->nAr[waveNumber][3];
//		std::complex<double> cEpsilon;
//		std::complex<double> cMu;
//		std::complex<double> cSolB, cSolGTotal;
//		for (auto e = elements.begin(); e != elements.end(); ++e) {
//			Integral_g gInts(&(*e));
//			
//			int nu = e->expansion[0];
//			int nv = e->expansion[1];
//			int nw = e->expansion[2];
//			//get the faces of the element that are connected to another element
//			
//			for (auto face = e->facet_indices.begin(); face != e->facet_indices.end(); ++face) {
//				int connectedTo; //the index of the other element
//				int connectedToType;
//				Facet this_facet = facets[*face];
//				int iFace = this_facet.index;
//				this_facet.element_indices.first == e->index ? connectedTo = this_facet.element_indices.second 
//					: connectedTo = this_facet.element_indices.first;
//				
//				//if not connected directly to another elements
//				if (connectedTo != 0) { //might need to -1, instead of 0
//					connectedToType = elements[connectedTo].materials.region;
//					if (connectedToType == 0 && e->materials.region == 1) { //might need to change 0 and 1 to be different
//						continue;
//					}
//				}
//				else {
//					continue;
//				}
//				int nglu = e->quadrature[0];
//				int nglv = e->quadrature[1];
//				int nglw = e->quadrature[2];
//				if (iFace == 1 || iFace == 2) nglu = 1;
//				else if (iFace == 3 || iFace == 4) nglv = 1;
//				else if (iFace == 5 || iFace == 6) nglw == 1;
//				//call new G-L integration params
//				std::vector<double> xglu, xglv, xglw, wglu, wglv, wglw;
//				functions::gaussk(nglu, xglu, wglu);
//				functions::gaussk(nglv, xglv, wglv);
//				functions::gaussk(nglw, xglw, wglw);
//				//modify the xglu,v,w vectors
//				switch (iFace){
//				case 1:
//					xglu[1] = -1;
//					break;
//				case 2:
//					xglu[1] = 1;
//					break;
//				case 3:
//					xglv[1] = -1;
//					break;
//				case 4:
//					xglv[1] = 1;
//					break;
//				case 5:
//					xglw[1] = -1;
//					break;
//				case 6:
//					xglw[1] = 1;
//					break;
//				}
//				unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw, 
//					e->geom_order, e->rs, e->rMatrix, e->auMatrix,e->avMatrix,
//					e->awMatrix,e->jacobian,e->nRs, integralType);
//				switch (indicatorBasisType) { //current set to only do kolundzija
//				case 1:
//					findPowersSens(nu, nv, nw, nglu, nglv, nglw, xglu, xglv, xglw,e->uPowers, 
//						e->vPowers, e->wPowers,e->fuPowers, e->fvPowers, 
//						e->fwPowers, e->fpuPowers, e->fpvPowers, e->fpwPowers);
//					break;
//				case 2:
//					//do legendre powers
//					break;
//				}
//				if (iHomCode == 1 && aCode == 1) {
//					cMu = std::complex<double>(1.0, 0.0);
//					cEpsilon = e->materials.epsr_list[0][0]; //might need to start at 1
//					e->MuRelInt(1, 1, 1, 1) = std::complex<double>(1.0, 0.0);
//					e->EpsRelInt(1, 1, 1, 1) = cEpsilon;
//				}
//				for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
//					int iCon = abs(vectorD[kMat]);
//					if (iCon == 0) continue;
//					int eh = eiUVWijk[kMat][1];
//					int iuvwh = eiUVWijk[kMat][2];
//					int ih = eiUVWijk[kMat][3];
//					int jh = eiUVWijk[kMat][4];
//					int kh = eiUVWijk[kMat][5];
//					gInts.set_e(&elements[eh - 1]);
//					cSolGTotal = 0;
//					gInts.findGInc_EpsOnly(scatter1, vectorD, iuvwh, ih, jh, kh, size1,k0List[freqNo], freqNo, cSolGTotal, asec2, hVector, asec3, val1);
//					cSensitivity += cSolGTotal*std::conj(cAlphaAdjoint[iCon]);
//					for (int lMat = kMat; lMat <= e->unknownsEnd; ++lMat) {
//						int jCon = abs(vectorD[lMat]);
//						if (jCon == 0) {
//							continue;
//						}
//						int el = eiUVWijk[lMat][1];
//						int iuvw = eiUVWijk[lMat][2];
//						int i = eiUVWijk[lMat][3];
//						int j = eiUVWijk[lMat][4];
//						int k = eiUVWijk[lMat][5];
//						cSolB = 0;
//// ------>fix			//findD
//						Integral_d d_integrator(&elements[el-1]);
//						d_integrator.findD(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolB, size1);
//						cSensitivity += cSolB*cAlphaForward[jCon] * std::conj(cAlphaAdjoint[iCon])*cEpsilon - std::complex<double>(1.0, 0.0);
//					}
//				}
//
//			}
//		}
//	}
//}