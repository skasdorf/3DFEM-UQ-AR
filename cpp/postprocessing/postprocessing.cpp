#include "postprocessing.h"

std::vector<std::complex<double>> postproc::findesc(Domain& dom, int freqNumber, std::vector<double> iR)
{
	std::vector<std::complex<double>> cEsc(4);
#pragma omp parallel for num_threads(8)
	for (auto e = dom.elements.begin(); e != dom.elements.end(); ++e) {
		Integral_g gInts(&(*e));
		int nu = e->expansion[0];
		int nv = e->expansion[1];
		int nw = e->expansion[2];

		int nglu = e->quadrature[0];
		int nglv = e->quadrature[1];
		int nglw = e->quadrature[2];
		std::vector<double> xglu, xglv, xglw;
		functions::gaussk(nglu, xglu, e->wglu);
		xglu[0] = -1.0;
		xglu[nglu + 1] = 1.0;
		functions::gaussk(nglv, xglv, e->wglv);
		xglv[0] = -1.0;
		xglv[nglv + 1] = 1.0;
		functions::gaussk(nglw, xglw, e->wglw);
		xglw[0] = -1.0;
		xglw[nglw + 1] = 1.0;
		/*e->jacobian = matrix3d<double>(nglu + 1, nglv + 1, nglw + 1);
		unitVectorsM::unitaryvectors(nglu, nglv, nglw, xglu, xglv, xglw,
			e->geom_order, e->rs, e->rMatrix, e->auMatrix, e->avMatrix, e->awMatrix,
			e->jacobian, e->nRs, 2);*/
		for (int kMat = e->unknownsStart; kMat <= e->unknownsEnd; ++kMat) {
			std::vector < std::complex<double> > cFF(4);
			int ihCon = abs(dom.vectorD[kMat]);
			if (ihCon == 0) continue;

			int eh = dom.eiUVWijk[kMat][1];
			int iuvwh = dom.eiUVWijk[kMat][2];
			int ih = dom.eiUVWijk[kMat][3];
			int jh = dom.eiUVWijk[kMat][4];
			int kh = dom.eiUVWijk[kMat][5];
			for (int face = 1; face <= 6; ++face) {
				if (((dom.elements[eh - 1]).abcList[face - 1] != 1) &&
					dom.elements[eh - 1].PML_boundary_list[face - 1] != 1) {
					continue;
				}
				if (((dom.elements[eh - 1]).PML_boundary_list[face - 1] == 1) && dom.elements[eh - 1].materials.region == 2) {
					continue;
				}
				gInts.findFF(iuvwh, ih, jh, kh, face, dom.scatter1.K0[freqNumber], iR, cFF);
				if (dom.vectorD[kMat] > 0) {
					#pragma omp critical
					cEsc = cEsc + cFF*dom.cAlpha[ihCon];
				}
				else {
					#pragma omp critical
					cEsc = cEsc +  cFF*(-dom.cAlpha[ihCon]);
				}
			}
		}
	}
	return cEsc;
}

double postproc::findrcs(Domain& dom, std::vector<std::complex<double>> cEsc, int waveNumber)
{
	double magEsc = abs(cEsc[1])*abs(cEsc[1]) + abs(cEsc[2])*abs(cEsc[2]) + abs(cEsc[3])*abs(cEsc[3]);
	auto cEteta = dom.scatter1.waveTheta[waveNumber];
	auto cEfi = dom.scatter1.wavePhi[waveNumber];
	double magEin = abs(cEteta)*abs(cEteta) + abs(cEfi)*abs(cEfi);
	double sigmaRCS = magEsc / magEin / 4.0 / constants::PI;
	return sigmaRCS;
}

std::vector<double> postproc::post(Domain& dom)
{
	std::vector<double> rcs_list;
	double tetaStep, fiStep;
	if (dom.scatter1.nTh != 1) {
		tetaStep = (dom.scatter1.thStop - dom.scatter1.thStart) / (dom.scatter1.nTh - 1);
	}
	else {
		tetaStep = 0.0;
	}
	if (dom.scatter1.nPh != 1) {
		fiStep = (dom.scatter1.phStop - dom.scatter1.phStart) / (dom.scatter1.nPh - 1);
	}
	else {
		fiStep = 0.0;
	}
	std::vector<double> tetas(dom.scatter1.nTh), fis(dom.scatter1.nPh);
	for (int i0 = 1; i0 <= dom.scatter1.nTh; ++i0) {
		tetas[i0 - 1] = dom.scatter1.thStart + (i0 - 1)*tetaStep;
	}
	for (int j0 = 1; j0 <= dom.scatter1.nPh; ++j0) {
		fis[j0 - 1] = dom.scatter1.phStart + (j0 - 1)*fiStep;
	}
	for (int i0 = 1; i0 <= dom.scatter1.nTh; ++i0) {
		for (int j0 = 1; j0 <= dom.scatter1.nPh; ++j0) {
			auto nSc = coordinates::findnAr(tetas[i0-1], fis[j0-1]);
			
			nSc = nSc*-1.0;
			//std::cout << "nSc: " << nSc[1] << " " << nSc[2] << " " << nSc[3] << std::endl;
			auto cEsc = findesc(dom, 1, nSc);
			//std::cout << "cEsc: " << cEsc[1] << " " << cEsc[2] << " " << cEsc[3] << std::endl;
			auto rcs = findrcs(dom, cEsc, 1);
		//	std::cout << "rcs: " << rcs << std::endl;
			rcs_list.push_back(rcs);
		}
	}
	return rcs_list;
}

std::vector<double> postproc::postplus(Domain& dom, std::vector<std::complex<double>>& errorterm) //RCS computation which includes adjoint error estimate
{
	std::vector<double> rcs_list;
	double tetaStep, fiStep;
	if (dom.scatter1.nTh != 1) {
		tetaStep = (dom.scatter1.thStop - dom.scatter1.thStart) / (dom.scatter1.nTh - 1);
	}
	else {
		tetaStep = 0.0;
	}
	if (dom.scatter1.nPh != 1) {
		fiStep = (dom.scatter1.phStop - dom.scatter1.phStart) / (dom.scatter1.nPh - 1);
	}
	else {
		fiStep = 0.0;
	}
	std::vector<double> tetas(dom.scatter1.nTh), fis(dom.scatter1.nPh);
	for (int i0 = 1; i0 <= dom.scatter1.nTh; ++i0) {
		tetas[i0 - 1] = dom.scatter1.thStart + (i0 - 1)*tetaStep;
	}
	for (int j0 = 1; j0 <= dom.scatter1.nPh; ++j0) {
		fis[j0 - 1] = dom.scatter1.phStart + (j0 - 1)*fiStep;
	}
	for (int i0 = 1; i0 <= dom.scatter1.nTh; ++i0) {
		for (int j0 = 1; j0 <= dom.scatter1.nPh; ++j0) {
			auto nSc = coordinates::findnAr(tetas[i0 - 1], fis[j0 - 1]);

			nSc = nSc * -1.0;
			//std::cout << "nSc: " << nSc[1] << " " << nSc[2] << " " << nSc[3] << std::endl;
			auto cEsc = findesc(dom, 1, nSc);
			cEsc = cEsc + errorterm;
			//std::cout << "cEsc: " << cEsc[1] << " " << cEsc[2] << " " << cEsc[3] << std::endl;
			auto rcs = findrcs(dom, cEsc, 1);
			//	std::cout << "rcs: " << rcs << std::endl;
			rcs_list.push_back(rcs);
		}
	}
	return rcs_list;
}