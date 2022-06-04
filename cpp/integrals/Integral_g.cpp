#define PI2 1.57079632679489661

#include "Integral_g.h"

#include <chrono>
inline void cdiv(const std::complex<double>& c1, const std::complex<double>& c2, std::complex<double>& c3) {
	double abs_val = (c2.real() * c2.real() + c2.imag() * c2.imag());
	//method 1
	//res += vals[i] / (abs_val)*(std::complex<double>(vals[i+1].real(), -vals[i+1].imag()));
	//method 2
	/*std::complex<double> temp = c1 / abs_val;
	c3.real(temp.real()*c2.real() + temp.imag()*c2.imag());
	c3.imag(temp.imag()*c2.real() - temp.real()*c2.imag());*/

	c3.real((c1.real() * c2.real() + c1.imag() * c2.imag()) / abs_val);
	c3.imag((c1.imag() * c2.real() - c1.real() * c2.imag()) / abs_val);
}
inline std::complex<double> cdiv(const std::complex<double>& c1, const std::complex<double>& c2) {
	double abs_val = (c2.real() * c2.real() + c2.imag() * c2.imag());
	//method 1
	//res += vals[i] / (abs_val)*(std::complex<double>(vals[i+1].real(), -vals[i+1].imag()));
	//method 2
	std::complex<double> temp = c1 / abs_val;
	return std::complex<double>((temp.real() * c2.real() + temp.imag() * c2.imag()),
		(temp.imag() * c2.real() - temp.real() * c2.imag()));
}
//void dotc(const std::vector<double>& v1, const std::vector < std::complex<double>>&v2, std::complex<double>& result) {
//	result.real(v1[1] * v2[1].real() + v1[2] * v2[2].real() + v1[3] * v2[3].real());
//	result.imag(v1[1] * v2[1].imag() + v1[2] * v2[2].imag() + v1[3] * v2[3].imag());
//}


void Integral_g::findFF(const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& face, const double& k0, const std::vector<double>& iR, std::vector<std::complex<double>>& cFF)
{
	int iuvwn, m, iuvw3, n;
	int sign1, sign2;
	double fuvw, rir;
	std::vector<double> asec(4), r(4), vec1(4), vec2(4), rfuvws(4);
	std::vector<std::complex<double>> cS(4);
	matrix2d<double> auvw(3, 3);

	switch (face) {
	case 1:
	case 2:
		iuvwn = 1;
		switch (face) {
		case 1:
			m = 0;
			if (iuvwh == 2) {
				sign1 = -1;
				iuvw3 = 3;
			}
			else {
				sign1 = 1;
				iuvw3 = 2;
			}
			sign2 = -1;
			break;
		case 2:
			m = e->quadrature[0] + 1;
			if (iuvwh == 2) {
				sign1 = 1;
				iuvw3 = 3;
			}
			else {
				sign1 = -1;
				iuvw3 = 2;
			}
			sign2 = 1;
			break;
		}
		for (n = 1; n <= e->quadrature[1]; ++n) {
			cS[1] = 0; cS[2] = 0; cS[3] = 0;
			for (int l = 1; l <= e->quadrature[2]; ++l) {
				auvw[1] = e->auMatrix(m, n, l);
				auvw[2] = e->avMatrix(m, n, l);
				auvw[3] = e->awMatrix(m, n, l);
				products::dot(e->rMatrix(m, n, l), iR, rir);
				functions::FINDSECONDARYUNIT(iuvwn, auvw[1], auvw[2], auvw[3], asec);
				ROTFSEC(iuvwh, m, n, l, ih, jh, kh, auvw[1], auvw[2], auvw[3], rfuvws);
				products::cross(asec, rfuvws, vec2);
				double jacobian = e->jacobian(m, n, l);
				if (iuvwh != iuvwn) {
					products::cross(auvw(iuvw3), iR, vec1);
					functions::FOFUVW(iuvwh, m, n, l, ih, jh, kh, e, fuvw);

					cS[1].real(cS[1].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[1] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[1] * sign1));
					cS[2].real(cS[2].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[2] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[2] * sign1));
					cS[3].real(cS[3].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[3] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[3] * sign1));

					cS[1].imag(cS[1].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[1] * sign1) + std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[2] * sign1) + std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[3] * sign1) + std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
				else {
					cS[1].real(cS[1].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].real(cS[2].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].real(cS[3].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[3] / jacobian * sign2));

					cS[1].imag(cS[1].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
			}
			cFF[1] += e->wglv[n] * cS[1];
			cFF[2] += e->wglv[n] * cS[2];
			cFF[3] += e->wglv[n] * cS[3];

		}
		break;
	case 3:
	case 4:
		iuvwn = 2;
		switch (face) {
		case 3:
			n = 0;
			if (iuvwh == 1) {
				sign1 = 1;
				iuvw3 = 3;
			}
			else {
				sign1 = -1;
				iuvw3 = 1;
			}
			sign2 = -1;
			break;
		case 4:
			n = e->quadrature[1] + 1;
			if (iuvwh == 1) {
				sign1 = -1;
				iuvw3 = 3;

			}
			else {
				sign1 = 1;
				iuvw3 = 1;
			}
			sign2 = 1;
			break;
		}
		for (m = 1; m <= e->quadrature[0]; ++m) {
			cS[1] = 0; cS[2] = 0; cS[3];
			for (int l = 1; l <= e->quadrature[2]; ++l) {
				auvw[1] = e->auMatrix(m, n, l);
				auvw[2] = e->avMatrix(m, n, l);
				auvw[3] = e->awMatrix(m, n, l);
				double jacobian = e->jacobian(m, n, l);
				products::dot(e->rMatrix(m, n, l), iR, rir);
				functions::FINDSECONDARYUNIT(iuvwn, auvw[1], auvw[2], auvw[3], asec);
				ROTFSEC(iuvwh, m, n, l, ih, jh, kh, auvw[1], auvw[2], auvw[3], rfuvws);
				products::cross(asec, rfuvws, vec2);
				if (iuvwh != iuvwn) {
					products::cross(auvw(iuvw3), iR, vec1);
					functions::FOFUVW(iuvwh, m, n, l, ih, jh, kh, e, fuvw);
					cS[1].real(cS[1].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[1] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[1] * sign1));
					cS[2].real(cS[2].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[2] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[2] * sign1));
					cS[3].real(cS[3].real() + e->wglw[l] * (std::cos(k0 * rir) * (vec2[3] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[3] * sign1));

					cS[1].imag(cS[1].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[1] * sign1) + std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[2] * sign1) + std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglw[l] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[3] * sign1) + std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
				else {
					cS[1].real(cS[1].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].real(cS[2].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].real(cS[3].real() + e->wglw[l] * (std::cos(k0 * rir) * vec2[3] / jacobian * sign2));

					cS[1].imag(cS[1].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglw[l] * (std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
			}
			cFF[1] += e->wglu[m] * cS[1];
			cFF[2] += e->wglu[m] * cS[2];
			cFF[3] += e->wglu[m] * cS[3];
		}
		break;
	case 5:
	case 6:
		int l;
		iuvwn = 3;
		switch (face) {
		case 5:
			l = 0;
			if (iuvwh == 1) {
				sign1 = -1;
				iuvw3 = 2;

			}
			else {
				sign1 = 1;
				iuvw3 = 1;
			}
			sign2 = -1;
			break;
		case 6:
			l = e->quadrature[2] + 1;
			if (iuvwh == 1) {
				sign1 = 1;
				iuvw3 = 2;
			}
			else {
				sign1 = -1;
				iuvw3 = 1;
			}
			sign2 = 1;
		}
		for (m = 1; m <= e->quadrature[0]; ++m) {
			cS[1] = 0; cS[2] = 0; cS[3] = 0;
			for (n = 1; n <= e->quadrature[1]; ++n) {
				auvw[1] = e->auMatrix(m, n, l);
				auvw[2] = e->avMatrix(m, n, l);
				auvw[3] = e->awMatrix(m, n, l);
				double jacobian = e->jacobian(m, n, l);
				products::dot(e->rMatrix(m, n, l), iR, rir);
				functions::FINDSECONDARYUNIT(iuvwn, auvw[1], auvw[2], auvw[3], asec);
				ROTFSEC(iuvwh, m, n, l, ih, jh, kh, auvw[1], auvw[2], auvw[3], rfuvws);
				products::cross(asec, rfuvws, vec2);
				if (iuvwh != iuvwn) {
					products::cross(auvw(iuvw3), iR, vec1);
					functions::FOFUVW(iuvwh, m, n, l, ih, jh, kh, e, fuvw);
					cS[1].real(cS[1].real() + e->wglv[n] * (std::cos(k0 * rir) * (vec2[1] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[1] * sign1));
					cS[2].real(cS[2].real() + e->wglv[n] * (std::cos(k0 * rir) * (vec2[2] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[2] * sign1));
					cS[3].real(cS[3].real() + e->wglv[n] * (std::cos(k0 * rir) * (vec2[3] / jacobian * sign2) - std::sin(k0 * rir) * k0 * fuvw * vec1[3] * sign1));

					cS[1].imag(cS[1].imag() + e->wglv[n] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[1] * sign1) + std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglv[n] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[2] * sign1) + std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglv[n] * (std::cos(k0 * rir) * (k0 * fuvw * vec1[3] * sign1) + std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
				else {
					cS[1].real(cS[1].real() + e->wglv[n] * (std::cos(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].real(cS[2].real() + e->wglv[n] * (std::cos(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].real(cS[3].real() + e->wglv[n] * (std::cos(k0 * rir) * vec2[3] / jacobian * sign2));

					cS[1].imag(cS[1].imag() + e->wglv[n] * (std::sin(k0 * rir) * vec2[1] / jacobian * sign2));
					cS[2].imag(cS[2].imag() + e->wglv[n] * (std::sin(k0 * rir) * vec2[2] / jacobian * sign2));
					cS[3].imag(cS[3].imag() + e->wglv[n] * (std::sin(k0 * rir) * vec2[3] / jacobian * sign2));
				}
			}
			cFF[1] += e->wglu[m] * cS[1];
			cFF[2] += e->wglu[m] * cS[2];
			cFF[3] += e->wglu[m] * cS[3];
		}
	}
}


void Integral_g::findGWave(Scatter* scatter1, const std::vector<int>& vectorD, const int& useAdjoint, const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
	const int& kMat, const int& iCon, const int& face) {


	cFF.resize(4);
	std::complex<double> cSolGtotal = 0.0, cSolG_eps = 0.0;
	std::vector<double> iR(4);
	std::complex<double> cj(0, 1);
	for (int waveNumber = 1; waveNumber <= scatter1->numberOfWaves; waveNumber++) {
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

		for (int freqNumber = 1; freqNumber <= scatter1->nF; freqNumber++) {
			//auto val1 = std::exp(-cj*scatter1->K0[freqNumber]);
			double val1 = -scatter1->K0[freqNumber];
			cFF[1] = 0; cFF[2] = 0; cFF[3] = 0;
			if (useAdjoint) {
				findFF(iuvwh, ih, jh, kh, face, scatter1->K0[freqNumber], iR, cFF);
				// cFF[3] is extracting the third component of the output scattered field vector
				// Ideally, there should be an input vector "w" such that we take cFF.dot(w) to generate the scalar instead
				// For selecting an arbitrary direction, but this is good for now.
				if (vectorD[kMat] > 0) {
					scatter1->cGr[freqNumber][iCon] = scatter1->cGr[freqNumber][iCon] + std::conj(this->cFF[3]);
				}
				else {
					scatter1->cGr[freqNumber][iCon] = scatter1->cGr[freqNumber][iCon] - std::conj(this->cFF[3]);
				}
			}
			else {

				findGInc(scatter1, vectorD, iuvwh, ih, jh, kh, size1, waveNumber, freqNumber, cSolGtotal, asec2, hVector, asec3, val1, cSolG_eps);

				if (cSolGtotal != 0.0 || cSolG_eps != 0.0) {

					if (vectorD[kMat] < 0)
						cSolG_eps *= -1.0;

					scatter1->cGr[freqNumber][iCon] = scatter1->cGr[freqNumber][iCon] - scatter1->K0[freqNumber] * scatter1->K0[freqNumber] * cSolGtotal;
					scatter1->cGr_eps[iCon] = scatter1->cGr_eps[iCon] - scatter1->K0[freqNumber] * scatter1->K0[freqNumber] * cSolG_eps;

					////////////////////////////////Frequency Perturbation//////////////////////////////////////////
					//e->cGr_eps_el[iCon] = -1.0*scatter1->K0[freqNumber] * 2.0 *cSolG_eps;
					//---------------------------------------------------------------------------------------------

					////////////////////////////////Material Perturbation//////////////////////////////////////////
					// This is suitable for the case when the material parameter is uniform in a cell
					e->cGr_eps_el[iCon] = -1.0 * scatter1->K0[freqNumber] * scatter1->K0[freqNumber] * cSolG_eps / (e->materials.epsr_list[0][0] - 1.0);
					//---------------------------------------------------------------------------------------------

				}

			}
		}
	}
}


void Integral_g::findGInc(Scatter* scatter1, const std::vector<int>& vectorD, const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
	const int& waveNumber, const int& freqNumber, std::complex<double>& cIntDInc, std::vector<std::complex<double>>& asec2, std::vector<std::complex<double>>& hVector, std::vector<std::complex<double>>& asec3, double val1, std::complex<double>& cIntDInc_eps) {

	std::vector<double> asec1(4), rfus1(4);
	std::vector<double> a1(4), a2(4);
	double f1;
	std::vector<std::complex<double>> f3(4), f4(4), asec2_t(4), asec3_t(4);

	std::complex<double> f2, f;
	std::complex<double> cEpsilon, cMu;

	std::complex<double> cj(0, 1);
	int m, n, l;
	std::complex<double> cSn, cSl, cSl_eps, cSn_eps;
	int coord;
	cIntDInc = 0.0;
	cIntDInc_eps = 0.0;
	for (m = 1; m <= e->wglu.size() - 1; m++) {
		cSn = 0.0;
		cSn_eps = 0.0;
		for (n = 1; n <= e->wglv.size() - 1; n++) {
			double mno1_part;
			double mno2_part;
			cSl = 0;
			cSl_eps = 0;
			double fuvw_part;
			if (iuvwh == 1) {
				fuvw_part = e->uPowers(m, ih) * e->fvPowers(n, jh);
				mno1_part = e->uPowers(m, ih) * e->fvPowers(n, jh);
				mno2_part = e->uPowers(m, ih) * e->fpvPowers(n, jh);
			}
			else if (iuvwh == 2) {
				fuvw_part = e->fuPowers(m, ih) * e->vPowers(n, jh);
				mno1_part = e->fpuPowers(m, ih) * e->vPowers(n, jh);
				mno2_part = e->fuPowers(m, ih) * e->vPowers(n, jh);
			}
			else {
				fuvw_part = e->fuPowers(m, ih) * e->fvPowers(n, jh);// *e->wPowers(l, kh);
				mno1_part = e->fuPowers(m, ih) * e->fpvPowers(n, jh);
				mno2_part = e->fpuPowers(m, ih) * e->fvPowers(n, jh);
			}
			for (l = 1; l <= e->wglw.size() - 1; l++) {
				double mno1, mno2;
				if (iuvwh == 1) {
					f1 = fuvw_part * e->fwPowers(l, kh);
					mno1 = mno1_part * e->fpwPowers(l, kh);
					mno2 = mno2_part * e->fwPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {

						a1[coord] = e->avMatrix(m, n, l, coord);
						a2[coord] = e->awMatrix(m, n, l, coord);
					}
				}
				else if (iuvwh == 2) {
					f1 = fuvw_part * e->fwPowers(l, kh);
					mno1 = mno1_part * e->fwPowers(l, kh);
					mno2 = mno2_part * e->fpwPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {
						a2[coord] = e->auMatrix(m, n, l, coord);
						a1[coord] = e->awMatrix(m, n, l, coord);
					}
				}
				else {
					f1 = fuvw_part * e->wPowers(l, kh);
					mno1 = mno1_part * e->wPowers(l, kh);
					mno2 = mno2_part * e->wPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {
						a1[coord] = e->auMatrix(m, n, l, coord);
						a2[coord] = e->avMatrix(m, n, l, coord);

					}
				}
				functions::FINDSECONDARYUNIT(iuvwh, a1, a2, asec1);
				double rn;
				products::dot(e->rMatrix(m, n, l), scatter1->nAr[waveNumber], rn);
				std::complex<double> val1_t;
				val1_t.imag(sin(val1 * rn));
				val1_t.real(sin(val1 * rn + PI2));

				std::complex<double> val2;
				cdiv(val1_t, scatter1->K0[freqNumber], val2);
				asec2_t[1] = asec2[1] * val1_t;
				asec2_t[2] = asec2[2] * val1_t;
				asec2_t[3] = asec2[3] * val1_t;

				if (e->materials.icode == 0) {
					if (e->materials.hcode == 0) {

						MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec2, m, n, l, f3);

						f3[1] = f3[1] - asec2_t[1];
						f3[2] = f3[2] - asec2_t[2];
						f3[3] = f3[3] - asec2_t[3];

					}
					else {

						MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec2_t, 1, 1, 1, f3);

						f3[1] = f3[1] - asec2_t[1];
						f3[2] = f3[2] - asec2_t[2];
						f3[3] = f3[3] - asec2_t[3];

					}
				}
				else {
					if (e->materials.hcode == 0) {

						f3[1] = asec2_t[1] * (e->EpsRelInt(1, m, n, l) - 1.0);
						f3[2] = asec2_t[2] * (e->EpsRelInt(1, m, n, l) - 1.0);
						f3[3] = asec2_t[3] * (e->EpsRelInt(1, m, n, l) - 1.0);

					}
					else {
						f3[1] = asec2_t[1] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);
						f3[2] = asec2_t[2] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);
						f3[3] = asec2_t[3] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);

					}
				}
				f = f3[1] * asec1[1] + f3[2] * asec1[2] + f3[3] * asec1[3];

				ROTFSEC(iuvwh, mno1, mno2, a1, a2, rfus1);
				asec3_t[1] = asec3[1] * val2;
				asec3_t[2] = asec3[2] * val2;
				asec3_t[3] = asec3[3] * val2;


				if (e->materials.icode == 0) {
					if (e->materials.hcode == 0) {

						MULTIPLYMURMATRIX(size1, e->MuRelIntInv, asec3, m, n, l, f4);

						f4[1] = f4[1] - hVector[1] * val2;
						f4[2] = f4[2] - hVector[2] * val2;
						f4[3] = f4[3] - hVector[3] * val2;


					}
					else {

						MULTIPLYMURMATRIX(size1, e->MuRelIntInv, asec3, 1, 1, 1, f4);

						f4[1] -= asec3_t[1];
						f4[2] -= asec3_t[2];
						f4[3] -= asec3_t[3];

					}
				}
				else {
					if (e->materials.hcode == 0) {
						cdiv(asec3_t[1], e->MuRelInt(1, m, n, l), f4[1]);
						cdiv(asec3_t[2], e->MuRelInt(1, m, n, l), f4[2]);
						cdiv(asec3_t[3], e->MuRelInt(1, m, n, l), f4[3]);
						f4[1] -= asec3_t[1];
						f4[2] -= asec3_t[2];
						f4[3] -= asec3_t[3];

					}
					else {


						std::complex<double> mui = e->MuRelInt(1, 1, 1, 1);
						cdiv(asec3_t[1], mui, f4[1]);
						cdiv(asec3_t[2], mui, f4[2]);
						cdiv(asec3_t[3], mui, f4[3]);
						f4[1] -= asec3_t[1];
						f4[2] -= asec3_t[2];
						f4[3] -= asec3_t[3];
					}
				}//else (aCode!=0)

				products::dotc(rfus1, f4, f2);

				cSl.real(cSl.real() + e->wglw[l] * (f2.real() + f1 * f.real()));
				cSl.imag(cSl.imag() + e->wglw[l] * (f2.imag() + f1 * f.imag()));
				cSl_eps.real(cSl_eps.real() + e->wglw[l] * (f1 * f.real()));
				cSl_eps.imag(cSl_eps.imag() + e->wglw[l] * (f1 * f.imag()));

			}//for l
			cSn = cSn + e->wglv[n] * cSl;
			cSn_eps = cSn_eps + e->wglv[n] * cSl_eps;

		}//for n
		cIntDInc += +e->wglu[m] * cSn;
		cIntDInc_eps += +e->wglu[m] * cSn_eps;
	}

}//findGInc



void Integral_g::MULTIPLYMURMATRIX(const int& size1, const matrix4d<std::complex<double>>& MuRelInv,
	const std::vector<std::complex<double>>& vector, const int& m, const int& n, const int& l, std::vector<std::complex<double>>& out) {
	//fix to explicitly set the real and imag parts, might speed things up
	switch (size1) {
	case (9):
		out[1] = vector[1] * MuRelInv(1, m, n, l) + vector[2] * MuRelInv(2, m, n, l) + vector[3] * MuRelInv(3, m, n, l);
		out[2] = vector[1] * MuRelInv(4, m, n, l) + vector[2] * MuRelInv(5, m, n, l) + vector[3] * MuRelInv(6, m, n, l);
		out[3] = vector[1] * MuRelInv(7, m, n, l) + vector[2] * MuRelInv(8, m, n, l) + vector[3] * MuRelInv(9, m, n, l);
		/*out[1].real(vector[1] * MuRelInv(1, m, n, l).real() + vector[2] * MuRelInv(2, m, n, l).real() + vector[3] * MuRelInv(3, m, n, l).real());
		out[1].imag(vector[1] * MuRelInv(1, m, n, l).imag() + vector[2] * MuRelInv(2, m, n, l).imag() + vector[3] * MuRelInv(3, m, n, l).imag());
		out[2].real(vector[1] * MuRelInv(4, m, n, l).real() + vector[2] * MuRelInv(5, m, n, l).real() + vector[3] * MuRelInv(6, m, n, l).real());
		out[2].imag(vector[1] * MuRelInv(4, m, n, l).imag() + vector[2] * MuRelInv(5, m, n, l).imag() + vector[3] * MuRelInv(6, m, n, l).imag());
		out[3].real(vector[1]*MuRelInv(7, m, n, l).real() + vector[2] * MuRelInv(8, m, n, l).real() + vector[3] * MuRelInv(9, m, n, l).real());
		out[3].imag(vector[1] * MuRelInv(7, m, n, l).imag() + vector[2] * MuRelInv(8, m, n, l).imag() + vector[3] * MuRelInv(9, m, n, l).imag());*/
		break;
	case 6:
		out[1] = vector[1] * MuRelInv(1, m, n, l) + vector[2] * MuRelInv(2, m, n, l) + vector[3] * MuRelInv(3, m, n, l);
		out[2] = vector[1] * MuRelInv(2, m, n, l) + vector[2] * MuRelInv(4, m, n, l) + vector[3] * MuRelInv(5, m, n, l);
		out[3] = vector[1] * MuRelInv(3, m, n, l) + vector[2] * MuRelInv(5, m, n, l) + vector[3] * MuRelInv(6, m, n, l);
		break;
	case 3:
		out[1] = vector[1] * MuRelInv(1, m, n, l);
		out[2] = vector[2] * MuRelInv(2, m, n, l);
		out[3] = vector[3] * MuRelInv(3, m, n, l);
		break;
	default:

		break;
	}
}
//fix this method to not be BS
void Integral_g::ROTFSEC(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw,
	std::vector<double>& rfs) {
	switch (iuvw) {
	case 1:
		ROTFUSEC(m, n, l, i, j, k, av, aw, rfs);
		break;
	case 2:
		ROTFVSEC(m, n, l, i, j, k, aw, au, rfs);
		break;
	case 3:
		ROTFWSEC(m, n, l, i, j, k, au, av, rfs);
		break;
	default:
		std::cout << "Invalid basis direction in ROTSEC!\n";
		break;
	}
}
//fix this method to not be BS
void Integral_g::ROTFUSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& av, const std::vector<double>& aw,
	std::vector<double>& rfus) {

	double mno1 = e->uPowers(m, i) * e->fvPowers(n, j) * e->fpwPowers(l, k);
	double mno2 = e->uPowers(m, i) * e->fpvPowers(n, j) * e->fwPowers(l, k);

	/*for (int cord = 1; cord <= 3; ++cord) {
		rfus[cord] = mno1 * av[cord] - mno2 * aw[cord];
	}*/
	rfus[1] = mno1 * av[1] - mno2 * aw[1];
	rfus[2] = mno1 * av[2] - mno2 * aw[2];
	rfus[3] = mno1 * av[3] - mno2 * aw[3];

}
//fix this method to not be BS
void Integral_g::ROTFVSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& aw, const std::vector<double>& au,
	std::vector<double>& rfvs) {
	double mno1 = e->fpuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	double mno2 = e->fuPowers(m, i) * e->vPowers(n, j) * e->fpwPowers(l, k);

	//for (int cord = 1; cord <= 3; ++cord) {
	//	rfvs[cord] = mno1 * aw[cord] - mno2 * au[cord];
	//}
	rfvs[1] = mno1 * aw[1] - mno2 * au[1];
	rfvs[2] = mno1 * aw[2] - mno2 * au[2];
	rfvs[3] = mno1 * aw[3] - mno2 * au[3];
}
//fix this method ot not be bs
void Integral_g::ROTFWSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
	std::vector<double>& rfws) {
	double mno1 = e->fuPowers(m, i) * e->fpvPowers(n, j) * e->wPowers(l, k);
	double mno2 = e->fpuPowers(m, i) * e->fvPowers(n, j) * e->wPowers(l, k);

	/*for (int cord = 1; cord <= 3; ++cord) {
		rfws[cord] = mno1 * au[cord] - mno2 * av[cord];
	}*/
	rfws[1] = mno1 * au[1] - mno2 * av[1];
	rfws[2] = mno1 * au[2] - mno2 * av[2];
	rfws[3] = mno1 * au[3] - mno2 * av[3];
}


//fix this method to not be BS
void Integral_g::ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& a1, const std::vector<double>& a2,
	std::vector<double>& rfs) {
	switch (iuvw) {
	case 1:
		ROTFUSEC(mno1, mno2, a1, a2, rfs);
		break;
	case 2:
		ROTFVSEC(mno1, mno2, a1, a2, rfs);
		break;
	case 3:
		ROTFWSEC(mno1, mno2, a1, a2, rfs);
		break;
	default:
		std::cout << "Invalid basis direction in ROTSEC!\n";
		break;
	}
}
//fix this method to not be BS
void Integral_g::ROTFUSEC(const double& mno1, const double& mno2, const std::vector<double>& av, const std::vector<double>& aw,
	std::vector<double>& rfus) {

	//double mno1 = e->uPowers(m, i) * e->fvPowers(n, j)  * e->fpwPowers(l, k);
	//double mno2 = e->uPowers(m, i) * e->fpvPowers(n, j) * e->fwPowers(l, k);

	/*for (int cord = 1; cord <= 3; ++cord) {
	rfus[cord] = mno1 * av[cord] - mno2 * aw[cord];
	}*/
	rfus[1] = mno1 * av[1] - mno2 * aw[1];
	rfus[2] = mno1 * av[2] - mno2 * aw[2];
	rfus[3] = mno1 * av[3] - mno2 * aw[3];

}
//fix this method to not be BS
void Integral_g::ROTFVSEC(const double& mno1, const double& mno2, const std::vector<double>& aw, const std::vector<double>& au,
	std::vector<double>& rfvs) {
	//double mno1 = e->fpuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	//double mno2 = e->fuPowers(m, i)  * e->vPowers(n, j) * e->fpwPowers(l, k);

	//for (int cord = 1; cord <= 3; ++cord) {
	//	rfvs[cord] = mno1 * aw[cord] - mno2 * au[cord];
	//}
	rfvs[1] = mno1 * aw[1] - mno2 * au[1];
	rfvs[2] = mno1 * aw[2] - mno2 * au[2];
	rfvs[3] = mno1 * aw[3] - mno2 * au[3];
}
//fix this method ot not be bs
void Integral_g::ROTFWSEC(const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
	std::vector<double>& rfws) {
	//double mno1 = e->fuPowers(m, i)  * e->fpvPowers(n, j) * e->wPowers(l, k);
	//double mno2 = e->fpuPowers(m, i) * e->fvPowers(n, j)  * e->wPowers(l, k);

	/*for (int cord = 1; cord <= 3; ++cord) {
	rfws[cord] = mno1 * au[cord] - mno2 * av[cord];
	}*/
	rfws[1] = mno1 * au[1] - mno2 * av[1];
	rfws[2] = mno1 * au[2] - mno2 * av[2];
	rfws[3] = mno1 * au[3] - mno2 * av[3];
}
void Integral_g::findGInc_EpsOnly(Scatter* scatter1, const std::vector<int>& vectorD, const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
	const int& waveNumber, const int& freqNumber, std::complex<double>& cIntDInc, std::vector<std::complex<double>>& asec2, std::vector<std::complex<double>>& hVector, std::vector<std::complex<double>>& asec3, double val1) {

	std::vector<double> asec1(4), rfus1(4);
	std::vector<double> a1(4), a2(4);
	double f1;
	std::vector<std::complex<double>> f3(4), f4(4), asec2_t(4), asec3_t(4);

	std::complex<double> f2, f;
	std::complex<double> cEpsilon, cMu;
	//std::vector<std::complex<double>> eVector(4), hVector(4), asec3(4);
	std::complex<double> cj(0, 1);
	int m, n, l;
	std::complex<double> cSn, cSl;
	int coord;
	cIntDInc = 0;
	for (m = 1; m <= e->wglu.size() - 1; m++) {
		cSn = 0;
		for (n = 1; n <= e->wglv.size() - 1; n++) {
			double mno1_part;
			double mno2_part;


			cSl = 0;
			double fuvw_part;
			if (iuvwh == 1) {
				fuvw_part = e->uPowers(m, ih) * e->fvPowers(n, jh);
				mno1_part = e->uPowers(m, ih) * e->fvPowers(n, jh);
				mno2_part = e->uPowers(m, ih) * e->fpvPowers(n, jh);
			}
			else if (iuvwh == 2) {
				fuvw_part = e->fuPowers(m, ih) * e->vPowers(n, jh);
				mno1_part = e->fpuPowers(m, ih) * e->vPowers(n, jh);
				mno2_part = e->fuPowers(m, ih) * e->vPowers(n, jh);
			}
			else {
				fuvw_part = e->fuPowers(m, ih) * e->fvPowers(n, jh);// *e->wPowers(l, kh);
				mno1_part = e->fuPowers(m, ih) * e->fpvPowers(n, jh);
				mno2_part = e->fpuPowers(m, ih) * e->fvPowers(n, jh);
			}
			for (l = 1; l <= e->wglw.size() - 1; l++) {
				double mno1, mno2;
				if (iuvwh == 1) {
					f1 = fuvw_part * e->fwPowers(l, kh);
					mno1 = mno1_part * e->fpwPowers(l, kh);
					mno2 = mno2_part * e->fwPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {

						a1[coord] = e->avMatrix(m, n, l, coord);
						a2[coord] = e->awMatrix(m, n, l, coord);
					}
				}
				else if (iuvwh == 2) {
					f1 = fuvw_part * e->fwPowers(l, kh);
					mno1 = mno1_part * e->fwPowers(l, kh);
					mno2 = mno2_part * e->fpwPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {
						a2[coord] = e->auMatrix(m, n, l, coord);
						a1[coord] = e->awMatrix(m, n, l, coord);
					}
				}
				else {
					f1 = fuvw_part * e->wPowers(l, kh);
					mno1 = mno1_part * e->wPowers(l, kh);
					mno2 = mno2_part * e->wPowers(l, kh);
					for (coord = 1; coord <= 3; coord++) {
						a1[coord] = e->auMatrix(m, n, l, coord);
						a2[coord] = e->avMatrix(m, n, l, coord);

					}
				}
				functions::FINDSECONDARYUNIT(iuvwh, a1, a2, asec1);
				double rn;
				products::dot(e->rMatrix(m, n, l), scatter1->nAr[waveNumber], rn);
				std::complex<double> val1_t;
				val1_t.imag(sin(val1 * rn));
				val1_t.real(sin(val1 * rn + PI2));

				//	auto val2 = val1 / scatter1->K0[freqNumber];
				std::complex<double> val2;
				cdiv(val1_t, scatter1->K0[freqNumber], val2);
				asec2_t[1] = asec2[1] * val1_t;
				asec2_t[2] = asec2[2] * val1_t;
				asec2_t[3] = asec2[3] * val1_t;

				if (e->materials.icode == 0) {
					if (e->materials.hcode == 0) {

						MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec2, m, n, l, f3);

						f3[1] = f3[1] - asec2_t[1];
						f3[2] = f3[2] - asec2_t[2];
						f3[3] = f3[3] - asec2_t[3];

					}
					else {

						MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec2_t, 1, 1, 1, f3);

						f3[1] = f3[1] - asec2_t[1];
						f3[2] = f3[2] - asec2_t[2];
						f3[3] = f3[3] - asec2_t[3];

					}
				}
				else {
					if (e->materials.hcode == 0) {

						f3[1] = asec2_t[1] * (e->EpsRelInt(1, m, n, l) - 1.0);
						f3[2] = asec2_t[2] * (e->EpsRelInt(1, m, n, l) - 1.0);
						f3[3] = asec2_t[3] * (e->EpsRelInt(1, m, n, l) - 1.0);

					}
					else {
						f3[1] = asec2_t[1] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);
						f3[2] = asec2_t[2] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);
						f3[3] = asec2_t[3] * (e->EpsRelInt(1, 1, 1, 1) - 1.0);

					}
				}
				f = f3[1] * asec1[1] + f3[2] * asec1[2] + f3[3] * asec1[3];

				//ROTFSEC(iuvwh, m, n, l, ih, jh, kh, au, av, aw, rfus1);
				/*ROTFSEC(iuvwh, mno1, mno2, a1, a2, rfus1);
				asec3_t[1] = asec3[1] * val2;
				asec3_t[2] = asec3[2] * val2;
				asec3_t[3] = asec3[3] * val2;*/

				cSl.real(cSl.real() + e->wglw[l] * (f1 * f.real()));
				cSl.imag(cSl.imag() + e->wglw[l] * (f1 * f.imag()));

			}//for l
			cSn = cSn + e->wglv[n] * cSl;

		}//for n
		cIntDInc = cIntDInc + e->wglu[m] * cSn;
	}

}//findGInc


