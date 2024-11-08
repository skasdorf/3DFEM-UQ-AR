
#include "Integral_d.h"
#include "../functions/integralFunctions.h"
#include "../structure/Element.h"
#include <chrono>

void Integral_d::findD(const int& iuvwh, const int& iuvw, const int& ih, const int& jh, const int& kh, const int& i, const int& j, const int& k, std::complex<double>& cIntD,
	const int& size1) {


	std::complex<double> cSn, cSl, f;
	std::vector<double> au(4), av(4), aw(4), asec1(4), asec2(4);

	std::vector<std::complex<double>> f3(4);

	double f1, f2;
	//!iuvwh, iuvw - Testing and basis function directions
	//!ih, jh, kh - indices i^, j^, k^
	//!i, j, k - indices i, j, k
	//!nglu, nglv, nglw - Orders of GL integrations
	//!xglu, xglv, xglw - Coordinates of GL integration
	//!wglu, wglv, wglw - Weighting Coefficients for GL integration
	//!cIntD - d i^ j^ k^ i j k




	int nglu = e->wglu.size() - 1;
	int nglv = e->wglv.size() - 1;
	int nglw = e->wglw.size() - 1;
	int iHomCode = e->materials.hcode;

	for (int m = 1; m <= nglu; ++m) {
		cSn = 0;
		for (int n = 1; n <= nglv; ++n) {
			cSl = 0;
			double f1_part, f2_part;
			if (iuvwh == 1) {
				f1_part = e->uPowers(m, ih) * e->fvPowers(n, jh);

			}
			else if (iuvwh == 2) {
				f1_part = e->fuPowers(m, ih) * e->vPowers(n, jh);
			}
			else {
				f1_part = e->fuPowers(m, ih) * e->fvPowers(n, jh);
			}
			if (iuvw == 1) {
				f2_part = e->uPowers(m, i) * e->fvPowers(n, j);

			}
			else if (iuvw == 2) {
				f2_part = e->fuPowers(m, i) * e->vPowers(n, j);
			}
			else {
				f2_part = e->fuPowers(m, i) * e->fvPowers(n, j);
			}
			for (int l = 1; l <= nglw; ++l) {
				for (int coord = 1; coord <= 3; ++coord) {
					au[coord] = e->auMatrix(m, n, l, coord);
					av[coord] = e->avMatrix(m, n, l, coord);
					aw[coord] = e->awMatrix(m, n, l, coord);
				}
				if (iuvwh == 1) {
					f1 = f1_part * e->fwPowers(l, kh);

				}
				else if (iuvwh == 2) {
					f1 = f1_part * e->fwPowers(l, kh);
				}
				else {
					f1 = f1_part * e->wPowers(l, kh);
				}
				if (iuvw == 1) {
					f2 = f2_part * e->fwPowers(l, k);

				}
				else if (iuvw == 2) {
					f2 = f2_part * e->fwPowers(l, k);
				}
				else {
					f2 = f2_part * e->wPowers(l, k);
				}

				functions::FINDSECONDARYUNIT(iuvw, au, av, aw, asec1);
				functions::FINDSECONDARYUNIT(iuvwh, au, av, aw, asec2);

				if (e->materials.icode == 0) {
					if (iHomCode == 0) {

						functions::MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec1, m, n, l, f3);


					}
					else {


						functions::MULTIPLYMURMATRIX(size1, e->EpsRelInt, asec1, m, n, l, f3);

					}
				}
				else {
					if (iHomCode == 0) {
						//	auto t16 = std::chrono::high_resolution_clock::now();
						std::vector<std::complex<double>> temp;
						for (auto val : asec1) {
							temp.push_back(val * e->EpsRelInt(1, m, n, l));
						}
						f3 = temp;

					}
					else {
						//	auto t16 = std::chrono::high_resolution_clock::now();
						f3[1] = asec1[1];
						f3[2] = asec1[2];
						f3[3] = asec1[3];

					}
				}

				products::dotc(asec2, f3, f);

				double temp = f1 * f2 * e->wglw[l] / e->jacobian(m, n, l);
				cSl.real(cSl.real() + temp * (f.real()));
				cSl.imag(cSl.imag() + temp * (f.imag()));

			}
			cSn = cSn + e->wglv[n] * cSl;
		}
		cIntD = cIntD + e->wglu[m] * cSn;
	}
}
