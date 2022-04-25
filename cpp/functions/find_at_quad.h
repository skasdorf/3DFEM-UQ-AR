#pragma once
void make_evals(int nglu, int nglv, int nglw) {
	for (auto e = elements.begin(); e != elements.end(); ++e) {
		
		std::vector<std::complex<double>> basisFunc(nglu*nglv*nglw), basisDeriv(nglu*nglv*nglw);
			for (int iUnkns = e->unknownsStart; iUnkns <= e->unknownsEnd; ++iUnkns) {
				int iuvw = eiUVWijk[iUnkns][2];
				int i = eiUVWijk[iUnkns][3];
				int j = eiUVWijk[iUnkns][4];
				int k = eiUVWijk[iUnkns][5];
				//iterate through quadrature points
				for (int m = 1; m <= nglu; ++m) {
					for (int n = 1; n <= nglv; ++n) {
						for (int l = 1; l <= nglw; ++l) {
							std::complex<double> f1;
							intcROTFSEC(iuvw, m, n, l, i, j, k, au, av, aw, rfs);
							rfs = cAlpha[iUnkns] * rfs;
							basisDeriv[(nglv*nglw*(m - 1)) + (nglw*(n - 1)) + l - 1] = rfs;
							functions::FOFUVW(iuvw, m, n, l, i, j, k, e, f1);
							f1 = cAlpha[iUnkns] * f1;
							basisFunc[(nglv*nglw*(m - 1)) + (nglw*(n - 1)) + l - 1] = f1;
						}
					}
				}
			}
	}
}