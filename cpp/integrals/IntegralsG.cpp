#include "Integral_g.h"

class IntegralsG {

public:
	//right now matrix reduction stuff is commented out for simplicity--re-add this later
	//double solverEpsilon = std::pow(10, -10); //solver epsilon -- anything smaller than this is considered zero
	std::vector<std::complex<double>> cFF;
	int nglu;
	int nglv;
	int nglw;
	std::vector<std::complex<double>>wglu, wglw, wglv;
	matrix4d<std::complex<double>> EpsRelInt;
	matrix4d<std::complex<double>> MuRelIntInv;
	matrix4d<std::complex<double>> MuRelInt;
	matrix2d<double> tuPowers, fvPowers, fwPowers, fpvPowers, fpwPowers, uPowers, vPowers, fuPowers, wPowers, fpuPowers;

	void findGWave(Domain &dom, int iuvwh, int iuvw, int ih, int jh, int kh, int i, int j, int k, int size1,
		matrix4d<double> auMatrix, matrix4d<double> avMatrix,
		matrix4d<double> awMatrix, std::complex<double> czero, int nglu, int nglv, int nglw, int aCode,
		int iHomCode, int kMat, int ihCon, matrix4d<double> rMatrix){

		
		cFF.resize(3);
		std::complex<double> cSolGtotal = 0;
		
		for (int waveNumber = 1; waveNumber <= dom.scatter1.nWaves; waveNumber++) {

			//get particular scattering values specific to (waveNumber)th wave
			//particular values start with w for distinction from lists for all waves i.e. wnAr is a particular nAr 
			std::complex<double> wcEtheta = dom.scatter1.waveTheta[waveNumber];
			std::complex<double> wcEphi = dom.scatter1.wavePhi[waveNumber];
			std::vector<double> wiPhiAr = dom.scatter1.iPhiAr[waveNumber];
			std::vector<double> wiThetaAr = dom.scatter1.iThetaAr[waveNumber];
			std::vector<double> wnAr = dom.scatter1.nAr[waveNumber];
			std::vector<double> wiR;
			wiR = wnAr;
			for (int i = 0; i < 3; i++) {
				wiR[i] = -1 * wiR[i]; //wiR = -wnAr
			}

			for (int freqNumber = 1; freqNumber <= dom.scatter1.nF; freqNumber++) {

				if (dom.sc.useAdjoint) {
					//cFF = findFF(iuvwh, ih, jh, kh, face, k0List[freqNumber], wiR); //ADD ME
					if (dom.vectorD[kMat] > 0) {
						dom.scatter1.cGr[ihCon][freqNumber] = dom.scatter1.cGr[ihCon][freqNumber] + std::conj(cFF[2]);
					}
					else {
						dom.scatter1.cGr[ihCon][freqNumber] = dom.scatter1.cGr[ihCon][freqNumber] - std::conj(cFF[2]);
					}
				}// if dom.sc.useAdjoint
				else { //not adjoint

					cSolGtotal = findGInc(dom, iuvwh, iuvw, ih,jh, kh,  i,  j, k,  size1,
						auMatrix, avMatrix, awMatrix,
						nglu, nglv, nglw, aCode, iHomCode, // . . .
						rMatrix, wnAr, wcEtheta,// . . .
						wcEphi, wiPhiAr,  wiThetaAr, dom.scatter1.K0[freqNumber]); 

					if (cSolGtotal != 0.0) {
						if (dom.vectorD[kMat] > 0) {
							dom.scatter1.cGr[ihCon][freqNumber] = dom.scatter1.cGr[ihCon][freqNumber] - dom.scatter1.K0[freqNumber]* dom.scatter1.K0[freqNumber]*cSolGtotal;
						}
						else {
							dom.scatter1.cGr[ihCon][freqNumber] = dom.scatter1.cGr[ihCon][freqNumber] + dom.scatter1.K0[freqNumber] * dom.scatter1.K0[freqNumber] * cSolGtotal;
						}
					}
				}// else
			}// for int freqNumber
		}// for int waveNumber
		return;
	}//findGWave


	std::complex<double> findGInc(Domain &dom, int iuvwh, int iuvw, int ih, int jh, int kh, int i, int j, int k, int size1,
		matrix4d<double> auMatrix, matrix4d<double> avMatrix, matrix4d<double> awMatrix,
		int nglu, int nglv, int nglw, int aCode, int iHomCode, // . . .
		matrix4d<double>rMatrix, std::vector<double> wnAr, std::complex<double> wcEtheta,// . . .
		std::complex<double> wcEphi, std::vector<double> wiPhiAr, std::vector<double> wiThetaAr, double k0) {

		std::complex<double> cIntDInc = 0;
		double jacobian;
		std::vector<double> au, av, aw;
		double f1;
		std::vector<std::complex<double>> f3, asec2, rfus1(3), asec1(3), f4;
		std::complex<double> f2, f;
		std::complex<double> cEpsilon, cMu;
		std::vector<std::complex<double>> eVector, hVector, asec3;
		std::complex<double> cj = std::sqrt(-1);
		int m, n, l;
		std::complex<double> cSn, cSl;
		int coord;
		
		for (m = 1; m <= nglu; m++) {
			cSn = 0;
			for (n = 1; n <= nglv; n++) {
				cSl = 0;
				for (l = 1; l <= nglw; l++) {
					for (coord = 1; coord <= 3; coord++) {
						au[coord] = auMatrix(m, n, l, coord);
						av[coord] = avMatrix(m, n, l, coord);
						aw[coord] = awMatrix(m, n, l, coord);
					}//for coord

					//The D Part
//ADD ME					//f1 = FOFUVW(iuvwh, m, n, l, ih, jh, kh); //ADD ME
//ADD ME					//asec1 = FINDSECONDARYUNIT(iuvwh, au, av, aw); ADD ME
					std::vector<double> r = rMatrix(m, n, l);
					double rn = products::dot(r, wnAr);
					asec2[0] = (wcEphi*wiPhiAr[0] + wcEtheta*wiThetaAr[0])*std::exp(-cj*k0*rn);
					asec2[1] = (wcEphi*wiPhiAr[1] + wcEtheta*wiThetaAr[1])*std::exp(-cj*k0*rn);
					asec2[2] = (wcEphi*wiPhiAr[2] + wcEtheta*wiThetaAr[2])*std::exp(-cj*k0*rn);

					if (aCode == 0) {
						if (iHomCode == 0) {
							std::vector<std::complex<double>> slicedEpsRel(size1 + 1);
							for (int var = 1; var <= size1; ++var) {
								slicedEpsRel[var] = EpsRelInt(var, m, n, l);
							}
							MULTIPLYMURMATRIX(size1, slicedEpsRel, asec2, f3);
							f3[0] = f3[0] - asec2[0];
							f3[1] = f3[1] - asec2[1];
							f3[2] = f3[2] - asec2[2];
						}
						else {
							std::vector<std::complex<double>> slicedEpsRel(size1 + 1);
							for (int var = 1; var <= size1; ++var) {
								slicedEpsRel[var] = EpsRelInt(var, 1, 1, 1);
							}
							MULTIPLYMURMATRIX(size1, slicedEpsRel, asec2, f3);
							f3[0] = f3[0] - asec2[0];
							f3[1] = f3[1] - asec2[1];
							f3[2] = f3[2] - asec2[2];
						}
					}
					else {
						if (iHomCode == 0) {
							f3[0] = asec2[0] * (EpsRelInt(1, m, n, l) - 1.0);
							f3[1] = asec2[1] * (EpsRelInt(1, m, n, l) - 1.0);
							f3[2] = asec2[2] * (EpsRelInt(1, m, n, l) - 1.0);
						}
						else {
							f3[0] = asec2[0] * (EpsRelInt(1, 1, 1, 1) - 1.0);
							f3[1] = asec2[1] * (EpsRelInt(1, 1, 1, 1) - 1.0);
							f3[2] = asec2[2] * (EpsRelInt(1, 1, 1, 1) - 1.0);
						}
					}//else (aCode!=0)
					f = products::dot(asec1, f3);

					//the C part
					ROTFSEC(iuvwh, m, n, l, ih, jh, kh, au, av, aw, rfus1);

					//From Maxwell's equations, magnitude of curl(E) is equal to magnitude of E times k, direction
					//equal to direction of H, phase shift of - j
					eVector[0] = wcEphi*wiPhiAr[0] + wcEtheta*wiThetaAr[0];
					eVector[1] = wcEphi*wiPhiAr[1] + wcEtheta*wiThetaAr[1];
					eVector[2] = wcEphi*wiPhiAr[2] + wcEtheta*wiThetaAr[2];
					std::vector <std::complex<double>> cmpxwnAr;
					cmpxwnAr[0] = wnAr[0];
					cmpxwnAr[1] = wnAr[1];
					cmpxwnAr[2] = wnAr[2];
					hVector = products::cross(cmpxwnAr, eVector);
					asec3[0] = -cj*hVector[0] * std::exp(-cj*k0*rn) / k0;
					asec3[1] = -cj*hVector[1] * std::exp(-cj*k0*rn) / k0;
					asec3[2] = -cj*hVector[2] * std::exp(-cj*k0*rn) / k0;

					if (aCode == 0) {
						if (iHomCode == 0) {
							std::vector<std::complex<double>> slicedMuRel(size1 + 1);
							for (int var = 1; var <= size1; ++var) {
								slicedMuRel[var] = MuRelIntInv(var, m, n, l);
							}
							MULTIPLYMURMATRIX(size1, slicedMuRel, asec3, f4);
							f4[0] = f4[0] - asec3[0];
							f4[1] = f4[1] - asec3[1];
							f4[2] = f4[2] - asec3[2];
						}
						else {
							std::vector<std::complex<double>> slicedMuRel(size1 + 1);
							for (int var = 1; var <= size1; ++var) {
								slicedMuRel[var] = MuRelIntInv(var, 1, 1, 1);
							}
							MULTIPLYMURMATRIX(size1, slicedMuRel, asec3, f4);
							f4[0] = f4[0] - asec3[0];
							f4[1] = f4[1] - asec3[1];
							f4[2] = f4[2] - asec3[2];
						}
					}
					else {
						if (iHomCode == 0) {
							f4[0] = asec3[0] / (MuRelInt(1, m, n, l) - asec3[0]);
							f4[1] = asec3[1] / (MuRelInt(1, m, n, l) - asec3[1]);
							f4[2] = asec3[2] / (MuRelInt(1, m, n, l) - asec3[2]);
						}
						else {
							f4[0] = asec3[0] / (MuRelInt(1, 1, 1, 1) - asec3[0]);
							f4[1] = asec3[1] / (MuRelInt(1, 1, 1, 1) - asec3[1]); 
							f4[2] = asec3[2] / (MuRelInt(1, 1, 1, 1) - asec3[2]);
						}
					}//else (aCode!=0)

					f2 = products::dot(rfus1, f4);

					cSl = cSl + wglw[l] * f2;
					cSl = cSl + wglw[l] * f1*f;
				}//for l
				cSn = cSn + wglv[n] * cSl;
			}//for n
			cIntDInc = cIntDInc + wglu[m] * cSn;
		}//for m
		return cIntDInc;
	}//findGInc


	void IntegralsG::MULTIPLYMURMATRIX(int size1, std::vector<std::complex<double>> MuRelInv,
		std::vector<std::complex<double>>  vector, std::vector<std::complex<double>> vectorNew) {

		switch (size1) {
		case 3:
			vectorNew[1] = vector[1] * MuRelInv[1];
			vectorNew[2] = vector[2] * MuRelInv[2];
			vectorNew[3] = vector[3] * MuRelInv[3];
			break;
		case 6:
			vectorNew[1] = vector[1] * MuRelInv[1] + vector[2] * MuRelInv[2] + vector[3] * MuRelInv[3];
			vectorNew[2] = vector[1] * MuRelInv[2] + vector[2] * MuRelInv[4] + vector[3] * MuRelInv[5];
			vectorNew[3] = vector[1] * MuRelInv[3] + vector[2] * MuRelInv[5] + vector[3] * MuRelInv[6];
			break;
		case 9:
			vectorNew[1] = vector[1] * MuRelInv[1] + vector[2] * MuRelInv[2] + vector[3] * MuRelInv[3];
			vectorNew[2] = vector[1] * MuRelInv[4] + vector[2] * MuRelInv[5] + vector[3] * MuRelInv[6];
			vectorNew[3] = vector[1] * MuRelInv[7] + vector[2] * MuRelInv[8] + vector[3] * MuRelInv[9];
			break;
		}
	}

	void IntegralsG::ROTFSEC(int iuvw, int m, int n, int l, int i, int j, int k, std::vector<double> au, std::vector<double> av,
		std::vector<double> aw, std::vector<std::complex<double>>  rfs) {
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

	void IntegralsG::ROTFUSEC(int m, int n, int l, int i, int j, int k, std::vector<double> av, std::vector<double> aw,
		std::vector<std::complex<double>>  rfus) {

		double mno1 = uPowers(m, i) * fvPowers(n, j)  * fpwPowers(l, k);
		double mno2 = uPowers(m, i) * fpvPowers(n, j) * fwPowers(l, k);

		for (int cord = 1; cord <= 3; ++cord) {
			rfus[cord] = mno1 * av[cord] - mno2 * aw[cord];
		}

	}

	void IntegralsG::ROTFVSEC(int m, int n, int l, int i, int j, int k, std::vector<double> aw, std::vector<double> au,
		std::vector<std::complex<double>>  rfvs) {
		double mno1 = fpuPowers(m, i) * vPowers(n, j) * fwPowers(l, k);
		double mno2 = fuPowers(m, i)  * vPowers(n, j) * fpwPowers(l, k);

		for (int cord = 1; cord <= 3; ++cord) {
			rfvs[cord] = mno1 * aw[cord] - mno2 * au[cord];
		}
	}

	void IntegralsG::ROTFWSEC(int m, int n, int l, int i, int j, int k, std::vector<double> au, std::vector<double> av,
		std::vector<std::complex<double>>  rfws) {
		double mno1 = fuPowers(m, i)  * fpvPowers(n, j) * wPowers(l, k);
		double mno2 = fpuPowers(m, i) * fvPowers(n, j)  * wPowers(l, k);

		for (int cord = 1; cord <= 3; ++cord) {
			rfws[cord] = mno1 * au[cord] - mno2 * av[cord];
		}
	}

};
