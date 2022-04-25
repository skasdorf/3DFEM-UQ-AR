//monostatic scattering class
//NEED TO FINISH waveSetup STILL
#ifndef SCATTER_H
#define SCATTER_H
#include "../functions/Coordinates.h"
#include <iostream>


class Scatter {

public:


	int nWaves;
	int nTh;
	int thStart;
	int thStop;
	int nPh;
	int phStart;
	int phStop;
	int nF;
	int fStart;
	int fStop;
	int numberOfWaves;

	std::vector<int> ThetaCount;
	std::vector<int> PhiCount;

	std::vector<double> PhiAr;
	std::vector<double> ThetaAr;
	std::vector<double> K0;
	std::vector<double> Frequency;

	std::vector<std::complex<double>> waveTheta;
	std::vector<std::complex<double>> wavePhi;
	std::vector<std::vector<double>> iThetaAr;
	std::vector<std::vector<double>> iPhiAr;
	std::vector<std::vector<double>> nAr;
	
	std::vector<std::vector<std::complex<double>>> cGr;
	std::vector<std::complex<double>> cGr_eps;
	void waveSetup(int matDimCon, int nfr, int fstart, int fstop, int numberOfWaves) {
		fStart = fstart;
		fStop = fstop;
		nF = nfr;
		this->nWaves = nTh*nPh; //number of incident waves tested
		ThetaCount.resize(nWaves + 1);
		PhiCount.resize(nWaves + 1);
		PhiAr.resize(nWaves + 1);
		ThetaAr.resize(nWaves + 1);
		K0.resize(nWaves + 1);
		Frequency.resize(nWaves + 1);
		//cGr.resize(matDimCon+1, std::vector<std::complex<double>>(nF+1));
		cGr.resize(nF + 1, std::vector<std::complex<double>>(matDimCon + 1));
		cGr_eps.resize(matDimCon + 1, 0.0);
		iThetaAr.resize(nWaves+1, std::vector<double>(4));
		iPhiAr.resize(nWaves+1, std::vector<double>(4));
		nAr.resize(nWaves+1, std::vector<double>(4));
		this->numberOfWaves = numberOfWaves;


		for (int waveNumber = 0; waveNumber <= nWaves; waveNumber++) {
			iThetaAr[waveNumber] = coordinates::findiPhi(PhiAr[waveNumber]);
			iPhiAr[waveNumber] = coordinates::findiTheta(ThetaAr[waveNumber], PhiAr[waveNumber]);
			nAr[waveNumber] = coordinates::findnAr(ThetaAr[waveNumber], PhiAr[waveNumber]);
		}

		waveTheta.resize(nWaves + 1);
		wavePhi.resize(nWaves + 1);
		
	}//void waveSetup
	void calcThetaPhi() {

		this->nWaves = nTh*nPh; //number of incident waves tested
		int iCount; //index of current incident wave

		//resize incident wave information vectors to nWaves
		//ThetaCount.resize(nWaves+1); //first entry is blank -- data start from [1]
		//PhiCount.resize(nWaves+1); //first entry is blank -- data start from [1]
		//ThetaAr.resize(nWaves+1); //first entry is blank -- data start from [1]
		//PhiAr.resize(nWaves+1); //first entry is blank -- data start from [1]
		//done resize incident wave information vectors to nWaves

		iCount = 0; //initialize iCount
		//make nPh*nth=nWaves sample points from thStart to thStop, phStart to phStop inclusive
		for (int thCount = 1; thCount <= nTh; ++thCount) {
			for (int phCount = 1; phCount <= nPh; ++phCount) {
				++iCount; //increment wave index
				ThetaCount[iCount] = thCount; //assign wave its designated theta index
				PhiCount[iCount] = phCount; //assign wave its designated phi index

				if (thCount == 1) {
					ThetaAr[iCount] = thStart; //first theta increment
				}
				else {
					ThetaAr[iCount] = thStart + (thCount - 1)*(thStop - thStart) / (nTh - 1); //next theta increment
				}
				if (phCount == 1) {
					PhiAr[iCount] = phStart; //first phi increment
				}
				else {
					PhiAr[iCount] = phStart + (phCount - 1)*(phStop - phStart) / (nPh - 1); //next phi increment
				} 

			} //for phCount
		} //for thCount
		//done make nPh*nth=nWaves sample points from thStart to thStop, phStart to phStop inclusive
		for (int waveNumber = 1; waveNumber <= nWaves; ++waveNumber) {
			iPhiAr[waveNumber] = coordinates::findiPhi(PhiAr[waveNumber]);
			iThetaAr[waveNumber] = coordinates::findiTheta(ThetaAr[waveNumber], PhiAr[waveNumber]);
			nAr[waveNumber] = coordinates::findnAr(ThetaAr[waveNumber], PhiAr[waveNumber]);
		}
		/*if (!functions::_debug_matrix("../exampleFiles/debug_data/theta_count_mat.txt", ThetaCount, 1)) {
			std::cout << "ERROR! ThetaCount does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! ThetaCount matches FORTRAN!" << std::endl;
		if (!functions::_debug_matrix("../exampleFiles/debug_data/phi_count_mat.txt", PhiCount, 1)) {
			std::cout << "ERROR! PhiCount does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! PhiCount matches FORTRAN!" << std::endl;
		if (!functions::_debug_matrix("../exampleFiles/debug_data/theta_ar_mat.txt", ThetaAr, 1)) {
			std::cout << "ERROR! ThetaAr does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! incWaveThetaAr matches FORTRAN!" << std::endl;
		if (!functions::_debug_matrix("../exampleFiles/debug_data/phi_ar_mat.txt", PhiAr, 1)) {
			std::cout << "ERROR! PhiAr does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! PhiAr matches FORTRAN!" << std::endl;*/
	} //void calcThetaPhi

	void calcFrequencies() {

		//Preallocates frequencies and k0 values
		//Used in faster filling of G vectors
		int freqNumber; //counter for allocating frequencies
		Frequency.resize(nF + 1); //first entry is blank -- data start from [1]
		K0.resize(nF + 1); //first entry is blank -- data start from [1]

		//make nF sample points from fStart to fStop inclusive
		Frequency[1] = fStart; //first frequency increment
		K0[1] = 2 * constants::PI*Frequency[1] * std::sqrt(constants::mu0*constants::eps0); //first k0 increment
		for (freqNumber = 2; freqNumber <= nF; freqNumber++) {
			Frequency[freqNumber] = fStart + (freqNumber - 1)*(fStop - fStart) / (nF - 1); //next frequency increment
			K0[freqNumber] = 2 * constants::PI*Frequency[freqNumber] * std::sqrt(constants::mu0*constants::eps0); //next k0 increment
		} //for freqNumber
		//done make nF sample points from fStart to fStop inclusive

		//if (!functions::_debug_matrix("../exampleFiles/debug_data/frequency_list_mat.txt", Frequency, 1)) {
		//	std::cout << "ERROR! Frequency does not match FORTRAN!" << std::endl;
		//}
		//else std::cout << "Success! Frequency matches FORTRAN!" << std::endl;

		/*if (!functions::_debug_matrix("../exampleFiles/debug_data/k0_list_mat.txt", K0, 1)) {
			std::cout << "ERROR! K0 does not match FORTRAN!" << std::endl;
		}
		else std::cout << "Success! K0 matches FORTRAN!" << std::endl;*/


	} //void calcFrequencies
};
#endif // !SCATTER_H