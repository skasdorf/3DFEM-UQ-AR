#//monostatic scattering class
//NEED TO FINISH waveSetup STILL
#ifndef SCATTER_H
#define SCATTER_H


#include "../functions/Coordinates.h"
class Scatter {

public:

	double thStart;
	double thStop;
	int nTh; //number of thetas tested
	double phStart;
	double phStop;
	int nPh; //number of phis tested
	double fStart;
	double fStop;
	int matDimCon;
	int waveNum; //should be equal to nWaves -- get rid of one of these later (waveNum is read in, nWaves is calculated here)
	double theta;
	double phi;
	std::complex<double> eTheta;
	std::complex<double> ePhi;
	int nF; //number of frequencies tested
	std::vector <int> incWaveThetaCount; //vector of thCount values -- first entry is blank -- data start from [1]
	std::vector <int> incWavePhiCount; //vector of phCount values -- first entry is blank -- data start from [1]
	std::vector <double> incWaveThetaAr; //vector of thAr values -- first entry is blank -- data start from [1]
	std::vector <double> incWavePhiAr; //vector of phAr values -- first entry is blank -- data start from [1]
	std::vector <double> incWaveFrequency; //vector of frequency values -- first entry is blank -- data start from [1]
	std::vector <double> incWaveK0; //vector of k0 values  -- first entry is blank -- data start from [1]
	std::vector <std::vector<double>> cGr; //matrix of right-hand side values -- matDimCon x nF -- data start from [1][1]
	std::vector <std::complex<double>> waveTheta;
	std::vector <std::complex<double>> wavePhi;
	std::vector <double> incWaveiThetaAr;
	std::vector <double> incWaveiPhiAr;
	std::vector <double> incWavenAr;

	void waveSetup() {
		int nWaves = nTh*nPh; //number of incident waves tested
		cGr.resize(matDimCon, std::vector<double>(nF));
		for (int waveNumber = 0; waveNumber <= nWaves; waveNumber++) {
			waveTheta[waveNumber] = eTheta;
			wavePhi[waveNumber] = ePhi;
			incWaveiThetaAr = coordinates::findiPhi(incWavePhiAr[1]);
			incWaveiPhiAr = coordinates::findiTheta(incWaveThetaAr[1], incWavePhiAr[1]);
			incWavenAr = coordinates::findnAr(incWaveThetaAr[1], incWavePhiAr[1]);
		}
	}

	void calcThetaPhi() {
		int nWaves = nTh*nPh; //number of incident waves tested
		int iCount; //index of current incident wave

					//resize incident wave information vectors to nWaves
		incWaveThetaCount.resize(nWaves + 1); //first entry is blank -- data start from [1]
		incWavePhiCount.resize(nWaves + 1); //first entry is blank -- data start from [1]
		incWaveThetaAr.resize(nWaves + 1); //first entry is blank -- data start from [1]
		incWavePhiAr.resize(nWaves + 1); //first entry is blank -- data start from [1]
										 //done resize incident wave information vectors to nWaves

		iCount = 0; //initialize iCount
					//make nPh*nth=nWaves sample points from thStart to thStop, phStart to phStop inclusive
		for (int thCount = 1; thCount <= nTh; ++thCount) {
			for (int phCount = 1; phCount <= nPh; ++phCount) {
				++iCount; //increment wave index
				incWaveThetaCount[iCount] = thCount; //assign wave its designated theta index
				incWavePhiCount[iCount] = phCount; //assign wave its designated phi index

				if (thCount == 1) {
					incWaveThetaAr[iCount] = thStart; //first theta increment
				}
				else {
					incWaveThetaAr[iCount] = thStart + (thCount - 1)*(thStop - thStart) / (nTh - 1); //next theta increment
				}
				if (phCount == 1) {
					incWavePhiAr[iCount] = phStart; //first phi increment
				}
				else {
					incWavePhiAr[iCount] = phStart + (phCount - 1)*(phStop - phStart) / (nPh - 1); //next phi increment
				}

			} //for phCount
		} //for thCount
		  //done make nPh*nth=nWaves sample points from thStart to thStop, phStart to phStop inclusive
	} //void calcThetaPhi

	void calcFrequencies() {
		//Preallocates frequencies and k0 values
		//Used in faster filling of G vectors
		int freqNumber; //counter for allocating frequencies
		incWaveFrequency.resize(nF + 1); //first entry is blank -- data start from [1]
		incWaveK0.resize(nF + 1); //first entry is blank -- data start from [1]

								  //make nF sample points from fStart to fStop inclusive
		incWaveFrequency[1] = fStart; //first frequency increment
		incWaveK0[1] = 2 * constants::PI*incWaveFrequency[1] * std::sqrt(constants::mu0*constants::eps0); //first k0 increment
		for (freqNumber = 2; freqNumber <= nF; freqNumber++) {
			incWaveFrequency[freqNumber] = fStart + (freqNumber - 1)*(fStop - fStart) / (nF - 1); //next frequency increment
			incWaveK0[freqNumber] = 2 * constants::PI*incWaveFrequency[freqNumber] * std::sqrt(constants::mu0*constants::eps0); //next k0 increment
		} //for freqNumber
		  //done make nF sample points from fStart to fStop inclusive

	} //void calcFrequencies
};
#endif // !SCATTER_H
