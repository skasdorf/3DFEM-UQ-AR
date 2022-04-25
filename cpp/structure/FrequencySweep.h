#pragma once
#include <vector>
#include <complex>
#include <math.h>
#include "../functions/unitVectorsM.h"
#include "../utility/constants.h"
//#include "../structure/Domain.h"
#include "../structure/Scatter.h"
class FrequencySweep {

public:
	//right now matrix reduction stuff is commented out for simplicity--re-add this later
	//double solverEpsilon = std::pow(10, -10); //solver epsilon -- anything smaller than this is considered zero

	void getCMatrix(Scatter scatter1, std::vector<std::complex<double>>&cFEMrSparse, // . . .
		std::vector<std::complex<double>>& cArSparse, std::vector<std::complex<double>>& cBrSparse, int nNz, int fCount) {
		cFEMrSparse.resize(nNz+1);
		//note that Reduce_Sparse has been combined with this
		//fCount is the current frequency to be tested
		//nNz is the number of nonzero entries in the FEM matrix -- this is the length of cFEMrSparse
		double k0 = scatter1.K0[fCount]; //a specific k0 value corresponding to the chosen frequency
		for (int i = 1; i <= nNz; i++) {
			//add A matrix entries to -k0*k0 weighted B matrix entries to form final c matrix (the matrix describing the linear system to be solved)
			cFEMrSparse[i] = cArSparse[i] - k0*k0*cBrSparse[i];
			//if (std::abs(cFEMrSparse[i]) <= solverEpsilon) {
			//	cFEMrSparse[i] = 0; //set the given entry of the matrix to zero if it is smaller than solver epsilon
			//}
		}//for int i
	}//getCMatrix



};
