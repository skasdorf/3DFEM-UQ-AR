#pragma once
#include "../structure/Domain.h"
#include "../utility/constants.h"
namespace postproc {
	std::vector<std::complex<double>> findesc(Domain& dom, int freqNumber, std::vector<double> iR);
	double findrcs(Domain& dom, std::vector<std::complex<double>> cEsc, int waveNumber);
	std::vector<double> post(Domain& dom);
	std::vector<double> postplus(Domain& dom, std::vector<std::complex<double>>& errorterm);
}