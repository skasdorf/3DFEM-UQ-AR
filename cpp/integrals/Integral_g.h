//#ifndef CPP_INTEGRAL_G_H
//#define CPP_INTEGRAL_G_H
#pragma once
#include "../utility/functions.h"
#include "../functions/integralFunctions.h"
#include "../functions/products.h"
#include <iostream>
#include <complex>
#include "../structure/Scatter.h"
#include "../structure/Element.h"
#include "../integrals/BasisEval.h"
//#include "../structure/Domain.h"


class Integral_g {
public:

	
		

		std::vector<std::complex<double>> cFF;
		
		Element *e;
		Integral_g(Element*e) {
			this->e = e;
		}
		void set_e(Element*e) {
			this->e = e;
		}

		
		void findFF(const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& face, const double& k0, const std::vector<double>& iR, std::vector < std::complex<double>>& cFF);

	void findGWave(Scatter*scatter1, const std::vector<int>& vectorD, const int& useAdjoint, const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
				    const int& kMat, const int& iCon, const int& face);

	void findGInc(Scatter*scatter1, const std::vector<int>& vectorD, const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
		 const int& waveNumber, const int& freqNumber, std::complex<double>& cIntDInc, std::vector<std::complex<double>>& asec2, std::vector<std::complex<double>>& hVector, std::vector<std::complex<double>>& asec3, double val1,std::complex<double>& cIntDInc_eps);

	void ROTFSEC(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw,
		 std::vector<double>&  rfs);

	void ROTFUSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double>&  rfus);

	void ROTFVSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& aw, const std::vector<double>& au,
		std::vector<double>&  rfvs);

	void ROTFWSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
		std::vector<double>& rfws);


	void ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& a1, const std::vector<double>& a2,
		 std::vector<double>&  rfs);

	void ROTFUSEC(const double& mno1, const double& mno2, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double>&  rfus);

	void ROTFVSEC(const double& mno1, const double& mno2, const std::vector<double>& aw, const std::vector<double>& au,
		std::vector<double>&  rfvs);

	void ROTFWSEC(const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
		std::vector<double>& rfws);

	//std::vector<std::complex<double>> MULTIPLYMURMATRIX(const int& size1, const std::vector<std::complex<double>>& MuRelInv, const std::vector<std::complex<double>>&);
	void MULTIPLYMURMATRIX(const int& size1, const matrix4d<std::complex<double>>& MuRelInv,
		const std::vector<std::complex<double>>&  vector, const int& m, const int& n, const int& l, std::vector<std::complex<double>>& out);
	void findGInc_EpsOnly(Scatter* scatter1, const std::vector<int>& vectorD,  const int& iuvwh, const int& ih, const int& jh, const int& kh, const int& size1,
		const int& waveNumber, const int& freqNumber, std::complex<double>& cIntDInc, std::vector<std::complex<double>>& asec2, std::vector<std::complex<double>>& hVector, std::vector<std::complex<double>>& asec3, double val1);
};

//#endif //CPP_INTEGRAL_G_H