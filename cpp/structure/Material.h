//material parameters
//Format:
/*
params format for each type specifier(here for eps) :
	0
	epsr
	1
	epsrxx epsryy epsrzz
	2
	epsrxx epsryx epsryy epsrzx epsrzy epsrzz
	3
	epsrxx epsrxy epsrxz epsryx epsryy epsryz epsrzx epsrzy epsrzz

*/
#ifndef MATERIAL_H
#define MATERIAL_H
#include <complex>
class Material {
public:
	/*unsigned int type_eps;
	unsigned int type_mu;
	std::vector<std::vector<std::complex<double>>> epsr;
	std::vector<std::vector<std::complex<double>>> mur;
	
	Material(unsigned int type_eps_spec, unsigned int type_mu_spec, std::vector < std::vector<std::complex<double>>> epsr_input, std::vector<std::vector<std::complex<double>>> mur_input) {
		type_eps = type_eps_spec;
		type_mu = type_mu_spec;
		epsr = epsr_input;
		mur = mur_input;

	}*/

	int hcode;
	int icode;
	int pmlcode;
	int region;
	//std::complex<double> epsr;
	//std::complex<double> mur;
	int Kuvw;
	int KuvwA;
	int sym;
	std::vector<std::vector<std::complex<double>>> epsr_list; //money
	std::vector<std::vector<std::complex<double>>> mur_list; //money!!

	Material(int hcode, int icode, int pmlcode,  int Kuvw, int sym, std::vector<std::vector<std::complex<double>>> epsr_list, std::vector<std::vector<std::complex<double>>> mur_list) {
		this->hcode = hcode;
		this->icode = icode;
		this->pmlcode = pmlcode;

		this->Kuvw = Kuvw;
		this->sym = sym;
		this->epsr_list = epsr_list;
		this->mur_list = mur_list;
	}
	Material() {
	}
	
};
#endif
