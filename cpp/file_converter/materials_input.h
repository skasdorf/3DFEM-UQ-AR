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
	std::complex<double> epsr;
	std::complex<double> mur;
	int Kuvw;
	int sym;
	std::vector<std::string> eps_list;
	std::vector<std::string> mu_list;
	//std::vector<std::vector<std::complex<double>>> epsr_mat;
	//std::vector<std::vector<std::complex<double>>> mur_mat;

	Material(int hcode, int icode, int pmlcode, std::complex<double> epsr, std::complex<double> mur, int Kuvw, int sym) {
		this->hcode = hcode;
		this->icode = icode;
		this->pmlcode = pmlcode;
		this->epsr = epsr;
		//this->epsr_mat = epsr_mat;
		//this->mur_mat = mur_mat;
	}
	Material() {
		hcode = NULL;
	}

};
#endif
