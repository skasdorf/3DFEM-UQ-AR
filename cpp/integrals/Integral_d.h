
#ifndef CPP_INTEGRAL_D_H
#define CPP_INTEGRAL_D_H

#include "../utility/functions.h"
//#include "../matricies/matrix4d.h"
//#include "../matricies/matrix3d.h"
//#include "../matricies/matrix2d.h"
#include "../functions/products.h"
#include "../structure/Element.h"
#include "../integrals/BasisEval.h"
#include <iostream>
#include <complex.h>
//void dotc(const std::vector<double>& v1, const std::vector < std::complex<double>>&v2, std::complex<double>& result) {
//	result.real(v1[1] * v2[1].real() + v1[2] * v2[2].real() + v1[3] * v2[3].real());
//	result.imag(v1[1] * v2[1].imag() + v1[2] * v2[2].imag() + v1[3] * v2[3].imag());
//}
class Integral_d {
public:
	Element *e;
	Integral_d(Element*e) {
		this->e = e;
	}
    void findD(const int& iuvwh, const int& iuvw, const int& ih, const int & jh, const int& kh, const int& i, 
		const int& j, const int& k, std::complex<double>& cIntD, const int& size1);
};


#endif //CPP_INTEGRAL_D_H
