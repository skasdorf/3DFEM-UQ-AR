
#ifndef CPP_INTEGRAL_C_H
#define CPP_INTEGRAL_C_H

#include "../utility/functions.h"
//#include "../matricies/matrix4d.h"
//#include "../matricies/matrix3d.h"
//#include "../matricies/matrix2d.h"
#include "../functions/products.h"
#include "../functions/integralFunctions.h"
#include <iostream>
#include <complex.h>
#include "../structure/Element.h"
#include "../integrals/BasisEval.h"
//#include "../structure/Point.h"
//template <class T>
//std::vector<T> operator*(const std::vector<T>& vec, const T& scalar) {
//	//assumes starts at 1
//	return{ 0, vec[1] * scalar, vec[2] * scalar, vec[3] * scalar };
//}
//template <class T>
//std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2) {
//	//assumes that the vectors both start at 1
//	return{ 0, vec1[1] + vec2[1], vec1[2] + vec2[2], vec1[3] + vec2[3] };
//}
//void dotc(const std::vector<double>& v1, const std::vector < std::complex<double>>&v2, std::complex<double>& result) {
//	result.real(v1[1] * v2[1].real() + v1[2] * v2[2].real() + v1[3] * v2[3].real());
//	result.imag(v1[1] * v2[1].imag() + v1[2] * v2[2].imag() + v1[3] * v2[3].imag());
//}
class Integral_c {
public:
	Element *e;
	Integral_c(Element* e);

    void findC(const int& iuvwh, const int& iuvw, const int& ih, const int& jh, const int& kh, const int& i, const int& j,
		const int& k, std::complex<double>& cIntC, const int& size1, BasisEval& eval);
	void ROTFUSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double> &rfus, BasisEval& eval);
	void ROTFVSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& aw, const std::vector<double>& au,
		std::vector<double> &rfvs, BasisEval& eval);
	void ROTFWSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
		std::vector<double> &rfws, BasisEval& eval);
	void ROTFSEC(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
		const std::vector<double>& aw, std::vector<double> &rfs, BasisEval& eval);
    void ROTFUSEC(const double& mno1, const double& mno2, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double> &rfus, BasisEval& eval);

    void ROTFVSEC(const double& mno1, const double& mno2, const std::vector<double>& aw, const std::vector<double>& au,
		std::vector<double> &rfvs, BasisEval& eval);

    void ROTFWSEC(const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
		std::vector<double> &rfws, BasisEval& eval);

    void ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
		const std::vector<double>& aw, std::vector<double> &rfs, BasisEval& eval);

	void ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& a1, const std::vector<double>& a2,
		std::vector<double> &rfs, BasisEval& eval);
};


#endif //CPP_INTEGRAL_C_H
