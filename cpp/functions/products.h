#ifndef PRODUCTS_H
#define PRODUCTS_H

#include "../utility/constants.h"
namespace products {
	inline void dotc(const std::vector<double>& v1, const std::vector < std::complex<double>>&v2, std::complex<double>& result) {
		result.real(v1[1] * v2[1].real() + v1[2] * v2[2].real() + v1[3] * v2[3].real());
		result.imag(v1[1] * v2[1].imag() + v1[2] * v2[2].imag() + v1[3] * v2[3].imag());
	}
	template<typename T, typename T2>
	inline void cross(const std::vector<T>& v1, const std::vector<T2>& v2, std::vector<T2>& v3) {
		//input and output vectors start at 1, NOT 0
		//std::vector<T2> v3(4);
		v3[1] = v1[2] * v2[3] - v2[2] * v1[3];
		v3[2] = v1[3] * v2[1] - v2[3] * v1[1];
		v3[3] = v1[1] * v2[2] - v2[1] * v1[2];
		
	}


	template<typename T, typename T2>
	inline void dot(const std::vector<T>& v1, const std::vector<T2>& v2, T2& result) {
	
	
			result = ((v1[1] * v2[1]) + (v1[2] * v2[2])) + (v1[3] * v2[3]);
		
		
	}
	
	//void dot(const std::vector<double>& v1, const std::vector<std::complex<double>> v2, std::complex<double>& result) {
	//	result.real(v1[1] * v2[1].real() + v1[2] * v2[2].real() + v1[3] * v2[3].real());
	//	result.imag(v1[1] * v2[1].imag() + v1[2] * v2[2].imag() + v1[3] * v2[3].imag());
	//}
	//void dot(const std::vector<double>& v1, const std::vector<double>& v2, double& result) {
	//	result = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
	//}
	/*void dot(const std::vector<std::complex<double>>& v1, const std::vector<std::complex<double>> v2, std::complex<double>& result) {
		result.real(v1[1].real() * v2[1].real() + v1[2].real() * v2[2].real() + v1[3].real() * v2[3].real()
			-v1[1].imag()*v2[1].imag() - v1[2].imag()*v2[2].imag() - v1[3].imag()*v2[3].imag());
		result.imag(v1[1].real() * v2[1].imag() + v1[2].real() * v2[2].imag() + v1[3].real() * v2[3].imag()
			+v1[1].imag()*v2[1].real() + v1[2].imag()*v2[2].real() + v1[3].imag()*v2[3].real());
	}*/
}

#endif