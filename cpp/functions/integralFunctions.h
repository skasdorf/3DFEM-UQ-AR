//
// Created by Blake Troksa on 6/5/18.
//

#ifndef CPP_INTEGRALFUNCTIONS_H
#define CPP_INTEGRALFUNCTIONS_H

#include "products.h"
#include "../utility/functions.h"
#include "../structure/Element.h"
#include <iostream>

namespace functions {
    void FINDSECONDARYUNIT(const int& iuvw, const std::vector<double>& a1, const std::vector<double>& a2,
                            std::vector<double>& asec);
	void FINDSECONDARYUNIT(const int& iuvw, const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double>& asec);
	//template <class T>
    void FOFUVW(const int& iuvw, const int &m, const int& n, const int& l, const int& i, const int& j, const int& k, Element* e, double& fuvw);
	void FOFUVW(const int& iuvw, const int &m, const int& n, const int& l, const int& i, const int& j, const int& k, Element* e, std::complex<double>& fuvw);

	//std::vector<std::complex<double>> MULTIPLYMURMATRIX(const int& size1, const std::vector<std::complex<double>>& MuRelInv, const std::vector<double>& vector);
	void MULTIPLYMURMATRIX(const int& size1, const matrix4d<std::complex<double>>& MuRelInv,
		const std::vector<double>& vector, const int& m, const int& n, const int& l, std::vector<std::complex<double>>& out);
};


#endif //CPP_INTEGRALFUNCTIONS_H
