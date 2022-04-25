//
// Created by Blake Troksa on 6/5/18.
//

#include "integralFunctions.h"
#include "../structure/Element.h"
#include <chrono>

namespace functions {

    void FINDSECONDARYUNIT(const int& iuvw, const std::vector<double>& a1, const std::vector<double>& a2,
                                                             std::vector<double>& asec) {
        //std::vector<double> asec;
        switch (iuvw) {
            case 1:
               // asec = products::cross(av, aw);
				products::cross(a1, a2, asec);
				break;
            case 2:
                //asec = products::cross(aw, au);
				products::cross(a1, a2, asec);
				break;
            case 3:
                //asec = products::cross(au, av);
				products::cross(a1, a2, asec);
				break;
            default:
                std::cerr << "Error: Invalid direction in FINDSECONDARYUNIT" << std::endl;
                exit(EXIT_FAILURE);
				break;
        }

       // return asec;
    }
	void FINDSECONDARYUNIT(const int& iuvw, const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw,
		std::vector<double>& asec) {
		//std::vector<double> asec;
		switch (iuvw) {
		case 1:
			// asec = products::cross(av, aw);
			products::cross(av, aw, asec);
			break;
		case 2:
			//asec = products::cross(aw, au);
			products::cross(aw, au, asec);
			break;
		case 3:
			//asec = products::cross(au, av);
			products::cross(au, av, asec);
			break;
		default:
			std::cerr << "Error: Invalid direction in FINDSECONDARYUNIT" << std::endl;
			exit(EXIT_FAILURE);
			break;
		}

		// return asec;
	}
	//template <class T>
    void FOFUVW(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, Element *e, double& fuvw) {

        //double fuvw;
        switch (iuvw) {
            case 1:
                fuvw = e->uPowers(m, i) * e->fvPowers(n, j) * e->fwPowers(l, k);
                break;
            case 2:
                fuvw = e->fuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
                break;
            case 3:
                fuvw = e->fuPowers(m, i) * e->fvPowers(n, j) * e->wPowers(l, k);
                break;
            default:
                std::cerr << "Invalid basis/testing direction in FOFUVW!";
                exit(8);
        }
       // return fuvw;
    }
	void FOFUVW(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, Element *e, std::complex<double>& fuvw) {

		//double fuvw;
		switch (iuvw) {
		case 1:
			fuvw = e->uPowers(m, i) * e->fvPowers(n, j) * e->fwPowers(l, k);
			break;
		case 2:
			fuvw = e->fuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
			break;
		case 3:
			fuvw = e->fuPowers(m, i) * e->fvPowers(n, j) * e->wPowers(l, k);
			break;
		default:
			std::cerr << "Invalid basis/testing direction in FOFUVW!";
			exit(8);
		}
		// return fuvw;
	}

	//void FOFUVW(const int& iuvw, const double& fuvw_part, double& fuvw) {

	//	//double fuvw;
	//	switch (iuvw) {
	//	case 1:
	//		fuvw = e->uPowers(m, i) * e->fvPowers(n, j) * e->fwPowers(l, k);
	//		break;
	//	case 2:
	//		fuvw = e->fuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	//		break;
	//	case 3:
	//		fuvw = e->fuPowers(m, i) * e->fvPowers(n, j) * e->wPowers(l, k);
	//		break;
	//	default:
	//		std::cerr << "Invalid basis/testing direction in FOFUVW!";
	//		exit(8);
	//	}
	//	// return fuvw;
	//}
	//void FOFUVW(const int& iuvw,const double& fuvw_part, std::complex<double>& fuvw) {

	//	//double fuvw;
	//	switch (iuvw) {
	//	case 1:
	//		fuvw = fuvw_part * e->fwPowers(l, k);
	//		break;
	//	case 2:
	//		fuvw = e->fuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	//		break;
	//	case 3:
	//		fuvw = e->fuPowers(m, i) * e->fvPowers(n, j) * e->wPowers(l, k);
	//		break;
	//	default:
	//		std::cerr << "Invalid basis/testing direction in FOFUVW!";
	//		exit(8);
	//	}
	//	// return fuvw;
	//}

	//void MULTIPLYMURMATRIX(const int& size1, const std::vector<std::complex<double>>& MuRelInv,
	//	const std::vector<double>& vector, std::vector<std::complex<double>>& out) {
	//	switch (size1) {
	//	case 9:
	//		out[1] = vector[1] * MuRelInv[1] + vector[2] * MuRelInv[2] + vector[3] * MuRelInv[3];
	//		out[
	//	}
	//	
	//	
	//	
	//	/*switch (size1) {
	//	case 9:
	//		return{ 0, vector[1] * MuRelInv[1] + vector[2] * MuRelInv[2] + vector[3] * MuRelInv[3],
	//			vector[1] * MuRelInv[4] + vector[2] * MuRelInv[5] + vector[3] * MuRelInv[6],
	//			vector[1] * MuRelInv[7] + vector[2] * MuRelInv[8] + vector[3] * MuRelInv[9] };
	//		break;

	//	case 6:
	//		return{ 0,vector[1] * MuRelInv[1] + vector[2] * MuRelInv[2] + vector[3] * MuRelInv[3],
	//			vector[1] * MuRelInv[2] + vector[2] * MuRelInv[4] + vector[3] * MuRelInv[5] ,
	//			vector[1] * MuRelInv[3] + vector[2] * MuRelInv[5] + vector[3] * MuRelInv[6] };
	//		break;
	//	case 3: 
	//		return{ 0,vector[1] * MuRelInv[1] , vector[2] * MuRelInv[2],vector[3] * MuRelInv[3] };
	//		break;
	//	default:
	//		return{ 0 };
	//		break;
	//	}*/

	//	
	//}
	void MULTIPLYMURMATRIX(const int& size1, const matrix4d<std::complex<double>>& MuRelInv,
		const std::vector<double>& vector, const int& m, const int& n, const int& l, std::vector<std::complex<double>>& out) {
		switch (size1) {
		case (9) :
			/*out[1] = vector[1] * MuRelInv(1, m, n, l) + vector[2] * MuRelInv(2, m, n, l) + vector[3] * MuRelInv(3, m, n, l);
			out[2] = vector[1] * MuRelInv(4, m, n, l) + vector[2] * MuRelInv(5, m, n, l) + vector[3] * MuRelInv(6, m, n, l);
			out[3] = vector[1] * MuRelInv(7, m, n, l) + vector[2] * MuRelInv(8, m, n, l) + vector[3] * MuRelInv(9, m, n, l);*/
			out[1].real(vector[1] * MuRelInv(1, m, n, l).real() + vector[2] * MuRelInv(2, m, n, l).real() + vector[3] * MuRelInv(3, m, n, l).real());
			out[1].imag(vector[1] * MuRelInv(1, m, n, l).imag() + vector[2] * MuRelInv(2, m, n, l).imag() + vector[3] * MuRelInv(3, m, n, l).imag());
			out[2].real(vector[1] * MuRelInv(4, m, n, l).real() + vector[2] * MuRelInv(5, m, n, l).real() + vector[3] * MuRelInv(6, m, n, l).real());
			out[2].imag(vector[1] * MuRelInv(4, m, n, l).imag() + vector[2] * MuRelInv(5, m, n, l).imag() + vector[3] * MuRelInv(6, m, n, l).imag());
			out[3].real(vector[1]*MuRelInv(7, m, n, l).real() + vector[2] * MuRelInv(8, m, n, l).real() + vector[3] * MuRelInv(9, m, n, l).real());
			out[3].imag(vector[1] * MuRelInv(7, m, n, l).imag() + vector[2] * MuRelInv(8, m, n, l).imag() + vector[3] * MuRelInv(9, m, n, l).imag());

			break;
		case 6:
			out[1] = vector[1] * MuRelInv(1, m, n, l) + vector[2] * MuRelInv(2, m, n, l) + vector[3] * MuRelInv(3, m, n, l);
			out[2] = vector[1] * MuRelInv(2, m, n, l) + vector[2] * MuRelInv(4, m, n, l) + vector[3] * MuRelInv(5, m, n, l);
			out[3] = vector[1] * MuRelInv(3, m, n, l) + vector[2] * MuRelInv(5, m, n, l) + vector[3] * MuRelInv(6, m, n, l);
			break;
		case 3:
			out[1] = vector[1] * MuRelInv(1, m, n, l);
			out[2] = vector[2] * MuRelInv(2, m, n, l);
			out[3] = vector[3] * MuRelInv(3, m, n, l);
			break;
		default:
			
			break;
		}
		/*switch (size1) {
		case 9:
			return{ 0, vector[1] * MuRelInv(1,m,n,l) + vector[2] * MuRelInv(2,m,n,l) + vector[3] * MuRelInv(3,m,n,l),
				vector[1] * MuRelInv(4,m,n,l) + vector[2] * MuRelInv(5,m,n,l) + vector[3] * MuRelInv(6,m,n,l),
				vector[1] * MuRelInv(7,m,n,l) + vector[2] * MuRelInv(8,m,n,l) + vector[3] * MuRelInv(9,m,n,l) };
			break;
		case 6:
			return{ 0, vector[1] * MuRelInv(1,m,n,l) + vector[2] * MuRelInv(2,m,n,l) + vector[3] * MuRelInv(3,m,n,l),
			vector[1] * MuRelInv(2,m,n,l) + vector[2] * MuRelInv(4,m,n,l) + vector[3] * MuRelInv(5,m,n,l),
			vector[1] * MuRelInv(3,m,n,l) + vector[2] * MuRelInv(5,m,n,l) + vector[3] * MuRelInv(6,m,n,l) };
		case 3:
			return{ 0, vector[1] * MuRelInv(1,m,n,l), vector[2] * MuRelInv(2,m,n,l), vector[3] * MuRelInv(3,m,n,l) };
			break;
		default:
			return{ 0 };
			break;
		}*/
	}
}