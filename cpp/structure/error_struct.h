#include <complex>
//e2: element that the val from e1 was merged to
//e1: element that had its val sent to e2
//i0, j0: location in e2 (of cAr,cBrPAK) where merge happened
//val_e2: val of the e2 portion before e1 part was added
//nNz: index of the nNz where the first occurence of the new happens
class Error_struct {
public:
	int e2;
	int e1;
	int i0;
	int j0;
	std::complex<double> cAr_val;
	std::complex<double> cBr_val;
	int nNz;
	Error_struct(int e2, int e1, int i0, int j0, std::complex<double> cAr_val, std::complex<double> cBr_val, int nNz) {
		this->e2 = e2;
		this->e1 = e1;
		this->i0 = i0;
		this->j0 = j0;
		this->cAr_val = cAr_val;
		this->cBr_val = cBr_val;
		this->nNz = nNz;
	}
};