
#include "Integral_c.h"


//Integral_c::Integral_c(matrix2d<double> fvPowers, matrix2d<double> fwPowers, matrix2d<double> fpvPowers, matrix2d<double> fpwPowers,
//                       matrix2d<double> uPowers, matrix2d<double> vPowers, matrix2d<double> fuPowers, matrix2d<double> wPowers,
//                       matrix2d<double> fpuPowers){
//    this->fvPowers = fvPowers;
//    this->fwPowers = fwPowers;
//    this->fpvPowers = fpvPowers;
//    this->fpwPowers = fpwPowers;
//    this->uPowers = uPowers;
//    this->vPowers = vPowers;
//    this->fuPowers = fuPowers;
//    this->wPowers = wPowers;
//    this->fpuPowers = fpuPowers;
//}


Integral_c::Integral_c(Element * e)
{
	this->e = e;
}
void Integral_c::ROTFUSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& av, const std::vector<double>& aw,
	std::vector<double> &rfus, BasisEval& eval) {

	double mno1 = e->uPowers(m, i) * e->fvPowers(n, j)  * e->fpwPowers(l, k);
	double mno2 = e->uPowers(m, i) * e->fpvPowers(n, j) * e->fwPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfus[cord] = mno1 * av[cord] - mno2 * aw[cord];
	}
	//rfus = av*mno1 + aw*(-mno2);

}

void Integral_c::ROTFVSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& aw, const std::vector<double>& au,
	std::vector<double> &rfvs, BasisEval& eval) {
	double mno1 = e->fpuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	double mno2 = e->fuPowers(m, i)  * e->vPowers(n, j) * e->fpwPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfvs[cord] = mno1 * aw[cord] - mno2 * au[cord];
	}
	//rfvs = aw*mno1 + au*(-mno2);
}

void Integral_c::ROTFWSEC(const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
	std::vector<double> &rfws, BasisEval& eval) {
	double mno1 = e->fuPowers(m, i)  * e->fpvPowers(n, j) * e->wPowers(l, k);
	double mno2 = e->fpuPowers(m, i) * e->fvPowers(n, j)  * e->wPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfws[cord] = mno1 * au[cord] - mno2 * av[cord];
	}
	//rfws = au*mno1 + av*(-mno2);
}

void Integral_c::ROTFSEC(const int& iuvw, const int& m, const int& n, const int& l, const int& i, const int& j, const int& k, const std::vector<double>& au, const std::vector<double>& av,
	const std::vector<double>& aw, std::vector<double> &rfs, BasisEval& eval) {
	switch (iuvw) {
	case 1:
		ROTFUSEC(m, n, l, i, j, k, av, aw, rfs, eval);
		break;
	case 2:
		ROTFVSEC(m, n, l, i, j, k, aw, au, rfs, eval);
		break;
	case 3:
		ROTFWSEC(m, n, l, i, j, k, au, av, rfs, eval);
		break;
	default:
		std::cerr << "Invalid basis direction in ROTSEC!" << std::endl;
		exit(9);
		break;
	}
}





void Integral_c::ROTFUSEC(const double& mno1, const double& mno2, const std::vector<double>& av, const std::vector<double>& aw,
	std::vector<double> &rfus, BasisEval& eval) {

	//double mno1 = e->uPowers(m, i) * e->fvPowers(n, j)  * e->fpwPowers(l, k);
	//double mno2 = e->uPowers(m, i) * e->fpvPowers(n, j) * e->fwPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfus[cord] = mno1 * av[cord] - mno2 * aw[cord];
	}
	//rfus = av*mno1 + aw*(-mno2);

}

void Integral_c::ROTFVSEC(const double& mno1, const double& mno2, const std::vector<double>& aw, const std::vector<double>& au,
	std::vector<double> &rfvs, BasisEval& eval) {
	//double mno1 = e->fpuPowers(m, i) * e->vPowers(n, j) * e->fwPowers(l, k);
	//double mno2 = e->fuPowers(m, i)  * e->vPowers(n, j) * e->fpwPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfvs[cord] = mno1 * aw[cord] - mno2 * au[cord];
	}
	//rfvs = aw*mno1 + au*(-mno2);
}

void Integral_c::ROTFWSEC(const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
	std::vector<double> &rfws, BasisEval& eval) {
	//double mno1 = e->fuPowers(m, i)  * e->fpvPowers(n, j) * e->wPowers(l, k);
	//double mno2 = e->fpuPowers(m, i) * e->fvPowers(n, j)  * e->wPowers(l, k);

	for (int cord = 1; cord <= 3; ++cord) {
		rfws[cord] = mno1 * au[cord] - mno2 * av[cord];
	}
	//rfws = au*mno1 + av*(-mno2);
}

void Integral_c::ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& au, const std::vector<double>& av,
	const std::vector<double>& aw, std::vector<double> &rfs, BasisEval& eval) {
	switch (iuvw) {
	case 1:
		ROTFUSEC(mno1,mno2, av, aw, rfs, eval);
		break;
	case 2:
		ROTFVSEC(mno1, mno2, aw, au, rfs, eval);
		break;
	case 3:
		ROTFWSEC(mno1, mno2, au, av, rfs, eval);
		break;
	default:
		std::cerr << "Invalid basis direction in ROTSEC!" << std::endl;
		exit(9);
		break;
	}
}

void Integral_c::ROTFSEC(const int& iuvw, const double& mno1, const double& mno2, const std::vector<double>& a1, const std::vector<double>& a2,
	std::vector<double> &rfs, BasisEval& eval) {
	switch (iuvw) {
	case 1:
		ROTFUSEC(mno1, mno2, a1, a2, rfs, eval);
		break;
	case 2:
		ROTFVSEC(mno1, mno2, a1, a2, rfs, eval);
		break;
	case 3:
		ROTFWSEC(mno1, mno2, a1, a2, rfs, eval);
		break;
	default:
		std::cerr << "Invalid basis direction in ROTSEC!" << std::endl;
		exit(9);
		break;
	}
}

void Integral_c::findC(const int& iuvwh, const int& iuvw, const int& ih, const int& jh, const int& kh, const int& i, const int& j,
	const int& k, std::complex<double>& cIntC, const int& size1, BasisEval& eval) {

    std::complex<double> cSn, cSl;
    std::vector<double> au(4), av(4), aw(4), rfus2(4), rfus1(4);

    std::vector<std::complex<double>> f1(4);
    std::complex<double> f2;
    int nglu = e->quadrature[0];
    int nglv = e->quadrature[1];
    int nglw = e->quadrature[2];
    int iHomCode = e->materials.hcode;

    cIntC = {0,0};
    for (int m = 1; m <= nglu; ++m) {
        cSn = {0,0};
		
        for (int n = 1; n <= nglv; ++n) {
            cSl = {0,0};
			
			double mno1;
			double mno2;
			double mno1_h;
			double mno2_h;
			//more efficient way of computing mno terms
			if (iuvw == 1) {
				mno1 = e->uPowers(m, i)*e->fvPowers(n, j);
				mno2 = e->uPowers(m, i)*e->fpvPowers(n, j);
				/*double mno1_u = e->uPowers(m, i)*e->fvPowers(n, j);
				double mno2_u = e->uPowers(m, i)*e->fpvPowers(n, j);*/
			}
			else if (iuvw == 2) {
				mno1 = e->fpuPowers(m, i)*e->vPowers(n, j);
				mno2 = e->fuPowers(m, i)*e->vPowers(n, j);
				/*double mno1_v = e->fpuPowers(m, i)*e->vPowers(n, j);
				double mno2_v = e->fuPowers(m, i)*e->vPowers(n, j);*/
			}
			else {
				mno1 = e->fuPowers(m, i)*e->fpvPowers(n, j);
				mno2 = e->fpuPowers(m, i)*e->fvPowers(n, j);
				//double mno1_w = e->fuPowers(m, i)*e->fpvPowers(n, j);
				//double mno2_w = e->fpuPowers(m, i)*e->fvPowers(n, j);
			}
			if (iuvwh == 1) {
				mno1_h = e->uPowers(m, ih)*e->fvPowers(n, jh);
				mno2_h = e->uPowers(m, ih)*e->fpvPowers(n, jh);
		/*		double mno1_un = e->uPowers(m, ih)*e->fvPowers(n, jh);
				double mno2_un = e->uPowers(m, ih)*e->fpvPowers(n, jh);*/
			}
			else if (iuvwh == 2) {
				mno1_h = e->fpuPowers(m, ih)*e->vPowers(n, jh);
				mno2_h = e->fuPowers(m, ih)*e->vPowers(n, jh);
				/*double mno1_vn = e->fpuPowers(m, ih)*e->vPowers(n, jh);
				double mno2_vn = e->fuPowers(m, ih)*e->vPowers(n, jh);*/
			}
			else {
				mno1_h = e->fuPowers(m, ih)*e->fpvPowers(n, jh);
				mno2_h = e->fpuPowers(m, ih)*e->fvPowers(n, jh);
				/*double mno1_wn = e->fuPowers(m, ih)*e->fpvPowers(n, jh);
				double mno2_wn = e->fpuPowers(m, ih)*e->fvPowers(n, jh);*/
			}
			//

            for (int l = 1; l <= nglw; ++l) {
                for (int coord = 1; coord <= 3; ++coord) { //fix this section
                   
                    aw[coord] = e->awMatrix(m, n, l, coord);
                }
				for (int coord = 1; coord <= 3; ++coord) { //fix this section
					au[coord] = e->auMatrix(m, n, l, coord);
					
				}
				for (int coord = 1; coord <= 3; ++coord) { //fix this section
					
					av[coord] = e->avMatrix(m, n, l, coord);
					
				}
				double mno1_l;
				double mno2_l;
				double mno1_lh;
				double mno2_lh;
				if (iuvw == 1) {
					mno1_l = mno1*e->fpwPowers(l, k);
					mno2_l = mno2*e->fwPowers(l, k);
	
				}
				else if (iuvw == 2) {
					mno1_l=	mno1*e->fwPowers(l, k);
					mno2_l = mno2*e->fpwPowers(l, k);
				
				}
				else {
					mno1_l = mno1*e->wPowers(l, k);
					mno2_l = mno2*e->wPowers(l, k);
				
				}
				if (iuvwh == 1) {
					mno1_lh = mno1_h*e->fpwPowers(l, kh);
					mno2_lh = mno2_h*e->fwPowers(l, kh);
				
				}
				else if (iuvwh == 2) {
					mno1_lh = mno1_h*e->fwPowers(l, kh);
					mno2_lh = mno2_h*e->fpwPowers(l, kh);

				}
				else {
					mno1_lh = mno1_h*e->wPowers(l, kh);
					mno2_lh = mno2_h*e->wPowers(l, kh);

				}
                
				ROTFSEC(iuvwh, mno1_lh, mno2_lh , au, av, aw, rfus1, eval);
				ROTFSEC(iuvw, mno1_l, mno2_l, au, av, aw, rfus2, eval);
                if (e->materials.icode == 0) {
                    if (iHomCode == 0) {
                       
					   functions::MULTIPLYMURMATRIX(size1, e->MuRelIntInv, rfus2, m, n, l, f1); //3.4% fix!!!!!!!!!!!
                    }
                    else {
                    
						functions::MULTIPLYMURMATRIX(size1, e->MuRelIntInv, rfus2, 1, 1, 1, f1);
                    }
                }
                else {
                    if (iHomCode == 0) {
                        std::vector<std::complex<double>> temp;
                        for(auto val : rfus2){
                            temp.push_back(val/e->MuRelInt(1, m, n, l));
                        }
                        f1 = temp;
                    }
                    else {
						f1[1] = rfus2[1];
						f1[2] = rfus2[2];
						f1[3] = rfus2[3];
                    }
                }

                products::dotc(rfus1, f1, f2); //2.04% fix!!!!!!!
             
				double temp = e->wglw[l] / e->jacobian(m, n, l);
				cSl.real(cSl.real() + temp*f2.real());
				cSl.imag(cSl.imag() + temp*f2.imag());
            }
            cSn = cSn + e->wglv[n] * cSl;
        }
        cIntC = cIntC + e->wglu[m] * cSn;
    }
}




