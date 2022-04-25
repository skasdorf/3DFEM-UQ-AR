

//#include "findIntegrals.h"
//
//
//findIntegrals::findIntegrals(Domain domain) {
//    dom = domain;
//}
//
//void findIntegrals::find_integrals(int eCNT, int iuvwh, int iuvw, int ih, int jh, int kh, int i, int j, int k,
//                                   std::complex<double>& cSolA, std::complex<double>& cSolB, std::complex<double>& cSolS,
//                                   std::complex<double>& cSolStotal, int size1, matrix4d<std::complex<double>>& MuRelInt,
//                                   matrix4d<std::complex<double>>& MuRelIntInv, matrix4d<std::complex<double>>& EpsRelInt){
//
//    Element element = dom.elements[eCNT];
//
//    Integral_c int_c = Integral_c();
//    int_c.findC(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolA, size1, element, MuRelInt, MuRelIntInv);
//
//    Integral_d int_d = Integral_d();
//    int_d.findD(iuvwh, iuvw, ih, jh, kh, i, j, k, cSolB, size1, element, EpsRelInt);
//}