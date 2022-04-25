//
//
//#include "sparsePacking.h"
//
//using dcomplex = std::complex<double>;
//
//namespace functions {
//
//    void SPARSE_PACKING(Domain& dom, std::vector<int>& iRow,
//                        std::vector<int>& jCol, int& nNZ, int& noAbcBCs, std::vector<dcomplex>& cArSparse, 
//		std::vector<dcomplex>& cBrSparse, std::vector<dcomplex>& cSrSparse) {
//		//noAbcBCs should be 0
//
//		//std::vector<Element>* elements = &(dom.elements);
//        int myElementCount = dom.elements.size();
//        std::vector<int> IsContained = std::vector<int>(myElementCount+1);
//
//        //*
////        *
////        *
////        *
////        *
////        *
////        *
////        *
////        *
//        //Need to assign these values to what they actuall should be
//     //   std::vector<dcomplex> cArSparse, cBrSparse, cSrSparse;
////        *
////        *
////        *
////        *
//
//        int e1, e2, i0, j0, i1, j1, lim1, D1, D2, C1, C2;
//        int index1, index2, k0, k1, k2, k3;
//        int FoundFlag;
//
//      //  for (int eCNT = 0; eCNT < myElementCount; ++eCNT) {
//        //    IsContained[e->index] = 1;
//       // }
//
//        for (auto elem = ++dom.elements.begin(); elem != dom.elements.end(); ++elem) {
//
//            //int eCNT = myElementsGlobal[elem];
//
//           // lim1 = dom.elements[eCNT].unknownsStart - 1;
//			lim1 = elem->unknownsStart - 1;
//            for (int i = elem->unknownsStart; i <= elem->unknownsEnd; ++i) {
//                i0 = i - lim1;
//                D1 = dom.vectorDUnique[i];
//                //C1 = The index of where the basis unknown goes globally
//                C1 = std::abs(dom.vectorD[i]);
//                //Skip if this unknown is zero
//                if (D1 == 0) {
//                    continue;
//                }
//
//                for (int j =elem->unknownsStart; j <= elem->unknownsEnd; ++j) {
//                    //How far from the beginning element we are for the testing function
//                    j0 = j - lim1;
//                    D2 = dom.vectorDUnique[j];
//                    //C2 = the index of where the testing unknown goes globally
//                    C2 = std::abs(dom.vectorD[j]);
//                    if (D2 != 0) {
//                        e2 = 0;
//
//                        //Index 1 and Index 2 tell us if the unknowns appear in more than one element
//                        index1 = dom.unknown_description[C1][1];
//                        index2 = dom.unknown_description[C2][1];
//
//                        //If both unknowns (basis and testing) appear in the global system at some point
//                        if ((index1 >= 1) && (index2 >= 1)) {
//
//                            //Whether or not we have found this entry on a different element
//                            FoundFlag = 0;
//                            //For all elements in which either unknown appears
//                            for (int k0 = 1; k0 <= index1; ++k0) {
//                                for (int k1 = 1; k1 <= index2; ++k1) {
//                                    //If both unknowns appear in a previous element, which is owned by this process,
//                                    //  those entries should be copied over to the previous element
//                                    //Note: If it corresponds to an element owned by another process, the entries will
//                                    //   automatically be added together by MUMPS (parallel) solver
//                                    if ((dom.unknown_description[C1][ 2 * k0] == dom.unknown_description[C2][ 2 * k1]) &&
//                                        (dom.unknown_description[C1][ 2 * k0] < elem->index) && (FoundFlag == 0)) //&&
//                                        /*(IsContained[unknown_description(C1, 2 * k0)] == 1))*/ {
//                                        //The element it should be copied to
//                                        e2 = dom.unknown_description[C1][ 2 * k0];
//                                        //The global disconnected unknown on that element (basis)
//                                        k2 = k0;
//                                        //The global disconnected unknown on that element (testing)
//                                        k3 = k1;
//                                        FoundFlag = 1;
//                                    }
//                                }
//                            }
//                        }
//
//
//                        //If we found somewhere else to put those entries, copy them over to the previous element
//                        if (e2 != 0) {
//
//                            i1 = dom.unknown_description[C1][ 2 * k2 + 1] - dom.elements[e2].unknownsStart + 1;
//                            j1 = dom.unknown_description[C2][ 2 * k3 + 1] - dom.elements[e2].unknownsStart + 1;
//
//                            dom.elements[e2].cArPAK(i1, j1) =
//                                    dom.elements[e2].cArPAK(i1, j1) + elem->cArPAK(i0, j0);
//                            dom.elements[e2].cBrPAK(i1, j1) =
//                                    dom.elements[e2].cBrPAK(i1, j1) + elem->cBrPAK(i0, j0);
//                            if (noAbcBCs > 0) {
//                                dom.elements[e2].cSrPAK(i1, j1) =
//                                        dom.elements[e2].cSrPAK(i1, j1) + elem->cSrPAK(i0, j0);
//                            }
//
//
//                            elem->cArPAK(i0, j0) = {0, 0};
//                            elem->cBrPAK(i0, j0) = {0, 0};
//                            if (noAbcBCs > 0) {
//                                elem->cSrPAK(i0, j0) = {0, 0};
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//
//        //        if(DEBUG .AND. id .EQ. 0) Then
//        //        Write(*,*) 'Finished Combining Elements, beginning Sparse Packing'
//        //FLUSH(6)
//        //End If
//
//
//        //for (int elem = 1; elem <= myElementCount; ++elem) {
//		for (auto elem = dom.elements.begin(); elem != dom.elements.end(); ++elem){
//          //  int eCNT = myElementsGlobal[elem];
//
//
//            //!WRITE(*,*) 'Packing Element No.',eCNT
//            //!FLUSH(6)
//
//            lim1 = elem->unknownsStart - 1;
//
//            for (int i = elem->unknownsStart; i <= elem->unknownsEnd; ++i) {
//                i0 = i - lim1;
//                D1 = abs(dom.vectorD[i]);
//                if (D1 == 0) { continue; }
//
//                for (int j = i; j <= elem->unknownsEnd; ++j) {
//                    j0 = j - lim1;
//                    D2 = abs(dom.vectorD[j]);
//                    if (D2 != 0) {
//
//                        dcomplex czero = {0,0};
//                        if ((elem->cArPAK(i0, j0) != czero) ||
//                            (elem->cBrPAK(i0, j0) != czero)) {
//
//                            if (D1 <= D2) {
//
//                                nNZ = nNZ + 1;
//                                cArSparse[nNZ] = elem->cArPAK(i0, j0);
//                                cBrSparse[nNZ] = elem->cBrPAK(i0, j0);
//                                if (noAbcBCs > 0) { cSrSparse[nNZ] = elem->cSrPAK(i0, j0); }
//
//                                iRow[nNZ] = D1;
//                                jCol[nNZ] = D2;
//
//                                if (D1 != D2) {
//                                    nNZ = nNZ + 1;
//                                    cArSparse[nNZ] = elem->cArPAK(j0, i0);
//                                    cBrSparse[nNZ] = elem->cBrPAK(j0, i0);
//                                    if (noAbcBCs > 0) { cSrSparse[nNZ] = elem->cSrPAK(j0, i0); }
//
//                                    iRow[nNZ] = D2;
//                                    jCol[nNZ] = D1;
//                                } else {
//
//                                    nNZ = nNZ + 1;
//                                    cArSparse[nNZ] = elem->cArPAK(i0, j0);
//                                    cBrSparse[nNZ] = elem->cBrPAK(i0, j0);
//                                    if (noAbcBCs > 0) { cSrSparse[nNZ] = elem->cSrPAK(i0, j0); }
//
//                                    iRow[nNZ] = D2;
//                                    jCol[nNZ] = D1;
//
//                                    if (D1 != D2) {
//                                        nNZ = nNZ + 1;
//                                        cArSparse[nNZ] = elem->cArPAK(j0, i0);
//                                        cBrSparse[nNZ] = elem->cBrPAK(j0, i0);
//                                        if (noAbcBCs > 0) cSrSparse[nNZ] = elem->cSrPAK(j0, i0);
//
//                                        iRow[nNZ] = D1;
//                                        jCol[nNZ] = D2;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        //    IF(DEBUG.AND.id.EQ.
//        //    0) Then
//        //    Write(*, *)
//        //    'Finished Sparse Packing'
//        //    End If
//
//    }
//}
