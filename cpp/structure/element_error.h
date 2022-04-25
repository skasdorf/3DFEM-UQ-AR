#pragma once
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include "../utility/functions.h"
#include <Eigen\Core>
#include "../structure/Domain.h"
#include "../functions/unitVectorsM.h"

#include "../structure/Point.h"
namespace qoi_error {
	/*bool is_higher(int& u_order, int& v_order, int& w_order,
		int& iuvwh, int& ih, int& jh, int& kh);*/
	std::complex<double> check_q(std::string cAlphaAdjointFile, std::string cGrFile);

	
}