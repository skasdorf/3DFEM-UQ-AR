//NEEDS COMMENTS FOR EXPLANATION
#pragma once
#ifndef CPP_COORDINATES_H
#define CPP_COORDINATES_H

#include <fstream>
#include <sstream>
#include "../utility/constants.h"
//#include "../functions/Coordinates.h"


namespace coordinates {
	double dSinD2(double angle);
	double dCosD2(double angle);
	std::vector<double> findiPhi(double phi);
	std::vector<double> findiTheta(double theta, double phi);
	std::vector<double> findnAr(double theta, double phi);
};


#endif //CPP_COORDINATES_H