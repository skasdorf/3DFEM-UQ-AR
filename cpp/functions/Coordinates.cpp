//NEEDS COMMENTS FOR EXPLANATION
#include "../functions/Coordinates.h"


namespace coordinates {

	double dSinD2(double angle) {
		double sinAngle;
		if (angle == 180.0) {
			sinAngle = 0.0;
		}
		else {
			sinAngle = std::sin(angle*constants::PI / 180.0);
		}
		return sinAngle;
	} //double dSinD2

	double dCosD2(double angle) {
		double cosAngle;
		if (angle == 90.0) {
			cosAngle = 0.;
		}
		else if (angle == 270.0) {
			cosAngle = 0.0;
		}
		else {
			cosAngle = std::cos(angle*constants::PI / 180.0);
		}
		return cosAngle;
	} //double dCosD2

	std::vector<double> findiPhi(double phi) {
		double cosPhi = dCosD2(phi);
		double sinPhi = dSinD2(phi);
		std::vector<double> iPhi(4);
		iPhi[1] = -sinPhi;
		iPhi[2] = cosPhi;
		iPhi[3] = 0.0;
		return iPhi;
	} //std::vector<double> findiPhi

	std::vector<double> findiTheta(double theta, double phi) {
		double cosTheta = dCosD2(theta);
		double sinTheta = dSinD2(theta);
		double cosPhi = dCosD2(phi);
		double sinPhi = dSinD2(phi);
		std::vector<double> iTheta(4);
		iTheta[1] = cosTheta*cosPhi;
		iTheta[2] = cosTheta*sinPhi;
		iTheta[3] = -sinTheta;
		return iTheta;
	} //std::vector<double> findiTheta

	std::vector<double> findnAr(double theta, double phi) {
		double cosTheta = dCosD2(theta);
		double sinTheta = dSinD2(theta);
		double cosPhi = dCosD2(phi);
		double sinPhi = dSinD2(phi);
		std::vector<double> nAr(4);
		nAr[1] = -sinTheta*cosPhi;
		nAr[2] = -sinTheta*sinPhi;
		nAr[3] = -cosTheta;
		return nAr;
	} //std::vector<double> findnAr

} //namespace coordinates

