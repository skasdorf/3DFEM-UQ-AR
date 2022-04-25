#ifndef BASISEVAL_H
#define BASISEVAL_H

#include "../matricies/matrix2d.h"
#include"../utility/functions.h"
#include <string>
#include <unordered_map>

using evaluation = std::vector<matrix2d<double>>;
using evalMap = std::unordered_map<std::string, evaluation>;

class BasisEval {
public: 
	//BasisEval(int indicatorBasisType) {
	//	this->basisType = indicatorBasisType; //basisType tells us what kind of basis functions we are using
	//}
	int basisType = 1; //LEGENDRE NUMBER ONE
	evalMap uMap;
	evalMap vMap;
	evalMap wMap;

	evaluation& get_u_samples(int nglu, int nu); 
	evaluation& get_v_samples(int nglu, int nu);
	evaluation& get_w_samples(int nglu, int nu);

	matrix2d<double> leg, leg2;
	void matrix_modifiy(matrix2d<double>& mat, int start_x, int end_x, int start_y, int end_y, double& scalar);
	void finduPowers(const int& nu, const int& nglu, const std::vector<double>& xglu, matrix2d<double>& uPowers, matrix2d<double>& fuPowers, matrix2d<double>& fpuPowers);
	void findvPowers(const int& nv, const int& nglv, const std::vector<double>& xglv, matrix2d<double>& vPowers, matrix2d<double>& fvPowers, matrix2d<double>& fpvPowers);
	void findwPowers(const int& nw, const int& nglw, const std::vector<double>& xglw, matrix2d<double>& wPowers, matrix2d<double>& fwPowers, matrix2d<double>& fpwPowers);

	void finduPowersLeg(const int& nu, const int& nglu, const std::vector<double>& xglu, matrix2d<double>& uPowers, matrix2d<double>& fuPowers, matrix2d<double>& fpuPowers);
	void findvPowersLeg(const int& nv, const int& nglv, const std::vector<double>& xglv, matrix2d<double>& vPowers, matrix2d<double>& fvPowers, matrix2d<double>& fpvPowers);
	void findwPowersLeg(const int& nw, const int& nglw, const std::vector<double>& xglw, matrix2d<double>& wPowers, matrix2d<double>& fwPowers, matrix2d<double>& fpwPowers);

	void setup_legendre();
};


#endif