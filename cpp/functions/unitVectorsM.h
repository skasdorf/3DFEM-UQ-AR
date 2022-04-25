#pragma once
//unitvectorsmeda
#ifndef UNITVECTORSM_H
#define UNITVECTORSM_H

//#include "../utility/constants.h"
#include "../functions/products.h"
namespace unitVectorsM {
	void unitaryvectors(const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv,
		const std::vector<double>& xglw, const int& elementOrders, const matrix2d<double>& rs, matrix4d<double>& rMatrix,
		matrix4d<double>& auMatrix, matrix4d<double>& avMatrix, matrix4d<double>& awMatrix, matrix3d<double>& jacobianMatrix, const int& nRs, const int& integralType);
	
	std::vector<double> findR(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);
	std::vector<double> findRA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);

	std::vector<double> findAU(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);
	std::vector<double> findAUA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);

	std::vector<double> findAV(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);
	std::vector<double> findAVA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);

	std::vector<double> findAW(const int& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);
	std::vector<double> findAWA(const std::vector<int>& elementOrders, const double& u, const double& v, const double& w, const matrix2d<double>& rs, const int& nRs);
	void find_jacobian(const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw, double& jacobian);
	void find_jacobian_surface(const std::vector<double>& au, const std::vector<double>& av, const std::vector<double>& aw, const int& face, double& result);

	
}
#endif //UNITVECTORSM_H
