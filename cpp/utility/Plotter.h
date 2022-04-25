//plotter class
#pragma once
#include "../matricies/matrix4d.h"
#include "../matricies/matrix2d.h"
#include "../matricies/matrix3d.h"
#include "../structure/Domain.h"
#include "../functions/unitVectorsM.h"
#include "../functions/products.h"
class PlotController {
public:
	int numberOfPlots;
	std::vector<int> element_num_for_plot;
	std::vector<double> u1, u2;
	std::vector<int> Nu;
	std::vector<double> v1, v2;
	std::vector<int> Nv;
	std::vector<double> w1, w2;
	std::vector<int> Nw;
	Domain* dom;

	PlotController(Domain* dom) {
		this->dom = dom;
	}
	void FUODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
		const int& k, const matrix2d<double>& uPowers, const matrix2d<double>& fvPowers, const matrix2d<double>& fwPowers,
		const int& nPu, const int& nPv, const int& nPw, double& fu);
	
void FVODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
	const int& k, const matrix2d<double>& fuPowers, const matrix2d<double>& vPowers, const matrix2d<double>& fwPowers,
	const int& nPu, const int& nPv, const int& nPw, double& fv);
void FWODUVWPLOT(const int& m, const int& n, const int& l, const int& i, const int& j,
	const int& k, const matrix2d<double>& fuPowers, const matrix2d<double>& fvPowers, const matrix2d<double>& wPowers,
	const int& nPu, const int& nPv, const int& nPw, double& fw);
void find_ausec(std::vector<double>& av, std::vector<double>& aw, std::vector<double>& ausec);
void find_avsec(std::vector<double>& aw, std::vector<double>& au, std::vector<double>& avsec);
void find_awsec(std::vector<double>& au, std::vector<double>& av, std::vector<double>& awsec);
void findPowers(const int& nu, const int& nv, const int& nw, const int& nglu, const int& nglv, const int& nglw, const std::vector<double>& xglu, const std::vector<double>& xglv, const std::vector<double>& xglw,
	matrix2d<double>& uPowers, matrix2d<double>& vPowers, matrix2d<double>& wPowers, matrix2d<double>& fuPowers, matrix2d<double>& fvPowers, matrix2d<double>& fwPowers,
	matrix2d<double>& fpuPowers, matrix2d<double>& fpvPowers, matrix2d<double>& fpwPowers);
void unitaryvectorsplot(int& nPu, int& nPv, int& nPw, std::vector<double>& xglu, std::vector<double>& xglv, std::vector<double>& xglw, int& elementOrders,
	matrix2d<double>& rs, matrix4d<double>& rMatrix, matrix4d<double>& auMatrix, matrix4d<double>& avMatrix, matrix4d<double>& awMatrix, matrix3d<double>& jacobianMatrix, const int& nRs);
void plotField(std::string mesh_name, bool plotAdjoint);

void plotterIO(std::string plot_file, std::string coeff_file, bool use_file);

};
