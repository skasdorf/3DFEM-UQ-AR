#include "../structure/Domain.h";
#include "../run_modes.h"

namespace Kriging {
	std::complex<double> sensitivity_to_epsr(std::string& file_name, std::vector<std::complex<double>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::complex<double> referenceFreq, std::complex<double>& gradient);
	std::complex<double> get_QoI(Domain& dom_forward, Domain& dom_adjoint);
	std::complex<double> get_QoI2(Domain& dom_forward, Domain& dom_adjoint);
	void monte_carlo_instance(std::string& file_name);
	void multi_HOPS_epsr(std::string& file_name);
	void multi_HOPS_multi_epsr(std::string& file_name);
	double normalFunction(double x, double mean, double std);
	//std::complex<double> taylorKrigingHigherOrder3(std::vector<double> variogramFit, std::vector<double> variogramFitFirstD, std::vector<double> variogramFitSecondD, std::vector<double> d, std::vector<std::complex<double>> xSample, std::vector<double>& weights, double xVal, std::vector<double> variogramFitGrad, std::vector<std::complex<double>> yVals, std::vector<std::complex<double>> yValsGrad);

	double trap_integral(double intervalMid, double mean, double std, double step);
	void sensitivity_to_multi_epsr(std::string& file_name, std::vector<std::vector<std::complex<double>>>& epsr_list, std::vector<std::complex<double>>& updated_qoi, std::vector<std::complex<double>>& reference_values);
}