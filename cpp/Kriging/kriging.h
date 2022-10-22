#include "../structure/Domain.h";
#include "../run_modes.h"

namespace Kriging {
	std::complex<long double> sensitivity_to_epsr(std::string& file_name, std::vector<std::complex<long double>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::complex<long double> referenceFreq, std::complex<long double>& gradient);
	std::complex<long double> get_QoI(Domain& dom_forward, Domain& dom_adjoint);
	std::complex<long double> get_QoI2(Domain& dom_forward, Domain& dom_adjoint);
	void monte_carlo_instance(std::string& file_name);
	void multi_HOPS_epsr(std::string& file_name);
	void multi_HOPS_multi_epsr(std::string& file_name);
	long double normalFunction(long double x, long double mean, long double std);
	long double trap_integral(long double intervalMid, long double mean, long double std, long double step);
	void sensitivity_to_multi_epsr(std::string& file_name, std::vector<std::vector<std::complex<long double>>>& epsr_list, std::vector<std::complex<long double>>& updated_qoi, std::vector<std::complex<long double>>& reference_values);
}