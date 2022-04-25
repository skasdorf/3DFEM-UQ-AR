
#pragma once
#ifndef CPP_FUNCTIONS_H
#define CPP_FUNCTIONS_H

#include <fstream>
#include <sstream>
#include "constants.h"

//#include "../functions/basis.h"


namespace functions {

    std::vector<std::string> split(const std::string &s, char delim);
    void MIMAISWAP(int *, int *);
	bool _debug_matrix2(std::string filename, std::vector<int>&, int);
	bool _debug_matrix(std::string filename, std::vector<int>&, int);
	bool _debug_matrix(std::string filename, std::vector<std::vector<int>>&, int);
	bool _debug_matrix(std::string filename, std::vector<std::vector<std::vector<int>>>&, int);
	bool _debug_matrix(std::string filename, std::vector<double>&, int);
	bool _debug_matrix(std::string filename, std::vector<std::vector<double>>&, int);
	bool _debug_matrix(std::string filename, std::vector<std::vector<std::vector<double>>>&, int);

	void gaussk(int, std::vector<double>&, std::vector<double>&);
	double trilagrangeoduvw(const int& m, const int& n, int& l, const int& i, const int& j, const int& k, const matrix2d<double>& fuPowersLagr, const matrix2d<double>& fvPowersLagr, 
		const matrix2d<double>& fwPowersLagr, const int& nglu, const int& nglv, const int& nglw);
	void find_mur_inv(const int& size1, const matrix4d<std::complex<double>>& MuRel, const int& nglu, const int& nglv, const int& nglw, matrix4d<std::complex<double>>& MuRelInv);

	void read_previous_solve(std::string mesh_name, std::vector<std::complex<double>>& cAlphaForLower, std::vector<std::complex<double>>& cAlphaAdjoint, std::vector<std::complex<double>>& cGrhigher);
		
};


#endif //CPP_FUNCTIONS_H