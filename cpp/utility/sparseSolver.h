#ifndef CPP_SOLVESYSTEM_H
#define CPP_SOLVESYSTEM_H
#include <Eigen/Core>
#include <Eigen/Sparse.h>
#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>

//#include "../eigen-eigen-5a0156e40feb/Eigen/Core.h"
//#include "../eigen-eigen-5a0156e40feb/Eigen/Sparse.h"
#include <vector>
typedef Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> SpMat;
typedef Eigen::PardisoLU<SpMat> Solver;
typedef Eigen::Triplet<std::complex<double>> Triplet;
typedef Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<std::complex<double>>> ItSolver;

class sparseSolver {
public:
	//Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> internalSolver;

	void doEverything(std::vector<std::complex<double>> &cFEMrSparse, std::vector<int>& iRow, std::vector<int>& jCol,
		int& nNZ, int& size, const std::vector<std::complex<double>>& cGr, SpMat& sparseMatrix, std::vector<std::complex<double>>& result);
	/*void setup(std::vector<std::complex<double>> &cFEMrSparse, std::vector<int>& iRow, std::vector<int>& jCol, int& nNZ, int& size);

	Eigen::SparseMatrix<std::complex<double>> makeMatrixEigenFormat(const std::vector<std::complex<double>>& cFEMrSparse, const std::vector<int>& iRow,
		const std::vector<int> & jCol,
		const int &nNZ, const int& size);

	void factorMatrix(Eigen::SparseMatrix<std::complex<double>>& stiffnessMatrix);
	
	Eigen::VectorXcd makeRHSEigenFormat(std::vector<std::complex<double>>& RHS);

	std::vector<std::complex<double>> solveSystem(std::vector<std::complex<double>>& RHS);*/

};

#endif //CPP_SOLVESYSTEM_H