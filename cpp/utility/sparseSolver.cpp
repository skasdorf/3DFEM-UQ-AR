
#define EIGEN_USE_MKL_ALL
#include "sparseSolver.h"

//#include "../eigen-eigen-5a0156e40feb/Eigen/Sparse.h"
//#include <Eigen/PardisoSupport>

void sparseSolver::doEverything(std::vector<std::complex<double>>& cFEMrSparse, std::vector<int>& iRow, std::vector<int>& jCol, 
	int& nNZ, int& size, const std::vector<std::complex<double>>& cGr, SpMat& sparseMatrix, std::vector<std::complex<double>>& result){
	//sets up the sparseSolver object
	//this->factorMatrix(makeMatrixEigenFormat(cFEMrSparse, iRow, jCol, nNZ, size));
	//SpMat sparseMatrix(size, size);
	sparseMatrix.reserve(nNZ);
	std::vector<Triplet> tripletList(nNZ);
	for (int i = 1; i <= nNZ; ++i) {
		tripletList.push_back(Triplet(iRow[i] - 1, jCol[i] - 1, cFEMrSparse[i]));
	}
	sparseMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	sparseMatrix.makeCompressed();
	//error function
	Eigen::VectorXcd newRHS(size);
	for (int i = 0; i < size; ++i) {
		newRHS[i] = cGr[i + 1];
	}
	//for pardiso
	Solver internalSolver;
	internalSolver.analyzePattern(sparseMatrix);
	internalSolver.factorize(sparseMatrix);
	//end for pardiso
	//for sparseQr
	//std::cout << "Doing solving..." << std::endl;
	////Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int>> internalSolver;
	//Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> internalSolver;
	//internalSolver.analyzePattern(sparseMatrix);
	//internalSolver.factorize(sparseMatrix);
	//end for sparseQr
	Eigen::VectorXcd x = internalSolver.solve(newRHS);
	for (int i = 0; i < size; ++i) {
		result[i + 1] = x[i];
	}

	//
	//FOR ITERATIVE SOLVER!
	//
	/*std::cout << "Solving matrix!" << std::endl;
	std::cout << "Num threads being used: " << Eigen::nbThreads();
	ItSolver internalSolver;
	internalSolver.compute(sparseMatrix);
	if (internalSolver.info() != Eigen::Success)
	{
		std::cout << internalSolver.info();
		throw std::exception("Error computing matrix");
	}
	Eigen::VectorXcd x = internalSolver.solve(newRHS);
	std::cout << "Iterative solver error (estimate): " << internalSolver.error() << std::endl;
	std::cout << "Iterative solver iterations: " << internalSolver.iterations() << std::endl;
	for (int i = 0; i < size; ++i) {
		result[i + 1] = x[i];
	}*/
}

