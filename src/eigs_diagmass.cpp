#include "eigs_diagmass.h"

void eigs_diagmass(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::SparseMatrix<double> & diagMass,
	const int numEigs,
	Eigen::VectorXd & eVal,
	Eigen::MatrixXd & eVec)
{
	using namespace Spectra;
	using namespace Eigen;

	// diagonal inverse of the sqrt of the mass matrix
  SparseMatrix<double> sqrtMinv;
  sqrtMinv.resize(diagMass.rows(), diagMass.cols());
  sqrtMinv = diagMass;
  for(int ii = 0; ii<sqrtMinv.outerSize(); ++ii)
  {
    for(typename SparseMatrix<double>::InnerIterator it (sqrtMinv,ii); it; ++it)
    {
      if(it.col() == it.row())
      {
        double v = it.value();
        assert(v != 0);
        v = (1.0)/sqrt(v);
        sqrtMinv.coeffRef(it.row(),it.col()) = v;
      }
    }
  }		

  // solve non-generalized eigen value problem
  double ncv = (3*numEigs > L.rows()) ? L.rows() : 3*numEigs;
  SparseMatrix<double> LHS(L.rows(), L.cols());
  LHS = sqrtMinv * L * sqrtMinv;
  SparseSymShiftSolve<double> op(LHS);
  SymEigsShiftSolver<double, LARGEST_MAGN, SparseSymShiftSolve<double>> 
  	eigs(&op, numEigs, ncv, -1e-6);
  eigs.init();
  eigs.compute();

  // get eigenvalue/vectors
	MatrixXd eVecTemp;
	eVal.resize(numEigs);
	if(eigs.info() == Spectra::SUCCESSFUL)
	{
		eVal = eigs.eigenvalues();
		eVecTemp = eigs.eigenvectors();
	}
	else
		std::cout << "eigs does not converge" << std::endl;

	// recover original eVec
	eVec.resize(eVecTemp.rows(), eVecTemp.cols());
	eVec = sqrtMinv * eVecTemp;

	// reorder the output
	eVal = eVal.reverse().eval();
	eVec = eVec.rowwise().reverse().eval();
}