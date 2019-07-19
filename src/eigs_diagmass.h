#ifndef EIGS_DIAGMASS_H
#define EIGS_DIAGMASS_H

#include <MatOp/SparseSymShiftSolve.h>
#include <SymEigsShiftSolver.h>
#include <iostream>

// Inputs:
//   L        system matrix 
//   diagMass mass matrix (only for diagonal mass matrix)
// 	 numEigs  number of eigenvalues/vectors to compute
// Outputs: 
//   eVal     eigenvalue from small to large (numEigs)
//   eVec     eigenvectors (n-by-numEigs)
void eigs_diagmass(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::SparseMatrix<double> & diagMass,
	const int numEigs,
	Eigen::VectorXd & eVal,
	Eigen::MatrixXd & eVec);
#endif