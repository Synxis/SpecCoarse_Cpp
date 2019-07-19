#ifndef STRENGTH_OF_CONNECTION_H
#define STRENGTH_OF_CONNECTION_H
#include <igl/find.h>
#include <igl/diag.h>

// Inputs:
//   L  SPSD system matrix (off diagonals are mostly negative)
//   M  diagonal mass matrix
// Outputs: 
//   SoC  symmetric strength of connection matrix
void strength_of_connection(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::SparseMatrix<double> & M,
	Eigen::SparseMatrix<double> & SoC);
#endif