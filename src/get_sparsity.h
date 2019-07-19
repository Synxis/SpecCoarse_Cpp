#ifndef GET_SPARSITY_H
#define GET_SPARSITY_H

#include <Eigen/Sparse>
#include <Eigen/Core>

// Inputs:
//   M  eigen sparse matrix
// Outputs: 
//   SM  the sparsity of matrix M (dtype double)
void get_sparsity(
	const Eigen::SparseMatrix<double> & M,
	Eigen::SparseMatrix<double> & SM);
#endif