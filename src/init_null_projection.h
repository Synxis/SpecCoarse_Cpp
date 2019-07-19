#ifndef INIT_NULL_PROJECTION_H
#define INIT_NULL_PROJECTION_H
#include <igl/sum.h>
#include <igl/find.h>

// Inputs:
//   S  sparsity matrix
//   P  projection matrix
//   nullVec null space vector (dense domain)
// Outputs: 
//   d null space vector in sparsity
//   invd2  row sum inverse ||d||^2 for projection 
void init_null_projection(
	const Eigen::SparseMatrix<double> & S,
	const Eigen::SparseMatrix<double> & P,
	const Eigen::VectorXd & nullVec,
	Eigen::SparseMatrix<double> & d,
	Eigen::VectorXd & invd2);
#endif