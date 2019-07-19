#ifndef PROJECTION_MAT_H
#define PROJECTION_MAT_H

#include <Eigen/Sparse>
#include <Eigen/Core>
// Inputs:
//   L  system matrix
//   clusterCenters  cluster centers (corresponding to rows of Gt)
// Outputs: 
//   ProjMat  projection matrix from dense to coarse
void projection_mat(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::VectorXi & clusterCenters,
	Eigen::SparseMatrix<double> & ProjMat);
#endif