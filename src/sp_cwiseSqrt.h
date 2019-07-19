#ifndef SP_CWISESQRT_H
#define SP_CWISESQRT_H

#include <Eigen/Sparse>

void sp_cwiseSqrt(
	const Eigen::SparseMatrix<double> & spM,
	Eigen::SparseMatrix<double> & spMCwiseSqrt);
#endif