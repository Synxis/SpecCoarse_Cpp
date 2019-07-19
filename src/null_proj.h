#ifndef NULL_PROJ_H
#define NULL_PROJ_H

#include <Eigen/Sparse>
#include <igl/sum.h>

void null_proj(
	const Eigen::SparseMatrix<double> & G,
    const Eigen::VectorXd & nullVec,
    const Eigen::SparseMatrix<double> & d,
    const Eigen::VectorXd & invd2,
	Eigen::SparseMatrix<double> & G_proj);
#endif