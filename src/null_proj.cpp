#include "null_proj.h"

void null_proj(
	const Eigen::SparseMatrix<double> & G,
    const Eigen::VectorXd & nullVec,
    const Eigen::SparseMatrix<double> & d,
    const Eigen::VectorXd & invd2,
	Eigen::SparseMatrix<double> & G_proj)
{
    using namespace Eigen;

    // compute d*G'
    SparseMatrix<double> dGMat;
    dGMat = d.cwiseProduct(G); 

    VectorXd dG;
    igl::sum(dGMat, 2, dG);

    // compute step size along the direction d
    VectorXd projStepSizeVec;
    projStepSizeVec = (nullVec - dG).array() * invd2.array();
    
    // compute stepSize * d
    SparseMatrix<double> step_d;
    step_d.resize(d.rows(), d.cols());
    step_d = d;
    for (int k=0; k<step_d.outerSize(); ++k) 
        for (SparseMatrix<double>::InnerIterator it(step_d,k); it; ++it)
            it.valueRef() *= projStepSizeVec(it.row());
    
    // compute progG = G + stepSize * d
    G_proj.resize(G.rows(), G.cols());
    G_proj = G + step_d;
}