#include "sp_cwiseSqrt.h"

void sp_cwiseSqrt(
	const Eigen::SparseMatrix<double> & spM,
	Eigen::SparseMatrix<double> & spMCwiseSqrt)
{
	using namespace Eigen;
	spMCwiseSqrt.resize(spM.rows(), spM.cols());
	spMCwiseSqrt = spM;
	for (int k=0; k<spMCwiseSqrt.outerSize(); ++k) 
		for (SparseMatrix<double>::InnerIterator it(spMCwiseSqrt,k); it; ++it){
			it.valueRef() = sqrt(it.valueRef());
		}
}