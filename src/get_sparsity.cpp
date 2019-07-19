#include "get_sparsity.h"

void get_sparsity(
	const Eigen::SparseMatrix<double> & M,
	Eigen::SparseMatrix<double> & SM)
{
    using namespace Eigen;
    std::vector<Triplet<double>> IJV;
	IJV.reserve(M.nonZeros());
    for(int ii = 0; ii<M.outerSize(); ++ii)
	{
		for(typename SparseMatrix<double>::InnerIterator it (M,ii); it; ++it)
		{
			assert(it.value() != 0);
            IJV.push_back(Triplet<double>(it.row(), it.col(), 1.0));
		}
	}	
	SM.resize(M.rows(), M.cols());
	SM.setFromTriplets(IJV.begin(),IJV.end());	  
}