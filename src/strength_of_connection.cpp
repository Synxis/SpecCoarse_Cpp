#include "strength_of_connection.h"

void strength_of_connection(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::SparseMatrix<double> & M,
	Eigen::SparseMatrix<double> & SoC)
{
	using namespace Eigen;
	using namespace std;

	VectorXd MVec = M.diagonal();
	SoC.resize(L.rows(), L.cols());

	// remove diagonal from L
	SoC = L;
	SoC.prune([](const int r, const int c, const double)->bool{return r!=c;});

	// construct the strength matrix
	for(int ii = 0; ii<SoC.outerSize(); ++ii)
	{
		for(typename SparseMatrix<double>::InnerIterator it (SoC,ii); it; ++it)
		{
			assert(it.value() != 0);
			int i = it.col();
			int j = it.row();
			it.valueRef() = -sqrt(MVec(i) + MVec(j)) / it.valueRef();
			if (it.valueRef() < 0)
				it.valueRef() = 0;
		}
	}		
}
