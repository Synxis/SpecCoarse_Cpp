#ifndef PRINT_SPMAT_H
#define PRINT_SPMAT_H

#include <Eigen/Sparse>
#include <Eigen/Core>
// Inputs:
//   mat eigen sparse matrix
void print_spmat(
	const Eigen::SparseMatrix<double> & mat)
{
	using namespace Eigen;
	using namespace std;

	int numPrint = 0;
	int numToPrint = 8;
	for (int k=0; k<mat.outerSize(); ++k)
	{
		for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
		{
			if (numPrint < numToPrint)
				cout << "(" << it.row() << "," << it.col() << "): " << it.value() << endl;
			else if (numPrint > mat.nonZeros() - numToPrint + 1)
				cout << "(" << it.row() << "," << it.col() << "): " << it.value() << endl;
			else if (numPrint == numToPrint)
				cout << "..." << endl;
			numPrint++;
		}
	}
	cout << "matrix size: (" << mat.rows() << "," << mat.cols() << ")" << endl;
	cout << "number of nonzeros: " << mat.nonZeros() << endl;
}
#endif