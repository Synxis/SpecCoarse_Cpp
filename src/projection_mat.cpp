#include "projection_mat.h"

void projection_mat(
	const Eigen::SparseMatrix<double> & L,
	const Eigen::VectorXi & clusterCenters,
	Eigen::SparseMatrix<double> & ProjMat)
{
	using namespace Eigen;
	std::vector<Triplet<double>> IJV;
	IJV.reserve(clusterCenters.size());
	for(int ii = 0; ii < clusterCenters.size(); ii++)
		IJV.push_back(Triplet<double>(ii, clusterCenters(ii), 1.0));
	ProjMat.resize(clusterCenters.size(), L.rows());
	ProjMat.setFromTriplets(IJV.begin(),IJV.end());
}