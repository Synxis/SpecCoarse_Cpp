#include "combinatorial_coarsening.h"

void combinatorial_coarsening(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::SparseMatrix<double> & K,
	Eigen::VectorXi & clusterCenters)
{
	using namespace Eigen;
	VectorXi clusterIdx;
	graph_kmedoids(SoC, seedsIdx, clusterIdx);
	igl::unique(clusterIdx, clusterCenters); 

	VectorXi idx;
	std::vector<Triplet<double>> IJV;
	IJV.reserve(clusterCenters.size());
	for(int ii = 0; ii < clusterCenters.size(); ii++)
	{	
		igl::find(clusterIdx.array() == clusterCenters(ii), idx); 
		for (int jj = 0; jj < idx.size(); jj++)
		{
			IJV.push_back(Triplet<double>(ii, idx(jj), 1.0));
		}
	}
	K.resize(clusterCenters.size(), clusterIdx.size());
	K.setFromTriplets(IJV.begin(),IJV.end());
	
}