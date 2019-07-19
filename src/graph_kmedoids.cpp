#include "graph_kmedoids.h"
// #include "matrix_bellmanford.h" // this one is slower

void graph_kmedoids(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::VectorXi & clusterIdx)
{
	using namespace Eigen;
	int maxIter = 20;
	clusterIdx.setConstant(SoC.rows(), -1);

	VectorXi centerIdx(seedsIdx.size());
	centerIdx << seedsIdx;

	VectorXi i, j;
	VectorXd dij;
	igl::find(SoC,i,j,dij);

	MatrixXi ij;
	{
		MatrixXi ijInit(i.size(), 2);
		ijInit << i, j;
		MatrixXi ijSort;
		MatrixXi sortIdx;
		igl::sort(ijInit, 2, true, ijSort, sortIdx);

		VectorXi uniIdxA2C, uniIdxC2A;
		igl::unique_rows(ijSort, ij, uniIdxA2C, uniIdxC2A);	
	}

	VectorXd dist, distFromBorder, distInCluster;
	VectorXi nearestCenter, isBorder, iCenter, jCenter, borderIdx, aggNodes, useless;
	MatrixXi ijBorder;
	VectorXd::Index maxIdx;
	for (int iter = 0; iter < maxIter; iter++)
	{
		// flow from centers to the entire mesh
		matrix_bellmanford_fast(SoC, centerIdx, i, j, dij, dist, nearestCenter);

		// extract boarder indices
		igl::slice(nearestCenter, ij.col(0), iCenter);
		igl::slice(nearestCenter, ij.col(1), jCenter);
		igl::find(iCenter.array() != jCenter.array(), isBorder);
		igl::slice(ij, isBorder, 1, ijBorder);
		borderIdx.resize(ijBorder.rows()*ijBorder.cols());
		borderIdx << ijBorder.col(0), ijBorder.col(1);
		igl::unique(borderIdx, borderIdx); // find unique component in ijVec

		// flow from borders to the entire mesh
		matrix_bellmanford_fast(SoC, borderIdx, i, j, dij, distFromBorder, useless);

		// set points that are far away to the boarders as new centers
		for (int ii = 0; ii < centerIdx.size(); ii++)
		{
			igl::find(nearestCenter.array() == centerIdx(ii), aggNodes); 
			igl::slice(distFromBorder, aggNodes, distInCluster);
			distInCluster.maxCoeff(&maxIdx); // get maximum index
			centerIdx.coeffRef(ii) = aggNodes(maxIdx);
		}

		// stopping criteria
		if ((clusterIdx.array() == nearestCenter.array()).all()) 
		{
			std::cout << "graph K-Medoids iteration: " << iter+1 << std::endl;
			break;
		}	
		else if (iter == maxIter-1)
		{
			clusterIdx << nearestCenter;
			std::cout << "graph K-Medoids iteration: " << maxIter << std::endl;
		}
		else clusterIdx << nearestCenter; // update the cluster index
	}
}

void matrix_bellmanford_fast(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	const Eigen::VectorXi & i,
	const Eigen::VectorXi & j,
	const Eigen::VectorXd & dij,
	Eigen::VectorXd & dist,
	Eigen::VectorXi & nearestSeeds)
{
	using namespace Eigen;

	int nV = SoC.rows();

	nearestSeeds.setConstant(nV, -1);
	igl::slice_into(seedsIdx, seedsIdx, nearestSeeds);

	dist.setConstant(nV, std::numeric_limits<double>::infinity());
 	igl::slice_into(VectorXd::Zero(seedsIdx.size()), seedsIdx, dist);

 	//  VectorXi i, j;
	// VectorXd dij;
	// igl::find(SoC,i,j,dij);

	VectorXd di(i.size()), dj(i.size());
	igl::slice(dist,i,di);
	igl::slice(dist,j,dj);

	VectorXi idx;
	while(true)
	{
		igl::find((di+dij).array() < dj.array(), idx); // MATLAB: idx = find(di+dij < dj);
		if (idx.size() > 0)
		{
			for (int ii = 0; ii < idx.size(); ii++)
			{
				dist.coeffRef(j(idx(ii))) = di.coeff(idx(ii)) + dij.coeff(idx(ii));
				nearestSeeds.coeffRef(j(idx(ii))) = nearestSeeds.coeff(i(idx(ii)));
			}
			igl::slice(dist,i,di);
			igl::slice(dist,j,dj);
		}
		else break;
	}

	assert((nearestSeeds.array() == -1).any() == false);
	assert((dist.array() == std::numeric_limits<double>::infinity()).any() == false);
}