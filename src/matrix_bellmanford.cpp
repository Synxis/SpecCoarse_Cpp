#include "matrix_bellmanford.h"

void matrix_bellmanford(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::VectorXd & dist,
	Eigen::VectorXi & nearestSeeds)
{
	using namespace Eigen;

	int nV = SoC.rows();

	nearestSeeds.setConstant(nV, -1);
	igl::slice_into(seedsIdx, seedsIdx, nearestSeeds);

	dist.setConstant(nV, std::numeric_limits<double>::infinity());
 	igl::slice_into(VectorXd::Zero(seedsIdx.size()), seedsIdx, dist);

	VectorXi i, j;
	VectorXd dij;
	igl::find(SoC,i,j,dij);

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