#include "init_null_projection.h"

void init_null_projection(
	const Eigen::SparseMatrix<double> & S,
	const Eigen::SparseMatrix<double> & P,
	const Eigen::VectorXd & nullVec,
	Eigen::SparseMatrix<double> & d,
	Eigen::VectorXd & invd2)
{
	using namespace Eigen;

	VectorXd nullVecCoarse = P * nullVec;

	VectorXi r, c;
	VectorXd _;
	igl::find(S,r,c,_);

	std::vector<Triplet<double>> IJV;
	IJV.reserve(r.size());
	for(int ii = 0; ii < r.size(); ii++)
			IJV.push_back(Triplet<double>(r(ii), c(ii), nullVecCoarse(c(ii))));
	d.resize(S.rows(), S.cols());
	d.setFromTriplets(IJV.begin(),IJV.end());

	SparseMatrix<double> d2 = d.cwiseProduct(d);
	VectorXd sumd2;
	igl::sum(d2,2,sumd2);
	invd2 = sumd2.array().inverse();  
}