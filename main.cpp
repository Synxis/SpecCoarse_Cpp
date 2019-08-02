// libigl includes
#include <igl/read_triangle_mesh.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/randperm.h>

// project code includes
#include "src/strength_of_connection.h"
#include "src/combinatorial_coarsening.h"
#include "src/sp_cwiseSqrt.h"
#include "src/sp_cwiseInverse.h"
#include "src/projection_mat.h"
#include "src/get_sparsity.h"
#include "src/eigs_diagmass.h"
#include "src/init_null_projection.h"
#include "src/null_proj.h"
#include "src/print_spmat.h"

// other includes
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <stdlib.h>

template<typename Derived>
void save_matrix_ascii(const std::string& filename, const Eigen::DenseBase<Derived>& mat, char sep = ',')
{
	std::ofstream out(filename);
	out.precision(14);
	for(Eigen::Index r = 0; r < mat.rows(); r++)
	{
		for(Eigen::Index c = 0; c < mat.cols(); c++)
			if(c == 0) out << mat(r, c);
			else out << sep << mat(r, c);
		out << '\n';
	}
	printf("Saved dense matrix in '%s'\n", filename.c_str());
}

template<typename Scalar, int Options, typename StorageIndex>
void save_matrix_ascii(const std::string& filename, const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat, char sep = '\t')
{
	using SpMat = Eigen::SparseMatrix<Scalar, Options, StorageIndex>;
	std::ofstream out(filename);
	out.precision(14);
	bool max_pos = false;
	for(int outer = 0; outer < mat.outerSize(); outer++)
		for(SpMat::InnerIterator it(mat, outer); it; ++it)
		{
			out << (it.row() + 1) << sep << (it.col() + 1) << sep << it.value() << '\n';
			if(it.row() + 1 == mat.rows() && it.col() + 1 == mat.cols()) max_pos = true;
		}
	if(!max_pos) out << mat.rows() << sep << mat.cols() << sep << "0.0\n";
	printf("Saved sparse matrix in '%s'\n", filename.c_str());
}

std::string meshPath;
int m; // number of vertices in the coarsened mesh
int k; // number of eigenvectors used in the optimization
int maxIter = 5000; // maximum iteration
double lr = 2e-2; // learning rate
int lrReduce = 0; // number of times to reduce learning rate

int main(int argc, char *argv[])
{
	using clock = std::chrono::high_resolution_clock;
	using namespace Eigen;
	using namespace std;
	srand(1);

	// input argument parsing
	if (argc != 6)
	{
		cout << "./specCoarsen_bin [meshPath] [m] [k] [lr] [lrReduce]" << endl;
		return 1;
	}
	else
	{
		meshPath = argv[1];
		m = atoi(argv[2]);
		k = atoi(argv[3]);
		lr = atof(argv[4]);
		lrReduce = atoi(argv[5]);
		cout << "=========================================================" << endl;
		cout << "meshPath: " << meshPath << endl;
		cout << "number of vertices in the coarsened mesh: " << m << endl;
		cout << "number of eigenvectors in use: " << k << endl;
		cout << "learning rate: " << lr << endl;
		cout << "number of learning rate reduction" << lrReduce << endl;
	}

	// load mesh
	MatrixXd V;
	MatrixXi F;
	{
		std::stringstream ss;
		igl::read_triangle_mesh(meshPath, V, F);
		assert(V.rows() > 0);

		// normalize to unit surface area
		Eigen::VectorXd dFA;
		igl::doublearea(V,F,dFA);
		double totalArea = dFA.sum() / 2;
		V /= sqrt(totalArea); // normalize to a unit area
		V *= 50; // increase numberics
		cout << "number of vertices in the original mesh: " << V.rows() << endl;
	}

	const clock::time_point TS = clock::now();

	// compute the Laplace operator and the mass matrix
	SparseMatrix<double> L, M, invM;
	int n; // number of rows/cols of L
	{ 
		igl::cotmatrix(V,F,L);
		igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
		igl::invert_diag(M,invM);
		L = (-L).eval();
		L.makeCompressed();
		n = L.rows();
	}

	// =========================
	// combinatorial coarsening
	// =========================
	// compute the strength of connections (aka operator distance)
	SparseMatrix<double> SoC;
	{
		strength_of_connection(L,M,SoC);
		SoC.makeCompressed();
		assert(SoC.nonZeros() > L.rows()); // SoC should have a lot of nonzeros. If not, probably the sign of L is wrong
	}

	// k-medoids clustering
	SparseMatrix<double> K; // assignment matrix
	SparseMatrix<double> P; // projection mateix
	VectorXi Cpts; // coarse vertex indices
	{
		// random seeds
		VectorXi seedsIdx;
		igl::randperm(L.rows(), seedsIdx);
		seedsIdx = seedsIdx.tail(m);

		// clustering
		combinatorial_coarsening(SoC, seedsIdx, K, Cpts);

		// construct projection matrix
		projection_mat(L, Cpts, P);
		P.makeCompressed();
	}
	
	// =========================
	// operator optimization
	// =========================
	// construct coarsened mass matrix
	SparseMatrix<double> Mc, invMc, sqrtMc;
	{
		Mc = K * M * K.transpose();
		Mc.makeCompressed();
		igl::invert_diag(Mc,invMc);
		sp_cwiseSqrt(Mc, sqrtMc);
	}

	// construct the sparsity matrix of G
	SparseMatrix<double> SG;
	{
		SparseMatrix<double> SL;
		get_sparsity(L, SL);

		SparseMatrix<double> Ac;
		get_sparsity(K * SL * K.transpose(), Ac);

		get_sparsity(K.transpose()*Ac, SG);
	}

	// compute test vectors (eigenvectors)
	MatrixXd U; 
	{
		VectorXd eVal;
		eigs_diagmass(L, M, k, eVal, U);
	}

	// pre-computation for optimization
	MatrixXd A(m,k), B(m,k); 
	MatrixXd AB(m,m), BB(m,m);
	{ 
		// (1) E = || Mc^(1/2) * (A - Mc^(-1) * G'LG * B) ||^2
		A = P * invM * L * U;
		B = P * U;
		AB = A * B.transpose();
		BB = B * B.transpose();
	}

	// initialize parameters for the nullspace projection
	VectorXd nullVec;       
	nullVec.setOnes(n);
	VectorXd nullVecCoarse = P * nullVec; 
	SparseMatrix<double> d, dGMat;
	VectorXd invd2;
	init_null_projection(SG, P, nullVec, d, invd2);
	SparseMatrix<double> deltaGMat(SG.rows(), SG.cols());
	VectorXd dG, projStepSizeVec;

	// place holders for gradient computation
	MatrixXd gradRHS(m, m), BBGLGinvMc(m, m);
	SparseMatrix<double, RowMajor> gradLHS(n, m);
	SparseMatrix<double> dEdG = SG;
	SparseMatrix<double> GLG, LG;
	GLG.resize(m, m);
	LG.resize(n, m);

	// optimization 
	SparseMatrix<double> G = K.transpose();
	SparseMatrix<double> Gbest = K.transpose();
	SparseMatrix<double> Gtemp = K.transpose();

	// for stopping criteria
	double objValBest = 1e10;
	double objVal;
	int stallCount = 0;
	int reduceLrCount = 0;

	// NADAM
	double timeStep = 0;
	double beta1 = 0.9;
	double beta2 = 0.9;
	double eps = 1e-8;
	SparseMatrix<double> mt(G.rows(), G.cols());     mt.setZero(); 
	SparseMatrix<double> nt(G.rows(), G.cols());     nt.setZero(); 
	SparseMatrix<double> mt_hat(G.rows(), G.cols()); mt_hat.setZero(); 
	SparseMatrix<double> nt_hat(G.rows(), G.cols()); nt_hat.setZero(); 
	SparseMatrix<double> nt_hatSqrt(G.rows(), G.cols()); nt_hatSqrt.setZero(); 
	SparseMatrix<double> nt_hatSqrtInv(G.rows(), G.cols()); nt_hatSqrtInv.setZero(); 

	// running the operator optimization (output "Gbest")
	for(int iter = 0; iter < maxIter; iter++)
	{
		// precompute GLG and LG
		LG = L * G;
		GLG = G.transpose() * LG;

		// compute objective function
		objVal = (sqrtMc * (A - invMc * GLG * B)).squaredNorm();
		if ((iter % 10) == 0) 
			cout << "iter: " << iter << ", objective value: " << objVal << endl;
		
		// compute gradient
		gradLHS = LG;
		BBGLGinvMc = BB*(GLG*invMc);
		gradRHS = - AB - AB.transpose() + BBGLGinvMc + BBGLGinvMc.transpose();
		// coeffRef is expensive, use it.valueRef instead
		for (int k=0; k<dEdG.outerSize(); ++k) 
			for (SparseMatrix<double>::InnerIterator it(dEdG,k); it; ++it){
				it.valueRef() = gradLHS.row(it.row()).dot(gradRHS.col(it.col()));
			}

		// NADAM update
		timeStep += 1;
		mt *= beta1;
		mt += (1 - beta1) * dEdG;
		nt *= beta2;
		nt += (1 - beta2) * dEdG.cwiseProduct(dEdG);
		mt_hat = mt / (1-pow(beta1, timeStep));
		nt_hat = nt / (1-pow(beta2, timeStep));
		nt_hatSqrt = nt_hat.cwiseSqrt();
		sp_cwiseInverse(nt_hatSqrt, nt_hatSqrtInv);
		dEdG = (beta1*mt_hat + (1-beta1)*dEdG/(1-pow(beta1,timeStep))).cwiseProduct(nt_hatSqrtInv);

		// update G
		Gtemp = G - lr * dEdG;
		null_proj(Gtemp, nullVec, d, invd2, G);

		// check feasibility
		double feasibility = (G * nullVecCoarse - nullVec).cwiseAbs().maxCoeff();
		if (feasibility > 1e-7)
			cout << "G is not feasible" << endl;

		// check stopping criteria
		if(objValBest < objVal)
		{
			stallCount += 1;
			if (stallCount > 10)
			{
				if (reduceLrCount < lrReduce)
				{
					reduceLrCount ++;
					cout << "reduce learning rate" << endl;
					lr *= 0.5;
					stallCount = 0;
				}
				else
				{
					cout << "iter " << iter << ", best objVal: " << objValBest << endl;
					break;
				}
			}
		}
		else // objVal < objValOld
		{
			Gbest = G;
			stallCount = 0;
			objValBest = objVal;
		}
	}

	const clock::time_point TE = clock::now();

	// =========================
	// evaluation: functional map
	// =========================
	// evaluation
	{
		SparseMatrix<double> Lc = Gbest.transpose() * L * Gbest;

		MatrixXd eVec, eVecc; 
		{
			VectorXd eVal, eValc;
			eigs_diagmass(L, M, k, eVal, eVec);
			eigs_diagmass(Lc, Mc, k, eValc, eVecc);
		}

		SparseMatrix<double> ptMap;
		{
			vector<Triplet<double>> IJV;
			IJV.reserve(m);
			for(int ii = 0; ii < m; ii++)
				IJV.push_back(Triplet<double>(Cpts(ii), ii, 1.0));
			ptMap.resize(L.rows(), m);
			ptMap.setFromTriplets(IJV.begin(),IJV.end());
		}
		MatrixXd fMap = eVecc.transpose() * Mc * ptMap.transpose() * eVec;

		save_matrix_ascii("result-Mc.txt", Mc);
		save_matrix_ascii("result-Lc.txt", Lc);
		save_matrix_ascii("result-P.txt", P);
		save_matrix_ascii("output_fmap.txt", fMap);
		std::cout << "Time: " << (std::chrono::duration_cast<std::chrono::milliseconds>(TE - TS).count() / 1000.0) << " s\n";
	}
}
