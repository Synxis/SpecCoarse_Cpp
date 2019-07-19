#include <igl/readOBJ.h>
// #include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/floor.h>
#include <igl/randperm.h>
#include <igl/slice.h>
#include <igl/sum.h>
#include <igl/doublearea.h>
#include <igl/diag.h>

#include <strength_of_connection.h>
#include <graph_kmeans.h>
#include <lloyd_aggregation.h>
#include <sparsity_mat.h>
#include <projection_mat.h>
#include <eigs_diagmass.h>
#include <init_null_projection.h>
#include <sp_cwiseInverse.h>
#include <sp_cwiseSqrt.h>
#include <null_proj.h>
#include <objFunc.h>

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <profc.h>
#include <stdlib.h>

#ifndef MESH_PATH
#define MESH_PATH "../../runtimeMeshes/"
#endif

std::string meshName;
int numNc, numTestVectors;
int maxIter = 3000;
double lr = 2e-2;

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;
	srand(1);

 	if (argc != 5)
 	{
 		cout << "./optPreserv_bin [meshName] [numNc] [numTestVectors] [lr]" << endl;
 		assert(false);
 	}
	else
	{
		meshName = argv[1];
		numNc = atoi(argv[2]);
		numTestVectors = atoi(argv[3]);
		lr = atof(argv[4]);
		cout << "===================" << endl;
		cout << "meshName: " << meshName << endl;
		cout << "numNc: " << numNc << endl;
		cout << "numTestVectors: " << numTestVectors << endl;
	}

	// load mesh
	MatrixXd V;
	MatrixXi F;
	{
		// PROFC_NODE("Read mesh")
		std::stringstream ss;
		ss << MESH_PATH + meshName + ".obj";
		igl::readOBJ(ss.str(), V, F);
		assert(V.rows() > 0);
		
		// normalize wrt surf area
		Eigen::VectorXd dFA;
		igl::doublearea(V,F,dFA);
		double totalArea = dFA.sum() / 2;
		V /= sqrt(totalArea); // normalize to a unit area
		V *= 50;
		cout << "numV: " << V.rows() << endl;
	}

	// compute system matrix mass matrix
	SparseMatrix<double> L, M, invM;
	{ 
		// PROFC_NODE("compute operator")
		igl::cotmatrix(V,F,L);
		igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
		igl::invert_diag(M,invM);
		L = (-L).eval();
		L.makeCompressed();
	}

	// compute strength of connections
	SparseMatrix<double> SoC;
	{
		strength_of_connection(L,M,SoC);
		SoC.makeCompressed();
		assert(SoC.nonZeros() > L.rows()); // SoC should have a lot of nonzeros. If not, probably the sign of L is wrong
	}

	// random seeds
	VectorXi seedsIdx;
	igl::randperm(L.rows(), seedsIdx);
	seedsIdx = seedsIdx.tail(numNc);

	// lloyd aggregation
	SparseMatrix<double> Gt; // tentative interpolation
	VectorXi Cpts; // center points
	{
		// PROFC_NODE("graph kmeans")
		lloyd_aggregation(SoC, seedsIdx, Gt, Cpts);
	}

	// Coarse grid mass matrix
	SparseMatrix<double> Mc, invMc, sqrtMc;
	{
		Mc = Gt.transpose() * M * Gt;
		Mc.makeCompressed();
		igl::invert_diag(Mc,invMc);
		sp_cwiseSqrt(Mc, sqrtMc);
	}

	// construct sparsity pattern
	SparseMatrix<double> S;
	sparsity_mat(SoC, Cpts, Gt, S);
	S.makeCompressed();

	// projection matrix
	SparseMatrix<double> P;
	projection_mat(L, Cpts, P);
	P.makeCompressed();

	// compute test vectors (eigenvectors)
	MatrixXd U; 
	{
		PROFC_NODE("eigs")
		VectorXd eVal;
		eigs_diagmass(L, M, numTestVectors, eVal, U);
	}

	// pre-computation for optimization
	MatrixXd A(numNc,numTestVectors), B(numNc,numTestVectors); 
	MatrixXd AB(numNc, numNc), BB(numNc, numNc);
	{ 
		// PROFC_NODE("initialize optimization")
		// (1) E = || Mc^(1/2) * (A - Mc^(-1) * G'LG * B) ||^2
		A = P * invM * L * U;
		B = P * U;
		AB = A * B.transpose();
		BB = B * B.transpose();
	}

	// initialize parameters for nullspace projection
	VectorXd nullVec;       
	nullVec.setOnes(L.rows());
	VectorXd nullVecCoarse = P * nullVec; 
	SparseMatrix<double> d, dGMat;
	VectorXd invd2;
	init_null_projection(S, P, nullVec, d, invd2);
	SparseMatrix<double> deltaGMat(Gt.rows(), Gt.cols());
	VectorXd dG, projStepSizeVec;

	// place holders for gradient computation
	MatrixXd gradRHS(numNc, numNc), BBGLGinvMc(numNc, numNc);
	SparseMatrix<double, RowMajor> gradLHS(L.rows(), numNc);
	SparseMatrix<double> dEdG = S;
	SparseMatrix<double> dEdG_dir;
	SparseMatrix<double> GLG, LG;
	GLG.resize(numNc, numNc);
	LG.resize(L.rows(), numNc);

	// optimization 
	SparseMatrix<double> G = Gt;
	SparseMatrix<double> Gbest = Gt;
	SparseMatrix<double> Gtemp = Gt;

	// for stopping criteria
	double objValBest = 1e10, objVal;
	int stallCount = 0;
	int reduceLrCount = 0;
	int numReduce = 0;

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
	
	time_t tstart; 
	tstart = clock();
	for(int iter = 0; iter < maxIter; iter++)
	{
		{ // precompute GLG and LG
			// PROFC_NODE("GLG & LG")
			LG = L * G;
			GLG = G.transpose() * LG;
		}

		{ // compute objective function
			// PROFC_NODE("compute objective")
			objVal = (sqrtMc * (A - invMc * GLG * B)).squaredNorm();
			// if ((iter % 1) == 0) 
			// 	cout << "iter: " << iter << ", objective value: " << objVal << endl;
		}
		
		{ // compute gradient
			// PROFC_NODE("compute gradient")
			gradLHS = LG;
			BBGLGinvMc = BB*(GLG*invMc);
			gradRHS = - AB - AB.transpose() + BBGLGinvMc + BBGLGinvMc.transpose();
			// coeffRef is expensive, use it.valueRef instead
			for (int k=0; k<dEdG.outerSize(); ++k) 
				for (SparseMatrix<double>::InnerIterator it(dEdG,k); it; ++it){
					it.valueRef() = gradLHS.row(it.row()).dot(gradRHS.col(it.col()));
				}
		}

		{ // NADAM update
			// PROFC_NODE("NADAM")
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
		}

		{ // update
			// PROFC_NODE("update G")
			Gtemp = G - lr * dEdG;
			null_proj(Gtemp, nullVec, d, invd2, G);
		}

		// { // check feasibility
		// 	PROFC_NODE("feasibility & stopping")
		// 	double feasibility = (G * nullVecCoarse - nullVec).cwiseAbs().maxCoeff();
		// 	if (feasibility > 1e-10)
		// 		cout << "G is not feasible" << endl;
		// }

		{ // check stop
			// PROFC_NODE("stop")
			{
				if(objValBest < objVal)
				{
					stallCount += 1;
					if (stallCount > 10)
					{
						if (reduceLrCount < numReduce)
						{
							reduceLrCount ++;
							cout << "reduce learning rate" << endl;
							lr *= 0.5;
							stallCount = 0;
							// mt.setZero(); 
						 	// nt.setZero(); 
							// mt_hat.setZero(); 
							// nt_hat.setZero(); 
							// nt_hatSqrt.setZero(); 
							// nt_hatSqrtInv.setZero(); 
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
		}
	}
	cout << float(clock() - tstart) / CLOCKS_PER_SEC << " sec."<< endl;

	// // evaluation
	// SparseMatrix<double> Lc = Gbest.transpose() * L * Gbest;

	// MatrixXd eVec, eVecc; 
	// int fMapSize = numTestVectors;
	// {
	// 	VectorXd eVal, eValc;
	// 	eigs_diagmass(L, M, fMapSize, eVal, eVec);
	// 	eigs_diagmass(Lc, Mc, fMapSize, eValc, eVecc);
	// }

	// SparseMatrix<double> ptMap;
	// {
	// 	vector<Triplet<double>> IJV;
	// 	IJV.reserve(numNc);
	// 	for(int ii = 0; ii < numNc; ii++)
	// 		IJV.push_back(Triplet<double>(Cpts(ii), ii, 1.0));
	// 	ptMap.resize(L.rows(), numNc);
	// 	ptMap.setFromTriplets(IJV.begin(),IJV.end());
	// }
	// MatrixXd fMap = eVecc.transpose() * Mc * ptMap.transpose() * eVec;

	// ofstream fMatFile;
	// std::stringstream outss;
	// outss << "../" + meshName + "_fmap.txt";
	// fMatFile.open (outss.str());
	// for (int ii = 0; ii < fMapSize; ii ++)
	// 	fMatFile << fMap.row(ii).head(fMapSize) << "\n";
	// fMatFile.close();

  // // Plot the mesh
  // igl::opengl::glfw::Viewer viewer;
  // viewer.data().set_mesh(V, F);

  // Eigen::MatrixXd C;
  // // igl::colormap(igl::COLOR_MAP_TYPE_JET,eVec.col(1),false,C);
  // viewer.data().set_colors(U.col(0));

  // // Eigen::MatrixXd points;
  // // igl::slice(V, Cpts, 1, points);
  // // RowVector3d ptColor;
  // // ptColor << 144./255., 210./255., 236./255.;
  // // viewer.data().set_points(points, ptColor);

  // viewer.launch();
}
