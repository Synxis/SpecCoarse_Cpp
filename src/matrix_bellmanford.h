#ifndef MATRIX_BELLMANFORD_H
#define MATRIX_BELLMANFORD_H
#include <igl/find.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

// Inputs:
//   SoC  SPSD strength of connection matrix
//   seedsIdx  seeds indices
// Outputs: 
//   dist  distance from seeds to all the other points
//   nearestSeeds  nearest seeds indices
void matrix_bellmanford(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::VectorXd & dist,
	Eigen::VectorXi & nearestSeeds);
#endif