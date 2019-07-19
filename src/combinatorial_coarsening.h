#ifndef COMBINATORIAL_COARSENING_H
#define COMBINATORIAL_COARSENING_H
#include <igl/find.h>
#include <igl/unique.h>
#include <vector>
#include "graph_kmedoids.h"

// Inputs:
//   SoC  symmetric strength of connections matrix
//   seedsIdx  indices of initial seeds
// Outputs: 
//   K  assignment matrix assigning original vertices to coarse vertices
//   clusterCenters  cluster centers (corresponding to rows of Gt)
void combinatorial_coarsening(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::SparseMatrix<double> & K,
	Eigen::VectorXi & clusterCenters);
#endif