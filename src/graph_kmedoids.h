#ifndef GRAPH_KMEDOIDS_H
#define GRAPH_KMEDOIDS_H
#include <igl/sort.h>
#include <igl/find.h>
#include <igl/unique_rows.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

// Inputs:
//   SoC  SPSD strength of connection matrix
//   seedsIdx  seeds indices
// Outputs: 
//   clusterIdx the cluster index for each point
void graph_kmedoids(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	Eigen::VectorXi & clusterIdx);

// this is only for the graph_kmedoids
void matrix_bellmanford_fast(
	const Eigen::SparseMatrix<double> & SoC,
	const Eigen::VectorXi & seedsIdx,
	const Eigen::VectorXi & i,
	const Eigen::VectorXi & j,
	const Eigen::VectorXd & dij,
	Eigen::VectorXd & dist,
	Eigen::VectorXi & nearestSeeds);
#endif