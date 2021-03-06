
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

#ifndef GENERATEPROBLEM_HPP
#define GENERATEPROBLEM_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include <vector>
#include <map> 
#include <cstdlib>

void GenerateProblem(SparseMatrix & A, Vector * b, Vector * x, std::map<int, std::vector<std::pair<int, double> > > & id_to_neighbors, std::map<int, double> & id_to_b);
#endif // GENERATEPROBLEM_HPP
