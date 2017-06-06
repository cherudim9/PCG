
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

/*!
 @file GenerateProblem_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#if defined(HPCG_DEBUG) || defined(HPCG_DETAILED_DEBUG)
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#endif
#include <cassert>
#include <iostream>
#include <algorithm>

#include "GenerateProblem_ref.hpp"

void GenerateProblem_ref(SparseMatrix & A, Vector * b, Vector * x, std::map<int, std::vector<std::pair<int, double> > > & id_to_neighbors, std::map<int, double> & id_to_b){
  //the ID of id_to_neighbors is 1-indexed, as is id_to_b

  local_int_t localNumberOfRows = A.geom->irow;
  assert(localNumberOfRows>0); // Throw an exception of the number of rows is less than zero (can happen if int overflow)
  local_int_t numberOfNonzerosPerRow = 50; // We are approximating a 27-point finite element/volume/difference 3D stencil

  global_int_t totalNumberOfRows = A.geom->grow; // Total number of grid points in mesh
  assert(totalNumberOfRows>0); // Throw an exception of the number of rows is less than zero (can happen if int overflow)

  // Allocate arrays that are of length localNumberOfRows
  char * nonzerosInRow = new char[localNumberOfRows];
  global_int_t ** mtxIndG = new global_int_t*[localNumberOfRows];
  local_int_t  ** mtxIndL = new local_int_t*[localNumberOfRows];
  double ** matrixValues = new double*[localNumberOfRows];
  double ** matrixDiagonal = new double*[localNumberOfRows];

  if (b!=0) InitializeVector(*b, localNumberOfRows);
  if (x!=0) InitializeVector(*x, localNumberOfRows);
  double * bv = 0;
  double * xv = 0;
  if (b!=0) bv = b->values; // Only compute exact solution if requested
  if (x!=0) xv = x->values; // Only compute exact solution if requested
  A.localToGlobalMap.resize(localNumberOfRows);

  // Use a parallel loop to do initial assignment:
  // distributes the physical placement of arrays of pointers across the memory system
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; ++i) {
    matrixValues[i] = 0;
    matrixDiagonal[i] = 0;
    mtxIndG[i] = 0;
    mtxIndL[i] = 0;
  }

#ifndef HPCG_CONTIGUOUS_ARRAYS
  // Now allocate the arrays pointed to
  for (local_int_t i=0; i< localNumberOfRows; ++i)
    mtxIndL[i] = new local_int_t[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< localNumberOfRows; ++i)
    matrixValues[i] = new double[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< localNumberOfRows; ++i)
   mtxIndG[i] = new global_int_t[numberOfNonzerosPerRow];

#else
  // Now allocate the arrays pointed to
  mtxIndL[0] = new local_int_t[localNumberOfRows * numberOfNonzerosPerRow];
  matrixValues[0] = new double[localNumberOfRows * numberOfNonzerosPerRow];
  mtxIndG[0] = new global_int_t[localNumberOfRows * numberOfNonzerosPerRow];

  for (local_int_t i=1; i< localNumberOfRows; ++i) {
  mtxIndL[i] = mtxIndL[0] + i * numberOfNonzerosPerRow;
  matrixValues[i] = matrixValues[0] + i * numberOfNonzerosPerRow;
  mtxIndG[i] = mtxIndG[0] + i * numberOfNonzerosPerRow;
  }
#endif

  local_int_t localNumberOfNonzeros = 0;
  // TODO:  This triply nested loop could be flattened or use nested parallelism
#ifndef HPCG_NO_OPENMP
  //#pragma omp parallel for
#endif
  for (local_int_t irow_id=0; irow_id<localNumberOfRows; irow_id++) {
    global_int_t grow_id = A.geom->irow_begin + irow_id;
    local_int_t currentLocalRow = irow_id;
    global_int_t currentGlobalRow = grow_id;
    //std::cout<<"id="<<grow_id<<std::endl;
    //std::cout<<"n_neighbors="<<id_to_neighbors[grow_id+1].size()<<std::endl;
#ifndef HPCG_NO_OPENMP
// C++ std::map is not threadsafe for writing
#pragma omp critical
#endif
    A.globalToLocalMap[currentGlobalRow] = currentLocalRow;
    A.localToGlobalMap[currentLocalRow] = currentGlobalRow;
#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout << " rank, globalRow, localRow = " << A.geom->rank << " " << currentGlobalRow << " " << A.globalToLocalMap[currentGlobalRow] << endl;
#endif
    char numberOfNonzerosInRow = 0;
    double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
    global_int_t * currentIndexPointerG = mtxIndG[currentLocalRow]; // Pointer to current index in current row

    int last_id = -1;

    for(int i1=0; i1<id_to_neighbors[grow_id+1].size(); i1++){
      int neighbor_id = id_to_neighbors[grow_id+1][i1].first - 1;
      
      //sanity check
      if (last_id >= neighbor_id){
	std::cout<< "wrong!!! neighbor not ordered" << std::endl;
      }

      double neighbor_rv = id_to_neighbors[grow_id+1][i1].second;
      if (neighbor_id == grow_id){
	matrixDiagonal[currentLocalRow] = currentValuePointer;
	*currentValuePointer++ = neighbor_rv;
      }else{
	*currentValuePointer++ = neighbor_rv;
      }
      *currentIndexPointerG++ = neighbor_id;
      numberOfNonzerosInRow++;
    }

    nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
#ifndef HPCG_NO_OPENMP
#pragma omp critical
#endif
    localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
    if (b!=0)      bv[currentLocalRow] = id_to_b[currentGlobalRow+1];
    if (x!=0)      xv[currentLocalRow] = 1.0;
  }

  std::cout     << "Process " << A.geom->rank << " of " << A.geom->size <<" has " << localNumberOfRows    << " rows."     << std::endl
		<< "Process " << A.geom->rank << " of " << A.geom->size <<" has " << localNumberOfNonzeros<< " nonzeros." <<std::endl;

  global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to sum all nonzeros
  MPI_Allreduce(&localNumberOfNonzeros, &totalNumberOfNonzeros, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif
  // If this assert fails, it most likely means that the global_int_t is set to int and should be set to long long
  // This assert is usually the first to fail as problem size increases beyond the 32-bit integer range.
  assert(totalNumberOfNonzeros>0); // Throw an exception of the number of nonzeros is less than zero (can happen if int overflow)

  A.title = 0;
  A.totalNumberOfRows = totalNumberOfRows;
  A.totalNumberOfNonzeros = totalNumberOfNonzeros;
  A.localNumberOfRows = localNumberOfRows;
  A.localNumberOfColumns = localNumberOfRows;
  A.localNumberOfNonzeros = localNumberOfNonzeros;
  A.nonzerosInRow = nonzerosInRow;
  A.mtxIndG = mtxIndG;
  A.mtxIndL = mtxIndL;
  A.matrixValues = matrixValues;
  A.matrixDiagonal = matrixDiagonal;

  return;
}
