
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
 @file ComputeProlongation_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "ComputeProlongation_ref.hpp"

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeProlongation_ref(const SparseMatrix & Af, Vector & xf) {

  double * xfv = xf.values;
  double * xcv = Af.mgData->xc->values;
  local_int_t ** c2f = Af.mgData->c2fOperator;
  local_int_t nc = Af.mgData->rc->localLength;

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
// TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i=0; i<nc; ++i){
    xfv[c2f[i][0]] += xcv[i];
    if (c2f[i][1] != -1){
      if  (c2f[i][1] < Af.localNumberOfRows)
	xfv[c2f[i][1]] += xcv[i];
      else{
	/*
#ifndef HPCG_NO_MPI
	double value = xcv[i];
	MPI_Send(&value, 1, MPI_DOUBLE, Af.geom->rank+1, 97, MPI_COMM_WORLD);
#endif
	*/
      }
    }
  }

  /*
#ifndef HPCG_NO_MPI
  if (Af.geom->irow_begin % 2 == 1){
    double value;
    MPI_Request request;
    MPI_Irecv(&value, 1, MPI_DOUBLE, Af.geom->rank-1, 97, MPI_COMM_WORLD, &request);
    MPI_Status status;
    if (MPI_Wait(&request, &status)){
      std::exit(-1);
    }
  }
#endif
  */

  return 0;
}
