
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
 @file GenerateGeometry.cpp

 HPCG routine
 */

#include <cmath>
#include <cstdlib>
#include <cassert>

#include "ComputeOptimalShapeXYZ.hpp"
#include "GenerateGeometry.hpp"

#ifdef HPCG_DEBUG
#include <fstream>
#include "hpcg.hpp"
using std::endl;

#endif

/*!
  Computes the factorization of the total number of processes into a
  3-dimensional process grid that is as close as possible to a cube. The
  quality of the factorization depends on the prime number structure of the
  total number of processes. It then stores this decompostion together with the
  parallel parameters of the run in the geometry data structure.

  @param[in]  size total number of MPI processes
  @param[in]  rank this process' rank among other MPI processes
  @param[in]  numThreads number of OpenMP threads in this process
  @param[in]  pz z-dimension processor ID where second zone of nz values start
  @param[in]  nx, ny, nz number of grid points for each local block in the x, y, and z dimensions, respectively
  @param[out] geom data structure that will store the above parameters and the factoring of total number of processes into three dimensions
*/
void GenerateGeometry(int size, int rank, int numThreads, int nrow,
  Geometry * geom)
{
  geom->size = size;
  geom->rank = rank;
  geom->numThreads = numThreads;
  
  geom->grow=nrow;
  geom->irow=nrow/size;
  if (rank<nrow%size){
    geom->irow++;
    geom->irow_begin=rank*(nrow/size)+rank;
  }else{
    geom->irow_begin=rank*(nrow/size)+nrow%size;
  }
  
  geom->irow_end=geom->irow_begin+geom->irow-1;

  return;
}
