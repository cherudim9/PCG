
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

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
const char* NULLDEVICE="nul";
#else
const char* NULLDEVICE="/dev/null";
#endif

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include "hpcg.hpp"


std::ofstream HPCG_fout; //!< output file stream for logging activities during HPCG run

static int
startswith(const char * s, const char * prefix) {
  size_t n = strlen( prefix );
  if (strncmp( s, prefix, n ))
    return 0;
  return 1;
}

int
HPCG_Init(int * argc_p, char ** *argv_p, HPCG_Params & params) {
  int argc = *argc_p;
  char ** argv = *argv_p;

#ifndef HPCG_NO_MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &params.comm_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &params.comm_size );
#else
  params.comm_rank = 0;
  params.comm_size = 1;
#endif

#ifdef HPCG_NO_OPENMP
  params.numThreads = 1;
#else
  #pragma omp parallel
  params.numThreads = omp_get_num_threads();
#endif

  time_t rawtime;
  time ( &rawtime );
  tm * ptm = localtime(&rawtime);
  char fname[80];
  sprintf( fname, "project%04d%02d%02dT%02d%02d%02d.txt",
      1900 + ptm->tm_year, ptm->tm_mon+1, ptm->tm_mday, ptm->tm_hour, ptm->tm_min, ptm->tm_sec );

  if (0 == params.comm_rank) {
    HPCG_fout.open(fname);
  } else {
#if defined(HPCG_DEBUG) || defined(HPCG_DETAILED_DEBUG)
    sprintf( fname, "project%04d%02d%02dT%02d%02d%02d_%d.txt",
        1900 + ptm->tm_year, ptm->tm_mon+1, ptm->tm_mday, ptm->tm_hour, ptm->tm_min, ptm->tm_sec, params.comm_rank );
    HPCG_fout.open(fname);
#else
    HPCG_fout.open(NULLDEVICE);
#endif
  }

  return 0;
}
