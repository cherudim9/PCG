
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
 @file GenerateProblem.cpp

 HPCG routine
 */

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#include <cassert>
#include "GenerateCoarseProblem.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"

void GenerateCoarseCoef(const SparseMatrix * Af, SparseMatrix * Ac, local_int_t nc, local_int_t * f2cOperator, local_int_t ** c2fOperator){

  local_int_t nf = Af->localNumberOfRows;

  local_int_t numberOfNonzerosPerRow = 400;
  global_int_t totalNumberOfRows = Ac->geom->grow;
  local_int_t * nonzerosInRow = new local_int_t[nc];
  global_int_t ** mtxIndG = new global_int_t * [nc];
  local_int_t ** mtxIndL = new local_int_t * [nc];
  double ** matrixValues = new double * [nc];
  double ** matrixDiagonal = new double * [nc];
  
  Ac->localToGlobalMap.resize(nc);

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for (local_int_t i=0; i<nc; ++i) {
    matrixValues[i] = 0;
    matrixDiagonal[i] = 0;
    mtxIndG[i] = 0;
    mtxIndL[i] = 0;
  }

  for (local_int_t i=0; i< nc; ++i)
    mtxIndL[i] = new local_int_t[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< nc; ++i)
    matrixValues[i] = new double[numberOfNonzerosPerRow];
  for (local_int_t i=0; i< nc; ++i)
    mtxIndG[i] = new global_int_t[numberOfNonzerosPerRow];

  local_int_t localNumberOfNonzeros = 0;

  //why not openMP???
  for(local_int_t ic=0; ic<nc; ic++){
    global_int_t global_ic = Ac->geom->irow_begin + ic;
#ifndef HPCG_NO_OPENMP
#pragma omp critical
#endif
    {
      Ac->globalToLocalMap[global_ic]=ic;
      Ac->localToGlobalMap[ic]=global_ic;
    }
    
    const global_int_t * fine0_mtxIndG = Af->mtxIndG[c2fOperator[ic][0]];
    const double * fine0_matrixValues = Af->matrixValues[c2fOperator[ic][0]];
    local_int_t fine0_nonzerosInRow = Af->nonzerosInRow[c2fOperator[ic][0]];
    const global_int_t * fine1_mtxIndG = (c2fOperator[ic][1]!=-1) ? Af->mtxIndG[c2fOperator[ic][1]]  : 0;
    const double * fine1_matrixValues = (c2fOperator[ic][1]!=-1) ? Af->matrixValues[c2fOperator[ic][1]] : 0;
    local_int_t fine1_nonzerosInRow = (c2fOperator[ic][1]!=-1) ? Af->nonzerosInRow[c2fOperator[ic][1]] : 0;
    
    global_int_t * currentIndexPointerG = mtxIndG[ic];
    double * currentValuePointer = matrixValues[ic];
    local_int_t numOfNonzerosInRow = 0;
    
    global_int_t fine0_col = 0, fine1_col = 0;

    while(fine0_col<fine0_nonzerosInRow || fine1_col<fine1_nonzerosInRow){
      if (fine0_col<fine0_nonzerosInRow && 
	  (fine1_col>=fine1_nonzerosInRow || *fine0_mtxIndG<=*fine1_mtxIndG)
	  ){
	local_int_t coarse0_mtxIndG = *fine0_mtxIndG/2;
	if (numOfNonzerosInRow==0 || *currentIndexPointerG!=coarse0_mtxIndG){
	  if (numOfNonzerosInRow!=0){
	    currentIndexPointerG++;
	    currentValuePointer++;
	  }
	  *currentIndexPointerG = coarse0_mtxIndG;
	  *currentValuePointer = *fine0_matrixValues;
	  numOfNonzerosInRow++;
	  if (coarse0_mtxIndG == global_ic){
	    matrixDiagonal[ic] = currentValuePointer;
	  }
	}else{
	  *currentValuePointer += *fine0_matrixValues;
	}
	fine0_col++;
	fine0_mtxIndG++;
	fine0_matrixValues++;
      }else{
	local_int_t coarse1_mtxIndG = *fine1_mtxIndG/2;
	if (numOfNonzerosInRow==0 || *currentIndexPointerG!=coarse1_mtxIndG){
	  if (numOfNonzerosInRow!=0){
	    currentIndexPointerG++;
	    currentValuePointer++;
	  }
	  *currentIndexPointerG = coarse1_mtxIndG;
	  *currentValuePointer = *fine1_matrixValues;
	  numOfNonzerosInRow++;
	  if (coarse1_mtxIndG == global_ic){
	    matrixDiagonal[ic] = currentValuePointer;
	  }
	}else{
	  *currentValuePointer += *fine1_matrixValues;
	}
	fine1_col++;
	fine1_mtxIndG++;
	fine1_matrixValues++;
      }
    }

    nonzerosInRow[ic] = numOfNonzerosInRow;;
#ifndef HPCG_NO_OPENMP
#pragma omp critical
#endif
    localNumberOfNonzeros += numOfNonzerosInRow;
  }

  global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
  MPI_Allreduce(&localNumberOfNonzeros, &totalNumberOfNonzeros, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif

  Ac->title = 0;
  Ac->totalNumberOfRows = totalNumberOfRows;
  Ac->totalNumberOfNonzeros = totalNumberOfNonzeros;
  Ac->localNumberOfRows = nc;
  Ac->localNumberOfColumns = nc;
  Ac->localNumberOfNonzeros = localNumberOfNonzeros;
  Ac->nonzerosInRow = nonzerosInRow; //***
  Ac->mtxIndG = mtxIndG;             //***
  Ac->mtxIndL = mtxIndL;
  Ac->matrixValues = matrixValues;   //***
  Ac->matrixDiagonal = matrixDiagonal;//***
 
  return;
}

void GenerateCoarseProblem(const SparseMatrix & Af) {

  local_int_t nf = Af.localNumberOfRows;
  local_int_t nc;

  Geometry * geomc = new Geometry();
  geomc->size = Af.geom->size;
  geomc->rank = Af.geom->rank;
  geomc->numThreads = Af.geom->numThreads;
  geomc->grow = Af.geom->grow/2;
  geomc->irow_begin = Af.geom->irow_begin / 2 + 
    ((Af.geom->irow_begin % 2 == 0) ? 0 : 1 );
  geomc->irow_end = Af.geom->irow_end / 2;
  nc = geomc->irow = geomc->irow_end - geomc->irow_begin + 1;

  local_int_t * f2cOperator = new local_int_t[nf];
  local_int_t ** c2fOperator = new local_int_t*[nc];

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for(local_int_t i=0; i<nf; i++)
    f2cOperator[i]=-1;

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i<nf; ++i) {
    int ic=(Af.geom->irow_begin+i)/2;
    if (ic>=geomc->irow_begin && ic<=geomc->irow_end){
      f2cOperator[i] = ic - geomc->irow_begin;
    }
  }

#ifndef HPCG_NO_OPENMP
#pragma omp parallel for
#endif
  for(local_int_t i=0; i<nc; i++){
    int if0 = (i + geomc->irow_begin) * 2 - Af.geom->irow_begin;
    int if1 = if0 + 1;
    c2fOperator[i] = new local_int_t[2];
    c2fOperator[i][0] = if0;
    if (Af.geom->irow_begin + if1 < Af.geom->grow && if1 < nf){
      c2fOperator[i][1] = if1;
    } else
      c2fOperator[i][1] = -1;
  }

  SparseMatrix * Ac = new SparseMatrix;
  InitializeSparseMatrix(*Ac, geomc);
  GenerateCoarseCoef(&Af, Ac, nc, f2cOperator, c2fOperator);
  SetupHalo(*Ac);
  Vector *rc = new Vector;
  Vector *xc = new Vector;
  Vector * Axf = new Vector;
  InitializeVector(*rc, Ac->localNumberOfRows);
  InitializeVector(*xc, Ac->localNumberOfColumns);
  InitializeVector(*Axf, Af.localNumberOfColumns);
  Af.Ac = Ac;
  MGData * mgData = new MGData;
  InitializeMGData(f2cOperator, c2fOperator, rc, xc, Axf, *mgData);
  Af.mgData = mgData;

  return;
}

