
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
 @file SetupHalo_ref.cpp

 HPCG routine
 */

#ifndef HPCG_NO_MPI
#include <mpi.h>
#include <map>
#include <set>
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#include <cassert>
#endif

#include "SetupHalo_ref.hpp"
#include "mytimer.hpp"
#include <iostream>

/*!
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/

#ifndef HPCG_NO_MPI
void ExchangeRowRange(SparseMatrix & A, int * & RowRange){
  int size, rank; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int MPI_MY_TAG = 99;
  MPI_Request * request = new MPI_Request[size];
  
  RowRange = new int[2*size];
  RowRange[2*rank]=A.geom->irow_begin;
  RowRange[2*rank+1]=A.geom->irow_end;

  for(int i=0; i<size; i++)
    if (i!=rank){
      MPI_Irecv(RowRange+2*i, 2, MPI_INT, i, MPI_MY_TAG, MPI_COMM_WORLD, request+i);
    }

  int * sendBuffer = RowRange+2*rank;
  for(int i=0; i<size; i++)
    if (i!=rank){
      MPI_Send(sendBuffer, 2, MPI_INT, i, MPI_MY_TAG, MPI_COMM_WORLD);
    }

  MPI_Status status;
  for(int i=0; i<size; i++)
    if (i!=rank)
      if (MPI_Wait(request+i, &status))
	std::exit(-1);

  delete []request;

  std::cout<<"Rank#"<<rank<<"   ";
  for(int i=0; i<size; i++)
    std::cout<<i<<":"<<RowRange[i*2]<<"->"<<RowRange[i*2+1]<<", ";
  std::cout<<std::endl;
}
#endif

int GetRankOfMatrixRow(Geometry * geom, int * RowRange, int curIndex){
  int size=geom->size;
  int rank=geom->rank;
  for(int i=0; i<size; i++)
    if (curIndex>=RowRange[2*i] && curIndex<=RowRange[2*i+1])
      return i;
  return -1;
}

void SetupHalo_ref(SparseMatrix & A) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  local_int_t  * nonzerosInRow = A.nonzerosInRow;
  global_int_t ** mtxIndG = A.mtxIndG;
  local_int_t ** mtxIndL = A.mtxIndL;

#ifdef HPCG_NO_MPI  // In the non-MPI case we simply copy global indices to local index storage
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    int cur_nnz = nonzerosInRow[i];
    for (int j=0; j<cur_nnz; j++) mtxIndL[i][j] = mtxIndG[i][j];
  }

#else // Run this section if compiling for MPI

  int * RowRange = 0;
  ExchangeRowRange(A, RowRange);
  assert(RowRange!=0);
  
  std::map< int, std::set< global_int_t> > sendList, receiveList;
  typedef std::map< int, std::set< global_int_t> >::iterator map_iter;
  typedef std::set<global_int_t>::iterator set_iter;
  std::map< local_int_t, local_int_t > externalToLocalMap;

  // TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    global_int_t currentGlobalRow = A.localToGlobalMap[i];
    for (int j=0; j<nonzerosInRow[i]; j++) {
      global_int_t curIndex = mtxIndG[i][j];
      //int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
      int rankIdOfColumnEntry = GetRankOfMatrixRow(A.geom, RowRange, curIndex);
      assert(rankIdOfColumnEntry!=-1);
#ifdef HPCG_DETAILED_DEBUG
      HPCG_fout << "rank, row , col, globalToLocalMap[col] = " << A.geom->rank << " " << currentGlobalRow << " "
          << curIndex << " " << A.globalToLocalMap[curIndex] << endl;
#endif
      if (A.geom->rank!=rankIdOfColumnEntry) {// If column index is not a row index, then it comes from another processor
        receiveList[rankIdOfColumnEntry].insert(curIndex);
        sendList[rankIdOfColumnEntry].insert(currentGlobalRow); // Matrix symmetry means we know the neighbor process wants my value
      }
    }
  }

  // Count number of matrix entries to send and receive
  local_int_t totalToBeSent = 0;
  for (map_iter curNeighbor = sendList.begin(); curNeighbor != sendList.end(); ++curNeighbor) {
    totalToBeSent += (curNeighbor->second).size();
  }
  local_int_t totalToBeReceived = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor) {
    totalToBeReceived += (curNeighbor->second).size();
  }

  int global_tbs = 0, global_tbr = 0;
  MPI_Allreduce(&totalToBeSent, &global_tbs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&totalToBeReceived, &global_tbr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (A.geom->rank == 0){
    int nonzeros = A.totalNumberOfNonzeros;
    std::cout << "totalToBeSent = " << global_tbs << " totalToBeReceived = " << global_tbr << std::endl;
    std:: cout << " ratio  = " << global_tbs + global_tbr << " / " << nonzeros << " = " << (global_tbs+global_tbr+0.0)/nonzeros*100.0 << "%"<<std::endl;
  }
  
  // Build the arrays and lists needed by the ExchangeHalo function.
  double * sendBuffer = new double[totalToBeSent];
  local_int_t * elementsToSend = new local_int_t[totalToBeSent];
  int * neighbors = new int[sendList.size()];
  local_int_t * receiveLength = new local_int_t[receiveList.size()];
  local_int_t * sendLength = new local_int_t[sendList.size()];
  int neighborCount = 0;
  local_int_t receiveEntryCount = 0;
  local_int_t sendEntryCount = 0;
  for (map_iter curNeighbor = receiveList.begin(); curNeighbor != receiveList.end(); ++curNeighbor, ++neighborCount) {
    int neighborId = curNeighbor->first; // rank of current neighbor we are processing
    neighbors[neighborCount] = neighborId; // store rank ID of current neighbor
    receiveLength[neighborCount] = receiveList[neighborId].size();
    sendLength[neighborCount] = sendList[neighborId].size(); // Get count if sends/receives
    for (set_iter i = receiveList[neighborId].begin(); i != receiveList[neighborId].end(); ++i, ++receiveEntryCount) {
      externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount; // The remote columns are indexed at end of internals
    }
    for (set_iter i = sendList[neighborId].begin(); i != sendList[neighborId].end(); ++i, ++sendEntryCount) {
      //if (geom.rank==1) HPCG_fout << "*i, globalToLocalMap[*i], sendEntryCount = " << *i << " " << A.globalToLocalMap[*i] << " " << sendEntryCount << endl;
      elementsToSend[sendEntryCount] = A.globalToLocalMap[*i]; // store local ids of entry to send
    }
  }

  // Convert matrix indices to local IDs
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< localNumberOfRows; i++) {
    for (int j=0; j<nonzerosInRow[i]; j++) {
      global_int_t curIndex = mtxIndG[i][j];
      int rankIdOfColumnEntry = GetRankOfMatrixRow(A.geom, RowRange, curIndex);
      if (A.geom->rank==rankIdOfColumnEntry) { // My column index, so convert to local index
        mtxIndL[i][j] = A.globalToLocalMap[curIndex];
      } else { // If column index is not a row index, then it comes from another processor
        mtxIndL[i][j] = externalToLocalMap[curIndex];
      }
    }
  }

  // Store contents in our matrix struct
  A.numberOfExternalValues = externalToLocalMap.size();
  A.localNumberOfColumns = A.localNumberOfRows + A.numberOfExternalValues;
  A.numberOfSendNeighbors = sendList.size();
  A.totalToBeSent = totalToBeSent;
  A.elementsToSend = elementsToSend;
  A.neighbors = neighbors;
  A.receiveLength = receiveLength;
  A.sendLength = sendLength;
  A.sendBuffer = sendBuffer;

  std::cout << " For rank " << A.geom->rank << " of " << A.geom->size << ", number of neighbors = " << A.numberOfSendNeighbors << std::endl;
  
  for (int i = 0; i < A.numberOfSendNeighbors; i++) {
    std::cout << "     rank " << A.geom->rank << " neighbor " << neighbors[i] << " send/recv length = " << sendLength[i] << "/" << receiveLength[i] << std::endl;
    /*
    for (local_int_t j = 0; j<sendLength[i]; ++j)
      std::cout << "       rank " << A.geom->rank << " elementsToSend[" << j << "] = " << elementsToSend[j] << std::endl;
    */
  }

#endif
// ifdef HPCG_NO_MPI

  return;
}
