#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <iostream>
#include <cstdlib>
#ifdef HPCG_DETAILED_DEBUG
using std::cin;
#endif
using std::endl;

#include <vector>

#include "hpcg.hpp"

#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "CheckProblem.hpp"
#include "ExchangeHalo.hpp"
#include "OptimizeProblem.hpp"
#include "WriteProblem.hpp"
#include "ReportResults.hpp"
#include "mytimer.hpp"
#include "ComputeSPMV_ref.hpp"
#include "ComputeMG_ref.hpp"
#include "ComputeResidual.hpp"
#include "CG.hpp"
#include "CG_ref.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "CGData.hpp"
#include "TestCG.hpp"
#include "TestSymmetry.hpp"
#include "TestNorms.hpp"

int main(int argc, char * argv[]) {

#ifndef HPCG_NO_MPI
  MPI_Init(&argc, &argv);
#endif

  HPCG_Params params;

  HPCG_Init(&argc, &argv, params);

  int size = params.comm_size, rank = params.comm_rank; // Number of MPI processes, My process ID

#ifdef HPCG_DETAILED_DEBUG
  if (size < 100 && rank==0) HPCG_fout << "Process "<<rank<<" of "<<size<<" is alive with " << params.numThreads << " threads." <<endl;

  if (rank==0) {
    char c;
    std::cout << "Press key to continue"<< std::endl;
    std::cin.get(c);
  }
#ifndef HPCG_NO_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  int ierr = 0;  // Used to check return codes on function calls

  /////////////////////////
  // Problem setup Phase //
  /////////////////////////

#ifdef HPCG_DEBUG
  double t1 = mytimer();
#endif

  //load data from file
  std::fstream data_stream(argv[1], std::ios_base::in);
  int net_cnt;
  data_stream>>net_cnt;
  int node_cnt;
  data_stream>>node_cnt;
  
  // Construct the geometry and linear system
  Geometry * geom = new Geometry;
  GenerateGeometry(size, rank, params.numThreads, node_cnt, geom);
  //  std::cout<<rank<<": irow="<<geom->irow<< " begin="<<geom->irow_begin<<" end="<<geom->irow_end<<std::endl;

  std::map<int, std::pair<int, std::pair<int, int> > > id_to_node;
  std::map<std::pair<int, std::pair<int,int> >, int> node_to_id;
  std::map<int, std::vector<std::pair<int, double> > > id_to_neighbors;
  std::map<int, double> id_to_b;
  std::map<int, double> id_to_voltage;
  for(int rem_node=node_cnt; rem_node>0; rem_node--){
    int id;
    data_stream>>id;
    int net, x, y;
    data_stream>>net>>x>>y;
    id_to_node[id]=std::make_pair(net,std::make_pair(x,y));
    node_to_id[std::make_pair(net,std::make_pair(x,y))]=id;
    int neighbor_cnt;
    data_stream>>neighbor_cnt;
    for(int rem_neighbor=neighbor_cnt; rem_neighbor>0; rem_neighbor--){
      int neighbor_id;
      double rv;
      data_stream>>neighbor_id>>rv;
      id_to_neighbors[id].push_back(std::make_pair(neighbor_id, rv));
    }
    double b;
    data_stream>>b;
    id_to_b[id]=b;
    double voltage;
    data_stream>>voltage;
    id_to_voltage[id]=voltage;
  }

  std::cout<<"node_cnt="<<node_cnt<<std::endl;
  std::cout<<"net_cnt="<<net_cnt<<std::endl;

  // Use this array for collecting timing information
  std::vector< double > times(10,0.0);

  double setup_time = mytimer();

  SparseMatrix A;
  InitializeSparseMatrix(A, geom);

  Vector b, x;
  GenerateProblem(A, &b, &x, id_to_neighbors, id_to_b);
  SetupHalo(A);

  int numberOfMgLevels = 0; // Number of levels including first
  SparseMatrix * curLevelMatrix = &A;
  for (int level = 1; level< numberOfMgLevels; ++level) {
    GenerateCoarseProblem(*curLevelMatrix);
    curLevelMatrix = curLevelMatrix->Ac; // Make the just-constructed coarse grid the next level
  }

  setup_time = mytimer() - setup_time; // Capture total time of setup
  times[9] = setup_time; // Save it for reporting

  curLevelMatrix = &A;
  Vector * curb = &b;
  Vector * curx = &x;
  for (int level = 0; level< numberOfMgLevels; ++level) {
    //CheckProblem(*curLevelMatrix, curb, curx, curxexact);
     curLevelMatrix = curLevelMatrix->Ac; // Make the nextcoarse grid the next level
     curb = 0; // No vectors after the top level
     curx = 0;
  }

  CGData data;
  InitializeSparseCGData(A, data);

  // Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;

#ifdef HPCG_DEBUG
  t1 = mytimer();
#endif

  int niters = 0;
  int totalNiters_ref = 0;
  double normr = 0.0;
  double normr0 = 0.0;
  int refMaxIters = 50;

  ZeroVector(x); // Zero out x
  double cg_time = mytimer();
  ierr = CG( A, data, b, x, 5000, 1e-4, niters, normr, normr0, &times[0], true);
  cg_time = mytimer() - cg_time;
  if (rank == 0)
    std::cout<<"Process #"<<rank<<" run time = " <<cg_time<<" secs."<<std::endl;
  if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
  if (rank == 0){
    std::cout<<"norm="<<normr<<std::endl;
    std::cout<<"niter="<<niters<<std::endl;
  }
  
  std::vector<int> error_node_vector;
  int error_node_cnt=0;
  for(int i=0; i<geom->irow; i++){
    int node_id = geom->irow_begin + i + 1;
    double dis = id_to_voltage[node_id] - x.values[i];
    dis = (dis>0.0) ? dis : -dis;
    if (dis>1e-2){
      error_node_cnt++;
      error_node_vector.push_back(node_id);
    }
  } 

/*
  for(int i = 0; i < error_node_vector.size(); i++){
    int node_id = error_node_vector[i];
    std::cout << id_to_node[node_id].first << "," << id_to_node[node_id].second.first << "," << id_to_node[node_id].second.second  << " : " << id_to_voltage[node_id] <<" " << x.values[node_id-1]<<std::endl;
  }
*/

  std::cout << "num of node with incorrect voltage = " << error_node_cnt << "/" << geom->irow << std::endl;

  // Clean up
  //DeleteMatrix(A); // This delete will recursively delete all coarse grid data
  DeleteCGData(data);
  DeleteVector(x);
  DeleteVector(b);

  HPCG_Finalize();

  // Finish up
#ifndef HPCG_NO_MPI
  MPI_Finalize();
#endif
  return 0;
}
