
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
 @file Geometry.hpp

 HPCG data structure for problem geometry
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

/*!
  This defines the type for integers that have local subdomain dimension.

  Define as "long long" when local problem dimension is > 2^31
*/
typedef int local_int_t;
//typedef long long local_int_t;

/*!
  This defines the type for integers that have global dimension

  Define as "long long" when global problem dimension is > 2^31
*/
//typedef int global_int_t;
typedef long long global_int_t;

// This macro should be defined if the global_int_t is not long long
// in order to stop complaints from non-C++11 compliant compilers.
//#define HPCG_NO_LONG_LONG

/*!
  This is a data structure to contain all processor geometry information
*/
struct Geometry_STRUCT {
  int size; //!< Number of MPI processes
  int rank; //!< This process' rank in the range [0 to size - 1]
  int numThreads; //!< This process' number of threads

  int irow;
  int irow_begin;
  int irow_end;
  
  int grow;

  local_int_t nx;   //!< Number of x-direction grid points for each local subdomain
  local_int_t ny;   //!< Number of y-direction grid points for each local subdomain
  local_int_t nz;   //!< Number of z-direction grid points for each local subdomain
  int npx;  //!< Number of processors in x-direction
  int npy;  //!< Number of processors in y-direction
  int npz;  //!< Number of processors in z-direction
  int pz; //!< partition ID of z-dimension process that starts the second region of nz values
  int npartz; //!< Number of partitions with varying nz values
  int * partz_ids; //!< Array of partition ids of processor in z-direction where new value of nz starts (valid values are 1 to npz)
  local_int_t * partz_nz; //!< Array of length npartz containing the nz values for each partition
  int ipx;  //!< Current rank's x location in the npx by npy by npz processor grid
  int ipy;  //!< Current rank's y location in the npx by npy by npz processor grid
  int ipz;  //!< Current rank's z location in the npx by npy by npz processor grid
  global_int_t gnx;  //!< Global number of x-direction grid points
  global_int_t gny;  //!< Global number of y-direction grid points
  global_int_t gnz;  //!< Global number of z-direction grid points
  global_int_t gix0;  //!< Base global x index for this rank in the npx by npy by npz processor grid
  global_int_t giy0;  //!< Base global y index for this rank in the npx by npy by npz processor grid
  global_int_t giz0;  //!< Base global z index for this rank in the npx by npy by npz processor grid

};
typedef struct Geometry_STRUCT Geometry;

/*!
  Returns the rank of the MPI process that is assigned the global row index
  given as the input argument.

  @param[in] geom  The description of the problem's geometry.
  @param[in] index The global row index

  @return Returns the MPI rank of the process assigned the row
*/
inline int ComputeRankOfMatrixRow(const Geometry & geom, global_int_t index) {
  int seg = geom.grow / geom.size;
  int par = (geom.grow % geom.size) * seg + (geom.grow % geom.size) - 1;
  int rank;
  if (index <= par)
    rank = index / (seg + 1);
  else
    rank = (index - geom.grow % geom.size) / seg;
  return rank;
}


/*!
 Destructor for geometry data.

 @param[inout] data the geometry data structure whose storage is deallocated
 */
inline void DeleteGeometry(Geometry & geom) {

  return;
}



#endif // GEOMETRY_HPP
