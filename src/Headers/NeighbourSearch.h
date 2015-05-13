//=================================================================================================
//  NeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching
//  data structures and algorithms.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _NEIGHBOUR_SEARCH_H_
#define _NEIGHBOUR_SEARCH_H_


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "MeshlessFV.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#endif
using namespace std;




//=================================================================================================
//  Class NeighbourSearch
/// \brief   NeighbourSearch class definition.
/// \details Class for creating the SPH neighbour search data structure, and for computing local
///          neighbour lists and calling SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    20/04/2015
//=================================================================================================
template <int ndim>
class NeighbourSearch
{
#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
#endif
 public:

  //-----------------------------------------------------------------------------------------------
  NeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                  SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    kernrange(kernrangeaux),
    kernrangesqd(kernrangeaux*kernrangeaux),
    box(boxaux),
    kernp(kernaux),
    kernfac(1.0),
    timing(timingaux) {};
  //NeighbourSearch() {};
  //NeighbourSearch(const NeighbourSearch &) {};
  virtual ~NeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int, const int,
                         const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int, const int,
                              const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *) {return -1;};
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *) {};
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *) {};
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const int, ParticleType<ndim> *);
  virtual void BuildGhostPrunedTree(const int, const DomainBox<ndim> &);
  virtual void BuildMpiGhostTree(bool, int, int, int, int, int,
                                 ParticleType<ndim> *, Hydrodynamics<ndim> *, FLOAT);
  virtual void CommunicatePrunedTrees(vector<int>&, int);
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Hydrodynamics<ndim> *hydro, int rank);
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int);
  virtual void InitialiseCellWorkCounters(void);
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                               Hydrodynamics<ndim> *);
  virtual void UpdateGravityExportList(int, int, int, ParticleType<ndim> *,
                                       Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UpdateHydroExportList(int, int, int, ParticleType<ndim> *,
                                     Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Hydrodynamics<ndim>* hydro, int rank);
#endif


  //-----------------------------------------------------------------------------------------------
  FLOAT kernrange;               ///< Kernel extent (in units of h)
  FLOAT kernrangesqd;            ///< Kernel extent (squared)


  //-----------------------------------------------------------------------------------------------
  CodeTiming *timing;                  ///< Pointer to code timing object
  DomainBox<ndim> *box;                ///< Pointer to simulation bounding box
  SmoothingKernel<ndim> *kernp;        ///< Pointer to SPH kernel object

  bool neibcheck;                      ///< Flag to verify neighbour lists
  FLOAT kernfac;                       ///< ..

};



//=================================================================================================
//  Class BruteForceSearch
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class BruteForceSearch : public virtual NeighbourSearch<ndim>
{
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  BruteForceSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                   SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~BruteForceSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int, const int,
                         const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int, const int,
                              const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *);
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *);
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *) {};
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const int, ParticleType<ndim> *) {};
  virtual void BuildGhostPrunedTree(const int, const DomainBox<ndim> &) {};
  virtual void BuildMpiGhostTree(bool, int, int, int, int, int,
                                 ParticleType<ndim> *, Hydrodynamics<ndim> *, FLOAT);
  virtual void CommunicatePrunedTrees(vector<int>&, int);
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Hydrodynamics<ndim> *hydro, int rank);
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int);
  virtual void InitialiseCellWorkCounters(void) {};
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                               Hydrodynamics<ndim> *);
  virtual void UpdateGravityExportList(int, int, int, ParticleType<ndim> *,
                                       Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UpdateHydroExportList(int, int, int, ParticleType<ndim> *,
                                     Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Hydrodynamics<ndim>* hydro, int rank);
#endif

};



//=================================================================================================
//  Class HydroTree
/// \brief   Class containing tree for efficient neighbour searching and gravity calculations.
/// \details ..
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class HydroTree : public virtual NeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using NeighbourSearch<ndim>::ids_active_particles;
  vector<int> N_imported_part_per_proc;
#endif
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::box;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  HydroTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
            DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  virtual ~HydroTree();


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int, const int,
                         const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int, const int,
                              const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *);
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *);
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *);
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const int, ParticleType<ndim> *);
  virtual void BuildGhostPrunedTree(const int, const DomainBox<ndim> &);
  virtual void BuildMpiGhostTree(bool, int, int, int, int, int,
                                 ParticleType<ndim> *, Hydrodynamics<ndim> *, FLOAT);
  virtual void CommunicatePrunedTrees(vector<int>&, int);
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Hydrodynamics<ndim> *hydro, int rank);
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int);
  virtual void InitialiseCellWorkCounters(void) {};
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                               Hydrodynamics<ndim> *);
  virtual void UpdateGravityExportList(int, int, int, ParticleType<ndim> *,
                                       Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UpdateHydroExportList(int, int, int, ParticleType<ndim> *,
                                     Hydrodynamics<ndim> *, Nbody<ndim> *);
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Hydrodynamics<ndim>* hydro, int rank);
#endif
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, int, int *, ParticleType<ndim> *, string);
#endif


  // Additional functions for binary tree neighbour search
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(const int);
  void DeallocateMemory(void);
  void ComputeCellMonopoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeCellQuadrupoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeFastMonopoleForces(int, int, TreeCell<ndim> *,
                                 TreeCell<ndim> &, ParticleType<ndim> *);


  // Additional variables for grid
  //-----------------------------------------------------------------------------------------------
  const int Nleafmax;                              ///< Max. number of particles per leaf cell
  const int Nmpi;                                  ///< No. of MPI processes
  const FLOAT thetamaxsqd;                         ///< Geometric opening angle squared
  const FLOAT invthetamaxsqd;                      ///< 1 / thetamaxsqd
  const FLOAT macerror;                            ///< Error tolerance for gravity tree-MAC
  const string gravity_mac;                        ///< Multipole-acceptance criteria for tree
  const string multipole;                          ///< Multipole-order for cell gravity

  // Class variables
  //-----------------------------------------------------------------------------------------------
  int Ncell;                                       ///< Current no. of grid cells
  int Ncellmax;                                    ///< Max. allowed no. of grid cells
  int Nlistmax;                                    ///< Max. length of neighbour list
  int Ntot;                                        ///< No. of current points in list
  int Ntotold;                                     ///< Prev. no. of particles
  int Ntotmax;                                     ///< Max. no. of points in list
  int Ntotmaxold;                                  ///< Old value of Ntotmax
  FLOAT hmax;                                      ///< Store hmax in the tree
  FLOAT theta;                                     ///< Geometric opening angle
  Tree<ndim,ParticleType,TreeCell> *tree;          ///< Pointer to main (local) tree
  Tree<ndim,ParticleType,TreeCell> *ghosttree;     ///< Pointer to tree containing ghosts
                                                   ///< on local domain

  bool allocated_buffer;                           ///< ..
  int Nthreads;                                    ///< ..
  int *Nneibmaxbuf;                                ///< ..
  int *Ngravcellmaxbuf;                            ///< ..
  int **activelistbuf;                             ///< ..
  int **levelneibbuf;                              ///< ..
  ParticleType<ndim> **neibpartbuf;                ///< Local copy of neighbouring ptcls
  ParticleType<ndim> **activepartbuf;              ///< Local copy of SPH particle
  TreeCell<ndim> **cellbuf;                        ///< ..

#ifdef MPI_PARALLEL
  static const int Nghostprunedmax = 27;           ///< ..
  int Nprunedcellmax;                              ///< ..
  int *Nghostpruned;                               ///< ..
  int *Ncellexport;                                ///< ..
  int *Npartexport;                                ///< ..
  TreeCell<ndim> ***cellexportlist;                ///< List of cells
  Tree<ndim,ParticleType,TreeCell> *mpighosttree;  ///< Pointer to tree containing
                                                   ///< ghosts from other MPI procs.
  Tree<ndim,ParticleType,TreeCell> **prunedtree;   ///< 'Pruned' tree for MPI nodes.
                                                   ///< i.e. only uses top levels
  Tree<ndim,ParticleType,TreeCell> ***ghostprunedtree; ///< Tree of periodic ghost cells created
                                                         ///< from pruned trees
#endif

};
#endif