//=============================================================================
//  KDRadiationTree.cpp
//  ...
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
//=============================================================================



#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include "KDRadiationTree.h"
#include "Precision.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=============================================================================
//  KDRadiationTree::KDRadiationTree()
/// Constructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
KDRadiationTree<ndim,ParticleType>::KDRadiationTree(int Nleafmaxaux)
{
  allocated_tree = false;
  ltot           = 0;
  Ntot           = 0;
  Ntotmax        = 0;
  Ntotmaxold     = 0;
  Nleafmax       = Nleafmaxaux;
#if defined _OPENMP
  Nthreads       = omp_get_max_threads();
#else
  Nthreads       = 1;
#endif
}



//=============================================================================
//  KDRadiationTree::~KDRadiationTree()
/// Destructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
KDRadiationTree<ndim,ParticleType>::~KDRadiationTree()
{
}



//=============================================================================
//  KDRadiationTree::AllocateTreeMemory
/// Allocate memory for KD-tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::AllocateMemory(void)
{
  int ithread;                      // Thread id number

  debug2("[KDRadiationTree::AllocateMemory]");

  if (!allocated_tree || Ntotmax > Ntotmaxold) {
    if (allocated_tree) DeallocateMemory();
    Ntotmax = max(Ntotmax,Ntot);

    ids = new int[Ntotmax];
    inext = new int[Ntotmax];
    radcell = new struct KDRadTreeCell<ndim>[Ncellmax];

    allocated_tree = true;
  }

  return;
}



//=============================================================================
//  KDRadiationTree::DeallocateMemory
/// Deallocates all KD-tree memory
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::DeallocateMemory(void)
{
  debug2("[KDRadiationTree::DeallocateMemory]");

  if (allocated_tree) {
    delete[] radcell;
    delete[] inext;
    delete[] ids;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  KDRadiationTree::BuildTree
/// Call all routines to build/re-build the KD-tree on the local node.
/// If OpenMP is activated, the local domain is partitioned into sub-trees 
/// in order to improve the scalability of building and stocking the tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::BuildTree
(int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 ParticleType<ndim> *partdata)      ///< Particle data array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  FLOAT bbmin[ndim];
  FLOAT bbmax[ndim];

  debug2("[KDRadiationTree::BuildTree]");
  //timing->StartTimingSection("BUILD_TREE",2);

  cout << "Building tree with " << Npart << " particles" << endl;

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  // Set no. of tree members to total number of SPH particles (inc. ghosts)
  ltot_old = ltot;
  Ntotold    = Ntot;
  Ntot       = Npart;
  Ntotmaxold = Ntotmax;
  Ntotmax    = max(Ntot,Ntotmax);
  
  // Compute the size of all tree-related arrays now we know number of points
  ComputeTreeSize();
  
  // Allocate (or reallocate if needed) all tree memory
  AllocateMemory();
  
  // If the number of levels in the tree has changed (due to destruction or  
  // creation of new particles) then re-create tree data structure 
  // including linked lists and cell pointers
  if (ltot != ltot_old) CreateTreeStructure();

  // Create bounding box of SPH particles
  for (k=0; k<ndim; k++) bbmin[k] = big_number;
  for (k=0; k<ndim; k++) bbmax[k] = -big_number;
  for (i=0; i<Ntot; i++) {
    for (k=0; k<ndim; k++) {
      if (partdata[i].r[k] + 2.0*partdata[i].h > bbmax[k])
	bbmax[k] = partdata[i].r[k] + 2.0*partdata[i].h;
      if (partdata[i].r[k] - 2.0*partdata[i].h < bbmin[k])
	bbmin[k] = partdata[i].r[k] - 2.0*partdata[i].h;	   
    }
  }
  

  // Set properties for root cell before constructing tree
  radcell[0].N = Ntot;
  radcell[0].ifirst = 0;
  radcell[0].ilast = Ntot - 1;
  for (k=0; k<ndim; k++) radcell[0].bbmin[k] = bbmin[k]; //-big_number;
  for (k=0; k<ndim; k++) radcell[0].bbmax[k] = bbmax[k]; //big_number;
  for (k=0; k<ndim; k++) radcell[0].cexit[0][k] = -1;
  for (k=0; k<ndim; k++) radcell[0].cexit[1][k] = -1;
  for (i=0; i<Ntot; i++) inext[i] = -1;
  
  // If number of particles remains unchanged, use old id list 
  // (nearly sorted list should be faster for quick select).
  if (Ntot != Ntotold)
    for (i=0; i<Ntot; i++) ids[i] = i;    
  
  // Recursively build tree from root node down
  DivideTreeCell(0,Ntot-1,partdata,radcell[0]);

  // Calculate more optimal cell quantities for speeding up ray walking on tree
  OptimiseTree();
  
  return;
}



//=============================================================================
//  KDRadiationTree::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of 
/// the KD tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::ComputeTreeSize(void)
{
  debug2("[KDRadiationTree::ComputeTreeSize]");

  // Calculate maximum level of tree that can contain max. no. of particles
  lmax = 0;
  while (Nleafmax*pow(2,lmax) < Ntotmax) {
    lmax++;
  };
  gmax = pow(2,lmax);
  Ncellmax = 2*gmax - 1;


  // Calculate level of tree that can contain all current particles
  ltot = 0;
  while (Nleafmax*pow(2,ltot) < Ntot) {
    ltot++;
  };
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;

  cout << "Nleafmax : " << Nleafmax << endl;
  cout << "No. of ptcls in tree  : " << Ntot << "   " << Ntotmax << endl;
  cout << "No. of grid-cells     : " << gtot << "   " << gmax << endl;
  cout << "No. of levels on tree : " << ltot << "   " << lmax << endl;
  cout << "No. of cells in tree  : " << Ncell << "   " << Ncellmax << endl; 

  return;
}



//=============================================================================
//  KDRadiationTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::CreateTreeStructure(void)
{
  int c;                            // Dummy id of tree-level, then tree-cell
  int g;                            // Dummy id of grid-cell
  int l;                            // Dummy id of level
  int *c2L;                         // Increment to second child-cell
  int *cNL;                         // Increment to next cell if cell unopened

  debug2("[KDRadiationTree::CreateTreeStructure]");

  // Allocate memory for local arrays
  c2L = new int[ltot + 1];
  cNL = new int[ltot + 1];

  // Set pointers to second child-cell (if opened) and next cell (if unopened)
  for (l=0; l<ltot; l++) {
    c2L[l] = pow(2,ltot - l);
    cNL[l] = 2*c2L[l] - 1;
  }

  // Zero tree cell variables
  for (c=0; c<Ncell; c++) {
    radcell[c].id = c;
    radcell[c].c1 = -1;
    radcell[c].c2 = -1;
    radcell[c].ifirst = -1;
    radcell[c].ilast = -1;
    radcell[c].N = 0;
  }
  g = 0;
  radcell[0].level = 0;

  // Loop over all cells and set all other pointers
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (radcell[c].level == ltot) {                    // If on leaf level
      radcell[c].cnext = c + 1;                        // id of next cell
    }
    else {
      radcell[c+1].level = radcell[c].level + 1;          // Level of 1st child
      radcell[c].c1 = c + 1;
      radcell[c].c2 = c + c2L[radcell[c].level];          // id of 2nd child
      radcell[radcell[c].c2].level = radcell[c].level + 1; // Level of 2nd child
      radcell[c].cnext = c + cNL[radcell[c].level];       // Next cell id
    }

  }
  //---------------------------------------------------------------------------


  // Free locally allocated memory
  delete[] cNL;
  delete[] c2L;

  return;
}



//=============================================================================
//  KDRadiationTree::DivideTreeCell
/// Recursive routine to divide a tree cell into two children cells.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::DivideTreeCell
(int ifirst,                        ///< [in] Aux. id of first particle in cell
 int ilast,                         ///< [in] Aux. id of last particle in cell
 ParticleType<ndim> *partdata,      ///< [in] Pointer to main SPH object
 KDRadTreeCell<ndim> &cell)         ///< [inout] Cell to be divided
{
  int i;                            // Aux. child cell counter
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int k_divide = 0;                 // Division dimension
  FLOAT rkmax = 0.0;                // Max. box size of all dimensions
  FLOAT rdivide;                    // Coordinate value at division


  // If cell is a leaf cell, do not divide further and set linked lists
  if (cell.level == ltot) {
    if (cell.N > 0) {
      for (j=cell.ifirst; j<cell.ilast; j++) inext[ids[j]] = ids[j+1];
      inext[ids[cell.ilast]] = -1;
      cell.ifirst = ids[cell.ifirst];
      cell.ilast = ids[cell.ilast];
    }
    else {
      cell.ifirst = -1;
      cell.ilast = -1;
    }
    StockCellProperties(cell,partdata);
    return;
  }

  // Determine dimension to split the cell along.
  // For now, simply split along direction of the bounding box's longest axis
  for (k=0; k<ndim; k++) {
    if (cell.bbmax[k] - cell.bbmin[k] > rkmax) {
      rkmax = cell.bbmax[k] - cell.bbmin[k];
      k_divide = k;
    }
  }
  cell.k_divide = k_divide;


  // Find median value along selected division plane and re-order array 
  // so particles reside on correct side of division
  rdivide = QuickSelect(cell.ifirst,cell.ilast,
			cell.ifirst+cell.N/2,k_divide,partdata);

  // Set properties of first child cell
  for (k=0; k<ndim; k++) radcell[cell.c1].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) radcell[cell.c1].bbmax[k] = cell.bbmax[k];
  for (k=0; k<ndim; k++) radcell[cell.c1].cexit[0][k] = cell.cexit[0][k];
  for (k=0; k<ndim; k++) radcell[cell.c1].cexit[1][k] = cell.cexit[1][k];
  radcell[cell.c1].N = cell.N/2;
  if (radcell[cell.c1].N != 0) {
    radcell[cell.c1].ifirst = ifirst;
    radcell[cell.c1].ilast = ifirst + cell.N/2 - 1;
  }

  // Set properties of second child cell
  for (k=0; k<ndim; k++) radcell[cell.c2].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) radcell[cell.c2].bbmax[k] = cell.bbmax[k];
  for (k=0; k<ndim; k++) radcell[cell.c2].cexit[0][k] = cell.cexit[0][k];
  for (k=0; k<ndim; k++) radcell[cell.c2].cexit[1][k] = cell.cexit[1][k];
  radcell[cell.c2].N = cell.N - radcell[cell.c1].N;
  if (radcell[cell.c2].N != 0) {
    radcell[cell.c2].ifirst = ifirst + cell.N/2;
    radcell[cell.c2].ilast = ilast;
  }
  assert(cell.N == radcell[cell.c1].N + radcell[cell.c2].N);

  /*if (radcell[cell.c1].N < radcell[cell.c2].N) {
    cout << "WTF?? : " << radcell[cell.c1].N << "   " << radcell[cell.c2].N << endl;
    cout << "p1 : " << radcell[cell.c1].ilast - radcell[cell.c1].ifirst + 1 << endl;
    cout << "p2 : " << radcell[cell.c2].ilast - radcell[cell.c2].ifirst + 1 << endl;
    exit(0);
    }*/

  // Set new cell boundaries depending on number of particles in cells
  if (radcell[cell.c1].N > 0 && radcell[cell.c2].N > 0) {
    radcell[cell.c1].bbmax[k_divide] = rdivide;
    radcell[cell.c2].bbmin[k_divide] = rdivide;
    radcell[cell.c1].cexit[1][k_divide] = cell.c2;
    radcell[cell.c2].cexit[0][k_divide] = cell.c1;
  }
  else if (radcell[cell.c2].N > 0) {
    radcell[cell.c1].bbmax[k_divide] = -big_number;
  }


  // Now divide the new child cells as a recursive function
#if defined _OPENMP
  if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel default(none) private(i) \
  shared(cell,ifirst,ilast,partdata) num_threads(2)
    {
#pragma omp for 
      for (i=0; i<2; i++) {
	if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,
                                   partdata,radcell[cell.c1]);
	else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,
                                        partdata,radcell[cell.c2]);
      }
#pragma omp barrier
    }
  }
  else {
    for (i=0; i<2; i++) {
      if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,
                                 partdata,radcell[cell.c1]);
      else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,
                                      partdata,radcell[cell.c2]);
    }
  }
#else
  for (i=0; i<2; i++) {
    if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,
                               partdata,radcell[cell.c1]);
    else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,
                                    partdata,radcell[cell.c2]);
  }
#endif

  if (cell.N != radcell[cell.c1].N + radcell[cell.c2].N) {
    cout << "Checking : " << cell.N << "   " << radcell[cell.c1].N 
	 << "    " << radcell[cell.c2].N << endl;
  }
  assert(cell.N == radcell[cell.c1].N + radcell[cell.c2].N);

  // Stock all cell properties once constructed
  StockCellProperties(cell,partdata);

  return;
}



//=============================================================================
//  KDRadiationTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct 
/// side of the division.  Uses the QuickSelect algorithm.
//=============================================================================
template <int ndim, template<int> class ParticleType>
FLOAT KDRadiationTree<ndim,ParticleType>::QuickSelect
(int left,                          ///< Left-most id of particle in array
 int right,                         ///< Right-most id of particle in array
 int jpivot,                        ///< Pivot/median point
 int k,                             ///< Dimension of sort
 ParticleType<ndim> *partdata)      ///< Pointer to main SPH object
{
  int i;                            // ..
  int j;                            // Aux.
  int jguess;                       // ..
  int jtemp;                        // ..
  FLOAT rpivot;                     // Position pivot for quick-select


  // Place all particles left or right of chosen pivot point.
  // Iterate until correct median pivot has been identified.
  //---------------------------------------------------------------------------
  do {

    // Make a guess of pivot value
    jguess = (left + right)/2;
    rpivot = partdata[ids[jguess]].r[k];

    // ..
    jtemp = ids[jguess];
    ids[jguess] = ids[right];
    ids[right] = jtemp;

    // ..
    jguess = left;


    //-------------------------------------------------------------------------
    for (j=left; j<right; j++) {
      assert(j < Ntot);
      if (partdata[ids[j]].r[k] < rpivot) {
	jtemp = ids[j];
	ids[j] = ids[jguess];
	ids[jguess] = jtemp;
        jguess++;
      }

    }
    //-------------------------------------------------------------------------


    // Move ?? particle from end of array to index i
    jtemp = ids[right];
    ids[right] = ids[jguess];
    ids[jguess] = jtemp;

    // jguess is lower than jpivot.  
    // Only need to search between jguess+1 and right
    if (jguess < jpivot) left = jguess + 1;

    // jguess is higher than jpivot.  
    // Only need to search between left and jguess-1
    else if (jguess > jpivot) right = jguess - 1;

  } while (jguess != jpivot);
  //---------------------------------------------------------------------------


  return rpivot;
}



//=============================================================================
//  KDRadiationTree::StockTree
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::StockTree
(KDRadTreeCell<ndim> &cell,         ///< Reference to cell to be stocked
 ParticleType<ndim> *partdata)      ///< SPH particle data array
{
  int i;                            // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.level != ltot) {
    if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel for default(none) private(i) \
  shared(cell,partdata) num_threads(2)
      for (i=0; i<2; i++) {
	if (i == 0) StockTree(radcell[cell.c1],partdata);
	else if (i == 1) StockTree(radcell[cell.c2],partdata);
      }
    }
    else {
      for (i=0; i<2; i++) {
	if (i == 0) StockTree(radcell[cell.c1],partdata);
	else if (i == 1) StockTree(radcell[cell.c2],partdata);
      }
    }
  }

  // Stock node once all children are stocked
  StockCellProperties(cell,partdata);

  return;
}



//=============================================================================
//  KDRadiationTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::StockCellProperties
(KDRadTreeCell<ndim> &cell,         ///< Reference to current tree cell
 ParticleType<ndim> *partdata)      ///< Particle data array
{
  int cc,ccc;                       // Cell counters
  int i;                            // Particle counter
  int iaux;                         // Aux. particle i.d. variable
  int j;                            // ..
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT mi;                         // Mass of particle i
  FLOAT p = 0.0;                    // ..
  FLOAT lambda = 0.0;               // ..
  KDRadTreeCell<ndim> &child1 = radcell[cell.c1];
  KDRadTreeCell<ndim> &child2 = radcell[cell.c2];


  // Zero all summation variables for all cells
  cell.uniform = false;
  cell.N = 0;
  cell.m = 0.0;
  cell.rho = 0.0;
  cell.temp = 0.0;
  cell.uphoton = 0.0;
  cell.opacity = 0.01;
  for (k=0; k<ndim; k++) cell.r[k] = 0.0;
  for (k=0; k<ndim; k++) cell.v[k] = 0.0;
  for (k=0; k<ndim; k++) cell.rcell[k] = 0.5*(cell.bbmax[k] + cell.bbmin[k]);

  // Calculate cell volume
  cell.volume = 1.0;
  for (k=0; k<ndim; k++) cell.volume *= (cell.bbmax[k] - cell.bbmin[k]);


  // If this is a leaf cell, sum over all particles
  //---------------------------------------------------------------------------
  if (cell.level == ltot) {
    
    // Loop over all particles in cell summing their contributions
    i = cell.ifirst;
    while (i != -1) {
      if (partdata[i].itype != dead) {
	cell.N++;
	cell.m += partdata[i].m;
        cell.rho += partdata[i].m*partdata[i].rho;
	for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
	for (k=0; k<ndim; k++) cell.v[k] += partdata[i].m*partdata[i].v[k];
	for (k=0; k<ndim; k++) {
	  if (partdata[i].r[k] < cell.bbmin[k]) 
            cell.bbmin[k] = partdata[i].r[k];
	  if (partdata[i].r[k] > cell.bbmax[k])
            cell.bbmax[k] = partdata[i].r[k];
	}
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
    
    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      cell.rho /= cell.m;
    }
    

  }
  // For non-leaf cells, sum together two children cells
  //---------------------------------------------------------------------------
  else {
    
    cell.N = child1.N + child2.N;
    if (cell.N > 0) {
      cell.m = child1.m + child2.m;
      cell.rho = (child1.m*child1.rho + child2.m*child2.rho)/cell.m;
      for (k=0; k<ndim; k++) cell.r[k] =
	(child1.m*child1.r[k] + child2.m*child2.r[k])/cell.m;
      for (k=0; k<ndim; k++) cell.v[k] =
	(child1.m*child1.v[k] + child2.m*child2.v[k])/cell.m;
    }
        
  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  KDRadiationTree::OptimiseTree
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void KDRadiationTree<ndim,ParticleType>::OptimiseTree(void)
{
  int c;                            // ..
  int c2;                           // ..
  int cexit;                        // ..
  int k;                            // ..
  int k_divide;                     // ..
  int level;                        // ..
  KDRadTreeCell<ndim> *cell;        // ..


  debug2("[KDRadiationTree::OptimiseTree]");


  // Loop over all cells in the tree
  //===========================================================================
  for (c=0; c<Ncell; c++) {

    cell = &(radcell[c]);
    level = cell->level;


    // Loop over all dimensions in the tree
    //-------------------------------------------------------------------------
    for (k=0; k<ndim; k++) {


      // First do left-hand side of cell
      //-----------------------------------------------------------------------
      cexit = radcell[c].cexit[0][k];

      while (radcell[cexit].level < level && cexit != -1) {
	k_divide = radcell[cexit].k_divide;
        c2 = radcell[cexit].c2;
        
	// First check if divide is in same direction, then use child directly
        if (k_divide == k && radcell[c2].N > 0)
          cexit = c2;
        else if (radcell[cexit+1].N > 0 && 
		 radcell[cexit+1].bbmax[k_divide] >= cell->bbmax[k_divide] &&
		 radcell[cexit+1].bbmin[k_divide] <= cell->bbmin[k_divide])
	  cexit = cexit + 1;
        else if (radcell[c2].N > 0 && 
		 radcell[c2].bbmax[k_divide] >= cell->bbmax[k_divide] &&
		 radcell[c2].bbmin[k_divide] <= cell->bbmin[k_divide])
	  cexit = c2;
	else 
	  break;
      };

      /*if (radcell[c].cexit[0][k] != cexit) {
	cout << "Found optimal exit cell : " << c << "   " 
	     << radcell[c].cexit[0][k] << "   " << cexit 
	     << "  k : " << k << endl;
	cout << "bbmin of cell  : " << cell->bbmin[0] << "   " 
	     << cell->bbmin[1] << "   " << cell->bbmin[2] << endl;
	cout << "bbmax of cell  : " << cell->bbmax[0] << "   " 
	     << cell->bbmax[1] << "   " << cell->bbmax[2] << endl;
	cout << "bbmin of cexit : " << radcell[cexit].bbmin[0] << "   " 
	     << radcell[cexit].bbmin[1] << "   " 
	     << radcell[cexit].bbmin[2] << endl;
	cout << "bbmax of cexit : " << radcell[cexit].bbmax[0] << "   " 
	     << radcell[cexit].bbmax[1] << "   " 
	     << radcell[cexit].bbmax[2] << endl;
	     }*/

      radcell[c].cexit[0][k] = cexit;


      // Now do right-hand side of cell
      //-----------------------------------------------------------------------
      cexit = radcell[c].cexit[1][k];

      // Loop down levels to find lowest cell that can be used for exit face
      while (radcell[cexit].level < level && cexit != -1) {
	k_divide = radcell[cexit].k_divide;
        c2 = radcell[cexit].c2;
        
	// First check if divide is in same direction, then use child directly
        if (k_divide == k && radcell[cexit+1].N > 0)
          cexit = cexit + 1;
        else if (k_divide == k && radcell[c2].N > 0)
          cexit = c2;
        else if (radcell[cexit+1].N > 0 && 
	 radcell[cexit+1].bbmax[k_divide] >= cell->bbmax[k_divide] &&
	 radcell[cexit+1].bbmin[k_divide] <= cell->bbmin[k_divide])
	  cexit = cexit + 1;
        else if (radcell[c2].N > 0 && 
	 radcell[c2].bbmax[k_divide] >= cell->bbmax[k_divide] &&
	 radcell[c2].bbmin[k_divide] <= cell->bbmin[k_divide])
	  cexit = c2;
	else 
	  break;
      };

      radcell[c].cexit[1][k] = cexit;


    }
    //-------------------------------------------------------------------------


  }
  //===========================================================================


  return;
}



//=============================================================================
//  KDRadiationTree::FindCell
/// Find i.d. of cell adjacent to current cell that the radiation packet is 
/// travelling into.
//=============================================================================
template <int ndim, template<int> class ParticleType>
int KDRadiationTree<ndim,ParticleType>::FindCell
(int cparent,                       ///< [in] i.d. of larger parent cell
 FLOAT rp[ndim])                    ///< [in] Position of point/ray
{
  int c = cparent;                  // Cell i.d.
  int c1;                           // i.d. of 1st cell child
  int k_divide;                     // Dimension of cell division

  // Walk back down through tree to bottom level
  //---------------------------------------------------------------------------
  while (radcell[c].level < ltot) {

#ifdef OUTPUT_ALL
    cout << "Searching for cell : " << radcell[c].level << "   " << ltot << endl;
#endif

    c1 = c + 1;
    k_divide = radcell[c].k_divide;

    // If point is left of divide, pick 1st child cell.  Else pick 2nd child.
    if (rp[k_divide] < radcell[c1].bbmax[k_divide])
      c = c1;
    else
      c = radcell[c].c2;

  };
  //---------------------------------------------------------------------------

  return c;
}




template class KDRadiationTree<1,GradhSphParticle>;
template class KDRadiationTree<2,GradhSphParticle>;
template class KDRadiationTree<3,GradhSphParticle>;
template class KDRadiationTree<1,SM2012SphParticle>;
template class KDRadiationTree<2,SM2012SphParticle>;
template class KDRadiationTree<3,SM2012SphParticle>;
template class KDRadiationTree<1,GodunovSphParticle>;
template class KDRadiationTree<2,GodunovSphParticle>;
template class KDRadiationTree<3,GodunovSphParticle>;
