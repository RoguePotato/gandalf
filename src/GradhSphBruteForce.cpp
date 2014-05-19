//=============================================================================
//  GradhSphBruteForce.cpp
//  Contains all routines for generating SPH neighbour lists using 
//  brute-force (i.e. direct summation over all particles).
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



#include <iostream>
#include <math.h>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "SphKernel.h"
#if defined MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=============================================================================
//  GradhSphBruteForce::GradhSphBruteForce
/// GradhSphBruteForce class constructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
GradhSphBruteForce<ndim,ParticleType>::GradhSphBruteForce
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  BruteForceSearch<ndim,ParticleType>(kernrangeaux,boxaux,kernaux,timingaux)
{
}



//=============================================================================
//  GradhSphBruteForce::~GradhSphBruteForce
/// GradhSphBruteForce class destructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
GradhSphBruteForce<ndim,ParticleType>::~GradhSphBruteForce()
{
}



//=============================================================================
//  GradhSphBruteForce::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphBruteForce<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                     // Particle and dimension counters
  int Nneib = 0;                    // No. of (non-dead) neighbours
  int okflag;                       // Flag valid smoothing length
  int *neiblist;                    // List of neighbours
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *drsqd;                     // Distance squared
  FLOAT *gpot;                      // Array of neib. grav. potentials
  FLOAT *m;                         // Array of neib. position vectors
  FLOAT *mu;                        // Array of neib. mass*u values
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphBruteForce::UpdateAllSphProperties]");

  // Store masses in separate array
  gpot = new FLOAT[Ntot];
  m = new FLOAT[Ntot];
  neiblist = new int[Ntot];
  for (i=0; i<Ntot; i++) {
    if (sphdata[i].itype == dead) continue;
    neiblist[Nneib] = i;
    gpot[Nneib] = sphdata[i].gpot;
    m[Nneib] = sphdata[i].m;
    Nneib++;
  }

  // Create parallel threads
  //===========================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,okflag,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nsph,Ntot,sph,sphdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //-------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nsph; i++) {

      // Skip over inactive particles
      if (!sphdata[i].active || sphdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle 
      // and all neighbours here
      //-----------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) { 
	j = neiblist[jj];
    	for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - rp[k];
    	drsqd[jj] = DotProduct(dr,dr,ndim);
      }
      //-----------------------------------------------------------------------

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,Nneib,big_number,m,mu,drsqd,
                             gpot,sphdata[i],nbody);
  
    }
    //-------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===========================================================================

  delete[] neiblist;
  delete[] m;
  delete[] gpot;

  return;
}



//=============================================================================
//  GradhSphBruteForce::UpdateAllSphHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphBruteForce<ndim,ParticleType>::UpdateAllSphHydroForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,k;                        // Particle and dimension counters
  int Nneib;                        // No. of neighbours
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Relative distance vector
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Gather kernel extent (squared)
  FLOAT hrangesqdj;                 // Scatter kernel extent (squared)
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *dr;                        // Array of neib. position vectors
  FLOAT *drmag;                     // Array of neib. distances
  FLOAT *invdrmag;                  // Array of neib. inverse distances
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphBruteForce::UpdateAllSphHydroForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      if (sphdata[j].itype == dead) continue;
      hrangesqdj = pow(kernfac*kernp->kernrange*sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i,Nneib,neiblist,drmag,invdrmag,dr,
			       sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  
  // Free all allocated memory
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



//=============================================================================
//  GradhSphBruteForce::UpdateAllSphForces
/// Empty function for now
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphBruteForce<ndim,ParticleType>::UpdateAllSphForces
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen,          ///< [in] ..
 Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
 Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphBruteForce::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*
      sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise
    // forces twice)
    Nneib = 0;
    for (j=0; j<Nsph; j++)
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=============================================================================
//  GradhSphBruteForce::UpdateAllSphGravForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphBruteForce<ndim,ParticleType>::UpdateAllSphGravForces
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen,          ///< [in] ..
 Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
 Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphBruteForce::UpdateAllSphGravForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];

  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*
      sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise
    // forces twice)
    Nneib = 0;
    for (j=0; j<Nsph; j++)
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



template class GradhSphBruteForce<1,GradhSphParticle>;
template class GradhSphBruteForce<2,GradhSphParticle>;
template class GradhSphBruteForce<3,GradhSphParticle>;
