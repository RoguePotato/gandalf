//=============================================================================
//  EnergyLeapfrogDKD.cpp
//  Contains functions for energy equation integration using the Leapfrog
//  drift-kick-drift integration scheme.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "EnergyEquation.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  EnergyLeapfrogDKD::EnergyLeapfrogDKD()
/// EnergyLeapfrogDKD class constructor
//=============================================================================
template <int ndim>
EnergyLeapfrogDKD<ndim>::EnergyLeapfrogDKD(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



//=============================================================================
//  EnergyLeapfrogDKD::~EnergyLeapfrogDKD()
/// EnergyLeapfrogDKD class destructor
//=============================================================================
template <int ndim>
EnergyLeapfrogDKD<ndim>::~EnergyLeapfrogDKD()
{
}



//=============================================================================
//  EnergyLeapfrogDKD::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to 
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
//=============================================================================
template <int ndim>
void EnergyLeapfrogDKD<ndim>::EnergyIntegration
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,  ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[EnergyLeapfrogDKD::EnergyIntegration]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,nstep)	\
     shared(n,Nsph,timestep,sphdata)
  for (i=0; i<Nsph; i++) {
    SphParticle<ndim> part = sphdata[i];
    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    part = sphdata[i];
    part.u = sphdata[i].u0 + part.dudt*dt;
  }
  //---------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  EnergyLeapfrogDKD::CorrectionTerms
/// Empty function.  No correction terms for Leapfrog DKD scheme
//=============================================================================
template <int ndim>
void EnergyLeapfrogDKD<ndim>::EnergyCorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,  ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  return;
}



//=============================================================================
//  EnergyLeapfrogDKD::EndTimestep
/// Record all important thermal quantities at the end of the step for the 
/// start of the new timestep.
//=============================================================================
template <int ndim>
void EnergyLeapfrogDKD<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 FLOAT timestep,                    ///< [in] Base timestep value
 SphParticle<ndim> *sphdata)  ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size

  debug2("[EnergyLeapfrogDKD::EndTimestep]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,nstep)	\
  shared(n,Nsph,sphdata,timestep)
  for (i=0; i<Nsph; i++) {
    SphParticle<ndim>& part = sphdata[i];
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (dn == nstep) {
      sphdata[i].u0 = part.u;
      sphdata[i].dudt0 = part.dudt;
    }
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  EnergyLeapfrogDKD::Timestep
/// Compute explicit timestep such that u cannot change by a large fraction 
/// in one step, i.e. dt = const*u/|dudt + epsilon| 
/// where epsilon is to prevent the denominator becoming zero.
//=============================================================================
template <int ndim>
DOUBLE EnergyLeapfrogDKD<ndim>::Timestep
(SphParticle<ndim> &part)           ///< [inout] SPH particle reference
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}



template class EnergyLeapfrogDKD<1>;
template class EnergyLeapfrogDKD<2>;
template class EnergyLeapfrogDKD<3>;


