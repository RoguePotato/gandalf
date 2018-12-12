//=================================================================================================
//  Sph.cpp
//  Contains important default routines for Sph class.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Parameters.h"
#include "EOS.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  Sph::Sph
/// Constructor for parent SPH class.  Initialises important variables and
/// sets important parameters using initialialisation lists.
//=================================================================================================
template <int ndim>
Sph<ndim>::Sph(int _hydro_forces, int _self_gravity, FLOAT _alpha_visc, FLOAT _beta_visc,
               FLOAT _h_fac, FLOAT _h_converge, aviscenum _avisc, acondenum _acond,
               tdaviscenum _tdavisc, string _gas_eos, string _KernelName, int _size_sph,
               SimUnits &units, Parameters *params):
  Hydrodynamics<ndim>(_hydro_forces, _self_gravity, _h_fac, _gas_eos,
                      _KernelName, _size_sph, units, params),
  size_sph_part(_size_sph),
  acond(_acond),
  avisc(_avisc),
  tdavisc(_tdavisc),
  alpha_visc(_alpha_visc),
  beta_visc(_beta_visc),
  h_converge(_h_converge),
  //create_sinks(0),
  fixed_sink_mass(0)
{
  Ngather = 0;
  hmin_sink = big_number;
  conservative_sph_star_gravity = params->intparams["conservative_sph_star_gravity"];

  // Set clump parameters.
  Nclump = 0;
  rho_unit = units.rho.outscale * units.rho.outcgs;
  clump_dens_step = params->floatparams["clump_dens_step"];
  clump_dens_min = pow10(params->floatparams["clump_dens_min"]) / rho_unit;
  clump_dens_max = pow10(params->floatparams["clump_dens_max"]) / rho_unit;
  clump_min_dist = params->floatparams["clump_min_dist"] / units.r.outscale;
  clump_min_star_dist = params->floatparams["clump_min_star_dist"] / units.r.outscale;
}



//=================================================================================================
//  Sph::InitialSmoothingLengthGuess
/// Perform initial guess of smoothing.  In the abscence of more sophisticated techniques, we guess
/// the smoothing length assuming a uniform density medium with the same volume and total mass.
//=================================================================================================
template <int ndim>
void Sph<ndim>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[ndim];                // Min. extent of bounding box
  FLOAT rmax[ndim];                // Max. extent of bounding box

  debug2("[Sph::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = (int) (2.0*kernp->kernrange*h_fac);
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro);
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro));
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = (int) (4.0*pi*pow(kernp->kernrange*h_fac,3)/3.0);
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0*volume*(FLOAT) Ngather)/(32.0*pi*(FLOAT) Nhydro),onethird);
  }
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    SphParticle<ndim>& part = GetSphParticlePointer(i);
    part.h = h_guess;
    part.hrangesqd = kernp->kernrangesqd*part.h*part.h;
  }

  return;
}

//=================================================================================================
//  Sph::ZeroAccelerations
/// Initialise key variables before force calculations
//=================================================================================================
template <int ndim>
void Sph<ndim>::ZeroAccelerations()
{
  for (int i=0; i< Nhydro; i++) {
    SphParticle<ndim>& part = GetSphParticlePointer(i);
    if (part.flags.check(active)) {
      part.levelneib = 0;
      part.div_v     = (FLOAT) 0.0;
      part.dudt      = (FLOAT) 0.0;
      part.gpot      = (FLOAT) 0.0;
      part.gpot_hydro= (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) part.atree[k] = (FLOAT) 0.0;
    }
  }
}

//=================================================================================================
//  Sph::ClumpFind
/// Performs the clump finding algorithm.
//=================================================================================================
template <int ndim>
void Sph<ndim>::ClumpFind(Nbody<ndim> *nbody) {
  assert(Nclump < 1024);

  for (int i = 0; i < Nhydro; ++i) {
    SphParticle<ndim> &part = GetSphParticlePointer(i);
    bool too_close = false;


    // Perform the density threshold check first. Reset it from being a clump if
    // the max density threshold is reached.
    if (part.rho < clump_dens_min || part.rho > clump_dens_max) {
      part.clump = 1E30;
      continue;
    }
    // Next perform a check see if the particle is too close to the star
    // (where densities may be a little high).
    if (Distance(part.r, nbody->stardata[0].r, 2) < clump_min_star_dist) {
      continue;
    }

    // Check if near other clumps, provided there exists at least one clump.
    for (int j = 0; j < Nclump && Nclump > 0; ++j) {
      SphParticle<ndim> &clump = GetSphParticlePointer(clumps[j]);

      // Don't let clumps be set when close to others.
      // FLOAT dist = max(clump_min_dist, 5.0 * part.h);
      if (Distance(part.r, clump.r, 2) < clump_min_dist || i == clumps[j]) {
        too_close = true;
        break;
      }
    }

    // All tests passed to label the particle a clump. We must also update the
    // clump array and the number of the clumps. Typically clumps should not
    // be forming ofter, therefore the clump_flag and clump_dens should be
    // find to update just once per timestep.
    if (!too_close) {
      clumps[Nclump++] = i;
      part.clump = part.rho;
    }
  }

  // Reset the clump flag.
  clump_flag = -1;

  // Subsequent density checks for output flag.
  for (int i = 0; i < Nclump; ++i) {
    SphParticle<ndim> &part = GetSphParticlePointer(clumps[i]);

    // FIXME: This stepping needs to be corrected. Probably want to round the
    // value to nearest density step.
    if (log10(part.rho * rho_unit) > log10(part.clump * rho_unit) + clump_dens_step) {
      part.clump = part.rho;

      // TODO: Could always pass through i as a clump_flag? Would require more
      // recording in the restart file. Output in restart file the clump
      // particle ID and its current number.
      clump_dens = part.clump;
      clump_flag = part.iorig;
      break;
    }
    // TODO: Else if the density down by some amount, stop it being a clump?
  } 
}

template class Sph<1>;
template class Sph<2>;
template class Sph<3>;
