// ============================================================================
// SphSnapshot.cpp
// ============================================================================

#include <ctime>
#include <cstdio>
#include <iostream>
#include "SphSnapshot.h"
#include "Sph.h"
#include "SphParticle.h"
#include "SphSimulation.h"
using namespace std;


// ============================================================================
// SphSnapshot::SphSnapshot
// ============================================================================
SphSnapshot::SphSnapshot(string auxfilename)
{
  allocated = false;
  nallocated = 0;
  Nsph = 0;
  t = 0.0;
  if (auxfilename != "") filename = auxfilename;
  LastUsed = time(NULL);
}



// ============================================================================
// SphSnapshot::~SphSnapshot
// ============================================================================
SphSnapshot::~SphSnapshot()
{
}



// ============================================================================
// SphSnapshot::AllocateBufferMemory
// ============================================================================
void SphSnapshot::AllocateBufferMemory(void)
{
  if (allocated) {
    if (Nsph > Nmax)
      DeallocateBufferMemory();
    else
      return;
  }

  if (ndim == 1) {
    x = new float[Nsph];
    vx = new float[Nsph];
    ax = new float[Nsph];
  }
  else if (ndim == 2) {
    x = new float[Nsph];
    y = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
  }
  else if (ndim == 3) {
    x = new float[Nsph];
    y = new float[Nsph];
    z = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    vz = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
    az = new float[Nsph];
  }
  
  m = new float[Nsph];
  h = new float[Nsph];
  rho = new float[Nsph];
  u = new float[Nsph];
  dudt = new float[Nsph];

  allocated = true;
  nallocated = 3*ndim + 5;
  Nmax = Nsph;

  return;
}



// ============================================================================
// SphSnapshot::DeallocateBufferMemory
// ============================================================================
void SphSnapshot::DeallocateBufferMemory(void)
{
  if (ndim == 1) {
    delete[] x;
    delete[] vx;
    delete[] ax;
  }
  else if (ndim == 2) {
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] ax;
    delete[] ay;
  }
  else if (ndim == 3) {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;
    delete[] ax;
    delete[] ay;
    delete[] az;
  }
  
  delete[] m;
  delete[] h;
  delete[] rho;
  delete[] u;
  delete[] dudt;

  allocated = false;
  nallocated = 0;

  return;
}



// ============================================================================
// SphSnapshot::CalculateMemoryUsage
// ============================================================================
int SphSnapshot::CalculateMemoryUsage(void)
{
  return Nsph*nallocated*sizeof(float);
}



// ============================================================================
// SphSnapshot::CopyDataFromSimulation
// ============================================================================
void SphSnapshot::CopyDataFromSimulation(int ndimaux, int Nsphaux, 
					 SphParticle *sphaux)
{
  ndim = ndimaux;
  Nsph = Nsphaux;

  AllocateBufferMemory();

  for (int i=0; i<Nsph; i++) {

    if (ndim == 1) {
      x[i] = (float) sphaux[i].r[0];
      vx[i] = (float) sphaux[i].v[0];
      ax[i] = (float) sphaux[i].a[0];
    }
    else if (ndim == 2) {
      x[i] = (float) sphaux[i].r[0];
      y[i] = (float) sphaux[i].r[1];
      vx[i] = (float) sphaux[i].v[0];
      vy[i] = (float) sphaux[i].v[1];
      ax[i] = (float) sphaux[i].a[0];
      ay[i] = (float) sphaux[i].a[1];
    }
    else if (ndim == 3) {
      x[i] = (float) sphaux[i].r[0];
      y[i] = (float) sphaux[i].r[1];
      z[i] = (float) sphaux[i].r[2];
      vx[i] = (float) sphaux[i].v[0];
      vy[i] = (float) sphaux[i].v[1];
      vz[i] = (float) sphaux[i].v[2];
      ax[i] = (float) sphaux[i].a[0];
      ay[i] = (float) sphaux[i].a[1];
      az[i] = (float) sphaux[i].a[2];
    }

    m[i] = (float) sphaux[i].m;
    h[i] = (float) sphaux[i].h;
    rho[i] = (float) sphaux[i].rho;
    u[i] = (float) sphaux[i].u;
    dudt[i] = (float) sphaux[i].dudt;

  }
  LastUsed = time(NULL);
  return;
}



// ============================================================================
// SphSnapshot::ExtractArray
// ============================================================================
void SphSnapshot::ExtractArray(string name, float** out_array, int* size_array)
{

  LastUsed = time(NULL);

  if (name == "x") *out_array = x;
  else if (name == "y") *out_array = y;
  else if (name == "z") *out_array = z;
  else if (name == "vx") *out_array = vx;
  else if (name == "vy") *out_array = vy;
  else if (name == "vz") *out_array = vz;
  else if (name == "ax") *out_array = ax;
  else if (name == "ay") *out_array = ay;
  else if (name == "az") *out_array = az;
  else if (name == "m") *out_array = m;
  else if (name == "h") *out_array = h;
  else if (name == "rho") *out_array = rho;
  else if (name == "u") *out_array = u;
  else if (name == "dudt") *out_array = dudt;
  else cout << "Warning: the selected array has not been recognized" << endl;

  *size_array = Nsph;

  return;
}



// ============================================================================
// SphSnapshot::ReadSnapshot
// ============================================================================
void SphSnapshot::ReadSnapshot(string format, SphSimulation *simulation) {

  simulation->ReadSnapshotFile(filename, format);
  CopyDataFromSimulation(simulation->simparams.intparams["ndim"],
			 simulation->sph->Nsph , simulation->sph->sphdata );
  t = simulation->t;

}