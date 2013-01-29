// ============================================================================
// SphSimulation.h
// ============================================================================


#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_


#include <map>
#include <string>
#include "Parameters.h"
#include "SimUnits.h"
#include "SphKernel.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "SphNeighbourSearch.h"
#include "SphIntegration.h"
using namespace std;


struct DomainBox {
  string x_boundary_lhs;
  string x_boundary_rhs;
  string y_boundary_lhs;
  string y_boundary_rhs;
  string z_boundary_lhs;
  string z_boundary_rhs;
  float boxmin[3];
  float boxmax[3];
  float boxsize[3];
  //float boxhalf[ndimmax];
  //float rmin[ndimmax];
  //float rmax[ndimmax];

};


// ============================================================================
// CLASS SphSimulation
// ============================================================================
class SphSimulation
{
 public:

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  SphSimulation();
  ~SphSimulation();

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  void GenerateIC(int);
  void Setup(void);
  void MainLoop(void);
  void Run(int,double);
  void ComputeBlockTimesteps(void);

  void SearchGhostParticles(void);
  void CreateGhostParticle(int,int,float,float);
  void CopyDataToGhosts(void);
  void CheckBoundaries(void);

  // Initial conditions routines
  // --------------------------------------------------------------------------
  void RandomBox(void);
  void ShockTube(void);

  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddRandomBox(int, float *, DomainBox);
  void AddRegularLattice(int, int *, float *, DomainBox);




#if !defined(FIXED_DIMENSIONS)
  int ndim;
  int vdim;
  int bdim;
#endif

  // Integer and physical Timestep counters
  // --------------------------------------------------------------------------
  int n;
  int Nsteps;
  double t;
  double timestep;

  // Name of parameters file and Parameters object that reads all data
  // --------------------------------------------------------------------------
  string paramfile;
  Parameters simparams;

  SimUnits simunits;
  Sph *sph;
  SphNeighbourSearch *sphneib;
  SphIntegration *sphint;
  SphSnapshot livesnap;
  DomainBox simbox;

};


#endif
