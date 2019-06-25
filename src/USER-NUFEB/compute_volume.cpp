/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_volume.h"
#include "atom.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

#define FOURTHIRDSPI 4.18879020478

/* ---------------------------------------------------------------------- */

ComputeVolume::ComputeVolume(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute nufeb/volume command");

  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeVolume::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  scalar = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    double r = atom->radius[i];
    scalar += FOURTHIRDSPI * r * r * r;
  }
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world); 
  return scalar;
}
