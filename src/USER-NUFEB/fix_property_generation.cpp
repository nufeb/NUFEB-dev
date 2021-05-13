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

#include "fix_property_generation.h"

#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyGeneration::FixPropertyGeneration(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix property/nufeb/generation command");

  create_attribute = 1;
  // use vprop if size_peratom_cols = 0
  size_peratom_cols = 0;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

void FixPropertyGeneration::init()
{
  for (int i = 0; i < atom->nlocal; i++) {
    vprop[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyGeneration::set_arrays(int j)
{
  vprop[j] = 0.0;
}

/* ----------------------------------------------------------------------
   update array values of two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyGeneration::update_arrays(int i, int j)
{
  vprop[j] = vprop[i] + 1;
}


