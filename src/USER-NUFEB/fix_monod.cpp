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

#include <cstdio>
#include <cstring>
#include <cmath>

#include "fix_monod.h"
#include "error.h"
#include "update.h"
#include "atom.h"
#include "grid.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixMonod::FixMonod(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  compute_flag = 1;
  reaction_flag = 1;
  growth_flag = 1;
  dt = 1.0;
}

/* ---------------------------------------------------------------------- */

int FixMonod::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compute") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	compute_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	compute_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "reaction") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	reaction_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	reaction_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "growth") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	growth_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	growth_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixMonod::init()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixMonod::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixMonod::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMonod::post_integrate()
{
  if (compute_flag)
    compute();
}


/* ---------------------------------------------------------------------- */

void FixMonod::update_atoms_coccus()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *biomass = atom->biomass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell = grid->cell(x[i]);
      const double density = rmass[i] /
    (four_thirds_pi * radius[i] * radius[i] * radius[i]);
      double growth = grid->growth[igroup][cell][0];
      // forward Euler to update rmass
      rmass[i] = rmass[i] * (1 + growth * dt);
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_mass[i] = 0;
      outer_radius[i] = radius[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonod::update_atoms_bacillus(AtomVecBacillus *&avec)
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *biomass = atom->biomass;

  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
      double acircle = MY_PI*atom->radius[i]*atom->radius[i];

      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];
      double length = bonus->length;

      double new_length;
      const int cell = grid->cell(x[i]);
      const double density = rmass[i] /	(vsphere + acircle * bonus->length);
      double growth = grid->growth[igroup][cell][0];
      // forward Eular to update rmass
      rmass[i] = rmass[i] * (1 + growth * dt);
      new_length = (rmass[i] - density * vsphere) / (density * acircle);
      bonus->length = new_length;
      // update coordinates of two poles
      double *pole1 = bonus->pole1;
      double *pole2 = bonus->pole2;

      pole1[0] *= new_length/length;
      pole1[1] *= new_length/length;
      pole1[2] *= new_length/length;
      pole2[0] *= new_length/length;
      pole2[1] *= new_length/length;
      pole2[2] *= new_length/length;
    }
  }
}
