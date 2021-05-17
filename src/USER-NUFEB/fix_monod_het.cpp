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
#include "fix_monod_het.h"
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixMonodHET::FixMonodHET(LAMMPS *lmp, int narg, char **arg) :
  FixMonod(lmp, narg, arg)
{
  if (narg < 11)
    error->all(FLERR, "Illegal fix nufeb/monod/het command");
  
  dynamic_group_allow = 1;

  isub = -1;
  io2 = -1;
  ino2 = -1;
  ino3 = -1;

  sub_affinity = 0.0;
  o2_affinity = 0.0;
  no2_affinity = 0.0;
  no3_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;
  eps_yield = 0.0;
  anoxic = 1.0;
  eps_dens = 1.0;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");
  sub_affinity = utils::numeric(FLERR,arg[4],true,lmp);

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate name");
  o2_affinity = utils::numeric(FLERR,arg[6],true,lmp);

  ino2 = grid->find(arg[7]);
  if (ino2 < 0)
    error->all(FLERR, "Can't find substrate name");
  no2_affinity = utils::numeric(FLERR,arg[8],true,lmp);

  ino3 = grid->find(arg[9]);
  if (ino3 < 0)
    error->all(FLERR, "Can't find substrate name");
  no3_affinity = utils::numeric(FLERR,arg[10],true,lmp);
  
  int iarg = 11;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maintain") == 0) {
      maintain = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "epsyield") == 0) {
      eps_yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "anoxic") == 0) {
      anoxic = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "epsdens") == 0) {
      eps_dens = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/monod/het command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodHET::compute()
{ 
  if (reaction_flag && growth_flag) {
    update_cells<1, 1>();
    update_atoms();
  } else if (reaction_flag && !growth_flag) {
    update_cells<1, 0>();
  } else if (!reaction_flag && growth_flag) {
    update_cells<0, 1>();
    update_atoms();
  }
}

/* ---------------------------------------------------------------------- */

template <int Reaction, int Growth>
void FixMonodHET::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = growth * conc[isub][i] / (sub_affinity + conc[isub][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
    double tmp2 = anoxic * growth * conc[isub][i] / (sub_affinity + conc[isub][i]) * conc[ino3][i] / (no3_affinity + conc[ino3][i]) * o2_affinity / (o2_affinity + conc[io2][i]);
    double tmp3 = anoxic * growth * conc[isub][i] / (sub_affinity + conc[isub][i]) * conc[ino2][i] / (no2_affinity + conc[ino2][i]) * o2_affinity / (o2_affinity + conc[io2][i]);
    double tmp4 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);
    double tmp5 = 1 / 2.86 * maintain * anoxic * conc[ino3][i] / (no3_affinity + conc[ino3][i]) * o2_affinity / (o2_affinity + conc[io2][i]);
    double tmp6 = 1 / 1.17 * maintain * anoxic * conc[ino2][i] / (no2_affinity + conc[ino2][i]) * o2_affinity / (o2_affinity + conc[io2][i]);

    if (Reaction && !(grid->mask[i] & GHOST_MASK)) {
      reac[isub][i] -= 1 / yield * (tmp1 + tmp2 + tmp3) * dens[igroup][i];
      reac[io2][i] -= (1 - yield - eps_yield) / yield * tmp1 * dens[igroup][i] + tmp4 * dens[igroup][i];
      reac[ino2][i] -= (1 - yield - eps_yield) / (1.17 * yield) * tmp3 * dens[igroup][i] + tmp6 * dens[igroup][i];
      reac[ino3][i] -= (1 - yield - eps_yield) / (2.86 * yield) * tmp2 * dens[igroup][i] + tmp5 * dens[igroup][i];
    }
  
    if (Growth) {
      double ***grow = grid->growth;
      grow[igroup][i][0] = tmp1 + tmp2 + tmp3 - tmp4 - tmp5 - tmp6 - decay;
      grow[igroup][i][1] = (eps_yield / yield) * (tmp1 + tmp2 + tmp3);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodHET::update_atoms()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *biomass = atom->biomass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;
  double ***growth = grid->growth;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell = grid->cell(x[i]);
      const double density = rmass[i] /
	(four_thirds_pi * radius[i] * radius[i] * radius[i]);
      double ratio = rmass[i] / biomass[i];
      // forward Euler to update biomass and rmass
      biomass[i] = biomass[i] * (1 + growth[igroup][cell][0] * dt);
      rmass[i] = rmass[i] * (1 + growth[igroup][cell][0] * dt * ratio);
      outer_mass[i] = four_thirds_pi *
	(outer_radius[i] * outer_radius[i] * outer_radius[i] -
	 radius[i] * radius[i] * radius[i]) *
	eps_dens + growth[igroup][cell][1] * rmass[i] * dt;
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_radius[i] = pow(three_quarters_pi *
			    (rmass[i] / density + outer_mass[i] / eps_dens),
			    third);
    }
  }
}
