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
#include "fix_monod_cyano.h"
#include "atom_vec_bacillus.h"
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixMonodCyano::FixMonodCyano(LAMMPS *lmp, int narg, char **arg) :
  FixMonod(lmp, narg, arg)
{
  avec = nullptr;
  avec = (AtomVecBacillus *) atom->style_match("bacillus");

  if (narg < 10)
    error->all(FLERR, "Illegal fix nufeb/monod/cyano command");

  dynamic_group_allow = 1;

  ilight = -1;
  ico2 = -1;
  igco2 = -1;
  isuc = -1;
  io2 = -1;

  light_affinity = 0.0;
  co2_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;
  suc_exp = 0.0;

  gco2_flag = 0;

  ilight = grid->find(arg[3]);
  if (ilight < 0)
    error->all(FLERR, "Can't find substrate(light) name");
  light_affinity = utils::numeric(FLERR,arg[4],true,lmp);
  if (light_affinity <= 0)
    error->all(FLERR, "light affinity must be greater than zero");

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate(o2) name");

  ico2 = grid->find(arg[6]);
  if (ico2 < 0)
    error->all(FLERR, "Can't find substrate(co2) name");
  co2_affinity = utils::numeric(FLERR,arg[7],true,lmp);
  if (co2_affinity <= 0)
    error->all(FLERR, "co2 affinity must be greater than zero");

  isuc = grid->find(arg[8]);
  if (isuc < 0)
    error->all(FLERR, "Can't find substrate(sucrose) name");

  igco2 = grid->find(arg[9]);
  if (igco2 < 0)
    error->all(FLERR, "Can't find substrate(gco2) name");

  int iarg = 10;
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
    } else if (strcmp(arg[iarg], "suc_exp") == 0) {
      suc_exp = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gco2_flag") == 0) {
      gco2_flag = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/monod/cyano command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodCyano::compute()
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
void FixMonodCyano::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    // cyanobacterial growth rate based on light(sub) and co2
    double tmp1 = growth * conc[ilight][i] / (light_affinity + conc[ilight][i]) * conc[ico2][i] / (co2_affinity + conc[ico2][i]);
    // sucrose export-induced growth reduction
    double tmp2 = 0.2 * tmp1 * suc_exp;
    double tmp3 = 4 * tmp1 * suc_exp;

    if (Reaction && !(grid->mask[i] & GHOST_MASK)) {
      // nutrient utilization
      reac[ilight][i] -= 1 / yield * (tmp1 + tmp3) * dens[igroup][i];
      reac[ico2][i] -= 1 / yield * (tmp1 + tmp3) * dens[igroup][i];
      reac[io2][i] -= 0.1 * maintain * dens[igroup][i];
      // oxygen evolution
      reac[io2][i] +=  (0.727 / yield) * (tmp1 + tmp3) * dens[igroup][i];
      // sucrose export
      reac[isuc][i] += 0.65 / yield * tmp3 * dens[igroup][i];

      // co2 dissolution
      if (gco2_flag == 1)
        reac[ico2][i] += (4.4e-6 * conc[igco2][i]) - (4.4e-6 * conc[ico2][i]);
    }

    if (Growth) {
      grid->growth[igroup][i][0] = tmp1 - tmp2 - decay - maintain;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodCyano::update_atoms()
{
  if (atom->coccus_flag) {
    update_atoms_coccus();
  } else {
    update_atoms_bacillus(avec);
  }
}
