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

#include "fix_gas_liquid.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "error.h"

#include "grid.h"
#include "grid_masks.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGasLiquid::FixGasLiquid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR,"Illegal fix nufeb/gas_liquid command");

  if (!grid->reactor_flag)
    error->all(FLERR,"Fix reactor requires nufeb/reactor grid style");

  compute_flag = 1;
  scalar_flag = 1;

  iliquid = -1;
  igas = -1;

  kga = 0.0;
  h = 1.0;
  temp = 1.0;
  reactor_vhead = 1.0;
  reactor_pres = 1.0;
  mw = 1.0;
  rg = 1.0;

  iliquid = grid->find(arg[3]);
  if (iliquid < 0)
    error->all(FLERR, "Can't find substrate(liquid) name");

  igas = grid->find(arg[4]);
  if (igas < 0)
    error->all(FLERR, "Can't find substrate(gas) name");

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kga") == 0) {
      kga = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "h") == 0) {
      h = force->numeric(FLERR, arg[iarg+1]);
      if (h <= 0)
	error->all(FLERR, "Henry's law solubility constant (H) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "temp") == 0) {
      temp = force->numeric(FLERR, arg[iarg+1]);
      if (temp <= 0)
	error->all(FLERR, "Temperature (temp) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_vhead") == 0) {
      reactor_vhead = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_vhead <= 0)
	error->all(FLERR, "Reactor headspace volume (reactor_vhead) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_pres") == 0) {
      reactor_pres = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_pres <= 0)
	error->all(FLERR, "Reactor headspace pressure (reactor_pres) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rg") == 0) {
      rg = force->numeric(FLERR, arg[iarg+1]);
      if (rg <= 0)
	error->all(FLERR, "Ideal gas constant must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "mw") == 0) {
      mw = force->numeric(FLERR, arg[iarg+1]);
      if (mw <= 0)
	error->all(FLERR, "Molar mass must be positive");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/gas_liquid command");
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixGasLiquid::modify_param(int narg, char **arg)
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
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

int FixGasLiquid::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGasLiquid::post_integrate()
{
  if (compute_flag)
    compute();
}

/* ---------------------------------------------------------------------- */

double FixGasLiquid::compute_scalar()
{
  double result = 0.0;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      result -= grid->reac[igas][i] * reactor_vhead;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  result /= reactor_pres;

  return result;
}


/* ---------------------------------------------------------------------- */

void FixGasLiquid::compute()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double p_g2l, n_l2g;
  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      p_g2l = kga * (grid->bulk[igas] - conc[iliquid][i]/(h * mw));
      n_l2g = -p_g2l * reactor_vhead / (rg * temp);
      // update reaction rates
      reac[igas][i] += p_g2l;
      reac[iliquid][i] += n_l2g * mw / vol;
    }
  }
}
