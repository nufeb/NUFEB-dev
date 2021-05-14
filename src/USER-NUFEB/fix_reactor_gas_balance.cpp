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

#include "fix_reactor_gas_balance.h"

#include <cstdio>
#include <cstring>
#include <cmath>

#include "error.h"
#include "grid.h"
#include "grid_masks.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "memory.h"

#include "fix_gas_liquid.h"
#include "modify.h"



using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReactorGasBalance::FixReactorGasBalance(LAMMPS *lmp, int narg, char **arg) :
  FixReactor(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal reactor/gas_balance command");

  igas = -1;

  q = 0.0;
  reactor_vhead = 1.0;

  nfix_gas_liquid = 0;
  fix_gas_liquid = NULL;

  igas = grid->find(arg[3]);
  if (igas < 0)
    error->all(FLERR, "Can't find substrate for reactor/gas_balance");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "reactor_vhead") == 0) {
      reactor_vhead = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_vhead <= 0)
	error->all(FLERR, "Reactor headspace volume (reactor_vhead) must be positive");
      iarg += 2;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixReactorGasBalance::~FixReactorGasBalance()
{
  delete [] fix_gas_liquid;
}

/* ---------------------------------------------------------------------- */

void FixReactorGasBalance::init()
{
  // allocate space for storing fixes
  fix_gas_liquid = new FixGasLiquid*[modify->nfix];

  // find fixes
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style, "nufeb/gas_liquid")) {
      fix_gas_liquid[nfix_gas_liquid++] = (FixGasLiquid *)modify->fix[i];
    }
  }
}

/* ----------------------------------------------------------------------
   return gas concentration in bulk
------------------------------------------------------------------------- */

double FixReactorGasBalance::compute_scalar()
{
  return grid->bulk[igas];
}

/* ----------------------------------------------------------------------
   update gas concentration in bulk
------------------------------------------------------------------------- */

void FixReactorGasBalance::compute()
{
  double **reac = grid->reac;
  double *bulk = grid->bulk;
  double sum_reac, ncells;

  for (int i = 0; i < nfix_gas_liquid; i++)
    q += fix_gas_liquid[i]->compute_scalar();

  sum_reac = 0;
  ncells = 0;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      sum_reac += reac[igas][i];
      ncells++;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum_reac, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &ncells, 1, MPI_INT, MPI_SUM, world);

  bulk[igas] += (sum_reac / ncells - (q / reactor_vhead * bulk[igas])) * update->dt;

}
