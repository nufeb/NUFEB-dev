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

#include "fix_reactor_solute_balance.h"

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

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReactorSoluteBalance::FixReactorSoluteBalance(LAMMPS *lmp, int narg, char **arg) :
  FixReactor(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix reactor/solute_balance command");

  iliq = -1;

  q = 0.0;
  rvol = 1.0;
  reactor_af = 0.0;
  domain_af = 1.0;
  inlet = 0.0;

  iliq = grid->find(arg[3]);
  if (iliq < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "q") == 0) {
      q = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_vol") == 0) {
      rvol = force->numeric(FLERR, arg[iarg+1]);
      if (rvol <= 0)
        error->all(FLERR, "Bioreactor volume must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_af") == 0) {
	reactor_af = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_af < 0)
        error->all(FLERR, "Reactor biofilm surface area cannot be negative");
      iarg += 2;
    } else if (strcmp(arg[iarg], "domain_af") == 0) {
      if (strcmp(arg[iarg+1],"xy") == 0) {
	domain_af = (domain->boxhi[0] - domain->boxlo[0]) * (domain->boxhi[1] - domain->boxlo[1]);
      } else if (strcmp(arg[iarg+1],"yz") == 0) {
	domain_af = (domain->boxhi[1] - domain->boxlo[1]) * (domain->boxhi[2] - domain->boxlo[2]);
      } else if (strcmp(arg[iarg+1],"xz") == 0) {
	domain_af = (domain->boxhi[0] - domain->boxlo[0]) * (domain->boxhi[2] - domain->boxlo[2]);
      } else
	 error->all(FLERR,"Illegal fix reactor/solute_balance command");
      iarg += 2;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixReactorSoluteBalance::init()
{
  inlet = grid->bulk[iliq];
}

/* ----------------------------------------------------------------------
   return solute concentration in bulk
------------------------------------------------------------------------- */

double FixReactorSoluteBalance::compute_scalar()
{
  return grid->bulk[iliq];
}

/* ----------------------------------------------------------------------
   update solute concentration in bulk
------------------------------------------------------------------------- */

void FixReactorSoluteBalance::compute()
{
  double **reac = grid->reac;
  double *bulk = grid->bulk;
  double sum_reac;
  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  sum_reac = 0;
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      sum_reac += reac[iliq][i];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum_reac, 1, MPI_DOUBLE, MPI_SUM, world);

  bulk[iliq] += ((q / rvol) * (inlet - bulk[iliq]) + ((reactor_af * sum_reac * vol) / (rvol * domain_af))) * update->dt;
}
