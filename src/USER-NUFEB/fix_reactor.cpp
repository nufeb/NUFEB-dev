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

#include "fix_reactor.h"

#include <string.h>
#include "error.h"
#include "update.h"
#include "grid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReactor::FixReactor(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!grid->reactor_flag)
    error->all(FLERR,"Fix reactor requires nufeb/reactor grid style");

  compute_flag = 1;
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

int FixReactor::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixReactor::modify_param(int narg, char **arg)
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

void FixReactor::post_integrate()
{
  if (compute_flag)
    compute();
}
