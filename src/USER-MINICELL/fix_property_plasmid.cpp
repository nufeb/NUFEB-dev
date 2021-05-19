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

#include "fix_property_plasmid.h"

#include <cstdlib>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "math_const.h"
#include "random_park.h"
#include "atom_vec_bacillus.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyPlasmid::FixPropertyPlasmid(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix nufeb/property/plasmid command");

  compute_flag = 1;
  scalar_flag = 1;

  replication = force->numeric(FLERR, arg[3]);
  seed = force->inumeric(FLERR, arg[4]);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  avec = NULL;
  // col 0 plasmid copy number, col 2 reaction time
  size_peratom_cols = 2;
  grow_arrays(atom->nmax);

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
}

/* ---------------------------------------------------------------------- */

void FixPropertyPlasmid::init()
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (random->uniform() > 0.2)
      aprop[i][0] = 2.0;
    else
      aprop[i][0] = 2.0;

    aprop[i][1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */

void FixPropertyPlasmid::compute()
{
  double propensity;
  double dt, t;
  double next = update->ntimestep * (1 + update->dt);
  double current = update->ntimestep * update->dt;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      while (aprop[i][1] < next) {
	if (aprop[i][1] != 0)
	  aprop[i][0]++;
	propensity = aprop[i][0] * replication;
	dt = -log(random->uniform())/propensity;
	aprop[i][1] += dt;
      }
    }
  }
}


/* ----------------------------------------------------------------------
   update plasmid copy numbers in two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyPlasmid::update_arrays(int i, int j)
{
  double vi, vj;
  int iplasmids, jplasmids;
  double vspherei, vspherej;
  double ratio;
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  vspherei = four_thirds_pi * atom->radius[i] * atom->radius[i] * atom->radius[i];
  vspherej = four_thirds_pi * atom->radius[j] * atom->radius[j] * atom->radius[j];

  if (atom->sphere_flag) {
    vi = vspherei;
    vj = vspherej;
  } else if (atom->bacillus_flag) {
    int ibac = atom->bacillus[i];
    int jbac = atom->bacillus[j];

    AtomVecBacillus::Bonus *ibouns = &avec->bonus[ibac];
    AtomVecBacillus::Bonus *jbouns = &avec->bonus[jbac];

    double acirclei = MY_PI * atom->radius[i] * atom->radius[i];
    double acirclej = MY_PI * atom->radius[j] * atom->radius[j];

    vi = vspherei + acirclei * ibouns->length;
    vj = vspherej + acirclej * jbouns->length;
  } else
    error->all(FLERR, "Fix nufeb/property/plasmid only supports coccus, sphere and bacillus atom styles");

  ratio = vi / (vi + vj);
  iplasmids = round(aprop[i][0] * ratio);
  aprop[i][0] = iplasmids;
  aprop[j][0] = 1 - iplasmids;
}

/* ----------------------------------------------------------------------
   compute average plasmid copy number
------------------------------------------------------------------------- */
double FixPropertyPlasmid::compute_scalar() {
  double result = 0.0;
  int n = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit)) {
      result += aprop[i][0];
      n++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_INT, MPI_SUM, world);

  if (n > 0) result /= n;
  else result = 0;

  return result;
}


