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

#include "fix_divide_bacillus.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_bacillus.h"
#include "error.h"
#include "force.h"
#include "lmptype.h"
#include "math_const.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "atom_masks.h"
#include "random_park.h"

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDivideBacillus::FixDivideBacillus(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Fix nufeb/monod/ecoli/wild requires "
      "atom style bacillus");

  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/divide/bacillus command");
  
  maxlength = force->numeric(FLERR, arg[3]);
  if (maxlength <= 0)
    error->all(FLERR, "Max division length cannot be less or equal to 0");
  seed = force->inumeric(FLERR, arg[4]);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

FixDivideBacillus::~FixDivideBacillus()
{
  delete random;
};

/* ---------------------------------------------------------------------- */

void FixDivideBacillus::compute()
{  
  int nlocal = atom->nlocal;
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      if (bonus->length >= maxlength) {
	double imass, ibiomass;
	double ilen, xp1[3], xp2[3];

//	double phiz = random->uniform() * 2e-8;

	double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
	double acircle = MY_PI*atom->radius[i]*atom->radius[i];
	double density = atom->rmass[i] / (vsphere + acircle * bonus->length);

	double oldx = atom->x[i][0];
	double oldy = atom->x[i][1];
	double oldz = atom->x[i][2];
	double old_len = bonus->length;

        imass = atom->rmass[i]/2;
        ibiomass = atom->biomass[i]/2;

        // conserve mass
        ilen = (imass / density - vsphere) / acircle;

	xp1[0] = xp1[1] = xp1[2] = 0.0;
	xp2[0] = xp2[1] = xp2[2] = 0.0;

	avec->get_pole_coords(i, xp1, xp2);

        // update daughter cell i
	double dl = (0.5*ilen + atom->radius[i]) / (0.5*old_len);
	atom->x[i][0] += (xp1[0] - oldx) * dl;
	atom->x[i][1] += (xp1[1] - oldy) * dl;
	atom->x[i][2] += (xp1[2] - oldz) * dl;

        atom->rmass[i] = imass;
        atom->biomass[i] = ibiomass;

        bonus->pole1[0] *= ilen / old_len;
        bonus->pole1[1] *= ilen / old_len;
        bonus->pole1[2] *= ilen / old_len;
        bonus->pole2[0] *= ilen / old_len;
        bonus->pole2[1] *= ilen / old_len;
        bonus->pole2[2] *= ilen / old_len;
        bonus->length = ilen;

        // create daughter cell j
        double *coord = new double[3];

	coord[0] = oldx + (xp2[0] - oldx) * dl;
	coord[1] = oldy + (xp2[1] - oldy) * dl;
	coord[2] = oldz + (xp2[2] - oldz) * dl;

        avec->create_atom(atom->type[i], coord);
        int j = atom->nlocal - 1;
        atom->bacillus[j] = 0;

        avec->set_bonus(j, bonus->pole1, bonus->diameter, bonus->quat, bonus->inertia);

        atom->tag[j] = 0;
        atom->mask[j] = atom->mask[i];
        atom->v[j][0] = atom->v[i][0];
        atom->v[j][1] = atom->v[i][1];
        atom->v[j][2] = atom->v[i][2];
	atom->f[j][0] = atom->f[i][0];
	atom->f[j][1] = atom->f[i][1];
	atom->f[j][2] = atom->f[i][2];
	atom->torque[j][0] = atom->torque[i][0];
	atom->torque[j][1] = atom->torque[i][1];
	atom->torque[j][2] = atom->torque[i][2];
	atom->angmom[j][0] = atom->angmom[i][0];
	atom->angmom[j][1] = atom->angmom[i][1];
	atom->angmom[j][2] = atom->angmom[i][2];
        atom->rmass[j] = imass;
        atom->biomass[j] = ibiomass;
        atom->radius[j] = atom->radius[i];

        modify->create_attribute(j);

        for (int m = 0; m < modify->nfix; m++)
          modify->fix[m]->update_arrays(i, j);

        delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
