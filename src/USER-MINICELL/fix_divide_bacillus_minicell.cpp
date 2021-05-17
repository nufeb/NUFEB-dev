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

#include "fix_divide_bacillus_minicell.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_bacillus.h"
#include "error.h"
#include "lmptype.h"
#include "math_const.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "domain.h"
#include "atom_masks.h"
#include "random_park.h"

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDivideBacillusMinicell::FixDivideBacillusMinicell(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Fix nufeb/divide/bacillus/minicell requires "
      "atom style bacillus");

  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/divide/bacillus/minicell command");
  
  imini = group->find(arg[3]);
  if (imini < 0)
    error->all(FLERR, "Can't find minicell group name");

  type = utils::inumeric(FLERR,arg[4],true,lmp);

  maxlength = utils::numeric(FLERR,arg[5],true,lmp);
  if (maxlength <= 0)
    error->all(FLERR, "Max division length must be greater than 0");
  prob = utils::numeric(FLERR,arg[6],true,lmp);

  if (prob < 0 || prob > 1)
    error->all(FLERR, "Minicell division probability must be between 0-1");
  seed = utils::numeric(FLERR,arg[7],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  maxradius = 0.0;
}

FixDivideBacillusMinicell::~FixDivideBacillusMinicell()
{
  delete random;
};

/* ---------------------------------------------------------------------- */

void FixDivideBacillusMinicell::compute()
{  
  int nlocal = atom->nlocal;
  int mini_mask = group->bitmask[imini];
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      if (bonus->length < maxlength) continue;

      double imass, jmass, ibiomass, jbiomass;
      double ilen, jlen, xp1[3], xp2[3];
      int n;

      double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
      double acircle = MY_PI*atom->radius[i]*atom->radius[i];
      double density = atom->rmass[i] / (vsphere + acircle * bonus->length);

      double oldx = atom->x[i][0];
      double oldy = atom->x[i][1];
      double oldz = atom->x[i][2];
      double old_len = bonus->length;

      xp1[0] = xp1[1] = xp1[2] = 0.0;
      xp2[0] = xp2[1] = xp2[2] = 0.0;

      avec->get_pole_coords(i, xp1, xp2);

      double prob_mini = random->uniform();
      // normal division
      if (prob_mini > prob) {
        imass = atom->rmass[i]/2;
        ibiomass = atom->biomass[i]/2;
        jmass = imass;
        jbiomass = ibiomass;
        // conserve mass
        ilen = (imass / density - vsphere) / acircle;;

        // update parent
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

        // create child
        double *coord = new double[3];

	coord[0] = oldx + (xp2[0] - oldx) * dl;
	coord[1] = oldy + (xp2[1] - oldy) * dl;
	coord[2] = oldz + (xp2[2] - oldz) * dl;

        avec->create_atom(atom->type[i], coord);
        n = atom->nlocal - 1;
        atom->bacillus[n] = 0;
        atom->mask[n] = atom->mask[i];
        avec->set_bonus(n, bonus->pole1, bonus->diameter, bonus->quat, bonus->inertia);

        delete[] coord;
      // abnormal division, generate one sphere (j) and one long rod (i)
      } else {
	double prob_pole = random->uniform();

	double *coord = new double[3];
	double idl, jdl, pole1[3];

	jmass = vsphere * density;
        jbiomass = atom->biomass[i] * jmass / atom->rmass[i];
        imass = atom->rmass[i] - jmass;
        ibiomass = atom->biomass[i] - jbiomass;
        ilen = (imass / density - vsphere) / acircle;
        idl = ilen / old_len;
	jdl = (ilen + 2*atom->radius[i]) / old_len;

	atom->rmass[i] = imass;
	atom->biomass[i] = ibiomass;
	bonus->length = ilen;
	bonus->pole1[0] *= ilen / old_len;
	bonus->pole1[1] *= ilen / old_len;
	bonus->pole1[2] *= ilen / old_len;
	bonus->pole2[0] *= ilen / old_len;
	bonus->pole2[1] *= ilen / old_len;
	bonus->pole2[2] *= ilen / old_len;

        if (prob_pole > 0.5) {
	  // update parent
	  atom->x[i][0] = xp1[0] + (oldx - xp1[0]) * idl;
	  atom->x[i][1] = xp1[1] + (oldy - xp1[1]) * idl;
	  atom->x[i][2] = xp1[2] + (oldz - xp1[2]) * idl;

	  // create child
	  coord[0] = oldx + (xp2[0] - oldx) * jdl;
	  coord[1] = oldy + (xp2[1] - oldy) * jdl;
	  coord[2] = oldz + (xp2[2] - oldz) * jdl;
        } else {
  	  // update parent
  	  atom->x[i][0] = xp2[0] + (oldx - xp2[0]) * idl;
  	  atom->x[i][1] = xp2[1] + (oldy - xp2[1]) * idl;
  	  atom->x[i][2] = xp2[2] + (oldz - xp2[2]) * idl;

  	  // create child
  	  coord[0] = oldx + (xp1[0] - oldx) * jdl;
  	  coord[1] = oldy + (xp1[1] - oldy) * jdl;
  	  coord[2] = oldz + (xp1[2] - oldz) * jdl;
        }

	avec->create_atom(type, coord);
	n = atom->nlocal - 1;
	atom->bacillus[n] = 0;
	atom->mask[n] = 1 | mini_mask;
	pole1[0] = pole1[1] = pole1[2] = 0;
	avec->set_bonus(n, pole1, bonus->diameter, bonus->quat, bonus->inertia);
	delete[] coord;
      }
      // set attributes for child
      atom->tag[n] = 0;
      atom->v[n][0] = atom->v[i][0];
      atom->v[n][1] = atom->v[i][1];
      atom->v[n][2] = atom->v[i][2];
      atom->f[n][0] = atom->f[i][0];
      atom->f[n][1] = atom->f[i][1];
      atom->f[n][2] = atom->f[i][2];
      atom->torque[n][0] = atom->torque[i][0];
      atom->torque[n][1] = atom->torque[i][1];
      atom->torque[n][2] = atom->torque[i][2];
      atom->angmom[n][0] = atom->angmom[i][0];
      atom->angmom[n][1] = atom->angmom[i][1];
      atom->angmom[n][2] = atom->angmom[i][2];
      atom->rmass[n] = jmass;
      atom->biomass[n] = jbiomass;
      atom->radius[n] = atom->radius[i];

      modify->create_attribute(n);
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

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixDivideBacillusMinicell::extract(const char *str, int &itype)
{
  if (strcmp(str,"radius") == 0) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->bacillus[i] >= 0) {
	double radius = atom->radius[i];
	if (radius > maxradius) maxradius = radius;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &maxradius, 1, MPI_DOUBLE, MPI_MAX, world);
    printf("maxradius=%e!\n",maxradius);
    return &maxradius;
  }
  return NULL;
}
