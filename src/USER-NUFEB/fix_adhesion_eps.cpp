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

#include "fix_adhesion_eps.h"

#include <math.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "group.h"
#include "pair.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DEFAULT,SQUARE};

/* ---------------------------------------------------------------------- */

FixEPSAdhesion::FixEPSAdhesion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/adhesion/eps command");

  virial_flag = 1;
  
  ieps = -1;
  ke = 0.0;
  disp = DEFAULT;

  ieps = group->find(arg[3]);
  if (ieps < 0)
    error->all(FLERR, "Can't find group in fix nufeb/adhesion/eps");
  ke = utils::numeric(FLERR,arg[4],true,lmp);

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "displace") == 0) {
      if (strcmp(arg[iarg+1], "default") == 0) {
	disp = DEFAULT;
      } else if (strcmp(arg[iarg+1], "square") == 0) {
	disp = SQUARE;
      } else {
	error->all(FLERR, "Illegal value for displace parameter in fix nufeb/adhesion/eps");
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal parameter in fix nufeb/adhesion/eps");
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixEPSAdhesion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdhesion::post_force(int vflag)
{
  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (disp == DEFAULT) {
    compute<0>(vflag);
  } else if (disp == SQUARE) {
    compute<1>(vflag);
  }
}

/* ---------------------------------------------------------------------- */

template <int DISP>
void FixEPSAdhesion::compute(int vflag)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double *outer_mass = atom->outer_mass;
  double *outer_radius = atom->outer_radius;
  double **f = atom->f;
  int newton_pair = force->newton_pair;
  
  NeighList *list = force->pair->list;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int epsmask = group->bitmask[ieps];

  for (int i = 0; i < 6; i++)
    virial[i] = 0.0;
  
  for (int ii = 0; ii < nlocal; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];

      double radi = radius[i];
      double oradi = outer_radius[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];

      double epsi = 0;
      if (mask[i] & epsmask)
	epsi = rmass[i];
      else
	epsi = outer_mass[i];

      for (int jj = 0; jj < jnum; jj++) {
	int j = jlist[jj];
	double delx = xtmp - x[j][0];
	double dely = ytmp - x[j][1];
	double delz = ztmp - x[j][2];
	double rsq = delx * delx + dely * dely + delz * delz;

	double radj = radius[j];
	double oradj = outer_radius[j];

	double epsj = 0;
	if (mask[j] & epsmask)
	  epsj = rmass[j];
	else
	  epsj = outer_mass[j];

	double radsum = radi + radj;
	double oradsum = oradi + oradj;
	double masssum = epsi + epsj;
	double r = sqrt(rsq);
	// double del = r - radsum;
	double del = r - 0.5 * (radsum + oradsum);
	double rinv = 1 / r;

	double ccel = 0;
	if (r > radsum && r < oradsum) {
	  if (DISP == DEFAULT)
	    ccel = -masssum * ke * del;
	  if (DISP == SQUARE)
	    ccel = -masssum * ke * (radsum / r) * (radsum / r);
	}

	double ccelx = delx * ccel * rinv;
	double ccely = dely * ccel * rinv;
	double ccelz = delz * ccel * rinv;

	f[i][0] += ccelx;
	f[i][1] += ccely;
	f[i][2] += ccelz;

	if (newton_pair || j < nlocal) {
	  f[j][0] -= ccelx;
	  f[j][1] -= ccely;
	  f[j][2] -= ccelz;
	}

	if (evflag) {
	  double v[6];
	  v[0] = delx*ccelx;
	  v[1] = dely*ccely;
	  v[2] = delz*ccelz;
	  v[3] = delx*ccely;
	  v[4] = delx*ccelz;
	  v[5] = dely*ccelz;
	  v_tally(i, v);
	}
      }
    }
  }
}
