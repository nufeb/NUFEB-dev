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

#include "fix_adhesion.h"

#include <math.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "group.h"
#include "memory.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdhesion::FixAdhesion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/adhesion command");

  virial_flag = 1;
  smin = 0.0;
  smax = 0.0;
  allocated = 0;

  ah = nullptr;

  smin = utils::numeric(FLERR,arg[3],true,lmp);
  smax = utils::numeric(FLERR,arg[4],true,lmp);
  if (smin > smax)
    error->all(FLERR, "Illegal value for smin/smax parameters in fix nufeb/adhesion");
}

/* ---------------------------------------------------------------------- */

FixAdhesion::~FixAdhesion()
{
  memory->destroy(ah);
}

/* ---------------------------------------------------------------------- */

int FixAdhesion::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ah") == 0) {
      if (!allocated) allocate();

      int ilo,ihi,jlo,jhi;
      utils::bounds(FLERR,arg[iarg+1],1,atom->ntypes,ilo,ihi,error);
      utils::bounds(FLERR,arg[iarg+2],1,atom->ntypes,jlo,jhi,error);

      double ah_one = utils::numeric(FLERR,arg[iarg+3],true,lmp);

      int count = 0;
      for (int i = ilo; i <= ihi; i++) {
        for (int j = MAX(jlo,i); j <= jhi; j++) {
          ah[i][j] = ah_one;
          ah[j][i] = ah_one;
          count++;
        }
      }
      if (count == 0) error->all(FLERR,"Incorrect args for ah coefficient");
      iarg += 4;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixAdhesion::init()
{
  if (!allocated) error->all(FLERR,"fix adhesion coeffs are not set");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void FixAdhesion::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(ah,atom->ntypes+1,atom->ntypes+1,"nufeb/cohesion:ah");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      ah[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixAdhesion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdhesion::post_force(int vflag)
{
  // energy and virial setup
  if (vflag) v_setup(vflag);
  else evflag = 0;

  compute<0>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int DISP>
void FixAdhesion::compute(int vflag)
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double *radius = atom->radius;
  double **f = atom->f;
  int newton_pair = force->newton_pair;
  
  NeighList *list = force->pair->list;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int i = 0; i < 6; i++)
    virial[i] = 0.0;
  
  for (int ii = 0; ii < nlocal; ii++) {
    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    int itype = type[i];
    double radi = radius[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;

      int jtype = type[j];
      double radj = radius[j];
      double radsum = radi + radj;
      double ccel = 0;

      if (rsq < (radsum + smax)*(radsum + smax)){
	double r = sqrt(rsq);
	double del = r - radsum;
	if (del > smin) {
	  double first = del*del+2*radi*del+2*radj*del;
	  double second = del*del+2*radi*del+2*radj*del+4*radi*radj;
	  ccel = -(ah[itype][jtype]*64*radi*radi*radi*radj*radj*radj*(del+radsum)) / (6.0*first*first*second*second);
	}
	else if (del >= 0 && del <= smin) {
	  double first = smin*smin+2*radi*smin+2*radj*smin;
	  double second = smin*smin+2*radi*smin+2*radj*smin+4*radi*radj;
	  ccel = -(ah[itype][jtype]*64*radi*radi*radi*radj*radj*radj*(del+radsum)) / (6.0*first*first*second*second);
	} else
	  ccel = 0;

	double rinv = 1/r;

	double ccelx = delx*ccel*rinv;
	double ccely = dely*ccel*rinv;
	double ccelz = delz*ccel*rinv;

	f[i][0] += ccelx;
	f[i][1] += ccely;
	f[i][2] += ccelz;

	if (newton_pair || j < nlocal){
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
