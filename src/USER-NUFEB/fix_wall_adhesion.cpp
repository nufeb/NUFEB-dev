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

#include <cmath>
#include <cstring>
#include "fix_wall_adhesion.h"
#include "atom.h"
#include "error.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallAdhesion::FixWallAdhesion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix nufeb/wall_adhesion command");

  ieps = group->find(arg[3]);
  if (ieps < 0)
    error->all(FLERR, "Can't find group name");
  eps_mask = group->bitmask[ieps];
  
  kn = utils::numeric(FLERR,arg[4],true,lmp);

  int iarg = 5;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],true,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],true,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],true,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],true,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = utils::numeric(FLERR,arg[iarg+1],true,lmp);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = utils::numeric(FLERR,arg[iarg+2],true,lmp);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = utils::numeric(FLERR,arg[iarg+1],true,lmp);
    iarg += 2;
  }
}

/* ---------------------------------------------------------------------- */

int FixWallAdhesion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallAdhesion::post_force(int vflag)
{
  double dx,dy,dz,del1,del2,delxy,delr,rsq;

  double wlo = lo;
  double whi = hi;

  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double eps_mass;
  double delta;
  double r, rinv, ccel, ccelx, ccely, ccelz;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & eps_mask) {
      eps_mass = rmass[i];
    }
    else {
      eps_mass = outer_mass[i];
    }

    if ((mask[i] & groupbit) && eps_mass > 0) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > outer_radius[i]) dz = cylradius;
        else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq <= 2*outer_radius[i]*outer_radius[i]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        delta = outer_radius[i] - r;
        ccel= delta*kn*eps_mass;
        ccelx = dx*ccel*rinv;
        ccely = dy*ccel*rinv;
        ccelz = dz*ccel*rinv;
        f[i][0] += ccelx;
        f[i][1] += ccely;
        f[i][2] += ccelz;
      }
    }
  }
}
