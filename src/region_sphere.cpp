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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "region_sphere.h"
#include "update.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

enum{CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

RegSphere::RegSphere(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-6,&arg[6]);

  xc = xscale*force->numeric(FLERR,arg[2]);
  yc = yscale*force->numeric(FLERR,arg[3]);
  zc = zscale*force->numeric(FLERR,arg[4]);

  rstr = NULL;
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    rstr = new char[n];
    strcpy(rstr,&arg[5][2]);
    radius = 0.0;
    rstyle = VARIABLE;
    varshape = 1;
    variable_check();
    shape_update();
  } else {
    radius = xscale*force->numeric(FLERR,arg[5]);
    rstyle = CONSTANT;
  }

  // error check

  if (radius < 0.0) error->all(FLERR,"Illegal region sphere command");

  // extent of sphere
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    extent_xlo = xc - radius;
    extent_xhi = xc + radius;
    extent_ylo = yc - radius;
    extent_yhi = yc + radius;
    extent_zlo = zc - radius;
    extent_zhi = zc + radius;
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
  tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegSphere::~RegSphere()
{
  delete [] rstr;
  delete [] contact;
}

/* ---------------------------------------------------------------------- */

void RegSphere::init()
{
  Region::init();
  if (rstr) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegSphere::inside(double x, double y, double z)
{
  double delx = x - xc;
  double dely = y - yc;
  double delz = z - zc;
  double r = sqrt(delx*delx + dely*dely + delz*delz);

  if (r <= radius) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from inner surface of sphere
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
   special case: no contact if x is at center of sphere
------------------------------------------------------------------------- */

int RegSphere::surface_interior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx*delx + dely*dely + delz*delz);
  if (r > radius || r == 0.0) return 0;

  double delta = radius - r;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx*(1.0-radius/r);
    contact[0].dely = dely*(1.0-radius/r);
    contact[0].delz = delz*(1.0-radius/r);
    contact[0].radius = -radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of sphere
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
------------------------------------------------------------------------- */

int RegSphere::surface_exterior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx*delx + dely*dely + delz*delz);
  if (r < radius) return 0;

  double delta = r - radius;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx*(1.0-radius/r);
    contact[0].dely = dely*(1.0-radius/r);
    contact[0].delz = delz*(1.0-radius/r);
    contact[0].radius = radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegSphere::shape_update()
{
  radius = xscale * input->variable->compute_equal(rvar);
  if (radius < 0.0)
    error->one(FLERR,"Variable evaluation in region gave bad value");
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegSphere::variable_check()
{
  rvar = input->variable->find(rstr);
  if (rvar < 0)
    error->all(FLERR,"Variable name for region sphere does not exist");
  if (!input->variable->equalstyle(rvar))
    error->all(FLERR,"Variable for region sphere is invalid style");
}


/* ----------------------------------------------------------------------
   Set values needed to calculate velocity due to shape changes.
   These values do not depend on the contact, so this function is
   called once per timestep by fix/wall/gran/region.

------------------------------------------------------------------------- */

void RegSphere::set_velocity_shape()
{
  xcenter[0] = xc;
  xcenter[1] = yc;
  xcenter[2] = zc;
  forward_transform(xcenter[0], xcenter[1], xcenter[2]);
  if (update->ntimestep > 0) rprev = prev[4];
  else rprev = radius;
  prev[4] = radius;
}



/* ----------------------------------------------------------------------
   add velocity due to shape change to wall velocity
------------------------------------------------------------------------- */

void RegSphere::velocity_contact_shape(double *vwall, double *xc)
{
  double delx, dely, delz; // Displacement of contact point in x,y,z

  delx = (xc[0] - xcenter[0])*(1 - rprev/radius);
  dely = (xc[1] - xcenter[1])*(1 - rprev/radius);
  delz = (xc[2] - xcenter[2])*(1 - rprev/radius);

  vwall[0] += delx/update->dt;
  vwall[1] += dely/update->dt;
  vwall[2] += delz/update->dt;
}

