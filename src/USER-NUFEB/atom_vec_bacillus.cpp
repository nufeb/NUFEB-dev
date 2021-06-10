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

#include "atom_vec_bacillus.h"

#include <cstdlib>
#include <math.h>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "memory.h"
#include "math_extra.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1.0e-7

enum{SPHERE,ROD};       // also in DumpImage
/* ---------------------------------------------------------------------- */

AtomVecBacillus::AtomVecBacillus(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = comm_f_only = 0;
  size_forward = 7;
  size_reverse = 6;
  size_border = 25;
  size_velocity = 6;
  size_data_atom = 8;
  size_data_vel = 7;
  size_data_bonus = 11;
  xcol_data = 5;

  atom->bacillus_flag = 1;
  atom->rmass_flag = 1;
  atom->biomass_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;
  atom->radius_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = NULL;
}

/* ---------------------------------------------------------------------- */

AtomVecBacillus::~AtomVecBacillus()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecBacillus::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  biomass = memory->grow(atom->biomass,nmax,"atom:biomass");
  angmom = memory->grow(atom->angmom,nmax,3,"atom:angmom");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
  bacillus = memory->grow(atom->bacillus,nmax,"atom:bacillus");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecBacillus::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  radius = atom->radius; rmass = atom->rmass; biomass = atom->biomass;
  angmom = atom->angmom; torque = atom->torque;
  bacillus = atom->bacillus;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecBacillus::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
   if delflag and atom J has bonus data, then delete it
------------------------------------------------------------------------- */

void AtomVecBacillus::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  radius[j] = radius[i];
  rmass[j] = rmass[i];
  biomass[j] = biomass[i];
  angmom[j][0] = angmom[i][0];
  angmom[j][1] = angmom[i][1];
  angmom[j][2] = angmom[i][2];

  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && bacillus[j] >= 0) {
    copy_bonus(nlocal_bonus-1,bacillus[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (bacillus[i] >= 0 && i != j) bonus[bacillus[i]].ilocal = j;
  bacillus[j] = bacillus[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset body that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecBacillus::copy_bonus(int i, int j)
{
  bacillus[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecBacillus::clear_bonus()
{
  nghost_bonus = 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_comm(int n, int *list, double *buf,
                          int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      if (bacillus[j] >= 0) {
        quat = bonus[bacillus[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      if (bacillus[j] >= 0) {
        quat = bonus[bacillus[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_comm_vel(int n, int *list, double *buf,
                              int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *quat;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      if (bacillus[j] >= 0) {
        quat = bonus[bacillus[j]].quat;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
      }
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = angmom[j][0];
      buf[m++] = angmom[j][1];
      buf[m++] = angmom[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        if (bacillus[j] >= 0) {
          quat = bonus[bacillus[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        if (bacillus[j] >= 0) {
          quat = bonus[bacillus[j]].quat;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
        }
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (bacillus[j] >= 0) {
      quat = bonus[bacillus[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (bacillus[i] >= 0) {
      quat = bonus[bacillus[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    if (bacillus[i] >= 0) {
      quat = bonus[bacillus[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    angmom[i][0] = buf[m++];
    angmom[i][1] = buf[m++];
    angmom[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    if (bacillus[i] >= 0) {
      quat = bonus[bacillus[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_border(int n, int *list, double *buf,
                            int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = biomass[j];
      if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[bacillus[j]].quat;
        inertia = bonus[bacillus[j]].inertia;
        pole1 = bonus[bacillus[j]].pole1;
        pole2 = bonus[bacillus[j]].pole2;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = pole1[0];
        buf[m++] = pole1[1];
        buf[m++] = pole1[2];
        buf[m++] = pole2[0];
        buf[m++] = pole2[1];
        buf[m++] = pole2[2];
        buf[m++] = bonus[bacillus[j]].length;
        buf[m++] = bonus[bacillus[j]].diameter;
      }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = biomass[j];
      if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[bacillus[j]].quat;
        inertia = bonus[bacillus[j]].inertia;
        pole1 = bonus[bacillus[j]].pole1;
        pole2 = bonus[bacillus[j]].pole2;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = pole1[0];
        buf[m++] = pole1[1];
        buf[m++] = pole1[2];
        buf[m++] = pole2[0];
        buf[m++] = pole2[1];
        buf[m++] = pole2[2];
        buf[m++] = bonus[bacillus[j]].length;
        buf[m++] = bonus[bacillus[j]].diameter;
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_border_vel(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = biomass[j];
      if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        quat = bonus[bacillus[j]].quat;
        inertia = bonus[bacillus[j]].inertia;
        pole1 = bonus[bacillus[j]].pole1;
        pole2 = bonus[bacillus[j]].pole2;
        buf[m++] = quat[0];
        buf[m++] = quat[1];
        buf[m++] = quat[2];
        buf[m++] = quat[3];
        buf[m++] = inertia[0];
        buf[m++] = inertia[1];
        buf[m++] = inertia[2];
        buf[m++] = pole1[0];
        buf[m++] = pole1[1];
        buf[m++] = pole1[2];
        buf[m++] = pole2[0];
        buf[m++] = pole2[1];
        buf[m++] = pole2[2];
        buf[m++] = bonus[bacillus[j]].length;
        buf[m++] = bonus[bacillus[j]].diameter;
      }
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = angmom[j][0];
      buf[m++] = angmom[j][1];
      buf[m++] = angmom[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = biomass[j];
        if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[bacillus[j]].quat;
          inertia = bonus[bacillus[j]].inertia;
          pole1 = bonus[bacillus[j]].pole1;
          pole2 = bonus[bacillus[j]].pole2;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
          buf[m++] = pole1[0];
          buf[m++] = pole1[1];
          buf[m++] = pole1[2];
          buf[m++] = pole2[0];
          buf[m++] = pole2[1];
          buf[m++] = pole2[2];
          buf[m++] = bonus[bacillus[j]].length;
          buf[m++] = bonus[bacillus[j]].diameter;
        }
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = biomass[j];
        if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
        else {
          buf[m++] = ubuf(1).d;
          quat = bonus[bacillus[j]].quat;
          inertia = bonus[bacillus[j]].inertia;
          pole1 = bonus[bacillus[j]].pole1;
          pole2 = bonus[bacillus[j]].pole2;
          buf[m++] = quat[0];
          buf[m++] = quat[1];
          buf[m++] = quat[2];
          buf[m++] = quat[3];
          buf[m++] = inertia[0];
          buf[m++] = inertia[1];
          buf[m++] = inertia[2];
          buf[m++] = pole1[0];
          buf[m++] = pole1[1];
          buf[m++] = pole1[2];
          buf[m++] = pole2[0];
          buf[m++] = pole2[1];
          buf[m++] = pole2[2];
          buf[m++] = bonus[bacillus[j]].length;
          buf[m++] = bonus[bacillus[j]].diameter;
        }
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = angmom[j][0];
        buf[m++] = angmom[j][1];
        buf[m++] = angmom[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
    buf[m++] = biomass[j];
    if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[bacillus[j]].quat;
      inertia = bonus[bacillus[j]].inertia;
      pole1 = bonus[bacillus[j]].pole1;
      pole2 = bonus[bacillus[j]].pole2;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
      buf[m++] = pole1[0];
      buf[m++] = pole1[1];
      buf[m++] = pole1[2];
      buf[m++] = pole2[0];
      buf[m++] = pole2[1];
      buf[m++] = pole2[2];
      buf[m++] = bonus[bacillus[j]].length;
      buf[m++] = bonus[bacillus[j]].diameter;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_border(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    biomass[i] = buf[m++];
    bacillus[i] = (int) ubuf(buf[m++]).i;
    if (bacillus[i] == 0) bacillus[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      pole1 = bonus[j].pole1;
      pole2 = bonus[j].pole2;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      pole1[0] = buf[m++];
      pole1[1] = buf[m++];
      pole1[2] = buf[m++];
      pole2[0] = buf[m++];
      pole2[1] = buf[m++];
      pole2[2] = buf[m++];
      bonus[j].length = buf[m++];
      bonus[j].diameter = buf[m++];
      bonus[j].ilocal = i;
      bacillus[i] = j;
      nghost_bonus++;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_border_vel(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    biomass[i] = buf[m++];
    bacillus[i] = (int) ubuf(buf[m++]).i;
    if (bacillus[i] == 0) bacillus[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      pole1 = bonus[j].pole1;
      pole2 = bonus[j].pole2;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      pole1[0] = buf[m++];
      pole1[1] = buf[m++];
      pole1[2] = buf[m++];
      pole2[0] = buf[m++];
      pole2[1] = buf[m++];
      pole2[2] = buf[m++];
      bonus[j].length = buf[m++];
      bonus[j].diameter = buf[m++];
      bonus[j].ilocal = i;
      bacillus[i] = j;
      nghost_bonus++;
    }
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    angmom[i][0] = buf[m++];
    angmom[i][1] = buf[m++];
    angmom[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    biomass[i] = buf[m++];
    bacillus[i] = (int) ubuf(buf[m++]).i;
    if (bacillus[i] == 0) bacillus[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();
      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      pole1 = bonus[j].pole1;
      pole2 = bonus[j].pole2;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      pole1[0] = buf[m++];
      pole1[1] = buf[m++];
      pole1[2] = buf[m++];
      pole2[0] = buf[m++];
      pole2[1] = buf[m++];
      pole2[2] = buf[m++];
      bonus[j].length = buf[m++];
      bonus[j].diameter = buf[m++];
      bonus[j].ilocal = i;
      bacillus[i] = j;
      nghost_bonus++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = biomass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (bacillus[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = bacillus[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    double *pole1= bonus[j].pole1;
    double *pole2= bonus[j].pole2;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = pole1[0];
    buf[m++] = pole1[1];
    buf[m++] = pole1[2];
    buf[m++] = pole2[0];
    buf[m++] = pole2[1];
    buf[m++] = pole2[2];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].diameter;
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  biomass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  bacillus[nlocal] = (int) ubuf(buf[m++]).i;
  if (bacillus[nlocal] == 0) bacillus[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    double *pole1= bonus[nlocal_bonus].pole1;
    double *pole2= bonus[nlocal_bonus].pole2;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    pole1[0] = buf[m++];
    pole1[1] = buf[m++];
    pole1[2] = buf[m++];
    pole2[0] = buf[m++];
    pole2[1] = buf[m++];
    pole2[2] = buf[m++];
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].diameter = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    bacillus[nlocal] = nlocal_bonus++;
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecBacillus::size_restart()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++)
    if (bacillus[i] >= 0) n += 33;
    else n += 18;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = biomass[i];
  buf[m++] = angmom[i][0];
  buf[m++] = angmom[i][1];
  buf[m++] = angmom[i][2];

  if (bacillus[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = bacillus[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    double *pole1= bonus[j].pole1;
    double *pole2= bonus[j].pole2;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = pole1[0];
    buf[m++] = pole1[1];
    buf[m++] = pole1[2];
    buf[m++] = pole2[0];
    buf[m++] = pole2[1];
    buf[m++] = pole2[2];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].diameter;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecBacillus::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  biomass[nlocal] = buf[m++];
  angmom[nlocal][0] = buf[m++];
  angmom[nlocal][1] = buf[m++];
  angmom[nlocal][2] = buf[m++];

  bacillus[nlocal] = (int) ubuf(buf[m++]).i;
  if (bacillus[nlocal] == 0) bacillus[nlocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    double *pole1= bonus[nlocal_bonus].pole1;
    double *pole2= bonus[nlocal_bonus].pole2;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    pole1[0] = buf[m++];
    pole1[1] = buf[m++];
    pole1[2] = buf[m++];
    pole2[0] = buf[m++];
    pole2[1] = buf[m++];
    pole2[2] = buf[m++];
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].diameter = buf[m++];
    bonus[nlocal_bonus].ilocal = nlocal;
    bacillus[nlocal] = nlocal_bonus++;
  }

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBacillus::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  radius[nlocal] = 0.5e-6;
  rmass[nlocal] = 1.0;
  biomass[nlocal] = 1.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;
  bacillus[nlocal] = -1;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBacillus::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  bacillus[nlocal] = atoi(values[2]);
  if (bacillus[nlocal] == 0) bacillus[nlocal] = -1;
  else if (bacillus[nlocal] == 1) bacillus[nlocal] = 0;
  else error->one(FLERR,"Invalid bacillusflag in Atoms section of data file");
  // use density to set mass here, correct later if bonus is defined
  rmass[nlocal] = atof(values[3]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Density must be greater than 0");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  angmom[nlocal][0] = 0.0;
  angmom[nlocal][1] = 0.0;
  angmom[nlocal][2] = 0.0;
  radius[nlocal] = 0.0;

  double ratio = atof(values[7]);
  if (ratio < 0 || ratio > 1)
    error->one(FLERR,"Biomass/Mass (dry/wet weight) ratio must be between 0-1");
  biomass[nlocal] = ratio;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecBacillus::data_atom_hybrid(int nlocal, char **values)
{
  bacillus[nlocal] = atoi(values[0]);
  if (bacillus[nlocal] == 0) bacillus[nlocal] = -1;
  else if (bacillus[nlocal] == 1) bacillus[nlocal] = 0;
  else error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rmass[nlocal] = atof(values[1]);
  if (rmass[nlocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  double ratio = atof(values[2]);
  if (ratio < 0 || ratio > 1)
    error->one(FLERR,"Biomass/Mass (dry/wet weight) ratio must be between 0-1");
  biomass[nlocal] = ratio;

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one body from Bacilli section of data file
------------------------------------------------------------------------- */

void AtomVecBacillus::data_atom_bonus(int m, char **values)
{
  if (bacillus[m])
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  // diagonalize inertia tensor

  double tensor[3][3];
  tensor[0][0] = atof(values[0]);
  tensor[1][1] = atof(values[1]);
  tensor[2][2] = atof(values[2]);
  tensor[0][1] = tensor[1][0] = atof(values[3]);
  tensor[0][2] = tensor[2][0] = atof(values[4]);
  tensor[1][2] = tensor[2][1] = atof(values[5]);

  double *inertia = bonus[nlocal_bonus].inertia;
  double evectors[3][3];
  int ierror = MathExtra::jacobi(tensor,inertia,evectors);
  if (ierror) error->one(FLERR,
                         "Insufficient Jacobi rotations for baccilus");

  // if any principal moment < scaled EPSILON, set to 0.0
  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON * max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON * max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON * max) inertia[2] = 0.0;
  //printf("i1=%e i2=%e i3=%e cut=%e \n", inertia[0],inertia[1],inertia[2],  EPSILON * max);

  // exyz_space = principal axes in space frame
  double ex_space[3],ey_space[3],ez_space[3];

  ex_space[0] = evectors[0][0];
  ex_space[1] = evectors[1][0];
  ex_space[2] = evectors[2][0];
  ey_space[0] = evectors[0][1];
  ey_space[1] = evectors[1][1];
  ey_space[2] = evectors[2][1];
  ez_space[0] = evectors[0][2];
  ez_space[1] = evectors[1][2];
  ez_space[2] = evectors[2][2];

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed
  double cross[3];
  MathExtra::cross3(ex_space,ey_space,cross);
  if (MathExtra::dot3(cross,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion
  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,bonus[nlocal_bonus].quat);

  double *pole1 = bonus[nlocal_bonus].pole1;
  double *pole2 = bonus[nlocal_bonus].pole2;
  double px = atof(values[6]);
  double py = atof(values[7]);
  double pz = atof(values[8]);

  pole1[0] = px;
  pole1[1] = py;
  pole1[2] = pz;

  double d = sqrt(px*px + py*py + pz*pz);

  bonus[nlocal_bonus].length = d * 2;

  pole2[0] = -px;
  pole2[1] = -py;
  pole2[2] = -pz;

  bonus[nlocal_bonus].diameter = atof(values[9]);
  if (bonus[nlocal_bonus].diameter < 0)
    error->one(FLERR, "Invalid diameter in Bacilli section of data file: diameter < 0");
  atom->radius[m] = bonus[nlocal_bonus].diameter * 0.5;

  // reset ellipsoid mass
  // previously stored density in rmass
  rmass[m] *= (4.0*MY_PI/3.0*
      atom->radius[m]*atom->radius[m]*atom->radius[m] +
      MY_PI*atom->radius[m]*atom->radius[m]*bonus[nlocal_bonus].length);
  biomass[m] = 1.0;

  bonus[nlocal_bonus].ilocal = m;
  bacillus[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   unpack one tri from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecBacillus::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  angmom[m][0] = atof(values[3]);
  angmom[m][1] = atof(values[4]);
  angmom[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one body in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecBacillus::data_vel_hybrid(int m, char **values)
{
  angmom[m][0] = atof(values[0]);
  angmom[m][1] = atof(values[1]);
  angmom[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBacillus::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    if (bacillus[i] < 0) buf[i][2] = ubuf(0).d;
    else buf[i][2] = ubuf(1).d;
    if (bacillus[i] < 0) buf[i][3] = rmass[i];
    else buf[i][3] = rmass[i] / (4.0*MY_PI/3.0*
	      atom->radius[i]*atom->radius[i]*atom->radius[i] +
	      MY_PI*atom->radius[i]*atom->radius[i]*bonus[bacillus[i]].length);
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = biomass[i];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_data_hybrid(int i, double *buf)
{
  if (bacillus[i] < 0) buf[0] = ubuf(0).d;
  else buf[0] = ubuf(1).d;
  if (bacillus[i] < 0) buf[1] = rmass[i];
  else buf[1] = rmass[i] / (4.0*MY_PI/3.0*
      atom->radius[i]*atom->radius[i]*atom->radius[i] +
      MY_PI*atom->radius[i]*atom->radius[i]*bonus[bacillus[i]].length);

  buf[2] = biomass[i];

  return 3;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBacillus::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,
            buf[i][3],buf[i][4],buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecBacillus::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d %-1.16e %-1.16e",(int) ubuf(buf[0]).i,buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecBacillus::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = angmom[i][0];
    buf[i][5] = angmom[i][1];
    buf[i][6] = angmom[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = angmom[i][0];
  buf[1] = angmom[i][1];
  buf[2] = angmom[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecBacillus::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecBacillus::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   reset quat orientation for atom M to quat_external
   called by Atom:add_molecule_atom()
------------------------------------------------------------------------- */

void AtomVecBacillus::set_quat(int m, double *quat_external)
{
  if (bacillus[m] < 0) error->one(FLERR,"Assigning quat to non-bacillus atom");
  double *quat = bonus[bacillus[m]].quat;
  quat[0] = quat_external[0]; quat[1] = quat_external[1];
  quat[2] = quat_external[2]; quat[3] = quat_external[3];
}

/* ----------------------------------------------------------------------
   set values in bonus data for bacillus m
------------------------------------------------------------------------- */

void AtomVecBacillus::set_bonus(int m, double *pole1, double diameter, double *quat, double *inertia)
{
  if (bacillus[m])
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  bonus[nlocal_bonus].pole1[0] = pole1[0];
  bonus[nlocal_bonus].pole1[1] = pole1[1];
  bonus[nlocal_bonus].pole1[2] = pole1[2];
  bonus[nlocal_bonus].pole2[0] = -pole1[0];
  bonus[nlocal_bonus].pole2[1] = -pole1[1];
  bonus[nlocal_bonus].pole2[2] = -pole1[2];

  bonus[nlocal_bonus].inertia[0] = inertia[0];
  bonus[nlocal_bonus].inertia[1] = inertia[1];
  bonus[nlocal_bonus].inertia[2] = inertia[2];

  bonus[nlocal_bonus].quat[0] = quat[0];
  bonus[nlocal_bonus].quat[1] = quat[1];
  bonus[nlocal_bonus].quat[2] = quat[2];
  bonus[nlocal_bonus].quat[3] = quat[3];

  double d = sqrt(pole1[0]*pole1[0] + pole1[1]*pole1[1] + pole1[2]*pole1[2]);

  bonus[nlocal_bonus].length = d * 2;
  bonus[nlocal_bonus].diameter = diameter;

  bonus[nlocal_bonus].ilocal = m;
  bacillus[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   get pole coordinate for bacillus m
------------------------------------------------------------------------- */
void AtomVecBacillus::get_pole_coords(int m, double *xp1, double *xp2)
{
  if (bacillus[m] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  double p[3][3];
  double *x;

  MathExtra::quat_to_mat(bonus[bacillus[m]].quat,p);
  MathExtra::matvec(p,bonus[bacillus[m]].pole1, xp1);
  MathExtra::matvec(p,bonus[bacillus[m]].pole2, xp2);

  x = atom->x[bonus[bacillus[m]].ilocal];

  xp1[0] += x[0];
  xp1[1] += x[1];
  xp1[2] += x[2];
  xp2[0] += x[0];
  xp2[1] += x[1];
  xp2[2] += x[2];
//
//  if (shift > 0) {
//    double dl = 2*shift/bonus[bacillus[m]].length;
//    xp1[0] += (xp1[0] - x[0]) * dl;
//    xp1[1] += (xp1[1] - x[1]) * dl;
//    xp1[2] += (xp1[2] - x[2]) * dl;
//    xp2[0] += (xp2[0] - x[0]) * dl;
//    xp2[1] += (xp2[1] - x[1]) * dl;
//    xp2[2] += (xp2[2] - x[2]) * dl;
//  }
}


/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecBacillus::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("biomass")) bytes += memory->usage(biomass,nmax);
  if (atom->memcheck("angmom")) bytes += memory->usage(angmom,nmax,3);
  if (atom->memcheck("torque")) bytes +=
                                  memory->usage(torque,nmax*comm->nthreads,3);
  if (atom->memcheck("bacillus")) bytes += memory->usage(bacillus,nmax);

  bytes += nmax_bonus*sizeof(Bonus);

  return bytes;
}
