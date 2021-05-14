/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "grid_vec_monod.h"
#include "grid.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include "grid_masks.h"
#include "comm.h"
#include "atom.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecMonod::GridVecMonod(LAMMPS *lmp) : GridVec(lmp)
{
  mask = NULL;
  conc = NULL;
  reac = NULL;
  dens = NULL;
  growth = NULL;
  grid->monod_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  if (n > nmax) {
    mask = memory->grow(grid->mask, n, "nufeb/monod:mask");
    conc = memory->grow(grid->conc, grid->nsubs, n, "nufeb/monod:conc");
    reac = memory->grow(grid->reac, grid->nsubs, n, "nufeb/monod:reac");
    dens = memory->grow(grid->dens, group->ngroup, n, "nufeb/monod:dens");
    growth = memory->grow(grid->growth, group->ngroup, n, 2, "nufeb/monod:grow");
    nmax = n;
    grid->nmax = nmax;

    grid->mask = mask;
    grid->conc = conc;
    grid->reac = reac;
    grid->dens = dens;
    grid->growth = growth;
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonod::pack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonod::pack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::set(int narg, char **arg)
{
  if (narg != 3 && narg != 9) error->all(FLERR, "Invalid grid_modify set command");
  int isub = grid->find(arg[1]);
  if (isub < 0) error->all(FLERR,"Cannot find substrate name");
  if (narg == 3) set_monod(isub, force->numeric(FLERR, arg[2]));
  else set_monod(isub, force->numeric(FLERR, arg[2]),
		force->numeric(FLERR, arg[3]), force->numeric(FLERR, arg[4]),
		force->numeric(FLERR, arg[5]), force->numeric(FLERR, arg[6]),
		force->numeric(FLERR, arg[7]), force->numeric(FLERR, arg[8]));
}


/* ---------------------------------------------------------------------- */

void GridVecMonod::set_monod(int sub, double domain)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK))
      conc[sub][i] = domain;
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::set_monod(int isub, double domain, double nx, double px,
		       double ny, double py, double nz, double pz)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK)) {
      if (mask[i] & X_NB_MASK) {
	conc[isub][i] = nx;
      } else if (mask[i] & X_PB_MASK) {
	conc[isub][i] = px;
      } else if (mask[i] & Y_NB_MASK) {
	conc[isub][i] = ny;
      } else if (mask[i] & Y_PB_MASK) {
	conc[isub][i] = py;
      } else if (mask[i] & Z_NB_MASK) {
	conc[isub][i] = nz;
      } else if (mask[i] & Z_PB_MASK) {
	conc[isub][i] = pz;
      } else {
	conc[isub][i] = domain;
      }
    }
    grid->reac[isub][i] = 0.0;
  }
}
