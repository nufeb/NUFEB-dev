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

#include "grid_vec_reactor.h"
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

GridVecReactor::GridVecReactor(LAMMPS *lmp) : GridVec(lmp)
{
  mask = NULL;
  conc = NULL;
  reac = NULL;
  dens = NULL;
  growth = NULL;
  bulk = NULL;
  grid->reactor_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GridVecReactor::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecReactor::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  if (n > nmax) {
    mask = memory->grow(grid->mask, n, "nufeb/reactor:mask");
    conc = memory->grow(grid->conc, grid->nsubs, n, "nufeb/reactor:conc");
    reac = memory->grow(grid->reac, grid->nsubs, n, "nufeb/reactor:reac");
    dens = memory->grow(grid->dens, group->ngroup, n, "nufeb/reactor:dens");
    growth = memory->grow(grid->growth, group->ngroup, n, 2, "nufeb/reactor:grow");
    bulk  = memory->grow(grid->bulk, grid->nsubs, "nufeb/reactor:bulk");
    nmax = n;
    grid->nmax = nmax;

    grid->mask = mask;
    grid->conc = conc;
    grid->reac = reac;
    grid->dens = dens;
    grid->growth = growth;
    grid->bulk = bulk;
  }
}

/* ---------------------------------------------------------------------- */

int GridVecReactor::pack_comm(int n, int *cells, double *buf)
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

void GridVecReactor::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int GridVecReactor::pack_exchange(int n, int *cells, double *buf)
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

void GridVecReactor::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridVecReactor::set(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Invalid grid_modify set command");
  int isub = grid->find(arg[1]);
  if (isub < 0) error->all(FLERR,"Cannot find substrate name");
  set_reactor(isub, force->numeric(FLERR, arg[2]), force->numeric(FLERR, arg[3]));
}

/* ---------------------------------------------------------------------- */

void GridVecReactor::set_reactor(int isub, double domain, double cbulk)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK)) {
      conc[isub][i] = domain;
    }
    grid->reac[isub][i] = 0.0;
  }
  bulk[isub] = cbulk;
}
