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

#include <cstdio>
#include <cstring>
#include <cmath>
#include "fix_diffusion_reaction.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "grid.h"
#include "grid_masks.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DIRICHLET,NEUMANN,PERIODIC,BULK};

/* ---------------------------------------------------------------------- */

FixDiffusionReaction::FixDiffusionReaction(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7)
    error->all(FLERR,"Illegal fix nufeb/diffusion_reaction command");

  compute_flag = 1;
  dynamic_group_allow = 1;
  scalar_flag = 1;

  closed_system = 0;

  ncells = 0;
  dt = 1.0;
  prev = nullptr;
  penult = nullptr;
  
  boundary[0] = boundary[1] = boundary[2] = boundary[3] =
  boundary[4] = boundary[5] = -1;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate for nufeb/diffusion_reaction");

  diff_coef = utils::numeric(FLERR,arg[4],true,lmp);

  int ndirichlet = 0;
  int nbulk = 0;

  for (int i = 0; i < 3; i++) {
    if ((arg[5+i][0] == 'p' && arg[5+i][1] != 'p') ||
	(arg[5+i][1] == 'p' && arg[5+i][0] != 'p'))
      error->all(FLERR, "Illegal boundary condition");

    for (int j = 0; j < 2; j++) {
      if (arg[5+i][j] == 'p') {
	boundary[2*i+j] = PERIODIC;
	grid->periodic[i] = 1;
      } else if (arg[5+i][j] == 'n') {
	boundary[2*i+j] = NEUMANN;
      } else if (arg[5+i][j] == 'd') {
	boundary[2*i+j] = DIRICHLET;
	ndirichlet++;
      } else if (arg[5+i][j] == 'b') {
	if (!grid->reactor_flag)
	  error->all(FLERR, "Illegal boundary condition");
	boundary[2*i+j] = BULK;
	nbulk++;
      } else {
	error->all(FLERR, "Illegal boundary condition");
      }
    }
  }

  if (!ndirichlet && !nbulk) closed_system = 1;

  if (narg < ndirichlet + 7)
    error->all(FLERR, "Not enough values for dirichlet boundaries");
  int iarg = 8;
  for (int i = 0; i < 6; i++) {
    if (boundary[i] == DIRICHLET) {
      dirichlet[i] = utils::numeric(FLERR,arg[iarg++],true,lmp);
    } else {
      dirichlet[i] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixDiffusionReaction::~FixDiffusionReaction()
{

  if (copymode) return;
  memory->destroy(prev);
  if (closed_system) memory->destroy(penult);
}

/* ---------------------------------------------------------------------- */

int FixDiffusionReaction::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixDiffusionReaction::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compute") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	compute_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	compute_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::init()
{
  ncells = grid->ncells;
  prev = memory->create(prev, ncells, "nufeb/diffusion_reaction:prev");
  for (int i = 0; i < ncells; i++)
    prev[i] = 0.0;
  if (closed_system) {
    penult = memory->create(penult, ncells, "nufeb/diffusion_reaction:penult");
    for (int i = 0; i < ncells; i++)
      penult[i] = 0.0;
  }
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::pre_force(int)
{
  if (compute_flag)
    compute_initial();
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::final_integrate()
{
  if (compute_flag)
    compute_final();
}

/* ---------------------------------------------------------------------- */

double FixDiffusionReaction::compute_scalar()
{
  double result = 0.0;
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      double res = fabs((grid->conc[isub][i] - prev[i]) / prev[i]);
      if (closed_system) {
        double res2 = fabs((prev[i] - penult[i]) / penult[i]);
        res = fabs(res - res2);
      }
      result = MAX(result, res);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, world);
  return result;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::compute_initial()
{
  if (ncells != grid->ncells) {
    ncells = grid->ncells;
    prev = memory->grow(prev, ncells, "nufeb/diffusion_reaction:prev");
    if (closed_system) penult = memory->grow(penult, ncells, "nufeb/diffusion_reaction:penult");
  }

  for (int i = 0; i < grid->ncells; i++) {
    // Dirichlet boundary conditions
    if (grid->mask[i] & X_NB_MASK && boundary[0] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[0];
    } else if (grid->mask[i] & X_PB_MASK && boundary[1] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[1];
    } else if (grid->mask[i] & Y_NB_MASK && boundary[2] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[2];
    } else if (grid->mask[i] & Y_PB_MASK && boundary[3] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[3];
    } else if (grid->mask[i] & Z_NB_MASK && boundary[4] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[4];
    } else if (grid->mask[i] & Z_PB_MASK && boundary[5] == DIRICHLET) {
      grid->conc[isub][i] = dirichlet[5];
    } else if (grid->mask[i] & X_NB_MASK && boundary[0] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & X_PB_MASK && boundary[1] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Y_NB_MASK && boundary[2] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Y_PB_MASK && boundary[3] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Z_NB_MASK && boundary[4] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Z_PB_MASK && boundary[5] == BULK) {
      grid->conc[isub][i] = grid->bulk[isub];
    }
    grid->reac[isub][i] = 0.0;
  }
}



/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::compute_final()
{
  int nx = grid->subbox[0];
  int nxy = grid->subbox[0] * grid->subbox[1];
  for (int i = 0; i < grid->ncells; i++) {
    // Neumann boundary conditions
    if (grid->mask[i] & X_NB_MASK && boundary[0] == NEUMANN) {
      grid->conc[isub][i] = grid->conc[isub][i+1];
    } else if (grid->mask[i] & X_PB_MASK && boundary[1] == NEUMANN) {
      grid->conc[isub][i] = grid->conc[isub][i-1];
    } else if (grid->mask[i] & Y_NB_MASK && boundary[2] == NEUMANN) {
      int py = i + nx;
      grid->conc[isub][i] = grid->conc[isub][py];
    } else if (grid->mask[i] & Y_PB_MASK && boundary[3] == NEUMANN) {
      int py = i - nx;
      grid->conc[isub][i] = grid->conc[isub][py];
    } else if (grid->mask[i] & Z_NB_MASK && boundary[4] == NEUMANN) {
      int pz = i + nxy;
      grid->conc[isub][i] = grid->conc[isub][pz];
    } else if (grid->mask[i] & Z_PB_MASK && boundary[5] == NEUMANN) {
      int pz = i - nxy;
      grid->conc[isub][i] = grid->conc[isub][pz];
    }
    if (closed_system) penult[i] = prev[i];
    prev[i] = grid->conc[isub][i];
  }

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      int nx = i - 1;
      int px = i + 1;
      int ny = i - grid->subbox[0];
      int py = i + grid->subbox[0];
      int nz = i - nxy;
      int pz = i + nxy;
      double dnx = diff_coef * (prev[i] - prev[nx]) / grid->cell_size;
      double dpx = diff_coef * (prev[px] - prev[i]) / grid->cell_size;
      double ddx = (dpx - dnx) / grid->cell_size;
      double dny = diff_coef * (prev[i] - prev[ny]) / grid->cell_size;
      double dpy = diff_coef * (prev[py] - prev[i]) / grid->cell_size;
      double ddy = (dpy - dny) / grid->cell_size;
      double dnz = diff_coef * (prev[i] - prev[nz]) / grid->cell_size;
      double dpz = diff_coef * (prev[pz] - prev[i]) / grid->cell_size;
      double ddz = (dpz - dnz) / grid->cell_size;
      // prevent negative concentrations
      grid->conc[isub][i] = MAX(0, prev[i] + dt * (ddx + ddy + ddz + grid->reac[isub][i]));
    }
  }
}

/* ----------------------------------------------------------------------
 Average substrate distribution before solving diffusion in closed system.
 ------------------------------------------------------------------------- */
void FixDiffusionReaction::closed_system_init()
{
  if (!closed_system) return;
  double sum = 0;
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      sum += grid->conc[isub][i];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
  sum /= (grid->box[0] * grid->box[1] * grid->box[2]);
  for (int i = 0; i < grid->ncells; i++) {
    grid->conc[isub][i] = sum;
  }
}

/* ----------------------------------------------------------------------
 Scaleup substrate concentrations in closed system based on residual value and
 biological timestep
 ------------------------------------------------------------------------- */
void FixDiffusionReaction::closed_system_scaleup(double biodt)
{
  if (!closed_system) return;
  for (int i = 0; i < grid->ncells; i++) {
    double res = grid->conc[isub][i] - prev[i];
    grid->conc[isub][i] += res / dt * biodt;
    grid->conc[isub][i] = MAX(0, grid->conc[isub][i]);
  }
}
