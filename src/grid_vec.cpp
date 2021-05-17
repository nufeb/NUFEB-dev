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

#include <cstring>
#include <cmath>
#include "grid_vec.h"
#include "grid.h"
#include "domain.h"
#include "error.h"
#include "grid_masks.h"
#include "comm.h"

using namespace LAMMPS_NS;

enum{DIRICHLET,NEUMANN,PERIODIC};

/* ---------------------------------------------------------------------- */

GridVec::GridVec(LAMMPS *lmp) : Pointers(lmp)
{
  kokkosable = 0;
  nargcopy = 0;
  argcopy = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

void GridVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  argcopy = new char*[nargcopy];
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(arg[i]) + 1;
    argcopy[i] = new char[n];
    strcpy(argcopy[i],arg[i]);
  }
}

/* ---------------------------------------------------------------------- */

void GridVec::process_args(int narg, char **arg)
{
  grid->nsubs = utils::inumeric(FLERR,arg[0],true,lmp);
  if (narg < grid->nsubs + 1)
    error->all(FLERR,"Missing substrate names in grid_style command");
  grid->sub_names = new char*[grid->nsubs];
  for (int i = 0; i < grid->nsubs; i++) {
    grid->sub_names[i] = new char[strlen(arg[i+1])+1];
    strcpy(grid->sub_names[i], arg[i+1]);
  }
  if (narg < grid->nsubs + 2)
    error->all(FLERR,"Missing cell size in grid_style command");
  grid->cell_size = utils::numeric(FLERR,arg[grid->nsubs+1],true,lmp);
}

/* ---------------------------------------------------------------------- */

void GridVec::init()
{
  if (lmp->kokkos != NULL && !kokkosable)
    error->all(FLERR,"KOKKOS package requires a kokkos enabled grid_style");

  const double small = 1e-12;
  grid->box[0] = static_cast<int>(domain->prd[0] / grid->cell_size + small);
  grid->box[1] = static_cast<int>(domain->prd[1] / grid->cell_size + small);
  grid->box[2] = static_cast<int>(domain->prd[2] / grid->cell_size + small);

  // check for incompatible sizes
  if (fabs(grid->cell_size * grid->box[0] - domain->prd[0]) > small)
    error->all(FLERR,"Grid cell size incompatible with simulation box x size.");
  if (fabs(grid->cell_size * grid->box[1] - domain->prd[1]) > small)
    error->all(FLERR,"Grid cell size incompatible with simulation box y size.");
  if (fabs(grid->cell_size * grid->box[2] - domain->prd[2]) > small)
    error->all(FLERR,"Grid cell size incompatible with simulation box z size.");

  // extend global grid size
  for (int i = 0; i < 3; i++)
    grid->extbox[i] = grid->box[i] + 2;

  // Fitting initial domain decomposition to the grid
  for (int i = 0; i < comm->procgrid[0]; i++) {
    int n = grid->box[0] * i * 1.0 / comm->procgrid[0];
    comm->xsplit[i] = (double) n / grid->box[0];
  }
  for (int i = 0; i < comm->procgrid[1]; i++) {
    int n = grid->box[1] * i * 1.0 / comm->procgrid[1];
    comm->ysplit[i] = (double) n / grid->box[1];
  }
  for (int i = 0; i < comm->procgrid[2]; i++) {
    int n = grid->box[2] * i * 1.0 / comm->procgrid[2];
    comm->zsplit[i] = (double) n / grid->box[2];
  }
  domain->set_local_box();
}

/* ---------------------------------------------------------------------- */

void GridVec::setup()
{
  const double small = 1e-12;
  for (int i = 0; i < 3; i++) {
    grid->sublo[i] = static_cast<int>((domain->sublo[i] - domain->boxlo[i]) /
				      grid->cell_size + small) - 1;
    grid->subhi[i] = static_cast<int>((domain->subhi[i] - domain->boxlo[i]) /
				      grid->cell_size + small) + 1;
    grid->subbox[i] = grid->subhi[i] - grid->sublo[i];
  }
  grid->ncells = grid->subbox[0] * grid->subbox[1] * grid->subbox[2];

  if (grid->ncells > grid->nmax) {
    grow(grid->ncells);
  }

  // setup mask
  int *mask = grid->mask;
  for (int z = 0; z < grid->subbox[2]; z++) {
    for (int y = 0; y < grid->subbox[1]; y++) {
      for (int x = 0; x < grid->subbox[0]; x++) {
	int m = 0;
	if (x == 0 || y == 0 || z == 0 ||
	    x == grid->subbox[0] - 1 ||
	    y == grid->subbox[1] - 1 ||
	    z == grid->subbox[2] - 1)
	  m |= GHOST_MASK;
	if ((x == 0 && y == 0 &&
	     grid->sublo[0] < 0 && grid->sublo[1] < 0) ||
	    (x == 0 && z == 0 &&
	     grid->sublo[0] < 0 && grid->sublo[2] < 0) ||
	    (x == 0 && y == grid->subbox[1] - 1 &&
	     grid->sublo[0] < 0 && grid->subhi[1] > grid->box[1]) ||
	    (x == 0 && z == grid->subbox[2] - 1 &&
	     grid->sublo[0] < 0 && grid->subhi[2] > grid->box[2]) ||
	    (x == grid->subbox[0] - 1 && y == 0 &&
	     grid->subhi[0] > grid->box[0] && grid->sublo[1] < 0) ||
	    (x == grid->subbox[0] - 1 && z == 0 &&
	     grid->subhi[0] > grid->box[0] && grid->sublo[2] < 0) ||
	    (x == grid->subbox[0] - 1 && y == grid->subbox[1] - 1 &&
	     grid->subhi[0] > grid->box[0] && grid->subhi[1] > grid->box[1]) ||
	    (x == grid->subbox[0] - 1 && z == grid->subbox[2] - 1 &&
	     grid->subhi[0] > grid->box[0] && grid->subhi[2] > grid->box[2]) ||
	    (y == 0 && z == 0 &&
	     grid->sublo[1] < 0 && grid->sublo[2] < 0) ||
	    (y == 0 && z == grid->subbox[2] - 1 &&
	     grid->sublo[1] < 0 && grid->subhi[2] > grid->box[2]) ||
	    (y == grid->subbox[1] - 1 && z == 0 &&
	     grid->subhi[1] > grid->box[1] && grid->sublo[2] < 0) ||
	    (y == grid->subbox[1] - 1 && z == grid->subbox[2] - 1 &&
	     grid->subhi[1] > grid->box[1] && grid->subhi[2] > grid->box[2]))
	  m |= CORNER_MASK;
	else {
	  if (grid->sublo[0] < 0 && x == 0)
	    m |= X_NB_MASK;
	  if (grid->subhi[0] > grid->box[0] && x == grid->subbox[0] - 1)
	    m |= X_PB_MASK;
	  if (grid->sublo[1] < 0 && y == 0)
	    m |= Y_NB_MASK;
	  if (grid->subhi[1] > grid->box[1] && y == grid->subbox[1] - 1)
	    m |= Y_PB_MASK;
	  if (grid->sublo[2] < 0 && z == 0)
	    m |= Z_NB_MASK;
	  if (grid->subhi[2] > grid->box[2] && z == grid->subbox[2] - 1)
	    m |= Z_PB_MASK;
	}
	mask[x + y * grid->subbox[0] + z * grid->subbox[0] * grid->subbox[1]] = m;
      }
    }
  }
}
