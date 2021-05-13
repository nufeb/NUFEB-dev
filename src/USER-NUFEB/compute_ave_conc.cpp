/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_ave_conc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"

#include "modify.h"
#include "grid.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeAveConc::ComputeAveConc(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute avgcon command");

  vector_flag = 1;
  extvector = 0;

  size_vector = grid->nsubs;
  memory->create(vector,size_vector,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputeAveConc::~ComputeAveConc()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeAveConc::compute_vector()
{
  invoked_vector = update->ntimestep;

  for(int isub = 0; isub < size_vector; isub++){
    double sum = 0;
    for (int i = 0; i < grid->ncells; i++) {
      if (!(grid->mask[i] & GHOST_MASK)) {
	sum += grid->conc[isub][i];
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    sum /= (grid->box[0] * grid->box[1] * grid->box[2]);
    vector[isub] = sum;
  }
}
