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
#include "fix_diffusion_reaction_kokkos.h"
#include "error.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DIRICHLET,NEUMANN,PERIODIC};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDiffusionReactionKokkos<DeviceType>::FixDiffusionReactionKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDiffusionReaction(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  closed_system = FixDiffusionReaction::closed_system;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDiffusionReactionKokkos<DeviceType>::~FixDiffusionReactionKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(prev);
  if (closed_system) memoryKK->destroy_kokkos(penult);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::init()
{
  FixDiffusionReaction::init();
  memory->destroy(prev);
  if (closed_system) memoryKK->destroy_kokkos(penult);
  grow_arrays(ncells);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_prev.template sync<LMPHostType>();
  memoryKK->grow_kokkos(k_prev, prev, nmax, "nufeb/diffusion_reaction:prev");
  d_prev = k_prev.template view<DeviceType>();
  k_prev.template modify<LMPHostType>();

  if (closed_system) {
    k_penult.template sync<LMPHostType>();
    memoryKK->grow_kokkos(k_penult, penult, nmax, "nufeb/diffusion_reaction:penult");
    d_penult = k_penult.template view<DeviceType>();
    k_penult.template modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double FixDiffusionReactionKokkos<DeviceType>::compute_scalar()
{
  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_conc = gridKK->k_conc.template view<DeviceType>();

  gridKK->sync(execution_space, GMASK_MASK);
  gridKK->sync(execution_space, CONC_MASK);

  copymode = 1;
  Functor f(this);
  double result = 0.0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<
    DeviceType,
    FixDiffusionReactionScalarTag>(0, grid->ncells), f, Kokkos::Max<double>(result));
  DeviceType().fence();
  copymode = 0;

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, world);
  return result;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::compute_initial()
{
  for (int i = 0; i < 6; i++) {
    boundary[i] = FixDiffusionReaction::boundary[i];
    dirichlet[i] = FixDiffusionReaction::dirichlet[i];
  }

  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_conc = gridKK->k_conc.template view<DeviceType>();
  d_reac = gridKK->k_reac.template view<DeviceType>();
  
  if (ncells < grid->ncells) {
    grow_arrays(grid->ncells);
    ncells = grid->ncells;
  }

  gridKK->sync(execution_space, GMASK_MASK | CONC_MASK);
  
  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDiffusionReactionInitialTag>(0, grid->ncells), f);
  copymode = 0;
  
  gridKK->modified(execution_space, CONC_MASK | REAC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::compute_final()
{
  cell_size = grid->cell_size;
  for (int i = 0; i < 3; i++)
    subbox[i] = grid->subbox[i];

  for (int i = 0; i < 6; i++)
    boundary[i] = FixDiffusionReaction::boundary[i];

  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_conc = gridKK->k_conc.template view<DeviceType>();
  d_reac = gridKK->k_reac.template view<DeviceType>();
  
  gridKK->sync(execution_space, GMASK_MASK);
  gridKK->sync(execution_space, REAC_MASK);
  
  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDiffusionReactionPreFinalTag>(0, grid->ncells), f);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDiffusionReactionFinalTag>(0, grid->ncells), f);
  copymode = 0;

  gridKK->modified(execution_space, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::closed_system_init()
{
  if (!closed_system) return;
  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_conc = gridKK->k_conc.template view<DeviceType>();

  gridKK->sync(execution_space, GMASK_MASK);
  gridKK->sync(execution_space, CONC_MASK);

  double result = 0.0;
  Kokkos::parallel_reduce(
    grid->ncells, LAMMPS_LAMBDA(int i, double& sum) {
      if (!(d_mask(i) & GHOST_MASK))
	sum += d_conc(isub, i);
    }, result);
  DeviceType().fence();

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);

  result /= (grid->box[0] * grid->box[1] * grid->box[2]);

  Kokkos::parallel_for(
    grid->ncells, LAMMPS_LAMBDA(int i) {
      d_conc(isub, i) = result;
  });

  gridKK->modified(execution_space, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDiffusionReactionKokkos<DeviceType>::closed_system_scaleup(double biodt)
{
  if (!closed_system) return;

  d_conc = gridKK->k_conc.template view<DeviceType>();

  gridKK->sync(execution_space, CONC_MASK);

  Kokkos::parallel_for(
    grid->ncells, LAMMPS_LAMBDA(int i) {
    double res = d_conc(isub, i) - d_prev(i);
    d_conc(isub, i) += res / dt * biodt;
    d_conc(isub, i) = MAX(0, d_conc(isub, i));
  });

  gridKK->modified(execution_space, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDiffusionReactionKokkos<DeviceType>::Functor::Functor(FixDiffusionReactionKokkos<DeviceType> *ptr):
  d_prev(ptr->d_prev), d_conc(ptr->d_conc),
  d_reac(ptr->d_reac), d_mask(ptr->d_mask), d_penult(ptr->d_penult),
  cell_size(ptr->cell_size), diff_coef(ptr->diff_coef), dt(ptr->dt),
  isub(ptr->isub), closed_system(ptr->closed_system)
{
  for (int i = 0; i < 3; i++)
    subbox[i] = ptr->subbox[i];
  for (int i = 0; i < 6; i++) {
    boundary[i] = ptr->boundary[i];
    dirichlet[i] = ptr->dirichlet[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDiffusionReactionKokkos<DeviceType>::Functor::operator()(FixDiffusionReactionInitialTag, int i) const
{
  // Dirichlet boundary conditions
  if (d_mask(i) & X_NB_MASK && boundary[0] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[0];
  } else if (d_mask(i) & X_PB_MASK && boundary[1] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[1];
  } else if (d_mask(i) & Y_NB_MASK && boundary[2] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[2];
  } else if (d_mask(i) & Y_PB_MASK && boundary[3] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[3];
  } else if (d_mask(i) & Z_NB_MASK && boundary[4] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[4];
  } else if (d_mask(i) & Z_PB_MASK && boundary[5] == DIRICHLET) {
    d_conc(isub, i) = dirichlet[5];
  }
  d_reac(isub, i) = 0.0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDiffusionReactionKokkos<DeviceType>::Functor::operator()(FixDiffusionReactionPreFinalTag, int i) const
{
  int nx = subbox[0];
  int nxy = subbox[0] * subbox[1];
  // Neumann boundary conditions
  if (d_mask(i) & X_NB_MASK && boundary[0] == NEUMANN) {
    d_conc(isub, i) = d_conc(isub, i+1);
  } else if (d_mask(i) & X_PB_MASK && boundary[1] == NEUMANN) {
    d_conc(isub, i) = d_conc(isub, i-1);
  } else if (d_mask(i) & Y_NB_MASK && boundary[2] == NEUMANN) {
    int py = i + nx;
    d_conc(isub, i) = d_conc(isub, py);
  } else if (d_mask(i) & Y_PB_MASK && boundary[3] == NEUMANN) {
    int py = i - nx;
    d_conc(isub, i) = d_conc(isub, py);
  } else if (d_mask(i) & Z_NB_MASK && boundary[4] == NEUMANN) {
    int pz = i + nxy;
    d_conc(isub, i) = d_conc(isub, pz);
  } else if (d_mask(i) & Z_PB_MASK && boundary[5] == NEUMANN) {
    int pz = i - nxy;
    d_conc(isub, i) = d_conc(isub, pz);
  }
  if (closed_system) d_penult(i) = d_prev(i);
  d_prev(i) = d_conc(isub, i);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDiffusionReactionKokkos<DeviceType>::Functor::operator()(FixDiffusionReactionFinalTag, int i) const
{
  int nxy = subbox[0] * subbox[1];
  if (!(d_mask(i) & GHOST_MASK)) {
    int nx = i - 1;
    int px = i + 1;
    int ny = i - subbox[0];
    int py = i + subbox[0];
    int nz = i - nxy;
    int pz = i + nxy;
    double dnx = diff_coef * (d_prev(i) - d_prev(nx)) / cell_size;
    double dpx = diff_coef * (d_prev(px) - d_prev(i)) / cell_size;
    double ddx = (dpx - dnx) / cell_size;
    double dny = diff_coef * (d_prev(i) - d_prev(ny)) / cell_size;
    double dpy = diff_coef * (d_prev(py) - d_prev(i)) / cell_size;
    double ddy = (dpy - dny) / cell_size;
    double dnz = diff_coef * (d_prev(i) - d_prev(nz)) / cell_size;
    double dpz = diff_coef * (d_prev(pz) - d_prev(i)) / cell_size;
    double ddz = (dpz - dnz) / cell_size;
    // prevent negative concentrations
    d_conc(isub, i) = MAX(0, d_prev(i) + dt * (ddx + ddy + ddz + d_reac(isub, i)));
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDiffusionReactionKokkos<DeviceType>::Functor::operator()(FixDiffusionReactionScalarTag, int i, double &max) const
{
  if (!(d_mask(i) & GHOST_MASK)) {
    double res = fabs((d_conc(isub, i) - d_prev(i)) / d_prev(i));
    if (closed_system) {
      double res2 = fabs((d_prev(i) - d_penult(i)) / d_penult(i));
      res = fabs(res - res2);
    }
    max = MAX(max, res);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDiffusionReactionKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDiffusionReactionKokkos<LMPHostType>;
#endif
}
