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
#include "fix_monod_cyano_kokkos.h"
#include "atom_kokkos.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixMonodCyanoKokkos<DeviceType>::FixMonodCyanoKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMonodCyano(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixMonodCyanoKokkos<DeviceType>::compute()
{ 
  if (reaction_flag && growth_flag) {
    update_cells<1, 1>();
    update_atoms();
  } else if (reaction_flag && !growth_flag) {
    update_cells<1, 0>();
  } else if (!reaction_flag && growth_flag) {
    update_cells<0, 1>();
    update_atoms();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int Reaction, int Growth>
void FixMonodCyanoKokkos<DeviceType>::update_cells()
{
  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_conc = gridKK->k_conc.template view<DeviceType>();
  d_reac = gridKK->k_reac.template view<DeviceType>();
  d_dens = gridKK->k_dens.template view<DeviceType>();
  d_growth = gridKK->k_growth.template view<DeviceType>();

  if (Reaction)
    gridKK->sync(execution_space, GMASK_MASK | CONC_MASK | DENS_MASK);
  else
    gridKK->sync(execution_space, GMASK_MASK | CONC_MASK);

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixMonodCyanoCellsTag<Reaction, Growth> >(0, grid->ncells), f);
  copymode = 0;

  if (Growth)
    gridKK->modified(execution_space, GROWTH_MASK);
  if (Reaction)
    gridKK->modified(execution_space, REAC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixMonodCyanoKokkos<DeviceType>::update_atoms()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;
  double ***growth = grid->growth;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  gridKK->sync(Host, GROWTH_MASK);
  
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell = grid->cell(x[i]);
      const double density = rmass[i] /
	(four_thirds_pi * radius[i] * radius[i] * radius[i]);

      rmass[i] = rmass[i] * (1 + growth[igroup][cell][0] * dt);
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_mass[i] = 0;
      outer_radius[i] = radius[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixMonodCyanoKokkos<DeviceType>::Functor::Functor(FixMonodCyanoKokkos<DeviceType> *ptr):
  igroup(ptr->igroup),
  ilight(ptr->ilight), ico2(ptr->ico2), igco2(ptr->igco2),
  isuc(ptr->isuc), io2(ptr->io2), suc_exp(ptr->suc_exp), gco2_flag(ptr->gco2_flag),
  light_affinity(ptr->light_affinity), co2_affinity(ptr->co2_affinity),
  growth(ptr->growth), yield(ptr->yield),
  maintain(ptr->maintain), decay(ptr->decay),
  d_mask(ptr->d_mask), d_conc(ptr->d_conc), d_reac(ptr->d_reac),
  d_dens(ptr->d_dens), d_growth(ptr->d_growth) {}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int Reaction, int Growth>
KOKKOS_INLINE_FUNCTION
void FixMonodCyanoKokkos<DeviceType>::Functor::operator()(FixMonodCyanoCellsTag<Reaction, Growth>, int i) const
{
  double tmp1 = growth * d_conc(ilight, i) / (light_affinity + d_conc(ilight, i)) *
      d_conc(ico2, i) / (co2_affinity + d_conc(ico2, i));
  double tmp2 = 0.2 * tmp1 * suc_exp;
  double tmp3 = 4 * tmp1 * suc_exp;

  if (Reaction && !(d_mask(i) & GHOST_MASK)) {
    // nutrient utilization
    d_reac(ilight, i) -= 1 / yield * (tmp1 + tmp3) * d_dens(igroup, i);
    d_reac(ico2, i) -= 1 / yield * (tmp1 + tmp3) * d_dens(igroup, i);
    d_reac(io2, i) -= 0.1 * maintain * d_dens(igroup, i);
    // oxygen evolution
    d_reac(io2, i) += (0.727 / yield) * (tmp1 + tmp3) * d_dens(igroup, i);
    // sucrose export
    d_reac(isuc, i) += (0.65 / yield) * tmp3 * d_dens(igroup, i);

    // co2 dissolution
    if (gco2_flag == 1)
      d_reac(ico2, i) += 4.4e-6 * d_conc(igco2, i) - (4.4e-6 * d_conc(ico2, i));
  }

  if (Growth) {
    d_growth(igroup, i, 0) = tmp1 - tmp2 - decay - maintain;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixMonodCyanoKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixMonodCyanoKokkos<LMPHostType>;
#endif
}
