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

#include "grid_vec_monod_kokkos.h"
#include "grid_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "grid_masks.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "group.h"
#include "lammps.h"
#include "kokkos.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecMonodKokkos::GridVecMonodKokkos(LAMMPS *lmp) : GridVecKokkos(lmp)
{
  mask = NULL;
  conc = NULL;
  reac = NULL;
  dens = NULL;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");
  
  if (n > nmax) {
    sync(Device, ALL_MASK);
    modified(Device, ALL_MASK);

    memoryKK->grow_kokkos(gridKK->k_mask, gridKK->mask, n, "nufeb/monod:mask");
    memoryKK->grow_kokkos(gridKK->k_conc, gridKK->conc, grid->nsubs, n, "nufeb/monod:conc");
    memoryKK->grow_kokkos(gridKK->k_reac, gridKK->reac, grid->nsubs, n, "nufeb/monod:reac");
    memoryKK->grow_kokkos(gridKK->k_dens, gridKK->dens, group->ngroup, n, "nufeb/monod:dens");
    memoryKK->grow_kokkos(gridKK->k_growth, gridKK->growth, group->ngroup, n, 2, "nufeb/monod:grow");
    nmax = n;
    grid->nmax = nmax;

    mask = gridKK->mask;
    d_mask = gridKK->k_mask.d_view;
    h_mask = gridKK->k_mask.h_view;
    conc = gridKK->conc;
    d_conc = gridKK->k_conc.d_view;
    h_conc = gridKK->k_conc.h_view;
    reac = gridKK->reac;
    d_reac = gridKK->k_reac.d_view;
    h_reac = gridKK->k_reac.h_view;
    dens = gridKK->dens;
    d_dens = gridKK->k_dens.d_view;
    h_dens = gridKK->k_dens.h_view;
    growth = gridKK->growth;
    d_growth = gridKK->k_growth.d_view;
    h_growth = gridKK->k_growth.h_view;

    sync(Host, ALL_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_comm(int n, int *cells, double *buf)
{
  sync(Host, CONC_MASK);
  
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_exchange(int n, int *cells, double *buf)
{
  sync(Host, CONC_MASK);
  
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct GridVecMonodKokkos_PackComm
{
  int _nsubs;
  
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_float_2d _conc;
  typename AT::t_int_1d _list;
  typename AT::t_xfloat_2d _buf;
  
  GridVecMonodKokkos_PackComm(
    int nsubs,
    const typename DAT::tdual_float_2d &conc,
    const DAT::tdual_int_1d &list,
    const DAT::tdual_xfloat_1d &buf):
    _nsubs(nsubs),
    _conc(conc.view<DeviceType>()),
    _list(list.view<DeviceType>())
  {
    _buf = typename AT::t_xfloat_2d_um(
      buf.view<DeviceType>().data(),
      buf.view<DeviceType>().extent(0)/nsubs,
      nsubs);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    int cell = _list(i);
    for (int j = 0; j < _nsubs; j++) {
      _buf(i, j) = _conc(j, cell);
    }
  }
};

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_comm_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{
  if (lmp->kokkos->forward_comm_on_host) {
    gridKK->sync(Host, CONC_MASK);
    struct GridVecMonodKokkos_PackComm<LMPHostType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPHostType>(first, last), f);
  } else {
    gridKK->sync(Device, CONC_MASK);
    struct GridVecMonodKokkos_PackComm<LMPDeviceType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
  }
  return (last-first)*grid->nsubs;
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
struct GridVecMonodKokkos_UnpackComm
{
  int _nsubs;
  
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_float_2d _conc;
  typename AT::t_int_1d _list;
  typename AT::t_xfloat_2d _buf;
  
  GridVecMonodKokkos_UnpackComm(
    int nsubs,
    const typename DAT::tdual_float_2d &conc,
    const DAT::tdual_int_1d &list,
    const DAT::tdual_xfloat_1d &buf):
    _nsubs(nsubs),
    _conc(conc.view<DeviceType>()),
    _list(list.view<DeviceType>())
  {
    _buf = typename AT::t_xfloat_2d_um(
      buf.view<DeviceType>().data(),
      buf.view<DeviceType>().extent(0)/nsubs,
      nsubs);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    int cell = _list(i);
    for (int j = 0; j < _nsubs; j++) {
      _conc(j, cell) = _buf(i, j);
    }
  }
};

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_comm_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{
  if (lmp->kokkos->forward_comm_on_host) {
    struct GridVecMonodKokkos_UnpackComm<LMPHostType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
    gridKK->modified(Host, CONC_MASK);
  } else {
    gridKK->sync(Device, CONC_MASK);
    struct GridVecMonodKokkos_UnpackComm<LMPDeviceType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
    gridKK->modified(Device, CONC_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_exchange_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_exchange_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::set(int narg, char **arg)
{
  sync(Host, GMASK_MASK);
  sync(Host, CONC_MASK);

  int isub = grid->find(arg[1]);
  if (isub < 0) error->all(FLERR,"Cannot find substrate name");
  if (narg == 3) set_monod(isub, force->numeric(FLERR, arg[2]));
  else set_monod(isub, force->numeric(FLERR, arg[2]),
		force->numeric(FLERR, arg[3]), force->numeric(FLERR, arg[4]),
		force->numeric(FLERR, arg[5]), force->numeric(FLERR, arg[6]),
		force->numeric(FLERR, arg[7]), force->numeric(FLERR, arg[8]));

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::set_monod(int sub, double domain)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK))
      conc[sub][i] = domain;
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::set_monod(int isub, double domain, double nx, double px,
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

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if ((mask & GMASK_MASK) && gridKK->k_mask.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(gridKK->k_mask, space);
    if ((mask & CONC_MASK) && gridKK->k_conc.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_conc, space);
    if ((mask & REAC_MASK) && gridKK->k_reac.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_reac, space);
    if ((mask & DENS_MASK) && gridKK->k_dens.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_dens, space);
    if ((mask & GROWTH_MASK) && gridKK->k_growth.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_3d>(gridKK->k_growth, space);
  } else {
    if ((mask & GMASK_MASK) && gridKK->k_mask.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(gridKK->k_mask, space);
    if ((mask & CONC_MASK) && gridKK->k_conc.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_conc, space);
    if ((mask & REAC_MASK) && gridKK->k_reac.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_reac, space);
    if ((mask & DENS_MASK) && gridKK->k_dens.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_dens, space);
    if ((mask & GROWTH_MASK) && gridKK->k_growth.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_3d>(gridKK->k_growth, space);
  }
}
