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

#ifdef GRID_CLASS

GridStyle(nufeb/monod/kk,GridVecMonodKokkos)
GridStyle(nufeb/monod/kk/device,GridVecMonodKokkos)
GridStyle(nufeb/monod/kk/host,GridVecMonodKokkos)

#else

#ifndef LMP_GRID_VEC_MONOD_KOKKOS_H
#define LMP_GRID_VEC_MONOD_KOKKOS_H

#include "grid_vec_kokkos.h"

namespace LAMMPS_NS {

class GridVecMonodKokkos : public GridVecKokkos {
 public:
  GridVecMonodKokkos(class LAMMPS *);
  ~GridVecMonodKokkos() {}
  void init();
  void grow(int);

  int pack_comm(int, int *, double *);
  void unpack_comm(int, int *, double *);
  int pack_exchange(int, int *, double *);
  void unpack_exchange(int, int *, double *);

  int pack_comm_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &);
  void unpack_comm_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &);
  int pack_exchange_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &);
  void unpack_exchange_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &);

  void set(int, char **);
  void set_monod(int, double);
  void set_monod(int, double, double, double, double, double, double, double);

  void sync(ExecutionSpace, unsigned int);
  void modified(ExecutionSpace, unsigned int);
  void sync_overlapping_device(ExecutionSpace, unsigned int);

 private:
  int *mask;
  double **conc;    // concentration
  double **reac;    // reaction rate
  double **dens;    // density
  double ***growth; // growth rate

  DAT::t_int_1d d_mask;
  HAT::t_int_1d h_mask;
  DAT::t_float_2d d_conc;
  HAT::t_float_2d h_conc;
  DAT::t_float_2d d_reac;
  HAT::t_float_2d h_reac;
  DAT::t_float_2d d_dens;
  HAT::t_float_2d h_dens;
  DAT::t_float_3d d_growth;
  HAT::t_float_3d h_growth;
};

}

#endif
#endif
