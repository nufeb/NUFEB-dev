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

#ifdef INTEGRATE_CLASS

IntegrateStyle(nufeb/kk,NufebRunKokkos)
IntegrateStyle(nufeb/kk/device,NufebRunKokkos)
IntegrateStyle(nufeb/kk/host,NufebRunKokkos)

#else

#ifndef LMP_NUFEB_RUN_KOKKOS_H
#define LMP_NUFEB_RUN_KOKKOS_H

#include "nufeb_run.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class NufebRunKokkos : public NufebRun {
 public:
  NufebRunKokkos(class LAMMPS *, int, char **);
  ~NufebRunKokkos() {}
  void init();
  void setup(int);
  void setup_minimal(int);
  void run(int);

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    f(i,0) += f_merge_copy(i,0);
    f(i,1) += f_merge_copy(i,1);
    f(i,2) += f_merge_copy(i,2);
  }

 protected:
  DAT::t_f_array f_merge_copy,f;

  void force_clear();
  void growth();
  void reactor();
  int diffusion();
  void disable_sync(class Fix *fix);
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
