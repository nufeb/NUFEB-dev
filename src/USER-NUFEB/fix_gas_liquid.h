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

#ifdef FIX_CLASS

FixStyle(nufeb/gas_liquid,FixGasLiquid)

#else

#ifndef LMP_FIX_GAS_LIQUID_H
#define LMP_FIX_GAS_LIQUID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGasLiquid : public Fix {
 public:
  int compute_flag;
  int iliquid;
  int igas;

  FixGasLiquid(class LAMMPS *, int, char **);
  ~FixGasLiquid() {}
  int modify_param(int, char **);
  int setmask();
  void post_integrate();
  void compute();

 protected:
  double kga;
  double h;
  double temp;
  double mw;
  double rg;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
