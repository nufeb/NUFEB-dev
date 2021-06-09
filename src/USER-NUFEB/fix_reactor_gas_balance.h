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

FixStyle(nufeb/reactor/gas_balance,FixReactorGasBalance)

#else

#ifndef LMP_FIX_REACTOR_GAS_BALANCE_H
#define LMP_FIX_REACTOR_GAS_BALANCE_H

#include "fix_reactor.h"

namespace LAMMPS_NS {

class FixReactorGasBalance : public FixReactor {
 public:

  FixReactorGasBalance(class LAMMPS *, int, char **);
  ~FixReactorGasBalance();
  double compute_scalar();
  virtual void init();
  virtual void compute();

 protected:
  int igas;

  double reactor_vhead;
  double reactor_pres;
  int nfix_gas_liquid;

  class FixGasLiquid **fix_gas_liquid;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
