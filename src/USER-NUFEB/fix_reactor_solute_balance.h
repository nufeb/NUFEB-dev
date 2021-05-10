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

FixStyle(nufeb/reactor/solute_balance,FixReactorSoluteBalance)

#else

#ifndef LMP_FIX_REACTOR_SOLUTE_BALANCE_H
#define LMP_FIX_REACTOR_SOLUTE_BALANCE_H

#include "fix_reactor.h"

namespace LAMMPS_NS {

class FixReactorSoluteBalance : public FixReactor {
 public:
  FixReactorSoluteBalance(class LAMMPS *, int, char **);
  virtual ~FixReactorSoluteBalance() {}
  double compute_scalar();
  virtual void init();
  virtual void compute();

 protected:
  int iliq;
  double inlet;

  double q;
  double rvol;
  double reactor_af;
  double domain_af;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
