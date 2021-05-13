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

#ifndef LMP_FIX_DIVIDE_H
#define LMP_FIX_DIVIDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDivide : public Fix {
 public:
  int compute_flag;
  
  FixDivide(class LAMMPS *, int, char **);
  virtual ~FixDivide() {};
  int modify_param(int, char **);
  int setmask();
  void post_integrate();
  void post_neighbor();
  virtual void compute() = 0;
};

}

#endif

/* ERROR/WARNING messages:
*/
