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

#ifndef FIX_AVE_TIME_H
#define FIX_AVE_TIME_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class LAMMPS *, int, char **);
  ~FixAveTime();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);

 private:
  int me;
  int nrepeat,nfreq,nvalid,irepeat,which,ifix;
  char *id;
  FILE *fp;

  int sflag,vflag,ave,nwindow,nsum;
  double scalar,*vector;
  int ncompute;
  class Compute **compute;
  class Fix *fix;

  int norm,iwindow,window_limit;
  double scalar_total;
  double *scalar_list;
  double *vector_total;
  double **vector_list;
};

}

#endif
