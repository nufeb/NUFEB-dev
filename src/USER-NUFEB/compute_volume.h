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

#ifdef COMPUTE_CLASS

ComputeStyle(nufeb/volume,ComputeVolume)

#else

#ifndef LMP_COMPUTE_VOLUME_H
#define LMP_COMPUTE_VOLUME_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVolume : public Compute {
 public:
  ComputeVolume(class LAMMPS *, int, char **);
  virtual ~ComputeVolume() {}
  virtual void init() {}
  virtual double compute_scalar();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
