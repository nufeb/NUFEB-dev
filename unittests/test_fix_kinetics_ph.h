#ifndef FIX_KINETICS_PH_TEST_H
#define FIX_KINETICS_PH_TEST_H

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "lammps.h"
#include "atom_vec.h"
#include "atom_vec_bio.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_thermo.h"
#include "fix_bio_kinetics_energy.h"
#include "fix_bio_kinetics_ph.h"
#include "bio.h"
#include "comm.h"
#include "input.h"
#include "mpi.h"
#include "modify.h"
#include "update.h"
#include "library.h"

class FixKineticsPHTest{
protected:
  FixKineticsPHTest();
  virtual ~FixKineticsPHTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
