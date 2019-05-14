#ifndef FIX_EPS_EXTRACT_TEST_H
#define FIX_EPS_EXTRACT_TEST_H

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "lammps.h"
#include "atom_vec.h"
#include "comm.h"
#include "input.h"
#include "mpi.h"
#include "atom_vec_bio.h"
#include <string>
#include <sstream>

class FixEPSExtractTest{
protected:
  FixEPSExtractTest();
  virtual ~FixEPSExtractTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
