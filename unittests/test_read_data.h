#ifndef READ_DATA_TEST_H
#define READ_DATA_TEST_H

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "lammps.h"
#include "atom_vec.h"
#include "atom_vec_bio.h"
#include "bio.h"
#include "comm.h"
#include "input.h"
#include "mpi.h"

class ReadDataTest{
protected:
  ReadDataTest();
  virtual ~ReadDataTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
