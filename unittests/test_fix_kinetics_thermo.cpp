/*
 * test_fix_kinetics_thermo.cpp
 *
 *  Created on: 22 Oct 2018
 *      Author: Bowen Li
 */

#include "test_fix_kinetics_thermo.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixKineticsThermoTest : public ::testing::Test {
  protected:
    FixKineticsThermoTest() {
      kinetics = NULL;

      char *argv[] = {"fix_kinetics", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      // Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixKineticsThermoTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 100 5.0e-7");

      lmp->input->one("read_data_bio inputs/thermo.in");
      // Create group
      lmp->input->one("group AOB type 1");
      // Define variables
      lmp->input->one("variable diffT equal 1e-4");
      lmp->input->one("variable layer equal -1");
      lmp->input->one("variable EPSdens equal 30");
      // set timestep and fix
      lmp->input->one("timestep 3600");
      lmp->input->one("fix k1 all kinetics 1 2 2 2 v_diffT v_layer");
      lmp->input->one("fix ki3 all kinetics/ph fix ph 7.5");

      int ifix = lmp->modify->find_fix("k1");
      if (ifix >= 0) {
        kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
      }
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    FixKinetics *kinetics;
  };

  /*
   *  Check metabolic energy
   */
 TEST_F(FixKineticsThermoTest, Var_Metabolic_Energy) {

   lmp->input->one("fix ki2 all kinetics/thermo");
   lmp->input->one("fix ki1 all kinetics/growth/energy v_EPSdens");

   // run lammps
   lmp->input->one("run 5");

   // loop over all grids which should have the same value
   for (int i = 0; i < 8; i++ ) {
     EXPECT_NEAR(-411.0335, kinetics->gibbs_cata[1][i], 1e-4);
     EXPECT_NEAR(234.6530, kinetics->gibbs_anab[1][i], 1e-4);
   }
 }

 /*
  *  Check gas liquid
  */
  TEST_F(FixKineticsThermoTest, Var_Gas_Liquid) {
    // TODO
  }

  /*
   *  Check dynamic yield
   */
   TEST_F(FixKineticsThermoTest, Var_Yield) {
     // TODO
   }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
