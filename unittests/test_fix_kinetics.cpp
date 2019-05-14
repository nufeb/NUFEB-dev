/*
 * test_fix_kinetics.cpp
 *
 *  Created on: 22 Oct 2018
 *      Author: Bowen Li
 */

#include "test_fix_kinetics.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixKineticsTest : public ::testing::Test {
  protected:
    FixKineticsTest() {
      kinetics = NULL;
      avec = NULL;

      char *argv[] = {"fix_kinetics", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      // Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixKineticsTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 100 5.0e-7");

      lmp->input->one("read_data_bio inputs/template.in");
      // Create group
      lmp->input->one("group AOB type 1");
      // Define variables
      lmp->input->one("variable diffT equal 1e-4");
      lmp->input->one("variable layer equal -1");
      lmp->input->one("variable EPSdens equal 30");

      avec = (AtomVecBio *) lmp->atom->style_match("bio");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    AtomVecBio *avec;
    FixKinetics *kinetics;
  };

  /*
   *  Check form activities after 1 step
   */
 TEST_F(FixKineticsTest, Compute_Activity) {
   lmp->input->one("timestep 3600");
   lmp->input->one("fix k1 all kinetics 1 2 2 2 v_diffT v_layer");
   lmp->input->one("fix ki2 all kinetics/thermo");
   lmp->input->one("fix ki3 all kinetics/ph fix");
   lmp->input->one("fix ki1 all kinetics/growth/energy v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   // run lammps
   lmp->input->one("run 1");

   EXPECT_NEAR(9.864638e-06, kinetics->activity[1][2][1], 1e-12);
   EXPECT_NEAR(1.672985e-24, kinetics->activity[2][1][1], 1e-26);
   EXPECT_NEAR(2.81e-04, kinetics->activity[3][1][1], 1e-6);
   EXPECT_NEAR(1.368558e-03, kinetics->activity[4][2][1], 1e-9);
 }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
