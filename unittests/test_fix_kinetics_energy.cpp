/*
 * test_fix_kinetics_energy.cpp
 *
 *  Created on: 22 Oct 2018
 *      Author: Bowen Li
 */

#include "test_fix_kinetics_energy.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixKineticsEnergyTest : public ::testing::Test {
  protected:
    FixKineticsEnergyTest() {
      kinetics = NULL;

      char *argv[] = {"fix_kinetics_energy", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      // Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixKineticsEnergyTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 100 5.0e-7");

      lmp->input->one("read_data_bio inputs/energy.in");
      // Create group
      lmp->input->one("group AOB type 1");
      // Define variables
      lmp->input->one("variable diffT equal 1e-4");
      lmp->input->one("variable layer equal -1");
      lmp->input->one("variable EPSdens equal 30");
      // set timestep and fix
      lmp->input->one("timestep 3600");
      lmp->input->one("fix k1 all kinetics 1 2 2 2 v_diffT v_layer");
      lmp->input->one("fix ki2 all kinetics/thermo");
      lmp->input->one("fix ki3 all kinetics/ph fix ph 7.5");
      lmp->input->one("fix ki1 all kinetics/growth/energy v_EPSdens");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    FixKinetics *kinetics;
  };

  /*
   *  Check growth and consumption rates
   */
   TEST_F(FixKineticsEnergyTest, Func_Growth) {
     FixKineticsEnergy *energy;

     int ifix = lmp->modify->find_fix("k1");
     if (ifix >= 0) {
       kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
     }

     ifix = lmp->modify->find_fix("ki1");
     if (ifix >= 0) {
       energy = static_cast<FixKineticsEnergy *>(lmp->modify->fix[ifix]);
     }

     // run lammps
     lmp->input->one("run 1");

     for (int i = 0; i < 8; i++ ) {
       // growth rate
       if (i==0 || i==7) EXPECT_NEAR(8.31679e-06, energy->growrate[1][i], 1e-11);
       else EXPECT_EQ(0, energy->growrate[1][i]);

       // reaction rate
       if (i==0)EXPECT_NEAR(-2.478869e-06, kinetics->nur[1][i], 1e-12);
       else if  (i==7)EXPECT_NEAR(-4.957737e-06, kinetics->nur[1][i], 1e-12);
       else EXPECT_EQ(0, kinetics->nur[1][i]);
     }

     lmp->input->one("run 9");

     // biomass
     EXPECT_NEAR(1.343139e-15, lmp->atom->rmass[0], 1e-21);
     EXPECT_NEAR(2.686278e-15, lmp->atom->rmass[1], 1e-21);
   }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
