/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_eps_extract.h"
#include <string>
#include <sstream>

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixEPSExtractTest : public ::testing::Test {
  protected:
    FixEPSExtractTest() {
      avec = NULL;

      char *argv[] = {"fix_division", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixEPSExtractTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");
      // Read input data
      lmp->input->one("read_data_bio inputs/eps_extract.in");
      // Create groups
      lmp->input->one("group HET type 1");
      lmp->input->one("group EPS type 2");
      // Define variables
      lmp->input->one("variable EPSdens equal 30");
      lmp->input->one("variable EPSratio equal 1.25");

      avec = (AtomVecBio *) lmp->atom->style_match("bio");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS *lmp;
    AtomVecBio *avec;
  };

  /*
   *  Check for mass, type etc after eps extraction
   */
 TEST_F(FixEPSExtractTest, EPS_Cell) {
   // Make sure the actual eps radio is greater than 1.25
   // Parent cell outer mass = 1.001e-15
   ASSERT_GT((avec->outer_radius[0]/lmp->atom->radius[0]), 1.25);

   char* prefix = "fix e1 HET eps_extract 1 v_EPSratio v_EPSdens ";

   // Test with 10 random cases
   for (int i = 0; i < 10; i++) {
     TearDown();
     SetUp();

     //initialize 10 random seeds
     int iseed = rand() % 1000 + 1;
     char str[4];
     char *seed;
     sprintf(str,"%d",iseed);
     seed = str;

     char arg[50];
     strcpy(arg, prefix);
     strcat(arg, str);

     lmp->input->one(arg);
     lmp->input->one("run 1");
     ASSERT_EQ(lmp->atom->natoms, 3) << "At " << i << " th iteration";
     // EPS type
     ASSERT_EQ(lmp->atom->type[2], 2);
     // EPS mass
     double mass = lmp->atom->rmass[2];
     ASSERT_TRUE((mass >= 0.39e-15) && (mass <= 0.61e-15))
     << "At " << i << " th iteration mass = " << mass;
     // Parent cell outer mass, note that id has changed from 0 to 1
     double omass = avec->outer_mass[1];
     ASSERT_TRUE((omass >= 0.39e-15) && (omass <= 0.61e-15))
     << "At " << i << " th iteration mass = " << omass;
   }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
