/*
 *
 *  Created on: 22 Oct 2018
 *      Author: Bowen Li
 */

#include "test_fix_division.h"
#include <string>
#include <sstream>

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixDivisionTest : public ::testing::Test {
  protected:
    FixDivisionTest() {
      avec = NULL;

      char *argv[] = {"fix_division", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixDivisionTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");

      lmp->input->one("read_data_bio inputs/division.in");

      // Create group
      lmp->input->one("group HET type 1");
      // Define variables
      lmp->input->one("variable EPSdens equal 30");
      // Define variables
      lmp->input->one("variable divDia equal 1.5631e-6");

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
   *  Check for daughter cells after division
   */
 TEST_F(FixDivisionTest, Daughter_Cell) {
   // Make sure initial mass is greater than divMass
   ASSERT_GT(lmp->atom->rmass[0], 2.77e-16);

   EXPECT_NEAR(1.256037e-13, avec->outer_mass[0], 1e-18);

   char* prefix = "fix d1 all divide 1 v_EPSdens v_divDia ";

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

     ASSERT_EQ(lmp->atom->natoms,2) << "At " << i << " th iteration";

     // Mass range between 1.2e-16 - 1.81e-16
     double mass1 = lmp->atom->rmass[0];
     ASSERT_TRUE((mass1 >= 1.2e-16) && (mass1 <= 1.8e-16))
     << "At " << i << " th iteration; mass = " << mass1;
     double mass2 = lmp->atom->rmass[1];
     ASSERT_TRUE((mass2 >= 1.2e-16) && (mass2 <= 1.8e-16))
     << "At " << i << " th iteration mass = " << mass2;

     // Outer mass range between 0.502e-13 - 0.754e-13
     double omass1 = avec->outer_mass[0];
     ASSERT_TRUE((omass1 >= 0.502e-13) && (omass1 <= 0.754e-13))
     << "At " << i << " th iteration; outer mass = " << omass1;
     double omass2 = avec->outer_mass[1];
     ASSERT_TRUE((omass2 >= 0.502e-13) && (omass2 <= 0.754e-13))
     << "At " << i << " th iteration outer mass = " << omass2;

     // Type
     ASSERT_TRUE(lmp->atom->type[0] == 1) << "At " << i << " th iteration";
     ASSERT_TRUE(lmp->atom->type[1] == 1) << "At " << i << " th iteration";

     // Growth Rate
     ASSERT_NEAR(6.9444e-5, avec->bio->q[1], 1e-6) << "At " << i << " th iteration";

   }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
