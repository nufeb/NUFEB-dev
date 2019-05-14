/*
 * test_read_data.cpp
 *
 *  Created on: 22 Oct 2018
 *      Author: Bowen Li
 */

#include "test_read_data.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class ReadDataTest : public ::testing::Test {
  protected:
    ReadDataTest() {
      char *argv[] = {"read_data_bio", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~ReadDataTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      avec = (AtomVecBio *) lmp->atom->style_match("bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");
      //Initialize simulation domain, type, atom and nutrients
      lmp->input->one("read_data_bio inputs/read_data.in");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    AtomVecBio *avec;
  };

 TEST_F(ReadDataTest, Coeff_Atom) {
   EXPECT_DOUBLE_EQ(5e-7, avec->outer_radius[0]);
   EXPECT_NEAR(1.37e-17, avec->outer_mass[0], 1e-18);
   EXPECT_DOUBLE_EQ(0.000069444, avec->bio->q[1]);
   EXPECT_NEAR(9.81e-18, lmp->atom->rmass[0], 1e-19);
 }

  TEST_F(ReadDataTest, Coeff_Type) {
    EXPECT_STREQ("het", avec->bio->tname[1]);
    EXPECT_DOUBLE_EQ(0.000069444, avec->bio->q[1]);
    EXPECT_DOUBLE_EQ(4e-3, avec->bio->ks[1][1]);
    EXPECT_DOUBLE_EQ(1e-2, avec->bio->ks[1][2]);
    EXPECT_DOUBLE_EQ(0.63, avec->bio->yield[1]);

    //Check Catabolism and Anabolim Coeffs
    EXPECT_DOUBLE_EQ(-1, avec->bio->cata_coeff[1][1]);
    EXPECT_DOUBLE_EQ(0, avec->bio->cata_coeff[1][2]);

    EXPECT_DOUBLE_EQ(0, avec->bio->anab_coeff[1][1]);
    EXPECT_DOUBLE_EQ(1.5, avec->bio->anab_coeff[1][2]);
  }

 TEST_F(ReadDataTest, Coeff_Nutrient) {
   EXPECT_TRUE(avec->bio->nnu == 2);

   EXPECT_STREQ("cod", avec->bio->nuname[1]);
   EXPECT_TRUE(avec->bio->nustate[1] == 0);
   EXPECT_DOUBLE_EQ(1e-9, avec->bio->diff_coeff[1]);

   EXPECT_STREQ("ac", avec->bio->nuname[2]);
   EXPECT_TRUE(avec->bio->nustate[2] == 1);
 }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
