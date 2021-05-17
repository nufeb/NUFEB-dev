/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_coccus.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecCoccus::AtomVecCoccus(LAMMPS *lmp) : AtomVec(lmp)
{
  mass_type = PER_ATOM;
  molecular = Atom::ATOMIC;

  atom->coccus_flag = 1;
  atom->sphere_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag =
    atom->torque_flag =  atom->outer_radius_flag  = atom->outer_mass_flag =
      atom->biomass_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "radius rmass omega torque outer_mass outer_radius biomass";
  fields_copy = (char *) "radius rmass omega outer_mass outer_radius biomass";
  fields_comm = (char *) "";
  fields_comm_vel = (char *) "omega";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "radius rmass outer_mass outer_radius biomass";
  fields_border_vel = (char *) "radius rmass omega outer_mass outer_radius biomass";
  fields_exchange = (char *) "radius rmass omega outer_mass outer_radius biomass";
  fields_restart = (char *) "radius rmass omega outer_mass outer_radius biomass";
  fields_create = (char *) "radius rmass omega outer_mass outer_radius biomass";
  fields_data_atom = (char *) "id type radius rmass x outer_radius biomass";
  fields_data_vel = (char *) "id v omega";
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecCoccus::process_args(int narg, char **arg)
{
  if (narg != 0 && narg != 1)
    error->all(FLERR,"Illegal atom_style coccus command");

  radvary = 1;
  if (narg == 1) {
    radvary = utils::numeric(FLERR,arg[0],true,lmp);
    if (radvary < 0 || radvary > 1)
      error->all(FLERR,"Illegal atom_style coccus command");
  }

  // dynamic particle properties must be communicated every step

  if (radvary) {
    fields_comm = (char *) "radius rmass outer_mass outer_radius";
    fields_comm_vel = (char *) "radius rmass omega outer_mass outer_radius";
  }

  // delay setting up of fields until now

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecCoccus::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag && radvary == 0)
        error->all(FLERR,"Fix adapt changes particle radii "
                   "but atom_style coccus is not dynamic");
    } else if (strcmp(modify->fix[i]->style,"nufeb/monod") == 0) {
	if (radvary == 0)
	  error->all(FLERR,"Fix nufeb/monod changes particle radii "
		     "but atom_style coccus is not dynamic");
    } else if (strcmp(modify->fix[i]->style,"nufeb/divide") == 0) {
	if (radvary == 0)
	  error->all(FLERR,"Fix nufeb/divide changes particle radii "
		     "but atom_style coccus is not dynamic");
    }
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecCoccus::grow_pointers()
{
  radius = atom->radius;
  rmass = atom->rmass;
  omega = atom->omega;
  biomass = atom->biomass;
  outer_radius = atom->outer_radius;
  outer_mass = atom->outer_mass;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecCoccus::create_atom_post(int ilocal)
{
  radius[ilocal] = 0.5e-6;
  rmass[ilocal] = 4.0*MY_PI/3.0 * 0.5e-6*0.5e-6*0.5e-6;
  biomass[ilocal] = rmass[ilocal];
  outer_radius[ilocal] = radius[ilocal];
  outer_mass[ilocal] = 0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecCoccus::data_atom_post(int ilocal)
{
  radius_one = 0.5 * atom->radius[ilocal];
  radius[ilocal] = radius_one;
  if (radius_one > 0.0)
    rmass[ilocal] *= 4.0*MY_PI/3.0 * radius_one*radius_one*radius_one;

  if (rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;

  outer_radius_one = 0.5 * atom->outer_radius[ilocal];
  outer_radius[ilocal] = outer_radius_one;
  if (outer_radius[ilocal] < radius[ilocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }

  outer_mass[ilocal] = (4.0*MY_PI/3.0)*
    ((outer_radius[ilocal]*outer_radius[ilocal]*outer_radius[ilocal])
     -(radius[ilocal]*radius[ilocal]*radius[ilocal])) * 30;

  biomass_one = atom->biomass[ilocal];
  if (biomass_one < 0 || biomass_one > 1)
    error->one(FLERR,"Biomass/Mass (dry/wet weight) ratio must be between 0-1");
  biomass[ilocal] = biomass_one*rmass[ilocal];
}


/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecCoccus::pack_data_pre(int ilocal)
{
  radius_one = radius[ilocal];
  rmass_one = rmass[ilocal];

  radius[ilocal] *= 2.0;
  if (radius_one != 0.0)
    rmass[ilocal] =
      rmass_one / (4.0*MY_PI/3.0 * radius_one*radius_one*radius_one);

  outer_radius_one = outer_radius[ilocal];
  outer_radius[ilocal] *= 2.0;

  biomass_one = biomass[ilocal];
  biomass[ilocal] = biomass_one/rmass_one;
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecCoccus::pack_data_post(int ilocal)
{
  radius[ilocal] = radius_one;
  rmass[ilocal] = rmass_one;
  outer_radius[ilocal] = outer_radius_one;
  biomass[ilocal] = biomass_one;
}
