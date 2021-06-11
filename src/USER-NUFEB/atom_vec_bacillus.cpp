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

#include "atom_vec_bacillus.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "memory.h"
#include "modify.h"
#include "math_const.h"
#include "math_extra.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1.0e-7

enum{SPHERE,ROD};       // also in DumpImage
/* ---------------------------------------------------------------------- */

AtomVecBacillus::AtomVecBacillus(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  bonus_flag = 1;

  // first 3 sizes do not include values from body itself
  // 1st,2nd body counts are added in process_args() via body style
  // 3rd body count is added in size_restart_bonus()
  // size_data_bonus is not used by Atom class for body style

  size_forward_bonus = 4;
  size_border_bonus = 16;
  size_restart_bonus_one = 16;
  size_data_bonus = 11;

  atom->bacillus_flag = 1;
  atom->rmass_flag = 1;
  atom->biomass_flag = 1;
  atom->angmom_flag = atom->torque_flag = 1;
  atom->radius_flag = 1;

  nlocal_bonus = nghost_bonus = nmax_bonus = 0;
  bonus = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) "radius rmass angmom torque bacillus biomass";
  fields_copy = (char *) "radius rmass angmom biomass";
  fields_comm = (char *) "rmass";
  fields_comm_vel = (char *) "angmom rmass";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "radius rmass biomass";
  fields_border_vel = (char *) "radius rmass angmom biomass";
  fields_exchange = (char *) "radius rmass angmom biomass";
  fields_restart = (char *) "radius rmass angmom biomass";
  fields_create = (char *) "radius rmass angmom bacillus biomass";
  fields_data_atom = (char *) "id type bacillus rmass x";
  fields_data_vel = (char *) "id v angmom";
}

/* ---------------------------------------------------------------------- */

AtomVecBacillus::~AtomVecBacillus()
{
  memory->sfree(bonus);
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecBacillus::process_args(int narg, char **arg)
{
  if (narg != 0 && narg != 1)
    error->all(FLERR,"Illegal atom_style bacillus command");

  radvary = 0;
  if (narg == 1) {
    radvary = utils::numeric(FLERR,arg[0],true,lmp);
    if (radvary < 0 || radvary > 1)
      error->all(FLERR,"Illegal atom_style bacillus command");
  }

  // dynamic particle properties must be communicated every step

  if (radvary) {
    fields_comm = (char *) "radius rmass";
    fields_comm_vel = (char *) "radius rmass angmom";
  }

  // delay setting up of fields until now

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag && radvary == 0)
        error->all(FLERR,"Fix adapt changes particle radii "
                   "but atom_style bacillus is not dynamic");
    }
}
/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecBacillus::grow_pointers()
{
  radius = atom->radius;
  bacillus = atom->bacillus;
  rmass = atom->rmass;
  biomass = atom->biomass;
  angmom = atom->angmom;
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecBacillus::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0)
    error->one(FLERR,"Per-processor system is too big");

  bonus = (Bonus *) memory->srealloc(bonus,nmax_bonus*sizeof(Bonus),
                                     "atom:bonus");
}

/* ----------------------------------------------------------------------
   copy atom I bonus info to atom J
------------------------------------------------------------------------- */

void AtomVecBacillus::copy_bonus(int i, int j, int delflag)
{
  // if deleting atom J via delflag and J has bonus data, then delete it

  if (delflag && bacillus[j] >= 0) {
    copy_bonus_all(nlocal_bonus-1,bacillus[j]);
    nlocal_bonus--;
  }

  // if atom I has bonus data, reset I's bonus.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's bonus data is already deleted

  if (bacillus[i] >= 0 && i != j) bonus[bacillus[i]].ilocal = j;
  bacillus[j] = bacillus[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset ellipsoid that points to I to now point to J
------------------------------------------------------------------------- */

void AtomVecBacillus::copy_bonus_all(int i, int j)
{
  bacillus[bonus[i].ilocal] = j;
  memcpy(&bonus[j],&bonus[i],sizeof(Bonus));
}

/* ----------------------------------------------------------------------
   clear ghost info in bonus data
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void AtomVecBacillus::clear_bonus()
{
  nghost_bonus = 0;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_comm_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (bacillus[j] >= 0) {
      quat = bonus[bacillus[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBacillus::unpack_comm_bonus(int n, int first, double *buf)
{
  int i,m,last;
  double *quat;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (bacillus[i] >= 0) {
      quat = bonus[bacillus[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::pack_border_bonus(int n, int *list, double *buf)
{
  int i,j,m;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (bacillus[j] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      quat = bonus[bacillus[j]].quat;
      inertia = bonus[bacillus[j]].inertia;
      pole1 = bonus[bacillus[j]].pole1;
      pole2 = bonus[bacillus[j]].pole2;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      buf[m++] = inertia[0];
      buf[m++] = inertia[1];
      buf[m++] = inertia[2];
      buf[m++] = pole1[0];
      buf[m++] = pole1[1];
      buf[m++] = pole1[2];
      buf[m++] = pole2[0];
      buf[m++] = pole2[1];
      buf[m++] = pole2[2];
      buf[m++] = bonus[bacillus[j]].length;
      buf[m++] = bonus[bacillus[j]].diameter;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_border_bonus(int n, int first, double *buf)
{
  int i,j,m,last;
  double *quat,*inertia,*pole1,*pole2;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    bacillus[i] = (int) ubuf(buf[m++]).i;
    if (bacillus[i] == 0) bacillus[i] = -1;
    else {
      j = nlocal_bonus + nghost_bonus;
      if (j == nmax_bonus) grow_bonus();

      quat = bonus[j].quat;
      inertia = bonus[j].inertia;
      pole1 = bonus[j].pole1;
      pole2 = bonus[j].pole2;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      inertia[0] = buf[m++];
      inertia[1] = buf[m++];
      inertia[2] = buf[m++];
      pole1[0] = buf[m++];
      pole1[1] = buf[m++];
      pole1[2] = buf[m++];
      pole2[0] = buf[m++];
      pole2[1] = buf[m++];
      pole2[2] = buf[m++];
      bonus[j].length = buf[m++];
      bonus[j].diameter = buf[m++];
      bonus[j].ilocal = i;
      bacillus[i] = j;
      nghost_bonus++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_exchange_bonus(int i, double *buf)
{
  int m = 0;

  if (bacillus[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = bacillus[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    double *pole1= bonus[j].pole1;
    double *pole2= bonus[j].pole2;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = pole1[0];
    buf[m++] = pole1[1];
    buf[m++] = pole1[2];
    buf[m++] = pole2[0];
    buf[m++] = pole2[1];
    buf[m++] = pole2[2];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].diameter;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBacillus::unpack_exchange_bonus(int ilocal, double *buf)
{
  int m = 0;

  bacillus[ilocal] = (int) ubuf(buf[m++]).i;
  if (bacillus[ilocal] == 0) bacillus[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    double *pole1= bonus[nlocal_bonus].pole1;
    double *pole2= bonus[nlocal_bonus].pole2;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    pole1[0] = buf[m++];
    pole1[1] = buf[m++];
    pole1[2] = buf[m++];
    pole2[0] = buf[m++];
    pole2[1] = buf[m++];
    pole2[2] = buf[m++];
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].diameter = buf[m++];

    bonus[nlocal_bonus].ilocal = ilocal;
    bacillus[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecBacillus::size_restart_bonus()
{
  int i;

  int n = 0;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) {
    if (bacillus[i] >= 0) n += size_restart_bonus_one;
    else n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including bonus data
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_restart_bonus(int i, double *buf)
{
  int m = 0;

  if (bacillus[i] < 0) buf[m++] = ubuf(0).d;
  else {
    buf[m++] = ubuf(1).d;
    int j = bacillus[i];
    double *quat = bonus[j].quat;
    double *inertia = bonus[j].inertia;
    double *pole1= bonus[j].pole1;
    double *pole2= bonus[j].pole2;
    buf[m++] = quat[0];
    buf[m++] = quat[1];
    buf[m++] = quat[2];
    buf[m++] = quat[3];
    buf[m++] = inertia[0];
    buf[m++] = inertia[1];
    buf[m++] = inertia[2];
    buf[m++] = pole1[0];
    buf[m++] = pole1[1];
    buf[m++] = pole1[2];
    buf[m++] = pole2[0];
    buf[m++] = pole2[1];
    buf[m++] = pole2[2];
    buf[m++] = bonus[j].length;
    buf[m++] = bonus[j].diameter;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including bonus data
------------------------------------------------------------------------- */

int AtomVecBacillus::unpack_restart_bonus(int ilocal, double *buf)
{
  int m = 0;

  bacillus[ilocal] = (int) ubuf(buf[m++]).i;
  if (bacillus[ilocal] == 0) bacillus[ilocal] = -1;
  else {
    if (nlocal_bonus == nmax_bonus) grow_bonus();
    double *quat = bonus[nlocal_bonus].quat;
    double *inertia = bonus[nlocal_bonus].inertia;
    double *pole1= bonus[nlocal_bonus].pole1;
    double *pole2= bonus[nlocal_bonus].pole2;
    quat[0] = buf[m++];
    quat[1] = buf[m++];
    quat[2] = buf[m++];
    quat[3] = buf[m++];
    inertia[0] = buf[m++];
    inertia[1] = buf[m++];
    inertia[2] = buf[m++];
    pole1[0] = buf[m++];
    pole1[1] = buf[m++];
    pole1[2] = buf[m++];
    pole2[0] = buf[m++];
    pole2[1] = buf[m++];
    pole2[2] = buf[m++];
    bonus[nlocal_bonus].length = buf[m++];
    bonus[nlocal_bonus].diameter = buf[m++];
    bonus[nlocal_bonus].ilocal = ilocal;
    bacillus[ilocal] = nlocal_bonus++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack one line from Bacilli section of data file
------------------------------------------------------------------------- */

void AtomVecBacillus::data_atom_bonus(int m, char **values)
{
  if (bacillus[m])
    error->one(FLERR,"Assigning bacillus parameters to non-ellipsoid atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  // diagonalize inertia tensor

  double tensor[3][3];
  tensor[0][0] = utils::numeric(FLERR,values[0],true,lmp);
  tensor[1][1] = utils::numeric(FLERR,values[1],true,lmp);
  tensor[2][2] = utils::numeric(FLERR,values[2],true,lmp);
  tensor[0][1] = tensor[1][0] = utils::numeric(FLERR,values[3],true,lmp);
  tensor[0][2] = tensor[2][0] = utils::numeric(FLERR,values[4],true,lmp);
  tensor[1][2] = tensor[2][1] = utils::numeric(FLERR,values[5],true,lmp);

  double *inertia = bonus[nlocal_bonus].inertia;
  double evectors[3][3];
  int ierror = MathExtra::jacobi(tensor,inertia,evectors);
  if (ierror) error->one(FLERR,
                         "Insufficient Jacobi rotations for bacillus");

  // if any principal moment < scaled EPSILON, set to 0.0
  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON * max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON * max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON * max) inertia[2] = 0.0;

  // exyz_space = principal axes in space frame
  double ex_space[3],ey_space[3],ez_space[3];

  ex_space[0] = evectors[0][0];
  ex_space[1] = evectors[1][0];
  ex_space[2] = evectors[2][0];
  ey_space[0] = evectors[0][1];
  ey_space[1] = evectors[1][1];
  ey_space[2] = evectors[2][1];
  ez_space[0] = evectors[0][2];
  ez_space[1] = evectors[1][2];
  ez_space[2] = evectors[2][2];

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed
  double cross[3];
  MathExtra::cross3(ex_space,ey_space,cross);
  if (MathExtra::dot3(cross,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion
  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,bonus[nlocal_bonus].quat);

  // initialise cooridnates os two cell poles
  double *pole1 = bonus[nlocal_bonus].pole1;
  double *pole2 = bonus[nlocal_bonus].pole2;
  double px = utils::numeric(FLERR,values[6],true,lmp);
  double py = utils::numeric(FLERR,values[7],true,lmp);
  double pz = utils::numeric(FLERR,values[8],true,lmp);

  pole1[0] = px;
  pole1[1] = py;
  pole1[2] = pz;

  double d = sqrt(px*px + py*py + pz*pz);

  bonus[nlocal_bonus].length = d * 2;

  pole2[0] = -px;
  pole2[1] = -py;
  pole2[2] = -pz;

  bonus[nlocal_bonus].diameter = utils::numeric(FLERR,values[9],true,lmp);
  if (bonus[nlocal_bonus].diameter < 0)
    error->one(FLERR, "Invalid diameter in Bacilli section of data file: diameter < 0");
  radius[m] = bonus[nlocal_bonus].diameter * 0.5;

  // reset bacillus mass
  // previously stored density in rmass
  if (radius[m] > 0.0)
  rmass[m] *= (4.0*MY_PI/3.0*radius[m]*radius[m]*radius[m] +
      MY_PI*radius[m]*radius[m]*bonus[nlocal_bonus].length);

  biomass[m] = 1.0;

  bonus[nlocal_bonus].ilocal = m;
  bacillus[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated bonus memory
------------------------------------------------------------------------- */

double AtomVecBacillus::memory_usage_bonus()
{
  double bytes = 0;
  bytes += nmax_bonus*sizeof(Bonus);
  return bytes;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBacillus::create_atom_post(int ilocal)
{
  radius[ilocal] = 0.5e-6;
  rmass[ilocal] *= (4.0*MY_PI/3.0*
      radius[ilocal]*radius[ilocal]*radius[ilocal] +
      MY_PI*radius[ilocal]*radius[ilocal]*bonus[nlocal_bonus].length);
  biomass[ilocal] = 1.0;
  bacillus[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBacillus::data_atom_post(int ilocal)
{
  bacillus_flag = bacillus[ilocal];
  if (bacillus_flag == 0) bacillus_flag = -1;
  else if (bacillus_flag == 1) bacillus_flag = 0;
  else error->one(FLERR,"Invalid bacillus flag in Atoms section of data file");
  bacillus[ilocal] = bacillus_flag;

  if (rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  angmom[ilocal][0] = 0.0;
  angmom[ilocal][1] = 0.0;
  angmom[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecBacillus::pack_data_pre(int ilocal)
{
  bacillus_flag = atom->bacillus[ilocal];
  rmass_one = atom->rmass[ilocal];

  if (bacillus_flag < 0) bacillus[ilocal] = 0;
  else bacillus[ilocal] = 1;

  if (bacillus_flag >= 0)
    rmass[ilocal] /= (4.0*MY_PI/3.0*radius[ilocal]*radius[ilocal]*radius[ilocal] +
	      MY_PI*radius[ilocal]*radius[ilocal]*bonus[bacillus[ilocal]].length);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecBacillus::pack_data_post(int ilocal)
{
  bacillus[ilocal] = bacillus_flag;
  rmass[ilocal] = rmass_one;
  biomass[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   pack bonus bacillus info for writing to data file
   if buf is nullptr, just return buffer size
------------------------------------------------------------------------- */

int AtomVecBacillus::pack_data_bonus(double *buf, int /*flag*/)
{
  int i,j,m;
  double p[3][3],pdiag[3][3],ispace[3][3];
  double *pole1 = bonus->pole1;
  double *inertia = bonus->inertia;
  double *quat = bonus->quat;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (bacillus[i] < 0) continue;
    if (buf) {
      buf[m++] = ubuf(tag[i]).d;
      j = bacillus[i];
      // 6 moments of inertia

      MathExtra::quat_to_mat(quat,p);
      MathExtra::times3_diag(p,inertia,pdiag);
      MathExtra::times3_transpose(pdiag,p,ispace);

      buf[m++] = ispace[0][0];
      buf[m++] = ispace[1][1];
      buf[m++] = ispace[2][2];
      buf[m++] = ispace[0][1];
      buf[m++] = ispace[0][2];
      buf[m++] = ispace[1][2];

      buf[m++] = pole1[0];
      buf[m++] = pole1[1];
      buf[m++] = pole1[2];

      buf[m++] = bonus->length;
    } else m += size_data_bonus;
  }

  return m;
}

/* ----------------------------------------------------------------------
   write bonus bacillus info to data file
------------------------------------------------------------------------- */

void AtomVecBacillus::write_data_bonus(FILE *fp, int n, double *buf, int /*flag*/)
{
  int i = 0;
  while (i < n) {
    fmt::print(fp,"{} {} {} {} {} {} {} {} {} {} {}\n",ubuf(buf[i]).i,
	       buf[i+1],buf[i+2],buf[i+3],buf[i+4],buf[i+5],buf[i+6],buf[i+7],
	       buf[i+8],buf[i+9],buf[i+10]);
    i += size_data_bonus;
  }
}

/* ----------------------------------------------------------------------
   set values in bonus data for bacillus m
------------------------------------------------------------------------- */

void AtomVecBacillus::set_bonus(int m, double *pole1, double diameter, double *quat, double *inertia)
{
  if (bacillus[m])
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  if (nlocal_bonus == nmax_bonus) grow_bonus();

  bonus[nlocal_bonus].pole1[0] = pole1[0];
  bonus[nlocal_bonus].pole1[1] = pole1[1];
  bonus[nlocal_bonus].pole1[2] = pole1[2];
  bonus[nlocal_bonus].pole2[0] = -pole1[0];
  bonus[nlocal_bonus].pole2[1] = -pole1[1];
  bonus[nlocal_bonus].pole2[2] = -pole1[2];

  bonus[nlocal_bonus].inertia[0] = inertia[0];
  bonus[nlocal_bonus].inertia[1] = inertia[1];
  bonus[nlocal_bonus].inertia[2] = inertia[2];

  bonus[nlocal_bonus].quat[0] = quat[0];
  bonus[nlocal_bonus].quat[1] = quat[1];
  bonus[nlocal_bonus].quat[2] = quat[2];
  bonus[nlocal_bonus].quat[3] = quat[3];

  double d = sqrt(pole1[0]*pole1[0] + pole1[1]*pole1[1] + pole1[2]*pole1[2]);

  bonus[nlocal_bonus].length = d * 2;
  bonus[nlocal_bonus].diameter = diameter;

  bonus[nlocal_bonus].ilocal = m;
  bacillus[m] = nlocal_bonus++;
}

/* ----------------------------------------------------------------------
   get pole coordinate for bacillus m
------------------------------------------------------------------------- */
void AtomVecBacillus::get_pole_coords(int m, double *xp1, double *xp2)
{
  if (bacillus[m] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  double p[3][3];
  double *x;

  MathExtra::quat_to_mat(bonus[bacillus[m]].quat,p);
  MathExtra::matvec(p,bonus[bacillus[m]].pole1, xp1);
  MathExtra::matvec(p,bonus[bacillus[m]].pole2, xp2);

  x = atom->x[bonus[bacillus[m]].ilocal];

  xp1[0] += x[0];
  xp1[1] += x[1];
  xp1[2] += x[2];
  xp2[0] += x[0];
  xp2[1] += x[1];
  xp2[2] += x[2];
}
