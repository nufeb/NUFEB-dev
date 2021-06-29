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

#include "fix_property_plasmid.h"

#include <cstdlib>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include "atom_vec_bacillus_ecoli.h"
#include "update.h"
#include "math_const.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define FILAMENT_VEL 0.1e-6;
#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixPropertyPlasmid::FixPropertyPlasmid(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg),
  xpm(nullptr), vpm(nullptr), filament(nullptr), xold(nullptr), reac_t(nullptr)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"fix nufeb/property/plasmid requires "
      "atom style bacillus");

  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/plasmid command");

  compute_flag = 1;
  scalar_flag = 1;

  replicate = 0.0;
  pmax = 5;
  pinit = 0;
  dia = 5e-7;
  acc = 1e-8;
  dt = 1;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "replicate") == 0) {
      replicate = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "max") == 0) {
      pmax = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "init") == 0) {
      pinit = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dia") == 0) {
      dia = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "acc") == 0) {
      acc = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }  else if (strcmp(arg[iarg], "dt") == 0) {
      dt = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/property/plasmid command");
    }
  }

  if (pinit > pmax) error->all(FLERR,"Illegal fix nufeb/property/plasmid command: "
      "initial plasmid cannot be more than maximum plasmid number");

  size_peratom_cols = 0;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */
FixPropertyPlasmid:: ~FixPropertyPlasmid() {
  memory->destroy(xpm);
  memory->destroy(vpm);
  memory->destroy(filament);
  memory->destroy(xold);
  memory->destroy(reac_t);

  delete random;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPropertyPlasmid::grow_arrays(int nmax)
{
  FixProperty::grow_arrays(nmax);
  memory->grow(xpm,nmax,pmax*3,"fix_nufeb/property/plasmid:xpm");
  memory->grow(xold,nmax,3,"fix_nufeb/property/plasmid:xold");
  memory->grow(vpm,nmax,pmax*3,"fix_nufeb/property/plasmid:vpm");
  memory->grow(reac_t,nmax,pmax,"fix_nufeb/property/plasmid:reac_t");
  memory->grow(filament,nmax,pmax,"fix_nufeb/property/plasmid:filament");
}

/* ---------------------------------------------------------------------- */

void FixPropertyPlasmid::init()
{
  for (int i = 0; i < atom->nlocal; i++) {
    AtomVecBacillus::Bonus *bouns = &avec->bonus[i];
    vprop[i] = pinit;
    // initialise plasmid position
    for (int j = 0; j < (int)vprop[i]; j++) {
      double xlimit = (bouns->length + bouns->diameter - dia) / 2;
      double x = -xlimit + 2*random->uniform()*xlimit;

      int jx0 = j*3;
      int jx1 = j*3+1;
      int jx2 = j*3+2;

      xpm[i][jx0] = x;
      xpm[i][jx1] = 0.0;
      xpm[i][jx2] = 0.0;

      vpm[i][jx0] = 0.0;
      vpm[i][jx1] = 0.0;
      vpm[i][jx2] = 0.0;

      reac_t[i][j] = -log(random->uniform())/replicate;;
      filament[i][j] = 0;
    }
  }
  dtsq = 0.5*dt*dt;
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */

void FixPropertyPlasmid::compute()
{
  int t = 0;
  while (t < update->dt) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit)
	motility(i);
    }
    t+=dt;
  }

  for (int i = 0; i < atom->nlocal; i++){
    replication(i);
  }

  // record x
  for (int i = 0; i < atom->nlocal; i++) {
    xold[i][0] = atom->x[i][0];
    xold[i][1] = atom->x[i][1];
    xold[i][2] = atom->x[i][2];
  }
}

/* ----------------------------------------------------------------------
   update plasmid position in atom i
------------------------------------------------------------------------- */

void FixPropertyPlasmid::motility(int i) {
  for (int m = 0; m < (int)vprop[i]; m++) {
    double pos[3];
    double xlimit[3];
    double ax, ay, az;

    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;

    ax = (2/sqrt(3)) * (-acc + 2*random->uniform()*acc);
    ay = (2/sqrt(3)) * (-acc + 2*random->uniform()*acc);
    az = (2/sqrt(3)) * (-acc + 2*random->uniform()*acc);

    if (filament[i][m] > 1) {
      filament[i][m]--;
      ax = 0.0;
      ay = 0.0;
      az = 0.0;
    } else if (filament[i][m] == 1) {
      filament[i][m]--;
      vpm[i][m0] = 0.0;
      vpm[i][m1] = 0.0;
      vpm[i][m2] = 0.0;
    } else {
      for (int n = 0; n < (int)vprop[i]; n++) {
        if (m == n) continue;

        int n0 = n*3;
        int n1 = n*3+1;
        int n2 = n*3+2;

        double rsq = (xpm[i][n2]-xpm[i][m2])*(xpm[i][n2]-xpm[i][m2])
  	  + (xpm[i][n1]-xpm[i][m1])*(xpm[i][n1]-xpm[i][m1])
  	  + (xpm[i][n0]-xpm[i][m0])*(xpm[i][n0]-xpm[i][m0]);

        if (rsq < dia*dia) {
	  filament[i][m] = 9;
	  filament[i][n] = 9;

	  if (xpm[i][m0] > xpm[i][n0]) {
	    vpm[i][m0] = FILAMENT_VEL;
	    vpm[i][n0] = -FILAMENT_VEL;
	  } else {
	    vpm[i][n0] = -FILAMENT_VEL;
	    vpm[i][m0] = FILAMENT_VEL;
	  }
	  vpm[i][m1] = 0.0;
	  vpm[i][m2] = 0.0;
	  vpm[i][n1] = 0.0;
	  vpm[i][n2] = 0.0;

	  ax = 0.0;
	  ay = 0.0;
	  az = 0.0;
        }
      }
    }

    pos[0] = xpm[i][m0] + vpm[i][m0]*dt + ax*dtsq;
    pos[1] = xpm[i][m1] + vpm[i][m1]*dt + ay*dtsq;
    pos[2] = xpm[i][m2] + vpm[i][m2]*dt + az*dtsq;

    vpm[i][m0] += ax*dt;
    vpm[i][m1] += ay*dt;
    vpm[i][m2] += az*dt;

    get_xlimit(xlimit, i, m);
    if (fabs(pos[0]) > xlimit[0]) {
      pos[0] = xpm[i][m0];
      vpm[i][m0] = 0.0;
    }
    if (fabs(pos[1]) > xlimit[1]) {
      pos[1] = xpm[i][m1];
      vpm[i][m1] = 0.0;
    }
    if (fabs(pos[2]) > xlimit[2]) {
      pos[2] = xpm[i][m2];
      vpm[i][m2] = 0.0;
    }

    xpm[i][m0] = pos[0];
    xpm[i][m1] = pos[1];
    xpm[i][m2] = pos[2];
  }
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */
void FixPropertyPlasmid::replication(int i) {
  double current_t = update->ntimestep * (1 + update->dt);
  double next_t = update->ntimestep * (1 + update->dt);
  double *ixpm = xpm[i];
  double *ivpm = vpm[i];

  for (int m = 0; m < (int)vprop[i]; m++) {
    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;

    while (reac_t[i][m] < next_t && (int)vprop[i] < pmax) {
      double ilimit[3];
      int n = vprop[i];
      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      double theta = random->uniform() * 2 * MY_PI;
      double phi = random->uniform() * (MY_PI);

      filament[i][n] = 0;
      ivpm[n0] = ivpm[m0];
      ivpm[n1] = ivpm[m1];
      ivpm[n2] = ivpm[m2];

      ixpm[n0] = ixpm[m0] + (dia * cos(theta) * sin(phi) * DELTA);
      ixpm[n1] = ixpm[m1] + (dia * sin(theta) * sin(phi) * DELTA);
      ixpm[n2] = ixpm[m2] + (dia * cos(phi) * DELTA);

      get_xlimit(ilimit, i, n);
      if (fabs(ixpm[n0]) > ilimit[0])
	ixpm[n0] = ilimit[0];
      if (fabs(ixpm[n1]) > ilimit[1])
	ixpm[n1] = ilimit[1];
      if (fabs(ixpm[n2]) > ilimit[2])
	ixpm[n2] = ilimit[2];

      vprop[i]++;
      reac_t[i][m] -= log(random->uniform())/replicate;
      reac_t[i][n] = reac_t[i][m];
    }
  }
}

/* ----------------------------------------------------------------------
   get coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_plasmid_coords(int i, int j, double *xp, double *x)
{
  xp[0] = x[0] + xpm[i][j*3];
  xp[1] = x[1] + xpm[i][j*3+1];
  xp[2] = x[2] + xpm[i][j*3+2];
}

/* ----------------------------------------------------------------------
   set related coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::set_plasmid_xpm(int i, int j, double *xp, double *x)
{
  xpm[i][j*3] = xp[0] - x[0];
  xpm[i][j*3+1] = xp[1] - x[1];
  xpm[i][j*3+2] = xp[2] - x[2];
}


/* ----------------------------------------------------------------------
   get maximum/minimum coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_xlimit(double *xlimit, int i, int j)
{
  AtomVecBacillus::Bonus *bouns = &avec->bonus[i];
  double d = bouns->length/2;
  double r = bouns->diameter/2;

  xlimit[0] = d + r - dia/2;
  if (fabs(xpm[i][j*3]) < d) xlimit[1] = r - dia/2;
  else xlimit[1] = sqrt(r*r - (fabs(xpm[i][j*3])-d)*(fabs(xpm[i][j*3])-d)) - dia/2;
  xlimit[2] = xlimit[1];
}

/* ----------------------------------------------------------------------
   update plasmid copy numbers in two daughter cells i, j
   called from fix_divide
------------------------------------------------------------------------- */
void FixPropertyPlasmid::update_arrays(int i, int j)
{
  int ibac = atom->bacillus[i];
  int jbac = atom->bacillus[j];
  double *ixpm = xpm[i];
  double *jxpm = xpm[j];
  double *ivpm = vpm[i];
  double *jvpm = vpm[j];
  int *dlist;
  memory->create(dlist,pmax,"fix nufeb/property/plasmid:dlist");

  AtomVecBacillus::Bonus *ibouns = &avec->bonus[ibac];
  AtomVecBacillus::Bonus *jbouns = &avec->bonus[jbac];

  vprop[j] = 0;

  for (int m = 0; m < (int)vprop[i]; m++) {
    double ilimit[3], jlimit[3], xp[3];
    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;

    get_plasmid_coords(i, m, xp, xold[i]);

    // plasmid transmission
    if (xp[0] < (atom->x[j][0]+atom->radius[j]+jbouns->length/2) &&
	xp[0] > (atom->x[j][0]-atom->radius[j]-jbouns->length/2)){
      // move plasmid k to cell j
      int n = (int)vprop[j];

      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      filament[j][n] = filament[i][m];
      jvpm[n0] = ivpm[m0];
      jvpm[n1] = ivpm[m1];
      jvpm[n2] = ivpm[m2];
      reac_t[j][n] = reac_t[i][m];

      set_plasmid_xpm(j, n, xp, atom->x[j]);
      get_xlimit(jlimit, j, n);
      if (fabs(jxpm[n0]) > jlimit[0])
	jxpm[n0] = jlimit[0];
      if (fabs(jxpm[n1]) > jlimit[1])
	jxpm[n1] = jlimit[1];
      if (fabs(jxpm[n2]) > jlimit[2])
	jxpm[n2] = jlimit[2];

      vprop[j]++;
      dlist[m] = 1;
    } else {
      dlist[m] = 0;
      set_plasmid_xpm(i, m, xp, atom->x[i]);
      get_xlimit(ilimit, i, m);
      if (fabs(jxpm[m0]) > ilimit[0])
	jxpm[m0] = ilimit[0];
      if (fabs(jxpm[m1]) > ilimit[1])
	jxpm[m1] = ilimit[1];
      if (fabs(jxpm[m2]) > ilimit[2])
	jxpm[m2] = ilimit[2];
    }
  }

  int k = 0;
  while (k < (int)vprop[i]) {
    // remove plasmid k from i
    if (dlist[k]) {
      copy_plasmid(i,k,(int)vprop[i]-1,1);
      dlist[k] = dlist[(int)vprop[i]-1];
      vprop[i]--;
    } else k++;
  }

  memory->destroy(dlist);
}


/* ----------------------------------------------------------------------
   compute average plasmid copy number
------------------------------------------------------------------------- */
double FixPropertyPlasmid::compute_scalar() {
  double result = 0.0;
  int n = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit)) {
      result += vprop[i];
      n++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_INT, MPI_SUM, world);

  if (n > 0) result /= n;
  else result = 0;

  return result;
}

/* ----------------------------------------------------------------------
   copy values within plasmid arrays
------------------------------------------------------------------------- */

void FixPropertyPlasmid::copy_plasmid(int i, int to, int from, int /*delflag*/)
{
  xpm[i][to*3] = xpm[i][from*3];
  xpm[i][to*3+1] = xpm[i][from*3+1];
  xpm[i][to*3+2] = xpm[i][from*3+2];

  vpm[i][to*3] = vpm[i][from*3];
  vpm[i][to*3+1] = vpm[i][from*3+1];
  vpm[i][to*3+2] = vpm[i][from*3+2];

  filament[i][to] = filament[i][from];
  reac_t[i][to] = reac_t[i][from];
}

/* ----------------------------------------------------------------------
   copy values within plasmid arrays
------------------------------------------------------------------------- */

void FixPropertyPlasmid::copy_arrays(int i, int j, int /*delflag*/)
{
  FixProperty::copy_arrays(i,j,0);

  for (int m = 0; m < pmax*3; m++) {
    xpm[j][m] = xpm[i][m];
    vpm[j][m] = vpm[i][m];
  }

  for (int m = 0; m < pmax; m++) {
    filament[j][m] = filament[i][m];
    reac_t[j][m] = reac_t[i][m];
  }

  for (int m = 0; m < 3; m++) {
    xold[j][m] = xold[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyPlasmid::pack_exchange(int i, double *buf)
{
  int m = FixProperty::pack_exchange(i,buf);

  for (int n = 0; n < pmax*3; n++) {
    buf[m++] = xpm[i][n];
    buf[m++] = vpm[i][n];
  }

  for (int n = 0; n < pmax; n++) {
    buf[m++] = filament[i][n];
    buf[m++] = reac_t[i][n];
  }

  buf[m++] = xold[i][0];
  buf[m++] = xold[i][1];
  buf[m++] = xold[i][2];

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyPlasmid::unpack_exchange(int nlocal, double *buf)
{
  int m = FixProperty::unpack_exchange(nlocal,buf);

  for (int n = 0; n < pmax*3; n++) {
    xpm[nlocal][n] = buf[m++];
    vpm[nlocal][n] = buf[m++];
  }

  for (int n = 0; n < pmax; n++) {
    filament[nlocal][n] = buf[m++];
    reac_t[nlocal][n] = buf[m++];
  }

  xold[nlocal][0] = buf[m++];
  xold[nlocal][1] = buf[m++];
  xold[nlocal][2] = buf[m++];

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyPlasmid::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = pmax*3 + 1;
  buf[m++] = vprop[i];

  for (int n = 0; n < pmax*3; n++)
    buf[m++] = xpm[i][n];

  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyPlasmid::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  vprop[nlocal] = extra[nlocal][m++];
  for (int n = 0; n < pmax*3; n++)
    xpm[nlocal][n] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyPlasmid::maxsize_restart()
{
  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyPlasmid::size_restart(int /*nlocal*/)
{
  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyPlasmid::memory_usage()
{
  double bytes;
  bytes += atom->nmax*sizeof(double);
  bytes += atom->nmax*pmax*sizeof(double);

  return bytes;
}
