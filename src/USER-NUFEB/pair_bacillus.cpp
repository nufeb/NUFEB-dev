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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_bacillus.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_bacillus.h"
#include "comm.h"
#include "force.h"
#include "fix.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathConst;

#define EPSILON 1e-30

enum {SPHERE,ROD};

/* ---------------------------------------------------------------------- */

PairBacillus::PairBacillus(LAMMPS *lmp) : Pair(lmp)
{
  nmax = 0;

  c_n = 0.1;
  c_t = 0.2;
  mu = 0.0;
  cutoff = 0.0;

  maxrad = nullptr;

  k_n = nullptr;
  k_na = nullptr;
}

/* ---------------------------------------------------------------------- */

PairBacillus::~PairBacillus()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
  }
}

/* ---------------------------------------------------------------------- */

void PairBacillus::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double rsq, leni, lenj, radi, radj;
  int ishape, jshape;
  double xtmp,ytmp,ztmp,delx,dely,delz, evdwl, facc[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  int *bacillus = atom->bacillus;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  AtomVecBacillus::Bonus *ibonus;
  AtomVecBacillus::Bonus *jbonus;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (bacillus[i] >= 0) {
      int index = atom->bacillus[i];
      ibonus = &avec->bonus[index];
      leni = ibonus->length/2;
      radi = ibonus->diameter/2;
      leni == 0 ? ishape = SPHERE:ishape = ROD;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (bacillus[i] < 0 || bacillus[j] < 0) continue;

      int index = atom->bacillus[j];
      jbonus = &avec->bonus[index];
      lenj = jbonus->length/2;
      radj = jbonus->diameter/2;
      lenj == 0 ? jshape = SPHERE:jshape = ROD;

      // no interaction
      double r = sqrt(rsq);
      if (r > radi+radj+leni+lenj+cutoff) continue;

      // sphere-sphere interaction

      if (ishape == SPHERE && jshape == SPHERE) {
        sphere_against_sphere(i, j, itype, jtype, delx, dely, delz,
                              rsq, v, f, radi, radj, evflag);
        continue;
      }

      // one of the two bacillus is a sphere
      if (jshape == SPHERE) {
        sphere_against_rod(i, j, itype, jtype, x, v, f, torque,
                            angmom, ibonus, evflag);
        continue;
      } else if (ishape == SPHERE) {
        sphere_against_rod(j, i, jtype, itype, x, v, f, torque,
                            angmom, jbonus, evflag);
        continue;
      }

      int contact = 0;
      Contact contact_list;

      // rod-rod interaction
      rod_against_rod(i, j, itype, jtype, x, v, f, torque, angmom, ibonus,
		      jbonus, contact, contact_list, evdwl, facc);

      if (contact > 0) {
        rescale_cohesive_forces(x, f, torque, contact_list, contact,
                                itype, jtype, facc);
      }

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBacillus::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_n,n+1,n+1,"pair:k_n");
  memory->create(k_na,n+1,n+1,"pair:k_na");
  memory->create(maxrad,n+1,"pair:maxerad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBacillus::settings(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal pair_style bacillus command");

  c_n = utils::numeric(FLERR,arg[0],true,lmp);
  c_t = utils::numeric(FLERR,arg[1],true,lmp);
  mu = utils::numeric(FLERR,arg[2],true,lmp);
  cutoff = utils::numeric(FLERR,arg[3],true,lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBacillus::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double k_n_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k_na_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBacillus::init_style()
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Pair bacillus requires "
                        "atom style bacillus");

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style bacillus requires "
               "newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair bacillus requires "
               "ghost atoms store velocity");

  neighbor->request(this);

  int i, itype;
  double radi, leni;
  int* bacillus = atom->bacillus;
  int* type = atom->type;
  int ntypes = atom->ntypes;
  int nlocal = atom->nlocal;

  double *merad = NULL;
  memory->create(merad,ntypes+1,"pair:merad");
  for (i = 1; i <= ntypes; i++)
    maxrad[i] = merad[i] = 0;

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  int idiv;
  for (idiv = 0; idiv < modify->nfix; idiv++)
    if (strcmp(modify->fix[idiv]->style,"nufeb/divide/bacillus/minicell") == 0) break;
  if (idiv == modify->nfix) idiv = -1;

  for (i = 1; i <= ntypes; i++) {
    merad[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[idep]->extract("radius",itype));
    }
    if (idiv >= 0) {
      itype = i;
      merad[i] =
	*((double *) modify->fix[idiv]->extract("radius",itype));
    }
  }

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if (bacillus[i] >= 0) {
      int index = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[index];
      leni = bonus->length/2;
      radi = bonus->diameter/2;
      if (leni+radi > merad[itype])
	merad[itype] = leni+radi;
    } else
      merad[itype] = 0;
  }

  MPI_Allreduce(&merad[1],&maxrad[1],ntypes,MPI_DOUBLE,MPI_MAX,world);

  memory->destroy(merad);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBacillus::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];

  return (maxrad[i]+maxrad[j]);
}


/* ----------------------------------------------------------------------
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBacillus::sphere_against_sphere(int i, int j,
  int itype, int jtype, double delx, double dely, double delz, double rsq,
  double** v, double** f, double rradi, double rradj, int evflag)
{
  double contact_dist;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair, energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  contact_dist = rradi + rradj;

  rij = sqrt(rsq);
  R = rij - contact_dist;

  energy = 0;
  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= 0) { // in contact
    // relative translational velocity

    vr1 = v[i][0] - v[j][0];
    vr2 = v[i][1] - v[j][1];
    vr3 = v[i][2] - v[j][2];

    // normal component

    rsqinv = 1.0/rsq;
    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // normal friction term at contact

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    // tangential friction term at contact,
    // excluding the tangential deformation term for now

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];
  }

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Interaction bt a rod (irod) and a sphere (jsphere)
---------------------------------------------------------------------- */

void PairBacillus::sphere_against_rod(int i, int j,
  int itype, int jtype, double** x, double** v, double** f, double** torque,
  double** angmom, AtomVecBacillus::Bonus *&ibonus, int evflag)
{
  int ni,nei,ifirst,iefirst,npi1,npi2;
  double xi1[3],xi2[3],vti[3],h[3],fn[3],ft[3],d,t;
  double delx,dely,delz,rsq,rij,rsqinv,R,fx,fy,fz,fpair,energy;
  double radi,radj,leni,contact_dist;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double *quat, *inertia;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  radi = atom->radius[i];
  radj = atom->radius[j];
  leni = ibonus->length/2;
  contact_dist = radi + radj;

  avec->get_pole_coords(i, xi1, xi2);

  // find shortest distance between i and j
  distance_bt_pt_rod(x[j], xi1, xi2, h, d, t);

  if (d > contact_dist + cutoff) return;
  if (t < 0 || t > 1) return;

  delx = h[0] - x[j][0];
  dely = h[1] - x[j][1];
  delz = h[2] - x[j][2];
  rsq = delx*delx + dely*dely + delz*delz;
  rij = sqrt(rsq);
  R = d - contact_dist;

  energy = 0;
  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  if (R <= 0) { // in contact

    // compute the velocity of the vertex in the space-fixed frame
    quat = ibonus->quat;
    inertia = ibonus->inertia;
    total_velocity(h, x[i], v[i], angmom[i],
		   inertia, quat, vti);

    // relative translational velocity

    vr1 = vti[0] - v[j][0];
    vr2 = vti[1] - v[j][1];
    vr3 = vti[2] - v[j][2];

    // normal component

    rsqinv = 1.0/rsq;
    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // normal friction term at contact

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    // tangential friction term at contact,
    // excluding the tangential deformation term

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];
  }

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;
  sum_torque(x[i], h, fx, fy, fz, torque[i]);

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                         energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's edges against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   contact_list = list of contacts
   num_contacts = number of contacts between i's edges and j's edges
   Return:

---------------------------------------------------------------------- */

void PairBacillus::rod_against_rod(int i, int j, int itype, int jtype,
				   double** x, double** v, double** f, double** torque,
				   double** angmom, AtomVecBacillus::Bonus *&ibonus,
				   AtomVecBacillus::Bonus *&jbonus, int &contact,
				   Contact &contact_list, double &evdwl, double* facc)
{
  double xi1[3],xi2[3],xpj1[3],xpj2[3];
  double r,t1,t2,h1[3],h2[3];
  double contact_dist, energy;

  avec->get_pole_coords(i, xi1, xi2);
  avec->get_pole_coords(j, xpj1, xpj2);

  contact_dist = (ibonus->diameter + jbonus->diameter)/2;

  int jflag = 1;
//  printf("  line1:\n p1 = (%e %e %e);\n p2 = (%e %e %e)\n \n"
//         "  line2:\n p1 = (%e %e %e);\n p2 = (%e %e %e)\n: "
//         "t1 = %f; t2 = %f; r = %e\n",
//    xi1[0], xi1[1], xi1[2], xi2[0], xi2[1], xi2[2],
//    xpj1[0], xpj1[1], xpj1[2], xpj2[0], xpj2[1], xpj2[2],
//    t1, t2, r);
  distance_bt_rods(xpj1, xpj2, xi1, xi2, h2, h1, t2, t1, r);

  // include the vertices for interactions
  if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1 &&
      r < contact_dist + cutoff) {

    pair_force_and_torque(j, i, h2, h1, r, contact_dist,
                          jtype, itype, x, v, f, torque,
			  angmom, jflag, energy, facc);

    if (r <= contact_dist) {
      // store the contact info
      contact_list.i = i;
      contact_list.j = j;
      contact_list.xi[0] = h1[0];
      contact_list.xi[1] = h1[1];
      contact_list.xi[2] = h1[2];
      contact_list.xj[0] = h2[0];
      contact_list.xj[1] = h2[1];
      contact_list.xj[2] = h2[2];
      contact_list.separation = r - contact_dist;
      contact = 1;
    }
  }

  evdwl += energy;
}


/* ----------------------------------------------------------------------
  Compute forces and torques between two bodies caused by the interaction
  between a pair of points on either bodies (similar to sphere-sphere)
------------------------------------------------------------------------- */

void PairBacillus::pair_force_and_torque(int i, int j,
                 double* pi, double* pj, double r, double contact_dist,
                 int itype, int jtype, double** x, double** v, double** f,
		 double** torque, double** angmom, int jflag, double& energy,
		 double* facc)
{
  double delx,dely,delz,R,fx,fy,fz,fpair;

  delx = pi[0] - pj[0];
  dely = pi[1] - pj[1];
  delz = pi[2] - pj[2];
  R = r - contact_dist;

  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/r;
  fy = dely*fpair/r;
  fz = delz*fpair/r;

  if (R <= 0) {

    // contact: accumulate normal and tangential contact force components

    contact_forces(i, j, pi, pj, delx, dely, delz, fx, fy, fz,
                   x, v, angmom, f, torque, facc);
  } else {

    // accumulate force and torque to both bodies directly

    f[i][0] += fx;
    f[i][1] += fy;
    f[i][2] += fz;
    sum_torque(x[i], pi, fx, fy, fz, torque[i]);

    facc[0] += fx; facc[1] += fy; facc[2] += fz;

    if (jflag) {
      f[j][0] -= fx;
      f[j][1] -= fy;
      f[j][2] -= fz;
      sum_torque(x[j], pj, -fx, -fy, -fz, torque[j]);
    }
  }
}

/* ----------------------------------------------------------------------
  Compute contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
  fx,fy,fz = unscaled cohesive forces
  fn = normal friction component
  ft = tangential friction component (-c_t * v_t)
------------------------------------------------------------------------- */

void PairBacillus::contact_forces(int i, int j, double *xi, double *xj,
				  double delx, double dely, double delz,
				  double fx, double fy, double fz, double** x,
				  double** v, double** angmom, double** f,
				  double** torque, double* facc)
{
  int ibonus,jbonus;
  double rsq,rsqinv,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3],vj[3];
  double *quat, *inertia;
  AtomVecBacillus::Bonus *bonus;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->bacillus[i];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xi, x[i], v[i], angmom[i],
                 inertia, quat, vi);

  // compute the velocity of the point on the edge in the space-fixed frame

  jbonus = atom->bacillus[j];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xj, x[j], v[j], angmom[j],
                 inertia, quat, vj);

  // vector pointing from the contact point on ibody to that on jbody

  rsq = delx*delx + dely*dely + delz*delz;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = vi[0] - vj[0];
  vr2 = vi[1] - vj[1];
  vr3 = vi[2] - vj[2];

  // normal component

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // normal friction term at contact

  fn[0] = -c_n * vn1;
  fn[1] = -c_n * vn2;
  fn[2] = -c_n * vn3;

  // tangential friction term at contact
  // excluding the tangential deformation term for now

  ft[0] = -c_t * vt1;
  ft[1] = -c_t * vt2;
  ft[2] = -c_t * vt3;

  // these are contact forces (F_n, F_t and F_ne) only
  // cohesive forces will be scaled by j_a after contact area is computed
  // mu * fne = tangential friction deformation during gross sliding
  // see Eq. 4, Fraige et al.

  fx = fn[0] + ft[0] + mu * fx;
  fy = fn[1] + ft[1] + mu * fy;
  fz = fn[2] + ft[2] + mu * fz;

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;
  sum_torque(x[i], xi, fx, fy, fz, torque[i]);

  f[j][0] -= fx;
  f[j][1] -= fy;
  f[j][2] -= fz;
  sum_torque(x[j], xj, -fx, -fy, -fz, torque[j]);

  facc[0] += fx; facc[1] += fy; facc[2] += fz;
}

/* ----------------------------------------------------------------------
 compute the shortest distance between sphere (point) and rod (line segments)
------------------------------------------------------------------------- */

void PairBacillus::distance_bt_pt_rod(const double* q,
     const double* xi1, const double* xi2, double* h, double& d, double& t)
{
  double vx = xi2[0] - xi1[0];
  double vy = xi2[1] - xi1[1];
  double vz = xi2[2] - xi1[2];

  double wx = q[0] - xi1[0];
  double wy = q[1] - xi1[1];
  double wz = q[2] - xi1[2];

  double c1 = dot(wx, wy, wz, vx, vy, vz);
  double c2 = dot(vx, vy, vz, vx, vy, vz);

  if (c1 <= 0) {
    t = 0;
    h[0] = xi1[0];
    h[1] = xi1[1];
    h[2] = xi1[2];

  } else if (c2 <= c1) {
    t = 1;
    h[0] = xi2[0];
    h[1] = xi2[1];
    h[2] = xi2[2];

  } else {
    t = c1 / c2;
    h[0] = xi1[0] + t * vx;
    h[1] = xi1[1] + t * vy;
    h[2] = xi1[2] + t * vz;
  }

  d = dist(q[0], q[1], q[2], h[0], h[1], h[2]);
}

/* ----------------------------------------------------------------------
 compute the shortest distance between two rods (line segments)
------------------------------------------------------------------------- */

void PairBacillus::distance_bt_rods(const double* x1,
		  const double* x2, const double* x3, const double* x4,
		  double* h1, double* h2, double& t1, double& t2, double& r)
{
    double ux = x2[0] - x1[0];
    double uy = x2[1] - x1[1];
    double uz = x2[2] - x1[2];

    double vx = x4[0] - x3[0];
    double vy = x4[1] - x3[1];
    double vz = x4[2] - x3[2];

    double wx = x1[0] - x3[0];
    double wy = x1[1] - x3[1];
    double wz = x1[2] - x3[2];

    double a = dot(ux, uy, uz, ux, uy, uz);         // always >= 0
    double b = dot(ux, uy, uz, vx, vy, vz);
    double c = dot(vx, vy, vz, vx, vy, vz);        // always >= 0
    double d = dot(ux, uy, uz, wx, wy, wz);
    double e = dot(vx, vy, vz, wx, wy, wz);
    double dd = a * c - b * b;
    double nt1, dt1 = dd;
    double nt2, dt2 = dd;

    // compute the line parameters of the two closest points
    if (dd < EPSILON) { // the lines are almost parallel
      nt1 = 0.0;
      dt1 = 1.0;         // to prevent possible division by 0.0 later
      nt2 = e;
      dt2 = c;
    } else {
      nt1 = (b*e - c*d);
      nt2 = (a*e - b*d);
      if (nt1 < 0.0) {
	nt1 = 0.0;
	nt2 = e;
	dt2 = c;
      } else if (nt1 > dt1) {
	nt1 = dt1;
	nt2 = e + b;
	dt2 = c;
      }
    }

    if (nt2 < 0.0) {
      nt2 = 0.0;

      if (-d < 0.0)
	nt1 = 0.0;
      else if (-d > a)
	nt1 = dt1;
      else {
	nt1 = -d;
	dt1 = a;
      }
    }
    else if (nt2 > dt2) {
      nt2 = dt2;
      // recompute sc for this edge
      if ((-d + b) < 0.0)
	nt1 = 0;
      else if ((-d + b) > a)
	nt1 = dt1;
      else {
	nt1 = (-d +  b);
	dt1 = a;
      }
    }

    // finally do the division to get t2 and t1
    t1 = (fabs(nt1) < EPSILON ? 0.0 : nt1 / dt1);
    t2 = (fabs(nt2) < EPSILON ? 0.0 : nt2 / dt2);

    double r1 = wx + t1*ux - t2*vx;
    double r2 = wy + t1*uy - t2*vy;
    double r3 = wz + t1*uz - t2*vz;

    r = norm(r1, r2, r3);   // return the closest distance

    h1[0] = x1[0] + (x2[0] - x1[0]) * t1;
    h1[1] = x1[1] + (x2[1] - x1[1]) * t1;
    h1[2] = x1[2] + (x2[2] - x1[2]) * t1;

    h2[0] = x3[0] + (x4[0] - x3[0]) * t2;
    h2[1] = x3[1] + (x4[1] - x3[1]) * t2;
    h2[2] = x3[2] + (x4[2] - x3[2]) * t2;
}

/* ----------------------------------------------------------------------
  Kernel force is model-dependent and can be derived for other styles
    here is the harmonic potential (linear piece-wise forces) in Wang et al.
------------------------------------------------------------------------- */

void PairBacillus::kernel_force(double R, int itype, int jtype,
  double& energy, double& fpair)
{
  double kn = k_n[itype][jtype];
  double kna = k_na[itype][jtype];
  double shift = kna * cutoff;
  double e = 0;
  if (R <= 0) {           // deformation occurs
    fpair = -kn * R - shift;
    e = (0.5 * kn * R + shift) * R;
  } else if (R <= cutoff) {   // not deforming but cohesive ranges overlap
    fpair = kna * R - shift;
    e = (-0.5 * kna * R + shift) * R;
  } else fpair = 0.0;
  energy += e;
}

/* ----------------------------------------------------------------------
  Rescale the forces and torques for all the contacts
------------------------------------------------------------------------- */

void PairBacillus::rescale_cohesive_forces(double** x,  double** f,
					   double** torque, Contact& contact_list,
					   int &num_contacts, int itype, int jtype,
					   double* facc)
{
  int m,ibody,jbody;
  double delx,dely,delz,fx,fy,fz,R,fpair,r;

  ibody = contact_list.i;
  jbody = contact_list.j;

  delx = contact_list.xi[0] - contact_list.xj[0];
  dely = contact_list.xi[1] - contact_list.xj[1];
  delz = contact_list.xi[2] - contact_list.xj[2];
  r = sqrt(delx*delx + dely*dely + delz*delz);
  R = contact_list.separation;

  double energy = 0;
  kernel_force(R, itype, jtype, energy, fpair);

  fx = delx*fpair/r;
  fy = dely*fpair/r;
  fz = delz*fpair/r;

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], contact_list.xi, fx, fy, fz, torque[ibody]);

  f[jbody][0] -= fx;
  f[jbody][1] -= fy;
  f[jbody][2] -= fz;
  sum_torque(x[jbody], contact_list.xj, -fx, -fy, -fz, torque[jbody]);

  facc[0] += fx; facc[1] += fy; facc[2] += fz;
}

/* ----------------------------------------------------------------------
  Accumulate torque to body from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

void PairBacillus::sum_torque(double* xm, double *x, double fx,
                                      double fy, double fz, double* torque)
{
  double rx = x[0] - xm[0];
  double ry = x[1] - xm[1];
  double rz = x[2] - xm[2];
  double tx = ry * fz - rz * fy;
  double ty = rz * fx - rx * fz;
  double tz = rx * fy - ry * fx;
  torque[0] += tx;
  torque[1] += ty;
  torque[2] += tz;
}


/* ----------------------------------------------------------------------
  Calculate the total velocity of a point (vertex, a point on an edge):
    vi = vcm + omega ^ (p - xcm)
------------------------------------------------------------------------- */

void PairBacillus::total_velocity(double* p, double *xcm,
  double* vcm, double *angmom, double *inertia, double *quat, double* vi)
{
  double r[3],omega[3],ex_space[3],ey_space[3],ez_space[3];
  r[0] = p[0] - xcm[0];
  r[1] = p[1] - xcm[1];
  r[2] = p[2] - xcm[2];
  MathExtra::q_to_exyz(quat,ex_space,ey_space,ez_space);
  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,
                             inertia,omega);
  vi[0] = omega[1]*r[2] - omega[2]*r[1] + vcm[0];
  vi[1] = omega[2]*r[0] - omega[0]*r[2] + vcm[1];
  vi[2] = omega[0]*r[1] - omega[1]*r[0] + vcm[2];
}
