/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include "fix_bio_kinetics.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"

#include "bio.h"
#include "atom_vec_bio.h"
#include "fix_bio_kinetics_ph.h"
#include "fix_bio_kinetics_thermo.h"
#include "fix_bio_kinetics_monod.h"
#include "fix_bio_immigration.h"
#include "fix_bio_kinetics_diffusion2.h"

#include "pointers.h"
#include "variable.h"
#include "modify.h"
#include "update.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  bio = avec->bio;

  if (narg != 13) error->all(FLERR,"Not enough arguments in fix kinetics command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics command: calling steps should be positive integer");

  var = new char*[6];
  ivar = new int[6];

  nx = atoi(arg[4]);
  ny = atoi(arg[5]);
  nz = atoi(arg[6]);

  for (int i = 0; i < 6; i++) {
    int n = strlen(&arg[7+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[7+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKinetics::~FixKinetics()
{
  int i;
  for (i = 0; i < 6; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  memory->destroy(gYield);
  memory->destroy(activity);
  memory->destroy(nuR);
  memory->destroy(nuS);
  memory->destroy(qGas);
  memory->destroy(DRGCat);
  memory->destroy(DRGAn);
  memory->destroy(kEq);
  memory->destroy(Sh);
  delete[] nuConv;;
}

/* ---------------------------------------------------------------------- */

int FixKinetics::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init()
{
  for (int n = 0; n < 6; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics is invalid style");
  }

  if (bio->nnus == 0)
    error->all(FLERR,"fix_kinetics requires # of Nutrients inputs");
  else if (bio->nuGCoeff == NULL)
    error->all(FLERR,"fix_kinetics requires Nutrient Energy inputs");
  else if (bio->iniS == NULL)
    error->all(FLERR,"fix_kinetics requires Nutrients inputs");

  // register fix kinetics with this class
  diffusion = NULL;
  monod = NULL;
  ph = NULL;
  thermo = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics/monod") == 0) {
      monod = static_cast<FixKineticsMonod *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"kinetics/diffusion2") == 0) {
      diffusion = static_cast<FixKineticsDiffusion2 *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"kinetics/ph") == 0) {
      ph = static_cast<FixKineticsPH *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"kinetics/thermo") == 0) {
      thermo = static_cast<FixKineticsThermo *>(lmp->modify->fix[j]);
    }
  }

  temp = input->variable->compute_equal(ivar[0]);
  rth = input->variable->compute_equal(ivar[1]);
  gVol = input->variable->compute_equal(ivar[2]);
  gasTrans = input->variable->compute_equal(ivar[3]);
  iph = input->variable->compute_equal(ivar[4]);
  diffT = input->variable->compute_equal(ivar[5]);

  ngrids = nx * ny * nz;

  nnus = bio->nnus;
  int ntypes = atom->ntypes;

  nuS = memory->create(nuS,nnus+1, ngrids, "kinetics:nuS");
  nuR = memory->create(nuR,nnus+1, ngrids, "kinetics:nuR");
  qGas = memory->create(qGas,nnus+1, ngrids, "kinetics:nuGas");
  gYield = memory->create(gYield,ntypes+1,ngrids,"kinetic:gYield");
  activity = memory->create(activity,nnus+1,5, ngrids,"kinetics:activity");
  DRGCat = memory->create(DRGCat,ntypes+1,ngrids,"kinetics:DRGCat");
  DRGAn = memory->create(DRGAn,ntypes+1,ngrids,"kinetics:DRGAn");
  kEq = memory->create(kEq,nnus+1,4,"kinetics:kEq");
  Sh = memory->create(Sh,ngrids,"kinetics:Sh");
  nuConv = new bool[nnus+1]();

  //initialize grid yield, inlet concentration, consumption
  for (int j = 0; j < ngrids; j++) {
    for (int i = 1; i <= ntypes; i++) {
      gYield[i][j] = bio->yield[i];
    }
    for (int i = 1; i <= nnus; i++) {
      nuS[i][j] = bio->iniS[i][0];
      nuR[i][j] = 0;
      qGas[i][j] = 0;
    }
  }

  init_keq();
  init_activity();
}

/* ---------------------------------------------------------------------- */

void FixKinetics::grow() {
  int ntypes = atom->ntypes;

  gYield = memory->grow(gYield,ntypes+1,ngrids,"kinetic:gYield");
  DRGCat = memory->grow(DRGCat,ntypes+1,ngrids,"kinetics:DRGCat");
  DRGAn = memory->grow(DRGAn,ntypes+1,ngrids,"kinetics:DRGAn");

  for (int j = 0; j < ngrids; j++) {
    for (int i = 1; i <= ntypes; i++) {
      gYield[i][j] = bio->yield[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init_keq()
{
  // water Kj/mol
  double dG0H2O = -237.18;
  int nnus = bio->nnus;
  double **nuGCoeff = bio->nuGCoeff;

  for (int i = 1; i < nnus+1; i++) {
    kEq[i][0] = exp((dG0H2O + nuGCoeff[i][0] - nuGCoeff[i][1]) / (-rth * temp));
    for (int j = 1; j < 4; j++) {
      double coeff = 0.0;

      if (nuGCoeff[i][j+1] > 10000) {
        coeff = j * 10001;
      } else {
        coeff = 0;
      }

      kEq[i][j] = exp((nuGCoeff[i][j+1] + coeff - nuGCoeff[i][j]) / (-rth * temp));
    }
  }
}

void FixKinetics::init_activity() {
  int nnus = bio->nnus;
  double *denm = memory->create(denm,nnus+1,"kinetics:denm");
  double gSh = pow(10, -iph);

  for (int k = 1; k < nnus+1; k++) {
    for (int j = 0; j < ngrids; j++) {
      double iniNuS = bio->iniS[k][0];
      Sh[j] = gSh;
      denm[k] = (1 + kEq[k][0]) * gSh * gSh * gSh + kEq[k][1] * gSh * gSh + kEq[k][2] * kEq[k][3] * gSh + kEq[k][3] * kEq[k][2] * kEq[k][1];
      if (denm[k] == 0) {
        lmp->error->all(FLERR,"denm returns a zero value");
      }
      // not hydrated form acitivity
      activity[k][0][j] = kEq[k][0] * nuS[k][j] * gSh * gSh * gSh / denm[k];
      // fully protonated form activity
      activity[k][1][j] = nuS[k][j] * gSh * gSh * gSh / denm[k];
      // 1st deprotonated form activity
      activity[k][2][j] = nuS[k][j] * gSh * gSh * kEq[k][1] / denm[k];
      // 2nd deprotonated form activity
      activity[k][3][j] = nuS[k][j] * gSh * kEq[k][1] * kEq[k][2] / denm[k];
      // 3rd deprotonated form activity
      activity[k][4][j] = nuS[k][j] * kEq[k][1] * kEq[k][2] * kEq[k][3] / denm[k];

      if (strcmp(bio->nuName[k], "h") == 0) {
        activity[k][1][j] = gSh;
      }
    }
  }

  memory->destroy(denm);
}

/* ---------------------------------------------------------------------- */

void FixKinetics::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  integration();
}

void FixKinetics::integration() {

  int iteration = 0;
  bool isConv = false;

  //initialization
  for (int i = 1; i <= nnus; i++) {
    if (strcmp(bio->nuName[i],"h2o")==0) nuConv[i] = true;
    else if (strcmp(bio->nuName[i],"h")==0) nuConv[i] = true;
    else nuConv[i] = false;
  }

  while (!isConv) {
    iteration++;
    isConv = true;

    if (ph != NULL) ph->solve_ph();
    else init_activity();
    if (thermo != NULL) thermo->thermo();
    if (monod != NULL) monod->growth(diffT);
    if (diffusion != NULL) nuConv = diffusion->diffusion(nuConv, iteration, diffT);
    else break;

    for (int i = 1; i <= nnus; i++){
      if (!nuConv[i]) {
        isConv = false;
        break;
      }
    }

    if (iteration > 10000) {
      isConv = true;
      for (int i = 1; i <= nnus; i++) {
        if (!nuConv[i]){
          nuConv[i] = true;
          printf( "%s  ", bio->nuName[i]);
        }
      }
    }
//    if (update->ntimestep > 20000 && iteration > 5000) {
//      error->all(FLERR,"# of iteration jumps to 5000, diffusion solver is unstable.");
//    }
  }

 printf( "number of iteration: %i \n", iteration);

 if (monod != NULL) monod->growth(update->dt * nevery);
}
