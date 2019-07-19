/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "dump_bio_sstl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <dirent.h>

#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "fix.h"
#include "modify.h"
#include "atom.h"
#include "bio.h"
#include "fix_bio_kinetics.h"
#include "math_const.h"
#include "pointers.h"

#include "compute_bio_height.h"
#include "compute_bio_ntypes.h"
#include "compute_bio_biomass.h"
#include "compute_bio_pressure.h"


struct stat sst = {0};

using namespace LAMMPS_NS;
using namespace MathConst;

// customize by adding keyword

/* ---------------------------------------------------------------------- */

DumpBioSSTL::DumpBioSSTL(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg <= 4) error->all(FLERR,"No dump sstl arguments specified");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal dump sstl command");

  nkeywords = narg - 4;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;
  keywords = NULL;
  fp = NULL;

  nus_flag = 0;
  mass_flag = 0;
  ntypes_flag = 0;
  height_flag = 0;
  pres_flag = 0;
  volfrac_flag = 0;

  mass_header = 0;
  type_header = 0;

  // customize for new sections
  keywords = (char **) memory->srealloc(keywords, nkeywords*sizeof(char *), "keywords");

  for (int i = 0; i < nkeywords; i++) {
    int n = strlen(arg[i+4]) + 2;
    keywords[i] = new char[n];
    strcpy(keywords[i], arg[i+4]);
  }
}

/* ---------------------------------------------------------------------- */

DumpBioSSTL::~DumpBioSSTL()
{
  for (int i = 0; i < nkeywords; i++) {
    delete [] keywords[i];
  }

  memory->sfree(keywords);
  if (volfrac_flag) memory->destroy(vol_frac);
}


/* ---------------------------------------------------------------------- */
void DumpBioSSTL::init_style()
{
  // register fix kinetics with this class
  kinetics = NULL;
  nfix = modify->nfix;
  ncompute = modify->ncompute;

  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required");

  bio = kinetics->bio;

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  int i = 0;
  int nnus = kinetics->bio->nnu;
  int ntypes = atom->ntypes;

  // create directory
  if (stat("./SSTL", &sst) == -1) {
      mkdir("./SSTL", 0700);
  }

  // create sstl graph model
  int len = 22;
  char path[len];
  strcpy(path, "./SSTL/sstl_model.txt");

  filename = path;
  fp = fopen(filename,"a");
  write_sstl_model();
  fclose(fp);

  while (i < nkeywords) {
    if (strcmp(keywords[i],"con") == 0) {
      nus_flag = 1;
    } else if (strcmp(keywords[i],"biomass") == 0) {
      mass_flag = 1;
    } else if (strcmp(keywords[i],"ave_height") == 0) {
      height_flag = 1;
    } else if (strcmp(keywords[i],"ntypes") == 0) {
      ntypes_flag = 1;
    } else if (strcmp(keywords[i],"pressure") == 0) {
      pres_flag = 1;
    } else if (strcmp(keywords[i],"vof") == 0) {
      volfrac_flag = 1;
      vol_frac = memory->create(vol_frac, atom->ntypes + 1, kinetics->ngrids, "dumpsstl:vol_frac");

      for (int j = 0; j < kinetics->ngrids; j++)
        for (int i = 0; i <= atom->ntypes; i++)
        	vol_frac[i][j] = 0;
    }

    else lmp->error->all(FLERR,"Undefined dump_sstl keyword");

    i++;
  }

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style,"ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      continue;
    }  else if (strcmp(modify->compute[j]->style,"ntypes") == 0) {
      ctype = static_cast<ComputeNufebNtypes *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"biomass") == 0) {
      cmass = static_cast<ComputeNufebBiomass *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"pressure") == 0) {
      cpres = static_cast<ComputeNufebPressure *>(lmp->modify->compute[j]);
      continue;
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write()
{
  if (update-> ntimestep == 0) return;

  if (ntypes_flag == 1) ctype->compute_vector();
  if (mass_flag == 1) cmass->compute_vector();
  if (height_flag == 1)  cheight->compute_scalar();
  if (pres_flag == 1) cpres->compute_scalar();

  int nnus = kinetics->bio->nnu;
  int ntypes = atom->ntypes;


  if (mass_flag == 1 && comm->me == 0) {
    int len = 35;
    char path[len];
    strcpy(path, "./SSTL/biomass.txt");

    filename = path;
    fp = fopen(filename,"a");
    write_biomass_data();
    fclose(fp);
  }

  if (height_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./SSTL/ave_height.txt");

    filename = path;
    fp = fopen(filename,"a");
    write_height_data();
    fclose(fp);
  }

  if (ntypes_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./SSTL/ntypes.txt");

    filename = path;
    fp = fopen(filename,"a");
    write_ntype_data();
    fclose(fp);
  }

  if (pres_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./SSTL/pressure.txt");

    filename = path;
    fp = fopen(filename,"a");
    write_pressure_data();
    fclose(fp);
  }

  if (nus_flag == 1) {
    for (int i = 1; i < nnus+1; i++) {
	    char *name = bio->nuname[i];
		int len = 30;
		len += strlen(name);
		char path[len];
		strcpy(path, "./SSTL/s_");
		strcat(path, name);
		strcat(path, ".txt");

		filename = path;
	    fp = fopen(filename,"a");
		write_nus_data(i);
		fclose(fp);
	  }
  }

  if (volfrac_flag == 1) {
    for (int i = 0; i < ntypes+1; i++) {
    	char *name;
	    if (i == 0) name = "all";
	    else name = bio->tname[i];
		int len = 30;
		len += strlen(name);
		char path[len];
		strcpy(path, "./SSTL/volf_");
		strcat(path, name);
		strcat(path, ".txt");

		filename = path;
	    fp = fopen(filename,"a");
	    update_volf();
	    write_volf_data(i);
		fclose(fp);
	  }
  }
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_header(bigint n)
{
}

void DumpBioSSTL::write_data(int n, double *mybuf)
{
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::pack(tagint *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_sstl_model()
{
  int nxyz = kinetics->nx*kinetics->ny*kinetics->nz;
  int ngrids = kinetics->ngrids;

  if (comm->me == 0) {
	  fprintf(fp, "LOCATIONS\n");
	  for(int i = 0; i < nxyz; i++){
		 fprintf(fp, "%i\n", i+1);
	  }

	  fprintf(fp, "EDGES\n");

	  // write edges
	  for(int i = 0; i < nxyz; i++){
		double gridx[3] = {0,0,0};
		get_global_gridx(i, gridx);

	    //int lhs = i - 1;   // x direction
	    int rhs = i + 1;  // x direction
	    //int bwd = i - kinetics->nx;  // y direction
	    int fwd = i + kinetics->nx;  // y direction
	    //int down = i - kinetics->nx * kinetics->ny; // z direction
	    int up = i + kinetics->nx * kinetics->ny;  // z direction
	    //printf("%i l=%i r=%i b=%i f=%i d=%i u=%i x=%e y=%e z=%e, %e \n", i, lhs, rhs, bwd, fwd, down, up, gridx[0], gridx[1], gridx[2], stepy);
		//if (gridx[0] - stepx > xlo) fprintf(fp, "%i %i 1\n",i, lhs);
		if (gridx[0] + stepx < xhi) fprintf(fp, "%i %i 1\n",i, rhs);
		//if (gridx[1] - stepy > ylo) fprintf(fp, "%i %i 1\n",i, bwd);
		if (gridx[1] + stepy < yhi) fprintf(fp, "%i %i 1\n",i, fwd);
		//if (gridx[2] - stepz > zlo) fprintf(fp, "%i %i 1\n",i, down);
		if (gridx[2] + stepz < zhi) fprintf(fp, "%i %i 1\n",i, up);
	  }
  }
}

int DumpBioSSTL::get_global_id(int i, double* gridx)
{
    int gxpos = gridx[0] / stepx;
    int gypos = gridx[1] / stepy;
    int gzpos = gridx[2] / stepz;

    int pos = gxpos + gypos * kinetics->nx + gzpos * kinetics->nx * kinetics->ny;

	return pos;
}

double* DumpBioSSTL::get_grid_x(int i, double* gridx) {

    int zpos = i / (kinetics->subn[0] * kinetics->subn[1]);
    int ypos = (i - zpos * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0];
    int xpos = i - zpos * (kinetics->subn[0] * kinetics->subn[1]) - ypos * kinetics->subn[0];

    xpos++;
    ypos++;
    zpos++;

    gridx[0] = kinetics->sublo[0] + xpos * stepx - stepx/2;
    gridx[1] = kinetics->sublo[1] + ypos * stepy - stepy/2;
    gridx[2] = kinetics->sublo[2] + zpos * stepz - stepz/2;

	return gridx;
}

double* DumpBioSSTL::get_global_gridx(int i, double* gridx) {

    int zpos = i / (kinetics->nx * kinetics->ny);
    int ypos = (i - zpos * (kinetics->nx * kinetics->ny)) / kinetics->nx;
    int xpos = i - zpos * (kinetics->nx * kinetics->ny) - ypos * kinetics->nx;

    xpos++;
    ypos++;
    zpos++;

    gridx[0] = xpos * stepx - stepx/2;
    gridx[1] = ypos * stepy - stepy/2;
    gridx[2] = zpos * stepz - stepz/2;

	return gridx;
}

void DumpBioSSTL::write_nus_data(int nuID)
{
  int ngrids = kinetics->ngrids;
  int nxyz = kinetics->nx*kinetics->ny*kinetics->nz;

  int *local_id = new int[ngrids];
  int *local_size = new int[comm->nprocs];

  int *full_id = new int[nxyz];
  double *full_s = new double[nxyz];

  for(int i = 0; i < kinetics->ngrids; i++){
    double gridx[3] = {0,0,0};
  	get_grid_x(i, gridx);
  	local_id[i] = get_global_id(i, gridx);
  }

  MPI_Gather(&ngrids, 1, MPI_INT, local_size, 1, MPI_INT, 0, world);

  // send grid id
  if(comm->me != 0) {
	  MPI_Request request;
	  MPI_Send(local_id, ngrids, MPI_INT, 0, comm->me, world);
  } else {
	 for(int i = 0; i < kinetics->ngrids; i++){
		 full_id[i] = local_id[i];
	 }

	 int start = ngrids;
	 for(int i = 1; i < comm->nprocs; i++){
		 MPI_Recv(&full_id[start], local_size[i], MPI_INT, i, i, world, MPI_STATUS_IGNORE);
		 start += local_size[i];
	 }
  }

  // send nus
  double *nus = kinetics->nus[nuID];
  if(comm->me != 0) {
	  MPI_Request request;
	  MPI_Send(nus, ngrids, MPI_DOUBLE, 0, comm->me, world);
  } else {
	 for(int i = 0; i < kinetics->ngrids; i++){
		 full_s[i] = nus[i];
	 }

	 int start = ngrids;
	 for(int i = 1; i < comm->nprocs; i++){
		 MPI_Recv(&full_s[start], local_size[i], MPI_DOUBLE, i, i, world, MPI_STATUS_IGNORE);
		 start += local_size[i];
	 }
  }

  if (comm->me == 0) {
	 fprintf(fp, "Time %f\n",get_time());
	 for(int i = 0; i < nxyz; i++){
		fprintf(fp, "%i, %e\n",full_id[i]+1, full_s[i]);
	 }
  }

  delete local_id;
  delete local_size;
  delete full_id;
  delete full_s;
}

void DumpBioSSTL::write_volf_data(int t)
{

  int ngrids = kinetics->ngrids;
  int nxyz = kinetics->nx*kinetics->ny*kinetics->nz;

  int *local_id = new int[ngrids];
  int *local_size = new int[comm->nprocs];

  int *full_id = new int[nxyz];
  double *full_vf = new double[nxyz];

  for(int i = 0; i < kinetics->ngrids; i++){
	double gridx[3] = {0,0,0};
	get_grid_x(i, gridx);
	local_id[i] = get_global_id(i, gridx);
  }

  MPI_Gather(&ngrids, 1, MPI_INT, local_size, 1, MPI_INT, 0, world);

  // send grid id
  if(comm->me != 0) {
	  MPI_Request request;
	  MPI_Send(local_id, ngrids, MPI_INT, 0, comm->me, world);
  } else {
	 for(int i = 0; i < kinetics->ngrids; i++){
		 full_id[i] = local_id[i];
	 }

	 int start = ngrids;
	 for(int i = 1; i < comm->nprocs; i++){
		 MPI_Recv(&full_id[start], local_size[i], MPI_INT, i, i, world, MPI_STATUS_IGNORE);
		 start += local_size[i];
	 }
  }

  // send volume fraction
  if(comm->me != 0) {
	  MPI_Request request;
	  MPI_Send(vol_frac[t], ngrids, MPI_DOUBLE, 0, comm->me, world);
  } else {
	 for(int i = 0; i < kinetics->ngrids; i++){
		 full_vf[i] = vol_frac[t][i];
	 }

	 int start = ngrids;
	 for(int i = 1; i < comm->nprocs; i++){
		 MPI_Recv(&full_vf[start], local_size[i], MPI_DOUBLE, i, i, world, MPI_STATUS_IGNORE);
		 start += local_size[i];
	 }
  }

  if (comm->me == 0) {
	 fprintf(fp, "Time %f\n",get_time());
	 for(int i = 0; i < nxyz; i++){
		fprintf(fp, "%i, %e\n",full_id[i]+1, full_vf[i]);
	 }
  }

  delete local_id;
  delete local_size;
  delete full_id;
  delete full_vf;
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_biomass_data()
{
  if (!mass_header) {
	fprintf(fp, "Time ");
    for(int i = 0; i < atom->ntypes+1; i++){
      if(i == 0) fprintf(fp, "Total ");
      else fprintf(fp, "%s ", kinetics->bio->tname[i]);
    }
    mass_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%f, ", get_time());

  for(int i = 0; i < atom->ntypes+1; i++){
    if(i != atom->ntypes) fprintf(fp, "%e, ", cmass->vector[i]);
    else fprintf(fp, "%e", cmass->vector[i]);
  }

  fprintf(fp, "\n");
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_height_data()
{
  fprintf(fp, "%f, %e \n", get_time(), cheight->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_pressure_data()
{
  fprintf(fp, "%f, %e \n", get_time(), cpres->scalar);
}


/* ---------------------------------------------------------------------- */

void DumpBioSSTL::write_ntype_data()
{
  if (!type_header) {
    fprintf(fp, "Time ");
    for(int i = 1; i < atom->ntypes+1; i++){
      fprintf(fp, "%s,", kinetics->bio->tname[i]);
    }
    type_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%f, ", get_time());

  for(int i = 1; i < atom->ntypes+1; i++){
	  if(i != atom->ntypes) fprintf(fp, "%i, ", (int)ctype->vector[i]);
	  else fprintf(fp, "%i", (int)ctype->vector[i]);
  }

  fprintf(fp, "\n");
}

double DumpBioSSTL::get_time()
{
	return (update->ntimestep*update->dt)/nevery;
}

void DumpBioSSTL::update_volf() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double vol = stepx * stepy * stepz;

  for (int i = 0; i <= atom->ntypes; i++) {
	for (int j = 0; j < kinetics->ngrids; j++) {
	  vol_frac[i][j] = 0;
	}
  }

  for (int i = 0; i < nlocal; i++) {
	int pos = kinetics->position(i);
	int t = atom->type[i];
	double volf = (4*(MY_PI*atom->radius[i]*atom->radius[i]*atom->radius[i])/3) / vol;
	vol_frac[t][pos] += volf;
	vol_frac[0][pos] += volf;
  }
}

bigint DumpBioSSTL::memory_usage() {
  return 0;
}



