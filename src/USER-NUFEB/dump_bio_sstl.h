/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(bio/sstl,DumpBioSSTL)

#else

#ifndef LMP_DUMP_BIO_SSTL_H
#define LMP_DUMP_BIO_SSTL_H

#include <mpi.h>
#include "dump.h"

namespace LAMMPS_NS {

class DumpBioSSTL : public Dump {
 public:
  DumpBioSSTL(LAMMPS *, int, char**);
  virtual ~DumpBioSSTL();

 protected:
  int nevery;                // dump frequency for output
  //gzFile gzFp;  // file pointer for the compressed output stream

  FILE *fp;                  // file to write dump to
  char *filename;            // user-specified file
  int nkeywords;

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  int nus_flag, mass_flag, pres_flag, height_flag, ntypes_flag, volfrac_flag, gridx_flag;
  int mass_header, type_header;

  double **vol_frac;

  double xlo,xhi,ylo,yhi,zlo,zhi;
  double stepx, stepy, stepz;
  int nx, ny, nz;               // # of grids in x, y and z

  char **keywords;
  class FixKinetics *kinetics;
  class BIO *bio;

  int ncompute;                  // # of Compute objects used by dump

  class ComputeNufebHeight *cheight;
  class ComputeNufebNtypes *ctype;
  class ComputeNufebBiomass *cmass;
  class ComputeNufebPressure *cpres;

  void write();
  void init_style();
  void write_header(bigint);
  void pack(tagint *);
  void write_data(int, double *);

  void write_nus_data(int);
  void write_biomass_data();
  void write_ntype_data();

  void write_pressure_data();
  void write_height_data();
  void write_volf_data(int);
  double get_gridx(int);
  void write_sstl_model();
  void write_gridx_data();

  double get_time();
  int get_global_id(int, double*);
  double* get_grid_x(int, double*);
  double* get_global_gridx(int, double*);
  void update_volf();

  bigint memory_usage();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump movie filename

The file produced by dump movie cannot be binary or compressed
and must be a single file for a single processor.

E: Support for writing movies not included

LAMMPS was not built with the -DLAMMPS_FFMPEG switch in the Makefile

E: Failed to open FFmpeg pipeline to file %s

The specified file cannot be opened.  Check that the path and name are
correct and writable and that the FFmpeg executable can be found and run.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
