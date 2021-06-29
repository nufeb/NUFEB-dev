/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "dump_hdf5.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"
#include "grid.h"
#include "grid_masks.h"

#include <regex>
#include <sstream>

using namespace LAMMPS_NS;

DumpHDF5::DumpHDF5(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg) {
  parse_fields(narg, arg);
  setup();

  if (!multifile)
    create_one_file();
}

DumpHDF5::~DumpHDF5() {}


void DumpHDF5::setup()
{
  for (int i = 0; i < 3; i++) {
    subdims[i] = grid->subbox[i] - 2;
    substart[i] = grid->sublo[i] + 1;
    dims[i] = grid->box[i];
  }
  ncells = subdims[0] * subdims[1] * subdims[2];
}

hid_t DumpHDF5::create_filespace_atom(bool oneperproc) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  return H5Screate_simple(1, &dims, NULL); 
}

hid_t DumpHDF5::create_filespace_grid(bool oneperproc) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[i] = subdims[i];
  } else {
    for (int i = 0; i < 3; i++)
      dims[i] = this->dims[i];
  }
  return H5Screate_simple(3, dims, NULL);
}

void DumpHDF5::write() {
  hid_t file;
  hid_t proplist = H5P_DEFAULT;
  std::string str(filename);
  int offset = 0;
  bool oneperproc = false;

  auto perc = str.find('%');
  if (perc != std::string::npos) {
    oneperproc = true;
    str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
  }

  if (!oneperproc) {
    proplist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(proplist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }

  if (multifile) {
    str = std::regex_replace(str, std::regex("\\*"), std::to_string(update->ntimestep));
    file = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, proplist);
  } else {
    file = H5Fopen(str.c_str(), H5F_ACC_RDWR, proplist);
  }

  H5Pclose(proplist);

  if (!oneperproc)
    MPI_Scan(&atom->nlocal, &offset, 1, MPI_INT, MPI_SUM, world);

  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "id") {
      if (multifile) {
	write_atoms_scalar(file, "id", H5T_NATIVE_INT, atom->tag, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "id", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "id/" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_INT, atom->tag, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "type") {
      if (multifile) {
	write_atoms_scalar(file, "type", H5T_NATIVE_INT, atom->type, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "type", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "type/" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_INT, atom->type, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "x") {
      if (multifile) {
	write_atoms_comp(file, "x", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "x", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "x/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "y") {
      if (multifile) {
	write_atoms_comp(file, "y", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "y", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "y/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "z") {
      if (multifile) {
	write_atoms_comp(file, "z", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "z", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "z/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "vx") {
      if (multifile) {
         write_atoms_comp(file, "vx", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "vx", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vx/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "vy") {
      if (multifile) {
	 write_atoms_comp(file, "vy", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "vy", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vy/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "vz") {
      if (multifile) {
	 write_atoms_comp(file, "vz", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "vz", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vz/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "fx") {
      if (multifile) {
	write_atoms_comp(file, "fx", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "fx", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fx/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "fy") {
      if (multifile) {
	write_atoms_comp(file, "fy", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "fy", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fy/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "fz") {
      if (multifile) {
	write_atoms_comp(file, "fz", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "fz", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fz/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "radius") {
      if (multifile) {
	write_atoms_scalar(file, "radius", H5T_NATIVE_DOUBLE, atom->radius, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "radius", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "radius/" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->radius, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "conc") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "concentration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "concentration", H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "concentration/" << grid->sub_names[i];
	double *conc = new double[ncells];
	trim_grids(grid->conc[i], conc);
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, conc, oneperproc);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, conc, oneperproc);
	  H5Gclose(subgroup);
	}
	delete conc;
      }
      H5Gclose(group);
    } else if (*it == "reac") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "reaction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "reaction", H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "reaction/" << grid->sub_names[i];
	double *reac = new double[ncells];
	trim_grids(grid->reac[i], reac);
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, reac, oneperproc);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, reac, oneperproc);
	  H5Gclose(subgroup);
	}
	delete reac;
      }
      H5Gclose(group);
    } else if (*it == "grow") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "growth", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "growth", H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "growth/" << grid->sub_names[i];
	double *grow = new double[ncells];
	trim_grids(grid->growth[i][0], grow);
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, grow, oneperproc);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, grow, oneperproc);
	  H5Gclose(subgroup);
	}
	delete grow;
      }
      H5Gclose(group);
    }
  }
  H5Fclose(file);
}

void DumpHDF5::create_one_file() {
  std::string str(filename);
  hid_t proplist = H5P_DEFAULT;
  bool oneperproc = false;

  auto perc = str.find('%');
  if (perc != std::string::npos) {
    str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
    oneperproc = true;
  }

  if (!oneperproc) {
    proplist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(proplist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }
  hid_t file = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, proplist);
  H5Pclose(proplist);

  int offset = 0;
  MPI_Scan(&atom->nlocal, &offset, 1, MPI_INT, MPI_SUM, world);
  // create data structure
  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "id") {
      hid_t group = H5Gcreate(file, "id", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "type") {
      hid_t group = H5Gcreate(file, "type", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "x") {
      hid_t group = H5Gcreate(file, "x", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "y") {
      hid_t group = H5Gcreate(file, "y", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "z") {
      hid_t group = H5Gcreate(file, "z", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vx") {
      hid_t group = H5Gcreate(file, "vx", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vy") {
      hid_t group = H5Gcreate(file, "vy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vz") {
      hid_t group = H5Gcreate(file, "vz", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fx") {
      hid_t group = H5Gcreate(file, "fx", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fy") {
      hid_t group = H5Gcreate(file, "fy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fz") {
      hid_t group = H5Gcreate(file, "fz", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "radius") {
      hid_t group = H5Gcreate(file, "radius", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "conc") {
      hid_t group = H5Gcreate(file, "concentration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "concentration/" << grid->sub_names[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "reac") {
      hid_t group = H5Gcreate(file, "reaction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "reaction/" << grid->sub_names[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "grow") {
      hid_t group = H5Gcreate(file, "growth", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 0; i < grid->nsubs; i++) {
	std::ostringstream oss;
	oss << "growth/" << grid->sub_names[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    }
  }
  H5Fclose(file);
}

int DumpHDF5::parse_fields(int narg, char **arg) {
  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    if (strcmp(arg[iarg],"id") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "type") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "x") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "y") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "z") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vx") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vy") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vz") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fx") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fy") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fz") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "radius") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "conc") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "reac") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "grow") == 0) {
      fields.push_back(arg[iarg]);
    }
  }
  return i;
}

template <class T>
herr_t DumpHDF5::write_atoms_scalar(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  hid_t filespace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const hsize_t count = atom->nlocal;
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    const hsize_t start = offset - atom->nlocal; 
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hid_t memspace = H5Screate_simple(1, &count, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpHDF5::write_atoms_comp(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset, int comp) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  hid_t filespace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    const hsize_t start = offset - atom->nlocal; 
    const hsize_t count = atom->nlocal;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount[2];
  memcount[0] = atom->nlocal;
  memcount[1] = 3;
  hid_t memspace = H5Screate_simple(2, memcount, NULL);
  hsize_t memstart[2];
  memstart[0] = 0;
  memstart[1] = comp;
  memcount[0] = atom->nlocal;
  memcount[1] = 1;
  hsize_t stride[2];
  stride[0] = 1;
  stride[1] = 3;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, stride, memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpHDF5::write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = subdims[i];
  } else {
    for (int i = 0; i < 3; i++) {
      dims[2 - i] = this->dims[i];
    }
  }
  hid_t filespace = H5Screate_simple(3, dims, NULL); 
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    hsize_t start[3];
    for (int i = 0; i < 3; i++)
      start[2 - i] = substart[i];
    hsize_t count[3];
    for (int i = 0; i < 3; i++)
      count[2 - i] = subdims[i];
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount = ncells;
  hid_t memspace = H5Screate_simple(1, &memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpHDF5::write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int n, int index) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = subdims[i];
  } else {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = this->dims[i];
  }
  hid_t filespace = H5Screate_simple(3, dims, NULL); 
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    hsize_t start[3];
    for (int i = 0; i < 3; i++)
      start[2 - i] = substart[i];
    hsize_t count[3];
    for (int i = 0; i < 3; i++)
      count[2 - i] = subdims[i];
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount[2];
  memcount[0] = n;
  memcount[1] = ncells;
  hid_t memspace = H5Screate_simple(2, memcount, NULL);
  hsize_t memstart[2];
  memstart[0] = index;
  memstart[1] = 0;
  memcount[0] = 1;
  memcount[1] = ncells;
  hsize_t stride[2];
  stride[0] = 1;
  stride[1] = 1;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, stride, memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
void DumpHDF5::trim_grids(T *bufin, T *bufout) {
  int n = 0;
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      bufout[n] = bufin[i];
      n++;
    }
  }
}
