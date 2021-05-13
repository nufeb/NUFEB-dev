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

#include "dump_grid_vtk.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"
#include "grid.h"
#include "grid_masks.h"
#include "domain.h"
#include "group.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLImageDataWriter.h>

#include <regex>
#include <sstream>

using namespace LAMMPS_NS;

DumpGridVTK::DumpGridVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  filewriter = 0;
  parse_fields(narg, arg);
}

DumpGridVTK::~DumpGridVTK() {}

void DumpGridVTK::init_style() {
  using std::placeholders::_1;

  if (!grid)
    error->all(FLERR, "No grid defined for dump grid/vtk");
  
  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "con") {
      packs.push_back(std::bind(&DumpGridVTK::pack_concentration, this, _1));
    } else if (*it == "rea") {
      packs.push_back(std::bind(&DumpGridVTK::pack_reaction, this, _1));
    } else if (*it == "den") {
      packs.push_back(std::bind(&DumpGridVTK::pack_density, this, _1));
    } else if (*it == "gro") {
      packs.push_back(std::bind(&DumpGridVTK::pack_growth, this, _1));
    }
  }
}

void DumpGridVTK::write() {
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(grid->subbox[0] - 1, grid->subbox[1] - 1, grid->subbox[2] - 1);
  image->SetSpacing(grid->cell_size, grid->cell_size, grid->cell_size);
  double origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = grid->cell_size * (grid->sublo[i] + 1) - domain->boxlo[i];
  image->SetOrigin(origin[0], origin[1], origin[2]);

  for (auto it = packs.begin(); it != packs.end(); ++it) {
    (*it)(image);
  }

  // filename must contain '%' and '*'
  std::string str(filename);
  auto perc = str.find('%');
  if (perc == std::string::npos) {
    error->all(FLERR, "dump grid filename must contain '%' special character");
  }
  auto star = str.find('*');
  if (star == std::string::npos) {
    error->all(FLERR, "dump grid filename must contain '*' special character");
  }

  str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
  str = std::regex_replace(str, std::regex("\\*"), std::to_string(update->ntimestep));

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(str.c_str());
  writer->SetInputData(image);
  writer->Write();
}

int DumpGridVTK::parse_fields(int narg, char **arg) {
  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    if (strcmp(arg[iarg],"con") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "rea") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "den") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "gro") == 0) {
      fields.push_back(arg[iarg]);
    }
  }
  return i;
}

void DumpGridVTK::pack_concentration(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "concentration", grid->conc, grid->sub_names, grid->nsubs);
}

void DumpGridVTK::pack_reaction(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "reaction", grid->reac, grid->sub_names, grid->nsubs);
}

void DumpGridVTK::pack_density(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "density", grid->dens, group->names, group->ngroup);
}

void DumpGridVTK::pack_growth(vtkSmartPointer<vtkImageData> image) {
  pack_tuple<2>(image, "growth", grid->growth, group->names, group->ngroup);
}

void DumpGridVTK::pack_tuple1(vtkSmartPointer<vtkImageData> image, const char *name, double *data) {
  vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK))
      array->InsertNextTuple1(data[i]);
  }
  image->GetCellData()->AddArray(array);
}

void DumpGridVTK::pack_tuple1(vtkSmartPointer<vtkImageData> image, const char *name, double **data, char **names, int count) {
  for (int n = 0; n < count; n++) {
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    std::ostringstream oss;
    oss << names[n] << " " << name;
    array->SetName(oss.str().c_str());
    array->SetNumberOfComponents(1);
    for (int i = 0; i < grid->ncells; i++) {
      if (!(grid->mask[i] & GHOST_MASK))
	array->InsertNextTuple1(data[n][i]);
    }
    image->GetCellData()->AddArray(array);
  }
}

template <int N>
void DumpGridVTK::pack_tuple(vtkSmartPointer<vtkImageData> image, const char *name, double ***data, char **names, int count) {
  for (int n = 0; n < count; n++) {
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    std::ostringstream oss;
    oss << names[n] << " " << name;
    array->SetName(oss.str().c_str());
    array->SetNumberOfComponents(N);
    for (int i = 0; i < grid->ncells; i++) {
      if (!(grid->mask[i] & GHOST_MASK)) {
	double tuple[N];
	for (int j = 0; j < N; j++)
	  tuple[j] = data[n][i][j];
	array->InsertNextTuple(tuple);
      }
    }
    image->GetCellData()->AddArray(array);
  }
}
