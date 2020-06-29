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

#include "force.h"
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include "style_bond.h"
#include "style_angle.h"
#include "style_dihedral.h"
#include "style_improper.h"
#include "style_pair.h"
#include "style_kspace.h"
#include "atom.h"
#include "comm.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_hybrid_overlay.h"
#include "bond.h"
#include "bond_hybrid.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"
#include "update.h"
#include "utils.h"
#include "fmt/format.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

Force::Force(LAMMPS *lmp) : Pointers(lmp)
{
  newton = newton_pair = newton_bond = 1;

  special_lj[0] = special_coul[0] = 1.0;
  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;
  special_extra = 0;

  dielectric = 1.0;
  qqr2e_lammps_real = 332.06371;          // these constants are toggled
  qqr2e_charmm_real = 332.0716;           // by new CHARMM pair styles

  pair = NULL;
  bond = NULL;
  angle = NULL;
  dihedral = NULL;
  improper = NULL;
  kspace = NULL;

  char *str = (char *) "none";
  int n = strlen(str) + 1;
  pair_style = new char[n];
  strcpy(pair_style,str);
  bond_style = new char[n];
  strcpy(bond_style,str);
  angle_style = new char[n];
  strcpy(angle_style,str);
  dihedral_style = new char[n];
  strcpy(dihedral_style,str);
  improper_style = new char[n];
  strcpy(improper_style,str);
  kspace_style = new char[n];
  strcpy(kspace_style,str);

  pair_restart = NULL;
  create_factories();
}

void _noopt Force::create_factories()
{
  // fill pair map with pair styles listed in style_pair.h

  pair_map = new PairCreatorMap();

#define PAIR_CLASS
#define PairStyle(key,Class) \
  (*pair_map)[#key] = &pair_creator<Class>;
#include "style_pair.h"
#undef PairStyle
#undef PAIR_CLASS

  bond_map = new BondCreatorMap();

#define BOND_CLASS
#define BondStyle(key,Class) \
  (*bond_map)[#key] = &bond_creator<Class>;
#include "style_bond.h"
#undef BondStyle
#undef BOND_CLASS

  angle_map = new AngleCreatorMap();

#define ANGLE_CLASS
#define AngleStyle(key,Class) \
  (*angle_map)[#key] = &angle_creator<Class>;
#include "style_angle.h"
#undef AngleStyle
#undef ANGLE_CLASS

  dihedral_map = new DihedralCreatorMap();

#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) \
  (*dihedral_map)[#key] = &dihedral_creator<Class>;
#include "style_dihedral.h"
#undef DihedralStyle
#undef DIHEDRAL_CLASS

  improper_map = new ImproperCreatorMap();

#define IMPROPER_CLASS
#define ImproperStyle(key,Class) \
  (*improper_map)[#key] = &improper_creator<Class>;
#include "style_improper.h"
#undef ImproperStyle
#undef IMPROPER_CLASS

  kspace_map = new KSpaceCreatorMap();

#define KSPACE_CLASS
#define KSpaceStyle(key,Class) \
  (*kspace_map)[#key] = &kspace_creator<Class>;
#include "style_kspace.h"
#undef KSpaceStyle
#undef KSPACE_CLASS
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
  delete [] pair_style;
  delete [] bond_style;
  delete [] angle_style;
  delete [] dihedral_style;
  delete [] improper_style;
  delete [] kspace_style;

  delete [] pair_restart;

  if (pair) delete pair;
  if (bond) delete bond;
  if (angle) delete angle;
  if (dihedral) delete dihedral;
  if (improper) delete improper;
  if (kspace) delete kspace;

  pair = NULL;
  bond = NULL;
  angle = NULL;
  dihedral = NULL;
  improper = NULL;
  kspace = NULL;

  delete pair_map;
  delete bond_map;
  delete angle_map;
  delete dihedral_map;
  delete improper_map;
  delete kspace_map;
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
  qqrd2e = qqr2e/dielectric;

  // check if pair style must be specified after restart
  if (pair_restart) {
    if (!pair)
      error->all(FLERR,fmt::format("Must re-specify non-restarted pair style "
                                   "({}) after read_restart", pair_restart));
  }

  if (kspace) kspace->init();         // kspace must come before pair
  if (pair) pair->init();             // so g_ewald is defined
  if (bond) bond->init();
  if (angle) angle->init();
  if (dihedral) dihedral->init();
  if (improper) improper->init();

  // print warnings if topology and force field are inconsistent

  if (comm->me == 0) {
    if (!bond && (atom->nbonds > 0)) {
      error->warning(FLERR,"Bonds are defined but no bond style is set");
      if ((special_lj[1] != 1.0) || (special_coul[1] != 1.0))
        error->warning(FLERR,"Likewise 1-2 special neighbor interactions != 1.0");
    }
    if (!angle && (atom->nangles > 0)) {
      error->warning(FLERR,"Angles are defined but no angle style is set");
      if ((special_lj[2] != 1.0) || (special_coul[2] != 1.0))
        error->warning(FLERR,"Likewise 1-3 special neighbor interactions != 1.0");
    }
    if (!dihedral && (atom->ndihedrals > 0)) {
      error->warning(FLERR,"Dihedrals are defined but no dihedral style is set");
      if ((special_lj[3] != 1.0) || (special_coul[3] != 1.0))
        error->warning(FLERR,"Likewise 1-4 special neighbor interactions != 1.0");
    }
    if (!improper && (atom->nimpropers > 0))
      error->warning(FLERR,"Impropers are defined but no improper style is set");
  }
}

/* ---------------------------------------------------------------------- */

void Force::setup()
{
  if (pair) pair->setup();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_pair(const std::string &style, int trysuffix)
{
  delete [] pair_style;
  if (pair) delete pair;
  if (pair_restart) delete [] pair_restart;
  pair_style = NULL;
  pair = NULL;
  pair_restart = NULL;

  int sflag;
  pair = new_pair(style,trysuffix,sflag);
  store_style(pair_style,style,sflag);
}

/* ----------------------------------------------------------------------
   generate a pair class
   if trysuffix = 1, try first with suffix1/2 appended
   return sflag = 0 for no suffix added, 1 or 2 for suffix1/2 added
------------------------------------------------------------------------- */

Pair *Force::new_pair(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator pair_creator = (*pair_map)[estyle];
        return pair_creator(lmp);
      }
    }
    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator pair_creator = (*pair_map)[estyle];
        return pair_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("pair",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per pair style in style_pair.h
------------------------------------------------------------------------- */

template <typename T>
Pair *Force::pair_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to Pair class if matches word or matches hybrid sub-style
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   if nsub > 0, match Nth hybrid sub-style
   return NULL if no match or if nsub=0 and multiple sub-styles match
------------------------------------------------------------------------- */

Pair *Force::pair_match(const std::string &word, int exact, int nsub)
{
  int iwhich,count;

  if (exact && (word == pair_style)) return pair;
  else if (!exact && utils::strmatch(pair_style,word)) return pair;
  else if (utils::strmatch(pair_style,"^hybrid")) {
    PairHybrid *hybrid = (PairHybrid *) pair;
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++)
      if ((exact && (word == hybrid->keywords[i])) ||
          (!exact && utils::strmatch(hybrid->keywords[i],word))) {
        iwhich = i;
        count++;
        if (nsub == count) return hybrid->styles[iwhich];
      }
    if (count == 1) return hybrid->styles[iwhich];
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   return style name of Pair class that matches Pair ptr
   called by Neighbor::print_neigh_info()
   return NULL if no match
------------------------------------------------------------------------- */

char *Force::pair_match_ptr(Pair *ptr)
{
  if (ptr == pair) return pair_style;

  if (utils::strmatch(pair_style,"^hybrid")) {
    PairHybrid *hybrid = (PairHybrid *) pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (ptr == hybrid->styles[i]) return hybrid->keywords[i];
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   create a bond style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_bond(const std::string &style, int trysuffix)
{
  delete [] bond_style;
  if (bond) delete bond;

  int sflag;
  bond = new_bond(style,trysuffix,sflag);
  store_style(bond_style,style,sflag);
}

/* ----------------------------------------------------------------------
   generate a bond class, fist with suffix appended
------------------------------------------------------------------------- */

Bond *Force::new_bond(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (bond_map->find(estyle) != bond_map->end()) {
        BondCreator bond_creator = (*bond_map)[estyle];
        return bond_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (bond_map->find(estyle) != bond_map->end()) {
        BondCreator bond_creator = (*bond_map)[estyle];
        return bond_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (bond_map->find(style) != bond_map->end()) {
    BondCreator bond_creator = (*bond_map)[style];
    return bond_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("bond",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per bond style in style_bond.h
------------------------------------------------------------------------- */

template <typename T>
Bond *Force::bond_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to current bond class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Bond *Force::bond_match(const std::string &style)
{
  if (style == bond_style) return bond;
  else if (strcmp(bond_style,"hybrid") == 0) {
    BondHybrid *hybrid = (BondHybrid *) bond;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create an angle style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_angle(const std::string &style, int trysuffix)
{
  delete [] angle_style;
  if (angle) delete angle;

  int sflag;
  angle = new_angle(style,trysuffix,sflag);
  store_style(angle_style,style,sflag);
}

/* ----------------------------------------------------------------------
   generate an angle class
------------------------------------------------------------------------- */

Angle *Force::new_angle(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (angle_map->find(estyle) != angle_map->end()) {
        AngleCreator angle_creator = (*angle_map)[estyle];
        return angle_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix;
      if (angle_map->find(estyle) != angle_map->end()) {
        AngleCreator angle_creator = (*angle_map)[estyle];
        return angle_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (angle_map->find(style) != angle_map->end()) {
    AngleCreator angle_creator = (*angle_map)[style];
    return angle_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("angle",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per angle style in style_angle.h
------------------------------------------------------------------------- */

template <typename T>
Angle *Force::angle_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to current angle class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Angle *Force::angle_match(const std::string &style)
{
  if (style == angle_style) return angle;
  else if (utils::strmatch(angle_style,"^hybrid")) {
    AngleHybrid *hybrid = (AngleHybrid *) angle;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create a dihedral style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_dihedral(const std::string &style, int trysuffix)
{
  delete [] dihedral_style;
  if (dihedral) delete dihedral;

  int sflag;
  dihedral = new_dihedral(style,trysuffix,sflag);
  store_style(dihedral_style,style,sflag);
}

/* ----------------------------------------------------------------------
   generate a dihedral class
------------------------------------------------------------------------- */

Dihedral *Force::new_dihedral(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (dihedral_map->find(estyle) != dihedral_map->end()) {
        DihedralCreator dihedral_creator = (*dihedral_map)[estyle];
        return dihedral_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (dihedral_map->find(estyle) != dihedral_map->end()) {
        DihedralCreator dihedral_creator = (*dihedral_map)[estyle];
        return dihedral_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (dihedral_map->find(style) != dihedral_map->end()) {
    DihedralCreator dihedral_creator = (*dihedral_map)[style];
    return dihedral_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("dihedral",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per dihedral style in style_dihedral.h
------------------------------------------------------------------------- */

template <typename T>
Dihedral *Force::dihedral_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to current angle class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Dihedral *Force::dihedral_match(const std::string &style)
{
  if (style == dihedral_style) return dihedral;
  else if (utils::strmatch(dihedral_style,"^hybrid")) {
    DihedralHybrid *hybrid = (DihedralHybrid *) dihedral;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create an improper style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_improper(const std::string &style, int trysuffix)
{
  delete [] improper_style;
  if (improper) delete improper;

  int sflag;
  improper = new_improper(style,trysuffix,sflag);
  store_style(improper_style,style,sflag);
}

/* ----------------------------------------------------------------------
   generate a improper class
------------------------------------------------------------------------- */

Improper *Force::new_improper(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (improper_map->find(estyle) != improper_map->end()) {
        ImproperCreator improper_creator = (*improper_map)[estyle];
        return improper_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 2;
      std::string estyle = style + "/" + lmp->suffix2;
      if (improper_map->find(estyle) != improper_map->end()) {
        ImproperCreator improper_creator = (*improper_map)[estyle];
        return improper_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (improper_map->find(style) != improper_map->end()) {
    ImproperCreator improper_creator = (*improper_map)[style];
    return improper_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("improper",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per improper style in style_improper.h
------------------------------------------------------------------------- */

template <typename T>
Improper *Force::improper_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to current improper class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Improper *Force::improper_match(const std::string &style)
{
  if (style == improper_style) return improper;
  else if (utils::strmatch(improper_style,"^hybrid")) {
    ImproperHybrid *hybrid = (ImproperHybrid *) improper;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (style == hybrid->keywords[i]) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   new kspace style
------------------------------------------------------------------------- */

void Force::create_kspace(const std::string &style, int trysuffix)
{
  delete [] kspace_style;
  if (kspace) delete kspace;

  int sflag;
  kspace = new_kspace(style,trysuffix,sflag);
  store_style(kspace_style,style,sflag);

  if (comm->style == 1 && !kspace_match("ewald",0))
    error->all(FLERR,
               "Cannot yet use KSpace solver with grid with comm style tiled");
}

/* ----------------------------------------------------------------------
   generate a kspace class
------------------------------------------------------------------------- */

KSpace *Force::new_kspace(const std::string &style, int trysuffix, int &sflag)
{
  if (trysuffix && lmp->suffix_enable) {
    if (lmp->suffix) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix;
      if (kspace_map->find(estyle) != kspace_map->end()) {
        KSpaceCreator kspace_creator = (*kspace_map)[estyle];
        return kspace_creator(lmp);
      }
    }

    if (lmp->suffix2) {
      sflag = 1;
      std::string estyle = style + "/" + lmp->suffix2;
      if (kspace_map->find(estyle) != kspace_map->end()) {
        KSpaceCreator kspace_creator = (*kspace_map)[estyle];
        return kspace_creator(lmp);
      }
    }
  }

  sflag = 0;
  if (style == "none") return NULL;
  if (kspace_map->find(style) != kspace_map->end()) {
    KSpaceCreator kspace_creator = (*kspace_map)[style];
    return kspace_creator(lmp);
  }

  error->all(FLERR,utils::check_packages_for_style("kspace",style,lmp));

  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per kspace style in style_kspace.h
------------------------------------------------------------------------- */

template <typename T>
KSpace *Force::kspace_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to Kspace class if matches word
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   return NULL if no match
------------------------------------------------------------------------- */

KSpace *Force::kspace_match(const std::string &word, int exact)
{
  if (exact && (word == kspace_style)) return kspace;
  else if (!exact && utils::strmatch(kspace_style,word)) return kspace;
  return NULL;
}

/* ----------------------------------------------------------------------
   store style name in str allocated here
   if sflag = 0, no suffix
   if sflag = 1/2, append suffix or suffix2 to style
------------------------------------------------------------------------- */

void Force::store_style(char *&str, const std::string &style, int sflag)
{
  std::string estyle = style;

  if (sflag == 1) estyle += std::string("/") + lmp->suffix;
  else if (sflag == 2) estyle += std::string("/") + lmp->suffix2;

  str = new char[estyle.size()+1];
  strcpy(str,estyle.c_str());
}

/* ----------------------------------------------------------------------
   set special bond values
------------------------------------------------------------------------- */

void Force::set_special(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal special_bonds command");

  // defaults, but do not reset special_extra

  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"amber") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.5;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 5.0/6.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"charmm") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 0.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"dreiding") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"fene") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 1.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 1.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"lj/coul") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = special_coul[1] = numeric(FLERR,arg[iarg+1]);
      special_lj[2] = special_coul[2] = numeric(FLERR,arg[iarg+2]);
      special_lj[3] = special_coul[3] = numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"lj") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = numeric(FLERR,arg[iarg+1]);
      special_lj[2] = numeric(FLERR,arg[iarg+2]);
      special_lj[3] = numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"coul") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_coul[1] = numeric(FLERR,arg[iarg+1]);
      special_coul[2] = numeric(FLERR,arg[iarg+2]);
      special_coul[3] = numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal special_bonds command");
      if (strcmp(arg[iarg+1],"no") == 0) special_angle = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) special_angle = 1;
      else error->all(FLERR,"Illegal special_bonds command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal special_bonds command");
      if (strcmp(arg[iarg+1],"no") == 0) special_dihedral = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) special_dihedral = 1;
      else error->all(FLERR,"Illegal special_bonds command");
      iarg += 2;
    } else error->all(FLERR,"Illegal special_bonds command");
  }

  for (int i = 1; i <= 3; i++)
    if (special_lj[i] < 0.0 || special_lj[i] > 1.0 ||
        special_coul[i] < 0.0 || special_coul[i] > 1.0)
      error->all(FLERR,"Illegal special_bonds command");
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = nmin to nmax,
     (3) i* = i to nmax, (4) *j = nmin to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Force::bounds(const char *file, int line, char *str,
                   int nmax, int &nlo, int &nhi, int nmin)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = nmin;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < nmin || nhi > nmax || nlo > nhi)
    error->all(file,line,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = nmin to nmax,
     (3) i* = i to nmax, (4) *j = nmin to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Force::boundsbig(const char *file, int line, char *str,
                      bigint nmax, bigint &nlo, bigint &nhi, bigint nmin)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = ATOBIGINT(str);
  } else if (strlen(str) == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = nmin;
    nhi = ATOBIGINT(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = ATOBIGINT(str);
    nhi = nmax;
  } else {
    nlo = ATOBIGINT(str);
    nhi = ATOBIGINT(ptr+1);
  }

  if (nlo < nmin || nhi > nmax || nlo > nhi)
    error->all(file,line,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double Force::numeric(const char *file, int line, char *str)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0)
    error->all(file,line,"Expected floating point parameter instead of"
               " NULL or empty string in input script or data file");

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    error->all(file,line,fmt::format("Expected floating point parameter "
               "instead of '{}' in input script or data file",str));
  }

  return atof(str);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int Force::inumeric(const char *file, int line, char *str)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0)
    error->all(file,line,"Expected integer parameter instead of "
               "NULL or empty string in input script or data file");

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,fmt::format("Expected integer parameter instead "
               "of '{}' in input script or data file",str));
  }

  return atoi(str);
}

/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

bigint Force::bnumeric(const char *file, int line, char *str)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0)
    error->all(file,line,"Expected integer parameter instead of "
               "NULL or empty string in input script or data file");

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,fmt::format("Expected integer parameter instead "
               "of '{}' in input script or data file",str));
  }

  return ATOBIGINT(str);
}

/* ----------------------------------------------------------------------
   read a tag integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

tagint Force::tnumeric(const char *file, int line, char *str)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0)
    error->all(file,line,"Expected integer parameter instead of "
               "NULL or empty string in input script or data file");

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,fmt::format("Expected integer parameter instead "
               "of '{}' in input script or data file",str));
  }

  return ATOTAGINT(str);
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name
   if fails, search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */

FILE *Force::open_potential(const char *name, int *auto_convert)
{
  std::string filepath = utils::get_potential_file_path(name);

  if(!filepath.empty()) {
    std::string unit_style = update->unit_style;
    std::string date       = utils::get_potential_date(filepath, "potential");
    std::string units      = utils::get_potential_units(filepath, "potential");

    if(!date.empty()) {
      utils::logmesg(lmp, fmt::format("Reading potential file {} "
                                      "with DATE: {}\n", name, date));
    }

    if (auto_convert == nullptr) {
      if (!units.empty() && (units != unit_style)) {
        error->one(FLERR, fmt::format("Potential file {} requires {} units "
                                      "but {} units are in use", name, units,
                                      unit_style));
        return nullptr;
      }
    } else {
      if (units.empty() || units == unit_style) {
        *auto_convert = utils::NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real")
            && (*auto_convert & utils::METAL2REAL)) {
          *auto_convert = utils::METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal")
            && (*auto_convert & utils::REAL2METAL)) {
          *auto_convert = utils::REAL2METAL;
        } else {
          error->one(FLERR, fmt::format("Potential file {} requires {} units "
                                        "but {} units are in use", name,
                                        units, unit_style));
          return nullptr;
        }
      }
      if (*auto_convert != utils::NOCONVERT)
        lmp->error->warning(FLERR, fmt::format("Converting potential file in "
                                               "{} units to {} units",
                                               units, unit_style));
    }
    return fopen(filepath.c_str(), "r");
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of force classes
------------------------------------------------------------------------- */

bigint Force::memory_usage()
{
  bigint bytes = 0;
  if (pair) bytes += static_cast<bigint> (pair->memory_usage());
  if (bond) bytes += static_cast<bigint> (bond->memory_usage());
  if (angle) bytes += static_cast<bigint> (angle->memory_usage());
  if (dihedral) bytes += static_cast<bigint> (dihedral->memory_usage());
  if (improper) bytes += static_cast<bigint> (improper->memory_usage());
  if (kspace) bytes += static_cast<bigint> (kspace->memory_usage());
  return bytes;
}
