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

#include <cstring>
#include <ctime>
#include "nufeb_run_kokkos.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "atom_masks.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "compute_pressure.h"
#include "compute_ke.h"

// NUFEB specific

#include "grid_kokkos.h"
#include "grid_masks.h"
#include "comm_grid.h"
#include "fix_density.h"
#include "fix_diffusion_reaction.h"
#include "fix_monod.h"
#include "fix_eps_extract.h"
#include "fix_divide.h"
#include "fix_death.h"
#include "fix_reactor.h"
#include "fix_gas_liquid.h"
#include "fix_property.h"
#include "compute_volume.h"

using namespace LAMMPS_NS;

template<class ViewA, class ViewB>
struct ForceAdder {
  ViewA a;
  ViewB b;
  ForceAdder(const ViewA& a_, const ViewB& b_):a(a_),b(b_) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    a(i,0) += b(i,0);
    a(i,1) += b(i,1);
    a(i,2) += b(i,2);
  }
};

/* ---------------------------------------------------------------------- */

template<class View>
struct Zero {
  View v;
  Zero(const View &v_):v(v_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    v(i,0) = 0;
    v(i,1) = 0;
    v(i,2) = 0;
  }
};

/* ---------------------------------------------------------------------- */

NufebRunKokkos::NufebRunKokkos(LAMMPS *lmp, int narg, char **arg) :
  NufebRun(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void NufebRunKokkos::init()
{
  // this is required because many places check for verlet style
  delete [] update->integrate_style;
  update->integrate_style = new char[13];
  strcpy(update->integrate_style, "verlet/nufeb\0");

  // create fix nufeb/density
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"nufeb_density";
  fixarg[1] = (char *)"all";
  fixarg[2] = (char *)"nufeb/density";
  modify->add_fix(3, fixarg, 1);
  delete [] fixarg;
  fix_density = (FixDensity *)modify->fix[modify->nfix-1];

  // allocate space for storing fixes
  fix_monod = new FixMonod*[modify->nfix];
  fix_diffusion = new FixDiffusionReaction*[modify->nfix];
  fix_eps_extract = new FixEPSExtract*[modify->nfix];
  fix_divide = new FixDivide*[modify->nfix];
  fix_death = new FixDeath*[modify->nfix];
  fix_reactor = new FixReactor*[modify->nfix];
  fix_gas_liquid = new FixGasLiquid*[modify->nfix];
  fix_property = new FixProperty*[modify->nfix];

  // find fixes
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style, "nufeb/monod")) {
      fix_monod[nfix_monod++] = (FixMonod *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/diffusion_reaction")) {
      fix_diffusion[nfix_diffusion++] = (FixDiffusionReaction *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/eps_extract")) {
      fix_eps_extract[nfix_eps_extract++] = (FixEPSExtract *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/divide")) {
      fix_divide[nfix_divide++] = (FixDivide *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/death")) {
      fix_death[nfix_death++] = (FixDeath *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/reactor")) {
      fix_reactor[nfix_reactor++] = (FixReactor *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/gas_liquid")) {
      fix_gas_liquid[nfix_gas_liquid++] = (FixGasLiquid *)modify->fix[i];
    } else if (strstr(modify->fix[i]->style, "nufeb/property")) {
      fix_property[nfix_property++] = (FixProperty *)modify->fix[i];
    }
  }

  // create compute volume
  char **volarg = new char*[3];
  volarg[0] = (char *)"nufeb_volume";
  volarg[1] = (char *)"all";
  volarg[2] = (char *)"nufeb/volume/kk";
  modify->add_compute(3, volarg, 1);
  delete [] volarg;
  comp_volume = (ComputeVolume *)modify->compute[modify->ncompute-1];

  // create compute pressure
  char **pressarg = new char*[6];
  pressarg[0] = (char *)"nufeb_pressure";
  pressarg[1] = (char *)"all";
  pressarg[2] = (char *)"pressure";
  pressarg[3] = (char *)"NULL";
  pressarg[4] = (char *)"pair";
  pressarg[5] = (char *)"fix";
  modify->add_compute(6, pressarg, 1);
  delete [] pressarg;
  comp_pressure = (ComputePressure *)modify->compute[modify->ncompute-1];

  // create compute ke
  char **kearg = new char*[3];
  kearg[0] = (char *)"nufeb_ke";
  kearg[1] = (char *)"all";
  kearg[2] = (char *)"ke/kk";
  modify->add_compute(3, kearg, 1);
  delete [] kearg;
  comp_ke = (ComputeKE *)modify->compute[modify->ncompute-1];

  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_fdotr_compute via sum over ghosts

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present for clearing force arrays

  int ifix = modify->find_fix("package_omp");
  if (ifix >= 0) external_force_clear = 1;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;

  // check if simulation box is orthogonal
  if (domain->triclinic)
    error->all(FLERR, "Triclinic simulation box not supported in nufeb run style");
  triclinic = 0;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void NufebRunKokkos::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fprintf(screen,"Setting up NUFEB run ...\n");
    if (flag) {
      fprintf(screen,"  Unit style    : %s\n", update->unit_style);
      fprintf(screen,"  Current step  : " BIGINT_FORMAT "\n",
              update->ntimestep);
      fprintf(screen,"  Time step     : %g\n", update->dt);
      timer->print_timeout(screen);
    }
  }

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atomKK->sync(Host,ALL_MASK);
  atomKK->modified(Host,ALL_MASK);
  atomKK->setup();
  modify->setup_pre_exchange();
  atomKK->sync(Host,ALL_MASK);
  atomKK->modified(Host,ALL_MASK);
  if (triclinic) domain->x2lamda(atomKK->nlocal);
  domain->pbc();
  atomKK->sync(Host,ALL_MASK);
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atomKK->sortfreq > 0) atomKK->sort();
  comm->borders();
  if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);
  atomKK->sync(Host,ALL_MASK);
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  atomKK->modified(Host,ALL_MASK);
  neighbor->build(1);
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;

  // NUFEB specific

  grid->setup();
  comm_grid->setup();

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    timer->stamp(Timer::PAIR);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
    }
    timer->stamp(Timer::BOND);
  }

  if(force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      force->kspace->compute(eflag,vflag);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
      timer->stamp(Timer::KSPACE);
    } else force->kspace->compute_dummy(eflag,vflag);
  }
  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  output->setup(flag);
  lmp->kokkos->auto_sync = 0;
  update->setupflag = 0;

  // NUFEB specific

  biodt = update->dt;
  
  // disable all fixes that will be called directly
  fix_density->compute_flag = 0;
  disable_sync(fix_density);
  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->compute_flag = 0;
  }
  for (int i = 0; i < nfix_diffusion; i++) {
    fix_diffusion[i]->compute_flag = 0;
  }
  for (int i = 0; i < nfix_eps_extract; i++) {
    fix_eps_extract[i]->compute_flag = 0;
    disable_sync(fix_eps_extract[i]);
  }
  for (int i = 0; i < nfix_divide; i++) {
    fix_divide[i]->compute_flag = 0;
    disable_sync(fix_divide[i]);
  }
  for (int i = 0; i < nfix_death; i++) {
    fix_death[i]->compute_flag = 0;
    disable_sync(fix_death[i]);
  }
  for (int i = 0; i < nfix_reactor; i++)
    fix_reactor[i]->compute_flag = 0;
  for (int i = 0; i < nfix_gas_liquid; i++)
    fix_gas_liquid[i]->compute_flag = 0;
  for (int i = 0; i < nfix_property; i++)
    fix_property[i]->compute_flag = 0;

  // compute density
  fix_density->compute();

  // run diffusion until it reaches steady state
  if (init_diff_flag) {
    int niter = diffusion();
    if (comm->me == 0)
      fprintf(screen, "Initial diffusion reaction converged in %d steps\n", niter);
  } else {
    if (comm->me == 0)
      fprintf(screen, "Initial diffusion reaction convergence disabled\n");
  }
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void NufebRunKokkos::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    atomKK->sync(Host,ALL_MASK);
    atomKK->modified(Host,ALL_MASK);

    modify->setup_pre_exchange();

    atomKK->sync(Host,ALL_MASK);
    atomKK->modified(Host,ALL_MASK);

    if (triclinic) domain->x2lamda(atomKK->nlocal);
    domain->pbc();

    atomKK->sync(Host,ALL_MASK);

    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);

    atomKK->sync(Host,ALL_MASK);

    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();

    atomKK->modified(Host,ALL_MASK);

    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    timer->stamp(Timer::PAIR);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);


  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
    }
    timer->stamp(Timer::BOND);
  }

  if(force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      force->kspace->compute(eflag,vflag);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
      timer->stamp(Timer::KSPACE);
    } else force->kspace->compute_dummy(eflag,vflag);
  }

  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  lmp->kokkos->auto_sync = 0;
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void NufebRunKokkos::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_post_neighbor = modify->n_post_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atomKK->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  f_merge_copy = DAT::t_f_array("nufeb:f_merge_copy",atomKK->k_f.extent(0));

  atomKK->sync(Device,ALL_MASK);

  timer->init_timeout();
  for (int i = 0; i < n; i++) {

    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }
    ntimestep = ++update->ntimestep;

    // needs to come before ev_set
    comp_pressure->addstep(ntimestep);

    ev_set(ntimestep);

    double t = get_time();
    gridKK->sync(Host,ALL_MASK);
    atomKK->sync(Host,ALL_MASK);
    growth();
    atomKK->modified(Host,ALL_MASK);
    atomKK->sync(Device,ALL_MASK);
    if (profile)
      fprintf(profile, "%d %e ", update->ntimestep, get_time()-t);
    
    update->dt = pairdt;
    reset_dt();
    
    t = get_time();
    double vol = comp_volume->compute_scalar();
    npair = 0;
    double press = 0.0;
    do {
      // initial time integration

      //ktimer.reset();
      timer->stamp();
      modify->initial_integrate(vflag);
      //time += ktimer.seconds();
      if (n_post_integrate) modify->post_integrate();
      timer->stamp(Timer::MODIFY);

      // regular communication vs neighbor list rebuild

      nflag = neighbor->decide();

      if (nflag == 0) {
	timer->stamp();
	comm->forward_comm();
	timer->stamp(Timer::COMM);
      } else {
	// added debug
	//atomKK->sync(Host,ALL_MASK);
	//atomKK->modified(Host,ALL_MASK);

	if (n_pre_exchange) {
	  timer->stamp();
	  modify->pre_exchange();
	  timer->stamp(Timer::MODIFY);
	}
	// debug
	//atomKK->sync(Host,ALL_MASK);
	//atomKK->modified(Host,ALL_MASK);
	if (triclinic) domain->x2lamda(atomKK->nlocal);
	domain->pbc();
	if (domain->box_change) {
	  domain->reset_box();
	  comm->setup();
	  if (neighbor->style) neighbor->setup_bins();
	}
	timer->stamp();

	// added debug
	//atomKK->sync(Device,ALL_MASK);
	//atomKK->modified(Device,ALL_MASK);

	comm->exchange();
	if (sortflag && ntimestep >= atomKK->nextsort) atomKK->sort();
	comm->borders();

	// added debug
	//atomKK->sync(Host,ALL_MASK);
	//atomKK->modified(Host,ALL_MASK);

	if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);

	timer->stamp(Timer::COMM);
	if (n_pre_neighbor) {
	  modify->pre_neighbor();
	  timer->stamp(Timer::MODIFY);
	}
	neighbor->build(1);
	timer->stamp(Timer::NEIGH);
    	if (n_post_neighbor) {
    	  modify->post_neighbor();
    	  timer->stamp(Timer::MODIFY);
    	}
      }

      // force computations
      // important for pair to come before bonded contributions
      // since some bonded potentials tally pairwise energy/virial
      // and Pair:ev_tally() needs to be called before any tallying

      force_clear();

      timer->stamp();

      if (n_pre_force) {
	modify->pre_force(vflag);
	timer->stamp(Timer::MODIFY);
      }

      bool execute_on_host = false;
      unsigned int datamask_read_device = 0;
      unsigned int datamask_modify_device = 0;
      unsigned int datamask_read_host = 0;

      if ( pair_compute_flag ) {
	if (force->pair->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->pair->datamask_read;
	  datamask_modify_device |= force->pair->datamask_modify;
	} else {
	  datamask_read_device   |= force->pair->datamask_read;
	  datamask_modify_device |= force->pair->datamask_modify;
	}
      }
      if ( atomKK->molecular && force->bond )  {
	if (force->bond->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->bond->datamask_read;
	  datamask_modify_device |= force->bond->datamask_modify;
	} else {
	  datamask_read_device   |= force->bond->datamask_read;
	  datamask_modify_device |= force->bond->datamask_modify;
	}
      }
      if ( atomKK->molecular && force->angle ) {
	if (force->angle->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->angle->datamask_read;
	  datamask_modify_device |= force->angle->datamask_modify;
	} else {
	  datamask_read_device   |= force->angle->datamask_read;
	  datamask_modify_device |= force->angle->datamask_modify;
	}
      }
      if ( atomKK->molecular && force->dihedral ) {
	if (force->dihedral->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->dihedral->datamask_read;
	  datamask_modify_device |= force->dihedral->datamask_modify;
	} else {
	  datamask_read_device   |= force->dihedral->datamask_read;
	  datamask_modify_device |= force->dihedral->datamask_modify;
	}
      }
      if ( atomKK->molecular && force->improper ) {
	if (force->improper->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->improper->datamask_read;
	  datamask_modify_device |= force->improper->datamask_modify;
	} else {
	  datamask_read_device   |= force->improper->datamask_read;
	  datamask_modify_device |= force->improper->datamask_modify;
	}
      }
      if ( kspace_compute_flag ) {
	if (force->kspace->execution_space==Host) {
	  execute_on_host  = true;
	  datamask_read_host   |= force->kspace->datamask_read;
	  datamask_modify_device |= force->kspace->datamask_modify;
	} else {
	  datamask_read_device   |= force->kspace->datamask_read;
	  datamask_modify_device |= force->kspace->datamask_modify;
	}
      }

      if (pair_compute_flag) {
	atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
	atomKK->sync(force->pair->execution_space,~(~force->pair->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	Kokkos::Impl::Timer ktimer;
	force->pair->compute(eflag,vflag);
	atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
	atomKK->modified(force->pair->execution_space,~(~force->pair->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	timer->stamp(Timer::PAIR);
      }

      if(execute_on_host) {
        if(pair_compute_flag && force->pair->datamask_modify!=(F_MASK | ENERGY_MASK | VIRIAL_MASK))
          Kokkos::fence();
        atomKK->sync_overlapping_device(Host,~(~datamask_read_host|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
        if(pair_compute_flag && force->pair->execution_space!=Host) {
          Kokkos::deep_copy(LMPHostType(),atomKK->k_f.h_view,0.0);
        }
      }

      if (atomKK->molecular) {
	if (force->bond) {
	  atomKK->sync(force->bond->execution_space,~(~force->bond->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	  force->bond->compute(eflag,vflag);
	  atomKK->modified(force->bond->execution_space,~(~force->bond->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	}
	if (force->angle) {
	  atomKK->sync(force->angle->execution_space,~(~force->angle->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	  force->angle->compute(eflag,vflag);
	  atomKK->modified(force->angle->execution_space,~(~force->angle->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	}
	if (force->dihedral) {
	  atomKK->sync(force->dihedral->execution_space,~(~force->dihedral->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	  force->dihedral->compute(eflag,vflag);
	  atomKK->modified(force->dihedral->execution_space,~(~force->dihedral->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	}
	if (force->improper) {
	  atomKK->sync(force->improper->execution_space,~(~force->improper->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	  force->improper->compute(eflag,vflag);
	  atomKK->modified(force->improper->execution_space,~(~force->improper->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	}
	timer->stamp(Timer::BOND);
      }

      if (kspace_compute_flag) {
	atomKK->sync(force->kspace->execution_space,~(~force->kspace->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	force->kspace->compute(eflag,vflag);
	atomKK->modified(force->kspace->execution_space,~(~force->kspace->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
	timer->stamp(Timer::KSPACE);
      }

      if(execute_on_host && !std::is_same<LMPHostType,LMPDeviceType>::value) {
	if(f_merge_copy.extent(0)<atomKK->k_f.extent(0)) {
	  f_merge_copy = DAT::t_f_array("VerletKokkos::f_merge_copy",atomKK->k_f.extent(0));
	}
	f = atomKK->k_f.d_view;
	Kokkos::deep_copy(LMPHostType(),f_merge_copy,atomKK->k_f.h_view);
	Kokkos::parallel_for(atomKK->k_f.extent(0),
			     ForceAdder<DAT::t_f_array,DAT::t_f_array>(atomKK->k_f.d_view,f_merge_copy));
	atomKK->k_f.clear_sync_state(); // special case
	atomKK->k_f.modify<LMPDeviceType>();
      }

      if (n_pre_reverse) {
	modify->pre_reverse(eflag,vflag);
	timer->stamp(Timer::MODIFY);
      }

      // reverse communication of forces

      if (force->newton) {
	Kokkos::fence();
	comm->reverse_comm();
	timer->stamp(Timer::COMM);
      }

      // force modifications, final time integration, diagnostics

      if (n_post_force) modify->post_force(vflag);
      modify->final_integrate();
      if (n_end_of_step) modify->end_of_step();
      timer->stamp(Timer::MODIFY);

      ++npair;
      
      press = comp_pressure->compute_scalar() * domain->xprd * domain->yprd * domain->zprd;
      press += comp_ke->compute_scalar();
      press /= 3.0 * vol;
    } while(fabs(press) > pairtol && ((pairmax > 0) ? npair < pairmax : true));
    if (profile)
      fprintf(profile, "%d %e ", npair, get_time()-t);
    if (comm->me == 0) fprintf(screen, "pair interaction: %d steps (pressure %e N/m2)\n", npair, press);

    // update densities

    t = get_time();
    fix_density->compute();
    if (profile)
      fprintf(profile, "%e ", get_time()-t);


    // run diffusion until it reaches steady state

    t = get_time();
    ndiff = diffusion();
    if (profile)
      fprintf(profile, "%d %e\n", ndiff, get_time()-t);
    if (comm->me == 0) fprintf(screen, "diffusion: %d steps\n", ndiff);

    // all output

    if (ntimestep == output->next) {
      atomKK->sync(Host,ALL_MASK);
      gridKK->sync(Host,ALL_MASK);
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void NufebRunKokkos::force_clear()
{
  if (external_force_clear) return;

  atomKK->k_f.clear_sync_state(); // ignore host forces/torques since device views
  atomKK->k_torque.clear_sync_state(); //   will be cleared below

  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  if (neighbor->includegroup == 0) {
    int nall = atomKK->nlocal;
    if (force->newton) nall += atomKK->nghost;

    Kokkos::parallel_for(nall, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_f.view<LMPDeviceType>()));
    atomKK->modified(Device,F_MASK);

    if (torqueflag) {
      Kokkos::parallel_for(nall, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_torque.view<LMPDeviceType>()));
      atomKK->modified(Device,TORQUE_MASK);
    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    Kokkos::parallel_for(atomKK->nfirst, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_f.view<LMPDeviceType>()));
    atomKK->modified(Device,F_MASK);

    if (torqueflag) {
      Kokkos::parallel_for(atomKK->nfirst, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_torque.view<LMPDeviceType>()));
      atomKK->modified(Device,TORQUE_MASK);
    }

    if (force->newton) {
      auto range = Kokkos::RangePolicy<LMPDeviceType>(atomKK->nlocal, atomKK->nlocal + atomKK->nghost);
      Kokkos::parallel_for(range, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_f.view<LMPDeviceType>()));
      atomKK->modified(Device,F_MASK);

      if (torqueflag) {
	Kokkos::parallel_for(range, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_torque.view<LMPDeviceType>()));
	atomKK->modified(Device,TORQUE_MASK);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void NufebRunKokkos::growth()
{
  // set biological dt

  update->dt = biodt;
  reset_dt();

  // setup monod growth fixes for growth

  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->reaction_flag = 0;
    fix_monod[i]->growth_flag = 1;
  }

  // grow atoms

  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->compute();
  }

  for (int i = 0; i < nfix_eps_extract; i++) {
    fix_eps_extract[i]->compute();
  }

  for (int i = 0; i < nfix_divide; i++) {
    fix_divide[i]->compute();
  }

  for (int i = 0; i < nfix_death; i++) {
    fix_death[i]->compute();
  }

  for (int i = 0; i < nfix_property; i++) {
    fix_property[i]->compute();
  }
}

/* ---------------------------------------------------------------------- */

int NufebRunKokkos::diffusion()
{
  // setup monod growth fixes for diffusion

  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->reaction_flag = 1;
    fix_monod[i]->growth_flag = 0;
  }

  // set diffusion dt

  update->dt = diffdt;
  reset_dt();
  
  for (int i = 0; i < nfix_diffusion; i++) {
    fix_diffusion[i]->closed_system_init();
  }

  int niter = 0;
  bool flag;
  bool converge[nfix_diffusion];
  for (int i = 0; i < nfix_diffusion; i++)
    converge[i] = false;
  do {
    timer->stamp();
    comm_grid->forward_comm();
    timer->stamp(Timer::COMM);

    flag = true;
    for (int i = 0; i < nfix_diffusion; i++) {
      if (!converge[i])
	fix_diffusion[i]->compute_initial();
    }
    
    for (int i = 0; i < nfix_monod; i++) {
      fix_monod[i]->compute();
    }

    for (int i = 0; i < nfix_gas_liquid; i++) {
      fix_gas_liquid[i]->compute();
    }

    for (int i = 0; i < nfix_diffusion; i++) {
      if (!converge[i]) {
	fix_diffusion[i]->compute_final();
	double res = fix_diffusion[i]->compute_scalar();
	if (res < difftol) converge[i] = true;
	if (!converge[i]) flag = false;
      }
    }

    timer->stamp(Timer::MODIFY);
    ++niter;

    if (diffmax > 0 && niter >= diffmax)
      flag = true;
    
  } while (!flag);

  for (int i = 0; i < nfix_diffusion; i++) {
    fix_diffusion[i]->closed_system_scaleup(biodt);
  }

  return niter;
}

/* ---------------------------------------------------------------------- */

void NufebRunKokkos::reactor()
{
  // set biological dt

  update->dt = biodt;
  reset_dt();

  // update reactor attributes

  for (int i = 0; i < nfix_reactor; i++) {
    fix_reactor[i]->compute();
  }
}

/* ---------------------------------------------------------------------- */

void NufebRunKokkos::disable_sync(Fix *fix)
{
  fix->datamask_read = EMPTY_MASK;
  fix->datamask_modify = EMPTY_MASK;
  fix->kokkosable = 1;
}
