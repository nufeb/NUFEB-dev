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
#include "comm_grid_kokkos.h"
#include "grid_kokkos.h"
#include "grid_vec_kokkos.h"
#include "comm_kokkos.h"
#include "memory_kokkos.h"
#include "intersect_list.h"
#include "domain.h"
#include "error.h"
#include "lammps.h"
#include "kokkos.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CommGridKokkos::CommGridKokkos(LAMMPS *lmp) : CommGrid(lmp) {}

/* ---------------------------------------------------------------------- */

CommGridKokkos::~CommGridKokkos()
{
  memoryKK->destroy_kokkos(k_recv_begin, recv_begin);
  memoryKK->destroy_kokkos(k_recv_end, recv_end);
  memoryKK->destroy_kokkos(k_send_begin, send_begin);
  memoryKK->destroy_kokkos(k_send_end, send_end);
  memoryKK->destroy_kokkos(k_recv_cells, recv_cells);
  memoryKK->destroy_kokkos(k_send_cells, send_cells);
  memoryKK->destroy_kokkos(k_buf_recv, buf_recv);
  memoryKK->destroy_kokkos(k_buf_send, buf_send);
  memoryKK->destroy_kokkos(k_recv_cells_self, recv_cells_self);
  memoryKK->destroy_kokkos(k_send_cells_self, send_cells_self);
  memoryKK->destroy_kokkos(k_buf_self, buf_self);
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::init()
{
  CommGrid::init();

  memory->destroy(recv_begin);
  memory->destroy(recv_end);
  memory->destroy(send_begin);
  memory->destroy(send_end);
  
  memoryKK->create_kokkos(k_recv_begin, recv_begin, comm->nprocs, "comm_grid:recv_begin");
  memoryKK->create_kokkos(k_recv_end, recv_end, comm->nprocs, "comm_grid:recv_end");
  memoryKK->create_kokkos(k_send_begin, send_begin, comm->nprocs, "comm_grid:send_begin");
  memoryKK->create_kokkos(k_send_end, send_end, comm->nprocs, "comm_grid:send_end");
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::setup()
{
  CommGrid::setup();

  k_recv_begin.modify<LMPHostType>();
  k_recv_begin.sync<LMPDeviceType>();
  k_recv_end.modify<LMPHostType>();
  k_recv_end.sync<LMPDeviceType>();
  k_recv_cells.modify<LMPHostType>();
  k_recv_cells.sync<LMPDeviceType>();
  k_send_begin.modify<LMPHostType>();
  k_send_begin.sync<LMPDeviceType>();
  k_send_end.modify<LMPHostType>();
  k_send_end.sync<LMPDeviceType>();
  k_send_cells.modify<LMPHostType>();
  k_send_cells.sync<LMPDeviceType>();
  k_recv_cells_self.modify<LMPHostType>();
  k_recv_cells_self.sync<LMPDeviceType>();
  k_send_cells_self.modify<LMPHostType>();
  k_send_cells_self.sync<LMPDeviceType>();
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::forward_comm()
{
  if (!lmp->kokkos->forward_comm_classic) {
    if (lmp->kokkos->forward_comm_on_host) forward_comm_device<LMPHostType>();
    else forward_comm_device<LMPDeviceType>();
    return;
  }

  gridKK->sync(Host, ALL_MASK); // maybe make it gvec->sync_forward instead
  CommGrid::forward_comm();
  gridKK->modified(Host, ALL_MASK); // maybe make it gvec->modified_forward instead
}

template <class DeviceType>
void CommGridKokkos::forward_comm_device()
{
  GridVecKokkos *gvec = (GridVecKokkos *)grid->gvec;
  
  for (int p = 0; p < nrecvproc; p++) {
    MPI_Irecv(k_buf_recv.view<DeviceType>().data() + recv_begin[p] * size_forward,
	      (recv_end[p] - recv_begin[p]) * size_forward,
	      MPI_DOUBLE, recvproc[p], 0, world, &requests[p]);
  }
  for (int p = 0; p < nsendproc; p++) {
    int n = gvec->pack_comm_kokkos(send_begin[p], send_end[p], k_send_cells, k_buf_send);
    DeviceType().fence();
    MPI_Send(k_buf_send.view<DeviceType>().data() + send_begin[p] * size_forward, n, MPI_DOUBLE, sendproc[p], 0, world);
  }
  MPI_Waitall(nrecvproc, requests, MPI_STATUS_IGNORE);
  for (int p = 0; p < nrecvproc; p++) {
    gvec->unpack_comm_kokkos(recv_begin[p], recv_end[p], k_recv_cells, k_buf_recv);
  }
  gvec->pack_comm_kokkos(0, nsend_self, k_send_cells_self, k_buf_self);
  gvec->unpack_comm_kokkos(0, nrecv_self, k_recv_cells_self, k_buf_self);
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::migrate()
{
  int boxlo[3*comm->nprocs];
  MPI_Allgather(grid->sublo, 3, MPI_INT, boxlo, 3, MPI_INT, world);

  int boxhi[3*comm->nprocs];
  MPI_Allgather(grid->subhi, 3, MPI_INT, boxhi, 3, MPI_INT, world);

  int newsublo[3];  // new subgrid lower bound
  int newsubhi[3];  // new subgrid upper bound
  int newsubbox[3]; // new subgrid box
  const double small = 1e-12;
  for (int i = 0; i < 3; i++) {
    newsublo[i] = static_cast<int>((domain->sublo[i] - domain->boxlo[i]) /
				   grid->cell_size + small) - 1;
    newsubhi[i] = static_cast<int>((domain->subhi[i] - domain->boxlo[i]) /
				   grid->cell_size + small) + 1;
    newsubbox[i] = newsubhi[i] - newsublo[i];
  }

  int newboxlo[3*comm->nprocs];
  MPI_Allgather(newsublo, 3, MPI_INT, newboxlo, 3, MPI_INT, world);

  int newboxhi[3*comm->nprocs];
  MPI_Allgather(newsubhi, 3, MPI_INT, newboxhi, 3, MPI_INT, world);
  
  int *recv_begin_ = memory->create(recv_begin_, comm->nprocs, "comm_grid:recv_begin_");
  int *send_begin_ = memory->create(send_begin_, comm->nprocs, "comm_grid:send_begin_");
  int *recv_end_ = memory->create(recv_end_, comm->nprocs, "comm_grid:recv_end_");
  int *send_end_ = memory->create(send_end_, comm->nprocs, "comm_grid:send_end_");
      
  IntersectList recvlist(lmp, comm->nprocs);
  IntersectList sendlist(lmp, comm->nprocs);
  
  int nrecv_ = 0;
  int nsend_ = 0;
  int nrecv_self_ = 0;
  int nsend_self_ = 0;
  int nrecvproc_ = 0;
  int nsendproc_ = 0;
  int lo[3];
  int hi[3];
  for (int p = 0; p < comm->nprocs; p++) {
    // receiving from other procs and self
    int n = intersect(newsublo, newsubhi, &boxlo[3*p], &boxhi[3*p],
		      -1, -1, 0, 0, 0, lo, hi, true);
    if (n > 0) {
      recvlist.add(p, lo, hi, n);
      recv_begin_[nrecvproc_] = nrecv_;
      if (comm->me == p) nrecv_self_ += n;
      else nrecv_ += n;
      recv_end_[nrecvproc_++] = nrecv_;
    }
    // sending to other procs and self
    n = intersect(grid->sublo, grid->subhi, &newboxlo[3*p], &newboxhi[3*p],
		  -1, -1, 0, 0, 0, lo, hi, true);
    if (n > 0) {
      sendlist.add(p, lo, hi, n);
      send_begin_[nsendproc_] = nsend_;
      if (comm->me == p) nsend_self_ += n;
      else nsend_ += n;
      send_end_[nsendproc_++] = nsend_;
    }
  }  

  if (nrecv_self_ != nsend_self_)
    error->all(FLERR, "Conflicting self send and recv sizes (this is possibly a bug)");

  int *recv_cells_ = memory->create(recv_cells_, nrecv_, "comm_grid:recv_cells_");
  int *send_cells_ = memory->create(send_cells_, nsend_, "comm_grid:send_cells_");
  int *recv_cells_self_ = memory->create(recv_cells_self_, nrecv_self_, "comm_grid:recv_cells_self_");
  int *send_cells_self_ = memory->create(send_cells_self_, nsend_self_, "comm_grid:send_cells_self_");
  
  int irecv = 0;
  int irecv_self = 0;
  for (int i = 0; i < recvlist.n; i++) {
    for (int z = recvlist.boxlo[3*i+2]; z < recvlist.boxhi[3*i+2]; z++) {
      for (int y = recvlist.boxlo[3*i+1]; y < recvlist.boxhi[3*i+1]; y++) {
	for (int x = recvlist.boxlo[3*i]; x < recvlist.boxhi[3*i]; x++) {
	  if (comm->me != recvlist.procs[i]) {
	    recv_cells_[irecv++] = x - newsublo[0] +
	      (y - newsublo[1]) * newsubbox[0] +
	      (z - newsublo[2]) * newsubbox[0] * newsubbox[1];
	  } else {
	    recv_cells_self_[irecv_self++] = x - newsublo[0] +
	      (y - newsublo[1]) * newsubbox[0] +
	      (z - newsublo[2]) * newsubbox[0] * newsubbox[1];
	  }
  	}
      }
    }
  }

  int isend = 0;
  int isend_self = 0;
  for (int i = 0; i < sendlist.n; i++) {
    for (int z = sendlist.boxlo[3*i+2]; z < sendlist.boxhi[3*i+2]; z++) {
      for (int y = sendlist.boxlo[3*i+1]; y < sendlist.boxhi[3*i+1]; y++) {
	for (int x = sendlist.boxlo[3*i]; x < sendlist.boxhi[3*i]; x++) {
	  if (comm->me != sendlist.procs[i]) {
	    send_cells_[isend++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	  } else {
	    send_cells_self_[isend_self++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	  }
  	}
      }
    }
  }

  double *buf_recv_ = memory->create(buf_recv_, nrecv_ * size_exchange, "comm_grid:buf_recv_");
  double *buf_send_ = memory->create(buf_send_, nsend_ * size_exchange, "comm_grid:buf_send_");
  double *buf_self_ = memory->create(buf_self_, nrecv_self_ * size_exchange, "comm_grid:buf_self_");
  
  int nrequest = 0;
  for (int p = 0; p < recvlist.n; p++) {
    if (comm->me != recvlist.procs[p]) {
      MPI_Irecv(&buf_recv_[recv_begin_[p] * size_exchange],
		(recv_end_[p] - recv_begin_[p]) * size_exchange,
		MPI_DOUBLE, recvlist.procs[p], 0, world, &requests[nrequest++]);
    }
  }
  for (int p = 0; p < sendlist.n; p++) {
    if (comm->me != sendlist.procs[p]) {
      int n = grid->gvec->pack_exchange(send_end_[p] - send_begin_[p],
					&send_cells_[send_begin_[p]],
					buf_send_);
      MPI_Send(buf_send_, n, MPI_DOUBLE, sendlist.procs[p], 0, world);
    }
  }
  int n = grid->gvec->pack_exchange(nsend_self_, send_cells_self_, buf_self_);
  // safe to call grid::setup() here because all grid data are in local buffers
  grid->setup();
  MPI_Waitall(nrequest, requests, MPI_STATUS_IGNORE);
  for (int p = 0; p < recvlist.n; p++) {
    grid->gvec->unpack_exchange(recv_end_[p] - recv_begin_[p],
				&recv_cells_[recv_begin_[p]],
				&buf_recv_[recv_begin_[p] * size_exchange]);
  }

  grid->gvec->unpack_exchange(nrecv_self_, recv_cells_self_, buf_self_);

  memory->destroy(recv_begin_);
  memory->destroy(send_begin_);
  memory->destroy(recv_end_);
  memory->destroy(send_end_);
  memory->destroy(recv_cells_);
  memory->destroy(send_cells_);
  memory->destroy(recv_cells_self_);
  memory->destroy(send_cells_self_);
  memory->destroy(buf_recv_);
  memory->destroy(buf_send_);
  memory->destroy(buf_self_);
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::grow_recv(int n)
{
  memoryKK->destroy_kokkos(k_recv_cells, recv_cells);
  memoryKK->create_kokkos(k_recv_cells, recv_cells, n, "comm_grid:recv_cells");
  memoryKK->destroy_kokkos(k_buf_recv, buf_recv);
  memoryKK->create_kokkos(k_buf_recv, buf_recv, n * max_size, "comm_grid:buf_recv");
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::grow_send(int n)
{
  memoryKK->destroy_kokkos(k_send_cells, send_cells);
  memoryKK->create_kokkos(k_send_cells, send_cells, n, "comm_grid:send_cells");
  memoryKK->destroy_kokkos(k_buf_send, buf_send);
  memoryKK->create_kokkos(k_buf_send, buf_send, n * max_size, "comm_grid:buf_send");
}

/* ---------------------------------------------------------------------- */

void CommGridKokkos::grow_self(int n)
{
  memoryKK->destroy_kokkos(k_recv_cells_self, recv_cells_self);
  memoryKK->create_kokkos(k_recv_cells_self, recv_cells_self, n, "comm_grid:recv_cells_self");
  memoryKK->destroy_kokkos(k_send_cells_self, send_cells_self);
  memoryKK->create_kokkos(k_send_cells_self, send_cells_self, n, "comm_grid:send_cells_self");
  memoryKK->destroy_kokkos(k_buf_self, buf_self);
  memoryKK->create_kokkos(k_buf_self, buf_self, n * max_size, "comm_grid:buf_self");
}
