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
#include "comm_grid.h"
#include "grid.h"
#include "grid_vec.h"
#include "comm.h"
#include "memory.h"
#include "intersect_list.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CommGrid::CommGrid(LAMMPS *lmp) : Pointers(lmp)
{
  size_forward = 0;
  size_exchange = 0;
  max_size = 0;

  nrecv = 0;
  nsend = 0;
  nrecvproc = 0;
  nsendproc = 0;
  recvproc = nullptr;
  sendproc = nullptr;
  recv_begin = nullptr;
  recv_end = nullptr;
  send_begin = nullptr;
  send_end = nullptr;
  recv_cells = nullptr;
  send_cells = nullptr;
  buf_recv = nullptr;
  buf_send = nullptr;

  nrecv_self = 0;
  nsend_self = 0;
  recv_cells_self = nullptr;
  send_cells_self = nullptr;
  buf_self = nullptr;
  
  requests = nullptr;
}

/* ---------------------------------------------------------------------- */

CommGrid::~CommGrid()
{
  memory->destroy(recvproc);
  memory->destroy(sendproc);
  memory->destroy(recv_begin);
  memory->destroy(recv_end);
  memory->destroy(send_begin);
  memory->destroy(send_end);
  memory->destroy(recv_cells);
  memory->destroy(send_cells);
  memory->destroy(buf_recv);
  memory->destroy(buf_send);
  memory->destroy(recv_cells_self);
  memory->destroy(send_cells_self);
  memory->destroy(buf_self);
  
  delete [] requests;
}

/* ---------------------------------------------------------------------- */

void CommGrid::init()
{
  size_forward = grid->gvec->size_forward;
  size_exchange = grid->gvec->size_exchange;
  max_size = MAX(size_forward, size_exchange);

  recvproc = memory->create(recvproc, comm->nprocs, "comm_grid:recvproc");
  sendproc = memory->create(sendproc, comm->nprocs, "comm_grid:sendproc");
  recv_begin = memory->create(recv_begin, comm->nprocs, "comm_grid:recv_begin");
  recv_end = memory->create(recv_end, comm->nprocs, "comm_grid:recv_end");
  send_begin = memory->create(send_begin, comm->nprocs, "comm_grid:send_begin");
  send_end = memory->create(send_end, comm->nprocs, "comm_grid:send_end");
}

/* ---------------------------------------------------------------------- */

void CommGrid::setup()
{
  nrecv = 0;
  nsend = 0;
  nrecvproc = 0;
  nsendproc = 0;
  nrecv_self = 0;
  nsend_self = 0;
  
  int boxlo[3*comm->nprocs];
  MPI_Allgather(grid->sublo, 3, MPI_INT, boxlo, 3, MPI_INT, world);

  int boxhi[3*comm->nprocs];
  MPI_Allgather(grid->subhi, 3, MPI_INT, boxhi, 3, MPI_INT, world);

  IntersectList recvlist(lmp, comm->nprocs);
  IntersectList sendlist(lmp, comm->nprocs);

  int lo[3];
  int hi[3];
  for (int p = 0; p < comm->nprocs; p++) {
    if (comm->me != p) {
      int n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
			0, -1, 0, 0, 0, lo, hi, false);
      if (n > 0) {
	recvproc[nrecvproc] = p;
	recv_begin[nrecvproc] = nrecv;
	nrecv += n;
	recv_end[nrecvproc++] = nrecv;
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, 0, 0, 0, lo, hi, false);
      if (n > 0) {
	sendproc[nsendproc] = p;
	send_begin[nsendproc] = nsend;
	nsend += n;
	send_end[nsendproc++] = nsend;
	sendlist.add(p, lo, hi, n);
      }
    }
    if (grid->periodic[0]) {
      int n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
			0, -1, -grid->box[0], 0, 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, grid->box[0], 0, 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    0, -1, grid->box[0], 0, 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, -grid->box[0], 0, 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
    }
    if (grid->periodic[1]) {
      int n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
			0, -1, 0, -grid->box[1], 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, 0, grid->box[1], 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    0, -1, 0, grid->box[1], 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, 0, -grid->box[1], 0, lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
    }
    if (grid->periodic[2]) {
      int n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
			0, -1, 0, 0, -grid->box[2], lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, 0, 0, grid->box[2], lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    0, -1, 0, 0, grid->box[2], lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nrecvproc > 0 && recvproc[nrecvproc-1] == p) {
	    recv_end[nrecvproc-1] += n;
	  } else {
	    recvproc[nrecvproc] = p;
	    recv_begin[nrecvproc] = nrecv;
	    recv_end[nrecvproc++] = nrecv + n;
	  }
	  nrecv += n;
	} else {
	  nrecv_self += n;
	}
	recvlist.add(p, lo, hi, n);
      }
      n = intersect(grid->sublo, grid->subhi, &boxlo[3*p], &boxhi[3*p],
		    -1, 0, 0, 0, -grid->box[2], lo, hi, false);
      if (n > 0) {
	if (comm->me != p) {
	  if (nsendproc > 0 && sendproc[nsendproc-1] == p) {
	    send_end[nsendproc-1] += n;
	  } else {
	    sendproc[nsendproc] = p;
	    send_begin[nsendproc] = nsend;
	    send_end[nsendproc++] = nsend + n;
	  }
	  nsend += n;
	} else {
	  nsend_self += n;
	}
	sendlist.add(p, lo, hi, n);
      }
    }
  }

  if (nrecv_self != nsend_self)
    error->all(FLERR, "Conflicting self send and recv sizes (this is possibly a bug)");

  grow_recv(nrecv);
  grow_send(nsend);
  grow_self(nrecv_self);
  
  int irecv = 0;
  int irecv_self = 0;
  for (int i = 0; i < recvlist.n; i++) {
    for (int z = recvlist.boxlo[3*i+2]; z < recvlist.boxhi[3*i+2]; z++) {
      for (int y = recvlist.boxlo[3*i+1]; y < recvlist.boxhi[3*i+1]; y++) {
	for (int x = recvlist.boxlo[3*i]; x < recvlist.boxhi[3*i]; x++) {
	  if (comm->me != recvlist.procs[i]) {
	    recv_cells[irecv++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	  } else {
	    recv_cells_self[irecv_self++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
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
	    send_cells[isend++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	  } else {
	    send_cells_self[isend_self++] = x - grid->sublo[0] +
	      (y - grid->sublo[1]) * grid->subbox[0] +
	      (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	  }
  	}
      }
    }
  }
  
  if (requests) delete [] requests;
  requests = new MPI_Request[nrecvproc];
}

/* ---------------------------------------------------------------------- */

void CommGrid::forward_comm()
{
  for (int p = 0; p < nrecvproc; p++) {
    MPI_Irecv(&buf_recv[recv_begin[p] * size_forward],
	      (recv_end[p] - recv_begin[p]) * size_forward,
	      MPI_DOUBLE, recvproc[p], 0, world, &requests[p]);
  }
  for (int p = 0; p < nsendproc; p++) {
    int n = grid->gvec->pack_comm(send_end[p] - send_begin[p],
				  &send_cells[send_begin[p]],
				  buf_send);
    MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[p], 0, world);
  }
  MPI_Waitall(nrecvproc, requests, MPI_STATUS_IGNORE);
  for (int p = 0; p < nrecvproc; p++) {
    grid->gvec->unpack_comm(recv_end[p] - recv_begin[p],
			    &recv_cells[recv_begin[p]],
			    &buf_recv[recv_begin[p] * size_forward]);
  }
  grid->gvec->pack_comm(nsend_self, send_cells_self, buf_self);
  grid->gvec->unpack_comm(nrecv_self, recv_cells_self, buf_self);
}

/* ---------------------------------------------------------------------- */

void CommGrid::migrate()
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

void CommGrid::grow_recv(int n)
{
  memory->destroy(recv_cells);
  recv_cells = memory->create(recv_cells, n, "comm_grid:recv_cells");
  memory->destroy(buf_recv);
  buf_recv = memory->create(buf_recv, n * max_size, "comm_grid:buf_recv");
}

/* ---------------------------------------------------------------------- */

void CommGrid::grow_send(int n)
{
  memory->destroy(send_cells);
  send_cells = memory->create(send_cells, n, "comm_grid:send_cells");
  memory->destroy(buf_send);
  buf_send = memory->create(buf_send, n * max_size, "comm_grid:buf_send");
}

/* ---------------------------------------------------------------------- */

void CommGrid::grow_self(int n)
{
  memory->destroy(recv_cells_self);
  recv_cells_self = memory->create(recv_cells_self, n, "comm_grid:recv_cells_self");
  memory->destroy(send_cells_self);
  send_cells_self = memory->create(send_cells_self, n, "comm_grid:send_cells_self");
  memory->destroy(buf_self);
  buf_self = memory->create(buf_self, n * max_size, "comm_grid:buf_self");
}

/* ---------------------------------------------------------------------- */

int CommGrid::intersect(int *lo1, int *hi1,
			int *lo2, int *hi2,
			int ext1, int ext2,
			int xshift, int yshift, int zshift,
			int *intlo, int *inthi, bool corner)
{
  int tmplo1[3];
  memcpy(tmplo1, lo1, 3*sizeof(int));
  int tmphi1[3];
  memcpy(tmphi1, hi1, 3*sizeof(int));
  int tmplo2[3];
  memcpy(tmplo2, lo2, 3*sizeof(int));
  int tmphi2[3];
  memcpy(tmphi2, hi2, 3*sizeof(int));

  tmplo2[0] += xshift;
  tmplo2[1] += yshift;
  tmplo2[2] += zshift;

  tmphi2[0] += xshift;
  tmphi2[1] += yshift;
  tmphi2[2] += zshift;

  int size[3];
  for (int d = 0; d < 3; d++) {
    tmplo1[d] -= ext1;
    tmphi1[d] += ext1;
    tmplo2[d] -= ext2;
    tmphi2[d] += ext2;
    intlo[d] = MAX(tmplo1[d], tmplo2[d]);
    inthi[d] = MIN(tmphi1[d], tmphi2[d]);
    size[d] = inthi[d] - intlo[d];
  }
  if (size[0] <= 0 || size[1] <= 0 || size[2] <= 0)
    return 0;
  return size[0] * size[1] * size[2];
}
