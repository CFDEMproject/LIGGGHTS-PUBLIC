/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author for triclinic: Andreas Aigner (JKU)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "fix_ave_euler.h"
#include "compute_stress_atom.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#define BIG 1000000000

#define INVOKED_PERATOM 8

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAveEuler::FixAveEuler(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  exec_every_(1),
  box_change_(false),
  cell_size_ideal_rel_(3.),
  cell_size_ideal_(0.),
  ncells_(0),
  ncells_max_(0),
  ncellptr_max_(0),
  cellhead_(NULL),
  cellptr_(NULL),
  center_(NULL),
  v_av_(NULL),
  vol_fr_(NULL),
  radius_(NULL),
  mass_(NULL),
  stress_(NULL),
  compute_stress_(NULL)
{
  // this fix produces a global array

  array_flag = 1;
  size_array_rows = BIG;
  size_array_cols = 3 + 1 + 3;

  triclinic_ = domain->triclinic;  

  // parse args
  if (narg < 6) error->all(FLERR,"Illegal fix ave/pic command");
  int iarg = 3;

  if(strcmp(arg[iarg++],"nevery"))
    error->fix_error(FLERR,this,"expecting keyword 'nevery'");
  exec_every_ = force->inumeric(FLERR,arg[iarg++]);
  if(exec_every_ < 1)
    error->fix_error(FLERR,this,"'nevery' > 0 required");
  nevery = exec_every_;

  if(strcmp(arg[iarg++],"cell_size_relative"))
    error->fix_error(FLERR,this,"expecting keyword 'cell_size_relative'");
  cell_size_ideal_rel_ = force->numeric(FLERR,arg[iarg++]);
  if(cell_size_ideal_rel_ < 3.)
    error->fix_error(FLERR,this,"'cell_size_relative' > 3 required");

  // add nfirst to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  bigint nfirst = (update->ntimestep/nevery)*nevery + nevery;
  modify->addstep_compute_all(nfirst);
}

/* ---------------------------------------------------------------------- */

FixAveEuler::~FixAveEuler()
{
  memory->destroy(cellhead_);
  memory->destroy(cellptr_);
  memory->destroy(center_);
  memory->destroy(v_av_);
  memory->destroy(vol_fr_);
  memory->destroy(radius_);
  memory->destroy(mass_);
  memory->destroy(stress_);
}

/* ---------------------------------------------------------------------- */

void FixAveEuler::post_create()
{
  //  register per particle properties
  if(!compute_stress_)
  {
        const char* arg[4];
        arg[0]="stress_faveu";
        arg[1]="all";
        arg[2]="stress/atom";
        arg[3]="pair";

        modify->add_compute(4,(char**)arg);
        compute_stress_ = static_cast<ComputeStressAtom*>(modify->compute[modify->find_compute(arg[0])]);
  }
}

/* ---------------------------------------------------------------------- */

int FixAveEuler::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveEuler::init()
{
  if(!atom->radius_flag)
    error->fix_error(FLERR,this,"requires atom attribute radius");
  if(!atom->rmass_flag)
    error->fix_error(FLERR,this,"requires atom attribute mass");

  // check if box constant
  box_change_ = false;
  if(domain->box_change)
    box_change_ = true;
}

/* ----------------------------------------------------------------------
   setup initial bins
   only does averaging if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveEuler::setup(int vflag)
{
    setup_bins();
    //end_of_step();
}

/* ----------------------------------------------------------------------
   setup 3d bins and their extent and coordinates
   called at setup() and when averaging occurs if box size changes
   similar to FixAveSpatial::setup_bins() and PairDSMC::init_style()

   bins are subbox - skin/2 so owned particles cannot move outside
   bins - so do not have to extrapolate
------------------------------------------------------------------------- */

void FixAveEuler::setup_bins()
{
    int ibin;

    // calc ideal cell size as multiple of max cutoff
    cell_size_ideal_ = cell_size_ideal_rel_ * (neighbor->cutneighmax-neighbor->skin);
    
    for(int dim = 0; dim < 3; dim++)
    {
        // get bounds
      if (triclinic_ == 0) {
        lo_[dim] = domain->sublo[dim];
        hi_[dim] = domain->subhi[dim];
      } else {
        lo_lamda_[dim] = domain->sublo_lamda[dim];
        hi_lamda_[dim] = domain->subhi_lamda[dim];

        cell_size_ideal_lamda_[dim] = cell_size_ideal_/domain->h[dim];
      }
    }
    if (triclinic_) {
      domain->lamda2x(lo_lamda_,lo_);
      domain->lamda2x(hi_lamda_,hi_);
    }

    for(int dim = 0; dim < 3; dim++) {
      // round down (makes cell size larger)
      // at least one cell
      if (triclinic_) {
      ncells_dim_[dim] = static_cast<int>((hi_lamda_[dim]-lo_lamda_[dim])/cell_size_ideal_lamda_[dim]);
      if (ncells_dim_[dim] < 1) {
        ncells_dim_[dim] = 1;
        error->warning(FLERR,"Number of cells for fix_ave_euler was less than 1");
      }
        cell_size_lamda_[dim] = (hi_lamda_[dim]-lo_lamda_[dim])/static_cast<double>(ncells_dim_[dim]);

        cell_size_[dim] = cell_size_lamda_[dim]*domain->h[dim];
      } else {
        ncells_dim_[dim] = static_cast<int>((hi_[dim]-lo_[dim])/cell_size_ideal_);
      if (ncells_dim_[dim] < 1) {
        ncells_dim_[dim] = 1;
        
        error->warning(FLERR,"Number of cells for fix_ave_euler was less than 1");
      }
        cell_size_[dim] = (hi_[dim]-lo_[dim])/static_cast<double>(ncells_dim_[dim]);
      }
    }

    for(int dim = 0; dim < 3; dim++)
    {
        cell_size_inv_[dim] = 1./cell_size_[dim];
      if (triclinic_) cell_size_lamda_inv_[dim] = 1./cell_size_lamda_[dim];

        // add 2 extra cells for ghosts
        // do not need to match subbox expaned by cutneighmax
        // cell size will always be larger than cutneighmax
        //ncells_dim_[dim] += 2;

        // do not change cell_size_ but adapt lo_/hi_
        //lo_[dim] -= cell_size_[dim];
        //hi_[dim] += cell_size_[dim];
    }

    ncells_ = ncells_dim_[0]*ncells_dim_[1]*ncells_dim_[2];
    
    cell_volume_ = cell_size_[0]*cell_size_[1]*cell_size_[2];

    // (re) allocate spatial bin arrays
    if (ncells_ > ncells_max_)
    {
        ncells_max_ = ncells_;
        memory->grow(cellhead_,ncells_max_,"ave/euler:cellhead_");
        memory->grow(center_,ncells_max_,3,"ave/euler:center_");
        memory->grow(v_av_,  ncells_max_,3,"ave/euler:v_av_");
        memory->grow(vol_fr_,ncells_max_,  "ave/euler:vol_fr_");
        memory->grow(radius_,ncells_max_,  "ave/euler:radius_");
        memory->grow(mass_,ncells_max_,    "ave/euler:mass_");
        memory->grow(stress_,ncells_max_,7,"ave/euler:stress_");
    }

    // calculate center corrdinates for cells
    for(int ix = 0; ix < ncells_dim_[0]; ix++)
    {
        for(int iy = 0; iy < ncells_dim_[1]; iy++)
        {
            for(int iz = 0; iz < ncells_dim_[2]; iz++)
            {
                ibin = iz*ncells_dim_[1]*ncells_dim_[0] + iy*ncells_dim_[0] + ix;

                if (triclinic_) {
                  center_[ibin][0] = lo_lamda_[0] + (static_cast<double>(ix)+0.5) * cell_size_lamda_[0];
                  center_[ibin][1] = lo_lamda_[1] + (static_cast<double>(iy)+0.5) * cell_size_lamda_[1];
                  center_[ibin][2] = lo_lamda_[2] + (static_cast<double>(iz)+0.5) * cell_size_lamda_[2];
                  domain->lamda2x(center_[ibin],center_[ibin]);

                } else {
                center_[ibin][0] = lo_[0] + (static_cast<double>(ix)+0.5) * cell_size_[0];
                center_[ibin][1] = lo_[1] + (static_cast<double>(iy)+0.5) * cell_size_[1];
                center_[ibin][2] = lo_[2] + (static_cast<double>(iz)+0.5) * cell_size_[2];
                }
            }
        }
    }

    // print to screen and log
    /*
    if (comm->me == 0)
    {
        if (screen) fprintf(screen,"Euler cell size on proc %d = %f (%d x %d x %d grid)\n",
            comm->me,cell_size_ideal_,ncells_dim_[0],ncells_dim_[1],ncells_dim_[2]);
        if (logfile) fprintf(logfile,"Euler cell size on proc %d = %f (%d x %d x %d grid)\n",
            comm->me,cell_size_ideal_,ncells_dim_[0],ncells_dim_[1],ncells_dim_[2]);
    }
    */
}

/* ---------------------------------------------------------------------- */

void FixAveEuler::end_of_step()
{
    
    // have to adapt grid if box changes
    if(box_change_)
        setup_bins();

    // bin atoms
    bin_atoms();

    // calculate Eulerian grid properties
    calculate_eu();
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
   this also implies we do not need to wrap around PBCs
   bin ghost atoms only if inside my grid
------------------------------------------------------------------------- */

void FixAveEuler::bin_atoms()
{
  int i,ibin;
  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < ncells_max_; i++)
    cellhead_[i] = -1;

  // re-alloc cellptr_ if necessary
  if(nall > ncellptr_max_)
  {
      ncellptr_max_ = nall;
      memory->grow(cellptr_,ncellptr_max_,"ave/pic:cellptr_");
  }

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list

  for (i = nall-1; i >= 0; i--)
  {
      if(! (mask[i] & groupbit)) continue;

      ibin = coord2bin(x[i]);

      // ghosts outside grid may return values ibin < 0 || ibin >= ncells_
      // lets ignore them
      if (ibin < 0 || ibin >= ncells_) {
        
        continue;
      }

      cellptr_[i] = cellhead_[ibin];
      cellhead_[ibin] = i;
  }
}

/* ----------------------------------------------------------------------
   map coord to grid, also return ix,iy,iz indices in each dim
------------------------------------------------------------------------- */

inline int FixAveEuler::coord2bin(double *x)
{
  int i,iCell[3];
  double float_iCell[3];

  if (triclinic_) {
    double tmp_x[3];
    domain->x2lamda(x,tmp_x);
    for (i=0;i<3;i++) {
      float_iCell[i] = (tmp_x[i]-lo_lamda_[i])*cell_size_lamda_inv_[i];
      iCell[i] = static_cast<int> (float_iCell[i] >= 0 ? float_iCell[i] : float_iCell[i]-1);
    }
  } else {
    for (i=0;i<3;i++) {
      float_iCell[i] = (x[i]-lo_[i])*cell_size_inv_[i];
      iCell[i] = static_cast<int> (float_iCell[i] >= 0 ? float_iCell[i] : float_iCell[i]-1);
    }
  }

  return iCell[2]*ncells_dim_[1]*ncells_dim_[0] + iCell[1]*ncells_dim_[0] + iCell[0];
}

/* ----------------------------------------------------------------------
   calculate Eulerian data, use interpolation function
------------------------------------------------------------------------- */

void FixAveEuler::calculate_eu()
{
    int ncount;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;

    double prefactor_vol_fr = 4./3.*M_PI/cell_volume_;
    double prefactor_stress = 1./cell_volume_;
    double vel_x_mass[3];

    // wrap compute with clear/add
    modify->clearstep_compute();

    // invoke compute if not previously invoked
    
    if (!(compute_stress_->invoked_flag & INVOKED_PERATOM)) {
        compute_stress_->compute_peratom();
        compute_stress_->invoked_flag |= INVOKED_PERATOM;
    }

    // forward comm per-particle stress from compute so neighs have it
    comm->forward_comm_compute(compute_stress_);

    // need to get pointer here since compute_peratom() may realloc
    double **stress_atom = compute_stress_->array_atom;

    // loop all binned particles
    // each particle can contribute to the cell that it has been binned
    // to plus its 26 neighs

    for(int icell = 0; icell < ncells_; icell++)
    {
        ncount = 0;
        vectorZeroize3D(v_av_[icell]);
        vol_fr_[icell] = 0.;
        radius_[icell] = 0.;
        mass_[icell] = 0.;
        vectorZeroizeN(stress_[icell],7);

        // add contributions of particle - v and volume fraction
        // v is favre-averaged (mass-averaged)
        // radius is number-averaged

        for(int iatom = cellhead_[icell/*+stencil*/]; iatom >= 0; iatom = cellptr_[iatom])
        {
            vectorScalarMult3D(v[iatom],rmass[iatom],vel_x_mass);
            vectorAdd3D(v_av_[icell],vel_x_mass,v_av_[icell]);
            vol_fr_[icell] += radius[iatom]*radius[iatom]*radius[iatom];
            radius_[icell] += radius[iatom];
            mass_[icell] += rmass[iatom];
            ncount++;
        }

        if(ncount) vectorScalarDiv3D(v_av_[icell],static_cast<double>(ncount)*mass_[icell]);
        if(ncount) radius_[icell]/=static_cast<double>(ncount);
        vol_fr_[icell] *= prefactor_vol_fr;

        // add contribution of particle - stress
        // need v before can calculate stress
        // stress is molecular diffusion + contact forces

        for(int iatom = cellhead_[icell/*+stencil*/]; iatom >= 0; iatom = cellptr_[iatom])
        {
            stress_[icell][1] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][0]-v_av_[icell][0]) + stress_atom[iatom][0];
            stress_[icell][2] += -rmass[iatom]*(v[iatom][1]-v_av_[icell][1])*(v[iatom][1]-v_av_[icell][1]) + stress_atom[iatom][1];
            stress_[icell][3] += -rmass[iatom]*(v[iatom][2]-v_av_[icell][2])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][2];
            stress_[icell][4] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][1]-v_av_[icell][1]) + stress_atom[iatom][3];
            stress_[icell][5] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][4];
            stress_[icell][6] += -rmass[iatom]*(v[iatom][1]-v_av_[icell][1])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][5]; 
            stress_[icell][0] = -0.333333333333333*(stress_[icell][1]+stress_[icell][2]+stress_[icell][3]);
        }
        vectorScalarMultN(7,stress_[icell],prefactor_stress);
    }

    // wrap with clear/add
    modify->addstep_compute(update->ntimestep + exec_every_);
}

/* ----------------------------------------------------------------------
   return I,J array value
   if I exceeds current bins, return 0.0 instead of generating an error
   column 1,2,3 = bin coords, next column = vol fr,
   remaining columns = vel, stress, radius
------------------------------------------------------------------------- */

double FixAveEuler::compute_array(int i, int j)
{
  if(i >= ncells_) return 0.0;

  else if(j < 3) return center_[i][j];
  else if(j == 3) return vol_fr_[i];
  else if(j < 7) return v_av_[i][j-4];
  else if(j == 7) return stress_[i][0];
  else if(j < 14) return stress_[i][j-8];
  else if(j < 15) return radius_[i];
  else return 0.0;
}
