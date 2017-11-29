/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "mpi_liggghts.h"
#include "fix_ave_euler.h"
#include "fix_multisphere.h"
#include "compute_stress_atom.h"
#include "math_extra_liggghts.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "neighbor.h"
#include "region.h"
#include "update.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#define BIG 1000000000

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAveEuler::FixAveEuler(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  parallel_(true),
  exec_every_(1),
  box_change_size_(false),
  box_change_domain_(false),
  cell_size_ideal_rel_(3.),
  cell_size_ideal_(0.),
  ncells_(0),
  ncells_max_(0),
  ncellptr_max_(0),
  cellhead_(NULL),
  cellptr_(NULL),
  idregion_(NULL),
  region_(NULL),
  center_(NULL),
  v_av_(NULL),
  vol_fr_(NULL),
  weight_(NULL),
  radius_(NULL),
  ncount_(NULL),
  mass_(NULL),
  stress_(NULL),
  compute_stress_(NULL),
  random_(0)
{
  // this fix produces a global array

  array_flag = 1;
  size_array_rows = BIG;
  size_array_cols = 3 + 1 + 3;

  triclinic_ = domain->triclinic;  

  // random number generator, seed is hardcoded
  random_ = new RanPark(lmp,"15485863");

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
  if(cell_size_ideal_rel_ < 1.)
    error->fix_error(FLERR,this,"'cell_size_relative' > 1 required");

  if(strcmp(arg[iarg++],"parallel"))
    error->fix_error(FLERR,this,"expecting keyword 'parallel'");
  if(strcmp(arg[iarg],"yes") == 0)
    parallel_ = true;
  else if(strcmp(arg[iarg],"no") == 0)
    parallel_ = false;
  else
    error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'parallel'");
  iarg++;

  while(iarg < narg)
  {
     if (strcmp(arg[iarg],"basevolume_region") == 0) {
       if (iarg+2 > narg) error->fix_error(FLERR,this,"");
       int iregion = domain->find_region(arg[iarg+1]);
       if (iregion == -1)
         error->fix_error(FLERR,this,"region ID does not exist");
       int n = strlen(arg[iarg+1]) + 1;
       idregion_ = new char[n];
       strcpy(idregion_,arg[iarg+1]);
       region_ = domain->regions[iregion];
       iarg += 2;
     } else {
       char *errmsg = new char[strlen(arg[iarg])+50];
       sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg]);
       error->fix_error(FLERR,this,errmsg);
       delete []errmsg;
     }
  }
}

/* ---------------------------------------------------------------------- */

FixAveEuler::~FixAveEuler()
{
  memory->destroy(cellhead_);
  memory->destroy(cellptr_);
  if(idregion_) delete []idregion_;
  memory->destroy(center_);
  memory->destroy(v_av_);
  memory->destroy(vol_fr_);
  memory->destroy(weight_);
  memory->destroy(radius_);
  memory->destroy(ncount_);
  memory->destroy(mass_);
  memory->destroy(stress_);
  if (random_) delete random_;
}

/* ---------------------------------------------------------------------- */

void FixAveEuler::post_create()
{
  //  stress computation, just for pairwise contribution
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

  // ensure that the compute is calculated in the first time step
  bigint nfirst = (update->ntimestep/nevery)*nevery + nevery;
  compute_stress_->addstep(nfirst);
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

  // does not work with MS
  FixMultisphere* fix_ms = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
  if(fix_ms)
      error->fix_error(FLERR,this,"does not work with multisphere");

  // check if box constant
  box_change_size_ = false;
  if(domain->box_change_size)
    box_change_size_ = true;
  box_change_domain_ = false;
  if(domain->box_change_domain)
    box_change_domain_ = true;

  if (region_)
  {
    int iregion = domain->find_region(idregion_);
    if (iregion == -1)
     error->fix_error(FLERR,this,"regions used my this command must not be deleted");
    region_ = domain->regions[iregion];
  }

  // error checks

  if (!parallel_ && 1 == domain->triclinic)
    error->fix_error(FLERR,this,"triclinic boxes only support 'parallel=yes'");
}

/* ----------------------------------------------------------------------
   setup initial bins
   only does averaging if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveEuler::setup(int vflag)
{
    setup_bins();
    end_of_step();
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
        lo_[dim] = parallel_ ? domain->sublo[dim] : domain->boxlo[dim];
        hi_[dim] = parallel_ ? domain->subhi[dim] : domain->boxhi[dim];
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

    // round down (makes cell size larger)
    // at least one cell
    for(int dim = 0; dim < 3; dim++)
    {
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
        memory->grow(weight_,ncells_max_,  "ave/euler:weight_");
        memory->grow(radius_,ncells_max_,  "ave/euler:radius_");
        memory->grow(ncount_,ncells_max_,    "ave/euler:ncount_");
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

    // calculate weight_[icell]
    if(!region_)
    {
        for(int icell = 0; icell < ncells_max_; icell++)
            weight_[icell] = 1.;
    }

    // MC calculation if region_ exists
    if(region_)
    {
        double x_try[3];
        int ibin;
        int ntry = ncells_ * ntry_per_cell(); // number of MC tries
        double contribution = 1./static_cast<double>(ntry_per_cell());  // contrib of each try

        for(int icell = 0; icell < ncells_max_; icell++)
            weight_[icell] = 0.;

        for(int itry = 0; itry < ntry; itry++)
        {
            x_try[0] = lo_[0]+(hi_[0]-lo_[0])*random_->uniform();
            x_try[1] = lo_[1]+(hi_[1]-lo_[1])*random_->uniform();
            x_try[2] = lo_[2]+(hi_[2]-lo_[2])*random_->uniform();
            if(region_->match(x_try[0],x_try[1],x_try[2]))
            {
                ibin = coord2bin(x_try);
                // only do this for points in my subbox
                if(ibin >= 0)
                    weight_[ibin] += contribution;
            }
        }

        // allreduce weights
        MPI_Sum_Vector(weight_,ncells_,world);

        // limit weight to 1
        for(int icell = 0; icell < ncells_max_; icell++)
            if(weight_[icell] > 1.) weight_[icell] = 1.;
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
    
    // have to adapt grid if box size changes
    if(box_change_size_ || (parallel_ && box_change_domain_))
    {
        if(region_)
            error->warning(FLERR,"Fix ave/euler using 'basevolume_region'"
                                "and changing simulation or load-balancing is huge over-head");
        setup_bins();
    }

    // bin atoms
    bin_atoms();

    // calculate Eulerian grid properties
    // performs allreduce if necessary
    calculate_eu();
}

/* ---------------------------------------------------------------------- */

int FixAveEuler::ncells_pack()
{
    // in parallel mode, each proc will pack
    if (parallel_)
        return ncells_;

    // in serial mode, only proc 0 will pack
    if(0 == comm->me)
        return ncells_;
    else
        return 0;
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
  // also use ghost atoms
  // skip if any atom is out of (sub) box

  for (i = nall-1; i >= 0; i--)
  {
      if(! (mask[i] & groupbit)) continue;

      ibin = coord2bin(x[i]);

      // particles outside grid may return values ibin < 0 || ibin >= ncells_
      // these are ignores
      
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

      // skip particles outside my subdomain
      
      if(x[i] <= domain->sublo[i] || x[i] >= domain->subhi[i])
        return -1;
      float_iCell[i] = (x[i]-lo_[i])*cell_size_inv_[i];
      iCell[i] = static_cast<int> (float_iCell[i]);
    }
  }

  return iCell[2]*ncells_dim_[1]*ncells_dim_[0] + iCell[1]*ncells_dim_[0] + iCell[0];
}

/* ----------------------------------------------------------------------
   calculate Eulerian data, use interpolation function
------------------------------------------------------------------------- */

void FixAveEuler::calculate_eu()
{
    //int ncount;
    double * const * const v = atom->v;
    double * const radius = atom->radius;
    double * const rmass = atom->rmass;

    double prefactor_vol_fr = 4./3.*M_PI/cell_volume_;
    double prefactor_stress = 1./cell_volume_;
    double vel_x_mass[3];
    #ifdef SUPERQUADRIC_ACTIVE_FLAG
    const double * const volume = atom->volume;
    const int superquadric_flag = atom->superquadric_flag;
    #endif

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
    // optionally plus its 26 neighs

    for(int icell = 0; icell < ncells_; icell++)
    {
        ncount_[icell] = 0;
        vectorZeroize3D(v_av_[icell]);
        vol_fr_[icell] = 0.;
        radius_[icell] = 0.;
        mass_[icell] = 0.;
        vectorZeroizeN(stress_[icell],7);

        // skip if no particles in cell
        
        if(-1 == cellhead_[icell])
            continue;

        // add contributions of particle - v and volume fraction
        // v is favre-averaged (mass-averaged)
        // radius is number-averaged

        for(int iatom = cellhead_[icell/*+stencil*/]; iatom >= 0; iatom = cellptr_[iatom])
        {
            vectorScalarMult3D(v[iatom],rmass[iatom],vel_x_mass);
            vectorAdd3D(v_av_[icell],vel_x_mass,v_av_[icell]);
            double r = radius[iatom];
            #ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(superquadric_flag)
                r = cbrt(0.75 * volume[iatom] / M_PI);
            #endif
            vol_fr_[icell] += r*r*r;
            radius_[icell] += r;
            mass_[icell] += rmass[iatom];
            ncount_[icell]++;
        }

    }

    // allreduce contributions so far if not parallel
    if(!parallel_ && ncells_ > 0)
    {
        MPI_Sum_Vector(&(v_av_[0][0]),3*ncells_,world);
        MPI_Sum_Vector(vol_fr_,ncells_,world);
        MPI_Sum_Vector(radius_,ncells_,world);
        MPI_Sum_Vector(mass_,ncells_,world);
        MPI_Sum_Vector(ncount_,ncells_,world);
    }

    // perform further calculations

    double eps_ntry = 1./static_cast<double>(ntry_per_cell());
    for(int icell = 0; icell < ncells_; icell++)
    {
        // calculate average vel and radius
        if(ncount_[icell]) vectorScalarDiv3D(v_av_[icell],mass_[icell]);
        if(ncount_[icell]) radius_[icell]/=static_cast<double>(ncount_[icell]);

        // calculate volume fraction
        //safety check, add an epsilon to weight if any particle ended up in that cell
        if(vol_fr_[icell] > 0. && MathExtraLiggghts::compDouble(weight_[icell],0.,1e-6))
           weight_[icell] = eps_ntry;
        if(weight_[icell] < eps_ntry )
            vol_fr_[icell] = 0.;
        else
            vol_fr_[icell] *= prefactor_vol_fr/weight_[icell];
        
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
        }
        stress_[icell][0] = -0.333333333333333*(stress_[icell][1]+stress_[icell][2]+stress_[icell][3]);
        if(weight_[icell] < eps_ntry)
            vectorZeroizeN(stress_[icell],7);
        else
            vectorScalarMultN(7,stress_[icell],prefactor_stress/weight_[icell]);
    }

    // allreduce stress if not parallel
    if(!parallel_ && ncells_ > 0)
    {
        MPI_Sum_Vector(&(stress_[0][0]),7*ncells_,world);

        // recalc pressure based on allreduced stress
        for(int icell = 0; icell < ncells_; icell++)
            stress_[icell][0] = -0.333333333333333*(stress_[icell][1]+stress_[icell][2]+stress_[icell][3]);
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
  else if(j < 14) return stress_[i][j-7];
  else if(j < 15) return radius_[i];
  else return 0.0;
}
