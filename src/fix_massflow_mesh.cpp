/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "fix_mesh_surface.h"
#include "fix_neighlist_mesh.h"
#include "fix_property_atom.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_massflow_mesh.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassflowMesh::FixMassflowMesh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fix_mesh_(0),
  fix_counter_(0),
  once_(false),
  mass_(0.),
  nparticles_(0),
  nparticles_last_(0.),
  mass_last_(0.),
  t_count_(0.),
  delta_t_(0.),
  fp_(0),
  reset_t_count_(true),
  havePointAtOutlet_(false),
  insideOut_(false),
  screenflag_(false),
  delete_atoms_(false),
  mass_deleted_(0.),
  nparticles_deleted_(0)
{
    vectorZeroize3D(nvec_);
    vectorZeroize3D(pref_);
    vectorZeroize3D(sidevec_);
    vectorZeroize3D(pointAtOutlet_);

    // parse args for this class

    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"vec_side") == 0) {
            if(narg < iarg+4)
                error->fix_error(FLERR,this,"not enough arguments for 'vec_side'");
            iarg++;
            sidevec_[0] = atof(arg[iarg++]);
            sidevec_[1] = atof(arg[iarg++]);
            sidevec_[2] = atof(arg[iarg++]);
            if(vectorMag3D(sidevec_) == 0.)
                error->fix_error(FLERR,this,"vec_side > 0 required");
            hasargs = true;
        } else if(strcmp(arg[iarg],"mesh") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'insert_stream'");
            iarg++;
            fix_mesh_ = static_cast<FixMeshSurface*>(modify->find_fix_id_style(arg[iarg++],"mesh/surface"));
            if(!fix_mesh_)
                error->fix_error(FLERR,this,"fix mesh ID does not exist");
            hasargs = true;
        } else if(strcmp(arg[iarg],"count") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'insert_stream'");
            iarg++;
            if(strcmp(arg[iarg],"once") == 0)
                once_ = true;
            else if(strcmp(arg[iarg],"multiple") == 0)
                once_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'once' or 'multiple' after 'count'");
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"point_at_outlet") == 0) {
            if(narg < iarg+4)
                error->fix_error(FLERR,this,"not enough arguments for 'point_at_outlet'");
            havePointAtOutlet_ = true;
            iarg++;
            pointAtOutlet_[0] = atof(arg[iarg++]);
            pointAtOutlet_[1] = atof(arg[iarg++]);
            pointAtOutlet_[2] = atof(arg[iarg++]);
            hasargs = true;
        } else if(strcmp(arg[iarg],"inside_out") == 0) {
            insideOut_ = true;
            iarg++;
            if(!havePointAtOutlet_)
                error->fix_error(FLERR,this,"the setting 'inside_out' has no meaning in case you do not use 'point_at_outlet'");
            hasargs = true;
        } else if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0)
          {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");

            char* filecurrent = new char[strlen(arg[iarg+1]) + 8];
            if (1 < comm->nprocs) //open a separate file for each processor
                 sprintf(filecurrent,"%s%s%d",arg[iarg+1],".",comm->me);
            else  //open one file for proc 0
                 sprintf(filecurrent,"%s",arg[iarg+1]);

            if (strcmp(arg[iarg],"file") == 0)
                fp_ = fopen(filecurrent,"w");
            else
                fp_ = fopen(filecurrent,"a");
            if (fp_ == NULL) {
               char str[128];
               sprintf(str,"Cannot open file %s",arg[iarg+1]);
                error->fix_error(FLERR,this,str);
            }
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"screen") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg+1],"yes") == 0) screenflag_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) screenflag_ = false;
            else error->all(FLERR,"Illegal fix print command");
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"delete_atoms") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg+1],"yes") == 0) delete_atoms_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) delete_atoms_ = false;
            else error->all(FLERR,"Illegal delete command");
            iarg += 2;
            hasargs = true;
        } else
            error->fix_error(FLERR,this,"unknown keyword");
    }

    if(fp_ && 1 < comm->nprocs && 0 == comm->me)
      fprintf(screen,"**FixMassflowMesh: > 1 process - "
                     " will write to multiple files\n");

    // error checks on necessary args

    if(!once_ && delete_atoms_)
           error->fix_error(FLERR,this,"using 'delete_atoms yes' requires 'count once'");
    if( !once_ && havePointAtOutlet_)
        error->fix_error(FLERR,this,"setting 'point_at_outlet' requires 'count once'");
    if( vectorMag3D(sidevec_)==0. && !havePointAtOutlet_)
        error->fix_error(FLERR,this,"expecting keyword 'vec_side'");
    if (!fix_mesh_)
        error->fix_error(FLERR,this,"expecting keyword 'mesh'");

    // get reference point on face
    // calculate normalvec

    fix_mesh_->triMesh()->node(0,0,pref_);
    fix_mesh_->triMesh()->surfaceNorm(0,nvec_);
    double dot = vectorDot3D(nvec_,sidevec_);

    if(fabs(dot) < 1e-6 && !havePointAtOutlet_ )
        error->fix_error(FLERR,this,"need to change 'vec_side', it is currently in or to close to the mesh plane");
    else if(dot < 0.)
        vectorScalarMult3D(nvec_,-1.);

    restart_global = 1;

    vector_flag = 1;
    size_vector = 4;
    global_freq = 1; 
}

/* ---------------------------------------------------------------------- */

FixMassflowMesh::~FixMassflowMesh()
{
   if(fp_) fclose(fp_);
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::post_create()
{
    // add per-particle count flag

    char **fixarg;
    fixarg=new char*[9];
    for (int kk=0;kk<9;kk++) fixarg[kk]=new char[50+strlen(id)];

    sprintf(fixarg[0],"massflow_%s",id);
    sprintf(fixid_,    "massflow_%s",id);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    sprintf(fixarg[3],"massflow_%s",id);
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="no";    
    fixarg[7]="no";    
    fixarg[8]="-1.";     
    modify->add_fix(9,fixarg);
    delete []fixarg;

    fix_counter_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid_,"property/atom","scalar",0,0,style));

    // add neighbor list

    fix_neighlist_ = fix_mesh_->createOtherNeighList(igroup,id);
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::pre_delete(bool unfixflag)
{
    if (unfixflag) modify->delete_fix(fixid_);
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixMassflowMesh::init()
{
    if (atom->rmass_flag == 0)
        error->fix_error(FLERR,this,"requires atoms have mass");

    if(modify->n_fixes_style("multisphere"))
        error->fix_error(FLERR,this,"does not support multi-sphere");

    if(delete_atoms_ && 0 == atom->map_style)
        error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::setup(int vflag)
{
    // check if face planar
    
    if(!fix_mesh_->triMesh()->isPlanar() && !havePointAtOutlet_)
       error->fix_error(FLERR,this,"requires a planar face mass flow measurement or using 'point_at_outlet'");
}

/* ---------------------------------------------------------------------- */

int FixMassflowMesh::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    if(delete_atoms_) mask |= PRE_EXCHANGE;
    return mask;
}

/* ----------------------------------------------------------------------
   evaluate mass which went thru face
   all nearby particles
  -1 = default value for counter, set if particle is not
       a neighbor of the TriMesh
   0 = was not on nvec_ last step and is thus eligible
   1 = was on nvec_ side last step
   2 = do not re-count, was counted already
------------------------------------------------------------------------- */

void FixMassflowMesh::post_integrate()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *counter = fix_counter_->vector_atom;
    double dot,delta[3];
    int *neighs, *nneighs;
    double mass_this = 0.;
    int nparticles_this = 0.;
    double deltan;
    int *tag = atom->tag;

    TriMesh *mesh = fix_mesh_->triMesh();
    int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

    fix_neighlist_->getPointers(neighs,nneighs);

    // update time for counter
    // also store values for last invokation
    t_count_ += update->dt;
    if(!reset_t_count_)
    {
        nparticles_last_ = nparticles_;
        mass_last_ = mass_;
        reset_t_count_ = true;
    }

    // loop owned and ghost triangles
    // count only if owned particle

    for(int iTri = 0; iTri < nTriAll; iTri++)
    {
        for(int iNeigh = 0; iNeigh < nneighs[iTri]; iNeigh++,neighs++)
        {
            int iPart = *neighs;

            // skip ghost particles
            if(iPart >= nlocal) continue;

            // in case of once_ == true, ignore
            // everything which has been already counted
            if(compDouble(counter[iPart],2.)) continue;

            if(havePointAtOutlet_)
            {
                //get the vector from the particle center
                //to the next triangle
                deltan = fix_mesh_->triMesh()->resolveTriSphereContact(iPart,iTri,radius[iPart],x[iPart],delta);
                if(deltan < radius[iPart])
                {
                    vectorSubtract3D(x[iPart],pointAtOutlet_,nvec_); //vector pointing to the particle location
                    dot = vectorDot3D(delta,nvec_);
                }
                else //particle is not overlapping with mesh, so continue
                    continue;

                if(insideOut_) dot =-dot;
            }
            else
            {
                vectorSubtract3D(x[iPart],pref_,delta);
                dot = vectorDot3D(delta,nvec_);
            }

            // first-time, just set 0 or 1
            if(compDouble(counter[iPart],-1.))
            {
                counter[iPart] = dot <= 0. ? 0. : 1.;
                continue;
            }

            // not first time
            if(dot > 0.)
            {
                if(compDouble(counter[iPart],0.))
                {
                    mass_this += rmass[iPart];
                    nparticles_this ++;
                    if(delete_atoms_)
                        atom_tags_delete_.push_back(atom->tag[iPart]);

                    if (screenflag_ && screen)
                        fprintf(screen," %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g \n ",
                                       tag[iPart],2.*radius[iPart]/force->cg(),
                                       x[iPart][0],x[iPart][1],x[iPart][2],
                                       v[iPart][0],v[iPart][1],v[iPart][2]);
                    if (fp_)
                    {
                        fprintf(fp_," %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g \n ",
                                   tag[iPart],2.*radius[iPart]/force->cg(),
                                   x[iPart][0],x[iPart][1],x[iPart][2],
                                   v[iPart][0],v[iPart][1],v[iPart][2]);
                        fflush(fp_);
                    }
                }

                counter[iPart] = once_ ? 2. : 1.;
            }
            else // dot <= 0
            {
                counter[iPart] = 0.;
            }

        }
    }

    MPI_Sum_Scalar(mass_this,world);
    MPI_Sum_Scalar(nparticles_this,world);

    mass_ += mass_this;
    nparticles_ += nparticles_this;

}

/* ----------------------------------------------------------------------
   perform particle deletion of marked particles
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixMassflowMesh::pre_exchange()
{
    int nlocal = atom->nlocal;
    int iPart;
    double *counter = fix_counter_->vector_atom;

    if (update->ntimestep != next_reneighbor) return;

    if (delete_atoms_)
    {
        double mass_deleted_this_ = 0.;
        int nparticles_deleted_this_ = 0.;

        // delete particles

        while (atom_tags_delete_.size() > 0)
        {
            iPart = atom->map(atom_tags_delete_[0]);

            mass_deleted_this_ += atom->rmass[iPart];
            nparticles_deleted_this_++;

            atom->avec->copy(atom->nlocal-1,iPart,1);
            atom->nlocal--;

            atom_tags_delete_.erase(atom_tags_delete_.begin());
        }

        atom_tags_delete_.clear();

        MPI_Sum_Scalar(mass_deleted_this_,world);
        MPI_Sum_Scalar(nparticles_deleted_this_,world);

        mass_deleted_ += mass_deleted_this_;
        nparticles_deleted_ += nparticles_deleted_this_;

        if(nparticles_deleted_this_)
        {
            int i;
            if (atom->molecular == 0) {
              int *tag = atom->tag;
              for (i = 0; i < atom->nlocal; i++) tag[i] = 0;
              atom->tag_extend();
            }

            if (atom->tag_enable) {
              if (atom->map_style) {
                atom->nghost = 0;
                atom->map_init();
                atom->map_set();
              }
            }
        }
    }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMassflowMesh::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = mass_;
  list[n++] = t_count_;
  list[n++] = mass_last_;
  list[n++] = static_cast<double>(nparticles_last_);
  list[n++] = mass_deleted_;
  list[n++] = static_cast<double>(nparticles_deleted_);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMassflowMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  mass_ = list[n++];
  t_count_ = list[n++];
  mass_last_ = list[n++];
  nparticles_last_ = static_cast<int>(list[n++]);
  mass_deleted_ = list[n++];
  nparticles_deleted_ = static_cast<int>(list[n++]);
}

/* ----------------------------------------------------------------------
   output mass
------------------------------------------------------------------------- */

double FixMassflowMesh::compute_vector(int index)
{
    if(reset_t_count_)
    {
        delta_t_ = t_count_;
        t_count_ = 0.;
        reset_t_count_ = false;
    }

    if(index == 0)
        return mass_;
    if(index == 1)
        return static_cast<double>(nparticles_);
    if(index == 2)
        return delta_t_ == 0. ? 0. : (mass_-mass_last_)/delta_t_;
    if(index == 3)
        return delta_t_ == 0. ? 0. : static_cast<double>(nparticles_-nparticles_last_)/delta_t_;
    if(index == 4)
        return mass_deleted_;
    if(index == 5)
        return static_cast<double>(nparticles_deleted_);

    return 0.;
}
