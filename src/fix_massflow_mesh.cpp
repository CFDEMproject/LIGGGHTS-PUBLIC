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

#include <cmath>
#include <stdlib.h>
#include <string.h>
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
#include "fix_multisphere.h"
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
  delete_atoms_(false),
  mass_deleted_(0.),
  nparticles_deleted_(0),
  once_(false),
  fix_orientation_(0),
  fix_counter_(0),
  fix_mesh_(0),
  fix_neighlist_(0),
  fix_volumeweight_ms_(0),
  havePointAtOutlet_(false),
  insideOut_(false),
  mass_(0.),
  nparticles_(0.),
  fix_property_(0),
  property_sum_(0.),
  screenflag_(false),
  fp_(0),
  writeTime_(false),
  mass_last_(0.),
  nparticles_last_(0.),
  t_count_(0.),
  delta_t_(0.),
  reset_t_count_(true)
{
    vectorZeroize3D(nvec_);
    vectorZeroize3D(pref_);
    vectorZeroize3D(sidevec_);
    vectorZeroize3D(pointAtOutlet_);

    // parse args for this class

    iarg_ = 3;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"vec_side") == 0) {
            if(narg < iarg_+4)
                error->fix_error(FLERR,this,"not enough arguments for 'vec_side'");
            iarg_++;
            sidevec_[0] = atof(arg[iarg_++]);
            sidevec_[1] = atof(arg[iarg_++]);
            sidevec_[2] = atof(arg[iarg_++]);
            if(vectorMag3D(sidevec_) == 0.)
                error->fix_error(FLERR,this,"vec_side > 0 required");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"mesh") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'insert_stream'");
            iarg_++;
            fix_mesh_ = static_cast<FixMeshSurface*>(modify->find_fix_id_style(arg[iarg_++],"mesh/surface"));
            if(!fix_mesh_)
                error->fix_error(FLERR,this,"fix mesh ID does not exist");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"sum_property") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sum_property'");
            iarg_++;
            fix_property_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(arg[iarg_++],"property/atom","scalar",0,0,style));
            hasargs = true;
        } else if(strcmp(arg[iarg_],"count") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'insert_stream'");
            iarg_++;
            if(strcmp(arg[iarg_],"once") == 0)
                once_ = true;
            else if(strcmp(arg[iarg_],"multiple") == 0)
                once_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'once' or 'multiple' after 'count'");
            iarg_++;
            hasargs = true;
        } else if( strcmp(arg[iarg_],"writeTime") == 0) {
            writeTime_ = true;
            iarg_++;
            hasargs = true;
        } else if(strcmp(arg[iarg_],"point_at_outlet") == 0) {
            if(narg < iarg_+4)
                error->fix_error(FLERR,this,"not enough arguments for 'point_at_outlet'");
            havePointAtOutlet_ = true;
            iarg_++;
            pointAtOutlet_[0] = atof(arg[iarg_++]);
            pointAtOutlet_[1] = atof(arg[iarg_++]);
            pointAtOutlet_[2] = atof(arg[iarg_++]);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"inside_out") == 0) {
            insideOut_ = true;
            iarg_++;
            if(!havePointAtOutlet_)
                error->fix_error(FLERR,this,"the setting 'inside_out' has no meaning in case you do not use 'point_at_outlet'");
            hasargs = true;
        } else if (strcmp(arg[iarg_],"file") == 0 || strcmp(arg[iarg_],"append") == 0)
          {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");

            char* filecurrent = new char[strlen(arg[iarg_+1]) + 8];
            if (1 < comm->nprocs) //open a separate file for each processor
                 sprintf(filecurrent,"%s%s%d",arg[iarg_+1],".",comm->me);
            else  //open one file for proc 0
                 sprintf(filecurrent,"%s",arg[iarg_+1]);

            if (strcmp(arg[iarg_],"file") == 0)
                fp_ = fopen(filecurrent,"w");
            else
                fp_ = fopen(filecurrent,"a");
            delete[] filecurrent;
            if (fp_ == NULL) {
               char str[512];
               sprintf(str,"Cannot open file %s",arg[iarg_+1]);
                error->fix_error(FLERR,this,str);
            }
            iarg_ += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"screen") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg_+1],"yes") == 0) screenflag_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) screenflag_ = false;
            else error->all(FLERR,"Illegal fix print command");
            iarg_ += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"delete_atoms") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg_+1],"yes") == 0) delete_atoms_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) delete_atoms_ = false;
            else error->all(FLERR,"Illegal delete command");
            iarg_ += 2;
            hasargs = true;
        } else if(strcmp(style,"massflow/mesh") == 0)
            error->fix_error(FLERR,this,"unknown keyword");
    }

    if(fp_ && 1 < comm->nprocs && 0 == comm->me)
      fprintf(screen,"**FixMassflowMesh: > 1 process - "
                     " will write to multiple files\n");
    if(fp_)
    {
        //write header
        fprintf(fp_,"# ID");

        if(writeTime_)
          fprintf(fp_," time ");

        fprintf(fp_," diameter x y z u v w");

        fprintf(fp_,"  (ex ey ez, color)\n");

        fflush(fp_);
    }

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
    setRefPoint();

    restart_global = 1;

    vector_flag = 1;
    size_vector = 6;
    if(fix_property_)
        size_vector = 7;
    global_freq = 1; 
}

/* ---------------------------------------------------------------------- */

FixMassflowMesh::~FixMassflowMesh()
{
    if(fp_)
        fclose(fp_);
    // do not delete fix_neighlist_, this will be handled by fix mesh/surface itself
    fix_neighlist_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::post_create()
{
    // add per-particle count flag

    const char * fixarg[9];

    sprintf(fixid_,"massflow_%s",id);
    fixarg[0]=fixid_;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=fixid_;
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="no";    
    fixarg[7]="no";    
    fixarg[8]="-1.";     
    modify->add_fix(9,const_cast<char**>(fixarg));

    fix_counter_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid_,"property/atom","scalar",0,0,style));

    // add neighbor list

    fix_neighlist_ = fix_mesh_->createOtherNeighList(igroup,id);

    fix_volumeweight_ms_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("volumeweight_ms","property/atom","scalar",0,0,style,false));

    // need to find multisphere here to be able to add property
    
    FixMultisphere *fix_ms = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(fix_ms)
    {
        Multisphere *ms = &fix_ms->data();
        char property_name[200];
        sprintf(property_name,"counter_ms_%s",id);

        if(!ms->prop().getElementProperty< ScalarContainer<int> >(static_cast<const char*>(property_name)))
        {
            (ms->prop().addElementProperty< ScalarContainer<int> >(static_cast<const char*>(property_name),"comm_exchange_borders","frame_invariant", "restart_yes"))->setDefaultValue(-1);
            
        }

        if(delete_atoms_)
            error->fix_error(FLERR,this,"can not use 'delete_atoms' with fix multisphere/*");
        if(!once_)
            error->fix_error(FLERR,this,"must use 'count once' with fix multisphere/*");
    }

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
    
    fix_volumeweight_ms_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("volumeweight_ms","property/atom","scalar",0,0,style,false));

    if (atom->rmass_flag == 0)
        error->fix_error(FLERR,this,"requires atoms have mass");

    if(delete_atoms_ && 1 != atom->map_style)
        error->fix_error(FLERR,this,"requires an atom map of type 'array', via an 'atom_modify map array' command");
/*
    if(!fix_ms_ && static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)))
        error->fix_error(FLERR,this,"fix multisphere must come before fix massflow/mesh in input script");*/
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
    int *mask = atom->mask;
    double *counter = fix_counter_->vector_atom;
    double dot,delta[3]={}, bary[3];
    double mass_this = 0.;
    double nparticles_this = 0.;
    double property_this = 0.;
    double deltan;
    int *tag = atom->tag;

    class FixPropertyAtom* fix_color=static_cast<FixPropertyAtom*>(modify->find_fix_property("color","property/atom","scalar",0,0,style,false));
    bool fixColFound = false;
    if (fix_color) fixColFound=true;

    fix_orientation_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("ex",
        "property/atom","vector",0,0,style,false));

    TriMesh *mesh = fix_mesh_->triMesh();
    int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

    // update reference point
    if (fix_mesh_->triMesh()->isMoving() || fix_mesh_->triMesh()->isDeforming()) {
        setRefPoint();
    }

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
        
        const std::vector<int> & neighborList = fix_neighlist_->get_contact_list(iTri);
        const int numneigh = neighborList.size();

        for(int iNeigh = 0; iNeigh < numneigh; iNeigh++)
        {
            const int iPart = neighborList[iNeigh];

            // skip ghost particles
            if(iPart >= nlocal)
                continue;

            // skip particles not in fix group
            if (!(mask[iPart] & groupbit))
                continue;

            // in case of once_ == true, ignore everything which has been already counted
            if(compDouble(counter[iPart],2.))
                continue;

            int barySign;
            deltan = fix_mesh_->triMesh()->resolveTriSphereContactBary(iPart,iTri,radius[iPart],x[iPart],delta,bary,barySign);

            if(havePointAtOutlet_)
            {
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

            // first-time, just set 0 or 1 depending on what side of the mesh
            if(compDouble(counter[iPart],-1.))
            {
                counter[iPart] = (dot <= 0.) ? 0. : 1.;
                
                continue;
            }

            // particle is now on nvec_ side
            if(dot > 0.  && 7 == barySign) 
            {
                //particle was not on nvec_ side before
                if((compDouble(counter[iPart],0.)) ) // compDouble(counter[iPart],0.))
                {
                    
                    mass_this += (fix_volumeweight_ms_ ? fix_volumeweight_ms_->vector_atom[iPart] : 1.) * rmass[iPart];
                    nparticles_this += (fix_volumeweight_ms_ ? fix_volumeweight_ms_->vector_atom[iPart] : 1.);

                    if(fix_property_)
                    {
                        
                        property_this += fix_property_->vector_atom[iPart];
                    }

                    if(delete_atoms_)
                    {
                        //reset counter to avoid problems with other fixes & mark to be deleted
                        counter[iPart] = -1.0;
                        atom_tags_delete_.push_back(atom->tag[iPart]);
                    }

                    if (screenflag_ && screen)
                        fprintf(screen," %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g \n ",
                                       tag[iPart],2.*radius[iPart]/force->cg(atom->type[iPart]),
                                       x[iPart][0],x[iPart][1],x[iPart][2],
                                       v[iPart][0],v[iPart][1],v[iPart][2]);
                    if(fp_)
                    {
                        fprintf(fp_,"%d", tag[iPart]);

                        if(writeTime_)
                            fprintf(fp_,"  %4.4g ", update->dt*update->ntimestep);

                        fprintf(fp_," %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g",
                                2.*radius[iPart]/force->cg(atom->type[iPart]),
                                   x[iPart][0],x[iPart][1],x[iPart][2],
                                v[iPart][0],v[iPart][1],v[iPart][2]);

                        if(fix_orientation_)
                        {
                            double **orientation = NULL;
                            orientation = fix_orientation_->array_atom;
                            fprintf(fp_,"    %4.4g %4.4g %4.4g ",
                                    orientation[iPart][0], orientation[iPart][1], orientation[iPart][2]);
                        }
                        if (fixColFound)
                            fprintf(fp_,"    %4.0g ", fix_color->vector_atom[iPart]);
                        fprintf(fp_,"\n");
                        fflush(fp_);
                    }
                }

                 if(!delete_atoms_) //only set if not marked for deletion
                    counter[iPart] = once_ ? 2. : 1.;
                
            }
            else if(dot <= 0.) // dot <= 0
            {
                counter[iPart] = 0.;
                
            }
        }
    }

    MPI_Sum_Scalar(mass_this,world);
    MPI_Sum_Scalar(nparticles_this,world);
    if(fix_property_) MPI_Sum_Scalar(property_this,world);

    mass_ += mass_this;
    nparticles_ += nparticles_this;
    property_sum_ += property_this;

}

/* ----------------------------------------------------------------------
   perform particle deletion of marked particles
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixMassflowMesh::pre_exchange()
{
    if (delete_atoms_)
    {
        double mass_deleted_this_ = 0.;
        double nparticles_deleted_this = 0.;
        int *atom_map_array = atom->get_map_array();

        // delete particles

        while (atom_tags_delete_.size() > 0)
        {
            int iPart = atom->map(atom_tags_delete_[0]);

            mass_deleted_this_ += atom->rmass[iPart];
            nparticles_deleted_this += (fix_volumeweight_ms_ ? fix_volumeweight_ms_->vector_atom[iPart] : 1.);

            atom->avec->copy(atom->nlocal-1,iPart,1);

            atom_map_array[atom->tag[atom->nlocal-1]] = iPart;

            atom->nlocal--;

            atom_tags_delete_.erase(atom_tags_delete_.begin());
        }

        atom_tags_delete_.clear();

        MPI_Sum_Scalar(mass_deleted_this_,world);
        MPI_Sum_Scalar(nparticles_deleted_this,world);

        mass_deleted_ += mass_deleted_this_;
        nparticles_deleted_ += nparticles_deleted_this;

        if(nparticles_deleted_this)
        {
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
  list[n++] = nparticles_last_;
  list[n++] = mass_deleted_;
  list[n++] = nparticles_deleted_;

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
  nparticles_last_ = list[n++];
  mass_deleted_ = list[n++];
  nparticles_deleted_ = list[n++];
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
        return nparticles_;
    if(index == 2)
        return delta_t_ == 0. ? 0. : (mass_-mass_last_)/delta_t_;
    if(index == 3)
        return delta_t_ == 0. ? 0. : (nparticles_-nparticles_last_)/delta_t_;
    if(index == 4)
        return mass_deleted_;
    if(index == 5)
        return nparticles_deleted_;
    if(index == 6 && fix_property_)
        return property_sum_;

    return 0.;
}

/* ----------------------------------------------------------------------
    get reference point on face
    calculate normalvec
------------------------------------------------------------------------- */

void FixMassflowMesh::setRefPoint()
{
    fix_mesh_->triMesh()->node(0,0,pref_);
    fix_mesh_->triMesh()->surfaceNorm(0,nvec_);
    double dot = vectorDot3D(nvec_,sidevec_);

    if(fabs(dot) < 1e-6 && !havePointAtOutlet_ )
        error->fix_error(FLERR,this,"need to change 'vec_side', it is currently in or to close to the mesh plane \n"
                                    "This error may be caused by a moving mesh command, since 'vec_side' is not moved with the mesh.");
    else if(dot < 0.)
        vectorScalarMult3D(nvec_,-1.);
}
