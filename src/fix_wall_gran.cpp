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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Philippe Seil (JKU Linz)
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_rigid.h"
#include "fix_mesh.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_contact_property_atom_wall.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "tri_mesh.h"
#include "primitive_wall.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"
#include "contact_interface.h"
#include "fix_property_global.h"
#include "domain_wedge.h"
#include <vector>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_extra_liggghts_nonspherical.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;
using namespace LIGGGHTS::Walls;
using namespace LIGGGHTS::ContactModels;

const double SMALL = 1e-12;

  // modes for conduction contact area calaculation
  // same as in fix_heat_gran_conduction.cpp

  enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
        CONDUCTION_CONTACT_AREA_CONSTANT,
        CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  impl(NULL),
  fix_sum_normal_force_(NULL)
{
    // wall/gran requires gran properties
    // sph not
    if (strncmp(style,"wall/gran",9) == 0 && (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag))
        error->fix_error(FLERR,this,"requires atom attributes radius, omega, torque");

    // defaults
    store_force_ = false;
    store_force_contact_ = false;
    store_force_contact_every_ = 1;
    store_force_contact_stress_ = false;
    stress_flag_ = false;
    n_FixMesh_ = 0;
    dnum_ = 0;

    r0_ = 0.;

    shear_ = 0;
    shearDim_ = shearAxis_ = -1;
    vectorZeroize3D(shearAxisVec_);

    atom_type_wall_ = 1; // will be overwritten during execution, but other fixes require a value here

    // initializations
    fix_wallforce_ = 0;
    fix_wallforce_contact_ = 0;
    fix_wallforce_contact_stress_ = 0;
    fix_store_multicontact_data_ = NULL;
    fix_rigid_ = NULL;
    heattransfer_flag_ = false;

    FixMesh_list_ = NULL;
    primitiveWall_ = NULL;
    fix_history_primitive_ = NULL;

    rebuildPrimitiveNeighlist_ = false;

    addflag_ = 0;
    cwl_ = NULL;

    computeflag_ = 1;

    meshwall_ = -1;

    track_energy_ = false;

    Temp_wall = -1.;
    fixed_contact_area_ = 0.;
    Q = Q_add = 0.;

    area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;

    // parse args
    //style = new char[strlen(arg[2])+2];
    //strcpy(style,arg[2]);

    iarg_ = 3;
    narg_ = narg;

    int nremaining = narg - 3;
    char ** remaining_args = &arg[3];

    // select model for wall/gran (not for SPH)
    if (strncmp(style,"wall/gran",9) == 0)
    {
        int64_t variant = Factory::instance().selectVariant("gran", nremaining, remaining_args,force->custom_contact_models);
        if (variant == -1)
            error->fix_error(FLERR, this, "Invalid model specified (check for typos and enable at least one model)");
        impl = Factory::instance().create("gran", variant, lmp, this);

        if(!impl)
            error->all(FLERR, "Internal errror");
    }

    iarg_ = narg - nremaining;

    bool allow_special_domain_periodic = false;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg_],"primitive") == 0) {
           iarg_++;
           meshwall_ = 0;

           if (meshwall_ == 1)
             error->fix_error(FLERR,this,"'mesh' and 'primitive' are incompatible, choose either of them");

           if (strcmp(arg[iarg_++],"type"))
             error->fix_error(FLERR,this,"expecting keyword 'type'");
           atom_type_wall_ = force->inumeric(FLERR,arg[iarg_++]);
           if (atom_type_wall_ < 1 || atom_type_wall_ > atom->ntypes)
             error->fix_error(FLERR,this,"1 <= type <= max type as defined in create_box'");

           char *wallstyle = arg[iarg_++];
           int nPrimitiveArgs = PRIMITIVE_WALL_DEFINITIONS::numArgsPrimitiveWall(wallstyle);
           
           if(narg-iarg_ < nPrimitiveArgs)
            error->fix_error(FLERR,this,"not enough arguments for primitive wall");

           double * argVec = new double[nPrimitiveArgs];
           for(int i=0;i<nPrimitiveArgs;i++)
           {
             
             argVec[i] = force->numeric(FLERR,arg[iarg_++]);
           }

           bool setflag = false;
           for(int w=0;w<(int)PRIMITIVE_WALL_DEFINITIONS::NUM_WTYPE;w++)
           {
             
             if(strcmp(wallstyle,PRIMITIVE_WALL_DEFINITIONS::wallString[w]) == 0)
             {
               primitiveWall_ = new PrimitiveWall(lmp,(PRIMITIVE_WALL_DEFINITIONS::WallType)w,nPrimitiveArgs,argVec);
               setflag = true;
               break;
             }
           }
           if(!setflag) error->fix_error(FLERR,this,"unknown primitive wall style");
           hasargs = true;
           delete[] argVec;
        } else if (strcmp(arg[iarg_],"mesh") == 0) {
           hasargs = true;
           meshwall_ = 1;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"track_energy") == 0) {
           hasargs = true;
           track_energy_ = true;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"store_force") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_contact_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_contact_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force_contact_'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact_every") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           store_force_contact_every_ = atoi(arg[iarg_+1]);
           if(store_force_contact_every_ <= 0)
                error->fix_error(FLERR,this,"store_force_contact_every must be > 0");
           store_force_contact_ = true;
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"allow_special_domain_periodic") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) allow_special_domain_periodic = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) allow_special_domain_periodic = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'allow_special_domain_periodic'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact_stress") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_contact_stress_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_contact_stress_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force_contact_stress_'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"n_meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'n_meshes'");
          if (iarg_+2 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          n_FixMesh_ = atoi(arg[iarg_+1]);
          if(n_FixMesh_ < 1)
              error->fix_error(FLERR,this,"'n_meshes' > 0 required");
          hasargs = true;
          iarg_ += 2;
        } else if (strcmp(arg[iarg_],"meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'meshes'");
          if(n_FixMesh_ == 0)
              error->fix_error(FLERR,this,"have to define 'n_meshes' before 'meshes'");
          if (narg < iarg_+1+n_FixMesh_)
              error->fix_error(FLERR,this,"not enough arguments");

          FixMesh_list_ = new FixMeshSurface*[n_FixMesh_];
          for(int i = 1; i <= n_FixMesh_; i++)
          {
              int f_i = modify->find_fix(arg[iarg_+i]);
              if (f_i == -1)
                  error->fix_error(FLERR,this,"could not find fix mesh id you provided");
              if (strncmp(modify->fix[f_i]->style,"mesh/surface",12))
                  error->fix_error(FLERR,this,"the fix belonging to the id you provided is not of type mesh");
              FixMesh_list_[i-1] = static_cast<FixMeshSurface*>(modify->fix[f_i]);

              if(FixMesh_list_[i-1]->trackStress())
                stress_flag_ = true;
              
          }
          hasargs = true;
          iarg_ += 1+n_FixMesh_;
        } else if (strcmp(arg[iarg_],"shear") == 0) {
          if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'shear'");
          if(!primitiveWall_)
            error->fix_error(FLERR,this,"have to define primitive wall before 'shear'. For mehs walls, please use fix move/mesh");

          if (strcmp(arg[iarg_+1],"x") == 0) shearDim_ = 0;
          else if (strcmp(arg[iarg_+1],"y") == 0) shearDim_ = 1;
          else if (strcmp(arg[iarg_+1],"z") == 0) shearDim_ = 2;
          else error->fix_error(FLERR,this,"illegal 'shear' dim");
          vshear_ = force->numeric(FLERR,arg[iarg_+2]);
          shear_ = 1;

          // update axis for cylinder etc if needed
          if(shearDim_ != primitiveWall_->axis())
          {
            shearAxis_ = primitiveWall_->axis();
            shearAxisVec_[shearAxis_] = vshear_;
          }

          hasargs = true;
          iarg_ += 3;
        } else if (strcmp(arg[iarg_],"temperature") == 0) {
            if (iarg_+1 >= narg)
              error->fix_error(FLERR,this,"not enough arguments for 'temperature'");
            Temp_wall = force->numeric(FLERR,arg[iarg_+1]);
            hasargs = true;
            iarg_ += 2;
        } else if(strcmp(arg[iarg_],"contact_area") == 0) {

          if(strcmp(arg[iarg_+1],"overlap") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
          else if(strcmp(arg[iarg_+1],"projection") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
          else if(strcmp(arg[iarg_+1],"constant") == 0)
          {
            if (iarg_+3 > narg)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
            fixed_contact_area_ = force->numeric(FLERR,arg[iarg_+2]);
            if (fixed_contact_area_ <= 0.)
                error->fix_error(FLERR,this,"'contact_area constant' value must be > 0");
            iarg_++;
          }
          else error->fix_error(FLERR,this,"expecting 'overlap', 'projection' or 'constant' after 'contact_area'");
          iarg_ += 2;
          hasargs = true;
        }
    }

    if(impl)
      impl->settings(narg - iarg_, &arg[iarg_], this);

    // error checks

    if(meshwall_ == -1 && primitiveWall_ == 0)
        error->fix_error(FLERR,this,"Need to use define style 'mesh' or 'primitive'");

    if(primitiveWall_ && (modify->n_fixes_style("particletemplate/convexhull") > 0 || modify->n_fixes_style("particletemplate/concave") > 0))
        error->fix_error(FLERR,this,"style 'primitive' is not compatible with convex or concave particles");

    if(meshwall_ == 1 && !FixMesh_list_)
        error->fix_error(FLERR,this,"Need to provide the number and a list of meshes by using 'n_meshes' and 'meshes'");

    if (modify->find_fix_style("continuum/weighted", 0))
        store_force_contact_stress_ = true;

    if (!allow_special_domain_periodic                               &&
        (domain->triclinic || dynamic_cast<DomainWedge*>(domain))     &&
        (domain->xperiodic || domain->yperiodic || domain->zperiodic) )
        error->fix_error(FLERR, this, "Triclinic or wedge domain is not allowed with periodic boundary conditions and meshes. This can be overridden by using the allow_special_domain_periodic option of fix wall/gran. In this case the user must ensure that meshes are sufficiently far away from periodic boundaries");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_create()
{
    if(strncmp(style,"wall/gran",9) != 0)
    {
      // case non-granular (sph)
      dnum_ = 0;
    }

    // register storage for wall force if required
    if(store_force_)
    {
          char *wallforce_name = new char[strlen(style)+1+6];
          strcpy(wallforce_name,"force_");
          strcat(wallforce_name,id);
          const char *fixarg[11];
          fixarg[0] = wallforce_name;
          fixarg[1] = "all";
          fixarg[2] = "property/atom";
          fixarg[3] = wallforce_name;
          fixarg[4] = "vector";
          fixarg[5] = "no";    // restart
          fixarg[6] = "no";    // communicate ghost
          fixarg[7] = "no";    // communicate rev
          fixarg[8] = "0.";
          fixarg[9] = "0.";
          fixarg[10] = "0.";
          modify->add_fix(11,const_cast<char**>(fixarg));
          fix_wallforce_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(wallforce_name,"property/atom","vector",3,0,style));
          delete []wallforce_name;
   }

   if(store_force_contact_ && 0 == meshwall_)
   {
        const char *fixarg[19];
        char fixid[200],ownid[200];
        sprintf(fixid,"contactforces_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]="all";
        fixarg[2]="contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]="6";
        fixarg[5]="fx";
        fixarg[6]="0";
        fixarg[7]="fy";
        fixarg[8]="0";
        fixarg[9]="fz";
        fixarg[10]="0";
        fixarg[11]="tx";
        fixarg[12]="0";
        fixarg[13]="ty";
        fixarg[14]="0";
        fixarg[15]="tz";
        fixarg[16]="0";
        fixarg[17]="primitive";
        fixarg[18]=ownid;
        modify->add_fix(19,const_cast<char**>(fixarg));
        fix_wallforce_contact_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
   }

   if(store_force_contact_stress_ && 0 == meshwall_)
   {
        const char *fixarg[25];
        char fixid[200],ownid[200];
        sprintf(fixid,"contactforces_stress_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]="all";
        fixarg[2]="contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]="9";
        fixarg[5]="fx";
        fixarg[6]="0";
        fixarg[7]="fy";
        fixarg[8]="0";
        fixarg[9]="fz";
        fixarg[10]="0";
        fixarg[11]="deltax";
        fixarg[12]="0";
        fixarg[13]="deltay";
        fixarg[14]="0";
        fixarg[15]="deltaz";
        fixarg[16]="0";
        fixarg[17]="vx";
        fixarg[18]="0";
        fixarg[19]="vy";
        fixarg[20]="0";
        fixarg[21]="vz";
        fixarg[22]="0";
        fixarg[23]="primitive";
        fixarg[24]=ownid;
        modify->add_fix(25,const_cast<char**>(fixarg));
        fix_wallforce_contact_stress_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
   }

   // create neighbor list for each mesh
   
   for(int i=0;i<n_FixMesh_;i++)
   {
       
       FixMesh_list_[i]->createWallNeighList(igroup);
       FixMesh_list_[i]->createContactHistory(dnum());

       if(store_force_contact_)
         FixMesh_list_[i]->createMeshforceContact();

       if(store_force_contact_stress_)
         FixMesh_list_[i]->createMeshforceContactStress();
   }

   // contact history for primitive wall
   if(meshwall_ == 0 && dnum_ > 0)
   {
          char *hist_name = new char[strlen(id)+1+10];
          strcpy(hist_name,"history_");
          strcat(hist_name,id);
          const char **fixarg = new const char*[8+dnum_];
          fixarg[0] = hist_name;
          fixarg[1] = "all";
          fixarg[2] = "property/atom";
          fixarg[3] = hist_name;
          if (dnum_ > 1)
              fixarg[4] = "vector";
          else
              fixarg[4] = "vector_one_entry";
          fixarg[5] = "yes";    // restart
          fixarg[6] = "no";    // communicate ghost
          fixarg[7] = "no";    // communicate rev
          for(int i = 8; i < 8+dnum_; i++)
              fixarg[i] = "0.";
          modify->add_fix(8+dnum_,const_cast<char**>(fixarg));
          fix_history_primitive_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name,"property/atom","vector",dnum_,0,style));
          delete []fixarg;
          delete []hist_name;
   }

}

/* ---------------------------------------------------------------------- */

void FixWallGran::pre_delete(bool unfixflag)
{
    if(unfixflag && store_force_)
        modify->delete_fix(fix_wallforce_->id);
    if(unfixflag && fix_history_primitive_)
        modify->delete_fix(fix_history_primitive_->id);

    if(unfixflag && store_force_contact_ && 0 == meshwall_)
        modify->delete_fix(fix_wallforce_contact_->id);

    if(unfixflag && store_force_contact_stress_ && 0 == meshwall_)
        modify->delete_fix(fix_wallforce_contact_stress_->id);

    if(unfixflag && cwl_)
       error->fix_error(FLERR,this,"need to uncompute the active compute wall/gran/local before unfixing the wall");

    if(unfixflag)
    {
       for(int i=0;i<n_FixMesh_;i++)
       {
           
           FixMesh_list_[i]->deleteWallNeighList();
           FixMesh_list_[i]->deleteContactHistory();

           if(store_force_contact_)
             FixMesh_list_[i]->deleteMeshforceContact();

           if(store_force_contact_stress_)
             FixMesh_list_[i]->deleteMeshforceContactStress();
       }
    }
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
    if(primitiveWall_ != 0) delete primitiveWall_;
    if(FixMesh_list_) delete []FixMesh_list_;
    delete impl;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::min_type()
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::max_type()
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

PrimitiveWall* FixWallGran::primitiveWall()
{ return primitiveWall_; }

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
    dt_ = update->dt;

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        // check if a fix rigid is registered - important for damp
        fix_rigid_ = static_cast<FixRigid*>(modify->find_fix_style_strict("rigid",0));

        if (strcmp(update->integrate_style,"respa") == 0)
          nlevels_respa_ = ((Respa *) update->integrate)->nlevels;

        if(impl)
          impl->init_granular();
        else
        {
          // init for derived classes
          init_granular();
        }

        // disallow more than one wall of non-primitive style
        
        if(is_mesh_wall())
        {
            int nfix = modify->n_fixes_style("wall/gran");
            for (int ifix = 0; ifix < nfix; ifix++)
            {
                FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
                if (fwg == this) continue;
                if (fwg->is_mesh_wall())
                    error->fix_error(FLERR,this,"More than one wall of type 'mesh' is not supported");
            }
        }
    }
    
    fix_sum_normal_force_ =
        static_cast<FixPropertyAtom*>
        (
            modify->find_fix_property("sum_normal_force_","property/atom","scalar",0,0,style, false)
        );
}

void FixWallGran::createMulticontactData()
{
    if(fix_store_multicontact_data_ == NULL && 0 == meshwall_)
    {
        // create a new per contact property which will contain the data for the computation according to Brodu et. al. 2016
        // surfPosIJ will contain the position of the contact surface ij, realtive to position i
        // normalForce will contain the normal component of the contact force
        const char *fixarg[17];
        char fixid[200], ownid[200];
        sprintf(fixid,"multicontactData_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]="all";
        fixarg[2]="contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]="4";
        fixarg[5]="surfPosIJ_x";
        fixarg[6]="0";
        fixarg[7]="surfPosIJ_y";
        fixarg[8]="0";
        fixarg[9]="surfPosIJ_z";
        fixarg[10]="0";
        fixarg[11]="normalForce";
        fixarg[12]="0";
        fixarg[13]="primitive";
        fixarg[14]=ownid;
        fixarg[15]="reset";
        fixarg[16]="no";
        modify->add_fix(17,const_cast<char**>(fixarg));
        fix_store_multicontact_data_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
    }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
    {
      pre_neighbor();
      pre_force(vflag);
      post_force(vflag);
    }
    else
    {
      ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa_-1);
      post_force_respa(vflag,nlevels_respa_-1,0);
      ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa_-1);
    }

    init_heattransfer();
}

/* ----------------------------------------------------------------------
   neighbor list via fix wall/gran is only relevant for primitive walls
------------------------------------------------------------------------- */

void FixWallGran::pre_neighbor()
{
    rebuildPrimitiveNeighlist_ = (primitiveWall_ != 0);
}

void FixWallGran::pre_force(int vflag)
{
    int nlocal = atom->nlocal;

    x_ = atom->x;
    radius_ = atom->radius;
    cutneighmax_ = neighbor->cutneighmax;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(atom->superquadric_flag) {
        quat_ = atom->quaternion;
        shape_ = atom->shape;
        blockiness_ = atom->blockiness;
    }
#endif

    // build neighlist for primitive walls
    
    if(rebuildPrimitiveNeighlist_)
      primitiveWall_->buildNeighList(radius_ ? neighbor->skin:(r0_+neighbor->skin),x_,radius_,nlocal);

    rebuildPrimitiveNeighlist_ = false;
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via verlet
------------------------------------------------------------------------- */

void FixWallGran::post_force(int vflag)
{
    computeflag_ = 1;
    shearupdate_ = 1;
    if (update->setupflag) shearupdate_ = 0;
    addflag_ = 0;

    post_force_wall(vflag);

}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via compute wall/gran
------------------------------------------------------------------------- */

void FixWallGran::post_force_pgl()
{
    computeflag_ = 0;
    shearupdate_ = 0;
    addflag_ = 1;

    post_force_wall(0);

    // call finalize for cwl
    if(cwl_)
    {
        
        cwl_->pair_finalize();
    }
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */

void FixWallGran::post_force_wall(int vflag)
{
  // set pointers and values appropriately
  nlocal_ = atom->nlocal;
  x_ = atom->x;
  f_ = atom->f;
  radius_ = atom->radius;
  rmass_ = atom->rmass;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  if(atom->superquadric_flag) {
    quat_ = atom->quaternion;
    shape_ = atom->shape;
    blockiness_ = atom->blockiness;
  }
#endif

  if(fix_rigid_)
  {
      body_ = fix_rigid_->body;
      masstotal_ = fix_rigid_->masstotal;
  }

  if(fix_wallforce_)
    wallforce_ = fix_wallforce_->array_atom;

  cutneighmax_ = neighbor->cutneighmax;

  if(nlocal_ && !radius_ && r0_ == 0.)
    error->fix_error(FLERR,this,"need either per-atom radius or r0_ being set");

  if(store_force_)
  {
      for(int i = 0; i < nlocal_; i++)
      {
          vectorZeroize3D(wallforce_[i]);
      }
  }
  if(meshwall_ == 1)
    post_force_mesh(vflag);
  else
    post_force_primitive(vflag);

  if(meshwall_ == 0 && store_force_contact_)
    fix_wallforce_contact_->do_forward_comm();

  if(meshwall_ == 0 && store_force_contact_stress_)
    fix_wallforce_contact_stress_->do_forward_comm();

  if(meshwall_ == 1 && store_force_contact_)
  {
    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
        FixMesh_list_[imesh]->meshforceContact()->do_forward_comm();
  }

  if(meshwall_ == 1 && store_force_contact_stress_)
  {
    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
        FixMesh_list_[imesh]->meshforceContactStress()->do_forward_comm();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
    if (ilevel == nlevels_respa_-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_mesh(int vflag)
{
    
    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    MultiVectorContainer<double,3,3> *vMeshC;
    double *radius = atom->radius;
    double ***vMesh = 0;
    int nlocal = atom->nlocal;
    int nTriAll, barysign = -1;
    const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;

    SurfacesIntersectData sidata;
    sidata.is_wall = true;

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // mark all contacts for delettion at this point
      
      if(fix_contact) fix_contact->markAllContacts();

      if(store_force_contact_)
        fix_wallforce_contact_ = FixMesh_list_[iMesh]->meshforceContact();

      if(store_force_contact_stress_)
        fix_wallforce_contact_stress_ = FixMesh_list_[iMesh]->meshforceContactStress();

      fix_store_multicontact_data_ = FixMesh_list_[iMesh]->meshMulticontactData();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = FixMesh_list_[iMesh]->meshNeighlist();

      // moving mesh
      vectorZeroize3D(v_wall);
      vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
      if(vMeshC)
        vMesh = vMeshC->begin();

      atom_type_wall_ = FixMesh_list_[iMesh]->atomTypeWall();

      // loop owned and ghost triangles
      for(int iTri = 0; iTri < nTriAll; iTri++)
      {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();
          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            
            const int iPart = neighborList[iCont];

            // do not handle ghost particles
            if (iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

            #ifdef SUPERQUADRIC_ACTIVE_FLAG
                if(atom->superquadric_flag) {
                  #ifdef LIGGGHTS_DEBUG
                    if(std::isnan(vectorMag3D(x_[iPart])))
                      error->fix_error(FLERR,this,"x_[iPart] is NaN!");
                    if(std::isnan(vectorMag4D(quat_[iPart])))
                      error->fix_error(FLERR,this,"quat_[iPart] is NaN!");
                  #endif

                  Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);

                  if(mesh->sphereTriangleIntersection(iTri, radius_[iPart], x_[iPart])) //check for Bounding Sphere-triangle intersection
                  {
                    deltan = mesh->resolveTriSuperquadricContact(iTri, delta, sidata.contact_point, particle, bary);
                    #ifdef LIGGGHTS_DEBUG
                        if(std::isnan(deltan))
                          error->fix_error(FLERR,this,"deltan is NaN!");
                        if(std::isnan(vectorMag3D(delta)))
                          error->fix_error(FLERR,this,"delta is NaN!");
                        if(std::isnan(vectorMag3D(sidata.contact_point)))
                          error->fix_error(FLERR,this,"sidata.contact_point is NaN!");
                    #endif
                  }
                  else
                    deltan = LARGE_TRIMESH;
                  sidata.is_non_spherical = true; //by default it is false
                } else {
                  sidata.radi = radius_ ? radius_[iPart] : r0_;
                  if (fix_store_multicontact_data_)
                  {
                      double * deltaData = NULL;
                      const bool contact = fix_store_multicontact_data_->haveContact(iPart, idTri, deltaData);
                      if (contact)
                          sidata.radi += deltaData[3];
                  }
                  deltan = mesh->resolveTriSphereContactBary(iPart, iTri, sidata.radi, x_[iPart], delta, bary, barysign, atom->shapetype_flag ? false : true);
                }
            #else
                sidata.radi = radius_ ? radius_[iPart] : r0_;
                if (fix_store_multicontact_data_)
                {
                    double * deltaData = NULL;
                    const bool contact = fix_store_multicontact_data_->haveContact(iPart, idTri, deltaData);
                    if (contact)
                        sidata.radi += deltaData[3];
                }
                
                deltan = mesh->resolveTriSphereContactBary(iPart, iTri, sidata.radi, x_[iPart], delta, bary, barysign, atom->shapetype_flag ? false : true);
            #endif
            
            if(deltan > cutneighmax_) continue;

            sidata.i = iPart;

            bool intersectflag = (deltan <= 0);

            sidata.mesh = mesh;

            if(atom->shapetype_flag)
            {
                
                sidata.j = iTri;
                fix_contact->handleContact(iPart,idTri,sidata.contact_history,intersectflag,false);
                if(vMeshC)
                {
                    for(int i = 0; i < 3; i++)
                        v_wall[i] = (bary[0]*vMesh[iTri][0][i] +
                                     bary[1]*vMesh[iTri][1][i] +
                                     bary[2]*vMesh[iTri][2][i] );
                }
                sidata.v_i = atom->v[iPart];
                sidata.omega_i = atom->omega[iPart];
                sidata.v_j = v_wall;
                sidata.shearupdate = shearupdate_;
                sidata.computeflag = computeflag_;
                intersectflag = impl->checkSurfaceIntersect(sidata);
                deltan = -sidata.deltan;
                
            }

            sidata.fix_mesh = FixMesh_list_[iMesh];

            if(deltan <= 0 || (radius && deltan < contactDistanceMultiplier*radius[iPart]))
            {
              
              if(!atom->shapetype_flag && fix_contact && ! fix_contact->handleContact(iPart,idTri,sidata.contact_history,intersectflag,7 == barysign)) continue;

              if(vMeshC && !atom->shapetype_flag)
              {
                for(int i = 0; i < 3; i++)
                    v_wall[i] = (bary[0]*vMesh[iTri][0][i] + bary[1]*vMesh[iTri][1][i] + bary[2]*vMesh[iTri][2][i]);
              }

              if(!sidata.is_non_spherical || atom->superquadric_flag)
                sidata.deltan   = -deltan;
              sidata.delta[0] = -delta[0];
              sidata.delta[1] = -delta[1];
              sidata.delta[2] = -delta[2];
              if(impl)
                impl->compute_force(this, sidata, intersectflag,v_wall,FixMesh_list_[iMesh],iMesh,mesh,iTri);
              else
              {
                sidata.r =  r0_ - sidata.deltan;
                compute_force(sidata, v_wall); // LEGACY CODE (SPH)
              }
            }
          }
      }

      // clean-up contacts
      
      if(fix_contact) fix_contact->cleanUpContacts();
    }
}

/* ----------------------------------------------------------------------
   post_force for primitive wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_primitive(int vflag)
{
  int *mask = atom->mask;

  SurfacesIntersectData sidata;
  sidata.is_wall = true;
  double *radius = atom->radius;
  const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;
  
  // contact properties
  double delta[3]={},deltan,rdist[3];
  double v_wall[] = {0.,0.,0.};
  double **c_history = 0;

  if(dnum() > 0)
    c_history = fix_history_primitive_->array_atom;

  // if shear, set velocity accordingly
  if (shear_) v_wall[shearDim_] = vshear_;

  // loop neighbor list
  int *neighborList;
  int nNeigh = primitiveWall_->getNeighbors(neighborList);

  for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++)
  {
    int iPart = *neighborList;

    if(!(mask[iPart] & groupbit)) continue;

    sidata.radi = radius_ ? radius_[iPart] : r0_;
    if (fix_store_multicontact_data_)
    {
        double * deltaData = NULL;
        const bool contact = fix_store_multicontact_data_->haveContact(iPart, 1, deltaData);
        if (contact)
            sidata.radi += deltaData[3];
    }
    deltan = primitiveWall_->resolveContact(x_[iPart], sidata.radi, delta);

    if(deltan>cutneighmax_) continue;

    if(deltan <= 0 || deltan < contactDistanceMultiplier*radius[iPart])
    {
      // spheres
      bool intersectflag = (deltan <= 0);

      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(atom->superquadric_flag) {
            double sphere_contact_point[3];
            vectorAdd3D(x_[iPart], delta, sphere_contact_point);
            double closestPoint[3], closestPointProjection[3], point_of_lowest_potential[3];
            #ifdef LIGGGHTS_DEBUG
              if(std::isnan(vectorMag3D(x_[iPart])))
                error->fix_error(FLERR,this,"x_[iPart] is NaN!");
              if(std::isnan(vectorMag4D(quat_[iPart])))
                error->fix_error(FLERR,this,"quat_[iPart] is NaN!");
            #endif

            Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);
            intersectflag = particle.plane_intersection(delta, sphere_contact_point, closestPoint, point_of_lowest_potential);
            #ifdef LIGGGHTS_DEBUG
                if(std::isnan(vectorMag3D(delta)))
                  error->fix_error(FLERR,this,"delta is NaN!");
                if(std::isnan(vectorMag3D(sphere_contact_point)))
                  error->fix_error(FLERR,this,"sphere_contact_point is NaN!");
                if(std::isnan(vectorMag3D(closestPoint)))
                  error->fix_error(FLERR,this,"closestPoint is NaN!");
                if(std::isnan(vectorMag3D(point_of_lowest_potential)))
                  error->fix_error(FLERR,this,"point_of_lowest_potential is NaN!");
            #endif
            deltan = -MathExtraLiggghtsNonspherical::point_wall_projection(delta, sphere_contact_point, closestPoint, closestPointProjection);

            #ifdef LIGGGHTS_DEBUG
                if(std::isnan(deltan))
                  error->fix_error(FLERR,this,"deltan is NaN!");
                if(std::isnan(vectorMag3D(delta)))
                  error->fix_error(FLERR,this,"delta is NaN!");
                if(std::isnan(vectorMag3D(closestPointProjection)))
                  error->fix_error(FLERR,this,"closestPointProjection is NaN!");
                if(std::isnan(vectorMag3D(closestPoint)))
                  error->fix_error(FLERR,this,"closestPoint is NaN!");
            #endif

            vectorCopy3D(closestPoint, sidata.contact_point);
            sidata.is_non_spherical = true; //by default it is false
      }
      #endif

      if(shear_ && shearAxis_ >= 0)
      {
          primitiveWall_->calcRadialDistance(x_[iPart],rdist);
          vectorCross3D(shearAxisVec_,rdist,v_wall);
          
      }
      sidata.i = iPart;
      sidata.contact_history = c_history ? c_history[iPart] : NULL;

      if(    (atom->superquadric_flag && deltan > 0.0)
          || (atom->shapetype_flag && !impl->checkSurfaceIntersect(sidata)))
      {
        if(c_history)
            vectorZeroizeN(c_history[iPart],dnum_);
        break;
      }

      if(!sidata.is_non_spherical || atom->superquadric_flag)
        sidata.deltan   = -deltan;
      sidata.delta[0] = -delta[0];
      sidata.delta[1] = -delta[1];
      sidata.delta[2] = -delta[2];

      if(impl)
        impl->compute_force(this, sidata, intersectflag, v_wall);
      else
      {
        sidata.r =  r0_ - sidata.deltan;
        compute_force(sidata, v_wall); // LEGACY CODE (SPH)
      }
    }
    
    else
    {
      if(c_history)
      {
         vectorZeroizeN(c_history[iPart],dnum_);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::compute_force(SurfacesIntersectData &, double *)
{
    
}

/* ---------------------------------------------------------------------- */

int FixWallGran::is_moving()
{
    if(is_mesh_wall())
    {
        for(int i = 0; i < n_FixMesh_; i++) {
            if(FixMesh_list_[i]->mesh()->isMoving())
               return 1;
        }
        return 0;
    }
    return shear_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int &nIntersect)
{
    
    if (!is_mesh_wall()) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(nIntersect);

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int &nIntersect)
{
    int local_nIntersect;
    int ncontacts = n_contacts_local(local_nIntersect);
    MPI_Sum_Scalar(ncontacts,world);
    MPI_Sum_Scalar(local_nIntersect,world);
    nIntersect = local_nIntersect;
    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int contact_groupbit, int &nIntersect)
{
    
    if (!is_mesh_wall()) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(contact_groupbit, nIntersect);

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int contact_groupbit, int &nIntersect)
{
    int local_nIntersect;
    int ncontacts = n_contacts_local(contact_groupbit, local_nIntersect);
    MPI_Sum_Scalar(ncontacts,world);
    MPI_Sum_Scalar(local_nIntersect,world);
    nIntersect = local_nIntersect;
    return ncontacts;
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixWallGran::register_compute_wall_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cwl_ != NULL)
     error->fix_error(FLERR,this,"Fix wall/gran allows only one compute of type wall/gran/local");
   cwl_ = ptr;
   dnum_compute = dnum_; //history values
}

void FixWallGran::unregister_compute_wall_local(ComputePairGranLocal *ptr)
{
   if(cwl_ != ptr)
     error->fix_error(FLERR,this,"Illegal situation in FixWallGran::unregister_compute_wall_local");
   cwl_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init_heattransfer()
{
    fppa_T = NULL;
    fppa_hf = NULL;
    fppa_htcw = NULL;
    deltan_ratio = NULL;

    // decide if heat transfer is to be calculated

    if (!is_mesh_wall() && Temp_wall < 0.) return;
    else if (is_mesh_wall())
    {
        int heatflag = 0;
        for(int imesh = 0; imesh < n_meshes(); imesh++)
        {
            heatflag = heatflag || mesh_list()[imesh]->mesh()->prop().getGlobalProperty<ScalarContainer<double> >("Temp") != NULL;
        }

        if(!heatflag) return;
    }

    // heat transfer is to be calculated - continue with initializations

    // set flag so addHeatFlux function is called
    heattransfer_flag_ = true;

    if(modify->n_fixes_style("heat/gran/conduction/fast") > 0)
        heattransfer_flag_ = false;

    // if(screen && comm->me == 0) fprintf(screen,"Initializing wall/gran heat transfer model\n");
    fppa_T = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,style));
    fppa_hf = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",1,0,style));
    fppa_htcw = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallHeattransferCoeff","property/atom","scalar",1,0,style,false));

    th_cond = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",0,0,style))->get_values();

    // if youngsModulusOriginal defined, get deltan_ratio
    Fix* ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style,false);
    // deltan_ratio is defined by heat transfer fix, see if there is one
    int n_htf = modify->n_fixes_style("heat/gran/conduction");

    // get deltan_ratio set by the heat transfer fix
    if(ymo_fix && n_htf) deltan_ratio = static_cast<FixPropertyGlobal*>(ymo_fix)->get_array_modified();
}

/* ---------------------------------------------------------------------- */

void FixWallGran::wall_temperature_unique(bool &has_temp,bool &temp_unique, double &temperature_unique)
{
    has_temp = false;
    temp_unique = true;
    temperature_unique = 0.;

    if(Temp_wall > 0)
    {
        has_temp = true;
        temp_unique = true;
        temperature_unique = Temp_wall;

    }

    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
    {
        
        if( (mesh_list()[imesh]->triMesh())->prop().getElementProperty<ScalarContainer<double> >("Temp"))
        {
            has_temp = true;
            temp_unique = false;
            
            return;
        }
        if( (mesh_list()[imesh]->triMesh())->prop().getGlobalProperty<ScalarContainer<double> >("Temp"))
        {
            has_temp = true;

            double Temp_mesh = (*((mesh_list()[imesh]->triMesh())->prop().getGlobalProperty< ScalarContainer<double> >("Temp")))(0);

            if(temperature_unique > 0. && Temp_mesh != temperature_unique)
            {
                temp_unique = false;
                return;
            }
            temperature_unique = Temp_mesh;
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::addHeatFlux(TriMesh *mesh,int ip, const double ri, double delta_n, double area_ratio)
{
    
    //r is the distance between the sphere center and wall
    double tcop, tcowall, hc, Acont=0.0, r;
    double reff_wall = ri;
    int itype = atom->type[ip];

    if(mesh)
    {
        ScalarContainer<double> *temp_ptr = mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp");
        if (!temp_ptr)
            return;

        Temp_wall = (*temp_ptr)(0);
    }

    double *Temp_p = fppa_T->vector_atom;
    double *heatflux = fppa_hf->vector_atom;

    if(CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    {
        
        if(deltan_ratio)
           delta_n *= deltan_ratio[itype-1][atom_type_wall_-1];

        r = ri - delta_n;

        Acont = (reff_wall*reff_wall-r*r)*M_PI*area_ratio; //contact area sphere-wall
    }
    else if (CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
        Acont = fixed_contact_area_;
    else if (CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    {
        Acont = M_PI*ri*ri;
    }

    tcop = th_cond[itype-1]; //types start at 1, array at 0
    tcowall = th_cond[atom_type_wall_-1];

    if ((fabs(tcop) < SMALL) || (fabs(tcowall) < SMALL)) hc = 0.;
    else hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(Acont);

    if(computeflag_)
    {
        double hf = (Temp_wall-Temp_p[ip]) * hc;
        heatflux[ip] += hf;
        Q_add += hf * update->dt;
        
        if(fppa_htcw)
            fppa_htcw->vector_atom[ip] = hc;
    }
    if(cwl_ && addflag_)
        cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
    
}
