/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
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
#include "fix_contact_history.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface_stress.h"
#include "primitive_wall.h"
#include "tri_mesh.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    // wall/gran requires gran properties
    // sph not
    if (strncmp(style,"wall/gran",9) == 0 && (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag))
        error->fix_error(FLERR,this,"requires atom attributes radius, omega, torque");

    // defaults
    store_force_ = false;
    stress_flag_ = false;
    n_FixMesh_ = 0;
    dnum_ = 0;
    skinDistance_ = 0.0;

    r0_ = 0.;

    shear_ = 0;
    shearDim_ = shearAxis_ = -1;
    vectorZeroize3D(shearAxisVec_);

    atom_type_wall_ = 1; // will be overwritten during execution, but other fixes require a value here

    // initializations
    pairgran_ = NULL;
    fix_wallforce_ = NULL;
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

    // parse args
    //style = new char[strlen(arg[2])+2];
    //strcpy(style,arg[2]);

    iarg_ = 3;
    narg_ = narg;

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
           atom_type_wall_ = atoi(arg[iarg_++]);
           if (atom_type_wall_ < 1 || atom_type_wall_ > atom->ntypes)
             error->fix_error(FLERR,this,"1 <= type <= max type as defined in create_box'");

           char *wallstyle = arg[iarg_++];
           int nPrimitiveArgs = PRIMITIVE_WALL_DEFINITIONS::numArgsPrimitiveWall(wallstyle);
           
           if(narg-iarg_ < nPrimitiveArgs)
            error->fix_error(FLERR,this,"not enough arguments for primitive wall");

           double * argVec = new double[nPrimitiveArgs];
           for(int i=0;i<nPrimitiveArgs;i++)
           {
             
             argVec[i] = force->numeric(arg[iarg_++]);
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
        } else if (strcmp(arg[iarg_],"store_force") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force'");
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
          vshear_ = atof(arg[iarg_+2]);
          shear_ = 1;

          // update axis for cylinder etc if needed
          if(shearDim_ != primitiveWall_->axis())
          {
            shearAxis_ = primitiveWall_->axis();
            shearAxisVec_[shearAxis_] = vshear_;
          }

          hasargs = true;
          iarg_ += 3;
        }
    }

    // error checks

    if(meshwall_ == -1 && primitiveWall_ == 0)
        error->fix_error(FLERR,this,"Need to use define style 'mesh' or 'primitive'");

    if(meshwall_ == 1 && !FixMesh_list_)
        error->fix_error(FLERR,this,"Need to provide the number and a list of meshes by using 'n_meshes' and 'meshes'");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_create()
{
    if(iarg_ < narg_)
        error->fix_error(FLERR,this,"invalid keyword or keyword(s) in wrong order");

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        //check if wall style and pair style fit together
        char *pairstyle = new char[strlen(style)-9+1];
        strcpy(pairstyle,&style[9]);

        PairGran *oldpair = pairgran_;
        pairgran_ = (PairGran*)force->pair_match(pairstyle,0);
        bool pair_changed = (pairgran_!=oldpair);

        delete []pairstyle;

        if (!pairgran_)
          error->fix_error(FLERR,this,"Fix wall/gran style and pair/gran style have to match, you have to define the wall after the pair style");

        if(dnum_ && dnum_ != pairgran_->dnum_pair())
          error->fix_error(FLERR,this,"Number of contact history values of pair style and wall style does not match");
        else dnum_ = pairgran_->dnum_pair();
    }
    // case non-granular (sph)
    else
        dnum_ = 0;

    // register storage for wall force if required
    if(store_force_)
    {
          char *wallforce_name = new char[strlen(style)+1+6];
          strcpy(wallforce_name,"force_");
          strcat(wallforce_name,id);
          char **fixarg = new char*[11];
          fixarg[0] = wallforce_name;
          fixarg[1] = (char *) "all";
          fixarg[2] = (char *) "property/atom";
          fixarg[3] = wallforce_name;
          fixarg[4] = (char *) "vector";
          fixarg[5] = (char *) "no";    // restart
          fixarg[6] = (char *) "no";    // communicate ghost
          fixarg[7] = (char *) "no";    // communicate rev
          fixarg[8] = (char *) "0.";
          fixarg[9] = (char *) "0.";
          fixarg[10] = (char *) "0.";
          modify->add_fix(11,fixarg);
          fix_wallforce_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(wallforce_name,"property/atom","vector",3,0,style));
          delete []fixarg;
          delete []wallforce_name;
   }

   // create neighbor list for each mesh
   // also create contact tracker
   for(int i=0;i<n_FixMesh_;i++)
   {
      FixMesh_list_[i]->createWallNeighList(igroup);
      //if(dnum()>0)FixMesh_list_[i]->createContactHistory(dnum());
      FixMesh_list_[i]->createContactHistory(dnum());
   }

   // contact history for primitive wall
   if(meshwall_ == 0 && dnum_ > 0)
   {
          char *hist_name = new char[strlen(style)+1+10];
          strcpy(hist_name,"history_");
          strcat(hist_name,id);
          char **fixarg = new char*[8+dnum_];
          fixarg[0] = hist_name;
          fixarg[1] = (char *) "all";
          fixarg[2] = (char *) "property/atom";
          fixarg[3] = hist_name;
          fixarg[4] = (char *) "vector";
          fixarg[5] = (char *) "yes";    // restart
          fixarg[6] = (char *) "no";    // communicate ghost
          fixarg[7] = (char *) "no";    // communicate rev
          for(int i = 8; i < 8+dnum_; i++)
              fixarg[i] = (char *) "0.";
          modify->add_fix(8+dnum_,fixarg);
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
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
    if(primitiveWall_ != 0) delete primitiveWall_;
    if(FixMesh_list_) delete []FixMesh_list_;
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

void FixWallGran::init()
{
    dt_ = update->dt;

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        //check if wall style and pair style fit together
        char *pairstyle = new char[strlen(style)-9+1];
        strcpy(pairstyle,&style[9]);

        PairGran *oldpair = pairgran_;
        pairgran_ = (PairGran*)force->pair_match(pairstyle,0);
        bool pair_changed = (pairgran_ != oldpair);

        delete []pairstyle;

        // prohibit changing pair style with wall active
        if(pair_changed)
        {
            if(dnum_ != pairgran_->dnum_pair())
                error->fix_error(FLERR,this,"Can not change to this pair style with fix wall/gran being active");
        }

        // re-initialize history if contact history was registered by fix with different # hist values
        for(int i=0;i<n_FixMesh_;i++)
            FixMesh_list_[i]->contactHistory()->reset_history(dnum_);

        // check if a fix rigid is registered - important for damp
        fix_rigid_ = pairgran_->fr_pair();

        if (strcmp(update->integrate_style,"respa") == 0)
          nlevels_respa_ = ((Respa *) update->integrate)->nlevels;

        // init for derived classes
        init_granular();
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
    double halfskin = neighbor->skin*0.5;
    int nlocal = atom->nlocal;

    x_ = atom->x;
    radius_ = atom->radius;
    cutneighmax_ = neighbor->cutneighmax;

    // build neighlist for primitive walls
    
    if(rebuildPrimitiveNeighlist_)
      primitiveWall_->buildNeighList(radius_ ? halfskin:(r0_+halfskin),x_,radius_,nlocal);

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
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */

void FixWallGran::post_force_wall(int vflag)
{
  // set pointers and values appropriately
  nlocal_ = atom->nlocal;
  int *mask = atom->mask;
  x_ = atom->x;
  f_ = atom->f;
  radius_ = atom->radius;
  rmass_ = atom->rmass;

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
    double force_old[3],force_wall[3],v_wall[3],bary[3];
    double delta[3],deltan;
    double *c_history = 0;
    MultiVectorContainer<double,3,3> *vMeshC;
    double ***vMesh;
    int nlocal = atom->nlocal;
    int nTriAll, iPart;

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistory *fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // mark all contacts for delettion at this point
      
      if(fix_contact) fix_contact->markAllContacts();

      int *neighborList  = 0, *numNeigh = 0;

      // get neighborList and numNeigh
      FixMesh_list_[iMesh]->meshNeighlist()->getPointers(neighborList, numNeigh);

      vectorZeroize3D(v_wall);
      vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

      atom_type_wall_ = FixMesh_list_[iMesh]->atomTypeWall();

      // moving mesh
      if(vMeshC)
      {
        vMesh = vMeshC->begin();

        // loop owned and ghost triangles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          for(int iCont = 0; iCont < numNeigh[iTri]; iCont++,neighborList++)
          {
            
            int iPart = *neighborList;

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

            deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_ ,x_[iPart],delta,bary);

            if(deltan > skinDistance_) //allow force calculation away from the wall
            {
              
            }
            else
            {
              
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,c_history)) continue;

              for(int i = 0; i < 3; i++)
                v_wall[i] = (bary[0]*vMesh[iTri][0][i] + bary[1]*vMesh[iTri][1][i] + bary[2]*vMesh[iTri][2][i]);

              post_force_eval_contact(iPart,deltan,delta,v_wall,c_history,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }

          }
        }
      }
      // non-moving mesh - do not calculate v_wall, use standard distance function
      else
      {
        // loop owned and ghost particles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          
          for(int iCont = 0; iCont < numNeigh[iTri]; iCont++,neighborList++)
          {
            int iPart = *neighborList;

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);
            deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);

            if(deltan > skinDistance_) //allow force calculation away from the wall
            {
              
            }
            else 
            {
              
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,c_history)) continue;
              post_force_eval_contact(iPart,deltan,delta,v_wall,c_history,iMesh,FixMesh_list_[iMesh],mesh,iTri);
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

  // contact properties
  double force_old[3],force_wall[3];
  double delta[3],deltan,rdist[3];
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

    deltan = primitiveWall_->resolveContact(x_[iPart],radius_?radius_[iPart]:r0_,delta);

    if(deltan > skinDistance_) //allow force calculation away from the wall
    {
      if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
    }
    else
    {
      if(shear_ && shearAxis_ >= 0)
      {
          primitiveWall_->calcRadialDistance(x_[iPart],rdist);
          vectorCross3D(shearAxisVec_,rdist,v_wall);
          
      }
      post_force_eval_contact(iPart,deltan,delta,v_wall,c_history?c_history[iPart]:0,NULL);
    }
  }
}

/* ----------------------------------------------------------------------
   actually calculate force, called for both mesh and primitive
------------------------------------------------------------------------- */

inline void FixWallGran::post_force_eval_contact(int iPart, double deltan, double *delta,
     double *v_wall, double *c_history, int iMesh, FixMeshSurface *fix_mesh, TriMesh *mesh, int iTri)
{

  double delr = (radius_ ? radius_[iPart] : r0_) + deltan;
  double rsqr = delr*delr;
  double dx = -delta[0];
  double dy = -delta[1];
  double dz = -delta[2];
  double mass = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];

  double force_old[3], f_pw[3];

  // if force should be stored - remember old force
  if(store_force_ || stress_flag_)
    vectorCopy3D(f_[iPart],force_old);

  // add to cwl
  if(cwl_ && addflag_)
  {
      double contactPoint[3];
      vectorAdd3D(x_[iPart],delta,contactPoint);
      cwl_->add_wall_1(iMesh,mesh->id(iTri),iPart,contactPoint);
  }

  // deltan > 0 in compute_force
  // but negative in distance algorithm
  compute_force(iPart, -deltan,rsqr,mass, dx, dy, dz, v_wall, c_history, 1.);

  // if force should be stored or evaluated
  if(store_force_ || stress_flag_)
  {
    vectorSubtract3D(f_[iPart],force_old,f_pw);

    if(store_force_)
        vectorAdd3D (wallforce_[iPart], f_pw, wallforce_[iPart]);

    if(stress_flag_ && fix_mesh->trackStress())
    {
        static_cast<FixMeshSurfaceStress*>(fix_mesh)->add_particle_contribution
        (
           iPart,f_pw,delta,iTri,v_wall
        );
    }
  }

  // add heat flux
  if(heattransfer_flag_)
    addHeatFlux(mesh,iPart,deltan,1.);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::is_moving()
{
    int i, flag;
    flag = 0;

    if (is_mesh_wall())
    {
        for(i = 0; i < n_FixMesh_; i++)
            if(FixMesh_list_[i]->mesh()->isMoving())
               flag = 1;
    }
    else flag = shear_;

    return flag;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local()
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts();

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all()
{
    int ncontacts = n_contacts_local();
    MPI_Sum_Scalar(ncontacts,world);
    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int contact_groupbit)
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(contact_groupbit);

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int contact_groupbit)
{
    int ncontacts = n_contacts_local(contact_groupbit);
    MPI_Sum_Scalar(ncontacts,world);
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
