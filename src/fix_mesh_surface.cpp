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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author, surface velocity rotation only:
   Evan Smuts (U Cape Town)
------------------------------------------------------------------------- */

#include "fix_mesh_surface.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include "error.h"
#include "group.h"
#include "force.h"
#include "bounding_box.h"
#include "input_mesh_tri.h"
#include "fix_contact_history_mesh.h"
#include "fix_contact_property_atom_wall.h"
#include "fix_neighlist_mesh.h"
#include "multi_node_mesh.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "style_mesh_module.h"
#include "mesh_module_stress.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON_V 0.00001

enum{NONE,CONSTANT,EQUAL};

FixMeshSurface::FixMeshSurface(LAMMPS *lmp, int narg, char **arg)
: FixMesh(lmp, narg, arg),
  fix_contact_history_mesh_(0),
  fix_mesh_neighlist_(0),
  fix_meshforce_contact_(NULL),
  fix_meshforce_contact_stress_(NULL),
  fix_mesh_multicontact_data_(NULL),
  velFlag_(false),
  vSurfStrX_(NULL),
  vSurfStrY_(NULL),
  vSurfStrZ_(NULL),
  vSurfVarX_(-1),
  vSurfVarY_(-1),
  vSurfVarZ_(-1),
  vSurfStyleX_(-1),
  vSurfStyleY_(-1),
  vSurfStyleZ_(-1),
  angVelFlag_(false),
  omegaStr_(NULL),
  omegaVar_(-1),
  omegaStyle_(-1),
  n_dump_active_(0),
  curvature_(0.),
  curvature_tolerant_(false),
  extrude_mesh_(false),
  extrusion_length_(0.0),
  extrusion_tri_count_(0),
  extrusion_tri_nodes_(NULL),
  extrusion_created_(false)
{
    // check if type has been read
    if(atom_type_mesh_ == -1)
       error->fix_error(FLERR,this,"expecting keyword 'type'");

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"surface_vel") == 0) {
            if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments");
            iarg_++;
            velFlag_ = true;
            if (strstr(arg[iarg_], "v_") == arg[iarg_])
            {
                int n = strlen(&arg[iarg_][2]) + 1;
                vSurfStrX_ = new char[n];
                strcpy(vSurfStrX_, &arg[iarg_++][2]);
            }
            else
            {
                vSurf_[0] = force->numeric(FLERR,arg[iarg_++]);
                vSurfStyleX_ = CONSTANT;
            }
            if (strstr(arg[iarg_], "v_") == arg[iarg_])
            {
                int n = strlen(&arg[iarg_][2]) + 1;
                vSurfStrY_ = new char[n];
                strcpy(vSurfStrY_, &arg[iarg_++][2]);
            }
            else
            {
                vSurf_[1] = force->numeric(FLERR,arg[iarg_++]);
                vSurfStyleY_ = CONSTANT;
            }
            if (strstr(arg[iarg_], "v_") == arg[iarg_])
            {
                int n = strlen(&arg[iarg_][2]) + 1;
                vSurfStrZ_ = new char[n];
                strcpy(vSurfStrZ_, &arg[iarg_++][2]);
            }
            else
            {
                vSurf_[2] = force->numeric(FLERR,arg[iarg_++]);
                vSurfStyleZ_ = CONSTANT;
            }
            hasargs = true;
        } else if (strcmp(arg[iarg_],"surface_ang_vel") == 0) {
            if (narg < iarg_+11) error->fix_error(FLERR,this,"not enough arguments");
            iarg_++;
            angVelFlag_ = true;
            if(strcmp(arg[iarg_++],"origin"))
                error->fix_error(FLERR,this,"expecting keyword 'origin' after 'rotation'");
            origin_[0] = force->numeric(FLERR,arg[iarg_++]);
            origin_[1] = force->numeric(FLERR,arg[iarg_++]);
            origin_[2] = force->numeric(FLERR,arg[iarg_++]);
            if(strcmp(arg[iarg_++],"axis"))
                error->fix_error(FLERR,this,"expecting keyword 'axis' after definition of 'origin'");
            axis_[0] = force->numeric(FLERR,arg[iarg_++]);
            axis_[1] = force->numeric(FLERR,arg[iarg_++]);
            axis_[2] = force->numeric(FLERR,arg[iarg_++]);
            if(vectorMag3D(axis_) < EPSILON_V)
                error->fix_error(FLERR,this,"'axis' vector must me non-zero");
            if(strcmp(arg[iarg_++],"omega"))
                error->fix_error(FLERR,this,"expecting keyword 'omega' after definition of 'axis'");
            // positive omega give anti-clockwise (CCW) rotation
            if (strstr(arg[iarg_], "v_") == arg[iarg_])
            {
                int n = strlen(&arg[iarg_][2]) + 1;
                omegaStr_ = new char[n];
                strcpy(omegaStr_, &arg[iarg_++][2]);
            }
            else
            {
                omegaSurf_ = force->numeric(FLERR,arg[iarg_++]);
                omegaStyle_ = CONSTANT;
            }
            hasargs = true;
      } else if (strcmp(arg[iarg_],"curvature") == 0) {
          if (narg < iarg_+2)
            error->fix_error(FLERR,this,"not enough arguments for 'curvature'");
          iarg_++;
          curvature_ = force->numeric(FLERR,arg[iarg_++]);
          if(curvature_ <= 0. || curvature_ > 60)
            error->fix_error(FLERR,this,"0° < curvature < 60° required");
          curvature_ = cos(curvature_*M_PI/180.);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"curvature_tolerant") == 0) {
          if (narg < iarg_+2)
            error->fix_error(FLERR,this,"not enough arguments for 'curvature_tolerant'");
          iarg_++;
          if(0 == strcmp(arg[iarg_],"yes"))
            curvature_tolerant_= true;
          else if(0 == strcmp(arg[iarg_],"no"))
            curvature_tolerant_= false;
          else
            error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'curvature_tolerant'");
          iarg_++;
          hasargs = true;
      } else if (strcmp(arg[iarg_], "extrude_planar") == 0) {
          if (narg < iarg_+2)
            error->fix_error(FLERR,this,"not enough arguments for 'extrude_planar'");
          iarg_++;
          extrude_mesh_ = true;
          extrusion_length_ = force->numeric(FLERR,arg[iarg_]);
          iarg_++;
          hasargs = true;
      } else if (strcmp(style,"mesh/surface") == 0) {
          char *errmsg = new char[strlen(arg[iarg_])+50];
          sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
          error->fix_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }

    // get list of available mesh modules
    std::map<std::string, MeshModuleCreator> *module_map = new std::map<std::string, MeshModuleCreator>();
    std::map<std::string, std::vector<std::string> > *module_restrict_map = new std::map<std::string, std::vector<std::string> >();
    #define MESHMODULE_CLASS
    #define MeshModuleStyle(key, Class) \
        (*module_map)[#key] = &meshmodule_creator<Class>;
    #define MeshModuleRestrict(key1, key2) \
        if (module_restrict_map->find(#key1) != module_restrict_map->end()) \
            (*module_restrict_map)[#key1].push_back(#key2); \
        else \
            module_restrict_map->insert( std::pair<std::string, std::vector<std::string> >(#key1, std::vector<std::string>(1,#key2))); \
        if (module_restrict_map->find(#key2) != module_restrict_map->end()) \
            (*module_restrict_map)[#key2].push_back(#key1); \
        else \
            module_restrict_map->insert( std::pair<std::string, std::vector<std::string> >(#key2, std::vector<std::string>(1,#key1)));
    #include "style_mesh_module.h"
    #undef MeshModuleStyle
    #undef MESHMODULE_CLASS

    // parse style name of fix
    std::string style(arg[2]);
    std::istringstream data(style);
    std::string line;
    char separator = '/';
    std::vector<std::string> module_list;

    bool is_planar = false;

    while (std::getline(data, line, separator))
    {
        if (line.compare("mesh") == 0 || line.compare("surface") == 0)
        { /* do nothing */ }
        else if (line.compare("planar") == 0) //support for the mesh/surface/planar style
        {
            if (active_mesh_modules.size() == 0)
                is_planar = true;
            else
                error->fix_error(FLERR,this,"Fix mesh/surface contains the keyword planar and others. The planar keyword is only allowed on its own, i.e. mesh/surface/planar");
        }
        else if (module_map->find(line) != module_map->end()) // mesh modules
        {
            if (is_planar)
                error->fix_error(FLERR,this,"Fix mesh/surface contains the keyword planar and others. The planar keyword is only allowed on its own, i.e. mesh/surface/planar");
            MeshModuleCreator mm_creator = (*module_map)[line];
            active_mesh_modules.insert(std::pair<std::string, MeshModule*>(line, mm_creator(lmp, iarg_, narg, arg, this)));
            mesh_module_order.push_back(line);

        }
        else
        {
            char *errmsg = new char[strlen(arg[2])+50];
            sprintf(errmsg,"unknown keyword in style definition: %s (%s)", arg[2], line.c_str());
            error->fix_error(FLERR,this,errmsg);
            delete []errmsg;
        }
    }

    // check if we have illegal module combinations
    std::map<std::string, MeshModule*>::iterator it;
    for(it = active_mesh_modules.begin(); it != active_mesh_modules.end(); it++)
    {
        // check if the name of this module appears in the restrict map
        std::map<std::string, std::vector<std::string> >::iterator cur_restrict_map = module_restrict_map->find(it->first);
        if (cur_restrict_map != module_restrict_map->end())
        {
            // if yes check the vector
            std::vector<std::string> *cur_restrict = &(cur_restrict_map->second);
            std::map<std::string, MeshModule*>::iterator it2 = it;
            // loop over remaining modlues and check if they are in the vector
            for (std::advance(it2, 1); it2 != active_mesh_modules.end(); it2++)
            {
                std::vector<std::string>::iterator other_module = std::find(cur_restrict->begin(), cur_restrict->end(), it2->first);
                if (other_module != cur_restrict->end()) {
                    char errormsg[200];
                    sprintf(errormsg, "Conflicting mesh modules loaded: %s and %s\n", cur_restrict_map->first.c_str(), other_module->c_str());
                    error->fix_error(FLERR, this, errormsg);
                }
            }
        }
    }

    // check if these modules want to write vectors
    for(it = active_mesh_modules.begin(); it != active_mesh_modules.end(); it++)
        size_vector += it->second->get_num_vector_components();
    if (size_vector > 0)
        vector_flag = 1;

    delete module_map;
    delete module_restrict_map;

    if (iarg_ < narg) {
        char *errmsg = new char[strlen(arg[iarg_])+50];
        sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
        error->fix_error(FLERR,this,errmsg);
        delete []errmsg;
    }
}

/* ---------------------------------------------------------------------- */

FixMeshSurface::~FixMeshSurface()
{
    delete [] vSurfStrX_;
    delete [] vSurfStrY_;
    delete [] vSurfStrZ_;
    delete [] omegaStr_;
    if (extrusion_tri_nodes_)
        delete [] extrusion_tri_nodes_;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::post_create()
{
    FixMesh::post_create();

    if(curvature_ > 0.)
        triMesh()->setCurvature(curvature_);

    if(curvature_tolerant_)
        triMesh()->setCurvatureTolerant(curvature_tolerant_);

    if(velFlag_ && angVelFlag_)
        error->fix_error(FLERR,this,"cannot use 'surface_vel' and 'surface_ang_vel' together");

    if(velFlag_)
        initVel();

    if(angVelFlag_)
        initAngVel();

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->post_create();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::post_create_pre_restart()
{
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->post_create_pre_restart();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::init()
{
    FixMesh::init();

    if (vSurfStrX_)
    {
        vSurfVarX_ = input->variable->find(vSurfStrX_);
        if (vSurfVarX_ < 0)
            error->all(FLERR, "Variable name for fix mesh/surface surfaceVelX does not exist");
        if (input->variable->equalstyle(vSurfVarX_))
            vSurfStyleX_ = EQUAL;
        else
            error->all(FLERR, "Variable for fix mesh/surface surfaceVelX has invalid style");
    }
    if (vSurfStrY_)
    {
        vSurfVarY_ = input->variable->find(vSurfStrY_);
        if (vSurfVarY_ < 0)
            error->all(FLERR, "Variable name for fix mesh/surface surfaceVelY does not exist");
        if (input->variable->equalstyle(vSurfVarY_))
            vSurfStyleY_ = EQUAL;
        else
            error->all(FLERR, "Variable for fix mesh/surface surfaceVelY has invalid style");
    }
    if (vSurfStrZ_)
    {
        vSurfVarZ_ = input->variable->find(vSurfStrZ_);
        if (vSurfVarZ_ < 0)
            error->all(FLERR, "Variable name for fix mesh/surface surfaceVelZ does not exist");
        if (input->variable->equalstyle(vSurfVarZ_))
            vSurfStyleZ_ = EQUAL;
        else
            error->all(FLERR, "Variable for fix mesh/surface surfaceVelZ has invalid style");
    }

    if (omegaStr_)
    {
        omegaVar_ = input->variable->find(omegaStr_);
        if (omegaVar_ < 0)
            error->all(FLERR, "Variable name for fix mesh/surface omega does not exist");
        if (input->variable->equalstyle(omegaVar_))
            omegaStyle_ = EQUAL;
        else
            error->all(FLERR, "Variable for fix mesh/surface omega has invalid style");
    }

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->init();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::setup(int vflag)
{
    FixMesh::setup(vflag);

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::pre_delete(bool unfixflag)
{
    if(unfixflag && n_dump_active_ > 0)
        error->fix_error(FLERR,this,"can not unfix while dump command is active on mesh");

    if(unfixflag && fix_contact_history_mesh_)
        error->fix_error(FLERR,this,"can not unfix while fix wall/gran command is active on mesh; need to unfix fix wall/gran first, then mesh");

    FixMesh::pre_delete(unfixflag);

    // contact tracker and neighlist are created via fix wall/gran
    deleteWallNeighList();
    deleteAllOtherNeighList();
    deleteContactHistory();
}

/* ---------------------------------------------------------------------- */

int FixMeshSurface::setmask()
{
    int mask = FixMesh::setmask();

    if(velFlag_ && (vSurfVarX_ >= 0 || vSurfVarY_ >= 0 || vSurfVarZ_ >= 0))
        mask |= PRE_FORCE;

    if(angVelFlag_ && omegaVar_ >= 0)
        mask |= PRE_FORCE;

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        mask |= active_mesh_modules[*it]->setmask();

    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::setup_pre_force(int vflag)
{
    FixMesh::setup_pre_force(vflag);

    // create neigh list for owned and local elements
    
    if(meshNeighlist())
        meshNeighlist()->initializeNeighlist();

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->setup_pre_force(vflag);

    if (extrude_mesh_ && !extrusion_created_)
    {
        if (update->nsteps > 0)
            error->all(FLERR, "Extrude mesh requires a run 0 after its definition");
        mesh()->extrudePlanarMesh(extrusion_length_, extrusion_tri_nodes_, extrusion_tri_count_);
        extrusion_created_ = true;
        can_create_mesh_ = true;
    }

}

/* ---------------------------------------------------------------------- */

void FixMeshSurface::pre_force(int vflag)
{
    FixMesh::pre_force(vflag);

    if(velFlag_ && (vSurfVarX_ >= 0 || vSurfVarY_ >= 0 || vSurfVarZ_ >= 0))
        setVel();

    if(angVelFlag_ && omegaVar_ >= 0)
        setAngVel();

    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->pre_force(vflag);
}

/* ----------------------------------------------------------------------
   initial integrate for mesh modules
------------------------------------------------------------------------- */

void FixMeshSurface::initial_integrate(int vflag)
{
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->initial_integrate(vflag);
}

/* ----------------------------------------------------------------------
   reverse comm for mesh
------------------------------------------------------------------------- */

void FixMeshSurface::final_integrate()
{
    // Mesh modules final integrate before reverse communication of properties
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->final_integrate_pre_comm();

    // Reverse communication of mesh properties happens in here
    FixMesh::final_integrate();

    // Default final integrate of mesh modules after reverse communication
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->final_integrate();
}

/* ----------------------------------------------------------------------
   End of step call to modules
------------------------------------------------------------------------- */

void FixMeshSurface::end_of_step()
{
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->end_of_step();
}

/* ----------------------------------------------------------------------
   Checks if any modules compute a vector
------------------------------------------------------------------------- */

double FixMeshSurface::compute_vector(int n)
{
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
    {
        const int num = active_mesh_modules[*it]->get_num_vector_components();
        if (n < num)
            return active_mesh_modules[*it]->compute_vector(n);
        else
            n -= num;
    }
    if (n == 0)
        error->fix_error(FLERR, this, "Internal error");

    return 0.0;
}

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMeshSurface::createWallNeighList(int igrp)
{
    if(fix_mesh_neighlist_) return;
    char *neighlist_name = new char[strlen(id)+1+20];
    sprintf(neighlist_name,"wall_neighlist_%s",id);

    const char *fixarg[4];
    fixarg[0]= neighlist_name;
    fixarg[1]= "all";
    fixarg[2]= "neighlist/mesh";
    fixarg[3]= id;
    modify->add_fix(4,const_cast<char**>(fixarg));

    fix_mesh_neighlist_ =
        static_cast<FixNeighlistMesh*>(modify->find_fix_id(neighlist_name));

    // fix added with "all", correct this now
    // groupbit identical to groupbit of this fix
    fix_mesh_neighlist_->igroup = igroup;
    fix_mesh_neighlist_->groupbit = group->bitmask[igroup];

    // neighlist build will use this: merging of wall and mesh groupbit
    fix_mesh_neighlist_->groupbit_wall_mesh = group->bitmask[igroup] | igrp;

    delete []neighlist_name;

    /*
    
    fix_mesh_neighlist_->pre_neighbor();
    fix_mesh_neighlist_->pre_force(0);
    */
}

void FixMeshSurface::deleteWallNeighList()
{
    
    if(fix_mesh_neighlist_) {
      modify->delete_fix(fix_mesh_neighlist_->id,true);
      fix_mesh_neighlist_ = NULL;
    }
}

/* ----------------------------------------------------------------------
   called from fix massflow/mesh out of post_create()
------------------------------------------------------------------------- */

class FixNeighlistMesh* FixMeshSurface::createOtherNeighList(int igrp,const char *nId)
{
    FixNeighlistMesh* neighlist;

    char *neighlist_name = new char[strlen(id)+1+20+strlen(nId)+1];
    sprintf(neighlist_name,"neighlist_%s_%s",nId,id);

    if(modify->find_fix_id(neighlist_name))
        error->fix_error(FLERR,this,"must not use the same mesh for fix massflow/mesh with same group");

    const char *fixarg[5];
    fixarg[0]= neighlist_name;
    fixarg[1]= "all";
    fixarg[2]= "neighlist/mesh";
    fixarg[3]= id;
    fixarg[4]= "other_yes";
    modify->add_fix(5,const_cast<char**>(fixarg));

    neighlist =
        static_cast<FixNeighlistMesh*>(modify->find_fix_id(neighlist_name));

    // fix added with "all", correct this now
    neighlist->igroup = igrp;
    neighlist->groupbit = group->bitmask[igrp];

    list_other_neighlist_[neighlist_name] = neighlist;

    delete []neighlist_name;

    return neighlist;
}

void FixMeshSurface::deleteOtherNeighList(const char *nId)
{
    char *neighlist_name = new char[strlen(id)+1+20+strlen(nId)+1];
    sprintf(neighlist_name,"neighlist_%s_%s",nId,id);

    std::map<std::string,FixNeighlistMesh*>::iterator it = list_other_neighlist_.find(neighlist_name);
    if (it != list_other_neighlist_.end()) {
      modify->delete_fix(it->second->id,true);
      list_other_neighlist_.erase(it);
    }
    
    delete []neighlist_name;
}

void FixMeshSurface::deleteAllOtherNeighList()
{
    for ( std::map<std::string,FixNeighlistMesh*>::iterator it = list_other_neighlist_.begin();
          it != list_other_neighlist_.end();
          ++ it) {
      modify->delete_fix(it->second->id,true);
    }
    list_other_neighlist_.clear();
}

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMeshSurface::createContactHistory(int dnum)
{
    if(fix_contact_history_mesh_) return;

    // create a contact tracker for the mesh
    char *contacthist_name = new char[strlen(id)+1+8];
    char *my_id = new char[strlen(id)+1];
    sprintf(contacthist_name,"tracker_%s",id);
    sprintf(my_id,"%s",id);

    char dnum_char[10];
    sprintf(dnum_char,"%d",dnum);

    const char *fixarg[5];
    fixarg[0] = contacthist_name;
    fixarg[1] = "all";
    fixarg[2] = "contacthistory/mesh";
    fixarg[3] = dnum_char;
    fixarg[4] = my_id;

    modify->add_fix(5,const_cast<char**>(fixarg));

    fix_contact_history_mesh_ = static_cast<FixContactHistoryMesh*>(modify->find_fix_id(contacthist_name));

    delete []contacthist_name;
    delete []my_id;
}

void FixMeshSurface::deleteContactHistory()
{
    // contact tracker and neighlist are created via fix wall/gran
    
    if(fix_contact_history_mesh_) {
      modify->delete_fix(fix_contact_history_mesh_->id,true);
      fix_contact_history_mesh_ = NULL;
    }
}

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMeshSurface::createMeshforceContact()
{
    if(fix_meshforce_contact_) return;

    const char *fixarg[19];
    char fixid[200],propertyid[200];
    sprintf(fixid,"contactforces_%s",id);
    sprintf(propertyid,"contactforces_%s",id);
    fixarg[0]=fixid;
    fixarg[1]="all";
    fixarg[2]="contactproperty/atom/wall";
    fixarg[3]=propertyid;
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
    fixarg[17]="mesh";
    fixarg[18]=this->id;
    modify->add_fix(19,const_cast<char**>(fixarg));
    fix_meshforce_contact_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
}

void FixMeshSurface::createMeshforceContactStress()
{
    if(fix_meshforce_contact_stress_) return;

    const char *fixarg[25];
    char fixid[200],propertyid[200];
    sprintf(fixid,"contactforces_stress_%s",id);
    sprintf(propertyid,"contactforces_stress_%s",id);
    fixarg[0]=fixid;
    fixarg[1]="all";
    fixarg[2]="contactproperty/atom/wall";
    fixarg[3]=propertyid;
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
    fixarg[23]="mesh";
    fixarg[24]=this->id;
    modify->add_fix(25,const_cast<char**>(fixarg));
    fix_meshforce_contact_stress_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
}

void FixMeshSurface::createMulticontactData()
{
    if(fix_mesh_multicontact_data_) return;

    // create a new per contact property which will contain the data for the computation according to Brodu et. al. 2016
    // surfPosIJ will contain the position of the contact surface ij, realtive to position i
    // normalForce will contain the normal component of the contact force
    const char *fixarg[17];
    char fixid[200];
    sprintf(fixid,"multicontactData_%s",id);
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
    fixarg[13]="mesh";
    fixarg[14]=id;
    fixarg[15]="reset";
    fixarg[16]="no";
    modify->add_fix(17,const_cast<char**>(fixarg));
    fix_mesh_multicontact_data_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
}

void FixMeshSurface::deleteMeshforceContact()
{
    // contact tracker and neighlist are created via fix wall/gran
    
    if(fix_meshforce_contact_) {
      modify->delete_fix(fix_meshforce_contact_->id,true);
      fix_meshforce_contact_ = NULL;
    }
}

void FixMeshSurface::deleteMeshforceContactStress()
{
    // contact tracker and neighlist are created via fix wall/gran
    
    if(fix_meshforce_contact_stress_) {
      modify->delete_fix(fix_meshforce_contact_stress_->id,true);
      fix_meshforce_contact_stress_ = NULL;
    }
}

void FixMeshSurface::deleteMeshMulticontactData()
{
    // contact tracker and neighlist are created via fix wall/gran

    if(fix_mesh_multicontact_data_) {
      modify->delete_fix(fix_mesh_multicontact_data_->id,true);
      fix_mesh_multicontact_data_ = NULL;
    }
}

/* ----------------------------------------------------------------------
   sets mesh velocity for conveyor model
------------------------------------------------------------------------- */

void FixMeshSurface::initVel()
{
    if (vSurfStyleX_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[0] = input->variable->compute_equal(vSurfVarX_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    if (vSurfStyleY_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[1] = input->variable->compute_equal(vSurfVarY_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    if (vSurfStyleZ_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[2] = input->variable->compute_equal(vSurfVarZ_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    double conv_vel[3];
    vectorCopy3D(vSurf_,conv_vel);

    // register new mesh property v
    mesh()->prop().addGlobalProperty< VectorContainer<double,3> >("v","comm_exchange_borders","frame_invariant","restart_no");
    mesh()->prop().setGlobalProperty< VectorContainer<double,3> >("v",conv_vel);

    // register new element property v - error if exists already via moving mesh
    mesh()->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_exchange_borders","frame_invariant","restart_no");

    setVel();
}

void FixMeshSurface::setVel()
{
    if (vSurfStyleX_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[0] = input->variable->compute_equal(vSurfVarX_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    if (vSurfStyleY_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[1] = input->variable->compute_equal(vSurfVarY_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    if (vSurfStyleZ_ == EQUAL)
    {
        modify->clearstep_compute();
        vSurf_[2] = input->variable->compute_equal(vSurfVarZ_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    double conv_vel[3];
    int size, nVec;
    double scp, tmp[3], facenormal[3], ***v_node;
    vectorCopy3D(vSurf_,conv_vel);
    double conv_vSurf_mag = vectorMag3D(conv_vel);

    size = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
    nVec = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->nVec();

    v_node = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();

    // set mesh velocity
    TriMesh *trimesh = triMesh();
    
    for (int i = 0; i < size; i++)
    {
        trimesh->surfaceNorm(i,facenormal);
        scp = vectorDot3D(conv_vel,facenormal);
        vectorScalarMult3D(facenormal,scp,tmp);
        for(int j = 0; j < nVec; j++)
        {
            vectorSubtract3D(conv_vel,tmp,v_node[i][j]);
            if(vectorMag3D(v_node[i][j]) > 0.)
            {
                vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
                vectorScalarMult3D(v_node[i][j],conv_vSurf_mag);
                
            }
        }
    }
}

/* ----------------------------------------------------------------------
   sets mesh velocity for rotation model
------------------------------------------------------------------------- */

void FixMeshSurface::initAngVel()
{
    // register new element property v - error if exists already via moving mesh
    mesh()->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_exchange_borders","frame_invariant","restart_no");

    setAngVel();
}

void FixMeshSurface::setAngVel()
{
    if (omegaStyle_ == EQUAL)
    {
        modify->clearstep_compute();
        omegaSurf_ = input->variable->compute_equal(omegaVar_);
        modify->addstep_compute(update->ntimestep + 1);
    }
    int size, nVec;
    double rot_origin[3],rot_axis[3],rot_omega;
    vectorCopy3D(origin_,rot_origin);
    vectorCopy3D(axis_,rot_axis);
    rot_omega = omegaSurf_;
    double tmp[3], scp, unitAxis[3], tangComp[3], Utang[3], surfaceV[3];
    double node[3],facenormal[3], magAxis, magUtang, ***v_node;

    size = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
    nVec = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->nVec();

    v_node = mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();

    // calculate unit vector of rotation axis
    magAxis = vectorMag3D(rot_axis);
    vectorScalarDiv3D(rot_axis,magAxis,unitAxis);

    TriMesh *trimesh = triMesh();
    for (int i = 0; i < size; i++)
    {
        trimesh->surfaceNorm(i,facenormal);
        // number of nodes per face (3)
        for(int j = 0; j < nVec; j++)
        {
            // position of node - origin of rotation (to get lever arm)
            trimesh->node(i,j,node);
            vectorSubtract3D(node,rot_origin,surfaceV);

            // lever arm X rotational axis = tangential component
            vectorCross3D(surfaceV,unitAxis,tangComp);

            // multiplying by omega scales the tangential component to give tangential velocity
            vectorScalarMult3D(tangComp,-rot_omega,Utang);

            // EPSILON is wall velocity, not rotational omega
            //if(vectorMag3D(Utang) < EPSILON_V)
            //    error->fix_error(FLERR,this,"Rotation velocity too low");

            scp = vectorDot3D(Utang, facenormal);
            vectorScalarMult3D(facenormal,scp,tmp);

            // removes components normal to wall
            vectorSubtract3D(Utang,tmp,v_node[i][j]);
            magUtang = vectorMag3D(Utang);

            if(vectorMag3D(v_node[i][j]) > 0.)
            {
               vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
               vectorScalarMult3D(v_node[i][j],magUtang);
            }
        }
    }
}

/* ----------------------------------------------------------------------
   Returns a pointer to a specific mesh module if it is loaded
------------------------------------------------------------------------- */

MeshModule* FixMeshSurface::get_module(std::string name)
{
    std::map<std::string, MeshModule*>::iterator it = active_mesh_modules.find(name);
    if (it != active_mesh_modules.end())
        return it->second;
    else
        return NULL;
}

/* ----------------------------------------------------------------------
   Creates a new mesh module
------------------------------------------------------------------------- */

template <typename T>
MeshModule *FixMeshSurface::meshmodule_creator(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh)
{
  return new T(lmp, iarg_, narg, arg, fix_mesh);
}

void FixMeshSurface::add_particle_contribution(int ip, double *frc, double *delta, int iTri, double *v_wall)
{
    std::vector<std::string>::iterator it;
    for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        active_mesh_modules[*it]->add_particle_contribution(ip, frc, delta, iTri, v_wall);
}

/* ----------------------------------------------------------------------
   Returns if the mesh tracks the stress.
   This function should not used within a performance critical code.
   Instead, save the value during an init/setup function.
------------------------------------------------------------------------- */

bool FixMeshSurface::trackStress()
{
    MeshModuleStress *mms = dynamic_cast<MeshModuleStress*>(get_module("stress"));
    if (mms)
      return mms->trackStress();
    else
      return false;
}

int FixMeshSurface::modify_param(int narg, char **arg)
{
    int ret = 0;
    std::string module(arg[0]);
    std::size_t slash = module.find_first_of('/');

    if (slash != std::string::npos)
    {
        MeshModule *mm_modify = get_module(module.substr(0,slash));
        if (!mm_modify)
        {
            std::string errmsg = "Could not find appropriate mesh module \"";
            errmsg.append(module.substr(0,slash));
            errmsg.append("\" in modify_param");
            error->fix_error(FLERR, this, errmsg.c_str());
        }
        ret = mm_modify->modify_param(narg, arg);
        return ret;
    }
    else
    {
        // Basic divide and conquer algorithm for modify_param in mesh modules. Note that this can cause issues if
        // two mesh modules use the same keyword.
        std::vector<std::string>::iterator it;
        for(it = mesh_module_order.begin(); it != mesh_module_order.end(); it++)
        {
            ret = active_mesh_modules[*it]->modify_param(narg, arg);
            if (ret > 0)
            {
                std::string warnmsg = "Using deprecated modify_param for attribute \"";
                warnmsg.append(module);
                warnmsg.append("\". Consider using ");
                warnmsg.append(*it);
                warnmsg.append("/");
                warnmsg.append(module);
                warnmsg.append(" in your input file to explicitly specify the mesh module the argument is used for.\n");
                error->warning(FLERR, warnmsg.c_str());
                return ret;
            }
        }
    }
    return ret;
}
