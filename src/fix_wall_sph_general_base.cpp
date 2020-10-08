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
    Andreas Eitzlmayr (TU Graz)

    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_sph_general_base.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "fix_neighlist_mesh.h"
#include "tri_mesh.h"
#include "mpi_liggghts.h"
#include "pair_sph.h"
#include <vector>
#include "fix_sph_velgrad.h"
#include "fix_sph_pressure.h"
//#include "fix_sph_density_sumconti.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LIGGGHTS::ContactModels;

/* ---------------------------------------------------------------------- */

FixWallSphGeneralBase::FixWallSphGeneralBase(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg)
{
  fix_wallContact_ = NULL;
  fix_wallForce_ = NULL;
  fix_pressure_ = NULL;
  fix_velgrad_ = NULL;
  fix_dvdx_ = NULL;
  fix_dvdy_ = NULL;
  fix_dvdz_ = NULL;

  fixName = arg[0];
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::post_create()
{
    FixWallGran::post_create();

    firstStep = 1;

    const char *fixarg[15];
    char * name = new char[12+strlen(fixName)+1];
    strcpy(name,"wallContact_");
    strcat(name,fixName);
    fixarg[0]=name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=name;
    fixarg[4]="vector";
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="yes";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fixarg[11]="0.";
    fixarg[12]="0.";
    fixarg[13]="0.";
    fixarg[14]="0.";
    fix_wallContact_ = modify->add_fix_property_atom(15,const_cast<char**>(fixarg),style);
    delete [] name;

    pairsph_ = (PairSph*)force->pair;
    pairStyle = pairsph_->returnPairStyle();
    viscosity = pairsph_->returnViscosity();
    kernel_id = pairsph_->sph_kernel_id();

    densityStyle = 0;
    rho0 = -1;

    char* fixID;
    int ifix_viscosity = -1;

    for (int ifix = 0; ifix < modify->nfix; ifix++)
    {
      if (strcmp("sph/density/continuity",modify->fix[ifix]->style) == 0) densityStyle = 1;
      if (strcmp("sph/density/summation",modify->fix[ifix]->style) == 0) densityStyle = 2;
/*      if (strcmp("sph/density/sumconti",modify->fix[ifix]->style) == 0) {
        fix_density_sumconti_ = static_cast<FixSphDensitySumconti*>(modify->fix[ifix]);
        densityStyle = 3;
      }*/
      if (strcmp("sph/pressure",modify->fix[ifix]->style) == 0) {
        fix_pressure_ = static_cast<FixSPHPressure*>(modify->fix[ifix]);
        rho0 = fix_pressure_->return_rho0();
      }
      if (strcmp("sph/velgrad",modify->fix[ifix]->style) == 0)
      {
        fix_velgrad_ = static_cast<FixSphVelgrad*>(modify->fix[ifix]);
        fix_dvdx_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdx","property/atom","vector",0,0,"FixWallSphGeneralBase",false));
        fix_dvdy_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdy","property/atom","vector",0,0,"FixWallSphGeneralBase",false));
        fix_dvdz_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dvdz","property/atom","vector",0,0,"FixWallSphGeneralBase",false));
        //dvdx_ = fix_dvdx_->array_atom;
        //dvdy_ = fix_dvdx_->array_atom;
        //dvdz_ = fix_dvdx_->array_atom;
      }

      fixID = modify->fix[ifix]->id;
      if (strcmp("viscosity",fixID) == 0) ifix_viscosity = ifix;
    }

    if (ifix_viscosity == -1) modelStyle = 1; // Newtonian
    else {  // non-Newtonian
      fix_visc_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("viscosity","property/atom","scalar",0,0,"FixWallGSphGeneralGap",false));
      modelStyle = 2;
    }

    if (rho0 == -1) error->fix_error(FLERR,this,"Requires a fix/sph/pressure also.");
    if (rho0 == 0) error->fix_error(FLERR,this,"Use style 'Tait' or 'relativ' in fix/sph/pressure and define parameter rho0.");

    char* fidforce = new char[30];
    strcpy(fidforce,"F_");
    strcat(fidforce,fixName);

    const char *fixargF[11];
    fixargF[0]=fidforce;
    fixargF[1]="all";
    fixargF[2]="property/atom";
    fixargF[3]=fidforce;
    fixargF[4]="vector";
    fixargF[5]="yes";
    fixargF[6]="yes";
    fixargF[7]="yes";
    fixargF[8]="0.";
    fixargF[9]="0.";
    fixargF[10]="0.";
    fix_wallForce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixargF),style);
    delete [] fidforce;

    fix_wallForce_->global_freq = 1;

    //wallContact_ = fix_wallContact_->array_atom;
    //wallForce_ = fix_wallForce_->array_atom;
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::pre_delete(bool unfixflag)
{
    FixWallGran::pre_delete(unfixflag);

    if(unfixflag && fix_wallContact_)
        modify->delete_fix(fix_wallContact_->id);

    if(unfixflag && fix_wallForce_)
        modify->delete_fix(fix_wallForce_->id);

    if(unfixflag && fix_pressure_)
        modify->delete_fix(fix_pressure_->id);

    if(unfixflag && fix_velgrad_)
        modify->delete_fix(fix_velgrad_->id);

    if(unfixflag && fix_dvdx_)
        modify->delete_fix(fix_dvdx_->id);

    if(unfixflag && fix_dvdy_)
        modify->delete_fix(fix_dvdy_->id);

    if(unfixflag && fix_dvdz_)
        modify->delete_fix(fix_dvdz_->id);

    if(unfixflag && fix_visc_)
        modify->delete_fix(fix_visc_->id);

//    if(unfixflag && fix_density_sumconti_)
//        modify->delete_fix(fix_density_sumconti_->id);
}

/* ---------------------------------------------------------------------- */

FixWallSphGeneralBase::~FixWallSphGeneralBase()
{

}

/* ---------------------------------------------------------------------- */

int FixWallSphGeneralBase::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::init()
{

}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::post_integrate()
{
  int i,nlocal = atom->nlocal;
  double **v = atom->v;
  double delx,dely,delz,rsq,r,mass,dtv;

/*  if (densityStyle == 3) {
    if (fix_density_sumconti_->returnSwitch() == 1) {
      densityStyle = 1; // switch to the continuity equation
    }
  }*/

  if ((densityStyle == 2)) // || (densityStyle == 3)) // with sph/density/summation
  {
    dtv = update->dt;

    for (i = 0; i < nlocal; i++) {
      if (wallContact_[i][0] != 0)
      {
        // wallContact is from last timestep; do update:
        delx = - wallContact_[i][1] + (v[i][0] - wallContact_[i][4])*dtv;
        dely = - wallContact_[i][2] + (v[i][1] - wallContact_[i][5])*dtv;
        delz = - wallContact_[i][3] + (v[i][2] - wallContact_[i][6])*dtv;
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);

        mass = rmass_ ? rmass_[i] : atom->mass[atom->type[i]];

        compute_density(i,r,mass);

        // reset wallContact_
        // wallContact_[i][0] = 0;
      }
    }
  } else {  // with sph/density/continuity
    /*for (i = 0; i < nlocal; i++) {
      // reset wallContact_
      wallContact_[i][0] = 0;
    }*/
  }
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::compute_density(int ip,double r,double mass)
{

}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::pre_force(int vflag)
{
  FixWallGran::pre_force(vflag);

 // Add wall contributions to velocity gradients

  int i,nlocal = atom->nlocal;
  double **v = atom->v;
  double delx,dely,delz,mass,dtv,v_wall[3];

  wallContact_ = fix_wallContact_->array_atom;

  if (firstStep == 1)
  {
    firstStep = 0; // do not calculate in the first time step (wall contact detection is later in post_force)
  }
  else
  {
    if (fix_velgrad_)
    {
      if (fix_velgrad_->return_velgrad_flag() == 1)
      {
        dtv = update->dt;

        for (i = 0; i < nlocal; i++)
        {
          if (wallContact_[i][0] != 0)
          {
            // wallContact is from last timestep; do update:
            v_wall[0] = wallContact_[i][4];
            v_wall[1] = wallContact_[i][5];
            v_wall[2] = wallContact_[i][6];

            delx = - wallContact_[i][1] + (v[i][0] - v_wall[0])*dtv;
            dely = - wallContact_[i][2] + (v[i][1] - v_wall[1])*dtv;
            delz = - wallContact_[i][3] + (v[i][2] - v_wall[2])*dtv;

            mass = rmass_ ? rmass_[i] : atom->mass[atom->type[i]];

            dvdx_ = fix_dvdx_->array_atom;
            dvdy_ = fix_dvdy_->array_atom;
            dvdz_ = fix_dvdz_->array_atom;

            compute_velgrad(i,delx,dely,delz,mass,v_wall);

            // reset wallContact_
            wallContact_[i][0] = 0;
          }
        }
      } else {
        for (i = 0; i < nlocal; i++) {
          // reset wallContact_
          wallContact_[i][0] = 0;
        }
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        // reset wallContact_
        wallContact_[i][0] = 0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralBase::compute_velgrad(int ip,double delx,double dely,double delz,double mass,double *vwall)
{

}

/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallSphGeneralBase::post_force_mesh(int vflag)
{
    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    double *c_history = 0;
    MultiVectorContainer<double,3,3> *vMeshC;
    double ***vMesh;
    int nlocal = atom->nlocal;
    int nTriAll, barysign = -1;
    double rSphere;

    wallContact_ = fix_wallContact_->array_atom;
    wallForce_ = fix_wallForce_->array_atom;

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // mark all contacts for delettion at this point

      if(fix_contact) fix_contact->markAllContacts();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = FixMesh_list_[iMesh]->meshNeighlist();

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
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();
          for(int iCont = 0; iCont < numneigh; iCont++)
          {

            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            //int idTri = mesh->id(iTri);

            rSphere = 0;
            deltan = mesh->resolveTriSphereContactBary(iPart,iTri,rSphere,x_[iPart],delta,bary,barysign);

            if (wallContact_[iPart][0] == 0)
            {
                wallContact_[iPart][0] = deltan;
                wallContact_[iPart][1] = delta[0];
                wallContact_[iPart][2] = delta[1];
                wallContact_[iPart][3] = delta[2];
                wallContact_[iPart][4] = (bary[0]*vMesh[iTri][0][0] + bary[1]*vMesh[iTri][1][0] + bary[2]*vMesh[iTri][2][0]);
                wallContact_[iPart][5] = (bary[0]*vMesh[iTri][0][1] + bary[1]*vMesh[iTri][1][1] + bary[2]*vMesh[iTri][2][1]);
                wallContact_[iPart][6] = (bary[0]*vMesh[iTri][0][2] + bary[1]*vMesh[iTri][1][2] + bary[2]*vMesh[iTri][2][2]);
            }
            else if (deltan < wallContact_[iPart][0])
            {
                wallContact_[iPart][0] = deltan;
                wallContact_[iPart][1] = delta[0];
                wallContact_[iPart][2] = delta[1];
                wallContact_[iPart][3] = delta[2];
                wallContact_[iPart][4] = (bary[0]*vMesh[iTri][0][0] + bary[1]*vMesh[iTri][1][0] + bary[2]*vMesh[iTri][2][0]);
                wallContact_[iPart][5] = (bary[0]*vMesh[iTri][0][1] + bary[1]*vMesh[iTri][1][1] + bary[2]*vMesh[iTri][2][1]);
                wallContact_[iPart][6] = (bary[0]*vMesh[iTri][0][2] + bary[1]*vMesh[iTri][1][2] + bary[2]*vMesh[iTri][2][2]);
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
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();
          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            //int idTri = mesh->id(iTri);
            rSphere = 0;
            deltan = mesh->resolveTriSphereContact(iPart,iTri,rSphere,x_[iPart],delta,barysign);

            if (wallContact_[iPart][0] == 0)
            {
                wallContact_[iPart][0] = deltan;
                wallContact_[iPart][1] = delta[0];
                wallContact_[iPart][2] = delta[1];
                wallContact_[iPart][3] = delta[2];
                wallContact_[iPart][4] = 0;
                wallContact_[iPart][5] = 0;
                wallContact_[iPart][6] = 0;
            }
            else if (deltan < wallContact_[iPart][0])
            {
                wallContact_[iPart][0] = deltan;
                wallContact_[iPart][1] = delta[0];
                wallContact_[iPart][2] = delta[1];
                wallContact_[iPart][3] = delta[2];
                wallContact_[iPart][4] = 0;
                wallContact_[iPart][5] = 0;
                wallContact_[iPart][6] = 0;
            }
          }
        }
      }

      // clean-up contacts
      if(fix_contact) fix_contact->cleanUpContacts();
    }

    for(int i = 0; i < nlocal; i++)
    {
        deltan = wallContact_[i][0];
        if (deltan!=0)
        {
          SurfacesIntersectData sidata;
          sidata.delta[0] = - wallContact_[i][1] / deltan; // normalize delta
          sidata.delta[1] = - wallContact_[i][2] / deltan;
          sidata.delta[2] = - wallContact_[i][3] / deltan;
          sidata.en[0] = sidata.delta[0];
          sidata.en[1] = sidata.delta[1];
          sidata.en[2] = sidata.delta[2];
          sidata.i = i;
          sidata.area_ratio = 1.0;
          sidata.deltan = deltan;
          sidata.contact_history = c_history;
          v_wall[0] = wallContact_[i][4];
          v_wall[1] = wallContact_[i][5];
          v_wall[2] = wallContact_[i][6];

          sidata.r = fabs(deltan);
          sidata.mi = rmass_ ? rmass_[i] : atom->mass[atom->type[i]];
          compute_force(sidata, v_wall);
        } else {
          // reset wall force for particles remote from wall
          wallForce_[i][0] = 0;
          wallForce_[i][1] = 0;
          wallForce_[i][2] = 0;
        }
    }
}

/* ----------------------------------------------------------------------
   post_force for primitive wall
------------------------------------------------------------------------- */

void FixWallSphGeneralBase::post_force_primitive(int vflag)
{
  error->fix_error(FLERR,this,"Primitive wall not yet implemented for sph");
}
