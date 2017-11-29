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
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_FIX_WALL_GRAN_BASE_H
#define LMP_FIX_WALL_GRAN_BASE_H

#include "fix_wall_gran.h"
#include "fix_contact_property_atom_wall.h"
#include "contact_interface.h"
#include "compute_pair_gran_local.h"
#include "tri_mesh.h"
#include "settings.h"
#include <string.h>
#include "force.h"
#include <stdlib.h>
#include "contact_models.h"
#include "granular_wall.h"
#include "fix_calculate_energy_wall.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_const.h"
#endif

namespace LIGGGHTS {
using namespace ContactModels;
namespace Walls {

template<typename ContactModel>
class Granular : private Pointers, public IGranularWall {
  ContactModel cmodel;
  FixWallGran * parent;
  int dissipation_offset_;
  FixCalculateWallEnergy *fix_dissipated_energy_;

public:
  Granular(LAMMPS * lmp, FixWallGran * parent, const int64_t hash) :
    Pointers(lmp),
    cmodel(lmp, parent,true /*is_wall*/, hash),
    parent(parent),
    dissipation_offset_(-1),
    fix_dissipated_energy_(NULL)
  {
  }

  virtual ~Granular()
  {
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
      force->registry.print_all(screen);
      fprintf(screen, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);

      fprintf(logfile, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
      force->registry.print_all(logfile);
      fprintf(logfile, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
    }
#endif
  }

  virtual void settings(int nargs, char ** args, IContactHistorySetup *hsetup) {
    Settings settings(lmp);
    cmodel.registerSettings(settings);
    bool success = settings.parseArguments(nargs, args);
    cmodel.postSettings(hsetup);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== WALL %s SETTINGS ====\n", parent->id);
      settings.print_all(screen);
      fprintf(screen, "==== WALL %s SETTINGS ====\n", parent->id);

      fprintf(logfile, "==== WALL %s SETTINGS ====\n", parent->id);
      settings.print_all(logfile);
      fprintf(logfile, "==== WALL %s SETTINGS ====\n", parent->id);
    }
#endif

    dissipation_offset_ = get_history_offset("dissipation_force");
    fix_dissipated_energy_ = static_cast<FixCalculateWallEnergy*>(modify->find_fix_style("calculate/wall_dissipated_energy", 0));
    if (dissipation_offset_ >= 0 && !fix_dissipated_energy_)
        error->one(FLERR, "Could not find fix calculate/wall_dissipated_energy");

    if(!success) {
      error->fix_error(FLERR, parent, settings.error_message.c_str());
    }
  }

  inline void force_update(double * const f, double * const torque,
      const ForceData & forces) {
    for (int coord = 0; coord < 3; coord++) {
      f[coord] += forces.delta_F[coord];
      torque[coord] += forces.delta_torque[coord];
    }
  }

  virtual bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
  {
    return cmodel.checkSurfaceIntersect(sidata);
  }

  virtual void compute_force(FixWallGran * wg, SurfacesIntersectData & sidata, bool intersectflag,double *vwall, class FixMeshSurface * fix_mesh = 0, int iMesh = 0, class TriMesh *mesh = 0,int iTri = 0)
  {
    const int ip = sidata.i;

    double *x = atom->x[ip];
    double *f = atom->f[ip];
    double *torque = atom->torque[ip];
    double *v = atom->v[ip];
    double *omega = atom->omega[ip];
    double mass = atom->rmass[ip];
    int *type = atom->type;

#ifdef LIGGGHTS_DEBUG
    if(std::isnan(vectorMag3D(x)))
      error->one(FLERR,"x is NaN!");
    if(std::isnan(vectorMag3D(f)))
      error->one(FLERR,"f is NaN!");
    if(std::isnan(vectorMag3D(torque)))
      error->one(FLERR,"torque is NaN!");
    if(std::isnan(vectorMag3D(omega)))
      error->one(FLERR,"omega is NaN!");
#endif

    // copy collision data to struct (compiler can figure out a better way to
    // interleave these stores with the double calculations above.
    ForceData i_forces;
    ForceData j_forces;
    memset(&i_forces, 0, sizeof(ForceData));
    memset(&j_forces, 0, sizeof(ForceData));
    sidata.v_i = v;
    sidata.v_j = vwall;
    sidata.omega_i = omega;

    sidata.r = sidata.radi - sidata.deltan; // sign corrected, because negative value is passed
    sidata.rsq = sidata.r*sidata.r;
    const double rinv = 1.0/sidata.r;
    sidata.rinv = rinv;
    sidata.area_ratio = 1.;

    //store type here as negative in case of primitive wall!
    sidata.j = mesh ? iTri : -wg->atom_type_wall();
    sidata.contact_flags = NULL;
    sidata.itype = type[ip];

    if(wg->fix_rigid() && wg->body(ip) >= 0)
      mass = wg->masstotal(wg->body(ip));

    sidata.meff = mass;
    sidata.mi = mass;

    sidata.computeflag = wg->computeflag();
    sidata.shearupdate = wg->shearupdate();
    sidata.jtype = wg->atom_type_wall();

    double force_old[3]={}, f_pw[3];

    // if force should be stored - remember old force
    if(wg->store_force() || fix_mesh)
        vectorCopy3D(f,force_old);

    // add to cwl
    if(wg->compute_wall_gran_local() && wg->addflag())
    {
      double contactPoint[3];
      vectorSubtract3D(x,sidata.delta,contactPoint);
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(atom->superquadric_flag) {
          vectorCopy3D(sidata.contact_point, contactPoint);
        }
      #endif
      wg->compute_wall_gran_local()->add_wall_1(iMesh,mesh->id(iTri),ip,contactPoint,vwall);
    }

    #ifdef SUPERQUADRIC_ACTIVE_FLAG
        double enx, eny, enz;
        if(atom->superquadric_flag) {
        #ifdef LIGGGHTS_DEBUG
          if(vectorMag3D(sidata.delta) == 0.0)
            error->one(FLERR, "vectorMag3D(sidata.delta) == 0.0");
        #endif
        const double delta_inv = 1.0 / vectorMag3D(sidata.delta);
          #ifdef LIGGGHTS_DEBUG
            if(std::isnan(delta_inv))
              error->one(FLERR, "delta_inv is NaN!");
          #endif
          enx = sidata.delta[0] * delta_inv;
          eny = sidata.delta[1] * delta_inv;
          enz = sidata.delta[2] * delta_inv;
          sidata.radi = cbrt(0.75 * atom->volume[ip] / M_PI);
          #ifdef LIGGGHTS_DEBUG
            if(std::isnan(sidata.radi))
              error->one(FLERR, "delta_inv is NaN!");
          #endif
        } else { // sphere case
          enx = sidata.delta[0] * rinv;
          eny = sidata.delta[1] * rinv;
          enz = sidata.delta[2] * rinv;
        }
    #else // sphere case
        const double enx = sidata.delta[0] * rinv;
        const double eny = sidata.delta[1] * rinv;
        const double enz = sidata.delta[2] * rinv;
    #endif

    sidata.radsum = sidata.radi;

#ifdef CONVEX_ACTIVE_FLAG
    if (!atom->shapetype_flag)
#endif
    {
        sidata.en[0] = enx;
        sidata.en[1] = eny;
        sidata.en[2] = enz;
    }

    double delta[3];
    if (dissipation_offset_ >= 0 && sidata.computeflag && sidata.shearupdate)
    {
        sidata.fix_mesh->triMesh()->get_global_vel(delta);
        vectorScalarMult3D(delta, update->dt);
        // displacement force from the previous time step
        double * const diss_force = &sidata.contact_history[dissipation_offset_];
        // in the scalar product with the current mesh delta
        const double dissipated_energy = vectorDot3D(delta, diss_force)*0.5;
        fix_dissipated_energy_->dissipate_energy(dissipated_energy, true);
        // reset the dissipation force
        diss_force[0] = 0.0;
        diss_force[1] = 0.0;
        diss_force[2] = 0.0;
    }

    if (intersectflag)
    {
       cmodel.surfacesIntersect(sidata, i_forces, j_forces);
       cmodel.endSurfacesIntersect(sidata, mesh, i_forces, j_forces);
       // if there is a surface touch, there will always be a force
       sidata.has_force_update = true;
    }
    // surfacesClose is not supported for convex particles
    else if (!atom->shapetype_flag)
    {
       // apply force update only if selected contact models have requested it
       sidata.has_force_update = false;
       cmodel.surfacesClose(sidata, i_forces, j_forces);
    }

    if (sidata.is_wall && dissipation_offset_ >= 0 && sidata.computeflag && sidata.shearupdate)
    {
        // new displacement force from this time step
        double * const diss_force = &sidata.contact_history[dissipation_offset_];
        // in the scalar product with the mesh displacement from earlier on
        const double dissipated_energy = vectorDot3D(delta, diss_force)*0.5;
        //printf("pdis %e %e\n", update->get_cur_time(), dissipated_energy);
        fix_dissipated_energy_->dissipate_energy(dissipated_energy, false);
    }

    if(sidata.computeflag)
    {
        if (sidata.has_force_update)
            force_update(f, torque, i_forces);
        // summation of f.n to compute a simplistic pressure
        if (wg->store_sum_normal_force())
        {
            double * const iforce = wg->get_sum_normal_force_ptr(ip);
            const double fDotN = vectorDot3D(i_forces.delta_F, sidata.en);
            *iforce += fDotN;
        }
    }

    if (wg->store_force_contact() && 0 == update->ntimestep % wg->store_force_contact_every())
    {
      wg->add_contactforce_wall(ip,i_forces,mesh?mesh->id(iTri):0);
    }

    if (wg->store_force_contact_stress())
      wg->add_contactforce_stress_wall(ip, i_forces, sidata.delta, vwall, mesh?mesh->id(iTri):0);

    if(wg->compute_wall_gran_local() && wg->addflag())
    {
          const double fx = i_forces.delta_F[0];
          const double fy = i_forces.delta_F[1];
          const double fz = i_forces.delta_F[2];
          const double tor1 = i_forces.delta_torque[0]*sidata.area_ratio;
          const double tor2 = i_forces.delta_torque[1]*sidata.area_ratio;
          const double tor3 = i_forces.delta_torque[2]*sidata.area_ratio;
          double normal[3];
          vectorCopy3D(sidata.en, normal);
          vectorNegate3D(normal);
          wg->compute_wall_gran_local()->add_wall_2(sidata.i,fx,fy,fz,tor1,tor2,tor3,sidata.contact_history,sidata.rsq, normal);
    }

    // add heat flux
    
    if(wg->heattransfer_flag())
       wg->addHeatFlux(mesh,ip,sidata.radi,sidata.deltan,1.);

    // if force should be stored or evaluated
    if(sidata.has_force_update && (wg->store_force() || fix_mesh) )
    {
        vectorSubtract3D(f,force_old,f_pw);

        if(wg->store_force())
            vectorAdd3D (wg->fix_wallforce()->array_atom[ip], f_pw, wg->fix_wallforce()->array_atom[ip]);

        if (fix_mesh)
        {
            double delta[3];
            delta[0] = -sidata.delta[0];
            delta[1] = -sidata.delta[1];
            delta[2] = -sidata.delta[2];
            fix_mesh->add_particle_contribution (ip,f_pw,delta,iTri,vwall);
        }
    }
  }

  int get_history_offset(const std::string hname)
  {
    return cmodel.get_history_offset(hname);
  }

  bool contact_match(const std::string mtype, const std::string model)
  {
    return cmodel.contact_match(mtype, model);
  }

};

}

}

#endif
