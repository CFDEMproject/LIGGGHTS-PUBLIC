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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifndef LMP_FIX_WALL_GRAN_BASE_H
#define LMP_FIX_WALL_GRAN_BASE_H

#include "fix_wall_gran.h"
#include "contact_interface.h"
#include "compute_pair_gran_local.h"
#include "settings.h"
#include "string.h"
#include "force.h"
#include <stdlib.h>
#include "contact_models.h"

namespace LAMMPS_NS {

using namespace ContactModels;

template<typename ContactModel>
class FixWallGranBase : public FixWallGran {
 public:
  FixWallGranBase(class LAMMPS * lmp, int narg, char **args) : FixWallGran(lmp, narg, args), cmodel(lmp, this) {
    // copy remaining arguments for later use to init contact model
    nfixargs = narg - iarg_;
    fixargs = new char*[nfixargs];

    for(int i = 0; i < nfixargs; i++)
    {
      fixargs[i] = new char[strlen(args[iarg_+i])+1];
      strcpy(fixargs[i], args[iarg_+i]);
    }
  }

  ~FixWallGranBase()
  {
    for(int i = 0; i < nfixargs; i++) delete [] fixargs[i];
    delete [] fixargs;
  }

 protected:
  ContactModel cmodel;
  int nfixargs;
  char ** fixargs;

  void init_granular() {
    cmodel.connectToProperties(force->registry);

    Settings settings(lmp);
    cmodel.registerSettings(settings);
    settings.registerDoubleSetting("temperature", Temp_wall, -1.0);
    bool success = settings.parseArguments(nfixargs, fixargs);

    if(!success) {
      error->fix_error(FLERR, this, settings.error_message.c_str());
    }
  }

  inline void force_update(double * const f, double * const torque,
      const ForceData & forces) {
    for (int coord = 0; coord < 3; coord++) {
      f[coord] += forces.delta_F[coord];
      torque[coord] += forces.delta_torque[coord];
    }
  }

  virtual void compute_force(CollisionData & cdata, double *vwall)
  {
    const int ip = cdata.i;

    double *f = atom->f[ip];
    double *torque = atom->torque[ip];
    double *v = atom->v[ip];
    double *omega = atom->omega[ip];
    double radius = atom->radius[ip];
    double mass = atom->rmass[ip];
    int *type = atom->type;

    if(fix_rigid_ && body_[ip] >= 0)
      mass = masstotal_[body_[ip]];

    const double r = cdata.r;
    const double rinv = 1.0/r;

    const double enx = cdata.delta[0] * rinv;
    const double eny = cdata.delta[1] * rinv;
    const double enz = cdata.delta[2] * rinv;

    // copy collision data to struct (compiler can figure out a better way to
    // interleave these stores with the double calculations above.
    ForceData i_forces;
    ForceData j_forces;
    cdata.v_i = v;
    cdata.v_j = vwall;
    cdata.omega_i = omega;
    cdata.en[0] = enx;
    cdata.en[1] = eny;
    cdata.en[2] = enz;
    cdata.computeflag = computeflag_;
    cdata.shearupdate = shearupdate_;
    cdata.i = ip;
    cdata.radi = radius;
    cdata.touch = NULL;
    cdata.itype = type[ip];
    cdata.jtype = atom_type_wall_;
    cdata.r = r;
    cdata.rinv = rinv;
    cdata.radsum = radius;
    cdata.mi = mass;

    cmodel.collision(cdata, i_forces, j_forces);

    if(computeflag_)
    {
      force_update(f, torque, i_forces);
    }

    if(cwl_ && addflag_)
      cwl_add_wall_2(cdata, i_forces);
  }

  void cwl_add_wall_2(CollisionData & cdata, ForceData & i_forces)
  {
    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0]*cdata.area_ratio;
    const double tor2 = i_forces.delta_torque[1]*cdata.area_ratio;
    const double tor3 = i_forces.delta_torque[2]*cdata.area_ratio;
    cwl_->add_wall_2(cdata.i,fx,fy,fz,tor1,tor2,tor3,cdata.contact_history,cdata.rsq);
  }
};

}

#endif
