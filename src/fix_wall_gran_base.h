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
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_FIX_WALL_GRAN_BASE_H
#define LMP_FIX_WALL_GRAN_BASE_H

#include "fix_wall_gran.h"
#include "fix_contact_property_atom_wall.h"
#include "contact_interface.h"
#include "compute_pair_gran_local.h"
#include "settings.h"
#include "string.h"
#include "force.h"
#include <stdlib.h>
#include "contact_models.h"
#include "granular_wall.h"

namespace LIGGGHTS {
using namespace ContactModels;
namespace Walls {

template<typename ContactModel>
class Granular : private Pointers, public IGranularWall {
  ContactModel cmodel;
  FixWallGran * parent;

public:
  Granular(LAMMPS * lmp, FixWallGran * parent) :
    Pointers(lmp),
    cmodel(lmp, parent),
    parent(parent)
  {
  }

  virtual ~Granular() {
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

  virtual void settings(int nargs, char ** args) {
    Settings settings(lmp);
    cmodel.registerSettings(settings);
    bool success = settings.parseArguments(nargs, args);

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

  virtual void compute_force(FixWallGran * wg, CollisionData & cdata, double *vwall)
  {
    const int ip = cdata.i;

    double *f = atom->f[ip];
    double *torque = atom->torque[ip];
    double *v = atom->v[ip];
    double *omega = atom->omega[ip];
    double radius = atom->radius[ip];
    double mass = atom->rmass[ip];
    int *type = atom->type;

    if(wg->fix_rigid() && wg->body(ip) >= 0)
      mass = wg->masstotal(wg->body(ip));

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
    cdata.i = ip;
    cdata.radi = radius;
    cdata.touch = NULL;
    cdata.itype = type[ip];

    cdata.r = r;
    cdata.rinv = rinv;
    cdata.radsum = radius;
    cdata.mi = mass;

    cmodel.collision(cdata, i_forces, j_forces);

    if(cdata.computeflag) {
      force_update(f, torque, i_forces);
    }

    if (wg->store_force_contact()) {
      wg->add_contactforce_wall(ip,i_forces);
    }

    if(wg->compute_pair_gran_local() && wg->addflag()) {
      wg->cwl_add_wall_2(cdata, i_forces);
    }
  }
};

}

}

#endif
