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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_FIX_WALL_GRAN_BASE_H
#define LMP_FIX_WALL_GRAN_BASE_H

#include "fix_wall_gran.h"
#include "fix_contact_property_atom_wall.h"
#include "contact_interface.h"
#include "compute_pair_gran_local.h"
#include "tri_mesh.h"
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

  virtual void compute_force(FixWallGran * wg, SurfacesIntersectData & sidata, double *vwall, class TriMesh *mesh = 0,int iTri = 0)
  {
    const int ip = sidata.i;

    double *f = atom->f[ip];
    double *torque = atom->torque[ip];
    double *v = atom->v[ip];
    double *omega = atom->omega[ip];
    double radius = atom->radius[ip];
    double mass = atom->rmass[ip];
    int *type = atom->type;

    if(wg->fix_rigid() && wg->body(ip) >= 0)
      mass = wg->masstotal(wg->body(ip));

    const double r = sidata.r;
    const double rinv = 1.0/r;

    const double enx = sidata.delta[0] * rinv;
    const double eny = sidata.delta[1] * rinv;
    const double enz = sidata.delta[2] * rinv;

    // copy collision data to struct (compiler can figure out a better way to
    // interleave these stores with the double calculations above.
    ForceData i_forces;
    ForceData j_forces;
    sidata.v_i = v;
    sidata.v_j = vwall;
    sidata.omega_i = omega;
    sidata.en[0] = enx;
    sidata.en[1] = eny;
    sidata.en[2] = enz;
    sidata.i = ip;
    sidata.radi = radius;
    sidata.contact_flags = NULL;
    sidata.itype = type[ip];

    sidata.r = r;
    sidata.rinv = rinv;
    sidata.radsum = radius;
    sidata.mi = mass;

    cmodel.surfacesIntersect(sidata, i_forces, j_forces);

    if(sidata.computeflag) {
      force_update(f, torque, i_forces);
    }

    if (wg->store_force_contact()) {
      wg->add_contactforce_wall(ip,i_forces,mesh?mesh->id(iTri):0);
    }

    if(wg->compute_pair_gran_local() && wg->addflag()) {
      wg->cwl_add_wall_2(sidata, i_forces);
    }
  }
};

}

}

#endif
