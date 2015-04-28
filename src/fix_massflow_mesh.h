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

#ifdef FIX_CLASS

FixStyle(massflow/mesh,FixMassflowMesh)

#else

#ifndef LMP_FIX_MASSFLOW_MESH_H
#define LMP_FIX_MASSFLOW_MESH_H

#include "fix.h"
#include <vector>

using namespace std;

namespace LAMMPS_NS {

class FixMassflowMesh : public Fix {

 public:

  FixMassflowMesh(class LAMMPS *lmp, int narg, char ** arg);
  ~FixMassflowMesh();

  void post_create();
  void pre_delete(bool unfixflag);

  void init();
  void setup(int vflag);
  int setmask();

  void post_integrate();
  void pre_exchange();

  void write_restart(FILE *fp);
  void restart(char *buf);

  double compute_vector(int index);

 protected:

  // in case particles counted should be deleted or transferred
  bool delete_atoms_;
  vector<int> atom_tags_delete_;
  double mass_deleted_;
  double nparticles_deleted_;

  // true if any given particle is
  // counted only once
  bool once_;

  int  iarg_;
  class FixPropertyAtom* fix_orientation_;

 protected:
  class FixPropertyAtom *fix_counter_;

 private:
  class FixMeshSurface *fix_mesh_;
  char fixid_[200];
  class FixNeighlistMesh *fix_neighlist_;
  double nvec_[3];
  double pref_[3];
  double sidevec_[3];

  bool   havePointAtOutlet_;
  bool   insideOut_;
  double pointAtOutlet_[3];

  // mass and particles which was counted
  double mass_;
  int nparticles_;

  // additional property to sum
  class FixPropertyAtom *fix_property_;
  double property_sum_;

  // data write
  bool screenflag_;
  FILE *fp_;
  bool writeTime_; //switch to write time to the outfile

  // data for particle and mass flow calculation
  double mass_last_;
  int nparticles_last_;
  double t_count_, delta_t_;
  bool reset_t_count_;

  class FixMultisphere* fix_ms_;
  class MultisphereParallel *ms_;
  class ScalarContainer<int> *ms_counter_;

}; //end class

}
#endif
#endif
