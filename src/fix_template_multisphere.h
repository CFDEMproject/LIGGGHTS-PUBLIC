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

FixStyle(particletemplate/multisphere,FixTemplateMultisphere)

#else

#ifndef LMP_FIX_TEMPLATE_MULTISPHERE_H
#define LMP_FIX_TEMPLATE_MULTISPHERE_H

#include "fix_template_multiplespheres.h"
#include "vector_liggghts.h"

namespace LAMMPS_NS {

class FixTemplateMultisphere : public FixTemplateMultiplespheres {
 public:
  FixTemplateMultisphere(class LAMMPS *, int, char **);
  ~FixTemplateMultisphere();
  virtual void post_create();

  void init();

  // called at single insertion
  virtual void randomize_single();

  // called at multi insertion
  void init_ptilist(int);
  void delete_ptilist();
  void randomize_ptilist(int ,int ,int);

  void finalize_insertion();

  int type() {return type_;}

  void fflag(bool *ff)
  { vectorCopy3D(fflag_,ff); }

  void tflag(bool *tf)
  { vectorCopy3D(tflag_,tf); }

 protected:

  void calc_volumeweight();
  void calc_inertia();
  void calc_eigensystem();
  void calc_displace_xcm_x_body();
  void print_info();

  // type of clump
  int type_;

  // flags
  bool mass_set_, moi_set_;
  int use_density_;

  // inertia of clump
  double moi_[3][3];      // 3x3 inertia tensor
  double inertia_[3]; // 3 principal components of inertia
  double ex_space_[3],ey_space_[3],ez_space_[3]; // principal axes in global coords

  // force and torque flags; DOFs are inactive if set false
  bool fflag_[3],tflag_[3];

  // sphere coordinates in ex,ey,ez frame
  double **displace_;

  // vector from center of mass (which is 0 0 0) to x_bound in body coordinates
  double xcm_to_xb_body_[3];

  // volume weight of each sphere
  // used for volume fraction calculation
  // 1 for spherical or non-overlapping multisphere
  // < 1 for overlapping multisphere
  double *volumeweight_;
};

}

#endif
#endif
