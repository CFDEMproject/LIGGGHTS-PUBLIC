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

FixStyle(particletemplate/multiplespheres,FixTemplateMultiplespheres)

#else

#ifndef LMP_FIX_TEMPLATE_MULTIPLESPHERES_H
#define LMP_FIX_TEMPLATE_MULTIPLESPHERES_H

#include "fix.h"
#include "fix_template_sphere.h"

namespace LAMMPS_NS {

class FixTemplateMultiplespheres : public FixTemplateSphere {
 public:

  FixTemplateMultiplespheres(class LAMMPS *, int, char **);
  virtual ~FixTemplateMultiplespheres();

  virtual void post_create();
  double max_r_bound();
  double max_rad();
  double min_rad();
  int maxtype();
  int mintype();
  int number_spheres();

  // single insertion
  virtual void randomize_single();

  // multi insertion
  virtual void init_ptilist(int);
  void randomize_ptilist(int ,int ,int);

  virtual void finalize_insertion() {}

  inline bool all_overlap_none()
  { return no_overlap; }

  inline bool all_overlap_atleast_one_slightly()
  { return overlap_slightly; }

 protected:

  // template calculations
  virtual void calc_bounding_sphere();
  virtual void calc_center_of_mass();
  virtual void check_overlap();
  virtual void print_info();

  // sqr distance from x_sphere[j] to xtest
  double dist_sqr(int j,double *xtest);

  // generate random point in bbox
  void generate_xtry(double *xtry);

  // number of spheres in template
  int nspheres;

  // coords of each sphere with respect to center of mass
  double **x_sphere;

  // radius of each sphere
  double *r_sphere;

  // scale factor if read from a file
  double scale_fact;

  // atom type might be variable if read from file
  int *atom_type_sphere;

  // bounding box
  double x_min[3], x_max[3];

  // bounding sphere - radius and coordinates with respect to com
  double r_bound;
  double x_bound[3];

  // radius of sphere with equal volume
  double r_equiv;

  // number of tries for mc
  int ntry;

  bool overlap_slightly;

  bool no_overlap;
};

}

#endif
#endif
