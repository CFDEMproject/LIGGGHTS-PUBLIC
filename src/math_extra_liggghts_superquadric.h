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
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_LIGGGHTS_SUPERQUADRIC_H
#define LMP_MATH_EXTRA_LIGGGHTS_SUPERQUADRIC_H

#include "math_extra_liggghts_nonspherical.h"
#include "superquadric.h"

using namespace LIGGGHTS::ContactModels;

namespace MathExtraLiggghtsNonspherical {

  void area_superquadric(const double *shape, const double *blockiness, double *result);
  double crossSectionalArea(const double* u_, const double* shape, const double *blockiness, const double* quat, int& index);
  void inertia_superquadric(const double *shape, const double *blockiness, const double density, double *inertia);
  void volume_superquadric(const double *shape, const double *blockiness, double *volume);
  void bounding_sphere_radius_superquadric(const double *shape, const double *blockiness, double *radius);

  double check_inequalities(const int i, const double *delta, const double *a, const double *b,
                          const double *A_delta, const double *B_delta, const double *C);
  bool obb_intersect(Superquadric *particle_i, Superquadric *particle_j);
  bool obb_intersect(Superquadric *particle_i, Superquadric *particle_j, unsigned int &inequality_to_start);

  double minimal_distance(Superquadric *particle1, Superquadric *particle2, const double *initial_point1, const double *initial_point2,
      double *result_point1, double *result_point2, bool *fail);

  void calc_contact_point(Superquadric *particle_i, Superquadric *particle_j,
      double ratio, const double *initial_point1, double *result_point, double &fi, double &fj, bool *fail, LAMMPS_NS::Error *error);

  bool capsules_intersect(Superquadric *particle_i, Superquadric *particle_j, double *capsule_contact_point);

  bool calc_contact_point_if_no_previous_point_avaialable(SurfacesIntersectData & sidata, Superquadric *particle_i, Superquadric *particle_j,
      double *contact_point, double &fi, double &fj,LAMMPS_NS::Error *error);
  bool calc_contact_point_using_prev_step(SurfacesIntersectData & sidata, Superquadric *particle_i, Superquadric *particle_j,
      double ratio, double dt, double *prev_step_point, double *contact_point, double &fi, double &fj, LAMMPS_NS::Error *error);
  void basic_overlap_algorithm(SurfacesIntersectData & sidata, Superquadric *particle_i, Superquadric *particle_j,
      double &alphai, double &alphaj, const double *contact_point, double *contact_point_i, double *contact_point_j);
  double extended_overlap_algorithm(Superquadric *particleA, Superquadric *particleB,
      double *en, double *const alpha_i, double *const alpha_j,
      const double *contact_point, double *contact_pointA, double *contact_pointB, double *delta);
  double common_normal(SurfacesIntersectData & sidata, Superquadric *particle_i, Superquadric *particle_j, bool particles_were_in_contact_flag,
      double *const contact_point_i_local, double *contact_point_j_local, double *contact_point_i, double *contact_point_j);
  double inverseMatrix4x4(const double *m, double *out);
  double determinant_4x4(double *mat);
#ifdef LIGGGHTS_DEBUG
  void printf_debug_data(Superquadric *particle_i, Superquadric *particle_j, double *initial_guess, LAMMPS_NS::Error *error);
#endif
};

#endif

