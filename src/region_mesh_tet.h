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

#ifdef REGION_CLASS

RegionStyle(mesh/tet,RegTetMesh)

#else

#ifndef LMP_REGION_TET_MESH_H
#define LMP_REGION_TET_MESH_H

#include "random_park.h"
#include "region.h"

namespace LAMMPS_NS {

class RegTetMesh : public Region {

 friend class InputMeshTet;

 public:

  RegTetMesh(class LAMMPS *, int, char **);
  ~RegTetMesh();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);

  void generate_random(double *pos,bool subdomain_flag);
  void generate_random_shrinkby_cut(double *pos,double cut,bool subdomain_flag);

  // volume calculation based on MC
  void volume_mc(int n_test,bool cutflag,double cut,double &vol_global,double &vol_local);

  void add_tet(double **n);
  int n_tet();
  double total_vol();
  double tet_vol(int i);
  double tet_acc_vol(int i);

  class TriMesh *get_tri_mesh()
  { return &tri_mesh; }

 protected:

   int is_inside_tet(int iTet,double *pos);
   bool nodesAreEqual(double *nodeToCheck1,double *nodeToCheck2,double precision);

   void grow_arrays();
   void set_extent_region();
   void set_extent_mesh();
   void build_neighs();
   void build_surface();
   double volume_of_tet(double* v0, double* v1, double* v2, double* v3);
   double volume_of_tet(int iTet);

   int mesh_randpos(double *pos);
   int  tet_rand_tri();

   char *filename;
   double scale_fact;
   double off_fact[3], rot_angle[3];

   int nTet,nTetMax;
   double ***node;
   double **center;
   double *rbound, rbound_max;

   int *n_face_neighs;
   int **face_neighs; 

   int **n_face_neighs_node;

   int *n_node_neighs;
   int **node_neighs; 

   int *n_surfaces;
   int **surfaces;

   double total_volume;
   double *volume;
   double *acc_volume;

   class BoundingBox &bounding_box_mesh;

   class RegionNeighborList &neighList;

   class TriMesh &tri_mesh;

   #include "region_mesh_tet_I.h"
};

}

#endif
#endif
