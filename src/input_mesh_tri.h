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

#ifndef LMP_INPUT_MESH_TRI_H
#define LMP_INPUT_MESH_TRI_H

#include "stdio.h"
#include "input.h"

namespace LAMMPS_NS {

class InputMeshTri : protected Input
{
  public:

    InputMeshTri(class LAMMPS *, int, char **);
    ~InputMeshTri();

    void meshtrifile(const char *,class TriMesh *,bool verbose);

  private:

    bool verbose_;

    void meshtrifile_vtk(class TriMesh *);
    void meshtrifile_stl(class TriMesh *);
    inline void addTriangle(class TriMesh *mesh,
         double *a, double *b, double *c,int lineNumber);

};

}

#endif
