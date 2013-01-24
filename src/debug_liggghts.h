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

#ifndef LMP_DEBUG_LIGGGHTS_H
#define LMP_DEBUG_LIGGGHTS_H

#include "lammps.h"
#include "string.h"
#include "stdlib.h"
#include "style_fix.h"
#include "vector_liggghts.h"
#include "container.h"

namespace LAMMPS_NS {

inline void __debug__(LAMMPS* lmp)
{
    //fprintf(lmp->screen,"test");

    for(int i = 0; i < lmp->modify->nfix; i++)
    {

        if(strcmp(lmp->modify->fix[i]->style,"move/mesh") == 0)
        {
            FixMoveMesh *fmm = static_cast<FixMoveMesh*>(lmp->modify->fix[i]);
            MultiVectorContainer<double,3,3> *v;
            v = fmm->mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
            if(v) printVec3D(lmp->screen,"VDEBUG",v->begin()[0][0]);

        }
    }
}

}

#endif
