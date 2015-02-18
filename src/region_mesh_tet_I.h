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

#ifndef LMP_REGION_MESH_TET_I_H
#define LMP_REGION_MESH_TET_I_H

/* ---------------------------------------------------------------------- */

inline void tet_randpos(int iTet,double *pos)
{
    double bary_coo[4];

    double s = random->uniform();
    double t = random->uniform();
    double u = random->uniform();

    if(s+t > 1.)
    {
        s = 1.-s;
        t = 1.-t;
    }
    if(t+u > 1.)
    {
        double tmp = u;
        u = 1.-s-t;
        t = 1.-tmp;
    }
    else if(s+t+u > 1.)
    {
        double tmp = u;
        u = s+t+u-1.;
        s = 1.-t-tmp;
    }
    bary_coo[0] = 1.-s-t-u;
    bary_coo[1] = s;
    bary_coo[2] = t;
    bary_coo[3] = u;

    bary_to_cart(iTet,bary_coo,pos);

}

/* ---------------------------------------------------------------------- */

inline void bary_to_cart(int iTet,double *bary_coo,double *pos)
{
    for(int i=0;i<3;i++)
       pos[i] = bary_coo[0] * node[iTet][0][i] + bary_coo[1] * node[iTet][1][i] + bary_coo[2] * node[iTet][2][i] + bary_coo[3] * node[iTet][3][i];
}

/* ---------------------------------------------------------------------- */

#endif
