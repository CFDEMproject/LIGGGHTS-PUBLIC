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
