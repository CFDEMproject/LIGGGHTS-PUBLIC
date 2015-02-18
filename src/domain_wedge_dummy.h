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

/* ----------------------------------------------------------------------
   Contributing authors:
   Stefan Amberger (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#ifndef DOMAIN_WEDGE_H
#define DOMAIN_WEDGE_H

#include "domain.h"

namespace LAMMPS_NS {

class DomainWedge : public Domain
{

  public:

    DomainWedge(class LAMMPS *lmp) : Domain(lmp) {};
    void set_domain(class RegWedge *rw) {}

    inline int index_axis()
    { return 0; }

    inline int index_phi()
    { return 0; }

    inline void n1(double *_n1)
    { }

    inline void n2(double *_n2)
    { }

    inline void center(double *_c)
    { }
};

}

#endif
