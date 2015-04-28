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

FixStyle(particletemplate/sphere,FixTemplateSphere)

#else

#ifndef LMP_FIX_TEMPLATE_SPHERE_H
#define LMP_FIX_TEMPLATE_SPHERE_H

#include "fix.h"
#include "probability_distribution.h"

namespace LAMMPS_NS {

class FixTemplateSphere : public Fix {
 public:

  FixTemplateSphere(class LAMMPS *, int, char **);
  ~FixTemplateSphere();

  // inherited from Fix
  virtual void post_create(){}
  virtual int setmask();
  void write_restart(FILE *);
  void restart(char *);

  // access to protected properties
  virtual double volexpect();           
  virtual double massexpect();          
  virtual double min_rad();
  virtual double max_rad();
  virtual double max_r_bound();
  virtual int number_spheres();
  virtual int maxtype();
  virtual int mintype();
  class Region *region();

  // single particle generation, used by fix pour/dev
  virtual void randomize_single();    
  class ParticleToInsert *pti;

  // many particle generation, used by fix insert commands
  virtual void init_ptilist(int);
  virtual void delete_ptilist();
  virtual void randomize_ptilist(int,int);
  int n_pti_max;
  class ParticleToInsert **pti_list;

  virtual void finalize_insertion() {}

 protected:

  int iarg;

  class Region *reg;
  class FixRegionVariable *reg_var;

  // random generator
  
  class RanPark *random_insertion;
  class RanPark *random_mc;
  int seed_insertion;
  int seed_mc;

  // properties of particle template
  int atom_type;
  class LMP_PROBABILITY_NS::PDF *pdf_radius;   
  class LMP_PROBABILITY_NS::PDF *pdf_density;
  double volume_expect;
  double mass_expect;
  double vol_limit;
};

}

#endif
#endif
