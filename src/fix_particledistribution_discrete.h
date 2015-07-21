
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

FixStyle(particledistribution/discrete,FixParticledistributionDiscrete)
FixStyle(particledistribution/discrete/numberbased,FixParticledistributionDiscrete)
FixStyle(particledistribution/discrete/massbased,FixParticledistributionDiscrete)

#else

#ifndef LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H
#define LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H

#include "fix.h"

enum{RAN_STYLE_CONSTANT_FPD,RAN_STYLE_UNIFORM_FPD,RAN_STYLE_GAUSSIAN_FPD};

namespace LAMMPS_NS {

class FixParticledistributionDiscrete : public Fix {
 public:
  friend class FixPourDev;
  FixParticledistributionDiscrete(class LAMMPS *, int, char **);
  ~FixParticledistributionDiscrete();
  void write_restart(FILE *);
  void restart(char *);

  int setmask();

  double vol_expect();
  double mass_expect();

  int max_type();
  int min_type();

  double min_rad(int);
  double max_rad(int);

  double min_rad()
  { return minrad; }
  double max_rad()
  { return maxrad; }
  double max_r_bound()
  { return maxrbound; }

  int max_nspheres();

  int random_init_single(int);         
  class Region* randomize_single();    

  void random_init_list(int);
  int randomize_list(int,int,int);     

  class ParticleToInsert *pti;
  class ParticleToInsert **pti_list;
  int n_pti, n_pti_max;

  void pre_insert(int n,class FixPropertyAtom *fp = 0,double val = 0.);
  int insert(int n);
  void finalize_insertion();

  inline int n_particletemplates()
  { return ntemplates; }
  inline class FixTemplateSphere** particletemplates()
  { return templates; }

  inline int dist_order(int i)
  { return (i >= ntemplates) ? (-1) : (distorder[i]); }

 protected:

  int ninserted;
  int ninsert;

  class RanPark *random;
  int seed;

  int iarg;

  bool mass_based;

  // particle templates
  int ntemplates;       
  double *distweight;   
  double *cumweight;    
  int *parttogen;       
  int *distorder;       
  class FixTemplateSphere **templates; 

  // mass and volume expectancy of this discrete distribution
  double volexpect;
  double massexpect;

  // min/maximum particle type to be inserted
  int maxtype;
  int mintype;

  // maximum number of spheres of all templates
  int maxnspheres;

  // maximum radius and bounding sphere radius of all templates
  double minrad,maxrad,maxrbound;
};

}

#endif
#endif
