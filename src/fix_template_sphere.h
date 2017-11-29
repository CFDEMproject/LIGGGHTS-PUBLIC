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

#include <vector>
#include <utility>
#include "fix.h"
#include "fix_property_atom.h"
#include "random_park.h"
#include "probability_distribution.h"
#include "vector_liggghts.h"

namespace PARTICLE_PACKING
{

class Sphere
{
public:
    Sphere() :
        pos_x(0.0),
        pos_y(0.0),
        pos_z(0.0),
        radius(0.0),
        density(0.0),
        is_local(false),
        id(0)
    {}

    Sphere(const double * const _x, const double _radius, const double _density, const int _id) :
        pos_x(_x[0]),
        pos_y(_x[1]),
        pos_z(_x[2]),
        radius(_radius),
        density(_density),
        is_local(false),
        id(_id)
    {}

    void move_particle(const double * const shift)
    {
        pos_x += shift[0];
        pos_y += shift[1];
        pos_z += shift[2];
    }

    inline void get_pos(double * const _x) const
    {
        _x[0] = pos_x;
        _x[1] = pos_y;
        _x[2] = pos_z;
    }

    void set_local()
    { is_local = true; }

    void unset_local()
    { is_local = false; }

    bool get_local() const
    { return is_local; }

    double get_radius() const
    { return radius; }

    double get_density() const
    { return density; }

    double get_volume() const
    { return radius*radius*radius*4.18879020478639098; }

    int get_id() const
    { return id; }

    size_t n_fix_properties() const
    { return fix_properties.size(); }

    size_t fix_property_nentries(const int i) const
    { return fix_property_values[i].size(); }

    LAMMPS_NS::FixPropertyAtom* get_fix_property(const int i) const
    { return fix_properties[i]; }

    double fix_property_value(const int i, const int j) const
    { return fix_property_values[i][j]; }

    void init_fix_properties(std::vector<std::pair<LAMMPS_NS::FixPropertyAtom*, int> > &fpa_list)
    {
        std::vector<std::pair<LAMMPS_NS::FixPropertyAtom*, int> >::iterator it = fpa_list.begin();
        size_t n = 0;
        for (; it != fpa_list.end(); it++)
        {
            if (it->first)
                n++;
        }
        fix_properties.resize(n);
        fix_property_values.resize(n);
        for (size_t i = 0, j = 0; i < fpa_list.size(); i++)
        {
            if (!fpa_list[i].first)
                continue;
            fix_properties[j] = fpa_list[i].first;
            size_t m = fpa_list[i].second;
            fix_property_values[j].resize(m);
            j++;
        }
    }

    void set_fix_property_values(const int i, const double * const data)
    {
        for (size_t j = 0; j < fix_property_values[i].size(); j++)
            fix_property_values[i][j] = data[j];
    }

protected:
    double pos_x, pos_y, pos_z, radius, density;
    bool is_local;
    int id;
    std::vector<LAMMPS_NS::FixPropertyAtom*> fix_properties;
    std::vector<std::vector<double> > fix_property_values;
};

}

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

  bool use_rad_for_cut_neigh_and_ghost()
  { return !relative; }

  bool is_relative()
  { return relative; }

  virtual double min_rad();
  virtual double max_rad();
  virtual double max_r_bound();
  virtual int number_spheres();
  virtual int maxtype();
  virtual int mintype();
  class Region *region();

  // single particle generation, used by fix insert/* commands
  virtual void randomize_single();    
  class ParticleToInsert *pti;

  // many particle generation, used by fix insert commands
  virtual void init_ptilist(int n_random_max, const bool enforce_single = false, FixPropertyAtom * const fix_release = NULL);
  virtual void delete_ptilist();
  virtual void randomize_ptilist(int,int,int distorder);
  virtual void direct_set_ptlist(const int i, const void * const data, const int distribution_groupbit, const int distorder);
  int n_pti_max;
  class ParticleToInsert **pti_list;

  // generates hash for identification
  virtual unsigned int generate_hash();

  virtual void finalize_insertion() {}

  inline int random_insertion_state()
  { return random_insertion->state(); }

  inline int random_mc_state()
  { return random_mc->state(); }

 protected:

  int iarg;

  class Region *reg;
  class FixRegionVariable *reg_var;

  // random generator
  
  class RanPark *random_insertion;
  class RanPark *random_mc;
  int seed_insertion;
  int seed_mc;
  int seed_orig;

  // properties of particle template
  int atom_type;
  class LMP_PROBABILITY_NS::PDF *pdf_radius;   
  class LMP_PROBABILITY_NS::PDF *pdf_density;
  double volume_expect;
  double mass_expect;
  double vol_limit;

  bool relative;

  void add_hash_value(const int value, unsigned int &start, unsigned int &hash);
  void add_hash_value(double value, unsigned int &start, unsigned int &hash);
};

}

#endif
#endif
