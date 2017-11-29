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

    Arno Mayrhofer (CFDEMresearch GmbH, Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contactproperty/atom,FixContactPropertyAtom)

#else

#ifndef LMP_FIX_CONTACTPROPERTY_ATOM_H
#define LMP_FIX_CONTACTPROPERTY_ATOM_H

#include "fix_contact_history.h"
#include "fix_property_atom.h"
#include "my_page.h"
#include <cmath>
#include "vector_liggghts.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "tri_mesh.h"

namespace LAMMPS_NS {

class FixContactPropertyAtom : public FixContactHistory {

 public:

  FixContactPropertyAtom(class LAMMPS *, int, char **);
  ~FixContactPropertyAtom();
  virtual void post_create();
  virtual int setmask();
  void init();
  //void set_arrays(int i);
  void setup_pre_exchange();
  void min_setup_pre_exchange();
  void pre_exchange();

  void setup_pre_force(int dummy);
  void min_setup_pre_force(int dummy);
  virtual void pre_force(int dummy);
  void min_pre_force(int dummy);

  void pre_neighbor();

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int unpack_exchange(int, double *);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_restart(int, int);
  void write_restart(FILE *fp);
  virtual void clear();

  void do_forward_comm();

  virtual class FixMeshSurface* getMesh() const
  { return NULL; }

  // inline access

  inline int has_partner(int i,int partner_id)
  {
      for(int ip = 0; ip < npartner_[i]; ip++)
      {
          if(partner_id == partner_[i][ip])
          {
            
            return ip;
          }
      }
      return -1;
  }

  void add_partner(int i, int partner_id, const double * const history)
  {
      
      partner_[i][npartner_[i]] = partner_id;
      vectorCopyN(history,&(contacthistory_[i][npartner_[i]*dnum_]),dnum_);
      //double *nneighs = fix_nneighs_full_->vector_atom;
      //int n = static_cast<int>(nneighs[i]);
      //printf("add_p: %p %d\n", &(contacthistory_[i][npartner_[i]*dnum_]), n);
      npartner_[i]++;
  }

  inline void update_partner(const int i, const int j, const double * const history)
  {
      if (j < npartner_[i])
          vectorCopyN(history, &(contacthistory_[i][j*dnum_]), dnum_);
  }

  inline int get_npartners(const int i)
  {
      return npartner_[i];
  }

 protected:

  class FixPropertyAtom *fix_nneighs_full_;

  bool build_neighlist_, reset_each_ts_;
};

}

#endif
#endif
