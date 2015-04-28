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

#ifndef LMP_CFD_DATACOUPLING_H
#define LMP_CFD_DATACOUPLING_H

#include "pointers.h"

namespace LAMMPS_NS {

class CfdDatacoupling : protected Pointers {
 public:

  CfdDatacoupling(class LAMMPS *lmp, int jarg, int narg, char **arg, class FixCfdCoupling* fc);
  ~CfdDatacoupling();

  int get_iarg() {return iarg_;}

  void add_pull_property(const char*, const char*);
  void add_push_property(const char*, const char*);

  virtual void pull(const char *name, const char *type, void *&ptr, const char *datatype);
  virtual void push(const char *name, const char *type, void *&ptr, const char *datatype);

  virtual void allocate_external(int    **&data, int len2,int len1,int    initvalue);
  virtual void allocate_external(double **&data, int len2,int len1,double initvalue);
  virtual void allocate_external(int    **&data, int len2,char *keyword,int initvalue);
  virtual void allocate_external(double **&data, int len2,char *keyword,double initvalue);

  void init();
  virtual void post_create() {}

  virtual bool error_push()
  { return true;}

  // exchange data with OF
  // does nothing in case of MPI coupling
  // for the MPI case, this is done withing the OF solver
  virtual void exchange() = 0;
  void check_datatransfer();

 protected:

  void grow_();

  // used to find properties
  virtual void* find_pull_property(const char *name, const char *type, int &len1, int &len2);
  virtual void* find_push_property(const char *name, const char *type, int &len1, int &len2);

  // data members

  bool liggghts_is_active;
  bool is_parallel;

  // max # of values stored in pull/push lists
  int nvalues_max_;

  // ------------------------------------
  // per atom or global values stored inchar
  // property fixes used in this fix
  // they are pulled from OF each coupling ts
  // ------------------------------------

  // number of values stored/used in this fix
  int npull_;
  // types can be scalar, vector
  char **pullnames_;
  char **pulltypes_;
  // flag used to check if transfer invoked - only if liggghts is not active
  int *pullinvoked_;

  // ------------------------------------
  // values stored in atom or a fix property
  // that are pushed to OF each coupling ts
  // ------------------------------------

  int npush_;
  char **pushnames_;
  char **pushtypes_;
  // flag used to check if transfer invoked - only if liggghts is not active
  int *pushinvoked_;

  int iarg_;
  class FixCfdCoupling *fc_;

  // reference to Properties class in PairGran
  class Properties *properties_;
};

}

#endif
