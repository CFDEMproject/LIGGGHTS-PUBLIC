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

FixStyle(property/atom/tracer/stream,FixPropertyAtomTracerStream)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_TRACER_STREAM_H
#define LMP_FIX_PROPERTY_ATOM_TRACER_STREAM_H

#include "fix_property_atom_tracer.h"
#include <vector>

namespace FixPropertyAtomTracerStreamAux {

    class Releasedata {
     public:
      int id, step;
      
      bool operator<(const Releasedata &rhs) const { return step < rhs.step; }
    };
}

using namespace std;
using FixPropertyAtomTracerStreamAux::Releasedata;

namespace LAMMPS_NS {

class FixPropertyAtomTracerStream : public FixPropertyAtomTracer {

 public:

  FixPropertyAtomTracerStream(class LAMMPS *, int, char **, bool parse = true);
  ~FixPropertyAtomTracerStream();

  void init();
  int setmask();

  void add_remove_packets();
  void mark_tracers(int ilo, int ihi);

 private:

  int construct_data(vector<Releasedata> data_c, int *&data);
  vector<Releasedata> construct_releasedata_all(int *data, int ndata);

  int n_marker_per_;
  int every_;

  vector<int> n_to_mark_;
  vector<int> mark_steps_;

  class FixInsertStream *fix_ins_stream_;

}; //end class

}
#endif
#endif
