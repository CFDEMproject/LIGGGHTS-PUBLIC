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
