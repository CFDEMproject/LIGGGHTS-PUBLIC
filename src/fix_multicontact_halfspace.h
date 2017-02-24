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

    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

/*
 * Multi contact model according to Brodu et al.: Multiple-contact discrete-element model for
 * simulating dense granular media. Phys. Rev. E, 2016
 */
#ifdef FIX_CLASS

FixStyle(multicontact/halfspace,FixMultiContactHalfSpace)

#else

#ifndef LMP_FIX_MULTICONTACT_HALFSPACE_H
#define LMP_FIX_MULTICONTACT_HALFSPACE_H

#include "fix.h"
#include "fix_contact_property_atom.h"

namespace LAMMPS_NS {

// The HistoryData class contains:
// 1.) char: either "p", "m", "w" for pair, mesh, wall primitive, respectively
// 2.) void*: pointer to the history fix or the property/atom fix for wall primitives
// 3.) int: offset of the sumDelta history

class HistoryData {
private:
    char type;
    void *fix_ptr;
    int offset;
public:
    HistoryData(const char c, void* const ptr, const int i) :
        type(c),
        fix_ptr(ptr),
        offset(i)
    { }
    const char get_type()
    { return type; }
    void* const get_ptr()
    { return fix_ptr; }
    const int get_offset()
    { return offset; }

    const int get_npartners(const int i);
    double * const get_data_ptr(const int i, const int j);
    void compute_surfPos(const int i, const int jj, const double * const* x, const double * const data_ptr, double * const surfPos_ij, double * const surfPos_ji, const double F_eps);
    double get_fn(const double * const data_ptr);
    void save_contact_property_atom(const int i, const int jj, const int * const tag, const double * const surfPos_ij, const double * const surfPos_ji, FixContactPropertyAtom *cpa);
};

class FixMultiContactHalfSpace : public Fix {
   public:
    FixMultiContactHalfSpace(class LAMMPS *, int, char **);
    ~FixMultiContactHalfSpace();
    void post_create();

    int setmask();
    void init();

    void setup_pre_force(int);
    void pre_force(int);

   private:

    // pointers to classes holding the data
    class PairGranProxy *pairgran_;

    // model properties
    double geometric_prefactor;

    // Youngs modulus and poisson ratio
    const double *Y, *nu;

    // lists of all offsets for the contact history
    // the HistoryData class contains:
    // 1.) char: either "p", "m", "w" for pair, mesh, wall primitive, respectively
    // 2.) void*: pointer to the history
    // 3.) int: offset of the sumDelta history
    std::vector<HistoryData> history_vector;

    // list of all fixes that contain contact property atom (wall) fixes
    std::vector<FixContactPropertyAtom*> contact_property_atom_vector;

};

}

#endif
#endif
