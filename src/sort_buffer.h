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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_SORT_BUFFER_H
#define LMP_SORT_BUFFER_H

#include "pointers.h"
#include "irregular.h"

namespace LAMMPS_NS
{

class SortBuffer : public Pointers
{
  public:
    SortBuffer(LAMMPS *lmp, const bool _sort_flag);
    ~SortBuffer();
    void sort(double * &buf, int &nme, int &maxbuf, const int _size_one, const int ntotal);
    int modify_param(const int narg, const char *const *const arg);
    bigint memory_usage(const int _size_one);
    void init(const int igroup);
    int *get_ids();
    void realloc_ids(const int nmax);

    bool sort_set()
    { return sort_flag; }

    int get_sortcol()
    { return sortcol; }

  private:
    bool sort_flag;             // 1 if sorted output
    int sortcol;                // 0 to sort on ID, 1-N on columns
    int sortcolm1;              // sortcol - 1
    int sortorder;              // ASCEND or DESCEND

    int size_one;

    int maxsort;                // size of bufsort, idsort and index
    double *bufsort;
    int *idsort;
    int *index;

    int maxids;                 // size of ids
    int *ids;                   // list of atom IDs, if sorting on IDs

    int maxproc;                // size of proclist
    int *proclist;

    Irregular *irregular;

    // reordering
    bool reorderflag;           // 1 if OK to reorder instead of sort
    int ntotal_reorder;         // # of atoms that must be in snapshot
    int nme_reorder;            // # of atoms I must own in snapshot
    int idlo;                   // lowest ID I own when reordering

    static int idcompare(const void *, const void *);
    static int bufcompare(const void *, const void *);
    static int bufcompare_reverse(const void *, const void *);
};

}

#endif
