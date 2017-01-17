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

#include "sort_buffer.h"

#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "force.h"
#include "memory.h"
#include "comm.h"
#include "group.h"
#include "atom.h"

#define BIG 1.0e20
#define IBIG 2147483647
#define EPSILON 1e-6

enum{ASCEND,DESCEND};

// This pointer to self is used for qsort
SortBuffer *sortptr;

SortBuffer::SortBuffer(LAMMPS *lmp, const bool _sort_flag) :
    Pointers(lmp),
    sort_flag(_sort_flag),
    sortcol(0),
    sortcolm1(0),
    sortorder(ASCEND),
    size_one(0),
    maxsort(0),
    bufsort(NULL),
    idsort(NULL),
    index(NULL),
    maxids(0),
    ids(NULL),
    maxproc(0),
    proclist(NULL),
    irregular(NULL),
    reorderflag(false),
    ntotal_reorder(0),
    nme_reorder(0),
    idlo(0)
{
}

SortBuffer::~SortBuffer()
{
    memory->destroy(bufsort);
    memory->destroy(idsort);
    memory->destroy(index);
    memory->destroy(ids);
    memory->destroy(proclist);
    delete irregular;
}

void SortBuffer::init(const int igroup)
{
    if (!sort_flag)
    {
        sortcol = 0;
        sortcolm1 = 0;
        sortorder = ASCEND;

        maxsort = 0;
        memory->destroy(bufsort);
        memory->destroy(idsort);
        memory->destroy(index);

        maxids = 0;
        memory->destroy(ids);

        maxproc = 0;
        memory->destroy(proclist);

        delete irregular;
    }
    else
    {
        if (comm->nprocs > 1 && irregular == NULL)
            irregular = new Irregular(lmp);

        // set reorderflag = 1 if can simply reorder local atoms rather than sort
        // criteria: sorting by ID, atom IDs are consecutive from 1 to Natoms
        //           min/max IDs of group match size of group
        // compute ntotal_reorder, nme_reorder, idlo/idhi to test against later

        bigint size = group->count(igroup);
        if (size > MAXSMALLINT)
            error->all(FLERR,"Too many atoms to dump sort");

        reorderflag = false;
        if (sortcol == 0 && atom->tag_consecutive())
        {
            int *tag = atom->tag;
            int *mask = atom->mask;
            int nlocal = atom->nlocal;

            int min = IBIG;
            int max = 0;
            for (int i = 0; i < nlocal; i++)
            {
                if (mask[i] & group->bitmask[igroup])
                {
                    min = MIN(min,tag[i]);
                    max = MAX(max,tag[i]);
                }
            }
            int minall,maxall;
            MPI_Allreduce(&min,&minall,1,MPI_INT,MPI_MIN,world);
            MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
            int isize = static_cast<int> (size);

            if (maxall-minall+1 == isize)
            {
                reorderflag = 1;
                double range = maxall-minall + EPSILON;
                const int me = comm->me;
                const int nprocs = comm->nprocs;
                idlo = static_cast<int> (range*me/nprocs + minall);
                int idhi = static_cast<int> (range*(me+1)/nprocs + minall);

                int lom1 = static_cast<int> ((idlo-1-minall)/range * nprocs);
                int lo = static_cast<int> ((idlo-minall)/range * nprocs);
                int him1 = static_cast<int> ((idhi-1-minall)/range * nprocs);
                int hi = static_cast<int> ((idhi-minall)/range * nprocs);
                if (me && me == lom1)
                    idlo--;
                else if (me && me != lo)
                    idlo++;
                if (me+1 == him1)
                    idhi--;
                else if (me+1 != hi)
                    idhi++;

                nme_reorder = idhi-idlo;
                ntotal_reorder = isize;
            }
        }
    }
}

int SortBuffer::modify_param(const int narg, const char *const *const arg)
{
    int iarg = 0;
    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg],"sort") == 0)
        {
            if (iarg+2 > narg)
                error->all(FLERR,"Illegal sort buffer command");
            if (strcmp(arg[iarg+1],"off") == 0)
                sort_flag = false;
            else if (strcmp(arg[iarg+1],"id") == 0)
            {
                sort_flag = true;
                sortcol = 0;
                sortorder = ASCEND;
            }
            else
            {
                sort_flag = true;
                sortcol = force->inumeric(FLERR,arg[iarg+1]);
                sortorder = ASCEND;
                if (sortcol == 0)
                    error->all(FLERR,"Illegal dump_modify command");
                if (sortcol < 0)
                {
                    sortorder = DESCEND;
                    sortcol = -sortcol;
                }
                sortcolm1 = sortcol - 1;
            }
            iarg += 2;
            hasargs = true;
        }
        else
            return iarg;
    }
    return iarg;
}

/* ---------------------------------------------------------------------- */

bigint SortBuffer::memory_usage(const int _size_one)
{
    size_one = _size_one;
    bigint bytes = 0;
    if (sort_flag)
    {
        if (sortcol == 0)
            bytes += memory->usage(ids,maxids);
        bytes += memory->usage(bufsort,size_one*maxsort);
        if (sortcol == 0)
            bytes += memory->usage(idsort,maxsort);
        bytes += memory->usage(index,maxsort);
        bytes += memory->usage(proclist,maxproc);
        if (irregular)
            bytes += irregular->memory_usage();
    }
    return bytes;
}

/* ---------------------------------------------------------------------- */

int * SortBuffer::get_ids()
{
    if (sort_flag && sortcol == 0)
        return ids;
    else
        return NULL;
}

/* ---------------------------------------------------------------------- */

void SortBuffer::realloc_ids(const int nmax)
{
    if (sort_flag && sortcol == 0 && nmax > maxids)
    {
        maxids = nmax;
        memory->destroy(ids);
        memory->create(ids,maxids,"dump:ids");
    }
}

/* ----------------------------------------------------------------------
   parallel sort of buf across all procs
   changes nme, reorders datums in buf, grows buf if necessary
------------------------------------------------------------------------- */

void SortBuffer::sort(double * &buf, int &nme, int &maxbuf, const int _size_one, const int ntotal)
{
    size_one = _size_one;
    if (!sort_flag)
        return;

    int i,iproc;
    double value;

    const int nprocs = comm->nprocs;

    // if single proc, swap ptrs to buf,ids <-> bufsort,idsort

    if (nprocs == 1) {
        if (nme > maxsort) {
            maxsort = nme;
            memory->destroy(bufsort);
            memory->create(bufsort,maxsort*size_one,"dump:bufsort");
            memory->destroy(index);
            memory->create(index,maxsort,"dump:index");
            if (sortcol == 0) {
                memory->destroy(idsort);
                memory->create(idsort,maxsort,"dump:idsort");
            }
        }

        double *dptr = buf;
        buf = bufsort;
        bufsort = dptr;

        if (sortcol == 0) {
            int *iptr = ids;
            ids = idsort;
            idsort = iptr;
        }

    // if multiple procs, exchange datums between procs via irregular

    } else {

        // grow proclist if necessary

        if (nme > maxproc) {
            maxproc = nme;
            memory->destroy(proclist);
            memory->create(proclist,maxproc,"dump:proclist");
        }

        // proclist[i] = which proc Ith datum will be sent to

        if (sortcol == 0) {
            int min = IBIG;
            int max = 0;
            for (i = 0; i < nme; i++) {
                min = MIN(min,ids[i]);
                max = MAX(max,ids[i]);
            }
            int minall,maxall;
            MPI_Allreduce(&min,&minall,1,MPI_INT,MPI_MIN,world);
            MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
            double range = maxall-minall + EPSILON;
            for (i = 0; i < nme; i++) {
                iproc = static_cast<int> ((ids[i]-minall)/range * nprocs);
                proclist[i] = iproc;
            }

        } else {
            double min = BIG;
            double max = -BIG;
            for (i = 0; i < nme; i++) {
                value = buf[i*size_one + sortcolm1];
                min = MIN(min,value);
                max = MAX(max,value);
            }
            double minall,maxall;
            MPI_Allreduce(&min,&minall,1,MPI_DOUBLE,MPI_MIN,world);
            MPI_Allreduce(&max,&maxall,1,MPI_DOUBLE,MPI_MAX,world);
            double range = maxall-minall + EPSILON*(maxall-minall);
            if (range == 0.0) range = EPSILON;
            for (i = 0; i < nme; i++) {
                value = buf[i*size_one + sortcolm1];
                iproc = static_cast<int> ((value-minall)/range * nprocs);
                proclist[i] = iproc;
            }
        }

        // create comm plan, grow recv bufs if necessary,
        // exchange datums, destroy plan
        // if sorting on atom IDs, exchange IDs also

        nme = irregular->create_data(nme,proclist);

        if (nme > maxsort) {
            maxsort = nme;
            memory->destroy(bufsort);
            memory->create(bufsort,maxsort*size_one,"dump:bufsort");
            memory->destroy(index);
            memory->create(index,maxsort,"dump:index");
            if (sortcol == 0) {
                memory->destroy(idsort);
                memory->create(idsort,maxsort,"dump:idsort");
            }
        }

        irregular->exchange_data((char *) buf,size_one*sizeof(double),
                                                         (char *) bufsort);
        if (sortcol == 0)
            irregular->exchange_data((char *) ids,sizeof(int),(char *) idsort);
        irregular->destroy_data();
    }

    // if reorder flag is set & total/per-proc counts match pre-computed values,
    // then create index directly from idsort
    // else quicksort of index using IDs or buf column as comparator

    if (reorderflag) {
        if (ntotal != ntotal_reorder) reorderflag = 0;
        int flag = 0;
        if (nme != nme_reorder) flag = 1;
        int flagall;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
        if (flagall) reorderflag = 0;

        if (reorderflag)
            for (i = 0; i < nme; i++)
                index[idsort[i]-idlo] = i;
    }

    if (!reorderflag) {
        sortptr = this;
        for (i = 0; i < nme; i++) index[i] = i;
        if (sortcol == 0) qsort(index,nme,sizeof(int),idcompare);
        else if (sortorder == ASCEND) qsort(index,nme,sizeof(int),bufcompare);
        else qsort(index,nme,sizeof(int),bufcompare_reverse);
    }

    // reset buf size and maxbuf to largest of any post-sort nme values
    // this insures proc 0 can receive everyone's info

    int nmax;
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);

    if (nmax > maxbuf) {
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*size_one,"dump:buf");
    }

    // copy data from bufsort to buf using index

    int nbytes = size_one*sizeof(double);
    for (i = 0; i < nme; i++)
        memcpy(&buf[i*size_one],&bufsort[index[i]*size_one],nbytes);
}

/* ----------------------------------------------------------------------
   compare two atom IDs
   called via qsort() in sort() method
   is a static method so access data via sortptr
------------------------------------------------------------------------- */

int SortBuffer::idcompare(const void *pi, const void *pj)
{
  int *idsort = sortptr->idsort;

  int i = *((int *) pi);
  int j = *((int *) pj);

  if (idsort[i] < idsort[j]) return -1;
  if (idsort[i] > idsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via sortptr
   sort in ASCENDing order
------------------------------------------------------------------------- */

int SortBuffer::bufcompare(const void *pi, const void *pj)
{
  double *bufsort = sortptr->bufsort;
  int size_one = sortptr->size_one;
  int sortcolm1 = sortptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] < bufsort[j]) return -1;
  if (bufsort[i] > bufsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer values with size_one stride
   called via qsort() in sort() method
   is a static method so access data via sortptr
   sort in DESCENDing order
------------------------------------------------------------------------- */

int SortBuffer::bufcompare_reverse(const void *pi, const void *pj)
{
  double *bufsort = sortptr->bufsort;
  int size_one = sortptr->size_one;
  int sortcolm1 = sortptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] > bufsort[j]) return -1;
  if (bufsort[i] < bufsort[j]) return 1;
  return 0;
}
