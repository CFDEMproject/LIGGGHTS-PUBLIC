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

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "fix_insert_stream.h"
#include "mpi_liggghts.h"
#include "fix_property_atom_tracer_stream.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtomTracerStream::FixPropertyAtomTracerStream(LAMMPS *lmp, int narg, char **arg,bool parse) :
  FixPropertyAtomTracer(lmp, narg, arg, false),
  n_marker_per_(-1),
  every_(-1),
  fix_ins_stream_(0)
{
    // error checks on settings from base class

    // parse args for this class

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"n_tracer") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'n_tracer'");
            iarg_++;
            n_marker_per_ = atoi(arg[iarg_++]);
            if(n_marker_per_ < 0)
                error->fix_error(FLERR,this,"n_tracer > 0 required");
            hasargs = true;
        }
        else if(strcmp(arg[iarg_],"insert_stream") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'insert_stream'");
            iarg_++;
            fix_ins_stream_ = static_cast<FixInsertStream*>(modify->find_fix_id_style(arg[iarg_++],"insert/stream"));
            if(!fix_ins_stream_)
                error->fix_error(FLERR,this,"insert_stream ID does not exist");
            fix_ins_stream_->register_tracer_callback(this);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"every") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'every'");
            iarg_++;
            if(strcmp(arg[iarg_],"once") == 0)
                every_ = 0;
            else
                every_ = atoi(arg[iarg_]);
            if(every_ < 0)
                error->fix_error(FLERR,this,"every > 0 required");
            iarg_++;
            hasargs = true;
        } else
        {
            //fprintf(screen,"%s\n",arg[iarg_]);
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    // disallow certain args from base class
    if (iregion_ > -1)
            error->fix_error(FLERR,this,"must not use keyword 'region_mark'");
    if (MARKER_HEAVISIDE == marker_style_)
            error->fix_error(FLERR,this,"must not use keyword 'marker_style'");

    // error checks on necessary args
    if (!fix_ins_stream_)
            error->fix_error(FLERR,this,"expecting keyword 'insert_stream'");
    if (step_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'mark_step'");
    if (n_marker_per_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'n_tracer'");
    if (every_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'every'");

    n_to_mark_.push_back(n_marker_per_);
    mark_steps_.push_back(step_);
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomTracerStream::~FixPropertyAtomTracerStream()
{
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixPropertyAtomTracerStream::init()
{
   
    if (atom->tag_enable == 0)
      error->fix_error(FLERR,this,"requires atoms have IDs");

    if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

    Fix *fix_ms = modify->find_fix_style("multisphere",0);
    if (fix_ms)
        error->warning(FLERR,"calculates the wrong mass in case of multisphere particles!");
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomTracerStream::setmask()
{
    
    int mask = FixPropertyAtom::setmask();
    mask |= INITIAL_INTEGRATE;
    return mask;
}

/* ----------------------------------------------------------------------
   update info on time vs marking
------------------------------------------------------------------------- */

void FixPropertyAtomTracerStream::add_remove_packets()
{
    int ts = update->ntimestep;
    int ins_every = fix_ins_stream_->ins_every();

    // in case of once - delete request if finished
    if(0 == every_)
    {
        if(mark_steps_.size() > 0 && n_to_mark_[0] == 0)
        {
            n_to_mark_.erase(n_to_mark_.begin());
            mark_steps_.erase(mark_steps_.begin());
        }
    }

    // in case of non-once
    if(every_ > 0)
    {
        // delete obsolete tracer packets from the past
        // warn if could not mark all requested tracers
        while(mark_steps_.size() > 0 && (mark_steps_[0] < ts-ins_every || n_to_mark_[0] == 0))
        {
            if(n_to_mark_[0] > 0)
                error->warning(FLERR,"Fix property/atom/tracer/stream: Not "
                               "all requested tracers could be marked");
            n_to_mark_.erase(n_to_mark_.begin());
            mark_steps_.erase(mark_steps_.begin());
        }

        // add packets
        for(int istep = step_+every_; istep < ts+ins_every; istep += every_)
        {
            if(istep > ts)
            {
                
                n_to_mark_.push_back(n_marker_per_);
                mark_steps_.push_back(istep);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomTracerStream::mark_tracers(int ilo, int ihi)
{
    std::vector<Releasedata> releasedata_local, releasedata_global;
    Releasedata one;
    FixPropertyAtom *fix_release = fix_ins_stream_->fix_prop_release();
    double **release_data = fix_release->array_atom;
    double *marker = this->vector_atom;
    int *mask = atom->mask;
    int s, step0;
    int release_step_index = fix_ins_stream_->release_step_index();

    add_remove_packets();

    step0 = mark_steps_[0];

    for(int i = ilo; i < ihi; i++)
    {
        // skip particles of wrong group
        if (!(mask[i] & groupbit))
            continue;

        s = static_cast<int>(release_data[i][release_step_index]);
        if (s >= step0)
        {
            one.id = atom->tag[i];
            one.step = s;
            releasedata_local.push_back(one);
        }
    }

    if(1 < comm->nprocs)
    {
        int *data = 0, *data_all = 0;
        int ndata = 0, ndata_all = 0;

        ndata = construct_data(releasedata_local,data);
        ndata_all = MPI_Allgather_Vector(data, ndata, data_all, world);
        releasedata_global = construct_releasedata_all(data_all,ndata_all);

        if(data) delete []data;
        if(data_all) delete []data_all;
    }
    else
        releasedata_global = releasedata_local;

    sort(releasedata_global.begin(), releasedata_global.end());

    int size_all = releasedata_global.size();
    int nmarked_this = 0;
    int ipacket = 0, n_packet = n_to_mark_.size();

    if(n_packet > 0)
    {
        for(int iall = 0; iall < size_all; iall++)
        {
            int step = releasedata_global[iall].step;

            if(step >= mark_steps_[ipacket])
            {
                nmarked_this++;
                n_to_mark_[ipacket]--;
                int ilocal = atom->map(releasedata_global[iall].id);
                
                if(ilocal >= 0)
                    marker[ilocal] = 1.;
            }

            if(n_to_mark_[ipacket] == 0)
                ipacket++;

            if(ipacket == n_packet) break;
        }
    }

    // add to counter of marked particles
    nmarked_ += nmarked_this;

}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomTracerStream::construct_data(std::vector<Releasedata> data_c, int *&data)
{
    int size = data_c.size();
    int datasize = 2*size;
    data = new int[datasize];

    for(int i = 0; i < size; i++)
    {
        data[2*i+0] = data_c[i].id;
        data[2*i+1] = data_c[i].step;
    }
    return datasize;
}

/* ---------------------------------------------------------------------- */

std::vector<Releasedata> FixPropertyAtomTracerStream::construct_releasedata_all(int *data, int ndata)
{
    std::vector<Releasedata> result;
    Releasedata r;

    for(int i = 0; i < ndata/2; i++)
    {
        r.id = data[i*2+0];
        r.step = data[i*2+1];
        result.push_back(r);
    }
    
    return result;
}

