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
    Arno Mayrhofer (DCS Computing GmbH, Linz)
    Stefan Radl (TU Graz)

    Copyright 2017 - DCS Computing GmbH, Linz
    Copyright 2016 - TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_multisphere_break.h"
#include "domain_wedge.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_template_multisphere.h"
#include "neighbor.h"
#include "fix_gravity.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "atom_vec.h"
#include "math_extra_liggghts.h"
#include "math_const.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum {LOOP_LOCAL,LOOP_ALL};
enum {NONE, VARIABLE, FIX};

/* ---------------------------------------------------------------------- */

FixMultisphereBreak::FixMultisphereBreak(LAMMPS *lmp, int narg, char **arg) :
  FixMultisphere(lmp, narg, arg),
  triggerFixName_(NULL),
  triggerFix_(NULL),
  triggerIdx_(-1),
  triggerArray_(NULL),
  maxatom_(0),
  triggerType_(NONE),
  triggerName_(NULL),
  triggerIndex_(-1),
  triggerThreshold_(0),
  triggerTimeStep_(0)
{
    int iarg=3;
    bool trigger_fixNameSet = false;
    bool trigger_nameSet = false;
    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        printf("iarg:%d \n", iarg);
        if(strcmp(arg[iarg],"trigger_threshold") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'trigger_threshold'");
            iarg++;
            triggerThreshold_ = atof(arg[iarg++]);
            printf("FixMultisphereBreak will use trigger_threshold: %g \n",
                    triggerThreshold_);
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"trigger_timeStep") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'trigger_timeStep'");
            iarg++;
            triggerTimeStep_ = atoi(arg[iarg++]);
            printf("FixMultisphereBreak will use trigger_timeStep: %d \n",
                    triggerTimeStep_);
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"trigger_name") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR, this, "not enough arguments for 'trigger_name'");
            iarg++;
            if (arg[iarg][0] == 'f')
                triggerType_ = FIX;
            else if (arg[iarg][0] == 'v')
                triggerType_ = VARIABLE;
            else
                error->fix_error(FLERR, this, "Require a fix with f_ or variable with v_");
            int n = strlen(arg[iarg]);
            char *suffix = new char[n];
            strcpy(suffix, &arg[iarg][2]);

            char *ptr = strchr(suffix,'[');
            if (ptr)
            {
                if (suffix[strlen(suffix)-1] != ']')
                    error->all(FLERR,"Illegal fix multisphere/break command");
                triggerIndex_ = atoi(ptr+1);
                *ptr = '\0';
            }
            else
                triggerIndex_ = 0;

            n = strlen(suffix) + 1;
            triggerName_ = new char[n];
            strcpy(triggerName_, suffix);
            delete [] suffix;
            printf("FixMultisphereBreak will use '%s' (length: %d) as trigger. \n",
                    triggerName_, n);
            trigger_nameSet = true;
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"trigger_fixName") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'trigger_fixName'");
            iarg++;
            int n = strlen(arg[iarg++]);
            triggerFixName_ = new char[n+1];
            strcpy(triggerFixName_,arg[iarg-1]);
            printf("FixMultisphereBreak will use fixPropertyAtom with name '%s' (length: %d) as trigger. \n",
                    triggerFixName_, n);
            trigger_fixNameSet = true;
            hasargs = true;
            error->warning(FLERR, "trigger_fixName is a deprecated argument and will be removed in future version, use trigger_name instead");
        }
        else
        {
            printf("WARNING from FixMultisphereBreak: Unknown keyword '%s'. This might be unproblematic in case the derived class handles the keyword correctly. \n", arg[iarg]);
            iarg++;
            hasargs = true;
        }
    }

    if(!trigger_fixNameSet && !trigger_nameSet)
        error->fix_error(FLERR,this,"No setting made for 'trigger_name'. You must make this setting!");
    if(trigger_fixNameSet && trigger_nameSet)
        error->fix_error(FLERR,this,"Setting made for 'trigger_name' and 'trigger_fixName' only one is allowed (preferably the former).");
}

/* ---------------------------------------------------------------------- */

FixMultisphereBreak::~FixMultisphereBreak()
{
    if (triggerName_ && triggerType_ == VARIABLE)
        memory->destroy(triggerArray_);
    if (triggerName_)
        delete[] triggerName_;

    if (triggerFixName_)
        delete[] triggerFixName_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FixMultisphereBreak::init()
{
    FixMultisphere::init();

    if (triggerName_)
    {
        if (triggerType_ == FIX)
        {
            triggerIdx_ = modify->find_fix(triggerName_);
            if (triggerIdx_ < 0)
                error->fix_error(FLERR,this,"fix with name set as trigger_name not found by fix multisphere/break. Ensure a fix with the proper name is available");
            if (triggerIndex_ == 0)
                triggerFix_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(triggerName_,"property/atom","scalar",1,0,style));
            else
                triggerFix_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(triggerName_,"property/atom","vector",triggerIndex_,0,style));
            if(!triggerFix_)
                error->fix_error(FLERR,this,"fix property/atom with name set as trigger_name not found by fix multisphere/break. Ensure a fix of style property/atom with the proper name is available");
        }
        else if (triggerType_ == VARIABLE)
        {
            triggerIdx_ = input->variable->find(triggerName_);
            if (triggerIdx_ < 0)
                error->fix_error(FLERR, this, "variable with name set as trigger_name not found by fix multisphere/break. Ensure a variable with the proper name is available");
            if (input->variable->atomstyle(triggerIdx_) == 0)
                error->fix_error(FLERR, this, "variable with name set as trigger_name is not of atom style in fix multisphere/break");
        }
    }
    else if (triggerFixName_)
    {
        triggerFix_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(triggerFixName_,"property/atom","scalar",1,0,style));
        if(!triggerFix_)
            error->fix_error(FLERR,this,"triggerFix not found by FixMultisphereBreak! Ensure a fix with the proper name is available!");
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FixMultisphereBreak::pre_neighbor()
{
    int nall = atom->nlocal + atom->nghost;
    double *corner_ghost = fix_corner_ghost_->vector_atom;
    vectorZeroizeN(corner_ghost,nall);

    fw_comm_flag_ = MS_COMM_FW_BODY;
    forward_comm();

    for(size_t irem = 0; irem < fix_remove_.size(); irem++)
        (fix_remove_[irem])->delete_bodies();

    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();
    multisphere_.remap_bodies(body_);
    rev_comm_flag_ = MS_COMM_REV_IMAGE;
    reverse_comm();
    multisphere_.exchange();

    multisphere_.calc_nbody_all();

    multisphere_.generate_map();

    // set deletion flag
    // if any deleted atoms, do re-neigh in 100 steps at latest to remove
    // remainder particles
    double   *delflag =   fix_delflag_->vector_atom;
    double *existflag = fix_existflag_->vector_atom;
    vectorZeroizeN(delflag,atom->nlocal+atom->nghost);
    vectorZeroizeN(existflag,atom->nlocal+atom->nghost);

    if(multisphere_.check_lost_atoms(body_,delflag,existflag,fix_volumeweight_ms_->vector_atom))
        next_reneighbor = update->ntimestep + 100;

    fix_delflag_->do_reverse_comm();
    fix_existflag_->do_reverse_comm();

    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();

    // DO NOT merge delflag and existflag, since we would like to keep atoms that are not in a body
//    int nlocal = atom->nlocal;
//    delflag =   fix_delflag_->vector_atom;
//    existflag = fix_existflag_->vector_atom;
//    for(int i = 0; i < nlocal; i++)
//    {
//            delflag[i] = (MathExtraLiggghts::compDouble(existflag[i],0.,1e-6)) ? 1. : delflag[i];
//    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FixMultisphereBreak::final_integrate()
{
    FixMultisphere::final_integrate();

    int ghostflag = LOOP_ALL;

    int step = update->ntimestep;
    const int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    double *existflag = fix_existflag_->vector_atom;

    double *density = atom->density;
    double *r = atom->radius;
    double *rmass = atom->rmass;

    int nloop = 0;
    int haveRemovedAtoms = 0;
    int haveRemovedAtomsDirectSelect = 0;
    double fourOverThreePi = 4.18879020479;

    if(ghostflag == LOOP_ALL)
        nloop = nlocal+nghost;
    else if(ghostflag == LOOP_LOCAL)
        nloop = nlocal;
    else
        error->all(FLERR,"Illegal call to FixMultisphereBreak::final_integrate()");

    if (triggerType_ == FIX)
        if (triggerName_ && update->ntimestep % (modify->fix[triggerIdx_]->peratom_freq))
            error->all(FLERR,"Fix used in fix multisphere/break not computed at compatible time");

    if (triggerName_ && triggerType_ == VARIABLE)
    {
        if (nlocal > maxatom_)
        {
            maxatom_ = atom->nmax;
            memory->destroy(triggerArray_);
            memory->create(triggerArray_, maxatom_, "multisphere/break:triggerArray_");
        }
        input->variable->compute_atom(triggerIdx_, igroup, triggerArray_, 1, 0);
    }
    else if ((triggerName_ && triggerType_ == FIX && triggerIndex_ == 0) || triggerFixName_)
        triggerArray_ = triggerFix_->vector_atom;
    const bool fixArray = (triggerName_ && triggerType_  == FIX && triggerIndex_ > 0);

    if(step > triggerTimeStep_)
    {
        for (int i = 0; i < nloop; i++)
        {
            double trigger = 0.;
            if (fixArray)
                trigger = triggerFix_->array_atom[i][triggerIndex_-1];
            else
                trigger = triggerArray_[i];

            if (body_[i] < 0 || trigger < triggerThreshold_)
                continue;
            int ibody = map(body_[i]);
            if (ibody < 0)
                continue;
            if(data().nrigid(ibody) < 2)
                continue;

//            printf("FixMultisphereBreak checks nrigid now for body %d/ibody: %d (n_body()/all: %d, %d) \n",
//                  body_[i], ibody, data().n_body(), data().n_body_all());

//            printf("FixMultisphereBreak resets body index for local atom-id %d with trigger %g, threshold: %g \n",
//                  i, trigger,triggerThreshold_);
            data().tagReset(ibody); //re-sets the tag of the body
            haveRemovedAtomsDirectSelect++;
        }

        //loop through atoms and mark for removal of body
        for (int i = 0; i < nloop; i++)
        {
            if (body_[i] < 0)
                continue;
            int ibody = map(body_[i]);
            if (ibody < 0)
                continue;
            if(data().nrigid(ibody) < 2)
                continue; //disregard bodies with only one or no atom

            if (tag(ibody) == -1)
                existflag[i] = 0; //mark for un-marking of body body_[i]   = -2;
        }

        //loop through atoms and ensure they have -2 in body if body is negatively tagged
        for (int i = 0; i < nloop; i++)
        {
            if (body_[i] < 0)
                continue;
            int ibody = map(body_[i]);
            if (ibody < 0)
                continue;
            if (tag(ibody) == -1 && existflag[i]==0)
            {
                body_[i]   = -1;
                if(data().nrigid(ibody)>1) //set zero --> removes the body!
                    data().nrigidReset(ibody, 0);

                if (rmass)
                    rmass[i] = r[i] * r[i] * r[i] * fourOverThreePi * density[i];

                haveRemovedAtoms++;
            }
        }

    }
    if(haveRemovedAtoms>0)
        FixMultisphere::add_body_finalize();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FixMultisphereBreak::calc_force(bool setupflag)
{

//  int ghostflag = LOOP_ALL;
//  int nlocal = atom->nlocal;
//  int nghost = atom->nghost;
//  int nloop = 0;

//  if(ghostflag == LOOP_ALL) nloop = nlocal+nghost;
//  else if(ghostflag == LOOP_LOCAL) nloop = nlocal;
//  else error->all(FLERR,"Illegal call to FixMultisphereBreak::final_integrate()");

  FixMultisphere::calc_force(setupflag);

  //add gravity to atoms not belonging to a body to ensure they move correctly
  if(false) //fix_gravity_)
  {
//      double grav[3];
//      fix_gravity_->get_gravity(grav);
//      if (rmass) {
//        for (int i = 0; i < nloop; i++)
//        {
//            if ( existflag[i]==0 )
//            {
//                massone = rmass[i];
/////                printf("massone: %g grav[2]: %g .\n",massone, grav[2]);
//                f[i][0] += massone*grav[0];
//                f[i][1] += massone*grav[1];
//                f[i][2] += massone*grav[2];
//            }
//        }
//      } else {
//        for (int i = 0; i < nloop; i++)
//        {
//            if ( existflag[i]==0 )
//            {
//                massone = mass[type[i]];
//                f[i][0] += massone*grav[0];
//                f[i][1] += massone*grav[1];
//                f[i][2] += massone*grav[2];
//            }
//        }
//      }
  }

}
