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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <cmath>
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force_implicit.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::FixCfdCouplingForceImplicit(LAMMPS *lmp, int narg, char **arg) :
    FixCfdCouplingForce(lmp,narg,arg),
    useCN_(false),
    CNalpha_(0.0),
    useAM_(false),
    CAddRhoFluid_(0.0),
    onePlusCAddRhoFluid_(1.0),
    fix_Ksl_(0),
    fix_uf_(0),
    fix_KslRotation_(0),
    fix_ex_(0),
    fix_KslExtra_(0)
{

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"CrankNicolson") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CrankNicholson'");
            iarg++;
            useCN_ = true;
            CNalpha_ = atof(arg[iarg]);
            if(CNalpha_<0 || CNalpha_>1)
                error->fix_error(FLERR,this,"incorrect choice for 'CrankNicholson': setting CNalpha_<0 or CNalpha_>1 is not appropriate");

            fprintf(screen,"cfd_coupling_foce_implicit will use Crank-Nicholson scheme with %f\n", CNalpha_);
            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"CAddRhoFluid") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CAddRhoFluid'");
            iarg++;
            useAM_ = true;
            CAddRhoFluid_        = atof(arg[iarg]);
            onePlusCAddRhoFluid_ = 1.0 + CAddRhoFluid_;
            fprintf(screen,"cfd_coupling_force_implicit will consider added mass with CAddRhoFluid = %f\n",
                    CAddRhoFluid_);
            iarg++;
        }
        else
           iarg++;
    }

  nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::~FixCfdCouplingForceImplicit()
{

}

/* ---------------------------------------------------------------------- */
int FixCfdCouplingForceImplicit::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_create()
{
    // do mother class init w/o dragforce
    FixCfdCouplingForce::post_create();

    // register Ksl
    if(!fix_Ksl_)
    {
        const char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Ksl";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register uf
    if(!fix_uf_)
    {
        const char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="uf";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
    if(!fix_KslRotation_)
    {
      const char* fixarg[11];
      fixarg[0]="KslRotation";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="KslRotation";
      fixarg[4]="vector"; // 1 vector per particle to be registered
      fixarg[5]="yes";    // restart
      fixarg[6]="no";     // communicate ghost
      fixarg[7]="no";     // communicate rev
      fixarg[8]= "0.";
      fixarg[9]= "0.";
      fixarg[10]="0.";
      fix_KslRotation_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
     }
    if(!fix_ex_)
    {
       const char* fixarg[11];
       fixarg[0]="ex";
       fixarg[1]="all";
       fixarg[2]="property/atom";
       fixarg[3]="ex";
       fixarg[4]="vector"; // 1 vector per particle to be registered
       fixarg[5]="yes";    // restart
       fixarg[6]="no";     // communicate ghost
       fixarg[7]="no";     // communicate rev
       fixarg[8]= "0.";
       fixarg[9]= "0.";
       fixarg[10]="0.";
       fix_ex_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
    if(!fix_KslExtra_)
    {
      const char* fixarg[11];
      fixarg[0]="KslExtra";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="KslExtra";
      fixarg[4]="vector"; // 1 vector per particle to be registered
      fixarg[5]="yes";    // restart
      fixarg[6]="no";     // communicate ghost
      fixarg[7]="no";     // communicate rev
      fixarg[8]= "0.";
      fixarg[9]= "0.";
      fixarg[10]="0.";
      fix_KslExtra_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_Ksl_) modify->delete_fix("Ksl");
    if(unfixflag && fix_uf_) modify->delete_fix("uf");
    if(unfixflag && fix_KslRotation_) modify->delete_fix("KslRotation");
    if(unfixflag && fix_KslExtra_) modify->delete_fix("KslExtra");
    if(unfixflag && fix_ex_) modify->delete_fix("ex");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::init()
{
    FixCfdCouplingForce::init();

    // values to come from OF
    fix_coupling_->add_pull_property("Ksl","scalar-atom");
    fix_coupling_->add_pull_property("uf","vector-atom");
    fix_coupling_->add_pull_property("KslRotation","vector-atom");
    fix_coupling_->add_pull_property("KslExtra","vector-atom");
    fix_coupling_->add_pull_property("ex","vector-atom");

    deltaT_ = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_force(int)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;
  double **dragforce = fix_dragforce_->array_atom;
  double frc[3];

  vectorZeroize3D(dragforce_total);
  vectorZeroize3D(hdtorque_total);

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        // calc force
        if(!useCN_)  //calculate drag force and add if not using Crank-Nicolson
        {
            vectorSubtract3D(uf[i],v[i],frc);
            vectorScalarMult3D(frc,Ksl[i]);

            vectorAdd3D(f[i],frc,f[i]);
            vectorAdd3D(dragforce_total,frc,dragforce_total);
        }

        // add other force
        vectorAdd3D(f[i],dragforce[i],f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::end_of_step()
{

  if(!useCN_) return; //return if CN not used

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;
  //double **dragforce = fix_dragforce_->array_atom;
  double KslMDeltaT, deltaU;
  double vN32[3];
  double frc[3];

  vectorZeroize3D(dragforce_total);

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
      if (rmass)  KslMDeltaT = Ksl[i]/(rmass[i]*onePlusCAddRhoFluid_)*deltaT_;
      else        KslMDeltaT = Ksl[i]/(mass[type[i]]*onePlusCAddRhoFluid_)*deltaT_;

        for(int dirI=0;dirI<3;dirI++)
        {
            //calculate new velocity
            vN32[dirI] = (  v[i][dirI]
                          + KslMDeltaT
                            *(   uf[i][dirI]
                              - (1.0-CNalpha_)*v[i][dirI]
                             )
                         )
                         /
                         (1.0+KslMDeltaT*CNalpha_);

            //calculate velocity difference and force
            deltaU    =  uf[i][dirI]
                           - (
                                (1.0-CNalpha_)*v[i][dirI]
                               +     CNalpha_ *vN32[dirI]
                             );
           frc[dirI] = Ksl[i] * deltaU;  //force required for the next time step

           //update the particle velocity
           //DCS IS NOW
           v[i][dirI] = vN32[dirI];  //update velocity for a half step!
           //RADL WAS
           //v[i][dirI] += KslMDeltaT/2.0 * deltaU;;  //update velocity for a half step!
        }

         // add force
        vectorAdd3D(f[i],frc,f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,frc,dragforce_total);
     }
  }
}
