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

    Andreas Eitzlmayr (TU Graz)

    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_sph_general_gap.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "modify.h"
#include "pair_sph_artvisc_tenscorr.h"

using namespace LAMMPS_NS;
using namespace LIGGGHTS::ContactModels;

/* ---------------------------------------------------------------------- */

FixWallSphGeneralGap::FixWallSphGeneralGap(LAMMPS *lmp, int narg, char **arg) :
  FixWallSphGeneralBase(lmp, narg, arg)
{
    if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments.");

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        if(strcmp(arg[iarg_++],"r0"))
            error->fix_error(FLERR,this,"illegal argument, expecting keyword 'r0'");
        r0 = force->numeric(FLERR,arg[iarg_++]);
        if(strcmp(arg[iarg_++],"D"))
            error->fix_error(FLERR,this,"illegal argument, expecting keyword 'D'");
        D  = force->numeric(FLERR,arg[iarg_++]);
        if(strcmp(arg[iarg_++],"gap"))
            error->fix_error(FLERR,this,"illegal argument, expecting keyword 'gap'");
        gapWidth  = force->numeric(FLERR,arg[iarg_++]);

        if (iarg_ < narg) {
          if(strcmp(arg[iarg_++],"vwall")==0) {
            StaticWall = 1;
            vwallX  = force->numeric(FLERR,arg[iarg_++]);
            vwallY  = force->numeric(FLERR,arg[iarg_++]);
            vwallZ  = force->numeric(FLERR,arg[iarg_++]);
          } else StaticWall = 0;
        }
    }

    if (r0 <= 0. || D < 0. )
      error->fix_error(FLERR,this,"values for r0 or D are invalid");

    set_r0(r0);

    mass_type = atom->avec->mass_type;
    int ntypes = atom->ntypes;
    fixName = arg[0];

    if (mass_type) {
      fppaSlType=static_cast<FixPropertyGlobal*>(modify->find_fix_property("sl","property/global","peratomtype",ntypes,0,force->pair_style));
    } else {
      fppaSl=static_cast<FixPropertyAtom*>(modify->find_fix_property("sl","property/atom","scalar",0,0,force->pair_style,false));
    }

    cs=static_cast<FixPropertyGlobal*>(modify->find_fix_property("speedOfSound","property/global","peratomtype",ntypes,0,force->pair_style));

    fix_wallCount_ = NULL;
    fix_wallContact2_ = NULL;
    fix_wallContact3_ = NULL;
    fix_wallForce2_ = NULL;
    fix_wallForce3_ = NULL;
    fix_usedGapmodel_ = NULL;
}

/* ---------------------------------------------------------------------- */

FixWallSphGeneralGap::~FixWallSphGeneralGap()
{

}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralGap::post_create()
{
    FixWallSphGeneralBase::post_create();
    char* fixID;
    char* myWallContact;
    char* myWallForce;

    myWallContact=new char[30];
    strcpy(myWallContact,"wallContact_");
    strcat(myWallContact,fixName);

    myWallForce=new char[30];
    strcpy(myWallForce,"F_");
    strcat(myWallForce,fixName);

    int ifix_wallCount = -1;
    int ifix_fgradP = -1;
    int ifix_usedGapmodel = -1;

    for (int ifix = 0; ifix < modify->nfix; ifix++)
    {
      fixID = modify->fix[ifix]->id;

      if (strcmp("wallCount",fixID) == 0) {
        ifix_wallCount = ifix;
      }

      if (strncmp("wallContact_",fixID,12) == 0) {
        if (strcmp(fixID,myWallContact) != 0) {

          if (!fix_wallContact2_)
            fix_wallContact2_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixID,"property/atom","vector",0,0,"FixWallGSphGeneralGap",false));

          else if (!fix_wallContact3_)
            fix_wallContact3_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixID,"property/atom","vector",0,0,"FixWallGSphGeneralGap",false));

          else error->fix_error(FLERR,this,"FixWallSphGeneralGap can currently not handle more than 3 meshes.");
        }
      }

      if (strncmp("F_",fixID,2) == 0) {
        if (strcmp(fixID,myWallForce) != 0) {

          if (!fix_wallForce2_)
            fix_wallForce2_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixID,"property/atom","vector",0,0,"FixWallGSphGeneralGap",false));

          else if (!fix_wallForce3_)
            fix_wallForce3_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixID,"property/atom","vector",0,0,"FixWallGSphGeneralGap",false));

          else error->fix_error(FLERR,this,"FixWallSphGeneralGap can currently not handle more than 3 meshes.");
        }
      }

      if (strcmp("fgradP",fixID) == 0) {
        ifix_fgradP = ifix;
      }

      if (strcmp("usedGapmodel",fixID) == 0) {
        ifix_usedGapmodel = ifix;
      }
    }

    if (ifix_wallCount == -1) {
      const char *fixarg[9];
      fixarg[0]="wallCount";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="wallCount";
      fixarg[4]="scalar";
      fixarg[5]="yes";
      fixarg[6]="yes";
      fixarg[7]="yes";
      fixarg[8]="0";
      fix_wallCount_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    } else {
      fix_wallCount_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallCount","property/atom","scalar",0,0,"FixWallGSphGeneralGap",false));
    }

    if (ifix_usedGapmodel == -1) {
      const char *fixarg[9];
      fixarg[0]="usedGapmodel";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="usedGapmodel";
      fixarg[4]="scalar";
      fixarg[5]="yes";
      fixarg[6]="yes";
      fixarg[7]="yes";
      fixarg[8]="0";
      fix_usedGapmodel_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    } else {
      fix_usedGapmodel_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("usedGapmodel","property/atom","scalar",0,0,"FixWallGSphGeneralGap",false));
    }

    if (ifix_fgradP == -1) error->fix_error(FLERR,this,"fix/wall/sph/general/gap needs a sph pairstyle.");
    else fix_fgradP_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("fgradP","property/atom","vector",0,0,"FixWallGSphGeneralGap",false));
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralGap::init()
{
    char* fixID;
    int ifix_integrity = -1;

    for (int ifix = 0; ifix < modify->nfix; ifix++)
    {
      fixID = modify->fix[ifix]->id;

      if (strcmp("int",fixID) == 0) {
        ifix_integrity = ifix;
      }
    }

     if (ifix_integrity == -1) error->fix_error(FLERR,this,"fix/wall/sph/general/gap needs a fix/sph/integrity.");
     else fix_integrity_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("int","property/atom","scalar",0,0,"FixWallGSphGeneralGap",false));
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralGap::pre_delete(bool unfixflag)
{
    FixWallSphGeneralBase::pre_delete(unfixflag);

    if(unfixflag && fix_wallCount_)
        modify->delete_fix(fix_wallCount_->id);

    if(unfixflag && fix_wallContact2_)
        modify->delete_fix(fix_wallContact2_->id);

    if(unfixflag && fix_wallContact3_)
        modify->delete_fix(fix_wallContact3_->id);

    if(unfixflag && fix_wallForce2_)
        modify->delete_fix(fix_wallForce2_->id);

    if(unfixflag && fix_wallForce3_)
        modify->delete_fix(fix_wallForce3_->id);

    if(unfixflag && fix_fgradP_)
        modify->delete_fix(fix_fgradP_->id);

    if(unfixflag && fix_integrity_)
        modify->delete_fix(fix_integrity_->id);

    if(unfixflag && fix_usedGapmodel_)
        modify->delete_fix(fix_usedGapmodel_->id);
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralGap::post_integrate()
{
    FixWallSphGeneralBase::post_integrate();

    // reset fix_wallCount_ and fix_usedGapmodel_
    int i,nlocal = atom->nlocal;
    wallCount_ = fix_wallCount_->vector_atom;
    usedGapmodel_ = fix_usedGapmodel_->vector_atom;

    for (i = 0; i < nlocal; i++) {
      wallCount_[i] = 0;
      usedGapmodel_[i] = 0;
    }
}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneralGap::compute_density(int ip,double r,double mass)
{
    double *rho = atom->rho;
    int  *type = atom->type;
    int itype;
    double h3,q,q2,q3,q4,q5,F,sli; // for pressure, density boundary contribution
    double sl_corr;

    // Get smoothing length sl:
    if (mass_type) {
      slType = fppaSlType->get_values();
      itype = type[ip];
      sli = slType[itype-1];
    } else {
      sl = fppaSl->vector_atom;
      sli = sl[ip];
    }

    h3 = sli*sli*sli;

    // correction for ratio smoothing length/particle spacing different from 1.2:
    sl_corr = h3 * rho0 / (1.728 * mass); // 1.2^3 = 1.728

    // wall coordinate:
    q = r / sli;
    q2 = q*q;
    q3 = q*q2;
    q4 = q2*q2;
    q5 = q2*q3;

    //if (kernel_id == 2) {
        // this polynomial was developed for the cubic spline kernel,
        // reasonable results with wendland and spiky kernel can be assumed
    if (q < 1.31) F = 0.13 * q5 - 0.7 * q4 + 1.14 * q3 - 0.072 * q2 - 1.33 * q + 0.863;
    else F = 0;
    F = F * sl_corr;
    //} else {
    //  error->fix_error(FLERR,this,"wall potentials for chosen kernel not available");
    //}

    // add wall contribution
    rho[ip] += mass * F / h3;
}

/* ----------------------------------------------------------------------
   compute velocity gradients
------------------------------------------------------------------------- */

void FixWallSphGeneralGap::compute_velgrad(int ip,double delx,double dely, double delz,double mass,double *vwall)
{
    double *rho = atom->rho;
    double **v = atom->vest;
    int  *type = atom->type;
    int itype;
    double rsq,r,h3,h4,q,q2,q3,q4,F,sli;
    double delvx,delvy,delvz,mF_rhoh4r;
    double sl_corr;

    // Get smoothing length sl:
    if (mass_type) {
      slType = fppaSlType->get_values();
      itype = type[ip];
      sli = slType[itype-1];
    } else {
      sl = fppaSl->vector_atom;
      sli = sl[ip];
    }

    h3 = sli * sli * sli;

    // correction for ratio smoothing length/particle spacing different from 1.2:
    sl_corr = h3 * rho0 / (1.728 * mass); // 1.2^3 = 1.728

    // wall distance:
    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    if (r > 0) { // avoid error if r = 0
      // wall coordinate:
      q = r / sli;

      if (q < 1.45) {
        if (StaticWall==1) {
          vwall[0] = vwallX;
          vwall[1] = vwallY;
          vwall[2] = vwallZ;
        }

        q2 = q*q;
        q3 = q*q2;
        q4 = q2*q2;

        //if (kernel_id == 2) {
            // this polynomial was developed for the cubic spline kernel,
            // reasonable results with wendland and spiky kernel can be assumed
        F = 1.37 - 0.607 * q4 + 2.59 * q3 - 3.09 * q2 - 0.059 * q;
        F = F * sl_corr;
        //} else {
        //  error->fix_error(FLERR,this,"wall potentials for chosen kernel not available");
        //}

        h4 = h3*sli;

        mF_rhoh4r = - mass * F / (rho[ip] * h4* r) * (1 + r0/r);
        delvx = vwall[0] - v[ip][0];
        delvy = vwall[1] - v[ip][1];
        delvz = vwall[2] - v[ip][2];

        // add wall contribution
        dvdx_[ip][0] += mF_rhoh4r * delvx * delx;
        dvdx_[ip][1] += mF_rhoh4r * delvy * delx;
        dvdx_[ip][2] += mF_rhoh4r * delvz * delx;
        dvdy_[ip][0] += mF_rhoh4r * delvx * dely;
        dvdy_[ip][1] += mF_rhoh4r * delvy * dely;
        dvdy_[ip][2] += mF_rhoh4r * delvz * dely;
        dvdz_[ip][0] += mF_rhoh4r * delvx * delz;
        dvdz_[ip][1] += mF_rhoh4r * delvy * delz;
        dvdz_[ip][2] += mF_rhoh4r * delvz * delz;
      }
    }
}

/* ----------------------------------------------------------------------
   compute force on particle
------------------------------------------------------------------------- */

void FixWallSphGeneralGap::compute_force(SurfacesIntersectData & sidata, double *vwall)
{
  const int ip = sidata.i;
  const double deltan = sidata.deltan;
  const double r = sidata.r;
  const double mass = sidata.mi;
  double * const c_history = sidata.contact_history;
  const double dx = sidata.delta[0];
  const double dy = sidata.delta[1];
  const double dz = sidata.delta[2];
  const double area_ratio = sidata.area_ratio;

    int otherWall;
    double wallDistance;
    double deltan2,r2=0.0,dx2=0.0,dy2=0.0,dz2=0.0,vwall2[3];

    wallCount_ = fix_wallCount_->vector_atom;
    wallCount_[ip] += 1;
    otherWall = 0;

    if (StaticWall==1) {
      vwall[0] += vwallX;
      vwall[1] += vwallY;
      vwall[2] += vwallZ;
    }

    if (wallCount_[ip] > 1) { // contact with more than 1 wall (maybe a gap)

      if (!fix_wallContact3_) {
        wallContact2_ = fix_wallContact2_->array_atom;

        otherWall = 2;
        wallDistance = deltan + wallContact2_[ip][0];
      } else {
        wallContact2_ = fix_wallContact2_->array_atom;
        wallContact3_ = fix_wallContact3_->array_atom;

        if (wallContact3_[ip][0] == 0) {
          otherWall = 2;
          wallDistance = deltan + wallContact2_[ip][0];
        } else if (wallContact2_[ip][0] == 0) {
          otherWall = 3;
          wallDistance = deltan + wallContact3_[ip][0];
        } else {
          if (wallContact2_[ip][0] < wallContact3_[ip][0]) {
            otherWall = 2;
            wallDistance = deltan + wallContact2_[ip][0];
          } else {
            otherWall = 3;
            wallDistance = deltan + wallContact3_[ip][0];
          }
        }
      }

      if (wallDistance > gapWidth) {  //no gap
        otherWall = 0;
      }
    }

    if (otherWall == 0) { // single wall model
      compute_force_eval_single(ip,deltan,r,mass,dx,dy,dz,vwall,c_history,area_ratio);
    } else { // use gap model
      if (otherWall == 2) {

        wallForce2_ = fix_wallForce2_->array_atom;

        deltan2 = wallContact2_[ip][0];

        dx2 = - wallContact2_[ip][1] / deltan2;
        dy2 = - wallContact2_[ip][2] / deltan2;
        dz2 = - wallContact2_[ip][3] / deltan2;

        vwall2[0] = wallContact2_[ip][4];
        vwall2[1] = wallContact2_[ip][5];
        vwall2[2] = wallContact2_[ip][6];

        r2 = fabs(deltan2);

      } else if (otherWall == 3) {

        wallForce2_ = fix_wallForce3_->array_atom;

        deltan2 = wallContact3_[ip][0];

        dx2 = - wallContact3_[ip][1] / deltan2;
        dy2 = - wallContact3_[ip][2] / deltan2;
        dz2 = - wallContact3_[ip][3] / deltan2;

        vwall2[0] = wallContact3_[ip][4];
        vwall2[1] = wallContact3_[ip][5];
        vwall2[2] = wallContact3_[ip][6];

        r2 = fabs(deltan2);
      }
      // in case of gap: calculate single wall before for wall contribution to density:
      compute_force_eval_single(ip,deltan,r,mass,dx,dy,dz,vwall,c_history,area_ratio);
      compute_force_eval_gap(ip,mass,r,dx,dy,dz,vwall,r2,dx2,dy2,dz2,vwall2);
    }
}

/* ------------------------------------------------------------------------- */

void FixWallSphGeneralGap::compute_force_eval_single(int ip,double deltan,double r,double mass,
                            double dx,double dy,double dz,double *vwall,
                            double *c_history,double area_ratio)
{
    // dx, dy, dz is normalized in case SPH

    double **f = atom->f;
    double *p = atom->p;
    double *rho = atom->rho;
    double **v = atom->vest;
    double *drho = atom->drho;
    int  *type = atom->type;
    double fwallRep=0.0;
    double frac; // for penetration force
    double h2,h3,h4,h5,q,q2,q3,q4,q5,F,vnormalabs,sli; // for pressure, density boundary contribution
    double Fnorm, Ftan,vrelative[3],vtan[3],vtanabs,fwallShear=0.0,factor,fwall[3]; // for viscous boundary contribution
    int itype;
    double rhoi,etai;
    double sl_corr;

    // Get smoothing length sl:
    if (mass_type) {
      slType = fppaSlType->get_values();
      itype = type[ip];
      sli = slType[itype-1];
    } else {
      sl = fppaSl->vector_atom;
      sli = sl[ip];
    }

    h2 = sli*sli;
    h3 = h2 * sli;

    // correction for ratio smoothing length/particle spacing different from 1.2:
    sl_corr = h3 * rho0 / (1.728 * mass); // 1.2^3 = 1.728

    // wall coordinate:

    q = r / sli;

    if (q < 1.45) {

      vrelative[0] = vwall[0] - v[ip][0];
      vrelative[1] = vwall[1] - v[ip][1];
      vrelative[2] = vwall[2] - v[ip][2];

      vnormalabs = vrelative[0]*dx + vrelative[1]*dy + vrelative[2]*dz;

      vtan[0] = vrelative[0] - vnormalabs * dx;
      vtan[1] = vrelative[1] - vnormalabs * dy;
      vtan[2] = vrelative[2] - vnormalabs * dz;

      vtanabs = sqrt(vtan[0]*vtan[0] + vtan[1]*vtan[1] + vtan[2]*vtan[2]);

      q2 = q*q;
      q3 = q*q2;
      q4 = q2*q2;
      q5 = q4*q;

      h4 = h2*h2;

      rhoi = rho[ip];

      if (modelStyle == 1) { // Newtonian
        etai = viscosity;
      } else {    // non-Newtonian
        visc_ = fix_visc_->vector_atom;
        etai = visc_[ip];
      }

      if (vtanabs > 0) {
        // viscous wall force

        if (pairStyle == 1) { // pairstyle sph/artVisc/tensCorr

          itype = type[ip];
          csValues = cs->get_values();
          factor = mass * mass * etai * (1 + r0/r) * csValues[itype-1] /(rhoi * h4);

          // boundary functions for artificial viscosity:
          // Ftan and Fnorm are polynomial fits for the tangential and normal component of
          // [h^4 / vtan * sum(mu*gradW)] as function of wall distance, sum over boundary particles
          // mu is the [mu_ab = h*v_ab*r_ab/(r_ab^2+eta^2)] of artificial viscosity

          //if (kernel_id == 2) { // cubicspline
              // these polynomials were developed for the cubic spline kernel,
              // reasonable results with wendland and spiky kernel can be assumed
          if (q < 1.22) Ftan = 0.0794 * q5 - 0.23 * q4 - 0.05 * q3 + 0.83 * q2 - 1.03 * q + 0.407;
          else Ftan = 0;
          Ftan = Ftan * sl_corr;

          if (q < 1.31) Fnorm = 0.165 * q5 - 0.791 * q4 + 1.328 * q3 - 0.711 * q2 - 0.311 * q + 0.335;
          else Fnorm = 0;
          Fnorm = Fnorm * sl_corr;
          //} else {
          //  error->fix_error(FLERR,this,"wall potentials for chosen kernel not available");
          //}

          fwallShear = Ftan * factor;
          fwallRep = Fnorm * factor * vtanabs;

        } else if (pairStyle == 2) { // pairstyle sph/morris/tensCorr

          h5 = h4*sli;

          //if (kernel_id == 2) { // cubicspline
              // this polynomial was developed for the cubic spline kernel,
              // reasonable results with wendland and spiky kernel were obtained (tested)
          if (q < 1.43) Ftan = -0.571 * q4 + 1.413 * q3 + 0.622 * q2 - 3.92 * q + 2.59;
          else Ftan = 0;
          Ftan = Ftan * sl_corr;
          //} else {
          //  error->fix_error(FLERR,this,"wall potentials for chosen kernel not available");
          //}

          fwallShear = Ftan * mass * mass * 2 * etai /(rhoi * rhoi * h5) * (1 + r0/r);
          fwallRep = 0;
        } else {
          error->fix_error(FLERR,this,"wall potentials for chosen pairstyle not available");
        }
      }

      // boundary contribution for pressure:
      // F is a fit for sum(gradW) over boundary particles as function of wall distance

      //if (kernel_id == 2) { // cubicspline
          // this polynomial was developed for the cubic spline kernel,
          // reasonable results with wendland and spiky kernel were obtained (tested)
      F = 1.37 - 0.607 * q4 + 2.59 * q3 - 3.09 * q2 - 0.059 * q;
      F = F * sl_corr;
      //} else {
      //  error->fix_error(FLERR,this,"wall potentials for chosen kernel not available");
      //}

      F /= h4;
      fwallRep += mass * mass * 2 * p[ip]/(rhoi*rhoi) * F;

      // boundary contribution for density with sph/density/continuity
      if (densityStyle == 1) drho[ip] += mass * vnormalabs * F;

      // add repulsive penetration force
      if (r < r0 && r != 0) {
        frac = r0/r;
        fwallRep += D * (frac - 1);
      }

      // set wall force
      fwall[0] = fwallRep * dx;
      fwall[1] = fwallRep * dy;
      fwall[2] = fwallRep * dz;

      if (vtanabs > 0) {
        fwall[0] += fwallShear * vtan[0];
        fwall[1] += fwallShear * vtan[1];
        fwall[2] += fwallShear * vtan[2];
      }

      f[ip][0] += fwall[0];
      f[ip][1] += fwall[1];
      f[ip][2] += fwall[2];

      // store wall force for particle ip:
      wallForce_[ip][0] = fwall[0];
      wallForce_[ip][1] = fwall[1];
      wallForce_[ip][2] = fwall[2];
    }
}

/* ------------------------------------------------------------------------- */

void FixWallSphGeneralGap::compute_force_eval_gap(int ip,double mass,double r1,
                            double dx1,double dy1,double dz1,double *vwall1,double r2,
                            double dx2,double dy2,double dz2,double *vwall2)
{
    double **v = atom->vest;
    double *rho = atom->rho;
    double **f = atom->f;
    double vPrel[3],vWrel[3],gradP[3];
    double wallDistance, r_rel,contactArea,prefac;
    double TwallP[3],TwallShear[3],fwall1[3],fwall2[3];
    double fwallRep1,fwallRep2,frac;
    int  *type = atom->type;
    int itype;
    double etai;
    double sli,particlesPerWidth;

    // Get smoothing length sl:
    if (mass_type) {
      slType = fppaSlType->get_values();
      itype = type[ip];
      sli = slType[itype-1];
    } else {
      sl = fppaSl->vector_atom;
      sli = sl[ip];
    }

    if (modelStyle == 1) { // Newtonian
      etai = viscosity;
    } else {    // non-Newtonian
      visc_ = fix_visc_->vector_atom;
      etai = visc_[ip];
    }

    wallDistance = r1 + r2;
    r_rel = r1/wallDistance;
    contactArea = mass / (rho[ip] * wallDistance);
    particlesPerWidth = wallDistance * 1.2 / sli;

    // relative velocity of wall2:
    vWrel[0] = vwall2[0] - vwall1[0];
    vWrel[1] = vwall2[1] - vwall1[1];
    vWrel[2] = vwall2[2] - vwall1[2];

    // relative particle velocity:
    vPrel[0] = v[ip][0] - vwall1[0];
    vPrel[1] = v[ip][1] - vwall1[1];
    vPrel[2] = v[ip][2] - vwall1[2];

    // apparent pressure gradient:
    prefac = 2 * etai / (r1 * r2) * (1 + 0.5 /(particlesPerWidth * particlesPerWidth));
    gradP[0] = prefac * (vWrel[0] * r_rel - vPrel[0]);
    gradP[1] = prefac * (vWrel[1] * r_rel - vPrel[1]);
    gradP[2] = prefac * (vWrel[2] * r_rel - vPrel[2]);

    // Wall forces:
    TwallP[0] = gradP[0] * wallDistance/2;
    TwallP[1] = gradP[1] * wallDistance/2;
    TwallP[2] = gradP[2] * wallDistance/2;

    TwallShear[0] = vWrel[0] / wallDistance * etai;
    TwallShear[1] = vWrel[1] / wallDistance * etai;
    TwallShear[2] = vWrel[2] / wallDistance * etai;

    fwall1[0] = (TwallP[0] - TwallShear[0]) * contactArea;
    fwall1[1] = (TwallP[1] - TwallShear[1]) * contactArea;
    fwall1[2] = (TwallP[2] - TwallShear[2]) * contactArea;

    fwall2[0] = (TwallP[0] + TwallShear[0]) * contactArea;
    fwall2[1] = (TwallP[1] + TwallShear[1]) * contactArea;
    fwall2[2] = (TwallP[2] + TwallShear[2]) * contactArea;

    // add repulsive penetration forces
    if (r1 < r0 && r1 != 0) {
      frac = r0/r1;
      fwallRep1 = D * (frac - 1);

      fwall1[0] += fwallRep1 * dx1;
      fwall1[1] += fwallRep1 * dy1;
      fwall1[2] += fwallRep1 * dz1;
    }

    if (r2 < r0 && r2 != 0) {
      frac = r0/r2;
      fwallRep2 = D * (frac - 1);

      fwall2[0] += fwallRep2 * dx2;
      fwall2[1] += fwallRep2 * dy2;
      fwall2[2] += fwallRep2 * dz2;
    }

    // store wall forces for particle ip:
    wallForce_[ip][0] = fwall1[0];
    wallForce_[ip][1] = fwall1[1];
    wallForce_[ip][2] = fwall1[2];

    wallForce2_[ip][0] = fwall2[0];
    wallForce2_[ip][1] = fwall2[1];
    wallForce2_[ip][2] = fwall2[2];

    // set particle force (overwrite force calculated by pairstyle):
    fgradP_ = fix_fgradP_->array_atom;
    integrity_ = fix_integrity_->vector_atom;
    usedGapmodel_ = fix_usedGapmodel_->vector_atom;

    f[ip][0] = fgradP_[ip][0] / integrity_[ip] + fwall1[0] + fwall2[0];
    f[ip][1] = fgradP_[ip][1] / integrity_[ip] + fwall1[1] + fwall2[1];
    f[ip][2] = fgradP_[ip][2] / integrity_[ip] + fwall1[2] + fwall2[2];

    usedGapmodel_[ip] = 1;
    /*///////// unit vector in direction of relative wall2 velocity:
    eVw[0] = vWrel[0] / vWmag;
    eVw[1] = vWrel[1] / vWmag;
    eVw[2] = vWrel[2] / vWmag;

    // relative particle velocity:
    vPrel[0] = v[ip][0] - vwall1[0];
    vPrel[1] = v[ip][1] - vwall1[1];
    vPrel[2] = v[ip][2] - vwall1[2];

    // particle velocity parallel to wall2 velocity:
    vPparaMag = vPrel[0]*eVw[0] + vPrel[1]*eVw[1] + vPrel[2]*eVw[2];
    vPpara[0] = eVw[0]*vPparaMag;
    vPpara[1] = eVw[1]*vPparaMag;
    vPpara[2] = eVw[2]*vPparaMag;

    // particle velocity normal to wall1 (should be negligible):
    vPnormMag = vPrel[0]*dx1 + vPrel[1]*dy1 + vPrel[2]*dz1;
    vPnorm[0] = dx1*vPnormMag;
    vPnorm[1] = dx2*vPnormMag;
    vPnorm[2] = dx3*vPnormMag;

    // particle velocity in cross direction:
    vPcross[0] = vPrel[0] - vPpara[0] - vPnorm[0];
    vPcross[1] = vPrel[1] - vPpara[1] - vPnorm[1];
    vPcross[2] = vPrel[2] - vPpara[2] - vPnorm[2];
    vPcrossMag = sqrt(vPcross[0]*vPcross[0] + vPcross[1]*vPcross[1] + vPcross[2]*vPcross[2]);

    // unit vector in cross direction:
    eCr[0] = vPcross[0] / vPcrossMag;
    eCr[1] = vPcross[1] / vPcrossMag;
    eCr[2] = vPcross[2] / vPcrossMag;

    // pressure gradient in direction of wall velocity:
    gradPpara = (vWmag * r_rel - vPparaMag) * 2 * viscosity / (r1 * r2);

    // Wall force in direction of wall velocity;
    fwallShear = vWmag / wallDistance * viscosity * contactArea;
    fwallPpara = gradPpara * wallDistance/2 * contactArea;

    fwallPara1 = fwallPpara - fwallShear;
    fwallPara2 = fwallPpara + fwallShear;

    // pressure gradient in cross direction:
    gradPcross = - vPcrossMag * 2 * viscosity / (r1 * r2);

    // Wall force in cross direction:
    fwallCross = gradPcross * wallDistance/2 * contactArea;

    // sum up both wall forces:
    fwall1[0] = fwallPara1 * eVw[0] + fwallCross * eCr[0];
    fwall1[1] = fwallPara1 * eVw[1] + fwallCross * eCr[1];
    fwall1[2] = fwallPara1 * eVw[2] + fwallCross * eCr[2];

    fwall2[0] = fwallPara2 * eVw[0] + fwallCross * eCr[0];
    fwall2[1] = fwallPara2 * eVw[1] + fwallCross * eCr[1];
    fwall2[2] = fwallPara2 * eVw[2] + fwallCross * eCr[2];*/
}

