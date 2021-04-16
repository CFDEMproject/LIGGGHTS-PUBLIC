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
#include "fix_wall_sph_general.h"
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

FixWallSphGeneral::FixWallSphGeneral(LAMMPS *lmp, int narg, char **arg) :
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

    if (mass_type) {
      fppaSlType=static_cast<FixPropertyGlobal*>(modify->find_fix_property("sl","property/global","peratomtype",ntypes,0,force->pair_style));
    } else {
      fppaSl=static_cast<FixPropertyAtom*>(modify->find_fix_property("sl","property/atom","scalar",0,0,force->pair_style,false));
    }

    cs=static_cast<FixPropertyGlobal*>(modify->find_fix_property("speedOfSound","property/global","peratomtype",ntypes,0,force->pair_style));
}

/* ---------------------------------------------------------------------- */

FixWallSphGeneral::~FixWallSphGeneral()
{

}

/* ---------------------------------------------------------------------- */

void FixWallSphGeneral::compute_density(int ip,double r,double mass)
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

void FixWallSphGeneral::compute_velgrad(int ip,double delx,double dely, double delz,double mass,double *vwall)
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

void FixWallSphGeneral::compute_force(SurfacesIntersectData & sidata, double *vwall)
{
    // dx, dy, dz is normalized in case SPH
    const int ip = sidata.i;
    const double r = sidata.r;
    const double dx = sidata.delta[0];
    const double dy = sidata.delta[1];
    const double dz = sidata.delta[2];
    const double mass = sidata.mi;

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

      if (StaticWall==1) {
        vwall[0] = vwallX;
        vwall[1] = vwallY;
        vwall[2] = vwallZ;
      }

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
