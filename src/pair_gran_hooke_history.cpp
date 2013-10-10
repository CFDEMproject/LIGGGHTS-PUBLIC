/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_history.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"
#include "compute_pair_gran_local.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::PairGranHookeHistory(LAMMPS *lmp) : PairGran(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum_pairgran = 3;

    Yeff = NULL;
    Geff = NULL;
    betaeff = NULL;
    veff = NULL;
    cohEnergyDens = NULL;
    coeffRestLog = NULL;
    coeffFrict = NULL;
    coeffRollFrict = NULL;
    coeffRollVisc = NULL;
    coeffMu = NULL;
    coeffRestMax = NULL;
    coeffStc = NULL;

    charVelflag = 1;

    force_off = false;
    sanity_checks = true;
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::~PairGranHookeHistory()
{
    memory->destroy(Yeff);
    memory->destroy(Geff);
    memory->destroy(betaeff);
    memory->destroy(veff);
    memory->destroy(cohEnergyDens);
    memory->destroy(coeffRestLog);
    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);
    memory->destroy(coeffRollVisc);

    memory->destroy(coeffMu);
    memory->destroy(coeffRestMax);
    memory->destroy(coeffStc);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value
    args[0] = (char *) "shearx";
    args[1] = (char *) "1";
    args[2] = (char *) "sheary";
    args[3] = (char *) "1";
    args[4] = (char *) "shearz";
    args[5] = (char *) "1";
    if (rollingflag == 2 || rollingflag == 3)
    {
      args[6] = (char *) "r_torquex_old"; 
      args[7] = (char *) "1";
      args[8] = (char *) "r_torquey_old"; 
      args[9] = (char *) "1";
      args[10] = (char *) "r_torquez_old"; 
      args[11] = (char *) "1";
    }
}

/* ---------------------------------------------------------------------- */

inline void PairGranHookeHistory::addCohesionForce(int &ip, int &jp,double &r, double &Fn_coh) 
{
    //r is the distance between the sphere's centeres
    double ri = atom->radius[ip];
    double rj = atom->radius[jp];
    double Acont;
    if(cohesionflag == 1)
     Acont = - M_PI/4 * ( (r-ri-rj)*(r+ri-rj)*(r-ri+rj)*(r+ri+rj) )/(r*r); //contact area of the two spheres
    else  Acont = M_PI * 2. * (2.*ri*rj/(ri+rj)) * (ri + rj - r);
    
    Fn_coh=cohEnergyDens[atom->type[ip]][atom->type[jp]]*Acont;
}

/* ---------------------------------------------------------------------- */

inline void PairGranHookeHistory::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu, double &vnnr) 
{
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double rj = atom->radius[jp];
    double ri = atom->radius[ip];
    double reff=ri*rj/(ri+rj);
    double stokes, coeffRestLogChosen;

    if (viscousflag)  {
       // Stokes Number from MW Schmeeckle (2001)
       stokes=meff*vnnr/(6.0*3.1416*coeffMu[itype][jtype]*reff*reff);
       // Empirical from Legendre (2006)
       coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
    } else {
       coeffRestLogChosen=coeffRestLog[itype][jtype];
    }

    kn = 16./15.*sqrt(reff)*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrt(reff)*Yeff[itype][jtype]),0.2);
    kt = kn;
    gamman=sqrt(4.*meff*kn/(1.+(M_PI/coeffRestLogChosen)*(M_PI/coeffRestLogChosen)));
    gammat=gamman;
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;

    return;
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::compute_force(int eflag, int vflag,int addflag)
{
    if     (rollingflag == 0) compute_force_eval<0>(eflag,vflag,addflag);
    else if(rollingflag == 1) compute_force_eval<1>(eflag,vflag,addflag);
    else if(rollingflag == 2) compute_force_eval<2>(eflag,vflag,addflag);
    else if(rollingflag == 3) compute_force_eval<3>(eflag,vflag,addflag);
}

/* ---------------------------------------------------------------------- */

template <int ROLLINGFRICTION>
void PairGranHookeHistory::compute_force_eval(int eflag, int vflag,int addflag)
{
  //calculated from the material properties 
  double kn,kt,kr,gamman,gammat,xmu,rmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3],dr_torque[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  double r_inertia,r_coef,r_torque_mag,r_torque_max,factor;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[dnum()*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;
        if (ROLLINGFRICTION == 2 || ROLLINGFRICTION == 3) {
          shear[3] = 0.0; // this is the r_torque_old
          shear[4] = 0.0; // this is the r_torque_old
          shear[5] = 0.0; // this is the r_torque_old
        }

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        double deltan=radsum-r;
        cri = radi-0.5*deltan;
        crj = radj-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          itype = type[i];
          jtype = type[j];
          mi = mass[itype];
          mj = mass[jtype];
        }
        if (fix_rigid)
        {
           if(body[i] >= 0) mi = masstotal[body[i]];
           if(body[j] >= 0) mj = masstotal[body[j]];
        }

        meff = mi*mj/(mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr);         

        // normal forces = Hookian contact + normal velocity damping

        damp = gamman*vnnr*rsqinv;  
        ccel = kn*(radsum-r)*rinv - damp;
        
        if (cohesionflag) { 
            addCohesionForce(i,j,r,Fn_coh);
            ccel-=Fn_coh*rinv;
        }

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

        shear = &allshear[dnum()*jj];

        if (shearupdate && computeflag)
        {
            shear[0] += vtr1*dt;
            shear[1] += vtr2*dt;
            shear[2] += vtr3*dt;

            // rotate shear displacements

            rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
            rsht *= rsqinv;
            shear[0] -= rsht*delx;
            shear[1] -= rsht*dely;
            shear[2] -= rsht*delz;
        }

        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +  shear[2]*shear[2]);

        // tangential forces = shear + tangential velocity damping

        fs1 = - (kt*shear[0]);
        fs2 = - (kt*shear[1]);
        fs3 = - (kt*shear[2]);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        // energy loss from sliding or damping
        if (fs > fn) {
            if (shrmag != 0.0) {
                fs1 *= fn/fs;
                fs2 *= fn/fs;
                fs3 *= fn/fs;
                shear[0] = -fs1/kt;
                shear[1] = -fs2/kt;
                shear[2] = -fs3/kt;
            }
            else fs1 = fs2 = fs3 = 0.0;
        }
        else
        {
            fs1 -= (gammat*vtr1);
            fs2 -= (gammat*vtr2);
            fs3 -= (gammat*vtr3);
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);

        // add rolling friction torque
        vectorZeroize3D(r_torque);
        if(ROLLINGFRICTION > 0)
        {
          if(ROLLINGFRICTION == 1)
          {
            vectorSubtract3D(omega[i],omega[j],wr_roll);
            wr_rollmag = vectorMag3D(wr_roll);

            if(wr_rollmag > 0.)
            {
                // calculate torque
                reff=radi*radj/(radi+radj);
                vectorScalarMult3D(wr_roll,rmu*kn*deltan*reff/wr_rollmag,r_torque);

                // remove normal (torsion) part of torque
                double rtorque_dot_delta = r_torque[0]*delx + r_torque[1]*dely + r_torque[2]*delz;
                r_torque_n[0] = delx * rtorque_dot_delta * rsqinv;
                r_torque_n[1] = dely * rtorque_dot_delta * rsqinv;
                r_torque_n[2] = delz * rtorque_dot_delta * rsqinv;
                vectorSubtract3D(r_torque,r_torque_n,r_torque);
            }
          }
          else // ROLLINGFRICTION == 2 || 3
          {
            double wr_roll_n[3],wr_roll_t[3];
            double r_inertia_red_i,r_inertia_red_j;

            itype = type[i];
            jtype = type[j];

            // relative rotational velocity
            vectorSubtract3D(omega[i],omega[j],wr_roll);

            // remove normal (torsion) part of relative rotation
            // use only tangential parts for rolling torque
            double wr_dot_delta = wr_roll[0]*delx+ wr_roll[1]*dely + wr_roll[2]*delz;
            wr_roll_n[0] = delx * wr_dot_delta * rsqinv;
            wr_roll_n[1] = dely * wr_dot_delta * rsqinv;
            wr_roll_n[2] = delz * wr_dot_delta * rsqinv;
            vectorSubtract3D(wr_roll,wr_roll_n,wr_roll_t);

            // spring
            reff=radi*radj/(radi+radj);
            if(ROLLINGFRICTION == 2)
              kr = 2.25*kn*rmu*rmu*reff*reff; 
            else
              kr = kt*reff*reff; 

            vectorScalarMult3D(wr_roll_t,update->dt*kr,dr_torque); 

            r_torque[0] = shear[3] + dr_torque[0];
            r_torque[1] = shear[4] + dr_torque[1];
            r_torque[2] = shear[5] + dr_torque[2];

            // limit max. torque
            r_torque_mag = vectorMag3D(r_torque);
            r_torque_max = fabs(ccel*r)*reff*rmu;
            if(r_torque_mag > r_torque_max)
            {
              factor = r_torque_max / r_torque_mag;

              r_torque[0] *= factor;
              r_torque[1] *= factor;
              r_torque[2] *= factor;

              // save rolling torque due to spring
              shear[3] = r_torque[0];
              shear[4] = r_torque[1];
              shear[5] = r_torque[2];

              // no damping / no dashpot in case of full mobilisation rolling angle
              // r_coef = 0.0;

            } else {

              // save rolling torque due to spring before adding damping torque
              shear[3] = r_torque[0];
              shear[4] = r_torque[1];
              shear[5] = r_torque[2];

              // dashpot only for the original epsd model
              if(ROLLINGFRICTION == 2)
              {
                // dashpot
                r_inertia_red_i = mi*radi*radi;
                r_inertia_red_j = mj*radj*radj;
                if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
                else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);

                r_coef = coeffRollVisc[itype][jtype] * 2 * sqrt(r_inertia*kr);

                // add damping torque
                r_torque[0] += r_coef*wr_roll_t[0];
                r_torque[1] += r_coef*wr_roll_t[1];
                r_torque[2] += r_coef*wr_roll_t[2];
              }
            }

          }

        }

        if(computeflag)
        {
            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;
            torque[i][0] -= cri*tor1 + r_torque[0];
            torque[i][1] -= cri*tor2 + r_torque[1];
            torque[i][2] -= cri*tor3 + r_torque[2];
        }

        if (j < nlocal && computeflag) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= crj*tor1 - r_torque[0];
          torque[j][1] -= crj*tor2 - r_torque[1];
          torque[j][2] -= crj*tor3 - r_torque[2];
        }

        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeHistory::settings(int narg, char **arg) 
{
    iarg_ = 0;

    // set defaults
    dampflag = 1;
    rollingflag = 0;
    cohesionflag = 0;
    viscousflag = 0;
    force_off = false;

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"cohesion") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'cohesion'");
            iarg_++;
            if(strcmp(arg[iarg_],"sjkr") == 0)
                cohesionflag = 1;
            else if(strcmp(arg[iarg_],"sjkr2") == 0)
                cohesionflag = 2;
            else if(strcmp(arg[iarg_],"off") == 0)
                cohesionflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'sjkr' or 'off' after keyword 'cohesion'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"force") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'force'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                force_off = false;
            else if(strcmp(arg[iarg_],"off") == 0)
                force_off = true;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'force'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"sanity_checks") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'sanity_checks'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                sanity_checks = true;
            else if(strcmp(arg[iarg_],"off") == 0)
                sanity_checks = false;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'sanity_checks'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"rolling_friction") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'rolling_friction'");
            iarg_++;
            if(strcmp(arg[iarg_],"cdt") == 0)
                rollingflag = 1;
            else if(strcmp(arg[iarg_],"epsd") == 0) {
                rollingflag = 2;
                dnum_pairgran = 6; 
            } else if(strcmp(arg[iarg_],"epsd2") == 0) {
                rollingflag = 3;
                dnum_pairgran = 6;
            } else if(strcmp(arg[iarg_],"off") == 0)
                rollingflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'cdt', 'epsd', 'epsd2' or 'off' after keyword 'rolling_friction'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"tangential_damping") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'tangential_damping'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                dampflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                dampflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'tangential_damping'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"viscous") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'viscous'");
            iarg_++;
            if(strcmp(arg[iarg_],"stokes") == 0)
                viscousflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                viscousflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'stokes' or 'off' after keyword 'viscous'");
            iarg_++;
            hasargs = true;
        } else if (force->pair_match("gran/hooke/history",1) || force->pair_match("gran/hertz/history",1))
            error->all(FLERR,"Illegal pair_style gran command, illegal keyword");
    }

    if(cohesionflag && domain->dimension!=3)
        error->all(FLERR,"Cohesion model valid for 3d simulations only");
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistory::init_granular()
{
  int max_type = mpg->max_type();

  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties
  
  Y1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,force->pair_style));
  v1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,force->pair_style));

  coeffRest1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  if(rollingflag == 2) // damping for original epsd model only
    coeffRollVisc1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingViscousDamping","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  if(viscousflag)
  {
    coeffMu1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("FluidViscosity","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    coeffRestMax1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("MaximumRestitution","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    coeffStc1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("CriticalStokes","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }

  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(charVelflag)
      charVel1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,force->pair_style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          double Yi=Y1->compute_vector(i-1);
          double Yj=Y1->compute_vector(j-1);
          double vi=v1->compute_vector(i-1);
          double vj=v1->compute_vector(j-1);
          double cor = coeffRest1->compute_array(i-1,j-1);

          // error checks on Y, v, e

          if(sanity_checks)
          {
              if(strcmp(update->unit_style,"si") == 0  && Yi < 5e6)
                 error->all(FLERR,"youngsModulus >= 5e6 required for SI units");
              if(strcmp(update->unit_style,"cgs") == 0 && Yi < 5e5)
                error->all(FLERR,"youngsModulus >= 5e5 required for CGS units");

              if(vi < 0. || vi > 0.5)
                error->all(FLERR,"0 <= poissonsRatio <= 0.5 required");

              if(cor <= 0.05 || cor > 1)
                error->all(FLERR,"0.05 < coefficientRestitution <= 1 required");
          }

          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);

          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

         if(viscousflag)
         {
           coeffMu[i][j] = coeffMu1->compute_array(i-1,j-1);
           coeffRestMax[i][j] = coeffRestMax1->compute_array(i-1,j-1);
           coeffStc[i][j] = coeffStc1->compute_array(i-1,j-1);

           // error check
           if(sanity_checks)
           {
               if(coeffRestMax[i][j] <= 0. || coeffRestMax[i][j] > 1)
                 error->all(FLERR,"0 < MaximumRestitution <= 1 required");
               if(coeffMu[i][j] <= 0.)
                 error->all(FLERR,"coeffMu > 0 required");
               if(coeffStc[i][j] <= 0.)
                 error->all(FLERR,"CriticalStokes > 0 required");
           }
         }

          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);
          if(rollingflag == 2) coeffRollVisc[i][j] = coeffRollVisc1->compute_array(i-1,j-1);
          
          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
          //omitting veff here

      }
  }

  if(charVelflag)
  {
      charVel = charVel1->compute_scalar();
      if(sanity_checks)
      {
            if(strcmp(update->unit_style,"si") == 0  && charVel < 1e-2)
                 error->all(FLERR,"characteristicVelocity >= 1e-2 required for SI units");
      }
  }

  // error checks on coarsegraining
  if((rollingflag || cohesionflag) && force->cg_active())
    error->cg(FLERR,"Granular model with rolling friction and / or cohesion");

  // error checks on coarsegraining
  if((rollingflag || cohesionflag) && force->cg_active())
    error->cg(FLERR,"Granular model with rolling friction and / or cohesion");
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeHistory::allocate_properties(int size)
{
    memory->destroy(Yeff);
    memory->destroy(Geff);
    memory->destroy(betaeff);
    memory->destroy(veff);
    memory->destroy(cohEnergyDens);
    memory->destroy(coeffRestLog);
    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);
    memory->destroy(coeffRollVisc);
    memory->destroy(coeffMu);
    memory->destroy(coeffRestMax);
    memory->destroy(coeffStc);

    memory->create(Yeff,size+1,size+1,"Yeff");
    memory->create(Geff,size+1,size+1,"Geff");
    memory->create(betaeff,size+1,size+1,"betaeff");
    memory->create(veff,size+1,size+1,"veff");
    memory->create(cohEnergyDens,size+1,size+1,"cohEnergyDens");
    memory->create(coeffRestLog,size+1,size+1,"coeffRestLog");
    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
    memory->create(coeffRollFrict,size+1,size+1,"coeffRollFrict");
    memory->create(coeffRollVisc,size+1,size+1,"coeffRollVisc");
    memory->create(coeffMu,size+1,size+1,"coeffMu");
    memory->create(coeffRestMax,size+1,size+1,"coeffRestMax");
    memory->create(coeffStc,size+1,size+1,"coeffStc");

}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairGranHookeHistory::write_restart_settings(FILE *fp) 
{
  int writeflag = dampflag + rollingflag * 2;
  fwrite(&writeflag,sizeof(int),1,fp);
  fwrite(&cohesionflag,sizeof(int),1,fp);
  fwrite(&viscousflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts 
------------------------------------------------------------------------- */

void PairGranHookeHistory::read_restart_settings(FILE *fp) 
{
  if (comm->me == 0) {
    int readflag;
    fread(&readflag,sizeof(int),1,fp);
    fread(&cohesionflag,sizeof(int),1,fp);
    fread(&viscousflag,sizeof(int),1,fp);
    dampflag = readflag & 1;
    rollingflag = readflag & 2;
  }
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
  MPI_Bcast(&cohesionflag,1,MPI_INT,0,world);
  MPI_Bcast(&rollingflag,1,MPI_INT,0,world);
  MPI_Bcast(&viscousflag,1,MPI_INT,0,world);
}
