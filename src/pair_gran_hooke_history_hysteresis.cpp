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
#include "pair_gran_hooke_history_hysteresis.h"
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

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryHysteresis::PairGranHookeHistoryHysteresis(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum_pairgran = 4;  

    kn2k2Max_ = NULL;
    kn2kc_ = NULL;
    phiF_ = NULL;

}

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryHysteresis::~PairGranHookeHistoryHysteresis()
{
  memory->destroy(kn2k2Max_);
  memory->destroy(kn2kc_);
  memory->destroy(phiF_);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryHysteresis::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value
    args[0] = (char *) "shearx";
    args[1] = (char *) "1";
    args[2] = (char *) "sheary";
    args[3] = (char *) "1";
    args[4] = (char *) "shearz";
    args[5] = (char *) "1";
    args[6] = (char *) "deltaMax"; 
    args[7] = (char *) "1";
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistoryHysteresis::compute_force(int eflag, int vflag,int addflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu,rmu; 

  double fHys,fTmp; 
  double k2,k2Max,kc,deltaMax,deltaMaxLim; 

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

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
        shear[3] = 0.0; 

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

        double mi,mj;
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

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr);         //modified C.K

        damp = gamman*vnnr*rsqinv;  

        itype = type[i];
        jtype = type[j];

        shear = &allshear[dnum()*jj];
        deltaMax = shear[3]; // the 4th value of the history array is deltaMax
        if (deltan > deltaMax) {
          shear[3] = deltan;
          deltaMax = deltan;
        }

        k2Max = kn * kn2k2Max_[itype][jtype]; 
        kc = kn * kn2kc_[itype][jtype]; 

        // k2 dependent on the maximum overlap
        // this accounts for an increasing stiffness with deformation
        deltaMaxLim =(k2Max/(k2Max-kn))*phiF_[itype][jtype]*2*radi*radj/(radi+radj);
        if (deltaMax >= deltaMaxLim) k2 = k2Max;
        else k2 = kn+(k2Max-kn)*deltaMax/deltaMaxLim;

        fTmp = k2*(deltan-deltaMax)+kn*deltaMax;//k2*(deltan-delta0);
        if (fTmp >= kn*deltan) {
          fHys = kn*deltan;
        } else {
          if (fTmp > -kc*deltan) {
            fHys = fTmp;
          } else fHys = -kc*deltan;
        }

        ccel = fHys*rinv - damp; //kn*(radsum-r)*rinv - damp;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

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
        if(rollingflag)
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

void PairGranHookeHistoryHysteresis::settings(int narg, char **arg) 
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
        if (strcmp(arg[iarg_],"force") == 0) {
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
        } else if (strcmp(arg[iarg_],"rolling_friction") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'rolling_friction'");
            iarg_++;
            if(strcmp(arg[iarg_],"cdt") == 0)
                rollingflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                rollingflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'cdt' or 'off' after keyword 'rolling_friction'");
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
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'dampflag'");
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

}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistoryHysteresis::init_granular()
{
  int max_type = mpg->max_type();

  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties
  
  Y1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,force->pair_style));
  v1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,force->pair_style));

  coeffRest1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  coeffKn2K2Max_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientMaxElasticStiffness","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  coeffKn2Kc_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientAdhesionStiffness","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  coeffPhiF_=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientPlasticityDepth","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(charVelflag) charVel1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,force->pair_style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          double Yi=Y1->compute_vector(i-1);
          double Yj=Y1->compute_vector(j-1);
          double vi=v1->compute_vector(i-1);
          double vj=v1->compute_vector(j-1);

          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);

          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

          kn2k2Max_[i][j] = coeffKn2K2Max_->compute_array(i-1,j-1);
          kn2kc_[i][j] = coeffKn2Kc_->compute_array(i-1,j-1);
          phiF_[i][j] = coeffPhiF_->compute_array(i-1,j-1);

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          //omitting veff here

      }
  }

  if(charVelflag) charVel = charVel1->compute_scalar();

  // error checks on coarsegraining
  if((rollingflag || cohesionflag) && force->cg_active())
    error->cg(FLERR,"Granular model with rolling friction and / or cohesion");
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeHistoryHysteresis::allocate_properties(int size)
{
    memory->destroy(Yeff);
    memory->destroy(Geff);
    memory->destroy(betaeff);
    memory->destroy(veff);
    memory->destroy(cohEnergyDens); 
    memory->destroy(coeffRestLog);
    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);
    memory->destroy(kn2k2Max_);
    memory->destroy(kn2kc_);
    memory->destroy(phiF_);
    memory->create(Yeff,size+1,size+1,"Yeff");
    memory->create(Geff,size+1,size+1,"Geff");
    memory->create(betaeff,size+1,size+1,"betaeff");
    memory->create(veff,size+1,size+1,"veff");
    memory->create(cohEnergyDens,size+1,size+1,"cohEnergyDens");
    memory->create(coeffRestLog,size+1,size+1,"coeffRestLog");
    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
    memory->create(coeffRollFrict,size+1,size+1,"coeffRollFrict");
    memory->create(kn2k2Max_,size+1,size+1,"kn2k2Max_");
    memory->create(kn2kc_,size+1,size+1,"kn2kc_");
    memory->create(phiF_,size+1,size+1,"phiF_");
}
