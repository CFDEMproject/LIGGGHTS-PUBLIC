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
#include "string.h"
#include "pair_gran_hooke.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "neigh_list.h"
#include "error.h"
#include "vector_liggghts.h"
#include "fix_rigid.h"
#include "compute_pair_gran_local.h"

using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

PairGranHooke::PairGranHooke(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
    no_virial_fdotr_compute = 0;

    //flag that we do not intend to use contact history
    history = 0;
    dnum_pairgran = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHooke::settings(int narg, char **arg) 
{
    PairGranHookeHistory::settings(narg,arg);

    if(rollingflag == 2)
        if (narg < iarg_+2) error->all(FLERR,"Cannot use 'epsd' with pair gran/hooke");
}

/* ---------------------------------------------------------------------- */

void PairGranHooke::compute_force(int eflag, int vflag,int addflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu,rmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3,cri,crj;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,ft,fs1,fs2,fs3;
  int *ilist,*jlist,*numneigh,**firstneigh;

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
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
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

      if (rsq < radsum*radsum) {
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
          mi=rmass[i];
          mj=rmass[j];
        } else {
          itype = type[i];
          jtype = type[j];
          mi=mass[itype];
          mj=mass[jtype];
        }
        if (fix_rigid)
        {
           if(body[i] >= 0) mi = masstotal[body[i]];
           if(body[j] >= 0) mj = masstotal[body[j]];
        }
        meff=mi*mj/(mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr);         

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

        // force normalization

        fn = xmu * fabs(ccel*r);
        fs = gammat*vrel;     
        if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
        else ft = 0.0;

        // tangential force due to tangential velocity damping

        fs1 = -ft*vtr1;
        fs2 = -ft*vtr2;
        fs3 = -ft*vtr3;

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

        if (computeflag && (newton_pair || j < nlocal)) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= crj*tor1 - r_torque[0];
          torque[j][1] -= crj*tor2 - r_torque[1];
          torque[j][2] -= crj*tor3 - r_torque[2];
        }
        
        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,NULL);

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair, 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}
