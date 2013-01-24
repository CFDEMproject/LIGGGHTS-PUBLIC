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
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "pair_gran_hertz_history.h"
#include "atom.h"
#include "force.h"
#include "fix_property_global.h"
#include "error.h"
#include "modify.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

#define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875

/* ---------------------------------------------------------------------- */

PairGranHertzHistory::PairGranHertzHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp)
{
    charVelflag = 0;
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */

inline void PairGranHertzHistory::deriveContactModelParams(int &ip, int &jp,double &meff, double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &vnnr) 
{
    
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double ri = atom->radius[ip];
    double rj = atom->radius[jp];
    double reff=ri*rj/(ri+rj);

    double sqrtval = sqrt(reff*deltan);

    double Sn=2.*Yeff[itype][jtype]*sqrtval;
    double St=8.*Geff[itype][jtype]*sqrtval;

    kn=4./3.*Yeff[itype][jtype]*sqrtval;
    kt=St;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
    gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    /* NP
    char *testmsg=new char[200];
    sprintf(testmsg,"Yeff=%f,reff=%f,deltan=%f, kn=%f, kt=%f, gamman=%f, gammat=%f, xmu=%f\n",Yeff, reff,deltan, kn,kt,gamman,gammat,xmu);
    error->warning(testmsg);
    delete []testmsg;*/

    return;
}
