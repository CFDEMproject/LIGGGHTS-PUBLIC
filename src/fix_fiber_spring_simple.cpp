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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2013-     TU Graz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include "force.h"
#include <stdlib.h>
#include <string.h>
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"
#include "fix_fiber_spring_simple.h"
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "pair_gran.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "properties.h"
#include "respa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

/* ---------------------------------------------------------------------- */

FixFiberSpringSimple::FixFiberSpringSimple(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix fiber/simpleSpring command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  group2 = NULL;

  Yeff = NULL;

  if (strcmp(arg[3],"tether") == 0) {
    if (narg != 9) error->all(FLERR,"Illegal fix spring command");
    k_spring = atof(arg[4]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[5],"NULL") == 0) xflag = 0;
    else xc = atof(arg[5]);
    if (strcmp(arg[6],"NULL") == 0) yflag = 0;
    else yc = atof(arg[6]);
    if (strcmp(arg[7],"NULL") == 0) zflag = 0;
    else zc = atof(arg[7]);
    r0 = atof(arg[8]);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix spring command");

  } else if (strcmp(arg[3],"couple") == 0) {
    if (narg != 10) error->all(FLERR,"Illegal fix spring command");

    int n = strlen(arg[4]) + 1;
    group2 = new char[n];
    strcpy(group2,arg[4]);
    igroup2 = group->find(arg[4]);
    if (igroup2 == -1)
      error->all(FLERR,"Fix spring couple group ID does not exist");
    if (igroup2 == igroup)
      error->all(FLERR,"Two groups cannot be the same in fix spring couple");
    group2bit = group->bitmask[igroup2];

    k_spring = atof(arg[5]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[6],"NULL") == 0) xflag = 0;
    else xc = atof(arg[6]);
    if (strcmp(arg[7],"NULL") == 0) yflag = 0;
    else yc = atof(arg[7]);
    if (strcmp(arg[8],"NULL") == 0) zflag = 0;
    else zc = atof(arg[8]);
    r0 = atof(arg[9]);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix spring command");

  } else error->all(FLERR,"Illegal fix spring command");

  ftotal[0] = ftotal[1] = ftotal[2] = ftotal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixFiberSpringSimple::~FixFiberSpringSimple()
{
  delete [] group2;
  memory->destroy(Yeff);
}

/* ---------------------------------------------------------------------- */

int FixFiberSpringSimple::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::init()
{
  // recheck that group 2 has not been deleted

  if (group2) {
    igroup2 = group->find(group2);
    if (igroup2 == -1)
      error->all(FLERR,"Fix spring couple group ID does not exist");
    group2bit = group->bitmask[igroup2];
  }

  //Get pointer to the fixes that have the material properties
  int max_type = 1;
  memory->destroy(Yeff);
  memory->create(Yeff,max_type+1,max_type+1,"Yeff");

  Y1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,force->pair_style));
  v1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,force->pair_style));
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
          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);

      }
  }
  charVel = charVel1->compute_scalar();

  if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix fiber/simpleSpring");
  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

  masstotal = group->mass(igroup);
//  if (styleflag == COUPLE) masstotal2 = group->mass(igroup2);

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::post_force(int vflag)
{

  //calculated from the material properties
  double deltaEquilibrium= 0.75;
  double snapDistance    = 0.90;
  double kn = 1000;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv;
  double meff,ccel;
  int *ilist,*jlist,*numneigh,**firstneigh;
  //int *contact_flag,**first_contact_flag;
  //double *all_contact_hist,**first_contact_hist;

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;

  //**************  LOOP OVER PARTICLES *********
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    //contact_flag = first_contact_flag[i];
    //all_contact_hist = first_contact_hist[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++)
    {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0]; //delta vector pointing from j to i
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);

      radj = radius[j];
      radsum = radi + radj;

      if (r <= radsum*snapDistance) //check if there is sufficient overlap
      {
//        fprintf(screen, "SNAPPED!! \n");

        //double deltan=radsum-r;
        rinv = 1.0/r;

        double mi,mj;
        if (rmass)
        {
          mi = rmass[i];
          mj = rmass[j];
        }
        else
        {
          itype = type[i];
          jtype = type[j];
          mi = mass[itype];
          mj = mass[jtype];
        }

        meff = mi*mj/(mi+mj);
        //deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu);
        deriveNormalStiffness(i,j,meff,kn);

        //Calculate the force and subtract to model bonding
        ccel = kn * radsum*(1.0-deltaEquilibrium) * rinv ;
//        ccel = kn*(radsum-r)*rinv;
        f[i][0] -= delx*ccel;
        f[i][1] -= dely*ccel;
        f[i][2] -= delz*ccel;

        if (j < nlocal)
        {
          f[j][0] += delx*ccel;
          f[j][1] += dely*ccel;
          f[j][2] += delz*ccel;
        }

//        fprintf(screen, "r: %3.3g, critDistance: %3.3g \n", r, radsum*snapDistance);
//        fprintf(screen, "NormalStiffness: %3.3g \n", kn);
//        fprintf(screen, "i: %2.1d f: %3.3g %3.3g %3.3g \n", i,  f[i][0],f[i][1],f[i][2]);

      }

    }
   }

}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFiberSpringSimple::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixFiberSpringSimple::compute_scalar()
{
  return espring;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixFiberSpringSimple::compute_vector(int n)
{
  return ftotal[n];
}

inline void FixFiberSpringSimple::deriveNormalStiffness(int &ip, int &jp,double &meff, double &kn)
{

    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double rj = atom->radius[jp];
    double ri = atom->radius[ip];
    double reff=ri*rj/(ri+rj);

    kn = 16./15.*sqrt(reff)*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrt(reff)*Yeff[itype][jtype]),0.2);

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;

    return;

}

