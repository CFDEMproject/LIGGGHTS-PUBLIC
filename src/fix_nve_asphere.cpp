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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h" 

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA_SPHEROID 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEAsphere::FixNVEAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
  updateRotation_(false), 
  fix_orientation_(NULL), 
  fix_shape_(NULL)
{ 
  if (narg < 3) error->all(FLERR,"Illegal fix nve/asphere command");

  // process extra keywords

  int iarg = 3;
  while (iarg < narg) { 
    if (strcmp(arg[iarg],"updateRotation") == 0)
    {
      updateRotation_=true;
      //printf("nve/asphere will update the rotation rate and orientation vector ex!\n");
      iarg += 1;
    } else error->all(FLERR,"Illegal fix nve/asphere command");
  }
}

/* ---------------------------------------------------------------------- */
void FixNVEAsphere::post_create()
{

   if(updateRotation_) 
   {
     const char *fixarg[11];

     if(fix_orientation_ == NULL)
     {
        //register orientation as property/atom
        fixarg[0]="ex";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="ex";
        fixarg[4]="vector";
        fixarg[5]="yes";    // restart
        fixarg[6]="yes";    // communicate ghost forward
        fixarg[7]="no";     // communicate ghost reverse
        fixarg[8]="1";
        fixarg[9]="0";
        fixarg[10]="0";
        modify->add_fix(11,const_cast<char**>(fixarg));
        fix_orientation_=static_cast<FixPropertyAtom*>(modify->find_fix_property("ex","property/atom","vector",0,0,style));
     }
     if(fix_shape_ == NULL)
     {
        //register orientation as property/atom
        fixarg[0]="shape";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="shape";
        fixarg[4]="vector";
        fixarg[5]="no";     // restart
        fixarg[6]="yes";    // communicate ghost forward
        fixarg[7]="no";     // communicate ghost reverse
        fixarg[8]="0";
        fixarg[9]="0";
        fixarg[10]="0";
        modify->add_fix(11,const_cast<char**>(fixarg));
        fix_shape_     = static_cast<FixPropertyAtom*>(modify->find_fix_property("shape","property/atom","vector",0,0,style));
     }
   }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute nve/asphere requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/asphere requires extended particles");

  FixNVE::init();
  fix_orientation_->do_forward_comm();
  fix_shape_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::initial_integrate(int vflag)
{
  double dtfm;
  double inertia[3],omega[3];
  double exone[3],eyone[3],ezone[3];  
  double *shape,*quat;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  //save rotation rate to array if necessary
  
  double **omegaParticles = NULL;
  double **orientation = NULL;
  double **shapeFix       = NULL;
  if(updateRotation_)
  {
      omegaParticles = atom->omega;
      if(fix_orientation_)
          orientation = fix_orientation_->array_atom;
      if(fix_shape_)
          shapeFix    = fix_shape_->array_atom;
  }

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      // principal moments of inertia

      shape = bonus[ellipsoid[i]].shape;
      quat = bonus[ellipsoid[i]].quat;

      if( (shape[0]<shape[1]) || (shape[0]<shape[2]) )
        error->one(FLERR,"Shape is not correctly specified. shape[0] must be the largest value!");

      //Moment of inertia in the Principal coordinate system (denoted as 'prime'), see http://mathworld.wolfram.com/Spheroid.html
      inertia[0] = INERTIA_SPHEROID*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = INERTIA_SPHEROID*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = INERTIA_SPHEROID*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);

      if(updateRotation_) 
      {
            omegaParticles[i][0]=omega[0];
            omegaParticles[i][1]=omega[1];
            omegaParticles[i][2]=omega[2];
            if(fix_orientation_)
            {
/*               g[0] = orientation[i][0] + dtv * (omega[1]*orientation[i][2]-omega[2]*orientation[i][1]);
               g[1] = orientation[i][1] + dtv * (omega[2]*orientation[i][0]-omega[0]*orientation[i][2]);
               g[2] = orientation[i][2] + dtv * (omega[0]*orientation[i][1]-omega[1]*orientation[i][0]);
               msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
               scale = 1.0/sqrt(msq);
               orientation[i][0] = g[0]*scale;
               orientation[i][1] = g[1]*scale;
               orientation[i][2] = g[2]*scale; */

               //Alternative calculation
               MathExtra::q_to_exyz(quat,exone,eyone,ezone);   
               orientation[i][0] = exone[0];
               orientation[i][1] = exone[1];
               orientation[i][2] = exone[2];
//               printf("exone[0]: %g,exone[1]: %g,exone[2]: %g. \n",
//                       exone[0],exone[1],exone[2]);
            }
            if(fix_shape_) 
            {
                shapeFix[i][0] = shape[0];
                shapeFix[i][1] = shape[1];
                shapeFix[i][2] = shape[2];
      		}
      }
    }
    fix_orientation_->do_forward_comm();
    fix_shape_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  AtomVecEllipsoid::Bonus *bonus = NULL;
  int *ellipsoid = NULL;
  double **omegaParticles = NULL;
  if(updateRotation_)
  {
      omegaParticles = atom->omega;
      ellipsoid = atom->ellipsoid;
      bonus = avec->bonus;
  }
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      if(updateRotation_) 
      {
            double inertia[3],omega[3];
            double *shape,*quat;

            shape = bonus[ellipsoid[i]].shape;
            quat = bonus[ellipsoid[i]].quat;

            inertia[0] = INERTIA_SPHEROID*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
            inertia[1] = INERTIA_SPHEROID*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
            inertia[2] = INERTIA_SPHEROID*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
            MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);

            omegaParticles[i][0]=omega[0];
            omegaParticles[i][1]=omega[1];
            omegaParticles[i][2]=omega[2];
      }
    }
}
