/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_nve_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h" //NP modified TUG

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEAsphere::FixNVEAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
  updateRotation_(false), //NP modified TUG
  fix_orientation_(0) //NP modified TUG
{ //NP modified TUG
  if (narg < 3) error->all(FLERR,"Illegal fix nve/asphere command");

  // process extra keywords

  int iarg = 3;
  while (iarg < narg) { //NP modified TUG
    if (strcmp(arg[iarg],"updateRotation") == 0) 
    {
      updateRotation_=true;
      //printf("nve/asphere will update the rotation rate and orientation vector ex!\n");
      iarg += 1;
    } else error->all(FLERR,"Illegal fix nve/asphere command");
  }

}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute nve/asphere requires atom style ellipsoid");

   if(updateRotation_) //NP modified TUG
   {
     fix_orientation_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("ex","property/atom","vector",0,0,style));
     if(fix_orientation_)
          printf("FOUND fix_orientation_!!\n");
   }

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
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::initial_integrate(int vflag)
{
  double dtfm;
  double inertia[3],omega[3];
  double g[3], msq, scale; //NP modified TUG
  double exone[3],eyone[3],ezone[3];  //NP modified TUG
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
  //NP modified TUG
  double **omegaParticles = NULL;
  double **orientation = NULL;
  if(updateRotation_)
  {
      omegaParticles = atom->omega;
      if(fix_orientation_)
          orientation = fix_orientation_->array_atom;
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

      inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);

      if(updateRotation_) //NP modified TUG
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
               MathExtra::q_to_exyz(quat,exone,eyone,ezone);   //NP modified TUG
               orientation[i][0] = exone[0];
               orientation[i][1] = exone[1];
               orientation[i][2] = exone[2];
//               printf("exone[0]: %g,exone[1]: %g,exone[2]: %g. \n",
//                       exone[0],exone[1],exone[2]);
            }
      }
    }
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

  //NP modified TUG begin
  AtomVecEllipsoid::Bonus *bonus = NULL;
  int *ellipsoid = NULL;
  double **omegaParticles = NULL;
  if(updateRotation_)
  {
      omegaParticles = atom->omega;
      ellipsoid = atom->ellipsoid;
      bonus = avec->bonus;
  }
  //NP modified TUG end

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

      if(updateRotation_) //NP modified TUG
      {
            double inertia[3],omega[3];
            double *shape,*quat;

            shape = bonus[ellipsoid[i]].shape;
            quat = bonus[ellipsoid[i]].quat;

            inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
            inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
            inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
            MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);

            omegaParticles[i][0]=omega[0];
            omegaParticles[i][1]=omega[1];
            omegaParticles[i][2]=omega[2];
      }
    }
}
