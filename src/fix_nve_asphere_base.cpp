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
    Alexander Podlozhnyuk, DCS Computing GmbH, Linz
    Christoph Kloss, DCS Computing GmbH, Linz

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <string.h>
#include "fix_nve_asphere_base.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"
#include "math_extra_liggghts_nonspherical.h"
#include "fix_property_atom.h"

/* ---------------------------------------------------------------------- */

FixNVEAsphereBase::FixNVEAsphereBase(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/superquadric command");

  time_integrate = 1;

  // process extra keywords

  integration_scheme = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"integration_scheme") == 0) {
      integration_scheme = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }
    else error->fix_error(FLERR,this,"unknown keyword");
  }

  hdtorque = NULL;
  orientation = NULL;
  ksl_rotation = NULL;
  couple_fix_id = -1;
  for(int ifix = 0; ifix < modify->nfix; ifix++) {
    if(strcmp(modify->fix[ifix]->style, "couple/cfd/force/implicit")==0) {
      couple_fix_id = ifix;
      break;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::init()
{
  FixNVE::init();

  // error checks might go here
}

/* ---------------------------------------------------------------------- */

//dynamic Euler equations for angular velocity in body's principal axes
void FixNVEAsphereBase::dynamic_euler(double *wbody, double *tbody, double *inertia, double *result)
{
  result[0] = tbody[0] / inertia[0] + wbody[1]*wbody[2]*((inertia[1] - inertia[2]) / inertia[0]);
  result[1] = tbody[1] / inertia[1] + wbody[2]*wbody[0]*((inertia[2] - inertia[0]) / inertia[1]);
  result[2] = tbody[2] / inertia[2] + wbody[0]*wbody[1]*((inertia[0] - inertia[1]) / inertia[2]);
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::integrate_dynamic_euler(double dt, double *wbody, double *tbody, double *inertia)
{
  double omega_der[3];
  double tol = 1e-12;
  if(LAMMPS_NS::vectorMag3D(wbody)*dt > 1.0)
    error->one(FLERR, "Timestep is too big for rotation integration!");
  dynamic_euler(wbody, tbody, inertia, omega_der);

  double omega_half_prev[3], delta[3];
  double omega_half[] = {0.0, 0.0, 0.0};
  while(1) {
    LAMMPS_NS::vectorCopy3D(omega_half, omega_half_prev);
    omega_half[0] = wbody[0] + 0.5*dt*omega_der[0];
    omega_half[1] = wbody[1] + 0.5*dt*omega_der[1];
    omega_half[2] = wbody[2] + 0.5*dt*omega_der[2];
    LAMMPS_NS::vectorSubtract3D(omega_half_prev, omega_half, delta);
    double omega_half_mag = LAMMPS_NS::vectorMag3D(omega_half);
    if(omega_half_mag > 0.0) {
      double eps = LAMMPS_NS::vectorMag3D(delta) / omega_half_mag;
      if(eps < tol)
        break;
      dynamic_euler(omega_half, tbody, inertia, omega_der);
    } else
      break;
  }
  wbody[0] += dt*omega_der[0];
  wbody[1] += dt*omega_der[1];
  wbody[2] += dt*omega_der[2];
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::integrate_quaternion(double dtq, double *quat, double *wbody)
{
  double q2[4];
  q2[0] = 1.0;
  q2[1] = -wbody[0] * dtq;
  q2[2] = -wbody[1] * dtq;
  q2[3] = -wbody[2] * dtq;
  MathExtraLiggghtsNonspherical::invquat(q2); //using implicit scheme

  double q_temp[4];
  MathExtra::quatquat(quat, q2, q_temp);
  vectorCopy4D(q_temp, quat);
  MathExtra::qnormalize(quat);
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::update_hdtorque(int i, double *rotation_matrix, double *omegaOld, double *omegaNew)
{
  if(ksl_rotation and hdtorque) {
    double deltaHydrotorquePrime[3], deltaHydrotorque[3];
    for(int dirI = 0; dirI < 3; dirI++)
      deltaHydrotorquePrime[dirI] = ksl_rotation[i][dirI]*(omegaOld[dirI]-omegaNew[dirI]);
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, deltaHydrotorquePrime, deltaHydrotorque);
    for(int dirI=0; dirI<3; dirI++)
      hdtorque[i][dirI] += deltaHydrotorque[dirI];
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::implicitRotationUpdate
(
    double deltaT, double* inertia,
    double *angMom, double *torque, double* KslRot,
    double *omegaNew,
    double *deltaHydrotorquePrime
)
{
      int index1_[3], index2_[3];
      index1_[0]=2;index1_[1]=0;index1_[2]=1;
      index2_[0]=1;index2_[1]=2;index2_[2]=0;
      double omegaOld[3], dtfm, KslMDeltaT;
      double omegaAngMomTerm;
      for(int dirI=0;dirI<3;dirI++)
      {
            if (inertia[dirI] == 0.0)
            {
              omegaOld[dirI] = 0.0;
              continue;
            }
            omegaOld[dirI] = angMom[dirI] / inertia[dirI];
      }
      for(int dirI=0;dirI<3;dirI++)
      {
            omegaAngMomTerm = omegaOld[index2_[dirI]]*angMom[index1_[dirI]]
                            - omegaOld[index1_[dirI]]*angMom[index2_[dirI]];

            dtfm           = deltaT
                           / ( inertia[dirI]);
            KslMDeltaT     = KslRot[dirI] * dtfm;

            //calculate new rotation velocity
            omegaNew[dirI] = (  omegaOld[dirI]
                              + dtfm
                              * (  torque[dirI] //this is the TOTAL torque!
                                 - omegaAngMomTerm
                                 + KslRot[dirI]*omegaOld[dirI]
                                )
                             )
                           /
                             ( 1.0+KslMDeltaT);

            //save increment in hdtorque to update torque since particle is now rotation with omegaNew
            deltaHydrotorquePrime[dirI] = KslRot[dirI]*(omegaOld[dirI]-omegaNew[dirI]);

            //update the angular momentum
            angMom[dirI] = omegaNew[dirI] * inertia[dirI];  //update velocity for a half step!
      }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::rotationUpdate(bool updateQuaternion)
{
  double omegaNew[3];
  double exone[3],eyone[3],ezone[3];
  double angMomPrime[3];
  double torquePrime[3];
  double omegaNewPrime[3];
  double deltaHydrotorquePrime[3];
  double deltaHydrotorque[3];
  double rot[3][3];

  // This is okey only for SPHERICAL particles
  // update 1/2 step for omega

  double **angmom  = atom->angmom;
  double **quat = atom->quaternion;
  double **inertia = atom->inertia;
  double **omega   = atom->omega;
  double **torque  = atom->torque;
  int    *mask     = atom->mask;
  int nlocal = atom->nlocal;

  double dtq = 0.50 * dtv;

  //save rotation rate to array if necessary
  orientation    = NULL;
  ksl_rotation    = NULL;
  hdtorque       = NULL;

  if(couple_fix_id > -1) {
      ksl_rotation = ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_KslRotation_->array_atom;
      hdtorque =     ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_hdtorque_->array_atom;
      orientation =  ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_ex_->array_atom;
    }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
    {
      MathExtra::quat_to_mat(quat[i],rot);     //rotation matrix

      //Compute all relevant quantities in the particle's coordinate system called 'prime'
      MathExtra::transpose_matvec(rot,angmom[i],angMomPrime);
      MathExtra::transpose_matvec(rot,torque[i],torquePrime);

/*      printf("\n***rotationUpdate(): angmom/prime: %g %g %g/%g %g %g, torque/prime: %g %g %g/%g %g %g \n\n",
             angmom[i][0],angmom[i][1],angmom[i][2],
             angMomPrime[0],angMomPrime[1],angMomPrime[2],
             torque[i][0],torque[i][1],torque[i][2],
             torquePrime[0],torquePrime[1],torquePrime[2]
             );

      printf("\n***rotationUpdate(): KslRotation: %g %g %g\n\n",
             KslRotation[i][0],KslRotation[i][1],KslRotation[i][2]
             );
*/
      //Implicit angmom and omegaNew update in the 'prime' coordinate system
      implicitRotationUpdate
      (
                dtf, inertia[i],
                angMomPrime, torquePrime, ksl_rotation[i],
                omegaNewPrime,
                deltaHydrotorquePrime
      );

      //Transform result to global coordinate system and update hydro torque (since it is the total!)
      MathExtra::matvec(rot,angMomPrime,angmom[i]);
      MathExtra::matvec(rot,omegaNewPrime,omegaNew);
      MathExtra::matvec(rot,deltaHydrotorquePrime,deltaHydrotorque); //save torque
      for(int dirI=0; dirI<3; dirI++)
          hdtorque[i][dirI] += deltaHydrotorque[dirI];

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion
      if(updateQuaternion)
        MathExtra::richardson(quat[i],angmom[i],omegaNew,inertia[i],dtq);

      omega[i][0]=omegaNew[0];
      omega[i][1]=omegaNew[1];
      omega[i][2]=omegaNew[2];

      if(ksl_rotation && updateQuaternion)
      {
          MathExtra::q_to_exyz(quat[i],exone,eyone,ezone);   
          orientation[i][0] = exone[0];
          orientation[i][1] = exone[1];
          orientation[i][2] = exone[2];
       }
    } //end of particle
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::initial_integrate(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double **angmom = atom->angmom;
  double **quat = atom->quaternion;
  double **inertia = atom->inertia;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  orientation    = NULL;
  ksl_rotation    = NULL;
  hdtorque       = NULL;

  if(couple_fix_id > -1) {
    ksl_rotation = ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_KslRotation_->array_atom;
    hdtorque =     ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_hdtorque_->array_atom;
    orientation =  ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_ex_->array_atom;
  }

  double tbody[3], rotation_matrix[9];
  double fquat[4], mbody[3], wbody[3], conjqm[4];

  double angMomPrime[3];
  double omegaNewPrime[3], omegaOldPrime[3];
  const double dtf2 = 2.0 * dtf;

  const double dtq = 0.5 * dtv;
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag3D(x[i])))
        error->fix_error(FLERR,this,"x[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(v[i])))
        error->fix_error(FLERR,this,"v[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(f[i])))
        error->fix_error(FLERR,this,"f[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag4D(quat[i])))
        error->fix_error(FLERR,this,"quat[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(angmom[i])))
        error->fix_error(FLERR,this,"angmom[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(omega[i])))
        error->fix_error(FLERR,this,"omega[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(torque[i])))
        error->fix_error(FLERR,this,"torque[i] is NaN!");
#endif
      MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);
      const double dtfm = dtf / rmass[i];
      MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom[i], angMomPrime);
      omegaOldPrime[0] = angMomPrime[0] / inertia[i][0];
      omegaOldPrime[1] = angMomPrime[1] / inertia[i][1];
      omegaOldPrime[2] = angMomPrime[2] / inertia[i][2];

      //update velocity by step t+1/2dt
      
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag3D(v[i])))
        error->fix_error(FLERR,this,"v[i] is NaN!");
#endif
      //update position by step t+dt
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag3D(x[i])))
        error->fix_error(FLERR,this,"x[i] is NaN!");
#endif

      if(integration_scheme == 0)
      {
        //update angular moment by step t+1/2dt
        angmom[i][0] += dtf * torque[i][0];
        angmom[i][1] += dtf * torque[i][1];
        angmom[i][2] += dtf * torque[i][2];
        //update angular velocity by step t+1/2dt
        MathExtra::mq_to_omega(angmom[i],quat[i],inertia[i],omega[i]);
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega[i], omegaNewPrime);
        update_hdtorque(i, rotation_matrix, omegaOldPrime, omegaNewPrime);
        //update quaternion by step t+dt
        MathExtra::richardson(quat[i],angmom[i],omega[i],inertia[i],dtq);
      }
      //symplectic scheme
      else if(integration_scheme == 1)
      {
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom[i],mbody);
        MathExtraLiggghtsNonspherical::calc_conjqm(quat[i],mbody,conjqm);

        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque[i], tbody);
        MathExtra::quatvec(quat[i], tbody, fquat);

        conjqm[0] += dtf2 * fquat[0];
        conjqm[1] += dtf2 * fquat[1];
        conjqm[2] += dtf2 * fquat[2];
        conjqm[3] += dtf2 * fquat[3];

        MathExtraLiggghtsNonspherical::no_squish_rotate(3,conjqm,quat[i],inertia[i],dtq);
        MathExtraLiggghtsNonspherical::no_squish_rotate(2,conjqm,quat[i],inertia[i],dtq);
        MathExtraLiggghtsNonspherical::no_squish_rotate(1,conjqm,quat[i],inertia[i],dtv);
        MathExtraLiggghtsNonspherical::no_squish_rotate(2,conjqm,quat[i],inertia[i],dtq);
        MathExtraLiggghtsNonspherical::no_squish_rotate(3,conjqm,quat[i],inertia[i],dtq);

        MathExtra::qnormalize(quat[i]);
        MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);

        MathExtra::invquatvec(quat[i],conjqm,mbody);
        mbody[0] *= 0.5;
        mbody[1] *= 0.5;
        mbody[2] *= 0.5;

        wbody[0] = mbody[0] / inertia[i][0];
        wbody[1] = mbody[1] / inertia[i][1];
        wbody[2] = mbody[2] / inertia[i][2];

        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom[i]);
        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega[i]);
      } else if (integration_scheme == 2) { //direct integration, 2nd order predictor-corrector
        double omega_half[3];
        MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix); //calculate rotation matrix from quaternion
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega[i],wbody); //angular velocity in body principal axes
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque[i], tbody); //torque in body principal axes
        omega_half[0] = wbody[0];
        omega_half[1] = wbody[1];
        omega_half[2] = wbody[2];
        integrate_dynamic_euler(dtf, wbody, tbody, inertia[i]); //calculate angular velocity at step t+dt
        omega_half[0] = 0.5*(wbody[0] + omega_half[0]); //angular velocity at step t+dt/2
        omega_half[1] = 0.5*(wbody[1] + omega_half[1]);
        omega_half[2] = 0.5*(wbody[2] + omega_half[2]);

        integrate_quaternion(dtq, quat[i], wbody); //calculate quaternion at step t+dt
        update_hdtorque(i, rotation_matrix, omegaOldPrime, wbody);

        mbody[0] = inertia[i][0]*wbody[0]; //calculate angular momentum at step t+dt from angular velocity in body's principal axes
        mbody[1] = inertia[i][1]*wbody[1];
        mbody[2] = inertia[i][2]*wbody[2];

        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom[i]); // angular momentum to global frame
        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega[i]); // angular velocity to global frame

      }  else if(integration_scheme == 3) { //woodem scheme
        double angmom_half[3], angmom_half_local[3], angmom_next_local[3];

        for(int j = 0; j < 3; j++)
          angmom_half[j] = angmom[i][j] + torque[i][j]*dtf;
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom_half, angmom_half_local);
        for(int j = 0; j < 3; j++)
          angmom[i][j] += torque[i][j]*dtf2;
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom[i], angmom_next_local);

        double omega_half_local[3], omega_next_local[3], q_der[4], quat_half[4];
        for(int j = 0; j < 3; j++) {
          omega_half_local[j] = angmom_half_local[j] / inertia[i][j];
          omega_next_local[j] = angmom_next_local[j] / inertia[i][j];
        }
        update_hdtorque(i, rotation_matrix, omegaOldPrime, omega_next_local);

        MathExtra::quatvec(quat[i], omega_half_local, q_der);
        for(int j = 0; j < 4; j++)
          quat_half[j] = quat[i][j] + q_der[j]*dtq*0.5;
        MathExtra::qnormalize(quat_half);

        MathExtra::quatvec(quat_half, omega_next_local, q_der);
        for(int j = 0; j < 4; j++)
          quat[i][j] += q_der[j]*dtq;
        MathExtra::qnormalize(quat[i]);

        MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);
        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, omega_next_local, omega[i]);
      } else if(integration_scheme == 4) {
        rotationUpdate(true);
      }
      else
        error->one(FLERR,"Invalid integration scheme! Please choose between 0, 1 (default), 2 or 3!");
#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag4D(quat[i])))
        error->fix_error(FLERR,this,"quat[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(angmom[i])))
        error->fix_error(FLERR,this,"angmom[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(omega[i])))
        error->fix_error(FLERR,this,"omega[i] is NaN!");
#endif
      if(orientation) {
        double exone[3], eyone[3], ezone[3];
        MathExtra::q_to_exyz(quat[i],exone,eyone,ezone);
        orientation[i][0] = exone[0];
        orientation[i][1] = exone[1];
        orientation[i][2] = exone[2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereBase::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double **angmom = atom->angmom;
  double **quat = atom->quaternion;
  double **inertia = atom->inertia;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double tbody[3], rotation_matrix[9];
  double fquat[4], mbody[3], wbody[3];
  double dtf2 = dtf * 2.0;
  double conjqm[4];
  double angMomPrime[3], omegaNewPrime[3], omegaOldPrime[3];

  ksl_rotation = 0;
  hdtorque = 0;
  orientation = 0;
  if(couple_fix_id > -1) {
    ksl_rotation = ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_KslRotation_->array_atom;
    hdtorque =     ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_hdtorque_->array_atom;
    orientation =  ((FixCfdCouplingForceImplicit*)modify->fix[couple_fix_id])->fix_ex_->array_atom;
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag3D(v[i])))
        error->fix_error(FLERR,this,"v[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(f[i])))
        error->fix_error(FLERR,this,"f[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag4D(quat[i])))
        error->fix_error(FLERR,this,"quat[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(angmom[i])))
        error->fix_error(FLERR,this,"angmom[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(omega[i])))
        error->fix_error(FLERR,this,"omega[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(torque[i])))
        error->fix_error(FLERR,this,"torque[i] is NaN!");
#endif
      const double dtfm = dtf / rmass[i];
      MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);
      MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom[i], angMomPrime);
      omegaOldPrime[0] = angMomPrime[0] / inertia[i][0];
      omegaOldPrime[1] = angMomPrime[1] / inertia[i][1];
      omegaOldPrime[2] = angMomPrime[2] / inertia[i][2];
      if(hdtorque) {
        vectorAdd3D(torque[i], hdtorque[i], torque[i]);
        vectorAdd3D(torque[i], hdtorque[i], torque[i]);
      }

      //update velocity by step t+dt
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag3D(v[i])))
        error->fix_error(FLERR,this,"v[i] is NaN!");
#endif

      if(integration_scheme == 0) {
        //update angular moment by step t+dt
        angmom[i][0] += dtf * torque[i][0];
        angmom[i][1] += dtf * torque[i][1];
        angmom[i][2] += dtf * torque[i][2];
        //update angular velocity by step t+dt
        MathExtra::mq_to_omega(angmom[i],quat[i],inertia[i],omega[i]);
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega[i], omegaNewPrime);
        update_hdtorque(i, rotation_matrix, omegaOldPrime, omegaNewPrime);
      } else if(integration_scheme == 1) {
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque[i], tbody);
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom[i], mbody);
        MathExtraLiggghtsNonspherical::calc_conjqm(quat[i], mbody, conjqm);

        MathExtra::quatvec(quat[i], tbody, fquat);

        conjqm[0] += dtf2 * fquat[0];
        conjqm[1] += dtf2 * fquat[1];
        conjqm[2] += dtf2 * fquat[2];
        conjqm[3] += dtf2 * fquat[3];

        MathExtra::invquatvec(quat[i], conjqm, mbody);
        mbody[0] *= 0.5;
        mbody[1] *= 0.5;
        mbody[2] *= 0.5;

        wbody[0] = mbody[0] / inertia[i][0];
        wbody[1] = mbody[1] / inertia[i][1];
        wbody[2] = mbody[2] / inertia[i][2];

        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom[i]);
        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega[i]);
        update_hdtorque(i, rotation_matrix, omegaOldPrime, wbody);
      } else if (integration_scheme == 2) {
        MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix); //calculate rotation matrix from quaternion
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega[i],wbody); //angular velocity in body principal axes
        MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque[i], tbody); //torque in body principal axes

        integrate_dynamic_euler(dtf, wbody, tbody, inertia[i]); //calculate angular velocity at step t+dt

        update_hdtorque(i, rotation_matrix, omegaOldPrime, wbody);

        mbody[0] = inertia[i][0]*wbody[0]; //calculate angular momentum at step t+dt from angular velocity in body's principal axes
        mbody[1] = inertia[i][1]*wbody[1];
        mbody[2] = inertia[i][2]*wbody[2];

        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom[i]); // angular momentum to global frame
        MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega[i]); // angular velocity to global frame
      } else if(integration_scheme == 4) {
        rotationUpdate(false);
      }

#ifdef LIGGGHTS_DEBUG
      if(std::isnan(LAMMPS_NS::vectorMag4D(quat[i])))
        error->fix_error(FLERR,this,"quat[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(angmom[i])))
        error->fix_error(FLERR,this,"angmom[i] is NaN!");
      if(std::isnan(LAMMPS_NS::vectorMag3D(omega[i])))
        error->fix_error(FLERR,this,"omega[i] is NaN!");
#endif
    }
}

