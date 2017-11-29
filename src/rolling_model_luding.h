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

   Contributing authors:
   Rahul Mohanty (University of Edinburgh, P&G)
------------------------------------------------------------------------- */
#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_LUDING,luding,6)
#else
#ifndef ROLLING_MODEL_LUDING_H_
#define ROLLING_MODEL_LUDING_H_
#include "contact_models.h"
#include "rolling_model_base.h"
#include <algorithm>
#include <cmath>
#include "math_vector.h"
#include "domain.h"
#include "math_extra_liggghts.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<>
  class RollingModel<ROLLING_LUDING> : public RollingModelBase
  {

  public:
    RollingModel(class LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) : RollingModelBase(lmp, hsetup, cmb), coeffRollFrict(NULL), coeffRollVisc(NULL), kR2kcMax(NULL),torsion(false)
    {
      history_offset = hsetup->add_history_value("r_torquex_old", "1");
      hsetup->add_history_value("r_torquey_old", "1");
      hsetup->add_history_value("r_torquez_old", "1");
      hsetup->add_history_value("r_tor_torquex_old", "1");
      hsetup->add_history_value("r_tor_torquey_old", "1");
      hsetup->add_history_value("r_tor_torquey_old", "1");
      kc_offset = cmb->get_history_offset("kc_offset");
      fo_offset = cmb->get_history_offset("fo_offset");
    }

    void registerSettings(Settings& settings) {
      settings.registerOnOff("torsion",torsion,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.registerProperty("coeffRollVisc", &MODEL_PARAMS::createCoeffRollVisc);
      registry.registerProperty("kR2kcMax", &MODEL_PARAMS::createCoeffRollingStiffness);

      registry.connect("coeffRollFrict", coeffRollFrict,"rolling_model luding");
      registry.connect("coeffRollVisc", coeffRollVisc,"rolling_model luding");
      registry.connect("kR2kcMax", kR2kcMax,"rolling_model luding");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"rolling model luding");
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      double r_torque[3],r_tor_torque[3];
      vectorZeroize3D(r_torque);
      vectorZeroize3D(r_tor_torque);

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_ROLLING_MODEL;

      double radi = sidata.radi;
      double radj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : (radi*radj/(radi+radj));

#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif

      if(sidata.is_wall) {
        const double wr1 = sidata.wr1;
        const double wr2 = sidata.wr2;
        const double wr3 = sidata.wr3;

        double r_inertia = 0.0; //pre-initialize to prevent compiler "warning"
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          const double rii = pointDistance(sidata.contact_point, atom->x[sidata.i]);
          const double omega_mag = sqrt(wr1*wr1 + wr2*wr2 + wr3*wr3);
          if(omega_mag != 0.0) {
            double er[3];
            er[0] = wr1 / omega_mag;
            er[1] = wr2 / omega_mag;
            er[2] = wr3 / omega_mag;
            const double Ix = atom->inertia[sidata.i][0];
            const double Iy = atom->inertia[sidata.i][1];
            const double Iz = atom->inertia[sidata.i][2];
            double inertia_tensor[9];
            double inertia_tensor_local[9] = { Ix, 0.0, 0.0,
                                               0.0, Iy, 0.0,
                                               0.0, 0.0, Iz };
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local, atom->quaternion[sidata.i], inertia_tensor);
            double temp[3];
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor, er, temp);
            double Ii = MathExtra::dot3(temp, er);
            r_inertia = Ii + sidata.mi*rii*rii;
          }
        } else {
          if (domain->dimension == 2) r_inertia = 1.5*sidata.mi*reff*reff;
          else  r_inertia = 1.4*sidata.mi*reff*reff;
        }
#else
        if (domain->dimension == 2) r_inertia = 1.5*sidata.mi*reff*reff;
        else  r_inertia = 1.4*sidata.mi*reff*reff;
#endif

        calcRollTorque(r_torque,sidata,reff,wr1,wr2,wr3,r_inertia);

        if(torsion){
          const int i = sidata.i;
          const double * const * const omega = atom->omega;
          double const* const normal = sidata.en;
          double crossOmegaI;
          crossOmegaI = vectorDot3D(normal,omega[i]);
          double omega_normal[3];
          vectorZeroize3D(omega_normal);
          omega_normal[0] = normal[0] * crossOmegaI;
          omega_normal[1] = normal[1] * crossOmegaI;
          omega_normal[2] = normal[2] * crossOmegaI;
          calcTorTorque(r_tor_torque,sidata,reff,omega_normal[0],omega_normal[1],omega_normal[2],r_inertia);
        }
      } else {
        double  wr_roll[3];
        const int i = sidata.i;
        const int j = sidata.j;
        const double * const * const omega = atom->omega;

        // relative rotational velocity
        vectorSubtract3D(omega[i],omega[j],wr_roll);
        double r_inertia = 0.0; //pre-initialize to prevent compiler "warning"
        double r_inertia_red_i, r_inertia_red_j;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          const double rii = pointDistance(sidata.contact_point, atom->x[i]);
          const double rjj = pointDistance(sidata.contact_point, atom->x[j]);
          const double omega_mag = vectorMag3D(wr_roll);
          if(omega_mag != 0.0) {
            double er[3];
            er[0] = wr_roll[0] / omega_mag;
            er[1] = wr_roll[1] / omega_mag;
            er[2] = wr_roll[2] / omega_mag;
            const double Ix_i = atom->inertia[i][0];
            const double Iy_i = atom->inertia[i][1];
            const double Iz_i = atom->inertia[i][2];

            const double Ix_j = atom->inertia[j][0];
            const double Iy_j = atom->inertia[j][1];
            const double Iz_j = atom->inertia[j][2];

            double inertia_tensor_i[9];
            double inertia_tensor_local_i[9] = { Ix_i, 0.0, 0.0,
                                                 0.0, Iy_i, 0.0,
                                                 0.0, 0.0, Iz_i };
            double inertia_tensor_j[9];
            double inertia_tensor_local_j[9] = { Ix_j, 0.0, 0.0,
                                                 0.0, Iy_j, 0.0,
                                                 0.0, 0.0, Iz_j };
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local_i, atom->quaternion[i], inertia_tensor_i);
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local_j, atom->quaternion[j], inertia_tensor_j);
            double temp[3];
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor_i, er, temp);
            double Ii = MathExtra::dot3(temp, er);
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor_j, er, temp);
            double Ij = MathExtra::dot3(temp, er);
            r_inertia_red_i = Ii + sidata.mi*rii*rii;
            r_inertia_red_j = Ij + sidata.mj*rjj*rjj;
            r_inertia = r_inertia_red_i*r_inertia_red_j / (r_inertia_red_i + r_inertia_red_j);
          }

        } else {
          r_inertia_red_i = sidata.mi*radi*radi;
          r_inertia_red_j = sidata.mj*radj*radj;
          if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
          else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
        }
#else
        r_inertia_red_i = sidata.mi*radi*radi;
        r_inertia_red_j= sidata.mj*radj*radj;
        if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
        else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
#endif

        calcRollTorque(r_torque,sidata,reff,wr_roll[0],wr_roll[1],wr_roll[2],r_inertia);

        if(torsion){
          const int i = sidata.i;
          const int j = sidata.j;
          const double * const * const omega = atom->omega;
          double const* const normal = sidata.en;
          double crossOmegaI;
          double crossOmegaJ;
          crossOmegaI = vectorDot3D(normal,omega[i]);
          crossOmegaJ = vectorDot3D(normal,omega[j]);
          double omega_normal[3];
          vectorZeroize3D(omega_normal);
          omega_normal[0] = normal[0] * (crossOmegaI - crossOmegaJ);
          omega_normal[1] = normal[1] * (crossOmegaI - crossOmegaJ);
          omega_normal[2] = normal[2] * (crossOmegaI - crossOmegaJ);
          calcTorTorque(r_torque,sidata,reff,omega_normal[0],omega_normal[1],omega_normal[2],r_inertia);
        }
      }

      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];
      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];

      if(torsion){
        i_forces.delta_torque[0] -= r_tor_torque[0];
        i_forces.delta_torque[1] -= r_tor_torque[1];
        i_forces.delta_torque[2] -= r_tor_torque[2];
        j_forces.delta_torque[0] += r_tor_torque[0];
        j_forces.delta_torque[1] += r_tor_torque[1];
        j_forces.delta_torque[2] += r_tor_torque[2];
      }
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_ROLLING_MODEL;
      double * const c_history = &scdata.contact_history[history_offset];
      c_history[0] = 0.0; // this is the r_torque_old
      c_history[1] = 0.0; // this is the r_torque_old
      c_history[2] = 0.0; // this is the r_torque_old
      double * const c_tor_history = &scdata.contact_history[history_offset];
      c_tor_history[0] = 0.0; // this is the r_tor_torque_old
      c_tor_history[1] = 0.0; // this is the r_tor_torque_old
      c_tor_history[2] = 0.0; // this is the r_tor_torque_old
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;
    double ** coeffRollVisc;
    double ** kR2kcMax;
    int history_offset;
    int kc_offset;
    int fo_offset;
    bool torsion;

    inline void calcRollTorque(double (&r_torque)[3],const SurfacesIntersectData & sidata,double reff,double wr1,double wr2,double wr3,double r_inertia) {
      double wr_n[3],wr_t[3],dr_torque[3];
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;

      const double enx = sidata.en[0];
      const double eny = sidata.en[1];
      const double enz = sidata.en[2];

      const double dt = update->dt;

      double * const c_history = &sidata.contact_history[history_offset]; // requires Style::TANGENTIAL == TANGENTIAL_HISTORY
      const double rmu= coeffRollFrict[itype][jtype];

      // remove normal (torsion) part of relative rotation
      // use only tangential parts for rolling torque
      const double wr_dot_delta = wr1*enx+ wr2*eny + wr3*enz;
      wr_n[0] = enx * wr_dot_delta;
      wr_n[1] = eny * wr_dot_delta;
      wr_n[2] = enz * wr_dot_delta;
      wr_t[0] = wr1 - wr_n[0];
      wr_t[1] = wr2 - wr_n[1];
      wr_t[2] = wr3 - wr_n[2];

      // spring
      const double k_r = sidata.kt;
      const double kr = k_r * kR2kcMax[itype][jtype]*reff*reff;  // tangential stiffness based on the ratio input

      const double kc = sidata.contact_history[kc_offset];
      const double f_adh = sidata.contact_history[fo_offset];

      vectorScalarMult3D(wr_t,dt*kr,dr_torque);

      r_torque[0] = c_history[0] + dr_torque[0];
      r_torque[1] = c_history[1] + dr_torque[1];
      r_torque[2] = c_history[2] + dr_torque[2];

      // limit max. torque
      const double r_torque_mag = vectorMag3D(r_torque);
      const double r_torque_max = rmu * fabs(sidata.Fn + kc * sidata.deltan-f_adh)*reff;   // Modified Columb limit
      if(r_torque_mag > r_torque_max)
      {
        const double factor = r_torque_max / r_torque_mag;
        r_torque[0] *= factor;
        r_torque[1] *= factor;
        r_torque[2] *= factor;
        // save rolling torque due to spring
        c_history[0] = r_torque[0];
        c_history[1] = r_torque[1];
        c_history[2] = r_torque[2];
        // no damping / no dashpot in case of full mobilisation rolling angle
      } else {
        // save rolling torque due to spring before adding damping torque
        c_history[0] = r_torque[0];
        c_history[1] = r_torque[1];
        c_history[2] = r_torque[2];
        // dashpot
        const double r_coef = sidata.gammat * coeffRollVisc[itype][jtype] * 2 * sqrt(r_inertia*kr);
        // add damping torque
        r_torque[0] += r_coef*wr_t[0];
        r_torque[1] += r_coef*wr_t[1];
        r_torque[2] += r_coef*wr_t[2];
      }
    }
    inline void calcTorTorque(double (&r_tor_torque)[3],const SurfacesIntersectData & sidata,double reff,double wr1,double wr2,double wr3,double r_inertia) {
      double wr_n[3],dr_tor_torque[3];

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;

      const double dt = update->dt;
      double * const c_tor_history = &sidata.contact_history[history_offset];
      const double rmu= coeffRollFrict[itype][jtype];   // rmu is used as torsion coefficient

      const double kc = sidata.contact_history[kc_offset];
      const double f_adh = sidata.contact_history[fo_offset];
      // use of normal component of omega for torque calculations
      wr_n[0] = wr1;
      wr_n[1] = wr2;
      wr_n[2] = wr3;

      // torsion spring
      const double k_r = sidata.kt;
      const double kr = k_r * kR2kcMax[itype][jtype]*reff*reff;  // rolling stiffness is used as torsion stiffness

      vectorScalarMult3D(wr_n,dt*kr,dr_tor_torque);

      r_tor_torque[0] = c_tor_history[0] + dr_tor_torque[0];
      r_tor_torque[1] = c_tor_history[1] + dr_tor_torque[1];
      r_tor_torque[2] = c_tor_history[2] + dr_tor_torque[2];

      // limit max. torsion torque
      const double r_torque_mag = vectorMag3D(r_tor_torque);
      const double r_torque_max = rmu * fabs(sidata.Fn + kc * sidata.deltan-f_adh)*reff;

      if(r_torque_mag > r_torque_max)
      {
        const double factor = r_torque_max / r_torque_mag;

        r_tor_torque[0] *= factor;
        r_tor_torque[1] *= factor;
        r_tor_torque[2] *= factor;
        // save torsion torque due to spring
        c_tor_history[0] = r_tor_torque[0];
        c_tor_history[1] = r_tor_torque[1];
        c_tor_history[2] = r_tor_torque[2];
        // no damping / no dashpot in case of full mobilisation torsion angle
      } else {
        // save torsion torque due to spring before adding damping torque
        c_tor_history[0] = r_tor_torque[0];
        c_tor_history[1] = r_tor_torque[1];
        c_tor_history[2] = r_tor_torque[2];
        // dashpot
        const double r_coef = sidata.gammat * coeffRollVisc[itype][jtype] * 2.0 * sqrt(r_inertia*kr);
        // add damping torque
        r_tor_torque[0] += r_coef*wr_n[0];
        r_tor_torque[1] += r_coef*wr_n[1];
        r_tor_torque[2] += r_coef*wr_n[2];
      }
    }
  };
}
}
#endif // ROLLING_MODEL_LUDING_H_
#endif
