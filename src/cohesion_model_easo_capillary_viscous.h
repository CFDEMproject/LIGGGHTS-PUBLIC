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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_EASO_CAPILLARY_VISCOUS,easo/capillary/viscous,8)
#else

#ifndef COHESION_MODEL_EASO_CAPILLARY_VISCOUS_H_
#define COHESION_MODEL_EASO_CAPILLARY_VISCOUS_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include <algorithm>
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"

namespace MODEL_PARAMS
{
    inline static ScalarProperty* createliquidContentInitialEaso(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* surfaceLiquidContentInitialScalar = MODEL_PARAMS::createScalarProperty(registry, "surfaceLiquidContentInitial", caller);
      return surfaceLiquidContentInitialScalar;
    }

    inline static ScalarProperty* createMinSeparationDistanceRatioEaso(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* minSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDistanceRatio", caller);
      return minSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createMaxSeparationDistanceRatioEaso(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* maxSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistanceRatio", caller);
      return maxSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createFluidViscosityEaso(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
      return fluidViscosityScalar;
    }
}

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_EASO_CAPILLARY_VISCOUS> : public CohesionModelBase {

  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      surfaceLiquidContentInitial(0.0),
      surfaceTension(0.0),
      contactAngle(0),
      minSeparationDistanceRatio(0.0),
      maxSeparationDistanceRatio(0.0),
      fluidViscosity(0.),
      history_offset(0),
      fix_surfaceliquidcontent(0),
      fix_liquidflux(0),
      fix_ste(0)
    {
      history_offset = hsetup->add_history_value("contflag", "0");
      
      if(cmb->is_wall())
        error->warning(FLERR,"Using cohesion model easo/capillary/viscous for walls only supports dry walls");
    }

    void registerSettings(Settings& settings)
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("surfaceLiquidContentInitial", &MODEL_PARAMS::createliquidContentInitialEaso);
      registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
      registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityEaso);
      registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngle);
      registry.registerProperty("minSeparationDistanceRatio", &MODEL_PARAMS::createMinSeparationDistanceRatioEaso);
      registry.registerProperty("maxSeparationDistanceRatio", &MODEL_PARAMS::createMaxSeparationDistanceRatioEaso);

      registry.connect("surfaceLiquidContentInitial", surfaceLiquidContentInitial,"cohesion_model easo/capillary/viscous");
      registry.connect("surfaceTension", surfaceTension,"cohesion_model easo/capillary/viscous");
      registry.connect("fluidViscosity", fluidViscosity,"cohesion_model easo/capillary/viscous");
      registry.connect("contactAngle", contactAngle,"cohesion_model easo/capillary/viscous");
      registry.connect("minSeparationDistanceRatio", minSeparationDistanceRatio,"cohesion_model easo/capillary/viscous");
      
      registry.connect("maxSeparationDistanceRatio", maxSeparationDistanceRatio,"cohesion_model easo/capillary/viscous");

      ln1overMinSeparationDistanceRatio = log(1./minSeparationDistanceRatio);

      fix_ste = modify->find_fix_scalar_transport_equation("liquidtransfer");
      if(!fix_ste)
      {
        
        char initstr[200];
        sprintf(initstr,"%e",surfaceLiquidContentInitial);

        const char * newarg[15];
        newarg[0] = "liquidtransfer";
        newarg[1] = "all";
        newarg[2] = "transportequation/scalar";
        newarg[3] = "equation_id";
        newarg[4] = "liquidtransfer";
        newarg[5] = "quantity";
        newarg[6] = "surfaceLiquidContent";
        newarg[7] = "default_value";
        newarg[8] = initstr;
        newarg[9] = "flux_quantity";
        newarg[10] = "liquidFlux";
        newarg[11] = "source_quantity";
        newarg[12] = "liquidSource";
        newarg[13] = "capacity_quantity";
        newarg[14] = "none";
        modify->add_fix(15,const_cast<char**>(newarg));
      }

      fix_surfaceliquidcontent = static_cast<FixPropertyAtom*>(modify->find_fix_property("surfaceLiquidContent","property/atom","scalar",0,0,"cohesion_model easo/capillary/viscous"));
      fix_liquidflux = static_cast<FixPropertyAtom*>(modify->find_fix_property("liquidFlux","property/atom","scalar",0,0,"cohesion_model easo/capillary/viscous"));
      fix_ste = modify->find_fix_scalar_transport_equation("liquidtransfer");

      if(!fix_surfaceliquidcontent || !fix_liquidflux || !fix_ste)
          error->all(FLERR,"internal error");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model easo/capillary/viscous");

      neighbor->register_contact_dist_factor(maxSeparationDistanceRatio*1.1); 
      if(maxSeparationDistanceRatio < 1.0)
            error->one(FLERR,"\n\ncohesion model easo/capillary/viscous requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int i = sidata.i;
      const int j = sidata.j;
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;

      const double radi = sidata.radi;
      const double radj = sidata.is_wall ? radi : sidata.radj;
      const double radsum = sidata.radsum;
      double const *surfaceLiquidContent = fix_surfaceliquidcontent->vector_atom;

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
      double * const contflag = &sidata.contact_history[history_offset];
      // store for noCollision
      contflag[0] = 1.0;

      const double volLi1000 = /* 4/3 * 1000 */ 1333.333333*M_PI*radi*radi*radi*surfaceLiquidContent[i];
      const double volLj1000 = /* 4/3 * 1000 */ sidata.is_wall ? 0.0 : 1333.333333*M_PI*radj*radj*radj*surfaceLiquidContent[j];
      const double volLiBond1000 = 0.5*volLi1000*(1.-sqrt(1.-radj*radj/(radsum*radsum)));
      const double volLjBond1000 = 0.5*volLj1000*(1.-sqrt(1.-radi*radi/(radsum*radsum)));
      const double volBond1000 = volLiBond1000+volLjBond1000;

      // skip if bond volume too small
      if(volBond1000 < 1e-14) return;

      const double rEff = radi*radj / (radi+radj);
      const double contactAngleEff = 0.5 * contactAngle[itype] * contactAngle[jtype];

      // capilar force
      // this is from Soulie et al, Intl. J Numerical and Analytical Methods in Geomechanics
      // 30 (2006), 213-228, Eqn. 13,14; separation distance = 0 in this case
      
      const double R2 = (radi >=radj) ? radi : radj;
      const double R2inv = 1./R2;
      const double volBondScaled = volBond1000*R2inv*0.001*R2inv*R2inv;
      const double Bparam = (-0.148*log(volBondScaled)-0.96)*contactAngleEff*contactAngleEff - 0.0082*log(volBondScaled) + 0.48;
      const double Cparam = 0.0018*log(volBondScaled)+0.078;
      const double Fcapilary = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Bparam)+Cparam);

      // viscous force
      // this is from Nase et al as cited in Shi and McCarthy, Powder Technology, 184 (2008), 65-75, Eqns 40,41
      const double stokesPreFactor = -6.*M_PI*fluidViscosity*rEff;
      const double FviscN = stokesPreFactor*sidata.vn/minSeparationDistanceRatio;
      const double FviscT_over_vt = (/* 8/15 */ 0.5333333*ln1overMinSeparationDistanceRatio + 0.9588) * stokesPreFactor;

      // tangential force components
      const double Ft1 = FviscT_over_vt*sidata.vtr1;
      const double Ft2 = FviscT_over_vt*sidata.vtr2;
      const double Ft3 = FviscT_over_vt*sidata.vtr3;

      // torques
      const double tor1 = sidata.en[1] * Ft3 - sidata.en[2] * Ft2;
      const double tor2 = sidata.en[2] * Ft1 - sidata.en[0] * Ft3;
      const double tor3 = sidata.en[0] * Ft2 - sidata.en[1] * Ft1;

      // add to fn, Ft
      if(tangentialReduce_) sidata.Fn += Fcapilary+FviscN;  
      //sidata.Ft += ...

      // apply normal and tangential force
      const double fx = (Fcapilary+FviscN) * sidata.en[0] + Ft1;
      const double fy = (Fcapilary+FviscN) * sidata.en[1] + Ft2;
      const double fz = (Fcapilary+FviscN) * sidata.en[2] + Ft3;

      // return resulting forces
      if(sidata.is_wall) {
        const double area_ratio = sidata.area_ratio;
        i_forces.delta_F[0] += fx * area_ratio;
        i_forces.delta_F[1] += fy * area_ratio;
        i_forces.delta_F[2] += fz * area_ratio;
        i_forces.delta_torque[0] += -sidata.cri * tor1 * area_ratio;
        i_forces.delta_torque[1] += -sidata.cri * tor2 * area_ratio;
        i_forces.delta_torque[2] += -sidata.cri * tor3 * area_ratio;
      } else {
        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;
        i_forces.delta_torque[0] += -sidata.cri * tor1;
        i_forces.delta_torque[1] += -sidata.cri * tor2;
        i_forces.delta_torque[2] += -sidata.cri * tor3;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
        j_forces.delta_torque[0] += -sidata.crj * tor1;
        j_forces.delta_torque[1] += -sidata.crj * tor2;
        j_forces.delta_torque[2] += -sidata.crj * tor3;
      }
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {

      double * const contflag = &scdata.contact_history[history_offset];

      // 3 cases: (i) no bridge present, (ii) bridge active, (iii) bridge breaks this step
      // for this model, bridge is created at contact and breaks at rupture distancy

      // case (i) no bridge
      if(!MathExtraLiggghts::compDouble(contflag[0],1.0,1e-6))
      {
          if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
          return;
      }

      // cases (i) and (ii)

      const int i = scdata.i;
      const int j = scdata.j;
      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double radi = scdata.radi;
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);
      const double radsum = scdata.radsum;
      const double dist = scdata.is_wall ? r - radi : r - (radi + radj);
      double const *surfaceLiquidContent = fix_surfaceliquidcontent->vector_atom;

      const double volLi1000 = /* 4/3 * 1000 */ 1333.333333*M_PI*radi*radi*radi*surfaceLiquidContent[i];
      const double volLj1000 = /* 4/3 * 1000 */ scdata.is_wall ? 0.0 : 1333.333333*M_PI*radj*radj*radj*surfaceLiquidContent[j];
      const double volLiBond1000 = 0.5*volLi1000*(1.-sqrt(1.-radj*radj/(radsum*radsum)));
      const double volLjBond1000 = 0.5*volLj1000*(1.-sqrt(1.-radi*radi/(radsum*radsum)));
      const double volBond1000 = volLiBond1000+volLjBond1000;

      const double rEff = radi*radj / (radi+radj);
      const double contactAngleEff = 0.5 * contactAngle[itype] * contactAngle[jtype];
      const double distMax = (1. + 0.5*contactAngleEff) * cbrt(volBond1000) *0.1;

      // check if liquid bridge exists
      bool bridge_active = false, bridge_breaks = false;

      if (dist > (maxSeparationDistanceRatio-1.0)*(radi+radj) && MathExtraLiggghts::compDouble(contflag[0],1.0,1e-6)) // in this case always break
      {
        bridge_breaks = true;
      }
      else if (dist < distMax && dist < (maxSeparationDistanceRatio-1.0)*(radi+radj) )
        bridge_active = true;
      else if(MathExtraLiggghts::compDouble(contflag[0],1.0,1e-6)) // can only break if exists
        bridge_breaks = true;

      // case (ii)
      if(bridge_active)
      {
          if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
          double **v = atom->v;

          // store for next step
          contflag[0] = 1.0;

          // skip if bond volume too small
          if(volBond1000 < 1e-14) return;

          // calculate forces, case no collision

          // capilary force
          // this is from Soulie et al, Intl. J Numerical and Analytical Methods in Geomechanics
          // 30 (2006), 213-228, Eqn. 13,14; separation distance = 0 in this case
          
          const double R2 = (radi >=radj) ? radi : radj;
          const double R2inv = 1./R2;
          const double volBondScaled = volBond1000*R2inv*0.001*R2inv*R2inv;
          const double Aparam = -1.1*pow((volBondScaled),-0.53);
          const double Bparam = (-0.148*log(volBondScaled)-0.96)*contactAngleEff*contactAngleEff - 0.0082*log(volBondScaled) + 0.48;
          const double Cparam = 0.0018*log(volBondScaled)+0.078;
          const double Fcapilary = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Aparam*dist/R2+Bparam)+Cparam);

          // calculate vn and vt since not in struct
          const double rinv = 1.0 / r;
          const double dx = scdata.delta[0];
          const double dy = scdata.delta[1];
          const double dz = scdata.delta[2];
          const double enx = dx * rinv;
          const double eny = dy * rinv;
          const double enz = dz * rinv;

          // relative translational velocity
          const double vr1 = v[i][0] - v[j][0];
          const double vr2 = v[i][1] - v[j][1];
          const double vr3 = v[i][2] - v[j][2];

          // normal component
          const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
          const double vn1 = vn * enx;
          const double vn2 = vn * eny;
          const double vn3 = vn * enz;

          // tangential component
          const double vt1 = vr1 - vn1;
          const double vt2 = vr2 - vn2;
          const double vt3 = vr3 - vn3;

          // relative rotational velocity
          double wr1, wr2, wr3;
          double const *omega_i = atom->omega[i];
          double const *omega_j = atom->omega[j];

          if(scdata.is_wall) {
            wr1 = radi * omega_i[0] * rinv;
            wr2 = radi * omega_i[1] * rinv;
            wr3 = radi * omega_i[2] * rinv;
          } else {
            wr1 = (radi * omega_i[0] + radj * omega_j[0]) * rinv;
            wr2 = (radi * omega_i[1] + radj * omega_j[1]) * rinv;
            wr3 = (radi * omega_i[2] + radj * omega_j[2]) * rinv;
          }

          // relative velocities
          const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
          const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
          const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

          // viscous force
          // this is from Nase et al as cited in Shi and McCarthy, Powder Technology, 184 (2008), 65-75, Eqns 40,41
          const double stokesPreFactor = -6.*M_PI*fluidViscosity*rEff;
          const double FviscN = stokesPreFactor*vn/std::max(minSeparationDistanceRatio,dist/rEff);
          const double FviscT_over_vt = (/* 8/15 */ 0.5333333*log(1./std::max(minSeparationDistanceRatio,dist/rEff)) + 0.9588) * stokesPreFactor;

          // tangential force components
          const double Ft1 = FviscT_over_vt*vtr1;
          const double Ft2 = FviscT_over_vt*vtr2;
          const double Ft3 = FviscT_over_vt*vtr3;

          // torques
          const double tor1 = eny * Ft3 - enz * Ft2;
          const double tor2 = enz * Ft1 - enx * Ft3;
          const double tor3 = enx * Ft2 - eny * Ft1;

          // add to fn, Ft
          //if(tangentialReduce_) scdata.Fn += Fcapilary+FviscN;
          //scdata.Ft += ...

          // apply normal and tangential force
          const double fx = (Fcapilary+FviscN) * enx + Ft1;
          const double fy = (Fcapilary+FviscN) * eny + Ft2;
          const double fz = (Fcapilary+FviscN) * enz + Ft3;

          scdata.has_force_update = true;

          // return resulting forces
          if(scdata.is_wall) {
            const double area_ratio = scdata.area_ratio;
            i_forces.delta_F[0] += fx * area_ratio;
            i_forces.delta_F[1] += fy * area_ratio;
            i_forces.delta_F[2] += fz * area_ratio;
            i_forces.delta_torque[0] += -radi * tor1 * area_ratio;
            i_forces.delta_torque[1] += -radi * tor2 * area_ratio;
            i_forces.delta_torque[2] += -radi * tor3 * area_ratio;
          } else {
            i_forces.delta_F[0] += fx;
            i_forces.delta_F[1] += fy;
            i_forces.delta_F[2] += fz;
            i_forces.delta_torque[0] += -radi * tor1; // using radius here, not contact radius
            i_forces.delta_torque[1] += -radi * tor2;
            i_forces.delta_torque[2] += -radi * tor3;

            j_forces.delta_F[0] -= fx;
            j_forces.delta_F[1] -= fy;
            j_forces.delta_F[2] -= fz;
            j_forces.delta_torque[0] += -radj * tor1; // using radius here, not contact radius
            j_forces.delta_torque[1] += -radj * tor2;
            j_forces.delta_torque[2] += -radj * tor3;
          }
      }
      // case (iii)
      else if(bridge_breaks)
      {
          if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

          // store for next step
          contflag[0] = 0.0;
          if (!scdata.is_wall)
          {
              // liquid transfer happens here
              // assume liquid distributes evenly
              double *liquidFlux = fix_liquidflux->vector_atom;
              
              const double invdt = 1./update->dt;
              const double rad_ratio = radj/radi;
              const double split_factor = 1.0/(1.0+rad_ratio*rad_ratio*rad_ratio);
              // liquid flux is in vol% per time
              liquidFlux[i] += invdt*(split_factor * volBond1000 - volLiBond1000) / (1333.333333*M_PI*radi*radi*radi) ;

              if (force->newton_pair || j < atom->nlocal)
                liquidFlux[j] += invdt*((1.-split_factor) * volBond1000 - volLjBond1000) / (1333.333333*M_PI*radj*radj*radj) ;
              
          }
      }
      // no else here, case (i) was already caught before
    }

  private:
    double surfaceLiquidContentInitial, surfaceTension, *contactAngle;
    double minSeparationDistanceRatio, maxSeparationDistanceRatio, fluidViscosity;
    double ln1overMinSeparationDistanceRatio;
    int history_offset;
    FixPropertyAtom *fix_surfaceliquidcontent;
    FixPropertyAtom *fix_liquidflux;
    FixScalarTransportEquation *fix_ste;
    bool tangentialReduce_;
  };
}
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
