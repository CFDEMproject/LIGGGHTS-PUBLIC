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
COHESION_MODEL(COHESION_WASHINO_CAPILLARY_VISCOUS,washino/capillary/viscous,7)
#else

#ifndef COHESION_MODEL_WASHINO_CAPILLARY_VISCOUS_H_
#define COHESION_MODEL_WASHINO_CAPILLARY_VISCOUS_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include <algorithm>
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "mesh_module_liquidtransfer.h"

namespace MODEL_PARAMS
{
    inline static ScalarProperty* createliquidContentInitialWashino(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* surfaceLiquidContentInitialScalar = MODEL_PARAMS::createScalarProperty(registry, "surfaceLiquidContentInitial", caller);
      return surfaceLiquidContentInitialScalar;
    }

    inline static ScalarProperty* createMinSeparationDistanceRatioWashino(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* minSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDistanceRatio", caller);
      return minSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createMaxSeparationDistanceRatioWashino(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* maxSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistanceRatio", caller);
      return maxSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createFluidViscosityWashino(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
      return fluidViscosityScalar;
    }

    inline static VectorProperty* createMaxLiquidContent(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
        return createPerTypeProperty(registry, "maxLiquidContent", caller, sanity_checks, 0.0, 1.0);
    }

    inline static ScalarProperty* createLbVolumeFraction(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
        return createScalarProperty(registry, "lbVolumeFraction", caller, sanity_checks, 0.0, 1.0);
    }

}

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_WASHINO_CAPILLARY_VISCOUS> : public CohesionModelBase {

  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      surfaceLiquidContentInitial(0.0),
      surfaceTension(0.0),
      contactAngle(0),
      minSeparationDistanceRatio(0.0),
      maxSeparationDistanceRatio(0.0),
      fluidViscosity(0.),
      ln1overMinSeparationDistanceRatio(0.0),
      maxLiquidContent(0),
      volumeFraction(0.05),
      history_offset(0),
      fix_surfaceliquidcontent(0),
      fix_liquidflux(0),
      fix_ste(0),
      limit_lqc_flag_(false),
      mod_lb_vol_flag_(false)
    {
      history_offset = hsetup->add_history_value("contflag", "0");
      
    }

    void registerSettings(Settings & settings)
    {
        settings.registerOnOff("limitLiquidContent", limit_lqc_flag_, false);
        settings.registerOnOff("modifyLbVolume", mod_lb_vol_flag_, false);
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("surfaceLiquidContentInitial", &MODEL_PARAMS::createliquidContentInitialWashino);
      registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
      registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityWashino);
      registry.registerProperty("contactAngle", &MODEL_PARAMS::createContactAngle);
      registry.registerProperty("minSeparationDistanceRatio", &MODEL_PARAMS::createMinSeparationDistanceRatioWashino);
      registry.registerProperty("maxSeparationDistanceRatio", &MODEL_PARAMS::createMaxSeparationDistanceRatioWashino);

      registry.connect("surfaceLiquidContentInitial", surfaceLiquidContentInitial,"cohesion_model washino/capillary/viscous");
      registry.connect("surfaceTension", surfaceTension,"cohesion_model washino/capillary/viscous");
      registry.connect("fluidViscosity", fluidViscosity,"cohesion_model washino/capillary/viscous");
      registry.connect("contactAngle", contactAngle,"cohesion_model washino/capillary/viscous");
      registry.connect("minSeparationDistanceRatio", minSeparationDistanceRatio,"cohesion_model washino/capillary/viscous");
      
      registry.connect("maxSeparationDistanceRatio", maxSeparationDistanceRatio,"cohesion_model washino/capillary/viscous");

      ln1overMinSeparationDistanceRatio = log(1./minSeparationDistanceRatio);

      // if limit_lqc_flag_ need additional material property
      if (limit_lqc_flag_) {
          registry.registerProperty("maxLiquidContent", &MODEL_PARAMS::createMaxLiquidContent);
          registry.connect("maxLiquidContent", maxLiquidContent, "cohesion_model washino/capillary/viscous");
      }

      // if mod_lb_vol_flag_ need additional property
      if (mod_lb_vol_flag_) {
          registry.registerProperty("lbVolumeFraction", &MODEL_PARAMS::createLbVolumeFraction);
          registry.connect("lbVolumeFraction", volumeFraction, "cohesion_model washino/capillary/viscous");
      }

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
        //fix_liquidcontent->set_all(liquidContentInitial);
      }

      //TODO-AM require initial liquid content of wall:
      //requires a thickness parameter

      fix_surfaceliquidcontent = static_cast<FixPropertyAtom*>(modify->find_fix_property("surfaceLiquidContent","property/atom","scalar",0,0,"cohesion_model washino/capillary/viscous"));
      fix_liquidflux = static_cast<FixPropertyAtom*>(modify->find_fix_property("liquidFlux","property/atom","scalar",0,0,"cohesion_model washino/capillary/viscous"));
      fix_ste = modify->find_fix_scalar_transport_equation("liquidtransfer");

      if(!fix_surfaceliquidcontent || !fix_liquidflux || !fix_ste)
          error->all(FLERR,"internal error");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model washino/capillary/viscous");

      if (limit_lqc_flag_) {
          const int max_type = registry.max_type();
          const double max_rad = registry.max_radius();
          const double min_rad = registry.min_radius();
          double max_dist_ratio = 0.0;
          for (int i = 1; i <= max_type; i++)
          {
              const double volL1000 = /* 2*4/3 * 1000 */ 2666.666666*M_PI*max_rad*max_rad*max_rad*maxLiquidContent[i];
              const double volBond1000 = (volL1000)*volumeFraction;
              const double contactAngleI = 0.5 * contactAngle[i] * contactAngle[i];
              const double distMax = (1. + 0.5*contactAngleI) * cbrt(volBond1000) * 0.1 /* 0.1*cbrt(1000)=1 */;
              max_dist_ratio = fmax(0.5*distMax/min_rad, max_dist_ratio); 
          }
          if (screen) fprintf(screen, "Warning: maxLiquidContent was specified, resulting in maxSeparationDistanceRatio being overwritten by %e (was %e)\n", 1.0 + max_dist_ratio, maxSeparationDistanceRatio);
          if (logfile) fprintf(logfile, "Warning: maxLiquidContent was specified, resulting in maxSeparationDistanceRatio being overwritten by %e (was %e)\n", 1.0 + max_dist_ratio, maxSeparationDistanceRatio);
          maxSeparationDistanceRatio = 1.0 + max_dist_ratio;
      }

      neighbor->register_contact_dist_factor(1.1*maxSeparationDistanceRatio);
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
      const double radj = sidata.radj;
      const double r = sidata.r;
      const double dist =  sidata.is_wall ? r - radi : r - (radi + radj);
      double const *surfaceLiquidContent = fix_surfaceliquidcontent->vector_atom;

      ScalarContainer<double> *liquidContPtr = sidata.is_wall && sidata.mesh ? sidata.mesh->prop().getElementProperty< ScalarContainer<double> >("LiquidContent") : NULL;

      MeshModuleLiquidTransfer *mm_liquid_transfer = liquidContPtr ? static_cast<MeshModuleLiquidTransfer*>((static_cast<FixMeshSurface*>(sidata.fix_mesh))->get_module("liquidtransfer")) : NULL;
      const double wallThickness = liquidContPtr ? mm_liquid_transfer->get_wall_thickness() : 0.0;

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
      double * const contflag = &sidata.contact_history[history_offset];
      // store for noCollision
      contflag[0] = 1.0;

      // limit maximum liquid content
      if (limit_lqc_flag_)
      {
          limitLiquidContent(i,itype);
          if (!sidata.is_wall)
            limitLiquidContent(j,jtype);
      }

      const double volLi1000 = /* 4/3 * 1000 */ 1333.333333*M_PI*radi*radi*radi*surfaceLiquidContent[i];
      const double volLj1000 = /* 4/3 * 1000 */ sidata.is_wall ?
        (liquidContPtr ? (*liquidContPtr)(sidata.j)*fmin(radi*radi*M_PI,sidata.mesh->areaElem(sidata.j))*1000.0*wallThickness : 0.0) :
        1333.333333*M_PI*radj*radj*radj*surfaceLiquidContent[j]; // TODO-wall liquid content at wall?
      const double volBond1000 = (volLi1000+volLj1000)*volumeFraction;

      // skip if bond volume too small
      if(volBond1000 < 1e-14) return;

      const double rEff = radi*radj / (radi+radj);
      const double contactAngleEff = 0.5 * contactAngle[itype] * contactAngle[jtype];

      // capilar force
      // this is from Rabinovich et al., Langmiur, 21 (2005), 10992-10997 - Eqn. A11
      // separation distance = 0 in this case
      const double Fcapilary = - 2.*M_PI*rEff*surfaceTension*cos(contactAngleEff)*(1.-dist*sqrt(1000.*M_PI*rEff/(2.*volBond1000)));

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
      // for this model, bridge is created and breaks at rupture distancy

      const int i = scdata.i;
      const int j = scdata.j;
      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double radi = scdata.radi;
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);
      const double dist =  scdata.is_wall ? r - radi : r - (radi + radj);
      double const *surfaceLiquidContent = fix_surfaceliquidcontent->vector_atom;

      ScalarContainer<double> *liquidContPtr = scdata.is_wall && scdata.mesh ? scdata.mesh->prop().getElementProperty< ScalarContainer<double> >("LiquidContent") : NULL;

      MeshModuleLiquidTransfer *mm_liquid_transfer = liquidContPtr ? static_cast<MeshModuleLiquidTransfer*>((static_cast<FixMeshSurface*>(scdata.fix_mesh))->get_module("liquidtransfer")) : NULL;
      const double wallThickness = liquidContPtr ? mm_liquid_transfer->get_wall_thickness() : 0.0;
      const double wallArea = liquidContPtr ? scdata.mesh->areaElem(scdata.j) : 0.0;

      // limit maximum liquid content
      if (limit_lqc_flag_)
      {
          limitLiquidContent(i,itype);
          if (!scdata.is_wall) limitLiquidContent(j,jtype);
      }

      const double volLi1000 = /* 4/3 * 1000 */ 1333.333333*M_PI*radi*radi*radi*surfaceLiquidContent[i];
      const double volLj1000 = /* 4/3 * 1000 */ scdata.is_wall ?
        (liquidContPtr ? (*liquidContPtr)(scdata.j)*fmin(radi*radi*M_PI, wallArea)*1000.0*wallThickness : 0.0) :
        1333.333333*M_PI*radj*radj*radj*surfaceLiquidContent[j]; // TODO-wall liquid content at wall?
      const double volBond1000 = (volLi1000+volLj1000)*volumeFraction;

      const double rEff = radi*radj / (radi+radj);
      const double contactAngleEff = 0.5 * contactAngle[itype] * contactAngle[jtype];
      const double distMax = (1. + 0.5*contactAngleEff) * cbrt(volBond1000) * 0.1 /* 0.1*cbrt(1000)=1 */;

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

          // store for next step
          contflag[0] = 1.0;

          // skip if bond volume too small
          if(volBond1000 < 1e-14) return;

          // calculate forces

          // capilary force, case no collision
          // this is from Langmiur, 21 (24), 2005 - Eqns. 19, 20, A3
          const double prefactor = -1. + sqrt(1.+2.*volBond1000/(M_PI*rEff*1000.*dist*dist));
          const double dSpSp =       0.5*dist*prefactor;
          const double alpha = sqrt(dist/rEff*prefactor);
          const double Fcapilary = - 2.*M_PI*rEff*surfaceTension* (cos(contactAngleEff) / (1. + dist/(2.*dSpSp)) + sin(alpha)*sin(alpha+contactAngleEff));

          // calculate vn and vt since not in struct
          const double rinv = 1.0 / r;
          const double dx = scdata.delta[0];
          const double dy = scdata.delta[1];
          const double dz = scdata.delta[2];
          const double enx = dx * rinv;
          const double eny = dy * rinv;
          const double enz = dz * rinv;

          // relative translational velocity
          const double vr1 = scdata.v_i[0] - scdata.v_j[0];
          const double vr2 = scdata.v_i[1] - scdata.v_j[1];
          const double vr3 = scdata.v_i[2] - scdata.v_j[2];

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
          double wr1 = 0.0, wr2 = 0.0, wr3 = 0.0;
          double const *omega_i = atom->omega[i];

          if(!scdata.is_wall) {
            double const *omega_j = atom->omega[j];
            wr1 = radj * omega_j[0];
            wr2 = radj * omega_j[1];
            wr3 = radj * omega_j[2];
          }
          wr1 = (radi * omega_i[0] + wr1)*rinv;
          wr2 = (radi * omega_i[1] + wr2)*rinv;
          wr3 = (radi * omega_i[2] + wr3)*rinv;

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

          // liquid transfer happens here
          if (!scdata.is_wall)
          {
              // assume liquid distributes evenly
              double *liquidFlux = fix_liquidflux->vector_atom;
              
              const double invdt = 1./update->dt;
              const double rad_ratio = radj/radi;
              const double split_factor = 1.0/(1.0+rad_ratio*rad_ratio*rad_ratio);
              // liquid flux is in vol% per time
              liquidFlux[i] += invdt*(split_factor * volBond1000 - volumeFraction*volLi1000) / (1333.333333*M_PI*radi*radi*radi) ;

              if (force->newton_pair || j < atom->nlocal)
                  liquidFlux[j] += invdt*((1.-split_factor) * volBond1000 - volumeFraction*volLj1000) / (1333.333333*M_PI*radj*radj*radj) ;
              
          }
          else if (liquidContPtr)
          {
              // assume liquid distributes evenly
              double *liquidFlux = fix_liquidflux->vector_atom;
              const double invdt = 1./update->dt;
              // liquid flux in vol per time
              double volFlux = invdt*(0.5 * volBond1000 - volumeFraction*volLi1000);
              // liquid flux in vol% per time for particle
              liquidFlux[i] += volFlux/(1333.333333*M_PI*radi*radi*radi);
              // TODO assume partArea < wallArea (overlap detection)
              const double wallLiquidFlux = -volFlux/(1000.0*wallThickness*wallArea);
              mm_liquid_transfer->add_liquid_flux(scdata.j, wallLiquidFlux, limit_lqc_flag_, limit_lqc_flag_ ? maxLiquidContent[scdata.jtype] : 0.0);
          }
      }
      // no else here, case (i) was already caught before
    }

  private: // private functions
    inline void limitLiquidContent(const int idx, const int itype)
    {
        double *surfaceLiquidContent = fix_surfaceliquidcontent->vector_atom;

        if (surfaceLiquidContent[idx] > maxLiquidContent[itype]) {
            
            surfaceLiquidContent[idx] = maxLiquidContent[itype];
        }
    }

  private:

    double surfaceLiquidContentInitial, surfaceTension, *contactAngle;
    double minSeparationDistanceRatio, maxSeparationDistanceRatio, fluidViscosity;
    double ln1overMinSeparationDistanceRatio, *maxLiquidContent;
    double volumeFraction;
    int history_offset;
    FixPropertyAtom *fix_surfaceliquidcontent;
    FixPropertyAtom *fix_liquidflux;
    FixScalarTransportEquation *fix_ste;

    bool limit_lqc_flag_;
    bool mod_lb_vol_flag_;

    bool tangentialReduce_;

    MeshModuleLiquidTransfer *mm_liquid_transfer;
  };
}
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
