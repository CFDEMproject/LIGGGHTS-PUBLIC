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
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_MULTICONTACT,multicontact,2)
#else
#ifndef SURFACE_MODEL_MULTICONTACT_H_
#define SURFACE_MODEL_MULTICONTACT_H_
#include "contact_models.h"
#include "surface_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "fix_multicontact_halfspace.h"

namespace LIGGGHTS {
namespace ContactModels
{
    template<>
    class SurfaceModel<SURFACE_MULTICONTACT> : public SurfaceModelBase
    {
      private:
        int delta_offset;
        class ContactModelBase *cmb;

      public:
        SurfaceModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) :
            SurfaceModelBase(lmp, hsetup, c),
            cmb(c)
        {
            if (cmb->is_wall()) {
                delta_offset = hsetup->add_history_value("surfPos_x", "0"); // position of wall contact
                hsetup->add_history_value("surfPos_y", "0");
                hsetup->add_history_value("surfPos_z", "0");
            } else {
                delta_offset = hsetup->add_history_value("radij", "0");     // radius(i) + expansion factor acting on particle i at contact ij
                hsetup->add_history_value("radji", "0");                    // radius(j) + expansion factor acting on particle j at contact ji
            }
            hsetup->add_history_value("fn", "0");                           // normal force
            cmb->add_history_offset("delta", delta_offset);
        }

        inline void registerSettings(Settings&) {}
        inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}
        inline void connectToProperties(PropertyRegistry&)
        {
            FixMultiContactHalfSpace *fix_mc = static_cast<FixMultiContactHalfSpace*>(modify->find_fix_style("multicontact/halfspace", 0));
            if (!fix_mc)
                error->one(FLERR, "No suitable multicontact fix found.");

            //set neighbor contact_distance_factor here, assume it's at most one radius away
            neighbor->register_contact_dist_factor(2.0);
        }

        inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
        {
            #ifdef SUPERQUADRIC_ACTIVE_FLAG
                sidata.is_non_spherical = false;
            #endif
            return true;
        }

        inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&)
        {
            const bool update_history = sidata.computeflag && sidata.shearupdate;
            if (update_history)
            {
                if (sidata.is_wall)
                {
                    sidata.contact_history[delta_offset  ] = -sidata.delta[0];
                    sidata.contact_history[delta_offset+1] = -sidata.delta[1];
                    sidata.contact_history[delta_offset+2] = -sidata.delta[2];
                }
                else
                {
                    sidata.contact_history[delta_offset  ] = sidata.radi;
                    sidata.contact_history[delta_offset+1] = sidata.radj;
                }
            }
            #ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(sidata.is_non_spherical)
              error->one(FLERR,"Using default surface model for non-spherical particles!");
            #endif
            const double enx = sidata.en[0];
            const double eny = sidata.en[1];
            const double enz = sidata.en[2];

            // relative translational velocity
            const double vr1 = sidata.v_i[0] - sidata.v_j[0];
            const double vr2 = sidata.v_i[1] - sidata.v_j[1];
            const double vr3 = sidata.v_i[2] - sidata.v_j[2];

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
            const double deltan = sidata.radsum - sidata.r;
            const double dx = sidata.delta[0];
            const double dy = sidata.delta[1];
            const double dz = sidata.delta[2];
            const double rinv = sidata.rinv;
            double wr1, wr2, wr3;

            if(sidata.is_wall) {
                // in case of wall contact, r is the contact radius
                const double cr = sidata.radi - 0.5*sidata.deltan;
                wr1 = cr * sidata.omega_i[0] * rinv;
                wr2 = cr * sidata.omega_i[1] * rinv;
                wr3 = cr * sidata.omega_i[2] * rinv;
                sidata.cri = cr;
            } else {
                const double cri = sidata.radi - 0.5 * deltan;
                const double crj = sidata.radj - 0.5 * deltan;
                wr1 = (cri * sidata.omega_i[0] + crj * sidata.omega_j[0]) * rinv;
                wr2 = (cri * sidata.omega_i[1] + crj * sidata.omega_j[1]) * rinv;
                wr3 = (cri * sidata.omega_i[2] + crj * sidata.omega_j[2]) * rinv;
                sidata.cri = cri;
                sidata.crj = crj;
            }

            // relative velocities
            const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
            const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
            const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

            sidata.vn = vn;
            sidata.deltan = deltan;
            sidata.wr1 = wr1;
            sidata.wr2 = wr2;
            sidata.wr3 = wr3;
            sidata.vtr1 = vtr1;
            sidata.vtr2 = vtr2;
            sidata.vtr3 = vtr3;
            sidata.P_diss = 0.;
        }

        inline void endSurfacesIntersect(SurfacesIntersectData& sidata, TriMesh * mesh, double * const i_forces)
        {
            // add the normal force to the contact history
            const int fn_offset = sidata.is_wall ? 3 : 2;
            const bool update_history = sidata.computeflag && sidata.shearupdate;
            if (update_history)
                sidata.contact_history[delta_offset+fn_offset] = fabs(vectorDot3D(i_forces, sidata.en));
        }

        inline void surfacesClose(SurfacesCloseData& sidata, ForceData&, ForceData&)
        {
            const bool update_history = sidata.computeflag && sidata.shearupdate;
            if (update_history)
            {
                if (sidata.is_wall) {
                    const double lenDelta = vectorMag3D(sidata.delta);
                    const double rinv = lenDelta > 1e-15 ? -1.0/lenDelta : 0.0;
                    sidata.contact_history[delta_offset  ] = sidata.delta[0]*sidata.radi*rinv;
                    sidata.contact_history[delta_offset+1] = sidata.delta[1]*sidata.radi*rinv;
                    sidata.contact_history[delta_offset+2] = sidata.delta[2]*sidata.radi*rinv;
                    sidata.contact_history[delta_offset+3] = 0.0;
                } else {
                    sidata.contact_history[delta_offset  ] = sidata.radi;
                    sidata.contact_history[delta_offset+1] = sidata.radj;
                    sidata.contact_history[delta_offset+2] = 0.0;
                }
            }
        }

        void beginPass(SurfacesIntersectData&, ForceData&, ForceData&) {}
        void endPass(SurfacesIntersectData&, ForceData&, ForceData&) {}
        inline void tally_pp(double,int,int,int) {}
        inline void tally_pw(double,int,int,int) {}
    };
}
}
#endif // SURFACE_MODEL_MULTICONTACT_H_
#endif
