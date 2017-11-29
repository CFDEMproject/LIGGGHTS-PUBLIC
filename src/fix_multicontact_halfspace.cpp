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

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_multicontact_halfspace.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "update.h"
#include "atom.h"
#include "force.h"
#include "pair_gran.h"
#include "pair_gran_proxy.h"
#include "fix_contact_history_mesh.h"
#include "fix_wall_gran.h"
#include "comm.h"
#include "neighbor.h"
#include "fix_property_atom.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include <list>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMultiContactHalfSpace::FixMultiContactHalfSpace(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    pairgran_(0),
    geometric_prefactor(1.125)
{
    nevery = 1;

    int iarg = 3;
    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg],"geometric_prefactor") == 0) {
            if (narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'geometric_prefactor'");
            iarg++;
            geometric_prefactor = force->numeric(FLERR,arg[iarg++]);
            if (geometric_prefactor <= 0.0)
                error->fix_error(FLERR,this,"geometric_prefactor > 0 required");
            hasargs = true;
        } else if(strcmp(style,"multicontact/halfspace") == 0) {
            char *errmsg = new char[strlen(arg[iarg])+50];
            sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg]);
            error->fix_error(FLERR,this,errmsg);
            delete []errmsg;
        }
    }

    if(!force->pair_match("gran", 0))
        error->fix_error(FLERR,this,"Please use a granular pair style before using this fix");

    int max_type = atom->get_properties()->max_type();
    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
}

/* ---------------------------------------------------------------------- */

FixMultiContactHalfSpace::~FixMultiContactHalfSpace()
{
    history_vector.clear();
    contact_property_atom_vector.clear();
}

/* ---------------------------------------------------------------------- */

void FixMultiContactHalfSpace::post_create()
{
    // request contactproperty atom (wall) fixes to be used
    static_cast<PairGran*>(force->pair_match("gran", 0))->do_store_multicontact_data();

    int nwalls = modify->n_fixes_style("wall/gran");
    for(int iwall = 0; iwall < nwalls; iwall++) {
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",iwall));

        if(fwg->is_mesh_wall())
        {
            int n_meshes = fwg->n_meshes();
            for(int imesh = 0; imesh < n_meshes; imesh++)
                (fwg->mesh_list())[imesh]->createMulticontactData();
        }
        else
            fwg->createMulticontactData();
    }
}

/* ---------------------------------------------------------------------- */

void FixMultiContactHalfSpace::init()
{
    if (force->pair == NULL)
        error->fix_error(FLERR,this,"No pair style is defined");

    pairgran_ = (PairGranProxy*)force->pair_match("gran",0);

    // check if surface model multicontact is specified
    // this automatically also creates the required contact history
    if (!pairgran_->contact_match("surface", "multicontact"))
        error->fix_error(FLERR,this,"Surface model is not multicontact");

    if (!pairgran_)
        error->fix_error(FLERR,this,"No valid granular pair style found");

    history_vector.clear();
    contact_property_atom_vector.clear();

    int sumDelta_offset_pair = pairgran_->get_history_offset("delta");
    if (sumDelta_offset_pair < 0)
        error->fix_error(FLERR,this,"Internal error: need delta history offset");
    history_vector.push_back(HistoryData('p', (void*) pairgran_, sumDelta_offset_pair));
    FixContactPropertyAtom *cpa = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("multicontactData_"));
    if (!cpa)
        error->fix_error(FLERR,this,"Internal error: no contactproperty/atom fix found for granular pair");
    contact_property_atom_vector.push_back(cpa);

    int nwalls = modify->n_fixes_style("wall/gran");
    for(int iwall = 0; iwall < nwalls; iwall++)
    {
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",iwall));

        if (!fwg->contact_match("surface", "multicontact"))
            error->fix_error(FLERR,this,"Surface model of wall is not multicontact");

        if(fwg->is_mesh_wall())
        {
            
            int n_meshes = fwg->n_meshes();
            for(int imesh = 0; imesh < n_meshes; imesh++)
            {
                char fixid[200];
                const int sumDelta_offset_wall = fwg->get_history_offset("delta");
                if(sumDelta_offset_wall < 0)
                    error->fix_error(FLERR,this,"Internal error: need delta history data (for wall meshes)");
                sprintf(fixid,"tracker_%s",(fwg->mesh_list())[imesh]->id);
                history_vector.push_back(HistoryData('m', (void*)modify->find_fix_id(fixid), sumDelta_offset_wall));
                sprintf(fixid,"multicontactData_%s",(fwg->mesh_list())[imesh]->id);
                cpa = static_cast<FixContactPropertyAtom*>(modify->find_fix_id(fixid));
                if (!cpa)
                    error->fix_error(FLERR,this,"Internal error: no contactproperty/atom fix found for mesh wall");
                contact_property_atom_vector.push_back(cpa);
            }
        }
        else
        {
            
            char fixid[200];
            const int sumDelta_offset_wall = fwg->get_history_offset("delta");
            if(sumDelta_offset_wall < 0)
                error->fix_error(FLERR,this,"Internal error: need delta history data (for wall primitives)");
            sprintf(fixid,"history_%s",fwg->id);
            history_vector.push_back(HistoryData('w', (void*)modify->find_fix_id(fixid), sumDelta_offset_wall));
            sprintf(fixid,"multicontactData_%s",fwg->id);
            cpa = static_cast<FixContactPropertyAtom*>(modify->find_fix_id(fixid));
            if (!cpa)
                error->fix_error(FLERR,this,"Internal error: no contactproperty/atom fix found for mesh wall");
            contact_property_atom_vector.push_back(cpa);
        }
    }

}

/* ---------------------------------------------------------------------- */

int FixMultiContactHalfSpace::setmask()
{
    int mask = 0;
    mask |= PRE_FORCE | MIN_PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMultiContactHalfSpace::setup_pre_force(int)
{
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixMultiContactHalfSpace::pre_force(int)
{
    // Clear all data in the FixContactPropertyAtom structure
    std::vector<FixContactPropertyAtom*>::iterator contact_property_atom = contact_property_atom_vector.begin();
    (*contact_property_atom)->clear();
    contact_property_atom++;
    for (; contact_property_atom < contact_property_atom_vector.end(); contact_property_atom++)
        static_cast<FixContactPropertyAtomWall*>(*contact_property_atom)->clear();

    int i,j;

    int nlocal = atom->nlocal;
    const double* const* x = atom->x;
    int *tag = atom->tag;
    int *type  = atom->type;
    int *mask = atom->mask;
    double *radius = atom->radius;

    // STEP 1:
    // Compute surface postion and get f.n and save it in a contactpropertyatom structure
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 0; i < nlocal; i++) {

        if (mask[i] & groupbit) {

            // small value for force epsilon
            const double Yi = Y[type[i]-1];
            const double F_eps = Yi*radius[i]*radius[i]*1e-14;

            // loop over all: pair and walls
            std::vector<HistoryData>::iterator history;
            for (history = history_vector.begin(), contact_property_atom = contact_property_atom_vector.begin();
                 history < history_vector.end() && contact_property_atom < contact_property_atom_vector.end();
                 history++, contact_property_atom++) {
                // loop over all contacts of particle i
                const int n = history->get_npartners(i);
                //for (j = 0; j < history->get_npartners(i); j++) {
                for (j = 0; j < n; j++) {

                    // get data pointer from contact history
                    double * const deltaData = history->get_data_ptr(i, j);
                    double surfPos_ij[4], surfPos_ji[4];

                    // compute surface position of contact ij and ji
                    history->compute_surfPos(i, j, x, deltaData, surfPos_ij, surfPos_ji, F_eps);
                    // get the normal force acting on the particle
                    double fn = history->get_fn(deltaData);

                    surfPos_ij[3] = fn;
                    surfPos_ji[3] = fn;
                    // save the newly computed variables in a contactpropertyatom fix that will be used throughout the time step
                    history->save_contact_property_atom(i, j, tag, surfPos_ij, surfPos_ji, *contact_property_atom);
                }
            }
        }
    }

    // STEP 2:
    // Compute delta_ij from contacts ik (k!=j)
    /////////////////////////////////////////////////////////
    for (i = 0; i < nlocal; i++) {

        if (mask[i] & groupbit) {

            const double Yi = Y[type[i]-1];
            const double nui = nu[type[i]-1];
            const double C1 = -(1.0+nui)/(2.0*M_PI*Yi)*geometric_prefactor;
            const double C2 = 3.0 - 4.0*nui;
            const double C3 = 1.0 - 2.0*nui;
            // small value for force epsilon
            const double F_eps = Yi*radius[i]*radius[i]*1e-14;

            std::list<double> delta_list;
            std::list<double*> surfPos;
            std::list<double> invMagSurfPos;
            std::list<double>::iterator delta_ij, delta_ik;
            std::list<double*>::iterator surfPos_ij, surfPos_ik;
            std::list<double>::iterator invMagSurfPos_ij, invMagSurfPos_ik;
            // loop over all contacts ij and get pointers to surfPos
            for (contact_property_atom = contact_property_atom_vector.begin(); contact_property_atom < contact_property_atom_vector.end(); contact_property_atom++) {
                for (j = 0; j < (*contact_property_atom)->get_npartners(i); j++) {
                        surfPos.push_back((*contact_property_atom)->contacthistory(i,j));
                        invMagSurfPos.push_back(1.0/fmax(vectorMag3D(surfPos.back()),1e-6*radius[i]));
                        delta_list.push_back(0.0);
                }
            }

            // loop over all contacts ij
            for (surfPos_ij = surfPos.begin(), invMagSurfPos_ij = invMagSurfPos.begin(), delta_ij = delta_list.begin();
                 *surfPos_ij != surfPos.back();
                 surfPos_ij++, invMagSurfPos_ij++, delta_ij++) {

                // loop over all contacts ik
                surfPos_ik = surfPos_ij;
                invMagSurfPos_ik = invMagSurfPos_ij;
                delta_ik = delta_ij;
                do {
                    surfPos_ik++;
                    invMagSurfPos_ik++;
                    delta_ik++;

                    const double f_ij = (*surfPos_ij)[3];
                    const double f_ik = (*surfPos_ik)[3];

                    if (f_ij > F_eps || f_ik > F_eps) { // check if one of the force magnitudes is larger than the threshold
                        double ukc[3];
                        vectorSubtract3D(*surfPos_ij, *surfPos_ik, ukc);
                        const double invDkc = 1.0/vectorMag3D(ukc);
                        if (invDkc < 1e10) {
                            const double nk_dot_nc = vectorDot3D(*surfPos_ij, *surfPos_ik)*(*invMagSurfPos_ij)*(*invMagSurfPos_ik);
                            const double ukc_dot_nc = -vectorDot3D(ukc, *surfPos_ij)*(*invMagSurfPos_ij)*invDkc;
                            const double ukc_dot_nk = -vectorDot3D(ukc, *surfPos_ik)*(*invMagSurfPos_ik)*invDkc;
                            double careful_expression = C2*nk_dot_nc + ukc_dot_nk*ukc_dot_nc;
                            if (f_ik > F_eps) {
                                double c3term = 0.0;
                                if (fabs(C3) > 1e-6) {
                                    const double divisor = 1.0 + ukc_dot_nk;
                                    if (divisor > 1e-6)
                                        c3term = -C3*(nk_dot_nc + ukc_dot_nc)/divisor;
                                }
                                *delta_ij += C1*f_ik * (careful_expression + c3term) * invDkc;
                            }
                            if (f_ij > F_eps) {
                                double c3term = 0.0;
                                if (fabs(C3) > 1e-6) {
                                    const double divisor = 1.0 - ukc_dot_nc;
                                    if (divisor > 1e-6)
                                        c3term = -C3*(nk_dot_nc - ukc_dot_nk)/divisor;
                                }
                                *delta_ik += C1*f_ij * (careful_expression + c3term) * invDkc;
                            }
                        }
                    }
                } while (*surfPos_ik != surfPos.back());
            }

            // loop over all contacts ij to save the delta values instead of the f.n values that were saved in the contact property atom structure before
            // this needs to be done after in order to avoid a race condition
            delta_ij = delta_list.begin();
            surfPos_ij = surfPos.begin();
            // loop over all contacts ij
            for (contact_property_atom = contact_property_atom_vector.begin(); contact_property_atom < contact_property_atom_vector.end(); contact_property_atom++) {
                for (j = 0; j < (*contact_property_atom)->get_npartners(i); j++) {
                    // get the pointer to the contact history
                    // replace f.n in the contactpropertyatom structure with delta_ij
                    (*surfPos_ij)[3] = *delta_ij;
                    // advance delta list iterator
                    delta_ij++;
                    surfPos_ij++;
                }
            }
            // clear the delta list so that it can be refilled for the next grain
            delta_list.clear();
        }
    }

    // STEP 3:
    // Communicate the contactpropertyatom data to all other processors
    /////////////////////////////////////////////////////////////////////////////////////////
    for (contact_property_atom = contact_property_atom_vector.begin(); contact_property_atom < contact_property_atom_vector.end(); contact_property_atom++)
        (*contact_property_atom)->do_forward_comm();
}

const int HistoryData::get_npartners(const int i)
{
    if (type == 'p') {
        NeighList* list = static_cast<PairGranProxy*>(fix_ptr)->list;
        return list->numneigh[i];
    } else if (type == 'm') {
        return static_cast<FixContactHistoryMesh*>(fix_ptr)->n_partner(i);
    } else { // == 'w'
        return 1; // primitives can only have one contact per particle
    }
}

double * const HistoryData::get_data_ptr(const int i, const int jj)
{
    if (type == 'p') {
        NeighList* list = static_cast<PairGranProxy*>(fix_ptr)->list;
        double * const all_contact_list = list->listgranhistory->firstdouble[i];
        return &all_contact_list[jj*list->listgranhistory->dnum + offset];
    } else if (type == 'm') {
        const int j = static_cast<FixContactHistoryMesh*>(fix_ptr)->get_contact(i, jj);
        return static_cast<FixContactHistoryMesh*>(fix_ptr)->contacthistory(i,j) + sizeof(double)*offset;
    } else { // == 'w'
        return static_cast<FixPropertyAtom*>(fix_ptr)->array_atom[i] + sizeof(double)*offset;
    }
}

void HistoryData::compute_surfPos(const int i, const int jj, const double * const *x, const double * const data_ptr, double * const surfPos_ij, double * const surfPos_ji, const double F_eps)
{
    if (type == 'p') {
        const NeighList * const list = static_cast<PairGranProxy*>(fix_ptr)->list;
        const int j = list->firstneigh[i][jj];
        // in the following the following computations are performed:
        // normal of contact
        // n = (x_j - x_i) / |(x_j - x_i)|
        // closest point of sphere_i to midpoint_j:
        // cp_i = x[i] + n * (r_i + delta_ij)
        // closest point of sphere_j to midpoint_i:
        // cp_j = x[j] - n * (r_j + delta_ji)
        // actual contact is at the midpoint:
        // midPoint = (cp_i + cp_j)/2
        // relative surface Position:
        // surfPos_ij = midPoint - x_i
        // surfPos_ji = midPoint - x_j
        //
        // the following code uses a few algebraic reformulations for optimization reasons
        
        double xij[3];
        vectorSubtract3D(x[j], x[i], &(xij[0]));
        const double length_xij = vectorMag3D(xij);
        if (data_ptr[2] > F_eps) {
            double tmp[3];
            vectorScalarMult3D(xij, (data_ptr[0] - data_ptr[1])*0.5/length_xij, tmp);
            vectorAddMultiple3D(tmp,  0.5, xij, surfPos_ij);
            vectorAddMultiple3D(tmp, -0.5, xij, surfPos_ji);
        } else {
            // no direct contact
            // surface position will be (r_i + delta_ij)*n
            vectorScalarMult3D(xij,  data_ptr[0]/length_xij, surfPos_ij);
            vectorScalarMult3D(xij, -data_ptr[1]/length_xij, surfPos_ji);
        }
    } else {
        vectorCopy3D(data_ptr, surfPos_ij);
    }
}

double HistoryData::get_fn(const double * const data_ptr)
{
    if (type == 'p')
        return data_ptr[2];
    else
        return data_ptr[3];
}

void HistoryData::save_contact_property_atom(const int i, const int jj, const int * const tag, const double * const surfPos_ij, const double * const surfPos_ji, FixContactPropertyAtom *cpa)
{
    if (type == 'p') {
        const NeighList * const list = static_cast<PairGranProxy*>(fix_ptr)->list;
        const int j = list->firstneigh[i][jj];
        cpa->add_partner(i, tag[j], surfPos_ij);
        cpa->add_partner(j, tag[i], surfPos_ji);
    } else if (type == 'm') {
        const FixContactHistoryMesh * const fix_history = static_cast<FixContactHistoryMesh*>(fix_ptr);
        const int idTri = fix_history->get_partner_idTri(i, jj);
        cpa->add_partner(i, idTri, surfPos_ij);
    } else {
        // only add if we have a surfpos != 0
        const double lenSurfPos = vectorMag3DSquared(surfPos_ij);
        if (lenSurfPos > 1e-14) {
            cpa->add_partner(i, 1, surfPos_ij);
        }
    }
}
