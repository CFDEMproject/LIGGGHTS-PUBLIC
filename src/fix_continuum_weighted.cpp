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
#include "fix_continuum_weighted.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "pair_gran.h"
#include "pair_gran_proxy.h"
#include "fix_wall_gran.h"
#include "force.h"
#include "fix_property_atom.h"
#include "fix_contact_property_atom.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContinuumWeighted::FixContinuumWeighted(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    kernel_radius_(0.0),
    kernel_sqRadius_(0.0),
    fix_stress_(NULL),
    fix_strain_(NULL),
    fix_cont_vars_(NULL),
    fix_contact_forces_(NULL),
    compute_stress(false),
    compute_strain(false),
    kernel_type(TOP_HAT)
{
    nevery = 1;

    int iarg = 3;
    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg],"kernel_radius") == 0) {
            if (narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'kernel_radius'");
            iarg++;
            kernel_radius_ = force->numeric(FLERR,arg[iarg++]);
            kernel_sqRadius_ = kernel_radius_*kernel_radius_;
            if (kernel_radius_ <= 0.0)
                error->fix_error(FLERR,this,"kernel_radius > 0 required");
            hasargs = true;
        } else if (strcmp(arg[iarg],"kernel_type") == 0) {
            if (narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'kernel_type'");
            iarg++;
            if (strcmp(arg[iarg],"Top_Hat") == 0)
                kernel_type = TOP_HAT;
            else if (strcmp(arg[iarg],"Gaussian") == 0)
                kernel_type = GAUSSIAN;
            else if (strcmp(arg[iarg],"Wendland") == 0)
                kernel_type = WENDLAND;
            else
                error->fix_error(FLERR,this,"Unknown kernel_type parameter");
            iarg++;
            hasargs = true;
        } else if (strcmp(arg[iarg],"compute") == 0) {
            if (narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'compute'");
            iarg++;
            if (strcmp(arg[iarg],"stress") == 0)
                compute_stress = true;
            else if (strcmp(arg[iarg],"strain") == 0)
                compute_strain = true;
            else if (strcmp(arg[iarg],"stress_strain") == 0) {
                compute_stress = true;
                compute_strain = true;
            } else
                error->fix_error(FLERR,this,"Unknown compute parameter");
            iarg++;
            hasargs = true;
        } else if (strcmp(style,"continuum/weighted") == 0) {
            char *errmsg = new char[strlen(arg[iarg])+50];
            sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg]);
            error->fix_error(FLERR,this,errmsg);
            delete []errmsg;
        }
    }

    if (kernel_sqRadius_ <= 0.0)
        error->fix_error(FLERR, this, "Please provide a kernel_radius > 0");

    if (!(compute_strain || compute_stress))
        error->fix_error(FLERR, this, "Please provide at least one compute target (stress, strain or stress_strain)");

    if (compute_strain && kernel_type == TOP_HAT)
        error->fix_error(FLERR, this, "Strain cannot be computed using the TOP_HAT kernel, please select another kernel");

    if(!force->pair_match("gran", 0))
        error->fix_error(FLERR,this,"Please use a granular pair style before using this fix");
    (static_cast<PairGran*>(force->pair_match("gran", 0)))->do_store_contact_forces_stress();
}

/* ---------------------------------------------------------------------- */

FixContinuumWeighted::~FixContinuumWeighted()
{ }

/* ---------------------------------------------------------------------- */

void FixContinuumWeighted::post_create()
{
    if(!fix_stress_) {
        const char * fixarg[17];
        fixarg[0]="stressTensor_";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="stressTensor_";
        fixarg[4]="vector"; 
        fixarg[5]="yes";    
        fixarg[6]="no";     
        fixarg[7]="yes";    
        fixarg[8]="1.0";    // sigma_xx
        fixarg[9]="0.0";    // sigma_xy
        fixarg[10]="0.0";   // sigma_xz
        fixarg[11]="0.0";   // sigma_yx
        fixarg[12]="1.0";   // sigma_yy
        fixarg[13]="0.0";   // sigma_yz
        fixarg[14]="0.0";   // sigma_zx
        fixarg[15]="0.0";   // sigma_zy
        fixarg[16]="1.0";   // sigma_zz
        fix_stress_ = modify->add_fix_property_atom(17,const_cast<char**>(fixarg),style);
    }
    if(!fix_strain_) {
        const char * fixarg[17];
        fixarg[0]="strainTensor_";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="strainTensor_";
        fixarg[4]="vector"; 
        fixarg[5]="yes";    
        fixarg[6]="no";     
        fixarg[7]="yes";    
        fixarg[8]="0.0";    // epsilon_xx
        fixarg[9]="0.0";    // epsilon_xy
        fixarg[10]="0.0";   // epsilon_xz
        fixarg[11]="0.0";   // epsilon_yx
        fixarg[12]="0.0";   // epsilon_yy
        fixarg[13]="0.0";   // epsilon_yz
        fixarg[14]="0.0";   // epsilon_zx
        fixarg[15]="0.0";   // epsilon_zy
        fixarg[16]="0.0";   // epsilon_zz
        fix_strain_ = modify->add_fix_property_atom(17,const_cast<char**>(fixarg),style);
    }
    if(!fix_cont_vars_) {
        const char * fixarg[15];
        fixarg[0]="cont_vars_";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="cont_vars_";
        fixarg[4]="vector"; 
        fixarg[5]="yes";    
        fixarg[6]="yes";    
        fixarg[7]="no";     
        fixarg[8]="0.0";    // p_x
        fixarg[9]="0.0";    // p_y
        fixarg[10]="0.0";   // p_z
        fixarg[11]="0.0";   // rho
        fixarg[12]="0.0";   // grad_rho_x
        fixarg[13]="0.0";   // grad_rho_y
        fixarg[14]="0.0";   // grad_rho_z
        fix_cont_vars_ = modify->add_fix_property_atom(15,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixContinuumWeighted::init()
{
    pairgran_ = (PairGran*)force->pair_match("gran",0);

    // contact forces between a pair
    fix_contact_forces_ = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("contactforces_stress_"));
    if(!fix_contact_forces_)
        error->fix_error(FLERR,this,"Internal error: need fix contactforces_stress_");

    // contact forces between particle & walls:
    fix_wall_contact_forces_vector_.clear();

    int nwalls = modify->n_fixes_style("wall/gran");
    for(int iwall = 0; iwall < nwalls; iwall++)
    {
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",iwall));

        if(!fwg->store_force_contact_stress())
            error->fix_error(FLERR,this,"Internal error: contact forces for stress computation are not stored (make sure this fix is added before the wall fixes or use the explicit store_force_contact_stress option in those fixes)");

        if(fwg->is_mesh_wall())
        {
            
            int n_meshes = fwg->n_meshes();
            for(int imesh = 0; imesh < n_meshes; imesh++)
            {
                char fixid[200];
                sprintf(fixid,"contactforces_stress_%s",(fwg->mesh_list())[imesh]->id);

                FixContactPropertyAtomWall *fix_contact_wall = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
                if(!fix_contact_wall)
                    error->fix_error(FLERR,this,"Internal error: need fix contactforces");

                fix_wall_contact_forces_vector_.push_back(fix_contact_wall);
            }
        }
        else
        {
            
            char fixid[200];
            sprintf(fixid,"contactforces_stress_%s",fwg->id);

            FixContactPropertyAtomWall *fix_contact_wall = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
            if(!fix_contact_wall)
                error->fix_error(FLERR,this,"Internal error: need fix contactforces");

            fix_wall_contact_forces_vector_.push_back(fix_contact_wall);
        }
    }

}

/* ---------------------------------------------------------------------- */

int FixContinuumWeighted::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {

template<>
double FixContinuumWeighted::weightingFunction<TOP_HAT>(const double r)
{
    // 1/ (4/3 pi r^3)
    const double vol_sphere = kernel_sqRadius_*kernel_radius_*4.188790204786391;
    return 1.0/vol_sphere;
}

template<>
double FixContinuumWeighted::gradWeightingFunction<TOP_HAT>(const double r)
{
    // this should not be used
    return 0.0;
}

template<>
double FixContinuumWeighted::weightingFunction<GAUSSIAN>(const double r)
{
    // w = exp(-r^2/kernel_sqRadius_*9)
    // the constant 9 is coming from 3^2 which is the cutoff of the gaussian (at 3 sigma)
    const double w = exp(-r*r*9.0/kernel_sqRadius_);
    // phi = w(r)/[4 pi \int_0^inf w(s) s^2 ds]
    // constant below is (pi (e^9 sqrt(pi) erf(3)-6))/(27 e^9)
    const double normalization = kernel_sqRadius_*kernel_radius_*0.20614365813686698838756664006357;
    return w/normalization;
}

template<>
double FixContinuumWeighted::gradWeightingFunction<GAUSSIAN>(const double r)
{
    // grad_w = dw/dr = -2/(kernel_sqRadius_/9) exp(-r^2/kernel_sqRadius_*9) * r
    // the multiplication with r is happening outside of this function (as x_ij is used)
    // the constant 9 is coming from 3^2 which is the cutoff of the gaussian (at 3 sigma)
    const double grad_w = exp(-r*r*9.0/kernel_sqRadius_)*(-18.0)/kernel_sqRadius_;
    // grad phi = grad_w(r)/[4 pi \int_0^inf w(s) s^2 ds]
    // constant below is (pi (e^9 sqrt(pi) erf(3)-6))/(27 e^9)
    const double normalization = kernel_sqRadius_*kernel_radius_*0.20614365813686698838756664006357;
    return grad_w/normalization;
}

template<>
double FixContinuumWeighted::weightingFunction<WENDLAND>(const double r)
{
    // q = 2 r / kernel_radius
    const double q = 2.0*r/kernel_radius_;
    // w = (1-q/2)^4 (1+2q)
    double w = (1.0 - q*0.5);
    w *= w;
    w *= w;
    w *= (1.0 + 2.0*q);
    // normalization: (2 pi R^3)/21, constant below is (2pi)/21
    const double normalization = kernel_sqRadius_*kernel_radius_*0.2991993003418850703;
    return w/normalization;
}

template<>
double FixContinuumWeighted::gradWeightingFunction<WENDLAND>(const double r)
{
    // q = 2 r / kernel_radius
    const double q = 2.0*r/kernel_radius_;
    // w = (q - 2)^3
    double w = (q - 2.0);
    w *= w*w;
    // normalization: (4 pi R^5)/105, constant below is 4pi/105
    const double normalization = kernel_sqRadius_*kernel_sqRadius_*kernel_radius_*0.1196797201367540281;
    return w/normalization;
}

}

inline double FixContinuumWeighted::compute_line_sphere_intersection(const double *const xij, const double *const xjk)
{
    const double dk = vectorMag3D(xjk); // distance r_j -> r_k
    // compute line sphere intersection (sphere with center xi and radius kernel_radius)
    const double nkj[3] = {-xjk[0]/dk, -xjk[1]/dk, -xjk[2]/dk};
    const double eDotx = -vectorDot3D(nkj, xij); // - coming from xij = x_i - x_j
    const double sqrtVal = eDotx*eDotx - vectorMag3DSquared(xij) + kernel_sqRadius_;
    if (sqrtVal > 1e-6)
    {
        // two intersection points of line spanned by r_j and r_k found
        const double dm = - eDotx - sqrt(sqrtVal); // signed distance from r_j to first intersection
        const double dp = - eDotx + sqrt(sqrtVal); // signed distance from r_j to second intersection
        if (dp > 0.0 && dm < dk) { // check if intersection interval intersects with r_j -> r_k
            const double lower = fmax(0.0, dm);
            const double upper = fmin(dk, dp);
            const double integralVal = integrate_phi(xij, nkj, lower, upper); // integral of the weighting function with radius kernel_sqRadius from r_j to r_k
                                                                              // note r=0 in the argument of the weighting function only works because of the heavyside function
            return integralVal;
        }
    }
    return 0.0;
}

inline double FixContinuumWeighted::integrate_phi(const double *const xij, const double *const nkj, const double a, const double b)
{
    const double length = b - a;
    double integral = 0.0;
    for (int i = 0; i <= 10; i++)
    {
        const double t = -(double)i/10.0*length;
        double r[3];
        vectorAddMultiple3D(xij, t, nkj, r);
        const double dist = vectorMag3D(r);
        integral += get_phi(dist);
    }
    return integral / 11.0 * length;
}

inline double FixContinuumWeighted::get_phi(const double r)
{
    if (kernel_type == TOP_HAT)
        return weightingFunction<TOP_HAT>(r);
    else if (kernel_type == GAUSSIAN)
        return weightingFunction<GAUSSIAN>(r);
    else if (kernel_type == WENDLAND)
        return weightingFunction<WENDLAND>(r);
    return -1.0;
}

/*
 * This function returns the norm of the gradient of phi divided by the distance r.
 * This allows to get the real gradient simply by multiplying the return value of
 * this function by the distance between two particles.
 */
inline double FixContinuumWeighted::get_grad_phi(const double r)
{
    if (kernel_type == TOP_HAT)
        return gradWeightingFunction<TOP_HAT>(r);
    else if (kernel_type == GAUSSIAN)
        return gradWeightingFunction<GAUSSIAN>(r);
    else if (kernel_type == WENDLAND)
        return gradWeightingFunction<WENDLAND>(r);
    return -1.0;
}

void FixContinuumWeighted::post_integrate()
{
    const double *const *const x = atom->x;
    const double *const *const v = atom->v;
    const double *const mass = atom->rmass;
    const int *const mask = atom->mask;
    const int nlocal = atom->nlocal;

    double **cont_vars = fix_cont_vars_->array_atom;
    double **strain = fix_strain_->array_atom;

    NeighList *list = pairgran_->list;
    const int inum = list->inum;
    const int *const ilist = list->ilist;
    const int *const numneigh = list->numneigh;
    const int *const *const firstneigh = list->firstneigh;

    for (int ii = 0; ii < nlocal; ii++)
        vectorZeroizeN(cont_vars[ii],7);
    if (compute_strain)
    {
        for (int i = 0; i < nlocal + atom->nghost; i++)
            vectorZeroizeN(strain[i], 9);
    }

    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;
        const double xi = x[i][0];
        const double yi = x[i][1];
        const double zi = x[i][2];
        const double vxi = v[i][0];
        const double vyi = v[i][1];
        const double vzi = v[i][2];
        const double mi = mass[i];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj <= jnum; jj++) {
            const int j = (jj == jnum ? i : jlist[jj]); // last j is i itself
            if (!(mask[j] & groupbit)) continue;

            // compute particle distance
            const double xij[3] = {xi - x[j][0], yi - x[j][1], zi - x[j][2]};
            const double sqDist = vectorMag3DSquared(xij);

            // assuming top hat kernel
            if (sqDist < kernel_sqRadius_) {
                const double mj = mass[j];
                const double dist = sqrt(sqDist);
                const double phi = get_phi(dist);
                const double rhoj = mj*phi;
                const double grad_phi = get_grad_phi(dist);
                const double vxj = v[j][0];
                const double vyj = v[j][1];
                const double vzj = v[j][2];
                cont_vars[i][0] += rhoj*vxj;
                cont_vars[i][1] += rhoj*vyj;
                cont_vars[i][2] += rhoj*vzj;
                cont_vars[i][3] += rhoj;
                if (compute_strain)
                {
                    const double gradrhoj = mj*grad_phi;
                    cont_vars[i][4] += gradrhoj*xij[0];
                    cont_vars[i][5] += gradrhoj*xij[1];
                    cont_vars[i][6] += gradrhoj*xij[2];
                    strain[i][0] += gradrhoj*xij[0]*vxj;
                    strain[i][1] += gradrhoj*xij[1]*vxj;
                    strain[i][2] += gradrhoj*xij[2]*vxj;
                    strain[i][3] += gradrhoj*xij[0]*vyj;
                    strain[i][4] += gradrhoj*xij[1]*vyj;
                    strain[i][5] += gradrhoj*xij[2]*vyj;
                    strain[i][6] += gradrhoj*xij[0]*vzj;
                    strain[i][7] += gradrhoj*xij[1]*vzj;
                    strain[i][8] += gradrhoj*xij[2]*vzj;
                }
                if (j < nlocal && j != i)
                {
                    const double rhoi = mi*phi;
                    cont_vars[j][0] += rhoi*vxi;
                    cont_vars[j][1] += rhoi*vyi;
                    cont_vars[j][2] += rhoi*vzi;
                    cont_vars[j][3] += rhoi;
                    if (compute_strain)
                    {
                        const double gradrhoi = mi*grad_phi;
                        cont_vars[j][4] -= gradrhoi*xij[0];
                        cont_vars[j][5] -= gradrhoi*xij[1];
                        cont_vars[j][6] -= gradrhoi*xij[2];
                        strain[j][0] -= gradrhoi*xij[0]*vxi;
                        strain[j][1] -= gradrhoi*xij[1]*vxi;
                        strain[j][2] -= gradrhoi*xij[2]*vxi;
                        strain[j][3] -= gradrhoi*xij[0]*vyi;
                        strain[j][4] -= gradrhoi*xij[1]*vyi;
                        strain[j][5] -= gradrhoi*xij[2]*vyi;
                        strain[j][6] -= gradrhoi*xij[0]*vzi;
                        strain[j][7] -= gradrhoi*xij[1]*vzi;
                        strain[j][8] -= gradrhoi*xij[2]*vzi;
                    }
                }
            }
        }

        // wall contribution for strain
        if (compute_strain) // TODO maybe also for stress vel avg?
        {
            std::vector<FixContactPropertyAtom *>::iterator it;
            for (it = fix_wall_contact_forces_vector_.begin(); it < fix_wall_contact_forces_vector_.end(); it++)
            {
                int n_contacts = (*it)->n_partner(i);
                for (int j = 0; j < n_contacts; j++)
                {
                    const double *const force_pos_ij = (*it)->contacthistory(i, j);
                    const double *const pos = &(force_pos_ij[3]);
                    const double *const vel = &(force_pos_ij[6]);
                    const double dist = 2.0*vectorMag3D(pos); // pos points to the contact points thus distance to center is twice that
                    if (dist < kernel_radius_) {
                        const double rhoi = mi*get_phi(dist);
                        const double gradrhoi = 2.0*mi*get_grad_phi(dist); // 2.0 again for same reason
                        fflush(stdout);
                        cont_vars[i][0] += rhoi*vel[0]; // vel of wall
                        cont_vars[i][1] += rhoi*vel[1];
                        cont_vars[i][2] += rhoi*vel[2];
                        cont_vars[i][3] += rhoi;
                        if (compute_strain)
                        {
                            cont_vars[i][4] += gradrhoi*pos[0];
                            cont_vars[i][5] += gradrhoi*pos[1];
                            cont_vars[i][6] += gradrhoi*pos[2];
                            strain[i][0] += gradrhoi*pos[0]*vel[0];
                            strain[i][1] += gradrhoi*pos[1]*vel[0];
                            strain[i][2] += gradrhoi*pos[2]*vel[0];
                            strain[i][3] += gradrhoi*pos[0]*vel[1];
                            strain[i][4] += gradrhoi*pos[1]*vel[1];
                            strain[i][5] += gradrhoi*pos[2]*vel[1];
                            strain[i][6] += gradrhoi*pos[0]*vel[2];
                            strain[i][7] += gradrhoi*pos[1]*vel[2];
                            strain[i][8] += gradrhoi*pos[2]*vel[2];
                        }
                    }
                }
            }
        }
    }

    // communicate results to ghosts
    fix_cont_vars_->do_forward_comm();

    // compute actual average velocity
    for (int i = 0; i < nlocal + atom->nghost; i++)
        vectorScalarDiv3D(cont_vars[i], cont_vars[i][3]);

    double **stress = fix_stress_->array_atom;

    const double dt = update->dt;

    const int *const tag = atom->tag;

    if (compute_stress)
    {
        for (int i = 0; i < nlocal + atom->nghost; i++)
            vectorZeroizeN(stress[i], 9);
    }
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;
        const double xi[3] = {x[i][0], x[i][1], x[i][2]};
        const double vi[3] = {v[i][0], v[i][1], v[i][2]};
        const double mi = mass[i];
        const int *jlist = firstneigh[i];
        const int jnum = numneigh[i];
        const int tag_i = tag[i];

        std::set<int> i_contacts;

        // standard neighbor loop
        for (int jj = -1; jj < jnum; jj++) {
            const int j = (jj == -1 ? i : jlist[jj]); // first j is i itself
            if (!(mask[j] & groupbit)) continue;

            // compute particle distance
            const double xj[3] = {x[j][0], x[j][1], x[j][2]};
            const double xij[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
            const double vj[3] = {v[j][0], v[j][1], v[j][2]};
            const double sqDist = vectorMag3DSquared(xij);

            // all kernels are scaled so that they are either 0 or can be cut off (Gaussian) at kernel_radius_
            if (compute_stress && sqDist < kernel_sqRadius_) {
                const double vpxj = vj[0]-cont_vars[i][0];
                const double vpyj = vj[1]-cont_vars[i][1];
                const double vpzj = vj[2]-cont_vars[i][2];
                const double mj = mass[j];
                const double phi = get_phi(sqDist);
                const double phi_mj = mj*phi;
                // j -> i
                stress[i][0] -= vpxj * vpxj * phi_mj;
                stress[i][1] -= vpxj * vpyj * phi_mj;
                stress[i][2] -= vpxj * vpzj * phi_mj;
                stress[i][3] -= vpyj * vpxj * phi_mj;
                stress[i][4] -= vpyj * vpyj * phi_mj;
                stress[i][5] -= vpyj * vpzj * phi_mj;
                stress[i][6] -= vpzj * vpxj * phi_mj;
                stress[i][7] -= vpzj * vpyj * phi_mj;
                stress[i][8] -= vpzj * vpzj * phi_mj;
                // i -> j
                if (j != i && j < nlocal) { // don't add the contribution i->i twice
                    const double vpxi = vi[0] - cont_vars[j][0];
                    const double vpyi = vi[1] - cont_vars[j][1];
                    const double vpzi = vi[2] - cont_vars[j][2];
                    const double phi_mi = mi*phi;
                    stress[j][0] -= vpxi * vpxi * phi_mi;
                    stress[j][1] -= vpxi * vpyi * phi_mi;
                    stress[j][2] -= vpxi * vpzi * phi_mi;
                    stress[j][3] -= vpyi * vpxi * phi_mi;
                    stress[j][4] -= vpyi * vpyi * phi_mi;
                    stress[j][5] -= vpyi * vpzi * phi_mi;
                    stress[j][6] -= vpzi * vpxi * phi_mi;
                    stress[j][7] -= vpzi * vpyi * phi_mi;
                    stress[j][8] -= vpzi * vpzi * phi_mi;
                }
            }

            // compute wall contribution
            if (compute_stress)
            {
                std::vector<FixContactPropertyAtom *>::iterator it;
                for (it = fix_wall_contact_forces_vector_.begin(); it < fix_wall_contact_forces_vector_.end(); it++)
                {
                    int n_contacts = (*it)->n_partner(j);
                    for (int k = 0; k < n_contacts; k++)
                    {
                        const double *const force_pos_jk = (*it)->contacthistory(j, k);
                        const double integralVal_ijk = compute_line_sphere_intersection(xij, &(force_pos_jk[3]));
                        stress[i][0] -= force_pos_jk[0]*force_pos_jk[3]*integralVal_ijk;
                        stress[i][1] -= force_pos_jk[0]*force_pos_jk[4]*integralVal_ijk;
                        stress[i][2] -= force_pos_jk[0]*force_pos_jk[5]*integralVal_ijk;
                        stress[i][3] -= force_pos_jk[1]*force_pos_jk[3]*integralVal_ijk;
                        stress[i][4] -= force_pos_jk[1]*force_pos_jk[4]*integralVal_ijk;
                        stress[i][5] -= force_pos_jk[1]*force_pos_jk[5]*integralVal_ijk;
                        stress[i][6] -= force_pos_jk[2]*force_pos_jk[3]*integralVal_ijk;
                        stress[i][7] -= force_pos_jk[2]*force_pos_jk[4]*integralVal_ijk;
                        stress[i][8] -= force_pos_jk[2]*force_pos_jk[5]*integralVal_ijk;
                    }
                    if (jj != -1)
                    {
                        int n_contacts = (*it)->n_partner(i);
                        for (int k = 0; k < n_contacts; k++)
                        {
                            const double *const force_pos_ik = (*it)->contacthistory(i, k);
                            const double xji[3] = {-xij[0], -xij[1], -xij[2]};
                            const double integralVal_jik = compute_line_sphere_intersection(xji, &(force_pos_ik[3]));
                            stress[j][0] -= force_pos_ik[0]*force_pos_ik[3]*integralVal_jik;
                            stress[j][1] -= force_pos_ik[0]*force_pos_ik[4]*integralVal_jik;
                            stress[j][2] -= force_pos_ik[0]*force_pos_ik[5]*integralVal_jik;
                            stress[j][3] -= force_pos_ik[1]*force_pos_ik[3]*integralVal_jik;
                            stress[j][4] -= force_pos_ik[1]*force_pos_ik[4]*integralVal_jik;
                            stress[j][5] -= force_pos_ik[1]*force_pos_ik[5]*integralVal_jik;
                            stress[j][6] -= force_pos_ik[2]*force_pos_ik[3]*integralVal_jik;
                            stress[j][7] -= force_pos_ik[2]*force_pos_ik[4]*integralVal_jik;
                            stress[j][8] -= force_pos_ik[2]*force_pos_ik[5]*integralVal_jik;
                        }
                    }
                }
            }

        } // end loop j neighbors

        // contactproperty atom loop
        for (int jj = 0; jj < fix_contact_forces_->get_npartners(i); jj++)
        {
            const double *const force_pos_ij = fix_contact_forces_->contacthistory(i, jj);
            const int j = (int)force_pos_ij[3];
            i_contacts.insert(tag[j]);
        }

        if (compute_stress)
        {
            // contactproperty atom loop
            for (int jj = 0; jj < fix_contact_forces_->get_npartners(i); jj++)
            {
                const double *const force_pos_ij = fix_contact_forces_->contacthistory(i, jj);
                const int j = (int)force_pos_ij[3];
                const double xji[3] = {x[j][0] - xi[0], x[j][1] - xi[1], x[j][2] - xi[2]};
                const double sqDistji = vectorMag3DSquared(xji);

                if (!compute_stress && sqDistji > kernel_sqRadius_)
                    continue;

                std::set<int> skip_contacts;
                for (int ll = 0; ll < fix_contact_forces_->get_npartners(j); ll++)
                {
                    const double *const force_pos_jl = fix_contact_forces_->contacthistory(j, ll);
                    const int l = (int)force_pos_jl[3];
                    const int tag_l = tag[l];
                    if (tag_l < tag_i && i_contacts.find(tag_l) != i_contacts.end())
                        skip_contacts.insert(tag_l);
                }

                for (int kk = 0; kk < fix_contact_forces_->get_npartners(i); kk++)
                {
                    const double *const force_pos_ik = fix_contact_forces_->contacthistory(i, kk);
                    const int k = (int)force_pos_ik[3];
                    const int tag_k = tag[k];
                    const double xik[3] = {xi[0] - x[k][0], xi[1] - x[k][1], xi[2] - x[k][2]};

                    if (skip_contacts.find(tag_k) != skip_contacts.end())
                        continue;

                    // stress contribution
                    // check if forces are actually present between those two particles, if not go to next contact
                    const double absForces = vectorMag3DSquared(force_pos_ik);
                    if (absForces < 1e-6)
                        continue;

                    const double integralVal_ijk = 0.5*compute_line_sphere_intersection(xji, xik);
                    if (integralVal_ijk > 1e-6)
                    {
                        stress[j][0] -= force_pos_ik[0]*xik[0]*integralVal_ijk; // - 1/2 f_ij (x) r_ij * integral
                        stress[j][1] -= force_pos_ik[0]*xik[1]*integralVal_ijk;
                        stress[j][2] -= force_pos_ik[0]*xik[2]*integralVal_ijk;
                        stress[j][3] -= force_pos_ik[1]*xik[0]*integralVal_ijk;
                        stress[j][4] -= force_pos_ik[1]*xik[1]*integralVal_ijk;
                        stress[j][5] -= force_pos_ik[1]*xik[2]*integralVal_ijk;
                        stress[j][6] -= force_pos_ik[2]*xik[0]*integralVal_ijk;
                        stress[j][7] -= force_pos_ik[2]*xik[1]*integralVal_ijk;
                        stress[j][8] -= force_pos_ik[2]*xik[2]*integralVal_ijk;
                    }
                }
            }
        }

    }

    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        //if(tag[i] == 87) printf("stress[i][0]: %e %d\n", stress[i][0], update->ntimestep);
        if (compute_strain)
        {
            const double rho = cont_vars[i][3];
            const double inv2Rho = 0.5/rho;
            const double pDivRho[3] = {cont_vars[i][0], cont_vars[i][1], cont_vars[i][2]}; // momentum density was already divided by rho before
            const double gradRho[3] = {cont_vars[i][4], cont_vars[i][5], cont_vars[i][6]};
            strain[i][0] = dt*(strain[i][0] - pDivRho[0]*gradRho[0])*inv2Rho;
            strain[i][1] = dt*(strain[i][1] - pDivRho[0]*gradRho[1])*inv2Rho;
            strain[i][2] = dt*(strain[i][2] - pDivRho[0]*gradRho[2])*inv2Rho;
            strain[i][3] = dt*(strain[i][3] - pDivRho[1]*gradRho[0])*inv2Rho;
            strain[i][4] = dt*(strain[i][4] - pDivRho[1]*gradRho[1])*inv2Rho;
            strain[i][5] = dt*(strain[i][5] - pDivRho[1]*gradRho[2])*inv2Rho;
            strain[i][6] = dt*(strain[i][6] - pDivRho[2]*gradRho[0])*inv2Rho;
            strain[i][7] = dt*(strain[i][7] - pDivRho[2]*gradRho[1])*inv2Rho;
            strain[i][8] = dt*(strain[i][8] - pDivRho[2]*gradRho[2])*inv2Rho;

        }
    }

    // reverse communicate results, i.e. add contributions from neib procs together
    // here neib procs only compute the contributions occuring from contacts with both
    // grains on another proc.
    if (compute_stress)
        fix_stress_->do_reverse_comm();
}
