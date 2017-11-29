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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
Thanks to Chris Stoltz (P&G) for providing
a Fortran version of the MC integrator
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_template_multisphere.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "vector_liggghts.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "fix_multisphere.h"
#include "memory.h"
#include "error.h"
#include "random_mars.h"
#include "random_park.h"
#include "fix_rigid.h"
#include "particleToInsert_multisphere.h"
#include "input_multisphere.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;
using namespace FixConst;

#define EPSILON 1.0e-7
#define TOLERANCE 1.0e-6

/* ---------------------------------------------------------------------- */

FixTemplateMultisphere::FixTemplateMultisphere(LAMMPS *lmp, int narg, char **arg) :
  FixTemplateMultiplespheres(lmp, narg, arg),
  mass_set_(false),
  moi_set_(false),
  use_density_(-1)
{
    delete pti;
    pti = new ParticleToInsertMultisphere(lmp,nspheres);

    memory->create(displace_,nspheres,3,"FixTemplateMultiplespheres:moi_");

    volumeweight_ = new double[nspheres];

    type_ = 0;

    mass_expect = 0;
    vectorZeroize3D(inertia_);
    vectorZeroize3D(ex_space_);
    vectorZeroize3D(ey_space_);
    vectorZeroize3D(ez_space_);
    vectorZeroize3D(moi_[0]);
    vectorZeroize3D(moi_[1]);
    vectorZeroize3D(moi_[2]);

    fflag_[0] = fflag_[1] = fflag_[2] = true;
    tflag_[0] = tflag_[1] = tflag_[2] = true;

    // parse args
    // should have the possibility to parse inertia tensor here
    // TODO: should parse cross section area matrix here

    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg],"type") == 0) {
            if (iarg+2 > narg)
                error->fix_error(FLERR,this,"not enough arguments for 'type'");
            iarg++;
            type_ = atoi(arg[iarg++]);
            hasargs = true;
        } else if(strcmp(arg[iarg],"mass") == 0) {
            if (iarg+2 > narg)
                error->fix_error(FLERR,this,"not enough arguments for 'mass'");
            iarg++;
            mass_expect = atof(arg[iarg++]);
            mass_set_ = true;
            hasargs = true;
        } else if(strcmp(arg[iarg],"use_volume") == 0) {
            iarg++;
            use_density_ = 0;
            hasargs = true;
        } else if(strcmp(arg[iarg],"use_density") == 0) {
            iarg++;
            use_density_ = 1;
            hasargs = true;
        } else if(strcmp(arg[iarg],"fflag") == 0) {
            if (iarg+4 > narg)
                error->fix_error(FLERR,this,"not enough arguments for 'fflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                fflag_[0] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                fflag_[0] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'fflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                fflag_[1] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                fflag_[1] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'fflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                fflag_[2] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                fflag_[2] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'fflag'");
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"tflag") == 0) {
            if (iarg+4 > narg)
                error->fix_error(FLERR,this,"not enough arguments for 'tflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                tflag_[0] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                tflag_[0] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'tflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                tflag_[1] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                tflag_[1] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'tflag'");
            iarg++;
            if(strcmp(arg[iarg],"on") == 0)
                tflag_[2] = true;
            else if(strcmp(arg[iarg],"off") == 0)
                tflag_[2] = false;
            else
                error->fix_error(FLERR,this,"expecting 'on or 'off' after 'tflag'");
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"inertia_tensor") == 0) {
            if (iarg+7 > narg)
                error->fix_error(FLERR,this,"not enough arguments for 'inertia_tensor'");
            iarg++;
            moi_[0][0] = atof(arg[iarg++]);
            moi_[0][1] = atof(arg[iarg++]);
            moi_[0][2] = atof(arg[iarg++]);
            moi_[1][1] = atof(arg[iarg++]);
            moi_[1][2] = atof(arg[iarg++]);
            moi_[2][2] = atof(arg[iarg++]);

            moi_[1][0] = moi_[0][1];
            moi_[2][0] = moi_[0][2];
            moi_[2][1] = moi_[1][2];
            moi_set_ = true;
            hasargs = true;
        } else if(strcmp(style,"particletemplate/multisphere") == 0)
            error->fix_error(FLERR,this,"unknown keyword");
    }

    // check if type has been defined
    if(!lmp->wb && type_ < 1)
        error->fix_error(FLERR,this,"have to provide a type >=1");

    if((mass_set_ && !moi_set_) || (!mass_set_ && moi_set_))
        error->fix_error(FLERR,this,"have to define either both 'mass' and 'inertia_tensor' or none of the two");

    if(mass_set_ && -1 == use_density_)
        error->fix_error(FLERR,this,"have to use either 'use_volume' or 'use_density' in case 'mass' and 'moi' are specified");
}

/* ---------------------------------------------------------------------- */

void FixTemplateMultisphere::post_create()
{
  // bounding sphere, center of mass calculations in mother class
  // calculate inertia_ tensor and its eigensystem

  // use density from input script and volume from MC; mass calculated
  if(!mass_set_)
  {
      // bounding sphere and center of mass

      FixTemplateMultiplespheres::post_create();
      calc_inertia();
      calc_eigensystem();
      calc_displace_xcm_x_body();
  }
  // use mass specified and volume_mc; density calculated
  else if(0 == use_density_)
  {
      // store mass read in from input script
      double mass = mass_expect;

      // bounding sphere and center of mass
      // calculates mass_expect, volume_expect, r_equiv
      FixTemplateMultiplespheres::post_create();

      // re-apply mass read in from input script
      // re-calc r_equiv without using density
      mass_expect = mass;
      r_equiv = pow(3.*volume_expect/(4.*M_PI),1./3.);
      pdf_density->set_params<RANDOM_CONSTANT>(mass_expect/volume_expect);
      
      // calc rest of stuff
      calc_eigensystem();
      calc_displace_xcm_x_body();
  }
  // use density specified in input script  and mass specified in input script ; volume calculated
  else
  {
      calc_bounding_sphere();

      check_overlap();

      calc_eigensystem();
      calc_displace_xcm_x_body();

      // calc expectancy values
      volume_expect = mass_expect / expectancy(pdf_density);
      r_equiv = pow(6.*mass_expect/(8.*expectancy(pdf_density)*M_PI),1./3.);
      
  }

  calc_volumeweight();

  print_info();
}

/* ---------------------------------------------------------------------- */

void FixTemplateMultisphere::init()
{
    // check consistency of type as defined in each command of this style
    // types must be from 1..ntemplates
    // is fulfilled if min and max is reached

    FixTemplateMultisphere *ftms_i,*ftms_j;

    int type_i;
    int type_min = 10000, type_max = 0;
    int n_fixes = modify->n_fixes_style_strict(style);

    for(int i = 0; i < n_fixes; i++)
    {
        ftms_i = static_cast<FixTemplateMultisphere*>(modify->find_fix_style_strict(style,i));
        type_i = ftms_i->type();

        if(type_i < type_min) type_min = type_i;
        if(type_i > type_max) type_max = type_i;

        // check if any other uses the same type
        for(int j = i+1; j < n_fixes; j++)
        {
            ftms_j = static_cast<FixTemplateMultisphere*>(modify->find_fix_style_strict(style,j));
            if(!lmp->wb && ftms_j != ftms_i && type_i == ftms_j->type())
                error->fix_error(FLERR,this,"multisphere template types have to be unique");
        }
    }

    // types are consecutive if no double usage and min/max is met
    if(!lmp->wb &&  (type_min != 1 || type_max != n_fixes))
        error->fix_error(FLERR,this,"multisphere template types have to be consecutive starting from 1");

    FixMultisphere *fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere", 0));
    if(fix_multisphere && (fix_multisphere->groupbit & groupbit) == 0)
        error->fix_error(FLERR,this,"Fix particletemplate/multisphere command and fix multisphere command are not compatible, must be member of the same group");
}

/* ----------------------------------------------------------------------
   calc volume weight for each sphere
   necessary for volume fraction on CFD side
------------------------------------------------------------------------- */

void FixTemplateMultisphere::calc_volumeweight()
{
    double x_try[3],distSqr,n_hits;

    bool *hits_j = new bool[nspheres];

    vectorZeroizeN(volumeweight_,nspheres);

    // volumeweight_ = 0.5 + 0.5 * hits only in j / hits in j

    int hits_tot = 0;
    for(int i = 0; i < ntry; i++)
    {
        generate_xtry(x_try);
        n_hits = 0;
        vectorInitializeN(hits_j,nspheres,false);

        for(int j = 0; j < nspheres; j++)
        {
            distSqr = dist_sqr(j,x_try);
            if(distSqr < r_sphere[j]*r_sphere[j])
            {
                hits_j[j] = true;
                n_hits++;
            }
        }

        for(int j = 0; j < nspheres; j++)
        {
            if(hits_j[j])
                volumeweight_[j] += 1./static_cast<double>(n_hits);
        }

        if (n_hits > 0)
            hits_tot++;
    }

    // calculate volume weight
    if (hits_tot > 0) {
        const double hits_tot_inv = 1./static_cast<double>(hits_tot);
        for(int j = 0; j < nspheres; j++)
            volumeweight_[j] *= hits_tot_inv;
    }

    delete []hits_j;
}

/* ----------------------------------------------------------------------
   calc inertia tensor
------------------------------------------------------------------------- */

void FixTemplateMultisphere::calc_inertia()
{
  double x_try[3],xcm[3],distSqr;
  bool alreadyChecked;

  for(int i = 0; i < 3; i++)
    vectorZeroize3D(moi_[i]);

  vectorZeroize3D(xcm);

  for(int i = 0; i < ntry; i++)
  {
      generate_xtry(x_try);

      alreadyChecked = false;
      for(int j = 0; j < nspheres; j++)
      {
          if (alreadyChecked) break;
          distSqr = dist_sqr(j,x_try);
          if(distSqr < r_sphere[j]*r_sphere[j])
          {
              moi_[0][0] +=  (x_try[1]-xcm[1])*(x_try[1]-xcm[1]) + (x_try[2]-xcm[2])*(x_try[2]-xcm[2]);
              moi_[0][1] -=  (x_try[0]-xcm[0])*(x_try[1]-xcm[1]);
              moi_[0][2] -=  (x_try[0]-xcm[0])*(x_try[2]-xcm[2]);
              moi_[1][0] -=  (x_try[1]-xcm[1])*(x_try[0]-xcm[0]);
              moi_[1][1] +=  (x_try[0]-xcm[0])*(x_try[0]-xcm[0]) + (x_try[2]-xcm[2])*(x_try[2]-xcm[2]);
              moi_[1][2] -=  (x_try[1]-xcm[1])*(x_try[2]-xcm[2]);
              moi_[2][0] -=  (x_try[2]-xcm[2])*(x_try[0]-xcm[0]);
              moi_[2][1] -=  (x_try[2]-xcm[2])*(x_try[1]-xcm[1]);
              moi_[2][2] +=  (x_try[0]-xcm[0])*(x_try[0]-xcm[0]) + (x_try[1]-xcm[1])*(x_try[1]-xcm[1]);
              alreadyChecked = true;
          }
      }
  }
  for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
          moi_[i][j] *= expectancy(pdf_density)/static_cast<double>(ntry)*(x_max[0]-x_min[0])*(x_max[1]-x_min[1])*(x_max[2]-x_min[2]);

  // check if tensor is symmetric enough(numerical errors)
  // not sure if this check is sufficient
  
  if(fabs(moi_[0][1]/moi_[1][0]-1.) > TOLERANCE) error->all(FLERR,"Fix particletemplate/multisphere:Error when calculating inertia_ tensor : Not enough accuracy. Boost ntry.");
  if(fabs(moi_[0][2]/moi_[2][0]-1.) > TOLERANCE) error->all(FLERR,"Fix particletemplate/multisphere:Error when calculating inertia_ tensor : Not enough accuracy. Boost ntry.");
  if(fabs(moi_[2][1]/moi_[1][2]-1.) > TOLERANCE) error->all(FLERR,"Fix particletemplate/multisphere:Error when calculating inertia_ tensor : Not enough accuracy. Boost ntry.");

  // make the tensor symmetric
  
  moi_[0][1] = (moi_[0][1]+moi_[1][0])/2.;
  moi_[1][0] = moi_[0][1];
  moi_[0][2] = (moi_[0][2]+moi_[2][0])/2.;
  moi_[2][0] = moi_[0][2];
  moi_[2][1] = (moi_[2][1]+moi_[2][0])/2.;
  moi_[1][2] = moi_[2][1];
}

/* ----------------------------------------------------------------------
   calc eigensystem of inertia_ tensor
------------------------------------------------------------------------- */

void FixTemplateMultisphere::calc_eigensystem()
{
  // following operations as in fix_rigid.cpp::init()

  //-----------------------------------------
  // calculate eigenvalues and eigenvectors
  //-----------------------------------------

  double evectors[3][3];
  
  int ierror = MathExtra::jacobi(moi_,static_cast<double*>(inertia_),evectors);
  if (ierror) error->fix_error(FLERR,this,"Insufficient Jacobi rotations for rigid body");

  ex_space_[0] = evectors[0][0];
  ex_space_[1] = evectors[1][0];
  ex_space_[2] = evectors[2][0];
  ey_space_[0] = evectors[0][1];
  ey_space_[1] = evectors[1][1];
  ey_space_[2] = evectors[2][1];
  ez_space_[0] = evectors[0][2];
  ez_space_[1] = evectors[1][2];
  ez_space_[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0
  
  double max;
  max = MAX(inertia_[0],inertia_[1]);
  max = MAX(max,inertia_[2]);

  if (inertia_[0] < EPSILON*max) inertia_[0] = 0.0;
  if (inertia_[1] < EPSILON*max) inertia_[1] = 0.0;
  if (inertia_[2] < EPSILON*max) inertia_[2] = 0.0;

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd evector if needed

  double ez[3];
  vectorCross3D(ex_space_,ey_space_,ez);
  double dot = vectorDot3D(ez,ez_space_);
  if (dot < 0.) vectorScalarMult3D(ez_space_,-1.);

}

/* ----------------------------------------------------------------------
   calc displace and xcm_to_xb_body_
------------------------------------------------------------------------- */

void FixTemplateMultisphere::calc_displace_xcm_x_body()
{
  // calculate displace_
  for(int i = 0; i < nspheres; i++)
  {
      
      // can do it like this because ex, ey, ez are orthogonal, using xcm = 0/0/0
      MathExtraLiggghts::cartesian_coosys_to_local_orthogonal(displace_[i],x_sphere[i], ex_space_, ey_space_, ez_space_,error);
      
  }

  // transform in body coodinates, could do this with
  // solve Mt*xcm_to_xb_body = xcm_to_xb (where xcm_to_xb == x_bound because xcm is 0 0 0)
  MathExtraLiggghts::cartesian_coosys_to_local_orthogonal(xcm_to_xb_body_,x_bound,ex_space_,ey_space_,ez_space_,error);

}

/* ---------------------------------------------------------------------- */

void FixTemplateMultisphere::print_info()
{
  if (logfile)
  {
    fprintf(logfile,"Finished calculating properties of template\n");
    fprintf(logfile,"   mass = %e, radius of bounding sphere = %e, radius of equivalent sphere = %e\n",mass_expect,r_bound,r_equiv);
    fprintf(logfile,"   center of mass = %e, %e, %e\n",0.,0.,0.);
    fprintf(logfile,"   center of bounding sphere in body coords = %e, %e, %e\n",xcm_to_xb_body_[0],xcm_to_xb_body_[1],xcm_to_xb_body_[2]);
    fprintf(logfile,"   Principal moments of inertia_: %e, %e, %e\n",inertia_[0],inertia_[1],inertia_[2]);
    fprintf(logfile,"     Eigenvector: %e, %e, %e\n",ex_space_[0],ex_space_[1],ex_space_[2]);
    fprintf(logfile,"     Eigenvector: %e, %e, %e\n",ey_space_[0],ey_space_[1],ey_space_[2]);
    fprintf(logfile,"     Eigenvector: %e, %e, %e\n",ez_space_[0],ez_space_[1],ez_space_[2]);
    fprintf(logfile,"     Inertia tensor: %e, %e, %e\n",moi_[0][0],moi_[1][0],moi_[2][0]);
    fprintf(logfile,"     Inertia tensor: %e, %e, %e\n",moi_[1][0],moi_[1][1],moi_[2][1]);
    fprintf(logfile,"     Inertia tensor: %e, %e, %e\n",moi_[2][0],moi_[2][1],moi_[2][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixTemplateMultisphere::~FixTemplateMultisphere()
{
    memory->destroy(displace_);
    delete []volumeweight_;

    delete pti;
    if(pti_list) delete_ptilist();
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultisphere::randomize_single()
{
  
  pti->nparticles = nspheres;
  pti->density_ins = expectancy(pdf_density);;
  pti->volume_ins = volume_expect;
  pti->mass_ins = mass_expect;
  pti->r_bound_ins = r_bound;
  vectorCopy3D(x_bound,pti->x_bound_ins);
  pti->atom_type = atom_type;

  ParticleToInsertMultisphere *pti_m = static_cast<ParticleToInsertMultisphere*>(pti);

  pti_m->type_ms = type_;

  for(int j = 0; j < nspheres; j++)
  {
      pti_m->radius_ins[j] = r_sphere[j];
      pti_m->volumeweight[j] = volumeweight_[j];
      vectorCopy3D(x_sphere[j],pti_m->x_ins[j]);
      vectorCopy3D(displace_[j],pti_m->displace[j]);
  }

  vectorCopy3D(inertia_,pti_m->inertia);
  vectorCopy3D(ex_space_,pti_m->ex_space);
  vectorCopy3D(ey_space_,pti_m->ey_space);
  vectorCopy3D(ez_space_,pti_m->ez_space);
  vectorCopy3D(fflag_,pti_m->fflag);
  vectorCopy3D(tflag_,pti_m->tflag);
  vectorCopy3D(xcm_to_xb_body_,pti_m->xcm_to_xbound);

  vectorZeroize3D(pti_m->xcm_ins);
  quatIdentity4D(pti_m->quat_ins);
  vectorZeroize3D(pti_m->v_ins);
  vectorZeroize3D(pti_m->omega_ins);

  pti->groupbit = groupbit;

   //pti->random_rotate(random->uniform(),random->uniform(),random->uniform());
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultisphere::init_ptilist(int n_random_max, const bool enforce_single, FixPropertyAtom * const fix_release)
{
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");

    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsertMultisphere(lmp,nspheres);

}

/* ----------------------------------------------------------------------*/

void FixTemplateMultisphere::delete_ptilist()
{
    
    if(n_pti_max == 0) return;

    for(int i = 0; i < n_pti_max; i++)
       delete pti_list[i];

    memory->sfree(pti_list);
    pti_list = NULL;
    n_pti_max = 0;
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultisphere::randomize_ptilist(int n_random,int distribution_groupbit,int distorder)
{
    
    for(int i = 0; i < n_random; i++)
    {
          
          ParticleToInsertMultisphere *pti_m = static_cast<ParticleToInsertMultisphere*>(pti_list[i]);

          pti_m->nparticles = nspheres;
          pti_m->density_ins = expectancy(pdf_density);
          pti_m->type_ms = type_;
          pti_m->volume_ins = volume_expect;
          pti_m->mass_ins = mass_expect;
          pti_m->r_bound_ins = r_bound;
          vectorCopy3D(x_bound,pti->x_bound_ins);
          pti_m->atom_type = atom_type;

          for(int j = 0; j < nspheres; j++)
          {
              pti_m->radius_ins[j] = r_sphere[j];
              pti_m->volumeweight[j] = volumeweight_[j];
              vectorCopy3D(x_sphere[j],pti_m->x_ins[j]);
              vectorCopy3D(displace_[j],pti_m->displace[j]);
          }

          vectorCopy3D(inertia_,pti_m->inertia);
          vectorCopy3D(ex_space_,pti_m->ex_space);
          vectorCopy3D(ey_space_,pti_m->ey_space);
          vectorCopy3D(ez_space_,pti_m->ez_space);
          vectorCopy3D(fflag_,pti_m->fflag);
          vectorCopy3D(tflag_,pti_m->tflag);
          vectorCopy3D(xcm_to_xb_body_,pti_m->xcm_to_xbound);

          vectorZeroize3D(pti_m->xcm_ins);
          quatIdentity4D(pti_m->quat_ins);
          vectorZeroize3D(pti_m->v_ins);
          vectorZeroize3D(pti_m->omega_ins);

          pti_m->groupbit = groupbit | distribution_groupbit; 

          pti_list[i]->distorder = distorder;
    }
}

/* ---------------------------------------------------------------------- */

void FixTemplateMultisphere::finalize_insertion()
{
    if(modify->n_fixes_style("multisphere") != 1)
        error->fix_error(FLERR,this,"Multi-sphere particle inserted: You have to use exactly one fix multisphere");

    static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0))->add_body_finalize();
}
