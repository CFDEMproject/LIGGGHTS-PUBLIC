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

#ifndef LMP_MULTISPHERE_I_H
#define LMP_MULTISPHERE_I_H

/* ---------------------------------------------------------------------- */

inline double Multisphere::max_r_bound()
{
    double max_r_bound = 0.;

    for(int i = 0; i < nbody_; i++)
        max_r_bound = std::max(max_r_bound,r_bound_(i));

    MPI_Max_Scalar(max_r_bound,world);

    return max_r_bound;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::copy_body(int from_local, int to_local)
{
    int tag_from =  id_(from_local);

    customValues_.copyElement(from_local, to_local);

    mapArray_[tag_from] = to_local;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::remove_body(int ilocal)
{
    
    mapArray_[id_(ilocal)] = -1;
    if(nbody_ > 1) mapArray_[id_(nbody_-1)] = ilocal;

    /*if(ilocal < nbody_-1)
        copy_body(nbody_-1,ilocal);*/
    customValues_.deleteElement(ilocal);

    nbody_--;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::calc_nbody_all()
{
   MPI_Sum_Scalar(nbody_,nbody_all_,world);
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::reset_forces(bool extflag)
{
    
    fcm_.setAll(nbody_,0.);
    torquecm_.setAll(nbody_,0.);
    if(extflag) dragforce_cm_.setAll(nbody_,0.);
    if(extflag) hdtorque_cm_.setAll(nbody_,0.);
}

/* ---------------------------------------------------------------------- */

inline int Multisphere::calc_n_steps(int, int body, double *p_ref, double *normalvec, double *v_normal)
{
    double pos_rel[3],dt,dist_normal;
    int ibody,timestep,n_steps;

    if(body < 0)
        return -1;

    ibody = map(body);

    timestep = update->ntimestep;
    dt = update->dt;

    if(ibody < 0)
        error->one(FLERR,"Illegal situation in FixMultisphere::calc_n_steps");

    if(start_step_(ibody) >= 0)
        return (start_step_(ibody) - timestep);

    vectorSubtract3D(p_ref,xcm_(ibody),pos_rel);
    dist_normal = vectorDot3D(pos_rel,normalvec);
    n_steps = static_cast<int>(dist_normal/(vectorMag3D(v_normal)*dt));
    start_step_(ibody) = n_steps + timestep;
    v_integrate_.set(ibody,v_normal);

    return n_steps;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::recalc_n_steps(double dt_ratio)
{
  
  for(int ibody = 0; ibody < nbody_; ibody++)
  {
        if(start_step_(ibody) > update->ntimestep)
        {
            
            start_step_(ibody) = static_cast<int>(update->ntimestep) + static_cast<int>(dt_ratio*(static_cast<double>(start_step_(ibody)) - static_cast<double>(update->ntimestep)));
            
        }
  }
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::release(int iatom,int body,double *v_toInsert,double *omega_toInsert)
{
    int ibody;

    if(body < 0)
        return;

    ibody = map(body);
    if(ibody < 0)
        return;

    bigint step = update->ntimestep;

    if(start_step_(ibody) >= 0 && step >= start_step_(ibody))
    {
        
        // set v and omega
        vcm_.set(ibody,v_toInsert);
        omega_.set(ibody,omega_toInsert);
        start_step_.set(ibody,-1);
    }
}

#endif
