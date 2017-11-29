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

#define DELTA 10000

#include "multisphere.h"
#include "domain.h"
#include "force.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "fix_heat_gran.h"
#include <cmath>
#include <algorithm>

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

Multisphere::Multisphere(LAMMPS *lmp) :
  Pointers(lmp),
  customValues_(*(new CustomValueTracker(lmp))),

  nbody_(0),
  nbody_all_(0),
  mapTagMax_(0), 
  mapArray_(0),

  id_ (*customValues_.addElementProperty< ScalarContainer<int> >("id_multisphere","comm_exchange_borders"/*ID does never change*/,"frame_invariant","restart_yes")),

  xcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("xcm","comm_exchange_borders","frame_invariant", "restart_yes")),
  vcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("vcm","comm_exchange_borders","frame_invariant", "restart_yes")),
  fcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("fcm","comm_none","frame_invariant", "restart_no")),
  torquecm_     (*customValues_.addElementProperty< VectorContainer<double,3> >("torque","comm_none","frame_invariant", "restart_no")),
  dragforce_cm_ (*customValues_.addElementProperty< VectorContainer<double,3> >("dragforce_cm","comm_none","frame_invariant", "restart_no")),
  hdtorque_cm_  (*customValues_.addElementProperty< VectorContainer<double,3> >("hdtorque_cm","comm_none","frame_invariant", "restart_no")),

  angmom_ (*customValues_.addElementProperty< VectorContainer<double,3> >("angmom","comm_exchange_borders","frame_invariant", "restart_yes")),
  omega_  (*customValues_.addElementProperty< VectorContainer<double,3> >("omega","comm_exchange_borders","frame_invariant", "restart_yes")),
  quat_   (*customValues_.addElementProperty< VectorContainer<double,4> >("quat","comm_exchange_borders","frame_invariant", "restart_yes")),

  atomtype_  (*customValues_.addElementProperty< ScalarContainer<int> >("atomtype","comm_exchange_borders","frame_invariant","restart_yes")),
  type_      (*customValues_.addElementProperty< ScalarContainer<int> >("clumptype","comm_exchange_borders","frame_invariant","restart_yes")),
  density_   (*customValues_.addElementProperty< ScalarContainer<double> >("density","comm_exchange_borders","frame_invariant","restart_yes")),
  masstotal_ (*customValues_.addElementProperty< ScalarContainer<double> >("masstotal","comm_exchange_borders","frame_invariant","restart_yes")),
  inertia_   (*customValues_.addElementProperty< VectorContainer<double,3> >("inertia","comm_exchange_borders","frame_invariant", "restart_yes")),
  ex_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ex_space","comm_exchange_borders","frame_invariant", "restart_yes")),
  ey_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ey_space","comm_exchange_borders","frame_invariant", "restart_yes")),
  ez_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ez_space","comm_exchange_borders","frame_invariant", "restart_yes")),

  nrigid_    (*customValues_.addElementProperty< ScalarContainer<int> >("nrigid","comm_exchange_borders","frame_invariant", "restart_yes")),

  imagebody_ (*customValues_.addElementProperty< ScalarContainer<tagint> >("imagebody","comm_exchange_borders","frame_invariant", "restart_yes")),
  remapflag_ (*customValues_.addElementProperty< VectorContainer<int,4> >("remapflag","comm_none","frame_invariant", "restart_no")),

  fflag_ (*customValues_.addElementProperty< VectorContainer<bool,3> >("fflag","comm_exchange_borders","frame_invariant", "restart_yes")),
  tflag_ (*customValues_.addElementProperty< VectorContainer<bool,3> >("tflag","comm_exchange_borders","frame_invariant", "restart_yes")),

  start_step_ (*customValues_.addElementProperty< ScalarContainer<int> >("start_step","comm_exchange_borders","frame_invariant", "restart_yes")),
  v_integrate_(*customValues_.addElementProperty< VectorContainer<double,3> >("v_integrate","comm_exchange_borders","frame_invariant", "restart_yes")),

  r_bound_       (*customValues_.addElementProperty< ScalarContainer<double> >("r_bound","comm_exchange_borders","frame_invariant", "restart_yes")),
  xcm_to_xbound_ (*customValues_.addElementProperty< VectorContainer<double,3> >("xcm_to_xbound","comm_exchange_borders","frame_invariant", "restart_yes")),

  temp_(*customValues_.addElementProperty< ScalarContainer<double> >("temp","comm_exchange_borders","frame_invariant","restart_yes")),
  temp_old_(*customValues_.addElementProperty< ScalarContainer<double> >("temp_old","comm_exchange_borders","frame_invariant","restart_yes"))
{

}

Multisphere::~Multisphere()
{
    delete &customValues_;

    // deallocate map memory if exists
    if(mapArray_) clear_map();
}

/* ----------------------------------------------------------------------
   add a new body
------------------------------------------------------------------------- */

void Multisphere::add_body(int nspheres, double *xcm_ins, double *xcm_to_xbound_ins,
               double r_bound_ins, double *v_ins, double *omega_ins, double mass_ins,double dens_ins,
               int atomtype_ins, int type_ins, double *inertia_ins,
               double *ex_space_ins, double *ey_space_ins, double *ez_space_ins,
               double **displace_ins, bool *fflag, bool *tflag, int start_step_ins,double *v_integrate_ins)
{
    
    int n = nbody_;

    customValues_.addUninitializedElement();

    // set initialize ID for element
    // ID starts from 0
    
    id_.set(n,-1);

    double zerovec[3] = {0.,0.,0.};
    double zerovec4[4] = {0.,0.,0.,0.};
    int zerovec4int[4] = {0,0,0,0};

    xcm_.set(n,xcm_ins);
    vcm_.set(n,v_ins);
    fcm_.set(n,zerovec);
    torquecm_.set(n,zerovec);
    dragforce_cm_.set(n,zerovec);
    hdtorque_cm_.set(n,zerovec);

    angmom_.set(n,zerovec);
    omega_.set(n,omega_ins);
    quat_.set(n,zerovec4);

    density_.set(n,dens_ins);
    atomtype_.set(n,atomtype_ins);
    type_.set(n,type_ins);
    masstotal_.set(n,mass_ins);
    inertia_.set(n,inertia_ins);
    ex_space_.set(n,ex_space_ins);
    ey_space_.set(n,ey_space_ins);
    ez_space_.set(n,ez_space_ins);

    nrigid_.set(n,nspheres);
    imagebody_.set(n,(IMGMAX << IMG2BITS) | (IMGMAX << IMGBITS) | IMGMAX);
    remapflag_.set(n,zerovec4int);

    fflag_.set(n,fflag);
    tflag_.set(n,tflag);

    start_step_.set(n,start_step_ins);
    if(v_integrate_ins)
        v_integrate_.set(n,v_integrate_ins);
    else
        v_integrate_.set(n,zerovec);

    r_bound_.set(n,r_bound_ins);
    xcm_to_xbound_.set(n,xcm_to_xbound_ins);

    // initialize the temperature with the initial value
    
    FixHeatGran *fix_heat = static_cast<FixHeatGran*>(modify->find_fix_style("heat/gran",0));
    if (fix_heat) {
        temp_.set(n,fix_heat->T0);
        temp_old_.set(n,fix_heat->T0);
    }
    else
    {
        temp_.set(n,0.);
        temp_old_.set(n,0.);
    }

    // calculate q and ang momentum

    MathExtra::exyz_to_q
    (
        ex_space_(n),ey_space_(n),ez_space_(n),
        quat_(n)
    );
    
    MathExtraLiggghts::angmom_from_omega
    (
        omega_(n),
        ex_space_(n),ey_space_(n),ez_space_(n),
        inertia_(n),
        angmom_(n)
    );

    // loop all non-initialized properties and set to their
    // default values
    int iProperty = 0;
    for(ContainerBase *cb = customValues_.getElementPropertyBase(iProperty); cb; cb = customValues_.getElementPropertyBase(++iProperty))
    {
        
        if(cb->useDefault())
        {
            
            cb->setToDefault(n);
            
        }
    }

    // increase local body counter
    nbody_++;

}

/* ----------------------------------------------------------------------
   remap bodies
------------------------------------------------------------------------- */

void Multisphere::remap_bodies(int *body)
{

  tagint original,oldimage,newimage;
  double xbnd[3],xbnd_old[3],xbnd_diff[3];

  // adjust body

  for (int ibody = 0; ibody < nbody_; ibody++) {
    original = imagebody_(ibody);

    x_bound(xbnd,ibody);
    vectorCopy3D(xbnd,xbnd_old);
    domain->remap(xbnd,imagebody_(ibody));
    vectorSubtract3D(xbnd,xbnd_old,xbnd_diff);
    vectorAdd3D(xcm_(ibody),xbnd_diff,xcm_(ibody));

    if (original == imagebody_(ibody)) remapflag_(ibody)[3] = 0;
    else {
      oldimage = original & IMGMASK;
      newimage = imagebody_(ibody) & IMGMASK;
      remapflag_(ibody)[0] = newimage - oldimage;
      oldimage = (original >> IMGBITS) & IMGMASK;
      newimage = (imagebody_(ibody) >> IMGBITS) & IMGMASK;
      remapflag_(ibody)[1] = newimage - oldimage;
      oldimage = original >> IMG2BITS;
      newimage = imagebody_(ibody) >> IMG2BITS;
      remapflag_(ibody)[2] = newimage - oldimage;
      remapflag_(ibody)[3] = 1;
    }
  }

  tagint *atomimage = atom->image;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  int ibody;
  tagint idim,otherdims;

  for (int i = 0; i < nlocal+nghost; i++)
  {

    if(body[i] < 0) continue;
    ibody = map(body[i]);

    if (ibody < 0) continue;
    if (remapflag_(ibody)[3] == 0) continue;

    if (remapflag_(ibody)[0]) {
      idim = atomimage[i] & IMGMASK;
      otherdims = atomimage[i] ^ idim;
      idim -= remapflag_(ibody)[0];
      idim &= IMGMASK;
      atomimage[i] = otherdims | idim;
    }
    if (remapflag_(ibody)[1]) {
      idim = (atomimage[i] >> IMGBITS) & IMGMASK;
      otherdims = atomimage[i] ^ (idim << IMGBITS);
      idim -= remapflag_(ibody)[1];
      idim &= IMGMASK;
      atomimage[i] = otherdims | (idim << IMGBITS);
    }
    if (remapflag_(ibody)[2]) {
      idim = atomimage[i] >> IMG2BITS;
      otherdims = atomimage[i] ^ (idim << IMG2BITS);
      idim -= remapflag_(ibody)[2];
      idim &= IMGMASK;
      atomimage[i] = otherdims | (idim << IMG2BITS);
    }
  }
}

/* ----------------------------------------------------------------------
   add unique ids to any body with id = -1
   new ids are grouped by proc and start after max current tag
   called after creating new atoms
------------------------------------------------------------------------- */

void Multisphere::id_extend_body_extend(int *body)
{
  int idmax, idmax_all;
  int nlocal = atom->nlocal;

  // calc total # of bodies
  calc_nbody_all();

  // return if no bodies are present
  if(nbody_all_ == 0)
    return;

  // idmax = max id for all bodies

  idmax = id_.max();
  MPI_Max_Scalar(idmax,idmax_all,this->world);

  // mapTagMax_ cannot get smaller - so ensure IDs are given only once
  
  mapTagMax_ = std::max(mapTagMax_,idmax_all);

  // noid = # of bodies I own with no id (id = -1)
  // noid_sum = # of total bodies on procs <= me with no tag
  // also check if # of atoms with no body id is consistent
  // nobody = # of atoms newly inserted with no body associated (body = -2)
  // nobody_check = # of atoms that should be newly inserted

  int noid = 0;
  int nobody = 0, nobody_first = 0, nobody_check = 0;
  for (int i = 0; i < nbody_; i++)
  {
    if (id_(i) == -1)
    {
        noid++;
        nobody_check += nrigid_(i);
    }
  }

  for(int i = 0; i < nlocal; i++)
  {
    if(body[i] == -2)
    {
       if(nobody == 0)
         nobody_first = i; 
       nobody++;
    }
  }

  if(nobody != nobody_check)
  {
    if(screen) fprintf(screen,"nobody: %d nobody_check: %d, nobody_first: %d. \n", nobody, nobody_check, nobody_first);
    error->one(FLERR,"Internal error: # of atoms with no associated body inconsistent");
  }

  int noid_sum;
  MPI_Scan(&noid,&noid_sum,1,MPI_INT,MPI_SUM,world);

  // itag = 1st new tag that my untagged bodies should use
  // give atoms body ID as well
  
  int itag = mapTagMax_ + noid_sum - noid + 1;
  for (int ibody = 0; ibody < nbody_; ibody++)
  {
    if (id_(ibody) == -1)
    {
        id_(ibody) = itag;
        
        if((nobody_first == nlocal-1) && ( nrigid_(ibody)>1 )) //allow body with a single atom
            error->one(FLERR,"Internal error: atom body id inconsistent: (nobody_first == nlocal-1) && ( nrigid_(ibody)>1 )");

        for(int iatom = nobody_first; iatom < nobody_first+nrigid_(ibody); iatom++)
        {
            if(body[iatom] != -2)
                error->one(FLERR,"Internal error: atom body id inconsistent");
            body[iatom] = itag;
            
        }
        nobody_first += nrigid_(ibody);

        // search for next particle with no body associated
        while(nobody_first < nlocal-1 && body[nobody_first] != -2)
            nobody_first++;

        itag++;
    }
  }
}

/* ----------------------------------------------------------------------
   clear and generate a global map for global-local lookup
------------------------------------------------------------------------- */

void Multisphere::clear_map()
{
    // deallocate old memory
    memory->destroy(mapArray_);
    mapArray_ = NULL;
}

void Multisphere::generate_map()
{
    int idmax, idmax_all;

    // deallocate old memory if exists
    if(mapArray_) clear_map();

    if(nbody_all_ == 0)
        return;

    // get max ID of all proc
    idmax = id_.max();
    MPI_Max_Scalar(idmax,idmax_all,world);
    mapTagMax_ = std::max(mapTagMax_,idmax_all);

    // alocate and initialize new array
    // IDs start at 1, have to go up to (inclusive) mapTagMax_
    
    memory->create(mapArray_,mapTagMax_+1,"Multisphere:mapArray_");
    for(int i = 0; i < mapTagMax_+1; i++)
        mapArray_[i] = -1;

    // build map
    for (int i = nbody_-1; i >= 0; i--)
    {
        
        mapArray_[id_(i)] = i;
    }
}

/* ----------------------------------------------------------------------
   check for lost atoms and bodies
------------------------------------------------------------------------- */

bool Multisphere::check_lost_atoms(int *body, double *atom_delflag, double *body_existflag, double *volumeweight)
{
    int body_tag,ibody,i;
    int nall = atom->nlocal + atom->nghost;
    int deleted = 0;

    int *nrigid_current = new int[nbody_];
    int *delflag = new int[nbody_];
    vectorZeroizeN(nrigid_current,nbody_);
    vectorZeroizeN(delflag,nbody_);

    //double _4pi_over_3 = 4.*M_PI/3.;

    for(i = 0; i < nall; i++)
    {
        body_tag = body[i];
        
        if(body_tag >= 0 && map(body_tag) >= 0 && domain->is_owned_or_first_ghost(i))
        {
            nrigid_current[map(body_tag)]++;
            
            body_existflag[i] = 1.;
        }
        else if (body_tag == -1)
            body_existflag[i] = 1.;
    }

    for(ibody = 0; ibody < nbody_; ibody++)
    {
        if(nrigid_current[ibody] > nrigid_(ibody))
        {
            
            error->one(FLERR,"Internal error in multisphere method");
        }
        if(nrigid_current[ibody] != nrigid_(ibody))
        {
            
            delflag[ibody] = 1;
            
        }
    }

    for(i = 0; i < nall; i++)
    {
        body_tag = body[i];

        if(body_tag < 0) continue;

        ibody = map(body_tag);

        if(ibody >= 0 && delflag[ibody])
        {
           
           atom_delflag[i] = 1.;
           body[i] = -1;
           deleted = 1;
           
           atom->rmass[i] *= volumeweight[i];
        }
    }

    ibody = 0;
    while(ibody < nbody_)
    {
        if(delflag[ibody] == 1)
        {
            
            delflag[ibody] = delflag[nbody_-1];
            remove_body(ibody);
            
        }
        else ibody++;
    }

    calc_nbody_all();

    MPI_Max_Scalar(deleted,world);

    delete []nrigid_current;
    delete []delflag;

    if(deleted == 1) return true;
    return false;
}

/* ----------------------------------------------------------------------
   restart - not available in PUBLIC
------------------------------------------------------------------------- */

void Multisphere::writeRestart(FILE *)
{
    error->one(FLERR,"Multisphere write_restart is not available in your version. See www.cfdem.com for details");
}

void Multisphere::restart(double *)
{
    error->one(FLERR,"Multisphere restart is not available in your version. See www.cfdem.com for details");
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   if don't recognize name, return NULL

   important: returns GLOBAL lengths (tag_max)
              because this is required by CfdDatacouplingMPI
------------------------------------------------------------------------- */

void *Multisphere::extract(const char *name, int &len1, int &len2)
{
    // scalars

    len1 = len2 = 1;

    if (strcmp(name,"nbody") == 0) return (void *) &nbody_;
    if (strcmp(name,"nbody_all") == 0) return (void *) &nbody_all_;

    // per-body properties

    len1 = mapTagMax_;
    ContainerBase *cb = customValues_.getElementPropertyBase(name);

    if(NULL == cb)
    {
        len1 = len2 = -1;
        return NULL;
        /*
        fprintf(screen,"ERROR: Property %s not found, which is required for multi-sphere\n",name);
        error->all(FLERR,"Property required for multi-sphere not found");
        */
    }

    len2 = cb->lenVec();
    if(cb->nVec() != 1)
       error->all(FLERR,"Internal error, cannot use multi-vector containers");
    return cb->begin_slow_dirty();
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   only return if data of type double**
   if don't recognize name or not double ** data, return NULL
------------------------------------------------------------------------- */

double** Multisphere::extract_double_vector(const char *name)
{
    VectorContainer<double,3> *cb = customValues_.getElementProperty<VectorContainer<double,3> >(name);
    if(!cb) return 0;
    return cb->begin();
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   only return if data of type double*
   if don't recognize name or not double * data, return NULL
------------------------------------------------------------------------- */

double* Multisphere::extract_double_scalar(const char *name)
{
    ScalarContainer<double> *cb = customValues_.getElementProperty<ScalarContainer<double> >(name);
    if(!cb) return 0;
    return cb->begin();
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
------------------------------------------------------------------------- */

double Multisphere::extract_ke()
{
  double ke = 0.0;
  double mvv2e = force->mvv2e;

  for (int i = 0; i < nbody_; i++)
    ke += masstotal_(i)*vectorMag3DSquared(vcm_(i));

  MPI_Sum_Scalar(ke,world);

  return 0.5*mvv2e*ke;
}

double Multisphere::extract_vave()
{
  double vave = 0.0;
  for (int i = 0; i < nbody_; i++)
    vave += vectorMag3D(vcm_(i));

  MPI_Sum_Scalar(vave,world);
  return vave/nbody_all_;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

double Multisphere::extract_rke()
{
  double wbody[3],rot[3][3];

  double rke = 0.0;
  for (int i = 0; i < nbody_; i++) {

    // wbody = angular velocity in body frame

    MathExtra::quat_to_mat(quat_(i),rot);
    MathExtra::transpose_matvec(rot,angmom_(i),wbody);
    if (inertia_(i)[0] == 0.0) wbody[0] = 0.0;
    else wbody[0] /= inertia_(i)[0];
    if (inertia_(i)[1] == 0.0) wbody[1] = 0.0;
    else wbody[1] /= inertia_(i)[1];
    if (inertia_(i)[2] == 0.0) wbody[2] = 0.0;
    else wbody[2] /= inertia_(i)[2];

    rke += inertia_(i)[0]*wbody[0]*wbody[0] +
      inertia_(i)[1]*wbody[1]*wbody[1] + inertia_(i)[2]*wbody[2]*wbody[2];
  }

  MPI_Sum_Scalar(rke,world);

  return 0.5*rke;
}

double Multisphere::extract_omega_ave()
{
  double omega_ave = 0.0;
  for (int i = 0; i < nbody_; i++)
    omega_ave += vectorMag3D(omega_(i));

  MPI_Sum_Scalar(omega_ave,world);
  return omega_ave/nbody_all_;
}
