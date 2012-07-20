/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_VOLUME_MESH_I_H
#define LMP_VOLUME_MESH_I_H

#define NTRY_MC_VOLUME_MESH_I_H 30000
#define NITER_MC_VOLUME_MESH_I_H 5
#define TOLERANCE_MC_VOLUME_MESH_I_H 0.05

/* ----------------------------------------------------------------------
   constructor(s), destructor
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
VolumeMesh<NUM_NODES,N_FACES>::VolumeMesh()
:   TrackingMesh<NUM_NODES>(),
    isInsertionMesh_(false),

    // TODO should keep volMeshSubdomain up-to-date more often for insertion faces
    volMesh_(*this->prop().template addMeshProperty   < ScalarContainer<double> >  ("volMesh", "comm_none","frame_trans_rot_invariant","restart_no",3)),

    vol_    (*this->prop().template addElementProperty< ScalarContainer<double> > ("vol",     "comm_none","frame_trans_rot_invariant", "restart_no",3)),
    volAcc_ (*this->prop().template addElementProperty< ScalarContainer<double> > ("volAcc",  "comm_none","frame_trans_rot_invariant", "restart_no",3))
{
    
    volMesh_.add(0.);
    volMesh_.add(0.);
    volMesh_.add(0.);
    volMesh_.add(0.);
}

template<int NUM_NODES,int N_FACES>
VolumeMesh<NUM_NODES,N_FACES>::~VolumeMesh()
{}

/* ----------------------------------------------------------------------
   set flag if used as insertion mesh
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::useAsInsertionMesh()
{
    isInsertionMesh_ = true;
}

/* ----------------------------------------------------------------------
   add / delete element
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::addElement(double **nodeToAdd)
{
    TrackingMesh<NUM_NODES>::addElement(nodeToAdd);

    calcVolPropertiesOfNewElement();
}

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::refreshOwned(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshOwned(setupFlag);
    // (re)calculate all properties for owned elements
    
    recalcLocalVolProperties();
}
template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::refreshGhosts(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshGhosts(setupFlag);

    recalcGhostVolProperties();
}

/* ----------------------------------------------------------------------
   recalculate properties of local elements
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::recalcLocalVolProperties()
{
    
    // volMeshGlobal [volMesh_(0)] and volMeshOwned [volMesh_(1)]
    // calculated here

    double areaAccOff;

    volMesh_(0) = 0.;
    volMesh_(1) = 0.;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
      vol(i) = calcVol(i);
      volAcc(i) = vol(i);
      if(i > 0) volAcc(i) += volAcc(i-1);

      // add to local volume
      volMesh_(1) += vol(i);
      
    }

    // mesh vol must be summed up
    MPI_Sum_Scalar(volMesh_(1),volMesh_(0),this->world);

}

/* ----------------------------------------------------------------------
   recalculate properties of ghost elements
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::recalcGhostVolProperties()
{
    double pos[3], areaCheck;
    int n_succ, n_iter;
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();

    // volMeshGhost [volMesh_(2)] and volMeshSubdomain [volMesh_(3)]
    // calculated here

    // accumulated vol includes owned and ghosts
    volMesh_(2) = 0.;
    for(int i = nlocal; i < nall; i++)
    {
      vol(i) = calcVol(i);
      volAcc(i) = vol(i);
      if(i > 0) volAcc(i) += volAcc(i-1);

      // add to ghost area
      volMesh_(2) += vol(i);
    }

    // calc area of owned and ghost elements in my subdomain
    
    volMesh_(3) = 0.;
    double volCheck = 0.;

    if(isInsertionMesh_)
    {
        n_succ = 0;
        n_iter = 0;

        // iterate long enough so MC has the desired tolerance
        while( (n_iter < NITER_MC_VOLUME_MESH_I_H) &&
               (fabs((volCheck-volMeshGlobal()))/volMeshGlobal() > TOLERANCE_MC_VOLUME_MESH_I_H) )
        {
            // only generate random positions if I have any mesh elements
            if(nall)
            {
                for(int i = 0; i < NTRY_MC_VOLUME_MESH_I_H; i++)
                {
                    // pick a random position on owned or ghost element
                    if((generateRandomOwnedGhost(pos) >= 0) && (this->domain->is_in_subdomain(pos)))
                        n_succ++;
                }
            }
            n_iter++;
            volMesh_(3) = static_cast<double>(n_succ)/static_cast<double>(NTRY_MC_VOLUME_MESH_I_H*n_iter) * (volMeshOwned()+volMeshGhost());

            MPI_Sum_Scalar(volMesh_(3),volCheck,this->world);
        }

        if(fabs((volCheck-volMeshGlobal()))/volMeshGlobal() > TOLERANCE_MC_VOLUME_MESH_I_H)
            this->error->all(FLERR,"Local mesh volume calculation failed, try boosting NITER_MC_VOLUME_MESH_I_H");

        // correct so sum of all owned vols is equal to global area
        volMesh_(3) *= volMeshGlobal()/volCheck;
    }

}

/* ----------------------------------------------------------------------
   generate a random Element by volAcc
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
int VolumeMesh<NUM_NODES,N_FACES>::randomOwnedGhostElement()
{
    
    if(!isInsertionMesh_) this->error->one(FLERR,"Illegal call for non-insertion mesh");
    double r = this->random_->uniform() * (volMeshOwned()+volMeshGhost());
    int nall = this->sizeLocal()+this->sizeGhost()-1;
    return searchElementByVolAcc(r,0,nall);
}

template<int NUM_NODES,int N_FACES>
int VolumeMesh<NUM_NODES,N_FACES>::searchElementByVolAcc(double vol,int lo, int hi)
{
    if( (lo < 1 || vol > volAcc(lo-1)) && (vol <= volAcc(lo)) )
        return lo;
    if( (hi < 1 || vol > volAcc(hi-1)) && (vol <= volAcc(hi)) )
        return hi;

    int mid = static_cast<int>((lo+hi)/2);
    if(vol > volAcc(mid))
        return searchElementByVolAcc(vol,mid,hi);
    else
        return searchElementByVolAcc(vol,lo,mid);
}

/* ----------------------------------------------------------------------
   calculate properties when adding new element
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::calcVolPropertiesOfNewElement()
{
    
    int n = MultiNodeMesh<NUM_NODES>::node.size()-1;

    // calc volume
    double vol_elem = calcVol(n);
    volMesh_(0) += vol_elem;
    vol_(n) = vol_elem;
    volAcc_(n) = vol_elem;
    if(n > 0) volAcc_(n) += volAcc_(n-1);
}

/* ----------------------------------------------------------------------
   build neighlist, generate mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::buildNeighbours()
{
    // iterate over all elems, over ghosts as well
    int nall = this->sizeLocal()+this->sizeGhost();

    // inititalize neigh topology - reset to default, ~n

    int neighs[NUM_NODES];
    for(int i=0;i<NUM_NODES;i++)
        neighs[i] = -1;

    for(int i = 0; i < nall; i++)
    {
        nNeighs_.set(i,0);
        neighElems_.set(i,neighs);
    }

    // build neigh topology, ~n*n/2
    for(int i = 0; i < nall; i++)
    {
      for(int j = i+1; j < nall; j++)
      {
        
        int iNode(0), jNode(0), iEdge(0), jEdge(0);
        if(!this->shareNode(i,j,iNode,jNode)) continue;

        this->error->all(FLERR,"TODO");
        //if(shareFace(i,j))
        {
            neighElems_(i)[nNeighs_(i)] = this->id(i);
            neighElems_(j)[nNeighs_(j)] = this->id(j);
            nNeighs_(i)++;
            nNeighs_(j)++;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   isInside etc
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
bool VolumeMesh<NUM_NODES,N_FACES>::isInside(double *p)
{
    // check subdomain
    if(!this->domain->is_in_subdomain(p)) return false;

    // check bbox
    if(!this->bbox_.isInside(p)) return false;

    int nall = this->size();

    // brute force
    for(int i = 0; i < nall; i++)
        if(isInside(i,p)) return true;

    return false;
}

/* ----------------------------------------------------------------------
   move, rotate, scale mesh
------------------------------------------------------------------------- */

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::move(double *vecTotal, double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecTotal,vecIncremental);
}

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::move(double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecIncremental);
}

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::scale(double factor)
{
    TrackingMesh<NUM_NODES>::scale(factor);
}

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::rotate(double *totalQ, double *dQ,double *totalDispl, double *dDisp)
{
    TrackingMesh<NUM_NODES>::rotate(totalQ,dQ,totalDispl,dDisp);

}

template<int NUM_NODES,int N_FACES>
void VolumeMesh<NUM_NODES,N_FACES>::rotate(double *dQ,double *dDispl)
{
    TrackingMesh<NUM_NODES>::rotate(dQ,dDispl);
}

#endif
