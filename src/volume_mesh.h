/*
 * VolumeMesh.h
 *
 *  Created on: 12.09.2011
 *      Author: phil
 */

#ifndef LMP_VOLUME_MESH_H_
#define LMP_VOLUME_MESH_H_

#include "tracking_mesh.h"
#include "container.h"

namespace LAMMPS_NS{

template<int NUM_NODES>
class VolumeMesh : public TrackingMesh<NUM_NODES>
{
  public:

    void calcVolPropertiesOfNewElement();
    bool isInside(double *p);
    void generateRandom(double *p);
    double volume();
    virtual void move(double *vec);

  protected:
    VolumeMesh();
    virtual ~VolumeMesh();

    virtual bool isInside(int nElem, double *p) =0;
    virtual double calcVol(int nElem) =0;
    virtual double calcCenter(int nElem) =0;

    int randomElement();
    virtual void generateRandom(int n, double *p) =0;

    virtual void rotate(double *totalQ, double *dQ,double *displacement);

    // mesh properties
    double volMesh;

    // per-element properties

    // add more properties as they are needed, such as
    // surface area, edge length ...
    ScalarContainer<double> &vol;
    ScalarContainer<double> &volAcc;
};

/* ----------------------------------------------------------------------
   constructor(s), destructor
------------------------------------------------------------------------- */

template<int NUM_NODES>
VolumeMesh<NUM_NODES>::VolumeMesh()
: TrackingMesh<NUM_NODES>(),
  vol_         (this->template addProperty< ScalarContainer<double> >("vol","comm_none","ref_trans_rot_invariant")),
  volAcc_      (this->template addProperty< ScalarContainer<double> >("volAcc","comm_none","ref_trans_rot_invariant")),
{

}

template<int NUM_NODES>
VolumeMesh<NUM_NODES>::~VolumeMesh()
{

}

/* ----------------------------------------------------------------------
   add / delete element
------------------------------------------------------------------------- */

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::addElement(double **nodeToAdd)
{
    TrackingMesh<NUM_NODES>::addElement(nodeToAdd);
    calcVolPropertiesOfNewElement();
}

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   calculate properties when adding new element
------------------------------------------------------------------------- */

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::calcVolPropertiesOfNewElement()
{
    double *vecTmp3 = create<double>(vecTmp3,3);

    int n = MultiNodeMesh<NUM_NODES>::node.size()-1;

    // calc volume
    double vol_elem = calcVol(n);
    volMesh += vol_elem;
    vol.set(n,vol_elem);
    volAcc.set(n,area_elem);
    if(n > 0) volAcc(n) += volAcc(n-1);

    destroy<double>(vecTmp3);
}

/* ----------------------------------------------------------------------
   re-calculate properties during simulation
------------------------------------------------------------------------- */

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::recalcVolProperties()
{
    recalcVol();
    recalcCenter();
}

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::recalcVol()
{
    volMesh = 0.;

    for(int i=0;i<size();i++){
      vol(i) = calcVol(i);
      volAcc(i) = vol(i);
      if(i > 0) volAcc(i) += volAcc(i-1);
      volMesh += area(i);
    }
}

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::recalcCenter()
{
    for(int i=0;i<center.size();i++){
      calcCenter(i, center(i));
    }
}

/* ----------------------------------------------------------------------
   isInside etc
------------------------------------------------------------------------- */

template<int NUM_NODES>
bool VolumeMesh<NUM_NODES>::isInside(double *p)
{
    // check subdomain
    if(!domain->is_in_subdomain(p)) return false;

    // check bbox
    if(!box_.isInside(p) return false;

    // brute force
    for(int i=0;i<size();i++)
        if(isInside(i,p)) return true;

    return false;
}

template<int NUM_NODES>
double VolumeMesh<NUM_NODES>::volume()
{
    return volMesh;
}

/* ----------------------------------------------------------------------
   random generation functions
------------------------------------------------------------------------- */

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::generateRandom(double *p)
{
    int n = randomElement();

}

template<int NUM_NODES>
int VolumeMesh<NUM_NODES>::randomElement()
{
    // primitive implementation, could do a binary search here
    double rd = volMesh * random->uniform();
    int chosen = 0;
    while (rd > volAcc(chosen) && chosen < size()-1)
        chosen++;
    return chosen;
}

/* ----------------------------------------------------------------------
   functions for altering the mesh on loading
------------------------------------------------------------------------- */

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::move(double *vec)
{
    TrackingMesh<NUM_NODES>::move(vec);
    // specific stuff to come here
    // center needs to be moved
    // other properties are translation-invariant
    for(int i=0;i<center.size();i++)
    {
      vectorAdd3D(center(i),vec,center(i));
    }
}

template<int NUM_NODES>
void SurfaceMesh<NUM_NODES>::scale(double factorX, double factorY, double factorZ)
{
    TrackingMesh<NUM_NODES>::scale(factorX,factorY,factorZ);
    areaMesh = 0.;
    for(int i=0;i<center.size();i++){
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcCenter(i, center(i));
      calcEdgeNormals(i, edgeNorm(i));
      area(i) = calcArea(i);
      areaMesh += area(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);
    }
}

template<int NUM_NODES>
void VolumeMesh<NUM_NODES>::rotate(double *totalQ, double *dQ,double *displacement)
{
      TrackingMesh<NUM_NODES>::rotate(totalQ, dQ, displacement);
}

} /* LAMMPS_NS */

#endif /* VOLUMEMESH_H_ */
