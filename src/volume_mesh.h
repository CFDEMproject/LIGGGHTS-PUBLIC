/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_VOLUME_MESH_H
#define LMP_VOLUME_MESH_H

#include "tracking_mesh.h"
#include "container.h"

namespace LAMMPS_NS{

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
class VolumeMesh : public TrackingMesh<NUM_NODES>
{
  public:

    bool addElement(double **nodeToAdd);

    void move(double *vecTotal, double *vecIncremental);
    void move(double *vecIncremental);
    void scale(double factor);

    bool isInside(double *p);

    virtual int generateRandomOwnedGhost(double *pos) = 0;
    virtual int generateRandomSubbox(double *pos) = 0;
    virtual int generateRandomSubboxWithin(double *pos,double delta) = 0;

    // public inline access

    // area of total mesh - all elements (all processes)
    
    inline double volMeshGlobal()
    { return volMesh_(0);}

    // area of owned elements
    inline double volMeshOwned()
    { return volMesh_(1);}

    // area of ghost elements
    inline double volMeshGhost()
    { return volMesh_(2);}

    // area of owned and ghost elements in my subdomain
    inline double volMeshSubdomain()
    { return volMesh_(3);}

  protected:

    VolumeMesh(LAMMPS *lmp);
    virtual ~VolumeMesh();

    void deleteElement(int n);

    void buildNeighbours();

    void refreshOwned(int setupFlag);
    void refreshGhosts(int setupFlag);

    inline void recalcLocalVolProperties();
    inline void recalcGhostVolProperties();

    void calcVolPropertiesOfNewElement();

    virtual bool shareFace(int i, int j, int &iFace, int &jFace) = 0;

    virtual bool isInside(int nElem, double *p) = 0;
    virtual double calcVol(int nElem) = 0;

    int randomOwnedGhostElement();

    void rotate(double *totalQ, double *dQ,double *totalDispl,double *dDispl);
    void rotate(double *dQ,double *dDispl);

    // inline access

    inline double&  vol(int i)
    { return (vol_)(i); }

    inline double& volAcc(int i)
    { return (volAcc_)(i); }

    inline double** faceNodes(int i)
    { return faceNodes_[i]; }

  private:

    void checkOrientation(int n);
    void calcFaceNormals(int n);

    int searchElementByVolAcc(double vol,int lo, int hi);

    // mesh properties

    ScalarContainer<double>& volMesh_; 

    // volume and accumulated volume for each element

    ScalarContainer<double> &vol_;
    ScalarContainer<double> &volAcc_;

    // faces for each element

    MultiVectorContainer<int,NUM_FACES,NUM_NODES_PER_FACE>& faceNodes_;
    MultiVectorContainer<double,NUM_FACES,3>& faceNormals_;
    VectorContainer<bool,NUM_FACES>& isBoundaryFace_;

    // neighbor topology for each element

    ScalarContainer<int>& nNeighs_;
    VectorContainer<int,NUM_FACES>& neighElems_;
};

// *************************************
#include "volume_mesh_I.h"
// *************************************

} /* LAMMPS_NS */

#endif /* VOLUMEMESH_H_ */
