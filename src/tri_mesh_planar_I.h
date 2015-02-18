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

#ifndef LMP_TRI_MESH_PLANAR_I_H
#define LMP_TRI_MESH_PLANAR_I_H

/* ----------------------------------------------------------------------
 generates a random point on the surface that lies within my subbox
 additionally, point is within a distance of delta from active edges
 so can generate a point within a distance delta from boundary on planar meshes
 tests against active edges
------------------------------------------------------------------------- */

inline int TriMeshPlanar::generateRandomOwnedGhostWithin(double *pos,double delta)
{
    int i, iSrf, iSrfCheck, iSrfCheckEdge, ntry;
    double dist;
    bool farEnoughAway;

    ntry = 0;

    do
    {
        ntry++;

        iSrf = generateRandomOwnedGhost(pos);
        farEnoughAway = true;
        
        if(minActiveEdgeDist_(iSrf) < delta)
        {
            
            i = 0;
            while(i < 2*NUM_NODES && nearestActiveEdgeID_(iSrf)[i] >= 0 && farEnoughAway)
            {
                iSrfCheck = TrackingMesh<NUM_NODES>::map
                (
                    nearestActiveEdgeID_(iSrf)[i]
                );
                iSrfCheckEdge = nearestActiveEdgeIndex_(iSrf)[i];

                dist = edgePointDist(iSrfCheck,iSrfCheckEdge,pos);

                if(dist < delta)
                    farEnoughAway = false;

                i++;
            }
        }
    }
    while(!farEnoughAway && ntry < MAXTRY_WITHIN);

    if(ntry == MAXTRY_WITHIN)
        error->warning(FLERR,"'All_in yes' did not succeed");

    return iSrf;
}

/* ----------------------------------------------------------------------
 search for position within the mesh - return ID, bary coords and normal dist
------------------------------------------------------------------------- */

inline bool TriMeshPlanar::locatePosition(double *pos,int &triID,double *bary,double &distance)
{
    double p[3], pProjected[3], temp[3];
    int nlocal = sizeLocal();
    vectorCopy3D(pos,p);

    // calc normal distance
    vectorSubtract3D(p,center_(0),temp);
    distance = vectorDot3D(temp,surfaceNorm(0));

    // calc projected point
    vectorScalarMult3D(surfaceNorm(0),distance,temp);
    vectorSubtract3D(p,temp,pProjected);

    for(int i = 0; i < nlocal; i++)
    {
        if(isInElement(pProjected,i))
        {
            triID = TrackingMesh<3>::id(i);
            vectorSubtract3D(pProjected,node_(i)[0],temp);
            MathExtraLiggghts::calcBaryTriCoords(temp,edgeVec(i),edgeLen(i),bary);
            
            return true;
        }
    }

    return false;
}

/* ----------------------------------------------------------------------
 calculate position from tri ID and bary coords
------------------------------------------------------------------------- */

inline bool TriMeshPlanar::constructPositionFromBary(int triID,double *bary,double *pos)
{
    int itri = TrackingMesh<3>::map(triID);

    if(itri < 0) return false;

    double tmp[3];
    vectorZeroize3D(pos);

    for(int i = 0; i < 3; i++)
    {
        vectorScalarMult3D(node_(itri)[i],bary[i],tmp);
        vectorAdd3D(pos,tmp,pos);
    }

    return true;
}

#endif
