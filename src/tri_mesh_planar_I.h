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
