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

#endif
