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

#include "tri_mesh_planar.h"

/* ----------------------------------------------------------------------
   constructor, destructor
------------------------------------------------------------------------- */

TriMeshPlanar::TriMeshPlanar(LAMMPS *lmp) :
    TriMesh(lmp),

    nearestActiveEdgeID_   (*this->prop().addElementProperty< VectorContainer<int,2*NUM_NODES> >("nearestActiveEdgeID","comm_exchange_borders","frame_invariant","restart_no")),
    nearestActiveEdgeIndex_(*this->prop().addElementProperty< VectorContainer<int,2*NUM_NODES> >("nearestActiveEdgeIndex","comm_exchange_borders","frame_invariant","restart_no")),
    minActiveEdgeDist_     (*this->prop().addElementProperty< ScalarContainer<double> >         ("minActiveEdgeDist","comm_exchange_borders","frame_invariant","restart_no"))
{
    
    doParallellization_ = false;;
}

TriMeshPlanar::~TriMeshPlanar()
{
}

/* ----------------------------------------------------------------------
   build up edge lists
------------------------------------------------------------------------- */

void TriMeshPlanar::postInitialSetup()
{
    // works for non-parallel mesh only
    if(!isPlanar())
        error->all(FLERR,"Face defined as planar face is in fact not planar. You might want to check the 'curvature' setting");

    if(this->isParallel())
        error->all(FLERR,"internal error");

    buildEdgeLists();
}

/* ----------------------------------------------------------------------
   build up edge lists
------------------------------------------------------------------------- */

void TriMeshPlanar::buildEdgeLists()
{
    double dist;
    const double INIT_LARGE = 1e8;
    int nall = this->sizeLocal()+this->sizeGhost();
    int nEdges, iLrg, dummy;

    // distances to the 6 edges that are stored
    double distEdges[2*NUM_NODES];

    // initialize containers
    nearestActiveEdgeID_.setAll(-1);
    nearestActiveEdgeIndex_.setAll(-1);
    minActiveEdgeDist_.setAll(0.);

    for(int iSrf = 0; iSrf < nall; iSrf++)
    {
        // initially no edges stored for iSrf
        // initialize distances with LARGE
        nEdges = 0;
        for(int i = 0; i < 2*NUM_NODES; i++)
            distEdges[i] = INIT_LARGE;

        for(int jSrf = 0; jSrf < nall; jSrf++)
        {
            
            if(n_active_edges(jSrf) == 0)
                continue;

            for(int iEdge = 0; iEdge < NUM_NODES; iEdge++)
            {
                for(int jEdge = 0; jEdge < NUM_NODES; jEdge++)
                {
                    
                    if(!edgeActive(jSrf,jEdge))
                        continue;

                    dist = edgeNodeDist(iSrf,iEdge,jSrf,jEdge);

                    if(nEdges < 2*NUM_NODES)
                    {
                        nearestActiveEdgeID_(iSrf)[nEdges] = TrackingMesh<NUM_NODES>::id(jSrf);
                        nearestActiveEdgeIndex_(iSrf)[nEdges] = jEdge;
                        distEdges[nEdges] = dist;
                        nEdges++;
                    }
                    
                    else
                    {
                        MathExtraLiggghts::max(distEdges,2*NUM_NODES,iLrg);
                        nearestActiveEdgeID_(iSrf)[iLrg] = TrackingMesh<NUM_NODES>::id(jSrf);
                        nearestActiveEdgeIndex_(iSrf)[iLrg] = jEdge;
                        distEdges[iLrg] = dist;
                    }
                }
            }
        }

        // store the minimum distance for quick check
        minActiveEdgeDist_(iSrf) = MathExtraLiggghts::min(distEdges,2*NUM_NODES,dummy);
    }
}
