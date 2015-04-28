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
