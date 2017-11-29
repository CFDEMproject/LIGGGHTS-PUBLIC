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

    Alexander Podlozhnyuk
    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_I_SUPERQUADRIC_H
#define LMP_TRI_MESH_I_SUPERQUADRIC_H

#ifdef SUPERQUADRIC_ACTIVE_FLAG

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSuperquadricContact(int nTri,
      double *delta, double *contactPoint, Superquadric particle)
  {
      double bary[3];
      return resolveTriSuperquadricContact(nTri, delta, contactPoint, particle, bary);
  }

//particle-triangle contact detection algorithm
  /* ---------------------------------------------------------------------- */
  inline double TriMesh::resolveTriSuperquadricContact(int nTri,
      double *delta, double *contactPoint, Superquadric particle, double *bary)
  {
    double **n = node_(nTri);
    double surfNorm[3];
    vectorCopy3D(SurfaceMeshBase::surfaceNorm(nTri), surfNorm);

    double pointOfMaximalPenetration[3], pointOfLowestPotential[3];
    bool particle_plane_intersection = particle.plane_intersection(surfNorm, n[0], pointOfMaximalPenetration, pointOfLowestPotential);
    if(particle_plane_intersection) {
      int flag = superquadricTriangleIntersection(nTri, pointOfLowestPotential, particle);
      double node0ToPoint[3];
      if(flag >= 0) {
        double f_node[3];
        f_node[0] = particle.shape_function_global(n[0]);
        f_node[1] = particle.shape_function_global(n[1]);
        f_node[2] = particle.shape_function_global(n[2]);

        double contactPointFace[3];
        double dFace = MathExtraLiggghtsNonspherical::point_wall_projection(surfNorm, n[0], pointOfMaximalPenetration, contactPointFace); //overlap distance
#ifdef LIGGGHTS_DEBUG
        if(std::isnan(dFace))
          this->error->one(FLERR,"dFace is NaN!");
#endif
        double dist[3];
        vectorSubtract3D(particle.center, n[0], dist);
        if(vectorDot3D(dist, surfNorm) > 0.0)
          vectorNegate3D(surfNorm); //overlap direction

        vectorSubtract3D(contactPointFace,n[0],node0ToPoint);
        double baryFace[3];
        MathExtraLiggghts::calcBaryTriCoords(node0ToPoint,edgeVec(nTri),edgeLen(nTri),baryFace);
        if(dFace < 0.001*particle.shape[0] or flag == 0) {
          vectorCopy3D(surfNorm, delta);
          vectorCopy3D(contactPointFace, contactPoint);
          vectorCopy3D(baryFace, bary);
          return -dFace;
        }
        else { //particle -edge/corner contact
          double closestPoints[3][3];
          double contactPointEdge[3];
          int iEdge, iCorner;
          bool edge_contact = false;
          bool corner_contact = false;
          double f_edge[3];
          f_edge[0] = particle.line_intersection(n[0], n[1], closestPoints[0]);
          f_edge[1] = particle.line_intersection(n[1], n[2], closestPoints[1]);
          f_edge[2] = particle.line_intersection(n[2], n[0], closestPoints[2]);
          if(f_edge[0] > 0.0 and f_edge[1] > 0.0 and f_edge[2] > 0.0) { //no particle-edge contact
            vectorCopy3D(surfNorm, delta);
            vectorCopy3D(contactPointFace, contactPoint);
            return -dFace;
          }

          if(f_edge[0] <= std::min(f_edge[1], f_edge[2]))
            iEdge = 0;
          else if(f_edge[1] <= std::min(f_edge[0], f_edge[2]))
            iEdge = 1;
          else
            iEdge = 2;

          int ip = iEdge;
          int ipp = (iEdge+1)%3;
          if(f_node[ip] <= f_edge[iEdge]) {
            iCorner = ip; //corner contact
            corner_contact = true;
          }
          else if(f_node[ipp] <= f_edge[iEdge]) {
            iCorner = ipp; //corner contact
            corner_contact = true;
          }
          else //edge contact
            edge_contact = true;

          if(edge_contact && !corner_contact) { //check if the edge is active
            if(edgeActive(nTri)[iEdge])
              LAMMPS_NS::vectorCopy3D(closestPoints[iEdge], contactPointEdge);
            else
              return LARGE_TRIMESH;
          } else if(!edge_contact && corner_contact){ //check if corner is active
            if(cornerActive(nTri)[iCorner])
              LAMMPS_NS::vectorCopy3D(n[iCorner], contactPointEdge);
            else
              return LARGE_TRIMESH;
          } else
            this->error->one(FLERR,"Internal error"); //should not be

          vectorSubtract3D(contactPointEdge,n[0],node0ToPoint);
          double baryEdge[3];
          MathExtraLiggghts::calcBaryTriCoords(node0ToPoint,edgeVec(nTri),edgeLen(nTri),baryEdge);

          particle.map_point(contactPointEdge, pointOfMaximalPenetration);
          double normal[3];
          particle.shape_function_gradient_global(contactPointEdge, normal);
          vectorNormalize3D(normal);
          if(edge_contact) {
            double direction[3];
            LAMMPS_NS::vectorSubtract3D(n[ip], n[ipp], direction);
            vectorNormalize3D(direction);
            double cosa = vectorDot3D(normal, direction);
            for(int i = 0; i < 3; i++)
              normal[i] = normal[i] - cosa*direction[i];
            vectorNormalize3D(normal);
          }
          double dEdge = particle.surface_line_intersection1(contactPointEdge, normal, pointOfMaximalPenetration);
#ifdef LIGGGHTS_DEBUG
          if(std::isnan(dEdge))
            this->error->one(FLERR,"dEdge is NaN!");
#endif
          if(dEdge < dFace) {
            vectorCopy3D(normal, delta);
            vectorCopy3D(contactPointEdge, contactPoint);
            vectorCopy3D(baryEdge, bary);
            return -dEdge; //return "overlap" distance
          } else {
            vectorCopy3D(surfNorm, delta);
            vectorCopy3D(contactPointFace, contactPoint);
            vectorCopy3D(baryFace, bary);
            return -dFace;
          }
        }
      } else
        return LARGE_TRIMESH;
    } else
      return LARGE_TRIMESH;
  }

//calculates distance from a point to a triangle
  inline double TriMesh::pointToTriangleDistance(int nTri, double *cSphere, double *delta, bool treatActiveFlag, double *bary)
  {
    double **n = node_(nTri);
    int obtuseAngleIndex = SurfaceMeshBase::obtuseAngleIndex(nTri);

    bary[0] = bary[1] = bary[2] = 0.;

    double node0ToSphereCenter[3];
    //double *surfNorm = SurfaceMeshBase::surfaceNorm(nTri);
    vectorSubtract3D(cSphere,n[0],node0ToSphereCenter);

    MathExtraLiggghts::calcBaryTriCoords(node0ToSphereCenter,edgeVec(nTri),edgeLen(nTri),bary);

    double invlen = 1./(2.*rBound_(nTri));
    int barySign = (bary[0] > -precision_trimesh()*invlen) + 2*(bary[1] > -precision_trimesh()*invlen) + 4*(bary[2] > -precision_trimesh()*invlen);

    double d(0.);

    switch(barySign)
    {
    case 1:
      d = resolveCornerContactBary(nTri,0,obtuseAngleIndex == 0,cSphere,delta,bary, treatActiveFlag);
      break;
    case 2:
      d = resolveCornerContactBary(nTri,1,obtuseAngleIndex == 1,cSphere,delta,bary, treatActiveFlag);
      break;
    case 3:
      d = resolveEdgeContactBary(nTri,0,cSphere,delta,bary, treatActiveFlag);
      break;
    case 4:
      d = resolveCornerContactBary(nTri,2,obtuseAngleIndex == 2,cSphere,delta,bary, treatActiveFlag);
      break;
    case 5:
      d = resolveEdgeContactBary(nTri,2,cSphere,delta,bary, treatActiveFlag);
      break;
    case 6:
      d = resolveEdgeContactBary(nTri,1,cSphere,delta,bary, treatActiveFlag);
      break;
    case 7: // face contact - all three barycentric coordinates are > 0
      d = resolveFaceContactBary(nTri,cSphere,node0ToSphereCenter,delta);
      break;
    default:
      this->error->one(FLERR,"Internal error");
      d = 1.; // doesn't exist, just to satisfy the compiler
      break;
    }
     // return distance - radius of the particle
    return d;
  }

  /* ---------------------------------------------------------------------- */

  inline bool TriMesh::sphereTriangleIntersection(int nTri, double rSphere, double *cSphere)
  {
    double delta[3], bary[3];
    return pointToTriangleDistance(nTri, cSphere, delta, false, bary) < rSphere;
  }

//particle-triangle intersection check (true/false)
  /* ---------------------------------------------------------------------- */
  inline int TriMesh::superquadricTriangleIntersection(int nTri, double *point_of_lowest_potential, Superquadric particle)
  {
    if(isInElement(point_of_lowest_potential,nTri))
      return 0;
    else {
      double **n = node_(nTri);
      if(particle.shape_function_global(n[0]) < 0.0 or
         particle.shape_function_global(n[1]) < 0.0 or
         particle.shape_function_global(n[2]) < 0.0)
        return 1;
      else {
        if(particle.edge_intersection(n[0], n[1]) or
           particle.edge_intersection(n[1], n[2]) or
           particle.edge_intersection(n[2], n[0]) )
          return 2;
        else
          return -1;
      }
    }
  }

#endif

#endif
