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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_I_H
#define LMP_TRI_MESH_I_H

#ifndef SMALL_TRIMESH
#define SMALL_TRIMESH (1.e-10)  
#endif
#define LARGE_TRIMESH 1000000

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSphereContact(int iPart,int nTri, double rSphere, double *cSphere, double *delta)
  {
    // this is the overlap algorithm, neighbor list build is
    // coded in resolveTriSphereNeighbuild

    double bary[3];
    return resolveTriSphereContactBary(iPart,nTri,rSphere,cSphere,delta,bary);
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSphereContactBary(int iPart, int nTri, double rSphere,
                                   double *cSphere, double *delta, double *bary)
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
      d = resolveCornerContactBary(nTri,0,obtuseAngleIndex == 0,cSphere,delta,bary);
      break;
    case 2: 
      d = resolveCornerContactBary(nTri,1,obtuseAngleIndex == 1,cSphere,delta,bary);
      break;
    case 3: 
      d = resolveEdgeContactBary(nTri,0,cSphere,delta,bary);
      break;
    case 4: 
      d = resolveCornerContactBary(nTri,2,obtuseAngleIndex == 2,cSphere,delta,bary);
      break;
    case 5: 
      d = resolveEdgeContactBary(nTri,2,cSphere,delta,bary);
      break;
    case 6: 
      d = resolveEdgeContactBary(nTri,1,cSphere,delta,bary);
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
    return d - rSphere;
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveEdgeContactBary(int iTri, int iEdge, double *p, double *delta, double *bary)
  {
      int ip = (iEdge+1)%3, ipp = (iEdge+2)%3;
      double nodeToP[3], d(1.);
      double **n = node_(iTri);

      vectorSubtract3D(p,n[iEdge],nodeToP);

      double distFromNode =  vectorDot3D(nodeToP,edgeVec(iTri)[iEdge]);

      if(distFromNode < -SMALL_TRIMESH){
        
        if(!cornerActive(iTri)[iEdge])
            return LARGE_TRIMESH;
        d = calcDist(p,n[iEdge],delta);
        bary[iEdge] = 1.; bary[ip] = 0.; bary[ipp] = 0.;
      }
      else if(distFromNode > edgeLen(iTri)[iEdge] + SMALL_TRIMESH){
        
        if(!cornerActive(iTri)[ip])
            return LARGE_TRIMESH;
        d = calcDist(p,n[ip],delta);
        bary[iEdge] = 0.; bary[ip] = 1.; bary[ipp] = 0.;
      }
      else{
        
        double closestPoint[3];

        if(!edgeActive(iTri)[iEdge])
            return LARGE_TRIMESH;

        vectorAddMultiple3D(n[iEdge],distFromNode,edgeVec(iTri)[iEdge],closestPoint);

        d = calcDist(p,closestPoint,delta);

        bary[ipp] = 0.;
        bary[iEdge] = 1. - distFromNode/edgeLen(iTri)[iEdge];
        bary[ip] = 1. - bary[iEdge];
      }

      return d;
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveCornerContactBary(int iTri, int iNode, bool obtuse,
                                                    double *p, double *delta, double *bary)
  {
      int ip = (iNode+1)%3, ipp = (iNode+2)%3;
      //double d(1.);
      double *n = node_(iTri)[iNode];

      if(obtuse){
        
        double **edge = edgeVec(iTri);
        double nodeToP[3], closestPoint[3];

        vectorSubtract3D(p,n,nodeToP);

        double distFromNode = vectorDot3D(nodeToP,edge[ipp]);
        if(distFromNode < SMALL_TRIMESH)
          {
            if(distFromNode > -edgeLen(iTri)[ipp]){
              
              if(!edgeActive(iTri)[ipp])
                return LARGE_TRIMESH;

              vectorAddMultiple3D(n,distFromNode,edge[ipp],closestPoint);

              bary[ip] = 0.;
              bary[iNode] = 1. + distFromNode/edgeLen(iTri)[ipp];
              bary[ipp] = 1. - bary[iNode];

              return calcDist(p,closestPoint,delta);
            } else{
              
              if(!cornerActive(iTri)[ipp])
                return LARGE_TRIMESH;

              bary[ipp] = 1.; bary[iNode] = bary[ip] = 0.;
              return calcDist(p,node_(iTri)[ipp],delta);
            }
          }

        distFromNode = vectorDot3D(nodeToP,edge[iNode]);
        if(distFromNode > -SMALL_TRIMESH)
          {
            if(distFromNode < edgeLen(iTri)[iNode]){
              
              if(!edgeActive(iTri)[iNode])
                return LARGE_TRIMESH;

              vectorAddMultiple3D(n,distFromNode,edge[ipp],closestPoint);

              bary[ipp] = 0.;
              bary[iNode] = 1. - distFromNode/edgeLen(iTri)[iNode];
              bary[ip] = 1. - bary[iNode];

              return calcDist(p,closestPoint,delta);
            } else{
              
              if(!cornerActive(iTri)[ip])
                return LARGE_TRIMESH;

              bary[ip] = 1.; bary[iNode] = bary[ipp] = 0.;
              return calcDist(p,node_(iTri)[ip],delta);

            }
          }
      }

      if(!cornerActive(iTri)[iNode])
          return LARGE_TRIMESH;

      bary[iNode] = 1.; bary[ip] = bary[ipp] = 0.;
      return calcDist(p,node_(iTri)[iNode],delta);
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveFaceContactBary(int iTri, double *p, double *node0ToSphereCenter, double *delta)
  {
      double *surfNorm = SurfaceMeshBase::surfaceNorm(iTri);

      double dNorm = vectorDot3D(surfNorm,node0ToSphereCenter);

      double csPlane[3], tmp[3];
      vectorScalarMult3D(surfNorm,dNorm,tmp);
      vectorSubtract3D(p,tmp,csPlane);

      return calcDist(p,csPlane,delta);
  }

  /* ---------------------------------------------------------------------- */

  inline bool TriMesh::resolveTriSphereNeighbuild(int nTri, double rSphere,
      double *cSphere, double treshold)
  {
    
    double maxDist = rSphere + treshold;

    double dNorm = fabs( calcDistToPlane(cSphere,
        SurfaceMeshBase::center_(nTri),SurfaceMeshBase::surfaceNorm(nTri)) );
    if(dNorm > maxDist) return false;

    double **node = MultiNodeMesh<3>::node_(nTri),
        **edgeNorm = SurfaceMeshBase::edgeNorm(nTri);

    // d_para^2 + d_norm^2 > maxDist^2 --> return false
    double dParaMax = maxDist*maxDist;// - dNorm*dNorm;

    for(int i=0;i<3;i++){
      double d = calcDistToPlane(cSphere,node[i],edgeNorm[i]);
      if(d>0 && d*d > dParaMax)
        return false;
    }

    /*
    for(int i=0;i<3;i++)
      if(calcDist(cSphere,node_[i]) > maxDist)
        return false;
    */
    return true;
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::calcDist(double *cs, double *closestPoint, double *delta)
  {
    vectorSubtract3D(closestPoint,cs,delta);
    return pointDistance(cs,closestPoint);
  }

  inline double TriMesh::calcDistToPlane(double *p, double *pPlane, double *nPlane)
  {
    double v[3];
    vectorSubtract3D(p,pPlane,v);
    // normal distance of sphere center_ to plane
    return vectorDot3D(nPlane,v);
  }

  /* ----------------------------------------------------------------------
   calculate area of a triangle
  ------------------------------------------------------------------------- */

  inline double TriMesh::calcArea(int n)
  {
    double *vecTmp3 = new double[3];

    vectorCross3D(SurfaceMeshBase::edgeVec(n)[0],
                    SurfaceMeshBase::edgeVec(n)[1],vecTmp3);

    // edgevecs are normalized so have to multiply with their lengths
    double area = 0.5*vectorMag3D(vecTmp3) * edgeLen(n)[0] * edgeLen(n)[1];
    delete[] vecTmp3;
    return area;
  }

  /* ----------------------------------------------------------------------
   check if point in triangle, within round-off
   from http://www.blackpawn.com/texts/pointinpoly/default.html
  ------------------------------------------------------------------------- */

  inline bool TriMesh::isInElement(double *pos,int i)
  {
    double v0[3],v1[3],v2[3];
    double dot00,dot01,dot02,dot11,dot12,invDenom,u,v;
    double ***node = node_.begin();

    vectorSubtract3D(node[i][2], node[i][0], v0);
    vectorSubtract3D(node[i][1], node[i][0], v1);
    vectorSubtract3D(pos,        node[i][0], v2);

    dot00 = vectorDot3D(v0, v0);
    dot01 = vectorDot3D(v0, v1);
    dot02 = vectorDot3D(v0, v2);
    dot11 = vectorDot3D(v1, v1);
    dot12 = vectorDot3D(v1, v2);

    invDenom = 1. / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    if((u > -SMALL_TRIMESH) && (v > -SMALL_TRIMESH) && (u + v < 1.+SMALL_TRIMESH))
        return true;
    else
        return false;
  }

  /* ----------------------------------------------------------------------
   generates a random point on the surface that lies within my subbox
  ------------------------------------------------------------------------- */

  inline int TriMesh::generateRandomSubbox(double *pos)
  {
      int index;
      do
      {
          index = generateRandomOwnedGhost(pos);
      }
      while(!domain->is_in_subdomain(pos));

      return index;
  }

  /* ----------------------------------------------------------------------
   generates a random point on the surface that lies on an owned or ghost element
  ------------------------------------------------------------------------- */

  inline int TriMesh::generateRandomOwnedGhost(double *pos)
  {
    double u,v, bary_0,bary_1,bary_2;
    double ***node = node_.begin();
    int nTri = sizeLocal() + sizeGhost();

    // step 1 - choose triangle
    int chosen = randomOwnedGhostElement();
    
    if(chosen >= nTri || chosen < 0)
    {
        
        error->one(FLERR,"TriMesh::generate_random error");
        return -1;
    }

    // step 2 - random bary coords
    
    do {
        u = random_->uniform();
        v = random_->uniform();
    } while( u+v > 1);

    bary_0 = 1. - u - v;
    bary_1 = v;
    bary_2 = u;

    pos[0] = bary_0 * node[chosen][0][0] + bary_1 * node[chosen][1][0] + bary_2 * node[chosen][2][0];
    pos[1] = bary_0 * node[chosen][0][1] + bary_1 * node[chosen][1][1] + bary_2 * node[chosen][2][1];
    pos[2] = bary_0 * node[chosen][0][2] + bary_1 * node[chosen][1][2] + bary_2 * node[chosen][2][2];

    return chosen;
  }

  /* ----------------------------------------------------------------------
   not implemented in thus class
  ------------------------------------------------------------------------- */

  inline int TriMesh::generateRandomOwnedGhostWithin(double *pos,double delta)
  {
      UNUSED(pos);
      UNUSED(delta);
      error->one(FLERR,"internal error");
      return 0;
  }

#endif
