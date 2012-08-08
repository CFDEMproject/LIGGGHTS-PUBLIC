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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_I_H
#define LMP_TRI_MESH_I_H

#define SMALL_TRIMESH 1.e-10

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSphereContact(int nTri, double rSphere, double *cSphere, double *delta)
  {
    // this is the overlap algorithm, neighbor list build is
    // coded in resolveTriSphereNeighbuild

    double tmp[3];

    // sphere-plane distance is coded explicitly because we need an intermediate result
    double triCenterToSphereCenter[3];
    double *surfNorm = SurfaceMesh<3>::surfaceNorm(nTri);
    vectorSubtract3D(cSphere,SurfaceMesh<3>::center_(nTri),triCenterToSphereCenter);

    // normal distance of sphere center to plane
    double dNorm = vectorDot3D(surfNorm,triCenterToSphereCenter);

    if(rSphere > 0. && (dNorm > rSphere || dNorm < -rSphere))
        return 1.;

    // else: go on with algorithm, calc projection of sphere center to plane
    double csPlane[3];
    vectorScalarMult3D(surfNorm,dNorm,tmp);
    vectorSubtract3D(cSphere,tmp,csPlane);

    double nodeToCsPlane[3];
    double **node = MultiNodeMesh<3>::node_(nTri),
       **edgeNorm = SurfaceMesh<3>::edgeNorm(nTri);
    int i;
    double distFromEdge(0.);
    for(i=0;i<3;i++){
      vectorSubtract3D(csPlane,node[i],nodeToCsPlane);
      distFromEdge = vectorDot3D(edgeNorm[i],nodeToCsPlane);
      
      if(distFromEdge > SMALL_TRIMESH) break;
    }

    if(i == 3) // then closest point on triangle is projection on surface
    {
        
        return (calcDist(cSphere,csPlane,delta) - rSphere);
    }

    double *edgeVec = SurfaceMesh<3>::edgeVec(nTri)[i];
    double distFromNode = vectorDot3D(nodeToCsPlane,edgeVec);

    if(distFromNode < 0.)
    {
      if(SurfaceMesh<3>::cornerActive(nTri)[i])
      {
          
          return calcDist(cSphere,node[i],delta) - rSphere;
      }
      else
          return 1.;
    }
    else if(distFromNode > edgeLen(nTri)[i])
    {
      if(SurfaceMesh<3>::cornerActive(nTri)[(i+1)%3])
      {
          
          return calcDist(cSphere,node[(i+1)%3],delta) - rSphere;
      }
      else
          return 1.;
    }

    if(!SurfaceMesh<3>::edgeActive(nTri)[i])
      return 1.;

    double edgeVecTmp[3];
    vectorScalarMult3D(edgeVec,distFromNode,edgeVecTmp);
    vectorAdd3D(node[i],edgeVecTmp,edgeVecTmp); // use edgeVec as contact point

    double d = calcDist(cSphere,edgeVecTmp,delta);
    return d - rSphere;
  }

  /* ---------------------------------------------------------------------- */

  inline double TriMesh::resolveTriSphereContactBary(int nTri, double rSphere,
                                   double *cSphere, double *delta, double *bary)
  {
    bool print = false;//(nTri == 25);
    // this is only the overlap algorithm, neighbor list build is
    // coded in resolveTriSphereNeighbuild

    bary[0] = bary[1] = bary[2] = 0.;
    double tmp[3];
    double **n = node_(nTri);

    // sphere-plane distance is coded explicitly because we need an intermediate result
    double node0ToSphereCenter[3];
    double *surfNorm = SurfaceMesh<3>::surfaceNorm(nTri);
    vectorSubtract3D(cSphere,n[0],node0ToSphereCenter);
    // normal distance of sphere center_ to plane
    double dNorm = vectorDot3D(surfNorm,node0ToSphereCenter);

    if(rSphere > 0. && (dNorm > rSphere || dNorm < -rSphere)) return 1.;

    // else: go on with algorithm, calc projection of sphere center_ to plane
    double csPlane[3];
    vectorScalarMult3D(surfNorm,dNorm,tmp);
    vectorSubtract3D(cSphere,tmp,csPlane);

    double node0ToCsPlane[3];
    vectorSubtract3D(csPlane,n[0],node0ToCsPlane);

    MathExtraLiggghts::calcBaryTriCoords(node0ToCsPlane,edgeVec(nTri),edgeLen(nTri),bary);

    int barySign = (bary[0] > -SMALL_TRIMESH) + 2*(bary[1] > -SMALL_TRIMESH) + 4*(bary[2] > -SMALL_TRIMESH);
/*
    if(print){
      printf("node_ ");
      for(int i=0;i<3;i++)
        printf("%f %f %f | ",n[i][0],n[i][1],n[i][2]);

      printf("%f %f %f\n",csPlane[0],csPlane[1],csPlane[2]);

      printf("iTri %d barySign %d bary %f %f %f \n",nTri,barySign,bary[0],bary[1],bary[2]);
    }
*/
    double d(0.);

    switch(barySign){
    case 1:
    case 2:
    case 3: // bary[2] < 0 --> edge contact on edge[0]
      d = resolveEdgeContact(nTri,0,cSphere,csPlane,delta,bary);
      break;
    case 4:
    case 6: // bary[0] < 0 --> edge contact on edge[1]
      d = resolveEdgeContact(nTri,1,cSphere,csPlane,delta,bary);
      break;
    case 5: // bary[1] < 0 --> edge contact on edge[2]
      d = resolveEdgeContact(nTri,2,cSphere,csPlane,delta,bary);
      break;
    case 7: // face contact - all three barycentric coordinates are > 0
      d = calcDist(cSphere,csPlane,delta);
      break;
    default:
      d = 1.; // doesn't exist, just to satisfy the compiler
      break;
    }

    return d == 1. ? d : d - rSphere;

  }

  /* ---------------------------------------------------------------------- */
    /*
     * p : sphere center_
     * pPlane : projection of p to triangle plane
     */
  inline double TriMesh::resolveEdgeContact(int iTri, int iEdge, double *p, double *pPlane, double *delta, double *bary)
  {
    bool print = false;//(iTri == 25);
    double tmp[3];
    double **n = node_(iTri);
    int ip = (iEdge+1)%3, ipp = (iEdge+2)%3;

    vectorSubtract3D(pPlane,n[iEdge],tmp);
    double d(0.), distFromNode = vectorDot3D(tmp,edgeVec(iTri)[iEdge]);

    if(distFromNode <= 0){
      if(!cornerActive(iTri)[iEdge]) d=1.;
      else{
        bary[iEdge] = 1.; bary[ip] = 0.; bary[ipp] = 0.;
        d = calcDist(p,node_(iTri)[iEdge],delta);
      }
    } else if(distFromNode >= edgeLen(iTri)[iEdge]){
      if(!cornerActive(iTri)[ip]) d=1.;
      else{
        bary[iEdge] = 0.; bary[ip] = 1.; bary[ipp] = 0.;
        d = calcDist(p,node_(iTri)[ip],delta);
      }
    } else{
      if(!edgeActive(iTri)[iEdge]) d=1.;
      else{
        bary[ipp] = 0.;
        bary[iEdge] = 1. - distFromNode/edgeLen(iTri)[iEdge];
        bary[ip] = 1. - bary[iEdge];
        vectorScalarMult3D(edgeVec(iTri)[iEdge],distFromNode,tmp);
        vectorAdd3D(tmp,n[iEdge],tmp);
        d = calcDist(p,tmp,delta);
      }
    }
    /*
        if(print)
          printf("triangle %d bary %f %f %f | distFromEdge %f edgeLen_ %f \n",
              iTri,bary[0],bary[1],bary[2],distFromNode,edgeLen_(iTri)[iEdge]);
    */
    return d;
  }

  /* ---------------------------------------------------------------------- */

  inline bool TriMesh::resolveTriSphereNeighbuild(int nTri, double rSphere,
      double *cSphere, double treshold)
  {
    
    double maxDist = rSphere + treshold;

    double dNorm = fabs( calcDistToPlane(cSphere,
        SurfaceMesh<3>::center_(nTri),SurfaceMesh<3>::surfaceNorm(nTri)) );
    if(dNorm > maxDist) return false;

    double **node = MultiNodeMesh<3>::node_(nTri),
        **edgeNorm = SurfaceMesh<3>::edgeNorm(nTri);

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

    vectorCross3D(SurfaceMesh<3>::edgeVec(n)[0],
                    SurfaceMesh<3>::edgeVec(n)[1],vecTmp3);

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

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    if((u > -SMALL_TRIMESH) && (v > -SMALL_TRIMESH) && (u + v < 1+SMALL_TRIMESH)) return true;
    return true;
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
   generates a random point on the surface that lies within my subbox
   additionally, point is within a distance of delta from active edges
   so can generate a point within a distance delta from boundary on planar meshes
   coarse approch: test only against active edges of chosen triangle and its neigbors
  ------------------------------------------------------------------------- */

  inline int TriMesh::generateRandomSubboxWithin(double *pos,double delta)
  {
      // TODO
      
      error->all(FLERR,"all_in 'yes' not yet implemented");
      return 0;
  }

  /* ----------------------------------------------------------------------
   generates a random point on the surface that lies on an owned or ghost element
  ------------------------------------------------------------------------- */

  inline int TriMesh::generateRandomOwnedGhost(double *pos)
  {
    double u,v, tmp, bary_0,bary_1,bary_2;
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

#endif
