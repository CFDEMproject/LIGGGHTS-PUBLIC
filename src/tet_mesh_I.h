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

#ifndef LMP_TET_MESH_I_H
#define LMP_TET_MESH_I_H

  /* ----------------------------------------------------------------------
   calculate volume of a tet
  ------------------------------------------------------------------------- */

  inline double TetMesh::calcVol(int n)
  {
     return calcTetVol(node_(n)[0],node_(n)[1],node_(n)[2],node_(n)[3]);
  }

  inline double TetMesh::calcTetVol(double* v0, double* v1, double* v2, double* v3)
  {
     double A[3];
     A[0] = v3[0] - v1[0];
     A[1] = v3[1] - v1[1];
     A[2] = v3[2] - v1[2];

     double B[3];
     B[0] = v2[0] - v1[0];
     B[1] = v2[1] - v1[1];
     B[2] = v2[2] - v1[2];

     double C[3];
     C[0] = v0[0] - v1[0];
     C[1] = v0[1] - v1[1];
     C[2] = v0[2] - v1[2];

     // cross product A x B
     double cp[] = {
        A[1]*B[2] - A[2]*B[1],
       -A[0]*B[2] + A[2]*B[0],
        A[0]*B[1] - A[1]*B[0]
     };

     // dot with C
     double volume = cp[0] * C[0] + cp[1] * C[1] + cp[2] * C[2];
     volume /= 6.;
     return volume;
  }

  /* ----------------------------------------------------------------------
   check if point if inside tet
  ------------------------------------------------------------------------- */

  inline bool TetMesh::isInside(int iTet,double *pos)
  {
      double ***node = node_.begin();
      double vol1,vol2,vol3,vol4;

      vol1 = calcTetVol(node[iTet][0], node[iTet][1], node[iTet][2], pos          );
      if(vol1 < 0.) return false;
      vol2 = calcTetVol(node[iTet][0], node[iTet][1], pos,           node[iTet][3]);
      if(vol2 < 0.) return false;
      vol3 = calcTetVol(node[iTet][0], pos,           node[iTet][2], node[iTet][3]);
      if(vol3 < 0.) return false;
      vol4 = calcTetVol(pos          , node[iTet][1], node[iTet][2], node[iTet][3]);
      if(vol4 < 0.) return false;

      return true;
  }

  /* ----------------------------------------------------------------------
   generates a random point on the surface that lies within my subbox
  ------------------------------------------------------------------------- */

  inline int TetMesh::generateRandomSubbox(double *pos)
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

  inline int TetMesh::generateRandomSubboxWithin(double *pos,double delta)
  {
      // TODO
      error->all(FLERR,"all_in 'yes' not yet implemented");
      return 0;
  }

  /* ----------------------------------------------------------------------
   generates a random point on the surface that lies on an owned or ghost element
  ------------------------------------------------------------------------- */

  inline int TetMesh::generateRandomOwnedGhost(double *pos)
  {
    double s,t,u,tmp,bary_coo[4];
    int nTet = sizeLocal() + sizeGhost();

    // step 1 - choose triangle
    int chosen = randomOwnedGhostElement();
    
    if(chosen >= nTet || chosen < 0)
    {
        
        error->one(FLERR,"TriMesh::generate_random error");
        return -1;
    }

    s = random_->uniform();
    t = random_->uniform();
    u = random_->uniform();

    if(s+t > 1.)
    {
        s = 1.-s;
        t = 1.-t;
    }
    if(t+u > 1.)
    {
        tmp = u;
        u = 1.-s-t;
        t = 1.-tmp;
    }
    else if(s+t+u > 1.)
    {
        tmp = u;
        u = s+t+u-1.;
        s = 1.-t-tmp;
    }

    bary_coo[0] = 1.-s-t-u;
    bary_coo[1] = s;
    bary_coo[2] = t;
    bary_coo[3] = u;

    baryToCart(chosen,bary_coo,pos);

    return chosen;
  }

  inline void TetMesh::baryToCart(int iTet,double *bary_coo,double *pos)
  {
    double ***node = node_.begin();
    for(int i=0;i<3;i++)
       pos[i] = bary_coo[0] * node[iTet][0][i] +
                bary_coo[1] * node[iTet][1][i] +
                bary_coo[2] * node[iTet][2][i] +
                bary_coo[3] * node[iTet][3][i];
  }

  /* ----------------------------------------------------------------------
   check if two elements share a face
  ------------------------------------------------------------------------- */

  bool TetMesh::shareFace(int i, int j, int &iFace, int &jFace)
  {

      //TODO: detect faces of volumes, facenormals etc
      this->error->all(FLERR,"END");
      return false;
  }

#endif
