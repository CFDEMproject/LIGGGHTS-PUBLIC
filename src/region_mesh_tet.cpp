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

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <algorithm>
#include "region_mesh_tet.h"
#include "lammps.h"
#include "bounding_box.h"
#include "tri_mesh.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "input_mesh_tet.h"

// include last to ensure correct macros
#include "domain_definitions.h"

#define DELTA_TET 1000

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegTetMesh::RegTetMesh(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg),
  bounding_box_mesh(*new BoundingBox(BIG,-BIG,BIG,-BIG,BIG,-BIG)),
  neighList(*new RegionNeighborList<interpolate_no>(lmp)),
  tri_mesh(*new TriMesh(lmp))
{
  if(narg < 14) error->all(FLERR,"Illegal region mesh/tet command");
  options(narg-14,&arg[14]);

  if(scaleflag) error->all(FLERR,"Lattice scaling not implemented for region mesh/tet, please use 'units box'");

  //TODO parse all_in yes
  //TODO disallow !interior and all_in yes

  if(strcmp(arg[2],"file"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'scale'");
  char *filename = arg[3];

  if(strcmp(arg[4],"scale"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'scale'");
  scale_fact = atof(arg[5]);
  if(strcmp(arg[6],"move"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'move'");
  off_fact[0] = atof(arg[7]);
  off_fact[1] = atof(arg[8]);
  off_fact[2] = atof(arg[9]);
  if(strcmp(arg[10],"rotate"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'rotate'");
  rot_angle[0] = atof(arg[11]);
  rot_angle[1] = atof(arg[12]);
  rot_angle[2] = atof(arg[13]);

  node = NULL;
  center = NULL;
  rbound = NULL;
  rbound_max = 0.;

  n_face_neighs = NULL;
  face_neighs = NULL;
  n_face_neighs_node = NULL;

  n_node_neighs = NULL;
  node_neighs = NULL;

  n_surfaces = NULL;
  surfaces = NULL;

  volume = NULL;
  acc_volume = NULL;
  nTet = 0;
  nTetMax = 0;
  total_volume = 0.;

  // manage input
  InputMeshTet *my_input = new InputMeshTet(lmp, 0, NULL);
  my_input->meshtetfile(filename,this,true);
  delete my_input;

  // extent of region and mesh

  set_extent_mesh();

  if (interior) {
    bboxflag = 1;
    set_extent_region();
    tri_mesh.useAsInsertionMesh(false);
    build_neighs();
    build_surface();
    tri_mesh.initialSetup();
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegTetMesh::~RegTetMesh()
{
  delete [] contact;

  delete &bounding_box_mesh;
  delete &neighList;
  delete &tri_mesh;

  memory->destroy(node);
  memory->destroy(center);
  memory->sfree(volume);
  memory->sfree(acc_volume);
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegTetMesh::inside(double x, double y, double z)
{
   double pos[3];
   pos[0] = x; pos[1] = y; pos[2] = z;

   // check subdomain
   if(!domain->is_in_subdomain(pos)) return 0;

   // check bbox, only if exists
   if(bboxflag)
   {
       if(pos[0] < extent_xlo || pos[0] > extent_xhi) return 0;
       if(pos[1] < extent_ylo || pos[1] > extent_yhi) return 0;
       if(pos[2] < extent_zlo || pos[2] > extent_zhi) return 0;
   }

   // brute force naive search
   int inside_mesh = 0;
   for(int i = 0; i < nTet; i++)
   {
        inside_mesh = inside_mesh + is_inside_tet(i,pos);
        if(inside_mesh) break;
   }

   //fprintf(screen,"checking pos %f %f %f, result %d; ntet %d\n",x,y,z,inside_mesh,nTet);

   return inside_mesh;
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::surface_interior(double *x, double cutoff)
{
  error->one(FLERR,"This feature is not available for tet mesh regions");
  return 0;
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::surface_exterior(double *x, double cutoff)
{
  error->one(FLERR,"This feature is not available for tet mesh regions");
  return 0;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::generate_random(double *pos,bool subdomain_flag)
{
    if(!interior) error->all(FLERR,"Impossible to generate random points on tet mesh region with side = out");

    int ntry = 0;

    do
    {
        mesh_randpos(pos);
        ntry++;
    }
    while(ntry < 10000 && (!subdomain_flag || !domain->is_in_subdomain(pos)));

    if(10000 == ntry)
        error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

// generates a random point within the region and has a min distance from surface
// i.e. generate random point in region "shrunk" by cut
void RegTetMesh::generate_random_shrinkby_cut(double *pos,double cut,bool subdomain_flag)
{
    int ntry = 0;
    bool is_near_surface = false;
    int barysign = -1;

    for(int i = 0; i < nTet; i++)
    {
        
    }

    do
    {
       ntry++;
       
       int iTetChosen = mesh_randpos(pos);
       is_near_surface = false;
       double delta[3];

       // check all surfaces of chosen tet
       for(int is = 0; is < n_surfaces[iTetChosen]; is++)
       {
         int iSurf = surfaces[iTetChosen][is];
         
         if(tri_mesh.resolveTriSphereContact(-1,iSurf,cut,pos,delta,barysign) < 0)
         {
            is_near_surface = true;
            break; 
          }
          
       }

       if(is_near_surface)
        continue; 

       // check all surfaces of face neigh tets
       for(int iFaceNeigh = 0; iFaceNeigh < n_face_neighs[iTetChosen]; iFaceNeigh++)
       {
           for(int is = 0; is < n_surfaces[face_neighs[iTetChosen][iFaceNeigh]]; is++)
           {
             int iSurf = surfaces[face_neighs[iTetChosen][iFaceNeigh]][is];
             
             if(tri_mesh.resolveTriSphereContact(-1,iSurf,cut,pos,delta,barysign) < 0)
             {
                is_near_surface = true;
                break; 
              }
              
           }
       }
       if(is_near_surface)
        continue; 

       // check all surfaces of node neigh tets
       for(int iNodeNeigh = 0; iNodeNeigh < n_node_neighs[iTetChosen]; iNodeNeigh++)
       {
           
           for(int is = 0; is < n_surfaces[node_neighs[iTetChosen][iNodeNeigh]]; is++)
           {
             int iSurf = surfaces[node_neighs[iTetChosen][iNodeNeigh]][is];
                
             if(tri_mesh.resolveTriSphereContact(-1,iSurf,cut,pos,delta,barysign) < 0)
             {
                is_near_surface = true;
                break; 
             }
             
           }
       }
    }
    // pos has to be within region, and within cut of region surface
    while(ntry < 10000 && (is_near_surface || (!subdomain_flag || !domain->is_in_subdomain(pos))));

    if(10000 == ntry)
    {
        error->one(FLERR,"internal error");
    }
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::add_tet(double **n)
{
    double ctr[3];

    if(nTet == nTetMax) grow_arrays();

    vectorZeroize3D(ctr);
    for(int i=0;i<4;i++)
    {
        vectorCopy3D(n[i],node[nTet][i]);
        vectorAdd3D(ctr,node[nTet][i],ctr);
    }
    vectorScalarDiv3D(ctr,4.);
    vectorCopy3D(ctr,center[nTet]);

    double vol = volume_of_tet(nTet);
    if(vol < 0.)
    {
        // flip nodes 0 and 3
        double node0[3];
        vectorCopy3D(node[nTet][0],node0);
        vectorCopy3D(node[nTet][3],node[nTet][0]);
        vectorCopy3D(node0,node[nTet][3]);
    }

    vol = volume_of_tet(nTet);
    if(vol < 0.) error->all(FLERR,"Fatal error: RegTetMesh::add_tet: vol < 0");

    volume[nTet] = vol;
    total_volume += volume[nTet];
    acc_volume[nTet] = volume[nTet];
    if(nTet > 0) acc_volume[nTet] += acc_volume[nTet-1];
    nTet++;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::build_neighs()
{
    neighList.reset();

    for(int i = 0; i < nTet; i++)
    {
        vectorZeroize3D(center[i]);

        for(int j=0;j<4;j++)
            vectorAdd3D(node[i][j],center[i],center[i]);
        vectorScalarDiv3D(center[i],4.);

        double rb = 0., vec[3];
        for(int j = 0; j < 4; j++)
        {
            vectorSubtract3D(center[i],node[i][j],vec);
            rb = std::max(rb,vectorMag3D(vec));
        }
        rbound[i] = rb;
        if(rb > rbound_max)
            rbound_max = rb;

    }

    neighList.setBoundingBox(bounding_box_mesh, rbound_max);

    for(int i = 0; i < nTet; i++)
    {
        n_face_neighs[i] = 0;
        n_node_neighs[i] = 0;
        vectorZeroizeN(n_face_neighs_node[i],4);
    }

    bool badMesh = false;
    for(int i = 0; i < nTet; i++)
    {
        std::vector<int> overlaps;
        neighList.hasOverlapWith(center[i], rbound[i],overlaps);

        for(size_t icontainer = 0; icontainer < overlaps.size(); icontainer++)
        {
            int iOverlap = overlaps[icontainer];

            if(iOverlap < 0 || iOverlap >= nTet)
               this->error->one(FLERR,"Mesh error: internal error");

            int nNodesEqual = 0;
            std::vector<int> iNodesInvolved, iOverlapNodesInvolved;

            for(int iNode = 0; iNode < 4; iNode++)
            {
                for(int iOverlapNode = 0; iOverlapNode < 4; iOverlapNode++)
                {
                    if(nodesAreEqual(node[i][iNode],node[iOverlap][iOverlapNode],1.e-12))
                    {
                        iNodesInvolved.push_back(iNode);
                        iOverlapNodesInvolved.push_back(iOverlapNode);
                        nNodesEqual++;
                    }
                }
            }

            if(0 == nNodesEqual)
                continue;
            
            else if(3 > nNodesEqual)
            {
                
                if(100 == n_node_neighs[i])
                    badMesh = true;
                else
                {
                    node_neighs[i][n_node_neighs[i]] = iOverlap;
                    n_node_neighs[i]++;
                }

                if(100 == n_node_neighs[iOverlap])
                    badMesh = true;
                else
                {
                    node_neighs[iOverlap][n_node_neighs[iOverlap]] = i;
                    n_node_neighs[iOverlap]++;
                }
            }
            
            else if(3 == nNodesEqual)
            {
                
                face_neighs[i][n_face_neighs[i]++] = iOverlap;
                face_neighs[iOverlap][n_face_neighs[iOverlap]++] = i;

                n_face_neighs_node[i][iNodesInvolved[0]]++;
                n_face_neighs_node[i][iNodesInvolved[1]]++;
                n_face_neighs_node[i][iNodesInvolved[2]]++;
                n_face_neighs_node[iOverlap][iOverlapNodesInvolved[0]]++;
                n_face_neighs_node[iOverlap][iOverlapNodesInvolved[1]]++;
                n_face_neighs_node[iOverlap][iOverlapNodesInvolved[2]]++;

                //TODO: build up the tri mesh here, link the two meshes
                //TODO: worst case: broad phase (wie auf zettel skizziert)
            }
            else
            {
                
                char errstr[256];

                sprintf(errstr,"duplicate elements %d and %d in tet for region %s",i,iOverlap,id);
                error->one(FLERR,errstr);
            }
        }

        neighList.insert(center[i], rbound[i],i);
    }
    if (badMesh)
        error->warningAll(FLERR,"Region mesh/tet: too many node neighbors, mesh is of bad quality; 'all_in yes' might not work correctly");

    for(int i = 0; i < nTet; i++)
    {
        
    }
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::build_surface()
{
    
    double **nodeTmp = create<double>(nodeTmp,3,3);

    int n_surf_elems = 0;
    for(int i = 0; i < nTet; i++)
    {
        n_surfaces[i] = 0;

        int dummy;

        // 4 neighs
        if(3 == n_face_neighs_node[i][0] && 3 == n_face_neighs_node[i][1] &&
           3 == n_face_neighs_node[i][2] && 3 == n_face_neighs_node[i][3])
        {
            if(! (4 == n_face_neighs[i]))
                error->one(FLERR,"assertion failed");
           continue;
        }

        // 3 neighs
        if(3 == n_face_neighs_node[i][0] || 3 == n_face_neighs_node[i][1] ||
           3 == n_face_neighs_node[i][2] || 3 == n_face_neighs_node[i][3])
        {
            if(!(3 == n_face_neighs[i]))
                error->one(FLERR,"assertion failed");
            
            int which;
            MathExtraLiggghts::max(n_face_neighs_node[i],4,which);
            for(int k=0;k<3;k++){
              nodeTmp[0][k] = node[i][(which+1)%4][k];
              nodeTmp[1][k] = node[i][(which+2)%4][k];
              nodeTmp[2][k] = node[i][(which+3)%4][k];
            }
            tri_mesh.addElement(nodeTmp,-1);
            surfaces[i][n_surfaces[i]] = n_surf_elems;
            n_surfaces[i]++;
            n_surf_elems++;
        }
        // 0 neighs
        else if(0 == MathExtraLiggghts::max(n_face_neighs_node[i],4,dummy))
        {
            if(!(0 == n_face_neighs[i]))
                error->one(FLERR,"assertion failed");

            // add 0 1 2 ,  1 2 3,  2 3 0,  3 0 1
            for(int iAdd = 0; iAdd < 4; iAdd++)
            {
                for(int k=0;k<3;k++){
                  nodeTmp[0][k] = node[i][(iAdd+0)%4][k];
                  nodeTmp[1][k] = node[i][(iAdd+1)%4][k];
                  nodeTmp[2][k] = node[i][(iAdd+2)%4][k];
                }
                tri_mesh.addElement(nodeTmp,-1);
                surfaces[i][n_surfaces[i]] = n_surf_elems;
                n_surfaces[i]++;
                n_surf_elems++;
            }
        }
        // 2 neighs
        else if( (2 == MathExtraLiggghts::max(n_face_neighs_node[i],4,dummy)) &&
                 (1 == MathExtraLiggghts::min(n_face_neighs_node[i],4,dummy)) &&
                 (6 == vectorSumN(n_face_neighs_node[i],4))
                )
        {
            if(! (2 == n_face_neighs[i]))
                error->one(FLERR,"assertion failed");

            int which_hi_1;
            MathExtraLiggghts::max(n_face_neighs_node[i],4,which_hi_1);
            int which_lo_1;
            MathExtraLiggghts::min(n_face_neighs_node[i],4,which_lo_1);
            int which_2 = (which_hi_1+1)%4;
            if(which_2 == which_lo_1)
                which_2 = (which_lo_1+1)%4;
            int which_3 = 6 - (which_hi_1+which_lo_1+which_2);

            int which_lo_2,which_hi_2;
            if(n_face_neighs_node[i][which_3] > n_face_neighs_node[i][which_2])
            {
                which_hi_2 = which_3;
                which_lo_2 = which_2;
            }
            else
            {
                which_hi_2 = which_2;
                which_lo_2 = which_3;
            }

            // add lo1 lo2 hi1   and lo1 lo2 hi2

            for(int k=0;k<3;k++){
              nodeTmp[0][k] = node[i][which_lo_1][k];
              nodeTmp[1][k] = node[i][which_lo_2][k];
              nodeTmp[2][k] = node[i][which_hi_1][k];
            }
            tri_mesh.addElement(nodeTmp,-1);
            surfaces[i][n_surfaces[i]] = n_surf_elems;
            n_surfaces[i]++;
            n_surf_elems++;
            for(int k=0;k<3;k++){
              nodeTmp[0][k] = node[i][which_lo_1][k];
              nodeTmp[1][k] = node[i][which_lo_2][k];
              nodeTmp[2][k] = node[i][which_hi_2][k];
            }
            tri_mesh.addElement(nodeTmp,-1);
            surfaces[i][n_surfaces[i]] = n_surf_elems;
            n_surfaces[i]++;
            n_surf_elems++;
        }
        // 1 neighs
        else if(0 == MathExtraLiggghts::min(n_face_neighs_node[i],4,dummy))
        {
            if(! (1 == n_face_neighs[i]))
                error->one(FLERR,"assertion failed");

            int which_lo = dummy;

            for(int iAdd = 0; iAdd < 3; iAdd++)
            {
                for(int k=0;k<3;k++){
                  nodeTmp[0][k] = node[i][(4+which_lo-2+iAdd)%4][k];
                  nodeTmp[1][k] = node[i][(4+which_lo-1+iAdd)%4][k];
                  nodeTmp[2][k] = node[i][(4+which_lo+0+iAdd)%4][k];
                }
                tri_mesh.addElement(nodeTmp,-1);
                surfaces[i][n_surfaces[i]] = n_surf_elems;
                n_surfaces[i]++;
                n_surf_elems++;
            }

        }
        else
            error->one(FLERR,"unrecognized");
    }
    destroy<double>(nodeTmp);
}

/* ---------------------------------------------------------------------- */

bool RegTetMesh::nodesAreEqual(double *nodeToCheck1,double *nodeToCheck2,double precision)
{
    for(int i=0;i<3;i++)
      if(!MathExtraLiggghts::compDouble(nodeToCheck1[i],nodeToCheck2[i],precision))
        return false;
    return true;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::grow_arrays()
{
    nTetMax += DELTA_TET;
    node = (double***)(memory->grow(node,nTetMax, 4 , 3, "vtk_tet_node"));
    center = (double**)(memory->grow(center,nTetMax, 3, "vtk_tet_center"));
    rbound = (double*)(memory->grow(rbound,nTetMax, "vtk_tet_rbound"));

    n_face_neighs = (int*)(memory->grow(n_face_neighs,nTetMax, "vtk_tet_n_face_neighs"));
    face_neighs = (int**)(memory->grow(face_neighs,nTetMax,4,"vtk_tet_face_neighs"));
    n_face_neighs_node = (int**)(memory->grow(n_face_neighs_node,nTetMax,4, "vtk_tet_n_face_neighs_node"));

    n_node_neighs = (int*)(memory->grow(n_node_neighs,nTetMax, "vtk_tet_n_node_neighs"));
    node_neighs = (int**)(memory->grow(node_neighs,nTetMax,100,"vtk_tet_node_neighs"));

    n_surfaces = (int*)(memory->grow(n_surfaces,nTetMax, "vtk_tet_n_surfaces"));
    surfaces = (int**)(memory->grow(surfaces,nTetMax,4,"vtk_tet_surfaces"));

    volume = (double*)(memory->srealloc(volume,nTetMax*sizeof(double),"vtk_tet_volume"));
    acc_volume = (double*)(memory->srealloc(acc_volume,nTetMax*sizeof(double),"vtk_tet_acc_volume"));
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::n_tet()
{
    return nTet;
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::total_vol()
{
    return total_volume;
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::tet_vol(int i)
{
    return volume[i];
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::tet_acc_vol(int i)
{
    return acc_volume[i];
}

/* ---------------------------------------------------------------------- */

inline double RegTetMesh::volume_of_tet(int iTet)
{
    return volume_of_tet(node[iTet][0],node[iTet][1],node[iTet][2],node[iTet][3]);
}

/* ---------------------------------------------------------------------- */

inline int RegTetMesh::is_inside_tet(int iTet,double *pos)
{
    double vol1,vol2,vol3,vol4;

    vol1 = volume_of_tet(node[iTet][0], node[iTet][1], node[iTet][2], pos          );
    vol2 = volume_of_tet(node[iTet][0], node[iTet][1], pos,           node[iTet][3]);
    vol3 = volume_of_tet(node[iTet][0], pos,           node[iTet][2], node[iTet][3]);
    vol4 = volume_of_tet(pos          , node[iTet][1], node[iTet][2], node[iTet][3]);

    if(vol1 > 0. && vol2 > 0. && vol3 > 0. && vol4 > 0.) return 1;
    return 0;
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::volume_of_tet(double* v0, double* v1, double* v2, double* v3)
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

/* ---------------------------------------------------------------------- */

inline void RegTetMesh::set_extent_mesh()
{
    for(int i = 0; i < nTet; i++)
        for(int j=0;j<4;j++)
            bounding_box_mesh.extendToContain(node[i][j]);
}

/* ---------------------------------------------------------------------- */

inline void RegTetMesh::set_extent_region()
{
    extent_xlo = extent_ylo = extent_zlo =  BIG;
    extent_xhi = extent_yhi = extent_zhi = -BIG;

    for(int i = 0; i < nTet; i++)
    {
        for(int j=0;j<4;j++)
        {
            if(node[i][j][0] < extent_xlo) extent_xlo = node[i][j][0];
            if(node[i][j][1] < extent_ylo) extent_ylo = node[i][j][1];
            if(node[i][j][2] < extent_zlo) extent_zlo = node[i][j][2];

            if(node[i][j][0] > extent_xhi) extent_xhi = node[i][j][0];
            if(node[i][j][1] > extent_yhi) extent_yhi = node[i][j][1];
            if(node[i][j][2] > extent_zhi) extent_zhi = node[i][j][2];
        }
    }
}

/* ---------------------------------------------------------------------- */

inline int RegTetMesh::mesh_randpos(double *pos)
{
    int iTriChosen = tet_rand_tri();
    tet_randpos(iTriChosen,pos);
    if(pos[0] == 0. && pos[1] == 0. && pos[2] == 0.)
        error->one(FLERR,"illegal RegTetMesh::mesh_randpos");
    
    return iTriChosen;
}

/* ---------------------------------------------------------------------- */

inline int RegTetMesh::tet_rand_tri()
{
    //SIMPLISTIC
    /*
    double rd = total_volume * random->uniform();
    int chosen = 0;
    while (rd > acc_volume[chosen] && chosen < nTet-1) chosen++;
    return chosen;
    */

    // FAST

    double rd = total_volume * random->uniform();

    int i = nTet/2;
    int imin = 0;
    int imax = nTet -1;
    int ntry = 0;

    // binary search
    do
    {
        if(rd < acc_volume[i] && ((i == 0) ? true : (rd > acc_volume[i-1])))
            return i;

        if(imax == imin)
            error->one(FLERR,"internal error");

        // must go up
        if(rd > acc_volume[i])
        {
            imin = i;
            i = (imax+imin) / 2;
            if (i == imin) i++;
        }
        // must go down
        else if (rd < acc_volume[i-1])
        {
            imax = i;
            i = (imax+imin) / 2;
            if (i == imin) i++;
            if (i == imax && i > 0) i--;
        }
        else
            error->one(FLERR,"internal error");

        ntry++;
    }
    while(ntry < 10000);

    error->one(FLERR,"internal error");
    return 0;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::volume_mc(int n_test,bool cutflag,double cut,double &vol_global,double &vol_local)
{
    double pos[3], volume_in_local = 0., vol_in_local_all;

    //error->all(FLERR,"end");

    //NO TODO: implementation for cutflag = true
    
    if(total_volume == 0.) error->all(FLERR,"mesh/tet region has zero volume, cannot continue");

    vol_global = total_volume; 

    for(int iTet = 0; iTet < nTet; iTet++)
    {
        
        for(int iNode = 0; iNode < 5; iNode++)
        {
            
            double weight = (iNode<4) ? (volume[iTet]*0.1) : (volume[iTet]*0.6);

            if(iNode<4)
                vectorCopy3D(node[iTet][iNode],pos);
            else
                vectorCopy3D(center[iTet],pos);

            if(!domain->is_in_domain(pos))
                error->one(FLERR,"mesh point outside simulation domain");

            // check if point is in subdomain
            if(domain->is_in_subdomain(pos))
                volume_in_local += weight;
        }
    }

    MPI_Sum_Scalar(volume_in_local,vol_in_local_all,world);
    if(vol_in_local_all < 1e-13)
        error->all(FLERR,"Unable to calculate region volume - are you operating on a 2d region?");

    // return calculated values
    vol_local  = volume_in_local;

    // sum of local volumes will not be equal to global volume because of
    // different random generator states - correct this now
    vol_local *= (vol_global/vol_in_local_all);

}
