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
    Richard Berger (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(neighlist/mesh,FixNeighlistMesh)

#else

#ifndef LMP_FIX_NEIGHLIST_MESH_H
#define LMP_FIX_NEIGHLIST_MESH_H

#include "fix.h"
#include "container.h"
#include <vector>
#include <algorithm>

namespace LAMMPS_NS
{

struct BinBoundary {
  int xlo;
  int xhi;
  int ylo;
  int yhi;
  int zlo;
  int zhi;
};

struct TriangleNeighlist {
  std::vector<int> contacts;
  std::vector<int> bins;
  BinBoundary boundary;
  int nchecked;

  TriangleNeighlist() : nchecked(0) {}
};

class FixNeighlistMesh : public Fix
{
  public:

    FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg);
    virtual
    ~FixNeighlistMesh();
    virtual int setmask();
    virtual void post_create();
    virtual void pre_delete(bool unfixflag);

    virtual void initializeNeighlist();
    virtual void setup_pre_force(int);
    virtual void min_setup_pre_force(int);

    virtual void pre_neighbor(); 
    virtual void pre_force(int vflag); 
    virtual void min_pre_force(int vflag); 

    virtual void post_run();

    const std::vector<int> & get_contact_list(int iTri) const {
      return triangles[iTri].contacts;
    }

    virtual int getSizeNumContacts();

    void enableTotalNumContacts(bool enable)
    {
      globalNumAllContacts_ = enable;
    }

    int getTotalNumContacts() { return numAllContacts_; }

    bool contactInList(int iTri, int iAtom)
    {
      std::vector<int> & neighbors = triangles[iTri].contacts;
      return std::find(neighbors.begin(), neighbors.end(), iAtom) != neighbors.end();
    }

    inline class FixPropertyAtom* fix_nneighs()
    { return fix_nneighs_; };

    // groupbit merged from groupbit of this fix and fix wall/gran (if exists)
    int groupbit_wall_mesh;

  protected:

    void handleTriangle(int iTri);
    void getBinBoundariesFromBoundingBox(class BoundingBox &b, int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax);
    void getBinBoundariesForTriangle(int iTri, int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax);

    class FixMeshSurface *caller_;
    class TriMesh *mesh_;

    class FixPropertyAtom *fix_nneighs_;
    char *fix_nneighs_name_;

    bool buildNeighList;

    std::vector<TriangleNeighlist> triangles;

    int numAllContacts_;
    bool globalNumAllContacts_;

    int mbinx,mbiny,mbinz,maxhead, *bins, *binhead;
    double skin;

    // max distance from triangle to particle center
    double distmax;

    double **x, *r;

    bool changingMesh;
    bool changingDomain;

    bigint last_bin_update;

    void generate_bin_list(size_t nall);

    class AtomVecEllipsoid *avec;

    bool otherList_;
};

} /* namespace LAMMPS_NS */

#endif /* FIX_MESH_NEIGHLIST_H_ */
#endif /* FIX_CLASS */
