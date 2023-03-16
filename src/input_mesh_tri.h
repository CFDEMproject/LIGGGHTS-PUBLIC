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

    Tóth János (MATE, Gödöllő)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_INPUT_MESH_TRI_H
#define LMP_INPUT_MESH_TRI_H

#include <stdio.h>
#include "input.h"
#include <fstream>

namespace LAMMPS_NS {

class InputMeshTri : protected Input
{
  public:
    enum GeneratedType {
        UNKNOWN = 0,
        CUBE,
        BOX,
        CYLINDER,
        PIPE,
        DISK,
        PLANE,
    };

    enum BoxMask {
        BOX_TOP = (1 << 0),
        BOX_BOTTOM = (1 << 1),
        BOX_FRONT = (1 << 2),
        BOX_BACK = (1 << 3),
        BOX_LEFT = (1 << 4),
        BOX_RIGHT = (1 << 5),
    };

    enum CylinderMask {
        CYL_TOP = (1 << 0),
        CYL_BOTTOM = (1 << 1),
        CYL_SIDE = (1 << 2),
    };

    enum PipeMask {
        PIPE_INNER_TOP = (1 << 0),
        PIPE_OUTER_TOP = (1 << 1),
        PIPE_INNER_BOTTOM = (1 << 2),
        PIPE_OUTER_BOTTOM = (1 << 3),
        PIPE_INNER_SIDE = (1 << 4),
        PIPE_OUTER_SIDE = (1 << 5),
    };

    struct GeneratorParameters {
        GeneratedType type;
        int mask;
        // cube/box: xsize, ysize, zsize
        // cylinder: radius, zsize
        // pipe: inner_radius, outer_radius, zsize
        // disk: radius
        // plane: size
        double dvalues[3];
        // cylinder/pipe/disk: number_of_segments
        int ivalues[1];
    };


    InputMeshTri(class LAMMPS *lmp, int narg, char **arg);
    ~InputMeshTri();

    void meshtrifile(const char *filename,class TriMesh *mesh,bool verbose,
                     const int size_exclusion_list, int *exclusion_list,
                     class Region *region);

    void meshgenerator(const GeneratorParameters * params, class TriMesh *mesh,bool verbose,
                     const int size_exclusion_list, int *exclusion_list,
                     class Region *region);

  private:

    void generate_box(const GeneratorParameters * params, class TriMesh *mesh, class Region *region);
    void generate_cylinder(const GeneratorParameters * params, class TriMesh *mesh, class Region *region);
    void generate_pipe(const GeneratorParameters * params, class TriMesh *mesh, class Region *region);
    void generate_disk(const GeneratorParameters * params, class TriMesh *mesh, class Region *region);
    void generate_plane(const GeneratorParameters * params, class TriMesh *mesh, class Region *region);

    void broadcast_and_add_triangle(class TriMesh *mesh, class Region *region, unsigned int *count, double v1x, double v1y, double v1z, double v2x, double  v2y, double v2z, double v3x, double v3y, double v3z);

    bool verbose_;
    int i_exclusion_list_;
    int size_exclusion_list_;
    int *exclusion_list_;

    void meshtrifile_vtk(class TriMesh *mesh,class Region *region);
    void meshtrifile_stl(class TriMesh *mesh,class Region *region, const char * filename);
    void meshtrifile_stl_binary(class TriMesh *, class Region *region, const char * filename);
    inline void addTriangle(class TriMesh *mesh, double *a, double *b, double *c,int lineNumber);

};

}

#endif
