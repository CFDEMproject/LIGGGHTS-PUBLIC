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

#ifdef LAMMPS_VTK

#ifdef DUMP_CLASS

DumpStyle(mesh/vtk,DumpMeshVTK)
DumpStyle(mesh/gran/VTK,DumpMeshVTK) 

#else

#ifndef LMP_DUMP_MESH_VTK_H
#define LMP_DUMP_MESH_VTK_H

#include "dump.h"
#include "dump_vtk.h"
#include "dump_mesh.h"
#include <vector>

namespace LAMMPS_NS {

class DumpMeshVTK : public Dump, public DumpVTK {

public:

    DumpMeshVTK(LAMMPS *, int, char**);
    virtual ~DumpMeshVTK();
    void init_style();

private:
    int modify_param(int, char **);
    void write_header(bigint ndump) {}
    int count();
    void pack(int *) {}

    void setFileCurrent();
    void write();
    void write_data(int, double *);

    char *filecurrent;
    DumpMesh *dumpMesh_;

    int vtk_file_format_;
    vtkSmartPointer<vtkMultiBlockDataSet> mbSet_;

    int dataMode_;
};

}

#endif
#endif
#endif //LAMMPS_VTK
