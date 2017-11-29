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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) 

#ifndef LMP_DUMP_MESH_H
#define LMP_DUMP_MESH_H

#include "pointers.h"
#include "container.h"
#include "tri_mesh.h"
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMPIController.h>

namespace LAMMPS_NS
{

class DumpMesh : public Pointers
{

  public:

    DumpMesh(LAMMPS *, int _nclusterprocs, int _multiproc, int _filewriter, int _fileproc, vtkMPIController *controller);
    virtual ~DumpMesh();
    int parse_parameters(const int narg, const char *const *const arg, std::list<std::string> keyword_list = std::list<std::string>());
    int init_style();
    int count();
    void prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet);

    bigint memory_usage()
    { return 0; }

  private:            // column labels

    std::list<TriMesh*> meshList_;
    int dump_what_;
    vtkSmartPointer<vtkMultiBlockDataSet> mbSet;
    vtkMPIController * localController;

    int n_calls_;
    int nclusterprocs;
    int multiproc;
    int filewriter;
    int fileproc;

    // buffer for data from all procs
    int n_all_, n_all_max_;
    double *buf_all_;

    // properties to be dumped
    // TODO: could make look-up more generic

    // stress
    class ScalarContainer<double> **sigma_n_, **sigma_t_;
    // wear
    class ScalarContainer<double> **wear_;
    // vel
    class MultiVectorContainer<double,3,3> **v_node_;
    // stresscomponents
    class VectorContainer<double,3> **f_node_;
    // temp
    class ScalarContainer<double> **T_;
    std::vector<bool> temp_per_element_;
    // min dist from active edge
    class ScalarContainer<double> **min_active_edge_dist_;
    // liquid content
    class ScalarContainer<double> **liquid_content_;

    // general implementation
    class ScalarContainer<double> ***scalar_containers_;
    char **scalar_container_names_;
    int n_scalar_containers_;
    class VectorContainer<double,3> ***vector_containers_;
    char **vector_container_names_;
    int n_vector_containers_;

    char **container_args_;
    int n_container_bases_;

    int modify_param(int, char **);
    void getGeneralRefs();

    void write_header(bigint ndump){}
    void getRefs();
    void pack(int *);
    void write_data(int, double *);
};

}

#endif
#endif // LAMMPS_VTK
