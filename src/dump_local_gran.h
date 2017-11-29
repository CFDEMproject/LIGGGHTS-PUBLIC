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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) 

#ifndef LMP_DUMP_LOCAL_GRAN_H
#define LMP_DUMP_LOCAL_GRAN_H

#include "pointers.h"
#include "sort_buffer.h"
#include <map>
#include <set>
#include <string>
#include <list>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkMultiBlockDataSet.h>

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;

namespace LAMMPS_NS {

/**
 * @brief DumpLocalGran class
 *              write gran bond data to vtk arrays.
 *
 * Similar to the DumpParticle class but uses the vtk library to copy data to vtk data types
 * In contrast to the DumpParticle class the attributes to be packed are stored in a std::map
 * to avoid duplicate entries and enforce correct ordering of vector components (except
 * for computes and fixes - these have to be given in the right order in the input script).
 * (Note: std::map elements are sorted by their keys.)
 */
class DumpLocalGran : public Pointers {
  public:
    DumpLocalGran(LAMMPS *lmp, int _igroup, int _nclusterprocs, int _multiproc, int _nevery, int _filewriter, int _fileproc);
    ~DumpLocalGran();

    void prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool use_poly_data = false);
    int modify_param(int, char **);
    int parse_parameters(const int narg, const char *const *const arg, bool optional_keyword = false, std::list<std::string> keyword_list = std::list<std::string>());
    int count();
    virtual void init_style();
    bigint memory_usage();

  protected:

    int nevery;             // dump frequency for output
    int nclusterprocs;      // number of procs that write to one file
    int multiproc;          // number of procs writing files
    int filewriter;         // 1 if this proc writes a file, else 0
    int fileproc;           // ID of proc in my cluster who writes to file
    int iregion;            // -1 if no region, else which region
    char *idregion;         // region ID

    int igroup;             // group id
    int groupbit;           // group mask
    int nchoose;            // # of selected local datums

    SortBuffer *sortBuffer;  // to allow sorting of data

    int maxbuf;             // max size of buffer
    double *buf;            // main buffer
    int size_one;           // number of doubles used per particle

    class ComputePairGranLocal *cpgl_;

    // private methods

    void pack(int *);

    int parse_fields(int, char **);
    int add_compute(char *);
    int add_fix(char *);
    int add_variable(char *);

    void write_data(int n, double *mybuf, vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool use_poly_data);

    void define_properties();
    typedef void (DumpLocalGran::*FnPtrPack)(int);
    
    std::map<int, FnPtrPack> pack_choice;   // ptrs to pack functions
    std::map<int, int> vtype;               // data type for each attribute
    std::map<int, std::string> name;        // label for each attribute
    std::set<int> vector_set;               // set of vector attributes; defines which are vectors

    // vtk data containers
    vtkSmartPointer<vtkPoints> points;      // list of points, 2 points for each line cell
    vtkSmartPointer<vtkCellArray> lineCells;// list of line cells
    std::map<int, vtkSmartPointer<vtkAbstractArray> > myarrays; // list of a list of arrays that is presents data for each atom (x, v,...)
                                                                // is then added to the point cells upon writing
    int n_calls_;

    void buf2arrays(int, double *); // transfer data from buf array to vtk arrays
    void reset_vtk_data_containers();

    // customize by adding a method prototype
    void pack_x1(int);
    void pack_x2(int);
    void pack_v1(int);
    void pack_v2(int);
    void pack_id1(int);
    void pack_id2(int);
    void pack_f(int);
    void pack_fn(int);
    void pack_ft(int);
    void pack_torque(int);
    void pack_torquen(int);
    void pack_torquet(int);
    void pack_area(int);
    void pack_delta(int);
    void pack_heat(int);
    void pack_contact_point(int);
    void pack_ms_id1(int);
    void pack_ms_id2(int);
};

}

#endif
#endif
