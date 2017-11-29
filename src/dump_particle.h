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

    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) 

#ifndef LMP_DUMP_PARTICLE_H
#define LMP_DUMP_PARTICLE_H

#include "pointers.h"
#include "sort_buffer.h"
#include <map>
#include <set>
#include <string>
#include <list>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyVertex.h>
#include <vtkMultiBlockDataSet.h>

class vtkAbstractArray;
class vtkRectilinearGrid;
class vtkUnstructuredGrid;

namespace LAMMPS_NS {

/**
 * @brief DumpParticle class
 *        write atom data to vtk files.
 *
 * Similar to the DumpCustom class but uses the vtk library to write data to vtk simple
 * legacy or xml format depending on the filename extension specified. (Since this
 * conflicts with the way binary output is specified, dump_modify allows to set the
 * binary flag for this dump command explicitly).
 * In contrast to DumpCustom class the attributes to be packed are stored in a std::map
 * to avoid duplicate entries and enforce correct ordering of vector components (except
 * for computes and fixes - these have to be given in the right order in the input script).
 * (Note: std::map elements are sorted by their keys.)
 * This dump command does not support compressed files, buffering or custom format strings,
 * multiproc is only supported by the xml formats, multifile option has to be used.
 */
class DumpParticle : public Pointers {
  public:
    DumpParticle(class LAMMPS *, int, int, int, int, int, int);
    virtual ~DumpParticle();
    int parse_parameters(int narg, const char *const *const arg, const bool pp_keyword_optional = false, std::list<std::string> keyword_list = std::list<std::string>());
    void prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool usePolyData = false);
    bigint memory_usage();
    virtual int modify_param(int, char **);
    virtual void init_style();
    int count();

  protected:
    int nevery;               // dump frequency for output
    int nclusterprocs;        // number of procs that write to one file
    int multiproc;            // number of procs writing files
    int filewriter;           // 1 if this proc writes a file, else 0
    int fileproc;             // ID of proc in my cluster who writes to file
    int iregion;              // -1 if no region, else which region
    char *idregion;           // region ID
    int igroup;               // group id
    int groupbit;             // group mask
    int nthresh;              // # of defined thresholds
    int *thresh_array;        // array to threshold on for each nthresh
    int *thresh_op;           // threshold operation for each nthresh
    double *thresh_value;     // threshold value for each nthresh

    int nchoose;              // # of selected atoms
    int maxlocal;             // size of atom selection and variable arrays
    int *choose;              // local indices of selected atoms
    double *dchoose;          // value for each atom to threshhold against
    int *clist;               // compressed list of indices of selected atoms

    int nfield;               // # of keywords listed by user
    int size_one;             // number of doubles used per particle

    std::map<int, int> field2index; // which compute,fix,variable calcs this field
    std::map<int, int> argindex;    // index into compute,fix scalar_atom,vector_atom
                                    // 0 for scalar_atom, 1-N for vector_atom values

    int ncompute;             // # of Compute objects used by dump
    char **id_compute;        // their IDs
    class Compute **compute;  // list of ptrs to the Compute objects

    int nfix;                 // # of Fix objects used by dump
    char **id_fix;            // their IDs
    class Fix **fix;          // list of ptrs to the Fix objects

    int nvariable;            // # of Variables used by dump
    char **id_variable;       // their names
    int *variable;            // list of indices for the Variables
    double **vbuf;            // local storage for variable evaluation

    int ntypes;               // # of atom types
    char **typenames;         // array of element names for each type

    int maxbuf;               // max size of buffer
    double *buf;              // array buffer

    // private methods

    void pack(int *);
    virtual void write_data(int, double *, vtkSmartPointer<vtkMultiBlockDataSet>, bool usePolyData);

    void identify_vectors();
    void identify_tensor();
    int add_compute(char *);
    int add_fix(char *);
    int add_variable(char *);

    void prepare_domain_data(vtkRectilinearGrid *);
    void prepare_domain_data_triclinic(vtkUnstructuredGrid *);
    void write_domain_vtk();
    void write_domain_vtk_triclinic();
    void write_domain_vtr();
    void write_domain_vtu_triclinic();

    typedef void (DumpParticle::*FnPtrPack)(int);
    std::map<int, FnPtrPack> pack_choice;  // ptrs to pack functions
    std::map<int, int> vtype;              // data type (INT, DOUBLE,...) for each type of entry for each atommyarrays
    std::map<int, std::string> name;       // attribute labels (e.g. "x", "v", ...) for each type of entry for each atommyarrays
    std::set<int> vector_set;              // set of vector attributes for each type of entry for each atommyarrays
    int current_pack_choice_key;

    // vtk data containers
    vtkSmartPointer<vtkPoints> points;                            // list of points, one point for each point cell
    vtkSmartPointer<vtkCellArray> pointsCells;                    // list of point cells
    std::map<int, vtkSmartPointer<vtkAbstractArray> > myarrays;   // list of a list of arrays that is presents data for each atom (x, v,...)
                                                                  // is then added to the point cells upon writing

    int n_calls_;
    double (*boxcorners)[3]; // corners of triclinic domain box

    bool convex_hull_detected;
    int convex_hull_max_n_tri;
    bool tensor_detected;

    double boxxlo,boxxhi;      // local copies of domain values
    double boxylo,boxyhi;      // lo/hi are bounding box for triclinic
    double boxzlo,boxzhi;

    SortBuffer *sortBuffer;

    void setFileCurrent();
    void buf2arrays(int, double *); // transfer data from buf array to vtk arrays
    void reset_vtk_data_containers();

    // customize by adding a method prototype
    void pack_compute(int);
    void pack_fix(int);
    void pack_variable(int);

    void pack_id(int);
    void pack_molecule(int);
    void pack_type(int);
    void pack_mass(int);

    void pack_x(int);
    void pack_y(int);
    void pack_z(int);
    void pack_points_convexhull(int);
    void pack_xs(int);
    void pack_ys(int);
    void pack_zs(int);
    void pack_xs_triclinic(int);
    void pack_ys_triclinic(int);
    void pack_zs_triclinic(int);
    void pack_xu(int);
    void pack_yu(int);
    void pack_zu(int);
    void pack_xu_triclinic(int);
    void pack_yu_triclinic(int);
    void pack_zu_triclinic(int);
    void pack_xsu(int);
    void pack_ysu(int);
    void pack_zsu(int);
    void pack_xsu_triclinic(int);
    void pack_ysu_triclinic(int);
    void pack_zsu_triclinic(int);
    void pack_ix(int);
    void pack_iy(int);
    void pack_iz(int);

    void pack_vx(int);
    void pack_vy(int);
    void pack_vz(int);
    void pack_fx(int);
    void pack_fy(int);
    void pack_fz(int);
    void pack_q(int);
    void pack_density(int); 
    void pack_p(int);       
    void pack_rho(int);     
    void pack_mux(int);
    void pack_muy(int);
    void pack_muz(int);
    void pack_mu(int);
    void pack_radius(int);
    void pack_diameter(int);

    void pack_omegax(int);
    void pack_omegay(int);
    void pack_omegaz(int);
    void pack_angmomx(int);
    void pack_angmomy(int);
    void pack_angmomz(int);
    void pack_tqx(int);
    void pack_tqy(int);
    void pack_tqz(int);
    void pack_spin(int);
    void pack_eradius(int);
    void pack_ervel(int);
    void pack_erforce(int);
    void pack_shapex(int); 
    void pack_shapey(int);
    void pack_shapez(int);
    void pack_blockiness1(int);
    void pack_blockiness2(int);
    void pack_quat1(int);
    void pack_quat2(int);
    void pack_quat3(int);
    void pack_quat4(int);
    void pack_tensor(int n);

    union ubuf {
        double   d;
        int64_t  i;
        ubuf(double arg) : d(arg) {}
        ubuf(int64_t arg) : i(arg) {}
        ubuf(int arg) : i(arg) {}
    };
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

E: Invalid attribute in dump custom command

Self-explantory.

E: Could not find dump custom compute ID

The compute ID needed by dump custom to compute a per-atom quantity
does not exist.

E: Could not find dump custom fix ID

Self-explanatory.

E: Dump custom and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump custom
needs them.

E: Could not find dump custom variable name

Self-explanatory.

E: Region ID for dump custom does not exist

Self-explanatory.

E: Threshhold for an atom property that isn't allocated

A dump threshhold has been requested on a quantity that is
not defined by the atom style used in this simulation.

E: Dumping an atom property that isn't allocated

The chosen atom style does not define the per-atom quantity being
dumped.

E: Dumping an atom quantity that isn't allocated

Only per-atom quantities that are defined for the atom style being
used are allowed.

E: Dump custom compute does not compute per-atom info

Self-explanatory.

E: Dump custom compute does not calculate per-atom vector

Self-explanatory.

E: Dump custom compute does not calculate per-atom array

Self-explanatory.

E: Dump custom compute vector is accessed out-of-range

Self-explanatory.

E: Dump custom fix does not compute per-atom info

Self-explanatory.

E: Dump custom fix does not compute per-atom vector

Self-explanatory.

E: Dump custom fix does not compute per-atom array

Self-explanatory.

E: Dump custom fix vector is accessed out-of-range

Self-explanatory.

E: Dump custom variable is not atom-style variable

Only atom-style variables generate per-atom quantities, needed for
dump output.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump_modify region ID does not exist

Self-explanatory.

E: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

E: Invalid attribute in dump modify command

Self-explantory.

E: Could not find dump modify compute ID

Self-explanatory.

E: Dump modify compute ID does not compute per-atom info

Self-explanatory.

E: Dump modify compute ID does not compute per-atom vector

Self-explanatory.

E: Dump modify compute ID does not compute per-atom array

Self-explanatory.

E: Dump modify compute ID vector is not large enough

Self-explanatory.

E: Could not find dump modify fix ID

Self-explanatory.

E: Dump modify fix ID does not compute per-atom info

Self-explanatory.

E: Dump modify fix ID does not compute per-atom vector

Self-explanatory.

E: Dump modify fix ID does not compute per-atom array

Self-explanatory.

E: Dump modify fix ID vector is not large enough

Self-explanatory.

E: Could not find dump modify variable name

Self-explanatory.

E: Dump modify variable is not atom-style variable

Self-explanatory.

E: Invalid dump_modify threshhold operator

Operator keyword used for threshold specification in not recognized.

*/
