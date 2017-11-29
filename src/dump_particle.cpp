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

#ifdef LAMMPS_VTK

#include <cmath>
#include "math_extra_liggghts.h"
#include <stdlib.h>
#include <string.h>
#include "dump_particle.h"

#ifdef CONVEX_ACTIVE_FLAG
#include "atom_vec_convexhull.h"
#endif

#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include <vector>
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkInformation.h>

// For compatibility with new VTK generic data arrays (VTK >= 7.0)
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace LAMMPS_NS;

// customize by
// * adding an enum constant (add vector components in consecutive order)
// * adding a pack_*(int) function for the value
// * adjusting parse_fields function to add the pack_* function to pack_choice
//   (in case of vectors, adjust identify_vectors as well)
// * adjusting thresh part in modify_param and count functions

enum{X,Y,Z, // required for vtk, must come first
     POINTS_CONVEXHULL, 
     ID,MOL,TYPE,ELEMENT,MASS,
     XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q, MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,SPIN,ERADIUS,ERVEL,ERFORCE,
     DENSITY, RHO, P, 
     VARIABLE,COMPUTE,FIX,
     SHAPEX, SHAPEY, SHAPEZ,
     QUAT1, QUAT2, QUAT3, QUAT4,
     EXTRA1, EXTRA2, TENSOR,
     BLOCKINESS1, BLOCKINESS2, 
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING,TENSOR_DOUBLE};      // same as in DumpCFG
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

/* ---------------------------------------------------------------------- */

DumpParticle::DumpParticle(LAMMPS *lmp, int _igroup, int _nclusterprocs, int _multiproc, int _nevery, int _filewriter, int _fileproc) :
    Pointers(lmp),
    nevery(_nevery),
    nclusterprocs(_nclusterprocs),
    multiproc(_multiproc),
    filewriter(_filewriter),
    fileproc(_fileproc),
    iregion(-1),
    idregion(NULL),
    igroup(_igroup),
    groupbit(group->bitmask[igroup]),
    nthresh(0),
    thresh_array(NULL),
    thresh_op(NULL),
    thresh_value(NULL),
    nchoose(0),
    maxlocal(0),
    choose(NULL),
    dchoose(NULL),
    clist(NULL),
    nfield(0),
    size_one(0),
    ncompute(0),
    id_compute(NULL),
    compute(NULL),
    nfix(0),
    id_fix(NULL),
    fix(NULL),
    nvariable(0),
    id_variable(NULL),
    variable(NULL),
    vbuf(NULL),
    ntypes(0),
    typenames(NULL),
    maxbuf(0),
    buf(NULL),
    current_pack_choice_key(0),
    n_calls_(0),
    convex_hull_detected(false),
    convex_hull_max_n_tri(0),
    tensor_detected(false),
    boxxlo(0.0),
    boxxhi(0.0),
    boxylo(0.0),
    boxyhi(0.0),
    boxzlo(0.0),
    boxzhi(0.0),
    sortBuffer(NULL)
{
#ifdef CONVEX_ACTIVE_FLAG
    if(dynamic_cast<AtomVecConvexHull*>(atom->avec))
    {
        convex_hull_detected = true;
        convex_hull_max_n_tri = dynamic_cast<AtomVecConvexHull*>(atom->avec)->get_ntri_max();
    }
#endif

    pack_choice.clear();
    vtype.clear();
    name.clear();
    myarrays.clear();
}

int DumpParticle::parse_parameters(const int narg, const char *const *const arg, const bool pp_keyword_optional, std::list<std::string> keyword_list)
{
    pack_choice[X] = &DumpParticle::pack_x;
    vtype[X] = DOUBLE;
    name[X] = "x";
    pack_choice[Y] = &DumpParticle::pack_y;
    vtype[Y] = DOUBLE;
    name[Y] = "y";
    pack_choice[Z] = &DumpParticle::pack_z;
    vtype[Z] = DOUBLE;
    name[Z] = "z";

    if(convex_hull_detected)
    {
        pack_choice[POINTS_CONVEXHULL] = &DumpParticle::pack_points_convexhull;
        vtype[POINTS_CONVEXHULL] = DOUBLE;
        name[POINTS_CONVEXHULL] = "points_convexhull";
    }
    int iarg = 0;

    // check if we have arguments at all, if not then just use positions
    if (iarg < narg)
    {
        if (pp_keyword_optional)
        {
            if (strcmp(arg[iarg], "particle_properties") == 0)
                iarg++;
        }
        else
        {
            if (strcmp(arg[iarg++], "particle_properties") != 0)
                error->all(FLERR, "Expected keyword 'particle_properties' in dump particle");
        }

        // customize by adding to if statement
        int i;
        for (; iarg < narg; iarg++) {
            i = iarg;

            if (strcmp(arg[iarg],"id") == 0) {
                pack_choice[ID] = &DumpParticle::pack_id;
                vtype[ID] = INT;
                name[ID] = arg[iarg];
            } else if (strcmp(arg[iarg],"mol") == 0 || strcmp(arg[iarg],"id_multisphere") == 0) {
                if (!atom->molecule_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[MOL] = &DumpParticle::pack_molecule;
                vtype[MOL] = INT;
                name[MOL] = arg[iarg];
            } else if (strcmp(arg[iarg],"type") == 0) {
                pack_choice[TYPE] = &DumpParticle::pack_type;
                vtype[TYPE] = INT;
                name[TYPE] =arg[iarg];
            } else if (strcmp(arg[iarg],"element") == 0) {
                pack_choice[ELEMENT] = &DumpParticle::pack_type;
                vtype[ELEMENT] = STRING;
                name[ELEMENT] = arg[iarg];
            } else if (strcmp(arg[iarg],"mass") == 0) {
                pack_choice[MASS] = &DumpParticle::pack_mass;
                vtype[MASS] = DOUBLE;
                name[MASS] = arg[iarg];

            } else if (strcmp(arg[iarg],"x") == 0) {
                // required property
            } else if (strcmp(arg[iarg],"y") == 0) {
                // required property
            } else if (strcmp(arg[iarg],"z") == 0) {
                // required property
            } else if (strcmp(arg[iarg],"xs") == 0) {
                if (domain->triclinic) pack_choice[XS] = &DumpParticle::pack_xs_triclinic;
                else pack_choice[XS] = &DumpParticle::pack_xs;
                vtype[XS] = DOUBLE;
                name[XS] = arg[iarg];
            } else if (strcmp(arg[iarg],"ys") == 0) {
                if (domain->triclinic) pack_choice[YS] = &DumpParticle::pack_ys_triclinic;
                else pack_choice[YS] = &DumpParticle::pack_ys;
                vtype[YS] = DOUBLE;
                name[YS] = arg[iarg];
            } else if (strcmp(arg[iarg],"zs") == 0) {
                if (domain->triclinic) pack_choice[ZS] = &DumpParticle::pack_zs_triclinic;
                else pack_choice[ZS] = &DumpParticle::pack_zs;
                vtype[ZS] = DOUBLE;
                name[ZS] = arg[iarg];
            } else if (strcmp(arg[iarg],"xu") == 0) {
                if (domain->triclinic) pack_choice[XU] = &DumpParticle::pack_xu_triclinic;
                else pack_choice[XU] = &DumpParticle::pack_xu;
                vtype[XU] = DOUBLE;
                name[XU] = arg[iarg];
            } else if (strcmp(arg[iarg],"yu") == 0) {
                if (domain->triclinic) pack_choice[YU] = &DumpParticle::pack_yu_triclinic;
                else pack_choice[YU] = &DumpParticle::pack_yu;
                vtype[YU] = DOUBLE;
                name[YU] = arg[iarg];
            } else if (strcmp(arg[iarg],"zu") == 0) {
                if (domain->triclinic) pack_choice[ZU] = &DumpParticle::pack_zu_triclinic;
                else pack_choice[ZU] = &DumpParticle::pack_zu;
                vtype[ZU] = DOUBLE;
                name[ZU] = arg[iarg];
            } else if (strcmp(arg[iarg],"xsu") == 0) {
                if (domain->triclinic) pack_choice[XSU] = &DumpParticle::pack_xsu_triclinic;
                else pack_choice[XSU] = &DumpParticle::pack_xsu;
                vtype[XSU] = DOUBLE;
                name[XSU] = arg[iarg];
            } else if (strcmp(arg[iarg],"ysu") == 0) {
                if (domain->triclinic) pack_choice[YSU] = &DumpParticle::pack_ysu_triclinic;
                else pack_choice[YSU] = &DumpParticle::pack_ysu;
                vtype[YSU] = DOUBLE;
                name[YSU] = arg[iarg];
            } else if (strcmp(arg[iarg],"zsu") == 0) {
                if (domain->triclinic) pack_choice[ZSU] = &DumpParticle::pack_zsu_triclinic;
                else pack_choice[ZSU] = &DumpParticle::pack_zsu;
                vtype[ZSU] = DOUBLE;
                name[ZSU] = arg[iarg];
            } else if (strcmp(arg[iarg],"ix") == 0) {
                pack_choice[IX] = &DumpParticle::pack_ix;
                vtype[IX] = INT;
                name[IX] = arg[iarg];
            } else if (strcmp(arg[iarg],"iy") == 0) {
                pack_choice[IY] = &DumpParticle::pack_iy;
                vtype[IY] = INT;
                name[IY] = arg[iarg];
            } else if (strcmp(arg[iarg],"iz") == 0) {
                pack_choice[IZ] = &DumpParticle::pack_iz;
                vtype[IZ] = INT;
                name[IZ] = arg[iarg];

            } else if (strcmp(arg[iarg],"vx") == 0) {
                pack_choice[VX] = &DumpParticle::pack_vx;
                vtype[VX] = DOUBLE;
                name[VX] = arg[iarg];
            } else if (strcmp(arg[iarg],"vy") == 0) {
                pack_choice[VY] = &DumpParticle::pack_vy;
                vtype[VY] = DOUBLE;
                name[VY] = arg[iarg];
            } else if (strcmp(arg[iarg],"vz") == 0) {
                pack_choice[VZ] = &DumpParticle::pack_vz;
                vtype[VZ] = DOUBLE;
                name[VZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"fx") == 0) {
                pack_choice[FX] = &DumpParticle::pack_fx;
                vtype[FX] = DOUBLE;
                name[FX] = arg[iarg];
            } else if (strcmp(arg[iarg],"fy") == 0) {
                pack_choice[FY] = &DumpParticle::pack_fy;
                vtype[FY] = DOUBLE;
                name[FY] = arg[iarg];
            } else if (strcmp(arg[iarg],"fz") == 0) {
                pack_choice[FZ] = &DumpParticle::pack_fz;
                vtype[FZ] = DOUBLE;
                name[FZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"q") == 0) {
                if (!atom->q_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[Q] = &DumpParticle::pack_q;
                vtype[Q] = DOUBLE;
                name[Q] = arg[iarg];
            } else if (strcmp(arg[iarg],"density") == 0) {
                if (!atom->density_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[DENSITY] = &DumpParticle::pack_density;
                vtype[DENSITY] = DOUBLE;
                name[DENSITY] = arg[iarg];
            } else if (strcmp(arg[iarg],"p") == 0) {
                if (!atom->p_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[P] = &DumpParticle::pack_p;
                vtype[P] = DOUBLE;
                name[P] = arg[iarg];
            } else if (strcmp(arg[iarg],"rho") == 0) {
                if (!atom->rho_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[RHO] = &DumpParticle::pack_rho;
                vtype[RHO] = DOUBLE;
                name[RHO] = arg[iarg];
            } else if (strcmp(arg[iarg],"mux") == 0) {
                if (!atom->mu_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[MUX] = &DumpParticle::pack_mux;
                vtype[MUX] = DOUBLE;
                name[MUX] = arg[iarg];
            } else if (strcmp(arg[iarg],"muy") == 0) {
                if (!atom->mu_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[MUY] = &DumpParticle::pack_muy;
                vtype[MUY] = DOUBLE;
                name[MUY] = arg[iarg];
            } else if (strcmp(arg[iarg],"muz") == 0) {
                if (!atom->mu_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[MUZ] = &DumpParticle::pack_muz;
                vtype[MUZ] = DOUBLE;
                name[MUZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"mu") == 0) {
                if (!atom->mu_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[MU] = &DumpParticle::pack_mu;
                vtype[MU] = DOUBLE;
                name[MU] = arg[iarg];
            } else if (strcmp(arg[iarg],"radius") == 0) {
                if (!atom->radius_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[RADIUS] = &DumpParticle::pack_radius;
                vtype[RADIUS] = DOUBLE;
                name[RADIUS] = arg[iarg];
            } else if (strcmp(arg[iarg],"diameter") == 0) {
                if (!atom->radius_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[DIAMETER] = &DumpParticle::pack_diameter;
                vtype[DIAMETER] = DOUBLE;
                name[DIAMETER] = arg[iarg];
            } else if (strcmp(arg[iarg],"omegax") == 0) {
                if (!atom->omega_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[OMEGAX] = &DumpParticle::pack_omegax;
                vtype[OMEGAX] = DOUBLE;
                name[OMEGAX] = arg[iarg];
            } else if (strcmp(arg[iarg],"omegay") == 0) {
                if (!atom->omega_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[OMEGAY] = &DumpParticle::pack_omegay;
                vtype[OMEGAY] = DOUBLE;
                name[OMEGAY] = arg[iarg];
            } else if (strcmp(arg[iarg],"omegaz") == 0) {
                if (!atom->omega_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[OMEGAZ] = &DumpParticle::pack_omegaz;
                vtype[OMEGAZ] = DOUBLE;
                name[OMEGAZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"angmomx") == 0) {
                if (!atom->angmom_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[ANGMOMX] = &DumpParticle::pack_angmomx;
                vtype[ANGMOMX] = DOUBLE;
                name[ANGMOMX] = arg[iarg];
            } else if (strcmp(arg[iarg],"angmomy") == 0) {
                if (!atom->angmom_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[ANGMOMY] = &DumpParticle::pack_angmomy;
                vtype[ANGMOMY] = DOUBLE;
                name[ANGMOMY] = arg[iarg];
            } else if (strcmp(arg[iarg],"angmomz") == 0) {
                if (!atom->angmom_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[ANGMOMZ] = &DumpParticle::pack_angmomz;
                vtype[ANGMOMZ] = DOUBLE;
                name[ANGMOMZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"tqx") == 0) {
                if (!atom->torque_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[TQX] = &DumpParticle::pack_tqx;
                vtype[TQX] = DOUBLE;
                name[TQX] = arg[iarg];
            } else if (strcmp(arg[iarg],"tqy") == 0) {
                if (!atom->torque_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[TQY] = &DumpParticle::pack_tqy;
                vtype[TQY] = DOUBLE;
                name[TQY] = arg[iarg];
            } else if (strcmp(arg[iarg],"tqz") == 0) {
                if (!atom->torque_flag)
                    error->all(FLERR,"Dumping an atom property that isn't allocated");
                pack_choice[TQZ] = &DumpParticle::pack_tqz;
                vtype[TQZ] = DOUBLE;
                name[TQZ] = arg[iarg];

            } else if (strcmp(arg[iarg],"spin") == 0) {
                if (!atom->spin_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[SPIN] = &DumpParticle::pack_spin;
                vtype[SPIN] = INT;
                name[SPIN] = arg[iarg];
            } else if (strcmp(arg[iarg],"eradius") == 0) {
                if (!atom->eradius_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[ERADIUS] = &DumpParticle::pack_eradius;
                vtype[ERADIUS] = DOUBLE;
                name[ERADIUS] = arg[iarg];
            } else if (strcmp(arg[iarg],"ervel") == 0) {
                if (!atom->ervel_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[ERVEL] = &DumpParticle::pack_ervel;
                vtype[ERVEL] = DOUBLE;
                name[ERVEL] = arg[iarg];
            } else if (strcmp(arg[iarg],"erforce") == 0) {
                if (!atom->erforce_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[ERFORCE] = &DumpParticle::pack_erforce;
                vtype[ERFORCE] = DOUBLE;
                name[ERFORCE] = arg[iarg];

            // compute value = c_ID
            // if no trailing [], then arg is set to 0, else arg is int between []

            } else if (strcmp(arg[iarg],"shapex") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[SHAPEX] = &DumpParticle::pack_shapex;
                vtype[SHAPEX] = DOUBLE;
                name[SHAPEX] = arg[iarg];
            } else if (strcmp(arg[iarg],"shapey") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[SHAPEY] = &DumpParticle::pack_shapey;
                vtype[SHAPEY] = DOUBLE;
                name[SHAPEY] = arg[iarg];
            } else if (strcmp(arg[iarg],"shapez") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[SHAPEZ] = &DumpParticle::pack_shapez;
                vtype[SHAPEZ] = DOUBLE;
                name[SHAPEZ] = arg[iarg];
            } else if (strcmp(arg[iarg],"blockiness1") == 0 or strcmp(arg[iarg],"roundness1") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                if(strcmp(arg[iarg],"roundness1") == 0)
                    error->warning(FLERR,"Keyword 'roundness1' will be deprecated in future, please use 'blockiness1' istead");
                pack_choice[BLOCKINESS1] = &DumpParticle::pack_blockiness1;
                vtype[BLOCKINESS1] = DOUBLE;
                name[BLOCKINESS1] = arg[iarg];
            } else if (strcmp(arg[iarg],"blockiness2") == 0 or strcmp(arg[iarg],"roundness2") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                if(strcmp(arg[iarg],"roundness2") == 0)
                    error->warning(FLERR,"Keyword 'roundness2' will be deprecated in future, please use 'blockiness2' istead");
                pack_choice[BLOCKINESS2] = &DumpParticle::pack_blockiness2;
                vtype[BLOCKINESS2] = DOUBLE;
                name[BLOCKINESS2] = arg[iarg];
            } else if (strcmp(arg[iarg],"quat1") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[QUAT1] = &DumpParticle::pack_quat1;
                vtype[QUAT1] = DOUBLE;
                name[QUAT1] = arg[iarg];
            } else if (strcmp(arg[iarg],"quat2") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[QUAT2] = &DumpParticle::pack_quat2;
                vtype[QUAT2] = DOUBLE;
                name[QUAT2] = arg[iarg];
            } else if (strcmp(arg[iarg],"quat3") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[QUAT3] = &DumpParticle::pack_quat3;
                vtype[QUAT3] = DOUBLE;
                name[QUAT3] = arg[iarg];
            } else if (strcmp(arg[iarg],"quat4") == 0) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,"Dumping an atom quantity that isn't allocated");
                pack_choice[QUAT4] = &DumpParticle::pack_quat4;
                vtype[QUAT4] = DOUBLE;
                name[QUAT4] = arg[iarg];
            } else if (strncmp(arg[iarg],"c_",2) == 0) {
                pack_choice[ATTRIBUTES+i] = &DumpParticle::pack_compute;
                vtype[ATTRIBUTES+i] = DOUBLE;

                int n = strlen(arg[iarg]);
                char *suffix = new char[n];
                strcpy(suffix,&arg[iarg][2]);

                char *ptr = strchr(suffix,'[');
                if (ptr)
                {
                    if (suffix[strlen(suffix)-1] != ']')
                        error->all(FLERR,"Invalid attribute in dump particle command");
                    argindex[ATTRIBUTES+i] = atoi(ptr+1);
                    *ptr = '\0';
                } else argindex[ATTRIBUTES+i] = 0;

                n = modify->find_compute(suffix);
                if (n < 0)
                    error->all(FLERR,"Could not find dump particle compute ID");
                if (modify->compute[n]->peratom_flag == 0)
                    error->all(FLERR,"Dump particle compute does not compute per-atom info");
                if (argindex[ATTRIBUTES+i] == 0 && modify->compute[n]->size_peratom_cols > 0)
                    error->all(FLERR, "Dump particle compute does not calculate per-atom vector");
                if (argindex[ATTRIBUTES+i] > 0 && modify->compute[n]->size_peratom_cols == 0)
                    error->all(FLERR, "Dump particle compute does not calculate per-atom array");
                if (argindex[ATTRIBUTES+i] > 0 &&
                        argindex[ATTRIBUTES+i] > modify->compute[n]->size_peratom_cols)
                    error->all(FLERR, "Dump particle compute vector is accessed out-of-range");

                field2index[ATTRIBUTES+i] = add_compute(suffix);
                name[ATTRIBUTES+i] = arg[iarg];
                delete [] suffix;

            // fix value = f_ID
            // if no trailing [], then arg is set to 0, else arg is between []

            } else if (strncmp(arg[iarg],"f_",2) == 0) {
                pack_choice[ATTRIBUTES+i] = &DumpParticle::pack_fix;
                vtype[ATTRIBUTES+i] = DOUBLE;

                int n = strlen(arg[iarg]);
                char *suffix = new char[n];
                strcpy(suffix,&arg[iarg][2]);

                char *ptr = strchr(suffix,'[');
                if (ptr) {
                    if (suffix[strlen(suffix)-1] != ']')
                        error->all(FLERR,"Invalid attribute in dump particle command");
                    argindex[ATTRIBUTES+i] = atoi(ptr+1);
                    *ptr = '\0';
                } else argindex[ATTRIBUTES+i] = 0;

                n = modify->find_fix(suffix);
                if (n < 0) error->all(FLERR,"Could not find dump particle fix ID");
                if (modify->fix[n]->peratom_flag == 0)
                    error->all(FLERR,"Dump particle fix does not compute per-atom info");
                if (argindex[ATTRIBUTES+i] == 0 && modify->fix[n]->size_peratom_cols > 0)
                    error->all(FLERR,"Dump particle fix does not compute per-atom vector");
                if (argindex[ATTRIBUTES+i] > 0 && modify->fix[n]->size_peratom_cols == 0)
                    error->all(FLERR,"Dump particle fix does not compute per-atom array");
                if (argindex[ATTRIBUTES+i] > 0 &&
                        argindex[ATTRIBUTES+i] > modify->fix[n]->size_peratom_cols)
                    error->all(FLERR,"Dump particle fix vector is accessed out-of-range");

                field2index[ATTRIBUTES+i] = add_fix(suffix);
                name[ATTRIBUTES+i] = arg[iarg];
                delete [] suffix;

            // variable value = v_name

            } else if (strncmp(arg[iarg],"v_",2) == 0) {
                pack_choice[ATTRIBUTES+i] = &DumpParticle::pack_variable;
                vtype[ATTRIBUTES+i] = DOUBLE;

                int n = strlen(arg[iarg]);
                char *suffix = new char[n];
                strcpy(suffix,&arg[iarg][2]);

                argindex[ATTRIBUTES+i] = 0;

                n = input->variable->find(suffix);
                if (n < 0) error->all(FLERR,"Could not find dump particle variable name");
                if (input->variable->atomstyle(n) == 0)
                    error->all(FLERR,"Dump particle variable is not atom-style variable");

                field2index[ATTRIBUTES+i] = add_variable(suffix);
                name[ATTRIBUTES+i] = suffix;
                delete [] suffix;

            }
            else
            {
                std::list<std::string>::iterator it;
                bool found_keyword = false;
                for (it = keyword_list.begin(); it != keyword_list.end(); it++)
                {
                    if (it->compare(arg[iarg]) == 0)
                    {
                        found_keyword = true;
                        break;
                    }
                }
                if (found_keyword)
                    break;
                else
                    error->all(FLERR, "Could not identify particle_property in dump/catalyst");
            }
        }
    }

    identify_vectors();
    identify_tensor();

    // TODO size_one for multiple tensors
    // TODO size_one changes for vectors?
    nfield = iarg - 1;
    size_one = pack_choice.size();
    if (tensor_detected)
        size_one += 8;

#if defined(NONSPHERICAL_ACTIVE_FLAG) && defined(CONVEX_ACTIVE_FLAG)
    if(convex_hull_detected)
    {
        int ntri_max = static_cast<AtomVecConvexHull*>(atom->avec)->get_ntri_max();
        size_one += (3*3*ntri_max);
    }
#endif

    if (filewriter)
        reset_vtk_data_containers();

    return iarg;
}

/* ---------------------------------------------------------------------- */

DumpParticle::~DumpParticle()
{
    if (idregion)
        delete [] idregion;
    memory->destroy(thresh_array);
    memory->destroy(thresh_op);
    memory->destroy(thresh_value);

    for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
    memory->sfree(id_compute);
    delete [] compute;

    for (int i = 0; i < nfix; i++) delete [] id_fix[i];
    memory->sfree(id_fix);
    delete [] fix;

    for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
    memory->sfree(id_variable);
    delete [] variable;
    for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
    delete [] vbuf;

    memory->destroy(choose);
    memory->destroy(dchoose);
    memory->destroy(clist);
    memory->destroy(buf);

    if (typenames) {
        for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
        delete [] typenames;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::init_style()
{
    // default for element names = C

    if (typenames == NULL) {
        typenames = new char*[ntypes+1];
        for (int itype = 1; itype <= ntypes; itype++) {
            typenames[itype] = new char[2];
            strcpy(typenames[itype],"C");
        }
    }

    // find current ptr for each compute,fix,variable
    // check that fix frequency is acceptable

    int icompute;
    for (int i = 0; i < ncompute; i++) {
        icompute = modify->find_compute(id_compute[i]);
        if (icompute < 0) error->all(FLERR,"Could not find dump particle compute ID");
        compute[i] = modify->compute[icompute];
    }

    int ifix;
    for (int i = 0; i < nfix; i++) {
        ifix = modify->find_fix(id_fix[i]);
        if (ifix < 0) error->all(FLERR,"Could not find dump particle fix ID");
        fix[i] = modify->fix[ifix];
        if (nevery % modify->fix[ifix]->peratom_freq)
            error->all(FLERR,"Dump particle and fix not computed at compatible times");
    }

    int ivariable;
    for (int i = 0; i < nvariable; i++) {
        ivariable = input->variable->find(id_variable[i]);
        if (ivariable < 0)
            error->all(FLERR,"Could not find dump particle variable name");
        variable[i] = ivariable;
    }

    // set index and check validity of region

    if (iregion >= 0) {
        iregion = domain->find_region(idregion);
        if (iregion == -1)
            error->all(FLERR,"Region ID for dump particle does not exist");
    }

#ifdef CONVEX_ACTIVE_FLAG
    if(dynamic_cast<AtomVecConvexHull*>(atom->avec))
    {
        convex_hull_detected = true;
        convex_hull_max_n_tri = dynamic_cast<AtomVecConvexHull*>(atom->avec)->get_ntri_max();
    }
    else
#endif
        convex_hull_detected = false;

    if (sortBuffer)
    {
        sortBuffer->init(igroup);

        if (sortBuffer->sort_set())
        {
            if (multiproc > 1)
                error->all(FLERR,
                           "Cannot dump sort when multiple procs write the dump file");
            if (sortBuffer->get_sortcol() == 0 && atom->tag_enable == 0)
                error->all(FLERR,"Cannot dump sort on atom IDs with no atom IDs defined");
            if (sortBuffer->get_sortcol() && sortBuffer->get_sortcol() > size_one)
                error->all(FLERR,"Dump sort column is invalid");
        }
    }
}

/* ---------------------------------------------------------------------- */

int DumpParticle::count()
{
    n_calls_ = 0;

    int i;

    // grow choose and variable vbuf arrays if needed

    int nlocal = atom->nlocal;
    if (nlocal > maxlocal) {
        maxlocal = atom->nmax;

        memory->destroy(choose);
        memory->destroy(dchoose);
        memory->destroy(clist);
        memory->create(choose,maxlocal,"dump:choose");
        memory->create(dchoose,maxlocal,"dump:dchoose");
        memory->create(clist,maxlocal,"dump:clist");

        for (i = 0; i < nvariable; i++) {
            memory->destroy(vbuf[i]);
            memory->create(vbuf[i],maxlocal,"dump:vbuf");
        }
    }

    // invoke Computes for per-atom quantities

    if (ncompute) {
        for (i = 0; i < ncompute; i++)
            if (!(compute[i]->invoked_flag & INVOKED_PERATOM)) {
                compute[i]->compute_peratom();
                compute[i]->invoked_flag |= INVOKED_PERATOM;
            }
    }

    // evaluate atom-style Variables for per-atom quantities

    if (nvariable)
        for (i = 0; i < nvariable; i++)
            input->variable->compute_atom(variable[i],igroup,vbuf[i],1,0);

    // choose all local atoms for output

    for (i = 0; i < nlocal; i++) choose[i] = 1;

    // un-choose if not in group

    if (igroup) {
        int *mask = atom->mask;
        for (i = 0; i < nlocal; i++)
            if (!(mask[i] & groupbit))
                choose[i] = 0;
    }

    // un-choose if not in region

    if (iregion >= 0) {
        Region *region = domain->regions[iregion];
        double **x = atom->x;
        for (i = 0; i < nlocal; i++)
            if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0)
                choose[i] = 0;
    }

    // un-choose if any threshold criterion isn't met

    if (nthresh) {
        double *ptr;
        double value;
        int nstride;
        int nlocal = atom->nlocal;

        for (int ithresh = 0; ithresh < nthresh; ithresh++) {

            // customize by adding to if statement

            if (thresh_array[ithresh] == ID) {
                int *tag = atom->tag;
                for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == MOL) {
                if (!atom->molecule_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                int *molecule = atom->molecule;
                for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == TYPE) {
                int *type = atom->type;
                for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ELEMENT) {
                int *type = atom->type;
                for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == MASS) {
                if (atom->rmass) {
                    ptr = atom->rmass;
                    nstride = 1;
                } else {
                    double *mass = atom->mass;
                    int *type = atom->type;
                    for (i = 0; i < nlocal; i++) dchoose[i] = mass[type[i]];
                    ptr = dchoose;
                    nstride = 1;
                }

            } else if (thresh_array[ithresh] == X) {
                ptr = &atom->x[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == Y) {
                ptr = &atom->x[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == Z) {
                ptr = &atom->x[0][2];
                nstride = 3;

            } else if (thresh_array[ithresh] == XS) {
                double **x = atom->x;
                double boxxlo = domain->boxlo[0];
                double invxprd = 1.0/domain->xprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (x[i][0] - boxxlo) * invxprd;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == YS) {
                double **x = atom->x;
                double boxylo = domain->boxlo[1];
                double invyprd = 1.0/domain->yprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (x[i][1] - boxylo) * invyprd;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ZS) {
                double **x = atom->x;
                double boxzlo = domain->boxlo[2];
                double invzprd = 1.0/domain->zprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (x[i][2] - boxzlo) * invzprd;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == XSTRI) {
                double **x = atom->x;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
                        h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == YSTRI) {
                double **x = atom->x;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
                        h_inv[3]*(x[i][2]-boxlo[2]);
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ZSTRI) {
                double **x = atom->x;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]);
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == XU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double xprd = domain->xprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = x[i][0] + ((image[i] & IMGMASK) - IMGMAX) * xprd;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == YU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double yprd = domain->yprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = x[i][1] +
                        ((image[i] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ZU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double zprd = domain->zprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = x[i][2] + ((image[i] >> IMG2BITS) - IMGMAX) * zprd;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == XUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *h = domain->h;
                int xbox,ybox,zbox;
                for (i = 0; i < nlocal; i++) {
                    xbox = (image[i] & IMGMASK) - IMGMAX;
                    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
                    zbox = (image[i] >> IMG2BITS) - IMGMAX;
                    dchoose[i] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
                }
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == YUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *h = domain->h;
                int ybox,zbox;
                for (i = 0; i < nlocal; i++) {
                    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
                    zbox = (image[i] >> IMG2BITS) - IMGMAX;
                    dchoose[i] = x[i][1] + h[1]*ybox + h[3]*zbox;
                }
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ZUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *h = domain->h;
                int zbox;
                for (i = 0; i < nlocal; i++) {
                    zbox = (image[i] >> IMG2BITS) - IMGMAX;
                    dchoose[i] = x[i][2] + h[2]*zbox;
                }
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == XSU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double boxxlo = domain->boxlo[0];
                double invxprd = 1.0/domain->xprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (x[i][0] - boxxlo) * invxprd +
                        (image[i] & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == YSU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double boxylo = domain->boxlo[1];
                double invyprd = 1.0/domain->yprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] =
                        (x[i][1] - boxylo) * invyprd +
                        (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == ZSU) {
                double **x = atom->x;
                tagint *image = atom->image;
                double boxzlo = domain->boxlo[2];
                double invzprd = 1.0/domain->zprd;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (x[i][2] - boxzlo) * invzprd +
                        (image[i] >> IMG2BITS) - IMGMAX;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == XSUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
                        h_inv[5]*(x[i][1]-boxlo[1]) +
                        h_inv[4]*(x[i][2]-boxlo[2]) +
                        (image[i] & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == YSUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
                        h_inv[3]*(x[i][2]-boxlo[2]) +
                        (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ZSUTRI) {
                double **x = atom->x;
                tagint *image = atom->image;
                double *boxlo = domain->boxlo;
                double *h_inv = domain->h_inv;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]) +
                        (image[i] >> IMG2BITS) - IMGMAX;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == IX) {
                tagint *image = atom->image;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (image[i] & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == IY) {
                tagint *image = atom->image;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == IZ) {
                tagint *image = atom->image;
                for (i = 0; i < nlocal; i++)
                    dchoose[i] = (image[i] >> IMG2BITS) - IMGMAX;
                ptr = dchoose;
                nstride = 1;

            } else if (thresh_array[ithresh] == VX) {
                ptr = &atom->v[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == VY) {
                ptr = &atom->v[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == VZ) {
                ptr = &atom->v[0][2];
                nstride = 3;
            } else if (thresh_array[ithresh] == FX) {
                ptr = &atom->f[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == FY) {
                ptr = &atom->f[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == FZ) {
                ptr = &atom->f[0][2];
                nstride = 3;

            } else if (thresh_array[ithresh] == Q) {
                if (!atom->q_flag)
                    error->all(FLERR,"Threshhold for an atom property that isn't allocated");
                ptr = atom->q;
                nstride = 1;
            } else if (thresh_array[ithresh] == P) { 
                if (!atom->p_flag)
                    error->all(FLERR,"Threshold for an atom property that isn't allocated");
                ptr = atom->p;
                nstride = 1;
            } else if (thresh_array[ithresh] == RHO) { 
                if (!atom->rho_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = atom->rho;
                nstride = 1;
            } else if (thresh_array[ithresh] == DENSITY) { 
                if (!atom->density_flag)
                    error->all(FLERR,"Threshold for an atom property that isn't allocated");
                ptr = atom->density;
                nstride = 1;
            } else if (thresh_array[ithresh] == MUX) {
                if (!atom->mu_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->mu[0][0];
                nstride = 4;
            } else if (thresh_array[ithresh] == MUY) {
                if (!atom->mu_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->mu[0][1];
                nstride = 4;
            } else if (thresh_array[ithresh] == MUZ) {
                if (!atom->mu_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->mu[0][2];
                nstride = 4;
            } else if (thresh_array[ithresh] == MU) {
                if (!atom->mu_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->mu[0][3];
                nstride = 4;

            } else if (thresh_array[ithresh] == RADIUS) {
                if (!atom->radius_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = atom->radius;
                nstride = 1;
            } else if (thresh_array[ithresh] == DIAMETER) {
                if (!atom->radius_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                double *radius = atom->radius;
                for (i = 0; i < nlocal; i++) dchoose[i] = 2.0*radius[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == OMEGAX) {
                if (!atom->omega_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->omega[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == OMEGAY) {
                if (!atom->omega_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->omega[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == OMEGAZ) {
                if (!atom->omega_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->omega[0][2];
                nstride = 3;
            } else if (thresh_array[ithresh] == ANGMOMX) {
                if (!atom->angmom_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->angmom[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == ANGMOMY) {
                if (!atom->angmom_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->angmom[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == ANGMOMZ) {
                if (!atom->angmom_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->angmom[0][2];
                nstride = 3;
            } else if (thresh_array[ithresh] == TQX) {
                if (!atom->torque_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->torque[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == TQY) {
                if (!atom->torque_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->torque[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == TQZ) {
                if (!atom->torque_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = &atom->torque[0][2];
                nstride = 3;

            } else if (thresh_array[ithresh] == SPIN) {
                if (!atom->spin_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                int *spin = atom->spin;
                for (i = 0; i < nlocal; i++) dchoose[i] = spin[i];
                ptr = dchoose;
                nstride = 1;
            } else if (thresh_array[ithresh] == ERADIUS) {
                if (!atom->eradius_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = atom->eradius;
                nstride = 1;
            } else if (thresh_array[ithresh] == ERVEL) {
                if (!atom->ervel_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = atom->ervel;
                nstride = 1;
            } else if (thresh_array[ithresh] == ERFORCE) {
                if (!atom->erforce_flag)
                    error->all(FLERR,
                                         "Threshold for an atom property that isn't allocated");
                ptr = atom->erforce;
                nstride = 1;

            } else if (thresh_array[ithresh] == COMPUTE) {
                i = ATTRIBUTES + nfield + ithresh;
                if (argindex[i] == 0) {
                    ptr = compute[field2index[i]]->vector_atom;
                    nstride = 1;
                } else {
                    ptr = &compute[field2index[i]]->array_atom[0][argindex[i]-1];
                    nstride = compute[field2index[i]]->size_peratom_cols;
                }

            } else if (thresh_array[ithresh] == FIX) {
                i = ATTRIBUTES + nfield + ithresh;
                if (argindex[i] == 0) {
                    ptr = fix[field2index[i]]->vector_atom;
                    nstride = 1;
                } else {
                    ptr = &fix[field2index[i]]->array_atom[0][argindex[i]-1];
                    nstride = fix[field2index[i]]->size_peratom_cols;
                }

            } else if (thresh_array[ithresh] == VARIABLE) {
                i = ATTRIBUTES + nfield + ithresh;
                ptr = vbuf[field2index[i]];
                nstride = 1;
            } else if (thresh_array[ithresh] == SHAPEX) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,
                        "Threshold for an atom property that isn't allocated");
                ptr = &atom->shape[0][0];
                nstride = 3;
            } else if (thresh_array[ithresh] == SHAPEY) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,
                        "Threshold for an atom property that isn't allocated");
                ptr = &atom->shape[0][1];
                nstride = 3;
            } else if (thresh_array[ithresh] == SHAPEZ) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,
                        "Threshold for an atom property that isn't allocated");
                ptr = &atom->shape[0][2];
                nstride = 3;
            } else if (thresh_array[ithresh] == BLOCKINESS1) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,
                        "Threshold for an atom property that isn't allocated");
                ptr = &atom->blockiness[0][0];
                nstride = 2;
            } else if (thresh_array[ithresh] == BLOCKINESS2) {
                if (!atom->superquadric_flag)
                    error->all(FLERR,
                        "Threshold for an atom property that isn't allocated");
                ptr = &atom->blockiness[0][1];
                nstride = 2;
            }

            // unselect atoms that don't meet threshold criterion

            value = thresh_value[ithresh];

            switch (thresh_op[ithresh]) {
            case LT:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr >= value) choose[i] = 0;
                break;
            case LE:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr > value) choose[i] = 0;
                break;
            case GT:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr <= value) choose[i] = 0;
                break;
            case GE:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr < value) choose[i] = 0;
                break;
            case EQ:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr != value) choose[i] = 0;
                break;
            case NEQ:
                for (i = 0; i < nlocal; i++, ptr += nstride)
                    if (choose[i] && *ptr == value) choose[i] = 0;
                break;
            }
        }
    }

    // compress choose flags into clist
    // nchoose = # of selected atoms
    // clist[i] = local index of each selected atom

    nchoose = 0;
    for (i = 0; i < nlocal; i++)
        if (choose[i]) clist[nchoose++] = i;

    return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpParticle::prepare_mbSet(vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool usePolyData)
{
    // simulation box bounds

    if (domain->triclinic == 0) {
        boxxlo = domain->boxlo[0];
        boxxhi = domain->boxhi[0];
        boxylo = domain->boxlo[1];
        boxyhi = domain->boxhi[1];
        boxzlo = domain->boxlo[2];
        boxzhi = domain->boxhi[2];
    } else {
        domain->box_corners();
        boxcorners = domain->corners;
    }

    // nme = # of dump lines this proc contributes to dump

    int nme = count();

    bigint ntotal = 0;
    bigint bnme = nme;
    MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

    int nmax;
    if (multiproc != comm->nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = nme;

    // ensure buf is sized for packing and communicating
    // use nmax to ensure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls

    if (nmax > maxbuf) {
        if ((bigint) nmax * size_one > MAXSMALLINT)
            error->all(FLERR,"Too much per-proc info for dump");
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*size_one,"dump:buf");
    }

    // ensure ids buffer is sized for sorting
    if (sortBuffer)
        sortBuffer->realloc_ids(nmax);

    // pack my data into buf
    // if sorting on IDs also request ID list from pack()
    // sort buf as needed

    if (sortBuffer)
        pack(sortBuffer->get_ids());
    else
        pack(NULL);
    if (sortBuffer)
        sortBuffer->sort(buf, nme, maxbuf, size_one, ntotal);

    // filewriter = 1 = this proc writes to file
    //   ping each proc in my cluster, receive its data, write data to file
    // else wait for ping from fileproc, send my data to fileproc

    int tmp,nlines;
    MPI_Status status;
    MPI_Request request;

    // comm and output buf of doubles

    if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
                MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,comm->me+iproc,0,world,&request);
                MPI_Send(&tmp,0,MPI_INT,comm->me+iproc,0,world);
                MPI_Wait(&request,&status);
                MPI_Get_count(&status,MPI_DOUBLE,&nlines);
                nlines /= size_one;
            } else nlines = nme;

            write_data(nlines, buf, mbSet, usePolyData);
        }
    } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
        MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack(int *ids)
{
    int n = 0;
    for (std::map<int,FnPtrPack>::iterator it=pack_choice.begin(); it!=pack_choice.end(); ++it, ++n) {
            current_pack_choice_key = it->first; // work-around for pack_compute, pack_fix, pack_variable
            (this->*(it->second))(n);
            if(current_pack_choice_key == TENSOR)
                n += 8;
            if(current_pack_choice_key == POINTS_CONVEXHULL)
                n += convex_hull_max_n_tri*3*3;
    }

    if (ids) {
        int *tag = atom->tag;
        for (int i = 0; i < nchoose; i++)
            ids[i] = tag[clist[i]];
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::buf2arrays(int n, double *mybuf)
{
    int ntri_max = 0;
#ifdef CONVEX_ACTIVE_FLAG
    if (convex_hull_detected)
        ntri_max = static_cast<AtomVecConvexHull*>(atom->avec)->get_ntri_max();
#endif
    // pid stores the ID(s) of the newly added point
    vtkIdType *pid;
    // pid size is the sum of:
    // 1.) maximum number of triangles of the convex hull * 3 (= for each vertex)
    // 2.) 1 (representing the position of the particle)
    const int pid_size = ntri_max*3 + 1;
    // safe upper bound to the number of points
    pid = new vtkIdType[pid_size];

    int ntri = 0;

    for (int iatom=0; iatom < n; ++iatom)
    {
        pid[0] = points->InsertNextPoint(mybuf[iatom*size_one],mybuf[iatom*size_one+1],mybuf[iatom*size_one+2]);
        int j=3; // 0,1,2 = x,y,z handled just above

        int npoints_extra = 0; 
        if(convex_hull_detected)
        {
            j += 3*3*ntri_max + 1;

            // each extra tri has 3 extra points
            ntri = mybuf[iatom*size_one+3];//(int) ubuf(mybuf[iatom*size_one+3]).i;
            npoints_extra = 3*mybuf[iatom*size_one+3];//(int) ubuf(mybuf[iatom*size_one+3]).i;

            for(int ipoint = 0; ipoint < npoints_extra; ipoint++)
            {
                if (1+ipoint >= pid_size)
                    error->all(FLERR,"Internal error: Overflow of pid array");

                pid[1+ipoint] = points->InsertNextPoint(mybuf[iatom*size_one+4+ipoint*3],
                                                        mybuf[iatom*size_one+5+ipoint*3],
                                                        mybuf[iatom*size_one+6+ipoint*3]);
            }
        }

        for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it)
        {
            vtkAbstractArray *paa = it->second;
            if (it->second->GetNumberOfComponents() == 3)
            {
                switch (vtype[it->first])
                {
                    case INT:
                    {
                        int iv3[3] = { static_cast<int>(mybuf[iatom*size_one+j  ]),
                                       static_cast<int>(mybuf[iatom*size_one+j+1]),
                                       static_cast<int>(mybuf[iatom*size_one+j+2]) };
                        vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
                        for(int ii = 0; ii < npoints_extra+1; ii++)
                            pia->InsertNextTupleValue(iv3);
                        break;
                    }
                    case DOUBLE:
                    {
                        vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
                        for(int ii = 0; ii < npoints_extra+1; ii++)
                            pda->InsertNextTupleValue(&mybuf[iatom*size_one+j]);
                        break;
                    }
                }
                j+=3;
            }
            else if (it->second->GetNumberOfComponents() == 9)
            {
                if(vtype[it->first] == TENSOR_DOUBLE)
                {
                    vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
                    for(int ii = 0; ii < npoints_extra+1; ii++)
                        pda->InsertNextTupleValue(&mybuf[iatom*size_one+j]);
                }
                else
                        error->all(FLERR,"Tensors of only double values are implemented!");
                j+=9;
            }
            else
            {
                switch (vtype[it->first])
                {
                    case INT:
                    {
                            vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
                            for(int ii = 0; ii < npoints_extra+1; ii++)
                                pia->InsertNextValue(mybuf[iatom*size_one+j]);
                            break;
                    }
                    case DOUBLE:
                    {
                            vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
                            for(int ii = 0; ii < npoints_extra+1; ii++)
                                pda->InsertNextValue(mybuf[iatom*size_one+j]);
                            break;
                    }
                    case STRING:
                    {
                            vtkStringArray *psa = static_cast<vtkStringArray*>(paa);
                            for(int ii = 0; ii < npoints_extra+1; ii++)
                                psa->InsertNextValue(typenames[static_cast<int>(mybuf[iatom*size_one+j])]);
                            break;
                    }
                }
                ++j;
            }
        }

        // 1 == npoints, not VTK_VERTEX (cell type), pid is the ID of the point added above
        if(!convex_hull_detected)
            pointsCells->InsertNextCell(1,pid);
        else
        {
            for(int itri = 0; itri < ntri; itri++)
            {
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

                triangle->GetPointIds()->SetNumberOfIds(3);
                for(int ipoint = 0; ipoint < 3; ipoint++)
                    triangle->GetPointIds()->SetId(ipoint,pid[1+itri*3+ipoint]);

                pointsCells->InsertNextCell(triangle);
            }

            /*
            vtkPolyVertex *polyvertex = vtkPolyVertex::New();

            polyvertex->GetPointIds()->SetNumberOfIds(npoints_extra+1);
            for(int ipoint = 0; ipoint < npoints_extra+1; ipoint++)
                    polyvertex->GetPointIds()->SetId(ipoint,pid[ipoint]);

            pointsCells->InsertNextCell(polyvertex);
            */
        }
    }

    delete[] pid;
}

/* ---------------------------------------------------------------------- */

void DumpParticle::prepare_domain_data(vtkRectilinearGrid *rgrid)
{
    vtkSmartPointer<vtkDoubleArray> xCoords =  vtkSmartPointer<vtkDoubleArray>::New();
    xCoords->InsertNextValue(boxxlo);
    xCoords->InsertNextValue(boxxhi);
    vtkSmartPointer<vtkDoubleArray> yCoords =  vtkSmartPointer<vtkDoubleArray>::New();
    yCoords->InsertNextValue(boxylo);
    yCoords->InsertNextValue(boxyhi);
    vtkSmartPointer<vtkDoubleArray> zCoords =  vtkSmartPointer<vtkDoubleArray>::New();
    zCoords->InsertNextValue(boxzlo);
    zCoords->InsertNextValue(boxzhi);

    rgrid->SetDimensions(2,2,2);
    rgrid->SetXCoordinates(xCoords);
    rgrid->SetYCoordinates(yCoords);
    rgrid->SetZCoordinates(zCoords);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::prepare_domain_data_triclinic(vtkUnstructuredGrid *hexahedronGrid)
{
    vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
    hexahedronPoints->SetNumberOfPoints(8);
    hexahedronPoints->InsertPoint(0, boxcorners[0][0], boxcorners[0][1], boxcorners[0][2]);
    hexahedronPoints->InsertPoint(1, boxcorners[1][0], boxcorners[1][1], boxcorners[1][2]);
    hexahedronPoints->InsertPoint(2, boxcorners[3][0], boxcorners[3][1], boxcorners[3][2]);
    hexahedronPoints->InsertPoint(3, boxcorners[2][0], boxcorners[2][1], boxcorners[2][2]);
    hexahedronPoints->InsertPoint(4, boxcorners[4][0], boxcorners[4][1], boxcorners[4][2]);
    hexahedronPoints->InsertPoint(5, boxcorners[5][0], boxcorners[5][1], boxcorners[5][2]);
    hexahedronPoints->InsertPoint(6, boxcorners[7][0], boxcorners[7][1], boxcorners[7][2]);
    hexahedronPoints->InsertPoint(7, boxcorners[6][0], boxcorners[6][1], boxcorners[6][2]);
    vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
    hexahedron->GetPointIds()->SetId(0, 0);
    hexahedron->GetPointIds()->SetId(1, 1);
    hexahedron->GetPointIds()->SetId(2, 2);
    hexahedron->GetPointIds()->SetId(3, 3);
    hexahedron->GetPointIds()->SetId(4, 4);
    hexahedron->GetPointIds()->SetId(5, 5);
    hexahedron->GetPointIds()->SetId(6, 6);
    hexahedron->GetPointIds()->SetId(7, 7);

    hexahedronGrid->Allocate(1, 1);
    hexahedronGrid->InsertNextCell(hexahedron->GetCellType(),
                                                                    hexahedron->GetPointIds());
    hexahedronGrid->SetPoints(hexahedronPoints);
}

/* ---------------------------------------------------------------------- */

void DumpParticle::write_data(int n, double *mybuf, vtkSmartPointer<vtkMultiBlockDataSet> mbSet, bool usePolyData)
{
    ++n_calls_;
    int cur_block = mbSet->GetNumberOfBlocks();

    buf2arrays(n, mybuf);

    if (n_calls_ < nclusterprocs)
        return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

    if (!usePolyData)
    {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

        unstructuredGrid->SetPoints(points);
        unstructuredGrid->SetCells(VTK_VERTEX, pointsCells);

        for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
            unstructuredGrid->GetPointData()->AddArray(it->second);
        }

        mbSet->SetBlock(cur_block++, unstructuredGrid);
    }
    else
    {
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

        polyData->SetPoints(points);
        if (!convex_hull_detected)
            polyData->SetVerts(pointsCells);
        else
            polyData->SetPolys(pointsCells);

        for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it)
            polyData->GetPointData()->AddArray(it->second);

        mbSet->SetBlock(cur_block++, polyData);
    }
    mbSet->GetMetaData(cur_block-1)->Set(mbSet->NAME(), "Particles");

    // DOMAIN DATA
    vtkSmartPointer<vtkDataSet> domainGrid;
    if (domain->triclinic)
    {
        domainGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        if (comm->me == 0)
            prepare_domain_data_triclinic(static_cast<vtkUnstructuredGrid*>(domainGrid.GetPointer()));
    }
    else
    {
        domainGrid = vtkSmartPointer<vtkRectilinearGrid>::New();
        if (comm->me == 0)
            prepare_domain_data(static_cast<vtkRectilinearGrid*>(domainGrid.GetPointer()));
    }
    mbSet->SetBlock(cur_block++, domainGrid);
    mbSet->GetMetaData(cur_block-1)->Set(mbSet->NAME(), "Domain");

    reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpParticle::reset_vtk_data_containers()
{
    points = vtkSmartPointer<vtkPoints>::New();
    pointsCells = vtkSmartPointer<vtkCellArray>::New();

    std::map<int,int>::iterator it=vtype.begin();
    ++it; ++it; ++it;

    if(convex_hull_detected)
        ++it;

    for (; it!=vtype.end(); ++it) {
        switch(vtype[it->first]) {
            case INT:
                myarrays[it->first] = vtkSmartPointer<vtkIntArray>::New();
                break;
            case DOUBLE:
                myarrays[it->first] = vtkSmartPointer<vtkDoubleArray>::New();
                break;
            case TENSOR_DOUBLE:
                myarrays[it->first] = vtkSmartPointer<vtkDoubleArray>::New();
                break;
            case STRING:
                myarrays[it->first] = vtkSmartPointer<vtkStringArray>::New();
                break;
        }

        if (vector_set.find(it->first) != vector_set.end()) {
            myarrays[it->first]->SetNumberOfComponents(3);
            myarrays[it->first]->SetName(name[it->first].c_str());
            ++it; ++it;
        } else if(vtype[it->first] == TENSOR_DOUBLE) {
            myarrays[it->first]->SetNumberOfComponents(9);
            myarrays[it->first]->SetName(name[it->first].c_str());
        } else {
            myarrays[it->first]->SetName(name[it->first].c_str());
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::identify_vectors()
{
    // detect vectors
    vector_set.insert(X); // required

    int vector3_starts[] = {XS, XU, XSU, IX, VX, FX, MUX, OMEGAX, ANGMOMX, TQX};
    int num_vector3_starts = sizeof(vector3_starts) / sizeof(int);

    for (int v3s = 0; v3s < num_vector3_starts; v3s++) {
        if(name.count(vector3_starts[v3s]  ) &&
             name.count(vector3_starts[v3s]+1) &&
             name.count(vector3_starts[v3s]+2) )
        {
            std::string vectorName = name[vector3_starts[v3s]];
            vectorName.erase(vectorName.find_first_of('x'));
            name[vector3_starts[v3s]] = vectorName;
            vector_set.insert(vector3_starts[v3s]);
        }
    }

    // compute and fix vectors
    for (std::map<int,std::string>::iterator it=name.begin(); it!=name.end(); ++it) {
        if (it->first < ATTRIBUTES) // neither fix nor compute
            continue;

        if(argindex[it->first] == 0) // single value
            continue;

        // assume components are grouped together and in correct order
        if(name.count(it->first + 1) && name.count(it->first + 2) ) { // more attributes?
            if(it->second.compare(0,it->second.length()-3,name[it->first + 1],0,it->second.length()-3) == 0  && // same attributes?
                 it->second.compare(0,it->second.length()-3,name[it->first + 2],0,it->second.length()-3) == 0 )
            {
                it->second.erase(it->second.length()-1);
                std::ostringstream oss;
                oss << "-" << argindex[it->first+2] << "]";
                it->second += oss.str();
                vector_set.insert(it->first);
                ++it; ++it;
            }
        }
    }
}

void DumpParticle::identify_tensor()
{
    if(pack_choice.count(QUAT1) > 0 &&
         pack_choice.count(QUAT2) > 0 &&
         pack_choice.count(QUAT3) > 0 &&
         pack_choice.count(QUAT4) > 0) {

             pack_choice.erase(QUAT1);
             pack_choice.erase(QUAT2);
             pack_choice.erase(QUAT3);
             pack_choice.erase(QUAT4);

             vtype.erase(QUAT1);
             vtype.erase(QUAT2);
             vtype.erase(QUAT3);
             vtype.erase(QUAT4);

             name.erase(QUAT1);
             name.erase(QUAT2);
             name.erase(QUAT3);
             name.erase(QUAT4);

             name[TENSOR] = "TENSOR";
             vtype[TENSOR] = TENSOR_DOUBLE;
             pack_choice[TENSOR] = &DumpParticle::pack_tensor;

             tensor_detected = true;
    }
}
/* ----------------------------------------------------------------------
     add Compute to list of Compute objects used by dump
     return index of where this Compute is in list
     if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParticle::add_compute(char *id)
{
    int icompute;
    for (icompute = 0; icompute < ncompute; icompute++)
        if (strcmp(id,id_compute[icompute]) == 0) break;
    if (icompute < ncompute) return icompute;

    id_compute = (char **)
        memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
    delete [] compute;
    compute = new Compute*[ncompute+1];

    int n = strlen(id) + 1;
    id_compute[ncompute] = new char[n];
    strcpy(id_compute[ncompute],id);
    ncompute++;
    return ncompute-1;
}

/* ----------------------------------------------------------------------
     add Fix to list of Fix objects used by dump
     return index of where this Fix is in list
     if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParticle::add_fix(char *id)
{
    int ifix;
    for (ifix = 0; ifix < nfix; ifix++)
        if (strcmp(id,id_fix[ifix]) == 0) break;
    if (ifix < nfix) return ifix;

    id_fix = (char **)
        memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
    delete [] fix;
    fix = new Fix*[nfix+1];

    int n = strlen(id) + 1;
    id_fix[nfix] = new char[n];
    strcpy(id_fix[nfix],id);
    nfix++;
    return nfix-1;
}

/* ----------------------------------------------------------------------
     add Variable to list of Variables used by dump
     return index of where this Variable is in list
     if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpParticle::add_variable(char *id)
{
    int ivariable;
    for (ivariable = 0; ivariable < nvariable; ivariable++)
        if (strcmp(id,id_variable[ivariable]) == 0) break;
    if (ivariable < nvariable) return ivariable;

    id_variable = (char **)
        memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                                         "dump:id_variable");
    delete [] variable;
    variable = new int[nvariable+1];
    delete [] vbuf;
    vbuf = new double*[nvariable+1];
    for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

    int n = strlen(id) + 1;
    id_variable[nvariable] = new char[n];
    strcpy(id_variable[nvariable],id);
    nvariable++;
    return nvariable-1;
}

/* ---------------------------------------------------------------------- */

int DumpParticle::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"region") == 0) {
        if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
        if (strcmp(arg[1],"none") == 0) iregion = -1;
        else {
            iregion = domain->find_region(arg[1]);
            if (iregion == -1)
                error->all(FLERR,"Dump_modify region ID does not exist");
            if (idregion)
                delete [] idregion;
            int n = strlen(arg[1]) + 1;
            idregion = new char[n];
            strcpy(idregion,arg[1]);
        }
        return 2;
    }

    if (strcmp(arg[0],"element") == 0) {
        if (narg < ntypes+1)
            error->all(FLERR,"Dump modify: number of element names do not match atom types");

        if (typenames) {
            for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
            delete [] typenames;
            typenames = NULL;
        }

        typenames = new char*[ntypes+1];
        for (int itype = 1; itype <= ntypes; itype++) {
            int n = strlen(arg[itype]) + 1;
            typenames[itype] = new char[n];
            strcpy(typenames[itype],arg[itype]);
        }
        return ntypes+1;
    }

    if (strcmp(arg[0],"thresh") == 0) {
        if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
        if (strcmp(arg[1],"none") == 0) {
            if (nthresh) {
                memory->destroy(thresh_array);
                memory->destroy(thresh_op);
                memory->destroy(thresh_value);
                thresh_array = NULL;
                thresh_op = NULL;
                thresh_value = NULL;
            }
            nthresh = 0;
            return 2;
        }

        if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

        // grow threshold arrays

        memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
        memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
        memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");

        // set attribute type of threshold
        // customize by adding to if statement

        if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
        else if (strcmp(arg[1],"mol") == 0 || strcmp(arg[1],"id_multisphere") == 0) thresh_array[nthresh] = MOL;
        else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
        else if (strcmp(arg[1],"mass") == 0) thresh_array[nthresh] = MASS;

        else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
        else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
        else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

        else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = XS;
        else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = XSTRI;
        else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = YS;
        else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = YSTRI;
        else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = ZS;
        else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = ZSTRI;

        else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = XU;
        else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = XUTRI;
        else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = YU;
        else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = YUTRI;
        else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = ZU;
        else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = ZUTRI;

        else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = XSU;
        else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = XSUTRI;
        else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = YSU;
        else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = YSUTRI;
        else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 0)
            thresh_array[nthresh] = ZSU;
        else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 1)
            thresh_array[nthresh] = ZSUTRI;

        else if (strcmp(arg[1],"ix") == 0) thresh_array[nthresh] = IX;
        else if (strcmp(arg[1],"iy") == 0) thresh_array[nthresh] = IY;
        else if (strcmp(arg[1],"iz") == 0) thresh_array[nthresh] = IZ;
        else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
        else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
        else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;
        else if (strcmp(arg[1],"fx") == 0) thresh_array[nthresh] = FX;
        else if (strcmp(arg[1],"fy") == 0) thresh_array[nthresh] = FY;
        else if (strcmp(arg[1],"fz") == 0) thresh_array[nthresh] = FZ;

        else if (strcmp(arg[1],"q") == 0) thresh_array[nthresh] = Q;
        else if (strcmp(arg[1],"density") == 0) thresh_array[nthresh] = DENSITY; 
        else if (strcmp(arg[1],"p") == 0) thresh_array[nthresh] = P; 
        else if (strcmp(arg[1],"rho") == 0) thresh_array[nthresh] = RHO; 
        else if (strcmp(arg[1],"mux") == 0) thresh_array[nthresh] = MUX;
        else if (strcmp(arg[1],"muy") == 0) thresh_array[nthresh] = MUY;
        else if (strcmp(arg[1],"muz") == 0) thresh_array[nthresh] = MUZ;
        else if (strcmp(arg[1],"mu") == 0) thresh_array[nthresh] = MU;

        else if (strcmp(arg[1],"radius") == 0) thresh_array[nthresh] = RADIUS;
        else if (strcmp(arg[1],"diameter") == 0) thresh_array[nthresh] = DIAMETER;
        else if (strcmp(arg[1],"omegax") == 0) thresh_array[nthresh] = OMEGAX;
        else if (strcmp(arg[1],"omegay") == 0) thresh_array[nthresh] = OMEGAY;
        else if (strcmp(arg[1],"omegaz") == 0) thresh_array[nthresh] = OMEGAZ;
        else if (strcmp(arg[1],"angmomx") == 0) thresh_array[nthresh] = ANGMOMX;
        else if (strcmp(arg[1],"angmomy") == 0) thresh_array[nthresh] = ANGMOMY;
        else if (strcmp(arg[1],"angmomz") == 0) thresh_array[nthresh] = ANGMOMZ;
        else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
        else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
        else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;

        else if (strcmp(arg[1],"spin") == 0) thresh_array[nthresh] = SPIN;
        else if (strcmp(arg[1],"eradius") == 0) thresh_array[nthresh] = ERADIUS;
        else if (strcmp(arg[1],"ervel") == 0) thresh_array[nthresh] = ERVEL;
        else if (strcmp(arg[1],"erforce") == 0) thresh_array[nthresh] = ERFORCE;
        else if (strcmp(arg[1],"shapex") == 0) thresh_array[nthresh] = SHAPEX;
        else if (strcmp(arg[1],"shapey") == 0) thresh_array[nthresh] = SHAPEY;
        else if (strcmp(arg[1],"shapez") == 0) thresh_array[nthresh] = SHAPEZ;
        else if (strcmp(arg[1],"blockiness1") == 0) thresh_array[nthresh] = BLOCKINESS1;
        else if (strcmp(arg[1],"blockiness2") == 0) thresh_array[nthresh] = BLOCKINESS2;
        else if (strcmp(arg[1],"quat1") == 0) thresh_array[nthresh] = QUAT1;
        else if (strcmp(arg[1],"quat2") == 0) thresh_array[nthresh] = QUAT2;
        else if (strcmp(arg[1],"quat3") == 0) thresh_array[nthresh] = QUAT3;
        else if (strcmp(arg[1],"quat4") == 0) thresh_array[nthresh] = QUAT4;

        // compute value = c_ID
        // if no trailing [], then arg is set to 0, else arg is between []

        else if (strncmp(arg[1],"c_",2) == 0) {
            thresh_array[nthresh] = COMPUTE;
            int n = strlen(arg[1]);
            char *suffix = new char[n];
            strcpy(suffix,&arg[1][2]);

            char *ptr = strchr(suffix,'[');
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']')
                    error->all(FLERR,"Invalid attribute in dump modify command");
                argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

            n = modify->find_compute(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump modify compute ID");

            if (modify->compute[n]->peratom_flag == 0)
                error->all(FLERR,
                                     "Dump modify compute ID does not compute per-atom info");
            if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
                    modify->compute[n]->size_peratom_cols > 0)
                error->all(FLERR,
                                     "Dump modify compute ID does not compute per-atom vector");
            if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
                    modify->compute[n]->size_peratom_cols == 0)
                error->all(FLERR,
                                     "Dump modify compute ID does not compute per-atom array");
            if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
                    argindex[ATTRIBUTES+nfield+nthresh] > modify->compute[n]->size_peratom_cols)
                error->all(FLERR,"Dump modify compute ID vector is not large enough");

            field2index[ATTRIBUTES+nfield+nthresh] = add_compute(suffix);
            delete [] suffix;

        // fix value = f_ID
        // if no trailing [], then arg is set to 0, else arg is between []

        } else if (strncmp(arg[1],"f_",2) == 0) {
            thresh_array[nthresh] = FIX;
            int n = strlen(arg[1]);
            char *suffix = new char[n];
            strcpy(suffix,&arg[1][2]);

            char *ptr = strchr(suffix,'[');
            if (ptr) {
                if (suffix[strlen(suffix)-1] != ']')
                    error->all(FLERR,"Invalid attribute in dump modify command");
                argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
                *ptr = '\0';
            } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

            n = modify->find_fix(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump modify fix ID");

            if (modify->fix[n]->peratom_flag == 0)
                error->all(FLERR,"Dump modify fix ID does not compute per-atom info");
            if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
                    modify->fix[n]->size_peratom_cols > 0)
                error->all(FLERR,"Dump modify fix ID does not compute per-atom vector");
            if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
                    modify->fix[n]->size_peratom_cols == 0)
                error->all(FLERR,"Dump modify fix ID does not compute per-atom array");
            if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
                    argindex[ATTRIBUTES+nfield+nthresh] > modify->fix[n]->size_peratom_cols)
                error->all(FLERR,"Dump modify fix ID vector is not large enough");

            field2index[ATTRIBUTES+nfield+nthresh] = add_fix(suffix);
            delete [] suffix;

        // variable value = v_ID

        } else if (strncmp(arg[1],"v_",2) == 0) {
            thresh_array[nthresh] = VARIABLE;
            int n = strlen(arg[1]);
            char *suffix = new char[n];
            strcpy(suffix,&arg[1][2]);

            argindex[ATTRIBUTES+nfield+nthresh] = 0;

            n = input->variable->find(suffix);
            if (n < 0) error->all(FLERR,"Could not find dump modify variable name");
            if (input->variable->atomstyle(n) == 0)
                error->all(FLERR,"Dump modify variable is not atom-style variable");

            field2index[ATTRIBUTES+nfield+nthresh] = add_variable(suffix);
            delete [] suffix;

        } else error->all(FLERR,"Invalid dump_modify threshold operator");

        // set operation type of threshold

        if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
        else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
        else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
        else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
        else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
        else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
        else error->all(FLERR,"Invalid dump_modify threshold operator");

        // set threshold value

        thresh_value[nthresh] = force->numeric(FLERR,arg[3]);

        nthresh++;
        return 4;
    }

    if (!sortBuffer)
        sortBuffer = new SortBuffer(lmp, false);

    return sortBuffer->modify_param(narg, arg);
}

/* ----------------------------------------------------------------------
     return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpParticle::memory_usage()
{
    int bytes = 0;
    bytes += memory->usage(choose,maxlocal);
    bytes += memory->usage(dchoose,maxlocal);
    bytes += memory->usage(clist,maxlocal);
    bytes += memory->usage(vbuf,nvariable,maxlocal);
    bytes += memory->usage(buf,maxbuf*size_one);
    if (sortBuffer)
        bytes += sortBuffer->memory_usage(size_one);
    return bytes;
}

/* ----------------------------------------------------------------------
     extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpParticle::pack_compute(int n)
{
    double *vector = compute[field2index[current_pack_choice_key]]->vector_atom;
    double **array = compute[field2index[current_pack_choice_key]]->array_atom;
    int index = argindex[current_pack_choice_key];

    if (index == 0) {
        for (int i = 0; i < nchoose; i++) {
            buf[n] = vector[clist[i]];
            n += size_one;
        }
    } else {
        index--;
        for (int i = 0; i < nchoose; i++) {
            buf[n] = array[clist[i]][index];
            n += size_one;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_fix(int n)
{
    double *vector = fix[field2index[current_pack_choice_key]]->vector_atom;
    double **array = fix[field2index[current_pack_choice_key]]->array_atom;
    int index = argindex[current_pack_choice_key];

    if (index == 0) {
        for (int i = 0; i < nchoose; i++) {
            buf[n] = vector[clist[i]];
            n += size_one;
        }
    } else {
        index--;
        for (int i = 0; i < nchoose; i++) {
            buf[n] = array[clist[i]][index];
            n += size_one;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_variable(int n)
{
    double *vector = vbuf[field2index[current_pack_choice_key]];

    for (int i = 0; i < nchoose; i++) {
        buf[n] = vector[clist[i]];
        n += size_one;
    }
}

/* ----------------------------------------------------------------------
     one method for every attribute dump particle can output
     the atom property is packed into buf starting at n with stride size_one
     customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpParticle::pack_id(int n)
{
    int *tag = atom->tag;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = tag[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_molecule(int n)
{
    int *molecule = atom->molecule;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = molecule[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_type(int n)
{
    int *type = atom->type;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = type[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_mass(int n)
{
    int *type = atom->type;
    double *mass = atom->mass;
    double *rmass = atom->rmass;

    if (rmass) {
        for (int i = 0; i < nchoose; i++) {
            buf[n] = rmass[clist[i]];
            n += size_one;
        }
    } else {
        for (int i = 0; i < nchoose; i++) {
            buf[n] = mass[type[clist[i]]];
            n += size_one;
        }
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_x(int n)
{
    double **x = atom->x;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = x[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_y(int n)
{
    double **x = atom->x;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = x[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_z(int n)
{
    double **x = atom->x;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = x[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_points_convexhull(int n)
{
#ifdef CONVEX_ACTIVE_FLAG
    int *shapetype = atom->shapetype;
    double **x = atom->x;
    double **quaternion = atom->quaternion;
    AtomVecConvexHull *avec = static_cast<AtomVecConvexHull*>(atom->avec);
    double point[3];

    for (int i = 0; i < nchoose; i++)
    {
        int shapetyp = shapetype[clist[i]];
        int ntri = avec->get_ntri(shapetyp);
        buf[n] = ntri; //ubuf(np).d;
        
        for(int itri = 0; itri < ntri; itri++)
        {
            for(int inode = 0; inode < 3; inode ++)
            {
                avec->get_tri_node(shapetyp,itri,inode,point);
                MathExtraLiggghts::vec_quat_rotate(point,quaternion[clist[i]]);
                vectorAdd3D(point,x[clist[i]],point);
                buf[n+itri*9+inode*3+1] = point[0];
                buf[n+itri*9+inode*3+2] = point[1];
                buf[n+itri*9+inode*3+3] = point[2];
            }
        }
        n += size_one;
    }
#endif
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xs(int n)
{
    double **x = atom->x;

    double boxxlo = domain->boxlo[0];
    double invxprd = 1.0/domain->xprd;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (x[clist[i]][0] - boxxlo) * invxprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ys(int n)
{
    double **x = atom->x;

    double boxylo = domain->boxlo[1];
    double invyprd = 1.0/domain->yprd;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (x[clist[i]][1] - boxylo) * invyprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zs(int n)
{
    double **x = atom->x;

    double boxzlo = domain->boxlo[2];
    double invzprd = 1.0/domain->zprd;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (x[clist[i]][2] - boxzlo) * invzprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xs_triclinic(int n)
{
    int j;
    double **x = atom->x;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
            h_inv[4]*(x[j][2]-boxlo[2]);
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ys_triclinic(int n)
{
    int j;
    double **x = atom->x;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]);
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zs_triclinic(int n)
{
    double **x = atom->x;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = h_inv[2]*(x[clist[i]][2]-boxlo[2]);
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double xprd = domain->xprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = x[j][0] + ((image[j] & IMGMASK) - IMGMAX) * xprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_yu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double yprd = domain->yprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = x[j][1] + ((image[j] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double zprd = domain->zprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = x[j][2] + ((image[j] >> IMG2BITS) - IMGMAX) * zprd;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *h = domain->h;
    int xbox,ybox,zbox;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        xbox = (image[j] & IMGMASK) - IMGMAX;
        ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[j] >> IMG2BITS) - IMGMAX;
        buf[n] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_yu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *h = domain->h;
    int ybox,zbox;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[j] >> IMG2BITS) - IMGMAX;
        buf[n] = x[j][1] + h[1]*ybox + h[3]*zbox;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *h = domain->h;
    int zbox;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        zbox = (image[j] >> IMG2BITS) - IMGMAX;
        buf[n] = x[j][2] + h[2]*zbox;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xsu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double boxxlo = domain->boxlo[0];
    double invxprd = 1.0/domain->xprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = (x[j][0] - boxxlo) * invxprd + (image[j] & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ysu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double boxylo = domain->boxlo[1];
    double invyprd = 1.0/domain->yprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = (x[j][1] - boxylo) * invyprd + (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zsu(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double boxzlo = domain->boxlo[2];
    double invzprd = 1.0/domain->zprd;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = (x[j][2] - boxzlo) * invzprd + (image[j] >> IMG2BITS) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_xsu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
            h_inv[4]*(x[j][2]-boxlo[2]) + (image[j] & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ysu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]) +
            (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_zsu_triclinic(int n)
{
    int j;
    double **x = atom->x;
    tagint *image = atom->image;

    double *boxlo = domain->boxlo;
    double *h_inv = domain->h_inv;

    for (int i = 0; i < nchoose; i++) {
        j = clist[i];
        buf[n] = h_inv[2]*(x[j][2]-boxlo[2]) + (image[j] >> IMG2BITS) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ix(int n)
{
    tagint *image = atom->image;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (image[clist[i]] & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_iy(int n)
{
    tagint *image = atom->image;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (image[clist[i]] >> IMGBITS & IMGMASK) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_iz(int n)
{
    tagint *image = atom->image;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = (image[clist[i]] >> IMG2BITS) - IMGMAX;
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_vx(int n)
{
    double **v = atom->v;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = v[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_vy(int n)
{
    double **v = atom->v;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = v[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_vz(int n)
{
    double **v = atom->v;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = v[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_fx(int n)
{
    double **f = atom->f;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = f[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_fy(int n)
{
    double **f = atom->f;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = f[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_fz(int n)
{
    double **f = atom->f;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = f[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_q(int n)
{
    double *q = atom->q;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = q[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_density(int n) 
{
    double *density = atom->density;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = density[clist[i]];
        n += size_one;
    }

}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_p(int n) 
{
    double *p = atom->p;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = p[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_rho(int n) 
{
    double *rho = atom->rho;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = rho[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_mux(int n)
{
    double **mu = atom->mu;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = mu[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_muy(int n)
{
    double **mu = atom->mu;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = mu[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_muz(int n)
{
    double **mu = atom->mu;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = mu[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_mu(int n)
{
    double **mu = atom->mu;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = mu[clist[i]][3];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_radius(int n)
{
    double *radius = atom->radius;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = radius[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_diameter(int n)
{
    double *radius = atom->radius;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = 2.0*radius[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_omegax(int n)
{
    double **omega = atom->omega;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = omega[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_omegay(int n)
{
    double **omega = atom->omega;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = omega[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_omegaz(int n)
{
    double **omega = atom->omega;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = omega[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_angmomx(int n)
{
    double **angmom = atom->angmom;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = angmom[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_angmomy(int n)
{
    double **angmom = atom->angmom;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = angmom[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_angmomz(int n)
{
    double **angmom = atom->angmom;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = angmom[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_tqx(int n)
{
    double **torque = atom->torque;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = torque[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_tqy(int n)
{
    double **torque = atom->torque;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = torque[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_tqz(int n)
{
    double **torque = atom->torque;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = torque[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_spin(int n)
{
    int *spin = atom->spin;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = spin[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_eradius(int n)
{
    double *eradius = atom->eradius;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = eradius[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_ervel(int n)
{
    double *ervel = atom->ervel;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = ervel[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_erforce(int n)
{
    double *erforce = atom->erforce;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = erforce[clist[i]];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_shapex(int n)
{
    double **shape = atom->shape;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = shape[clist[i]][0];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_shapey(int n)
{
    double **shape = atom->shape;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = shape[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_shapez(int n)
{
    double **shape = atom->shape;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = shape[clist[i]][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_blockiness1(int n)
{
    double **blockiness = atom->blockiness;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = blockiness[clist[i]][0];
        n += size_one;
    }
}
/* ---------------------------------------------------------------------- */

void DumpParticle::pack_blockiness2(int n)
{
    double **blockiness = atom->blockiness;

    for (int i = 0; i < nchoose; i++) {
        buf[n] = blockiness[clist[i]][1];
        n += size_one;
    }
}
//dump rotation matrix
void DumpParticle::pack_tensor(int n)
{
    double **quaternion = atom->quaternion;
    for (int i = 0; i < nchoose; i++) {
        double mat[3][3];
        MathExtra::quat_to_mat(quaternion[i], mat);
        buf[n]  =  mat[0][0];
        buf[n+1] = mat[0][1];
        buf[n+2] = mat[0][2];
        buf[n+3] = mat[1][0];
        buf[n+4] = mat[1][1];
        buf[n+5] = mat[1][2];
        buf[n+6] = mat[2][0];
        buf[n+7] = mat[2][1];
        buf[n+8] = mat[2][2];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_quat1(int n)
{
    double **quaternion = atom->quaternion;
    for (int i = 0; i < nchoose; i++) {
        buf[n] = quaternion[clist[i]][0];
        n += size_one;
    }
}
/* ---------------------------------------------------------------------- */

void DumpParticle::pack_quat2(int n)
{
    double **quaternion = atom->quaternion;
    for (int i = 0; i < nchoose; i++) {
        buf[n] = quaternion[clist[i]][1];
        n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpParticle::pack_quat3(int n)
{
    double **quaternion = atom->quaternion;
    for (int i = 0; i < nchoose; i++) {
        buf[n] = quaternion[clist[i]][2];
        n += size_one;
    }
}
/* ---------------------------------------------------------------------- */

void DumpParticle::pack_quat4(int n)
{
    double **quaternion = atom->quaternion;
    for (int i = 0; i < nchoose; i++) {
        buf[n] = quaternion[clist[i]][3];
        n += size_one;
    }
}

#endif // LAMMPS_VTK
