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

    Copyright 2017-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#ifndef LMP_DUMP_VTK_H
#define LMP_DUMP_VTK_H

#include "lammps.h"
#include <vtkSmartPointer.h>
#include <vtkXMLWriter.h>
#include <vtkDataWriter.h>
#include <vtkDataObject.h>
#include <vtkAlgorithmOutput.h>
#include <vtkMPIController.h>
#include <list>
#include <string>

namespace LAMMPS_NS
{

namespace VTK_FILE_FORMATS
{
// file formats
// serial need to come first
enum
{
    VTK,
    VTP,
    VTU,
    VTI,
    VTR,
    VTM,
    PVTP,
    PVTU,
    PVTI,
    PVTR,
    VTK_INVALID
};

// number of serial vtk file types
const int vtk_serial_file_types = 6;
}; // namespace VTK_FILE_FORMATS

class DumpVTK
{
public:
    DumpVTK(LAMMPS *lmp);

    int modify_param(int narg, char **arg);

    void setVtkWriterOptions(vtkSmartPointer<vtkXMLWriter> writer);
    void setVtkWriterOptions(vtkSmartPointer<vtkDataWriter> writer);

    void write_vtp(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename);
    void write_vtu(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename);
    void write_vti(vtkSmartPointer<vtkAlgorithmOutput> data, const int vtk_file_format, const char * const filename);
    void write_vtr(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename);

    void write_vtk_poly(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label = NULL);
    void write_vtk_unstructured_grid(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label = NULL);
    void write_vtk_rectilinear_grid(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label = NULL);

    vtkMPIController *getLocalController();

    void setFileCurrent(char * &filecurrent, char * const filename, const int multifile, const int padflag);
    int identify_file_type(char * const filename, std::list<int> &allowed_extensions, char * const style, int &multiproc, int &nclusterprocs, int &filewriter, int &fileproc, MPI_Comm &world, MPI_Comm &clustercomm);

private:

    void type_error(std::string msg, char * const style, std::list<int> &allowed_extensions);

    // compressors
    enum
    {
        VTK_COMP_ZLIB,
#if VTK_MAJOR_VERSION >= 8
        VTK_COMP_LZ4,
#endif
        VTK_COMP_NONE
    };

    LAMMPS *lmp_;
    int vtk_compressor_;
    bool binary_;

    char * filesuffixes[VTK_FILE_FORMATS::VTK_INVALID];
};

}; // namespace

#endif // LMP_DUMP_VTK_H

#endif // LAMMPS_VTK
