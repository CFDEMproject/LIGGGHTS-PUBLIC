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

#include "dump_vtk.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

DumpVTK::DumpVTK(LAMMPS *lmp) :
    lmp_(lmp),
    vtk_compressor_(VTK_COMP_NONE),
    binary_(false)
{}

int DumpVTK::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"binary") == 0) {
        if (narg < 2)
            lmp_->error->all(FLERR,"Illegal dump_modify command [binary]");
        if (strcmp(arg[1],"yes") == 0)
            binary_ = true;
        else if (strcmp(arg[1],"no") == 0)
            binary_ = false;
        else
            lmp_->error->all(FLERR,"Illegal dump_modify command [binary]");
        return 2;
    }

    if (strcmp(arg[0],"compressor") == 0)
    {
        if (narg < 2)
            lmp_->error->all(FLERR,"Illegal dump_modify command [compressor]");

        if      (strcmp(arg[1],"zlib") == 0)
            vtk_compressor_ = VTK_COMP_ZLIB;
        else if (strcmp(arg[1],"lz4") == 0)
#if VTK_MAJOR_VERSION >= 8
            vtk_compressor_ = VTK_COMP_LZ4;
#else
            lmp_->error->all(FLERR, "Lz4 compressor is only available for VTK >= 8");
#endif
        else if (strcmp(arg[1],"none") == 0)
            vtk_compressor_ = VTK_COMP_NONE;
        else
            lmp_->error->all(FLERR,"Illegal dump_modify command [compressor]");

        // set binary on if compressor is used
        if (vtk_compressor_ != VTK_COMP_NONE && !binary_)
        {
            lmp_->error->warning(FLERR, "Vtk dump will switch to binary writing as compressor is used");
            binary_ = true;
        }
        return 2;
    }

    return 0;
}

void DumpVTK::setVtkWriterOptions(vtkSmartPointer<vtkDataWriter> writer)
{
    if (vtk_compressor_ != VTK_COMP_NONE && lmp_->comm->me == 0)
        lmp_->error->warning(FLERR, "Vtk compressor enabled but data format does not support compression. To avoid this message do not use the *.vtk file ending");

    if (binary_)
        writer->SetFileTypeToBinary();
    else
        writer->SetFileTypeToASCII();
}

void DumpVTK::setVtkWriterOptions(vtkSmartPointer<vtkXMLWriter> writer)
{
    if (binary_)
        writer->SetDataModeToBinary();
    else
        writer->SetDataModeToAscii();

    switch (vtk_compressor_)
    {
    case VTK_COMP_ZLIB:
        writer->SetCompressorTypeToZLib();
        break;
#if VTK_MAJOR_VERSION >= 8
    case VTK_COMP_LZ4:
        writer->SetCompressorTypeToLZ4();
        break;
#endif
    default:
        writer->SetCompressorTypeToNone();
        break;
    }
}

#endif // LAMMPS_VTK
