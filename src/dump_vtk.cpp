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
#include "update.h"
#include "universe.h"
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkRectilinearGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLPRectilinearGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPImageDataWriter.h>
#include <sstream>

namespace LAMMPS_NS
{

DumpVTK::DumpVTK(LAMMPS *lmp) :
    lmp_(lmp),
    vtk_compressor_(VTK_COMP_NONE),
    binary_(false)
{
    filesuffixes[0] = (char*) ".vtk";
    filesuffixes[1] = (char*) ".vtp";
    filesuffixes[2] = (char*) ".vtu";
    filesuffixes[3] = (char*) ".vti";
    filesuffixes[4] = (char*) ".vtr";
    filesuffixes[5] = (char*) ".vtm";
    filesuffixes[6] = (char*) ".pvtp";
    filesuffixes[7] = (char*) ".pvtu";
    filesuffixes[8] = (char*) ".pvti";
    filesuffixes[9] = (char*) ".pvtr";
}

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

void DumpVTK::write_vtp(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename)
{
    if (vtk_file_format == VTK_FILE_FORMATS::PVTP)
    {
        vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
        pwriter->SetFileName(filename);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(data);
#else
        pwriter->SetInputData(data);
#endif

        pwriter->SetNumberOfPieces(lmp_->comm->nprocs);
        pwriter->SetStartPiece(lmp_->comm->me);
        pwriter->SetEndPiece(lmp_->comm->me);
        pwriter->Write();
    }
    else if (vtk_file_format == VTK_FILE_FORMATS::VTP)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
        writer->SetInput(data);
#else
        writer->SetInputData(data);
#endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vtu(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename)
{
    if (vtk_file_format == VTK_FILE_FORMATS::PVTU)
    {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        pwriter->SetFileName(filename);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(data);
#else
        pwriter->SetInputData(data);
#endif

        pwriter->SetNumberOfPieces(lmp_->comm->nprocs);
        pwriter->SetStartPiece(lmp_->comm->me);
        pwriter->SetEndPiece(lmp_->comm->me);
        pwriter->Write();
    }
    else if (vtk_file_format == VTK_FILE_FORMATS::VTU)
    {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

#if VTK_MAJOR_VERSION < 6
        writer->SetInput(data);
#else
        writer->SetInputData(data);
#endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vti(vtkSmartPointer<vtkAlgorithmOutput> data, const int vtk_file_format, const char * const filename)
{
    if (vtk_file_format == VTK_FILE_FORMATS::PVTI)
    {
        vtkSmartPointer<vtkXMLPImageDataWriter> pwriter = vtkSmartPointer<vtkXMLPImageDataWriter>::New();
        pwriter->SetFileName(filename);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

        pwriter->SetInputConnection(data);

        pwriter->SetNumberOfPieces(lmp_->comm->nprocs);
        pwriter->SetStartPiece(lmp_->comm->me);
        pwriter->SetEndPiece(lmp_->comm->me);
        pwriter->SetWriteSummaryFile(lmp_->comm->me == 0 ? 1 : 0);
        pwriter->Write();
    }
    else if (vtk_file_format == VTK_FILE_FORMATS::VTI)
    {
        vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

        writer->SetInputConnection(data);

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vtr(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename)
{
    if (vtk_file_format == VTK_FILE_FORMATS::PVTR)
    {
        vtkSmartPointer<vtkXMLPRectilinearGridWriter> pwriter = vtkSmartPointer<vtkXMLPRectilinearGridWriter>::New();
        pwriter->SetFileName(filename);

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));

        #if VTK_MAJOR_VERSION <= 5
        pwriter->SetInputConnection(data->GetProducerPort());
        #else
        pwriter->SetInputData(data);
        #endif

        pwriter->SetNumberOfPieces(lmp_->comm->nprocs);
        pwriter->SetStartPiece(lmp_->comm->me);
        pwriter->SetEndPiece(lmp_->comm->me);
        pwriter->Write();
    }
    else if (vtk_file_format == VTK_FILE_FORMATS::VTR)
    {
        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(writer));

        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(data->GetProducerPort());
        #else
        writer->SetInputData(data);
        #endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vtk_unstructured_grid(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label)
{
    if (vtk_file_format == VTK_FILE_FORMATS::VTK)
    {
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

        if(label)
            writer->SetHeader(label);
        else
            writer->SetHeader("Generated by LIGGGHTS");
        setVtkWriterOptions(vtkDataWriter::SafeDownCast(writer));

        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(data->GetProducerPort());
        #else
        writer->SetInputData(data);
        #endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vtk_rectilinear_grid(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label)
{
    if (vtk_file_format == VTK_FILE_FORMATS::VTK)
    {
        vtkSmartPointer<vtkRectilinearGridWriter> writer = vtkSmartPointer<vtkRectilinearGridWriter>::New();

        if(label)
            writer->SetHeader(label);
        else
            writer->SetHeader("Generated by LIGGGHTS");
        setVtkWriterOptions(vtkDataWriter::SafeDownCast(writer));

        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(data->GetProducerPort());
        #else
        writer->SetInputData(data);
        #endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

void DumpVTK::write_vtk_poly(vtkSmartPointer<vtkDataObject> data, const int vtk_file_format, const char * const filename, char * const label)
{
    if (vtk_file_format == VTK_FILE_FORMATS::VTK)
    {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

        if(label)
            writer->SetHeader(label);
        else
            writer->SetHeader("Generated by LIGGGHTS");
        setVtkWriterOptions(vtkDataWriter::SafeDownCast(writer));

        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(data->GetProducerPort());
        #else
        writer->SetInputData(data);
        #endif

        writer->SetFileName(filename);
        writer->Write();
    }
    else
        lmp_->error->all(FLERR, "Internal error");
}

vtkMPIController * DumpVTK::getLocalController()
{
    vtkMPIController *vtkGlobalController = static_cast<vtkMPIController*>(vtkMultiProcessController::GetGlobalController());
    if (!vtkGlobalController)
        lmp_->error->all(FLERR, "Global VTK MPI Controller not found");

    if (lmp_->universe->existflag == 0)
        return vtkGlobalController;
    else
    {
        vtkMPIController *vtkLocalController = vtkGlobalController->PartitionController(lmp_->universe->iworld, 0);
        if (!vtkLocalController)
            lmp_->error->all(FLERR, "Local VTK MPI Controller not found");
        return vtkLocalController;
    }
}

void DumpVTK::setFileCurrent(char * &filecurrent, char * const filename, const int multifile, const int padflag)
{
    delete [] filecurrent;
    filecurrent = NULL;

    if (multifile == 0)
    {
        // contains no '*' -> simply copy
        filecurrent = new char[strlen(filename) + 1];
        strcpy(filecurrent, filename);
    }
    else
    {
        // contains '*' -> replace with time step
        filecurrent = new char[strlen(filename) + 16];
        char *ptr = strchr(filename,'*');
        *ptr = '\0';
        if (padflag == 0)
        {
            sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
                    filename,lmp_->update->ntimestep,ptr+1);
        }
        else
        {
            char bif[8],pad[16];
            strcpy(bif,BIGINT_FORMAT);
            sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
            sprintf(filecurrent,pad,filename,lmp_->update->ntimestep,ptr+1);
        }
        *ptr = '*';
    }
}

int DumpVTK::identify_file_type(char * const filename, std::list<int> &allowed_extensions, char * const style, int &multiproc, int &nclusterprocs, int &filewriter, int &fileproc, MPI_Comm &world, MPI_Comm &clustercomm)
{
    // ensure no old % format is used
    // this is set in dump.cpp
    if (multiproc)
        type_error("It is no longer allowed to enable parallel writing by setting the \% character, please see the documentation for help.", style, allowed_extensions);
    // find last dot
    char *suffix = strrchr(filename, '.');
    if (strlen(suffix) == 5)
    {
        multiproc = 1;
        nclusterprocs = 1;
        filewriter = 1;
        fileproc = lmp_->comm->me;
        MPI_Comm_split(world,lmp_->comm->me,0,&clustercomm);
        std::list<int>::iterator it = allowed_extensions.begin();
        for (; it != allowed_extensions.end(); it++)
        {
            if (*it >= VTK_FILE_FORMATS::vtk_serial_file_types && strcmp(suffix, filesuffixes[*it]) == 0)
                return *it;
        }
        type_error("Could not find allowed filetype for parallel writing.", style, allowed_extensions);
    }
    else if (strlen(suffix) == 4)
    {
        std::list<int>::iterator it = allowed_extensions.begin();
        for (; it != allowed_extensions.end(); it++)
        {
            if (*it < VTK_FILE_FORMATS::vtk_serial_file_types && strcmp(suffix, filesuffixes[*it]) == 0)
                return *it;
        }
        type_error("Could not find allowed filetype for serial writing.", style, allowed_extensions);
    }
    else
        type_error("Could not find allowed filetype for writing of VTK file.", style, allowed_extensions);
    return -1;
}

void DumpVTK::type_error(std::string msg, char * const style, std::list<int> &allowed_extensions)
{
    std::stringstream ss;
    ss << "dump " << std::string(style) << ": " << msg << " Allowed file extensions for this dump style are:";
    std::list<int>::iterator it = allowed_extensions.begin();
    for (; it != allowed_extensions.end(); it++)
        ss << " " << std::string(filesuffixes[*it]);
    lmp_->error->all(FLERR, ss.str().c_str());
}

}

#endif // LAMMPS_VTK
