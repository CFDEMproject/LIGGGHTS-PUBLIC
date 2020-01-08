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

#ifdef LAMMPS_VTK

#include <string.h>
#include "dump_mesh_vtm.h"
#include "tri_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkMPIController.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMeshVTM::DumpMeshVTM(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    filecurrent(NULL),
    multiname_ex(NULL),
    dumpMesh(NULL)
{
    if (narg < 7)
        error->all(FLERR,"Illegal dump mesh/vtm command");

    //INFO: CURRENTLY ONLY PROC 0 writes

    format_default = NULL;

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPICommunicatorOpaqueComm vtkWorldOpaqueComm(&world);
        vtkMPICommunicator * vtkWorldComm = vtkMPICommunicator::New();
        vtkWorldComm->InitializeExternal(&vtkWorldOpaqueComm);
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->SetCommunicator(vtkWorldComm);
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
    vtkMPIController * controller = getLocalController();

    int ioptional = 5;
    dumpMesh = new DumpMesh(lmp, nclusterprocs, multiproc, filewriter, fileproc, controller);
    ioptional += dumpMesh->parse_parameters(narg-ioptional, &(arg[ioptional]));

    if (ioptional < narg)
        error->all(FLERR,"Invalid attribute in dump mesh/vtm command");

    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }
}

/* ---------------------------------------------------------------------- */

DumpMeshVTM::~DumpMeshVTM()
{
    if (filecurrent)
        delete [] filecurrent;
    delete dumpMesh;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::init_style()
{
    size_one = dumpMesh->init_style();
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write_header(bigint ndump)
{
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTM::count()
{
    return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::pack(int *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::setFileCurrent() {
    delete [] filecurrent;
    filecurrent = NULL;

    char *filestar = filename;
    if (multiproc)
    {
        if (multiproc > 1) // if dump_modify fileper or nfile was used
        {
            delete [] multiname_ex;
            multiname_ex = NULL;
            char *ptr = strchr(filename,'%');
            if (ptr)
            {
                int id;
                if (me + nclusterprocs == nprocs) // last filewriter
                    id = multiproc -1;
                else
                    id = me/nclusterprocs;
                multiname_ex = new char[strlen(filename) + 16];
                *ptr = '\0';
                sprintf(multiname_ex,"%s_%d%s",filename,id,ptr+1);
                *ptr = '%';
            }
        } // else multiname_ex built in constructor is OK
        filestar = multiname_ex;
    }

    if (multifile == 0)
    {
        filecurrent = new char[strlen(filestar) + 1];
        strcpy(filecurrent, filestar);
    }
    else
    {
        filecurrent = new char[strlen(filestar) + 16];
        char *ptr = strchr(filestar,'*');
        *ptr = '\0';
        if (padflag == 0)
        {
            sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
                    filestar,update->ntimestep,ptr+1);
        }
        else
        {
            char bif[8],pad[16];
            strcpy(bif,BIGINT_FORMAT);
            sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
            sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
        }
        *ptr = '*';
    }
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write()
{
    write_data(0, NULL);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write_data(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkMultiBlockDataSet> mbSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    dumpMesh->prepare_mbSet(mbSet);

    if (!filewriter)
        return;

    vtkSmartPointer<vtkXMLMultiBlockDataWriter> mbWriter = vtkXMLMultiBlockDataWriter::New();
    mbWriter->SetFileName(filecurrent);
    setVtkWriterOptions(vtkXMLWriter::SafeDownCast(mbWriter));
#if VTK_MAJOR_VERSION < 6
    mbWriter->SetInput(mbSet);
#else
    mbWriter->SetInputData(mbSet);
#endif
    mbWriter->Write();
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTM::modify_param(int narg, char **arg)
{
    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;

    return 0;
}

#endif // LAMMPS_VTK
