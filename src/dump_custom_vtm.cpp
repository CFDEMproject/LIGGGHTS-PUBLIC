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

#ifdef LAMMPS_VTK

#include <cmath>
#include "math_extra_liggghts.h"
#include <stdlib.h>
#include <string>
#include <cstddef>
#include "dump_custom_vtm.h"

#ifdef CONVEX_ACTIVE_FLAG
#include "atom_vec_convexhull.h"
#endif

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "region.h"
#include "update.h"
#include "variable.h"
#include <algorithm>
#include <vector>
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>
#include <vtkInformation.h>
#include <vtkMPIController.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>

#ifdef _WIN32
    #include "dirent.h"
#else
    #include <dirent.h>
#endif

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

/* ---------------------------------------------------------------------- */

DumpCustomVTM::DumpCustomVTM(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    nevery(0),
    filecurrent(NULL),
    parallelfilecurrent(NULL),
    multiname_ex(NULL),
    dumpParticle(NULL),
    dumpMesh(NULL),
    write_pv_config(false),
    write_pvd_file(false),
    pvd_ntimesteps(0)
{
    if (narg <= 7)
        error->all(FLERR,"dump custom/vtm lacks arguments");

    clearstep = 1;

    nevery = force->inumeric(FLERR,arg[3]);

    // TODO Remove the following lines
    // Single proc writing due to a paraview bug
    if (me != 0) filewriter = 0;
    fileproc = 0;
    multiproc = 0;
    nclusterprocs = nprocs;
    // TODO end remove

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

    bool hasargs = true;
    int iarg = 5;
    // all keywords that start a new dump class
    keyword_list.push_back("particle_properties");
    keyword_list.push_back("meshes");
    keyword_list.push_back("local_gran");
    keyword_list.push_back("write_pv_config");
    keyword_list.push_back("write_pvd_file");
    // loop over all arguments
    while (iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg], "particle_properties") == 0)
        {
            if (!dumpParticle)
                dumpParticle = new DumpParticle(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc);
            else
                error->all(FLERR, "Internal error (multiple particle_properties keywords?)");
            iarg += dumpParticle->parse_parameters(narg-iarg, &(arg[iarg]), false, keyword_list);
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "meshes") == 0)
        {
            if (!dumpMesh)
                dumpMesh = new DumpMesh(lmp, nclusterprocs, multiproc, filewriter, fileproc, controller);
            else
                error->all(FLERR, "Internal error (multiple meshes keywords?)");
            iarg += dumpMesh->parse_parameters(narg-iarg, &(arg[iarg]), keyword_list);
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "local_gran") == 0)
        {
            dumpLocalGranList.push_back(DumpLocalGran(lmp, igroup, nclusterprocs, multiproc, nevery, filewriter, fileproc));
            iarg += dumpLocalGranList.back().parse_parameters(narg-iarg, &(arg[iarg]), false, keyword_list);
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "write_pv_config") == 0)
        {
            if (++iarg >= narg)
                error->all(FLERR, "write_pv_config in dump custom/vtm requires one argument");
            if (strcmp(arg[iarg], "yes") == 0)
                write_pv_config = true;
            else if (strcmp(arg[iarg], "no") == 0)
                write_pv_config = false;
            else
                error->all(FLERR, "write_pv_config in dump custom/vtm requires either \"yes\" or \"no\" as argument");
            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "write_pvd_file") == 0)
        {
            if (++iarg >= narg)
                error->all(FLERR, "write_pvd_file in dump custom/vtm requires one argument");
            if (strcmp(arg[iarg], "yes") == 0)
                write_pvd_file = true;
            else if (strcmp(arg[iarg], "no") == 0)
                write_pvd_file = false;
            else
                error->all(FLERR, "write_pvd_file in dump custom/vtm requires either \"yes\" or \"no\" as argument");
            iarg++;
            hasargs = true;
        }
        else
        {
            error->all(FLERR, "Unknown keyword in dump custom/vtm");
        }

    }

    if (iarg < narg)
        error->all(FLERR,"Invalid attribute in dump custom/vtm command");

    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }

    for (int i = 0; i < output->ndump; i++)
    {
        if (strcmp(output->dump[i]->style, "custom/vtm") == 0)
            error->all(FLERR, "Only one dump custom/vtm is allowed");
    }

    if (write_pvd_file && comm->me == 0)
    {
        std::string fname(filename);

        DIR *dir;
        struct dirent *ent;
        std::size_t found = fname.find_last_of("/\\");
        std::string dirname = "";
        if (found != std::string::npos)
            dirname = fname.substr(0,found+1);
        //while (true)
        //{
        //    found = dirname.find("/");
        //    if (found == std::string::npos)
        //        break;
        //    dirname.replace(found, 1, "\\");
        //}
        if ((dir = opendir (dirname.c_str())) != NULL)
        {
            std::string basename = fname.substr(found+1);
            std::size_t last_spec_char = basename.find_last_of("*%");
            std::size_t first_spec_char = basename.find_first_of("*%");
            std::string prefix = basename.substr(0, first_spec_char);
            std::string postfix = basename.substr(last_spec_char+1);
            std::size_t star_char = basename.find("*");
            if (star_char == std::string::npos)
                error->one(FLERR, "Could not find * in filename");
            std::string starpostfix = basename.substr(star_char+1);
            // print all the files and directories within directory
            while ((ent = readdir (dir)) != NULL)
            {
                std::string curfile = std::string(ent->d_name);
                if (curfile.find(prefix) != 0)
                    continue;
                if (curfile.find(postfix) != curfile.length() - postfix.length())
                    continue;
                std::size_t end_nr = curfile.find(starpostfix);
                std::string nr = curfile.substr(star_char, end_nr - star_char);
                int ts = atoi(nr.c_str());
                // add each entry to a <int, std::string> pair list
                ts_files.push_back(std::pair<int, std::string>(ts, curfile));
            }
            closedir (dir);
            // sort the pair list by the int
            std::sort(ts_files.begin(), ts_files.end());
            // open pvd file
            std::fstream pvdFile;
            std::string pvdFname = dirname;
            pvdFname.append("liggghts_simulation.pvd");
            pvdFile.open(pvdFname.c_str(), std::fstream::out);
            // write header
            pvdFile << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
            pvdFile << "<Collection>" << endl;
            // write pvd file
            for (pvd_ntimesteps = 0; pvd_ntimesteps < ts_files.size(); pvd_ntimesteps++)
            {
                if (ts_files[pvd_ntimesteps].first >= update->ntimestep)
                    break;
                pvdFile << "<DataSet timestep=\"" << pvd_ntimesteps << "\" part=\"0\" file=\"" << ts_files[pvd_ntimesteps].second << "\"/>" << endl;
            }
            // delete unused entries
            if (pvd_ntimesteps < ts_files.size())
                ts_files.erase(ts_files.begin()+pvd_ntimesteps, ts_files.end());
            // write footer
            pvdFile << "</Collection>" << endl;
            pvdFile << "</VTKFile>" << endl;
            // close pvd file
            pvdFile.close();
        }
        else
        {
            /* could not open directory */
            error->one(FLERR, "Could not open directory in custom/vtm");
        }
    }
}

/* ---------------------------------------------------------------------- */

DumpCustomVTM::~DumpCustomVTM()
{
    delete [] filecurrent;
    delete [] parallelfilecurrent;
    delete [] multiname_ex;

    delete dumpParticle;

    if (dumpMesh)
        delete dumpMesh;

    dumpLocalGranList.clear();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTM::init_style()
{
    dumpParticle->init_style();

    if (dumpMesh)
        dumpMesh->init_style();

    std::list<DumpLocalGran>::iterator it;
    for (it = dumpLocalGranList.begin(); it != dumpLocalGranList.end(); it++)
        it->init_style();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTM::write_header(bigint /*ndump*/)
{
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTM::count()
{
    return 0;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTM::setFileCurrent() {
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

    // filename of parallel file
    if (multiproc)
    {
        delete [] parallelfilecurrent;
        parallelfilecurrent = NULL;

        // remove '%' character
        // -> string length stays the same
        // p is not added to filename as we are writing multiblock data
        char *ptr = strchr(filename,'%');
        filestar = new char[strlen(filename)];
        *ptr = '\0';
        sprintf(filestar,"%s%s",filename,ptr+1);
        *ptr = '%';
        ptr = strrchr(filestar,'.');
        ptr++;
        *ptr++='v';
        *ptr++='t';
        *ptr++='m';
        *ptr++= 0;

        if (multifile == 0)
        {
            parallelfilecurrent = new char[strlen(filestar) + 1];
            strcpy(parallelfilecurrent, filestar);
        }
        else
        {
            parallelfilecurrent = new char[strlen(filestar) + 16];
            char *ptr = strchr(filestar,'*');
            *ptr = '\0';
            if (padflag == 0)
                sprintf(parallelfilecurrent,"%s" BIGINT_FORMAT "%s",
                        filestar,update->ntimestep,ptr+1);
            else
            {
                char bif[8],pad[16];
                strcpy(bif,BIGINT_FORMAT);
                sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
                sprintf(parallelfilecurrent,pad,filestar,update->ntimestep,ptr+1);
            }
            *ptr = '*';
        }
        delete [] filestar;
        filestar = NULL;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTM::write()
{
    setFileCurrent();

    // multi block data set that will contain everything
    vtkSmartPointer<vtkMultiBlockDataSet> mbSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // mesh data
    if (dumpMesh)
        dumpMesh->prepare_mbSet(mbSet);

    std::list<DumpLocalGran>::iterator it;
    for (it = dumpLocalGranList.begin(); it != dumpLocalGranList.end(); it++)
        it->prepare_mbSet(mbSet, true);

    // particle & domain data
    
    dumpParticle->prepare_mbSet(mbSet);

    if (!filewriter)
        return;

    if (multiproc)
    {
        vtkSmartPointer<vtkXMLPMultiBlockDataWriter> pwriter = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(pwriter));
        pwriter->SetFileName(parallelfilecurrent);

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(mbSet);
#else
        pwriter->SetInputData(mbSet);
#endif

        pwriter->Write();
    }
    else if (comm->me == 0 || comm->nprocs == 1)
    {
        vtkSmartPointer<vtkXMLMultiBlockDataWriter> mbWriter = vtkXMLMultiBlockDataWriter::New();
        setVtkWriterOptions(vtkXMLWriter::SafeDownCast(mbWriter));
        mbWriter->SetFileName(filecurrent);
#if VTK_MAJOR_VERSION < 6
        mbWriter->SetInput(mbSet);
#else
        mbWriter->SetInputData(mbSet);
#endif
        mbWriter->Write();
    }

    // write out paraview configuration file for python script (paraview_generic_display.py)
    if (write_pv_config && comm->me == 0)
    {
        std::fstream configFile;
        std::string path = std::string(filecurrent);
        std::size_t endpath = path.find_last_of("/\\");
        std::string config_fname = "pv_config.txt";
        if (endpath != std::string::npos)
            config_fname = path.substr(0,endpath+1).append(config_fname);
        configFile.open(config_fname.c_str(), std::fstream::out);
        int nlocalgran = dumpLocalGranList.size();
        int nmeshes = mbSet->GetNumberOfBlocks() - nlocalgran - 2; // -2 = particles & domain
        configFile << "# This is a config file for paraview_generic_display.py" << endl;
        configFile << "# This file was created automatically using LIGGGHTS and dump custom/vtm" << endl;
        configFile << "########################" << endl;
        configFile << "# DO NOT EDIT THIS FILE" << endl;
        configFile << "########################" << endl;
        configFile << "nmeshes: " << nmeshes << endl;
        configFile << "nlocalgran: " << nlocalgran << endl;
        if (domain->triclinic)
        {
            configFile << "bbxmin: " << domain->boxlo_bound[0] << endl;
            configFile << "bbxmax: " << domain->boxhi_bound[0] << endl;
            configFile << "bbymin: " << domain->boxlo_bound[1] << endl;
            configFile << "bbymax: " << domain->boxhi_bound[1] << endl;
            configFile << "bbzmin: " << domain->boxlo_bound[2] << endl;
            configFile << "bbzmax: " << domain->boxhi_bound[2] << endl;
        }
        else
        {
            configFile << "bbxmin: " << domain->boxlo[0] << endl;
            configFile << "bbxmax: " << domain->boxhi[0] << endl;
            configFile << "bbymin: " << domain->boxlo[1] << endl;
            configFile << "bbymax: " << domain->boxhi[1] << endl;
            configFile << "bbzmin: " << domain->boxlo[2] << endl;
            configFile << "bbzmax: " << domain->boxhi[2] << endl;
        }
        configFile.close();
    }

    if (write_pvd_file && comm->me == 0)
    {
        std::string dirname = "";
        std::string fname(filename);
        std::size_t found = fname.find_last_of("/\\");
        if (found != std::string::npos)
            dirname = fname.substr(0,found+1);
        // add a new entry to a <int, std::string> pair list
        std::string curfile = multiproc ? parallelfilecurrent : filecurrent;
        found = curfile.find_last_of("/\\");
        if (found != std::string::npos)
            curfile = curfile.substr(found+1);
        ts_files.push_back(std::pair<int, std::string>(update->ntimestep, curfile));
        // open pvd file
        std::fstream pvdFile;
        std::string pvdFname = dirname;
        pvdFname.append("liggghts_simulation.pvd");
        pvdFile.open(pvdFname.c_str(), std::fstream::out);
        // write header
        pvdFile << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
        pvdFile << "<Collection>" << endl;
        // write pvd file
        for (unsigned int i = 0; i < ts_files.size(); i++)
        {
            pvdFile << "<DataSet timestep=\"" << i << "\" part=\"0\" file=\"" << ts_files[i].second << "\"/>" << endl;
        }
        // write footer
        pvdFile << "</Collection>" << endl;
        pvdFile << "</VTKFile>" << endl;
        // close pvd file
        pvdFile.close();
        pvd_ntimesteps++;
    }
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTM::modify_param(int narg, char **arg)
{
    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;

    return dumpParticle->modify_param(narg, arg);
}

/* ----------------------------------------------------------------------
     return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpCustomVTM::memory_usage()
{
    bigint bytes = Dump::memory_usage();
    bytes += dumpParticle->memory_usage();
    if (dumpMesh)
        bytes += dumpMesh->memory_usage();
    std::list<DumpLocalGran>::iterator it;
    for (it = dumpLocalGranList.begin(); it != dumpLocalGranList.end(); it++)
        bytes += it->memory_usage();
    return bytes;
}

#endif
