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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include "dump_decomposition_vtk.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"

// include last to ensure correct macros
#include "domain_definitions.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::DumpDecompositionVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5)
    error->all(FLERR,"Illegal dump decomposition command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc

  format_default = NULL;

  //number of properties written out in one line with buff
  size_one=1;  //dont use buff

  lasttimestep=-1;

  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;

  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::~DumpDecompositionVTK()
{
  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::init_style()
{
  if (domain->triclinic == 1)
    error->all(FLERR,"Can not perform dump decomposition for triclinic box");
  if (binary)
    error->all(FLERR,"Can not perform dump decomposition in binary mode");

  // default format not needed

  delete [] format;
  format = new char[150];

  // setup function ptrs

  header_choice = &DumpDecompositionVTK::header_item;
  pack_choice = &DumpDecompositionVTK::pack_item;
  write_choice = &DumpDecompositionVTK::write_item;

  // open single file, one time only

  if (multifile == 0) openfile();

  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;
  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump decomposition' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::count()
{
  if (comm->me!=0) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::pack(int *ids)
{
   (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::header_item(bigint ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/gran/VTK export\nASCII\n");
}

void DumpDecompositionVTK::footer_item()
{
  return;

}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::pack_item()
{
  
  xdata[0] = -BIG;
  if(comm->myloc[0] == 0) xdata[0] = domain->sublo[0];
  for(int i = 0; i < comm->procgrid[0]; i++)
  {
      xdata[i+1] = -BIG;
      if(comm->myloc[0] == i) xdata[i+1] = domain->subhi[0];
  }

  ydata[0] = -BIG;
  if(comm->myloc[1] == 0) ydata[0] = domain->sublo[1];
  for(int i = 0; i < comm->procgrid[1]; i++)
  {
      ydata[i+1] = -BIG;
      if(comm->myloc[1] == i) ydata[i+1] = domain->subhi[1];
  }

  zdata[0] = -BIG;
  if(comm->myloc[2] == 0) zdata[0] = domain->sublo[2];
  for(int i = 0; i < comm->procgrid[2]; i++)
  {
      zdata[i+1] = -BIG;
      if(comm->myloc[2] == i) zdata[i+1] = domain->subhi[2];
  }

  MPI_Allreduce(xdata,xdata_all,len[0],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(ydata,ydata_all,len[1],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(zdata,zdata_all,len[2],MPI_DOUBLE,MPI_MAX,world);

  return;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_item(int n, double *mybuf)
{
  
  if (comm->me!=0) return;

  //ensure it is only written once in multi-proc (work-around)
  if(lasttimestep==update->ntimestep)return;
  lasttimestep=update->ntimestep;

  //write the data
  fprintf(fp,"DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\n",len[0],len[1],len[2]);

  fprintf(fp,"X_COORDINATES %d float\n",len[0]);
  for (int i = 0; i < len[0]; i++)
     fprintf(fp,"%f ",xdata_all[i]);
  fprintf(fp,"\n");

  fprintf(fp,"Y_COORDINATES %d float\n",len[1]);
  for (int i = 0; i < len[1]; i++)
     fprintf(fp,"%f ",ydata_all[i]);
  fprintf(fp,"\n");

  fprintf(fp,"Z_COORDINATES %d float\n",len[2]);
  for (int i = 0; i < len[2]; i++)
     fprintf(fp,"%f ",zdata_all[i]);
  fprintf(fp,"\n");

  //footer not needed
}
