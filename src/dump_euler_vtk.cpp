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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "string.h"
#include "dump_euler_vtk.h"
#include "fix_ave_euler.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpEulerVTK::DumpEulerVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  fix_euler_(0),
  n_calls_(0),
  n_all_(0),
  n_all_max_(0),
  buf_all_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump pic/vtk command");

  // CURRENTLY ONLY PROC 0 writes

  format_default = NULL;
}

/* ---------------------------------------------------------------------- */

DumpEulerVTK::~DumpEulerVTK()
{
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::init_style()
{
  fix_euler_ = static_cast<FixAveEuler*>(modify->find_fix_style("ave/euler",0));
  if(!fix_euler_)
    error->all(FLERR,"Illegal dump euler/vtk command, need a fix ave/euler");

  // multifile=1;             // 0 = one big file, 1 = one file per timestep
  // multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multifile != 1)
    error->all(FLERR,"You should use a filename like 'dump*.vtk' for the 'dump euler/vtk' command to produce one file per time-step");
  if (multiproc != 0)
    error->all(FLERR,"Your 'dump euler/vtk' command is writing one file per processor, where all the files contain the same data");

//  if (domain->triclinic == 1)
//    error->all(FLERR,"Can not dump VTK files for triclinic box");
  if (binary)
    error->all(FLERR,"Can not dump VTK files in binary mode");

  // node center (3), av vel (3), volume fraction, stress, radius
  size_one = 9;

  delete [] format;
}

/* ---------------------------------------------------------------------- */

int DumpEulerVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump euler/vtk' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::write_header(bigint ndump)
{
  write_header_ascii(ndump);
}

void DumpEulerVTK::write_header_ascii(bigint ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/VTK export\nASCII\n");
}

/* ---------------------------------------------------------------------- */

int DumpEulerVTK::count()
{
  n_calls_ = 0;
  n_all_ = 0;
  return fix_euler_->ncells_pack();
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::pack(int *ids)
{
  int m = 0;

  // have to stick with this order (all per-element props)
  // as multiple procs pack

  int ncells = fix_euler_->ncells_pack();

  for(int i = 0; i < ncells; i++)
  {
    buf[m++] = fix_euler_->cell_center(i,0);
    buf[m++] = fix_euler_->cell_center(i,1);
    buf[m++] = fix_euler_->cell_center(i,2);

    buf[m++] = fix_euler_->cell_v_av(i,0);
    buf[m++] = fix_euler_->cell_v_av(i,1);
    buf[m++] = fix_euler_->cell_v_av(i,2);

    buf[m++] = fix_euler_->cell_vol_fr(i);
    buf[m++] = fix_euler_->cell_radius(i);
    buf[m++] = fix_euler_->cell_pressure(i);
  }
  return ;
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::write_data(int n, double *mybuf)
{
    //only proc 0 writes
    if (comm->me != 0) return;

    n_calls_++;

    // grow buffer if necessary
    if(n_all_+n*size_one > n_all_max_)
    {
        n_all_max_ = n_all_ + n*size_one;
        memory->grow(buf_all_,n_all_max_,"DumpEulerVTK:buf_all_");
    }

    // copy to buffer
    vectorCopyN(mybuf,&(buf_all_[n_all_]),n*size_one);
    n_all_ += n*size_one;

    // write on last call
    if(n_calls_ == comm->nprocs)
        write_data_ascii(n_all_/size_one,buf_all_);
}

void DumpEulerVTK::write_data_ascii(int n, double *mybuf)
{

  int m, buf_pos;

  // n is the number of elements

  // write point data
  fprintf(fp,"DATASET POLYDATA\nPOINTS %d float\n",n);
  m = 0;
  buf_pos = 0;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += size_one ;
  }
  buf_pos += 3;

  // write polygon data
  fprintf(fp,"VERTICES %d %d\n",n,2*n);
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%d %d\n",1,i);
  }

  // write point data header
  fprintf(fp,"POINT_DATA %d\n",n);

  // write cell data

  fprintf(fp,"VECTORS v_avg float\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
     fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
     m += size_one;
  }
  buf_pos += 3;

  fprintf(fp,"SCALARS volumefraction float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  fprintf(fp,"SCALARS radius float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  fprintf(fp,"SCALARS pressure float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  // footer not needed
  // if would be needed, would do like in dump stl
}
