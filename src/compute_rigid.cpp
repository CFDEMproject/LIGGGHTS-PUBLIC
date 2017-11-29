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

#include <cmath>
#include <string.h>
#include "compute_rigid.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "comm.h"
#include "container.h"
#include "fix_multisphere.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRigid::ComputeRigid(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg),
  is_single_(false),
  id_single_(-1),
  multisphere_(NULL),
  property_(NULL),
  len_(0)
{
  if (narg < iarg+2)
      error->compute_error(FLERR,this,"expecting at least 2 arguments");

  if(strcmp("all",group->names[igroup]))
      error->compute_error(FLERR,this,"must use group 'all'");

  if(strstr(style,"single"))
  {
      is_single_ = true;
      if(strcmp(arg[iarg++],"id"))
          error->compute_error(FLERR,this,"expecting 'id'");
      id_single_ = atoi(arg[iarg++]);
  }
  
  else
  {
      local_flag = 1;
      array_local = 0;
      vector_local = 0;
  }

  if(modify->n_fixes_style("multisphere") != 1)
      error->compute_error(FLERR,this,"defining exactly one fix multisphere is required");
  multisphere_ = & (static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0))->data());

  // parse property
  if(strcmp(arg[iarg++],"property"))
      error->compute_error(FLERR,this,"expecting keyword 'property'");
  property_ = multisphere_->prop().getElementPropertyBase(arg[iarg++]);

  // set sizes for single
  vector = 0;
  if(is_single_)
  {
      if(0 == property_->lenVec())
        error->compute_error(FLERR,this,"property has length of 0");
      
      else if(1 == property_->lenVec())
      {
          scalar_flag = 1;
      }
      
      else
      {
          vector_flag = 1;
          size_vector = property_->lenVec();
          vector = new double [size_vector];

          if(property_->isIntData())
            error->compute_error(FLERR,this,"int vectors currently not supported");
      }
  }

  // check arg validity
  if(!property_)
    error->compute_error(FLERR,this,"illegal property name used");

  update_pointers();
}

/* ---------------------------------------------------------------------- */

ComputeRigid::~ComputeRigid()
{
    if(vector) delete []vector;
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::update_pointers()
{
    
    if(is_single_)
    {
        int ibody = multisphere_->map(id_single_);
        if(vector_flag)
        {
            vectorZeroizeN(vector,size_vector);
            double **ptr = (double**) property_->begin_slow_dirty();
            if(ibody >= 0 && property_->isDoubleData())
                vectorAddN(vector,ptr[ibody],size_vector);
            MPI_Sum_Vector(vector,size_vector,world);
        }
        
        else
        {
            scalar = 0.;
            if(ibody >= 0 && property_->isDoubleData())
            {
                double *ptr = (double*) property_->begin_slow_dirty();
                scalar += ptr[ibody];
            }
            else if(property_->isIntData())
            {
                int *ptr = (int*) property_->begin_slow_dirty();
                scalar += static_cast<double>(ptr[ibody]);
            }
            MPI_Sum_Scalar(scalar,world);
        }
    }
    else
    {
        size_local_rows = multisphere_->n_body();

        if(property_->lenVec() > 1)
        {
            size_local_cols = property_->lenVec();

            if(len_ < size_local_rows)
            {
                len_ = size_local_rows;
                memory->grow(array_local,len_,size_local_cols,"compute/rigid:array_local");
            }

            if(property_->isDoubleData())
            {
                double **ptr = (double**) property_->begin_slow_dirty();
                for(int i = 0; i < size_local_rows; i++)
                    for(int j = 0; j < size_local_cols; j++)
                        array_local[i][j] = static_cast<double>(ptr[i][j]);
            }
            else if(property_->isIntData())
            {
                int **ptr = (int**) property_->begin_slow_dirty();
                for(int i = 0; i < size_local_rows; i++)
                    for(int j = 0; j < size_local_cols; j++)
                        array_local[i][j] = static_cast<double>(ptr[i][j]);
            }
        }
        else
        {
            size_local_cols  = 0;

            if(len_ < size_local_rows)
            {
                len_ = size_local_rows;
                memory->grow(vector_local,len_,"compute/rigid:vector_local");
            }

            if(property_->isDoubleData())
            {
                double *ptr = (double*) property_->begin_slow_dirty();
                for(int i = 0; i < size_local_rows; i++)
                    vector_local[i] = static_cast<double>(ptr[i]);
            }
            else if(property_->isIntData())
            {
                int *ptr = (int*) property_->begin_slow_dirty();
                for(int i = 0; i < size_local_rows; i++)
                    vector_local[i] = static_cast<double>(ptr[i]);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::init()
{
  
  update_pointers();
}

/* ---------------------------------------------------------------------- */

double ComputeRigid::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  update_pointers();

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::compute_vector()
{
  invoked_vector = update->ntimestep;

  update_pointers();
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::compute_local()
{
  invoked_local = update->ntimestep;

  update_pointers();
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeRigid::memory_usage()
{
  if(size_local_cols == 0)
    return size_local_rows;
  return size_local_rows*size_local_cols;
}
