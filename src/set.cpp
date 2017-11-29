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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include "math_extra.h"
#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"
#include "atom_vec_superquadric.h"
#endif
#include <stdlib.h>
#include <string.h>
#include "set.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "math_extra.h"
#include "fix_multisphere.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "modify.h" 
#include "fix_property_atom.h" 
#include "sph_kernels.h" 
#include "fix_sph.h" 
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{ATOM_SELECT,MOL_SELECT,TYPE_SELECT,GROUP_SELECT,REGION_SELECT};
enum{TYPE,TYPE_FRACTION,MOLECULE,X,Y,Z,CHARGE,MASS,SHAPE,LENGTH,
     DIPOLE,DIPOLE_RANDOM,QUAT,QUAT_RANDOM, QUAT_DIRECT, THETA,ANGMOM,
     DIAMETER,DENSITY,VOLUME,IMAGE,BOND,ANGLE,DIHEDRAL,IMPROPER,
     MESO_E,MESO_CV,MESO_RHO,INAME,DNAME,
     VX,VY,VZ,OMEGAX,OMEGAY,OMEGAZ,PROPERTYPERATOM,BLOCKINESS,
     ASPECTRATIO, INERTIAX, INERTIAY, INERTIAZ,
     SMD_MASS_DENSITY, SMD_CONTACT_RADIUS}; 

#define BIG INT_MAX

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Set command with no atoms existing");
  if (narg < 3) error->all(FLERR,"Illegal set command");

  int n_ms = modify->n_fixes_style("multisphere");
  if(n_ms > 0 && !static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0))->allow_group_and_set())
    error->all(FLERR,"Set command may not be used together with fix multisphere");

  // style and ID info

  if (strcmp(arg[0],"atom") == 0) style = ATOM_SELECT;
  else if (strcmp(arg[0],"mol") == 0) style = MOL_SELECT;
  else if (strcmp(arg[0],"type") == 0) style = TYPE_SELECT;
  else if (strcmp(arg[0],"group") == 0) style = GROUP_SELECT;
  else if (strcmp(arg[0],"region") == 0) style = REGION_SELECT;
  else error->all(FLERR,"Illegal set command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id,arg[1]);
  select = NULL;
  selection(atom->nlocal);

  add    = 0; 
  until  = 1; 

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen,"Setting atom values ...\n");

  int allcount,origarg;

  int iarg = 2;
  while (iarg < narg) {
    varflag = varflag1 = varflag2 = varflag3 = varflag4 = 0;
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      set(TYPE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"type/fraction") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      newtype = force->inumeric(FLERR,arg[iarg+1]);
      fraction = force->numeric(FLERR,arg[iarg+2]);
      ivalue = force->inumeric(FLERR,arg[iarg+3]);
      if (newtype <= 0 || newtype > atom->ntypes)
        error->all(FLERR,"Invalid value in set command");
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR,"Invalid value in set command");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      setrandom(TYPE_FRACTION);
      iarg += 4;

    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (!atom->molecule_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MOLECULE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(X);
      iarg += 2;

    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(Y);
      iarg += 2;

    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(Z);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vx") == 0) {  
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VX);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VZ);
      iarg += 2;
    } else if (strcmp(arg[iarg],"omegax") == 0) { 
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(OMEGAX);
      iarg += 2;
    } else if (strcmp(arg[iarg],"omegay") == 0) {  
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(OMEGAY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {  
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(OMEGAZ);
      iarg += 2;

    } else if (strcmp(arg[iarg],"inertiax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if(!atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(INERTIAX);
      iarg += 2;
    } else if (strcmp(arg[iarg],"inertiay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if(!atom->superquadric_flag)
          error->all(FLERR,"Cannot set this attribute for this atom style");
      set(INERTIAY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"inertiaz") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if(!atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(INERTIAZ);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->q_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(CHARGE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rmass_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MASS);
      iarg += 2;

    } else if (strcmp(arg[iarg],"shape") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->ellipsoid_flag && !atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(SHAPE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"blockiness") == 0 or strcmp(arg[iarg],"roundness") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (!atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if(strcmp(arg[iarg],"roundness") == 0)
        error->warning(FLERR,"Keyword 'roundness' will be deprecated in future, please use blockiness istead");
      set(BLOCKINESS);
      iarg += 3;
    } else if (strcmp(arg[iarg],"aspect_ratio") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (!atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(ASPECTRATIO);
      iarg += 3;
    } else if (strcmp(arg[iarg],"length") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->line_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(LENGTH);
      iarg += 2;

    } else if (strcmp(arg[iarg],"dipole") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->mu_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(DIPOLE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"dipole/random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      dvalue = force->numeric(FLERR,arg[iarg+2]);
      if (!atom->mu_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      if (dvalue <= 0.0)
        error->all(FLERR,"Invalid dipole length in set command");
      setrandom(DIPOLE_RANDOM);
      iarg += 3;

    } else if (strcmp(arg[iarg],"quat") == 0 || strcmp(arg[iarg],"quat_direct") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) varparse(arg[iarg+4],4);
      else wvalue = force->numeric(FLERR,arg[iarg+4]);
      if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if(strcmp(arg[iarg],"quat") == 0)
      set(QUAT);
      else
        set(QUAT_DIRECT);
      iarg += 5;

    } else if (strcmp(arg[iarg],"quat/random") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      setrandom(QUAT_RANDOM);
      iarg += 2;

    } else if (strcmp(arg[iarg],"theta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        dvalue = force->numeric(FLERR,arg[iarg+1]);
      dvalue *= MY_PI/180.0;
      }
      if (!atom->line_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(THETA);
      iarg += 2;

    } else if (strcmp(arg[iarg],"angmom") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->ellipsoid_flag && !atom->tri_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(ANGMOM);
      iarg += 4;

    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->radius_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (dvalue < 0.0) error->all(FLERR,"Invalid diameter in set command");
      set(DIAMETER);
      iarg += 2;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rmass_flag && !atom->density_flag) 
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (dvalue <= 0.0) error->all(FLERR,"Invalid density in set command");
      set(DENSITY);
      iarg += 2;

    } else if (strcmp(arg[iarg],"volume") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->vfrac_flag && !atom->superquadric_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (dvalue <= 0.0) error->all(FLERR,"Invalid volume in set command");
      set(VOLUME);
      iarg += 2;

    } else if (strcmp(arg[iarg],"image") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1],"NULL") != 0) {
        ximageflag = 1;
        ximage = force->inumeric(FLERR,arg[iarg+1]);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0) {
        yimageflag = 1;
        yimage = force->inumeric(FLERR,arg[iarg+2]);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0) {
        zimageflag = 1;
        zimage = force->inumeric(FLERR,arg[iarg+3]);
      }
      if (ximageflag && ximage && !domain->xperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (yimageflag && yimage && !domain->yperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (zimageflag && zimage && !domain->zperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      set(IMAGE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->bonds_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nbondtypes)
        error->all(FLERR,"Invalid value in set command");
      topology(BOND);
      iarg += 2;

    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nangletypes)
        error->all(FLERR,"Invalid value in set command");
      topology(ANGLE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->ndihedraltypes)
        error->all(FLERR,"Invalid value in set command");
      topology(DIHEDRAL);
      iarg += 2;

    } else if (strcmp(arg[iarg],"improper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->impropers_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nimpropertypes)
        error->all(FLERR,"Invalid value in set command");
      topology(IMPROPER);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso_e") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->e_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MESO_E);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso_cv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->cv_flag)
            error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MESO_CV);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso_rho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rho_flag)
        error->all(FLERR,"Cannot set meso_rho for this atom style");
      set(MESO_RHO);
      iarg += 2;
    } else if (strcmp(arg[iarg],"smd/mass/density") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
          if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
          else dvalue = force->numeric(FLERR,arg[iarg+1]);
          if (!atom->smd_flag)
            error->all(FLERR,"Cannot set smd/mass/density for this atom style");
          set(SMD_MASS_DENSITY);
          iarg += 2;

    } else if (strcmp(arg[iarg],"smd/contact/radius") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
          if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
          else dvalue = force->numeric(FLERR,arg[iarg+1]);
          if (!atom->smd_flag)
            error->all(FLERR,"Cannot set smd/contact/radius "
                       "for this atom style");
          set(SMD_CONTACT_RADIUS);
          iarg += 2;

    } else if (strstr(arg[iarg],"i_") == arg[iarg]) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      int flag;
      index_custom = atom->find_custom(&arg[iarg][2],flag);
      if (index_custom < 0 || flag != 0)
        error->all(FLERR,"Set command integer vector does not exist");
      set(INAME);
      iarg += 2;

    } else if (strstr(arg[iarg],"d_") == arg[iarg]) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      int flag;
      index_custom = atom->find_custom(&arg[iarg][2],flag);
      if (index_custom < 0 || flag != 1)
        error->all(FLERR,"Set command floating point vector does not exist");
      set(DNAME);
      iarg += 2;
    } else if (strcmp(arg[iarg],"add") == 0){  
      if (iarg+1 > narg)
        error->all(FLERR,"Illegal set command for add");
      if(strcmp(arg[iarg+1],"yes") == 0)
      add = 1;
      else if(strcmp(arg[iarg+1],"no") == 0)
      add = 0;
      else error->all(FLERR,"Illegal 'add' option called");
      iarg +=2;
    } else if (strcmp(arg[iarg],"until") == 0){ 
     if (iarg+1 > narg)
        error->all(FLERR,"Illegal set command for until");
      until = atof(arg[iarg+1]);
      if (until <= 0.0)
        error->all(FLERR,"Illegal 'until' option called. Please set keyword value >0");
      iarg +=2;
    } else if (strncmp(arg[iarg],"property/atom",13) == 0) { 
      if (iarg+1 > narg)
        error->all(FLERR,"Illegal set command for property/atom");
      //find the fix (there should be only one fix with the same variablename, this is ensured by the fix itself)
      int n = strlen(arg[iarg+1]) + 1;
      char* variablename = new char[n];
      strcpy(variablename,arg[iarg+1]);
      updFix = NULL;
      for (int ifix = 0; ifix < (lmp->modify->nfix); ifix++){
        if ((strncmp(modify->fix[ifix]->style,"property/atom",13) == 0) && (strcmp(((FixPropertyAtom*)(modify->fix[ifix]))->variablename,variablename)==0) ){
            updFix=(FixPropertyAtom*)(lmp->modify->fix[ifix]);
        }
      }
      delete []variablename;
      if (updFix==NULL)
        error->all(FLERR,"Could not identify the per-atom property you want to set");
      nUpdValues=updFix->nvalues;

      if (nUpdValues != (narg-iarg-2) )
        error->all(FLERR,"The number of values for the set property/atom does not match the number needed");

      //get to update values
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2])
      {
        if (nUpdValues != 1 && nUpdValues != 3)
            error->all(FLERR,"Set command for property/atom and variables does only work with scalar and 3-vector quantities");
        varparse(arg[iarg+2],1);
        if(nUpdValues == 3)
        {
            varparse(arg[iarg+3],2);
            varparse(arg[iarg+4],3);
        }
        updValues = 0;
      }
      else
      {
        updValues = new double[nUpdValues];
        for(int j=0;j<nUpdValues ;j++)
          updValues[j]=atof(arg[iarg+1+1+j]);
      }

      set(PROPERTYPERATOM);
      if(updValues) delete []updValues;
      iarg += (2+nUpdValues);
    } else if (strcmp(arg[iarg],"sphkernel") == 0) { 
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");

      // check uniqueness of kernel IDs

      int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
      if(flag < 0) error->all(FLERR, "Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

      // get kernel id

      dvalue = SPH_KERNEL_NS::sph_kernel_id(arg[iarg+1]);
      if(dvalue < 0) error->all(FLERR, "Illegal pair_style sph command, unknown sph kernel");

      // set kernel_id in all sph fixes

      if (comm->me == 0 && screen) {
        fprintf(screen,"Setting undefined fix_sph kernel IDs ...\n");
        fprintf(screen,"  Sph styles with undefined kernel_id found: \n");
      }
      for (int ifix = 0; ifix < modify->nfix; ifix++)
      {
        if (strstr(modify->fix[ifix]->style,"sph") && dynamic_cast<FixSph*>(modify->fix[ifix]) ) {
          if (((FixSph *)(modify->fix[ifix]))->kernel_flag && ((FixSph *)(modify->fix[ifix]))->get_kernel_id() < 0) {
            if (comm->me == 0 && screen) fprintf(screen,"  Fix style = %s\n",modify->fix[ifix]->style);
            ((FixSph *)(modify->fix[ifix]))->set_kernel_id(dvalue);
            count++;
          }
        }
      }

      iarg += 2;
    } else error->all(FLERR,"Illegal set command");

    // statistics

    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);

    if (comm->me == 0) {
      if (screen) fprintf(screen,"  %d settings made for %s\n",
                          allcount,arg[origarg]);
      if (logfile) fprintf(logfile,"  %d settings made for %s\n",
                           allcount,arg[origarg]);
    }
  }

  // free local memory

  delete [] id;
  delete [] select;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, MOLECULE, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set::selection(int n)
{
  delete [] select;
  select = new int[n];
  int nlo,nhi;

  if (style == ATOM_SELECT) {
    if (atom->tag_enable == 0)
      error->all(FLERR,"Cannot use set atom with no atom IDs defined");
    force->bounds(id,BIG,nlo,nhi);

    int *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlo && tag[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == MOL_SELECT) {
    if (atom->molecule_flag == 0)
      error->all(FLERR,"Cannot use set mol with no molecule IDs defined");
    else force->bounds(id,BIG,nlo,nhi,0);

    int *molecule = atom->molecule;
    for (int i = 0; i < n; i++)
      if (molecule[i] >= nlo && molecule[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    force->bounds(id,atom->ntypes,nlo,nhi);

    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR,"Could not find set group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all(FLERR,"Set region ID does not exist");

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
        select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set owned atom properties directly
   either scalar or per-atom values from atom-style variable(s)
------------------------------------------------------------------------- */

void Set::set(int keyword)
{
  // evaluate atom-style variable(s) if necessary

  vec1 = vec2 = vec3 = vec4 = NULL;

  if (varflag) {
    int nlocal = atom->nlocal;
    if (varflag1) {
      memory->create(vec1,nlocal,"set:vec1");
      input->variable->compute_atom(ivar1,0,vec1,1,0);
    }
    if (varflag2) {
      memory->create(vec2,nlocal,"set:vec2");
      input->variable->compute_atom(ivar2,0,vec2,1,0);
    }
    if (varflag3) {
      memory->create(vec3,nlocal,"set:vec3");
      input->variable->compute_atom(ivar3,0,vec3,1,0);
    }
    if (varflag4) {
      memory->create(vec4,nlocal,"set:vec4");
      input->variable->compute_atom(ivar4,0,vec4,1,0);
    }
  }

  // loop over selected atoms

  AtomVecEllipsoid *avec_ellipsoid =
    (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecLine *avec_line = (AtomVecLine *) atom->style_match("line");
  #ifdef SUPERQUADRIC_ACTIVE_FLAG
  AtomVecSuperquadric *avec_superquadric = (AtomVecSuperquadric *) atom->style_match("superquadric");
  #endif

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    // overwrite dvalue, ivalue, xyzw value if variables defined
    // else the input script scalar value remains in place

    if (varflag1) {
      dvalue = xvalue = vec1[i];
      ivalue = static_cast<int> (dvalue);
    }
    if (varflag2) yvalue = vec2[i];
    if (varflag3) zvalue = vec3[i];
    if (varflag4) wvalue = vec4[i];

    // set values in per-atom arrays
    // error check here in case atom-style variables generated bogus value

    if (keyword == TYPE) {
      if (ivalue <= 0 || ivalue > atom->ntypes)
        error->one(FLERR,"Invalid value in set command");
      atom->type[i] = ivalue;
    }
    else if (keyword == MOLECULE) atom->molecule[i] = ivalue;
    else if (keyword == X) atom->x[i][0] = dvalue;
    else if (keyword == Y) atom->x[i][1] = dvalue;
    else if (keyword == Z) atom->x[i][2] = dvalue;
    else if (keyword == VX) atom->v[i][0] = dvalue; 
    else if (keyword == VY) atom->v[i][1] = dvalue;
    else if (keyword == VZ) atom->v[i][2] = dvalue;
    #ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (keyword == INERTIAX) atom->inertia[i][0] = dvalue;
    else if (keyword == INERTIAY) atom->inertia[i][1] = dvalue;
    else if (keyword == INERTIAZ) atom->inertia[i][2] = dvalue;
    else if (keyword == OMEGAX && atom->superquadric_flag) {
        atom->omega[i][0] = dvalue;
        MathExtraLiggghtsNonspherical::omega_to_angmom(atom->quaternion[i], atom->omega[i], atom->inertia[i],atom->angmom[i]);
    }
    #endif
    else if (keyword == OMEGAX) atom->omega[i][0] = dvalue;  
    #ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (keyword == OMEGAY && atom->superquadric_flag) {
        atom->omega[i][1] = dvalue;
        MathExtraLiggghtsNonspherical::omega_to_angmom(atom->quaternion[i], atom->omega[i], atom->inertia[i],atom->angmom[i]);
    }
    #endif
    else if (keyword == OMEGAY) atom->omega[i][1] = dvalue;  
    #ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (keyword == OMEGAZ && atom->superquadric_flag) {
        atom->omega[i][2] = dvalue;
        MathExtraLiggghtsNonspherical::omega_to_angmom(atom->quaternion[i], atom->omega[i], atom->inertia[i],atom->angmom[i]);
    }
    #endif
    else if (keyword == OMEGAZ) atom->omega[i][2] = dvalue;  
    else if (keyword == CHARGE) atom->q[i] = dvalue;
    else if (keyword == MASS) {
              if (dvalue <= 0.0) error->one(FLERR,"Invalid mass in set command");
        atom->rmass[i] = dvalue;
    }
    else if (keyword == DIAMETER) {
       if (dvalue < 0.0) error->one(FLERR,"Invalid diameter in set command");
        atom->radius[i] = 0.5 * dvalue * force->cg(atom->type[i]);
        
        if(atom->rmass_flag && atom->density_flag && atom->density[i] > 0.)
        {
          if(atom->superquadric_flag) {
            double vol = MY_PI/6.0 * dvalue * dvalue * dvalue;
            atom->volume[i] = vol;
          }
          else {
          if (domain->dimension == 2)
            atom->rmass[i] = MY_PI * atom->radius[i]*atom->radius[i] * atom->density[i];
          else
            atom->rmass[i] = 4.0*MY_PI/3.0 * atom->radius[i]*atom->radius[i]*atom->radius[i] * atom->density[i];
          }
        }
    }
    else if (keyword == VOLUME) {
      if (dvalue <= 0.0) error->one(FLERR,"Invalid volume in set command");
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if (avec_superquadric)
        atom->volume[i] = dvalue;
      else
        atom->vfrac[i] = dvalue;
#else
        atom->vfrac[i] = dvalue;
#endif
    }
    else if (keyword == MESO_E) atom->e[i] = dvalue;
    else if (keyword == MESO_CV) atom->cv[i] = dvalue;
    else if (keyword == MESO_RHO) atom->rho[i] = dvalue;

    else if (keyword == SMD_MASS_DENSITY) { 
      // set mass from volume and supplied mass density
      atom->rmass[i] = atom->vfrac[i] * dvalue;
    }
    else if (keyword == SMD_CONTACT_RADIUS) atom->contact_radius[i] = dvalue;
    // set shape of ellipsoidal particle

    else if (keyword == SHAPE) {
      if (xvalue < 0.0 || yvalue < 0.0 || zvalue < 0.0)
        error->one(FLERR,"Invalid shape in set command");
      if (xvalue > 0.0 || yvalue > 0.0 || zvalue > 0.0) {
        if (xvalue == 0.0 || yvalue == 0.0 || zvalue == 0.0)
          error->one(FLERR,"Invalid shape in set command");
      }
      if(avec_ellipsoid)
        avec_ellipsoid->set_shape(i,0.5*xvalue,0.5*yvalue,0.5*zvalue);
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      else if(avec_superquadric) {
        atom->shape[i][0] = xvalue;
        atom->shape[i][1] = yvalue;
        atom->shape[i][2] = zvalue;
        MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(atom->shape[i], atom->blockiness[i], atom->radius+i); //re-calculate bounding sphere radius
        MathExtraLiggghtsNonspherical::volume_superquadric(atom->shape[i], atom->blockiness[i], atom->volume+i); //re-calculate volume
        MathExtraLiggghtsNonspherical::area_superquadric(atom->shape[i], atom->blockiness[i], atom->area+i);  //re-calculate surface area
        atom->rmass[i] = atom->volume[i] * atom->density[i]; //re-calculate mass
        MathExtraLiggghtsNonspherical::inertia_superquadric(atom->shape[i], atom->blockiness[i], atom->density[i], atom->inertia[i]); //re-calculate inertia tensor

      }
      #endif // SUPERQUADRIC_ACTIVE_FLAG
      else
        error->one(FLERR,"Cannot set shape for this type of atom");
    }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (keyword == ASPECTRATIO) {
      if(avec_superquadric) {
        if (xvalue < 0.0 || yvalue < 0.0)
          error->one(FLERR,"Invalid aspect ratio in set command");
        if (xvalue > 0.0 || yvalue > 0.0) {
          if (xvalue == 0.0 || yvalue == 0.0)
            error->one(FLERR,"Invalid aspect ratio in set command");
        }
        double k1 = xvalue;
        double k2 = yvalue;

        double shape_[3]={1.0, 1.0, 1.0};
        double f_;
        double vol = atom->volume[i];
        MathExtraLiggghtsNonspherical::volume_superquadric(shape_, atom->blockiness[i], &f_);
        atom->shape[i][0] = cbrt(atom->volume[i] / (k1*k2*f_));
        atom->shape[i][1] = k1*atom->shape[i][0];
        atom->shape[i][2] = k2*atom->shape[i][0];

        MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(atom->shape[i], atom->blockiness[i], atom->radius+i); //re-calculate bounding sphere radius
        MathExtraLiggghtsNonspherical::area_superquadric(atom->shape[i], atom->blockiness[i], atom->area+i);  //re-calculate surface area
        MathExtraLiggghtsNonspherical::inertia_superquadric(atom->shape[i], atom->blockiness[i], atom->density[i], atom->inertia[i]); //re-calculate inertia tensor
      }
     else
      error->one(FLERR,"Cannot set shape for this type of atom");
    }
#endif // SUPERQUADRIC_ACTIVE_FLAG
   //set roundness parameters for superquadric *
    else if (keyword == BLOCKINESS) {
      if (xvalue < 2.0 || yvalue < 2.0)
        error->one(FLERR,"Invalid blockiness (<2) in set command");
      if(0) {}
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      else if(avec_superquadric) {
        atom->blockiness[i][0] = xvalue;
        atom->blockiness[i][1] = yvalue;
        MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(atom->shape[i], atom->blockiness[i], atom->radius+i); //re-calculate bounding sphere radius
        MathExtraLiggghtsNonspherical::volume_superquadric(atom->shape[i], atom->blockiness[i], atom->volume+i); //re-calculate volume
        MathExtraLiggghtsNonspherical::area_superquadric(atom->shape[i], atom->blockiness[i], atom->area+i); //re-calculate surface area
        atom->rmass[i] = atom->density[i] * atom->volume[i]; //re-calculate mass
        MathExtraLiggghtsNonspherical::inertia_superquadric(atom->shape[i], atom->blockiness[i], atom->density[i], atom->inertia[i]); //re-calculate inertia tensor

      }
      #endif // SUPERQUADRIC_ACTIVE_FLAG
      else
        error->one(FLERR,"Cannot set shape for this type of atom");
    }
    // set desired per-atom property
    else if (keyword == PROPERTYPERATOM) { 

        // if fix was just created, its default values have not been set,
        // so have to add a run 0 to call setup
        if(updFix->just_created)
            error->one(FLERR,"May not use the set command right after fix property/atom without a prior run. Add a 'run 0' between fix property/atom and set");

        if(add == 1)
        {
              currentTimestep = update->ntimestep;
              if (currentTimestep >= until)
                continue;
        }
        if(varflag)
        {
            
            if (updFix->data_style)
            {
                updFix->array_atom[i][0] = xvalue;
                updFix->array_atom[i][1] = yvalue;
                updFix->array_atom[i][2] = zvalue;
            }
            else updFix->vector_atom[i]=dvalue;
        }
        else
        {
            if (updFix->data_style)
            {
                  for (int m = 0; m < nUpdValues; m++)
                     updFix->array_atom[i][m] = updValues[m];
            }
            else updFix->vector_atom[i]=updValues[0];
        }
    }

    // set length of line particle

    else if (keyword == LENGTH) {
      if (dvalue < 0.0) error->one(FLERR,"Invalid length in set command");
      avec_line->set_length(i,dvalue);
    }

    // set rmass via density
    // if radius > 0.0, treat as sphere
    // if shape > 0.0, treat as ellipsoid
    // if length > 0.0, treat as line
    // if area > 0.0, treat as tri
    // else set rmass to density directly

    else if (keyword == DENSITY) {
      if (dvalue <= 0.0) error->one(FLERR,"Invalid density in set command");
      if (atom->radius_flag && atom->radius[i] > 0.0)
      {
          atom->density[i] = dvalue;
          if(!atom->superquadric_flag) {
              if (domain->dimension == 2)
                atom->rmass[i] = MY_PI * atom->radius[i]*atom->radius[i] * atom->density[i]; 
              else
                atom->rmass[i] = 4.0*MY_PI/3.0 * atom->radius[i]*atom->radius[i]*atom->radius[i] * atom->density[i]; 
          }
          else {
            #ifdef SUPERQUADRIC_ACTIVE_FLAG
            if (domain->dimension == 3) {
              atom->rmass[i] = atom->density[i] * atom->volume[i];
              MathExtraLiggghtsNonspherical::inertia_superquadric(atom->shape[i], atom->blockiness[i], atom->density[i], atom->inertia[i]);
            }
            else
              error->one(FLERR,"Superquadrics are implemented only in 3D");
            #endif // SUPERQUADRIC_ACTIVE_FLAG
          }
      }
      else if (atom->density_flag)
        atom->density[i] = dvalue;
      else if (atom->ellipsoid_flag && atom->ellipsoid[i] >= 0) {
        double *shape = avec_ellipsoid->bonus[atom->ellipsoid[i]].shape;
        atom->rmass[i] = 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2] * dvalue;
      } else if (atom->line_flag && atom->line[i] >= 0) {
        double length = avec_line->bonus[atom->line[i]].length;
        atom->rmass[i] = length * dvalue;
      }
       else atom->rmass[i] = dvalue;
    }

    // set dipole moment

    else if (keyword == DIPOLE) {
      double **mu = atom->mu;
      mu[i][0] = xvalue;
      mu[i][1] = yvalue;
      mu[i][2] = zvalue;
      mu[i][3] = sqrt(mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] +
                      mu[i][2]*mu[i][2]);
    }

    // set quaternion orientation of ellipsoid or tri particle or superquadric

    else if (keyword == QUAT || keyword == QUAT_DIRECT ) {
      double *quat = NULL;
      if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
        quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      else if (avec_superquadric)
        quat = avec_superquadric->return_quat_ptr(i);
      #endif
      else
        error->one(FLERR,"Cannot set quaternion for atom that has none");

      if(keyword == QUAT) {
      double theta2 = MY_PI2 * wvalue/180.0;
      double sintheta2 = sin(theta2);
      quat[0] = cos(theta2);
      quat[1] = xvalue * sintheta2;
      quat[2] = yvalue * sintheta2;
      quat[3] = zvalue * sintheta2;
      } else {
        //direct quaternion components setup
        quat[0] = xvalue;
        quat[1] = yvalue;
        quat[2] = zvalue;
        quat[3] = wvalue;
      }
      MathExtra::qnormalize(quat);
    }

    // set theta of line particle

    else if (keyword == THETA) {
      if (atom->line[i] < 0)
        error->one(FLERR,"Cannot set theta for atom that is not a line");
      avec_line->bonus[atom->line[i]].theta = dvalue;
    }

    // set angmom of ellipsoidal or tri particle

    else if (keyword == ANGMOM) {
      atom->angmom[i][0] = xvalue;
      atom->angmom[i][1] = yvalue;
      atom->angmom[i][2] = zvalue;
    }

    // reset any or all of 3 image flags

    else if (keyword == IMAGE) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if (ximageflag) xbox = ximage;
      if (yimageflag) ybox = yimage;
      if (zimageflag) zbox = zimage;
      atom->image[i] = ((tagint) (xbox + IMGMAX) & IMGMASK) |
        (((tagint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((tagint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }

    // set value for custom integer or double vector

    else if (keyword == INAME) {
      atom->ivector[index_custom][i] = ivalue;
    }

    else if (keyword == DNAME) {
      atom->dvector[index_custom][i] = dvalue;
    }

    count++;
  }

  // clear up per-atom memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
  memory->destroy(vec4);
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom tag
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set::setrandom(int keyword)
{
  int i;

  AtomVecEllipsoid *avec_ellipsoid =
    (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  atom->style_match("line"); // DEAD CODE?
  #ifdef SUPERQUADRIC_ACTIVE_FLAG
  AtomVecSuperquadric *avec_superquadric = (AtomVecSuperquadric *) atom->style_match("superquadric");
  #endif

  RanPark *random = new RanPark(lmp, "12345787");
  double **x = atom->x;
  int seed = ivalue;

  // set fraction of atom types to newtype

  if (keyword == TYPE_FRACTION) {
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (select[i]) {
        random->reset(seed,x[i]);
        if (random->uniform() > fraction) continue;
        atom->type[i] = newtype;
        count++;
      }

  // set dipole moments to random orientations in 3d or 2d
  // dipole length is determined by dipole type array

  } else if (keyword == DIPOLE_RANDOM) {
    double **mu = atom->mu;
    int nlocal = atom->nlocal;

    double msq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          mu[i][0] = random->uniform() - 0.5;
          mu[i][1] = random->uniform() - 0.5;
          mu[i][2] = random->uniform() - 0.5;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][2] *= scale;
          mu[i][3] = dvalue;
          count++;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          mu[i][0] = random->uniform() - 0.5;
          mu[i][1] = random->uniform() - 0.5;
          mu[i][2] = 0.0;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][3] = dvalue;
          count++;
        }
    }

  // set quaternions to random orientations in 3d or 2d

  } else if (keyword == QUAT_RANDOM) {
    
    int nlocal = atom->nlocal;
    double *quat = NULL;

    if (domain->dimension == 3) {
      double s,t1,t2,theta1,theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          #ifdef SUPERQUADRIC_ACTIVE_FLAG
          else if (avec_superquadric)
            quat = avec_superquadric->return_quat_ptr(i);
          #endif
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          random->reset(seed,x[i]);
          s = random->uniform();
          t1 = sqrt(1.0-s);
          t2 = sqrt(s);
          theta1 = 2.0*MY_PI*random->uniform();
          theta2 = 2.0*MY_PI*random->uniform();
          quat[0] = cos(theta2)*t2;
          quat[1] = sin(theta1)*t1;
          quat[2] = cos(theta1)*t1;
          quat[3] = sin(theta2)*t2;
          count++;
        }

    } else {
      double theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          random->reset(seed,x[i]);
          theta2 = MY_PI*random->uniform();
          quat[0] = cos(theta2);
          quat[1] = 0.0;
          quat[2] = 0.0;
          quat[3] = sin(theta2);
          count++;
        }
    }
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

void Set::topology(int keyword)
{
  int m,atom1,atom2,atom3,atom4;

  // border swap to acquire ghost atom info
  // enforce PBC before in case atoms are outside box
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen) fprintf(screen,"  system init for set ...\n");
  lmp->init();

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // select both owned and ghost atoms

  selection(atom->nlocal + atom->nghost);

  // for BOND, each of 2 atoms must be in group

  if (keyword == BOND) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++) {
        atom1 = atom->map(atom->bond_atom[i][m]);
        if (atom1 == -1) error->one(FLERR,"Bond atom missing in set command");
        if (select[i] && select[atom1]) {
          atom->bond_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for ANGLE, each of 3 atoms must be in group

  if (keyword == ANGLE) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++) {
        atom1 = atom->map(atom->angle_atom1[i][m]);
        atom2 = atom->map(atom->angle_atom2[i][m]);
        atom3 = atom->map(atom->angle_atom3[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Angle atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3]) {
          atom->angle_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for DIHEDRAL, each of 4 atoms must be in group

  if (keyword == DIHEDRAL) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++) {
        atom1 = atom->map(atom->dihedral_atom1[i][m]);
        atom2 = atom->map(atom->dihedral_atom2[i][m]);
        atom3 = atom->map(atom->dihedral_atom3[i][m]);
        atom4 = atom->map(atom->dihedral_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Dihedral atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->dihedral_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for IMPROPER, each of 4 atoms must be in group

  if (keyword == IMPROPER) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++) {
        atom1 = atom->map(atom->improper_atom1[i][m]);
        atom2 = atom->map(atom->improper_atom2[i][m]);
        atom3 = atom->map(atom->improper_atom3[i][m]);
        atom4 = atom->map(atom->improper_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Improper atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->improper_type[i][m] = ivalue;
          count++;
        }
      }
  }
}

/* ---------------------------------------------------------------------- */

void Set::varparse(char *name, int m)
{
  varflag = 1;

  name = &name[2];
  int n = strlen(name) + 1;
  char *str = new char[n];
  strcpy(str,name);

  int ivar = input->variable->find(str);
  delete [] str;

  if (ivar < 0)
    error->all(FLERR,"Variable name for set command does not exist");
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable for set command is invalid style");

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  } else if (m == 4) {
    varflag4 = 1; ivar4 = ivar;
  }
}
