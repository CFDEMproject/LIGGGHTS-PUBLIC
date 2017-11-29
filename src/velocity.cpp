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

#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "velocity.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "random_park.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "fix_multisphere.h"

using namespace LAMMPS_NS;

enum{SET,SETANGULAR,RAMP,ZERO};
enum{ALL,LOCAL,GEOM};
enum{NONE,CONSTANT,EQUAL,ATOM};

#define WARMUP 100
#define SMALL  0.001

/* ---------------------------------------------------------------------- */

Velocity::Velocity(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Velocity::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal velocity command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Velocity command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Velocity command with no atoms existing");

  // atom masses must all be set

  atom->check_mass();

  // identify group

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find velocity group ID");
  groupbit = group->bitmask[igroup];

  // identify style

  if (strcmp(arg[1],"set") == 0) style = SET;
  else if (strcmp(arg[1],"setAngular")==0) style = SETANGULAR;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"zero") == 0) style = ZERO;
  else error->all(FLERR,"Illegal velocity command");

  // set defaults

  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 1;
  rotation_flag = 0;
  loop_flag = ALL;
  scale_flag = 0; 
  rfix = -1;

  // multispheres
  fix_ms_ =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));

  // read options from end of input line
  // change defaults as options specify

  if (style == SET) options(narg-5,&arg[5]);
  else if (style == SETANGULAR) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == ZERO) options(narg-3,&arg[3]);

  // initialize velocities based on style
  // create() invoked differently, so can be called externally

  if (style == SET) set(narg-2,&arg[2]);
  else if (style == SETANGULAR) setAngular(narg-2,&arg[2]);
  else if (style == RAMP) ramp(narg-2,&arg[2]);
  else if (style == ZERO) zero(narg-2,&arg[2]);
}

/* ----------------------------------------------------------------------
   initialization of defaults before calling velocity methods externaly
------------------------------------------------------------------------- */

void Velocity::init_external(const char *extgroup)
{
  igroup = group->find(extgroup);
  if (igroup == -1) error->all(FLERR,"Could not find velocity group ID");
  groupbit = group->bitmask[igroup];

  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 1;
  rotation_flag = 0;
  loop_flag = ALL;
  scale_flag = 0; 
}

/* ---------------------------------------------------------------------- */

void Velocity::set(int narg, char **arg)
{
  int xstyle,ystyle,zstyle,varflag;
  double vx=0.0,vy=0.0,vz=0.0;
  char *xstr,*ystr,*zstr;
  int xvar=0,yvar=0,zvar=0;

  // parse 3 args

  xstyle = ystyle = zstyle = CONSTANT;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[0],"v_") == arg[0]) {
    int n = strlen(&arg[0][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[0][2]);
  } else if (strcmp(arg[0],"NULL") == 0) xstyle = NONE;
  else vx = force->numeric(FLERR,arg[0]);

  if (strstr(arg[1],"v_") == arg[1]) {
    int n = strlen(&arg[1][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[1][2]);
  } else if (strcmp(arg[1],"NULL") == 0) ystyle = NONE;
  else vy = force->numeric(FLERR,arg[1]);

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[2][2]);
  } else if (strcmp(arg[2],"NULL") == 0) zstyle = NONE;
  else vz = force->numeric(FLERR,arg[2]);

  // set and apply scale factors

  xscale = yscale = zscale = 1.0;

  if (xstyle && !xstr) {
    if (scale_flag) xscale = domain->lattice->xlattice;
    vx *= xscale;
  }
  if (ystyle && !ystr) {
    if (scale_flag) yscale = domain->lattice->ylattice;
    vy *= yscale;
  }
  if (zstyle && !zstr) {
    if (scale_flag) zscale = domain->lattice->zlattice;
    vz *= zscale;
  }

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  // error check for 2d models

  if (domain->dimension == 2) {
    if (zstyle == CONSTANT && vz != 0.0)
      error->all(FLERR,"Cannot set non-zero z velocity for 2d simulation");
    if (zstyle == EQUAL || zstyle == ATOM)
      error->all(FLERR,"Cannot set variable z velocity for 2d simulation");
  }

  // allocate vfield array if necessary

  double **vfield = NULL;
  if (varflag == ATOM) memory->create(vfield,atom->nlocal,3,"velocity:vfield");

  // set velocities via constants

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle) v[i][0] = vx;
          if (ystyle) v[i][1] = vy;
          if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle) v[i][0] += vx;
          if (ystyle) v[i][1] += vy;
          if (zstyle) v[i][2] += vz;
        }
      }
    }

  // set velocities via variables

  } else {
    if (xstyle == EQUAL) vx = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && vfield)
      input->variable->compute_atom(xvar,igroup,&vfield[0][0],3,0);
    if (ystyle == EQUAL) vy = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM && vfield)
      input->variable->compute_atom(yvar,igroup,&vfield[0][1],3,0);
    if (zstyle == EQUAL) vz = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM && vfield)
      input->variable->compute_atom(zvar,igroup,&vfield[0][2],3,0);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle == ATOM) v[i][0] = vfield[i][0];
          else if (xstyle) v[i][0] = vx;
          if (ystyle == ATOM) v[i][1] = vfield[i][1];
          else if (ystyle) v[i][1] = vy;
          if (zstyle == ATOM) v[i][2] = vfield[i][2];
          else if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle == ATOM) v[i][0] += vfield[i][0];
          else if (xstyle) v[i][0] += vx;
          if (ystyle == ATOM) v[i][1] += vfield[i][1];
          else if (ystyle) v[i][1] += vy;
          if (zstyle == ATOM) v[i][2] += vfield[i][2];
          else if (zstyle) v[i][2] += vz;
        }
      }
  }

  // clean up

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy(vfield);
}

/* ----------------------------------------------------------------------
   apply the angular momentum
------------------------------------------------------------------------- */

void Velocity::setAngular(int narg, char **arg)
{
  int xstyle,ystyle,zstyle,varflag;
  double vx=0.0,vy=0.0,vz=0.0;
  char *xstr,*ystr,*zstr;
  int xvar=0,yvar=0,zvar=0;

  // parse 3 args

  xstyle = ystyle = zstyle = CONSTANT;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[0],"v_") == arg[0]) {
    int n = strlen(&arg[0][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[0][2]);
  } else if (strcmp(arg[0],"NULL") == 0) xstyle = NONE;
  else vx = force->numeric(FLERR,arg[0]);

  if (strstr(arg[1],"v_") == arg[1]) {
    int n = strlen(&arg[1][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[1][2]);
  } else if (strcmp(arg[1],"NULL") == 0) ystyle = NONE;
  else vy = force->numeric(FLERR,arg[1]);

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[2][2]);
  } else if (strcmp(arg[2],"NULL") == 0) zstyle = NONE;
  else vz = force->numeric(FLERR,arg[2]);

  // set and apply scale factors

  xscale = yscale = zscale = 1.0;

  if (xstyle && !xstr) {
    if (scale_flag) xscale = domain->lattice->xlattice;
    vx *= xscale;
  }
  if (ystyle && !ystr) {
    if (scale_flag) yscale = domain->lattice->ylattice;
    vy *= yscale;
  }
  if (zstyle && !zstr) {
    if (scale_flag) zscale = domain->lattice->zlattice;
    vz *= zscale;
  }

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  // error check for 2d models

  if (domain->dimension == 2) {
    if (zstyle == CONSTANT && vz != 0.0)
      error->all(FLERR,"Cannot set non-zero z velocity for 2d simulation");
    if (zstyle == EQUAL || zstyle == ATOM)
      error->all(FLERR,"Cannot set variable z velocity for 2d simulation");
  }

  // other error checks
  if (!atom->angmom_flag)
      error->all(FLERR,"setAngular requires the atom property 'angmom' as defined by atom_style ellipsoid.");

  // allocate vfield array if necessary

  double **vfield = NULL;
  if (varflag == ATOM) memory->create(vfield,atom->nlocal,3,"velocity:vfield");

  // set velocities via constants

  double **angmom  = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle) angmom[i][0] = vx;
          if (ystyle) angmom[i][1] = vy;
          if (zstyle) angmom[i][2] = vz;
        } else {
          if (xstyle) angmom[i][0] += vx;
          if (ystyle) angmom[i][1] += vy;
          if (zstyle) angmom[i][2] += vz;
        }
      }
    }

  // set velocities via variables

  } else {
    if (xstyle == EQUAL) vx = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && vfield)
      input->variable->compute_atom(xvar,igroup,&vfield[0][0],3,0);
    if (ystyle == EQUAL) vy = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM && vfield)
      input->variable->compute_atom(yvar,igroup,&vfield[0][1],3,0);
    if (zstyle == EQUAL) vz = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM && vfield)
      input->variable->compute_atom(zvar,igroup,&vfield[0][2],3,0);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle == ATOM) angmom[i][0] = vfield[i][0];
          else if (xstyle) angmom[i][0] = vx;
          if (ystyle == ATOM) angmom[i][1] = vfield[i][1];
          else if (ystyle) angmom[i][1] = vy;
          if (zstyle == ATOM) angmom[i][2] = vfield[i][2];
          else if (zstyle) angmom[i][2] = vz;
        } else {
          if (xstyle == ATOM) angmom[i][0] += vfield[i][0];
          else if (xstyle) angmom[i][0] += vx;
          if (ystyle == ATOM) angmom[i][1] += vfield[i][1];
          else if (ystyle) angmom[i][1] += vy;
          if (zstyle == ATOM) angmom[i][2] += vfield[i][2];
          else if (zstyle) angmom[i][2] += vz;
        }
      }
  }

  // clean up

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy(vfield);
}

/* ----------------------------------------------------------------------
   apply a ramped set of velocities
------------------------------------------------------------------------- */

void Velocity::ramp(int narg, char **arg)
{
  // set scale factors

  if (scale_flag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // parse args

  int v_dim = 0;
  if (strcmp(arg[0],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[0],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[0],"vz") == 0) v_dim = 2;
  else error->all(FLERR,"Illegal velocity command");

  if (v_dim == 2 && domain->dimension == 2)
    error->all(FLERR,"Velocity ramp in z for a 2d problem");

  double v_lo,v_hi;
  if (v_dim == 0) {
    v_lo = xscale*force->numeric(FLERR,arg[1]);
    v_hi = xscale*force->numeric(FLERR,arg[2]);
  } else if (v_dim == 1) {
    v_lo = yscale*force->numeric(FLERR,arg[1]);
    v_hi = yscale*force->numeric(FLERR,arg[2]);
  } else if (v_dim == 2) {
    v_lo = zscale*force->numeric(FLERR,arg[1]);
    v_hi = zscale*force->numeric(FLERR,arg[2]);
  }

  int coord_dim = 0;
  if (strcmp(arg[3],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[3],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[3],"z") == 0) coord_dim = 2;
  else error->all(FLERR,"Illegal velocity command");

  double coord_lo,coord_hi;
  if (coord_dim == 0) {
    coord_lo = xscale*force->numeric(FLERR,arg[4]);
    coord_hi = xscale*force->numeric(FLERR,arg[5]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*force->numeric(FLERR,arg[4]);
    coord_hi = yscale*force->numeric(FLERR,arg[5]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*force->numeric(FLERR,arg[4]);
    coord_hi = zscale*force->numeric(FLERR,arg[5]);
  }

  // vramp = ramped velocity component for v_dim
  // add or set based on sum_flag

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double fraction,vramp;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      if (sum_flag) v[i][v_dim] += vramp;
      else v[i][v_dim] = vramp;
    }
}

/* ----------------------------------------------------------------------
   zero linear or angular momentum of a group
   if using rigid/small requires init of entire system since
      its methods perform forward/reverse comm,
      comm::init needs neighbor::init needs pair::init needs kspace::init, etc
      also requires setup_pre_neighbor call to setup bodies
------------------------------------------------------------------------- */

void Velocity::zero(int narg, char **arg)
{
  if (strcmp(arg[0],"linear") == 0) {
    if (rfix < 0) zero_momentum();
    else {
      if (strcmp(modify->fix[rfix]->style,"rigid/small") == 0) {
        lmp->init();
        modify->fix[rfix]->setup_pre_neighbor();
        modify->fix[rfix]->zero_momentum();
      } else if (strstr(modify->fix[rfix]->style,"rigid")) {
        modify->fix[rfix]->zero_momentum();
      } else error->all(FLERR,"Velocity rigid used with non-rigid fix-ID");
    }

  } else if (strcmp(arg[0],"angular") == 0) {
    if (rfix < 0) zero_rotation();
    else {
      if (strcmp(modify->fix[rfix]->style,"rigid/small") == 0) {
        lmp->init();
        modify->fix[rfix]->setup_pre_neighbor();
        modify->fix[rfix]->zero_rotation();
      } else if (strstr(modify->fix[rfix]->style,"rigid")) {
        modify->fix[rfix]->zero_rotation();
      } else error->all(FLERR,"Velocity rigid used with non-rigid fix-ID");
    }

  } else if (strcmp(arg[0],"angularIndividual") == 0) {
    if (rfix < 0) 
        zero_rotation_individual();
    else
        error->all(FLERR,"angularIndividual not used correctly");
  } else error->all(FLERR,"Illegal velocity command");
}

/* ----------------------------------------------------------------------
   zero the linear momentum of a group of atoms by adjusting v by -Vcm
------------------------------------------------------------------------- */

void Velocity::zero_momentum()
{
  // cannot have no atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Cannot zero momentum of no atoms");

  // compute velocity of center-of-mass of group

  double masstotal = group->mass(igroup);
  double vcm[3];
  group->vcm(igroup,masstotal,vcm);

  // adjust velocities by vcm to zero linear momentum

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double zerovec[3] = {0., 0., 0.};

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vcm[0];
      v[i][1] -= vcm[1];
      v[i][2] -= vcm[2];

      if (fix_ms_)
        fix_ms_->set_v_body_from_atom_index(i,zerovec);

    }
}

/* ----------------------------------------------------------------------
   zero the angular velocity of individual atoms
------------------------------------------------------------------------- */
void Velocity::zero_rotation_individual()
{

  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      omega[i][0] = 0;
      omega[i][1] = 0;
      omega[i][2] = 0;
    }
}

/* ----------------------------------------------------------------------
   zero the angular momentum of a group of atoms by adjusting v by -(w x r)
------------------------------------------------------------------------- */

void Velocity::zero_rotation()
{
  int i;

  // cannot have no atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Cannot zero momentum of no atoms");

  // compute omega (angular velocity) of group around center-of-mass

  double xcm[3],angmom[3],inertia[3][3],omega[3];
  double masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  group->angmom(igroup,xcm,angmom);
  group->inertia(igroup,xcm,inertia);
  group->omega(angmom,inertia,omega);

  // adjust velocities to zero omega
  // vnew_i = v_i - w x r_i
  // must use unwrapped coords to compute r_i correctly

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];

  double zerovec[3] = {0., 0., 0.};

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      v[i][0] -= omega[1]*dz - omega[2]*dy;
      v[i][1] -= omega[2]*dx - omega[0]*dz;
      v[i][2] -= omega[0]*dy - omega[1]*dx;

      if (fix_ms_)
        fix_ms_->set_omega_body_from_atom_index(i,zerovec);

    }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of velocity input line
------------------------------------------------------------------------- */

void Velocity::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal velocity command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dist") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"uniform") == 0) dist_flag = 0;
      else if (strcmp(arg[iarg+1],"gaussian") == 0) dist_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) sum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) sum_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) momentum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) momentum_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rot") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) rotation_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) rotation_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"all") == 0) loop_flag = ALL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop_flag = LOCAL;
      else if (strcmp(arg[iarg+1],"geom") == 0) loop_flag = GEOM;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rigid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      rfix = modify->find_fix(arg[iarg+1]);
      if (rfix < 0) error->all(FLERR,"Fix ID for velocity does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"box") == 0) scale_flag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scale_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else error->all(FLERR,"Illegal velocity command");
  }
}
