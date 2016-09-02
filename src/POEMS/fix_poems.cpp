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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   FixPOEMS is a LAMMPS interface to the POEMS coupled multi-segment simulator
   POEMS authors: Rudranarayan Mukherjee (mukher@rpi.edu)
                  Kurt Anderson (anderk5@rpi.edu)


    POEMS and the POEMS fix has been re-worked by Stefan Radl and 
    Mingqiu WU (TU Graz) to be integrated with LIGGGHTS-TUG's fibre modules
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "workspace.h"
#include "fix_poems.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "force.h"
#include "output.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXSEGMENT 2         // currently 2 since only linear chains allowed
#define DELTA 128
#define TOLERANCE 1.0e-1
#define EPSILON 1.0e-7
#define MAXJACOBI 50

/* ----------------------------------------------------------------------
   define rigid segments and joints, initiate POEMS
------------------------------------------------------------------------- */

FixPOEMS::FixPOEMS(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  mydebug=0;
  int i,j,isegment;

  time_integrate = 1;
  rigid_flag = 1;
  virial_flag = 1;
  model_switch_flag  = 0;

  MPI_Comm_rank(world,&me);

  // perform initial allocation of atom-based arrays
  // register with atom class

  natom2segment = NULL;
  atom2segment = NULL;
  displace = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // initialize each atom to belong to no rigid segments

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) natom2segment[i] = 0;

  // create an atom map if one doesn't exist already
  // readfile() and jointbuild() use GLOBAL atom IDs!!

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_style = 1;
    atom->map_init();
    atom->map_set();
  }

  // parse command-line args
  // set natom2segment, atom2segment for all atoms and nsegment = # of rigid segments
  // atoms must also be in fix group to be in a segment

  if (narg < 4) error->all(FLERR,"Illegal fix poems command");

  // group = arg has list of groups

  if (strcmp(arg[3],"group") == 0) {
    nsegment = narg-6;
    if (nsegment <= 0) error->all(FLERR,"Illegal fix poems command");
    if (strcmp(arg[narg-2],"modelType") == 0)
    if (strcmp(arg[narg-1],"SFJ") == 0){
       model_switch_flag = 1;
    }else 
    {
       model_switch_flag = 0;
    }
    
    int *igroups = new int[nsegment];
    for (isegment = 0; isegment < nsegment; isegment++) {
      igroups[isegment] = group->find(arg[isegment+4]);
      if (igroups[isegment] == -1) 
	    error->all(FLERR,"Could not find fix poems group ID");
    }

    int *mask = atom->mask;

    for (int i = 0; i < nlocal; i++) 
   {
      if (mask[i] & groupbit)
	for (isegment = 0; isegment < nsegment; isegment++)
	  if (mask[i] & group->bitmask[igroups[isegment]]) {
	    if (natom2segment[i] < MAXSEGMENT) atom2segment[i][natom2segment[i]] = isegment;
	    natom2segment[i]++;
	  }
   }

   delete [] igroups;
    
  // file = read segments from file
  // file read doesn't pay attention to fix group,
  //   so after read, reset natom2segment = 0 if atom is not in fix group

  } else if (strcmp(arg[3],"file") == 0) {

    readfile(arg[4]);
    if (strcmp(arg[narg-2],"modelType") == 0)
    if (strcmp(arg[narg-1],"SFJ") == 0){
       model_switch_flag = 1;
    }else 
    {
       model_switch_flag = 0;
    }
    int *mask = atom->mask;
    for (int i = 0; i < nlocal; i++)
       if (!(mask[i] & groupbit)) natom2segment[i] = 1;
  // each molecule in fix group is a rigid segment
  // maxmol = largest molecule #
  // ncount = # of atoms in each molecule (have to sum across procs)
  // nsegment = # of non-zero ncount values
  // use nall as incremented ptr to set atom2segment[] values for each atom

  } else if (strcmp(arg[3],"molecule") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal fix poems command");
    if (atom->molecular == 0)
      error->all(FLERR,"Must use a molecular atom style with fix poems molecule");

    int *mask = atom->mask;
    int *molecule = atom->molecule;
    int nlocal = atom->nlocal;

    int maxmol = -1;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) maxmol = MAX(maxmol,molecule[i]);

    int itmp;
    MPI_Allreduce(&maxmol,&itmp,1,MPI_INT,MPI_MAX,world);
    maxmol = itmp + 1;

    int *ncount = new int[maxmol];
    for (i = 0; i < maxmol; i++) ncount[i] = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) ncount[molecule[i]]++;

    int *nall = new int[maxmol];
    MPI_Allreduce(ncount,nall,maxmol,MPI_INT,MPI_SUM,world);

    nsegment = 0;
    for (i = 0; i < maxmol; i++)
      if (nall[i]) nall[i] = nsegment++;
      else nall[i] = -1;

    for (i = 0; i < nlocal; i++) {
      natom2segment[i] = 0;
      if (mask[i] & groupbit) {
	    natom2segment[i] = 1;
	    atom2segment[i][0] = nall[molecule[i]];
      }
    }
  
    delete [] ncount;
    delete [] nall;

  } else error->all(FLERR,"Illegal fix poems command");

  // error if no segments
  // error if any atom in too many segments

  if (nsegment == 0) error->all(FLERR,"No rigid segments defined");

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (natom2segment[i] > MAXSEGMENT) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"Atom in too many rigid segments - boost MAXSEGMENT");

  //create FixPropertyAtom
  fix_xcm = NULL;
  fix_orientation = NULL;
  
  // create all nsegment-length arrays
  nrigid = new int[nsegment];
  masstotal = new double[nsegment];
  memory->create(xcm,nsegment,3,"poems:xcm");
  memory->create(vcm,nsegment,3,"poems:vcm");
  memory->create(fcm,nsegment,3,"poems:fcm");
  memory->create(inertia,nsegment,3,"poems:inertia");
  memory->create(ex_space,nsegment,3,"poems:ex_space");
  memory->create(ey_space,nsegment,3,"poems:ey_space");
  memory->create(ez_space,nsegment,3,"poems:ez_space");
  memory->create(angmom,nsegment,3,"poems:angmom");
  memory->create(omega,nsegment,3,"poems:omega");
  memory->create(torque,nsegment,3,"poems:torque");

  memory->create(sum,nsegment,6,"poems:sum");
  memory->create(all,nsegment,6,"poems:all");
  
  if (model_switch_flag == 1){
     fprintf(screen, "YOUR ARE USING MODEL TYPE (SphericalFlexiblejoint): %x \n",
     model_switch_flag);
  }
  else{
     fprintf(screen, "YOUR ARE USING MODEL TYPE (SphericalJoint): %x \n",
     model_switch_flag);
  }
  //fprintf(screen, "modelType: %d \n",
  //            model_switch_flag);

  
  // nrigid[n] = # of atoms in Nth rigid segment
  // double count joint atoms as being in multiple segments
  // error if one or zero atoms
 // if (mydebug) fprintf(screen, "nsegment: %d, nclocal: %d, nrigid: %d\n",nsegment,nlocal,nrigid);     
  int *ncount = new int[nsegment];
  for (isegment = 0; isegment < nsegment; isegment++) ncount[isegment] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < natom2segment[i]; j++)
    {
       fprintf(screen, "me: %d, nlocal: %d. particle with i: %d is part of joint j: %d\n",
                me, nlocal, 
                i, j);        
       ncount[j]++;
    }

  }

  MPI_Allreduce(ncount,nrigid,nsegment,MPI_INT,MPI_SUM,world);
  delete [] ncount;


  for (isegment = 0; isegment < nsegment; isegment++)
//FIXME: Allow 1 atom in each segment - Previously atoms needed to be part of multiple segments to define a joint
//    if (nrigid[isegment] <= 1) error->all(FLERR,"One or zero atoms in rigid segment"); 
    //if (nrigid[isegment] <= 0) error->all(FLERR,"Zero atoms in rigid segment");
      fprintf(screen, "nsegment: %d, nclocal: %d, nrigid: %d\n",nsegment,nlocal,nrigid[isegment]);   
  // build list of joint connections and check for cycles and trees
  jointbuild();
  
  // delete temporary atom map
//  fprintf(screen, "delete temporary atom map..\n");
//  if (mapflag) 
//  {
//    atom->map_delete();
//    atom->map_style = 0;
//  }

  // create POEMS instance 
  poems = new Workspace;
  
  // print statistics
  int nsum = 0;
  for (isegment = 0; isegment < nsegment; isegment++) nsum += nrigid[isegment];
  
  if (me == 0) {
    if (screen)
      fprintf(screen,"%d clusters, %d segments, %d joints, %d atoms\n",
	      ncluster,nsegment,njoint,nsum);
    if (logfile)
      fprintf(logfile,"%d clusters, %d segments, %d joints, %d atoms\n",
	      ncluster,nsegment,njoint,nsum);
  }
  if(mydebug) fprintf(screen, "POEMS initialized!\n");
}

/* ----------------------------------------------------------------------
   free all memory for rigid segments, joints, and POEMS
------------------------------------------------------------------------- */

FixPOEMS::~FixPOEMS()
{

  // if atom class still exists:
  //   unregister this fix so atom class doesn't invoke it any more

  if (atom) atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(natom2segment);
  memory->destroy(atom2segment);
  memory->destroy(displace);

  // delete nsegment-length arrays

  delete [] nrigid;
  delete [] masstotal;
  memory->destroy(xcm);
  memory->destroy(vcm);
  memory->destroy(fcm);
  memory->destroy(inertia);
  memory->destroy(ex_space);
  memory->destroy(ey_space);
  memory->destroy(ez_space);
  memory->destroy(angmom);
  memory->destroy(omega);
  memory->destroy(torque);

  memory->destroy(sum);
  memory->destroy(all);

  // delete joint arrays

  memory->destroy(jointsegment);
  memory->destroy(xjoint);
  delete [] freelist;

  // delete POEMS object

  delete poems;
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::post_create()
{

  if(mydebug) fprintf(screen, "POEMS::post_create()!\n");
  // register fixes for quantities to be saved to disk
  // see fix_property_atom.cpp for meaning of fixargs 
  if(!fix_xcm)
  {
        char* fixarg[11];
        fixarg[0]="xcm";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="xcm";
        fixarg[4]="vector";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
	    fixarg[9]="0.";
	    fixarg[10]="0.";
	    fix_xcm = modify->add_fix_property_atom(11,fixarg,style);
  }
  if(!fix_orientation)
  {
        char* fixarg[11];
        fixarg[0]="ex";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="ex";
        fixarg[4]="vector";
        fixarg[5]="no";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.";
	    fixarg[9]="0.";
	    fixarg[10]="0.";
	    fix_orientation = modify->add_fix_property_atom(11,fixarg,style);
  }

}

/* ---------------------------------------------------------------------- */

void FixPOEMS::updatePtrs()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int i, isegment;


  //Set FixPropertyAtom for each atom in segment
  for (i = 0; i < nlocal; i++) 
  {
    if (natom2segment[i]) 
    {
      isegment = atom2segment[i][0];
      fix_xcm->array_atom[i][0] = xcm[isegment][0];
      fix_xcm->array_atom[i][1] = xcm[isegment][1];
      fix_xcm->array_atom[i][2] = xcm[isegment][2];

     //Search joint that belongs to this segment
     //Need to refresh first
     if(njoint)
     {
           fix_orientation->array_atom[i][0] = ex_space[isegment][0] ;
           fix_orientation->array_atom[i][1] = ex_space[isegment][1] ;
           fix_orientation->array_atom[i][2] = ex_space[isegment][2] ;
 

//                fprintf(screen, "i: %d ex: %g %g %g \n", 
//                        i,
//                        fix_orientation->array_atom[i][0],
//                        fix_orientation->array_atom[i][1],
//                        fix_orientation->array_atom[i][2]); 
     }
    }
  }
}
/* ---------------------------------------------------------------------- */

int FixPOEMS::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  mask |= POST_FORCE_RESPA;  
  return mask;
}

int FixPOEMS::model_switch_value()
{
  // to communicate the choice of Spherical Flexible Joint Model
  return model_switch_flag;
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::init()
{

  int i,isegment;

  // warn if more than one POEMS fix

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"poems") == 0) count++;
  if (count > 1 && comm->me == 0) error->warning(FLERR,"More than one fix poems");

  // error if npt,nph fix comes before rigid fix

  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) break;
    if (strcmp(modify->fix[i]->style,"nph") == 0) break;
  }
  if (i < modify->nfix) {
    for (int j = i; j < modify->nfix; j++)
      if (strcmp(modify->fix[j]->style,"poems") == 0)
	error->all(FLERR,"POEMS fix must come before NPT/NPH fix");
  }

  // timestep info

  dtv = update->dt;  
  dtf = 0.5 * update->dt * force->ftm2v;  
  dthalf = 0.5 * update->dt;

  // rRESPA info

  if (strstr(update->integrate_style,"respa")) {
    step_respa = ((Respa *) update->integrate)->step;
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
  }

  // compute masstotal & center-of-mass xcm of each rigid segment
  // only count joint atoms in 1st segment

  int *type = atom->type;
  int *image = atom->image;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int xbox,ybox,zbox;
  double massone;

  for (isegment = 0; isegment < nsegment; isegment++)
    for (i = 0; i < 6; i++) sum[isegment][i] = 0.0; //reset

  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i]) {
      isegment = atom2segment[i][0];
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }
      sum[isegment][0] += (x[i][0] + xbox*xprd) * massone;
      sum[isegment][1] += (x[i][1] + ybox*yprd) * massone;
      sum[isegment][2] += (x[i][2] + zbox*zprd) * massone;
      sum[isegment][3] += massone;
      sum[isegment][4] += massone *
	                  (  v[i][0]*v[i][0] 
                       + v[i][1]*v[i][1] 
                       + v[i][2]*v[i][2]
                      );
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nsegment,MPI_DOUBLE,MPI_SUM,world);

  total_ke = 0.0;
  for (isegment = 0; isegment < nsegment; isegment++) {
    masstotal[isegment] = all[isegment][3];
    xcm[isegment][0] = all[isegment][0]/masstotal[isegment];
    xcm[isegment][1] = all[isegment][1]/masstotal[isegment];
    xcm[isegment][2] = all[isegment][2]/masstotal[isegment];
    total_ke += 0.5 * all[isegment][4];
  }

  // compute 6 moments of inertia of each segment
  // only count joint atoms in 1st segment
  // dx,dy,dz = coords relative to center-of-mass

  double dx,dy,dz;

  for (isegment = 0; isegment < nsegment; isegment++)
    for (i = 0; i < 6; i++) sum[isegment][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i]) {
      isegment = atom2segment[i][0];

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[isegment][0];
      dy = x[i][1] + ybox*yprd - xcm[isegment][1];
      dz = x[i][2] + zbox*zprd - xcm[isegment][2];
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }
    //HARDCODE: CYLINDER-MOMENT OF INERTIA
	//Cylinder in x-direction, each segment consists of cylinders

	double radius=0.05;
	double length=1.0;
      sum[isegment][0] += 1.0/2.0 * massone * radius * radius;
      sum[isegment][1] += 1.0/12.0*massone * (3.0 * radius * radius + length * length);
      sum[isegment][2] += 1.0/12.0*massone * (3.0 * radius * radius + length * length);
      sum[isegment][3] -= 0;
      sum[isegment][4] -= 0;
      sum[isegment][5] -= 0;

    }
  }

  MPI_Allreduce(sum[0],all[0],6*nsegment,MPI_DOUBLE,MPI_SUM,world);

  // inertia = 3 eigenvalues = principal moments of inertia
  // ex_space,ey_space,ez_space = 3 eigenvectors = principal axes of rigid segment

  double **tensor,**evectors;
  memory->create(tensor,3,3,"fix_rigid:tensor");
  memory->create(evectors,3,3,"fix_rigid:evectors");

  int ierror;
  double ez0,ez1,ez2;

  for (isegment = 0; isegment < nsegment; isegment++) {
    tensor[0][0] = all[isegment][0];
    tensor[1][1] = all[isegment][1];
    tensor[2][2] = all[isegment][2];
    tensor[0][1] = tensor[1][0] = all[isegment][3];
    tensor[1][2] = tensor[2][1] = all[isegment][4];
    tensor[0][2] = tensor[2][0] = all[isegment][5];
  
    ierror = jacobi(tensor,inertia[isegment],evectors);
    if (ierror) error->all(FLERR,"Insufficient Jacobi rotations for POEMS segment");

    ex_space[isegment][0] = evectors[0][0];
    ex_space[isegment][1] = evectors[1][0];
    ex_space[isegment][2] = evectors[2][0];
    
    ey_space[isegment][0] = evectors[0][1];
    ey_space[isegment][1] = evectors[1][1];
    ey_space[isegment][2] = evectors[2][1];
    
    ez_space[isegment][0] = evectors[0][2];
    ez_space[isegment][1] = evectors[1][2];
    ez_space[isegment][2] = evectors[2][2];
    
    // if any principal moment < scaled EPSILON, error
    // this is b/c POEMS cannot yet handle degenerate segments
  
    double max;
    max = MAX(inertia[isegment][0],inertia[isegment][1]);
    max = MAX(max,inertia[isegment][2]);
    
//    fprintf(screen, "max: %.3g; isegment: %d; inertia[isegment][i] %g %g %g \n", 
//                max,  
//                isegment,
//                inertia[isegment][0], inertia[isegment][1], inertia[isegment][2]);
  
    if (inertia[isegment][0] < EPSILON*max ||
	inertia[isegment][1] < EPSILON*max ||
	inertia[isegment][2] < EPSILON*max)
      error->all(FLERR,"Rigid segment has degenerate moment of inertia");

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd evector if needed
  
    ez0 = ex_space[isegment][1]*ey_space[isegment][2] -
      ex_space[isegment][2]*ey_space[isegment][1];
    ez1 = ex_space[isegment][2]*ey_space[isegment][0] -
      ex_space[isegment][0]*ey_space[isegment][2];
    ez2 = ex_space[isegment][0]*ey_space[isegment][1] -
      ex_space[isegment][1]*ey_space[isegment][0];
  
    if (ez0*ez_space[isegment][0] + ez1*ez_space[isegment][1] + 
	ez2*ez_space[isegment][2] < 0.0) {
      ez_space[isegment][0] = -ez_space[isegment][0];
      ez_space[isegment][1] = -ez_space[isegment][1];
      ez_space[isegment][2] = -ez_space[isegment][2];
    }
  }

  // free temporary memory
  
  memory->destroy(tensor);
  memory->destroy(evectors);

  // displace = initial atom coords in basis of principal axes
  // only set joint atoms relative to 1st segment
  // set displace = 0.0 for atoms not in any rigid segment

  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i]) {
      isegment = atom2segment[i][0];

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[isegment][0];
      dy = x[i][1] + ybox*yprd - xcm[isegment][1];
      dz = x[i][2] + zbox*zprd - xcm[isegment][2];
      
      displace[i][0] = dx*ex_space[isegment][0] + dy*ex_space[isegment][1] +
	dz*ex_space[isegment][2];
      displace[i][1] = dx*ey_space[isegment][0] + dy*ey_space[isegment][1] +
	dz*ey_space[isegment][2];
      displace[i][2] = dx*ez_space[isegment][0] + dy*ez_space[isegment][1] +
	dz*ez_space[isegment][2];
    } else displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
  }  

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // only count joint atoms in 1st segment
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // (ddx,ddy,ddz) is projection of atom within rigid segment onto principal axes
  // 6 moments use  (ddx,ddy,ddz) displacements from principal axes

  for (isegment = 0; isegment < nsegment; isegment++)
    for (i = 0; i < 6; i++) sum[isegment][i] = 0.0;

  double ddx,ddy,ddz;

  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i]) {
      isegment = atom2segment[i][0];

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[isegment][0];
      dy = x[i][1] + ybox*yprd - xcm[isegment][1];
      dz = x[i][2] + zbox*zprd - xcm[isegment][2];
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }

      ddx = dx*ex_space[isegment][0] + dy*ex_space[isegment][1] + 
	dz*ex_space[isegment][2];
      ddy = dx*ey_space[isegment][0] + dy*ey_space[isegment][1] +
	dz*ey_space[isegment][2];
      ddz = dx*ez_space[isegment][0] + dy*ez_space[isegment][1] +
	dz*ez_space[isegment][2];

//HARDCODE
	float radius=0.05;
	float length=1.0;
      sum[isegment][0] += 1./2. * massone * radius * radius;
      sum[isegment][1] += 1./12.*massone * (3. * radius * radius + length * length);
      sum[isegment][2] += 1./12.*massone * (3. * radius * radius + length * length);
      sum[isegment][3] -= 0.;
      sum[isegment][4] -= 0.;
      sum[isegment][5] -= 0.;
/*
      sum[isegment][0] += massone * (ddy*ddy + ddz*ddz);
      sum[isegment][1] += massone * (ddx*ddx + ddz*ddz);
      sum[isegment][2] += massone * (ddx*ddx + ddy*ddy);
      sum[isegment][3] -= massone * ddx*ddy;
      sum[isegment][4] -= massone * ddy*ddz;
      sum[isegment][5] -= massone * ddx*ddz;
*/
    }
  }
  
  MPI_Allreduce(sum[0],all[0],6*nsegment,MPI_DOUBLE,MPI_SUM,world);
  
  for (isegment = 0; isegment < nsegment; isegment++) {
    if (fabs(all[isegment][0]-inertia[isegment][0]) > TOLERANCE || 
	fabs(all[isegment][1]-inertia[isegment][1]) > TOLERANCE ||
	fabs(all[isegment][2]-inertia[isegment][2]) > TOLERANCE)
      error->all(FLERR,"Bad principal moments - type 1");
    if (fabs(all[isegment][3]) > TOLERANCE || 
	fabs(all[isegment][4]) > TOLERANCE ||
	fabs(all[isegment][5]) > TOLERANCE)
      error->all(FLERR,"Bad principal moments - type 2");
  }
  //fprintf(screen,"%g inertia\n",
  //	      inertia[1][1]);
  // find fix and assign values
  fix_xcm = static_cast<FixPropertyAtom*>(modify->find_fix_property("xcm","property/atom","vector",0,0,style));
  fix_orientation= static_cast<FixPropertyAtom*>(modify->find_fix_property("ex","property/atom","vector",0,0,style));

}

/* ----------------------------------------------------------------------
   compute initial rigid segment info
   make setup call to POEMS
------------------------------------------------------------------------- */

void FixPOEMS::setup(int vflag)
{

  int i,n,isegment;

  // vcm = velocity of center-of-mass of each rigid segment
  // angmom = angular momentum of each rigid segment
  // only count joint atoms in 1st segment

  int *type = atom->type;
  int *image = atom->image;
  double *rmass = atom->rmass; 
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int xbox,ybox,zbox;
  double massone,dx,dy,dz;

  for (isegment = 0; isegment < nsegment; isegment++)
    for (i = 0; i < 6; i++) sum[isegment][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i]) {
      isegment = atom2segment[i][0];
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }
       //fprintf(screen, "massone: %g \n",massone);
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[isegment][0];
      dy = x[i][1] + ybox*yprd - xcm[isegment][1];
      dz = x[i][2] + zbox*zprd - xcm[isegment][2];

      sum[isegment][0] += v[i][0] * massone;
      sum[isegment][1] += v[i][1] * massone;
      sum[isegment][2] += v[i][2] * massone;
      sum[isegment][3] += dy * massone*v[i][2] - dz * massone*v[i][1];
      sum[isegment][4] += dz * massone*v[i][0] - dx * massone*v[i][2];
      sum[isegment][5] += dx * massone*v[i][1] - dy * massone*v[i][0]; 
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nsegment,MPI_DOUBLE,MPI_SUM,world);

  for (isegment = 0; isegment < nsegment; isegment++) {
    vcm[isegment][0] = all[isegment][0]/masstotal[isegment];
    vcm[isegment][1] = all[isegment][1]/masstotal[isegment];
    vcm[isegment][2] = all[isegment][2]/masstotal[isegment];
    angmom[isegment][0] = all[isegment][3];
    angmom[isegment][1] = all[isegment][4];
    angmom[isegment][2] = all[isegment][5];  
  }

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set omega from angmom & orientation of rigid segment

  for (isegment = 0; isegment < nsegment; isegment++)
    omega_from_mq(angmom[isegment],ex_space[isegment],ey_space[isegment],
		  ez_space[isegment],inertia[isegment],omega[isegment]);
  set_v();

  // guestimate virial as 2x the set_v contribution

  if (vflag_global)
    for (n = 0; n < 6; n++) virial[n] *= 2.0;
  if (vflag_atom) {
    for (i = 0; i < nlocal; i++)
      for (n = 0; n < 6; n++)
	vatom[i][n] *= 2.0;
  }

  // use post_force() to compute initial fcm & torque

  post_force(vflag);

  // setup for POEMS

  poems->MakeSystem(nsegment,masstotal,inertia,xcm,vcm,omega,
		    ex_space,ey_space,ez_space,
		    njoint,jointsegment,xjoint,nfree,freelist,
		    dthalf,dtv,force->ftm2v,total_ke,model_switch_flag);
   if(mydebug){
   int currAtom=0;
   fprintf(screen, "masstotal: %g, currAtom: %d dthalf %g,dtv %g ,force->ftm2v %g , total_ke %g\n",
		     masstotal[currAtom], currAtom,
             dthalf,
             dtv,
             force->ftm2v,
             total_ke);

   }
  //update fixes to report fibre data
  updatePtrs();

}

/* ----------------------------------------------------------------------
   update vcm,omega by 1/2 step and xcm,orientation by full step
   set x,v of segment atoms accordingly
   ---------------------------------------------------------------------- */

void FixPOEMS::initial_integrate(int vflag)
{

  if(mydebug) fprintf(screen, "POEMS::initial_integrate()!\n");

  // perform POEMS integration

   poems->LobattoOne(xcm,vcm,omega,torque,fcm,ex_space,ey_space,ez_space);
/*
   int currAtom=1;
   fprintf(screen, "currAtom: %d xcm %g %g %g,vcm %g %g %g ,omega %g %g %g, torque  %g %g %g, fcm  %g %g %g\n",
            currAtom,
             xcm[currAtom][0], xcm[currAtom][1], xcm[currAtom][2],
             vcm[currAtom][0], vcm[currAtom][1], vcm[currAtom][2],
             omega[currAtom][0], omega[currAtom][1], omega[currAtom][2],
             torque[currAtom][0],torque[currAtom][1],torque[currAtom][2],
             fcm[currAtom][0],fcm[currAtom][1],fcm[currAtom][2]);

   fprintf(screen, "currAtom: %d e_space_1 %g %g %g,e_space_2 %g %g %g ,e_space_3 %g %g %g \n",
            currAtom,
            ex_space[currAtom][0],ey_space[currAtom][0],ez_space[currAtom][0],
            ex_space[currAtom][1],ey_space[currAtom][1],ez_space[currAtom][1],
            ex_space[currAtom][2],ey_space[currAtom][2],ez_space[currAtom][2]);*/

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords and velocities of atoms in rigid segments
  set_xv();
}

/* ----------------------------------------------------------------------
   compute fcm,torque on each rigid segment
   only count joint atoms in 1st segment
------------------------------------------------------------------------- */

void FixPOEMS::post_force(int vflag)
{
  int i,isegment;
  int xbox,ybox,zbox;
  double dx,dy,dz;

  int *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  for (isegment = 0; isegment < nsegment; isegment++)
    for (i = 0; i < 6; i++) sum[isegment][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) 
  {
      isegment = i;    //exploit fact that each atom is a segment

      sum[isegment][0] += f[i][0];
      sum[isegment][1] += f[i][1];
      sum[isegment][2] += f[i][2];
      
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[isegment][0];
      dy = x[i][1] + ybox*yprd - xcm[isegment][1];
      dz = x[i][2] + zbox*zprd - xcm[isegment][2];
    
      sum[isegment][3] += dy*f[i][2] - dz*f[i][1];
      sum[isegment][4] += dz*f[i][0] - dx*f[i][2];
      sum[isegment][5] += dx*f[i][1] - dy*f[i][0];

//      fprintf(screen, "sum[%d]: %g %g %g %g %g %g \n", 
//                 isegment,
//                 sum[isegment][0],sum[isegment][1],sum[isegment][2],sum[isegment][3],sum[isegment][4],sum[isegment][5]);
  }
  
  MPI_Allreduce(sum[0],all[0],6*nsegment,MPI_DOUBLE,MPI_SUM,world);

  for (isegment = 0; isegment < nsegment; isegment++) {
    fcm[isegment][0] = all[isegment][0];
    fcm[isegment][1] = all[isegment][1];
    fcm[isegment][2] = all[isegment][2];
    torque[isegment][0] = all[isegment][3];
    torque[isegment][1] = all[isegment][4];
    torque[isegment][2] = all[isegment][5];
  }
}

/* ----------------------------------------------------------------------
   update vcm,omega by last 1/2 step
   set v of segment atoms accordingly
------------------------------------------------------------------------- */

void FixPOEMS::final_integrate()
{

  // perform POEMS integration
  poems->LobattoTwo(vcm,omega,torque,fcm);

  // set velocities of atoms in rigid segments
  // virial is already setup from initial_integrate

  set_v();
  
  //update fixes to report fibre data
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ----------------------------------------------------------------------
   remap xcm of each rigid segment back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   if don't do this, then atoms of a segment which drifts far away
     from a triclinic box will be remapped back into box
     with huge displacements when the box tilt changes via set_x() 
   NOTE: cannot do this by changing xcm of each segment in cluster
         or even 1st segment in cluster
	 b/c POEMS library does not see xcm but only sets xcm
	 so remap needs to be coordinated with POEMS library
	 thus this routine does nothing for now
------------------------------------------------------------------------- */

void FixPOEMS::pre_neighbor() {}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_poems for atoms in igroup 
------------------------------------------------------------------------- */

int FixPOEMS::dof(int igroup)
{

  int groupbit = group->bitmask[igroup];

  // ncount = # of atoms in each rigid segment that are also in group
  // only count joint atoms as part of first segment

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *ncount = new int[nsegment];
  for (int isegment = 0; isegment < nsegment; isegment++) ncount[isegment] = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (natom2segment[i]) ncount[atom2segment[i][0]]++;

  int *nall = new int[nsegment];
  MPI_Allreduce(ncount,nall,nsegment,MPI_INT,MPI_SUM,world);

  // remove 3N - 6 dof for each rigid segment if at least 2 atoms are in igroup

  int n = 0;
  for (int isegment = 0; isegment < nsegment; isegment++)
    if (nall[isegment] > 2) n += 3*nall[isegment] - 6;

  // subtract 3 additional dof for each joint if atom is also in igroup

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (natom2segment[i] > 1 && (mask[i] & groupbit)) m += 3*(natom2segment[i]-1);
  int mall;
  MPI_Allreduce(&m,&mall,1,MPI_INT,MPI_SUM,world);
  n += mall;

  // delete local memory

  delete [] ncount;
  delete [] nall;

  return n;
}

/* ----------------------------------------------------------------------
   adjust xcm of each cluster due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
   NOTE: cannot do this by changing xcm of each segment in cluster
         or even 1st segment in cluster
	 b/c POEMS library does not see xcm but only sets xcm
	 so deform needs to be coordinated with POEMS library
	 thus this routine does nothing for now
------------------------------------------------------------------------- */

void FixPOEMS::deform(int flag) {}

/* ---------------------------------------------------------------------- */

void FixPOEMS::readfile(char *file)
{


  FILE *fp;

  if (me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix poems file %s",file);
      error->one(FLERR,str);
    }
  }

  nsegment = 0;
  nsystem = 0; 
  char *line = NULL;
  int maxline = 0;
  char *ptr;
  int nlocal = atom->nlocal;
  int i,id,nlen;

  while (1) {
    if (me == 0) nlen = readline(fp,&line,&maxline);
    MPI_Bcast(&nlen,1,MPI_INT,0,world);
    if (nlen == 0) break;
    MPI_Bcast(line,nlen,MPI_CHAR,0,world);

    ptr = strtok(line," ,\t\n\0");
    if (ptr == NULL || ptr[0] == '#') continue;
    ptr = strtok(NULL," ,\t\n\0");

    while (ptr = strtok(NULL," ,\t\n\0")) { //FIXME?
      id = atoi(ptr);
      i = atom->map(id);
      printf("id: %d, local id: %d \n", id, i); 

      if (i < 0 || i >= nlocal) continue;
      if (natom2segment[i] < MAXSEGMENT){ 
      atom2segment[i][natom2segment[i]] = i;
      nsegment++;
      //printf("atom2segment]: %d \n", atom2segment[i][natom2segment[i]]); 
      }else 
      {
          printf("natom2segment: %d \n", natom2segment[i]);
          error->one(FLERR,"natom2segment[i] >= MAXSEGMENT. Cannot proceed.\n");
          natom2segment[i]++;
      }
      printf("natom2segment: %d \n", natom2segment[i]);
    }
    printf("Filled system %d. \n", nsystem);
    nsystem++;
  }

  memory->destroy(line);
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPOEMS::readline(FILE *fp, char **pline, int *pmaxline)
{


  int n = 0;
  char *line = *pline;
  int maxline = *pmaxline;

  while (1) {
    if (n+1 >= maxline) {
      maxline += DELTA;
      memory->grow(line,maxline,"fix_poems:line");
    }
    if (fgets(&line[n],maxline-n,fp) == NULL) {
      n = 0;
      break;
    }
    n = strlen(line);
    if (n < maxline-1 || line[n-1] == '\n') break;
  }

  *pmaxline = maxline;
  *pline = line;
  return n;
}

/* ----------------------------------------------------------------------
   build list of joints and error check for cycles and trees
------------------------------------------------------------------------- */

void FixPOEMS::jointbuild()
{

  int i,j;

  //WE DONT WANT JOINT ATOMS; RATHER JOINTS!
  // convert atom2segment into list of joint atoms on this proc
  // local_cpu_joint = # of joint atoms in this proc
  // an atom in N rigid segments, infers N-1 joints between 1st segment and others
  // mylist = [0],[1] = 2 segment indices, [2] = global ID of joint atom
  double **x = atom->x;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
 
  int local_cpu_joint = nlocal - 1; //TODO: ISSUE: number of joints on this CPU


  for (i = 0; i < nlocal; i++) {
    if (natom2segment[i] <= 0) continue; //FIXME: changed from 1
    local_cpu_joint += natom2segment[i]-1;
  }


  // jlist = mylist concatenated across all procs via MPI_Allgatherv
  MPI_Allreduce(&local_cpu_joint,&njoint,1,MPI_INT,MPI_SUM,world);
  int **jlist = NULL;
  if (njoint) memory->create(jlist,njoint,3,"poems:jlist");

  int nprocs;
  MPI_Comm_size(world,&nprocs);

  int *recvcounts = new int[nprocs];
  int tmp = 3*local_cpu_joint;
  MPI_Allgather(&tmp,1,MPI_INT,recvcounts,1,MPI_INT,world);
  fprintf(screen, "nprocess=: %i, tmp=:%i, localCPUjoint=:%i,njoint=:%i \n",nprocs,tmp,local_cpu_joint,njoint);
  int *displs = new int[nprocs];
  displs[0] = 0;
  for (i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  delete [] recvcounts;
  delete [] displs;

  // warning if no joints

  if (njoint == 0 && me == 0)
    error->warning(FLERR,"No joints between rigid segments, use fix rigid instead");

  // sort joint list in ascending order by segment indices
  // check for loops in joint connections between rigid segments
  // check for trees = same segment in more than 2 joints
  //sortlist(njoint,jlist);


  // allocate and setup joint arrays
  // jointsegment stores segment indices from 1 to nsegment to pass to POEMS
  // each proc sets xjointTemp if it owns joint atom
  // MPI_Allreduce gives all procs the xjoint coords

  jointsegment = NULL;
  xjoint = NULL;
  double **xjointTemp = NULL;
  if (njoint) {
    memory->create(jointsegment,njoint,2,"poems:jointsegment");
    memory->create(xjoint,njoint,3,"poems:xjoint");
    memory->create(xjointTemp,njoint,3,"poems:xjointTemp");
  }


  for (i = 0; i < njoint; i++) 
  {
    
    //Create joints between atoms stored in nlocal (exactly in the middle)
    //Joint position
    xjointTemp[i][0]=(x[i][0]+x[i+1][0])/2;
    xjointTemp[i][1]=(x[i][1]+x[i+1][1])/2;
    xjointTemp[i][2]=(x[i][2]+x[i+1][2])/2;
    //xjointTemp[i][0]=(x[0][0]+x[1][0])/2;
    //xjointTemp[i][1]=(x[0][1]+x[1][1])/2;
    //xjointTemp[i][2]=(x[0][2]+x[1][2])/2;

    jointsegment[i][0] = i+1; //TODO: BIG ISSUE!!
    jointsegment[i][1] = i+2;  

  }
  
  if (mydebug) fprintf(screen, "me: %d, joint list complete. \n");

  if (njoint)  
    MPI_Allreduce(xjointTemp[0],xjoint[0],3*njoint,MPI_DOUBLE,MPI_SUM,world);

  // compute freelist of nfree single unconnected segments
  // POEMS could do this itself
  //HARDCODED
  nfree = 0;
  freelist = NULL;
  ncluster = 1;

  // free memory local to this routine
  memory->destroy(xjointTemp);

  for (i = 0; i < njoint; i++) 
  {
    fprintf(screen, "me: %d, xjoint[%d]: %g %g %g \n",
              me,
              i, 
              xjoint[i][0],xjoint[i][1],xjoint[i][2]);
  }
  

}

/* ----------------------------------------------------------------------
  sort joint list (Numerical Recipes shell sort)
  sort criterion: sort on 1st segment, if equal sort on 2nd segment
------------------------------------------------------------------------- */

void FixPOEMS::sortlist(int n, int **list)
{

  int i,j,v0,v1,v2,flag;

  int inc = 1;
  while (inc <= n) inc = 3*inc + 1;

  do {
    inc /= 3;
    for (i = inc+1; i <= n; i++) {
      v0 = list[i-1][0];
      v1 = list[i-1][1];
      v2 = list[i-1][2];
      j = i;
      flag = 0;
      if (list[j-inc-1][0] > v0 || 
	  (list[j-inc-1][0] == v0 && list[j-inc-1][1] > v1)) flag = 1;
      while (flag) {
	list[j-1][0] = list[j-inc-1][0];
	list[j-1][1] = list[j-inc-1][1];
	list[j-1][2] = list[j-inc-1][2];
	j -= inc;
	if (j <= inc) break;
	flag = 0;
	if (list[j-inc-1][0] > v0 || 
	    (list[j-inc-1][0] == v0 && list[j-inc-1][1] > v1)) flag = 1;
      }
      list[j-1][0] = v0;
      list[j-1][1] = v1;
      list[j-1][2] = v2;
    }
  } while (inc > 1);
}

/* ----------------------------------------------------------------------
  check for cycles in list of joint connections between rigid segments
  treat as graph: vertex = segment, edge = joint between 2 segments
------------------------------------------------------------------------- */

int FixPOEMS::loopcheck(int nvert, int nedge, int **elist)
{

  int i,j,k;

  // ecount[i] = # of vertices connected to vertex i via edge
  // elistfull[i][*] = list of vertices connected to vertex i

  int *ecount = new int[nvert];
  for (i = 0; i < nvert; i++) ecount[i] = 0;
  for (i = 0; i < nedge; i++) {
    ecount[elist[i][0]]++;
    ecount[elist[i][1]]++;
  }

  int emax = 0;
  for (i = 0; i < nvert; i++) emax = MAX(emax,ecount[i]);
  
  int **elistfull;
  memory->create(elistfull,nvert,emax,"poems:elistfull");
  for (i = 0; i < nvert; i++) ecount[i] = 0;
  for (i = 0; i < nedge; i++) {
    elistfull[elist[i][0]][ecount[elist[i][0]]++] = elist[i][1];
    elistfull[elist[i][1]][ecount[elist[i][1]]++] = elist[i][0];
  }

  // cycle detection algorithm
  // mark = 0/1 marking of each vertex, all initially unmarked
  // outer while loop:
  //   if all vertices are marked, no cycles, exit loop
  //   push an unmarked vertex on stack and mark it, parent is -1
  //   while stack is not empty:
  //     pop vertex I from stack
  //     loop over vertices J connected to I via edge
  //       if J is parent (vertex that pushed I on stack), skip it
  //       else if J is marked, a cycle is found, return 1
  //       else push J on stack and mark it, parent is I
  //   increment ncluster each time stack empties since that is new cluster

  int *parent = new int[nvert];
  int *mark = new int[nvert];
  for (i = 0; i < nvert; i++) mark[i] = 0;

  int nstack = 0;
  int *stack = new int[nvert];
  ncluster = 0;

  while (1) {
    for (i = 0; i < nvert; i++)
      if (mark[i] == 0) break;
    if (i == nvert) break;
    stack[nstack++] = i;
    mark[i] = 1;
    parent[i] = -1;

    while (nstack) {
      i = stack[--nstack];
      for (k = 0; k < ecount[i]; k++) {
	j = elistfull[i][k];
	if (j == parent[i]) continue;
	if (mark[j]) return 1;
	stack[nstack++] = j;
	mark[j] = 1;
	parent[j] = i;
      }
    }
    ncluster++;
  }

  // free memory local to this routine

  delete [] ecount;
  memory->destroy(elistfull);
  delete [] parent;
  delete [] mark;
  delete [] stack;

  return 0;
}

/* ----------------------------------------------------------------------
   compute evalues and evectors of 3x3 real symmetric matrix
   based on Jacobi rotations
   adapted from Numerical Recipes jacobi() function
------------------------------------------------------------------------- */

int FixPOEMS::jacobi(double **matrix, double *evalues, double **evectors)
{
  int i,j,k;
  double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) evectors[i][j] = 0.0;
    evectors[i][i] = 1.0;
  }
  for (i = 0; i < 3; i++) {
    b[i] = evalues[i] = matrix[i][i];
    z[i] = 0.0;
  }
  
  for (int iter = 1; iter <= MAXJACOBI; iter++) {
    sm = 0.0;
    for (i = 0; i < 2; i++)
      for (j = i+1; j < 3; j++)
	sm += fabs(matrix[i][j]);
    if (sm == 0.0) return 0;
    
    if (iter < 4) tresh = 0.2*sm/(3*3);
    else tresh = 0.0;
    
    for (i = 0; i < 2; i++) {
      for (j = i+1; j < 3; j++) {
	g = 100.0*fabs(matrix[i][j]);
	if (iter > 4 && fabs(evalues[i])+g == fabs(evalues[i])
	    && fabs(evalues[j])+g == fabs(evalues[j]))
	  matrix[i][j] = 0.0;
	else if (fabs(matrix[i][j]) > tresh) {
	  h = evalues[j]-evalues[i];
	  if (fabs(h)+g == fabs(h)) t = (matrix[i][j])/h;
	  else {
	    theta = 0.5*h/(matrix[i][j]);
	    t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt(1.0+t*t);
	  s = t*c;
	  tau = s/(1.0+c);
	  h = t*matrix[i][j];
	  z[i] -= h;
	  z[j] += h;
	  evalues[i] -= h;
	  evalues[j] += h;
	  matrix[i][j] = 0.0;
	  for (k = 0; k < i; k++) rotate(matrix,k,i,k,j,s,tau);
	  for (k = i+1; k < j; k++) rotate(matrix,i,k,k,j,s,tau);
	  for (k = j+1; k < 3; k++) rotate(matrix,i,k,j,k,s,tau);
	  for (k = 0; k < 3; k++) rotate(evectors,k,i,k,j,s,tau);
	}
      }
    }
    
    for (i = 0; i < 3; i++) {
      evalues[i] = b[i] += z[i];
      z[i] = 0.0;
    }
  }
  return 1;
}

/* ----------------------------------------------------------------------
   perform a single Jacobi rotation
------------------------------------------------------------------------- */

void FixPOEMS::rotate(double **matrix, int i, int j, int k, int l,
		      double s, double tau)
{
  double g = matrix[i][j];
  double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in segment frame
   set wbody component to 0.0 if inertia component is 0.0
     otherwise segment can spin easily around that axis
   project space-frame angular momentum onto segment axes
     and divide by principal moments
------------------------------------------------------------------------- */

void FixPOEMS::omega_from_mq(double *m, double *ex, double *ey, double *ez,
			     double *inertia, double *w)
{
  double wbody[3];

  if (inertia[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / inertia[0];
  if (inertia[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / inertia[1];
  if (inertia[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / inertia[2];

  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
}
/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid segment
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixPOEMS::set_xv()
{
  int isegment;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double vr[6];

  int *image = atom->image;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass; 
  double *mass = atom->mass; 
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  // set x and v of each atom
  // only set joint atoms for 1st rigid segment they belong to

  for (int i = 0; i < nlocal; i++) {
    if (natom2segment[i] == 0) continue;
    isegment = atom2segment[i][0];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    // save old positions and velocities for virial

    if (evflag) {
      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on segment orientation
    // v = vcm + omega around center-of-mass

    x[i][0] = ex_space[isegment][0]*displace[i][0] +
      ey_space[isegment][0]*displace[i][1] + 
      ez_space[isegment][0]*displace[i][2];
    x[i][1] = ex_space[isegment][1]*displace[i][0] +
      ey_space[isegment][1]*displace[i][1] + 
      ez_space[isegment][1]*displace[i][2];
    x[i][2] = ex_space[isegment][2]*displace[i][0] +
      ey_space[isegment][2]*displace[i][1] + 
      ez_space[isegment][2]*displace[i][2];

    v[i][0] = omega[isegment][1]*x[i][2] - omega[isegment][2]*x[i][1] +
      vcm[isegment][0];
    v[i][1] = omega[isegment][2]*x[i][0] - omega[isegment][0]*x[i][2] +
      vcm[isegment][1];
    v[i][2] = omega[isegment][0]*x[i][1] - omega[isegment][1]*x[i][0] +
      vcm[isegment][2];
    
    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    //fprintf(screen, "oriential: %g %g %g \n",
    //          ex_space[isegment][0],ex_space[isegment][1],ex_space[isegment][2]);

    x[i][0] += xcm[isegment][0] - xbox*xprd;
    x[i][1] += xcm[isegment][1] - ybox*yprd;
    x[i][2] += xcm[isegment][2] - zbox*zprd;

    // virial = unwrapped coords dotted into segment constraint force
    // segment constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to segment
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2]; 

      vr[0] = 0.5*fc0*x0;
      vr[1] = 0.5*fc1*x1;
      vr[2] = 0.5*fc2*x2;
      vr[3] = 0.5*fc1*x0;
      vr[4] = 0.5*fc2*x0;
      vr[5] = 0.5*fc2*x1;

      v_tally(1,&i,1.0,vr);
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid segment
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixPOEMS::set_v()
{
  int isegment;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double vr[6];

  double *rmass = atom->rmass; 
  double *mass = atom->mass; 
  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // set v of each atom
  // only set joint atoms for 1st rigid segment they belong to

  for (int i = 0; i < nlocal; i++) {
    if (natom2segment[i] == 0) continue;
    isegment = atom2segment[i][0];

    dx = ex_space[isegment][0]*displace[i][0] +
      ey_space[isegment][0]*displace[i][1] + 
      ez_space[isegment][0]*displace[i][2];
    dy = ex_space[isegment][1]*displace[i][0] +
      ey_space[isegment][1]*displace[i][1] + 
      ez_space[isegment][1]*displace[i][2];
    dz = ex_space[isegment][2]*displace[i][0] +
      ey_space[isegment][2]*displace[i][1] + 
      ez_space[isegment][2]*displace[i][2];

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[isegment][1]*dz - omega[isegment][2]*dy + vcm[isegment][0];
    v[i][1] = omega[isegment][2]*dx - omega[isegment][0]*dz + vcm[isegment][1];
    v[i][2] = omega[isegment][0]*dy - omega[isegment][1]*dx + vcm[isegment][2];

    // virial = unwrapped coords dotted into segment constraint force
    // segment constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to segment
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) 
     {
        massone = rmass[i];
     }
      else
     {
        massone = mass[type[i]];		
      }
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2]; 

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;

      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      vr[0] = 0.5*fc0*x0;
      vr[1] = 0.5*fc1*x1;
      vr[2] = 0.5*fc2*x2;
      vr[3] = 0.5*fc1*x0;
      vr[4] = 0.5*fc2*x0;
      vr[5] = 0.5*fc2*x1;

      v_tally(1,&i,1.0,vr);
    }
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixPOEMS::grow_arrays(int nmax)
{
  memory->grow(natom2segment,nmax,"fix_poems:natom2segment");
  memory->grow(atom2segment,nmax,MAXSEGMENT,"fix_poems:atom2segment");
  memory->grow(displace,nmax,3,"fix_poems:displace");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixPOEMS::copy_arrays(int i, int j)
{
  natom2segment[j] = natom2segment[i];
  for (int k = 0; k < natom2segment[j]; k++) atom2segment[j][k] = atom2segment[i][k];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixPOEMS::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*MAXSEGMENT * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixPOEMS::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = static_cast<double> (natom2segment[i]);
  for (int j = 0; j < natom2segment[i]; j++) 
    buf[m++] = static_cast<double> (atom2segment[i][j]);
  buf[m++] = displace[i][0];
  buf[m++] = displace[i][1];
  buf[m++] = displace[i][2];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixPOEMS::unpack_exchange(int nlocal, double *buf)
{

  int m = 0;
  natom2segment[nlocal] = static_cast<int> (buf[m++]);
  for (int i = 0; i < natom2segment[nlocal]; i++)
    atom2segment[nlocal][i] = static_cast<int> (buf[m++]);
  displace[nlocal][0] = buf[m++];
  displace[nlocal][1] = buf[m++];
  displace[nlocal][2] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixPOEMS::reset_dt()
{
  dtv = update->dt;  
  dtf = 0.5 * update->dt * force->ftm2v;  
  dthalf = 0.5 * update->dt;
}
