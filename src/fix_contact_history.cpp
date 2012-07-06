/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "fix_contact_history.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "force.h"
#include "update.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA_MAXTOUCH_PAIR 15
#define DELTA_MAXTOUCH_MESH 3

/* ---------------------------------------------------------------------- */

FixContactHistory::FixContactHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //parse args
  if(narg < 6)
    error->fix_error(FLERR,this,"not enough parameters");

  int iarg = 3;

  if(strcmp(arg[iarg],"pair") == 0)
    is_pair = true;
  else if(strcmp(arg[iarg],"mesh") == 0)
    is_pair = false;
  else
    error->fix_error(FLERR,this,"expecting keyword 'pair' or 'mesh'");
  iarg++;

  // for mesh style, read mesh
  if(!is_pair)
  {
      Fix *f = modify->find_fix_id(arg[iarg++]);
      if(!f || strncmp(f->style,"mesh/surface",12) )
        error->fix_error(FLERR,this,"wrong ID for fix mesh/surface");
      mesh_ = (static_cast<FixMeshSurface*>(f))->triMesh();
  }

  //read dnum
  dnum = atoi(arg[iarg++]);
  
  if(dnum < 1)
    error->all(FLERR,"dnum must be >=1 in fix contacthistory");
  if(dnum > 10)
    error->warning(FLERR,"dnum >10 in fix contacthistory - are you really sure you intend this?");

  newtonflag = NULL;
  history_id = NULL;
  if(is_pair)
  {
      //read newtonflag
      if(narg-iarg < 2*dnum)
        error->all(FLERR,"Illegal fix contacthistory command - not enough parameters (need to specify an id and a newtonflag for each dnum)");

      newtonflag = new int[dnum];
      //history_id = new char*[dnum];
      history_id = (char**) memory->smalloc((dnum)*sizeof(char*),"FixContactHistory:history_id");

      for(int i = 0 ; i < dnum; i++)
      {
        
        history_id[i] = new char[strlen(arg[iarg])+1];
        strcpy(history_id[i],arg[iarg++]);
        newtonflag[i] = atoi(arg[iarg++]);
        if(newtonflag[i] != 0 && newtonflag[i] != 1)
            error->all(FLERR,"Illegal fix history command - newtonflag must be either 0 or 1");
      }
  }

  restart_peratom = 1;
  restart_global = 1; 
  create_attribute = 1;

  // set time_depend so that history will be preserved correctly
  // across multiple runs via laststep setting in granular pair styles
  time_depend = 1;

  if(is_pair)
    maxtouch = DELTA_MAXTOUCH_PAIR;
  else
    maxtouch = DELTA_MAXTOUCH_MESH;

  // perform initial allocation of atom-based arrays
  // register with atom class

  atom->add_callback(0);
  atom->add_callback(1);

  npartner = NULL;
  partner = NULL;
  contacthistory = NULL;
  delflag = NULL;
  pair_gran = NULL;

  // initialize npartner to 0 so neighbor list creation is OK the 1st time
  if(atom->nmax > 0)
  {
      grow_arrays(atom->nmax);
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++) npartner[i] = 0;
  }

  laststep = -1;
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::post_create()
{
    if( ( (strcmp(style,"contacthistory") == 0) && ! is_pair ) ||
        ( (strcmp(style,"contacthistory/mesh") == 0) &&  is_pair ) )
        error->fix_error(FLERR,this,"illegal call");
}

/* ---------------------------------------------------------------------- */

FixContactHistory::~FixContactHistory()
{
    // unregister this fix so atom class doesn't invoke it any more
    atom->delete_callback(id,0);
    atom->delete_callback(id,1);

    // delete locally stored arrays
    memory->destroy(npartner);
    memory->destroy(partner);
    memory->destroy(contacthistory);
    memory->destroy(delflag);

    if(is_pair)
    {
        delete[]newtonflag;

        for(int i = 0; i < dnum; i++)
            delete [](history_id[i]);

        memory->sfree(history_id);
    }
}

/* ---------------------------------------------------------------------- */

int FixContactHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::init()
{
    if (atom->tag_enable == 0)
      error->fix_error(FLERR,this,"using contact history requires atoms have IDs");

    if(is_pair)
    {
        if(!force->pair_match("gran", 0))
            error->fix_error(FLERR,this,"Please use a granular pair style for fix contacthistory");
        pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
    }
    else if(!mesh_)
        error->fix_error(FLERR,this,"illegal");
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::setup_pre_exchange()
{
  if(laststep > 0) pre_exchange();
}

/* ----------------------------------------------------------------------
   pre_exchange - decide if pair or mesh
------------------------------------------------------------------------- */

void FixContactHistory::pre_exchange()
{
    
    if(is_pair)
        pre_exchange_pair();

    check_grow();

    laststep = update->ntimestep;
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be exchanged with atoms
------------------------------------------------------------------------- */

void FixContactHistory::pre_exchange_pair()
{
  int i,j,ii,jj,m,inum,jnum,d;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *hist,*allhist,**firsthist;

  // zero npartners for all current atoms
  // do not do this if inum = 0, since will wipe out npartner counts
  // inum is 0 when pre_exchange() is called from integrate->setup()
  //   before first run when no neighbor lists yet exist
  //   partner info may have just been read from restart file
  //     and won't be used until granular neighbor lists are built
  //   if nothing read from restart file, constructor sets npartner = 0

  NeighList *list = pair_gran->list;
  int nlocal = atom->nlocal;

  if (list->inum)
    for (i = 0; i < nlocal; i++)
      npartner[i] = 0;

  // copy contact info from neighbor list atoms to atom arrays

  int *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listgranhistory->firstneigh;
  firsthist = list->listgranhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allhist = firsthist[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        hist = &allhist[dnum*jj]; 
        j = jlist[jj];
            j &= NEIGHMASK;

        if (npartner[i] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 
        m = npartner[i];
        partner[i][m] = tag[j];
        for (d = 0; d < dnum; d++) {
            contacthistory[i][m][d] = hist[d]; 
        }
        npartner[i]++;
        if (j < nlocal) {
            if (npartner[j] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 

            m = npartner[j];
            partner[j][m] = tag[i];
            for (d = 0; d < dnum; d++) {
               if(newtonflag[d]) contacthistory[j][m][d] = -hist[d]; 
               else              contacthistory[j][m][d] =  hist[d]; 
            }
          npartner[j]++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   look if other procs grew, and if so, do so too
------------------------------------------------------------------------- */

void FixContactHistory::check_grow()
{
  int maxtouch_all;
  MPI_Allreduce(&maxtouch,&maxtouch_all,1,MPI_INT,MPI_MAX,world);
  while (maxtouch<maxtouch_all) grow_arrays_maxtouch(atom->nmax);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixContactHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*maxtouch * sizeof(int);  
  bytes += nmax*maxtouch * dnum * sizeof(double); 
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistory::grow_arrays(int nmax)
{
  
  memory->grow(npartner,nmax,"contacthistory:npartner");
  memory->grow(partner,nmax,maxtouch,"contacthistory:partner");
  memory->grow(contacthistory,nmax,maxtouch,dnum,"contact_history:contacthistory");
  memory->grow(delflag,nmax,maxtouch,"contact_history:delflag");
}

/* ----------------------------------------------------------------------
   grow local atom-based arrays in case maxtouch is too small
------------------------------------------------------------------------- */

void FixContactHistory::grow_arrays_maxtouch(int nmax)
{
  if(comm->me==0)
  {
      if(screen) fprintf(screen,  "INFO: more than %d touching neighbor %s found, growing contact history.\n",maxtouch,is_pair?"atoms":"mesh elements");
      if(logfile) fprintf(logfile,"INFO: more than %d touching neighbor %s found, growing contact history.\n",maxtouch,is_pair?"atoms":"mesh elements");
  }

  int delta;
  if(is_pair)
    delta = DELTA_MAXTOUCH_PAIR;
  else
    delta = DELTA_MAXTOUCH_MESH;

  int **partner_g;
  memory->create(partner_g,nmax,maxtouch+delta,"contacthistory:partner_g");
  double ***contacthistory_g;
  memory->create(contacthistory_g,nmax,maxtouch+delta,dnum,"contacthistory:contacthistory_g");
  bool **delflag_g;
  memory->create(delflag_g,nmax,maxtouch+delta,"contacthistory:delflag_g");

  for (int i = 0; i < nmax; i++)
  {
      for (int j = 0; j < maxtouch; j++)
      {
          partner_g[i][j] = partner[i][j];
          delflag_g[i][j] = delflag[i][j];
          for (int k = 0 ; k < dnum; k++)
            contacthistory_g[i][j][k] = contacthistory[i][j][k];
      }

  }

  maxtouch += delta;

  int **h1; double ***h2; bool **h3;
  h1 = partner;
  h2 = contacthistory;
  h3 = delflag;
  partner = partner_g;
  contacthistory = contacthistory_g;
  delflag = delflag_g;
  memory->destroy(h1);
  memory->destroy(h2);
  memory->destroy(h3);
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays from i to j
------------------------------------------------------------------------- */

void FixContactHistory::copy_arrays(int i, int j)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    delflag[j][m] = delflag[i][m];
    for (int d = 0; d < dnum; d++) {
      contacthistory[j][m][d] = contacthistory[i][m][d];
    }
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixContactHistory::set_arrays(int i)
{
    
    npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = static_cast<double>(delflag[i][n]);
    for (int d = 0; d < dnum; d++) {
      buf[m++] = contacthistory[i][n][d];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    delflag[nlocal][n] = static_cast<bool> (buf[m++]);
    for (int d = 0; d < dnum; d++) {
      contacthistory[nlocal][n][d] = buf[m++];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixContactHistory::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = static_cast<double>(dnum);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  pre_exchange();
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixContactHistory::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  int pack_dnum = static_cast<int> (list[n++]);
  if(pack_dnum != dnum)
    error->all(FLERR,"Saved simulation state used different contact history model - can not restart");

}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixContactHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = (dnum+1)*npartner[i] + 2;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    for (int d = 0; d < dnum; d++) {
      buf[m++] = contacthistory[i][n][d];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  int d;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    delflag[nlocal][n] = false;
    for (d = 0; d < dnum; d++) {
      contacthistory[nlocal][n][d] = extra[nlocal][m++];
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixContactHistory::maxsize_restart()
{
  return (dnum+1)*maxtouch + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixContactHistory::size_restart(int nlocal)
{
  return (dnum+1)*npartner[nlocal] + 2;
}

/* ----------------------------------------------------------------------
   clear all entries where deleteFlag is set
------------------------------------------------------------------------- */

void FixContactHistory::cleanUpContacts()
{
    int nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
    {
        int j = 0;
        while(j < npartner[i])
        {
            // delete and copy values
            if(delflag[i][j])
            {
                // swap with last contact if more than 1
                if(npartner[i] > 1)
                {
                    partner[i][j] = partner[i][npartner[i]-1];
                    delflag[i][j] = delflag[i][npartner[i]-1];

                    if(dnum)
                        for(int d = 0; d < dnum; d++)
                            contacthistory[i][j][d] = contacthistory[i][npartner[i]-1][d];
                }

                partner[i][npartner[i]-1] = -1;
                delflag[i][npartner[i]-1] = false;

                if(dnum)
                    for(int d = 0; d < dnum; d++)
                        contacthistory[i][j][d] = 0.;

                npartner[i]--;
            }
            else j++;
        }
    }
}
