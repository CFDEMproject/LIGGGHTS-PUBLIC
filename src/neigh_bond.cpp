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

#include "neighbor.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BONDDELTA 10000

// bondlist, anglelist, dihedrallist, improperlist
//   no longer store atom->map() of the bond partners
// instead store domain->closest_image() of the bond partners of atom I
// this enables distances between list atoms to be calculated
//   w/out invoking domain->minimium_image(), e.g. in bond->compute()

/* ---------------------------------------------------------------------- */

void Neighbor::bond_all()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  double ***bond_hist = atom->bond_hist;  
  int *tag = atom->tag;
  int newton_bond = force->newton_bond;
  int n_bondhist = atom->n_bondhist;  

  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
        char str[128];
        sprintf(str,
                "Bond atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
                tag[i],bond_atom[i][m],me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      if (newton_bond || i < atom1) {
        if (nbondlist == maxbond) {
          maxbond += BONDDELTA;
          memory->grow(bondlist,maxbond,4,"neighbor:bondlist");  
          if(atom->n_bondhist)
            memory->grow(bondhistlist,maxbond,atom->n_bondhist,"neighbor:bondhistlist");  
        }
        bondlist[nbondlist][0] = i;
        bondlist[nbondlist][1] = atom1;
        bondlist[nbondlist][2] = bond_type[i][m];
        bondlist[nbondlist][3] = 0; 
        if(n_bondhist) { 
            for(int j = 0; j < n_bondhist; j++)
            {
                bondhistlist[nbondlist][j] = bond_hist[i][m][j];
                
            }
        }
        nbondlist++;
      }
    }
  if (cluster_check) bond_check();
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_partial()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  double ***bond_hist = atom->bond_hist; 
  int *tag = atom->tag;
  int newton_bond = force->newton_bond;
  int n_bondhist = atom->n_bondhist; 

  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_type[i][m] <= 0) continue;
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
        char str[128];
        sprintf(str,
                "Bond atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
                tag[i],bond_atom[i][m],me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      if (newton_bond || i < atom1) {
        if (nbondlist == maxbond) {
          maxbond += BONDDELTA;
          memory->grow(bondlist,maxbond,3,"neighbor:bondlist");
        }
        bondlist[nbondlist][0] = i;
        bondlist[nbondlist][1] = atom1;
        bondlist[nbondlist][2] = bond_type[i][m];
        bondlist[nbondlist][3] = 0; 
        if(n_bondhist) { 
            for(int j = 0; j < n_bondhist; j++)
                bondhistlist[nbondlist][j] = bond_hist[i][m][j];
        }
        nbondlist++;
      }
    }
  if (cluster_check) bond_check();
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_check()
{
  int i,j;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  for (int m = 0; m < nbondlist; m++) {
    i = bondlist[m][0];
    j = bondlist[m][1];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Bond extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_all()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        char str[128];
        sprintf(str,
                "Angle atoms %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
        if (nanglelist == maxangle) {
          maxangle += BONDDELTA;
          memory->grow(anglelist,maxangle,4,"neighbor:anglelist");
        }
        anglelist[nanglelist][0] = atom1;
        anglelist[nanglelist][1] = atom2;
        anglelist[nanglelist][2] = atom3;
        anglelist[nanglelist][3] = angle_type[i][m];
        nanglelist++;
      }
    }
  if (cluster_check) angle_check();
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_partial()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      if (angle_type[i][m] <= 0) continue;
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        char str[128];
        sprintf(str,
                "Angle atoms %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
        if (nanglelist == maxangle) {
          maxangle += BONDDELTA;
          memory->grow(anglelist,maxangle,4,"neighbor:anglelist");
        }
        anglelist[nanglelist][0] = atom1;
        anglelist[nanglelist][1] = atom2;
        anglelist[nanglelist][2] = atom3;
        anglelist[nanglelist][3] = angle_type[i][m];
        nanglelist++;
      }
    }
  if (cluster_check) angle_check();
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_check()
{
  int i,j,k;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  // check all 3 distances
  // in case angle potential computes any of them

  for (int m = 0; m < nanglelist; m++) {
    i = anglelist[m][0];
    j = anglelist[m][1];
    k = anglelist[m][2];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Angle extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        char str[128];
        sprintf(str,
                "Dihedral atoms %d %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                dihedral_atom1[i][m],dihedral_atom2[i][m],
                dihedral_atom3[i][m],dihedral_atom4[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (ndihedrallist == maxdihedral) {
          maxdihedral += BONDDELTA;
          memory->grow(dihedrallist,maxdihedral,5,"neighbor:dihedrallist");
        }
        dihedrallist[ndihedrallist][0] = atom1;
        dihedrallist[ndihedrallist][1] = atom2;
        dihedrallist[ndihedrallist][2] = atom3;
        dihedrallist[ndihedrallist][3] = atom4;
        dihedrallist[ndihedrallist][4] = dihedral_type[i][m];
        ndihedrallist++;
      }
    }
  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      if (dihedral_type[i][m] <= 0) continue;
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        char str[128];
        sprintf(str,
                "Dihedral atoms %d %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                dihedral_atom1[i][m],dihedral_atom2[i][m],
                dihedral_atom3[i][m],dihedral_atom4[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (ndihedrallist == maxdihedral) {
          maxdihedral += BONDDELTA;
          memory->grow(dihedrallist,maxdihedral,5,"neighbor:dihedrallist");
        }
        dihedrallist[ndihedrallist][0] = atom1;
        dihedrallist[ndihedrallist][1] = atom2;
        dihedrallist[ndihedrallist][2] = atom3;
        dihedrallist[ndihedrallist][3] = atom4;
        dihedrallist[ndihedrallist][4] = dihedral_type[i][m];
        ndihedrallist++;
      }
    }
  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_check(int nlist, int **list)
{
  int i,j,k,l;
  double dx,dy,dz,dxstart,dystart,dzstart;

  double **x = atom->x;
  int flag = 0;

  // check all 6 distances
  // in case dihedral/improper potential computes any of them

  for (int m = 0; m < nlist; m++) {
    i = list[m][0];
    j = list[m][1];
    k = list[m][2];
    l = list[m][3];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[l][0];
    dystart = dy = x[i][1] - x[l][1];
    dzstart = dz = x[i][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[l][0];
    dystart = dy = x[j][1] - x[l][1];
    dzstart = dz = x[j][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[k][0] - x[l][0];
    dystart = dy = x[k][1] - x[l][1];
    dzstart = dz = x[k][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Dihedral/improper extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        char str[128];
        sprintf(str,
                "Improper atoms %d %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                improper_atom1[i][m],improper_atom2[i][m],
                improper_atom3[i][m],improper_atom4[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain-> closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (nimproperlist == maximproper) {
          maximproper += BONDDELTA;
          memory->grow(improperlist,maximproper,5,"neighbor:improperlist");
        }
        improperlist[nimproperlist][0] = atom1;
        improperlist[nimproperlist][1] = atom2;
        improperlist[nimproperlist][2] = atom3;
        improperlist[nimproperlist][3] = atom4;
        improperlist[nimproperlist][4] = improper_type[i][m];
        nimproperlist++;
      }
    }
  if (cluster_check) dihedral_check(nimproperlist,improperlist);
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      if (improper_type[i][m] <= 0) continue;
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        char str[128];
        sprintf(str,
                "Improper atoms %d %d %d %d missing on proc %d at step "
                BIGINT_FORMAT,
                improper_atom1[i][m],improper_atom2[i][m],
                improper_atom3[i][m],improper_atom4[i][m],
                me,update->ntimestep);
        error->one(FLERR,str);
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      if (newton_bond ||
          (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
        if (nimproperlist == maximproper) {
          maximproper += BONDDELTA;
          memory->grow(improperlist,maximproper,5,"neighbor:improperlist");
        }
        improperlist[nimproperlist][0] = atom1;
        improperlist[nimproperlist][1] = atom2;
        improperlist[nimproperlist][2] = atom3;
        improperlist[nimproperlist][3] = atom4;
        improperlist[nimproperlist][4] = improper_type[i][m];
        nimproperlist++;
      }
    }
  if (cluster_check) dihedral_check(nimproperlist,improperlist);
}
