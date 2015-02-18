/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program - DROPPED*/

#ifdef __cplusplus
//extern "C" {
#endif

int liggghts_get_maxtag(void *ptr);

void* locate_coupling_fix(void *ptr);
void data_liggghts_to_of(char *name,char *type,void *ptr,void *&data,char *datatype);
void data_of_to_liggghts(char *name,char *type,void *ptr,void *data,char *datatype);
void update_rm(void *ptr);
void check_datatransfer(void *ptr);

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr);
void allocate_external_int(int    **&data, int len2,char *,  int    initvalue,void *ptr);

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr);
void allocate_external_double(double **&data, int len2,char *,  double initvalue,void *ptr);

#ifdef __cplusplus
//}
#endif
