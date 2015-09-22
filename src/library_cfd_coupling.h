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
int liggghts_get_maxtag_ms(void *ptr);
int liggghts_get_ntypes_ms(void *ptr);
double* liggghts_get_vclump_ms(void *ptr);
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
