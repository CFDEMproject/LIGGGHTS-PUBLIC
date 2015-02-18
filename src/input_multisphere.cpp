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

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "style_command.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "force.h"
#include "pair.h"
#include "min.h"
#include "modify.h"
#include "compute.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "special.h"
#include "variable.h"
#include "error.h"
#include "memory.h"
#include "input_multisphere.h"

using namespace LAMMPS_NS;

#define MAXLINE 2048
#define DELTA 4

InputMultisphere::InputMultisphere(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv)
{}

InputMultisphere::~InputMultisphere()
{}

/* ----------------------------------------------------------------------
   process clump file
------------------------------------------------------------------------- */

int InputMultisphere::clmpfile(double **xclmp,double *rclmp,int *atomtypeclmp,int nclmps)
{
  int n,m;
  int iClmp = 0;

  while (1) {
    // read a line from input script
    // n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line

    if (me == 0) {
      m = 0;
      while (1) {
        if (maxline-m < 2) reallocate(line,maxline,0);
        if (fgets(&line[m],maxline-m,nonlammps_file) == NULL) {
          if (m) n = strlen(line) + 1;
          else n = 0;
          break;
        }
        m = strlen(line);
        if (line[m-1] != '\n') continue;

        m--;
        while (m >= 0 && isspace(line[m])) m--;
        if (m < 0 || line[m] != '&') {
          line[m+1] = '\0';
          n = m+2;
          break;
        }
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      break;
    }

    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    //parse one line from the clump file
    parse_nonlammps();

    //skip empty lines
    if(narg == 0)
    {
        if (me == 0) fprintf(screen,"Note: Skipping empty line or comment line in clump file\n");
        continue;
    }

    if(iClmp >= nclmps)
        error->all(FLERR,"Number of clumps in file larger than number specified");

    if((0 == atomtypeclmp) && (narg < 4))
        error->all(FLERR,"Not enough arguments in one line of clump file, need to specify "
                         "[xcoo ycoo zcoo radius] in each line");

    if(atomtypeclmp && (narg < 5))
        error->all(FLERR,"Not enough arguments in one line of clump file, need to specify "
                         "[xcoo ycoo zcoo radius type] in each line");

    rclmp[iClmp] = atof(arg[3]);
    if(atomtypeclmp)
        atomtypeclmp[iClmp] = atoi(arg[4]);

    for(int j = 0; j < 3; j++)
       xclmp[iClmp][j] = atof(arg[j]);

    iClmp++;
  }

  return iClmp;
}

/* ----------------------------------------------------------------------
   process all input from file
------------------------------------------------------------------------- */

void InputMultisphere::clmpfile(const char *filename, double **xclmp,double *rclmp,int *atomtypeclmp,int nclmps)
{
  if (me == 0)
  {
    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL)
    {
      char str[128];
      sprintf(str,"Cannot open clump file %s",filename);
      error->one(FLERR,str);
    }
  }
  else nonlammps_file = NULL;

  if(clmpfile(xclmp,rclmp,atomtypeclmp,nclmps) != nclmps)
    error->all(FLERR,"Number of clumps in file does not match number of clumps that were specified");

  if(nonlammps_file) fclose(nonlammps_file);

}
