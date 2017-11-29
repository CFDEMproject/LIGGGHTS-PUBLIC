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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctype.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include <cmath>
#include "vector_liggghts.h"
#include "input_mesh_tet.h"
#include "region_mesh_tet.h"

using namespace LAMMPS_NS;

#define MAXLINE 2048
#define DELTA 4

InputMeshTet::InputMeshTet(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
verbose_(false)
{}

InputMeshTet::~InputMeshTet()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshTet::meshtetfile(const char *filename, class RegTetMesh *mesh, bool verbose)
{
  verbose_ = verbose;

  if(strlen(filename) < 5)
    error->all(FLERR,"Illegal command, file name too short for input of tet mesh");

  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set stl___file

  if (me == 0) {

    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[512];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(FLERR,str);
    }
  } else nonlammps_file = NULL;

  meshtetfile_vtk(mesh);

  if(nonlammps_file) fclose(nonlammps_file);

}

/* ----------------------------------------------------------------------
   process VTK file
------------------------------------------------------------------------- */

void InputMeshTet::meshtetfile_vtk(class RegTetMesh *mesh)
{
  int n,m;

  double phix = (mesh->rot_angle[0])*M_PI/180.;
  double phiy = (mesh->rot_angle[1])*M_PI/180.;
  double phiz = (mesh->rot_angle[2])*M_PI/180.;

  int flag_outside = 0;

  double **points = NULL;
  int ipoint = 0, npoints = 0;
  double vert_before_rot[3], vert_after_rot[3];

  int **cells = NULL;
  int icell = 0, ncells = 0;

  int ntets = 0;
  int iLine = 0;

  int nLines = 0;

  int flag_other_than_tet = 0;

  bool allPointsRead = false, allCellsRead = false;
  int lastPointLine = 0, lastCellLine = 0;

  if (!domain->box_exist)
      error->all(FLERR,"Cannot load mesh before simulation box is defined.");

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

    // lines start with 1 (not 0)
    nLines++;

    //parse one line from the file
    parse_nonlammps();

    // skip empty lines
    if(narg == 0){
         if (me == 0 && verbose_)
            fprintf(screen,"Note: Skipping empty line in VTK mesh file\n");
      continue;
    }

    //increase line counter
    iLine++;

    if(iLine < 2) continue;

    if(iLine == 2)
    {
        if(strcmp(arg[0],"ASCII")) error->all(FLERR,"Expecting ASCII VTK mesh file, cannot continue");
        continue;
    }

    if(iLine == 3)
    {
        if(strcmp(arg[0],"DATASET") || strcmp(arg[1],"UNSTRUCTURED_GRID")) error->all(FLERR,"Expecting ASCII VTK unstructured grid mesh file, cannot continue");
        continue;
    }

    if(iLine == 4)
    {
        if(strcmp(arg[0],"POINTS")) error->all(FLERR,"Expecting 'POINTS' section in ASCII VTK mesh file, cannot continue");
        npoints = atoi(arg[1]);

        memory->create(points,npoints,3,"input_mesh:points");
        continue;
    }

    if (!allPointsRead)
    {
        if(narg%3 != 0) error->all(FLERR,"Expecting 3 values for each point in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        for (int p=0;p<narg;p=p+3) {
            //read the vertex, translate and scale it
            for (int j=0;j<3;j++) vert_before_rot[j]=(atof(arg[p+j])+(mesh->off_fact[j]))*(mesh->scale_fact);

            //rotate the vertex
            vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
            vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
            vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

            if (!domain->is_in_domain(vert_after_rot))
                flag_outside = 1;

            //store the vertex
            vectorCopy3D(vert_after_rot,points[ipoint]);
            ipoint++;
        }
        if (ipoint == npoints) {
            allPointsRead = true;
            lastPointLine = iLine;
        }

        continue;
    }

    if(allPointsRead && iLine == lastPointLine+1)
    {
        if(strcmp(arg[0],"CELLS")) error->all(FLERR,"Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);

        memory->create(cells,ncells,4,"input_mesh:cells");
        continue;
    }

    //copy data of all which have 4 values - can be tet, quad, poly_line or triangle_strip
    if(!allCellsRead)
    {
        if(narg == 5)
        {
            for (int j=0;j<4;j++) cells[icell][j] = atoi(arg[1+j]);

        }
        else
        {
            cells[icell][0] = -1;
        }

        icell++;
        if (icell == ncells) {
            allCellsRead = true;
            lastCellLine = iLine;
        }

        continue;
    }

    if(allCellsRead && iLine == lastCellLine+1)
    {
        if(strcmp(arg[0],"CELL_TYPES")) error->all(FLERR,"Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))  error->all(FLERR,"Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        icell = 0;

        continue;
    }

    //only take tetraeders (cell type 10 according to VTK standard) - count them
    if(allCellsRead && iLine <= lastCellLine+1+ncells)
    {
        if(strcmp(arg[0],"10"))
        {
            cells[icell][0] = -1; //remove if not a tet
            flag_other_than_tet = 1;
        }
        else ntets++;
        icell++;
        continue;
    }
}

  //must throw an error here since regions may not extend outside box
  if(flag_outside)
    error->all(FLERR,"VTK mesh file is incompatible with simulation box: One or more vertices outside simulation box");

  //now that everything is parsed, write the data into the mesh
  while(mesh->nTetMax < ntets) mesh->grow_arrays();

  double **tetnodes;
  memory->create(tetnodes,4,3,"input_mesh:tetnodes");

  for(int i = 0; i < ncells; i++)
  {
      if(cells[i][0] == -1) continue;
      for (int j = 0; j < 4; j++)
        vectorCopy3D(points[cells[i][j]],tetnodes[j]);

      mesh->add_tet(tetnodes);
  }

  if(flag_other_than_tet && comm->me == 0)
    error->warning(FLERR,"VTK file contains other types than tetrahedra - only tets are currently supported in LIGGGHTS, other cells are discarded");

  if(ncells == 0)
    error->all(FLERR,"VTK mesh file containing no tet cell - cannnot continue");

  memory->destroy(tetnodes);
  memory->destroy(points);
  memory->destroy(cells);
}

