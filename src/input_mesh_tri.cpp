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
#include "memory.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "math.h"
#include "vector_liggghts.h"
#include "input_mesh_tri.h"
#include "tri_mesh.h"

using namespace LAMMPS_NS;

InputMeshTri::InputMeshTri(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
verbose_(false)
{}

InputMeshTri::~InputMeshTri()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile(const char *filename, class TriMesh *mesh,bool verbose)
{
  verbose_ = verbose;
  if(strlen(filename) < 5)
    error->all(FLERR,"Illegal command, file name too short for input of triangular mesh");
  const char *ext = &(filename[strlen(filename)-3]);

  // read file
  // case STL file or VTK file

  bool is_stl = (strcmp(ext,"stl") == 0) || (strcmp(ext,"STL") == 0);
  bool is_vtk = (strcmp(ext,"vtk") == 0) || (strcmp(ext,"VTK") == 0);

  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set stl___file

  if (me == 0)
  {
    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[128];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(FLERR,str);
    }
  } else nonlammps_file = NULL;

  if(is_stl)
  {
      if (comm->me == 0) fprintf(screen,"\nReading STL file '%s' \n",filename);
      meshtrifile_stl(mesh);
  }
  else if(is_vtk)
  {
      if (comm->me == 0) fprintf(screen,"\nReading VTK file '%s' \n",filename);
      meshtrifile_vtk(mesh);
  }
  else error->all(FLERR,"Illegal command, need either an STL file or a VTK file as input for triangular mesh.");

  if(nonlammps_file) fclose(nonlammps_file);
}

/* ----------------------------------------------------------------------
   process VTK file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_vtk(class TriMesh *mesh)
{
  int n,m;

  double **points = NULL;
  int ipoint = 0,npoints = 0;

  int **cells = NULL, *lines = NULL;
  int icell = 0,ncells = 0;

  int ntris = 0;
  int iLine = 0;

  int nLines = 0;

  while (1)
  {
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

    // parse one line from the file
    parse_nonlammps();
    // skip empty lines
    if(narg == 0){
         if (me == 0 && verbose_)
            fprintf(screen,"Note: Skipping empty line in VTK mesh file\n");
      continue;
    }

    //increase line counter
    iLine++;

    if(iLine < 3) continue;

    if(iLine == 3)
    {
        if(strcmp(arg[0],"ASCII"))
            error->all(FLERR,"Expecting ASCII VTK mesh file, cannot continue");
        continue;
    }

    if(iLine == 4)
    {
        if(strcmp(arg[0],"DATASET UNSTRUCTURED_GRID"))
            error->all(FLERR,"Expecting ASCII VTK unstructured grid mesh file, cannot continue");
        continue;
    }

    if(iLine == 5)
    {
        if(strcmp(arg[0],"POINTS"))
            error->all(FLERR,"Expecting 'POINTS' section in ASCII VTK mesh file, cannot continue");
        npoints = atoi(arg[1]);
        memory->create<double>(points,npoints,3,"input_mesh:points");
        continue;
    }

    if(iLine <= 5+npoints)
    {
        if(narg != 3)
            error->all(FLERR,"Expecting 3 values for each point in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        points[ipoint][0] = atof(arg[0]);
        points[ipoint][1] = atof(arg[1]);
        points[ipoint][2] = atof(arg[2]);
        ipoint++;
        continue;
    }

    if(iLine == 6+npoints)
    {
        if(strcmp(arg[0],"CELLS"))
            error->all(FLERR,"Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);
        memory->create<int>(cells,ncells,3,"input_mesh:cells");
        memory->create<int>(lines,ncells,"input_mesh:lines");
        continue;
    }

    //copy data of all which have 3 values - can be tri, polygon etc
    if(iLine <= 6+npoints+ncells)
    {
        if(narg == 4)
            for (int j=0;j<3;j++) cells[icell][j] = atoi(arg[1+j]);
        else
            cells[icell][0] = -1;

        lines[icell] = nLines;

        icell++;
        continue;
    }

    if(iLine == 7+npoints+ncells)
    {
        if(strcmp(arg[0],"CELL_TYPES"))
            error->all(FLERR,"Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))
            error->all(FLERR,"Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
         icell = 0;
        continue;
    }

    //only take triangles (cell type 5 according to VTK standard) - count them
    if(iLine <= 7+npoints+2*ncells)
    {
        if(strcmp(arg[0],"5")) cells[icell][0] = -1; //remove if not a tet
        else ntris++;
        icell++;
        continue;
    }

  }

  //now that everything is parsed, write the data into the mesh
  for(int i = 0; i < ncells; i++)
  {
      if(cells[i][0] == -1) continue;
      addTriangle(mesh,points[cells[i][0]],points[cells[i][1]],points[cells[i][2]],lines[i]);
  }

  memory->destroy<double>(points);
  memory->destroy<int>(cells);
}

/* ----------------------------------------------------------------------
   process STL file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_stl(class TriMesh *mesh)
{
  int n,m;
  int iVertex = 0;
  double vertices[3][3];
  bool insideSolidObject = false;
  bool insideFacet = false;
  bool insideOuterLoop = false;

  int nLines = 0, nLinesTri = 0;

  while (1)
  {
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

    // parse one line from the stl file
    parse_nonlammps();

    // skip empty lines
    if(narg==0){
         if (me == 0 && verbose_)
            fprintf(screen,"Note: Skipping empty line in STL file\n");
      continue;
    }

    // detect begin and end of a solid object, facet and vertices
    if (strcmp(arg[0],"solid") == 0)
    {
      if (insideSolidObject)
        error->all(FLERR,"Corrupt or unknown STL file: New solid object begins without closing prior solid object.");
      insideSolidObject=true;
      if (me == 0 && verbose_){
         fprintf(screen,"Solid body detected in STL file\n");
       }
    }
    else if (strcmp(arg[0],"endsolid") == 0)
    {
       if (!insideSolidObject)
         error->all(FLERR,"Corrupt or unknown STL file: End of solid object found, but no begin.");
       insideSolidObject=false;
       if (me == 0 && verbose_) {
         fprintf(screen,"End of solid body detected in STL file.\n");
       }
    }

    // detect begin and end of a facet within a solids object
    else if (strcmp(arg[0],"facet") == 0)
    {
      if (insideFacet)
        error->all(FLERR,"Corrupt or unknown STL file: New facet begins without closing prior facet.");
      if (!insideSolidObject)
        error->all(FLERR,"Corrupt or unknown STL file: New facet begins outside solid object.");
      insideFacet = true;

      nLinesTri = nLines;

      // check for keyword normal belonging to facet
      if (strcmp(arg[1],"normal") != 0)
        error->all(FLERR,"Corrupt or unknown STL file: Facet normal not defined.");

      // do not import facet normal (is calculated later)
    }
    else if (strcmp(arg[0],"endfacet") == 0)
    {
       if (!insideFacet)
         error->all(FLERR,"Corrupt or unknown STL file: End of facet found, but no begin.");
       insideFacet = false;
       if (iVertex != 3)
         error->all(FLERR,"Corrupt or unknown STL file: Number of vertices not equal to three (no triangle).");

      // add triangle to mesh
      //printVec3D(screen,"vertex",vertices[0]);
      //printVec3D(screen,"vertex",vertices[1]);
      //printVec3D(screen,"vertex",vertices[2]);
      addTriangle(mesh,vertices[0],vertices[1],vertices[2],nLinesTri);

       if (me == 0) {
         //fprintf(screen,"  End of facet detected in in solid body.\n");
       }
    }

    //detect begin and end of an outer loop within a facet
    else if (strcmp(arg[0],"outer") == 0)
    {
      if (insideOuterLoop)
        error->all(FLERR,"Corrupt or unknown STL file: New outer loop begins without closing prior outer loop.");
      if (!insideFacet)
        error->all(FLERR,"Corrupt or unknown STL file: New outer loop begins outside facet.");
      insideOuterLoop = true;
      iVertex = 0;

      if (me == 0){
         //fprintf(screen,"    Outer loop detected in facet.\n");
       }
    }
    else if (strcmp(arg[0],"endloop") == 0)
    {
       if (!insideOuterLoop)
         error->all(FLERR,"Corrupt or unknown STL file: End of outer loop found, but no begin.");
       insideOuterLoop=false;
       if (me == 0) {
         //fprintf(screen,"    End of outer loop detected in facet.\n");
       }
    }

    else if (strcmp(arg[0],"vertex") == 0)
    {
       if (!insideOuterLoop)
         error->all(FLERR,"Corrupt or unknown STL file: Vertex found outside a loop.");

       if (me == 0) {
         //fprintf(screen,"      Vertex found.\n");
       }

      // read the vertex
      for (int j=0;j<3;j++)
        vertices[iVertex][j]=atof(arg[1+j]);

      iVertex++;
      if (iVertex > 3)
         error->all(FLERR,"Corrupt or unknown STL file: Can not have more than 3 vertices "
                          "in a facet (only triangular meshes supported).");
    }
  }
}

/* ----------------------------------------------------------------------
   add a triangle to the mesh
------------------------------------------------------------------------- */

void InputMeshTri::addTriangle(TriMesh *mesh,double *a, double *b, double *c,int lineNumber)
{
    double **nodeTmp = create<double>(nodeTmp,3,3);
    for(int i=0;i<3;i++){
      nodeTmp[0][i] = a[i];
      nodeTmp[1][i] = b[i];
      nodeTmp[2][i] = c[i];
    }
    mesh->addElement(nodeTmp,lineNumber);
    destroy<double>(nodeTmp);
}
