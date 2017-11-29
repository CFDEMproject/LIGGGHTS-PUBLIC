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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Arno Mayrhofer (CFDEMresearch GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ctype.h"
#include "memory.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "region.h"
#include "domain.h"
#include <cmath>
#include "vector_liggghts.h"
#include "input_mesh_tri.h"
#include "tri_mesh.h"

using namespace LAMMPS_NS;

InputMeshTri::InputMeshTri(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
verbose_(false),
i_exclusion_list_(0),
size_exclusion_list_(0),
exclusion_list_(0)
{}

InputMeshTri::~InputMeshTri()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile(const char *filename, class TriMesh *mesh,bool verbose,
                               const int size_exclusion_list, int *exclusion_list,
                               class Region *region)
{
  verbose_ = verbose;
  size_exclusion_list_ = size_exclusion_list;
  exclusion_list_ = exclusion_list;
  
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
      char str[512];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(FLERR,str);
    }
  } else nonlammps_file = NULL;

  if(is_stl)
  {
      if (comm->me == 0) fprintf(screen,"\nReading STL file '%s' (mesh processing step 1/3) \n",filename);
      meshtrifile_stl(mesh,region,filename);
      
  }
  else if(is_vtk)
  {
      if (comm->me == 0) fprintf(screen,"\nReading VTK file '%s' (mesh processing step 1/3) \n",filename);
      meshtrifile_vtk(mesh,region);
  }
  else error->all(FLERR,"Illegal command, need either an STL file or a VTK file as input for triangular mesh.");

  if(nonlammps_file) fclose(nonlammps_file);
}

/* ----------------------------------------------------------------------
   process VTK file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_vtk(class TriMesh *mesh,class Region *region)
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

    if (0 == (nLines % 100000) && comm->me == 0)
        fprintf(screen,"   successfully read a chunk of 100000 elements in VTK file\n");

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
      if(size_exclusion_list_ > 0 && lines[i] == exclusion_list_[i_exclusion_list_])
      {
         if(i_exclusion_list_ < size_exclusion_list_-1)
            i_exclusion_list_++;
         continue;
      }
      if(!region || ( region->match(points[cells[i][0]]) && region->match(points[cells[i][1]]) && region->match(points[cells[i][2]]) ) )
      addTriangle(mesh,points[cells[i][0]],points[cells[i][1]],points[cells[i][2]],lines[i]);
  }

  memory->destroy<double>(points);
  memory->destroy<int>(cells);
}

/* ----------------------------------------------------------------------
   process STL file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_stl(class TriMesh *mesh,class Region *region, const char *filename)
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

    if (0 == (nLines % 100000) && comm->me == 0)
        fprintf(screen,"   successfully read a chunk of 100000 elements in STL file '%s'\n",filename);

    // parse one line from the stl file

    parse_nonlammps();

    // skip empty lines
    if(narg==0){
         if (me == 0 && verbose_)
            fprintf(screen,"Note: Skipping empty line in STL file\n");
      continue;
    }

    if (strcmp(arg[0],"solid") != 0 && nLines == 1)
    {
        if (me == 0)
        {
            fclose(nonlammps_file);
            if (verbose_)
                fprintf(screen,"Note: solid keyword not found, assuming binary stl file\n");
        }
        nonlammps_file = NULL;
        meshtrifile_stl_binary(mesh, region, filename);
        break;
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

      // nLinesTri is the line
      
      if(size_exclusion_list_ > 0 && nLinesTri == exclusion_list_[i_exclusion_list_])
      {
         
         if(i_exclusion_list_ < size_exclusion_list_-1)
         {
            i_exclusion_list_++;
            
            while((exclusion_list_[i_exclusion_list_-1] == exclusion_list_[i_exclusion_list_]) && (i_exclusion_list_ < size_exclusion_list_-1))
                i_exclusion_list_++;
         }
      }
      else if(!region || ( region->match(vertices[0]) && region->match(vertices[1]) && region->match(vertices[2]) ) )
      {
          addTriangle(mesh,vertices[0],vertices[1],vertices[2],nLinesTri);
      }

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

void InputMeshTri::meshtrifile_stl_binary(class TriMesh *mesh, class Region *region, const char *filename)
{
    unsigned int num_of_facets;
    std::ifstream stl_file;

    if (me == 0) {
        // open file for reading
        stl_file.open(filename, std::ifstream::in | std::ifstream::binary);

        // read 80 byte header into nirvana
        for (int i=0; i<20; i++){
          float dum;
          stl_file.read((char *)&dum, sizeof(float));
        }

        // read number of triangles
        stl_file.read((char *)&num_of_facets, sizeof(int));
    }

    // communicate number of triangles
    MPI_Bcast(&num_of_facets,1,MPI_INT,0,world);

    unsigned int count = 0;
    int nElems = 0;
    float tri_data[12];
    while(1) {
        if (me == 0) {
            // read one triangle dataset (normal + 3 vertices)
            for (int i=0; i<12; i++)
              stl_file.read((char *)&tri_data[i], sizeof(float));
            // read triangle attribute into nirvana
            short dum;
            stl_file.read((char *)&dum, sizeof(short));
            // error handling
            if (!stl_file)
                error->one(FLERR,"Corrupt STL file: Error in reading binary STL file.");
        }
        // communicate triangle data
        MPI_Bcast(&tri_data,12,MPI_FLOAT,0,world);
        // increase triangle counter
        count++;
        if(size_exclusion_list_ > 0 && count == (unsigned int)exclusion_list_[i_exclusion_list_]) {
            if(i_exclusion_list_ < size_exclusion_list_-1)
                i_exclusion_list_++;
        }
        else
        {
            double vert[9];
            // copy float vertices into double variables (ignore normal at the beginning of tri_data)
            for (int i=0; i<9; i++)
                vert[i] = (double) tri_data[i+3];
            if(!region || ( region->match(&vert[0]) && region->match(&vert[3]) && region->match(&vert[6]) ) )
            {
                addTriangle(mesh, &vert[0], &vert[3], &vert[6], count);
                nElems++;
                if (0 == (nElems % 100000) && comm->me == 0)
                    fprintf(screen,"   successfully read a chunk of 100000 elements in STL file '%s'\n",filename);
            }
        }

        if (count == num_of_facets)
            break;
    }

    if (me == 0)
        stl_file.close();
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
