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

    Stefan Amberger (JKU Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "region_wedge.h"
#include "error.h"
#include "domain.h"
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "math_extra.h"
#include "update.h"
#include "vector_liggghts.h"

// include last to ensure correct macros
#include "domain_definitions.h"

using namespace LAMMPS_NS;

/* -----------------------------------------------------------------------------
 OVERVIEW:
 constructor
 parses command line arguments
 checks validity conditions for input

 SETS THE FOLLOWING VARIABLES:
 extent_ylo and 'hi' for all directions (bounding box)
 -----------------------------------------------------------------------------*/

RegWedge::RegWedge(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  if(narg < 16)
    error->all(FLERR,"Illegal region wegde command, not enough arguments");

  options(narg-16, &arg[16]);

  // define helper-data
  pi_half = M_PI*0.5;

  int iarg = 2;

  if(strcmp(arg[iarg++],"axis"))
    error->all(FLERR,"Illegal region wegde command, expecting keyword 'axis'");

  // axis needs to be defined, otherwise: error
  axis = arg[iarg++][0];
  if (axis != 'x' && axis != 'y' && axis != 'z')
    error->all(FLERR,"Illegal region wedge command, expecting 'x', 'y', or 'z'");

  if (axis == 'x') {
    if(strcmp(arg[iarg++],"center"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'center'");
    c1 = yscale*atof(arg[iarg++]);
    c2 = zscale*atof(arg[iarg++]);
    if(strcmp(arg[iarg++],"radius"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'radius'");
    radius = yscale*atof(arg[iarg++]);
  } else if (axis == 'y') {
    if(strcmp(arg[iarg++],"center"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'center'");
    c1 = xscale*atof(arg[iarg++]);
    c2 = zscale*atof(arg[iarg++]);
    if(strcmp(arg[iarg++],"radius"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'radius'");
    radius = xscale*atof(arg[iarg++]);
  } else if (axis == 'z') {
    if(strcmp(arg[iarg++],"center"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'center'");
    c1 = xscale*atof(arg[iarg++]);
    c2 = yscale*atof(arg[iarg++]);
    if(strcmp(arg[iarg++],"radius"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'radius'");
    radius = xscale*atof(arg[iarg++]);
  }

  if(strcmp(arg[iarg++],"bounds"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'bounds'");

  // set 'lo' and 'hi' to be extent of cylinder
  // 'lo'
  if (strcmp(arg[iarg],"INF") == 0 || strcmp(arg[iarg],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[iarg++],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[0];
      else lo = domain->boxlo_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[iarg++],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[1];
      else lo = domain->boxlo_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[iarg++],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[2];
      else lo = domain->boxlo_bound[2];
    }
  } else {
    if (axis == 'x') lo = xscale*atof(arg[iarg++]);
    if (axis == 'y') lo = yscale*atof(arg[iarg++]);
    if (axis == 'z') lo = zscale*atof(arg[iarg++]);
  }

  // 'hi'
  if (strcmp(arg[iarg],"INF") == 0 || strcmp(arg[iarg],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[iarg++],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[0];
      else hi = domain->boxhi_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[iarg++],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[1];
      else hi = domain->boxhi_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[iarg++],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[2];
      else hi = domain->boxhi_bound[2];
    }
  } else {
    if (axis == 'x') hi = xscale*atof(arg[iarg++]);
    if (axis == 'y') hi = yscale*atof(arg[iarg++]);
    if (axis == 'z') hi = zscale*atof(arg[iarg++]);
    //debug
  }

  if(strcmp(arg[iarg++],"angle0"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'angle0'");
  angle1 = atof(arg[iarg++])*M_PI/180.0;
  cosang1 = cos(angle1);
  sinang1 = sin(angle1);

  if(strcmp(arg[iarg++],"angle"))
        error->all(FLERR,"Illegal region wegde command, expecting keyword 'angle'");
  dang = atof(arg[iarg++])*M_PI/180.0;
  cosdang = cos(dang);
  sindang = sin(dang);
  cosmdang = cos(-dang);
  sinmdang = sin(-dang);

  //debug

  if(iarg != 16)
    error->all(FLERR,"Illegal region wedge implementation, adapt options() call");
  //  fprintf(screen,"debug wedge: region id wedge %c %f %f %f %f %f %f %f\n",
  //       axis,c1,c2,radius,lo,hi,angle1,dang);

  // error checks
  if (radius <= 0.0) error->all(FLERR,"Illegal region cylinder command");
  if (dang < 5.0*M_PI/180.0) error->all(FLERR,"Wedge too flat. Wedge-angle has "
                                        "to be >= 5.0 degrees");
  if (dang > 180.0) error->all(FLERR, "Maximum wedge-angle of 180 "
                                          "deg exceeded");

  // calculate helper variables
  onedivr = 1.0/radius;
  angle2 = angle1 + dang;
  cosang2 = cos(angle2);
  sinang2 = sin(angle2);

  // extent of cylinder
  if (interior) {

    // bounding box is calculated as follows:
    // the bounding box is instantiated to encompass the origin and the two
    // bounding angles. then, in a loop from angle1 to angle2 all multiples
    // of pi_half are incorporated in the bounding box, in case the wedge does
    // not reside in one quadrant only.

    // a ... horizontal direction from viewpoint of axis (cos part of angle)
    // b ... vertical direction from viewpoint of axis (sin part of angle)

    double amin, amax;
    double bmin, bmax;

    amin = amax = 0.0;
    bmin = bmax = 0.0;

    // find min and max in both directions for end- and start-angle
    amin = MIN(amin, cosang1);
    amin = MIN(amin, cosang2);
    amax = MAX(amax, cosang1);
    amax = MAX(amax, cosang2);

    bmin = MIN(bmin, sinang1);
    bmin = MIN(bmin, sinang2);
    bmax = MAX(bmax, sinang1);
    bmax = MAX(bmax, sinang2);

    // sin and cos functions are only monotonous within the four quadrants.
    // if the wedge crossees quadrands, these functions might have a maximum
    // there, which should be a limit of the bounding box.

    double phi, sinphi, cosphi;
    int n = 0;
    // could be, that angle1 is -190 degrees
    while (static_cast<double>(n)*pi_half > angle1)
      n--;

    // since the wedge can have a maximum of 180 degrees
    // and we start smaller then angle1 i=0 to i=2 suffices
    for (int i = 0; i < 3; i++, n += pi_half){
      phi = static_cast<double>(n)*pi_half;

      if (angle1 <= phi && phi <= angle1 + dang){
      sinphi = sin(phi);
      cosphi = cos(phi);

      amin = MIN(amin, cosphi);
      amax = MAX(amax, cosphi);
      bmin = MIN(bmin, sinphi);
      bmax = MAX(bmax, sinphi);
      }
    }

    // adjust for radius
    amin *= radius;
    amax *= radius;
    bmin *= radius;
    bmax *= radius;

    bboxflag = 1;
    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 + amin;
      extent_yhi = c1 + amax;
      extent_zlo = c2 + bmin;
      extent_zhi = c2 + bmax;
    }
    if (axis == 'y') {
      extent_xlo = c1 + bmin;
      extent_xhi = c1 + bmax;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 + amin;
      extent_zhi = c2 + amax;
    }
    if (axis == 'z') {
      extent_xlo = c1 + amin;
      extent_xhi = c1 + amax;
      extent_ylo = c2 + bmin;
      extent_yhi = c2 + bmax;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else bboxflag = 0;

  // calculate normal vectors of the angular planes of the wedge
  double vec[2];
  // for angle1
  vec[0] = cosang1;
  vec[1] = sinang1;
  normal1[0] = vec[1];
  normal1[1] = -vec[0];
  // for angle2
  vec[0] = cosang2;
  vec[1] = sinang2;
  normal2[0] = -vec[1];
  normal2[1] = vec[0];

  // particle could contact cylinder surface and 2 ends and
  // 2 sides of the wedge
  cmax = 5;
  contact = new Contact[cmax];

}

/* -----------------------------------------------------------------------------
 destructor
 frees all memory allocated by this class
 -----------------------------------------------------------------------------*/

RegWedge::~RegWedge(){
  delete [] contact;
}

/* -----------------------------------------------------------------------------
 inside = 1 if x,y,z is inside or on surface
 inside = 0 if x,y,z is ouside and not on surface
 -----------------------------------------------------------------------------*/
int RegWedge::inside(double x, double y, double z){
  double lohi=0.0, distsq, sp1, sp2;
  double del[2]={};

  if (axis == 'x'){
    lohi = x;
    del[0] = y - c1;
    del[1] = z - c2;
  }
  else if(axis == 'y'){
    lohi = y;
    del[0] = z - c1;
    del[1] = x - c2;
  }
  else if(axis == 'z'){
    lohi = z;
    del[0] = x - c1;
    del[1] = y - c2;
  }

  // check if (x,y,z) between 'lo' and 'hi'
  if (lohi < lo || lohi > hi) return 0;

  // check radius
  distsq = del[0]*del[0] + del[1]*del[1];
  if (distsq > radius*radius) return 0;

  // check angle
  sp1 = vectorDot2D(normal1,del);
  sp2 = vectorDot2D(normal2,del);

  if (dang <= M_PI){
    if (sp1 > 0 || sp2 > 0) return 0;
  } else{
    if (sp1 > 0 && sp2 > 0) return 0;
  }

  return 1;
}

/* -----------------------------------------------------------------------------
 TODO
 -----------------------------------------------------------------------------*/
int RegWedge::surface_exterior(double *x, double cutoff){

  error->all(FLERR,"surface_exterior not implemented");

  return 0;
}

/* -----------------------------------------------------------------------------
 contact if 0 <= x < cutoff from one or more inner surfaces of cylinder
 can be one contact for each of 3 cylinder surfaces
 no contact if outside (possible if called from union/intersect)
 delxyz = vector from nearest point on cylinder to x
 special case: no contact with curved surf if x is on center axis
 ---------------------------------------------------------------------------- */
int RegWedge::surface_interior(double *x, double cutoff){
  int n = 0; // n ... number of contacts
  double contacts[5][3];  // general contacts - for actual contacts indices need
                          // to be changed depending on which axis the wedge
                          // is aligned to
  double delta;           // ... distance from surface to particle
  double lohi = 0.0;      // ... coord of particle in dimension of axis
  double del[2] = {};     // ... vector from center to particle
  double delxyz[2];       // ... vector from nearest point on surface to particle
  double rr;              // ... distance from point to center

  // assure that point is inside wedge
  if (!inside(x[0],x[1],x[2]))
      return n;

  // set parameters that depend upon what the axis is aligned to
  if (axis == 'x'){
    lohi = x[0];
    del[0] = x[1]-c1;
    del[1] = x[2]-c2;
  }
  else if (axis == 'y'){
    lohi = x[1];
    del[0] = x[2]-c1;
    del[1] = x[0]-c2;
  }
  else if (axis == 'z'){
    lohi = x[2];
    del[0] = x[0]-c1;
    del[1] = x[1]-c2;
  }
  rr = sqrt(del[0]*del[0]+del[1]*del[1]);

  // #1 upper wall
  delta = hi-lohi;
  if (delta < cutoff){
    contact[n].r = delta;
    contacts[n][0] = -delta;
    contacts[n][1] = 0.0;
    contacts[n][2] = 0.0;
    n++;
  }

  // #2 lower wall
  delta = lohi-lo;
  if (delta < cutoff){
    contact[n].r = delta;
    contacts[n][0] = delta;
    contacts[n][1] = 0.0;
    contacts[n][2] = 0.0;
    n++;
  }

  // #3 outer cylinder wall
  delta = radius - rr;
  if (delta < cutoff && rr > 0.0) {
    contact[n].r = delta;
    contacts[n][0] = 0.0;

    // delxyz is vector from center to particle w/ length radius
    snormalize2(radius, del, delxyz);
    delxyz[0] -= del[0];
    delxyz[1] -= del[1];

    // contacts is vector from nearest point on surface to particle
    contacts[n][1] = -delxyz[0];
    contacts[n][2] = -delxyz[1];
    n++;
  }

  // TODO what if dang > 180 --> wrong
  // #4 first plane of wedge
  delta = vectorDot2D(del,normal1);
  if (fabs(delta) < cutoff){
    contact[n].r = fabs(delta);
    contacts[n][0] = 0.0;
    contacts[n][1] = delta*normal1[0];
    contacts[n][2] = delta*normal1[1];
    n++;
  }

  // TODO what if dang > 180 --> wrong
  // #5 second plane of wedge
  delta = vectorDot2D(del,normal2);
  if (fabs(delta) < cutoff){
    contact[n].r = fabs(delta);
    contacts[n][0] = 0.0;
    contacts[n][1] = delta*normal2[0];
    contacts[n][2] = delta*normal2[1];
    n++;
  }

  // change indices according to axis in use
  if (axis == 'x'){
    for (int i = 0; i<n; i++){
      contact[i].delx = contacts[i][0];
      contact[i].dely = contacts[i][1];
      contact[i].delz = contacts[i][2];
    }
  }
  else if (axis == 'y') {
    for (int i = 0; i<n; i++){
      contact[i].delx = contacts[i][2];
      contact[i].dely = contacts[i][0];
      contact[i].delz = contacts[i][1];
    }
  }
  else if (axis == 'z') {
    for (int i = 0; i<n; i++){
      contact[i].delx = contacts[i][1];
      contact[i].dely = contacts[i][2];
      contact[i].delz = contacts[i][0];
    }
  }
  else {
    error->all(FLERR, "axis does not match 'x', 'y' or 'z'");
  }

  // debug
  // printContacts(x,n);

  return n;
}

/* -----------------------------------------------------------------------------
 helper functions
 -----------------------------------------------------------------------------*/

void RegWedge::snormalize2(const double length, const double *v, double *ans){
  double scale = length/sqrt(v[0]*v[0] + v[1]*v[1]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
}

void RegWedge::printRegion(){
  printProperty("interior", interior);
  printProperty("c1", c1);
  printProperty("c2", c2);
  printProperty("lo", lo);
  printProperty("hi", hi);
  printProperty("dang", dang);
  printProperty("radius", radius);
  printProperty("inside(1,1,1)", inside(1.0, 1.0, 1.0));
}

void RegWedge::printProperty(const char *name, double val){
  printf("%s:    %f\n",name,val);
}

void RegWedge::printContacts(double *x, int n){
  for (int i = 0; i<n; i++){
    printf("step " BIGINT_FORMAT " Contact %i\n",update->ntimestep,i);
    printf("\tx\t: %f\t%f\t%f\n",x[0],x[1],x[2]);
    printf("\tr\t: %f\n\tdx\t: %f\n\tdy\t: %f\n\tdz\t: %f\n",contact[i].r,contact[i].delx,contact[i].dely,contact[i].delz);
  }
  if (n>0)
    printf("\n");
}
