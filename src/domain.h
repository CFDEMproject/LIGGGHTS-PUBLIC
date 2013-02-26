/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include "pointers.h"
#include "error.h" 
#include "comm.h" 
#include "vector_liggghts.h" 

#define SMALL_DMBRDR 1.0e-8 

namespace LAMMPS_NS {

class Domain : protected Pointers {
 public:
  int box_exist;                         // 0 = not yet created, 1 = exists
  int dimension;                         // 2 = 2d, 3 = 3d
  int nonperiodic;                       // 0 = periodic in all 3 dims
                                         // 1 = periodic or fixed in all 6
                                         // 2 = shrink-wrap in any of 6
  int xperiodic,yperiodic,zperiodic;     // 0 = non-periodic, 1 = periodic
  int periodicity[3];                    // xyz periodicity as array

  int boundary[3][2];                    // settings for 6 boundaries
                                         // 0 = periodic
                                         // 1 = fixed non-periodic
                                         // 2 = shrink-wrap non-periodic
                                         // 3 = shrink-wrap non-per w/ min

  int triclinic;                         // 0 = orthog box, 1 = triclinic

                                         // orthogonal box
  double xprd,yprd,zprd;                 // global box dimensions
  double xprd_half,yprd_half,zprd_half;  // half dimensions
  double prd[3];                         // array form of dimensions
  double prd_half[3];                    // array form of half dimensions

                                         // triclinic box
                                         // xprd,xprd_half,prd,prd_half =
                                         // same as if untilted
  double prd_lamda[3];                   // lamda box = (1,1,1)
  double prd_half_lamda[3];              // lamda half box = (0.5,0.5,0.5)

  double boxlo[3],boxhi[3];              // orthogonal box global bounds

                                         // triclinic box
                                         // boxlo/hi = same as if untilted
  double boxlo_lamda[3],boxhi_lamda[3];  // lamda box = (0,1)
  double boxlo_bound[3],boxhi_bound[3];  // bounding box of tilted domain
  double corners[8][3];                  // 8 corner points

                                         // orthogonal box & triclinic box
  double minxlo,minxhi;                  // minimum size of global box
  double minylo,minyhi;                  // when shrink-wrapping
  double minzlo,minzhi;                  // tri only possible for non-skew dims

                                         // orthogonal box
  double sublo[3],subhi[3];              // sub-box bounds on this proc

                                         // triclinic box
                                         // sublo/hi = undefined
  double sublo_lamda[3],subhi_lamda[3];  // bounds of subbox in lamda

                                         // triclinic box
  double xy,xz,yz;                       // 3 tilt factors
  double h[6],h_inv[6];                           // shape matrix in Voigt notation
  double h_rate[6],h_ratelo[3];          // rate of box size/shape change

  int box_change;                 // 1 if box bounds ever change, 0 if fixed
  int deform_flag;                // 1 if fix deform exist, else 0
  int deform_vremap;              // 1 if fix deform remaps v, else 0
  int deform_groupbit;            // atom group to perform v remap for

  class Lattice *lattice;                  // user-defined lattice

  int nregion;                             // # of defined Regions
  int maxregion;                           // max # list can hold
  class Region **regions;                  // list of defined Regions

  Domain(class LAMMPS *);
  virtual ~Domain();
  virtual void init();
  void set_initial_box();
  virtual void set_global_box();
  virtual void set_lamda_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  int minimum_image_check(double, double, double);
  void minimum_image(double &, double &, double &);
  void minimum_image(double *);
  void closest_image(const double * const, const double * const,
                     double * const);
  void remap(double *, int &);
  void remap(double *);
  void remap_near(double *, double *);
  void unmap(double *, int);
  void unmap(double *, int, double *);
  void image_flip(int, int, int);
  void set_lattice(int, char **);
  void add_region(int, char **);
  void delete_region(int, char **);
  int find_region(char *);
  void set_boundary(int, char **, int);
  void print_box(const char *);

  virtual void lamda2x(int);
  virtual void x2lamda(int);
  void lamda2x(double *, double *);
  void x2lamda(double *, double *);
  void x2lamda(double *, double *, double *, double *);
  void bbox(double *, double *, double *, double *);
  void box_corners();

  inline int is_in_domain(double* pos); 
  inline int is_in_subdomain(double* pos); 
  inline int is_in_extended_subdomain(double* pos); 
  inline double dist_subbox_borders(double* pos); 
  int is_periodic_ghost(int i); 
  int is_periodic_ghost_of_owned(int i); 

 private:
  double small[3];                  // fractions of box lengths
};

/* ----------------------------------------------------------------------
   check if coordinate in domain, subdomain or extended subdomain
   need to test with <= and >= for domain, and < and >= for subdomain
   inlined for performance
------------------------------------------------------------------------- */

inline int Domain::is_in_domain(double* pos) 
{
    if
    (
        pos[0] >= boxlo[0] && pos[0] <= boxhi[0] &&
        pos[1] >= boxlo[1] && pos[1] <= boxhi[1] &&
        pos[2] >= boxlo[2] && pos[2] <= boxhi[2]
    )   return 1;
    return 0;
}

inline int Domain::is_in_subdomain(double* pos) 
{
    double checkhi[3];

    // need to test with and < and >= for subdomain
    // but need to ensure elements right at boxhi
    // are tested correctly
    if(subhi[0] == boxhi[0])
        checkhi[0] = boxhi[0] + SMALL_DMBRDR;
    else
        checkhi[0] = subhi[0];

    if(subhi[1] == boxhi[1])
        checkhi[1] = boxhi[1] + SMALL_DMBRDR;
    else
        checkhi[1] = subhi[1];

    if(subhi[2] == boxhi[2])
        checkhi[2] = boxhi[2] + SMALL_DMBRDR;
    else
        checkhi[2] = subhi[2];

    if ( pos[0] >= sublo[0] && pos[0] < checkhi[0] &&
         pos[1] >= sublo[1] && pos[1] < checkhi[1] &&
         pos[2] >= sublo[2] && pos[2] < checkhi[2])
        return 1;
    return 0;
}

inline int Domain::is_in_extended_subdomain(double* pos) 
{
    // called on insertion
    // yields true if particle would be in subdomain after box extension
    
    if (is_in_subdomain(pos))
        return 1;
    else if (dimension == 2)
        error->all(FLERR,"Domain::is_in_extended_subdomain() not implemented for 2d");
    else 
    {
        bool flag = true;
        for(int idim = 0; idim < 3; idim++)
        {
            
            if (comm->procgrid[idim] == 1) {}
            else if(comm->myloc[idim] == comm->procgrid[idim]-1)
                flag = flag && (pos[idim] >= sublo[idim]);
            else if(comm->myloc[idim] == 0)
                flag = flag && (pos[idim] <= subhi[idim]);
            
            else
                flag = flag && (pos[idim] >= sublo[idim] && pos[idim] < subhi[idim]);
        }
        if(flag) return 1;
        return 0;
    }
    return 0;
}

inline double Domain::dist_subbox_borders(double* pos) 
{
    double deltalo[3], deltahi[3], dlo, dhi;

    vectorSubtract3D(sublo,pos,deltalo);
    vectorSubtract3D(subhi,pos,deltahi);
    vectorAbs3D(deltalo);
    vectorAbs3D(deltahi);

    dlo = vectorMin3D(deltalo);
    dhi = vectorMin3D(deltahi);

    if(dlo < dhi)
        return dlo;
    return dhi;
}

}

#endif

/* ERROR/WARNING messages:

E: Box bounds are invalid

The box boundaries specified in the read_data file are invalid.  The
lo value must be less than the hi value for all 3 dimensions.

E: Cannot skew triclinic box in z for 2d simulation

Self-explanatory.

E: Triclinic box skew is too large

The displacement in a skewed direction must be less than half the box
length in that dimension.  E.g. the xy tilt must be between -half and
+half of the x box length.

E: Illegal simulation box

The lower bound of the simulation box is greater than the upper bound.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Reuse of region ID

A region ID cannot be used twice.

E: Invalid region style

The choice of region style is unknown.

E: Delete region ID does not exist

Self-explanatory.

E: Both sides of boundary must be periodic

Cannot specify a boundary as periodic only on the lo or hi side.  Must
be periodic on both sides.

*/
