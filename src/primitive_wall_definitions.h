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

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_PRIMITIVE_WALL_DEFINITIONS
#define LMP_PRIMITIVE_WALL_DEFINITIONS

/*
 * Necessary steps to add new primitive walls:
 * (1) add an enum for your primitive to WallType, but insert it before NUM_WTYPE
 * (2) add a string that you want to use in your input script to wallString and
 *     the number of arguments the wall requires to numArgs
 * (3) implement distance and neighbor list build functions
 * (4) add them to the switch statements in chooseContactTemplate() and
 *     chooseNeighlistTemplate() located at the bottom of this file
 */

namespace LAMMPS_NS
{
  namespace PRIMITIVE_WALL_DEFINITIONS
  {

    enum WallType
    {
        XPLANE,
        YPLANE,
        ZPLANE,
        XCYLINDER,
        YCYLINDER,
        ZCYLINDER,
        NUM_WTYPE
    };

    static const char *wallString[] =
    {
        "xplane",
        "yplane",
        "zplane",
        "xcylinder",
        "ycylinder",
        "zcylinder"
    };

    static int numArgs[] =
    {
        1,
        1,
        1,
        3,
        3,
        3
    };

    /*
     * templates for the different wall primitives
     * each wall type needs to be coded in a template specialization
     * default implementation: no contact
     */
    template<WallType W>
    double resolveContactTemplate(double *x, double r, double *delta, double *param) {return 1.;}

    /*
     * default neighbor list template returns true --> if neighbor list distance function is not
     * coded explicitly, every particle will be added
     */
    template<WallType W>
    bool resolveNeighlistTemplate(double *x, double r, double treshold, double *param) {return true;}

    /*
     * declaration of choosing functions
     * definitions need to be after ALL definitions of resolveContactTemplate and resolveNeighlistTemplate
     */
    inline double chooseContactTemplate(double *x, double r, double *delta, double *param, WallType wType);
    inline bool chooseNeighlistTemplate(double *x, double r, double treshold, double *param, WallType wType);

/* ---------------------------------------------------------------------- */

    /*
     * x,y,z planes can be handled by a template with dimension as template parameter
     */

    template<int dim>
    struct Dim
    {
      static const int x = dim, y = (dim+1)%3, z = (dim+2)%3;
    };

    template<int dim>
    struct Plane : public Dim<dim>
    {
      typedef Dim<dim> d;
      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        if(pos[d::x] > *param){
          delta[d::x] = *param - pos[d::x]; delta[d::y] = 0.; delta[d::z] = 0.;
          return pos[d::x] - *param - r;
        } else{
          delta[d::x] = *param - pos[d::x]; delta[d::y] = 0.; delta[d::z] = 0.;
          return *param - pos[d::x] - r;
        }
      }
      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dMax = r + treshold;
        double dist = pos[d::x] - *param;
        double absdist = (dist > 0.0) ? dist : -dist;
        return (absdist <= dMax);
      }
    };

/* ---------------------------------------------------------------------- */
    /*
     * same holds for x,y,z cylinders
     * param[0] = radius
     * param[1] = first coordinate of center
     * param[2] = second coordinate of center
     */
    template<int dim>
    struct Cylinder : public Dim<dim>
    {
    public:

      typedef Dim<dim> d;
      static double calcRadialDistance(double *pos, double *param, double &dy, double &dz)
      {
        dy = pos[d::y]-param[1];
        dz = pos[d::z]-param[2];
        return sqrt(dy*dy+dz*dz);
      }

      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        double dx, dy,dz, fact;
        double dist = calcRadialDistance(pos,param,dy,dz);
        if(dist > *param){
          dx = dist - *param - r;
          fact = (dist - *param) / dist;
          delta[d::x] = 0.; delta[d::y] = -dy*fact; delta[d::z] = -dz*fact;

        } else{
          dx = *param - dist - r;
          fact = (*param - dist) / dist;
          delta[d::x] = 0.; delta[d::y] = +dy*fact; delta[d::z] = +dz*fact;
        }
        return dx;
      }
      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dy,dz;
        double dMax = r + treshold;
        double dist = calcRadialDistance(pos,param,dy,dz) - *param;
        return (dMax < dist || -dMax < dist);
      }

    };

/* ---------------------------------------------------------------------- */

    /*
     * functions to choose the correct template
     */

/* ---------------------------------------------------------------------- */

    inline double chooseContactTemplate(double *x, double r, double *delta, double *param, WallType wType)
    {
      //TODO: find a way to create switch statement automatically
      switch(wType){
      case XPLANE:
        return Plane<0>::resolveContact(x,r,delta,param);
      case YPLANE:
        return Plane<1>::resolveContact(x,r,delta,param);
      case ZPLANE:
        return Plane<2>::resolveContact(x,r,delta,param);
      case XCYLINDER:
        return Cylinder<0>::resolveContact(x,r,delta,param);
      case YCYLINDER:
        return Cylinder<1>::resolveContact(x,r,delta,param);
      case ZCYLINDER:
        return Cylinder<2>::resolveContact(x,r,delta,param);

      default: // default: no contact
        return 1.;
      }
    }

    inline bool chooseNeighlistTemplate(double *x, double r, double treshold, double *param, WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XPLANE:
        return Plane<0>::resolveNeighlist(x,r,treshold,param);
      case YPLANE:
        return Plane<1>::resolveNeighlist(x,r,treshold,param);
      case ZPLANE:
        return Plane<2>::resolveNeighlist(x,r,treshold,param);
      case XCYLINDER:
        return Cylinder<0>::resolveNeighlist(x,r,treshold,param);
      case YCYLINDER:
        return Cylinder<1>::resolveNeighlist(x,r,treshold,param);
      case ZCYLINDER:
        return Cylinder<2>::resolveNeighlist(x,r,treshold,param);

      default: // default value: every particle will be added to neighbor list
        return true;
      }
    }

    inline int chooseAxis(WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XCYLINDER:
        return 0;
      case YCYLINDER:
        return 1;
      case ZCYLINDER:
        return 2;

      default:
        return -1;
      }
    }

    inline double chooseCalcRadialDistance(double *pos, double *param, double &dx, double &dy, double &dz, WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XCYLINDER:
        dx=0.;
        return Cylinder<0>::calcRadialDistance(pos,param,dy,dz);
      case YCYLINDER:
        dy = 0.;
        return Cylinder<1>::calcRadialDistance(pos,param,dx,dz);
      case ZCYLINDER:
        dz = 0.;
        return Cylinder<2>::calcRadialDistance(pos,param,dx,dy);

      default:
        dx = dy = dz = 0.;
        return -1.;
      }
    }

    inline int chooseNumArgs(char *style)
    {
        for(int i = 0; i < NUM_WTYPE; i++)
            if(strcmp(style,wallString[i]) == 0)
                return numArgs[i];
        return 0;
    }

    inline int numArgsPrimitiveWall(char *style)
    {
        return chooseNumArgs(style);
    }
  }
}

#endif /* PRIMITIVEWALLDEFINITIONS_H_ */
