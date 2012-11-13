/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_PRIMITIVE_WALL
#define LMP_PRIMITIVE_WALL

#include "container.h"
#include "primitive_wall_definitions.h"

namespace LAMMPS_NS
{

  /*
   * class PrimitiveWall only holds logic for resolving contacts, neighbor list and contact tracking
   */

  class PrimitiveWall : protected Pointers
  {
      public:
        PrimitiveWall(LAMMPS *lmp,PRIMITIVE_WALL_DEFINITIONS::WallType wType_, int nParam_, double *param_);
        virtual
        ~PrimitiveWall();

        int getNeighbors(int *&contactPtr);
        void handleContact(int iPart,double *&c_history);
        void handleNoContact(int iPart);
        void setContactHistorySize(int nPart);

        void buildNeighList(double neighCutoff, double **x, double *r, int nPart);

        double resolveContact(double *x, double r, double *delta);
        bool resolveNeighlist(double *x, double r, double treshold);

        int axis();
        double calcRadialDistance(double *pos, double *distvec);

      private:
        ScalarContainer<int> neighlist;
        ScalarContainer<int> contactHashmap;
        VectorContainer<double,3> history;
        PRIMITIVE_WALL_DEFINITIONS::WallType wType;

        double *param;
        int nParam;

  };

  /*
   * implementation of class PrimitiveWall starts here
   */

  PrimitiveWall::PrimitiveWall(LAMMPS *lmp,PRIMITIVE_WALL_DEFINITIONS::WallType wType_, int nParam_, double *param_)
  : Pointers(lmp), wType(wType_), nParam(nParam_)
  {
    param = new double[nParam];
    for(int i=0;i<nParam;i++)
      param[i] = param_[i];
  }

  PrimitiveWall::~PrimitiveWall()
  {
    delete []param;
  }

  double PrimitiveWall::resolveContact(double *x, double r, double *delta)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseContactTemplate(x, r, delta, param, wType);
  }

  int PrimitiveWall::axis()
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseAxis(wType);
  }

  double PrimitiveWall::calcRadialDistance(double *pos, double *distvec)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseCalcRadialDistance(pos, param, distvec[0],distvec[1],distvec[2], wType);
  }

  bool PrimitiveWall::resolveNeighlist(double *x, double r, double treshold)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseNeighlistTemplate(x,r,treshold,param,wType);
  }

  int PrimitiveWall::getNeighbors(int *&contactPtr)
  {
    contactPtr = neighlist.begin();
    return neighlist.size();
  }

  void PrimitiveWall::handleContact(int iPart, double *&c_history)
  {
    int index = contactHashmap(iPart);
    static double zeroVec[] = {0.,0.,0.};
    if(index != -1){
      c_history = history(index);
      //printf("contact to particle %d found\n",iPart);
    } else{
      contactHashmap(iPart) = history.size();
      history.add(zeroVec);
      c_history = history(history.size()-1);
      //printf("contact to particle %d added\n",iPart);
    }
  }

  void PrimitiveWall::handleNoContact(int iPart)
  {
    int index = contactHashmap(iPart);
    //printf("history-size: %d | index: %d ",history.size(),index);
    if(index == -1) return;
    //printf("contact to particle %d deleted\n",iPart);
    //printf("*** history-size: %d | index: %d \n",history.size(),index);
      contactHashmap(iPart) = -1;
      history.del(index);
  }

  void PrimitiveWall::setContactHistorySize(int nPart)
  {
    int currentSize = contactHashmap.size();
    if(nPart > currentSize){
      int delta = nPart - currentSize;
      contactHashmap.addUninitialized(delta);
      for(int i=0;i<delta;i++)
        contactHashmap(currentSize+i) = -1;
    }

  }
  void PrimitiveWall::buildNeighList(double treshold, double **x, double *r, int nPart)
  {
    neighlist.empty();
    for(int iPart=0;iPart<nPart;iPart++)
    {
      
      if(resolveNeighlist(x[iPart],r?r[iPart]:0.,treshold))
        neighlist.add(iPart);
    }
  }
} /* namespace LAMMPS_NS */
#endif /* PRIMITIVEWALL_H_ */
