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

#ifndef LMP_MESH_MOVER_H
#define LMP_MESH_MOVER_H

#include "tri_mesh.h"
#include "force.h"

namespace LAMMPS_NS
{
  class MeshMover : protected Pointers
  {
      public:
        MeshMover(LAMMPS * lmp,AbstractMesh *_mesh) :
            Pointers(lmp), mesh_(_mesh), isFirst_(false) {}
        virtual ~MeshMover() {}

        virtual void initial_integrate(double dT,double dt) = 0;
        virtual void final_integrate(double dT,double dt) {};
        virtual void pre_delete() = 0;
        inline bool isFirst()
        { return isFirst_; }

        double ***get_nodes()
        {
            return mesh_->nodePtr();
        }

        double ***get_v()
        {
            double ***ptr;
            if(mesh_->numNodes() == 3)
                ptr = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();
            else if(mesh_->numNodes() == 4)
                return mesh_->prop().getElementProperty<MultiVectorContainer<double,4,3> >("v")->begin();
            if(!ptr)
                error->one(FLERR,"Illegal call to MeshMover::get_v");
            return ptr;
        }

        AbstractMesh *mesh_;
        bool isFirst_;
  };

  MeshMover* createMeshMover(LAMMPS *lmp,AbstractMesh *mesh,char **arg, int narg);

  /* ---------------------------------------------------------------------- */

  class MeshMoverLinear : public MeshMover{

      public:
        MeshMoverLinear(LAMMPS *lmp,AbstractMesh *_mesh,double vx, double vy, double vz);
        virtual ~MeshMoverLinear();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        double vel_[3];
  };

  /* ----------------------------------------------------------------------- */

   class MeshMoverLinearVariable : public MeshMover{

      public:
        MeshMoverLinearVariable(LAMMPS *lmp,AbstractMesh *_mesh,char* var1, char* var2, char* var3);
        virtual ~MeshMoverLinearVariable();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        char *var1str,*var2str,*var3str;
        int myvar1,myvar2,myvar3;
        double dX_[3];
        double vel_[3];
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverWiggle : public MeshMover{

      public:
        MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh,
                        double ax, double ay, double az,
                        double T);
        virtual ~MeshMoverWiggle();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        double amplitude[3],omega;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRotate : public MeshMover {

      public:
        MeshMoverRotate(LAMMPS *lmp,AbstractMesh *_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            double T);
        virtual ~MeshMoverRotate();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        double axis[3], p[3], omega;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRotateVariable : public MeshMover {

      public:
        MeshMoverRotateVariable(LAMMPS *lmp,AbstractMesh *_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            char* var1);
        virtual ~MeshMoverRotateVariable();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        char *var1str;
        int myvar1,myvar2,myvar3;
        double axis[3], p[3], omega, totalPhi;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRiggle : public MeshMover {

      public:
        MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            double T, double ampl);
        virtual ~MeshMoverRiggle();

        void initial_integrate(double dT,double dt);
        void final_integrate(double dT,double dt) {}
        void pre_delete();

      private:
        double axis[3], p[3], omega, amplitude;
  };

   /* ----------------------------------------------------------------------
    Mesh Mover
    ------------------------------------------------------------------------- */

    inline MeshMover* createMeshMover(LAMMPS *lmp,AbstractMesh *mesh,char **arg, int narg)
    {
        if(narg < 1) return 0;

        char *name = arg[0];
        if(strcmp(name,"linear") == 0){
          if(narg < 4) return 0;

          return new MeshMoverLinear(lmp,mesh,
                          lmp->force->numeric(arg[1]),
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]));
        } else if(strcmp(name,"linear/variable") == 0){
          if(narg < 4) return 0;

          return new MeshMoverLinearVariable(lmp,mesh,
                          arg[1],
                          arg[2],
                          arg[3]);
        } else if(strcmp(name,"rotate") == 0){
          if(narg < 11) return 0;
          else
          {
            if(strcmp("origin",arg[1]))
                return 0;
            if(strcmp("axis",arg[5]))
                return 0;
            if(strcmp("period",arg[9]))
                return 0;

            return new MeshMoverRotate(lmp,mesh,
                          // origin
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          // axis
                          lmp->force->numeric(arg[6]),
                          lmp->force->numeric(arg[7]),
                          lmp->force->numeric(arg[8]),
                          // period
                          lmp->force->numeric(arg[10]));
          }
	    } else if(strcmp(name,"rotate/variable") == 0){
          if(narg < 11) return 0;
          else
          {
            if(strcmp("origin",arg[1]))
                return 0;
            if(strcmp("axis",arg[5]))
                return 0;
            if(strcmp("period",arg[9]))
                return 0;

            return new MeshMoverRotateVariable(lmp,mesh,
                          // origin
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          // axis
                          lmp->force->numeric(arg[6]),
                          lmp->force->numeric(arg[7]),
                          lmp->force->numeric(arg[8]),
                          // variable name for OMEGA (because var could be zero !)
                          arg[10]);
          }
        } else if(strcmp(name,"wiggle") == 0){
          if(narg < 7) return 0;
          else
          {
            if(strcmp("amplitude",arg[1]))
                return 0;
            if(strcmp("period",arg[5]))
                return 0;

            return new MeshMoverWiggle(lmp,mesh,
                          //amplitude
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          //period
                          lmp->force->numeric(arg[6]));
          }
        } else if(strcmp(name,"riggle") == 0){
          if(narg < 13) return 0;
          else
          {
            if(strcmp("origin",arg[1]))
                return 0;
            if(strcmp("axis",arg[5]))
                return 0;
            if(strcmp("period",arg[9]))
                return 0;
            if(strcmp("amplitude",arg[11]))
                return 0;

            return new MeshMoverRiggle(lmp,mesh,
                          // origin
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          // axis
                          lmp->force->numeric(arg[6]),
                          lmp->force->numeric(arg[7]),
                          lmp->force->numeric(arg[8]),
                          // period
                          lmp->force->numeric(arg[10]),
                          // amplitude
                          lmp->force->numeric(arg[12]));
          }
        }

        return 0;
    }

} /* LAMMPS_NS */

#endif /* MESHMOVER_H_ */
