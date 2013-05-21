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
   Niels Dallinger (TU Chemnitz, viblin and vibrot)
------------------------------------------------------------------------- */

#ifndef LMP_MESH_MOVER_H
#define LMP_MESH_MOVER_H

#include "tri_mesh.h"
#include "fix_move_mesh.h"
#include "force.h"

namespace LAMMPS_NS
{
  class MeshMover : protected Pointers
  {
      public:

        MeshMover(LAMMPS * lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh) :
            Pointers(lmp),
            mesh_(_mesh),
            fix_move_mesh_(_fix_move_mesh),
            isFirst_(false)
        {}

        virtual ~MeshMover()
        {}

        virtual void pre_delete() = 0;
        virtual void setup() {};

        virtual void initial_integrate(double dTAbs,double dTSetup,double dt) = 0;
        virtual void final_integrate(double dTAbs,double dTSetup,double dt) {};

        inline bool isFirst()
        { return isFirst_; }

        virtual int n_restart()
        { return 0; }

        virtual void write_restart(double *buf) {}
        virtual void read_restart(double *buf) {}

        void add_reference_point(double *point)
        { fix_move_mesh_->add_reference_point(point); }

        void get_reference_point(double *point)
        { fix_move_mesh_->get_reference_point(point); }

        double ***get_nodes()
        { return mesh_->nodePtr(); }

        double ***get_v()
        {
            double ***ptr = NULL;
            if(mesh_->numNodes() == 3)
                ptr = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin();
            else if(mesh_->numNodes() == 4)
                return mesh_->prop().getElementProperty<MultiVectorContainer<double,4,3> >("v")->begin();
            if(!ptr)
                error->one(FLERR,"Illegal call to MeshMover::get_v");
            return ptr;
        }

        AbstractMesh *mesh_;
        FixMoveMesh *fix_move_mesh_;
        bool isFirst_;
  };

  MeshMover* createMeshMover(LAMMPS *lmp,AbstractMesh *mesh,char **arg, int narg);

  /* ---------------------------------------------------------------------- */

  class MeshMoverLinear : public MeshMover{

      public:

        MeshMoverLinear(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                        double vx, double vy, double vz);
        virtual ~MeshMoverLinear();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:

        double vel_[3];
  };

  /* ----------------------------------------------------------------------- */

   class MeshMoverLinearVariable : public MeshMover{

      public:

        MeshMoverLinearVariable(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                                char* var1, char* var2, char* var3);
        virtual ~MeshMoverLinearVariable();

        void pre_delete();
        void setup();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}

        int n_restart();
        void write_restart(double *buf);
        void read_restart(double *buf);

      private:

        char *var1str_,*var2str_,*var3str_;
        int myvar1_,myvar2_,myvar3_;
        double dX_[3];
        double vel_[3];
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverWiggle : public MeshMover{

      public:

        MeshMoverWiggle(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                        double ax, double ay, double az,
                        double T);
        virtual ~MeshMoverWiggle();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:

        double amplitude_[3],omega_;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRotate : public MeshMover {

      public:

        MeshMoverRotate(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            double T);
        virtual ~MeshMoverRotate();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:

        double axis_[3], point_[3], omega_;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRotateVariable : public MeshMover {

      public:

        MeshMoverRotateVariable(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            char* var1);
        virtual ~MeshMoverRotateVariable();

        void pre_delete();
        void setup();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}

        int n_restart();
        void write_restart(double *buf);
        void read_restart(double *buf);

      private:

        char *var1str_;
        int myvar1_,myvar2_,myvar3_;
        double axis_[3], point_[3], omega_, totalPhi_;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverRiggle : public MeshMover {

      public:

        MeshMoverRiggle(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            double T, double ampl);
        virtual ~MeshMoverRiggle();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:

        double axis_[3], point_[3], omega_, amplitude_;
  };

  /* ---------------------------------------------------------------------- */

  class MeshMoverVibRot : public MeshMover {

      public:
        MeshMoverVibRot(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            double px, double py,double pz,
                            double axisX, double axisY, double axisZ,
                            int order, double amplitude[10], double phase[10],
                            double T);
        virtual ~MeshMoverVibRot();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:
        double axis_[3], ord, ampl[10], phi[10], p_[3], omega_;

  };

 /* ---------------------------------------------------------------------- */
  class MeshMoverVibLin : public MeshMover {

      public:
        MeshMoverVibLin(LAMMPS *lmp,AbstractMesh *_mesh,FixMoveMesh *_fix_move_mesh,
                            double axisX, double axisY, double axisZ,
                            int order, double amplitude[10], double phase[10],
                            double T);
        virtual ~MeshMoverVibLin();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}
        void pre_delete();

      private:
        double axis_[3], ord, omega_, ampl[10], phi[10];

  };

   /* ----------------------------------------------------------------------
    Mesh Mover
    ------------------------------------------------------------------------- */

    inline MeshMover* createMeshMover(LAMMPS *lmp,AbstractMesh *mesh,FixMoveMesh *fix_mm,char **arg, int narg)
    {
        if(narg < 1) return 0;

        char *name = arg[0];
        if(strcmp(name,"linear") == 0){
          if(narg < 4) return 0;

          return new MeshMoverLinear(lmp,mesh,fix_mm,
                          lmp->force->numeric(arg[1]),
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]));
        } else if(strcmp(name,"linear/variable") == 0){
          if(narg < 4) return 0;

          return new MeshMoverLinearVariable(lmp,mesh,fix_mm,
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

            return new MeshMoverRotate(lmp,mesh,fix_mm,
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
            if(strcmp("omega",arg[9]))
                return 0;

            return new MeshMoverRotateVariable(lmp,mesh,fix_mm,
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

            return new MeshMoverWiggle(lmp,mesh,fix_mm,
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

            return new MeshMoverRiggle(lmp,mesh,fix_mm,
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
         else if(strcmp(name,"viblin") == 0){
            int order = lmp->force->numeric(arg[6]);
            if(narg < 10+2*order) return 0;
          else
          {
            if(strcmp("axis",arg[1]))
                return 0;
            if(strcmp("order",arg[5]))
                return 0;
            if(strcmp("amplitude",arg[7]))
                return 0;
            if(strcmp("phase",arg[8+order]))
                return 0;
            if(strcmp("period",arg[9+2*order]))
                return 0;
            double pha[10];
            double amp[10];
            // creating array of amplitude and phase
            for (int zv=0; zv<order; zv++) {
                //amplitude
                amp[zv] = lmp->force->numeric(arg[8+zv]);
                // angle of phase
                pha[zv] = lmp->force->numeric(arg[9+order+zv]);
               }

            return new MeshMoverVibLin(lmp,mesh,fix_mm,
                          // direction
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          // order
                          lmp->force->numeric(arg[6]),
                          // amplitudes
                          amp,
                          // phases
                          pha,
                          // periode
                          lmp->force->numeric(arg[10+2*order]));
          }
        }
        else if(strcmp(name,"vibrot") == 0){
             int order = lmp->force->numeric(arg[10]);
             if (narg < 14+2*order) return 0;
          else
          {
            if(strcmp("origin",arg[1]))
                return 0;
            if(strcmp("axis",arg[5]))
                return 0;
            if(strcmp("order",arg[9]))
                return 0;
            if(strcmp("amplitude",arg[11]))
                return 0;
            if(strcmp("phase",arg[12+order]))
                return 0;
            if(strcmp("period",arg[13+2*order]))
                return 0;
            double pha[10];
            double amp[10];
            // creating array of amplitude and phase
            for (int zv=0; zv<order; zv++) {
                //amplitude
                amp[zv] = lmp->force->numeric(arg[12+zv]);
                // angle of phase
                pha[zv] = lmp->force->numeric(arg[13+order+zv]);
               }
            return new MeshMoverVibRot(lmp,mesh,fix_mm,
                          // origin px py pz
                          lmp->force->numeric(arg[2]),
                          lmp->force->numeric(arg[3]),
                          lmp->force->numeric(arg[4]),
                          // axis axisx axisy axisz
                          lmp->force->numeric(arg[6]),
                          lmp->force->numeric(arg[7]),
                          lmp->force->numeric(arg[8]),
                          // order
                          lmp->force->numeric(arg[10]),
                          // amplitudes
                          amp,
                          // phases
                          pha,
                          // periode
                          lmp->force->numeric(arg[14+2*order]));
           }
        }
        return 0;
    }

} /* LAMMPS_NS */

#endif /* MESHMOVER_H_ */
