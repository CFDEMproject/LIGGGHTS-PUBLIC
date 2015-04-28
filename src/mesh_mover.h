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

    Christoph Kloss (DCS Computing GmbH, Linz, JKU Linz)
    Philippe Seil (JKU Linz)
    Niels Dallinger (TU Chemnitz, viblin and vibrot)
    Christian Richter (OVGU Magdeburg, linear variable and rotate/variable)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
    Copyright 2013      TU Chemnitz
    Copyroght 2013      OVGU Magdeburg
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

        virtual void post_create() = 0;
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
        void post_create();

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
        void post_create();
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
        void post_create();

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
        void post_create();

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
        void post_create();
        void setup();

        void initial_integrate(double dTAbs,double dTSetup,double dt);
        void final_integrate(double dTAbs,double dTSetup,double dt) {}

        int n_restart();
        void write_restart(double *buf);
        void read_restart(double *buf);

      private:

        char *var1str_;
        int myvar1_;
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
        void post_create();

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
        void post_create();

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
        void post_create();

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
                          lmp->force->numeric(FLERR,arg[1]),
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]));
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
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          // axis
                          lmp->force->numeric(FLERR,arg[6]),
                          lmp->force->numeric(FLERR,arg[7]),
                          lmp->force->numeric(FLERR,arg[8]),
                          // period
                          lmp->force->numeric(FLERR,arg[10]));
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
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          // axis
                          lmp->force->numeric(FLERR,arg[6]),
                          lmp->force->numeric(FLERR,arg[7]),
                          lmp->force->numeric(FLERR,arg[8]),
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
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          //period
                          lmp->force->numeric(FLERR,arg[6]));
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
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          // axis
                          lmp->force->numeric(FLERR,arg[6]),
                          lmp->force->numeric(FLERR,arg[7]),
                          lmp->force->numeric(FLERR,arg[8]),
                          // period
                          lmp->force->numeric(FLERR,arg[10]),
                          // amplitude
                          lmp->force->numeric(FLERR,arg[12]));
          }
        }
         else if(strcmp(name,"viblin") == 0){
            int order = lmp->force->numeric(FLERR,arg[6]);
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
                amp[zv] = lmp->force->numeric(FLERR,arg[8+zv]);
                // angle of phase
                pha[zv] = lmp->force->numeric(FLERR,arg[9+order+zv]);
               }

            return new MeshMoverVibLin(lmp,mesh,fix_mm,
                          // direction
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          // order
                          lmp->force->numeric(FLERR,arg[6]),
                          // amplitudes
                          amp,
                          // phases
                          pha,
                          // periode
                          lmp->force->numeric(FLERR,arg[10+2*order]));
          }
        }
        else if(strcmp(name,"vibrot") == 0){
             int order = lmp->force->numeric(FLERR,arg[10]);
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
                amp[zv] = lmp->force->numeric(FLERR,arg[12+zv]);
                // angle of phase
                pha[zv] = lmp->force->numeric(FLERR,arg[13+order+zv]);
               }
            return new MeshMoverVibRot(lmp,mesh,fix_mm,
                          // origin px py pz
                          lmp->force->numeric(FLERR,arg[2]),
                          lmp->force->numeric(FLERR,arg[3]),
                          lmp->force->numeric(FLERR,arg[4]),
                          // axis axisx axisy axisz
                          lmp->force->numeric(FLERR,arg[6]),
                          lmp->force->numeric(FLERR,arg[7]),
                          lmp->force->numeric(FLERR,arg[8]),
                          // order
                          lmp->force->numeric(FLERR,arg[10]),
                          // amplitudes
                          amp,
                          // phases
                          pha,
                          // periode
                          lmp->force->numeric(FLERR,arg[14+2*order]));
           }
        }
        return 0;
    }

} /* LAMMPS_NS */

#endif /* MESHMOVER_H_ */
