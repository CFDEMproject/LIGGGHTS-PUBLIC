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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "math.h"
#include "error.h"
#include "fix_check_timestep_gran.h"
#include "pair_gran.h"
#include "properties.h"
#include "fix_property_global.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_wall_gran.h"
#include "fix_mesh_surface.h"
#include "neighbor.h"
#include "mpi_liggghts.h"
#include "property_registry.h"
#include "global_properties.h"
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MODEL_PARAMS;

#define BIG 1000000.

/* ---------------------------------------------------------------------- */

FixCheckTimestepGran::FixCheckTimestepGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix check/timestep/gran command, not enough arguments");

  nevery = atoi(arg[3]);
  fraction_rayleigh_lim = atof(arg[4]);
  fraction_hertz_lim = atof(arg[5]);

  int iarg = 6;

  warnflag = true;
  if(iarg < narg){
      if (narg < 8) error->all(FLERR,"Illegal fix check/timestep/gran command, not enough arguments");
      if(strcmp(arg[iarg++],"warn")!=0) error->all(FLERR,"Illegal fix check/timestep/gran command, use keyword 'warn'");
      if(strcmp(arg[iarg++],"no")==0) warnflag=false;
  }

  vector_flag = 1;
  size_vector = 3;
  global_freq = nevery;
  extvector = 1;

  fraction_rayleigh = fraction_hertz = fraction_skin = 0.;
  Yeff = NULL;
}

/* ---------------------------------------------------------------------- */

int FixCheckTimestepGran::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::init()
{
  //some error checks
  if(!atom->radius_flag || !atom->density_flag)
    error->all(FLERR,"Fix check/timestep/gran can only be used together with atom style sphere");

  pg = (PairGran*)force->pair_match("gran",1);
  if(!pg) pg = (PairGran*)force->pair_match("gran/omp",1);

  if (!pg)
    error->all(FLERR,"Fix check/timestep/gran can only be used together with: gran");

  properties = atom->get_properties();
  int max_type = properties->max_type();

  fwg = NULL;
  for (int i = 0; i < modify->n_fixes_style("wall/gran"); i++)
      if(static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",i))->is_mesh_wall())
        fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",i));

  Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style));
  nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style));

  if(!Y || !nu)
    error->all(FLERR,"Fix check/timestep/gran only works with a pair style that defines youngsModulus and poissonsRatio");

  force->registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
  force->registry.connect("Yeff", Yeff,this->style);
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::end_of_step()
{
    calc_rayleigh_hertz_estims();

    double skin = neighbor->skin;
    double dt = update->dt;

    fraction_rayleigh = dt/rayleigh_time;
    fraction_hertz = dt/hertz_time;
    fraction_skin = (vmax * dt) / neighbor->skin;

    if(warnflag&&comm->me==0)
    {
        if(fraction_skin >= 0.5)
        {
            if(screen)  fprintf(screen ,"WARNING: time step too large or skin too small - particles may travel a distance of %f per time-step, but half skin is %f\n",vmax*dt,0.5*skin);
            if(logfile) fprintf(logfile,"WARNING: time step too large or skin too small - particles may travel a distance of %f per time-step, but half skin is %f\n",vmax*dt,0.5*skin);
        }

        if(vmax * dt > r_min)
        {
            if(screen)  fprintf(screen  ,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
            if(logfile)  fprintf(logfile,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
        }

        if(fraction_rayleigh > fraction_rayleigh_lim)
        {
            if(screen) fprintf(screen,  "WARNING: time-step is %f %% of rayleigh time\n",fraction_rayleigh*100.);
            if(logfile) fprintf(logfile,"WARNING: time-step is %f %% of rayleigh time\n",fraction_rayleigh*100.);
        }
        if(fraction_hertz > fraction_hertz_lim)
        {
            if(screen) fprintf(screen,  "WARNING: time-step is %f %% of hertz time\n",fraction_hertz*100.);
            if(logfile) fprintf(logfile,"WARNING: time-step is  %f %% of hertz time\n",fraction_hertz*100.);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::calc_rayleigh_hertz_estims()
{
  double **v = atom->v;
  double *density = atom->density;
  double *r = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int max_type = properties->max_type();

  //check rayleigh time and vmax of particles
  rayleigh_time = BIG;
  r_min = BIG;
  vmax = 0;

  double vmag;
  double rayleigh_time_i;

  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        double shear_mod = Y->get_values()[type[i]-1]/(2.*(nu->get_values()[type[i]-1]+1.));
        rayleigh_time_i = M_PI*r[i]*sqrt(density[i]/shear_mod)/(0.1631*nu->get_values()[type[i]-1]+0.8766);
        if(rayleigh_time_i < rayleigh_time) rayleigh_time = rayleigh_time_i;

        vmag = sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
        if(vmag > vmax) vmax=vmag;

        if(r[i] < r_min) r_min = r[i];
    }
  }

  MPI_Min_Scalar(r_min,world);
  MPI_Max_Scalar(vmax,world);
  MPI_Min_Scalar(rayleigh_time,world);

  // get vmax of geometry
  FixMeshSurface ** mesh_list;
  TriMesh * mesh;
  double *v_node;
  double vmax_mesh=0.;

  if(fwg)
  {
      mesh_list = fwg->mesh_list();
      for(int imesh = 0; imesh < fwg->n_meshes(); imesh++)
      {
          mesh = (mesh_list[imesh])->triMesh();
          if(mesh->isMoving())
          {
              // check if perElementProperty 'v' exists
              if (mesh->prop().getElementPropertyIndex("v") == -1)
                  error->one(FLERR,"Internal error - mesh has no perElementProperty 'v' \n");
              // loop local elements only
              for(int itri=0;itri<mesh->sizeLocal();itri++)
                  for(int inode=0;inode<3;inode++)
                  {
                      v_node = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->begin()[itri][inode];
                      vmag = vectorMag3D(v_node);
                      if(vmag>vmax_mesh) vmax_mesh=vmag;
                  }
          }
      }
  }

 MPI_Max_Scalar(vmax_mesh,world);

  // decide vmax - either particle-particle or particle-mesh contact
  vmax = std::max(2.*vmax,vmax+vmax_mesh);

  // check estimation for hertz time
  // this is not exact...
  // loop over all material comibinations
  //  loop all particles
  //     test collision of particle with itself
  double hertz_time_min = 1000000.;
  double hertz_time_i,meff,reff;

  for(int ti = 1; ti < max_type+1; ti++)
  {
      for(int tj =  ti; tj < max_type+1; tj++)
      {
          const double Eeff = Yeff[ti][tj];

          for(int i = 0; i < nlocal; i++)
          {
            if (mask[i] & groupbit)
            {
                if(type[i]!=ti || type[i]!=tj) continue;
                meff = 4.*r[i]*r[i]*r[i]*M_PI/3.*density[i];
                reff = r[i]/2.;
                hertz_time_i = 2.87*pow(meff*meff/(reff*Eeff*Eeff*vmax),0.2);
                if(hertz_time_i<hertz_time_min)
                    hertz_time_min=hertz_time_i;
            }
          }
      }
  }

 MPI_Min_Scalar(hertz_time_min,world);
  hertz_time = hertz_time_min;
}

/* ----------------------------------------------------------------------
   return fractions of rayleigh/hertz time-step
------------------------------------------------------------------------- */

double FixCheckTimestepGran::compute_vector(int n)
{
  if(n == 0)      return fraction_rayleigh;
  else if(n == 1) return fraction_hertz;
  else if(n == 2) return fraction_skin;
  return 0.;
}
