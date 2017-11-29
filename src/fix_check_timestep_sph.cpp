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

    Andreas Aigner (JKU Linz))

    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include <cmath>
#include "error.h"
#include "fix_check_timestep_sph.h"
#include "fix_property_global.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_wall_gran.h"
#include "fix_mesh_surface.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "mpi_liggghts.h"
#include "sph_kernels.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1000000.

/* ---------------------------------------------------------------------- */

FixCheckTimestepSph::FixCheckTimestepSph(LAMMPS *lmp, int narg, char **arg) :
  FixSph(lmp, narg, arg)
{
  int iarg = 3;
  if (narg < iarg+2) error->all(FLERR,"Illegal fix check/timestep/sph command, not enough arguments");

  nevery = atoi(arg[iarg]);
  fraction_courant_lim=atof(arg[iarg+1]);

  iarg += 2;

  warnflag = true;
  if(iarg < narg){
      if (narg < iarg+2) error->all(FLERR,"Illegal fix check/timestep/sph command, not enough arguments");
      if(strcmp(arg[iarg],"warn")!=0) error->all(FLERR,"Illegal fix check/timestep/sph command, use keyword 'warn'");
      if(strcmp(arg[iarg+1],"no")==0) warnflag=false;
  }

  iarg += 2;

  vector_flag = 1;
  size_vector = 2;
  global_freq = nevery;
  extvector = 1;

  fraction_courant = fraction_skin = 0.;
}

/* ---------------------------------------------------------------------- */

int FixCheckTimestepSph::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepSph::init()
{
  FixSph::init();

  //some error checks
  if(!atom->rho_flag) error->all(FLERR,"Fix check/timestep/sph can only be used together with a sph atom style");

/*
  fwg = NULL;
  for (int i = 0; i < modify->n_fixes_style("wall/gran"); i++)
      if(static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",i))->is_mesh_wall())
        fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",i));
*/

  int ntype = atom->ntypes;
  cs = static_cast<FixPropertyGlobal*>(modify->find_fix_property("speedOfSound","property/global","peratomtype",ntype,0,"check/sph/timestep"));

  if(!cs) error->all(FLERR,"Fix check/timestep/sph only works with a pair style that defines speedOfSound");

}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepSph::end_of_step()
{
    calc_courant_estims();

    double skin = neighbor->skin;
    double dt = update->dt;

    fraction_courant = dt/courant_time;
    fraction_skin = (vmax * dt) / skin;

    if(warnflag&&comm->me==0)
    {
        if(fraction_skin > 0.1)
        {
            if(screen)  fprintf(screen ,"WARNING: time step too large or skin too small - particles may travel a relative distance of %g per time-step, but 0.1 * skin is %g\n",vmax*dt,0.1*skin);
            if(logfile) fprintf(logfile,"WARNING: time step too large or skin too small - particles may travel a relative distance of %g per time-step, but 0.1 * skin is %g\n",vmax*dt,0.1*skin);
        }

        if(fraction_courant > fraction_courant_lim)
        {
            if(screen) fprintf(screen,  "WARNING: time-step is %f %% of courant time\n",fraction_courant*100.);
            if(logfile) fprintf(logfile,"WARNING: time-step is %f %% of courant time\n",fraction_courant*100.);
        }

/*
        if(vmax * dt > r_min)
        {
            if(screen)  fprintf(screen  ,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
            if(logfile)  fprintf(logfile,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
        }
*/
    }
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepSph::calc_courant_estims()
{
  //template function for using per atom or per atomtype smoothing length
  if (mass_type) calc_courant_estims_eval<1>();
  else calc_courant_estims_eval<0>();

}

/* ----------------------------------------------------------------------
   return fraction of courant time-step
------------------------------------------------------------------------- */

double FixCheckTimestepSph::compute_vector(int n)
{
  if(n == 0)      return fraction_courant;
  else if(n == 1) return fraction_skin;
  return 0.;
}

/* ---------------------------------------------------------------------- */

template <int MASSFLAG>
void FixCheckTimestepSph::calc_courant_estims_eval()
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double delvx,delvy,delvz,mu;
  double sli,slj,slCom,cut;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double vmag,courant_time_one;
  double cmean;
  int j_maxmu = 0;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  updatePtrs(); // get sl

  vmax = -1;
  mumax = -1;
  courant_time = BIG;

  // calculate minimum courant time step
  //inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii=0; ii<nlocal;ii++){
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (MASSFLAG) {
      itype = type[i];
      sli = sl[itype-1];
    } else {
      sli = sl[i];
    }

    vmag = sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
    if(vmag > vmax) vmax=vmag;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (!(mask[j] & groupbit)) continue;

      if (MASSFLAG) {
        jtype = type[j];
        slCom = slComType[itype][jtype];
      } else {
        slj = sl[j];
        slCom = interpDist(sli,slj);
      }

      cut = slCom*kernel_cut;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut*cut) {

        delvx = v[i][0] - v[j][0];
        delvy = v[i][1] - v[j][1];
        delvz = v[i][2] - v[j][2];

        mu = slCom * (delvx*delx+delvy*dely+delvz*delz) / (rsq);

        if (mu > mumax) {
          mumax = mu;
          j_maxmu = j;
        }

      }
    } //neighbor loop

    cmean = 0.5*(cs->get_values()[type[i]-1]+cs->get_values()[type[j_maxmu]-1]);
    courant_time_one = sli / (cmean + mumax);
    courant_time = MIN(courant_time,courant_time_one);

  }

  MPI_Max_Scalar(vmax,world);
  MPI_Max_Scalar(courant_time,world);
  /*
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
  vmax = fmax(2.*vmax,vmax+vmax_mesh);

*/

}
