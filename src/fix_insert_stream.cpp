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
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2015 JKU Linz
------------------------------------------------------------------------- */
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "fix_insert_stream.h"
#include "fix_mesh_surface.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_atom_tracer_stream.h"
#include "fix_particledistribution_discrete.h"
#include "fix_multisphere.h"
#include "fix_template_sphere.h"
#include "particleToInsert.h"
#include "tri_mesh_planar.h"

enum{FACE_NONE,FACE_MESH,FACE_CIRCLE};

using namespace LAMMPS_NS;
using namespace FixConst;

#define FIX_INSERT_NTRY_SUBBOX 500
#define FIX_INSERT_STREAM_TINY 1e-14

/* ---------------------------------------------------------------------- */

FixInsertStream::FixInsertStream(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg),
  recalc_release_ms(false),
  dt_ratio(0.),
  save_template_(false),
  fix_template_(NULL)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    
    if (strcmp(arg[iarg],"insertion_face") == 0)
    {
      
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int f_i = modify->find_fix(arg[iarg+1]);
      if (f_i == -1) error->fix_error(FLERR,this,"Could not find fix mesh/surface id you provided");
      if (strncmp(modify->fix[f_i]->style,"mesh",4))
        error->fix_error(FLERR,this,"The fix belonging to the id you provided is not of type mesh");
      ins_face = (static_cast<FixMeshSurface*>(modify->fix[f_i]))->triMesh();
      ins_face->useAsInsertionMesh(false);
      face_style = FACE_MESH;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"extrude_length") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      extrude_length = atof(arg[iarg+1]);
      if(extrude_length < 0. ) error->fix_error(FLERR,this,"invalid extrude_length");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"duration") == 0 || strcmp(arg[iarg],"duration_time") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      if(strcmp(arg[iarg],"duration_time") == 0)
          duration = static_cast<int>(atof(arg[iarg+1])/update->dt);
      else
        duration = atoi(arg[iarg+1]);
      if(duration < 1 ) error->fix_error(FLERR,this,"'duration' can not be < 1");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"parallel") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      if(strcmp("yes",arg[iarg+1]) == 0)
        parallel = true;
      else if(strcmp("no",arg[iarg+1]) == 0)
        parallel = false;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'parallel'");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"ntry_mc") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      ntry_mc = atoi(arg[iarg+1]);
      if(ntry_mc < 1000) error->fix_error(FLERR,this,"ntry_mc must be > 1000");
      iarg += 2;
      hasargs = true;
    }
    else if (strcmp(arg[iarg], "save_template") == 0)
    {
        if (iarg+2 > narg)
            error->fix_error(FLERR,this,"not enough arguments");

        if(strcmp("yes",arg[iarg+1]) == 0)
            save_template_ = true;
        else if(strcmp("no",arg[iarg+1]) == 0)
            save_template_ = false;
        else
            error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'save_template'");
        iarg += 2;
        hasargs = true;
    }
    else if (0 == strcmp(style,"insert/stream")) 
      error->fix_error(FLERR,this,"unknown keyword or wrong keyword order");
  }

  fix_release = NULL;
  i_am_integrator = false;

  tracer = NULL;
  ntracer = 0;

  ins_fraction = 0.;
  do_ins_fraction_calc = true;

  nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixInsertStream::~FixInsertStream()
{
    if(tracer) delete []tracer;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::post_create()
{
    FixInsert::post_create();

    // only register property if I am the first fix/insert/stream in the simulation
    
    if(modify->n_fixes_style(style) == 1)
    {
        const char* fixarg[22];
        fixarg[0]="release_fix_insert_stream";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="release_fix_insert_stream";
        fixarg[4]="vector"; 
        fixarg[5]="yes";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fixarg[11]="0.";
        fixarg[12]="0.";
        fixarg[13]="0.";
        fixarg[14]="0.";
        fixarg[15]="0.";
        fixarg[16]="0.";
        fixarg[17]="0.";
        fixarg[18]="0.";
        fixarg[19]="0.";
        fixarg[20]="0.";
        fixarg[21]="0.";
        modify->add_fix_property_atom(22,const_cast<char**>(fixarg),style);

        fix_release = static_cast<FixPropertyAtom*>(modify->find_fix_property("release_fix_insert_stream","property/atom","vector",14,0,style));
        if(!fix_release) error->fix_error(FLERR,this,"Internal error in fix insert/stream");

        if( modify->fix_restart_in_progress())
            recalc_release_restart();
    }

    if (save_template_)
    {
        fix_template_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("insertion_template_", "property/atom", "scalar", 1, 0, style, false));
        if (!fix_template_)
        {
            const char *fixarg[9];
            fixarg[0] = "insertion_template_";
            fixarg[1] = "all";
            fixarg[2] = "property/atom";
            fixarg[3] = "insertion_template_";
            fixarg[4] = "scalar";
            fixarg[5] = "yes"; // restart
            fixarg[6] = "yes"; // ghost
            fixarg[7] = "no";  // reverse
            fixarg[8] = "-1.0";
            fix_template_ = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), style);
        }
        fix_distribution->save_templates(fix_template_);
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::pre_delete(bool unfixflag)
{
    // delete if I am the last fix of this style to be deleted
    if(unfixflag && modify->n_fixes_style(style) == 1)
        modify->delete_fix("release_fix_insert_stream");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init_defaults()
{
    face_style = FACE_NONE;
    extrude_length = 0.;

    extrude_length_min = extrude_length_max = 0.;

    duration = 0;

    parallel = false;

    ntry_mc = 100000;

    vel_normal_to_face = false;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::register_tracer_callback(FixPropertyAtomTracerStream* tr)
{
    // just return if I already have this callback
    for(int i = 0; i < ntracer; i++)
        if(tracer[i] == tr) return;

    FixPropertyAtomTracerStream** tracer_new = new FixPropertyAtomTracerStream*[ntracer+1];

    for(int i = 0; i < ntracer; i++)
        tracer_new[i] = tracer[i];

    tracer_new[ntracer] = tr;
    ntracer++;
    delete []tracer;
    tracer = tracer_new;
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixInsertStream::calc_insertion_properties()
{
    double dt,dot,extrude_vec[3],t1[3],t2[3];

    // error check on insertion face
    if(face_style == FACE_NONE)
        error->fix_error(FLERR,this,"must define an insertion face");

    // check properties of insertion face
    if(face_style == FACE_MESH)
    {
        // check if face planar
        if(!ins_face->isPlanar())
            error->fix_error(FLERR,this,"command requires a planar face for insertion");

        if(all_in_flag)
        {
            if(!dynamic_cast<TriMeshPlanar*>(ins_face))
                error->fix_error(FLERR,this,"using all_in yes requires you to use a fix mesh/surface/planar");
        }

        // get normal vector of face 0
        ins_face->surfaceNorm(0,normalvec);

        // flip normal vector so dot product with v_insert is > 0
        dot = vectorDot3D(v_insert,normalvec);
        if(dot < 0) vectorScalarMult3D(normalvec,-1.);

        // calc v normal
        dot = vectorDot3D(v_insert,normalvec);
        vectorCopy3D(normalvec,v_normal);
        vectorScalarMult3D(v_normal,dot);

        double diff[3];
        vectorSubtract3D(v_insert,v_normal,diff);

        if(vectorMag3DSquared(diff) < 1e-6)
            vel_normal_to_face = true;
        else
            vel_normal_to_face = false;

        // error check on v normal
        if(vectorMag3D(v_normal) < 1.e-3)
          error->fix_error(FLERR,this,"insertion velocity projected on face normal is < 1e-3");

        // get reference point on face
        ins_face->node(0,0,p_ref);
    }
    else error->fix_error(FLERR,this,"FixInsertStream::calc_insertion_properties(): Implementation missing");

    // error check on insertion velocity
    if(vectorMag3D(v_insert) < 1e-5)
        error->fix_error(FLERR,this,"insertion velocity too low");

    // further error-checks
    if(insert_every == -1 && extrude_length == 0.)
      error->fix_error(FLERR,this,"must define either 'insert_every' or 'extrude_length'");
    if(insert_every > -1 && extrude_length > 0.)
      error->fix_error(FLERR,this,"must not provide both 'insert_every' and 'extrude_length'");
    if(extrude_length > 0. && duration > 0)
      error->fix_error(FLERR,this,"must not provide both 'extrude_length' and 'duration'");

    dt = update->dt;

    // if extrude_length given, calculate insert_every
    if(insert_every == -1)
    {
        // no duration allowed here (checked before)

        if(extrude_length < 3.*max_r_bound() && (all_in_flag || check_ol_flag))
            error->fix_error(FLERR,this,"'extrude_length' is too small");
        // add TINY for resolving round-off
        insert_every = static_cast<int>((extrude_length+FIX_INSERT_STREAM_TINY)/(dt*vectorMag3D(v_normal)));
        
        if(insert_every == 0)
          error->fix_error(FLERR,this,"insertion velocity too high or extrude_length too low");
    }
    // if insert_every given, calculate extrude_length
    // take into account duration can be != insert_every
    else
    {
        if(insert_every < 1) error->fix_error(FLERR,this,"'insert_every' must be > 0");

        // duration = insert_every by default (if already > 0, defined directly)
        if(duration == 0) duration = insert_every;
        else if (duration > insert_every) error->fix_error(FLERR,this,"'duration' > 'insert_every' not allowed");

        extrude_length = static_cast<double>(duration) * dt * vectorMag3D(v_normal);
        
        if(extrude_length < 2.*max_r_bound())
          error->fix_error(FLERR,this,"'insert_every' or 'vel' is too small, or (bounding) radius of inserted particles too large");
    }

    // ninsert - if ninsert not defined directly, calculate it
    if(ninsert == 0 && ninsert_exists)
    {
        if(massinsert/fix_distribution->mass_expect() > 2.e9)
           error->fix_error(FLERR,this,"you are attempting to insert more than 2e9 particles. Reduce the mass to be inserted or increase the particle diameter");

        if(massinsert > 0.) ninsert = static_cast<int>((massinsert+FIX_INSERT_STREAM_TINY) / fix_distribution->mass_expect());
        else error->fix_error(FLERR,this,"must define either 'nparticles' or 'mass'");
    }

    // flow rate
    if(nflowrate == 0.)
    {
        if(massflowrate == 0.) error->fix_error(FLERR,this,"must define either 'massrate' or 'particlerate'");
        nflowrate = massflowrate / fix_distribution->mass_expect();
    }
    else massflowrate = nflowrate * fix_distribution->mass_expect();

    // ninsert_per and massinsert
    ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);
    if(ninsert_exists) massinsert = static_cast<double>(ninsert) * fix_distribution->mass_expect();

    // calculate bounding box of extruded face
    if(face_style == FACE_MESH)
    {
        // get bounding box for face
        ins_face->getGlobalBoundingBox().getBoxBounds(ins_vol_xmin,ins_vol_xmax);

        // get bounding box for extruded face - store in t1,t2
        vectorScalarMult3D(normalvec,-extrude_length,extrude_vec);
        vectorAdd3D(ins_vol_xmin,extrude_vec,t1);
        vectorAdd3D(ins_vol_xmax,extrude_vec,t2);

        // take min and max
        vectorComponentMin3D(ins_vol_xmin,t1,ins_vol_xmin);
        vectorComponentMax3D(ins_vol_xmax,t2,ins_vol_xmax);

    }
    else error->fix_error(FLERR,this,"Missing implementation in calc_insertion_properties()");

    extrude_length_min = 0.;
    extrude_length_max = extrude_length;
}

/* ---------------------------------------------------------------------- */

int FixInsertStream::setmask()
{
    int mask = FixInsert::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init()
{
    
    FixInsert::init();

    if(fix_multisphere && v_randomSetting != RANDOM_CONSTANT)
        error->fix_error(FLERR,this,"Currently only fix insert/stream with multisphere particles only supports constant velocity");

    fix_release = static_cast<FixPropertyAtom*>(modify->find_fix_property("release_fix_insert_stream","property/atom","vector",5,0,style));
    if(!fix_release) error->fix_error(FLERR,this,"Internal error if fix insert/stream");
    fix_release->set_internal();

    i_am_integrator = modify->i_am_first_of_style(this);

    // error check on insertion face
    if(face_style == FACE_NONE)
        error->fix_error(FLERR,this,"must define an insertion face");

    if(ins_face->isMoving() || ins_face->isScaling())
        error->fix_error(FLERR,this,"cannot translate, rotate, scale mesh which is used for particle insertion");

    if(recalc_release_ms)
    {
        recalc_release_ms = false;
        
        if(fix_multisphere && dt_ratio > 0.)
            fix_multisphere->data().recalc_n_steps(dt_ratio);
        
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::setup_pre_exchange()
{

}

/* ---------------------------------------------------------------------- */

double FixInsertStream::insertion_fraction()
{
    
    // have to re-calculate insertion fraction for my subbox
    // in case subdomains of simulation box are changing
    
    if(domain->box_change || do_ins_fraction_calc || ins_face->isMoving())
        calc_ins_fraction();

    return ins_fraction;
}

/* ----------------------------------------------------------------------
   calculate insertion fraction for my subbox
   has to be called at initialization and before every insertion in case
   box is changing
------------------------------------------------------------------------- */

void FixInsertStream::calc_ins_fraction()
{
    
    do_ins_fraction_calc = false;

    double pos[3], boxedgevec[3], dot;
    int n_in_local = 0, n_test = ntry_mc;

    for(int i = 0; i < n_test; i++)
    {
        generate_random_global(pos);

        if(domain->is_in_subdomain(pos))
           n_in_local++;
    }

    ins_fraction = static_cast<double>(n_in_local)/static_cast<double>(n_test);

    // also calculate min and max extrusion
    // this can speed up insertion if extrusion volume extends across multiple procs

    if(parallel)
    {
        extrude_length_min = extrude_length;
        extrude_length_max = 0.;

        for(int ix = 0; ix < 2; ix++)
            for(int iy = 0; iy < 2; iy++)
                for(int iz = 0; iz < 2; iz++)
                {
                    vectorConstruct3D
                    (
                        boxedgevec,
                        (ix == 0 ? domain->sublo[0] : domain->subhi[0]) - p_ref[0],
                        (iy == 0 ? domain->sublo[1] : domain->subhi[1]) - p_ref[1],
                        (iz == 0 ? domain->sublo[2] : domain->subhi[2]) - p_ref[2]
                    );

                    dot = -vectorDot3D(boxedgevec,normalvec);
                    
                    if(dot > 0. && dot < extrude_length)
                    {
                        extrude_length_max = std::max(extrude_length_max,dot);
                        extrude_length_min = std::min(extrude_length_min,dot);
                    }
                    else if(dot < 0.)
                        extrude_length_min = 0.;
                    else if(dot >= extrude_length)
                        extrude_length_max = extrude_length;
                }
        if(extrude_length_min == extrude_length)
            extrude_length_min = 0.;
        if(extrude_length_max == 0.)
            extrude_length_max = extrude_length;
    }

    double ins_fraction_all;
    MPI_Sum_Scalar(ins_fraction,ins_fraction_all,world);
    if(ins_fraction_all < 0.9 || ins_fraction_all > 1.1)
        error->fix_error(FLERR,this,"insertion volume could not be distributed properly in parallel. "
                                     "Bad decomposition or insertion face extrusion is too small or outside domain");
}

/* ---------------------------------------------------------------------- */

bool FixInsertStream::pre_insert()
{
    if((!domain->is_in_domain(ins_vol_xmin) || !domain->is_in_domain(ins_vol_xmax)) && comm->me == 0)
      error->warning(FLERR,"Fix insert/stream: Extruded insertion face extends outside domain, may not insert all particles correctly");

    return true;
}

/* ---------------------------------------------------------------------- */

inline int FixInsertStream::is_nearby(int i)
{
    double pos_rel[3], pos_projected[3], t[3];
    double **x = atom->x;

    vectorSubtract3D(x[i],p_ref,pos_rel);
    double dist_normal = vectorDot3D(pos_rel,normalvec);

    // on wrong side of extrusion
    if(dist_normal > maxrad) return 0;

    // on right side of extrusion, but too far away
    // 3*maxrad b/c extrude_length+rad is max extrusion for overlapcheck yes
    if(dist_normal < -(extrude_length + 3.*maxrad)) return 0;

    // on right side of extrusion, within extrude_length
    // check if projection is on face or not

    vectorScalarMult3D(normalvec,dist_normal,t);
    vectorAdd3D(x[i],t,pos_projected);

    //TODO also should check if projected point is NEAR surface

    return ins_face->isOnSurface(pos_projected);
}

/* ---------------------------------------------------------------------- */

BoundingBox FixInsertStream::getBoundingBox() {
  BoundingBox bb = ins_face->getGlobalBoundingBox();

  const double cut = 3*maxrad;
  const double delta = -(extrude_length + 2*cut);
  bb.extrude(delta, normalvec);
  bb.shrinkToSubbox(domain->sublo, domain->subhi);

  const double extend = 3*maxrad /*cut*/ + 2.*fix_distribution->max_r_bound(); 
  bb.extendByDelta(extend);

  return bb;
}

/* ----------------------------------------------------------------------
   generate random positions on insertion face
   extrude by random length in negative face normal direction
     currently only implemented for all_in_flag = 0
     since would be tedious to check/implement otherwise
------------------------------------------------------------------------- */

inline void FixInsertStream::generate_random(double *pos, double rad)
{
    double r, ext[3];

    // generate random position on the mesh
    
    if(all_in_flag)
        ins_face->generateRandomOwnedGhostWithin(pos,rad);
        
    else
        ins_face->generateRandomOwnedGhost(pos);
        
    // extrude the position
    
    if(check_ol_flag)
        r = -1.*(random->uniform()*(extrude_length_max         ) + rad + extrude_length_min);
    else
        r = -1.*(random->uniform()*(extrude_length_max - 2.*rad) + rad + extrude_length_min);

    vectorScalarMult3D(normalvec,r,ext);
    vectorAdd3D(pos,ext,pos);

}

/* ----------------------------------------------------------------------
   generate random positions on shallow copy insertion face
   extrude by random length in negative face normal direction
     currently only implemented for all_in_flag = 0
     since would be tedious to check/implement otherwise
------------------------------------------------------------------------- */

inline void FixInsertStream::generate_random_global(double *pos)
{
    double r, ext[3];

    // generate random position on the mesh
    
    ins_face->generateRandomOwnedGhost(pos);

    // extrude the position
    r = -1.*(random->uniform()*extrude_length);
    vectorScalarMult3D(normalvec,r,ext);
    vectorAdd3D(pos,ext,pos);

}

/* ----------------------------------------------------------------------
   generate random positions within extruded face
   perform overlap check via xnear if requested
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertStream::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
    ninserted_this_local = ninserted_spheres_this_local = 0;
    mass_inserted_this_local = 0.;

    int nins;
    double pos[3];
    ParticleToInsert *pti;

    double omega_tmp[] = {0.,0.,0.};

    int ntry = 0;
    int maxtry = ninsert_this_local * maxattempt;

    // no overlap check
    // insert with v_normal, no omega
    if(!check_ol_flag)
    {
        for(int itotal = 0; itotal < ninsert_this_local; itotal++)
        {
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rad_to_insert = pti->r_bound_ins;

            do
            {
                generate_random(pos,rad_to_insert);
                ntry++;
            }
                        
            while(ntry < maxtry && (!domain->is_in_subdomain(pos)));

            if(ntry < maxtry)
            {
                // randomize quat here

                if(quat_random_)
                        MathExtraLiggghts::random_unit_quat(random,quat_insert);

                nins = pti->set_x_v_omega(pos,v_normal,omega_tmp,quat_insert);

                ninserted_spheres_this_local += nins;
                mass_inserted_this_local += pti->mass_ins;
                ninserted_this_local++;
            }
        }
    }
    // overlap check
    // account for maxattempt
    // pti checks against xnear and adds self contributions
    else
    {
        while(ntry < maxtry && ninserted_this_local < ninsert_this_local)
        {
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rad_to_insert = pti->r_bound_ins;

            nins = 0;
            while(nins == 0 && ntry < maxtry)
            {
                do
                {
                    generate_random(pos,rad_to_insert);
                    ntry++;

                }
                while(ntry < maxtry && ((!domain->is_in_subdomain(pos)) || (domain->dist_subbox_borders(pos) < rad_to_insert)));

                if(ntry < maxtry)
                {
                    // randomize quat here
                    if(quat_random_)
                        MathExtraLiggghts::random_unit_quat(random,quat_insert);
                    
                    nins = pti->check_near_set_x_v_omega(pos,v_normal,omega_tmp,quat_insert,neighList);
                }
            }

            if(nins > 0)
            {
                ninserted_spheres_this_local += nins;
                mass_inserted_this_local += pti->mass_ins;
                ninserted_this_local++;
            }
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::finalize_insertion(int ninserted_spheres_this_local)
{
    // nins particles have been inserted on this proc, set initial position, insertion step and release step according to pos
    
    int n_steps = -1;
    int step = update->ntimestep;
    int ilo = atom->nlocal - ninserted_spheres_this_local;
    int ihi = atom->nlocal;

    double pos_rel[3], dist_normal;
    double **x = atom->x;
    double dt = update->dt;

    double **release_data = fix_release->array_atom;

    Multisphere *multisphere = NULL;
    if(fix_multisphere) multisphere = &fix_multisphere->data();

    for(int i = ilo; i < ihi; i++)
    {
        
        if(multisphere)
            n_steps = fix_multisphere->calc_n_steps(i,p_ref,normalvec,v_normal);
        if(!multisphere || n_steps == -1)
        {
            vectorSubtract3D(p_ref,x[i],pos_rel);
            dist_normal = vectorDot3D(pos_rel,normalvec);
            n_steps = static_cast<int>((dist_normal+FIX_INSERT_STREAM_TINY)/(vectorMag3D(v_normal)*dt));
        }

        // first 3 values is original position to integrate
        vectorCopy3D(x[i],release_data[i]);

        // 4th value is insertion step
        release_data[i][3] = static_cast<double>(step);

        // 5th value is step to release
        release_data[i][4] = static_cast<double>(step + n_steps);

        // 6-8th value is integration velocity
        vectorCopy3D(v_normal,&release_data[i][5]);

        // set initial conditions
        // randomize vel, omega here
        double v_toInsert[3],omega_toInsert[3];

        vectorCopy3D(v_insert,v_toInsert);
        vectorCopy3D(omega_insert,omega_toInsert);

        // could randomize vel, omega here
        generate_random_velocity(v_toInsert);

        // 9-11th value is velocity, 12-14 is omega
        vectorCopy3D(v_toInsert,&release_data[i][8]);
        vectorCopy3D(omega_toInsert,&release_data[i][11]);

    }

    if(ntracer)
        for(int i = 0; i < ntracer; i++)
            tracer[i]->mark_tracers(ilo,ihi);
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::end_of_step()
{
    int r_step, i_step;

    int step = update->ntimestep;
    int nlocal = atom->nlocal;
    double **release_data = fix_release->array_atom;
    double time_elapsed, dist_elapsed[3], v_integrate[3], *v_toInsert, *omega_toInsert;
    double dt = update->dt;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int *mask = atom->mask;

    // only one fix handles the integration
    if(!i_am_integrator) return;

    for(int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if(MathExtraLiggghts::compDouble(release_data[i][3],0.,1.e-13))
                continue;

            i_step = static_cast<int>(release_data[i][3]+FIX_INSERT_STREAM_TINY);
            r_step = static_cast<int>(release_data[i][4]+FIX_INSERT_STREAM_TINY);
            vectorCopy3D(&release_data[i][5],v_integrate);

            if(step > r_step) continue;
            else if (r_step == step)
            {
                
                if(fix_multisphere && fix_multisphere->belongs_to(i) >= 0)
                {
                    
                    v_toInsert = &release_data[i][8];
                    omega_toInsert = &release_data[i][11];
                    fix_multisphere->release(i,v_toInsert,omega_toInsert);
                    continue;
                }

                // integrate with constant vel and set v,omega

                time_elapsed = (step - i_step) * dt;

                // particle moves with v_integrate
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                double *x_ins = release_data[i];

                // set x,v,omega
                // zero out force, torque

                vectorAdd3D(x_ins,dist_elapsed,x[i]);

                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);

                v_toInsert = &release_data[i][8];
                omega_toInsert = &release_data[i][11];

                vectorCopy3D(v_toInsert,v[i]);
                vectorCopy3D(omega_toInsert,omega[i]);

            }
            // step < r_step, only true for inserted particles
            //   b/c r_step is 0 for all other particles
            // integrate with constant vel
            else
            {
                
                time_elapsed = (step - i_step) * dt;

                // particle moves with v_integrate
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                double *x_ins = release_data[i];

                // set x,v,omega
                vectorAdd3D(x_ins,dist_elapsed,x[i]);
                
                vectorCopy3D(v_integrate,v[i]);
                vectorZeroize3D(omega[i]);

                // zero out force, torque
                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);
            }
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::reset_timestep(bigint newstep,bigint oldstep)
{
    FixInsert::reset_timestep(newstep,oldstep);
    
    reset_releasedata(newstep,oldstep);
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::reset_releasedata(bigint newstep,bigint oldstep)
{
  
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **release_data = fix_release->array_atom;

  for(int i = 0; i < nlocal; i++)
  {
        
        // first 3 values is original position to integrate
        vectorCopy3D(x[i],release_data[i]);

        // 4th value is insertion step
        release_data[i][3] -= static_cast<double>(oldstep-newstep);

        // 5th value is step to release
        release_data[i][4] -= static_cast<double>(oldstep-newstep);

        // 6-8th value is integration velocity
        vectorCopy3D(v_normal,&release_data[i][5]);
  }
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::recalc_release_restart()
{
  
  const int nlocal = atom->nlocal;
  double **x = atom->x;
  double **release_data = fix_release->array_atom;
  const double dt = update->dt;
  double change_ratio = -1.;

  for(int i = 0; i < nlocal; i++)
  {
        if(release_data[i][4] > update->ntimestep)
        {
            double dx[3];
            vectorSubtract3D(x[i],release_data[i],dx);
            const double dt_old = (vectorMag3D(dx)) / (vectorMag3D(&release_data[i][5]) * (static_cast<double>(update->ntimestep) - release_data[i][3]));

            change_ratio = dt_old / dt;
            
            release_data[i][4] = static_cast<double>(update->ntimestep) + static_cast<double>(static_cast<int>(change_ratio*(release_data[i][4] - static_cast<double>(update->ntimestep))));
            
        }
  }

  recalc_release_ms = true;
  dt_ratio = change_ratio;
}
