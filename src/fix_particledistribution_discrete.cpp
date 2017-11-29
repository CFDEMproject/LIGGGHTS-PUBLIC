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

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_template_sphere.h"
#include "fix_template_multiplespheres.h"
#include "fix_particledistribution_discrete.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fix_template_(NULL)
{
  restart_global = 1;

  mass_based = true;

  if(strstr(arg[2],"numberbased"))
    mass_based = false;

  // random number generator, same for all procs

  if (narg < 7)
    error->fix_error(FLERR,this,"not enough arguments");
  random = new RanPark(lmp, arg[3], true);
  seed = random->getSeed();
  ntemplates = atoi(arg[4]);
  if(ntemplates < 1)
    error->fix_error(FLERR,this,"illegal number of templates");

  templates = new FixTemplateSphere*[ntemplates];
  distweight = new double[ntemplates];
  cumweight = new double[ntemplates];
  parttogen = new int[ntemplates];
  distorder = new int[ntemplates];

  iarg = 5;

  int itemp=0;

  if(narg != iarg+2*ntemplates)
    error->fix_error(FLERR,this,"# of templates does not match # of arguments");

  // parse further args
  do {
    if(itemp == ntemplates) break;
    if(narg < iarg+1)
        error->fix_error(FLERR,this,"not enough arguments");
    int ifix = modify->find_fix(arg[iarg]);

    if(ifix < 0)
        error->fix_error(FLERR,this,"invalid ID for fix particletemplate provided");

    if(strncmp(modify->fix[ifix]->style,"particletemplate/",16))
        error->fix_error(FLERR,this,"fix is not of type particletemplate");

    templates[itemp] = static_cast<FixTemplateSphere*>(modify->fix[ifix]);
    distweight[itemp] = atof(arg[iarg+1]);
    if (distweight[itemp] < 0) error->fix_error(FLERR,this,"invalid weight");
    itemp++;
    iarg += 2;
  } while (iarg < narg);

  // check for double use of template which is not allowed
  for(int i = 0; i < ntemplates; i++)
      for(int j = 0; j < i; j++)
        if(templates[i] == templates[j])
            error->fix_error(FLERR,this,"cannot use the same template twice");

  // normalize distribution
  double weightsum = 0;
  for(int i = 0; i < ntemplates; i++)
    weightsum += distweight[i];

  if(comm->me == 0 && fabs(weightsum-1.) > 0.00001)
    error->warning(FLERR,"particledistribution/discrete: sum of distribution weights != 1, normalizing distribution");

  for(int i = 0; i < ntemplates; i++)
    distweight[i]/=weightsum;

  if(mass_based && comm->me == 0 && screen)
  {
      fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on mass%%:\n",this->id);
      for(int i = 0; i < ntemplates; i++)
        fprintf(screen,"    %s: d=%e (max. bounding sphere) mass%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
  }

  // convert distribution from mass% to number%
  // do not do if already number-based
  if(mass_based)
  {
      for(int i=0;i<ntemplates; i++)
        distweight[i]=distweight[i]/templates[i]->massexpect();

      weightsum=0;
      for(int i=0;i<ntemplates; i++)
        weightsum+=distweight[i];

      for(int i=0;i<ntemplates; i++)
        distweight[i]/=weightsum;
  }

  if(comm->me == 0 && screen)
  {
      fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on number%%:\n",this->id);
      for(int i = 0; i < ntemplates; i++)
        fprintf(screen,"    %s: d=%e (max. bounding sphere) number%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
  }

  cumweight[0] = distweight[0];
  for(int i = 1; i < ntemplates; i++)
    cumweight[i] = distweight[i]+cumweight[i-1];

  volexpect=0.;massexpect=0.;

  for(int i = 0; i < ntemplates; i++)
  {
      volexpect  += templates[i]->volexpect()  * distweight[i];
      massexpect += templates[i]->massexpect() * distweight[i];
  }

  //get min/maxtype
  maxtype = 0;
  mintype = 10000;
  for(int i = 0; i < ntemplates; i++)
  {
    if(templates[i]->maxtype() > maxtype)
      maxtype = templates[i]->maxtype();
    if(templates[i]->mintype() < mintype)
      mintype = templates[i]->mintype();
  }

  // check which template has the most spheres
  maxnspheres = 0;
  for(int i = 0; i < ntemplates;i++)
    if(templates[i]->number_spheres() > maxnspheres)
      maxnspheres=templates[i]->number_spheres();

  // sort the distributions by insertion volume (in descending order)
  // use bubble sort
  for(int i = 0; i < ntemplates; i++)
    distorder[i]=i;

  bool swaped;
  int n = ntemplates;
  do
  {
      swaped = false;
      for(int i = 0; i < ntemplates-1; i++)
      {
          if(templates[distorder[i]]->volexpect() < templates[distorder[i+1]]->volexpect())
          {
            //swap
            int tmp = distorder[i];
            distorder[i] = distorder[i+1];
            distorder[i+1] = tmp;
            swaped = true;
          }
      }
      n--;
  } while(swaped && n > 0);

  pti = templates[distorder[0]]->pti;

  pti_list = NULL;
  n_pti = n_pti_max = 0;

  //calc max radius and bounding sphere radius

  maxrad = maxrbound = 0.;
  minrad = 1000.;

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->max_r_bound() > maxrbound)
        maxrbound = templates[i]->max_r_bound();

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->max_rad() > maxrad)
        maxrad = templates[i]->max_rad();

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->min_rad() < minrad)
        minrad = templates[i]->min_rad();

}

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::~FixParticledistributionDiscrete()
{
    delete []templates;
    delete []distweight;
    delete []cumweight;
    delete []parttogen;
    delete []distorder;
    if(pti_list) delete []pti_list;
    delete random;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::setmask()
{
    int mask = 0;
    return mask;
}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_single() commands
   typically called once per insertion step
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::random_init_single(int ntotal)
{
    ninsert = ntotal;
    ninserted = 0;

    for(int i = 0; i < ntemplates; i++)
       parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i] + random->uniform());

    ninsert = 0;
    for(int i = 0; i < ntemplates; i++)
        ninsert += parttogen[i];
    return ninsert;
}

/* ----------------------------------------------------------------------
   request one template to generate one pti
------------------------------------------------------------------------- */

Region* FixParticledistributionDiscrete::randomize_single()
{
    if(ntemplates == 1){
         templates[0]->randomize_single();
         
         return templates[0]->region(); 
    }

    //choose a template from the discrete distribution, beginning from large to small particles
    int chosen = 0;
    int chosendist = distorder[chosen];
    int ntoinsert = parttogen[chosendist];
    while(ninserted >= ntoinsert && chosen < ntemplates-1)
    {
        chosen++;
        chosendist = distorder[chosen];
        ntoinsert += parttogen[chosendist];
    }

    templates[chosendist]->randomize_single();

    pti = templates[chosendist]->pti;

    ninserted++;

    return templates[chosendist]->region();

}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_list() command
   also prepares templates
       - deletes their old lists if present and allocates new lists
   typically only called once before first insertion step

   allocates for max # particles

   can be called by multiple fix insert commands, so check first if max #
   particles to be inserted is exceeded and only re-allocate in this case
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::random_init_list(int ntotal)
{
    int parttogen_max_i, n_pti_max_requested;
    int nprocs = comm->nprocs;

    ntotal += 2 * ntemplates;

    // number of requested pti
    n_pti_max_requested = 0;

    for(int i = 0; i < ntemplates; i++)
    {
        parttogen_max_i = static_cast<int>(static_cast<double>(ntotal) * distweight[i] + static_cast<double>(1.01)*(ntemplates+nprocs));
        n_pti_max_requested += parttogen_max_i;

        // re-allocated if need more ptis in this template than allocated so far
        if(parttogen_max_i > templates[i]->n_pti_max)
        {
            templates[i]->delete_ptilist();
            templates[i]->init_ptilist(parttogen_max_i);
        }
    }

    // re-allocate if need more total ptis in distribution than allocated so far
    if(n_pti_max_requested > n_pti_max)
    {
        n_pti_max = n_pti_max_requested;
        if(pti_list) delete []pti_list;
        pti_list = new ParticleToInsert*[n_pti_max];
        
    }

}

void FixParticledistributionDiscrete::direct_init_list(const int * const parttogen, FixPropertyAtom * const fix_release)
{
    int n_pti_max_requested = 0;
    for (int i = 0; i < ntemplates; i++)
    {
        if (parttogen[i] > templates[i]->n_pti_max)
        {
            templates[i]->delete_ptilist();
            templates[i]->init_ptilist(parttogen[i], true, fix_release);
        }
        n_pti_max_requested += parttogen[i];
    }

    // re-allocate if need more total ptis in distribution than allocated so far
    if(n_pti_max_requested > n_pti_max)
    {
        n_pti_max = n_pti_max_requested;
        if(pti_list) delete []pti_list;
        pti_list = new ParticleToInsert*[n_pti_max];
        
    }
}

/* ----------------------------------------------------------------------
   tell all templates to generate their pti_list, wire their pti_list to
   the list in this fix. returns number of particles to be inserted.
   typically called once per insertion step

   for exact_number = 1, truncate distribution so to exactly meet
                               requested # particles
   for exact_number = 0, use random gen to fulfill distribution
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::randomize_list(int ntotal,int insert_groupbit,int exact_number)
{
    if(ntotal > n_pti_max)
    {
        
        error->one(FLERR,"Faulty implementation: FixParticledistributionDiscrete::randomize_list() called for more particles than defined in random_init_list()");
    }

    ninsert = ntotal;
    ninserted = 0;

    // use random generator so long-time average of insertion will represent distribution correctly
    if(exact_number == 0)
    {

        for(int i = 0; i < ntemplates; i++)
        {
           parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i] + random->uniform());
           
        }
    }
    // truncate distribution so # particles to insert is met exactly
    else
    {
        int ninsert_truncated = 0, j;
        double *remainder = new double[ntemplates], rsum, r;

        // distribute particles and calculate remainder
        for(int i = 0; i < ntemplates; i++)
        {
           parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i]);
           ninsert_truncated += parttogen[i];
           remainder[i] = static_cast<double>(ninsert) * distweight[i] - static_cast<double>(parttogen[i]);
           
        }

        int ninsert_gap = ninsert - ninsert_truncated;

        // distribute remaining ninsert_gap particles
        for(int i = 0; i < ninsert_gap; i++)
        {
            r = random->uniform() * static_cast<double>(ninsert_gap);
            j = 0;
            rsum = remainder[0];

            while(rsum < r && j < (ntemplates-1) )
            {
                j++;
                rsum += remainder[j];
            }
            
            parttogen[j]++;
        }

        delete []remainder;
    }

    // count total particle number to be inserted, let templates generate a pti_list
    ninsert = 0;
    for(int i = 0; i < ntemplates; i++)
    {
        ninsert += parttogen[i];
        templates[i]->randomize_ptilist(parttogen[i],groupbit | insert_groupbit,distorder[i]);
    }

    // wire lists, make sure in correct order (large to small particles)

    n_pti = 0;
    for(int i = 0; i < ntemplates; i++)
    {
        int chosendist = distorder[i];
        for (int j = 0; j < parttogen[chosendist]; j++)
        {
            pti_list[n_pti + j] = templates[chosendist]->pti_list[j];
        }
        n_pti += parttogen[chosendist];
    }

    if(n_pti != ninsert)
        error->one(FLERR,"Internal error in FixParticledistributionDiscrete::randomize_list");

    ninserted = ninsert;
    return ninsert;
}

/* ---------------------------------------------------------------------- */

void FixParticledistributionDiscrete::direct_set_ptlist(const int itemplate, const int i, const void * const data, const int distribution_groupbit)
{
    templates[itemplate]->direct_set_ptlist(i, data, distribution_groupbit | groupbit, distorder[itemplate]);
}

/* ---------------------------------------------------------------------- */

int FixParticledistributionDiscrete::update_ptlist_pointer(const int * ext_parttogen)
{
    n_pti = 0;
    ninsert = 0;
    for (int i = 0; i < ntemplates; i++)
    {
        ninsert += ext_parttogen[i];
        int chosendist = distorder[i];
        parttogen[chosendist] = ext_parttogen[chosendist];
        for (int j = 0; j < parttogen[chosendist]; j++)
        {
            pti_list[n_pti + j] = templates[chosendist]->pti_list[j];
        }
        n_pti += parttogen[chosendist];
    }

    if(n_pti != ninsert)
        error->one(FLERR,"Internal error in FixParticledistributionDiscrete::update_ptlist_ptr");

    ninserted = ninsert;
    return ninsert;
}

/* ----------------------------------------------------------------------
   preparations before insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::pre_insert(int n,class FixPropertyAtom *fp,double val)
{
    // allow fixes to e.g. update some pointers before set_arrays is called
    // set_arrays called in ParticleToInsert::insert()

    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    for (int j = 0; j < nfix; j++)
        if (fix[j]->create_attribute) fix[j]->pre_set_arrays();

    // set fix property as desired by fix insert
    // loop to n, not n_pti
    if(fp)
    {
        
        for(int i = 0; i < ntemplates; i++)
        {
            if( dynamic_cast<FixTemplateMultiplespheres*>(templates[i]) &&
                dynamic_cast<FixTemplateMultiplespheres*>(templates[i])->is_bonded()
              )
              error->one(FLERR,"'bonded' and setting values for a fix property upon insertion can not be used together");
        }

        for(int i = 0; i < n; i++)
        {
            if (!pti_list[i]->fix_property)
            {
                pti_list[i]->fix_property = new FixPropertyAtom*[1];
                if (pti_list[i]->fix_property_value)
                    error->one(FLERR, "Internal error (fix property pti list)");
                pti_list[i]->fix_property_value = new double*[1];
                pti_list[i]->fix_property_value[0] = new double[1];
                if (pti_list[i]->fix_property_nentry)
                    error->one(FLERR, "Internal error (fix property pti list)");
                pti_list[i]->fix_property_nentry = new int[1];
            }
            pti_list[i]->fix_property[0] = fp;
            pti_list[i]->fix_property_value[0][0] = val;
            pti_list[i]->n_fix_property = 1;
            pti_list[i]->fix_property_nentry[0] = 1;
        }
    }

    for(int i = 0; i < n; i++)
        pti_list[i]->setFixTemplate(fix_template_);
}

/* ----------------------------------------------------------------------
   set particle properties - only pti needs to know which properties to set
   loop to n, not n_pti, since not all particles may have been inserted
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::insert(int n)
{
    int ninserted_spheres_local = 0;
    for(int i = 0; i < n; i++)
    {
        
        ninserted_spheres_local += pti_list[i]->insert();
    }
    return ninserted_spheres_local;
}

/* ----------------------------------------------------------------------
   wrap up insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::finalize_insertion()
{
    for(int i = 0; i < ntemplates; i++)
        templates[i]->finalize_insertion();
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::vol_expect()
{
    return volexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::mass_expect()
{
    return massexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::min_rad(int type)
{
    //get minrad
    double minrad_type = 1000.;
    for(int i = 0; i < ntemplates;i++)
    {
      if(
            (type >= templates[i]->mintype() && type <= templates[i]->maxtype()) &&
            (templates[i]->min_rad() < minrad_type)
        )
        minrad_type = templates[i]->min_rad();
    }

    return minrad_type;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::max_rad(int type)
{
    //get maxrad
    double maxrad_type = 0.;
    for(int i = 0; i < ntemplates;i++)
    {
      
      if(!templates[i]->use_rad_for_cut_neigh_and_ghost())
        continue;

      if(
          (type >= templates[i]->mintype() && type <= templates[i]->maxtype()) &&
          (templates[i]->max_rad() > maxrad_type)
        )
        maxrad_type = templates[i]->max_rad();
    }

    return maxrad_type;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_type()
{
    return maxtype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::min_type()
{
    return mintype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_nspheres()
{
    return maxnspheres;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]) + comm->me;

  random->reset(seed);
}

/* ----------------------------------------------------------------------
   generate a hash for unique identification
------------------------------------------------------------------------- */
unsigned int FixParticledistributionDiscrete::generate_hash()
{
    unsigned int hash = 0;
    unsigned int start = seed*420001; // it's magic
    add_hash_value(ntemplates, start, hash);
    for (int i=0; i<ntemplates; i++)
    {
        add_hash_value(distweight[i], start, hash);
        add_hash_value(distorder[i], start, hash);
        add_hash_value((int)templates[i]->generate_hash(), start, hash);
    }
    add_hash_value(maxtype, start, hash);
    add_hash_value(mintype, start, hash);
    add_hash_value(volexpect, start, hash);
    add_hash_value(massexpect, start, hash);
    add_hash_value(maxnspheres, start, hash);
    return hash;
}

void FixParticledistributionDiscrete::add_hash_value(const int value, unsigned int &start, unsigned int &hash)
{
    if (value >= 0)
        hash = hash*start + (unsigned int)value;
    else
        hash = hash*start + (unsigned int)(-value) + INT_MAX;
    start = start*seed;
}

void FixParticledistributionDiscrete::add_hash_value(double value, unsigned int &start, unsigned int &hash)
{
    int ivalue = 0;
    if (value < 0)
        value *= -1.0;
    if (value > 1e-50)
    {
        while (value > 1e6)
            value *= 1e-6;

        while (value < 1e0)
            value *= 1e6;
        int high = (int) value;
        double remainder = (value - (double)high)*1e6;
        int low = (int) remainder;
        ivalue = high + low;
    }
    add_hash_value(ivalue, start, hash);
}
