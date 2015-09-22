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

    Richard Berger (JKU Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "comm.h"
#include "error.h"

#include "pair_gran_proxy.h"
#include "granular_pair_style.h"

using namespace LAMMPS_NS;
using namespace LIGGGHTS::PairStyles;

PairGranProxy::PairGranProxy(LAMMPS * lmp) : PairGran(lmp), impl(NULL)
{
}

PairGranProxy::~PairGranProxy()
{
  delete impl;
}

void PairGranProxy::settings(int nargs, char ** args)
{
  delete impl;
  int64_t variant = Factory::instance().selectVariant("gran", nargs, args);
  impl = Factory::instance().create("gran", variant, lmp, this);

  if(impl) {
    impl->settings(nargs, args);
  } else {
    
    error->one(FLERR, "unknown contact model");
  }
}

void PairGranProxy::init_granular()
{
  impl->init_granular();
}

void PairGranProxy::write_restart_settings(FILE * fp)
{
  impl->write_restart_settings(fp);
}

void PairGranProxy::read_restart_settings(FILE * fp)
{
  int me = comm->me;

  int64_t selected = -1;
  if(me == 0){
    // read model hashcode, but reset file pointer afterwards.
    // this way read_restart_settings can still read the hashcode (sanity check)
    fread(&selected, sizeof(int64_t), 1, fp);
    fseek(fp, -sizeof(int64_t), SEEK_CUR);
  }
  MPI_Bcast(&selected,8,MPI_CHAR,0,world);

  impl = Factory::instance().create("gran", selected, lmp, this);

  if(impl) {
    impl->read_restart_settings(fp);
  } else {
    error->one(FLERR, "unknown contact model");
  }
}

void PairGranProxy::compute_force(int eflag, int vflag, int addflag)
{
  impl->compute_force(this, eflag, vflag, addflag);
}

double PairGranProxy::stressStrainExponent()
{
  return impl->stressStrainExponent();
}

int64_t PairGranProxy::hashcode() {
  return impl->hashcode();
}
