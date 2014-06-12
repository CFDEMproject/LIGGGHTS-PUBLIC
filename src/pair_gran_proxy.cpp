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
