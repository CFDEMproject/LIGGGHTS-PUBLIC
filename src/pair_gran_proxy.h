#ifdef PAIR_CLASS
PairStyle(gran,PairGranProxy)
#else

#ifndef PAIR_GRAN_PROXY_H
#define PAIR_GRAN_PROXY_H

#include "pair_gran.h"
#include "granular_pair_style.h"

namespace LAMMPS_NS
{
class PairGranProxy : public PairGran
{
  LIGGGHTS::PairStyles::IGranularPairStyle * impl;

public:
  PairGranProxy(LAMMPS * lmp);
  virtual ~PairGranProxy();

  virtual void settings(int nargs, char ** args);
  virtual void init_granular();
  virtual void write_restart_settings(FILE * fp);
  virtual void read_restart_settings(FILE * fp);
  virtual void compute_force(int eflag, int vflag, int addflag);

  virtual double stressStrainExponent();
  virtual int64_t hashcode();
};
}

#endif // PAIR_GRAN_PROXY_H

#endif
