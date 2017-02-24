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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_HISTOGRAM_H
#define LMP_HISTOGRAM_H

#include "fix.h"
#include "modify.h"
#include "container.h"
#include "fix_property_atom.h"
#include "atom.h"
#include <vector>

namespace LAMMPS_NS {

template<typename T>
class Histogram : protected Pointers
{
  public:

    Histogram(LAMMPS *lmp,int _type,const char *_id, T _min_val,T _max_val,int _n_bins);
    ~Histogram();
    inline void add_fix(FixPropertyAtom *_fix_histogram_counter);

    inline void reset_counters();

    inline void add_atom_value(int iatom,T value);
    inline void add_global_value(T value);

    inline void tally_global_values();
    inline double get_global_value(int ibin);

    inline int n_bins()
    { return n_bins_; }

    inline const char *id()
    { return id_; }

    inline ScalarContainer<int> & global_count_this()
    { return global_count_this_; }

    inline FixPropertyAtom* fix_histogram_counter()
    { return fix_histogram_counter_; }

    inline int size_restart();
    inline int write_restart(double *buf);
    inline int read_restart(double *buf);

  private:

    int type_;
    char *id_;

    T min_val_,max_val_;
    int n_bins_;
    T inv_distance;

    FixPropertyAtom *fix_histogram_counter_;
    ScalarContainer<int> global_count_this_;
    ScalarContainer<int> global_count_;
};

/* ---------------------------------------------------------------------- */

template<typename T>
Histogram<T>::Histogram(LAMMPS *lmp, int _type,const char *_id,
                        T _min_val,T _max_val,int _n_bins)
    : Pointers(lmp),
    type_(_type),
    id_(0),
    min_val_(_min_val),
    max_val_(_max_val),
    n_bins_(_n_bins),
    fix_histogram_counter_(0)
{
    id_ = new char[strlen(_id)+1];
    strcpy(id_,_id);

    // hack
    if(sizeof(double) == sizeof(T))
    {
        double distance = (static_cast<double>(max_val_) - static_cast<double>(min_val_)) / static_cast<double>(n_bins_);
        inv_distance = 1./distance;
    }

    reset_counters();
}

/* ---------------------------------------------------------------------- */

template<typename T>
Histogram<T>::~Histogram()
{
    delete []id_;
}

/* ---------------------------------------------------------------------- */

template<typename T>
void Histogram<T>::add_fix(FixPropertyAtom *_fix_histogram_counter)
{
    fix_histogram_counter_ = _fix_histogram_counter;
}

/* ---------------------------------------------------------------------- */

template<typename T>
void Histogram<T>::reset_counters()
{
    global_count_.clearContainer();
    global_count_this_.clearContainer();

    for(int i = 0; i < n_bins_; i++)
    {
        global_count_.add(0);
        global_count_this_.add(0);
    }

    if(fix_histogram_counter_)
    {
        int nall = atom->nlocal + atom->nghost;
        vectorZeroizeN(&(fix_histogram_counter_->array_atom[0][0]),nall*n_bins_);
    }
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void Histogram<T>::add_atom_value(int iatom,T value)
{
    error->one(FLERR,"Illegal call");
}

/* ---------------------------------------------------------------------- */

template<>
inline void Histogram<double>::add_atom_value(int iatom,double value)
{
    int ibin = static_cast<int>((value - min_val_) * inv_distance);
    
    if(ibin > -1 && ibin < n_bins_)
        fix_histogram_counter_->array_atom[iatom][ibin]++;
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void Histogram<T>::add_global_value(T value)
{
    error->one(FLERR,"Illegal call, need template specialization");
}

/* ---------------------------------------------------------------------- */

template<>
inline void Histogram<double>::add_global_value(double value)
{
    int ibin = static_cast<int>((value - min_val_) * inv_distance);
    
    if(ibin > -1 && ibin < n_bins_)
        global_count_this_(ibin)++;
}

/* ---------------------------------------------------------------------- */

template<typename T>
void Histogram<T>::tally_global_values()
{
    for(int i = 0; i < n_bins_; i++)
    {
        global_count_(i) += global_count_this_(i);
        global_count_this_(i) = 0;
    }
}

/* ---------------------------------------------------------------------- */

template<typename T>
double Histogram<T>::get_global_value(int ibin)
{
    
    if(ibin < n_bins_)
        return static_cast<double>(global_count_(ibin));

    return 0.;
}

/* ---------------------------------------------------------------------- */

template<typename T>
int Histogram<T>::size_restart()
{
    int size = 4;
    size += global_count_.size();

    return size;
}

/* ---------------------------------------------------------------------- */

template<typename T>
int Histogram<T>::write_restart(double *buf)
{
    int m = 0;

    buf[m++] = static_cast<double>(type_);
    buf[m++] = static_cast<double>(min_val_);
    buf[m++] = static_cast<double>(max_val_);
    buf[m++] = static_cast<double>(n_bins_);
    m += global_count_.pushToBuffer_plain(&(buf[m]));

    return m;
}

/* ---------------------------------------------------------------------- */

template<typename T>
int Histogram<T>::read_restart(double *buf)
{
    int m = 0;

    int type_re = static_cast<int>(buf[m++]);
    if(type_re != type_)
        error->one(FLERR,"Histogram re-started with different properties");

    T min_val_re = static_cast<T>(buf[m++]);
    if(min_val_re != min_val_)
        error->one(FLERR,"Histogram re-started with different properties");

    T max_val_re = static_cast<T>(buf[m++]);
    if(max_val_re != max_val_)
        error->one(FLERR,"Histogram re-started with different properties");

    int n_bins_re = static_cast<int>(buf[m++]);
    if(n_bins_re != n_bins_)
        error->one(FLERR,"Histogram re-started with different properties");

    m += global_count_.pullFromBuffer_plain(&(buf[m]));

    return m;
}

};

#endif
