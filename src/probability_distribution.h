/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_PROBABILITY_DISTRIBUTION_H
#define LMP_PROBABILITY_DISTRIBUTION_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "random_park.h"
#include "error.h"
#include "pointers.h"

enum{RANDOM_CONSTANT,RANDOM_UNIFORM,RANDOM_GAUSSIAN,RANDOM_LOGNORMAL};

namespace LMP_PROBABILITY_NS {
  class PDF
  {
    public:
      PDF(LAMMPS_NS::Error *error)
      {
          mu_ = sigma_ = min_ = max_ = 0.;
          h1_ = h2_ = 0.;
          mass_shift_ = 0;
          this->error = error;
      }
      ~PDF(){}

      int rand_style_;

      double mu_,sigma_;
      double min_,max_;

      // helper
      double h1_,h2_;

      // if 1, pdf is shifted from number to mass based pdf
      int mass_shift_;

      LAMMPS_NS::Error *error;

      inline int rand_style()
      { return rand_style_; }

      inline void set_min_max(double min,double max)
      {
          min_ = min;
          max_ = max;
      }

      inline void activate_mass_shift()
      { mass_shift_ = 1; }

      template<int RAND_STYLE> void set_params(double)
      { error->all(FLERR,"Faulty usage of Probability::set_params"); }

      template<int RAND_STYLE> void set_params(double,double)
      { error->all(FLERR,"Faulty usage of Probability::set_params"); }
  };

  inline double pdf_max(PDF *pdf)
  {
      return pdf->max_;
  }

  inline double pdf_min(PDF *pdf)
  {
      return pdf->min_;
  }

  template <int RAND_STYLE> inline double expectancy_value(PDF *pdf)
  {
      pdf->error->all(FLERR,"Faulty usage of Probability::expectancy");
      return 0.;
  }

  template <int RAND_STYLE> inline double cubic_expectancy_value(PDF *pdf)
  {
      pdf->error->all(FLERR,"Faulty usage of Probability::volume_expectancy");
      return 0.;
  }

  template <int RAND_STYLE> inline double rand_value(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
      pdf->error->all(FLERR,"Faulty usage of Probability::rand");
      return 0.;
  }

  //------------------------------------------------------------------------------
  // CONSTANT
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_CONSTANT>(double val)
  {
         rand_style_ = RANDOM_CONSTANT;
         mu_ = val;
         set_min_max(mu_,mu_);
  }

  template<> inline double cubic_expectancy_value<RANDOM_CONSTANT>(PDF *pdf)
  {

     return pdf->mu_*pdf->mu_*pdf->mu_;
  }

  template<> inline double expectancy_value<RANDOM_CONSTANT>(PDF *pdf)
  {

     return pdf->mu_;
  }

  template<> inline double rand_value<RANDOM_CONSTANT>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     return pdf->mu_;
  }

  //------------------------------------------------------------------------------
  // UNIFORM
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_UNIFORM>(double min, double max)
  {
         rand_style_ = RANDOM_UNIFORM;
         set_min_max(min,max);
         h1_ = 2./(1./(min_*min_)-1./(max_*max_));
         h2_ = h1_/(2.*min_*min_);
  }

  template<> inline double cubic_expectancy_value<RANDOM_UNIFORM>(PDF *pdf)
  {
     if(pdf->mass_shift_)
        pdf->error->all(FLERR,"mass distribution not implemented for uniform");

     return 0.25*(pdf->max_*pdf->max_*pdf->max_+
                  pdf->max_*pdf->max_*pdf->min_+
                  pdf->max_*pdf->min_*pdf->min_+
                  pdf->min_*pdf->min_*pdf->min_);
  }

  template<> inline double expectancy_value<RANDOM_UNIFORM>(PDF *pdf)
  {
     if(!pdf->mass_shift_)
        return 0.5 * (pdf->min_ + pdf->max_);
     else
        return sqrt(pdf->h1_/(2.*(pdf->h2_-0.5)));
  }

  template<> inline double rand_value<RANDOM_UNIFORM>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     double rn =  rp->uniform();
     if(!pdf->mass_shift_)
        return (pdf->min_) + rn * (pdf->max_ - pdf->min_);
     else
        return sqrt(pdf->h1_/(2.*(pdf->h2_-rn)));
  }

  //------------------------------------------------------------------------------
  // GAUSSIAN
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_GAUSSIAN>(double mu, double sigma)
  {
         rand_style_ = RANDOM_GAUSSIAN;
         mu_ = mu;
         sigma_ = sigma;

         // set min-max to +- 3 sigma (99.73% of all values)
         set_min_max(mu_-3.*sigma_, mu_+3.*sigma_);

         if(min_ < 0.)
            error->all(FLERR,"Probablity distribution: mu-3*sigma < 0, please increase mu or decrease sigma");
  }

  template<> inline double cubic_expectancy_value<RANDOM_GAUSSIAN>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for gaussian");
     return pdf->mu_*(pdf->mu_*pdf->mu_+3*pdf->sigma_*pdf->sigma_);
  }

  template<> inline double expectancy_value<RANDOM_GAUSSIAN>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for gaussian");
     return pdf->mu_;
  }

  template<> inline double rand_value<RANDOM_GAUSSIAN>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for gaussian");
     double value;
     do
     {
         value = pdf->mu_ + rp->gaussian() * pdf->sigma_;
     } while (value < pdf->min_ || value > pdf->max_);
     return value;
  }

  //------------------------------------------------------------------------------
  // LOGNORMAL
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_LOGNORMAL>(double mu, double sigma)
  {
      error->all(FLERR,"lognormal distribution currently deactivated");

      rand_style_ = RANDOM_LOGNORMAL;
      mu_ = mu;
      sigma_ = sigma;

      // also here, take +- 3 sigma as min/max
      // change in expectancy considered negligable
      double min =  exp(mu_ - 3. * sigma_);
      double max =  exp(mu_ + 3. * sigma_);
      set_min_max(min, max);
      if(min_ < 0.)
            error->all(FLERR,"Probablity distribution: exp(mu-3*sigma) < 0, please increase mu or decrease sigma");
  }

  template<> inline double cubic_expectancy_value<RANDOM_LOGNORMAL>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for lognormal");
     return exp(3.*pdf->mu_+4.5*pdf->sigma_*pdf->sigma_);
  }

  template<> inline double expectancy_value<RANDOM_LOGNORMAL>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for lognormal");
     return exp(pdf->mu_ + 0.5 * pdf->sigma_ * pdf->sigma_);
  }

  template<> inline double rand_value<RANDOM_LOGNORMAL>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     if(pdf->mass_shift_) pdf->error->all(FLERR,"mass distribution not implemented for lognormal");
     double value;
     do
     {
         value = exp(pdf->mu_ + rp->gaussian() * pdf->sigma_);
     } while (value < pdf->min_ || value > pdf->max_);
     return value;
  }

  //------------------------------------------------------------------------------
  // MASTER FUNCTIONS
  //------------------------------------------------------------------------------

  inline double expectancy(PDF *pdf)
  {
      if(pdf->rand_style_ == RANDOM_CONSTANT) return expectancy_value<RANDOM_CONSTANT>(pdf);
      else if(pdf->rand_style_ == RANDOM_UNIFORM) return expectancy_value<RANDOM_UNIFORM>(pdf);
      else if(pdf->rand_style_ == RANDOM_GAUSSIAN) return expectancy_value<RANDOM_GAUSSIAN>(pdf);
      else if(pdf->rand_style_ == RANDOM_LOGNORMAL) return expectancy_value<RANDOM_LOGNORMAL>(pdf);
      else pdf->error->all(FLERR,"Faulty implemantation in Probability::expectancy");
      return 0.;
  }

  inline double cubic_expectancy(PDF *pdf)
  {
      if(pdf->rand_style_ == RANDOM_CONSTANT) return cubic_expectancy_value<RANDOM_CONSTANT>(pdf);
      else if(pdf->rand_style_ == RANDOM_UNIFORM) return cubic_expectancy_value<RANDOM_UNIFORM>(pdf);
      else if(pdf->rand_style_ == RANDOM_GAUSSIAN) return cubic_expectancy_value<RANDOM_GAUSSIAN>(pdf);
      else if(pdf->rand_style_ == RANDOM_LOGNORMAL) return cubic_expectancy_value<RANDOM_LOGNORMAL>(pdf);
      else pdf->error->all(FLERR,"Faulty implemantation in Probability::expectancy");
      return 0.;
  }

  inline double rand(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
      if(pdf->rand_style_ == RANDOM_CONSTANT) return rand_value<RANDOM_CONSTANT>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_UNIFORM) return rand_value<RANDOM_UNIFORM>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_GAUSSIAN) return rand_value<RANDOM_GAUSSIAN>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_LOGNORMAL) return rand_value<RANDOM_LOGNORMAL>(pdf,rp);
      else pdf->error->all(FLERR,"Faulty implemantation in Probability::rand");
      return 0.;
  }

};

#endif
