/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_MEMORY_NS_H
#define LMP_MEMORY_NS_H

#include "stdlib.h"
#include "lmptype.h"

using namespace LAMMPS_NS;

namespace LAMMPS_MEMORY_NS {

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *create(TYPE *&array, int n)
    {
      bigint nbytes = ((bigint) sizeof(TYPE)) * n;
      array = (TYPE *) malloc(nbytes);
      return array;
    }

/* ----------------------------------------------------------------------
   grow or shrink 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    TYPE *grow(TYPE *&array, int n)
    {
      if (array == NULL) return create(array,n);

      bigint nbytes = ((bigint) sizeof(TYPE)) * n;
      array = (TYPE *) realloc(array,nbytes);
      return array;
    }

/* ----------------------------------------------------------------------
   destroy a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
    void destroy(TYPE *array)
    {
      free(array);
    }

  /* ----------------------------------------------------------------------
     create a 2d array
  ------------------------------------------------------------------------- */

    template <typename TYPE>
      TYPE **create(TYPE **&array, int n1, int n2)
      {
        bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
        TYPE *data = (TYPE *) malloc(nbytes);
        nbytes = ((bigint) sizeof(TYPE *)) * n1;
        array = (TYPE **) malloc(nbytes);

        bigint n = 0;
        for (int i = 0; i < n1; i++) {
          array[i] = &data[n];
          n += n2;
        }
        return array;
      }

  /* ----------------------------------------------------------------------
     grow or shrink 1st dim of a 2d array
     last dim must stay the same
  ------------------------------------------------------------------------- */

    template <typename TYPE>
      TYPE **grow(TYPE **&array, int n1, int n2)
      {
        if (array == NULL) return create(array,n1,n2);

        bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
        TYPE *data = (TYPE *) realloc(array[0],nbytes);
        nbytes = ((bigint) sizeof(TYPE *)) * n1;
        array = (TYPE **) realloc(array,nbytes);

        bigint n = 0;
        for (int i = 0; i < n1; i++) {
          array[i] = &data[n];
          n += n2;
        }
        return array;
      }

  /* ----------------------------------------------------------------------
     destroy a 2d array
  ------------------------------------------------------------------------- */

    template <typename TYPE>
      void destroy(TYPE **array)
      {
        if (array == NULL) return;
        free(array[0]);
        free(array);
      }

    /* ----------------------------------------------------------------------
       create a 3d array
    ------------------------------------------------------------------------- */

      template <typename TYPE>
        TYPE ***create(TYPE ***&array, int n1, int n2, int n3)
        {
          bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
          TYPE *data = (TYPE *) malloc(nbytes);
          nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
          TYPE **plane = (TYPE **) malloc(nbytes);
          nbytes = ((bigint) sizeof(TYPE **)) * n1;
          array = (TYPE ***) malloc(nbytes);

          int i,j;
          bigint m;
          bigint n = 0;
          for (i = 0; i < n1; i++) {
            m = ((bigint) i) * n2;
            array[i] = &plane[m];
            for (j = 0; j < n2; j++) {
              plane[m+j] = &data[n];
              n += n3;
            }
          }
          return array;
        }

    /* ----------------------------------------------------------------------
       grow or shrink 1st dim of a 3d array
       last 2 dims must stay the same
    ------------------------------------------------------------------------- */

      template <typename TYPE>
        TYPE ***grow(TYPE ***&array, int n1, int n2, int n3)
        {
          if (array == NULL) return create(array,n1,n2,n3);

          bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
          TYPE *data = (TYPE *) realloc(array[0][0],nbytes);
          nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
          TYPE **plane = (TYPE **) realloc(array[0],nbytes);
          nbytes = ((bigint) sizeof(TYPE **)) * n1;
          array = (TYPE ***) realloc(array,nbytes);

          int i,j;
          bigint m;
          bigint n = 0;
          for (i = 0; i < n1; i++) {
            m = ((bigint) i) * n2;
            array[i] = &plane[m];
            for (j = 0; j < n2; j++) {
              plane[m+j] = &data[n];
              n += n3;
            }
          }
          return array;
        }

    /* ----------------------------------------------------------------------
       destroy a 3d array
    ------------------------------------------------------------------------- */

      template <typename TYPE>
        void destroy(TYPE ***array)
        {
          if (array == NULL) return;
          free(array[0][0]);
          free(array[0]);
          free(array);
        }
} /* LAMMPS_NS */

#endif /* MEMORY_H_ */
