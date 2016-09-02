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

#ifndef LMP_MULTISPHERE_PARALLEL_I_H
#define LMP_MULTISPHERE_PARALLEL_I_H

/* ---------------------------------------------------------------------- */

inline int MultisphereParallel::pack_exchange_rigid(int i, double *buf)
{
  
  int m = 1; 
  double xbound[3];
  bool dummy = false;

  // calculate xbound in global coo sys
  MathExtraLiggghts::local_coosys_to_cartesian
  (
    xbound,xcm_to_xbound_(i),
    ex_space_(i),ey_space_(i),ez_space_(i)
  );

  vectorAdd3D(xcm_(i),xbound,xbound);

  // have to pack xbound first because exchange() tests against first 3 values in buffer
  vectorToBuf3D(xbound,buf,m);
  
  m += customValues_.pushElemToBuffer(i,&(buf[m]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);

  buf[0] = m;
  
  return m;
}

/* ---------------------------------------------------------------------- */

inline int MultisphereParallel::unpack_exchange_rigid(double *buf)
{
  double xbound[3];
  bool dummy = false;
  int m = 0;

  int nvalues = buf[m++];

  bufToVector3D(xbound,buf,m);
  m += customValues_.popElemFromBuffer(&(buf[m]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);

  nbody_++;

  return nvalues;
}

#endif
