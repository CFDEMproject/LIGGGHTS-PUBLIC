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

#ifndef LMP_FIX_MULTISPHERE_COMM_I_H_
#define LMP_FIX_MULTISPHERE_COMM_I_H_

      // restart
      int pack_restart(int i, double *buf);
      void unpack_restart(int nlocal, int nth);

      // communication
      int pack_exchange(int i, double *buf);
      int unpack_exchange(int nlocal, double *buf);

      void forward_comm();
      void reverse_comm();

      int pack_comm(int, int*, double*, int, int*);
      int pack_comm_body(int, int*, double*, int, int*);
      int pack_comm_image_displace(int, int*, double*, int, int*);
      int pack_comm_v_omega(int, int*, double*, int, int*);
      int pack_comm_f_torque(int, int*, double*, int, int*);
      int pack_comm_temp(int, int*, double*, int, int*);

      void unpack_comm(int, int, double*);
      void unpack_comm_body(int, int, double*);
      void unpack_comm_image_displace(int, int, double*);
      void unpack_comm_v_omega(int, int, double*);
      void unpack_comm_f_torque(int, int, double*);
      void unpack_comm_temp(int, int, double*);

      int pack_reverse_comm(int, int, double*);
      int pack_reverse_comm_x_v_omega(int, int, double*);
      int pack_reverse_comm_v_omega(int, int, double*);
      int pack_reverse_comm_image(int n, int first, double *buf);
      int pack_reverse_comm_displace(int n, int first, double *buf);
      int pack_reverse_comm_temp(int n, int first, double *buf);
      void unpack_reverse_comm(int, int*, double*);
      void unpack_reverse_comm_x_v_omega(int, int*, double*);
      void unpack_reverse_comm_v_omega(int, int*, double*);
      void unpack_reverse_comm_image(int n, int *list, double *buf);
      void unpack_reverse_comm_displace(int n, int *list, double *buf);
      void unpack_reverse_comm_temp(int n, int *list, double *buf);

#endif
