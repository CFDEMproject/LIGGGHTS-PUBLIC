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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ERROR/WARNING messages:

E: Bond atoms %d %d missing on proc %d at step %ld

The 2nd atom needed to compute a particular bond is missing on this
processor.  Typically this is because the pairwise cutoff is set too
short or the bond has blown apart and an atom is too far away.

E: Bond extent > half of periodic box length

This error was detected by the neigh_modify check yes setting.  It is
an error because the bond atoms are so far apart it is ambiguous how
it should be defined.

E: Angle atoms %d %d %d missing on proc %d at step %ld

One or more of 3 atoms needed to compute a particular angle are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the angle has blown apart and an atom is
too far away.

E: Angle extent > half of periodic box length

This error was detected by the neigh_modify check yes setting.  It is
an error because the angle atoms are so far apart it is ambiguous how
it should be defined.

E: Dihedral atoms %d %d %d %d missing on proc %d at step %ld

One or more of 4 atoms needed to compute a particular dihedral are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the dihedral has blown apart and an atom is
too far away.

E: Dihedral/improper extent > half of periodic box length

This error was detected by the neigh_modify check yes setting.  It is
an error because the dihedral atoms are so far apart it is ambiguous
how it should be defined.

E: Improper atoms %d %d %d %d missing on proc %d at step %ld

One or more of 4 atoms needed to compute a particular improper are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the improper has blown apart and an atom is
too far away.

*/
