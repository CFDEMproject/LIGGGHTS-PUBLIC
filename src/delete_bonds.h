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

#ifdef COMMAND_CLASS

CommandStyle(delete_bonds,DeleteBonds)

#else

#ifndef LMP_DELETE_BONDS_H
#define LMP_DELETE_BONDS_H

#include "pointers.h"

namespace LAMMPS_NS {

class DeleteBonds : protected Pointers {
 public:
  DeleteBonds(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Delete_bonds command before simulation box is defined

The delete_bonds command cannot be used before a read_data,
read_restart, or create_box command.

E: Delete_bonds command with no atoms existing

No atoms are yet defined so the delete_bonds command cannot be used.

E: Cannot use delete_bonds with non-molecular system

Your choice of atom style does not have bonds.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot find delete_bonds group ID

Group ID used in the delete_bonds command does not exist.

E: Bond atom missing in delete_bonds

The delete_bonds command cannot find one or more atoms in a particular
bond on a particular processor.  The pairwise cutoff is too short or
the atoms are too far apart to make a valid bond.

E: Angle atom missing in delete_bonds

The delete_bonds command cannot find one or more atoms in a particular
angle on a particular processor.  The pairwise cutoff is too short or
the atoms are too far apart to make a valid angle.

E: Dihedral atom missing in delete_bonds

The delete_bonds command cannot find one or more atoms in a particular
dihedral on a particular processor.  The pairwise cutoff is too short
or the atoms are too far apart to make a valid dihedral.

E: Improper atom missing in delete_bonds

The delete_bonds command cannot find one or more atoms in a particular
improper on a particular processor.  The pairwise cutoff is too short
or the atoms are too far apart to make a valid improper.

*/
