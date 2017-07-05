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

    Andreas Aigner (DCS Computing GmbH)

    Copyright 2014-2015 JKU Linz
    Copyright 2015-     DCS Computing GmbH
------------------------------------------------------------------------- */

#ifndef REGION_NEIGHBOR_LIST_BASE_H
#define REGION_NEIGHBOR_LIST_BASE_H

namespace LAMMPS_NS {

/**
 * @brief Interface class for \a RegionNeighborList
 *
 * Allows to interacte with \a RegionNeighborList via a non-template interface.
 *
 * This interface is further used for \a IRegionNeighborFieldList. The design
 * requires to reimplement / redirect all here defined funtions in
 * \a RegionNeighborFieldList. Therefore, keep in small and clean!
 */

class IRegionNeighborList
{
public:
    virtual ~IRegionNeighborList() {}

    virtual int mbins() const = 0;
    virtual void getDimensions(int *dims) const = 0;
    virtual void getLocalDimensions(int *dims) const = 0;
    virtual void getBinSize(double *binsize) const = 0;
    virtual void getOrigin(double *origin) const = 0;

    virtual int getSizeOne() const = 0;
};

}

#endif // REGION_NEIGHBOR_LIST_BASE_H

