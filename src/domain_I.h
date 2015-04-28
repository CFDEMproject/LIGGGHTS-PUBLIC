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

#ifndef LMP_DOMAIN_I_H
#define LMP_DOMAIN_I_H

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   check if coordinate in domain, subdomain or extended subdomain
   need to test with <= and >= for domain, and < and >= for subdomain
   inlined for performance
------------------------------------------------------------------------- */

inline int Domain::is_in_domain(double* pos) 
{
    if(is_wedge)
        return is_in_domain_wedge(pos);

    if
    (
        pos[0] >= boxlo[0] && pos[0] <= boxhi[0] &&
        pos[1] >= boxlo[1] && pos[1] <= boxhi[1] &&
        pos[2] >= boxlo[2] && pos[2] <= boxhi[2]
    )   return 1;
    return 0;
}

inline int Domain::is_in_subdomain(double* pos) 
{
    if(is_wedge)
        return is_in_subdomain_wedge(pos);

    double checkhi[3];

    // need to test with and < and >= for subdomain
    // but need to ensure elements right at boxhi
    // are tested correctly
    if(subhi[0] == boxhi[0])
        checkhi[0] = boxhi[0] + SMALL_DMBRDR;
    else
        checkhi[0] = subhi[0];

    if(subhi[1] == boxhi[1])
        checkhi[1] = boxhi[1] + SMALL_DMBRDR;
    else
        checkhi[1] = subhi[1];

    if(subhi[2] == boxhi[2])
        checkhi[2] = boxhi[2] + SMALL_DMBRDR;
    else
        checkhi[2] = subhi[2];

    if ( pos[0] >= sublo[0] && pos[0] < checkhi[0] &&
         pos[1] >= sublo[1] && pos[1] < checkhi[1] &&
         pos[2] >= sublo[2] && pos[2] < checkhi[2])
        return 1;
    return 0;
}

inline int Domain::is_in_extended_subdomain(double* pos) 
{
    if(is_wedge)
        return is_in_extended_subdomain_wedge(pos);

    // called on insertion
    // yields true if particle would be in subdomain after box extension
    
    if (is_in_subdomain(pos))
        return 1;
    else if (dimension == 2)
        error->all(FLERR,"Domain::is_in_extended_subdomain() not implemented for 2d");
    else 
    {
        bool flag = true;
        for(int idim = 0; idim < 3; idim++)
        {
            
            if (comm->procgrid[idim] == 1) {}
            else if(comm->myloc[idim] == comm->procgrid[idim]-1)
                flag = flag && (pos[idim] >= sublo[idim]);
            else if(comm->myloc[idim] == 0)
                flag = flag && (pos[idim] <= subhi[idim]);
            
            else
                flag = flag && (pos[idim] >= sublo[idim] && pos[idim] < subhi[idim]);
        }
        if(flag) return 1;
        return 0;
    }
    return 0;
}

/* ----------------------------------------------------------------------
   check distance from borders of subbox
------------------------------------------------------------------------- */

inline double Domain::dist_subbox_borders(double* pos) 
{
    if(is_wedge)
        return dist_subbox_borders_wedge(pos);

    double deltalo[3], deltahi[3], dlo, dhi;

    vectorSubtract3D(sublo,pos,deltalo);
    vectorSubtract3D(subhi,pos,deltahi);
    vectorAbs3D(deltalo);
    vectorAbs3D(deltahi);

    dlo = vectorMin3D(deltalo);
    dhi = vectorMin3D(deltahi);

    if(dlo < dhi)
        return dlo;
    return dhi;
}

/* ----------------------------------------------------------------------
   return smallest extent ob subbox
------------------------------------------------------------------------- */

inline void Domain::min_subbox_extent(double &min_extent,int &dim) 
{
    if(is_wedge)
        error->one(FLERR,"missing implementation");

    double delta[3];
    vectorSubtract3D(subhi,sublo,delta);
    
    min_extent = vectorMin3D(delta,dim);
}

/* ----------------------------------------------------------------------
   domain check - not used very often, so not inlined
------------------------------------------------------------------------- */

inline int Domain::is_periodic_ghost(int i) 
{
    if(is_wedge)
        return is_periodic_ghost_wedge(i);

    int idim;
    int nlocal = atom->nlocal;
    double *x = atom->x[i];
    double halfskin = 0.5*neighbor->skin;

    if(i < nlocal) return 0;
    else
    {
        for(idim = 0; idim < 3; idim++)
             if ((x[idim] < (boxlo[idim]+halfskin) || x[idim] > (boxhi[idim]-halfskin)) && periodicity[idim])
             
                return 1;
    }
    return 0;
}

/* ----------------------------------------------------------------------
   check if atom is the true representation of the particle on this subdomain
   used when tallying stats across owned and ghost particles
------------------------------------------------------------------------- */

inline bool Domain::is_owned_or_first_ghost(int i) 
{
    if(!atom->tag_enable)
        error->one(FLERR,"The current simulation setup requires atoms to have tags");
    if(0 == atom->map_style)
        error->one(FLERR,"The current simulation setup requires an 'atom_modify map' command to allocate an atom map");

    if(i == atom->map(atom->tag[i]))
        return true;
    return false;
}

#endif
