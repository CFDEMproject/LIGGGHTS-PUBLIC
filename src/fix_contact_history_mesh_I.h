/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_CONTACT_HISTORY_MESH_I_H
#define LMP_CONTACT_HISTORY_MESH_I_H

  /* ---------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::handleContact(int iP, int idTri, double *&history)
  {
    
    // check if contact with iTri was there before
    // if so, set history to correct location and return
    if(haveContact(iP,idTri,history))
      return true;

    // else new contact - add contact if did not calculate contact with coplanar neighbor already
    
    if(coplanarContactAlready(iP,idTri))
        // did not add new contact
        return false;
    else
    {
        addNewTriContactToExistingParticle(iP,idTri,history);

        // check if one of the contacts of previous steps is coplanar with iTri
        
        // if so, copy history
        // also check if this contact has delflag = false, i.e. has been executed already
        // this step. If so, signalize not to execute this contact (return false)
        checkCoplanarContactHistory(iP,idTri,history);
        return true;
    }
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::swap(int ilocal,int ineigh, int jneigh, bool delflag_swap)
  {
      
      int id_temp;
      bool delflag_temp;

      id_temp                  = partner_[ilocal][ineigh];
      partner_[ilocal][ineigh] = partner_[ilocal][jneigh];
      partner_[ilocal][jneigh] = id_temp;

      vectorCopyN(&(contacthistory_[ilocal][ineigh*dnum_]),swap_,                                   dnum_);
      vectorCopyN(&(contacthistory_[ilocal][jneigh*dnum_]),&(contacthistory_[ilocal][ineigh*dnum_]),dnum_);
      vectorCopyN(swap_,                                   &(contacthistory_[ilocal][jneigh*dnum_]),dnum_);

      if(delflag_swap)
      {
          delflag_temp             = delflag_[ilocal][ineigh];
          delflag_[ilocal][ineigh] = delflag_[ilocal][jneigh];
          delflag_[ilocal][jneigh] = delflag_temp;
      }
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::haveContact(int iP, int idTri, double *&history)
  {
    int *tri = partner_[iP];
    double *nn = fix_nneighs_->vector_atom;
    int nneighs = static_cast<int>(round(nn[iP]));

    for(int i = 0; i < nneighs; i++)
    {
        if(tri[i] == idTri)
        {
            if(dnum_ > 0) history = &(contacthistory_[iP][i*dnum_]);
            delflag_[iP][i] = false;
            return true;
        }
    }
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::coplanarContactAlready(int iP, int idTri)
  {
    double *nn = fix_nneighs_->vector_atom;
    int nneighs = static_cast<int>(round(nn[iP]));
    for(int i = 0; i < nneighs; i++)
    {
      
      int idPartnerTri = partner_[iP][i];

      if(idPartnerTri >= 0 && idPartnerTri != idTri && mesh_->map(idPartnerTri) >= 0 && mesh_->areCoplanarNodeNeighs(idPartnerTri,idTri))
      {
        
        // other coplanar contact handled already - do not handle this contact
        if(!delflag_[iP][i]) return true;
      }
    }

    // no coplanar contact found - handle this contact
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::checkCoplanarContactHistory(int iP, int idTri, double *&history)
  {
    int *tri = partner_[iP];
    double *nn = fix_nneighs_->vector_atom;
    int nneighs = static_cast<int>(round(nn[iP]));

    for(int i = 0; i < nneighs; i++)
    {
      
      if(tri[i] >= 0 && tri[i] != idTri && mesh_->map(tri[i]) >= 0 && mesh_->areCoplanarNodeNeighs(tri[i],idTri))
      {
          
          // copy contact history
          if(dnum_ > 0) vectorCopyN(&(contacthistory_[iP][i*dnum_]),history,dnum_);
          
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::addNewTriContactToExistingParticle(int iP, int idTri, double *&history)
  {
      
      int nneighs = static_cast<int>(round(fix_nneighs_->vector_atom[iP]));
      int iContact = -1;

      if(-1 == idTri)
        error->one(FLERR,"internal error");

      if(npartner_[iP] >= nneighs)
      {
        
        error->one(FLERR,"internal error");
      }

      for(int ineigh = 0; ineigh < nneighs; ineigh++)
      {
          if(-1 == partner_[iP][ineigh])
          {
              iContact = ineigh;
              break;
          }
      }

      if(iContact >= nneighs)
        error->one(FLERR,"internal error");

      partner_[iP][iContact] = idTri;
      delflag_[iP][iContact] = false;

      if(dnum_ > 0)
      {
          history = &(contacthistory_[iP][iContact*dnum_]);
          vectorZeroizeN(history,dnum_);
      }
      else
          history = 0;

      npartner_[iP]++;

  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistoryMesh::n_contacts()
  {
    int ncontacts = 0, nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
           ncontacts += npartner_[i];
    return ncontacts;
  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistoryMesh::n_contacts(int contact_groupbit)
  {
    int ncontacts = 0, nlocal = atom->nlocal;
    int *mask = atom->mask;

    for(int i = 0; i < nlocal; i++)
        if(mask[i] & contact_groupbit)
           ncontacts += npartner_[i];
    return ncontacts;
  }

#endif
