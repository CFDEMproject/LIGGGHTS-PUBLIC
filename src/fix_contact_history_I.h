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

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#ifndef LMP_CONTACT_HISTORY_I_H
#define LMP_CONTACT_HISTORY_I_H

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistory::handleContact(int iP, int idTri, double *&history)
  {
    // check if contact with iTri was there before
    // if so, set history to correct location and return
    if(haveContact(iP,idTri,history))
      return true;

    // else new contact - add contact if did not calculate contact with coplanar neighbor already
    
    if(!coplanarContactAlready(iP,idTri))
    {
        addNewTriContactToExistingParticle(iP,idTri,history);

        // check if one of the contacts of previous steps is coplanar with iTri
        
        // if so, copy history
        // also check if this contact has delflag = false, i.e. has been executed already
        // this step. If so, signalize not to execute this contact (return false)
        checkCoplanarContactHistory(iP,idTri,history);
        return true;
    }

    // did not add new contact
    return false;
  }

  /* ----------------------------------------------------------------------
     mark all contacts for deletion
  ------------------------------------------------------------------------- */

  inline void FixContactHistory::markAllContacts()
  {
      int nlocal = atom->nlocal;

      for(int i = 0; i < nlocal; i++)
          for(int j = 0; j < npartner[i]; j++)
              delflag[i][j] = true;
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistory::haveContact(int iP, int idTri, double *&history)
  {
    int *tri = partner[iP];

    for(int i = 0; i < npartner[iP]; i++)
    {
        if(tri[i] == idTri)
        {
            if(dnum > 0) history = contacthistory[iP][i];
            delflag[iP][i] = false;
            return true;
        }
    }
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistory::coplanarContactAlready(int iP, int idTri)
  {
    int *tri = partner[iP];
    for(int i = 0; i < npartner[iP]; i++)
    {
      
      if(tri[i] != idTri && mesh_->map(tri[i]) >= 0 && mesh_->areCoplanarNodeNeighs(tri[i],idTri))
      {
        
        // other coplanar contact handled already - do not handle this contact
        if(!delflag[iP][i]) return true;
      }
    }

    // no coplanar contact found - handle this contact
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistory::checkCoplanarContactHistory(int iP, int idTri, double *&history)
  {
    int *tri = partner[iP];
    for(int i = 0; i < npartner[iP]; i++)
    {
      
      if(tri[i] != idTri && mesh_->map(tri[i]) >= 0 && mesh_->areCoplanarNodeNeighs(tri[i],idTri))
      {
          
          // copy contact history
          if(dnum > 0) vectorCopyN(contacthistory[iP][i],history,dnum);
          
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistory::addNewTriContactToExistingParticle(int iP, int idTri, double *&history)
  {
      int numCont = npartner[iP];
      if(numCont == maxtouch)
        grow_arrays_maxtouch(atom->nmax);

      partner[iP][numCont] = idTri;
      delflag[iP][numCont] = false;

      if(dnum > 0)
        history = contacthistory[iP][numCont];
      for(int i = 0; i < dnum; i++)
        history[i] = 0.;

      npartner[iP]++;

  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistory::n_contacts()
  {
    int ncontacts = 0, nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
           ncontacts += npartner[i];
    return ncontacts;
  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistory::n_contacts(int contact_groupbit)
  {
    int ncontacts = 0, nlocal = atom->nlocal;
    int *mask = atom->mask;

    for(int i = 0; i < nlocal; i++)
        if(mask[i] & contact_groupbit)
           ncontacts += npartner[i];
    return ncontacts;
  }
#endif
