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

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(file,CfdDatacouplingFile)

#else

#ifndef LMP_CFD_DATACOUPLING_FILE_H
#define LMP_CFD_DATACOUPLING_FILE_H

#include "cfd_datacoupling.h"

namespace LAMMPS_NS {

class CfdDatacouplingFile : public CfdDatacoupling {
 public:
  CfdDatacouplingFile(class LAMMPS *, int, int, char **,class FixCfdCoupling* fc);
  ~CfdDatacouplingFile();
  friend class FixTempFromFile;

  void pull(const char *, const char *, void *&, const char *);
  void push(const char *, const char *, void *&, const char *);
  virtual void post_create();

  void exchange();

  private:
   char* filepath;
   
   int append;

   bool firstexec;
   int t0;

   char * getFilePath(const char *name, bool flag);
   void op_complete(const char *name);
   void writeVectorData(const char *name,  double ** field);
   void writeScalarData(const char *name,  double * field);
   void writeGlobalVectorData(const char *name,  double * field,int);
   void writeGlobalArrayData(const char *name,  double ** field,int,int);

   void readVectorData(const char *name,  double ** field);
   void readScalarData(const char *name,  double * field);
   void readGlobalVectorData(const char* name, double *field, int &len);
   void readGlobalArrayData(const char *name, double ** field, int &len1, int &len2);

};

}

#endif
#endif
