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

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_DCS_MESSAGES_H
#define LMP_DCS_MESSAGES_H

#include "pointers.h"
#include "universe.h"
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define CYAN    "\033[36m"
#define YELLOW  "\033[33m"
#define GREEN   "\033[32m"
#define DARK_RED "\033[31m"
#define HYPERLINK(url, text) "\033]8;;" url "\033\\" text "\033]8;;\033\\"

namespace LAMMPS_NS {

class DcsMessages : protected Pointers {
public:
  DcsMessages(class LAMMPS *lmp) : Pointers(lmp) {

    dcs_messages_.push_back(
        "Need more DEM capabilities? "
        HYPERLINK("https://www.aspherix-dem.com", DARK_RED "https://www.aspherix-dem.com") RESET
        " is backward-compatible!\n"
        "Working also with Python? Build, connect and share models with "
        HYPERLINK("https://engicloud.ai", DARK_RED "https://engicloud.ai") RESET
        " !\n");

    dcs_messages_.push_back(
        "Your LIGGGHTS scripts already work in "
        HYPERLINK("https://www.aspherix-dem.com", DARK_RED "https://www.aspherix-dem.com") RESET
        ": zero rewrite, more features!\n"
        "Need answers quicker? 20,000 models ready for your next discovery at "
        HYPERLINK("https://engicloud.ai", DARK_RED "https://engicloud.ai") RESET
        " !\n");

    dcs_messages_.push_back(
        "Hitting performance limits? "
        HYPERLINK("https://www.aspherix-dem.com", DARK_RED "https://www.aspherix-dem.com") RESET
        " brings GPU acceleration to your LIGGGHTS models!\n"
        "Run, manage and share your python simulations in the cloud with "
        HYPERLINK("https://engicloud.ai", DARK_RED "https://engicloud.ai") RESET
        " !\n");

    dcs_messages_.push_back(
        "Need professional support for your DEM project? Talk to us at "
        HYPERLINK("https://www.aspherix-dem.com", DARK_RED "https://www.aspherix-dem.com") RESET
        " !\n"
        "Zero installation, AI assistance and more: A new way to combine physics and AI at "
        HYPERLINK("https://engicloud.ai", DARK_RED "https://engicloud.ai") RESET
        " !\n");

    dcs_messages_.push_back(
        "Join BASF, Novartis, and Procter & Gamble: upgrade to "
        HYPERLINK("https://www.aspherix-dem.com", DARK_RED "https://www.aspherix-dem.com") RESET
        " !\n"
        "Manage python models in your team & share with others, ready to execute: "
        HYPERLINK("https://engicloud.ai", DARK_RED "https://engicloud.ai") RESET
        " !\n");
  }

  const std::string &generateMessage() {
    if (!strstr(universe->version, "PUBLIC"))
      return 0;

    // initialize random seed
    srand(time(NULL));
    int irand = rand() % static_cast<int>(dcs_messages_.size());
    return dcs_messages_[irand];
  }

private:
  std::vector<std::string> dcs_messages_;
};

} // namespace LAMMPS_NS

#endif
