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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2017-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "unistd.h"
#include <stdlib.h>
#include <string>

namespace LAMMPS_NS
{

namespace SignalHandler
{

bool request_write_restart(false);
bool request_quit(false);
bool enable_restart_writing(true);

void int_handler(int signum)
{
    if (request_quit)
    {
        if (request_write_restart)
        {
            std::string warning("\nRestart file was not written yet.");
            const int len = warning.length();
#ifndef  _WINDOWS
            write(STDOUT_FILENO, warning.c_str(), len);
#endif
        }

        std::string error("\nSecond SIGINT/SIGTERM caught - Force quitting.\n");
        const int len = error.length();
#ifndef  _WINDOWS
        write(STDERR_FILENO, error.c_str(), len);

#endif
        _Exit(EXIT_FAILURE);
    }
    else
    {
        std::string warning;
        if (enable_restart_writing)
        {
            warning = "\nSIGINT/SIGTERM caught - Writing restart on next occasion and quitting after that.\n";
            request_write_restart = true;
        }
        else
            warning = "\nSIGINT/SIGTERM caught - quitting on first occasion.\n";
        const int len = warning.length();
#ifndef  _WINDOWS
        write(STDOUT_FILENO, warning.c_str(), len);

#endif
        request_quit = true;
    }
}

void usr1_handler(int signum)
{
    std::string warning;
    if (enable_restart_writing)
    {
        warning = "\nSIGUSR1 caught - Writting restart file on next occasion.\n";
        request_write_restart = true;
    }
    else
        warning = "\nSIGUSR1 caught - No action performed.\n";
    const int len = warning.length();
#ifndef  _WINDOWS
    write(STDOUT_FILENO, warning.c_str(), len);
#endif
}

};
};
