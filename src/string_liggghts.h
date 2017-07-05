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

#ifndef LMP_STRING_LIGGGHTS_H
#define LMP_STRING_LIGGGHTS_H

namespace LAMMPS_NS {

//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

/**
 * @brief str_ends_with
 * @param str
 * @param suffix
 * @return
 */
inline bool strEndWith(const char * str, const char * suffix) {

  if( str == NULL || suffix == NULL )
    return false;

  const size_t str_len = strlen(str);
  const size_t suffix_len = strlen(suffix);

  if(suffix_len > str_len)
    return false;

  return strncmp( str + str_len - suffix_len, suffix, suffix_len ) == 0;
}

inline void strtrim(char* str) {
  int start = 0; // number of leading spaces
  char* buffer = str;
  while (*str && *str++ == ' ')
      ++start;
  while (*str++); // move to end of string
  int end = str - buffer - 1;
  while (end > 0 && buffer[end - 1] == ' ')
      --end; // backup over trailing spaces
  buffer[end] = 0; // remove trailing spaces
  if (end <= start || start == 0)
      return; // exit if no leading spaces or string is now empty
  str = buffer + start;
  while ((*buffer++ = *str++));  // remove leading spaces: K&R
}

} // LAMMPS_NS

#endif
