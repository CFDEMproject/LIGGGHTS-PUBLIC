/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "sys/stat.h"
#include "input.h"
#include "style_command.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "force.h"
#include "pair.h"
#include "min.h"
#include "modify.h"
#include "compute.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "special.h"
#include "variable.h"
#include "accelerator_cuda.h"
#include "error.h"
#include "memory.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace LAMMPS_NS;

//#define MAXLINE 2048
#define MAXLINE 8192 
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(LAMMPS *lmp, int argc, char **argv) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  work = new char[MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

  variable = new Variable(lmp);

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 0;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-var") == 0 || strcmp(argv[iarg],"-v") == 0) {
      int jarg = iarg+3;
      while (jarg < argc && argv[jarg][0] != '-') jarg++;
      variable->set(argv[iarg+1],jarg-iarg-2,&argv[iarg+2]);
      iarg = jarg;
    } else if (strcmp(argv[iarg],"-echo") == 0 ||
               strcmp(argv[iarg],"-e") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
    } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  delete variable;
  delete [] line;
  delete [] copy;
  delete [] work;
  if (labelstr) delete [] labelstr;
  memory->sfree(arg);
  memory->sfree(infiles);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int m,n;

  while (1) {

    // read a line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = length of line including str terminator, 0 if end of file
    // m = position of last printable char in line or -1 if blank line

    if (me == 0) {
      m = 0;
      while (1) {
        if (fgets(&line[m],MAXLINE-m,infile) == NULL) n = 0;
        else n = strlen(line) + 1;
        if (n == 0) break;
        m = n-2;
        while (m >= 0 && isspace(line[m])) m--;
        if (m < 0 || line[m] != '&') break;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all(FLERR,"Label wasn't found in input script");
      if (me == 0) {
        if (infile != stdin) fclose(infile);
        nfile--;
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long

    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(FLERR,str);
    }

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s",line);
      if (echo_log && logfile) fprintf(logfile,"%s",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == NULL) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char str[MAXLINE];
      sprintf(str,"Unknown command: %s",line);
      error->all(FLERR,str);
    }
  }
}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void Input::file(const char *filename)
{
  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set infile, infiles[0]

  if (me == 0) {
    if (nfile > 1)
      error->one(FLERR,"Another input script is already being processed");
    if (infile != stdin) fclose(infile);
    infile = fopen(filename,"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",filename);
      error->one(FLERR,str);
    }
    infiles[0] = infile;
  } else infile = NULL;

  file();
}

/* ----------------------------------------------------------------------
   parse the command in single and execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(const char *single)
{
  strcpy(line,single);

  // echo the command unless scanning for label

  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s\n",line);
    if (echo_log && logfile) fprintf(logfile,"%s\n",line);
  }

  // parse the line
  // if no command, just return NULL

  parse();
  if (command == NULL) return NULL;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return NULL;

  // execute the command and return its name

  if (execute_command()) {
    char str[MAXLINE];
    sprintf(str,"Unknown command: %s",line);
    error->all(FLERR,str);
  }

  return command;
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by replacing it with 0
  // do not strip # inside single/double quotes

  char quote = '\0';
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && !quote) {
      *ptr = '\0';
      break;
    }
    if (*ptr == quote) quote = '\0';
    else if (*ptr == '"' || *ptr == '\'') quote = *ptr;
    ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy,1);

  // command = 1st arg

  char *next;
  command = nextword(copy,&next);
  if (command == NULL) return;

  // point arg[] at each subsequent arg
  // nextword() inserts zeroes in copy to delimit args
  // nextword() treats text between single/double quotes as one arg

  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double quotes as one arg
     matching quote must be followed by whitespace char if not end of string
     strip quotes from returned word
   return ptr to start of word
   return next = ptr after word or NULL if word ended with 0
   return NULL if no word in string
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return NULL;

  if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    if (stop[1] && !isspace(stop[1]))
      error->all(FLERR,"Input line quote not followed by whitespace");
    start++;
  } else stop = &start[strcspn(start," \t\n\v\f\r")];

  if (*stop == '\0') *next = NULL;
  else *next = stop+1;
  *stop = '\0';
  return start;
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str and return it
   str assumed to be long enough to hold expanded version
   print updated string if flag is set and not searching for label
------------------------------------------------------------------------- */

void Input::substitute(char *str, int flag)
{
  // use work[] as scratch space to expand str, then copy back to str
  // do not replace $ inside single/double quotes
  // var = pts at variable name, ended by NULL
  //   if $ is followed by '{', trailing '}' becomes NULL
  //   else $x becomes x followed by NULL
  // beyond = pts at text following variable

  char *var,*value,*beyond;
  char quote = '\0';
  char *ptr = str;

  while (*ptr) {
    if (*ptr == '$' && !quote) {
      if (*(ptr+1) == '{') {
        var = ptr+2;
        int i = 0;
        while (var[i] != '\0' && var[i] != '}') i++;
        if (var[i] == '\0') error->one(FLERR,"Invalid variable name");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
      } else {
        var = ptr;
        var[0] = var[1];
        var[1] = '\0';
        beyond = ptr + strlen(var) + 1;
      }
      value = variable->retrieve(var);
      if (value == NULL) error->one(FLERR,"Substitution for illegal variable");

      *ptr = '\0';
      strcpy(work,str);
      if (strlen(work)+strlen(value) >= MAXLINE)
        error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,value);
      if (strlen(work)+strlen(beyond) >= MAXLINE)
        error->one(FLERR,"Input line too long after variable substitution");
      strcat(work,beyond);
      strcpy(str,work);
      ptr += strlen(value);
      if (flag && me == 0 && label_active == 0) {
        if (echo_screen && screen) fprintf(screen,"%s",str);
        if (echo_log && logfile) fprintf(logfile,"%s",str);
      }
      continue;
    }
    if (*ptr == quote) quote = '\0';
    else if (*ptr == '"' || *ptr == '\'') quote = *ptr;
    ptr++;
  }
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"echo")) echo();
  else if (!strcmp(command,"if")) ifthenelse();
  else if (!strcmp(command,"include")) include();
  else if (!strcmp(command,"jump")) jump();
  else if (!strcmp(command,"label")) label();
  else if (!strcmp(command,"log")) log();
  else if (!strcmp(command,"next")) next_command();
  else if (!strcmp(command,"partition")) partition();
  else if (!strcmp(command,"print")) print();
  else if (!strcmp(command,"quit")) quit();
  else if (!strcmp(command,"shell")) shell();
  else if (!strcmp(command,"variable")) variable_command();

  else if (!strcmp(command,"angle_coeff")) angle_coeff();
  else if (!strcmp(command,"angle_style")) angle_style();
  else if (!strcmp(command,"atom_modify")) atom_modify();
  else if (!strcmp(command,"atom_style")) atom_style();
  else if (!strcmp(command,"bond_coeff")) bond_coeff();
  else if (!strcmp(command,"bond_style")) bond_style();
  else if (!strcmp(command,"boundary")) boundary();
  else if (!strcmp(command,"communicate")) communicate();
  else if (!strcmp(command,"compute")) compute();
  else if (!strcmp(command,"compute_modify")) compute_modify();
  else if (!strcmp(command,"dielectric")) dielectric();
  else if (!strcmp(command,"dihedral_coeff")) dihedral_coeff();
  else if (!strcmp(command,"dihedral_style")) dihedral_style();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"dump")) dump();
  else if (!strcmp(command,"dump_modify")) dump_modify();
  else if (!strcmp(command,"fix")) fix();
  else if (!strcmp(command,"fix_modify")) fix_modify();
  else if (!strcmp(command,"group")) group_command();
  else if (!strcmp(command,"improper_coeff")) improper_coeff();
  else if (!strcmp(command,"improper_style")) improper_style();
  else if (!strcmp(command,"kspace_modify")) kspace_modify();
  else if (!strcmp(command,"kspace_style")) kspace_style();
  else if (!strcmp(command,"lattice")) lattice();
  else if (!strcmp(command,"mass")) mass();
  else if (!strcmp(command,"min_modify")) min_modify();
  else if (!strcmp(command,"min_style")) min_style();
  else if (!strcmp(command,"neigh_modify")) neigh_modify();
  else if (!strcmp(command,"neighbor")) neighbor_command();
  else if (!strcmp(command,"newton")) newton();
  else if (!strcmp(command,"package")) package();
  else if (!strcmp(command,"pair_coeff")) pair_coeff();
  else if (!strcmp(command,"pair_modify")) pair_modify();
  else if (!strcmp(command,"pair_style")) pair_style();
  else if (!strcmp(command,"pair_write")) pair_write();
  else if (!strcmp(command,"processors")) processors();
  else if (!strcmp(command,"region")) region();
  else if (!strcmp(command,"reset_timestep")) reset_timestep();
  else if (!strcmp(command,"restart")) restart();
  else if (!strcmp(command,"run_style")) run_style();
  else if (!strcmp(command,"special_bonds")) special_bonds();
  else if (!strcmp(command,"suffix")) suffix();
  else if (!strcmp(command,"thermo")) thermo();
  else if (!strcmp(command,"thermo_modify")) thermo_modify();
  else if (!strcmp(command,"thermo_style")) thermo_style();
  else if (!strcmp(command,"timestep")) timestep();
  else if (!strcmp(command,"uncompute")) uncompute();
  else if (!strcmp(command,"undump")) undump();
  else if (!strcmp(command,"unfix")) unfix();
  else if (!strcmp(command,"units")) units();

  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // check if command is added via style.h

  if (0) return 0;      // dummy line to enable else-if macro expansion

#define COMMAND_CLASS
#define CommandStyle(key,Class)         \
  else if (strcmp(command,#key) == 0) { \
    Class key(lmp);                        \
    key.command(narg,arg);              \
    return 0;                           \
  }
#include "style_command.h"
#undef COMMAND_CLASS

  // unrecognized command

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR,"Illegal clear command");
  lmp->destroy();
  lmp->create();
  lmp->post_create();
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all(FLERR,"Illegal echo command");

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all(FLERR,"Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg < 3) error->all(FLERR,"Illegal if command");

  // substitute for variables in Boolean expression for "if"
  // in case expression was enclosed in quotes
  // must substitute on copy of arg else will step on subsequent args

  char *scopy = new char[MAXLINE];
  strcpy(scopy,arg[0]);
  substitute(scopy,0);

  // evaluate Boolean expression for "if"
  
  double btest;

  if(strncmp(scopy,"property_",9) == 0)
  {
      if(strlen(scopy) < 11)
        error->all(FLERR,"Illegal if command, length of argument too short");
      if(modify->find_fix_property(&(scopy[9]),"property/global","any",0,0,"if",false))
        btest = 1.0;
      else
        btest = 0.0;
  }
  
  else
    btest = variable->evaluate_boolean(scopy);

  // bound "then" commands

  if (strcmp(arg[1],"then") != 0) error->all(FLERR,"Illegal if command");

  int first = 2;
  int iarg = first;
  while (iarg < narg &&
         (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
    iarg++;
  int last = iarg-1;

  // execute "then" commands
  // make copies of all arg string commands
  // required because re-parsing a command via one() will wipe out args

  if (btest != 0.0) {
    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR,"Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }

    for (int i = 0; i < ncommands; i++) one(commands[i]);

    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
    delete [] scopy;

    return;
  }

  // done if no "elif" or "else"

  if (iarg == narg) {
    delete [] scopy;
    return;
  }

  // check "elif" or "else" until find commands to execute
  // substitute for variables and evaluate Boolean expression for "elif"
  // must substitute on copy of arg else will step on subsequent args
  // bound and execute "elif" or "else" commands

  while (1) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal if command");
    if (strcmp(arg[iarg],"elif") == 0) {
      strcpy(scopy,arg[iarg+1]);
      substitute(scopy,0);
      btest = variable->evaluate_boolean(scopy);
      first = iarg+2;
    } else {
      btest = 1.0;
      first = iarg+1;
    }

    iarg = first;
    while (iarg < narg &&
           (strcmp(arg[iarg],"elif") != 0 && strcmp(arg[iarg],"else") != 0))
      iarg++;
    last = iarg-1;

    if (btest == 0.0) continue;

    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR,"Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR,"Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }

    // execute the list of commands

    for (int i = 0; i < ncommands; i++) one(commands[i]);

    // clean up

    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
    delete [] scopy;

    return;
  }
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all(FLERR,"Illegal include command");

  if (me == 0) {
    if (nfile == maxfile) {
      maxfile++;
      infiles = (FILE **)
        memory->srealloc(infiles,maxfile*sizeof(FILE *),"input:infiles");
    }
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(FLERR,str);
    }
    infiles[nfile++] = infile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (strcmp(arg[0],"SELF") == 0) rewind(infile);
    else {
      if (infile != stdin) fclose(infile);
      infile = fopen(arg[0],"r");
      if (infile == NULL) {
        char str[128];
        sprintf(str,"Cannot open input script %s",arg[0]);
        error->one(FLERR,str);
      }
      infiles[nfile-1] = infile;
    }
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr,arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all(FLERR,"Illegal label command");
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg != 1) error->all(FLERR,"Illegal log command");

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = NULL;
    else {
      logfile = fopen(arg[0],"w");
      if (logfile == NULL) {
        char str[128];
        sprintf(str,"Cannot open logfile %s",arg[0]);
        error->one(FLERR,str);
      }
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::partition()
{
  if (narg < 3) error->all(FLERR,"Illegal partition command");

  int yesflag;
  if (strcmp(arg[0],"yes") == 0) yesflag = 1;
  else if (strcmp(arg[0],"no") == 0) yesflag = 0;
  else error->all(FLERR,"Illegal partition command");

  int ilo,ihi;
  force->bounds(arg[1],universe->nworlds,ilo,ihi);

  // copy original line to copy, since will use strtok() on it
  // ptr = start of 4th word

  strcpy(copy,line);
  copy[strlen(copy)-1] = '\0';
  char *ptr = strtok(copy," \t\n\r\f");
  ptr = strtok(NULL," \t\n\r\f");
  ptr = strtok(NULL," \t\n\r\f");
  ptr += strlen(ptr) + 1;
  ptr += strspn(ptr," \t\n\r\f");

  // execute the remaining command line on requested partitions

  if (yesflag) {
    if (universe->iworld+1 >= ilo && universe->iworld+1 <= ihi) one(ptr);
  } else {
    if (universe->iworld+1 < ilo || universe->iworld+1 > ihi) one(ptr);
  }
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg != 1) error->all(FLERR,"Illegal print command");

  // substitute for $ variables (no printing) and print arg

  substitute(arg[0],0);
  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",arg[0]);
    if (logfile) fprintf(logfile,"%s\n",arg[0]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::quit()
{
  if (narg) error->all(FLERR,"Illegal quit command");
  error->done();
}

/* ---------------------------------------------------------------------- */

void Input::shell()
{
  if (narg < 1) error->all(FLERR,"Illegal shell command");

  if (strcmp(arg[0],"cd") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal shell command");
    chdir(arg[1]);

  } else if (strcmp(arg[0],"mkdir") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
#if !defined(WINDOWS) && !defined(__MINGW32_VERSION)
    if (me == 0)
      for (int i = 1; i < narg; i++)
        mkdir(arg[i], S_IRWXU | S_IRGRP | S_IXGRP);
#endif

  } else if (strcmp(arg[0],"mv") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal shell command");
    if (me == 0) rename(arg[1],arg[2]);

  } else if (strcmp(arg[0],"rm") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
    if (me == 0)
      for (int i = 1; i < narg; i++) unlink(arg[i]);

  } else if (strcmp(arg[0],"rmdir") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal shell command");
    if (me == 0)
      for (int i = 1; i < narg; i++) rmdir(arg[i]);

  // use work to concat args back into one string separated by spaces
  // invoke string in shell via system()

  } else {
    strcpy(work,arg[0]);
    for (int i = 1; i < narg; i++) {
      strcat(work," ");
      strcat(work,arg[i]);
    }
    if (me == 0) system(work);
  }
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg,arg);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one function for each LAMMPS-specific input script command
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::angle_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Angle_coeff command before simulation box is defined");
  if (force->angle == NULL)
    error->all(FLERR,"Angle_coeff command before angle_style is defined");
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Angle_coeff command when no angles allowed");
  force->angle->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::angle_style()
{
  if (narg < 1) error->all(FLERR,"Illegal angle_style command");
  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Angle_style command when no angles allowed");
  force->create_angle(arg[0],lmp->suffix);
  if (force->angle) force->angle->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::atom_modify()
{
  atom->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::atom_style()
{
  if (narg < 1) error->all(FLERR,"Illegal atom_style command");
  if (domain->box_exist)
    error->all(FLERR,"Atom_style command after simulation box is defined");
  atom->create_avec(arg[0],narg-1,&arg[1],lmp->suffix);
}

/* ---------------------------------------------------------------------- */

void Input::bond_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Bond_coeff command before simulation box is defined");
  if (force->bond == NULL)
    error->all(FLERR,"Bond_coeff command before bond_style is defined");
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Bond_coeff command when no bonds allowed");
  force->bond->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::bond_style()
{
  if (narg < 1) error->all(FLERR,"Illegal bond_style command");
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Bond_style command when no bonds allowed");
  force->create_bond(arg[0],lmp->suffix);
  if (force->bond) force->bond->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  if (domain->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");
  domain->set_boundary(narg,arg,0);
}

/* ---------------------------------------------------------------------- */

void Input::communicate()
{
  comm->set(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::compute()
{
  modify->add_compute(narg,arg,lmp->suffix);
}

/* ---------------------------------------------------------------------- */

void Input::compute_modify()
{
  modify->modify_compute(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dielectric()
{
  if (narg != 1) error->all(FLERR,"Illegal dielectric command");
  force->dielectric = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Dihedral_coeff command before simulation box is defined");
  if (force->dihedral == NULL)
    error->all(FLERR,"Dihedral_coeff command before dihedral_style is defined");
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Dihedral_coeff command when no dihedrals allowed");
  force->dihedral->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_style()
{
  if (narg < 1) error->all(FLERR,"Illegal dihedral_style command");
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,"Dihedral_style command when no dihedrals allowed");
  force->create_dihedral(arg[0],lmp->suffix);
  if (force->dihedral) force->dihedral->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (narg != 1) error->all(FLERR,"Illegal dimension command");
  if (domain->box_exist)
    error->all(FLERR,"Dimension command after simulation box is defined");
  domain->dimension = atoi(arg[0]);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR,"Illegal dimension command");

  // must reset default extra_dof of all computes
  // since some were created before dimension command is encountered

  for (int i = 0; i < modify->ncompute; i++)
    modify->compute[i]->reset_extra_dof();
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  output->modify_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
  modify->add_fix(narg,arg,lmp->suffix);
}

/* ---------------------------------------------------------------------- */

void Input::fix_modify()
{
  modify->modify_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::group_command()
{
  group->assign(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::improper_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Improper_coeff command before simulation box is defined");
  if (force->improper == NULL)
    error->all(FLERR,"Improper_coeff command before improper_style is defined");
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Improper_coeff command when no impropers allowed");
  force->improper->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::improper_style()
{
  if (narg < 1) error->all(FLERR,"Illegal improper_style command");
  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Improper_style command when no impropers allowed");
  force->create_improper(arg[0],lmp->suffix);
  if (force->improper) force->improper->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::kspace_modify()
{
  if (force->kspace == NULL)
    error->all(FLERR,"KSpace style has not yet been set");
  force->kspace->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::kspace_style()
{
  force->create_kspace(narg,arg,lmp->suffix);
}

/* ---------------------------------------------------------------------- */

void Input::lattice()
{
  domain->set_lattice(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::mass()
{
  if (narg != 2) error->all(FLERR,"Illegal mass command");
  if (domain->box_exist == 0)
    error->all(FLERR,"Mass command before simulation box is defined");
  atom->set_mass(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_modify()
{
  update->minimize->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Min_style command before simulation box is defined");
  update->create_minimize(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
  neighbor->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
  neighbor->set(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::newton()
{
  int newton_pair,newton_bond;

  if (narg == 1) {
    if (strcmp(arg[0],"off") == 0) newton_pair = newton_bond = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair = newton_bond = 1;
    else error->all(FLERR,"Illegal newton command");
  } else if (narg == 2) {
    if (strcmp(arg[0],"off") == 0) newton_pair = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair= 1;
    else error->all(FLERR,"Illegal newton command");
    if (strcmp(arg[1],"off") == 0) newton_bond = 0;
    else if (strcmp(arg[1],"on") == 0) newton_bond = 1;
    else error->all(FLERR,"Illegal newton command");
  } else error->all(FLERR,"Illegal newton command");

  force->newton_pair = newton_pair;

  if (newton_bond == 0) {
    if (domain->box_exist && force->newton_bond == 1)
      error->all(FLERR,"Newton bond change after simulation box is defined");
    force->newton_bond = 0;
  } else {
    if (domain->box_exist && force->newton_bond == 0)
      error->all(FLERR,"Newton bond change after simulation box is defined");
    force->newton_bond = 1;
  }

  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;
}

/* ---------------------------------------------------------------------- */

void Input::package()
{
  if (domain->box_exist)
    error->all(FLERR,"Package command after simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal package command");

  if (strcmp(arg[0],"cuda") == 0) {
    if (!lmp->cuda)
      error->all(FLERR,"Package cuda command without USER-CUDA installed");
    lmp->cuda->accelerator(narg-1,&arg[1]);

  } else if (strcmp(arg[0],"gpu") == 0) {
    char **fixarg = new char*[2+narg];
    fixarg[0] = (char *) "package_gpu";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "GPU";
    for (int i = 1; i < narg; i++) fixarg[i+2] = arg[i];
    modify->allow_early_fix = 1;
    modify->add_fix(2+narg,fixarg,NULL);
    modify->allow_early_fix = 0;
    delete [] fixarg;

  } else if (strcmp(arg[0],"omp") == 0) {
    char **fixarg = new char*[2+narg];
    fixarg[0] = (char *) "package_omp";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "OMP";
    for (int i = 1; i < narg; i++) fixarg[i+2] = arg[i];
    modify->allow_early_fix = 1;
    modify->add_fix(2+narg,fixarg,NULL);
    modify->allow_early_fix = 0;
    delete [] fixarg;

  } else error->all(FLERR,"Illegal package command");
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Pair_coeff command before simulation box is defined");
  if (force->pair == NULL)
    error->all(FLERR,"Pair_coeff command before pair_style is defined");
  force->pair->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::pair_modify()
{
  if (force->pair == NULL)
    error->all(FLERR,"Pair_modify command before pair_style is defined");
  force->pair->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   if old pair style exists and new style is same, just change settings
   else create new pair class
------------------------------------------------------------------------- */

void Input::pair_style()
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");
  if (force->pair && strcmp(arg[0],force->pair_style) == 0) {
    force->pair->settings(narg-1,&arg[1]);
    return;
  }
  force->create_pair(arg[0],lmp->suffix);
  if (force->pair) force->pair->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::pair_write()
{
  if (force->pair == NULL)
    error->all(FLERR,"Pair_write command before pair_style is defined");
  force->pair->write_file(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::processors()
{
  if (domain->box_exist)
    error->all(FLERR,"Processors command after simulation box is defined");
  comm->set_processors(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::region()
{
  domain->add_region(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::reset_timestep()
{
  update->reset_timestep(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::restart()
{
  output->create_restart(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::run_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Run_style command before simulation box is defined");
  update->create_integrate(narg,arg,lmp->suffix);
}

/* ---------------------------------------------------------------------- */

void Input::special_bonds()
{
  // store 1-3,1-4 and dihedral/extra flag values before change
  // change in 1-2 coeffs will not change the special list

  double lj2 = force->special_lj[2];
  double lj3 = force->special_lj[3];
  double coul2 = force->special_coul[2];
  double coul3 = force->special_coul[3];
  int angle = force->special_angle;
  int dihedral = force->special_dihedral;
  int extra = force->special_extra;

  force->set_special(narg,arg);

  // if simulation box defined and saved values changed, redo special list

  if (domain->box_exist && atom->molecular) {
    if (lj2 != force->special_lj[2] || lj3 != force->special_lj[3] ||
        coul2 != force->special_coul[2] || coul3 != force->special_coul[3] ||
        angle != force->special_angle ||
        dihedral != force->special_dihedral ||
        extra != force->special_extra) {
      Special special(lmp);
      special.build();
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::suffix()
{
  if (narg != 1) error->all(FLERR,"Illegal suffix command");

  if (strcmp(arg[0],"off") == 0) lmp->suffix_enable = 0;
  else if (strcmp(arg[0],"on") == 0) lmp->suffix_enable = 1;
  else {
    delete [] lmp->suffix;
    int n = strlen(arg[0]) + 1;
    lmp->suffix = new char[n];
    strcpy(lmp->suffix,arg[0]);
    lmp->suffix_enable = 1;
  }
}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
  if (narg != 1) error->all(FLERR,"Illegal thermo command");
  int n = atoi(arg[0]);
  if (n < 0) error->all(FLERR,"Illegal thermo command");
  output->thermo_every = n;
}

/* ---------------------------------------------------------------------- */

void Input::thermo_modify()
{
  output->thermo->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
  output->create_thermo(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all(FLERR,"Illegal timestep command");
  update->dt = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::uncompute()
{
  if (narg != 1) error->all(FLERR,"Illegal uncompute command");
  modify->delete_compute(arg[0],true); 
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
  if (narg != 1) error->all(FLERR,"Illegal undump command");
  output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
  if (narg != 1) error->all(FLERR,"Illegal unfix command");
  modify->delete_fix(arg[0],true); 
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
  if (narg != 1) error->all(FLERR,"Illegal units command");
  if (domain->box_exist)
    error->all(FLERR,"Units command after simulation box is defined");
  update->set_units(arg[0]);
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse_nonlammps()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by replacing it with 0
  // do not strip # inside single/double quotes

  char quote = '\0';
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && !quote) {
      *ptr = '\0';
      break;
    }
    if (*ptr == quote) quote = '\0';
    else if (*ptr == '"' || *ptr == '\'') quote = *ptr;
    ptr++;
  }

  // point arg[] at each arg
  // nextword() inserts zeroes in copy to delimit args
  // nextword() treats text between single/double quotes as one arg

  narg = 0;
  ptr = copy;
  char *next;

  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}
