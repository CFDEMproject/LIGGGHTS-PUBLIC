#!/bin/bash

#will compile
# - the poems library
# - the fedora_fpic verions of LIGGGHTS including the fix_poems

#WARNING: will overwrite fix_poems.* in the src dir!!!

nProc=3
currDir=$PWD

#check if PaSCal dir exists
if [[ -n "$PASCAL_SRC_DIR" ]] && [[ "$1" == "yes-PASCAL" ]] ; then
  if [ -d "$PASCAL_SRC_DIR" ]; then
    echo "*******************"
    echo "PaScal src directory is set and exists. Be sure that PaScal is ALREADY compiled as library! Continuing compilation..."
    echo "Will continue compilation now with PaScal as libary..."
    make yes-PASCAL
    echo "*******************"
    sleep 2
  else
    echo "*******************"
    echo "PaScal src directory is set but DOES NOT exist exists."
    echo "WARNING: Will NOT compile with PaScal!"
    echo "*******************"
    sleep 2
  fi
else
    echo "*******************"
    echo "WARNING: PaScal src directory IS NOT SET or you have NOT Set 'yes-PASCAL'. Skipping compilation of PaScal..."
    echo "WARNING: Will NOT compile with PaScal, i.e., will uninstall PASCAL package!"
    echo "*******************"
    make no-PASCAL
    sleep 2
fi


make -j $nProc fedora_fpic
make makelib
make -f Makefile.lib fedora_fpic
