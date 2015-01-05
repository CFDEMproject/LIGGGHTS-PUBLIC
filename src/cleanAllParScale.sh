#!/bin/bash

#will clean
# - the poems library
# - the fedora_fpic verions of LIGGGHTS
# - uninstall the poems package (and take care of an edit of fix_poems.*)

currDir=$PWD
cp *pascal* PASCAL

cp Makefile.package.empty Makefile.package

make no-PASCAL
