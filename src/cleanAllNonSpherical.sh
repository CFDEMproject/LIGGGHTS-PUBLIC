#!/bin/bash

#will clean
# - the poems library
# - the fedora_fpic verions of LIGGGHTS
# - uninstall the poems package (and take care of an edit of fix_poems.*)

currDir=$PWD

cd ../lib/poems
make -f Makefile.g++ clean

cd $currDir
#make clean-all
cp fix_poems.* POEMS
cp fix_*asphere.* ASPHERE
cp fix_*Asphere.* ASPHERE
cp pair_gayberne* ASPHERE
cp *pascal* PASCAL
cp *_smd_* smd_* USER-SMD

make no-POEMS
make no-DIPOLE
make no-ASPHERE
make no-PASCAL
make no-USER-SMD
rm make.log

cp Makefile.package.empty Makefile.package
rm lmp_*
rm *.so *.a
