# Make.sh = update Makefile.lib, Makefile.shlib, Makefile.list
#           or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.shlib
#         sh Make.sh Makefile.list

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

style () {
  # modified C.K. create version_liggghts.h
  builddate=`date +%Y-%m-%d-%H:%M:%S`
  wai=`whoami`
  vers=`cat version_liggghts.txt`
  bra=`cat version_liggghts_branch.txt`
  githash=`git log -1 --format="%H"`
  echo "#define LIGGGHTS_VERSION \"$bra $vers, compiled $builddate by $wai, git commit $githash\"" > version_liggghts.h

  list=`grep -sl $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    if (test ! -e style_$3.h) then
      touch style_$3.h
    elif (test "`cat style_$3.h`" != "") then
    rm -f style_$3.h
    touch style_$3.h
      rm -f Obj_*/$4.d
      if (test $5) then
        rm -f Obj_*/$5.d
      fi
      rm -f Obj_*/lammps.d
    fi
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files
# called by "make machine"
# col 1 = string to search for
# col 2 = search in *.h files starting with this name
# col 3 = prefix of style file
# col 4

if (test $1 = "style") then

  style ANGLE_CLASS     angle_      angle      force
  style ATOM_CLASS      atom_vec_   atom       atom      atom_vec_hybrid
  style BODY_CLASS      body_       body       atom_vec_body
  style BOND_CLASS      bond_       bond       force
  style COMMAND_CLASS   ""          command    input
  style COMPUTE_CLASS   compute_    compute    modify    modify_cuda
  style DIHEDRAL_CLASS  dihedral_   dihedral   force
  style DUMP_CLASS      dump_       dump       output
  style FIX_CLASS       fix_        fix        modify
  style IMPROPER_CLASS  improper_   improper   force
  style INTEGRATE_CLASS ""          integrate  update
  style KSPACE_CLASS    ""          kspace     force
  style MINIMIZE_CLASS  min_        minimize   update
  style PAIR_CLASS      pair_       pair       force
  style SURFACE_MODEL    surface_model_     surface_model     force
  style NORMAL_MODEL     normal_model_      normal_model      force
  style TANGENTIAL_MODEL tangential_model_  tangential_model  force
  style COHESION_MODEL   cohesion_model_    cohesion_model    force
  style ROLLING_MODEL    rolling_model_     rolling_model     force
  style READER_CLASS    reader_     reader     read_dump
  style REGION_CLASS    region_     region     domain
  style CFD_DATACOUPLING_CLASS      cfd_datacoupling_  cfd_datacoupling  fix_cfd_coupling
  style CFD_REGIONMODEL_CLASS       cfd_regionmodel_  cfd_regionmodel  fix_cfd_coupling
  style LB_CLASS        ""          lb
  style SPH_KERNEL_CLASS  sph_kernel_  sph_kernel  pair_sph-fix_sph
elif (test $1 = "models") then
  sed_ex="sed -E" # BSD sed
  sed --version 2>&1 | grep -i gnu &> /dev/null
  [ $? -eq 0 ] && sed_ex="sed -r" # GNU sed

  surface_models=`grep -s -E '^SURFACE_MODEL' surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  normal_models=`grep -s -E '^NORMAL_MODEL' normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  tangential_models=`grep -s -E '^TANGENTIAL_MODEL' tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  cohesion_models=`grep -s -E '^COHESION_MODEL' cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  rolling_models=`grep -s -E '^ROLLING_MODEL' rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`

  # check for duplicate constants
  sm_duplicates=`grep -s -E '^SURFACE_MODEL' surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  nm_duplicates=`grep -s -E '^NORMAL_MODEL' normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  tm_duplicates=`grep -s -E '^TANGENTIAL_MODEL' tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  cm_duplicates=`grep -s -E '^COHESION_MODEL' cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  rm_duplicates=`grep -s -E '^ROLLING_MODEL' rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`

  if [ -n "$sm_duplicates" ]; then echo "ERROR: duplicate surface model identifiers:"; echo $sm_duplicates; exit -1; fi
  if [ -n "$nm_duplicates" ]; then echo "ERROR: duplicate normal model identifiers:"; echo $nm_duplicates; exit -1; fi
  if [ -n "$tm_duplicates" ]; then echo "ERROR: duplicate tangential model identifiers:"; echo $tm_duplicates; exit -1; fi
  if [ -n "$cm_duplicates" ]; then echo "ERROR: duplicate cohesion model identifiers:"; echo $cm_duplicates; exit -1; fi
  if [ -n "$rm_duplicates" ]; then echo "ERROR: duplicate rolling model identifiers:"; echo $rm_duplicates; exit -1; fi

  stylefile=style_contact_model.h
  tmpfile=style_contact_model.tmp
  filteredfile=style_contact_model_filtered.tmp

  if (test -e $tmpfile) then
    rm -f $tmpfile
  fi

  if (test -e $filteredfile) then
    rm -f $filteredfile
  fi

  cohesion_models="COHESION_OFF $cohesion_models"
  rolling_models="ROLLING_OFF $rolling_models"

  # build all model combinations
  for surf in $surface_models; do
    for norm in $normal_models; do
      for tang in $tangential_models; do
        for coh in $cohesion_models; do
          for roll in $rolling_models; do
            echo "GRAN_MODEL($norm, $tang, $coh, $roll, $surf)" >> $tmpfile
          done
        done
      done
    done
  done

  if (test -e style_contact_model.blacklist) then
    grep -v -f style_contact_model.blacklist $tmpfile > $filteredfile
    rm $tmpfile
  else
    mv $tmpfile $filteredfile
  fi


  if (test ! -e $filteredfile) then
    rm -f $stylefile
    touch $stylefile
  elif (test ! -e $stylefile) then
    mv $filteredfile $stylefile
    rm -f Obj_*/force.d
    rm -f Obj_*/modify.d
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief $stylefile $filteredfile`" != "") then
    mv $filteredfile $stylefile
    rm -f Obj_*/force.d
    rm -f Obj_*/modify.d
    rm -f Obj_*/lammps.d
  else
    rm -f $filteredfile
  fi

# edit Makefile.lib, for creating non-shared lib
# called by "make makelib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.lib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.lib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.lib
# edit Makefile.shlib, for creating shared lib
# called by "make makeshlib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.shlib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.shlib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.shlib

# edit Makefile.list
# called by "make makelist"
# use current list of *.cpp and *.h files in src dir

elif (test $1 = "Makefile.list") then

  list=`ls -1 *.cpp | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.list
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.list

fi
