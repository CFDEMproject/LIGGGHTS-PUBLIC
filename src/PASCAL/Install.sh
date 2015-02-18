# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
  action $file
done

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then
  echo "...updating Makefile.package"
  if (test -e ../Makefile.package) then
    # remove old paths
    sed -i -e 's/[^ \t]*PASCAL[^ \t]* //g' ../Makefile.package # global for all occurences of PASCAL
    sed -i -e 's/[^ \t]*pasc[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*sundial[^ \t]* //g' ../Makefile.package # global
    sed -i -e 's/[^ \t]*chemkinreader[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*Qt5[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*hdf5[^ \t]* //g' ../Makefile.package # global
    sed -i -e 's/[^ \t]*boost[^ \t]* //' ../Makefile.package
    # add new paths
    sed -i -e 's|^PKG_INC =[ \t]*|&-I$(PASCAL_SRC_DIR) |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L$(PASCAL_SRC_DIR) -L$(PASCAL_INST_DIR)/lib64 -L$(PASCAL_INST_DIR)/lib -L$(PASCAL_QT5_DIR)/lib -L$(PASCAL_HDF5_DIR)/lib -L$(PASCAL_THIRDPARTY_DIR)/chemkinReader/src |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lpasc_fedora_fpic -lsundials_cvode -lsundials_nvecserial -lQt5Core -lhdf5 -lhdf5_hl -lhdf5_cpp -lchemkinreader -lboost_regex |' ../Makefile.package
  fi

#  if (test -e ../Makefile.package.settings) then
#    sed -i -e '/^include.*pascal.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
#    sed -i -e '4 i \include $(PASCAL_SRC_DIR) ' ../Makefile.package.settings
#  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*PASCAL[^ \t]* //g' ../Makefile.package # global for all occurences of PASCAL
    sed -i -e 's/[^ \t]*pasc[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*sundial[^ \t]* //g' ../Makefile.package # global
    sed -i -e 's/[^ \t]*chemkinreader[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*Qt5[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*hdf5[^ \t]* //g' ../Makefile.package # global
    sed -i -e 's/[^ \t]*boost[^ \t]* //' ../Makefile.package
  fi

#  if (test -e ../Makefile.package.settings) then
#    sed -i -e '/^include.*pascal.*$/d' ../Makefile.package.settings
#  fi

fi
