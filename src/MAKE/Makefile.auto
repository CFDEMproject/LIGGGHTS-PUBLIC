# auto = Recommended automatic makefile (settings in MAKE/Makefile.user or MAKE/Makefile.user_${postfix})
#
# Attempts to detect MPI, VTK, etc. automatically
#
# WARNING:
# This file should NOT be modified by the user. Instead see the Makefile.user file in the present
# directory of this file and set your variables there. Should this file not exist, simply run
# $ make auto
# once and it will be created for you including default values.
#
# Should you have any issues with this makefile please send us a bug report on www.cfdem.com so
# that we can improve it.

# Subtargets
EMPTY=
SUBTARGET=
ifneq (${postfix}, ${EMPTY})
    SUBTARGET = _${postfix}
endif

TMP := $(shell ls ../MAKE/Makefile.user${SUBTARGET} 2> /dev/null 1> /dev/null && echo 1 || echo 0)
ifeq ($(TMP), 0)
    # Use default Makefile.user
    ifeq (${postfix}, ${EMPTY})
        TMP := $(shell tail -n +11 ../MAKE/Makefile.user_default > ../MAKE/Makefile.user || echo -1)
    # Use Makefile.user_${postfix}
    else
        # Check if Makefile.user exists, if yes clone this one, otherwise Makefile.user_default
        TMP := $(shell ls ../MAKE/Makefile.user 2> /dev/null 1> /dev/null && echo 1 || echo 0)
        ifneq ($(TMP), 0)
            TMP := $(shell tail -n +11 ../MAKE/Makefile.user > ../MAKE/Makefile.user${SUBTARGET} || echo -1)
        else
            TMP := $(shell tail -n +11 ../MAKE/Makefile.user_default > ../MAKE/Makefile.user${SUBTARGET} || echo -1)
        endif
    endif
    ifeq ($(TMP), -1)
        $(error Could not create Makefile.user${SUBTARGET})
    else
        $(info )
        $(info >>> WARNING <<<)
        $(info Created MAKE/Makefile.user${SUBTARGET} for the first time. Please check that the options in that file are set correctly.)
        $(info )
    endif
endif

sinclude ../MAKE/Makefile.user${SUBTARGET}

# Version consistency check
# The following line needs to be identical in Makefile.user_default
MAKEFILE_USER_DEFAULT_VERSION = "1.1"
ifeq ($(origin MAKEFILE_USER_VERSION), undefined)
    $(info Could not find MAKEFILE_USER_VERSION in MAKE/Makefile.user${SUBTARGET}. No version consistency check will be performed and there might be errors. Please remove [backup] this file and run "make auto" again to generate a new one.)
else
    ifneq ($(MAKEFILE_USER_VERSION), $(MAKEFILE_USER_DEFAULT_VERSION))
        $(info $(MAKEFILE_USER_VERSION) $(MAKEFILE_USER_DEFAULT_VERSION))
        $(error Makefile.auto has detected an outdated MAKE/Makefile.user${SUBTARGET}. Please remove [backup] this file and run "make auto" again to generate a new one.)
    endif
endif

# OS detection
ifeq ($(OS),Windows_NT)
    MINGW := 1
else
    MINGW := 0
endif

ECHO := echo
TMPFILE = /tmp/liggghtscompile.tmp
TMP = $(shell $(ECHO) "\n" > $(TMPFILE) && cat $(TMPFILE))
EMPTY =
ifneq ($(TMP), $(EMPTY))
    ECHO := $(ECHO) -e
    TMP = $(shell $(ECHO) "\n" > $(TMPFILE) && cat $(TMPFILE))
    ifneq ($(TMP), $(EMPTY))
        $(error Could not create empty file with \n entry using echo and echo -e)
    endif
endif

# Default settings should ../Makefile.user${SUBTARGET} not exist for some reason or not contain the minimum
# variables
USE_MPI ?= "ON"
USE_VTK ?= "ON"
USE_SUPERQUADRICS ?= "OFF"
USE_CATALYST ?= "OFF"
USE_JPG ?= "OFF"
USE_FPIC ?= "ON"
USE_DEBUG ?= "OFF"
USE_PROFILE ?= "OFF"
USE_CONVEX ?= "OFF"
USE_GZIP ?= "OFF"
USE_XDR ?= "OFF"
AUTOINSTALL_CONVEX ?= "OFF"

# Default variables
ifeq ($(MINGW), 0)
    # for Linux
    MPI_CCFLAGS_USR ?=-funroll-loops -fstrict-aliasing -Wall -Wno-unused-result
    MPI_LDFLAGS_USR ?=-lstdc++
    MPI_DEPFLAGS_USR ?=-M
    CCFLAGS_USR ?=-funroll-loops -fstrict-aliasing -Wall -Wno-unused-result
    LDFLAGS_USR ?=-lstdc++
    DEPFLAGS_USR ?=-M
    ARCHIVE_USR ?= ar
    ARFLAGS_USR ?= -rcsv
    SHLIBFLAGS_USR ?= -shared
    SIZE_USR ?= size
else
    # for Windows
    ifeq ($(USE_VTK), "ON") # disable VTK
        USE_VTK = "OFF"
        $(info Warning: VTK has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_SUPERQUADRICS), "ON") # disable SUPERQUADRICS
        USE_SUPERQUADRICS = "OFF"
        $(info Warning: SUPERQUADRICS has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_CATALYST), "ON") # disable CATALYST
        USE_CATALYST = "OFF"
        $(info Warning: VTK has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_MPI), "ON") # disable MPI
        USE_MPI = "OFF"
        $(info Warning: MPI has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_JPG), "ON") # disable JPG
        USE_JPG = "OFF"
        $(info Warning: JPG has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_CONVEX), "ON") # disable CONVEX
        USE_CONVEX = "OFF"
        $(info Warning: CONVEX has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_DEBUG), "ON") # disable DEBUG
        USE_DEBUG = "OFF"
        $(info Warning: DEBUG has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_PROFILE), "ON") # disable PROFILE
        USE_PROFILE = "OFF"
        $(info Warning: DEBUG has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_GZIP), "ON") # disable GZIP
        USE_GZIP = "OFF"
        $(info Warning: GZIP has been disabled due to compilation on Windows)
    endif
    ifeq ($(USE_XDR), "OFF") # enable XDR
        USE_XDR = "ON"
        $(info Warning: XDR has been enabled due to compilation on Windows)
    endif
    CCFLAGS_USR ?=-I../STUBS -O3 -march=i686 -mtune=generic -mfpmath=387 -mpc64 -fno-exceptions -fno-rtti -ffast-math -funroll-loops -fstrict-aliasing -Wall -W -Wno-uninitialized
    LDFLAGS_USR ?=-O -lwsock32
    DEPFLAGS_USR ?=-M
    ARCHIVE_USR ?= ar
    ARFLAGS_USR ?= -rcsv
    SHLIBFLAGS_USR ?=
    SIZE_USR ?= i686-pc-mingw32-size
    OBJ += mpi.o
endif

ifeq ($(USE_CATALYST), "ON")
    ifneq ($(USE_VTK), "ON")
        USE_VTK = "ON"
        $(info Information: USE_VTK was not activated but USE_CATALYST overrode this setting)
    endif
endif

# Read previously used variables
OPTIONS_FILE=../MAKE/Makefile.auto.options${SUBTARGET}
sinclude $(OPTIONS_FILE)
# Set them to be equal the current variables if they have not been set
OPTIONS_USE_MPI           ?= $(USE_MPI)
OPTIONS_USE_VTK           ?= $(USE_VTK)
OPTIONS_USE_SUPERQUADRICS ?= $(USE_SUPERQUADRICS)
OPTIONS_USE_CATALYST      ?= $(USE_CATALYST)
OPTIONS_USE_JPG           ?= $(USE_JPG)
OPTIONS_USE_FPIC          ?= $(USE_FPIC)
OPTIONS_USE_DEBUG         ?= $(USE_DEBUG)
OPTIONS_USE_PROFILE       ?= $(USE_PROFILE)
OPTIONS_USE_CONVEX        ?= $(USE_CONVEX)
OPTIONS_USE_GZIP          ?= $(USE_GZIP)
OPTIONS_USE_XDR           ?= $(USE_XDR)
# Check if they have changed
ifneq ($(OPTIONS_USE_MPI), $(USE_MPI))
    $(error variable USE_MPI has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_VTK), $(USE_VTK))
    $(error variable USE_VTK has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_SUPERQUADRICS), $(USE_SUPERQUADRICS))
    $(error variable USE_SUPERQUADRICS has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_CATALYST), $(USE_CATALYST))
    $(error variable USE_CATALYST has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_JPG), $(USE_JPG))
    $(error variable USE_JPG has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_FPIC), $(USE_FPIC))
    $(error variable USE_FPIC has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_DEBUG), $(USE_DEBUG))
    $(error variable USE_DEBUG has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_PROFILE), $(USE_PROFILE))
    $(error variable USE_PROFILE has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_CONVEX), $(USE_CONVEX))
    $(error variable USE_CONVEX has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_GZIP), $(USE_GZIP))
    $(error variable USE_GZIP has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
ifneq ($(OPTIONS_USE_XDR), $(USE_XDR))
    $(error variable USE_XDR has changed. Please run "make clean-all" or "make clean-auto" before compiling again.)
endif
# Write out current variables
$(shell echo '# Automatically generated file from Makefile.auto.\n' > $(OPTIONS_FILE))
$(shell echo '# DO NO EDIT THIS FILE\n' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_MPI = $(USE_MPI)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_VTK = $(USE_VTK)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_SUPERQUADRICS = $(USE_SUPERQUADRICS)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_CATALYST = $(USE_CATALYST)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_JPG = $(USE_JPG)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_FPIC = $(USE_FPIC)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_PROFILE = $(USE_PROFILE)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_CONVEX = $(USE_CONVEX)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_GZIP = $(USE_GZIP)' >> $(OPTIONS_FILE))
$(shell echo 'OPTIONS_USE_XDR = $(USE_XDR)' >> $(OPTIONS_FILE))

# Sources
LIBCCD_GIT ?= https://github.com/danfis/libccd.git
LIBCCD_ZIP ?= https://github.com/danfis/libccd/archive/master.zip

# Lib path in which 3rd party libraries are installed
LIB_PATH ?= ../../lib

# Initialization of EXTRA_* variables
# All -I include paths
EXTRA_INC=
# All -L library paths
EXTRA_LIB=
# All -l libraries
EXTRA_ADDLIBS=

# Debug settings
#
ifeq ($(USE_DEBUG), "ON")
    OPT_LVL = -Og -g
    LMP_INC += -DLIGGGHTS_DEBUG
else
    ifeq ($(USE_DEBUG), "NOOPT")
        OPT_LVL = -O0 -g
        LMP_INC += -DLIGGGHTS_DEBUG
    else
        ifeq ($(USE_DEBUG), "OFF")
            ifeq ($(MINGW), 0)
                OPT_LVL = -O2
            else
                OPT_LVL = # for mingw the optimization is already set in CCFLAGS_USR
            endif
        else
            $(error 'USE_DEBUG needs to be either "ON", "NOOPT" or "OFF"')
        endif
    endif
endif

# Profiling settings
#
ifeq ($(USE_PROFILE), "ON")
    PROF_FLAG = -pg
else
    PROF_FLAG =
endif

# MPI settings
#
ifeq ($(USE_MPI), "ON")
    # Test if MPICXX is set and if yes, then check if this executable exists
    ifneq (, $(MPICXX))
        ifeq (, $(shell which $(MPICXX)))
            $(error "Could not find executable set in variable MPICXX")
        endif
    # MPICXX is not set
    else
        # Check if MPICXX_USR is set
        ifneq (, $(MPICXX_USR))
            ifeq (, $(shell which $(MPICXX_USR)))
                $(error "Could not find executable set in MPICXX_USR")
            else
                # found it, so set MPICXX accordingly
                MPICXX = $(MPICXX_USR)
            endif
        # if MPICXX_USR is not set then try to find mpicxx or mpic++
        else
            ifeq (, $(shell which mpicxx))
                ifeq (, $(shell which mpic++))
                    $(error 'Could not find a suitable mpi compiler (mpicxx or mpic++). Please specify one in the Makefile.user${SUBTARGET} using the MPICXX_USR variable or set USE_MPI to "OFF"')
                else
                    MPICXX = mpic++
                endif
            else
                MPICXX = mpicxx
            endif
        endif
    endif
    # Detect openmpi mvapich mpich2 etc.
    ifneq ($(origin MPI_INC_USR), undefined)
        MPI_INC ?= $(MPI_INC_USR)
    endif
    TMP_INC =
    ifneq ($(origin MPI_INC), undefined)
        TMP_INC = -I$(MPI_INC)
    endif
    # We assume that the compiler supports #pragma message
    TMP := $(shell $(ECHO) '\#include <mpi.h> \n \#if defined(MPICH) \n \#pragma message "MPICH" \n \#elif defined(OPEN_MPI) \n \#pragma message "OpenMPI" \n \#else \n \#pragma message "Unknown" \n \#endif' > $(TMPFILE) && $(MPICXX) $(OPT_LVL) $(PROF_FLAG) $(TMP_INC) -xc++ -E $(TMPFILE) 2> /dev/null | grep pragma | grep -m 1 message || echo -1)
    # See if compilation has worked out
    ifeq ($(TMP), -1)
        # Maybe it failed because of the optimization as -Og is not known
        ifeq ($(USE_DEBUG), "ON")
            TMP := $(shell $(ECHO) '\#include <mpi.h> \n \#if defined(MPICH) \n \#pragma message "MPICH" \n \#elif defined(OPEN_MPI) \n \#pragma message "OpenMPI" \n \#else \n \#pragma message "Unknown" \n \#endif' > $(TMPFILE) && $(MPICXX) -O0 -g $(PROF_FLAG) $(TMP_INC) -xc++ -E $(TMPFILE) 2> /dev/null | grep pragma | grep -m 1 message || echo -1)
            ifeq ($(TMP), -1)
                $(error 'Could not compile a simple MPI example (testing with -Og and -O0). Test was done with MPI_INC="$(TMP_INC)" and MPICXX="$(MPICXX)"')
            else
                OPT_LVL = -O0 -g
                $(info 'Could not use -Og compilation flag, replacing with -O0')
            endif
        else
            $(error 'Could not compile a simple MPI example. Test was done with MPI_INC="$(TMP_INC)" and MPICXX="$(MPICXX)"')
        endif
    endif
    MPI_TYPE := $(word $(words $(TMP)), $(TMP))
    ifeq ($(MPI_TYPE), MPICH)
        ifneq ($(origin MPI_ADDLIBS), undefined)
            MPI_ADDLIBS += -lmpich -lmpichcxx
        else
            ifneq ($(origin MPI_ADDLIBS_USR), undefined)
                MPI_ADDLIBS_USR += -lmpich -lmpichcxx
            else
                MPI_ADDLIBS = -lmpich -lmpichcxx
            endif
        endif
    else
        ifeq ($(MPI_TYPE), OpenMPI)
            # Here we could do some OpenMPI related stuff
        endif
    endif

    # Now that we have found an mpi compiler check that we can actually do something with it
    # 0.) Prepare the environment

    ifneq ($(origin MPI_LIB_USR), undefined)
        MPI_LIB ?= $(MPI_LIB_USR)
    endif
    MPI_ADDLIBS ?= $(MPI_ADDLIBS_USR)
    MPI_CCFLAGS ?= $(MPI_CCFLAGS_USR)
    MPI_LDFLAGS ?= $(MPI_LDFLAGS_USR)
    MPI_CCFLAGS += $(OPT_LVL) $(PROF_FLAG)
    MPI_LDFLAGS += $(OPT_LVL) $(PROF_FLAG)
    # Linking test
    TMP := $(shell $(ECHO) 'main(){}' > $(TMPFILE) && $(MPICXX) -xc++ $(MPI_LDFLAGS) -L$(MPI_LIB) $(MPI_ADDLIBS) -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile and link a simple MPI example. Test was done with MPI_LIB=$(MPI_LIB)')
    endif
    # Finally set the compiler environment
    CXX = $(MPICXX)
    CCFLAGS = $(MPI_CCFLAGS)
    LDFLAGS = $(MPI_LDFLAGS)
    DEPFLAGS ?= $(MPI_DEPFLAGS_USR)
    ifneq ($(origin MPI_INC), undefined)
        EXTRA_INC += -I$(MPI_INC)
    endif
    ifneq ($(origin MPI_LIB), undefined)
        EXTRA_LIB += -L$(MPI_LIB)
    endif
    EXTRA_ADDLIBS += $(MPI_ADDLIBS)

# If USE_MPI is not set, check that we have a C++ compiler available
else
    # Test if MPICXX is set and if yes, then check if this executable exists
    ifneq (, $(CXX))
        ifeq (, $(shell which $(CXX)))
            $(error "Could not find executable set in variable CXX")
        endif
    # CXX is not set
    else
        # Check if CXX_USR is set
        ifneq (, $(CXX_USR))
            ifeq (, $(shell which $(CXX_USR)))
                $(error "Could not find executable set in CXX_USR")
            else
                # found it, so set MPICXX accordingly
                CXX = $(CXX_USR)
            endif
        else
            # if CXX is not set then try to find g++ on linux
            ifeq ($(MINGW), 0)
                ifeq (, $(shell which g++))
                    $(error 'Could not find a suitable c++ compiler (g++). Please specify one in the Makefile.user${SUBTARGET} using the CXX_USR variable')
                else
                    CXX = g++
                endif
            else # or the appropriate mingw compiler
                ifeq (, $(shell which i686-pc-mingw32-g++))
                    $(error 'Could not find a suitable c++ compiler (i686-pc-mingw32-g++). Please specify one in the Makefile.user${SUBTARGET} using the CXX_USR variable')
                else
                    CXX = i686-pc-mingw32-g++
                endif
            endif
        endif
    endif

    # Now that we have found a c++ compiler check that we can actually do something with it
    # 0.) Prepare the environment
    CCFLAGS ?= $(CCFLAGS_USR)
    LDFLAGS ?= $(LDFLAGS_USR)
    DEPFLAGS ?= $(DEPFLAGS_USR)
    CCFLAGS += $(OPT_LVL)
    LDFLAGS += $(OPT_LVL)
    EXTRA_INC += -I../STUBS
    EXTRA_LIB += -L../STUBS
    EXTRA_ADDLIBS += -lmpi_stubs
    # 1.) Attempt compilation
    TMP := $(shell $(ECHO) '\#include <mpi.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_INC) $(EXTRA_LIB) $(EXTRA_ADDLIBS) -xc++ $(LDFLAGS) $(CCFLAGS) -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile a simple c++ example. Please make sure that you have run "make stubs" before compiling LIGGGHTS itself. Test was done with CXX=$(CXX), EXTRA_INC=$(EXTRA_INC), EXTRA_LIB=$(EXTRA_LIB) and EXTRA_ADDLIBS=$(EXTRA_ADDLIBS).')
    endif
endif

# Check for c++11 support
C11SUPPORT = 0
TMP := $(shell $(ECHO) 'main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_INC) $(EXTRA_LIB) $(EXTRA_ADDLIBS) -xc++ $(LDFLAGS) $(CCFLAGS) -std=c++11 -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
ifeq ($(TMP),0)
    C11SUPPORT = 1
endif

# Position independent code
#
ifeq ($(USE_FPIC), "ON")
    SHFLAGS = -fPIC
    CCFLAGS += -fPIC
else
    ifeq ($(USE_FPIC), "OFF")
        SHFLAGS =
    else
        $(error 'USE_FPIC needs to be either "ON" or "OFF"')
    endif
endif

# Check for gzip
#
ifeq ($(USE_GZIP), "ON")
    ifneq (, $(shell which gzip))
        LMP_INC += -DLAMMPS_GZIP
    else
        $(error 'USE_GZIP is set to "ON" but the gzip executable cannot be found. Ensure that it is available in your PATH variable or set USE_GZIP to "OFF"')
    endif
endif

ifeq ($(USE_XDR), "ON")
    LMP_INC += -DLAMMPS_XDR
endif

# Default for building shared libraries
SHLIBFLAGS = -shared

# Library compilation flags
ARCHIVE ?= $(ARCHIVE_USR)
ARFLAGS ?= $(ARFLAGS_USR)
SHFLAGS ?= $(SHFLAGS_USR)
SHLIBFLAGS ?= $(SHLIBFLAGS_USR)

# Size
SIZE ?= $(SIZE_USR)

# Linker (for now equal to CXX) if not set explicitly
LINK ?= $(CXX)

# CATALYST for in-situ ParaView visualization
#
ifeq ($(USE_CATALYST), "ON")
    TMP := $(shell cp ../MAKE/CMakeLists_Catalyst.txt CMakeLists.txt && mkdir -p .tmp_catalyst_build && cd .tmp_catalyst_build && cd - && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not create .tmp_catalyst_build and create CMakeLists.txt')
    endif
    # All vtk includes in LIGGGHTS (to be inserted into catalyst_dummy.cpp)
    #grep -rh "#include" *{h,cpp} | grep "<vtk"
    ifeq (, $(shell which cmake))
        $(error 'Could not find cmake. Please ensure that it is in your PATH environment variable')
    endif
    PV_BUILD_DIR ?= $(PV_BUILD_DIR_USR)
    TMP := $(shell cd .tmp_catalyst_build && cmake -DParaView_DIR="$(PV_BUILD_DIR)" .. 1> /dev/null && cd - && echo 0 || echo -1)
    TMP2 := $(shell rm -rf CMakeLists.txt CMakeCache.txt CMakeFiles)
    ifeq ($(TMP), -1)
        $(error 'Could not generate makefiles using cmake. Make sure PV_BUILD_DIR="$(PV_BUILD_DIR)" is set correctly')
    endif
    # TODO check if CXX_INCLUDES is set, if yes save it somewhere temporarily (same with other flags in flags.make)
    sinclude .tmp_catalyst_build/CMakeFiles/CatalystAdaptor.dir/flags.make
    CATALYST_INC = $(CXX_INCLUDES)
    # Next the link_raw needs to remove the dummy crap and put this into the normal linking flags
    CATALYST_LIB := $(shell cat .tmp_catalyst_build/CMakeFiles/CatalystDummy.dir/link.txt | sed s/.*libCatalystAdaptor.a// || echo -1)
    ifeq ($(CATALYST_LIB), -1)
        $(error 'Could not get linking flags from cmake files')
    endif
    # prepare the environment
    ifneq ($(origin VTK_INC_USR), undefined)
        $(info Information: VTK_INC_USR was set, but USE_CATALYST overrides this)
    endif
    ifneq ($(origin VTK_LIB_USR), undefined)
        $(info Information: VTK_LIB_USR was set, but USE_CATALYST overrides this)
    endif
    VTK_INC_USR = $(CATALYST_INC)
    VTK_LIB_USR = $(CATALYST_LIB)
    # Set the USE_CATALYST_FLAG
    LMP_INC += -DUSE_CATALYST_FLAG
endif

# VTK library, OPTIONAL
# INC = path for VTK header files
# LIB = path for VTK library
# ADDLIBS = name of VTK library

ifeq ($(USE_VTK), "ON")
    # Use user settings if available
    # The following ifneq are required as otherwise VTK_{INC,LIB} will be set (but are empty) if VTK_{INC,LIB}_USR is undefined. This will conflict with the following ?= assignments further below.
    ifneq ($(origin VTK_INC_USR), undefined)
        VTK_INC ?= $(VTK_INC_USR)
    endif
    ifneq ($(origin VTK_LIB_USR), undefined)
        VTK_LIB ?= $(VTK_LIB_USR)
    endif
    VTK_ADDLIBS ?= $(VTK_ADDLIBS_USR)
    VTK_INC ?= -I$(shell ls -d /usr/include/vtk* 2> /dev/null | tail -n 1)
    ifeq ($(VTK_INC),-I)
        VTK_INC =
    endif
    # note we assume here that our compiler supports #pragma message
    VTK_TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n \#define XSTR(x) STR(x) \n \#define STR(x) \#x \n \#pragma message XSTR(VTK_MAJOR_VERSION)' > $(TMPFILE) && $(CXX) -Wno-deprecated -E $(VTK_INC) -xc++ $(TMPFILE) 2> /dev/null | grep "pragma" | grep "message" || echo -1)
    ifeq ($(VTK_TMP), -1)
        $(error Could not obtain VTK_MAJOR_VERSION)
    endif
    VTK_MAJOR_VERSION := $(patsubst "%",%,$(word $(words $(VTK_TMP)),$(VTK_TMP)))
    VTK_TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n \#define XSTR(x) STR(x) \n \#define STR(x) \#x \n \#pragma message XSTR(VTK_MINOR_VERSION)' > $(TMPFILE) && $(CXX) -Wno-deprecated -E $(VTK_INC) -xc++ $(TMPFILE) 2> /dev/null | grep "pragma" | grep "message" || echo -1)
    ifeq ($(VTK_TMP), -1)
        $(error Could not obtain VTK_MINOR_VERSION)
    endif
    VTK_MINOR_VERSION := $(patsubst "%",%,$(word $(words $(VTK_TMP)),$(VTK_TMP)))
    VTK_VERSION=$(VTK_MAJOR_VERSION).$(VTK_MINOR_VERSION)

    # check if the vtk compiled with ParaView is going to be used
    USE_PV_VTK = "OFF"
    ifeq ($(USE_CATALYST), "ON")
        # VTK_LIB is set already
        FILTER = $(foreach v,$(2),$(if $(findstring $(1),$(v)),$(v),))
        ifeq ($(VTK_MAJOR_VERSION),5)
            VTK_APPENDIX = $(strip $(call FILTER,libvtkCommon, $(VTK_LIB)))
        else
            ifneq (,$(filter $(VTK_MAJOR_VERSION),6 7 8))
                VTK_APPENDIX = $(strip $(call FILTER,libvtkParallelMPI, $(VTK_LIB)))
            else
                $(error 'Unsupported VTK version $(VTK_VERSION)')
            endif
        endif
        USE_PV_VTK = $(findstring -pv, $(VTK_APPENDIX))
        ifeq ($(USE_PV_VTK), -pv)
            USE_PV_VTK = "ON"
            # strip the path, library name and ending
            VTK_APPENDIX := $(subst libvtkParallelMPI-,,$(patsubst %.so,%,$(basename $(notdir $(VTK_APPENDIX)))))
            VTK_LIB += -L$(PV_BUILD_DIR)/lib
        endif
    endif

    ifeq ($(USE_PV_VTK), "OFF") # uses system vtk
        ifeq ($(origin VTK_APPENDIX_USR), undefined)
            VTK_APPENDIX = $(VTK_VERSION)
        else
            VTK_APPENDIX = $(VTK_APPENDIX_USR)
        endif
    endif

    ifeq ($(VTK_MAJOR_VERSION),5)
        ifneq ($(origin VTK_APPENDIX_USR), undefined)
            VTK_APPENDIX_5 = -$(VTK_APPENDIX)
        endif
        ifeq ($(USE_PV_VTK), "ON")
            VTK_APPENDIX_5 = -$(VTK_APPENDIX)
        endif
        VTK_LIB_SET = 0
        ifneq ($(origin VTK_LIB), undefined)
            VTK_LIB_SET = 1
        else
            VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkCommon$(VTK_APPENDIX_5).so' | tail -n 1))
            ifeq ($(VTK_LIB),-L)
                VTK_LIB =
            endif
        endif
        TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(VTK_LIB) $(VTK_INC) $(EXTRA_ADDLIBS) -lvtkCommon$(VTK_APPENDIX_5) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
        ifeq ($(TMP), -1)
            ifeq ($(VTK_LIB_SET), 0)
                VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkCommon.so' | tail -n 1))
                ifeq ($(VTK_LIB),-L)
                    VTK_LIB =
                endif
            endif
            TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(VTK_LIB) $(VTK_INC) $(EXTRA_ADDLIBS) -lvtkCommon $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
            ifeq ($(TMP), -1)
                $(error 'Could not determine suitable appendix of VTK library with VTK_INC="$(VTK_INC)", VTK_LIB="$(VTK_LIB)" and VTK_APPENDIX="$(VTK_APPENDIX)"')
            else
                # everything worked, so no appendix needed
                VTK_APPENDIX_5 = $(EMPTY)
            endif
        endif
        VTK_ADDLIBS +=  -lvtkCommon$(VTK_APPENDIX_5) \
                        -lvtkFiltering$(VTK_APPENDIX_5) \
                        -lvtkIO$(VTK_APPENDIX_5) \
                        -lvtkParallel$(VTK_APPENDIX_5) \
                        -lvtkGraphics$(VTK_APPENDIX_5)
    else
        VTK_APPENDIX := -$(VTK_APPENDIX)
        VTK_LIB_SET = 0
        ifneq ($(origin VTK_LIB), undefined)
            VTK_LIB_SET = 1
        else
            VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkParallelMPI$(VTK_APPENDIX).so' 2> /dev/null | tail -n 1))
            ifeq ($(VTK_LIB),-L)
                VTK_LIB =
            endif
        endif
        TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(VTK_LIB) $(VTK_INC) $(EXTRA_ADDLIBS) -lvtksys$(VTK_APPENDIX) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
        ifeq ($(TMP), -1)
            ifeq ($(VTK_LIB_SET), 0)
                VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkParallelMPI.so' 2> /dev/null | tail -n 1))
                ifeq ($(VTK_LIB),-L)
                    VTK_LIB =
                endif
            endif
            TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(VTK_LIB) $(VTK_INC) $(EXTRA_ADDLIBS) -lvtksys $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
            ifeq ($(TMP), -1)
                $(error 'Could not determine suitable appendix of VTK library with VTK_INC="$(VTK_INC)", VTK_LIB="$(VTK_LIB)" and VTK_APPENDIX="$(VTK_APPENDIX)"')
            else
                # everything worked, so no appendix needed
                VTK_APPENDIX = $(EMPTY)
            endif
        endif
        VTKSET = 0
        ifeq ($(VTK_MAJOR_VERSION),6)
            VTKSET = 1
            VTK_ADDLIBS +=  -lvtkCommonCore$(VTK_APPENDIX) \
                            -lvtkIOCore$(VTK_APPENDIX) \
                            -lvtkIOXML$(VTK_APPENDIX) \
                            -lvtkIOLegacy$(VTK_APPENDIX) \
                            -lvtkIOImage$(VTK_APPENDIX) \
                            -lvtkCommonDataModel$(VTK_APPENDIX) \
                            -lvtkParallelCore$(VTK_APPENDIX) \
                            -lvtkParallelMPI$(VTK_APPENDIX) \
                            -lvtkCommonExecutionModel$(VTK_APPENDIX) \
                            -lvtkFiltersCore$(VTK_APPENDIX) \
                            -lvtksys$(VTK_APPENDIX) \
                            -lvtkCommonMisc$(VTK_APPENDIX) \
                            -lvtkCommonTransforms$(VTK_APPENDIX) \
                            -lvtkCommonMath$(VTK_APPENDIX) \
                            -lvtkIOXMLParser$(VTK_APPENDIX) \
                            -lvtkCommonSystem$(VTK_APPENDIX) \
                            -lvtkDICOMParser$(VTK_APPENDIX) \
                            -lvtkmetaio$(VTK_APPENDIX)
            # For VTK >= 6.2 we need to use lvtkIOParallelXML
            MINOR_GE_2 := $(shell [ $(VTK_MINOR_VERSION) -ge 2 ] && echo true)
            ifeq ($(MINOR_GE_2),true)
                VTK_ADDLIBS += -lvtkIOParallelXML$(VTK_APPENDIX)
            else
                VTK_ADDLIBS += -lvtkIOParallel$(VTK_APPENDIX)
            endif
            ifeq ($(VTK_LIB_SET), 0)
                VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkParallelMPI$(VTK_APPENDIX).so' 2> /dev/null | tail -n 1))
                ifeq ($(VTK_LIB),-L)
                    VTK_LIB =
                endif
            endif
        endif
        ifneq (, $(filter $(VTK_MAJOR_VERSION),7 8))
            VTKSET = 1
            VTK_ADDLIBS +=  -lvtkCommonCore$(VTK_APPENDIX) \
                            -lvtkIOCore$(VTK_APPENDIX) \
                            -lvtkIOXML$(VTK_APPENDIX) \
                            -lvtkIOLegacy$(VTK_APPENDIX) \
                            -lvtkIOImage$(VTK_APPENDIX) \
                            -lvtkCommonDataModel$(VTK_APPENDIX) \
                            -lvtkIOParallelXML$(VTK_APPENDIX) \
                            -lvtkParallelCore$(VTK_APPENDIX) \
                            -lvtkParallelMPI$(VTK_APPENDIX) \
                            -lvtkCommonExecutionModel$(VTK_APPENDIX) \
                            -lvtkFiltersCore$(VTK_APPENDIX) \
                            -lvtksys$(VTK_APPENDIX) \
                            -lvtkCommonMisc$(VTK_APPENDIX) \
                            -lvtkCommonTransforms$(VTK_APPENDIX) \
                            -lvtkCommonMath$(VTK_APPENDIX) \
                            -lvtkIOXMLParser$(VTK_APPENDIX) \
                            -lvtkCommonSystem$(VTK_APPENDIX) \
                            -lvtkDICOMParser$(VTK_APPENDIX) \
                            -lvtkmetaio$(VTK_APPENDIX)
            ifeq ($(VTK_LIB_SET), 0)
                VTK_LIB := -L$(dir $(shell find /usr/lib* -name 'libvtkParallelMPI$(VTK_APPENDIX).so' 2> /dev/null | tail -n 1))
                ifeq ($(VTK_LIB),-L)
                    VTK_LIB =
                endif
            endif
            ifeq ($(VTK_MAJOR_VERSION),8)
                ifeq ($(C11SUPPORT),1)
                    CCFLAGS += -std=c++11
                else
                    $(error VTK version $(VTK_VERSION) requires support for -std=c++11 in the compiler flags)
                endif
            endif
        endif
        ifeq ($(VTKSET),0)
            TMP:=$(error Unsupported VTK version $(VTK_VERSION))
        endif
    endif

    # Test the settings
    TMP := $(shell $(ECHO) '\#include <vtkVersion.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(VTK_LIB) $(VTK_INC) $(EXTRA_ADDLIBS) $(VTK_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile VTK example with VTK_INC="$(VTK_INC)", VTK_LIB="$(VTK_LIB)" and VTK_ADDLIBS="$(VTK_ADDLIBS)"')
    endif
    # Prepare environment
    LMP_INC += -DLAMMPS_VTK
    EXTRA_ADDLIBS += $(VTK_ADDLIBS)
    EXTRA_LIB += $(VTK_LIB)
    EXTRA_INC += $(VTK_INC)

    # Write out VTK Data so that it can be used by CFDEMcoupling without detecting it again
    $(shell echo 'OPTIONS_VTK_ADDLIBS = $(VTK_ADDLIBS)' >> $(OPTIONS_FILE))
    $(shell echo 'OPTIONS_VTK_LIB = $(VTK_LIB)' >> $(OPTIONS_FILE))
    $(shell echo 'OPTIONS_VTK_INC = $(VTK_INC)' >> $(OPTIONS_FILE))
endif

# Superquadric particles
#
ifeq ($(USE_SUPERQUADRICS), "ON")
    BOOST_INC ?= $(BOOST_INC_USR)
    # Include test
    TMP := $(shell $(ECHO) '\#include "boost/math/special_functions/beta.hpp" \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(BOOST_INC) $(EXTRA_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile boost example with BOOST_INC="$(BOOST_INC)" as boost/math/special_functions/beta.hpp could not be found')
    endif
    # Prepare environment
    LMP_INC += -DSUPERQUADRIC_ACTIVE_FLAG -DNONSPHERICAL_ACTIVE_FLAG
    EXTRA_INC += $(BOOST_INC)
endif

# JPG output for sphere renderings
#
ifeq ($(USE_JPG), "ON")
    JPG_INC ?= $(JPG_INC_USR)
    JPG_LIB ?= $(JPG_LIB_USR)
    JPG_ADDLIBS += -ljpeg
    # Test settings
    # Include test
    TMP := $(shell $(ECHO) '\#include <libjpeg.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(JPG_INC) $(EXTRA_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile jpg example with JPG_INC="$(JPG_INC)"')
    endif
    # Linking test
    TMP := $(shell $(ECHO) '\#include <libjpeg.h> \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(JPG_LIB) $(JPG_INC) $(EXTRA_ADDLIBS) $(JPG_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile and link jpg example with JPG_INC="$(JPG_INC)", JPG_LIB="$(JPG_LIB)" and JPG_ADDLIBS="$(JPG_ADDLIBS)"')
    endif
    # Prepare environment
    LMP_INC += -DLAMMPS_JPEG
    EXTRA_INC += $(JPG_INC)
    EXTRA_LIB += $(JPG_LIB)
    EXTRA_ADDLIBS += $(JPG_ADDLIBS)
endif

# Polyhedral particles (premium only)
#
ifeq ($(USE_CONVEX), "ON")
    ifneq ($(origin CONVEX_INC_USR), undefined)
        CONVEX_INC ?= $(CONVEX_INC_USR)
    endif
    ifneq ($(origin CONVEX_LIB_USR), undefined)
        CONVEX_LIB ?= $(CONVEX_LIB_USR)
    endif
    # Defaults
    CONVEX_INC ?= -I$(LIB_PATH)/libccd/src
    CONVEX_LIB ?= -L$(LIB_PATH)/libccd/src
    # This should really not be required to be set
    CONVEX_ADDLIBS += $(CONVEX_ADDLIBS_USR)
    # This is the ccd library
    CONVEX_ADDLIBS += -lccd
    # Test settings
    # Link test
    $(info $(CXX) $(EXTRA_LIB) $(CONVEX_LIB) $(CONVEX_INC) $(EXTRA_ADDLIBS) $(CONVEX_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE))
    TMP := $(shell $(ECHO) '\#include "ccd/ccd.h" \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(CONVEX_LIB) $(CONVEX_INC) $(EXTRA_ADDLIBS) $(CONVEX_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    # Automatic download and compilation if AUTODOWNLOAD_CONVEX is set
    ifeq ($(TMP), -1)
        ifeq ($(AUTOINSTALL_CONVEX), "ON")
            # Check if sources provide the convex surface model
            ifneq ($(wildcard ../surface_model_convex.h),)
                $(error 'It appears that you do not posses the source files for the convex hull surface model. Note that this is available only in the premium version of LIGGGHTS.')
            endif
            $(info 'Could not compile Convex (libccd) example. As AUTOINSTALL_CONVEX is set to "ON". libccd will now be automatically downloaded to ../lib/libccd')
            # Check if $(LIB_PATH)/libccd/src folder is available
            TMP := $(shell ls $(LIB_PATH)/libccd/src && echo 0 || echo -1)
            # If not download libccd
            ifeq ($(TMP), -1)
                # Create $(LIB_PATH) folder
                TMP := $(shell mkdir -p $(LIB_PATH) && echo 0 || echo -1)
                ifeq ($(TMP), -1)
                    $(error 'Could not create $(LIB_PATH) directory')
                endif
                # Check for git
                ifneq (, $(shell which git))
                    # Download git repository
                    TMP := $(shell git clone $(LIBCCD_GIT) $(LIB_PATH)/libccd && echo 0 || echo -1)
                    ifeq ($(TMP), -1)
                        $(error 'Could not perform "git clone" of $(LIBCCD_GIT) into $(LIB_PATH)/libccd')
                    endif
                # Git not available
                else
                    # Check for wget
                    ifeq (,$(shell which wget))
                        $(error "Could not find either git or wget to download libccd. Please install either of the two in order to allow the auto-installation of libccd.")
                    endif
                    # Download zip file
                    TMP := $(shell wget $(LIBCCD_ZIP) -P $(LIB_PATH) &> /dev/null && echo 0 || echo -1)
                    ifeq ($(TMP), -1)
                        $(error 'Failed to download libccd using wget. Please install libccd manually')
                    endif
                    # Check if unzip is available
                    ifeq (, $(shell which unzip))
                        $(error 'Could not find unzip. Please install it to allow the auto-installation of libccd.')
                    endif
                    # Unzip file
                    TMP := $(shell unzip $(LIB_PATH)/$(notdir $(LIBCCD_ZIP)) -d $(LIB_PATH) && echo 0 || echo -1)
                    ifeq ($(TMP), -1)
                        $(error 'Could not unzip $(notdir $(LIBCCD_ZIP)) in $(LIB_PATH)')
                    endif
                endif
            endif
            # At this stage we now have libccd downloaded. Next we need to compile it
            OBJDIR := $(PWD)
            TMP := $(shell $(ECHO) '\#!/bin/bash \n cd $(OBJDIR)/$(LIB_PATH)/libccd/src \n make PREFIX=$(PWD)/../../ USE_DOUBLE=yes &> /dev/null' > $(TMPFILE))
            TMP := $(shell bash $(TMPFILE) && echo 0 || echo -1)
            ifeq ($(TMP), -1)
                $(error 'Compilation of libccd failed. Please install it manually')
            endif
        else
            $(error 'Could not compile Convex (libccd) example with CONVEX_INC="$(CONVEX_INC)"')
        endif
    endif
    # Include test
    TMP := $(shell $(ECHO) '\#include "ccd/ccd.h" \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(CONVEX_INC) $(EXTRA_ADDLIBS) $(CCFLAGS) -xc++ -E $(TMPFILE) &> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile Convex (libccd) example with CONVEX_INC="$(CONVEX_INC)"')
    endif
    # Link test
    TMP := $(shell $(ECHO) '\#include "ccd/ccd.h" \n main(){}' > $(TMPFILE) && $(CXX) $(EXTRA_LIB) $(CONVEX_LIB) $(CONVEX_INC) $(EXTRA_ADDLIBS) $(CONVEX_ADDLIBS) $(CCFLAGS) -xc++ -o /dev/null $(TMPFILE) 2> /dev/null && echo 0 || echo -1)
    ifeq ($(TMP), -1)
        $(error 'Could not compile and link Convex (libccd) example with CONVEX_INC="$(CONVEX_INC)", CONVEX_LIB="$(CONVEX_LIB)" and CONVEX_ADDLIBS="$(CONVEX_ADDLIBS)"')
    endif
    # Prepare environment
    LMP_INC += -DNONSPHERICAL_ACTIVE_FLAG -DCONVEX_ACTIVE_FLAG
    EXTRA_INC += $(CONVEX_INC)
    EXTRA_LIB += $(CONVEX_LIB)
    EXTRA_ADDLIBS += $(CONVEX_ADDLIBS)
endif

# Add pre-processor defines to includes
EXTRA_INC += $(LMP_INC)

# For mingw redefine the target
ifeq ($(MINGW), 1)
    TARGET = $(EXE).exe
else
    TARGET = $(EXE)
endif

$(shell rm $(TMPFILE))

# ---------------------------------------------------------------------
# build rules and dependencies
# no need to edit this section

include	Makefile.package.settings
include	Makefile.package

EXTRA_INC += $(PKG_INC) $(PKG_SYSINC)
EXTRA_LIB += $(PKG_LIB) $(PKG_SYSLIB)
EXTRA_PATH += $(PKG_PATH) $(PKG_SYSPATH)
EXTRA_ADDLIBS += $(PKG_ADDLIBS) $(PKG_SYSADDLIBS)

# Path to src files

ifeq ($(MINGW), 0)
    vpath %.cpp ..
    vpath %.h ..
endif

# Link target

ifeq ($(MINGW), 1)
$(EXE): $(TARGET) ../MAKE/mingw_cross.nsis
	makensis ../MAKE/mingw_cross.nsis
	(cd ..; zip -0 lammps-icms-win.zip lammps-icms-win.exe)
	touch $(EXE)
endif

$(TARGET):	$(OBJ)
	$(LINK) $(LDFLAGS) $(OBJ) $(EXTRA_PATH) $(EXTRA_LIB) $(EXTRA_ADDLIBS) -o $(EXE)
	$(SIZE) $(EXE)

# Library targets

lib:	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(EXE) $(OBJ)

shlib:	$(OBJ)
	$(CXX) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) $(EXTRA_PATH) -o $(EXE) \
        $(OBJ) $(EXTRA_ADDLIBS) $(EXTRA_LIB)

# Compilation rules

%.o:%.cpp
	$(CXX) $(CCFLAGS) $(SHFLAGS) $(EXTRA_PATH) $(EXTRA_INC) -c $<

%.d:%.cpp
	$(CXX) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

ifeq ($(MINGW), 1)
mpi.o:  ../STUBS/mpi.c ../STUBS/mpi.h
	$(CXX) $(CCFLAGS) $(EXTRA_INC) -c $<
endif

DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)
