INSTRUCTIONS FOR COMPILING LIGGGHTS WITH VISUAL STUDIO 2010/2012 
(Ultimate, Professional or Express Versions)

To compile LIGGGHTS open the LIGGGHTS_VS2012 Solution

The LIGGGHTS project has configurations to compile either with MPI 
support or with MPI stubs. *

To compile with MPI:

1.  Install MPICH for Windows (http://www.mpich.org/), validate the 
    corresponding include and lib directories in the project properties 
    of LIGGGHTS: LIGGGHTS/Properties/Configuration Properties/VC++ Directories **

2.  Compile LIGGGHTS using Debug or Release configurations from the
    provided projects (use x64 for 64bit binary)

To compile with MPI STUBS
   
1.  Compile STUBS.vcproj 

2.  Compile LIGGGHTS using Debug_STUBS or Release_STUBS configurations
from the provided project (use x64 for 64bit binary)


* For Visual Studio versions prior to 2012  the Platform Toolset setting has to
be adjusted for each project. This setting can be changed under:
Properties/Configuration Properties/General/General/Platform Toolset

** Note: Depending on your MPICH Installation (32bit, 64bit) and your Windows 
Installation (32bit, 64bit) you might have to change the 32bit paths from
c:\Program Files (x86)\ to C:\Program Files\. The supplied solution was created
on a 64bit system.