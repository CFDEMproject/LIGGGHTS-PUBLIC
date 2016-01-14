//#**************************************************************
//#
//# filename:             ioincludes.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						2008
//# description:          replaces old iostreams and stdio 
//#												
//# remarks:						  
//#
//# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian
//# Center of Competence in Mechatronics GmbH, Institute of Technical Mechanics at the 
//# Johannes Kepler Universitaet Linz, Austria. All rights reserved.
//#
//# This file is part of HotInt.
//# HotInt is free software: you can redistribute it and/or modify it under the terms of 
//# the HOTINT license. See folder 'licenses' for more details.
//#
//# bug reports are welcome!!!
//# WWW:		www.hotint.org
//# email:	bug_reports@hotint.org or support@hotint.org
//#**********************************************finite_element_definition****************


#ifndef IOINCLUDES__H
#define IOINCLUDES__H

//comment the following line out for old Visual Studio.NET
#define my_new_stdiostream

//use the following flag, in case that you use the old visual studio 6
//#ifndef NEW_VISUAL_STUDIO 
//#define sprintf_s sprintf
//#endif

#ifdef my_new_stdiostream
 #include <iostream>
 #include <ostream>
 #include <istream>
 #include <fstream>
 #include <cstdlib>
 #include <cstdio>
 using namespace std;
#else
 

 #include <iostream.h>
 #include <fstream.h>
 #include <stdlib.h>
 #include <stdio.h> 
#endif



#endif
