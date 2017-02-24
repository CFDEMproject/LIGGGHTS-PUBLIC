//#***************************************************************************************
//# filename:     preprocessor_includes.h
//#
//# author:				Johannes Gerstmayr, Yuri Vetyukov
//# 
//# generated:      
//# description:  
//#                       
//# comments:      
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
//#***************************************************************************************

#ifndef  __PREPROCESSOR_INCLUDES_H__
#define __PREPROCESSOR_INCLUDES_H__



//EDIT THE FOLLOWING: choose right one from the list (for checkin to CVS, choose: __DEVELOPER_VERSION__):

//#define __RELEASE_VERSION__
//#define __COMPANY_VERSION__
#define __DEVELOPER_VERSION__    //(default for CVS checkin)




//GENERALLY DON'T EDIT THE FOLLOWING flags for respective version:

#ifdef __RELEASE_VERSION__
	#define __EXCLUDE_EXPERIMENTAL_OBJECTS__
#endif

#ifdef __COMPANY_VERSION__
	#define __EXCLUDE_EXPERIMENTAL_MENU_ITEMS__
	#define __EXCLUDE_EXPERIMENTAL_OBJECTS__
	// define further company specific flags
	// ... 
	// which you shall not commit to CVS!
#endif

#ifdef __DEVELOPER_VERSION__
	// define whatever you like
	// #define __ASSERT_IN_RELEASE_MODE__    (uncomment for assert information in release mode, commented by default)
	// ..
	// but comment out before CVS-commit !
#endif



#endif  /*__PREPROCESSOR_INCLUDES_H__*/