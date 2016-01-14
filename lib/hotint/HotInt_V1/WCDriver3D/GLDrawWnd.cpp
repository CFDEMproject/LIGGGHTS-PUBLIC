//#**************************************************************
//# filename:             GLDrawWnd.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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
 


#include "stdafx.h"

//#include <iostream.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <fstream.h>

#define my_new_stdiostream

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>


//#include <math.h>

#include "savewindowbitmap.h"
#include "GLDrawWnd.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include "resource.h"

#define ROTATION_TIMER_ID 1
#define ROTATION_TIMER_DELAY 20
#define SCALING_ANIMATE_TIMER_ID 2

const int MAX_SIZE = 10000000;
const float MAX_FOVY = 160.0f;				// Maximal viewing angle

const int ToolBarButtons = 10;
const int ToolBarButtonSize = 24;

const int FrameNumberDigits = 6;


const float NoSavedModelviewMatrix = 1.8945e10f;



//   glEnable(GL_POLYGON_STIPPLE);
//   glPolygonStipple(stippleMask[0]);  /* 0% opaqueness */
//   glPolygonStipple(stippleMask[8]);  /* 50% opaqueness */
//   glPolygonStipple(stippleMask[16]); /* 100% opaqueness */

const GLubyte stippleMask[17][128] =
{
  /* NOTE: 0% opaqueness is faster to set and probably faster to render with:
	glDisable(GL_POLYGON_STIPPLE);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); */
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},

  {0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},

  {0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x00, 0x00, 0x00, 0x00},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x44, 0x44, 0x44, 0x44, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x11, 0x11, 0x11, 0x11},

  {0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55},

  {0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0x55, 0x55},

  {0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xee, 0xee, 0xee, 0xee, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55},

  {0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xbb, 0xbb, 0xbb, 0xbb, 0x55, 0x55, 0x55, 0x55},

  {0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55},

  {0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x55, 0x55, 0x55, 0x55},

  {0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xdd, 0xdd, 0xdd, 0xdd, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77},

  {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x77, 0x77, 0x77, 0x77},

  /* NOTE: 100% opaqueness is faster to set and probably faster to render with:
        glDisable(GL_POLYGON_STIPPLE); */
  {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
};



/////////////////////////////////////////////////////////////////////////////
// CGLDrawWnd

CToolBar CGLDrawWnd::ToolBar;

//GLsizei CGLDrawWnd::SelBufSize = 20;
const int MAX_N_TEXTURES = 1; //one is needed to pass texture //for storing textures: max. tested value: 136000; //approx. 2^17
const int MAX_N_FE_COLOR_TEXTURES = 8; 


CGLDrawWnd::CGLDrawWnd() :
m_pDC(NULL),
AxesPosition(0),
bFittingTheView(false),
bShowLegend(1),
bRotationTimerRunning(false)
{
	//bLighting = TRUE;
	prohibit_redraw = 0;
	StartingPoint.x = -1000000000;
	StartingPoint.y = -1000000000;
	Action = ActionNone;

	DetailLevel = 2.8f;
	TranslationStepCoeff = .01f;
	RotationStep = 0.5f;
	DetailLevelStepCoeff = 0.005f; // 0.00025
	PerspectiveStepCoeff = 0.01f;
	AspectRatio = 1;
	CenterPoint.x = CenterPoint.y = 0;
	CenterOffset.x = CenterOffset.y = CenterOffset.z = 0;
	SCENE_OFFSET_COEFF = 1.5f;

	saved_modelview_matrix[0] = NoSavedModelviewMatrix;

	char buf[2000];
	GetCurrentDirectory(2000,buf);
	ProgramDirectory = buf;

	BackgroundR = 0;
	BackgroundG = 0;
	BackgroundB = 0.3f;

	ScrTextR = 1.0f;
	ScrTextG = 1.0f;
	ScrTextB = 0.7f;

	//textures:
	gltexturecnt = 0; //counter for OpenGL textures
	gltexturecnt_warned = 0;
	nmaxtextures = MAX_N_TEXTURES;
	glstoretextures = 0; //do not store textures
}

void CGLDrawWnd::SetWCDI(WCDInterface * p_wcdi)
{
	pWCDI = p_wcdi;
	MaxSceneCoord = pWCDI->GetSceneMaxAbsCoordinate();
	bNoRotationState = !pWCDI->AllowRotation();
	DialogFramesRecording.SetWCDI(pWCDI);  
}

CGLDrawWnd::~CGLDrawWnd()
{
	glDeleteTextures(MAX_N_TEXTURES, texNames);
	delete[] texNames;

	glDeleteTextures(MAX_N_FE_COLOR_TEXTURES, FEcolortexNames);
	delete[] FEcolortexNames;

	/*	if(hEnhMetaFile)
	DeleteEnhMetaFile(hEnhMetaFile);*/
}


BEGIN_MESSAGE_MAP(CGLDrawWnd, CWnd)
	//{{AFX_MSG_MAP(CGLDrawWnd)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_LBUTTONDBLCLK()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_COMMAND(ID_BUTTON_NO_ROTATION, OnButtonNoRotation)
	ON_COMMAND(ID_BUTTON_STANDARD_VIEW_XY, OnButtonStandardViewXy)
	ON_COMMAND(ID_BUTTON_STANDARD_VIEW_XYZ, OnButtonStandardViewXyz)
	ON_COMMAND(ID_BUTTON_STANDARD_VIEW_XZ, OnButtonStandardViewXz)
	ON_COMMAND(ID_BUTTON_STANDARD_VIEW_YZ, OnButtonStandardViewYz)
	ON_COMMAND(ID_BUTTON_STANDARD_VIEW_ROTATE, OnButtonStandardViewRotate)
	ON_COMMAND(ID_BUTTON_FIT, OnButtonFit)
	ON_COMMAND(ID_BUTTON_MOVE_AXES, OnButtonMoveAxes)
	ON_WM_TIMER()
	ON_COMMAND(ID_BUTTON_SAVE_IMAGE, OnButtonSaveImage)
	ON_COMMAND(ID_BUTTON_FRAMES_RECORDING, OnButtonFramesRecording)
	//}}AFX_MSG_MAP
	ON_MESSAGE(WM_REDRAW,OnRedraw)
	ON_NOTIFY_EX( TTN_NEEDTEXT, 0, ToolBarToolTipsSupport )
//	ON_WM_KEYUP()
ON_WM_MOUSEWHEEL()
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CGLDrawWnd message handlers

int CGLDrawWnd::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CWnd::OnCreate(lpCreateStruct) == -1)
		return -1;

	m_pDC = new CClientDC(this);
	ASSERT(m_pDC != NULL);

	Init();

	if (!ToolBar.CreateEx(GetParent(), TBSTYLE_FLAT | TBSTYLE_TOOLTIPS, WS_CHILD | WS_VISIBLE | CBRS_TOP
		| CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC  )
		||
		!ToolBar.LoadToolBar(IDR_TOOLBAR_GL))
	{
		ASSERT(FALSE);
	}
	ToolBar.SetOwner(this);
	EnableToolTips(TRUE);


	return 0;
}

BOOL CGLDrawWnd::ToolBarToolTipsSupport( UINT id, NMHDR * pTTTStruct, LRESULT * pResult )
{
	TOOLTIPTEXT *pTTT = (TOOLTIPTEXT *)pTTTStruct;
	UINT nID =pTTTStruct->idFrom;
	if(nID)
	{
		pTTT->lpszText = MAKEINTRESOURCE(nID);
		pTTT->hinst = AfxGetResourceHandle();
		return(TRUE);
	}
	return(FALSE);

}

void CGLDrawWnd::Init()
{
	PIXELFORMATDESCRIPTOR pfd;
	int         n;

	if (!SetupPixelFormat())
		return;

	n = ::GetPixelFormat(m_pDC->GetSafeHdc());
	::DescribePixelFormat(m_pDC->GetSafeHdc(), n, sizeof(pfd), &pfd);

	hrc = wglCreateContext(m_pDC->GetSafeHdc());
	wglMakeCurrent(m_pDC->GetSafeHdc(), hrc);

	CRect r;
	GetWindowRect(&r);
	glViewport(0, 0, r.Width(), r.Height());
	AspectRatio = (GLdouble)r.Width()/r.Height();

	//setup open GL parameters
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);

	SelectObject(m_pDC->GetSafeHdc(), GetStockObject (SYSTEM_FONT)); 
	wglUseFontBitmaps (m_pDC->GetSafeHdc(), 0, 255, 0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	if(saved_modelview_matrix[0] == NoSavedModelviewMatrix)
		OnButtonStandardViewXyz();
	else
	{
		glLoadIdentity();
		glMultMatrixf(saved_modelview_matrix);
	}

	gltexturecnt = 0;
	texNames = new unsigned int[MAX_N_TEXTURES];
	glGenTextures(MAX_N_TEXTURES, texNames);

	FEcolortexNames = new unsigned int[MAX_N_FE_COLOR_TEXTURES];
	glGenTextures(MAX_N_FE_COLOR_TEXTURES, FEcolortexNames);

	ResetOpenGLParam();

	SetPerspective();

	if(bNoRotationState)
		OnButtonStandardViewXy();

	GL_SCENE_LIST = glGenLists(1);

	if(pWCDI->GetIOption(221))	//$!DR 2013-12-18 do not compute or draw anything if no hotint window is shown
	{
		this->MaxSceneCoord = pWCDI->GetSceneMaxAbsCoordinate(); // (AD) recompute bounding box on init 
		Redraw();

		if(saved_modelview_matrix[0] == NoSavedModelviewMatrix)
			OnButtonFit();
	}

 }

void CGLDrawWnd::ContentsChanged(int forcefit) //Redraw and recompute scene size
{
	if(pWCDI->GetIOption(221))	//$!DR 2013-12-18 do not compute or draw anything if no hotint window is shown
	{
		float newmax = pWCDI->GetSceneMaxAbsCoordinate();

		if (newmax > 1.1*MaxSceneCoord || newmax < 0.5*MaxSceneCoord || forcefit)
		{
			char str[100];
			//sprintf_s(str, "old max=%g, new max=%g", MaxSceneCoord, newmax);
			MaxSceneCoord = newmax;

			ButtonFit();
			//Redraw();
		}
		else
		{
			Redraw();
		}
	}
}

void CGLDrawWnd::ResetOpenGLParam()
{
	//+++++++++++++++++

	glEnable(GL_LINE_SMOOTH);

//OpenGL 2.0 or higher: ??? does not work!, need to tile triangles!
//#define GL_PHONG_WIN                      0x80EA
//#define GL_PHONG_HINT_WIN                 0x80EB
//#define GL_SHADOW_AMBIENT_SGIX            0x80BF

	if (pWCDI->GetIOption(207))
	{
		glShadeModel(GL_SMOOTH);
	}
	else
		glShadeModel(GL_FLAT);
 
	glDisable (GL_LIGHTING);
 
	glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//++++++++++++++++++++++++++



	glLineWidth(1.0f); 
	glPointSize(4.0f);


  GLfloat lightamb   = (float)pWCDI->GetDOption(207);
  GLfloat lightdiff  = (float)pWCDI->GetDOption(208);
  GLfloat lightspec  = (float)pWCDI->GetDOption(209);
  GLfloat lightamb2  = (float)pWCDI->GetDOption(210);
  GLfloat lightdiff2 = (float)pWCDI->GetDOption(211);
  GLfloat lightspec2 = (float)pWCDI->GetDOption(212);
	GLfloat mat_shininess = (float)pWCDI->GetDOption(219);

  GLfloat mat_spec_col[4];					// = { 1, 1, 1, 1 };
	mat_spec_col[0] = (float)pWCDI->GetDOption(220);
	mat_spec_col[1] = (float)pWCDI->GetDOption(220);
	mat_spec_col[2] = (float)pWCDI->GetDOption(220);
	mat_spec_col[3] = 1.f;

  GLfloat light_position1[4];				// = { 1, 1, -3, 0 }; //position x,y,z and 0 for directional source or 1 for exact light source
	light_position1[0] = (float)pWCDI->GetDOption(213);
	light_position1[1] = (float)pWCDI->GetDOption(214);
	light_position1[2] = (float)pWCDI->GetDOption(215);
	light_position1[3] = (float)pWCDI->GetIOption(210);

	GLfloat light_position2[4];				// = { 0, 3, 2, 0 };
	light_position2[0] = (float)pWCDI->GetDOption(216);
	light_position2[1] = (float)pWCDI->GetDOption(217);
	light_position2[2] = (float)pWCDI->GetDOption(218);
	light_position2[3] = (float)pWCDI->GetIOption(211);

  GLfloat vals[3];
  vals[0] = vals[1] = vals[2] = lightamb;
  glLightfv(GL_LIGHT0, GL_AMBIENT, vals);
  vals[0] = vals[1] = vals[2] = lightdiff;
  glLightfv(GL_LIGHT0, GL_DIFFUSE, vals);
  vals[0] = vals[1] = vals[2] = lightspec;
  glLightfv(GL_LIGHT0, GL_SPECULAR, vals);

  glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
  
  vals[0] = vals[1] = vals[2] = lightamb2;
  glLightfv(GL_LIGHT1, GL_AMBIENT, vals);
  vals[0] = vals[1] = vals[2] = lightdiff2;
  glLightfv(GL_LIGHT1, GL_DIFFUSE, vals);
  vals[0] = vals[1] = vals[2] = lightspec2;
  glLightfv(GL_LIGHT1, GL_SPECULAR, vals);

  glLightfv(GL_LIGHT1, GL_POSITION, light_position2);


	glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);
	//glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
  //glMaterialfv (GL_FRONT, GL_SPECULAR, mat_spec_col);

	glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 0);
	GLint locviewer = 0;
  glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, locviewer);

  if (pWCDI->GetIOption(208)) glEnable (GL_LIGHT0);
	else glDisable (GL_LIGHT0);
  if (pWCDI->GetIOption(209)) glEnable (GL_LIGHT1);
  else glDisable (GL_LIGHT1);

	
	switch (pWCDI->GetIOption(213))
	{
	case 1:
		glCullFace(GL_FRONT);
		glEnable(GL_CULL_FACE);
		break;
	case 2:
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
		break;
	case 3:
		glCullFace(GL_FRONT_AND_BACK);
		glEnable(GL_CULL_FACE);
		break;
	default:
		glDisable(GL_CULL_FACE);
	}

	CreateFEColorTexture(1024, 0, 0); //color (iso/no iso)
	//CreateFEColorTexture(256, 1, 1); //grey

}

void CGLDrawWnd::CreateFEColorTexture (int ncols, int texnum, int grey)
{
	int i;

	if (ncols < 2) ncols = 2;

	colortexture = new GLubyte[4*ncols+8];

	const double colp[][3] = //color
	{
		//{ 1, 0, 0 }, //NETGEN
		//{ 1, 1, 0 },
		//{ 0, 1, 0 },
		//{ 0, 1, 1 },
		//{ 0, 0, 1 },
		//{ 1, 0, 1 },
		//{ 1, 0, 0 },
		//{ 0.25, 0.25, 0.25 },
		{ 0.1, 0.1, 0.1 },
		{ 0.1, 0.1, 0.9 },
		{ 0.1, 0.9, 0.9 },
		{ 0.1, 0.9, 0.1 },
		{ 0.9, 0.9, 0.1 },
		{ 0.9, 0.1, 0.1 },
		{ 0.9, 0.9, 0.9 },
		{ 0.9, 0.9, 0.9 }
	};

	const double colpg[][3] = //grey
	{
		{ 0.1, 0.1, 0.1 },
		{ 0.2, 0.2, 0.2 },
		{ 0.8, 0.8, 0.8 },
		{ 0.8, 0.8, 0.8 },
		{ 0.8, 0.8, 0.95}
	};

	double niso = pWCDI->GetIOption(103);

	for (i = 0; i < ncols; i++)
	{
		double value = ((double)i)/ (ncols-1.);


		if (niso != 33 && niso != 0)
		{
			value = ((int)((niso)*value+0*0.5)); //for isolines ...
			int niso2 = (int)niso-1; if (niso2 == 0) niso2 = 1;
			value = value/niso2+1e-12;
		}

		if (pWCDI->GetIOption(105))
			value *= 1.;
		else
			value *= 4.;

		int iv = int(value);
		double r = value - iv;

		iv += 1;
		if (i == 0) { iv = 0; r=0;}
		if (i == ncols-1) {iv = 5; r=1;}


		GLdouble col[3];
		int j;
		for (j = 0; j < 3; j++)
		{
			if (!pWCDI->GetIOption(105))
				col[j] = (1-r) * colp[iv][j] + r * colp[iv+1][j];
			else
				col[j] = (1-r) * colpg[iv][j] + r * colpg[iv+1][j];
		}

		//glColor3d(col[0], col[1], col[2]); //??????
		//char str[256];
		//sprintf_s(str, "R=%.3f, G=%.3f, B=%.3f\n", col[0], col[1],col[2]);
		//pWCDI->GetUserInterface()->AddText(str);

		colortexture[4*i] = GLubyte (255 * col[0]);
		colortexture[4*i+1] = GLubyte (255 * col[1]);
		colortexture[4*i+2] = GLubyte (255 * col[2]);
		colortexture[4*i+3] = GLubyte(255);
	}


	glBindTexture (GL_TEXTURE_1D, FEcolortexNames[texnum]);

	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexImage1D (GL_TEXTURE_1D, 0, GL_RGBA, ncols, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);
	glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

	delete [] colortexture;

	
	//cout << "linear = " << linear << endl;
	//glBindTexture (GL_TEXTURE_1D, coltexname);
	//if (linear)
	//{
	//glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//}
	//else
	//{
	//glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//}
}



BOOL CGLDrawWnd::SetupPixelFormat()
{
	static PIXELFORMATDESCRIPTOR pfd = 
	{
		sizeof(PIXELFORMATDESCRIPTOR),  // size of this pfd
			1,                              // version number
			PFD_DRAW_TO_WINDOW |            // support window
			PFD_SUPPORT_OPENGL |          // support OpenGL
			PFD_DOUBLEBUFFER,             // double buffered
			PFD_TYPE_RGBA,                  // RGBA type
			24,                             // 24-bit color depth
			0, 0, 0, 0, 0, 0,               // color bits ignored
			0,                              // no alpha buffer
			0,                              // shift bit ignored
			0,                              // no accumulation buffer
			0, 0, 0, 0,                     // accum bits ignored
			32,                             // 32-bit z-buffer
			0,                              // no stencil buffer
			0,                              // no auxiliary buffer
			PFD_MAIN_PLANE,                 // main layer
			0,                              // reserved
			0, 0, 0                         // layer masks ignored
	};
	int pixelformat;

	if ( (pixelformat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd)) == 0 )
	{
		MessageBox("ChoosePixelFormat failed");
		return FALSE;
	}

	if (SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, &pfd) == FALSE)
	{
		MessageBox("SetPixelFormat failed");
		return FALSE;
	}

	return TRUE;
}

void CGLDrawWnd::StopOpenGL()
{
	wglMakeCurrent(NULL,NULL);
	wglDeleteContext(hrc);
}

// the structured text vertical position on the screen
float CGLDrawWnd::GetStructTextVertPos(int nLineNo)
{
	float height = (float)(m_pDC->GetTextExtent("a").cy);
	CRect r;
	GetWindowRect(&r);
	return 1-2*height*(nLineNo+1)/r.Height();
}

// the structured text horizontal position on the screen
float CGLDrawWnd::GetStructTextHorPos(int nXPos,const CString & text)
{
	if(nXPos < 0)
		return -1.0f;		// left side
	float width = (float)(m_pDC->GetTextExtent(text).cx);
	CRect r;
	GetWindowRect(&r);
	if(nXPos > 0)
		return 1.0f - 2*width/r.Width();	// right side
	return -width/r.Width();	// right side
}

// here the texts with fixed positions are displayed inthe OpenGL scene
void CGLDrawWnd::PrintTextsFixed()
{
	if(bFittingTheView)
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	::glColor3f(ScrTextR, ScrTextG, ScrTextB);
	glListBase(0);
	//glTranslated(0.f,0.f,4.f);

	POSITION pos = Texts.GetHeadPosition();
	while(pos)
	{
		const TextPortion & tp = Texts.GetNext(pos);
		switch(tp.type)
		{
		case TextPortion::Text3D:
			continue;
		case TextPortion::Text2D:
			glRasterPos3f(tp.x, tp.y, -1.0f);
			break;
		case TextPortion::TextStruct:
			glRasterPos3f(GetStructTextHorPos(tp.nXPos,tp.text),GetStructTextVertPos(tp.nLineNo),-1.0f);
			break;
		default: ASSERT(FALSE);
		}
		glCallLists(tp.text.GetLength(), GL_UNSIGNED_BYTE, tp.text);
	}
	
	if(bNoRotationState)
	{
		static const CString NoRotation = "No rotation";
		glRasterPos3f(-1.f, -1.f, -1.f);
		glCallLists(NoRotation.GetLength(), GL_UNSIGNED_BYTE, NoRotation);
	}
	if(bRotationTimerRunning)
	{
		static const CString AutoRotation = "Automatic rotation";
		glRasterPos3f(GetStructTextHorPos(0,AutoRotation), -1.0f, 0.0f);
		glCallLists(AutoRotation.GetLength(), GL_UNSIGNED_BYTE, AutoRotation);
	}
	if(DialogFramesRecording.m_bCheckRecordFrames && DialogFramesRecording.m_bShowFrameNumbers)
	{
		CString FramesRecording;
		FramesRecording.Format("frame %d",DialogFramesRecording.m_nFrameCounter);
		glRasterPos3f(GetStructTextHorPos(1,FramesRecording), -1.0f, 0.0f);
		glCallLists(FramesRecording.GetLength(), GL_UNSIGNED_BYTE, FramesRecording);
	}

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void CGLDrawWnd::OnPaint() 
{
	SetPerspective();

	CPaintDC dc(this); // device context for painting

	glClearColor(BackgroundR, BackgroundG, BackgroundB, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// the text should not be affected by the lighting effects
	glDisable(GL_LIGHTING);
	PrintTextsFixed();

	DrawScene();

	glFinish();
	SwapBuffers(wglGetCurrentDC());
}

// the main painting function
void CGLDrawWnd::DrawScene()
{
	float m[16];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	ResetOpenGLParam();
	glPopMatrix();


	// first we draw the axes
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPerspective((GLdouble)10,AspectRatio,6,10);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	glLoadIdentity();
	float f1 = (1.27f*(float)AspectRatio-0.25f);
	float f2 = 0.85f;
	glTranslated(-0.5f*f1*f2+0.01f,-0.5f*f2-0.02f,-8.f+1.5f);
	switch(AxesPosition)
	{
	case 1: glTranslated(1.f*f1*f2,0,0); break;
	case 2: glTranslated(1.f*f1*f2,1.0f*f2,0); break;
	case 3: glTranslated(0,1.0f*f2,0); break;
	case 4: glTranslated(0.5f*f1*f2,0.5f*f2,0); break;
	} 
	glMultMatrixf(m);
	ChooseAxesColor();
	if(AxesPosition != 5 && !bFittingTheView)
	{
		glBegin(GL_LINES);
		glVertex3f(0,0,0);
		glVertex3f(0.1f,0,0);
		glVertex3f(0,0,0);
		glVertex3f(0,0.1f,0);
		glVertex3f(0,0,0);
		glVertex3f(0,0,-0.1f);
		::glEnd();
		// names of the axes
		::glColor3f(ScrTextR, ScrTextG, ScrTextB);
		glListBase(0);
		glRasterPos3f(0.13f, 0.0f, 0.0f);
		glCallLists(1, GL_UNSIGNED_BYTE, "x");
		glRasterPos3f(0.0f, 0.13f, 0.0f);
		glCallLists(1, GL_UNSIGNED_BYTE, "y");
		glRasterPos3f(0.0f, 0.0f, -0.13f);
		glCallLists(1, GL_UNSIGNED_BYTE, "z");
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);


	// getting prepared to paint the rest
	glMatrixMode(GL_MODELVIEW);
	

	glPushMatrix();
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	glLoadIdentity();
	glTranslated(CenterPoint.x,CenterPoint.y,(-MaxSceneCoord)*(1+SCENE_OFFSET_COEFF));
	//glTranslated(CenterPoint.x-CenterOffset.x,CenterPoint.y+CenterOffset.y,(-MaxSceneCoord+CenterOffset.z)*(1+SCENE_OFFSET_COEFF));

	glMultMatrixf(m);
	glTranslated(-CenterOffset.x,-CenterOffset.y,CenterOffset.z);

	if(pWCDI->GetIOption(206))
		glEnable(GL_LIGHTING);

	glCallList(GL_SCENE_LIST);
	glDisable(GL_LIGHTING);

	if(pWCDI->GetIOption(157))
	{
		if(pWCDI->GetIOption(149)) //activate clipping plane 1
		{
			GLdouble eqn[4] = {
				pWCDI->GetDOption(125), // x-component of normal
				pWCDI->GetDOption(126), // y-component of normal
				pWCDI->GetDOption(127), // z-component of normal
				pWCDI->GetDOption(128)  // distance from origin
			};
			glClipPlane(GL_CLIP_PLANE0, eqn);
			glEnable(GL_CLIP_PLANE0);
		}
		else
		{
			glDisable(GL_CLIP_PLANE0);
		}

		if(pWCDI->GetIOption(156)) //activate clipping plane 2
		{
			GLdouble eqn[4] = {
				pWCDI->GetDOption(129), // x-component of normal
				pWCDI->GetDOption(130), // y-component of normal
				pWCDI->GetDOption(131), // z-component of normal
				pWCDI->GetDOption(132)  // distance from origin
			};
			glClipPlane(GL_CLIP_PLANE1, eqn);
			glEnable(GL_CLIP_PLANE1);
		}
		else
		{
			glDisable(GL_CLIP_PLANE1);
		}
	}

	glPopMatrix();

	// and now the text3D, draw in front of text:
	glPushMatrix();
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	glLoadIdentity();
	glTranslated(CenterPoint.x,CenterPoint.y,(-MaxSceneCoord)*(1+SCENE_OFFSET_COEFF));

	if (pWCDI->GetIOption(122)) glTranslated(0.f,0.f,-1.f); //option: bring text to front

	glMultMatrixf(m);
	glTranslated(-CenterOffset.x,-CenterOffset.y,CenterOffset.z);

	if(!bFittingTheView)
	{
		::glColor3f(ScrTextR, ScrTextG, ScrTextB);
		glListBase(0);
		POSITION pos = Texts.GetHeadPosition();
		while(pos)
		{
			const TextPortion & tp = Texts.GetNext(pos);
			if(tp.type == TextPortion::Text3D)
			{
				glRasterPos3f(tp.x, tp.y, tp.z);
				glCallLists(tp.text.GetLength(), GL_UNSIGNED_BYTE, tp.text);
			}
		}
	}

	glPopMatrix();

}

// creating the perspective matrix
void CGLDrawWnd::SetPerspective()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//$!PG 2011-3-15:[
	//$ PG 2011-3-15: commented the following line, in order to avoid the scene to be rescaled during load/time iterations
	//MaxSceneCoord = pWCDI->GetSceneMaxAbsCoordinate();
	//$!PG 2011-3-15:]
	gluPerspective((GLdouble)(MAX_FOVY/DetailLevel), AspectRatio,
		MaxSceneCoord*SCENE_OFFSET_COEFF, MaxSceneCoord*(2+SCENE_OFFSET_COEFF));
	//orig:
	//gluPerspective((GLdouble)(MAX_FOVY/DetailLevel), AspectRatio,
	//	MaxSceneCoord*SCENE_OFFSET_COEFF, MaxSceneCoord*(2+SCENE_OFFSET_COEFF));

	glMatrixMode(GL_MODELVIEW);


}

LRESULT CGLDrawWnd::OnRedraw(WPARAM, LPARAM)
{
	Redraw();

	// Eat spurious WM_REDRAW messages
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_REDRAW, WM_REDRAW, PM_REMOVE));

	return 0;
}

extern int vertexcnt;
extern int gltrigcnt;
extern int glquadcnt;
extern int gllinecnt;

// rebuild the scene
void CGLDrawWnd::Redraw()
{
	if (prohibit_redraw) return;

	Texts.RemoveAll();

	//gluLookAt(0,0,10,0,0,5,0,1,0);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//texture:
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	gltexturecnt = 0;
	texturemode = 1-glstoretextures;
	glNewList(GL_SCENE_LIST,GL_COMPILE); //COMPILE, only for textures

	pWCDI->RenderScene(this);

	glEndList();

/*
	if (gltexturecnt > 100 || !glstoretextures)
	{
		//COMPILE, only for textures
		//do not include textures in call-list ...
		glCallList(GL_SCENE_LIST);
		gltexturecnt = 0;
		texturemode = 0;
		glNewList(GL_SCENE_LIST,GL_COMPILE); //COMPILE remaining objects
		pWCDI->RenderScene(this);
		glEndList();
	}
*/

	/*
	char str[100];
	sprintf_s(str, "vertexcnt=%d", vertexcnt);
	PrintText2D(-0.99,-0.05,str);
	sprintf_s(str, "linecnt=%d", gllinecnt);
	PrintText2D(-0.99,-0.10,str);
	sprintf_s(str, "trigcnt=%d", gltrigcnt);
	PrintText2D(-0.99,-0.15,str);
	sprintf_s(str, "quadcnt=%d", glquadcnt);
	PrintText2D(-0.99,-0.20,str);
	int textcnt = Texts.GetCount();
	sprintf_s(str, "textcnt=%d", textcnt);
	PrintText2D(-0.99,-0.25,str);
	*/

	vertexcnt = 0;
	gltrigcnt = 0;
	glquadcnt = 0;
	gllinecnt = 0;
	
	SetPerspective();

	RedrawWindow();

	PerformVideoFramesRecord();	
}

// the size of the window has changed
void CGLDrawWnd::OnSize(UINT nType, int cx, int cy) 
{
	CWnd::OnSize(nType, cx, cy);

	if(cy > 0) 
	{    
		glViewport(0, 0, cx, cy);
		AspectRatio = (GLdouble)cx/cy;
		Redraw();
	}

	CRect r;
	GetParent()->GetWindowRect(&r);

	//ToolBar.MoveWindow(r.Width()-ToolBarButtons*ToolBarButtonSize-5-2,-2,r.Width(),ToolBarButtonSize-2); //!AD: 2012-07-25 in Win7 the last button is clipped by window frame
	ToolBar.MoveWindow(r.Width()-ToolBarButtons*ToolBarButtonSize-5-12,-2,r.Width(),ToolBarButtonSize-2); //!AD: 2012-07-25 looks nicer for Win7 
}

BOOL CGLDrawWnd::OnEraseBkgnd(CDC* pDC) 
{
	return 1;		//Prevents flickering

	//return CView::OnEraseBkgnd(pDC);
}

void CGLDrawWnd::OnDestroy() 
{
	StopOpenGL();
	CWnd::OnDestroy();	
	ToolBar.ShowWindow(SW_HIDE);
	if(m_pDC)
		delete m_pDC;
}

//1.5 is standard value, smaller for higher distorted perspective!
void CGLDrawWnd::SetSceneOffsetCoeff(float x) {SCENE_OFFSET_COEFF = x;}

void CGLDrawWnd::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default

	//$ AD 2013-5-10: [ new code to catch an Element object in the GL Window

	// define structure of name buffer
	GLuint buffer[32];
	GLubyte* p;

	
	// switch to selection mode
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);


	//$ AD 2013-5-10: end ]

	RedrawWindow();


	CWnd::OnLButtonDblClk(nFlags, point);
}

// begin action - BoxZoom
void CGLDrawWnd::OnLButtonDown(UINT nFlags, CPoint point) 
{
	Action = ActionMove; //if no button is pressed
	if(!(nFlags & MK_CONTROL) && (nFlags & MK_SHIFT))
		Action = ActionSelect;
	
	CurrentPoint = StartingPoint = point;
	SetCapture();


	CWnd::OnLButtonDown(nFlags, point);
}

// end action -- always zoom
void CGLDrawWnd::OnLButtonUp(UINT nFlags, CPoint point) 
{
	if (!prohibit_redraw && Action == ActionSelect/*(StartingPoint.x != -1000000000)*/) //catch, if file dialog in OnFileOpenMBS is double-clicked!!!
	{
		ReleaseCapture();

		MakeZoom(StartingPoint,point);
		//	else
		//SelectRegion(StartingPoint,point,(nFlags & MK_SHIFT) != 0,(nFlags & MK_CONTROL) != 0);
		RedrawWindow();
	}
	else
	{
		ReleaseCapture();
	}
	Action = ActionNone;

	CWnd::OnLButtonUp(nFlags, point);
}

// the current action must be processed
void CGLDrawWnd::OnMouseMove(UINT nFlags, CPoint point) 
{
	if(	(Action == ActionRotate) ||
		(Action == ActionZoom) ||
		(Action == ActionPerspective) ||
		(Action == ActionMove) )
	{
		float m[16];
		glGetFloatv(GL_MODELVIEW_MATRIX,m);
		glLoadIdentity();
		switch(Action)
		{
		case ActionMove:
			{
				CRect r;
				GetWindowRect(&r);
				float c2 = 290.f/(float)r.Height(); //depends on zoom and perspective!!!
				float coeff = TranslationStepCoeff*MaxSceneCoord * (SCENE_OFFSET_COEFF+1.f) *c2/DetailLevel;
				CenterPoint.x += (point.x-CurrentPoint.x)*coeff;
				CenterPoint.y += -(point.y-CurrentPoint.y)*coeff;
			}
			break;
		case ActionZoom:
			{
			//DetailLevel *= 1 - (point.y-CurrentPoint.y)*DetailLevelStepCoeff/* *MaxSceneCoord*/;

			float fact = (float)(point.y-CurrentPoint.y);
			if (fact < 0)
				DetailLevel *= 1.f + fabsf(fact)*DetailLevelStepCoeff/* *MaxSceneCoord*/;
			else
				DetailLevel *= 1.f/(1.f + fabsf(fact)*DetailLevelStepCoeff)/* *MaxSceneCoord*/;

			/*
			DetailLevel -= 1;
			DetailLevel *= 1 - (point.y-CurrentPoint.y)*DetailLevelStepCoeff*MaxSceneCoord;
			DetailLevel += 1;*/
			DetailLevel = max(DetailLevel,0.9f);
			SetPerspective();
			//Redraw();		// can be used to remove defects while zooming
			//  if BeginLinesOnSolid() are used in the scene,
			//  but leads to flickering
			break;
			}
		case ActionPerspective:
			SCENE_OFFSET_COEFF *= 1 - (CurrentPoint.y - point.y)*PerspectiveStepCoeff/* *MaxSceneCoord*/;
			SCENE_OFFSET_COEFF = max(SCENE_OFFSET_COEFF,0.001f);
			SCENE_OFFSET_COEFF = min(SCENE_OFFSET_COEFF,100.f);
			SetPerspective();
			break;
		case ActionRotate:
			glRotatef	(	(point.x-CurrentPoint.x)*RotationStep,
				0.0f,
				1.0f,
				0.0f					);
			glRotatef	(	(point.y-CurrentPoint.y)*RotationStep,
				1.0f,
				0.0f,
				0.0f					);
			break;					
		}
		glMultMatrixf(m);
		RedrawWindow();
		PerformVideoFramesRecord();
	}
	if(Action == ActionSelect)
	{
		CRect r1(StartingPoint,point);
		CRect r2(StartingPoint,CurrentPoint);
		r1.NormalizeRect();
		r2.NormalizeRect();
		m_pDC->DrawDragRect(r1,CSize(1,1),r2,CSize(1,1));
	}
	CurrentPoint = point;


	CWnd::OnMouseMove(nFlags, point);
}

void CGLDrawWnd::OnMouseWheelGLDW(UINT nFlags, short zDelta, CPoint pt)
{
	//zDelta: Indicates distance rotated. The zDelta value is expressed in multiples or divisions of WHEEL_DELTA, which is 120. 
	//A value less than zero indicates rotating back (toward the user) while a value greater than zero indicates rotating forward (away from the user). The user can reverse this response by changing the Wheel setting in the mouse software. 
	//See the Remarks for more information about this parameter. 
	
	//float m[16];
	//glGetFloatv(GL_MODELVIEW_MATRIX,m);
	//glLoadIdentity();

	double fact = (double)zDelta/4.; //one wheel step is equal to moving the mouse 30 pts

	if (fact > 0)
		DetailLevel *= 1. + fabs(fact)*DetailLevelStepCoeff/* *MaxSceneCoord*/;
	else
		DetailLevel *= 1./(1. + fabs(fact)*DetailLevelStepCoeff)/* *MaxSceneCoord*/;

	DetailLevel = max(DetailLevel,0.9f);
	//SetPerspective();
	Redraw();		// can be used to remove defects while zooming

	//glMultMatrixf(m);
	//RedrawWindow();
	PerformVideoFramesRecord();
}

// begin action depending on the buttons currently pressed
void CGLDrawWnd::OnRButtonDown(UINT nFlags, CPoint point) 
{
	SetCapture();
	CurrentPoint = StartingPoint = point;
	Action = ActionRotate;
	if( ((nFlags & MK_CONTROL) && !(nFlags & MK_SHIFT)) || bNoRotationState)
		Action = ActionMove;
	if(!(nFlags & MK_CONTROL) && (nFlags & MK_SHIFT))
		Action = ActionZoom;
	if((nFlags & MK_CONTROL) && (nFlags & MK_SHIFT))
		Action = ActionPerspective;


	CWnd::OnRButtonDown(nFlags, point);
}

// end action
void CGLDrawWnd::OnRButtonUp(UINT nFlags, CPoint point) 
{
	ReleaseCapture();
	if(Action == ActionZoom)
	{
		Redraw();
	}
	/*
	else if(Action == ActionRotate)
	{
		double dist = sqrt((point.x-StartingPoint.x)*(point.x-StartingPoint.x)+(point.y-StartingPoint.y)*(point.y-StartingPoint.y));
		if (dist < 1)
		{
			DetailLevel *= 0.5;
			Redraw();
		}
	}*/
	Action = ActionNone;


	CWnd::OnRButtonUp(nFlags, point);
}

/*void CGLDrawWnd::SelectRegion(const CPoint & p1, const CPoint & p2, bool bShift, bool bCtrl)
{
GLuint * buf;
bool bSizeEnough;

do
{
bSizeEnough = true;
buf = new GLuint[SelBufSize];
Selection.RemoveAll();

glSelectBuffer( SelBufSize-1, buf );
glRenderMode( GL_SELECT );

GLfloat mproj[16];
glMatrixMode( GL_PROJECTION );
glGetFloatv( GL_PROJECTION_MATRIX, mproj );
glPushMatrix();

GLint vp[4];
glGetIntegerv( GL_VIEWPORT, vp );
glLoadIdentity();

int x = (p1.x+p2.x)/2;
int y = (p1.y+p2.y)/2;
int width = abs(p1.x-p2.x);
int height = abs(p1.y-p2.y);
if(width < 3)
width = 3;
if(height < 3)
height = 3;
gluPickMatrix( x, vp[3]-y, width, height, vp );
glMultMatrixf( mproj );
glMatrixMode( GL_MODELVIEW );

float m[16];
glPushMatrix();
glGetFloatv(GL_MODELVIEW_MATRIX,m);
glLoadIdentity();
glTranslated(CenterPoint.x,CenterPoint.y,-m_fViewDepth);
glMultMatrixf(m);
glCallList(GL_LIST);
glFinish();
glPopMatrix();

glMatrixMode( GL_PROJECTION );
glPopMatrix();
glMatrixMode( GL_MODELVIEW );
GLint nRecords = glRenderMode( GL_RENDER );

if(nRecords == -1)		// не хватило буфера
{
bSizeEnough = false;
delete[] buf;
SelBufSize = (3*SelBufSize)/2;
}
else
{
int i=0, pos = 0;
for( i=0; i < nRecords; i++ )
{
for(unsigned int k = pos+3; k < pos+3+buf[pos]; k++)
Selection.AddTail(buf[k]);
pos += 3+buf[pos];
}
}
} while(!bSizeEnough);

pGLDFB->Select(Selection,bShift,bCtrl);

delete[] buf;
}*/


// handling the toolbar buttons

void CGLDrawWnd::OnButtonNoRotation() 
{
	bNoRotationState = !bNoRotationState;
	RedrawWindow();
}

void CGLDrawWnd::OnButtonStandardViewXy() 
{
	OnButtonStandardViewXz();
	glRotated(-90,1,0,0);
	RedrawWindow();
}

void CGLDrawWnd::OnButtonStandardViewXyz() 
{
	OnButtonStandardViewYz();
	glRotated(-90,0,0,1);
	glRotated(-90,1,0,0);
	glRotated(180,0,1,0);

	double rot[4];
	rot[1] = 0;
	rot[2] = 0;
	rot[3] = 0;
	rot[pWCDI->GetIOption(200)] = 1;
	glRotated(pWCDI->GetDOption(200),rot[1],rot[2],rot[3]);
	rot[pWCDI->GetIOption(200)] = 0;
	rot[pWCDI->GetIOption(201)] = 1;
	glRotated(pWCDI->GetDOption(201),rot[1],rot[2],rot[3]);
	rot[pWCDI->GetIOption(201)] = 0;
	rot[pWCDI->GetIOption(202)] = 1;
	glRotated(pWCDI->GetDOption(202),rot[1],rot[2],rot[3]);
	rot[pWCDI->GetIOption(202)] = 0;

	RedrawWindow();
}

void CGLDrawWnd::OnButtonStandardViewXz() 
{
	OnButtonStandardViewYz();
	glRotated(90,0,0,1);
	RedrawWindow();
}

void CGLDrawWnd::OnButtonStandardViewYz() 
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	float m[16];
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	m[0] = -1;
	glMultMatrixf(m);
	glRotated(90,1,0,0);
	glRotated(90,0,0,1);
	RedrawWindow();
}

void CGLDrawWnd::SetAnimateScalingTimer(int flag)
{
	old_scaling_factor = pWCDI->GetDOption(105);
	if(flag)
	{
		old_scaling_factor = pWCDI->GetDOption(105);
		virtual_animate_time = 0;
		SetTimer(SCALING_ANIMATE_TIMER_ID,ROTATION_TIMER_DELAY,NULL);
	}
	else
	{
		KillTimer(SCALING_ANIMATE_TIMER_ID);
		RedrawWindow();
		pWCDI->GetDOption(105) = old_scaling_factor;
	}
}


void CGLDrawWnd::OnButtonStandardViewRotate() 
{
	bNoRotationState = false;
	if(!bRotationTimerRunning)
	{
		bRotationTimerRunning = true;
		SetTimer(ROTATION_TIMER_ID,ROTATION_TIMER_DELAY,NULL);
	}
	else
	{
		bRotationTimerRunning = false;
		KillTimer(ROTATION_TIMER_ID);
		RedrawWindow();
	}
}

void CGLDrawWnd::MakeZoom(const CPoint & p1, const CPoint & p2)
{
	int MidX = (p1.x+p2.x)/2;
	int MidY = (p1.y+p2.y)/2;
	CRect r;
	GetWindowRect(r);
	CPoint CentPoint = r.CenterPoint();
	ScreenToClient(&CentPoint);
	int dx = CentPoint.x - MidX;
	int dy = CentPoint.y - MidY;

	//char str[256];
	//sprintf_s(str, "centerx=%d, centery=%d, px=%d, py=%d\n", CentPoint.x, CentPoint.y, MidX, MidY);
	//pWCDI->GetUserInterface()->AddText(str);


	//old: float coeff = MaxSceneCoord/DetailLevel;
	float c2 = 0.01f*290.f/(float)r.Height(); //depends on zoom and perspective!!!
	float coeff = MaxSceneCoord * (SCENE_OFFSET_COEFF+1.f) *c2/DetailLevel;

	// the following coefficients adjust the accuracy of zooming
	CenterPoint.x += dx*coeff;		//old: /r.Width()*7.8f
	CenterPoint.y += -dy*coeff;	//old: /r.Height()*6.5f
	int width = abs(p1.x-p2.x);
	int height = abs(p1.y-p2.y);
	float Magnification;
	if(width < 5 || height < 5)
		Magnification = 2;
	else
		Magnification = min(((float)r.Width())/width,((float)r.Height())/height);
	if (Magnification > 100) Magnification = 100;
	DetailLevel *= Magnification;

	Redraw();
}


void CGLDrawWnd::OnButtonFit() 
{
	MaxSceneCoord = pWCDI->GetSceneMaxAbsCoordinate();

	ButtonFit();
}

void CGLDrawWnd::ButtonFit() 
{
	//pWCDI->GetUserInterface()->AddText("button fit\n");

	// the scene needs to fit the screen
	CenterPoint.x = CenterPoint.y = 0;
	DetailLevel = 20.*1.5f;//1.5f;
	SCENE_OFFSET_COEFF = 10.f; //1.5

	SetPerspective();
	//RedrawWindow();

	// now we paint the scene virtually and find the bounds
	// First we fill the FeedBackBuffer in
	GLfloat *buf;    
	int size = 3000;
	GLint nValues;
	bFittingTheView = true;
	do
	{
		buf = new GLfloat[ size ];
		glFeedbackBuffer( size, GL_3D_COLOR, buf );
		glRenderMode( GL_FEEDBACK );
		DrawScene();
		nValues = glRenderMode( GL_RENDER );
		if( nValues < 0 )
		{
			delete [] buf;
			size *= 2;
		}
	}   while( size <= MAX_SIZE   &&   nValues < 0 );
	bFittingTheView = false;
	if( nValues < 0 )
	{
		AfxMessageBox( "The scene is too complicated for screen fit!", MB_ICONSTOP );
		return;
	}	
	// And now we are finding the borders of the scene
	// For simplicity we take into account only the points
	int i=0;
	int dim = 3 + 4;
	bool bPointConsidered = false;
	CRect r;
	GetWindowRect(r);
	int WindowHeight = r.Height();
	while( i < nValues )   switch( ( GLuint )buf[i] )
	{
				case GL_POINT_TOKEN:
					{
						i++;
						int x = (int)buf[i], y = WindowHeight - (int)buf[i+1];
						if(bPointConsidered)
						{
							r.left = min(r.left,x);
							r.right = max(r.right,x);
							r.top = max(r.top,y);
							r.bottom = min(r.bottom,y);
						}
						else
						{
							r.SetRect(x-2,y+2,x+2,y-2);
							bPointConsidered = true;
						}
						i += dim;
						break;
					}
				case GL_LINE_TOKEN:
				case GL_LINE_RESET_TOKEN:
					{
						i++;
						int x = (int)buf[i], y = WindowHeight - (int)buf[i+1];
						if(bPointConsidered)
						{
							r.left = min(r.left,x);
							r.right = max(r.right,x);
							r.top = max(r.top,y);
							r.bottom = min(r.bottom,y);
						}
						else
						{
							r.SetRect(x-2,y+2,x+2,y-2);
							bPointConsidered = true;
						}
						i += 2*dim;
						break;
					}
				case GL_BITMAP_TOKEN: break;
					ASSERT( FALSE );
				case GL_DRAW_PIXEL_TOKEN:
					ASSERT( FALSE );
				case GL_COPY_PIXEL_TOKEN:
					ASSERT( FALSE );
				case GL_PASS_THROUGH_TOKEN:
					ASSERT( FALSE );
				case GL_POLYGON_TOKEN:
					{
						int nVertices = ( int )buf[ ++i ];
						i++;
						int x = (int)buf[i], y = WindowHeight - (int)buf[i+1];
						if(bPointConsidered)
						{
							r.left = min(r.left,x);
							r.right = max(r.right,x);
							r.top = max(r.top,y);
							r.bottom = min(r.bottom,y);
						}
						else
						{
							r.SetRect(x-2,y+2,x+2,y-2);
							bPointConsidered = true;
						}
						i += nVertices * dim;
						break;
					}
				default:
					ASSERT( FALSE );
	}
	delete [] buf;

	r.NormalizeRect();
	r.InflateRect(10,10,10,10); //originally: r.InflateRect(50,50,50,50);
	// And finally we make a BoxZoom for the determined screen coordinates
	
	//char str[256];
	//sprintf_s(str, "MakeZoom: left=%d, top=%d, right=%d, bottom=%d\n", r.left, r.top, r.right, r.bottom);
	//pWCDI->GetUserInterface()->AddText(str); MaxSceneCoord = 1;

	MakeZoom(r.TopLeft(),r.BottomRight());
	// DEBUG
	/*CClientDC dc(this);
	dc.LineTo(r.left,r.top);
	dc.LineTo(r.left,r.bottom);
	dc.LineTo(r.right,r.bottom);
	dc.LineTo(r.right,r.top);
	dc.LineTo(r.left,r.top);*/
}

void CGLDrawWnd::OnButtonMoveAxes() 
{
	AxesPosition = (AxesPosition + 1) % 6;
	RedrawWindow();
}


void CGLDrawWnd::PrintText2D(float x, float y, const char * text)
{
	PrintText3D(x,y,0,text);
	Texts.GetTail().type = TextPortion::Text2D;
}


void CGLDrawWnd::PrintText3D(float x, float y, float z, const char * text)
{
	TextPortion tp;
	tp.text = text;
	tp.x = x;
	tp.y = y;
	tp.z = -z;
	tp.type = TextPortion::Text3D;
	Texts.AddTail(tp);
}

void CGLDrawWnd::PrintTextStruct(int nLineNo, int nXPos, const char * text)
{
	TextPortion tp;
	tp.text = text;
	tp.nLineNo = nLineNo;
	tp.nXPos = nXPos;
	tp.type = TextPortion::TextStruct;
	Texts.AddTail(tp);
}


// the scene may be rotating automatically with the timer, or a funny image may be animated

void CGLDrawWnd::OnTimer(UINT_PTR nIDEvent) 
{
	if(nIDEvent == ROTATION_TIMER_ID)
	{
		float m[16];
		glMatrixMode(GL_MODELVIEW);
		glGetFloatv(GL_MODELVIEW_MATRIX,m);
		glLoadIdentity();
		glRotated(0.7,0,1,0.1);
		glMultMatrixf(m);
		RedrawWindow();
		PerformVideoFramesRecord();
	}
	if(nIDEvent == SCALING_ANIMATE_TIMER_ID)
	{
		double t = virtual_animate_time;
		virtual_animate_time += 0.02;
		pWCDI->GetDOption(105) = old_scaling_factor * cos(2.*3.1415926535897932*t);

		//check if animation cycle is only once:
		if (pWCDI->GetIOption(153) && virtual_animate_time > 1.)
		{
			KillTimer(SCALING_ANIMATE_TIMER_ID);

		}
		else
		{
			Redraw();
			PerformVideoFramesRecord();
		}
	}
	CWnd::OnTimer(nIDEvent);
}

void CGLDrawWnd::OnButtonSaveImage() 
{
	// save a single image (screenshot)
	mystr path = pWCDI->GetTOption(120);
	mystr file = pWCDI->GetTOption(121);
	CString FileName_noext = CString(path) + "\\" + file;

	int radio_format = pWCDI->GetIOption(167);
	switch(radio_format)
	{
	case 0:	SaveWindowBitmap(this, FileName_noext + ".jpg", Gdiplus::ImageFormatJPEG); break;
	case 1:	SaveWindowBitmap(this, FileName_noext + ".png", Gdiplus::ImageFormatPNG); break;
	case 2:	SaveWindowBitmap(this, FileName_noext + ".bmp", Gdiplus::ImageFormatBMP); break;
	default: SaveWindowBitmap(this, FileName_noext + ".jpg", Gdiplus::ImageFormatJPEG); break;
	}
}

void CGLDrawWnd::OnButtonFramesRecording() 
{
	CRect r;
	GetWindowRect(&r);
	DialogFramesRecording.SetWindowDimensions(r);
	BOOL bOldStatusRecordFrames = DialogFramesRecording.m_bCheckRecordFrames;
	BOOL bOldStatusShowFrameNumber = DialogFramesRecording.m_bShowFrameNumbers;
  DialogFramesRecording.DoModal();
	if(		bOldStatusRecordFrames != DialogFramesRecording.m_bCheckRecordFrames ||
		bOldStatusShowFrameNumber != DialogFramesRecording.m_bShowFrameNumbers)
	 	RedrawWindow();
}

// here the logics is implemented to work with the automatically saved frames of the animation
void CGLDrawWnd::PerformVideoFramesRecord()
{
	if(!DialogFramesRecording.m_bCheckRecordFrames)
		return;

	static int counter = 0;
	if(++counter < DialogFramesRecording.m_nRecordEachFrameOf)
		return;
	
	//!AD: 2013-02-08   fix: prevent repeated saving of the same frame
	static double lastdrawtime = 0.;
	double drawtime = pWCDI->GetActualDrawTime();
	if (abs(drawtime-lastdrawtime) < abs(drawtime/1e8) )
	{
		return; // no more repeated draw of the same time point but allow for decreasing time value...
	}
	lastdrawtime = drawtime;

	counter = 0;

	CString FrameNumber;
	FrameNumber.Format("%d",DialogFramesRecording.m_nFrameCounter);
	DialogFramesRecording.m_nFrameCounter++;
	while(FrameNumber.GetLength() < FrameNumberDigits)
		FrameNumber = CString('0') + FrameNumber;

	// save a single image (screenshot)
	mystr path = pWCDI->GetTOption(122);
	mystr file = pWCDI->GetTOption(123);
	CString FileName_noext = CString(path) + "\\" + file + FrameNumber;
	
	int radio_format = pWCDI->GetIOption(167);
	switch(radio_format)
	{
		case 0:	SaveWindowBitmap(this, FileName_noext + ".jpg", Gdiplus::ImageFormatJPEG); break;
		case 1:	SaveWindowBitmap(this, FileName_noext + ".png", Gdiplus::ImageFormatPNG); break;
		case 2:	SaveWindowBitmap(this, FileName_noext + ".bmp", Gdiplus::ImageFormatBMP); break;
		default: SaveWindowBitmap(this, FileName_noext + ".jpg", Gdiplus::ImageFormatJPEG); break;
	}

return;

//AD 2012-11-09: new routine ends here...
// ONLY OLD CODE BELOW
//	FileName.Format("%d",DialogFramesRecording.m_nFrameCounter++);
//	while(FileName.GetLength() < FrameNumberDigits)
//		FileName = CString('0') + FileName;
//	FileName = DialogFramesRecording.m_strPathToImageFiles + "image" + FileName;
//
//	if(!SaveWindowBitmap(this,FileName + ".bmp",Gdiplus::ImageFormatJPEG))
//	{
//		AfxMessageBox("Errors happened while saving the image.\nMake sure the path for the image files exists.",
//			MB_OK | MB_ICONSTOP);
//		DialogFramesRecording.m_bCheckRecordFrames = FALSE;
//		RedrawWindow();
//	}
//	else
//		if(DialogFramesRecording.m_bProcessImage)
//		{
//			SetProgramDirectoryAsCurrent();
//
//			const char* str_bat = "process_image.bat";
//			const char* str_txt = "process_image.txt";
//
//#ifdef my_new_stdiostream
//			int opt = ios::in;
//			ifstream* ifs = new ifstream(str_bat, opt);
//			int isgood = ifs->good() && !ifs->fail();
//			ifs->close();
//			delete ifs;
//
//			ifstream* ifs2 = new ifstream(str_txt, opt);
//			int isgood2 = ifs2->good() && !ifs2->fail();
//			ifs2->close();
//			delete ifs2;
//#else
//			int opt = ios::in|ios::nocreate;
//			ifstream* ifs = new ifstream(str_bat, opt, filebuf::sh_read);
//			int isgood = ifs->good() && !ifs->fail();
//			ifs->close();
//			delete ifs;
//
//			ifstream* ifs2 = new ifstream(str_txt, opt, filebuf::sh_read);
//			int isgood2 = ifs2->good() && !ifs2->fail();
//			ifs2->close();
//			delete ifs2;
//#endif
//
//
//			if (!isgood) 
//			{
//				if (isgood2)
//				{
//					//copy file!
//					const int buflen = 4096;
//					char buf[buflen+1];
//					char ch;
//
//#ifdef my_new_stdiostream
//					ifstream* ifile = new ifstream(str_txt, opt);
//#else
//					ifstream* ifile = new ifstream(str_txt, opt, filebuf::sh_read);
//#endif
//
//					int i=0;
//					while (!ifile->eof() && i < buflen)
//					{
//						ifile->get(ch);
//						if (!ifile->eof() && ch != (char)EOF)
//						{
//							buf[i] = ch;
//							i++;
//						}
//					} 
//					buf[i] = (char)0;
//
//					ifile->close();
//					delete ifile;
//
//					if (i == buflen) 
//					{
//						AfxMessageBox("The file for processing images 'process_image.txt'\nis too long (> 4096 kb)!\nThe file has been recreated!");
//						isgood2=0;
//					}
//					else
//					{
//						ofstream ofs(str_bat);
//						ofs << buf;
//					}
//				}
//				if (!isgood2)
//				{
//					//file does not exist --> generate new file!
//					int opt = ios::out;
//					ofstream* ofs = new ofstream();
//#ifdef my_new_stdiostream
//					ofs->open(str_bat, opt);
//#else
//					ofs->open(str_bat, opt, filebuf::sh_read);
//#endif
//					(*ofs) << "@echo off\ndel %1.jpg\nC:\\Programme\\ImageMagick\\convert -quality 100 %1.bmp %1.jpg\ndel %1.bmp\n";
//					ofs->close();
//					delete ofs;
//				}
//			}
//
//
//			ShellExecute(
//				NULL,
//				"open",
//				str_bat,
//				FileName,
//				".",
//				SW_HIDE
//				//				SW_SHOWMINIMIZED
//				);
//		}
}

void CGLDrawWnd::Configuration2EDC(ElementDataContainer& edc)
{
	const int model_view_mat_len = 16;
	float m[model_view_mat_len];
	double v[model_view_mat_len];
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	for(int i = 0; i < model_view_mat_len; i++)
	{
		v[i] = (double)m[i];
	}

	edc.TreeSetVectorC("HOTINTConfiguration.OpenGL.model_view_matrix",v,model_view_mat_len,"OpenGL model view matrix");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.Center_point.xpos",(double)CenterPoint.x,"Centerpoint x-position");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.Center_point.ypos",(double)CenterPoint.y,"Centerpoint y-position");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.Center_offset.xpos",(double)CenterOffset.x,"Centeroffset x-position");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.Center_offset.ypos",(double)CenterOffset.y,"Centeroffset y-position");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.Center_offset.zpos",(double)CenterOffset.z,"Centeroffset z-position");

	edc.TreeSetIntC("HOTINTConfiguration.axes_position",AxesPosition,"position of axes: left, bottom, right, top, no axes");
	edc.TreeSetIntC("HOTINTConfiguration.lock_rotation",bNoRotationState,"lock rotation of model (for 2D models)");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.scene_offset",(double)SCENE_OFFSET_COEFF,"Scene offset");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.zoom_factor",(double)DetailLevel,"zoom factor of 3D view");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.zoom_factor_mouse_inc",(double)DetailLevelStepCoeff,"increment of zoom factor per mousemove of 3D view on mousemove");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.translation_mouse_inc",(double)TranslationStepCoeff,"increment of translation per mousemove of 3D view on mousemove");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.rotation_mouse_inc",(double)RotationStep,"increment of rotation per mousemove of 3D view on mousemove");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.perspective_mouse_inc",(double)PerspectiveStepCoeff,"increment of perspective per mousemove of 3D view on mousemove");

	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.aspect_ratio",(double)AspectRatio,"current aspect ratio in OpenGL view");
	edc.TreeSetDoubleC("HOTINTConfiguration.3DView.maximum_scene_coordinates",(double)MaxSceneCoord,"stored maximum scene coordinates in OpenGL view");

	//saved_modelview_matrix[0] = NoSavedModelviewMatrix;
	saved_modelview_matrix[0] = NoSavedModelviewMatrix; //for loading, restore

	DialogFramesRecording.Configuration2EDC(edc);
}

void CGLDrawWnd::EDC2Configuration(const ElementDataContainer& edc)
{
	const int model_view_mat_len = 16;

	double* vptr;
	int len;
	edc.TreeGetVector("HOTINTConfiguration.OpenGL.model_view_matrix",&vptr,len);
	if (len != model_view_mat_len)
	{
		AfxMessageBox("Warning: HOTINTConfiguration.OpenGL.model_view_matrix must have length 16! Loaded matrix is ignored!");
	}
	else
	{
		for(int i = 0; i < model_view_mat_len; i++)
		{
			saved_modelview_matrix[i] = (float)(vptr[i]);
		}
	}

	CenterPoint.x = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.Center_point.xpos");
	CenterPoint.y = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.Center_point.ypos");
	
	CenterOffset.x = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.Center_offset.xpos");
	CenterOffset.y = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.Center_offset.ypos");
	CenterOffset.z = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.Center_offset.zpos");
	
	AxesPosition = edc.TreeGetInt("HOTINTConfiguration.axes_position");
	bNoRotationState = (bool)edc.TreeGetInt("HOTINTConfiguration.lock_rotation");

	SCENE_OFFSET_COEFF = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.scene_offset");
	DetailLevel = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.zoom_factor");

	DetailLevelStepCoeff = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.zoom_factor_mouse_inc");
	TranslationStepCoeff = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.translation_mouse_inc");
	RotationStep = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.rotation_mouse_inc");
	PerspectiveStepCoeff = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.perspective_mouse_inc");

	AspectRatio = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.aspect_ratio");
	MaxSceneCoord = (float)edc.TreeGetDouble("HOTINTConfiguration.3DView.maximum_scene_coordinates");

	DialogFramesRecording.EDC2Configuration(edc);

}

//will be erased:
void CGLDrawWnd::SaveConfig(CArchive & ar)
{
	float m[16];
	glGetFloatv(GL_MODELVIEW_MATRIX,m);
	for(int i = 0; i < 16; i++)
		ar << m[i];

	ar << CenterPoint.x << CenterPoint.y << DetailLevel;// << bLighting; //lighting moved to ioption124
	ar << SCENE_OFFSET_COEFF;
	ar << AxesPosition;
	ar << bNoRotationState;

	ar << DetailLevel;
	ar << TranslationStepCoeff;
	ar << RotationStep;
	ar << DetailLevelStepCoeff; 
	ar << PerspectiveStepCoeff;
	ar << AspectRatio;
	ar << CenterOffset.x;
	ar << CenterOffset.y;
	ar << CenterOffset.z;
	ar << MaxSceneCoord;

	saved_modelview_matrix[0] = NoSavedModelviewMatrix;

	DialogFramesRecording.Serialize(ar);
}

//will be erased:
void CGLDrawWnd::LoadConfig(CArchive & ar)
{
	for(int i = 0; i < 16; i++)
		ar >> saved_modelview_matrix[i];

	ar >> CenterPoint.x >> CenterPoint.y >> DetailLevel;// >> bLighting;

	ar >> SCENE_OFFSET_COEFF;
	ar >> AxesPosition;
	ar >> bNoRotationState;

	ar >> DetailLevel;
	ar >> TranslationStepCoeff;
	ar >> RotationStep;
	ar >> DetailLevelStepCoeff; 
	ar >> PerspectiveStepCoeff;
	ar >> AspectRatio;
	ar >> CenterOffset.x;
	ar >> CenterOffset.y;
	ar >> CenterOffset.z;
	ar >> MaxSceneCoord;


	DialogFramesRecording.Serialize(ar);
}

void CGLDrawWnd::ChooseAxesColor()
{
	if( (BackgroundR + BackgroundG + BackgroundB) / 3 < 0.5f )
		::glColor3f(0.9f, 0.9f, 0.9f);
	else
		::glColor3f(0.1f, 0.1f, 0.1f);
}


void CGLDrawWnd::BeginLinesOnSolid()
{
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	float fact = 1.0001;
	gluPerspective((GLdouble)(MAX_FOVY/DetailLevel), AspectRatio,
		MaxSceneCoord*SCENE_OFFSET_COEFF, MaxSceneCoord*(fact*2.+SCENE_OFFSET_COEFF));
	//gluPerspective((GLdouble)(MAX_FOVY/DetailLevel), AspectRatio,
	//	MaxSceneCoord*SCENE_OFFSET_COEFF, MaxSceneCoord*(2.0001+SCENE_OFFSET_COEFF));
	glMatrixMode(GL_MODELVIEW);

	glDepthMask( FALSE );
}

void CGLDrawWnd::EndLinesOnSolid()
{
	glMatrixMode( GL_PROJECTION );
	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
	glDepthMask( TRUE );
}



