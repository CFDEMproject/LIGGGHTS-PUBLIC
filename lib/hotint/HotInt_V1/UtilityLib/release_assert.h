//#**************************************************************
//#
//# filename:             release_assert.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						October 2010
//# description:					Definition of the preprocessor makro "release_assert"
//# remarks:		          In contrast to the makro "assert" from <assert.h>, "release_assert" will
//#                       stay active also in release mode, if the flag "__ASSERT_IN_RELEASE_MODE__"
//#                       is set. If not, "release_assert" is equivalent to "assert".
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
//#**************************************************************

#ifndef  __RELEASE_ASSERT_H__
#define  __RELEASE_ASSERT_H__


#include <assert.h>
#include "preprocessor_includes.h"


#undef release_assert

#ifndef  __ASSERT_IN_RELEASE_MODE__

#define release_assert(_Expression) (assert(_Expression))

#else  /*__ASSERT_IN_RELEASE_MODE__*/

#ifndef  NDEBUG

#define release_assert(_Expression) (assert(_Expression))

#else  /*NDEBUG*/

#ifdef  __cplusplus
extern "C" {
#endif

_CRTIMP void __cdecl _wassert(__in_z const wchar_t * _Message, __in_z const wchar_t *_File, __in unsigned _Line);

#ifdef  __cplusplus
}
#endif

#define release_assert(_Expression) (void)( (!!(_Expression)) || (_wassert(_CRT_WIDE(#_Expression), _CRT_WIDE(__FILE__), __LINE__), 0) )

#endif  /* NDEBUG */

#endif  /* __ASSERT_IN_RELEASE_MODE__ */




#endif  /*__RELEASE_ASSERT_H__*/