//#***************************************************************************************
//# filename:     MBSElementsAndModels.cpp
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

// MBSElementsAndModels.cpp : Definiert den Einstiegspunkt fuer die DLL-Anwendung.
//

#include <windows.h>
#include "mbs_interface.h"
#include "mbsmodelslibrary.h"
#include "mbsobjectfactory.h"

#ifdef _MANAGED
#pragma managed(push, off)
#endif

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
		{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
    }
    return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif

extern "C"
{
	__declspec(dllexport) MBSModelsLibraryInterface * GetModelsLibrary()
	{
		return &modelsLibrary;
	}
	__declspec(dllexport) MBSObjectFactoryInterface * GetObjectFactory()
	{
		return &objectFactory;
	}
}
