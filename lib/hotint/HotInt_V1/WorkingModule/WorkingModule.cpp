//#**************************************************************
//# filename:							WorkingModule.cpp
//#
//# author:               Yury Vetyukov
//#
//# generated:						2003
//# description:          main file of the DLL
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
//#***************************************************************************************
 
#include "stdafx.h"
#include "windows.h"
#include "..\mbs_interface\ElementsAndModelsLibraryInterface.h"

HMODULE hElementsAndModelsModule = NULL;
MBSModelsLibraryInterface * pModelsLibrary;
MBSObjectFactoryInterface * pObjectFactory;

BOOL LoadElementsAndModelsModule()
{
	hElementsAndModelsModule = LoadLibraryA(ELEMENTS_AND_MODELS_DLL_NAME);
	if(hElementsAndModelsModule == NULL)
	{
		MessageBox(NULL, ELEMENTS_AND_MODELS_DLL_NAME, "Could not load the library", 0);
		return FALSE;
	}
	typedef MBSModelsLibraryInterface * (*GetMBSModelsLibraryInterfacePFN)();
	GetMBSModelsLibraryInterfacePFN models_pfn =
		(GetMBSModelsLibraryInterfacePFN)GetProcAddress(hElementsAndModelsModule, MODELS_LIBRARY_ACCESS_FUNCTION_NAME);
	typedef MBSObjectFactoryInterface * (*GetMBSObjectFactoryInterfacePFN)();
	GetMBSObjectFactoryInterfacePFN objects_pfn =
		(GetMBSObjectFactoryInterfacePFN)GetProcAddress(hElementsAndModelsModule, OBJECT_FACTORY_ACCESS_FUNCTION_NAME);
	if(!models_pfn || !objects_pfn)
	{
		FreeLibrary(hElementsAndModelsModule);	
		MessageBox(NULL, ELEMENTS_AND_MODELS_DLL_NAME, "Incorrect library version loaded", 0);
		return FALSE;
	}
	pModelsLibrary = (*models_pfn)();
	pObjectFactory = (*objects_pfn)();
	return TRUE;
}

BOOL UnloadElementsAndModelsModule()
{
	if(hElementsAndModelsModule != NULL)
		FreeLibrary(hElementsAndModelsModule);
	return TRUE;
}

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
			return LoadElementsAndModelsModule();
		case DLL_PROCESS_DETACH:
			return UnloadElementsAndModelsModule();
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
			break;
    }
    return TRUE;
}

#include "mbs.h"
#include "script_parser.h"

MultiBodySystem * TIMBS;

extern TArray<double> TMtspent;
extern TArray<double> TMtstart;
extern UserOutputInterface * global_uo;
#ifdef gencnt
extern int *genvec, *genmat;
#endif

extern "C"
{
	__declspec(dllexport) WCDInterface * CreateWCDObject()
	{
		TIMBS = new MultiBodySystem();
#ifdef gencnt
		pModelsLibrary->SetGlobalVariablesPointers(global_uo, &TMtspent, &TMtstart, genvec, genmat);
#else
		pModelsLibrary->SetGlobalVariablesPointers(global_uo, &TMtspent, &TMtstart);
#endif
		pObjectFactory->SetMBS(TIMBS);
		TIMBS->SetModelsLibrary(pModelsLibrary);
		TIMBS->EDCParser().SetObjectFactory(pObjectFactory);
		return TIMBS;
	}
}