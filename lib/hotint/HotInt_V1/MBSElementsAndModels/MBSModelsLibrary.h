//#***************************************************************************************
//# filename:     MBSModelsLibrary.h
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

#pragma once
#include "elementsandmodelslibraryinterface.h"

class MBSModelsLibrary : public MBSModelsLibraryInterface
{
	TArray<MBSModelInterface*> models;

public:
	MBSModelsLibrary(void);
	~MBSModelsLibrary(void);

	// this function should be used by model objects, which implement MBSModelInterface, to register themselves
	static void AddModel(MBSModelInterface * model);

	virtual int GetModelsCount() { return models.Length(); }
	virtual MBSModelInterface * GetModelInterface(int nModel) { return models(nModel); }

	// global variables need to be shared with the kernel module;
	// here the pointers to them are set from the kernel to the client module
	virtual void SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel
		);
#ifdef gencnt
	virtual void SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel,
		int *global_genvec,
		int *global_genmat
		);
#endif
};

// there is a global models library object created after all model
extern MBSModelsLibrary modelsLibrary;