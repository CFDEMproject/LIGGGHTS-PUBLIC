//#***************************************************************************************
//# filename:     MBSModelsLibrary.cpp
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

#include "MBSModelsLibrary.h"

// these pointers to global variables need to be set from the kernel
UserOutputInterface * global_uo;
extern TArray<double> * ptrTMtspent;
extern TArray<double> * ptrTMtstart;

#ifdef gencnt 
int *genvec;
int *genmat; 
#endif

void MBSModelsLibrary::SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel
		)
{
	global_uo = global_uo_kernel;
	ptrTMtspent = ptrTMtspent_kernel;
	ptrTMtstart = ptrTMtstart_kernel;
}
#ifdef gencnt
void MBSModelsLibrary::SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel,
		int *global_genvec,
		int *global_genmat
		)
{
	global_uo = global_uo_kernel;
	ptrTMtspent = ptrTMtspent_kernel;
	ptrTMtstart = ptrTMtstart_kernel;
	genvec = global_genvec;
	genmat = global_genmat;
}
#endif

MBSModelsLibrary::MBSModelsLibrary(void)
{
}

MBSModelsLibrary::~MBSModelsLibrary(void)
{
	for(int i = 1; i <= GetModelsCount(); i++)
		delete GetModelInterface(i);
}

void MBSModelsLibrary::AddModel(MBSModelInterface * model)
{
	modelsLibrary.models.Add(model);
}

// there is a global models library object, which needs to be created before the particular model objects
#pragma warning(disable:4074)
#pragma init_seg(compiler)
MBSModelsLibrary modelsLibrary;

// first we implement the old functionality with function pointers
struct ModelDataWrapperOldFnPointers : MBSModelsLibraryInterface::MBSModelInterface
{
	mystr modelName;
	mystr modelDescription;
	int(*function_ptr)(MBS* mbs);
	int(*function_ptr_init_modeldata)(MBS* mbs);

	virtual mystr & GetMBSModelName() { return modelName; }
	virtual mystr & GetMBSModelDescription() { return modelDescription; }
	virtual int CreateMBSModel(MBS * mbs) { return function_ptr(mbs); }
	virtual bool HasMBSModelInitData() { return function_ptr_init_modeldata != NULL; }
	virtual int InitializeMBSModelData(MBS * mbs) { return function_ptr_init_modeldata(mbs); }
};

void ModelFunctionAutoRegistration(int(*function_ptr)(MBS* mbs), const char* functionName, const char* description, int option,
																	 int(*function_ptr_init_modeldata)(MBS* mbs))
{
	ModelDataWrapperOldFnPointers * model = new ModelDataWrapperOldFnPointers;
	model->function_ptr = function_ptr;
	model->function_ptr_init_modeldata = function_ptr_init_modeldata;
	model->modelName = functionName;
	model->modelDescription = description;

	modelsLibrary.AddModel(model);
}