//#***************************************************************************************
//# filename:     ElementsAndModelsLibraryInterface.h
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

#include "mbs_interface.h"

// instances of the objects (object factory and models library)
// are created and destroyed within the dll automatically
// when the dll is loaded/unloaded

#define ELEMENTS_AND_MODELS_DLL_NAME "MBSElementsAndModels.dll"
#define OBJECT_FACTORY_ACCESS_FUNCTION_NAME "GetObjectFactory"		// no arguments
#define MODELS_LIBRARY_ACCESS_FUNCTION_NAME "GetModelsLibrary"		// no arguments

// manages, creates and adds to MBS different classes of objects:
// Elements
// Sensors
// Nodes
// Loads
// Materials
// GeomElements
// within each class different types may be available
// types are identified by int type_id
// classes are defined by

enum MBSObjectFactoryClass
{
	OFCElement = 1,
	OFCSensor,
	OFCNode,
	OFCLoad,
	OFCMaterial,
	OFCBeamProperties, //$ DR+PG 2013-01-21
	OFCGeomElement,
	OFCMaxVal	// DR: do not delete, always has to be the last entry! If you want to add new object classes, do it before (!) this entry!
};

typedef enum {TAEBody = 1, TAEflexible=2, TAEconstraint=4, TAE2D = 8, TAEspecial_connector=16, TAEinput_output = 32, TAENotInRelease=64,TAE1D = 128} TAddElementType;
//$ DR 2013-01-21: added TAENotInRelease
//$ DR 2013-06-19: added TAE1D

struct MBSObjectFactoryInterface
{
	// these functions add objects of a given class and type to mbs
	// the number in the corresponding array in mbs (numerical identificator) is returned
	// or -1 if no object could be added
	virtual int AddObject(MBSObjectFactoryClass objectClass, int objectTypeId) = 0;
	// available range types: maximal value of objectType for a given class
	virtual int GetAvailableTypesCount(MBSObjectFactoryClass objectClass) = 0;
	// name of an object type id and object type id for a given name (if exists, otherwise -1)
	virtual const mystr & GetTypeName(MBSObjectFactoryClass objectClass, int objectTypeId) = 0;
	virtual int GetObjectTypeId(MBSObjectFactoryClass objectClass, const mystr & typeName) = 0;
	// description for the documentation and for the user interface
	virtual const mystr & GetTypeDescription(MBSObjectFactoryClass objectClass, int objectTypeId) = 0;
	// flags - characteristic properties of the type (e.g. constraint or not, 2D/3D for elements)
	virtual int GetTypeFlags(MBSObjectFactoryClass objectClass, int objectTypeId) = 0;
	// example of the type in skript language 
	virtual mystr GetTypeExample(MBSObjectFactoryClass objectClass, int objectTypeId) = 0;
	// if this function returns 1, then all experimental elements and menus, etc. are removed from Hotint
	virtual int ExcludeExperimentalObjects() = 0;

	// an object factory needs mbs to be able to create objects and to add them to the collections
	virtual void SetMBS(MBS * mbs) = 0;
};

// manages models
struct MBSModelsLibraryInterface
{
	struct MBSModelInterface
	{
		virtual ~MBSModelInterface() {}
		virtual mystr & GetMBSModelName() = 0;
		virtual mystr & GetMBSModelDescription() = 0;
		virtual int CreateMBSModel(MBS * mbs) = 0;
		// the two functions below need to be overridden in case the model can re-initialize itself with different model data
		virtual bool HasMBSModelInitData() { return false; }
		virtual int InitializeMBSModelData(MBS * mbs) { return 0; }
	};
	virtual int GetModelsCount() = 0;
	virtual MBSModelInterface * GetModelInterface(int nModel) = 0;		// 1-based: nModel = 1 .. GetModelsCount()

	// global variables need to be shared with the kernel module;
	// here the pointers to them are set from the kernel to the client module
	virtual void SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel
		) = 0;
#ifdef gencnt
	virtual void SetGlobalVariablesPointers(
		UserOutputInterface * global_uo_kernel,
		TArray<double> * ptrTMtspent_kernel,
		TArray<double> * ptrTMtstart_kernel,
		int *global_genvec,
		int *global_genmat
		) = 0;
#endif
};