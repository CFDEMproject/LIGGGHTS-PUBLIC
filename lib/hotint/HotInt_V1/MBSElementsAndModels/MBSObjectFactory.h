//#***************************************************************************************
//# filename:     MBSObjectFactory.h
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

class MBSObjectFactory : public MBSObjectFactoryInterface
{
	struct ObjectTypeInfo
	{
		mystr name;						// name of the object, e.g. Rigid3D
		int flags;						// flags - characteristic properties of the type (e.g. constraint or not, 2D/3D for elements)
		mystr description;		// tex description for the docu
		mystr example;				// example file for the docu

		void operator=(const ObjectTypeInfo & oti)
		{
			name = oti.name;
			description = oti.description;
			flags = oti.flags;
			example = oti.example;
		}
		ObjectTypeInfo() : flags(0) {}
		ObjectTypeInfo(const ObjectTypeInfo & oti)
		{
			*this = oti;
		}
		ObjectTypeInfo(char * name, int flags, char * example, char * description);
	};

	// first index (outer array) - class
	// second index (inner array) - type id
	TArrayDynamic<TArrayDynamic<ObjectTypeInfo>> objectTypeInfos;

	MBS * mbs;

	// adding objects of particular classes
	int AddElement(int objectTypeId);
	int AddSensor(int objectTypeId);
	int AddLoad(int objectTypeId);
	int AddMaterial(int objectTypeId);
	int AddBeamProperties(int objectTypeId);
	int AddNode(int objectTypeId);
	int AddGeomElement(int objectTypeId);

	// initialization of object type infos is performed in constructor;
	// first the automatically generated part is called
	void AddObjectInfos_Auto();
	// this function adds a new object type info entry
	void AddObjectInfo(MBSObjectFactoryClass objectClass, char * name, int flag, char * example, char * description);

	MBS * GetMBS() { return mbs; }

public:
	MBSObjectFactory(void);
	~MBSObjectFactory(void);
	virtual void SetMBS(MBS * mbs) { this->mbs = mbs; }

	// these functions add objects of a given class and type to mbs
	// the number in the corresponding array in mbs (numerical identificator) is returned
	// or -1 if no object could be added
	virtual int AddObject(MBSObjectFactoryClass objectClass, int objectTypeId);
	// available range types: maximal value of objectType for a given class
	virtual int GetAvailableTypesCount(MBSObjectFactoryClass objectClass) { return objectTypeInfos(objectClass).Length(); }
	// name of an object type id and object type id for a given name (if exists, otherwise -1)
	virtual const mystr & GetTypeName(MBSObjectFactoryClass objectClass, int objectTypeId) { return objectTypeInfos(objectClass)(objectTypeId).name; }
	virtual int GetObjectTypeId(MBSObjectFactoryClass objectClass, const mystr & typeName);
	// description for the documentation and for the user interface
	virtual const mystr & GetTypeDescription(MBSObjectFactoryClass objectClass, int objectTypeId)  { return objectTypeInfos(objectClass)(objectTypeId).description; }
	// flags - characteristic properties of the type (e.g. constraint or not, 2D/3D for elements)
	virtual int GetTypeFlags(MBSObjectFactoryClass objectClass, int objectTypeId)  { return objectTypeInfos(objectClass)(objectTypeId).flags; }
	// example of the type in skript language 
	virtual mystr GetTypeExample(MBSObjectFactoryClass objectClass, int objectTypeId)  { return objectTypeInfos(objectClass)(objectTypeId).example; }
	// if this function returns 1, then all experimental elements are removed from Hotint
	virtual int ExcludeExperimentalObjects() { 
#ifdef __EXCLUDE_EXPERIMENTAL_OBJECTS__
		return 1;
#else
		return 0; 
#endif
	}

};

// there is a global object factory object
extern MBSObjectFactory objectFactory;