//#***************************************************************************************
//# filename:     MBSObjectFactory.cpp
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

#include "MBSObjectFactory.h"
#include "mbs_interface.h"
#include "element.h"
#include "material.h"
#include "sensorsspecific.h"
#include "geomelements.h"
#include "node.h"

// there is a global object factory instance, which needs to be created before the particular model objects
#pragma warning(disable:4074)
#pragma init_seg(compiler)
MBSObjectFactory objectFactory;

MBSObjectFactory::MBSObjectFactory(void)
{
	// all objects, that are added in this function are available for the parser
	objectTypeInfos.SetLen(OFCMaxVal-1);	

	// ATTENTION: the order of the objects of a class MUST coincide with the numbers in AddElement, AddLoad,...

	// add automatically generated elements
	AddObjectInfos_Auto();

	// add other objects manually
	AddObjectInfo(OFCSensor, "FVElementSensor", 0, "FVElementSensor.txt", "");		// flag not used yet, therefore = 0
	AddObjectInfo(OFCSensor, "ElementSensor", 0, "ElementSensor.txt", "");			// flag not used yet, therefore = 0
	AddObjectInfo(OFCSensor, "LoadSensor", 0, "LoadSensor.txt", "");				// flag not used yet, therefore = 0
	AddObjectInfo(OFCSensor, "MultipleSensor", 0, "", "");		// flag not used yet, therefore = 0
	AddObjectInfo(OFCSensor, "SystemSensor", 0, "SystemSensor.txt", "");				// flag not used yet, therefore = 0

	AddObjectInfo(OFCLoad, "GCLoad", 0, "GCLoad.txt", "");					// flag not used yet, therefore = 0
	AddObjectInfo(OFCLoad, "BodyLoad", 0, "", "");				// flag not used yet, therefore = 0
	AddObjectInfo(OFCLoad, "ForceVector2D", 0, "ForceVector2D.txt", "");		
	AddObjectInfo(OFCLoad, "Moment2D", TAENotInRelease, "", "");				
	AddObjectInfo(OFCLoad, "ForceVector3D", 0, "ForceVector3D.txt", "");		// flag not used yet, therefore = 0
	AddObjectInfo(OFCLoad, "MomentVector3D", 0, "", "");	// flag not used yet, therefore = 0
	AddObjectInfo(OFCLoad, "Gravity", 0, "Gravity.txt", "");					// flag not used yet, therefore = 0
	AddObjectInfo(OFCLoad, "SurfacePressure", 0, "SurfacePressure.txt", "");					// flag not used yet, therefore = 0


	AddObjectInfo(OFCGeomElement, "GeomMesh3D", 0, "", "");				// flag not used yet, therefore = 0
	AddObjectInfo(OFCGeomElement, "GeomCylinder3D", 0, "", "");		// flag not used yet, therefore = 0
	AddObjectInfo(OFCGeomElement, "GeomSphere3D", 0, "", "");			// flag not used yet, therefore = 0
	AddObjectInfo(OFCGeomElement, "GeomCube3D", 0, "", "");				// flag not used yet, therefore = 0
	AddObjectInfo(OFCGeomElement, "GeomOrthoCube3D", 0, "", "");	// flag not used yet, therefore = 0

	
	AddObjectInfo(OFCMaterial, "Material", 0, "Material.txt", "");				
	AddObjectInfo(OFCBeamProperties, "Beam2DProperties", TAENotInRelease, "", "");
	AddObjectInfo(OFCBeamProperties, "Beam3DProperties", 0, "Beam3DProperties.txt", "");

	AddObjectInfo(OFCNode, "Node3D", TAENotInRelease, "Node3D.txt", "");
	AddObjectInfo(OFCNode, "Node3DS1rot1", 0, "Node3DS1rot1.txt", "");
	AddObjectInfo(OFCNode, "Node3DS2S3", 0, "", "");
	//AddObjectInfo(OFCNode, "NodeLinearBeam3DTest", 0, "", "");
	AddObjectInfo(OFCNode, "Node3DRxyz", 0, "", "");
	AddObjectInfo(OFCNode, "Node2DS2", TAENotInRelease, "", "");
	AddObjectInfo(OFCNode, "Node3DR123", 0, "", "");
	AddObjectInfo(OFCNode, "Node3DS1S2", 0, "", "");
}

MBSObjectFactory::~MBSObjectFactory(void)
{
}

MBSObjectFactory::ObjectTypeInfo::ObjectTypeInfo(char * name, int flags, char * example, char * description)
{
	this->name = name;
	this->flags = flags;
	this->example = example;
	this->description = description;
}

int MBSObjectFactory::AddObject(MBSObjectFactoryClass objectClass, int objectTypeId)
{
	switch(objectClass)
	{
	case OFCElement: return AddElement(objectTypeId);		// AddElement is generated automatically
	case OFCSensor: return AddSensor(objectTypeId);
	case OFCNode: return AddNode(objectTypeId); 
	case OFCLoad: return AddLoad(objectTypeId);
	case OFCMaterial:	return AddMaterial(objectTypeId);
	case OFCBeamProperties:	return AddBeamProperties(objectTypeId);
	case OFCGeomElement: return AddGeomElement(objectTypeId);
	}
	assert(0);
	return -1;
}

int MBSObjectFactory::GetObjectTypeId(MBSObjectFactoryClass objectClass, const mystr & typeName)
{
	// a simple search needs to be done
	TArrayDynamic<ObjectTypeInfo> & classTypeInfos = objectTypeInfos(objectClass);
	for(int i = 1; i <= classTypeInfos.Length(); i++)
	{
		if(classTypeInfos(i).name == typeName)
			return i;
	}
	return -1;
}

void MBSObjectFactory::AddObjectInfo(MBSObjectFactoryClass objectClass, char * name, int flag, char * example, char * description)
{
	objectTypeInfos(objectClass).Add(ObjectTypeInfo(name, flag, example, description));
}

// below particular objects are created and added to mbs

int MBSObjectFactory::AddSensor(int objectTypeId)
{
	// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //sensor type
	{
	case 1: //FieldVariableElementSensor
		{
			FieldVariableElementSensor sens(mbs);
			sens.SetFVESPos3D(1,FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_x),Vector3D(0));
			return GetMBS()->AddSensor(&sens);
		}
	case 2: //SpecialValueElementSensor
		{
			SingleElementDataSensor sens(mbs);
			sens.SetSingleElementDataSensor(1,"");
			return GetMBS()->AddSensor(&sens);
		}
	case 3: //LoadSensor
		{
			LoadSensor sens(mbs);
			sens.SetLoadSensor(1);
			return GetMBS()->AddSensor(&sens);
		}
	case 4: //MultipleSensor
		{
			TArray<int> sensors;
			MultipleSensor sens(mbs);
			sens.SetMultipleSensor(sensors,"maximum");
			return GetMBS()->AddSensor(&sens);
		}
	case 5: //SystemSensor
		{
			SystemSensor sens(mbs);
			sens.SetSystemSensor(SystemSensor::None);
			return GetMBS()->AddSensor(&sens);
		}
	}
	assert(0);
	return -1;
}

int MBSObjectFactory::AddLoad(int objectTypeId)
{
	MBSLoad load;
	load.SetMBS(mbs);

	// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //load type
	{
	case 1: //GCLoad
		{
			load.SetGCLoad(0,1);
			return GetMBS()->AddLoad(load);
		}
	case 2: //BodyLoad
		{
			load.SetBodyLoad(0,1);
			return GetMBS()->AddLoad(load);
		}
	case 3: //ForceVector2D
		{
			load.SetForceVector2D(Vector2D(),Vector2D()); 
			return GetMBS()->AddLoad(load);
		}
	case 4: //Moment2D
		{
			load.SetMomentVector2D(0,Vector2D()); 
			return GetMBS()->AddLoad(load);
		}
	case 5: //ForceVector3D
		{
			load.SetForceVector3D(Vector3D(),Vector3D()); 
			return GetMBS()->AddLoad(load);
		}
	case 6: //MomentVector3D
		{
			load.SetMomentVector3D(0,Vector3D()); 
			return GetMBS()->AddLoad(load);
		}
	case 7: //Gravity
		{
			load.SetGravity(9.81,1); 
			return GetMBS()->AddLoad(load);
		}
	case 8: //SurfacePressure
		{
			load.SetSurfacePressure(0.,1); 
			return GetMBS()->AddLoad(load);
		}
	}
	assert(0);
	return -1;
}

int MBSObjectFactory::AddMaterial(int objectTypeId)
{
		// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //load type
	{
	case 1: //Material
		{
			Material mat(mbs);
			return GetMBS()->AddMaterial(mat);
		}
	}
	assert(0);
	return -1;
}

int MBSObjectFactory::AddBeamProperties(int objectTypeId)
{
		// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //load type
	{
	case 1: //Beam2DProperties
		{
			Beam2DProperties m2dp(mbs);
			return GetMBS()->AddMaterial(m2dp);
		}
	case 2: //Beam3DProperties
		{
			Beam3DProperties m3dp(mbs);
			return GetMBS()->AddMaterial(m3dp);
		}
	}
	assert(0);
	return -1;
}

int MBSObjectFactory::AddNode(int objectTypeId)
{
		// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //load type
	{
	case 1: //Node3D
		{
			Node3D n(mbs);
			return GetMBS()->AddNode(&n);
		}
	case 2: //ANCFNodeS1rot1_3D
		{
			ANCFNodeS1rot1_3D n(mbs);
			return GetMBS()->AddNode(&n);
		}
	case 3: //ANCFNodeS2S3
		{
			ANCFNodeS2S3_3D n(mbs);
			return GetMBS()->AddNode(&n);
		}
	//case 4: //NodeRot
	//	{
	//		NodeBeam3D n(mbs);
	//		return GetMBS()->AddNode(&n);
	//	}
	case 4: //node for Beam3D
		{
			Node3DRxyz n(mbs);
			return GetMBS()->AddNode(&n);
		}
	case 5: //ANCFNodeS2_2D
		{
			ANCFNodeS2_2D n(mbs);
			return GetMBS()->AddNode(&n);
		}
	case 6: // node for RotorBeamXAxis
		{
			Node3DR123 n(mbs);
			return GetMBS()->AddNode(&n);
		}
	case 7: //ANCFNodeS1S2
		{
			ANCFNodeS1S2_3D n(mbs);
			return GetMBS()->AddNode(&n);
		}
	}
	assert(0);
	return -1;
}

// taken from Default_ParserObjectValues
Vector3D defaultbodycol(0.2,0.2,0.8);
int defaultgeomtile = 16;

int MBSObjectFactory::AddGeomElement(int objectTypeId)
{
	// the order of this list MUST coincide with the order in MBSObjectFactory()
	switch(objectTypeId) //element type
	{
	case 1: //GeomMesh3D
		{
			GeomMesh3D ge(GetMBS(), 0, defaultbodycol);
			return GetMBS()->Add(ge);
		}
	case 2: //GeomCyl3D
		{
			GeomZyl3D ge(GetMBS(), Vector3D(0.,0.,0.), Vector3D(0.,0.,0.), 0, defaultgeomtile, defaultbodycol);
			return GetMBS()->Add(ge);
		}
	case 3: //GeomSphere3D
		{
			GeomSphere3D ge(GetMBS(), 0, Vector3D(0.,0.,0.), 0, defaultgeomtile, defaultbodycol);
			ge.SetDrawParam(Vector3D(0.001,16,0.));
			return GetMBS()->Add(ge);
		}
	case 4: //GeomCube3D
		{
			GeomCube3D ge(GetMBS(), Vector3D(0.,0.,0.), Vector3D(1.,1.,1.), defaultbodycol);
			return GetMBS()->Add(ge);
		}
	case 5: //GeomOrthoCube3D
		{
			GeomOrthoCube3D ge(GetMBS(), Vector3D(0.,0.,0.), Vector3D(1.,1.,1.), Matrix3D(1.), defaultbodycol);
			return GetMBS()->Add(ge);
		}
	}
	assert(0);
	return -1;
}