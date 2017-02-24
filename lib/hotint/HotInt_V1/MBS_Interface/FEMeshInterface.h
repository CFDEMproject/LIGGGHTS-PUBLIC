//#***************************************************************************************
//# filename:     FEMeshInterface.h
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

#include "finite_element_definitions.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ struct FEMesh_Settings
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ This structure and corresponding functions were added by YV 03.11.10
struct FEMesh_Settings
{
	//$ LA 2011-2-16: this status will be given to the elements when transfered to mbs
	// at the moment only set for 3D finite elements within FEMesh::GetPtrToMBSElement_new
	// if not set: default value is GNS_NonlinearLargeStrain
	// for FFRF elements internaly set to GNS_Linear
	GeometricNonlinearityStatus FEMesh_geometricNonlinearityStatus;
	// when "true" finite elements are created,
	// this flags determines the choice between normal finite elements (derived from FiniteElementXD)
	// and FFRF elements (derived from FiniteElement3DFFRF)
	bool generate_ffrf_elements;
	// this cms number will be passed to the created elements in the Set-function
	int cms_element_number;
	// domain to be used for the elements; by default is set to -1, which means
	// that the domain ids locally stored inside the FEElement structure should be used
	int domain;
	int flag_use_plate2D;  // 0 .. use quadrilateral, 1 .. use Plate2Dlin, 2 .. use Plate2Dlinhom
	double thickness2D; // thickness of 2D finite elements
	int flag_mbselement_warnings; // warnings for nonexisting elements can be turned off -this is used for Plate2dlinhom, funciton GetPtrToMBSElement_new is overloaded in derived class FEMesh3D_STA
	FEMesh_Settings() : FEMesh_geometricNonlinearityStatus(GNS_NonlinearLargeStrain), generate_ffrf_elements(false), cms_element_number(0), domain(-1), flag_use_plate2D(0), thickness2D(1), flag_mbselement_warnings(1) {}

	double& Thickness2D() { return thickness2D;}
	int& MBSElementWarnings() { return flag_mbselement_warnings; }
};

struct FEMeshInterface
{
	virtual int NElements() const = 0;
	virtual int NPoints() const = 0;
	virtual const Vector3D GetPoint3D(int i) const = 0;
	virtual Element * GetPtrToMBSElement_new(int elem_num) const = 0;
	virtual FEMesh_Settings & Settings() const = 0;
};