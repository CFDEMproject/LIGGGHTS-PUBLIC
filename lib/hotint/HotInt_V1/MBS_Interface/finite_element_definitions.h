//#**************************************************************
//#
//# filename:             finite_element_definitions.h
//#
//# author:               YV
//#
//# generated:						august 2011
//# description:          general definitions of types, which are relevant for finite elements in hotint
//#
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
//#**************************************************************

#pragma once

/// All existing types of finite elements should be listed here:
/// under this "keyword" they will be identified and will obtain an integration rule.
enum TFiniteElementType
{
	// 2D
	TFE_Quadrilateral = 1,
	TFE_Triangle,
	TFE_ThinPlate,					// ANCFThinPlate3D
	TFE_ThinBeam2D,					// ANCFBeamBE2D
	TFE_Beam2D,							// BeamShear2D
	// 3D
	TFE_Hexahedral,
	TFE_Tetrahedral,
	TFE_Prism,							//$ EK 2013-03-05 : FE prismatic type
	TFE_Pyramid,            //$ EK 2013-03-05 : FE pyramid type
	TFE_Quadrilateral_custom, 	//$ EK 2012-08-08 : FE Quad type with extended/custom functionality 
	TFE_Hexahedral_custom,			//$ EK 2012-09-14 : FE Hex type with extended/custom functionality 
	TFE_Tetrahedral_custom      //$ EK 2013-03-05 : FE Tet type with extended/custom functionality 
};

/// Possible statuses of geometric nonlinearity considered.
enum GeometricNonlinearityStatus
{
	GNS_Linear,								///< no geometrically nonlinear effects
	GNS_NonlinearSmallStrain,	///< large rotations/displacements, but local strains are small
	GNS_NonlinearLargeStrain	///< general case: arbitrary deformations
};