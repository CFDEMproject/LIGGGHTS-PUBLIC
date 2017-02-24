//#**************************************************************
//# filename:             ANCFBeamShear3DGetdeidDiffVar.h
//#
//# author:               Gruber, Nachbagauer
//#
//# generated:						
//# description:          
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
 

TArray<Vector3D> GetdeidDiffVar(int i, double dr1dy, double dr2dy, double dr3dy, double dr1dz, double dr2dz, double dr3dz) const
{
 TArray<Vector3D> val(12);
 switch(i)
 { 
 case 1:
 
val(1) = Vector3D(0,-dr3dz,dr2dz);

val(2) = Vector3D(dr3dz,0,-dr1dz);

val(3) = Vector3D(-dr2dz,dr1dz,0);

val(4) = Vector3D(0,dr3dy,-dr2dy);

val(5) = Vector3D(-dr3dy,0,dr1dy);

val(6) = Vector3D(dr2dy,-dr1dy,0);

 break; 
 case 2:
 
val(1) = Vector3D(Power(dr2dz,2) + Power(dr3dz,2),-(dr1dz*dr2dz),-(dr1dz*dr3dz));

val(2) = Vector3D(-(dr1dz*dr2dz),Power(dr1dz,2) + Power(dr3dz,2),-(dr2dz*dr3dz));

val(3) = Vector3D(-(dr1dz*dr3dz),-(dr2dz*dr3dz),Power(dr1dz,2) + Power(dr2dz,2));

val(4) = Vector3D(-(dr2dy*dr2dz) - dr3dy*dr3dz,2*dr1dz*dr2dy - dr1dy*dr2dz,2*dr1dz*dr3dy - dr1dy*dr3dz);

val(5) = Vector3D(-(dr1dz*dr2dy) + 2*dr1dy*dr2dz,-(dr1dy*dr1dz) - dr3dy*dr3dz,2*dr2dz*dr3dy - dr2dy*dr3dz);

val(6) = Vector3D(-(dr1dz*dr3dy) + 2*dr1dy*dr3dz,-(dr2dz*dr3dy) + 2*dr2dy*dr3dz,-(dr1dy*dr1dz) - dr2dy*dr2dz);

 break; 
 case 3:
 
val(1) = Vector3D(0,0,0);

val(2) = Vector3D(0,0,0);

val(3) = Vector3D(0,0,0);

val(4) = Vector3D(1,0,0);

val(5) = Vector3D(0,1,0);

val(6) = Vector3D(0,0,1);

 break; 
 default:
 break; 
 } 
 return val; 
}