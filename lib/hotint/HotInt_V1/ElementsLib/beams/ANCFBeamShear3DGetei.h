//#**************************************************************
//# filename:             ANCFBeamShear3DGetei.h
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
 
Vector3D Getei(int i, double dr1dy, double dr2dy, double dr3dy, double dr1dz, double dr2dz, double dr3dz) const
{
 switch(i)
 { 
 case 1:
 return Vector3D(-(dr2dz*dr3dy) + dr2dy*dr3dz,dr1dz*dr3dy - dr1dy*dr3dz,-(dr1dz*dr2dy) + dr1dy*dr2dz); 
 break; 
 case 2:
 return Vector3D(-(dr1dz*(dr2dy*dr2dz + dr3dy*dr3dz)) + dr1dy*(Power(dr2dz,2) + Power(dr3dz,2)),Power(dr1dz,2)*dr2dy - dr1dy*dr1dz*dr2dz + dr3dz*(-(dr2dz*dr3dy) + dr2dy*dr3dz),Power(dr1dz,2)*dr3dy - dr1dy*dr1dz*dr3dz + dr2dz*(dr2dz*dr3dy - dr2dy*dr3dz)); 
 break; 
 case 3:
 return Vector3D(dr1dz,dr2dz,dr3dz); 
 break;
 default:
 break;
 } 
}

Vector3D GeteiP(int i, double dr1dy, double dr2dy, double dr3dy, double dr1dz, double dr2dz, double dr3dz, double drP1dy, double drP2dy, double drP3dy, double drP1dz, double drP2dz, double drP3dz) const
{
 switch(i)
 { 
 case 1:
 return Vector3D(-(drP2dz*dr3dy+dr2dz*drP3dy) + drP2dy*dr3dz + dr2dy*drP3dz, drP1dz*dr3dy + dr1dz*drP3dy - drP1dy*dr3dz - dr1dy*drP3dz,-(drP1dz*dr2dy + dr1dz*drP2dy) + drP1dy*dr2dz + dr1dy*drP2dz); 
 break; 
 case 2:
 return Vector3D(-(drP1dz*(dr2dy*dr2dz + dr3dy*dr3dz) + dr1dz*(drP2dy*dr2dz + dr2dy*drP2dz + drP3dy*dr3dz + dr3dy*drP3dz)) + drP1dy*(Power(dr2dz,2) + Power(dr3dz,2)) + dr1dy*(2*dr2dz*drP2dz + 2*dr3dz*drP3dz),2*dr1dz*drP1dz*dr2dy + Power(dr1dz,2)*drP2dy - drP1dy*dr1dz*dr2dz - dr1dy*drP1dz*dr2dz - dr1dy*dr1dz*drP2dz + drP3dz*(-(dr2dz*dr3dy) + dr2dy*dr3dz) + dr3dz*(-(drP2dz*dr3dy) - (dr2dz*drP3dy) + drP2dy*dr3dz + dr2dy*drP3dz),2*dr1dz*drP1dz*dr3dy + Power(dr1dz,2)*drP3dy - drP1dy*dr1dz*dr3dz - dr1dy*drP1dz*dr3dz - dr1dy*dr1dz*drP3dz + drP2dz*(dr2dz*dr3dy - dr2dy*dr3dz)  + dr2dz*(drP2dz*dr3dy + dr2dz*drP3dy - drP2dy*dr3dz - dr2dy*drP3dz)); 
 break; 
 case 3:
 return Vector3D(drP1dz,drP2dz,drP3dz); 
 break;
 default:
 break;
 }
}