//#**************************************************************
//# filename:             ANCFBeamShear3DGetdeidxi.h
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
 
Vector3D Getdeidxi(int i, double dr1dy, double dr2dy, double dr3dy, double dr1dz, double dr2dz, double dr3dz, double ddr1dydxi, double ddr2dydxi, double ddr3dydxi, double ddr1dzdxi, double ddr2dzdxi, double ddr3dzdxi) const
{
 switch(i)
 { 
 case 1:
 return Vector3D(ddr3dzdxi*dr2dy - ddr3dydxi*dr2dz - ddr2dzdxi*dr3dy + ddr2dydxi*dr3dz,-(ddr3dzdxi*dr1dy) + ddr3dydxi*dr1dz + ddr1dzdxi*dr3dy - ddr1dydxi*dr3dz,ddr2dzdxi*dr1dy - ddr2dydxi*dr1dz - ddr1dzdxi*dr2dy + ddr1dydxi*dr2dz); 
 break; 
 case 2:
 return Vector3D(-(dr1dz*(ddr2dzdxi*dr2dy + ddr2dydxi*dr2dz + ddr3dzdxi*dr3dy + ddr3dydxi*dr3dz)) + 2*dr1dy*(ddr2dzdxi*dr2dz + ddr3dzdxi*dr3dz) - ddr1dzdxi*(dr2dy*dr2dz + dr3dy*dr3dz) + ddr1dydxi*(Power(dr2dz,2) + Power(dr3dz,2)),ddr2dydxi*Power(dr1dz,2) - ddr1dzdxi*dr1dy*dr2dz - dr1dz*(ddr2dzdxi*dr1dy - 2*ddr1dzdxi*dr2dy + ddr1dydxi*dr2dz) - ddr3dzdxi*dr2dz*dr3dy + 2*ddr3dzdxi*dr2dy*dr3dz - ddr3dydxi*dr2dz*dr3dz - ddr2dzdxi*dr3dy*dr3dz + ddr2dydxi*Power(dr3dz,2),ddr3dydxi*Power(dr1dz,2) - ddr3dzdxi*dr2dy*dr2dz + ddr3dydxi*Power(dr2dz,2) + 2*ddr2dzdxi*dr2dz*dr3dy - ddr1dzdxi*dr1dy*dr3dz - ddr2dzdxi*dr2dy*dr3dz - ddr2dydxi*dr2dz*dr3dz - dr1dz*(ddr3dzdxi*dr1dy - 2*ddr1dzdxi*dr3dy + ddr1dydxi*dr3dz)); 
 break; 
 case 3:
 return Vector3D(ddr1dzdxi,ddr2dzdxi,ddr3dzdxi); 
 break;
 default:
 break;
 } 
}