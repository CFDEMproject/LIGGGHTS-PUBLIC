//#**************************************************************
//# filename:             ANCFBeamShear3DGetddeidxidDiffVar.h
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
 

TArray<Vector3D> GetddeidxidDiffVar(int i, double dr1dy, double dr2dy, double dr3dy, double dr1dz, double dr2dz, double dr3dz, double ddr1dydxi, double ddr2dydxi, double ddr3dydxi, double ddr1dzdxi, double ddr2dzdxi, double ddr3dzdxi) const
{
 TArray<Vector3D> val(12);
 switch(i)
 { 
 case 1:
 
val(1) = Vector3D(0,-ddr3dzdxi,ddr2dzdxi);

val(2) = Vector3D(ddr3dzdxi,0,-ddr1dzdxi);

val(3) = Vector3D(-ddr2dzdxi,ddr1dzdxi,0);

val(4) = Vector3D(0,ddr3dydxi,-ddr2dydxi);

val(5) = Vector3D(-ddr3dydxi,0,ddr1dydxi);

val(6) = Vector3D(ddr2dydxi,-ddr1dydxi,0);

val(7) = Vector3D(0,-dr3dz,dr2dz);

val(8) = Vector3D(dr3dz,0,-dr1dz);

val(9) = Vector3D(-dr2dz,dr1dz,0);

val(10) = Vector3D(0,dr3dy,-dr2dy);

val(11) = Vector3D(-dr3dy,0,dr1dy);

val(12) = Vector3D(dr2dy,-dr1dy,0);

 break; 
 case 2:
 
val(1) = Vector3D(2*(ddr2dzdxi*dr2dz + ddr3dzdxi*dr3dz),-(ddr2dzdxi*dr1dz) - ddr1dzdxi*dr2dz,-(ddr3dzdxi*dr1dz) - ddr1dzdxi*dr3dz);

val(2) = Vector3D(-(ddr2dzdxi*dr1dz) - ddr1dzdxi*dr2dz,2*(ddr1dzdxi*dr1dz + ddr3dzdxi*dr3dz),-(ddr3dzdxi*dr2dz) - ddr2dzdxi*dr3dz);

val(3) = Vector3D(-(ddr3dzdxi*dr1dz) - ddr1dzdxi*dr3dz,-(ddr3dzdxi*dr2dz) - ddr2dzdxi*dr3dz,2*(ddr1dzdxi*dr1dz + ddr2dzdxi*dr2dz));

val(4) = Vector3D(-(ddr2dzdxi*dr2dy) - ddr2dydxi*dr2dz - ddr3dzdxi*dr3dy - ddr3dydxi*dr3dz,-(ddr2dzdxi*dr1dy) + 2*ddr2dydxi*dr1dz + 2*ddr1dzdxi*dr2dy - ddr1dydxi*dr2dz,-(ddr3dzdxi*dr1dy) + 2*ddr3dydxi*dr1dz + 2*ddr1dzdxi*dr3dy - ddr1dydxi*dr3dz);

val(5) = Vector3D(2*ddr2dzdxi*dr1dy - ddr2dydxi*dr1dz - ddr1dzdxi*dr2dy + 2*ddr1dydxi*dr2dz,-(ddr1dzdxi*dr1dy) - ddr1dydxi*dr1dz - ddr3dzdxi*dr3dy - ddr3dydxi*dr3dz,-(ddr3dzdxi*dr2dy) + 2*ddr3dydxi*dr2dz + 2*ddr2dzdxi*dr3dy - ddr2dydxi*dr3dz);

val(6) = Vector3D(2*ddr3dzdxi*dr1dy - ddr3dydxi*dr1dz - ddr1dzdxi*dr3dy + 2*ddr1dydxi*dr3dz,2*ddr3dzdxi*dr2dy - ddr3dydxi*dr2dz - ddr2dzdxi*dr3dy + 2*ddr2dydxi*dr3dz,-(ddr1dzdxi*dr1dy) - ddr1dydxi*dr1dz - ddr2dzdxi*dr2dy - ddr2dydxi*dr2dz);

val(7) = Vector3D(Power(dr2dz,2) + Power(dr3dz,2),-(dr1dz*dr2dz),-(dr1dz*dr3dz));

val(8) = Vector3D(-(dr1dz*dr2dz),Power(dr1dz,2) + Power(dr3dz,2),-(dr2dz*dr3dz));

val(9) = Vector3D(-(dr1dz*dr3dz),-(dr2dz*dr3dz),Power(dr1dz,2) + Power(dr2dz,2));

val(10) = Vector3D(-(dr2dy*dr2dz) - dr3dy*dr3dz,2*dr1dz*dr2dy - dr1dy*dr2dz,2*dr1dz*dr3dy - dr1dy*dr3dz);

val(11) = Vector3D(-(dr1dz*dr2dy) + 2*dr1dy*dr2dz,-(dr1dy*dr1dz) - dr3dy*dr3dz,2*dr2dz*dr3dy - dr2dy*dr3dz);

val(12) = Vector3D(-(dr1dz*dr3dy) + 2*dr1dy*dr3dz,-(dr2dz*dr3dy) + 2*dr2dy*dr3dz,-(dr1dy*dr1dz) - dr2dy*dr2dz);

 break; 
 case 3:
 
val(1) = Vector3D(0,0,0);

val(2) = Vector3D(0,0,0);

val(3) = Vector3D(0,0,0);

val(4) = Vector3D(0,0,0);

val(5) = Vector3D(0,0,0);

val(6) = Vector3D(0,0,0);

val(7) = Vector3D(0,0,0);

val(8) = Vector3D(0,0,0);

val(9) = Vector3D(0,0,0);

val(10) = Vector3D(1,0,0);

val(11) = Vector3D(0,1,0);

val(12) = Vector3D(0,0,1);

 break; 
 default:
 break; 
 } 
 return val; 
}
