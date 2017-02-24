//#**************************************************************
//#
//# filename:             Beam2DFFRF.cpp
//#
//# author:               Dibold Markus
//#
//# generated:						April 2006
//# description:          2D-plane stress plate element for floating frame of reference or absolute coordinates
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
//#***************************************************************************************
 
#include "rigid2d.h"
#include "femathhelperfunctions.h"
//#include "material.h"
#include "Node.h"
#include "referenceFrame2D.h"
#include "Beam2DFFRF.h"
//#include "graphicsconstants.h"
//#include "elementdataaccess.h"
//#include "solversettings_auto.h"
//#include "sensors\sensors.h"



void Beam2D::BuildDSMatrices() 
{
};

//ploc = -1 .. +1
double Beam2D::GetS0(double ploc, int shape) const
{
	if (ploc == -1 || ploc == 1) return 0;

	//if (GetMBS()->GetTime() > 0.0001) global_uo << "s";
	double r = ploc;

	if (ffrfmode==1)
	{
		//global_uo << "Eigenmodes, NS = " << NS() << "\n";

		switch(shape) //Eigenformen
		{
		case  1: return r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		case  2: return	-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		case  3: return	cosh(0.2365020372431352E1*r+0.2365020372431352E1)-1.0*cos(0.2365020372431352E1*r+0.2365020372431352E1)-0.9825022145762381*sinh(0.2365020372431352E1*r+0.2365020372431352E1)+0.9825022145762381*sin(0.2365020372431352E1*r+0.2365020372431352E1);
		case  4: return	cosh(0.3926602312047919E1*r+0.3926602312047919E1)-1.0*cos(0.3926602312047919E1*r+0.3926602312047919E1)-0.1000777311907269E1*sinh(0.3926602312047919E1*r+0.3926602312047919E1)+0.1000777311907269E1*sin(0.3926602312047919E1*r+0.3926602312047919E1);
		case  5: return	cosh(0.5497803919000835E1*r+0.5497803919000835E1)-1.0*cos(0.5497803919000835E1*r+0.5497803919000835E1)-0.9999664501254088*sinh(0.5497803919000835E1*r+0.5497803919000835E1)+0.9999664501254088*sin(0.5497803919000835E1*r+0.5497803919000835E1);
		case  6: return	cosh(0.7068582745628732E1*r+0.7068582745628732E1)-1.0*cos(0.7068582745628732E1*r+0.7068582745628732E1)-0.1000001449897656E1*sinh(0.7068582745628732E1*r+0.7068582745628732E1)+0.1000001449897656E1*sin(0.7068582745628732E1*r+0.7068582745628732E1);
		case  7: return	cosh(0.8639379828699741E1*r+0.8639379828699741E1)-1.0*cos(0.8639379828699741E1*r+0.8639379828699741E1)-0.9999999373443833*sinh(0.8639379828699741E1*r+0.8639379828699741E1)+0.9999999373443833*sin(0.8639379828699741E1*r+0.8639379828699741E1);
		case  8: return	cosh(0.1021017612281303E2*r+0.1021017612281303E2)-1.0*cos(0.1021017612281303E2*r+0.1021017612281303E2)-0.1000000002707595E1*sinh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1000000002707595E1*sin(0.1021017612281303E2*r+0.1021017612281303E2);
		default: return	0;
		}
	}
	else if (ffrfmode==2)
	{
		//global_uo << "approximated Eigenmodes 1, NS = " << NS() << "\n";

		switch(shape) //approximierte Eigenmoden
		{
		case  1: return r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		case  2: return	-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		case  3: return	0.1588077761E1-0.3176155521E1*r*r+0.1588077761E1*r*r*r*r;
		case  4: return	-0.5690568867E1*r+0.136095169E2*r*r*r-0.101473272E2*r*r*r*r*r+0.2228379167E1*r*r*r*r*r*r*r;
		case  5: return	-0.1405998458E1+0.2087060599E2*r*r-0.4843760064E2*r*r*r*r+0.3988737715E2*r*r*r*r*r*r-0.1091438404E2*r*r*r*r*r*r*r*r;
		case  6: return	0.9978018898E1*r-0.8304579388E2*r*r*r+0.2033363063E3*r*r*r*r*r-0.2264050907E3*r*r*r*r*r*r*r+0.1250943438E3*r*r*r*r*r*r*r*r*r-0.2895778429E2*pow(r,11.0);
		case  7: return	0.1414567736E1-0.5236059149E2*r*r+0.3200418442E3*r*r*r*r-0.7539913233E3*r*r*r*r*r*r+0.8705952324E3*r*r*r*r*r*r*r*r-0.5060683207E3*pow(r,10.0)+0.1203685911E3*pow(r,12.0);
		case  8: return	-0.1444444331E2*r+0.2510192916E3*r*r*r-0.1305754014E4*r*r*r*r*r+0.3200277873E4*r*r*r*r*r*r*r-0.4418988734E4*r*r*r*r*r*r*r*r*r+0.3634255529E4*pow(r,11.0)-0.1688890926E4*pow(r,13.0)+0.3425254233E3*pow(r,15.0);
		default: return	0;
		}
	}
	else if (ffrfmode==3)
	{
		//global_uo << "approximated Eigenmodes 2, NS = " << NS() << "\n";

		switch(shape) //approximierte Eigenmoden 2, Gerstm
		{
		case  1: return r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		case  2: return	-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		case  3: return	r*r*r*r-2.0*r*r+1.0;
		case  4: return	-r*r*r*r*r-r+2.0*r*r*r;
		case  5: return	-1.0+118.0/9.0*r*r+100.0/9.0*r*r*r*r*r*r-209.0/9.0*r*r*r*r;
		case  6: return	22016.0/3675.0*r-386048.0/11025.0*r*r*r+573952.0/11025.0*r*r*r*r*r-253952.0/11025.0*r*r*r*r*r*r*r;
		case  7: return	0;
		case  8: return	0;
		default: return	0;
		}
	}
	else if (ffrfmode==4)
	{
		//global_uo << "approximated Eigenmodes 1, NS = " << NS() << "\n";

		switch(shape) //approximierte Eigenmoden wie unter ffrfmode=2, nicht vereinfacht
		{
		case  1: return r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		case  2: return	-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		case  3: return	0.1588056235093594E1+0.2350484118576435E-1*r-0.3176112470187187E1*r*
										r-0.470096823715287E-1*r*r*r+0.1588056235093594E1*r*r*r*r+0.2350484118576435E-1
										*r*r*r*r*r;
		case  4: return	-0.6318476614386524E-3-0.5690568523498349E1*r+0.4861797896147454E-2*
										r*r+0.1360951607891885E2*r*r*r-0.782805280797895E-2*r*r*r*r
										-0.1014732658734265E2*r*r*r*r*r+0.3598102573270149E-2*r*r*r*r*r*r+
										0.2228379031922149E1*r*r*r*r*r*r*r;
		case  5: return	-0.1405998464277214E1+0.2087061002379698E2*r*r+
										0.6630567074313145E-13*r*r*r-0.48437608258235E2*r*r*r*r-0.1939551704933784E-12*
										r*r*r*r*r+0.3988738030218793E2*r*r*r*r*r*r+0.1889933287573624E-12*r*r*r*r*r*r*r
										-0.1091438360347269E2*r*r*r*r*r*r*r*r-0.6134382900711549E-13*r*r*r*r*r*r*r*r*r;
		case  6: return	-0.1285708230034643E-13+0.9978024577309942E1*r+
										0.4175591936322908E-12*r*r-0.8304589607938295E2*r*r*r-0.5078706583809296E-11*r*
										r*r*r+0.2033372322328371E3*r*r*r*r*r+0.1197797657429938E-10*r*r*r*r*r*r
										-0.2264073872733016E3*r*r*r*r*r*r*r-0.103257847601983E-10*r*r*r*r*r*r*r*r+
										0.125096539279074E3*r*r*r*r*r*r*r*r*r+0.3021812658376271E-11*pow(r,10.0)
										-0.289585127365365E2*pow(r,11.0);
		case  7: return	0.1414567511330691E1-0.5236069201914458E2*r*r+0.2660660977140584E-10
										*r*r*r+0.3200430116862467E3*r*r*r*r-0.3905843414857162E-9*r*r*r*r*r
										-0.7539969802118095E3*r*r*r*r*r*r+0.2099775437281532E-8*r*r*r*r*r*r*r+
										0.8706062489318276E3*r*r*r*r*r*r*r*r-0.4368369832462425E-8*r*r*r*r*r*r*r*r*r
										-0.5060775489454747E3*pow(r,10.0)+0.3866717670166089E-8*pow(r,11.0)+
										0.1203713930470238E3*pow(r,12.0)-0.1234145543270886E-8*pow(r,13.0);
		case  8: return	0.8630221265383098E-12-0.1444419679559728E2*r-0.7344540036389513E-10
										*r*r+0.251006768545024E3*r*r*r+0.2127669376223108E-8*r*r*r*r
										-0.1305463183338323E4*r*r*r*r*r-0.1824390736979506E-7*r*r*r*r*r*r+
										0.3198238771759708E4*r*r*r*r*r*r*r+0.6779662183255282E-7*r*r*r*r*r*r*r*r
										-0.4412513068171108E4*r*r*r*r*r*r*r*r*r-0.1158566211230247E-6*pow(r,10.0)+
										0.3624133763745404E4*pow(r,11.0)+0.9109529059375317E-7*pow(r,12.0)
										-0.1681298727025664E4*pow(r,13.0)-0.2684647093147202E-7*pow(r,14.0)+
										0.340339871280556E3*pow(r,15.0);
		default: return	0;
		}
	}
}

//ploc = -1 .. +1
void Beam2D::GetS0(Vector& sf, double ploc) const
{
	sf.SetLen(FlexDOF());

	if (ploc == -1 || ploc == 1)
	{
		sf.SetAll(0); return;
	}
	//if (GetMBS()->GetTime() > 0.0001) global_uo << "S";

	double r = ploc;
	if (ffrfmode==1)
	{
		//Eigenformen
		sf(1)=r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		sf(2)=-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		if (FlexDOF() >= 3) sf(3)=cosh(0.2365020372431352E1*r+0.2365020372431352E1)-1.0*cos(0.2365020372431352E1*r+0.2365020372431352E1)-0.9825022145762381*sinh(0.2365020372431352E1*r+0.2365020372431352E1)+0.9825022145762381*sin(0.2365020372431352E1*r+0.2365020372431352E1);
		if (FlexDOF() >= 4) sf(4)=cosh(0.3926602312047919E1*r+0.3926602312047919E1)-1.0*cos(0.3926602312047919E1*r+0.3926602312047919E1)-0.1000777311907269E1*sinh(0.3926602312047919E1*r+0.3926602312047919E1)+0.1000777311907269E1*sin(0.3926602312047919E1*r+0.3926602312047919E1);
		if (FlexDOF() >= 5) sf(5)=cosh(0.5497803919000835E1*r+0.5497803919000835E1)-1.0*cos(0.5497803919000835E1*r+0.5497803919000835E1)-0.9999664501254088*sinh(0.5497803919000835E1*r+0.5497803919000835E1)+0.9999664501254088*sin(0.5497803919000835E1*r+0.5497803919000835E1);
		if (FlexDOF() >= 6) sf(6)=cosh(0.7068582745628732E1*r+0.7068582745628732E1)-1.0*cos(0.7068582745628732E1*r+0.7068582745628732E1)-0.1000001449897656E1*sinh(0.7068582745628732E1*r+0.7068582745628732E1)+0.1000001449897656E1*sin(0.7068582745628732E1*r+0.7068582745628732E1);
		if (FlexDOF() >= 7) sf(7)=cosh(0.8639379828699741E1*r+0.8639379828699741E1)-1.0*cos(0.8639379828699741E1*r+0.8639379828699741E1)-0.9999999373443833*sinh(0.8639379828699741E1*r+0.8639379828699741E1)+0.9999999373443833*sin(0.8639379828699741E1*r+0.8639379828699741E1);
		if (FlexDOF() >= 8) sf(8)=cosh(0.1021017612281303E2*r+0.1021017612281303E2)-1.0*cos(0.1021017612281303E2*r+0.1021017612281303E2)-0.1000000002707595E1*sinh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1000000002707595E1*sin(0.1021017612281303E2*r+0.1021017612281303E2);
	}
	else if (ffrfmode==2)
	{
		//approximierte Eigenmoden
		sf(1)=r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		sf(2)=-pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		if (FlexDOF() >= 3) sf(3)=0.1588077761E1-0.3176155521E1*r*r+0.1588077761E1*r*r*r*r;
		if (FlexDOF() >= 4) sf(4)=-0.5690568867E1*r+0.136095169E2*r*r*r-0.101473272E2*r*r*r*r*r+0.2228379167E1*r*r*r*r*r*r*r;
		if (FlexDOF() >= 5) sf(5)=-0.1405998458E1+0.2087060599E2*r*r-0.4843760064E2*r*r*r*r+0.3988737715E2*r*r*r*r*r*r-0.1091438404E2*r*r*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6)=0.9978018898E1*r-0.8304579388E2*r*r*r+0.2033363063E3*r*r*r*r*r-0.2264050907E3*r*r*r*r*r*r*r+0.1250943438E3*r*r*r*r*r*r*r*r*r-0.2895778429E2*pow(r,11.0);
		if (FlexDOF() >= 7) sf(7)=0.1414567736E1-0.5236059149E2*r*r+0.3200418442E3*r*r*r*r-0.7539913233E3*r*r*r*r*r*r+0.8705952324E3*r*r*r*r*r*r*r*r-0.5060683207E3*pow(r,10.0)+0.1203685911E3*pow(r,12.0);
		if (FlexDOF() >= 8) sf(8)=-0.1444444331E2*r+0.2510192916E3*r*r*r-0.1305754014E4*r*r*r*r*r+0.3200277873E4*r*r*r*r*r*r*r-0.4418988734E4*r*r*r*r*r*r*r*r*r+0.3634255529E4*pow(r,11.0)-0.1688890926E4*pow(r,13.0)+0.3425254233E3*pow(r,15.0);
	}
	else if (ffrfmode==3)
	{
		//approximierte Eigenmoden 2, Gerstm
		sf(1)= r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		sf(2)= -pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		if (FlexDOF() >= 3) sf(3)= r*r*r*r-2.0*r*r+1.0;
		if (FlexDOF() >= 4) sf(4)= -r*r*r*r*r-r+2.0*r*r*r;
		if (FlexDOF() >= 5) sf(5)= -1.0+118.0/9.0*r*r+100.0/9.0*r*r*r*r*r*r-209.0/9.0*r*r*r*r;
		if (FlexDOF() >= 6) sf(6)= 22016.0/3675.0*r-386048.0/11025.0*r*r*r+573952.0/11025.0*r*r*r*r*r-253952.0/11025.0*r*r*r*r*r*r*r;
		if (FlexDOF() >= 7) sf(7)= 0;
		if (FlexDOF() >= 8) sf(8)= 0;
	}
	else if (ffrfmode==4)
	{
		//approximierte Eigenmoden wie unter ffrfmode=2, nicht vereinfacht
		sf(1)= r+1.0-pow(r+1.0,2.0)+pow(r+1.0,3.0)/4.0;
		sf(2)= -pow(r+1.0,2.0)/2.0+pow(r+1.0,3.0)/4.0;
		if (FlexDOF() >= 3) sf(3)= 0.1588056235093594E1+0.2350484118576435E-1*r-0.3176112470187187E1*r*
r-0.470096823715287E-1*r*r*r+0.1588056235093594E1*r*r*r*r+0.2350484118576435E-1
*r*r*r*r*r;
		if (FlexDOF() >= 4) sf(4)= -0.6318476614386524E-3-0.5690568523498349E1*r+0.4861797896147454E-2*
r*r+0.1360951607891885E2*r*r*r-0.782805280797895E-2*r*r*r*r
-0.1014732658734265E2*r*r*r*r*r+0.3598102573270149E-2*r*r*r*r*r*r+
0.2228379031922149E1*r*r*r*r*r*r*r;
		if (FlexDOF() >= 5) sf(5)= -0.1405998464277214E1+0.2087061002379698E2*r*r+
0.6630567074313145E-13*r*r*r-0.48437608258235E2*r*r*r*r-0.1939551704933784E-12*
r*r*r*r*r+0.3988738030218793E2*r*r*r*r*r*r+0.1889933287573624E-12*r*r*r*r*r*r*r
-0.1091438360347269E2*r*r*r*r*r*r*r*r-0.6134382900711549E-13*r*r*r*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6)= -0.1285708230034643E-13+0.9978024577309942E1*r+
0.4175591936322908E-12*r*r-0.8304589607938295E2*r*r*r-0.5078706583809296E-11*r*
r*r*r+0.2033372322328371E3*r*r*r*r*r+0.1197797657429938E-10*r*r*r*r*r*r
-0.2264073872733016E3*r*r*r*r*r*r*r-0.103257847601983E-10*r*r*r*r*r*r*r*r+
0.125096539279074E3*r*r*r*r*r*r*r*r*r+0.3021812658376271E-11*pow(r,10.0)
-0.289585127365365E2*pow(r,11.0);
		if (FlexDOF() >= 7) sf(7)= 0.1414567511330691E1-0.5236069201914458E2*r*r+0.2660660977140584E-10
*r*r*r+0.3200430116862467E3*r*r*r*r-0.3905843414857162E-9*r*r*r*r*r
-0.7539969802118095E3*r*r*r*r*r*r+0.2099775437281532E-8*r*r*r*r*r*r*r+
0.8706062489318276E3*r*r*r*r*r*r*r*r-0.4368369832462425E-8*r*r*r*r*r*r*r*r*r
-0.5060775489454747E3*pow(r,10.0)+0.3866717670166089E-8*pow(r,11.0)+
0.1203713930470238E3*pow(r,12.0)-0.1234145543270886E-8*pow(r,13.0);
		if (FlexDOF() >= 8) sf(8)= 0.8630221265383098E-12-0.1444419679559728E2*r-0.7344540036389513E-10
*r*r+0.251006768545024E3*r*r*r+0.2127669376223108E-8*r*r*r*r
-0.1305463183338323E4*r*r*r*r*r-0.1824390736979506E-7*r*r*r*r*r*r+
0.3198238771759708E4*r*r*r*r*r*r*r+0.6779662183255282E-7*r*r*r*r*r*r*r*r
-0.4412513068171108E4*r*r*r*r*r*r*r*r*r-0.1158566211230247E-6*pow(r,10.0)+
0.3624133763745404E4*pow(r,11.0)+0.9109529059375317E-7*pow(r,12.0)
-0.1681298727025664E4*pow(r,13.0)-0.2684647093147202E-7*pow(r,14.0)+
0.340339871280556E3*pow(r,15.0);
	}
}


//ploc = -1 .. +1
double Beam2D::GetDS0(double ploc, int shape) const
{

	if (ploc == -1) 
	{
		if (shape == 1) return 1;
		else return 0;
	} 
	else if (ploc == 1) 
	{
		if (shape == 2) return 1;
		else return 0;
	}

	double r = ploc;
	if (ffrfmode==1)
	{
		switch(shape) //Eigenformen
		{//D_sf/D_r:
		case 1: return -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		case 2: return -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		case 3: return 0.2365020372431352E1*sinh(0.2365020372431352E1*r+0.2365020372431352E1)+0.2365020372431352E1*sin(0.2365020372431352E1*r+0.2365020372431352E1)-0.2323637753431723E1*cosh(0.2365020372431352E1*r+0.2365020372431352E1)+0.2323637753431723E1*cos(0.2365020372431352E1*r+0.2365020372431352E1);
		case 4: return 0.3926602312047919E1*sinh(0.3926602312047919E1*r+0.3926602312047919E1)+0.3926602312047919E1*sin(0.3926602312047919E1*r+0.3926602312047919E1)-0.3929654506780184E1*cosh(0.3926602312047919E1*r+0.3926602312047919E1)+0.3929654506780184E1*cos(0.3926602312047919E1*r+0.3926602312047919E1);
		case 5: return 0.5497803919000835E1*sinh(0.5497803919000835E1*r+0.5497803919000835E1)+0.5497803919000835E1*sin(0.5497803919000835E1*r+0.5497803919000835E1)-0.5497619468368826E1*cosh(0.5497803919000835E1*r+0.5497803919000835E1)+0.5497619468368826E1*cos(0.5497803919000835E1*r+0.5497803919000835E1);
		case 6: return 0.7068582745628732E1*sinh(0.7068582745628732E1*r+0.7068582745628732E1)+0.7068582745628732E1*sin(0.7068582745628732E1*r+0.7068582745628732E1)-0.706859299435029E1*cosh(0.7068582745628732E1*r+0.7068582745628732E1)+0.706859299435029E1*cos(0.7068582745628732E1*r+0.7068582745628732E1);
		case 7: return 0.8639379828699741E1*sinh(0.8639379828699741E1*r+0.8639379828699741E1)+0.8639379828699741E1*sin(0.8639379828699741E1*r+0.8639379828699741E1)-0.863937928739407E1*cosh(0.8639379828699741E1*r+0.8639379828699741E1)+0.863937928739407E1*cos(0.8639379828699741E1*r+0.8639379828699741E1);
		case 8: return 0.1021017612281303E2*sinh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1021017612281303E2*sin(0.1021017612281303E2*r+0.1021017612281303E2)-0.1021017615045805E2*cosh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1021017615045805E2*cos(0.1021017612281303E2*r+0.1021017612281303E2);
		default: return 0;
		} 
	}
	else if (ffrfmode==2)
	{
		switch(shape) //approximierte Eigenformen
		{//D_sf/D_r:
		case 1: return -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		case 2: return -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		case 3: return -0.6352311042E1*r+0.6352311044E1*r*r*r;
		case 4: return -0.5690568867E1+0.408285507E2*r*r-0.50736636E2*r*r*r*r+0.15598654169E2*r*r*r*r*r*r;
		case 5: return 0.4174121198E2*r-0.19375040256E3*r*r*r+0.2393242629E3*r*r*r*r*r-0.8731507232E2*r*r*r*r*r*r*r;
		case 6: return 0.9978018898E1-0.24913738164E3*r*r+0.10166815315E4*r*r*r*r-0.15848356349E4*r*r*r*r*r*r+0.11258490942E4*r*r*r*r*r*r*r*r-0.31853562719E3*pow(r,10.0);
		case 7: return -0.10472118298E3*r+0.12801673768E4*r*r*r-0.45239479398E4*r*r*r*r*r+0.69647618592E4*r*r*r*r*r*r*r-0.5060683207E4*r*r*r*r*r*r*r*r*r+0.14444230932E4*pow(r,11.0);
		case 8: return -0.1444444331E2+0.7530578748E3*r*r-0.652877007E4*r*r*r*r+0.22401945111E5*r*r*r*r*r*r-0.39770898606E5*r*r*r*r*r*r*r*r+0.39976810819E5*pow(r,10.0)-0.21955582038E5*pow(r,12.0)+0.51378813495E4*pow(r,14.0);
		default: return 0;
		} 
	}
	else if (ffrfmode==3)
	{
		switch(shape) //approximierte Eigenformen 2, Gerstm
		{//D_sf/D_r:
		case 1: return -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		case 2: return -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		case 3: return 4.0*r*r*r-4.0*r;
		case 4: return -5.0*r*r*r*r-1.0+6.0*r*r;
		case 5: return 236.0/9.0*r+200.0/3.0*r*r*r*r*r-836.0/9.0*r*r*r;
		case 6: return 22016.0/3675.0-386048.0/3675.0*r*r+573952.0/2205.0*r*r*r*r-253952.0/1575.0*r*r*r*r*r*r;
		case 7: return 0;
		case 8: return 0;
		default: return 0;
		} 
	}
	else if (ffrfmode==4)
	{
		switch(shape) //approximierte Eigenformen wie unter ffrfmode=2, nicht vereinfacht
		{//D_sf/D_r:
		case 1: return -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		case 2: return -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		case 3: return 0.2350484118576435E-1-0.6352224940374374E1*r-0.1410290471145861*r*r+
0.6352224940374374E1*r*r*r+0.1175242059288218*r*r*r*r;
		case 4: return -0.5690568523498349E1+0.9723595792294907E-2*r+0.4082854823675654E2*r
*r-0.313122112319158E-1*r*r*r-0.5073663293671323E2*r*r*r*r+
0.2158861543962089E-1*r*r*r*r*r+0.1559865322345504E2*r*r*r*r*r*r;
		case 5: return 0.4174122004759395E2*r+0.1989170122293943E-12*r*r-0.19375043303294E3
*r*r*r-0.9697758524668919E-12*r*r*r*r+0.2393242818131276E3*r*r*r*r*r+
0.1322953301301537E-11*r*r*r*r*r*r-0.8731506882778152E2*r*r*r*r*r*r*r
-0.5520944610640394E-12*r*r*r*r*r*r*r*r;
		case 6: return 0.9978024577309942E1+0.8351183872645816E-12*r-0.2491376882381488E3*r
*r-0.2031482633523719E-10*r*r*r+0.1016686161164186E4*r*r*r*r+
0.7186785944579626E-10*r*r*r*r*r-0.1584851710913111E4*r*r*r*r*r*r
-0.8260627808158637E-10*r*r*r*r*r*r*r+0.1125868853511666E4*r*r*r*r*r*r*r*r+
0.3021812658376271E-10*r*r*r*r*r*r*r*r*r-0.3185436401019016E3*pow(r,10.0);
		case 7: return -0.1047213840382892E3*r+0.7981982931421753E-10*r*r+
0.1280172046744987E4*r*r*r-0.1952921707428581E-8*r*r*r*r-0.4523981881270857E4*r
*r*r*r*r+0.1469842806097072E-7*r*r*r*r*r*r+0.6964849991454621E4*r*r*r*r*r*r*r
-0.3931532849216182E-7*r*r*r*r*r*r*r*r-0.5060775489454747E4*r*r*r*r*r*r*r*r*r+
0.4253389437182698E-7*pow(r,10.0)+0.1444456716564286E4*pow(r,11.0)
-0.1604389206252152E-7*pow(r,12.0);
		case 8: return -0.1444419679559728E2-0.1468908007277903E-9*r+0.7530203056350721E3*r
*r+0.8510677504892433E-8*r*r*r-0.6527315916691617E4*r*r*r*r
-0.1094634442187704E-6*r*r*r*r*r+0.2238767140231796E5*r*r*r*r*r*r+
0.5423729746604226E-6*r*r*r*r*r*r*r-0.3971261761353997E5*r*r*r*r*r*r*r*r
-0.1158566211230247E-5*r*r*r*r*r*r*r*r*r+0.3986547140119945E5*pow(r,10.0)+
0.1093143487125038E-5*pow(r,11.0)-0.2185688345133363E5*pow(r,12.0)
-0.3758505930406082E-6*pow(r,13.0)+0.510509806920834E4*pow(r,14.0);
		default: return 0;
		} 
	}
}

//ploc = -1 .. +1
void Beam2D::GetDSMatrix0(Vector& sf, double ploc) const
{
	sf.SetLen(FlexDOF());
	double r = ploc;

	if (ffrfmode==1)
	{	
		//Eigenformen
		sf(1) = -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		sf(2) = -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		if (FlexDOF() >= 3) sf(3) =0.2365020372431352E1*sinh(0.2365020372431352E1*r+0.2365020372431352E1)+0.2365020372431352E1*sin(0.2365020372431352E1*r+0.2365020372431352E1)-0.2323637753431723E1*cosh(0.2365020372431352E1*r+0.2365020372431352E1)+0.2323637753431723E1*cos(0.2365020372431352E1*r+0.2365020372431352E1); 
		if (FlexDOF() >= 4) sf(4) =0.3926602312047919E1*sinh(0.3926602312047919E1*r+0.3926602312047919E1)+0.3926602312047919E1*sin(0.3926602312047919E1*r+0.3926602312047919E1)-0.3929654506780184E1*cosh(0.3926602312047919E1*r+0.3926602312047919E1)+0.3929654506780184E1*cos(0.3926602312047919E1*r+0.3926602312047919E1); 
		if (FlexDOF() >= 5) sf(5) =0.5497803919000835E1*sinh(0.5497803919000835E1*r+0.5497803919000835E1)+0.5497803919000835E1*sin(0.5497803919000835E1*r+0.5497803919000835E1)-0.5497619468368826E1*cosh(0.5497803919000835E1*r+0.5497803919000835E1)+0.5497619468368826E1*cos(0.5497803919000835E1*r+0.5497803919000835E1); 
		if (FlexDOF() >= 6) sf(6) =0.7068582745628732E1*sinh(0.7068582745628732E1*r+0.7068582745628732E1)+0.7068582745628732E1*sin(0.7068582745628732E1*r+0.7068582745628732E1)-0.706859299435029E1*cosh(0.7068582745628732E1*r+0.7068582745628732E1)+0.706859299435029E1*cos(0.7068582745628732E1*r+0.7068582745628732E1); 
		if (FlexDOF() >= 7) sf(7) =0.8639379828699741E1*sinh(0.8639379828699741E1*r+0.8639379828699741E1)+0.8639379828699741E1*sin(0.8639379828699741E1*r+0.8639379828699741E1)-0.863937928739407E1*cosh(0.8639379828699741E1*r+0.8639379828699741E1)+0.863937928739407E1*cos(0.8639379828699741E1*r+0.8639379828699741E1); 
		if (FlexDOF() >= 8) sf(8) =0.1021017612281303E2*sinh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1021017612281303E2*sin(0.1021017612281303E2*r+0.1021017612281303E2)-0.1021017615045805E2*cosh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1021017615045805E2*cos(0.1021017612281303E2*r+0.1021017612281303E2); 
	}
	else if (ffrfmode==2)
	{
		//approximierte Eigenformen
		sf(1) = -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		sf(2) = -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		if (FlexDOF() >= 3) sf(3) =-0.6352311042E1*r+0.6352311044E1*r*r*r;
		if (FlexDOF() >= 4) sf(4) =-0.5690568867E1+0.408285507E2*r*r-0.50736636E2*r*r*r*r+0.15598654169E2*r*r*r*r*r*r;
		if (FlexDOF() >= 5) sf(5) =0.4174121198E2*r-0.19375040256E3*r*r*r+0.2393242629E3*r*r*r*r*r-0.8731507232E2*r*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6) =0.9978018898E1-0.24913738164E3*r*r+0.10166815315E4*r*r*r*r-0.15848356349E4*r*r*r*r*r*r+0.11258490942E4*r*r*r*r*r*r*r*r-0.31853562719E3*pow(r,10.0);
		if (FlexDOF() >= 7) sf(7) =-0.10472118298E3*r+0.12801673768E4*r*r*r-0.45239479398E4*r*r*r*r*r+0.69647618592E4*r*r*r*r*r*r*r-0.5060683207E4*r*r*r*r*r*r*r*r*r+0.14444230932E4*pow(r,11.0);
		if (FlexDOF() >= 8) sf(8) =-0.1444444331E2+0.7530578748E3*r*r-0.652877007E4*r*r*r*r+0.22401945111E5*r*r*r*r*r*r-0.39770898606E5*r*r*r*r*r*r*r*r+0.39976810819E5*pow(r,10.0)-0.21955582038E5*pow(r,12.0)+0.51378813495E4*pow(r,14.0);
	}
	else if (ffrfmode==3)
	{
		//approximierte Eigenformen 2, Gerstm
		sf(1) = -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		sf(2) = -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		if (FlexDOF() >= 3) sf(3) = 4.0*r*r*r-4.0*r;
		if (FlexDOF() >= 4) sf(4) =	-5.0*r*r*r*r-1.0+6.0*r*r;
		if (FlexDOF() >= 5) sf(5) =	236.0/9.0*r+200.0/3.0*r*r*r*r*r-836.0/9.0*r*r*r;
		if (FlexDOF() >= 6) sf(6) =	22016.0/3675.0-386048.0/3675.0*r*r+573952.0/2205.0*r*r*r*r-253952.0/1575.0*r*r*r*r*r*r;
		if (FlexDOF() >= 7) sf(7) =	0;
		if (FlexDOF() >= 8) sf(8) =	0;
	}
  else if (ffrfmode==4)
	{
		//approximierte Eigenformen wie unter ffrfmode=2, nicht vereinfacht
		sf(1) = -1.0-2.0*r+3.0/4.0*pow(r+1.0,2.0);
		sf(2) = -r-1.0+3.0/4.0*pow(r+1.0,2.0);
		if (FlexDOF() >= 3) sf(3) = 0.2350484118576435E-1-0.6352224940374374E1*r-0.1410290471145861*r*r+
0.6352224940374374E1*r*r*r+0.1175242059288218*r*r*r*r;
		if (FlexDOF() >= 4) sf(4) =	-0.5690568523498349E1+0.9723595792294907E-2*r+0.4082854823675654E2*r
*r-0.313122112319158E-1*r*r*r-0.5073663293671323E2*r*r*r*r+
0.2158861543962089E-1*r*r*r*r*r+0.1559865322345504E2*r*r*r*r*r*r;
		if (FlexDOF() >= 5) sf(5) =	0.4174122004759395E2*r+0.1989170122293943E-12*r*r-0.19375043303294E3
*r*r*r-0.9697758524668919E-12*r*r*r*r+0.2393242818131276E3*r*r*r*r*r+
0.1322953301301537E-11*r*r*r*r*r*r-0.8731506882778152E2*r*r*r*r*r*r*r
-0.5520944610640394E-12*r*r*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6) =	0.9978024577309942E1+0.8351183872645816E-12*r-0.2491376882381488E3*r
*r-0.2031482633523719E-10*r*r*r+0.1016686161164186E4*r*r*r*r+
0.7186785944579626E-10*r*r*r*r*r-0.1584851710913111E4*r*r*r*r*r*r
-0.8260627808158637E-10*r*r*r*r*r*r*r+0.1125868853511666E4*r*r*r*r*r*r*r*r+
0.3021812658376271E-10*r*r*r*r*r*r*r*r*r-0.3185436401019016E3*pow(r,10.0);
		if (FlexDOF() >= 7) sf(7) =	-0.1047213840382892E3*r+0.7981982931421753E-10*r*r+
0.1280172046744987E4*r*r*r-0.1952921707428581E-8*r*r*r*r-0.4523981881270857E4*r
*r*r*r*r+0.1469842806097072E-7*r*r*r*r*r*r+0.6964849991454621E4*r*r*r*r*r*r*r
-0.3931532849216182E-7*r*r*r*r*r*r*r*r-0.5060775489454747E4*r*r*r*r*r*r*r*r*r+
0.4253389437182698E-7*pow(r,10.0)+0.1444456716564286E4*pow(r,11.0)
-0.1604389206252152E-7*pow(r,12.0);
		if (FlexDOF() >= 8) sf(8) =	-0.1444419679559728E2-0.1468908007277903E-9*r+0.7530203056350721E3*r
*r+0.8510677504892433E-8*r*r*r-0.6527315916691617E4*r*r*r*r
-0.1094634442187704E-6*r*r*r*r*r+0.2238767140231796E5*r*r*r*r*r*r+
0.5423729746604226E-6*r*r*r*r*r*r*r-0.3971261761353997E5*r*r*r*r*r*r*r*r
-0.1158566211230247E-5*r*r*r*r*r*r*r*r*r+0.3986547140119945E5*pow(r,10.0)+
0.1093143487125038E-5*pow(r,11.0)-0.2185688345133363E5*pow(r,12.0)
-0.3758505930406082E-6*pow(r,13.0)+0.510509806920834E4*pow(r,14.0);
	}
}

//ploc = -1 .. +1
void Beam2D::GetDDSMatrix0(Vector& sf, double ploc) const
{
	sf.SetLen(FlexDOF());
	double r = ploc;

	if (ffrfmode==1)
	{
		sf(1) = -1.0/2.0+3.0/2.0*r; 
		sf(2) = 1.0/2.0+3.0/2.0*r;								
		if (FlexDOF() >= 3) sf(3) = 0.5593321362015331E1*cosh(0.2365020372431352E1*r+0.2365020372431352E1)+0.5593321362015331E1*cos(0.2365020372431352E1*r+0.2365020372431352E1)-0.5495450625016643E1*sinh(0.2365020372431352E1*r+0.2365020372431352E1)-0.5495450625016643E1*sin(0.2365020372431352E1*r+0.2365020372431352E1);
		if (FlexDOF() >= 4) sf(4) = 0.1541820571698006E2*cosh(0.3926602312047919E1*r+0.3926602312047919E1)+0.1541820571698006E2*cos(0.3926602312047919E1*r+0.3926602312047919E1)-0.1543019047187259E2*sinh(0.3926602312047919E1*r+0.3926602312047919E1)-0.1543019047187259E2*sin(0.3926602312047919E1*r+0.3926602312047919E1);
		if (FlexDOF() >= 5) sf(5) = 0.3022584793178094E2*cosh(0.5497803919000835E1*r+0.5497803919000835E1)+0.3022584793178094E2*cos(0.5497803919000835E1*r+0.5497803919000835E1)-0.3022483385837342E2*sinh(0.5497803919000835E1*r+0.5497803919000835E1)-0.3022483385837342E2*sin(0.5497803919000835E1*r+0.5497803919000835E1);
		if (FlexDOF() >= 6) sf(6) = 0.4996486203180022E2*cosh(0.7068582745628732E1*r+0.7068582745628732E1)+0.4996486203180022E2*cos(0.7068582745628732E1*r+0.7068582745628732E1)-0.4996493447573659E2*sinh(0.7068582745628732E1*r+0.7068582745628732E1)-0.4996493447573659E2*sin(0.7068582745628732E1*r+0.7068582745628732E1);
		if (FlexDOF() >= 7) sf(7) = 0.7463888382454396E2*cosh(0.8639379828699741E1*r+0.8639379828699741E1)+0.7463888382454396E2*cos(0.8639379828699741E1*r+0.8639379828699741E1)-0.7463887914799867E2*sinh(0.8639379828699741E1*r+0.8639379828699741E1)-0.7463887914799867E2*sin(0.8639379828699741E1*r+0.8639379828699741E1);
		if (FlexDOF() >= 8) sf(8) = 0.1042476964588613E3*cosh(0.1021017612281303E2*r+0.1021017612281303E2)+0.1042476964588613E3*cos(0.1021017612281303E2*r+0.1021017612281303E2)-0.1042476967411219E3*sinh(0.1021017612281303E2*r+0.1021017612281303E2)-0.1042476967411219E3*sin(0.1021017612281303E2*r+0.1021017612281303E2);
	}
	else if (ffrfmode==2)
	{
		//approximierte Eigenformen
		sf(1) = -1.0/2.0+3.0/2.0*r; 
		sf(2) = 1.0/2.0+3.0/2.0*r;								
		if (FlexDOF() >= 3) sf(3) = -0.6352311042E1+0.19056933132E2*r*r;
		if (FlexDOF() >= 4) sf(4) = 0.816571014E2*r-0.202946544E3*r*r*r+0.93591925014E2*r*r*r*r*r;
		if (FlexDOF() >= 5) sf(5) = 0.4174121198E2-0.58125120768E3*r*r+0.11966213145E4*r*r*r*r-0.61120550624E3*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6) = -0.49827476328E3*r+0.4066726126E4*r*r*r-0.95090138094E4*r*r*r*r*r+0.90067927536E4*r*r*r*r*r*r*r-0.31853562719E4*r*r*r*r*r*r*r*r*r;
		if (FlexDOF() >= 7) sf(7) = -0.10472118298E3+0.38405021304E4*r*r-0.22619739699E5*r*r*r*r+0.487533330144E5*r*r*r*r*r*r-0.45546148863E5*r*r*r*r*r*r*r*r+0.158886540252E5*pow(r,10.0);
		if (FlexDOF() >= 8) sf(8) = 0.15061157496E4*r-0.2611508028E5*r*r*r+0.134411670666E6*r*r*r*r*r-0.318167188848E6*r*r*r*r*r*r*r+0.39976810819E6*r*r*r*r*r*r*r*r*r-0.263466984456E6*pow(r,11.0)+0.71930338893E5*pow(r,13.0);
	}
	else if (ffrfmode==3)
	{
		//approximierte Eigenformen 2, Gerstm
		sf(1) = -1.0/2.0+3.0/2.0*r;
		sf(2) = 1.0/2.0+3.0/2.0*r;							
		if (FlexDOF() >= 3) sf(3) = 12.0*r*r-4.0;
		if (FlexDOF() >= 4) sf(4) = -20.0*r*r*r+12.0*r;
		if (FlexDOF() >= 5) sf(5) = 236.0/9.0+1000.0/3.0*r*r*r*r-836.0/3.0*r*r;
		if (FlexDOF() >= 6) sf(6) = -772096.0/3675.0*r+2295808.0/2205.0*r*r*r-507904.0/525.0*r*r*r*r*r;
		if (FlexDOF() >= 7) sf(7) = 0;
		if (FlexDOF() >= 8) sf(8) = 0;
	}
	else if (ffrfmode==4)
	{
		//approximierte Eigenformen wie unter ffrfmode=2, nicht vereinfacht 
		sf(1) = -1.0/2.0+3.0/2.0*r;
		sf(2) = 1.0/2.0+3.0/2.0*r;							
		if (FlexDOF() >= 3) sf(3) = -0.6352224940374374E1-0.2820580942291722*r+0.1905667482112312E2*r*r+
0.470096823715287*r*r*r;
		if (FlexDOF() >= 4) sf(4) = 0.9723595792294907E-2+0.8165709647351308E2*r-0.939366336957474E-1*r*
r-0.2029465317468529E3*r*r*r+0.1079430771981045*r*r*r*r+0.9359191934073026E2*r*
r*r*r*r;
		if (FlexDOF() >= 5) sf(5) = 0.4174122004759395E2+0.3978340244587887E-12*r-0.58125129909882E3*r*r
-0.3879103409867568E-11*r*r*r+0.1196621409065638E4*r*r*r*r+
0.7937719807809222E-11*r*r*r*r*r-0.6112054817944707E3*r*r*r*r*r*r
-0.4416755688512316E-11*r*r*r*r*r*r*r;
		if (FlexDOF() >= 6) sf(6) = 0.8351183872645816E-12-0.4982753764762977E3*r-0.6094447900571156E-10
*r*r+0.4066744644656742E4*r*r*r+0.3593392972289813E-9*r*r*r*r
-0.9509110265478668E4*r*r*r*r*r-0.5782439465711046E-9*r*r*r*r*r*r+
0.9006950828093331E4*r*r*r*r*r*r*r+0.2719631392538644E-9*r*r*r*r*r*r*r*r
-0.3185436401019016E4*r*r*r*r*r*r*r*r*r;
		if (FlexDOF() >= 7) sf(7) = -0.1047213840382892E3+0.1596396586284351E-9*r+0.384051614023496E4*r*
r-0.7811686829714325E-8*r*r*r-0.2261990940635429E5*r*r*r*r+
0.8819056836582435E-7*r*r*r*r*r+0.4875394994018235E5*r*r*r*r*r*r
-0.3145226279372946E-6*r*r*r*r*r*r*r-0.4554697940509272E5*r*r*r*r*r*r*r*r+
0.4253389437182698E-6*r*r*r*r*r*r*r*r*r+0.1588902388220714E5*pow(r,10.0)
-0.1925267047502582E-6*pow(r,11.0);
		if (FlexDOF() >= 8) sf(8) = -0.1468908007277903E-9+0.1506040611270144E4*r+0.255320325146773E-7*r
*r-0.2610926366676647E5*r*r*r-0.5473172210938518E-6*r*r*r*r+
0.1343260284139078E6*r*r*r*r*r+0.3796610822622958E-5*r*r*r*r*r*r
-0.3177009409083198E6*r*r*r*r*r*r*r-0.1042709590107222E-4*r*r*r*r*r*r*r*r+
0.3986547140119945E6*r*r*r*r*r*r*r*r*r+0.1202457835837542E-4*pow(r,10.0)
-0.2622826014160036E6*pow(r,11.0)-0.4886057709527907E-5*pow(r,12.0)+
0.7147137296891676E5*pow(r,13.0);
	}
}

//ploc = -L/2 .. +L/2
Vector2D Beam2D::GetPos2D(const Vector2D& p_loc) const
{
	return GetPos2Drel(p_loc);
}

//ploc = -L/2 .. +L/2
Vector2D Beam2D::GetVel2D(const Vector2D& p_loc) const
{
	return GetVel2Drel(p_loc);
}

//ploc = -L/2 .. +L/2
Vector2D Beam2D::GetPos2DD(const Vector2D& p_loc) const
{
	return GetPos2DrelD(p_loc);
}

//ploc = -L/2 .. +L/2
Vector2D Beam2D::GetVel2DD(const Vector2D& p_loc) const
{
	return GetVel2DrelD(p_loc);
}
/////////////////////////////////////////////////////////////////////
////////////RELATIVE POSITIONS AND VELOCITIES////////////////////////
/////////////////////////////////////////////////////////////////////

//ploc = -L/2 .. +L/2
//relative position!
Vector2D Beam2D::GetPos2Drel(const Vector2D& p_loc) const
{
	double p0 = 2.*p_loc.X()/beaml;
	if (p0 == -1 || p0 == 1) return Vector2D(p_loc.X(), 0);
	double w = 0;

	for (int j = 1; j <= FlexDOF(); j++)
	{
		w += GetS0(p0, j) * XG(j);
	}

	return Vector2D(p_loc.X(), -w);
};

//ploc = -L/2 .. +L/2
//relative velocity!
Vector2D Beam2D::GetVel2Drel(const Vector2D& p_loc) const
{
	double p0 = 2.*p_loc.X()/beaml;
	if (p0 == -1 || p0 == 1) return Vector2D(0, 0);
	double wp = 0;

	for (int j = 1; j <= FlexDOF(); j++)
	{
		wp += GetS0(p0, j) * XGP(j);
	}

	return Vector2D(0, -wp);
};

//ploc = -L/2 .. +L/2
//relative position!
Vector2D Beam2D::GetPos2DrelD(const Vector2D& p_loc) const
{
	double w = 0;

	for (int j = 1; j <= FlexDOF(); j++)
	{
		w += GetS0(2.*p_loc.X()/beaml, j) * XGD(j);
	}

	return Vector2D(p_loc.X(), -w);
};

//ploc = -L/2 .. +L/2
//relative velocity!
Vector2D Beam2D::GetVel2DrelD(const Vector2D& p_loc) const
{
	double wp = 0;

	for (int j = 1; j <= FlexDOF(); j++)
	{
		wp += GetS0(2.*p_loc.X()/beaml, j) * XGPD(j);
	}

	return Vector2D(0, -wp);
}







//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//floating frame of reference formulation: for FFRF elements

//fill in sos x sos components, m might be larger
//compute stiffness matrix
void Beam2D::StiffnessMatrix(Matrix& m) 
{
	int ns = FlexDOF();

	m.SetSize(ns, ns);
	m.SetAll(0);

	Matrix B(1,ns);
	Vector temp;
	Matrix temp2;

	GetIntegrationRule(x1,w1,orderxy);

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		GetDDSMatrix0(temp, x1(i1)); //d^2 w / dr^2 --> mult. with 4/L^2
		temp *= 4./Sqr(beaml);
		B.SetRowVec(temp, 1);

		double fact = w1(i1) * EI * beaml * 0.5;
		m += fact * (B.GetTp() * B);
	}

	//UO() << "Klin=" << m << "\n";
}


void Beam2D::GetH(Matrix& H) 
{
	if (Hmatrix.Getrows() == FlexDOF())
	{
		H = Hmatrix;
		return;
	}
	else
	{
		int ns = FlexDOF();

		H.SetSize(ns, 1);
		H.SetAll(0);

		Vector temp(ns);
		Matrix mtemp(ns,1);

		GetIntegrationRule(x1,w1,orderxyH);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			GetS0(temp, x1(i1));
			mtemp.SetColVec(temp, 1);

			double fact = w1(i1) * beaml * 0.5;
			H += fact * mtemp;
		}

		Hmatrix = H;
		//UO() << "H=" << H << "\n";
	}
}


void Beam2D::EvalM(Matrix& m, double t) 
{
	if (massmatrix.Getcols() == FlexDOF())
	{
		m = massmatrix;
		return;
	}
	else
	{
		//...............
		UO() << "Error: wrong EvalM called\n";

		massmatrix = m;

	}
};


void Beam2D::EvalF2(Vector& f, double t) 
{
	//..........................
}; 

void Beam2D::DrawElement() 
{
	mbs->SetColor(col);

	int tile = 32;
	double eps = 1e-8;
	Vector2D p1, p2, p1p, p2p, n1, n2;
	for (int i = 1; i <= tile; i++)
	{
		double x1 = ((double)(i-1.)/(double)tile-0.5)*beaml;
		double x2 = ((double)(i   )/(double)tile-0.5)*beaml;
		double x1p = x1 + eps;
		double x2p = x2 + eps;

		p1 = GetPos2DD(Vector2D(x1,0));
		p2 = GetPos2DD(Vector2D(x2,0));
		p1p = GetPos2DD(Vector2D(x1p,0));
		p2p = GetPos2DD(Vector2D(x2p,0));
		n1 = (p1p-p1);
		n2 = (p2p-p2);
		n1.Normalize();
		n2.Normalize();
		n1 = Vector2D(-n1.Y(), n1.X());
		n2 = Vector2D(-n2.Y(), n2.X());

		GetMBS()->DrawQuad(ToP3D(p1-0.5*beamh*n1),ToP3D(p2-0.5*beamh*n2),ToP3D(p2+0.5*beamh*n2),ToP3D(p1+0.5*beamh*n1));
	}

};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/////////////////////////////////////
//2-node beam element with variable number of internal shape functions,  no axial deformation:
/////////////////////////////////////
Beam2DFFRF::Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, double L0, double H0, double T0, 
											 double rhoi, double EmIi, const Vector3D& coli, int isCMSi):Beam2D(mbsi),
											 K(),	SbarS(), Ibar12S(), I1S()
{
	//Add FFRF reference
	FFRFind = FFRFindex;

	beaml = L0;
	beamh = H0;
	beamt = T0;

	EI = EmIi;
	//UO() << "Bending Stiffness = " << EI << "\n";
	rho = rhoi;
	rhoA = rho*H0*T0;

	col = coli;
	if (isCMSi) AddType(TCMS);

	nfm = nfmi;
	ffrfmode = 1;

	mass = rhoA * beaml; //for EvalM()

	xg = Vector(SOS()); //for intial conditions
	xgd = xg; //for initial drawing

	x_init = Vector(2*SOS());
	x_init.SetAll(0);

	Vector rf_xinit = GetMBS()->GetElement(FFRFindex).GetXInit();
	x_init(SOSowned()+1) = rf_xinit(1);
	x_init(SOSowned()+2) = rf_xinit(2);
	x_init(SOSowned()+3) = rf_xinit(3);
	x_init(SOS()+SOSowned()+1) = rf_xinit(3+1);
	x_init(SOS()+SOSowned()+2) = rf_xinit(3+2);
	x_init(SOS()+SOSowned()+3) = rf_xinit(3+3);


	orderxy = 17-2; //9+2+6; //only for nonlinear
	orderxyM = 17-2; //6+2+6;
	orderxyH = 17-2; //4+2+6;

}


Beam2DFFRF::Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, int ffrfmodei, double L0, double H0, double T0, 
double rhoi, double EmIi, const Vector3D& coli, int isCMSi):Beam2D(mbsi),
K(),	SbarS(), Ibar12S(), I1S()
{
//Add FFRF reference
FFRFind = FFRFindex;

beaml = L0;
beamh = H0;
beamt = T0;

EI = EmIi;
//UO() << "Bending Stiffness = " << EI << "\n";
rho = rhoi;
rhoA = rho*H0*T0;

col = coli;
if (isCMSi) AddType(TCMS);

nfm = nfmi;
ffrfmode = ffrfmodei;

mass = rhoA * beaml; //for EvalM()

xg = Vector(SOS()); //for intial conditions
xgd = xg; //for initial drawing

x_init = Vector(2*SOS());
x_init.SetAll(0);

Vector rf_xinit = GetMBS()->GetElement(FFRFindex).GetXInit();
x_init(SOSowned()+1) = rf_xinit(1);
x_init(SOSowned()+2) = rf_xinit(2);
x_init(SOSowned()+3) = rf_xinit(3);
x_init(SOS()+SOSowned()+1) = rf_xinit(3+1);
x_init(SOS()+SOSowned()+2) = rf_xinit(3+2);
x_init(SOS()+SOSowned()+3) = rf_xinit(3+3);


orderxy = 17-2; //9+2+6; //only for nonlinear
orderxyM = 17-2; //6+2+6;
orderxyH = 17-2; //4+2+6;

}

Beam2DFFRF::Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, int ffrfmodei, double L0, double H0, double T0, 
double H20, double T20, double rhoi, double EmIi, const Vector3D& coli, int isCMSi):Beam2D(mbsi),
K(),	SbarS(), Ibar12S(), I1S()
{
	//to be used for hollow cross sections
//Add FFRF reference
FFRFind = FFRFindex;

beaml = L0;
beamh = H0;
beamt = T0;

beamhi = H20;
beamti = T20;

EI = EmIi;
//UO() << "Bending Stiffness = " << EI << "\n";
rho = rhoi;
rhoA = rho*(H0*T0-H20*T20);

col = coli;
if (isCMSi) AddType(TCMS);

nfm = nfmi;
ffrfmode = ffrfmodei;

mass = rhoA * beaml; //for EvalM()

xg = Vector(SOS()); //for intial conditions
xgd = xg; //for initial drawing

x_init = Vector(2*SOS());
x_init.SetAll(0);

Vector rf_xinit = GetMBS()->GetElement(FFRFindex).GetXInit();
x_init(SOSowned()+1) = rf_xinit(1);
x_init(SOSowned()+2) = rf_xinit(2);
x_init(SOSowned()+3) = rf_xinit(3);
x_init(SOS()+SOSowned()+1) = rf_xinit(3+1);
x_init(SOS()+SOSowned()+2) = rf_xinit(3+2);
x_init(SOS()+SOSowned()+3) = rf_xinit(3+3);


orderxy = 17-2; //9+2+6; //only for nonlinear
orderxyM = 17-2; //6+2+6;
orderxyH = 17-2; //4+2+6;
}


void Beam2DFFRF::Initialize() 
{
	Body2D::Initialize();
	int ns = FlexDOF();

	//Store beaml for testval
	testval = beaml;

	//compute mass of element:
	mass = beaml * rhoA;
	StiffnessMatrix(K);

	//+++++++++++++++++++++++++++++++++++++++++
	//compute mass matrix --> stored in "massmatrix"
	Matrix tmp;
	EvalMff(tmp,0);

	//stored functions:
	GetSbar(SbarS);
	GetIbarkl(1,2,Ibar12S);
	GetI1(I1S);
	I1122S = (GetIkl(1,1)+GetIkl(2,2));

	//char str[200];
	//sprintf(str, "%.24g, %.24g, %.24g, %.24g", SbarS(1), SbarS(2), SbarS(3), SbarS(4));

	//UO() << "I1S=" << I1S << "\nSbarS=" << str << "\n";

	//Vector2D tttA(0,0);
	//UO() << "InitPos = " << GetPos2D(tttA) << "\n";

}

void Beam2DFFRF::LinkToElements()
{
	//order of DOF in LTG: 0 .. 1*FlexDOF() Xref Yref phi

	if (SOSowned() != SOS())
	{
		TArray<int> storeltg(FlexDOF()*2);
		for (int i=1; i <= FlexDOF()*2; i++)
		{
			storeltg.Add(LTG(i));
		}
		LTGreset();

		//q_f:
		for (int i=1; i <= FlexDOF(); i++)
		{
			AddLTG(storeltg(i));
		}

		//reference frame:
		for (int i = 1; i <= ReferenceFrame().SOS(); i++)
		{
			//UO() << "RF-ltg" << i << "=" << ReferenceFrame().LTG(i) << "\n";
			AddLTG(ReferenceFrame().LTG(i));
		}

		//q_f_dot:
		for (int i=1; i <= FlexDOF(); i++)
		{
			AddLTG(storeltg(i+FlexDOF()));
		}

		//reference frame_dot
		for (int i = 1; i <= ReferenceFrame().SOS(); i++)
		{
			AddLTG(ReferenceFrame().LTG(i+ReferenceFrame().SOS()));
		}
	}
}


//insert all 9 blocks of mass matrix
void Beam2DFFRF::EvalM(Matrix& m, double t)
{
	//UO() << "Mff=" << massmatrix << "\n";
	//UO() << "S~=" << Sbar_tilde << "\n";

	static Vector temp;
	static Vector temp2;
	//static Matrix mtemp;

	int off = FlexDOF(); //offset where rigid body entries start!
	//m.SetAll(0);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRR = Unit(2,2) * mass
	m(1+off,1+off) = GetMass();
	m(2+off,2+off) = GetMass();
	m(1+off,2+off) = 0.;
	m(2+off,1+off) = 0.;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRtheta = A_theta*[I1+Sbar*qf]
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.);

	for (int j = 1; j <= FlexDOF(); j++)
	{
		pt.Y() += SbarS(j)*XG(j);
	}

	I1 += pt;
	Matrix3D ADphi = ReferenceFrame().GetRotMatrixDphi2D();
	I1 = ADphi*I1;
	//UO() << "I1=" << I1 << "\n";
	m(1+off,3+off) = I1(1);
	m(2+off,3+off) = I1(2);
	m(3+off,1+off) = I1(1);
	m(3+off,2+off) = I1(2);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_Rf: = A*Sbar //very important term

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();

	for (int j = 1; j <= FlexDOF(); j++)
	{
		Vector2D v(0,SbarS(j)); //stimmt mit rho*GetH() überein!!!
		v = A*v;
		m(j,1+off) = v(1);
		m(j,2+off) = v(2);
		m(1+off,j) = v(1);
		m(2+off,j) = v(2);
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_theta_rr: I_11+I_22
	m(3+off,3+off) = I1122S; //(GetIkl(1,1)+GetIkl(2,2));

	//m_theta_theta_ff: q_f^T*m_ff*q_f
	double mthth = 0;
	for (int i=1; i <= FlexDOF(); i++) xg(i) = XG(i);

	for (int i=1; i <= FlexDOF(); i++)
	{
		for (int j=1; j <= FlexDOF(); j++)
		{
			mthth += xg(i)*massmatrix(i,j)*xg(j);
		}
	}
	m(3+off,3+off) += mthth;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_ff
	for (int i=1; i <= FlexDOF(); i++)
	{
		for (int j=1; j <= FlexDOF(); j++)
		{
			m(i,j) = massmatrix(i,j);
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_f:
	temp.SetLen(FlexDOF());
	temp = Ibar12S;

	for (int i=1; i <= FlexDOF(); i++) 
	{
		m(i,off+3) = temp(i);
		m(off+3,i) = temp(i);
	}

	//UO() << "Massmatrix = " << m << "\n";
}


void Beam2DFFRF::EvalMff(Matrix& m, double t) 
{
	int ns = FlexDOF();
	if (massmatrix.Getcols() == ns)
	{
		m = massmatrix;
		return;
	}
	else
	{
		int ns = FlexDOF();

		m.SetSize(ns, ns);
		m.SetAll(0);

		Matrix H(1,ns);
		Vector temp;

		GetIntegrationRule(x1,w1,orderxyM);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			GetS0(temp, x1(i1));
			H.SetRowVec(temp, 1);

			double fact = w1(i1) * rhoA * beaml * 0.5;
			m += fact * (H.GetTp() * H);
		}

		massmatrix = m;
		//UO() << "Mff=" << m << "\n";
	}
};


//insert quadratic velocity vector
void Beam2DFFRF::EvalF2(Vector& f, double t)
{
	Body2D::EvalF2(f,t);

	TMStartTimer(22);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//add linear stiffness terms (negative on right hand side!!!!)
	for (int i=1; i <= FlexDOF(); i++)
	{
		for (int j=1; j <= FlexDOF(); j++)
		{
			f(i) -= K(i,j)*XG(j);
		}
	}
	TMStopTimer(22);


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//quadratic velocity vector, Shabana p.229:

	static Vector temp; 
	static Vector Qvf;
	Qvf.SetLen(FlexDOF());
	Qvf.SetAll(0);
	double Qvtheta;
	Vector2D QvR(0.,0.);

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();
	Matrix3D Atheta = ReferenceFrame().GetRotMatrixDphi2D();

	double theta = ReferenceFrame().GetAngle2D();
	double thetap = ReferenceFrame().GetAngle2DP();


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_R:
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.); //temporary point2D variable
	for (int j = 1; j <= FlexDOF(); j++)
	{
		pt.Y() += SbarS(j)*XG(j);
	}

	I1 += pt;
	QvR = Sqr(thetap)*(A*I1);

	pt=Vector2D(0.,0.);
	for (int j = 1; j <= FlexDOF(); j++)
	{
		pt.Y() += SbarS(j)*XGP(j);
	}
	QvR += (-2.*thetap)*(Atheta*pt);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_theta:
	temp.SetLen(FlexDOF());
	temp.SetAll(0);

	double mthth = 0;
	xg.SetLen(FlexDOF());

	for (int i=1; i <= FlexDOF(); i++) xg(i) = XG(i);
	Mult(massmatrix,xg,temp);

	//Ibar_0 = Ibar_11+Ibar_22 = 0

	for (int i=1; i <= FlexDOF(); i++)
	{
		mthth += XGP(i)*temp(i);
	}
	Qvtheta = -2.*thetap * mthth; //does not influence up very much


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Qv_f: //temp is still m_ff*q_f+Ibar_0

	for (int i=1; i <= FlexDOF(); i++)
	{
		Qvf(i) = Sqr(thetap)*temp(i);
	}

	//Sbar_tilde = 0

	//if (!GetMBS()->IsJacobianComputation())
	f(FlexDOF()+1) += QvR(1);
	f(FlexDOF()+2) += QvR(2);
	f(FlexDOF()+3) += Qvtheta;

	//fill in velocity vector terms (positive on right hand side):
	for (int i=1; i <= FlexDOF(); i++)
	{
		f(i) += Qvf(i);
	}
	//UO() << "F2 = " << f << "\n";

	//mass damping ... ?
	if (this->GetMassDamping() != 0) UO() << "Error: Mass damping in FFRF not possible!!!\n";
}

//->only for volumeloads (gravity ...); u = r_0 + A (u_bar + S*q_f)
void Beam2DFFRF::GetIntDuDq(Matrix& dudq) 
{
	//UO() << "Not yet implemented\n";
	dudq.FillWithZeros();

	int off = FlexDOF();
	//double rho = GetRho(); //DR 2013-02-04 deleted rho from class element

	//same as mRR/rho
	dudq(1+off,1) = GetMass()/rho;
	dudq(2+off,2) = GetMass()/rho;
	//UO() << "vol=" <<GetMass()/rho << "\n";

	//same asmRtheta/rho = A_theta*[I1+Sbar*qf]/rho
	Vector2D I1;
	I1(1) = I1S(1);
	I1(2) = I1S(2);

	Vector2D pt(0.,0.);

	for (int j = 1; j <= FlexDOF(); j++)
	{
		pt.Y() += SbarS(j)*XG(j);
	}

	I1 += pt;
	Matrix3D ADphi = ReferenceFrame().GetRotMatrixDphi2D();
	//UO() << "Adphi=" << ADphi << "\n";
	I1 = ADphi*I1;
	dudq(3+off,1) = I1(1)/rho;
	dudq(3+off,2) = I1(2)/rho;

	//same as m_Rf = A*Sbar:

	Matrix3D A = ReferenceFrame().GetRotMatrix2D();

	for (int j = 1; j <= FlexDOF(); j++)
	{
		Vector2D v(0,SbarS(j)); //stimmt mit rho*GetH() überein!!!
		v = A*v;
		dudq(j,1) = v(1)/rho;
		dudq(j,2) = v(2)/rho;
	}
	//UO() << "dudq=" << dudq << "\n";

}

//compute d/dq( u = r_0 + A (u_bar + S*q_f) )
void Beam2DFFRF::GetdPosdqT(const Vector2D& ploc, Matrix& d)
{
	//UO() << "Not yet implemented\n";

	Vector2D prel = GetPos2Drel(ploc);
	//double phi = ReferenceFrame().GetAngle2D();
	//UO() << "ploc = " << ploc << "; prel = " << prel << "\n";

	//dprel/dq
	d.SetSize(SOS(),2);
	d.FillWithZeros();

	double p0 = (2./beaml)*ploc.X();
	if (!(p0 == -1 || p0 == 1))
	{
		GetS0(SV, p0);

		for (int i = 1; i <= FlexDOF(); i++)
		{
			Matrix3D A = ReferenceFrame().GetRotMatrix2D();
			Vector2D tmp = A*Vector2D(0.,-SV(i));
			d(i,1)=tmp.X();
			d(i,2)=tmp.Y();
		}
	}

	//d_pref/dq:
	d(FlexDOF()+1,1) = 1;
	d(FlexDOF()+2,2) = 1;

	//d_A/dq*prel:
	Vector2D rotdphi_prel = ReferenceFrame().GetRotMatrixDphi2D()*prel;
	//UO() << "GetdPos:rotdphi=" << rotdphi_prel << "\n";
	d(FlexDOF()+3,1) = rotdphi_prel.X();
	d(FlexDOF()+3,2) = rotdphi_prel.Y();

	//UO() << "dpdq=" << d << "\n";
}



void Beam2DFFRF::GetSbar(Vector& Sbar) 
{
	//compute int_V (rho*x_k * S_l) d_V of element:
	//int dim = 1; 

	//compute int_V (rho * S) d_V of element:
	int ns = FlexDOF();

	Sbar.SetLen(ns);
	Sbar.SetAll(0);

	GetIntegrationRule(x1,w1,orderxyH);

	for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	{
		GetS0(SV,x1(i1));
		for (int i=1; i<=ns; i++)
		{
			Sbar(i) += - rhoA * w1(i1) * beaml * 0.5 * SV(i); //minus because y = -w(x)
		}
	}
}


//Shabana p. 209-211, inertia terms, 2D
void Beam2DFFRF::GetI1(Vector& I1) 
{
	//I1 = int(rho * u0bar)dV
	I1.SetLen(2);
	I1(1) = 0; //0 if x=0, y=0 is center of mass
	I1(2) = 0;
};

double Beam2DFFRF::GetIkl(int k, int l)
{
	//compute 'mass moment of inertia' of element:
	double Ikl = 0;

	if (k == 1 && l == 1)
	{
		Ikl = rhoA * Cub(beaml) / 12.;
	}

	return Ikl;
}

void Beam2DFFRF::GetIbarkl(int k, int l, Vector& I1)
{
	int ns = FlexDOF();
	I1.SetLen(ns);
	I1.SetAll(0);

	if (k == 1 && l == 2)
	{
		GetIntegrationRule(x1,w1,orderxyH); //max quadratic
		SV.SetLen(ns);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			GetS0(SV, x1(i1));

			for (int i = 1; i <= ns; i++)
			{
				I1(i) += - rhoA * w1(i1) * beaml*0.5*x1(i1) * beaml*0.5*SV(i); //minus because y = -w(x)
			}
		}
	}

}


//**********Calculation of Positions********//

Vector2D Beam2DFFRF::GetPos2D(const Vector2D& p_loc) const
{
	Vector2D rpos = ReferenceFrame().GetRefPos2D();
	return rpos+ReferenceFrame().GetRotMatrix2D()*GetPos2Drel(p_loc);
}

//-1..+1 based!!!
Vector2D Beam2DFFRF::GetVel2D(const Vector2D& p_loc) const
{
	return ReferenceFrame().GetRefVel2D()+ReferenceFrame().GetRotMatrix2D()*GetVel2Drel(p_loc)+
		ReferenceFrame().GetRotMatrix2DP()*GetPos2Drel(p_loc);
}

//-1..+1 based!!!
Vector2D Beam2DFFRF::GetPos2DD(const Vector2D& p_loc) const
{
	Vector2D rpos = ReferenceFrame().GetRefPos2DD();
	return rpos + ReferenceFrame().GetRotMatrix2DD()*GetPos2DrelD(p_loc);
}

//-1..+1 based!!!
Vector2D Beam2DFFRF::GetVel2DD(const Vector2D& p_loc) const
{
	return ReferenceFrame().GetRefVel2DD()+ReferenceFrame().GetRotMatrix2DD()*GetVel2DrelD(p_loc)+
		ReferenceFrame().GetRotMatrix2DPD()*GetPos2DrelD(p_loc);
}

//for roational constraint
void Beam2DFFRF::GetdAngle2DdqT(const Vector2D& ploc, Matrix& d)
{
	d.SetSize(SOS(),1);

	double p = 2.*ploc.X()/beaml;
	//double wx = GetWx(ploc.X());
	double fact = 1.;//1./pow(1.+sqr(wx),0.5)-wx/pow(1.+sqr(wx),1.5);
	for (int j = 1; j <= FlexDOF(); j++)
	{
		d(j, 1) = -fact*GetDS0(p, j) * 2./beaml;
	}

	d(FlexDOF()+1, 1) = 0;
	d(FlexDOF()+2, 1) = 0;
	d(FlexDOF()+3, 1) = 1;
}

double Beam2DFFRF::GetAngle2D(const Vector2D& ploc) const
{
	double p = 2.*ploc.X()/beaml;


	if (p == -1) 
	{
		return XG(FlexDOF()+3) - XG(1) * 2./beaml;
	} 
	else if (p == 1) 
	{
		return XG(FlexDOF()+3) - XG(2) * 2./beaml;
	}

	double w_x = 0; 
	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_x -= GetDS0(p, j) * XG(j) * 2./beaml;
	}
	//double factwx = 1./sqrt(1+Sqr(w_x));
	//w_x*=factwx;

	return XG(FlexDOF()+3) + w_x;
}

double Beam2DFFRF::GetAngle2DP(const Vector2D& ploc) const
{
	double p = 2.*ploc.X()/beaml;
	double w_x = 0;
	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_x -= GetDS0(p, j) * XGP(j) * 2./beaml;
	}

	//double factwx = 1./sqrt(1+Sqr(w_x));
	//w_x*=factwx;

	return XGP(FlexDOF()+3) + w_x;
}

double Beam2DFFRF::GetAngle2DD(const Vector2D& ploc) const
{
	double p = 2.*ploc.X()/beaml;
	double w_x = 0;

	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_x -= GetDS0(p, j) * XGD(j) * 2./beaml;
	}

	return XGD(FlexDOF()+3) + w_x;
}


double Beam2DFFRF::GetW(double ploc)
{
	double p = 2.*ploc/beaml;
	double w = 0;
	for (int j = 1; j <= FlexDOF(); j++)
	{
		w -= GetS0(p, j) * XG(j) * 2./beaml;
	}
	return w;
}


double Beam2DFFRF::GetWx(double ploc)
{
	double p = 2.*ploc/beaml;
	double w_x = 0;
	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_x -= GetDS0(p, j) * XG(j) * 2./beaml;
	}
	return w_x;
}

double Beam2DFFRF::GetWxx(double ploc)
{
  double p = 2.*ploc/beaml;

	double w_xx = 0;
	Vector temp;

	GetDDSMatrix0(temp,ploc);
	temp *= 4./Sqr(beaml);

	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_xx +=  temp(j)* XG(j);
	}
	return w_xx;
}

double Beam2DFFRF::GetWxxP(double ploc)
{
  double p = 2.*ploc/beaml;

	double w_xx = 0;
	Vector temp;

	GetDDSMatrix0(temp,ploc);
	temp *= 4./Sqr(beaml);

	for (int j = 1; j <= FlexDOF(); j++)
	{
		w_xx +=  temp(j)* XGP(j);
	}
	return w_xx;
}

double Beam2DFFRF::GetPower()
{
	//UO() << " Power FFRF \n";
	int ns = FlexDOF();
	double Udot = 0;


	int n = 10;
	double w = 1./(double)n*beaml;
	//for (int i1 = 1; i1 <= x1.GetLen(); i1++)
	for (int i1 = 1; i1 <= n; i1++)
	{ 
		double x = beaml * (((double)i1-0.5)/(double)n-0.5);
		
		double wxx = GetWxx(x);
		double wxxp = GetWxxP(x);
		Udot += fabs(EI * w * wxx*wxxp); 
	}
	return Udot;
}