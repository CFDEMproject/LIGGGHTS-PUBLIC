//#**************************************************************
//# filename:             BEBeamGetKappa1.h
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
 

double GetKappa1(double dR1,double dR2,double dR3,double ddR1,double ddR2,double ddR3,double Ang,double dAng,double D1,double D2,double D3,double dD1,double dD2,double dD3) const {
 return (Power(D3,2)*(ddR2*dR1*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - ddR1*dR2*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + dAng*(Power(dR1,2) + Power(dR2,2))*(Power(dR1,2) + Power(dR2,2) + Power(dR3,2))) + Power(D1,2)*(ddR3*dR1*dR2*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - ddR2*dR1*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + dAng*(Power(dR2,2) + Power(dR3,2))*(Power(dR1,2) + Power(dR2,2) + Power(dR3,2))) + D1*(D3*ddR2*Power(dR1,2)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - D2*ddR3*Power(dR1,2)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - D3*ddR1*dR1*dR2*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + D2*ddR3*Power(dR2,2)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + D2*ddR1*dR1*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - D2*ddR2*dR2*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + D3*ddR3*dR2*dR3*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - D3*ddR2*Power(dR3,2)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - 2*dAng*dR1*(D2*dR2 + D3*dR3)*(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) - dD3*dR2*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5) + dD2*dR3*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5)) + D2*(dD3*dR1*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5) - dD1*dR3*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5) + D2*(dAng*Power(dR2,2)*(Power(dR1,2) + Power(dR3,2)) + dAng*Power(Power(dR1,2) + Power(dR3,2),2) + dR2*(-(ddR3*dR1) + ddR1*dR3)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)))) + D3*(-(dD2*dR1*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5)) + dD1*dR2*Power(Power(dR1,2) + Power(dR2,2) + Power(dR3,2),1.5) + D2*(-2*dAng*Power(dR2,3)*dR3 - ddR1*Power(dR2,2)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + dR3*(-(ddR3*dR1) + ddR1*dR3)*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2)) + dR2*(-2*dAng*dR3*(Power(dR1,2) + Power(dR3,2)) + ddR2*dR1*Sqrt(Power(dR1,2) + Power(dR2,2) + Power(dR3,2))))))/((Power(dR1,2) + Power(dR2,2) + Power(dR3,2))*(-2*D1*D2*dR1*dR2 + Power(D3,2)*(Power(dR1,2) + Power(dR2,2)) - 2*D3*(D1*dR1 + D2*dR2)*dR3 + Power(D2,2)*(Power(dR1,2) + Power(dR3,2)) + Power(D1,2)*(Power(dR2,2) + Power(dR3,2))));}