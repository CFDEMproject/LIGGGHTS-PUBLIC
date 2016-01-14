//#**************************************************************
//#
//# filename:             FE3DHexTet.cpp
//#
//# author:               Gerstmayr Johannes, YV
//#
//# generated:						October 2010
//# description:          3D-HexahedralGeneric and 3D-TetrahedralGeneric finite elements
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
 
#include "finiteElement3D.h"
#include "Material.h"
#include "node.h"
#include "solversettings_auto.h"
#include "femathhelperfunctions.h"
#include "rigid3d.h"
#include "rigid3dkardan.h"
#include "femeshinterface.h"
#include "referenceframe3d.h"
#include "basecmselement.h"
#include "finiteelement3dffrf.h"
#include "fe3dhextet.h"

template class HexahedralGeneric<FiniteElement3D>;
template class HexahedralGeneric<FiniteElement3DFFRF>;
template class TetrahedralGeneric<FiniteElement3D>;
template class TetrahedralGeneric<FiniteElement3DFFRF>;
template class PrismGeneric<FiniteElement3D>;
template class PyramidGeneric<FiniteElement3D>;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW HEXAHEDRAL    NEW HEXAHEDRAL    NEW HEXAHEDRAL    NEW HEXAHEDRAL    NEW HEXAHEDRAL    NEW HEXAHEDRAL
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//nodal coordinates and initial velocities of points in 3D; 
//order == 1 --> linear, order == 2 --> quadratic, etc.
template<class FiniteElement3DBase>
HexahedralGeneric<FiniteElement3DBase>::HexahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
											 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelememti): 
FiniteElement3DBase(mbsi)
{
	SetHexahedral(bodyindi, xc, vc, rhoi, Emi, nui, coli, CMSelememti);
};

template<class FiniteElement3DBase>
HexahedralGeneric<FiniteElement3DBase>::HexahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
											 int materialnumi, const Vector3D& coli, int CMSelementi): 
FiniteElement3DBase(mbsi)
{
	SetHexahedral(bodyindi, xc, vc, materialnumi, coli, CMSelementi);
};

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
    											double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi)
{

	if (xc.Length() != 8 && xc.Length() != 20) {GetMBS()->UO() << "Warning: HexahedralGeneric, number of nodes in constructor is invalid!!\n";}

	TArray<int> nodelist;

	SetFFRFElement(CMSelementi);
	for (int i=1; i <= xc.Length(); i++)
	{
		Node n1(3, bodyindi, xc(i));
		nodelist.Add(AddBodyNode(n1)); // adds the node either to the reference frame or to the bode and checks if node exists
	}

	SetHexahedral(bodyindi, xc, vc, nodelist, rhoi, Emi, nui, coli, CMSelementi);
}

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
    											int material_num, const Vector3D& coli, int CMSelementi)
{
	if (xc.Length() != 8 && xc.Length() != 20) {GetMBS()->UO() << "Warning: HexahedralGeneric, number of nodes in constructor is invalid!!\n";}

	TArray<int> nodelist;
	SetFFRFElement(CMSelementi);

	for (int i=1; i <= xc.Length(); i++)
	{
		Node n1(3, bodyindi, xc(i));
		nodelist.Add(AddBodyNode(n1)); // adds the node either to the reference frame or to the bode and checks if node exists
	}
	//UO() << "hexnodes=" << nodelist << "\n";

	SetHexahedral(bodyindi, xc, vc, nodelist, material_num, coli, CMSelementi);
}

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, const TArray<int>& nodelist,
    											double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi)
{
	//initial positions and velocities must be given in nodes; values of older files are ignored!
	Material mat(GetMBS(), rhoi, Emi, nui, 0);
	int material_num = GetMBS()->AddMaterial(mat);

	SetHexahedral(bodyindi, xc, vc, nodelist, material_num, coli, CMSelementi);
};

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, const TArray<int>& nodelist,
    											int material_num, const Vector3D& coli, int CMSelementi)
{
	//initial positions and velocities must be given in nodes; values of older files are ignored!
	SetHexahedral(bodyindi, nodelist, material_num, coli, CMSelementi);
};

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::SetHexahedral(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	if (nodelist.Length() != 8 && nodelist.Length() != 20) {GetMBS()->UO() << "Warning: HexahedralGeneric, number of nodes in constructor is invalid!!\n";}

	FiniteElement3DBase::SetFiniteElement3D(bodyindi, nodelist, material_num, coli, CMSelementi);
};

template<class FiniteElement3DBase>
double HexahedralGeneric<FiniteElement3DBase>::GetS0(const Vector3D& ploc, int shape) const
{
	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();

	//order of nodes:
	//           7      8
	//           +------+
	//          /|     /|
	// Z     5 +------+6|     
	// ^       | |    | |
	// | Y     |3+- - |-+ 4
	// |/      |/     |/
	// --->X 1 +------+ 2     
	//

	if (NNodes() == 8)
	{
		switch(shape)
		{
		case 1: return (1.0-r)*(1.0-s)*(1.0-t)/8.0; break;
		case 2: return (1.0+r)*(1.0-s)*(1.0-t)/8.0; break;
		case 3: return (1.0-r)*(1.0+s)*(1.0-t)/8.0; break;
		case 4: return (1.0+r)*(1.0+s)*(1.0-t)/8.0; break;
		case 5: return (1.0-r)*(1.0-s)*(1.0+t)/8.0; break;
		case 6: return (1.0+r)*(1.0-s)*(1.0+t)/8.0;	break;
		case 7: return (1.0-r)*(1.0+s)*(1.0+t)/8.0; break;
		case 8: return (1.0+r)*(1.0+s)*(1.0+t)/8.0;	break;
		}
	}
	else
	{
		//with case switch:
		switch(shape)
		{
		case 1: return (1.0-r)*(1.0-s)*(1.0-t)*(-r-s-t-2.0)/8.0; break;
		case 2: return (1.0+r)*(1.0-s)*(1.0-t)*(r-s-t-2.0)/8.0; break;
		case 3: return (1.0-r)*(1.0+s)*(1.0-t)*(-r+s-t-2.0)/8.0; break;
		case 4: return (1.0+r)*(1.0+s)*(1.0-t)*(r+s-t-2.0)/8.0; break;
		case 5: return (1.0-r)*(1.0-s)*(1.0+t)*(-r-s+t-2.0)/8.0; break;
		case 6: return (1.0+r)*(1.0-s)*(1.0+t)*(r-s+t-2.0)/8.0;	break;
		case 7: return (1.0-r)*(1.0+s)*(1.0+t)*(-r+s+t-2.0)/8.0; break;
		case 8: return (1.0+r)*(1.0+s)*(1.0+t)*(r+s+t-2.0)/8.0;	break;
		case 9: return (1.0-r*r)*(1.0-s)*(1.0-t)/4.0; break;
		case 10: return (1.0-r*r)*(1.0+s)*(1.0-t)/4.0; break;
		case 11: return (1.0-r*r)*(1.0-s)*(1.0+t)/4.0; break;
		case 12: return (1.0-r*r)*(1.0+s)*(1.0+t)/4.0; break;
		case 13: return (1.0-s*s)*(1.0-r)*(1.0-t)/4.0; break;
		case 14: return (1.0-s*s)*(1.0+r)*(1.0-t)/4.0; break;
		case 15: return (1.0-s*s)*(1.0-r)*(1.0+t)/4.0; break;
		case 16: return (1.0-s*s)*(1.0+r)*(1.0+t)/4.0; break;
		case 17: return (1.0-t*t)*(1.0-r)*(1.0-s)/4.0; break;
		case 18: return (1.0-t*t)*(1.0+r)*(1.0-s)/4.0; break;
		case 19: return (1.0-t*t)*(1.0-r)*(1.0+s)/4.0; break;
		case 20: return (1.0-t*t)*(1.0+r)*(1.0+s)/4.0; break;
		}
	}
	assert(0 && "Wrong number of nodes or wrong index.");
	return 0;
}

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const
{
	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();
	sm.SetSize(Dim(),NS());

	if (NNodes() == 8)
	{
		sm(1,1) = -(1.0-s)*(1.0-t)/8.0;
		sm(1,2) = (1.0-s)*(1.0-t)/8.0;
		sm(1,3) = -(1.0+s)*(1.0-t)/8.0;
		sm(1,4) = (1.0+s)*(1.0-t)/8.0;
		sm(1,5) = -(1.0-s)*(1.0+t)/8.0;
		sm(1,6) = (1.0-s)*(1.0+t)/8.0;
		sm(1,7) = -(1.0+s)*(1.0+t)/8.0;
		sm(1,8) = (1.0+s)*(1.0+t)/8.0;
		sm(2,1) = -(1.0-r)*(1.0-t)/8.0;
		sm(2,2) = -(1.0+r)*(1.0-t)/8.0;
		sm(2,3) = (1.0-r)*(1.0-t)/8.0;
		sm(2,4) = (1.0+r)*(1.0-t)/8.0;
		sm(2,5) = -(1.0-r)*(1.0+t)/8.0;
		sm(2,6) = -(1.0+r)*(1.0+t)/8.0;
		sm(2,7) = (1.0-r)*(1.0+t)/8.0;
		sm(2,8) = (1.0+r)*(1.0+t)/8.0;
		sm(3,1) = -(1.0-r)*(1.0-s)/8.0;
		sm(3,2) = -(1.0+r)*(1.0-s)/8.0;
		sm(3,3) = -(1.0-r)*(1.0+s)/8.0;
		sm(3,4) = -(1.0+r)*(1.0+s)/8.0;
		sm(3,5) = (1.0-r)*(1.0-s)/8.0;
		sm(3,6) = (1.0+r)*(1.0-s)/8.0;
		sm(3,7) = (1.0-r)*(1.0+s)/8.0;
		sm(3,8) = (1.0+r)*(1.0+s)/8.0;
	}
	else
	{
		sm(1,1) = -(1.0-s)*(1.0-t)*(-r-s-t-2.0)/8.0-(1.0-r)*(1.0-s)*(1.0-t)/8.0;
		sm(1,2) = (1.0-s)*(1.0-t)*(r-s-t-2.0)/8.0+(1.0+r)*(1.0-s)*(1.0-t)/8.0;
		sm(1,3) = -(1.0+s)*(1.0-t)*(-r+s-t-2.0)/8.0-(1.0-r)*(1.0+s)*(1.0-t)/8.0;
		sm(1,4) = (1.0+s)*(1.0-t)*(r+s-t-2.0)/8.0+(1.0+r)*(1.0+s)*(1.0-t)/8.0;
		sm(1,5) = -(1.0-s)*(1.0+t)*(-r-s+t-2.0)/8.0-(1.0-r)*(1.0-s)*(1.0+t)/8.0;
		sm(1,6) = (1.0-s)*(1.0+t)*(r-s+t-2.0)/8.0+(1.0+r)*(1.0-s)*(1.0+t)/8.0;
		sm(1,7) = -(1.0+s)*(1.0+t)*(-r+s+t-2.0)/8.0-(1.0-r)*(1.0+s)*(1.0+t)/8.0;
		sm(1,8) = (1.0+s)*(1.0+t)*(r+s+t-2.0)/8.0+(1.0+r)*(1.0+s)*(1.0+t)/8.0;
		sm(1,9) = -r*(1.0-s)*(1.0-t)/2.0;
		sm(1,10) = -r*(1.0+s)*(1.0-t)/2.0;
		sm(1,11) = -r*(1.0-s)*(1.0+t)/2.0;
		sm(1,12) = -r*(1.0+s)*(1.0+t)/2.0;
		sm(1,13) = -(1.0-s*s)*(1.0-t)/4.0;
		sm(1,14) = (1.0-s*s)*(1.0-t)/4.0;
		sm(1,15) = -(1.0-s*s)*(1.0+t)/4.0;
		sm(1,16) = (1.0-s*s)*(1.0+t)/4.0;
		sm(1,17) = -(1.0-t*t)*(1.0-s)/4.0;
		sm(1,18) = (1.0-t*t)*(1.0-s)/4.0;
		sm(1,19) = -(1.0-t*t)*(1.0+s)/4.0;
		sm(1,20) = (1.0-t*t)*(1.0+s)/4.0;
		sm(2,1) = -(1.0-r)*(1.0-t)*(-r-s-t-2.0)/8.0-(1.0-r)*(1.0-s)*(1.0-t)/8.0;
		sm(2,2) = -(1.0+r)*(1.0-t)*(r-s-t-2.0)/8.0-(1.0+r)*(1.0-s)*(1.0-t)/8.0;
		sm(2,3) = (1.0-r)*(1.0-t)*(-r+s-t-2.0)/8.0+(1.0-r)*(1.0+s)*(1.0-t)/8.0;
		sm(2,4) = (1.0+r)*(1.0-t)*(r+s-t-2.0)/8.0+(1.0+r)*(1.0+s)*(1.0-t)/8.0;
		sm(2,5) = -(1.0-r)*(1.0+t)*(-r-s+t-2.0)/8.0-(1.0-r)*(1.0-s)*(1.0+t)/8.0;
		sm(2,6) = -(1.0+r)*(1.0+t)*(r-s+t-2.0)/8.0-(1.0+r)*(1.0-s)*(1.0+t)/8.0;
		sm(2,7) = (1.0-r)*(1.0+t)*(-r+s+t-2.0)/8.0+(1.0-r)*(1.0+s)*(1.0+t)/8.0;
		sm(2,8) = (1.0+r)*(1.0+t)*(r+s+t-2.0)/8.0+(1.0+r)*(1.0+s)*(1.0+t)/8.0;
		sm(2,9) = -(1.0-r*r)*(1.0-t)/4.0;
		sm(2,10) = (1.0-r*r)*(1.0-t)/4.0;
		sm(2,11) = -(1.0-r*r)*(1.0+t)/4.0;
		sm(2,12) = (1.0-r*r)*(1.0+t)/4.0;
		sm(2,13) = -s*(1.0-r)*(1.0-t)/2.0;
		sm(2,14) = -s*(1.0+r)*(1.0-t)/2.0;
		sm(2,15) = -s*(1.0-r)*(1.0+t)/2.0;
		sm(2,16) = -s*(1.0+r)*(1.0+t)/2.0;
		sm(2,17) = -(1.0-t*t)*(1.0-r)/4.0;
		sm(2,18) = -(1.0-t*t)*(1.0+r)/4.0;
		sm(2,19) = (1.0-t*t)*(1.0-r)/4.0;
		sm(2,20) = (1.0-t*t)*(1.0+r)/4.0;
		sm(3,1) = -(1.0-r)*(1.0-s)*(-r-s-t-2.0)/8.0-(1.0-r)*(1.0-s)*(1.0-t)/8.0;
		sm(3,2) = -(1.0+r)*(1.0-s)*(r-s-t-2.0)/8.0-(1.0+r)*(1.0-s)*(1.0-t)/8.0;
		sm(3,3) = -(1.0-r)*(1.0+s)*(-r+s-t-2.0)/8.0-(1.0-r)*(1.0+s)*(1.0-t)/8.0;
		sm(3,4) = -(1.0+r)*(1.0+s)*(r+s-t-2.0)/8.0-(1.0+r)*(1.0+s)*(1.0-t)/8.0;
		sm(3,5) = (1.0-r)*(1.0-s)*(-r-s+t-2.0)/8.0+(1.0-r)*(1.0-s)*(1.0+t)/8.0;
		sm(3,6) = (1.0+r)*(1.0-s)*(r-s+t-2.0)/8.0+(1.0+r)*(1.0-s)*(1.0+t)/8.0;
		sm(3,7) = (1.0-r)*(1.0+s)*(-r+s+t-2.0)/8.0+(1.0-r)*(1.0+s)*(1.0+t)/8.0;
		sm(3,8) = (1.0+r)*(1.0+s)*(r+s+t-2.0)/8.0+(1.0+r)*(1.0+s)*(1.0+t)/8.0;
		sm(3,9) = -(1.0-r*r)*(1.0-s)/4.0;
		sm(3,10) = -(1.0-r*r)*(1.0+s)/4.0;
		sm(3,11) = (1.0-r*r)*(1.0-s)/4.0;
		sm(3,12) = (1.0-r*r)*(1.0+s)/4.0;
		sm(3,13) = -(1.0-s*s)*(1.0-r)/4.0;
		sm(3,14) = -(1.0-s*s)*(1.0+r)/4.0;
		sm(3,15) = (1.0-s*s)*(1.0-r)/4.0;
		sm(3,16) = (1.0-s*s)*(1.0+r)/4.0;
		sm(3,17) = -t*(1.0-r)*(1.0-s)/2.0;
		sm(3,18) = -t*(1.0+r)*(1.0-s)/2.0;
		sm(3,19) = -t*(1.0-r)*(1.0+s)/2.0;
		sm(3,20) = -t*(1.0+r)*(1.0+s)/2.0;
	}

}

// optimized by YV
template<class FiniteElement3DBase>
const FEFace & HexahedralGeneric<FiniteElement3DBase>::GetLocalFace(int i) const 
{
	static FEFace localFaces[] = {
		FEFace(1,5,7,3),
		FEFace(2,4,8,6),
		FEFace(1,2,6,5),
		FEFace(3,7,8,4),
		FEFace(1,3,4,2),
		FEFace(5,6,8,7),
		FEFace(0,0,0,0)
	};
	if(i < 1 || i > 6)
		i = 7;
	return localFaces[i-1];
}

//compute local face coordinates for drawing or for accessing face surface
//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const
{
	refpos = Vector3D(0.,0.,0.);
	v1 = Vector3D(1.,0.,0.);
	v2 = Vector3D(0.,2.,0.);
	v3 = Vector3D(0.);
	//example for hexahedral:
	switch(face)
	{
	case 1: 
		{ //left
			refpos = Vector3D(-1., 1.,-1.);
			v1 = Vector3D(0,-2.,0);
			v2 = Vector3D(0,0,2.);
			break;
		}
	case 2:
		{ //right
			refpos = Vector3D( 1.,-1.,-1.);
			v1 = Vector3D(0,2.,0);
			v2 = Vector3D(0,0,2.);
			break;
		}
	case 3:
		{ //bottom
			refpos = Vector3D(-1.,-1.,-1.);
			v1 = Vector3D(2.,0,0);
			v2 = Vector3D(0,0,2.);
			break;
		}
	case 4:
		{ //top
			refpos = Vector3D(-1., 1., 1.);
			v1 = Vector3D(2.,0,0);
			v2 = Vector3D(0,0,-2.);
			break;
		}
	case 5:
		{ //back
			refpos = Vector3D(-1., 1.,-1.);
			v1 = Vector3D(2.,0,0);
			v2 = Vector3D(0,-2.,0);
			break;
		}
	case 6:
		{ //front
			refpos = Vector3D(-1.,-1., 1.);
			v1 = Vector3D(2.,0,0);
			v2 = Vector3D(0,2.,0);
			break;
		}
	default: break;
	}
}

template<class FiniteElement3DBase>
Vector3D HexahedralGeneric<FiniteElement3DBase>::GetNodeLocPos(int i) const
{
	//generate local node position hierarchically:

	if (i <= 8)
	{
		Vector3D lp(-1.,-1.,-1.);
		int j = i-1;

		if (j&1) lp.X() += 2;
		if (j&2) lp.Y() += 2;
		if (j&4) lp.Z() += 2;
		return lp;
	}
	else
	{
		if (i <= 12) //X=0-midnode
		{
			Vector3D lp(0.,-1.,-1.);
			int j = i-9;
			if (j&1) lp.Y() += 2;
			if (j&2) lp.Z() += 2;
			return lp;
		}
		else if (i <= 16) //Y=0-midnode
		{
			Vector3D lp(-1.,0.,-1.);
			int j = i-13;
			if (j&1) lp.X() += 2;
			if (j&2) lp.Z() += 2;
			return lp;
		}
		else //Z=0-midnode
		{
			Vector3D lp(-1.,-1.,0.);
			int j = i-17;
			if (j&1) lp.X() += 2;
			if (j&2) lp.Y() += 2;
			return lp;
		}
	}
}

template<class FiniteElement3DBase>
int HexahedralGeneric<FiniteElement3DBase>::GetActualInterpolationOrder() const
{
	switch(NNodes())
	{
	case 8:		return 3;
	case 20:	return 4;
	}
	assert(0);
	return 0;
}

template<class FiniteElement3DBase>
void HexahedralGeneric<FiniteElement3DBase>::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Hexahedral);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
		ruleOrder = 2;
	else
	{
		if(integrationRule.settings.interpolationOrder == 3)
		{
			if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
				ruleOrder = 3;
			else
				ruleOrder = 2;
		}
		else
		{
			if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			{
				if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
					ruleOrder = 4;
				else
					ruleOrder = 3;
			}
			else
			{
				if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
					ruleOrder = 7;
				else
					ruleOrder = 4;
			}
		}
	}
	assert(ruleOrder != 0);

	IntegrationRule::DefineIntegrationRuleCube(integrationRule, ruleOrder);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW TET     NEW TET     NEW TET     NEW TET     NEW TET     NEW TET     NEW TET
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<class FiniteElement3DBase>
TetrahedralGeneric<FiniteElement3DBase>::TetrahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
											 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi): 
FiniteElement3DBase(mbsi)
{
	SetTetrahedral(bodyindi, xc, vc, rhoi, Emi, nui, coli);
};

template<class FiniteElement3DBase>
void TetrahedralGeneric<FiniteElement3DBase>::SetTetrahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
    											double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi)
{
	if (xc.Length() != 4 && xc.Length() != 10) {GetMBS()->UO() << "Warning: TetrahedralGeneric, number of nodes in constructor is invalid!!\n";}

	TArray<int> nodelist;

	SetFFRFElement(CMSelementi);
	for (int i=1; i <= xc.Length(); i++)
	{
		Node n1(3, bodyindi, xc(i));
		nodelist.Add(AddBodyNode(n1)); // adds the node either to the reference frame or to the bode and checks if node exists
	}

	Material mat(GetMBS(), rhoi, Emi, nui, 0);
	int material_num = GetMBS()->AddMaterial(mat);

	SetTetrahedral(bodyindi, nodelist, material_num, coli, CMSelementi);
}

template<class FiniteElement3DBase>
void TetrahedralGeneric<FiniteElement3DBase>::SetTetrahedral(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	if (nodelist.Length() != 4 && nodelist.Length() != 10) {GetMBS()->UO() << "Warning: TetrahedralGeneric, number of nodes in constructor is invalid!!\n";}

	FiniteElement3DBase::SetFiniteElement3D(bodyindi, nodelist, material_num, coli, CMSelementi);
};

template<class FiniteElement3DBase>
double TetrahedralGeneric<FiniteElement3DBase>::GetS0(const Vector3D& ploc, int shape) const
{
	//see also getnodelocnode in .h file
	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();

	//order of nodes:
	//                          
	//       t                  
	//      4+ _  9             
	//       |\  +_   s         
	//Z      | \    +3               
	//^    10+  +8 /|               
	//| Y    |  +\  +6               
	//|/     | /7 \ |                 
	//--->X 1+--+--+2 r              
	//          5

	if (NNodes() == 4)
	{
		switch(shape)
		{
		case 1: return 1. - r - s - t; break;
		case 2: return r; break;
		case 3: return s; break;
		case 4: return t; break;

		default: return 0;
		}
	}
	else //nnodes==10
	{

		switch(shape)
		{
		case 1: return 1.0-r-s-t-0.2E1*r*(1.0-r-s-t)-0.2E1*s*(1.0-r-s-t)-0.2E1*t*(1.0-r-s-t);
			break; case 2: return r-0.2E1*r*(1.0-r-s-t)-0.2E1*r*s-0.2E1*r*t;
			break; case 3: return s-0.2E1*r*s-0.2E1*s*(1.0-r-s-t)-0.2E1*s*t;
			break; case 4: return t-0.2E1*r*t-0.2E1*s*t-0.2E1*t*(1.0-r-s-t);
			break; case 5: return 4.0*r*(1.0-r-s-t);
			break; case 6: return 4.0*r*s;
			break; case 7: return 4.0*s*(1.0-r-s-t);
			break; case 8: return 4.0*r*t;
			break; case 9: return 4.0*s*t;
			break; case 10: return 4.0*t*(1.0-r-s-t);
			break; 
		default: return 0;
		}

	}

}

template<class FiniteElement3DBase>
void TetrahedralGeneric<FiniteElement3DBase>::GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const
{

	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();
	sm.SetSize(Dim(),NS());

	if (NNodes() == 4)
	{
		sm(1,1) = -1.0;
		sm(1,2) = 1.0;
		sm(1,3) = 0.0;
		sm(1,4) = 0.0;
		sm(2,1) = -1.0;
		sm(2,2) = 0.0;
		sm(2,3) = 1.0;
		sm(2,4) = 0.0;
		sm(3,1) = -1.0;
		sm(3,2) = 0.0;
		sm(3,3) = 0.0;
		sm(3,4) = 1.0;
	}
	else
	{
		sm(1,1) = -0.3E1+0.4E1*r+0.4E1*s+0.4E1*t;
		sm(1,2) = -0.1E1+0.4E1*r;
		sm(1,3) = 0.0;
		sm(1,4) = 0.0;
		sm(1,5) = 4.0-8.0*r-4.0*s-4.0*t;
		sm(1,6) = 4.0*s;
		sm(1,7) = -4.0*s;
		sm(1,8) = 4.0*t;
		sm(1,9) = 0.0;
		sm(1,10) = -4.0*t;
		sm(2,1) = -0.3E1+0.4E1*r+0.4E1*s+0.4E1*t;
		sm(2,2) = 0.0;
		sm(2,3) = -0.1E1+0.4E1*s;
		sm(2,4) = 0.0;
		sm(2,5) = -4.0*r;
		sm(2,6) = 4.0*r;
		sm(2,7) = 4.0-4.0*r-8.0*s-4.0*t;
		sm(2,8) = 0.0;
		sm(2,9) = 4.0*t;
		sm(2,10) = -4.0*t;
		sm(3,1) = -0.3E1+0.4E1*r+0.4E1*s+0.4E1*t;
		sm(3,2) = 0.0;
		sm(3,3) = 0.0;
		sm(3,4) = -0.1E1+0.4E1*t;
		sm(3,5) = -4.0*r;
		sm(3,6) = 0.0;
		sm(3,7) = -4.0*s;
		sm(3,8) = 4.0*r;
		sm(3,9) = 4.0*s;
		sm(3,10) = 4.0-4.0*r-4.0*s-8.0*t;

	}

}

// optimized by YV
template<class FiniteElement3DBase>
const FEFace & TetrahedralGeneric<FiniteElement3DBase>::GetLocalFace(int i) const 
{
	static FEFace localFaces[] = {
		FEFace(1,3,2),
		FEFace(1,2,4),
		FEFace(1,4,3),
		FEFace(2,3,4),
		FEFace(0,0,0)
	};
	if(i < 1 || i > 4)
		i = 5;
	return localFaces[i-1];
}

//compute local face coordinates for drawing or for accessing face surface
//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
template<class FiniteElement3DBase>
void TetrahedralGeneric<FiniteElement3DBase>::GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const
{
	refpos = Vector3D(0.,0.,0.);
	//before was: refpos = Vector3D(s2,s2,s2);
	v1 = Vector3D(1.,0.,0.);
	v2 = Vector3D(0.,2.,0.);
	v3 = Vector3D(0.);
	double s3 = 1; //other method for shrinking

	//tetrahedral:
	switch(face)
	{
	case 1:
		{ //bottom s-r
			v1 = Vector3D(0.,s3,0.);
			v2 = Vector3D(s3,0.,0.);
			break;
		}
	case 2:
		{ //side r-t
			v1 = Vector3D(s3,0.,0.);
			v2 = Vector3D(0.,0.,s3);
			break;
		}
	case 3:
		{ //side t-s
			v1 = Vector3D(0.,0.,s3);
			v2 = Vector3D(0.,s3,0.);
			break;
		}
	case 4:
		{ //side rs-t
			v1 = Vector3D(s3,0.,0.);
			v2 = Vector3D(0.,0.,s3);
			v3 = Vector3D(0.,s3,0.);
			break;
		}
	default: break;
	}
}

template<class FiniteElement3DBase>
Vector3D TetrahedralGeneric<FiniteElement3DBase>::GetNodeLocPos(int i) const
{
		//acc. to Bathe (new german version), Fig. 5.13, p.438
		//generate local node position hierarchically:
		switch(i)
		{
		//linear:
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 1., 0.); break;
		case 4: return Vector3D(0., 0., 1.); break;
		//quadratic:
		case 5: return Vector3D(0.5, 0. , 0. ); break;
		case 6: return Vector3D(0.5, 0.5, 0. ); break;
		case 7: return Vector3D(0. , 0.5, 0. ); break;
		case 8: return Vector3D(0.5, 0. , 0.5); break;
		case 9: return Vector3D(0. , 0.5, 0.5); break;
		case 10:return Vector3D(0. , 0. , 0.5); break;
		}

		return Vector3D(0.);
}

template<class FiniteElement3DBase>
int TetrahedralGeneric<FiniteElement3DBase>::GetActualInterpolationOrder() const
{
	switch(NNodes())
	{
	case 4:		return 1;
	case 10:	return 2;
	}
	assert(0);
	return 0;
}

template<class FiniteElement3DBase>
void TetrahedralGeneric<FiniteElement3DBase>::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Tetrahedral);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
		ruleOrder = 2;
	else
	{
		if(integrationRule.settings.interpolationOrder == 1)
		{
			if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
				ruleOrder = 3;
			else
				ruleOrder = 2;
		}
		else
		{
			if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			{
				if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
					ruleOrder = 4;
				else
					ruleOrder = 3;
			}
			else
			{
				if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Stiffness)
					ruleOrder = 5;
				else
					ruleOrder = 4;
			}
		}
	}
	assert(ruleOrder != 0);

	//$AP2012-05-25 bug fixed, correct element type chosen for integration rule
	//IntegrationRule::DefineIntegrationRuleCube(integrationRule, ruleOrder);
	IntegrationRule::DefineIntegrationRuleTetrahedron(integrationRule, ruleOrder);

}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW PRISM    NEW PRISM     NEW PRISM     NEW PRISM     NEW PRISM     NEW PRISM 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$EK 2013-03-04: created class PrismGeneric
template<class FiniteElement3DBase>
PrismGeneric<FiniteElement3DBase>::PrismGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
											 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi): 
FiniteElement3DBase(mbsi)
{
	SetPrism(bodyindi, xc, vc, rhoi, Emi, nui, coli);
};

template<class FiniteElement3DBase>
void PrismGeneric<FiniteElement3DBase>::SetPrism(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
    											double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi)
{
	if (xc.Length() != 6 && xc.Length() != 15) {GetMBS()->UO() << "Warning: PrismGeneric, number of nodes in constructor is invalid!!\n";}

	TArray<int> nodelist;

	//Floating frame reference finite element .... ????EK????
	SetFFRFElement(CMSelementi);
	for (int i=1; i <= xc.Length(); i++)
	{
		Node n1(3, bodyindi, xc(i));
		nodelist.Add(AddBodyNode(n1)); // adds the node either to the reference frame or to the bode and checks if node exists
	}

	Material mat(GetMBS(), rhoi, Emi, nui, 0);
	int material_num = GetMBS()->AddMaterial(mat);

	SetPrism(bodyindi, nodelist, material_num, coli, CMSelementi);
}

template<class FiniteElement3DBase>
void PrismGeneric<FiniteElement3DBase>::SetPrism(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	if (nodelist.Length() != 6 && nodelist.Length() != 15) {GetMBS()->UO() << "Warning: PrismGeneric, number of nodes in constructor is invalid!!\n";}

	FiniteElement3DBase::SetFiniteElement3D(bodyindi, nodelist, material_num, coli, CMSelementi);
};

template<class FiniteElement3DBase>
double PrismGeneric<FiniteElement3DBase>::GetS0(const Vector3D& ploc, int shape) const
{
	//see also getnodelocnode in .h file
	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();

  //order of nodes:
	//          6+
	//          /|\
	//       15+ | +14
	//       t/  |13\
	//      4+---+---+5
  //       |   +12 |
  //       |   |   |
  //       |   |   |
	//       |   |   |
	//     10+   |s  +11
	//Z      |  3+   |
	//^      |  / \  |
	//| Y    |9+  8+ |
	//|/     |/     \|
	//--->X 1+---+---+2 r
  //           7


		

	if (NNodes() == 6)
	{
		switch(shape)
		{
		case 1: return (1. - t)*(1. - r - s); break;
		case 2: return (1. - t)*r; break;
		case 3: return (1. - t)*s; break;
		case 4: return t*(1. - r - s); break;
		case 5: return t*r; break;
		case 6: return t*s; break;
		default: return 0;
		}
	}
	else //nnodes==15
	{
		switch(shape)
		{
		case 1: return (1 - r - s)*(1 - 2*r - 2*s)*(1-t)*(1-2*t); break;
		case 2: return (2*r-1)*r*(1-t)*(1-2*t); break;
		case 3: return (2*s-1)*s*(1-t)*(1-2*t); break;
		case 4: return (1-r-s)*(1-2*r-2*s)*t*(2*t-1); break;
		case 5: return (2*r-1)*r*t*(2*t-1); break;
		case 6: return (2*s-1)*s*t*(2*t-1); break;
		case 7: return 4 * r * (1 - r - s) * (1 - t)*(1-2*t); //node 7
		case 8: return 4 * r*s*(1 - t)*(1-2*t); //node 8
		case 9: return 4 * s * (1 - r - s) * (1 - t)*(1-2*t); //node 9
		case 10: return 4*(1-r-s)*t*(1-t); //node 10
	  case 11: return 4*r*t*(1-t); //node 11
		case 12: return 4*s*t*(1-t); //node 12
		//case 10: return 4*(1-r-s)*(1 - 2*r - 2*s)*t*(1-t); //node 10
		//case 11: return 4*r*(2*r-1)*t*(1-t); //node 11
		//case 12: return 4*s*(2*s-1)*t*(1-t); //node 12
		
		case 13: return 4 * (1-r-s)*r*t*(2*t-1); //node 13
		case 14: return 4 * r*s*t*(2*t-1); //node 14
		case 15: return 4 * (1-r-s)*s*t*(2*t-1); //node 15
		default: return 0;
		}
	}

}

template<class FiniteElement3DBase>
void PrismGeneric<FiniteElement3DBase>::GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const
{

	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();
	sm.SetSize(Dim(),NS());

	if (NNodes() == 6)
	{
		sm(1,1) = t - 1.;
		sm(1,2) = 1. - t;
		sm(1,3) = 0.;
		sm(1,4) = -t;
		sm(1,5) = t;
		sm(1,6) = 0.;
		sm(2,1) = t - 1.;
		sm(2,2) = 0.0;
		sm(2,3) = 1. - t;
		sm(2,4) = -t;
		sm(2,5) = 0.;
		sm(2,6) = t;
		sm(3,1) = r + s - 1.;
		sm(3,2) = -r;
		sm(3,3) = -s;
		sm(3,4) = 1. - r - s;
		sm(3,5) = r;
		sm(3,6) = s;
	}
	else
	{
		//sm(1,1) = (-3+4*r+4*s)*(1+t*(2*t-3)); //(1 - r - s)*(1 - 2*r - 2*s)*(1-t)*(1-2*t)  
		sm(1,1) =  -(1 - 2*r - 2*s)*(1-t)*(1-2*t) -2*(1 - r - s)*(1-t)*(1-2*t);
		sm(1,2) = (4*r-1)*(1-t)*(1-2*t); 
		sm(1,3) = 0.0;
		sm(1,4) = (-3+4*r+4*s)*t*(2*t-1);
		sm(1,5) = (4*r-1)*t*(2*t-1);
		sm(1,6) = 0;
		sm(1,7) = 4*(1-2*r-s)*(1-t)*(1-2*t);
		sm(1,8) = 4*s*(1-t)*(1-2*t);
		sm(1,9) = -4*s*(1-t)*(1-2*t);
		sm(1,10) = -4*t*(1-t);
		sm(1,11) = 4*t*(1-t);
		sm(1,12) = 0; 
		//sm(1,10) = 4*(-3+4*r+4*s)*t*(1-t);
		//sm(1,11) = 4*(4*r-1)*t*(1-t);
		//sm(1,12) = 0;

		sm(1,13) = 4*(1-2*r-s)*t*(2*t-1);
		sm(1,14) = 4*r*t*(2*t-1);
		sm(1,15) = -4*s*t*(2*t-1);

		//sm(2,1) = (-3+4*r+4*s)*(1+t*(2*t-3)); //(1 - r - s)*(1 - 2*r - 2*s)*(1-t)*(1-2*t)  
		sm(2,1) =  -(1 - 2*r - 2*s)*(1-t)*(1-2*t) -2*(1 - r - s)*(1-t)*(1-2*t);
		sm(2,2) = 0;
		sm(2,3) = (4*s-1)*(1-t)*(1-2*t);
		sm(2,4) = (-3+4*r+4*s)*t*(2*t-1);
		sm(2,5) = 0;
		sm(2,6) = (4*s-1)*t*(2*t-1);
		sm(2,7) = -4*r*(1-t)*(1-2*t);
		sm(2,8) = 4*r*(1-t)*(1-2*t);
		sm(2,9) = 4*(1-r-2*s)*(1-t)*(1-2*t);
		sm(2,10) = -4*t*(1-t);
		sm(2,11) = 0;
		sm(2,12) = 4*t*(1-t); 
    //sm(2,10) = 4*(-3+4*r+4*s)*t*(1-t);
		//sm(2,11) = 0;
		//sm(2,12) = 4*(4*s-1)*t*(1-t); 

		sm(2,13) = -4*r*t*(2*t-1);
		sm(2,14) =4*r*t*(2*t-1);
		sm(2,15) = 4*(1-r-2*s)*t*(2*t-1);

		sm(3,1) = (1-r-s)*(1-2*r-2*s)*(4*t-3);
		sm(3,2) = (2*r-1)*r*(4*t-3);
		sm(3,3) = (2*s-1)*s*(4*t-3);
		sm(3,4) = (1-r-s)*(1-2*r-2*s)*(4*t-1);
		sm(3,5) = (2*r-1)*r*(4*t-1);
		sm(3,6) = (2*s-1)*s*(4*t-1);
		sm(3,7) = 4*(1-r-s)*r*(4*t-3);
		sm(3,8) = 4*r*s*(4*t-3);
		sm(3,9) = 4*(1-r-s)*s*(4*t-3);
		sm(3,10) = 4*(1-r-s)*(1-2*t);
		sm(3,11) = 4*r*(1-2*t);
		sm(3,12) = 4*s*(1-2*t); 

		//sm(3,10) = 4*(1-r-s)*(1 - 2*r - 2*s)*(1-2*t);
		//sm(3,11) = 4*r*(2*r-1)*(1-2*t);
		//sm(3,12) = 4*s*(2*s-1)*(1-2*t); 

		sm(3,13) = 4*(1-r-s)*r*(4*t-1);
		sm(3,14) = 4*r*s*(4*t-1);	
		sm(3,15) = 4*(1-r-s)*s*(4*t-1);	
	}

}

// optimized by YV
template<class FiniteElement3DBase>
const FEFace & PrismGeneric<FiniteElement3DBase>::GetLocalFace(int i) const 
{
	static FEFace localFaces[] = {
		FEFace(1,2,5,4),
		FEFace(1,4,6,3),
		FEFace(2,3,6,5),
		FEFace(1,3,2),
		FEFace(4,6,5),
		FEFace(0,0,0)
	};
	if(i < 1 || i > 5)
		i = 6;
	return localFaces[i-1];
}

//compute local face coordinates for drawing or for accessing face surface
//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
template<class FiniteElement3DBase>
void PrismGeneric<FiniteElement3DBase>::GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const
{
	refpos = Vector3D(0.,0.,0.);
	//before was: refpos = Vector3D(s2,s2,s2);
	v1 = Vector3D(1.,0.,0.);
	v2 = Vector3D(0.,1.,0.);
	v3 = Vector3D(0.);
	double s3 = 1; //other method for shrinking

	//$EK  //tetrahedral:
	//prism
	switch(face)
	{
	case 1:
		{ //bottom r-t
			v1 = Vector3D(s3,0.,0.);
			v2 = Vector3D(0.,0.,s3);
			break;
		}
	case 2:
		{ //side s-t
			//v1 = Vector3D(0.,s3,0.);
			//v2 = Vector3D(0.,0.,s3);
			v1 = Vector3D(0, 0, 1.);
			v2 = Vector3D(0.,1.,0.);
			break;
		}
	case 3:
		{ //side rs-t
			refpos = Vector3D(1.,0.,0.);
			v1 = Vector3D(-s3,s3,0.);
			v2 = Vector3D(0.,0.,s3);
			//v3 = Vector3D(0.,s3,0.);
			break;
		}		
	case 4:
		{ //side r-s; bottom
			//v1 = Vector3D(s3,0.,0.);
			//v2 = Vector3D(0.,s3,0.);
			v1 = Vector3D(0, 1.,0.);
			v2 = Vector3D(1.,0.,0.);
			break;
		}
	case 5:
		{ //side r-s, top
			refpos = Vector3D(0.,0.,1.);
			//v1 = Vector3D(s3,0.,0.);
			//v2 = Vector3D(0.,s3,0.);
			v1 = Vector3D(0, 1.,0.);
			v2 = Vector3D(1.,0.,0.);
			break;
		}
	default: break;
	}
}



template<class FiniteElement3DBase>
Vector3D PrismGeneric<FiniteElement3DBase>::GetNodeLocPos(int i) const
{
		//acc. to Bathe (new german version), Fig. 5.13, p.438
		//generate local node position hierarchically:
		switch(i)
		{
		//linear:
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 1., 0.); break;
		case 4: return Vector3D(0., 0., 1.); break;
		case 5: return Vector3D(1., 0., 1.); break;
		case 6: return Vector3D(0., 1., 1.); break;
		//quadratic:
		case 7: return Vector3D(0.5, 0. , 0. ); break;
		case 8: return Vector3D(0.5, 0.5, 0. ); break;
		case 9: return Vector3D(0., 0.5, 0. ); break;
		case 10: return Vector3D(0., 0., 0.5); break;
		case 11: return Vector3D(1., 0., 0.5); break;
		case 12: return Vector3D(0., 1., 0.5); break;
		case 13: return Vector3D(0.5, 0., 1.); break;
		case 14: return Vector3D(0.5, 0.5, 1.); break;
		case 15: return Vector3D(0., 0.5, 1.); break;
		}
		return Vector3D(0.);
}

template<class FiniteElement3DBase>
int PrismGeneric<FiniteElement3DBase>::GetActualInterpolationOrder() const
{
	//$EK NEEDS TO BE CHECKED
	switch(NNodes())
	{
	case 6:		return 2;
	case 15:	return 4;
	}
	assert(0);
	return 0;
}

template<class FiniteElement3DBase>
void PrismGeneric<FiniteElement3DBase>::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Prism);

	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
	{
		//if (integrationRule.settings.interpolationOrder == 2)
		//	ruleOrder = 1;
		//else
		//	ruleOrder = 2;
		ruleOrder = integrationRule.settings.interpolationOrder;
	}
	else if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Mass)
	{
		ruleOrder = integrationRule.settings.interpolationOrder;
	}
	else //stiffness
	{
		if(integrationRule.settings.interpolationOrder == 1)
		{
			if (integrationRule.settings.geometricNonlinearityStatus == GNS_Linear)
			{	ruleOrder = 1; }
			else
			{				
				if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
				{	ruleOrder = 3; }
				else
				{	ruleOrder = 2; }
			}
		}
		else
		{
			if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			{	ruleOrder = 4; }
			else
			{	ruleOrder = 3; }
		}
	}
	assert(ruleOrder != 0);

	//$AP2012-05-25 bug fixed, correct element type chosen for integration rule
	//IntegrationRule::DefineIntegrationRuleCube(integrationRule, ruleOrder);
	IntegrationRule::DefineIntegrationRulePrism(integrationRule, ruleOrder);

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW PYRAMID  NEW PYRAMID   NEW PYRAMID   NEW PYRAMID   NEW PYRAMID   NEW PYRAMID 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$EK 2013-03-04: created class PyramidGeneric
template<class FiniteElement3DBase>
PyramidGeneric<FiniteElement3DBase>::PyramidGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
											 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi): 
FiniteElement3DBase(mbsi)
{
	SetPyramid(bodyindi, xc, vc, rhoi, Emi, nui, coli);
};

template<class FiniteElement3DBase>
void PyramidGeneric<FiniteElement3DBase>::SetPyramid(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
    											double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi)
{
	if (xc.Length() != 5 && xc.Length() != 13) {GetMBS()->UO() << "Warning: PyramidGeneric, number of nodes in constructor is invalid!!\n";}

	TArray<int> nodelist;

	//Floating frame reference finite element .... ????EK????
	SetFFRFElement(CMSelementi);
	for (int i=1; i <= xc.Length(); i++)
	{
		Node n1(3, bodyindi, xc(i));
		nodelist.Add(AddBodyNode(n1)); // adds the node either to the reference frame or to the bode and checks if node exists
	}

	Material mat(GetMBS(), rhoi, Emi, nui, 0);
	int material_num = GetMBS()->AddMaterial(mat);

	SetPyramid(bodyindi, nodelist, material_num, coli, CMSelementi);
}

template<class FiniteElement3DBase>
void PyramidGeneric<FiniteElement3DBase>::SetPyramid(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi)
{
	if (nodelist.Length() != 5 && nodelist.Length() != 13) {GetMBS()->UO() << "Warning: PyramidGeneric, number of nodes in constructor is invalid!!\n";}

	FiniteElement3DBase::SetFiniteElement3D(bodyindi, nodelist, material_num, coli, CMSelementi);
};

// $EK 2013-03-06 only first order elements are implemented, second order missing
template<class FiniteElement3DBase>
double PyramidGeneric<FiniteElement3DBase>::GetS0(const Vector3D& ploc, int shape) const
{
	//see also getnodelocnode in .h file
	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();

	 //order of nodes:
	//     5 +...
	//       |\\.\...    13
	//       | \ \.  \..+
	//       |12+  \     \...
	//       |   \  \        \..
	//       |    \s \    7     \ 
	//		 10+    3+--+---+------+4
	//       |    /  11\        /
	//Z      |   /      \      /
	//^      | 8+        \   9+
	//| Y    | /          \  /
	//|/     |/            \/
	//--->X 1+------+------+2  r
	//              6

	//according to NGSolve
	if (t == 1.) t -= 1e-10;
	if (NNodes() == 5 || NNodes() == 13) 
	{
		switch(shape)
		{
		// $EK possibility for basis functions - problem with Area Loads!!
		//case 1: return (1. - r)*(1- s)*(1-t); break;
		//case 2: return r*(1- s)*(1-t); break;
		//case 3: return (1. - r)*s*(1-t); break;
		//case 4: return r*s*(1-t); break;		
		//case 5: return t; break;
		// $EK basis functions of NGSolve
		case 1: return (1 - r -t)*(1 - s - t)/(1-t); break;
		case 2: return r*(1 - s - t)/(1-t); break;
		case 3: return s*(1 - r - t)/(1-t); break;
		case 4: return r*s/(1-t); break;
		case 5: return t; break;
		default: return 0;
		}
	}
	if (NNodes() == 13) //nnodes==13
	{
//-----------------------------------------------------------------------------
	//NGSolve: Higher order pyramids... problem -> is not a partition of unity, 
	//which is needed in HOTINT for drawing reasons ($EK)
		/*
		    // horizontal edge dofs 
    for (int i = 0; i < 4; i++)
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
	  INT<2> e = GetEdgeSort (i, vnums);	  

	  Tx xi = sigma[e[1]]-sigma[e[0]]; 
	  Tx lam_e = lambda[e[0]]+lambda[e[1]];
	  Tx bub = 0.25 * lam_e * (1 - xi*xi)*(1-z)*(1-z);
	  
	  LegendrePolynomial::
	    EvalScaledMult (p-2, xi*(1-z), 1-z, bub, shape.Addr(ii));
	  ii += p-1;
	}
    
    // vertical edges
    for (int i = 4; i < 8; i++) 
      if (order_edge[i] >= 2)
	{
	  int p = order_edge[i];
	  INT<2> e = GetEdgeSort (i, vnums);	  

	  Tx xi = lambda3d[e[1]]-lambda3d[e[0]]; 
	  Tx lam_e = lambda3d[e[0]]+lambda3d[e[1]];
	  Tx bub = 0.25 * (lam_e*lam_e-xi*xi);
	  
	  LegendrePolynomial::
	    EvalScaledMult (p-2, xi, lam_e, bub, shape.Addr(ii));
	  ii += p-1;
	}
		*/
//-----------------------------------------------------------------------------
		//$EK not implemented so far !!!!!!!!!!!!!!!!!!!!!!!!!!
		switch(shape)
		{
		case 5: return 0; break;
		case 6: return 0; break;
		case 7: return 0; break;
		case 8: return 0; break;
		case 9: return 0; break;
		case 10: return 0; break;
		case 11: return 0; break;
		case 12: return 0; break;
		case 13: return 0; break;
		}
		return 0;
	}

}

template<class FiniteElement3DBase>
void PyramidGeneric<FiniteElement3DBase>::GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const
{

	double r = ploc.X();
	double s = ploc.Y();
	double t = ploc.Z();
	sm.SetSize(Dim(),NS());
	if (t == 1.) t -= 1e-10;	
	if (NNodes() == 5)
	{
		/*sm(1,1) = -(1- s)*(1-t);
		sm(1,2) = (1- s)*(1-t);
		sm(1,3) = -s*(1-t);
		sm(1,4) = s*(1-t);
		sm(1,5) = 0;
		sm(2,1) = -(1. - r)*(1-t);
		sm(2,2) = -r*(1-t);
		sm(2,3) = (1. - r)*(1-t);
		sm(2,4) = r*(1-t);
		sm(2,5) = 0.;
		sm(3,1) = -(1. - r)*(1- s);
		sm(3,2) = -r*(1- s);
		sm(3,3) = -(1. - r)*s;
		sm(3,4) = -r*s;
		sm(3,5) = 1;*/
		
		// according to the basis functions of NGSolve
		sm(1,1) = -(1 - s - t)/(1-t);
		sm(1,2) = (1 - s - t)/(1-t);
		sm(1,3) = -s/(1-t);
		sm(1,4) = s/(1-t);
		sm(1,5) = 0;
		sm(2,1) = -(1 - r -t)/(1-t);
		sm(2,2) = -r/(1-t);
		sm(2,3) = (1 - r - t)/(1-t);
		sm(2,4) = r/(1-t);
		sm(2,5) = 0.;
		sm(3,1) = (  ( -(1 - r -t) - (1 - s - t))*(1-t) + (1 - r -t)*(1 - s - t))/sqr(1-t);
		sm(3,2) = (-r*(1-t) + r*(1 - s - t) )/sqr(1-t);
		sm(3,3) = (-s*(1-t) + s*(1 - r - t))/sqr(1-t);
		sm(3,4) = r*s/sqr(1-t);
		sm(3,5) = 1;

	}
	else
	{
		//$EK not implemented so far !!!!!!!!!!!!!!!!!!!!!!!!!!
		sm = 0.; 
		
	}

}

// optimized by YV
template<class FiniteElement3DBase>
const FEFace & PyramidGeneric<FiniteElement3DBase>::GetLocalFace(int i) const 
{
	static FEFace localFaces[] = {
		FEFace(1,2,4,3),
		FEFace(1,2,5),
		FEFace(2,4,5),
		FEFace(3,5,4),
		FEFace(1,5,3),
		FEFace(0,0,0)
	};
	if(i < 1 || i > 5)
		i = 6;
	return localFaces[i-1];
}

//compute local face coordinates for drawing or for accessing face surface
//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
template<class FiniteElement3DBase>
void PyramidGeneric<FiniteElement3DBase>::GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const
{
	refpos = Vector3D(0.,0.,0.);
	//before was: refpos = Vector3D(s2,s2,s2);
	v1 = Vector3D(1.,0.,0.);
	v2 = Vector3D(0.,1.,0.);
	v3 = Vector3D(0.);

	//$EK - pyramid:
	switch(face)
	{
	case 1:
		{ //bottom r-s
			v1 = Vector3D(1,0.,0.);
			v2 = Vector3D(0.,1,0);
			break;
		}
	case 2:
		{ //side r-t
			v1 = Vector3D(1,0.,0.);
			v2 = Vector3D(0.,0.,1.);
			break;
		}
	case 3:
		{ //side rs-t
			refpos = Vector3D(1.,0.,0.);
			v1 = Vector3D(0.,1.,0.);
			v2 = Vector3D(-1.,0.,1.);
			//v3 = Vector3D(0.,s3,0.);
			break;
		}		
	case 4:
		{ //side r-t; back
			refpos = Vector3D(0.,1.,0.);
			v1 = Vector3D(0.,-1.,1.);
			v2 = Vector3D(1.,0.,0.);
			break;
		}
	case 5:
		{ //side r-s, top
			v1 = Vector3D(0.,0.,1.);
			v2 = Vector3D(0.,1.,0.);
			break;
		}
	default: break;
	}
}



template<class FiniteElement3DBase>
Vector3D PyramidGeneric<FiniteElement3DBase>::GetNodeLocPos(int i) const
{
		//acc. to Bathe (new german version), Fig. 5.13, p.438
		//generate local node position hierarchically:
		switch(i)
		{
		//linear:
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 1., 0.); break;
		case 4: return Vector3D(1., 1., 0.); break;
		case 5: return Vector3D(0., 0., 1.); break;
		//quadratic:
		case 6: return Vector3D(0.5, 0., 0.); break;
		case 7: return Vector3D(0.5, 1. , 0. ); break;
		case 8: return Vector3D(0., 0.5, 0. ); break;
		case 9: return Vector3D(1., 0.5, 0. ); break;
		case 10: return Vector3D(0., 0., 0.5); break;
		case 11: return Vector3D(0.5, 0., 0.5); break;
		case 12: return Vector3D(0., 0.5, 0.5); break;
		case 13: return Vector3D(0.5, 0.5, 0.5); break;
		}
		return Vector3D(0.);
}

template<class FiniteElement3DBase>
int PyramidGeneric<FiniteElement3DBase>::GetActualInterpolationOrder() const
{
	//$EK NEEDS TO BE CHECKED
	switch(NNodes())
	{
	case 5:		return 3;
	case 13:	return 4;
	}
	assert(0);
	return 0;
}

template<class FiniteElement3DBase>
void PyramidGeneric<FiniteElement3DBase>::DefineIntegrationRule(IntegrationRule & integrationRule)
{
	assert(integrationRule.settings.elementType == TFE_Pyramid);
	int ruleOrder = 0;
	// here the particular rule order will be chosen depending on the settings
	if(integrationRule.settings.integratedValueType == IntegrationRule::IVT_Load)
	{		
		ruleOrder = integrationRule.settings.interpolationOrder;
	}
	else if (integrationRule.settings.integratedValueType == IntegrationRule::IVT_Mass)
	{
		ruleOrder = integrationRule.settings.interpolationOrder;
	}
	else //stiffness
	{
		if(integrationRule.settings.interpolationOrder == 3)
		{
			if (integrationRule.settings.geometricNonlinearityStatus == GNS_Linear)
			{	ruleOrder = 1; }
			else
			{				
				if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
				{	ruleOrder = 3; }
				else
				{	ruleOrder = 2; }
			}
		}
		else
		{
			if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearLargeStrain)
			{	ruleOrder = 4; }
			else if(integrationRule.settings.geometricNonlinearityStatus == GNS_NonlinearSmallStrain)
			{	ruleOrder = 3; }
			else
			{	ruleOrder = 2; }
		}
	}
	assert(ruleOrder != 0);

	//$EK 2013-03-04 integration rule for Pyramid
	//IntegrationRule::DefineIntegrationRuleCube(integrationRule, ruleOrder);
	IntegrationRule::DefineIntegrationRulePyramid(integrationRule, ruleOrder);
}
