//#**************************************************************
//#
//# filename:             referenceframe3D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						10. Juli  2007
//# description:          Reference frame element for FFRF formulation
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
 
#include "element.h"
#include "body3d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "node.h"
#include "rigid3d.h"
#include "rigid3dkardan.h"
#include "femeshinterface.h"
#include "rendercontext.h"
#include "referenceframe3d.h"
#include "elementdataaccess.h"


		// Set-Function,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
template <class RIGID>
void ReferenceFrame3D<RIGID>::SetReferenceFrame(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip,
																								const Vector3D& sizei, const Vector3D& coli)
{
	ComputeInitialConditions(p, v, phi, phip, x_init);

	//////// HACK AP
	//double xinit_double[14] = {0.9771035473772125, 0.8119080368041125, 0.8149451665816472, 0.4345721941229590, 0.4196667118043127, 0.4654793137759562, 0.6468043502831272, 0.1137455299827921, 0.8447456427769386, 0.7828729837346822, -0.05939125082671293, -0.4017941539650536, 0.3501751878606030, 0.04859291907932162};
	//for (int i=1; i<=14; i++)
	//	x_init(i) = xinit_double[i-1];

	col = coli;
	size = sizei;

	InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
	nodes.SetLen(0);
	isACRS = 0;
}


	// SetFunction for Rigid3DKardan,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
// nimodesi .. number of internal modes
// sizei ..... size of frame (for drawing)
// coli ...... color (for drawing, currently not used)
// rs ........ Kardan angle rotation sequence
template<>
void ReferenceFrame3D<Rigid3DKardan>::SetReferenceFrameKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip,
																															const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs)
{
	rotationsSequence = rs;
	ComputeInitialConditions(p, v, phi, phip, x_init);

	col = coli;
	size = sizei;

	InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
	nodes.SetLen(0);
	isACRS = 0;
}

template<>
void ReferenceFrame3D<Rigid3D>::SetReferenceFrameKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip,
																															const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs)
{
	GetMBS()->UO().InstantMessageText("SetReferenceFrameKardan called for ReferenceFrame3D<Rigid3D>!\n");
}

template <class RIGID>
void ReferenceFrame3D<RIGID>::AddFEMesh(const FEMeshInterface& femesh)
{
		// add nodes from FEMeshInterface
		for (int n=1; n<=femesh.NPoints(); n++)
		{
//			Node node(3, 1, femesh.Point(n)); //$ AD 2011-02-24: capsuled Nodes in FEMesh, Point(i) function split in Get & Set
			Node node(3, 1, femesh.GetPoint3D(n));
			int nn = mbs->AddNode(&node);
			AddNode(&mbs->GetNode(nn));
		}
		GetMBS()->UO() << "number of nodes " << femesh.NPoints() << "\n";

		for (int el=1; el<=femesh.NElements(); el++)
		{
			//JG, AD: use function from FEMesh to convert FEElements to MBSElements

			//*YV{
			femesh.Settings().domain = 1;
			femesh.Settings().generate_ffrf_elements = true;
			femesh.Settings().cms_element_number = GetOwnNum();
			Element * element_ptr = femesh.GetPtrToMBSElement_new(el); //get new element, force domain=1, 
			//*YV}
			if (element_ptr)
			{
				// reference-frame: add element to mbs
				int elnr = GetMBS()->AddElement(element_ptr);
				AddFFRFElement(elnr);

				delete element_ptr; //delete element, has been generated with new in "FEMesh::GetPtrToMBSElement_new"!
			}
		}
}


template <class RIGID>
void ReferenceFrame3D<RIGID>::EvalM(Matrix& m, double t) 
{
	m.FillWithZeros();
	/*
	UO() << "LTG=";
	for (int i=1; i <= 2*SOS(); i++)
	{
		UO() << LTG(i) << ",";
	}
	UO() << "\n";
	*/

	//RIGID::EvalM(m,t);
	//UO() << "RF EvalM\n";
	//everything is done by the elements
}; 

template <class RIGID>
void ReferenceFrame3D<RIGID>::EvalF2(Vector& f, double t) 
{
	//RIGID::EvalF2(f,t);
	
	Element::EvalF2(f,t); //not equivalent to RIGID !!!

	AddEPCqTterms(f);
	
}; 



template <class RIGID>
int ReferenceFrame3D<RIGID>::AddNode(Node* n)
{
	double tol = 1e-8*GetMBS()->CharacteristicLength();

	//optimized:

	IVector items;

	Box3D box(n->Pos(),n->Pos());
	box.Increase(tol); //disabled in example for witteven adams comparison Milano conference 2007

	searchtree.GetItemsInBox(box,items);
	for (int i = 1; i <= items.Length(); i++)
	{
		//if (Dist(n->Pos(),nodes(items(i))->Pos()) <= tol) return items(i);
		if (Dist(n->Pos(),nodes(items(i))->Pos()) <= tol && (n->GetBodyInd() == nodes(items(i))->GetBodyInd()) ) return items(i); //$ AD FENodes have domain
	}
	Node* nc = n->GetCopy();
	int index = nodes.Add(nc);
	searchtree.AddItem(Box3D(n->Pos(),n->Pos()), index);

	nc->NodeNum() = index;

	return index;

}

template <class RIGID>
int ReferenceFrame3D<RIGID>::GetNode(Vector3D p_loc) const
{
	double tol = 1e-8*GetMBS()->CharacteristicLength();
/*
	for (int i = 1; i <=nodes.Length(); i++)
	{
		if (Dist(p_loc,nodes(i)->Pos2D()) <= tol) return i;
	}
	GetMBS()->UO() << "ERROR: node not found:" << p_loc << " !!!\n";
	return 0;
*/

	IVector items;
	Vector3D p(p_loc);

	Box3D box(p, p);
	box.Increase(tol);
	searchtree.GetItemsInBox(box,items);
	for (int i = 1; i <= items.Length(); i++)
	{
		Vector3D p2(nodes(items(i))->Pos());
		if (Dist(p,p2) <= tol) return items(i);
	}
	GetMBS()->UO() << "ERROR: node not found:" << p_loc << " !!!\n";
	return 0;
}

template <class RIGID>
void ReferenceFrame3D<RIGID>::DrawElement() 
{

	mbs->SetColor(col);


	double lx = 0.5*size.X(); double ly = 0.5*size.Y(); double lz = 0.5*size.Z();
	//Vector3D p8(GetPosD(Vector3D(-lx,-ly,-lz)));
	//Vector3D p7(GetPosD(Vector3D(-lx,-ly, lz)));
	//Vector3D p6(GetPosD(Vector3D( lx,-ly,-lz)));
	//Vector3D p5(GetPosD(Vector3D( lx,-ly, lz)));
	//Vector3D p4(GetPosD(Vector3D(-lx, ly,-lz)));
	//Vector3D p3(GetPosD(Vector3D(-lx, ly, lz)));
	//Vector3D p2(GetPosD(Vector3D( lx, ly,-lz)));
	//Vector3D p1(GetPosD(Vector3D( lx, ly, lz)));
	Vector3D p0(GetPosD(Vector3D( 0.,0.,0.)));
	Vector3D px(GetPosD(Vector3D( lx,0.,0.)));
	Vector3D py(GetPosD(Vector3D( 0.,ly,0.)));
	Vector3D pz(GetPosD(Vector3D( 0.,0.,lz)));

	if (mbs->GetIOption(133)) //draw body outline
	{
		double th = 0.03*lx; // 1+mbs->GetIOption(115); //rigid body line thickness
		double hs = 0.09*lx;
		Vector3D color_black(0.,0.,0.);
		mbs->MyDrawArrow(p0, px, color_black, th, hs, 8);
		mbs->MyDrawArrow(p0, py, color_black, th, hs, 8);
		mbs->MyDrawArrow(p0, pz, color_black, th, hs, 8);
		// Text:
		char textx[10], texty[10], textz[10];
		sprintf(textx,"x%d", GetOwnNum());
		Vector3D pxt(GetPosD(Vector3D( lx,0.,0.1*lz)));
		GetMBS()->GetRC()->PrintText3D((float)pxt.X(),(float)pxt.Y(),(float)pxt.Z(),textx);
		sprintf(texty,"y%d", GetOwnNum());
		Vector3D pyt(GetPosD(Vector3D( 0.,ly,0.1*lz)));
		GetMBS()->GetRC()->PrintText3D((float)pyt.X(),(float)pyt.Y(),(float)pyt.Z(),texty);
		sprintf(textz,"z%d", GetOwnNum());
		Vector3D pzt(GetPosD(Vector3D( 0.1*lx,0.1*ly,lz)));
		GetMBS()->GetRC()->PrintText3D((float)pzt.X(),(float)pzt.Y(),(float)pzt.Z(),textz);
		//mbs->MyDrawLine(p1,p3, th);
		//mbs->MyDrawLine(p3,p4, th);
		//mbs->MyDrawLine(p4,p2, th);
		//mbs->MyDrawLine(p2,p1, th);

		//mbs->MyDrawLine(p5,p7, th);
		//mbs->MyDrawLine(p7,p8, th);
		//mbs->MyDrawLine(p8,p6, th);
		//mbs->MyDrawLine(p6,p5, th);

		//mbs->MyDrawLine(p1,p5, th);
		//mbs->MyDrawLine(p3,p7, th);
		//mbs->MyDrawLine(p4,p8, th);
		//mbs->MyDrawLine(p2,p6, th);
	}

};

template <class RIGID>
void ReferenceFrame3D<RIGID>::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	RIGID::GetElementData(edc);

	ElementData ed;

	ed.SetInt(isACRS, "IsACRS"); edc.Add(ed);
}

template <class RIGID>
int ReferenceFrame3D<RIGID>::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = RIGID::SetElementData(edc);

	GetElemDataInt(mbs, edc, "IsACRS", isACRS, 1);
	return rv;
}


template class ReferenceFrame3D<Rigid3D>;
template class ReferenceFrame3D<Rigid3DKardan>;