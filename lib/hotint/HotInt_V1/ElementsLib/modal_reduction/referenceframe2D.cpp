//#**************************************************************
//#
//# filename:             referenceframe2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						20. April  2006
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
#include "body2d.h"
#include "femathhelperfunctions.h"
#include "material.h"
#include "node.h"
#include "referenceframe2d.h"
#include "elementdataaccess.h"

void ReferenceFrame2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body2D::GetElementData(edc);

	//IVector FFRFelements; //elements connected to ReferenceFrame
	//SearchTree searchtree; //for optimized node fill
	//TArray<Node*> nodes;
	//int resortconstraint;  //activates resorting of the DOF of the reference frame into the constraint part
	//int isACRS;  //absolute coordinates reduced strain --> the frame rotation and translation is not taken into account in GetPos2D() etc.


	ElementData ed;
	SetElemDataIVector(edc, FFRFelements, "FFRF_elements");

	ed.SetBool(draw_frame, "Draw_Frame"); edc.Add(ed);
}

int ReferenceFrame2D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body2D::SetElementData(edc);

	GetElemDataIVector(mbs, edc, "FFRF_elements", FFRFelements, 0);

	GetElemDataBool(mbs, edc, "Draw_Frame", draw_frame, 0);

	return rv;
}



void ReferenceFrame2D::EvalM(Matrix& m, double t) 
{
	//UO() << "RF EvalM\n";
	//everything is done by the elements
}; 

void ReferenceFrame2D::EvalF2(Vector& f, double t) 
{
	Element::EvalF2(f,t);
	//UO() << "frame-f=" << f << "\n";

	//UO() << "RF EvalF2\n";
	//everything is done by the elements

}; 



int ReferenceFrame2D::AddNode(Node* n)
{
	n->SetAuxNode();

	double tol = 1e-8*GetMBS()->CharacteristicLength();

	//standard:
	/*
	for (int i = 1; i <=nodes.Length(); i++)
	{
		if (Dist(n->Pos(),nodes(i)->Pos()) <= tol) return i;
	}
	Node* nc = n->GetCopy();
	return nodes.Add(nc);
	*/

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
	nc->NodeNum() = index; //store node number in reference configuration (local == global) for const int& access of virtual function nodenum in referenceframe2D

	searchtree.AddItem(Box3D(n->Pos(),n->Pos()), index);

	return index;
	
}

int ReferenceFrame2D::GetNode(Vector2D p_loc) const
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
	Vector3D p(p_loc.X(),p_loc.Y(),0.);

	Box3D box(p, p);
	box.Increase(tol);
	searchtree.GetItemsInBox(box,items);
	for (int i = 1; i <= items.Length(); i++)
	{
		Vector3D p2(nodes(items(i))->Pos().X(), nodes(items(i))->Pos().Y(), 0.);
		if (Dist(p,p2) <= tol) return items(i);
	}
	GetMBS()->UO() << "ERROR: node not found:" << p_loc << " !!!\n";
	return 0;
}

void ReferenceFrame2D::DrawElement() 
{
	if (draw_frame)
	{
		mbs->SetColor(col);

		//UO() << "XGD(1)=" << XGD(1) << ", XGD(2)=" << XGD(2) << "\n";

		Vector3D p1(GetPosD(Vector3D( 0.5*size.X(), 0.5*size.Y(),0.)));
		Vector3D p2(GetPosD(Vector3D(-0.5*size.X(), 0.5*size.Y(),0.)));
		Vector3D p3(GetPosD(Vector3D(-0.5*size.X(),-0.5*size.Y(),0.)));
		Vector3D p4(GetPosD(Vector3D( 0.5*size.X(),-0.5*size.Y(),0.)));

		/* //slidercrank:
		double x1 = 0.173205;
		double y1 = -0.1;
		double x2 = 0.926795;
		double y2 = 0.1;

		Vector3D p1(GetPosD(Vector3D(x2, y2,0.)));
		Vector3D p2(GetPosD(Vector3D(x1, y2,0.)));
		Vector3D p3(GetPosD(Vector3D(x1, y1,0.)));
		Vector3D p4(GetPosD(Vector3D(x2, y1,0.)));
		*/
		double th = 1;
		mbs->MyDrawLine(p1,p2,th);
		mbs->MyDrawLine(p2,p3,th);
		mbs->MyDrawLine(p3,p4,th);
		mbs->MyDrawLine(p4,p1,th);
	}
	
};

