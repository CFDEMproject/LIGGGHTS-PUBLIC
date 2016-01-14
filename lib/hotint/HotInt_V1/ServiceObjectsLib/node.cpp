//#**************************************************************
//#
//# filename:             node.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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


#include "element.h"
#include "node.h"
#include "elementdataaccess.h"
#include "rendercontext.h"

void Node::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	mystr sname = "";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ed.SetText("General", "Node_type"); ed.SetLocked(1); edc.Add(ed); //if different node types are used, insert here!
	ed.SetInt(nodenum, "Node_number"); ed.SetLocked(1); edc.Add(ed);
	ed.SetText(nodename.c_str(), "Node_name"); edc.Add(ed);

	ed.SetInt(sos, "Node_DOF_size"); ed.SetToolTipText("Set number of nodal coordinates"); edc.Add(ed);
	ed.SetInt(body, "Body_number"); ed.SetToolTipText("Set body number where node belongs to otherwise zero"); edc.Add(ed);

	ed.SetVector3D(pos.X(), pos.Y(), pos.Z(), "Node_position"); edc.Add(ed);
	SetElemDataIVector(edc, ltg, "Node_DOF_numbers"); edc.Last().SetLocked(1);
	edc.Last().SetToolTipText("Local to global DOF reference numbers, position and velocities");
}

int Node::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = 1;

	body = 0;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	GetElemDataText(0, edc, "Node_name", nodename, 0);

	GetElemDataInt(0, edc, "Node_DOF_size", sos, 1);
	GetElemDataInt(0, edc, "Body_number", body, 0);

	GetElemDataVector3D(0, edc, "Node_position", pos, 1);
	//GetElemDataIVector(0, edc, "Node_DOF_numbers", ltg, 1);
	//if (ltg.Length() != 2*sos) GetMBS()->UO() << "Error in node " << nodenum << ": Node DOF size and length of Node DOF numbers do not match!!!\n";

	return rv;
}

int Node::Dim() const
{
	if (sos <= 3) return sos;
	if (elements.Length() != 0)
	{
		int elem = elements(1);
		if (elem > 0 && elem <= GetMBS()->NE())
			return GetMBS()->GetElement(elem).Dim();
	}
	return 3; //default value
}

void Node::DrawNode()
{
	int res = GetMBS()->GetIOption(147);
	double size = GetMBS()->GetDOption(124);

	if (res == 0) return;

	if (size != 0)
	{
		mbs->SetColor(GetCol());
		mbs->DrawSphere(GetPosD(),0.5*size,res);
	}
	if (GetMBS()->GetIOption(148)) //show node numbers
	{
		Vector3D p = GetPosD();

		char text[40];
		sprintf(text,"%d", NodeNum());

		GetMBS()->GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
	}

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%% Generic Node %%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void GenericNode::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{

	GetElementDataAuto(edc);
}

int GenericNode::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = 	SetElementDataAuto(edc);

	return rv;
}


int GenericNode::AddNodeDOF(NodeDOFType ndt)
{
	// first add NodeDOFType to list
	int i = doftypes.Add(ndt);
	// than increase the respective value of sos, es, is
	// THIS is just an idea how it could look like
	if (ndt & NDT_sos) sos++;
	else if (ndt & NDT_es) es++;
	else if (ndt & NDT_is) is++;
	else mbs->UO() << "invalid NodeDOFType: DOF-order missing\n";
	
	return i;
}

int GenericNode::Dim() const
{
	int dim = 0;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_position || doftypes(i) & NDT_displacement)
			dim++;
	}
	return dim;
}

Vector3D GenericNode::GetDisplacement() const
{
	Vector3D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_displacement)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XG(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XG(i);
			else if (doftypes(i) & NDT_comp_3) p(3) = XG(i);
		}
		else if (doftypes(i) & NDT_position)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XG(i) - pos(1);
			else if (doftypes(i) & NDT_comp_2) p(2) = XG(i) - pos(2);
			else if (doftypes(i) & NDT_comp_3) p(3) = XG(i) - pos(3);
		}
	}
	return p;
}

Vector3D GenericNode::GetDisplacementD() const
{
	Vector3D p;
	double fact = mbs->GetDOption(105);

	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_displacement)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = fact*XGD(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = fact*XGD(i);
			else if (doftypes(i) & NDT_comp_3) p(3) = fact*XGD(i);
		}
		else if (doftypes(i) & NDT_position)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = fact*(XGD(i) - pos(1));
			else if (doftypes(i) & NDT_comp_2) p(2) = fact*(XGD(i) - pos(2));
			else if (doftypes(i) & NDT_comp_3) p(3) = fact*(XGD(i) - pos(3));
		}
	}
	return p;
}

Vector2D GenericNode::GetDisplacement2D() const
{
	Vector2D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_position)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XG(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XG(i);
		}
		else if (doftypes(i) & NDT_displacement)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XG(i) - pos(1);
			else if (doftypes(i) & NDT_comp_2) p(2) = XG(i) - pos(2);
		}
	}
	return p;
}

Vector2D GenericNode::GetDisplacement2DD() const
{
	Vector2D p;
	double fact = mbs->GetDOption(105);

	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_displacement)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = fact*XGD(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = fact*XGD(i);
		}
		else if (doftypes(i) & NDT_position)
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = fact*(XGD(i) - pos(1));
			else if (doftypes(i) & NDT_comp_2) p(2) = fact*(XGD(i) - pos(2));
		}
	}
	return p;
}

Vector3D GenericNode::GetVel() const
{
	Vector3D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_sos && (doftypes(i) & NDT_displacement || doftypes(i) & NDT_position))
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XGP(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XGP(i);
			else if (doftypes(i) & NDT_comp_3) p(3) = XGP(i);
		}
	}
	return p;
}

Vector3D GenericNode::GetVelD() const
{
	Vector3D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_sos && (doftypes(i) & NDT_displacement || doftypes(i) & NDT_position))
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XGPD(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XGPD(i);
			else if (doftypes(i) & NDT_comp_3) p(3) = XGPD(i);
		}
	}
	return p;
}

Vector2D GenericNode::GetVel2D() const
{
	Vector2D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_sos && (doftypes(i) & NDT_displacement || doftypes(i) & NDT_position))
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XGP(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XGP(i);
		}
	}
	return p;
}

Vector2D GenericNode::GetVel2DD() const
{
	Vector2D p;
	for (int i = 1; i <= doftypes.Length(); i++)
	{
		if (doftypes(i) & NDT_sos && (doftypes(i) & NDT_displacement || doftypes(i) & NDT_position))
		{
			if		(doftypes(i) & NDT_comp_1) p(1) = XGPD(i);
			else if (doftypes(i) & NDT_comp_2) p(2) = XGPD(i);
		}
	}
	return p;
}

void GenericNode::DrawNode()
{
	int res = mbs->GetIOption(147);
	double size = mbs->GetDOption(124);
	
	if (res == 0) return;

	if (size != 0)
	{
		mbs->SetColor(col);
		mbs->DrawSphere(GetPosD(),0.5*size,res);
	}
	if (mbs->GetIOption(148)) //show node numbers
	{
		Vector3D p = GetPosD();

		char text[40];
		sprintf(text,"%d", NodeNum());

		mbs->GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
	}
}

void Node2D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	GenericNode::GetElementData(edc);
}

int Node2D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = GenericNode::SetElementData(edc);

	return rv;
}

void Node3D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	GenericNode::GetElementData(edc);
}

int Node3D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = GenericNode::SetElementData(edc);

	return rv;
}

//void NodeBeam3D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
//{
//	GenericNode::GetElementData(edc);
//}
//
//int NodeBeam3D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
//{
//	int rv = GenericNode::SetElementData(edc);
//
//	return rv;
//}

void ANCFNodeS1rot1_3D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	Node3D::GetElementData(edc);
}

int ANCFNodeS1rot1_3D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = Node3D::SetElementData(edc);

	return rv;
}

void ANCFNodeS2S3_3D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	Node3D::GetElementData(edc);
}

int ANCFNodeS2S3_3D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = Node3D::SetElementData(edc);

	return rv;
}

void ANCFNodeS2_2D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	GenericNode::GetElementData(edc);
}

int ANCFNodeS2_2D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = GenericNode::SetElementData(edc);

	return rv;
}

void ANCFNodeS1S2_3D::GetElementData(ElementDataContainer& edc) 		//fill in all material data
{
	Node3D::GetElementData(edc);
}

int ANCFNodeS1S2_3D::SetElementData(ElementDataContainer& edc) //set material data according to ElementDataContainer
{
	int rv = Node3D::SetElementData(edc);

	return rv;
}