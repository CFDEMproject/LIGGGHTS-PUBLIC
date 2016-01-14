//#**************************************************************
//#
//# filename:             model_tools.cpp
//#
//# author:               Peter Gruber
//#
//# generated:						December 2012
//# description:          Modeling-functionality for automated node and element generation
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
//#**************************************************************


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "mbs_interface.h"
#include "element.h"
#include "material.h"
#include "body3d.h"
#include "ANCFBeam3DTorsion.h"
#include "node.h"

#include "model_tools.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++    Automated Node Generation     +++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Funktion mit StartPunkt, Orientierung und Länge, templatisiert für alle Knoten!!!!==>PG!!!!!! JG
IVector GenerateNodesOnStraightLine(MBS* mbs, const ANCFNodeS1rot1_3D& initial_node, const ANCFNodeS1rot1_3D& last_node, int number_of_nodes)
{
	int number_of_segments = number_of_nodes - 1;
	assert(number_of_segments>0 && "GenerateNodesOnArc: number of nodes > 1 required");
	assert( (initial_node.GetRefAngles() - last_node.GetRefAngles()).Norm() < 1e-14);   // rotation of initial and last node has to be equal

	Vector3D initial_node_pos = initial_node.Pos();
	Vector3D last_node_pos = last_node.Pos();
	Vector3D vector_to_next_node(last_node_pos - initial_node_pos); 
	
	Vector3D ref_kardan_angles(initial_node.GetRefAngles());

	double beam_length = vector_to_next_node.Norm();
	
	vector_to_next_node *= 1./number_of_segments;

	// create nodes
	IVector node_array(number_of_nodes);
	Vector3D pos(initial_node_pos);
	for (int i=0; i<=number_of_segments; i++)
	{
		ANCFNodeS1rot1_3D node(mbs);
		node.SetANCFNodeS1rot1_3D(pos,ref_kardan_angles);
		node_array(i+1) = mbs->AddNode(&node);
		pos += vector_to_next_node;
	}
	
	return node_array;
}

IVector GenerateNodesOnArc(MBS* mbs, const ANCFNodeS1rot1_3D& initial_node, const ANCFNodeS1rot1_3D& last_node, int number_of_nodes)
{
	int number_of_segments = number_of_nodes - 1;
	assert(number_of_segments>0 && "GenerateNodesOnArc: number of nodes > 1 required");

	Vector3D initial_node_e1, last_node_e1;
	Vector3D initial_node_e2, last_node_e2;
	Vector3D initial_node_e3, last_node_e3;
	for (int i=1; i<=3; i++)
	{
		initial_node_e1(i) = initial_node.GetLocalFrame(i,1);
		initial_node_e2(i) = initial_node.GetLocalFrame(i,2);
		initial_node_e3(i) = initial_node.GetLocalFrame(i,3);
		last_node_e1(i) = last_node.GetLocalFrame(i,1);
		last_node_e2(i) = last_node.GetLocalFrame(i,2);
		last_node_e3(i) = last_node.GetLocalFrame(i,3);
	}

	double delta_e3 = (initial_node_e3 - last_node_e3).Norm();
	assert(delta_e3 >= 1e-14 && "GenerateNodesOnArc: use GenerateNodesOnStraightLine(...) for generation of nodes on a straight line");

	double radius = (initial_node.Pos() - last_node.Pos()).Norm() / delta_e3;
	Vector3D center = initial_node.Pos() - radius*initial_node_e3;

	// test for consistency: 
	assert((last_node_e2-initial_node_e2).Norm() < 1e-14);  // by convention: e2 must be normal to the plane in which the arc lies

	// some renaming, in order to make the generation algorithm independent from interface (to be substituted for acceleration of computation)
	Vector3D tangential_vector_at_initial_node(initial_node_e1);
	Vector3D normal(initial_node_e2);
	Vector3D vector_to_initial_node(initial_node_e3);
	Vector3D vector_to_last_node(last_node_e3);
			
	double delta_phi = acos(vector_to_initial_node*vector_to_last_node)/number_of_segments;

	//create nodes
	IVector node_array(number_of_nodes);
	for (int i=0; i<=number_of_segments; i++)
	{
		Vector3D unit_vector_to_node_i(cos(i*delta_phi)*vector_to_initial_node + sin(i*delta_phi)*tangential_vector_at_initial_node);
		Vector3D tangential_vector_at_node_i(normal.Cross(unit_vector_to_node_i));

		Vector3D pos(center+radius*unit_vector_to_node_i);
		
		Matrix3D local_frame;
		for (int j=1; j<=3; j++) 
		{
			local_frame(j,1)=tangential_vector_at_node_i(j);
			local_frame(j,2)=normal(j);
			local_frame(j,3)=unit_vector_to_node_i(j);
		}
		
		Vector3D ref_angles;
		RotMatToKardanAngles(local_frame, ref_angles);
		ANCFNodeS1rot1_3D node(mbs); node.SetANCFNodeS1rot1_3D(pos, ref_angles);
		node_array(i+1) = mbs->AddNode(&node);
	}

	return node_array;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++    Automated Beam Generation     +++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVector GenerateANCFBeam3DTorsionBeam(MBS* mbs, const IVector& node_array, int matnr, const Vector3D& color)
{
	int number_of_elements = node_array.Length() - 1;
	assert(number_of_elements>0);

	//create elements
	IVector	element_array(number_of_elements);

	for (int i=1; i<=number_of_elements; i++)
	{
		ANCFBeam3DTorsion beam_element(mbs);
		beam_element.SetANCFBeam3DTorsion(node_array(i),node_array(i+1), matnr, color, 1, 1);
		element_array(i) = mbs->AddElement(&beam_element);
	}

	return element_array;
}

void ApplyLoadToElements(MBS* mbs, const IVector& element_array, const MBSLoad& load)
{
	for (int i=1; i<=element_array.Length(); i++)
	{
		int idx = element_array(i);
		assert(idx > 0);
		mbs->GetElementPtr(idx)->AddLoad(load);
	}
}