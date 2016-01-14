//#**************************************************************
//#
//# filename:             model_tools.h
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


#ifndef MODEL_TOOLS__H
#define MODEL_TOOLS__H


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++    Automated Node Generation     +++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IVector GenerateNodesOnStraightLine(MBS* mbs, const ANCFNodeS1rot1_3D& initial_node, const ANCFNodeS1rot1_3D& last_node, int number_of_nodes);
IVector GenerateNodesOnArc(MBS* mbs, const ANCFNodeS1rot1_3D& initial_node, const ANCFNodeS1rot1_3D& last_node, int number_of_nodes);



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++    Automated Beam Generation     +++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVector GenerateANCFBeam3DTorsionBeam(MBS* mbs, const IVector& node_array, int matnr, const Vector3D& color);
void ApplyLoadToElements(MBS* mbs, const IVector& element_array, const MBSLoad& load);

#endif //AUTOMATED_NODE_AND_ELEMENT_GENERATION__H
