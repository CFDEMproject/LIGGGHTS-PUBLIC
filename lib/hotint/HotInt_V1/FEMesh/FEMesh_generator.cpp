//#**************************************************************
//#
//# filename:             FEMesh_generator.cpp
//#
//# author:               Dorninger Alexander 
//#                       Daniel Reischl ( original version of GenerateCylinder )
//#                       Larissa Aigner ( original version of GenerateHexahedral deformed )
//#
//# generated:						November 2011
//# description:          Finite Element mesh
//# remarks:						  * original file split
//#                       * FEMesh_aux contains       variety of helper classes, 
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
 

#include "mbs_interface.h"
#include "stepsettings.h"
#include "FEMesh.h"
#include "myfile.h"
#include "femathhelperfunctions.h"


//#include "windows.h" //for shell execute
#include <direct.h>  //for getcwd
#include  <io.h>     //file operation _access
#include  <stdio.h>
#include  <stdlib.h>
#include "FEMesh.h"
#include "FEMesh_aux.h"

#include "mbsload.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ class FEMesh_Generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: encapsulate all kinds of generation functions

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ access functions - placed here since class FEMesh is used
//TODO: proper forward declaration & back to .h file
FEMesh* FEMesh_Generator::GetMesh() { return mesh; }

int FEMesh_Generator::NN() { return GetMesh()->NN(); }
Vector3D FEMesh_Generator::GetPoint3D(int i) { return GetMesh()->GetPoint3D(i); }
int FEMesh_Generator::AddNodeCheck(Vector3D pos, int domain, double tol) { return GetMesh()->AddNodeCheck(pos, domain, tol); }
int FEMesh_Generator::AddNodeCheck(Vector2D pos, int domain, double tol) { return GetMesh()->AddNodeCheck(pos, domain, tol); }

FEElement& FEMesh_Generator::GetElement(int i) {	return GetMesh()->GetElement(i); }
int FEMesh_Generator::AddElement(FEElement& fep) { return GetMesh()->AddElement(fep); }
TArray<FEMesh_Face>& FEMesh_Generator::Faces() {	return GetMesh()->Faces(); }

Box3D& FEMesh_Generator::GetNodeBox() { return GetMesh()->GetNodeBox(); } 
SearchTree& FEMesh_Generator::NodesTree() { return GetMesh()->NodesTree(); }
void FEMesh_Generator::ComputeNodeBox(IVector& subset = IVector(0)) { GetMesh()->ComputeNodeBox(subset); }
void FEMesh_Generator::ResizeNodeSearchTree(Box3D addbox, int addnodes) { GetMesh()->ResizeNodeSearchTree(addbox,addnodes); }

int FEMesh_Generator::NLoads() { return GetMesh()->NLoads(); }
int FEMesh_Generator::AddLoad(FEMesh_Load& loadi) { return GetMesh()->AddLoad(loadi); }

void FEMesh_Generator::LinearToQuadratic(IVector& subset = IVector(0)) {	GetMesh()->LinearToQuadratic(subset); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// generation of entire blocks 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Hexahedral Blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// extern calls with full set of parameters

// generates a block of hexaherals elements ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		MBSLoad& bodyload, int bodyloadactive,
																		IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();

	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
	CreateHexBlock_Hexes(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Hexes_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Hexes_AllFaces(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);
	
	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}

// generates a block of prism elements (hexahedrals split into 2 prisms each) ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock_Prisms(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		       IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		       IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       MBSLoad& bodyload, int bodyloadactive,
																		       IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		       int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();
	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
  CreateHexBlock_Prisms(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Prisms_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Prisms_AllFaces(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);

	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}

//$EK 2013-03-04
// generates a block of prism elements (hexahedrals split into 2 prisms each) ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock_Prisms_new(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		       IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		       IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       MBSLoad& bodyload, int bodyloadactive,
																		       IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		       int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();
	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
  CreateHexBlock_Prisms_new(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Prisms_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Prisms_AllFaces_new(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);

	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}

// generates a block of pyramid elements (hexahedrals split into 3 pyramids each) ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock_Pyrams(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		       IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		       IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       MBSLoad& bodyload, int bodyloadactive,
																		       IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		       int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();
	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
  CreateHexBlock_Pyrams(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Pyrams_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Pyrams_AllFaces(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);

	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}
//$EK 2012-03-04
// generates a block of reals pyramid elements (hexahedrals split into 3 pyramids each) ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock_Pyramids_new(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		       IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		       IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       MBSLoad& bodyload, int bodyloadactive,
																		       IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		       int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();
	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
  CreateHexBlock_Pyramids_new(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Pyrams_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Pyramids_AllFaces_new(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);

	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}


// generates a block of tetrahedral elements (hexahedrals split into 5 tetrahedrals each) ( including nodes, constraints, ... )
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock_Tetras(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		       IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces,			
																		       IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       MBSLoad& bodyload, int bodyloadactive,
																		       IVector& loadfaces_active, Vector& loadfaces_loadstrength,
																		       int order)
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();
	TArray<IVector*> blockouterelements; 

	CreateHexBlock_Nodes(blocknodes, block, divisions_xyz, bodynr);
  CreateHexBlock_Tetras(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Tetras_ComputeOuterElements(blockouterelements,divisions_xyz);
	CreateHexBlock_Tetras_AllFaces(blockfaces,blockelements,divisions_xyz);
	if(order==2) LinearToQuadratic(blockelements);

	if(constraintfaces_constainttypes.Length())
		CreateHexBlock_AllFaceConstraints(constraintfaces_constainttypes,constraintfaces_directions,blockfaces);
	if (bodyloadactive)
		CreateHexBlock_BodyLoads(bodyload, blockelements);
	if(loadfaces_active.Length())
		CreateHexBlock_LoadFaces(loadfaces_active,loadfaces_loadstrength,blockfaces,block);

	return 0;
}

// generates a DEFORMED hexahedral block of hexaherals elements (return values, no boundary conditions)
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock(TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces)			
{
	ReleaseArray_TemplatePtr(blockfaces);
	blockfaces.SetLen(6);
	for(int i=1; i <= 6; i++) blockfaces(i) = new IVector();

	CreateHexBlock_Nodes(blocknodes, corners, divisions_xyz, bodynr);
	CreateHexBlock_Hexes(blockelements, blocknodes, divisions_xyz, bodynr, matnr, color);
	CreateHexBlock_Hexes_AllFaces(blockfaces,blockelements,divisions_xyz);	

	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// abbreviated calls

// generates a block of hexaherals elements (no return values, no boundary conditions)
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color)
{
	return GenerateHexahedralBlock(block, divisions_xyz, bodynr, matnr, color, IVector(), IVector(), TArray<IVector*>(), 
		                      IVector(), IVector(), MBSLoad(), 0, IVector(), Vector(), 1);
}

// generates a block of hexaherals elements (return values, no boundary conditions)
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces)			
{
	return GenerateHexahedralBlock(block, divisions_xyz, bodynr, matnr, color, blockelements, blocknodes, blockfaces, 
		                      IVector(), IVector(), MBSLoad(), 0, IVector(), Vector(), 1);
}

// generates a DEFORMED hexahedral block of hexaherals elements (return values, no boundary conditions)
// this may only make sense in a Mesh3D
int FEMesh_Generator::GenerateHexahedralBlock(TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr, int matnr, Vector3D color)           
{
	return GenerateHexahedralBlock(corners, divisions_xyz, bodynr, matnr, color, IVector(0), IVector(0), TArray<IVector*>(0));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create nodes

// creates nodes for hexahedral block - adds nodes to mesh.points
int FEMesh_Generator::CreateHexBlock_Nodes(IVector& blocknodes, Box3D block, int3 divisions_xyz, int bodynr, int order)
{
// resize nodestree
	ComputeNodeBox();
	Box3D newbox(GetNodeBox());
	if(! (newbox.IsIn(block.PMax()) && newbox.IsIn(block.PMin())) )
	{
		newbox.Add(block);
		int additional_nodes = (divisions_xyz(1)+1)*(divisions_xyz(2)+1)*(divisions_xyz(3)+1);
		int divs = (int) pow( (NN()+additional_nodes)*0.1 ,1/3.) +1;
		NodesTree().ResetSearchTree(divs, divs, divs, newbox);
		for(int i=1; i<=NN(); i++) 
			NodesTree().AddItem(Box3D(GetPoint3D(i),1e-8),i);
	}

	blocknodes.Flush();
// create all nodes for the block, do not make any double nodes, save nodenumbers of the nodes of the block in IVector
	for (int k=0; k <= divisions_xyz(3); k++)
	{
		double fz = (double)(k)/(double)divisions_xyz(3);
	  for (int j=0; j <= divisions_xyz(2); j++)
		{
			double fy = (double)(j)/(double)divisions_xyz(2);
			for (int i=0; i <= divisions_xyz(1); i++)
			{
				double fx = (double)(i)/(double)divisions_xyz(1);
				Vector3D pos = block.PMin() + Vector3D( block.SizeX()*fx, block.SizeY()*fy, block.SizeZ()*fz );
				
				blocknodes.Add(AddNodeCheck(pos,bodynr,1E-10));		
			}
		}
	}
	return blocknodes.Length();
}

// creates nodes for DEFORMED hexahedral block - adds nodes to mesh.points
int FEMesh_Generator::CreateHexBlock_Nodes(IVector& blocknodes, TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr, int order)
{
// resize nodestree
	ComputeNodeBox();
	Box3D newbox(GetNodeBox());
	int contained = 1;
	for(int i=1; i<=corners.Length(); i++)
	{
		if(!newbox.IsIn(corners.Get(i)))
		{
			newbox.Add(corners.Get(i));
			contained = 0;
		}
	}
	if(!contained)
	{
		int additional_nodes = (divisions_xyz(1)+1)*(divisions_xyz(2)+1)*(divisions_xyz(3)+1);
		int divs = (int) pow( (NN()+additional_nodes)*0.1 ,1/3.) +1;
		NodesTree().ResetSearchTree(divs, divs, divs, newbox);
		for(int i=1; i<=NN(); i++) 
			NodesTree().AddItem(Box3D(GetPoint3D(i),1e-8),i);
	}

	blocknodes.Flush();
	// create all nodes for the hexahedral with 'corners', do not make any double nodes, save nodenumbers of the nodes of the block in IVector
	for (int k=0; k <= divisions_xyz(3); k++)
	{
		double fz = (double)(k)/(double)divisions_xyz(3);
	  for (int j=0; j <= divisions_xyz(2); j++)
		{
			double fy = (double)(j)/(double)divisions_xyz(2);
			for (int i=0; i <= divisions_xyz(1); i++)
			{
				double fx = (double)(i)/(double)divisions_xyz(1);
				Vector3D pos_xy0 = corners.Get(1) + (corners.Get(2)-corners.Get(1))*fx +
												  (corners.Get(3)+(corners.Get(4)-corners.Get(3))*fx - (corners.Get(1)+(corners.Get(2)-corners.Get(1))*fx))*fy;
				Vector3D pos_xyz = corners.Get(5) + (corners.Get(6)-corners.Get(5))*fx +
												  (corners.Get(7)+(corners.Get(8)-corners.Get(7))*fx - (corners.Get(5)+(corners.Get(6)-corners.Get(5))*fx))*fy;
				Vector3D pos = pos_xy0 + (pos_xyz - pos_xy0)*fz;
											 				
				blocknodes.Add(AddNodeCheck(pos,bodynr,1E-10));		
			}
		}
	}
	return blocknodes.Length();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create elements

// returns Nodelist of a single (hexahedral) element of the block
int FEMesh_Generator::CreateHexBlock_HexNodes(IVector& points, IVector& blocknodes, int firstnode, int dx, int dy, int dz, int order)
{
// node numbering in HotInt sequence:
	if(order == 1)
	{
		points.SetLen(8);
		points(1) = blocknodes( firstnode                              );
		points(2) = blocknodes( firstnode                          + 1 );
		points(3) = blocknodes( firstnode                 + (dx+1)     );
		points(4) = blocknodes( firstnode                 + (dx+1) + 1 );
		points(5) = blocknodes( firstnode + (dy+1)*(dx+1)              );
		points(6) = blocknodes( firstnode + (dy+1)*(dx+1)          + 1 );
		points(7) = blocknodes( firstnode + (dy+1)*(dx+1) + (dx+1)     );
		points(8) = blocknodes( firstnode + (dy+1)*(dx+1) + (dx+1) + 1 );
		return 8;
	}
	return points.Length();
}

// creates Hexahedrals from nodelist - adds FEHex to mesh.elements
int FEMesh_Generator::CreateHexBlock_Hexes(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(8);

	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 1 hex
				FEHex fehex(points);
				fehex.MaterialNum() = matnr;
				fehex.Domain() = bodynr;
				fehex.Color() = color;
				blockelements.Add(AddElement(fehex));
			}
		}
	}
	return blockelements.Length();
}

// creates Prisms from nodelist - adds 2 FEHex to mesh.elements
int FEMesh_Generator::CreateHexBlock_Prisms(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(8);
	TArray<IVector*> prisms(2); 
	for(int i=1; i<=2; i++) prisms(i) = new IVector(0);
	TArray<Vector3D> colorblend; colorblend.Set2(Vector3D(1.0,1.0,1.0),Vector3D(0.0,0.0,0.0));

	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 2prisms
				prisms(1)->SetXN(8, points(1),points(2),points(3),points(3),points(5),points(6),points(7),points(7));
				prisms(2)->SetXN(8, points(4),points(3),points(2),points(2),points(8),points(7),points(6),points(6)); 
// swap every 2nd prism ? -> not necessary
				for(int i=1; i<=2; i++)
				{
					FEHex fepri(*prisms(i));
					fepri.MaterialNum() = matnr;
					fepri.Domain() = bodynr;
					fepri.Color() = (color + colorblend(i)) * 0.5;
					blockelements.Add(AddElement(fepri));
				}
			}
		}
	}
	for(int i=1; i<=2; i++) delete prisms(i);
	return blockelements.Length();
}

// $EK 2013-03-04
// creates Prisms from nodelist - adds 2 FEPrism to mesh.elements
int FEMesh_Generator::CreateHexBlock_Prisms_new(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(6);
	TArray<IVector*> prisms(2); 
	for(int i=1; i<=2; i++) prisms(i) = new IVector(0);
	TArray<Vector3D> colorblend; colorblend.Set2(Vector3D(1.0,1.0,1.0),Vector3D(0.0,0.0,0.0));

	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 2prisms
				//prisms(1)->SetXN(8, points(1),points(2),points(3),points(3),points(5),points(6),points(7),points(7));
				//prisms(2)->SetXN(8, points(4),points(3),points(2),points(2),points(8),points(7),points(6),points(6)); 
				prisms(1)->SetXN(6, points(1),points(2),points(3),points(5),points(6),points(7));
				prisms(2)->SetXN(6, points(4),points(3),points(2),points(8),points(7),points(6));
// swap every 2nd prism ? -> not necessary
				for(int i=1; i<=2; i++)
				{
					//FEHex fepri(*prisms(i));
					FEPrism fepri(*prisms(i));
					fepri.MaterialNum() = matnr;
					fepri.Domain() = bodynr;
					fepri.Color() = (color + colorblend(i)) * 0.5;
					blockelements.Add(AddElement(fepri));
				}
			}
		}
	}
	for(int i=1; i<=2; i++) delete prisms(i);
	return blockelements.Length();
}

// creates Pyramids from nodelist - adds 3 FEHex to mesh.elements
int FEMesh_Generator::CreateHexBlock_Pyrams(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(8);
	TArray<IVector*> pyrams(3); 
	for(int i=1; i<=3; i++) pyrams(i) = new IVector(0);
	TArray<Vector3D> colorblend; colorblend.Set3(Vector3D(1.0,0.0,0.0),Vector3D(0.0,1.0,0.0),Vector3D(0.0,0.0,1.0));
	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 3 pyramids
				if ( (ix+iy+iz)%2 == 1)
				{
					pyrams(1)->SetXN(8, points(1),points(2),points(3),points(4),points(5),points(5),points(5),points(5));
					pyrams(2)->SetXN(8, points(6),points(8),points(2),points(4),points(5),points(5),points(5),points(5)); 
					pyrams(3)->SetXN(8, points(7),points(3),points(8),points(4),points(5),points(5),points(5),points(5)); 
				}
				else // every 2nd hex swapped for continuation
				{
					pyrams(1)->SetXN(8, points(8),points(7),points(6),points(5),points(4),points(4),points(4),points(4));
					pyrams(2)->SetXN(8, points(3),points(1),points(7),points(5),points(4),points(4),points(4),points(4));
					pyrams(3)->SetXN(8, points(2),points(6),points(1),points(5),points(4),points(4),points(4),points(4));
				
				}
				for(int i=1; i<=3; i++)
				{
					FEHex fepyr(*pyrams(i));
					fepyr.MaterialNum() = matnr;
					fepyr.Domain() = bodynr;
					fepyr.Color() = (color + colorblend(i)) * 0.5;
					blockelements.Add(AddElement(fepyr));
				}
			}
		}
	}
	for(int i=1; i<=3; i++) delete pyrams(i);
	return blockelements.Length();
}

//$EK 2013-03-04 added for FEPyramids....
// creates Pyramids from nodelist - adds 3 FEPyramid to mesh.elements
int FEMesh_Generator::CreateHexBlock_Pyramids_new(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(5);
	TArray<IVector*> pyrams(3); 
	for(int i=1; i<=3; i++) pyrams(i) = new IVector(0);
	TArray<Vector3D> colorblend; colorblend.Set3(Vector3D(1.0,0.0,0.0),Vector3D(0.0,1.0,0.0),Vector3D(0.0,0.0,1.0));
	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 3 pyramids
				if ( (ix+iy+iz)%2 == 1)
				{
					//pyrams(1)->SetXN(8, points(1),points(2),points(3),points(4),points(5),points(5),points(5),points(5));
					//pyrams(2)->SetXN(8, points(6),points(8),points(2),points(4),points(5),points(5),points(5),points(5)); 
					//pyrams(3)->SetXN(8, points(7),points(3),points(8),points(4),points(5),points(5),points(5),points(5)); 
					pyrams(1)->SetXN(5, points(1),points(2),points(3),points(4),points(5));
					pyrams(2)->SetXN(5, points(6),points(8),points(2),points(4),points(5)); 
					pyrams(3)->SetXN(5, points(7),points(3),points(8),points(4),points(5)); 

				}
				else // every 2nd hex swapped for continuation
				{
					pyrams(1)->SetXN(5, points(8),points(7),points(6),points(5),points(4));
					pyrams(2)->SetXN(5, points(3),points(1),points(7),points(5),points(4));
					pyrams(3)->SetXN(5, points(2),points(6),points(1),points(5),points(4));
				
				}
				for(int i=1; i<=3; i++)
				{
					FEPyramid fepyr(*pyrams(i));
					fepyr.MaterialNum() = matnr;
					fepyr.Domain() = bodynr;
					fepyr.Color() = (color + colorblend(i)) * 0.5;
					blockelements.Add(AddElement(fepyr));
				}
			}
		}
	}
	for(int i=1; i<=3; i++) delete pyrams(i);
	return blockelements.Length();
}


// creates Tetrahedralss from nodelist - adds 5 FETet to mesh.elements
int FEMesh_Generator::CreateHexBlock_Tetras(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(8);
	TArray<IVector*> tets(5);
	for(int i=1; i<=5; i++) tets(i) = new IVector(0);
	TArray<Vector3D> colorblend; colorblend.SetXN(5,Vector3D(1.0,0.0,0.0),Vector3D(0.0,1.0,0.0),Vector3D(0.0,0.0,1.0),Vector3D(1.0,1.0,1.0),Vector3D(0.0,0.0,0.0));

	int firstnode;
	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);

 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				firstnode = (iz-1)*(dy+1)*(dx+1) + (iy-1)*(dx+1) + (ix-1)*1 + 1;
				CreateHexBlock_HexNodes(points,blocknodes,firstnode,dx,dy,dz,order);
// 5 tetrahedrals 
				if ( (ix+iy+iz)%2 == 1)
				{
					tets(1)->SetXN(5,points(1),points(2),points(3),points(5));
					tets(2)->SetXN(5,points(4),points(3),points(2),points(8));
					tets(3)->SetXN(5,points(6),points(5),points(8),points(2));
					tets(4)->SetXN(5,points(7),points(5),points(8),points(3));
					tets(5)->SetXN(5,points(2),points(5),points(3),points(8));
				}
				else // every 2nd hex swapped for continuation
				{
					tets(1)->SetXN(5,points(8),points(7),points(6),points(4));
					tets(2)->SetXN(5,points(5),points(6),points(7),points(1));
					tets(3)->SetXN(5,points(3),points(4),points(1),points(7));
					tets(4)->SetXN(5,points(2),points(4),points(1),points(6));
					tets(5)->SetXN(5,points(7),points(4),points(6),points(1));
				}

				for(int i=1; i<=5; i++)
				{
					FETet fetet(*tets(i));
					fetet.MaterialNum() = matnr;
					fetet.Domain() = bodynr;
					fetet.Color() = (color + colorblend(i)) * 0.5;
					blockelements.Add(AddElement(fetet));
				}
			}
		}
	}
	for(int i=1; i<=5; i++) delete tets(i) ;
	return blockelements.Length();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Cylindrical Blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// extern calls with full set of parameters

// generates a 1/4 disc of hexaherals elements, with or without hole in the center
// this may only make sense in a Mesh3D, use MirrorMeshAtPlane to get the rest of the cylinder
int FEMesh_Generator::GenerateDisc(double radius, double radius_hole, double height, int divisions_angle, int divisions_radial, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
{ 
	IVector blockelements_tmp, blocknodes_tmp;
	blockelements.Flush();
	blocknodes.Flush();
	double radius_1, radius_2;
	if(!divisions_radial) assert(0);

	if(!radius_hole)						// disc withot hole
	{
		radius_hole = radius/(divisions_radial);
		GenerateCylinder(radius_hole, radius_hole/2.0, height, divisions_angle, divisions_height, pos, bodynr, matnr, color, blockelements_tmp, blocknodes_tmp);
		blockelements.Merge(blockelements_tmp);
		blocknodes.Merge(blocknodes_tmp);
		divisions_radial = divisions_radial-1;
	}

	double radius_incr = (radius-radius_hole)*(1.0/divisions_radial);
	int j=2;
	int k=1;
	radius_2 = radius_hole;

	if(divisions_radial)
	{
		for (int i=1; i<=divisions_radial;i++)
		{
			radius_1 = radius_2;
			radius_2 = radius_1 + radius_incr;	
			if(i==k)
			{
				GenerateHollowCylinderDoubleDiv(radius_2, radius_1,height, divisions_angle,divisions_height, pos, bodynr, matnr, color, blockelements_tmp, blocknodes_tmp);
				blockelements.Merge(blockelements_tmp);
				blocknodes.Merge(blocknodes_tmp);
				divisions_angle=divisions_angle*2;
				j = j+1;
				k = k+j;
			}
			else
			{
				GenerateHollowCylinder(radius_2, radius_1,height, divisions_angle,divisions_height, pos, bodynr, matnr, color, blockelements_tmp, blocknodes_tmp);
				blockelements.Merge(blockelements_tmp);
				blocknodes.Merge(blocknodes_tmp);
			}
		}
	}
	RemoveRedundantEntries(blocknodes);
	return 0;
}

// generates a 1/4 cylinder of hexaherals elements 
// this may only make sense in a Mesh3D, use MirrorMeshAtPlane to get the rest of the cylinder
int FEMesh_Generator::GenerateCylinder(double radius, double small_radius, double height, int divisions_angle, int divisions_height, Vector3D pos, 
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
{
	/*
	Input:
	radius:			the radius of the cylinder
	small_radius:	the radius where the mesh is changing from meshing a square to meshing a circle
	*/

	Box3D addbox(pos, pos+Vector3D(radius, radius, height));
	int additional_nodes = (3*divisions_angle*divisions_angle + 3*divisions_angle+ 1 )*(divisions_height+1); // dH*nN (see below)
	ResizeNodeSearchTree(addbox, additional_nodes);
	blocknodes.Flush();
	TArray<double> x,y;

	CreateCylinder_Mesh2D(radius, small_radius, divisions_angle,x,y);
	Create3DNodesOutOf2DCoordinates(x, y, height, divisions_height, pos, bodynr, blocknodes);
	CreateCylinder_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);
	
	return 0;
}

// generates a 1/4 hollow cylinder of hexaherals elements 
// this may only make sense in a Mesh3D, use MirrorMeshAtPlane to get the rest of the cylinder
int FEMesh_Generator::GenerateHollowCylinder(double radius, double radius_hole, double height, int divisions_angle, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
{
	Box3D addbox(pos, pos+Vector3D(radius, radius, height));
	int additional_nodes = (4*divisions_angle+2)*(divisions_height+1); // dH*nN (see below)
	ResizeNodeSearchTree(addbox, additional_nodes);
	blocknodes.Flush();
	TArray<double> x,y;

	CreateHollowCylinder_Mesh2D(radius, radius_hole, divisions_angle,x,y);
	Create3DNodesOutOf2DCoordinates(x, y, height, divisions_height, pos, bodynr, blocknodes);
	CreateHollowCylinder_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);

	return 0;
}

int FEMesh_Generator::GenerateCylinder(Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r, int divisions_angle, int divisions_height, 
																			 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
{
	//				x-axis					
	//					|								
	//					|							
	//					|			 p1l		p1r	
	//					|			/________\		
	//					| 		|#########|
	//					|			|#########|
	//					|	p2l	|#########| p2r 
	//					|		 \|#########|/
	//					|
	//	z-axis	-.-.-.-.-.-.-.-.-.-.-.-.-.-.->

	// if p2l.X==p2r.X == 0, than a full cylinder is generated
	// if both, p2l.X AND p2r.X, are != 0 than a hollow cylinder is generated
	// it is possible to define a distorted cylinder: p1l.X != p1r.X, p1l.Z != p2l.Z, ...

	Box3D addbox(Vector3D(0.,0.,p2l.Z()), Vector3D(p1r.X(),p1r.X(),p1r.Z()));
	int additional_nodes;
	blocknodes.Flush();
	TArray<double> x,y;
	int hollowCylinder = 0;

	if(p2l.X() && p2r.X()) // hollow cylinder
	{
		additional_nodes = (4*divisions_angle+2)*(divisions_height+1); 
		hollowCylinder = 1;
	}
	else if((p2l.X()==0) && (p2r.X()==0))	//solid cylinder
	{
		additional_nodes = (3*divisions_angle*divisions_angle + 3*divisions_angle+ 1 )*(divisions_height+1);
	}
	else	// bad points
	{
		return -1;
	}
	ResizeNodeSearchTree(addbox, additional_nodes);

	if(hollowCylinder) // hollow cylinder
	{
		CreateHollowCylinder_Mesh2D(p1l.X(), p2l.X(), divisions_angle,x,y);
		Create3DNodesOutOf2DCoordinates(x, y, p2r.Z()-p2l.Z(), divisions_height,bodynr,blocknodes,p1l,p1r,p2l,p2r);
		CreateHollowCylinder_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);
	}
	else	//solid cylinder
	{
		CreateCylinder_Mesh2D(p1l.X(), p1l.X()/2., divisions_angle,x,y);
		Create3DNodesOutOf2DCoordinates(x, y, p2r.Z()-p2l.Z(), divisions_height,bodynr,blocknodes,p1l,p1r,p2l,p2r);
		CreateCylinder_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);
	}
	return 0;
}


// generates an 1/4 hollow cylinder of hexaherals elements with increasing number of divisions/rad
// this may only make sense in a Mesh3D, use MirrorMeshAtPlane to get the rest of the cylinder
int FEMesh_Generator::GenerateHollowCylinderDoubleDiv(double radius, double radius_hole, double height, int divisions_angle, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
{
	Box3D addbox(pos, pos+Vector3D(radius, radius, height));
	int additional_nodes = (9*divisions_angle+2)*(divisions_height+1); // dH*nN (see below)
	ResizeNodeSearchTree(addbox, additional_nodes);
	blocknodes.Flush();
	TArray<double> x,y;

	CreateHollowCylinderDoubleDiv_Mesh2D(radius, radius_hole, divisions_angle,x,y);
	Create3DNodesOutOf2DCoordinates(x, y, height, divisions_height, pos, bodynr, blocknodes);
	CreateHollowCylinderDoubleDiv_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);
	
	return 0;
}

	int FEMesh_Generator::GenerateHollowCylinderDoubleDiv(Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r, int divisions_angle, int divisions_height, 
								int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes)
	{
	//				x-axis					
	//					|								
	//					|							
	//					|			 p1l		p1r	
	//					|			/________\		
	//					| 		|#########|
	//					|			|#########|
	//					|	p2l	|#########| p2r 
	//					|		 \|#########|/
	//					|
	//	z-axis	-.-.-.-.-.-.-.-.-.-.-.-.-.-.->

	// it is possible to define a distorted cylinder: p1l.X != p1r.X, p1l.Z != p2l.Z, ...

	Box3D addbox(Vector3D(0.,0.,p2l.Z()), Vector3D(p1r.X(),p1r.X(),p1r.Z()));
	int additional_nodes; 
	blocknodes.Flush();
	TArray<double> x,y;

	if(p2l.X() && p2r.X()) // hollow cylinder
	{
		additional_nodes = (9*divisions_angle+2)*(divisions_height+1);
	}
	else	// bad points
	{
		return -1;
	}
	ResizeNodeSearchTree(addbox, additional_nodes);

	CreateHollowCylinderDoubleDiv_Mesh2D(p1l.X(), p2l.X(), divisions_angle,x,y);
	Create3DNodesOutOf2DCoordinates(x, y, p2r.Z()-p2l.Z(), divisions_height,bodynr,blocknodes,p1l,p1r,p2l,p2r);
	CreateHollowCylinderDoubleDiv_Hexes(blockelements, blocknodes, divisions_angle, divisions_height, bodynr, matnr, color);

	return 0;
	}

// generates a rotor based on a mesh2d by rotating the mesh2d around the specified rotation axis
//$ DR 2013-03-28 added this function
int FEMesh_Generator::GenerateRotorOutOfMesh2D(FEMesh* mesh2d, int rotation_axis_number, int angular_segments, double final_angle_deg)
{
	FEMesh2D *meshP = (FEMesh2D *)mesh2d;

  // strategy
	// **********************************************************************
	// 1: create the nodes for each angle/plane by means of a rotation matrix
	// **********************************************************************
	TArrayDynamic<FEMesh_Set> nodes_per_plane;

	// loop over all required planes
	for(int i=1; i<=angular_segments+1; i++)
	{
		FEMesh_Set nodes(TSetNodes);

		// compute the rotation matrix for this plane		
		double angle = (MY_PI/180.) * (final_angle_deg/angular_segments) * (i-1.);
		Matrix3D rot;
		if (rotation_axis_number == 1) rot = RotMatrix1(angle);
		else if (rotation_axis_number == 2) rot = RotMatrix2(angle);
		else /*if (rotation_axis_number == 3)*/ rot = RotMatrix3(angle);


		for(int j=1; j<=meshP->NNodes(); j++)
		{
			FEMesh_Node& node2d = meshP->GetNode(j);
			Vector3D nodecoord = node2d.GetCoords3D();
			Vector3D nodecoord_angle = rot * nodecoord;
			int nodenr = AddNodeCheck(nodecoord_angle, node2d.Domain());
			nodes.Add(nodenr);
		}
		nodes_per_plane.Add(nodes);
	}
	// All Nodes for sweep are now created...

	// ***********************************************************************************
	// 2: loop throug the elements of the 2d mesh and create the corresponding 3d elmeents
	// ***********************************************************************************
	// outer loop: all elements
	for(int i=1;i<=meshP->NElements(); i++)
	{
		// inner loop: all corresponding elements
		for(int j=1; j<=angular_segments; j++)
		{
			if(meshP->GetElement(i).Type() == TFETrig)
			{
				FETrig& the_2d_element = (FETrig&) meshP->GetElement(i);
				int planenumber1 = j;
				int planenumber2 = j+1; 

				FEMesh_Set& nodes_p1 = nodes_per_plane(planenumber1);
				FEMesh_Set& nodes_p2 = nodes_per_plane(planenumber2);

				// check if it is a tet, pyramid or prism
				int flag_nodesOnAxis=0;
				int NnodesOnAxis=0;

				if(nodes_p1.Node(the_2d_element.GetNode(1))==nodes_p2.Node(the_2d_element.GetNode(1)))
				{
					flag_nodesOnAxis += 1;
					NnodesOnAxis++;
				}
				if(nodes_p1.Node(the_2d_element.GetNode(2))==nodes_p2.Node(the_2d_element.GetNode(2)))
				{
					flag_nodesOnAxis += 2;
					NnodesOnAxis++;
				}
				if(nodes_p1.Node(the_2d_element.GetNode(3))==nodes_p2.Node(the_2d_element.GetNode(3)))
				{
					flag_nodesOnAxis += 4;
					NnodesOnAxis++;
				}

				// create FE elements
				if(NnodesOnAxis == 0)	// no point on axis --> prism
				{
					TArray<int> points_prism;
					points_prism.Add( nodes_p1.Node( the_2d_element.GetNode(1) ) );
					points_prism.Add( nodes_p1.Node( the_2d_element.GetNode(2) ) );
					points_prism.Add( nodes_p1.Node( the_2d_element.GetNode(3) ) );

					points_prism.Add( nodes_p2.Node( the_2d_element.GetNode(1) ) );
					points_prism.Add( nodes_p2.Node( the_2d_element.GetNode(2) ) );
					points_prism.Add( nodes_p2.Node( the_2d_element.GetNode(3) ) );

					FEPrism* the_prism = (FEPrism*) mesh->MakePrism(points_prism, the_2d_element.MaterialNum());
					the_prism->Domain() = the_2d_element.Domain();
					mesh->AddElement(the_prism);
				}
				else if(NnodesOnAxis == 1)	// 1 point on axis --> pyramid
				{
					TArray<int> points_pyr;

					if(flag_nodesOnAxis & 1)					// first node on axis 
					{ 
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(2) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(3) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(2) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(3) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(1) ) );	
					}
					else if(flag_nodesOnAxis & 2)			// second node on axis 
					{	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(1) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(3) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(1) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(3) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(2) ) );	
					}
					else															// third node on axis 
					{	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(1) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(2) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(1) ) );	
						points_pyr.Add( nodes_p2.Node( the_2d_element.GetNode(2) ) );	
						points_pyr.Add( nodes_p1.Node( the_2d_element.GetNode(3) ) );	
					}

					FEPyramid* the_pyr = (FEPyramid*)mesh->MakePyramid(points_pyr, the_2d_element.MaterialNum());
					the_pyr->Domain() = the_2d_element.Domain();
					mesh->AddElement(the_pyr);
				}
				else if(NnodesOnAxis == 2)	// 2 points on axis --> tetrahedral
				{
					TArray<int> points_tet;
					points_tet.Add( nodes_p1.Node( the_2d_element.GetNode(1) ) );
					points_tet.Add( nodes_p1.Node( the_2d_element.GetNode(2) ) );
					points_tet.Add( nodes_p1.Node( the_2d_element.GetNode(3) ) );

					if(!(flag_nodesOnAxis & 1))				// first node not on axis 
					{ points_tet.Add( nodes_p2.Node( the_2d_element.GetNode(1) ) );	}
					else if(!(flag_nodesOnAxis & 2))	// second node not on axis 
					{	points_tet.Add( nodes_p2.Node( the_2d_element.GetNode(2) ) );	}
					else															// third node not on axis 
					{	points_tet.Add( nodes_p2.Node( the_2d_element.GetNode(3) ) );	}

					FETet* the_tet = (FETet*)mesh->MakeTet(points_tet, the_2d_element.MaterialNum());
					the_tet->Domain() = the_2d_element.Domain();
					mesh->AddElement(the_tet);
				}
			}
			else /*if(the_2d_element.Type() == TFEQuad)*/
			{

			}
		}
	}

	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// abbreviated calls


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create nodes

// create 3D nodes based on 2D mesh. Used by functions GenerateCylinder, GenerateHollowCylinder, GenerateHollowCylinderDoubleDiv 	- adds nodes to mesh.points
// mesh stored in x and y is duplicated and scaled (divisions_height+1)-times beginning with z-coordinate pos.Z up to pos.Z+height
// origin is shifted to (pos.X,pos.Y, pos.Z)
//
// if the MathFunction* are provided:
// it is assumed that the MathFunctions describe the surface of a (hollow) rotor
// the coordinates are interpolated between the values of the MathFunctions
// "left" indicate scaling factor at positions with z = pos.Z, "right" means positions at pos.Z+height
// WARNING: do not use deformation in both directions, inner-outer AND left-right. This leads to a wrong geometry
// If the surfaces are just linear functions, use the other call of Create3DNodesOutOf2DCoordinates

int FEMesh_Generator::Create3DNodesOutOf2DCoordinates(TArray<double>& x, TArray<double>& y, double height, int divisions_height, Vector3D pos, int bodynr, IVector& blocknodes,MathFunction* outer, MathFunction* inner, MathFunction* left,MathFunction* right)
{
	double	h	= height;
	int		dh	= divisions_height;
	double	h0	= h/dh;		// increment for the height
	double	h_tmp;

	double	x_offset=pos.X();
	double	y_offset=pos.Y();
	double	z_offset=pos.Z();
	int nN = x.Length();
	TArray<Vector3D> nodes_tmp;
	Vector3D pos_tmp;

	// generate the 3D-coordinates of a straight geometry
	for (int j=0; j<=dh; j++)
	{
		h_tmp=j*h0;
		for (int i=1; i<=nN; i++)
		{
			pos_tmp = Vector3D(x(i), y(i), h_tmp)+pos;
			nodes_tmp.Add(pos_tmp);
		}
	}

	// distort the temporary mesh in radial direction
	if(outer != NULL)	
	{
		MathFunction* mf_out= new MathFunction;
		MathFunction* mf_in= new MathFunction;

		mf_out = outer;

		if(inner != NULL) {	mf_in = inner; }
		else 							{	mf_in->SetConstant(0.);	}


		// find the maximum radius
		double ref_out=0;
		for (int i=1; i<=nN; i++)
		{
			if(x(i) > ref_out) {ref_out = x(i);}
		}

		// find the minimum radius
		double ref_in = nodes_tmp(1).X();		//x-coordinate of first node corresponds to radius of hole

		Distort3DMesh(nodes_tmp,3,1,mf_out,ref_out,mf_in,ref_in,2);
	}

	// distort the temporary mesh in axial direction
	if(left != NULL)
	{
		double ref_left = nodes_tmp(1).Z();			// axial coordinate of first node
		double ref_right = nodes_tmp(nodes_tmp.Length()).Z();	// axial coordinate of last node
		//Distort3DMesh(nodes_tmp,1,3,right,ref_right,left,ref_left,0,2);
		Distort3DMesh(nodes_tmp,1,3,right,ref_right,left,ref_left,0,2);
	}

	// Create nodes
	for (int i=1; i<=nodes_tmp.Length(); i++)
	{
		//nodes_tmp(i) = nodes_tmp(i) + Vector3D(x_offset,y_offset,0);
		blocknodes.Add(AddNodeCheck(nodes_tmp(i),bodynr,1E-10));
	}

	return 0;
}

int FEMesh_Generator::Create3DNodesOutOf2DCoordinates(TArray<double>& x, TArray<double>& y, double height, int divisions_height, int bodynr,IVector& blocknodes,Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r)
{
	double	h	= height;
	int		dh	= divisions_height;
	double	h0	= h/dh;		// increment for the height
	double	h_tmp;

	Vector3D pos(0,0,p2l.Z());

	double	x_offset=pos.X();
	double	y_offset=pos.Y();
	double	z_offset=pos.Z();
	int nN = x.Length();
	TArray<Vector3D> nodes_tmp;
	Vector3D pos_tmp;

	// generate the 3D-coordinates of a straight geometry
	for (int j=0; j<=dh; j++)
	{
		h_tmp=j*h0;
		for (int i=1; i<=nN; i++)
		{
			pos_tmp = Vector3D(x(i), y(i), h_tmp)+pos;
			nodes_tmp.Add(pos_tmp);
		}
	}

	// distort the temporary mesh in radial direction

		MathFunction* mf_out= new MathFunction;
		MathFunction* mf_in= new MathFunction;

		//mf_out->SetPiecewise(Vector(p1l.Z(),p1r.Z()),Vector(p1l.X(),p1r.X()),1);		// linear interpolation of outer surface
		mf_out->SetPiecewise(Vector(p2l.Z(),p2r.Z()),Vector(p1l.X(),p1r.X()),1);		// linear interpolation of pseudo outer surface, used for first distortion step
		mf_in->SetPiecewise(Vector(p2l.Z(),p2r.Z()),Vector(p2l.X(),p2r.X()),1);			// linear interpolation of inner surface

		// find the maximum radius
		double ref_out=0;
		for (int i=1; i<=nN; i++)
		{
			if(x(i) > ref_out) {ref_out = x(i);}
		}

		// find the minimum radius
		double ref_in = nodes_tmp(1).X();		//x-coordinate of first node corresponds to radius of hole

		Distort3DMesh(nodes_tmp,3,1,mf_out,ref_out,mf_in,ref_in,2);


	// distort the temporary mesh in axial direction

		if((p1l.Z()!=p2l.Z())||(p1r.Z()!=p2r.Z()))
		{

			MathFunction* right= new MathFunction;
			MathFunction* ref_right= new MathFunction;
			MathFunction* left= new MathFunction;
			MathFunction* ref_left= new MathFunction;

			if(p2l.X()==p2r.X())	//straight inner surface
			{
				if(p1l.X()==p1r.X())	//straigth outer surface + straight inner surface
				{
					ref_left->SetConstant(p2l.Z());		
					left->SetPiecewise(Vector(p2l.X(),p1l.X()),Vector(p2l.Z(),p1l.Z()),1);		
					ref_right->SetConstant(p2r.Z());	
					right->SetPiecewise(Vector(p2r.X(),p1r.X()),Vector(p2r.Z(),p1r.Z()),1);		
				}
				else
				{
					if(p1l.X()>p1r.X())	//straight inner surface + decreasing outer radius
					{
						ref_left->SetConstant(p2l.Z());		
						left->SetPiecewise(Vector(p2l.X(),p1l.X()),Vector(p2l.Z(),p1l.Z()),1);	
						ref_right->SetPiecewise(Vector(p2r.X(),p1r.X(),p1l.X()),Vector(p2r.Z(),p2r.Z(),p2l.Z()),1);	
						right->SetPiecewise(Vector(p2r.X(),p1r.X(),p1l.X()),Vector(p2r.Z(),p1r.Z(),p1l.Z()),1);	
					}
					else		//straight inner surface + increasing outer radius
					{
						ref_left->SetPiecewise(Vector(p2l.X(),p1l.X(),p1r.X()),Vector(p2l.Z(),p2l.Z(),p2r.Z()),1);	
						left->SetPiecewise(Vector(p2l.X(),p1l.X(),p1r.X()),Vector(p2l.Z(),p1l.Z(),p1r.Z()),1);		
						ref_right->SetConstant(p2r.Z());	
						right->SetPiecewise(Vector(p2r.X(),p1r.X()),Vector(p2r.Z(),p1r.Z()),1);		
					}
				}
			}
			else		// non constant inner surface
			{
				if(p1l.X()==p1r.X())	//straight outer surface
				{
					if(p2l.X()>p2r.X()) //straight outer surface + decreasing inner radius
					{
						ref_left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X()),Vector(p2r.Z(),p2l.Z(),p2l.Z()),1);		
						left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X()),Vector(p2r.Z(),p2l.Z(),p1l.Z()),1);	
						ref_right->SetConstant(p2r.Z());
						right->SetPiecewise(Vector(p2r.X(),p1r.X()),Vector(p2r.Z(),p1r.Z()),1);
					}
					else //straight outer surface + increasing inner radius
					{
						ref_left->SetConstant(p2l.Z());	
						left->SetPiecewise(Vector(p2l.X(),p1l.X()),Vector(p2l.Z(),p1l.Z()),1);	
						ref_right->SetPiecewise(Vector(p2l.X(),p2r.X(),p1r.X()),Vector(p2l.Z(),p2r.Z(),p2r.Z()),1);
						right->SetPiecewise(Vector(p2l.X(),p2r.X(),p1r.X()),Vector(p2l.Z(),p2r.Z(),p1r.Z()),1);
					}
				}
				else	// neither outer nor inner surface straight
				{
					if(p1l.X()<p1r.X())	//increasing outer radius
					{
						if(p2l.X()<p2r.X())	//increasing inner radius
						{
							ref_left->SetPiecewise(Vector(p2l.X(),p1l.X(),p1r.X()),Vector(p2l.Z(),p2l.Z(),p2r.Z()),1);		
							left->SetPiecewise(Vector(p2l.X(),p1l.X(),p1r.X()),Vector(p2l.Z(),p1l.Z(),p1r.Z()),1);		
							ref_right->SetPiecewise(Vector(p2l.X(),p2r.X(),p1r.X()),Vector(p2l.Z(),p2r.Z(),p2r.Z()),1);
							right->SetPiecewise(Vector(p2l.X(),p2r.X(),p2r.X()),Vector(p2l.Z(),p2r.Z(),p1r.Z()),1);
						}
						else	//decreasing inner radius
						{
							ref_left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X(),p1r.X()),Vector(p2r.Z(),p2l.Z(),p2l.Z(),p2r.Z()),1);		
							left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X(),p1r.X()),Vector(p2r.Z(),p2l.Z(),p1l.Z(),p1r.Z()),1);		
							ref_right->SetConstant(p2r.Z());
							right->SetPiecewise(Vector(p2r.X(),p1r.X()),Vector(p2r.Z(),p1r.Z()),1);
						}
					}
					else //decreasing outer radius
					{
						if(p2l.X()<p2r.X())	//increasing inner radius
						{
							ref_left->SetConstant(p2l.Z());
							left->SetPiecewise(Vector(p2l.X(),p1l.X()),Vector(p2l.Z(),p1l.Z()),1);	
							ref_right->SetPiecewise(Vector(p2l.X(),p2r.X(),p1r.X(),p1l.X()),Vector(p2l.Z(),p2r.Z(),p2r.Z(),p2l.Z()),1);
							right->SetPiecewise(Vector(p2l.X(),p2r.X(),p1r.X(),p1l.X()),Vector(p2l.Z(),p2r.Z(),p1r.Z(),p1l.Z()),1);
						}
						else	//decreasing inner radius
						{
							ref_left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X()),Vector(p2r.Z(),p2l.Z(),p2l.Z()),1);		
							left->SetPiecewise(Vector(p2r.X(),p2l.X(),p1l.X()),Vector(p2r.Z(),p2l.Z(),p1l.Z()),1);	
							ref_right->SetPiecewise(Vector(p2r.X(),p1r.X(),p1l.X()),Vector(p2r.Z(),p2r.Z(),p2l.Z()),1);
							right->SetPiecewise(Vector(p2r.X(),p1r.X(),p1l.X()),Vector(p2r.Z(),p1r.Z(),p1l.Z()),1);
						}
					}
				}
			}

			//for(int i=1;i<=11;i++)
			//{
			//	double t=i/10.;
			//	
			//	GetMBS()->UO()<< "ref_left("<<t<<")="<<ref_left->Evaluate(t)<<"\n";
			//	GetMBS()->UO()<< "left("<<t<<")="<<left->Evaluate(t)<<"\n";
			//	GetMBS()->UO()<< "ref_right("<<t<<")="<<ref_right->Evaluate(t)<<"\n";
			//	GetMBS()->UO()<< "right("<<t<<")="<<right->Evaluate(t)<<"\n";
			//}

			Distort3DMesh(nodes_tmp,1,3,right,ref_right,left,ref_left,0,2);
		}

	// Create nodes
	for (int i=1; i<=nodes_tmp.Length(); i++)
	{
		//nodes_tmp(i) = nodes_tmp(i) + Vector3D(x_offset,y_offset,0);
		blocknodes.Add(AddNodeCheck(nodes_tmp(i),bodynr,1E-10));
	}

	return 0;
}
// The coordinate "deformed_coord" of nodes is changed w.r.t. MathFunctions
// The coordinate is interpolated linearly between the values of outer and inner
// ref_outer and ref_inner are the reference values (=values before distortion)

// if deformed_coord2!=0 than radial deformation is performed
// if independent_coord2!=0: new_value = f(r)*old_value, with r = sqrt(independent_coord^2+independent_coord2^2)

void FEMesh_Generator::Distort3DMesh(TArray<Vector3D>& nodes, int independent_coord, int deformed_coord, MathFunction* outer, double ref_outer, MathFunction* inner, double ref_inner, int deformed_coord2, int independent_coord2)
{
	MathFunction* MFref_outer= new MathFunction;
	MathFunction* MFref_inner= new MathFunction;
	MFref_outer->SetConstant(ref_outer);
	MFref_inner->SetConstant(ref_inner);
	Distort3DMesh(nodes, independent_coord, deformed_coord,outer, MFref_outer, inner,MFref_inner, deformed_coord2, independent_coord2);
}

void FEMesh_Generator::Distort3DMesh(TArray<Vector3D>& nodes, int independent_coord, int deformed_coord, MathFunction* outer, MathFunction* ref_outer, MathFunction* inner, MathFunction* ref_inner, int deformed_coord2, int independent_coord2)
{
	int j;
	double t,value, old_value;
	double r_out, r_in;
	Vector3D tmp;

	double fo;	// fo... final value outside, 
	double fi;	// fi... final value inside

	// deform mesh
	for (int i=1; i<=nodes.Length(); i++)
	{
		tmp=nodes(i);

		// get value of independent variable of point
		if(independent_coord2)
		{
			t=sqrt(tmp(independent_coord)*tmp(independent_coord)+tmp(independent_coord2)*tmp(independent_coord2));	//radius
		}
		else
		{
			t=tmp(independent_coord);										
		}

		// get the reference value
		if(deformed_coord2)
		{
			old_value = sqrt(tmp(deformed_coord)*tmp(deformed_coord)+tmp(deformed_coord2)*tmp(deformed_coord2));
		}
		else
		{
			old_value = tmp(deformed_coord);
		}

		// get new value
		fo = 0;
		fi = 0;

		if(outer != NULL)	{ fo = outer->Evaluate(t);}		// get the value of the outer MathFunction
		if(inner != NULL)	{ fi = inner->Evaluate(t);}		// get the value of the inner MathFunction

		r_out = ref_outer->Evaluate(t);			// get the reference value of the outer MathFunction
		r_in = ref_inner->Evaluate(t);			// get the reference value of the inner MathFunction

		if(r_out!=r_in)
		{
			value = fi + (fo-fi)/(r_out-r_in)*(old_value-r_in);	//interpolation
		}
		else
		{
			value = fo;
		}

		if(deformed_coord2)
		{
			if(old_value!=0)	// otherwise coordinates remain unchanged for radial shearing
			{
				double factor = value/old_value;
				tmp(deformed_coord)	=factor*tmp(deformed_coord);
				tmp(deformed_coord2)=factor*tmp(deformed_coord2);
			}
		}
		else
		{
			tmp(deformed_coord) = value;
		}
		nodes(i)=tmp;
	}
}

// creates nodes for cylinder - does not add nodes to mesh.points - just to TArrays
int FEMesh_Generator::CreateCylinder_Mesh2D(double radius, double small_radius, int divisions_angle, TArray<double>& x, TArray<double>& y)
{
	// init
	int		d = divisions_angle;
	double	r = radius;
	double	r1= small_radius;
	double	l = r1/d;     // size of a small square
	double sqrt2 = sqrt((double)2);
	double	r0=	r1*sqrt2; // diameter of the square-area
	int		n  = d+1;			
	int		nN = n*n+(d+n)*d;	// number of nodes in a 1/4 circle
	
	double x_tmp, y_tmp;

	// xy-coordinates of the squares
	for (int i=0; i<=d; i++)
	{
		y_tmp = i*l;
		for (int j=0; j<=d; j++)
		{
			x_tmp=j*l;
			x.Add(x_tmp/*+x_offset*/);
			y.Add(y_tmp/*+y_offset*/);
		}
	}

	double r0_tmp, r1_tmp, alpha_max, r2_tmp, x0, alpha_tmp;
	double tmp1,tmp2;
	int j;
	int m=1;

	// xy-coordinates of the nodes on the circles
	for (int i=1; i<=d; i++)
	{
		r0_tmp=r0+((r-r0)/d)*i;
		r1_tmp=r1+((r-r1)/d)*i;
		tmp1 =-r0_tmp * (sqrt2 * r1_tmp - r0_tmp) / (-r0_tmp*r0_tmp + r0_tmp * sqrt2 * r1_tmp - r1_tmp*r1_tmp);
		tmp2=-(r0_tmp * sqrt2 - r1_tmp) * r1_tmp / (-r0_tmp*r0_tmp + r0_tmp * sqrt2 * r1_tmp - r1_tmp*r1_tmp);
		// alpha_max=atan2(-r0_tmp * (sqrt(0.2e1) * r1_tmp - r0_tmp) / (-r0_tmp ^ 2 + r0_tmp * sqrt(0.2e1) * r1_tmp - r1_tmp ^ 2), -(r0_tmp * sqrt(0.2e1) - r1_tmp) * r1_tmp / (-r0_tmp ^ 2 + r0_tmp * sqrt(0.2e1) * r1_tmp - r1_tmp ^ 2));
		alpha_max=atan2(tmp1,tmp2);
		r2_tmp=r0_tmp*sin(MY_PI/4)/sin(alpha_max);
		x0=r2_tmp-r1_tmp;
		for (j=0; j<=d; j++)
		{
			alpha_tmp=j*alpha_max/d;
			x_tmp=r2_tmp*cos(alpha_tmp)-x0;
			y_tmp=r2_tmp*sin(alpha_tmp);
			x.Add(x_tmp/*+x_offset*/);
			y.Add(y_tmp/*+y_offset*/);
			m++;
		}
		j=d;
		while (j>0)
		{
			j=j-1;
			alpha_tmp=j*alpha_max/d;
			x_tmp=r2_tmp*sin(alpha_tmp);
			y_tmp=r2_tmp*cos(alpha_tmp)-x0;
			x.Add(x_tmp/*+x_offset*/);
			y.Add(y_tmp/*+y_offset*/);
		}
	}

	return x.Length()-nN;		// should be 0
}

// creates nodes for an hollow cylinder - does not add nodes to mesh.points - just to TArrays
int FEMesh_Generator::CreateHollowCylinder_Mesh2D(double radius, double radius_hole, int divisions_angle, TArray<double>& x, TArray<double>& y)
{
	// init
	int	d  = divisions_angle;
	double	r2 = radius;
	double	r1= radius_hole;

	//double sqrt2 = sqrt((double)2);
	//int		n  = d+1;			
	int		nN = 2*(2*d+1);	// number of nodes in a 1/4 circle

	double x_tmp, y_tmp, r_tmp, phi_tmp;

	// xy-coordinates of the hole (circle 1)
	r_tmp=r1;
	for (int j=0; j<=2*d; j++)
	{
		phi_tmp=j*MY_PI/(4*d);
		x_tmp=r_tmp*cos(phi_tmp);
		y_tmp=r_tmp*sin(phi_tmp);
		x.Add(x_tmp);
		y.Add(y_tmp);
	}

	// xy-coordinates of the surface (circle 3)
	r_tmp=r2;
	for (int j=0; j<=2*d; j++)
	{
		phi_tmp=j*MY_PI/(4*d);
		x_tmp=r_tmp*cos(phi_tmp);
		y_tmp=r_tmp*sin(phi_tmp);
		x.Add(x_tmp);
		y.Add(y_tmp);
	}

	return x.Length()-nN;		// should be 0
}

// creates nodes for an hollow cylinder with increasing number of divisions/rad - does not add nodes to mesh.points - just to TArrays
int FEMesh_Generator::CreateHollowCylinderDoubleDiv_Mesh2D(double radius, double radius_hole, int divisions_angle, TArray<double>& x, TArray<double>& y)
{
	// init
	int		d  = divisions_angle;
	double	r3 = radius;	// radius of the surface
	double	r1 = radius_hole;
	double	r2 = (r3-r1)/2+r1;
	int		f  = 5;			// factor coresp. to the distance of the nodes at the connection-circle

	int		n  = 2*d+1;	// number of nodes for 1/4 circle 1
	int		n2 = 3*d;   // number of nodes for connection circle
	int		nN = n+n2+4*d+1;	// number of nodes in a 1/4 circle

	double x_tmp, y_tmp, r_tmp, phi_tmp, d_phi, phi_0;

	// xy-coordinates of the hole (circle 1)
	r_tmp=r1;
	for (int j=0; j<=2*d; j++)
	{
		phi_tmp=j*MY_PI/(4*d);
		x_tmp=r_tmp*cos(phi_tmp);
		y_tmp=r_tmp*sin(phi_tmp);
		x.Add(x_tmp/*+x_offset*/);
		y.Add(y_tmp/*+y_offset*/);
	}

	// xy-coordinates of the connection ring (circle 2)

	r_tmp = r2;
	d_phi = MY_PI/(2*d*f);
	phi_0 = MY_PI/(4*d);
	for (int j=1; j<=d; j++)
	{
		for (int i=-1; i<=1; i++)
		{
			phi_tmp=phi_0+i*d_phi;
			x_tmp=r_tmp*cos(phi_tmp);
			y_tmp=r_tmp*sin(phi_tmp);
			x.Add(x_tmp/*+x_offset*/);
			y.Add(y_tmp/*+y_offset*/);
		}
		phi_0=phi_0+MY_PI/(2*d);
	}

	// xy-coordinates of the surface (circle 3)
	r_tmp=r3;
	for (int j=0; j<=4*d; j++)
	{
		phi_tmp=j*MY_PI/(8*d);
		x_tmp=r_tmp*cos(phi_tmp);
		y_tmp=r_tmp*sin(phi_tmp);
		x.Add(x_tmp/*+x_offset*/);
		y.Add(y_tmp/*+y_offset*/);
	}

	return x.Length()-nN;		// should be 0
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create elements

// creates Hexaherals for cylinder from nodelist - adds FEHex to mesh.elements 
int FEMesh_Generator::CreateCylinder_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color)
{
	blockelements.Flush();

	int k=0; int j=0; 
	int		d = divisions_angle;
	int		dh= divisions_height;
	int		n = d+1;	
	TArray<int4> E;	// Elements in xy-plane
	int4 el2D;		// node numbers of the element in xy-plane
	
	// Elements of the square
	for (j=1; j<=d;j++)
	{
		for (int i=1; i<=d;i++)
		{
			k=i+(j-1)*n;
			//E(m,:)=[m k k+1 (k+1)+n k+n];
			el2D(1) = k;
			el2D(2) = k+1;
			el2D(3) = k+1+n;
			el2D(4) = k+n;
			E.Add(el2D);
		}
	}

	// Elements of the connection area
	k=0;
	for (j=1; j<=d;j++)
	{
		//E(m,:)=[m j*n n^2+j n^2+j+1 (j+1)*n];
		el2D(1) = j*n;
		el2D(2) = n*n+j;
		el2D(3) = n*n+j+1;
		el2D(4) = (j+1)*n;
		E.Add(el2D);
	}
	for (j=0; j<=d-1;j++)
	{
		//E(m,:)=[m n^2-j n*(n+1)+j n*(n+1)+j+1 n^2-j-1];
		el2D(1) = n*n-j;
		el2D(2) = n*(n+1)+j;
		el2D(3) = n*(n+1)+j+1;
		el2D(4) = n*n-j-1;
		E.Add(el2D);
	}

	// Elements of the circles
	k=n*n;
	for (int i=1; i<=d-1;i++)
	{
		for (j=1; j<=2*d;j++)
		{
			//E(m,:)=[m k+j k+2*n+j-1 k+2*n+j k+j+1];
			el2D(1) = k+j;
			el2D(2) = k+2*n+j-1;
			el2D(3) = k+2*n+j;
			el2D(4) = k+j+1;
			E.Add(el2D);
		}
		k=k+2*d+1;
	}

	int numberElements=E.Length();

	// Elements of the cylinder --> creating 3D-elements
	int nN = n*n+(d+n)*d;	// number of nodes in a 1/4 circle
	IVector points(8);
	for (int i=1; i<=dh;i++)
	{
		for (j=1; j<=numberElements;j++)
		{
			//V(m,:)=[m E(j,2)+(i-1)*nN E(j,3)+(i-1)*nN E(j,4)+(i-1)*nN E(j,5)+(i-1)*nN E(j,2)+i*nN E(j,3)+i*nN E(j,4)+i*nN E(j,5)+i*nN];
			points(1) = blocknodes(E(j)(1)+(i-1)*nN	);
			points(2) = blocknodes(E(j)(2)+(i-1)*nN	);
			points(4) = blocknodes(E(j)(3)+(i-1)*nN	);
			points(3) = blocknodes(E(j)(4)+(i-1)*nN	);
			points(5) = blocknodes(E(j)(1)+i*nN		);
			points(6) = blocknodes(E(j)(2)+i*nN		);
			points(8) = blocknodes(E(j)(3)+i*nN		);
			points(7) = blocknodes(E(j)(4)+i*nN		);


			FEHex fehex(points);
			fehex.MaterialNum() = matnr;
			fehex.Domain() = bodynr;
			fehex.Color() = color;
			blockelements.Add(AddElement(fehex));
		}
	}

	return blockelements.Length();
}

// creates Hexaherals for an hollow cylinder from nodelist - adds FEHex to mesh.elements 
int FEMesh_Generator::CreateHollowCylinder_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color)
{
	blockelements.Flush();

	int k=0; int j=0; 
	int		d = divisions_angle;
	int		dh= divisions_height;
	int		n = d+1;	
	TArray<int4> E;	// Elements in xy-plane
	int4 el2D;		// node numbers of the element in xy-plane
	

	// Elements of the circles
	k=0;

	for (j=1; j<=2*d;j++)
	{
		//E(m,:)=[m k+j k+2*n+j-1 k+2*n+j k+j+1];
		el2D(1) = k+j;
		el2D(2) = k+2*n+j-1;
		el2D(3) = k+2*n+j;
		el2D(4) = k+j+1;
		E.Add(el2D);
	}

	int numberElements=E.Length();

	// Elements of the cylinder --> creating 3D-elements
	int nN = 2*(2*d+1);	// number of nodes in a 1/4 circle
	IVector points(8);
	for (int i=1; i<=dh;i++)
	{
		for (j=1; j<=numberElements;j++)
		{
			//V(m,:)=[m E(j,2)+(i-1)*nN E(j,3)+(i-1)*nN E(j,4)+(i-1)*nN E(j,5)+(i-1)*nN E(j,2)+i*nN E(j,3)+i*nN E(j,4)+i*nN E(j,5)+i*nN];
			points(1) = blocknodes( E(j)(1)+(i-1)*nN	);
			points(2) = blocknodes( E(j)(2)+(i-1)*nN	);
			points(4) = blocknodes( E(j)(3)+(i-1)*nN	);
			points(3) = blocknodes( E(j)(4)+(i-1)*nN	);
			points(5) = blocknodes( E(j)(1)+i*nN		);
			points(6) = blocknodes( E(j)(2)+i*nN		);
			points(8) = blocknodes( E(j)(3)+i*nN		);
			points(7) = blocknodes( E(j)(4)+i*nN		);

			FEHex fehex(points);
			fehex.MaterialNum() = matnr;
			fehex.Domain() = bodynr;
			fehex.Color() = color;
			blockelements.Add(AddElement(fehex));
		}
	}

	return blockelements.Length();
}

// creates Hexaherals for an hollow cylinder with increasing number of divisions/rad from nodelist - adds FEHex to mesh.elements 
int FEMesh_Generator::CreateHollowCylinderDoubleDiv_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color)
{
	blockelements.Flush();

	int k=0; int j=0; 
	int		d = divisions_angle;
	int		dh= divisions_height;
	int		n = 2*d+1;	// number of nodes on the 1/4 circle 1
	int		n2 = 3*d;	// number of nodes for 1/4 connection circle


	TArray<int4> E;	// Elements in xy-plane
	int4 el2D;		// node numbers of the element in xy-plane
	

	// Elements of ring 1
	int k2=n;
	for (int i=1; i<=2*d;i++)
	{
		k=i;   
		//E(m,:)=[m k k+k2 k+k2+1 k+1];
		el2D(1) = k;
		el2D(2) = k+k2;
		el2D(3) = k+k2+1;
		el2D(4) = k+1;
		E.Add(el2D);
		if((i+1)%2) k2=k2+1;
	}

	// Elements of ring 2
	k=n;
	for (int i=1; i<=2*d;i++)
	{
		k=k+1;   
		//E(m,:)=[m k k+k2 k+k2+1 k+1];
		el2D(1) = k;
		el2D(2) = k+k2;
		el2D(3) = k+k2+1;
		el2D(4) = k+1;
		E.Add(el2D);
		if((i+1)%2)
		{
			k2=k2+1;
			k=k+1;
		}
	}

	// Elements of connection 1
	k = -1;
	int k1 = n-2;
	k2 = n+n2-3;
	for (int i=1; i<=d;i++)
	{
		k=k+2;
		k1=k1+3;
		k2=k2+4;
		//E(m,:)=[m k k2 k2+1 k1];
		el2D(1) = k;
		el2D(2) = k2;
		el2D(3) = k2+1;
		el2D(4) = k1;
		E.Add(el2D);
	}

	// Elements of connection 2
	k = 1;
	k1 = n;
	k2 = n+n2;
	for (int i=1; i<=d;i++)
	{
		k=k+2;
		k1=k1+3;
		k2=k2+4;
		//E(m,:)=[m k k1 k2 k2+1];
		el2D(1) = k;
		el2D(2) = k1;
		el2D(3) = k2;
		el2D(4) = k2+1;
		E.Add(el2D);
	}

	int numberElements=E.Length();

	// Elements of the cylinder --> creating 3D-elements
	int		nN = n+n2+4*d+1;	// number of nodes in a 1/4 circle
	IVector points(8);
	for (int i=1; i<=dh;i++)
	{
		for (j=1; j<=numberElements;j++)
		{
			//V(m,:)=[m E(j,2)+(i-1)*nN E(j,3)+(i-1)*nN E(j,4)+(i-1)*nN E(j,5)+(i-1)*nN E(j,2)+i*nN E(j,3)+i*nN E(j,4)+i*nN E(j,5)+i*nN];
			points(1) = blocknodes( E(j)(1)+(i-1)*nN	);
			points(2) = blocknodes( E(j)(2)+(i-1)*nN	);
			points(4) = blocknodes( E(j)(3)+(i-1)*nN	);
			points(3) = blocknodes( E(j)(4)+(i-1)*nN	);
			points(5) = blocknodes( E(j)(1)+i*nN		);
			points(6) = blocknodes( E(j)(2)+i*nN		);
			points(8) = blocknodes( E(j)(3)+i*nN		);
			points(7) = blocknodes( E(j)(4)+i*nN		);

			FEHex fehex(points);
			fehex.MaterialNum() = matnr;
			fehex.Domain() = bodynr;
			fehex.Color() = color;
			blockelements.Add(AddElement(fehex));
		}
	}

	return blockelements.Length();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Quadrilateral Blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// abbreviated calls

// generates a block of quadrilateral elements (no return values, no boundary conditions)
int FEMesh_Generator::GenerateQuadrilateralBlock(Box2D block, double thickness, int2 divisions_xy, int bodynr, int matnr, Vector3D color)
{
	IVector blocknodes;
	IVector blockelements;
	int order = 1;
  
	CreateQuadBlock_Nodes(blocknodes, block, divisions_xy, bodynr, order);
	CreateQuadBlock_Quads(blockelements, blocknodes, thickness, divisions_xy, bodynr, matnr, color, order );

	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create nodes

// creates nodes for quadrilateral block - adds nodes to mesh.points 
int FEMesh_Generator::CreateQuadBlock_Nodes(IVector& blocknodes, Box2D block, int2 divisions_xy, int bodynr, int order)
{
	// resize nodestree
	ComputeNodeBox();
	Box3D newbox(GetNodeBox());
	if(! (newbox.IsIn(block.PMax().MakeV3D()) && newbox.IsIn(block.PMin().MakeV3D())) )
	{
		newbox.Add(block.PMax().MakeV3D());
		newbox.Add(block.PMin().MakeV3D());
		int additional_nodes = (divisions_xy(1)+1)*(divisions_xy(2)+1);
		int divs = (int) pow( (NN()+additional_nodes)*0.1 ,1/3.) +1;
		NodesTree().ResetSearchTree(divs, divs, divs, newbox);
		for(int i=1; i<=NN(); i++) 
			NodesTree().AddItem(Box3D(GetPoint3D(i),1e-8),i);
	}

	blocknodes.Flush();
	// create all nodes for the quadrilateral with 'corners', do not make any double nodes, save nodenumbers of the nodes of the block in IVector
	for (int j=0; j <= divisions_xy(2); j++)
	{
		double fy = (double)(j)/(double)divisions_xy(2);
		for (int i=0; i <= divisions_xy(1); i++)
		{
			double fx = (double)(i)/(double)divisions_xy(1);
			Vector2D pos = block.PMin() + Vector2D( block.SizeX()*fx, block.SizeY()*fy );

			blocknodes.Add(AddNodeCheck(pos,bodynr,1E-10));		
		}
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functions to create elements

// returns Nodelist of a single (quadrilateral) element of the block
int FEMesh_Generator::CreateQuadBlock_QuadNodes(IVector& points, IVector& blocknodes, int firstnode, int dx, int dy, int order)
{
	// node numbering in HotInt sequence:
	if(order == 1)
	{
		points.SetLen(8);
		points(1) = blocknodes( firstnode              );
		points(2) = blocknodes( firstnode          + 1 );
		points(3) = blocknodes( firstnode + (dx+1) + 1 );
		points(4) = blocknodes( firstnode + (dx+1)     );
		return 4;
	}
	return points.Length();
}

// creates Quadrilaterals from nodelist - adds FEQuadx to mesh.elements
int FEMesh_Generator::CreateQuadBlock_Quads(IVector& blockelements, IVector& blocknodes, double thickness, int2 divisions_xy, int bodynr, int matnr, Vector3D color, int order)
{
	blockelements.Flush();
	IVector points(4);
	
	int firstnode;
	int dx=divisions_xy(1);
	int dy=divisions_xy(2);

	for (int iy=1; iy <= dy; iy++)
	{
		for (int ix=1; ix <= dx; ix++)
		{
			firstnode = (iy-1)*(dx+1) + (ix-1)*1 + 1;
			CreateQuadBlock_QuadNodes(points,blocknodes,firstnode,dx,dy,order);
			// 1 quad
			FEQuad fequad(points);
			fequad.Thickness() = thickness;
			fequad.MaterialNum() = matnr;
			fequad.Domain() = bodynr;
			fequad.Color() = color;
			blockelements.Add(AddElement(fequad));
		}
	}
	return blockelements.Length();
}








//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//TODO: kill off all these functions... outerfaces etc should be handeld in class SegmentedBlock
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// create outer faces - adds FEQuad to mesh.faces
int FEMesh_Generator::CreateHexBlock_Hexes_OneFace(IVector& outerfaces, IVector& subset_elements, int side)
{
	outerfaces.Flush();

	for(int i=1; i <= subset_elements.Length(); i++)
	{
		int elnr = subset_elements(i);
		int4 tempnodes = GetElement(elnr).GetSideNodeNum4(side);
		FEMesh_Face tempface(tempnodes, elnr, side);
		outerfaces.Add(Faces().Add(tempface));
	}
	return outerfaces.Length();
}

// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
int FEMesh_Generator::CreateHexBlock_Hexes_AllFaces( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Hexes_ComputeOuterElements(outerelems,divisions_xyz);
	
	for(int i=1; i <= 6; i++)
	{
		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements(outerelems(i)->Get(j));
			int4 tempnodes = GetElement(elnr).GetSideNodeNum4(i);
			FEMesh_Face tempface(tempnodes, elnr, i);
			outerfaces(i)->Add(Faces().Add(tempface));
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
int FEMesh_Generator::CreateHexBlock_Prisms_AllFaces( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Prisms_ComputeOuterElements(outerelems,divisions_xyz);
	
	for(int i=1; i <= 6; i++)
	{
		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements( outerelems(i)->Get(j) );
			int sidenr = i;
			int4 tempnodes;

			if( i<5 )   // quads for sides 1,3: first prism,; hex sides 1,3 == prism1 sides 1,3
			{           // quads for sides 2,4: second prism; hex sides 2,4 == prism2 sides 1,3
				if ((i==2)||(i==4)) sidenr--;

				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				FEMesh_Face tempface(tempnodes, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
      if ((i==5)||(i==6)) // trigs for sides 5,6: both prisms; hex sides 5,6 == prism1 sides 5,6 {4,1,2} + prism2 sides 5,6 {4,1,2}
			{
				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr); // element is still hex
				int3 pts3(tempnodes(4),tempnodes(1),tempnodes(2));
				FEMesh_Face tempface(pts3, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

// $EK 2013-03-04
// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
int FEMesh_Generator::CreateHexBlock_Prisms_AllFaces_new( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Prisms_ComputeOuterElements(outerelems,divisions_xyz);
	
	for(int i=1; i <= 6; i++)
	{
		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements( outerelems(i)->Get(j) );
			int sidenr;
			switch (i)
			{
				case 1: case 2:
					sidenr = 2; break;
				case 3: case 4: 
					sidenr = 1; break;
				case 5:  
					sidenr = 4; break;
				case 6: 
					sidenr = 5; break;
			}
			int4 tempnodes;
			//order is -x, +x, -y, +y, -z, +z
			if (sidenr < 4)   // quads for sides 1,3: first prism,; hex sides 1,3 == prism1 sides 1,3
			{           // quads for sides 2,4: second prism; hex sides 2,4 == prism2 sides 1,3
				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				FEMesh_Face tempface(tempnodes, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
			else // trigs for sides 5,6: 
			{
				int3 pts3 = GetElement(elnr).GetSideNodeNum3(sidenr); // element is still hex
				//int3 pts3(tempnodes(4),tempnodes(1),tempnodes(2));
				FEMesh_Face tempface(pts3, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces 
int FEMesh_Generator::CreateHexBlock_Pyrams_AllFaces( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Pyrams_ComputeOuterElements(outerelems,divisions_xyz);

	for(int i=1; i <= 6; i++)
	{
		int face_has_trigs; // 1 if face of hexahedral#1 has 2 trigs: faces(1,3,6), 0 if face is quad;
		switch(i)
		{
			case 1: case 3: case 6: face_has_trigs = 1; break;     
			case 2: case 4: case 5: face_has_trigs = 0; break;
			default: break;
		}

		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements( outerelems(i)->Get(j) ); // ("global)" element number within the mesh !
			int sidenr;
			int4 tempnodes;

//			int hexnr = (elemnr-1)/3 +1;
			int hexnr = (outerelems(i)->Get(j)-1)/3 +1;
			int ix = (hexnr-1)%divisions_xyz(1) +1;
			int iy = ((hexnr-1-(ix-1))/divisions_xyz(1))%divisions_xyz(2) +1;
			int iz = ((hexnr-1-(ix-1)-(divisions_xyz(1)*(iy-1))))/(divisions_xyz(1)*divisions_xyz(2)) +1;

			int swapped = (ix + iy + iz -1)%2;
			int trigs = (face_has_trigs + swapped)%2; 
	
			if( trigs == 0) // quad
			{
				switch(i)
				{
					case 1: case 2: sidenr = 5; break; // pyram 2 side 5
					case 3: case 4: sidenr = 5; break; // pyram 3 side 5
					case 5: case 6: sidenr = 5; break; // pyram 1 side 5
					default: break;
				}
				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				FEMesh_Face tempface(tempnodes, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
			else // trigs
			{
				switch(i)
				{
					case 1: case 2: if((elnr%3 == 1)) sidenr = 1; else sidenr = 3; break; // T1: pyram1 side 1 or T2: pyram3 side 3
					case 3: case 4: if((elnr%3 == 1)) sidenr = 3; else sidenr = 1; break; // T1: pyram1 side 3 or T2: pyram2 side 1
					case 5: case 6: if((elnr%3 == 2)) sidenr = 3; else sidenr = 1; break; // T1: pyram2 side 3 or T2: pyram3 side 1
					default: break;
				}
				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				int3 pts3(tempnodes(4),tempnodes(1),tempnodes(2));
				FEMesh_Face tempface(pts3, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

//$EK 2013-03-04 adapted for real pyramids
// creates all faces of block - 6 individual lists - adds FEQuad, FETrig to mesh.faces   
int FEMesh_Generator::CreateHexBlock_Pyramids_AllFaces_new( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Pyrams_ComputeOuterElements(outerelems,divisions_xyz);

	for(int i=1; i <= 6; i++)
	{
		int face_has_trigs; // 1 if face of hexahedral#1 has 2 trigs: faces(1,3,6), 0 if face is quad;
		//1 -x, 2 +x, 3 -y, 4 +y, 5 -z, 6 +z
		switch(i)
		{
			case 1: case 3: case 6: face_has_trigs = 1; break;     
			case 2: case 4: case 5: face_has_trigs = 0; break;
			default: break;
		}

		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements( outerelems(i)->Get(j) ); // ("global)" element number within the mesh !
			int elnr_rel = elnr - GetMesh()->NElements() + 3*(divisions_xyz(1)*divisions_xyz(2)*divisions_xyz(3));
			int sidenr;
			int4 tempnodes;

//			int hexnr = (elemnr-1)/3 +1;
			int hexnr = (outerelems(i)->Get(j)-1)/3 +1;
			int ix = (hexnr-1)%divisions_xyz(1) +1;
			int iy = ((hexnr-1-(ix-1))/divisions_xyz(1))%divisions_xyz(2) +1;
			int iz = ((hexnr-1-(ix-1)-(divisions_xyz(1)*(iy-1))))/(divisions_xyz(1)*divisions_xyz(2)) +1;

			int swapped = (ix + iy + iz -1)%2;
			int trigs = (face_has_trigs + swapped)%2; 
	
			if( trigs == 0) // quad
			{
				sidenr = 1;
				tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				FEMesh_Face tempface(tempnodes, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
			else // trigs
			{
				switch(i)
				{
					case 1: case 2: if((elnr_rel%3 == 1)) sidenr = 5; else sidenr = 2; break; // T1: pyram1 side 1 or T2: pyram3 side 3
					case 3: case 4: if((elnr_rel%3 == 1)) sidenr = 2; else sidenr = 5; break; // T1: pyram1 side 3 or T2: pyram2 side 1
					case 5: case 6: if((elnr_rel%3 == 2)) sidenr = 2; else sidenr = 5; break; // T1: pyram2 side 3 or T2: pyram3 side 1
					default: break;
				}

				//tempnodes = GetElement(elnr).GetSideNodeNum4(sidenr);
				//int3 pts3(tempnodes(4),tempnodes(1),tempnodes(2));
				int3 pts3 = GetElement(elnr).GetSideNodeNum3(sidenr);
				FEMesh_Face tempface(pts3, elnr, sidenr);
				outerfaces(i)->Add(Faces().Add(tempface));
			}
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
int FEMesh_Generator::CreateHexBlock_Tetras_AllFaces( TArray<IVector*>& outerfaces, IVector& blockelements, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerfaces);
	outerfaces.SetLen(6);
	for(int i=1; i <= 6; i++) outerfaces(i) = new IVector();
	TArray<IVector*> outerelems;

	CreateHexBlock_Tetras_ComputeOuterElements(outerelems,divisions_xyz);

	for(int i=1; i <= 6; i++)
	{
		for(int j=1; j <= outerelems(i)->Length(); j++)
		{
			int elnr = blockelements( outerelems(i)->Get(j) );
			int sidenr;
			int3 tempnodes3;
  
			int trignr = (j -1)%2 +1;
			int hexnr = (outerelems(i)->Get(j)-1)/5 +1;
			int ix = (hexnr-1)%divisions_xyz(1) +1;
			int iy = ((hexnr-1-(ix-1))/divisions_xyz(1))%divisions_xyz(2) +1;
			int iz = ((hexnr-1-(ix-1)-(divisions_xyz(1)*(iy-1))))/(divisions_xyz(1)*divisions_xyz(2)) +1;

			int swapped = (ix +iy +iz - 1) % 2; 

			int ibar = i;
			if(swapped)
			{
				if(i%2) ibar = i+1; // 1->2, ...
				else ibar = i-1;    // 2->1, ... 
			}

			switch(ibar)
			{
			case 1: if ((j%2) == 1) sidenr = 4; else sidenr = 2; break; // T1: tet1, side4 and T2: tet4, side2
			case 2: sidenr = 4; break; // T1: tet2, side4 and T2: tet3, side4
			case 3: sidenr = 2; break; // T1: tet1, side2 and T2: tet3, side2
			case 4: if ((j%2) == 1) sidenr = 2; else sidenr = 4; break; // T1: tet2, side2 and T2: tet4, side4
			case 5: sidenr = 1; break; // T1: tet1, side1 and T2: tet2, side1
			case 6: sidenr = 1; break; // T1: tet3, side1 and T2: tet4, side1
			default: break;
			}

			tempnodes3 = GetElement(elnr).GetSideNodeNum3(sidenr);
			FEMesh_Face tempface(tempnodes3, elnr, sidenr);
			outerfaces(i)->Add(Faces().Add(tempface));
		}
	}
	ReleaseArray_TemplatePtr(outerelems);

	int totalfaces=0;
	for(int i=1; i <= 6; i++) totalfaces += outerfaces(i)->Length();
	return totalfaces;
}

// constraints for one blockface - not implemented yet
int FEMesh_Generator::CreateHexBlock_OneFaceConstraints(int type, int direction, IVector& subset_outerfaces)
{
// this does not work yet, add more loadtypes in mesh
	return 0;
}

// constraints for all blockfaces - not implemented yet
int FEMesh_Generator::CreateHexBlock_AllFaceConstraints(IVector& types, IVector& directions, TArray<IVector*>& blockfaces)
{
  for(int i=1; i<=6; i++) CreateHexBlock_OneFaceConstraints(types(i), directions(i), *blockfaces(i));
	return 0;
}

// apply a bodyload on all elements of the block
int FEMesh_Generator::CreateHexBlock_BodyLoads(MBSLoad& bodyload, IVector& blockelements)
{
	int loads_before = NLoads();
	
	int gravity_flag;
	double val;
	int dir;
 
	if (bodyload.LoadType() == Tgravity) { gravity_flag = 1; bodyload.GetGravity(val,dir);}
	else { gravity_flag = 0; bodyload.GetBodyLoad(val,dir); }

	for(int i=1; i <= blockelements.Length(); i++)
	{
		Vector3D forcevector(0.);
		if (dir == 1) forcevector.X() = val;
		if (dir == 2) forcevector.Y() = val;
		if (dir == 3) forcevector.Z() = val;
		AddLoad(FEMesh_BodyLoad(blockelements(i), forcevector, gravity_flag));
	}
	return NLoads()-loads_before;
}

// apply face-normal area loads to outer faces
int FEMesh_Generator::CreateHexBlock_LoadFaces(IVector& loadfaces_active, Vector& loadfaces_loadstrength, TArray<IVector*>& blockfaces, Box3D block)
{
	int loads_before = NLoads();
	
	for(int i=1; i<=6; i++)
	{
		if(loadfaces_active(i))
		{
// compute Load per unitarea of Blockface
			Vector3D lv;	
			if((i==1) || (i==2)) lv = Vector3D(loadfaces_loadstrength(i) / (block.SizeY()*block.SizeZ()), 0.0, 0.0);
			if((i==3) || (i==4)) lv = Vector3D(0.0, loadfaces_loadstrength(i) / (block.SizeX()*block.SizeZ()), 0.0);
			if((i==5) || (i==6)) lv = Vector3D(0.0, 0.0, loadfaces_loadstrength(i) / (block.SizeX()*block.SizeY()));
			
			if (i%2 == 0) lv *= -1;

			for(int j=1; j<=blockfaces(i)->Length(); j++)
			{
				//AddLoad(TFaceLoad,blockfaces(i)->Get(j),lv);
				AddLoad(FEMesh_FaceLoad(blockfaces(i)->Get(j),lv));
			}
		}
	}
	return NLoads()-loads_before;
}

// find all outer elements - 6 individual lists
int FEMesh_Generator::CreateHexBlock_Hexes_ComputeOuterElements(TArray<IVector*>& outerelems, int3 divisions_xyz)
{	
	ReleaseArray_TemplatePtr(outerelems);
	outerelems.SetLen(6);
	for(int i=1; i <= 6; i++) outerelems(i) = new IVector();

	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);
	 
 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				int elnr = (iz-1)*dy*dx + (iy-1)*dx + (ix-1)*1 + 1; // (block) element number of hex
				if(ix ==  1) outerelems(1)->Add(elnr);
				if(ix == dx) outerelems(2)->Add(elnr);
				if(iy ==  1) outerelems(3)->Add(elnr);
				if(iy == dy) outerelems(4)->Add(elnr);
				if(iz ==  1) outerelems(5)->Add(elnr);
				if(iz == dz) outerelems(6)->Add(elnr);
			}
		}
	}

	IVector allouter(0);
	for(int i=1; i <= 6; i++) allouter = allouter.Union(*outerelems(i));
	return allouter.Length();
}

// find all outer elements - 6 individual lists
int FEMesh_Generator::CreateHexBlock_Prisms_ComputeOuterElements(TArray<IVector*>& outerelems, int3 divisions_xyz)
{	
	ReleaseArray_TemplatePtr(outerelems);
	outerelems.SetLen(6);
	for(int i=1; i <= 6; i++) outerelems(i) = new IVector();

	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);
	 
 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				int elnr = 2*((iz-1)*dy*dx + (iy-1)*dx + (ix-1)*1) + 1; // (block) element number of prism with node#1 (== xmin & ymin )
				if(ix ==  1) outerelems(1)->Add(elnr);
				if(ix == dx) outerelems(2)->Add(elnr+1);
				if(iy ==  1) outerelems(3)->Add(elnr);
				if(iy == dy) outerelems(4)->Add(elnr+1);
				if(iz ==  1) { outerelems(5)->Add(elnr); outerelems(5)->Add(elnr+1); }
				if(iz == dz) { outerelems(6)->Add(elnr); outerelems(6)->Add(elnr+1); }
			}
		}
	}

	IVector allouter(0);
	for(int i=1; i <= 6; i++) allouter = allouter.Union(*outerelems(i));
	return allouter.Length();
}

// find all outer elements - 6 individual lists
int FEMesh_Generator::CreateHexBlock_Pyrams_ComputeOuterElements(TArray<IVector*>& outerelems, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerelems);
	outerelems.SetLen(6);
	for(int i=1; i <= 6; i++) outerelems(i) = new IVector();

	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);
	 
 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				int elnr = 3*((iz-1)*dy*dx + (iy-1)*dx + (ix-1)*1) + 1; // (block) element number of pyram with node#1 (== xmin & ymin )
				if( (ix+iy+iz)%2 == 1)
				{
					if(ix ==  1) { outerelems(1)->Add(elnr); outerelems(1)->Add(elnr+2); }
					if(ix == dx) outerelems(2)->Add(elnr+1);
					if(iy ==  1) { outerelems(3)->Add(elnr); outerelems(3)->Add(elnr+1); }
					if(iy == dy) outerelems(4)->Add(elnr+2);
					if(iz ==  1) outerelems(5)->Add(elnr);
					if(iz == dz) { outerelems(6)->Add(elnr+1); outerelems(6)->Add(elnr+2); }
				}
				else 
				{
					if(ix == dx) { outerelems(2)->Add(elnr); outerelems(2)->Add(elnr+2); }
					if(ix ==  1) outerelems(1)->Add(elnr+1);
					if(iy == dy) { outerelems(4)->Add(elnr); outerelems(4)->Add(elnr+1); }
					if(iy ==  1) outerelems(3)->Add(elnr+2);
					if(iz == dz) outerelems(6)->Add(elnr);
					if(iz ==  1) { outerelems(5)->Add(elnr+1); outerelems(5)->Add(elnr+2); }
				}
			}
		}
	}

	IVector allouter(0);
	for(int i=1; i <= 6; i++) allouter = allouter.Union(*outerelems(i));
	return allouter.Length();
}

// find all outer elements - 6 individual lists
int FEMesh_Generator::CreateHexBlock_Tetras_ComputeOuterElements(TArray<IVector*>& outerelems, int3 divisions_xyz)
{
	ReleaseArray_TemplatePtr(outerelems);
	outerelems.SetLen(6);
	for(int i=1; i <= 6; i++) outerelems(i) = new IVector();

	int dx=divisions_xyz(1);
	int dy=divisions_xyz(2);
	int dz=divisions_xyz(3);
	 
 	for (int iz=1; iz <= dz; iz++)
	{
	  for (int iy=1; iy <= dy; iy++)
		{
			for (int ix=1; ix <= dx; ix++)
			{
				int elnr = 5*((iz-1)*dy*dx + (iy-1)*dx + (ix-1)*1) + 1; // (block) element number of tet with node#1 (== xmin & ymin )
				if( (ix+iy+iz)%2 == 1)
				{
					if(ix ==  1) { outerelems(1)->Add(elnr); outerelems(1)->Add(elnr+3); } //1+4
					if(ix == dx) { outerelems(2)->Add(elnr+1); outerelems(2)->Add(elnr+2); } //2+3
					if(iy ==  1) { outerelems(3)->Add(elnr); outerelems(3)->Add(elnr+2); } //1+3
					if(iy == dy) { outerelems(4)->Add(elnr+1); outerelems(4)->Add(elnr+3); } //2+4
					if(iz ==  1) { outerelems(5)->Add(elnr); outerelems(5)->Add(elnr+1); } //1+2
					if(iz == dz) { outerelems(6)->Add(elnr+2); outerelems(6)->Add(elnr+3); } //3+4
				}
				else 
				{
					if(ix == dx) { outerelems(2)->Add(elnr); outerelems(2)->Add(elnr+3); } //1+4
					if(ix ==  1) { outerelems(1)->Add(elnr+1); outerelems(1)->Add(elnr+2); } //2+3
					if(iy == dy) { outerelems(4)->Add(elnr); outerelems(4)->Add(elnr+2); } //1+3
					if(iy ==  1) { outerelems(3)->Add(elnr+1); outerelems(3)->Add(elnr+3); } //2+4
					if(iz == dz) { outerelems(6)->Add(elnr); outerelems(6)->Add(elnr+1); } //1+2
					if(iz ==  1) { outerelems(5)->Add(elnr+2); outerelems(5)->Add(elnr+3); } //3+4
				}
			}
		}
	}

	IVector allouter(0);
	for(int i=1; i <= 6; i++) allouter = allouter.Union(*outerelems(i));
	return allouter.Length();
} 
