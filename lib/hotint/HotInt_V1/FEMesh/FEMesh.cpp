//#**************************************************************
//#
//# filename:             FEMesh.cpp
//#
//# author:               Gerstmayr Johannes
//#                       Dorninger Alexander
//#
//# generated:						April 2010
//# description:          Finite Element mesh
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
 

#include "mbs_interface.h"
//#include "windows.h" //for shell execute
#include <direct.h>  //for getcwd
#include  <io.h>     //file operation _access
#include  <stdio.h>
#include  <stdlib.h>
#include "stepsettings.h"
#include "element.h"
#include "material.h"
#include "FEMesh.h"
#include "myfile.h"
#include "node.h"
#include "femathhelperfunctions.h"


// here we have to explicitly include element headers
#include "FE3DHexTet.h"
#include "rigid3d.h"
#include "FiniteElement3DFFRF.h"
#include "rigid2d.h"
#include "FiniteElement2D.h"
#include "FE2DTriQuad.h"
#include "referenceframe2d.h"
#include "plate2d.h"
#include "kinematicpairs.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ (abstract base) class FEMesh:  derive FEMesh2d, FEMesh3d
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ access and manipulation of arrays
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mesh data - nodes

// creates node if not already present 
int FEMesh::AddNodeCheck(Vector3D pos, int domain, double tol)
{
	CurSel().SetType(TSetNodes);
	GetNodesInBox(CurSel(),Box3D(pos,tol));
	
	for(int i=1; i<= CurSel().NNodes(); i++)
	{
		if (domain == GetNodeDomain(CurSel(i)))
			return CurSel(i);
	}
	int n = AddPoint(pos,domain);
	nodestree.AddItem(Box3D(GetPoint3D(n),1E-10),n);
	return n;
}

// creates node if not already present 
int FEMesh::AddNodeCheck(Vector2D pos, int domain, double tol)
{
	return AddNodeCheck(Vector3D(pos.X(),pos.Y(),0.0), domain, tol);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// materials

void FEMesh::GetMaterial(int matnum, double& rho, double& Em, double& nu) const // returns parameters of a specified material
{
	rho = material_data(matnum)->Density();
	Em  = material_data(matnum)->YoungsModulus();
	nu  = material_data(matnum)->PoissonRatio();
}

void FEMesh::GetMaterial(int matnum, double& rho, double& Em, double& nu, Vector3D& col) const // returns parameters of a specified material
{
	rho = material_data(matnum)->Density();
	Em  = material_data(matnum)->YoungsModulus();
	nu  = material_data(matnum)->PoissonRatio();
	col = material_data(matnum)->GetMaterialColor();
}

int FEMesh::AddMaterial(double rho, double Em, double nu, Vector3D color) 
{
	Material* pmat = new Material(mbs);

	pmat->Density() = rho;
	pmat->YoungsModulus() = Em;
	pmat->PoissonRatio() = nu;
	pmat->SetMaterialColor(color);
	
	material_data.Add(pmat);
	return material_data.Length();
}

int FEMesh::AddMaterial(Material& mat) 
{
	Material* pmat = mat.GetCopy();
	material_data.Add(pmat);
	return material_data.Length();
}

void FEMesh::ReplaceMaterial(int i, Material& mat) 
{
	Material* old = material_data(i);
	material_data(i) = mat.GetCopy();
	delete old;
}

void FEMesh::ReplaceMaterial(int i, Material* p_mat)
{
	Material* old = material_data(i);
	material_data(i) = p_mat->GetCopy();
	delete old;
}

// add material to the material list(s), if it does not exist:
int FEMesh::AddElasticMaterial(double rho, double Em, double nu)
	{
		for (int i=1; i <= NMaterials(); i++)
		{
			if (rho == material_data(i)->Density() && Em == material_data(i)->YoungsModulus() && nu == material_data(i)->PoissonRatio())
			{
				return i;
			}
		}
		return AddMaterial(rho,Em,nu);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ initial conditions

//add initial rotational velocity field (w)x(r'-r) to all nodes
void FEMesh::AddInitialRotationalVelocity(const Vector3D& omega, const Vector3D& point_at_rot_axis)
{
	if (NN() != init_vel.Length()) InitializeInitVelocities();

	Vector3D w = omega;
	Vector3D x0 = point_at_rot_axis;
	Vector3D v,r;

	for (int i=1; i<= NN(); i++)
	{
		r = GetPoint3D(i);
		v = GetInitialNodalVelocity3D(i) + w.Cross(r - x0);
		SetInitialNodalVelocity(i, v);
	}
}

//add initial translational velocity 'vel' to all nodes
void FEMesh::AddInitialTranslationalVelocity(const Vector3D& vel)
{
	if (NN() != init_vel.Length()) InitializeInitVelocities();
	
	for (int i=1; i<= NN(); i++)
	{
		SetInitialNodalVelocity(i, GetInitialNodalVelocity3D(i) + vel);
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ quick search and mapping
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// computes box that contains all nodes
void FEMesh::ComputeNodeBox(IVector& subset)
{
	nodebox.Clear();
	if (subset.Length() == 0)	
	{ 
		subset.SetLen(NPoints()); 
		for (int i=1; i<=NPoints(); i++) subset(i)=i; 
	}
	for (int i=1; i<=subset.Length(); i++) nodebox.Add(GetPoint3D(subset(i)));
}

// initializes searchtree (sets bins, adds nodes) with current nodes in the mesh
void FEMesh::InitializeNodeSearchTree()
{
// prepare searchtree
	double tol = 1E-10;

	this->ComputeNodeBox();
	Box3D streebox(nodebox);
	streebox.InflateFactor(1.0+tol*100.0);
	//streebox.Increase(tol*100.0);

	int3 divs = GetSearchTreeDivisions(streebox, NN());
	nodestree.ResetSearchTree(divs(1), divs(2), divs(3), streebox);

	//nodestree.ResetSearchTree(bins_x, bins_y, bins_z,streebox);
	for(int i=1; i<=NN(); i++) 
		nodestree.AddItem(Box3D(GetPoint3D(i),tol),i);
}

// resizes the Searchtree (adds new region) - used in GenerateBlock
void FEMesh::ResizeNodeSearchTree(Box3D addbox, int addnodes)
{
	ComputeNodeBox();
	Box3D box(GetNodeBox()); // box containing all nodes already in mesh
	
	if(! (box.IsIn(addbox.PMax()) && box.IsIn(addbox.PMin())) ) // !(2x point is in box) --> Addbox is in old box, no update needed
	{
		box.Add(addbox); // resize
		int3 divs = GetSearchTreeDivisions(box, NN()+addnodes);
		nodestree.ResetSearchTree(divs(1), divs(2), divs(3), box);
		for(int i=1; i<=NN(); i++) 
			nodestree.AddItem(Box3D(GetPoint3D(i),1e-8),i);
	}
}

// compute divisions for the Searchtree, cells should be roughly cubes
int3 FEMesh::GetSearchTreeDivisions(Box3D& box, int nnodes)
{
	int nodespercell = 10; // hardcoded value... is there an optimum cellsize?
	double celllen;
	int bins_x ,bins_y, bins_z;

	if (box.SizeZ() == 0.0) // 2D case: z=0
	{
		double boxarea = box.SizeX()*box.SizeY();
		double cellarea = boxarea / NN() * nodespercell;
		celllen = pow(cellarea,1./2.);
		bins_x = Maximum(1,(int)( box.SizeX()/celllen )+1);
		bins_y = Maximum(1,(int)( box.SizeY()/celllen )+1);
		bins_z = 1;
		box.Increase(0.0,0.0,0.5); // set sizeZ to 1 // needed ?
	}
	else // 3D case
	{
		double boxvol = box.SizeX()*box.SizeY()*box.SizeZ();
		double cellvol = boxvol / NN() * nodespercell;
		celllen = pow(cellvol,1./3.);

		bins_x = Maximum(1,(int)( box.SizeX()/celllen )+1);
		bins_y = Maximum(1,(int)( box.SizeY()/celllen )+1);
		bins_z = Maximum(1,(int)( box.SizeZ()/celllen )+1);
	}
	return int3(bins_x, bins_y, bins_z);
}

// generates a list that holds a list of elements the node is part for each node
void FEMesh::ComputeNodesToElements()
{
	int i;
	if (nodes_to_elements.Length() != NN())
	{
		ReleaseArray_TemplatePtr(nodes_to_elements);
		nodes_to_elements.SetLen(NN());

		for (i=1; i <= NN(); i++)
		{
			nodes_to_elements(i) = new IVector(10);
		}
	}
	for (i=1; i <= NN(); i++)
	{
		nodes_to_elements(i)->SetLen(0);
	}
	for (i=1; i <= elements.Length(); i++)
	{
		for (int j=1; j <= elements(i)->NNodes(); j++)
		{
			int nn = elements(i)->GetNode(j);
			if (nodes_to_elements(nn)->Last() != i) 
				nodes_to_elements(nn)->Add(i); //element could have double nodes
		}
	}
}

// access node_to_element_list, automatic recomppute if number of entries differs from number of nodes  
TArray<TArray<int>*>& FEMesh::NodesToElementList(int flag_autorecompute)
{
 	if(nodes_to_elements.Length() != NN ())
	{
		ComputeNodesToElements();
	}
	return nodes_to_elements;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ indirect - statistics - these quantities are computed
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// checks all used bodynumbers - updates mesh_register
int FEMesh::CountBodies(IVector& subset)
{
	if(subset.Length() == 0)
	{
		subset = NaturalNumbers(NElements());
	}

	mesh_register.Bodies().Flush();
	for(int i=1; i<=subset.Length(); i++)
	{
		mesh_register.AddIfNotExists( GetElement(subset(i)).Domain(), TSetBodies );
	}

	return mesh_register.NBodies();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ import functions - read from external files
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// general import functions

// loads Nodes data from various sources
int FEMesh::LoadNodes(const mystr & filename, const mystr & format)
{
	// *** retvalues: -2: header not found | -1: #expected!=#read | n: nodes read ***
	CMFile file_nodes(filename,TFMread);
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_warn) << "No Filename for Nodesfile specified !\n";
		return -2;
	}
	if (!file_nodes.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename.c_str();
		mbs->UO().InstantMessageText("Could not open Nodes File");
		return -2;
	}

// Data: nodenumber | x-coord | y-coord | z-coord  (Mesh3D) or nodenumber | x-coord | y-coord  (Mesh2D)
	Vector3D npos;
	TArray<int> nodenumber;
	TArray<Vector3D> nodeposition;

	if (format==mystr("ANSYS") || format==mystr("HOTINT"))
	{
// format ANSYS
// Header:  [NODES]
// Header2: total number of nodes (int)  
// Data: nodenumber | x-coord | y-coord | z-coord   or nudenumber | x-coord | y-coord
	  mystr blocktag("[NODES]");
	  mystr buffer("");
		bool foundblock = false;      
		int totalnumberofnodes = 0;

		mbs->UO(UO_LVL_ext) << mystr("Looking for Header [NODES] in File: ") + filename + mystr(" (from ANSYS)\n");
// find and read Blockheader
		while( !file_nodes.EndOfFile() )
		{
			file_nodes.RWSoleStrBrackets(buffer);
			if (file_nodes.EndOfFile()) { mbs->UO(UO_LVL_err) << mystr("No header found in ") + filename + mystr("\n"); return -2; }
			if (blocktag.Compare(buffer)) foundblock = true;
			if (!foundblock) continue;
// read 2nd line of Blockheader
			totalnumberofnodes = (int) file_nodes.GetRWDouble();
			break;
		}
		mbs->UO(UO_LVL_ext) << mystr(totalnumberofnodes) + mystr(" nodes expected. Reading ...");
// prepare intermediate arrays
//		nodenumber.SetLen(totalnumberofnodes); nodenumber.SetAll(0);
//		nodeposition.SetLen(totalnumberofnodes); nodeposition.SetAll(Vector3D());

		while (!file_nodes.EndOfFile())
		{
			file_nodes.RWSoleStrBrackets(buffer); // read first entry of new line
			if (file_nodes.EndOfFile()) break;	  // end of file -> stop reading nodes
			if ( buffer[0] == '[') break;         // different tag -> stop reading nodes
			nodenumber.Add((int)buffer.MakeDouble());   

			file_nodes.RWDouble( npos.X());
			file_nodes.RWDouble( npos.Y());
			if (GetMeshDim() == 3) file_nodes.RWDouble( npos.Z());

			nodeposition.Add(npos);
		}
		mbs->UO(UO_LVL_ext) << mystr(nodenumber.Length()) + mystr(" nodes read.\n");
	} // end format ANSYS
	//else if (format==mystr("NETGEN"))
	//{
	//	// format NETGEN
	//	mystr blocktag("points");
	//	mystr buffer("");
	//	mystr mystr_file("");
	//	int pos=0;

	//	bool foundblock = false;      
	//	int totalnumberofnodes = 0;

	//	mbs->UO(UO_LVL_ext) << mystr("reading file: ") + filename + mystr(" (from NETGEN)\n");
	//	filename.RWF(mystr_file);
	//	mbs->UO(UO_LVL_ext) << mystr("Looking for Header 'points' in buffer. \n");

	//	// find and read Blockheader
	//	while( mystr_file.GetUntil(pos,'\n',buffer,1))
	//	{
	//		if(pos==-1) { mbs->UO(UO_LVL_err) << mystr("No header found in ") + filename + mystr("\n"); return -2; }
	//		if(buffer.Compare(blocktag)) foundblock = true;
	//		if (!foundblock) continue;
	//		// read 2nd line of Blockheader
	//		totalnumberofnodes = (int) file_nodes.GetRWDouble();
	//		break;
	//	}
	//	mbs->UO(UO_LVL_ext) << mystr(totalnumberofnodes) + mystr(" nodes expected. Reading ...");
	//	// prepare intermediate arrays
	//	//		nodenumber.SetLen(totalnumberofnodes); nodenumber.SetAll(0);
	//	//		nodeposition.SetLen(totalnumberofnodes); nodeposition.SetAll(Vector3D());

	//	//while (!file_nodes.EndOfFile())
	//	//{
	//	//	file_nodes.RWSoleStrBrackets(buffer); // read first entry of new line
	//	//	if (file_nodes.EndOfFile()) break;	  // end of file -> stop reading nodes
	//	//	if ( buffer[0] == '[') break;         // different tag -> stop reading nodes
	//	//	nodenumber.Add((int)buffer.MakeDouble());   

	//	//	file_nodes.RWDouble( npos.X());
	//	//	file_nodes.RWDouble( npos.Y());
	//	//	if (GetMeshDim() == 3) file_nodes.RWDouble( npos.Z());

	//	//	nodeposition.Add(npos);
	//	//}
	//	//mbs->UO(UO_LVL_ext) << mystr(nodenumber.Length()) + mystr(" nodes read.\n");		


	//}// end format NETGEN
	else 
	{
		mbs->UO(UO_LVL_err) << mystr("No valid format specified for LoadNodes");
		return -2;
	}

// check highest nodenumber 
	int highestnodenumber = 0;
//	for (int i=1; nodenumber.Length(); i++)
	//	if(nodenumber(i) > highestnodenumber) highestnodenumber = nodenumber(i);
	
// copy to Mesh
	for (int i=1; i<=nodenumber.Length(); i++)
	{
		AddPoint(nodeposition(i)); // AddPoint checks meshdimension, default domain 1 !
	}
	ComputeNodeBox();
	InitializeNodeSearchTree();
	if (NPoints() == nodenumber.Length()) return nodenumber.Length();
	else return -1;
}

// loads Elements data from various sources, writes them into chosen Array (which must be defined in FEMesh)
int FEMesh::LoadElements(const mystr& filename, const mystr& format, const mystr& writeto)
{
// *** retvalues: -2: header not found | -1: #expected!=#read | n: elements read
	CMFile file_elements(filename,TFMread);
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_warn) << "No Filename for Elementsfile specified !\n";
		return -2;
	}
	if (!file_elements.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename.c_str();
		mbs->UO().InstantMessageText("Could not open Elements File");
		return -2;
	}

// Data: elementnumber | e.dimension | e.order | nodes used :<=N | node1 | ... | nodeN | e.material
	TArray<int> elementnumber;
	TArray<int> elementdimension;
	TArray<int> elementorder;
	TArray<int> nodesused;
	TMatrix<int> nodelist;
	TArray<int> elementmaterial;

	if(format==mystr("ANSYS") || format==mystr("HOTINT"))
	{
// format ANSYS
// Header:  [ELEMENTS]  or  [FACES]
// Header2: total number of elements (int) | columns used for nodes (int) := N
// Data: elementnumber | e.dimension | e.order | nodes used :<=N | node1 | ... | nodeN | e.material

		mystr blocktag;
		if(writeto==mystr("elements")) blocktag = mystr("[ELEMENTS]");
		if(writeto==mystr("faces")) blocktag = mystr("[FACES]");

		mystr buffer("");
		bool foundblock = false;      
		int totalnumberofelements = 0;
		int columnsfornodes = 0;
		double dummy = 0.0;

		mbs->UO(UO_LVL_ext) << mystr("Looking for Header ") + blocktag + mystr(" in File: ") + filename + mystr(" (from ANSYS)\n");
// find and read Blockheader
		while (!file_elements.EndOfFile())
		{
			file_elements.RWSoleStrBrackets(buffer);
			if (file_elements.EndOfFile()) { mbs->UO(UO_LVL_err) << mystr("no header found in ") + filename.c_str() + mystr("\n"); return -2; }
			if (blocktag.Compare(buffer)) foundblock = true;
			if (!foundblock) continue;
	// read 2nd line of Blockheader
			totalnumberofelements = (int) file_elements.GetRWDouble();
			columnsfornodes = (int) file_elements.GetRWDouble();
			break;
		}
		mbs->UO(UO_LVL_ext) << mystr(totalnumberofelements) + mystr(" elements expected. Reading ...");

		nodelist.SetDim(totalnumberofelements,columnsfornodes);

		while(!file_elements.EndOfFile())
		{
			file_elements.RWSoleStrBrackets(buffer); // read first entry of new line
			if (file_elements.EndOfFile()) break;	   // end of file -> stop reading elements
			if (buffer[0] == '[') break;             // different tag -> stop reading elements
			int i = (int)buffer.MakeDouble(); 
			elementnumber.Add(i);
		
			elementdimension.Add((int) file_elements.GetRWDouble());
			elementorder.Add((int) file_elements.GetRWDouble());
			nodesused.Add((int) file_elements.GetRWDouble());
			for (int j=1; j<=columnsfornodes; j++)
			{
				nodelist(i,j) = (int)file_elements.GetRWDouble();
			}
			elementmaterial.Add( (int)file_elements.GetRWDouble());		
		}
		mbs->UO(UO_LVL_ext) << elementnumber.Length() <<" elements read.\n";
	}// end format ANSYS

// copy to Elements Array
	if(writeto==mystr("elements")) 
	{
		TArray<int> points;
		for (int i=1; i<=elementnumber.Length(); i++)
		{
			points.SetLen(0); // get nodes for the element in TArray
			for (int j=1; j<=nodelist.NCols(i); j++)
			{
				points.Add(nodelist(i,j));
			}
			if(format == mystr("ANSYS"))
			{
				MapNodesAnsys2Hotint(elementdimension(i),elementorder(i),nodesused(i),points); // Ansys node order is different than Hotint node order 
			}
			// 3D elements:
			// Hex quad and its degenerate elements ( save indices of degenerate elements in "current selection"
			if ( ((elementdimension(i)==3) && (elementorder(i)==2) && (nodesused(i)==20)) ||  //Hex quad
				((elementdimension(i)==3) && (elementorder(i)==2) && (nodesused(i)==15)) ||  //Prism quad
				((elementdimension(i)==3) && (elementorder(i)==2) && (nodesused(i)==13)) )	  //Pyram quad
			{
				Elements().Add(MakeHexquad(points,elementmaterial(i)));
				if (nodesused(i) == 15) 
				{
					mbs->UO(UO_LVL_dbg1 /*UO_LVL_warn*/) << mystr("Warning: Element #") + mystr(i) + mystr(" is Prism\n");
					CurSel().Add(i,TSetElements);
				}
				if (nodesused(i) == 13) 
				{
					mbs->UO(UO_LVL_dbg1 /*UO_LVL_warn*/) << mystr("Warning: Element #") + mystr(i) + mystr(" is Pyramid\n");
					CurSel().Add(i,TSetElements);	
				}
			}
			// Tet quad
			else if ((elementdimension(i)==3) && (elementorder(i)==2) && (nodesused(i)==10)) //Tet quad
			{
				Elements().Add(MakeTetquad(points,elementmaterial(i)));
			}
			// Hex and its degenerate elements
			else if ( ((elementdimension(i)==3) && (elementorder(i)==1) && (nodesused(i)==8)) ||  //Hex
				((elementdimension(i)==3) && (elementorder(i)==1) && (nodesused(i)==6)) )   //Prism
			{
				Elements().Add(MakeHex(points,elementmaterial(i)));
				if (nodesused(i) == 6) 
				{
					mbs->UO(UO_LVL_dbg1 /*UO_LVL_warn*/) << "Warning: Element #" << i << " is Prism\n";
					CurSel().Add(i,TSetElements);
				}
			}
			// 2D elements:
			// Quad quad
			else if ( ((elementdimension(i)==2) && (elementorder(i)==2) && (nodesused(i)==8)) )   //Quad quad
			{
			// !!! additional node in the middle of the element is required !!! FEMesh has only 9 node element !!!
				int center_nodenr = AddCenterNode(points);
				Elements().Add(MakeQuadquad(points,elementmaterial(i)));
			}
			// Trig quad
			else if ( ((elementdimension(i)==2) && (elementorder(i)==2) && (nodesused(i)==6)) )   //Trig quad
			{
				Elements().Add(MakeTrigquad(points,elementmaterial(i)));
			}
			// Quad
			else if ( ((elementdimension(i)==2) && (elementorder(i)==1) && (nodesused(i)==4)) )   //Quad
			{
				Elements().Add(MakeQuad(points,elementmaterial(i)));
			}
			// Trig
			else if ( ((elementdimension(i)==2) && (elementorder(i)==1) && (nodesused(i)==3)) )   //Trig
			{
				Elements().Add(MakeTrig(points,elementmaterial(i)));
			}
			// 1D elements:
			// LineQuad
			else if ( ((elementdimension(i)==1) && (elementorder(i)==2) && (nodesused(i)==3)) )		//Line quad
			{
				Elements().Add(MakeLinequad(points,elementmaterial(i)));
			}
			// Line
			else if (((elementdimension(i)==1) && (elementorder(i)==1) && (nodesused(i)==2)) )	  //Line
			{
				Elements().Add(MakeLine(points,elementmaterial(i)));
			}
			else
				mbs->UO(UO_LVL_err) << "Element #" << i << " has no implemented/known type\n";
			// add here and in MapNodesAnsys2Hotint(..)
		}
		if (Elements().Length() == elementnumber.Length())
			return elementnumber.Length();
		else 
			return -1;
	}

// copy to FEFaces Array
	if(writeto==mystr("faces"))
	{
		TArray<int> points;
		for (int i=1; i<=elementnumber.Length(); i++)
		{
			points.SetLen(0); // get nodes for the face 
			for (int j=1; j<=nodelist.NCols(i); j++)
			{
				points.Add(nodelist(i,j));
			}
			if(format == mystr("ANSYS"))
				MapNodesAnsys2Hotint(elementdimension(i),elementorder(i),nodesused(i),points); // Ansys node order is different than Hotint node order 
			// 2D elements:
			// Quad quad
			if ( ((elementdimension(i)==2) && (elementorder(i)==2) && (nodesused(i)==8)) )   //Quad quad
			{
				Faces().Add(FEMesh_Face(int4(points(1),points(2),points(3),points(4)), 0, 0));
			}
			// Trig quad
			else if ( ((elementdimension(i)==2) && (elementorder(i)==2) && (nodesused(i)==6)) )   //Trig quad
			{
				Faces().Add(FEMesh_Face(int4(points(1),points(2),points(3),0), 0, 0));
			}
			// Quad
			else if ( ((elementdimension(i)==2) && (elementorder(i)==1) && (nodesused(i)==4)) )   //Quad
			{
				Faces().Add(FEMesh_Face(int4(points(1),points(2),points(3),points(4)), 0, 0));
			}
			// Trig
			else if ( ((elementdimension(i)==2) && (elementorder(i)==1) && (nodesused(i)==3)) )   //Trig
			{
				Faces().Add(FEMesh_Face(int4(points(1),points(2),points(3),0), 0, 0));
			}
			else
				// 3D elements:
				// 1D elements:
				// other elements:
				mbs->UO(UO_LVL_err) << "Element #" << i << " has no implemented/known type\n";
		}
		if (Faces().Length() == elementnumber.Length()) 
			return elementnumber.Length();
		else 
			return -1;
	}
	return -1;
}

// loads Materials data from File - replaces any previous entries in materials_* arrays
int FEMesh::LoadMaterials(const mystr & filename)
{
	// *** retvalues: -2: header not found | -1: #expected!=#read | n: materials read
	// Header:  [MATERIALS]
	// Header2: total number of materials (int) 
	// Data: materialnumber | density | youngs_modulus | poisson_ratio
	CMFile file_materials(filename,TFMread);
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_err) << "No Filename for Materials specified !\n";
		return -2;
	}
	if (!file_materials.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename.c_str();
		mbs->UO().InstantMessageText("Could not open Materials File");
		return -2;
	}
	mystr blocktag("[MATERIALS]");
	mystr buffer("");
	bool foundblock = false;      
	int totalnumberofmaterials = 0;

	mbs->UO(UO_LVL_all) << "Looking for Header [MATERIALS] in File: " << filename.c_str() << " (from ANSYS)\n";
// find and read Blockheader
	while (!file_materials.EndOfFile())
	{
		file_materials.RWSoleStrBrackets(buffer);
		if (file_materials.EndOfFile()) return -2;		
		if (blocktag.Compare(buffer)) foundblock = true;
		if (!foundblock) continue;
// read 2nd line of Blockheader
		totalnumberofmaterials = (int)file_materials.GetRWDouble();
		break;
	}
	if(GetMBS()->UO().GetGlobalMessageLevel()==UO_LVL_dbg1)	//$ DR 2011-09-15
	{
		mbs->UO(UO_LVL_dbg1) << totalnumberofmaterials <<" materials expected. Reading ...";
	}
	int numberofcolumns=0;

	TArray<int> materialnumber; 

	mystr tag_density("[density]");
	mystr tag_youngs("[youngs_modulus]");
	mystr tag_poisson("[poisson_ratio]");
	mystr tag_colorRGB("[drawing_color]");
	mystr tag_name("[name]");

	// Set all materials to defaults
	material_data.Flush();


	//material_data.SetLen(totalnumberofmaterials);  
	for(int i = 1; i <= totalnumberofmaterials; i++)
	{
		material_data.Add(new Material(this->mbs));
	}

// Data: MaterialNumber (int)  |  Number of Parameters (int) |   {n * ParameterName (str)}  |  {n* ParameterValue}
// Number and Order of parameters can change for each material ! 
	while (!file_materials.EndOfFile())
	{
		file_materials.RWSoleStrBrackets(buffer);
		if (file_materials.EndOfFile()) break;	
		if (buffer[0] == '[') break;

		materialnumber.Add((int)buffer.MakeDouble());
		if (materialnumber.Length() > totalnumberofmaterials)
		{
			mbs->UO(UO_LVL_err) << mystr("WARNING from FEMesh::LoadMaterials : Material file contains more Materials then stated in Header ( ")
													 + mystr(totalnumberofmaterials) +mystr(" ), no further Materials are read! \n");
			break;
		}
		numberofcolumns = (int)file_materials.GetRWDouble();
// read all ParameterNames
		TArray<mystr*> colnames; 
		colnames.SetLen(numberofcolumns); 
		for (int i = 1; i <= numberofcolumns; i++)
		{
			file_materials.RWSoleStrBrackets(buffer);
			colnames(i) = new mystr(buffer);
		}
// read all ParameterValues to material_* arrays
		for (int i = 1; i <= numberofcolumns; i++)
		{
			if (*colnames(i) == tag_density) 
			{
				material_data(materialnumber.Last())->Density() = file_materials.GetRWDouble();
			}
			if (*colnames(i) == tag_youngs)
			{
				material_data(materialnumber.Last())->YoungsModulus() = file_materials.GetRWDouble();
			}
			if (*colnames(i) == tag_poisson)
			{
				material_data(materialnumber.Last())->PoissonRatio() = file_materials.GetRWDouble();
			}
			if (*colnames(i) == tag_colorRGB)
			{
					mystr str_rgb;
					Vector3D mat_color;

					file_materials.RWSoleStr(str_rgb);
					mat_color.X() = str_rgb.MakeDouble();
					file_materials.RWSoleStr(str_rgb);
					mat_color.Y() = str_rgb.MakeDouble();
					file_materials.RWSoleStr(str_rgb);
					mat_color.Z() = str_rgb.MakeDouble();
					material_data(materialnumber.Last())->SetMaterialColor(mat_color);
			}
			if (*colnames(i) == tag_name)
			{
				mystr mat_name;
				file_materials.RWSoleStr(mat_name);
				material_data(materialnumber.Last())->GetMaterialName() = mat_name;
			}
		}
		if(GetMBS()->UO().GetGlobalMessageLevel()==UO_LVL_dbg1)	//$ DR 2011-09-15
		{
		mbs->UO(UO_LVL_dbg1) << mystr("material ") + mystr(materialnumber.Last()) 
			+ mystr(": density= ") + mystr(material_data(materialnumber.Last())->Density()) 
			+ mystr(", youngs= ") + mystr(material_data(materialnumber.Last())->YoungsModulus())
			+ mystr(", poisson= ") + mystr(material_data(materialnumber.Last())->PoissonRatio()) 
			+ mystr(", color= ") + (material_data(materialnumber.Last())->GetMaterialColor()).MakeString() + mystr("\n");
		}

	}
	if (totalnumberofmaterials == materialnumber.Length()) return materialnumber.Length();
	else return -1;
}

int FEMesh::LoadMaterialsEDC(const mystr & filename)		// loads Materials data from File, assuming EDC structure there - replaces any previous entries in materials_* arrays //$ DR 2012-10
{
	GetMBS()->UO(UO_LVL_0) << "read " << filename << "\n";

	ElementDataContainer edc_file;
	ElementDataContainer* edcP_Mat;
	ElementData* edP;
	int index;

	// read parameters from file and store it in ElementDataContainer
	int rv = GetMBS()->File2EDC(filename, &edc_file);
	if(!rv)return rv; 

	int totalnumberofmaterials =edc_file.Length();
	mbs->UO(UO_LVL_dbg1) << totalnumberofmaterials <<" materials expected. Reading ...";

	// Set all materials to defaults
	material_data.Flush();


	//material_data.SetLen(totalnumberofmaterials);  
	for(int i = 1; i <= totalnumberofmaterials; i++)
	{
		material_data.Add(new Material(this->mbs));
	}

	for(int i=1; i<=totalnumberofmaterials; i++)
	{
		index = 0;
		index  = edc_file.Find(mystr("Material")+mystr(i));
		if(index)
		{
			edcP_Mat = edc_file.GetPtr(index)->GetEDC();
			edP = edcP_Mat->TreeFind("Mechanics.density");
			if(edP){ material_data(i)->Density() = edcP_Mat->TreeGetDouble("Mechanics.density");}
			edP = edcP_Mat->TreeFind("Mechanics.youngs_modulus");
			if(edP){ material_data(i)->YoungsModulus() = edcP_Mat->TreeGetDouble("Mechanics.youngs_modulus");}
			edP = edcP_Mat->TreeFind("Mechanics.poisson_ratio");
			if(edP){ material_data(i)->PoissonRatio() = edcP_Mat->TreeGetDouble("Mechanics.poisson_ratio");}
			edP = edcP_Mat->TreeFind("Graphics.color");
			if(edP)
			{ 
				Vector3D vtemp;
				edcP_Mat->TreeGetVector3D("Graphics.color",vtemp.X(),vtemp.Y(),vtemp.Z());
				material_data(i)->SetMaterialColor(vtemp);
			}
			edP = edcP_Mat->TreeFind("name");
			if(edP){ material_data(i)->GetMaterialName() = edcP_Mat->TreeGetString("name");}
		}
		else
		{
			return -1;
		}
	}

	return totalnumberofmaterials;
}

// loads Loads data - Not Implemented yet 
int FEMesh::LoadLoads(const mystr & filename, const mystr & format)
{
	mbs->UO(UO_LVL_warn) <<" ********** \n WARNING: function LoadLoads renamed to LoadLoads_FORMAT1 on 11-04-11 \n please change in your model file a.s.a.p. \n ********** \n";
	mbs->UO().InstantMessageText("this function is not implemented yet"); 
	return 0;
}

// loads Constraints data - Not Implemented yet 
int FEMesh::LoadConstraints(const mystr & filename, const mystr & format)
{
	mbs->UO(UO_LVL_warn) <<" ********** \n WARNING: function LoadLoads renamed to LoadLoads_FORMAT1 on 11-04-11 \n please change in your model file a.s.a.p. \n ********** \n";
	mbs->UO().InstantMessageText("this function is not implemented yet"); 
	return 0;
}

// loads a mapping that combines several faces to areas 
int FEMesh::LoadFaces2Areas(const mystr & filename)
{
	// *** retvalues: -2: header not found | -1: #expected!=#read | n: areas read
	// Header:  [AELEM]
	// Header2: total number of areas(int) 
	// Data: AreaNumber(int) | NumberOfFaces(int) | Face#1 ... Face#n
	CMFile file_areas(filename,TFMread);
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_err) << "No Filename for Faces2Areas specified !\n";
		return -2;
	}
	if (!file_areas.IsGood())
	{
		mbs->UO(UO_LVL_err) << mystr("Could not open File: ") + filename + mystr("\n");
		mbs->UO().InstantMessageText("Could not open Faces2Areas File");
		return -2;
	}

	mystr blocktag("[AELEM]");
	mystr buffer("");
	bool foundblock = false;      
	int totalnumberofareas = 0;

	mbs->UO(UO_LVL_all) << mystr("Looking for Header [AELEM] in File: ") + filename + mystr(" (from ANSYS)\n");
// find and read Blockheader
	while (!file_areas.EndOfFile())
	{
		file_areas.RWSoleStrBrackets(buffer);
		if (file_areas.EndOfFile()) { mbs->UO(UO_LVL_err) << mystr("no header found in ") + filename + mystr("\n"); return -2; }	
		if (blocktag.Compare(buffer)) foundblock = true;
		if (!foundblock) continue;
// read 2nd line of Blockheader
		totalnumberofareas = (int)file_areas.GetRWDouble();
		break;
	}
	if(GetMBS()->UO().GetGlobalMessageLevel()==UO_LVL_dbg1)	//$ DR 2011-09-15
	{
		mbs->UO(UO_LVL_dbg1) << mystr(totalnumberofareas) + mystr(" aresa expected. Reading ...");
	}
	TArray<int> areanumber;
	Areas().Flush();
	Areas().SetLen(totalnumberofareas);
	Areas().SetAll(FEMesh_Set(TSetFaces));

	while(!file_areas.EndOfFile())
	{
		file_areas.RWSoleStrBrackets(buffer);		// read first entry of new line
		if (file_areas.EndOfFile()) break;			// end of file -> stop reading elements
		if (buffer[0] == '[') break;            // different tag -> stop reading elements
		int i = (int)buffer.MakeDouble(); 
		areanumber.Add(i);
		
		int facesinarea = (int)file_areas.GetRWDouble();
		for (int j=1; j<=facesinarea; j++)
		{
			GetArea(areanumber.Last()).Faces().Add(	(int)file_areas.GetRWDouble() );	
		}
	}
	mbs->UO(UO_LVL_ext) << areanumber.Length() <<" areas read.\n";

	if (totalnumberofareas == areanumber.Length()) return areanumber.Length();
	else return -1;
}

// load data generated with Ansys2Hotint.mac - BlockTag: [NODES]
int FEMesh::LoadAnsys2HotintNODES(const mystr & filename)
{
	return LoadNodes(filename,"ANSYS");
}

// load data generated with Ansys2Hotint.mac - BlockTag: [ELEMENTS]
int FEMesh::LoadAnsys2HotintELEMENTS(const mystr & filename)
{
	return LoadElements(filename,mystr("ANSYS"),mystr("elements"));
}

// load data generated with Ansys2Hotint.mac - BlockTag: [FACES]
int FEMesh::LoadAnsys2HotintFACES(const mystr & filename)
{
	return LoadElements(filename,mystr("ANSYS"),mystr("faces"));
}

// loads Loads data - Text file is in a specific format
int FEMesh::LoadLoads_FORMAT1(const mystr & filename)
{
	// *** retvalues: -2: header not found | -1: #expected!=#read | n: materials read
	// Header:  [LOADS]
	// Header2: total number of loads(int) | Units(string)
	// Data: LoadType(string) | LoadTarget#(int) | X-Component LoadVector(float) | Y-Component(float) | Z-Component(float)
	CMFile file_loads(filename,TFMread);
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_err) << "No Filename for Loadsfile specified !\n";
		return -2;
	}
	if (!file_loads.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename.c_str() << "\n";
		mbs->UO().InstantMessageText("Could not open Loads File");
		return -2;
	}
	mystr blocktag("[LOADS]");
	mystr buffer("");
	bool foundblock = false;      
	int totalnumberofloads = 0;
	mystr units("");

	mbs->UO(UO_LVL_ext) << "Looking for Header [LOADS] in File: " << filename.c_str() << " (from ANSYS)\n";
// find and read Blockheader
	while (!file_loads.EndOfFile())
	{
		file_loads.RWSoleStrBrackets(buffer);
		if (file_loads.EndOfFile()) { mbs->UO(UO_LVL_err) << mystr("no header found in ") + filename + mystr("\n"); return -2; }
		if (blocktag.Compare(buffer)) foundblock = true;
		if (!foundblock) continue;
// read 2nd line of Blockheader
		totalnumberofloads = (int)file_loads.GetRWDouble();
		file_loads.RWSoleStrBrackets(units);
		break;
	}
	mbs->UO(UO_LVL_ext) << totalnumberofloads <<" loads expected. Reading ...";

// Data: LoadType(string) | LoadTarget#(int) | X-Component LoadVector(float) | Y-Component(float) | Z-Component(float)
	int loadtarget;
	Vector3D loadvector;
	int nr_of_loadsteps=0;
	TArray<StepSettings> loadsteps;
// implemented LoadType identifier -> any changes in load type numbering: change in this routine, in AddLoadsToMBS() and very likely in mirrormesh, linear2quadratic,...
	mystr loadtype_areaload("AreaSurfaceLoad");					// LoadType TAreaLoad
	mystr loadtype_areaconstraint("AreaConstraint");    // LoadType TAreaConstraint
	mystr loadtype_nodalconstraint("NodalConstraint");  // LoadType TNodalConstraint
	mystr loadtype_bodyload("BodyLoad");                // LoadType TBodyLoad
	mystr loadtype_faceload("FaceLoad");							  // LoadType TFaceLoad
	mystr loadtype_faceconstraint("FaceConstraint");    // LoadType TFaceConstraint
	mystr loadtype_areacontact("AreaContact");					// LoadType TAreaContact

	while (!file_loads.EndOfFile())
	{
		file_loads.RWSoleStrBrackets(buffer); // first entry of new line
		if (file_loads.EndOfFile()) break;	  // end of file -> stop reading loads
		if (buffer[0] == '[') break;          // different tag -> stop reading loads
		if ((int)buffer.MakeDouble() != 0) continue;  // skip any (integer) numbers at the beginning of the line - optional number of load

//$ AD 2011-04-11 new: loads + constraints		
		if ( buffer == loadtype_areacontact )  //$ LA 2011-06-08: special case for area_contact to read 3D-stiffness
		{
			loadtarget = (int) file_loads.GetRWDouble();
			int otherarea = (int) file_loads.GetRWDouble();
			file_loads.RWDouble( loadvector.X());
			file_loads.RWDouble( loadvector.Y());
			file_loads.RWDouble( loadvector.Z());
			double damping = file_loads.GetRWDouble();
			
			nr_of_loadsteps = (int) file_loads.GetRWDouble();
			loadsteps.Flush();
			loadsteps.SetLen(nr_of_loadsteps);
			for(int i=1; i <= nr_of_loadsteps; i++)
			{
				double lf; 
				file_loads.RWDouble(lf);
				loadsteps(i).SetStepSettings(lf, TCSRLinear);
			}

			// content of LV:=  stiffnessvector
			int penaltyflag = (loadvector.X() != 0 || loadvector.Y() != 0 || loadvector.Z() != 0);
			AddConstraint(FEMesh_AreaContact(loadtarget, otherarea, penaltyflag, loadvector, Vector3D(damping)));		
		}
		else
		{
			loadtarget = (int) file_loads.GetRWDouble();
			file_loads.RWDouble( loadvector.X());
			file_loads.RWDouble( loadvector.Y());
			file_loads.RWDouble( loadvector.Z());

			nr_of_loadsteps = (int) file_loads.GetRWDouble();
			loadsteps.Flush();
			loadsteps.SetLen(nr_of_loadsteps);
			for(int i=1; i <= nr_of_loadsteps; i++)
			{
				double lf; 
				file_loads.RWDouble(lf);
				loadsteps(i).SetStepSettings(lf, TCSRLinear);
			}
		}

		if (loadtype_areaload == buffer)	
		{
// content of LV := Vector - Load Vector for a unit area
			AddLoad(FEMesh_AreaLoad(loadtarget, loadvector, loadsteps));
		}
		else if (loadtype_faceload == buffer) 
		{
// content of LV := Vector - Load Vector for a unit area
			AddLoad(FEMesh_FaceLoad(loadtarget, loadvector, loadsteps));
		}
		else if (loadtype_bodyload == buffer)  
		{
// content of LV:=  x-> strength of load, y-> axis, z->flag (0=bodyload, 1=gravityload)
			Vector3D forcevector(0.);
			double strength = loadvector.X();
			int axis = (int) loadvector.Y();
			int gravity = (int) loadvector.Z();
			if (axis == 1) forcevector.X() = strength;
			if (axis == 2) forcevector.Y() = strength;
			if (axis == 3) forcevector.Z() = strength;

			AddLoad(FEMesh_BodyLoad(loadtarget, forcevector, gravity, loadsteps));
		}
		else if (loadtype_areaconstraint == buffer)
		{
// content of LoadVector:=  x-> constrained axis , y-> springstiffness, z->damping	
			int penaltyflag = (loadvector.Y() != 0);
		
			AddConstraint(FEMesh_AreaConstraint(loadtarget, (int) loadvector.X(), penaltyflag, Vector3D(loadvector.Y()), Vector3D(loadvector.Z())));
		}
		else if (loadtype_faceconstraint == buffer)
		{
// content of LoadVector:=  x-> constrained axis , y-> springstiffness, z->damping	
			int penaltyflag = (loadvector.Y() != 0);

			AddConstraint(FEMesh_FaceConstraint(loadtarget, (int) loadvector.X(), penaltyflag, Vector3D(loadvector.Y()), Vector3D(loadvector.Z())));
		}
		else if (loadtype_nodalconstraint == buffer)
		{
// content of LoadVector:=  x-> constrained axis , y-> springstiffness, z->damping	
			int penaltyflag = (loadvector.Y() != 0);
			
			AddConstraint(FEMesh_NodeConstraint(loadtarget, (int) loadvector.X(), penaltyflag, Vector3D(loadvector.Y()), Vector3D(loadvector.Z())));
		}

		else if (loadtype_areacontact == buffer)
		{
			// already treated when read (see above)
//// content of LV:=  x-> 2nd area number , y-> springstiffness, z->damping
//			int penaltyflag = (loadvector.Y() != 0);
//			int otherarea = (int) loadvector.X();
//			AddConstraint(FEMesh_AreaContact(loadtarget, otherarea, penaltyflag, Vector3D(loadvector.Y()), Vector3D(loadvector.Z())));
		}			
		else ; // unknown LoadType

	} // end while

	mbs->UO(UO_LVL_ext) << load_data.Length() <<" Loads read.\n";
	mbs->UO(UO_LVL_ext) << constraint_data.Length() <<" Constraints read.\n";

	if (totalnumberofloads == (load_data.Length()+constraint_data.Length()) ) return totalnumberofloads;
	else return -1;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// import helper functions

//constant mapping arrays ... which HOTINT node number is the node read from ANSYS file...
// 3D elements:
	const int map_hex[] = {1, 2, 4, 3, 5, 6, 8, 7};
	const int map_prism[] = {1, 2, 3, 3, 5, 6, 7, 7}; // mapped to Hex
	const int map_pyram[] = {1, 2, 4, 3, 5, 5, 5, 5}; // mapped to Hex
	const int map_tet[] = {1, 2, 3, 4};

//	const int map_hexquad[] = {4, 1, 3, 2, 8, 5, 7, 6, 12, 10, 16, 14, 11, 9, 15, 13, 20, 17, 19, 18}; // rotated 90° for some reason... identical to table below
//	const int map_hexquad[] = {2, 4, 3, 1, 6, 8, 7, 5, 14, 10, 13, 9, 16, 12, 15, 11, 18, 20, 19, 17}; // checked inverse of table (as works for tetquad) -> nonsense
  const int map_hexquad[] = {1, 2, 4, 3, 5, 6, 8, 7, 9, 11, 13, 15, 12, 10, 16, 14, 17, 18, 20, 19}; // this would be non-rotated...
//	const int map_prismquad[] = {4, 1, 3, 2, 8, 5, 7, 6, 12, 10, 16, 14, 11, 9, 15, 13, 20, 17, 19, 18}; // mapped to Hexquad
	const int map_prismquad[] = {1, 2, 4, 3, 5, 6, 8, 7, 9, 11, 13, 15, 12, 10, 16, 14, 17, 18, 20, 19}; // this would be non-rotated...
//	const int map_pyramquad[] = {4, 1, 3, 2, 8, 5, 7, 6, 12, 10, 16, 14, 11, 9, 15, 13, 20, 17, 19, 18}; // mapped to Hexquad
	const int map_pyramquad[] = {1, 2, 4, 3, 5, 6, 8, 7, 9, 11, 13, 15, 12, 10, 16, 14, 17, 18, 20, 19}; // this would be non-rotated...
//	const int map_tetquad[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	const int map_tetquad[] = {1, 2, 3, 4, 5, 6, 7, 9, 10, 8}; // this gives a correctly drawn tetquad
//	const int map_tetquad[] = {1, 2, 3, 4, 5, 6, 7, 10, 8, 9}; // this should be correct from table below
	const int map_tetquad20[] = {1, 2, 3, 5, 9, 10, 12, 19, 17, 18}; // mapped to Hexquad

// 2D elements:
	const int map_quad[] = {1, 2, 3, 4};
	const int map_trig[] = {1, 2, 3};
	const int map_quadquad[] = {1, 2, 3, 4, 5, 6, 7, 8};
	const int map_trigquad[] = {1, 2, 3, 4, 5, 6};

// map Node enumeration of Ansys to Hotint, Files generated with Ansys2Hotint.mac macro
int FEMesh::MapNodesAnsys2Hotint(int elementdimension, int elementorder, int nodesused, TArray<int>& nodeshotint)
{
	if (elementdimension != 3) return nodesused;
// Assumptions: 
// *** retvalue: -2: element type unknown | -1: consistency error on multiple nodes | 0: element type not implemented yet | n: length of Node-Array
	TArray<int> nodesansys(nodeshotint); // TArray passes in: Ansys, passes out: Hotint
	nodeshotint.Init();

	if ((elementdimension == 3) && (elementorder == 1) && (nodesused == 8)) //Hex
	{ //e.g.: Mesh200(keyopt=10), Solid45, Solid185, 
		for (int i = 0; i <= 7; i++) nodeshotint.Add(nodesansys(map_hex[i]));
	}
	else if ((elementdimension == 3) && (elementorder == 1) && (nodesused == 6)) // Prism
	{ //
		for (int i = 0; i < 8; i++) nodeshotint.Add(nodesansys(map_prism[i]));
	}
	else if ((elementdimension == 3) && (elementorder == 1) && (nodesused == 5)) // Pyram
	{
		for (int i = 0; i < 8; i++) nodeshotint.Add(nodesansys(map_pyram[i]));
	}
	else if ((elementdimension == 3) && (elementorder == 1) && (nodesused == 4)) // Tet
	{
		for (int i = 0; i < 4; i++) nodeshotint.Add(nodesansys(map_tet[i]));
	}

	else if ((elementdimension == 3) && (elementorder == 2) && (nodesused == 20)) // Hexquad
	{ //e.g.: Solid95, Solid186, 
		for (int i = 0; i < 20; i++) nodeshotint.Add(nodesansys(map_hexquad[i]));
	}	
	else if ((elementdimension == 3) && (elementorder == 2) && (nodesused == 15)) // PrismQuad
	{ // degenerate Hex - e.g.: Solid95, Solid186, 
		for (int i = 0; i < 20; i++) nodeshotint.Add(nodesansys(map_prismquad[i]));
	}
	else if ((elementdimension == 3) && (elementorder == 2) && (nodesused == 13)) // PyramQuad
	{ // degenerate Hex - e.g.: Solid95, Solid186, 
		for (int i = 0; i < 20; i++) nodeshotint.Add(nodesansys(map_pyramquad[i]));
	}
	else if ((elementdimension == 3) && (elementorder == 2) && (nodesused == 10)) // TetQuad
	{ // check if true tetquad or just degenerated hexquad
		if ((nodesansys.Length() == 20) && (nodesansys(11) != 0) && (nodesansys(12) != 0) && (nodesansys(13) != 0) && (nodesansys(14) != 0) && (nodesansys(15) != 0)) //...
		{// degenerate Hex - e.g.: Solid95, Solid186, 
			for (int i = 0; i < 10; i++) nodeshotint.Add(nodesansys(map_tetquad20[i]));
		}
		else 
		{ //true Tet - e.g.: Solid187, 
			for (int i = 0; i < 10; i++) nodeshotint.Add(nodesansys(map_tetquad[i]));
		}
	}
	else if (elementdimension == 2) // do not change sorting order of 2d elements
	{
		nodeshotint = nodesansys;
	}
	else
	{
		nodeshotint.CopyFrom(nodesansys);
		mbs->UO(UO_LVL_err) << "NodeMapping failed! FEType not implemented.\n";
	}
	return nodeshotint.Length();
}

// changes node order (hotint) to be node order (ansys)
int FEMesh::MapNodesHotint2Ansys(int elementdimension, int elementorder, int nodesused, TArray<int>& nodesansys)
{
	if (elementdimension != 3) return nodesused;
// Assumptions: 
// *** retvalue: -2: element type unknown | -1: consistency error on multiple nodes | 0: element type not implemented yet | n: length of Node-Array
	TArray<int> nodeshotint(nodesansys); // TArray passes in: Hotint, passes out: Ansys
	nodesansys.Init();

// get the mapping function Ansys->Hotint
	NaturalNumbers indices(nodesused);
	//TArray<int> indices(nodeshotint);
	//for(int i=1; i<=indices.Length(); i++)
	//{
	//	indices(i) = i;
	//}
	MapNodesAnsys2Hotint(elementdimension, elementorder, nodesused, indices);

	for(int i=1; i<=nodesused; i++)
	{
		nodesansys(indices(i)) = nodeshotint(i);
	}

	return nodesansys.Length();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ forward data to MBS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// forward by entity type

void FEMesh::AddPointsToMBS(int bodyindex)
{
	mbsnodenumbers.SetLen(NN());
	mbsnodenumbers.SetAll(0);

	// initialize searchtree used by mbs->AddBodyNode
	ComputeNodeBox();
	int meshdim = GetMeshDim();
	int3 divs;
	if(meshdim == 3) // meshdim == 3 
	{
		int div = (int) pow(NN()*0.1,1/3.) +1;
		divs.Set(div,div,div);
	}
	else // meshdim == 2
	{
		int div = (int) pow(NN()*0.1,1/2.) +1;
		divs.Set(div,div,1);
	}
	add2mbstree.ResetSearchTree(divs(1),divs(2),divs(3),GetNodeBox());

	// add all nodes
	for (int i = 1; i <= NN(); i++)
	{
		AddSingleNodeToMBS(i, add2mbstree);
	}
}

// adds a single node to the MBS, node dimension is computed from existing elements if not forced
int FEMesh::AddSingleNodeToMBS(int i, SearchTree& addtree, int forced_nodedim)
{
	// node dimension
	int nodedim = forced_nodedim;
	if(!nodedim)
	{
		if(Is3DNode(i)) 
		{
			nodedim = 3;
		}
		else if (Is2DNode(i))
		{
			nodedim = 2;
		}
		else // mixed dimension node - used in 2D and 3D finite elements
		{
			mbs->UO(UO_LVL_warn) << mystr("Warning: Node ") + mystr(i) + mystr(" is used by 3D and 2D elements! use diffent body numbers and constraint! \n");
			nodedim = 3;
		}
	}
	FEMesh_Node& fenode = GetNode(i);
  // mbs node object & initia values
	Node n(nodedim, fenode.Domain(), fenode.GetCoords3D()); // dimension, bodynumber, position
	Vector3D init_d(0.);
	Vector3D init_v(0.);
	if(init_disp.Length()) init_d = GetInitialNodalDisplacement3D(i);
	if(init_vel.Length()) init_v = GetInitialNodalVelocity3D(i);

	if(nodedim == 3) // nodedim == 3 
	{
		n.SetX_Init(init_d, init_v); // initial conditions 3D
	}
	else // nodedim == 2
	{
		n.SetX_Init(init_d.MakeV2D(), init_v.MakeV2D()); // initial conditions 2D
	}
	// add to MBS
	mbsnodenumbers(i)=(mbs->AddBodyNode(&n,addtree));
	return mbsnodenumbers(i);
}

//pass material data from internal array to mbs
void FEMesh::AddMaterialsToMBS()
{
	mbsmaterialnumbers.SetLen(NMaterials());
	mbsmaterialnumbers.SetAll(0);

	for (int i=1; i <= material_data.Length(); i++)
	{
		if(this->GetMeshDim() == 2)
		{
			material_data(i) -> SetType( (TMBSMaterial) (material_data(i)->GetType() | TMat2D) );
		}
		mbsmaterialnumbers(i)=mbs->AddMaterial(*material_data(i));
	}
}

//pass element data from internal array to mbs, !nodes and materials must be already added!
void FEMesh::AddElementsToMBS()
{
	mbselementnumbers.SetLen(NElements());
	mbselementnumbers.SetAll(0);
	SearchTree& addtree = add2mbstree;

	for (int i=1; i<=NElements(); i++)
	{
		FEElement& feelem = GetElement(i);
		//if(feelem.Type() == 2 && feelem.Order() == 2)   //!AD:?Temporary? hack for Plate2Dquad, FEElement 8 nodes -> MBSElement 9 nodes
		//{
		//	FEQuadquad& feqq = (FEQuadquad&) feelem;
		//	Vector2D center = feelem.GetCenterPoint(this->points).MakeV2D();
		//	int additional_node_number = AddPoint(center,feelem.Domain());
		//	int mbs_additional_node_number = AddSingleNodeToMBS(additional_node_number, addtree, 2); // force nodedimension 2
		//	feqq.SurplusNode() = additional_node_number;
		//}
		Element* element_ptr = GetPtrToMBSElement_new(i); //get new element
		if (element_ptr)
		{
			mbselementnumbers(i)=mbs->AddElement(element_ptr);
			delete element_ptr; //delete element!
		}
	}
}

void FEMesh::AddLoadsToMBS()
{
	// this is intended for loads from an external file!!
	// Loads are stored in 3 arrays: 
	// (1) load_type (defines, what the contents other two arrays really mean)
	// (2) number of area to apply to (could be number of element or node too)
	// (3) load vector 3D (direction and strength, could also hold axis and spring stiffness)

// LOADS
	for (int i=1; i <= load_data.Length(); i++)
	{
		//int loadtype,loadtarget;
		//Vector3D loadvector;
//old:		GetLoad(i,loadtype,loadtarget,loadvector); 
//new:
		FEMesh_Load& load = GetLoad(i);

		if(load.Type() == TAreaLoad)
		{
// AreaLoad: external force applied to defined surface (area == array of faces <-> selection of surfaceelements) 
// content of LV := Vector - Load Vector for a unit area
			AddLoadToMBS_TAreaLoad((FEMesh_AreaLoad&) load);
		}

		else if(load.Type() == TFaceLoad)
		{
// FaceLoad: external force applied to a face 
// content of LV := Vector - Load Vector for a unit area
// use for linear elements only
			AddLoadToMBS_TFaceLoad((FEMesh_FaceLoad&) load);
		}

		else if(load.Type() == TBodyLoad) 
		{
// BodyLoad: 
// content of LV:=  x-> strength of load, y-> axis, z->flag (0=bodyload, 1=gravityload)
			AddLoadToMBS_TBodyLoad((FEMesh_BodyLoad&) load);
		}
		else
		{
			mbs->UO(UO_LVL_warn) <<  "Load #" << i << "'s Type is not implemented, Load NOT added to MBS\n";
		}
	}
}

//pass constraint data from internal array to mbs, !elements must be already added!
void FEMesh::AddConstraintsToMBS()
{
// !!! to prevent that constraints are set multiple times, a list is computed before the constraints are added to the MBS !!!

// nodal Constraints
// list containing the corresponding FEMesh_Constaints for all three directions for all nodes
	TArray<int3> nodeconstraintlist;
	nodeconstraintlist.SetLen(NN());
	nodeconstraintlist.SetAll(int3(0,0,0)); 

// spherical joints
// list containing the corresponding FEMesh_Constaints for all nodes
	TArray<TArray<int>*> sphericaljointlist;
  sphericaljointlist.SetLen(NN());
	for(int i=1; i <= NN(); i++)
	{
		sphericaljointlist(i) = new TArray<int>;
	}

// first loop: filter for redundant entries - mark all nodes 
	for(int i=1; i <= constraint_data.Length(); i++)
	{
		//new:
		FEMesh_Constraint& constraint = GetConstraint(i);

		if(constraint.Type() == TAreaConstraint) 
		{
// AreaConstraint: external constraint applied to defined surface (surfaceelement) 
// content of LoadVector:=  x-> constrained axis , y-> springstiffness, z->damping	
			AddConstraintToMBS_TAreaConstraint( (FEMesh_AreaConstraint&) constraint, nodeconstraintlist, i);
		}

		else if (constraint.Type() == TFaceConstraint)
		{
// FaceConstraint: external constraint applied to defined face (surfaceelement) 
// content of LoadVector:=  x-> constrained axis , y-> springstiffness, z->damping	
			AddConstraintToMBS_TFaceConstraint( (FEMesh_FaceConstraint&) constraint, nodeconstraintlist, i);
		}

		else if(constraint.Type() == TNodalConstraint) // NodalConstraint
		{
// NodalConstraint: external constraint applied to defined node 
// content of LV:=  x-> constrained axis , y-> springstiffness, z->damping
			AddConstraintToMBS_TNodalConstraint( (FEMesh_NodeConstraint&) constraint, nodeconstraintlist, i);
		}
		else if (constraint.Type() == TAreaContact)
		{
// TAreaContact: contact between two areas 
// content of LV:=  x-> 2nd area number , y-> springstiffness, z->damping
// all nodes of 1st area (lta) are constrained to a 2nd area (lv.x) with a spherical joint
			AddConstraintToMBS_TAreaContact_Filter( (FEMesh_AreaContact&) constraint, sphericaljointlist, i);
		}
		else
		{
			mbs->UO(UO_LVL_warn) <<  "Constraint #" << i << "'s Type is not implemented, Constraint NOT added to MBS\n";
		}
	}

// Add all Nodal Constraints to MBS
	for(int i=1; i <= NN(); i++)
	{
		for(int j=1; j <= 3; j++)
		{
			if (nodeconstraintlist(i)(j))
			{
				FEMesh_Constraint& c = GetConstraint(nodeconstraintlist(i)(j));   // get FEConstraint for this node j-direction
				AddNodeConstraintToMBS(i, j, c.GetPenalty(), c.Stiffness(), c.Damping(), c.Steps());
			}
		}
	}
// Add all Spherical Joints to MBS
	for(int i=1; i <= constraint_data.Length(); i++)
	{
		//new:
		FEMesh_Constraint& constraint = GetConstraint(i);
	
		if (constraint.Type() == TAreaContact)
		{
// TAreaContact: contact between two areas 
// content of LV:=  x-> 2nd area number , y-> springstiffness, z->damping
// all nodes of 1st area (lta) are constrained to a 2nd area (lv.x) with a spherical joint
			AddConstraintToMBS_TAreaContact( (FEMesh_AreaContact&) constraint, sphericaljointlist, i);

		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper functions - elements

// returns new MBSElement of FEelement "elem_num"
// called by: FEMesh::AddElementsToMBS & CMSElement::AddFEMesh 
// uses settings (changed by YV, 03.11.10)
//$ EK 2012-08-09: type of finite element would have been defined in FeMesh_aux.cpp (TFEType)
Element* FEMesh::GetPtrToMBSElement_new(int elem_num) const
{
// 
	int i = elem_num;

	FEElement* e = elements(i);
	TArray<int> nodelist(20);

	double thickness2D = Settings().Thickness2D();
	if(e->Thickness() > 0.)
		thickness2D = e->Thickness(); // thickness is determined by element itself

//JG 2011-03-21: in case of ffrf elements, nodes are not available!!!
	int matnum;
	if (settings.generate_ffrf_elements) // FEMesh -> FFRF-Element
	{
		for (int j=1; j<=e->NNodes(); j++) // nodenumbers as in FEMesh
		{
			nodelist.Add(e->GetNode(j));
		}
		matnum = e->MaterialNum(); // materials from FEMesh
	}
	else // FEMesh -> MBS
	{
		for (int j=1; j<=e->NNodes(); j++)
		//	nodelist.Add(e->GetNode(j));
		{ 
			int meshnodenum = e->GetNode(j);
			int mbsnodenum = GetMBSNodeNumber(meshnodenum);
			nodelist.Add(mbsnodenum); // nodenumbers as in MBS
		}
		//if(e->Type()==2 && e->Order()==2) //!AD:?Temporary? hack for Plate2Dquad, FEElement 8 nodes -> MBSElement 9 nodes
		//{
		//	int meshnodenum = ((FEQuadquad*) e)->SurplusNode();
		//	int mbsnodenum = GetMBSNodeNumber(meshnodenum);
		//	nodelist.Add(mbsnodenum); // nodenumbers as in MBS
		//}
		matnum = GetMBSMaterialNumber(e->MaterialNum()); // materialnumbers as in MBS
	}

	int domain = settings.domain;
// domain == -1 -> use from element
	if(domain == -1)
		domain = e->Domain();
	if(e->Color()==Vector3D(-1.0,0.0,0.0)) //$ LA 2010-10: default color -RED -1.0,0.0,0.0, replace with material color
	{
		if (material_data.Length() >= e->MaterialNum()) // check if materialdata exists
			e->Color() = material_data(e->MaterialNum())->GetMaterialColor(); // overwrite with material.color
	}

//3D - common base class FEElement3D (base::Setbase function)
// creating 3D finite elements - rewritten by YV, 03.11.2010
// $EK 2013-03-04 extended to TFEPrism and TFEPyramid
	if ((e->Type() ==	TFETet || e->Type() == TFEHex || e->Type() == TFEPrism || e->Type() == TFEPyramid) && (e->Order() <= 2))
	{
		FiniteElement3D * fe3D;
		if(e->Type() ==	TFETet)
		{
			if(settings.generate_ffrf_elements)
			{
				TetrahedralFFRF * p = new TetrahedralFFRF(mbs);
				p->SetTetrahedral(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				//$ LA 2011-2-18: Geometric nonlinearity status is not set for FFRF elements, since it is already set to GNS_Linear within FiniteElement3DFFRF::SetFFRF(int isffrf)
			}
			else
			{
				Tetrahedral * p = new Tetrahedral(mbs);
				p->SetTetrahedral(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				fe3D->SetGeometricNonlinearityStatus(Settings().FEMesh_geometricNonlinearityStatus); //$ LA 2011-2-18: Sets geometric nonlinearity status as defined in FEMesh settings (default GNS_NonlinearLargeStrain)
			}
		}
		//$EK 2013-03-04
		else if(e->Type() == TFEPrism)
		{
			if(settings.generate_ffrf_elements)
			{
				mbs->UO(UO_LVL_err) << "**error: FFRF not implemented so far for Prismatic Elements!\n";
			}
			else
			{
				Prism * p = new Prism(mbs);
				p->SetPrism(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				fe3D->SetGeometricNonlinearityStatus(Settings().FEMesh_geometricNonlinearityStatus); //$ LA 2011-2-18: Sets geometric nonlinearity status as defined in FEMesh settings (default GNS_NonlinearLargeStrain)
			}
		}
		//$EK 2013-03-04
		else if(e->Type() == TFEPyramid)
		{
			if(settings.generate_ffrf_elements)
			{
				mbs->UO(UO_LVL_err) << "**error: FFRF not implemented so far for Pyramids!\n";
			}
			else
			{
				Pyramid * p = new Pyramid(mbs);
				p->SetPyramid(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				fe3D->SetGeometricNonlinearityStatus(Settings().FEMesh_geometricNonlinearityStatus); //$ LA 2011-2-18: Sets geometric nonlinearity status as defined in FEMesh settings (default GNS_NonlinearLargeStrain)
			}
		}
		else
		{
			if(settings.generate_ffrf_elements)
			{
				HexahedralFFRF * p = new HexahedralFFRF(mbs);
				p->SetHexahedral(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				//$ LA 2011-2-18: Geometric nonlinearity status is not set for FFRF elements, since it is already set to GNS_Linear within FiniteElement3DFFRF::SetFFRF(int isffrf)
			}
			else
			{
				Hexahedral * p = new Hexahedral(mbs);
				p->SetHexahedral(domain, nodelist, matnum, e->Color(), settings.cms_element_number);
				fe3D = p;
				fe3D->SetGeometricNonlinearityStatus(Settings().FEMesh_geometricNonlinearityStatus); //$ LA 2011-2-18: Sets geometric nonlinearity status as defined in FEMesh settings (default GNS_NonlinearLargeStrain)
			}
		}
		return fe3D;
	}
	//else if ((e->Type() == 3) && (e->Order() == 2))
	//{
	//	elem = new Tetrahedral(mbs);
	//	tet.SetTetrahedral(domain, nodelist, matnum, e->Color(),CMSElementNr);
	//	mbselementnumbers(i)=mbs->AddElement(&tet);
	//}
	//else if ((e->Type() == 4) && (e->Order() == 2))
	//{
	//	elem = new Hexahedral(mbs);
	//	hex.SetHexahedral(domain, nodelist, matnum, e->Color(),CMSElementNr);
	//	mbselementnumbers(i)=mbs->AddElement(&hex);
	//}
	//2D - common base class Plate2D (base::Setbase function)

// creation of the elements below should be updated to the new 2D finite elements, YV

	else if ((e->Type() == TFEQuad) && (e->Order() == 1))
	{
		if(settings.flag_use_plate2D == 0)
		{
			Quadrilateral* quad = new Quadrilateral(mbs);
			quad->SetQuadrilateral(domain, nodelist, matnum, thickness2D, e->Color());
			return quad;
		}
		else if(settings.flag_use_plate2D == 1)
		{
			Plate2Dlin* quad = new Plate2Dlin(mbs);
			quad->SetPlate2Dlin(domain, nodelist, matnum, thickness2D, e->Color(), settings.cms_element_number);
			return quad;
		}
	}
	else if ((e->Type() == TFETrig) && (e->Order() == 1))
	{
		Trig2Dlin* trig = new Trig2Dlin(mbs);
		trig->SetTrig2Dlin(domain, nodelist, matnum, thickness2D, e->Color(), settings.cms_element_number);
		return trig;
	}
	else if ((e->Type() == TFEQuad) && (e->Order() == 2))
	{
		if(settings.flag_use_plate2D == 0)
		{
			Quadrilateral* quadquad = new Quadrilateral(mbs);
			quadquad->SetQuadrilateral(domain, nodelist, matnum, thickness2D, e->Color());
			return quadquad;
		}
		else if (settings.flag_use_plate2D == 1)
		{
			Plate2Dquad* quadquad = new Plate2Dquad(mbs);
			quadquad->SetPlate2Dquad(domain, nodelist, matnum, thickness2D, e->Color(), settings.cms_element_number);
			return quadquad;
		}
	}
	else if ((e->Type() == TFETrig) && (e->Order() == 2))
	{
		Trig2Dquad* trigquad = new Trig2Dquad(mbs);
		trigquad->SetTrig2Dquad(domain, nodelist, matnum, thickness2D, e->Color(), settings.cms_element_number);
		return trigquad;
	}
  if(Settings().MBSElementWarnings())
	{
		mbs->UO(UO_LVL_err) << "Error in FEMesh::GetPtrToMBSElement_new, could not convert FEElement\n";
	}
	return 0; //should not happen
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper functions - loads

// Processes TAreaLoad type Load and adds to MBS
void FEMesh::AddLoadToMBS_TAreaLoad(FEMesh_AreaLoad& load)
{
	int areanumber = load.GetElement();
	Vector3D loadvector = load.GetLoadVector();
	TArray<StepSettings> loadsteps = load.GetSteps();
	int flag_pressure = load.GetFlagPressure();
// load is applied to the 3 or 4 nodes of the element
// use for linear elements only
	IVector& facelist = GetArea(areanumber).Faces();
	TArray<double> node_weight; node_weight.SetLen(NN()); node_weight.SetAll(0.);
	TArray<int> elemnrs; elemnrs.SetLen(NN()); elemnrs.SetAll(0);
	TArray<int> local_node_nrs; local_node_nrs.SetLen(NN()); local_node_nrs.SetAll(0);

	double area = 0.;
	if (flag_pressure == 0)
	{
		for (int i = 1; i <= facelist.Length(); ++i)
			area += GetFaceArea(GetFace(facelist.Get(i)));
	}
	else	
		area = 1.;

	ComputeNodeWeight(facelist, node_weight, elemnrs, local_node_nrs);
	for(int i=1; i <= NN(); i++)
	{
		if(node_weight(i) > 0.)
		{
			AddNodalLoadToMBS(i, loadvector*(node_weight(i)/area), elemnrs(i), local_node_nrs(i), loadsteps);
		}
	}
}

// Processes TFaceLoad type Load and adds to MBS
void FEMesh::AddLoadToMBS_TFaceLoad(FEMesh_FaceLoad& load)
{
	int facenumber = load.GetElement();
	Vector3D loadvector = load.GetLoadVector();
	TArray<StepSettings> loadsteps = load.GetSteps();
// load is applied to the 3 or 4 nodes of the element
// use for linear elements only
	TArray<int> facelist(0);
	facelist.Add(facenumber);
	TArray<double> nodes_weight; nodes_weight.SetLen(NN()); nodes_weight.SetAll(0.);
	TArray<int> elemnrs; elemnrs.SetLen(NN()); elemnrs.SetAll(0);
	TArray<int> local_node_nrs; local_node_nrs.SetLen(NN()); local_node_nrs.SetAll(0);

	ComputeNodeWeight(facelist, nodes_weight, elemnrs, local_node_nrs);

	FEMesh_Set facenodes(TSetNodes);
	GetNodesOfFaces(facenodes, facenumber);
	
	for(int j=1; j <= facenodes.NNodes(); j++)
	{
		int i = facenodes(j);
 		if(nodes_weight(i) >= 0.)
		{
			AddNodalLoadToMBS(i, loadvector*nodes_weight(i), elemnrs(i), local_node_nrs(i), loadsteps);
		}
	}
}

// Processes TBodyLoad type Load and adds to MBS
void FEMesh::AddLoadToMBS_TBodyLoad(FEMesh_BodyLoad& load)
{
	int bodynumber = load.GetElement();
	Vector3D loadvector = load.GetLoadVector();
	TArray<StepSettings> loadsteps = load.GetSteps();
	int flag_gravity = load.GetGravity();
// applies a load to a single body
	MBSLoad mbs_load;
	mbs_load.SetLoadSteps(loadsteps);
	if (loadvector.X() != 0.0)
	{
		if(!flag_gravity) mbs_load.SetBodyLoad(loadvector.X(), 1);
		else mbs_load.SetGravity(loadvector.X(), 1);
		mbs->GetElement(mbselementnumbers(bodynumber)).AddLoad(mbs_load);
	}
	if (loadvector.Y() != 0.0)
	{
		if(!flag_gravity) mbs_load.SetBodyLoad(loadvector.Y(), 2);
		else mbs_load.SetGravity(loadvector.Y(), 2);
		mbs->GetElement(mbselementnumbers(bodynumber)).AddLoad(mbs_load);
	}
	if (loadvector.Z() != 0.0)
	{
		if(!flag_gravity) mbs_load.SetBodyLoad(loadvector.Z(), 3);
		else mbs_load.SetGravity(loadvector.Z(), 3);
		mbs->GetElement(mbselementnumbers(bodynumber)).AddLoad(mbs_load);
	}
}

// computes the relative weight of each node in a TAreaLoad or TFaceLoad constraint, saves (1 possible) element and sidenumber for each node
void FEMesh::ComputeNodeWeight(TArray<int>& facelist, TArray<double>& nodes_weight, TArray<int>& elemnrs, TArray<int>& local_node_nr)
{
	for(int i=1; i <= facelist.Length(); i++) // compute weight for all nodes of all faces
	{
		int elemnr = GetElementOfFace(facelist(i));
		int sidenr = GetElementSideOfFace(facelist(i));
		int elemorder = GetElement(elemnr).Order();

		double area = GetElementFaceArea(elemnr, sidenr);
		double elemorder_factor = 1./elemorder;

		int4 nodes;
		GetElement(elemnr).GetSideNodeNum(sidenr,nodes); // global node numbers
		int len = nodes.RemoveAllRedundantEntries(); // redundant
		
		int4 nodes_side;
		GetElement(elemnr).GetSide(sidenr,nodes_side); // local node numbers - does not recognize degenerated hex (prisms & pyramids)
		
		int4 nodes_local;
		for(int j=1; j <= len; j++)      // create a mapping to local node numbers (local DOFs)
		{
			for(int k=1; k <= 4; k++)
			{
				int ng = nodes(j);
				int nl = GetElement(elemnr).GetNode(nodes_side(k)); 
				if ( ng == nl )
				{
					nodes_local(j) = nodes_side(k);
				}
			}
		}
		Vector3D center;
		if (elemorder == 2)
			center = GetFaceCenterPoint(nodes);
		//double sum_of_interior_angles = MY_PI * (len-2.);    //(180.0 * (len-2.));

		for(int j=1; j <= len; j++) // individual node - remember element and local node number, add weight
		{
// angular factor - angle at a node of the side
			////double interior_angle = ComputeAngleAtNode(j, nodes);
			////int node_global = nodes(j);
			////int node_local = nodes_local(j);

			////elemnrs(node_global) = elemnr;
			////local_node_nr(node_global) = node_local;
			////nodes_weight(node_global) += (area * interior_angle) / sum_of_interior_angles; // add weight to global node number

// area factor - area of deltoid 
			double nodearea = ComputeNodeArea(j, nodes);
			int node_global = nodes(j);
			int node_local = nodes_local(j);

			elemnrs(node_global) = elemnr;
			local_node_nr(node_global) = node_local;
			// $EK 2013-04-03 for quadratic elements intermediate points have to be weigthed
			//nodes_weight(node_global) += nodearea ; // add weight to global node number
			nodes_weight(node_global) += elemorder_factor*nodearea ; // add weight to global node number
			if (elemorder == 2)
			{
				int internode_loc = GetElement(elemnr).GetIntermediateNode(nodes_local(j), nodes_local(j%len+1));
				int internode = GetElement(elemnr).GetNode(internode_loc);
				if (!internode) 
				{
					mbs->UO(UO_LVL_warn) << "**** WARNING: GetIntermediateNodeNum not implemented for element!\n";
					continue;
				}
				Vector3D p1 = GetPoint3D(nodes(j));
				Vector3D p2 = GetPoint3D(internode);
				Vector3D p3 = GetPoint3D(nodes(j%len+1));
				p1 = (p1 + p2)*0.5;
				p3 = (p2 + p3)*0.5;
				nodes_weight(internode) += fabs(::Area(center, p1, p3));
				elemnrs(internode) = elemnr;
				local_node_nr(internode) = internode_loc;
			}
		}
	}
}

// adds the nodal load to the MBS 
void FEMesh::AddNodalLoadToMBS(int nodenumber, Vector3D& loadvector, int elemnr, int local_node_nr, TArray<StepSettings> & loadsteps)
{
	MBSLoad load;
	if (loadvector.X() != 0.0)
	{
		load.SetGCLoad(loadvector.X(), (local_node_nr-1)*3 +1 );
		load.SetLoadSteps(loadsteps);
		mbs->GetElement(elemnr).AddLoad(load);
	}
	if (loadvector.Y() != 0.0)
	{
		load.SetGCLoad(loadvector.Y(), (local_node_nr-1)*3 +2 );
		load.SetLoadSteps(loadsteps);
		mbs->GetElement(elemnr).AddLoad(load);
	}
	if (loadvector.Z() != 0.0)
	{
		load.SetGCLoad(loadvector.Z(), (local_node_nr-1)*3 +3 );
		load.SetLoadSteps(loadsteps);
		mbs->GetElement(elemnr).AddLoad(load);
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper functions - constraints

// Processes TAreaConstraint type Constraint and fill nodeconstraint list
void FEMesh::AddConstraintToMBS_TAreaConstraint(FEMesh_AreaConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr)
{
	int areanumber = constraint.GetElement();
	int direction = constraint.GetAxis();
	int penalty = constraint.GetPenalty();
	Vector3D stiffness = constraint.GetStiffness();
	Vector3D damping = constraint.GetDamping();
// constraints are applied to all nodes of the defined area
	FEMesh_Set nodelist(TSetNodes);
	GetNodesOfArea(nodelist,areanumber);
	::RemoveRedundantEntries(nodelist.Nodes());

	for(int i=1; i <= nodelist.NNodes(); i++)
	{
		nodeconstraintlist(nodelist(i))(direction) = nr; // enter all nodal constraints of this area constraint to list
	}
}

// Processes TFaceConstraint type Constraint and fill nodeconstraint list
void FEMesh::AddConstraintToMBS_TFaceConstraint(FEMesh_FaceConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr)
{
	int facenumber = constraint.GetElement();
	int direction = constraint.GetAxis();
	int penalty = constraint.GetPenalty();
	Vector3D stiffness = constraint.GetStiffness();
	Vector3D damping = constraint.GetDamping();
// constraint applied to all nodes of a single face

	FEMesh_Set facenodes(TSetFaces);
	GetNodesOfFaces(facenodes, facenumber);
	
	for(int i=1; i <= facenodes.NNodes(); i++)
	{
		nodeconstraintlist(facenodes(i))(direction) = nr; // enter all nodal constraints of this face constraint to list
	}
}

// Processes TNodalConstraint type Constraint and fill nodeconstraint list
void FEMesh::AddConstraintToMBS_TNodalConstraint(FEMesh_NodeConstraint& constraint, TArray<int3>& nodeconstraintlist, int nr)
{
	int nodenumber = constraint.GetElement();
	int direction = constraint.GetAxis();
	int penalty = constraint.GetPenalty();	
	Vector3D stiffness = constraint.GetStiffness();
	Vector3D damping = constraint.GetDamping();
// constraint applied to a single node		
// some variables for the constructor that are not in the current "load vector"
// will be included in FEConstraint class ...

	nodeconstraintlist(nodenumber)(direction) = nr; // enter nodal constraint to list
}

// adds the nodal constraint to the MBS
void FEMesh::AddNodeConstraintToMBS(int nodenumber, int direction, int penaltyflag, Vector3D& stiffness, Vector3D& damping, TArray<StepSettings>& steps)
{
	double ground = 0.;
	int velocityconstraint = 0;

	NodalConstraint nc(mbs);
	if( (direction == 1) || (direction == 2) || (direction == 3)) // valid directions 1=x, 2=y, 3=z
	{
		if(!penaltyflag)
		{
			nc.SetNodalConstraint(nodenumber, direction, ground);// lagrange  //, drawsize);  //$ LA 2011-4-22: drawsize uses global value
		}
		else
		{
			nc.SetNodalConstraintSpringDamper(nodenumber, direction, ground, -1., stiffness.X(), damping.X()); // penalty //$ LA 2011-4-22: drawsize uses global value
		}
		nc.SetConstraintSteps(steps);
		mbs->AddElement(&nc);		
	}
	else 
	{
		mbs->UO(UO_LVL_warn) <<  "Constraint's direction is invalid! Set 1=x, 2=y or 3=z. Constraint NOT added to MBS\n";
	}
}

// Processes TAreaContact type Constraint and fill spherical joint list
void FEMesh::AddConstraintToMBS_TAreaContact_Filter(FEMesh_AreaContact& contact, TArray<TArray<int>*>& sphericaljointlist, int nr)
{
	int sendingarea = contact.GetElement();
	int targetarea = contact.GetTarget();
	int penaltyflag = contact.GetPenalty();
	Vector3D stiffness = contact.GetStiffness();
	Vector3D damping = contact.GetDamping();
// all nodes of "sending" area  are constrained to a "target" area with a spherical joint

// list of nodes
	FEMesh_Set nodes_area1(TSetNodes);
	GetNodesOfArea(nodes_area1, sendingarea);

// mark the node with number of Constraint
	for(int i=1; i <= nodes_area1.NNodes(); i++)
	{
		sphericaljointlist(nodes_area1(i))->Add(nr);
	}
}

// Processes TAreaConstraint type Constraint and adds to MBS
void FEMesh::AddConstraintToMBS_TAreaContact(FEMesh_AreaContact& contact, TArray<TArray<int>*>& sphericaljointlist, int nr)
{
	int sendingarea = contact.GetElement();
	int targetarea = contact.GetTarget();
	int penaltyflag = contact.GetPenalty();
	Vector3D stiffness = contact.GetStiffness();
	Vector3D damping = contact.GetDamping();
	TArray<StepSettings> steps = contact.GetSteps();
// all nodes of "sending" area  are constrained to a "target" area with a spherical joint
// it is assumed that the nodes are exactly ON face elements of targetarea

// list of nodes
	FEMesh_Set nodes_area1(TSetNodes),  nodes_area2(TSetNodes);
	GetNodesOfArea(nodes_area1, sendingarea);
	GetNodesOfArea(nodes_area2, sendingarea);

// all global positions 
	TArray<Vector3D> global_positions;
	for(int i=1; i <= nodes_area1.NNodes(); i++)
	{
		for(int j=1; j <= sphericaljointlist(nodes_area1(i))->Length(); j++)
		{
// filter out nodes that are associated with other area contacts
			if(sphericaljointlist(nodes_area1(i))->Get(j) == nr)
			{
				global_positions.Add(GetPoint3D(nodes_area1(i)));
			}
		}
	}

// all elements of area1 - via faces
	FEMesh_Set elems_area1(TSetElements);
	FEMesh_Set faces_area1(TSetFaces);
	faces_area1 = GetArea(sendingarea);
	for(int i=1; i <= faces_area1.NFaces(); i++)
	{
		elems_area1.Add(GetElementOfFace(faces_area1(i)));
	}
	::RemoveRedundantEntries(elems_area1.Elements(),1);
// all elements of area2
	FEMesh_Set elems_area2(TSetElements);
	FEMesh_Set faces_area2(TSetFaces);
	faces_area2 = GetArea(targetarea);
	for(int i=1; i <= faces_area2.NFaces(); i++)
	{
		elems_area2.Add(GetElementOfFace(faces_area2(i)));
	}
	::RemoveRedundantEntries(elems_area2.Elements(),1);

// ELNR1, LC1
	TArray<int> elnr1(global_positions.Length());
	TArray<Vector3D> lc1(global_positions.Length());
	for(int i=1; i <= global_positions.Length(); i++)
	{
		elnr1(i) = GetElementAndLocalCoordOfGlobalPosition(global_positions(i),lc1(i),elems_area1.Elements());
	}
// ELNR2, LC2
	TArray<int> elnr2(global_positions.Length());
	TArray<Vector3D> lc2(global_positions.Length());
	for(int i=1; i <= global_positions.Length(); i++)
	{
		elnr2(i) = GetElementAndLocalCoordOfGlobalPosition(global_positions(i),lc2(i),elems_area2.Elements());
	}

	int rv_elemnr;  // debug only
// create spherical joints
	for(int i=1; i <= global_positions.Length(); i++)
	{
		if( elnr2(i) != 0 )
		{
			if(penaltyflag)
			{
				SphericalJoint sj(mbs,elnr1(i),elnr2(i),lc1(i),lc2(i),stiffness,-1,Vector3D(1.,1.,1.));
				sj.SetConstraintSteps(steps);
				rv_elemnr = mbs->AddElement(&sj);
			}
			else
			{
				SphericalJoint sj(mbs,elnr1(i),elnr2(i),lc1(i),lc2(i),-1,Vector3D(1.,1.,1.));
				sj.SetConstraintSteps(steps);
				rv_elemnr = mbs->AddElement(&sj);	
			}
			
		}	
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ export functions - export in specified format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// export surface in STL format - can be read to GeomMesh3D object

void FEMesh::ExportSurfaceTrigsToSTL(char* file, int binary)
{
	TArray<Vector3D> trigs;
	SurfaceTrigsAsPointSequenceTArray(trigs);
	int n = trigs.Length()/3;

	if (!binary)
	{
		ofstream fout(file);

		fout << "solid\n";
		for (int i=1; i <= n; i++)
		{
			Vector3D normal;
			Vector3D p[3];
			p[0] = trigs(i*3-2);
			p[1] = trigs(i*3-1);
			p[2] = trigs(i*3-0);
			Normal3D(p[0],p[1],p[2], normal);
			fout << "facet ";
			fout << "normal " << normal.X() << " " << normal.Y() << " " << normal.Z() << "\n";
			fout << "outer loop\n";
			for (int j=0; j < 3; j++)
			{
				fout << "vertex " << p[j].X() << " " << p[j].Y() << " " << p[j].Z() << "\n";				
			}
			fout << "endloop\n";
			fout << "endfacet\n";
		}
		fout << "endsolid\n";
	}
	else
	{
		CMFile outfile(file, TFMwrite, 1);
		if (sizeof(int) != 4 || sizeof(float) != 4) mbs->UO(UO_LVL_err) << mystr("for stl-binary compatibility only use 32 bit compilation!!!\n");

		//specific settings for stl-binary format - compare GeomMesh3D::ReadSTLMesh
		const int namelen = 80; //length of name of header in file
		const int nospaces = 2; //number of spaces after a triangle

		char buf[namelen+1]; buf[namelen]='\n';
		mystr str(buf);
		//str.SetLen(80);
		//str.SetAll(' ');

		//write leading 80 characters
		outfile.RWchars(namelen, str);
		//write 4-byte int number of facets
		outfile.RWbinaryInt(n);

		float f;
		char spaces[nospaces];

		//write triangles:
		int i, j;
		for (i = 1; i <= n; i++)
		{
			Vector3D normal;
			Vector3D p[3];
			p[0] = trigs(i*3-2);
			p[1] = trigs(i*3-1);
			p[2] = trigs(i*3-0);
			Normal3D(p[0],p[1],p[2], normal);

			float f;
			//write 3 floats for triangle normal
			f = normal.X(); outfile.RWbinaryFloat(f);
			f = normal.Y(); outfile.RWbinaryFloat(f);
			f = normal.Z(); outfile.RWbinaryFloat(f);

			//write 3 x 3 floats for (x,y,z) coordinates of triangle points
			for (j = 0; j < 3; j++)
			{
				f = p[j].X(); outfile.RWbinaryFloat(f);
				f = p[j].Y(); outfile.RWbinaryFloat(f);
				f = p[j].Z(); outfile.RWbinaryFloat(f);
			} 

			//write two bytes for attributes (usually unused)
			outfile.RWchar(spaces[0]); 
			outfile.RWchar(spaces[1]);
		}
	}
}

// create a list of all surfaceelement as trigs ( nodenumbers )
int FEMesh::SurfaceAsTrigs(TArray<int3>& surfacetrigs)
{
	surfacetrigs.Flush();
	if(Surface().N() == 0) { mbs->UO(UO_LVL_warn) << "No surfaceelements identified, List of surface trigs empty !\n"; return 0; }

	for(int i=1; i <= Surface().N(); i++)
	{
		FEMesh_Face& feface = GetFace(Surface()(i));
		if(feface.NNodes() == 4) // quad
		{
			surfacetrigs.Add(int3(feface.GetNode(1),feface.GetNode(2),feface.GetNode(3)));
			surfacetrigs.Add(int3(feface.GetNode(1),feface.GetNode(3),feface.GetNode(4)));
		}
		else if(feface.NNodes() == 3) // trig
		{
			surfacetrigs.Add(int3(feface.GetNode(1),feface.GetNode(2),feface.GetNode(3)));
		}	
	}
	return surfacetrigs.Length();
}

// creates a list of surface trigs defined by 3 consecutive vector3d
int FEMesh::SurfaceTrigsAsPointSequenceTArray(TArray<Vector3D>& pointsequence)
{
	TArray<int3> surfacetrigs;
	SurfaceAsTrigs(surfacetrigs);

	pointsequence.Flush();
	for (int i=1; i <= surfacetrigs.Length(); i++)
	{
		int3 strig = surfacetrigs(i);
		pointsequence.Add( GetPoint3D(strig(1)) );
		pointsequence.Add( GetPoint3D(strig(2)) );
		pointsequence.Add( GetPoint3D(strig(3)) );
	}
	return pointsequence.Length();
}

// creates an array of surface trigs defined by 3 consecutive vector3d
void FEMesh::SurfaceTrigsAsPointSequence(Vector3D** points, int& len)
{
	TArray<int3> surfacetrigs;
	SurfaceAsTrigs(surfacetrigs);

	int ntrigs = surfacetrigs.Length();
	len = 3 * ntrigs;

	Vector3D* sequence = new Vector3D[len];

	for (int i=1; i <= ntrigs; i++)
	{
		int3 strig = surfacetrigs(i);
		sequence[3*i-3] = GetPoint3D(strig(1));
		sequence[3*i-2] = GetPoint3D(strig(2));
		sequence[3*i-1] = GetPoint3D(strig(3));
	}
	*points = sequence;
	return;
}

// creates an array of surface trigs defined by 9 consecutive doubles (3 points x 3 coordinates)
void FEMesh::SurfaceTrigsAsDoubleSequence(double** doubles, int& len)
{
	TArray<int3> surfacetrigs;
	SurfaceAsTrigs(surfacetrigs);
	int ntrigs = surfacetrigs.Length();
	len = 9 * ntrigs;

	double* sequence = new double[len];

	for (int i=1; i <= ntrigs; i++)
	{
		int3 strig = surfacetrigs(i);
		for (int j=1; j<= 3; j++)
		{
			Vector3D pt = GetPoint3D(strig(j));
			//9*(i-1)+3*(j-1)+a == 9i+3j-12+a
			sequence[9*i+3*j-12] = pt.X();
			sequence[9*i+3*j-11] = pt.Y();
			sequence[9*i+3*j-10] = pt.Z();
		}
	}
	*doubles = sequence;
	return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// various export formats

//$ YV 2011-03-17: Export of elements, nodes and material properties to a simple text file
// the format assumes that all elements have equal number of nodes;
// the function returns true on success;
// file format:
// first line contains (divided by spaces):
// 1. dimensionality of the model (2 or 3);
// 2. number of nodes;
// 3. number of nodes per element;
// 4. number of elements;
// 5. number of materials;
// then follow the definitions of nodes by 2 or 3 coordinates (according to the dimensionality of the model);
// then follow the definitions of elements: 1-based numbers of elements and in the end 1-based number of material;
// then follow the definitions of materials: E, nu, rho
bool FEMesh::ExportAsText(const char * filename)
{
	if(NElements() == 0)
		return false;
	const FEElement & e = GetElement(1);
	int nNodesPerElement = e.NNodes();
	ofstream f;
	f.open(filename);
	if(!f.is_open())
		return false;
	if(GetMeshDim() != 2 && GetMeshDim() != 3)
		return false;
	f << GetMeshDim() << " " << NN() << " " << nNodesPerElement << " " << NElements() << " " << NMaterials() << "\n";
	for(int i = 1; i <= NN(); i++)
	{
		if(GetMeshDim() == 2)
			f << GetPoint2D(i).X() << " " << GetPoint2D(i).Y();
		else
			f << GetPoint3D(i).X() << " " << GetPoint3D(i).Y() << " " << GetPoint3D(i).Z();
		f << "\n";
	}
	for(int i = 1; i <= NElements(); i++)
	{
		const FEElement & e = GetElement(i);
		if(e.NNodes() != nNodesPerElement)
			return false;
		for(int j = 1; j <= nNodesPerElement; j++)
			f << e.GetNode(j) << " ";
		f << e.MaterialNum() << "\n";
	}
	for(int i = 1; i <= NMaterials(); i++)
	{
		const Material & m = GetMaterial(i);
		f << m.YoungsModulus() << " " << m.PoissonRatio() << " " << m.Density() << "\n";
	}

	return true;
}

//$ AD 2012-06-01: Export nodes and elements
bool FEMesh::ExportNodesAndElements(const char * filename, const mystr &format)
{
	CMFile outfile(filename,TFMwrite);

	// NODES
	outfile.RWSoleStr(mystr("[NODES]\n") + mystr(NNodes()));
	for(int i=1; i<=NNodes(); i++)
	{
		FEMesh_Node& the_node = GetNode(i);
		mystr line = mystr(i) + mystr("\t") + mystr(the_node.X()) + mystr("\t") + mystr(the_node.Y()) + mystr("\t") +  mystr(the_node.Z());
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

	// ELEMENTS - differences in node order for ANSYS and HOTINT
	int maxnnr = 0;
	for(int i=1; i<=NElements(); i++)
	{
		if(GetElement(i).NNodes() > maxnnr)
			maxnnr = GetElement(i).NNodes();
	}
	outfile.RWSoleStr(mystr("[ELEMENTS]\n") + mystr(NElements()) + mystr("\t") + mystr(maxnnr));
	for(int i=1; i<=NElements(); i++)
	{
		FEElement& the_element = GetElement(i);
		mystr line = mystr(i) + mystr("\t") + mystr(the_element.Dim())  + mystr("\t") + mystr(the_element.Order()) + mystr("\t") + mystr(the_element.NNodes()) + mystr("\t");

		TArray<int> nodenumbers(20);
		for(int j=1; j<=maxnnr; j++)
		{
			nodenumbers(j) = the_element.GetNode(j);
		}

		if( format == mystr("ANSYS") )
		{
			MapNodesHotint2Ansys(the_element.Dim(), the_element.Order(), the_element.NNodes(), nodenumbers);
		}

		for(int j=1; j<=maxnnr; j++)
		{
			if(j<=the_element.NNodes())
			{
				line = line + mystr(nodenumbers(j)) + mystr("\t");
			}
			else
			{
				line = line + mystr("0\t");
			}
		}
		line = line + mystr(the_element.MaterialNum());
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

// FACES
	outfile.RWSoleStr(mystr("[FACES]\n") + mystr(NFaces()) + mystr("\t") + mystr(maxnnr));
	for(int i=1; i<=NFaces(); i++)
	{
		FEMesh_Face& the_face = GetFace(i);
		FEElement& the_element = GetElement(the_face.GetElement());
		mystr line = mystr(i) + mystr("\t") + mystr(the_element.Dim()-1)  + mystr("\t") + mystr(the_element.Order()) + mystr("\t") + mystr(the_face.NNodes()) + mystr("\t");
		for(int j=1; j<=maxnnr; j++)
		{
			if(j<=the_face.NNodes())
			{
				line = line + mystr(the_face.GetNode(j)) + mystr("\t");
			}
			else
			{
				line = line + mystr("0\t");
			}
		}
		line = line + mystr(the_element.MaterialNum());
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

// ELEMENTS MATERIAL
	outfile.RWSoleStr(mystr("[VMAT]\n") + mystr(NMaterials()));
	for(int i=1; i<=NMaterials(); i++)
	{
		mystr line = mystr(i) + mystr("\t") + mystr(GetMaterial(i).GetMaterialName());
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

// AREA FACES
	outfile.RWSoleStr(mystr("[AELEM]\n") + mystr(NAreas()));
	for(int i=1; i<=NAreas(); i++)
	{
		FEMesh_Set the_area = GetArea(i);
		mystr line = mystr(i) + mystr("\t") + mystr(the_area.NFaces());
		for(int j=1; j<=the_area.NFaces(); j++)
		{
			line = line + mystr("\t") + mystr(the_area.Face(j)) ;
		}
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

// AREA MATERIAL
	outfile.RWSoleStr(mystr("[AMAT]\n") + mystr(NAreas()));
	for(int i=1; i<=NAreas(); i++)
	{
		mystr line = mystr(i) + mystr("\t") + mystr(GetArea(i).Name());
		outfile.RWSoleStr(line);
	}
	outfile.RWSoleStr(mystr("\n"));

// AREA LOADS
	outfile.RWSoleStr(mystr("[LOADS]\n") + mystr(NConstraints()+NLoads()) + mystr("\t{N,m}"));
	for(int i=1; i<=NConstraints(); i++)
	{
		FEMesh_Constraint& the_constraint = GetConstraint(i);	
		mystr line;
		if(the_constraint.Type() == TAreaConstraint) 
		{
			line = mystr("AreaConstraint\t"); 
		}
		else if (the_constraint.Type() == TFaceConstraint)
		{
			line = mystr("FaceConstraint\t");	
		}
		else if(the_constraint.Type() == TNodalConstraint) 
		{
			line = mystr("NodalConstraint\t");	
		}
		else if (the_constraint.Type() == TAreaContact)
		{
			line = mystr("AreaContact\t");	
		}
		else line = mystr("Contact\t");
		line = line + mystr(the_constraint.GetElement()) + mystr("\t") + mystr(the_constraint.GetAxis()) + mystr("\t");
		if(the_constraint.GetPenalty() == 1)
		{
			line = line +	mystr(the_constraint.GetStiffness().X()) + mystr("\t") + mystr(the_constraint.GetDamping().X()) + mystr("\t");
		}
		else
		{
			line = line + mystr("0\t0\t");
		}
		line = line + mystr(the_constraint.GetSteps().Length());
		for(int j=1; j<=the_constraint.GetSteps().Length(); j++)
		{
			StepSettings& the_step = (StepSettings&) the_constraint.GetSteps()(j);
			line = line + mystr("\t") + mystr(the_step.LoadFactor());
		}
		outfile.RWSoleStr(line);
	}
	for(int i=1; i<=NLoads(); i++)
	{
		FEMesh_Load& the_load = GetLoad(i);
		mystr line;
		if(the_load.Type() == TAreaLoad)
		{
			line = mystr("AreaSurfaceLoad\t");
		}
		else if(the_load.Type() == TFaceLoad)
		{
			line = mystr("FaceLoad\t");
		}
		else if(the_load.Type() == TBodyLoad) 
		{
			line = mystr("BodyLoad\t");
		}
		line = line + mystr(the_load.GetElement()) + mystr("\t");
		line = line + mystr(the_load.GetLoadVector().X()) + mystr("\t");
		line = line + mystr(the_load.GetLoadVector().Y()) + mystr("\t");
		line = line + mystr(the_load.GetLoadVector().Z()) + mystr("\t");
		line = line + mystr(the_load.GetSteps().Length());
		for(int j=1; j<=the_load.GetSteps().Length(); j++)
		{
			StepSettings& the_step = (StepSettings&) the_load.GetSteps()(j);
			line = line + mystr("\t") + mystr(the_step.LoadFactor());
		}
		outfile.RWSoleStr(line);
	}
	return 0;
}

//$ LA 2011-06-09: Export of elements, nodes and material properties to an Abaqus input text file
bool FEMesh::ExportToAbaqus(const char * filename)
{
	IVector material_mask;
	//Info: **xyz is ignored from Abaqus, treated as comment
	//      *xyz means that xyz is treated as keyword from Abaqus

	int i;
	ofstream abaqus_input;
	abaqus_input.open(filename);

	// catch general failures --> return false
	if(!abaqus_input.is_open())
		return false;
	if(GetMeshDim() != 2 && GetMeshDim() != 3)
		return false;

	// set precision for doubles
	abaqus_input.precision(8);

	abaqus_input << "*Heading" << endl;
	abaqus_input << " " << filename << endl;

	// Parts --> Bodys
	// only osed for naming, details are contained in instances (see assembly)
	abaqus_input << "**" << endl;
	abaqus_input << "**PARTS" << endl;
	abaqus_input << "**" << endl;
	int bodynr;
	
	CountBodies();
	IVector& bodynrs = mesh_register.Bodies();;
	for(int i=1; i<=bodynrs.Length(); i++)
	{
		abaqus_input << "*Part, name=body" << bodynrs(i) << endl;
		abaqus_input << "*End Part" << endl;
	}

	// Assembly
	abaqus_input << "**" << endl;
	abaqus_input << "**ASSEMBLY" << endl;
	abaqus_input << "**" << endl;
	abaqus_input << "*Assembly, name=HOTINTEXPORT" << endl;


	// Instances
	for(int i=1; i<=bodynrs.Length(); i++)  // if each body has only one material
	{
		abaqus_input << "*Instance, name=body" << bodynrs(i) << ", part=body" << bodynrs(i) << endl;
		abaqus_input << "0.0,0.0,0.0" << endl;                                 // position shift

		// write nodes of body;
		abaqus_input << "*Node, nset=n_body" << bodynrs(i) << endl;

		FEMesh_Set nodenrs(TSetNodes);
		GetNodesOfBody(nodenrs, bodynrs(i));
		for(int j=1; j<=nodenrs.NNodes(); j++)
		{
			Vector3D coord = GetPoint3D(nodenrs(j));
			abaqus_input << nodenrs(j) << ", " << coord.X() << ", " << coord.Y()  << ", " << coord.Z() << endl;
		}

		// write elements of body;
		FEMesh_Set elemnrs(TSetElements);
		GetElementsOfBody(elemnrs, bodynrs(i));         // this contains all elements of the body

// manual splitting into f.t.t.b. 4 element types
		IVector elems_hex_o1(0), elems_hex_o2(0), elems_prism_o1, elems_prism_o2;
		for(int j=1; j<=elemnrs.NElements(); j++)
		{
			if( (ElementIsHex(elemnrs(j))) && (GetElement(elemnrs(j)).Order() == 1) ) elems_hex_o1.Add(elemnrs(j));
			else if( (ElementIsHex(elemnrs(j))) && (GetElement(elemnrs(j)).Order() == 2) ) elems_hex_o2.Add(elemnrs(j));
			else if( (ElementIsPrism(elemnrs(j))) && (GetElement(elemnrs(j)).Order() == 1) ) elems_prism_o1.Add(elemnrs(j));
			else if( (ElementIsPrism(elemnrs(j))) && (GetElement(elemnrs(j)).Order() == 2) ) elems_prism_o2.Add(elemnrs(j));
			else GetMBS()->UO() << mystr("Elementtype not yet implemented in export function") + "\n";
		}
		for(int j=1; j<=4; j++)
		{
			IVector* p_elems;		
			mystr buffer("*Element, ");
			char letter = 'A';
			if(j==1) { p_elems = &elems_hex_o1; buffer += "type=C3D8R\n"; }
			else if(j==2) { p_elems = &elems_hex_o2; buffer += "type=C3D20\n"; }
			else if(j==3) { p_elems = &elems_prism_o1; buffer += "type=C3D6R\n"; }
			else if(j==4) { p_elems = &elems_prism_o2; buffer += "type=C3D15\n"; }
			else GetMBS()->UO() << mystr("Elementtype not yet implemented in export function") + "\n";

			if(p_elems->Length()) 
			{
				abaqus_input << buffer;
		//		abaqus_input << "elset=el_body" << bodynrs(i) << letter++ << endl;   
			}
			for(int k=1; k<=p_elems->Length(); k++)
			{
			  abaqus_input << p_elems->Get(k);
				FEMesh_Set nodenumbers(TSetElements); 
				GetNodesOfElement(nodenumbers, p_elems->Get(k));
				RemoveRedundantEntries(nodenumbers.Nodes());
				
// vgl Abaqus manual & QuadNodes routine
const int map_hexquad[] = {1, 2, 4, 3,  5, 6, 8, 7,  9, 14, 10, 13,  11, 16, 12, 15,  17, 18, 20, 19}; // hotint hex to abaqus hex
const int map_prismquad[] = {1, 2, 3,  5, 6, 7,  9, 14, 13,  11, 16, 15,  17, 18, 19};				// hotint degenerated hex to abaqus prism
				buffer.SetLength(0);
				if(j==1) // Hex O1: sequence 1,2,4,3, 5,6,8,7; 
				{
					for(int l=0; l<8; l++)
					{
						buffer += ", " + mystr(nodenumbers(map_hexquad[l]));
					}
				}
				else if(j==2) // Hex O2: sequence 1,2,4,3, 5,6,8,7, 9,14,10,13, 11,16,12,15, 17,18,20,19;
				{
					for(int l=0; l<20; l++)
					{
						buffer += ", " + mystr(nodenumbers(map_hexquad[l]));
					}
				}

				else if(j==3) // Prism O1: sequence 1,2,3, 5,6,7; (hex sequence & skip all l%4==2)
				{
					for(int l=0; l<6; l++)
					{
						buffer += ", " + mystr(nodenumbers(map_prismquad[l]));
					}
				}
				else if(j==4) // Prism O2: sequence 1,2,3, 5,6,7, 9,14,13, 11,16,15, 17,18,19; (hex sequence & skip all l%4==2)
				{
					for(int l=0; l<6; l++)
					{
						buffer += ", " + mystr(nodenumbers(map_prismquad[l]));
					}
				}
				abaqus_input << buffer << endl; 
			}
		}

// manual splitting into material sections
		material_mask.SetLen(NMaterials());
		material_mask.SetAll(0);

		for(int j=1; j<=NMaterials(); j++)
		{
			IVector elemmat(0);
			for(int k=1; k<=elemnrs.NElements(); k++)
			{
				if(GetElement(elemnrs(k)).MaterialNum()==j)
					elemmat.Add(elemnrs(k));
			}
			if(elemmat.Length())
			{
				material_mask(j) = 1;
				abaqus_input << "*Elset, elset=el_body" <<  bodynrs(i)<< "M" << j;
const int entries_per_line = 16;				
				for(int l=1; l<=elemmat.Length(); l++)
				{
					if(l%entries_per_line == 1)
					{
						abaqus_input << endl << elemmat(l);
					}
					else
					{
						abaqus_input << ", " << elemmat(l);
					}
				}
				abaqus_input << endl;
			}
		}

		// finish instance
		for(int j=1; j<=NMaterials(); j++)
		{
			if(material_mask(j))
				abaqus_input << "*Solid Section, elset=el_body" <<  bodynrs(i)<< "M" << j << ", material=" << GetMaterial(j).GetMaterialName() << "Nr" << j << endl;
		}	
		abaqus_input << "1.," << endl;
		abaqus_input << "*End Instance" << endl;

		//// write node sets of instance;   /e.g surface sets, ...

		//// write element sets of instance;
	}

	// write surfaces (dependent on elementsets);

	abaqus_input << "*End Assembly" << endl;

	// materials
	for( i=1; i<=NMaterials(); i++ )
	{
		if(material_mask(i))
		{
			abaqus_input << "*Material, name=" << GetMaterial(i).GetMaterialName() << "Nr" << i << endl;
			abaqus_input << "*Density" << endl;
			abaqus_input << GetMaterial(i).Density() << "," << endl;
			if(GetMaterial(i).GetType()&TMatFlexGeneral /*|| GetMaterial(i).GetType()&TMatFlexGeneral2D *//*|| GetMaterial(i).GetType()&TMatPlast*/) //$ JG: changed TMat2D, //$ PG 2013-8-1: TMatPlast has been deleted in class Material
			{
				abaqus_input << "*Elastic" << endl;
				abaqus_input << GetMaterial(i).YoungsModulus() << "," << GetMaterial(i).PoissonRatio() << endl;
			}
			else 
			{
				////GetMBS()->UO() << "Material type not yet implemented" << endl;
			}
		}
	}



	//** 
	//** INTERACTION PROPERTIES
	//** 
	//*Surface Interaction, name=Blech
	//1.,
	//*Friction, slip tolerance=0.005
	// 0.22,
	//*Surface Behavior, pressure-overclosure=EXPONENTIAL
	//0.012, 80.
	//
	//*Surface Interaction, name=contact1
	//1
	//*Friction, slip tolerance=1e-11
	//2.00000000E-001
	//*Surface Behavior,pressure-overclosure=linear
	//1.02400000E+005

	//** 
	//** INTERACTIONS
	//** 
	//*Contact Pair, interaction=Blech, type=SURFACE TO SURFACE, tracking=STATE
	//GH.Blech, Blech.Unten_Hi
	//*Contact Pair, interaction=contact1, type=SURFACE TO SURFACE
	//s_1_o, s_2_u
	//*Contact Pair, interaction=contact1, type=SURFACE TO SURFACE
	//s_2_o, s_3_u
	//*Contact Pair, interaction=contact1, type=SURFACE TO SURFACE
	//s_3_o, s_4_u




	//** ----------------------------------------------------------------
	//** 
	//** STEP: Testbsp 
	//** 
	//*Step, name=Testbsp, inc=100000
	//Testbsp
	//*Static
	//0.001, 1., 1e-8, 1.
	//** 
	//** BOUNDARY CONDITIONS
	//** 
	//**Fixierung Seitenränder 
	//*Boundary
	//n_a1_l, 1, 3
	//n_a1_r, 1, 1
	//** 
	//** LOADS
	//** 
	//**verteilter Druck
	//*Dsload
	//s_4_po, P, 1e-5
	//** 
	//** OUTPUT REQUESTS
	//** 
	//*Restart, write, overlay, frequency=1
	//** 
	//**FiELD OUTPUT
	//** 
	//*Output, field
	//*Contact Output
	//CSTRESS
	//*Node Output
	//COORD, U
	//*Element Output
	//S 
	//*Output, history, frequency=1
	//*End Step

	abaqus_input.close();
	return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ mesh opetations 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// elements

// find the highest ocurring material number in a subset of elements
int FEMesh::FindHighestMaterialNumber(FEMesh_Set subset_elements)
{
	if (subset_elements.NElements() == 0) subset_elements = FEMesh_Set(TSetElements,NaturalNumbers(this->NElements()));
	int highestmatnr = 0;
	for(int i=1; i<=subset_elements.NElements(); i++)
	{
		if (highestmatnr < GetElement(subset_elements(i)).MaterialNum()) highestmatnr = GetElement(subset_elements(i)).MaterialNum();
	}
	return highestmatnr;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// faces (and areas)

//finds elementnumber and sidenumber for a subset of faces
int FEMesh::MapFacesToElements(FEMesh_Set& subset_faces, int trybothcycles)
{
	bool useall = false;
	if(subset_faces.GetType() != TSetFaces)
    useall = true;
	if(subset_faces.GetFaces().Length() == 0)
    useall = true;
	if(useall)
	{
		IntVectorWithLengthN<1> temp;
		temp.SetArithmetic(1,1,NFaces());
		subset_faces.Faces() = temp;
		subset_faces.SetType(TSetFaces);
	}
	
	for(int i=1; i<= subset_faces.NFaces(); i++)
	{
		FindElementAndSideForFEFace(subset_faces(i), trybothcycles);
	}
	return 0;
}

//finds elementnumber and sidenumber for a subset of faces
int FEMesh::MapAreasToElements(TArray<int> & areastomap, int trybothcycles)
{
	FEMesh_Set allareas(TSetFaces,mystr("areas to map"));
	for(int i = 1; i<= areastomap.Length(); i++)
	{
		FEMesh_Set area(GetArea(areastomap(i)));
		allareas += area;
	}
	this->MapFacesToElements(allareas, trybothcycles);
	return 0;
}

//finds elementnumber and sidenumber for a FEMesh_Face with only nodes known
int FEMesh::FindElementAndSideForFEFace(int facenumber, int trybothcycles)
{
// MAJOR CHANGE: "face" is loaded from file: list of nodes
// NEW: face is loaded into FEMesh_face, nodes only
//      routine assigns element and side to the 
// rename the routine to sth like: "assign element and side to feface"

	FEMesh_Face& feface = GetFace(facenumber);
	int len = feface.RemoveRedundant();

	if(nodes_to_elements.Length() != NN() )
		ComputeNodesToElements();

	// possible elements - intersection of all ElementsOfNode: 
	FEMesh_Set parents(TSetElements, GetElementsOfNode(feface.GetNode(1)), mystr("possible elements for face"));
	for(int i=2; i<=len; i++)
	{
		TArray<int>& elemnodes = GetElementsOfNode(feface.GetNode(i));
		FEMesh_Set nodeelements(TSetElements, elemnodes);
		parents &= nodeelements;
	}

	// only few candidates left now... FEMesh_Face comparison
	int4 tempnodes;
	FEMesh_Face tempface;
	bool found = false;
	for(int i=1; i<= parents.N(); i++)
	{
		int elnr = parents(i);
		for(int j=1; j<= GetElement(elnr).NSides(); j++)
		{
			int templen = GetElement(elnr).GetSideNodeNum(j,tempnodes);
			tempface.SetFEMesh_Face(tempnodes, elnr, j);

			if( feface.IsCyclicEqual(tempface) )
			{
				feface.SetElement(elnr);
				feface.SetSide(j);
				found = true;
				break;
			}
			if(trybothcycles)
			{
				tempface.InvertRotationSense();
				if( feface.IsCyclicEqual(tempface) )
				{
					feface.SetElement(elnr);
					feface.SetSide(j);
					found = true;
					mbs->UO(UO_LVL_warn) << mystr("Warning: Face #") + mystr(facenumber) + mystr(" was inverted (trybothcycles == true) \n");

					break;
				}
			}
		}
	}
	if(!found)
	{
		mbs->UO(UO_LVL_warn) << mystr("Warning: Face #") + mystr(facenumber) + mystr(" could not be mapped to any existing element\n");
		return 0;
	}
	return 1;
}

// create faces from a list of nodes (can be restricted to subset of elements)
int FEMesh::FacesFromNodes(TArray<int>& facenrs, TArray<int>& nodenrs, TArray<int>& subset_validelements)
{
	facenrs.Flush();

	RemoveRedundantEntries(nodenrs,1); // !sort! array of node 

	for(int i=1; i<=nodenrs.Length(); i++)
	{
		int nodenr = nodenrs(i);
		TArray<int>& elements_of_node = GetElementsOfNode(nodenr);

		if(subset_validelements.Length() != 0)
		{
			elements_of_node = elements_of_node.Intersection(subset_validelements); // only elements from subset are allowed
		}

		//GetElementSides containing nodes
		for(int j=1; j<=elements_of_node.Length(); j++)
		{
			FEElement& elem = GetElement(elements_of_node(j));
			for(int k=1; k<= elem.NSides(); k++)
			{
				//int4 side = elem.GetSideNodeNum4(k);
				int4 side;
				//$EK 2013-03-04 in case of mixed surface elements (prism)
				if (k <= elem.NSides4())
					side = elem.GetSideNodeNum4(k);
        		if(side == int4(0,0,0,0) || k > elem.NSides4()) 
				{
					side = elem.GetSideNodeNum3(k);
					side(4) = side(1);
				}
// only interested, if loopnode is lowest node in the face
// face has already added in previous loop iteration otherwise
				int4 unsortedcopyofside = side;
				side.Sort();
				if(side(1) == nodenr)
				{
					int redlen = side.RemoveAllRedundantEntries();

					// check if all nodes of the face are in nodenrs array
					int not_found = 0;
					for(int l=2; l<=redlen; l++)
					{
						if(! nodenrs.Find(side(l)) )
							not_found = 1;
					}
					if(not_found) // do not add face if not all the nodes are in the nodenrs array
						continue;

					FEMesh_Face feface(unsortedcopyofside, elements_of_node(j), k);
					if(feface.NNodes() == 4 || feface.NNodes() == 3)
					{
						int newfacenr = Faces().Add(feface);
						facenrs.Add(newfacenr);
					}
				}
			}
		}
	}
	return facenrs.Length();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// surface
//creates all faces of the surfce in list FACES, saves according numbers in set SURFACE, optional filter NOT WORKING YET
void FEMesh::FindSurface(FEMesh_Set& subset_of_any_type)
{
// TO DO: Convert to subset of elements
// for the time being - use ALL elements
	IntVectorWithLengthN<1> allelems;
	allelems.SetArithmetic(1,1,NElements());
	subset_of_any_type.SetFEMesh_Set(TSetElements,allelems,mystr("subset elements"));
// SUBSET OF ELEMENTS - alias
	FEMesh_Set& subset = subset_of_any_type;
	int ne_global = NElements();
	int ne_subset = subset.N();

// refresh NodeToElement List
	if (NN() != nodes_to_elements.Length()) 
	{ 
		ComputeNodesToElements();
	}

// highest material number
	int nmat=0;
	for (int i=1; i <= ne_subset; i++)
	{ 
		if (nmat < GetElement(subset(i)).MaterialNum() ) nmat = GetElement(subset(i)).MaterialNum(); 
	}

// prepare arrays for all elements - assuming 6 faces per element 
// entries will contain material number of neighbor
// easier implementation if arrays are prepared for ALL elements, not just subset
	TArray<IVector*> element_to_faceflags(ne_global);
	for (int i=1; i <= ne_global; i++)
	{
		element_to_faceflags(i) = new IVector(6);
	}

// some variables used for comparison of degeneratred faces...
	TIntX& faceX = int4();
	TIntX* face;
	int4 face4;
	int3 face3;
	int2 face2;
	TArray<int> neighbours(20);

	//OLD: outer faces are marked as 0, other have element number of neighbour element
	//NEW: outer faces are marked as 0, other have material number of neighbour element
	for (int i = 1; i <= ne_global; i++)
	{
		// $EK 2013-02-18: the following is only necessary if subset is a real subset of elements
		if (subset.NElements() < ne_global)
		{
			if (! subset.Elements().Find(i) )
				continue; // skip all elements that are not in the subset 
		}
		GetNeighbourElements(i,neighbours);
		// $EK 2013-02-18: the following is only necessary if subset is a real subset of elements
		if (subset.NElements() < ne_global)
			neighbours &= subset.Elements();

		for (int j=1; j<=GetElement(i).NSides(); j++)
		{
			faceX.Set(0,0,0,0);
			int len = GetElement(i).GetSideNodeNum(j,faceX);
			int redlen = faceX.RemoveAllRedundantEntries();
			if(redlen == 4) 
			{ 
				face4.Set(faceX(1),faceX(2),faceX(3),faceX(4)); 
				face = &face4;
			} 
			else if(redlen == 3) 
			{ 
				face3.Set(faceX(1),faceX(2),faceX(3)); 
				face = &face3; 
			} 
			else if((redlen == 2) || (redlen == 1)) 
			{ 
				face2.Set(faceX(1),faceX(2)); 
				face = &face2; 
			} 
			else mbs->UO(UO_LVL_err) << "FindSurfaceElements: returned element side has more nodes then variable faceX can take!\n";

			int neighb = FaceInNeighbour(*face, neighbours);
			int found_in_subset = subset.Elements().Find(neighb);

			if (found_in_subset) 
				element_to_faceflags(i)->Add(GetElement(neighb).MaterialNum()); //NEW
			else 
			{ // no neighbour -> is surfaceelement: check dimension ?
				element_to_faceflags(i)->Add(0); // 0 for true surfaceelement
			}
		}
	}

	//outer faces are generated
	int quadcount = 0;
	int trigcount = 0;
	int linecount = 0;
	FEMesh_Face feface;

	for (int i = 1; i <= ne_global; i++)
	{
		if (! subset.Elements().Find(i) )
			continue; // skip all elements that are not in the subset 

		for (int j=1; j<=GetElement(i).NSides(); j++)
		{
			if (element_to_faceflags(i)->Get(j) == 0)  // entry 0 for surfaceelements
			{
				faceX.Set(0,0,0,0);
				int len = GetElement(i).GetSideNodeNum(j,faceX);
				int redlen = faceX.RemoveAllRedundantEntries();

				if(redlen == 4) 
				{ 
					quadcount++;
				} 
				else if(redlen == 3) 
				{ 
					trigcount++;
				}
				else if((redlen == 2) || (redlen == 1)) 
				{ 
					linecount++;
				}
				else mbs->UO(UO_LVL_err) << "FindSurfaceElements: returned element side has more nodes then variable faceX can take!\n";

				feface.SetNodes(faceX);
				feface.SetElement(i);
				feface.SetSide(j);
				int facenr =faces.Add(feface);
				Surface() << facenr;
			}
		}
	}
}

//get a list of elements which are neighbours of this element
void FEMesh::GetNeighbourElements(int elem, TArray<int>& neighbours)
{
	neighbours.SetLen(0);
	int nn=0;
	int el=0;

	for (int i=1; i<=GetElement(elem).NNodes(); i++)      
	{
		nn = GetElement(elem).GetNode(i);
		for (int j=1; j<=nodes_to_elements(nn)->Length(); j++)
		{
			el = nodes_to_elements(nn)->Get(j);
			if (!neighbours.Find(el) && el != elem) neighbours.Add(el);
		}
	}
}

//get a list of elements which are side-neighbours of this element (subset of node neighbors)
void FEMesh::GetNeighbourElements_Side(int elem, TArray<int> &neighbours)
{
	IVector nodeneighbours(0);
	neighbours.SetLen(0);

	GetNeighbourElements(elem,nodeneighbours);
	
	TIntX& faceX = int4();
	TIntX* face;
	int4 face4;
	int3 face3;
	int2 face2;

	for (int j=1; j<=GetElement(elem).NSides(); j++)
	{
		int len = GetElement(elem).GetSideNodeNum(j,faceX);
		int redlen = faceX.RemoveAllRedundantEntries();
		if(redlen == 4) { face4.Set(faceX(1),faceX(2),faceX(3),faceX(4)); face = &face4; } 
		else if(redlen == 3) { face3.Set(faceX(1),faceX(2),faceX(3)); face = &face3; } 
		else if((redlen == 2) || (redlen == 1)) { face2.Set(faceX(1),faceX(2)); face = &face2; } 
		else mbs->UO(UO_LVL_err) << "GetNeighbourElements_Side: returned element side has more nodes then variable faceX can take!\n";

		int neighb = FaceInNeighbour(*face, nodeneighbours);
		if (neighb) neighbours.Add(neighb); 
	}
}

//returns neighbour element number for the face or 0
int FEMesh::FaceInNeighbour(TIntX& face, const TArray<int>& neighbours)
{
	face.Invert();
	
	TIntX& ofaceX = int4();
	TIntX* oface;
	int4 oface4;
	int3 oface3;
	int2 oface2;
	
	for (int i=1; i <= neighbours.Length(); i++)
	{
		for (int j=1; j<=GetElement(neighbours(i)).NSides(); j++)
		{
			ofaceX.Set(0,0,0,0);
			int len = GetElement(neighbours(i)).GetSideNodeNum(j,ofaceX);
			int redlen = ofaceX.RemoveAllRedundantEntries();
			if(redlen == 4) 
			{ 
				oface4.Set(ofaceX(1),ofaceX(2),ofaceX(3),ofaceX(4)); 
				oface = &oface4; 
			} 
			else if(redlen == 3) 
			{ 
				oface3.Set(ofaceX(1),ofaceX(2),ofaceX(3)); 
				oface = &oface3; 
			} 
			else if((redlen == 2) || (redlen == 1)) 
			{ 
				oface2.Set(ofaceX(1),ofaceX(2)); 
				oface = &oface2; 
			} 
			else mbs->UO(UO_LVL_err) << "FindSurfaceElements: returned element side has more nodes then variable face can take!\n";
		
			if (face.IsCyclicEqual(*oface)) 
			{ face.Invert(); return neighbours(i); }
		}
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mirror

//adds new nodes, elements, faces, areas and new loads/constraints to internal arrays 	//calls FindSurface() unless flag is set to 0
void FEMesh::MirrorMeshAtPlane(Vector3D nplane, double cplane, double tol, int flag_recalc_surface)
{
	TArray<int> mirrornodes; 
	TArray<int> mirrorelements;
	TArray<int> mirrorfaces;
	TArray<int> mirrorareas;

	DuplicateNodesAtPlane(mirrornodes, nplane, cplane);
	DuplicateExistingElements(mirrorelements, mirrornodes);
	DuplicateExistingFaces(mirrorfaces, mirrornodes, mirrorelements);
	DuplicateExistingAreas(mirrorareas, mirrorfaces);

	DuplicateExistingLoads(mirrornodes, mirrorelements, mirrorfaces, mirrorareas, nplane);
	DuplicateExistingConstraints(mirrornodes, mirrorelements, mirrorfaces, mirrorareas, nplane);
	ComputeNodeBox();
	ComputeNodesToElements();
	if(flag_recalc_surface)
	{
		FindSurface();
	}
}

//adds new nodes, elements, faces, areas and new loads/constraints to internal arrays 	//calls FindSurface() unless flag is set to 0
void FEMesh::MirrorMeshAtPlane(Vector2D nplane, double cplane, double tol, int flag_recalc_surface)
{
	MirrorMeshAtPlane(nplane.MakeV3D(), cplane, tol, flag_recalc_surface);
}

// duplicates Nodes for mirrormesh
void FEMesh::DuplicateNodesAtPlane(IVector& mirrornode, Vector3D nplane, double cplane, double tol)
{
	double d;
	int nn = NN(); // length before mirrornodes are added 
	int d_pos=0;
	int d_neg=0;
	Vector3D n0(nplane); n0.Normalize();

	for(int i=1; i<=nn; i++) 
	{
		d = DistanceFromPlane(i,nplane,cplane);
		if(fabs(d)<tol) mirrornode.Add(i); // at plane, dont add 
		else
		{
			if (d>0) d_pos++;
			else d_neg++;

			Vector3D mirrorpos = GetPoint3D(i) - 2.0 * d * n0;
			int newnodenumber = AddNodeCheck(mirrorpos,GetNodeDomain(i));
			mirrornode.Add(newnodenumber); 
		}
	}
	if( d_pos && d_neg)
		mbs->UO().InstantMessageText("Warning: Mesh has nodes on both sides of the mirrorplane,\nintersecting elements are not checked by this routine.");
}

// duplicates Elements for mirrormesh
void FEMesh::DuplicateExistingElements(IVector& mirrorelement, const IVector& mirrornodes)
{
	mirrorelement.Flush();
  int ne = NElements();

	for(int i=1; i<=ne; i++)
	{
		FEElement* melem = elements(i)->GetCopy();
		for(int j=1; j<=melem->NNodes(); j++) 
		{
			melem->SetNode(j,mirrornodes(melem->GetNode(j))); // nodenumbers of the duplicate
		}
		melem->Invert();
		int newelementnumber = AddElement(*melem);
		mirrorelement.Add(newelementnumber);
		delete melem;
	}
}

	// duplicates Faces for mirrormesh - assume Nodes&Elements already duplicated - only faces that have a mapping to a element/side
void FEMesh::DuplicateExistingFaces(IVector& mirrorfaces, const IVector& mirrornodes, const IVector& mirrorelements)
{
	mirrorfaces.Flush();
	int nf = NFaces();

	for(int i=1; i<=nf; i++)
	{
		FEMesh_Face& originalface = GetFace(i);
		FEMesh_Face copyface;
		for(int j=1; j<=originalface.NNodes(); j++)
		{
			copyface.SetNode(j, mirrornodes(originalface.GetNode(j)));    
		}

		if(originalface.GetElement() && originalface.GetSide())
		{
			copyface.SetElement(mirrorelements(originalface.GetElement()));
			copyface.SetSide(originalface.GetSide());

			int newfacenr = Faces().Add(copyface);
			mirrorfaces.Add(newfacenr);
		}
	}
}

// duplicate Areas for mirrormesh
void FEMesh::DuplicateExistingAreas(IVector& mirrorareas, const IVector&mirrorfaces)
{
	mirrorareas.Flush();
	int na = NAreas();

	for(int i=1; i<=na; i++)
	{
		TArray<int>& oldarea_facenumbers = GetArea(i).Faces();
		int nfaces = oldarea_facenumbers.Length();
		FEMesh_Set newarea;
		for(int j=1; j<=nfaces; j++)
		{
			newarea.Add(mirrorfaces(oldarea_facenumbers(j)));
		}
		int newareanumber = AddArea(newarea);
		mirrorareas.Add(newareanumber);
	}
}

// duplicates existing Loads - assume nodes, elements, faces and areas duplicated
void FEMesh::DuplicateExistingLoads(const IVector& mirrornodes, const IVector& mirrorelements, const IVector& mirrorfaces, const IVector& mirrorareas, const Vector3D& nplane)
{
	int nl = NLoads();
	
	for(int i=1; i <= nl; i++)
	{
		FEMesh_Load& load = GetLoad(i);
		int elem = load.GetElement();
		Vector3D lvm = MirrorAtPlane3D(load.GetLoadVector(), nplane, 0.0);

		if(load.GetType() == TAreaLoad)
		{
			FEMesh_AreaLoad newload;
			newload.SetElement( mirrorareas(elem) );  
			newload.SetLoadVector( lvm ); // mirror the direction of force at the mirrorplane as well
			AddLoad(newload);
		}
		else if(load.GetType() == TFaceLoad)
		{
			FEMesh_FaceLoad newload;
			newload.SetElement( mirrorfaces(elem) );
			newload.SetLoadVector( lvm ); // mirror the direction of force at the mirrorplane as well
			AddLoad(newload);
		}
		else if(load.GetType() == TBodyLoad)
		{
			FEMesh_BodyLoad newload( (FEMesh_BodyLoad&) load);
			newload.SetGravity( ((FEMesh_BodyLoad&) load).GetGravity() );
			newload.SetElement( mirrorelements(elem) );
			if(newload.GetGravity())
			{
				newload.SetLoadVector( ((FEMesh_BodyLoad&) load).GetLoadVector() ); // body load direction is not mirrored ?
			}
			else
			{
					newload.SetLoadVector( ((FEMesh_BodyLoad&) load).GetLoadVector() ); // gravity is NEVER mirrored 
			}
			AddLoad(newload);
		}
	}
}

// duplicates existing Constraints  - assume nodes, elements, faces and areas duplicated
void FEMesh::DuplicateExistingConstraints(const IVector& mirrornodes, const IVector& mirrorelements, const IVector& mirrorfaces, const IVector& mirrorareas, const Vector3D& nplane)
{
	int nc = NConstraints();
	
	for(int i=1; i <= nc; i++)
	{
		FEMesh_Constraint& constraint = GetConstraint(i);

		if(constraint.GetType() == TAreaConstraint)
		{
			FEMesh_AreaConstraint newconstraint( (FEMesh_AreaConstraint&) constraint);
			newconstraint.SetElement( mirrorareas(newconstraint.GetElement()) );
			AddConstraint(newconstraint);
		}
		else if(constraint.GetType() == TFaceConstraint)
		{
			FEMesh_FaceConstraint newconstraint( (FEMesh_FaceConstraint&) constraint);
			newconstraint.SetElement( mirrorfaces(newconstraint.GetElement()) );
			AddConstraint(newconstraint);
		}
		else if(constraint.GetType() == TNodalConstraint)
		{
			FEMesh_NodeConstraint newconstraint( (FEMesh_NodeConstraint&) constraint);
			newconstraint.SetElement( mirrornodes(newconstraint.GetElement()) );
			AddConstraint(newconstraint);
		}
		else if(constraint.GetType() == TAreaContact)
		{
			FEMesh_AreaContact newconstraint( (FEMesh_AreaContact&) constraint);
			newconstraint.SetElement( mirrorareas(newconstraint.GetElement()) );
			newconstraint.SetTarget( mirrorareas(newconstraint.GetTarget()) );
			AddConstraint(newconstraint);
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// refine

// ATTENTION: works with hexahedrals only so far 
// refine the mesh with a factor 2 in all directions
void FEMesh::RefineMesh2(TArray<int>& subset)
{
	IVector allrv(0);
	if (subset.Length() == 0) 
	{
		for(int i=1; i<=NElements(); i++) 
			subset.Add(i);
	}

	int ne = NElements();
  for(int i=1; i<= ne; i++)
	{
		if(!subset.Find(i))
			continue;
		FEElement& elem(GetElement(i));
		if(elem.Order() != 1) // this subroutine is only for linear elements
		{
			mbs->UO(UO_LVL_ext) << mystr("RefineMesh only for linear elements: Element #") + mystr(i) + mystr(" is not linear, no changes for this element \n");
			break;
		}
	
		// identify the face numbers of faces for this element // areas to modify
		IVector facenrs;	
		IVector areanrs;
		for(int j=1; j <= NFaces(); j++)
		{
			FEMesh_Face& feface = GetFace(j);
			if(feface.GetElement() == i) // found face on this element
			{
				facenrs.Add(j);
				for(int k=1; k <= areas.Length(); k++)
				{
					for(int l=1; l <= areas(k).NFaces(); l++)
					{
						if(areas(k).Faces().Get(l) == j) // found face in area
						{
							areanrs.Add(k);
						}
					}
				}
			}
		}
		::RemoveRedundantEntries(facenrs,1);
		::RemoveRedundantEntries(areanrs,1);

		// compute the new nodes
		IVector refnodes;
		for(int j=1; j <= elem.NNodes(); j++)
		{
			refnodes.Add(elem.GetNode(j));
		}
		RefinedNodes(refnodes);

		// make new elements
		IVector newelements; // remember the indices of the new elements
		for(int iz = 1; iz <= 2; iz++)
		{
			for(int iy = 1; iy <= 2; iy++)
			{
				for(int ix = 1; ix <= 2; ix++)
				{
					IVector points;
					int offset = (ix-1)*1 + (iy-1)*3 + (iz-1)*9;    // offset for first node
					int index = (ix-1)*1 + (iy-1)*2 + (iz-1)*4 +1;  // index 1..8 for new hex
					
					points.Add(refnodes(offset + 1));
					points.Add(refnodes(offset + 2));
					points.Add(refnodes(offset + 4));
					points.Add(refnodes(offset + 5));
					points.Add(refnodes(offset + 10));
					points.Add(refnodes(offset + 11));
					points.Add(refnodes(offset + 13));
					points.Add(refnodes(offset + 14));

					FEElement* newelem = MakeHex(points,elem.MaterialNum());
					newelem->Color() = elem.Color();
					if(offset == 0)
					{
						elements(i) = newelem; 
						newelements(index) = i;
					}
					else
					{
						newelements(index) = AddElement(newelem);
					}
				}
			}
			allrv.Append(newelements);
		}

		// make new faces and entries in surfaceparents / area
		for(int j=1; j <= facenrs.Length(); j++)
		{
			int oldfacenr = facenrs(j);
			int oldelem = GetFace(oldfacenr).GetElement();
			int oldside = GetFace(oldfacenr).GetSide();

			IVector newfacenrs;
			IVector pickedelems;

			int4 tempnodes = elements(i)->GetSide4(oldside);
			pickedelems.Set4(tempnodes(1),tempnodes(2),tempnodes(3),tempnodes(4)); // select the daughter elements containing the new faces

			for(int k=1; k <= 4; k++) // create 4 new faces for each old face
			{
				int elnr = newelements(pickedelems(k));
				tempnodes = GetElement(elnr).GetSideNodeNum4(oldside); // nodenumbers of the new face

				FEMesh_Face newface(tempnodes, elnr, oldside);
				int newfacenr = Faces().Add(newface);
				newfacenrs.Add(newfacenr);
			}

			for(int k=1; k <= areanrs.Length(); k++)
			{
				IVector* looparea = & areas(areanrs(k)).Faces();
				for(int l=1; l <= looparea->Length(); l++)
				{
					if( looparea->Get(l) == oldfacenr)
					{
						looparea->Erase(l);
						looparea->Add(newfacenrs(1));
						looparea->Add(newfacenrs(2));
						looparea->Add(newfacenrs(3));
						looparea->Add(newfacenrs(4));
					}
				}
			}
		}
	}
	subset = allrv; // return value: all refined elements
}

// calculate refined node positions(27), create nodes, return these nodes numbers 
// !! ASSUME LINEAR HEXAHEDRALS , FACTOR 2 refinement!! 
void FEMesh::RefinedNodes(IVector& points)
{
	int domain = this->GetNodeDomain(points(1));
	Vector3D p1(GetPoint3D(points(1)));
	Vector3D p2(GetPoint3D(points(2)));
	Vector3D p3(GetPoint3D(points(3)));
	Vector3D p4(GetPoint3D(points(4)));
	Vector3D p5(GetPoint3D(points(5)));
	Vector3D p6(GetPoint3D(points(6)));
	Vector3D p7(GetPoint3D(points(7)));
	Vector3D p8(GetPoint3D(points(8)));
	points.Flush();

	points.Add(AddNodeCheck(p1,domain));
	points.Add(AddNodeCheck((p1+p2)*0.5,domain));
	points.Add(AddNodeCheck(p2,domain));
	
	points.Add(AddNodeCheck((p1+p3)*0.5,domain));
	points.Add(AddNodeCheck((p1+p2+p3+p4)*0.25,domain));
	points.Add(AddNodeCheck((p2+p4)*0.5,domain));

	points.Add(AddNodeCheck(p3,domain));
	points.Add(AddNodeCheck((p3+p4)*0.5,domain));
	points.Add(AddNodeCheck(p4,domain));

	points.Add(AddNodeCheck((p1+p5)*0.5,domain));
	points.Add(AddNodeCheck((p1+p2+p5+p6)*0.25,domain));
	points.Add(AddNodeCheck((p2+p6)*0.5,domain));

	points.Add(AddNodeCheck((p1+p3+p5+p7)*0.25,domain));
	points.Add(AddNodeCheck((p1+p2+p3+p4+p5+p6+p7+p8)*0.125,domain));
	points.Add(AddNodeCheck((p2+p4+p6+p8)*0.25,domain));

	points.Add(AddNodeCheck((p3+p7)*0.5,domain));
	points.Add(AddNodeCheck((p3+p4+p7+p8)*0.25,domain));
	points.Add(AddNodeCheck((p4+p8)*0.5,domain));

	points.Add(AddNodeCheck(p5,domain));
	points.Add(AddNodeCheck((p5+p6)*0.5,domain));
	points.Add(AddNodeCheck(p6,domain));
	
	points.Add(AddNodeCheck((p5+p7)*0.5,domain));
	points.Add(AddNodeCheck((p5+p6+p7+p8)*0.25,domain));
	points.Add(AddNodeCheck((p6+p8)*0.5,domain));

	points.Add(AddNodeCheck(p7,domain));
	points.Add(AddNodeCheck((p7+p8)*0.5,domain));
	points.Add(AddNodeCheck(p8,domain));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// linear <--> quadratic

// make all linear elements quadratic - consistency !NOT! conserved unless all elements are linear in the beginning 
void FEMesh::LinearToQuadratic(IVector& subset)
{
	if (subset.Length() == 0) 
	{
		for(int i=1; i<=NElements(); i++) //*JG: changed from old code "i<NElements()" to <=; otherwise not all elements are converted!
			subset.Add(i);
	}
  IVector newnodes(20);
	for(int i=1; i<=subset.Length(); i++)
	{
		FEElement* linelem = elements(subset(i));

		if(linelem->Order() != 1) // this subroutine is only for linear elements
		{
			mbs->UO(UO_LVL_ext) << mystr("LinearToQuadratic: Element #") + mystr(subset(i)) + mystr(" is not linear, no changes for this element \n");
			break;
		}
		newnodes.Flush();
		QuadNodes(newnodes,linelem);
		Vector3D color = linelem->Color();

		if(linelem->Type() == TFELine) elements(subset(i)) = MakeLinequad(newnodes,linelem->MaterialNum());
		if(linelem->Type() == TFETrig) elements(subset(i)) = MakeTrigquad(newnodes,linelem->MaterialNum());
		if(linelem->Type() == TFEQuad) elements(subset(i)) = MakeQuadquad(newnodes,linelem->MaterialNum());
		if(linelem->Type() == TFETet) elements(subset(i)) = MakeTetquad(newnodes,linelem->MaterialNum());
		if(linelem->Type() == TFEHex) elements(subset(i)) = MakeHexquad(newnodes,linelem->MaterialNum());
		// $EK 2013-03-05 added functions for PrismQuad and PyramidQuad
		if(linelem->Type() == TFEPrism) elements(subset(i)) = MakePrismquad(newnodes,linelem->MaterialNum());
		if(linelem->Type() == TFEPyramid) 
		{
			mbs->UO(UO_LVL_warn) << "WARNING in FEMesh::LinearToQuadratic - pyramids are not implemented with quadratic shape functions!\n";
			elements(subset(i)) = MakePyramidquad(newnodes,linelem->MaterialNum());
		}
		delete linelem;	 // deletes old element
		elements(subset(i))->Color() = color;
	}
}

// make all quadratic elements linear- consistency !NOT! conserved unless all elements are quadratic in the beginning 
void FEMesh::QuadraticToLinear(IVector& subset, int removeunused)
{
	if (subset.Length() == 0) 
	{
		for(int i=1; i<=NElements(); i++)
			subset.Add(i);
	}
  IVector newnodes(20);
	for(int i=1; i<=subset.Length(); i++)
	{
		FEElement* quadelem = elements(subset(i));
		if(quadelem->Order() != 2) // this subroutine is only for quadratic elements
		{
			mbs->UO(UO_LVL_ext) << mystr("LinearToQuadratic: Element #") + mystr(subset(i)) + mystr(" is not quadratic, no changes for this element \n");
			break;
		}
		newnodes.Flush();
		LinearNodes(newnodes,quadelem);
		Vector3D color = quadelem->Color();
		
		if(quadelem->Type() == 0) elements(subset(i)) = MakeLine(newnodes,quadelem->MaterialNum());
		if(quadelem->Type() == 1) elements(subset(i)) = MakeTrig(newnodes,quadelem->MaterialNum());
		if(quadelem->Type() == 2) elements(subset(i)) = MakeQuad(newnodes,quadelem->MaterialNum());
		if(quadelem->Type() == 3) elements(subset(i)) = MakeTet(newnodes,quadelem->MaterialNum());
		if(quadelem->Type() == 4) elements(subset(i)) = MakeHex(newnodes,quadelem->MaterialNum());
		delete quadelem;	 // deletes old element
		elements(subset(i))->Color() = color;
	}
	if (removeunused)
	{
		this->RemoveUnusedNodes();
	}
}

// returns nodelist for the quadratic element, adds intermediate nodes 
void FEMesh::QuadNodes(IVector& newpoints, FEElement* linelem)
{
	IVector oldpoints;
	for (int i=1; i<=linelem->NNodes(); i++)
		oldpoints.Add(linelem->GetNode(i));

	newpoints = oldpoints;

	if((linelem->Type() == TFELine) && (newpoints.Length() == 2)) 
	{
		newpoints.Add(AddIntermediateNode(linelem->GetNode(1),linelem->GetNode(2)));
	}

	if((linelem->Type() == TFETrig) && (newpoints.Length() == 3)) 
	{
		const int inter_trig1[] = { 1, 2, 3};
		const int inter_trig2[] = { 2, 3, 1};
		for(int i=1; i<=3; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_trig1[i-1]),linelem->GetNode(inter_trig2[i-1])));
	}
	if((linelem->Type() == TFEQuad) && (newpoints.Length() == 4)) //!AD: 2012-03-07 changed to 9-node element
	{
		const int inter_quad1[] = { 1, 2, 3, 4, 1};
		const int inter_quad2[] = { 2, 3, 4, 1, 3};
		for(int i=1; i<=5; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_quad1[i-1]),linelem->GetNode(inter_quad2[i-1])));
	}
	if((linelem->Type() == TFETet) && (newpoints.Length() == 4)) 
	{
		//const int inter_tet1[] = { 1, 2, 3, 3, 1, 2}; // wrong
		//const int inter_tet2[] = { 2, 3, 1, 4, 4, 4}; // wrong
		const int inter_tet1[] = { 1, 2, 3, 2, 3, 1};
		const int inter_tet2[] = { 2, 3, 1, 4, 4, 4};
		for(int i=1; i<=6; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_tet1[i-1]),linelem->GetNode(inter_tet2[i-1])));
	}
	if((linelem->Type() == TFEHex) && (newpoints.Length() == 8)) 
	{
		const int inter_hex1[] = { 1, 3, 5, 7, 1, 2, 5, 6, 1, 2, 3, 4};
		const int inter_hex2[] = { 2, 4, 6, 8, 3, 4, 7, 8, 5, 6, 7, 8};
		for(int i=1; i<=12; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_hex1[i-1]),linelem->GetNode(inter_hex2[i-1])));
	}
	// $EK 2013-03-06 added for quadratic prisms
	if((linelem->Type() == TFEPrism) && (newpoints.Length() == 6)) 
	{
		const int inter_prism1[] = { 1, 2, 3, 1, 2, 3, 4, 5, 6};
		const int inter_prism2[] = { 2, 3, 1, 4, 5, 6, 5, 6, 4};
		for(int i=1; i<=9; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_prism1[i-1]),linelem->GetNode(inter_prism2[i-1])));
	}
	// $EK 2013-03-06 added for quadratic pyramids
	if((linelem->Type() == TFEPyramid) && (newpoints.Length() == 5)) 
	{
		const int inter_pyramid1[] = { 1, 2, 4, 3, 1, 2, 3, 4};
		const int inter_pyramid2[] = { 2, 4, 3, 1, 5, 5, 5, 5};
		for(int i=1; i<=8; i++) 
			newpoints.Add(AddIntermediateNode(linelem->GetNode(inter_pyramid1[i-1]),linelem->GetNode(inter_pyramid2[i-1])));
	}
	return; 
}

// returns nodelist for the linear element 
void FEMesh::LinearNodes(IVector& newpoints, FEElement* quadelem)
{
	newpoints.Flush();
	if(quadelem->Type() == 0)
	{
		for (int i=1; i<=2; i++) newpoints.Add(quadelem->GetNode(i));
	}
	if(quadelem->Type() == 1)
	{
		for (int i=1; i<=3; i++) newpoints.Add(quadelem->GetNode(i));
	}
	if(quadelem->Type() == 2)
	{
		for (int i=1; i<=4; i++) newpoints.Add(quadelem->GetNode(i));
	}
		if(quadelem->Type() == 3)
	{
		for (int i=1; i<=4; i++) newpoints.Add(quadelem->GetNode(i));
	}
	if(quadelem->Type() == 4)
	{
		for (int i=1; i<=8; i++) newpoints.Add(quadelem->GetNode(i));
	}
	return;
}

// creates intermediate node between two existing (global number) nodes, adds to nodetree
int FEMesh::AddIntermediateNode(int i, int j)
{
	Vector3D newpos = 0.5* (GetPoint3D(i)+GetPoint3D(j));
	int n = AddNodeCheck( newpos, GetNodeDomain(i) );
	return n;
}

// creates a node in the center of all nodes registered in array points, adds the index of this node to points
int FEMesh::AddCenterNode(TArray<int>& points)
{
	Vector3D newpos(0.,0.,0.);
	for(int i=1; i<=points.Length(); i++)
	{
		newpos += GetPoint3D(points(i));
	}
	newpos *= (1./points.Length());
	int n = AddNodeCheck( newpos, GetNodeDomain(points(1)) );
	points.Add(n);
	return n;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// automatic joints

//join all nodes of slave to faces of master, returns number of AreaContact
int FEMesh::JoinAtPlane(FEMesh_Set& master_elements, FEMesh_Set& slave_elements, Vector3D point, Vector3D normal)
{
	FEMesh_Set master_nn(TSetNodes), master_nn_plane(TSetNodes), master_fn(TSetFaces);
	FEMesh_Set slave_nn(TSetNodes),  slave_nn_plane(TSetNodes),  slave_fn(TSetFaces);

	double d = (normal * point) / normal.Norm();

	GetNodesOfElements(master_nn, master_elements);
	GetNodesOnPlane(master_nn_plane, normal, d, master_nn);
	FacesFromNodes(master_fn.Faces(), master_nn_plane.Nodes(), master_elements.Elements());
	int area_master = AddArea(master_fn);

	GetNodesOfElements(slave_elements, slave_nn);
	GetNodesOnPlane(slave_nn_plane, normal, d, slave_nn);
	FacesFromNodes(slave_fn.Faces(), slave_nn_plane.Nodes(), slave_elements.Elements());
	int area_slave = AddArea(slave_fn);

	return AddConstraint(FEMesh_AreaContact(area_slave, area_master, 0, 0, 0));
}

//join all nodes of slave to faces of master, returns number of AreaContact
int FEMesh::JoinBodiesAtPlane(int master_body, int slave_body, Vector3D point, Vector3D normal)
{
	FEMesh_Set master_en(TSetElements);
	FEMesh_Set slave_en(TSetElements);

	GetElementsOfBody(master_en, master_body );
	GetElementsOfBody(slave_en,  slave_body );

	return JoinAtPlane(master_en, slave_en, point, normal);
}

//join all nodes of slave to faces of master, returns number of AreaContact
int FEMesh::JoinMaterialAtPlane(int master_material, int slave_material, Vector3D point, Vector3D normal)
{
	FEMesh_Set master_en(TSetElements);
	FEMesh_Set slave_en(TSetElements);

	GetElementsOfMaterial(master_en, master_material);
	GetElementsOfMaterial(slave_en, slave_material);

	return JoinAtPlane(master_en, slave_en, point, normal);
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ operations on sets and single entities
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// elements

// changes the domain of a set of elements, adds additional nodes, also updates faces if they are associated with elements
void FEMesh::SetElementDomain(FEMesh_Set& in_set_elements, int in_int_domain_number)
{
	// Step 1: identify nodes in all the chosen elements and create nodes for new domain if required
	FEMesh_Set nodes_in_elements(TSetNodes);
  GetNodesOfElements(nodes_in_elements, in_set_elements);

	FEMesh_Set nodes_in_new_domain(TSetNodes);
	for(int i=1; i<=nodes_in_elements.NNodes(); i++)
	{
		FEMesh_Node& the_oldnode = GetNode(nodes_in_elements(i));
		int newnodenr = AddNodeCheck(the_oldnode.Coords3D(), in_int_domain_number);
		nodes_in_new_domain.Add(newnodenr);
	}

	// Step 2: identify faces associated with chosen elements and exchange old node numbers for the new ones
	for(int i=1; i<=NFaces(); i++)
	{
		FEMesh_Face& the_face = GetFace(i);
		if( in_set_elements.Elements().Find(the_face.Element()) ) // determine if the node is on an element that changed domain
		{
			for(int j=1; j<=the_face.NNodes(); j++)
			{
				int found = nodes_in_elements.Nodes().Find(the_face.Node(j));
				if( found )
				{
					the_face.Node(j) = nodes_in_new_domain(found);
				}
			}
		}
	}

	// Step 3: for the chosen elements exchange old node numbers for the new ones and domainnumber in the course as well
	for(int i=1; i<=in_set_elements.NElements(); i++)
	{
		FEElement& the_element = GetElement(in_set_elements(i));
		for(int j=1; j<=the_element.NNodes(); j++)
		{

			int found = nodes_in_elements.Nodes().Find(the_element.GetNode(j));
			if( found )
			{
				int newnnr = nodes_in_new_domain(found);
				the_element.SetNode(j,newnnr);
			}
			else
			{
				mbs->UO() << mystr("ERROR: routine FEMesh::SetElementDomain failed to find node");
				// this should not happen
			}
		}
		the_element.Domain() = in_int_domain_number;
	}

	// Step 4: cleanup
	RemoveUnusedNodes();
}

// changes the elementcolor for a set of elements
void FEMesh::SetElementColor(FEMesh_Set& in_set_elements, Vector3D in_V3D_new_color)
{
	for(int i=1; i<=in_set_elements.NElements(); i++)
	{
		GetElement(i).Color() = in_V3D_new_color;
	}
}

// splits 1 prism into 3 hexehedrals - NOT IMPLEMENTED YET
void FEMesh::PrismTo3Hexes(int elemnr)
{
// assume element is prism
	FEElement& original = GetElement(elemnr);

// identify the triangle-faces
	int triangleside = 0;
	for(int i=1; i <= 6; i++)
	{
		int4 face = original.GetSide4(i);
		int len = face.RemoveAllRedundantEntries();
		if(len == 3)
		{
			triangleside = i; // triangle identified
			break; // break loop
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ filter functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// general

// determine if the Node is used in 3D elements only
int FEMesh::Is3DNode(int nodenumber)
{
	TArray<int>& elementlist = GetElementsOfNode(nodenumber);
	for(int i=1; i<= elementlist.Length(); i++)
	{
		FEElement& elem = GetElement(elementlist(i));
		if(elem.Dim() != 3)
			return 0;
	}
	return 3;
}

// determine if the Node is used in 2D elements only
int FEMesh::Is2DNode(int nodenumber)
{
	TArray<int>& elementlist = GetElementsOfNode(nodenumber);
	for(int i=1; i<= elementlist.Length(); i++)
	{
		FEElement& elem = GetElement(elementlist(i));
		if(elem.Dim() != 2)
			return 0;
	}
	return 2;
}

// determines if the Finite Element element is really a hexahedral
int FEMesh::ElementIsHex(int elemnr)
{
	FEElement& elem = GetElement(elemnr);
	FEMesh_Set nodenumbers(TSetNodes);
	GetNodesOfElement(nodenumbers, elemnr);

	if( (elem.Order() == 1) && (nodenumbers.NNodes() == 8) ) // linear hex: 8 nodes
		return 1;
	if( (elem.Order() == 2) && (nodenumbers.NNodes() == 20) ) // quadratic hex: 20 nodes
		return 1;
	return 0;
}

// determines if the Finite Element element is a prism
int FEMesh::ElementIsPrism(int elemnr)
{
	FEElement& elem = GetElement(elemnr);
	FEMesh_Set nodenumbers(TSetNodes);
	GetNodesOfElement(nodenumbers, elemnr);

	if( (elem.Order() == 1) && (nodenumbers.NNodes() == 6) ) // linear prism: 6 nodes
		return 1;
	if( (elem.Order() == 2) && (nodenumbers.NNodes() == 15) ) // quadratic prism: 15 nodes
		return 1;
	return 0;
}

// determines if the Finite Element element is a pyramid
int FEMesh::ElementIsPyram(int elemnr)
{
	FEElement& elem = GetElement(elemnr);
	FEMesh_Set nodenumbers(TSetNodes);
	GetNodesOfElement(nodenumbers, elemnr);

	if( (elem.Order() == 1) && (nodenumbers.NNodes() == 5) ) // linear pyram: 5 nodes
		return 1;
	if( (elem.Order() == 2) && (nodenumbers.NNodes() == 13) ) // quadratic pyram: 13 nodes
		return 1;
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// convert sets  functions - use FEMesh_Set only
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// computes a set of nodes for a given set of elements, subset of valid nodes can be included
void FEMesh::ComputeNodesOfElements(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_elements, FEMesh_Set& subset_nodes)
{
	rv_set_nodes.Nodes().Flush();
	TArray<int> tempnodes;

	for(int i=1; i<=in_set_elements.NElements(); i++)
	{
		FEElement& elem = GetElement(in_set_elements.Elements()(i));
		tempnodes.Flush();
		for(int j=1; j<= elem.NNodes(); j++)
		{
			tempnodes.Add(elem.GetNode(j));
		}
		rv_set_nodes.Nodes() += tempnodes;
	}

	if(subset_nodes.NNodes())
	{
		rv_set_nodes.Nodes() &= subset_nodes.Nodes();
	}
}

// computes a set of ALL elements for a given set of nodes, subset of valid elements can be included 
void FEMesh::ComputeElementsOfNodes(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();

	for(int i=1; i<=in_set_nodes.NNodes(); i++)
	{
		rv_set_elements.Elements() += GetElementsOfNode( in_set_nodes.Nodes()(i) );
	}

	if(subset_elements.NElements())
	{
		rv_set_elements.Elements() &= subset_elements.Elements();
	}
}

// computes a set of elements that can be built with a given set of nodes, subset of valid elements can be included 
void FEMesh::ComputeElementsFromNodes(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();

	// can loop over subset here
	if(subset_elements.NElements() == 0)
	{
		subset_elements.Elements() = NaturalNumbers(NElements());
	}

	for(int i=1; i<=subset_elements.NFaces(); i++)
	{
		int elemnr = subset_elements.Elements()(i);
		FEElement& elem = GetElement(elemnr);
		int nodesfound = 0;
		for(int j=1; j<=elem.NNodes(); j++)
		{
			if( in_set_nodes.Nodes().Find(elem.GetNode(j)) ) 
			{
				nodesfound++;
			}
		}
		if( nodesfound == elem.NNodes() )  // ALL nodes found 
		{
			rv_set_elements.Elements().Add(elemnr);
		}
	}
}

// computes a set of nodes for a given set of faces, subset of valid nodes can be included
void FEMesh::ComputeNodesOfFaces(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_faces, FEMesh_Set& subset_nodes)
{
	rv_set_nodes.Nodes().Flush();
	TArray<int> tempnodes;

	for(int i=1; i<=in_set_faces.NFaces(); i++)
	{
		FEMesh_Face& face = GetFace(in_set_faces.Faces()(i));
		tempnodes.Flush();
		for(int j=1; j<= face.NNodes(); j++)
		{
			tempnodes.Add(face.GetNode(j));
		}
		rv_set_nodes.Nodes() += tempnodes;
	}

	if(subset_nodes.NNodes())
	{
		rv_set_nodes.Nodes() &= subset_nodes.Nodes();
	}
}

// computes a set of ALL faces for a given set of nodes, subset of valid faces can be included
void FEMesh::ComputeFacesOfNodes(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_faces)
{
	rv_set_faces.Faces().Flush();

	// can loop over subset here
	if(subset_faces.NFaces() == 0)
	{
		subset_faces.Faces() = NaturalNumbers(NFaces());
	}
	
	for(int i=1; i<=subset_faces.NFaces(); i++)
	{
		int facenr = subset_faces.Faces()(i);
		FEMesh_Face& face = GetFace(facenr);
		for(int j=1; j<=face.NNodes(); j++)
		{
			if( in_set_nodes.Nodes().Find(face.GetNode(j)) ) 
			{
				rv_set_faces.Faces().Add(facenr); // ANY node found
				continue;
			}
		}
	}
}

// computes a set of faces that can be built a given set of nodes, subset of valid faces can be included
void FEMesh::ComputeFacesFromNodes(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_faces)
{
	rv_set_faces.Faces().Flush();

	// can loop over subset here
	if(subset_faces.NFaces() == 0)
	{
		subset_faces.Faces() = NaturalNumbers(NFaces());
	}

	for(int i=1; i<=subset_faces.NFaces(); i++)
	{
		int facenr = subset_faces.Faces()(i);
		FEMesh_Face& face = GetFace(facenr);
		int nodesfound = 0;
		for(int j=1; j<=face.NNodes(); j++)
		{
			if( in_set_nodes.Nodes().Find(face.GetNode(j)) ) 
			{
				nodesfound++;
			}
		}
		if( nodesfound == face.NNodes() )  // ALL nodes found 
		{
			rv_set_faces.Faces().Add(facenr);
		}
	}
}

// computes a set of faces for a given set of elements, subset of valid faces can be included
void FEMesh::ComputeFacesOfElements(FEMesh_Set& rv_set_faces, FEMesh_Set& in_set_elements, FEMesh_Set& subset_faces)
{
	rv_set_faces.Faces().Flush();

	// can loop over subset here
	if(subset_faces.NFaces() == 0)
	{
		subset_faces.Faces() = NaturalNumbers(NFaces());
	}

	for(int i=1; i<=subset_faces.NFaces(); i++)
	{
		int facenr = subset_faces.Faces()(i);
		FEMesh_Face& face = GetFace(facenr);

		if( in_set_elements.Elements().Find(face.Element()) )
		{
			rv_set_faces.Faces().Add(facenr);
		}
	}
}

// computes a set of elements for a given set of faces, subset of valis elements can be included
void FEMesh::ComputeElementsOfFaces(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_faces, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();

	for(int i=1; i<=in_set_faces.NFaces(); i++)
	{
		FEMesh_Face& face = GetFace(in_set_faces.Faces()(i));
		rv_set_elements.Elements().AddIfNotExists( face.Element() ); // Add each element only once
	}

	if(subset_elements.NElements())
	{
		rv_set_elements.Elements() &= subset_elements.Elements();
	}
}

// computes a set of all materials used in a given set of elements, subset of valid materials can be included
void FEMesh::ComputeMaterialsOfElements(FEMesh_Set& rv_set_materials, FEMesh_Set& in_set_elements, FEMesh_Set& subset_materials)
{
	rv_set_materials.Materials().Flush();

	for(int i=1; i<=in_set_elements.NElements(); i++)
	{
		FEElement& elem = GetElement(in_set_elements.Elements()(i));
		rv_set_materials.Materials().AddIfNotExists( elem.MaterialNum() );
	}

	if(subset_materials.NMaterials())
	{
		rv_set_materials.Materials() &= subset_materials.Materials();
	}
}

// computes a set of all elements with materials in a given set, subset of valid elements can be included
void FEMesh::ComputeElementsOfMaterials(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_materials, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();

// can loop over subset here
	if(subset_elements.NElements() == 0)
	{
		subset_elements.Elements() = NaturalNumbers(NElements()); 
	}

	for(int i=1; i<=subset_elements.NElements(); i++)
	{
		FEElement& elem = GetElement( subset_elements.Elements()(i) );
		if( in_set_materials.Materials().Find( elem.MaterialNum() ) )
		{
			rv_set_elements.Elements().Add(i);
		}
	}
}

// computes a set of all bodies used in a given set of elements, subset of valid bodies can be included
void FEMesh::ComputeBodiesOfElements(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_elements, FEMesh_Set& subset_bodies)
{
	rv_set_bodies.Bodies().Flush();

	for(int i=1; i<=in_set_elements.NElements(); i++)
	{
		FEElement& elem = GetElement(in_set_elements.Elements()(i));
		rv_set_bodies.Bodies().AddIfNotExists( elem.Domain() );
	}

	if(subset_bodies.NBodies())
	{
		rv_set_bodies.Bodies() &= subset_bodies.Bodies();
	}
}

// computes a set of all elements with bodynumbers in a given set, subset of valid elements can be included
void FEMesh::ComputeElementsOfBodies(FEMesh_Set& rv_set_elements, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();

// can loop over subset here
	if(subset_elements.NElements() == 0)
	{
		subset_elements.Elements() = NaturalNumbers(NElements()); 
	}

	for(int i=1; i<=subset_elements.NElements(); i++)
	{
		FEElement& elem = GetElement( subset_elements.Elements()(i) );
		if( in_set_bodies.Bodies().Find( elem.Domain() ) )
		{
			rv_set_elements.Elements().Add(i);
		}
	}
}

// computes a set of all bodies used in a given set of nodes, subset of valid bodies can be included
void FEMesh::ComputeBodiesOfNodes(FEMesh_Set& rv_set_bodies, FEMesh_Set& in_set_nodes, FEMesh_Set& subset_bodies) 
{
	rv_set_bodies.Bodies().Flush();

	for(int i=1; i<=in_set_nodes.NNodes(); i++)
	{
		FEMesh_Node& node = GetNode(in_set_nodes.Nodes()(i));
		rv_set_bodies.Bodies().AddIfNotExists( node.Domain() );
	}

	if(subset_bodies.NBodies())
	{
		rv_set_bodies.Bodies() &= subset_bodies.Bodies();
	}
}

// computes a set of all nodes with bodynumbers in a given set of bodies
void FEMesh::ComputeNodesOfBodies(FEMesh_Set& rv_set_nodes, FEMesh_Set& in_set_bodies, FEMesh_Set& subset_nodes) 
{
	rv_set_nodes.Nodes().Flush();

// can loop over subset here
	if(subset_nodes.NNodes() == 0)
	{
		subset_nodes.Nodes() = NaturalNumbers(NNodes()); 
	}

	for(int i=1; i<=subset_nodes.NNodes(); i++)
	{
		FEMesh_Node& node = GetNode( subset_nodes.Nodes()(i) );
		if( in_set_bodies.Bodies().Find( node.Domain() ) )
		{
			rv_set_nodes.Nodes().Add(i);
		}
	}
}

// convert a set of any type into a set of desired type
// flag to force clean set otherwise old entries in arrays are kept 
void FEMesh::ConvertSet(FEMesh_Set& the_set, TSetType typei, int flag_flushall)
{
	TSetType the_type = the_set.GetType();
	FEMesh_Set converted(typei);   // the single type result of the convert
	FEMesh_Set temp_elements(TSetElements);
	
	switch (typei)  // type of the converted set
	{
	case TSetNodes: // ComputeNodesOf***
		switch (the_type)
		{
		case TSetNodes: return;
		case TSetElements: ComputeNodesOfElements(converted, the_set); break;
		case TSetFaces: ComputeNodesOfFaces(converted, the_set); break;
		case TSetBodies: ComputeNodesOfBodies(converted, the_set); break;
		case TSetMaterials: ComputeNodesOfMaterials(converted, the_set); break;
		default: break;
		}
		break;
	case TSetElements: // ComputeElementsOf***
		switch (the_type)
		{
		case TSetNodes: ComputeElementsOfNodes(converted, the_set); break;
		case TSetElements: return; break;
		case TSetFaces: ComputeElementsOfFaces(converted, the_set); break;
		case TSetBodies: ComputeElementsOfBodies(converted, the_set); break;
		case TSetMaterials: ComputeElementsOfMaterials(converted, the_set); break;
		default: break;
		}
		break;
	case TSetFaces: // ComputeFacesOf***
		switch (the_type)
		{
		case TSetNodes: ComputeFacesOfNodes(converted, the_set); break;
		case TSetElements: ComputeFacesOfElements(converted, the_set); break;
		case TSetFaces: return; break;
		case TSetBodies: ComputeFacesOfBodies(converted, the_set); break;
		case TSetMaterials: ComputeFacesOfMaterials(converted, the_set); break;
		default: break;
		}
		break;
	case TSetBodies: // ComputeBodiesOf***
		switch (the_type)
		{
		case TSetNodes: ComputeBodiesOfNodes(converted, the_set); break;
		case TSetElements: ComputeBodiesOfElements(converted, the_set); break;
		case TSetFaces: ComputeBodiesOfFaces(converted, the_set); break;
		case TSetBodies: return; break;
		case TSetMaterials: ComputeBodiesOfMaterials(converted, the_set); break;
		default: break;
		}
		break;
	case TSetMaterials: // ComputeMaterialsOf***
		switch (the_type)
		{
		case TSetNodes: ComputeMaterialsOfNodes(converted, the_set); break;
		case TSetElements: ComputeMaterialsOfElements(converted, the_set); break;
		case TSetFaces: ComputeMaterialsOfFaces(converted, the_set); break;
		case TSetBodies: ComputeMaterialsOfBodies(converted, the_set); break;
		case TSetMaterials: return; break;
		default: break;
		}
		break;
	default:
		return;
	}
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ geometry functions - 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// nodes

//get (first) nodenumber for a given position, 0 if no node is found
int FEMesh::GetNodeAtPos(Vector3D pos, FEMesh_Set& subset_nodes, double tol)
{
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set
	
	FEMesh_Set tree_nodes(TSetNodes);
	NodesTree().GetItemsInBox(Box3D(pos,tol),tree_nodes.Nodes());
	if (subset_nodes.NNodes() != 0) tree_nodes &= subset_nodes;

	for (int i=1; i<=tree_nodes.NNodes(); i++)
	{
		if(DistanceFromPoint(tree_nodes(i),pos) < tol ) return tree_nodes(i);
	}
	return 0;
}

// get nodes on a plane (as HNF-equation)
void FEMesh::GetNodesOnPlane(FEMesh_Set& rv_set_nodes, Vector3D nplane, double cplane, FEMesh_Set& subset_nodes , double tol)
{
	rv_set_nodes.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set

	// can loop over subset here
	if(subset_nodes.NNodes() == 0)
	{
		subset_nodes.Nodes() = NaturalNumbers(NNodes()); 
	}

	for (int i=1; i<=subset_nodes.NNodes(); i++)
	{
		double d = DistanceFromPlane(subset_nodes(i),nplane,cplane);
		if( (-tol < d) && (d < tol)) rv_set_nodes.Add(subset_nodes(i));
	}
}

//get a list of nodes that have minimal vaule of distance to plane (min d in HNF-equation)
void FEMesh::GetNodesMinRespVect(FEMesh_Set& rv_set_nodes, Vector3D nplane, FEMesh_Set& subset_nodes, double tol)
{
	rv_set_nodes.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set

	// can loop over subset here
	if(subset_nodes.NNodes() == 0)
	{
		subset_nodes.Nodes() = NaturalNumbers(NNodes()); 
	}
	if(subset_nodes.NNodes() == 0) return; // this line prevents crashes when mesh has no nodes...
	
	double dmin = DistanceFromPlane(subset_nodes(1),nplane,0.0); // first node 

	for (int i=2; i<=subset_nodes.N(); i++)
	{
		double d = DistanceFromPlane(subset_nodes(i),nplane,0.0);
		
		if(((dmin-tol) < d) && (d < (dmin+tol)))
		{
			rv_set_nodes.Add(subset_nodes(i));
		}
		else if(d < (dmin-tol))
		{
			rv_set_nodes.Nodes().Flush();
			rv_set_nodes.Add(subset_nodes(i));
			dmin = d;
		}
	}
}

// sort nodes by their position / coordinate (ascending)
int FEMesh::SortNodesByCoordinate(FEMesh_Set& set_to_be_sorted, int coordinate)
{
  TArray<int> original_nodes(set_to_be_sorted.Nodes());
	TArray<double> coord_vals;
	for (int i=1; i <= set_to_be_sorted.NNodes(); i++)
	{
		double val = GetPoint3D(set_to_be_sorted.Node(i))(coordinate);
		coord_vals.Add(val);
	}
	NaturalNumbers order(set_to_be_sorted.NNodes());

	QuicksortDouble(coord_vals,order);
	
	for(int i=1; i<=set_to_be_sorted.NNodes(); i++)
	{
		set_to_be_sorted.Node(i) = original_nodes(order(i));
	}

	return set_to_be_sorted.NNodes();
}

// get all nodes in a box
void FEMesh::GetNodesInBox(FEMesh_Set& rv_set_nodes, Box3D& box, FEMesh_Set& subset_nodes, double tol)
{
	rv_set_nodes.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set
	
	// $EK 2012-12-12 taking care of tolerance
	box.Increase(tol);

	FEMesh_Set tree_nodes(TSetNodes);
	NodesTree().GetItemsInBox(box, tree_nodes.Nodes());
	RemoveRedundantEntries(tree_nodes.Nodes());
	if (subset_nodes.NNodes() != 0) tree_nodes &= subset_nodes;

	for (int i=1; i<=tree_nodes.N(); i++)
	{ 
		if(box.IsIn(GetPoint3D(tree_nodes(i)))) rv_set_nodes.Add(tree_nodes(i));
	}
}

// computes nodes that are on the circle/ on the disc 
void FEMesh::GetNodesOnCircle(FEMesh_Set& rv_set_nodes_sorted, Vector3D center, Vector3D normal, double radius, int flag_DISC, double radius_hole,FEMesh_Set& subset_nodes, int get_every_X_node, double tol)
{	
  rv_set_nodes_sorted.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set
	
	// can loop over subset here
	if(subset_nodes.NNodes() == 0)
	{
		subset_nodes.Nodes() = NaturalNumbers(NNodes()); 
	}
	if(subset_nodes.NNodes() == 0) return;				// this line prevents crashes when mesh has no nodes...
	if(!normal.Norm()) {assert(0);}								// n must not be (0.0,0.0,0.0)

	FEMesh_Set nodes_temp(TSetNodes);

	Vector3D P,v;
	for (int i=1; i<=subset_nodes.N(); i++)
	{
		P = GetPoint3D(subset_nodes(i));						// get Point of subset
		v = P - center;															// get vector from M to P in subset		
		if(!flag_DISC)
		{
			if( abs(v.Norm()-radius) < tol){					// distance MP = r
				if(abs(v*normal) < tol)				          // vector MP is perpendicular to n
					nodes_temp.Add(subset_nodes(i));			// add node to temporaty list
			}
		}
		else //$AD: entire area
		{
			if( (v.Norm()-radius) < tol){							// distance MP <= r
				if( (radius_hole - v.Norm()) < tol){		// distance MP >= r_hole	//$ DR 2012-08-20
					if(abs(v*normal) < tol)									// vector MP is perpendicular to n
						nodes_temp.Add(subset_nodes(i));			// add node to temporary list
				}
			}
		}
	}
	
	// sort nodes by angle
	if(nodes_temp.N() > 2)
	{	// at least 3 nodes found
		SortNodesOnCircle(nodes_temp,center,normal,radius);

		if(get_every_X_node==1) 
		{
			rv_set_nodes_sorted = nodes_temp;					// return all N found nodes
		}				
		else																				// return only N/get_every_X_node nodes
		{
			for(int i=1; i<=(int)(nodes_temp.N()/get_every_X_node);i++)
			{
				rv_set_nodes_sorted.Add(nodes_temp(i*get_every_X_node));				//nodelist.Add(nodes_neg(i*get_every_X_node));
			}
		}
	}
}

// sorts the nodes by angle, first node defined 0°
void FEMesh::SortNodesOnCircle(FEMesh_Set& rv_set_nodes, Vector3D center, Vector3D normal, double radius)
{
	if(rv_set_nodes.N()>2)
	{	
		TArray<double> cosphi_list, sign_sin, cosphi_list_pos, cosphi_list_neg;
		TArray<int> nodes_pos, nodes_neg;
		double cosphi;
		Vector3D cross,P,v;

		Vector3D P1 = GetPoint3D(rv_set_nodes(1));		// first node is reference (0°)
		Vector3D v1 = P1 - center;										// vector of center to reference point
		cosphi_list.Add(1.0);													// cosphi of reference node = 1
		sign_sin.Add(1.0);														// sign of sinus(phi)
		for(int i=2;i<=rv_set_nodes.N();i++)
		{
			P = GetPoint3D(rv_set_nodes(i));						// get Point P
			v = P - center;															// get vector from M to P	
			cosphi = (v*v1)/((v.Norm())*(v1.Norm()));	  // cos(phi)
			cosphi_list.Add(cosphi);
			cross = v1.Cross(v);
			if(normal(1))				{sign_sin.Add((cross(1)/normal(1)));}	// get sign of sin(phi)
			else if(normal(2))	{sign_sin.Add((cross(2)/normal(2)));}
			else								{sign_sin.Add((cross(3)/normal(3)));}
		}

		QuicksortDouble(sign_sin,cosphi_list,rv_set_nodes.Nodes());		// nodes are sorted according to sign of sin(phi) 

		for(int i=1;i<=rv_set_nodes.N();i++)
		{
			if(sign_sin(i)<0.0){
				nodes_neg.Add(rv_set_nodes(i));
				cosphi_list_neg.Add(cosphi_list(i));
			}
			else
			{
				nodes_pos.Add(rv_set_nodes(i));
				cosphi_list_pos.Add(-cosphi_list(i));
			}
		}

		QuicksortDouble(cosphi_list_neg,nodes_neg);						// all nodes with neg. sin(phi) are sorted
		QuicksortDouble(cosphi_list_pos,nodes_pos);						// all nodes with pos. sin(phi) are sorted
		rv_set_nodes.Nodes() = nodes_neg;
		rv_set_nodes.Nodes() += nodes_pos;
	}
}

// computes nodes that are on the cylinder shell/ in the cylinder
void FEMesh::GetNodesOnCylinder(FEMesh_Set& rv_set_nodes_sorted, Vector3D center_bottom, Vector3D center_top, double radius, int flag_FULL, FEMesh_Set& subset_nodes, int get_every_X_node, double tol)
{
  rv_set_nodes_sorted.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set
	
	// can loop over subset here
	if(subset_nodes.NNodes() == 0)
	{
		subset_nodes.Nodes() = NaturalNumbers(NNodes()); 
	}
	if(subset_nodes.NNodes() == 0) return;				// this line prevents crashes when mesh has no nodes...
	Vector3D axis = center_top - center_bottom;
	if(!axis.Norm()) {assert(0);}									// n must not be (0.0,0.0,0.0)

	FEMesh_Set nodes_temp(TSetNodes);

	double dist, tmpd;
	Vector3D P,v, tmp;
	for (int i=1; i<=subset_nodes.N(); i++)
	{
		P = GetPoint3D(subset_nodes(i));						// get Point of subset
		v = P - center_bottom;											// get vector from M1 to P in subset	

		tmp = v.Cross(axis);
		if(axis.Norm()>0) 
		{
			dist = tmp.Norm() / axis.Norm();					// distance axis-P = r
		}
		else dist = 0;

		if(!flag_FULL)
		{
			if( abs(dist-radius) < tol){							// distance axis-P = r
				tmpd = v*axis;
				dist = tmpd / axis.Norm();
				if( ((dist+tol)>0) && ((dist-tol)< axis.Norm()) )						// axial position of P is between M1 and M2
					nodes_temp.Add(subset_nodes(i));			// add node to list
			}
		}
		else //$AD: entire volume
		{
			if( (dist-radius) < tol){									// distance axis-P <= r
				tmpd = v*axis;
				dist = tmpd / axis.Norm();
				if( ((dist+tol)>0) && ((dist-tol)< axis.Norm()) )						// axial position of P is between M1 and M2
					nodes_temp.Add(subset_nodes(i)); // add node to list
			}
		}
	}

	SortNodesOnCylinder(nodes_temp, center_bottom, center_top, radius);
	//rv_set_nodes_sorted = nodes_temp;	//$ DR 2012-08-13 replaced this line by the following lines (get_every_X_node)

	if(get_every_X_node==1) 
	{
		rv_set_nodes_sorted = nodes_temp;					// return all N found nodes
	}				
	else																				// return only N/get_every_X_node nodes
	{
		for(int i=1; i<=(int)(nodes_temp.N()/get_every_X_node);i++)
		{
			rv_set_nodes_sorted.Nodes().Add(nodes_temp(i*get_every_X_node));		
		}
	}

}

// sorts the nodes by angle and axial position
void FEMesh::SortNodesOnCylinder(FEMesh_Set& rv_set_nodes, Vector3D center_bottom, Vector3D center_top, double radius)
{
	Vector3D axis = center_top - center_bottom;

	// sort nodes according to axial position
	if(!axis.Norm()) {assert(0);}		// n must not be (0.0,0.0,0.0)

	double dist, tmpd;
	Vector3D P,v;
	TArray<double> axial;

	for (int i=1; i<=rv_set_nodes.N(); i++)
	{
		P = GetPoint3D(rv_set_nodes(i));			// get Point of subset
		v = P - center_bottom;								// get vector from M1 to P 	
		
		tmpd = v*axis;
		dist = tmpd / axis.Norm();
		axial.Add(dist);
	}
	QuicksortDouble(axial,rv_set_nodes.Nodes());	

	// detect circles
	TArray<double> diff_axial;
	double max=0.0;

	diff_axial.Add(0.0);
	for (int i=2; i<=axial.Length(); i++)
	{
		diff_axial.Add(axial(i)-axial(i-1));
	}
	
	for (int i=1; i<=axial.Length(); i++)
	{
		if(max <= diff_axial(i)) max = diff_axial(i);
	}

	// sort nodes on circles
	FEMesh_Set nodelist_tmp(TSetNodes), nodelist_sorted(TSetNodes);
	Vector3D M_circle;
	axis.Normalize();
	
	for (int i=1; i<=rv_set_nodes.N() ; i++)
	{
		if ((diff_axial(i) > 0.5*max) || (i == rv_set_nodes.N()))
		{
			if (i == rv_set_nodes.N())				// add last node in list
			{
				nodelist_tmp.Add(rv_set_nodes(i));
			}
			M_circle = center_bottom + axial(i-1) * axis;		// center of circle
			SortNodesOnCircle(nodelist_tmp, M_circle, axis, radius);	// sort nodes on circle
			
			nodelist_sorted += nodelist_tmp;

			nodelist_tmp.Nodes().Flush();
			nodelist_tmp.Add(rv_set_nodes(i));		// add first node of new circle
		}
		else
		{
			nodelist_tmp.Add(rv_set_nodes(i));
		}
	}
	rv_set_nodes = nodelist_sorted;
}

// identify all nodes in a cylindrical shell (node is in shell)
int FEMesh::GetNodesInCylinderShell(FEMesh_Set& rv_set_nodes, Vector3D base, Vector3D axis, double ri, double ro, double hb, double ht, FEMesh_Set& subset_nodes, double tol)
{
	rv_set_nodes.Nodes().Flush();
	ConvertSet(subset_nodes, TSetNodes); // automatically convert the subset to Nodes-Set

	Vector3D a0 = axis;
	a0.Normalize();
	Vector3D m1 = base + hb*a0;
	Vector3D m2 = base + ht*a0;

	FEMesh_Set tree_nodes(TSetNodes);
	Box3D surroundingbox(m1,m2);
	surroundingbox.Increase(ro);
	this->NodesTree().GetItemsInBox(surroundingbox,tree_nodes.Nodes());
	if (subset_nodes.NNodes() != 0) tree_nodes &= subset_nodes;

	FEMesh_Set outer_cyl(TSetNodes), inner_cyl(TSetNodes), inner_rad(TSetNodes);

	GetNodesOnCylinder(outer_cyl, m1, m2, ro, 1, tree_nodes);
	rv_set_nodes = outer_cyl;

	if( ri > tol )
	{
		GetNodesOnCylinder(inner_cyl, m1, m2, ri, 1, tree_nodes);
		GetNodesOnCylinder(inner_rad, m1, m2, ri, 0, tree_nodes);
		rv_set_nodes -= inner_cyl;
		rv_set_nodes += outer_cyl;
	}
	return rv_set_nodes.N();
}

void FEMesh::ScaleMesh(double scaling_factor)	// coordinates of the nodes are scaled: new = scaling_factor*old
{
	Vector3D tmp;
	double max, min;
	for(int i=1; i<=points.Length();i++)
	{
		tmp=points(i).GetCoords3D();
		points(i).SetCoords3D(scaling_factor*tmp);

		if(i==1)
		{
			max = tmp.X(); min = max;
		}
		else
		{
			if(tmp.X()>max) {max=tmp.X();}
			if(tmp.Y()>max) {max=tmp.Y();}
			if(tmp.Z()>max) {max=tmp.Z();}
			if(tmp.X()<min) {min=tmp.X();}
			if(tmp.Y()<min) {min=tmp.Y();}
			if(tmp.Z()<min) {min=tmp.Z();}
		}
	}

//$ YV 2013-01-12: no SolverUO here; another solution is needed if necessary	
//	if(GetMBS()->SolverUO().GetGlobalMessageLevel() >= UO_LVL_dbg1)
	{
		mbs->UO()<< "size of mesh AFTER FEMesh::ScaleMesh: max. value = "<< scaling_factor*max << ", min. value = " << scaling_factor*min <<"\n";
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// elements

// identify all elements in a box
int FEMesh::GetElementsInBox(FEMesh_Set& rv_set_elements, Box3D& box, FEMesh_Set& subset_elements, double tol)
{
	rv_set_elements.Elements().Flush();
	ConvertSet(subset_elements, TSetElements); // automatically convert the subset to Elements-Set

		// can loop over subset here
	if(subset_elements.NElements() == 0)
	{
		subset_elements.Elements() = NaturalNumbers(NElements()); 
	}
	if(subset_elements.NElements() == 0) return 0;				// this line prevents crashes when mesh has no elements...

	for(int i=1; i<=subset_elements.N(); i++)
	{
		Vector3D center = GetCenterOfElement(subset_elements(i));
		if(box.IsIn(center))
			rv_set_elements.Add(subset_elements(i));
	}
	return rv_set_elements.N();
}

// identify all elements in a cylindrical shell (center point is in cylinder)
int FEMesh::GetElementsInCylinderShell(FEMesh_Set& rv_set_elements, Vector3D base, Vector3D axis, double ri, double ro, double hb, double ht, FEMesh_Set& subset_elements, double tol)
{
	rv_set_elements.Elements().Flush();
	ConvertSet(subset_elements, TSetElements); // automatically convert the subset to Elements-Set

	// can loop over subset here
	if(subset_elements.NElements() == 0)
	{
		subset_elements.Elements() = NaturalNumbers(NElements()); 
	}
	if(subset_elements.NElements() == 0) return 0;				// this line prevents crashes when mesh has no elements...
	
	Vector3D a0 = axis;
	a0.Normalize();
	Vector3D m1 = base + hb*a0;
	Vector3D m2 = base + ht*a0;
	
	for(int i=1; i<=subset_elements.N(); i++)
	{
		Vector3D center = GetCenterOfElement(subset_elements(i));
		if(IsPointInCylinder(center, m1, m2, ri, ro))
			rv_set_elements.Add(i);
	}
	return rv_set_elements.N();
}

// checks if a point is in a cylinder shell
int FEMesh::IsPointInCylinder(Vector3D point, Vector3D m1, Vector3D m2, double ri, double ro, double tol)
{
	// distance of a point from a line defined by point & direction || two points
	// d = |(m2-m1)x(m1-p)| / |m2-m1| == |(m2-m1)_0 x (m1-p)|
	Vector3D a0 = m2 - m1;
	double height = a0.Norm();
	a0.Normalize();																	// (m2-m1)_0    unit vector					
  Vector3D v = point - m1;				  							// (m1-p)       vector to point
	Vector3D tmp = v.Cross(a0);										  // ()x()        cross product
	double radial = tmp.Norm();											// ||           vector norm

	double radial2 = ::DistToLine(m1, m2, point);

  if (radial+tol > ro)  return 0;									// distance from axis is too large
	if ((ri > tol) && (radial-tol < ri))  return 0;	// distance from axis is too small

	double axial = v*a0;
	if (axial+tol < 0.) return 0;										// below bottom plane of cylinder
	if (axial-tol > height) return 0;               // above top plane of cylinder
	
	return 1;
}

// returns a surrounding Box3D for a FElement
Box3D FEMesh::GetBox3DofElement(int elemnr)
{
  Vector3D vec = GetPoint3D(GetElement(elemnr).GetNode(1));
	Box3D surround(vec,MESH_STD_TOL);

	for(int i=2; i<=GetElement(elemnr).NNodes(); i++)
	{
		surround.Add(GetPoint3D(GetElement(elemnr).GetNode(i)));
	}
	return surround;
}

// returns a surrounding Box3D for a Set of FElement
Box3D FEMesh::GetBox3DofElements(FEMesh_Set& in_set_elements)
{
	Box3D surround(GetBox3DofElement(in_set_elements(1)));

	for(int i=2; i<=in_set_elements.NElements(); i++)
	{
		surround.Add(GetBox3DofElement(in_set_elements(i)));
	}
	return surround;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// nodes and elements
// finds elements containing the given global node, also store the local node number in Nodes-Array of set
void FEMesh::GetElementsAndLocalNodesOfGlobalNode(FEMesh_Set& rv_set_elements, int in_int_globalnodenumber, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();
	rv_set_elements.Nodes().Flush();
  ConvertSet(subset_elements,TSetElements);

	FEMesh_Set elements_of_globalnode(TSetElements,GetElementsOfNode(in_int_globalnodenumber));
	if(subset_elements.N() != 0)
	{
		elements_of_globalnode &= subset_elements;
	}
	for(int i=1; i<= elements_of_globalnode.N(); i++)
	{
		FEElement& elem = GetElement(elements_of_globalnode(i));
		for(int j=1; j<= elem.NNodes(); j++)
		{
			if(elem.GetNode(j) == in_int_globalnodenumber)
			{
				rv_set_elements.Add(elements_of_globalnode(i),TSetElements); // add the element number to elements array
				rv_set_elements.Add(j,TSetNodes);                            // add the local node number to node array
			}
		}
	}
}

// finds elements containing the given global node, also stores the local coordinates in Vector3D-Array
void FEMesh::GetElementsAndLocalCoordsOfGlobalNode(FEMesh_Set& rv_set_elements, TArray<Vector3D>& rv_array_localcoords, int in_int_globalnodenumber, FEMesh_Set& subset_elements)
{
	Vector3D globPos = GetPoint3D(in_int_globalnodenumber);
	GetElementsAndLocalCoordsOfGlobalPosition(rv_set_elements, rv_array_localcoords, globPos, subset_elements);
}

// finds elements containing the given global position, also stores the local coordinates in Vector3D-Array
void FEMesh::GetElementsAndLocalCoordsOfGlobalPosition(FEMesh_Set& rv_set_elements, TArray<Vector3D>& rv_array_localcoords, Vector3D& in_V3D_globalposition, FEMesh_Set& subset_elements)
{
	rv_set_elements.Elements().Flush();
	rv_array_localcoords.Flush();
  ConvertSet(subset_elements,TSetElements);
	
	if(subset_elements.NElements() == 0)
	{
		subset_elements.SetArray(NaturalNumbers(NElements()));
	}

	for(int i=1; i<=subset_elements.NElements(); i++)
	{
		Box3D box = GetBox3DofElement((subset_elements(i)));
		if(box.IsIn(in_V3D_globalposition))
		{
			rv_set_elements.Add(subset_elements(i));
		}
	}

	for(int i=1; i<= rv_set_elements.NElements(); i++)
	{
		Vector3D lc = GetLocalCoord(rv_set_elements(i),in_V3D_globalposition);
		rv_array_localcoords.Add(lc);
	}
}

// returns first (elementr/localnodenr) pair for a given global node number
int2 FEMesh::GetFirstElementAndLocalNodeOfGlobalNode(int node)
{
  if (nodes_to_elements.Length() != NN()) ComputeNodesToElements();
  int elemnr = nodes_to_elements(node)->Get(1);
	FEElement& elem = this->GetElement(elemnr);
	for (int j=1; j<=elem.NNodes(); j++)
	{
		if(elem.GetNode(j) == node)
		{
			return (int2(elemnr,j));
		}
	}
	return int2();
}

// fills a list with the element nodes global positions
int FEMesh::GetGlobalPositionsOfElementNodes(int elemnr, TArray<Vector3D>& points)
{
	FEElement& elem = GetElement(elemnr);
	points.Flush();
	points.SetLen(elem.NNodes());

	for(int i=1; i <= elem.NNodes(); i++)
	{
		points(i) = GetPoint3D(elem.GetNode(i));
	}
	return points.Length();
}

// returns element number and local coordinate of a global position
int FEMesh::GetElementAndLocalCoordOfGlobalPosition(Vector3D& global, Vector3D& local, IVector& subset)
{
// assume search in all elements if not specified
	if(subset.Length() == 0)
	{
		for(int i=1; i<=NElements(); i++)
			subset.Add(i);
	}
	for(int i=1; i <= subset.Length(); i++)
	{
		FEElement& elem = GetElement(subset(i));	

		int is_on_side = elem.PointIsOnSide(global, GetNodesArray());

		if (is_on_side)
		{
			local = elem.GetLocalPosOnSide(global, GetNodesArray());
			return subset(i);
		}
	}
	return 0; // not found
}

// returns center of an element (global coords)
Vector3D FEMesh::GetCenterOfElement(FEElement& el)
{
	return el.GetCenterPoint(this->GetNodesArray());
}

// calculates local position in respect to a chosen element of a global vector
// ATTENTION: restricted to a global position on the side of the element
Vector3D FEMesh::GetLocalCoord(FEElement& el, Vector3D& globalcoord)
{
// new: FEElement calculates - 
	Vector3D localcoords;
	localcoords = el.GetLocalPosOnSide(globalcoord, GetNodesArray());
	return localcoords;
}

// calculates local position of a global node in the given element
Vector3D FEMesh::GetNodeLocalCoord(FEElement& el, int nodenumber_global)
{
	Vector3D localcoords;
	localcoords = el.GetGlobalNodeLocalCoord(nodenumber_global);
	return localcoords;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// faces

// returns Area or Length of a Side of an Finite Element
double FEMesh::GetFaceArea(FEMesh_Face& feface) 
{
	if (feface.NNodes() == 2) // 2 nodes -> "face" is a line 
	{
		return ::Dist(GetPoint2D(feface.GetNode(1)),GetPoint2D(feface.GetNode(2)));
	}
	else if (feface.NNodes() == 3) // 3 nodes -> face is trig
	{
		return ::Area(GetPoint3D(feface.GetNode(1)),GetPoint3D(feface.GetNode(2)),GetPoint3D(feface.GetNode(3)));
	}
	else if (feface.NNodes() == 4) // 4 nodes -> face is quad == 2 trigs
	{
		return ::Area(GetPoint3D(feface.GetNode(1)),GetPoint3D(feface.GetNode(2)),GetPoint3D(feface.GetNode(3)),GetPoint3D(feface.GetNode(4)));
	}
	return 0.0;
}

// returns Area or Length of a Side of an Finite Element
double FEMesh::GetElementFaceArea(int elem,int side) 
{
	int4 face;
	int len = elements(elem)->GetSideNodeNum(side,face);

	if(len == 2) // 2 nodes -> "face" is line ...
	{
		return ::Dist(GetPoint2D(face(1)),GetPoint2D(face(2)));
	}
	else if (len == 3) // 3 nodes -> face is trig
	{
		return ::Area(GetPoint3D(face(1)),GetPoint3D(face(2)),GetPoint3D(face(3)));
	}
	else if (len == 4) // 4 nodes -> face is quad == 2 trigs
	{
		return ::Area(GetPoint3D(face(1)),GetPoint3D(face(2)),GetPoint3D(face(3)),GetPoint3D(face(4)));
	}
	else mbs->UO().InstantMessageText("ERROR, Element does not have requested Face!\n");
	return 0.0;
}

// calculate the interior angle at a node for a given face
double FEMesh::ComputeAngleAtNode(int nodenumber, int4& face)
{
// always return same weight for all nodes - not weighted by angle
	int leng = face.GetLen();
	double sum_of_interior_angles = MY_PI * (leng-2.);    //(180.0 * (len-2.));
	return sum_of_interior_angles / leng;

// calculate the angle
	Vector3D thispoint, prevpoint, nextpoint;
	
	thispoint = GetPoint3D(face(nodenumber));
	int len = face.GetLen(); // this includes the "0" entries which are sorted to be at the end of the intX

	for(int i=1; i <= len; i++) // first nonzero entry prior to node 
	{
		if(face.GetMod(nodenumber+len-i) > 0)
		{
			prevpoint = GetPoint3D(face.GetMod(nodenumber+len-i));
			break;
		}
	}
	
	for(int i=1; i <= len; i++) // first nonzero entry after node 
	{
		if(face.GetMod(nodenumber+i) > 0)
		{
			nextpoint = GetPoint3D(face.GetMod(nodenumber+i));
			break;
		}
	}

	Vector3D a = prevpoint - thispoint;
	Vector3D b = nextpoint - thispoint; 

	double angle_rad = ::VectorAngle(a,b);
	return angle_rad;
}

// Computes part of the faces ares that is associated with the node (4-gon: node, sidecenter, facecenter, sidecenter)
double FEMesh::ComputeNodeArea(int nodenumber, int4& face)
{
	Vector3D centerpoint = GetFaceCenterPoint(face);

	Vector3D thispoint = GetPoint3D(face(nodenumber));
	Vector3D prevpoint = GetPrevNodePos(nodenumber, face);
	Vector3D nextpoint = GetNextNodePos(nodenumber, face);

	Vector3D h1 = (thispoint+prevpoint)*0.5;
	Vector3D h2 = (thispoint+nextpoint)*0.5;

	double area = fabs(::Area(thispoint, h1, centerpoint, h2));
	return area;
}

// returns number of valid entries in face and fills p1..p4 with nodepositions
int FEMesh::GetFaceNodePositions(int4& face, Vector3D& p1, Vector3D& p2, Vector3D& p3,  Vector3D& p4)
{
	int len = face.RemoveAllRedundantEntries();

	if (face(1) > 0) p1 = GetPoint3D(face(1));
	if (face(2) > 0) p2 = GetPoint3D(face(2));
	if (face(3) > 0) p3 = GetPoint3D(face(3));
	if (face(4) > 0) p4 = GetPoint3D(face(4));
	return len;
}

// returns the center point of the face ( Center of Gravity )
Vector3D FEMesh::GetFaceCenterPoint(int4& face)
{
	Vector3D p1,p2,p3,p4;
	int len = GetFaceNodePositions(face, p1, p2, p3, p4); 
	TArray<int3> trigs;
	TArray<Vector3D> points; 
	points.Set4( p1,p2,p3,p4 );

	if (len <= 1) return p1;                       // point
	else if (len == 2) return (p1+p2)*0.5;         // line
	else if (len == 3)                             // trig
	{
		return (p1+p2+p3)*(1./3.);   // this computation is WRONG (but ComputeCenterOfMassTrigs is wrong too)
	} 
	else                                           // quad
	{
		return (p1+p2+p3+p4)*0.25;   // this computation is WRONG (but ComputeCenterOfMassTrigs is wrong too)
	}
}

// returns the position of the (cyclic) previous node	
Vector3D FEMesh::GetPrevNodePos(int thisnodenr, int4& face)
{
	int len = face.GetLen(); // full length
	for(int i=1; i <= len; i++) // first nonzero entry prior to node 
	{
		if(face.GetMod(thisnodenr+len-i) > 0)
		{
			return GetPoint3D(face.GetMod(thisnodenr+len-i)); // returns "thisnode" in last loop cycle
		}
	}
	return 0;
}

// returns the position of the (cyclic) next node
Vector3D FEMesh::GetNextNodePos(int thisnodenr, int4& face)
{
	int len = face.GetLen(); // full length
	for(int i=1; i <= len; i++) // first nonzero entry after node 
	{
		if(face.GetMod(thisnodenr+i) > 0)
		{
			 return GetPoint3D(face.GetMod(thisnodenr+i));
		}
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// general - vector - distance

// returns 1 if all points are in one plane 
int FEMesh::IsPlanar(Vector3D p1,Vector3D p2,Vector3D p3,Vector3D p4)
{
	Vector3D d2 = p2-p1;
	Vector3D d3 = p3-p1;
	Vector3D d4 = p4-p1;

	Vector3D c23 = d2.Cross(d3);
  double dummy = c23 * d4;
	
	if (fabs(dummy) < 1E-10) return 1;
	else return 0;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ cleanup
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// returns list of nodes new numbers (... a(n) = 0 : delete node)
void FEMesh::FindUnusedNodes(IVector& nodelist)
{
// returns "0" for unused, "newnum" for used nodes
	nodelist.SetLen(NN());
	nodelist.SetAll(0);
	int newnumber=0;

	for(int i=1; i<= NN(); i++)
	{
		if(nodes_to_elements(i)->Length() != 0) 
		{
			newnumber++;
			nodelist(i) = newnumber;
		}
	}
}

// remove nodes that are not used in any element
int FEMesh::RemoveUnusedNodes()
{
 	ComputeNodesToElements();
	IVector nodenums;
	FindUnusedNodes(nodenums);
	int rv = RemoveNodes(nodenums);
	InitializeNodeSearchTree();
 	ComputeNodesToElements();
	return rv;
}


// returns list of nodes new numbers (... a(n+1) <= a(n) : merge node)
void FEMesh::FindDoubleNodes(IVector& nodelist, double tol)
{
// returns "newnum" for new nodes, "i<newnum" for double nodes
	nodelist.SetLen(NN());
	nodelist.SetAll(0);
	int newnumber=0;
	
	InitializeNodeSearchTree();
	for (int i=1; i <= NN(); i++)
	{
		if (nodelist(i) == 0) // not maked as identical in previous loop
		{
			newnumber++;
			nodelist(i) = newnumber;
			IVector subset;
			nodestree.GetItemsInBox(Box3D(GetPoint3D(i),tol),subset);
			for (int j=1; j<=subset.Length(); j++) 
			{	
				if (subset(j)>i) // search upwards
				{
					if (DistanceFromPoint(subset(j),GetPoint3D(i)) < tol) 
						nodelist(subset(j)) = newnumber; 
				}
			}
		}
	}
}


// merge nodes that occupy the same position (with tolerance)
int FEMesh::RemoveDoubleNodes(double tol)
{
	ComputeNodesToElements();
	IVector nodenums;
	FindDoubleNodes(nodenums,tol);
	int rv = RemoveNodes(nodenums);
	InitializeNodeSearchTree();
	return rv;
}

// actually removes the nodes, changes nodenumbers in elements accordingly
int FEMesh::RemoveNodes(IVector& newnodenumbers)
{
	// loop through elements and faces, enter new nodenumbers there
	for (int i=1; i<=elements.Length(); i++)
	{
		for (int j=1; j<=elements(i)->NNodes(); j++)
		{
			elements(i)->SetNode(j,newnodenumbers(elements(i)->GetNode(j)));
		}
	}
	for (int i=1; i<=NFaces(); i++)
	{
		for (int j=1; j<=GetFace(i).NNodes(); j++)
		{
			GetFace(i).SetNode(j,newnodenumbers(GetFace(i).GetNode(j)));
		}
	}
	return this->EraseManyNodes(newnodenumbers);
}

int FEMesh::RemoveUnusedFaces()
{
	// check if faces have removed nodes, erase those faces
	int nf = NFaces();
	TArray<int> mask;
	mask.SetLen(nf);
	mask.SetAll(0);
	for (int i=1; i<=nf; i++)
	{
		FEMesh_Face& feface = GetFace(i);
		int oldlen = feface.NNodes();
		int newlen = feface.RemoveRedundant();
		if (oldlen == newlen)
			mask(i)=0; // keep face
		else
			mask(i)=1; // mark for delete
	}
	Faces().EraseMany(mask);
	return nf-NFaces();

}

// erase single face from the dataarray
void FEMesh::EraseFace(int i)
{
	this->faces.Erase(i);
	InitializeNodeSearchTree();
}
// erase single element from the dataarray
void FEMesh::EraseElem(int i)
{
	this->elements.Erase(i);
	InitializeNodeSearchTree();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CLEANUP DONE UP TO THIS LINE - CONTINUE CLEANUP BELOW !!! 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//////// fills a list (elementr/localnodenr) pairs for a given global node number
//////void FEMesh::GetElementsAndLocalNodeOfGlobalNode(int node, TArray<int2>& list)
//////{
//////	list.Flush();
//////	if (nodes_to_elements.Length() != NN()) ComputeNodesToElements();
//////
//////	for (int i=1; i<=nodes_to_elements(node)->Length(); i++)
//////	{
//////		int elemnr = nodes_to_elements(node)->Get(i);
//////		FEElement& elem = this->GetElement(elemnr);
//////
//////		for (int j=1; j<=elem.NNodes(); j++)
//////		{
//////			if(elem.GetNode(j) == node)
//////			{
//////				list.Add(int2(elemnr,j));
//////				break;
//////			}
//////		}
//////	}
//////}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ mesh processing
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ mesh processing 
//initial nodal conditions 
/*
//add initial rotational velocity field (w)x(r'-r) to all nodes
void FEMesh::AddInitialRotationalVelocity(const Vector3D& omega, const Vector3D& point_at_rot_axis)
{
	if (NN() != GetInitVels3D()->Length()) InitializeInitVelocities();

	Vector3D w = omega;
	Vector3D x0 = point_at_rot_axis;
	Vector3D v,r;

	for (int i=1; i<= NN(); i++)
	{
		r = GetPoint3D(i);
		v = GetInitialNodalVelocity3D(i) + w.Cross(r - x0);
		SetInitialNodalVelocity(i, v);
	}
}

//add initial rotational velocity field (w)x(r'-r) to all nodes, where w is in z-direction
void FEMesh::AddInitialRotationalVelocity(double omegaz, const Vector2D& point_at_rot_axis)
{
	if (NN() != GetInitVels2D()->Length()) InitializeInitVelocities();

	double w = omegaz;
	Vector2D x0 = point_at_rot_axis.X();
	Vector2D v,r,n;
	
	for (int i=1; i<= NN(); i++)
	{
		r = GetPoint2D(i);
		n = Vector2D( -w*(r-x0).Y(), w*(r-x0).X());
		v = GetInitialNodalVelocity2D(i)+n;
		SetInitialNodalVelocity(i, v);
	}
}

//add initial translational velocity 'vel' to all nodes
void FEMesh::AddInitialTranslationalVelocity(const Vector3D& vel)
{
	if (NN() != GetInitVels3D()->Length()) InitializeInitVelocities();
	
	for (int i=1; i<= NN(); i++)
	{
		SetInitialNodalVelocity(i, GetInitialNodalVelocity3D(i) + vel);
	}
}
//add initial translational velocity 'vel' to all nodes
void FEMesh::AddInitialTranslationalVelocity(const Vector2D& vel)
{
	if (NN() != GetInitVels2D()->Length()) InitializeInitVelocities();
	for (int i=1; i<= NN(); i++)
	{
		SetInitialNodalVelocity(i, GetInitialNodalVelocity2D(i) + vel);
	}
}
*/
//+ mesh processing 
//get elements and local nodenumers for global nodenumbers




//+ mesh processing 
//get elements and local position for global position


////// searches for global position in (a subset of) elements and returns list of elements and corresponding local position 
////void FEMesh::GetElementsAndLocalCoordOfGlobalPosition(Vector3D& pos, TArray<int>& elems, TArray<Vector3D>& loccoords, TArray<int>& subset)
////{
////// assume search in all elements if not specified
////	if(subset.Length() == 0)
////	{
////		for(int i=1; i<=NElements(); i++)
////			subset.Add(i);
////	}
////
////	for(int i=1; i<=subset.Length(); i++)
////	{
////		Box3D box = GetBox3DofElement((subset(i)));
////		if(box.IsIn(pos))
////		{
////			elems.Add(subset(i));
////		}
////	}
////
////	for(int i=1; i<=elems.Length(); i++)
////	{
////		Vector3D lc = GetLocalCoord(elems(i),pos);
////		loccoords.Add(lc);
////	}
////}






//+ mesh processing 
//surfaces & interfaces


//computes list of all surfaceelements, fills arrays surfaceelements and surfaceparents, same for interfaces
//void FEMesh::FindSurfaceElements()



//////finds the corresponding surfaceelement for all faces from input file - uses surfaceelements
////void FEMesh::MapFacesToSurfaceElements(int trybothcycles)
////{
////
////	int sel=surfaceelements.Length();
////	if (sel == 0)
////	{
////			mbs->UO(UO_LVL_warn) << "no surfaceelements ! run FindSurfaceElements()\n";
////			return;
////	}
////
////	int fl=faces.Length();
////	if (fl == 0) 
////	{
////		mbs->UO(UO_LVL_warn) << "no FACES file loaded ! use surfaceelements as faces\n";
////		face_to_surfaceelement.SetLen(sel);
////		for (int i=1; i<=sel; i++) face_to_surfaceelement(i) = i;
////	}
////
////	int nnse=0;
////	int nnf=0;
////	face_to_surfaceelement.SetLen(fl);
////	face_to_surfaceelement.SetAll(0);
////
////	TIntX* surelem;
////	TIntX* face;
////	int4 surf4,face4;
////	int3 surf3,face3;
////	int2 surf2,face2;
////
////	for (int i=1; i <= fl; i++)
////	{
////		face_to_surfaceelement(i) = FindSurfaceElement4Face(i,trybothcycles);
////	}
////}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ derived class FEMesh2D:   base class FEMesh 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//pass node data from internal array to mbs
//void FEMesh2D::AddPointsToMBS(int bodyindex)
//{
//	mbsnodenumbers.SetLen(NN());
//	mbsnodenumbers.SetAll(0);
//	
//	// temporal searchtree used by mbs->AddBodyNode
//	ComputeNodeBox();
//	int div = (int) pow(NN()*0.1,1/2.) +1;
//	SearchTree addtree(div,div,1,GetNodeBox());
//
//	for (int i = 1; i <= NN(); i++)
//	{
//		FEMesh_Node fenode = GetNode(i);
//		Node n(2, fenode.Domain(), fenode.GetCoords3D()); // position
//			
//		Vector2D init_d(0.);
//		Vector2D init_v(0.);
//		if(init_disp.Length()) init_d = GetInitialNodalDisplacement2D(i);
//		if(init_vel.Length()) init_v = GetInitialNodalVelocity2D(i);
//	
//		n.SetX_Init(init_d, init_v); // initial conditions
//		mbsnodenumbers(i)=(mbs->AddBodyNode(&n,addtree));
//	}
//}


// special function to load a Mesh created in Netgen
void FEMesh2D::LoadNetgenMesh2D(const mystr& filename, int type, int order, int domain_number) //type: 1==trig, 2==quad; order: 1=linear, 2=quadratic
{
	Reset();
	int np;
	int nelements;

	ifstream file(filename.c_str());
	if(!file.is_open())	//$ DR 2013-03-01 added this error handling
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: could not open file " << filename << "\n";
		return;
	}

	file >> np;
	GetMBS()->UO() << "npoints=" << np << "\n";

	TArray<Vector2D> newpoints;

	for (int i=1; i <= np; i++)
	{
		Vector2D p;
		file >> p.X();
		file >> p.Y();
		newpoints.Add(p);
	}

	file >> nelements;
	GetMBS()->UO() << "nelements=" << nelements << "\n";

	if (order == 1)
	{
		for (int i=1; i <= nelements; i++)
		{
			int domain;
			int p1, p2, p3;
			file >> domain;
			file >> p1;
			file >> p2;
			file >> p3;
			FETrig trig(p1, p2, p3);
			trig.MaterialNum() = domain;	//$!DR 2012-07-09
			//trig.Domain() = 1;	//$ DR 2013-04-11 old code
			if(domain_number) trig.Domain() = domain_number;
			else trig.Domain() = domain;	//$!DR 2013-04-11
			AddElement(trig);
		}
	}
	else if (order == 2)
	{
		for (int i=1; i <= nelements; i++)
		{
			int domain;
			int p1[6];
			file >> domain;
			for (int j = 0; j < 3; j++)
			{
				file >> p1[j];
			}
			
			file >> p1[4];
			file >> p1[5];
			file >> p1[3];
			FETrigquad trig(p1);
			//trig.MaterialNum() = 1;	//$!DR 2012-06-20 HACK
			trig.MaterialNum() = domain;	//$!DR 2012-07-09
			//trig.Domain() = 1;	//$ DR 2013-04-11 old code
			if(domain_number) trig.Domain() = domain_number;
			else trig.Domain() = domain;	//$!DR 2013-04-11
			AddElement(trig);
		}
	}

	//points = newpoints;

	//sort points for better band structure:
	TArray<int> index;
	TArray<int> index2;
	TArray<double> xval;

	for (int i=1; i <= newpoints.Length(); i++)
	{
		xval.Add(newpoints(i).X());
		index.Add(i);
	}

	QuicksortDouble(xval, index); //sort xval, index has reference to original index of xval

	index2.SetLen(index.Length());
	for (int i = 1; i <= index.Length(); i++)
	{
		index2(index(i)) = i;
	}

	
	points.SetLen(newpoints.Length());
	for (int i=1; i <= newpoints.Length(); i++)
	{
		points(i) = newpoints(index(i));
	}

	for (int i = 1; i <= elements.Length(); i++)
	{
		for (int j = 1; j <= elements(i)->NNodes(); j++)
		{
			elements(i)->SetNode(j, index2(elements(i)->GetNode(j)));
		}
	}

}

void FEMesh2D::ComputeBoundary()
{
	//compute boundary edges:
	boundaryedges.SetLen(0);

	GetMBS()->UO() << "nelem=" << elements.Length() << "\n";

	for (int i=1; i <= elements.Length(); i++)
	{
		for (int j = 1; j <= elements(i)->NSides(); j++)
		{
			int2 edge(elements(i)->GetSideNodeNum2(j));
			edge.Swap();
			//GetMBS()->UO() << "edge=" << edge(1) << ", " << edge.(2) << "\n";

			int found = 0;
			for (int i1=1; i1 <= elements.Length(); i1++)
			{ 
				if (i != i1)
				{
					//GetMBS()->UO() << "  edge_" << i1 << ":";

					for (int j1 = 1; j1 <= elements(i1)->NSides2(); j1++)
					{
						int2 edge2(elements(i1)->GetSideNodeNum2(j1));
					
						//GetMBS()->UO() << "(" << edge2(1) << ", " << edge2.(2) << "),";
						//GetMBS()->UO() << "  edge2=" << edge2(1) << ", " << edge2.(2) << "\n";
						if ((edge.Get(1) == edge2.Get(1) && edge.Get(2) == edge2.Get(2)) || 
							(edge.Get(1) == edge2.Get(2) && edge.Get(2) == edge2.Get(1))) found = i1;
					}
					//GetMBS()->UO() << "\n";
				}
				if (found) break;
			}
			if (!found /*&& i < found*/)
			{
				edge.Swap();
				boundaryedges.Add(edge);
			}
		}
	}
	GetMBS()->UO() << "n-boundaryedge=" << boundaryedges.Length() << "\n";
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ derived class FEMesh3D:   base class FEMesh 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//type: 1==tet; order: 1=linear, 2=quadratic
//netgenmode=1: netgen from export mesh (neutral), netgenmode=0: from ANSYS
//void FEMesh3D::LoadNetgenMesh3D(const mystr& filename, int type, int order, double rho0, double Em0, double nu0, int netgenmode) 
void FEMesh3D::LoadNetgenMesh3D(const mystr& filename, int order, int mat_num, int bodynum) //$ DR 2011-03-07
{
	//int netgen = netgenmode;
	Reset();
	//NETGEN neutral file format:
	//number of points
	//X1 Y1 Z1
	//...
	//number of tets
	//domain np1 np2 np3 np4
	//...
	//number of surface elements
	//surface number np1 np2 np3
	//...
	//...

	//$ DR 2012-04-12: read domain number
	if(mat_num)
	{
		GetMBS()->UO() << "LoadNetgenMesh3D: material number and bodynum manually set to " << mat_num << " and " << bodynum <<" for whole mesh.\n";
	}

	int np;
	int nelements;

	ifstream file(filename.c_str());

	file >> np;
	GetMBS()->UO() << "npoints=" << np << "\n";

	TArray<Vector3D> newpoints;

	for (int i=1; i <= np; i++)
	{
		Vector3D p;
		file >> p.X();
		file >> p.Y();
		file >> p.Z();
		newpoints.Add(p);
	}
	/*
	for (int i=1; i <= np; i++)
	{
		GetMBS()->UO() << "p" << i << "=" << newpoints(i) << "\n";
	}*/

	file >> nelements;
	//nelements = 1;
	GetMBS()->UO() << "nelements=" << nelements << "\n";

	if (order == 1)
	{
		for (int i=1; i <= nelements; i++)
		{
			int domain;
			int p1, p2, p3, p4;
			file >> domain;	
			file >> p1;
			file >> p2;
			file >> p3;
			file >> p4;
			FETet tet(p1, p2, p3, p4);

			if(mat_num) //$ DR 2012-04-12: use domain number
			{
				tet.Domain() = bodynum;
				tet.MaterialNum() = mat_num;
			}
			else
			{
				tet.Domain() = domain;
				tet.MaterialNum() = domain;
			}
			AddElement(tet);
		}
	}
	else if (order == 2)
	{
		for (int i=1; i <= nelements; i++)
		{
			double domain;
			int p1[10];
			double p2[10];
			file >> domain;				

			//NETGEN --> HOTINT
			//NETGEN: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
			//HOTINT: 1, 3, 2, 4, 7, 5,10, 6, 9, 8

			file >> p1[1-1];
			file >> p1[3-1];
			file >> p1[2-1];
			file >> p1[4-1];
			file >> p1[7-1];
			file >> p1[5-1];
			file >> p1[10-1];
			file >> p1[6-1];
			file >> p1[9-1];
			file >> p1[8-1];

			FETetquad tet(p1);

			if(mat_num) //$ DR 2012-04-12: use domain number
			{
				tet.Domain() = bodynum;
				tet.MaterialNum() = mat_num;
			}
			else
			{
				tet.Domain() = domain;
				tet.MaterialNum() = domain;
			}
			AddElement(tet);
		}
	}

	 //do not read surfaceelements at the moment
	int nsurfaceelements;
	file >> nsurfaceelements;
	GetMBS()->UO() << "nsurfelements=" << nsurfaceelements << "\n";

//$ AD 2011-11-09: read surfaceelements as FEMesh_Faces..., 
//                 ignoring order = 2..., ignoring surfaceindex...
	for (int i=1; i<= nsurfaceelements; i++)
	{
		int surfaceindex, dummy;
		int3 p;
		file >> surfaceindex;
		file >> p(1);
		file >> p(2);
		file >> p(3);
		if(order == 2)
		{
			file >> dummy;
			file >> dummy;
			file >> dummy;
		}
		FEMesh_Face feface(p, 0, 0);
		int facenr = Faces().Add(feface);
		Surface().Add(facenr);
	}

	//////if (order == 1)
	//////{
	//////	for (int i=1; i <= nsurfaceelements; i++)
	//////	{
	//////		int surfaceindex;
	//////		int p1, p2, p3;
	//////		file >> surfaceindex;
	//////		file >> p1;
	//////		file >> p2;
	//////		file >> p3;
	//////		FETrig trig(p1, p2, p3);
	//////		trig.Domain() = surfaceindex;

	//////		FEElement* elem = trig.GetCopy();
	//////		surfaceelements.Add(elem);
	//////	}
	//////}
	//////else if (order == 2)
	//////{
	//////	for (int i=1; i <= nsurfaceelements; i++)
	//////	{
	//////		int surfaceindex;
	//////		file >> surfaceindex;
	//////		int p1[6];
	//////		for (int j = 0; j < 3; j++)
	//////		{
	//////			file >> p1[j];
	//////		}
	//////		
	//////		file >> p1[4];
	//////		file >> p1[5];
	//////		file >> p1[3];
	//////		FETrigquad trig(p1);

	//////		trig.Domain() = surfaceindex;

	//////		FEElement* elem = trig.GetCopy();
	//////		surfaceelements.Add(elem);
	//////	}
	//////}

	points.Flush();
	for(int i=1; i<= newpoints.Length(); i++)
		points.Add(newpoints(i));
	return;

	//sort points for better band structure:
	TArray<int> index;
	TArray<int> index2;
	TArray<double> xval;

	for (int i=1; i <= newpoints.Length(); i++)
	{
		xval.Add(newpoints(i).X());
		index.Add(i);
	}

	QuicksortDouble(xval, index); //sort xval, index has reference to original index of xval

	index2.SetLen(index.Length());
	for (int i = 1; i <= index.Length(); i++)
	{
		index2(index(i)) = i;
	}

	points.SetLen(newpoints.Length());
	for (int i=1; i <= newpoints.Length(); i++)
	{
		//points(i) = newpoints(index(i));
		points(i) = FEMesh_Node(newpoints(index(i)),bodynum);
		//GetMBS()->UO() << "p" << i << "=" << points(i) << "\n";
	}

	for (int i = 1; i <= elements.Length(); i++)
	{
		for (int j = 1; j <= elements(i)->NNodes(); j++)
		{
			elements(i)->SetNode(j, index2(elements(i)->GetNode(j)));
			elements(i)->MaterialNum() = mat_num;		//$ DR 2011-03-07:[
			elements(i)->Domain() = bodynum;				//$ DR 2011-03-07:]
		}
	}
}

// Load mesh from ANSYS node and element files
//type: 1==tet; 2== hex; order: 1=linear, 2=quadratic
// so far, only linear hexes are supported
void FEMesh3D:: LoadAnsysNodesAndElements(const string& filename_nodes, const string& filename_elements, int type, int order, int bodynum)
{
	if (type != 2) 
	{
		GetMBS()->UO() << "Error in FEMesh3D::LoadANSYSMesh3D: element type " << type << " not supported, only hexes (type 2) work!\n";
		return;
	}
	if (order != 1)
	{
		GetMBS()->UO() << "Error in FEMesh3D::LoadANSYSMesh3D: order " << order << " not supported, only support order 1!\n";
		return;
	}

	// ANSYS element list:
	//n1 n2 n3 n4 n5 n6 n7 n8 matnum dummy dummy dummy dummy elementnumber

	int np = 0;
	int nelements = 0;

	ifstream file_nodes(filename_nodes.c_str());
	ifstream file_elements(filename_elements.c_str());

	TArray<Vector3D> newpoints;	
	TArray<int> nodenums;
	int dummy, readint;
	while( !file_nodes.eof() )
	{
		file_nodes >> readint;
		// test if end of file was reached but for empty lines..
		if (file_nodes.eof()) break;
		nodenums.Add(readint);
		if (readint > np) np = readint;
		Vector3D p;
		file_nodes >> p.X();
		file_nodes >> p.Y();
		file_nodes >> p.Z();

		newpoints.Add(p);
		char ch;
		file_nodes.get(ch);
	}

	Quicksort(nodenums, newpoints);

	points.Flush();
	for(int i=1; i<=newpoints.Length(); i++)
		points.Add(FEMesh_Node(newpoints(i)));
	
	np = newpoints.Length();

	int matnum, elnum = 0;
	Vector3D color(1,0,0);
	while (!file_elements.eof())
	{
		TArray<int> p(8);
		file_elements >> p(1);
		// test if end of file was reached but for empty lines..
		if (file_elements.eof()) break;
		file_elements >> p(2);
		file_elements >> p(4);
		file_elements >> p(3);
		file_elements >> p(5);
		file_elements >> p(6);
		file_elements >> p(8);
		file_elements >> p(7);

		file_elements >> matnum;
		file_elements >> dummy;
		file_elements >> dummy;
		file_elements >> dummy;
		file_elements >> dummy;
		file_elements >> elnum;

		FEHex hexel(p);
		hexel.MaterialNum() = matnum;
		hexel.Domain() = bodynum;

		FEElement* elem = hexel.GetCopy();
		elements.Add(elem);
	}
}

//load files generated with Tool "Gene.mac" from MN
void FEMesh3D::LoadAnsysNodesAndElementsWithCGandMaterial(const mystr& filename_nodes, const mystr& filename_elements)
{
	// ANSYS node list (conversion of MN):
	//nodenumber_i n_i_x n_i_y n_i_z
	//
	Reset();

	ifstream file_nodes(filename_nodes.c_str());

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//READ NODES
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int np = 0; //node count

	double readint;
	while( !file_nodes.eof() )
	{
		file_nodes >> readint; //nodenumber
		// test if end of file was reached but for empty lines..
		if (file_nodes.eof()) break;

		if ((int)readint > np) np = (int)readint;
		Vector3D p;
		file_nodes >> p.X();
		file_nodes >> p.Y();
		file_nodes >> p.Z();

		points.Add(p);
	}

	mbs->UO(UO_LVL_ext) << "read " << np << " nodes from ANSYS\n";
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//READ ELEMENTS
	//
	// structure:
	// Elnum CGx CGy CGz nnodes node_1 node_2 ... node_n volume rho Em nu
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int nelements = 0;
	ifstream file_elements(filename_elements.c_str());
	TArray<int> nodenumbers; //ANSYS
	TArray<int> nodenumbers_hi; //HOTINT

	int warned = 0;

	int cnt = 0;
	while( !file_elements.eof() /*&& cnt++ < 100*/)
	{
		file_elements >> readint; //element number
		// test if end of file was reached but for empty lines..
		if (file_elements.eof()) break;

		if ((int)readint > nelements) nelements = (int)readint;
		
		Vector3D CG;
		file_elements >> CG.X();
		file_elements >> CG.Y();
		file_elements >> CG.Z();

		file_elements >> readint; //number of nodes for element
		int nn = (int)readint;
		nodenumbers.SetLen(0);
		nodenumbers_hi.SetLen(0);

		if (nn > 10) nn = 20; //all elements except tet have 20 nodes ==> everything transformed to Hexahedrals == ANSYS Element 186 with options!

		for (int i=1; i <= nn; i++)
		{
			file_elements >> readint; //node numbers
			nodenumbers.Add((int)readint);
		}
		for (int i=nn+1; i <= 20; i++) //each element has 20 node numbers, some are empty
		{
			file_elements >> readint; //read zero entries
		}
		//mbs->UO(UO_LVL_dbg1) << "element has " << nn << "nodes = " << nodenumbers << "\n";


		double volume, rho, Em, nu;
		file_elements >> volume;
		file_elements >> rho;
		file_elements >> Em;
		file_elements >> nu;

		elementCG.Add(CG);
		elementVol.Add(volume);
		int matnum = AddElasticMaterial(rho, Em, nu);
		
		if (nn == 8)
		{
			nodenumbers_hi.Add(nodenumbers( 4)); //1
			nodenumbers_hi.Add(nodenumbers( 1)); //2
			nodenumbers_hi.Add(nodenumbers( 3)); //3
			nodenumbers_hi.Add(nodenumbers( 2)); //4
			nodenumbers_hi.Add(nodenumbers( 8)); //5
			nodenumbers_hi.Add(nodenumbers( 5)); //6
			nodenumbers_hi.Add(nodenumbers( 7)); //7
			nodenumbers_hi.Add(nodenumbers( 6)); //8

			FEHex hex(nodenumbers_hi);
			hex.MaterialNum() = matnum;
			hex.Domain() = matnum;

			AddElement(hex);
		}
		else if (nn == 20)
		{
	//order of nodes for linear hex:
	//           7      8
	//           +------+
	//          /|     /|
	// Z     5 +------+6|     
	// ^       | |    | |
	// | Y     |3+- - |-+ 4
	// |/      |/     |/
	// --->X 1 +------+ 2     
	//
			nodenumbers_hi.Add(nodenumbers( 4)); //1
			nodenumbers_hi.Add(nodenumbers( 1)); //2
			nodenumbers_hi.Add(nodenumbers( 3)); //3
			nodenumbers_hi.Add(nodenumbers( 2)); //4
			nodenumbers_hi.Add(nodenumbers( 8)); //5
			nodenumbers_hi.Add(nodenumbers( 5)); //6
			nodenumbers_hi.Add(nodenumbers( 7)); //7
			nodenumbers_hi.Add(nodenumbers( 6)); //8
			nodenumbers_hi.Add(nodenumbers(12)); //9
			nodenumbers_hi.Add(nodenumbers(10)); //10
			nodenumbers_hi.Add(nodenumbers(16)); //11
			nodenumbers_hi.Add(nodenumbers(14)); //12
			nodenumbers_hi.Add(nodenumbers(11)); //13
			nodenumbers_hi.Add(nodenumbers( 9)); //14
			nodenumbers_hi.Add(nodenumbers(15)); //15
			nodenumbers_hi.Add(nodenumbers(13)); //16
			nodenumbers_hi.Add(nodenumbers(20)); //17
			nodenumbers_hi.Add(nodenumbers(17)); //18
			nodenumbers_hi.Add(nodenumbers(19)); //19
			nodenumbers_hi.Add(nodenumbers(18)); //20

			FEHexquad hexquad(nodenumbers_hi);
			hexquad.MaterialNum() = matnum;
			hexquad.Domain() = matnum;

			AddElement(hexquad);
		}
		else if (nn == 10)
		{
			nodenumbers_hi.Add(nodenumbers( 1)); //1
			nodenumbers_hi.Add(nodenumbers( 2)); //2
			nodenumbers_hi.Add(nodenumbers( 3)); //3
			nodenumbers_hi.Add(nodenumbers( 4)); //4
			nodenumbers_hi.Add(nodenumbers( 5)); //5
			nodenumbers_hi.Add(nodenumbers( 6)); //6
			nodenumbers_hi.Add(nodenumbers( 7)); //7
			nodenumbers_hi.Add(nodenumbers(10)); //8
			nodenumbers_hi.Add(nodenumbers( 8)); //9
			nodenumbers_hi.Add(nodenumbers( 9)); //10

			FETetquad tetquad(nodenumbers_hi.GetDataPtr());
			tetquad.MaterialNum() = matnum;
			tetquad.Domain() = matnum;

			AddElement(tetquad);
		}
		else
		{
			if (!warned)
			{
				mbs->UO().InstantMessageText(mystr("Number of nodes is ")+mystr(nn)+mystr(" which is a currently not implemented element\nElement ignored!"));
				warned = 1;
			}
		}
	}
}


	// transform the whole mesh
	// x_new = translation + rotation * x_old
	// ATTENTION: only mesh point list is transformed, nodes already added to mbs are not changed
void FEMesh3D:: Transform(const Vector3D& translation, const Matrix3D& rotation, const IVector& subset)
{
	if (subset.Length() == 0)
	{
		for (int i=1; i<=points.Length(); i++)
		{
			Vector3D help(points(i).GetCoords3D());
			points(i).SetCoords3D(translation + rotation*help);
		}
	}
	else
	{
		for (int i=1; i<=subset.Length(); i++)
		{
			Vector3D help(points(subset(i)).GetCoords3D());
			points(subset(i)).SetCoords3D(translation + rotation*help);
		}
	}
}

	// deform the whole Mesh, by nonlinear shearing
	// Coordinate deformed_coord of all points will be changed w.r.t. MathFunction final_shape
	// ATTENTION: only mesh point list is transformed, nodes already added to mbs are not changed
void FEMesh3D::Distort(int independent_coord, int deformed_coord, TArray<MathFunction*> final_shape, TArray<double> ref_values, int deformed_coord2, const IVector& subset_nodes)
{
	// final_shape(i) of final mesh corresponds to ref_values(i) of original mesh
	// ref_values have to be sorted ascending
	// 
	// if shearing should be performed radial, deformed_coord2 is used to set the 2nd deformed coordinate
	int N,j;
	double t,value, old_value;

	Vector old_values;
	Vector3D tmp;
	IVector subset;

	// initialize Vectors
	if (subset_nodes.Length() == 0)
	{ 
		N = points.Length();
		subset.SetLen(N);
		for(int i=1; i<=N; i++)
		{
			subset(i) = i;
		}
	}
	else 
	{
		N = subset_nodes.Length();
		subset.SetLen(N);
		for(int i=1; i<=N; i++)
		{
			subset(i) = subset_nodes(i);
		}
	}

	// get a copy of the old values of the dependent variables
	old_values.SetLen(N);
	for(int i=1; i<=N; i++)
	{
		tmp=points(subset(i)).GetCoords3D();		// get coordinates of old node
		if(deformed_coord2)
		{
			old_values(i) = sqrt(tmp(deformed_coord)*tmp(deformed_coord)+tmp(deformed_coord2)*tmp(deformed_coord2));
		}
		else
		{
			old_values(i) = tmp(deformed_coord);		// copy value of dependent variable of point
		}
	}

	// deform mesh
	int i_out, i_in, move_point;
	for (int i=1; i<=N; i++)
	{
		move_point = 1;
		tmp=points(subset(i)).GetCoords3D();
		t=tmp(independent_coord);										// get value of independent variable of point
		old_value = old_values(i);

		// find valid MathFunction
		for(i_in=0; i_in<= ref_values.Length()-1; i_in++)
		{
			if(old_value < ref_values(i_in+1))
			{
				break;
			}
		}
		if(i_in==0)
		{
			mbs->UO(UO_LVL_warn) << "WARNING in FEMesh3D::NonLinShear: node is not moved, because the deformable coordinate of the node is smaller than the min. ref_value!\n";
			move_point = 0;
		}

		if(old_value == ref_values(i_in)) {i_out = i_in;}		// point is exactly on MathFunction i_in
		else
		{
			if(i_in == ref_values.Length())	// point is outside max. ref_value
			{
				mbs->UO(UO_LVL_warn) << "WARNING in FEMesh3D::NonLinShear: node is not moved, because the deformable coordinate of the node is greater than the max. ref_value!\n";
				move_point = 0;
			}
			else
			{
				i_out = i_in + 1;
			}
		}

		if(move_point)
		{
			// get new value
			double fo=0;	// fo... final value outside, 
			double fi=0;	// fi... final value inside

			if(final_shape(i_out) != NULL)	{ fo = final_shape(i_out)->Evaluate(t);}		// get the value of the outer MathFunction

			if(i_in == i_out) {	value = fo;}	// point is exactly on MathFunction i_in
			else															// interpolation
			{
				if((i_in!=0) && (final_shape(i_in) != NULL)) { fi = final_shape(i_in)->Evaluate(t);} // get the value of the inner MathFunction

				if(fo<=fi)		// outer MathF is smaller than or equal to inner MathF
				{
					mbs->UO(UO_LVL_err) << "ERROR in FEMesh3D::NonLinShear: final_shape(" << i_in << ") >= final_shape(" << i_out << ")! \n";
					return;
				}
				else 
				{
					value = fi + (fo-fi)/(ref_values(i_out)-ref_values(i_in))*(old_value-ref_values(i_in));
				}
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
			points(subset(i)).SetCoords3D(tmp);	// set deformed variable of point
		}
	}
}
		
void FEMesh3D::AddPointsToReferenceFrame(int frameind, int bodyindex)
{
	Rigid3D* rf = (Rigid3D*)GetMBS()->GetElementPtr(frameind);

	for (int i = 1; i <= points.Length(); i++)
	{
		Node n(3, bodyindex, points(i).GetCoords3D());
		int nn = rf->AddNode(&n);
	}
}

void FEMesh3D::DrawElements(const TArray<Vector3D>& node_disp)
{
	TArray<Vector3D> pts; // coordinates only
	for(int i=1; i<=points.Length(); i++)
		pts.Add(points(i).GetCoords3D());

	for (int i=1; i <= NElements(); i++)
	{
		int dom = GetElement(i).Domain();
		Vector3D col(0.,0.,0.8);

//		if (dom <= material_color.Length()) col = material_color(dom);
		if (dom <= material_data.Length()) col = material_data(dom)->GetMaterialColor();

		mbs->SetColor(col);

	  GetElement(i).DrawElement(mbs, pts, node_disp);	
	}
}


























//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ THE SCRAP YARD - these functions sould be eliminated a.s.a.p.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// returns sorted arrays of cutting planes with all box borders and all prefilled entries within limits of overall box
void FindAllCuttingPlanes(TArray<Box3D>& boxes, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut)
{	

//size of overall Box
	Box3D overall;
	overall = boxes(1);
	for(int i=2; i<= boxes.Length(); i++)
	{
		overall.Add(boxes(i));
	}

// cutting planes due to the boxes
	for(int i=1; i<=boxes.Length(); i++)
	{
		xcut.Add(boxes(i).PMin().X());
		xcut.Add(boxes(i).PMax().X());
		ycut.Add(boxes(i).PMin().Y());
		ycut.Add(boxes(i).PMax().Y());
		zcut.Add(boxes(i).PMin().Z());
		zcut.Add(boxes(i).PMax().Z());
	}

	RemoveRedundantEntries(xcut,1);
	RemoveRedundantEntries(ycut,1);
	RemoveRedundantEntries(zcut,1);

// remove all values out of overall box, only cutting planes 
	while (xcut(1) < overall.PMin().X()) xcut.Erase(1);
	while (xcut.Last() > overall.PMax().X()) xcut.Erase(xcut.Length());
	while (ycut(1) < overall.PMin().Y())ycut.Erase(1);
	while (ycut.Last() > overall.PMax().Y())ycut.Erase(ycut.Length());
	while (zcut(1) < overall.PMin().Z())zcut.Erase(1);
	while (zcut.Last() > overall.PMax().Z())zcut.Erase(zcut.Length());
}

// remove all cutting planes within error tolerance of an other cutting plane & change boxes
void ApplyErrorTolerances(TArray<Box3D>& boxes, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut, double tol = 1e-10)
{
	for(int i=2; i<=xcut.Length(); i++)
	{
		if( abs(xcut(i-1) - xcut(i) ) < tol )
		{
			for(int j=1; j<= boxes.Length(); j++)
			{
				Vector3D& pmin = (Vector3D&) boxes(j).PMin(); // nonconst PMin()...
				Vector3D& pmax = (Vector3D&) boxes(j).PMax();

				if( abs( pmin.X() - xcut(i-1) ) < tol ) 
					pmin.X() = xcut(i-1);
				if( abs( pmax.X() - xcut(i-1) ) < tol ) 
					pmax.X() = xcut(i-1);
			
			}
			xcut.Erase(i);
			i--;
		}
	}
	for(int i=2; i<=ycut.Length(); i++)
	{
		if( abs(ycut(i-1) - ycut(i) ) < tol )
		{
			for(int j=1; j<= boxes.Length(); j++)
			{
				Vector3D& pmin = (Vector3D&) boxes(j).PMin(); // nonconst PMin()...
				Vector3D& pmax = (Vector3D&) boxes(j).PMax();

				if( abs( pmin.Y() - ycut(i-1) ) < tol ) 
					pmin.Y() = ycut(i-1);
				if( abs( pmax.Y() - ycut(i-1) ) < tol ) 
					pmax.Y() = ycut(i-1);
			}
			ycut.Erase(i);
			i--;
		}
	}
	for(int i=2; i<=zcut.Length(); i++)
	{
		if( abs(zcut(i-1) - zcut(i) ) < tol )
		{
			for(int j=1; j<= boxes.Length(); j++)
			{
				Vector3D& pmin = (Vector3D&) boxes(j).PMin(); // nonconst PMin()...
				Vector3D& pmax = (Vector3D&) boxes(j).PMax();

				if( abs( pmin.Z() - zcut(i-1) ) < tol ) 
					pmin.Z() = zcut(i-1);
				if( abs( pmax.Z() - zcut(i-1) ) < tol ) 
					pmax.Z() = zcut(i-1);
			}
			zcut.Erase(i);
			i--;
		}
	}
}

// computes field which block occupies the volume (the latter block always gets priority)
void ComputeOccupationArray(TArray<Box3D>& boxes, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut, TArray<int>& occupation)
{
	int xcl = xcut.Length();
	int ycl = ycut.Length();
	int zcl = zcut.Length();

	occupation.Flush();
	occupation.SetLen((xcl-1)*(ycl-1)*(zcl-1));
	occupation.SetAll(0);

	//size of overall Box
	Box3D overall;
	overall = boxes(1);
	for(int i=2; i<= boxes.Length(); i++)
	{
		overall.Add(boxes(i));
	}

	for(int i=1; i <= boxes.Length(); i++) // all blocks
	{
		for(int iz = 1; iz <= zcl-1; iz++)
		{
			for(int iy = 1; iy <= ycl-1; iy++)
			{
				for(int ix = 1; ix <= xcl-1; ix++)
				{
					Vector3D center;
					center.X() = (xcut(ix) + xcut(ix+1)) * 0.5;
					center.Y() = (ycut(iy) + ycut(iy+1)) * 0.5;
					center.Z() = (zcut(iz) + zcut(iz+1)) * 0.5;

					////if(ix == 1) center.X() = (overall.PMin().X() + xcut(1)) * 0.5; // first segment
					////else if (ix == (xcl-1)) center.X() = (xcut.Last() + overall.PMax().X()) * 0.5; // last segment
					////else center.X() = (xcut(ix-1) + xcut(ix)) * 0.5;

					////if(iy == 1) center.Y() = (overall.PMin().Y() + ycut(1)) * 0.5;
					////else if (iy == (ycl+1)) center.Y() = (ycut.Last() + overall.PMax().Y()) * 0.5;
					////else center.Y() = (ycut(iy-1) + ycut(iy)) * 0.5;

					////if(iz == 1) center.Z() = (overall.PMin().Z() + zcut(1)) * 0.5;
					////else if (iz == (zcl+1)) center.Z() = (zcut.Last() + overall.PMax().Z()) * 0.5;
					////else center.Z() = (zcut(iz-1) + zcut(iz)) * 0.5;

					if( boxes(i).IsIn(center) ) // occupy volume section with block
					{
						int idx = 1 + (ix-1) + (iy-1)*(xcl-1) + (iz-1)*(xcl-1)*(ycl-1);
						////int idx = 1 + (ix-1) + (iy-1)*(xcl+1) + (iz-1)*(xcl+1)*(ycl+1);
						occupation( idx ) = i;
					}
				}
			}
		}
	}
}

// writes computed sectioning to the blocks
void SetBlockData(TArray<SegmentedHexBlock*>& blocks, TArray<double>& xcut, TArray<double>& ycut, TArray<double>& zcut, TArray<int>& occupation)
{
	int xcl = xcut.Length();
	int ycl = ycut.Length();
	int zcl = zcut.Length();

	//size of overall Box
	Box3D overall;
	overall = blocks(1)->GetBox();
	for(int i=2; i<= blocks.Length(); i++)
	{
		overall.Add(blocks(i)->GetBox());
	}
	
	for(int i=1; i <= blocks.Length(); i++)
	{
// compute indices range of this block in set
		int xfirst = xcut.Find(blocks(i)->XMin());
		int xlast = xcut.Find(blocks(i)->XMax());

		int yfirst = ycut.Find(blocks(i)->YMin());
		int ylast = ycut.Find(blocks(i)->YMax());

		int zfirst = zcut.Find(blocks(i)->ZMin());
		int zlast = zcut.Find(blocks(i)->ZMax());

		blocks(i)->XCut().CopyFrom(xcut,xfirst,xlast);
		blocks(i)->YCut().CopyFrom(ycut,yfirst,ylast);
		blocks(i)->ZCut().CopyFrom(zcut,zfirst,zlast);
		blocks(i)->SetAllOccupations(0);
// loop over all sections of the current block
		for(int iz = zfirst; iz <= zlast; iz++)
		{
			for(int iy = yfirst; iy <= ylast; iy++)
			{
				for(int ix = xfirst; ix <= xlast; ix++)
				{
					int idx =	1 + (ix-1) + (iy-1)*(xcl-1) + (iz-1)*(xcl-1)*(ycl-1)	;
					int occu = occupation( idx ); 
					if(occu == i)
					{
						blocks(i)->Occupation(ix-xfirst+1, iy-yfirst+1, iz-zfirst+1) = occu; // set occupation
					}
				}
			}
		}
	}
}


// cuts a set of blocks - compute cutting planes, check priority for occupation, fill data (cutting planes, populaiton array) in HexBlock-class
void CutBlocks(TArray<SegmentedHexBlock*>& blocks, TArray<double> xcut, TArray<double> ycut, TArray<double> zcut)
{
	TArray<Box3D> boxes;
	for(int i=1; i<= blocks.Length(); i++)
	{
		boxes.Add(blocks(i)->GetBox());
	}
	//TArray<double> xcut;
	//TArray<double> ycut;
	//TArray<double> zcut;
	TArray<int> occu;

	// hard coded additional cutting planes for symmetry
	//xcut.Add(0.);
	//ycut.Add(0.);
	FindAllCuttingPlanes(boxes, xcut, ycut, zcut);
	ApplyErrorTolerances(boxes, xcut, ycut, zcut, 1e-10);
	ComputeOccupationArray(boxes, xcut, ycut, zcut, occu);
	SetBlockData(blocks, xcut, ycut, zcut, occu);
}


