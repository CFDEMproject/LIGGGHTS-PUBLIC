//#**************************************************************
//#
//# filename:             FEMesh_aux.h
//#
//# author:               Gerstmayr Johannes
//#                       Dorninger Alexander
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
 

#ifndef FEMESH_AUX__H
#define FEMESH_AUX__H

#include "..\WCDriver3D\PlotToolDlg_aux.h" // Color Palette
#include "geomelements.h"
#include "mathfunc.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ HELPER FUNCTIONS without class binding: 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// MEMORYCLEANUP: Flush TArray<T*>, TMatrix<T*>

// delete all data content of a TArray<T*> and Flush()
template <class T>
void ReleaseArray_TemplatePtr(TArray<T*>& release_me)
{ 
	T* a;
	for (int i = 1; i <= release_me.Length(); i++) 
	{
		a = release_me.Get(i);
		if (a != NULL) 
		{
			delete a; 
			a = NULL;
		}
	}
	release_me.Flush();
};

// delete all data content of a TMatrix<T*> and Flush()
template <class T>
void ReleaseMatrix_TemplatePtr(TMatrix<TArray<T*>*>& release_me)
{
	for (int i=1; i<=release_me.NLines(); i++)
	{
		for (int j=1; j<=release_me.NCols(); j++)
		{
			ReleaseArray_TemplatePtr(release_me.Get(i,j));
		}
	}
	release_me.Flush();
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ (base) class FEElement:  derive linear: {FEHex, FETet, FEQuad, FETrig, FELine} and the respective quadratic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodenumbers, functions to generate faces, domain, material, color
//
typedef enum { TFELine = 0, TFETrig = 1, TFEQuad = 2, TFETet = 3, TFEHex = 4, TFEPrism = 5, TFEPyramid = 6 } TFEType;
typedef enum { TFE0D = 0, TFE1D = 1, TFE2D = 2, TFE3D = 3} TFEDim;
typedef enum { TFElinear = 1, TFEquadratic = 2} TFEOrder;
// classes to process finite elements within FEMesh
 class FEElement
{
	static const int type = 0;
	static const int dim = 0;
	static const int order = 0;
	static const int nnodes = 0;
	static const int nadditionalnodes = 0;
	static const int nsides = 0;

public:
	FEElement() { color = Vector3D(-1.0,0.0,0.0); } // default color -red
	virtual FEElement* GetCopy() const { return 0; }
	FEElement(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FEElement& e)
	{
		*this = e;                      
		domain = e.domain;	
		rho = e.rho;
		em = e.em;
		nu = e.nu;
		materialnum = e.materialnum;
		color = e.color;
	}

	virtual int Type() const {return type;} // 0==line, 1==trig, 2==quad, 3==tet, 4==hex
	virtual int Dim() const {return dim;} 
	virtual int Order() const {return order;} //1==linear, 2==quadratic
	virtual int NNodes() const {return nnodes;}
	virtual const int& Domain() const {return domain;}
	virtual int& Domain() {return domain;}
	virtual int GetNode(int i) const {return 0;}
	virtual void SetNode(int i, int pnum) {;};
	virtual int NSides() const {return nsides;}
	virtual int NSides2() const {return 0;}
	virtual int NSides3() const {return 0;}
	virtual int NSides4() const {return 0;}
	virtual int2 GetSide2(int i) const {return int2(0,0);}
	virtual int3 GetSide3(int i) const {return int3(0,0,0);}
	virtual int4 GetSide4(int i) const {return int4(0,0,0,0);}
	virtual int GetSide(int i, TIntX& retval) const { retval = int4(0,0,0,0); return 0;}
	virtual int2 GetSideNodeNum2(int i) const {return int2(0,0);}
	virtual int3 GetSideNodeNum3(int i) const {return int3(0,0,0);}
	virtual int4 GetSideNodeNum4(int i) const {return int4(0,0,0,0);}
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = int4(0,0,0,0); return 0;}
	// $EK 2013-03-05 returns the node number lying between nodes i and j
	// needed for quadratic elements
	virtual int GetIntermediateNodeNum(int i, int j)
	{ return GetNode(GetIntermediateNode(i,j));	}
	virtual int GetIntermediateNode(int i, int j)
	{ return 0;	}
	virtual void Swap(int* i, int* j) {int dummy = *i; *i = *j; *j = dummy;}
	virtual void Invert() {;};

	virtual const double& Rho() const {return rho;}
	virtual double& Rho() {return rho;}
	virtual const double& Em() const {return em;}
	virtual double& Em() {return em;}
	virtual const double& Nu() const {return nu;}
	virtual double& Nu() {return nu;}
	virtual const double& Thickness() const { return thickness; }
	virtual double& Thickness() { return thickness; }
	virtual const int& MaterialNum() const {return materialnum;}
	virtual int& MaterialNum() {return materialnum;}
	virtual const Vector3D& Color() const {return color;}
	virtual Vector3D& Color() {return color;}

	virtual Vector3D GetGlobalNodeLocalCoord(int nodenumber_global); // returns local coordinates of a (global) node, -1 if node is not in the element
	virtual Vector3D GetNodeLocalCoord(int nodenumber); // returns local coordinates of a (local) node,  if node is not in the element return (0.,0.,0.)
	virtual Vector3D GetCenterPoint(const TArray<class FEMesh_Node>& nodepos); //returns the coordinates of the center point of the finite element

	virtual int PointIsOnSide(const Vector3D& globalpos, const TArray<class FEMesh_Node>& nodes); 	// returns sidenumber of finite element if global point is on that side - returns only first side
	virtual Vector3D GetLocalPosOnSide(const Vector3D& globalpos, const TArray<class FEMesh_Node>& nodes); 	// !!HEX ONLY!!: returns local coordinates of a given global position - assuming that this point is on a side on the finite element
	virtual Vector3D GetLocalPos(const Vector3D& globalpos, const TArray<class FEMesh_Node>& nodepos); // !! NOT IMPLEMENTED YET!! returns local coordinates of a given global position

	virtual void DrawElement(MBS* mbs, const TArray<Vector3D>& node_ref_pos, const TArray<Vector3D>& node_disp); // called by DrawSystem

private:
	int domain; //or surface index
	double rho, em, nu;
	double thickness; // for 2D Elements, individual thickness can overwrite mesh.setting
	int materialnum;
	Vector3D color;
	virtual void Abstractor() = 0;					// (AD) use this to make abstract base class
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +  3D Finite Elements:   Tet, Tetquad, Hex, Hexquad
// TO DO: Prism, Pyramid
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class FETet:public FEElement
{
	static const int type =   TFETet;
	static const int dim =    TFE3D;
	static const int order =  TFElinear;
	static const int nnodes = 4;
	static const int nsides = 4;
public:
	FETet():FEElement() {};
	FETet(int p1, int p2, int p3, int p4):FEElement() {node[0] = p1; node[1] = p2; node[2] = p3; node[3] = p4;};
	FETet(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};

	virtual FETet* GetCopy() const
	{
		FETet* ed = new FETet(*this);
		return ed;
	}
	FETet(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FETet& e)
	{
		FEElement::CopyFrom(e);
		node[0] = e.node[0];
		node[1] = e.node[1];
		node[2] = e.node[2];
		node[3] = e.node[3];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return nsides;}
	virtual int3 GetSide3(int i) const 
	{
		switch(i)
		{
		case 1: return int3(1,3,2); break;
		case 2: return int3(1,2,4); break;
		case 3: return int3(2,3,4); break;
		case 4: return int3(3,1,4); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide3(i); return 3; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum3(i); return 3; }

	virtual void Invert() {Swap(&node[2-1],&node[3-1]);}
private:
	int node[nnodes];
	void Abstractor() {;};
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class FETetquad:public FEElement
{
	static const int type =   TFETet;
	static const int dim =    TFE3D;
  static const int order =  TFEquadratic;
	static const int nnodes = 10;
	static const int nsides = 4;
public:
	FETetquad():FEElement() {};
	FETetquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FETetquad(int* p1):FEElement() 
	{
		for (int i = 0; i < 10; i++) node[i] = p1[i];
	};

	virtual FETetquad* GetCopy() const
	{
		FETetquad* ed = new FETetquad(*this);
		return ed;
	}
	FETetquad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FETetquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < 10; i++) node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return nsides;}
	virtual int3 GetSide3(int i) const //maybe extend for all 16 triangulated surfaces?
	{
		switch(i)
		{
		case 1: return int3(1,3,2); break;
		case 2: return int3(1,2,4); break;
		case 3: return int3(2,3,4); break;
		case 4: return int3(3,1,4); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide3(i); return 3; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum3(i); return 3; }
	// $EK 2013-03-04 needed for area loads in case of quadratics
	virtual int GetIntermediateNode(int i, int j)
	{ 
	//order of nodes (from TetrahedralGeneric<FiniteElement3DBase>::GetS0 (FEHexTet.cpp)):
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
		static const int nodeneighbors[] = {0, 5, 7, 10, 
													5, 0, 6, 8, 
													7, 6, 0, 9, 
													10, 8, 9, 0};
		return nodeneighbors[(i-1)*4+j-1];
	}

	virtual void Invert(){Swap(&node[2-1],&node[3-1]); Swap(&node[5-1],&node[7-1]); Swap(&node[8-1],&node[9-1]);}
private:
	int node[nnodes];
	void Abstractor() {;};
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$EK 2012-11-22
class FEPrism:public FEElement
{
	static const int type =   TFEPrism;
	static const int dim =    TFE3D;
	static const int order =  TFElinear;
	static const int nnodes = 6;
	static const int nsides = 5;
public:
	FEPrism():FEElement() {};
	FEPrism(int p1, int p2, int p3, int p4, int p5, int p6):FEElement() 
		{node[0] = p1; node[1] = p2; node[2] = p3; node[3] = p4; node[4] = p5; node[5] = p6;};
	FEPrism(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FEPrism(const FEElement& e) {CopyFrom(e);}

	virtual FEPrism* GetCopy() const
	{
		FEPrism* ed = new FEPrism(*this);
		return ed;
	}

	virtual void CopyFrom(const FEPrism& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < nnodes; ++i)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return 2;}
	virtual int3 GetSide3(int i) const
	{
		switch(i)
		{
		case 4: return int3(1,3,2); break;
		case 5: return int3(4,6,5); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}

	virtual int NSides4() const {return 3;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,2,5,4); break;
		case 2: return int4(1,4,6,3); break;
		case 3: return int4(2,3,6,5); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const 
	{ 
		if (i <= 3)
			retval = GetSide4(i);
		else
			retval = GetSide3(i); 
		return (i<= 3) ? 4 : 3; 
	}
	virtual int GetSideNodeNum(int i, TIntX& retval) const 
	{ 
		if (i <= 3)
			retval = GetSideNodeNum4(i); 
		else
			retval = GetSideNodeNum3(i); 
		return (i<= 3) ? 4 : 3; 
	}
	
	virtual void Invert() 
	{
		Swap(&node[1-1],&node[4-1]);
		Swap(&node[2-1],&node[5-1]); 
		Swap(&node[3-1],&node[6-1]);
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class FEPrismquad ... quadratic prismatic Finite Element
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$EK 2012-11-27
class FEPrismquad:public FEElement
{
	static const int type =   TFEPrism;
	static const int dim =    TFE3D;
	static const int order =  TFEquadratic;
	static const int nnodes = 15;
	static const int nsides = 5;
public:
	FEPrismquad():FEElement() {};
	FEPrismquad(int* p1):FEElement() 
	{
		for (int i = 0; i < NNodes(); ++i)
			node[i] = p1[i];
	};
	FEPrismquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FEPrismquad(const FEElement& e) {CopyFrom(e);}

	virtual FEPrismquad* GetCopy() const
	{
		FEPrismquad* ed = new FEPrismquad(*this);
		return ed;
	}

	virtual void CopyFrom(const FEPrismquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < NNodes(); ++i)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return 2;}
	virtual int3 GetSide3(int i) const
	{
		switch(i)
		{
		case 4: return int3(1,3,2); break;
		case 5: return int3(4,6,5); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}

	virtual int NSides4() const {return 3;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,2,5,4); break;
		case 2: return int4(1,4,6,3); break;
		case 3: return int4(2,3,6,5); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const 
	{ 
		if (i <= 3)
			retval = GetSide4(i);
		else
			retval = GetSide3(i); 
		return (i<= 3) ? 4 : 3; 
	}
	virtual int GetSideNodeNum(int i, TIntX& retval) const 
	{ 
		if (i <= 3)
			retval = GetSideNodeNum4(i); 
		else
			retval = GetSideNodeNum3(i); 
		return (i<= 3) ? 4 : 3; 
	}

	//virtual int GetSideReal(int i, TArray<int>& nodenum)
	//{
	//	nodenum.SetLen(i<=3 ? 8 : 6);
	//	//nodenum
	//	switch(i)
	//	{
	//	case 1: nodenum.Set6(1,7,2,11,5,13); nodenum(7) = 4; nodenum(8) =10; break;
	//	case 2: nodenum.Set6(1,10,4,15,6,12); nodenum(7) = 3; nodenum(8) = 9; break;
	//	case 3: nodenum.Set6(2,8,3,12,6,14); nodenum(7) = 5; nodenum(8) = 11; break;
	//	case 4: nodenum.Set6(1,9,3,8,2,7); break;
	//	case 5: nodenum.Set6(4,15,6,14,5,13); break;
	//	default: assert(0); return 0;
	//	}
	//	return nodenum.Length();
	//}

	/*virtual int GetSideNodeNumReal(int i, TArray<int>& nodenum)
	{
		TArray<int> loc_nodes;
		GetSideReal(i, loc_nodes);
		nodenum.SetLen(loc_nodes.Length());
		for (int i = 1; i <= nodenum.Length(); ++i)
			nodenum(i) = GetNode(loc_nodes(i));
		return nodenum.Length();
	}*/
	virtual int GetIntermediateNode(int i, int j)
	{
		//order of nodes (from PrismGeneric<FiniteElement3DBase>::GetS0 (FEHexTet.cpp)):
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
		static const int nodeneighbors[] = {0, 7, 9, 10, 0, 0,
													7, 0, 8, 0, 11, 0,
													9, 8, 0, 0, 0, 12,
													10, 0, 0, 0, 13, 15,
													0, 11, 0, 13, 0, 14,
													0, 0, 12, 15, 14, 0};
		return nodeneighbors[(i-1)*6+j-1];
	}

	virtual void Invert() 
	{
		Swap(&node[1-1],&node[4-1]);
		Swap(&node[2-1],&node[5-1]); 
		Swap(&node[3-1],&node[6-1]);
		Swap(&node[7-1],&node[13-1]);
		Swap(&node[8-1],&node[14-1]); 
		Swap(&node[9-1],&node[15-1]);
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class FEPyramidquad ... quadratic pyramid Finite Element
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$EK 2013-03-05 - not tested/finished implementation
class FEPyramidquad:public FEElement
{
	static const int type =   TFEPyramid;
	static const int dim =    TFE3D;
	static const int order =  TFEquadratic;
	static const int nnodes = 13;
	static const int nsides = 5;
public:
	FEPyramidquad():FEElement() {};
	FEPyramidquad(int* p1):FEElement() 
	{
		for (int i = 0; i < NNodes(); ++i)
			node[i] = p1[i];
	};
	FEPyramidquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FEPyramidquad(const FEElement& e) {CopyFrom(e);}

	virtual FEPyramidquad* GetCopy() const
	{
		FEPyramidquad* ed = new FEPyramidquad(*this);
		return ed;
	}

	virtual void CopyFrom(const FEPyramidquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < NNodes(); ++i)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return 4;}
	virtual int3 GetSide3(int i) const
	{
		switch(i)
		{
		case 2: return int3(1,2,5); break;
		case 3: return int3(2,4,5); break;
		case 4: return int3(3,5,4); break;
		case 5: return int3(1,5,3); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}

	virtual int NSides4() const {return 1;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,2,4,3); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const 
	{ 
		if (i == 1)
			retval = GetSide4(i);
		else
			retval = GetSide3(i); 
		return (i == 1) ? 4 : 3; 
	}
	virtual int GetSideNodeNum(int i, TIntX& retval) const 
	{ 
		if (i == 1)
			retval = GetSideNodeNum4(i); 
		else
			retval = GetSideNodeNum3(i); 
		return (i == 1) ? 4 : 3; 
	}

	virtual int GetIntermediateNode(int i, int j)
	{
		//order of nodes (from PyramidGeneric<FiniteElement3DBase>::GetS0 (FEHexTet.cpp)):
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
		static const int nodeneighbors[] = {0,  6,  8,  0,  10, 
																				6,  0,  0,  9,  11, 
																				8,  0,  0,  7,  12, 
																				0,  9,  7,  0,  13,
																				10, 11, 12, 13, 0};
		return nodeneighbors[(i-1)*5+j-1];
	}

	virtual void Invert() 
	{
		//what does invert for a pyramid mean??? EK
		//assert(0);
		/*Swap(&node[1-1],&node[4-1]);
		Swap(&node[2-1],&node[5-1]); 
		Swap(&node[3-1],&node[6-1]);*/
		Swap(&node[1-1],&node[4-1]);
		Swap(&node[10-1], &node[13-1]);
		Swap(&node[7-1], &node[8-1]);
		Swap(&node[9-1], &node[6-1]);
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$EK 2012-11-22
class FEPyramid:public FEElement
{
	static const int type =   TFEPyramid;
	static const int dim =    TFE3D;
	static const int order =  TFElinear;
	static const int nnodes = 5;
	static const int nsides = 5;
public:
	FEPyramid():FEElement() {};
	FEPyramid(int p1, int p2, int p3, int p4, int p5):FEElement() 
		{node[0] = p1; node[1] = p2; node[2] = p3; node[3] = p4; node[4] = p5; };
	FEPyramid(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FEPyramid(const FEElement& e) {CopyFrom(e);}

	virtual FEPyramid* GetCopy() const
	{
		FEPyramid* ed = new FEPyramid(*this);
		return ed;
	}

	virtual void CopyFrom(const FEPyramid& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < nnodes; ++i)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides3() const {return 4;}
	virtual int3 GetSide3(int i) const
	{
		switch(i)
		{
		case 2: return int3(1,2,5); break;
		case 3: return int3(2,4,5); break;
		case 4: return int3(3,5,4); break;
		case 5: return int3(1,5,3); break;
		default: assert(0); return int3(0,0,0);
		}
	}
	virtual int3 GetSideNodeNum3(int i) const {return int3(GetNode(GetSide3(i).Get(1)),GetNode(GetSide3(i).Get(2)),GetNode(GetSide3(i).Get(3)));}

	virtual int NSides4() const {return 1;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,2,4,3); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const 
	{ 
		if (i == 1)
			retval = GetSide4(i);
		else
			retval = GetSide3(i); 
		return (i == 1) ? 4 : 3; 
	}
	virtual int GetSideNodeNum(int i, TIntX& retval) const 
	{ 
		if (i == 1)
			retval = GetSideNodeNum4(i); 
		else
			retval = GetSideNodeNum3(i); 
		return (i == 1) ? 4 : 3; 
	}
	
	virtual void Invert() 
	{
		//what does invert for a pyramid mean??? EK
		//assert(0);
		/*Swap(&node[1-1],&node[4-1]);
		Swap(&node[2-1],&node[5-1]); 
		Swap(&node[3-1],&node[6-1]);*/
		Swap(&node[1-1],&node[4-1]);
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class FEHex:public FEElement
{
	static const int type =   TFEHex;
	static const int dim =    TFE3D;
	static const int order =  TFElinear;
	static const int nnodes = 8;
	static const int nsides = 6;
public:
	FEHex():FEElement() {};
	FEHex(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};

	virtual FEHex* GetCopy() const
	{
		FEHex* ed = new FEHex(*this);
		return ed;
	}
	FEHex(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FEHex& e)
	{
		FEElement::CopyFrom(e);
		for (int i=0; i<NNodes(); i++)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides4() const {return nsides;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,5,7,3); break;
		case 2: return int4(2,4,8,6); break;
		case 3: return int4(1,2,6,5); break;
		case 4: return int4(3,7,8,4); break;
		case 5: return int4(1,3,4,2); break;
		case 6: return int4(5,6,8,7); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide4(i); return 4; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum4(i); return 4; }

	virtual void Invert() 
	{
		Swap(&node[1-1],&node[5-1]);
		Swap(&node[2-1],&node[6-1]); 
		Swap(&node[3-1],&node[7-1]);
		Swap(&node[4-1],&node[8-1]);
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};

class FEHexquad:public FEElement
{
	static const int type =   TFEHex;
	static const int dim =    TFE3D;
	static const int order =  TFEquadratic;
	static const int nnodes = 20;
	static const int nsides = 6;
public:
	FEHexquad():FEElement() {};
	FEHexquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};

	virtual FEHexquad* GetCopy() const
	{
		FEHexquad* ed = new FEHexquad(*this);
		return ed;
	}
	FEHexquad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FEHexquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i=0; i<NNodes(); i++)
			node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides4() const {return nsides;}
	virtual int4 GetSide4(int i) const 
	{
		switch(i)
		{
		case 1: return int4(1,5,7,3); break;
		case 2: return int4(2,4,8,6); break;
		case 3: return int4(1,2,6,5); break;
		case 4: return int4(3,7,8,4); break;
		case 5: return int4(1,3,4,2); break;
		case 6: return int4(5,6,8,7); break;
		default: assert(0); return int4(0,0,0,0);
		}
	}
	virtual int4 GetSideNodeNum4(int i) const {return int4(GetNode(GetSide4(i).Get(1)),GetNode(GetSide4(i).Get(2)),GetNode(GetSide4(i).Get(3)),GetNode(GetSide4(i).Get(4)));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide4(i); return 4; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum4(i); return 4; }

// $EK 2013-03-04 returns the quadratic local node for hexahedral element
// needed for quadratic elements -> Area Loads
	virtual int GetIntermediateNode(int i, int j)
	{
	//order of nodes (compare with HexahedralGeneric<FiniteElement3DBase>::GetS0 (FEHexTet.cpp)):	
	//         7      12     8
	//         +------+------+
	//        /|            /|
  //     15+ |        16 + |
	//      /  +19        /  +20
	//    5+------+------+6  |
	//     |   | 11  10  |   |
	//     |  3+------+--|---+4          
	//	 17+  /        18+  / 
	//     | +13         | +14
	//     |/            |/
  //     +------+------+
  //     1      9      2
		//node neighbors have to be symmetric
		static const int nodeneighbors[] = {0,  9, 13,  0, 17,  0,  0,  0,
																				9,  0,  0, 14,  0, 18,  0,  0,
																				13, 0,  0, 10,  0,  0, 19,  0,
																				0, 14, 10,  0,  0,  0,  0, 20,
																				17, 0,  0,  0,  0, 11, 15,  0,
																				0, 18,  0,  0, 11,  0,  0, 16,
																				0,  0, 19,  0, 15,  0,  0, 12,
																				0,  0,  0, 20,  0, 16, 12,  0};
		return nodeneighbors[(i-1)*8+j-1];
	}

	virtual void Invert() 
	{
		Swap(&node[1-1],&node[5-1]); Swap(&node[2-1],&node[6-1]); Swap(&node[3-1],&node[7-1]); Swap(&node[4-1],&node[8-1]); // corners
		Swap(&node[9-1],&node[11-1]); Swap(&node[10-1],&node[12-1]); // x
		Swap(&node[13-1],&node[15-1]); Swap(&node[14-1],&node[16-1]); // y
	}
private:
	int node[nnodes];
	void Abstractor() {;};
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +  2D Finite Elements:   Trig, Trigquad, Quad, Quadquad
class FETrig:public FEElement
{
	static const int type =   TFETrig;
	static const int dim =    TFE2D;
	static const int order =  TFElinear;
	static const int nnodes = 3;
	static const int nsides = 3;
public:
	FETrig():FEElement() {};
	FETrig(int p1, int p2, int p3):FEElement() {node[0] = p1; node[1] = p2; node[2] = p3;};
	FETrig(int3 points):FEElement() {node[0] = points(1); node[1] = points(2); node[2] = points(3);};
	FETrig(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	virtual FETrig* GetCopy() const
	{
		FETrig* ed = new FETrig(*this);
		return ed;
	}
	FETrig(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FETrig& e)
	{
		FEElement::CopyFrom(e);
		node[0] = e.node[0];
		node[1] = e.node[1];
		node[2] = e.node[2];
	}
public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides2() const {return nsides;}
	virtual int2 GetSide2(int i) const {return int2(i,i%nnodes+1);}
	virtual int2 GetSideNodeNum2(int i) const {return int2(GetNode(i),GetNode(i+1));}	
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide2(i); return 2; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum2(i); return 2; }

	virtual void Invert() {Swap( &node[2-1], &node[3-1] );}
private:
	int node[nnodes];
	void Abstractor() {;};
};

class FETrigquad:public FEElement
{	
	static const int type =   TFETrig;
	static const int dim =    TFE2D;
	static const int order =  TFEquadratic;
	static const int nnodes = 6;
	static const int nsides = 3;
public:
	FETrigquad():FEElement() {};
	FETrigquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FETrigquad(int* p1):FEElement() 
	{
		for (int i = 0; i < 6; i++) node[i] = p1[i];
	};

	virtual FETrigquad* GetCopy() const
	{
		FETrigquad* ed = new FETrigquad(*this);
		return ed;
	}
	FETrigquad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FETrigquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < 6; i++) node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides2() const {return nsides;}
	virtual int2 GetSide2(int i) const {return int2(i,i%nnodes+1);}
	virtual int2 GetSideNodeNum2(int i) const {return int2(GetNode(GetSide2(i).Get(1)),GetNode(GetSide2(i).Get(2)));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide2(i); return 2; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum2(i); return 2; }

	virtual void Invert() {Swap(&node[2-1],&node[3-1]); Swap(&node[4-1],&node[6-1]);}
private:
	int node[nnodes];
	void Abstractor() {;};
};

class FEQuad:public FEElement
{
	static const int type =   TFEQuad;
	static const int dim =    TFE2D;
	static const int order =  TFElinear;
	static const int nnodes = 4;
	static const int nsides = 4;
public:
	FEQuad():FEElement() {};
	FEQuad(int p1, int p2, int p3, int p4):FEElement() {node[0] = p1; node[1] = p2; node[2] = p3; node[3] = p4;};
	FEQuad(int4 points):FEElement() {node[0] = points(1); node[1] = points(2); node[2] = points(3); node[3] = points(4);};
	FEQuad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	virtual FEQuad* GetCopy() const
	{
		FEQuad* ed = new FEQuad(*this);
		return ed;
	}
	FEQuad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FEQuad& e)
	{
		FEElement::CopyFrom(e);
		node[0] = e.node[0];
		node[1] = e.node[1];
		node[2] = e.node[2];
		node[3] = e.node[3];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides2() const {return nsides;}
	virtual int2 GetSide2(int i) const {return int2(i,i%nsides+1);}
	virtual int2 GetSideNodeNum2(int i) const {return int2(GetNode(i),GetNode(i+1));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide2(i); return 2;  }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum2(i); return 2; }

	virtual void Invert() {Swap(&node[2-1],&node[4-1]);}
private:
	int node[nnodes];
	void Abstractor() {;};
};

class FEQuadquad:public FEElement
{
	static const int type =   TFEQuad;
	static const int dim =    TFE2D;
	static const int order =  TFEquadratic;
	static const int nnodes = 9; 
	static const int nsides = 4;
public:
	FEQuadquad():FEElement() {};
	FEQuadquad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	FEQuadquad(int* p1):FEElement() 
	{
		for (int i = 0; i < 8; i++) node[i] = p1[i];
	};

	virtual FEQuadquad* GetCopy() const
	{
		FEQuadquad* ed = new FEQuadquad(*this);
		return ed;
	}
	FEQuadquad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FEQuadquad& e)
	{
		FEElement::CopyFrom(e);
		for (int i = 0; i < 8; i++) node[i] = e.node[i];
	}

public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

	virtual int NSides2() const {return nsides;}
	virtual int2 GetSide2(int i) const {return int2(i,i%nsides+1);}
	virtual int2 GetSideNodeNum2(int i) const {return int2(GetNode(i),GetNode(i+1));}
	virtual int GetSide(int i, TIntX& retval) const { retval = GetSide2(i); return 2; }
	virtual int GetSideNodeNum(int i, TIntX& retval) const { retval = GetSideNodeNum2(i); return 2; }

	virtual void Invert() {Swap(&node[2-1],&node[4-1]); Swap(&node[5-1],&node[8-1]); Swap(&node[6-1],&node[7-1]);}
private:
	int node[nnodes];	
	void Abstractor() {;};
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +  1D Finite Elements:   Line, Linequad
class FELine:public FEElement
{
	static const int type =   TFELine;
	static const int dim =    TFE1D;
	static const int order =  TFElinear;
	static const int nnodes = 2;
	static const int nsides = 0;
public:
	FELine():FEElement() {};
	FELine(int p1, int p2):FEElement() {node[0] = p1; node[1] = p2; };
	FELine(int2 points):FEElement() {node[0] = points(1); node[1] = points(2); };
	FELine(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	virtual FELine* GetCopy() const
	{
		FELine* ed = new FELine(*this);
		return ed;
	}
	FELine(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FELine& e)
	{
		FEElement::CopyFrom(e);
		node[0] = e.node[0];
		node[1] = e.node[1];
	}
public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

private:
	int node[nnodes]; 
	void Abstractor() {;};
};

class FELinequad:public FEElement
{
	static const int type =   TFELine;
	static const int dim =    TFE1D;
	static const int order =  TFEquadratic;
	static const int nnodes = 3;
	static const int nsides = 0;
public:
	FELinequad():FEElement() {};
	FELinequad(int p1, int p2, int p3):FEElement() {node[0] = p1; node[1] = p2; node[2] = p3;};
	FELinequad(int2 points):FEElement() {node[0] = points(1); node[1] = points(2); node[2] = points(3);};
	FELinequad(TArray<int>& points):FEElement() 
	{
		for (int i=0; i<NNodes(); i++)
			node[i] = points(i+1);
	};
	virtual FELinequad* GetCopy() const
	{
		FELinequad* ed = new FELinequad(*this);
		return ed;
	}
	FELinequad(const FEElement& e) {CopyFrom(e);}
	virtual void CopyFrom(const FELinequad& e)
	{
		FEElement::CopyFrom(e);
		node[0] = e.node[0];
		node[1] = e.node[1];
		node[2] = e.node[2];
	}
public:
	virtual int Type() const {return type;} 
	virtual int Dim() const {return dim;}
	virtual int Order() const {return order;} 
	virtual int NNodes() const {return nnodes;}
	virtual int NSides() const {return nsides;}
	virtual int GetNode(int i) const {return node[(i-1)%nnodes];}
	virtual void SetNode(int i, int pnum) {node[(i-1)%nnodes] = pnum;}

private:
	int node[nnodes];
	void Abstractor() {;};
};

//+ end FEElement and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class FENode:  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: position, domain number, of a Node in the FEMesh
class FEMesh_Node
{
public: // lifecycle
	FEMesh_Node():coord(0.,0.,0.),domain(0) 
	{	
	}
	FEMesh_Node(FEMesh_Node& other)
	{
		CopyFrom(other);
	}
	~FEMesh_Node()
	{
	}

public: // lifecycle II
	FEMesh_Node(Vector3D coordi, int domaini = 1)
	{	
		SetFEMesh_Node(coordi ,domaini);
	}
	FEMesh_Node(Vector2D coordi, int domaini = 1)
	{
		SetFEMesh_Node(coordi.MakeV3D() ,domaini);
	}
	void CopyFrom(FEMesh_Node& other)
	{
		coord = other.coord;
		domain = other.domain;
	}
	void SetFEMesh_Node(Vector3D coordi, int domaini = 1)
	{	
		coord = coordi;
		domain = domaini;
	}

public: // access
	virtual const Vector3D& GetCoords3D() const { return coord; }
	virtual void SetCoords3D(const Vector3D coordi) { coord = coordi; }
	virtual void SetCoords3D(const Vector2D coordi) { coord = Vector3D(coordi.X(),coordi.Y(),0.0); }
	virtual Vector3D& Coords3D() { return coord; } 

	virtual const Vector2D GetCoords2D() const { return coord.MakeV2D(); }
	virtual void SetCoords2D(const Vector3D coordi) { coord = Vector3D(coordi.X(),coordi.Y(),0.0); }
	virtual void SetCoords2D(const Vector2D coordi) { coord = Vector3D(coordi.X(),coordi.Y(),0.0); }

	virtual const int& Domain() const { return domain; }
	virtual int& Domain() { return domain; }

	virtual double X() { return coord.X(); }
	virtual double Y() { return coord.Y(); }
	virtual double Z() { return coord.Z(); }

public: // functionality

private: // variables
	Vector3D coord;
	int domain;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class FEMesh_Face:  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: list of nodes, element/side of a face in the FEMesh
class FEMesh_Face
{
public: // lifecycle
	FEMesh_Face():nodelist(0,0,0,0),len(0),element(0),side(0) 
	{	
	}
	FEMesh_Face(FEMesh_Face& other)
	{
		CopyFrom(other);
	}
	~FEMesh_Face()
	{
	}

public: // lifecycle II
	FEMesh_Face(TIntX& nodelisti, int elementi, int sidei)
	{	
		SetFEMesh_Face(nodelisti, elementi, sidei);
	}
	void CopyFrom(FEMesh_Face& other)
	{
		nodelist = other.nodelist;
		len = other.len;
		element = other.element;
		side = other.side;
	}
	void SetFEMesh_Face(TIntX& nodelisti, int elementi, int sidei)
	{	
		nodelist = nodelisti;
		int i = nodelisti.GetLen();
		// $EK 2013-02-13 initialize remaining part of nodelist
		for (int j = i+1; j <= nodelist.GetLen(); ++j)
			nodelist(j) = 0;
		//the length of nodelist is determined properly
		RemoveRedundant();
		element = elementi;
		side = sidei;
	}

public: // access
	virtual const TIntX& GetNodes() const { return nodelist; }
	virtual void SetNodes(TIntX& nodelisti) { nodelist = nodelisti; }
	virtual TIntX& Nodes() { return nodelist; }

	virtual const int GetNode(int i) const { return nodelist(i); }
	virtual void SetNode(int i, int nn) { nodelist(i) = nn; }
	virtual int& Node(int i) {return nodelist(i); }

	virtual const int NNodes() const { return len; }
//	virtual const int GetLen() const { return len; }

	virtual const int GetElement() const { return element; }
	virtual void SetElement(int elementi) { element = elementi; }
	virtual int& Element() { return element; }
	
	virtual const int GetSide() const { return side; }
	virtual void SetSide(int sidei) { side = sidei; }
	virtual int& Side() { return side; }

public: // functionality
	virtual int RemoveRedundant() { len = nodelist.RemoveAllRedundantEntries(); return len; }

	virtual int IsCyclicEqual(FEMesh_Face& other) // uses the comparison routines of intX, typecasts in here...
	{
// assume that any redundant entries were removed previously --> len is the correct length
		if(NNodes() != other.NNodes())
			return 0; // length not equal, can never be same	
		if(len == 3)
		{
			return ((int3)GetNodes()).IsCyclicEqual( (int3) other.GetNodes() ); 
		}
		else if(len == 2)
		{
			return ((int2)GetNodes()).IsCyclicEqual( (int2) other.GetNodes() ); 
		}
		else //(len == 4)
		{
			return GetNodes().IsCyclicEqual( other.GetNodes() ); 
		}
	}

	virtual void InvertRotationSense()
	{
// assume that any redundant entries were removed previously --> len is the correct length
		if(len == 2) 
		{
			int2 tmp = (int2) GetNodes();
			tmp.Invert();
			SetNodes( tmp );
		}
		if(len == 3)
		{
			int3 tmp = (int3) GetNodes();
			tmp.Invert();
			SetNodes( tmp );
		}
		else // len == 4
			Nodes().Invert();
	}
private: // variables
	int4 nodelist;     // corner-NODES of the face - entries set to ZERO when not used
	int len;           // number of nodes used for face
	int element;       // number of the element
	int side;          // number of the side of the face within the element
};

typedef enum {TAreaLoad = 1, TAreaConstraint = 2, TNodalConstraint = 3, TBodyLoad = 4, TFaceLoad = 5, TFaceConstraint = 6, TAreaContact = 7} TLoadType;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class FEMesh_Load: base class for all Loads
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: !type-proxy!, element number, load vector,  of a Load in the FEMesh
// type is set in derived classes
class FEMesh_Load
{
public: // lifecycle
	FEMesh_Load():type(0),element(0),loadvector(0),steps(0) 
	{	
	}
	FEMesh_Load(FEMesh_Load& other)
	{
		CopyFrom(other);
	}
	~FEMesh_Load()
	{
	}

public: // lifecycle II
	FEMesh_Load(int typei, int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi = TArray<StepSettings>(0))
	{	
		SetFEMesh_Load(typei, elementi, loadvectori, stepsi);
	}
	void CopyFrom(FEMesh_Load& other)
	{
		type = other.type;
		element = other.element;
		loadvector = other.loadvector;
		steps.CopyFrom(other.steps);
	}
	virtual FEMesh_Load* GetCopy() { return 0; } // override in all derived classes
	void SetFEMesh_Load(int typei, int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi)
	{	
		type = typei;
		element = elementi;
		loadvector = loadvectori;
		steps.CopyFrom(stepsi);
	}

public: // access
	virtual const int GetType() const { return type; }
	virtual void SetType(int typei) { type = typei; }
	virtual int& Type() { return type; }

	virtual const int GetElement() const { return element; }
	virtual void SetElement(int elementi) { element = elementi; }
	virtual int& Element() { return element; }
	
	virtual const Vector3D GetLoadVector() const { return loadvector; }
	virtual void SetLoadVector(Vector3D loadvectori) { loadvector = loadvectori; }
	virtual Vector3D& LoadVector() { return loadvector; }

	virtual const TArray<StepSettings> GetSteps() const { return steps; }
	virtual void SetSteps(TArray<StepSettings> &stepsi) { steps.CopyFrom(stepsi); }
	virtual TArray<StepSettings>& Steps() { return steps; }

public: // functionality

private: // variables
  int type;            // proxy for load type
	int element;         // number of the *** the Load is applied on, can be elementnr, facenr, bodynr 
	Vector3D loadvector; // strength of the load
	TArray<StepSettings> steps; // 
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Load on a Body, can be set to gravity or manual...
class FEMesh_BodyLoad : public FEMesh_Load 
{
public: // lifecycle
	FEMesh_BodyLoad():FEMesh_Load(),flaggravity(0) { SetType(TBodyLoad); }
	FEMesh_BodyLoad(FEMesh_BodyLoad& other) { CopyFrom(other); }
	~FEMesh_BodyLoad() {}
public: // lifecycle II
	FEMesh_BodyLoad(int elementi, Vector3D loadvectori, int flaggravityi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) {	SetFEMesh_BodyLoad(elementi, loadvectori, flaggravityi, stepsi); }
	void CopyFrom(FEMesh_BodyLoad& other) { FEMesh_Load::CopyFrom(other); flaggravity = other.flaggravity; }
	FEMesh_BodyLoad* GetCopy()  { FEMesh_BodyLoad* l = new FEMesh_BodyLoad(); l->CopyFrom(*this); return l; }
	void SetFEMesh_BodyLoad(int elementi, Vector3D loadvectori, int flaggravityi, TArray<StepSettings>& stepsi) { FEMesh_Load::SetFEMesh_Load(TBodyLoad, elementi, loadvectori, stepsi); flaggravity = flaggravityi; }
public: // access
	virtual const int GetGravity() const { return flaggravity; }
	virtual void SetGravity(int flaggravityi) { flaggravity = flaggravityi; }
	virtual int& Gravity() { return flaggravity; }
private: // own variables
	int flaggravity;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Load applied on an area (will be split into loads on faces and further to loads on individual nodes)
class FEMesh_AreaLoad : public FEMesh_Load  
{
public: // lifecycle
	FEMesh_AreaLoad():FEMesh_Load() { SetType(TAreaLoad); }
	FEMesh_AreaLoad(FEMesh_AreaLoad& other) { CopyFrom(other); }
	~FEMesh_AreaLoad() {}
public: // lifecycle II
	FEMesh_AreaLoad(int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi = TArray<StepSettings>(0), int pressure = 1) {	SetFEMesh_AreaLoad(elementi, loadvectori, stepsi, pressure); }
	void CopyFrom(FEMesh_AreaLoad& other) { FEMesh_Load::CopyFrom(other); flag_pressure = other.flag_pressure; }
	FEMesh_AreaLoad* GetCopy() { FEMesh_AreaLoad* l = new FEMesh_AreaLoad(); l->CopyFrom(*this); return l; }
	void SetFEMesh_AreaLoad(int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi, int pressure) { FEMesh_Load::SetFEMesh_Load(TAreaLoad, elementi, loadvectori, stepsi); flag_pressure= pressure; }
public: // access
	virtual const int GetFlagPressure() const { return flag_pressure; }
	virtual void SetFlagPressure(int flag_pressurei) { flag_pressure = flag_pressurei; }
	virtual int& FlagPressure() { return flag_pressure; }
private: // own variables
	int flag_pressure;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Load applied on a face (will be split into loads on individual nodes)
class FEMesh_FaceLoad : public FEMesh_Load  
{
public: // lifecycle
	FEMesh_FaceLoad():FEMesh_Load() { SetType(TFaceLoad); }
	FEMesh_FaceLoad(FEMesh_FaceLoad& other) { CopyFrom(other); }
	~FEMesh_FaceLoad() {}
public: // lifecycle II
	FEMesh_FaceLoad(int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) {	SetFEMesh_FaceLoad(elementi, loadvectori, stepsi); }
	void CopyFrom(FEMesh_FaceLoad& other) { FEMesh_Load::CopyFrom(other); }
	FEMesh_FaceLoad* GetCopy() { FEMesh_FaceLoad* l = new FEMesh_FaceLoad(); l->CopyFrom(*this); return l; }
	void SetFEMesh_FaceLoad(int elementi, Vector3D loadvectori, TArray<StepSettings>& stepsi) { FEMesh_Load::SetFEMesh_Load(TFaceLoad, elementi, loadvectori, stepsi); }
public: // access
private: // own variables
};

//+ end FEMesh_Load and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class FEMesh_Constraint: base class for all Constraints
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: !type-proxy!, element number, penalty/lagrange flag, stiffness, damping   of a Constraint in the FEMesh
class FEMesh_Constraint
{
public: // lifecycle
	FEMesh_Constraint():type(0), element(0),axis(0),penalty(0),stiffness(0.,0.,0.),damping(0.,0.,0.),steps(0)
	{	
	}
	FEMesh_Constraint(FEMesh_Constraint& other)
	{
		CopyFrom(other);
	}
	~FEMesh_Constraint()
	{
	}

public: // lifecycle II
	FEMesh_Constraint(int typei, int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0))
	{	
		SetFEMesh_Constraint(typei, elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi);
	}
	void CopyFrom(FEMesh_Constraint& other)
	{
		type = other.type;
		element = other.element;
		axis = other.axis;
		penalty = other.penalty;
		stiffness = other.stiffness;
		damping = other.damping;
		steps.CopyFrom(other.steps);
	}
	virtual FEMesh_Constraint* GetCopy() { return 0; } // override in all derived classes
	void SetFEMesh_Constraint(int typei, int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0))
	{	
		type = typei;
		element = elementi;
		axis = axisi;
		penalty = penaltyi;
		stiffness = stiffnessi;
		damping = dampingi;
		steps.CopyFrom(stepsi);
	}

public: // access
	virtual const int GetType() const { return type; }
	virtual void SetType(int typei) { type = typei; }
	virtual int& Type() { return type; }

	virtual const int GetElement() const { return element; }
	virtual void SetElement(int elementi) { element = elementi; }
	virtual int& Element() { return element; }
	
	virtual const int GetAxis() const { return axis; }
	virtual void SetAxis(int axisi) { axis = axisi; }
	virtual int& Axis() { return axis; }

	virtual const int GetPenalty() const { return penalty; }
	virtual void SetPenalty(int penaltyi) { penalty = penaltyi; }
	virtual int& Penalty() { return penalty; }

	virtual const Vector3D GetStiffness() const { return stiffness; }
	virtual void SetStiffness(Vector3D stiffnessi) { stiffness = stiffnessi; }
	virtual void SetStiffness(double stiffnessi) { stiffness = Vector3D(stiffnessi,stiffnessi,stiffnessi); }
	virtual Vector3D& Stiffness() { return stiffness; }

	virtual const Vector3D GetDamping() const { return damping; }
	virtual void SetDamping(Vector3D dampingi) { damping = dampingi; }
	virtual void SetDamping(double dampingi) { damping = Vector3D(dampingi,dampingi,dampingi); }
	virtual Vector3D& Damping() { return damping; }

	virtual const TArray<StepSettings> GetSteps() const { return steps; }
	virtual void SetSteps(TArray<StepSettings> &stepsi) { steps.CopyFrom(stepsi); }
	virtual TArray<StepSettings>& Steps() { return steps; }

public: // functionality

private: // variables
	int type; 
	int element; 
	int axis;
  int penalty;
	Vector3D stiffness;
	Vector3D damping;
	TArray<StepSettings> steps;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Contstraint applied to an Area (will be split into constraint on individual nodes)
class FEMesh_AreaConstraint : public FEMesh_Constraint  
{
public: // lifecycle
	FEMesh_AreaConstraint():FEMesh_Constraint() { SetType(TAreaConstraint); }
	FEMesh_AreaConstraint(FEMesh_AreaConstraint& other) { CopyFrom(other); }
	~FEMesh_AreaConstraint() {}
public: // lifecycle II
	FEMesh_AreaConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_AreaConstraint(elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
	FEMesh_AreaConstraint(int elementi, int axisi, int penaltyi, double stiffnessi, double dampingi) { SetFEMesh_AreaConstraint(elementi, axisi, penaltyi, Vector3D(stiffnessi), Vector3D(dampingi)); }
	void CopyFrom(FEMesh_AreaConstraint& other) { FEMesh_Constraint::CopyFrom(other); }
	FEMesh_AreaConstraint* GetCopy() { FEMesh_AreaConstraint* c = new FEMesh_AreaConstraint(); c->CopyFrom(*this); return c; } 
	void SetFEMesh_AreaConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_Constraint(TAreaConstraint, elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
public: // access
private: // own variables
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Contstraint applied to an Face (will be split into constraint on individual nodes)
class FEMesh_FaceConstraint : public FEMesh_Constraint  
{
public: // lifecycle
	FEMesh_FaceConstraint():FEMesh_Constraint() { SetType(TFaceConstraint); }
	FEMesh_FaceConstraint(FEMesh_FaceConstraint& other) { CopyFrom(other); }
	~FEMesh_FaceConstraint() {}
public: // lifecycle II
	FEMesh_FaceConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_FaceConstraint(elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
	FEMesh_FaceConstraint(int elementi, int axisi, int penaltyi, double stiffnessi, double dampingi) { SetFEMesh_FaceConstraint(elementi, axisi, penaltyi, Vector3D(stiffnessi), Vector3D(dampingi)); }
	void CopyFrom(FEMesh_FaceConstraint& other) { FEMesh_Constraint::CopyFrom(other); }
	FEMesh_FaceConstraint* GetCopy() { FEMesh_FaceConstraint* c = new FEMesh_FaceConstraint(); c->CopyFrom(*this); return c; } 
	void SetFEMesh_FaceConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_Constraint(TFaceConstraint, elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
public: // access
private: // own variables
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Constraint applied to a single Node
class FEMesh_NodeConstraint : public FEMesh_Constraint  
{
public: // lifecycle
	FEMesh_NodeConstraint():FEMesh_Constraint() { SetType(TNodalConstraint); }
	FEMesh_NodeConstraint(FEMesh_NodeConstraint& other) { CopyFrom(other); }
	~FEMesh_NodeConstraint() {}
public: // lifecycle II
	FEMesh_NodeConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_NodeConstraint(elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
	FEMesh_NodeConstraint(int elementi, int axisi, int penaltyi, double stiffnessi, double dampingi) { SetFEMesh_NodeConstraint(elementi, axisi, penaltyi, Vector3D(stiffnessi), Vector3D(dampingi)); }
	void CopyFrom(FEMesh_NodeConstraint& other) { FEMesh_Constraint::CopyFrom(other); }
	FEMesh_NodeConstraint* GetCopy() { FEMesh_NodeConstraint* c = new FEMesh_NodeConstraint(); c->CopyFrom(*this); return c; } 
	void SetFEMesh_NodeConstraint(int elementi, int axisi, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_Constraint(TNodalConstraint, elementi, axisi, penaltyi, stiffnessi, dampingi, stepsi); }
public: // access
private: // own variables
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Contact between two areas - all nodes from "sending" area constrained to "target" area with spherical joints
class FEMesh_AreaContact : public FEMesh_Constraint  
{
public: // lifecycle
	FEMesh_AreaContact():FEMesh_Constraint(),target(0) { SetType(TAreaContact); }
	FEMesh_AreaContact(FEMesh_AreaContact& other) { CopyFrom(other); }
	~FEMesh_AreaContact() {}
public: // lifecycle II
	FEMesh_AreaContact(int elementi, int targeti, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_AreaContact(elementi, targeti, penaltyi, stiffnessi, dampingi, stepsi); }
	FEMesh_AreaContact(int elementi, int targeti, int penaltyi, double stiffnessi, double dampingi) { SetFEMesh_AreaContact(elementi, targeti, penaltyi, Vector3D(stiffnessi), Vector3D(dampingi)); }
	void CopyFrom(FEMesh_AreaContact& other) { FEMesh_Constraint::CopyFrom(other); target = other.target; }
	FEMesh_AreaContact* GetCopy() { FEMesh_AreaContact* c = new FEMesh_AreaContact(); c->CopyFrom(*this); return c; } 
	void SetFEMesh_AreaContact(int elementi, int targeti, int penaltyi, Vector3D& stiffnessi, Vector3D& dampingi, TArray<StepSettings>& stepsi = TArray<StepSettings>(0)) { SetFEMesh_Constraint(TAreaContact, elementi, 0/*axis=0*/, penaltyi, stiffnessi, dampingi, stepsi); target = targeti; }
public: // access
	virtual const int GetTarget() const { return target; }
	virtual void SetTarget(int targeti) { target = targeti; }
	virtual int& Target() { return target; }
private: // own variables
	int target;
};

//+ end FEMesh_Constraints and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#define SetTypesAsP2_not
#ifdef SetTypesAsP2
// SetTypes defined as powers of 2 - easier to define and use mixed entity sets ( but that is not used by convert functions )
typedef enum { TSetUnknown = 0, TSetNodes = 1, TSetElements = 2, TSetFaces = 4, TSetBodies = 8, TSetMaterials = 16} TSetType;
#else
// SetTypes defined as integers - asserting the a set has a MAIN entity, other may or may not be present
typedef enum { TSetUnknown = 0, TSetElements = 1, TSetNodes = 2, TSetFaces = 3, TSetBodies = 4, TSetMaterials = 5} TSetType;
#endif
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class Set: manage all kinds of sets of meshelements
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: type, several Integer-Arrays
class FEMesh_Set
{
public: // lifecycle
	FEMesh_Set(): nodenrs(0), elementnrs(0), facenrs(0), bodynrs(0), materialnrs(0), type(TSetUnknown), name("") {;}
	FEMesh_Set(const FEMesh_Set& other) { CopyFrom(other); }
	~FEMesh_Set() {}
public: // lifecycle II - often used
	FEMesh_Set(TSetType typei, TArray<int> listi, mystr namei) { SetFEMesh_Set(typei, listi, namei); } // normal initialization
	FEMesh_Set(TSetType typei, TArray<int> listi) { SetFEMesh_Set(typei, listi, mystr("")); }          // lacking name
	FEMesh_Set(TSetType typei, mystr namei) { SetFEMesh_Set(typei, TArray<int>(0), namei); }           // lacking data
	FEMesh_Set(TSetType typei) { SetFEMesh_Set(typei); }
	FEMesh_Set(mystr namei) { SetFEMesh_Set(TSetUnknown, TArray<int>(0), namei); }                     // only name, used for inline TArray<T>::Find()
	void CopyFrom(const FEMesh_Set& other) 
	{ 
		type = other.type;
		if(other.name.Length())
			name = other.name;
		else
			name = mystr("");
		nodenrs.CopyFrom(other.nodenrs);
		elementnrs.CopyFrom(other.elementnrs);
		facenrs.CopyFrom(other.facenrs);
		bodynrs.CopyFrom(other.bodynrs);
		materialnrs.CopyFrom(other.materialnrs);
	}
	FEMesh_Set* GetCopy() { FEMesh_Set* s = new FEMesh_Set(); s->CopyFrom(*this); return s; }
	void SetFEMesh_Set(TSetType typei, TArray<int> listi = TArray<int>(0), mystr namei = mystr(""))
	{
		type = typei;
		FlushArrays(); // clear all arrays ??? 
		if(type==TSetNodes) nodenrs.CopyFrom(listi);
		if(type==TSetElements) elementnrs.CopyFrom(listi);
		if(type==TSetFaces) facenrs.CopyFrom(listi);
		if(type==TSetBodies) bodynrs.CopyFrom(listi);
		if(type==TSetMaterials) materialnrs.CopyFrom(listi);
		name = namei;
	}
	void FlushArrays() { nodenrs.Flush(); elementnrs.Flush(); facenrs.Flush(); bodynrs.Flush(); materialnrs.Flush(); }
  
public: // access I - individual access on members
	virtual int& Node(int i) {return nodenrs(i);}
	virtual int& Element(int i) {return elementnrs(i);}
	virtual int& Face(int i) {return facenrs(i);}
	virtual int& Body(int i) {return bodynrs(i);}
	virtual int& Material(int i) {return materialnrs(i);}
	virtual const TArray<int>& GetNodes() const { return nodenrs; }
	virtual const TArray<int>& GetElements() const { return elementnrs; }
	virtual const TArray<int>& GetFaces() const { return facenrs; }
	virtual const TArray<int>& GetBodies() const { return bodynrs; }
	virtual const TArray<int>& GetMaterials() const { return materialnrs; }
	virtual void SetNodes(TArray<int>& nodenrsi) { nodenrs = nodenrsi; }
	virtual void SetElements(TArray<int>& elementnrsi) { elementnrs = elementnrsi; }
	virtual void SetFaces(TArray<int>& facenrsi) { facenrs = facenrsi; }
	virtual void SetBodies(TArray<int>& bodynrsi) { bodynrs = bodynrsi; }
	virtual void SetMaterials(TArray<int>& materialnrsi) { materialnrs = materialnrsi; }
	virtual TArray<int>& Nodes() { return nodenrs;}
	virtual TArray<int>& Elements() { return elementnrs;}
	virtual TArray<int>& Faces() { return facenrs;}
	virtual TArray<int>& Bodies() { return bodynrs;}
	virtual TArray<int>& Materials() { return materialnrs;}
	virtual int NNodes() { return nodenrs.Length(); }
	virtual int NElements() { return elementnrs.Length(); }
	virtual int NFaces() { return facenrs.Length(); }
	virtual int NBodies() { return bodynrs.Length(); }
	virtual int NMaterials() { return materialnrs.Length(); }

	virtual const TSetType GetType() const { return type; }
	virtual void SetType(TSetType typei) { type = typei; }
	virtual TSetType& Type() { return type; }

	virtual const mystr GetName() const { return name; }
	virtual void SetName(mystr namei) { name = namei; }
	virtual mystr& Name() { return name; }

public: // access II - access with a selection flag (typei)
	virtual const TArray<int>& GetArray() const { return GetArray(GetType()); }
	virtual const TArray<int>& GetArray(TSetType typei) const 
	{
		switch (typei)
		{
		case TSetNodes: return GetNodes(); break;
		case TSetElements: return GetElements(); break;
		case TSetFaces: return GetFaces(); break;
		case TSetBodies: return GetBodies(); break;
		case TSetMaterials: return GetMaterials(); break;
		default: return *(new TArray<int>(0)); break;
		}
	}
	virtual void SetArray(TArray<int>& arrayi) { SetArray(arrayi, GetType()); }
  virtual void SetArray(TArray<int>& arrayi, TSetType typei)
	{
		switch (typei)
		{
		case TSetNodes: SetNodes(arrayi); break;
		case TSetElements: SetElements(arrayi); break;
		case TSetFaces: SetFaces(arrayi); break;
		case TSetBodies: SetBodies(arrayi); break;
		case TSetMaterials: SetMaterials(arrayi); break;
		default: break;
		}
	}
	virtual TArray<int>& Array() { return Array(GetType()); }
  virtual TArray<int>& Array(TSetType typei)
	{
		switch (typei)
		{
		case TSetNodes: return Nodes(); break;
		case TSetElements: return Elements(); break;
		case TSetFaces: return Faces(); break;
		case TSetBodies: return Bodies(); break;
		case TSetMaterials: return Materials(); break;
		default: return *(new TArray<int>(0)); break;
		}
	}
	virtual int& operator() (int i) { return Array(GetType())(i); }
	virtual int& operator() (int i, TSetType typei) { return Array(typei)(i); }

	virtual int N() { return Array(GetType()).Length(); }
	virtual int N(TSetType typei) { return Array(typei).Length(); }
	virtual int Length() {return N();}

public: // operations
  // add a single integer ( "Add", "Append", "<<" )
	int Add(int i) { return Array().Add(i); }
	int Append(int i) { return Array().Add(i); }
	int operator<< (int i) { return Array().Add(i); }
	int AddIfNotExists(int i) { return Array().AddIfNotExists(i); }

	int Add(int i, TSetType typei) { return Array(typei).Add(i); }
	int Append(int i, TSetType typei) { return Array(typei).Add(i); } 
	int AddIfNotExists(int i, TSetType typei) { return Array(typei).AddIfNotExists(i); }

	int SetNaturalIfEmpty(int length)
	{
		if (N() == 0)
		{
			Array() = NaturalNumbers(length);
		}
	}

	// copy/replace only "data" content ( keep name and type )
	friend FEMesh_Set& operator<< (FEMesh_Set& target, const FEMesh_Set& other) 
	{
		target.nodenrs = other.nodenrs;
		target.elementnrs = other.elementnrs;
		target.facenrs = other.facenrs;
		target.bodynrs = other.bodynrs;
		target.materialnrs = other.materialnrs;
		return target;
	}

	friend FEMesh_Set operator+ (const FEMesh_Set& a, const FEMesh_Set& b)  //set union, assignment - use c = a+b or c << a+b
	{
		FEMesh_Set c(a);
		c += b;
		return c;
	}
	void operator+= (const FEMesh_Set& other) //set union, self
	{
// keep type from this set
// keep name of this set
		nodenrs += other.nodenrs;
		elementnrs += other.elementnrs;
		facenrs += other.facenrs;
		bodynrs += other.bodynrs;
		materialnrs += other.materialnrs;
	}
	void operator+= (const TArray<int>& single) // set union for single array
	{
		Array()+= single;
	}

	friend FEMesh_Set operator& (const FEMesh_Set& a, const FEMesh_Set& b) //set intersection, assignment - use c = a&b or c << a&b
	{
		FEMesh_Set c(a);
		c &= b;
		return c;
	}
	void operator&= (const FEMesh_Set& other) //set intersection, self
	{
// keep type from this set
// keep name of this set
		nodenrs &= other.nodenrs;
		elementnrs &= other.elementnrs;
		facenrs &= other.facenrs;
		bodynrs &= other.bodynrs;
		materialnrs &= other.materialnrs;
	}

	friend FEMesh_Set operator- (const FEMesh_Set& a, const FEMesh_Set& b) //set difference, assignment - use c = a-b or c << a-b
	{
		FEMesh_Set c(a);
		c -= b;
		return c;
	}
	void operator-= (const FEMesh_Set& other) //set difference, self
	{
// keep type from this set
// keep name of this set
		nodenrs -= other.nodenrs;
		elementnrs -= other.elementnrs;
		facenrs -= other.facenrs;
		bodynrs -= other.bodynrs;
		materialnrs -= other.materialnrs;
	}

	bool operator== (const FEMesh_Set& other) // equality operator for TArray<T>::Find(const T& t)
	{
// to find in list the name has to match
		return name.Compare(other.name);
	}

private: // own variables
	TArray<int> nodenrs;
	TArray<int> elementnrs;
	TArray<int> facenrs;
	TArray<int> bodynrs;
	TArray<int> materialnrs;
  TSetType type;                 // defines the (main) type of the set entity - primary used array, others may be nonempty  
	mystr name;	                   // assign a name to the set 
};

//+ end FEMesh_Set
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  classes to create  cylindrical geometries from a 2d crossection (2d polygon)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class LineSegment
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains (consistent to netgen .in2d file format):  start point, end point, spline point, left&right domain, local refinement,     boundary contition number
class LineSegment : public GeomLine2D
{
public: // lifecycle
	LineSegment():GeomLine2D(),flag_spline(0),domain_left(0),domain_right(0),boundarycondition(0),refinement(0) { } // empty constructor required fpr TArray !!!
	LineSegment(LineSegment& other) { CopyFrom(other); }
	~LineSegment() { }
public: // lifecycle II
	void Reset() 
	{ 
		locpoints2D.Flush();
		P1() = 0;
		P2() = 0;
		P3() = 0;
		flag_spline = 0; 
		Domain1() = 0; 
		Domain2() = 0; 
		BC() = 0; 
		Ref() = 0; 
	}
	LineSegment(MBS* mbsi) { mbs = mbsi; }
	LineSegment(MBS* mbsi, int p1, int p2, int dom1, int dom2, int boundary, int ref) { SetLineSegment(p1,p2,dom1,dom2,boundary,ref); mbs = mbsi; }
	LineSegment(MBS* mbsi, int p1, int p2, int pspline, int dom1, int dom2, int boundary, int ref) { SetLineSegment(p1,p2,pspline,dom1,dom2,boundary,ref); mbs = mbsi; }
	LineSegment* GetCopy() { return (new LineSegment(*this)); }
	void CopyFrom(LineSegment& other)
	{
		GeomLine2D::CopyFrom(other);
		P1() = other.P1();
		P2() = other.P2();
		P3() = other.P3();
		flag_spline = other.IsSpline();
		Domain1() = other.Domain1();
		Domain2() = other.Domain2();
		BC() = other.BC();
		Ref() = other.Ref();
	}
	void SetLineSegment(int p1, int p2, int dom1, int dom2, int boundary, int ref) // 2 points - no default values possible for bc&ref due to ambiguities
	{	P1() = p1; P2() = p2; flag_spline = 0; Domain1() = dom1; Domain2() = dom2; BC() =  boundary; Ref() = ref; }
	void SetLineSegment(int p1, int p2, int pspline, int dom1, int dom2, int boundary, int ref) // 3 points - no default values possible for bc&ref due to ambiguities
	{	P1() = p1; P2() = p2; P3() = pspline; flag_spline = 1; Domain1() = dom1; Domain2() = dom2; BC() =  boundary; Ref() = ref; }

	void SetMBS(MBS* mbsi) { mbs = mbsi; }
public: // access
	int& P1() { return point1; }
	int& P2() { return point2; }
	int& P3() { return pointspline; }
	int IsSpline() { return flag_spline; }
	int& Domain1() { return domain_left; }         // "inner" domain
	int& Domain2() { return domain_right; }        // "outer" domain
	int& Ref() { return refinement; }
	int& BC() { return boundarycondition; }
public: // drawing routines
	void DrawYourself()
	{
		GeomLine2D::DrawYourself(); // !!! this routine is NEVER CALLED  !!!
	}
	void SetDrawPoints(Vector2D p1, Vector2D p2, Vector2D pspline = Vector2D(0.,0.))
	{
		locpoints2D.Flush();
		AddLocPoint2D(p1);
		AddLocPoint2D(p2);
		AddLocPoint2D(pspline);
	}
	void SetDrawColor(MyColorList* p_palette=NULL, int i=-1)
	{
    if(p_palette == NULL)
		{
			p_palette = new class MyColorList(); // default palet
		}
		if(i == -1) 
		{
			SetCol(p_palette->GetCol3D(BC())); // draw color depends on bondary condition number
		}
		else 
		{
			SetCol(p_palette->GetCol3D(i)); // draw color is chosen by call
		}
	}

public:
private: // own variables
// save coordinates in GeomLine::locpoints2D 
	int point1, point2, pointspline;           // global point numbers
	int flag_spline;
	int domain_left, domain_right;
	int boundarycondition;
	int refinement;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class ClosedPolygon
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: list of points (numbers), list of line sections(numbers)
// lists are stored in counterclockwise sequence, a closed polygon has only ONE domain 
class ClosedPolygon : public GeomPolygon2D
{
public: // lifecycle
	ClosedPolygon():GeomPolygon2D(),pointnumbers(0),linesegmentnumbers(0),pointsOuterEdge(0),pointsInnerEdge(0) { } // empty constructor required for TArray !!!
	ClosedPolygon(ClosedPolygon& other):GeomPolygon2D(other) { CopyFrom(other); }
	~ClosedPolygon() { }
public: // lifecycle II
	void Reset() 
	{ 
		PointNrs().Flush(); 
		LineNrs().Flush(); 
	}
	ClosedPolygon(MBS* mbsi):GeomPolygon2D() { mbs = mbsi; }
	ClosedPolygon(MBS* mbsi, TArray<int>& pointnumbers, TArray<int>& linesegmentnumbers):GeomPolygon2D() { SetClosedPolygon(pointnumbers,linesegmentnumbers); mbs = mbsi; }
	ClosedPolygon* GetCopy() { ClosedPolygon* pcp = new ClosedPolygon(*this); return pcp; }
	void CopyFrom(ClosedPolygon& other)	
	{	
		GeomElement::CopyFrom(other); 
		PointNrs().CopyFrom(other.PointNrs()); 
		LineNrs().CopyFrom(other.LineNrs());	
		OuterEdge().CopyFrom(other.OuterEdge());
		InnerEdge().CopyFrom(other.InnerEdge());
	}
	void SetClosedPolygon(TArray<int>& pointnumbersi, TArray<int>& linesegmentnumbersi) { PointNrs().CopyFrom(pointnumbersi); LineNrs().CopyFrom(linesegmentnumbersi); }
	void SetMBS(MBS* mbsi) { mbs = mbsi; }
public: // insert & remove Points / LineSegments

public: // access
	TArray<int>& PointNrs() { return pointnumbers; }
	int& PointNr(int i) { return PointNrs()(i); }
	int NPoints() { return PointNrs().Length(); }
	TArray<int>& LineNrs() { return linesegmentnumbers; }
	int& LineNr(int i) { return LineNrs()(i); }
	int NLines() { return LineNrs().Length(); }
	int& Domain() { return domainnumber; }
	TArray<int>& OuterEdge() { return pointsOuterEdge; }
	TArray<int>& InnerEdge() { return pointsInnerEdge; }
	int& PointNrOuterEdge(int i) { return OuterEdge()(i); }
	int& PointNrInnerEdge(int i) { return InnerEdge()(i); }


public: // base class overloaded
	virtual const char* GetElementSpec() const {return "ClosedPolygon";}
	virtual TGeomElement GetType() const {return TGeomPolygon2D;}

public: // drawing routines
	void DrawYourself() 
	{
		GeomPolygon2D::DrawYourself(); // !!! this routine is NEVER CALLED  !!!
	}
	int operator<< (int n) { return Rotate(n); }
	int operator>> (int n) { return Rotate(-n); }
	int Rotate(int n); // moves the first point of the polygon n steps in counter-clockwise direction
	
	void SetDrawPoints(TArray<Vector2D>& loccoordsi)
	{
		locpoints2D.Flush();
		for(int i=1; i<= loccoordsi.Length(); i++)
		{
			AddLocPoint2D(loccoordsi(i));
		}
	}
	void SetDrawColor(MyColorList* p_palette=NULL, int i=-1)
	{
    if(p_palette == NULL)
		{
			p_palette = new class MyColorList(); // default palette
		}
		Vector3D color; // help debug
		if(i == -1) 
		{
			color = p_palette->GetCol3D(Domain()); // draw color depends on domain number
		}
		else 
		{
			color = p_palette->GetCol3D(i); // draw color is set in function call
		}
		SetCol(color);
	}

private: // own variables
	TArray<int> pointnumbers;          // global point numbers
	TArray<int> linesegmentnumbers;    // global linesegment numbers
	int domainnumber;                  // domainnumber of the enclosed area

	TArray<int> pointsOuterEdge;			// global point numbers of outer edge with increasing x-coordinate
																		// it is possible that pointsOuterEdge(1)= pointsInnerEdge(1)
	TArray<int> pointsInnerEdge;			// global point numbers of outer edge with increasing x-coordinate	
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class Crosssection - convert a 2D crossetion defined by polygons into a quad-mesh and cylindrical parts
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: list of points, list of line sections, list of closed polygons
//
// TODO: routines for: input from file, draw, generate 4-gons,  
class Crosssection
{
public: // lifecycle
	Crosssection() { }  // empty constructor required fpr TArray !!!
	Crosssection(Crosssection& other) { CopyFrom(other); }
	~Crosssection() { }
public: // lifecycle II 
	void Reset() { Points().Flush(); Lines().Flush(); Polygons().Flush(); Boundaries().Reset();}
	Crosssection(MBS* mbsi) { mbs = mbsi; } 
	void CopyFrom(Crosssection& other)
	{	
		mbs = other.mbs;
		Points().CopyFrom(other.Points());
		Lines().CopyFrom(other.Lines());
		Polygons().CopyFrom(other.Polygons());
		Boundaries().CopyFrom(other.Boundaries());	//$ DR 2012-08-16
	}

public: // file import 
	int ImportFromIN2D(mystr& filename, int read_boundaries=0); 
	int ReadBoundaryConditions(mystr& filename); //$ DR 2012-08-16

public: // automatic polygon identification
	int FindDomainPolygons(); // find closed polygons, each containing only a single domain 
	int HighestDomainNumber(); // search highest domainnumber in all lines			//$ DR 2012-08-20 moved from protected to public
protected:
	int FilterLineSegmentsForDomainNumber(int domainnumber, TArray<int>& domain_on_left_side, TArray<int>& domain_on_right_side, TArray<int>& domain_on_both_sides = TArray<int>(0)); // finds all LineSegments that border a domain, use these to find closed polygons
  int FindContinuation(int lastpoint, TArray<int>& candidates, int flag_leftright); // Find a continuation for the open polygon, returns number of next line
  int ResortPolygonCCW(int i); // make sure that the polygon encloses the domain on the left ( counter clock wise orientation ) - ATTETNION: DO THIS AT AN EARLY STAGE, WILL NOT WORK FOR FINISHED MESHING

public:	
	void GetPolygonOuterAndInnerEdge();

public: // Splitting I: Add Points
	int AddPointOnLine(Vector2D splitpoint, int linesegmentnumber = 0); // Additional Point - splits Linesegment in two - updates Polygons
protected:
	int FindLineSegment(Vector2D& splitpoint); // Find the Segment the point is on
	int IsOnLineSegment(Vector2D& splitpoint, int linenr); // checks if the point is on the given line segment
  int SplitLineSegment(int linesegmentnumber, int splitpointnumber); // splits one linesegment into two at the splitpoint
	int StitchPolygon(int polygonnr, int line_index, int newlinenumber, int splitpointnumber); // closes hole in polygon

public: // turn splines into several linear segments
 	int SplinesToLinear(int segments_manual = -1); // split the segments into linear sections, if no segmentation is defined, use the "refinement" entry
	int SplitSpLine(int linesegment_number, int segments); // split the line into equally long segments, checks if the line segment is a spline and moves the points accordingly

protected:
	virtual Vector2D InterpolateLinear(double factor, Vector2D& p1, Vector2D& p2);
	virtual Vector2D InterpolateQuadraticBezier(double factor, Vector2D& p1, Vector2D& p2, Vector2D& ps);


public: // Splitting II: bud a new polygon
	int BudPolygon(int p1_nr, int p2_nr); // bud a polygon where points are known
	//////int BudPolygon(int polygonnumber, int linenr); // bud a polygon where linesegments are known

public: // Splitting III: Cutting along a specified vector
	int CutAligned(Vector2D guide, double angular_variance_deg = 1e-3); // cut all polygons along the guide vector with allowed angular variance
	int CutAligned_ContinueLines(Vector2D guide, double angular_variance_deg = 1e-3); // try to continue existing lines
	int CutAligned_AtEachPoint(Vector2D guide, double angular_variance_deg = 1e-3); // try to cut at each point
protected:
	int IsInnerLineOfPolygon(int polygonnr, int startingpointnr, int targetlinenr, Vector2D endpoint); 
  int GetBudlineEndPointCoords(int startpointnr, int targetlinenr, Vector2D guide, Vector2D& rv_endpoint_coords, double angular_variance_rad = 1e-3 / 180. * MY_PI);
	int LineSegmentExistsInPolygon(int startpoint, int endpoint, int polygonnumber);
protected: 
	double TwoPointVectorAngleToReference(int p1_idx, int p2_idx, Vector2D reference = Vector2D(0.,1.)); // returns the angle in rad of a vector specified by 2 points to the specified reference
	double LineSegmentAngleToReference(int linesegmentnumber, Vector2D reference = Vector2D(0.,1.)); // returns the angle in rad of the linesegment to the specified reference
  Vector2D CutStraightLineWithLineSegment_LambdaMu(int pointnumber, Vector2D guide, int linesegmentnumber); // returns parameters lambda and mu of the intersection point on the linesegment
	Vector2D CutTwoLineSegments_LambdaMu(int linesegmentnumber1, int linesegmentnumber2); // returns the parameters lambda and mu of the intersection point on the second linesegment
protected: 
	int CommonPointOfLineSegments(int linenr1, int linenr2);
  int CommonPolygonOfLineSegments(int linenr1, int linenr2);
  int CommonPolygonOfPointAndLineSegment(int pointnr, int linenr);

public: // Merging


public: // Drawing preparation
	int AssingDrawPointsAndColorToPolygon(int polygonnumber); // assign values for MBS::Draw

public: // export
	int ExportEdgesForNetgen(TArray<Vector3D>& leftpoints, TArray<Vector3D>& rightpoints, TArray<int>& materials, TArray<int>& refinements, int draw_edges=0); // export the edges to create the 3D rotatory geometry in Netgen
	int SimplifyConesToCylinders(TArray<Vector3D>& leftpoints, TArray<Vector3D>& rightpoints); //$ DR 2012-10-22: simplification of the geometry: (cutted) cones are converted to (cutted) cylinders
	int IsOuterEdge(int polygonnumber, int linesegmentnumber); // determines if the linesegment is an "outer" edge of the polygon
	int IsInnerEdge(int polygonnumber, int linesegmentnumber); // determines if the linesegment is an "inner" edge of the polygon - line at y=0 (rotation axis) does not count as inner edge
  

public: // access
	TArray<Vector2D>& Points() { return points; }
	Vector2D& Point(int i) { return Points()(i); }
	int NPoints() { return Points().Length(); }
	TArrayDynamic<LineSegment>& Lines() { return linesegments; }
	LineSegment& Line(int i) { return Lines()(i); }
	int NLines() { return Lines().Length(); }
	TArrayDynamic<ClosedPolygon>& Polygons() { return polygons; }
	ClosedPolygon& Polygon(int i) { return Polygons()(i); }
	int NPolygons() { return Polygons().Length(); }
	MBS* GetMBS() { return mbs; }
	ElementDataContainer& Boundaries() {return boundaries;}	//$ DR 2012-08-16

private: // own variables
	TArray<Vector2D> points;                   // point coordinates
	TArrayDynamic<LineSegment> linesegments;
	TArrayDynamic<ClosedPolygon> polygons;
	MBS* mbs;
	class MyColorList mycolors;
	ElementDataContainer boundaries; //$ DR 2012-08-16
};

// enum for directions (of blocks)
typedef enum { Txmin = 1, Txmax = 2, Tymin = 3, Tymax = 4, Tzmin = 5, Tzmax = 6 } THexDirections;
typedef enum { Tinnner = 1, Touter = 2, Tbottom = 5, Ttop = 6  } TCylDirections;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class SegmentedBlock: 3D segmented hexahedral or cylinder (cutting planes)
//+                               2D segmented quadrilateral 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodes, elements, cutting-planes, divisions (of segments), occupation (of segments)
//
// a BLOCK is a volume defined by geometry parameters
// the SEGMENTEDBLOCK defines fixed cutting planes and occupation for all segments
// the BlockFaces sum over the Faces of the corresponding segments
//
// TODO: fix stand-alone generation...
class SegmentedBlock
{
public: // lifecycle
	SegmentedBlock():blocknodes(TSetNodes),blockelements(TSetElements) { }
	SegmentedBlock(SegmentedBlock& other) { CopyFrom(other); }
	~SegmentedBlock() { }
public: // lifecycle II
	void Reset() // set defaults
	{ 
		QMin() = Vector3D(0.); 
		QMax() = Vector3D(1.); 
		Thickness() = 0.;
		BlockDivs() = int4(1,1,1,0);
		MatNr() = 1; 
		Domain() = 1;
		Color() = Vector3D(-1.,0.,0.); // default Color -RED will be replaced by material color !!
		CutQ1().Flush(); 
		CutQ2().Flush(); 
		CutQ3().Flush(); 
		OccuArray().Flush(); 
		DivsArray().Flush();
		BlockNodes().FlushArrays();
		BlockElements().FlushArrays(); 
	}
	SegmentedBlock(Vector3D q_min, Vector3D q_max, int matnr, int domain) { SetSegmentedBlock( q_min, q_max, matnr, domain); }
	SegmentedBlock* GetCopy() { return (new SegmentedBlock(*this)); }
	void CopyFrom(SegmentedBlock& other)
	{
	  QMin() = other.QMin();
		QMax() = other.QMax();
		Thickness() = other.Thickness();
		BlockDivs() = other.BlockDivs();
		MatNr() = other.MatNr();
		Domain() = other.Domain();
		Color() = other.Color();

		CutQ1().CopyFrom(other.CutQ1());
		CutQ2().CopyFrom(other.CutQ2());
		CutQ3().CopyFrom(other.CutQ3());
		OccuArray().CopyFrom(other.OccuArray());
		DivsArray().CopyFrom(other.DivsArray());
	
		BlockNodes().CopyFrom(other.BlockNodes());
		BlockElements().CopyFrom(other.BlockElements());
		ResetBlockFaces(other.BlockFaces().Length());
		for(int i=1; i<= BlockFaces().Length(); i++)
		{
			BlockFace(i).CopyFrom(other.BlockFace(i));
		}
	}
	void SetSegmentedBlock( Vector3D q_mini, Vector3D q_maxi, int matnri, int domaini, Vector3D color = Vector3D(1.,0.,0.) )
	{
		QMin() = q_mini;
		QMax() = q_maxi;
	  CutQ1().Add(QMin().X());
		CutQ1().Add(QMax().X());
	  CutQ2().Add(QMin().Y());
		CutQ2().Add(QMax().Y());
	  CutQ3().Add(QMin().Z());
		CutQ3().Add(QMax().Z());
		MatNr() = matnri;
		Domain() = domaini;
		Color() = color;
	}
	void ResetBlockFaces(int len = 0) 
	{ 
		if (!len) len = BlockFaces().Length();
		BlockFaces().Flush();
		BlockFaces().SetLen(len);
		BlockFaces().SetAll(FEMesh_Set(TSetFaces));
	}
public: // segments
	void SetAllOccupations(int i)
	{
		OccuArray().Flush();
		OccuArray().SetLen((CutQ1().Length()-1)*(CutQ2().Length()-1)*(CutQ3().Length()-1));
		OccuArray().SetAll(i);
	}
	void SetAllDivisions(int4 i)
	{
		DivsArray().Flush();
		DivsArray().SetLen((CutQ1().Length()-1)*(CutQ2().Length()-1)*(CutQ3().Length()-1));
		DivsArray().SetAll(i);
	}
public: // access
	Vector3D& QMin() { return q_min; }              // additional in derived class
	Vector3D& QMax() { return q_max; }							// additional in derived class
	double& Thickness() { return thickness; }
	TArray<int4>& DivsArray() { return divisions; } 
	int4& BlockDivs() { return blockdivisions; }
	int& MatNr() { return matnr; }
	int& Domain() { return domain; }
	Vector3D& Color() { return color; }

	TArray<double>& CutQ1() { return cut_q1; }      // additional in derived class
	TArray<double>& CutQ2() { return cut_q2; }      // additional in derived class
	TArray<double>& CutQ3() { return cut_q3; }      // additional in derived class
	TArray<int>& OccuArray() { return occuptaion; }
	FEMesh_Set& BlockNodes() { return blocknodes; }
	FEMesh_Set& BlockElements() { return blockelements; }
	FEMesh_Set& BlockFace(int i) { return blockfaces(i); }
	TArrayDynamic<FEMesh_Set>& BlockFaces() { return blockfaces; }
private: // own variables
// geometry 
	Vector3D q_min, q_max;                 // set of coordinates: [x,y,z] or [r,phi,h] --> have own access functions in derived class. DEFAULT [0,0,0 - 1,1,1]
	double thickness;                      // for 2D: element thickness in 3rd dimension. DEFAULT 0. - mesh.setting will be used
// mesh
	int4 blockdivisions;                   // mesh divisions according to coordinates  --> have own access functions in derived class
// elements
	int matnr;                             // material number (defined in mesh)
	int domain;                            // domain number == body number
  Vector3D color;                        // assigned color for the elements

// segmentation 
  TArray<double> cut_q1, cut_q2, cut_q3; // cutting "planes"                         --> have own access functions in derived class
	TArray<int> occuptaion;                // list whether the segments are created 
	TArray<int4> divisions;                // list containing all divisions of the segment
// statistics
	FEMesh_Set blocknodes;
	FEMesh_Set blockelements;
	TArrayDynamic<FEMesh_Set> blockfaces;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  segmented hexahedral
class SegmentedHexBlock: public SegmentedBlock
{
public: // lifecycle
	SegmentedHexBlock() { Reset(); }
	SegmentedHexBlock(SegmentedHexBlock& other) { CopyFrom(other); }
	~SegmentedHexBlock() { }
public: // lifecycle II
	void Reset() { SegmentedBlock::Reset(); SegmentedBlock::ResetBlockFaces(6); }
	SegmentedHexBlock(Box3D box, int matnr, int domain) { Reset(); SetSegmentedHexBlock(box, matnr, domain); }
	SegmentedHexBlock* GetCopy() { return (new SegmentedHexBlock(*this)); }
	void CopyFrom(SegmentedHexBlock& other)	{	SegmentedBlock::CopyFrom(other); }
	void SetSegmentedHexBlock(Box3D boxi, int matnri, int domaini) { SetSegmentedBlock(boxi.PMin(), boxi.PMax(), matnri, domaini); }
public: // access
	Vector3D& Corner1() { return QMin(); }
	Vector3D& Corner2() { return QMax(); }
	TArray<double>& XCut() { return CutQ1(); }
	TArray<double>& YCut() { return CutQ2(); }
	TArray<double>& ZCut() { return CutQ3(); }
public: // additional access
	Box3D GetBox() { return Box3D(QMin(),QMax()); }
	double& XMin() { return QMin().X(); }
	double& XMax() { return QMax().X(); }
	double& YMin() { return QMin().Y(); }
	double& YMax() { return QMax().Y(); }
	double& ZMin() { return QMin().Z(); }
	double& ZMax() { return QMax().Z(); }
	double& XCut(int i) { return XCut()(i); }
	double& YCut(int i) { return YCut()(i); }
	double& ZCut(int i) { return ZCut()(i); }
	int& XDivs(int i) { return DivsArray()(i)(1); }
	int& YDivs(int i) { return DivsArray()(i)(2); }
	int& ZDivs(int i) { return DivsArray()(i)(3); }
	int& XDivs_Block() { return BlockDivs()(1); }
	int& YDivs_Block() { return BlockDivs()(2); }
	int& ZDivs_Block() { return BlockDivs()(3); }
public: // segments
	int& Occupation(int x, int y, int z) { return OccuArray()(IofXYZ(x,y,z)); }
	int4& Division(int x, int y, int z) { return DivsArray()(IofXYZ(x,y,z)); }
	int IofXYZ(int x, int y, int z) // maps x,y,z counts to index count(array populate)
	{
		int divx = XCut().Length()-1; // number of intervals = number of cutting planes plus 1
		int divy = YCut().Length()-1;
		int divz = ZCut().Length()-1;
	  int i = 1 + (x-1) + (y-1)*divx + (z-1)*divx*divy; // compute linear index
		return i;
	}
public: // Read/Generate
	void ReadFromEDC(ElementDataContainer* edc);  // read a single hex block from edc !defined format!
	void GenerateIn(class FEMesh* p_mesh, Vector3D meshsize); // creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
  void RefreshLists(class FEMesh* p_mesh, int flag_faces = 0); // recomputes the lists of elements, nodes, etc for the block

private: // own variables
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  segmented cylinder 
class SegmentedCylBlock: public SegmentedBlock
{
public: // lifecycle
	SegmentedCylBlock() { Reset(); }
	SegmentedCylBlock(SegmentedCylBlock& other) { CopyFrom(other); }
	~SegmentedCylBlock() {}
public: // lifecycle II
	void Reset() { SegmentedBlock::Reset();  SegmentedBlock::ResetBlockFaces(6); AngRef() = 0; }
	SegmentedCylBlock(double ri, double ro, double hb, double ht, int matnr, int domain) { Reset(); SetSegmentedCylBlock(ri, ro, hb, ht, matnr, domain); }
	SegmentedCylBlock* GetCopy() { return (new SegmentedCylBlock(*this)); }
	void CopyFrom(SegmentedCylBlock& other) { SegmentedBlock::CopyFrom(other); AngRef() = other.AngRef();	}
	void SetSegmentedCylBlock(double rii, double roi, double hbi, double hti, int matnri, int domaini) { SetSegmentedBlock(Vector3D(rii,0.,hbi), Vector3D(roi,0.,hti), matnri, domaini);	}
public: // access
	TArray<double>& RadCut() { return CutQ1(); }
	TArray<double>& AxiCut() { return CutQ3(); }
public: // additional access
	double& RI() { return QMin().X(); }
	double& RO() { return QMax().X(); }
	double& HB() { return QMin().Z(); }
	double& HT() { return QMax().Z(); }
	double& RadCut(int i) { return RadCut()(i); }
	double& AxiCut(int i) { return AxiCut()(i); }
	int& RadDivs(int i) { return DivsArray()(i)(1); }
	int& TanDivs(int i) { return DivsArray()(i)(2); }
	int& AxiDivs(int i) { return DivsArray()(i)(3); }
	int& Refine(int i) { return DivsArray()(i)(4); }
	int& RadDivs_Block() { return BlockDivs()(1); }
	int& TanDivs_Block() { return BlockDivs()(2); }
	int& AxiDivs_Block() { return BlockDivs()(3); }
	int& AngRef() { return BlockDivs()(4); }
public: // segments
	int& Occupation(int r, int h) {	return OccuArray()(IofRH(r,h));	}
	int4& Division(int r, int h) { return DivsArray()(IofRH(r,h)); }
	int& RadDivs(int r, int h) { return DivsArray()(IofRH(r,h))(1); }
	int& TanDivs(int r, int h) { return DivsArray()(IofRH(r,h))(2); }
	int& AxiDivs(int r, int h) { return DivsArray()(IofRH(r,h))(3); }
	int& Refine(int r, int h) { return DivsArray()(IofRH(r,h))(4); }
	int IofRH(int r, int h) // maps r,h counts to index count(array populate)
	{
		int divr = RadCut().Length()-1; // number of intervals = number of cutting planes plus 1
		int divh = AxiCut().Length()-1;
	  int i = 1 + (r-1) + (h-1)*divr; // compute linear index
		return i;
	}
public: // Read/Generate
	void ReadFromEDC(ElementDataContainer* edc);  // read a single cyl block from edc !defined format!
	void GenerateIn(class FEMesh* themesh, Vector3D meshsize); // creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
	void RefreshLists(class FEMesh* p_mesh, Vector3D& base, Vector3D& axis, int flag_faces=0); 

private: // own variables
	int angref;                                   // flag for angular refinement, # of refinement steps
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  segmented quadrilateral
class SegmentedQuadBlock: public SegmentedHexBlock
{
public: // lifecycle
	SegmentedQuadBlock() { Reset(); }
	SegmentedQuadBlock(SegmentedQuadBlock& other) { CopyFrom(other); }
	~SegmentedQuadBlock() { }
public: // lifecycle II
	void Reset() { SegmentedBlock::Reset(); SegmentedBlock::ResetBlockFaces(4); ZMin() = -MESH_STD_TOL; ZMax() = MESH_STD_TOL; }
	SegmentedQuadBlock(Box3D box, int matnr, int domain) { Reset(); SetSegmentedQuadBlock(box, matnr, domain); }
	SegmentedQuadBlock* GetCopy() { return (new SegmentedQuadBlock(*this)); }
	void CopyFrom(SegmentedQuadBlock& other)	{	SegmentedBlock::CopyFrom(other); }
	void SetSegmentedQuadBlock(Box3D boxi, int matnri, int domaini) { SetSegmentedBlock(boxi.PMin(), boxi.PMax(), matnri, domaini); }
public: // access - from SegmentedHexBlock, ATTENTION Corners() are Vector3D
public: // additional access - from SegmentedHexBlock, ATTENTION Corners() are Vector3D
public: // segments
	int& Occupation(int x, int y) { return OccuArray()(IofXY(x,y)); }
	int4& Division(int x, int y) { return DivsArray()(IofXY(x,y)); }
	int IofXY(int x, int y) // maps x,y,z counts to index count(array populate)
	{
		return IofXYZ(x,y,1);
	}
public: // Read/Generate
	void ReadFromEDC(ElementDataContainer* edc);  // read a single hex block from edc !defined format!
	void GenerateIn(class FEMesh* p_mesh, Vector3D meshsize); // creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
  void RefreshLists(class FEMesh* p_mesh, int flag_faces = 0); // recomputes the lists of elements, nodes, etc for the block
private: // own variables
};

//+  end SegmentedBlock and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class MeshedPart: 3D several hexahedral or cylinder (SegmentedBlocks)with consistent mesh
//+                           2D several quadilateral (SegmentedBlocks) with consistent mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodes, elements, blocks, cutting-planes, divisions (of segments), occupation (of segments),
//
// a MESHEDPART contains a list of blocks parameters common to all blocks
// functions of MESHEDPART guarante a consistent mesh for all blocks
class MeshedPart
{
public: // lifecycle
	MeshedPart() { }
	MeshedPart(MeshedPart& other) { CopyFrom(other); }
	~MeshedPart() { }
public: // lifecycle II
  void Reset()
	{
	  ResetBlocks();
		MeshWidth() = Vector3D(1.,1.,1.);
		CutQ1().Flush();
		CutQ2().Flush();
		CutQ3().Flush();
		
		PartElements().FlushArrays();
		PartNodes().FlushArrays();
		Name() = "";

		Translation() = Vector3D(0.,0.,0.);
		Rotation().SetDiag(1.);
	}
	MeshedPart(FEMesh* p_mesh) { MeshPtr() = p_mesh; }
	MeshedPart* GetCopy() { return (new MeshedPart(*this)); }
	void CopyFrom(MeshedPart& other)
	{
		ResetBlocks(other.Blocks().Length());
		for(int i=1; i<=other.Blocks().Length(); i++)
		{
			this->_Block(i).CopyFrom(other._Block(i));
		}
		MeshWidth() = other.MeshWidth();
		CutQ1().CopyFrom(other.CutQ1());
		CutQ2().CopyFrom(other.CutQ2());
		CutQ3().CopyFrom(other.CutQ3());
		MeshPtr() = (class FEMesh*) other.MeshPtr();
		PartElements().CopyFrom(other.PartElements());
		PartNodes().CopyFrom(other.PartNodes());
		ResetPartFaces(other.PartFaces().Length());
		for(int i=1; i<=other.PartFaces().Length(); i++)
		{
			this->PartFace(i).CopyFrom(other.PartFace(i));
		}
		Name() = other.Name();
		Translation() = other.Translation();
		Rotation() = other.Rotation();
	}
	void SetMeshedPart(TArray<SegmentedBlock*> blocksi) 
	{ 
		ResetBlocks(blocksi.Length());
		for(int i=1; i<=blocksi.Length(); i++)
		{
			_Block(i).CopyFrom(*blocksi(i));
		}
	}
	void ResetBlocks(int len = 0)
	{
		ReleaseArray_TemplatePtr(blocks); 
		Blocks().SetLen(len);	
		for(int i=1; i <= len; i++) 
			Blocks()(i) = new SegmentedBlock; 
	}
	void ResetPartFaces(int len = 0) 
	{ 
		if (!len) len = PartFaces().Length();
		PartFaces().Flush();
		PartFaces().SetLen(len);	
		PartFaces().SetAll(FEMesh_Set(TSetFaces));
	}
public: // functions I - override in derived classes
	virtual void ReadFromEDC(ElementDataContainer* edc) {} // read a sub-edc (part) containing several hexblocks !defined format!
	virtual void Cut() {}                                  // cuts the blocks in the list - prepares consistend meshing of the part
	virtual void FindCuttingPlanes() {}                    // computes all cutting planes (sorted), arrays include overall limits
	virtual void ComputeOccupationArray() {}               // computes which block occupies the volume-segment (the latter block always gets priority)
	virtual void ComputeDivisionsArray() {}                // computes the divisions for the blocks
	virtual void SetBlockData() {}                         // writes computed sectioning to the blocks
	virtual void Generate() {}                             // creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
	virtual void RefreshLists(int flag_faces = 0) {}       // recomputes the lists of elements, nodes, etc for the entire part

public: // functions II
  Box3D GetSurroundingBox()
	{
		if(! Blocks().Length()) return Box3D();
		Box3D overall(GetBox(1));
		for(int i=2; i<=Blocks().Length(); i++)
		{
			overall.Add(_Block(i).QMin());
			overall.Add(_Block(i).QMax());
		}
		return overall;
	}
	Box3D GetBox(int i) { return Box3D(_Block(i).QMin(), _Block(i).QMax()); }
public: // ... 
	void AddCuttingPlane(int i, double val)         // adds a cutting plane for the part (global cutting plane)
	{
		if(i==1) CutQ1().Add(val);
		else if(i==2) CutQ2().Add(val);
		else if(i==3) CutQ3().Add(val);
	}
	int AddBlock(SegmentedBlock& blocki)
	{
		return Blocks().Add(blocki.GetCopy());
	}
public: // access
	TArray<SegmentedBlock*>& Blocks() { return blocks; }
	SegmentedBlock& _Block(int i) { return *(blocks(i)); } // type cast version in derived classes
	Vector3D& MeshWidth() { return mesh_width; }         // additional in derived class
	double& MeshWidth(int i) { return mesh_width(i); }   // additional in derived class
	TArray<double>& CutQ1() { return cut_q1; }      // additional in derived class
	TArray<double>& CutQ2() { return cut_q2; }      // additional in derived class
	TArray<double>& CutQ3() { return cut_q3; }      // additional in derived class
	TArray<int>& OccuArray() { return occupation; }
	int& Occupation(int i) { return OccuArray()(i); }
	TArray<int4>& DivsArray() { return divisions; }
	int4& Division(int i) { return DivsArray()(i); }
	(class FEMesh*)& MeshPtr() { return p_mesh; }
	FEMesh_Set& PartElements() { return elements; }
	FEMesh_Set& PartNodes() { return nodes; }
	FEMesh_Set& PartFace(int i) { return partfaces(i); }
	TArray<FEMesh_Set>& PartFaces() { return partfaces; }
	mystr& Name() { return name; }
	Vector3D& Base() { return base; }
	Vector3D& Axis() { return axis; }
	Vector3D& Translation() { return translation; }
	Matrix3D& Rotation() { return rotation; }
protected: // access, restricted

private: // own variables
  TArray<SegmentedBlock*> blocks;                 // array of Blocks
	Vector3D mesh_width;                            // desired size for the elements [x,y,z] or [radial, tangential, axial]
  TArray<double> cut_q1, cut_q2, cut_q3;          // cutting "planes"            
	TArray<int> occupation;                         // occupation array for the entire part
	TArray<int4> divisions;                         // divisions [x,y,z, ] or [radial, tangential, axial, refinement]
	class FEMesh* p_mesh;                           
	FEMesh_Set elements;
	FEMesh_Set nodes;
	TArrayDynamic<FEMesh_Set> partfaces;            // TODO: outer faces for the parts...
	mystr name;
	Vector3D base;
	Vector3D axis;
	Vector3D translation;
	Matrix3D rotation;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  several hexahedral 
class MeshedHexPart: public MeshedPart
{
public: // lifecycle
	MeshedHexPart() { Reset(); }
	MeshedHexPart(MeshedHexPart& other) { CopyFrom(other); }
	~MeshedHexPart() { Reset(); }
public: // lifecycle II
	void Reset() { MeshedPart::Reset(); }
	MeshedHexPart(FEMesh* p_mesh) { MeshPtr() = p_mesh; }
	MeshedHexPart(TArray<SegmentedHexBlock*>& blocks, FEMesh* p_mesh = NULL ) { SetMeshedHexPart(blocks); MeshPtr() = p_mesh; }
	MeshedHexPart* GetCopy() { return (new MeshedHexPart(*this)); }
	void CopyFrom(MeshedHexPart& other) { MeshedPart::CopyFrom(other); }
	void SetMeshedHexPart(TArray<SegmentedHexBlock*>& blocksi) // why does base class function not work ?
	{
		ResetBlocks(blocksi.Length());
		for(int i=1; i<=blocksi.Length(); i++)
		{
			Blocks()(i) = blocksi(i)->GetCopy();
		}
	}
public: // functions
  void ReadFromEDC(ElementDataContainer* edc);  // read a sub-edc (part) containing several hexblocks !defined format!
	void Cut();                                   // cuts the blocks in the list - prepares consistend meshing of the part
	void FindCuttingPlanes();                     // computes all cutting planes (sorted), arrays include overall limits
	void ComputeOccupationArray();                // computes which block occupies the volume-segment (the latter block always gets priority)
  void ComputeDivisionsArray();                 // computes the divisions for the blocks
	void SetBlockData();                          // writes computed sectioning to the blocks
	void Generate();                              // creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
	void RefreshLists(int flag_faces = 0);        // recomputes the lists of elements, nodes, etc for the entire part
public: // access  
	SegmentedHexBlock& Block(int i) { return (SegmentedHexBlock&) _Block(i); }
	SegmentedHexBlock& operator() (int i) { return (SegmentedHexBlock&) _Block(i); }
	double& MeshX() { return MeshWidth(1); }
	double& MeshY() { return MeshWidth(2); }
	double& MeshZ() { return MeshWidth(3); }
	TArray<double>& XCut() { return CutQ1(); }
	TArray<double>& YCut() { return CutQ2(); }
	TArray<double>& ZCut() { return CutQ3(); }
	double& XCut(int i) { return XCut()(i); }
	double& YCut(int i) { return YCut()(i); }
	double& ZCut(int i) { return ZCut()(i); }
private: // own variables
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  several cylinder
class MeshedCylPart: public MeshedPart
{
public: // lifecycle
	MeshedCylPart() { Reset(); }
	MeshedCylPart(MeshedCylPart& other) { CopyFrom(other); }
	~MeshedCylPart() { Reset(); }
public: // lifecycle II
	void Reset() { MeshedPart::Reset(); Quadrants() = 1; }
	MeshedCylPart(FEMesh* p_mesh) { MeshPtr()  = p_mesh; }
	MeshedCylPart(TArray<SegmentedCylBlock*>& blocks, FEMesh* p_mesh = NULL) { SetMeshedCylPart(blocks); MeshPtr() = p_mesh; }
	MeshedCylPart* GetCopy() { return (new MeshedCylPart(*this)); }
	void CopyFrom(MeshedCylPart& other) { MeshedPart::CopyFrom(other); Quadrants() = other.Quadrants(); }
	void SetMeshedCylPart(TArray<SegmentedCylBlock*>& blocksi) // why does base class function not work ?
	{
		ResetBlocks(blocksi.Length());
		for(int i=1; i<=blocksi.Length(); i++)
		{
			Blocks()(i) = blocksi(i)->GetCopy();
		}
	}
public: // functions
	void ReadFromEDC(ElementDataContainer* edc);  // read a sub-edc (part) containing several cylblocks !defined format!
	void Cut();                                   // cuts the blocks in the list - prepares consistend meshing of the part
	void FindCuttingPlanes();                     // computes all cutting planes (sorted), arrays include overall limits
	void ComputeOccupationArray();                // computes which block occupies the volume-segment (the latter block always gets priority)
  void ComputeDivisionsArray();                 // computes the divisions for the blocks
	void SetBlockData();                          // writes computed sectioning to the blocks
	void Generate();                              // creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
	void RefreshLists(int flag_faces = 0);        // recomputes the lists of elements, nodes, etc for the entire part
public: // access  
	SegmentedCylBlock& Block(int i) { return (SegmentedCylBlock&) _Block(i); }
	SegmentedCylBlock& operator() (int i) { return (SegmentedCylBlock&) _Block(i); }
	double& MeshRad() { return MeshWidth(1); }
	double& MeshTan() { return MeshWidth(2); }
	double& MeshAxi() { return MeshWidth(3); }
	TArray<double>& RadCut() { return CutQ1(); }
	TArray<double>& AxiCut() { return CutQ3(); }
	double& RadCut(int i) { return RadCut()(i); }
	double& AxiCut(int i) { return AxiCut()(i); }
	int& Quadrants() { return quadrants; }

protected: 
	TArray<int>& AngSeg() { return anglesegments; }
	int& AngSeg(int i) { return AngSeg()(i); }
	TArray<int> anglesegments; // list of number of angular elements (45°) for each radial cutting plane

private: // own variables
	int quadrants;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  several quads 
class MeshedQuadPart: public MeshedHexPart
{
public: // lifecycle
	MeshedQuadPart() { Reset(); }
	MeshedQuadPart(MeshedQuadPart& other) { CopyFrom(other); }
	~MeshedQuadPart() { Reset(); }
public: // lifecycle II
	void Reset() { MeshedPart::Reset(); }
	MeshedQuadPart(FEMesh* p_mesh) { MeshPtr() = p_mesh; }
	MeshedQuadPart(TArray<SegmentedQuadBlock*>& blocks, FEMesh* p_mesh = NULL ) { SetMeshedQuadPart(blocks); MeshPtr() = p_mesh; }
	MeshedQuadPart* GetCopy() { return (new MeshedQuadPart(*this)); }
	void CopyFrom(MeshedQuadPart& other) { MeshedPart::CopyFrom(other); }
	void SetMeshedQuadPart(TArray<SegmentedQuadBlock*>& blocksi) // why does base class function not work ?
	{
		ResetBlocks(blocksi.Length());
		for(int i=1; i<=blocksi.Length(); i++)
		{
			Blocks()(i) = blocksi(i)->GetCopy();
		}
	}
public: // functions
  void ReadFromEDC(ElementDataContainer* edc);  // read a sub-edc (part) containing several hexblocks !defined format!
	void Cut();                                   // cuts the blocks in the list - prepares consistend meshing of the part
	void FindCuttingPlanes();                     // computes all cutting planes (sorted), arrays include overall limits
	void ComputeOccupationArray();                // computes which block occupies the volume-segment (the latter block always gets priority)
  void ComputeDivisionsArray();                 // computes the divisions for the blocks
	void SetBlockData();                          // writes computed sectioning to the blocks
	void Generate();                              // creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
	void RefreshLists(int flag_faces = 0);        // recomputes the lists of elements, nodes, etc for the entire part
public: // access  
	SegmentedQuadBlock& Block(int i) { return (SegmentedQuadBlock&) _Block(i); }
	SegmentedQuadBlock& operator() (int i) { return (SegmentedQuadBlock&) _Block(i); }
private: // own variables
};


//+  end MeshedPart and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ class FEMesh_Generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: encapsulate all kinds of generation functions
class FEMesh_Generator
{
public: // lifecycle
	FEMesh_Generator():mbs(0),mesh(0) {}
	FEMesh_Generator(FEMesh_Generator& other) { CopyFrom(other); }
	~FEMesh_Generator() {}
public: // lifecycle II
	FEMesh_Generator(MBS* mbsi, FEMesh* meshi)
	{
		SetFEMesh_Generator(mbsi,meshi);
	}
	void SetFEMesh_Generator(MBS* mbsi, FEMesh* meshi)
	{
		mbs = mbsi;
		mesh = meshi;
	}
	void CopyFrom(FEMesh_Generator& other)
	{
		mbs = other.mbs;
		mesh = other.mesh;
	}

public: // access
	MBS* GetMBS() { return mbs; }
	FEMesh* GetMesh(); 

public: // access on members of the FEMesh class
	virtual int NN();
	virtual Vector3D GetPoint3D(int i);
	virtual int AddNodeCheck(Vector3D pos, int domain = 1, double tol=1E-10);
	virtual int AddNodeCheck(Vector2D pos, int domain = 1, double tol=1E-10);

	virtual FEElement& GetElement(int i);
	virtual int AddElement(FEElement& fep);

	virtual TArray<FEMesh_Face>& Faces();

	virtual Box3D& GetNodeBox();
	virtual SearchTree& NodesTree();
	virtual void ComputeNodeBox(IVector& subset);
	virtual void ResizeNodeSearchTree(Box3D addbox, int addnodes);
// to be removed soon....
	virtual int NLoads();
	virtual int AddLoad(FEMesh_Load& loadi);
	virtual void LinearToQuadratic(IVector& subset);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Hexahedral Blocks

// extern calls with full set of parameters
public: 
// creates nodes, elements, outerfaces and add them to FEMesh,
// returns the respective numbers of (nodes/elements/faces) in FEMesh in ret_*
// 6 faces: 1=xmin, 2=xmax, 3=ymin, 4=ymax, 5=zmin, 6=zmax
// 3 directions: 1=x, 2=y, 4=z -> 0..7 combinations for xyz on/off
	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		  IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																			MBSLoad& bodyload, int bodyloadactive,
		  																IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			int order=1);

	// generates a block of prism elements (hexahedrals split into 2 prisms each) ( including nodes, constraints, ... )
	virtual int GenerateHexahedralBlock_Prisms(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			       IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		         IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       	 MBSLoad& bodyload, int bodyloadactive,
		  																       IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			       int order=1);
	//$EK - 2013-03-04
	// generates a block of real prism elements (hexahedrals split into 2 prisms each) ( including nodes, constraints, ... )
	virtual int GenerateHexahedralBlock_Prisms_new(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			       IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		         IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		       	 MBSLoad& bodyload, int bodyloadactive,
		  																       IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			       int order=1);



// generates a block of pyramid elements (hexahedrals split into 3 pyramids each) ( including nodes, constraints, ... )
	virtual int GenerateHexahedralBlock_Pyrams(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			       IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		         IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		         MBSLoad& bodyload, int bodyloadactive,
		  																       IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			       int order=1);
	//$EK - 2013-03-04
	// generates a block of real pyramid elements (hexahedrals split into 3 pyramids each) ( including nodes, constraints, ... )
	virtual int GenerateHexahedralBlock_Pyramids_new(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			       IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		         IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		         MBSLoad& bodyload, int bodyloadactive,
		  																       IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			       int order=1);

// generates a block of tetrahedral elements (hexahedrals split into 5 tetrahedrals each) ( including nodes, constraints, ... )
	virtual int GenerateHexahedralBlock_Tetras(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,
																			       IVector& ret_elementnumbers, IVector& ret_nodenumbers, TArray<IVector*>& ret_faces,
																		         IVector& constraintfaces_constainttypes, IVector& constraintfaces_directions,
																		         MBSLoad& bodyload, int bodyloadactive,
		  																       IVector& loadfaces_active, Vector& loadfaces_loadstrength,														
																			       int order=1);

// generates a DEFORMED block of hexaherals elements (return values, no boundary conditions)
	virtual int GenerateHexahedralBlock(TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces);			

// abbreviated calls
public: 
// generates a block of hexaherals elements (no return values, no boundary conditions)
	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color);

// generates a block of hexaherals elements (return values, no boundary conditions)
	virtual int GenerateHexahedralBlock(Box3D block, int3 divisions_xyz, int bodynr, int matnr, Vector3D color,           
																		IVector& blockelements, IVector& blocknodes, TArray<IVector*>& blockfaces);			

// generates a DEFORMED block of hexaherals elements (return values, no boundary conditions)
	virtual int GenerateHexahedralBlock(TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr, int matnr, Vector3D color);

// functions to create nodes
protected: 
// creates nodes for hexahedral block - adds nodes to mesh.points
  virtual int CreateHexBlock_Nodes(IVector& blocknodes, Box3D block, int3 divisions_xyz, int bodynr=1, int order=1);
// creates nodes for DEFORMED hexahedral block - adds nodes to mesh.points
	virtual int CreateHexBlock_Nodes(IVector& blocknodes, TArray<Vector3D>& corners, int3 divisions_xyz, int bodynr=1, int order=1);

// functions to create elements
protected: 
// returns Nodelist of a single (hexahedral) element of the block
  virtual int CreateHexBlock_HexNodes(IVector& points, IVector& blocknodes, int firstnode, int dx, int dy, int dz, int order=1);
// creates Hexahedrals from nodelist - adds FEHex to mesh.elements
	virtual int CreateHexBlock_Hexes(IVector& ret_elements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);	
// creates Prisms from nodelist - adds 2 FEHex to mesh.elements
	virtual int CreateHexBlock_Prisms(IVector& ret_elements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);
	// $EK 2013-03-04 creates Prisms from nodelist - adds 2 FEPrism to mesh.elements 
	virtual int CreateHexBlock_Prisms_new(IVector& ret_elements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);
	// creates Pyramids from nodelist - adds 3 FEHex to mesh.elements
	virtual int CreateHexBlock_Pyrams(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);
	//$EK 2013-03-04 added for real FEPyramids
	// creates Pyramids from nodelist - adds 3 FEPyramids to mesh.elements
	virtual int CreateHexBlock_Pyramids_new(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);
// creates Pyramids from nodelist - adds 5 FETet to mesh.elements
	virtual int CreateHexBlock_Tetras(IVector& blockelements, IVector& blocknodes, int3 divisions_xyz, int bodynr, int matnr, Vector3D color, int order=1);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Cylindrical Blocks

// extern calls with full set of parameters
public: 
// generates a 1/4 disc of hexaherals elements, with or without hole in the center
	virtual int GenerateDisc(double radius, double radius_hole, double height, int divisions_angle, int divisions_radial, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes);

	// generates a 1/4 cylinder of hexaherals elements 
	virtual int GenerateCylinder(double radius, double small_radius, double height, int divisions_angle, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes);

// generates a 1/4 hollow cylinder of hexaherals elements 
	virtual int GenerateHollowCylinder(double radius, double radius_hole, double height, int divisions_angle, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes);

// generates a 1/4 (hollow) cylinder of hexaherals elements 
	virtual int GenerateCylinder(Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r, int divisions_angle, int divisions_height, 
								int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes);

// generates a 1/4 hollow cylinder of hexaherals elements with increasing number of divisions/rad 
	virtual int GenerateHollowCylinderDoubleDiv(double radius, double radius_hole, double height, int divisions_angle, int divisions_height, Vector3D pos,
							 int bodynr, int matnr, Vector3D color,
							 IVector& blockelements, IVector& blocknodes);

// generates a 1/4 hollow cylinder of hexaherals elements with increasing number of divisions/rad 
	virtual int GenerateHollowCylinderDoubleDiv(Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r, int divisions_angle, int divisions_height, 
								int bodynr, int matnr, Vector3D color, IVector& blockelements, IVector& blocknodes);

// generates a rotor based on a mesh2d by rotating the mesh2d around the specified rotation axis
	virtual int GenerateRotorOutOfMesh2D(FEMesh* mesh2d, int rotation_axis_number, int angular_segments, double final_angle_deg);

// abbreviated calls
public: 

// functions to create nodes
protected: 
// creates nodes for cylinder 
  virtual int CreateCylinder_Mesh2D(double radius, double small_radius, int divisions_angle, TArray<double>& x, TArray<double>& y);
// creates nodes for an hollow cylinder 
  virtual int CreateHollowCylinder_Mesh2D(double radius, double radius_hole, int divisions_angle, TArray<double>& x, TArray<double>& y);
// creates nodes for an hollow cylinder with increasing number of divisions/rad
  virtual int CreateHollowCylinderDoubleDiv_Mesh2D(double radius, double radius_hole, int divisions_angle, TArray<double>& x, TArray<double>& y);

// create 3D nodes based on 2D mesh. Used by functions GenerateCylinder, GenerateHollowCylinder, GenerateHollowCylinderDoubleDiv
//	- adds nodes to mesh.points
	virtual int Create3DNodesOutOf2DCoordinates(TArray<double>& x, TArray<double>& y, double height, int divisions_height, Vector3D pos, int bodynr, IVector& blocknodes,MathFunction* outer=NULL, MathFunction* inner=NULL, MathFunction* left = NULL,MathFunction* right=NULL);
	virtual int Create3DNodesOutOf2DCoordinates(TArray<double>& x, TArray<double>& y, double height, int divisions_height, int bodynr,IVector& blocknodes,Vector3D p1l,Vector3D p1r,Vector3D p2l,Vector3D p2r);
	virtual void Distort3DMesh(TArray<Vector3D>& nodes, int independent_coord, int deformed_coord, MathFunction* outer, double ref_outer, MathFunction* inner, double ref_inner, int deformed_coord2=0, int independent_coord2=0);
	virtual void Distort3DMesh(TArray<Vector3D>& nodes, int independent_coord, int deformed_coord, MathFunction* outer, MathFunction* ref_outer, MathFunction* inner, MathFunction* ref_inner, int deformed_coord2=0, int independent_coord2=0);


// functions to create elements
protected: 
// creates Hexaherals for a cylinder from nodelist - adds FEHex to mesh.elements 
	virtual int CreateCylinder_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color);
// creates Hexaherals for an hollow cylinder from nodelist - adds FEHex to mesh.elements 
	virtual int CreateHollowCylinder_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color);
// creates Hexaherals for an hollow cylinder with increasing number of divisions/rad from nodelist - adds FEHex to mesh.elements 
	virtual int CreateHollowCylinderDoubleDiv_Hexes(IVector& blockelements, IVector& blocknodes, int divisions_angle, int divisions_height, int bodynr, int matnr, Vector3D color);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ Functions for Quadrilateral Blocks

// abbreviated calls
public: 
// generates a block of quadrilateral elements (no return values, no boundary conditions)
	virtual int GenerateQuadrilateralBlock(Box2D block, double thickness, int2 divisions_xy, int bodynr, int matnr, Vector3D color);
	virtual int GenerateQuadrilateralBlock(Box2D block, int2 divisions_xy, int bodynr, int matnr, Vector3D color) { return GenerateQuadrilateralBlock(block, 0., divisions_xy, bodynr, matnr, color); }


// functions to create nodes
protected: 
// creates nodes for quadrilateral block - adds nodes to mesh.points 
	virtual int CreateQuadBlock_Nodes(IVector& blocknodes, Box2D block, int2 divisions_xy, int bodynr, int order);

// functions to create elements
protected: 
// returns Nodelist of a single (quadrilateral) element of the block
	virtual int CreateQuadBlock_QuadNodes(IVector& points, IVector& blocknodes, int firstnode, int dx, int dy, int order);

// creates Quadrilaterals from nodelist - adds FEQuadx to mesh.elements
	virtual int CreateQuadBlock_Quads(IVector& blockelements, IVector& blocknodes, double thickness, int2 divisions_xy, int bodynr, int matnr, Vector3D color, int order);




//TODO: kill off all these functions... outerfaces etc should be handeld in class SegmentedBlock
 // outer faces and constraints applied on them 
protected:
	// create outer faces - adds FEQuad to mesh.faces
	virtual int CreateHexBlock_Hexes_OneFace(IVector& ret_faces, IVector& subset_elements, int side);
// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
	virtual int CreateHexBlock_Hexes_AllFaces(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
	virtual int CreateHexBlock_Prisms_AllFaces(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
	// $EK 2013-03-04 - for real Prismatic elements
	// creates all faces of block - 6 individual lists - adds FEQuad, FETrig to mesh.faces   
	virtual int CreateHexBlock_Prisms_AllFaces_new(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces   
	virtual int CreateHexBlock_Pyrams_AllFaces(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
	// $EK  2013-03-04 - for real Pyramids
	// creates all faces of block - 6 individual lists - adds FEQuad, FETrig to mesh.faces   
	virtual int CreateHexBlock_Pyramids_AllFaces_new(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
	// creates all faces of block - 6 individual lists - adds FEQuad to mesh.faces  
	virtual int CreateHexBlock_Tetras_AllFaces(TArray<IVector*>& ret_faces, IVector& blockelements, int3 divisions_xyz);
	
// constraints for one blockface - not implemented yet
  virtual int CreateHexBlock_OneFaceConstraints(int type, int direction, IVector& subset_outerfaces);
// constraints for all blockfaces - not implemented yet
  virtual int CreateHexBlock_AllFaceConstraints(IVector& types, IVector& directions, TArray<IVector*>& blockfaces);

// apply a bodyload on all elements of the block
	virtual int CreateHexBlock_BodyLoads(MBSLoad& bodyload, IVector& blockelements);

// apply face-normal area loads to outer faces
	virtual int CreateHexBlock_LoadFaces(IVector& loadfaces_active, Vector& loadfaces_loadstrength, TArray<IVector*>& blockfaces, Box3D block);

// find all outer elements - 6 individual lists
	virtual int CreateHexBlock_Hexes_ComputeOuterElements(TArray<IVector*>& ret_outerelems, int3 divisions_xyz);
// find all outer elements - 6 individual lists
	virtual int CreateHexBlock_Prisms_ComputeOuterElements(TArray<IVector*>& ret_outerelems, int3 divisions_xyz);
// find all outer elements - 6 individual lists
	virtual int CreateHexBlock_Pyrams_ComputeOuterElements(TArray<IVector*>& ret_outerelems, int3 divisions_xyz);
// find all outer elements - 6 individual lists
	virtual int CreateHexBlock_Tetras_ComputeOuterElements(TArray<IVector*>& ret_outerelems, int3 divisions_xyz);




private: // variables
	MBS* mbs;
	FEMesh* mesh;

};


#endif //FEMESH_AUX__H