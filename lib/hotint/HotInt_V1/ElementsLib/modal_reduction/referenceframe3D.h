//#**************************************************************
//#
//# filename:             ReferenceFrame3D.h
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
 
#ifndef REFERENCEFRAME3D__H
#define REFERENCEFRAME3D__H

class Rigid3DKardan;

//class BaseReferenceFrame3D: public Element
//{
//public:
//	//ReferenceFrame3D():Element() {mbs = NULL;};
//	BaseReferenceFrame3D(MBS* mbsi): Element(mbsi)
//	{ ;	};
//
//	virtual void Initialize() = 0;
//
//	virtual void InitializeSearchtree(const Vector3D& size1, const Vector3D& size2, int ix, int iy, int iz) = 0;
//	virtual const char* GetElementSpec() const {return "BaseReferenceFrame3D";}
//
//	virtual void EvalM(Matrix& m, double t) = 0; 
//	virtual void EvalF2(Vector& f, double t) = 0;
//
//	virtual int NFFRFElements() const = 0;
//	virtual int AddFFRFElement(int elnum) = 0;
//	virtual int GetFFRFElementNum(int i) const = 0;
//	virtual const Element& GetFFRFElement(int i) const  = 0;
//	virtual Element& GetFFRFElement(int i)  = 0;
//
//	virtual void SetResortConstraint(int i)  = 0;
//	virtual int IsResortConstraint() const  = 0;
//
//	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	//CMS
//	virtual void SetIsACRS(int i=1)  = 0;
//	virtual int IsACRS() const  = 0;
//	virtual int IsCMS() const  = 0;
//	virtual const double& GetXactFull(int i) const  = 0;
//	virtual double& GetXactFull(int i)  = 0;
//	virtual const double& GetDrawValueFull(int i) const  = 0;
//
//	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	//most functions identical to RIGID!
//	//
//	virtual int FlexDOF() const  = 0;
//	virtual int IsRigid() const  = 0;
//	virtual int SOSRigid() const  = 0;
//
//	//end: from RIGID
//	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//	// zero MBS nodes
//	virtual int NNodes() const  = 0;
//	// own nodes, which do not belong to MBS
//	virtual int NCMSNodes() const  = 0;
//	virtual int AddNode(Node* n) = 0;
//	virtual Node& GetNode(int i)  = 0;
//	virtual const Node& GetNode(int i) const  = 0;
//	int GetNode(Vector3D p_loc) const = 0;
//
//	virtual Vector3D GetNodePos(int i) const  = 0;
//
//	virtual Vector3D GetNodeVel(int i) const  = 0;
//
//	virtual Vector3D GetNodePosD(int i) const = 0;
//
//	virtual void DrawElement() = 0;;
//
//	virtual void PrintToAbaqus(ofstream& os)  = 0;
//
//	virtual void GetElementData(ElementDataContainer& edc) = 0;; 		//fill in all element data
//	virtual int SetElementData(ElementDataContainer& edc) = 0;; //set element data according to ElementDataContainer
//};

template <class RIGID>
class ReferenceFrame3D: public RIGID
{
public:
	//ReferenceFrame3D():Element() {mbs = NULL;};
	ReferenceFrame3D(MBS* mbsi):RIGID(mbsi), nodes(), searchtree(), FFRFelements(), isACRS(0), resortconstraint(1)
	{
		elementname = GetElementSpec();
	};
	ReferenceFrame3D(const ReferenceFrame3D& e):RIGID(e.mbs), nodes(), searchtree(), isACRS(0), resortconstraint(1) {CopyFrom(e);};
	//To be overwritten in derived class:

	ReferenceFrame3D(MBS* mbsi, const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, 
		const Vector3D& sizei, const Vector3D& coli):RIGID(mbsi), nodes(), FFRFelements(), searchtree()
	{
		elementname = GetElementSpec();
		SetReferenceFrame(p, v, phi, phip, sizei, coli);
	};

	// Set-Function,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	void SetReferenceFrame(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip,
		const Vector3D& sizei, const Vector3D& coli);

	// SetFunction for Rigid3DKardan,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	// rs ........ Kardan angle rotation sequence
	void SetReferenceFrameKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip,
		const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs);

	virtual Element* GetCopy()
	{
		Element* ec = new ReferenceFrame3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		RIGID::CopyFrom(e);
		const ReferenceFrame3D& ce = (const ReferenceFrame3D&)e;
		size = ce.size;

		FFRFelements.SetLen(0);
		for (int i=1; i<=ce.FFRFelements.Length(); i++)
		{
			FFRFelements.Add(ce.FFRFelements(i));
		}
		searchtree = ce.searchtree;

		nodes.SetLen(0);
		for (int i=1; i<=ce.nodes.Length(); i++)
		{
			nodes.Add(ce.nodes(i));
		}	

		isACRS = ce.isACRS;
		resortconstraint = ce.resortconstraint;

	}

	virtual void AddFEMesh(const FEMeshInterface& femesh);

	virtual void Initialize() 
	{
		RIGID::Initialize();
	};

	virtual void InitializeSearchtree(const Vector3D& size1, const Vector3D& size2, int ix, int iy, int iz)
	{
		searchtree = SearchTree(ix, iy, iz, Box3D(size1, size2)); //define number of cubes and maximum size of search tree
	}
	virtual const char* GetElementSpec() const {return "ReferenceFrame3D";}

	virtual void EvalM(Matrix& m, double t); 
	virtual void EvalF2(Vector& f, double t);

	virtual int NFFRFElements() const {return FFRFelements.Length();}
	virtual int AddFFRFElement(int elnum) {return FFRFelements.Add(elnum);}
	virtual int GetFFRFElementNum(int i) const {return FFRFelements(i);}
	virtual const Element& GetFFRFElement(int i) const {return GetMBS()->GetElement(FFRFelements(i));}
	virtual Element& GetFFRFElement(int i) {return GetMBS()->GetElement(FFRFelements(i));}

	virtual void SetResortConstraint(int i) {resortconstraint = i;}
	virtual int IsResortConstraint() const {return resortconstraint;}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//CMS
	// AP: flag ACRS is replaced by GCMS-Element
	//virtual void SetIsACRS(int i=1) {isACRS = i;}
	//virtual int IsACRS() const {return isACRS;}
	virtual int IsCMS() const {return 0;}  // AP: returns 1 if frame is BaseCMSElement, CMSElement, GCMSElement
	virtual int IsGCMS() const {return 0;} // AP: returns 1 if frame is GCMSElement
	virtual const double& GetXactFull(int i) const {assert(0); static double d; return d;};
	virtual double& GetXactFull(int i) {assert(0); static double d; return d;}; //should not be used ...
	virtual const double& GetDrawValueFull(int i) const {assert(0); static double d; return d;};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//most functions identical to RIGID!
	//
	virtual int FlexDOF() const {return 0;}
	virtual int IsRigid() const {return 0;} //default value
	virtual int SOSRigid() const {return RIGID::SOS();}  // number of degrees of freedom from the rigid

	//end: from RIGID
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// zero MBS nodes
	virtual int NNodes() const {return 0;}
	// own nodes, which do not belong to MBS
	virtual int NCMSNodes() const {return nodes.Length();}
	virtual int AddNode(Node* n); 
	virtual Node& GetNode(int i) {return *nodes(i);}
	virtual const Node& GetNode(int i) const {return *nodes(i);}
	int GetNode(Vector3D p_loc) const;

	virtual Vector3D GetNodePos(int i) const 
	{
		Vector3D prel = Vector3D(
			nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),
			nodes(i)->Pos().Y()+XG(nodes(i)->Get(2)),
			nodes(i)->Pos().Z()+XG(nodes(i)->Get(3)));
		return GetRefPos() + GetRotMatrix()*prel;
	}

	virtual Vector3D GetNodeVel(int i) const 
	{
		Vector3D prel = Vector3D(
			nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),
			nodes(i)->Pos().Y()+XG(nodes(i)->Get(2)),
			nodes(i)->Pos().Z()+XG(nodes(i)->Get(3)));
		Vector3D vrel = Vector3D(XGP(nodes(i)->Get(1)), XGP(nodes(i)->Get(2)), XGP(nodes(i)->Get(3)));
		return GetRefVel() + GetRotMatrix()*vrel + GetRotMatrixP()*prel;
	}

	virtual Vector3D GetNodePosD(int i) const 
	{
		Vector3D prel = Vector3D(
			nodes(i)->Pos().X()+XGD(nodes(i)->Get(1)),
			nodes(i)->Pos().Y()+XGD(nodes(i)->Get(2)),
			nodes(i)->Pos().Z()+XGD(nodes(i)->Get(3)));
		return GetRefPosD() + GetRotMatrixD()*prel;
	}

	virtual void DrawElement();

	virtual void PrintToAbaqus(ofstream& os) {};

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

protected:
	//Vector3D size;    //for graphics; do not draw if size=Vector3D(0,0,0)
	IVector FFRFelements; //elements connected to ReferenceFrame
	SearchTree searchtree; //for optimized node fill
	TArray<Node*> nodes;
	int resortconstraint;  //activates resorting of the DOF of the reference frame into the constraint part
	int isACRS;  //absolute coordinates reduced strain --> the frame rotation and translation is not taken into account in GetPos2D() etc.
};




#endif

