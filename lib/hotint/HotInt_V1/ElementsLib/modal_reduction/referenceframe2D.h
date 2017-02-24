//#**************************************************************
//#
//# filename:             referenceframe2D.h
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
 
#ifndef REFERENCEFRAME2D__H
#define REFERENCEFRAME2D__H


class ReferenceFrame2D: public Body2D
{
public:
	//ReferenceFrame2D():Element() {mbs = NULL;};
	ReferenceFrame2D(MBS* mbsi):Body2D(mbsi), nodes(), searchtree(), isACRS(0), resortconstraint(0)
	{
	};
	ReferenceFrame2D(const ReferenceFrame2D& e):Body2D(e.mbs), nodes(), searchtree(), isACRS(0), resortconstraint(0) {CopyFrom(e);};
	//To be overwritten in derived class:

	ReferenceFrame2D(MBS* mbsi, const Vector2D& p, const Vector2D& v, double phi, double phip, 
		const Vector3D& sizei, const Vector3D& coli):Body2D(mbsi), nodes(), FFRFelements(), searchtree()
	{
		x_init.SetLen(6);
		x_init(1) = p.X(); x_init(2) = p.Y();
		x_init(3) = phi;

		x_init(4) = v.X(); x_init(5) = v.Y();
		x_init(6) = phip;
		
		col = coli;
		size = sizei;

		resortconstraint = 0; //in old solver version this was activated (1)

		InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
		nodes.SetLen(0);
		isACRS = 0;
		draw_frame = 0;
	};

	virtual Element* GetCopy()
	{
		Element* ec = new ReferenceFrame2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const ReferenceFrame2D& ce = (const ReferenceFrame2D&)e;
		//size = ce.size;

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

		draw_frame = ce.draw_frame;
	}

	virtual const char* GetElementSpec() const {return "ReferenceFrame2D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer


	virtual void Initialize() 
	{
		Body2D::Initialize();
	};

	virtual void InitializeSearchtree(const Vector3D& size1, const Vector3D& size2, int ix, int iy, int iz)
	{
		searchtree = SearchTree(ix, iy, iz, Box3D(size1, size2)); //define number of cubes and maximum size of search tree
	}

	virtual void EvalM(Matrix& m, double t); 
	virtual void EvalF2(Vector& f, double t);

	virtual int IS() const {return 0;};  //implicit (algebraic) size
	virtual int SOS() const {return 3;}; //size of second order equations, len(u)
	virtual int SOSowned() const {return SOS();}; //size of second order equations, len(u)

	virtual int SOSowned_RS() const {return SOS()-3*resortconstraint;}; //size of second order equations, len(u)
	virtual int IS_RS() const {return 3*resortconstraint;}; //size of second order equations, len(u)

	virtual int FlexDOF() const {return 0;}

	virtual int Dim() const {return 2;} //default value
	virtual int IsRigid() const {return 0;} //default value

	virtual int NFFRFElements() const {return FFRFelements.Length();}
	virtual int AddFFRFElement(int elnum) {return FFRFelements.Add(elnum);}
	virtual int GetFFRFElementNum(int i) const {return FFRFelements(i);}
	virtual const Element& GetFFRFElement(int i) const {return GetMBS()->GetElement(FFRFelements(i));}
	virtual Element& GetFFRFElement(int i) {return GetMBS()->GetElement(FFRFelements(i));}

	virtual void SetResortConstraint(int i) {resortconstraint = i;}
	virtual int IsResortConstraint() const {return resortconstraint;}

	virtual const int& DrawFrame() const {return draw_frame;}
	virtual int& DrawFrame() {return draw_frame;}

	virtual void SetIsACRS(int i=1) {isACRS = i;}
	virtual int IsACRS() const {return isACRS;}
	virtual const double& GetXactFull(int i) const {assert(0); static double d; return d;};
	virtual double& GetXactFull(int i) {assert(0); static double d; return d;}; //should not be used ...
	virtual const double& GetDrawValueFull(int i) const {assert(0); static double d; return d;};

	//get angle in 2D
	virtual double GetAngle2D() const {return XG(3);};
	virtual double GetAngle2DP() const {return XGP(3);};
	virtual double GetAngle2DD() const {return XGD(3);};
	virtual double GetAngle2DPD() const {return XGPD(3);};

	virtual Vector2D GetRefPos2D() const 
	{
		return Vector2D(XG(1),XG(2));
	};
	virtual Vector2D GetRefPos2DD() const 
	{
		return Vector2D(XGD(1),XGD(2));
	};
	virtual Vector2D GetRefVel2D() const {return Vector2D(XGP(1),XGP(2));};
	virtual Vector2D GetRefVel2DD() const {return Vector2D(XGPD(1),XGPD(2));};

	//get the rotated and translated position of a local point at the body
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const
  {
		return GetRotMatrix2D()*p_loc+GetRefPos2D();
	};
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const
  {
		return GetRotMatrix2DP()*p_loc+GetRefVel2D();
	};

	//get the rotated and translated position of a local point at the body
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const
  {
		return GetRotMatrix2DD()*p_loc+GetRefPos2DD();
	};

	//for loads:
	virtual void GetIntDuxDq(Vector& dudq) {};
	virtual void GetIntDuyDq(Vector& dudq) {};
	virtual void GetIntDuzDq(Vector& dudq) {};
	virtual void ApplyDprefdq(Vector& f, const Vector2D& x) {};
	virtual void ApplyDrotrefdq(Vector& f, const double& x) {};

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi)
	{
		dpdqi.SetSize(3,2);
		dpdqi.SetAll(0.);
		dpdqi(1,1) = 1.;
		dpdqi(2,2) = 1.;
		Vector2D dpdq3;
		dpdq3 = GetRotMatrixDphi2D() * ploc;
		dpdqi(3,1) = dpdq3.X();
		dpdqi(3,2) = dpdq3.Y();
	}

	virtual Matrix3D GetRotMatrix2D() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(GetAngle2D()); 
		double sinphi = sin(GetAngle2D()); 
		rot(1,1) = cosphi;
		rot(1,2) =-sinphi;
		rot(2,1) = sinphi;
		rot(2,2) = cosphi;

		return rot;
	}
	virtual Matrix3D GetRotMatrix2DP() const
	{
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(GetAngle2D()); 
		double sinphi = sin(GetAngle2D()); 
		double phip = GetAngle2DP(); 
		rot(1,1) =-phip*sinphi;
		rot(1,2) =-phip*cosphi;
		rot(2,1) = phip*cosphi;
		rot(2,2) =-phip*sinphi;

		return rot;
	}
	virtual Matrix3D GetRotMatrix() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double cosphi = cos(GetAngle2D()); 
		double sinphi = sin(GetAngle2D()); 
		rot.Get0(0,0) = cosphi;
		rot.Get0(0,1) =-sinphi;
		rot.Get0(0,2) = 0;
		rot.Get0(1,0) = sinphi;
		rot.Get0(1,1) = cosphi;
		rot.Get0(1,2) = 0;
		rot.Get0(2,0) = 0;
		rot.Get0(2,1) = 0;
		rot.Get0(2,2) = 1;

		return rot;
	}
	virtual Matrix3D GetRotMatrixP() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double cosphi = cos(GetAngle2D()); 
		double sinphi = sin(GetAngle2D()); 
		double phip = GetAngle2DP(); 
		rot.Get0(0,0) =-phip*sinphi;
		rot.Get0(0,1) =-phip*cosphi;
		rot.Get0(0,2) = 0;
		rot.Get0(1,0) = phip*cosphi;
		rot.Get0(1,1) =-phip*sinphi;
		rot.Get0(1,2) = 0;
		rot.Get0(2,0) = 0;
		rot.Get0(2,1) = 0;
		rot.Get0(2,2) = 0;

		return rot;
	}

	//drawing matrices:
	virtual Matrix3D GetRotMatrix2DD() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double phi = GetAngle2DD();
		rot(1,1) = cos(phi);
		rot(1,2) =-sin(phi);
		rot(2,1) = sin(phi);
		rot(2,2) = cos(phi);

		return rot;
	}
	virtual Matrix3D GetRotMatrix2DPD() const
	{
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(GetAngle2DD()); 
		double sinphi = sin(GetAngle2DD()); 
		double phip = GetAngle2DPD(); 
		rot(1,1) =-phip*sinphi;
		rot(1,2) =-phip*cosphi;
		rot(2,1) = phip*cosphi;
		rot(2,2) =-phip*sinphi;

		return rot;
	}
	virtual Matrix3D GetRotMatrixD() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double phi = GetAngle2DD(); 
		rot.Get0(0,0) = cos(phi);
		rot.Get0(0,1) =-sin(phi);
		rot.Get0(0,2) = 0;
		rot.Get0(1,0) = sin(phi);
		rot.Get0(1,1) = cos(phi);
		rot.Get0(1,2) = 0;
		rot.Get0(2,0) = 0;
		rot.Get0(2,1) = 0;
		rot.Get0(2,2) = 1;

		return rot;
	}
	virtual Matrix3D GetRotMatrixPD() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double phi  = GetAngle2DD(); 
		double phip = GetAngle2DPD(); 
		rot.Get0(0,0) =-phip*sin(phi);
		rot.Get0(0,1) =-phip*cos(phi);
		rot.Get0(0,2) = 0;
		rot.Get0(1,0) = phip*cos(phi);
		rot.Get0(1,1) =-phip*sin(phi);
		rot.Get0(1,2) = 0;
		rot.Get0(2,0) = 0;
		rot.Get0(2,1) = 0;
		rot.Get0(2,2) = 0;

		return rot;
	}

	virtual Matrix3D GetRotMatrixDphi2D() const
	{ 
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(GetAngle2D()); 
		double sinphi = sin(GetAngle2D()); 
		rot(1,1) =-sinphi;
		rot(1,2) =-cosphi;
		rot(2,1) = cosphi;
		rot(2,2) =-sinphi;

		return rot;
	}

	virtual const int& NodeNum(int i) const {return nodes(i)->NodeNum();}
	virtual int& NodeNum(int i) {return nodes(i)->NodeNum();}

	virtual int NNodes() const {return nodes.Length();}
	virtual int AddNode(Node* n); 
	virtual Node& GetNode(int i) {return *nodes(i);}
	virtual const Node& GetNode(int i) const {return *nodes(i);}
	int GetNode(Vector2D p_loc) const;

	virtual Vector2D GetNodePos2D(int i) const 
	{
		Vector2D prel = Vector2D(nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),nodes(i)->Pos().Y()+XG(nodes(i)->Get(2)));
		return GetRefPos2D() + GetRotMatrix2D()*prel;
	}

	virtual Vector2D GetNodeVel2D(int i) const 
	{
		Vector2D prel = Vector2D(nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),nodes(i)->Pos().Y()+XG(nodes(i)->Get(2)));
		Vector2D vrel = Vector2D(XGP(nodes(i)->Get(1)),XGP(nodes(i)->Get(2)));
		return GetRefVel2D() + GetRotMatrix2D()*vrel + GetRotMatrix2DP()*prel;
	}

	virtual Vector2D GetNodePos2DD(int i) const 
	{
		Vector2D prel = Vector2D(nodes(i)->Pos().X()+XGD(nodes(i)->Get(1)),nodes(i)->Pos().Y()+XGD(nodes(i)->Get(2)));
		return GetRefPos2DD() + GetRotMatrix2DD()*prel;
	}

	virtual void DrawElement();


	virtual void PrintToAbaqus(ofstream& os) {};

protected:
	//Vector3D size;    //for graphics; do not draw if size=Vector3D(0,0,0)
	IVector FFRFelements; //elements connected to ReferenceFrame
	SearchTree searchtree; //for optimized node fill
	TArray<Node*> nodes;
	int resortconstraint;  //activates resorting of the DOF of the reference frame into the constraint part
	int isACRS;  //absolute coordinates reduced strain --> the frame rotation and translation is not taken into account in GetPos2D() etc.
	int draw_frame;
};




#endif

