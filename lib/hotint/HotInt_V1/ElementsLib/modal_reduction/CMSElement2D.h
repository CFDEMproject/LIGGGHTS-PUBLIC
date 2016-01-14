//#**************************************************************
//#
//# filename:             CMSElement2D.h
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
 
#ifndef CMSElement2D__H
#define CMSElement2D__H


//this element contains the whole structure of the body (finite element numbers, nodes) in this class.
//all nodes have local indices in the CMSElement
class CMSElement2D: public ReferenceFrame2D
{
public:
	//CMSElement2D():Element() {mbs = NULL;};
	CMSElement2D(MBS* mbsi):ReferenceFrame2D(mbsi), 
		boundarynode(), boundarydofelem(), boundarydof()
	{
	};
	CMSElement2D(const CMSElement2D& e):ReferenceFrame2D(e.mbs),
		boundarynode(), boundarydofelem(), boundarydof()  {CopyFrom(e);};
	//To be overwritten in derived class:

	CMSElement2D(MBS* mbsi, const Vector2D& p, const Vector2D& v, double phi, double phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli):ReferenceFrame2D(mbsi),
		boundarynode(), boundarydofelem(), boundarydof()
	{
		x_init.SetLen(6); //temporary, change after modal synthesis
		x_init(1) = p.X(); x_init(2) = p.Y();
		x_init(3) = phi;

		x_init(4) = v.X(); x_init(5) = v.Y();
		x_init(6) = phip;

		col = coli;
		size = sizei;

		nbmodes = 0;
		nimodes = nimodesi;

		InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
		boundarydofelem.SetLen(0);
		boundarynode.SetLen(0);

	};

	virtual Element* GetCopy()
	{
		Element* ec = new CMSElement2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		ReferenceFrame2D::CopyFrom(e);
		const CMSElement2D& ce = (const CMSElement2D&)e;

		nbmodes = ce.nbmodes;
		nimodes = ce.nimodes;

		boundarydofelem.SetLen(0);
		for (int i=1; i<=ce.boundarydofelem.Length(); i++)
		{
			boundarydofelem.Add(ce.boundarydofelem(i));
		}	

		boundarydof.SetLen(0);
		for (int i=1; i<=ce.boundarydof.Length(); i++)
		{
			boundarydof.Add(ce.boundarydof(i));
		}	

		boundarynode.SetLen(0);
		for (int i=1; i<=ce.boundarynode.Length(); i++)
		{
			boundarynode.Add(ce.boundarynode(i));
		}	

		M = ce.M;
		K = ce.K;
		Msparse = ce.Msparse;
		Ksparse = ce.Ksparse;
		Mr = ce.Mr;
		Kr = ce.Kr;
		Hr = ce.Hr;

		Phi_CB = ce.Phi_CB;


	}

	virtual void Initialize();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//must be called before Assemble!!!!!
	//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	virtual void DoModalAnalysis(const TArray<int2>& fixednodes);

	//final matrices of FFRF-formulation:
	virtual void EvalM(Matrix& m, double t); 
	//stiffness+quadratic velocity vector:
	virtual void EvalF2(Vector& f, double t);

	virtual int FastStiffnessMatrix() const;
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components only of stiffness matrix, m might be larger

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual void GetI1(Vector& I1); //Shabana p. 209-211
	virtual double GetIkl(int k, int l);
	virtual void GetIbarkl(int k, int l, Vector& I1);
	virtual void GetSbar(Matrix& Sbar);
	virtual void GetSbarkl(int k, int l, Matrix& Sbar);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual int SOS() const {return NModes()+3;}; //size of second order equations, len(u)
	virtual int SOSowned() const {return SOS();}; //size of second order equations, len(u)
	virtual int SOSowned_RS() const {return SOS()-3;}; //size of second order equations, len(u)
	virtual int IS_RS() const {return 3;}; //size of second order equations, len(u)


	virtual int SOSFull() const 
	{
		int sosfull = 0;
		//nodal dof:
		for (int i=1; i <= nodes.Length(); i++) sosfull+=nodes(i)->SOS();

		//element dof:
		for (int i=1; i <= NFFRFElements(); i++) sosfull+=GetFFRFElement(i).SOSowned();

		return sosfull;
	}

	virtual const Element& GetFFRFElement(int i) const {return GetMBS()->GetAuxElement(FFRFelements(i));}
	virtual Element& GetFFRFElement(int i) {return GetMBS()->GetAuxElement(FFRFelements(i));}

	virtual int FlexDOF() const {return NModes();}
	virtual int NModes() const {return NBModes()+NIModes();}

	virtual int NBModes() const {return nbmodes;} //boundary modes
	virtual int NIModes() const {return nimodes;} //internal modes

	virtual void SetNIModes(int i) {nimodes = i;} //internal modes



	virtual void AddBoundaryDOF(int elem, int locelemdof) {boundarydofelem.Add(int2(elem,locelemdof));}
	virtual void AddBoundaryNode(int node) {boundarynode.Add(node);}


	virtual const double& GetXactFull(int i) const;
	virtual double& GetXactFull(int i); //should not be used ...
	virtual const double& GetDrawValueFull(int i) const;

	//get angle in 2D
	virtual double GetAngle2D() const { return XG(3+NModes());};
	virtual double GetAngle2DP() const { return XGP(3+NModes());};
	virtual double GetAngle2DD() const { return XGD(3+NModes());};
	virtual double GetAngle2DPD() const { return XGPD(3+NModes());};

	virtual Vector2D GetRefPos2D() const { return Vector2D(XG(1+NModes()),XG(2+NModes()));};
	virtual Vector2D GetRefPos2DD() const { return Vector2D(XGD(1+NModes()),XGD(2+NModes()));};
	virtual Vector2D GetRefVel2D() const { return Vector2D(XGP(1+NModes()),XGP(2+NModes()));};
	virtual Vector2D GetRefVel2DD() const { return Vector2D(XGPD(1+NModes()),XGPD(2+NModes()));};
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi); 
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);
	virtual void AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f);   // f += dpdq*lambda

	virtual void GetIntDuDq(Matrix& dudq); //in fact it is DuDq Transposed

	virtual void GetH(Matrix& H)
	{
		H = Hr;
	}

	virtual Box3D GetElementBox() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBox());
		}
		return box;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBoxD());
		}
		return box;
	}

private:


	//CMS:
	int nbmodes;
	int nimodes;

	Vector xg;

	Matrix Phi_CB; //transformation matrix reduced to full dof
	TArray<int2> boundarydofelem; //element number and local dof which are fixed
	TArray<int> boundarydof; //mark dof with 1 which are fixed
	TArray<int> boundarynode;     //node number which is fixed

	Matrix M, K, Mr, Kr, Hr; //full and reduced system matrices
	SparseMatrix Msparse, Ksparse;
	double volume; //dor intdudq

	//double xactfullvalue; //dummy variable, for non-const call of XActFull()
	//double xactfullvalueD; //dummy variable, for non-const call of XActFull()

	///+++++++++++++++++++++
	//store matrices for FFRF:
	Matrix Sbar_tilde;
	Matrix SbarS;
	SparseMatrix Sbar_tildeSM;
	SparseMatrix SbarSM;
	Vector Ibar11S,Ibar21S,Ibar12S,Ibar22S;
	Vector I1S;
	double I1122S;
	//++++++++++++++++++++++
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class ACRSElement2D: public ReferenceFrame2D
{
public:
	//ACRSElement2D():Element() {mbs = NULL;};
	ACRSElement2D(MBS* mbsi):ReferenceFrame2D(mbsi),
		boundarynode(), boundarydof()
	{
	};
	ACRSElement2D(const ACRSElement2D& e):ReferenceFrame2D(e.mbs),
		boundarynode(), boundarydof()  {CopyFrom(e);};
	//To be overwritten in derived class:

	ACRSElement2D(MBS* mbsi, const Vector3D& searchtreesizei, const Vector3D& coli, const int nimodesi = 0):
	ReferenceFrame2D(mbsi), boundarynode(), boundarydof()
	{
		nbmodes = 0;
		nimodes = 2*nimodesi; //take modes + orthogonal (rotated) modes
		isCMS = 0;
		if (nimodesi != 0) isCMS = 1;

		useSparseMK = 1;
		usesparsepre = 1;
		SetIsACRS(1);

		x_init.SetLen(0); //temporary, change after assembly

		col = coli;
		size = searchtreesizei;


		InitializeSearchtree(-0.1*size,1.1*size,10,10,10);

		boundarynode.SetLen(0);

		sosfull = 0;

		initphi = 0;
		initrot.SetSize(2,2);
		initrot.SetAll(0);
	};

	virtual Element* GetCopy()
	{
		Element* ec = new ACRSElement2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		ReferenceFrame2D::CopyFrom(e);
		const ACRSElement2D& ce = (const ACRSElement2D&)e;

		nbmodes = ce.nbmodes;
		nimodes = ce.nimodes;



		boundarydof.SetLen(0);
		for (int i=1; i<=ce.boundarydof.Length(); i++)
		{
			boundarydof.Add(ce.boundarydof(i));
		}	

		boundarynode.SetLen(0);
		for (int i=1; i<=ce.boundarynode.Length(); i++)
		{
			boundarynode.Add(ce.boundarynode(i));
		}	

		//M = ce.M;
		//K = ce.K;
		Msparse = ce.Msparse;
		Ksparse = ce.Ksparse;
		Mr = ce.Mr;
		Kr = ce.Kr;
		Hr = ce.Hr;

		Phi_CB = ce.Phi_CB;

		isCMS = ce.isCMS;
		sosfull = ce.sosfull;
		refnode1 = ce.refnode1;
		refnode2 = ce.refnode2;

		useSparseMK = ce.useSparseMK;
		usesparsepre = ce.usesparsepre;

		newmode = ce.newmode;
		reduced1 = ce.reduced1;
		reduced2 = ce.reduced2;

		initphi = ce.initphi;
		initrot = ce.initrot;
	}

	virtual void Initialize();

	//compute reduced matrices:
	virtual void FinishAssembly();
	virtual void DoModalAnalysis();


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//absolute coordinate reduced strain formulation
	//different mass and stiffness matrix and variation of position or integrals
	virtual void SetIsCMS(int i=1) {isCMS = i;}
	virtual int IsCMS() const {return isCMS;}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//sparse operations for ACRS
	virtual int ElementBandwidth() const;
	virtual void SetUseSparseMK(int i=1) {useSparseMK = i;} //only for CMS element
	virtual int UseSparseM() const {return useSparseMK;}
	virtual int UseSparseK() const {return useSparseMK;}
	virtual int UseSparseMK() const {return useSparseMK;} //only for CMS element
	virtual void AddMSparse(SparseMatrix& m, double t);  //add sparse matrix into full system matrix
	virtual void AddKSparse(SparseMatrix& m, double t);  //add sparse matrix into full system matrix
	virtual int TransformJacApply() const;
	virtual void ApplyTransform(const Vector& v, Vector& Av, int mode); //compute Av=A^T*v in mode==0 and Av=A*v in mode==1

	//define position of 2 nodes which define the frame, boundary node must be added as well!!!
	virtual void SetReferenceNodes(int node1, int node2) {refnode1 = node1; refnode2 = node2; };
	virtual int RefNode1() const {return refnode1;}
	virtual int RefNode2() const {return refnode2;}
	virtual void ApplyRotation(const Matrix3D& rot, Vector& u) const;
	virtual void ComputeURigid(const Matrix3D& A, const Vector2D& uref1, Vector& urigid) const;
	virtual void ComputeDAijDq(int i, int j, Vector& v);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//must be called before Assemble!!!!!
	//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	virtual void DoModalAnalysis(int ninternalmodes) {};

	//constant stiffness matrix
	virtual void EvalM(Matrix& m, double t); 
	//linear and nonlinear part
	virtual void EvalF2(Vector& f, double t);

	virtual int FastStiffnessMatrix() const;
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components only of stiffness matrix, m might be larger


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual int SOS() const {if (isCMS) return NModes(); else return sosfull;}; //size of second order equations, len(u)
	virtual int SOSowned() const {return SOS();}; //size of second order equations, len(u)
	virtual int SOSowned_RS() const {return SOSowned();}; //for resort, the implicit element variables belong to SOS2
	virtual int IS_RS() const {return IS();};  //implicit (algebraic) size


	virtual int SOSFull() const
	{
		int sosfull = 0;
		//nodal dof:
		for (int i=1; i <= nodes.Length(); i++) sosfull+=nodes(i)->SOS();

		//element dof:
		for (int i=1; i <= NFFRFElements(); i++) sosfull+=GetFFRFElement(i).SOSowned();

		return sosfull;
	}

	virtual const Element& GetFFRFElement(int i) const {return GetMBS()->GetAuxElement(FFRFelements(i));}
	virtual Element& GetFFRFElement(int i) {return GetMBS()->GetAuxElement(FFRFelements(i));}

	virtual int FlexDOF() const {return SOS();}
	virtual int NModes() const {return NBModes()+NIModes();}

	virtual int NBModes() const {return nbmodes;} //boundary modes
	virtual int NIModes() const {return nimodes;} //internal modes

	virtual void SetNIModes(int i) {nimodes = i;} //internal modes

	virtual void AddBoundaryNode(int node) {boundarynode.Add(node);}


	virtual const double& GetXactFull(int i) const;
	virtual double& GetXactFull(int i); //should not be used ...
	virtual const double& GetDrawValueFull(int i) const;

	//get angle in 2D

	virtual double GetAngle2D() const 
	{ 
		Vector2D p1 = GetNodePos2D(RefNode1());
		Vector2D p2 = GetNodePos2D(RefNode2());
		double Lact = (p1-p2).Norm();
		double s = 1./Lact*(p2.Y()-p1.Y());
		double c = 1./Lact*(p2.X()-p1.X());
		return atan2(s,c)-initphi;
	};
	virtual double GetAngle2DP() const { assert(0); return 0;};
	virtual double GetAngle2DD() const 
	{ 
		Vector2D p1 = GetNodePos2DD(RefNode1());
		Vector2D p2 = GetNodePos2DD(RefNode2());
		double Lact = (p1-p2).Norm();
		double s = 1./Lact*(p2.Y()-p1.Y());
		double c = 1./Lact*(p2.X()-p1.X());
		return atan2(s,c)-initphi;
	};
	virtual double GetAngle2DPD() const { assert(0); return 0;};

	virtual Matrix3D GetRotMatrix2D() const
	{
		Matrix3D rot;
		rot.SetSize(Dim(),Dim());
		Vector2D p1 = GetNodePos2D(RefNode1());
		Vector2D p2 = GetNodePos2D(RefNode2());
		double Lact = (p1-p2).Norm();
		double s = 1./Lact*(p2.Y()-p1.Y());
		double c = 1./Lact*(p2.X()-p1.X());

		rot(1,1) = c;
		rot(1,2) =-s;
		rot(2,1) = s;
		rot(2,2) = c;

		rot = initrot*rot;

		return rot;
	}

	virtual Matrix3D GetRotMatrix2DD() const
	{
		Matrix3D rot;
		rot.SetSize(Dim(),Dim());
		Vector2D p1 = GetNodePos2DD(RefNode1());
		Vector2D p2 = GetNodePos2DD(RefNode2());
		double Lact = (p1-p2).Norm();
		double s = 1./Lact*(p2.Y()-p1.Y());
		double c = 1./Lact*(p2.X()-p1.X());

		rot(1,1) = c;
		rot(1,2) =-s;
		rot(2,1) = s;
		rot(2,2) = c;

		rot = initrot*rot;

		return rot;
	}

	virtual Vector2D GetRefPos2D() const { return GetNodePos2D(RefNode1());};
	virtual Vector2D GetRefPos2DD() const {	return GetNodePos2DD(RefNode1());	};
	virtual Vector2D GetRefVel2D() const { return GetNodeVel2D(RefNode1());};
	virtual Vector2D GetRefVel2DD() const { return GetNodeVel2DD(RefNode1());};
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual Vector2D GetNodePosInit2D(int i) const 
	{
		return Vector2D(nodes(i)->Pos().X(),nodes(i)->Pos().Y());
	}
	virtual Vector2D GetNodeDisp2D(int i) const 
	{
		//**
		//only for boundary nodes!!!
		if (newmode)
		{
      if ((i != RefNode1() || !reduced1) && (i != RefNode2() || !reduced2))
				return Vector2D(XG(nodes(i)->Get(1))-XG(nodes(i)->Get(2)+1),XG(nodes(i)->Get(1)+1)+XG(nodes(i)->Get(2)));
			else
				return Vector2D(XG(nodes(i)->Get(1)),XG(nodes(i)->Get(1)+1));
		}
		else
			return Vector2D(XG(nodes(i)->Get(1)),XG(nodes(i)->Get(2)));
	}

	virtual Vector2D GetNodePos2D(int i) const 
	{
		//**
		//only for boundary nodes!!!
		if (newmode)
		{
      if ((i != RefNode1() || !reduced1) && (i != RefNode2() || !reduced2))
	      return Vector2D(nodes(i)->Pos().X()+XG(nodes(i)->Get(1))-XG(nodes(i)->Get(2)+1),
					nodes(i)->Pos().Y()+XG(nodes(i)->Get(1)+1)+XG(nodes(i)->Get(2)));
			else
			  return Vector2D(nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),
					nodes(i)->Pos().Y()+XG(nodes(i)->Get(1)+1));
		}
		else
			return Vector2D(nodes(i)->Pos().X()+XG(nodes(i)->Get(1)),nodes(i)->Pos().Y()+XG(nodes(i)->Get(2)));
	}

	virtual Vector2D GetNodeVel2D(int i) const 
	{
		//**
		//only for boundary nodes!!!
		if (newmode)
		{
      if ((i != RefNode1() || !reduced1) && (i != RefNode2() || !reduced2))
	      return Vector2D(XGP(nodes(i)->Get(1))-XGP(nodes(i)->Get(2)+1), XGP(nodes(i)->Get(1)+1)+XGP(nodes(i)->Get(2)));
			else
				return Vector2D(XGP(nodes(i)->Get(1)),XGP(nodes(i)->Get(1)+1));
		}
		else
			return Vector2D(XGP(nodes(i)->Get(1)),XGP(nodes(i)->Get(2)));
	}

	virtual Vector2D GetNodePos2DD(int i) const 
	{
		//**
		//only for boundary nodes!!!
		if (newmode)
		{
      if ((i != RefNode1() || !reduced1) && (i != RefNode2() || !reduced2))
	      return Vector2D(nodes(i)->Pos().X()+XGD(nodes(i)->Get(1))-XGD(nodes(i)->Get(2)+1),
					nodes(i)->Pos().Y()+XGD(nodes(i)->Get(1)+1)+XGD(nodes(i)->Get(2)));
			else
			  return Vector2D(nodes(i)->Pos().X()+XGD(nodes(i)->Get(1)), nodes(i)->Pos().Y()+XGD(nodes(i)->Get(1)+1));
		}
		else
			return Vector2D(nodes(i)->Pos().X()+XGD(nodes(i)->Get(1)),nodes(i)->Pos().Y()+XGD(nodes(i)->Get(2)));
	}

	virtual Vector2D GetNodeVel2DD(int i) const 
	{
		//**
		//only for boundary nodes!!!
		if (newmode)
		{
      if ((i != RefNode1() || !reduced1) && (i != RefNode2() || !reduced2))
	      return Vector2D(XGPD(nodes(i)->Get(1))-XGPD(nodes(i)->Get(2)+1),XGPD(nodes(i)->Get(1)+1)+XGPD(nodes(i)->Get(2)));
			else
			  return Vector2D(XGPD(nodes(i)->Get(1)), XGPD(nodes(i)->Get(1)+1));
		}
		else
			return Vector2D(XGPD(nodes(i)->Get(1)),XGPD(nodes(i)->Get(2)));
	}

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi);
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);
	virtual void AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f);   // f += dpdq*lambda

	virtual void GetIntDuDq(Matrix& dudq) //in fact it is DuDq Transposed
	{
		GetH(dudq);
	}

	virtual void GetH(Matrix& H)
	{
		H = Hr;
	}

	virtual Box3D GetElementBox() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBox());
		}
		return box;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D box;
		for (int i=1; i <= NFFRFElements(); i++) 
		{
			box.Add(GetFFRFElement(i).GetBoundingBoxD());
		}
		return box;
	}

	virtual void DrawElement();

	virtual void PrintToAbaqus(ofstream& os);

private:

	//++++++++++++++++++++++
	//ACRS:
	int refnode1;
	int refnode2;

	//++++++++++++++++++++++
	//CMS:
	int nbmodes;
	int nimodes;
	int sosfull;
	int isCMS;

	int useSparseMK;
	int usesparsepre;
	//++++++++++++++++++++++
	int newmode;
	int reduced1;
	int reduced2;
	int usematlab; //use matlab for eigencomputation
	//++++++++++++++++++++++

	Matrix3D initrot;
	double initphi;

	Vector xg;

	Matrix Phi_CB; //transformation matrix reduced to full dof
	TArray<int> boundarydof;      //mark dof with 1 which are fixed
	TArray<int> boundarynode;     //node number which is fixed

	//Matrix M, K; //not used, only temporary
	Matrix Mr, Kr, Hr; //full and reduced system matrices
	SparseMatrix Msparse, Ksparse;
	double volume; //dor intdudq

	//double xactfullvalue; //dummy variable, for non-const call of XActFull()
	//double xactfullvalueD; //dummy variable, for non-const call of XActFull()

	///+++++++++++++++++++++
	//store matrices for FFRF:
	/*Matrix Sbar_tilde;
	Matrix SbarS;
	SparseMatrix Sbar_tildeSM;
	SparseMatrix SbarSM;
	Vector Ibar11S,Ibar21S,Ibar12S,Ibar22S;
	Vector I1S;
	double I1122S;*/
	//++++++++++++++++++++++
};





#endif

