//#**************************************************************
//#
//# filename:             CMSElement.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						11. July 2007
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
 
#ifndef CMSElement3D__H
#define CMSElement3D__H



//this element contains the whole structure of the body (finite element numbers, nodes) in this class.
//all nodes have local indices in the CMSElement
template<class RIGID>
class CMSElement: public BaseCMSElement<RIGID>
{
public:

	// Constructor, TCMSflag is set
	CMSElement(MBS* mbsi):BaseCMSElement<RIGID>(mbsi)
	{
		type = (TMBSElement)(type|TCMSflag);
	};

	// Copy constructor
	CMSElement(const CMSElement& e):BaseCMSElement<RIGID>(e.mbs)
		{CopyFrom(e);};

	// Constructor,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	CMSElement(MBS* mbsi, const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli) : BaseCMSElement<RIGID>(mbsi)
	{
		type = (TMBSElement)(type|TCMSflag);
		solverparameters.DefaultInitialize();
		EVfile = mystr("");
		SetCMSElement(p, v, phi, phip, nimodesi, sizei, coli);
	};

	// Set-Function,
	// p ..... initial position
	// v ..... initial velocity
	// phi ... initial angle
	// phip .. initial angular velocity
	// nimodesi .. number of internal modes
	// sizei ..... size of frame (for drawing)
	// coli ...... color (for drawing, currently not used)
	void SetCMSElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
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
	void SetCMSElementKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs);


	// destructor
	~CMSElement()
	{
	}


	virtual Element* GetCopy()
	{
		Element* ec = new CMSElement(*this);
		return ec;
	}
	
	virtual void CopyFrom(const Element& e)
	{
		BaseCMSElement<RIGID>::CopyFrom(e);
		const CMSElement& ce = (const CMSElement&)e;
	}

	// initialization of  FFRF matrices is done here
	virtual void Initialize();
	virtual int PerformNodeCheck() const {return 0;}
	virtual const char* GetElementSpec() const {return "CMSElement";}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//must be called before Assemble!!!!!
	//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, transformation matrix, modal matrices
	virtual void DoModalAnalysis(const TArray<int2>& fixednodes);
	virtual void InitializeFFRFMatrices();
	//final matrices of FFRF-formulation:
	// compute mass matrix of CMS element (not constant!)
	virtual void EvalM(Matrix& m, double t); 
	// stiffness+quadratic velocity vector:
	virtual void EvalF2(Vector& f, double t);
	// quadratic velocity vector
	virtual void AddQuadraticVelocityVector(Vector& f, double t, Vector& temp);
	// for quadratic velocity vector
	virtual void ComputeIbarThetaTheta(Matrix3D& Ibartt, Matrix3D& uklbar);
	virtual void ComputeIbarThetaThetaP(Matrix3D& Ibarttp, Matrix3D& uklbar, Matrix3D& uklbarp);
	virtual void ComputeIbarThetaF(Vector& IbarThetaF0, Vector& IbarThetaF1, Vector& IbarThetaF2, Vector& xg, Vector& temp);
		

	// can be set to 2 if the stiffness matrix routine works fine (quadratic velocity vector missing currently)
	virtual int FastStiffnessMatrix() const;
	// returns Kr
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components only of stiffness matrix, m might be larger

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// FFRF matricestemplate <class RIGID>
	// ltg-map from sosfull x nmodes to sos_red x nmodes
	virtual void GetLocalMatrix(const Matrix& globalM, const TArray<int>& ltg, int sos, Matrix& localM)
	{
		localM.SetSize(sos, globalM.Getcols());
		for (int j=1; j<=sos; j++)
		{
			for (int k=1; k<=globalM.Getcols(); k++)
			{
				localM(j,k) = globalM(ltg.Get(j),k);
			}
		}
	}
	virtual void GetI1(Vector& I1); //Shabana p. 209-211
	virtual double GetIkl(int k, int l);
	virtual void GetIbarkl(int k, int l, Vector& I1);
	virtual void GetSbar(Matrix& Sbar);
	virtual void GetSbarkl(int k, int l, Matrix& Sbar);

	// scalar Integral  Int( rho u_k u_l ) dV
	virtual double GetIntRhoUkUl(int k, int l, const Vector& xgloc, Vector& temp, int use_mthetatheta_full);
	// FFRF: compute integrals d/dt int_V \rho ubar_k * ubar_l dV
	// attention: d/dt int_V \rho ubar_k * ubar_l dV = uklbarp(k,l) + uklbarp(l,k)
	virtual double GetIntRhoUkUlP(int k, int l, const Vector& xgloc, const Vector& xgploc, Vector& temp);
	// Vector3D S_tbar = I1 + Sbar * qf
	virtual void GetStbar(Vector3D& S_tbar, const Vector& xg, int use_Stbar_full=1)
	{
		for (int i = 1; i <= Dim(); i++)
		{
			S_tbar(i) = I1S(i);
			if (use_Stbar_full)
			{
				for (int j = 1; j <= NModes(); j++)
				{
					S_tbar(i) += SbarS(i,j)*xg(j); //this term is important for jacobian
				}
			}
		}
	}



	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual int SOS() const {return NModes()+SOSRigid();}; //size of second order equations, len(u), number of modes + 7 rigid body dofs
	virtual int SOSowned() const {return SOS();}; //size of second order dofs, len(u), which are not taken from external mbs nodes (all dofs in this case)
	virtual int SOSowned_RS() const {return SOS()-SOSRigid();}; //size of second order equations for resorting???
	virtual int IS_RS() const {return SOSRigid();}; // implicit size  for resorting???


	// number of flexible degrees of freedom
	virtual int FlexDOF() const {return NModes();}
	virtual int NModes() const {return NBModes()+NIModes();}

	virtual int IsCMS() const {return 1;}
	virtual int IsGCMS() const {return 0;}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position / velocity of node in the CMS element
	virtual Vector3D GetNodePos(int i) const ;
	virtual Vector3D GetNodeLocPos(int i) const ;
	virtual Vector3D GetNodeVel(int i) const ;
	virtual Vector3D GetNodePosD(int i) const;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// d = d pos / dq(ploc)
	// GetdPosdqt: only for rigid body position of rigid frame!!
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d); 
	// GetdPosdqt: only for rigid body rotation of rigid frame!!
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);
	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	// differentiation of global position of Node node with respect to REDUCED degrees of freedom!
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi);
	virtual void AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f);   // f += dpdq*lambda

	// int du/dq dV for the full cms element, nonlinear, depends on xg
	virtual void GetIntDuDq(Matrix& dudq); //in fact it is DuDq Transposed
	// int rho du/dq dV for the full cms element, nonlinear, depends on xg
	virtual void GetIntRhoDuDq(Matrix& rhodudq); //in fact it is DuDq Transposed
	// Phi_CB^T H, where H is the assembled H-matrix of the FFRFelements
	virtual void GetH(Matrix& H)
	{
		H = Hr;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//overwrite Rigid3D access functions for Euler parameters--> new ordering, Euler parameters are after mode-coordinates!!!
	virtual int GetIndRefPos(int i) const {return NModes()+i;}
	virtual int GetIndBeta(int i) const {return NModes()+3+i;}
	virtual Vector3D GetRefPos() const 
	{
		return Vector3D(XG(GetIndRefPos(1)),XG(GetIndRefPos(2)),XG(GetIndRefPos(3)));
	}
	virtual Vector3D GetRefVel() const 
	{
		return Vector3D(XGP(GetIndRefPos(1)),XGP(GetIndRefPos(2)),XGP(GetIndRefPos(3)));
	}
	virtual Vector3D GetRefPosD() const 
	{
		return Vector3D(XGD(GetIndRefPos(1)),XGD(GetIndRefPos(2)),XGD(GetIndRefPos(3)));
	}
	virtual Vector3D GetRefVelD() const 
	{
		return Vector3D(XGPD(GetIndRefPos(1)),XGPD(GetIndRefPos(2)),XGPD(GetIndRefPos(3)));
	}
	virtual Vector3D GetRefPosInit() const 
	{
		return Vector3D(x_init(GetIndRefPos(1)),x_init(GetIndRefPos(2)),x_init(GetIndRefPos(3)));
	}

	////the following overwritten functions lead to correct GetGbar, GetRotMatrix and similar !!!!
	//// should be replaced by functions with Vector& beta below (cmp. Rigid3D)
	//virtual void GetBetaP(double& betap0, double& betap1, double& betap2, double& betap3) const;
	//virtual void GetBeta(double& beta0, double& beta1, double& beta2, double& beta3) const;
	//virtual void GetBetaD(double& beta0, double& beta1, double& beta2, double& beta3) const;
	//virtual void GetBetaPD(double& betap0, double& betap1, double& betap2, double& betap3) const;
	//
	//new functions, should be used and tested
	virtual void GetBetaP(Vector& betap) const;
	virtual void GetBeta(Vector& beta) const;
	virtual void GetBetaD(Vector& beta) const;
	virtual void GetBetaPD(Vector& betap) const;


	virtual double GetLagrangeMultEP() const 
	{
		int ind = (NModes()+7)*2+1;
		return XG(ind);
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// AH: get angular momementum with respect to a reference point p_ref
	virtual Vector3D GetAngularMomentum(const Vector3D& p_ref) const;
	virtual double GetKineticEnergy();
	virtual double GetPotentialEnergy();

	virtual void SetDiagonalInternalDamping(const Vector& DrI) {Dr = DrI;}
	virtual const Vector& GetDiagonalInternalDamping() {return Dr;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer


	// under construction - create a matlab/simulink c-file for this element
	virtual void MatlabExportF2(Vector& f, double t);

private:

	Vector xg;

	///+++++++++++++++++++++
	//store matrices for FFRF, these are already the reduced matrices!!! --> DENSE Matrices
	Vector I1S;
	Matrix SbarS;
	double IklS[3][3];
	Vector IbarklS[3][3];
	Matrix SbarklS[3][3]; //reduced matrix! 
	//++++++++++++++++++++++
};


template <class RIGID>
class CMSSurface
{
public:
	CMSSurface(MBS* mbsi, CMSElement<RIGID>* cmsi)
	{
		mbs = mbsi;
		cms = cmsi;
	}

	virtual MBS* GetMBS() {return mbs;}
	virtual int AddNode(int nodenum, const Vector3D& locnode)
	{
		locnodes.Add(locnode);
		nodenums.Add(nodenum);
		return locnodes.Length();
	}
	virtual int NNodes() const {return locnodes.Length();}

	virtual const Vector3D& NodePos(int i) const {return locnodes(i);}

	//compute rigid body surface constraint transformation and return 1 if success and 0 if failed
	//rb_modenum specifies where to insert in transformation matrix PhiRB
	virtual int ComputeRigidBodyModes(int rb_modenum, Vector3D& midpoint, SparseMatrix& modes)
	{

	}

	//compute rigidbody_DOF
	//return 0 if not possible
	virtual int SelectRigidBodyDOF()
	{
		if (NNodes() < 3) return 0;

		double maxd = 1e30;
		Vector3D p(maxd, maxd, maxd);
		int np1 = 0;

		for (int i=1; i <= NNodes(); i++)
		{
			Vector3D pi = NodePos(i);
			if (pi.X() <= p.X() && pi.Y() <= p.Y() && pi.Z() <= p.Z()) 
			{
				p = pi;
				np1 = i;
			}
		}

		Vector3D p1 = p;
		p=Vector3D(-maxd, -maxd, -maxd);
		int np2 = 0;

		for (int i=1; i <= NNodes(); i++)
		{
			Vector3D pi = NodePos(i);
			if (pi.X() >= p.X() && pi.Y() >= p.Y() && pi.Z() >= p.Z()) 
			{
				p = pi;
				np2 = i;
			}
		}

		Vector3D p2 = p;
		double maxnorm = 0;
		Vector3D v0 = p1-p2;
		int np3 = 0;

		for (int i=1; i <= NNodes(); i++)
		{
			Vector3D pi = NodePos(i);
			Vector3D vn = pi-p1;
			if ((v0.Cross(vn)).Norm() > maxnorm) 
			{
				p = pi;
				np3 = i;
			}
		}

		GetMBS()->UO() << "nodes=" << np1 << "," << np2 << "," << np3 << "\n";

		



	}

private:
	TArray<Vector3D> locnodes;	//local node position of surface
	TArray<int> nodenums;				//node numbers in CMS
	CMSElement<RIGID>* cms;						//link to CMS element
	MBS* mbs;										//link to MBS system
	int rigidbody_DOF[6];				//selected rigid body DOF

	Vector3D midpoint;					//mid or reference point of surface
};


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   NodalConstraintCMS
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Constraint for a Node of an CMS element
//leads to linear constraint equations
template <class RIGID>
class NodalConstraintCMS: public NodalConstraint
{
public:
	// Constraint needs multibody system MBS and CMS element!
	NodalConstraintCMS(MBS* mbsi, int cms_element):NodalConstraint(mbsi)
	{	
		cms_el = cms_element;
	};

	// Copy constructor
	NodalConstraintCMS(const NodalConstraintCMS& ct): NodalConstraint(ct.mbs)
	{
		CopyFrom(ct);
	};

	virtual Element* GetCopy()
	{
		Element* ec = new NodalConstraintCMS(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)
	{
		NodalConstraint::CopyFrom(e);
		const NodalConstraintCMS& ce = (const NodalConstraintCMS&)e;
		cms_el = ce.cms_el;
	}

	virtual const char* GetElementSpec() const {return "NodalConstraintCMS";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//evaluate constraints: len(z)
	virtual void EvalG(Vector& f, double t);
	virtual void EvalF2(Vector& f, double t);

	virtual BaseCMSElement<RIGID>& GetCMSElement() { return dynamic_cast<BaseCMSElement<RIGID>&> (mbs->GetElement(cms_el)); }
	virtual const BaseCMSElement<RIGID>& GetCMSElement() const { return dynamic_cast<const BaseCMSElement<RIGID>&> (mbs->GetElement(cms_el)); }

	virtual Node& GetNode(int i) {return GetCMSElement().GetNode(NodeNum(i));} //get local node number i
	virtual const Node& GetNode(int i) const {return GetCMSElement().GetNode(NodeNum(i));} //get local node number i

	virtual int SOS() const {return GetCMSElement().SOS();};  // explicit size, number degrees of freedom of the CMS element, since all modes of CMS element can affect node position
	virtual int SOSowned() const {return 0;}  // number of explicit equations added by element
	virtual int IS() const { if (UsePenaltyFormulation()) return 0; else return loccoords.Length();};  

	virtual int Dim() const {return 3;}  
	// one node which belongs to CMS element
	virtual int NCMSNodes() const {return 1; };
	// zero MBS-nodes
	virtual int NNodes() const {return 0;}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return GetCMSElement().GetNodePos(NodeNum(1));
	}
	
	virtual void DrawElement();
	virtual void LinkToElements();

	// do not change global init vector due to CMS constraint
	virtual void SetGlobalInitConditions(Vector& x_glob)
	{
		;
	}
	virtual int IsLocalNodeConstraint() const {return 0;}

protected:
	// number of CMS element in MBS
	int cms_el;
};


#endif

