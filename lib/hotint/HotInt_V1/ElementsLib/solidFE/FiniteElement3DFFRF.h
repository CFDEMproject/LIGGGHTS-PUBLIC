//#**************************************************************
//#
//# filename:             FiniteElement3DFFRF.h
//#
//# author:               Gerstmayr Johannes, Sinwell Astrid, YV
//#
//# generated:						October 2010
//# description:          functionality of 3D finite elements related to FFRF & CMS
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
 

#pragma once

#include "FiniteElement3D.h"
#include "FE3DHexTet.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    class FFRFData
///  stores necessary information for floating frame of reference formulation elements
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// needs to be generalized to 2D case
class FFRFData
{
public:

	FFRFData() : mbs(NULL), CMSelnum(0), I1S(0), SbarS(0)
	{ ; }

	FFRFData(MBS* mbsi) : mbs(mbsi), CMSelnum(0), I1S(0), SbarS(0)
	{ ; }

	FFRFData(const FFRFData& ff)
	{
		mbs = ff.mbs;
		CMSelnum = ff.CMSelnum;
		I1S = ff.I1S;
		SbarS = ff.SbarS;
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				IklS[i][j] = ff.IklS[i][j];
				IbarklS[i][j] = ff.IbarklS[i][j];
				SbarklS[i][j] = ff.SbarklS[i][j];
			}
			SbarTilde[i] = ff.SbarTilde[i];
		}
	}

	void Initialize()
	{ ; }	

	virtual Rigid3D& GetReferenceFrame() {return dynamic_cast<Rigid3D&>(GetMBS()->GetElement(CMSelnum));}
	virtual const Rigid3D& GetReferenceFrame() const {return dynamic_cast<const Rigid3D&>(GetMBS()->GetElement(CMSelnum));}


	MBS* GetMBS() {return mbs;}
	const MBS* GetMBS() const {return mbs;}

	MBS* mbs;
	int CMSelnum; //CMSelnum!=0: element number of Reference frame; CMSelnum==0: classical finite element
	//int isCMS; //1...is modally reduced element in CMS (AuxElement), 0...classical FFRF element
	Vector I1S;
	Matrix SbarS;
	double IklS[3][3];
	Vector IbarklS[3][3];
	Matrix SbarklS[3][3]; //reduced matrix! 
	Matrix SbarTilde[3]; //SbarTilde23, SbarTilde31, SbarTilde12 = //Sbar23-Sbar32, Sbar31-Sbar13, Sbar12-Sbar21
};


/// this class implements particular functionality of FFRF and CMS elements
class FiniteElement3DFFRF : public FiniteElement3D
{
public:
	// default constructor
	FiniteElement3DFFRF(MBS * mbs) : FiniteElement3D(mbs), ffrfdata(NULL) {}
	~FiniteElement3DFFRF()
	{
		if(ffrfdata != NULL)
			delete ffrfdata;
	}

	/// Filling in the data from an existing element
	virtual void CopyFrom(const Element& e);

	// late initialization
	virtual void Initialize();
	virtual void LinkToElements(); ///< For a given mesh topology, computes the mapping between the local and the global degrees of freedom

	virtual void SetFFRFElement(int CMSelementI);

	virtual int IsFFRF() const { return 1; } // {return ffrfdata!=0;}
	virtual int IsCMS() const { return type & TCMS; }
	virtual int IsGCMS() const { return type & TGCMS; }
	virtual int FastStiffnessMatrix() const { return 0; }

	// initialize FFRFData
	virtual void SetFFRF(int isffrf);
	// underlying reference frame
	const Rigid3D& ReferenceFrame() const { return ffrfdata->GetReferenceFrame(); }
	Rigid3D& ReferenceFrame() { return ffrfdata->GetReferenceFrame(); }
	// shortcuts to the properties of the reference frame
	virtual Vector3D RefFramePos() const {return ReferenceFrame().GetRefPos();}
	virtual Vector3D RefFrameVel() const {return ReferenceFrame().GetRefVel();}
	virtual Vector3D RefFramePosD() const {return ReferenceFrame().GetRefPosD();}
	virtual Vector3D RefFramePosInit() const {return ReferenceFrame().GetRefPosInit();}
	virtual Vector3D RefFrameVelD() const {return ReferenceFrame().GetRefVelD();}
	virtual Matrix3D RefFrameRot() const {return ReferenceFrame().GetRotMatrix();}
	virtual Matrix3D RefFrameRotP() const {return ReferenceFrame().GetRotMatrixP();}
	virtual Matrix3D RefFrameRotD() const {return ReferenceFrame().GetRotMatrixD();}
	virtual Matrix3D RefFrameRotPD() const {return ReferenceFrame().GetRotMatrixPD();}

	virtual const Element& GetElement(int globind) const 
	{
		if (!IsCMS())
			return GetMBS()->GetElement(globind);
		else
			return GetMBS()->GetAuxElement(globind);
	}
	virtual Element* GetElementPtr(int globind)
	{
		if (!IsCMS())
			return GetMBS()->GetElementPtr(globind);
		else
			return GetMBS()->GetAuxElementPtr(globind);
	}
	virtual const Node& GetNode(int i) const 
	{
		if (!IsCMS())
			return GetMBS()->GetNode(NodeNum(i));
		else
			return ReferenceFrame().GetNode(NodeNum(i));
	}
	virtual Node& GetNode(int i) 
	{
		if (!IsCMS())
			return GetMBS()->GetNode(NodeNum(i));
		else
			return ReferenceFrame().GetNode(NodeNum(i));
	}

	virtual int SOS() const { if (IsGCMS()) {return FlexDOF();} else {return FlexDOF()+FFRFDim();} } ///<Size of K and M (stiffness matrix and mass matrix)

	// number of degrees of freedom of the Floating Frame
	virtual int FFRFDim() const { if (IsGCMS()) {return 0;} else {return ReferenceFrame().SOSRigid();} }
	// number of translational degrees of freedom for the Floating Frame
	virtual int FFRFTranslationDim() const { if (IsGCMS()) {return 0;} else {return Dim();} }
	// number of rotational degrees of freedom for the Floating Frame
	virtual int FFRFRotationDim() const { if (IsGCMS()) {return 0;} else {return ReferenceFrame().NRotParam();} }
	// size of K and M for not reduced system, number of total degrees of freedom
	virtual int XGLength() const
	{
		if (IsCMS())
			return FlexDOF();
		else
			return ReferenceFrame().SOSRigid() + FlexDOF();
	}

	virtual const double& GetXact(int i) const;
	virtual const double& GetDrawValue(int iloc) const;

	// AP: Was macht diese Funktion??
	virtual int PerformNodeCheck() const { return IsCMS()==0; }

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);

	// precomputed matrices - inertial properties and coupling with the rigid body degrees of freedom of the reference frame
	virtual void GetI1(Vector& I1); //Shabana p. 209-211
	virtual void GetSbar(Matrix& Sbar);
	virtual double GetIkl(int k, int l);
	virtual void GetIbarkl(int k, int l, Vector& I1);
	virtual void GetSbarkl(int k, int l, Matrix& Sbar);
	// useful for FFRF
	virtual double GetIntRhoUkUl(int k, int l, const Vector& xgloc);
	virtual void GetIntRhoUkUlMat(Matrix3D& mat, const Vector& xgloc);
	virtual double GetIntRhoUkUlP(int k, int l, const Vector& xgloc, const Vector& xglocp);

	// quadratic velocity vector - geometric stiffening effect
	virtual void AddQuadraticVelocityVector(Vector& fadd, double t);
  virtual void GetIbarThetaF(Matrix& Ibar_theta_f, const Vector& xgloc);
  virtual void GetIbarThetaTheta(Matrix3D& Ibar_theta_theta, const Vector& xgloc);
  virtual void GetIbarThetaThetaP(Matrix3D& Ibar_theta_theta, const Vector& xgloc, const Vector& xgploc);

	// absolute nodal positions
	virtual Vector3D GetNodePos(int i) const;
	virtual Vector3D GetNodePosD(int i) const;
	virtual Vector3D GetNodeVel(int i) const;

	// absolute positions inside the volume 
	// FFRF/CMS: with the accounting of the motion of the reference frame
	// GCMS: using the absolute coordinates as in standard FiniteElement3D
	virtual Vector3D GetPos(const Vector3D& p_loc) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const;
	virtual Vector3D GetPosD(const Vector3D& p_loc, int use_magnification) const;
	virtual Vector3D GetVelD(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacement(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const;

	virtual void EvalM(Matrix& m, double t);

	// for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq); //in fact it is DuDq Transposed

	virtual int AddBodyNode(Node & n);

protected:
	FFRFData* ffrfdata;
};

typedef HexahedralGeneric<FiniteElement3DFFRF> HexahedralFFRF;
typedef TetrahedralGeneric<FiniteElement3DFFRF> TetrahedralFFRF;