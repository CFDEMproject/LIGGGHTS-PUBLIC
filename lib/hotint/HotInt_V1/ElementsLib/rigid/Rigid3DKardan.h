//#**************************************************************
//#
//# filename:             Rigid3DKardan.h
//#
//# author:               Ludwig Rafael, Gerstmayr Johannes, Yury Vetyukov
//#
//# generated:						10.5.2010
//# description:          3D Element Library - rigid body with Kardan angles
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
 
#ifndef RIGID3DKARDAN__H
#define RIGID3DKARDAN__H
//rigid cube using Kardan angles
class Rigid3DKardan: public Rigid3D //$EDC$[beginclass,classname=Rigid3DKardan,parentclassname=Rigid3D,addelementtype=TAEBody,addelementtypename=Rigid3DKardan,
//texdescription="A rigid body in 3D, implemented with bryant angles.",texdescriptionDOF="The first 3 degrees of freedom are those describing the position. The rotation is parameterized with 3 bryant angles.",
//texdescriptionGeometry="The center of gravity, S, is defined by the vector initial\_position, which is in global coordinates, see figure \ref{Rigid3Dfigure2}. The rotation of the body-fixed local coordinate system w.r.t. the global coordiante system is defined by the Matrix initial\_rotation. \\In order to define the position of a point P of the element, e.g. for connectors or sensors, the local coordinate system is used. The reference point is the center of mass, S, so the values of the local coordinates can be positive or negative.",
//example="Rigid3DKardan.txt",figure="Rigid3D,Rigid3DKardan"]
{
public:
	// these are the types of angles:
	// axis of the first rotation, axis of the second rotation, and of the third one.
	enum RotationsSequence
	{
		xyz,	// Cardan angles for x axis
		zxy,	// Cardan angles for z axis - default
		zxz		// Euler angles
	};
	// here we can set the default rotations sequence,
	// which will be used if the actual one is not provided in the constructor
	static void SetDefaultRotationsSequence(RotationsSequence rs) { defaultRotationsSequence = rs; }

protected:
	RotationsSequence rotationsSequence;
	static RotationsSequence defaultRotationsSequence;

public:
	//Body3D():Element() {mbs = NULL;};
	Rigid3DKardan(MBS* mbsi);
	Rigid3DKardan(const Rigid3DKardan& e):Rigid3D(e.mbs) {CopyFrom(e);};

	//x0 - initial position and velocity,
	// phi0 - initial Euler angles and vector of angular velocity
	Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
		const Vector3D& si, const Vector3D& coli, RotationsSequence rotSeq);

	Rigid3DKardan(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi,
		const Vector3D& si, const Vector3D& coli);

	virtual void ElementDefaultConstructorInitialization(); //$ DR 2012-12: ElementDefaultConstructorInitialization added

	//compute initial vector from Euler angles:
	void ComputeInitialConditions(const Vector3D& xp, const Vector3D& vp, 
		const Vector3D& phi, const Vector3D& phip, Vector& xinit);

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);
	virtual Element* GetCopy();

	virtual const char* GetElementSpec() const {return "Rigid3DKardan";}

	virtual void GetElementData(ElementDataContainer& edc);
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter


	virtual const Matrix3D& GetRotInertia() const {return Iphi;}
// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}
	
	virtual void EvalF(Vector& f, double t) {}  

	virtual void EvalG(Vector& f, double t) {}
	virtual void EvalF2(Vector& f, double t); 

	virtual int SOS() const {return 6;}; //size of second order equations, len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	virtual int SOSowned_RS() const {return 6;}; //for resort, the algebraic element variables belong to SOS2
	virtual int IS_RS() const {return 0;};  //algebraic (implicit) size


	virtual int IsRigid() const {return 1;} //default value

//deprecated
//	virtual void GetBetaP(double& betap0, double& betap1, double& betap2) const; //GetBetaP(Vector& betap)
//	virtual void GetBeta(double& beta0, double& beta1, double& beta2) const;
//	virtual void GetBetaD(double& beta0, double& beta1, double& beta2) const;
//	virtual void GetBetaPD(double& betap0, double& betap1, double& betap2) const;
	//new
	virtual int GetIndBeta(int i) const {return 3+i;}
	virtual double GetLagrangeMultEP() const {assert(0); return XG(13);}//

	virtual Matrix3D GetG_Kardan(double beta0, double beta1, double beta2)	const;
  	virtual Matrix3D GetGT_Kardan(double beta0, double beta1, double beta2) const;
	virtual Matrix3D GetG() const;
	virtual Matrix3D GetG(Vector& beta) const;
	virtual Matrix3D GetGT() const;
	virtual Matrix3D GetGbar() const;
	virtual Matrix3D GetGbarT() const;
	virtual Matrix3D GetGbarp() const;
	virtual Matrix3D GetGbarpT() const;

	// derivatives of the matrix Gbar with respect to the rotational degrees of freedom;
	// GetDGbarDBeta3 === 0
	Matrix3D GetDGbarDBeta1() const;
	Matrix3D GetDGbarDBeta2() const;

	// Compute d omega / d theta = [ d Gbar/dtheta1 thetap , d Gbar/dtheta2 thetap , 0 ]
	virtual Matrix3D GetDOmegaDTheta(Vector& betap) const;

	//virtual void Rigid3DKardan::PostprocessingStep(); // test function
	//rot matrix A describes the transformation of the local(body) to the
	//global coordinate system, such that r_g=R_g+A_bg*u_b
	virtual int NRotParam() const { return 3; }

	virtual Matrix3D ComputeRotMatrix(const Vector& beta) const;

	virtual Matrix3D ComputeRotMatrixP(const Vector& beta, const Vector& betap) const;

	virtual Vector3D GetPosD(const Vector3D& p_loc) const;

// under construction - create a matlab/simulink c-file for this element
	virtual void MatlabExportF2();//(Vector& f, double t);

protected:
	//mechanical:
	//Matrix3D Iphi;

	//Vector3D betap; //do not use for graphics!!!!
}; //$EDC$[endclass,Rigid3DKardan]
#endif
