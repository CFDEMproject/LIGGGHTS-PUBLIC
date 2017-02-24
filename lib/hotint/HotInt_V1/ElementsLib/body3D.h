//#**************************************************************
//# filename:             body3D.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
//# comments:
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

#include "element.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Body3D   Body3D   Body3D   Body3D   Body3D   Body3D   Body3D   Body3D   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const int show_orientation=0;
const int max_rigid3D_coordinates = 7; //maximum coordinates in any rigid body; this is to avoid static vectors instead of ConstVector and ConstMatrix

class Body3D: public Element //$EDC$[beginclass,classname=Body3D,parentclassname=Element]
{
public:
	//Body3D():Element() {mbs = NULL;};
	Body3D(MBS* mbsi):Element(mbsi) {type = TBody;};
	Body3D(const Body3D& e):Element(e.mbs) {CopyFrom(e);};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Body3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Element::CopyFrom(e);
		const Body3D& ce = (const Body3D&)e;
		size = ce.size;
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void Initialize() 
	{
		Element::Initialize();
	};
	virtual const Vector3D& GetSize() const {return size;}
	
	virtual void EvalF(Vector& f, double t) {};  
	virtual void EvalG(Vector& f, double t) {}; 
	virtual void EvalM(Matrix& m, double t) {}; 
	virtual void EvalF2(Vector& f, double t) 
	{
		Element::EvalF2(f,t);
	}; 
	virtual int FastStiffnessMatrix() const {return 0;}
	virtual void StiffnessMatrix(Matrix& m) {}; //fill in sos x sos components, m might be larger

	virtual int IS() const {return 0;};  //implicit (algebraic) size
	virtual int SOS() const {return 0;}; //size of second order equations, len(u)
	virtual int ES() const {return 0;};  //size of first order explicit equations

	virtual int Dim() const {return 3;} //default value
	virtual int IsRigid() const {return 0;} //default value
	virtual int ElementBandwidth() const {return SS();}

	//# Timeint specific derived functions: for discontinuities
	virtual double PostNewtonStep(double t) {return 0;};
	virtual void PostprocessingStep() {};
	virtual double GetError() const {return Element::GetError();};

	virtual Vector3D GetRefPosD() const;
	//get the rotated and translated position of a local point at the body
	//virtual Vector3D GetPos_dc(const Vector3D& ploc, TComputeDrawInitFlag flag) const;


	virtual Vector3D GetPos(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacement(const Vector3D& p_loc) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const;
	virtual double GetAngleAroundAxis(const Vector3D& p_loc, const Vector3D& axis) const;
	virtual double GetAngleAroundLocalAxis(const Vector3D& p_loc, const Vector3D& axis) const;
	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const 
	{
		mbs->UO() << "Error: Angular velocity for element type of element " << GetOwnNum() << " not yet implemented!\n";
		return Vector3D(0);
	}

	//for loads:
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x) {UO() << "ERROR:called Body3D::EmptyFunction\n";};
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x) {UO() << "ERROR:called Body3D::EmptyFunction\n";};
	virtual void GetIntDuDq(Matrix& dudq) {UO() << "ERROR:called Body3D::EmptyFunction\n";}; //in fact it is DuDq Transposed
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0) {UO() << "ERROR:called Body3D::EmptyFunction\n";};

	virtual void GetDuxDq(Vector& dudq) {UO() << "ERROR:called Body3D::EmptyFunction\n";};
	virtual void GetDuyDq(Vector& dudq) {UO() << "ERROR:called Body3D::EmptyFunction\n";};
	virtual void GetDuzDq(Vector& dudq) {UO() << "ERROR:called Body3D::EmptyFunction\n";};
	//return the derivative of a position (local at ploc) with respect to all coordinates q, for ext. force
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdPosdqT\n";};
	virtual void GetdAngVeldqpT(const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdAngVeldqT\n";};
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi) {assert(0);}; //get matrix with derivatives
	virtual void AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f) {assert(0);};   // f += dpdq*lambda

	//return the derivative of the global roations (x/y/z) at ploc with respect to all coordinates q, for ext. moment
	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdRotdqT\n";};
	//return the derivative d(Rot*vloc)/dq
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdRotvdqT\n";};
	//return the derivative dP(x,y,z)/dx
	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) {UO() << "ERROR:called Body3D::GetdPosdx\n";};

	//rot matrix A describes the transformation of the local(body) to the
	//global coordinate system, such that r_g=R_g+A_bg*u_b
	virtual Matrix3D GetRotMatrix() const {return Matrix3D(1);}
	virtual Matrix3D GetRotMatrixP() const {return Matrix3D();}
	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const {return GetRotMatrix();}
	virtual Matrix3D GetRotMatrixP(const Vector3D& ploc) const {return GetRotMatrixP();}
	virtual Matrix3D GetRotMatrixD(const Vector3D& ploc) const {return GetRotMatrixD();}

	virtual Vector3D GetRotMatv(const Vector3D& ploc, const Vector3D& v) {return GetRotMatrix(ploc)*v;}
	virtual Vector3D GetRotMatPv(const Vector3D& ploc, const Vector3D& v) {return GetRotMatrixP(ploc)*v;}

	//functions for drawing:
	virtual Vector3D GetPosD(const Vector3D& p_loc) const;
	virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const;
	virtual Vector3D GetVelD(const Vector3D& p_loc) const;

	//drawing matrices:
	virtual Matrix3D GetRotMatrixD() const {return Matrix3D();}
	virtual Matrix3D GetRotMatrixPD() const {return Matrix3D();}

	virtual Box3D GetElementBox() const
	{
		Box3D b;
		b.Add(GetPos(Vector3D(-0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPos(Vector3D( 0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPos(Vector3D(-0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPos(Vector3D( 0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPos(Vector3D(-0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPos(Vector3D( 0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPos(Vector3D(-0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPos(Vector3D( 0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));
		return b;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b;
		b.Add(GetPosD(Vector3D(-0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPosD(Vector3D( 0.5*size.X(),-0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPosD(Vector3D(-0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPosD(Vector3D( 0.5*size.X(), 0.5*size.Y(),-0.5*size.Z())));
		b.Add(GetPosD(Vector3D(-0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPosD(Vector3D( 0.5*size.X(),-0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPosD(Vector3D(-0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));
		b.Add(GetPosD(Vector3D( 0.5*size.X(), 0.5*size.Y(), 0.5*size.Z())));
		return b;
	}

	virtual void DrawElement();

protected:
	Vector3D size; 
}; //$EDC$[endclass,Body3D]