//#**************************************************************
//# filename:             body2D.h
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
//   Body2D   Body2D   Body2D   Body2D   Body2D   Body2D   Body2D   Body2D  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class Body2D: public Element //$EDC$[beginclass,classname=Body2D,parentclassname=Element]
{
public:
	//Body2D():Element() {mbs = NULL;};
	Body2D(MBS* mbsi):Element(mbsi) //, pref3D(0.,0.,0.), rotref3D(1.,1.,1.) 
	{
		Body2D::ElementDefaultConstructorInitialization();
		//type = TBody;
	};
	Body2D(const Body2D& e):Element(e.mbs) {CopyFrom(e);};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Body2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Element::CopyFrom(e);
		const Body2D& ce = (const Body2D&)e;
		pref3D = ce.pref3D;
		rotref3D = ce.rotref3D;
		size = ce.size;
	}

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization()
	{
		type = TBody;
		pref3D = Vector3D(0.,0.,0.);
		rotref3D = Matrix3D(1.);
		size = Vector3D(0.,0.,0.);
	}

	virtual void Initialize() 
	{
		Element::Initialize();
	};

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for load/save/edit:
	virtual const char* GetElementSpec() const {return "Body2D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	virtual const Vector3D& GetSize() const {return size;}

	virtual Box3D GetElementBox() const	//$ DR 2013-09-23
	{
		return Box3D(ToP3D(GetPos2D(Vector2D(-0.5*GetSize().X())))+pref3D,ToP3D(GetPos2D(Vector2D(0.5*GetSize().X())))+pref3D); //$ DR 2013-09-23
	}

	virtual Box3D GetElementBoxD() const //$ DR 2013-09-23
	{
		return Box3D(ToP3D(GetPos2DD(Vector2D(-0.5*GetSize().X())))+pref3D,ToP3D(GetPos2DD(Vector2D(0.5*GetSize().X())))+pref3D); //$ DR 2013-09-23
	}

	virtual void EvalF(Vector& f, double t) {};  
	virtual void EvalG(Vector& f, double t) {}; 
	virtual void EvalM(Matrix& m, double t) {}; 
	virtual void EvalF2(Vector& f, double t) 
	{
		Element::EvalF2(f,t);
	}; 


	virtual int IS() const {return 0;};  //implicit (algebraic) size
	virtual int SOS() const {return 0;}; //size of second order equations, len(u)
	virtual int ES() const {return 0;};  //size of first order explicit equations
	virtual int SS() const {return 2*SOS()+ES()+IS();};  //system size

	//# Timeint specific derived functions: for discontinuities
	virtual double PostNewtonStep(double t) {return 0;};
	virtual void PostprocessingStep() {};
	virtual double GetError() const 
	{
		return Element::GetError();
	};

	virtual int Dim() const {return 2;} //default value
	virtual int IsRigid() const {return 0;} //default value

	//get angle in 2D
	virtual double GetAngle2D() const {return XG(3);};
	virtual double GetAngle2D(const Vector2D& p_loc) const {return XG(3);};
	virtual double GetAngle2DP(const Vector2D& p_loc) const {return XGP(3);};

	virtual double GetCurvature2D(const Vector2D& p_loc) const {return 0.;};
	virtual double GetCurvature2DP(const Vector2D& p_loc) const {return 0.;};
	virtual double GetCurvature2DD(const Vector2D& p_loc) const {return 0.;};


	virtual Vector2D GetRefPos2D() const {return Vector2D(XG(1),XG(2));};
	//virtual Vector2D GetRefPos2DD() const {return Vector2D(XGD(1),XGD(2));};
	virtual Vector2D GetRefPos2DD() const 
	{
		if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.)
		{
			Vector2D p = Vector2D(XGD(1),XGD(2));
			Vector2D p0 = GetRefPosInit2D();
			return (p - p0)*GetMBS()->GetDOption(105) + p0;
		}
		else
		{
			return Vector2D(XGD(1),XGD(2));
		}
	};
	virtual Vector2D GetRefVel2D() const {return Vector2D(XGP(1),XGP(2));};
	virtual Vector2D GetRefVel2DD() const {return Vector2D(XGPD(1),XGPD(2));};


	virtual Vector3D GetRefPos() const
	{
		Vector2D p = GetRefPos2D();
		return ToP3D(Vector3D(p.X(), p.Y(), 0.));
	}
	virtual Vector3D GetRefPosD() const
	{
		Vector2D p = GetRefPos2DD();
		return ToP3D(Vector3D(p.X(), p.Y(), 0.));
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//transform 2D into 3D points for drawing and constraints:
	virtual Vector3D GetPos(const Vector3D& p_loc) const
	{
		Vector2D p = GetPos2D(Vector2D(p_loc.X(), p_loc.Y()));
		return ToP3D(Vector3D(p.X(), p.Y(), p_loc.Z()));
	}
	virtual Vector3D GetPosD(const Vector3D& p_loc) const
	{
		Vector2D p = GetPos2DD(Vector2D(p_loc.X(), p_loc.Y()));
		return ToP3D(Vector3D(p.X(), p.Y(), p_loc.Z()));
	}
	virtual Vector3D GetVel(const Vector3D& p_loc) const
	{
		Vector2D p = GetVel2D(Vector2D(p_loc.X(), p_loc.Y()));
		return ToV3D(p);
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
	//if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
	//{
	//	Matrix3D A = GetRotMatrix2DD();
	//	Matrix3D A0= GetRotMatrixInit2D();
	//	double fact = GetMBS()->GetDOption(105);

	//	return (fact*A-(fact-1.)*A0)*p_loc+GetRefPos2DD();
	//}
	//else
	//{
		return GetRotMatrix2DD()*p_loc+GetRefPos2DD();
	//}
	};

	//get angle in 2D
	virtual double GetAngle2DD() const {return XGD(3);};

	//for loads:
	virtual void GetIntDuxDq(Vector& dudq) {};
	virtual void GetIntDuyDq(Vector& dudq) {};
	virtual void GetIntDuzDq(Vector& dudq) {};
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0) {UO() << "ERROR:called Body2D::EmptyFunction\n";};
	virtual void ApplyDprefdq(Vector& f, const Vector2D& x) {};
	virtual void ApplyDrotrefdq(Vector& f, const double& x) {};

	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& dpdqi) {assert(0);};

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi) {assert(0);};
	virtual void GetdRotvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d) {assert(0);}; 
	virtual void GetdPosdx(const Vector2D& ploc, Vector2D& dpdx) {assert(0);}; 
	virtual void GetdAngle2Ddx(const Vector2D& ploc, double& dphidx) {assert(0);}; 

	virtual void GetNodedPosdqT(int node, Matrix& dpdqi) {assert(0);}; //get matrix with derivatives
	virtual void AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f) {assert(0);};   // f += dpdq*lambda

	virtual Matrix3D GetRotMatrix2D() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(XG(3)); 
		double sinphi = sin(XG(3)); 
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
		double cosphi = cos(XG(3)); 
		double sinphi = sin(XG(3)); 
		double phip = XGP(3); 
		rot(1,1) =-phip*sinphi;
		rot(1,2) =-phip*cosphi;
		rot(2,1) = phip*cosphi;
		rot(2,2) =-phip*sinphi;

		return rot;
	}
	//$ SW 2013-10-18: added function GetRotMatrix2DPP
	virtual Matrix2D GetRotMatrix2DPP() const
	{
		Matrix2D rot;
    rot.SetSize(Dim(),Dim());
		double cosphi = cos(XG(3)); 
		double sinphi = sin(XG(3)); 
		double phip = XGP(3); 
		double phip2 = phip*phip;
		double phipp = XGPP(3);

		rot(1,1) =-phipp*sinphi-phip2*cosphi;
		rot(1,2) =-phipp*cosphi+phip2*sinphi;
		rot(2,1) = phipp*cosphi-phip2*sinphi;
		rot(2,2) =-phipp*sinphi-phip2*cosphi;
		return rot;
	}

	virtual Matrix3D GetRotMatrix() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double cosphi = cos(XG(3)); 
		double sinphi = sin(XG(3)); 
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
		double phi  = XG(3); 
		double phip = XGP(3); 
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

	//drawing matrices:
	virtual Matrix3D GetRotMatrix2DD() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());

		double phi = XGD(3);
		if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
		{
			double fact = GetMBS()->GetDOption(105);
			phi = fact*XGD(3) - (fact-1.)*x_init(3);
		}

		rot(1,1) = cos(phi);
		rot(1,2) =-sin(phi);
		rot(2,1) = sin(phi);
		rot(2,2) = cos(phi);

		return rot;
	}
	////get initial rotation matrix; for drawing and computation:
	//virtual Matrix3D GetRotMatrixInit2D() const
	//{
	//	//assumes ordering of DOFS: x,y,phi
	//	Matrix3D rot;
 //   rot.SetSize(Dim(),Dim());
	//	double phi = x_init(3);
	//	rot(1,1) = cos(phi);
	//	rot(1,2) =-sin(phi);
	//	rot(2,1) = sin(phi);
	//	rot(2,2) = cos(phi);

	//	return rot;
	//}
	//virtual Matrix3D GetRotMatrixInit() const
	//{
	//	//assumes ordering of DOFS: x,y,phi
	//	Matrix3D rot(0.);
	//	double phi = x_init(3);
	//	rot(1,1) = cos(phi);
	//	rot(1,2) =-sin(phi);
	//	rot(2,1) = sin(phi);
	//	rot(2,2) = cos(phi);

	//	return rot;
	//}
	//	
	virtual Matrix3D GetRotMatrixD() const
	{
		//assumes ordering of DOFS: x,y,phi
		Matrix3D rot;
		double phi = XGD(3); 

		if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
		{
			double fact = GetMBS()->GetDOption(105);
			phi = fact*XGD(3) - (fact-1.)*x_init(3);
		}

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
		double phi  = XGD(3); 
		double phip = XGPD(3); 
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

	//transform 2D points
	virtual Vector3D ToP3D(const Vector2D& p) const 
	{
		return Vector3D(p.X(),p.Y(),0)*rotref3D+pref3D;
	}
	//set offset for 3D drawing and computation
	virtual void SetPRef3D(const Vector3D& p) {pref3D = p;}
	//set rotation for 3D transformation
	virtual void SetRotRef3D(const Matrix3D& rot) {rotref3D = rot;}

	//for velocities
	virtual Vector3D ToV3D(const Vector2D& v) const 
	{
		return Vector3D(v.X(),v.Y(),0)*rotref3D;
	}
	//transform 2D points (including virtual z-ccordinate), makes no sense for velocities
	virtual Vector3D ToP3D(const Vector3D& p) const 
	{
		return Vector3D(p.X(),p.Y(),p.Z())*rotref3D+pref3D;
	}

	virtual void DrawElement() 
	{
		Element::DrawElement();
	};

protected:
	Vector3D pref3D;   //$EDC$[varaccess,EDCvarname="reference_position",EDCfolder="Geometry",tooltiptext="Reference point for transformation of planar objects to 3D; p = [X, Y, Z]"]
	Matrix3D rotref3D; //$EDC$[varaccess,EDCvarname="rotation_matrix",EDCfolder="Geometry",tooltiptext="Rotation matrix for transformation of planar objects to 3D"]
	Vector3D size;     //$EDC$[varaccess,EDCvarname="general_size",EDCfolder="Geometry",tooltiptext="Dimensions of a regular cube [L_x, L_y, L_z]"]
	//Vector3D pref3D;
	//Matrix3D rotref3D;
	//Vector3D size;
}; //$EDC$[endclass,Body2D]