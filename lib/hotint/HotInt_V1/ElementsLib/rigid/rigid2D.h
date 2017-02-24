//#**************************************************************
//#
//# filename:             element2D.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
 
#ifndef ELEMENT2D__H
#define ELEMENT2D__H

#include "element.h"
#include "body2d.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Mass2D  Mass2D  Mass2D  Mass2D  Mass2D  Mass2D  Mass2D  Mass2D  Mass2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class Mass2D: public Body2D //$EDC$[beginclass,classname=Mass2D,parentclassname=Body2D,addelementtype=TAEBody+TAE2D,addelementtypename=Mass2D,
//texdescription="A point mass in two dimensions with 2 position coordinates. The computation of the dynamics of the point mass is extremely simple, thus the Mass2D can be used for many body simulations (e.g. particles).",
//texdescriptionDOF="2 degrees of freedom: the position in 2 coordinates",
//texdescriptionLimitations="The mass has no rotations, thus external moments can not be applied. The transformation of local to global coordinates is based on a translation, i.e., the global mass position is added to the local coordinates.", 
//texdescriptionEquations="
//\begin{equation}
//m \ddot{\boldsymbol{x}} = \boldsymbol{F}	
//\end{equation}",
//example="mass2D.txt",
//figure="Mass2D"]
{
public:
	//Mass2D():Element() {mbs = NULL;};
	Mass2D(MBS* mbsi):Body2D(mbsi) {ElementDefaultConstructorInitialization();};
	Mass2D(const Mass2D& e):Body2D(e.mbs) {CopyFrom(e);};
	Mass2D(MBS* mbsi, const Vector& x0, double radi, double massi, const Vector3D& coli):
	  Body2D(mbsi)
	{
		ElementDefaultConstructorInitialization();
		mass=massi;
		size.X()=radi;
		size.Y()=radi;
		size.Z()=radi;
		//drawres = 5;

		x_init = x0;
		//rho = mass / (4./3. * MY_PI * Cub(radi)); 
		col = coli;
		//elementname = GetElementSpec();
	};

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization();		//$ DR 2012-07: moved to cpp-file

	//$ DR 2012-07: CheckConsistency added
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Mass2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const Mass2D& ce = (const Mass2D&)e;
		drawres = ce.drawres;
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for load/save/edit:
	virtual const char* GetElementSpec() const {return "Mass2D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const // $ MSax 2013-08-05 : added
	{
		TKinematicsAccessFunctions kaf = Body2D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_position_2D+TKAF_velocity_2D+TKAF_D_pos_D_q_2D);
	}

	virtual void Initialize() 
	{
		Body2D::Initialize();
	};

	virtual void SetDrawRes(int res) {drawres = res;}

	virtual const Vector3D& GetSize() const {return size;}
	
	virtual void EvalM(Matrix& m, double t) 
	{
		m(1,1) = mass; m(2,2) = mass;
	}; 
	virtual void EvalF2(Vector& f, double t) 
	{
		Body2D::EvalF2(f,t);

		if (GetMassDamping() != 0)
		{
			f(1) -= GetMassDamping()*mass*XGP(1);
			f(2) -= GetMassDamping()*mass*XGP(2);
		}

	}; 
	//virtual void EvalMinvF2(Vector& f, double t) 
	//{
	//	EvalF2(f, t);

	//	double minv = 1;
	//	if (mass != 0) minv = 1./mass;

	//	f *= minv;
	//}; 

	virtual Matrix3D GetRotMatrix() const // $ MSax 2013-08-08 : added
	{
		return Matrix3D(1.);
	}

	virtual Matrix3D GetRotMatrix2D() const // $ MSax 2013-08-08 : added
	{
		Matrix3D rot;
    rot.SetSize(Dim(),Dim());
		rot(1,1) = 1.;
		rot(1,2) = 0;
		rot(2,1) = 0;
		rot(2,2) = 1.;

		return rot;
	}

	virtual Matrix3D GetRotMatrixD() const // $ MSax 2013-08-08 : added
	{
		return Matrix3D(1.);
	}

	virtual Matrix3D GetRotMatrixP() const // $ MSax 2013-08-09 : added
	{
		return Matrix3D(0.);
	}

	virtual int SOS() const {return 2;}; //size of second order equations, len(u)

	virtual int Dim() const {return 2;} //default value
	virtual int IsRigid() const {return 1;} //default value

	virtual const double& Rho() const {assert(0); double dummy = -1.; return dummy; /*return rho;*/} //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	virtual double& Rho() {assert(0); double dummy = -1.; return dummy; /*return rho;*/} //$ DR 2013-02-04 deleted rho from class element, do not use it here!

	//$ SW 2013-10-18: added function GetAcceleration
	//compute acceleration vector (after Timestep is successfully computed)
	virtual Vector2D GetAcceleration(const Vector2D& ploc=Vector2D(0.)) const
	{
		//the only point in a Mass2D is the center of gravity
		assert(ploc.X() == 0. && ploc.Y() == 0); 
		return Vector2D(XGPP(1),XGPP(2));
	}

	//get the rotated and translated position of a local point at the body
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const
  {
		return GetRefPos2D()+p_loc;
	};
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const
  {
		return GetRefVel2D();
	};

	//for loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector2D& x)
	{
		f(1) = x(1);
		f(2) = x(2);
	}
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector2D& x)
	{
		f(1) = 0;
		f(2) = 0;
	}
	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetDuxDq(Vector& dudq)
	{
		dudq.SetAll(0);
		dudq(1) = 1;
	}
	virtual void GetDuyDq(Vector& dudq)
	{
		dudq.SetAll(0);
		dudq(2) = 1;
	}
	//return the derivative of a position (local at ploc) with respect to all coordinates q
	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d)
	{
		d.SetSize(2,2);
		d(1,1)=1; d(1,2)=0;
		d(2,1)=0; d(2,2)=1;
	}
	//return the derivative d(Rot*vloc)/dq
	virtual void GetdRotvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d) 
	{
		d.SetSize(2,2);
		d.SetAll(0);
	};
	//return the derivative dP(x,y,z)/dx
	virtual void GetdPosdx(const Vector2D& ploc, Vector2D& dpdx) 
	{
		dpdx = Vector2D(1,0);
	};

	//functions for drawing:
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const
  {
		return GetRefPos2DD()+p_loc;
	};
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const
  {
		return GetRefVel2DD();
	};

	virtual void DrawElement();
	
	//for visualization and sensoring of displacements and velocities:
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);

protected:
	int drawres; //$EDC$[varaccess,EDCvarname="drawing_tiling",EDCfolder="Graphics",tooltiptext="tiling of circle/sphere to represent Mass2D;
	//the drawing_tiling should be set small in order to improve efficiency,
	//but large for nice graphical representations"]
	
  //EDC double size.X() //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="drawing radius of mass"] //DO NOT MODIFY THIS COMMENT!!!!

	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_position",EDCfolder="Initialization",vecstart=1,vecend=2,tooltiptext="initial values for position [x,y]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="initial_velocity",EDCfolder="Initialization",vecstart=3,vecend=4,tooltiptext="initial values for velocity [vx,vy]"]

	//EDC double mass; //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="total mass of point mass"] 

	//$ DR 2013-01-11 remove unused options:
	//EDC Vector3D size;     //$EDC$[varaccess,EDCvarname="general_size",EDCfolder="Geometry",remove]


}; //$EDC$[endclass,Mass2D]


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   Rigid2D   Rigid2D   Rigid2D   Rigid2D   Rigid2D   Rigid2D   Rigid2D   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// coords(sos=3): x,y,phi
class Rigid2D: public Body2D //$EDC$[beginclass,classname=Rigid2D,parentclassname=Body2D,addelementtype=TAEBody+TAE2D,addelementtypename=Rigid2D,
//texdescription="A rigid body in 2D.",texdescriptionDOF="The first 2 degrees of freedom are those describing the position in the xy-plane. The rotation around the local z-axis is parameterized with the third degree of freedom.",
//texdescriptionGeometry="The center of gravity, S, is defined by the vector initial\_position, which is in global coordinates. The rotation of the body-fixed local coordinate system w.r.t. the global coordiante system is defined by the variable initial\_rotation. \\In order to define the position of a point P of the element, e.g. for connectors or sensors, the local coordinate system is used. The reference point is the center of mass, S, so the values of the local coordinates can be positive or negative.",
//example="Rigid2D.txt",figure="Rigid2D,Rigid2D"]
{
public:
	Rigid2D(MBS* mbsi):Body2D(mbsi) {ElementDefaultConstructorInitialization(); };
	Rigid2D(MBS* mbsi, const Vector& x0, double rhoi, double Vi, double Ip,
		const Vector3D& si, const Vector3D& coli):
	  Body2D(mbsi)
	{
		SetRigid2D(x0,rhoi,Vi,Ip,si,coli);
	};

	Rigid2D(MBS* mbsi, const Vector& x0, double rhoi,
		const Vector3D& si, const Vector3D& coli):
	  Body2D(mbsi)
	{
		SetRigid2D(x0,rhoi,si,coli);
	};

	void SetRigid2D(const Vector& x0, double rhoi,
		const Vector3D& si, const Vector3D& coli)
	{
		double vol = si.X()*si.Y()*si.Z();
		SetRigid2D(x0,rhoi,vol,1./12.*vol*rhoi*(Sqr(si.X())+Sqr(si.Y())),si,coli);
	};

	void SetRigid2D(const Vector& x0, double rhoi, double Vi, double Ip,
		const Vector3D& si, const Vector3D& coli)
	{
		Iphi = Ip;
		size = si;
		mass = rhoi*Vi;
		x_init = x0;
		rho = rhoi;
		col = coli;
		//mbs->uout << "Rigid2D: m=" << mass << ", I=" << Iphi << "\n";
		elementname = GetElementSpec();
	};
	Rigid2D(const Rigid2D& e):Body2D(e.mbs) 
	{
		ElementDefaultConstructorInitialization(); 
		CopyFrom(e);
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Rigid2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const Rigid2D& ce = (const Rigid2D&)e;
		Iphi = ce.Iphi;
		// size = ce.size; //$RE 2013-08-20: already done in Body2D::CopyFrom(e) --> commented out
		rho = ce.rho; //DR 2013-02-04 deleted rho from class element
	}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for load/save/edit:
	virtual const char* GetElementSpec() const {return "Rigid2D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual void ElementDefaultConstructorInitialization();
 	virtual const double& Rho() const {return rho;}
	virtual double& Rho() {return rho;}

	virtual int IsRigid() const {return 1;} //default value

	virtual void Initialize() 
	{
		Body2D::Initialize();
	};
	
	//evaluate mass matrix
	virtual void EvalM(Matrix& m, double t) 
	{
		//m.FillWithZeros();
		m(1,1) = mass;
		m(2,2) = mass;
		m(3,3) = Iphi;
	}; 
	//second order equations: M \ddot u = F2(u,\dot u,t), len(u)
	//constraints set in Element!!!
	virtual void EvalF2(Vector& f, double t) 
	{
		Body2D::EvalF2(f,t);

		if (damping_m)
		{
			f(1) -= damping_m*mass*XGP(1);
			f(2) -= damping_m*mass*XGP(2);
			f(3) -= damping_m*Iphi*XGP(3);
		}
	}; 
	
	virtual int IS() const {return 0;};  //implicit (algebraic) size
	virtual int SOS() const {return 3;}; //size of second order equations, len(u)
	virtual int ES() const {return 0;};  //size of first order explicit equations
	virtual int SS() const {return 2*SOS()+ES()+IS();};  //system size


	//redundant with Body2D::GetPos()????? ++++++++++++++++++++++++++++
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
	virtual Vector3D GetVelD(const Vector3D& p_loc) const
	{
		Vector2D p = GetVel2DD(Vector2D(p_loc.X(), p_loc.Y()));
		return ToV3D(p);
	}

	//$ SW 2013-10-18: Added function GetAcceleration
	//compute acceleration vector (after Timestep is successfully computed)
	virtual Vector2D GetAcceleration(const Vector2D& ploc) const
	{
		Vector2D Rpp(XGPP(1),XGPP(2));
		Matrix2D App = GetRotMatrix2DPP();
		
		return Rpp + App*ploc;
	}

	//2013-21-08:[ similar to Rigid3D
	//virtual Vector3D GetAngularVel(const Vector3D& p_loc) const  //$ DR 2013-08-26 removed, does not make sense
	//{		
		//return GetRotMatrix()*GetAngularVelLocal(); 
	//}
	//Vector3D GetAngularVelLocal() const		//$ DR 2013-08-26 removed, does not make sense
	//{
	//	return Vector3D(0., 0., XG(6));
	//}
    //2013-21-08:] 
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		return GetRefPosD();
	}

	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		if (idof == 1) return ToP3D(Vector2D(1.,0.));
		else if (idof == 2) return ToP3D(Vector2D(0.,1.));
		else if (idof == 3) return ToP3D(Vector3D(0.,0.,2.));
		return Vector3D(0.,0.,0.);
	}

	//virtual Box3D GetElementBox() const	//$ DR 2013-09-23 moved to body2d
	//{
		//return Box3D(ToP3D(GetPos2D(Vector2D(-0.5*GetSize().X(),0))),
			//ToP3D(GetPos2D(Vector2D(0.5*GetSize().X(),0))));
	//}

	//virtual Box3D GetElementBoxD() const	//$ DR 2013-09-23 moved to body2d
	//{
		//return Box3D(ToP3D(GetPos2DD(Vector2D(-0.5*GetSize().X(),0))),
			//ToP3D(GetPos2DD(Vector2D(0.5*GetSize().X(),0)))+Vector3D(0,0,size.Z()));
	//}

	//---------------------------------------------------------------
	//for visualization of displacements and velocities:
	//---------------------------------------------------------------
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);
	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		TKinematicsAccessFunctions kaf = Body2D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_position_2D+TKAF_velocity_2D+TKAF_D_pos_D_q_2D); // $ MSax 2013-08-26 : added dposdq2d
	}
	virtual void DrawElement();
	

		//for body loads:
	virtual void GetDuxDq(Vector& dudq)
	{
		dudq(1) = 1; dudq(2) = 0; dudq(3) = 0;
	}
	virtual void GetDuyDq(Vector& dudq)
	{
		dudq(1) = 0; dudq(2) = 1; dudq(3) = 0;
	}
	virtual void ApplyDprefdq(Vector& f, const Vector2D& x) 
	{
		f(1) = x.X();
		f(2) = x.Y();
		f(3) = 0;
	};
	virtual void ApplyDrotrefdq(Vector& f, const double& x) 
	{
		f(1) = 0; f(2) = 0;
		f(3) = x;
	};

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi);
	virtual void GetdRotvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& dpdqi);
	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& dpdqi);


protected:
	//mechanical:

	//double Iphi;  principal inertia; automatic access not possible because used in special case only
	double rho;  //DR 2013-02-04 deleted rho from class element
	
	//RE 2013-08:[ access for Get/SetElementDataAuto
	double Iphi; //$EDC$[varaccess,EDCvarname="moment_of_inertia",EDCfolder="Physics",tooltiptext="[I_ZZ]"]	
	//EDC double mass;								//$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="mass of the body in kg"]
	//EDC Vector3D size;							//$EDC$[varaccess,EDCvarname="body_dimensions",EDCfolder="Graphics",tooltiptext="Dimensions of a regular cube [L_x, L_y, (L_z)]"]
	
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=1, vecend=2,EDCvarname="initial_position",EDCfolder="Initialization",tooltiptext="[X, Y]"]
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=4, vecend=5,EDCvarname="initial_velocity",EDCfolder="Initialization",tooltiptext="[vX, vY]"]
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=3, vecend=3,EDCvarname="initial_rotation",EDCfolder="Initialization",tooltiptext="rotation in rad"]
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=6, vecend=6,EDCvarname="initial_angular_velocity",EDCfolder="Initialization",tooltiptext="Angular velocity in rad/s"]
	
	// remove unused options:
	//EDC Vector3D size;							//$EDC$[varaccess,remove,EDCvarname="general_size",EDCfolder="Geometry"]
	//RE 2013-08:]

};//$EDC$[endclass,Rigid2D]





#endif