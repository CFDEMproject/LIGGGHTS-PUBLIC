//#**************************************************************
//#
//# filename:             element3D.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          3D Element Library
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
 
#ifndef ELEMENT3D__H
#define ELEMENT3D__H

#include "body3d.h"
#include "material.h"
#include "Node.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Mass3D  Mass3D  Mass3D  Mass3D  Mass3D  Mass3D  Mass3D  Mass3D  Mass3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Mass3D: public Body3D //$EDC$[beginclass,classname=Mass3D,parentclassname=Body3D,addelementtype=TAEBody,addelementtypename=Mass3D,texdescription="A point mass in three dimensions with 3 position coordinates. The computation of the dynamics of the point mass is extremely simple, thus the Mass3D can be used for many body simulations (e.g. particles).",
//texdescriptionDOF="3 degrees of freedom: the position in 3 coordinates",
//texdescriptionLimitations="The mass has no rotations, thus external moments can not be applied. The transformation of local to global coordinates is based on a translation, e.g. the global mass position is added to the local coordinates.",
//example="AddElement.txt"]
{
public:
	//Mass3D():Element() {mbs = NULL;};
	Mass3D(MBS* mbsi):Body3D(mbsi) {ElementDefaultConstructorInitialization();}; //$ DR 2012-07: ElementDefaultConstructorInitialization for CEDCParser
	Mass3D(const Mass3D& e):Body3D(e.mbs) {CopyFrom(e);};
	Mass3D(MBS* mbsi, const Vector& x0, double radi, double massi, const Vector3D& coli):
	  Body3D(mbsi)
	{
		ElementDefaultConstructorInitialization(); //$ DR 2012-07: ElementDefaultConstructorInitialization for CEDCParser
		mass=massi;
		size.X()=radi;
		size.Y()=radi;
		size.Z()=radi;
		//drawres = 6; //$ DR 2012-07: moved to ElementDefaultConstructorInitialization

		x_init = x0;
		double rho = mass / (4./3. * MY_PI * Cub(radi));
		col = coli;


		Material mat(GetMBS());
		mat.SetMaterialRigid(rho);

		int material_num = GetMBS()->AddMaterial(mat);
		SetMaterialNum(material_num);

		//elementname = GetElementSpec(); //$ DR 2012-07: moved to ElementDefaultConstructorInitialization
	};

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization();		//$ DR 2012-07: ElementDefaultConstructorInitialization added

	//$ DR 2012-07: CheckConsistency added
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Mass3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body3D::CopyFrom(e);
		const Mass3D& ce = (const Mass3D&)e;
		drawres = ce.drawres;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	virtual const char* GetElementSpec() const {return "Mass3D";}

	virtual void Initialize() 
	{
		Body3D::Initialize();
	};

	virtual void SetDrawRes(int res) {drawres = res;}

	virtual void SetMass(double m) {mass = m;}	//$ DR 2013-02-05 added

	virtual const Vector3D& GetSize() const {return size;}
	
	virtual void EvalM(Matrix& m, double t) 
	{
		m(1,1) = mass; m(2,2) = mass; m(3,3) = mass; 
	}; 
	virtual void EvalF2(Vector& f, double t) 
	{
		Body3D::EvalF2(f,t);

		if (GetMassDamping() != 0)
		{
			f(1) -= GetMassDamping()*mass*XGP(1);
			f(2) -= GetMassDamping()*mass*XGP(2);
			f(3) -= GetMassDamping()*mass*XGP(3);
		}

	}; 
	virtual void EvalMinvF2(Vector& f, double t) 
	{
		EvalF2(f, t);

		double minv = 1;
		if (mass != 0) minv = 1./mass;

		f *= minv;
	}; 

	virtual int SOS() const {return 3;}; //size of K and M

	virtual int Dim() const {return 3;} //default value
	virtual int IsRigid() const {return 1;} //default value

	//get the rotated and translated position of a local point at the body
	virtual Vector3D GetPos(const Vector3D& p_loc) const
  {
		return GetRefPos()+p_loc;
	};
	virtual Vector3D GetVel(const Vector3D& p_loc) const
  {
		return GetRefVel();
	};

	virtual Vector3D GetRefConfPos(const Vector3D& ploc) const { return GetRefPosInit();}		//$ DR+PG 2013-05-15

	//compute acceleration vector (after Timestep is successfully computed)
	virtual Vector3D GetAcceleration(const Vector3D& ploc=Vector3D(0.)) const //$ RL 2011-3-14: added for use in MBSSensor.
	{
		                                                             // _ 
		assert(ploc.X() == 0. && ploc.Y() == 0. && ploc.Z() == 0. ); // u = 0 --> Term 2 == Term 3 == 0
		//Term 1  +    Term 2       + Term 3
		//..   ..      _    _   _        _      _    ..
		//ri = Ri + Ai[w x (w x 0) + Ai(alpha x 0) = Ri
		//return Vector3D(GetMBS()->GetAcceleration(LTG(1 + SOS())), GetMBS()->GetAcceleration(LTG(2 + SOS())), GetMBS()->GetAcceleration(LTG(3 + SOS()))); // d²/dt(R) // $ MSax 2013-07-16 : removed
		return Vector3D(XGPP(1),XGPP(2),XGPP(3)); // d²/dt(R) // $ MSax 2013-07-16 : added
	}

	//functions for drawing:
	virtual Vector3D GetPosD(const Vector3D& p_loc) const
  {
		return GetRefPosD()+p_loc;
	};
	virtual Vector3D GetVelD(const Vector3D& p_loc) const
  {
		return GetRefVelD();
	};

	//for loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x)
	{
		f(1) = x(1);
		f(2) = x(2);
		f(3) = x(3);
	}
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x)
	{
		f(1) = 0;
		f(2) = 0;
		f(3) = 0;
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
	virtual void GetDuzDq(Vector& dudq)
	{
		dudq.SetAll(0);
		dudq(3) = 1;
	}
	//return the derivative of a position (local at ploc) with respect to all coordinates q
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d)
	{
		d.SetSize(3,3);
		d(1,1)=1; d(1,2)=0; d(1,3)=0;
		d(2,1)=0; d(2,2)=1; d(2,3)=0;
		d(3,1)=0; d(3,2)=0; d(3,3)=1;
	}
	//return the derivative d(Rot*vloc)/dq
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d) 
	{
		d.SetSize(3,3);
		d.SetAll(0);
	};
	//return the derivative dP(x,y,z)/dx
	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) 
	{
		dpdx = Vector3D(1,0,0);
	};

	virtual void DrawElement();
	virtual void DrawElementFunc(const Vector3D& pos, const Vector3D& refpos); //this function does the drawing
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_ref_conf_position+TKAF_D_pos_D_x+TKAF_D_pos_D_q+TKAF_D_rot_v_D_q+TKAF_acceleration);
	}

protected:
	int drawres; //$EDC$[varaccess,EDCvarname="drawing_tiling",EDCfolder="Graphics",tooltiptext="tiling of circle/sphere to represent Sphere"]
  //EDC double size.X() //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="drawing radius of mass"] 
	//EDC Vector3D x_init //$EDC$[varaccess,EDCvarname="initial_position",EDCfolder="Initialization",vecstart=1,tooltiptext="coordinates for initial position of mass [X Y Z]"]
	//EDC Vector3D x_init //$EDC$[varaccess,EDCvarname="initial_velocity",EDCfolder="Initialization",vecstart=4,tooltiptext="coordinates for initial velocity of mass [X Y Z]"]

	//EDC double mass; //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="total mass of point mass"] 

}; //$EDC$[endclass,Mass3D]


class NodalMass3D: public Mass3D //$EDC$[beginclass,classname=NodalMass3D,parentclassname=Mass3D]
{
public:
	//Mass3D():Element() {mbs = NULL;};
	NodalMass3D(MBS* mbsi):Mass3D(mbsi) {};
	NodalMass3D(const NodalMass3D& e):Mass3D(e.mbs) {CopyFrom(e);};
	NodalMass3D(MBS* mbsi, int nodenumi, double radi, double massi, const Vector3D& coli):
	  Mass3D(mbsi)
	{
		nodenum = nodenumi;
		mass=massi;
		size.X()=radi;
		size.Y()=radi;
		size.Z()=radi;
		drawres = 6;

		x_init = GetNode(1).X_Init();
		if (x_init.Length() != 6) //if node has no initialization, assume zeros
		{
			x_init = Vector(6);
		}

		double rho = mass / (4./3. * MY_PI * Cub(radi));
		col = coli;

		Material mat(GetMBS());				// DR 2011-05-10: changed due to new constructor of material
		mat.SetMaterialRigid(rho);		// DR 2011-05-10: changed due to new constructor of material
		int material_num = GetMBS()->AddMaterial(mat);
		SetMaterialNum(material_num);

		elementname = GetElementSpec();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new NodalMass3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Mass3D::CopyFrom(e);
		const NodalMass3D& ce = (const NodalMass3D&)e;
		drawres = ce.drawres;
		nodenum = ce.nodenum;
	}

	virtual void ElementDefaultConstructorInitialization() // $ MSax 2013-02-19 : added
	{
		Mass3D::ElementDefaultConstructorInitialization();
		nodenum = 1;
		drawres = 6;
	}

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Mass3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf-TKAF_ref_conf_position+TKAF_node_ref_conf_position);
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	virtual const char* GetElementSpec() const {return "NodalMass3D";}

	virtual int SOS() const {return 3;}; //size of K and M
	virtual int SOSowned() const {return 0;}; //size of second order degrees of freedom added by element; DOF from node!

	virtual void LinkToElements(); //this function maps the degrees of freedom of the node to NodalMass3D
	virtual void Initialize();
	//virtual void SetGlobalInitConditions(Vector& x_glob) {}//do not add initial conditions!
	//virtual const Vector& GetXInit() const {return GetNode(1).X_Init();}

	//nodal functions:
	virtual int NNodes() const {return 1; };
	virtual const int& NodeNum(int i) const {return nodenum;}
	virtual int& NodeNum(int i) {return nodenum;}

	virtual const Node& GetNode(int i) const {return GetMBS()->GetNode(NodeNum(i));} //get local node number i
	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));} //get local node number i
	virtual Vector3D GetNodeLocPos(int i) const {return Vector3D(0.,0.,0.);} //local position of node

	virtual Vector3D GetNodePos(int i) const {return GetRefPos();}

	virtual Vector3D GetNodePosD(int i) const {return GetRefPosD();}

	virtual Vector3D GetNodeVel(int i) const {return GetRefVel();}

	virtual Vector3D GetNodeRefPos(int i) const
	{
#ifdef _DEBUG
		assert(i==1);
#endif
		return GetRefPos();
	}
	virtual Vector3D GetRefPos() const 
	{
		return GetMBS()->GetNode(NodeNum(1)).GetPos();
	}
	virtual Vector3D GetRefVel() const 
	{
		return GetMBS()->GetNode(NodeNum(1)).GetVel();
	}
	virtual Vector3D GetRefPosD() const 
	{
		return GetMBS()->GetNode(NodeNum(1)).GetPosD();
	}
	virtual Vector3D GetRefVelD() const 
	{
		return GetMBS()->GetNode(NodeNum(1)).GetVelD();
	}
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi) 
	{
		dpdqi.SetSize(3,3);
		dpdqi.FillWithZeros();

		dpdqi(1,1) = 1; //dnx/dnx
		dpdqi(2,2) = 1; //dny/dny
		dpdqi(3,3) = 1; //dny/dny
	};

	virtual void DrawElement();

protected:
	int nodenum; //$EDC$[varaccess,EDCvarname="node_number",EDCfolder="",tooltiptext="node number to which the mass refers"]
	//access functions, which are referring to Node entries, because Set/GetElementData calls parent function of "Body3D" and not of "Mass3D"
	//EDC int drawres; //!EDC$[varaccess,EDCvarname="drawing_tiling",EDCfolder="Graphics",tooltiptext="tiling of circle/sphere to represent Sphere"] // MSax 2013-02-20: deleted beacause already added in base class to edc
	//EDC double mass; //!EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="total mass of point mass"] //DO NOT MODIFY THIS COMMENT!!!! // MSax 2013-02-20: deleted beacause already added in base class to edc
  //EDC double size.X() //!EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="drawing radius of mass"] //DO NOT MODIFY THIS COMMENT!!!! // MSax 2013-02-20: deleted beacause already added in base class to edc
}; //$EDC$[endclass,NodalMass3D]


//MSax is responsible for NodalDiskMass3D
class NodalDiskMass3D: public NodalMass3D //$EDC$[beginclass,classname=NodalDiskMass3D,parentclassname=NodalMass3D,addelementtype=TAEBody,addelementtypename=NodalDiskMass3D,texdescription="This is a disk mass for the purpose of rotordynamics applications and should be used together with the RotorBeamXAxis element.",texdescriptionNode="The DOF of the disk element are stored in a node. To create a new disk element the user has to define a 'Node3DR123' node. This node type has 6 DOF. The first 3 DOF describe the node displacement ($x,y,z$) w.r.t local rotor element coordinate system, the last 3 DOF are angles of rotation ($\phi_x, \phi_y, \phi_z$) w.r.t local rotor element coordinate system. The rotation about the local x-axis is considered as large, the rotations about the local y and z-axes are considered as small (linearized angles).",example="NodalDiskMass3D.txt",figure="NodalDiskMass3D"]
{
public:
	NodalDiskMass3D(MBS* mbsi):NodalMass3D(mbsi) {ElementDefaultConstructorInitialization();}; // $ MSax 2013-04-17 : added ElementDefaultConstructorInitialization for CEDCParser
	NodalDiskMass3D(const NodalDiskMass3D& e):NodalMass3D(e.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(e);
	};
	NodalDiskMass3D(MBS* mbsi, int nodenumi, double radi, double thicknessi, double massi, const Vector3D& Ipi, const Vector3D& coli, int full_mass_matrixi = 1):
	  NodalMass3D(mbsi)
	{
		nodenum = nodenumi;
		mass=massi;
		Ip=Ipi;
		size.X()=thicknessi;
		size.Y()=radi;
		size.Z()=radi;
		drawres = 6;
		full_mass_matrix = full_mass_matrixi;

		x_init = GetNode(1).X_Init();
		if (x_init.Length() != 12) //if node has no initialization, assume zeros
		{
			x_init = Vector(12);
		}

		col = coli;

		elementname = GetElementSpec();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new NodalDiskMass3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		NodalMass3D::CopyFrom(e);
		const NodalDiskMass3D& ce = (const NodalDiskMass3D&)e;
		Ip = ce.Ip;
		full_mass_matrix = ce.full_mass_matrix;
	}

	virtual void ElementDefaultConstructorInitialization() // $ MSax 2013-02-19: added
	{
		NodalMass3D::ElementDefaultConstructorInitialization();
		x_init = Vector(SS()); //zero initialized
		Ip = Vector3D(1.);
		full_mass_matrix = 1;
	}
	virtual int IsRigid() const {return 0;} // $ MSax 2013-04-25 added because GetdPosdqT should be used by MBSLoad
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		return TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_angular_velocity+TKAF_D_pos_D_q+TKAF_D_rot_D_q);
	}

	virtual void EvalM(Matrix& m, double t) 
	{
		m(1,1) = mass; m(2,2) = mass; m(3,3) = mass; m(4,4) = Ip(1); m(5,5) = Ip(2); m(6,6) = Ip(3); 

		if (full_mass_matrix)
		{
			m(5,4) = XG(6)*Ip(1); 
			m(4,5) = m(5,4);
			m(6,4) = -1*XG(5)*Ip(1);
			m(4,6) = m(6,4);
		}
	}; 

	virtual void EvalF2(Vector& f, double t) 
	{
		Body3D::EvalF2(f,t);

		// gyro terms
		f(5) -= XGP(6)*Ip(1)*XGP(4);
		f(6) += XGP(5)*Ip(1)*XGP(4);
	}; 

	virtual void GyroscopicMatrix(SparseMatrix& locgy)const  // $ MSax 2013-07-25 : added
	{
		locgy.SetSize(SOS(),SOS());
		locgy.FillWithZeros();
		locgy(5,6) = -1.*Ip(1)*XGP(4);
		locgy(6,5) = Ip(1)*XGP(4);
	};

	virtual void EvalMinvF2(Vector& f, double t) 
	{
		Element::EvalMinvF2(f,t); // slow version
	}; 

	virtual void GetRotMatrixFromXGD(Matrix3D& rot) const
	{
		double theta = XGD(4); double phi_y = XGD(5); double phi_z = XGD(6);
		GetRotMatrixFromAngles(theta, phi_y, phi_z, rot);
	}

	virtual void GetRotMatrixFromXG(Matrix3D& rot) const
	{
		double theta = XG(4); double phi_y = XG(5); double phi_z = XG(6);
		GetRotMatrixFromAngles(theta, phi_y, phi_z, rot);
	}

	virtual void GetRotMatrixFromAngles(double theta, double phi_y, double phi_z, Matrix3D& rot) const
	{
		rot(1,1) = 1;        rot(1,2) = -1*cos(theta)*phi_z+sin(theta)*phi_y;       rot(1,3) = sin(theta)*phi_z+cos(theta)*phi_y;
		rot(2,1) = phi_z;    rot(2,2) = cos(theta);                                 rot(2,3) = -1*sin(theta);
		rot(3,1) = -1*phi_y; rot(3,2) = sin(theta);                                 rot(3,3) = cos(theta);
	}

	virtual void GetRotMatrixFromXGP(Matrix3D& rot) const // $ MSax 2013-07-11: added
	{
		double theta = XG(4); double phi_y = XG(5); double phi_z = XG(6); 		double thetaP = XGP(4); double phi_yP = XGP(5); double phi_zP = XGP(6);
		rot(1,1) = 0;						rot(1,2) = -1*(-sin(theta)*thetaP*phi_z+cos(theta)*phi_zP)+ cos(theta)*thetaP*phi_y+sin(theta)*phi_yP;		rot(1,3) = cos(theta)*thetaP*phi_z+sin(theta)*phi_zP-sin(theta)*thetaP*phi_y+cos(theta)*phi_yP;
		rot(2,1) = phi_zP;			rot(2,2) = -1*sin(theta)*thetaP;																																					rot(2,3) = -1*cos(theta)*thetaP;
		rot(3,1) = -1*phi_yP;		rot(3,2) = cos(theta)*thetaP;																																							rot(3,3) = -sin(theta)*thetaP;
	}

	virtual void GetRotMatrixFromXGPD(Matrix3D& rot) const //$ SW 2013-10-23: added
	{
		double theta = XGD(4); double phi_y = XGD(5); double phi_z = XGD(6); 		double thetaP = XGPD(4); double phi_yP = XGPD(5); double phi_zP = XGPD(6);
		rot(1,1) = 0;						rot(1,2) = -1*(-sin(theta)*thetaP*phi_z+cos(theta)*phi_zP)+ cos(theta)*thetaP*phi_y+sin(theta)*phi_yP;		rot(1,3) = cos(theta)*thetaP*phi_z+sin(theta)*phi_zP-sin(theta)*thetaP*phi_y+cos(theta)*phi_yP;
		rot(2,1) = phi_zP;			rot(2,2) = -1*sin(theta)*thetaP;																																					rot(2,3) = -1*cos(theta)*thetaP;
		rot(3,1) = -1*phi_yP;		rot(3,2) = cos(theta)*thetaP;																																							rot(3,3) = -sin(theta)*thetaP;
	}

	virtual void GetRotMatrixFromXGPP(Matrix3D& rot) const //$ SW 2013-10-21: added
	{
		double theta = XG(4);     double phi_y = XG(5);     double phi_z = XG(6);
		double thetaP = XGP(4);   double phi_yP = XGP(5);   double phi_zP = XGP(6);
		double thetaPP = XGPP(4); double phi_yPP = XGPP(5); double phi_zPP = XGPP(6);
		double thetaP2 = thetaP*thetaP; double phi_yP2 = phi_yP*phi_yP; double phi_zP2 = phi_zP*phi_zP;

		rot(1,1) = 0;						rot(1,2) = -1*(-cos(theta)*thetaP2*phi_z -sin(theta)*thetaPP*phi_z -sin(theta)*thetaP*phi_zP  -sin(theta)*thetaP*phi_zP + cos(theta)*phi_zPP)   - sin(theta)*thetaP2*phi_y + cos(theta)*thetaPP*phi_y + cos(theta)*thetaP*phi_yP + cos(theta)*thetaP*phi_yP + sin(theta)*phi_yPP;		rot(1,3) = -1*sin(theta)*thetaP2*phi_z +  cos(theta)*thetaPP*phi_z +  cos(theta)*thetaP*phi_zP +   cos(theta)*thetaP*phi_zP + sin(theta)*phi_zPP   - cos(theta)*thetaP2*phi_y - sin(theta)*thetaPP*phi_y - sin(theta)*thetaP*phi_yP  - sin(theta)*thetaP*phi_yP + cos(theta)*phi_yPP;
		rot(2,1) = phi_zPP;			rot(2,2) = -1*cos(theta)*thetaP2 - 1*sin(theta)*thetaPP;																									rot(2,3) = sin(theta)*thetaP2 - cos(theta)*thetaPP;
		rot(3,1) = -1*phi_yPP;	rot(3,2) = -sin(theta)*thetaP2 + cos(theta)*thetaPP;																											rot(3,3) = -1*cos(theta)*thetaP2 - sin(theta)*thetaPP;
	}

	virtual Matrix3D GetRotMatrixD(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRotMatrixFromXGD(rot);
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		return n.GetLocalFrame()*rot;
	}

	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRotMatrixFromXG(rot);
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		return n.GetLocalFrame()*rot;
	}

	virtual Matrix3D GetRotMatrixP(const Vector3D& ploc) const // $ MSax 2013-07-11: added
	{
		Matrix3D rot;
		GetRotMatrixFromXGP(rot);
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		return n.GetLocalFrame()*rot;
	}

	virtual Matrix3D GetRotMatrixPD(const Vector3D& ploc) const //$ SW 2013-10-23: added
	{
		Matrix3D rot;
		GetRotMatrixFromXGPD(rot);
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		return n.GetLocalFrame()*rot;
	}

	virtual Matrix3D GetRotMatrixPP(const Vector3D& ploc) const // $ SW 2013-10-23: added
	{
		Matrix3D rot;
		GetRotMatrixFromXGPP(rot);
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		return n.GetLocalFrame()*rot;
	}

	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const  // $ MSax 2013-07-11: added
	{
		Matrix3D omega_skew = (GetRotMatrixP(p_loc)*GetRotMatrix(p_loc).GetTp());
		return Vector3D(-omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2));
	}

	virtual Vector3D GetAngularVelD(const Vector3D& p_loc) const  //$ SW 2013-10-25: added
	{
		Matrix3D omega_skew = (GetRotMatrixPD(p_loc)*GetRotMatrixD(p_loc).GetTp());
		return Vector3D(-omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2));
	}

	//for loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x)
	{
		UO() << "Warning: Not implemented function ApplyDprefdq!\n";
		assert(0);
	}

	//return the derivative of a position (local at ploc) with respect to all coordinates q, needed if a point mass is added to the disk
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d)
	{
		d.SetSize(6,3);
		d.SetAll(0);

		Matrix tmp;

		GetdRotvdqTLoc(ploc, ploc, tmp);

		tmp(1,1)=1;
		tmp(2,2)=1;
		tmp(3,3)=1;

		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		Matrix A0T = (Matrix)n.GetLocalFrame().GetTp();

		Mult(tmp,A0T,d);
	}

	//for loads:
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x)
	{
		UO() << "Warning: Not implemented function ApplyDrotrefdq!\n";
		assert(0);
	}

	//return the derivative d(Rot*vloc)/dq
	virtual void GetdRotvdqTLoc(const Vector3D& vloc, const Vector3D& ploc, Matrix& d) 
	{
		d.SetSize(6,3);
		d.SetAll(0);

		d(4,1) = (XG(6)*sin(XG(4))+XG(5)*cos(XG(4)))*vloc.Y()+(XG(6)*cos(XG(4))-XG(5)*sin(XG(4)))*vloc.Z();
		d(4,2) = -sin(XG(4))*vloc.Y()-cos(XG(4))*vloc.Z();
		d(4,3) = cos(XG(4))*vloc.Y()-sin(XG(4))*vloc.Z();
		d(5,1) = sin(XG(4))*vloc.Y()+cos(XG(4))*vloc.Z();
		d(6,1) = -cos(XG(4))*vloc.Y()+sin(XG(4))*vloc.Z();
		d(5,3) = -vloc.X();
		d(6,2) = vloc.X();
	};

	//return the derivative d(Rot*vloc)/dq
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d) 
	{
		Matrix tmp;
		GetdRotvdqTLoc(vloc, ploc, tmp);

		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		Matrix A0T = (Matrix)n.GetLocalFrame().GetTp();

		Mult(tmp,A0T,d);
	};

	//return the derivative of the global roations (x/y/z) at ploc with respect to all coordinates q, for ext. moment
	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d) 
	{
		d.SetSize(6,3);
		d.SetAll(0);

		d(4,1) = 1;          d(4,2) = XG(6);    d(4,3) = -1*XG(5);
		d(5,1) = -1*XG(6);   d(5,2) = 1; 
		d(6,1) = XG(5);                         d(6,3) = 1;

		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		Matrix A0 = (Matrix)n.GetLocalFrame();
		TransformMatrix(A0,d);
	};

	virtual void TransformMatrix(Matrix3D rot,Matrix& d) //transform dudq, etc. matrices by means of rot0 reference rotation
	{
	//apply reference rotation:
	for (int i=1; i <= d.Getrows(); i++)
	{
		Vector3D v = rot*Vector3D(d(i,1),d(i,2),d(i,3));
		d(i,1) = v.X();
		d(i,2) = v.Y();
		d(i,3) = v.Z();
	} 
}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	virtual const char* GetElementSpec() const {return "NodalDiskMass3D";}

	virtual int SOS() const {return 6;}; //size of K and M
	virtual int SOSowned() const {return 0;}; //size of second order degrees of freedom added by element; DOF from node!

	virtual void LinkToElements(); //this function maps the degrees of freedom of the node to NodalDiskMass3D
	virtual void Initialize();

	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> &variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

	//$ SW 2013-10-21: added
	virtual Vector3D GetAcceleration(const Vector3D& ploc=Vector3D(0.)) const
	{
		Vector3D Rpp(XGPP(1),XGPP(2),XGPP(3));
		Matrix3D App = GetRotMatrixPP(ploc);
		return Rpp+App*ploc;
	}
	
	//$ SW 2013-10-22: added
	virtual Vector3D GetVel(const Vector3D& ploc=Vector3D(0.)) const
	{
		Vector3D Rp(XGP(1),XGP(2),XGP(3));
		Matrix3D Ap = GetRotMatrixP(ploc);
		return Rp+Ap*ploc;
	}

	//$ SW 2013-10-23: added
	virtual Vector3D GetVelD(const Vector3D& ploc=Vector3D(0.)) const
	{
		Vector3D RpD(XGPD(1),XGPD(2),XGPD(3));
		Matrix3D ApD = GetRotMatrixPD(ploc);
		return RpD+ApD*ploc;
	}

	virtual Vector3D GetPosD(const Vector3D& p_loc) const
  {
		if (GetMBS()->GetIOption(151)) // scale rigid body displacements
		{
			Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
			Matrix3D rot;
			GetRotMatrixFromAngles(XGD(4),XGD(5),XGD(6), rot);
			Matrix3D rotInit;
			GetRotMatrixFromAngles(XGD(4), GetXInit(5), GetXInit(6), rotInit);

			Vector3D unscaled_position = GetRefPosD()+n.GetLocalFrame()*rot*p_loc;
			Vector3D initial_position_rotated = GetRefPosD() + n.GetLocalFrame()*rotInit*p_loc;
			return initial_position_rotated + GetMBS()->GetDOption(105)*(unscaled_position - initial_position_rotated); // GetMBS()->GetDOption(105) ==> scaling factor: BodiesDeformationScaleFactor
		}
		else
		{
			Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
			Matrix3D rot;
			GetRotMatrixFromXGD(rot);
			return GetRefPosD()+n.GetLocalFrame()*rot*p_loc;
		}
	};

	virtual Vector3D GetPos(const Vector3D& p_loc) const
  {
		Node3DR123& n = (Node3DR123&)mbs->GetNode(nodenum);
		Matrix3D rot;
		GetRotMatrixFromXG(rot);
		return GetRefPos()+n.GetLocalFrame()*rot*p_loc;
	};

	virtual void GetNodedPosdqT(int node, Matrix& dpdqi) 
	{
		dpdqi.SetSize(6,3);
		dpdqi.FillWithZeros();

		dpdqi(1,1) = 1; //dnx/dnx
		dpdqi(2,2) = 1; //dny/dny
		dpdqi(3,3) = 1; //dny/dny
	};

	virtual void DrawElement();

private:
	int full_mass_matrix; //$EDC$[varaccess,int_bool,EDCvarname="full_mass_matrix",EDCfolder="Physics",tooltiptext="set to 1 if influence of tilted mass should be considered in the mass matrix"]
	Vector3D Ip; //$EDC$[varaccess,EDCvarname="moment_of_inertia",EDCfolder="Physics",tooltiptext="moments of inertia of the disk"]
	//EDC double size.X() //$EDC$[varaccess,remove,EDCvarname="radius",EDCfolder="Graphics"]
  //EDC double size.X() //$EDC$[varaccess,EDCvarname="thickness",EDCfolder="Graphics",tooltiptext="drawing thickness of disk mass"]
  //EDC double size.Y() //$EDC$[varaccess,EDCvarname="radius",EDCfolder="Graphics",tooltiptext="drawing radius of disk mass"]

	//EDC double mass; //$EDC$[varaccess,remove,EDCvarname="mass",EDCfolder="Physics"]
	//EDC double mass; //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="total mass of disk"] 

	//EDC Vector3D x_init //$EDC$[varaccess,remove,EDCvarname="initial_position",EDCfolder="Initialization"]
	//EDC Vector3D x_init //$EDC$[varaccess,remove,EDCvarname="initial_velocity",EDCfolder="Initialization"]
}; //$EDC$[endclass,NodalDiskMass3D]


//rigid cube
class Rigid3D: public Body3D //$EDC$[beginclass,classname=Rigid3D,parentclassname=Body3D,addelementtype=TAEBody,addelementtypename=Rigid3D,
//texdescription="A rigid body in 3D.",texdescriptionDOF="The first 3 degrees of freedom are those describing the position. The rotation is parameterized with 4 degrees of freedom and one additional algebraic equation.",
//texdescriptionGeometry="The center of gravity, S, is defined by the vector initial\_position, which is in global coordinates, see figure \ref{Rigid3Dfigure2}. The rotation of the body-fixed local coordinate system w.r.t. the global coordiante system is defined by the Matrix initial\_rotation. \\In order to define the position of a point P of the element, e.g. for connectors or sensors, the local coordinate system is used. The reference point is the center of mass, S, so the values of the local coordinates can be positive or negative.",
//example="Rigid3D.txt",figure="Rigid3D,Rigid3D",figure="Rigid3Dcoordinates,local and global coordinate system for a Rigid3D"]
{
public:
	//Body3D():Element() {mbs = NULL;};
	Rigid3D(MBS* mbsi):Body3D(mbsi) 
	{	
		ElementDefaultConstructorInitialization(); //$ DR 2012-07: ElementDefaultConstructorInitialization added
	};
	
	Rigid3D(const Rigid3D& e):Body3D(e.mbs) 
	{
		ElementDefaultConstructorInitialization(); //$ DR 2012-07: ElementDefaultConstructorInitialization added
		CopyFrom(e);
	};

	//x0 ... pos0 and vel0, phi0 ... initial euler angles, 3 quaternions_p
	Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi,
		const Vector3D& si, const Vector3D& coli);

	virtual void ElementDefaultConstructorInitialization(); //$ DR 2012-07: ElementDefaultConstructorInitialization added

	//$ DR 2012-07: CheckConsistency added
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute


	// set functions
	//x0 ... pos0 and vel0, phi0 ... initial euler angles, 3 quaternions_p
	virtual void SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	virtual void SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
		const Vector3D& si, const Vector3D& coli);

	virtual void SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi,
		const Vector3D& si, const Vector3D& coli);


	//compute initial vector from Euler angles:
	void ComputeInitialConditions(const Vector3D& xp, const Vector3D& vp, 
		const Vector3D& phi, const Vector3D& phip, Vector& xinit);

	//To be overwritten in derived class:
	virtual Element* GetCopy();
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor

	virtual int SetElementDataMassAndInertia(ElementDataContainer& edc); //$ DR 2012-12-13 sub function of SetElementData, used for Rigid3D and Rigid3DKardan

	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
	virtual const char* GetElementSpec() const {return "Rigid3D";}
	virtual void Initialize();

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_ref_conf_position+TKAF_ref_position+TKAF_D_pos_D_x+TKAF_D_pos_D_q+TKAF_angular_velocity+TKAF_D_ang_vel_D_q+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_rot_D_q+TKAF_D_rot_v_D_q+TKAF_acceleration);
	}

	virtual const Matrix3D& GetRotInertia() const {return Iphi;}
// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}

	virtual void EvalF(Vector& f, double t) {};  

	virtual void EvalG(Vector& f, double t);
	virtual void EvalM(Matrix& m, double t); 
	virtual void EvalF2(Vector& f, double t); 

	virtual int SOS() const {return 7;}; //size of second order equations, len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 1;};  //implicit (algebraic) size

	virtual int SOSowned_RS() const {return 8;}; //for resort, the algebraic element variables belong to SOS2
	virtual int IS_RS() const {return 0;};  //algebraic (implicit) size

	virtual int IsRigid() const {return 1;} //default value

	//deprecated
	//virtual void GetBetaP(double& betap0, double& betap1, double& betap2, double& betap3) const;
	//virtual void GetBeta(double& beta0, double& beta1, double& beta2, double& beta3) const;
	//virtual void GetBetaD(double& beta0, double& beta1, double& beta2, double& beta3) const;
	//virtual void GetBetaPD(double& betap0, double& betap1, double& betap2, double& betap3) const;
	//new
	virtual void GetBetaP(Vector& betap) const;
	virtual void GetBeta(Vector& beta) const;
	virtual void GetBetaD(Vector& beta) const;
	virtual void GetBetaInitD(Vector& beta) const;
	virtual void GetBetaPD(Vector& betap) const;
	
	virtual int GetIndBeta(int i) const {return 3+i;}
	virtual double GetLagrangeMultEP() const {return XG(15);}
	virtual void GetBetaPP(Vector& betapp) const//$ RL 2011-3-15: acceleration of generalized rotational dofs
	{	
		for(int i=1; i<=NRotParam();i++)
		{
			//betapp(i) = GetMBS()->GetAcceleration(LTG(i+3+SOS())); // $ MSax 2013-07-16 : removed
			betapp(i) = XGPP(i+3); // $ MSax 2013-07-16 : added
		}
	}

	////compute acceleration vector (after Timestep is successfully computed)
	virtual Vector3D GetAcceleration(const Vector3D& ploc) const//$ RL 2011-3-14: added for use in MBSSensor. Attention: Use this function only at end of computed time step (otherwise not tested).
	{
		//Term 1  +    Term 2       + Term 3
		//..   ..       _    _   _          _      _    ..
		//ri = Ri + Ai( w x (w x ui)) + Ai(alpha x u) = Ri 
		// Vector3D Rpp(GetMBS()->GetAcceleration(LTG(1 + SOS())), GetMBS()->GetAcceleration(LTG(2 + SOS())), GetMBS()->GetAcceleration(LTG(3 + SOS()))); // d²/dt²(R)  // $ MSax 2013-07-16 : removed
		Vector3D Rpp(XGPP(1),XGPP(2),XGPP(3)); // d²/dt²(R)  // $ MSax 2013-07-16 : added
		Matrix3D A = GetRotMatrix();

		ConstVector<4> betap(NRotParam());
		GetBetaP(betap);
		ConstVector<4> betapp(NRotParam());
		GetBetaPP(betapp);

		Vector3D omega_bar = GetGbar()*betap;
		Vector3D alpha_bar = GetGbar()*betapp; // Gbar.betapp (Shabana 1994)
		return Rpp + A * (omega_bar.Cross(omega_bar.Cross(ploc))+alpha_bar.Cross(ploc));

		//global vectors (alternative):
		//Vector3D omega = GetG()*betap;
		//Vector3D alpha = GetG()*betapp;
		//Vector3D u = A * ploc;
		//return Rpp + omega.Cross(omega.Cross(u))+alpha.Cross(u);
	}
	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const;
	virtual Vector3D GetAngularVelLocal() const;
	virtual Matrix3D GetG() const;
	virtual Matrix3D GetG(Vector& beta) const{assert(0 && "Use only in child class.");return Matrix3D(0.);}
	virtual Matrix3D GetGT() const;
	virtual Matrix3D GetGbar() const;
	virtual Matrix3D GetGbarT() const;
	virtual Matrix3D GetGbarp() const;
	virtual Matrix3D GetGbarpT() const;

	// Compute d omega / d theta = - dot Gbar
	virtual Matrix3D GetDOmegaDTheta(Vector& betap) const;
	//C_q^T*\lambda for Euler parameter equation:
	virtual void AddEPCqTterms(Vector& f);

	virtual Vector3D GetRefConfPos(const Vector3D& ploc) const { return GetRefPosInit()+GetRotMatrixInit()*ploc;}		//$ DR+PG 2013-05-15

	//rot matrix A describes the transformation of the local(body) to the
	//global coordinate system, such that r_g=R_g+A_bg*u_b
	virtual int NRotParam() const{return 4;};
	
	virtual Matrix3D GetRotMatrix() const
	{
		ConstVector<4> beta(NRotParam());
		GetBeta(beta);		
		return ComputeRotMatrix(beta);
	}
	virtual Matrix3D GetRotMatrixP() const
	{
		ConstVector<4> betap(NRotParam());
		GetBetaP(betap);
		ConstVector<4> beta(NRotParam());
		GetBeta(beta);

		return ComputeRotMatrixP(beta, betap);
	}

	//drawing matrices:
	virtual Matrix3D GetRotMatrixD() const
	{
		ConstVector<4> beta(NRotParam());
		GetBetaD(beta);
		return ComputeRotMatrix(beta);
	}
	//drawing matrices:
	virtual Matrix3D GetRotMatrixInit() const
	{
		ConstVector<4> beta(NRotParam());
		GetBetaInitD(beta);
		return ComputeRotMatrix(beta);
	}
	virtual Matrix3D GetRotMatrixPD() const
  {
		ConstVector<4> betap(NRotParam());
		GetBetaPD(betap);
		ConstVector<4> beta(NRotParam());
		GetBetaD(beta);

		return ComputeRotMatrixP(beta, betap);
		//return ComputeRotMatrixP(XGD(4), XGD(5), XGD(6), XGD(7), XGPD(4), XGPD(5), XGPD(6), XGPD(7));
	}

	virtual Matrix3D ComputeRotMatrix(const Vector& beta) const;

	virtual Matrix3D ComputeRotMatrixP(const Vector& beta, const Vector& betap) const;

	//for body loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x);

		//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x);

	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetDuxDq(Vector& dudq);
 
	virtual void GetDuyDq(Vector& dudq);

	virtual void GetDuzDq(Vector& dudq);

	virtual void GetH3T(const Vector3D& vloc, Matrix3D& d);

	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);

	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetdAngVeldqpT(const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx);

	virtual void DrawElement();

	virtual double GetSpecialSensorValue(int nr, double time) const;

	// virtual functions for reference frame
	virtual int SOSRigid() const { mbs->UO() << "Rigid3D::SOSRigid() called, should be called for FFRFElement only!!\n"; return SOS();}  // number of degrees of freedom from the rigid
	virtual int AddNode(Node* n) { mbs->UO() << "Rigid3D::AddNode(Node* n) called, shoulb be called for FFRFElement only!!\n"; return 0;} // add node to the floating frame
	virtual int IsCMS() const { mbs->UO() << "Rigid3D::IsCMS() called, should be called for FFRFElement only!!\n"; return 0;}
	virtual int IsGCMS() const { mbs->UO() << "Rigid3D::IsGCMS() called, should be called for FFRFElement only!!\n"; return 0;}

	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);

	virtual double GetVolume() const {return volume;}	//$ DR 2013-06-04

protected:
	//mechanical:
	//$ DR 2012-07: access for Get/SetElementDataAuto[
	Matrix3D Iphi;										//$EDC$[varaccess,EDCvarname="moment_of_inertia",EDCfolder="Physics",tooltiptext="[I_XX,I_XY,I_XZ; ...]"]
	double volume;										//$EDC$[varaccess,EDCvarname="volume",EDCfolder="Physics",tooltiptext="volume of the body in m*m*m"]
	Vector betap;											//do not use betap for graphics!!!!

	
	//EDC Vector3D size;							//$EDC$[varaccess,EDCvarname="body_dimensions",EDCfolder="Graphics",tooltiptext="Dimensions of a regular cube [L_x, L_y, L_z] in m"]
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=1, vecend=3,EDCvarname="initial_position",EDCfolder="Initialization",tooltiptext="[X, Y, Z]"]
	//EDC Vector x_init;							//$EDC$[varaccess,vecstart=NRotParam()+4, vecend=NRotParam()+6,EDCvarname="initial_velocity",EDCfolder="Initialization",tooltiptext="[X, Y, Z]"]
	//EDC double mass;								//$EDC$[varaccess,EDCvarname="mass",EDCfolder="Physics",tooltiptext="mass of the body in kg"]	


	//$ DR 2012-07: access for Get/SetElementDataAuto]

}; //$EDC$[endclass,Rigid3D]



#endif

