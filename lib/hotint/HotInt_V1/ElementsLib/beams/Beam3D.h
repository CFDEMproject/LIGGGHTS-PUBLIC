//#**************************************************************
//#
//# filename:             Beam3D.h
//#
//# author:               Martin Saxinger, PG
//#
//# generated:						26. April 2012
//# description:          
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
 
#ifndef Beam3D__H 
#define Beam3D__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//const int nmtorsion = 2; //number of modes/shape functions for torsional deformation
//const int nmaxial = 2; //number of modes/shape functions for axial deformation

// Beam element for rotor simulation, beam axis in z-direction, nodal coordinates
class Beam3D: public Body3D //$EDC$[beginclass,classname=Beam3D,parentclassname=Body3D,addelementtype=TAEBody,addelementtypename=LinearBeam3D,texdescription="The Beam3D element is
// a three dimensional elastic beam element which is aligned along the local x axis. It provides a decoupled calculation of bending in y and z direction, axial deformation in x direction
// and torsion about the x axis. Shear deformation is not considered. The decoupled calculation is a simplification of the real, nonlinear problem, but for small deformations the results coincidence highly with the exact solution.",texdescriptionDOF="Bending is described by 4 DOF, the number of DOF for axial deformation as well as torsion is 2. These 12 DOF are stored in two nodes i and j. The DOF vector of the LinearBeam3D read as follows 
//\begin{equation}
//\mathbf{q}^{(i)} = [\mathbf{q}^{(i)},\mathbf{q}^{(j)}] = [x^{(i)},y^{(i)},z^{(i)},\phi_x^{(i)},\phi_y^{(i)},\phi_z^{(i)},x^{(j)},y^{(j)},z^{(j)},\phi_x^{(j)},\phi_y^{(j)},\phi_z^{(j)}]^T.
//\end{equation}",
// texdescriptionGeometry="The beam geometry is fully defined by 2 'Node3DRxyz' nodes and a 'Beam3DProperties' material element. Beam length and orientation is specified due to node positions and the orientation of the first node. The beam cross section size is defined due to the material. See figure \ref{LinearBeam3Dfigure1} for more details. In order to define the position of point P of the element, e.g. for connectors or sensors, the local coordinate system is used. The origin of the local coordinate system is the center of gravity of the beam, $p_0$ is the vector to the center of gravity.",texdescriptionLimitations="Shear deformation is not considered. The decoupled calculation is a simplification of the real, nonlinear problem, but for small deformations the results coincidence highly with the exact solution.",texdescriptionNode="To create a new beam element the user has to define two 'Node3DRxyz' nodes i and j. Every node of this type has 6 DOF. The first 3 DOF describe the node displacement ($x, y, z$) w.r.t global coordinate system, the last 3 DOF are angles of rotation ($\phi_x, \phi_y, \phi_z$) w.r.t global coordinate system. All angles are considered as small (linearized angles).
// The reference positions of the nodes define the beam ends at initial configuration and so the length of the beam. The beam orientation is defined due to reference rot angles of node i. The advantage of using nodes with global DOF is the possibility to discretize a beam element into small beams easily without needing complicated constraint conditions. The beam elements do not even have to be aligned along a straight line.
// If using the same node number for the boundary point of the adjoint beams, beam elements are constrained automatically, see figure \ref{LinearBeam3Dfigure2}.",figure="beam3d,LinearBeam3D - Geometry",figure="beam3ddisc,LinearBeam3D - Nodes",example="LinearBeam3D.txt"]
{
public:
	//typedef Rigid3DKardan RIGID;
	Beam3D(MBS* mbsi) : Body3D(mbsi) 
	{	
		ElementDefaultConstructorInitialization();
	};

	Beam3D(const Beam3D& e):Body3D(e.mbs), x1(), w1(), q0() {CopyFrom(e);};

	void SetBeam3D(Vector3D p0i, Matrix3D rot0i, int n1i, int n2i, int materialnumi, 
                   		 Vector3D si, int axialdeformationI = 1, Vector3D coli = Vector3D(0.,0.,1.), Vector q0i = Vector(12));
	
	// standard set function with oriented nodes (default set function for script language)
	void SetBeam3D(int n1nr, int n2nr, int matnr);
	
	virtual void ElementDefaultConstructorInitialization();

	virtual Element* GetCopy()
	{
		Element* ec = new Beam3D(*this);
		return ec;
	}
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor()
	{
	}

	virtual void Initialize() 
	{
		Body3D::Initialize();
		CalcMK_Matrices();
	}

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_acceleration+TKAF_angular_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_pos_D_q+TKAF_D_pos_D_x+TKAF_D_rot_D_q+TKAF_D_rot_v_D_q+TKAF_int_D_u_D_q);
	}

	virtual void CalcMK_Matrices();
	virtual void LinkToElements();
	virtual int ElementDOFtoNodeDOF(int ind, int nodedim) const;

	virtual const char* GetElementSpec() const {return "LinearBeam3D";}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual mystr GetElementCoordinatesTexDescription();
	virtual void GetElementNodeImage(mystr& name, mystr& caption);

	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual int NNodes() const {return 2;}
	virtual int SOS() const {return FlexDOF();} // 
	virtual int SOSowned() const {return 0;}; //anz. eigene Freiheitsgrade (bei uns sind alle FG in den Knoten vorhanden) len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 1-axialdeformation;};  //implicit (algebraic) size
	inline int NSPos() const {return 4;} //cubic --> nmodes
	inline int NSTor() const {return 2;} //linear --> nmtorsion
	inline int NSAx() const {return 2;} //linear --> nmaxial

	virtual int FlexDOF() const {return 2*NSPos() + NSTor() + NSAx()*useAllDOF;} //number of flexible dof

	virtual int UseGlobalDOF() const {return 1;}

	virtual int& NodeNum(int i) //this int can be changed by user
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else 
		{
			GetMBS()->UO() << "Error in Beam3D::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	virtual double GetLagrangeMultEP() const {return XG(2*SOS()+1);}

	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));}

	virtual int AxialDeformation() const {return axialdeformation;}
	virtual int UseAllDOF() const {return useAllDOF;}

	virtual int IsRigid() const {return 0;}

	//$ MSax 2013-4-5:[
	virtual const double XGL(int iloc) const;
	virtual const double XGPL(int iloc) const;
	virtual const double XGDL(int iloc) const;
	virtual const double XGPDL(int iloc) const;
	virtual const double XGPPL(int iloc) const;
	virtual const double q0l(int iloc) const;

	virtual Matrix3D GetRefRotN1() const;
	virtual Matrix3D GetRefRotN2() const;

	//$ MSax 2013-4-5:]

	//x0 is normalized to -1..+1
	virtual double GetS0w(double x0, int shape) const;  //shape functions for bending
	virtual double GetS0wx(double x0, int shape) const; 
	virtual double GetS0wxx(double x0, int shape) const;
	virtual double GetS0wxxx(double x0, int shape) const; 
	virtual double GetS0v(double x0, int shape) const;  //shape functions for bending
	virtual double GetS0vx(double x0, int shape) const; 
	virtual double GetS0vxx(double x0, int shape) const;
	virtual double GetS0vxxx(double x0, int shape) const; 
	virtual double GetS0theta(double x0, int shape) const; //shape functions for torsion
	virtual double GetS0thetax(double x0, int shape) const;
	virtual double GetS0u(double x0, int shape) const;  //shape functions for axial deformation
	virtual double GetS0ux(double x0, int shape) const;
	virtual double GetS0beta(double x0, int shape) const;  //shape functions for shear deformation
	virtual double GetS0betax(double x0, int shape) const;

	virtual double GetU(double x0) const;  //axial deformation
	virtual double GetV(double x0) const;  //deflection in Y-direction
	virtual double GetW(double x0) const;	 //deflection in Z-direction
	virtual double GetUx(double x0) const; //d(U)/dx
	virtual double GetVx(double x0) const; //d(V)/dx
	virtual double GetWx(double x0) const; //d(W)/dx
	virtual double GetVxx(double x0) const; //d(V)/dxx
	virtual double GetWxx(double x0) const; //d(W)/dxx
	virtual double GetVxxx(double x0) const; //d(V)/dxx
	virtual double GetWxxx(double x0) const; //d(W)/dxx
	virtual double GetTheta(double x0) const; //Theta (theta = rotation of cross section, about beam axis)
	virtual double GetThetax(double x0) const; //d(Theta)/dx
	virtual double GetAngleX(double x0) const; //(linearized) rotation about local x-axis
	virtual double GetAngleY(double x0) const; //(linearized) rotation about local y-axis
	virtual double GetAngleZ(double x0) const; //(linearized) rotation about local z-axis

	virtual double GetdAngleXdx(double x0) const;
	virtual double GetdAngleYdx(double x0) const;
	virtual double GetdAngleZdx(double x0) const;

	virtual double GetUP(double x0) const;  //axial deformation
	virtual double GetVP(double x0) const;  //deflection in Y-direction
	virtual double GetWP(double x0) const;	 //deflection in Z-direction
	virtual double GetUxP(double x0) const; //d(U)/dx
	virtual double GetVxP(double x0) const; //d(V)/dx
	virtual double GetWxP(double x0) const; //d(W)/dx
	virtual double GetVxxP(double x0) const; //d(V)/dxx
	virtual double GetWxxP(double x0) const; //d(W)/dxx
	virtual double GetThetaP(double x0) const; //Theta (rotation of cross section)
	virtual double GetThetaxP(double x0) const; //d(Theta)/dx
	virtual double GetAngleXP(double x0) const; //(linearized) rotation about local x-axis
	virtual double GetAngleYP(double x0) const; //(linearized) rotation about local y-axis
	virtual double GetAngleZP(double x0) const; //(linearized) rotation about local z-axis

	/***********************************************************************/
	// $ MSax 2012-12-14: added accelerations
	virtual double GetUPP(double x0) const;  //axial deformation
	virtual double GetVPP(double x0) const;  //deflection in Y-direction
	virtual double GetWPP(double x0) const;	 //deflection in Z-direction

	virtual double GetAngleXPP(double x0) const; //(linearized) rotation about local x-axis
	virtual double GetAngleYPP(double x0) const; //(linearized) rotation about local y-axis
	virtual double GetAngleZPP(double x0) const; //(linearized) rotation about local z-axis

	/***********************************************************************/

	virtual double GetUD(double x0) const;  //axial deformation
	virtual double GetVD(double x0) const;  //deflection in Y-direction
	virtual double GetWD(double x0) const;	 //deflection in Z-direction
	virtual double GetUxD(double x0) const; //d(U)/dx
	virtual double GetVxD(double x0) const; //d(V)/dx
	virtual double GetWxD(double x0) const; //d(W)/dx
	virtual double GetVxxD(double x0) const; //d(V)/dxx
	virtual double GetWxxD(double x0) const; //d(W)/dxx
	virtual double GetVxxxD(double x0) const; //d(V)/dxx
	virtual double GetWxxxD(double x0) const; //d(W)/dxx
	virtual double GetThetaD(double x0) const; //Theta (rotation of cross section)
	virtual double GetThetaxD(double x0) const; //d(Theta)/dx
	virtual double GetAngleXD(double x0) const; //(linearized) rotation about local x-axis
	virtual double GetAngleYD(double x0) const; //(linearized) rotation about local y-axis
	virtual double GetAngleZD(double x0) const; //(linearized) rotation about local z-axis

	virtual double GetAngleXInit(double x0) const; //(linearized) rotation about local x-axis, in beam-coordinates (reference configuration)
	virtual double GetAngleYInit(double x0) const; //(linearized) rotation about local y-axis, in beam-coordinates (reference configuration)
	virtual double GetAngleZInit(double x0) const; //(linearized) rotation about local z-axis, in beam-coordinates (reference configuration)

	virtual double GetUInit(double x0) const;  //axial deformation (reference configuration)     //MS:  mit q0 oder XInit() oder beides addiert
	virtual double GetVInit(double x0) const;  //deflection in Y-direction (reference configuration)
	virtual double GetWInit(double x0) const;	 //deflection in Z-direction (reference configuration)

	virtual double GetVxInit(double x0) const; //d(V)/dx (reference configuration)
	virtual double GetWxInit(double x0) const; //d(W)/dx (reference configuration)
	virtual double GetThetaInit(double x0) const; //Theta (rotation of cross section) (reference configuration)

	virtual void GetRot(double x, Matrix3D& rot) const;  //linearized rotation matrix at global position x (-size.X()/2..+size.X()/2)
	virtual void GetdRotdx(double x, Matrix3D& rot) const; 
	virtual void GetRotP(double x, Matrix3D& rot) const; 
	virtual void GetRotPP(double x, Matrix3D& rot) const; 
	virtual void GetRotD(double x, Matrix3D& rot) const; 

	// $ MSax 2012-12-14: added accelerations
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	virtual double GetUPD(double x0) const;  //axial deformation
	virtual double GetVPD(double x0) const;  //deflection in Y-direction
	virtual double GetWPD(double x0) const;	 //deflection in Z-direction

	virtual double GetAngleXPD(double x0) const; //(linearized) rotation about local x-axis
	virtual double GetAngleYPD(double x0) const; //(linearized) rotation about local y-axis
	virtual double GetAngleZPD(double x0) const; //(linearized) rotation about local z-axis

	virtual double GetVxPD(double x0) const; //d(V)/dx
	virtual double GetWxPD(double x0) const; //d(W)/dx
	virtual double GetThetaPD(double x0) const; //Theta (rotation of cross section)

	virtual void GetRotPD(double x, Matrix3D& rot) const; 

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual void GetRotInit(double x, Matrix3D& rot) const; //linearized rotation matrix (reference configuration) at global position x (-size.X()/2..+size.X()/2)

	virtual Vector3D GetNodeLocPos(int i) const; //local position of node $ SW + MSax 2013-08-27: added

	virtual Vector3D GetPos(double x) const;
	virtual Vector3D GetPos(const Vector3D& p_loc) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const;

	virtual Vector3D GetAcceleration(const Vector3D& p_loc) const;

	virtual Vector3D GetPosD(const Vector3D& p_loc) const;
	virtual Vector3D GetVelD(const Vector3D& p_loc) const;

	virtual Vector3D GetInitRefPos() const {return p0;};

	virtual Vector3D GetRefPos() const {return GetPos(0.);}
	virtual Vector3D GetRefPosD() const {return GetPosD(Vector3D(0.,0.,0.));};

	virtual Vector3D GetDOFPosD(int idof) const; //returns postion of i-th DOF
	virtual Vector3D GetDOFDirD(int idof) const; //returns direction of action of i-th DOF

	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);

	virtual void EvalF2(Vector& f, double t); 

	//C_q^T*\lambda for Euler parameter equation:
	virtual void AddEPCqTterms(Vector& f);

	virtual void EvalG(Vector& f, double t); // ui-uj=0 

	virtual void TransformationToNodeDOF(Matrix& d);

	virtual int FastStiffnessMatrix() const {return 0;}

	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRot(ploc.X(), rot);
		return rot;
	}

	virtual Matrix3D GetRotMatrixP(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRotP(ploc.X(), rot);
		return rot;
	}

	virtual Matrix3D GetRotMatrixPD(const Vector3D& ploc) const //$ SW 2013-10-24: added
	{
		Matrix3D rot;
		GetRotPD(ploc.X(), rot);
		return rot;
	}

	virtual Matrix3D GetRotMatrixPP(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRotPP(ploc.X(), rot);
		return rot;
	}

	virtual Matrix3D GetRotMatrixD(const Vector3D& ploc) const 
	{
		Matrix3D rot;
		GetRotD(ploc.X(), rot);
		return rot;
	}

	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const 
	{
		Matrix3D omega_skew = /*rot0*/(GetRotMatrixP(p_loc)*GetRotMatrix(p_loc).GetTp()); // $ MSax 2013-04-19 : removed multiplication with rot0 from left because it's wrong: omega_skew: rot0*rotp*rotp^T*rot0^T

		return Vector3D(-omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2)); // $ MSax 2013-04-19 : changed sign from x coordinate of vector 3d
	}
	
	//$ SW 2013-10-24: added
	virtual Vector3D GetAngularVelD(const Vector3D& p_loc) const 
	{
		//see comments above
		Matrix3D omega_skew = (GetRotMatrixPD(p_loc)*GetRotMatrixD(p_loc).GetTp());
		return Vector3D(-omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2));
	}

	virtual void TransformMatrix(Matrix3D rot,Matrix& d); //transform dudq, etc. matrices by means of rot0 reference rotation

		//for body loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";

	}
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}
	virtual void ApplyDrotdq(Vector& f, const Vector3D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}
	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq);

	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);
	virtual void GetdRotvdqTLoc(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi); //$ SW 2013-08-28: added

	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx);

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, int local_node_number, bool flagD) //$ SW 2013-08-28: added
	{
		return GetFieldVariableValue(fvd,GetNodeLocPos(local_node_number),flagD);
	}

	virtual void DrawElement();

protected:
	//global DOF
	static const int beam3d_uN1 = 1;
	static const int beam3d_vN1 = 2;
	static const int beam3d_wN1 = 3;
	static const int beam3d_phiXN1 = 4;
	static const int beam3d_phiYN1 = 5;
	static const int beam3d_phiZN1 = 6;
	static const int beam3d_uN2 = 7;
	static const int beam3d_vN2 = 8;
	static const int beam3d_wN2 = 9;
	static const int beam3d_phiXN2 = 10;
	static const int beam3d_phiYN2 = 11;
	static const int beam3d_phiZN2 = 12;

	int beam3d_xg_v_list[4];
	int beam3d_xg_w_list[4];
	int beam3d_xg_theta_list[2];
	int beam3d_xg_u_list[2];
	int beam3d_xg_list[12];
	int beam3d_xg_inverse_list[12];

	//mechanical:
	//EDC Vector3D size;	//!EDC$[varaccess,minval=0,EDCvarname="body_dimensions",EDCfolder="Geometry",tooltiptext="dimensions of the beam. [L_x (length), L_y (width), L_z (height)]"]

	int n1;	//$EDC$[varaccess,minval=0,EDCvarname="node_1",EDCfolder="Geometry",tooltiptext="number of Node 1"] 
	int n2;	//$EDC$[varaccess,minval=0,EDCvarname="node_2",EDCfolder="Geometry",tooltiptext="number of Node 2"] 
	Vector3D p0;   //!EDC$[varaccess,EDCvarname="reference_point_p0",EDCfolder="Geometry",tooltiptext="reference point of center of gravity."]
	Matrix3D rot0; //!EDC$[varaccess,EDCvarname="reference_orientation_rot0",EDCfolder="Geometry",tooltiptext="rotation of reference (initial) configuration."]
	Vector q0; //initial vector of generalized coordinates

	int axialdeformation; //$EDC$[varaccess,int_bool,EDCvarname="axial_deformation",EDCfolder="Physics",tooltiptext="include effect of axial deformation"]
	
	int useAllDOF;	// is always 1, but will become important for future beams

	//integration points
	Vector x1,w1;
	
	Matrix Kv, Kw, Ktheta, Ku, Mv, Mw, Mtheta, Mu, M_with_u, M_without_u, M_with_u_global_dof;//, transformation_node_DOF;
	//EDC int materialnum; //$EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number which contains the main material properties of the beam"] 

}; //$EDC$[endclass,Beam3D]


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// Beam element for rotor simulation, beam axis in z-direction, nodal coordinates
class RotorBeamXAxis: public Beam3D//$EDC$[beginclass,classname=RotorBeamXAxis,parentclassname=Beam3D,addelementtype=TAEBody,addelementtypename=RotorBeamXAxis,texdescription="The RotorBeamXAxis element is
// a three dimensional elastic rotor beam element. It has exact the same characteristics and properties as the LinearBeam3D element except two differences. The first difference is that for a rotor element it is necessary to enable big rotation about the rotor axis instead of the small rotation of the LinearBeam3D. The second difference is that all element DOF are stored w.r.t. local beam coordinate system.",texdescriptionDOF="Bending is described by 4 DOF, the number of DOF for axial deformation as well as torsion is 2. These 12 DOF are stored in two nodes i and j. The DOF vector of the LinearBeam3D read as follows 
//\begin{equation}
//\mathbf{q}^{(i)} = [\mathbf{q}^{(i)},\mathbf{q}^{(j)}] = [x^{(i)},y^{(i)},z^{(i)},\phi_x^{(i)},\phi_y^{(i)},\phi_z^{(i)},x^{(j)},y^{(j)},z^{(j)},\phi_x^{(j)},\phi_y^{(j)},\phi_z^{(j)}]^T.
//\end{equation}",
// texdescriptionGeometry="The rotor beam geometry is fully defined by 2 'Node3DR123' nodes and a 'Beam3DProperties' material element. Beam length and orientation is specified due to node positions and the beam cross section size due to the material. The rotor beam has a circular cross section.",texdescriptionNode="To create a new rotor beam element the user has to define two 'Node3DR123' nodes i and j. Every node of this type has 6 DOF. The first 3 DOF describe the node displacement ($x, y, z$) w.r.t local rotor element coordinate system, the last 3 DOF are angles of rotation ($\phi_x, \phi_y, \phi_z$) w.r.t local rotor element coordinate system. The rotation about the local x-axis is considered as large, the rotations about the local y and z-axes are considered as small (linearized angles).
// The reference positions of the nodes define the beam ends at initial configuration and so the length of the beam. The beam orientation is defined due to reference rot angles of node i.",example="RotorBeamXAxis.txt",figure="RotorBeamXAxis"]
{
public:
	RotorBeamXAxis(MBS* mbsi):Beam3D(mbsi) 
	{
		elementname = GetElementSpec();
	};

	void SetRotorBeamXAxis(Vector3D p0i, Matrix3D rot0i, int n1i, int n2i, int materialnumi, 
                   		 Vector2D dimi, int axialdeformationI = 1, Vector3D coli = Vector3D(0.,0.,1.), Vector q0i = Vector(12));

	// standard set function with oriented nodes (default set function for script language)
	void SetRotorBeamXAxis(int n1nr, int n2nr, int matnr);

  RotorBeamXAxis(const RotorBeamXAxis& e):Beam3D(e.mbs) 
	{
		elementname = GetElementSpec();

		CopyFrom(e);
	};

	virtual Element* GetCopy()
	{
		Element* ec = new RotorBeamXAxis(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)
	{
		Beam3D::CopyFrom(e);
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute

	virtual int UseGlobalDOF() const {return 0;}

	virtual Matrix3D GetRefRotN1() const; // only used for global to local DOF transformation
	virtual Matrix3D GetRefRotN2() const; // only used for global to local DOF transformation

	virtual const char* GetElementSpec() const {return "RotorBeamXAxis";}

	virtual void GetRot(double x, Matrix3D& rot) const;  //linearized rotation matrix at global position x (-lx/2..+lx/2)
	virtual void GetRotP(double x, Matrix3D& rot) const; //linearized rotation matrix at global position x (-lx/2..+lx/2)
	virtual void GetRotD(double x, Matrix3D& rot) const; //linearized rotation matrix at global position x (-lx/2..+lx/2)
	virtual void GetRotPD(double x, Matrix3D& rot) const; 
	virtual void GetRotInit(double x, Matrix3D& rot) const; //linearized rotation matrix (reference configuration) at global position x (-lx/2..+lx/2)

	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);
	virtual void GetdRotvdqTLoc(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);

	virtual void DrawElement();

	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx);
	virtual void GetNodedPosdqT(int node, Matrix& dpdqi); //$ SW: 2013-08-28: added

	virtual void GetRotTheta(double x, Matrix3D& rot) const;
	virtual void GetRotThetaP(double x, Matrix3D& rot) const;
	virtual void GetRotThetaD(double x, Matrix3D& rot) const;
	virtual void GetRotThetaPD(double x, Matrix3D& rot) const;

	virtual Vector3D GetPosD(const Vector3D& p_loc) const;

	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

private:
  //EDC Vector3D size; //$EDC$[varaccess,remove,EDCvarname="body_dimensions",EDCfolder="Geometry"]
}; //$EDC$[endclass,RotorBeamXAxis]

#endif