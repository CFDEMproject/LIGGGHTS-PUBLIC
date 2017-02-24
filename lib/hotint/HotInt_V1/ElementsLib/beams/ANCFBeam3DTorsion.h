//#**************************************************************
//# filename:             ANCFBeam3DTorsion.h
//#
//# author:               PG & KN
//#
//# generated:						2012
//# description:          3D ANCF beam element with torsion, without shear deformation (BE beam theory)
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
 

//#include "../modelslib/models/mynode.h"

#ifndef ANCFBeam3DTorsion__H
#define ANCFBeam3DTorsion__H

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeam3DTorsion    ANCFBeam3DTorsion    ANCFBeam3DTorsion    ANCFBeam3DTorsion 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ANCFBeam3DTorsion: public Body3D  //$EDC$[beginclass,classname=ANCFBeam3DTorsion,parentclassname=Body3D,addelementtype=TAEBody,addelementtypename=ANCFBeam3DTorsion,
//texdescription="ANCFBeam3DTorsion is a Bernoulli-Euler beam finite element in ANCF (Absolute Nodal Coordinate Formulation) capable of large axial, bending, and torsional deformations.",
//texdescriptionNode="The element operates with two Nodes of type $\mathtt{Node3DS1rot1}$, each of which are located at either tip of the beam element. The integer values $\mathtt{Geometry.node\_number1}$ and $\mathtt{Geometry.node\_number2}$ refer to the index of the nodes in the multibody system. Each of these Nodes is instantiated by the user with a position and a rotation (kardan angles), and provides a frame $\left(\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3\right)$ (which is measured in the global frame of the multibody system) for the instantiation of the beam elemtent: At each node, the slope of the beam axis $\mathbf{r}'$ is identical with $\mathbf{e}_1$, and the director is defined as $\mathbf{e}_3$.",
//figure="ANCFBeam3DTorsion_Geometry, The geometry of the element is defined by nodal values for (a) the axial position, (b) the axial slope vector, and (c) the torsional angle of the cross section around the beam axis. This angle is measured with respect to a reference direction in the global frame (director). Between the nodal values, the axial position is interpolated cubically, the axial slope is interpolated quadratically, and the torsional angle of the cross section (around the beam axis) as well as the director are interpolated linearly.",
//texdescriptionGeometry="The geometry of the element is defined by the nodal values for axial position $\mathbf{r}$, the axial slope vector $\mathbf{r}'$, and the torsional angle $\theta$ of the cross section around the beam axis, see Fig.~\ref{ANCFBeam3DTorsionfigure2}. This angle is measured with respect to a reference direction in the global frame (director). Between the nodal values, the axial position is interpolated cubically, the axial slope is interpolated quadratically, and the torsional angle of the cross section (around the beam axis) as well as the director are interpolated linearly.",
//figure="ANCFBeam3DTorsion_DOFs, Ordering of the generalized coordinates.",
//texdescriptionDOF="The element affects $14$ degrees of freedom (generalized coordinates) in total, which are $7$ degrees of freedom per node, i.e., at each node we have: the axial displacement $\mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the change of the axial slope $\mathbf{u}' = \mathbf{r}'-\mathbf{r}'_0$, and the change of the torsional angle $\theta-\theta_0$. Each quantity with index $0$ here confers to the reference configuration. The element wise ordering of the degrees of freedom is displayed in Fig.~\ref{ANCFBeam3DTorsionfigure2}.",
//example="ANCFBeam3DTorsion.txt",
//reference="K. Nachbagauer, P. Gruber, Yu. Vetyukov, J. Gerstmayr. A spatial thin beam finite element based on the absolute nodal coordinate formulation without singularities. Proceedings of the ASME 2011 International Design Engineering Technical Conferences, Computers and Information in Engineering Conference IDETC/CIE 2011, Paper No. DETC2011/MSNDC-47732, Washington, DC, USA, 2011.",
//reference="P. Gruber, K. Nachbagauer, Yu. Vetyukov, J. Gerstmayr. A novel director-based Bernoulli-Euler beam finite element in absolute nodal coordinate formulation free of geometric singularities. Mechanical Science, 2013 (to appear).",
//texdescriptionComments="For details on the element, such as the definition of the elastic forces and the kinetic terms, see \cite{ANCFBeam3DTorsionreference1,ANCFBeam3DTorsionreference2}."]
{
public:
	ANCFBeam3DTorsion(MBS* mbsi):Body3D(mbsi), massmatrix(), Hmatrix(), x1(), w1(), q0() 
	{
		ElementDefaultConstructor();
	}
	ANCFBeam3DTorsion(const ANCFBeam3DTorsion& e):Body3D(e.mbs), massmatrix(), Hmatrix(), x1(), w1(), q0() 
	{
		CopyFrom(e);
	};

	// standard set function
	//
	// xc1: position and axial slope at node 1 (size 6 Vector)
	// theta1: torsional angle at node 1
	// director1: arbitrary vector of inertial system at node 1. must not point into direction of axial slope vector. local frame is defined by projection of director into normal plane of axial slope vector, and rotation around torsional angle theta.
	// n1i: global node index of node 1
	// xc2, theta2, director2, n2i: same at node 2
	// materialnumi: global index of material to be used
	// coli: body color of this element
	// do_update_directorsi: bool flag, if 1, then directors are updated after every time or load step, so that they keep perpendicular to the slope vectors at the nodes
	// kinematic_computation_modei: int flag, 0 ... exact terms, 5th order gaussian integration (slow);
	//                                        1 ... exact terms, low order (1st order lobatto) integration (fast);
	//                                        2 ... approximate mass matrix (torsional terms approximated), no quadratic velocity vector (fastest)
	void SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		const Vector3D& director1, const Vector3D& director2, int n1i, int n2i, int materialnumi,
		const Vector3D& coli, const int do_update_directorsi, const int kinematic_computation_modei);

	// standard set functions with oriented nodes
	void SetANCFBeam3DTorsion(int n1nr, int n2nr, int matnr, const Vector3D& coli = Vector3D(0.2,0.2,0.8), const int do_update_directorsi = 1, const int kinematic_computation_modei = 0);  // default set function for script language
	
	// alternative set function with explicit size
	void SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		const Vector3D& director1, const Vector3D& director2, int n1i, int n2i, int materialnumi, const Vector3D& si, const Vector3D& coli, const int do_update_directorsi = 1, const int kinematic_computation_modei=0);

	// alternative set function: determines directors automatically (see description of standard set function)
	// to be used with caution: torsional angle theta is measured with respect to these automatically generated directors!!
	void SetANCFBeam3DTorsion(const Vector& xc1, const Vector& xc2, const double& theta1, const double& theta2,
		int n1i, int n2i, int materialnumi, const Vector3D& si, const Vector3D& coli, const int do_update_directorsi = 1, const int kinematic_computation_modei=0);
	

private:
	// determine one standard director (for the whole element) by the knowledge of the axial slopes in the nodes
	// if slopes are significantly distinct, then the director is chosen to be the cross product of both
	// else
	//   if axial direction is not identical with 3rd kartesian direction then 3rd kartesian direction is chosen for director
	//   else 2nd kartesian direction is chosen for director
	Vector3D DetermineStandardDirector(Vector3D slopevector1, Vector3D slopevector2) const;
	double CalculateElementLength() const;

public:
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeam3DTorsion(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "ANCFBeam3DTorsion";}

	virtual void CopyFrom(const Element& e)
	{
		Body3D::CopyFrom(e);
		const ANCFBeam3DTorsion& ce = (const ANCFBeam3DTorsion&)e;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		kinematic_computation_mode = ce.kinematic_computation_mode;

		//integration points
		x1 = ce.x1;
		w1 = ce.w1;
		intorder_mass = ce.intorder_mass;
		intorder_axial_strain = ce.intorder_axial_strain;
		intorder_curvature = ce.intorder_curvature;

		q0 = ce.q0;  //initial values
		n1 = ce.n1; n2 = ce.n2;
		
		do_update_directors = ce.do_update_directors;

		materialnum = ce.materialnum;
	}

	virtual void Initialize();
	virtual void LinkToElements();

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter


	//nnodes =2
  virtual int NSPos() const {return 2*NNodes();}//= 4 r1,r1',r2,r2'
	virtual int NSRot() const {return NNodes();}//= 2 //SOSRot()=2, angle is only one component
	virtual int SOSPos() const {return Dim()*NSPos();}//=12
	virtual int SOS() const {return SOSPos()+NSRot();}//= 14   //Second Order Size (SO-DOFS owned by FE-nodes)
	virtual int SOSowned() const {return 0;}; //SO-DOFS owned by the element itself (i.e., not by nodes or other elements)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	//virtual int DataS() const {return 0;}	//Data size for non-state variables: for this element (Vector3D director1 + Vector3D director2)
	virtual int DataS() const {return 6;}	  //2*3 director components (left node director + right node director are saved by XData(), between the nodes director is interpolated linearly! should be better kind of interpolation later on!!!

	virtual int NNodes() const {return 2;} //Get-Fkt, we need a Set-function

	virtual const int& NodeNum(int i) const //this int stays const
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else 
		{
			mbs->UO() << "Error in ANCFBeam3DTorsion::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}
	virtual int& NodeNum(int i) //this int can be changed by user
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else 
		{
			mbs->UO() << "Error in ANCFBeam3DTorsion::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}
	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));}

	virtual void SetMaterialNum(int matnr)
	{
		materialnum = matnr;
	}
	
	//rename the functions in Mathematica-Input
	static double Cos(double phi)
	{
		return cos(phi);
	}
	static double Sin(double phi)
	{
		return sin(phi);
	}
	static double Power(double x, double y)
	{
		if (y == 2) return x*x;
		if (y == 3) return x*x*x;
		else return pow(x,y);
	}
	static double Sqrt(double x)
	{
		//$ YV 2012-12-11: commented out
		//if (x < 0) mbs->UO() << "Sqrt(" << x << ")\n";
		return sqrt(x);
	}

#include "BEBeamGetKappa1.h"
#include "BEBeamGetKappa2.h"
#include "BEBeamGetKappa3.h"

#include "BEBeamGetDKappa1.h"
#include "BEBeamGetDKappa2.h"
#include "BEBeamGetDKappa3.h"

#include "BEBeamdRotdr.h"
#include "BeamBEdRotdrT.h"

	void SetBeamParameters(double beamEIi, double beamEAi, double beamRhoAi)
	{
		GetMaterial().BeamEIy() = beamEIi;
		GetMaterial().BeamEA() = beamEAi;
		GetMaterial().BeamRhoA() = beamRhoAi;
	}

	virtual int IsRigid() const {return 0;}

	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}

//----------------
//Shapefunctions
//----------------
	//Shapefunctions for position -> length(sf)=4
	virtual void GetShapesPos(Vector& sf, double xi) const;
	virtual double GetSFPos(int i, double xi) const;
	virtual void GetShapesPosx(Vector& sf, double xi) const;
	virtual double GetSFPosx(int i, double xi) const;
	virtual void GetShapesPosxx(Vector& sf, double xi) const;
	virtual double GetSFPosxx(int i, double xi) const;
	virtual void GetShapesPosxxx(Vector& sf, double xi) const;
	virtual double GetSFPosxxx(int i, double xi) const;

	//Shapefunctions for rotation -> length(sf)=2
	virtual void GetShapesRot(Vector& sf, double xi) const;
	virtual double GetSFRot(int i, double xi) const;
	virtual void GetShapesRotx(Vector& sf, double xi) const;
	virtual double GetSFRotx(int i, double xi) const;

	//----------------
	//GetPos
	//----------------
	//r, r', r'', r''' (flagD for Visualization)
	virtual Vector3D GetPos(double xi, const Vector& xg) const;
	virtual Vector3D GetPos(double xi, int flagD=0) const;
	virtual Vector3D GetPosx(double xi, const Vector& xg) const;
	virtual Vector3D GetPosx(double xi, int flagD=0) const;
	virtual Vector3D GetPosxx(double xi, const Vector& xg) const;
	virtual Vector3D GetPosxx(double xi, int flagD=0) const;
	virtual Vector3D GetPosxxx(double xi, const Vector& xg) const;
  virtual Vector3D GetPosxxx(double xi, int flagD=0) const;
	virtual Vector3D GetPosxP(double xi, int flagD=0) const;

	//Get actual position of relative point p_loc in range [-lx/2..lx/2, etc.]
	virtual Vector3D GetPos(const Vector3D & p0) const; //$ JG

	// just for drawing...
	virtual Vector3D GetPosD(const Vector3D & p0) const;
	virtual Vector3D GetRefPosD() const;
	virtual Vector3D GetDisplacement(const Vector3D& p0, int flagD) const;
	virtual Vector3D GetVel(double xi, int flagD) const;
	virtual Vector3D GetVel(const Vector3D& p_loc) const;
	virtual Vector3D GetVelD(const Vector3D& p_loc) const;
	virtual Vector3D GetNodePos(int node_idx) const; //returns position of i-th node
	virtual Vector3D GetNodePosD(int node_idx) const; //returns position of i-th node (draw mode)
	virtual Vector3D GetDOFPosD(int idof) const; //returns position of i-th DOF
	virtual Vector3D GetDOFDirD(int idof) const; //returns direction of action of i-th DOF


	//----------------
	//GetRot
	//----------------
	//theta, theta' (flagD for Visualization)
	virtual double GetRot(double xi, const Vector& xg) const;  //faster - call XG() avoided
	virtual double GetRot(double xi, int flagD=0) const;
	virtual double GetRotx(double xi, const Vector& xg) const;  //faster - call XG() avoided
	virtual double GetRotx(double xi, int flagD=0) const;
	virtual double GetRotP(double xi, int flagD=0) const;

	//----------------
	//GetRotMatrix & derivative
	//----------------
	virtual Matrix3D GetRotMatrix(double xi, int flagD=0) const;
	virtual Matrix3D GetRotMatrixP(double xi, int flagD=0) const;
	virtual Matrix3D GetRotMatrix(const Vector3D& p_loc) const {return GetRotMatrix(p_loc(1), 0);};
	virtual Matrix3D GetRotMatrixD(const Vector3D& p_loc) const {return GetRotMatrix(p_loc(1), 1);};
	virtual Matrix3D GetRotMatrixP(const Vector3D& p_loc) const {return GetRotMatrixP(p_loc(1), 0);};
	virtual Matrix3D GetRotMatrixPD(const Vector3D& p_loc) const {return GetRotMatrixP(p_loc(1), 1);};
	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT);
	//alternative variant, where you can decide whether to use A or A^T for multiplication with vloc
	void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT, int use_transposed);

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//ONLY for beam axis ==> PG CHANGE!!!!!
	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Integration of du/dq = SFPos over axis
	virtual void GetIntDuDq(Matrix& dudq);

	//----------------
	//Init Position
	//----------------
	virtual Vector3D GetInitPos3D(double xi) const;
	virtual Vector3D GetInitPosx3D(double xi) const;
	virtual double GetInitRot3D(double xi) const;
	virtual double GetInitRotx3D(double xi) const;
	virtual Matrix3D GetRotMatrixInit() const { return GetInitRotMatrix3D(0.); };
	virtual Matrix3D GetInitRotMatrix3D(double xi) const;

	//1 director per node (vectors in 3D), interpolated linearly between nodes
	virtual Vector3D GetDirector(const double& xi, int flagD=0) const;    //PG: now (fast solution) linear interpolation between left and right directors (in kartesian coordinates)!
	virtual Vector3D GetDirectorx(const double& xi, int flagD=0) const;   //derivative of director wrt xi at xi
	virtual Vector3D GetInitDirector(const double& xi) const;
	void SetDirector1(const Vector3D& directori);
	void SetDirector2(const Vector3D& directori);
	Vector3D GetDirector1(int flagd=0) const;
	Vector3D GetDirector2(int flagd=0) const;
	int UpdateDirectors(const Vector3D& r1dx, const Vector3D& r2dx, const double theta1, const double theta2);

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode) const
	{
		TKinematicsAccessFunctions kaf = Body3D::GetKinematicsAccessFunctions(mode);
		return TKinematicsAccessFunctions(kaf+TKAF_node_position+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_pos_D_q+TKAF_D_rot_D_q+TKAF_D_rot_v_D_q+TKAF_int_D_u_D_q);
	}

	//Gamma = r' - A r0', Gamma1 = (|r'|-|r0'|)/|r0'|
	virtual void GetGamma1(const double& xi, double& Gamma1, int flagD=0) const;
	virtual void GetDeltaGamma1(const double& xi, Vector& DeltaGamma1) const;   //derivative of Gamma wrt q at xi
	
	//Kappa = k - A k0 = Kappa_i(r', r'', theta, theta') e_i
	virtual void GetKappa1(const double& xi, double& Kappa1, int flagD=0) const;
	virtual void GetKappa2(const double& xi, double& Kappa2, int flagD=0) const;
	virtual void GetKappa3(const double& xi, double& Kappa3, int flagD=0) const;
	virtual void GetDeltaKappa1(const double& xi, Vector& DeltaKappa1) const;   //derivative of Kappa1 wrt q at xi
	virtual void GetDeltaKappa2(const double& xi, Vector& DeltaKappa2) const;   //derivative of Kappa1 wrt q at xi
	virtual void GetDeltaKappa3(const double& xi, Vector& DeltaKappa3) const;   //derivative of Kappa1 wrt q at xi


	void GetQuadraticVelocityVector(const double& xi, Vector& qvv) const; //eqs.(47)-(48), xi in [-lx/2, lx/2]
	void GetQuadraticVelocityVectorAtNode(Vector& f, int n) const; // used for 2nd order lobatto integration (IPs are in nodes)
	void GetM(const double& xi, Matrix& m) const; //eq.(49), xi in [-lx/2, lx/2]
	void GetdMdqk(const double& xi, const int k, Matrix& dMdqk) const; //eq.(51), xi in [-lx/2, lx/2], k in {1,..,14}
	void GetL0(const double& xi, Matrix& L0) const; //eq.(42), xi in [-lx/2, lx/2]
	void Getdeidq(const double& xi, const int i, Matrix& deidq) const; //eq.(66)-(68), xi in [-lx/2, lx/2], i in {2,3}
	void Getdeidposx(const double& xi, const int i, Matrix& deidposx) const; //eq.(66)-(68), xi in [-lx/2, lx/2], i in {2,3}
	void Getdei0dposx(const double& xi, const int i, Matrix& dei0dposx) const; //eq.(69)-(76) + (55)-(56), xi in [-lx/2, lx/2], i in {2,3}
	void Getdei0hatdposx(const double& xi, const int i, Matrix& dei0hatdposx) const; //eq.(69)-(76), xi in [-lx/2, lx/2], i in {2,3}
	void Getdeidrot(const double& xi, const int i, Vector& deidrot) const; //eq.(68)+(80), xi in [-lx/2, lx/2], i in {2,3}
	void Getei0(const double& xi, const int i, Vector& ei0) const; //eq.(59)-(62) + (54), xi in [-lx/2, lx/2], i in {2,3}
	void Getei0hat(const double& xi, const int i, Vector& ei0hat) const; //eq.(59)-(62), xi in [-lx/2, lx/2], i in {2,3}

	void GetdL0dqk(const double& xi, const int k, Matrix& dL0dqk) const; //eq.(53), xi in [-lx/2, lx/2], k in {1,..,14}
	void Getddeidqdqk(const double& xi, const int i, const int k, Matrix& ddeidqdqk) const; //eq.(81), xi in [-lx/2, lx/2], i in {2,3}, k in {1,..,14}
	void Getddeidposxdposx(const double& xi, const int i, Matrix& ddeidposxdposx) const; //eq.(79)+(91), xi in [-lx/2, lx/2], i in {2,3}
	void Getddeidposxdrot(const double& xi, const int i, Matrix& ddeidposxdrot) const; //eq.(87)+(92), xi in [-lx/2, lx/2], i in {2,3}
	void Getddeidrotdrot(const double& xi, const int i, Vector& ddeidrotdrot) const; //eq.(88)+(93), xi in [-lx/2, lx/2], i in {2,3}
	void Getddei0dposxdposx(const double& xi, const int i, Matrix& ddei0dposxdposx) const; //eq.(57)-(58)+(79)+(91), xi in [-lx/2, lx/2], i in {2,3}
	void Getddei0hatdposxdposx(const double& xi, const int i, Matrix& ddei0hatdposxdposx) const; //eq.(80)-(84), xi in [-lx/2, lx/2], i in {2,3}

	void Getde1dposx(const double& xi, Matrix& de1dposx) const; //xi in [-lx/2, lx/2], e1 = posx/|posx|
	void Getdde1dposxdposx(const double& xi, Matrix& dde1dposxdposx) const; //xi in [-lx/2, lx/2], e1 = posx/|posx|
	void Getde30hatde1(const double& xi, Matrix& de30hatde1) const; //xi in [-lx/2, lx/2], e1 = posx/|posx|
	void Getdde30hatde1de1(const double& xi, Matrix& dde30hatde1de1) const; //xi in [-lx/2, lx/2], e1 = posx/|posx|, dde30hatde1de1 has 3 rows (according to de30hat) and 3*3 columns (according to de1)

	void GetFOverAbsF(const Vector& f, Vector& f_over_abs_f) const; //eq.(54)  computes (f/|f|)', where f maps a scalar to a vector
	void GetdFOverAbsFdx(const Vector& f, const Vector& dfdx, Vector& d_f_over_abs_f_prime_dx) const; //eq.(55)-(56)  computes (f/|f|)', where f maps a scalar to a vector   //move to linalg.h/cpp
	void GetddFOverAbsFdxdx(const Vector& f, const Vector& dfdx, const Vector& ddfdxdx, Vector& dd_f_over_abs_f_dxdx) const; //eq.(57)-(58)  computes (f/|f|)'', where f maps a scalar to a vector   //move to linalg.h/cpp
	void GetddFOverAbsFdxdy(const Vector& f, const Vector& dfdx, const Vector& dfdy, const Vector& ddfdxdy, Vector& dd_f_over_abs_f_dxdy) const; //eq.(57)-(58)  computes (f/|f|)'', where f maps a scalar to a vector   //move to linalg.h/cpp

	// for faster calculation: fixed drdx1, ..
	void Getdeidposx_and_deidtheta ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 , c37 , c38 , c39 , c40 , c41 , c42 , c43 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = d3*drdx1-d1*drdx3 ; c4 = 2*d3*c3-2*d2*c2 ; c5 = sqrt(pow(c3,2)+pow(c2,2)+pow(c1,2)) ; c6 = 1/pow(c5,3) ; c7 = cos(theta) ; c8 = pow(drdx1,2) ; c9 = pow(drdx2,2) ; c10 = pow(drdx3,2) ; c11 = sqrt(c10+c9+c8) ; c12 = 1/c11 ; c13 = d3*drdx3*c12+d2*drdx2*c12+d1*drdx1*c12 ; c14 = d1-drdx1*c12*c13 ; c15 = 1/pow(c11,3) ; c16 = d1*c12-d3*drdx1*drdx3*c15-d2*drdx1*drdx2*c15-d1*c8*c15 ; c17 = -c12*c13 ; c18 = c17+c8*c15*c13-drdx1*c12*c16 ; c19 = drdx1*drdx2*c15*c13 ; c20 = c19-drdx2*c12*c16 ; c21 = d2-drdx2*c12*c13 ; c22 = drdx1*drdx3*c15*c13 ; c23 = c22-drdx3*c12*c16 ; c24 = d3-drdx3*c12*c13 ; c25 = 2*c23*c24+2*c20*c21+2*c18*c14 ; c26 = sqrt(pow(c24,2)+pow(c21,2)+pow(c14,2)) ; c27 = 1/pow(c26,3) ; c28 = sin(theta) ; c29 = 1/c26 ; c30 = 1/c5 ; c31 = 2*d1*c2-2*d3*c1 ; c32 = d2*c12-d3*drdx2*drdx3*c15-d2*c9*c15-d1*drdx1*drdx2*c15 ; c33 = c19-drdx1*c12*c32 ; c34 = c17+c9*c15*c13-drdx2*c12*c32 ; c35 = drdx2*drdx3*c15*c13 ; c36 = c35-drdx3*c12*c32 ; c37 = 2*c36*c24+2*c34*c21+2*c33*c14 ; c38 = 2*d2*c1-2*d1*c3 ; c39 = d3*c12-d3*c10*c15-d2*drdx2*drdx3*c15-d1*drdx1*drdx3*c15 ; c40 = c22-drdx1*c12*c39 ; c41 = c35-drdx2*c12*c39 ; c42 = c17+c10*c15*c13-drdx3*c12*c39 ; c43 = 2*c42*c24+2*c41*c21+2*c40*c14 ; ret.SetXN(24,c18*c29*c28-c14*c25*c27*c28/2-c1*c4*c6*c7/2,c20*c29*c28-c21*c25*c27*c28/2+d3*c30*c7-c3*c4*c6*c7/2,c23*c29*c28-c24*c25*c27*c28/2-d2*c30*c7-c2*c4*c6*c7/2,c33*c29*c28-c14*c37*c27*c28/2-d3*c30*c7-c1*c31*c6*c7/2,c34*c29*c28-c21*c37*c27*c28/2-c3*c31*c6*c7/2,c36*c29*c28-c24*c37*c27*c28/2+d1*c30*c7-c2*c31*c6*c7/2,c40*c29*c28-c14*c43*c27*c28/2+d2*c30*c7-c1*c38*c6*c7/2,c41*c29*c28-c21*c43*c27*c28/2-d1*c30*c7-c3*c38*c6*c7/2,c42*c29*c28-c24*c43*c27*c28/2-c2*c38*c6*c7/2,c14*c29*c7-c1*c30*c28,c21*c29*c7-c3*c30*c28,c24*c29*c7-c2*c30*c28,c1*c4*c6*c28/2+c18*c29*c7-c14*c25*c27*c7/2,-d3*c30*c28+c3*c4*c6*c28/2+c20*c29*c7-c21*c25*c27*c7/2,d2*c30*c28+c2*c4*c6*c28/2+c23*c29*c7-c24*c25*c27*c7/2,d3*c30*c28+c1*c31*c6*c28/2+c33*c29*c7-c14*c37*c27*c7/2,c3*c31*c6*c28/2+c34*c29*c7-c21*c37*c27*c7/2,-d1*c30*c28+c2*c31*c6*c28/2+c36*c29*c7-c24*c37*c27*c7/2,-d2*c30*c28+c1*c38*c6*c28/2+c40*c29*c7-c14*c43*c27*c7/2,d1*c30*c28+c3*c38*c6*c28/2+c41*c29*c7-c21*c43*c27*c7/2,c2*c38*c6*c28/2+c42*c29*c7-c24*c43*c27*c7/2,-c14*c29*c28-c1*c30*c7,-c21*c29*c28-c3*c30*c7,-c24*c29*c28-c2*c30*c7);} 
	void Getddeidposx1dposx1 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = d3*drdx1-d1*drdx3 ; c4 = 2*d3*c3-2*d2*c2 ; c5 = pow(c4,2) ; c6 = sqrt(pow(c3,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 2*pow(d3,2)+2*pow(d2,2) ; c10 = 1/pow(c6,3) ; c11 = pow(drdx1,2) ; c12 = sqrt(pow(drdx3,2)+pow(drdx2,2)+c11) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d1*c13-d3*drdx1*drdx3*c16-d2*drdx1*drdx2*c16-d1*c11*c16 ; c18 = -c13*c14+c11*c16*c14-drdx1*c13*c17 ; c19 = drdx1*drdx2*c16*c14-drdx2*c13*c17 ; c20 = d2-drdx2*c13*c14 ; c21 = drdx1*drdx3*c16*c14-drdx3*c13*c17 ; c22 = d3-drdx3*c13*c14 ; c23 = 2*c21*c22+2*c19*c20+2*c18*c15 ; c24 = pow(c23,2) ; c25 = sqrt(pow(c22,2)+pow(c20,2)+pow(c15,2)) ; c26 = 1/pow(c25,5) ; c27 = sin(theta) ; c28 = 1/pow(c25,3) ; c29 = pow(drdx1,3) ; c30 = 1/pow(c12,5) ; c31 = -d3*drdx3*c16-d2*drdx2*c16-3*d1*drdx1*c16+3*d3*c11*drdx3*c30+3*d2*c11*drdx2*c30+3*d1*c29*c30 ; c32 = 3*drdx1*c16*c14-3*c29*c30*c14-2*c13*c17+2*c11*c16*c17-drdx1*c13*c31 ; c33 = drdx2*c16*c14-3*c11*drdx2*c30*c14+2*drdx1*drdx2*c16*c17-drdx2*c13*c31 ; c34 = drdx3*c16*c14-3*c11*drdx3*c30*c14+2*drdx1*drdx3*c16*c17-drdx3*c13*c31 ; c35 = 2*c34*c22+2*pow(c21,2)+2*c33*c20+2*pow(c19,2)+2*pow(c18,2)+2*c32*c15 ; c36 = 1/c25 ; ret.SetXN(6,c32*c36*c27-c15*c35*c28*c27/2-c18*c23*c28*c27+3*c15*c24*c26*c27/4-c9*c1*c10*c8/2+3*c1*c5*c7*c8/4,c33*c36*c27-c20*c35*c28*c27/2-c19*c23*c28*c27+3*c20*c24*c26*c27/4-d3*c4*c10*c8-c9*c3*c10*c8/2+3*c3*c5*c7*c8/4,c34*c36*c27-c22*c35*c28*c27/2-c21*c23*c28*c27+3*c22*c24*c26*c27/4+d2*c4*c10*c8-c9*c2*c10*c8/2+3*c2*c5*c7*c8/4,c9*c1*c10*c27/2-3*c1*c5*c7*c27/4+c32*c36*c8-c15*c35*c28*c8/2-c18*c23*c28*c8+3*c15*c24*c26*c8/4,d3*c4*c10*c27+c9*c3*c10*c27/2-3*c3*c5*c7*c27/4+c33*c36*c8-c20*c35*c28*c8/2-c19*c23*c28*c8+3*c20*c24*c26*c8/4,-d2*c4*c10*c27+c9*c2*c10*c27/2-3*c2*c5*c7*c27/4+c34*c36*c8-c22*c35*c28*c8/2-c21*c23*c28*c8+3*c22*c24*c26*c8/4); } 
	void Getddeidposx2dposx1 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 , c37 , c38 , c39 , c40 , c41 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = d3*drdx1-d1*drdx3 ; c4 = 2*d3*c3-2*d2*c2 ; c5 = 2*d1*c2-2*d3*c1 ; c6 = sqrt(pow(c3,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 1/pow(c6,3) ; c10 = pow(drdx1,2) ; c11 = pow(drdx2,2) ; c12 = sqrt(pow(drdx3,2)+c11+c10) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d1*c13-d3*drdx1*drdx3*c16-d2*drdx1*drdx2*c16-d1*c10*c16 ; c18 = -c13*c14 ; c19 = c18+c10*c16*c14-drdx1*c13*c17 ; c20 = drdx1*drdx2*c16*c14 ; c21 = c20-drdx2*c13*c17 ; c22 = d2-drdx2*c13*c14 ; c23 = drdx1*drdx3*c16*c14-drdx3*c13*c17 ; c24 = d3-drdx3*c13*c14 ; c25 = 2*c23*c24+2*c21*c22+2*c19*c15 ; c26 = d2*c13-d3*drdx2*drdx3*c16-d2*c11*c16-d1*drdx1*drdx2*c16 ; c27 = c20-drdx1*c13*c26 ; c28 = c18+c11*c16*c14-drdx2*c13*c26 ; c29 = drdx2*drdx3*c16*c14-drdx3*c13*c26 ; c30 = 2*c29*c24+2*c28*c22+2*c27*c15 ; c31 = sqrt(pow(c24,2)+pow(c22,2)+pow(c15,2)) ; c32 = 1/pow(c31,5) ; c33 = sin(theta) ; c34 = 1/pow(c12,5) ; c35 = -d1*drdx2*c16-d2*drdx1*c16+3*d3*drdx1*drdx2*drdx3*c34+3*d2*drdx1*c11*c34+3*d1*c10*drdx2*c34 ; c36 = drdx2*c16*c14-3*c10*drdx2*c34*c14-c13*c26+c10*c16*c26+drdx1*drdx2*c16*c17-drdx1*c13*c35 ; c37 = drdx1*c16*c14-3*drdx1*c11*c34*c14+drdx1*drdx2*c16*c26-c13*c17+c11*c16*c17-drdx2*c13*c35 ; c38 = -3*drdx1*drdx2*drdx3*c34*c14+drdx1*drdx3*c16*c26+drdx2*drdx3*c16*c17-drdx3*c13*c35 ; c39 = 2*c38*c24+2*c37*c22+2*c36*c15+2*c21*c28+2*c27*c19+2*c23*c29 ; c40 = 1/pow(c31,3) ; c41 = 1/c31 ; ret.SetXN(6,c36*c41*c33-c19*c30*c40*c33/2-c27*c25*c40*c33/2-c15*c39*c40*c33/2+3*c15*c25*c30*c32*c33/4+d3*c4*c9*c8/2+d1*d2*c1*c9*c8+3*c1*c4*c5*c7*c8/4,c37*c41*c33-c21*c30*c40*c33/2-c28*c25*c40*c33/2-c22*c39*c40*c33/2+3*c22*c25*c30*c32*c33/4-d3*c5*c9*c8/2+d1*d2*c3*c9*c8+3*c3*c4*c5*c7*c8/4,c38*c41*c33-c23*c30*c40*c33/2-c29*c25*c40*c33/2-c24*c39*c40*c33/2+3*c24*c25*c30*c32*c33/4+d2*c5*c9*c8/2-d1*c4*c9*c8/2+d1*d2*c2*c9*c8+3*c2*c4*c5*c7*c8/4,-d3*c4*c9*c33/2-d1*d2*c1*c9*c33-3*c1*c4*c5*c7*c33/4+c36*c41*c8-c19*c30*c40*c8/2-c27*c25*c40*c8/2-c15*c39*c40*c8/2+3*c15*c25*c30*c32*c8/4,d3*c5*c9*c33/2-d1*d2*c3*c9*c33-3*c3*c4*c5*c7*c33/4+c37*c41*c8-c21*c30*c40*c8/2-c28*c25*c40*c8/2-c22*c39*c40*c8/2+3*c22*c25*c30*c32*c8/4,-d2*c5*c9*c33/2+d1*c4*c9*c33/2-d1*d2*c2*c9*c33-3*c2*c4*c5*c7*c33/4+c38*c41*c8-c23*c30*c40*c8/2-c29*c25*c40*c8/2-c24*c39*c40*c8/2+3*c24*c25*c30*c32*c8/4); } 
	void Getddeidposx3dposx1 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 , c37 , c38 , c39 , c40 , c41 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = d3*drdx1-d1*drdx3 ; c4 = 2*d3*c3-2*d2*c2 ; c5 = 2*d2*c1-2*d1*c3 ; c6 = sqrt(pow(c3,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 1/pow(c6,3) ; c10 = pow(drdx1,2) ; c11 = pow(drdx3,2) ; c12 = sqrt(pow(drdx2,2)+c11+c10) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d1*c13-d3*drdx1*drdx3*c16-d2*drdx1*drdx2*c16-d1*c10*c16 ; c18 = -c13*c14 ; c19 = c18+c10*c16*c14-drdx1*c13*c17 ; c20 = drdx1*drdx2*c16*c14-drdx2*c13*c17 ; c21 = d2-drdx2*c13*c14 ; c22 = drdx1*drdx3*c16*c14 ; c23 = c22-drdx3*c13*c17 ; c24 = d3-drdx3*c13*c14 ; c25 = 2*c23*c24+2*c20*c21+2*c19*c15 ; c26 = d3*c13-d3*c11*c16-d2*drdx2*drdx3*c16-d1*drdx1*drdx3*c16 ; c27 = c22-drdx1*c13*c26 ; c28 = drdx2*drdx3*c16*c14-drdx2*c13*c26 ; c29 = c18+c11*c16*c14-drdx3*c13*c26 ; c30 = 2*c29*c24+2*c28*c21+2*c27*c15 ; c31 = sqrt(pow(c24,2)+pow(c21,2)+pow(c15,2)) ; c32 = 1/pow(c31,5) ; c33 = sin(theta) ; c34 = 1/pow(c12,5) ; c35 = -d1*drdx3*c16-d3*drdx1*c16+3*d3*drdx1*c11*c34+3*d2*drdx1*drdx2*drdx3*c34+3*d1*c10*drdx3*c34 ; c36 = drdx3*c16*c14-3*c10*drdx3*c34*c14-c13*c26+c10*c16*c26+drdx1*drdx3*c16*c17-drdx1*c13*c35 ; c37 = -3*drdx1*drdx2*drdx3*c34*c14+drdx1*drdx2*c16*c26+drdx2*drdx3*c16*c17-drdx2*c13*c35 ; c38 = drdx1*c16*c14-3*drdx1*c11*c34*c14+drdx1*drdx3*c16*c26-c13*c17+c11*c16*c17-drdx3*c13*c35 ; c39 = 2*c38*c24+2*c37*c21+2*c36*c15+2*c23*c29+2*c27*c19+2*c20*c28 ; c40 = 1/pow(c31,3) ; c41 = 1/c31 ; ret.SetXN(6,c36*c41*c33-c19*c30*c40*c33/2-c27*c25*c40*c33/2-c15*c39*c40*c33/2+3*c15*c25*c30*c32*c33/4-d2*c4*c9*c8/2+d1*d3*c1*c9*c8+3*c1*c4*c5*c7*c8/4,c37*c41*c33-c20*c30*c40*c33/2-c28*c25*c40*c33/2-c21*c39*c40*c33/2+3*c21*c25*c30*c32*c33/4-d3*c5*c9*c8/2+d1*c4*c9*c8/2+d1*d3*c3*c9*c8+3*c3*c4*c5*c7*c8/4,c38*c41*c33-c23*c30*c40*c33/2-c29*c25*c40*c33/2-c24*c39*c40*c33/2+3*c24*c25*c30*c32*c33/4+d2*c5*c9*c8/2+d1*d3*c2*c9*c8+3*c2*c4*c5*c7*c8/4,d2*c4*c9*c33/2-d1*d3*c1*c9*c33-3*c1*c4*c5*c7*c33/4+c36*c41*c8-c19*c30*c40*c8/2-c27*c25*c40*c8/2-c15*c39*c40*c8/2+3*c15*c25*c30*c32*c8/4,d3*c5*c9*c33/2-d1*c4*c9*c33/2-d1*d3*c3*c9*c33-3*c3*c4*c5*c7*c33/4+c37*c41*c8-c20*c30*c40*c8/2-c28*c25*c40*c8/2-c21*c39*c40*c8/2+3*c21*c25*c30*c32*c8/4,-d2*c5*c9*c33/2-d1*d3*c2*c9*c33-3*c2*c4*c5*c7*c33/4+c38*c41*c8-c23*c30*c40*c8/2-c29*c25*c40*c8/2-c24*c39*c40*c8/2+3*c24*c25*c30*c32*c8/4); } 
	void Getddeidposx2dposx2 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = 2*d1*c2-2*d3*c1 ; c4 = pow(c3,2) ; c5 = d3*drdx1-d1*drdx3 ; c6 = sqrt(pow(c5,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 2*pow(d3,2)+2*pow(d1,2) ; c10 = 1/pow(c6,3) ; c11 = pow(drdx2,2) ; c12 = sqrt(pow(drdx3,2)+pow(drdx1,2)+c11) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d2*c13-d3*drdx2*drdx3*c16-d2*c11*c16-d1*drdx1*drdx2*c16 ; c18 = drdx1*drdx2*c16*c14-drdx1*c13*c17 ; c19 = -c13*c14+c11*c16*c14-drdx2*c13*c17 ; c20 = d2-drdx2*c13*c14 ; c21 = drdx2*drdx3*c16*c14-drdx3*c13*c17 ; c22 = d3-drdx3*c13*c14 ; c23 = 2*c21*c22+2*c19*c20+2*c18*c15 ; c24 = pow(c23,2) ; c25 = sqrt(pow(c22,2)+pow(c20,2)+pow(c15,2)) ; c26 = 1/pow(c25,5) ; c27 = sin(theta) ; c28 = 1/pow(c25,3) ; c29 = 1/pow(c12,5) ; c30 = pow(drdx2,3) ; c31 = -d3*drdx3*c16-3*d2*drdx2*c16-d1*drdx1*c16+3*d3*c11*drdx3*c29+3*d2*c30*c29+3*d1*drdx1*c11*c29 ; c32 = drdx1*c16*c14-3*drdx1*c11*c29*c14+2*drdx1*drdx2*c16*c17-drdx1*c13*c31 ; c33 = 3*drdx2*c16*c14-3*c30*c29*c14-2*c13*c17+2*c11*c16*c17-drdx2*c13*c31 ; c34 = drdx3*c16*c14-3*c11*drdx3*c29*c14+2*drdx2*drdx3*c16*c17-drdx3*c13*c31 ; c35 = 2*c34*c22+2*pow(c21,2)+2*c33*c20+2*pow(c19,2)+2*pow(c18,2)+2*c32*c15 ; c36 = 1/c25 ; ret.SetXN(6,c32*c36*c27-c15*c35*c28*c27/2-c18*c23*c28*c27+3*c15*c24*c26*c27/4+d3*c3*c10*c8-c9*c1*c10*c8/2+3*c1*c4*c7*c8/4,c33*c36*c27-c20*c35*c28*c27/2-c19*c23*c28*c27+3*c20*c24*c26*c27/4-c9*c5*c10*c8/2+3*c5*c4*c7*c8/4,c34*c36*c27-c22*c35*c28*c27/2-c21*c23*c28*c27+3*c22*c24*c26*c27/4-d1*c3*c10*c8-c9*c2*c10*c8/2+3*c2*c4*c7*c8/4,-d3*c3*c10*c27+c9*c1*c10*c27/2-3*c1*c4*c7*c27/4+c32*c36*c8-c15*c35*c28*c8/2-c18*c23*c28*c8+3*c15*c24*c26*c8/4,c9*c5*c10*c27/2-3*c5*c4*c7*c27/4+c33*c36*c8-c20*c35*c28*c8/2-c19*c23*c28*c8+3*c20*c24*c26*c8/4,d1*c3*c10*c27+c9*c2*c10*c27/2-3*c2*c4*c7*c27/4+c34*c36*c8-c22*c35*c28*c8/2-c21*c23*c28*c8+3*c22*c24*c26*c8/4); } 
	void Getddeidposx3dposx2 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 , c37 , c38 , c39 , c40 , c41 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d3*drdx1-d1*drdx3 ; c3 = 2*d2*c1-2*d1*c2 ; c4 = d1*drdx2-d2*drdx1 ; c5 = 2*d1*c4-2*d3*c1 ; c6 = sqrt(pow(c4,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 1/pow(c6,3) ; c10 = pow(drdx2,2) ; c11 = pow(drdx3,2) ; c12 = sqrt(pow(drdx1,2)+c11+c10) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d2*c13-d3*drdx2*drdx3*c16-d2*c10*c16-d1*drdx1*drdx2*c16 ; c18 = drdx1*drdx2*c16*c14-drdx1*c13*c17 ; c19 = -c13*c14 ; c20 = c19+c10*c16*c14-drdx2*c13*c17 ; c21 = d2-drdx2*c13*c14 ; c22 = drdx2*drdx3*c16*c14 ; c23 = c22-drdx3*c13*c17 ; c24 = d3-drdx3*c13*c14 ; c25 = 2*c23*c24+2*c20*c21+2*c18*c15 ; c26 = d3*c13-d3*c11*c16-d2*drdx2*drdx3*c16-d1*drdx1*drdx3*c16 ; c27 = drdx1*drdx3*c16*c14-drdx1*c13*c26 ; c28 = c22-drdx2*c13*c26 ; c29 = c19+c11*c16*c14-drdx3*c13*c26 ; c30 = 2*c29*c24+2*c28*c21+2*c27*c15 ; c31 = sqrt(pow(c24,2)+pow(c21,2)+pow(c15,2)) ; c32 = 1/pow(c31,5) ; c33 = sin(theta) ; c34 = 1/pow(c12,5) ; c35 = -d2*drdx3*c16-d3*drdx2*c16+3*d3*drdx2*c11*c34+3*d2*c10*drdx3*c34+3*d1*drdx1*drdx2*drdx3*c34 ; c36 = -3*drdx1*drdx2*drdx3*c34*c14+drdx1*drdx2*c16*c26+drdx1*drdx3*c16*c17-drdx1*c13*c35 ; c37 = drdx3*c16*c14-3*c10*drdx3*c34*c14-c13*c26+c10*c16*c26+drdx2*drdx3*c16*c17-drdx2*c13*c35 ; c38 = drdx2*c16*c14-3*drdx2*c11*c34*c14+drdx2*drdx3*c16*c26-c13*c17+c11*c16*c17-drdx3*c13*c35 ; c39 = 2*c38*c24+2*c37*c21+2*c36*c15+2*c23*c29+2*c28*c20+2*c18*c27 ; c40 = 1/pow(c31,3) ; c41 = 1/c31 ; ret.SetXN(6,c36*c41*c33-c18*c30*c40*c33/2-c27*c25*c40*c33/2-c15*c39*c40*c33/2+3*c15*c25*c30*c32*c33/4-d2*c5*c9*c8/2+d3*c3*c9*c8/2+d2*d3*c1*c9*c8+3*c1*c3*c5*c7*c8/4,c37*c41*c33-c20*c30*c40*c33/2-c28*c25*c40*c33/2-c21*c39*c40*c33/2+3*c21*c25*c30*c32*c33/4+d1*c5*c9*c8/2+d2*d3*c2*c9*c8+3*c2*c3*c5*c7*c8/4,c38*c41*c33-c23*c30*c40*c33/2-c29*c25*c40*c33/2-c24*c39*c40*c33/2+3*c24*c25*c30*c32*c33/4-d1*c3*c9*c8/2+d2*d3*c4*c9*c8+3*c4*c3*c5*c7*c8/4,d2*c5*c9*c33/2-d3*c3*c9*c33/2-d2*d3*c1*c9*c33-3*c1*c3*c5*c7*c33/4+c36*c41*c8-c18*c30*c40*c8/2-c27*c25*c40*c8/2-c15*c39*c40*c8/2+3*c15*c25*c30*c32*c8/4,-d1*c5*c9*c33/2-d2*d3*c2*c9*c33-3*c2*c3*c5*c7*c33/4+c37*c41*c8-c20*c30*c40*c8/2-c28*c25*c40*c8/2-c21*c39*c40*c8/2+3*c21*c25*c30*c32*c8/4,d1*c3*c9*c33/2-d2*d3*c4*c9*c33-3*c4*c3*c5*c7*c33/4+c38*c41*c8-c23*c30*c40*c8/2-c29*c25*c40*c8/2-c24*c39*c40*c8/2+3*c24*c25*c30*c32*c8/4); } 
	void Getddeidposx3dposx3 ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d3*drdx1-d1*drdx3 ; c3 = 2*d2*c1-2*d1*c2 ; c4 = pow(c3,2) ; c5 = d1*drdx2-d2*drdx1 ; c6 = sqrt(pow(c5,2)+pow(c2,2)+pow(c1,2)) ; c7 = 1/pow(c6,5) ; c8 = cos(theta) ; c9 = 2*pow(d2,2)+2*pow(d1,2) ; c10 = 1/pow(c6,3) ; c11 = pow(drdx3,2) ; c12 = sqrt(pow(drdx2,2)+pow(drdx1,2)+c11) ; c13 = 1/c12 ; c14 = d3*drdx3*c13+d2*drdx2*c13+d1*drdx1*c13 ; c15 = d1-drdx1*c13*c14 ; c16 = 1/pow(c12,3) ; c17 = d3*c13-d3*c11*c16-d2*drdx2*drdx3*c16-d1*drdx1*drdx3*c16 ; c18 = drdx1*drdx3*c16*c14-drdx1*c13*c17 ; c19 = drdx2*drdx3*c16*c14-drdx2*c13*c17 ; c20 = d2-drdx2*c13*c14 ; c21 = -c13*c14+c11*c16*c14-drdx3*c13*c17 ; c22 = d3-drdx3*c13*c14 ; c23 = 2*c21*c22+2*c19*c20+2*c18*c15 ; c24 = pow(c23,2) ; c25 = sqrt(pow(c22,2)+pow(c20,2)+pow(c15,2)) ; c26 = 1/pow(c25,5) ; c27 = sin(theta) ; c28 = 1/pow(c25,3) ; c29 = 1/pow(c12,5) ; c30 = pow(drdx3,3) ; c31 = -3*d3*drdx3*c16-d2*drdx2*c16-d1*drdx1*c16+3*d3*c30*c29+3*d2*drdx2*c11*c29+3*d1*drdx1*c11*c29 ; c32 = drdx1*c16*c14-3*drdx1*c11*c29*c14+2*drdx1*drdx3*c16*c17-drdx1*c13*c31 ; c33 = drdx2*c16*c14-3*drdx2*c11*c29*c14+2*drdx2*drdx3*c16*c17-drdx2*c13*c31 ; c34 = 3*drdx3*c16*c14-3*c30*c29*c14-2*c13*c17+2*c11*c16*c17-drdx3*c13*c31 ; c35 = 2*c34*c22+2*pow(c21,2)+2*c33*c20+2*pow(c19,2)+2*pow(c18,2)+2*c32*c15 ; c36 = 1/c25 ; ret.SetXN(6,c32*c36*c27-c15*c35*c28*c27/2-c18*c23*c28*c27+3*c15*c24*c26*c27/4-d2*c3*c10*c8-c9*c1*c10*c8/2+3*c1*c4*c7*c8/4,c33*c36*c27-c20*c35*c28*c27/2-c19*c23*c28*c27+3*c20*c24*c26*c27/4+d1*c3*c10*c8-c9*c2*c10*c8/2+3*c2*c4*c7*c8/4,c34*c36*c27-c22*c35*c28*c27/2-c21*c23*c28*c27+3*c22*c24*c26*c27/4-c9*c5*c10*c8/2+3*c5*c4*c7*c8/4,d2*c3*c10*c27+c9*c1*c10*c27/2-3*c1*c4*c7*c27/4+c32*c36*c8-c15*c35*c28*c8/2-c18*c23*c28*c8+3*c15*c24*c26*c8/4,-d1*c3*c10*c27+c9*c2*c10*c27/2-3*c2*c4*c7*c27/4+c33*c36*c8-c20*c35*c28*c8/2-c19*c23*c28*c8+3*c20*c24*c26*c8/4,c9*c5*c10*c27/2-3*c5*c4*c7*c27/4+c34*c36*c8-c22*c35*c28*c8/2-c21*c23*c28*c8+3*c22*c24*c26*c8/4); } 
	void Getddeidthetadposx ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 , c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 , c28 , c29 , c30 , c31 , c32 , c33 , c34 , c35 , c36 , c37 , c38 , c39 , c40 , c41 , c42 , c43 ; c1 = pow(drdx1,2) ; c2 = pow(drdx2,2) ; c3 = pow(drdx3,2) ; c4 = sqrt(c3+c2+c1) ; c5 = 1/c4 ; c6 = d3*drdx3*c5+d2*drdx2*c5+d1*drdx1*c5 ; c7 = d1-drdx1*c5*c6 ; c8 = 1/pow(c4,3) ; c9 = d1*c5-d3*drdx1*drdx3*c8-d2*drdx1*drdx2*c8-d1*c1*c8 ; c10 = -c5*c6 ; c11 = c10+c1*c8*c6-drdx1*c5*c9 ; c12 = drdx1*drdx2*c8*c6 ; c13 = c12-drdx2*c5*c9 ; c14 = d2-drdx2*c5*c6 ; c15 = drdx1*drdx3*c8*c6 ; c16 = c15-drdx3*c5*c9 ; c17 = d3-drdx3*c5*c6 ; c18 = 2*c16*c17+2*c13*c14+2*c11*c7 ; c19 = sqrt(pow(c7,2)+pow(c17,2)+pow(c14,2)) ; c20 = 1/pow(c19,3) ; c21 = cos(theta) ; c22 = 1/c19 ; c23 = d2*drdx3-d3*drdx2 ; c24 = d1*drdx2-d2*drdx1 ; c25 = d3*drdx1-d1*drdx3 ; c26 = 2*d3*c25-2*d2*c24 ; c27 = sqrt(pow(c25,2)+pow(c24,2)+pow(c23,2)) ; c28 = 1/pow(c27,3) ; c29 = sin(theta) ; c30 = 1/c27 ; c31 = d2*c5-d3*drdx2*drdx3*c8-d2*c2*c8-d1*drdx1*drdx2*c8 ; c32 = c12-drdx1*c5*c31 ; c33 = c10+c2*c8*c6-drdx2*c5*c31 ; c34 = drdx2*drdx3*c8*c6 ; c35 = c34-drdx3*c5*c31 ; c36 = 2*c35*c17+2*c33*c14+2*c32*c7 ; c37 = 2*d1*c24-2*d3*c23 ; c38 = d3*c5-d3*c3*c8-d2*drdx2*drdx3*c8-d1*drdx1*drdx3*c8 ; c39 = c15-drdx1*c5*c38 ; c40 = c34-drdx2*c5*c38 ; c41 = c10+c3*c8*c6-drdx3*c5*c38 ; c42 = 2*c41*c17+2*c40*c14+2*c39*c7 ; c43 = 2*d2*c23-2*d1*c25 ; ret.SetXN(18,c23*c26*c28*c29/2+c11*c22*c21-c7*c18*c20*c21/2,-d3*c30*c29+c25*c26*c28*c29/2+c13*c22*c21-c14*c18*c20*c21/2,d2*c30*c29+c24*c26*c28*c29/2+c16*c22*c21-c17*c18*c20*c21/2,d3*c30*c29+c23*c37*c28*c29/2+c32*c22*c21-c7*c36*c20*c21/2,c25*c37*c28*c29/2+c33*c22*c21-c14*c36*c20*c21/2,-d1*c30*c29+c24*c37*c28*c29/2+c35*c22*c21-c17*c36*c20*c21/2,-d2*c30*c29+c23*c43*c28*c29/2+c39*c22*c21-c7*c42*c20*c21/2,d1*c30*c29+c25*c43*c28*c29/2+c40*c22*c21-c14*c42*c20*c21/2,c24*c43*c28*c29/2+c41*c22*c21-c17*c42*c20*c21/2,-c11*c22*c29+c7*c18*c20*c29/2+c23*c26*c28*c21/2,-c13*c22*c29+c14*c18*c20*c29/2-d3*c30*c21+c25*c26*c28*c21/2,-c16*c22*c29+c17*c18*c20*c29/2+d2*c30*c21+c24*c26*c28*c21/2,-c32*c22*c29+c7*c36*c20*c29/2+d3*c30*c21+c23*c37*c28*c21/2,-c33*c22*c29+c14*c36*c20*c29/2+c25*c37*c28*c21/2,-c35*c22*c29+c17*c36*c20*c29/2-d1*c30*c21+c24*c37*c28*c21/2,-c39*c22*c29+c7*c42*c20*c29/2-d2*c30*c21+c23*c43*c28*c21/2,-c40*c22*c29+c14*c42*c20*c29/2+d1*c30*c21+c25*c43*c28*c21/2,-c41*c22*c29+c17*c42*c20*c29/2+c24*c43*c28*c21/2); } 
	void Getddeidthetadtheta ( TArray<double>& ret, double drdx1, double drdx2, double drdx3, double theta, double d1, double d2, double d3 ) const { double c1 , c2 , c3 , c4 , c5 , c6 , c7 , c8 , c9 , c10 , c11 , c12 ; c1 = d2*drdx3-d3*drdx2 ; c2 = d1*drdx2-d2*drdx1 ; c3 = d3*drdx1-d1*drdx3 ; c4 = 1/sqrt(pow(c3,2)+pow(c2,2)+pow(c1,2)) ; c5 = cos(theta) ; c6 = 1/sqrt(pow(drdx3,2)+pow(drdx2,2)+pow(drdx1,2)) ; c7 = d3*drdx3*c6+d2*drdx2*c6+d1*drdx1*c6 ; c8 = d1-drdx1*c6*c7 ; c9 = d2-drdx2*c6*c7 ; c10 = d3-drdx3*c6*c7 ; c11 = 1/sqrt(pow(c9,2)+pow(c8,2)+pow(c10,2)) ; c12 = sin(theta) ; ret.SetXN(6,-c8*c11*c12-c1*c4*c5,-c9*c11*c12-c3*c4*c5,-c10*c11*c12-c2*c4*c5,c1*c4*c12-c8*c11*c5,c3*c4*c12-c9*c11*c5,c2*c4*c12-c10*c11*c5); } 





	virtual void GetCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XG(i);
	}

	virtual void GetCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGP(i);
	}

	virtual void GetDrawCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGD(i);
	}

	virtual void GetDrawCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGPD(i);
	}
	
	virtual void EvalM(Matrix& m, double t);   // siehe ANCFBeam3D...
	virtual void EvalF2(Vector& f, double t);
	virtual void EvalQVV(Vector& f, double t);
	virtual double PostNewtonStep(double t);

	//for visualization
	virtual Vector3D GetPos3D0D(const Vector3D& p_loc) const ;
	virtual Vector3D GetPos3D0D(const Vector3D& p_loc, double defscale) const ;

	virtual void DrawElement();
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

	virtual void StartTimeStep();

protected:

	void ElementDefaultConstructor()
	{
		SetElementName(GetElementSpec());

		// make sure, that all members are set to a value, which are declared 'readonly' for the EDCConverter (see declaration of members below)
		intorder_mass = 4;
		intorder_axial_strain = 9;
		intorder_curvature = 5;

		kinematic_computation_mode = 0;

		// make sure, to set the size of all dynamical objects, which are used by EDCConverter
		q0.SetLen(SOS());
		q0.SetAll(0.);
		size = Vector3D(1,0.1,0.1);

		// some other default values of members, e.g. global node numbers
		n1=1;
		n2=2;

		x_init.SetLen(2*SOS());
		x_init.SetAll(0);

		Vector3D director;
		Vector3D slopevector1(1,0,0);   /// PG - get vaules from nodes
		Vector3D slopevector2(1,0,0);
		director = DetermineStandardDirector(slopevector1, slopevector2);
		Vector datainit(director.X(), director.Y(), director.Z(), director.X(), director.Y(), director.Z());
		datainit(1)=1.;
		SetDataInit(datainit);

		do_update_directors = 0; //$ DR added
	}


	//for faster calculation:
	//mass and stiffness matrix
	Matrix massmatrix, Hmatrix; //M = int(rho*((S)^T).S, dV,V); H = int(S,dV,V)	
	//integration points and weights
	Vector x1,w1;

	// node numbers (global)
	int n1;	//$EDC$[varaccess,EDCvarname="node_number1",EDCfolder="Geometry",tooltiptext="global number of node 1 (left), node must already exist"] 
	int n2; //$EDC$[varaccess,EDCvarname="node_number2",EDCfolder="Geometry",tooltiptext="global number of node 2 (right), node must already exist"] 

	//update directors in post newton step, s.t. singularities during the deformation process are avoided.
	int do_update_directors;        //$EDC$[varaccess,EDCvarname="update_directors",EDCfolder="Geometry",int_bool,tooltiptext="update directors during calculation"]

	//kinematic_computation_mode =   0 ... exact terms, 5th order gaussian integration (slow);
	//                               1 ... exact terms, low order (1st order lobatto) integration (fast);
	//                               2 ... approximate mass matrix (torsional terms approximated), no quadratic velocity vector (fastest) --- see also paper by dmitrochenko
	// must be zero (exact terms + 5th order gaussian integration) for the user (models-cpp programmers: use only, if you precisely know what you are doing)
	int kinematic_computation_mode; //$EDC$[varaccess,EDCvarname="kinematic_computation_mode",EDCfolder="Computation",tooltiptext="0 .. exact kinematic terms + 5th order gaussian integration (slow), 1 .. exact terms + 1st order lobatto integration (fast), 2 .. constant mass matrix approximation (fastest)"]
	
	//order of numerical integration along beam axis
	int intorder_mass;          //$EDC$[varaccess,EDCvarname="mass",EDCfolder="Computation.IntegrationOrder",tooltiptext="integration order for mass terms"]
	int intorder_axial_strain;  //$EDC$[varaccess,EDCvarname="axial_strain",EDCfolder="Computation.IntegrationOrder",tooltiptext="integration order for work of axial strain"]
	int intorder_curvature;     //$EDC$[varaccess,EDCvarname="curvature",EDCfolder="Computation.IntegrationOrder",tooltiptext="integration order for work of curvature"]
	
	//reference configuration
	Vector q0;

	//EDC int materialnum; //$EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number which contains the main material properties of the beam"] 
  //EDC double size.X(); //$EDC$[varaccess,readonly,EDCvarname="beam_length",EDCfolder="Info",tooltiptext="length of the beam, calculated by numerical integration of |r'|"] 
}; //$EDC$[endclass,ANCFBeam3DTorsion]



#endif
