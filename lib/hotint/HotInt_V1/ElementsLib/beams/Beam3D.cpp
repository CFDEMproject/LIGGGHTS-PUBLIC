//#**************************************************************
//#
//# filename:             Beam3D.cpp
//#
//# author:               Saxinger Martin, Gerstmayr Johannes
//#
//# generated:						10. Juli 2012
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
 
#include "body3d.h"
#include "Beam3D.h"
#include "elementDataAccess.h"
#include "Material.h"
#include "node.h"
#include "femathHelperFunctions.h"
#include "graphicsConstants.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D Beam3D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Beam3D::SetBeam3D(Vector3D p0i, Matrix3D rot0i, int n1i, int n2i, int materialnumi, 
                   		 Vector3D si, int axialdeformationI, Vector3D coli, Vector q0i)
{
	n1=n1i; // node numbers
	n2=n2i;

	axialdeformation = axialdeformationI; //axial deformation
	
	size = si; // size of beam
	
	materialnum = materialnumi;
	
	col = coli;

	//vector of dofs (reference configuration)
	q0 = q0i;				// only used for drawing
	rot0 = rot0i;   // rigid body deformation
	p0 = p0i;       // rigid body deformation

	x_init = Vector(SS()); //Position & Velocity initial conditions initialized with zeros
	Vector x_init_1 = GetNode(1).X_Init();
	Vector x_init_2 = GetNode(2).X_Init();
	
	if (x_init_1.Length() != SS()*.5 || x_init_2.Length() != SS()*.5) //if node has no initialization, assume zeros
	{
	}
	else
	{
		int node_sos = x_init_1.Length()/2;

		// $ MSax 2013-07-29 : added new code
		x_init.Copy(x_init_1,1,1,node_sos);
		x_init.Copy(x_init_2,1,1+node_sos,node_sos);
		x_init.Copy(x_init_1,1+node_sos,1+2*node_sos,node_sos);
		x_init.Copy(x_init_2,1+node_sos,1+3*node_sos,node_sos);
	}
}

void Beam3D::ElementDefaultConstructorInitialization()
{	
	elementname = GetElementSpec();
	axialdeformation = 1;    // 1|(0) (no) axial deformation
	useAllDOF=1;

	n1=1;
	n2=2;
	size.X()=1;
	size.Y()=0.1;
	size.Z()=0.1;
	materialnum=1;
	q0 = Vector(12);
	rot0 = Matrix3D(1.);
	p0 = Vector3D(0.);

	beam3d_xg_v_list[0] = beam3d_vN1; beam3d_xg_v_list[1] = beam3d_phiZN1; beam3d_xg_v_list[2] = beam3d_vN2; beam3d_xg_v_list[3] = beam3d_phiZN2;
	beam3d_xg_w_list[0] = beam3d_wN1; beam3d_xg_w_list[1] = beam3d_phiYN1; beam3d_xg_w_list[2] = beam3d_wN2; beam3d_xg_w_list[3] = beam3d_phiYN2;
	beam3d_xg_theta_list[0] = beam3d_phiXN1; beam3d_xg_theta_list[1] = beam3d_phiXN2;
	beam3d_xg_u_list[0] = beam3d_uN1; beam3d_xg_u_list[1] = beam3d_uN2;
	for (int i=0; i<NSPos(); i++) 
	{
		beam3d_xg_list[i] = beam3d_xg_v_list[i];
		beam3d_xg_list[i+NSPos()] = beam3d_xg_w_list[i];
	}
	for (int i=0; i<NSTor(); i++)
	{
		beam3d_xg_list[i+2*NSPos()] = beam3d_xg_theta_list[i];
	}
	for (int i=0; i<NSAx(); i++)
	{
		beam3d_xg_list[i+2*NSPos()+NSTor()] = beam3d_xg_u_list[i];
	}
	for (int i=0; i<FlexDOF(); i++)
	{
		for (int j=0; j<FlexDOF(); j++)
		{
			if (beam3d_xg_list[j] == i+1) beam3d_xg_inverse_list[i] = j+1;
		}
	}
}

void Beam3D::CopyFrom(const Element& e)
{
	Body3D::CopyFrom(e);
	const Beam3D& ce = (const Beam3D&)e;
	size = ce.size;

	p0 = ce.p0;
	rot0 = ce.rot0;
	q0 = ce.q0;

	axialdeformation = ce.axialdeformation;
	useAllDOF = ce.useAllDOF;


	//integration points
	x1 = ce.x1;
	w1 = ce.w1;

	n1 = ce.n1;
	n2 = ce.n2;

	for (int i=0; i<NSPos(); i++)
	{
		beam3d_xg_v_list[i] = ce.beam3d_xg_v_list[i];
		beam3d_xg_w_list[i] = ce.beam3d_xg_w_list[i];
	}
	for (int i=0; i<NSTor(); i++)
	{
		beam3d_xg_theta_list[i] = ce.beam3d_xg_theta_list[i];
	}
	for (int i=0; i<NSAx(); i++)
	{
		beam3d_xg_u_list[i] = ce.beam3d_xg_u_list[i];
	}
	for (int i=0; i<FlexDOF(); i++)
	{
		beam3d_xg_list[i] = ce.beam3d_xg_list[i];
	}
	for (int i=0; i<FlexDOF(); i++)
	{
		beam3d_xg_inverse_list[i] = ce.beam3d_xg_inverse_list[i];
	}
}

void Beam3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);
}

int Beam3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);
	SetBeam3D(n1, n2, materialnum);
	x_init = Vector(SS()); //Position & Velocity initial conditions initialized with zeros
	return rv;
}

// standard set function with oriented nodes (default set function for script language)
void Beam3D::SetBeam3D(int n1nr, int n2nr, int matnr)
{
	Node3DRxyz& n1 = (Node3DRxyz&)mbs->GetNode(n1nr);
	Node3DRxyz& n2 = (Node3DRxyz&)mbs->GetNode(n2nr);
	assert((n1.GetTypeName()).CStrCompare("Node3DRxyz") && (n2.GetTypeName()).CStrCompare("Node3DRxyz"));

	Matrix3D f1 = n1.GetLocalFrame();
	Matrix3D f2 = n2.GetLocalFrame();

	if ((f1-f2).AbsNorm() > 1e-14)
	{
		UO(UO_LVL_warn) << "Warning LinearBeam3D: Local frame of node 1 does not coincide with local frame of node 2.";
	}

	//assert((f1-f2).AbsNorm() < 1e-14); // $ MSax 2013-04-15 : removed
	
	rot0 = f1;
	p0 = (n1.Pos()+n2.Pos())*0.5;

	size.X() = (n1.Pos()-n2.Pos()).Norm(); //CalculateElementLength(n1,n2);
	Material mat = GetMaterial();
	if (mat.IsMaterialOfBeamWithRectangularCrossSection())
	{
		size.Y() = GetMaterial().GetBeamThicknessY();
		size.Z() = GetMaterial().GetBeamThicknessZ();
	}
	else
	{
		assert (GetMaterial().IsMaterialOfBeamWithRectangularCrossSection());   //$ PG 2012-12-18: (as yet) a rectangular cross section (which is specified by the material) is required in ANCFBeam3DTorsion
	}

	SetBeam3D(p0, rot0, n1nr, n2nr, matnr, size, axialdeformation, col);
}

mystr Beam3D::GetElementCoordinatesTexDescription()
{
	mystr str = "1-4 & Second order ODE & Bending in local y - direction. \\newline 1. $v\_i$ : Displacement of node i in local y - direction. \\newline 2. $\\psi\_i$ : Local linearized angle of node i about z - axis. \\newline 3. $v\_j$ : Displacement of node j in local y - direction. \\newline 4. $\\psi\_j$ : Local linearized angle of node j about z - axis.\\\\ \\hline";
	str += "5-8 & Second order ODE & Bending in local w - direction. \\newline 5. $w\_i$ : Displacement of node i in local z - direction. \\newline 6. $\\phi\_i$ : Local linearized angle of node i about y - axis. \\newline 7. $w\_j$ : Displacement of node j in local z - direction. \\newline 8. $\\phi\_j$ : Local linearized angle of node j about y - axis.\\\\ \\hline";
	str += "9-10 & Second order ODE & Torsion about local x - axis. \\newline 9. $\\theta\_i$ : Torsion about local x - axis at position of node i. \\newline 10. $\\theta\_j$ : Torsion about local x - axis at position of node j. \\\\ \\hline";
	str += "11-12 & Second order ODE & Displacement in local x - axis. \\newline 11. $u\_i$ : Displacement of node i in local x - direction. \\newline 12. $u\_j$ : Displacement of node j in local x - direction. \\\\ \\hline";
	str += "13-24 & Second order ODE & Time derivatives of the coordinates 1-12. \\\\ \\hline";
	str += "25 & Algebraic equation & Only used if axial deformation is not allowed. \\newline Constraint equation: $C=u\_i-u\_j=0$. \\newline 25. $\\lambda$ is the Lagrange multiplier for the constraint and equates the internal axial force of the beam. \\\\ \\hline";
	return str;
}

void Beam3D::GetElementNodeImage(mystr& name, mystr& caption)
{
	name = "beam3dCoordinates";
	caption = "Beam3D element as line model.";
}

int Beam3D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Body3D::GetAvailableSpecialValues(available_variables);
	
	available_variables.Add(ReadWriteElementDataVariableType("Internal.acceleration",SOS(),0,0.,mystr("accelerations of the element. range: 1-") + mystr(SOS()) )) ;

	GetAvailableSpecialValuesAuto(available_variables);
	return 0;
}

int Beam3D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine

	int rv = Element::ReadSingleElementData(RWdata);

	if(strcmp(RWdata.variable_name.c_str(), "Internal.acceleration") == 0)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= SOS()) //range check
		{
			RWdata.value = XGPP(RWdata.comp1); return 1; 
		}
		else return -2; 
	}

	if (rv == 1) return 1;

	// manual things to read  

	return ReadSingleElementDataAuto(RWdata);
}

int Beam3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	Node& node1 = mbs->GetNode(n1);
	Node& node2 = mbs->GetNode(n1);

	if((node1.GetTypeName()).CStrCompare("Node3DRxyz") != 1 || (node2.GetTypeName()).CStrCompare("Node3DRxyz") != 1 )
	{
		rv = 2;
		errorstr = mystr("Beam3D: Wrong node type detected!\n");
	}

	if(size.X()*size.Y()*size.Z() == 0.)
	{
		rv = 2;
		errorstr = mystr("Beam3D: size.X(), size.Y() and size.Z() must be greater than zero!\n");
	}

	return rv;
}

int Beam3D::ElementDOFtoNodeDOF(int ind, int nodedim) const
{
	// used for [ui,vi,wi,thetai,phiyi,phizi] sequence of node DOF
	const int nodenum6DOF[] = {1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2};
	const int DOFnum6DOF[] = {1,2,3,4,5,6,1,2,3,4,5,6,7,8,9,10,11,12,7,8,9,10,11,12};
	int node;
	int DOF;

	node = nodenum6DOF[ind-1];
	DOF = DOFnum6DOF[ind-1];

	const Node& node1 = GetMBS()->GetNode(n1);
	const Node& node2 = GetMBS()->GetNode(n2);

	if (node == 1) {
		return node1.Get(DOF); // returns global index of local DOF
	} else if (node == 2) {
		return node2.Get(DOF);
	} else return 0;
}

void Beam3D::LinkToElements()
{
	//Vector of degrees of freedom: [ vi, vi', vj, vj', wi, wi', wj, wj', thetai, thetaj, (ui), (uj) ]
	//= 12 DOFs
	//i=left node, j=right node
	//v,w,u = position, v',w',u' = slope, theta = rotation

	TArray<int> storeltg(IS()); // store old ltg list temporary
	for (int i=1; i <= IS(); i++)
	{
		storeltg.Add(LTG(i));
	}
	LTGreset(); // reset ltg list

	int nodedim;
	//UO() << "Link nodes to elements";
	const Node& node1 = GetMBS()->GetNode(n1);
	const Node& node2 = GetMBS()->GetNode(n2);
	
	nodedim = 6;

	// positions:
	for (int i=1; i<=12; i++)
	{
		AddLTG(ElementDOFtoNodeDOF(i, nodedim));
	}

	//same for velocities:
	for (int i=1+nodedim*2; i<=12+nodedim*2; i++)
	{
		AddLTG(ElementDOFtoNodeDOF(i,nodedim));
	}

	// remaining implicit size dofs
	for (int i=1; i <= IS(); i++)
	{
		AddLTG(storeltg(i));
	}
}

//$ MSax 2013-4-8:[
const double Beam3D::XGL(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D disp_1 = GetRefRotN1().GetTp()*Vector3D(XG(1),XG(2),XG(3));
		return disp_1(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D rot_1 = GetRefRotN1().GetTp()*Vector3D(XG(4),XG(5),XG(6));
		return rot_1(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D disp_2 = GetRefRotN2().GetTp()*Vector3D(XG(7),XG(8),XG(9));
		return disp_2(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D rot_2 = GetRefRotN2().GetTp()*Vector3D(XG(10),XG(11),XG(12));
		return rot_2(iloc-9);
	}
}

const double Beam3D::XGPL(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D disp_1 = GetRefRotN1().GetTp()*Vector3D(XGP(1),XGP(2),XGP(3));
		return disp_1(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D rot_1 = GetRefRotN1().GetTp()*Vector3D(XGP(4),XGP(5),XGP(6));
		return rot_1(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D disp_2 = GetRefRotN2().GetTp()*Vector3D(XGP(7),XGP(8),XGP(9));
		return disp_2(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D rot_2 = GetRefRotN2().GetTp()*Vector3D(XGP(10),XGP(11),XGP(12));
		return rot_2(iloc-9);
	}
}

const double Beam3D::XGDL(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D disp_1 = GetRefRotN1().GetTp()*Vector3D(XGD(1),XGD(2),XGD(3));
		return disp_1(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D rot_1 = GetRefRotN1().GetTp()*Vector3D(XGD(4),XGD(5),XGD(6));
		return rot_1(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D disp_2 = GetRefRotN2().GetTp()*Vector3D(XGD(7),XGD(8),XGD(9));
		return disp_2(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D rot_2 = GetRefRotN2().GetTp()*Vector3D(XGD(10),XGD(11),XGD(12));
		return rot_2(iloc-9);
	}
}

const double Beam3D::XGPDL(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D disp_1 = GetRefRotN1().GetTp()*Vector3D(XGPD(1),XGPD(2),XGPD(3));
		return disp_1(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D rot_1 = GetRefRotN1().GetTp()*Vector3D(XGPD(4),XGPD(5),XGPD(6));
		return rot_1(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D disp_2 = GetRefRotN2().GetTp()*Vector3D(XGPD(7),XGPD(8),XGPD(9));
		return disp_2(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D rot_2 = GetRefRotN2().GetTp()*Vector3D(XGPD(10),XGPD(11),XGPD(12));
		return rot_2(iloc-9);
	}
}

Matrix3D Beam3D::GetRefRotN1() const
{
	Node3DRxyz& node_1 = (Node3DRxyz&)mbs->GetNode(n1);
	if(UseGlobalDOF())
	{
		return node_1.GetLocalFrame();
	}
	return Matrix3D(1.);
}
Matrix3D Beam3D::GetRefRotN2() const
{
	// $ MSax 2013-04-15 changed: Only rotation of node 1 is used. This should be changed in future version of beam
	//Node3DRxyz& node_2 = (Node3DRxyz&)mbs->GetNode(n2);
	//if(UseGlobalDOF())
	//{
	//	return node_2.GetLocalFrame();
	//}
	//return Matrix3D(1.);
	return GetRefRotN1();
}

//$ MSax 2013-4-8:]

// -1<=x0<=1
// cubic interpolation
double Beam3D::GetS0w(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return 1.0/2.0-3.0/4.0*x0+Cub(x0)/4.0;							//q1=w (x0=-1)
	case 2: return -size.X()*(1.0-x0-Sqr(x0)+Cub(x0))/8.0;			//q2=w'(x0=-1)
	case 3: return 1.0/2.0+3.0/4.0*x0-Cub(x0)/4.0;							//q3=w (x0= 1)
	case 4: return -size.X()*(1.0+x0)*(1.0+x0)*(-1.0+x0)/8.0;	//q4=w'(x0= 1)
	default: return 0;
	}
}

double Beam3D::GetS0wx(double x0, int shape) const //d(S0w)/dx0 * dx0/dx
{
	switch(shape)
	{
	case 1: return (2./size.X()) * (-3.0/4.0+3.0/4.0*x0*x0);
	case 2: return -(2./size.X()) * (size.X()*(-1.0-2.0*x0+3.0*x0*x0)/8.0);			
	case 3: return (2./size.X()) * (3.0/4.0-3.0/4.0*x0*x0);							
	case 4: return -(2./size.X()) * (size.X()*(1.0+x0)*(-1.0+x0)/4.0+size.X()*(1.0+x0)*(1.0+x0)/8.0);	
	default: return 0;
	}
}

double Beam3D::GetS0wxx(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return (4./Sqr(size.X())) * (3.0/2.0*x0          );
	case 2: return -(4./Sqr(size.X())) * (size.X()*(-1.0+3.0*x0)/4.0);
	case 3: return (4./Sqr(size.X())) * (-3.0/2.0*x0         );
	case 4: return -(4./Sqr(size.X())) * (size.X()*( 1.0+3.0*x0)/4.0);
	default: return 0;
	}
}

double Beam3D::GetS0wxxx(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return (8./Cub(size.X())) * (3.0/2.0   );							
	case 2: return -(8./Cub(size.X())) * (size.X()*3.0/4.0);			
	case 3: return (8./Cub(size.X())) * (-3.0/2.0  );							
	case 4: return -(8./Cub(size.X())) * (size.X()*3.0/4.0);
	default: return 0;
	}
}

double Beam3D::GetS0v(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return GetS0w(x0, 1);						
	case 2: return -GetS0w(x0, 2);			
	case 3: return GetS0w(x0, 3);						
	case 4: return -GetS0w(x0, 4);
	default: return 0;
	}
}

double Beam3D::GetS0vx(double x0, int shape) const //d(S0w)/dx0 * dx0/dx
{
	switch(shape)
	{
	case 1: return GetS0wx(x0, 1);
	case 2: return -GetS0wx(x0, 2);		
	case 3: return GetS0wx(x0, 3);			
	case 4: return -GetS0wx(x0, 4);
	default: return 0;
	}
}

double Beam3D::GetS0vxx(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return GetS0wxx(x0, 1);
	case 2: return -GetS0wxx(x0, 2);
	case 3: return GetS0wxx(x0, 3);
	case 4: return -GetS0wxx(x0, 4);
	default: return 0;
	}
}

double Beam3D::GetS0vxxx(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return GetS0wxxx(x0, 1);						
	case 2: return -GetS0wxxx(x0, 2);		
	case 3: return GetS0wxxx(x0, 3);							
	case 4: return -GetS0wxxx(x0, 4);
	default: return 0;
	}
}

double Beam3D::GetS0u(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return 0.5*(1.-x0);	//q1=u(x0=-1)
	case 2: return 0.5*(1.+x0);	//q2=u(x0= 1)
	default: return 0;
	}
}

double Beam3D::GetS0ux(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return (2./size.X()) * (-0.5);	
	case 2: return (2./size.X()) * ( 0.5);	
	default: return 0;
	}
}

double Beam3D::GetS0theta(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return 0.5*(1.-x0);	//q1=theta(x0=-1)
	case 2: return 0.5*(1.+x0);	//q2=theta(x0= 1)
	default: return 0;
	}
}

double Beam3D::GetS0thetax(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return (2./size.X()) * (-0.5);	
	case 2: return (2./size.X()) * ( 0.5);	
	default: return 0;
	}
}

double Beam3D::GetS0beta(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return -0.5*(1.-x0)*x0;	  //q1=beta(x0=-1)
	case 2: return  (1.-x0)*(x0-1.);	//q2=beta(x0= 0)
	case 3: return 0.5*(1.+x0)*x0;	  //q3=beta(x0= 1)
	default: return 0;
	}
}

double Beam3D::GetS0betax(double x0, int shape) const
{
	switch(shape)
	{
	case 1: return (2./size.X()) * (x0-0.5  );	     
	case 2: return (2./size.X()) * (2.-2.*x0);	
	case 3: return (2./size.X()) * (x0+0.5  );	
	default: return 0;
	}
}

double Beam3D::GetU(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*XGL(beam3d_xg_u_list[i-1]);
	return u;
}


double Beam3D::GetV(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*XGL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetW(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0w(x0,i)*XGL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetUx(double x0) const //d(U)/dx
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0ux(x0,i)*XGL(beam3d_xg_v_list[i-1]);
	return u;
}

double Beam3D::GetVx(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*XGL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWx(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wx(x0,i)*XGL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetVxx(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vxx(x0,i)*XGL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWxx(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wxx(x0,i)*XGL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetVxxx(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vxxx(x0,i)*XGL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWxxx(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wxxx(x0,i)*XGL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetTheta(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*XGL(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetThetax(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0thetax(x0,i)*XGL(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetAngleX(double x0) const //(linearized) rotation about local x-axis
{
	return GetTheta(x0);
}

double Beam3D::GetAngleY(double x0) const //(linearized) rotation about local y-axis
{
	return -GetWx(x0);
}

double Beam3D::GetAngleZ(double x0) const //(linearized) rotation about local z-axis
{
	return GetVx(x0);
}

double Beam3D::GetdAngleXdx(double x0) const
{
	return GetThetax(x0);
}

double Beam3D::GetdAngleYdx(double x0) const
{
	return -GetWxx(x0);
}

double Beam3D::GetdAngleZdx(double x0) const
{
	return GetVxx(x0);
}

double Beam3D::GetUP(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*XGPL(beam3d_xg_u_list[i-1]);
	return u;
}

double Beam3D::GetVP(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*XGPL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWP(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0w(x0,i)*XGPL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetUxP(double x0) const //d(U)/dx
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0ux(x0,i)*XGPL(beam3d_xg_u_list[i-1]);
	return u;
}

double Beam3D::GetVxP(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*XGPL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWxP(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wx(x0,i)*XGPL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetVxxP(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vxx(x0,i)*XGPL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWxxP(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wxx(x0,i)*XGPL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetThetaP(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*XGPL(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetThetaxP(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0thetax(x0,i)*XGPL(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetAngleXP(double x0) const //(linearized) rotation about local x-axis
{
	return GetThetaP(x0);
}

double Beam3D::GetAngleYP(double x0) const //(linearized) rotation about local y-axis
{
	return -GetWxP(x0);
}

double Beam3D::GetAngleZP(double x0) const //(linearized) rotation about local z-axis
{
	return GetVxP(x0);
}

const double Beam3D::XGPPL(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D acc_1 = GetRefRotN1().GetTp()*Vector3D(XGPP(1),XGPP(2),XGPP(3));  // $ MSax 2013-07-16 : changed from GetMBS()->GetAcceleration(i+SOS()) to XGPP(i)
		return acc_1(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D acc_1 = GetRefRotN1().GetTp()*Vector3D(XGPP(4),XGPP(5),XGPP(6));  // $ MSax 2013-07-16 : changed from GetMBS()->GetAcceleration(i+SOS()) to XGPP(i)
		return acc_1(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D acc_2 = GetRefRotN2().GetTp()*Vector3D(XGPP(7),XGPP(8),XGPP(9));  // $ MSax 2013-07-16 : changed from GetMBS()->GetAcceleration(i+SOS()) to XGPP(i)
		return acc_2(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D acc_2 = GetRefRotN2().GetTp()*Vector3D(XGPP(10),XGPP(11),XGPP(12));  // $ MSax 2013-07-16 : changed from GetMBS()->GetAcceleration(i+SOS()) to XGPP(i)
		return acc_2(iloc-9);
	}
}

const double Beam3D::q0l(int iloc) const
{
	if(iloc >= 1 && iloc <=3)
	{
		Vector3D vec = GetRefRotN1().GetTp()*Vector3D(q0(1),q0(2),q0(3));
		return vec(iloc);
	} 
	else if (iloc >= 4 && iloc <=6)
	{
		Vector3D vec = GetRefRotN1().GetTp()*Vector3D(q0(4),q0(5),q0(6));
		return vec(iloc-3);
	} 
	else if (iloc >= 7 && iloc <=9)
	{
		Vector3D vec = GetRefRotN2().GetTp()*Vector3D(q0(7),q0(8),q0(9));
		return vec(iloc-6);
	} 
	else if (iloc >= 10 && iloc <=12)
	{
		Vector3D vec = GetRefRotN2().GetTp()*Vector3D(q0(10),q0(11),q0(12));
		return vec(iloc-9);
	}
}

double Beam3D::GetUPP(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*XGPPL(beam3d_xg_u_list[i-1]);
	return u;
}

double Beam3D::GetVPP(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*XGPPL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWPP(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0w(x0,i)*XGPPL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetAngleXPP(double x0) const //(linearized) rotation about local x-axis
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*XGPPL(beam3d_xg_theta_list[i-1]);
	return theta;}

double Beam3D::GetAngleYPP(double x0) const //(linearized) rotation about local y-axis
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += -1*GetS0wx(x0,i)*XGPPL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetAngleZPP(double x0) const //(linearized) rotation about local z-axis
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*XGPPL(beam3d_xg_v_list[i-1]);
	return v;
}

void Beam3D::GetRotPP(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleXPP(x0),GetAngleYPP(x0),GetAngleZPP(x0));
	rot = rot0*rot;
}

double Beam3D::GetUD(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*(XGDL(beam3d_xg_u_list[i-1])+q0l(beam3d_xg_u_list[i-1]));
	return u;
}

double Beam3D::GetVD(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*(XGDL(beam3d_xg_v_list[i-1])+q0l(beam3d_xg_v_list[i-1]));
	return v;
}

double Beam3D::GetWD(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++) 
		w += GetS0w(x0,i)*(XGDL(beam3d_xg_w_list[i-1])+q0l(beam3d_xg_w_list[i-1]));
	return w;
}

double Beam3D::GetUPD(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*XGPDL(beam3d_xg_u_list[i-1]);
	return u;
}

double Beam3D::GetVPD(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*XGPDL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetWPD(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0w(x0,i)*XGPDL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetAngleXPD(double x0) const //(linearized) rotation about local x-axis
{
	return GetThetaPD(x0);
}

double Beam3D::GetAngleYPD(double x0) const //(linearized) rotation about local y-axis
{
	return -GetWxPD(x0);
}

double Beam3D::GetAngleZPD(double x0) const //(linearized) rotation about local z-axis
{
	return GetVxPD(x0);
}

double Beam3D::GetThetaPD(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*XGPDL(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetWxPD(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wx(x0,i)*XGPDL(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetVxPD(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*XGPDL(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetUxD(double x0) const //d(U)/dx
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0ux(x0,i)*(XGDL(beam3d_xg_u_list[i-1])+q0l(beam3d_xg_u_list[i-1]));
	return u;
}

double Beam3D::GetVxD(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*(XGDL(beam3d_xg_v_list[i-1])+q0l(beam3d_xg_v_list[i-1]));
	return v;
}

double Beam3D::GetWxD(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wx(x0,i)*(XGDL(beam3d_xg_w_list[i-1])+q0l(beam3d_xg_w_list[i-1]));
	return w;
}

double Beam3D::GetVxxD(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vxx(x0,i)*(XGDL(beam3d_xg_v_list[i-1])+q0l(beam3d_xg_v_list[i-1]));
	return v;
}

double Beam3D::GetWxxD(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wxx(x0,i)*(XGDL(beam3d_xg_w_list[i-1])+q0l(beam3d_xg_w_list[i-1]));
	return w;
}

double Beam3D::GetVxxxD(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vxxx(x0,i)*(XGDL(beam3d_xg_v_list[i-1])+q0l(beam3d_xg_v_list[i-1]));
	return v;
}

double Beam3D::GetWxxxD(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wxxx(x0,i)*(XGDL(beam3d_xg_w_list[i-1])+q0l(beam3d_xg_w_list[i-1]));
	return w;
}

double Beam3D::GetThetaD(double x0) const 
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*(XGDL(beam3d_xg_theta_list[i-1])+q0l(beam3d_xg_theta_list[i-1]));
	return theta;
}

double Beam3D::GetThetaxD(double x0) const 
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0thetax(x0,i)*(XGDL(beam3d_xg_theta_list[i-1])+q0l(beam3d_xg_theta_list[i-1]));
	return theta;
}

double Beam3D::GetAngleXD(double x0) const //(linearized) rotation about local x-axis, in beam-coordinates
{
	return GetThetaD(x0);
}

double Beam3D::GetAngleYD(double x0) const //(linearized) rotation about local y-axis, in beam-coordinates
{
	return -GetWxD(x0);
}

double Beam3D::GetAngleZD(double x0) const //(linearized) rotation about local z-axis, in beam-coordinates
{
	return GetVxD(x0);
}

double Beam3D::GetAngleXInit(double x0) const //(linearized) rotation about local x-axis, in beam-coordinates
{
	return GetThetaInit(x0);
}

double Beam3D::GetAngleYInit(double x0) const //(linearized) rotation about local y-axis, in beam-coordinates
{
	return -GetWxInit(x0);
}

double Beam3D::GetAngleZInit(double x0) const //(linearized) rotation about local z-axis, in beam-coordinates
{
	return GetVxInit(x0);
}

double Beam3D::GetWxInit(double x0) const //d(W)/dx
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0wx(x0,i)*q0l(beam3d_xg_w_list[i-1]);
	return w;
}

double Beam3D::GetVxInit(double x0) const //d(V)/dx
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0vx(x0,i)*q0l(beam3d_xg_v_list[i-1]);
	return v;
}

double Beam3D::GetThetaInit(double x0) const //d(W)/dx
{
	double theta = 0;
	for (int i=1; i <= NSTor(); i++)
		theta += GetS0theta(x0,i)*q0l(beam3d_xg_theta_list[i-1]);
	return theta;
}

double Beam3D::GetUInit(double x0) const  //axial deformation
{
	double u = 0;
	for (int i=1; i <= NSAx()*useAllDOF; i++)
		u += GetS0u(x0,i)*q0l(beam3d_xg_u_list[i-1]);
	return u;
}
double Beam3D::GetVInit(double x0) const  //deflection in Y-direction
{
	double v = 0;
	for (int i=1; i <= 4; i++)
		v += GetS0v(x0,i)*q0l(beam3d_xg_v_list[i-1]);
	return v;
}
double Beam3D::GetWInit(double x0) const	 //deflection in Z-direction
{
	double w = 0;
	for (int i=1; i <= 4; i++)
		w += GetS0w(x0,i)*q0l(beam3d_xg_w_list[i-1]);
	return w;
}

void Beam3D::GetRot(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleX(x0),GetAngleY(x0),GetAngleZ(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	rot = rot0*rot;
}

void Beam3D::GetdRotdx(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetdAngleXdx(x0),GetdAngleYdx(x0),GetdAngleZdx(x0));
	rot = rot0*rot;
}

void Beam3D::GetRotP(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleXP(x0),GetAngleYP(x0),GetAngleZP(x0));
	rot = rot0*rot;
}

// $ MSax 2012-12-14
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Beam3D::GetRotPD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleXPD(x0),GetAngleYPD(x0),GetAngleZPD(x0));
	rot = rot0*rot;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Beam3D::GetRotD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleXD(x0),GetAngleYD(x0),GetAngleZD(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	rot = rot0*rot;
}

void Beam3D::GetRotInit(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(GetAngleXInit(x0),GetAngleYInit(x0),GetAngleZInit(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	rot = rot0*rot;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Beam3D::GetdPosdx(const Vector3D& ploc, Vector3D& dpdx)
{
	double x0 = ploc.X()/(0.5*size.X());  //x0 is normalized to -1..+1
	Matrix3D rotx;
	GetdRotdx(ploc.X(), rotx);
	dpdx = rot0*Vector3D(1+GetUx(x0),GetVx(x0),GetWx(x0))+rotx*Vector3D(0.,ploc.Y(),ploc.Z());
};

Vector3D Beam3D::GetNodeLocPos(int i) const //local position of node $ SW + MSax 2013-08-27: added
{
	if(i>2 || i<1)
	{
		GetMBS()->UO(UO_LVL_warn) << "WARNING: Beam3D: GetNodeLocPos: wrong local node number detected.";
		return Vector3D(0.);
	}
	
	Node3DRxyz& no1 = (Node3DRxyz&)mbs->GetNode(n1);
	Node3DRxyz& no2 = (Node3DRxyz&)mbs->GetNode(n2);
	double halfLength = (no1.Pos()-no2.Pos()).Norm()/2;
	if(i==1)
	{
		return Vector3D(-halfLength,0,0);
	}
	else
	{
		return Vector3D(halfLength,0,0);
	}
}

Vector3D Beam3D::GetPos(double x) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	return p0 + rot0*Vector3D(x+GetU(x0),GetV(x0),GetW(x0));
};

//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D Beam3D::GetPos(const Vector3D& p_loc) const
{
	double x0 = p_loc.X()*2./size.X();
	Matrix3D rot;
	GetRot(p_loc.X(),rot);
	return p0 + rot0*Vector3D(p_loc.X()+GetU(x0),GetV(x0),GetW(x0))
		 + rot*Vector3D(0,p_loc.Y(),p_loc.Z());
};

//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D Beam3D::GetVel(const Vector3D& p_loc) const
{
	double x0 = p_loc.X()*2./size.X();
	Matrix3D rotp;
	GetRotP(p_loc.X(),rotp);
	
	return rot0*Vector3D(GetUP(x0),GetVP(x0),GetWP(x0)) 
		+ rotp*Vector3D(0,p_loc.Y(),p_loc.Z());
};

//compute acceleration vector (after Timestep is successfully computed)
Vector3D Beam3D::GetAcceleration(const Vector3D& p_loc) const //$ MS 2012-11-5: added for use in MBSSensor.
{
	double x0 = p_loc.X()*2./size.X();
	Matrix3D rotpp;
	GetRotPP(p_loc.X(),rotpp);

	return rot0*Vector3D(GetUPP(x0),GetVPP(x0),GetWPP(x0)) + rotpp*Vector3D(0,p_loc.Y(),p_loc.Z());
}

//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D Beam3D::GetPosD(const Vector3D& p_loc) const
{
	double x0 = p_loc.X()*2./size.X();
	Matrix3D rotD;
	GetRotD(p_loc.X(),rotD);
	Matrix3D rotInit;
	GetRotInit(p_loc.X(),rotInit);
	Vector3D unscaled_position = p0 + rot0*Vector3D(p_loc.X()+GetUD(x0),GetVD(x0),GetWD(x0)) + rotD*Vector3D(0,p_loc.Y(),p_loc.Z());
	Vector3D initial_position = p0 + rot0*Vector3D(p_loc.X()+GetUInit(x0),GetVInit(x0),GetWInit(x0)) + rotInit*Vector3D(0,p_loc.Y(),p_loc.Z());
	return initial_position + GetMBS()->GetDOption(105)*(unscaled_position - initial_position);
};

// $ MSax 2012-12-14: added
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D Beam3D::GetVelD(const Vector3D& p_loc) const
{
		double x0 = p_loc.X()*2./size.X();
		Matrix3D rotp;
		GetRotPD(p_loc.X(),rotp);
		
		return rot0*Vector3D(GetUPD(x0),GetVPD(x0),GetWPD(x0)) 
			+ rotp*Vector3D(0,p_loc.Y(),p_loc.Z());
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// v0,v'0,vL,v'L (y-direction)
	// w0,w'0,wL,w'L (z-direction)
	// theta0, thetaL (torsion; rotation about x-axis)
	//optional:
	// u0, uL (axial deformation; optional)
//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D Beam3D::GetDOFPosD(int idof) const
{
	// $ MSax 2013-03-04: changed for new dof order

	if (idof == beam3d_vN1 || idof == beam3d_phiZN1 || idof == beam3d_wN1 || idof == beam3d_phiYN1 || idof == beam3d_phiXN1 || idof == beam3d_uN1) // node 1
	{
		return GetPosD(Vector3D(-0.5*size.X(),0.,0.));
	}
	else // node 2
	{
		return GetPosD(Vector3D(0.5*size.X(),0.,0.));
	}
};

Vector3D Beam3D::GetDOFDirD(int idof) const
{
	// $ MSax 2013-03-04: changed for new dof order

	if (idof == beam3d_vN1 || idof == beam3d_vN2) // vi or vj
	{
		return rot0*Vector3D(0.,1.,0.);
	}
	else if (idof == beam3d_wN1 || idof == beam3d_wN2) // wi or wj
	{
		return rot0*Vector3D(0.,0.,1.);
	}
	else if (idof == beam3d_uN1 || idof == beam3d_uN2) // ui or uj
	{
		return rot0*Vector3D(1.,0.,0.);
	}
	else if (idof == beam3d_phiZN1 || idof == beam3d_phiZN2) // vi' or vj'
	{
		return rot0*Vector3D(0.,0.,2.);
	}
	else if (idof == beam3d_phiYN1 || idof == beam3d_phiYN2) // wi' or wj'
	{
		return rot0*Vector3D(0.,2.,0.);
	}
	else if (idof == beam3d_phiXN1 || idof == beam3d_phiXN2) // thetai or thetaj
	{
		return rot0*Vector3D(2.,0.,0.);
	}
	return Vector3D(0.,0.,0.);
};

void Beam3D::GetH(Matrix& H) 
{
		int dim = Dim();

		H.SetSize(SOS(),dim);
		H.SetAll(0);

		GetIntegrationRule(x1,w1,7); //3x1x1 !!!!!

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			double x0 = x1(i1);
			double fact = w1(i1)*size.X()/2.*GetMaterial().BeamRhoA()/GetMaterial().Density(); // changed: SM 05122012 (integrate du/dq*dV)

			for (int i=1; i <= NSAx()*useAllDOF; i++) //U
				H(2*(4)+NSTor()+i,1) += fact*GetS0u(x0,i);

			for (int i=1; i <= 4; i++) //V
				H(i,2) += fact*GetS0v(x0,i);

			for (int i=1; i <= 4; i++) //W
				H(i+4,3) += fact*GetS0w(x0,i);
		}
}

void Beam3D::CalcMK_Matrices()
{
		double fakt;

		Kv = Matrix(4,4);
		Kw = Matrix(4,4);
		Ktheta = Matrix(2,2);
		Ku = Matrix(2,2);
		Mv = Matrix(4,4);
		Mw = Matrix(4,4);
		Mtheta = Matrix(2,2);
		Mu = Matrix(2,2);
		M_with_u = Matrix(12,12);
		M_without_u = Matrix(10,10);

		// Kv
		fakt = (-1.)*GetMaterial().BeamEIz()/Cub(size.X());
		Kv(1,1) = fakt*12.; Kv(2,1) = fakt*6.*size.X(); Kv(3,1) = fakt*(-12.); Kv(4,1) = fakt*6.*size.X(); Kv(1,2) = fakt*6.*size.X(); Kv(2,2) = fakt*4.*size.X()*size.X(); Kv(3,2) = fakt*(-6.)*size.X(); Kv(4,2) = fakt*2.*size.X()*size.X();
		Kv(1,3) = fakt*(-12.); Kv(2,3) = fakt*(-6.)*size.X(); Kv(3,3) = fakt*12.; Kv(4,3) = fakt*(-6.)*size.X(); Kv(1,4) = fakt*6.*size.X(); Kv(2,4) = fakt*2.*size.X()*size.X(); Kv(3,4) = fakt*(-6.)*size.X(); Kv(4,4) = fakt*4.*size.X()*size.X();
		
		// Kw
		fakt = (-1.)*GetMaterial().BeamEIy()/Cub(size.X());
		Kw(1,1) = fakt*12.; Kw(2,1) = fakt*(-6.)*size.X(); Kw(3,1) = fakt*(-12.); Kw(4,1) = fakt*(-6.)*size.X(); 
		Kw(1,2) = fakt*(-6.)*size.X(); Kw(2,2) = fakt*4.*size.X()*size.X(); Kw(3,2) = fakt*6.*size.X(); Kw(4,2) = fakt*2.*size.X()*size.X();
		Kw(1,3) = fakt*(-12.); Kw(2,3) = fakt*6.*size.X(); Kw(3,3) = fakt*12.; Kw(4,3) = fakt*6.*size.X(); 
		Kw(1,4) = fakt*(-6.)*size.X(); Kw(2,4) = fakt*2.*size.X()*size.X(); Kw(3,4) = fakt*6.*size.X(); Kw(4,4) = fakt*4.*size.X()*size.X();

		// Ktheta
		fakt = (-1.)*GetMaterial().BeamGJkx()/size.X();
		Ktheta(1,1) = fakt; Ktheta(2,1) = (-1.)*fakt; Ktheta(1,2) = (-1.)*fakt; Ktheta(2,2) = fakt; 

		// Ku
		fakt = (-1.)*GetMaterial().BeamEA()/size.X();
		Ku(1,1) = fakt; Ku(2,1) = (-1.)*fakt; Ku(1,2) = (-1.)*fakt; Ku(2,2) = fakt;


		// Mv
		fakt = GetMaterial().BeamRhoA()*size.X()/420.;
		Mv(1,1) = fakt*156.; Mv(2,1) = fakt*22.*size.X(); Mv(3,1) = fakt*54.; Mv(4,1) = fakt*(-13.)*size.X(); Mv(1,2) = fakt*22.*size.X(); Mv(2,2) = fakt*4.*size.X()*size.X(); Mv(3,2) = fakt*13.*size.X(); Mv(4,2) = fakt*(-3.)*size.X()*size.X();
		Mv(1,3) = fakt*54.; Mv(2,3) = fakt*13.*size.X(); Mv(3,3) = fakt*156.; Mv(4,3) = fakt*(-22.)*size.X(); Mv(1,4) = fakt*(-13.)*size.X(); Mv(2,4) = fakt*(-3.)*size.X()*size.X(); Mv(3,4) = fakt*(-22.)*size.X(); Mv(4,4) = fakt*4.*size.X()*size.X();

		// Mw
		Mw(1,1) = fakt*156.; Mw(2,1) = fakt*(-22.)*size.X(); Mw(3,1) = fakt*54.; Mw(4,1) = fakt*13.*size.X(); Mw(1,2) = fakt*(-22.)*size.X(); Mw(2,2) = fakt*4.*size.X()*size.X(); Mw(3,2) = fakt*(-13.)*size.X(); Mw(4,2) = fakt*(-3.)*size.X()*size.X();
		Mw(1,3) = fakt*54.; Mw(2,3) = fakt*(-13.)*size.X(); Mw(3,3) = fakt*156.; Mw(4,3) = fakt*22.*size.X(); Mw(1,4) = fakt*13.*size.X(); Mw(2,4) = fakt*(-3.)*size.X()*size.X(); Mw(3,4) = fakt*22.*size.X(); Mw(4,4) = fakt*4.*size.X()*size.X();
		
		// Mtheta
		fakt = GetMaterial().BeamRhoIx()*size.X();
		Mtheta(1,1) = fakt*1./3.; Mtheta(2,1) = fakt*1./6.; Mtheta(1,2) = fakt*1./6.; Mtheta(2,2) = fakt*1./3.; 

		// Mu
		fakt = GetMaterial().BeamRhoA()*size.X();
		Mu(1,1) = fakt*1./3.; Mu(2,1) = fakt*1./6.; Mu(1,2) = fakt*1./6.; Mu(2,2) = fakt*1./3.; 

		M_without_u.InsertMatrix(1,1,Mv);
		M_without_u.InsertMatrix(NSPos()+1,NSPos()+1,Mw);
		M_without_u.InsertMatrix(2*(NSPos())+1,2*(NSPos())+1,Mtheta);

		M_with_u.InsertMatrix(1,1,M_without_u);
		M_with_u.InsertMatrix(2*(NSPos())+NSTor()+1,2*(NSPos())+NSTor()+1,Mu);

		Matrix tr(FlexDOF(),FlexDOF());
		tr.SetDiagMatrix(1);

		TransformationToNodeDOF(tr);
		M_with_u_global_dof = tr*M_with_u*tr.GetTp();
}

void Beam3D::TransformationToNodeDOF(Matrix& d)
{
	int col = d.Getcols();
	int row = d.Getrows();

	Matrix res; 
	res.SetSize(row,col);

	//sort DOF
	for (int i=1; i <= FlexDOF(); i++)
	{
		for (int j=1; j <= col; j++)
		{
			res(i,j) = d(beam3d_xg_inverse_list[i-1],j);
		}
	}

	Matrix tmp;
	tmp.SetSize(3,col);
	Matrix d_part;
	Matrix rot = (Matrix)GetRefRotN1();

	// mult 1
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= col; j++)
		{
			tmp(i,j) = res(i,j);
		}
	}

	Mult(rot,tmp,d_part);
	d.SetSubmatrix(d_part,1,1);

	// mult 2
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= col; j++)
		{
			tmp(i,j) = res(i+3,j);
		}
	}

	Mult(rot,tmp,d_part);
	d.SetSubmatrix(d_part,4,1);

	// mult 3
	rot = (Matrix)GetRefRotN2();

	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= col; j++)
		{
			tmp(i,j) = res(i+6,j);
		}
	}

	Mult(rot,tmp,d_part);
	d.SetSubmatrix(d_part,7,1);

	// mult 4
	for (int i=1; i <= 3; i++)
	{
		for (int j=1; j <= col; j++)
		{
			tmp(i,j) = res(i+9,j);
		}
	}

	Mult(rot,tmp,d_part);
	d.SetSubmatrix(d_part,10,1);
}

void Beam3D::EvalM(Matrix& m, double t) 
{
	m = M_with_u_global_dof;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Beam3D::EvalG(Vector& f, double t) // ui-uj=0 
{
	if (MaxIndex()==3) // displacement
	{
		f(1)=XGL(beam3d_uN1)-XGL(beam3d_uN2);
	}
	else
	{
		f(1)=XGPL(beam3d_uN1)-XGPL(beam3d_uN2);
	}
} 

void Beam3D::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);
	TMStartTimer(22);


	ConstVector<12> fadd(12);

	ConstVector<4> xv(XGL(beam3d_vN1),XGL(beam3d_phiZN1),XGL(beam3d_vN2),XGL(beam3d_phiZN2));
	ConstVector<4> xw(XGL(beam3d_wN1),XGL(beam3d_phiYN1),XGL(beam3d_wN2),XGL(beam3d_phiYN2));
	ConstVector<2> xtheta(XGL(beam3d_phiXN1),XGL(beam3d_phiXN2));
	ConstVector<2> xu(XGL(beam3d_uN1),XGL(beam3d_uN2));

	Vector fv, fw, ftheta, fu;
	fv.LinkWith(fadd,1,4);
	fw.LinkWith(fadd,5,4);
	ftheta.LinkWith(fadd,9,2);
	fu.LinkWith(fadd,11,2);

	Mult(Kv,xv,fv);
	Mult(Kw,xw,fw);
	Mult(Ktheta,xtheta,ftheta);
	Mult(Ku,xu,fu);

	AddEPCqTterms(fadd);

	Matrix tr(FlexDOF(),FlexDOF());
	tr.SetDiagMatrix(1);

	TransformationToNodeDOF(tr);

	f += tr*fadd;

	//Add C_q^T terms
	
	TMStopTimer(22);

}

//C_q^T*\lambda for Euler parameter equation:
void Beam3D::AddEPCqTterms(Vector& f)
{
	if(axialdeformation == 0)
	{
		double lambda = GetLagrangeMultEP();

		f(11) -= lambda;
		f(12) += lambda;
	}
}

void Beam3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class

	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_torsion);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	variables.Add(FieldVariableDescriptor(FieldVariableDescriptor::FVT_beam_force_transversal, FieldVariableDescriptor::FVCI_y));
	variables.Add(FieldVariableDescriptor(FieldVariableDescriptor::FVT_beam_force_transversal, FieldVariableDescriptor::FVCI_z));
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_torsional);
	variables.Add(FieldVariableDescriptor(FieldVariableDescriptor::FVT_beam_moment_bending, FieldVariableDescriptor::FVCI_y));
	variables.Add(FieldVariableDescriptor(FieldVariableDescriptor::FVT_beam_moment_bending, FieldVariableDescriptor::FVCI_z));
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration,FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_bryant_angle, FieldVariableDescriptor::FVCI_z); //$ DR 2013-01-29 added
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity, FieldVariableDescriptor::FVCI_z);  //$ DR 2013-01-29 added
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity_local_basis, FieldVariableDescriptor::FVCI_z);  //$ DR 2013-01-29 added
}

double Beam3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	// local_position is local position in [-size.X()/2., +size.X()/2.] x [-size.Y()/2., +size.Y()/2.] x [-size.Z()/2., +size.Z()/2.]
	double beamEA = GetMaterial().BeamEA();
	double beamEIz = GetMaterial().BeamEIz();
	double beamEIy = GetMaterial().BeamEIy();
	double beamGJkx = GetMaterial().BeamGJkx();
	if(!flagD)
	{
		//ploc is from -1 .. +1
		//x0 is from -1 .. +1
		double x0 = local_position.X()*2./size.X();

		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:	
			return fvd.GetComponent(GetPos(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			{
				Vector3D pr0 = GetInitRefPos();
				return fvd.GetComponent(GetPos(local_position)-pr0-rot0*local_position);
			}
		case FieldVariableDescriptor::FVT_beam_torsion: return (GetAngleX(1)-GetAngleX(-1))/size.X();
		case FieldVariableDescriptor::FVT_beam_force_axial: return beamEA * GetUx(x0);
		case FieldVariableDescriptor::FVT_beam_force_transversal:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
				{
					return -beamEIz * GetVxxx(x0); //Qy=-dMz/dx; see Mang, page 196, p.176, p.178
				}
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
				{
					return -beamEIy * GetWxxx(x0); //Qz= dMy/dx; see Mang, page 196, p.176, p.178

				}
			}
		case FieldVariableDescriptor::FVT_beam_moment_torsional: return beamGJkx * GetThetax(x0);	//Mx = theta_x * GJ *kx
		case FieldVariableDescriptor::FVT_beam_moment_bending:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
					return beamEIy * GetWxx(x0);			//My = -EIy * w''; see Mang, page 196, p.176
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
					return beamEIz * GetVxx(x0);			//Mz =  EIz * v''; see Mang, page 196, p.178
			}
		case FieldVariableDescriptor::FVT_acceleration:
			{
				return fvd.GetComponent(GetAcceleration(local_position));
			}
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2013-01-29 added
			return fvd.GetComponent(GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2013-01-29 added
			return fvd.GetComponent(GetRotMatrix(local_position).GetTp()*GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2013-01-29 added
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrix(local_position), phi);
				return fvd.GetComponent(phi);
			}
		}
	} 
	else 
	{
		//ploc is from -1 .. +1
		//x0 is from -1 .. +1
		double x0 = local_position.X()*2./size.X();

		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:		
			return fvd.GetComponent(GetPosD(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVelD(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			{
				Vector3D pr0 = GetInitRefPos();
				return fvd.GetComponent(GetPosD(local_position)-pr0-rot0*local_position);
			}
		case FieldVariableDescriptor::FVT_beam_force_axial: return beamEA * GetUxD(x0);
		case FieldVariableDescriptor::FVT_beam_torsion: return (GetAngleXD(1)-GetAngleXD(-1))/size.X();
		case FieldVariableDescriptor::FVT_beam_force_transversal:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
				{
					return -beamEIz * GetVxxxD(x0); //Qy=-dMz/dx; see Mang, page 196, p.176, p.178
				}
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
				{
					return -beamEIy * GetWxxxD(x0); //Qz= dMy/dx; see Mang, page 196, p.176, p.178

				}
			}
		case FieldVariableDescriptor::FVT_beam_moment_torsional: return beamGJkx * GetThetaxD(x0);	//Mx = theta_x * GJ *kx
		case FieldVariableDescriptor::FVT_beam_moment_bending:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
					return beamEIy * GetWxxD(x0);			//My = -EIy * w''; see Mang, page 196, p.176
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
					return beamEIz * GetVxxD(x0);			//Mz =  EIz * v''; see Mang, page 196, p.178
			}
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2013-01-29 added
			return fvd.GetComponent(GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2013-01-29 added
			return fvd.GetComponent(GetRotMatrixD(local_position).GetTp()*GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2013-01-29 added
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrixD(local_position), phi);
				return fvd.GetComponent(phi);
			}
		}
	}

	// if not implemented here, call parent class
	//FVT_position, FVT_velocity: implemented in Element
	return Body3D::GetFieldVariableValue(fvd, local_position, flagD);
}

void Beam3D::DrawElement() 
{
	mbs->SetColor(col);

	double lx1 = 0.5*size.X(); double ly1 = 0.5*size.Y()*GetMBS()->GetMagnifyYZ(); double lz1 = 0.5*size.Z()*GetMBS()->GetMagnifyYZ();

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	}
	if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}

	//linemode = 0;
	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double tilex = GetMBS()->GetIOption(137);
	double tiley = 1;//GetMBS()->GetIOption(138); //no tiling in transverse direction

	TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
	TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
	double v=0;

	for (int side=1; side <= 6; side++)
	{
		points.SetLen(0); vals.SetLen(0);
		Vector3D p0, vx, vy;
		int tileyn = (int)tiley;
		int tilexn = (int)tilex;

		switch(side)
		{
		case 1:
		{ //bottom
			p0 = Vector3D(-lx1,-ly1,-lz1);
			vx = Vector3D(2.*lx1/tilexn,0,0);
			vy = Vector3D(0,0,2.*lz1/tileyn); break;
		}
		case 2:
		{ //top
			p0 = Vector3D(-lx1, ly1, lz1);
			vx = Vector3D(2.*lx1/tilexn,0,0);
			vy = Vector3D(0,0,-2.*lz1/tileyn); break;
		}
		case 3:
		{ //front
			p0 = Vector3D(-lx1,-ly1, lz1);
			vx = Vector3D(2.*lx1/tilexn,0,0);
			vy = Vector3D(0,2.*ly1/tileyn,0); break;
		}
		case 4:
		{ //back
			p0 = Vector3D(-lx1, ly1,-lz1);
			vx = Vector3D(2.*lx1/tilexn,0,0);
			vy = Vector3D(0,-2.*ly1/tileyn,0); break;
		}
		case 5:
		{ //left
			p0 = Vector3D(-lx1, ly1,-lz1);
			vx = Vector3D(0,-2.*ly1/tilexn,0);
			vy = Vector3D(0,0,2.*lz1/tileyn); break;
		}
		case 6:
		{ //right
			p0 = Vector3D( lx1,-ly1,-lz1);
			vx = Vector3D(0,2.*ly1/tilexn,0);
			vy = Vector3D(0,0,2.*lz1/tileyn); break;
		}
		default: ;
		}

		for (double iy = 0; iy <= tileyn+1e-10; iy++)
		{
			for (double ix = 0; ix <= tilexn+1e-10; ix++)
			{
				Vector3D ploc = (p0+ix*vx+iy*vy);
				Vector3D p(GetPosD(ploc));
				points.Add(p);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(),ploc, true);

				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
	}

};

void Beam3D::TransformMatrix(Matrix3D rot,Matrix& d) //transform dudq, etc. matrices by means of rot0 reference rotation
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

void Beam3D::GetIntDuDq(Matrix& dudq)
{
	GetH(dudq);
	TransformMatrix(rot0,dudq);
	TransformationToNodeDOF(dudq);
}


void Beam3D::GetdRotvdqTLoc(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	//compute d(phi)/dq (= dRotdq):
	d.SetSize(SOS(),Dim());
	d.SetAll(0);

	double x0 = ploc.X()*2./size.X();
	for (int i=1; i <= NSTor(); i++) //rotation about x-axis
		d(2*(NSPos())+i,1) += GetS0theta(x0,i);

	for (int i=1; i <= NSPos(); i++) //rotation about y-axis
		d((NSPos())+i,2) += -GetS0wx(x0,i);

	for (int i=1; i <= NSPos(); i++) //rotation about z-axis
		d(i,3) += GetS0vx(x0,i);

	//compute d(phi)/dq x v
	for (int i=1; i <= SOS(); i++)
	{
		Vector3D phi(d(i,1),d(i,2),d(i,3));
		Vector3D phi_x_v = rot0*(phi.Cross(vloc));
		d(i,1) = phi_x_v.X();
		d(i,2) = phi_x_v.Y();
		d(i,3) = phi_x_v.Z();
	}
};

void Beam3D::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	GetdRotvdqTLoc(vloc, ploc,d);
	TransformationToNodeDOF(d);
};

void Beam3D::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	//ploc=[-size.X()/2..size.X()/2]
	d.SetSize(SOS(),Dim());
	d.SetAll(0);

	double x0 = ploc.X()*2./size.X();

	//rotations:
	GetdRotvdqTLoc(Vector3D(0.,ploc.Y(),ploc.Z()), ploc, d);

	//deformation of beam axis
	for (int i=1; i <= NSAx()*useAllDOF; i++) //U
	{
		Vector3D v = rot0*Vector3D(GetS0u(x0,i),0.,0.);
		d(2*(NSPos())+NSTor()+i,1) += v.X();
		d(2*(NSPos())+NSTor()+i,2) += v.Y();
		d(2*(NSPos())+NSTor()+i,3) += v.Z();
	}

	for (int i=1; i <= NSPos(); i++) //V
	{
		Vector3D v = rot0*Vector3D(0.,GetS0v(x0,i),0.);
		d(i,1) += v.X();
		d(i,2) += v.Y();
		d(i,3) += v.Z();
	}

	for (int i=1; i <= NSPos(); i++) //W
	{
		Vector3D v = rot0*Vector3D(0.,0.,GetS0w(x0,i));
		d(i+NSPos(),1) += v.X();
		d(i+NSPos(),2) += v.Y();
		d(i+NSPos(),3) += v.Z();
	}

	TransformationToNodeDOF(d);
}

//$ SW 2013-08-28: added
void Beam3D::GetNodedPosdqT(int node, Matrix& dpdqi){
	if (node == 1)
	{
		assert(SOS() == 12);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
		dpdqi(1,1) = 1;
		dpdqi(2,2) = 1;
		dpdqi(3,3) = 1;
	}
	else if (node == 2)
	{
		assert(SOS() == 12);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
		dpdqi(7,1) = 1;
		dpdqi(8,2) = 1;
		dpdqi(9,3) = 1;
	}
	else {
		GetMBS()->UO(UO_LVL_warn) << "WARNING: Beam3D: GetNodePosdqT: wrong local node number detected.";
		assert(0);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
	}
}

void Beam3D::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	d.SetSize(SOS(),Dim());
	d.SetAll(0);

	double x0 = ploc.X()*2./size.X();
	for (int i=1; i <= NSTor(); i++) //rotation about x-axis
		d(2*(NSPos())+i,1) += GetS0theta(x0,i);

	for (int i=1; i <= NSPos(); i++) //rotation about y-axis
		d((NSPos())+i,2) += -GetS0wx(x0,i);

	for (int i=1; i <= NSPos(); i++) //rotation about z-axis
		d(i,3) += GetS0vx(x0,i);

	TransformMatrix(rot0,d);
	TransformationToNodeDOF(d);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis RotorBeamXAxis 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void RotorBeamXAxis::SetRotorBeamXAxis(Vector3D p0i, Matrix3D rot0i, int n1i, int n2i, int materialnumi, 
                   		 Vector2D dim, int axialdeformationI, Vector3D coli, Vector q0i)
{
	Vector3D si = Vector3D(dim.X(),dim.Y()*2,dim.Y()*2); // dim : (length, diameter), si : (lx, ly, lz) // $ MSax 2013-04-16 added factor 2 ==> radius to diameter
	SetBeam3D(p0i, rot0i, n1i, n2i, materialnumi, si, axialdeformationI, coli, q0i);
}


// standard set function with oriented nodes (default set function for script language)
void RotorBeamXAxis::SetRotorBeamXAxis(int n1nr, int n2nr, int matnr)
{
	Node3DR123& n1 = (Node3DR123&)mbs->GetNode(n1nr);
	Node3DR123& n2 = (Node3DR123&)mbs->GetNode(n2nr);
	Matrix3D f1 = n1.GetLocalFrame();
	Matrix3D f2 = n2.GetLocalFrame();

	assert((f1-f2).AbsNorm() < 1e-14);
	
	rot0 = f1;
	p0 = (n1.Pos()+n2.Pos())*0.5;

	Vector2D dim;
	dim.X() = (n1.Pos()-n2.Pos()).Norm(); //CalculateElementLength(n1,n2);
	Material mat = GetMaterial();
	if (mat.IsMaterialOfBeamWithCircularCrossSection())
	{
		dim.Y() = GetMaterial().GetBeamRadius();
	}
	else
	{
		assert (GetMaterial().IsMaterialOfBeamWithCircularCrossSection());   //$ PG 2012-12-18: (as yet) a rectangular cross section (which is specified by the material) is required in ANCFBeam3DTorsion
	}

	SetRotorBeamXAxis(p0, rot0, n1nr, n2nr, matnr, dim, axialdeformation, col);
}

int RotorBeamXAxis::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	Node& node1 = mbs->GetNode(n1);
	Node& node2 = mbs->GetNode(n1);

	if((node1.GetTypeName()).CStrCompare("Node3DR123") != 1 || (node2.GetTypeName()).CStrCompare("Node3DR123") != 1 )
	{
		rv = 2;
		errorstr = mystr("RotorBeamXAxis: Wrong node type detected!\n");
	}

	if(size.X()*size.Y()*size.Z() == 0.)
	{
		rv = 2;
		errorstr = mystr("RotorBeamXAxis: length and radius must be greater than zero!\n");
	}

	return rv;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RotorBeamXAxis::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);
}

int RotorBeamXAxis::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);
	SetRotorBeamXAxis(n1, n2, materialnum);
	return rv;
}

int RotorBeamXAxis::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Body3D::GetAvailableSpecialValues(available_variables);

	available_variables.Add(ReadWriteElementDataVariableType("Internal.acceleration",SOS(),0,0.,mystr("accelerations of the element. range: 1-") + mystr(SOS()) )) ;

	GetAvailableSpecialValuesAuto(available_variables);
	return 0;
}

int RotorBeamXAxis::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Element::ReadSingleElementData(RWdata);

	if(strcmp(RWdata.variable_name.c_str(), "Internal.acceleration") == 0)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= SOS()) //range check
		{
			RWdata.value = XGPP(RWdata.comp1); return 1; 
		}
		else return -2; 
	}

	if (rv == 1) return 1;

	// manual things to read  

	return ReadSingleElementDataAuto(RWdata);
}

Matrix3D RotorBeamXAxis::GetRefRotN1() const
{
	Node3DR123& node_1 = (Node3DR123&)mbs->GetNode(n1);
	if(UseGlobalDOF())
	{
		return node_1.GetLocalFrame();
	}
	return Matrix3D(1.);
}
Matrix3D RotorBeamXAxis::GetRefRotN2() const
{
	return GetRefRotN1();
}

void RotorBeamXAxis::GetRot(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(0.,GetAngleY(x0),GetAngleZ(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	Matrix3D rotTheta;
	GetRotTheta(x, rotTheta);
	rot = rot0*rot*rotTheta;
}

void RotorBeamXAxis::GetRotP(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1

	rot.SetSkew(0.,GetAngleY(x0),GetAngleZ(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	
	Matrix3D rotP;
	rotP.SetSkew(0.,GetAngleYP(x0),GetAngleZP(x0));

	Matrix3D rotTheta;
	GetRotTheta(x, rotTheta);
	Matrix3D rotThetaP;
	GetRotThetaP(x, rotThetaP);
	rot = rot0*(rot*rotThetaP+rotP*rotTheta);
}

void RotorBeamXAxis::GetRotPD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1

	rot.SetSkew(0.,GetAngleYD(x0),GetAngleZD(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	
	Matrix3D rotP;
	rotP.SetSkew(0.,GetAngleYPD(x0),GetAngleZPD(x0));

	Matrix3D rotThetaD;
	GetRotThetaD(x, rotThetaD);
	Matrix3D rotThetaPD;
	GetRotThetaPD(x, rotThetaPD);
	rot = rot0*(rot*rotThetaPD+rotP*rotThetaD);
}

void RotorBeamXAxis::GetRotD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(0.,GetAngleYD(x0),GetAngleZD(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	Matrix3D rotThetaD;
	GetRotThetaD(x, rotThetaD);
	rot = rot0*rot*rotThetaD;
}

void RotorBeamXAxis::GetRotInit(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	rot.SetSkew(0.,GetAngleYInit(x0),GetAngleZInit(x0));
	rot(1,1)+=1; rot(2,2)+=1; rot(3,3)+=1; 
	
	double thetaInit = /*GetAngleXInit(x0)+ $ MSax 2013-07-12 : removed*/(GetAngleXD(-1)+GetAngleXD(1))/2.; // $ MSax 2013-06-14 : added average angle for proper scaling
	Matrix3D rotThetaInit(1.,0.,0.,0.,cos(thetaInit),-sin(thetaInit),0.,sin(thetaInit),cos(thetaInit));
	rot = rot0*rot*rotThetaInit;
}

void RotorBeamXAxis::GetdRotvdqTLoc(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	// $ MSax 2013-04-24 rotations about y- and z- axes are considered as small
	d.SetSize(SOS(),Dim());
	d.SetAll(0);
	double x0 = ploc.X()*2./size.X();

	for (int i=1; i <= NSPos(); i++)
	{
		d(i,1) = -GetS0vx(x0,i)*cos(GetTheta(x0))*vloc.Y()+GetS0vx(x0,i)*sin(GetTheta(x0))*vloc.Z();		
	}
	for (int i=1; i <= NSPos(); i++)
	{
		d(NSPos()+i,1) = -GetS0wx(x0,i)*cos(GetTheta(x0))*vloc.Z()-GetS0wx(x0,i)*sin(GetTheta(x0))*vloc.Y();		
	}
	for (int i=1; i <= NSPos(); i++)
	{
		d(i,2) = GetS0vx(x0,i)*vloc.X();		
	}
	for (int i=1; i <= NSTor(); i++)
	{
		d(2*NSPos()+i,2) = -sin(GetTheta(x0))*GetS0theta(x0,i)*vloc.Y()-cos(GetTheta(x0))*GetS0theta(x0,i)*vloc.Z();		
	}
	for (int i=1; i <= NSPos(); i++)
	{
		d(NSPos()+i,3) = GetS0wx(x0,i)*vloc.X();		
	}
	for (int i=1; i <= NSTor(); i++)
	{
		d(2*NSPos()+i,3) = (sin(GetTheta(x0))*GetS0theta(x0,i)*GetAngleZ(x0)+cos(GetTheta(x0))*GetS0theta(x0,i)*GetAngleY(x0))*vloc.Y()+(cos(GetTheta(x0))*GetS0theta(x0,i)*GetAngleZ(x0)-sin(GetTheta(x0))*GetS0theta(x0,i)*GetAngleY(x0))*vloc.Z();		
	}

	TransformMatrix(rot0,d);
};

void RotorBeamXAxis::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	GetdRotvdqTLoc(vloc, ploc, d);
	TransformationToNodeDOF(d);
};

void RotorBeamXAxis::DrawElement() 
{
	mbs->SetColor(col);

	double lx1 = 0.5*size.X(); double ly1 = 0.5*size.Y()*GetMBS()->GetMagnifyYZ(); double lz1 = 0.5*size.Z()*GetMBS()->GetMagnifyYZ();

	double radius = ly1;
	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	double thickness = GetMBS()->GetDOption(102);
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	}
	if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		colormode = 1;

	double tilex = GetMBS()->GetIOption(137);
	double tiler = Sqr(GetMBS()->GetIOption(138)); //no tiling in transverse direction

	TArray<Vector3D> points((int)(tilex+1)*(int)(tiler+1));
	TArray<double> vals((int)(tilex+1)*(int)(tiler+1));
	double v=0;

	for (int side=1; side <= 3; side++)
	{
		points.SetLen(0); vals.SetLen(0);
		Vector3D px, vx, vr;
		int tilern = (int)tiler;
		int tilexn = (int)tilex;

		switch(side)
		{
		case 1:
		{ //cylinder
			px = Vector3D(-lx1,0.,0.);
			vx = Vector3D(2.*lx1/tilexn,0,0);
			break;
		}
		case 2:
		{ //left side
			px = Vector3D(-lx1,0., 0.);
			break;
		}
		case 3:
		{ //right side
			px = Vector3D(lx1,0., 0.);
			break;
		}
		default: ;
		}

		double del_phi = 2.*MY_PI/tilern;
		for (double ir = 0; ir <= tilern+1e-10; ir++)
		{
			for (double ix = 0; ix <= tilexn+1e-10; ix++)
			{
				Vector3D ploc;
				if (side==1) 
				{
					ploc = px+ix*vx+Vector3D(0.,-radius*cos(del_phi*ir),-radius*sin(del_phi*ir));
				} else 
				{
					ploc = px+Vector3D(0.,-radius/tilexn*ix*cos(del_phi*ir),-radius/tilexn*ix*sin(del_phi*ir));
				}
				Vector3D p(GetPosD(ploc));
				points.Add(p);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(),ploc, true);
				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tilern+1,colormode,0);

		for (double ir = 0; ir <= tilern+1e-10; ir++) // draw lines at last for better visibility
		{
			if (linemode != 0)
			{
				mbs->MyDrawLine(GetPosD(px+Vector3D(0.,-radius*cos(del_phi*ir),-radius*sin(del_phi*ir))),GetPosD(px+Vector3D(0.,-radius*cos(del_phi*(ir+1.)),-radius*sin(del_phi*(ir+1.)))),thickness);
			}
		}
	}
};

void RotorBeamXAxis::GetdPosdx(const Vector3D& ploc, Vector3D& dpdx)
{
	UO() << "Error: RotorBeamXAxis::Not implemented function GetdPosdx!\n";
	assert(0);
};


//$ SW 2013-08-28: added
void RotorBeamXAxis::GetNodedPosdqT(int node, Matrix& dpdqi){
	if (node == 1)
	{
		assert(SOS() == 12);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
		//dpdqi.SetSubmatrix(GetRotMatrix(GetNodeLocPos(node)),1,1);  exact solution
		dpdqi.SetSubmatrix(rot0.GetTp(),1,1);
	}
	else if (node == 2)
	{
		assert(SOS() == 12);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
		//dpdqi.SetSubmatrix(GetRotMatrix(GetNodeLocPos(node)),7,1);
		dpdqi.SetSubmatrix(rot0.GetTp(),7,1);
	}
	else {
		UO() << "Error: RotorBeamXAxis::GetNodedPosdqT: local node number is " 
			<< node << " but should be 1 or 2";
		assert(0);
		dpdqi.SetSize(SOS(),Dim());
		dpdqi.SetAll(0);
	}
}

void RotorBeamXAxis::GetRotTheta(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	double theta = GetAngleX(x0);
	rot.SetAll(0);
	rot(1,1) = 1;
	rot(2,2) = cos(theta);
	rot(2,3) = -sin(theta);
	rot(3,2) = sin(theta);
	rot(3,3) = cos(theta);
}

void RotorBeamXAxis::GetRotThetaP(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	double theta = GetAngleX(x0);
	double thetaP = GetAngleXP(x0);

	rot.SetAll(0);
	rot(2,2) = -sin(theta)*thetaP;
	rot(2,3) = -cos(theta)*thetaP;
	rot(3,2) = cos(theta)*thetaP;
	rot(3,3) = -sin(theta)*thetaP;
}

void RotorBeamXAxis::GetRotThetaD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	double theta = GetAngleXD(x0);
	rot.SetAll(0);
	rot(1,1) = 1;
	rot(2,2) = cos(theta);
	rot(2,3) = -sin(theta);
	rot(3,2) = sin(theta);
	rot(3,3) = cos(theta);
}

void RotorBeamXAxis::GetRotThetaPD(double x, Matrix3D& rot) const
{
	double x0 = x*2./size.X(); //x0 is normalized to -1..+1
	double theta = GetAngleXD(x0);
	double thetaP = GetAngleXPD(x0);

	rot.SetAll(0);
	rot(2,2) = -sin(theta)*thetaP;
	rot(2,3) = -cos(theta)*thetaP;
	rot(3,2) = cos(theta)*thetaP;
	rot(3,3) = -sin(theta)*thetaP;
}

//p_loc is from -size.X()*0.5..size.X()*0.5, etc.
Vector3D RotorBeamXAxis::GetPosD(const Vector3D& p_loc) const
{
	double x0 = p_loc.X()*2./size.X();
	Matrix3D rotD;
	GetRotD(p_loc.X(),rotD);
	Matrix3D rotInit;
	GetRotInit(p_loc.X(),rotInit);
	Vector3D unscaled_position = p0 + rot0*Vector3D(p_loc.X()+GetUD(x0),GetVD(x0),GetWD(x0)) + rotD*Vector3D(0,p_loc.Y(),p_loc.Z());
	Vector3D initial_position_rotated = p0 + rot0*Vector3D(p_loc.X()+GetUInit(x0),GetVInit(x0),GetWInit(x0)) + rotInit*Vector3D(0,p_loc.Y(),p_loc.Z());
	return initial_position_rotated + GetMBS()->GetDOption(105)*(unscaled_position - initial_position_rotated);
};

double RotorBeamXAxis::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD) // $ MSax 2013-06-14 : added because there are some differences compared to Beam3D, see comments
{
	// local_position is local position in [-size.X()/2., +size.X()/2.] x [-size.Y()/2., +size.Y()/2.] x [-size.Z()/2., +size.Z()/2.]
	double beamEA = GetMaterial().BeamEA();
	double beamEIz = GetMaterial().BeamEIz();
	double beamEIy = GetMaterial().BeamEIy();
	double beamGJkx = GetMaterial().BeamGJkx();
	if(!flagD)
	{
		////ploc is from -1 .. +1
		//double x0 = local_position.X();
		//x0 is from -1 .. +1
		double x0 = local_position.X()*2./size.X();

		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:	
			return fvd.GetComponent(GetPos(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			{
				Vector3D pr0 = GetInitRefPos();
				
				Matrix3D rot; // $ MSax 2013-06-14 : added matrix rot
				double theta = (GetAngleX(-1)+GetAngleX(1))/2.;
				rot.SetAll(0);
				rot(1,1) = 1;
				rot(2,2) = cos(theta);
				rot(2,3) = -sin(theta);
				rot(3,2) = sin(theta);
				rot(3,3) = cos(theta);

				return fvd.GetComponent(GetPos(local_position)-pr0-rot0*rot*local_position); // $ MSax 2013-06-14 : added 
			}
		case FieldVariableDescriptor::FVT_beam_torsion: return (GetAngleX(1)-GetAngleX(-1))/size.X();
		case FieldVariableDescriptor::FVT_beam_force_axial: return beamEA * GetUx(x0);
		case FieldVariableDescriptor::FVT_beam_force_transversal:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
				{
					return -beamEIz * GetVxxx(x0); //Qy=-dMz/dx; see Mang, page 196, p.176, p.178
				}
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
				{
					return -beamEIy * GetWxxx(x0); //Qz= dMy/dx; see Mang, page 196, p.176, p.178

				}
			}
		case FieldVariableDescriptor::FVT_beam_moment_torsional: return beamGJkx * GetThetax(x0);	//Mx = theta_x * GJ *kx
		case FieldVariableDescriptor::FVT_beam_moment_bending:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
					return beamEIy * GetWxx(x0);			//My = -EIy * w''; see Mang, page 196, p.176
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
					return beamEIz * GetVxx(x0);			//Mz =  EIz * v''; see Mang, page 196, p.178
			}
		case FieldVariableDescriptor::FVT_acceleration:
			{
				return fvd.GetComponent(GetAcceleration(local_position));
			}
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2013-01-29 added
			return fvd.GetComponent(GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2013-01-29 added
			return fvd.GetComponent(GetRotMatrix(local_position).GetTp()*GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2013-01-29 added
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrix(local_position), phi);
				return fvd.GetComponent(phi);
			}
		}
	} 
	else 
	{
		////ploc is from -1 .. +1
		//double x0 = local_position.X();
		//x0 is from -1 .. +1
		double x0 = local_position.X()*2./size.X();

		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:		
			return fvd.GetComponent(GetPosD(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVelD(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			{
				Vector3D pr0 = GetInitRefPos();

				Matrix3D rot; // $ MSax 2013-06-14 : added matrix rot
				double theta = (GetAngleXD(-1)+GetAngleXD(1))/2.;
				rot.SetAll(0);
				rot(1,1) = 1;
				rot(2,2) = cos(theta);
				rot(2,3) = -sin(theta);
				rot(3,2) = sin(theta);
				rot(3,3) = cos(theta);

				return fvd.GetComponent(GetPosD(local_position)-pr0-rot0*rot*local_position); // $ MSax 2013-06-14 : added 
			}
		case FieldVariableDescriptor::FVT_beam_force_axial: return beamEA * GetUxD(x0);
		case FieldVariableDescriptor::FVT_beam_torsion: return (GetAngleXD(1)-GetAngleXD(-1))/size.X();
		case FieldVariableDescriptor::FVT_beam_force_transversal:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
				{
					return -beamEIz * GetVxxxD(x0); //Qy=-dMz/dx; see Mang, page 196, p.176, p.178
				}
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
				{
					return -beamEIy * GetWxxxD(x0); //Qz= dMy/dx; see Mang, page 196, p.176, p.178

				}
			}
		case FieldVariableDescriptor::FVT_beam_moment_torsional: return beamGJkx * GetThetaxD(x0);	//Mx = theta_x * GJ *kx
		case FieldVariableDescriptor::FVT_beam_moment_bending:
			{
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_y)
					return beamEIy * GetWxxD(x0);			//My = -EIy * w''; see Mang, page 196, p.176
				if(fvd.ComponentIndex1() == FieldVariableDescriptor::FVCI_z)
					return beamEIz * GetVxxD(x0);			//Mz =  EIz * v''; see Mang, page 196, p.178
			}
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2013-01-29 added
			return fvd.GetComponent(GetAngularVelD(local_position));//$ SW 2013-10-25: changed to GetAngularVelD
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2013-01-29 added
			return fvd.GetComponent(GetRotMatrixD(local_position).GetTp()*GetAngularVelD(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2013-01-29 added
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrixD(local_position), phi);
				return fvd.GetComponent(phi);
			}
		}
	}

	// if not implemented here, call parent class
	//FVT_position, FVT_velocity: implemented in Element
	return Body3D::GetFieldVariableValue(fvd, local_position, flagD);
}