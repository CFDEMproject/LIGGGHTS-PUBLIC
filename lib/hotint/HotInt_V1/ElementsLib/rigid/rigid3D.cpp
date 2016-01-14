//#**************************************************************
//#
//# filename:             rigid3D.cpp
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
 
#include "element.h"
#include "material.h"
#include "node.h"
#include "solversettings_auto.h"
#include "options_class_auto.h"
#include "rigid3d.h"
#include "graphicsconstants.h"
#include "rendercontext.h"
#include "elementdataaccess.h"
#include "geomelements.h"


double Body3D::GetAngleAroundAxis(const Vector3D& p_loc, const Vector3D& axis) const
{
	Vector3D laxis = GetRotMatrix(p_loc).GetTp()*axis;

	return GetAngleAroundLocalAxis(p_loc, laxis);
}

double Body3D::GetAngleAroundLocalAxis(const Vector3D& p_loc, const Vector3D& axis) const
{
	//local/global axis ....
	Vector3D t,n, tglob, nglob;
	Vector3D rot = axis;

	rot.SetNormalBasis(t,n);
	tglob = GetRotMatrix(p_loc)*t;
	nglob = GetRotMatrix(p_loc)*n;

	//Vector2D t1(1,0);
	Vector2D t2(tglob*t, tglob*n);
	double len = t2.Norm();
	if (len == 0) 
	{
		Vector2D n2(nglob*t, nglob*n);
		double len = n2.Norm();
		if (len == 0) return 0;

		n2 *= 1./len;

		double phi = atan2(n2.Y(), n2.X());
		if (phi > MY_PI) phi -= 2.*MY_PI;
		if (phi <-MY_PI) phi += 2.*MY_PI;

		return phi;
	}

	t2 *= 1./len;

	double phi = atan2(t2.Y(), t2.X());
	if (phi > MY_PI) phi -= 2.*MY_PI;
	if (phi <-MY_PI) phi += 2.*MY_PI;

	return phi;
}

Vector3D Body3D::GetRefPosD() const 
{
	if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.)
	{
		Vector3D p = Vector3D(XGD(1),XGD(2),XGD(3));
		Vector3D p0 = GetRefPosInit();
		return (p - p0)*GetMBS()->GetDOption(105) + p0;
	}
	else
	{
		return Vector3D(XGD(1),XGD(2),XGD(3));
	}
}

//Vector3D Body3D::GetPos_dc(const Vector3D& ploc, TComputeDrawInitFlag flag) const
//{
//	if (flag&TCD_draw && GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
//	{
//		Matrix3D A = GetRotMatrix_dc(TCD_draw);
//		Matrix3D A0= GetRotMatrix_dc(TCD_initial_values);
//		double fact = GetMBS()->GetDOption(105);
//
//		return (fact*A-(fact-1.)*A0)*p_loc+GetRefPos_dc(TCD_draw);
//	}
//	else
//	{
//	  return GetRotMatrix_dc(flag)*p_loc+GetRefPos_dc();
//	}
//}

Vector3D Body3D::GetPos(const Vector3D& p_loc) const
{
	return GetRotMatrix()*p_loc+GetRefPos();
};

Vector3D Body3D::GetDisplacement(const Vector3D& p_loc) const
{
	return GetRotMatrix()*p_loc + GetRefPos() - (GetRotMatrixInit()*p_loc + GetRefPosInit());
};

Vector3D Body3D::GetVel(const Vector3D& p_loc) const
{
	return GetRotMatrixP()*p_loc+GetRefVel();
};

//functions for drawing:
Vector3D Body3D::GetPosD(const Vector3D& p_loc) const
{
	if (GetMBS()->GetIOption(151) && IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
	{
		Matrix3D A = GetRotMatrixD();
		Matrix3D A0= GetRotMatrixInit();
		double fact = GetMBS()->GetDOption(105);

		return (fact*A-(fact-1.)*A0)*p_loc+GetRefPosD();
	}
	else
	{
		return GetRotMatrixD()*p_loc+GetRefPosD();
	}
};

Vector3D Body3D::GetDisplacementD(const Vector3D& p_loc) const
{
	return GetRotMatrixD()*p_loc + GetRefPosD() - (GetRotMatrixInit()*p_loc + GetRefPosInit());
};

Vector3D Body3D::GetVelD(const Vector3D& p_loc) const
{
	return GetRotMatrixPD()*p_loc+GetRefVelD();
};

void Body3D::DrawElement() 
{
	Element::DrawElement();
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++          MASS3D            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$ DR 2012-07: CheckConsistency added
int Mass3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!mass)
	{
		rv = 1;
		errorstr = mystr("No mass is defined for the Mass3D!\n");
	}

	return rv;
}

void Mass3D::ElementDefaultConstructorInitialization() //$ DR 2012-07: ElementDefaultConstructorInitialization added
{
	drawres = 6; 
	x_init = Vector(SS()); //zero initialized
	elementname = GetElementSpec();
	size.X() = 0.1;						// DR just some default value
}

void Mass3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);
}

int Mass3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);
	return rv;
}


void Mass3D::DrawElement() 
{
	Body3D::DrawElement();

	Vector3D pos = GetRefPosD();
	Vector3D refpos = GetRefPosInit(); //$JG2012-02-21: old: Vector3D(GetXInit()(1),GetXInit()(2),GetXInit()(3));

	DrawElementFunc(pos, refpos);
}

void Mass3D::DrawElementFunc(const Vector3D& pos, const Vector3D& refpos) //this function does the drawing
{
	//double factor = GetMBS()->GetDOption(105); //magnification of displacements
	//if (!GetMBS()->GetIOption(151)) factor = 1;


	double drawrad = size.X();
	double factor = mbs->GetOptions()->PostProcOptions()->BodiesDeformationScaleFactor();
	if (!GetMBS()->GetIOption(151)) factor = 1;

	if (type & TParticle)
	{
		drawrad *= mbs->GetOptions()->PostProcOptions()->BodiesParticlesDrawSizeFactor();
		factor = mbs->GetOptions()->PostProcOptions()->BodiesParticlesDisplacementScaleFactor();
	}

	Vector3D p = factor * (pos - refpos) + refpos;


	//double drawrad = size.X();
	//if (type & TParticle) drawrad *= mbs->GetDOption(155);    //GetDOption(155) corresponds to PostProcOptions.Bodies.Particles.draw_size_factor

	if (mbs->GetIOption(134)) //draw sphere
	{

		if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL)
		{
			double v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), Vector3D(0.), true);
			mbs->DrawColorSphere(v, p, drawrad, drawres);
		}
		else
		{
			mbs->SetColor(col);
			mbs->DrawSphere(p, drawrad, drawres);
		}
	}

	if (mbs->GetIOption(133)) //only draw a point for outline drawing
	{
		mbs->SetColor(colgrey1);
		mbs->GetRC()->glBeginPoints();
		mbs->GetRC()->ChoosePointSize((float)mbs->GetIOption(117));
		mbs->GetRC()->glVertex(p.GetVecPtr());
		mbs->GetRC()->glEnd();
	}
}


double Mass3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{	
	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:		
			return fvd.GetComponent(GetPosD(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVelD(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPosD(local_position) - Vector3D(x_init(1), x_init(2), x_init(3)));
		case FieldVariableDescriptor::FVT_acceleration: //$ MS 2013-9-4: Added
			return fvd.GetComponent(GetAcceleration(local_position));
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:	
			return fvd.GetComponent(GetPos(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetPos(local_position) - Vector3D(x_init(1), x_init(2), x_init(3)));
		case FieldVariableDescriptor::FVT_acceleration:  //$ MS 2013-9-4: Added
			return fvd.GetComponent(GetAcceleration(local_position));

		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}


void Mass3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	
	// add some specific ones of this class
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration, FieldVariableDescriptor::FVCI_z); //$ MS 2013-9-4: Added //$ SW 2013-10-21:added FieldVariableDescriptor::FVCI_z
	
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++         NODALMASS3D            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NodalMass3D::LinkToElements()
{
	LTGreset();

	// Node positions
	for (int j = 1; j <= NNodes(); j++)
	{
		const Node& node = GetNode(j);
		//Position:
		for (int i=1; i <= node.SOS(); i++)
		{
			AddLTG(node.Get(i));
		}
	}

	// Node velocities
	for (int j = 1; j <= NNodes(); j++)
	{
		const Node& node = GetNode(j);
		//velocity:
		for (int i=1; i <= node.SOS(); i++)
		{
			AddLTG(node.Get(i+node.SOS()));
		}
	}
}

void NodalMass3D::Initialize() 
{
	Mass3D::Initialize();

	x_init = GetNode(1).X_Init();
	if (x_init.Length() != 6) 
	{
		x_init = Vector(6);
	//	GetMBS()->UO().InstantMessageText("Error: NodalMass3D: Node does not have initialization!");
	}
};


void NodalMass3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);

	//int dind = edc.Find("Density");
	//if (dind) edc.Delete(dind);
	ElementData* edf = edc.TreeFind("Physics.density");
	if (edf) edf->SetLocked(1);

	//ElementData ed;

	//ed.SetDouble(mass, "Mass"); edc.Add(ed);
	//ed.SetDouble(size.X(), "Radius"); edc.Add(ed);

	//ed.SetInt(nodenum, "Node_number"); edc.Add(ed);

	//ed.SetInt(drawres, "Draw_resolution", 1, 1000); ed.SetToolTipText("Number of segments used to approximate circumference of sphere"); edc.Add(ed);
}

int NodalMass3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);

	GetElemDataDouble(GetMBS(), edc, "Mass", mass, 1);

	double r = size.X();
	//GetElemDataDouble(GetMBS(), edc, "Radius", size.X(), 1);
	//size = r;
	//if (r != 0) 
	//{
	//	double rho = mass/ (4./3. * MY_PI * Cub(r));
	//	GetMaterial().Density() = rho;
	//}

	//GetElemDataInt(GetMBS(), edc, "Draw_resolution", drawres, 1);
	//GetElemDataInt(GetMBS(), edc, "Node_number", nodenum, 1);


	return rv;
}

void NodalMass3D::DrawElement() 
{
	Body3D::DrawElement();

	Vector3D pos = GetRefPosD();
	Vector3D refpos = GetMBS()->GetNode(NodeNum(1)).RefConfPos();

	DrawElementFunc(pos, refpos);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++         NODALDISKMASS3D            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void NodalDiskMass3D::LinkToElements()
{
	LTGreset();

	// Node positions
	for (int j = 1; j <= NNodes(); j++)
	{
		const Node& node = GetNode(j);
		//Position:
		for (int i=1; i <= node.SOS(); i++)
		{
			AddLTG(node.Get(i));
		}
	}

	// Node velocities
	for (int j = 1; j <= NNodes(); j++)
	{
		const Node& node = GetNode(j);
		//velocity:
		for (int i=1; i <= node.SOS(); i++)
		{
			AddLTG(node.Get(i+node.SOS()));
		}
	}
}

void NodalDiskMass3D::Initialize() 
{
	Body3D::Initialize();

	x_init = GetNode(1).X_Init();
	if (x_init.Length() != 12) 
	{
		x_init = Vector(12);
	}
};

//$ SW 2013-11-20: added
void NodalDiskMass3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> &variables)
{
	// add all the field variables of the parent class
	NodalMass3D::GetAvailableFieldVariables(variables);
	
	// add some specific ones of this class
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity, FieldVariableDescriptor::FVCI_z);
}

//$ SW 2013-11-20: added
double NodalDiskMass3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{	
	if (!flagD)
	{
		if (fvd.VariableType() == FieldVariableDescriptor::FVT_angular_velocity)
		{
			return fvd.GetComponent(GetAngularVel(local_position));
		}
	}
	else
	{
		if (fvd.VariableType() == FieldVariableDescriptor::FVT_angular_velocity)
		{
			return fvd.GetComponent(GetAngularVelD(local_position));
		}
	}
	return NodalMass3D::GetFieldVariableValue(fvd,local_position,flagD);
}

void NodalDiskMass3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Mass3D::GetElementData(edc);
}

int NodalDiskMass3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Mass3D::SetElementData(edc);
	x_init = GetNode(1).X_Init();
	if (x_init.Length() != 12) 
	{
		x_init = Vector(12);
	}
	return rv;
}

int NodalDiskMass3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = NodalMass3D::CheckConsistency(errorstr);
	if (rv) return rv;

	Node& n = mbs->GetNode(nodenum);

	if((n.GetTypeName()).CStrCompare("Node3DR123") != 1)
	{
		rv = 2;
		errorstr = mystr("NodalDiskMass3D: Wrong node type detected. You must use 'Node3DR123'!\n");
	}

	return rv;
}

void NodalDiskMass3D::DrawElement() 
{
	Body3D::DrawElement();

	double drawthickness = size.X();
	double drawrad = size.Y();

	if (mbs->GetIOption(134)) //draw cylinder
	{
		mbs->SetColor(col);
		mbs->DrawZyl(GetPosD(Vector3D(drawthickness*0.5,0,0)),GetPosD(Vector3D(-1*drawthickness*0.5,0,0)),drawrad,drawres);
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++         RIGID3D            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Rigid3D::Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
								 const Vector3D& si, const Vector3D& coli): Body3D(mbsi)
{
	SetRigid3D(x0, phi0, rhoi, Vi, Ip, si, coli);
};

Rigid3D::Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
								 const Vector3D& si, const Vector3D& coli):Body3D(mbsi)
{
	SetRigid3D(x0, phi0, rhoi, Vi, Ip, si, coli);
};

Rigid3D::Rigid3D(MBS* mbsi, const Vector& x0, const Vector& phi0, double rhoi,
								 const Vector3D& si, const Vector3D& coli):
Body3D(mbsi)
{
	SetRigid3D(x0, phi0, rhoi, si, coli);
}

//$ DR 2012-07: ElementDefaultConstructorInitialization added
void Rigid3D::ElementDefaultConstructorInitialization()
{
	
	mass = 1;
	size = Vector3D(1);
	volume = 1;

	Iphi = Matrix3D(1);
	Iphi(1,1)= 1./12.*mass*(Sqr(size.Y())+Sqr(size.Z()));
	Iphi(2,2)= 1./12.*mass*(Sqr(size.X())+Sqr(size.Z()));
	Iphi(3,3)= 1./12.*mass*(Sqr(size.X())+Sqr(size.Y()));

	x_init.SetLen(15);
	x_init.SetAll(0);

	elementname = GetElementSpec();
}

//$ DR 2012-07: CheckConsistency added
int Rigid3D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Element::CheckConsistency(errorstr);
	if (rv) return rv;

	if(!mass)
	{
		rv = 1;
		errorstr = mystr("No mass is defined for the Rigid3D!\n");
	}

	// some check for Iphi?
	return rv;
}

void Rigid3D::SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Vector3D& Ip,
								 const Vector3D& si, const Vector3D& coli)
{
	Iphi = Matrix3D(Ip.X(),Ip.Y(),Ip.Z());

	mass=rhoi*Vi;
	size = si;

	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));

	ComputeInitialConditions(xp, vp, phi, phip, x_init);

	//rho = rhoi;
	col = coli;

	Material mat(GetMBS());
	mat.SetMaterialRigid(rhoi);
	int material_num = GetMBS()->AddMaterial(mat);
	SetMaterialNum(material_num);

	//elementname = GetElementSpec(); //$ DR 2012-07: moved to ElementDefaultConstructorInitialization
}


void Rigid3D::SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi, double Vi, const Matrix3D& Ip,
								 const Vector3D& si, const Vector3D& coli)
{
	Iphi = Ip;

	mass=rhoi*Vi;
	size = si;

	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));

	ComputeInitialConditions(xp, vp, phi, phip, x_init);

	//rho = rhoi;
	col = coli;
	Material mat(GetMBS());
	mat.SetMaterialRigid(rhoi);
	int material_num = GetMBS()->AddMaterial(mat);
	SetMaterialNum(material_num);

	//elementname = GetElementSpec(); //$ DR 2012-07: moved to ElementDefaultConstructorInitialization
}

void Rigid3D::SetRigid3D(const Vector& x0, const Vector& phi0, double rhoi,
								 const Vector3D& si, const Vector3D& coli)
{
	size = si;
	volume = si.X()*si.Y()*si.Z();
	mass = volume*rhoi;

	Iphi.SetAll(0);
	Iphi(1,1)= 1./12.*mass*(Sqr(si.Y())+Sqr(si.Z()));
	Iphi(2,2)= 1./12.*mass*(Sqr(si.X())+Sqr(si.Z()));
	Iphi(3,3)= 1./12.*mass*(Sqr(si.X())+Sqr(si.Y()));

	Vector3D xp(x0(1),x0(2),x0(3));
	Vector3D vp(x0(4),x0(5),x0(6));
	Vector3D phi(phi0(1),phi0(2),phi0(3));
	Vector3D phip(phi0(4),phi0(5),phi0(6));

	ComputeInitialConditions(xp, vp, phi, phip, x_init);

	//rho = rhoi;
	col = coli;
	Material mat(GetMBS());
	mat.SetMaterialRigid(rhoi);
	int material_num = GetMBS()->AddMaterial(mat);
	SetMaterialNum(material_num);

	//elementname = GetElementSpec(); //$ DR 2012-07: moved to ElementDefaultConstructorInitialization
};


void Rigid3D::ComputeInitialConditions(const Vector3D& xp, const Vector3D& vp, 
																			 const Vector3D& phi, const Vector3D& phip, Vector& xinit)
{
	xinit.SetLen(SS()); //15+GGL
	xinit(1) = xp.X(); xinit(2) = xp.Y(); xinit(3) = xp.Z();
	double b0, b1, b2, b3;
	RotMatToQuaternions(ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z()),b0, b1, b2, b3);
	//RotMatToQuaternions(ComputeRotMatrixWithKardanAngles(phi.X(),phi.Y(),phi.Z()),b0, b1, b2, b3);

	xinit(4) = b0; xinit(5) = b1;  xinit(6) = b2;    xinit(7) = b3;

	//xinit(8) = vp.X(); xinit(9) = vp.Y(); xinit(10) = vp.Z();
	xinit(NRotParam()+4) = vp.X(); xinit(NRotParam()+5) = vp.Y(); xinit(NRotParam()+6) = vp.Z();

	//set beta_p = 0, for case of singularity (should not happen!):
	//xinit(11) = 0;
	//xinit(12) = 0;
	//xinit(13) = 0;
	//xinit(14) = 0;
	for(int i=NRotParam()+7; i<=2*NRotParam()+6;i++)
	{
		xinit(i) = 0.;
	}

	//initial values for time-wise derivative quaternions beta_t:
	double beta0 = 2*xinit(4); double beta1 = 2*xinit(5); double beta2 = 2*xinit(6); double beta3 = 2*xinit(7); 
	//the betas are already multiplied with 2, compared to Shabana
	Matrix G(4,4);
	G(1,1) = -beta1; G(1,2) =  beta0; G(1,3) = -beta3; G(1,4) =  beta2;
	G(2,1) = -beta2; G(2,2) =  beta3; G(2,3) =  beta0; G(2,4) = -beta1;
	G(3,1) = -beta3; G(3,2) = -beta2; G(3,3) =  beta1; G(3,4) =  beta0;
	G(4,1) =  beta0; G(4,2) =  beta1; G(4,3) =  beta2; G(4,4) =  beta3; //this is the time-derivative of the quaternion condition | |^2 = 1

	Vector f(phip.X(),phip.Y(),phip.Z(),0.); //omega and zero
	Vector q(4); 
	int rv = G.Solve(f, q);
	if (!rv) {GetMBS()->UO() << "Rigid3D:Initialization: could not determine initial Euler parameter velocities due to singularity!!!\n";}
	else
	{
		xinit(11) = q(1);
		xinit(12) = q(2);
		xinit(13) = q(3);
		xinit(14) = q(4); 
	}

	xinit(15) = 0; //for constraint
	//if (GGLStabilisation()) xinit(16) = 0; //for constraint
}

Element* Rigid3D::GetCopy()
{
	Element* ec = new Rigid3D(*this);
	return ec;
}

void Rigid3D::CopyFrom(const Element& e)
{
	Body3D::CopyFrom(e);
	const Rigid3D& ce = (const Rigid3D&)e;
	Iphi = ce.Iphi;
	volume = ce.volume;
}

void Rigid3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);

	ElementData ed;

	//$ DR 2013-01-30 removed this option, see change log 379
	//ed.SetInt(1,"inertia_mode",1,3); ed.SetToolTipText("1.. directly entered mass (m) and inertia values (I) are used. 2..m and I are computed from density and body dimensions. 3..m and I are computed from density and GeomElement.");  edc.TreeAdd("Physics",ed);

	Vector3D phi;
	QuaternionsToKardanAngles(x_init(4),x_init(5),x_init(6),x_init(7), phi);
	
	ed.SetVector3D(phi.X(), phi.Y(), phi.Z(), "initial_rotation"); 
	ed.SetToolTipText(mystr("3 consecutive rotations (global rotation axes): [rot3_X, rot2_Y, rot1_Z] in rad")); edc.TreeAdd("Initialization",ed);


	Vector3D phip;
	QuaternionsPToAngularVelocities(x_init(4),x_init(5),x_init(6),x_init(7), 
			x_init(11),x_init(12),x_init(13),x_init(14),phip);


	ed.SetVector3D(phip.X(), phip.Y(), phip.Z(), "initial_angular_velocity"); 
	ed.SetToolTipText(mystr("Angular velocity vector in global coordinates: [ang_X, ang_Y, ang_Z] in rad/s")); edc.TreeAdd("Initialization",ed);

	////Vector3D xp(x_init(1),x_init(2),x_init(3));	//$ DR 2012-08-22: not used at all in GetElementData

	//QuaternionsToEulerAngles(x_init(4),x_init(5),x_init(6),x_init(7), phi);
	//	
	//QuaternionsPToAngularVelocities(x_init(4),x_init(5),x_init(6),x_init(7), 
	//	x_init(11),x_init(12),x_init(13),x_init(14),phip);
	//}
	//else
	//{		
	//	// Rotation matrix --> Quaternions --> Euler Angles
	//	//Matrix3D A = GetRotMatrix();
	//	Matrix3D A;
	//	if(ltg.Length())	{	A = GetRotMatrix();	}		
	//	else {	A = Matrix3D(1.);}	//$ DR 2012-12-13: necessary if GetElementData is called before assemble, ltg is not initialized at this moment

	//	double b0, b1, b2, b3;
	//	RotMatToQuaternions(A, b0, b1,b2, b3);
	//	QuaternionsToEulerAngles(b0, b1,b2, b3, phi);
	//	Matrix3D G = GetG(Vector(x_init(4), x_init(5), x_init(6)));
	// 	phip = G*Vector3D(x_init(10), x_init(11), x_init(12)); //omega = G beta_p
	//	
	//}	
	//Vector3D vp(x_init(NRotParam()+4),x_init(NRotParam()+5),x_init(NRotParam()+6));

	//if (!GetMBS()->IsLoadSaveMode() && GetMBS()->GetIOption(140) == 1) //RotXYZ
	//{
	//	QuaternionsToRot1Rot2Rot3Angles(x_init(4),x_init(5),x_init(6),x_init(7), phi);
	//	phi *= frad;
	//	ed.SetVector3D(phi.X(), phi.Y(), phi.Z(), "initial_rotation_XYZ"); 
	//	ed.SetToolTipText((mystr("3 consecutive rotations (global rotation axes): [rot3_X, rot2_Y, rot1_Z] ")+GetRotUnitStr(GetMBS()->GetIOption(120))).c_str()); edc.Add(ed);
	//}
	//else if (!GetMBS()->IsLoadSaveMode() && GetMBS()->GetIOption(140) == 2) //Euler parameters
	//{
	//	Vector v(4);
	//	v(1) = x_init(4); v(2) = x_init(5); v(3) = x_init(6); v(4) = x_init(7); 
	//	ed.SetVector(v.GetVecPtr(), 4, "initial_euler_parameters"); 
	//	ed.SetToolTipText("Euler parameters beta0, beta1, beta2, beta3"); edc.Add(ed);
	//}
	//else //Euler rotations
	//{
	//	phi *= frad;
	//	ed.SetVector3D(phi.X(), phi.Y(), phi.Z(), "initial_euler_angles"); 
	//	ed.SetToolTipText((mystr("3 consecutive rotations (global rotation axes): [rot3_Z, rot2_X, rot1_Z] ")+GetRotUnitStr(GetMBS()->GetIOption(120))).c_str()); edc.Add(ed);
	//}

	//phip *= frad;
	//ed.SetVector3D(phip.X(), phip.Y(), phip.Z(), "initial_angular_velocity"); 
	//ed.SetToolTipText((mystr("Angular velocity vector in global coordinates: [ang_X, ang_Y, ang_Z] ")+GetRotUnitStr(GetMBS()->GetIOption(120))).c_str()); edc.Add(ed);
}

int Rigid3D::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Body3D::SetElementData(edc);

	Vector3D xp, vp, phi, phip;

	xp=Vector3D(x_init(1),x_init(2),x_init(3)); //$ DR 2012-07: is now in Get/SetElementDataAuto]
	vp=Vector3D(x_init(NRotParam()+4),x_init(NRotParam()+5),x_init(NRotParam()+6)); //$ DR 2012-07: is now in Get/SetElementDataAuto]

	GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_rotation", phi, 0);

	GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_angular_velocity", phip, 0);

	//$!DR 2013-01-28: phi are now kardan angles!
	Matrix3D rot = ComputeRotMatrixWithKardanAngles(phi.X(),phi.Y(),phi.Z());
	RotMatToEulerAngles(rot, phi);

	ComputeInitialConditions(xp, vp, phi, phip, x_init);

	//SetElementDataMassAndInertia(edc); //$ DR 2013-01-30 removed this option, see change log 379

	return rv;
}

//$ DR 2012-12-13 sub function of SetElementData, used for Rigid3D and Rigid3DKardan
int Rigid3D::SetElementDataMassAndInertia(ElementDataContainer& edc)
{
	int massmode = 1; // "Use_inertia_values" is default!

	GetElemDataInt(GetMBS(), edc, "Physics.inertia_mode", massmode,1);

	if (massmode == 2)
	{
		//compute inertia from size:
		mass = size.X()*size.Y()*size.Z()*Rho();

		Iphi.SetAll(0);
		Iphi(1,1)= 1./12.*mass*(Sqr(size.Y())+Sqr(size.Z()));
		Iphi(2,2)= 1./12.*mass*(Sqr(size.X())+Sqr(size.Z()));
		Iphi(3,3)= 1./12.*mass*(Sqr(size.X())+Sqr(size.Y()));
	}
	else if (massmode == 3)
	{
		if (NGeomElements() == 0) 
		{
			UO() << "Warning: No GeomElements, mass and inertia can not be computed from GeomElement!!!\n";
			massmode = 1;
		}
		else if (NGeomElements() > 1) 
		{
			UO() << "Warning: More than one GeomElement, mass and inertia is only computed from first GeomElement!!!\n";
		}
		if (NGeomElements() >= 1)
		{
			mass = Rho() * GetGeomElement(1)->ComputeVolume();
			Iphi = GetGeomElement(1)->ComputeMassMomentOfInertia(Rho());
			Vector3D center_of_mass = GetGeomElement(1)->ComputeCenterOfMass();

			if (center_of_mass.Norm() > 1e-14)
			{
				GetGeomElement(1)->Translate(-1*center_of_mass);
				UO() << "Warning: The center of mass of the GeomElement was\nCOM=" << center_of_mass << "\nand has been corrected to (0,0,0). The coordinates have been translated by t=" << -1.*center_of_mass << "!!!\n";
			}
		}
	}

	return 1;
}

void Rigid3D::Initialize() 
{
	Body3D::Initialize();
};

void Rigid3D::EvalG(Vector& f, double t) 
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);

	if (MaxIndex() == 3)
	{	//position level:
		f(1) = Sqr(beta(1))+Sqr(beta(2))+Sqr(beta(3))+Sqr(beta(4))-1.;             //equ1: sum(beta_i,i=0..3)-1=0 ...algebraic equation
	} 
	else
	{	//velocity level:
		ConstVector<4> betap(NRotParam());
		GetBetaP(betap);
		f(1) = (2*beta(1)*betap(1)+2*beta(2)*betap(2)+2*beta(3)*betap(3)+2*beta(4)*betap(4));

		//if (GGLStabilisation())
		//{
		//	f(2) = Sqr(beta(1))+Sqr(beta(2))+Sqr(beta(3))+Sqr(beta(4))-1.;             //equ1: sum(beta_i,i=0..3)-1=0 ...algebraic equation
		//}
	}
}; 

void Rigid3D::EvalM(Matrix& m, double t) 
{
	//compare Shabana 1998, p. 158, eq. 3.156/157 and p.159, eq. 3.164
	//center of mass must be at origin(0,0,0) !!!!!!!!!!!!!!!!!!!!

	m(1,1) = mass; m(2,2) = mass; m(3,3) = mass; 

	Matrix3D mtheta = (GetGbarT()*Iphi)*GetGbar();

	for (int i=0; i<NRotParam(); i++)
		for (int j=0; j<NRotParam(); j++)
			m(4+i,4+j) = mtheta.Get0(i,j);
};

void Rigid3D::EvalF2(Vector& f, double t) 
{
	Body3D::EvalF2(f,t);

	//compare Shabana 1998, p. 158, eq. 3.156/157 and p.159, eq. 3.164
	//quadratic velocity vector (Shabana Computational Dynamics 1994, p. 396 or p. 414 in newer version)
	//alternative representation: -Gbar.GetTp()*omegabar.Cross(Iphi*omegabar) = -2 Gbarp.GetTp()*Iphi*omegabar ... faster???

	//center of mass must be at origin(0,0,0) !!!!!!!!!!!!!!!!!!!!


	ConstVector<4> betap(NRotParam());
	GetBetaP(betap);


	if (1)
	{
		Matrix3D Gbar = GetGbar();
		Vector3D omegabar = Gbar*betap; //--->Mult is faster

		Vector3D temp = (omegabar.Cross(Iphi*omegabar));
		Mult(Gbar.GetTp(),temp,betap);

	}
	else
	{
		//alternative:
		Matrix3D Gbar = GetGbar();
		Vector3D omegabar = Gbar*betap; //--->Mult is faster, alternative Gbarp*beta???
		Matrix3D GbarpTp = GetGbarpT();
		omegabar *= 2;
		Vector3D temp = Iphi*omegabar;

		Mult(GbarpTp,temp,betap);  //betap = GbarpTp * temp;
	}

	for(int i = 1;i<=NRotParam();i++)
	{
		f(i+3) -= betap(i); //betap is now result of the above evaluation!
    }

	//f(4) -= betap(1); //betap is now result of the above evaluation!
	//f(5) -= betap(2);
	//f(6) -= betap(3);
	//f(7) -= betap(4);



	//Add C_q^T terms
	AddEPCqTterms(f);

}; 

Vector3D Rigid3D::GetAngularVel(const Vector3D& p_loc) const
{
	//local position is ignored in rigid body
	Matrix3D G = GetG();
	ConstVector<4> betap(NRotParam());

	GetBetaP(betap);
	Vector3D omega = G*betap;

	return omega;
}

Vector3D Rigid3D::GetAngularVelLocal() const
{
	//local position is ignored in rigid body
	Matrix3D Gbar = GetGbar();
	
	//double mem_betap[4];
	//Vector betap;
	//betap.LinkWith(&mem_betap[0], 4);
	
	////static Vector betap;
	////betap.SetLen(4);

	ConstVector<4> betap(NRotParam());
	GetBetaP(betap);
	Vector3D omega = Gbar*betap;

	return omega;
}

void Rigid3D::GetBetaP(Vector& betap) const
{
	for(int i=1; i<=NRotParam();i++)
	{
		betap(i) = XGP(i+3);	
	}	


//	betap(1) = XGP(4); betap(2) = XGP(5); betap(3) = XGP(6); betap(4) = XGP(7);
//	if (GGLstab && 0)
//	{
//		//add stabilisation terms: d(g)/dq * eta, eta is additional Lagrange Multiplier
//		double eta = XG(16);
//
//		betap0 += 2*XG(4)*eta; //negative sign because dCqTLambda is put on rhs
//		betap1 += 2*XG(5)*eta;
//		betap2 += 2*XG(6)*eta;
//		betap3 += 2*XG(7)*eta;
//	}
}

void Rigid3D::GetBeta(Vector& beta) const
{
	for(int i=1; i<=NRotParam();i++)
	{
		beta(i) = XG(i+3);	
	}	

	//beta(1) = XG(4); beta(2) = XG(5); beta(3) = XG(6); beta(4) = XG(7);
	//if (GGLstab && 1)
	//{
	//	//add stabilisation terms: d(g)/dq * eta, eta is additional Lagrange Multiplier
	//	double eta = -GetMBS()->GetStepSize()*2.*XG(16);

	//	beta0 += beta0*eta; //negative sign because dCqTLambda is put on rhs
	//	beta1 += beta1*eta;
	//	beta2 += beta2*eta;
	//	beta3 += beta3*eta;
	//}
}

void Rigid3D::GetBetaD(Vector& beta) const
{
	for(int i=1; i<=NRotParam();i++)
	{
		beta(i) = XGD(i+3);
	}
}

void Rigid3D::GetBetaInitD(Vector& beta) const
{
	for(int i=1; i<=NRotParam();i++)
	{
		beta(i) = x_init(i+3);
	}
}

void Rigid3D::GetBetaPD(Vector& betap) const
{
	for(int i=1; i<=NRotParam();i++)
	{
		betap(i) = XGPD(i+3);	
	}
	//betap(1) = XGPD(4); betap(2) = XGPD(5); betap(3) = XGPD(6); betap(4) = XGPD(7);
}


Matrix3D Rigid3D::GetG() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);
	beta*=2;// beta1*=2; beta2*=2; beta3*=2;
	//the betas are already multiplied with 2, compared to Shabana
	
	/*//construct consistent beta:
	double fact = 2./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
	beta0 *= fact;
	beta1 *= fact;
	beta2 *= fact;
	beta3 *= fact;*/

	//return Matrix3D(
	//	-beta1, beta0,-beta3, beta2,
	//	-beta2, beta3, beta0,-beta1,
	//	-beta3,-beta2, beta1, beta0);

		return Matrix3D(
		-beta(2), beta(1),-beta(4), beta(3),
		-beta(3), beta(4), beta(1),-beta(2),
		-beta(4),-beta(3), beta(2), beta(1));
}

Matrix3D Rigid3D::GetGT() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);
	beta*=2;
	//the betas are already multiplied with 2, compared to Shabana

	/*//construct consistent beta:
	double fact = 2./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
	beta0 *= fact;
	beta1 *= fact;
	beta2 *= fact;
	beta3 *= fact;*/

	Matrix3D GT;
	GT.Set43(
		-beta(2),-beta(3),-beta(4),
		beta(1), beta(4),-beta(3),
		-beta(4), beta(1), beta(2),
		beta(3),-beta(2), beta(1));
	return GT;
}

Matrix3D Rigid3D::GetGbar() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);
	beta*=2;


	//construct consistent beta:
	/*double fact = 2./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
	beta0 *= fact;
	beta1 *= fact;
	beta2 *= fact;
	beta3 *= fact;*/


	//the betas are already multiplied with 2, compared to Shabana
	return Matrix3D(
		-beta(2), beta(1), beta(4),-beta(3),
		-beta(3),-beta(4), beta(1), beta(2),
		-beta(4), beta(3),-beta(2), beta(1));
}

Matrix3D Rigid3D::GetGbarT() const
{
	ConstVector<4> beta(NRotParam());
	GetBeta(beta);
	beta*=2;

	////construct consistent beta:
	//double fact = 2./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
	//beta0 *= fact;
	//beta1 *= fact;
	//beta2 *= fact;
	//beta3 *= fact;


	//the betas are already multiplied with 2, compared to Shabana
	Matrix3D Gbar;
	Gbar.Set43(
		-beta(2),-beta(3),-beta(4),
		beta(1),-beta(4), beta(3),
		beta(4), beta(1),-beta(2),
		-beta(3), beta(2), beta(1));
	return Gbar;
}

Matrix3D Rigid3D::GetGbarp() const
{
	ConstVector<4> betap(NRotParam());
	GetBetaP(betap);
	betap*=2;

	return Matrix3D(
		-betap(2), betap(1), betap(4),-betap(3),
		-betap(3),-betap(4), betap(1), betap(2),
		-betap(4), betap(3),-betap(2), betap(1));
}

Matrix3D Rigid3D::GetGbarpT() const
{
	ConstVector<4> betap(NRotParam());
	GetBetaP(betap);
	betap*=2;

	//the betaps are already multiplied with 2, compared to Shabana
	Matrix3D Gbar;
	Gbar.Set43(
		-betap(2),-betap(3),-betap(4),
		betap(1),-betap(4), betap(3),
		betap(4), betap(1),-betap(2),
		-betap(3), betap(2), betap(1));
	return Gbar;
}

// d omega / d theta = - dot Gbar
Matrix3D Rigid3D::GetDOmegaDTheta(Vector& betap) const
{
	return -2*Matrix3D(
		-betap(2), betap(1), betap(4),-betap(3),
		-betap(3),-betap(4), betap(1), betap(2),
		-betap(4), betap(3),-betap(2), betap(1));
}


//C_q^T*\lambda for Euler parameter equation:
void Rigid3D::AddEPCqTterms(Vector& f)
{
	if(NRotParam() == 4)
	{
		ConstVector<4> beta(NRotParam());
		GetBeta(beta);
	
		double v = 2*GetLagrangeMultEP();//*GetMass();
	
		f(GetIndBeta(1)) -= v*beta(1);
		f(GetIndBeta(2)) -= v*beta(2);
		f(GetIndBeta(3)) -= v*beta(3);
		f(GetIndBeta(4)) -= v*beta(4);
	}
}

// old:
//Matrix3D Rigid3D::ComputeRotMatrix(const double& beta0, const double& beta1,
//																	 const double& beta2, const double& beta3) const
//{
//	double fact;
//	if (GGLStabilisation()) fact = 1;
//	else fact = 1./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
//
//	double b0 = beta0*fact;
//	double b1 = beta1*fact;
//	double b2 = beta2*fact;
//	double b3 = beta3*fact;
//
//	Matrix3D rot;
//
//	rot.Get0(0,0)= -2.0*b3*b3-2.0*b2*b2+1.0;
//	rot.Get0(0,1)= -2.0*b3*b0+2.0*b2*b1;
//	rot.Get0(0,2)= 2.0*b3*b1+2.0*b2*b0;
//	rot.Get0(1,0)= 2.0*b3*b0+2.0*b2*b1;
//	rot.Get0(1,1)= -2.0*b3*b3-2.0*b1*b1+1.0;
//	rot.Get0(1,2)= 2.0*b3*b2-2.0*b1*b0;
//	rot.Get0(2,0)= -2.0*b2*b0+2.0*b3*b1;
//	rot.Get0(2,1)= 2.0*b3*b2+2.0*b1*b0;
//	rot.Get0(2,2)= -2.0*b2*b2-2.0*b1*b1+1.0;
//	return rot;
//}
//
//Matrix3D Rigid3D::ComputeRotMatrixP(const double& beta0, const double& beta1,
//																		const double& beta2, const double& beta3, const double& betap0, const double& betap1,
//																		const double& betap2, const double& betap3) const
//{
//
//	double fact = 1.;//1./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
//	double b0 = beta0*fact;
//	double b1 = beta1*fact;
//	double b2 = beta2*fact;
//	double b3 = beta3*fact;
//	/*
//	double b0 = beta0;
//	double b1 = beta1;
//	double b2 = beta2;
//	double b3 = beta3;
//	*/
//	Matrix3D rotp;
//	rotp.Get0(0,0)= -4.0*b3*betap3-4.0*b2*betap2;
//	rotp.Get0(0,1)= -2.0*betap3*b0-2.0*b3*betap0+2.0*betap2*b1+2.0*b2*betap1;
//	rotp.Get0(0,2)= 2.0*betap3*b1+2.0*b3*betap1+2.0*betap2*b0+2.0*b2*betap0;
//	rotp.Get0(1,0)= 2.0*betap3*b0+2.0*b3*betap0+2.0*betap2*b1+2.0*b2*betap1;
//	rotp.Get0(1,1)= -4.0*b3*betap3-4.0*b1*betap1;
//	rotp.Get0(1,2)= 2.0*betap3*b2+2.0*b3*betap2-2.0*betap1*b0-2.0*b1*betap0;
//	rotp.Get0(2,0)= -2.0*betap2*b0-2.0*b2*betap0+2.0*betap3*b1+2.0*b3*betap1;
//	rotp.Get0(2,1)= 2.0*betap3*b2+2.0*b3*betap2+2.0*betap1*b0+2.0*b1*betap0;
//	rotp.Get0(2,2)= -4.0*b2*betap2-4.0*b1*betap1;
//	return rotp;
//}

//new:
Matrix3D Rigid3D::ComputeRotMatrix(const Vector& beta) const
{
	double fact;
	//if (GGLStabilisation()) fact = 1;
	//else
  fact = 1./sqrt(beta(1)*beta(1)+beta(2)*beta(2)+beta(3)*beta(3)+beta(4)*beta(4)); //JG2013-01-21: maybe this can be erased?

	double b0 = beta(1)*fact;
	double b1 = beta(2)*fact;
	double b2 = beta(3)*fact;
	double b3 = beta(4)*fact;

	Matrix3D rot;

	rot.Get0(0,0)= -2.0*b3*b3-2.0*b2*b2+1.0;
	rot.Get0(0,1)= -2.0*b3*b0+2.0*b2*b1;
	rot.Get0(0,2)= 2.0*b3*b1+2.0*b2*b0;
	rot.Get0(1,0)= 2.0*b3*b0+2.0*b2*b1;
	rot.Get0(1,1)= -2.0*b3*b3-2.0*b1*b1+1.0;
	rot.Get0(1,2)= 2.0*b3*b2-2.0*b1*b0;
	rot.Get0(2,0)= -2.0*b2*b0+2.0*b3*b1;
	rot.Get0(2,1)= 2.0*b3*b2+2.0*b1*b0;
	rot.Get0(2,2)= -2.0*b2*b2-2.0*b1*b1+1.0;
	return rot;
}

Matrix3D Rigid3D::ComputeRotMatrixP(const Vector& beta, const Vector& betap) const
{

	double fact = 1.;//1./sqrt(beta0*beta0+beta1*beta1+beta2*beta2+beta3*beta3);
	double b0 = beta(1)*fact;
	double b1 = beta(2)*fact;
	double b2 = beta(3)*fact;
	double b3 = beta(4)*fact;

	Matrix3D rotp;
	rotp.Get0(0,0)= -4.0*b3*betap(4)-4.0*b2*betap(3);
	rotp.Get0(0,1)= -2.0*betap(4)*b0-2.0*b3*betap(1)+2.0*betap(3)*b1+2.0*b2*betap(2);
	rotp.Get0(0,2)= 2.0*betap(4)*b1+2.0*b3*betap(2)+2.0*betap(3)*b0+2.0*b2*betap(1);
	rotp.Get0(1,0)= 2.0*betap(4)*b0+2.0*b3*betap(1)+2.0*betap(3)*b1+2.0*b2*betap(2);
	rotp.Get0(1,1)= -4.0*b3*betap(4)-4.0*b1*betap(2);
	rotp.Get0(1,2)= 2.0*betap(4)*b2+2.0*b3*betap(3)-2.0*betap(2)*b0-2.0*b1*betap(1);
	rotp.Get0(2,0)= -2.0*betap(3)*b0-2.0*b2*betap(1)+2.0*betap(4)*b1+2.0*b3*betap(2);
	rotp.Get0(2,1)= 2.0*betap(4)*b2+2.0*b3*betap(3)+2.0*betap(2)*b0+2.0*b1*betap(1);
	rotp.Get0(2,2)= -4.0*b2*betap(3)-4.0*b1*betap(2);
	return rotp;


}

//for body loads:
//Computes f = d p_ref/d q * x
void Rigid3D::ApplyDprefdq(Vector& f, const Vector3D& x)
{
	f(1) = x(1);
	f(2) = x(2);
	f(3) = x(3);
	for(int i=4; i<=3+NRotParam();i++)
	{
		f(i) = 0.;	
	}
	//f(4) = 0;
	//f(5) = 0;
	//f(6) = 0;
	//f(7) = 0; //change for quaternions!!!
}
//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
void Rigid3D::ApplyDrotrefdq(Vector& f, const Vector3D& x)
{
	f(1) = 0;
	f(2) = 0;
	f(3) = 0;
	static Vector fh; fh.LinkWith(f,4,NRotParam()); //change for quaternions!!!
	//static Vector fh; fh.SetLen(4);
	Mult(GetGT(),x,fh);
	/*f(4) = fh(1);
	f(5) = fh(2);
	f(6) = fh(3);
	f(7) = fh(4);*/
}
//only displacements, rotations make no sense, even in rigid body
//->only for volumeloads (gravity ...)
void Rigid3D::GetDuxDq(Vector& dudq)
{
	dudq.SetAll(0);
	dudq(1) = 1;
}
void Rigid3D::GetDuyDq(Vector& dudq)
{
	dudq.SetAll(0);
	dudq(2) = 1;
}

void Rigid3D::GetDuzDq(Vector& dudq)
{
	dudq.SetAll(0);
	dudq(3) = 1;
}

void Rigid3D::GetH3T(const Vector3D& vloc, Matrix3D& d)
{
	//realize d=(skew(pglob)^T*G)^T = G^T*skew(A*ploc)
	Matrix3D m;
	//m.SetSkew(GetRotMatrix()*vloc);
	//d = GetGT()*m;
	m.SetSkew(vloc);
	d = GetGbarT()*m*GetRotMatrix().GetTp();
}

void Rigid3D::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	// -> simple Sparse matrix mit set(i,j)/get(ii,j)/m*v, SetSubMat/AddSubmat, not ordered!!!!, Clear, dim
	d.SetSize(3+NRotParam(),3);
	for (int i=1; i<=3; i++)
		for (int j=1; j<=3; j++)
			d(i,j) = 0;

	Matrix3D H;
	GetH3T(vloc,H);
	d.SetSubmatrix(H,4,1);
}

void Rigid3D::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	// A screw(ploc) A^T = screw(pglob) (Shabana 1994, p.372)
	//
	// L^T =  | I                        |  =   |   I              |
	//        | (-A screw(ploc) Gbar)^T  |      | G^T screw(pglob) |
	//
	//  (-A screw(ploc) Gbar)^T = (-A screw(ploc) A^T G)^T = (-screw(pglob) G)^T = G^T screw(pglob)
	//	G^T*skew(A*ploc) = -G	(-A screw(ploc) Gbar)^T
	d.SetSize(3+NRotParam(),3);
	d(1,1)=1; d(1,2)=0; d(1,3)=0;
	d(2,1)=0; d(2,2)=1; d(2,3)=0;
	d(3,1)=0; d(3,2)=0; d(3,3)=1;
	Matrix3D G;
	GetH3T(ploc, G);

	d.SetSubmatrix(G,4,1);
}

void Rigid3D::GetdAngVeldqpT(const Vector3D& ploc, Matrix& d)
{
	d.SetSize(3+NRotParam(),3);
	d(1,1)=0; d(1,2)=0; d(1,3)=0;
	d(2,1)=0; d(2,2)=0; d(2,3)=0;
	d(3,1)=0; d(3,2)=0; d(3,3)=0;

	d.SetSubmatrix(GetGT(),4,1);
}

void Rigid3D::GetdRotdqT(const Vector3D &ploc, Matrix &d)
{
	d.SetSize(3+NRotParam(),3);
	d(1,1)=0; d(1,2)=0; d(1,3)=0;
	d(2,1)=0; d(2,2)=0; d(2,3)=0;
	d(3,1)=0; d(3,2)=0; d(3,3)=0;

	d.SetSubmatrix(GetGT(),4,1);
}

void Rigid3D::GetdPosdx(const Vector3D& ploc, Vector3D& dpdx)
{
	dpdx = GetRotMatrix()*Vector3D(1,0,0);
}

//$ DR 2012-11-02 deprecated, use ReadSingleElementData instead
double Rigid3D::GetSpecialSensorValue(int nr, double time) const
{
	if (1)
	{
		/*ConstVector<4> theta_c, thetap_c;
		Vector3D theta, thetap;
		this->GetBeta(theta_c);
 
		//Vector3D omega = GetAngularVel();
		Vector3D omegabar = GetAngularVelLocal();
		Vector3D L = GetRotMatrix() * (Iphi*omegabar);

		if (nr == 1)
			return L(1);
		if (nr == 2)
			return L(2);
		if (nr == 3)
			return L(3);*/
		if (nr == 4)
			return GetAngularMomentum().X();
		if (nr == 5)
			return GetAngularMomentum().Y();
		if (nr == 6)
			return GetAngularMomentum().Z();
	}
	if (nr == 7)
	{
		return const_cast<Rigid3D*>(this)->GetKineticEnergy();
	}
	if (nr == 8)
	{
		return const_cast<Rigid3D*>(this)->GetPotentialEnergy();
	}
	
	return 0.;
}

int Rigid3D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Body3D::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	Rigid3D::GetAvailableSpecialValuesAuto(available_variables);

	// Manual entries for this class
	return 0;
}

int Rigid3D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Element::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read  
	if( RWdata.variable_name == mystr("Internal.angular_momentum") )
	{
		if (RWdata.comp1 ==1)
		{
			RWdata.value = GetAngularMomentum().X(); return 1; 
		}
		else if (RWdata.comp1 ==2)
		{
			RWdata.value = GetAngularMomentum().Y(); return 1; 
		}
		else if (RWdata.comp1 ==3)
		{
			RWdata.value = GetAngularMomentum().Z(); return 1; 
		}
		else return -2; 
	}
	else if ( RWdata.variable_name == mystr("Internal.kinetic_energy") )
	{
		RWdata.value = const_cast<Rigid3D*>(this)->GetKineticEnergy(); return 1;
	}
	else if ( RWdata.variable_name == mystr("Internal.potential_energy") )
	{
		RWdata.value = const_cast<Rigid3D*>(this)->GetPotentialEnergy(); return 1;
	}

	return ReadSingleElementDataAuto(RWdata);
}

void Rigid3D::DrawElement() 
{
	Body3D::DrawElement();

	////draw local frame at origin:
	//if (GetMBS()->GetIOption(125))
	//{
	//	double s = GetMBS()->GetDOption(104);

	//	GetMBS()->ChooseColor(0.3f,0.3f,0.3f);

	//	Vector3D v1(-0.2*s, 0,0);
	//	Vector3D v2( s, 0,0);
	//	Vector3D v3( 0,-0.2*s,0);
	//	Vector3D v4( 0, s,0);
	//	Vector3D v5( 0,0,-0.2*s);
	//	Vector3D v6( 0,0, s);
	//	v1 = GetPosD(v1);
	//	v2 = GetPosD(v2);
	//	v3 = GetPosD(v3);
	//	v4 = GetPosD(v4);
	//	v5 = GetPosD(v5);
	//	v6 = GetPosD(v6);
	//	double d = GetMBS()->GetDOption(114);
	//	GetMBS()->MyDrawLine(v1,v2,d);
	//	GetMBS()->MyDrawLine(v3,v4,d);
	//	GetMBS()->MyDrawLine(v5,v6,d);

	//	char str[20];
	//	sprintf(str, "X%d", GetOwnNum());
	//	GetMBS()->GetRC()->PrintText3D((float)v2.X(), (float)v2.Y(), (float)v2.Z(), str);
	//	sprintf(str, "Y%d", GetOwnNum());
	//	GetMBS()->GetRC()->PrintText3D((float)v4.X(), (float)v4.Y(), (float)v4.Z(), str);
	//	sprintf(str, "Z%d", GetOwnNum());
	//	GetMBS()->GetRC()->PrintText3D((float)v6.X(), (float)v6.Y(), (float)v6.Z(), str);
	//}

	mbs->SetColor(col);

	//mbs->uout << "Rigid3D: m=" << mass << ", I=" << Iphi << ", Size=" << size << "\n";
	//GetMBS()->UO() << "r_posD=" << GetPosD(Vector3D(1,1,1)) << "\n";

	double lx = 0.5*size.X(); double ly = 0.5*size.Y(); double lz = 0.5*size.Z();
	Vector3D p8(GetPosD(Vector3D(-lx,-ly,-lz)));
	Vector3D p7(GetPosD(Vector3D(-lx,-ly, lz)));
	Vector3D p6(GetPosD(Vector3D( lx,-ly,-lz)));
	Vector3D p5(GetPosD(Vector3D( lx,-ly, lz)));
	Vector3D p4(GetPosD(Vector3D(-lx, ly,-lz)));
	Vector3D p3(GetPosD(Vector3D(-lx, ly, lz)));
	Vector3D p2(GetPosD(Vector3D( lx, ly,-lz)));
	Vector3D p1(GetPosD(Vector3D( lx, ly, lz)));

	if (mbs->GetIOption(134)) //draw body faces
	{
		mbs->DrawHex(p1,p2,p3,p4,p5,p6,p7,p8);
	}

	if (mbs->GetIOption(133)) //draw body outline
	{
		double th = 1+mbs->GetDOption(115); //rigid body line thickness
		mbs->MyDrawLine(p1,p3, th);
		mbs->MyDrawLine(p3,p4, th);
		mbs->MyDrawLine(p4,p2, th);
		mbs->MyDrawLine(p2,p1, th);

		mbs->MyDrawLine(p5,p7, th);
		mbs->MyDrawLine(p7,p8, th);
		mbs->MyDrawLine(p8,p6, th);
		mbs->MyDrawLine(p6,p5, th);

		mbs->MyDrawLine(p1,p5, th);
		mbs->MyDrawLine(p3,p7, th);
		mbs->MyDrawLine(p4,p8, th);
		mbs->MyDrawLine(p2,p6, th);
	}
};

//$ DR 2012-10 added
double Rigid3D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{	
	// special imlementation of specific class
	if (flagD)
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:		
			return fvd.GetComponent(GetPosD(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVelD(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetDisplacementD(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2012-12-12 added according to JG
			return fvd.GetComponent(GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2012-12-12 added according to JG
			return fvd.GetComponent(GetRotMatrixD().GetTp()*GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2012-12-12 added according to JG
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrixD(), phi);
				return fvd.GetComponent(phi);
			}

			case FieldVariableDescriptor::FVT_acceleration: //$ MS 2013-9-4: Added
			return fvd.GetComponent(GetAcceleration(local_position));
		}
	}
	else
	{
		switch (fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_position:	
			return fvd.GetComponent(GetPos(local_position));
		case FieldVariableDescriptor::FVT_velocity:
			return fvd.GetComponent(GetVel(local_position));
		case FieldVariableDescriptor::FVT_displacement:
			return fvd.GetComponent(GetDisplacement(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity:		//$ DR 2012-12-12 added according to JG
			return fvd.GetComponent(GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_angular_velocity_local_basis: //$ DR 2012-12-12 added according to JG
			return fvd.GetComponent(GetRotMatrix().GetTp()*GetAngularVel(local_position));
		case FieldVariableDescriptor::FVT_bryant_angle:			//$ DR 2012-12-12 added according to JG
			{
				Vector3D phi;
				RotMatToKardanAngles(GetRotMatrix(), phi);
				return fvd.GetComponent(phi);
			}
			case FieldVariableDescriptor::FVT_acceleration: //$ MS 2013-9-4: Added
			return fvd.GetComponent(GetAcceleration(local_position));
		}
	}
	return FIELD_VARIABLE_NO_VALUE;
}

//$ DR 2012-10 added
void Rigid3D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// add all the field variables of the parent class
	Body3D::GetAvailableFieldVariables(variables);
	//FVT_position,FVT_displacement, FVT_velocity: implemented in Element
	// possibility to remove entries does not exist yet, just do not use the Function of the parent class

	// add some specific ones of this class
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_bryant_angle, FieldVariableDescriptor::FVCI_z); //$ DR 2012-12-12 added according to JG
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity, FieldVariableDescriptor::FVCI_z); //$ DR 2012-12-12 added according to JG
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_angular_velocity_local_basis, FieldVariableDescriptor::FVCI_z); //$ DR 2012-12-12 added according to JG
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_acceleration,FieldVariableDescriptor::FVCI_z); //$ MS 2013-9-4: Added
	//	//FVT_displacement_local_basis,						// vector, 3 components in the local basis
	//	//FVT_velocity_local_basis,								// vector, 3 components in the local basis
	//	//FVT_acceleration_local_basis,						// vector, 3 components in the local basis
}
