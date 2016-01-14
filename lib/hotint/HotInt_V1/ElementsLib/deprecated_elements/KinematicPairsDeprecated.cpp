//#**************************************************************
//#
//# filename:             KinematicPairsDeprecated.cpp
//#
//# author:               Gerstmayr Johannes, Reischl Daniel
//#
//# generated:						29. August 2013

//# description:          This class contains only DEPRECATED kinematic pairs. This file exists only for backward compatibilty. DO NOT USE THE CLASSES FOR NEW MODELS.
//#												Kinematic pairs can be divided in lower (surface contact) and higher (point or line contact) kinematic pairs. 
//#												In this file the lower kinematic pairs and the rigid joints (all d.o.f. are constrained) are provided.
//#												Other implemented constraints are in SpecialConstraints.h
//#                       
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
 
//# Rigid:				all d.o.f. are constrained, e.g. welding

//# lower kinematic pairs:

//# Spherical:		all translatory d.o.f. are constrained
//# Rotary:				all translatory and 2 rotatory d.o.f. are constrained, only rotation about 1 axis possible
//# Cylindrical:	like Rotary, but additionally translation along rotational-axis possible
//# Translatory:	all rotary and 2 translatory d.o.f are constrained
//# Plane:				not implemented
//# Skrew-type:		not implemented

#include "kinematicpairsDeprecated.h"

#include "element.h"
#include "constraint.h"
#include "elementdataaccess.h"
#include "graphicsconstants.h"
//#include "rendercontext.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical			
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical: CylindricalJointOLD 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CylindricalJointOLD::GetElementData(ElementDataContainer& edc) 		// data to edc
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt(GetDrawSizeResolution(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);

	ed.SetDouble(GetDrawSizeAxisLength(), "Draw_axis_length"); ed.SetToolTipText("Drawing length along axis of joint"); edc.Add(ed);

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Global_joint_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} 
	else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Global_joint_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
}

int CylindricalJointOLD::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{																																			// data from edc to variable
	int rv = Constraint::SetElementData(edc);
	
	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0);
	SetDrawSizeScalar(tmp);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) SetDrawSizeResolution(dd);
	GetElemDataDouble(GetMBS(), edc, "Draw_axis_length", tmp);
	SetDrawSizeAxisLength(tmp);

	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis", loccoords(2));
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis", loccoords(3));
	}

	return rv;
}


void CylindricalJointOLD::Initialize() 
{
	if (elements.Length()==1)
	{
		Vector3D lpos = loccoords(1);
		Vector3D vrot = loccoords(2);
		Vector3D ln1, lt1;
		vrot.SetNormalBasis(ln1,lt1);

		Matrix3D RT=GetBody3D(1).GetRotMatrix(lpos).GetTp();
		//UO() << "RT=" << RT << "\n";
		//UO() << "vrot=" << vrot << "\n";

		ln1 = RT*ln1;
		lt1 = RT*lt1;

		loccoords(3)=ln1;
		loccoords(4)=lt1;
	} else
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		Vector3D vrot = loccoords(3);
		Vector3D vrot2 = loccoords(3);
		Matrix3D RT=GetBody3D(1).GetRotMatrix(lp1).GetTp();
		vrot = RT*vrot;
		vrot.Normalize();
		loccoords(6) = vrot;

		RT=GetBody3D(2).GetRotMatrix(lp2).GetTp();
		vrot2 = RT*vrot2;

		Vector3D ln2, lt2;
		vrot2.SetNormalBasis(ln2,lt2);

		loccoords(4) = ln2;
		loccoords(5) = lt2;
	}
};

void CylindricalJointOLD::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: CylindricalJointOLD::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos = loccoords(1);
			const Vector3D& vrot = loccoords(2);
			const Vector3D& ln = loccoords(3);
			const Vector3D& lt = loccoords(4);
			Vector3D v = GetBody3D(1).GetPos(lpos)-p_global;
			//UO() << "pg=" << p_global << ", p1=" << GetBody3D(1).GetPos(lpos) << "\n";
			Matrix3D A=GetBody3D(1).GetRotMatrix(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(1) = v*(A*ln); //constrain relative motion normal to rotation axis 
			f(2) = v*(A*lt); 

			f(3) = vrot*(A*ln);
			f(4) = vrot*(A*lt);
			//UO() << "f=" << f << "\n";
		}
		else
			if (elements.Length()==2) 
			{
				const Vector3D& lp1 = loccoords(1);
				const Vector3D& lp2 = loccoords(2);
				const Vector3D& lr1 = loccoords(6);
				const Vector3D& ln2 = loccoords(4);
				const Vector3D& lt2 = loccoords(5);
				Vector3D v = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
				//vrot1^T*A(1)^T*A(2)*ln2=0, vrot1^T*A(1)^T*A(2)*lt2=0
				Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);

				f(1) = v*(A2*ln2); //constrain relative motion normal to rotation axis 
				f(2) = v*(A2*lt2); 
				v = A1*lr1; 
				f(3) = v*(A2*ln2);
				f(4) = v*(A2*lt2);
			}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos = loccoords(1);
			const Vector3D& vrot = loccoords(2);
			const Vector3D& ln = loccoords(3);
			const Vector3D& lt = loccoords(4);

			Vector3D x = GetBody3D(1).GetPos(lpos)-p_global;
			Vector3D v = GetBody3D(1).GetVel(lpos);

			Matrix3D A=GetBody3D(1).GetRotMatrix(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D Ap=GetBody3D(1).GetRotMatrixP(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(1) = v*(A*ln) + x*(Ap*ln); //constrain relative motion normal to rotation axis
			f(2) = v*(A*lt) + x*(Ap*lt); 

			f(3) = vrot*(Ap*ln);
			f(4) = vrot*(Ap*lt);
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			const Vector3D& lr1 = loccoords(6);
			const Vector3D& ln2 = loccoords(4);
			const Vector3D& lt2 = loccoords(5);
			Vector3D x = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
			Vector3D v = GetBody3D(1).GetVel(lp1)-GetBody3D(2).GetVel(lp2);
			
			//vrot1^T*A(1)^T*A(2)*ln2=0, vrot1^T*A(1)^T*A(2)*lt2=0
			Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A1p=GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2p=GetBody3D(2).GetRotMatrixP(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(1) = v*(A2*ln2) + x*(A2p*ln2); //constrain relative motion normal to rotation axis 
			f(2) = v*(A2*lt2) + x*(A2p*lt2); 

			v = A1*lr1; 
			Vector3D vp = A1p*lr1; 
			f(3) = v*(A2p*ln2)+vp*(A2*ln2);
			f(4) = v*(A2p*lt2)+vp*(A2*lt2);
			//UO() << "fg=" << f << "\n";
		}
	}
};

void CylindricalJointOLD::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc

	hmat.SetSize(f.Length(),4);

	if (elements.Length() == 1) //--> ignore locelemind, only first element!
	{ //lpos,ln,lt
		const Vector3D& lpos = loccoords(1);
		const Vector3D& lvj = loccoords(2);
		const Vector3D& lv1i = loccoords(3);
		const Vector3D& lv2i = loccoords(4);

		Vector3D rpij = GetBody3D(1).GetPos(lpos)-p_global;
		Matrix3D Ai = GetBody3D(1).GetRotMatrix(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);

		Vector3D v1i = Ai*lv1i;
		Vector3D v2i = Ai*lv2i;

		//acc. to Comp.DynShabana 1998, p. 434, Hj=0

		GetBody3D(1).GetdPosdqT(lpos,dpdq);
		Mult(dpdq,v1i,hvec);
		hmat.SetColVec(hvec,1);
		Mult(dpdq,v2i,hvec);
		hmat.SetColVec(hvec,2);

		GetBody3D(1).GetdRotvdqT(lv1i,lpos,dpdq);
		Mult(dpdq,lvj,hvec); //hvec *= -1 ???
		hmat.SetColVec(hvec,3);
		Mult(dpdq,rpij,hvec);
		hmat.AddColVec(1,hvec);

		GetBody3D(1).GetdRotvdqT(lv2i,lpos,dpdq);
		Mult(dpdq,lvj,hvec);
		hmat.SetColVec(hvec,4);
		Mult(dpdq,rpij,hvec);
		hmat.AddColVec(2,hvec);

	}
	else
	{
		//acc. to Comp.DynShabana 1998, p. 434, but here, body i and j are INTERCHANGED!!!!

		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lvj = loccoords(6); //body1=j, rotation axis
		const Vector3D& lv1i = loccoords(4); //body2=i, normal vector1
		const Vector3D& lv2i = loccoords(5); //body2=i, normal vector2

		Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
			
		Matrix3D Aj=GetBody3D(1).GetRotMatrix(lp1);
		Matrix3D Ai=GetBody3D(2).GetRotMatrix(lp2);

		Vector3D v1i = Ai*lv1i;
		Vector3D v2i = Ai*lv2i;
		Vector3D vj = Aj*lvj;

		if (locelemind==2) //bodyi = body2
		{
			GetBody3D(2).GetdPosdqT(lp2,dpdq);
			Mult(dpdq,v1i,hvec);
			hmat.SetColVec(hvec,1);
			Mult(dpdq,v2i,hvec);
			hmat.SetColVec(hvec,2);

			GetBody3D(2).GetdRotvdqT(lv1i,lp2,dpdq);
			Mult(dpdq,vj,hvec);
			hmat.SetColVec(hvec,3);
			Mult(dpdq,rpij,hvec);
			hmat.AddColVec(1,hvec);

			GetBody3D(2).GetdRotvdqT(lv2i,lp2,dpdq);
			Mult(dpdq,vj,hvec);
			hmat.SetColVec(hvec,4);
			Mult(dpdq,rpij,hvec);
			hmat.AddColVec(2,hvec);
		}
		else //body j: = body1
		{
			GetBody3D(1).GetdPosdqT(lp1,dpdq);
			Mult(dpdq,v1i,hvec); hvec *= -1;
			hmat.SetColVec(hvec,1);
			Mult(dpdq,v2i,hvec); hvec *= -1;
			hmat.SetColVec(hvec,2);

			GetBody3D(1).GetdRotvdqT(lvj,lp1,dpdq);
			Mult(dpdq,v1i,hvec);
			hmat.SetColVec(hvec,3);

			GetBody3D(1).GetdRotvdqT(lvj,lp1,dpdq);
			Mult(dpdq,v2i,hvec);
			hmat.SetColVec(hvec,4);
		}
	}
	//UO() << "hmat=" << hmat;
	hvec.SetLen(4);
	hvec(1) = XG(1);
	hvec(2) = XG(2);
	hvec(3) = XG(3);
	hvec(4) = XG(4);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*hvec(1)+hmat(i,2)*hvec(2)+hmat(i,3)*hvec(3)+hmat(i,4)*hvec(4);
	}
};

Vector3D CylindricalJointOLD::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void CylindricalJointOLD::DrawElement() 
// drawing options
// draw_dim.X()==diameter
// draw_dim.Y()==axis length
// draw_dim.Z()==draw resolution
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	int res = (int)GetDrawSizeResolution();//$ DR 2011-04-21: old code: int res = (int)draw_dim.Z();
	//if (res == 0) res = 16;

	if (elements.Length()==1) 
	{
		Vector3D rot;
		rot = loccoords(2);
		rot.Normalize();
		rot *= 0.5*GetDrawSizeAxisLength();					//$ DR 2011-04-21: old code: draw_dim.Y();
		Vector3D p = GetBody3D(1).GetPosD(loccoords(1));
		mbs->DrawZyl(p+rot,p-rot,0.5*GetDrawSizeScalar(),res);	//$ DR 2011-04-21: old code: mbs->DrawZyl(p+rot,p-rot,0.5*draw_dim.X(),res);
		mbs->SetColor(colgrey2);
		mbs->DrawZyl(p+1.5*rot,p-1.5*rot,0.25*GetDrawSizeScalar(),res); //$ DR 2011-04-21: old code: mbs->DrawZyl(p+1.5*rot,p-1.5*rot,0.25*draw_dim.X(),res);
	} else
	{
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		Vector3D p2 = GetBody3D(2).GetPosD(loccoords(2));

		Vector3D rot1, rot2;
		rot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*loccoords(6);
		rot1.Normalize();
		rot1 *= 0.5*GetDrawSizeAxisLength();					//$ DR 2011-04-21: old code: draw_dim.Y();
		rot2 = GetBody3D(2).GetRotMatrixD(loccoords(2))*(loccoords(4).Cross(loccoords(5)));
		rot2.Normalize();
		rot2 *= 0.5*GetDrawSizeAxisLength();					//$ DR 2011-04-21: old code: draw_dim.Y();

		mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*GetDrawSizeScalar(),res); //$ DR 2011-04-21: old code: 0.5*draw_dim.X()
		mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*GetDrawSizeScalar(),res); //$ DR 2011-04-21: old code: 0.5*draw_dim.X()
		mbs->SetColor(colgrey2);
		mbs->DrawZyl(p2+1.5*rot2,p2-1.5*rot2,0.25*GetDrawSizeScalar(),res);	//$ DR 2011-04-21: old code: 0.25*draw_dim.X()
	}
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Translatory	Translatory	Translatory	Translatory	Translatory	Translatory				
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Translatory: PrismaticJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PrismaticJointOLD::PrismaticJointOLD(MBS* mbsi, int en1, 
															 const Vector3D& lp1, const Vector3D& gp2,
															 const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
{	
	x_init = Vector(SS());
	GetCol() = coli;
	//draw_dim = ddim;
	SetDrawSizeScalar(ddim.X());
	AddElement(en1);
	loccoords.Add(lp1);
	loccoords.Add(gp2);
	loccoords.Add(Vector3D(0));//will be filled in initialize!!!
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	//UO() << "prismatic joint changed, check sliding axis!\n";
	elementname = GetElementSpec();
};
//relative motion of body and ground along local body1 axis laxis1, lp1 and gp2 may be equal!!
//relative rotation around axis lp1-lp2 is constrained
PrismaticJointOLD::PrismaticJointOLD(MBS* mbsi, int en1, 
															 const Vector3D& lp1, const Vector3D& gp2, const Vector3D& laxis1,
															 const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
{	
	x_init = Vector(SS());
	GetCol() = coli;
	//draw_dim = ddim;
	SetDrawSizeScalar(ddim.X());
	AddElement(en1);
	loccoords.Add(lp1);
	loccoords.Add(gp2);
	loccoords.Add(laxis1);
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));

	elementname = GetElementSpec();
};
//relative motion of two bodies along axis lp1-lp2, lp1 and lp2 must never be equal!!
//relative rotation around axis lp1-lp2 is constrained
PrismaticJointOLD::PrismaticJointOLD(MBS* mbsi, int en1, int en2, 
															 const Vector3D& lp1, const Vector3D& lp2,
															 const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
{	
	x_init = Vector(SS());
	GetCol() = coli;
	//draw_dim = ddim;
	SetDrawSizeScalar(ddim.X());
	AddElement(en1);
	AddElement(en2);
	loccoords.Add(lp1);
	loccoords.Add(lp2);
	loccoords.Add(Vector3D(0)); //will be filled in initialize!!!
//	UO() << "prismatic joint changed, check sliding axis!\n";
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));

	elementname = GetElementSpec();
};
//relative motion of two bodies along local (body1) axis laxis1, lp1 and lp2 may be equal!!
//relative rotation around axis lp1-lp2 is constrained
PrismaticJointOLD::PrismaticJointOLD(MBS* mbsi, int en1, int en2, 
															 const Vector3D& lp1, const Vector3D& lp2, const Vector3D& laxis1, 
															 const Vector3D& ddim, const Vector3D& coli):Constraint(mbsi), loccoords(), dpdq(), hmat(), hvec() 
{	
	x_init = Vector(SS());
	GetCol() = coli;
	//draw_dim = ddim;
	SetDrawSizeScalar(ddim.X());
	AddElement(en1);
	AddElement(en2);
	loccoords.Add(lp1);
	loccoords.Add(lp2);
	loccoords.Add(laxis1);
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));
	loccoords.Add(Vector3D(0));

	elementname = GetElementSpec();
};



void PrismaticJointOLD::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Local_sliding_axis"); edc.Get(edc.Length()).SetToolTipText("Direction of sliding axis, in local body coordinates [X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Local_sliding_axis1"); edc.Get(edc.Length()).SetToolTipText("Direction of sliding axis, in local body coordinates [X, Y, Z]");
	}
}

int PrismaticJointOLD::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size",tmp, 0);
	SetDrawSizeScalar(tmp);
	
	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Local_sliding_axis", loccoords(3));
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Local_sliding_axis1", loccoords(3));
	}

	return rv;
}




void PrismaticJointOLD::Initialize() 
{
	//find orthogonal vectors:
	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	
	if (loccoords(3).Norm() == 0) 
	{
		if (elements.Length() == 1) 
			loccoords(3) = GetBody3D(1).GetRotMatrix(lp1).GetTp()*(GetBody3D(1).GetPos(lp1) - lp2); //lp2 = gp2!!!!
		else
			loccoords(3) = GetBody3D(1).GetRotMatrix(lp1).GetTp()*(GetBody3D(1).GetPos(lp1) - GetBody3D(2).GetPos(lp2));
	}
	
	const Vector3D& laxis1 = loccoords(3); //local body1 axis of the sliding direction


	Vector3D vi1,vi2,vi3;
	Vector3D vj1,vj3;
	Vector3D gp1 = GetBody3D(1).GetPos(lp1);
	Vector3D gp2 = lp2; //for ground joints, local = global
	if (elements.Length() == 2) 
		gp2 = GetBody3D(2).GetPos(lp2);

	//vi3 = gp1 - gp2; //global sliding axis
	vi3 = GetBody3D(1).GetRotMatrix(lp1) * laxis1; //global sliding axis
	//UO() << "Prismatic joint axis =" << vi3 << "\n";

	vj3 = vi3;
	vi3.SetNormalBasis(vi1,vi2);
	vj1 = vi1;

	Matrix3D RTi=GetBody3D(1).GetRotMatrix(lp1).GetTp();
	Matrix3D RTj(1); //Identity matrix
	if (elements.Length() == 2) 
		RTj = GetBody3D(2).GetRotMatrix(lp2).GetTp();
	vi1 = RTi*vi1;
	vi2 = RTi*vi2;
	vi3 = RTi*vi3;
	vj1 = RTj*vj1;
	vj3 = RTj*vj3;

	loccoords(4) = vi1; //like ln2 in UniversalJoint
	loccoords(5) = vi2;	//like lt2 in UniversalJoint
	loccoords(6) = vj1; 
	loccoords(7) = vj3; //like vrot in UniversalJoint

};

void PrismaticJointOLD::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: PrismaticJointOLD::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& gp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(4);
			const Vector3D& lvi2 = loccoords(5);
			const Vector3D& vj1 = loccoords(6);
			const Vector3D& vj3 = loccoords(7);

			Vector3D rpij = GetBody3D(1).GetPos(lp1)-gp2;

			Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A*lvi1;
			Vector3D gvi2 = A*lvi2;

			f(1) = vj3 * gvi1;
			f(2) = vj3 * gvi2;
			f(3) = gvi1 * rpij;
			f(4) = gvi2 * rpij;
			f(5) = gvi2 * vj1;
			//UO() << "f=" << f << "\n";
		}
		else
			if (elements.Length()==2) 
			{
				const Vector3D& lp1 = loccoords(1);
				const Vector3D& lp2 = loccoords(2);
				const Vector3D& lvi1 = loccoords(4);
				const Vector3D& lvi2 = loccoords(5);
				const Vector3D& lvj1 = loccoords(6);
				const Vector3D& lvj3 = loccoords(7);

				Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);

				Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Vector3D gvi1 = A1*lvi1;
				Vector3D gvi2 = A1*lvi2;
				Vector3D gvj3 = A2*lvj3;

				f(1) = gvj3 * gvi1;
				f(2) = gvj3 * gvi2;
				f(3) = gvi1 * rpij;
				f(4) = gvi2 * rpij;
				f(5) = gvi2 * A2*lvj1;
				//UO() << "f=" << f << "\n";
			}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& gp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(4);
			const Vector3D& lvi2 = loccoords(5);
			const Vector3D& vj1 = loccoords(6);
			const Vector3D& vj3 = loccoords(7);

			Vector3D rpij = GetBody3D(1).GetPos(lp1) - gp2;
			Vector3D vpij = GetBody3D(1).GetVel(lp1);

			Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D Ap= GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A*lvi1;
			Vector3D gvi2 = A*lvi2;
			Vector3D vgvi1 = Ap*lvi1;
			Vector3D vgvi2 = Ap*lvi2;

			f(1) = vj3*vgvi1;
			f(2) = vj3*vgvi2;
			f(3) = vgvi1 * rpij + gvi1 * vpij;
			f(4) = vgvi2 * rpij + gvi2 * vpij;
			f(5) = vgvi2 * vj1;
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(4);
			const Vector3D& lvi2 = loccoords(5);
			const Vector3D& lvj1 = loccoords(6);
			const Vector3D& lvj3 = loccoords(7);

			Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
			Vector3D rpijp = GetBody3D(1).GetVel(lp1)-GetBody3D(2).GetVel(lp2);

			Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A1p= GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2p= GetBody3D(2).GetRotMatrixP(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A1*lvi1;
			Vector3D gvi2 = A1*lvi2;
			Vector3D gvj3 = A2*lvj3;
			Vector3D gvi1p = A1p*lvi1;
			Vector3D gvi2p = A1p*lvi2;
			Vector3D gvj3p = A2p*lvj3;

			f(1) = gvj3 * gvi1p + gvj3p * gvi1;
			f(2) = gvj3 * gvi2p + gvj3p * gvi2;
			f(3) = gvi1 * rpijp + gvi1p * rpij;
			f(4) = gvi2 * rpijp + gvi2p * rpij;
			f(5) = gvi2 * A2p*lvj1 + gvi2p * A2*lvj1;
		}
	}
};

void PrismaticJointOLD::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	hmat.SetSize(f.Length(),5);

	if (elements.Length() == 1) //--> ignore locelemind, only first element!
	{
		const Vector3D& lp1 =  loccoords(1);
		const Vector3D& gp2 =  loccoords(2);
		const Vector3D& lvi1 = loccoords(4);
		const Vector3D& lvi2 = loccoords(5);
		const Vector3D& vj1 =  loccoords(6);
		const Vector3D& vj3 =  loccoords(7);

		Vector3D rpij = GetBody3D(1).GetPos(lp1) - gp2;
		Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A*lvi1;
		Vector3D gvi2 = A*lvi2;

		GetBody3D(1).GetdRotvdqT(lvi1,lp1,dpdq); //Hi1
		Mult(dpdq,vj3,hvec);
		hmat.SetColVec(hvec,1);
		Mult(dpdq,rpij,hvec);
		hmat.SetColVec(hvec,3);

		GetBody3D(1).GetdRotvdqT(lvi2,lp1,dpdq); //Hi2
		Mult(dpdq,vj3,hvec);
		hmat.SetColVec(hvec,2);
		Mult(dpdq,rpij,hvec);
		hmat.SetColVec(hvec,4);
		Mult(dpdq,vj1,hvec);
		hmat.SetColVec(hvec,5);

		GetBody3D(1).GetdPosdqT(lp1,dpdq); //Hip
		Mult(dpdq,gvi1,hvec);
		hmat.AddColVec(3,hvec);
		Mult(dpdq,gvi2,hvec);
		hmat.AddColVec(4,hvec);
	}
	else
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lvi1 = loccoords(4);
		const Vector3D& lvi2 = loccoords(5);
		const Vector3D& lvj1 = loccoords(6);
		const Vector3D& lvj3 = loccoords(7);

		if (locelemind==1)
		{
			Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);

			Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A1*lvi1;
			Vector3D gvi2 = A1*lvi2;
			Vector3D gvj1 = A2*lvj1;
			Vector3D gvj3 = A2*lvj3;


			GetBody3D(1).GetdRotvdqT(lvi1,lp1,dpdq); //Hi1
			Mult(dpdq,gvj3,hvec);
			hmat.SetColVec(hvec,1);
			Mult(dpdq,rpij,hvec);
			hmat.SetColVec(hvec,3);

			GetBody3D(1).GetdRotvdqT(lvi2,lp1,dpdq); //Hi2
			Mult(dpdq,gvj3,hvec);
			hmat.SetColVec(hvec,2);
			Mult(dpdq,rpij,hvec);
			hmat.SetColVec(hvec,4);
			Mult(dpdq,gvj1,hvec);
			hmat.SetColVec(hvec,5);

			GetBody3D(1).GetdPosdqT(lp1,dpdq); //Hip
			Mult(dpdq,gvi1,hvec);
			hmat.AddColVec(3,hvec);
			Mult(dpdq,gvi2,hvec);
			hmat.AddColVec(4,hvec);
			//UO() << "hmat1=" << hmat << "\n";
		}
		else
		{
			Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A1*lvi1;
			Vector3D gvi2 = A1*lvi2;


			GetBody3D(2).GetdRotvdqT(lvj3,lp2,dpdq); //Hi1
			Mult(dpdq,gvi1,hvec);
			hmat.SetColVec(hvec,1);
			Mult(dpdq,gvi2,hvec);
			hmat.SetColVec(hvec,2);

			GetBody3D(2).GetdRotvdqT(lvj1,lp2,dpdq); //Hi2
			Mult(dpdq,gvi2,hvec);
			hmat.SetColVec(hvec,5);

			GetBody3D(2).GetdPosdqT(lp2,dpdq); //Hip
			Mult(dpdq,gvi1,hvec); hvec *= -1;
			hmat.SetColVec(hvec,3);
			Mult(dpdq,gvi2,hvec); hvec *= -1;
			hmat.SetColVec(hvec,4);
			//UO() << "hmat2=" << hmat << "\n";
		}
	}
	//UO() << "hmat=" << hmat;
	hvec.SetLen(5);
	hvec(1) = XG(1);
	hvec(2) = XG(2);
	hvec(3) = XG(3);
	hvec(4) = XG(4);
	hvec(5) = XG(5);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*hvec(1)+hmat(i,2)*hvec(2)+hmat(i,3)*hvec(3)+hmat(i,4)*hvec(4)+hmat(i,5)*hvec(5);
	}
};

Vector3D PrismaticJointOLD::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void PrismaticJointOLD::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	//...
	if (GetDrawSizeScalar() == 0) return;
	if (elements.Length()==1) 
	{
		Vector3D rot1, rot2;
		rot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*(loccoords(3).Cross(loccoords(4)));
		rot1.Normalize();
		rot1 *= 0.5*GetDrawSizeScalar();
		rot2 = loccoords(6);
		rot2.Normalize();
		rot2 *= 0.5*GetDrawSizeScalar();
		Vector3D p = GetBody3D(1).GetPosD(loccoords(1));
		
		mbs->MyDrawLine(p+rot1,p-rot1,2);
		mbs->MyDrawLine(p+1.5*rot2,p-1.5*rot2,2);
		//mbs->DrawZyl(p+rot1,p-rot1,0.5*draw_dim.X(),20);
		//mbs->SetColor(colgrey2);
		//mbs->DrawZyl(p+1.5*rot2,p-1.5*rot2,0.07*draw_dim.X(),8);
	} else
	{

		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		Vector3D p2 = GetBody3D(2).GetPosD(loccoords(2));

		Vector3D rot1, rot2;
		rot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*(loccoords(3).Cross(loccoords(4)));
		rot1.Normalize();
		rot1 *= 0.5*GetDrawSizeScalar();
		rot2 = GetBody3D(2).GetRotMatrixD(loccoords(2))*(loccoords(6));
		rot2.Normalize();
		rot2 *= 0.5*GetDrawSizeScalar();

		mbs->MyDrawLine(p1+rot1,p1-rot1,2);
		mbs->MyDrawLine(p2+1.5*rot2,p2-1.5*rot2,2);
		/*
		mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*draw_dim.X(),24);
		mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*draw_dim.X(),24);
		mbs->SetColor(colgrey2);
		mbs->DrawZyl(p2+1.15*rot2,p2-1.15*rot2,0.4*draw_dim.X(),24);
		*/
	}
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rotary		Rotary		Rotary		Rotary		Rotary		Rotary		Rotary		
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rotary: RevoluteJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RevoluteJointOLD::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt(GetDrawSizeResolution(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);

	ed.SetDouble(GetDrawSizeAxisLength(), "Draw_axis_length"); ed.SetToolTipText("Drawing length along axis of joint"); edc.Add(ed);

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Global_joint_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(3), "Global_joint_axis"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
}

int RevoluteJointOLD::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0);
	SetDrawSizeScalar(tmp);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) SetDrawSizeResolution(dd);
	GetElemDataDouble(GetMBS(), edc, "Draw_axis_length", tmp);
	SetDrawSizeAxisLength(tmp);

	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis", loccoords(2));
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_axis", loccoords(3));
	}

	return rv;
}


void RevoluteJointOLD::Initialize() 
{
	if (elements.Length()==1)
	{
		Vector3D lpos = loccoords(1);
		Vector3D vrot = loccoords(2);
		Vector3D ln1, lt1;
		vrot.SetNormalBasis(ln1,lt1);

		Matrix3D RT=GetBody3D(1).GetRotMatrix(lpos).GetTp();

		ln1 = RT*ln1;
		lt1 = RT*lt1;

		loccoords(3)=ln1;
		loccoords(4)=lt1;
	} else
	{
		const Vector3D& lp1 = loccoords(1);  //local position vector 1
		const Vector3D& lp2 = loccoords(2);  //local position vector 2
		Vector3D vrot = loccoords(3);        //rotation axis 1 (unit vector)
		Vector3D vrot2 = loccoords(3);       //rotation axis 2 (unit vector)
		Matrix3D RT=GetBody3D(1).GetRotMatrix(lp1).GetTp(); // A: from inertial to body fixed coordinates => u = A^t . ubar
		vrot = RT*vrot;
		vrot.Normalize();
		loccoords(6) = vrot; 

		RT=GetBody3D(2).GetRotMatrix(lp2).GetTp(); 
		vrot2 = RT*vrot2;   // u = A^t . ubar  //?J Ist Transponiere RT richtig?
		vrot2.Normalize();

		Vector3D ln2, lt2;
		vrot2.SetNormalBasis(ln2,lt2);

		loccoords(4) = ln2; 
		loccoords(5) = lt2; 

	}
};

void RevoluteJointOLD::EvalG(Vector& f, double t) 
{
	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: RevoluteJointOLD::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos = loccoords(1);
			const Vector3D& vrot = loccoords(2);
			const Vector3D& ln = loccoords(3);
			const Vector3D& lt = loccoords(4);
			Vector3D v = GetBody3D(1).GetPos(lpos)-p_global;
			//UO() << "pg=" << p_global << ", p1=" << GetBody3D(1).GetPos(lpos) << "\n";
			f(1) = v(1); f(2) = v(2); f(3) = v(3);
			//vrot^T*A(1)*ln=0, vrot^T*A(1)*lt=0
			Matrix3D A=GetBody3D(1).GetRotMatrix(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(4) = vrot*(A*ln);
			f(5) = vrot*(A*lt);
			//UO() << "f=" << f << "\n";
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			const Vector3D& lr1 = loccoords(6);
			const Vector3D& ln2 = loccoords(4);
			const Vector3D& lt2 = loccoords(5);
			Vector3D v = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);
			//vrot1^T*A(1)^T*A(2)*ln2=0, vrot1^T*A(1)^T*A(2)*lt2=0
			Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			v = A1*lr1; 
			f(4) = v*(A2*ln2);
			f(5) = v*(A2*lt2);
		}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			const Vector3D& lpos = loccoords(1);
			const Vector3D& vrot = loccoords(2);
			const Vector3D& ln = loccoords(3);
			const Vector3D& lt = loccoords(4);

			Vector3D v = GetBody3D(1).GetVel(lpos);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);
			//vrot^T*A(1)*ln=0, vrot^T*A(1)*lt=0
			Matrix3D A=GetBody3D(1).GetRotMatrixP(lpos); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			f(4) = vrot*(A*ln);
			f(5) = vrot*(A*lt);
			/*			UO() << "Ap=\n" << A;
			UO() << "f =" << f << "\n";
			UO() << "vrot =" << vrot << "\n";
			UO() << "vrot =" << vrot << "\n";
			UO() << "ln   =" << ln << "\n";
			UO() << "lt   =" << lt << "\n";*/
		}
		else if (elements.Length()==2) 
		{
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			const Vector3D& lr1 = loccoords(6);
			const Vector3D& ln2 = loccoords(4);
			const Vector3D& lt2 = loccoords(5);
			Vector3D v = GetBody3D(1).GetVel(lp1)-GetBody3D(2).GetVel(lp2);
			f(1) = v(1); f(2) = v(2); f(3) = v(3);
			//vrot1^T*A(1)^T*A(2)*ln2=0, vrot1^T*A(1)^T*A(2)*lt2=0
			Matrix3D A1=GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2=GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A1p=GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2p=GetBody3D(2).GetRotMatrixP(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			v = A1*lr1; 
			Vector3D vp = A1p*lr1; 
			f(4) = v*(A2p*ln2)+vp*(A2*ln2);
			f(5) = v*(A2p*lt2)+vp*(A2*lt2);
			//UO() << "fg=" << f << "\n";
		}
	}
};

void RevoluteJointOLD::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc

	hmat.SetSize(f.Length(),5);

	if (elements.Length() == 1) //--> ignore locelemind, only first element!
	{ //lpos,ln,lt
		const Vector3D& lpos = loccoords(1);
		const Vector3D& vrot = loccoords(2);
		const Vector3D& ln = loccoords(3);
		const Vector3D& lt = loccoords(4);

		GetBody3D(1).GetdPosdqT(lpos,dpdq);
		hmat.SetSubmatrix(dpdq,1,1,-1);

		GetBody3D(1).GetdRotvdqT(ln,lpos,dpdq);
		Mult(dpdq,vrot,hvec); hvec*=-1;
		hmat.SetColVec(hvec,4);

		GetBody3D(1).GetdRotvdqT(lt,lpos,dpdq);
		Mult(dpdq,vrot,hvec); hvec*=-1;
		hmat.SetColVec(hvec,5);
		//UO() << "hmat=" << hmat;
	}
	else
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lr1 = loccoords(6);
		const Vector3D& ln2 = loccoords(4);
		const Vector3D& lt2 = loccoords(5);

		double sign = 1;
		if (locelemind==2) sign = -1;

		//TMStartTimer(20);
		GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
		hmat.SetSubmatrix(dpdq,1,1, sign);
		//TMStopTimer(20);

		//TMStartTimer(21);
		if (locelemind==1)
		{
			Vector3D vn2_glob = GetBody3D(2).GetRotMatrix(lp2)*ln2;
			Vector3D vt2_glob = GetBody3D(2).GetRotMatrix(lp2)*lt2;
			GetBody3D(1).GetdRotvdqT(lr1,lp1,dpdq);
			Mult(dpdq,vn2_glob,hvec);
			hmat.SetColVec(hvec,4);

			Mult(dpdq,vt2_glob,hvec);
			hmat.SetColVec(hvec,5);
		}
		else
		{
			
			Vector3D vrglob = GetBody3D(1).GetRotMatrix(lp1)*lr1;
			/*
			if (GetMBS()->GetTime() > 0.303) 
			{
			double xxxx=5;
			}
			double n1 = vrglob.Norm();
			double n2 = ln2.Norm();
			double n3 = lt2.Norm();*/

			GetBody3D(2).GetdRotvdqT(ln2,lp2,dpdq);
			Mult(dpdq,vrglob,hvec);
			hmat.SetColVec(hvec,4);

			GetBody3D(2).GetdRotvdqT(lt2,lp2,dpdq);
			Mult(dpdq,vrglob,hvec);
			hmat.SetColVec(hvec,5);
		}
		//TMStopTimer(21);
	}
	//UO() << "hmat=" << hmat;
	hvec.SetLen(5);
	hvec(1) = -XG(1);
	hvec(2) = -XG(2);
	hvec(3) = -XG(3);
	hvec(4) = -XG(4);
	hvec(5) = -XG(5);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*hvec(1)+hmat(i,2)*hvec(2)+hmat(i,3)*hvec(3)+hmat(i,4)*hvec(4)+hmat(i,5)*hvec(5);
	}
};

Vector3D RevoluteJointOLD::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void RevoluteJointOLD::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	int res = GetDrawSizeResolution();

	if (elements.Length()==1) 
	{
		Vector3D rot;
		rot = loccoords(2);
		rot.Normalize();
		rot *= 0.5*GetDrawSizeAxisLength();
		Vector3D p = GetBody3D(1).GetPosD(loccoords(1));
		mbs->DrawZyl(p+rot,p-rot,0.5*GetDrawSizeScalar(),res);
		mbs->SetColor(colgrey2);
		mbs->DrawZyl(p+1.5*rot,p-1.5*rot,0.07*GetDrawSizeScalar(),res);
	} else
	{
		Vector3D p1 = GetBody3D(1).GetPosD(loccoords(1));
		Vector3D p2 = GetBody3D(2).GetPosD(loccoords(2));

		Vector3D rot1, rot2;
		rot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*loccoords(6);
		rot1.Normalize();
		rot1 *= 0.5*GetDrawSizeAxisLength();
		rot2 = GetBody3D(2).GetRotMatrixD(loccoords(2))*(loccoords(4).Cross(loccoords(5)));
		rot2.Normalize();
		rot2 *= 0.5*GetDrawSizeAxisLength();

		mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*GetDrawSizeScalar(),res);
		mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*GetDrawSizeScalar(),res);
		mbs->SetColor(colgrey2);
		mbs->DrawZyl(p2+1.2*rot2,p2-1.2*rot2,0.4*GetDrawSizeScalar(),res);
	}
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	Rigid	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Rigid: RigidJointOLD
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RigidJointOLD::GetElementData(ElementDataContainer& edc) 		// data to edc
{
	Constraint::GetElementData(edc);

	//$ DR 2012-10: is already in GetSetElementDataAuto
	//ElementData ed;
	//ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);	

	//if (elements.Length()==1)
	//{
	//	SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	//	SetElemDataVector3D(edc, loccoords(2), "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	//} else
	//{
	//	SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	//	SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	//}
}

int RigidJointOLD::SetElementData(ElementDataContainer& edc) // set element data according to ElementDataContainer
{																																// data from edc to variable
	int rv = Constraint::SetElementData(edc);

	//draw_dim = 0;
	//double tmp;
	//GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0); //$ DR 2012-10: is already in GetSetElementDataAuto
	//SetDrawSizeScalar(tmp);

	//if (elements.Length()==1)
	//{
	//	GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
	//	GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", loccoords(2));
	//} else
	//{
	//	GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
	//	GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	//}

	return rv;
}

void RigidJointOLD::ElementDefaultConstructorInitialization()
{
	x_init = Vector(SS());
	loccoords.SetLen(7);
	loccoords.SetAll(0);
	elementname = GetElementSpec();
	elements.Set2(1,0);
}

void RigidJointOLD::Initialize() 
{
	//find orthogonal vectors:
	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	Vector3D vi1,vi2,vi3;
	Vector3D vj1,vj3;
	
	//Vector3D gp1 = GetBody3D(1).GetPos(lp1);
	//Vector3D gp2 = lp2; //for ground joints, local = global

	//if (elements.Length() == 2) 
	//	gp2 = GetBody3D(2).GetPos(lp2);

	vi3 = Vector3D(1.,0.,0.);
	vj3 = vi3;
	vi3.SetNormalBasis(vi1,vi2);
	vj1 = vi1;

	Matrix3D RTi=GetBody3D(1).GetRotMatrix(lp1).GetTp();
	Matrix3D RTj(1); //Identity matrix
	//if (elements.Length() == 2) 
	if (elements(2))	//$ DR 2012-10 elements.Length is now always 2 
		RTj = GetBody3D(2).GetRotMatrix(lp2).GetTp();
	vi1 = RTi*vi1;
	vi2 = RTi*vi2;
	vi3 = RTi*vi3;
	vj1 = RTj*vj1;
	vj3 = RTj*vj3;

	loccoords(3) = vi1; //body i
	loccoords(4) = vi2;	//body i
	loccoords(5) = vi3;	//body i
	loccoords(6) = vj1; //body j
	loccoords(7) = vj3; //body j

};

void RigidJointOLD::EvalG(Vector& f, double t) 
{
	//$ DR 2012-10 elements.Length is now always 2
	//if (elements.Length()<1 || elements.Length()>2)
	//{
	//	mbs->uout << "ERROR: RigidJointOLD::EvalG, number of elements != 1 or 2\n"; return;
	//}

	if (MaxIndex()==3)
	{
		//if (elements.Length()==1)
		if (elements(2)==0)	//$ DR 2012-10 elements.Length is now always 2
		{
			//ground joint:
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& gp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(3);
			const Vector3D& lvi2 = loccoords(4);
			const Vector3D& lvi3 = loccoords(5);
			const Vector3D& lvj1 = loccoords(6);
			const Vector3D& lvj3 = loccoords(7);

			Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A*lvi1;
			Vector3D gvi2 = A*lvi2;
			Vector3D v = GetBody3D(1).GetPos(lp1) - gp2;

			f(1) = v.X(); //translation
			f(2) = v.Y(); //translation
			f(3) = v.Z(); //translation

			f(4) = lvj3 * gvi1; //rotation
			f(5) = lvj3 * gvi2; //rotation
			f(6) = lvj1 * gvi2; //rotation
		}
		else
			//if (elements.Length()==2) 
			if (elements(2))	//$ DR 2012-10 elements.Length is now always 2
			{
				//body-body joint:
				const Vector3D& lp1 = loccoords(1);
				const Vector3D& lp2 = loccoords(2);
				const Vector3D& lvi1 = loccoords(3);
				const Vector3D& lvi2 = loccoords(4);
				const Vector3D& lvi3 = loccoords(5);
				const Vector3D& lvj1 = loccoords(6);
				const Vector3D& lvj3 = loccoords(7);

				Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
				Vector3D gvi1 = A1*lvi1;
				Vector3D gvi2 = A1*lvi2;
				Vector3D gvj1 = A2*lvj1;
				Vector3D gvj3 = A2*lvj3;
				Vector3D v = GetBody3D(1).GetPos(lp1) - GetBody3D(2).GetPos(lp2);

				f(1) = v.X(); //translation
				f(2) = v.Y(); //translation
				f(3) = v.Z(); //translation

				f(4) = gvj3 * gvi1; //rotation
				f(5) = gvj3 * gvi2; //rotation
				f(6) = gvj1 * gvi2; //rotation

			}
	}
	else if (MaxIndex()<=2)
	{
		//if (elements.Length()==1)
		if (elements(2)==0)	//$ DR 2012-10 elements.Length is now always 2
		{
			//ground joint, velocity:
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& gp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(3);
			const Vector3D& lvi2 = loccoords(4);
			const Vector3D& lvi3 = loccoords(5);
			const Vector3D& lvj1 = loccoords(6);
			const Vector3D& lvj3 = loccoords(7);

			Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D Ap= GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A*lvi1;
			Vector3D gvi2 = A*lvi2;
			Vector3D gvi1p = Ap*lvi1;
			Vector3D gvi2p = Ap*lvi2;
			//Vector3D v = GetBody3D(1).GetPos(lp1) - gp2;
			Vector3D vp = GetBody3D(1).GetVel(lp1);

			f(1) = vp.X(); //translation
			f(2) = vp.Y(); //translation
			f(3) = vp.Z(); //translation

			f(4) = lvj3 * gvi1p; //rotation
			f(5) = lvj3 * gvi2p; //rotation
			f(6) = lvj1 * gvi2p; //rotation
		}
		//else if (elements.Length()==2) 
		else if (elements(2))	//$ DR 2012-10 elements.Length is now always 2
		{
			//body-body joint, velocity:
			const Vector3D& lp1 = loccoords(1);
			const Vector3D& lp2 = loccoords(2);
			const Vector3D& lvi1 = loccoords(3);
			const Vector3D& lvi2 = loccoords(4);
			const Vector3D& lvi3 = loccoords(5);
			const Vector3D& lvj1 = loccoords(6);
			const Vector3D& lvj3 = loccoords(7);

			Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A1p= GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Matrix3D A2p= GetBody3D(2).GetRotMatrixP(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
			Vector3D gvi1 = A1*lvi1;
			Vector3D gvi2 = A1*lvi2;
			Vector3D gvj1 = A2*lvj1;
			Vector3D gvj3 = A2*lvj3;
			Vector3D gvi1p = A1p*lvi1;
			Vector3D gvi2p = A1p*lvi2;
			Vector3D gvj1p = A2p*lvj1;
			Vector3D gvj3p = A2p*lvj3;
			Vector3D vp = GetBody3D(1).GetVel(lp1) - GetBody3D(2).GetVel(lp2);

			f(1) = vp.X(); //translation
			f(2) = vp.Y(); //translation
			f(3) = vp.Z(); //translation

			f(4) = gvj3 * gvi1p + gvj3p * gvi1; //rotation
			f(5) = gvj3 * gvi2p + gvj3p * gvi2; //rotation
			f(6) = gvj1 * gvi2p + gvj1p * gvi2; //rotation
		}
	}
};

void RigidJointOLD::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	//if (elements.Length()==2)
	//{
	//	if (locelemind == 1)
	//		mbs->UO() << "A1:" << Matrix3D(0,0,1,0,1,0,-1,0,0)*GetBody3D(1).GetRotMatrix(loccoords(1)) << "\n";
	//	else
	//		mbs->UO() << "A2:" << Matrix3D(0,1,0,-1,0,0,0,0,1)*GetBody3D(2).GetRotMatrix(loccoords(2)) << "\n";
	//}

	hmat.SetSize(f.Length(),3);

	//if (elements.Length() == 1) //--> ignore locelemind, only first element!
	if (elements(2)==0)	//$ DR 2012-10 elements.Length is now always 2
	{
		//ground joint, velocity:
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& gp2 = loccoords(2);
		const Vector3D& lvi1 = loccoords(3);
		const Vector3D& lvi2 = loccoords(4);
		const Vector3D& lvi3 = loccoords(5);
		const Vector3D& lvj1 = loccoords(6);
		const Vector3D& lvj3 = loccoords(7);

		Matrix3D A = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A*lvi1;
		Vector3D gvi2 = A*lvi2;
		Vector3D v = GetBody3D(1).GetPos(lp1) - gp2;

		//for position constraints:
		GetBody3D(1).GetdPosdqT(lp1,dpdq);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
		}

		//for rotational constraints:
		GetBody3D(1).GetdRotvdqT(lvi1,lp1,dpdq); //Hi1
		Mult(dpdq,lvj3,hvec);
		hmat.SetColVec(hvec,1);

		GetBody3D(1).GetdRotvdqT(lvi2,lp1,dpdq); //Hi1
		Mult(dpdq,lvj3,hvec);
		hmat.SetColVec(hvec,2);
		Mult(dpdq,lvj1,hvec);
		hmat.SetColVec(hvec,3);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= (hmat(i,1)*XG(4)+hmat(i,2)*XG(5)+hmat(i,3)*XG(6));
		}
	}
	else
	{
		//body-body joint:
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lvi1 = loccoords(3);
		const Vector3D& lvi2 = loccoords(4);
		const Vector3D& lvi3 = loccoords(5);
		const Vector3D& lvj1 = loccoords(6);
		const Vector3D& lvj3 = loccoords(7);

		Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Matrix3D A2 = GetBody3D(2).GetRotMatrix(lp2); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A1*lvi1;
		Vector3D gvi2 = A1*lvi2;
		Vector3D gvj1 = A2*lvj1;
		Vector3D gvj3 = A2*lvj3;
		Vector3D v = GetBody3D(1).GetPos(lp1) - GetBody3D(2).GetPos(lp2);

		if (locelemind==1)
		{
			//body 1:
			//for position constraints:
			GetBody3D(1).GetdPosdqT(lp1,dpdq);

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= (dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
			}

			//for rotational constraints:
			GetBody3D(1).GetdRotvdqT(lvi1,lp1,dpdq); //Hi1
			Mult(dpdq,gvj3,hvec);
			hmat.SetColVec(hvec,1);

			GetBody3D(1).GetdRotvdqT(lvi2,lp1,dpdq); //Hi1
			Mult(dpdq,gvj3,hvec);
			hmat.SetColVec(hvec,2);
			Mult(dpdq,gvj1,hvec);
			hmat.SetColVec(hvec,3);

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= (hmat(i,1)*XG(4)+hmat(i,2)*XG(5)+hmat(i,3)*XG(6));
			}

		}
		else
		{
			//body 2
			//for position constraints:
			GetBody3D(2).GetdPosdqT(lp2,dpdq);

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= -(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
			}

			//for rotational constraints:
			GetBody3D(2).GetdRotvdqT(lvj3,lp2,dpdq); //Hi1
			Mult(dpdq,gvi1,hvec);
			hmat.SetColVec(hvec,1);
			Mult(dpdq,gvi2,hvec);
			hmat.SetColVec(hvec,2);

			GetBody3D(2).GetdRotvdqT(lvj1,lp2,dpdq); //Hi1
			Mult(dpdq,gvi2,hvec);
			hmat.SetColVec(hvec,3);

			for (int i=1; i <= f.Length(); i++)
			{
				f(i) -= (hmat(i,1)*XG(4)+hmat(i,2)*XG(5)+hmat(i,3)*XG(6));
			}
		}
	}
};

Vector3D RigidJointOLD::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void RigidJointOLD::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	if (GetDrawSizeScalar() == 0) return;

	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	const Vector3D& lvi1 = loccoords(3);
	const Vector3D& lvi2 = loccoords(4);
	const Vector3D& lvi3 = loccoords(5);
	const Vector3D& lvj1 = loccoords(6);
	const Vector3D& lvj3 = loccoords(7);

	//int tt = (int)draw_dim.Z(); //tiling of cylinders
	Vector3D pi;
	Vector3D pj;

	Vector3D roti1, roti2, roti3;
	Vector3D rotj1, rotj2, rotj3;

	roti1 = GetBody3D(1).GetRotMatrixD(lp1)*lvi1;
	roti1 *= 1.*GetDrawSizeScalar();
	roti2 = GetBody3D(1).GetRotMatrixD(lp1)*lvi2;
	roti2 *= 1.*GetDrawSizeScalar();
	roti3 = GetBody3D(1).GetRotMatrixD(lp1)*lvi3;
	roti3 *= 1.*GetDrawSizeScalar();

	//if (elements.Length()==1) 
	if (elements(2)==0)	//$ DR 2012-10 elements.Length is now always 2
	{
		pi = GetBody3D(1).GetPosD(lp1);
		pj = lp2;

		rotj1 = lvj1;
		rotj3 = lvj3;
		rotj2 = lvj3.Cross(lvj1);
	} else
	{
		pi = GetBody3D(1).GetPosD(lp1);
		pj = GetBody3D(2).GetPosD(lp2);

		rotj1 = GetBody3D(2).GetRotMatrixD(lp2)*lvj1;
		rotj2 = GetBody3D(2).GetRotMatrixD(lp2)*(lvj3.Cross(lvj1));
		rotj3 = GetBody3D(2).GetRotMatrixD(lp2)*lvj3;
	}
	roti1 *= 1.05;
	rotj1 *= 1.00*GetDrawSizeScalar();
	rotj2 *= 0.95*GetDrawSizeScalar();
	rotj3 *= 1.05*GetDrawSizeScalar();



	pi -= 0.5*(roti1+roti2+roti3);
	pj -= 0.5*(rotj1+rotj2+rotj3);

	mbs->DrawCube(pi, roti1, roti2, roti3);
	mbs->SetColor(colgrey2);
	mbs->DrawCube(pj, rotj1, rotj2, rotj3);
};
