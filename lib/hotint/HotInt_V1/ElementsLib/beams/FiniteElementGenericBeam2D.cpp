//#**************************************************************
//# filename:             FiniteElementGenericBeam2D.cpp
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
 


#include "element.h"
#include "material.h"
#include "Node.h"
#include "body2d.h"
#include "femathhelperfunctions.h"
#include "FiniteElementGenericBeam2D.h"

FiniteElementGenericBeam2D::~FiniteElementGenericBeam2D(void) {}

void FiniteElementGenericBeam2D::SetFiniteElementGenericBeam2D(
	int bodyind,
	const Vector& xc1, 
	const Vector& xc2, 
	int n1, 
	int n2, 
	int material_num,
	const Vector3D& asize,
	const Vector3D& color)
{
	TArray<int>node_list(2); node_list.Add(n1); node_list.Add(n2);		// generate list of node ids

	// SetFiniteElementGeneric(bodyind, node_list, material_num, color, 0);
	//AP 22-02-2013 uncommented the following line by YV, since the element is for nonlinear large strain
	//SetGeometricNonlinearityStatus(GNS_NonlinearSmallStrain);		// YV

	this->bodyind = bodyind;
	this->nodes = node_list;
	this->elementname = GetElementSpec();

	SetMaterialNum(material_num);
	col = color;

	//this->lx = size.X();
	//this->ly = size.Y();
	//this->lz = size.Z();
	this->size = asize;

	this->geometricNonlinearityStatus = GNS_NonlinearLargeStrain;

	q0 = (xc1.Append(xc2)).Append(Vector(SOS()));
	x_init = Vector(2*FlexDOF());
	x_init.SetAll(0);

	plasticip_x = 0;
	plasticip_y = 0;

	xgReferenceState.SetLen(FlexDOF());
	xgReferenceState.SetAll(0);
	xgInit.SetXGProvider(this, &xgReferenceState);
	xgCompute.SetXGProvider(this, &xgReferenceState);
	xgDraw.SetXGProvider(this, &xgReferenceState);
};

void FiniteElementGenericBeam2D::PreAssemble()
{
	IntegrationRule::IntegrationRuleSettings settings(
		GetElementType(), GetActualInterpolationOrder(), IntegrationRule::IVT_Stiffness, GetGeometricNonlinearityStatus(), 0
		);
	//$AP 2013-01-11 GetIntegrationRulesLibrary is now member of the element, 
	// I did not test this, only compile!!
	integrationRuleStiffness = /*mbs->*/GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);
	settings.integratedValueType = IntegrationRule::IVT_Mass;
	integrationRuleMass = /*mbs->*/GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);
	settings.integratedValueType = IntegrationRule::IVT_Load;
	integrationRuleLoad = /*mbs->*/GetIntegrationRulesLibrary()->GetIntegrationRule(settings, this);

	assert(
		integrationRuleStiffness != NULL &&
		integrationRuleMass != NULL &&
		integrationRuleLoad != NULL
		);

	// pre-compute shape functions in integration points
	// it turns out that this does not really improve things

	/*IntegrationPointsIterator ip(integrationRuleStiffness);
	for (ip; !ip.IsEnd(); ++ip)
	{
	if (ip.GetVectorData(1) == 0)
	{			
	Vector * sx = ip.GetVectorData(ip.AddVectorData());
	Vector * sxx = ip.GetVectorData(ip.AddVectorData());

	sx->SetLen(NS());
	sxx->SetLen(NS());

	for (int j = 1; j <= NS(); j++)
	{
	double x = ip.Point2D().X();
	(*sx)(j) = GetS0x(x, j);
	(*sxx)(j) = GetS0xx(x, j);
	}
	}
	}*/
};

void FiniteElementGenericBeam2D::CopyFrom(const Element &e)
{
	FiniteElementGeneric<Body2D>::CopyFrom(e);
	const FiniteElementGenericBeam2D& ce = (const FiniteElementGenericBeam2D&)e;
	this->q0 = ce.q0;
	//this->lx = ce.lx;
	//this->ly = ce.ly;
	//this->lz = ce.lz;
	this->plasticip_x = ce.plasticip_x;
	this->plasticip_y = ce.plasticip_y;
	xgReferenceState = ce.xgReferenceState;
};

Vector2D FiniteElementGenericBeam2D::GetRefPos2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(Vector2D(x0,0), i);
		p(1) += s * q0(2*i-1);
		p(2) += s * q0(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetRefPos2D(const Vector2D& p_loc) const
{
	Vector2D p0 = GetRefPos2D(p_loc.X());
	return p0 + GetLy()*0.5 * p_loc.Y() * GetRefInplaneUnitVector2D(p_loc.X());
};

Vector2D FiniteElementGenericBeam2D::GetDisplacement2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(Vector2D(x0,0), i);
		p(1) += s * XG(2*i-1);
		p(2) += s * XG(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetDisplacement2DD(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(Vector2D(x0,0), i);
		p(1) += s * XGD(2*i-1);
		p(2) += s * XGD(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetDisplacementx2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(Vector2D(x0,0), i) * 2./GetLx();
		p(1) += sx * XG(2*i-1);
		p(2) += sx * XG(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetDisplacementx2DD(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(Vector2D(x0,0), i) * 2./GetLx();
		p(1) += sx * XGD(2*i-1);
		p(2) += sx * XGD(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetDisplacementxx2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sxx = GetS0xx(Vector2D(x0,0), i) * 4./(GetLx()*GetLx());
		p(1) += sxx * XG(2*i-1);
		p(2) += sxx * XG(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetDisplacementxx2DD(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sxx = GetS0xx(Vector2D(x0,0), i) * 4./(GetLx()*GetLx());
		p(1) += sxx * XGD(2*i-1);
		p(2) += sxx * XGD(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetRefPosx2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(Vector2D(x0,0), i) * 2./GetLx();	// don't forget transformation to unit element, change to element jacobian in future
		p(1) += sx * q0(2*i-1);
		p(2) += sx * q0(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetRefPosxx2D(const double& x0) const
{
	Vector2D p(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sxx = GetS0xx(Vector2D(x0,0), i) * 4./(GetLx()*GetLx());	// don't forget transformation to unit element, change to element jacobian in future
		p(1) += sxx * q0(2*i-1);
		p(2) += sxx * q0(2*i  );
	}
	return p;
};

Vector2D FiniteElementGenericBeam2D::GetVel2D(const double& x0) const
{
	Vector2D v(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double s = GetS0(Vector2D(x0,0), i);
		v(1) += s * XGP(2*i-1);
		v(2) += s * XGP(2*i  );
	}
	return v;
}

Vector2D FiniteElementGenericBeam2D::GetVel2D(const Vector2D& p_loc) const
{
	return GetVel2D(p_loc.X()) + p_loc.Y() * GetInplaneUnitVectorP2D(p_loc.X());
}

Vector2D FiniteElementGenericBeam2D::GetVelx2D(const double& x0) const
{
	Vector2D v(0.,0.);
	for (int i = 1; i <= NS(); i++)
	{
		double sx = GetS0x(Vector2D(x0,0), i) * 2./GetLx();
		v(1) += sx * XGP(2*i-1);
		v(2) += sx * XGP(2*i  );
	}
	return v;
}

void FiniteElementGenericBeam2D::ComputeStrainMatrix2D(Vector2D ploc, Matrix2D& strain, const XGProvider& xg) const
{
	mbs->UO() << "Error in FiniteElementGenericBeam2D::ComputeStrainMatrix2D: not implemented for base class!\n";
	assert(0);
}

#pragma region draw routines
void FiniteElementGenericBeam2D::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement,
		FieldVariableDescriptor::FVCI_y);
	//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
	//	FieldVariableDescriptor::FVCI_y, true);
	//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_curvature);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_bending);
}

double FiniteElementGenericBeam2D::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D& p_loc, bool flagD)
{
	// displacements
	if(fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
		return fvd.GetComponent(GetDisplacement2D(p_loc.X()));

	double kappa, eps0;

	if (flagD)
	{
		kappa = GetKappaD(p_loc.X());
		eps0 = GetEpsAxialD(p_loc.X()); 
	}
	else
	{
		kappa = GetKappa(p_loc.X());
		eps0 = GetEpsAxial(p_loc.X()); 
	}
	double EI = GetBeamEIy();
	double EA = GetBeamEA();

	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_beam_axial_extension: return eps0;
	case FieldVariableDescriptor::FVT_beam_curvature: return kappa;
	case FieldVariableDescriptor::FVT_beam_force_axial: return EA * eps0;
	case FieldVariableDescriptor::FVT_beam_moment_bending: return -EI * kappa;		//not true for prescribed case!
	default: return FIELD_VARIABLE_NO_VALUE;
	}
};


void FiniteElementGenericBeam2D::DrawElement()
{
	mbs->SetColor(col);

	//// draw elements in front of origin, at z = 0.01 (AH: was 0.1 once)
	Vector3D offset(0, 0, 0.01);

	// length and magnified thickness of element
	double lx1 = GetLx(); double ly1 = GetLy()*GetMBS()->GetMagnifyYZ();

	// get draw mode from Visualization window
	//0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	int linemode = 1; 
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

	double def_scale = GetMBS()->GetDOption(105); //deformation scaling

	int drawlinemode = GetMBS()->GetIOption(117); // Draw Flat Elements, in this case, function is plotted normal to element axis

	if ( (1||colormode) && !drawlinemode) // plot coloured elements
	{
		lx1 /= GetLx();
		ly1 /= GetLy();

		double modeval = 0;
		int xgset = 0;

		double tilex = GetMBS()->GetIOption(137);
		double tiley = GetMBS()->GetIOption(138);

		TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
		TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
		double v=0;

		points.SetLen(0); vals.SetLen(0);
		Vector2D p0, p0_mag, vx, vy, vy_mag;
		int tileyn = (int)tiley;
		int tilexn = (int)tilex;

		p0 = Vector2D(-lx1,-1.);
		vx = Vector2D(2.*lx1/tilexn,0);
		vy = Vector2D(0,2./tileyn);

		p0_mag = Vector2D(-lx1,-ly1);
		vy_mag = Vector2D(0,2.*ly1/tileyn);

		for (double iy = 0; iy <= tileyn+1e-10; iy++)
		{
			for (double ix = 0; ix <= tilexn+1e-10; ix++)
			{
				Vector2D ploc = (p0+ix*vx+iy*vy);
				Vector2D ploc_mag = (p0_mag+ix*vx+iy*vy_mag);
				// Vector2D pg = GetPos2D0D(ploc_mag, def_scale); (AH)
				Vector2D pg = GetPos2DD(ploc_mag, def_scale);

				Vector3D p(ToP3D(pg)+offset);
				points.Add(p);
				if (colormode)
					v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
				vals.Add(v);
			}
		}
		mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
	}
	else if (drawlinemode) // plot function on center line over beam axis
	{
		double scalex = lx1/GetLx();
		double scaley = ly1/GetLy();

		Vector3D col_line(0.8,0,0);

		int tilex = GetMBS()->GetIOption(137);

		TArray<Vector3D> points(tilex+1), normals(tilex+1), drawpoints(tilex+1);
		TArray<double> vals(tilex+1);
		double v=0;

		points.SetLen(0); vals.SetLen(0); normals.SetLen(0); drawpoints.SetLen(0);
		Vector2D p0, vx, vy;

		p0 = Vector2D(-scalex,0);
		vx = Vector2D(2.*scalex/tilex,0);
		vy = Vector2D(0.,0.);

		int cnt = 0;
		for (double ix = 0.; ix <= tilex+1e-10; ix++)
		{
			// point on the middle line
			Vector2D ploc = (p0+ix*vx);
			Vector2D pg = GetPos2DD(ploc, def_scale);

			Vector3D p(ToP3D(pg)+offset);
			points.Add(p);

			// normal direction (numerical differetiation for the moment)
			Vector2D ploc_eps = ploc + 1e-5*vx;
			Vector2D pg_eps = GetPos2DD(ploc_eps, def_scale);

			Vector3D normal = Vector3D(-pg_eps.Y()+pg.Y(), pg_eps.X()-pg.X(), 0.);
			double scale = ly1/normal.Norm();
			normal *= scale;
			normals.Add(normal);

			// computation of value
			if (colormode)
				v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc, true);
			vals.Add(v);

			MBS* constmbs = const_cast<MBS*>(mbs);
			constmbs->UpdateFEMinMaxCol(v);
		}

		for (int cnt=1; cnt<=tilex+1; cnt++)
		{
			// scale normal according to minimum/maximum
			double scale = Maximum(fabs(GetMBS()->GetTImincol()), fabs(GetMBS()->GetTImaxcol()));
			if (scale==0) scale = 1.;
			// point on the function line
			drawpoints.Add(points(cnt)+(vals(cnt)/scale)*normals(cnt));

			// draw arrow for value of function, in normal direction to middle line, colored matching legend
			double diff = mbs->GetFEmaxcol()-mbs->GetFEmincol();
			Vector3D fecolor = mbs->FEColor( (vals(cnt)-mbs->GetFEmincol())/diff );
			mbs->MyDrawArrow(points(cnt), drawpoints(cnt), fecolor);
			if (cnt>1)
			{
				// red line for plotting function value
				mbs->MyDrawLine(drawpoints(cnt-1), drawpoints(cnt), ly1*1e2, col_line);
				// middle line of element
				mbs->MyDrawLine(points(cnt-1), points(cnt), ly1*5e1);
				// upper and lower line of element, unscaled
				mbs->MyDrawLine(points(cnt-1)+(0.5/scaley)*normals(cnt-1), points(cnt)+(0.5/scaley)*normals(cnt), ly1*5e1);
				mbs->MyDrawLine(points(cnt-1)-(0.5/scaley)*normals(cnt-1), points(cnt)-(0.5/scaley)*normals(cnt), ly1*5e1);
			}

		}

		TArray<Vector3D> nodepoints(4);
		nodepoints(1) = ToP3D(GetPos2DD(Vector2D(-scalex, 1.), def_scale)) + offset;
		nodepoints(2) = ToP3D(GetPos2DD(Vector2D(-scalex,-1.), def_scale)) + offset;
		nodepoints(3) = ToP3D(GetPos2DD(Vector2D( scalex,-1.), def_scale)) + offset;
		nodepoints(4) = ToP3D(GetPos2DD(Vector2D( scalex, 1.), def_scale)) + offset;
		//mbs->DrawPolygon(nodepoints, 1, ly*ly1*1e-1);
		// left and right line of element
		mbs->MyDrawLine(nodepoints(1), nodepoints(2), ly1*4e1);
		mbs->MyDrawLine(nodepoints(3), nodepoints(4), ly1*4e1);
	}
	/*else	// AH: imho, this won't ever be used again, keep it commented out though
	{
	double tiling = 24;

	Vector2D p1,p2,p3,p4;
	mbs->SetColor(col);

	for (double i = 0; i < tiling; i++)
	{
	double l1 = -lx1*0.5+lx1*i/tiling;
	double l2 = -lx1*0.5+lx1*(i+1)/tiling;
	if (i == 0)
	{
	p1 = GetPos2D0D(Vector2D(l1/(0.5*lx), 1.), def_scale);
	p2 = GetPos2D0D(Vector2D(l1/(0.5*lx),-1.), def_scale);
	}
	else
	{
	p1 = p4;
	p2 = p3;
	}
	p3 = GetPos2DD(Vector2D(l2,-ly*0.5));
	p4 = GetPos2DD(Vector2D(l2, ly*0.5));
	mbs->DrawQuad(ToP3D(p4),ToP3D(p3),ToP3D(p2),ToP3D(p1));
	}
	}*/
};
#pragma endregion

#pragma region plasticity
void FiniteElementGenericBeam2D:: DataToPlasticVariables(Vector& plastic_variables, int ip_x, int ip_y) const
{
	plastic_variables.SetLen(NPlasticParameters());
	if (plasticip_x == 0)
	{
		plastic_variables.SetAll(0.);
		return;
	}
	int offset = 0;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		plastic_variables(i) = XData(offset+plasticip_x*(ip_y-1)+ip_x);
		offset += NPlasticGridpoints();
	}
}
void FiniteElementGenericBeam2D:: DataToPlasticVariablesLastTimeStep(Vector & plastic_variables, int ip_x, int ip_y) const
{
	plastic_variables.SetLen(NPlasticParameters());
	if (plasticip_x == 0)
	{
		plastic_variables.SetAll(0.);
		return;
	}
	int idx = plasticip_x*(ip_y-1)+ip_x;
	int offset = 0;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		plastic_variables(i) = GetMBS()->GetLastDataVector()(LTGdata(offset+idx));
		offset += NPlasticGridpoints();
	}
}
void FiniteElementGenericBeam2D:: DataToPlasticVariablesD(Vector & plastic_variables, int ip_x, int ip_y) const
{
	plastic_variables.SetLen(NPlasticParameters());
	if (plasticip_x == 0)
	{
		plastic_variables.SetAll(0.);
		return;
	}
	int offset = 0;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		plastic_variables(i) = XDataD(offset+plasticip_x*(ip_y-1)+ip_x);
		offset += NPlasticGridpoints();
	}
}
// interpolate plastic variables to ploc
void FiniteElementGenericBeam2D:: DataToPlasticVariables(Vector & plastic_variables, Vector2D& ploc)
{
	plastic_variables.SetLen(NPlasticParameters());
	if (plasticip_x == 0)
	{
		plastic_variables.SetAll(0.);
		return;
	}
	int offset = 0;
	Matrix plasticstrain_quantity;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		plasticstrain_quantity.LinkWith(plasticip_y, plasticip_x, &XData(offset+1));
		plastic_variables(i) = plasticstrain_quantity.Interpolate(ploc);
		offset += NPlasticGridpoints();
	}
}
void FiniteElementGenericBeam2D:: DataToPlasticVariablesD(Vector & plastic_variables, const Vector2D ploc)
{
	plastic_variables.SetLen(NPlasticParameters());
	if (plasticip_x == 0)
	{
		plastic_variables.SetAll(0.);
		return;
	}
	int offset = 0;
	Matrix plasticstrain_quantity;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		plasticstrain_quantity.LinkWith(plasticip_y, plasticip_x, &XDataD(offset+1));
		plastic_variables(i) = plasticstrain_quantity.Interpolate(ploc);
		offset += NPlasticGridpoints();
	}
}

// interpolate plastic strain only, not hardening parameter and yield function
void FiniteElementGenericBeam2D:: DataToPlasticStrain(Vector & plasticstrain, Vector2D& ploc)
{
	// -2 due to missing hardening parameter and yield function
	plasticstrain.SetLen(NPlasticParameters()-2);
	if (plasticip_x == 0)
	{
		plasticstrain.SetAll(0.);
		return;
	}
	int offset = 0;
	Matrix plasticstrain_quantity;
	for (int i=1; i<=NPlasticParameters()-2; i++)
	{
		plasticstrain_quantity.LinkWith(plasticip_y, plasticip_x, &XData(offset+1));
		plasticstrain(i) = plasticstrain_quantity.Interpolate(ploc);
		offset += NPlasticGridpoints();
	}
}


// set plastic variables to internal XData
void FiniteElementGenericBeam2D:: PlasticVariablesToData(const Vector & plastic_variables, int ip_x, int ip_y)
{
	int offset = 0;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		XData(offset+plasticip_x*(ip_y-1)+ip_x) = plastic_variables(i);
		offset += NPlasticGridpoints();
	}
}
void FiniteElementGenericBeam2D:: PlasticVariablesToDataD(const Vector & plastic_variables, int ip_x, int ip_y)
{
	int offset = 0;
	for (int i=1; i<=NPlasticParameters(); i++)
	{
		XDataD(offset+plasticip_x*(ip_y-1)+ip_x) = plastic_variables(i);
		offset += NPlasticGridpoints();
	}
}

double FiniteElementGenericBeam2D::PostNewtonStep(double t)
{
	// no PostNewtonStep for elastic materials
	if (!GetMaterial().IsInelasticMaterial())
	{
		return 0;
	}
	// Plasticity...
	double error = 0;
	// loop over all grid points
	double hx = 2./(plasticip_x-1.);
	double hy = 2./(plasticip_y-1.);
	XGProviderCached<FE2DmaxDOF> xg(xgCompute, SOS());

	for (int i=1; i<=plasticip_y; i++)
	{
		for (int j=1; j<=plasticip_x; j++)
		{
			// grid points are ordered as in matrix
			Vector2D ploc(-1.+(j-1)*hx, 1.-(i-1)*hy);
			// strain matrix
			Matrix2D strain;
			ComputeStrainMatrix2D(ploc, strain, xg);
			// plasticVariables in struct;
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars(NPlasticParameters());
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars_store(NPlasticParameters());
			DataToPlasticVariables(plasticvars, j, i);
			plasticvars_store = plasticvars;
			ConstVector<MAX_PLASTIC_QUANTITIES> plasticvars_last(NPlasticParameters());
			DataToPlasticVariablesLastTimeStep(plasticvars_last, j, i);

			plasticvars = plasticvars_last;

			// recompute strain due to new idea...
			double sigma22_ref = Em()*(strain(2,2)-GetPlasticThicknessStrain(ploc(1)));
			double old_strain22 = strain(2,2);
			strain(2,2) = sigma22_ref/Em()+plasticvars(2);
			if (GetMaterial().IsPlaneStrain())
			{
				GetMaterial().ComputeInelasticVariablesFromStrainPlaneStrain2D(plasticvars, strain);
			}
			else // plane stress
			{
				GetMaterial().ComputeInelasticVariablesFromStrainPlaneStress2D(plasticvars, strain);
			}

			// add to error the difference between the stored plasticvars_store and the new plasticvars (only strain components
			// scale by length of element and by number of integration points
			double delta_epsplast_laststep = 0;
			double delta_epsplast_thisstep = 0;
			for (int ii=1; ii<=NPlasticParameters()-2; ii++) 
			{
				delta_epsplast_laststep += fabs(plasticvars(ii)-plasticvars_last(ii)); 
				delta_epsplast_thisstep += fabs(plasticvars(ii)-plasticvars_store(ii)); 
			}
			// if delta_epsplast > 0, then yield function has to be zero (compatibility condition)
			if (plasticvars(NPlasticParameters()) < - GetMBS()->DiscontinuousAccuracy() && delta_epsplast_laststep > GetMBS()->DiscontinuousAccuracy()/Em())
			{
				error += size(1)*delta_epsplast_laststep	/(plasticip_x*plasticip_y);
				plasticvars = plasticvars_last;
			}
			else // otherwise, add error quantity and use plastic variables returned in vector
			{
				error += size(1)*delta_epsplast_thisstep	/(plasticip_x*plasticip_y);
			}

			// map vector back to plasticvars and then to XData
			PlasticVariablesToData(plasticvars, j, i);


		}
	}

	return error;
}


// link XData corresponding to the plastic quantity quantitynr with Matrix mat
// 1 .. epsp_xx
// 2 .. epsp_yy
// 3 .. epsp_xy
// 4 .. hardeningparameter
// 5 .. yield function
void FiniteElementGenericBeam2D::GetPlasticQuantityMatrix(int quantitynr, const Matrix& mat, int flagD) const
{
	if (!GetMaterial().IsInelasticMaterial()) return;
	int firstel = (quantitynr-1)*NPlasticGridpoints()+1;
		FiniteElementGenericBeam2D& constthis = const_cast<FiniteElementGenericBeam2D&>(*this);
		Matrix& constmat = const_cast<Matrix&>(mat);
	if (!flagD)
	{
		constmat.LinkWith(plasticip_y, plasticip_x, &(constthis.XData(firstel)) );
	}
	else if (flagD)
	{
		constmat.LinkWith(plasticip_y, plasticip_x, &(constthis.XDataD(firstel)) );
	}
}
void FiniteElementGenericBeam2D::GetPlasticQuantityMatrix(int quantitynr, Matrix& mat, int flagD)
{
	if (!GetMaterial().IsInelasticMaterial()) return;
	int firstel = (quantitynr-1)*NPlasticGridpoints()+1;
	if (!flagD)
	{
		mat.LinkWith(plasticip_y, plasticip_x, &(XData(firstel)) );
	}
	else if (flagD)
	{
		mat.LinkWith(plasticip_y, plasticip_x, &(XDataD(firstel)) );
	}
}
// routines for retrieving integral values of plastic strain
// plastic axial strain = 1/h int_{-h/2}^h/2 epsp_xx dy
double FiniteElementGenericBeam2D:: GetPlasticAxialStrain(double x0, int flagD)
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	Matrix epsp_xx;
	GetPlasticQuantityMatrix(1, epsp_xx, flagD);
	double plasticaxstrain = IntegralMean(x0, epsp_xx);
	return plasticaxstrain;
}
// plastic material curvature = 12/h^3 int_{-h/2}^h/2 y epsp_xx dy
double FiniteElementGenericBeam2D:: GetPlasticKappa(double x0, int flagD)
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	double kappa = 0;
	double hy = size(2)/(plasticip_y-1.);
	// Matrix contains plastic strain matrix
	Matrix epsp_xx;
	GetPlasticQuantityMatrix(1, epsp_xx, flagD);
	// Vector plasticstrains contains plastic strains in gridpoints over height interpolated at x-coord x0
	ConstVector<MAX_PLAST_GRIDP_Y> epsxx_vec(plasticip_y);
	epsp_xx.InterpolateX(x0, epsxx_vec);

	// Integrate plasticstrain*y over height - exact formula in each plastic cell
	for (int cell=1; cell<plasticip_y; cell++)
	{
		double y1 = hy*(cell-1)-0.5*size(2);
		double y2 = hy*(cell)-0.5*size(2);
		double eps1 = epsxx_vec(plasticip_y-(cell-1));
		double eps2 = epsxx_vec(plasticip_y-cell);
		kappa += 2./Cub(size(2)) * ( eps2*(2*y2*y2-y1*y1) + (eps1-eps2)*y1*y2 + eps1*(y2*y2-2*y1*y1) );
		//kappa2 += 12./(Cub(size(2)))*0.5*(y2*eps2+y1*eps1)*hy;
	}
	return -kappa;
}
// plastic shear strain = 1/h int_{-h/2}^h/2 epsp_xy dy
double FiniteElementGenericBeam2D::  GetPlasticShearStrain(double x0, int flagD)
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	Matrix epsp_xy;
	GetPlasticQuantityMatrix(3, epsp_xy, flagD);
	double plasticshearstrain = IntegralMean(x0, epsp_xy);
	return plasticshearstrain;
}
// plastic thickness strain = 1/h int_{-h/2}^h/2 epsp_yy dy
double FiniteElementGenericBeam2D::  GetPlasticThicknessStrain(double x0, int flagD)
{
	if (!GetMaterial().IsInelasticMaterial()) return 0.;
	Matrix epsp_yy;
	GetPlasticQuantityMatrix(2, epsp_yy, flagD);
	double plasticthicknessstrain = IntegralMean(x0, epsp_yy);
	return plasticthicknessstrain;
}

// help routine to integrate some plastic strain over height
// return 1/h int_{-h/2}^h/2 epsp dy
double FiniteElementGenericBeam2D::  IntegralMean(double x0, Matrix& plasticstrain_mat)
{
	int ny = plasticip_y;
	// Vector plasticstrain_vec contains plastic strains in gridpoints over height interpolated at x-coord x0
	ConstVector<MAX_PLAST_GRIDP_Y> plasticstrain_vec(ny);
	plasticstrain_mat.InterpolateX(x0, plasticstrain_vec);

	// Integrate plasticstrain over height - trapezoidal rule
	double eps_ax = 0;
	for (int cell=1; cell<ny; cell++)
	{
		eps_ax += 0.5*(plasticstrain_vec(ny-(cell-1))+plasticstrain_vec(ny-cell));
	}
	return eps_ax/(ny-1.);
}

// help routine to evaluate some plastic strain on center line
// return plastic strain evaluated at y=0
double FiniteElementGenericBeam2D::  CenterLineValue(double x0, Matrix& plasticstrain_mat)
{
	int ny = plasticip_y;
	// Vector plasticstrain_vec contains plastic strains in gridpoints over height interpolated at x-coord x0
	ConstVector<MAX_PLAST_GRIDP_Y> plasticstrain_vec(ny);
	plasticstrain_mat.InterpolateX(x0, plasticstrain_vec);

	int odd_gridpoints = plasticip_y%2;
	// odd number of gridpoints -> use center value
	if (odd_gridpoints)
	{
		return plasticstrain_vec(plasticip_y/2+1);
	}
	else // even number of gridpoints -> use mean value of the two center values
	{
		return 0.5*(plasticstrain_vec(plasticip_y/2) + plasticstrain_vec(plasticip_y/2+1));
	}

}

#pragma endregion