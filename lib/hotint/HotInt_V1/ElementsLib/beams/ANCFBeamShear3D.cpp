//#**************************************************************
//# filename:             ANCFBeamShear3D.cpp
//#
//# author:               Karin & Astrid
//#
//# generated:						Nov. 2010
//# description:          3D ANCF beam element with shear deformation
//# comments:							formulation is based on the following paper: 
//#												"A 3D Shear Deformable Finite Element Based on the Absolute Nodal Coordinate Formulation"
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
#include "femathhelperfunctions.h"
#include "material.h"
#include "ANCFBeamShear3D.h"
#include "Node.h"
//#include "graphicsconstants.h"
#include "elementdataaccess.h"
//#include "solversettings_auto.h"




#pragma region ANCFBeamShear3DGeneric

double ANCFBeamShear3DGeneric::CalculateElementLength() const
{
	//to obtain stress free reference configuration: same numerical integration as in axial deformation part of EvalF2()
	int int_order_axial = (BeamFormulation == 4) ? 3 : order_axial;   // must be set equally to integration order in axial deformation part of strain energy (EvalF2()) for each BeamFormulation, respectively!!!
	Vector xT, wT;
	GetIntegrationRule(xT, wT, int_order_axial); 

	double val = 0.;
	for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
	{
		double x = xT(i1);		//x in unit configuration [-1 , 1]
		double dxidx = GetLx()*0.5;	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])
		double xi = x*dxidx;	//xi in scaled configuration [-lx/2 , lx/2]
		Vector3D p0(xi, 0., 0.);

		// integrate r,xi along xi in [-lx/2, +lx/2]
		Vector sfx(NS(), 1); GetShapesx(sfx, p0);
		Matrix Sfx(Dim(), NS()*Dim());
		for (int i=1; i<=NS(); i++)
		{
			for (int j=1; j<=Dim(); j++)
			{
				Sfx(j, (i-1)*Dim() + j) = sfx(i);
			}
		}
		Vector posx = Sfx*q0;
		double norm_posx = posx.GetNorm();
		val += norm_posx*wT(i1)*dxidx;
	}

	return val;
}


// default set function
void ANCFBeamShear3DGeneric::SetANCFBeamShear3DGeneric(int materialnumi, const Vector3D& coli)
{
	//set reference configuration
	int idx = 1;
	for(int i=1; i <= nnodes; i++)   
	{
		ANCFNodeS2S3_3D& n((ANCFNodeS2S3_3D&)GetNode(i));
		q0(idx++) = n.Pos().X();
		q0(idx++) = n.Pos().Y();
		q0(idx++) = n.Pos().Z();
		q0(idx++) = n.GetRefSlope2().X();
		q0(idx++) = n.GetRefSlope2().Y();
		q0(idx++) = n.GetRefSlope2().Z();
		q0(idx++) = n.GetRefSlope3().X();
		q0(idx++) = n.GetRefSlope3().Y();
		q0(idx++) = n.GetRefSlope3().Z();
	}

	materialnum = materialnumi;

	double lx = CalculateElementLength();
	double ly, lz;
	Beam3DProperties& bp((Beam3DProperties&)GetMaterial());
	if (bp.IsMaterialOfBeamWithRectangularCrossSection())
	{
		ly = bp.GetBeamThicknessY();
		lz = bp.GetBeamThicknessZ();
	}
	else
	{
		mbs->UO(UO_LVL_err) << "WARNING: Only rectangular cross section implemented for ANCFBeamShear3D at the moment!\n";
	}
	size.Set(lx,ly,lz);

	xg.SetLen(SOS());
	xg.SetAll(0.);
	
	col = coli;
	BuildDSMatrices();

	penalty.SetLen(9);
	penalty.SetAll(1.);

	x_init.SetLen(2*SOS());
	x_init.SetAll(0.);

	SetElementName(GetElementSpec());
}

// deprecated set function!!!
void ANCFBeamShear3DGeneric::SetANCFBeamShear3DGeneric(int materialnumi, const Vector3D& si, const Vector3D& coli)
{
	size = si;

	xg.SetLen(SOS());
	xg.SetAll(0.);
	
	//lx = size.X(); ly = size.Y(); lz = size.Z(); //lx is parameter to identify reference configuration
	materialnum = materialnumi;
	//mass = lx*ly*lz*GetMaterial().Density();
	col = coli;
	BuildDSMatrices();

	penalty.SetLen(9);
	penalty.SetAll(1.);

	x_init.SetLen(2*SOS());
	x_init.SetAll(0.);

	SetElementName(GetElementSpec());
}

void ANCFBeamShear3DGeneric::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Body3D::GetElementData(edc);
}

//user has changed data ==> write new data and parameters to element
int ANCFBeamShear3DGeneric::SetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	int rv = Body3D::SetElementData(edc);
	return rv;
}

void ANCFBeamShear3DGeneric::BuildDSMatrices()
{
};


void ANCFBeamShear3DGeneric::Gradu(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const
{
	gradu.SetSize(3,3);
	gradu.SetAll(0);

	int dim = Dim();
	int l;
	for (int j = 1; j <= dim; j++) 
	{
		for (int i = 1; i <= NS(); i++)
		{
			l = (i-1)*dim+j;
			gradu(j,1) += GetSFx(i,ploc)*u(l);
			gradu(j,2) += GetSFy(i,ploc)*u(l);
			gradu(j,3) += GetSFz(i,ploc)*u(l);
		}
	}
}

//----------------
//Velocity
//----------------

//flagD = 0 für Berechnungen, flagD = 1 for Visualization

Vector3D ANCFBeamShear3DGeneric::GetVel3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++) //i=1,2,3
	{
		for (int j = 1; j <= NS(); j++)
		{	
			if(flagD==0)
				p(i) += GetSF(j,p_loc)*XGP((j-1)*Dim()+i);  //XGP=actual solution vector

			else
				p(i) += GetSF(j,p_loc)*XGPD((j-1)*Dim()+i);
		}
	}
	return p;
};

//----------------
//Displacement
//----------------
//only for visualization
Vector3D ANCFBeamShear3DGeneric::GetDisplacementD(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSF(j,p_loc)*(XGD((j-1)*Dim()+i));
		}
	}
	return p;
};

Vector3D ANCFBeamShear3DGeneric::GetDisplacement(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NS(); j++)
		{
			p(i) += GetSF(j,p_loc)*(XG((j-1)*Dim()+i));
		}
	}
	return p;
};

//Vector3D ANCFBeamShear3DGeneric::GetDisplacement3D(const Vector3D& p_loc, const Vector& xg) const
//{
//	Vector3D p(0.,0.,0.);
//	for (int i = 1; i <= Dim(); i++)
//	{
//		for (int j = 1; j <= NS(); j++)
//		{
//			p(i) += GetSF(j,p_loc)*(xg((j-1)*Dim()+i));
//		}
//	}
//	return p;
//};
//

void ANCFBeamShear3DGeneric::GetIntDuDq(Matrix& dudq)
{
	//integrate du/dq over beam volume   => (crossection area) * (integration over axis)
	//u(q) = r(q) - r0  ==>  du/dq = dr/dq

	dudq.SetSize(NS()*Dim(), Dim());
	dudq.SetAll(0);

	ConstVector<4> w1, x1;//max. 4 integration points -> order 7
	GetIntegrationRule(x1, w1, 3);
	double fact = GetLx()*0.5*GetLy()*GetLz();									//crosssection area times determinant for line integral

	for (int i1=1; i1<=x1.GetLen(); i1++)
	{
		double scaled_weight = w1(i1) * fact;			//scale weight
		double xi = 0.5*GetLx()*x1(i1);
		Vector3D ploc(xi, 0., 0.);

		for (int j=1; j<=NS(); j++)
		{
			double sf_value = GetSF(j, ploc);			//value returned by i'th (positional) shape function at point xi
			double block_j = Dim()*(j-1);						//addresses (j-1)st block in dudq
			for (int k=1; k<=Dim(); k++)
			{
				dudq(k+block_j, k) += scaled_weight*sf_value;
			}
		}
	}
}

//----------------
//GetPos
//p_loc in [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [-Lz/2,Lz/2]
//----------------

Vector3D ANCFBeamShear3DGeneric::GetPosD(const Vector3D& p_loc) const 
{ 
	return GetPos3D(p_loc, 1); 
}

//Vector3D ANCFBeamShear3DGeneric::GetPosD(const Vector3D& p_loc) const
//// p_loc in [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [-Lz/2,Lz/2]
//{
//	//mbs->UO() << "posy=" << GetPosy3D(p_loc, 1) << "\n";
//	//mbs->UO() << "posz=" << GetPosz3D(p_loc, 1) << "\n";
//	//mbs->UO() << "ploc=" << p_loc << "\n";
//	Vector3D p0(p_loc.X(),0.,0.);
//
//	return GetPos3D(p_loc, 1) + p_loc.Y()*GetPosy3D(p_loc, 1) + p_loc.Z()*GetPosz3D(p_loc, 1);
//}

//GetPos2D berechnet Vektor Positionsvektor r im deformierten Element
//GetPos2D computes r in deformed element
//xg = generalized coordinates (= our unknowns q)
Vector3D ANCFBeamShear3DGeneric::GetPos3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0.);
	//TMStartTimer(24);
	Vector sf(NS(), 1);
	GetShapes(sf, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	//TMStopTimer(24);
	return p;
};

//GetPos2D von oben, aber falls im Aufruf als zweite Eingabe kein Vektor, sondern ein Integer steht,
//wird diese Fkt für Grafik aufgerufen; globale Variable XG wird für Rechnung verwendet anstatt Eingabe xg
//xg=XG wird einmal gesetzt, weshalb GetPos3D mit xg schneller ist, wohingegen GetPos3D mit XG jedes Mal aufs globale XG zugreift und das dauert lange

//same GetPos2D function as above, but with the second parameter flag D the function for the visualization is activated
Vector3D ANCFBeamShear3DGeneric::GetPos3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapes(sf, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sf(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{ 
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sf(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
};

//Ableitungen von GetPos entsprechen unseren Richtungsvektoren
//GetPosx2D (Ableitung in Richtung der Balkenachse: dr/dx=r,x)
//GetPosx2D is the derivative with respect to the direction of the beam centerline
Vector3D ANCFBeamShear3DGeneric::GetPosx3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0.);
	Vector sfx(NS(), 1);
	GetShapesx(sfx, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sfx(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosx3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfx(NS(), 1);
	GetShapesx(sfx, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfx(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{ 
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfx(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
}


void ANCFBeamShear3DGeneric::GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) 
{
	dpdx = GetPosx3D(ploc);
};


//GetPosy2D (Ableitung in Richtung des Querschnitts: dr/dy=r,y)
//GetPosy2D is the derivative with respect to the direction of the cross section
Vector3D ANCFBeamShear3DGeneric::GetPosy3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0.);
	Vector sfy(NS(), 1);
	GetShapesy(sfy, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sfy(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosy3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfy(NS(), 1);
	GetShapesy(sfy, p_loc);
	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfy(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfy(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosyP3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfy(NS(), 1);
	GetShapesy(sfy, p_loc);
	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfy(j)*XGP(offset + i);
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfy(j)*XGPD(offset + i);
			}
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosz3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0);
	Vector sfz(NS(), 1);
	GetShapesz(sfz, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sfz(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosz3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfz(NS(), 1);
	GetShapesz(sfz, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfz(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfz(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPoszP3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfz(NS(), 1);
	GetShapesz(sfz, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfz(j)*XGP(offset + i);
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfz(j)*XGPD(offset + i);
			}
		}
	}
	return p;
}


//GetPosxy2D: d2r/dydx=r,xy
//is used in DeltaThetax
Vector3D ANCFBeamShear3DGeneric::GetPosxy3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0.);
	Vector sfxy(NS(), 1);
	GetShapesxy(sfxy, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sfxy(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosxy3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfxy(NS(), 1);
	GetShapesxy(sfxy, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfxy(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfxy(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosxz3D(const Vector3D& p_loc, const Vector& xg) const
{
	Vector3D p(0.,0.,0.);
	Vector sfxz(NS(), 1);
	GetShapesxz(sfxz, p_loc);

	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sfxz(j)*(xg(offset + i) + q0(offset + i)); //r,y=S,y*u+S,y*q_0		
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetPosxz3D(const Vector3D& p_loc, int flagD) const
{
	Vector3D p(0.,0.,0.);
	Vector sfxz(NS(), 1);
	GetShapesxz(sfxz, p_loc);

	if (!flagD)
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfxz(j)*(XG(offset + i) + q0(offset + i));
			}
		}
	}
	else
	{
		for (int j = 1; j <= NS(); j++)
		{
			int offset = (j-1)*Dim();
			for (int i = 1; i <= Dim(); i++)
			{
				p(i) += sfxz(j)*(XGD(offset + i) + q0(offset + i));
			}
		}
	}
	return p;
}

//computes position vector in the reference element: r_0 
//p_loc in [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [-Lz/2,Lz/2]
Vector3D ANCFBeamShear3DGeneric::GetInitPos3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapes(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
};

Vector3D ANCFBeamShear3DGeneric::GetInitPosx3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapesx(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
};

Vector3D ANCFBeamShear3DGeneric::GetInitPosy3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapesy(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetInitPosz3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapesz(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetInitPosxy3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapesxy(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
}

Vector3D ANCFBeamShear3DGeneric::GetInitPosxz3D(const Vector3D& p_loc) const
{
	Vector3D p(0.,0.,0.);
	Vector sf(NS(), 1);
	GetShapesxz(sf, p_loc);
	for (int j = 1; j <= NS(); j++)
	{
		int offset = (j-1)*Dim();
		for (int i = 1; i <= Dim(); i++)
		{
			p(i) += sf(j)*q0(offset + i);
		}
	}
	return p;
}

//for visualization
Vector3D ANCFBeamShear3DGeneric::GetPos3D0D(const Vector3D& p_loc) const 
{
	Vector3D plocscaled;
	plocscaled(1)=p_loc(1)*GetLx()/2.;
	plocscaled(2)=p_loc(2)*GetLy()/2.;
	plocscaled(3)=p_loc(3)*GetLz()/2.;
	return GetPos3D(plocscaled,1); 
};

Vector3D ANCFBeamShear3DGeneric::GetPos3D0D(const Vector3D& p_loc, double defscale) const 
{
	Vector3D plocscaled;
	plocscaled(1)=p_loc(1)*GetLx()/2.;
	plocscaled(2)=p_loc(2)*GetLy()/2.;
	plocscaled(3)=p_loc(3)*GetLz()/2.;
	return GetInitPos3D(plocscaled) + defscale*GetDisplacementD(plocscaled); 

	GetNodePos(1);
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Vector3D ANCFBeamShear3DGeneric::GetNodeLocPos(int i) const //returns position of i-th node
{
	switch(i)
	{
	case 1:	return Vector3D(-GetLx()*.5, 0, 0);
	case 2:	return Vector3D(+GetLx()*.5, 0, 0);
	case 3:	assert(nnodes >= 3); return Vector3D(0.);
	}
	assert(0);
	return Vector3D(0.);
}


Vector3D ANCFBeamShear3DGeneric::GetNodePos(int i) const //returns position of i-th node
{
	switch(i)
	{
	case 1:	return GetPos(Vector3D(-GetLx()*.5, 0, 0));
	case 2:	return GetPos(Vector3D(+GetLx()*.5, 0, 0));
	case 3:	assert(nnodes >= 3); return GetPos(Vector3D(0.));
	}
	assert(0);
	return Vector3D(0.);
}

Vector3D ANCFBeamShear3DGeneric::GetNodePosD(int i) const //returns position of i-th node (draw mode)
{
	switch(i)
	{
	case 1:	return GetPosD(Vector3D(-GetLx()*.5, 0, 0));
	case 2:	return GetPosD(Vector3D(+GetLx()*.5, 0, 0));
	case 3:	assert(nnodes >= 3); return GetPosD(Vector3D(0.));
	}
	assert(0);
	return Vector3D(0.);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector3D ANCFBeamShear3DGeneric::GetDOFDirD(int idof) const //returns direction of action of i-th DOF
{
	if ((idof >= 1 && idof <= 3) ||
		(idof >= 10 && idof <= 12) ||
		(idof >= 19 && idof <= 21))
	{
		Vector3D v(0.);
		int component = (idof-1)%3+1;
		v(component) = 1;
		return v;
	}
	//for other components, the direction can not be easily identified and is therefore set to zero (default)

	return Vector3D(0.,0.,0.);
}

Vector3D ANCFBeamShear3DGeneric::GetDOFPosD(int idof) const //returns position of i-th DOF
{
	if (idof <= 9) //left node: 1,...,9
	{
		return GetPos3D(Vector3D(-0.5*GetLx(),0.,0.),1);//1 für flagD
	}
	else if (idof <= 18) //right node: 10,...,18
	{
		return GetPos3D(Vector3D(0.5*GetLx(),0.,0.),1);
	}
	else if (idof <= 27) //middle node: 19,...,27
	{
		return GetPos3D(Vector3D(0.,0.,0.),1);//middle node: (0,0)
	}
	else
	{
		return Vector3D(0, 0, 0);
	}
}


Matrix3D ANCFBeamShear3DGeneric::GetRotMatrix(const Vector3D& p0, int on_reference_element, int flagD) const //p0(1) in [-lx/2,lx/2]
{
	//ROTATION TENSOR: 
	//A = [ e1 | e2 | e3 ]

	//e1 = e1^{\bar}/|e1^{\bar}| with e1^{\bar} = drdy x drdz
	//e3 = e3^{\bar}/|e3^{\bar}| with e3^{\bar} = drdz
	//e2 = e2^{\bar}/|e2^{\bar}| with e2^{\bar} = drdz x (drdy x drdz) = e3^{\bar} x e1^{\bar}
	
	//Attention: Getei(i) returns ei^{\bar} 

  Vector3D drdy, drdz;

	if (on_reference_element)
	{
		drdy = GetInitPosy3D(p0);
		drdz = GetInitPosz3D(p0);
	}
	else
	{
		drdy = GetPosy3D(p0, flagD);
		drdz = GetPosz3D(p0, flagD);
	}
	
	TArray<Vector3D> e(Dim());
	for (int i=1; i<=Dim(); i++)
	{
		e(i) = Getei(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		e(i).Normalize();
	}

	Matrix3D A(e(1).X(),e(2).X(),e(3).X(),
					   e(1).Y(),e(2).Y(),e(3).Y(),
					   e(1).Z(),e(2).Z(),e(3).Z());
	return A;
}

Matrix3D ANCFBeamShear3DGeneric::GetRotMatrixP(const Vector3D& p0, int flagD) const
{
	//Time derivative of rotation tensor: 
	//AP = [ e1P | e2P | e3P ]

	//e1 = e1^{\bar}/|e1^{\bar}| with eP1^{\bar} = drdyP x drdz + drdy x drdzP
	//e3 = e3^{\bar}/|e3^{\bar}| with eP3^{\bar} = drdzP
	//e2 = e2^{\bar}/|e2^{\bar}| with eP2^{\bar} = drdzP x (drdy x drdz) + drdz x (drdyP x drdz + drdy x drdzP) = eP3^{\bar} x e1^{\bar} + e3^{\bar} x eP1^{\bar}
	
	//Attention: Getei(i) returns ei^{\bar}, and GetePi(i) returns ePi^{\bar} 

  Vector3D drdy, drdz, drPdy, drPdz;

	drdy = GetPosy3D(p0, flagD);
	drdz = GetPosz3D(p0, flagD);
	drPdy = GetPosyP3D(p0, flagD);
	drPdz = GetPoszP3D(p0, flagD);
	
	Matrix3D AP;
	for (int i=1; i<=3; i++)
	{
		Vector3D ei = Getei(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		double ei_norm = ei.Norm();
		Vector3D eiP = GeteiP(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z(), drPdy.X(), drPdy.Y(), drPdy.Z(), drPdz.X(), drPdz.Y(), drPdz.Z());
		double ei_eiP = ei*eiP;
		eiP *= (1./ei_norm);
		eiP += (-ei_eiP/pow(ei_norm,3))*ei;	
		for (int j=1; j<=3; j++)
		{
			AP(j,i) = eiP(j);
		}
	}

	return AP;
}

Vector3D ANCFBeamShear3DGeneric::GetTwistAndCurvature(const Vector3D& ploc, int on_reference_element, int flagD) const 
// ploc.X() in [ -Lx/2 , Lx/2 ]
// calculate vector of twist and curvature k = 1/2 sum_i e_i x e_i'
// or, for the reference element (see flag) k_0 = 1/2 sum_i e0_i x e0_i'
// NOTE: both vectors are returned with their components w.r.t. the local frame basis (e_i)!
{
	if(IsStraightBeamInReferenceConfiguration() && on_reference_element)  // for speedup
	{
		return Vector3D();   //std constructor initializes with zero
	}

	Vector3D drdy, drdz, ddrdxdy, ddrdxdz;
	if (on_reference_element)
	{
		drdy = GetInitPosy3D(ploc);
		drdz = GetInitPosz3D(ploc);
		ddrdxdy = GetInitPosxy3D(ploc);
		ddrdxdz = GetInitPosxz3D(ploc);
	}
	else
	{
		drdy = GetPosy3D(ploc, flagD);
		drdz = GetPosz3D(ploc, flagD);
		ddrdxdy = GetPosxy3D(ploc, flagD);
		ddrdxdz = GetPosxz3D(ploc, flagD);
	}
	
	Vector3D ei;     // returns ei^{\bar}, ei = ei^{\bar} / |ei^{\bar}|
                                 //e1 = e1^{\bar}/|e1^{\bar}| with e1^{\bar} = drdy x drdz
	                               //e3 = e3^{\bar}/|e3^{\bar}| with e3^{\bar} = drdz
	                               //e2 = e2^{\bar}/|e2^{\bar}| with e2^{\bar} = drdz x (drdy x drdz) = e3^{\bar} x e1^{\bar}
	Vector3D deidxi; // returns d ei^{\bar} / d xi

	Vector3D k;
	for (int i=1; i<=Dim(); i++)
	{
		ei = Getei(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		deidxi = Getdeidxi(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z(), ddrdxdy.X(), ddrdxdy.Y(), ddrdxdy.Z(), ddrdxdz.X(), ddrdxdz.Y(), ddrdxdz.Z());

		k += ( 0.5 / ei.Norm2() ) * ei.Cross(deidxi);  //k = 1/2 sum_i e_i x e_i' -> simplified written using e_i^{\bar}-Vectors (see the paper)
	}

	Matrix3D AT = GetRotMatrix(ploc, on_reference_element, flagD).GetTp();   //in order to represent the vector k in the local frame basis

	return AT * k;
}

Vector3D ANCFBeamShear3DGeneric::GetAngularVel(const Vector3D& p_loc) const
{
	Matrix retmat;
	Vector xp(SOS());
	for (int j=1; j<=SOS(); j++)
	{
		xp(j)=XGP(j);
	}
	Matrix3D A = GetRotMatrix(p_loc);
	Matrix3D omega_skew;
	for (int i=1; i<=Dim(); i++)
	{
		Vector3D ei(A(1,i),A(2,i),A(3,i));
		GetdRotvdqT(ei, p_loc, retmat, 1);
		Vector wi(3);
		Mult(xp, retmat, wi);
		for (int j=1; j<=Dim(); j++)
		{
			omega_skew(i,j) = wi(j);
		}
	}
	//Matrix3D omega_skew = GetRotMatrixP(p_loc)*GetRotMatrix(p_loc).GetTp();
	//omega_skew * A = RotMatrixP

	return Vector3D(omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2));
}

void ANCFBeamShear3DGeneric::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& drotvdqT, int use_transposed_rotmatrix) const
{
	// rot = (rot11 rot12 rot13
	//        rot21 rot22 rot23
	//        rot31 rot32 rot33)
	// Basis vectors are columns of rotation tensor
	// rotT.v = (rot11*v1 + rot21*v2 + rot31*v3, rot12*v1 + rot22*v2 + rot32*v3, rot13*v1 + rot23*v2 + rot33*v3)
	// p = (r,y , r,z) = ( r1dy , r2dy , r3dy , r1dz , r2dz , r3dz )
	// drotTv/dq = drotTv/dp . dp/dq

	double xi = ploc.X();		// ploc.X() in [ -Lx/2 , Lx/2 ]

	//rotation tensor only depends on r,y and r,z
	Vector3D drdy = GetPosy3D(xi);
	Vector3D drdz = GetPosz3D(xi);

	Vector3D e1(Getei(1, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
	Vector3D e2(Getei(2, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
	Vector3D e3(Getei(3, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
	TArray<Vector3D> de1dp(GetdeidDiffVar(1, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
	TArray<Vector3D> de2dp(GetdeidDiffVar(2, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
	TArray<Vector3D> de3dp(GetdeidDiffVar(3, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));

// calculate derivatives of normalized ei
	double e1norm = e1.Norm();
	double e2norm = e2.Norm();
	double e3norm = e3.Norm();
	for (int j=1; j<=de1dp.Length(); j++)
	{
		de1dp(j) = de1dp(j)*(1./e1norm) - e1*((e1*de1dp(j))/Cub(e1norm));
		de2dp(j) = de2dp(j)*(1./e2norm) - e2*((e2*de2dp(j))/Cub(e2norm));
		de3dp(j) = de3dp(j)*(1./e3norm) - e3*((e3*de3dp(j))/Cub(e3norm));
	}

  Matrix drotvdp(Dim(), 6);//3x6
	if (use_transposed_rotmatrix)    //use transposed A: A^T*v
	{
		for (int h=1; h<=6; h++)//h...DiffVar
		{
			drotvdp(1,h) = de1dp(h)*vloc;
			drotvdp(2,h) = de2dp(h)*vloc;
			drotvdp(3,h) = de3dp(h)*vloc;
		}
	}
	else //dont use transposed A: A*v
	{
		for (int h=1; h<=6; h++) //h...DiffVar
		{
			drotvdp(1,h) = de1dp(h).X()*vloc.X() + de2dp(h).X()*vloc.Y() + de3dp(h).X()*vloc.Z();
			drotvdp(2,h) = de1dp(h).Y()*vloc.X() + de2dp(h).Y()*vloc.Y() + de3dp(h).Y()*vloc.Z();
			drotvdp(3,h) = de1dp(h).Z()*vloc.X() + de2dp(h).Z()*vloc.Y() + de3dp(h).Z()*vloc.Z();
		}
	}
	//h = 1,2,3 -> SFy
	//h = 4,5,6 -> SFz
	drotvdqT.SetSize(SOS(),Dim());
	drotvdqT.SetAll(0.);
	Vector sfy(NS(), 1); GetShapesy(sfy, ploc);
	Vector sfz(NS(), 1); GetShapesz(sfz, ploc);

	for (int k=1; k<=Dim(); k++)
	{
		for (int i=1; i<=3; i++)
			for (int j=1; j<=NS(); j++)
				drotvdqT((j-1)*Dim()+i, k) = drotvdp(k,i)*sfy(j) + drotvdp(k,Dim()+i)*sfz(j);
	}
}

void ANCFBeamShear3DGeneric::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	//ploc in [-lx/2,lx/2] x [-ly/2,ly/2] x [-lz/2,lz/2]
	////////////////////////////////////////////////////////////
	//KN&PG's approach   \cite{hodges06} 
	//let e_i denote the i-th base vector of the actual system (i-th column of the rotation matrix)
	//\delta \theta = \sum_{i=1}^3  e_i (\delta e_j . e_k)    where j = i%3+1, k = j%3+1
	//(\grad_q \theta^T)_{l,m} = (e_i)_m (e_k)_n (partial (e_j)_n/ \partial q_l\)  
	//note: for this method the call of routine GetdRotvdqT has to sum up (e_j)_n v_n instead of (e_j)_n v_j, i.e. use_transposed_rotmatrix = 1
	////////////////////////////////////////////////////////////

	d.SetSize(SOS(),Dim());
	d.FillWithZeros();

	double xi = ploc.X();		// ploc.X() in [ -Lx/2 , Lx/2 ]

	//rotation tensor only depends on r,y and r,z
	Vector3D drdy = GetPosy3D(xi);
	Vector3D drdz = GetPosz3D(xi);

	//-----------------------------------------
	// new version
	//-----------------------------------------
	Vector3D ei;
	Vector3D ej;
	Vector3D ek;
	TArray<Vector3D> dejdp(6);
	TArray<Vector3D> dekdp(6);
	TArray<Vector3D> deidp(6);
	ConstVector<27> helpy;
	helpy.SetAll(0.);


	bool define_theta_by_a_transposed = false;
	bool derivatives_wrt_rotated_vars = false;

	for (int i=1; i<=Dim(); i++)	//i=1,2,3
	{
		int j = i%3 + 1;		//i=1->j=2->k=3, i=2->j=3->k=1, i=3->j=1->k=2
		int k = j%3 + 1;

		// e_i = ebar_i	/ |ebar_i|
		ei = Getei(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		ej = Getei(j, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		ek = Getei(k, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());

		if (define_theta_by_a_transposed)
		{
			double memo;
			memo = ei(j);
			ei(j) = ej(i);
			ej(i) = memo;
			memo = ej(k);
			ej(k) = ek(j);
			ek(j) = memo;
			memo = ek(i);
			ek(i) = ei(k);
			ei(k) = memo;
		}

		ei.Normalize();
		ek.Normalize();
		double ej_norm = ej.Norm();
		double ej_norm3 = Cub(ej_norm);
		ej.Normalize();

		// GetdeidDiffVar(j,...) returns debar_j_dp
		deidp = GetdeidDiffVar(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		dejdp = GetdeidDiffVar(j, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		dekdp = GetdeidDiffVar(k, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());

		// calculate all n=2*3=6 partial derivatives of ej (=normalized ebar_j) by using the derivatives of ebar_j
		for (int n=1; n<=dejdp.Length(); n++)
		{
			if (define_theta_by_a_transposed)
			{
				dejdp(n)(i) = deidp(n)(j);
				dejdp(n)(k) = dekdp(n)(j);
			}

			//dejdp(n) = dejdp(n)*(1./ej_norm) - ej*((ej*dejdp(n))/ej_norm3);
			dejdp(n) *= 1./ej_norm;
			dejdp(n) -= ej*(ej*dejdp(n));
		}

		if (!derivatives_wrt_rotated_vars)
		{
			//assembling of a help-vector saving at the (m*Dim()+n)th position the value
			//(ek^T * dejdpy(n))*SFy(m) + (ek^T * dejdpz(n))*SFz(m)
			for (int n=1; n<=Dim(); n++)
			{
				double fac_y = ek*dejdp(n);
				double fac_z = ek*dejdp(n+Dim());
				for (int m=1; m<=NS(); m++)
				{
					helpy((m-1)*Dim() + n) = fac_y*GetSFy(m, ploc) + fac_z*GetSFz(m, ploc);
				}
			}
		}
		else
		{
			//assembling of a help-vector saving at the (m*Dim()+n)th position the value
			//(ek^T * dejdpy(n))*A^T*SFy(m) + (ek^T * dejdpz(n))*A^T*SFz(m)

			Vector3D vec_y, vec_z;
			for (int n=1; n<=Dim(); n++)
			{
				vec_y(n) = ek*dejdp(n);
				vec_z(n) = ek*dejdp(n+Dim());
			}

			Matrix3D A = GetRotMatrix(ploc);

			vec_y = vec_y*A.GetTp();
			vec_z = vec_z*A.GetTp();


			for (int n=1; n<=Dim(); n++)
			{
				for (int m=1; m<=NS(); m++)
				{
					helpy((m-1)*Dim() + n) = vec_y(n)*GetSFy(m, ploc) + vec_z(n)*GetSFz(m, ploc);
				}
			}
		}

		//dimension of output matrix d is SOS() times Dim()   (transpose of dAdq)
		for (int m=1; m<=SOS(); m++)
		{
			for (int n=1; n<=Dim(); n++)
			{
				d(m,n) += helpy(m)*ei(n);
			}
		}
	}


	////-----------------------------------------	
	////old version:
	////-----------------------------------------
	////
	//Matrix3D A = GetRotMatrix(ploc);
	////Matrix3D I = A.GetTp()*A;
	////mbs->UO() << I.DoubleContract(I)-3. << "\n";

	//Matrix drotvdqT(SOS(),Dim());
	//for (int i=1; i<=Dim(); i++)//i=1,2,3
	//{
	//	int j = i%3 + 1;//i=1->j=2->k=3, i=2->j=3->k=1, i=3->j=1->k=2
	//	int k = j%3 + 1;
	//	Vector3D ei(A(1,i), A(2,i), A(3,i)); //i-th column of A
	//	GetdRotvdqT(Vector3D(A(1,k), A(2,k), A(3,k)), ploc, drotvdqT, 1);  //k-th column of A (=v) * dRotdq -> dRotvdq
	//	for (int l=1; l<=SOS(); l++)//SOS=Dim*NS=3*9=27
	//		for (int m=1; m<=Dim(); m++)//3
	//			d(l,m) += drotvdqT(l,j)*ei(m);
	//}
}


//----------------
//Massmatrix M
//----------------
void ANCFBeamShear3DGeneric::EvalM(Matrix& m, double t)
{
	if (massmatrix.Getcols() == SOS())
	{
		m = massmatrix;
		return;
	}
	else
	{
		int dim = Dim(); 
		int ns = NS();
		Vector SV;
		SV.SetLen(ns);
		Matrix HL(ns*dim,dim);
		Vector x1,x2,x3,w1,w2,w3; //integration points and according weights

		//GetIntegrationRule => Gauß points on interval [-1,+1]
		GetIntegrationRule(x1,w1,7); //optimal: 6x3x3, 6x1x1 geht auch; (GetIntegrationRule-Order == 6 complies with == 7)
		GetIntegrationRule(x2,w2,3);
		GetIntegrationRule(x3,w3,3);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				for (int i3=1; i3<=x3.GetLen(); i3++)
				{
					double x = x1(i1)*0.5*GetLx();  //x1 in [-1,+1], scaled rect. reference element: \xi
					double y = x2(i2)*0.5*GetLy();
					double z = x3(i3)*0.5*GetLz();
					Vector3D ploc(x,y,z);
					Vector3D rx0 = GetInitPosx3D(ploc);
					double rxn = rx0.Norm(); //curvature in reference element
					double jacdet = 0.5*GetLx()*0.5*GetLy()*0.5*GetLz() * rxn; //element transformation

					GetShapes(SV,ploc);
					for (int i=0; i<ns; i++)
					{
						for (int j=1; j<=dim; j++)
						{
							HL(i*dim+j,j)=SV(i+1);
						}
					}
					m += fabs (jacdet) * Rho() * w1(i1)*w2(i2)*w3(i3) * (HL*HL.GetTp());
				}
			}
		}
	}  

};

void ANCFBeamShear3DGeneric::GetJacobi(Matrix3D& jac, const Vector3D& ploc, const Vector& xg0) const
{
	jac.SetSize(3,3);
	jac.SetAll(0.);
	for (int i = 1; i <= Dim(); i++)
	{ 
		for (int k=1; k <= NS(); k++)
		{ 
			jac(i,1) += GetSFx(k, ploc)*xg0((k-1)*Dim()+i);
			jac(i,2) += GetSFy(k, ploc)*xg0((k-1)*Dim()+i);
			jac(i,3) += GetSFz(k, ploc)*xg0((k-1)*Dim()+i);
		}
	}
}

void ANCFBeamShear3DGeneric::EvalF2(Vector& f, double t)
{
	Body3D::EvalF2(f,t);
	//TMStartTimer(22);  //CPU timing
	Vector f_elastic_forces(SOS());	 //element residual - components are set to 0

	if (InContinuumMechanicsFormulation())
	{
		EvalF2CM(f_elastic_forces, t);
	}
	else
	{
		EvalF2SM(f_elastic_forces, t);
	}

	f -= f_elastic_forces;  //f=f-fadd, (Ku=f -> Residual=f-Ku, Ku=fadd)
	//TMStopTimer(22);
}

//residual of variational formulation in structural mechanics formulation
void ANCFBeamShear3DGeneric::EvalF2SM(Vector& fadd, double t)
{
	if (!InContinuumMechanicsFormulation())
		//Reissner-Simo-VuQuoc Formulation
	{

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//integration rule must be same as in GetFieldVariableValue(...) !!!! JG2013-09
		int order_x_bend = 3;						//3		 
		int order_x_thickness_eps = 4;  //4
		int order_x_axial_eps = 3;      //3
		int order_x_shear = 3;					//3
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		double beamGJ = GetMaterial().BeamGJkx();
		double beamEIy = GetMaterial().BeamEIy();
		double beamEIz = GetMaterial().BeamEIz();
		double beamEA = GetMaterial().BeamEA();
		double beamGAky = GetMaterial().BeamGAky();
		double beamGAkz = GetMaterial().BeamGAkz();
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		double beamGAkyz = BeamGAkyz(beamGAky, beamGAkz);   // ACHTUNG: Literaturstudie notwendig - dieser Faktor wird für den letzten Term bzgl. Dickendeformation verwendet (GAkyz*2E_yz*2deltaE_yz).
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		//UO() << "beamEIy=" << beamEIy << "\n";
		//UO() << "beamEIz=" << beamEIz << "\n";

		Vector faddadd(SOS());
		Vector3D gk(0.);								//container for Gamma, Kappa, and Strain E=(E_yy,E_zz,E_yz) respectively
		Matrix delta_gk(SOS(),Dim());		//container for DeltaGamma, and DeltaKappa, and DeltaStrain respectively
		Vector wT,xT;										//integration points and weights
		GetIntegrationRule(xT,wT,order_x_bend);
		
		// add kappa terms (torsion + 2*bending)
		//TMStartTimer(23);  //CPU timing %52%
		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*GetLx()*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*GetLx();	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact1 = penalty(1)*beamGJ*wT(i1)*det;	// beamGJ  = torsional stiffness (including torsional correction factor, depending on the shape of the cross section)
			double fact2 = penalty(2)*beamEIy*wT(i1)*det;	// beamEIy = bending stiffness (in e2=y direction))
			double fact3 = penalty(3)*beamEIz*wT(i1)*det;	// beamEIz = bending stiffness (in e3=z direction)

			GetDeltaKappa(xi, delta_gk, gk);
			//mbs->UO() << "int-point" << i1 << "\n";
			//if(i1 == 2)
			//{
			//	//mbs->UO() << "test\n";
			//	ofstream ofs; ofs.open("..\\..\\output\\FWF\\kappa_bei_GC12_Biegung_um_y_linelem_ip3.txt");
			//	ofs << "gk = " << gk << "\n" << "delta_gk = " << delta_gk << "\n";
			//	ofs.close();
			//}
			gk.Scale(1./fact1, 1./fact2, 1./fact3);			
			
			Mult(delta_gk, gk, faddadd);
			fadd += faddadd;
		}

		// add gamma terms (axial strains only)
		GetIntegrationRule(xT,wT,order_x_axial_eps); 
		//GetIntegrationRuleLobatto(xT,wT,2); 

		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*GetLx()*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*GetLx();	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact1 = penalty(4)*beamEA*wT(i1)*det;		  // beamEA  = axial stiffness 
			//integration of axial term only!!!
			double fact2 = 0.; //Gamma2 and Gamma3 is integrated seperately
			double fact3 = 0.;		

			//GetGamma(xi, gk);
			GetDeltaGamma(xi, delta_gk, gk);
			gk.Scale(1./fact1, 1./fact2, 1./fact3);	
			//UO() << "gk=" << gk << "\n";
			Mult(delta_gk, gk, faddadd);

			fadd += faddadd;
		}

		// add gamma terms (shear terms only)
		GetIntegrationRule(xT,wT,order_x_shear);
		//GetIntegrationRuleLobatto(xT,wT,order_x_shear);
		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*GetLx()*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*GetLx();	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact1 = 0.;		//integration of shear terms only!!!
			double fact2 = penalty(5)*beamGAky*wT(i1)*det;		
			double fact3 = penalty(6)*beamGAkz*wT(i1)*det;		
			
			GetDeltaGamma(xi, delta_gk, gk);
			gk.Scale(1./fact1, 1./fact2, 1./fact3);	
			//UO() << "gk=" << gk << "\n";
			Mult(delta_gk, gk, faddadd);
			fadd += faddadd;
		}
		//TMStopTimer(24);  //CPU timing

		// add thickness terms (E_yy E_zz 2*E_yz)
		//TMStartTimer(25);  //CPU timing %7.5%
		//GetIntegrationRule(xT,wT,order_x_thickness_eps);
		GetIntegrationRuleLobatto(xT,wT,order_x_thickness_eps);
		for (int i1=1; i1<=xT.GetLen(); i1++)  //i1-loop over all integration points
		{
			double x = xT(i1);		//x in unit configuration [-1 , 1]
			double xi = x*GetLx()*0.5;	//xi in scaled configuration [-lx/2 , lx/2]
			double det = 0.5*GetLx();	//determinant of the jacobian (from element transformation [-1 , 1] -> [-lx/2 , lx/2])

			double fact1 = penalty(7)*beamEA*wT(i1)*det;
			double fact2 = penalty(8)*beamEA*wT(i1)*det;
			double fact3 = penalty(9)*beamGAkyz*wT(i1)*det;		// ACHTUNG: beamGAkyz: Literaturstudie noch notwendig!!!!!!!
			
			GetDeltaThicknessStrain(xi, delta_gk, gk);
			gk.Scale(1./fact1, 1./fact2, 1./fact3);	

			Mult(delta_gk, gk, faddadd);
			fadd += faddadd;
		}
		//TMStopTimer(25);  //CPU timing
	}
}

//residual of variational formulation in continuum mechanics formulation
void ANCFBeamShear3DGeneric::EvalF2CM(Vector& fadd, double t)
{
	int dim = Dim();//=3
	int ns = NS();//=9
	int sos = SOS();//9*3

	ConstVector<27> temp;//KN-todo: 27 ersetzen durch ConstInt MaxDOF ...
	temp.SetLen(SOS());
	temp.SetAll(0);

	double A = GetMaterial().BeamA();					//cross section area
	double EA = GetMaterial().BeamEA();				//Em*area
	double GAky = GetMaterial().BeamGAky();		
	double GAkz = GetMaterial().BeamGAkz();		
	double G = Em()/(2.+2.*Nu());

	SetComputeCoordinates();//xg(i) = XG(i);

	Matrix3D strain, piola1, F;
	strain.SetAll(0);
	piola1.SetAll(0);

	int poissoncorrection = 0; //1==reduced integration of poisson part
	if (BeamFormulation == 2) poissoncorrection = 1;

	static Vector x1,x2,x3,w1,w2,w3;//KN-todo: ersetzen durch ConstVector

	//Vector3D ks(1.,GAky/(G*A),GAkz/(G*A)); //vector of shear correction factors
	//Matrix3D R = GetRotMatrix(Vector3D(0.,0.,0.), 1).GetTp(); //WARNING:for curved beams, this must be evaluated at every point of integration!
	//Vector3D ksR = R*ks;
	//ksR.X() = fabs(ksR.X());
	//ksR.Y() = fabs(ksR.Y());
	//ksR.Z() = fabs(ksR.Z());
	//UO() << "rot=" << R << "\n";
	//UO() << "ks=" << ks << "\n";
	//UO() << "ksR=" << ksR << "\n";

	for (int kk=1; kk <= 1+poissoncorrection; kk++)
	{
		ConstMatrix<36> Dm;
		Dm.SetSize(6,6);
		Dm.SetAll(0.);
		if (poissoncorrection)
		{
			if (BeamFormulation == 2)//locking compensated
			{
				if (kk == 1)//D^0 (equ. (25))
				{
					//poisson ratio zero:
					Dm(1,1) = Em();
					Dm(2,2) = Em();//reduced thickness stiffness 0.5*
					Dm(3,3) = Em();//reduced thickness stiffness

					//Dm(4,4) = G*5./6.;//1.;//*ksR.X();//statt GAks/A; //gamma_yz
					//Dm(5,5) = G*5./6.;//*ksR.Y();           //gamma_xz
					//Dm(6,6) = G*5./6.;//*ksR.Z();//event. weglassen -> konstant über Querschnitt //gamma_xy
					//UO() << "CMF (locking-free)" << "\n";
					Dm(4,4) = G;//statt GAks/A; //gamma_yz
					Dm(5,5) = GAky/A;           //gamma_xz
					Dm(6,6) = GAkz/A;//event. weglassen -> konstant über Querschnitt //gamma_xy

					//x_i ... integration points, w_i ... integration weights
					GetIntegrationRule(x1,w1,order_axial);//for dx
					GetIntegrationRule(x2,w2,order_crosssectional);//for dy
					GetIntegrationRule(x3,w3,order_crosssectional);//for dz
				}
				else if (kk == 2)//D^v (equ. (26))
				{
					//$ KN 2011-10-28: help1, help2 mit Paper verglichen -> passt
					double help1 = (Em()* (1-Nu()) / ((1+Nu())*(1-2.*Nu()))) - Em();
					//help1 = (2.*Em()*Nu()*Nu())/((1+Nu())*(1-2.*Nu()))) ... equ.(26)
					Dm(1,1) = help1;
					Dm(2,2) = help1;
					Dm(3,3) = help1;
					
					double help2 = Em()*Nu()/((1+Nu())*(1-2.*Nu()));
					Dm(1,2) = help2;
					Dm(2,1) = help2;
					Dm(1,3) = help2;
					Dm(3,1) = help2;
					Dm(2,3) = help2;
					Dm(3,2) = help2;

					GetIntegrationRule(x1,w1,order_axial);//x_i ... integration points, w_i ... integration weights
					GetIntegrationRule(x2,w2,1);//only 1 int. point along cross section//GaußIntegration
					GetIntegrationRule(x3,w3,1);
				}
			}
		}
		else if(BeamFormulation == 1)//cont. mech. approach which shows locking phenomenon
		{
			//x_i ... integration points, w_i ... integration weights
			GetIntegrationRule(x1,w1,order_axial);//for dx
			GetIntegrationRule(x2,w2,order_crosssectional);//for dy//Gauß!!!
			GetIntegrationRule(x3,w3,order_crosssectional);//for dz

			double help3 = Em()* (1-Nu()) / ((1+Nu())*(1-2*Nu()));
			double help4 = Em()*Nu() / ((1+Nu())*(1-2*Nu()));

			Dm(1,1)=help3;
			Dm(1,2)=help4;
			Dm(1,3)=help4;

			Dm(2,1)=help4;
			Dm(2,2)=help3;
			Dm(2,3)=help4;

			Dm(3,1)=help4;
			Dm(3,2)=help4;
			Dm(3,3)=help3;

			Dm(4,4) = G;//statt GAks/A; //gamma_yz
			Dm(5,5) = GAky/A;           //gamma_xz
			Dm(6,6) = GAkz/A;//event. weglassen -> konstant über Querschnitt //gamma_xy
			/*UO() << "false!\n";*/
			//Dm(4,4)=G*5./6.;
			//Dm(5,5)=G*5./6.;
			//Dm(6,6)=G*5./6.;
			//Dm(4,4)=10*Em();//G*5./6.;
			//Dm(5,5)=10*Em();//G*5./6.;
			//Dm(6,6)=10*Em();//G*5./6.;
			//UO() << "D=" << Dm << "\n";
		}
		else if(BeamFormulation == 3)// only D^0 (equ. (25))
		{
			Dm(1,1) = Em();
			Dm(2,2) = Em();
			Dm(3,3) = Em();
			//Dm(4,4) = G*ksR.X();//statt GAks/A; //gamma_yz
			//Dm(5,5) = G*ksR.Y();           //gamma_xz
			//Dm(6,6) = G*ksR.Z();//event. weglassen -> konstant über Querschnitt //gamma_xy
			/*UO() << "false!\n";*/
			Dm(4,4) = G*5./6.;//statt GAks/A;
			Dm(5,5) = G*5./6.;
			Dm(6,6) = G*5./6.;

			//x_i ... integration points, w_i ... integration weights
			GetIntegrationRule(x1,w1,order_axial);//for dx
			GetIntegrationRule(x2,w2,order_crosssectional);//for dy
			GetIntegrationRule(x3,w3,order_crosssectional);//for dz
		}
		//mbs->UO()<<"Dm="<<Dm<<"\n";

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				for ( int i3=1; i3<=x3.GetLen(); i3++)
				{
					int i,j,k;
					Vector3D p(x1(i1)*0.5*GetLx(), x2(i2)*0.5*GetLy(), x3(i3)*0.5*GetLz());

					Matrix3D jac;
					GetJacobi(jac,p,q0);//p(1) in [-lx/2,lx/2]
					double jacdet = jac.Det();
					//mbs->UO()<<"jacdet = "<<jacdet<<"\n";
					//mbs->UO()<<"0.5*lx*0.5*ly*0.5*lz = "<<0.5*lx*0.5*ly*0.5*lz<<"\n";

					Matrix3D jacinv;
					jac.GetInverse(jacinv);
					jacinv = jacinv.GetTp();
					//mbs->UO() << "jacinv = " << jacinv << "\n";

					static Matrix grad0, grad;//am Einheitselement
					grad0.SetSize(Dim(),NS());
					for(i=1; i<=NS(); i++)
					{
						grad0(1,i)=GetSFx(i,p);//3x9-matrix
						grad0(2,i)=GetSFy(i,p);
						grad0(3,i)=GetSFz(i,p);
					}

					Mult(jacinv, grad0, grad);//grad0*jacinv
					F.SetAll(0);//muss in jedem Integrationsschritt zurück auf 0 gesetzt werden!!!!
					for (int j = 1; j <= dim; j++) 
					{
						for (int i = 1; i <= NS(); i++)
						{ 
							int l = (i-1)*dim+j;
							F(j,1) += grad(1,i)*xg(l);//3x9-matrix*9-vector
							F(j,2) += grad(2,i)*xg(l);
							F(j,3) += grad(3,i)*xg(l);
						}
					}
					//F = F.GetTp();//FTranspose?
					if(calclinear == 0)
					{
						F(1,1) += 1;//Positionsgrad (nicht Grad(u), sondern Grad(r))
						F(2,2) += 1;
						F(3,3) += 1;
						//E = 1/2*(F^T F-I)
						strain = 0.5*(F.GetTp()*F);//eps=0.5*(F.GetTp()*F-Id)
						strain(1,1) -= 0.5;
						strain(2,2) -= 0.5;
						strain(3,3) -= 0.5;
					}
					else if(calclinear == 1)
					{
						strain = 0.5*(F + F.GetTp());
					}

					ConstVector<6> s3(strain(1,1), strain(2,2), strain(3,3), 2.*strain(2,3), 2.*strain(1,3), 2.*strain(1,2));
					ConstVector<6> stress;
					stress = Dm * s3;//Dm=6x6-Matrix, //S = D : E

					//generate 2. Piola-Kirchhoff-tensor, 2.PK = symm. tensor 
					piola1(1,1)=stress(1);
					piola1(2,2)=stress(2);
					piola1(3,3)=stress(3);

					piola1(2,3)=stress(4);
					piola1(3,2)=stress(4);
					piola1(1,3)=stress(5);
					piola1(3,1)=stress(5);
					piola1(1,2)=stress(6);
					piola1(2,1)=stress(6);

					if(calclinear == 0)
					{
						piola1 = F * piola1;//P = F * S, 1.PK = F * 2.PK
					}

					//grad = grad.GetTp();//gradTranspose?
					for (int j=1; j <= dim; j++)
					{
						for (int i = 0; i < ns; i++)
						{
							temp(dim*i+j) = grad(1,i+1)*piola1(j,1) + grad(2,i+1)*piola1(j,2) + grad(3,i+1)*piola1(j,3);
						}
					}
					//integrate(P : dF/dq) 
					fadd.MultAdd(fabs(jacdet) * 0.5*GetLx()*0.5*GetLy()*0.5*GetLz() * w1(i1)*w2(i2)*w3(i3),temp);
		
				}
			}
		}	
	}
};

//----------------------------------------------------------------------------
//Gamma & Kappa & ThicknessStrain  +  Derivatives of those with respect to dof
//----------------------------------------------------------------------------
void ANCFBeamShear3DGeneric::GetGamma(const double& xi, Vector3D& gamma, int flagD)		// -Lx/2 <= xi <= Lx/2
{
	//gamma1 = t1*drdx/|e1| - 1
	//gamma2 = t2*drdxi
	//gamma3 = t3*drdxi
	//
	//gamma = (gamma1, gamma2, gamma3)
	Vector3D p0(xi, 0., 0.);     //xi in [-Lx/2, Lx/2]

	Vector3D drdx, drdy, drdz, e1, e2, e3;
	drdx = GetPosx3D(p0, flagD);
	drdy = GetPosy3D(p0, flagD);
	drdz = GetPosz3D(p0, flagD);
	e1 = Getei(1, 
		drdy.X(), drdy.Y(), drdy.Z(), 
		drdz.X(), drdz.Y(), drdz.Z());
	e2 = Getei(2,
		drdy.X(), drdy.Y(), drdy.Z(),
		drdz.X(), drdz.Y(), drdz.Z());
	e3 = Getei(3,
		drdy.X(), drdy.Y(), drdy.Z(),
		drdz.X(), drdz.Y(), drdz.Z());
	e1.Normalize();
	e2.Normalize();
	e3.Normalize();

	gamma.Set(
		e1*drdx - 1.,
		e2*drdx,
		e3*drdx);
}

void ANCFBeamShear3DGeneric::GetDeltaGamma(const double& xi, Matrix& dgammadq, Vector3D& gamma)		// -Lx/2 <= xi <= Lx/2
{
	//dgamma1dq = ( ... )   see paper
	//dgamma2dq = ( r,y.S,x + r,x.S,y - (r,x.r,y)/(r,y.r,y).r,y.S,y ) / (norm(r,y))
	//dgamma3dq = ( r,z.S,x + r,x.S,z - (r,x.r,z)/(r,z.r,z).r,z.S,z ) / (norm(r,z))
	//
	//dgammadq = (dgamma1dq, dgamma2dq, dgamma3dq)

	Vector3D p0(xi, 0., 0.);     //xi in [-Lx/2, Lx/2]
	Vector sfx(NS(), 1); GetShapesx(sfx, p0);
	Vector sfy(NS(), 1); GetShapesy(sfy, p0);
	Vector sfz(NS(), 1); GetShapesz(sfz, p0);

	dgammadq.SetSize(SOS(), Dim());
	dgammadq.SetAll(0.);

	Vector3D drdx, drdy, drdz;
	drdx = GetPosx3D(p0);
	drdy = GetPosy3D(p0);
	drdz = GetPosz3D(p0);

	// help variables
	double norm_drdy_squ = drdy.Norm2();
	double norm_drdy_lin = sqrt(norm_drdy_squ);
	double norm_drdy_cub = norm_drdy_lin * norm_drdy_squ;
	double norm_drdz_squ = drdz.Norm2();
	double norm_drdz_lin = sqrt(norm_drdz_squ);
	double norm_drdz_cub = norm_drdz_lin * norm_drdz_squ;
	
	for (int i=1; i<=Dim(); i++)
	{
		Vector3D ei = Getei(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z());
		TArray<Vector3D> deidp(GetdeidDiffVar(i, drdy.X(), drdy.Y(), drdy.Z(), drdz.X(), drdz.Y(), drdz.Z()));
		double norm_ei_squ = ei.Norm2();
		double norm_ei_lin = sqrt(norm_ei_squ);
		double norm_ei_cub = norm_ei_lin * norm_ei_squ;

		gamma(i) = (ei*drdx)/norm_ei_lin;
		Vector3D v = (1/norm_ei_lin)*drdx - (gamma(i)/norm_ei_squ)*ei;
		//double facy = v * (de1dp(1) + de1dp(2) + de1dp(3));
		//double facz = v * (de1dp(4) + de1dp(5) + de1dp(6));

		//assembling of dgamma1dq
		for (int j=1; j<=NS(); j++)
		{
			for (int k=1; k<=Dim(); k++)
			{
				double rowidx = (j-1)*Dim() + k;
				dgammadq(rowidx, i) = ((1/norm_ei_lin)*ei(k))*sfx(j);
				dgammadq(rowidx, i) += (v*deidp(k))*sfy(j);
				dgammadq(rowidx, i) += (v*deidp(k+Dim()))*sfz(j);			
			}
		}
	}
	gamma(1) -= 1;
}


void ANCFBeamShear3DGeneric::GetKappa(const double& xi, Vector3D& kappa, int flagD)		// -Lx/2 <= xi <= Lx/2
{
	Vector3D p0 = Vector3D(xi, 0., 0.);
	kappa = GetTwistAndCurvature(p0, 0, flagD) - GetTwistAndCurvature(p0, 1, flagD);
}


void ANCFBeamShear3DGeneric::GetDeltaKappa(const double& xi, Matrix& dkappadq, Vector3D& kappa)   //derivative of vector Kappa wrt vector q at xi (which is between -Lx/2 and Lx/2)
{
	//TMStartTimer(28);
	Vector3D p0 = Vector3D(xi, 0, 0);

	Vector sfy(NS(), 1);  GetShapesy(sfy, p0);
	Vector sfz(NS(), 1);  GetShapesz(sfz, p0);
	Vector sfxy(NS(), 1); GetShapesxy(sfxy, p0);
	Vector sfxz(NS(), 1); GetShapesxz(sfxz, p0);


	Vector3D drdy = GetPosy3D(p0);
	Vector3D drdz = GetPosz3D(p0);
	Vector3D ddrdxdy = GetPosxy3D(p0);
	Vector3D ddrdxdz = GetPosxz3D(p0);

	// e1 = drdy x drdz, e3 = drdz, e2 = e3 x e1 = -e1 x e3
	Vector3D ei;
	Vector3D deidxi;
	TArray<Vector3D> deidp(6);
	TArray<Vector3D> ddeidxidp(12);

	// ei     depends on p = (dr1dy dr2dy dr3dy dr1dz dr2dz dr3dz),
	// deidxi depends on p = (dr1dy dr2dy dr3dy dr1dz dr2dz dr3dz ddr1dxdy ddr2dxdy ddr3dxdy ddr1dxdz ddr2dxdz ddr3dxdz)
	// delta kappa = delta k - delta(Rot*k_0)

	// calculation of delta k:  in the way  dkdq = dkdp * dpdq

	TArray<Vector3D> dkdp(12);    // Vector3D() sets components to 0!
	for (int i=1; i<=Dim(); i++) 
	{
		ei = Getei(i, 
			drdy.X(), drdy.Y(), drdy.Z(), 
			drdz.X(), drdz.Y(), drdz.Z());
		deidxi = Getdeidxi(i, 
			drdy.X(), drdy.Y(), drdy.Z(), 
			drdz.X(), drdz.Y(), drdz.Z(), 
			ddrdxdy.X(), ddrdxdy.Y(), ddrdxdy.Z(), 
			ddrdxdz.X(), ddrdxdz.Y(), ddrdxdz.Z());
		deidp = GetdeidDiffVar(i, 
			drdy.X(), drdy.Y(), drdy.Z(), 
			drdz.X(), drdz.Y(), drdz.Z());
		ddeidxidp = GetddeidxidDiffVar(i, 
			drdy.X(), drdy.Y(), drdy.Z(), 
			drdz.X(), drdz.Y(), drdz.Z(), 
			ddrdxdy.X(), ddrdxdy.Y(), ddrdxdy.Z(), 
			ddrdxdz.X(), ddrdxdz.Y(), ddrdxdz.Z());

		double einormsquare = ei*ei;

		for (int j=1; j<=12; j++)
		{
			if (j<=6) 
			{
				dkdp(j) -= ((ei*deidp(j)) / Sqr(einormsquare)) * ei.Cross(deidxi);
				dkdp(j) += (0.5/einormsquare) * deidp(j).Cross(deidxi);
			}

			dkdp(j) += (0.5/einormsquare) * ei.Cross(ddeidxidp(j));
		}
	}

	// delta k = [ d k / d p ] . [ d p / dq ]
	dkappadq.SetAll(0.);
	Vector3D rowvec;

	for(int k=1; k<=NS(); k++)
	{
		for(int j=1; j<=Dim(); j++) 
		{
			rowvec =
				dkdp(j  )*sfy(k) + 
				dkdp(j+3)*sfz(k) + 
				dkdp(j+6)*sfxy(k) + 
				dkdp(j+9)*sfxz(k); 
			int rowidx = j + (k-1)*Dim();

			dkappadq(rowidx, 1) = rowvec.X();
			dkappadq(rowidx, 2) = rowvec.Y();
			dkappadq(rowidx, 3) = rowvec.Z();
		}
	}

	// [kappa]_e = A^T (k - k_0)              [..] be the representation in the local frame basis
	// d [kappa]_e / d q = d A^T /dq (k - k_0) + A^T dk/dq

	Matrix helper = dkappadq;
	Matrix3D A = GetRotMatrix(p0);
	Mult(helper, A, dkappadq); //dkappadq = helper * A;

	//Matrix3D AT = A;//.GetTp();
	//for(int i=1; i<=NS()*Dim(); i++)
	//{
	//	Vector3D v;
	//	v(1)=helper(i,1);
	//	v(2)=helper(i,2);
	//	v(3)=helper(i,3);
	//	v = AT*v;
	//	dkappadq(i,1) = v(1);
	//	dkappadq(i,2) = v(2);
	//	dkappadq(i,3) = v(3);
	//}

//TMStopTimer(29);

//TMStartTimer(30);
	kappa = GetTwistAndCurvature(p0, 0) - GetTwistAndCurvature(p0, 1);

//TMStopTimer(30);
//TMStartTimer(20);
	GetdRotvdqT(A*kappa, p0, helper, 1);
//TMStopTimer(20);
	dkappadq += helper;

	//kappa = A.GetTp()*kappa;
}


void ANCFBeamShear3DGeneric::GetThicknessStrain(const double& xi, Vector3D& thicknessstrain, int flagD)		// -Lx/2 <= xi <= Lx/2
{
	//thicknessstrain_yy = 1/2 * (drdy*drdy - 1)
	//thicknessstrain_zz = 1/2 * (drdz*drdz - 1)
	//thicknessstrain_yz_twice = 2 * thicknessstrain_yz = drdy*drdz
	//
	//thicknessstrain = (thicknessstrain_yy, thicknessstrain_zz, thicknessstrain_yz_twice)

	Vector3D p0(xi, 0., 0.);     //xi in [-Lx/2, Lx/2]

	Vector3D drdy = GetPosy3D(p0, flagD);
	Vector3D drdz = GetPosz3D(p0, flagD);

	thicknessstrain.Set(
		(drdy*drdy - 1.)/2.,
		(drdz*drdz - 1.)/2.,
		drdy*drdz );
}

void ANCFBeamShear3DGeneric::GetDeltaThicknessStrain(const double& xi, Matrix& dthicknessstraindq, Vector3D& thicknessstrain)		// -Lx/2 <= xi <= Lx/2
{
	//dthicknessstrainyydq = r,y.S,y
	//dthicknessstrainzzdq = r,z.S,z
	//dthicknessstrainyztwicedq = r,y.S,z + r,z.S,y
	//
	//dthicknessstraindq = (dthicknessstrainyydq, dthicknessstrainzzdq, dthicknessstrainyztwicedq)

	Vector3D p0(xi, 0., 0.);     //xi in [-Lx/2, Lx/2]
	Vector sfy(NS(), 1); GetShapesy(sfy, p0);
	Vector sfz(NS(), 1); GetShapesz(sfz, p0);

	dthicknessstraindq.SetSize(SOS(), Dim());
	dthicknessstraindq.SetAll(0.);

	Vector3D drdy = GetPosy3D(p0);
	Vector3D drdz = GetPosz3D(p0);

	//set thickness strain
	thicknessstrain.Set(
		(drdy*drdy - 1.)/2.,
		(drdz*drdz - 1.)/2.,
		drdy*drdz );

	//assembling of dthicknessstraindq
	for (int i=1; i<=NS(); i++)
	{
		for (int j=1; j<=Dim(); j++)
		{
			double rowidx = (i-1)*Dim() + j;
			dthicknessstraindq(rowidx, 1) += drdy(j)*sfy(i);
			dthicknessstraindq(rowidx, 2) += drdz(j)*sfz(i);
			dthicknessstraindq(rowidx, 3) += drdy(j)*sfz(i) + drdz(j)*sfy(i);
		}
	}
}


//for plasticity:
double ANCFBeamShear3DGeneric::PostNewtonStep(double t)		// changes plastic strains, returns yieldfunction(sigma)/Em
{
	return 0;
}


//---------------------------------------------------------------
//for Displacement/Stress/Strain/Stress Resultants-Visualization:
//---------------------------------------------------------------
double ANCFBeamShear3DGeneric::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{		

	//ploc_scaled = local position: -l/2 <= ploc_scaled <= +l/2
	Vector3D p0 = local_position;//-l/2 <= p0 <= +l/2

	xgd.SetLen(SOS());
	GetDrawCoordinates(xgd);  //xgd=XGD

	//Position:
	if (fvd.VariableType() == FieldVariableDescriptor::FVT_position)
	{
		Vector3D u;
		if(flagD) 
			u = GetPosD(p0);
		else
			u = GetPos(p0);

		return fvd.GetComponent(u);
	}

	//Displacement:
	if (fvd.VariableType() == FieldVariableDescriptor::FVT_displacement)
	{
		Vector3D u;
		if(flagD) 
			u = GetDisplacementD(p0);
		else
			u = GetDisplacement(p0);

		return fvd.GetComponent(u);
	}

	//Velocity:
	if (fvd.VariableType() == FieldVariableDescriptor::FVT_velocity)
	{
		Vector3D u;
		u = GetVel3D(p0, flagD);

		return fvd.GetComponent(u);
	}

	//Stress and Strain:
	if (InContinuumMechanicsFormulation())
	{
		Matrix3D strain(0), F, piola1;
		ConstMatrix<36> Dm;
		Dm.SetSize(6,6);
		Dm.SetAll(0);
		Vector3D rx, ry, rz;
		rx=GetPosx3D(p0,flagD);
		ry=GetPosy3D(p0,flagD);
		rz=GetPosz3D(p0,flagD);

		//strain:
		F.SetAll(0); //$JG: WRONG: this is the deformation gradient in local coordinates
		for (int j = 1; j <= Dim(); j++) 
		{
			F(j,1) = rx(j);//grad*xg
			F(j,2) = ry(j);
			F(j,3) = rz(j);
		}
		if(calclinear == 0)
		{
			//F(1,1) += 1;//Positionsgrad (nicht Grad(u), sondern Grad(r))//ux->rx: +1, hier haben wir schon rx, kein +1
			//F(2,2) += 1;
			//F(3,3) += 1;
			//E = 1/2*(F^T F-I)
			strain = 0.5*(F.GetTp()*F);//eps=0.5*(F.GetTp()*F-Id)
			strain(1,1) -= 0.5;
			strain(2,2) -= 0.5;
			strain(3,3) -= 0.5;
		}
		else if(calclinear == 1)
		{
			F(1,1) -= 1;//rx->ux: -1
			F(2,2) -= 1;
			F(3,3) -= 1;
			strain = 0.5*(F + F.GetTp());
		}
		//stress:
		double help3 = Em()* (1-Nu()) / ((1+Nu())*(1-2*Nu()));
		double help4 = Em()*Nu() / ((1+Nu())*(1-2*Nu()));
		double GAky = GetMaterial().BeamGAky();		
		double GAkz = GetMaterial().BeamGAkz();		
		double G = Em()/(2.+2.*Nu());
		double A = GetMaterial().BeamA();

		Dm(1,1)=help3;
		Dm(1,2)=help4;
		Dm(1,3)=help4;

		Dm(2,1)=help4;
		Dm(2,2)=help3;
		Dm(2,3)=help4;

		Dm(3,1)=help4;
		Dm(3,2)=help4;
		Dm(3,3)=help3;

		Dm(4,4)=G;
		Dm(5,5)=G;
		Dm(6,6)=G;

		ConstVector<6> s3(strain(1,1), strain(2,2), strain(3,3), 2.*strain(2,3), 2.*strain(1,3), 2.*strain(1,2));
		ConstVector<6> stress;
		stress = Dm * s3;//Dm=6x6-Matrix, //S = D : E

		//generate 2. Piola-Kirchhoff-tensor, 2.PK = symm. tensor 
		piola1(1,1)=stress(1);
		piola1(2,2)=stress(2);
		piola1(3,3)=stress(3);

		piola1(2,3)=stress(4);
		piola1(3,2)=stress(4);
		piola1(1,3)=stress(5);
		piola1(3,1)=stress(5);
		piola1(1,2)=stress(6);
		piola1(2,1)=stress(6);

		switch(fvd.VariableType())
		{
		case FieldVariableDescriptor::FVT_stress_mises:							return piola1.Mises();
		case FieldVariableDescriptor::FVT_stress:										return fvd.GetComponent(piola1);
		case FieldVariableDescriptor::FVT_total_strain:							return fvd.GetComponent(strain);
		}
	}
	else		// SM-Formulation
	{
		double xi = p0.X();
		Vector3D gamma; GetGamma(xi, gamma, flagD);
		Vector3D kappa; GetKappa(xi, kappa, flagD);
		Vector3D thicknessstrain; GetThicknessStrain(xi, thicknessstrain, flagD);

		//compute interpolated values for gamma, otherwise gamma has oscillatory solution in plot!!!! JG2013-09
		//integration rule must be same as in EvalF2SM() !!!!
		int order_x_bend = 3;						//3		 
		int order_x_thickness_eps = 4;  //4
		int order_x_axial_eps = 3;      //3
		int order_x_shear = 3;					//3

		Vector wT,xT;										//integration points and weights
		GetIntegrationRule(xT,wT,order_x_axial_eps);
		
		double x1 = GetLx()*0.5*xT(1);
		double x2 = GetLx()*0.5*xT(2);
		Vector3D gamma_1; GetGamma(x1, gamma_1, flagD);
//		Vector3D kappa_1; GetKappa(x1, kappa_1, flagD);
		Vector3D gamma_2; GetGamma(x2, gamma_2, flagD);
//		Vector3D kappa_2; GetKappa(x2, kappa_2, flagD);

		int use_interpolated_values = 1;
		//compute interpolated values:
		if (use_interpolated_values)
		{
			gamma = 1./(x2-x1)*((x2-xi)*gamma_1+(xi-x1)*gamma_2); //linear interpolation of all gamma values
		}
		
		double beamEA = GetMaterial().BeamEA();
		double beamGAky = GetMaterial().BeamGAky();
		double beamGAkz = GetMaterial().BeamGAkz();
		double beamGAkyz = BeamGAkyz(beamGAky, beamGAkz);      
		double beamGJkx = GetMaterial().BeamGJkx();
		double beamEIy = GetMaterial().BeamEIy();
		double beamEIz = GetMaterial().BeamEIz();

		Vector3D moment(kappa);
		moment.X() *= beamGJkx;
		moment.Y() *= beamEIy;
		moment.Z() *= beamEIz;

		Vector3D beam_force(gamma);
		beam_force.X() *= beamEA;
		beam_force.Y() *= beamGAky;
		beam_force.Z() *= beamGAkz;

		/*
		//strains and stresses are computed in WRONG WAY JG2013-09 !!!
		//the part of bending strain has been forgotten!!!

		// for sm-formulation the total strain contains:
		//  Gamma_1 Gamma_2 Gamma_3
		//  Gamma_2  E_yy    E_yz
		//  Gamma_3  E_yz    E_zz
		// ... same for the stress

		Matrix3D strains(
			gamma.X(), gamma.Y(), gamma.Z(),
			gamma.Y(), thicknessstrain.X(), 0.5*thicknessstrain.Z(),
			gamma.Z(), 0.5*thicknessstrain.Z(), thicknessstrain.Y()
			);

		Matrix3D stresses(
			beamEA*gamma.X(), beamGAky*gamma.Y(), beamGAkz*gamma.Z(),
			beamGAky*gamma.Y(), beamEA*thicknessstrain.X(), beamGAkyz*0.5*thicknessstrain.Z(),
			beamGAkz*gamma.Z(), beamGAkyz*0.5*thicknessstrain.Z(), beamEA*thicknessstrain.Y()
			);
		*/
		
		switch(fvd.VariableType())
		{
			//beam_curvature is defined as 2-component vector! JG2013-09 
			//case FieldVariableDescriptor::FVT_beam_curvature: return fvd.GetComponent(Vector3D(kappa.X(), kappa.Y(), kappa.Z()));

			//strains and stresses are computed in WRONG WAY JG2013-09 !!!
			//case FieldVariableDescriptor::FVT_total_strain:	  return fvd.GetComponent(strains);
			//case FieldVariableDescriptor::FVT_stress:         return fvd.GetComponent(stresses);

			case FieldVariableDescriptor::FVT_beam_torsion: return kappa.X();
			case FieldVariableDescriptor::FVT_beam_curvature: return fvd.GetComponent(kappa);

			case FieldVariableDescriptor::FVT_beam_moment_torsional: return moment.X();
			case FieldVariableDescriptor::FVT_beam_moment_bending: return fvd.GetComponent(moment);

			case FieldVariableDescriptor::FVT_beam_axial_extension: return gamma.X();
			case FieldVariableDescriptor::FVT_beam_shear: return fvd.GetComponent(gamma);

			case FieldVariableDescriptor::FVT_beam_force_axial: return beam_force.X();
			case FieldVariableDescriptor::FVT_beam_force_transversal: return fvd.GetComponent(beam_force);
		}
	}

	return FIELD_VARIABLE_NO_VALUE;
}


void ANCFBeamShear3DGeneric::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_z);
	FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement, FieldVariableDescriptor::FVCI_z);
	if (InContinuumMechanicsFormulation())
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain,
			FieldVariableDescriptor::FVCI_z, true);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress,
			FieldVariableDescriptor::FVCI_z, true);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress_mises);
	}
	else
	{
		
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_curvature, FieldVariableDescriptor::FVCI_z);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_torsion);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_bending, FieldVariableDescriptor::FVCI_z);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment_torsional);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_shear, FieldVariableDescriptor::FVCI_z);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_axial_extension);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_transversal, FieldVariableDescriptor::FVCI_z);
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_force_axial);

		//this is the global beam moment! JG2013-09:
		//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_beam_moment, FieldVariableDescriptor::FVCI_z);
		// for sm-formulation the total strain contains:
		//  Gamma_1 Gamma_2 Gamma_3
		//  Gamma_2  E_yy    E_yz
		//  Gamma_3  E_yz    E_zz
		// ... same for the stress
		
		//strains and stresses are computed in WRONG WAY JG2013-09 !!!
		//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_total_strain, FieldVariableDescriptor::FVCI_z, true);   
		//FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_stress, FieldVariableDescriptor::FVCI_z, true);

	}
}

//ACHTUNG: nur kopiert von ANCFBeam3D und ComputeStress-Aufrufe auskommentiert
void ANCFBeamShear3DGeneric::DrawElement()
{
	mbs->SetColor(col);

	double lx1 = 1; double ly1 = 1/**GetMBS()->GetMagnifyYZ()*/; double lz1 = 1/**GetMBS()->GetMagnifyYZ()*/; //GetMagnifyYZ() leads to wrong stress evaluation
	double def_scale = GetMBS()->GetDOption(105); //deformation scaling

	int linemode = 1; //0=no lines, 1=outline+color, 2=outline, 3=elementline+color, 4=elementline
	if (GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 2;
	}
	if (!GetMBS()->GetIOption(110) && GetMBS()->GetIOption(111))
	{
		linemode = 0;
	}
	if (!GetMBS()->GetIOption(110) && !GetMBS()->GetIOption(111))
	{
		linemode = 4;
	}

	int colormode = 0;
	if (GetMBS()->GetActualPostProcessingFieldVariable() != NULL) colormode = 1;
	//if (GetMBS()->GetDrawMode() != 0) colormode = 1;//Yuri
	else
	{
		colormode = 0; 
	}

	if (colormode)
	{
		double modeval = 0;
		int xgset = 0;

		double tilex = GetMBS()->GetIOption(137);
		double tiley = GetMBS()->GetIOption(138);
		TArray<Vector3D> points((int)(tilex+1)*(int)(tiley+1));
		TArray<double> vals((int)(tilex+1)*(int)(tiley+1));
		double v=0;

		for (int side=1; side <= 6; side++)
		{
			points.SetLen(0); vals.SetLen(0);
			Vector3D p0, vx, vy;
			int tileyn = (int)tiley;
			int tilexn = (int)tilex;

			if (side == 1)
			{ //bottom
				p0 = Vector3D(-lx1,-ly1,-lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 2)
			{ //top
				p0 = Vector3D(-lx1, ly1, lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,0,-2.*lz1/tileyn);
			}
			else if (side == 3)
			{ //front
				p0 = Vector3D(-lx1,-ly1, lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,2.*ly1/tileyn,0);
			}
			else if (side == 4)
			{ //back
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(2.*lx1/tilexn,0,0);
				vy = Vector3D(0,-2.*ly1/tileyn,0);
			}
			else if (side == 5)
			{ //left
				p0 = Vector3D(-lx1, ly1,-lz1);
				vx = Vector3D(0,-2.*ly1/tilexn,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}
			else if (side == 6)
			{ //right
				p0 = Vector3D( lx1,-ly1,-lz1);
				vx = Vector3D(0,2.*ly1/tilexn,0);
				vy = Vector3D(0,0,2.*lz1/tileyn);
			}

			for (double iy = 0; iy <= tileyn+1e-10; iy++)
			{
				for (double ix = 0; ix <= tilexn+1e-10; ix++)
				{
					Vector3D ploc = (p0+ix*vx+iy*vy);//-1<=p0<=+1
					Vector3D ploc_scaled(ploc(1)*GetLx()*0.5,ploc(2)*GetLy()*0.5,ploc(3)*GetLz()*0.5);//-l/2<=ploc_scaled<=+l/2
					Vector3D pg = GetPos3D0D(ploc, def_scale);
					//Vector3D p(pg);//p=pg
					points.Add(pg);
						if (colormode)
							v = GetFieldVariableValue(*GetMBS()->GetActualPostProcessingFieldVariable(), ploc_scaled, true);
						vals.Add(v);

					//if (colormode) ComputeStressD(ploc,type,comp1,comp2,v);//old_KN
					//vals.Add(v);
				}
			}
			mbs->DrawColorQuads(points,vals,(int)tilexn+1,(int)tileyn+1,colormode,linemode);
			//mbs->UO()<<"points"<<points(1)<<"\n";
		}
	}
	else
	{
		int drawgrid = GetMBS()->GetIOption(110);
		double thickness = GetMBS()->GetDOption(102);
		double tiling = GetMBS()->GetIOption(137);

		mbs->SetDrawlines(0);
		mbs->SetDrawlinesH(0);
		Vector3D p1,p2,p3,p4,p5,p6,p7,p8;
		for (double i = 0; i < tiling; i++)
		{
			double l1 = -lx1+2*lx1*i/tiling;
			double l2 = -lx1+2*lx1*(i+1)/tiling;
			if (i == 0)
			{
				//p8 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1,-lz1)));
				//p7 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1, lz1)));
				//p4 = Vector3D(GetPos3D0D(Vector3D(l1, ly1,-lz1))); 
				//p3 = Vector3D(GetPos3D0D(Vector3D(l1, ly1, lz1)));
				p8 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1,-lz1),def_scale));
				p7 = Vector3D(GetPos3D0D(Vector3D(l1,-ly1, lz1),def_scale));
				p4 = Vector3D(GetPos3D0D(Vector3D(l1, ly1,-lz1),def_scale)); 
				p3 = Vector3D(GetPos3D0D(Vector3D(l1, ly1, lz1),def_scale));
				if (drawgrid)
				{
					mbs->MyDrawLine(p8,p7,thickness);
					mbs->MyDrawLine(p7,p3,thickness);
					mbs->MyDrawLine(p4,p3,thickness);
					mbs->MyDrawLine(p4,p8,thickness);
				}
			}
			else
			{
				p8 = p6;
				p7 = p5;
				p4 = p2;
				p3 = p1;
			}
			//p6 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1,-lz1)));
			//p5 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1, lz1)));
			//p2 = Vector3D(GetPos3D0D(Vector3D(l2, ly1,-lz1)));
			//p1 = Vector3D(GetPos3D0D(Vector3D(l2, ly1, lz1)));
			p6 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1,-lz1),def_scale));
			p5 = Vector3D(GetPos3D0D(Vector3D(l2,-ly1, lz1),def_scale));
			p2 = Vector3D(GetPos3D0D(Vector3D(l2, ly1,-lz1),def_scale));
			p1 = Vector3D(GetPos3D0D(Vector3D(l2, ly1, lz1),def_scale));
			int dout = 0;
			if (i == 0 || i == tiling-1) dout = 1;
			if (drawgrid) //draw mesh
			{
				mbs->MyDrawLine(p6,p5,thickness);
				mbs->MyDrawLine(p5,p1,thickness);
				mbs->MyDrawLine(p2,p1,thickness);
				mbs->MyDrawLine(p2,p6,thickness);
				mbs->MyDrawLine(p6,p8,thickness);
				mbs->MyDrawLine(p5,p7,thickness);
				mbs->MyDrawLine(p2,p4,thickness);
				mbs->MyDrawLine(p1,p3,thickness);
			}
			if (GetMBS()->GetIOption(111)) //draw solid
			{
				mbs->DrawHex(p1,p2,p3,p4,p5,p6,p7,p8,dout);
			}
		}
		mbs->SetDrawlines(0);

		//// draw slope vectors
		//Vector3D position, slope1, slope2;
		//double drawlength = max(ly,lz)*2.;
		//double linethickness = max(2.,0.04);
		//Vector3D color(0., 0.5, 0.);
		//
		//position = GetPos3D0D(Vector3D(-1., 0., 0.));   //left node
		//slope1(1) = XGD(4) + q0(4);
		//slope1(2) = XGD(5) + q0(5);
		//slope1(3) = XGD(6) + q0(6);
		//slope1.Normalize();
		//slope1 *= drawlength;
		//mbs->MyDrawLine(position, position+slope1, linethickness, color);
		//slope2(1) = XGD(7) + q0(7);
		//slope2(2) = XGD(8) + q0(8);
		//slope2(3) = XGD(9) + q0(9);
		//slope2.Normalize();
		//slope2 *= drawlength;
		//mbs->MyDrawLine(position, position+slope2, linethickness, color);

		//position = GetPos3D0D(Vector3D(1., 0., 0.));   //right node
		//slope1(1) = XGD(13) + q0(13);
		//slope1(2) = XGD(14) + q0(14);
		//slope1(3) = XGD(15) + q0(15);
		//slope1.Normalize();
		//slope1 *= drawlength;
		//mbs->MyDrawLine(position, position+slope1, linethickness, color);
		//slope2(1) = XGD(16) + q0(16);
		//slope2(2) = XGD(17) + q0(17);
		//slope2(3) = XGD(18) + q0(18);
		//slope2.Normalize();
		//slope2 *= drawlength;
		//mbs->MyDrawLine(position, position+slope2, linethickness, color);

		//if (nnodes == 3)
		//{
		//	position = GetPos3D0D(Vector3D(0., 0., 0.));   //middle node
		//	slope1(1) = XGD(22) + q0(22);
		//	slope1(2) = XGD(23) + q0(23);
		//	slope1(3) = XGD(24) + q0(24);
		//	slope1.Normalize();
		//	slope1 *= drawlength;
		//	mbs->MyDrawLine(position, position+slope1, linethickness, color);
		//	slope2(1) = XGD(25) + q0(25);
		//	slope2(2) = XGD(26) + q0(26);
		//	slope2(3) = XGD(27) + q0(27);
		//	slope2.Normalize();
		//	slope2 *= drawlength;
		//	mbs->MyDrawLine(position, position+slope2, linethickness, color);
		//}
	}
};

#pragma endregion

#pragma region ANCFBeamShear3DLinear

void ANCFBeamShear3DLinear::ElementDefaultConstructorInitialization() 
{
	// set nodes
	nnodes=2;
	n1=1; n2=2;
	
	//set reference configuration
	q0.SetLen(SOS());
	q0.SetAll(0.);
	size = Vector3D(1,0.1,0.1);
	penalty.SetLen(9);
	penalty.SetAll(1);
	x_init.SetLen(2*SOS());
	x_init.SetAll(0.);

	xg.SetLen(SOS());
	xg.SetAll(0.);
	
	is_straight_beam_in_reference_configuration = 0;
}

// default set function!!!
void ANCFBeamShear3DLinear::SetANCFBeamShear3DLinear(int n1i, int n2i, int materialnumi, const Vector3D& coli)
{
	ElementDefaultConstructorInitialization();

	// set nodes
	n1=n1i; n2=n2i;

	// set common stuff for ancfbeamshear3d elements
	SetANCFBeamShear3DGeneric(materialnumi, coli);

	//used for computation speedup: test if beam element is straight in reference configuration
	is_straight_beam_in_reference_configuration = 1; // linear (2-noded) element is always straight, only cross section may vary
};

// deprecated set function!!!
void ANCFBeamShear3DLinear::SetANCFBeamShear3DLinear(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
																				 const Vector3D& si, const Vector3D& coli)
{
	ElementDefaultConstructorInitialization();
	// set nodes
	//nnodes=2;		// is now in ElementDefaultConstructorInitialization
	n1=n1i; n2=n2i;
	
	//set reference configuration
	//q0.SetLen(SOS()); // is now in ElementDefaultConstructorInitialization
	//q0.SetAll(0.); // is now in ElementDefaultConstructorInitialization
	int i;
	int ndofs_per_node = SOS()/nnodes;    //each of the nodes owns ndofs_per_node degrees of freedom
	for(int i=1;i<=ndofs_per_node;i++)   
	{
		q0(i) = xc1(i);
		q0(i+ndofs_per_node) = xc2(i);
	}

	// set common stuff for ancfbeamshear3d elements
	SetANCFBeamShear3DGeneric(materialnumi, si, coli);

	//used for computation speedup: test if beam element is straight in reference configuration
	is_straight_beam_in_reference_configuration = 1; // linear (2-noded) element is always straight, only cross section may vary
};

void ANCFBeamShear3DLinear::LinkToElements()
{
	if (SOSowned() == 0)
	{
		//UO() << "Link nodes to elements in ANCFBeamShear3DLinear\n";
		const Node& node1 = GetMBS()->GetNode(n1);
		const Node& node2 = GetMBS()->GetNode(n2);
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i));
		}
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i+node1.SOS()));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i+node2.SOS()));
		}
	}
}

//-------------------------------------------------------
//Shapefunctions
//
//defined on scaled rectangular element: ploc=(\xi,\eta,\zeta)
//-L/2<\xi<+L/2, -H/2<\eta<+H/2, -B/2<\zeta<+B/2
//-------------------------------------------------------


void ANCFBeamShear3DLinear::GetShapes(Vector& sf, const Vector3D& ploc) const 
//compute all Shapefunctions at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = ploc(1)/GetLx();
	double sf1 = .5 - fact;
	double sf4 = .5 + fact;
	sf(1) = sf1;
	sf(2) = sf1*ploc(2);
	sf(3) = sf1*ploc(3);
	sf(4) = sf4;
	sf(5) = sf4*ploc(2);
	sf(6) = sf4*ploc(3);
}

double ANCFBeamShear3DLinear::GetSF(int i, const Vector3D& ploc) const
//computes i-th Shapefunction at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return (.5 - ploc(1)/GetLx());
	case 2: return (.5 - ploc(1)/GetLx())*ploc(2);
	case 3: return (.5 - ploc(1)/GetLx())*ploc(3);
	case 4: return (.5 + ploc(1)/GetLx());
	case 5: return (.5 + ploc(1)/GetLx())*ploc(2);
	case 6: return (.5 + ploc(1)/GetLx())*ploc(3);
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DLinear::GetShapesx(Vector& sf, const Vector3D& ploc) const 
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = 1./GetLx();
	sf(1) = -fact;
	sf(2) = -fact*ploc(2);
	sf(3) = -fact*ploc(3);
	sf(4) = fact;
	sf(5) = fact*ploc(2);
	sf(6) = fact*ploc(3);
}

double ANCFBeamShear3DLinear::GetSFx(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(1)=xi at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return -1./GetLx();
	case 2: return -ploc(2)/GetLx();
	case 3: return -ploc(3)/GetLx();
	case 4: return 1./GetLx();
	case 5: return ploc(2)/GetLx();
	case 6: return ploc(3)/GetLx();
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DLinear::GetShapesy(Vector& sf, const Vector3D& ploc) const 
//compute derivative of all Shapefunctions w.r.t. ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = ploc(1)/GetLx();
	sf(1) = 0.;
	sf(2) = .5 - fact;
	sf(3) = 0.;
	sf(4) = 0.;
	sf(5) = .5 + fact;
	sf(6) = 0.;
}

double ANCFBeamShear3DLinear::GetSFy(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return 0.;
	case 2: return .5 - ploc(1)/GetLx();
	case 3: return 0.;
	case 4: return 0.;
	case 5: return .5 + ploc(1)/GetLx();
	case 6: return 0.;
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DLinear::GetShapesz(Vector& sf, const Vector3D& ploc) const 
//compute derivative of all Shapefunctions w.r.t. ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = ploc(1)/GetLx();
	sf(1) = 0.;
	sf(2) = 0.;
	sf(3) = .5 - fact;
	sf(4) = 0.;
	sf(5) = 0.;
	sf(6) = .5 + fact;
}

double ANCFBeamShear3DLinear::GetSFz(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return 0.;
	case 2: return 0.;
	case 3: return .5 - ploc(1)/GetLx();
	case 4: return 0.;
	case 5: return 0.;
	case 6: return .5 + ploc(1)/GetLx();
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DLinear::GetShapesxy(Vector& sf, const Vector3D& ploc) const 
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi and ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = 1./GetLx();
	sf(1) = 0.;
	sf(2) = -fact;
	sf(3) = 0.;
	sf(4) = 0.;
	sf(5) = fact;
	sf(6) = 0.;
}

double ANCFBeamShear3DLinear::GetSFxy(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. plox(1)=xi and ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return 0.;
	case 2: return -1./GetLx();
	case 3: return 0.;
	case 4: return 0.;
	case 5: return 1./GetLx();
	case 6: return 0.;
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DLinear::GetShapesxz(Vector& sf, const Vector3D& ploc) const 
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi and ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());
	double fact = 1./GetLx();
	sf(1) = 0.;
	sf(2) = 0.;
	sf(3) = -fact;
	sf(4) = 0.;
	sf(5) = 0.;
	sf(6) = fact;
}

double ANCFBeamShear3DLinear::GetSFxz(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. plox(1)=xi and ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return 0.;
	case 2: return 0.;
	case 3: return -1./GetLx();
	case 4: return 0.;
	case 5: return 0.;
	case 6: return 1./GetLx();
	default: mbs->UO()<<"only 6 Shapefunctions\n"; return 0.;
	}
	return 0.;
}
#pragma endregion



#pragma region ANCFBeamShear3DQuadratic

void ANCFBeamShear3DQuadratic::ElementDefaultConstructorInitialization() 
{
	// set nodes
	nnodes=3;
	n1=1; n2=2; n3 = 3;
	
	//set reference configuration
	q0.SetLen(SOS());
	q0.SetAll(0.);

	//lx = 1; 
	//ly = 0.1;
	//lz = 0.1;
	size = Vector3D(1,0.1,0.1);

	penalty.SetLen(9);
	penalty.SetAll(1);

	x_init.SetLen(2*SOS());
	x_init.SetAll(0.);

	xg.SetLen(SOS());
	xg.SetAll(0.);
	
	is_straight_beam_in_reference_configuration = 0;
}

// default set function
void ANCFBeamShear3DQuadratic::SetANCFBeamShear3DQuadratic(int n1i, int n2i, int n3i, int materialnumi, const Vector3D& coli)
{
	ElementDefaultConstructorInitialization();
	// set nodes
	//nnodes=3;	// is now in ElementDefaultConstructorInitialization()
	n1=n1i; n2=n2i, n3=n3i;
	
	// set common stuff for ancfbeamshear3d elements
	SetANCFBeamShear3DGeneric(materialnumi, coli);

	//used for computation speedup: test if beam element is straight in reference configuration
	Vector3D pos_left(q0(1), q0(2), q0(3));
	Vector3D pos_right(q0(10), q0(11), q0(12));
	Vector3D pos_mid(q0(19), q0(20), q0(21));
	Vector3D v1 = pos_right - pos_mid;
	Vector3D v2 = pos_mid - pos_left;
	v1.Normalize();
	v2.Normalize();
	is_straight_beam_in_reference_configuration = (v1*v2 > 1-1e-14);   // if both sections of the beam are colinear, then this scalar product gives numerically 1
};

// deprecated set function!!!
void ANCFBeamShear3DQuadratic::SetANCFBeamShear3DQuadratic(const Vector& xc1, const Vector& xc2, const Vector& xc3, int n1i, int n2i, int n3i, int materialnumi,
																				 const Vector3D& si, const Vector3D& coli)
{
	ElementDefaultConstructorInitialization();
	// set nodes
	//nnodes=3;	// is now in ElementDefaultConstructorInitialization()
	n1=n1i; n2=n2i, n3=n3i;
	
	//set reference configuration
	//q0.SetLen(SOS()); // is now in ElementDefaultConstructorInitialization()
	//q0.SetAll(0.); // is now in ElementDefaultConstructorInitialization()
	int i;
	int ndofs_per_node = SOS()/nnodes;    //each of the nodes owns ndofs_per_node degrees of freedom
	for(int i=1;i<=ndofs_per_node;i++)   
	{
		q0(i) = xc1(i);
		q0(i+ndofs_per_node) = xc2(i);
		q0(i+2*ndofs_per_node) = xc3(i);
	}

	//used for computation speedup: test if beam element is straight in reference configuration
	Vector3D pos_left(q0(1), q0(2), q0(3));
	Vector3D pos_right(q0(10), q0(11), q0(12));
	Vector3D pos_mid(q0(10), q0(11), q0(12));
	Vector3D v1 = pos_right - pos_mid;
	Vector3D v2 = pos_mid - pos_left;
	v1.Normalize();
	v2.Normalize();
	is_straight_beam_in_reference_configuration = (v1*v2 > 1-1e-14);   // if both sections of the beam are colinear, then this scalar product gives numerically 1

	// set common stuff for ancfbeamshear3d elements
	SetANCFBeamShear3DGeneric(materialnumi, si, coli);
};

void ANCFBeamShear3DQuadratic::LinkToElements()
{
	if (SOSowned() == 0)
	{
		//UO() << "Link nodes to elements in ANCFBeamShear3DQuadratic\n";
		const Node& node1 = GetMBS()->GetNode(n1);
		const Node& node2 = GetMBS()->GetNode(n2);
		const Node& node3 = GetMBS()->GetNode(n3);
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i));
		}
		for (int i=1; i <= node3.SOS(); i++)
		{
			AddLTG(node3.Get(i));
		}
		for (int i=1; i <= node1.SOS(); i++)
		{
			AddLTG(node1.Get(i+node1.SOS()));
		}
		for (int i=1; i <= node2.SOS(); i++)
		{
			AddLTG(node2.Get(i+node2.SOS()));
		}
		for (int i=1; i <= node3.SOS(); i++)
		{
			AddLTG(node3.Get(i+node3.SOS()));
		}
	}
}

//-------------------------------------------------------
//Shapefunctions
//
//defined on scaled rectangular element: ploc=(\xi,\eta,\zeta)
//-L/2<\xi<+L/2, -H/2<\eta<+H/2, -B/2<\zeta<+B/2
//-------------------------------------------------------


void ANCFBeamShear3DQuadratic::GetShapes(Vector& sf, const Vector3D& ploc) const
//compute all Shapefunctions at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact=2./(GetLx()*GetLx());
	sf(1)=-fact*(GetLx()/2.-ploc(1))*ploc(1);
	sf(2)=sf(1)*ploc(2);
	sf(3)=sf(1)*ploc(3);
	sf(4)=fact*(GetLx()/2.+ploc(1))*ploc(1);
	sf(5)=sf(4)*ploc(2);
	sf(6)=sf(4)*ploc(3);
	sf(7)=-2.*fact*(GetLx()/2.+ploc(1))*(ploc(1)-GetLx()/2.);
	sf(8)=sf(7)*ploc(2);
	sf(9)=sf(7)*ploc(3);
}

double ANCFBeamShear3DQuadratic::GetSF(int i, const Vector3D& ploc) const
//computes i-th Shapefunction at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{ 
	switch(i)
	{
	case 1: return -2./(GetLx()*GetLx())*ploc(1)*(GetLx()/2.-ploc(1));
	case 2: return -2./(GetLx()*GetLx())*ploc(1)*ploc(2)*(GetLx()/2.-ploc(1));
	case 3: return -2./(GetLx()*GetLx())*ploc(1)*ploc(3)*(GetLx()/2.-ploc(1));
	case 4: return 2./(GetLx()*GetLx())*ploc(1)*(GetLx()/2.+ploc(1));
	case 5: return 2./(GetLx()*GetLx())*ploc(1)*ploc(2)*(GetLx()/2.+ploc(1));
	case 6: return 2./(GetLx()*GetLx())*ploc(1)*ploc(3)*(GetLx()/2.+ploc(1));
	case 7: return -4./(GetLx()*GetLx())*(ploc(1)-GetLx()/2.)*(ploc(1)+GetLx()/2.);
	case 8: return -4./(GetLx()*GetLx())*ploc(2)*(ploc(1)-GetLx()/2.)*(ploc(1)+GetLx()/2.);
	case 9: return -4./(GetLx()*GetLx())*ploc(3)*(ploc(1)-GetLx()/2.)*(ploc(1)+GetLx()/2.);
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DQuadratic::GetShapesx(Vector& sf, const Vector3D& ploc) const
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact = 4.*ploc(1)/(GetLx()*GetLx());
	double sf1 = -1./GetLx()+fact;
	double sf4 = 1./GetLx()+fact;
	double sf7 = -2.*fact;
	sf(1) = sf1;
	sf(2) = sf1*ploc(2);
	sf(3) = sf1*ploc(3);
	sf(4) = sf4;
	sf(5) = sf4*ploc(2);
	sf(6) = sf4*ploc(3);
	sf(7) = sf7;
	sf(8) = sf7*ploc(2);
	sf(9) = sf7*ploc(3);
}

double ANCFBeamShear3DQuadratic::GetSFx(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(1)=xi at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{
	switch(i)
	{
	case 1: return -1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 2: return -ploc(2)/GetLx()+4.*ploc(1)*ploc(2)/(GetLx()*GetLx());
	case 3: return -ploc(3)/GetLx()+4.*ploc(1)*ploc(3)/(GetLx()*GetLx());
	case 4: return 1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 5: return ploc(2)/GetLx()+4.*ploc(1)*ploc(2)/(GetLx()*GetLx());
	case 6: return ploc(3)/GetLx()+4.*ploc(1)*ploc(3)/(GetLx()*GetLx());
	case 7: return -8.*ploc(1)/(GetLx()*GetLx());
	case 8: return -8.*ploc(1)*ploc(2)/(GetLx()*GetLx());
	case 9: return -8.*ploc(1)*ploc(3)/(GetLx()*GetLx());
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DQuadratic::GetShapesy(Vector& sf, const Vector3D& ploc) const
//compute derivative of all Shapefunctions w.r.t. ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact1 = ploc(1)/GetLx();
	double fact2 = 2.*fact1*fact1;
	sf(1) = 0;
	sf(2) = -fact1 + fact2;
	sf(3) = 0;
	sf(4) = 0;
	sf(5) = fact1 + fact2;
	sf(6) = 0;
	sf(7) = 0;
	sf(8) = 1. - 2.*fact2;
	sf(9) = 0;
}

double ANCFBeamShear3DQuadratic::GetSFy(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{
	switch(i)
	{
	case 1: return 0;
	case 2: return -ploc(1)/GetLx()+2.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	case 3: return 0;
	case 4: return 0;
	case 5: return ploc(1)/GetLx()+2.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	case 6: return 0;
	case 7: return 0;
	case 8: return 1.-4.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	case 9: return 0;
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DQuadratic::GetShapesxy(Vector& sf, const Vector3D& ploc) const
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi and ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact1 = 1./GetLx();
	double fact2 = 4.*ploc(1)*fact1*fact1;
	sf(1) = 0;
	sf(2) = -fact1 + fact2;
	sf(3) = 0;
	sf(4) = 0;
	sf(5) = fact1 + fact2;
	sf(6) = 0;
	sf(7) = 0;
	sf(8) = -2.*fact2;
	sf(9) = 0;
}

double ANCFBeamShear3DQuadratic::GetSFxy(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(1)=xi and ploc(2)=eta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{
	switch(i)
	{
	case 1: return 0;
	case 2: return -1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 3: return 0;
	case 4: return 0;
	case 5: return 1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 6: return 0;
	case 7: return 0;
	case 8: return -8.*ploc(1)/(GetLx()*GetLx());
	case 9: return 0;
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DQuadratic::GetShapesz(Vector& sf, const Vector3D& ploc) const
//compute derivative of all Shapefunctions w.r.t. ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact1 = ploc(1)/GetLx();
	double fact2 = 2.*fact1*fact1;
	sf(1) = 0;
	sf(2) = 0;
	sf(3) = -fact1 + fact2;
	sf(4) = 0;
	sf(5) = 0;
	sf(6) = fact1 + fact2;
	sf(7) = 0;
	sf(8) = 0;
	sf(9) = 1. - 2.*fact2;
}

double ANCFBeamShear3DQuadratic::GetSFz(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2] 
{
	switch(i)
	{
	case 1: return 0;
	case 2: return 0;
	case 3: return -ploc(1)/GetLx()+2.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	case 4: return 0;
	case 5: return 0;
	case 6: return ploc(1)/GetLx()+2.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	case 7: return 0;
	case 8: return 0;
	case 9: return 1.-4.*ploc(1)*ploc(1)/(GetLx()*GetLx());
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

void ANCFBeamShear3DQuadratic::GetShapesxz(Vector& sf, const Vector3D& ploc) const
//compute derivative of all Shapefunctions w.r.t. ploc(1)=xi and ploc(3)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{ 
	sf.SetLen(NS());  //NS=3*nnodes, vector sf has 9 entries
	double fact1 = 1./GetLx();
	double fact2 = 4.*ploc(1)*fact1*fact1;
	sf(1) = 0;
	sf(2) = 0;
	sf(3) = -fact1 + fact2;
	sf(4) = 0;
	sf(5) = 0;
	sf(6) = fact1 + fact2;
	sf(7) = 0;
	sf(8) = 0;
	sf(9) = -2.*fact2;
}

double ANCFBeamShear3DQuadratic::GetSFxz(int i, const Vector3D& ploc) const
//computes derivative of i-th Shapefunction w.r.t. ploc(1)=xi and ploc(2)=zeta at ploc=(xi,eta,zeta) in [-L/2,L/2] x [-h/2,h/2] x [-w/2,w/2]
{
	switch(i)
	{
	case 1: return 0;
	case 2: return 0;
	case 3: return -1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 4: return 0;
	case 5: return 0;
	case 6: return 1./GetLx()+4.*ploc(1)/(GetLx()*GetLx());
	case 7: return 0;
	case 8: return 0;
	case 9: return -8.*ploc(1)/(GetLx()*GetLx());
	default: mbs->UO()<<"only 9 Shapefunctions\n"; return 0.;
	}
	return 0.;
}

#pragma endregion
