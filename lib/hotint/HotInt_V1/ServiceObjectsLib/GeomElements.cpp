//#**************************************************************
//#
//# filename:             drawtools.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						May 2005
//# description:          Spatial objects 3D
//#                       
//# remarks:						  
//#
///# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian
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
//#**************************************************************

#include "geomElements.h"
#include "elementDataAccess.h"
#include "element.h"
#include "rigid2d.h"		// is needed for ToP3D in case of a 2D body
#include "body3d.h"		// is needed for GetRotMatrixD(localPosition)
#include "..\mbs_interface\renderContext.h"
#include "myfile.h"
#include "graphicsConstants.h"
#include "mystring.h"

void GeomElement::ElementDefaultConstructorInitialization()
{
	col = Vector3D(0.2,0.2,0.8);
	elnum = 0; 
	objdata = 0; 
	draw_dim = Vector3D(1,1,1);
	transparency = -1;
	translated = Vector3D(0,0,0); 
	name = GetElementSpec();
	drawstyle = TDSDrawOutline | TDSFillAreas;
	pointsize = 0.1;
	linethickness = 1.;
}

void GeomElement::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	ed.SetText(GetElementSpec(), "geom_element_type"); ed.SetToolTipText("specification of GeomElement type. Once the element is added to the mbs, you MUST NOT change this type anymore!"); edc.Add(ed);

	ed.SetText(GetName(), "name"); ed.SetToolTipText("name of the GeomElement"); edc.Add(ed);
	ed.SetInt(GetElnum(), "reference_element_number"); ed.SetToolTipText("0 ... ground, otherwise insert number of existing element"); edc.Add(ed);

	ed.SetVector3D(col.X(),col.Y(),col.Z(),"RGB_color");	ed.SetToolTipText("[red, green, blue], range = 0..1");	edc.TreeAdd("Graphics",ed);
	//SetElemDataVector3D(edc, col, "RGB_color");	edc.Last().SetToolTipText("[red, green, blue], range = 0..1");
	ed.SetDouble(transparency, "transparency"); ed.SetToolTipText("transparency [0..1], 0=transparent, 1=solid, set -1 if global transparency is used"); edc.TreeAdd("Graphics",ed);
	
	ed.SetInt(drawstyle, "drawstyle"); ed.SetToolTipText("+1: draw outline, +2 fill area, +4 highlight points, +8 colored: outline"); edc.TreeAdd("Graphics",ed);
	ed.SetDouble(pointsize, "pointsize"); ed.SetToolTipText("size for highlighted points [m]"); edc.TreeAdd("Graphics",ed);
	ed.SetDouble(linethickness, "linethickness"); ed.SetToolTipText("thickness of lines [pts]"); edc.TreeAdd("Graphics",ed);

	//SetElemDataVector3D(edc, translated, "Geometry_translated"); edc.Last().SetLocked(1); edc.Last().SetToolTipText("Geometry has been translated probably in order that the center of mass=[0,0,0]");
}

int GeomElement::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = 1;

	mystr type_old = GetElementSpec();
	mystr type_new = edc.TreeGetString("geom_element_type");
	if(!type_new.Compare(type_old))
	{
		GetMBS()->UO().InstantMessageText("ERROR: You MUST NOT change the type of the GeomElement!");
		return 0;
	}
	GetElemDataText(mbs, edc, "name", GetName(), 1);
	GetElemDataInt(mbs, edc, "reference_element_number", elnum, 1);

	edc.TreeGetVector3D("Graphics.RGB_color",col.X(),col.Y(),col.Z());


	double val=transparency;
	GetElemDataDouble(mbs, edc, "transparency", val, 0);
	edc.TreeGetDouble("Graphics.transparency",val);
	transparency = (float)val;

	edc.TreeGetDouble("Graphics.drawstyle",drawstyle);
	edc.TreeGetDouble("Graphics.pointsize",pointsize);
	edc.TreeGetDouble("Graphics.linethickness",linethickness);

	return rv;
}

void GeomElement::SetGeomElementColor()
{
	if (GetElnum() != 0)
	{
		if (col.X() == -1 && transparency == -1)
		{
			GetMBS()->SetColor(GetElement().GetCol());
		}
		else if (col.X() != -1 && transparency == -1)
		{
			GetMBS()->SetColor(col);
		}
		else if (col.X() == -1 && transparency != -1)
		{
			Vector3D c = GetElement().GetCol();
			mbs->GetRC()->glColor4f((float)c[0],(float)c[1],(float)c[2], transparency);
		}
		else
		{
			mbs->GetRC()->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);
		}
	}
	else
	{
		if (col.X() != -1)
		{
			if (transparency == -1)
			{
				GetMBS()->SetColor(col);
			}
			else
			{
				mbs->GetRC()->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);
			}
		}
	}
}


void GeomElement::SetDrawPoints()
{
	if (elnum)
	{
		for (int i=1; i <= locpoints.Length(); i++)
		{
			ptd(i) = GetElement().GetPosD(locpoints(i));
		}
		//GetMBS()->UO() << "in\n";
	}	
	else if (objdata)
	{
		int r = (int)(GetMBS()->GetDrawTime()/GetMBS()->GetObjDataStepSize());
		//global_uo << "r=" << r << ", objdata=" << GetObjdata() << " " <<
		//	GetMBS()->GetObjData(r, 3*(objdata-1)+1) << " " << GetMBS()->GetObjData(r, 3*(objdata-1)+2) << "\n";

		for (int i=1; i <= locpoints.Length(); i++)
		{

			ptd(i).X() = GetMBS()->GetObjData(r, 3*(objdata-1)+1) + locpoints(i).X();
			ptd(i).Y() = GetMBS()->GetObjData(r, 3*(objdata-1)+2) + locpoints(i).Y();
			ptd(i).Z() = GetMBS()->GetObjData(r, 3*(objdata-1)+3) + locpoints(i).Z();
		}
	}
	else CopyPtD();
}

void GeomElement::SetComputePoints()
{
	if (elnum != NULL)
	{
		for (int i=1; i <= locpoints.Length(); i++)
		{
			pt(i) = GetElement().GetPos(locpoints(i));
		}
	}	
	else CopyPt();
}

//2D: +++++++++++++++++++++++++++++++++++++
void GeomElement::CopyPt2D()
{
	for (int i=1; i <= locpoints2D.Length(); i++)
	{
		pt2D(i) = locpoints2D(i);
	}
}

void GeomElement::CopyPtD2D()
{
	for (int i=1; i <= locpoints2D.Length(); i++)
	{
		ptd2D(i) = locpoints2D(i);
	}
}

void GeomElement::SetDrawPoints2D()
{
	if (elnum != 0)
	{
		for (int i=1; i <= locpoints2D.Length(); i++)
		{
			ptd2D(i) = GetElement().GetPos2DD(locpoints2D(i));
		}
	}	
	else CopyPtD2D();
}

void GeomElement::SetComputePoints2D()
{
	if (elnum)
	{
		for (int i=1; i <= locpoints2D.Length(); i++)
		{
			pt2D(i) = GetElement().GetPos2D(locpoints2D(i));
		}
	}	
	else CopyPt2D();
}


//actual (global) position:
Vector3D GeomElement::GetTPoint(int i) const 
{
	if (elnum == 0) return locpoints(i);
	else return GetElement().GetPos(locpoints(i));
}

Vector3D GeomElement::GetTPointD(int i) const 
{
	if (elnum == 0) return locpoints(i);
	else return GetElement().GetPosD(locpoints(i));
}

Vector2D GeomElement::GetTPoint2D(int i) const 
{
	if (elnum == 0) return locpoints2D(i);
	else 
	{
		if (nodenums.Length() >= i) return GetElement().GetNodePos2D(nodenums(i));
		else return GetElement().GetPos2D(locpoints2D(i));
	}
}

Vector2D GeomElement::GetTPointD2D(int i) const 
{
	if (elnum == 0) return locpoints2D(i);
	else return GetElement().GetPos2DD(locpoints2D(i));
}

Vector3D GeomElement::ToP3D(const Vector2D& p) const
{
	if (elnum == 0) return Vector3D(p.X(), p.Y(), 0);
	else return ((Body2D&)GetElement()).ToP3D(p);
}

Vector2D GeomElement::GetPos2D(const Vector2D& ploc) const
{
	if (elnum == 0) return ploc;
	else return GetElement().GetPos2D(ploc);
}

Vector2D GeomElement::GetVel2D(const Vector2D& ploc) const
{
	if (elnum == 0) return Vector2D(0,0);
	else return GetElement().GetVel2D(ploc);
}

Vector3D GeomElement::GetPos(const Vector3D& ploc) const
{
	if (elnum == 0) return ploc;
	else return GetElement().GetPos(ploc);
}

Vector3D GeomElement::GetVel(const Vector3D& ploc) const
{
	if (elnum == 0) return Vector3D(0.,0.,0.);
	else return GetElement().GetVel(ploc);
}

const Element& GeomElement::GetElement() const 
{
	return mbs->GetElement(elnum);
}

MBS* GeomElement::GetMBS() const 
{
	return mbs;
}





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GeomText2D::GeomText2D(MBS* mbsi, int elnumi, const char* textI, int rowI, int alignmentI, const Vector3D& coli):GeomElement()
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	row = rowI;
	alignment = alignmentI;
	text = new char[strlen(textI)+1];
	strcpy(text, textI);
}

void GeomText2D::DrawYourself()
{
	mbs->GetRC()->PrintTextStruct(row, alignment-1, text);
}








//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GeomLine2D::GeomLine2D(MBS* mbsi, int elnumi, const Vector2D& p1, const Vector2D& p2, const Vector3D& coli):GeomElement()
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	locpoints2D.Add(p1);
	locpoints2D.Add(p2);

	CopyPt2D();
	CopyPtD2D();
}

GeomLine2D::GeomLine2D(MBS* mbsi, int elnumi, int n1, int n2, const Vector3D& coli):GeomElement()
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	Vector2D p1, p2;
	if (elnumi != 0) 
	{
		p1 = GetElement().GetNodeLocPos2D(n1);
		p2 = GetElement().GetNodeLocPos2D(n2);
	}
	else GetMBS()->UO() << "ERROR: GeomLine2D: elnum==0 not possible for nodal coordinates\n";

	locpoints2D.Add(p1);
	locpoints2D.Add(p2);
	nodenums.Add(n1);
	nodenums.Add(n2);

	CopyPt2D();
	CopyPtD2D();
}

//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
double GeomLine2D::GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp)
{
	ind  = 1; //only one index possible in line
	Vector2D p1 = GetTPoint2D(1);
	Vector2D p2 = GetTPoint2D(2);
	double gap = RighthandMinDist2D(p1,p2,p,pp);
	return gap;
}

//get local position on element from global position and index (for polygon etc.)
Vector2D GeomLine2D::GetLocPos(int ind, const Vector2D& pglob) const
{
	Vector2D mp1 = GetTPoint2D(1);
	Vector2D mp2 = GetTPoint2D(2);

	//compute local position by linear interpolation
	double d1 = Dist(pglob, mp1);
	double d2 = Dist(pglob, mp2);
	double d = Dist(mp1,mp2);
	if (d < 1e-16) {d = 1; d1 = 0; d2 = 1;}

	return (1.-d1/d)*GetLocPoint2D(1)+(1.-d2/d)*GetLocPoint2D(2);
}

//get (deformed) normalized normal vector (outwards) at locpos
Vector2D GeomLine2D::GetNormal(int ind, const Vector2D& ploc) const
{
	Vector2D mp1 = GetTPoint2D(1);
	Vector2D mp2 = GetTPoint2D(2);

	Vector2D tt = mp2 - mp1; //tangential direction
	tt.Normalize();
	return Vector2D(tt.Y(),-tt.X()); //normal direction, outwards
}

void GeomLine2D::DrawYourself()
{
	//if (elnum == 0)
	{
		SetDrawPoints2D();

		Vector2D p1 = ptd2D(1);
		Vector2D p2 = ptd2D(2);

		Vector3D p1r = ToP3D(ptd2D(1));
		Vector3D p2r = ToP3D(ptd2D(2));

		//mbs->SetColor(col);
		//mbs->DrawZyl(p1r, p2r, 0.3*10*draw_dim.X(), 12);
		//mbs->MyDrawLine(p1r, p2r, draw_dim.X(), Vector3D(0.1,0.8,0.8));
		mbs->MyDrawLine(p1r, p2r, draw_dim.X(), col);

		if ( DrawStyle() & TDSArrowHeadLines ) //draw arrowhead
		{
// AD: changed .X() to .Y() and .Z() since drawdim.X() is supposed to be in [pts] but .Y() and .Z() supposed to be in [m]
			Vector2D v = p2-p1;
			v.Normalize();
			Vector2D n(v.Y(),-v.X());

			//Draw vectors ...
			double dd = draw_dim.Y();
			mbs->DrawZyl(p2r, ToP3D(p2-2*dd*v+0.5*dd*n), 0.4*draw_dim.Z(), 4);
			mbs->DrawZyl(p2r, ToP3D(p2-2*dd*v-0.5*dd*n), 0.4*draw_dim.Z(), 4);
		}
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


GeomPolygon2D::GeomPolygon2D(MBS* mbsi, int elnumi, const TArray<Vector2D>& p, const Vector3D& coli):GeomElement()
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	locpoints2D = p;

	CopyPt2D();
	CopyPtD2D();
}

////!AD:OLD
//void GeomPolygon2D::DrawYourself()
//{
//	SetDrawPoints2D();
//
//	ptd.SetLen(ptd2D.Length());
//
//	for (int i=1; i <= ptd2D.Length(); i++)
//		ptd(i) = ToP3D(ptd2D(i));
//
//	mbs->SetColor(col);
//	mbs->DrawPolygon(ptd,GetMBS()->GetIOption(110));
//}

//!AD:New-fancy
void GeomPolygon2D::DrawYourself()
{
	SetDrawPoints2D();

	ptd.SetLen(ptd2D.Length());

	for (int i=1; i <= ptd2D.Length(); i++)
		ptd(i) = ToP3D(ptd2D(i));

	mbs->SetColor(col);
	
	if( DrawStyle() & TDSFillAreas )
	{
		mbs->DrawPolygon(ptd,GetMBS()->GetIOption(110));
	}

	if( DrawStyle() & TDSDrawOutline )
	{
		mbs->DrawPolygonOutline(ptd,LineThickness());
	}

	if( DrawStyle() & TDSColoredLines )
	{
		Vector3D old_linecolor = mbs->ColLine();
		mbs->ColLine() = col;
		mbs->DrawPolygonOutline(ptd,LineThickness());
		mbs->ColLine() = old_linecolor;
	}

	if( DrawStyle() & TDSHighlightPoints )
	{
		for(int i=1; i<=ptd.Length(); i++)
		{
			mbs->DrawSphere(ptd(i),PointSize(),4);
		}
	}
}

int GeomPolygon2D::ReadFromFile(mystr& filename)
{
// reads file in format:
// optional:       % comments are skipped
// 1nd line:       point1.X   point1.Y ...            we assume 2 entries per line
// ...
	CMFile file_points(filename,TFMread);
	mystr filestring;
	if (!file_points.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename.c_str();
		return 0;
	}
	file_points.RWF(filestring); // entire file read to string

// parse block
	ptd2D.Flush();
	int rv;
	int pos=0;
	mystr line;
	int line_count = 0; // count data lines
	
	while (pos!=-1)
	{
		rv = filestring.GetUntil(pos,'\n',line,1);
		if(!line.Left(1).Compare(mystr("%")) && pos!=-1)
		{
			line_count++;
			int pos_line = 0;
			if(line.Length() == 1) pos_line=-1; // skip last line.... 

			while (pos_line!=-1)
			{
				mystr word;
				word = line.GetWord(pos_line,1);
				double x = word.MakeDouble();
				word = line.GetWord(pos_line,1);
				double y = word.MakeDouble();
				ptd2D.Add(Vector2D(x,y));
				break;
			}
		}
	}
	return ptd2D.Length();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//full circle:
GeomCircle2D::GeomCircle2D(MBS* mbsi, int elnumi, const Vector2D& p1, double rad, const Vector3D& coli):GeomElement()
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	locpoints2D.Add(p1);
	r = rad;

	CopyPt2D();
	CopyPtD2D();
}

//only segment:
GeomCircle2D::GeomCircle2D(MBS* mbsi, int elnumi, const Vector2D& p1, const Vector2D& pseg1, const Vector2D& pseg2, double rad, const Vector3D& coli)
{
	mbs = mbsi;
	elnum = elnumi;
	col = coli;

	locpoints2D.Add(p1);
	locpoints2D.Add(pseg1);
	locpoints2D.Add(pseg2);
	r = rad;

	CopyPt2D();
	CopyPtD2D();
}

//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
double GeomCircle2D::GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp)
{
	ind  = 1; //only one index possible in line
	double gap = 0;

	if (locpoints2D.Length() == 1)
	{
		Vector2D pc = GetTPoint2D(1);
		pp = p-pc;
		gap = pp.Norm();
		if (gap != 0) pp *= (r/gap);
		pp += pc;
		gap -= r;
	}
	else
	{
		GetMBS()->UO() << "ERROR: GeomCircle2D::not yet implemented!\n";
	}
	return gap;
}

//get local position on element from global position and index (for polygon etc.)
Vector2D GeomCircle2D::GetLocPos(int ind, const Vector2D& pglob) const
{
	if (GetElnum() == 0)
	{
		return pglob;
	}
	else
	{
		//only for rigid bodies!!!
		if (!GetElement().IsRigid()) GetMBS()->UO() << "ERROR: GeomCircle2D: only rigid bodies allowed!\n";

		Vector2D v = pglob - GetElement().GetRefPos2D();
		return ((Body2D&)GetElement()).GetRotMatrix2D().GetTp() * v;
	}
}

//get (deformed) normalized normal vector (outwards) at locpos
Vector2D GeomCircle2D::GetNormal(int ind, const Vector2D& ploc) const
{

	//Vector between circle point and position at circle ....
	Vector2D pc = GetTPoint2D(1);
	Vector2D p = GetPos2D(ploc);
	Vector2D v = p-pc;
	v.Normalize();
	return v;
}

void GeomCircle2D::DrawYourself()
{
	double ns = 12;
	if (draw_dim.Y() != 0) ns = draw_dim.Y();

	//draw points must be Vector3D
	ptd.SetLen((int)ns);
	Vector2D pc = GetTPointD2D(1);

	double phi;
	for (int i=0; i < ns; i++)
	{
		phi = (double)i/ns*2.*MY_PI;
		ptd(i+1) = ToP3D(pc + Vector2D(sin(phi)*r,cos(phi)*r));
	}

	////$ PG 2011-12-15:[ testing, colorplot fieldvariables of "mother" finite element
	//bool colormode = false;
	//FieldVariableDescriptor * fvdptr = GetMBS()->GetActualPostProcessingFieldVariable();
	//if (fvdptr != NULL)
	//{
	//	colormode = true;
	//}

	//if (0) //(colormode)
	//{
	//	Vector2D lp = locpoints2D(1);
	//	Mass2D* el2D = (Mass2D*) (&GetElement());
	//	mbs->SetColor(el2D->GetFieldVariableValue(*fvdptr, lp, true));
	//}
	//else
	//{
	//	mbs->SetColor(col);
	//}
	////$ PG 2011-12-15:]

	mbs->SetColor(col);

	mbs->DrawPolygon(ptd,1,0.5);

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
double GeomTrig3D::GetNearestPoint(const Vector3D& p, int& ind, Vector3D& pp)
{
	ind  = 1; //only one index possible in Trig!!!
	double gap;

	Vector3D p1 = GetTPoint(1);
	Vector3D p2 = GetTPoint(2);
	Vector3D p3 = GetTPoint(3);
	gap = MinDistTP(p1,p2,p3,p,pp);

	Vector3D n;
	Normal3D(p1,p2,p3,n);
	if (n * (p-pp) < 0) gap = -gap;

	return gap;
}

//get local position on element from global position and index (for polygon etc.)
Vector3D GeomTrig3D::GetLocPos(int ind, const Vector3D& pglob) const
{
	if (GetElnum() == 0)
	{
		return pglob;
	}
	else
	{
		//only for rigid bodies!!!
		if (!GetElement().IsRigid()) GetMBS()->UO() << "ERROR: GeomTrig3D: only rigid bodies allowed!\n";

		Vector3D v = pglob - GetElement().GetRefPos();
		return GetElement().GetRotMatrix().GetTp() * v;
	}
}

//get (deformed) normalized normal vector (outwards) at locpos
Vector3D GeomTrig3D::GetNormal(int ind, const Vector3D& ploc) const
{

	//Vector between circle point and position at circle ....
	Vector3D n;
	Normal3D(GetTPoint(1),GetTPoint(2),GetTPoint(3),n);
	return n;
}

Vector3D GeomTrig3D::GetTangent(int ind, const Vector3D& ploc) const
{
	Vector3D v = GetTPoint(2) - GetTPoint(1);
	v.Normalize();
	return v;
}

void GeomTrig3D::DrawYourself()
{

	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();

	if (transparency != 0)
	{
		SetGeomElementColor();
		pCurrentRC->glBeginTriangles();
		//pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

		Vector3D n;
		Normal3D(ptd(1),ptd(2),ptd(3),n);
		pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);

		pCurrentRC->glVertex((float)ptd(1).X(),(float)ptd(1).Y(),(float)ptd(1).Z());
		pCurrentRC->glVertex((float)ptd(2).X(),(float)ptd(2).Y(),(float)ptd(2).Z());
		pCurrentRC->glVertex((float)ptd(3).X(),(float)ptd(3).Y(),(float)ptd(3).Z());

		pCurrentRC->glEnd();
	}

	if (0)
	{
		mbs->ChooseColor(0.f,0.f,0.f);
		pCurrentRC->ChooseLineThickness(1);
		pCurrentRC->glBeginLines();
		pCurrentRC->glVertex((float)ptd(1).X(),(float)ptd(1).Y(),(float)ptd(1).Z());
		pCurrentRC->glVertex((float)ptd(2).X(),(float)ptd(2).Y(),(float)ptd(2).Z());
		pCurrentRC->glVertex((float)ptd(2).X(),(float)ptd(2).Y(),(float)ptd(2).Z());
		pCurrentRC->glVertex((float)ptd(3).X(),(float)ptd(3).Y(),(float)ptd(3).Z());
		pCurrentRC->glVertex((float)ptd(3).X(),(float)ptd(3).Y(),(float)ptd(3).Z());
		pCurrentRC->glVertex((float)ptd(1).X(),(float)ptd(1).Y(),(float)ptd(1).Z());
		pCurrentRC->glEnd();
	}

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GeomQuad3D::DrawYourself()
{

	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();


	SetGeomElementColor();
	//pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

	Vector3D n;
	Normal3D(ptd(1), ptd(2), ptd(4), n);
	pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);

	pCurrentRC->glBeginQuads();
	pCurrentRC->glVertex((float)ptd(4).X(),(float)ptd(4).Y(),(float)ptd(4).Z());
	pCurrentRC->glVertex((float)ptd(3).X(),(float)ptd(3).Y(),(float)ptd(3).Z());
	pCurrentRC->glVertex((float)ptd(2).X(),(float)ptd(2).Y(),(float)ptd(2).Z());
	pCurrentRC->glVertex((float)ptd(1).X(),(float)ptd(1).Y(),(float)ptd(1).Z());
	pCurrentRC->glEnd();

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomCube3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GeomElement::GetElementData(edc);

	ElementData ed;
	
	//depricated: use GeomOrthoCube if center point and size is needed
	//Box3D b;
	//for (int i=1; i <= 8; i++) b.Add(locpoints(i));
	//Vector3D center = b.Center();
	//Vector3D size = b.PMax()-b.PMin();

	//SetElemDataVector3D(edc, center, "center_point"); edc.Last().SetToolTipText("Center point in global coordinates");
	//SetElemDataVector3D(edc, size, "size"); edc.Last().SetToolTipText("Dimension of cube in X, Y and Z-direction");
	//ed.SetBool(0, "use_center_point_and_size"); ed.SetToolTipText("Use center point of cube and size to define cube."); edc.Add(ed);

	//#bottom: 1-2-4-3
	//#top:    5-6-8-7

	ed.SetVector3D(locpoints(1).X(),locpoints(1).Y(),locpoints(1).Z(), "point1"); ed.SetToolTipText("Bottom point 1 of bottom points: 1-2-4-3");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(2).X(),locpoints(2).Y(),locpoints(2).Z(), "point2"); ed.SetToolTipText("Bottom point 2 of bottom points: 1-2-4-3");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(3).X(),locpoints(3).Y(),locpoints(3).Z(), "point3"); ed.SetToolTipText("Bottom point 3 of bottom points: 1-2-4-3");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(4).X(),locpoints(4).Y(),locpoints(4).Z(), "point4"); ed.SetToolTipText("Bottom point 4 of bottom points: 1-2-4-3");	edc.TreeAdd("Geometry",ed);

	ed.SetVector3D(locpoints(5).X(),locpoints(5).Y(),locpoints(5).Z(), "point5"); ed.SetToolTipText("Bottom point 5 of bottom points: 5-6-8-7");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(6).X(),locpoints(6).Y(),locpoints(6).Z(), "point6"); ed.SetToolTipText("Bottom point 6 of bottom points: 5-6-8-7");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(7).X(),locpoints(7).Y(),locpoints(7).Z(), "point7"); ed.SetToolTipText("Bottom point 7 of bottom points: 5-6-8-7");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(8).X(),locpoints(8).Y(),locpoints(8).Z(), "point8"); ed.SetToolTipText("Bottom point 8 of bottom points: 5-6-8-7");	edc.TreeAdd("Geometry",ed);
}

int GeomCube3D::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = GeomElement::SetElementData(edc);

	//depricated: use GeomOrthoCube if center point and size is needed
	//int use_center_point_and_size = 0;
	//GetElemDataBool(mbs, edc, "use_center_point_and_size", use_center_point_and_size, 0);
	//
	//if (use_center_point_and_size)
	//{
	//	Vector3D center, size;
	//	GetElemDataVector3D(mbs, edc, "center_point", center, 1);
	//	GetElemDataVector3D(mbs, edc, "size", size, 1);

	//	Vector3D size2 = 0.5*size;
	//	SetLocPoints(center, size2);
	//}
	//else
	//{		
		GetElemDataVector3D(mbs, edc, "Geometry.point1", locpoints(1), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point2", locpoints(2), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point3", locpoints(3), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point4", locpoints(4), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point5", locpoints(5), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point6", locpoints(6), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point7", locpoints(7), 1);
		GetElemDataVector3D(mbs, edc, "Geometry.point8", locpoints(8), 1);
	//}

	return rv;
}

void GeomCube3D::Translate(const Vector3D& translate)
{
	translated += translate;
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) += translate;
	}
}

void GeomCube3D::Rotate(const Matrix3D& rotate)
{
  for (int i=1; i <= locpoints.Length(); i++)
	{
		Vector3D help(locpoints(i));
		locpoints(i) = rotate*help;
	}
}

void GeomCube3D::GetTrigs(TArray<int3>& trigs) const
{
	trigs.SetLen(0);

	const int ind[36]={
		  1,3,4,
			4,2,1,
			5,6,8,
			8,7,5,
			1,2,6,
			6,5,1,
			2,4,8,
			8,6,2,
			4,3,7,
			7,8,4,
			3,1,5,
			5,7,3};

	for (int face=0; face < 12; face++)
	{
		trigs.Add(int3(ind[face*3+0], ind[face*3+1], ind[face*3+2]));
	}	
}

double GeomCube3D::ComputeVolume() const
{
	TArray<int3> trigs;
	GetTrigs(trigs);
	//return ComputeVolumeTrigs(trigs, locpoints);
	//$ DR 2013-01-31, see log 382, computation moved to MultibodySystem

	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(1.0,"density"); edc.Add(ed);	// just a dummy value
	Matrix3D t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);
	return edc_rv.TreeGetDouble("double");
}

Vector3D GeomCube3D::ComputeCenterOfMass() const
{
	TArray<int3> trigs;
	GetTrigs(trigs);
	//return ComputeCenterOfMassTrigs(trigs, locpoints);
	//$ DR 2013-01-31, see log 382, computation moved to MultibodySystem

	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(1.0,"density"); edc.Add(ed);	// just a dummy value
	Matrix3D t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);
	return edc_rv.TreeGetDouble("center_of_mass");
}

Matrix3D GeomCube3D::ComputeMassMomentOfInertia(double rho) const
{
	TArray<int3> trigs;
	GetTrigs(trigs);
	//return ComputeMassMomentOfInertiaTrigs(rho, trigs, locpoints);

	//$ DR 2013-01-31, see log 382, computation moved to MultibodySystem
	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(rho,"density"); edc.Add(ed);	
	Matrix3D t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);

	double * IMp;
	int size = 3;
	edc_rv.TreeGetMatrix("moment_of_inertia",&IMp,size,size);
	Matrix IMmat(size,size,IMp);
	return IMmat;
}


void GeomCube3D::DrawYourself()
{
	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();

	const int ind[24]=
	{1,3,4,2,
	5,6,8,7,
	1,2,6,5,
	2,4,8,6,
	4,3,7,8,
	3,1,5,7
	};

	SetGeomElementColor();

	Vector3D n;
	if (GetMBS()->GetIOption(134)) //draw faces
	{
		for (int face=0; face < 6; face++)
		{

			pCurrentRC->glBeginQuads();

			//pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

			Vector3D& p1 = ptd(ind[face*4+0]);
			Vector3D& p2 = ptd(ind[face*4+1]);
			Vector3D& p3 = ptd(ind[face*4+2]);
			Vector3D& p4 = ptd(ind[face*4+3]);

			Normal3D(p1,p2,p4,n);
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p1.X(),(float)p1.Y(),(float)p1.Z());
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p2.X(),(float)p2.Y(),(float)p2.Z());
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p3.X(),(float)p3.Y(),(float)p3.Z());
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glVertex((float)p4.X(),(float)p4.Y(),(float)p4.Z());

			pCurrentRC->glEnd();
		}
	}

	if (GetMBS()->GetIOption(133)) //draw lines
	{
		double th = 1+mbs->GetDOption(115); //rigid body line thickness
		GetMBS()->MyDrawLine(ptd(1), ptd(3), th);
		GetMBS()->MyDrawLine(ptd(3), ptd(4), th);
		GetMBS()->MyDrawLine(ptd(4), ptd(2), th);
		GetMBS()->MyDrawLine(ptd(2), ptd(1), th);

		GetMBS()->MyDrawLine(ptd(5), ptd(6), th);
		GetMBS()->MyDrawLine(ptd(6), ptd(8), th);
		GetMBS()->MyDrawLine(ptd(8), ptd(7), th);
		GetMBS()->MyDrawLine(ptd(7), ptd(5), th);

		
		GetMBS()->MyDrawLine(ptd(1), ptd(5), th);
		GetMBS()->MyDrawLine(ptd(3), ptd(7), th);
		GetMBS()->MyDrawLine(ptd(4), ptd(8), th);
		GetMBS()->MyDrawLine(ptd(2), ptd(6), th);
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomOrthoCube3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GeomElement::GetElementData(edc);
	ElementData ed;
	ed.SetVector3D(locpoints(9).X(),locpoints(9).Y(),locpoints(9).Z(), "center_point"); ed.SetToolTipText("Center point in global coordinates");	edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(size.X(),size.Y(),size.Z(), "size"); ed.SetToolTipText("Dimension of cube in X, Y and Z-direction");	edc.TreeAdd("Geometry",ed);
	ed.SetMatrix(A.GetMatPtr(),A.Getrows(), A.Getcols(),"rotation_matrix"); ed.SetToolTipText("The rotation matrix defines the orientation of the cube (global_point = matrix . local_point)."); edc.TreeAdd("Geometry",ed);


}

int GeomOrthoCube3D::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = GeomElement::SetElementData(edc);
	GetElemDataVector3D(mbs, edc, "Geometry.center_point", locpoints(9), 1);
	GetElemDataVector3D(mbs, edc, "Geometry.size", size, 1);
	Matrix tmp;
	GetElemDataMatrix(mbs, edc, "Geometry.rotation_matrix", tmp, 1);
	A = Matrix3D(tmp(1,1),tmp(1,2),tmp(1,3), tmp(2,1), tmp(2,2), tmp(2,3), tmp(3,1), tmp(3,2), tmp(3,3)); 
	SetLocPoints(locpoints(9), size, A);
	return rv;
}

void GeomOrthoCube3D::Rotate(const Matrix3D& rotate)
{
  for (int i=1; i <= 8; i++) // same as GeomOrthoCube3D, Length()-1 is used, rotation relative to center point
	{
		Vector3D help(locpoints(i));
		locpoints(i) = rotate*help;
	}
}
void GeomOrthoCube3D::SetLocPoints(const Vector3D center_point, const Vector3D sizei, const Matrix3D A)
{
	Vector3D size2 = 0.5*sizei;
	GeomCube3D::SetLocPoints(Vector3D(0.), size2);
	Rotate(A);
	Translate(center_point);
	size = sizei;
	locpoints(9) = center_point;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomZyl3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GeomElement::GetElementData(edc);

	ElementData ed;

	ed.SetInt(tile, "draw_resolution"); ed.SetToolTipText("Number of quadrangles to draw the cylinder surface"); edc.TreeAdd("Graphics",ed);
	ed.SetDouble(rz, "radius"); ed.SetToolTipText("radius of the cylinder"); edc.TreeAdd("Geometry",ed);
	ed.SetDouble(ri, "radius_hole"); ed.SetToolTipText("inner radius of the cylinder (0 if full cylinder)"); edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(1).X(),locpoints(1).Y(),locpoints(1).Z(),"axis_point1"); ed.SetToolTipText("point on axis of rotation"); edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(2).X(),locpoints(2).Y(),locpoints(2).Z(),"axis_point2"); ed.SetToolTipText("point on axis of rotation"); edc.TreeAdd("Geometry",ed);

}

int GeomZyl3D::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = GeomElement::SetElementData(edc);
	tile = 16;

	GetElemDataInt(mbs, edc, "Graphics.draw_resolution", tile, 0);
	GetElemDataDouble(mbs, edc, "Geometry.radius", rz, 1);
	GetElemDataDouble(mbs, edc, "Geometry.radius_hole", ri, 1);
	GetElemDataVector3D(mbs, edc, "Geometry.axis_point1", locpoints(1), 1);
	GetElemDataVector3D(mbs, edc, "Geometry.axis_point2", locpoints(2), 1);


	return rv;
}

void GeomZyl3D::Translate(const Vector3D& translate)
{
	translated += translate;
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) += translate;
	}
}

double GeomZyl3D::ComputeVolume() const
{
	return (Sqr(rz)-Sqr(ri))*MY_PI*Dist(locpoints(1), locpoints(2));
}

Vector3D GeomZyl3D::ComputeCenterOfMass() const
{
	return 0.5*(locpoints(1) + locpoints(2));
}

Matrix3D GeomZyl3D::ComputeMassMomentOfInertia(double rho) const
{
	Vector3D n, t1, t2;
	n = locpoints(1) - locpoints(2);
	double length = n.Norm();

	n.Normalize();
	n.SetNormalBasis(t1,t2);

	double m = rho*ComputeVolume();

	Matrix3D A;
	A.Set(n, t1, t2);
	Matrix3D Iphi(0.5*m*(Sqr(rz)+Sqr(ri)), 1./12.*m*(3*(Sqr(rz)+Sqr(ri))+Sqr(length)), 1./12.*m*(3*(Sqr(rz)+Sqr(ri))+Sqr(length))); //MomentofInertia about rotation axis (component 1) and other axis (component 2 and 3)

	return A*Iphi*A.GetTp();
}


void GeomZyl3D::DrawYourself()
{
	RenderContext* pCurrentRC = mbs->GetRC();

	if (!GetMBS()->GetIOption(134) && !GetMBS()->GetIOption(133) && GetElnum() != 0) return;

	SetDrawPoints();

	double tn = 1+mbs->GetDOption(115); //rigid body line thickness

	Vector3D p1 = ptd(1);
	Vector3D p2 = ptd(2);
	Vector3D n;
	double r = rz;//; 

	double i, phi, dphi;
	double dtile = tile;
	if (dtile == 0) {dtile = 1;}

	Vector3D p3=p1+Vector3D(1,1,1);
	Vector3D v21=p2-p1;
	Vector3D v=v21;
	v.Normalize();
	double np1=v*p1;

	if (v.X() != 0) {p3.X() = 1./v.X()*(np1-p3.Y()*v.Y()-p3.Z()*v.Z());}
	else if (v.Y() != 0) {p3.Y() = 1./v.Y()*(np1-p3.X()*v.X()-p3.Z()*v.Z());}
	else if (v.Z() != 0) {p3.Z() = 1./v.Z()*(np1-p3.Y()*v.Y()-p3.X()*v.X());}
	else {return;}

	Vector3D n1 = p3-p1;
	n1.Normalize();
	Vector3D n2 = v.Cross(n1);
	n2.Normalize();

	for (i = 0; i < dtile; i++)
	{
		SetGeomElementColor();

		phi = i*2.*MY_PI/dtile;
		dphi = 2.*MY_PI/dtile;

		Vector3D pz1=r*cos(phi)*n1+r*sin(phi)*n2+p1;
		Vector3D pz2=pz1+v21;
		Vector3D pz4=r*cos(phi+dphi)*n1+r*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+v21;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		if (GetMBS()->GetIOption(134) || GetElnum() == 0) //draw faces
		{
			pCurrentRC->glBeginQuads();
			pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
			pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
			pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
			pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
			pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
			pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
			pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
			pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
			pCurrentRC->glEnd();

			if (ri == 0.) 
			{//full cylinder
				pCurrentRC->glBeginTriangles();
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);

				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
				pCurrentRC->glEnd();
			}
			else
			{ // hollow cylinder
				Vector3D pz1i=ri*cos(phi)*n1+ri*sin(phi)*n2+p1;
				Vector3D pz2i=pz1i+v21;
				Vector3D pz4i=ri*cos(phi+dphi)*n1+ri*sin(phi+dphi)*n2+p1;
				Vector3D pz3i=pz4i+v21;
				Vector3D nz1i = pz1i-p1;
				Vector3D nz2i = pz4i-p1;
				nz1i.Normalize();
				nz2i.Normalize();

				//inner area
				pCurrentRC->glBeginQuads();
				pCurrentRC->glNormal(-1.f*(float)nz1i[0],-1.f*(float)nz1i[1],-1.f*(float)nz1i[2]);
				pCurrentRC->glVertex((float)pz1i[0],(float)pz1i[1],(float)pz1i[2]);
				pCurrentRC->glNormal(-1.f*(float)nz2i[0],-1.f*(float)nz2i[1],-1.f*(float)nz2i[2]);
				pCurrentRC->glVertex((float)pz4i[0],(float)pz4i[1],(float)pz4i[2]);
				pCurrentRC->glNormal(-1.f*(float)nz2i[0],-1.f*(float)nz2i[1],-1.f*(float)nz2i[2]);
				pCurrentRC->glVertex((float)pz3i[0],(float)pz3i[1],(float)pz3i[2]);
				pCurrentRC->glNormal(-1.f*(float)nz1i[0],-1.f*(float)nz1i[1],-1.f*(float)nz1i[2]);
				pCurrentRC->glVertex((float)pz2i[0],(float)pz2i[1],(float)pz2i[2]);
				pCurrentRC->glEnd();

				//sides

				pCurrentRC->glBeginQuads();
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz1i[0],(float)pz1i[1],(float)pz1i[2]);
				pCurrentRC->glNormal(-1.f*(float)v[0],-1.f*(float)v[1],-1.f*(float)v[2]);
				pCurrentRC->glVertex((float)pz4i[0],(float)pz4i[1],(float)pz4i[2]);

				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz3i[0],(float)pz3i[1],(float)pz3i[2]);
				pCurrentRC->glNormal((float)v[0],(float)v[1],(float)v[2]);
				pCurrentRC->glVertex((float)pz2i[0],(float)pz2i[1],(float)pz2i[2]);
				pCurrentRC->glEnd();

				if (GetMBS()->GetIOption(133) || (GetElnum() == 0 && mbs->GetIOption(134)==0)) //draw lines
				{
					GetMBS()->MyDrawLine(pz1i, pz4i, tn);
					GetMBS()->MyDrawLine(pz2i, pz3i, tn);
				}
			}

		}

		if (GetMBS()->GetIOption(133) || (GetElnum() == 0 && mbs->GetIOption(134)==0)) //draw lines
		{
			GetMBS()->MyDrawLine(pz1, pz4, tn);
			GetMBS()->MyDrawLine(pz2, pz3, tn);

			if (i==0 || i==tile/4 || i==tile/2 || i==3*tile/4)
			{
				//GetMBS()->MyDrawLine(pz1, pz2, tn); //JG 2013-01-18, because lines can rotate around axis of cylinder, bad graphics. ...
			}
		}
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GeomPipe3D::DrawYourself()
{

	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();

	Vector3D p1 = ptd(1);
	Vector3D p2 = ptd(2);
	Vector3D va = p2-p1;
	Vector3D n1 = ptd(3) - p1;

	Vector3D n0 = va;
	n0.Normalize();
	Vector3D n2 = n0.Cross(n1);

	double dtile = tile;

	pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

	double phitot = phimax-phimin;
	double phi, dphi;


	for (int i = 0; i <= tile-1; i++)
	{
		phi = i*phitot/dtile+phimin;
		dphi = phitot/dtile;

		Vector3D pz1=ro*cos(phi)*n1+ro*sin(phi)*n2+p1;
		Vector3D pz2=pz1+va;
		Vector3D pz4=ro*cos(phi+dphi)*n1+ro*sin(phi+dphi)*n2+p1;
		Vector3D pz3=pz4+va;
		Vector3D nz1 = pz1-p1;
		Vector3D nz2 = pz4-p1;
		nz1.Normalize();
		nz2.Normalize();

		pCurrentRC->glBeginQuads();
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
		pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
		pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
		pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
		pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
		pCurrentRC->glEnd();

		Vector3D pz1i=ri*cos(phi)*n1+ri*sin(phi)*n2+p1;
		Vector3D pz2i=pz1i+va;
		Vector3D pz4i=ri*cos(phi+dphi)*n1+ri*sin(phi+dphi)*n2+p1;
		Vector3D pz3i=pz4i+va;

		pCurrentRC->glBeginQuads();
		pCurrentRC->glNormal(-(float)nz1[0],-(float)nz1[1],-(float)nz1[2]);
		pCurrentRC->glVertex((float)pz1i[0],(float)pz1i[1],(float)pz1i[2]);
		pCurrentRC->glNormal(-(float)nz2[0],-(float)nz2[1],-(float)nz2[2]);
		pCurrentRC->glVertex((float)pz4i[0],(float)pz4i[1],(float)pz4i[2]);
		pCurrentRC->glNormal(-(float)nz2[0],-(float)nz2[1],-(float)nz2[2]);
		pCurrentRC->glVertex((float)pz3i[0],(float)pz3i[1],(float)pz3i[2]);
		pCurrentRC->glNormal(-(float)nz1[0],-(float)nz1[1],-(float)nz1[2]);
		pCurrentRC->glVertex((float)pz2i[0],(float)pz2i[1],(float)pz2i[2]);
		pCurrentRC->glEnd();

		if (1)
		{
			//side faces
			Vector3D n = -1*va;
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glBeginQuads();
			pCurrentRC->glVertex((float)pz1[0], (float)pz1[1], (float)pz1[2]);
			pCurrentRC->glVertex((float)pz4[0], (float)pz4[1], (float)pz4[2]);
			pCurrentRC->glVertex((float)pz4i[0],(float)pz4i[1],(float)pz4i[2]);
			pCurrentRC->glVertex((float)pz1i[0],(float)pz1i[1],(float)pz1i[2]);
			pCurrentRC->glEnd();

			n = va;
			pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
			pCurrentRC->glBeginQuads();
			pCurrentRC->glVertex((float)pz2[0], (float)pz2[1], (float)pz2[2]);
			pCurrentRC->glVertex((float)pz3[0], (float)pz3[1], (float)pz3[2]);
			pCurrentRC->glVertex((float)pz3i[0],(float)pz3i[1],(float)pz3i[2]);
			pCurrentRC->glVertex((float)pz2i[0],(float)pz2i[1],(float)pz2i[2]);
			pCurrentRC->glEnd();
		}


		if (phitot < 2.*MY_PI && drawendface)
		{
			if (i == 0)
			{
				Vector3D n;
				Normal3D(pz1,pz2,pz1i,n);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glBeginQuads();
				pCurrentRC->glVertex((float)pz1[0], (float)pz1[1], (float)pz1[2]);
				pCurrentRC->glVertex((float)pz2[0], (float)pz2[1], (float)pz2[2]);
				pCurrentRC->glVertex((float)pz2i[0],(float)pz2i[1],(float)pz2i[2]);
				pCurrentRC->glVertex((float)pz1i[0],(float)pz1i[1],(float)pz1i[2]);
				pCurrentRC->glEnd();
			}
			if (i == tile-1)
			{
				Vector3D n;
				Normal3D(pz3,pz4,pz3i,n);
				pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
				pCurrentRC->glBeginQuads();
				pCurrentRC->glVertex((float)pz4[0], (float)pz4[1], (float)pz4[2]);
				pCurrentRC->glVertex((float)pz3[0], (float)pz3[1], (float)pz3[2]);
				pCurrentRC->glVertex((float)pz3i[0],(float)pz3i[1],(float)pz3i[2]);
				pCurrentRC->glVertex((float)pz4i[0],(float)pz4i[1],(float)pz4i[2]);
				pCurrentRC->glEnd();
			}
		}
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

GeomRotObject3D::GeomRotObject3D(MBS* mbsi, const Vector3D& p0, TArray<Vector2D> segs, 
																 TArray<Vector3D> cols, const Matrix3D& rotmat, int tilei):GeomElement(mbsi)
{
	elnum = 0;
	mbs = mbsi;
	tile = tilei;
	if (tile < 2) tile = 2;

	segments = segs;
	colorseg = cols;
	if (colorseg.Length() < segments.Length())
	{
		Vector3D c = colorseg(1);
		colorseg.SetLen(segments.Length());
		for (int i=1; i <= segments.Length(); i++)
			colorseg(i) = c;
	}

	Matrix3D A=rotmat;

	int np = tile*segs.Length();
	locpoints.SetLen(np);

	for (int j = 1; j <= segs.Length(); j++)
	{
		for (int i = 1; i <= tile; i++)
		{
			double td = tile;
			double nd = np;
			double id = i;

			double phi1 = (id-1.)/td*2.*MY_PI;
			//double phi2 = (id   )/td*2.*MY_PI;
			Vector2D p1 = segments(j);
			//Vector2D p2 = segments(j+1);

			Vector3D p;
			p.X() = sin(phi1)*p1.Y();
			p.Y() = cos(phi1)*p1.Y();
			p.Z() = p1.X();
			locpoints(i+(j-1)*tile) = A*p + p0;
		}
	}

	ptd.SetLen(np);
	pt.SetLen(np);

	CopyPt();
	CopyPtD();
}

void GeomRotObject3D::CopyFrom(const GeomElement& e)
{
	GeomElement::CopyFrom(e);
	const GeomRotObject3D& ce = (const GeomRotObject3D&)e;

	tile =ce.tile;
	rz = ce.rz;
	colorseg = ce.colorseg;
	segments = ce.segments;

}

void GeomRotObject3D::DrawYourself()
{
	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();


	Vector3D p1, p2, p3, p4;
	Vector3D p1n, p2n, p3n, p4n; //next quad
	Vector3D p1l, p2l, p3l, p4l; //last quad
	Vector3D n, nnext, nlast;

	int np = segments.Length();

	for (int j = 1; j < np; j++)
	{
		for (int i = 1; i <= tile; i++)
		{
			int inext = i%tile+1;
			int inext2= (i+1)%tile+1;
			int ilast = i-1; if (ilast < 1) ilast = tile;

			p1 = ptd(i+(j-1)*tile);
			p2 = ptd(i+(j-0)*tile);
			p3 = ptd(inext+(j-0)*tile);
			p4 = ptd(inext+(j-1)*tile);

			p1n = ptd(inext +(j-1)*tile);
			p2n = ptd(inext +(j-0)*tile);
			p3n = ptd(inext2+(j-0)*tile);
			p4n = ptd(inext2+(j-1)*tile);

			p1l = ptd(ilast +(j-1)*tile);
			p2l = ptd(ilast +(j-0)*tile);
			p3l = ptd(i+(j-0)*tile);
			p4l = ptd(i+(j-1)*tile);

			Normal3D(p1, p2, p3, n);
			Normal3D(p1l, p2l, p3l, nlast);
			Normal3D(p1n, p2n, p3n, nnext);

			Vector3D n1 = 0.5*(nlast + n);
			Vector3D n2 = 0.5*(nnext + n);

			col = colorseg(j);

			//if (i == 1) col *= 0.5;
			pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

			pCurrentRC->glBeginQuads();
			pCurrentRC->glNormal((float)n1[0],(float)n1[1],(float)n1[2]);
			pCurrentRC->glVertex((float)p1[0],(float)p1[1],(float)p1[2]);
			pCurrentRC->glNormal((float)n2[0],(float)n2[1],(float)n2[2]);
			pCurrentRC->glVertex((float)p4[0],(float)p4[1],(float)p4[2]);
			pCurrentRC->glNormal((float)n2[0],(float)n2[1],(float)n2[2]);
			pCurrentRC->glVertex((float)p3[0],(float)p3[1],(float)p3[2]);
			pCurrentRC->glNormal((float)n1[0],(float)n1[1],(float)n1[2]);
			pCurrentRC->glVertex((float)p2[0],(float)p2[1],(float)p2[2]);
			pCurrentRC->glEnd();


		}
	}
}







//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GeomSphere3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GeomElement::GetElementData(edc);

	ElementData ed;

	ed.SetInt(draw_dim.Y(), "draw_resolution"); ed.SetToolTipText("Number of quadrangles to draw the sphere"); edc.TreeAdd("Graphics",ed);
	ed.SetDouble(rz, "radius"); ed.SetToolTipText("radius of the sphere"); edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(locpoints(1).X(),locpoints(1).Y(),locpoints(1).Z(),"center_point"); ed.SetToolTipText("center point of the sphere"); edc.TreeAdd("Geometry",ed);

}

int GeomSphere3D::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = GeomElement::SetElementData(edc);
	draw_dim.Y() = 16;

	GetElemDataDouble(mbs, edc, "Graphics.draw_resolution", draw_dim.Y(), 0);
	GetElemDataDouble(mbs, edc, "Geometry.radius", rz, 1);
	GetElemDataVector3D(mbs, edc, "Geometry.center_point", locpoints(1), 1);

	return rv;
}

void GeomSphere3D::Translate(const Vector3D& translate)
{
	translated += translate;
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) += translate;
	}
}

double GeomSphere3D::ComputeVolume() const
{
	return 4./3.*Cub(rz)*MY_PI;
}

Vector3D GeomSphere3D::ComputeCenterOfMass() const
{
	return locpoints(1);
}

Matrix3D GeomSphere3D::ComputeMassMomentOfInertia(double rho) const
{
	//See Shabana, Computational Dynamics 1994, p. 381
	double m = rho*ComputeVolume();

	return Matrix3D(2./5.*m*Sqr(rz));
}


//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
double GeomSphere3D::GetNearestPoint(const Vector3D& p, int& ind, Vector3D& pp)
{
	ind  = 1; //only one index possible in Sphere
	double gap = 0;

	Vector3D pc = GetTPoint(1);		//center of sphere, global; ACCELERATE: getrefpos directly?
	pp = p-pc;										//vector from center pc of sphere to p
	gap = pp.Norm();							//distance of pc and p
	if (gap != 0) pp *= (rz/gap);	//pp=relative vector from center to surface of sphere
	pp += pc;											//pp=global projected point at surface of sphere
	gap -= rz;										//distance between p and projected point pp

	return gap;
}

//get local position on element from global position and index (for polygon etc.)
Vector3D GeomSphere3D::GetLocPos(int ind, const Vector3D& pglob) const
{
	if (GetElnum() == 0)
	{
		return pglob;
	}
	else
	{
		//only for rigid bodies!!!
		if (!GetElement().IsRigid()) GetMBS()->UO() << "ERROR: GeomTrig3D: only rigid bodies allowed!\n";

		Vector3D v = pglob - GetElement().GetRefPos();
		//return GetBody3D().GetRotMatrix().GetTp() * v;
		return v; //only for Mass3D !!!
	}
}

//get (deformed) normalized normal vector (outwards) at locpos
Vector3D GeomSphere3D::GetNormal(int ind, const Vector3D& ploc) const
{
	return Vector3D(1.,1.,1.); //not used!!!
}

//get (deformed) normalized unique tangent vector at locpos
Vector3D GeomSphere3D::GetTangent(int ind, const Vector3D& ploc) const
{
	Vector3D n = GetNormal(ind, ploc);
	Vector3D t1,t2;
	n.SetNormalBasis(t1, t2);
	return t1;
}


void GeomSphere3D::DrawYourself()
{
	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();

	Vector3D pk = ptd(1);

	double i, j, phi, dphi;
	if (draw_dim.Y() == 0) {draw_dim.Y() = 1;}

	SetGeomElementColor();


	if (mbs->GetIOption(134)) //draw faces
	{
		double fill = 1;
		for (j = 0; j < draw_dim.Y()*fill; j++)
		{
			for (i = 0; i < draw_dim.Y(); i++)
			{
				phi = i*2.*MY_PI/draw_dim.Y();
				dphi = 2.*MY_PI/draw_dim.Y();
				double phik1 = j*MY_PI/draw_dim.Y();
				double phik2 = (j+1.)*MY_PI/draw_dim.Y();
				double xpos1 = rz*cos(phik1);
				double xpos2 = rz*cos(phik2);
				double ypos1 = rz*sin(phik1);
				double ypos2 = rz*sin(phik2);

				Vector3D n1(0.,1.,0.);
				Vector3D n2(0.,0.,1.);
				Vector3D vr1(xpos1,0,0);
				Vector3D vr2(xpos2,0,0);

				Vector3D pz1=ypos1*cos(phi)*n1+ypos1*sin(phi)*n2+pk+vr1;
				Vector3D pz2=ypos2*cos(phi)*n1+ypos2*sin(phi)*n2+pk+vr2;
				Vector3D pz4=ypos1*cos(phi+dphi)*n1+ypos1*sin(phi+dphi)*n2+pk+vr1;
				Vector3D pz3=ypos2*cos(phi+dphi)*n1+ypos2*sin(phi+dphi)*n2+pk+vr2;
				Vector3D nz1 = pz1-pk;
				Vector3D nz2 = pz2-pk;
				Vector3D nz3 = pz3-pk;
				Vector3D nz4 = pz4-pk;
				nz1.Normalize();
				nz2.Normalize();
				nz3.Normalize();
				nz4.Normalize();

				pCurrentRC->glBeginQuads();
				pCurrentRC->glNormal((float)nz1[0],(float)nz1[1],(float)nz1[2]);
				pCurrentRC->glVertex((float)pz1[0],(float)pz1[1],(float)pz1[2]);
				pCurrentRC->glNormal((float)nz4[0],(float)nz4[1],(float)nz4[2]);
				pCurrentRC->glVertex((float)pz4[0],(float)pz4[1],(float)pz4[2]);
				pCurrentRC->glNormal((float)nz3[0],(float)nz3[1],(float)nz3[2]);
				pCurrentRC->glVertex((float)pz3[0],(float)pz3[1],(float)pz3[2]);
				pCurrentRC->glNormal((float)nz2[0],(float)nz2[1],(float)nz2[2]);
				pCurrentRC->glVertex((float)pz2[0],(float)pz2[1],(float)pz2[2]);
				pCurrentRC->glEnd();

			}
		}
	}

	if (GetMBS()->GetIOption(133) || (GetElnum() == 0 && GetMBS()->GetIOption(134)==0))
	{
		GetMBS()->SetColor(colgrey1);
		GetMBS()->GetRC()->glBeginPoints();
		GetMBS()->GetRC()->ChoosePointSize(GetMBS()->GetIOption(117));
		GetMBS()->GetRC()->glVertex(pk.GetVecPtr());
		GetMBS()->GetRC()->glEnd();
	}

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomShip3D::DrawYourself()
{
	//Center is center of mass in element = (0,0,0)
	RenderContext* pCurrentRC = mbs->GetRC();

	float transp = transparency;

	Vector3D pref(0.,0.,0.);
	Matrix3D A;
	A.SetDiag(1);

	if (GetElnum())
	{
		pref = GetElement().GetRefPosD();
		A = GetElement().GetRotMatrixD();
	}

	double l = size.X();
	double b = size.Y();
	double h = size.Z();
	double lbug = l*0.2;

	//Center in shape is rear bottom center
	TArray<Vector3D> pshape;

	pshape.Add(Vector3D(0.,0.,0.));
	pshape.Add(Vector3D(0.,0.5*b,0.));
	pshape.Add(Vector3D(l-lbug,0.5*b,0.));

	//bug:
	for (int i=0; i <= tile-1; i++)
	{
		double x = (double)i/(double)tile;
		pshape.Add(Vector3D(l-lbug+lbug*x,0.5*b*(1.-Sqr(x)),0.));
	}

	pshape.Add(Vector3D(l,0.,0.));

	//double and mirror shape points with respect to y-coordinate:
	int nshape = pshape.Length();
	for (int i=1; i <= nshape; i++)
	{
		Vector3D v = pshape(nshape-i+1);
		v.Y() *= -1;
		pshape.Add(v);
	}

	int tileh = tile;
	TArray<double> hshape;
	for (int i=1; i <= tileh; i++)
	{
		double x = (double)(i-1.)/(double)(tileh-1.);
		x *= x;
		double v = (cos(x*MY_PI*0.5+MY_PI*0.1)-cos(MY_PI*0.6))/(cos(MY_PI*0.1)-cos(MY_PI*0.6));
		hshape.Add(v);
	}


	Vector3D p[5];
	pCurrentRC->glColor4f(0.6f,0.6f,0.8f, transp);
	for (int i=1; i <= nshape-1; i++)
	{
		p[1] = pshape(i);
		p[2] = pshape(i+1);
		p[3] = pshape(2*nshape-i);
		p[4] = pshape(2*nshape-i+1);

		for (int j = 1; j <= 4; j++)
		{
			p[j].Z() += h;
		}

		for (int j = 1; j <= 4; j++)
		{
			p[j] = pref+A*(p[j]-com);
		}
		mbs->DrawQuad(p[1], p[2], p[3], p[4]);
	}

	pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transp);

	double xfact = 0.06;
	for (int i=1; i <= 2*nshape; i++)
	{
		for (int k=1; k <= hshape.Length()-1; k++)
		{
			double z1 = h*(1.-(double)(k-1)/(double)(hshape.Length()-1));
			double z2 = h*(1.-(double)(k  )/(double)(hshape.Length()-1));

			p[1] = pshape(i+1); //top
			p[2] = pshape(i);   //top
			p[3] = pshape(i);		//bottom
			p[4] = pshape(i+1); //bottom

			for (int j = 1; j <= 2; j++)
			{
				p[j].X() -= xfact*(1-hshape(k))*p[j].X();
				p[j].Y() *= hshape(k);
			}
			for (int j = 3; j <= 4; j++)
			{
				p[j].X() -= xfact*(1-hshape(k+1))*p[j].X();
				p[j].Y() *= hshape(k+1);
			}

			p[1].Z() += z1;
			p[2].Z() += z1;
			p[3].Z() += z2;
			p[4].Z() += z2;

			for (int j = 1; j <= 4; j++)
			{
				p[j] = pref+A*(p[j]-com);
			}
			mbs->DrawQuad(p[4], p[3], p[2], p[1]);
			//mbs->DrawQuad(p[1], p[2], p[3], p[4]);
		}
	}

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Vector3D GeomGround::FStreet(double x, int mode, double sx, double sy, double sz) const
{
	Vector3D f;
	double fz;
	if (x < 0.8) fz = -(1-Sqr((x-0.8)/0.8))*sz;
	else fz = -sz;

	if (mode == 1)
	{
		if (x < 0.25) f = Vector3D(sx*x,0,fz);
		else if (x < 0.75) 
		{
			//double fy = (-1.+3.*x-1./4.*Cub(4*x-2));
			double xx = 2*(x-0.25);
			double fy = 10*Cub(xx)-15*To4(xx)+6*To5(xx); //f' and f'' 0

			f = Vector3D(sx*x, sy*fy,fz);
		}
		else f = Vector3D(sx*x, sy ,fz);
	}
	if (mode == 2)
	{
		if (x < 0.2) f = Vector3D(sx*x,0,fz);
		else if (x < 1) 
		{
			//double fy = (-1.+3.*x-1./4.*Cub(4*x-2));
			double xx = (x-0.2)/0.8;
			double fy = 1-cos((xx*2*MY_PI)*2); //f' and f'' 0

			f = Vector3D(sx*x, sy*fy,fz);
		}
		else f = Vector3D(sx*x, 0 ,fz);
	}
	return f;
}


void GeomGround::SetStreet()
{

	double sizex=160; //dimension of street midline
	double sizey=2;
	double sizez=-0.1;  //height over ground
	double streetwidth = 1; //street size
	double streetdec = 30;
	double groundwidth = 8;
	double resx = 120;
	double resyg = 16; //resolution of one side of ground
	double resys = 2;  //resolution of street
	double ry = resys+2*resyg;
	int mode = 2;
	Vector3D streetcol(0.6,0.52,0.);
	//Vector3D streetcol(0.7,0.7,0.7);
	Vector3D groundcol1(0.1,0.7,0.);
	Vector3D groundcol2(0.6,0.5,0.);
	Vector3D quadcol(0.7,0.7,0.7);

	//generate street line:
	TArray<Vector3D> streetline;
	for (double i = 0; i <= resx+1; i++)
	{
		double x = i/resx;
		streetline.Add(FStreet(x,mode,sizex, sizey, streetdec));
	}

	//generate points:
	for (double i = 1; i <= resx+1; i++)
	{
		double x = (i-1)/resx;
		double fz;
		if (x < 0.8) fz = -(1-Sqr((x-0.8)/0.8))*streetdec;
		else fz = -streetdec;

		for (double j = 0; j <= resyg-1; j++)
		{
			Vector3D v = streetline((int)i+1)-streetline((int)i);
			Vector3D n(-v.Y(),v.X(),0);
			n.Normalize();
			double y = -0.5*streetwidth;
			Vector3D side(i/resx*sizex,-groundwidth,fz);

			Vector3D p((resyg-j)/resyg*side+j/resyg*(y*n+streetline((int)i)));
			//p.Z() = -(1-Sqr(j/resyg))*sizez;
			//p.Z() = -(1-(j/resyg))*sizez;
			p.Z() = -Sqr(1-(j/resyg))*sizez+fz;
			locpoints.Add(p);
		}
		for (double j = 0; j <= resys-1; j++)
		{
			Vector3D v = streetline((int)i+1)-streetline((int)i);
			Vector3D n(-v.Y(),v.X(),0);
			n.Normalize();
			double y = 0.5*streetwidth*(j/resys*2-1);
			locpoints.Add(streetline((int)i)+y*n);
		}
		for (double j = 0; j <= resyg; j++)
		{
			Vector3D v = streetline((int)i+1)-streetline((int)i);
			Vector3D n(-v.Y(),v.X(),0);
			n.Normalize();
			double y = +0.5*streetwidth;
			Vector3D side(i/resx*sizex,sizey+groundwidth,-sizez+fz);

			Vector3D p(j/resyg*side+(resyg-j)/resyg*(y*n+streetline((int)i)));
			//p.Z() = -(1-Sqr((resyg-j)/resyg))*sizez;
			//p.Z() = -(1-((resyg-j)/resyg))*sizez;
			p.Z() = -Sqr(1-((resyg-j)/resyg))*sizez+fz;
			locpoints.Add(p);
		}
	}

	//generate quads:
	for (double i = 0; i < resx; i++)
	{
		for (double j = 1; j <= ry; j++)
		{
			quads.Add(int4(((int)ry+1)*(int)i+(int)j,((int)ry+1)*(int)i+(int)j+1,
				((int)ry+1)*((int)i+1)+(int)j+1,((int)ry+1)*((int)i+1)+(int)j));
			if (j <= resyg)
				cols.Add(j/resyg*groundcol2 + (resyg-j)/resyg*groundcol1);
			else if (j <= resyg+resys)
				cols.Add(streetcol);
			else
				cols.Add((j-resys-resyg)/resyg*groundcol1 + (ry-j)/resyg*groundcol2);
		}
	}

	//add object:
	{
		double objy = groundwidth*0.5;
		double objx = groundwidth*0.1;
		double x = 0.16;
		double xoff = sizex*x;
		double objz = 2;
		double zoff = -(1-Sqr((x-0.8)/0.8))*streetdec;

		int i = locpoints.Length();
		int qi = quads.Length();
		locpoints.Add(Vector3D(-objx+xoff,-0.5*objy,0*objz+zoff));
		locpoints.Add(Vector3D( objx+xoff,-0.5*objy,0*objz+zoff));
		locpoints.Add(Vector3D(-objx+xoff, 0.5*objy,0*objz+zoff));
		locpoints.Add(Vector3D( objx+xoff, 0.5*objy,0*objz+zoff));
		locpoints.Add(Vector3D(-objx+xoff,-0.5*objy,1*objz+zoff));
		locpoints.Add(Vector3D( objx+xoff,-0.5*objy,1*objz+zoff));
		locpoints.Add(Vector3D(-objx+xoff, 0.5*objy,1*objz+zoff));
		locpoints.Add(Vector3D( objx+xoff, 0.5*objy,1*objz+zoff));
		quads.Add(int4(i+1,i+3,i+4,i+2));
		quads.Add(int4(i+5,i+6,i+8,i+7));
		quads.Add(int4(i+1,i+2,i+6,i+5));
		quads.Add(int4(i+2,i+4,i+8,i+6));
		quads.Add(int4(i+4,i+3,i+7,i+8));
		quads.Add(int4(i+3,i+1,i+5,i+7));
		for (int j = qi+1; j <= qi+6; j++)
		{
			cols.Add(quadcol);
			quads(j).Invert();
		}
	}

	//Model is ready ....
	CopyPt();
	CopyPtD();

	//build bounding boxes:
	Box3D box;
	boxes.Flush();
	for (int i=1; i <= quads.Length(); i++)
	{
		Box3D b;
		b.Add(ptd(quads(i).Get(1)));
		b.Add(ptd(quads(i).Get(2)));
		b.Add(ptd(quads(i).Get(3)));
		b.Add(ptd(quads(i).Get(4)));
		boxes.Add(b);
		box.Add(b);
	}

	tree = SearchTree(50,50,1, box);
	for (int i=1; i <= quads.Length(); i++)
	{
		Box3D b;
		b.Add(ptd(quads(i).Get(1)));
		b.Add(ptd(quads(i).Get(2)));
		b.Add(ptd(quads(i).Get(3)));
		b.Add(ptd(quads(i).Get(4)));

		tree.AddItem(b,i);
	}

	mbs->UO() << "generated " << quads.Length() << " quads\n";


}

void GeomGround::SetFlat(const Vector3D& dim, int resx, int resy)
{

	double sizex = dim.X();
	double sizey = dim.Y();
	double sizez = dim.Z();  //height over ground

	for (double i = 1; i <= resx+1; i++)
	{
		double x = (i-1.)/(double)resx*sizex;

		for (double j = 1; j <= resy+1; j++)
		{
			double y = (j-1.)/(double)resy*sizey;
			locpoints.Add(Vector3D(x-sizex*0.5,y-sizey*0.5,sizez));
		}
	}

	//generate quads:
	for (int i = 0; i < resx; i++)
	{
		for (int j = 1; j <= resy; j++)
		{
			quads.Add(int4((resy+1)*i+j,(resy+1)*i+j+1, (resy+1)*(i+1)+j+1,(resy+1)*(i+1)+j));
			double a = (double)i/resx;
			double b = (double)j/resy;
			cols.Add(Vector3D(0.4+0.2*(1-Sqr(a)),0.3+0.2*b,0.1));
		}
	}


	//Model is ready ....
	CopyPt();
	CopyPtD();

	//build bounding boxes:
	Box3D box;
	boxes.Flush();
	for (int i=1; i <= quads.Length(); i++)
	{
		Box3D b;
		b.Add(ptd(quads(i).Get(1)));
		b.Add(ptd(quads(i).Get(2)));
		b.Add(ptd(quads(i).Get(3)));
		b.Add(ptd(quads(i).Get(4)));
		boxes.Add(b);
		box.Add(b);
	}

	tree = SearchTree(20,20,1, box);
	for (int i=1; i <= quads.Length(); i++)
	{
		Box3D b;
		b.Add(ptd(quads(i).Get(1)));
		b.Add(ptd(quads(i).Get(2)));
		b.Add(ptd(quads(i).Get(3)));
		b.Add(ptd(quads(i).Get(4)));

		tree.AddItem(b,i);
	}

	mbs->UO() << "generated " << quads.Length() << " quads\n";


}


void GeomGround::DrawYourself()
{
	RenderContext* pCurrentRC = mbs->GetRC();

	int drawlines = 1;
	Vector3D n;
	for (int i=1; i <= quads.Length(); i++)
	{
		pCurrentRC->glBeginQuads();

		pCurrentRC->glColor4f((float)cols(i)[0],(float)cols(i)[1],(float)cols(i)[2], transparency);
		//pCurrentRC->glColor3f(0.4f,0.4f,0.05f);
		Normal3D(ptd(quads(i).Get(4)),ptd(quads(i).Get(2)),ptd(quads(i).Get(1)),n);
		pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);

		pCurrentRC->glVertex((float)ptd(quads(i).Get(4)).X(),(float)ptd(quads(i).Get(4)).Y(),(float)ptd(quads(i).Get(4)).Z());
		pCurrentRC->glVertex((float)ptd(quads(i).Get(3)).X(),(float)ptd(quads(i).Get(3)).Y(),(float)ptd(quads(i).Get(3)).Z());
		pCurrentRC->glVertex((float)ptd(quads(i).Get(2)).X(),(float)ptd(quads(i).Get(2)).Y(),(float)ptd(quads(i).Get(2)).Z());
		pCurrentRC->glVertex((float)ptd(quads(i).Get(1)).X(),(float)ptd(quads(i).Get(1)).Y(),(float)ptd(quads(i).Get(1)).Z());
		pCurrentRC->glEnd();

		if (drawlines)
		{
			pCurrentRC->ChooseLineThickness(1);
			pCurrentRC->glColor3f(0.2f,0.2f,0.2f);
			pCurrentRC->glBeginLineStrip();
			pCurrentRC->glVertex((float)ptd(quads(i).Get(4)).X(),(float)ptd(quads(i).Get(4)).Y(),(float)ptd(quads(i).Get(4)).Z());
			pCurrentRC->glVertex((float)ptd(quads(i).Get(3)).X(),(float)ptd(quads(i).Get(3)).Y(),(float)ptd(quads(i).Get(3)).Z());
			pCurrentRC->glVertex((float)ptd(quads(i).Get(2)).X(),(float)ptd(quads(i).Get(2)).Y(),(float)ptd(quads(i).Get(2)).Z());
			pCurrentRC->glVertex((float)ptd(quads(i).Get(1)).X(),(float)ptd(quads(i).Get(1)).Y(),(float)ptd(quads(i).Get(1)).Z());
			pCurrentRC->glVertex((float)ptd(quads(i).Get(4)).X(),(float)ptd(quads(i).Get(4)).Y(),(float)ptd(quads(i).Get(4)).Z());
			pCurrentRC->glEnd();
		}

	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++           GEOMMesh                 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void GeomMesh3D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	GeomElement::GetElementData(edc);

	ElementData ed;

	ed.SetBool(smooth_drawing, "smooth_drawing"); ed.SetToolTipText("Draw smooth interpolation of surface"); edc.TreeAdd("Graphics",ed);

	double angle = drawedge_angle;
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) angle *= 180./MY_PI;
	ed.SetDouble(angle, "draw_edge_angle"); ed.SetToolTipText((mystr("Minimum angle between two triangles that defines an edge ")+GetRotUnitStr(GetMBS()->GetIOption(120))).c_str()); edc.TreeAdd("Graphics",ed);
	//ed.SetDouble(stl_tolerance, "STL_tolerance"); ed.SetToolTipText("Minimum angle between two triangles that defines an edge"); edc.Add(ed);

	Vector3D resize(1.,1.,1.);
	Vector3D rotate(0.,0.,0.);
	Vector3D translate(0.,0.,0.);
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) rotate *= 180./MY_PI; //not necessary, but for completeness!

	ed.SetVector3D(resize.X(), resize.Y(), resize.Z(),"transform_scale"); ed.SetToolTipText("Resize GeomElement in X, Y and Z direction [sX, sY, sZ]"); edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(rotate.X(), rotate.Y(), rotate.Z(),"transform_rotation"); ed.SetToolTipText("Resize GeomElement in X, Y and Z direction [sX, sY, sZ]"); edc.TreeAdd("Geometry",ed);
	ed.SetVector3D(translate.X(), translate.Y(), translate.Z(),"transform_position"); ed.SetToolTipText("Translate GeomElement in X, Y and Z direction [tX, tY, tZ]"); edc.TreeAdd("Geometry",ed);



	int rows = trigs.Length();
	int cols = 3;
	if (rows == 0) cols = 0;
	Matrix m(rows, cols);
	for (int i = 1; i <= rows; i++)
	{
		m(i,1) = trigs(i).Get(1);
		m(i,2) = trigs(i).Get(2);
		m(i,3) = trigs(i).Get(3);
	}
	ed.SetMatrix(m.GetMatPtr(),rows, cols, "triangles"); ed.SetToolTipText("Fill in point numbers of each triangle: p1, p2, p3; p4, p5, p6 ..."); edc.TreeAdd("MeshData",ed);



	rows = locpoints.Length();
	if (rows == 0) cols = 0;
	else cols = 3;
	m = Matrix(rows, cols);
	for (int i = 1; i <= rows; i++)
	{
		m(i,1) = locpoints(i).X();
		m(i,2) = locpoints(i).Y();
		m(i,3) = locpoints(i).Z();
	}

	ed.SetMatrix(m.GetMatPtr(),rows, cols, "points"); ed.SetToolTipText("Fill in point coordinates: X1, Y1, Z1; X2, Y2, Z2 ..."); edc.TreeAdd("MeshData",ed);

}

int GeomMesh3D::SetElementData(const ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = GeomElement::SetElementData(edc);

	smooth_drawing = 1;
	drawedge_angle = MY_PI;

	GetElemDataBool(mbs, edc, "Graphics.smooth_drawing", smooth_drawing, 0);
	GetElemDataDouble(mbs, edc, "Graphics.draw_edge_angle", drawedge_angle, 0);
	if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) drawedge_angle *= MY_PI/180.;

	Matrix m;
	GetElemDataMatrix(mbs, edc, "MeshData.points", m, 1);

	locpoints.SetLen(m.Getrows());

	for (int i = 1; i <= locpoints.Length(); i++)
	{
		locpoints(i).X() = m(i,1);
		locpoints(i).Y() = m(i,2);
		locpoints(i).Z() = m(i,3);
	}

	//GetMBS()->UO() << "n-points=" << locpoints.Length() << "\n";

	GetElemDataMatrix(mbs, edc, "MeshData.triangles", m, 1);

	if (locpoints.Length() == 0) trigs.SetLen(0);
	else
	{
		trigs.SetLen(m.Getrows());

		for (int i = 1; i <= trigs.Length(); i++)
		{
			for (int j = 1; j <= 3; j++)
			{
				trigs(i).Get(j) = 1;
				if ((int)m(i,j) < 1) GetMBS()->UO() << "Error in GeomElement: Triangle point index < 1! Indices start with 1!\n";
				else if ((int)m(i,j) > locpoints.Length()) GetMBS()->UO() << "Error in GeomElement: Triangle point index " << (int)m(i,j) << " > maximum (" << locpoints.Length() << ")!\n";
				else trigs(i).Get(j) = (int)m(i,j);
			}
		}
		//GetMBS()->UO() << "n-trigs=" << trigs.Length() << "\n";


		Vector3D resize(1.,1.,1.);
		Vector3D rotate(0.,0.,0.);
		Vector3D translate(0.,0.,0.);

		GetElemDataVector3D(mbs, edc, "Geometry.transform_scale", resize, 0);
		GetElemDataVector3D(mbs, edc, "Geometry.transform_rotation", rotate, 0); 
    	GetElemDataVector3D(mbs, edc, "Geometry.transform_position", translate, 0); 

		if (GetMBS()->GetIOption(120) && !GetMBS()->IsLoadSaveMode()) rotate *= MY_PI/180.;

		if (!(resize.X() == 1. && resize.Y() == 1. && resize.Z() == 1. && rotate.Norm() == 0. && translate.Norm() == 0.))
		{
			Matrix3D rot = RotMatrix3(rotate.Z())*RotMatrix2(rotate.Y())*RotMatrix1(rotate.X());
			SetGeomModification(resize, rot, translate);
		}

		FinishMesh(); //needs to be called to compute edges, smooth normals, etc.!
	}


	return rv;
}


void GeomMesh3D::Init()
{
	drawedge_angle = MY_PI*2./10.;
	smooth_drawing = 1;
	//stl_tolerance = 1e-4;

	locpoints.SetLen(0);
	ptd.SetLen(0);
	pt.SetLen(0);
}

void GeomMesh3D::CopyFrom(const GeomElement& e)
{
	GeomElement::CopyFrom(e);
	const GeomMesh3D& ce = (const GeomMesh3D&)e;

	normals = ce.normals;
	drawnormals = ce.drawnormals;
	trigs = ce.trigs;

	points_to_trigs.SetLen(ce.points_to_trigs.Length());
	for (int i=1; i <= ce.points_to_trigs.Length(); i++)
	{
		points_to_trigs(i) = new IVector(*ce.points_to_trigs(i));
	}

	edges = ce.edges;

	points_to_edges.SetLen(ce.points_to_edges.Length());
	for (int i=1; i <= ce.points_to_edges.Length(); i++)
	{
		points_to_edges(i) = new IVector(*ce.points_to_edges(i));
	}

	smooth_drawing = ce.smooth_drawing;
	drawedge_angle = ce.drawedge_angle;
	//stl_tolerance	= ce.stl_tolerance;
}

GeomMesh3D::~GeomMesh3D()
{
	for (int i=1; i <= points_to_trigs.Length(); i++)
	{
		delete points_to_trigs(i);
	}
	for (int i=1; i <= points_to_edges.Length(); i++)
	{
		delete points_to_edges(i);
	}
}

int GeomMesh3D::NTrigs() const 
{
	return trigs.Length();
}

void GeomMesh3D::GetTrig0(int i, Vector3D& p1, Vector3D& p2, Vector3D& p3, Vector3D& normal)
{
	p1 = locpoints(trigs(i).Get(1));
	p2 = locpoints(trigs(i).Get(2));
	p3 = locpoints(trigs(i).Get(3));
	if (smooth_drawing)
	{
		normal += normals(i*3-2);
		normal += normals(i*3-1);
		normal += normals(i*3);
		normal.Normalize();
	}
	else
	{
		normal = normals(i);
	}
}

void GeomMesh3D::ClearMesh() 
{
	locpoints.SetLen(0);
	ptd.SetLen(0);
	pt.SetLen(0);

	trigs.SetLen(0);
	normals.SetLen(0);
	drawnormals.SetLen(0);

	for (int i=1; i <= points_to_trigs.Length(); i++)
	{
		delete points_to_trigs(i);
	}
	for (int i=1; i <= points_to_edges.Length(); i++)
	{
		delete points_to_edges(i);
	}
	points_to_edges.SetLen(0);
	edges.SetLen(0);

	Init();
}
/*
int GeomMesh3D::AddPoint(const Vector3D& p1, double eps)
{
	double eps2 = Sqr(eps);
	for (int i=1; i <= locpoints.Length(); i++)
	{
		if (Dist2(p1, locpoints(i)) <= eps2)
		{
			return i;
			}
			}
			return locpoints.Add(p1);
			}*/

int GeomMesh3D::AddTrig(int n1, int n2, int n3, const Vector3D& normal)
{
	if (smooth_drawing)
	{
		normals.Add(normal);
		normals.Add(normal);
		normals.Add(normal);
	}
	else
	{
		normals.Add(normal);
	}

	return trigs.Add(int3(n1,n2,n3));
}

void GeomMesh3D::ReadSTLMesh(char* file, int binary)
{
	//$ 2013-01-14 reading from file moved to File2EDC
	ElementDataContainer edc;
	mbs->File2EDC(file,&edc);

	Matrix m;
	GetElemDataMatrix(mbs,edc,"triangles",m);
	for(int i=1; i<=m.Getrows();i++)
	{
		trigs.Add(int3(m(i,1),m(i,2),m(i,3)));
	}

	GetElemDataMatrix(mbs,edc,"points",m);
	for(int i=1; i<=m.Getrows();i++)
	{
		locpoints.Add(Vector3D(m(i,1),m(i,2),m(i,3)));
	}

	GetElemDataMatrix(mbs,edc,"normals",m);
	for(int i=1; i<=m.Getrows();i++)
	{
		normals.Add(Vector3D(m(i,1),m(i,2),m(i,3)));
		if (smooth_drawing)
		{
			normals.Add(Vector3D(m(i,1),m(i,2),m(i,3)));
			normals.Add(Vector3D(m(i,1),m(i,2),m(i,3)));
		}
	}
}

//void GeomMesh3D::ReadFEMesh3D(FEMesh3D* femesh3d)
void GeomMesh3D::ReadFromPointSequence(Vector3D* v3dseq, int n)
{
  TArray<Vector3D>* v3darr = new TArray<Vector3D>;
	for(int i=1; i<=n; i++)
	{
		v3darr->Add(v3dseq[i-1]);
	}
	ReadFromPointArray(v3darr);
	delete v3darr;
}

void GeomMesh3D::ReadFromDoubleSequence(double* dseq, int n)
{
  TArray<Vector3D>* v3darr = new TArray<Vector3D>;
	for(int i=1; i<=n; i+=3)
	{
		Vector3D v3d(dseq[i-1],dseq[i],dseq[i+1]);
		v3darr->Add(v3d);
	}
	ReadFromPointArray(v3darr);
	delete v3darr;
}

void GeomMesh3D::ReadFromPointArray(TArray<Vector3D>* v3darr)
{
	int n = v3darr->Length();

// create empty searchtree structure with bounding box for all points, 
// create mapping to reduced point list in pnums, 
// add unique points to searchtree and GeomMesh::locpoints
	IVector pnums; 
	pnums.SetLen(n); 
	Box3D box(Vector3D(0.0,0.0,0.0),Vector3D(0.0,0.0,0.0));
	for (int i=1; i <= n ; i++) 
		box.Add(v3darr->Get(i)); // compute bounding box only
	
	SearchTree st(20,20,20, box); // is initialized empty !
	IVector items;
	double tol = 1E-10;

	for (int i=1; i <= n; i++)
	{
		Box3D b(v3darr->Get(i), tol);
		st.GetItemsInBox(b, items);
		int found = 0;
		for (int j=1; j <= items.Length(); j++)
		{
			if (Dist2(v3darr->Get(i), locpoints(items(j))) <= tol) 
			{
				found = items(j);
				break;
			}
		}
		if (found)
		{
			pnums(i) = found;
		}
		else
		{
			pnums(i) = locpoints.Add(v3darr->Get(i));
			st.AddItem(Box3D(v3darr->Get(i),v3darr->Get(i)), pnums(i));
		}
	}

	for (int i=1; i <= (n/3); i++)
	{
		Vector3D normal;
		Normal3D(v3darr->Get(3*i-2), v3darr->Get(3*i-1), v3darr->Get(3*i),normal);
		AddTrig(pnums(3*i-2), pnums(3*i-1), pnums(3*i), normal);
	}
	mbs->UO(UO_LVL_ext) << mystr("GeomMesh3D found ") + mystr(locpoints.Length()) + mystr(" different points in ") + mystr(trigs.Length()) + mystr(" triangles \n");
}

//get number of triangle neighbouring at 'side'
//side = 1,2,3; side 1= node 1 and 2, side 2=node 2 and 3, ...
int GeomMesh3D::GetNeighbourTrig(int trig, int side, int& ntrig, int& nside) 
{
	int n1 = trigs(trig).Get(side);
	int n2 = trigs(trig).GetMod(side + 1);

	int2 edge1(n1,n2); 
	edge1.Sort(); //Need to sort, two triangles have same edge with reverted indices!!!!!

	const IVector& ptrigs = *points_to_trigs(n1);

	for (int i=1; i <= ptrigs.Length(); i++)
	{
		if (ptrigs(i) != trig)
		{
			int3 t = trigs(ptrigs(i));
			for (int j = 1; j <= 3; j++)
			{
				int2 edge2(t.GetMod(j), t.GetMod(j+1));
				edge2.Sort();

				if (edge1 == edge2)
				{
					ntrig = ptrigs(i);
					nside = j;
					return 1;
				}
			}
		}
	}
	ntrig = 0;
	nside = 0;

	return 0;
}

Vector3D GeomMesh3D::GeometryNormal(int trig1) const
{
	const Vector3D& p1 = locpoints(trigs(trig1).Get(1));
	const Vector3D& p2 = locpoints(trigs(trig1).Get(2));
	const Vector3D& p3 = locpoints(trigs(trig1).Get(3));

	Vector3D n = (p2-p1).Cross(p3-p1);
	n.Normalize();
	return n;
}

double GeomMesh3D::GetAngle(int trig1, int trig2) const
{
	return VectorAngle(GeometryNormal(trig1), GeometryNormal(trig2));
}

double GeomMesh3D::GetArea(int trig) const
{
	const Vector3D& p1 = locpoints(trigs(trig).Get(1));
	const Vector3D& p2 = locpoints(trigs(trig).Get(2));
	const Vector3D& p3 = locpoints(trigs(trig).Get(3));

	return Area(p1, p2, p3);
}

int GeomMesh3D::IsEdge(int2 edge) const //check if two point indices represent an edge, return edge number or zero
{
	edge.Sort();
	const IVector& iv = *points_to_edges(edge.Get(1));
	for (int j = 1; j <= iv.Length(); j++)
	{
		if (edge == edges(iv(j))) return iv(j);
	}
	return 0;
}

void GeomMesh3D::FinishMesh() //call after all triangles have been added
{
	//compute points_to_trigs
	points_to_trigs.SetLen(locpoints.Length());
	for (int i=1; i <= points_to_trigs.Length(); i++)
	{
		points_to_trigs(i) = new IVector(8); //well above average, not too many resizes
		points_to_trigs(i)->SetLen(0);
	}
	for (int i=1; i <= trigs.Length(); i++)
	{
		for (int j=1; j <= 3; j++)
		{
			points_to_trigs(trigs(i).Get(j))->Add(i);
		}
	}

	points_to_edges.SetLen(locpoints.Length());
	for (int i=1; i <= points_to_trigs.Length(); i++)
	{
		points_to_edges(i) = new IVector(2); //edges are seldom
		points_to_edges(i)->SetLen(0);
	}


	if (smooth_drawing) //each triangle has 3 own normals!!!
	{
		normals.SetLen(trigs.Length()*3);
		for (int i=1; i <= trigs.Length(); i++)
		{
			normals(i*3-2) = GeometryNormal(i);
			normals(i*3-1) = GeometryNormal(i);
			normals(i*3  ) = GeometryNormal(i);
		}
	}
	else //one normal per triangle
	{
		normals.SetLen(trigs.Length());
		for (int i=1; i <= normals.Length(); i++)
		{
			normals(i) = GeometryNormal(i);
		}
	}

	//GetMBS()->UO() << "number trigs=" << trigs.Length() << "\n";
	//GetMBS()->UO() << "number points=" << locpoints.Length() << "\n";

	edges.SetLen(0);

	//find neighbours, compute edges and compute smooth normals
	int nt, nside;
	int ninconsistent = 0;
	for (int i=1; i <= trigs.Length(); i++)
	{
		for (int j = 1; j <= 3; j++)
		{
			if (GetNeighbourTrig(i, j, nt, nside))
			{
				//GetMBS()->UO() << "Trig " << i << ": neighbor = " << nt << "\n";

				if (nt > i) //do not add edges twice!!!
				{
					double angle = GetAngle(nt, i);
					//GetMBS()->UO() << "Trig " << i << ": angle = " << angle << "\n";
					if (fabs(angle) > drawedge_angle)
					{
						int2 edge(trigs(i).GetMod(j),trigs(i).GetMod(j+1));
						edge.Sort();
						int edgenum = edges.Add(edge);
						points_to_edges(edge.Get(1))->Add(edgenum);
						points_to_edges(edge.Get(2))->Add(edgenum);

//						GetMBS()->UO() << "Edge (" << edge.Get(1) << "," << edge.Get(2) << "\n";

						//i=trigs.Length();
					}
				}
			}
			else
			{
				ninconsistent++;
			}
		}
	}

	//GetMBS()->UO() << "Edges added:" << edges.Length() << "\n";

	if (ninconsistent)
	{
		GetMBS()->UO() << "Error in STL mesh: " << ninconsistent << " triangles found with no neighbours!!!\n";
		return;
	}

	if (smooth_drawing)
	{
		IVector average_trigs;

		for (int i=1; i <= points_to_trigs.Length(); i++)
		{
			//ckeck point number i
			//check all triangles at point i:
			const IVector& iv = *points_to_trigs(i);
			const IVector& ive = *points_to_edges(i);

			if (ive.Length() != 0)
			{
				for (int j = 1; j <= ive.Length(); j++)
				{
					average_trigs.SetLen(0);
					int2 edge = edges(ive.Get(j));
					if (edge.Get(2) == i) edge.Swap();

					//find starting triangle with edge:
					int k=1;
					int found = 0;
					while (k <= iv.Length() && !found)
					{
						int p1 = trigs(iv(k)).Find(i);

						if (trigs(iv(k)).GetMod(p1+1) == edge.Get(2)) found = iv(k);
						else k++;
					}

					if (found)
					{
						average_trigs.Add(found);
						found = 0;
						//find further triangles until next edge

						int seg = trigs(average_trigs.Last()).Find(i) + 2;
						seg = (seg-1)%3 + 1;

						while (!found)
						{
							int ntrig, nseg;
							GetNeighbourTrig(average_trigs.Last(), seg, //trig, seg
								ntrig, nseg);

							average_trigs.Add(ntrig);
							seg = (nseg+2-1)%3 + 1;

							int2 isedge(trigs(ntrig).Get(nseg), trigs(ntrig).GetMod(nseg+2));
							if (IsEdge(isedge)) found = 1;
						}
					}
					else
					{
						ninconsistent++;
						/*
						GetMBS()->UO() << "point " << i << ": edge (" << edge.Get(1) << ", " << edge.Get(2) << ")\n";
						for (int ii=1; ii <= iv.Length(); ii++)
						{
							GetMBS()->UO() << "  trig " << iv(ii) << ":" << trigs(iv(ii)).Get(1) << ", "
								 << trigs(iv(ii)).Get(2) << ", " << trigs(iv(ii)).Get(3) << "\n";
						}
						Vector3D asdf(0.,0.,0.);

						*/
					}

					Vector3D anormal(0.,0.,0.);
					for (int k = 1; k <= average_trigs.Length(); k++)
					{
						int ti = average_trigs(k);
						double A = GetArea(ti);
						Vector3D normal = GeometryNormal(ti);
						anormal += A*normal;
					}
					anormal.Normalize();
					for (int k = 1; k <= average_trigs.Length(); k++)
					{ 
						int ti = average_trigs(k);
						int locp = trigs(ti).Find(i);
						normals(ti*3-3+locp) = anormal;
					}
				}
			}
			else
			{
				Vector3D anormal(0.,0.,0.);
				for (int k = 1; k <= iv.Length(); k++)
				{
					int ti = iv(k);
					int locp = trigs(ti).Find(i);
					double A = GetArea(ti);
					/*
					A = 1;
					double dist = Dist(locpoints(trigs(ti).GetMod(locp)),locpoints(trigs(ti).GetMod(locp+1)))
									+ Dist(locpoints(trigs(ti).GetMod(locp)),locpoints(trigs(ti).GetMod(locp+2)));
					if (dist != 0) A = 1./dist;
					*/

					Vector3D normal = GeometryNormal(ti);
					anormal += A*normal;
				}
				anormal.Normalize();
				for (int k = 1; k <= iv.Length(); k++)
				{ 
					int ti = iv(k);
					int locp = trigs(ti).Find(i);
					normals(ti*3-3+locp) = anormal;
				}
			}
		}
	}

	if (ninconsistent)
	{
		GetMBS()->UO() << "Error in STL mesh: " << ninconsistent << " triangles found with inconsistent neighbours!!!\n";
		return;
	}
}

void GeomMesh3D::SetGeomModification(const Vector3D& resize, const Matrix3D& rotate, const Vector3D& translate)
{
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i).X() *= resize.X();
		locpoints(i).Y() *= resize.Y();
		locpoints(i).Z() *= resize.Z();

		locpoints(i) = translate + rotate*locpoints(i);
	}
  for (int i=1; i <= normals.Length(); i++)
	{
		normals(i) = rotate*normals(i);
	}
}

void GeomMesh3D::Translate(const Vector3D& translate)
{
	translated += translate;
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) += translate;
	}
}

void GeomMesh3D::Rotation(const Matrix3D &rotate)
{
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) = rotate*locpoints(i);
	}
  for (int i=1; i <= normals.Length(); i++)
	{
		normals(i) = rotate*normals(i);
	}

}

void GeomMesh3D::Stretch(const double factor)
{
  for (int i=1; i <= locpoints.Length(); i++)
	{
		locpoints(i) *= factor;
	}
}

double GeomMesh3D::ComputeVolume() const
{
	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(1.0,"density"); edc.Add(ed);	// just a dummy value
	Matrix t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);
	return edc_rv.TreeGetDouble("volume");

	//$ DR 2013-01-31, old code, see log 382, computation moved to MultibodySystem
	////see Maple file "maple/integration_triangle_volume.mws" and paper 'computing moments of piecewise polynomial surfaces'
	//double vol = 0;
	//for (int i=1; i <= trigs.Length(); i++)
	//{
	//	const Vector3D& p1 = locpoints(trigs(i).Get(1));
	//	const Vector3D& p2 = locpoints(trigs(i).Get(2));
	//	const Vector3D& p3 = locpoints(trigs(i).Get(3));

	//	double v1x = p2.X()-p1.X();
	//	double v1y = p2.Y()-p1.Y();
	//	double v1z = p2.Z()-p1.Z();

	//	double v2x = p3.X()-p1.X();
	//	double v2y = p3.Y()-p1.Y();
	//	double v2z = p3.Z()-p1.Z();
	//	double pz = p1.Z();

	//	vol += 1./6.*v1z*(v1x*v2y-v2x*v1y)+1./6.*v2z*(v1x*v2y-v2x*v1y)+1./2.*pz*(v1x*v2y-v2x*v1y);
	//}

	//return vol;
}

Vector3D GeomMesh3D::ComputeCenterOfMass() const
{
	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(1.0,"density"); edc.Add(ed);	// just a dummy value
	Matrix t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);
	Vector3D com;
	edc_rv.TreeGetVector3D("center_of_mass", com.X(), com.Y(), com.Z());
	return com;

	//$ DR 2013-01-31, old code, see log 382, computation moved to MultibodySystem
	////see Maple file "maple/integration_triangle_volume.mws" and paper 'computing moments of piecewise polynomial surfaces'
	//Vector3D center(0.,0.,0.);
	//double vol = ComputeVolume();

	//for (int i=1; i <= trigs.Length(); i++)
	//{
	//	const Vector3D& p1 = locpoints(trigs(i).Get(1));
	//	const Vector3D& p2 = locpoints(trigs(i).Get(2));
	//	const Vector3D& p3 = locpoints(trigs(i).Get(3));

	//	double v1x = p2.X()-p1.X();
	//	double v1y = p2.Y()-p1.Y();
	//	double v1z = p2.Z()-p1.Z();

	//	double v2x = p3.X()-p1.X();
	//	double v2y = p3.Y()-p1.Y();
	//	double v2z = p3.Z()-p1.Z();
	//	double px = p1.X();
	//	double py = p1.Y();
	//	double pz = p1.Z();

 //   center.X() += v1x*v1z*(v1x*v2y-v2x*v1y)/12.0+(v2x*v1z+v1x*v2z)*(v1x*v2y-v2x*v1y)/24.0+v2x*v2z*(v1x*v2y-v2x*v1y)/12.0+(px*v1z+v1x*pz)*(v1x*v2y-v2x*v1y)/6.0+(px*v2z+v2x*pz)*(v1x*v2y-v2x*v1y)/6.0+px*pz*(v1x*v2y-v2x*v1y)/2.0;
 //   center.Y() += v1y*v1z*(v1x*v2y-v2x*v1y)/12.0+(v2y*v1z+v1y*v2z)*(v1x*v2y-v2x*v1y)/24.0+v2y*v2z*(v1x*v2y-v2x*v1y)/12.0+(py*v1z+v1y*pz)*(v1x*v2y-v2x*v1y)/6.0+(py*v2z+v2y*pz)*(v1x*v2y-v2x*v1y)/6.0+py*pz*(v1x*v2y-v2x*v1y)/2.0;
 //   center.Z() += v1z*v1z*(v1x*v2y-v2x*v1y)/24.0+v2z*v1z*(v1x*v2y-v2x*v1y)/24.0+v2z*v2z*(v1x*v2y-v2x*v1y)/24.0+pz*v1z*(v1x*v2y-v2x*v1y)/6.0+pz*v2z*(v1x*v2y-v2x*v1y)/6.0+pz*pz*(v1x*v2y-v2x*v1y)/4.0;
	//}

	//if (vol != 0) center *= 1./vol;

	//return center;
}

Matrix3D GeomMesh3D::ComputeMassMomentOfInertia(double rho) const
{
	ElementDataContainer edc, edc_rv;
	ElementData ed;

	ed.SetDouble(rho,"density"); edc.Add(ed);
	Matrix t,pts;
	t.SetSize(trigs.Length(),3);
	pts.SetSize(locpoints.Length(),3);

	for (int i=1; i <= trigs.Length(); i++)
	{
		t(i,1)= trigs(i).Get(1);
		t(i,2)= trigs(i).Get(2);
		t(i,3)= trigs(i).Get(3);
	}

	for (int i=1; i <= locpoints.Length(); i++)
	{
		pts(i,1)= locpoints(i).X();
		pts(i,2)= locpoints(i).Y();
		pts(i,3)= locpoints(i).Z();
	}

	ed.SetMatrix(pts.GetMatPtr(),pts.Getrows(),3,"points"); edc.TreeAdd("MeshData",ed);
	ed.SetMatrix(t.GetMatPtr(),t.Getrows(),3,"triangles"); edc.TreeAdd("MeshData",ed);
	GetMBS()->ComputeInertia(&edc,&edc_rv);

	double * IMp;
	int size = 3;
	edc_rv.TreeGetMatrix("moment_of_inertia",&IMp,size,size);
	Matrix IMmat(size,size,IMp);
	return IMmat;

	//$ DR 2013-01-31, old code, see log 382, computation moved to MultibodySystem
	////see Maple file "maple/integration_triangle_volume.mws" and paper 'computing moments of piecewise polynomial surfaces'
	//Matrix3D IM(0.);

	//double im11=0, im12=0, im13=0, im22=0, im23=0, im33=0;

	//for (int i=1; i <= trigs.Length(); i++)
	//{
	//	const Vector3D& p1 = locpoints(trigs(i).Get(1));
	//	const Vector3D& p2 = locpoints(trigs(i).Get(2));
	//	const Vector3D& p3 = locpoints(trigs(i).Get(3));

	//	double v1x = p2.X()-p1.X();
	//	double v1y = p2.Y()-p1.Y();
	//	double v1z = p2.Z()-p1.Z();

	//	double v2x = p3.X()-p1.X();
	//	double v2y = p3.Y()-p1.Y();
	//	double v2z = p3.Z()-p1.Z();
	//	double px = p1.X();
	//	double py = p1.Y();
	//	double pz = p1.Z();

	//	im11 += v1x*v1x*v1z*(v1x*v2y-v2x*v1y)/20.0+px*px*pz*(v1x*v2y-v2x*v1y)/2.0+(
	//		2.0*px*v1x*v1z+v1x*v1x*pz)*(v1x*v2y-v2x*v1y)/12.0+(2.0*v2x*v1x*v1z+v1x*v1x*v2z)
	//		*(v1x*v2y-v2x*v1y)/60.0+(2.0*px*v2x*v1z+2.0*px*v1x*v2z+2.0*v2x*v1x*pz)*(v1x*v2y
	//		-v2x*v1y)/24.0+(v2x*v2x*v1z+2.0*v2x*v1x*v2z)*(v1x*v2y-v2x*v1y)/60.0+(px*px*v1z+
	//		2.0*px*v1x*pz)*(v1x*v2y-v2x*v1y)/6.0+(px*px*v2z+2.0*px*v2x*pz)*(v1x*v2y-v2x*v1y
	//		)/6.0+(2.0*px*v2x*v2z+v2x*v2x*pz)*(v1x*v2y-v2x*v1y)/12.0+v2x*v2x*v2z*(v1x*v2y-
	//		v2x*v1y)/20.0;
	//	double MapleGenVar1 = ((px*v1y+v1x*py)*v1z+v1x*v1y*pz)*(v1x*v2y-v2x*v1y)/12.0+((
	//		v2x*v1y+v1x*v2y)*v1z+v1x*v1y*v2z)*(v1x*v2y-v2x*v1y)/60.0+((px*v2y+v2x*py)*v1z+(
	//		px*v1y+v1x*py)*v2z+(v2x*v1y+v1x*v2y)*pz)*(v1x*v2y-v2x*v1y)/24.0+(v2x*v2y*v1z+(
	//		v2x*v1y+v1x*v2y)*v2z)*(v1x*v2y-v2x*v1y)/60.0+(px*py*v2z+(px*v2y+v2x*py)*pz)*(
	//		v1x*v2y-v2x*v1y)/6.0;
	//	im12 += MapleGenVar1+((px*v2y+v2x*py)*v2z+v2x*v2y*pz)*(v1x*v2y-v2x*v1y)/12.0
	//		+px*py*pz*(v1x*v2y-v2x*v1y)/2.0+(px*py*v1z+(px*v1y+v1x*py)*pz)*(v1x*v2y-v2x*v1y
	//		)/6.0+v1x*v1y*v1z*(v1x*v2y-v2x*v1y)/20.0+v2x*v2y*v2z*(v1x*v2y-v2x*v1y)/20.0;
	//	im13 += px*pz*pz*(v1x*v2y-v2x*v1y)/4.0+v1x*v1z*v1z*(v1x*v2y-v2x*v1y)/40.0+
	//		v2x*v2z*v2z*(v1x*v2y-v2x*v1y)/40.0+(px*v1z*v1z+2.0*v1x*pz*v1z)*(v1x*v2y-v2x*v1y
	//		)/24.0+(v2x*v1z*v1z+2.0*v1x*v2z*v1z)*(v1x*v2y-v2x*v1y)/120.0+(2.0*(px*v2z+v2x*
	//		pz)*v1z+2.0*v1x*pz*v2z)*(v1x*v2y-v2x*v1y)/48.0+(2.0*v2x*v2z*v1z+v1x*v2z*v2z)*(
	//		v1x*v2y-v2x*v1y)/120.0+(2.0*px*pz*v2z+v2x*pz*pz)*(v1x*v2y-v2x*v1y)/12.0+(px*v2z
	//		*v2z+2.0*v2x*pz*v2z)*(v1x*v2y-v2x*v1y)/24.0+(2.0*px*pz*v1z+v1x*pz*pz)*(v1x*v2y-
	//		v2x*v1y)/12.0;
	//	im22 += py*py*pz*(v1x*v2y-v2x*v1y)/2.0+v2y*v2y*v2z*(v1x*v2y-v2x*v1y)/20.0+
	//		v1y*v1y*v1z*(v1x*v2y-v2x*v1y)/20.0+(2.0*py*v1y*v1z+v1y*v1y*pz)*(v1x*v2y-v2x*v1y
	//		)/12.0+(2.0*v2y*v1y*v1z+v1y*v1y*v2z)*(v1x*v2y-v2x*v1y)/60.0+(2.0*py*v2y*v1z+2.0
	//		*py*v1y*v2z+2.0*v2y*v1y*pz)*(v1x*v2y-v2x*v1y)/24.0+(v2y*v2y*v1z+2.0*v2y*v1y*v2z
	//		)*(v1x*v2y-v2x*v1y)/60.0+(py*py*v1z+2.0*py*v1y*pz)*(v1x*v2y-v2x*v1y)/6.0+(py*py
	//		*v2z+2.0*py*v2y*pz)*(v1x*v2y-v2x*v1y)/6.0+(2.0*py*v2y*v2z+v2y*v2y*pz)*(v1x*v2y-
	//		v2x*v1y)/12.0;
	//	im23 += v2y*v2z*v2z*(v1x*v2y-v2x*v1y)/40.0+py*pz*pz*(v1x*v2y-v2x*v1y)/4.0+
	//		v1y*v1z*v1z*(v1x*v2y-v2x*v1y)/40.0+(py*v1z*v1z+2.0*v1y*pz*v1z)*(v1x*v2y-v2x*v1y
	//		)/24.0+(v2y*v1z*v1z+2.0*v1y*v2z*v1z)*(v1x*v2y-v2x*v1y)/120.0+(2.0*(py*v2z+v2y*
	//		pz)*v1z+2.0*v1y*pz*v2z)*(v1x*v2y-v2x*v1y)/48.0+(2.0*v2y*v2z*v1z+v1y*v2z*v2z)*(
	//		v1x*v2y-v2x*v1y)/120.0+(2.0*py*pz*v2z+v2y*pz*pz)*(v1x*v2y-v2x*v1y)/12.0+(py*v2z
	//		*v2z+2.0*v2y*pz*v2z)*(v1x*v2y-v2x*v1y)/24.0+(2.0*py*pz*v1z+v1y*pz*pz)*(v1x*v2y-
	//		v2x*v1y)/12.0;
	//	im33 += v2z*v1z*v1z*(v1x*v2y-v2x*v1y)/60.0+v1z*v1z*v1z*(v1x*v2y-v2x*v1y)/
	//		60.0+pz*v2z*v1z*(v1x*v2y-v2x*v1y)/12.0+pz*pz*v1z*(v1x*v2y-v2x*v1y)/6.0+v2z*v2z*
	//		v1z*(v1x*v2y-v2x*v1y)/60.0+pz*pz*v2z*(v1x*v2y-v2x*v1y)/6.0+pz*v2z*v2z*(v1x*v2y-
	//		v2x*v1y)/12.0+pz*pz*pz*(v1x*v2y-v2x*v1y)/6.0+v2z*v2z*v2z*(v1x*v2y-v2x*v1y)/60.0
	//		+pz*v1z*v1z*(v1x*v2y-v2x*v1y)/12.0;
	//}

	//IM = rho * Matrix3D(im22+im33,-im12,-im13, -im12, im11+im33, -im23, -im13, -im23, im11+im22);

	//return IM;
}

int GeomMesh3D::PointIsInside(Vector3D p, double& dist) //return 1 if point is inside geometry; return dist to nearest triangle point
{
	int inside = 0;
	dist = 1e100;
	
	//point is inside if it cuts in one direction of a line triangles x many times && x is odd
	Vector3D line(0,0,1);
	Vector3D p1,p2,p3; //triangular points
	Vector3D n; //normal

	Vector3D pp; //projected point

	TArray<Vector3D> cutpoints;
	TArray<Vector3D> normals;

	for (int i=1; i <= NTrigs(); i++)
	{
		GetTrig0(i, p1, p2, p3, n);
		//here we could transform points to global points:
		//p1 = GetPos(p1); //this does not work, because bodys get their XG function working only after Assemble and Initialize ==> could be changed, by assigning to XG in case that system is not initialized, the x_init value ...
		//p2 = GetPos(p2);
		//p3 = GetPos(p3);

		//compute minimum distance to every triangle:
		pp = p;
		double actdist = MinDistTP(p1,p2,p3, pp);
		dist = min(dist, actdist);

		//project line in plane and compute distance
		Normal3D(p1,p2,p3,n); //for safety, if normal is not correct!
		
		pp = p;
		int rv = IntersectLineWithPlane(p1, n, line, pp);
		
		if (rv)
		{
			double md = MinDistTP(p1,p2,p3, pp);
			if (md < 1e-14)
			{
				int found = 0;
				for (int j=1; j<=cutpoints.Length(); j++)
				{
					if ((cutpoints(j) - pp).Norm() < 1e-14) 
					{
						found = j;
					}
				}
				if(!found || (normals(found)*line)*(n*line) < 0)
				{	
					cutpoints.Add(pp);
					normals.Add(n);
				}
			}
			
		}
	}

	// test cutpoints on edges of triangles

	//count number of cuts in one direction:
	int ncuts_in_one_direction = 0;
	for (int i=1; i <= cutpoints.Length(); i++)
	{
		if ((cutpoints(i)-p)*line > 0) ncuts_in_one_direction++;
	}

	if (ncuts_in_one_direction % 2 == 1) {return 1;}
	else {return 0;}
}

int GeomMesh3D::CutInsideOutsidePoints(TArray<Vector3D>& points, double min_wall_dist, int inside)
{
	TArray<Vector3D> pnew;

	for (int i=1; i<=points.Length(); i++)
	{
		double dist = 0;
		int point_inside = PointIsInside(points(i), dist);

		if ((point_inside && inside == 1) || (!point_inside && inside == 0))
		{
			if (dist >= min_wall_dist)
			{
				pnew.Add(points(i));
			}
		}
	}

	int n_cut_points = points.Length()-pnew.Length();

	points = pnew;

	return n_cut_points;
}


void GeomMesh3D::DrawYourself()
{
	if (locpoints.Length() == 0) return ;
	if (!GetMBS()->GetIOption(134) && !GetMBS()->GetIOption(133) && GetElnum() != 0) return;

	RenderContext* pCurrentRC = mbs->GetRC();
	SetDrawPoints();

	int old_transp_mode = GetMBS()->GetTransparency();
	int use_cutting_planes = mbs->UseCuttingPlanes();

	if (elnum)
	{
		drawnormals.SetLen(normals.Length());

		if (GetElement().Dim() == 3)
		{
			if (GetElement().IsRigid())
			{

				Matrix3D A = GetElement().GetRotMatrixD();

				//scale rotation:
				if (GetMBS()->GetIOption(151) && GetElement().IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
				{
					Matrix3D A0 = GetElement().GetRotMatrixInit();
					double fact = GetMBS()->GetDOption(105);

					A = (fact*A-(fact-1.)*A0);
				}

				Vector3D p0 = GetElement().GetRefPosD(); //refpos is already scaled

				for (int i=1; i <= locpoints.Length(); i++)
				{
					ptd(i) = p0 + A * locpoints(i);
				}
				for (int i=1; i <= normals.Length(); i++)
				{
					drawnormals(i) = A * normals(i);
				}
			}
			else
			{
				for (int i=1; i <= locpoints.Length(); i++)
				{
					ptd(i) = GetElement().GetPosD(locpoints(i));
				}
				for (int i=1; i <= normals.Length(); i++)
				{
					drawnormals(i) = ((Body3D&)GetElement()).GetRotMatrixD(locpoints(i)) * normals(i);
				}
			}
		}
		else //2D-Body:
		{
			if (GetElement().IsRigid())
			{

				Matrix3D A = GetElement().GetRotMatrixD();

				////scale rotation:
				//if (GetMBS()->GetIOption(151) && GetBody2D().IsRigid() && GetMBS()->GetDOption(105) != 1.) //for deformation scaling
				//{
				//	Matrix3D A0 = GetBody2D().GetRotMatrixInit();
				//	double fact = GetMBS()->GetDOption(105);

				//	A = (fact*A-(fact-1.)*A0);
				//}

				Vector3D p0 = GetElement().GetRefPosD(); //refpos is already scaled

				for (int i=1; i <= locpoints.Length(); i++)
				{
					ptd(i) = p0 + A * locpoints(i);
				}
				for (int i=1; i <= normals.Length(); i++)
				{
					drawnormals(i) = A * normals(i);
				}
			}
			else
			{
				mbs->UO() << "warning: drawtools, geommesh, draw yourself: GetBody2D().GetRotMatrixD() does not scale\n";
			}
		}
	}	
	else
	{
		CopyPtD();
		drawnormals = normals;
	}

	SetGeomElementColor();

	double thickness = 2+GetMBS()->GetDOption(119); //-->maybe check if rigid body/constraint/background object


	if (!(!GetMBS()->GetIOption(134) && GetElnum() != 0))
	{

		int supersmooth = GetMBS()->GetIOption(132)*smooth_drawing;
		if (!supersmooth) pCurrentRC->glBeginTriangles();
		//if (!supersmooth) pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

		for (int i=1; i <= trigs.Length(); i++)
		{
			const int3& t = trigs(i);
			Vector3D c = 1./3.*(ptd(t(1))+ptd(t(2))+ptd(t(3)));

			if (!use_cutting_planes || mbs->CuttingPlanesAllow(c))
			{

				if (smooth_drawing) //3 normals per triangle
				{
					if (!supersmooth)
					{
						for (int j=1; j <= 3; j++)
						{
							pCurrentRC->glNormal(drawnormals((i-1)*3+j).GetVecPtr());
							pCurrentRC->glVertex(ptd(t(j)).GetVecPtr());
						}
					}
					else
					{ //super-smooth:
						pCurrentRC->glBeginTriangleFan();
						//pCurrentRC->glColor4f((float)col[0],(float)col[1],(float)col[2], transparency);

						//			Vector3D c = 1./3.*(ptd(t(1))+ptd(t(2))+ptd(t(3))); // already calculated
						Vector3D n = Vector3D(drawnormals((i-1)*3+1)+drawnormals((i-1)*3+2)+drawnormals((i-1)*3+3));
						n.Normalize();

						pCurrentRC->glNormal(n.GetVecPtr());
						pCurrentRC->glVertex(c.GetVecPtr());

						for (int j=1; j <= 3; j++)
						{
							pCurrentRC->glNormal(drawnormals((i-1)*3+j).GetVecPtr());
							pCurrentRC->glVertex(ptd(t(j)).GetVecPtr());
							if (j!=3)
							{
								Vector3D n(0.5*drawnormals((i-1)*3+j)+0.5*drawnormals((i-1)*3+j+1));
								n.Normalize();
								pCurrentRC->glNormal(n.GetVecPtr());
								pCurrentRC->glVertex((0.5*ptd(t(j))+0.5*ptd(t(j+1))).GetVecPtr());
							}
							else
							{
								Vector3D n(0.5*drawnormals((i-1)*3+j)+0.5*drawnormals((i-1)*3+1));
								n.Normalize();
								pCurrentRC->glNormal(n.GetVecPtr());
								pCurrentRC->glVertex((0.5*ptd(t(j))+0.5*ptd(t(1))).GetVecPtr());
							}
						}
						pCurrentRC->glNormal(drawnormals((i-1)*3+1).GetVecPtr());
						pCurrentRC->glVertex(ptd(t(1)).GetVecPtr());

						pCurrentRC->glEnd();

					}
				}
				else	//1 normal per triangle
				{
					const Vector3D& n = drawnormals(i);
					pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					pCurrentRC->glVertex((float)ptd(t(1)).X(),(float)ptd(t(1)).Y(),(float)ptd(t(1)).Z());
					pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					pCurrentRC->glVertex((float)ptd(t(2)).X(),(float)ptd(t(2)).Y(),(float)ptd(t(2)).Z());
					pCurrentRC->glNormal((float)n[0],(float)n[1],(float)n[2]);
					pCurrentRC->glVertex((float)ptd(t(3)).X(),(float)ptd(t(3)).Y(),(float)ptd(t(3)).Z());
				}
			}
		}
		if (!supersmooth) pCurrentRC->glEnd();
	}

	int advancedlinemode = 1; //draw nices lines (for presenation), but no redraw during zooming

	GetMBS()->SetTransparency(0);

	if (advancedlinemode) pCurrentRC->BeginLinesOnSolid();

	//draw normals:
	if (0)
	{
		for (int i=1; i <= trigs.Length(); i++)
		{
			double len = 0.05;
			const int3& t = trigs(i);
			for (int j=1; j <= 3; j++)
			{
				Vector3D n, t1, t2;
				n = drawnormals((i-1)*3+j);
				n.SetNormalBasis(t1,t2);
				GetMBS()->MyDrawLine(ptd(t(j)),ptd(t(j))+len*n,thickness);

				GetMBS()->MyDrawLine(ptd(t(j))+len*n, ptd(t(j))+len*0.9*n+(len*0.05)*t1,thickness);
				GetMBS()->MyDrawLine(ptd(t(j))+len*n, ptd(t(j))+len*0.9*n+(-len*0.05)*t1,thickness);
			}
		}
	}

	//draw edges:
	if (!(!GetMBS()->GetIOption(133) /*&& GetElnum() != 0*/))
	{
		for (int i=1; i <= edges.Length(); i++)
		{
			GetMBS()->MyDrawLine(ptd(edges(i).Get(1)),ptd(edges(i).Get(2)),thickness);
		}
	}
	if (advancedlinemode) pCurrentRC->EndLinesOnSolid();

	GetMBS()->SetTransparency(old_transp_mode);

	if (0) //trig numbers
	{
		for (int i=1; i <= trigs.Length(); i++)
		{
			const int3& t = trigs(i);
			Vector3D p = 1./3.*(ptd(t(1))+ptd(t(2))+ptd(t(3)));

			char text[40];
			sprintf(text,"t%d", i);

			GetMBS()->GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
		}
	}

}
