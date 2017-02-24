//#**************************************************************
//#
//# filename:             GeomElements.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						Mai 2005
//# description:          Drawing tools
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
//#**************************************************************

#ifndef DRAWTOOLS__H
#define DRAWTOOLS__H

#include "mbs_interface.h"

//GeomElement
//GeomQuad3D
//GeomZyl3D
//GeomPipe3D
//GeomRotObject3D
//GeomCube3D --> communication changed
//GeomSphere3D
//GeomGround
//GeomCircle2D
//GeomPolygon2D
//GeomTrig3D
//GeomMesh3D

typedef enum {TGeomGeneral = 0, TGeomQuad3D = 1, TGeomZyl3D = 2, TGeomPipe3D = 3, TGeomRotObject3D = 4, TGeomCube3D = 5, TGeomSphere3D = 6, 
              TGeomGround = 7, TGeomCircle2D = 8, TGeomPolygon2D = 9, TGeomText2D = 10, TGeomLine2D = 11, TGeomTrig3D = 12, TGeomMesh3D = 13, TGeomOrthoCube3D = 14} TGeomElement;

typedef enum { TDSDrawOutline = 1,     // GeomPolygon2D - does not scale zoom [pts]
							 TDSFillAreas = 2,       // GeomPolygon2D - disable when using non-convex polygons...
							 TDSHighlightPoints = 4, // GeomPolygon2D scales zoom [m]
							 TDSColoredLines = 8,    // GeomPolygon2D - does not scale zoom [pts]
							 TDSArrowHeadLines = 16, // GeomLine2D - scales zoom [m]
} TDrawStyle;

//either fixed to ground (element==0) or fixed to element with pointer to element
class GeomElement
{
public:
	GeomElement():locpoints(), ptd(), pt(), locpoints2D(), ptd2D(), pt2D(), nodenums() 
	{
		col = Vector3D(0,0,0); elnum = 0; objdata = 0; mbs = 0; draw_dim = Vector3D(0.001,4,0); 
		drawstyle = TDSDrawOutline | TDSFillAreas; pointsize = 0.1; linethickness = 1.;
		transparency = 1; translated = Vector3D(0,0,0); 
		name = GetElementSpec();
	}
	GeomElement(const Vector3D& coli):locpoints(), ptd(), pt(), locpoints2D(), ptd2D(), pt2D(), nodenums() 
	{
		col = coli; elnum = 0; objdata = 0; mbs = 0; draw_dim = Vector3D(0.001,4,0); 
		drawstyle = TDSDrawOutline | TDSFillAreas; pointsize = 0.1; linethickness = 1.;
		transparency = 1; translated = Vector3D(0,0,0); 
		name = GetElementSpec();
	}
	GeomElement(MBS* mbsI):locpoints(), ptd(), pt(), locpoints2D(), ptd2D(), pt2D(), nodenums() 
	{
		mbs = mbsI;
		ElementDefaultConstructorInitialization();
	}
	virtual void ElementDefaultConstructorInitialization();
	
	virtual void CopyFrom(const GeomElement& e)
	{
		col = e.col;
		locpoints = e.locpoints;
		ptd = e.ptd;
		pt = e.pt;
		locpoints2D = e.locpoints2D;
		ptd2D = e.ptd2D;
		pt2D = e.pt2D;
		nodenums = e.nodenums;
		elnum = e.elnum;
		objdata = e.objdata;
		mbs = e.mbs;
		draw_dim = e.draw_dim;
		transparency = e.transparency;
		name = e.name;
		translated = e.translated;
		drawstyle = e.drawstyle;
		pointsize = e.pointsize;
		linethickness = e.linethickness;
	}

	GeomElement(const GeomElement& e): locpoints(), ptd(), pt(), locpoints2D(), ptd2D(), pt2D(), nodenums()
	{
		CopyFrom(e);
	}

	GeomElement& operator=(const GeomElement& e) 
	{
		if (this == &e) {return *this;}
		CopyFrom(e);
		return *this;
	}

	virtual ~GeomElement()
	{
		locpoints.Flush();
		ptd.Flush();
		pt.Flush();
		locpoints2D.Flush();
		ptd2D.Flush();
		pt2D.Flush();
		nodenums.Flush();
	}

	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomElement(*this);
		return ed;
	}

	virtual const mystr& GetName() const {return name;}
	virtual mystr& GetName() {return name;}
	virtual void SetName(const char* nameI) {name = nameI;}
	virtual const char* GetElementSpec() const {return "GeomElement";}
	virtual TGeomElement GetType() const {return TGeomGeneral;}
	
	virtual int& DrawStyle() { return drawstyle; }
	virtual double& PointSize() { return pointsize; }
	virtual double& LineThickness() { return linethickness; }

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void Translate(const Vector3D& translate) {assert(0);};
	virtual double ComputeVolume() const {assert(0); return 0;};
	virtual Vector3D ComputeCenterOfMass() const {assert(0); return Vector3D(0.,0.,0.);};
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const {assert(0); return Matrix3D(0.);};

	virtual Box3D GetBoundingBoxD() const 
	{
		Box3D b;
		AddBoundingBoxD(b);
		return b;
	}

	virtual Box3D GetBoundingBox() const 
	{
		Box3D b;
		AddBoundingBox(b);
		return b;
	}

	virtual Box2D GetBoundingBoxD2D() const 
	{
		Box2D b;
		AddBoundingBoxD2D(b);
		return b;
	}

	virtual Box2D GetBoundingBox2D() const 
	{
		Box2D b;
		AddBoundingBox2D(b);
		return b;
	}


	//3D: +++++++++++++++++++++++++++++++++++++
	virtual void CopyPt()
	{
		for (int i=1; i <= locpoints.Length(); i++)
		{
			pt(i) = locpoints(i);
		}
	}

	virtual void CopyPtD()
	{
		for (int i=1; i <= locpoints.Length(); i++)
		{
			ptd(i) = locpoints(i);
		}
	}

	virtual void SetDrawPoints();
	virtual void SetComputePoints();

	//2D: +++++++++++++++++++++++++++++++++++++
	virtual void CopyPt2D();
	virtual void CopyPtD2D();
	virtual void SetDrawPoints2D();
	virtual void SetComputePoints2D();

	virtual void DrawYourself() {}; //DrawYourself(WCDInterface::RenderContext* pCurrentRC)


	//to transform locpoints into deformed or rotated configuration:
	virtual const TArray<Vector3D> & GetPoints() {return locpoints;}

	virtual TArray<Vector3D> & GetPointsTransformedD() {return ptd;}
	virtual const TArray<Vector3D> & GetPointsTransformedD() const {return ptd;}

	virtual TArray<Vector3D> & GetPointsTransformed() {return pt;}
	virtual const TArray<Vector3D> & GetPointsTransformed() const {return pt;}

	virtual const TArray<Vector2D> & GetPoints2D() {return locpoints2D;}

	virtual TArray<Vector2D> & GetPointsTransformedD2D() {return ptd2D;}
	virtual const TArray<Vector2D> & GetPointsTransformedD2D() const {return ptd2D;}

	virtual TArray<Vector2D> & GetPointsTransformed2D() {return pt2D;}
	virtual const TArray<Vector2D> & GetPointsTransformed2D() const {return pt2D;}


	virtual const Vector3D& PT(int i) const {return pt(i);}
	virtual const Vector3D& PTD(int i) const {return ptd(i);}

	virtual const Vector2D& PT2D(int i) const {return pt2D(i);}
	virtual const Vector2D& PTD2D(int i) const {return ptd2D(i);}

	//e.g. for cutting plane and others
	virtual Vector3D GetRefPosD(/*int ind, const Vector3D& pglob*/) const
	{
		return GetBoundingBoxD().Center();
	}

	//local position
	virtual Vector3D GetLocPoint(int i) const 
	{
		return locpoints(i);
	}

	virtual Vector2D GetLocPoint2D(int i) const 
	{
		return locpoints2D(i);
	}

	virtual int AddLocPoint(const Vector3D& p) {return locpoints.Add(p);}
	virtual int AddLocPoint2D(const Vector2D& p) {return locpoints2D.Add(p);}


	//actual (global) position:
	virtual Vector3D GetTPoint(int i) const;
	virtual Vector3D GetTPointD(int i) const;
	virtual Vector2D GetTPoint2D(int i) const;
	virtual Vector2D GetTPointD2D(int i) const;
	virtual Vector3D ToP3D(const Vector2D& p) const;
	virtual Vector2D GetPos2D(const Vector2D& ploc) const;
	virtual Vector2D GetVel2D(const Vector2D& ploc) const;

	virtual Vector3D GetPos(const Vector3D& ploc) const;
	virtual Vector3D GetVel(const Vector3D& ploc) const;

	virtual int NP() const {return locpoints.Length() + locpoints2D.Length();} //one of the two arrays must be zero!!!

	virtual const Element& GetElement() const;
	//$ YV 2012-12-10: commented out
	/*virtual const Body2D& GetBody2D() const;
	virtual Body2D& GetBody2D();
	virtual const Body3D& GetBody3D() const;
	virtual Body3D& GetBody3D();*/
	virtual MBS* GetMBS() const;

	virtual int GetElnum() const {return elnum;}
	virtual void SetElnum(int en) {elnum = en;}
	virtual int GetObjdata() const {return objdata;}
	virtual void SetObjectdata(int od) {objdata = od;}

	virtual int GetNodeNum(int i) const 
	{
		if (nodenums.Length() == locpoints.Length()) return nodenums(i);
		else return 0;
	}

	virtual void SetDrawParam(const Vector3D& v) {draw_dim = v;}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//functions to be set by each GeomElement:

	//add bounding boxes to existing box
	virtual void AddBoundingBoxD(Box3D& box) const 
	{
		if (locpoints2D.Length() != 0)
		{
			for (int i=1; i <= locpoints2D.Length(); i++) box.Add(ToP3D(GetTPointD2D(i)));
		}
		else
			for (int i=1; i <= locpoints.Length(); i++) box.Add(GetTPointD(i));
	}

	virtual void AddBoundingBox(Box3D& box) const 
	{
		if (locpoints2D.Length() != 0)
		{
			for (int i=1; i <= locpoints2D.Length(); i++) box.Add(ToP3D(GetTPoint2D(i)));
		}
		else
			for (int i=1; i <= locpoints.Length(); i++) box.Add(GetTPoint(i));
	}

	virtual void AddBoundingBoxD2D(Box2D& box) const 
	{
		for (int i=1; i <= locpoints2D.Length(); i++) box.Add(GetTPointD2D(i));
	}

	virtual void AddBoundingBox2D(Box2D& box) const 
	{
		for (int i=1; i <= locpoints2D.Length(); i++) box.Add(GetTPoint2D(i));
	}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp) {assert(0); return 0;}

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector2D GetLocPos(int ind, const Vector2D& pglob) const {assert(0); return Vector2D(0.,0.);};

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector2D GetNormal(int ind, const Vector2D& ploc) const {assert(0); return Vector2D(0.,0.);};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//3D:
	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector3D& p, int& ind, Vector3D& pp) {assert(0); return 0;}

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector3D GetLocPos(int ind, const Vector3D& pglob) const {assert(0); return Vector3D(0.,0.,0.);};

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector3D GetNormal(int ind, const Vector3D& ploc) const {assert(0); return Vector3D(0.,0.,0.);};

	//get (deformed) normalized unique tangent vector at locpos
	virtual Vector3D GetTangent(int ind, const Vector3D& ploc) const {assert(0); return Vector3D(0.,0.,0.);};

	virtual const Vector3D& GetCol() const {return col;}
	virtual void SetCol(const Vector3D& colI) {col = colI;}
	virtual void SetGeomElementColor();
	virtual float GetTransparency() const {return transparency;}
	virtual void SetTransparency(float transparencyI) {transparency = transparencyI;}

protected:
	mystr name;								 //name for load/save, dialogs
	Vector3D col;
	float transparency;				 //object transparency; 0=transparent, 1=solid
	Vector3D draw_dim;				 //(linewidth, #segments for curved shapes, 0)
	MBS* mbs;
	int elnum;								 //element number
	int objdata;							 //position refers to objectdata in MBS
	TArray<int> nodenums;			 //2D and 3D node number; 0 if does not exist
	Vector3D translated;			 //Geometry has been translated (e.g. for center of mass)
	
	int drawstyle;						 // drawstyle according to enum
	double pointsize;          // drawsize of highlighted points
	double linethickness;      // thickness of colored lines

	//3D:
	TArray<Vector3D> locpoints;//local points
	TArray<Vector3D> ptd;		   //transformed points
	TArray<Vector3D> pt;       //transformed points for computation
	//2D:
	TArray<Vector2D> locpoints2D; //local points
	TArray<Vector2D> ptd2D;    //transformed points
	TArray<Vector2D> pt2D;     //transformed points for computation
};




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class GeomText2D:public GeomElement
{
public:
	GeomText2D():GeomElement()
	{
		locpoints2D.SetLen(0);
		ptd2D.SetLen(0);
		pt2D.SetLen(0);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomText2D(*this);
		return ed;
	}
	GeomText2D(const GeomText2D& e):GeomElement(e) {CopyFrom(e);}
	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomText2D& ce = (const GeomText2D&)e;

		row = ce.row;
		alignment = ce.alignment;
		text = new char[strlen(ce.text)+1];
		strcpy(text, ce.text);
	}

	GeomText2D(MBS* mbsi, int elnumi, const char* textI, int rowI, int alignmentI, const Vector3D& coli);

	virtual ~GeomText2D()
	{
		delete[] text;
		GeomElement::~GeomElement();
	}
	virtual const char* GetElementSpec() const {return "GeomText2D";}
	virtual TGeomElement GetType() const {return TGeomText2D;}

	virtual void DrawYourself();

private:
	int alignment;
	int row;
	char* text;
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class GeomLine2D:public GeomElement
{
public:
	GeomLine2D():GeomElement()
	{
		locpoints2D.SetLen(2);
		ptd2D.SetLen(2);
		pt2D.SetLen(2);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomLine2D(*this);
		return ed;
	}
	GeomLine2D(const GeomLine2D& e):GeomElement(e) {CopyFrom(e);}

	GeomLine2D(MBS* mbsi, int elnumi, const Vector2D& p1, const Vector2D& p2, const Vector3D& coli);

	GeomLine2D(MBS* mbsi, int elnumi, int n1, int n2, const Vector3D& coli);

	virtual const char* GetElementSpec() const {return "GeomLine2D";}
	virtual TGeomElement GetType() const {return TGeomLine2D;}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp);

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector2D GetLocPos(int ind, const Vector2D& pglob) const;

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector2D GetNormal(int ind, const Vector2D& ploc) const;
	virtual void DrawYourself();

	//private:
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class GeomPolygon2D:public GeomElement
{
public:
	GeomPolygon2D():GeomElement()
	{}

	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomPolygon2D(*this);
		return ed;
	}
	GeomPolygon2D(const GeomPolygon2D& e):GeomElement(e) {CopyFrom(e);}

	GeomPolygon2D(MBS* mbsi, int elnumi, const TArray<Vector2D>& p, const Vector3D& coli);
	//GeomPolygon2D(MBS* mbsi, int elnumi, int n1, int n2, const Vector3D& coli);

	virtual const char* GetElementSpec() const {return "GeomPolygon2D";}
	virtual TGeomElement GetType() const {return TGeomPolygon2D;}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp) {assert(0); return 0;};

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector2D GetLocPos(int ind, const Vector2D& pglob) const {assert(0); return Vector2D(0,0);};

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector2D GetNormal(int ind, const Vector2D& ploc) const {assert(0); return Vector2D(0,0);}
	virtual void DrawYourself();
	virtual int ReadFromFile(mystr& filename);
	//private:
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class GeomCircle2D:public GeomElement
{
public:
	GeomCircle2D():GeomElement()
	{
		locpoints2D.SetLen(2);
		ptd2D.SetLen(2);
		pt2D.SetLen(2);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomCircle2D(*this);
		return ed;
	}
	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomCircle2D& ce = (const GeomCircle2D&)e;

		r = ce.r;
	}

	GeomCircle2D(const GeomCircle2D& e):GeomElement(e) {CopyFrom(e);}

	//full circle:
	GeomCircle2D(MBS* mbsi, int elnumi, const Vector2D& p1, double rad, const Vector3D& coli);

	//only segment:
	GeomCircle2D(MBS* mbsi, int elnumi, const Vector2D& p1, const Vector2D& pseg1, const Vector2D& pseg2, double rad, const Vector3D& coli);

	virtual const char* GetElementSpec() const {return "GeomCircle2D";}
	virtual TGeomElement GetType() const {return TGeomCircle2D;}

	//add bounding boxes to existing box
	virtual void AddBoundingBoxD(Box3D& box) const 
	{
		Vector2D p(GetTPointD2D(1));
		box.Add(Box3D(ToP3D(p+Vector2D(-r,-r)),ToP3D(p+Vector2D(r,r))));
	}

	virtual void AddBoundingBox(Box3D& box) const 
	{
		Vector2D p(GetTPoint2D(1));
		box.Add(Box3D(ToP3D(p+Vector2D(-r,-r)),ToP3D(p+Vector2D(r,r))));
	}

	virtual void AddBoundingBoxD2D(Box2D& box) const 
	{
		Vector2D p(GetTPointD2D(1));
		box.Add(Box2D(p+Vector2D(-r,-r),p+Vector2D(r,r)));
	}

	virtual void AddBoundingBox2D(Box2D& box) const 
	{
		Vector2D p(GetTPoint2D(1));
		box.Add(Box2D(p+Vector2D(-r,-r),p+Vector2D(r,r)));
	}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector2D& p, int& ind, Vector2D& pp);

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector2D GetLocPos(int ind, const Vector2D& pglob) const;

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector2D GetNormal(int ind, const Vector2D& ploc) const;
	virtual void DrawYourself();

private:
	double r; //radius of circle
};







//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//             3D GeomQuad
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomQuad3D:public GeomElement
{
public:
	GeomQuad3D():GeomElement()
	{
		locpoints.SetLen(4);
		ptd.SetLen(4);
		pt.SetLen(4);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomQuad3D(*this);
		return ed;
	}
	GeomQuad3D(const GeomQuad3D& e):GeomElement(e) {CopyFrom(e);}

	GeomQuad3D(MBS* mbsi, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, const Vector3D& coli):GeomElement()
	{
		mbs = mbsi;
		col = coli;
		elnum = 0;

		locpoints.SetLen(4);
		ptd.SetLen(4);
		pt.SetLen(4);

		locpoints(1) = p1;
		locpoints(2) = p2;
		locpoints(3) = p3;
		locpoints(4) = p4;
		CopyPt();
		CopyPtD();
	}
	virtual const char* GetElementSpec() const {return "GeomQuad3D";}
	virtual TGeomElement GetType() const {return TGeomQuad3D;}

	virtual void DrawYourself();


	//private:
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//             3D GeomTrig
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomTrig3D:public GeomElement
{
public:
	GeomTrig3D():GeomElement()
	{
		locpoints.SetLen(3);
		ptd.SetLen(3);
		pt.SetLen(3);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomTrig3D(*this);
		return ed;
	}
	GeomTrig3D(const GeomTrig3D& e):GeomElement(e) {CopyFrom(e);}

	GeomTrig3D(MBS* mbsI, int elnumI, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& colI):GeomElement()
	{
		mbs = mbsI;
		col = colI;
		elnum = elnumI;

		locpoints.SetLen(3);
		ptd.SetLen(3);
		pt.SetLen(3);

		locpoints(1) = p1;
		locpoints(2) = p2;
		locpoints(3) = p3;
		CopyPt();
		CopyPtD();
	}

	virtual const char* GetElementSpec() const {return "GeomTrig3D";}
	virtual TGeomElement GetType() const {return TGeomTrig3D;}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector3D& p, int& ind, Vector3D& pp);

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector3D GetLocPos(int ind, const Vector3D& pglob) const;

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector3D GetNormal(int ind, const Vector3D& ploc) const;

	//get (deformed) normalized unique tangent vector at locpos
	virtual Vector3D GetTangent(int ind, const Vector3D& ploc) const;

	virtual void DrawYourself();


	//private:
};





//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomCube3D:public GeomElement
{
public:
	//binary ordering of locpoints: X=0/1: i+=1, Y=0/1: i+=2, Z=0/1: i+=4
	//#outer surface right: 2-4-8-6
	//#outer surface left:  1-3-7-5
	//#bottom: 1-2-4-3
	//#top: 5-6-8-7
	GeomCube3D():GeomElement()
	{
		locpoints.SetLen(8);
		ptd.SetLen(8);
		pt.SetLen(8);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomCube3D(*this);
		return ed;
	}
	GeomCube3D(const GeomCube3D& e):GeomElement(e) {CopyFrom(e);}

	GeomCube3D(MBS* mbsi, const Vector3D& pmid, const Vector3D& size, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		Vector3D size2 = 0.5*size;
		SetLocPoints(pmid, size2);
		CopyPt();
		CopyPtD();
	}

	//rectangular ground-shape, 4 individual heights
	GeomCube3D(MBS* mbsi, const Vector3D& p0, double xsize, double ysize, double z1, double z2, double z3, double z4, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		locpoints(1) = p0+Vector3D(0.   ,0.   ,0.);
		locpoints(2) = p0+Vector3D(xsize,0.   ,0.);
		locpoints(3) = p0+Vector3D(0.   ,ysize,0.);
		locpoints(4) = p0+Vector3D(xsize,ysize,0.);
		locpoints(5) = p0+Vector3D(0.   ,0.   ,z1);
		locpoints(6) = p0+Vector3D(xsize,0.   ,z2);
		locpoints(7) = p0+Vector3D(0.   ,ysize,z3);
		locpoints(8) = p0+Vector3D(xsize,ysize,z4);
		CopyPt();
		CopyPtD();
	}

	//arbitrary ground-shape, one height
	GeomCube3D(MBS* mbsi, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, double z, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		locpoints(1) = p1;
		locpoints(2) = p2;
		locpoints(3) = p3;
		locpoints(4) = p4;
		locpoints(5) = p1+Vector3D(0.,0.,z);
		locpoints(6) = p2+Vector3D(0.,0.,z);
		locpoints(7) = p3+Vector3D(0.,0.,z);
		locpoints(8) = p4+Vector3D(0.,0.,z);
		CopyPt();
		CopyPtD();
	}

	GeomCube3D(MBS* mbsi, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, 
		const Vector3D& p5, const Vector3D& p6, const Vector3D& p7, const Vector3D& p8, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		locpoints(1) = p1;
		locpoints(2) = p2;
		locpoints(3) = p3;
		locpoints(4) = p4;
		locpoints(5) = p5;
		locpoints(6) = p6;
		locpoints(7) = p7;
		locpoints(8) = p8;
		CopyPt();
		CopyPtD();
	}

	virtual const char* GetElementSpec() const {return "GeomCube3D";}
	virtual TGeomElement GetType() const {return TGeomCube3D;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void GetTrigs(TArray<int3>& trigs) const;
	virtual void SetLocPoints(const Vector3D center_point, const Vector3D size)
	{
		locpoints(1) = center_point+Vector3D(-size.X(),-size.Y(),-size.Z());
		locpoints(2) = center_point+Vector3D( size.X(),-size.Y(),-size.Z());
		locpoints(3) = center_point+Vector3D(-size.X(), size.Y(),-size.Z());
		locpoints(4) = center_point+Vector3D( size.X(), size.Y(),-size.Z());
		locpoints(5) = center_point+Vector3D(-size.X(),-size.Y(), size.Z());
		locpoints(6) = center_point+Vector3D( size.X(),-size.Y(), size.Z());
		locpoints(7) = center_point+Vector3D(-size.X(), size.Y(), size.Z());
		locpoints(8) = center_point+Vector3D( size.X(), size.Y(), size.Z());
	}
	virtual void Translate(const Vector3D& translate);
	virtual void Rotate(const Matrix3D& rotate);
	virtual double ComputeVolume() const;
	virtual Vector3D ComputeCenterOfMass() const;
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const;

	virtual void DrawYourself();

	//private:
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// orthogonal cube with stored center point and rotation matrix 
class GeomOrthoCube3D:public GeomCube3D
{
public:
	//binary ordering of locpoints: X=0/1: i+=1, Y=0/1: i+=2, Z=0/1: i+=4
	//#outer surface right: 2-4-8-6
	//#outer surface left:  1-3-7-5
	//#bottom: 1-2-4-3
	//#top: 5-6-8-7
	//#center point stored
	GeomOrthoCube3D():GeomCube3D()
	{
		locpoints.SetLen(9);
		ptd.SetLen(9);
		pt.SetLen(9);
		A = Matrix3D(1.);
		size = Vector3D(1.);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomOrthoCube3D(*this);
		
		return ed;
	}
	GeomOrthoCube3D(const GeomOrthoCube3D& e):GeomCube3D(e) {CopyFrom(e);A = e.A;size = e.size;}

	GeomOrthoCube3D(MBS* mbsi, const Vector3D& pmid, const Vector3D& sizei, const Matrix3D& Ai, const Vector3D& coli):GeomCube3D(mbsi, pmid, sizei, coli)
	{
		locpoints.SetLen(9);
		ptd.SetLen(9);
		pt.SetLen(9);		
		A = Ai;		
		size = sizei; 
		SetLocPoints(pmid, size, A);
	}

	virtual const char* GetElementSpec() const {return "GeomOrthoCube3D";}
	virtual TGeomElement GetType() const {return TGeomOrthoCube3D;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void SetLocPoints(const Vector3D center_point, const Vector3D sizei, const Matrix3D A);
	virtual void Translate(const Vector3D& translate){GeomCube3D::Translate(translate);};
	virtual void Rotate(const Matrix3D& rotate);
	virtual double ComputeVolume() const{return size(1)*size(2)*size(3);};
	virtual Vector3D ComputeCenterOfMass() const{return locpoints(9);};
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const{return GeomCube3D::ComputeMassMomentOfInertia(rho);};

	virtual void DrawYourself(){GeomCube3D::DrawYourself();}

	private:
		Matrix3D A;
		Vector3D size;
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ SW 2013-11-13: extended GeomZyl3D such that it can now draw hollow cylinders
class GeomZyl3D:public GeomElement
{ 
public:
	GeomZyl3D():GeomElement()
	{
		locpoints.SetLen(2);
		ptd.SetLen(2);
		pt.SetLen(2);
		ri = 0.; 
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomZyl3D(*this);
		return ed;
	}
	GeomZyl3D(const GeomZyl3D& e):GeomElement(e) {CopyFrom(e);}

	GeomZyl3D(MBS* mbsi, const Vector3D& pz1, const Vector3D& pz2, double rzi, int tilei, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		tile = tilei;
		rz = rzi;
		locpoints(1) = pz1;
		locpoints(2) = pz2;
		CopyPt();
		CopyPtD();
		ri = 0.;
	}

	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomZyl3D& ce = (const GeomZyl3D&)e;

		tile =ce.tile;
		rz = ce.rz;
		ri = ce.ri;
	}

	virtual const char* GetElementSpec() const {return "GeomCylinder3D";}
	virtual TGeomElement GetType() const {return TGeomZyl3D;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void Translate(const Vector3D& translate);
	virtual double ComputeVolume() const;
	virtual Vector3D ComputeCenterOfMass() const;
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const;

	virtual void DrawYourself();

private:
	int tile;
	double rz;
	double ri;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//A full or partial pipe
//draw pipe from phimin to phimax with respect to Vector vrad
class GeomPipe3D:public GeomElement
{
public:
	GeomPipe3D():GeomElement()
	{
		locpoints.SetLen(3);
		ptd.SetLen(3);
		pt.SetLen(3);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomPipe3D(*this);
		return ed;
	}
	GeomPipe3D(const GeomPipe3D& e):GeomElement(e) {CopyFrom(e);}

	GeomPipe3D(MBS* mbsi, const Vector3D& pz0, const Vector3D& vaxis, const Vector3D& vrad, double riI, double roI, int tileI, 
		const Vector3D& coli, double phiminI=0, double phimaxI=2.*MY_PI, int drawendfaceI=1):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = coli;
		tile = tileI;
		if (tile == 0) tile = 1;

		ri = riI;
		ro = roI;
		phimin = phiminI;
		phimax = phimaxI;
		drawendface = drawendfaceI;

		Vector3D va = vaxis;
		Vector3D vr = vrad;
		vr.Normalize();

		locpoints(1) = pz0;
		locpoints(2) = pz0 + va;
		locpoints(3) = pz0 + vr;
		CopyPt();
		CopyPtD();
	}

	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomPipe3D& ce = (const GeomPipe3D&)e;

		tile =ce.tile;
		ri = ce.ri;
		ro = ce.ro;
		phimin = ce.phimin;
		phimax = ce.phimax;
		drawendface = ce.drawendface;
	}

	virtual const char* GetElementSpec() const {return "GeomPipe3D";}
	virtual TGeomElement GetType() const {return TGeomPipe3D;}

	virtual void DrawYourself();

private:
	int tile;
	double ri, ro; //inner and outer radius
	double phimin, phimax;
	int drawendface;
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomSphere3D:public GeomElement
{
public:
	GeomSphere3D():GeomElement()
	{
		locpoints.SetLen(1);
		ptd.SetLen(1);
		pt.SetLen(1);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomSphere3D(*this);
		return ed;
	}
	GeomSphere3D(const GeomSphere3D& e):GeomElement(e) {CopyFrom(e);}

	GeomSphere3D(MBS* mbsi, int elnumI, const Vector3D& pm, double rzi, int tilei, const Vector3D& coli):GeomElement(mbsi)
	{
		elnum = elnumI;
		mbs = mbsi;
		col = coli;
		draw_dim.Y() = tilei;
		rz = rzi;

		locpoints.SetLen(1);
		ptd.SetLen(1);
		pt.SetLen(1);
		locpoints(1) = pm;

		CopyPt();
		CopyPtD();
	}

	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomSphere3D& ce = (const GeomSphere3D&)e;

		rz = ce.rz;
	}

	virtual const char* GetElementSpec() const {return "GeomSphere3D";}
	virtual TGeomElement GetType() const {return TGeomSphere3D;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void Translate(const Vector3D& translate);
	virtual double ComputeVolume() const;
	virtual Vector3D ComputeCenterOfMass() const;
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const;

	//add bounding boxes to existing box
	virtual void AddBoundingBoxD(Box3D& box) const 
	{
		Vector3D p(GetTPointD(1));
		box.Add(Box3D(p+Vector3D(-rz,-rz,-rz),p+Vector3D(rz,rz,rz)));
	}

	virtual void AddBoundingBox(Box3D& box) const 
	{
		Vector3D p(GetTPoint(1));
		box.Add(Box3D(p+Vector3D(-rz,-rz,-rz),p+Vector3D(rz,rz,rz)));
	}

	//Get nearest global point pp on element from global point p, ind is the index of found segment, returns the signed gap/penetration
	virtual double GetNearestPoint(const Vector3D& p, int& ind, Vector3D& pp);

	//get local position on element from global position and index (for polygon etc.)
	virtual Vector3D GetLocPos(int ind, const Vector3D& pglob) const;

	//get (deformed) normalized normal vector (outwards) at locpos
	virtual Vector3D GetNormal(int ind, const Vector3D& ploc) const;

	//get (deformed) normalized unique tangent vector at locpos
	virtual Vector3D GetTangent(int ind, const Vector3D& ploc) const;

	virtual void DrawYourself();

	virtual const double& GetRadius() const {return rz;}

private:
	double rz;
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//rotation object
//a 2D curve builds a rotation object
//the standard rotation axis is z-axis
class GeomRotObject3D:public GeomElement
{
public:
	GeomRotObject3D():GeomElement()
	{
		locpoints.SetLen(0);
		ptd.SetLen(0);
		pt.SetLen(0);
		segments.SetLen(0);
		colorseg.SetLen(0);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomRotObject3D(*this);
		return ed;
	}
	GeomRotObject3D(const GeomRotObject3D& e):GeomElement(e) {CopyFrom(e);}

	GeomRotObject3D(MBS* mbsi, const Vector3D& p0, TArray<Vector2D> segs, TArray<Vector3D> cols, const Matrix3D& rotmat, int tilei);

	virtual const char* GetElementSpec() const {return "GeomRotObject3D";}
	virtual TGeomElement GetType() const {return TGeomRotObject3D;}

	virtual void CopyFrom(const GeomElement& e);

	virtual void DrawYourself();

private:
	int tile;
	double rz;
	TArray<Vector2D> segments;
	TArray<Vector3D> colorseg;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//triangular surface mesh for drawing and contact:
class GeomMesh3D:public GeomElement
{
public:
	GeomMesh3D():GeomElement(), normals(), drawnormals(), trigs(), points_to_trigs(), points_to_edges(), edges() 
	{
		Init();
	}
	GeomMesh3D(MBS* mbsi, int elnumi, const Vector3D& coli):GeomElement(mbsi), normals(), drawnormals(), trigs(), points_to_trigs(), edges() 
	{
		Init();
		mbs = mbsi;
		elnum = elnumi;
		col = coli;
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomMesh3D(*this);
		return ed;
	}
	GeomMesh3D(const GeomMesh3D& e):GeomElement(e) {CopyFrom(e);}
	virtual void CopyFrom(const GeomElement& e);
	virtual ~GeomMesh3D();

	virtual const char* GetElementSpec() const {return "GeomMesh3D";}
	virtual TGeomElement GetType() const {return TGeomMesh3D;}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(const ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void Init();

	virtual int NTrigs() const;
	virtual void GetTrig0(int i, Vector3D& p1, Vector3D& p2, Vector3D& p3, Vector3D& normal);
	virtual int AddTrig(int n1, int n2, int n3, const Vector3D& normal);

	//virtual int AddPoint(const Vector3D& p1, double eps = 1e-16); //returns point number if exists, or adds new point
	virtual void ClearMesh();
	virtual Vector3D GeometryNormal(int trig1) const;
	virtual double GetAngle(int trig1, int trig2) const;
	virtual double GetArea(int trig) const;
	virtual int IsEdge(int2 edge) const; //check if two point indices represent an edge, return edge number or zero

	virtual void ReadSTLMesh(char* file, int binary=0);
//	virtual void ReadFEMesh3D(class FEMesh3D* femesh);
	virtual void ReadFromPointArray(TArray<Vector3D>* v3darr);
	virtual void ReadFromPointSequence(Vector3D* v3dseq, int n);
	virtual void ReadFromDoubleSequence(double* dseq, int n);


	//virtual void SetSTLtolerance(double tol) {stl_tolerance = tol;}
	virtual void SetDrawEdgeAngle(double edge_angle) {drawedge_angle = edge_angle;}
	virtual void SetSmoothDrawing(int smooth_draw = 1) {smooth_drawing = smooth_draw;}

	virtual int GetNeighbourTrig(int trig, int side, int& ntrig, int& nside); //get number of triangle neighbouring at 'side'
	virtual void FinishMesh(); //call after all triangles have been added

	virtual void SetGeomModification(const Vector3D& resize, const Matrix3D& rotate, const Vector3D& translate);

	virtual void Translate(const Vector3D& translate);
	virtual void Rotation(const Matrix3D &rotate);
	virtual void Stretch(const double factor);
	virtual double ComputeVolume() const;
	virtual Vector3D ComputeCenterOfMass() const;
	virtual Matrix3D ComputeMassMomentOfInertia(double rho = 1) const;

	virtual int PointIsInside(Vector3D p, double& dist); //return 1 if point is inside geometry; return dist to nearest triangle point;
	virtual int PointIsInside(Vector3D p) {double dist = 0; return PointIsInside(p, dist);}; //return 1 if point is inside geometry; return dist to nearest triangle point

	//cut points which are outside (inside) of mesh geometry; if inside=1, then cut inside points; returns number of cut points
	//points are allowed to have distance min_wall_dist to any triangle
	virtual int CutInsideOutsidePoints(TArray<Vector3D>& points, double min_wall_dist, int inside = 1);

	virtual void DrawYourself();

private:
	TArray<Vector3D> normals;						//Normals of triangles
	TArray<Vector3D> drawnormals;				//Normals of triangles, for drawing only: drawnormals = A*normals
	TArray<int3> trigs;									//Triangles to points
	TArray<IVector*> points_to_trigs;   //triangles that belong to a point
	TArray<IVector*> points_to_edges;   //edges that belong to a point
	TArray<int2> edges;									//draw edges

	int smooth_drawing;									//if 1, each triangle has 3 own normals, otherwise only one
	double drawedge_angle;							// if > MY_PI, do not draw edges!
	//double stl_tolerance;								//tolerance for equal points when reading
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomGround:public GeomElement
{
public:
	GeomGround():GeomElement()
	{
		col = Vector3D(0,0,0);
		locpoints.SetLen(0);
		ptd.SetLen(0);
		pt.SetLen(0);
		quads.SetLen(0);
		boxes.SetLen(0);
	}
	GeomGround(MBS* mbsi):GeomElement(mbsi)
	{
		elnum = 0;
		mbs = mbsi;
		col = Vector3D(0,0,0);
		locpoints.SetLen(0);
		ptd.SetLen(0);
		pt.SetLen(0);
		quads.SetLen(0);
		boxes.SetLen(0);
	}
	virtual ~GeomGround()
	{
		quads.Flush();
		cols.Flush();
		boxes.Flush();
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomGround(*this);
		return ed;
	}

	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomGround& ce = (const GeomGround&)e;

		quads.Flush();
		for (int i=1; i <= quads.Length(); i++)
			quads.Add(ce.quads(i));

		cols.Flush();
		for (int i=1; i <= cols.Length(); i++)
			cols.Add(ce.cols(i));

		boxes.Flush();
		for (int i=1; i <= boxes.Length(); i++)
			boxes.Add(ce.boxes(i));

		tree = ce.tree;
	}
	GeomGround(const GeomGround& e):GeomElement(e) {CopyFrom(e);}

	virtual Vector3D FStreet(double x, int mode, double sx, double sy, double sz) const;

	virtual void SetStreet();
	virtual void SetFlat(const Vector3D& dim, int resx, int resy);

	virtual void DrawYourself();

	virtual int NQuads() const {return quads.Length();}
	virtual int4 GetQuad(int i) const {return quads(i);}

	virtual const Box3D& Boxes(int i) const {return boxes(i);}

	virtual void GetQuadsInBox(const Box3D& b, TArray<int>& items)
	{
		tree.GetItemsInBox(b,items);
	}

private:
	TArray<int4> quads;
	TArray<Vector3D> cols;
	TArray<Box3D> boxes;
	SearchTree tree;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeomShip3D:public GeomElement
{
public:
	GeomShip3D():GeomElement()
	{
		locpoints.SetLen(2);
		ptd.SetLen(2);
		pt.SetLen(2);
	}
	virtual GeomElement* GetCopy() const
	{
		GeomElement* ed = new GeomShip3D(*this);
		return ed;
	}
	GeomShip3D(const GeomShip3D& e):GeomElement(e) {CopyFrom(e);}

	GeomShip3D(MBS* mbsi, int elnumi, const Vector3D& sizei, const Vector3D& comi, double tilei, const Vector3D& coli): GeomElement(mbsi)
	{
		elnum = elnumi;
		mbs = mbsi;
		col = coli;

		size = sizei;
		com = comi;
		tile = (int)tilei;

		locpoints.SetLen(2);
		ptd.SetLen(2);
		pt.SetLen(2);

		locpoints(1) = size-com;
		locpoints(2) = -1*(size-com);
		CopyPt();
		CopyPtD();

	}

	virtual void CopyFrom(const GeomElement& e)
	{
		GeomElement::CopyFrom(e);
		const GeomShip3D& ce = (const GeomShip3D&)e;

		tile =ce.tile;
		size = ce.size;
		com = ce.com;
	}

	virtual void DrawYourself();

private:
	Vector3D size; //dimensions of ship
	int tile;
	Vector3D com; //center of mass

};

#endif