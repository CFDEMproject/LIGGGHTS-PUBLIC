//#**************************************************************
//#
//# filename:             csg_geometry.h
//#
//# author:               Gerstmayr Johannes, Peter Gruber, Daniel Reischl
//#
//# generated:						Feb 2012
//# description:          Collection of data structures and methods, which provide help with 
//#                       reading & writing CSG geometry files.
//# 
//# remarks:						  CSG means 'Constructive Solid Geometry', e.g. '*.geo' files in NETGEN
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

#pragma once

#include "mbs_interface.h"



class GEOObject
{
public:
	GEOObject() {};

	virtual void Translate(const Vector3D& v)
	{
		for (int i=1; i <= points.Length(); i++)
		{
			points(i) += v;
		}
	}

	virtual void Rotate(double ang, int axis)
	{
		Matrix3D rot;
		if (axis == 1) rot = RotMatrix1(ang);
		else if (axis == 2) rot = RotMatrix2(ang);
		else if (axis == 3) rot = RotMatrix3(ang);

		for (int i=1; i <= points.Length(); i++)
		{
			points(i) = rot*points(i);
		}
		for (int i=1; i <= normals.Length(); i++)
		{
			normals(i) = rot*normals(i);
		}
	}

	virtual void Export(ostream& os, const mystr& name) {}

	virtual void WritePoint(ostream& os, int i)
	{
		os << points(i).X() << ", " << points(i).Y() << ", " << points(i).Z();
	}
	virtual void WriteNormal(ostream& os, int i)
	{
		os << normals(i).X() << ", " << normals(i).Y() << ", " << normals(i).Z();
	}

protected:
	TArray<Vector3D> points;
	TArray<Vector3D> normals;
};

class GEOCube: public GEOObject
{
public:
	GEOCube() {points.SetLen(6); normals.SetLen(6);};

	virtual void SetOrthoCube(const Vector3D& pmid, const Vector3D& size)
	{
		Vector3D p1 = pmid - 0.5*size;
		Vector3D p2 = pmid + 0.5*size;

		points(1) = p1; points(2) = p1; points(3) = p1;
		points(4) = p2; points(5) = p2; points(6) = p2;

		normals(1) = Vector3D(0,0,-1);
		normals(2) = Vector3D(0,-1,0);
		normals(3) = Vector3D(-1,0,0);
		normals(4) = Vector3D(0,0,1);
		normals(5) = Vector3D(0,1,0);
		normals(6) = Vector3D(1,0,0);
	}

	virtual void SetCube(const TArray<Vector3D>& corners)
	{
// numbering taken from FEHex.GetSide4
		points(1) = corners(1); Normal3D( corners(1), corners(5), corners(7), normals(1));
		points(2) = corners(2); Normal3D( corners(2), corners(4), corners(8), normals(2));
		points(3) = corners(1); Normal3D( corners(1), corners(2), corners(6), normals(3));
		points(4) = corners(3); Normal3D( corners(3), corners(7), corners(8), normals(4));
		points(5) = corners(1); Normal3D( corners(1), corners(3), corners(4), normals(5));
		points(6) = corners(5); Normal3D( corners(5), corners(6), corners(8), normals(6));
	}

	virtual void Export(ostream& os, const mystr& name)
	{
		os << "solid " << name << " = ";// << GetNetgenString();
		for (int i = 1; i <= 6; i++)
		{
			if (i != 1) os << "         and ";
			os << "plane(" << points(i).X() << ", " << points(i).Y() << ", " << points(i).Z() << "; ";
			os << normals(i).X() << ", " << normals(i).Y() << ", " << normals(i).Z() << ")";
			if (i == 6) os << ";\n\n";
			else os << "\n";
		}
	}

	//virtual mystr GetNetgenString()
	//{
	//	mystr buffer("");
	//  for (int i = 1; i <= 6; i++)
	//	{
	//		if (i != 1) buffer += "         and ";
	//		buffer += "plane(" + mystr(points(i).X()) + ", " + mystr(points(i).Y()) + ", " + mystr(points(i).Z()) + "; ";
	//		buffer += mystr(normals(i).X()) + ", " + mystr(normals(i).Y()) + ", " + mystr(normals(i).Z()) + ")";
	//		if (i == 6) buffer += ";\n\n";
	//		else buffer += "\n";
	//	}
	//	return buffer;
	//}
};

class GEOTorus: public GEOObject
{
public:
	GEOTorus() {points.SetLen(1); normals.SetLen(1);};

	virtual void SetTorus(const Vector3D& p1, const Vector3D& n1, double r_out, double r_in)
	{
		points(1) = p1; 
		normals(1) = n1;
		rad_outer = r_out;
		rad_inner = r_in;
	}

	virtual void Export(ostream& os, const mystr& name)
	{
		os << "solid " << name << " = ";

		Vector3D n1 = points(1)-points(2); n1.Normalize();
		Vector3D n2 = points(2)-points(1); n2.Normalize();

		os << "torus("; WritePoint(os,1); 
		os << "; "; WriteNormal(os,1);
		os << ";" << rad_outer << ";" << rad_inner << ");\n\n";
	}
private:
	double rad_outer, rad_inner;
};

class GEOCyl: public GEOObject
{
public:
	GEOCyl() {points.SetLen(2); normals.SetLen(0);};

	virtual void SetCyl(const Vector3D& p1, const Vector3D& p2, double r)
	{
		points(1) = p1; 
		points(2) = p2; 
		rad = r;
	}

	virtual void Export(ostream& os, const mystr& name)
	{
		os << "solid " << name << " = ";

		Vector3D n1 = points(1)-points(2); n1.Normalize();
		Vector3D n2 = points(2)-points(1); n2.Normalize();

		os << "cylinder("; WritePoint(os,1); 
		os << "; "; WritePoint(os,2); os << ";" << rad << ") and\n";
		os << "         plane("; WritePoint(os,1);
		os << "; " << n1.X() << ", " << n1.Y() << ", " << n1.Z() << ") and\n";
		os << "         plane("; WritePoint(os,2);
		os << "; " << n2.X() << ", " << n2.Y() << ", " << n2.Z() << ");\n\n";

	}
private:
	double rad;
};

class GEOSphere: public GEOObject
{
public:
	GEOSphere() {points.SetLen(1);};

	virtual void SetSphere(const Vector3D& p1, double r)
	{
		points(1) = p1; 
		rad = r;
	}

	virtual void Export(ostream& os, const mystr& name)
	{
		os << "solid " << name << " = ";
		os << "sphere("; WritePoint(os,1); 
		os << ";" << rad << ");\n\n";
	}
private:
	double rad;
};


//$ DR 2012-02-24
// writes and executes bat-file to mesh geometry, defined in geo-file, in NETGEN
bool CreateMeshWithNETGEN(const char * geoFilename,const char * meshFilename,const char * NETGENPath,const char * WorkDir,const char * NETGENMeshFileType, const char* coarsity="coarse");

//$ DR 2012-04
// generates geo-file, for use in NETGEN, to define geometry of a rotor
bool ExportAsGeoFileCylinders(const char * filename, TArray<Vector3D> PointsLeft, TArray<Vector3D> PointsRight, TArray<int> Material, double maxh=1, const char * rotor_name="SimpleRotor", const char * title="geo-file for rotor, generated in HOTINT", const char * comment="some additional comments");

//$ DR 2012-04
bool ExportAsMeshedCylPartEDCtxtFile(const char * filename, TArray<Vector3D> PointsLeft, TArray<Vector3D> PointsRight, TArray<int> Material, const char * rotor_name="SimpleRotor", const char * title="autogenerated part of a model file for HOTINT", const char * comment="some additional comments");

//$ AD 2012-22
// generates an entry for a geo-file to define a hexahedral block
mystr GetGeoFileEntryForCube( TArray<Vector3D>& points, int number );