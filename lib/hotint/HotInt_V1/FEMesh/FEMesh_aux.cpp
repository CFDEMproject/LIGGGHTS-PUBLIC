//#**************************************************************
//#
//# filename:             FEMesh_aux.cpp
//#
//# author:               Gerstmayr Johannes
//#                       Dorninger Alexander
//#
//# generated:						November 2011
//# description:          Finite Element mesh
//# remarks:						  * original file split
//#                       * FEMesh_aux contains       variety of helper classes, 
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
 
#include "mbs_interface.h"
#include "stepsettings.h"
#include "FEMesh.h"
#include "myfile.h"


//#include "windows.h" //for shell execute
#include <direct.h>  //for getcwd
#include  <io.h>     //file operation _access
#include  <stdio.h>
#include  <stdlib.h>
#include "FEMesh_aux.h"

#include "mbsload.h"
#include "femathhelperfunctions.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ (base) class FEElement:  derive linear: {FEHex, FETet, FEQuad, FETrig, FELine} and the respective quadratic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodenumbers, functions to generate faces, domain, material, color
//
// classes to process finite elements within FEMesh
// TODO: 
// * remove functions:
// * implement functions: GetLocalPos
// THINK ABOUT:
// * default return local position (0,0,0) COULD be confusing
// * connect to mbs thus no meed to pass list of nodes for geometry

// returns local coordinates of a (global) node, -1 if node is not in the element
Vector3D FEElement::GetGlobalNodeLocalCoord(int nodenumber_global)
{
	int local = -1;
	for(int i=1; i<= NNodes(); i++)
	{
		if(this->GetNode(i) == nodenumber_global)
		{
			local = i;
			break;
		}
	}
	return GetNodeLocalCoord(local);
}

// returns local coordinates of a (local) node,  if node is not in the element return (0.,0.,0.)
Vector3D FEElement::GetNodeLocalCoord(int nodenumber_local)
{
	Vector3D loccoord(0.,0.,0.);

	switch(Type())
	{
	case TFELine: // Line
		switch(nodenumber_local)
		{
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 0.5, 0.); break;
		default:return Vector3D(0.);
		}
		break;
	case TFETrig: // Trigonal
		switch(nodenumber_local)
		{
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 1., 0.); break;
		case 4: return Vector3D(0.5, 0. , 0. ); break;
		case 5: return Vector3D(0.5, 0.5, 0. ); break;
		case 6: return Vector3D(0. , 0.5, 0. ); break;
		default: return Vector3D(0.);
		}
		break;
	case TFEQuad: // Quadrilateral
		//loccoord = Quadrilateral::GetNodeLocPos(nodenumber); // can not use this function (non static member function...)
	  switch(nodenumber_local)
		{ //1..4 for linear element, 1..9 for quadratic element
		case 1: return Vector3D( 1., 1., 0.); break;
		case 2: return Vector3D(-1., 1., 0.); break;
		case 3: return Vector3D(-1.,-1., 0.); break;
		case 4: return Vector3D( 1.,-1., 0.); break;
		case 5: return Vector3D( 0., 1., 0.); break;
		case 6: return Vector3D(-1., 0., 0.); break;
		case 7: return Vector3D( 0.,-1., 0.); break;
		case 8: return Vector3D( 1., 0., 0.); break;
		case 9: return Vector3D( 0., 0., 0.); break;
		default: return Vector3D(0.);
		}
		break;
	case TFETet: // Tetrahedral
		//loccoord = Tetrahedral::GetNodeLocPos(nodenumber); // can not use this function (non static member function...)
		switch(nodenumber_local)
		{
		//linear:
		case 1: return Vector3D(0., 0., 0.); break;
		case 2: return Vector3D(1., 0., 0.); break;
		case 3: return Vector3D(0., 1., 0.); break;
		case 4: return Vector3D(0., 0., 1.); break;
		//quadratic:
		case 5: return Vector3D(0.5, 0. , 0. ); break;
		case 6: return Vector3D(0.5, 0.5, 0. ); break;
		case 7: return Vector3D(0. , 0.5, 0. ); break;
		case 8: return Vector3D(0.5, 0. , 0.5); break;
		case 9: return Vector3D(0. , 0.5, 0.5); break;
		case 10:return Vector3D(0. , 0. , 0.5); break;
		default: return Vector3D(0.);
		}
		break;
	case TFEHex: // Hexahedral
		//loccoord = Hexahedral::GetNodeLocPos(nodenumber); // can not use this function (non static member function...)
			//generate local node position hierarchically:
		int i = nodenumber_local;
		if (i <= 8)
		{
			Vector3D lp(-1.,-1.,-1.);
			int j = i-1;

			if (j&1) lp.X() += 2;
			if (j&2) lp.Y() += 2;
			if (j&4) lp.Z() += 2;
			return lp;
		}
		else
		{
			if (i <= 12) //X=0-midnode
			{
				Vector3D lp(0.,-1.,-1.);
				int j = i-9;
				if (j&1) lp.Y() += 2;
				if (j&2) lp.Z() += 2;
				return lp;
			}
			else if (i <= 16) //Y=0-midnode
			{
				Vector3D lp(-1.,0.,-1.);
				int j = i-13;
				if (j&1) lp.X() += 2;
				if (j&2) lp.Z() += 2;
				return lp;
			}
			else //Z=0-midnode
			{
				Vector3D lp(-1.,-1.,0.);
				int j = i-17;
				if (j&1) lp.X() += 2;
				if (j&2) lp.Y() += 2;
				return lp;
			}
		}
		break;
	}
	return loccoord;
}

//returns the coordinates of the center point of the finite element
Vector3D FEElement::GetCenterPoint(const TArray<FEMesh_Node>& nodes)
{
	Vector3D center(0.,0.,0.);
	for(int i=1; i<=NNodes(); i++)
	{
		center += nodes(GetNode(i)).GetCoords3D();
	}
	center *= 1./NNodes();
	return center;
}

// returns sidenumber of finite element if global point is on that side - returns only first side
int FEElement::PointIsOnSide(const Vector3D& globalpos, const TArray<FEMesh_Node>& nodes)
{
	double dist = 0.;
	int4 sidenodes;
	if(Type() == TFEHex)
	{
		for(int i=1; i<=NSides(); i++)
		{
			sidenodes = GetSideNodeNum4(i);
			dist = ::DistToQuad(nodes(sidenodes(1)).GetCoords3D(),nodes(sidenodes(2)).GetCoords3D(),
				                  nodes(sidenodes(3)).GetCoords3D(),nodes(sidenodes(4)).GetCoords3D(), globalpos);
			if( abs(dist) < MESH_STD_TOL)
			{
				return i;
			}
		}
	}
	return 0;
}

// !!HEX ONLY!!: returns local coordinates of a given global position - assuming that this point is on a side on the finite element
Vector3D FEElement::GetLocalPosOnSide(const Vector3D& globalpos, const TArray<FEMesh_Node>& nodes)
{
	int is_on_side = PointIsOnSide(globalpos,nodes);
	Vector3D localpos;

	if(Type() == TFEHex)
	{
		if(is_on_side != 0)
		{
			int4 sidenodes = GetSideNodeNum4(is_on_side);
// interpolation 
			double y_loc, x_loc;
			Vector3D p1 = nodes(sidenodes(1)).GetCoords3D();
			Vector3D p2 = nodes(sidenodes(2)).GetCoords3D();
			Vector3D p3 = nodes(sidenodes(3)).GetCoords3D();
			Vector3D p4 = nodes(sidenodes(4)).GetCoords3D();

			Vector3D edge12 = p2-p1; //nodepos(sidenodes(2)) - nodepos(sidenodes(1));
			Vector3D edge43 = p3-p4; //nodepos(sidenodes(3)) - nodepos(sidenodes(4));

			Vector3D edge14 = p4-p1; //nodepos(sidenodes(4)) - nodepos(sidenodes(1));
			Vector3D edge23 = p3-p2; //nodepos(sidenodes(3)) - nodepos(sidenodes(2));
			
			double dot_12_43 = abs( edge12 * edge43 / (edge12.Norm() * edge43.Norm()) );
			double dot_14_23 = abs( edge14 * edge23 / (edge14.Norm() * edge23.Norm()) );

			if( abs(dot_12_43-1.) < MESH_STD_TOL ) // edge12 is par. to edge 43
			{
				// parallelogram or trapezoid
				Vector3D pt_at_23;
				Vector3D pt_at_14;
				Vector2D lm;
				int one;
				one = ::IntersectLines3D( globalpos, edge12,     p4, edge14,    pt_at_14,    lm);
				assert(one==1);
				one = ::IntersectLines3D( globalpos, edge12,     p3, edge23,    pt_at_23,    lm);
				assert(one==1);

				double ly1 = (p1-pt_at_14).Norm();
				double ly2 = (p4-pt_at_14).Norm();
				y_loc = (ly1 - ly2) / (ly1 + ly2); // from -1 @p1(ly1 = 0) to +1 @p4(ly2 = 0)
			
				double lx1 = (globalpos-pt_at_14).Norm();
				double lx2 = (globalpos-pt_at_23).Norm();
				x_loc = (lx1 - lx2) / (lx1 + lx2); // from -1 @p14 to +1 @p23
			}
			else if( abs(dot_14_23-1.) < MESH_STD_TOL)
			{
				// parallelogram or trapezoid
				Vector3D pt_at_12;
				Vector3D pt_at_43;
				Vector2D lm;
				int one;
				one = ::IntersectLines3D( globalpos, edge23,     p1, edge12,    pt_at_12,    lm);
				assert(one==1);
				one = ::IntersectLines3D( globalpos, edge23,     p4, edge43,    pt_at_43,    lm);
				assert(one==1);

				double lx1 = (p1-pt_at_12).Norm();
				double lx2 = (p2-pt_at_12).Norm();
				x_loc = (lx1 - lx2) / (lx1 + lx2); // from -1 @p1 to +1 @p2

				double ly1 = (globalpos-pt_at_12).Norm();
				double ly2 = (globalpos-pt_at_43).Norm();
				y_loc = (ly1 - ly2) / (ly1 + ly2); // from -1 @p12 to +1 @p43
			}
			else
			{
				// oblique
				Vector3D inter;
				Vector3D pt_at_12;
				Vector3D pt_at_43;
				Vector2D lm;
				int one;
				one = ::IntersectLines3D( p1, edge14,    p2, edge23,   inter,    lm);
				assert(one==1);

				Vector3D dir = inter-globalpos;
				one = ::IntersectLines3D( globalpos, dir,     p1, edge12,    pt_at_12,    lm);
				assert(one==1);
				one = ::IntersectLines3D( globalpos, dir,     p4, edge43,    pt_at_43,    lm);
				assert(one==1);

				double lx1 = (p1-pt_at_12).Norm();
				double lx2 = (p2-pt_at_12).Norm();
				x_loc = (lx1 - lx2) / (lx1 + lx2); // from -1 @p1 to +1 @p2

				double ly1 = (globalpos-pt_at_12).Norm();
				double ly2 = (globalpos-pt_at_43).Norm();
				y_loc = (ly1 - ly2) / (ly1 + ly2); // from -1 @p12 to +1 @p43
			}

			switch(is_on_side)
			{
			case 1: localpos = Vector3D( -1., y_loc, x_loc); break; // (1,5,7,3) !! rotate x: 1->5 || z
			case 2: localpos = Vector3D(  1., x_loc, y_loc); break; // (2,4,8,6) 
			
			case 3: localpos = Vector3D(x_loc, -1., y_loc); break; // (1,2,6,5)
			case 4: localpos = Vector3D(y_loc,  1., x_loc); break; // (3,7,8,4) !! rotate x: 3->7 || z
			
			case 5: localpos = Vector3D(y_loc, x_loc, -1.); break; // (1,3,4,2) !! rotate x: 1->3 || y
			case 6: localpos = Vector3D(x_loc, y_loc,  1.); break; // (5,6,8,7)
			}
		}
		else
		{
			localpos = GetLocalPos(globalpos, nodes);
		}
	}
  return localpos;
}


// !! NOT IMPLEMENTED YET!! returns local coordinates of a given global position
Vector3D FEElement::GetLocalPos(const Vector3D& globalpos, const TArray<FEMesh_Node>& nodes)
{
	Vector3D localpos;
	// not implemented yet
	// 
	return localpos;
}

//called by DrawSystem
void FEElement::DrawElement(MBS* mbs, const TArray<Vector3D>& node_ref_pos, const TArray<Vector3D>& node_disp)
{
  for (int i=1; i <= NSides(); i++)
	{
		int ntrig = 1; //number of triangles
		int3 side_trig[2]; //first drawing triangle

		if (Type() == 4) //Hex
		{
			ntrig = 2;
			int4 side = GetSideNodeNum4(i);
			side_trig[0] = int3(side(1),side(2),side(3));
			side_trig[1] = int3(side(1),side(3),side(4));
		}
		else if (Type() == 3) //Tet
		{
			ntrig = 1;
			side_trig[0] = GetSideNodeNum3(i);
		}
		else if (Type() == 5) //Prism
		{
			ntrig = 1;
			if (NSides() <= NSides3()) 
			{
				side_trig[0] = GetSideNodeNum3(i);
			}
			else
			{
				ntrig = 2;

				int4 side = GetSideNodeNum4(i-NSides3());
				side_trig[0] = int3(side(1),side(2),side(3));
				side_trig[1] = int3(side(1),side(3),side(4));
			}
		}

		for (int j = 1; j <= ntrig; j++)
		{
			int3 trig = side_trig[j-1];

			Vector3D v[3];
			for (int k = 1; k <= 3; k++)
			{
				int nn = trig.Get(k);
				v[k-1] = node_ref_pos(nn) + node_disp(nn);
			}

			mbs->DrawTrig(v[0],v[1],v[2]);
		}
	}
}
//+ end FEElement and derived
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class ClosedPolygon
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: list of points (numbers), list of line sections(numbers)
// lists are stored in counterclockwise sequence, a closed polygon has only ONE domain 

int ClosedPolygon::Rotate(int n)
{
	PointNrs().Erase(PointNrs().Length());
  if(n>0)
	{
		for(int i=1; i<=n; i++)
		{
			int lastpoint = PointNrs().Last();
			PointNrs().Insert(1,lastpoint);
			PointNrs().Erase(PointNrs().Length());
			int lastline = LineNrs().Last();
			LineNrs().Insert(1,lastline);
			LineNrs().Erase(LineNrs().Length());
		}	
	}
	else
	{
	  for(int i=1; i<=(-n); i++)
		{
			PointNrs().Add(PointNr(1));
			PointNrs().Erase(PointNrs().Length());
			LineNrs().Add(LineNr(1));
			LineNrs().Erase(LineNrs().Length());
		}
	}
	PointNrs().Add(PointNr(1));
	return PointNr(1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  class Crossection
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: list of points, list of line sections, list of closed polygons

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Import from file
// read_boundaries = 1 looks for matching entries in the edc Boundaries, which has to be filled first
int Crosssection::ImportFromIN2D(mystr& filename, int read_boundaries)
{
	MBS* mbs = GetMBS();
	CMFile file_in2d(filename,TFMread);
	
	if (!filename.Length())
	{
		mbs->UO(UO_LVL_err) << "No Filename for .in2d geometry specified !\n";
		return -2;
	}
	if (!file_in2d.IsGood())
	{
		mbs->UO(UO_LVL_err) << "Could not open File: " << filename;
		mbs->UO().InstantMessageText("Could not open .in2d geometry file");
		return -2;
	}

	mystr entire, line, word;
	file_in2d.RWF(entire);
	int pos=0;

	int foundheader = false;
	int grading = 0;
	int nr_points = 0;
	int nr_splines = 0;

	int in2d_version = 0;

// IN2D HEADER: 
// line 1 : splinecurves2d or splinedcurves2dv2
// line 2 : parameter grading
// find and read header
	while(pos!=-1)
	{
		entire.GetUntilEOL(pos, '\n', line, 1);
		if(pos==-1) { mbs->UO(UO_LVL_err) << mystr("EOF: no header found in ") + filename + mystr("\n"); return -2; }
		if(line[0]=='#') continue;       // skip comment line
		if(line.Length()==0) continue;   // skip empty line
		if(in2d_version == 2 && ( line[0] < '0' || line[0] > '9' )) continue; // skip entries "points" etc. from version 2

		if (line.Compare("splinecurves2dv2\n") || line.Compare("splinecurves2dv2 \n") ) 
		{
			foundheader = true;
			in2d_version = 2;
			mbs->UO(UO_LVL_ext) << mystr("In2D V2 - Header identified in file ") + filename +  mystr("\n");
			break;
		}
		else if (line.Compare("splinecurves2d\n") || line.Compare("splinecurves2d \n")) 
		{
			foundheader = true;
			in2d_version = 1;
			mbs->UO(UO_LVL_ext) << mystr("In2D - Header identified in file ") + filename +  mystr("\n");
			break;
		}
		if (!foundheader) continue;      // skip garbage
	}
	while(pos!=-1)
	{
		entire.GetUntilEOL(pos, '\n', line, 1);
		if(pos==-1) { mbs->UO(UO_LVL_err) << mystr("EOF: no entry grading found in ") + filename + mystr("\n"); return -2; }
		if (line.IsValidNumber())
		{
			grading = line.MakeInt();
			break;
		}
	}

// IN2D DATA - Points
// VERSION 1:
// line 1 : number of points
// line 2+: point data { x y ref }
// VERSION 2:
// line 1: "POINTS"
// line 2+: { pointnumber x y }
	int pointsread = 0;  
	int skip_once = 0;
	while(pos!=-1)
	{
// exit loop version 1: reached number of points as specified in first line
		if(in2d_version == 1 && pointsread != 0 && pointsread >= nr_points ) 
			break;             

		entire.GetUntilEOL(pos, '\n', line, 1);
		if(pos==-1) { mbs->UO(UO_LVL_err) << mystr("EOF: unexpected end of file in point data") + filename + mystr("\n"); return -2; }
		if(line[0]=='#') continue;       // skip comment line
		if(line.Compare("\n")) continue;   // skip empty line
		if(in2d_version == 2 && ( line[0] < '0' || line[0] > '9' && skip_once==0)) {skip_once++; continue;} // skip entries "points" etc. from version 2

//exit loop version 2: reached a text entry
		if(in2d_version == 2 && ( line[0] < '0' || line[0] > '9' && skip_once!=0)) 
			break;   

		if(line.IsValidNumber() && nr_points == 0) // first line of block 
		{
			nr_points = line.MakeInt();
			pointsread = 0;
		}
		else              // parse line - Assert
		{
		  Vector2D p;
			int pos_line = 0;
//  point number
			if(in2d_version == 2)
			{
				word = line.GetWord(pos_line);
				if(!word.IsValidNumber()) {	mbs->UO(UO_LVL_err) << mystr("FORMAT:line cannot be parsed: ") + line + mystr("\n"); return -2; }
				int nr = word.MakeInt();
			}
//  x coordinate
			word = line.GetWord(pos_line);
			if(!word.IsValidNumber()) {	mbs->UO(UO_LVL_err) << mystr("FORMAT:line cannot be parsed: ") + line + mystr("\n"); return -2; }
			p.X() = word.MakeDouble();
//  y coordinate
			word = line.GetWord(pos_line);
			if(!word.IsValidNumber()) {	mbs->UO(UO_LVL_err) << mystr("FORMAT:line cannot be parsed: ") + line + mystr("\n"); return -2; }
			p.Y() = word.MakeDouble();
//  refinement
			if(in2d_version == 1)
			{
				// skip for now
			}

// Add to Array
			Points().Add(p);
			pointsread++;
		}
	}

// IN2D DATA - (spline) segments
// VERSION 1:
// line 1 : number of segments
// line 2+: segment data { domainleft domainright nrpoints(2|3) point1 point2 (point3)  }
// VERSION 2:
// line 1: "SEGMENTS"
// line 2+: segment data { domainleft domainright nrpoints(2|3) point1 point2 (point3) optional: refinement,bc }
	int splinesread = 0;
	skip_once = 0;
	while(pos!=-1)
	{
// exit loop version 1: reached number of points as specified in first line
		if(in2d_version == 1 && splinesread!= 0 && splinesread >= nr_splines) 
			break;             

		entire.GetUntilEOL(pos, '\n', line, 1);
		if(pos==-1) { mbs->UO(UO_LVL_err) << mystr("EOF: unexpected end of file in segment data") + filename + mystr("\n"); return -2; }
		if(line[0]=='#') continue;       // skip comment line
		if(line.Compare("\n")) continue;   // skip empty line
		if(in2d_version == 2 && ( line[0] < '0' || line[0] > '9' && skip_once==0)) {skip_once++; continue;} // skip entries "points" etc. from version 2

//exit loop version 2: reached a text entry
		if(in2d_version == 2 && ( line[0] < '0' || line[0] > '9' && skip_once!=0)) 
			break;    

		if(line.IsValidNumber() && nr_splines == 0) // first line of block 
		{
			nr_splines = line.MakeInt();
			splinesread = 0;
		}
		else             // parse line 
		{
			int pos_line = 0;
			int domain1, domain2, p1,p2,psp, bc=0,ref=0;
			word = line.GetWord(pos_line);      domain1 = word.MakeInt();
			word = line.GetWord(pos_line);		  domain2 = word.MakeInt();
			word = line.GetWord(pos_line);      int pts = word.MakeInt();
			if(pts<2 || pts >3) { mbs->UO(UO_LVL_err) << mystr("FORMAT: spline segment may have '2' or '3' as 3rd parameter") + mystr("\n") + line + mystr("\n"); return -2; }
			word = line.GetWord(pos_line);      p1 = word.MakeInt();
			word = line.GetWord(pos_line);      p2 = word.MakeInt();
			if(pts == 3)
			{
				word = line.GetWord(pos_line);    psp = word.MakeInt();
			}
// DONE WITH NON-OPTIONAL CONTENT !
      while(pos_line!=-1)
			{
				word = line.GetWord(pos_line);
				if(pos_line!=-1)
				{
					if(word.IsValidNumber())
					{
						ref = word.MakeInt();					
					}
					else if(word[0]=='-')
					{
						pos_line++;
						word = line.GetWord(pos_line);
						// old version:
						//if(word[0] == 'F') bc = 1;
						//if(word[0] == 'K') bc = 2;
						//if(word[0] == 'S') bc = 3;

						if(read_boundaries)	//$ DR 2012-11-20 added this flag
						{
							// new versionCrosssection::ImportFromIN2d
							if(Boundaries().Length())
							{
								bc = Boundaries().Find(word);	// name of bc is mapped to integer of bc
								if(!bc)
								{
									mbs->UO(UO_LVL_err) << "Could not find the boundary condition '" << word << "' written in "<< filename << "in the edc of boundary conditions!\n";
								}
							}
							else
							{
								mbs->UO(UO_LVL_err) << "No boundary conditions defined before Crosssection::ImportFromIN2D was called.";
							}
						}
					}
				}
			}
			// skip entry refinement
						
// Add to Array
			if(pts == 2) { Lines().Add(LineSegment(mbs, p1, p2, domain1, domain2, bc, ref)); }
			if(pts == 3) { Lines().Add(LineSegment(mbs, p1, psp, p2, domain1, domain2, bc, ref)); }
			splinesread++;
		}
	}

  return 1;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Read the boundary conditions from file
int Crosssection::ReadBoundaryConditions(mystr& filename) //$ DR 2012-08-16
{
	return GetMBS()->File2EDC(filename, &Boundaries());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find all closed polygons
int Crosssection::FindDomainPolygons() // find closed polygons, each containing only a single domain
{
	Polygons().Flush();
	// find highest domain number in all linesegments
	int highest = HighestDomainNumber();

	for(int i=1; i<=highest; i++)
	{
		TArray<int> left, right; // Line segments available for closed polygon
		FilterLineSegmentsForDomainNumber(i, left, right);

		while(left.Length() + right.Length() != 0) // have still lines available for given domain number
		{
			TArray<int> pointsequence;
			TArray<int> linesequence;
// pick starting point
			if(left.Length() != 0)
			{
				linesequence.Add( left(1) );
				pointsequence.Add( Line(left(1)).P1() );
				pointsequence.Add( Line(left(1)).P2() );
				left.Erase(1);
			}
			else
			{
				linesequence.Add( right(1) );
				pointsequence.Add( Line(right(1)).P2() );
				pointsequence.Add( Line(right(1)).P1() );
				right.Erase(1);
			}

			while(pointsequence(1) != pointsequence.Last()) // while not closed 
			{
				int next = 0;
				next = FindContinuation(pointsequence.Last(), left, 1);
				if(next)
				{
					linesequence.Add(left(next));
					pointsequence.Add(Line(left(next)).P2());
					left.Erase(next);					
				}
				else
				{
					next = FindContinuation(pointsequence.Last(), right, 2);
					if(next)
					{
						linesequence.Add(right(next));
						pointsequence.Add(Line(right(next)).P1());
						right.Erase(next);					
					}
				}
				if(!next) continue; // no continuation found !!!
			} // end while
// have closed polygon here
			ClosedPolygon the_polygon(mbs, pointsequence, linesequence);
			the_polygon.Domain() = i;
			Polygons().Add(the_polygon);
		}
	}
	
	for(int i=1; i<=Polygons().Length(); i++)
	{
		ResortPolygonCCW(i);
	}

	// get inner + outer edges
	// jeweils IVector



	return Polygons().Length();
}

void Crosssection::GetPolygonOuterAndInnerEdge()
{
	// ATTENTION: it is assumed that the ClosedPolygons are already sorted counterclockwise --> ResortPolygonCCW
	// the meaning of counterclockwise: assume that the x-axis shows to the right and y-axis up.

	// result: poly->InnerEdge and poly->OuterEdge contain the numbers of the points
	// the meaning of "inner": small y-coordinate

	ClosedPolygon* poly;
	for(int i=1; i<=Polygons().Length(); i++)
	{
		// find starting point (=point with minimum x-coordinate)
		poly = &Polygon(i);
		int nrMinX=1;
		double min = Point(poly->PointNr(1)).X();
		poly->InnerEdge().Flush();
		poly->OuterEdge().Flush();

//$ YV 2013-01-12: no SolverUO here; another solution is needed if necessary
//		if(GetMBS()->SolverUO().GetGlobalMessageLevel() >= UO_LVL_dbg1)
		{
			for(int j=1; j<=poly->NPoints()-1; j++)
			{
				GetMBS()->UO(UO_LVL_ext) << "Polygons("<<i<<"): point "<<j<<", node="<< poly->PointNr(j) <<" ("<<Point(poly->PointNr(j)).X()<<","<<Point(poly->PointNr(j)).Y()<<")\n";
			}
		}

		for(int j=1; j<=poly->NPoints()-1; j++)
		{
			if(Point(poly->PointNr(j)).X() <= min)	// ###
			{
				min = Point(poly->PointNr(j)).X();
				nrMinX = j;
			}
		}		// "min" is the minimum x-coordinate of the poly, 
				// nrMinX is the corresp. nr of the point

		if(nrMinX==(poly->NPoints()-1))	// nrMinX is found as the last point, so it can also be (one of) the first node(s)
		{
			for (int j=1; j<=(poly->NPoints()-1);j++)
			{
				if(Point(poly->PointNr(j)).X() > min)
				{
					nrMinX = j-1;
					j = poly->NPoints(); // to break the for loop
				}
			}
			if(nrMinX==0) {nrMinX=poly->NPoints()-1;}	// if it is really the last point
		}

		// find the last point of the inner edge. it is the last one before increasing x-coordinate again
		int lastInner = nrMinX;
		poly->InnerEdge().Add(poly->PointNr(nrMinX));
		for(int j=nrMinX+1;j<=poly->NPoints()-1;j++)		//from starting point to end of polygon
		{
			lastInner = j;
			if(Point(poly->PointNr(j)).X() >= Point(poly->InnerEdge().Last()).X())	// x-coordinates are increasing
			{
				poly->InnerEdge().Add(poly->PointNr(j));
			}
			else
			{			
				lastInner = j-1;
				break;
			}
		}

		if(nrMinX!=1)
		{
			for(int j=1;j<=nrMinX-1;j++)	//from begin of polygon to starting point
			{
				if(Point(poly->PointNr(j)).X() >= Point(poly->InnerEdge().Last()).X())	// x-coordinates are increasing
				{
					poly->InnerEdge().Add(poly->PointNr(j));
					lastInner = j;
				}
				else
				{	
					break;
				}
			}
		}
		 //all the points of the inner edge are now stored ascending

//$ YV 2013-01-12: no SolverUO here; another solution is needed if necessary
//		if(GetMBS()->SolverUO().GetGlobalMessageLevel() >= UO_LVL_dbg1)
		{
			GetMBS()->UO(UO_LVL_ext) << "Polygons("<<i<<"): nrMinX="<<nrMinX <<", lastInner = "<<lastInner<<"\n";
		}

		// store the points of the outer edge
		if(lastInner < nrMinX)
		{
			for(int j=nrMinX;j>=lastInner;j--)
			{
					poly->OuterEdge().Add(poly->PointNr(j));
			}
		}
		else
		{
			for(int j=nrMinX;j>=1;j--)
			{
					poly->OuterEdge().Add(poly->PointNr(j));
			}
			for(int j=poly->NPoints()-1;j>=lastInner;j--)
			{
					poly->OuterEdge().Add(poly->PointNr(j));
			}
		}

	}

}

// search highest domainnumber in all lines
int Crosssection::HighestDomainNumber()
{
	int highest = 0;
	for(int i=1; i<= Lines().Length(); i++)
	{
		if(Line(i).Domain1() > highest) highest = Line(i).Domain1();
		if(Line(i).Domain2() > highest) highest = Line(i).Domain2();
	}
	return highest;
}

// finds all LineSegments that border a domain, use these to find closed polygons
int Crosssection::FilterLineSegmentsForDomainNumber(int domainnumber, TArray<int>& domain_on_left_side, TArray<int>& domain_on_right_side, TArray<int>& domain_on_both_sides)
{
	domain_on_left_side.Flush();
	domain_on_right_side.Flush();
	domain_on_both_sides.Flush();

	for(int i=1; i<=Lines().Length(); i++)
	{
		if(Line(i).Domain1() == domainnumber && Line(i).Domain2() != domainnumber) domain_on_left_side.Add(i);
		if(Line(i).Domain1() != domainnumber && Line(i).Domain2() == domainnumber) domain_on_right_side.Add(i);
		if(Line(i).Domain1() == domainnumber && Line(i).Domain2() == domainnumber) domain_on_both_sides.Add(i);
	}

	return domain_on_left_side.Length() + domain_on_right_side.Length();
}

// Find a continuation for the open polygon, returns number of next line
int Crosssection::FindContinuation(int lastpoint, TArray<int>& candidates, int flag_leftright)
{
	for(int i=1; i<=candidates.Length(); i++)
	{
		if(flag_leftright == 1 && lastpoint == Line(candidates(i)).P1()) return i;
		if(flag_leftright == 2 && lastpoint == Line(candidates(i)).P2()) return i;
	}
	return 0;
}

int Crosssection::ResortPolygonCCW(int i)
{
	ClosedPolygon& the_polygon = Polygon(i);
	int domain = the_polygon.Domain();
// look at order of two points in respect to first line segment for ALL lines up to the first line where left and right border are
	int counter = 1;
	while( Line(the_polygon.LineNr(counter)).Domain1() == Line(the_polygon.LineNr(counter)).Domain2() )
		counter++;

	LineSegment& the_line = Line(the_polygon.LineNr(counter));
	int lp1nr = the_line.P1();
	int lp2nr = the_line.P2();

// cw -> resort
	if( ( lp1nr == the_polygon.PointNr(counter) && lp2nr == the_polygon.PointNr(counter+1) && the_line.Domain2() == the_polygon.Domain() ) ||
		  ( lp1nr == the_polygon.PointNr(counter+1) && lp2nr == the_polygon.PointNr(counter) && the_line.Domain1() == the_polygon.Domain() ) )
	{
		TArray<int> copy_of_pointsequence(the_polygon.PointNrs());
		TArray<int> copy_of_linesequence(the_polygon.LineNrs());

		for(int i=1; i<= the_polygon.NPoints(); i++)
		{
			the_polygon.PointNr(i) = copy_of_pointsequence(the_polygon.NPoints()+1-i);
		}
		for(int i=1; i<= the_polygon.NLines(); i++)
		{
			the_polygon.LineNr(i) = copy_of_linesequence(the_polygon.NLines()+1-i);
		}
	}
	else
	{
	  ;
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Additional Point - splits Linesegment in two - updates Polygons
int Crosssection::AddPointOnLine(Vector2D splitpoint, int linesegmentnumber)
{
// check if point is valid: *point does not exist *point is ON a linesegment * linesegmentnumber is correct
	int splitpointnumber = Points().Find(splitpoint);
	if(splitpointnumber) return splitpointnumber; // new point is corner -> done
  if(linesegmentnumber == 0) linesegmentnumber = FindLineSegment(splitpoint); // no segmentnumber defined, find
	if(linesegmentnumber == 0) return 0; // new point is not on any linesegment ( after search )
	if(!IsOnLineSegment(splitpoint, linesegmentnumber)) return 0; // point is not on specified
	
	splitpointnumber = Points().Add(splitpoint); // Add the Point

	int newlinenumber = SplitLineSegment(linesegmentnumber,splitpointnumber);// split the line in two

// close the polygons again
// in each of the affected polygons the loop is opened - the linesegment is aready the shortened, with the new point on one side
// 
  for(int i=1; i<=Polygons().Length(); i++)
	{
	  ClosedPolygon& the_polygon = Polygon(i);

		int ls_idx = the_polygon.LineNrs().Find(linesegmentnumber);
		if(ls_idx) // this polygon has to be changed
		{
			StitchPolygon(i,ls_idx, newlinenumber, splitpointnumber);
		}
	}
	return splitpointnumber;
}

// Find the Segment the point is on
int Crosssection::FindLineSegment(Vector2D& splitpoint)
{
	for(int i=1; i<=Lines().Length(); i++)
	{
		if(IsOnLineSegment(splitpoint, i))
			return i;
	}
  return 0;
}

// checks if the point is on the given line segment
int Crosssection::IsOnLineSegment(Vector2D& splitpoint, int linenr)
{
// linear
	Vector2D p1 = Point( Line(linenr).P1() );
	Vector2D p2 = Point( Line(linenr).P2() );
	Vector2D v12 = p2-p1;
  
	if( (v12.X() == 0.) )
	{	// horizontal line
	  if (abs(p1.X()-splitpoint.X())>MESH_STD_TOL) return 0; // off a vertical line
		double dy = (splitpoint.Y() - p1.Y()) / v12.Y();
		if( (dy<=0.) || (dy>=1.) ) return 0;
		return linenr;
	}

	if( (v12.Y() == 0.) )
	{ // vertical line 
		if(abs(p1.Y()-splitpoint.Y())>MESH_STD_TOL) return 0; // off a horizontal line
		double dx = (splitpoint.X() - p1.X()) / v12.X();
		if( (dx<=0.) || (dx>=1.) ) return 0;
		return linenr;
	}

	double dx = (splitpoint.X() - p1.X()) / v12.X();
	double dy = (splitpoint.Y() - p1.Y()) / v12.Y();

	if( abs(dx-dy)>MESH_STD_TOL ) return 0; 
	if( (dx<=0.) || (dx>=1.) ) return 0;
	if( (dy<=0.) || (dy>=1.) ) return 0;
	return linenr;

	// TODO: SPLINE
}

// splits one linesegment into two at the splitpoint
int Crosssection::SplitLineSegment(int linesegmentnumber, int splitpointnumber)
{
// ALWAYS ALWAYS ALWAYS ( function ClosedPolygon::Stitch reqires that ) ALWAYS ALWAYS ALWAYS
// old line P1 -> Pnew
// new line Pnew -> P2
	LineSegment& oldline(Line(linesegmentnumber));
	LineSegment newline(oldline);
	oldline.P2() = splitpointnumber;
	newline.P1() = splitpointnumber;
	newline.P3() = 0; // new lines are always LINEAR, not spline
	return Lines().Add(newline);
}

// closes hole in polygon
int Crosssection::StitchPolygon(int polygonnr, int line_index, int newlinenumber, int splitpointnumber)
{
  ClosedPolygon& the_polygon = Polygon(polygonnr);
	int p1_nr = Line(the_polygon.LineNr(line_index)).P1();
	int p2_nr = Line(the_polygon.LineNr(line_index)).P2();
	int p1_idx = the_polygon.PointNrs().Find(p1_nr);
	int p2_idx = the_polygon.PointNrs().Find(p2_nr);

	int nr_of_points = the_polygon.NPoints() - 1; // number of points in polygon for modulo operation. use -1 since first and last point are the same !
	
	if(p2_nr == splitpointnumber)
	{
		int oldendpoint_nr = Line(newlinenumber).P2();
		int oldendpoint_idx = the_polygon.PointNrs().Find(oldendpoint_nr);
		int diff = (oldendpoint_idx - p1_idx + nr_of_points) % nr_of_points; // valid results: +1 and -1
		if(diff == 1)
		{
			the_polygon.PointNrs().Insert(p1_idx+1, splitpointnumber);
			the_polygon.LineNrs().Insert(line_index+1, newlinenumber); 
		}
		else // diff==-1
		{
			the_polygon.PointNrs().Insert(p1_idx, splitpointnumber);
			the_polygon.LineNrs().Insert(line_index, newlinenumber); 
		}
	}
	else //(p1_nr == splitpointnumber)
	{
		int oldendpoint_nr = Line(newlinenumber).P2();
		int oldendpoint_idx = the_polygon.PointNrs().Find(oldendpoint_nr);
		int diff = (oldendpoint_idx - p2_idx + nr_of_points) % nr_of_points; // valid results: +1 and -1

		//if( (oldendpoint_idx-1)&nr_of_points > (p2_idx-1)%nr_of_points )
		if(diff == 1)
		{
			the_polygon.PointNrs().Insert(p2_idx+1, splitpointnumber);
			the_polygon.LineNrs().Insert(line_index+1, newlinenumber); 	
		}
		else // diff==-1
		{
			the_polygon.PointNrs().Insert(p2_idx, splitpointnumber);
			the_polygon.LineNrs().Insert(line_index, newlinenumber); 	
		}
	}
	return 0;
}

// split the segments into linear sections, if no segmentation is defined, use the "refinement" entry
int Crosssection::SplinesToLinear(int segments_manual) 
{
	int nlines = Lines().Length();
  for (int i=1; i<=nlines; i++)
	{
		if( Line(i).IsSpline() )
		{
			if( segments_manual <= 0) 
			{
				SplitSpLine(i, Line(i).Ref());
			}
			else
			{
				SplitSpLine(i, segments_manual);
			}
		}
	}
	return 0;
}

// split the line into equally long segments, checks if the line segment is a spline and moves the points accordingly
int Crosssection::SplitSpLine(int linesegmentnumber, int segments)
{
	LineSegment& the_line = Line(linesegmentnumber);
	Vector2D p1 = Point(the_line.P1());
	Vector2D p2 = Point(the_line.P2());
  Vector2D ps;// = (p1+p2)*0.5;
	if(the_line.IsSpline())
	{
		ps = Point(the_line.P3());
	}

	TArray<int> inter_point_nrs;
	TArray<Vector2D> inter_points_true_coords;
// split as if it were all linear
	for(int i=1; i<=segments-1; i++)
	{
		double interpolation_factor = 1.-(double)i/(double)segments; // work backwards to have same line segment number !!!
		Vector2D splitpoint_coord_linear = InterpolateLinear(interpolation_factor, p1, p2);
		int newpoint_nr = AddPointOnLine(splitpoint_coord_linear, linesegmentnumber); 
		inter_point_nrs.Add(newpoint_nr);
		
		if(the_line.IsSpline())
		{
			Vector2D splitpoint_coord_quadratic_bezier = InterpolateQuadraticBezier(interpolation_factor, p1, p2, ps);
			inter_points_true_coords.Add(splitpoint_coord_quadratic_bezier);
		}
		else
		{
			inter_points_true_coords.Add(splitpoint_coord_linear);
		}
	}
// correction of point coordinates
	for(int i=1; i<=inter_point_nrs.Length(); i++)
	{
		Point(inter_point_nrs(i)) = inter_points_true_coords(i);
	}
	return 0;
}

Vector2D Crosssection::InterpolateLinear(double factor, Vector2D& p1, Vector2D& p2)
{
	double& t = factor;
	return (1-t)*p1 + t*p2;
}

Vector2D Crosssection::InterpolateQuadraticBezier(double factor, Vector2D& p1, Vector2D& p2, Vector2D& ps)
{
	double& t = factor;
  return (1-t)*((1-t)*p1 + t*ps) + t*((1-t)*ps + t*p2);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Bud a new Polygon - additional linesegment, create 2 polygons
int Crosssection::BudPolygon(int p1_nr, int p2_nr)
{
// check if point is valid: *line_segment does not exist *EXACTLY ONE polygon where both points are in
	for (int i=1; i<=Lines().Length(); i++) 
	{ 
		if( ((Line(i).P1() == p1_nr) && (Line(i).P2() == p1_nr)) || ((Line(i).P1() == p2_nr) && (Line(i).P2() == p1_nr)) ) 
			return 0; // error: line exists
	}
	int polygonnr = 0;
	for (int i=1; i<=Polygons().Length(); i++)
	{
		if(Polygon(i).PointNrs().Find(p1_nr) && Polygon(i).PointNrs().Find(p2_nr))
		{
			if(polygonnr) return 0; // error: polygon not unique
			else polygonnr = i;
		}
	}

	ClosedPolygon& the_polygon = Polygon(polygonnr);
// make the additional linesegment
	LineSegment divisionline(mbs, p1_nr, p2_nr, the_polygon.Domain(), the_polygon.Domain(), 1, 1);
	int divisionlinenumber = Lines().Add(divisionline);

// make new Polygon
	TArray<int> budpoints, budlines;
	int p1_idx = the_polygon.PointNrs().Find(p1_nr);
	int p2_idx = the_polygon.PointNrs().Find(p2_nr);
	int min_idx = (p1_idx < p2_idx ? p1_idx : p2_idx);
	int max_idx = (p1_idx < p2_idx ? p2_idx : p1_idx);

// built arrays for new polygon
	for(int i=min_idx; i<=max_idx; i++)
	{
	  budpoints.Add(the_polygon.PointNr(i));
   	if(i!=max_idx)
		{
			budlines.Add(the_polygon.LineNr(i));
		}
		else
		{
			budlines.Add(divisionlinenumber);
		}
	}
	budpoints.Add(the_polygon.PointNr(min_idx));
	ClosedPolygon bud(mbs,budpoints, budlines);
	bud.Domain() = the_polygon.Domain();

// remove from original polygon
	for(int i=max_idx-1; i>=min_idx+1; i--)
	{
		the_polygon.PointNrs().Erase(i);
		the_polygon.LineNrs().Erase(i);
	}
	the_polygon.LineNrs().Erase(min_idx);
	the_polygon.LineNrs().Insert(min_idx, divisionlinenumber);
	
	return Polygons().Add(bud);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
#define iterative
#ifdef iterative

int Crosssection::CutAligned(Vector2D guide, double angular_variance_deg)
{
	int i1 = 0; //CutAligned_ContinueLines(guide,angular_variance_deg);
	int i2 = CutAligned_AtEachPoint(guide,angular_variance_deg);
	return i1+i2;
}

int Crosssection::CutAligned_ContinueLines(Vector2D guide, double angular_variance_deg)
{
	double angular_variance_rad = angular_variance_deg / 180. * MY_PI;
	guide.Normalize();
	int budcount = 0;

	// try to continue every line
	for(int i=1; i<=NLines(); i++)
	{
		int startline_nr = i;
		double ang_dev = TwoPointVectorAngleToReference(Line(i).P1(), Line(i).P2(), guide);
		if(abs(ang_dev) <= angular_variance_rad)
		{
			int nlines = NLines(); // dont cut more then once
			for(int j=1; j<=nlines; j++)
			{
				int targetline_nr = j;
				// CHECK I: no common point 
				if( CommonPointOfLineSegments(i,j) == 0 ) 
				{
					// dont try to cut a line that is within angular variance of guide
					double ang_dev_target = TwoPointVectorAngleToReference(Line(j).P1(), Line(j).P2(), guide);
					if(abs(ang_dev_target) > angular_variance_rad) 
					{				
						Vector2D lambda_mu = CutTwoLineSegments_LambdaMu(i,j);
						double mu = lambda_mu.Y();
						if(0.<=mu && mu<=1.)
						{
							// starting point for budline 
							double lambda = lambda_mu.X();
							int startpoint_nr;
							if(lambda<0.5)
								startpoint_nr = Line(i).P1();
							else
								startpoint_nr = Line(i).P2();

							// CHECK II: the "nearest point" and the target line are in the same polygon
							int commonpolygon = CommonPolygonOfPointAndLineSegment(startpoint_nr,j);
							if(commonpolygon)
							{
								ClosedPolygon& the_polygon = Polygon(commonpolygon);	

								// Compute the end point of the line
								Vector2D endpoint_coords;
								int endpoint = GetBudlineEndPointCoords(startpoint_nr, targetline_nr, guide, endpoint_coords, angular_variance_rad);

								// CHECK III: budline must be an inner line of the polygon
								if(IsInnerLineOfPolygon(commonpolygon, startpoint_nr, targetline_nr, endpoint_coords)) 
								{
									if(endpoint == 0)
									{
										AddPointOnLine(endpoint_coords, targetline_nr);
										endpoint = Points().Length();
									}		
									// CHECK IV: line does not exist yet
									if(!LineSegmentExistsInPolygon(startpoint_nr, endpoint, commonpolygon))
									{
										BudPolygon(startpoint_nr, endpoint);
										budcount++;
									}
								}
							}
						}
					}
				}
			}
		}	
	}
	return budcount;
}

int Crosssection::CutAligned_AtEachPoint(Vector2D guide, double angular_variance_deg)
{
	double angular_variance_rad = angular_variance_deg / 180. * MY_PI;
	guide.Normalize();
	int budcount = 0;

	// try to continue every point, even new ones... 
	int npoints = NPoints();
	for(int i=1; i<=NPoints(); i++)
	{
		int startpoint_nr = i;
		int nlines = NLines(); // dont cut more then once
		for(int j=1; j<=nlines; j++)
		{
			int targetline_nr = j;
			// CHECK I: no common point 
			if( Line(j).P1() != i && Line(j).P2() != i ) 
			{
				// dont try to cut a line that is within angular variance of guide
				double ang_dev_target = TwoPointVectorAngleToReference(Line(j).P1(), Line(j).P2(), guide);
				if(abs(ang_dev_target) > angular_variance_rad) 
				{
					Vector2D lambda_mu = CutStraightLineWithLineSegment_LambdaMu(startpoint_nr, guide, targetline_nr);
					double mu = lambda_mu.Y();
					if(0.<=mu && mu<=1.)
					{
						// CHECK II: the "nearest point" and the target line are in the same polygon
						int commonpolygon = CommonPolygonOfPointAndLineSegment(startpoint_nr, targetline_nr);
						if(commonpolygon)
						{
							ClosedPolygon& the_polygon = Polygon(commonpolygon);	

							// HACK AD - manaully filtering the lines not "inside"
							//if( startpoint_nr == 4 || startpoint_nr == 8 || startpoint_nr == 11 || startpoint_nr == 15 || startpoint_nr == 20 || startpoint_nr == 22 )
							//	continue;
							// END HACKS
						
							// Compute the end point of the line
							Vector2D endpoint_coords;
							int endpoint_nr = GetBudlineEndPointCoords(startpoint_nr, targetline_nr, guide, endpoint_coords, angular_variance_rad);

							// CHECK III: budline must be an inner line of the polygon
							int innerlineflag = IsInnerLineOfPolygon(commonpolygon, startpoint_nr, targetline_nr, endpoint_coords);
							if(innerlineflag) 
							{
								if(endpoint_nr == 0)
								{
									AddPointOnLine(endpoint_coords, targetline_nr);
									endpoint_nr = Points().Length();
								}		
								// CHECK IV: line does not exist yet
								if(!LineSegmentExistsInPolygon(startpoint_nr, endpoint_nr, commonpolygon))
								{
									BudPolygon(startpoint_nr, endpoint_nr);
									budcount++;
								}
							}
						}				
					}
				}
			}
		}	
	}
	return budcount;
}
#else
int Crosssection::CutAligned(Vector2D guide, double angular_variance_deg)
{
  double angular_variance_rad = angular_variance_deg / 180. * MY_PI;
	guide.Normalize();
	Vector2D n_guide(-guide.Y(), guide.X());
	
	TArray<int2> budlines;
	for(int i=1; i<=Lines().Length(); i++)
	{
    double ang_dev = TwoPointVectorAngleToReference(Line(i).P1(), Line(i).P2(), guide);
		if(abs(ang_dev) <= angular_variance_rad)
		{
// only search if lines of higher index are cut ...
			int nlines = Lines().Length(); // dont cut more then once
			for(int j=1; j<=nlines; j++)
			{
// no common point
				if( (Line(i).P1() != Line(j).P1()) && (Line(i).P1() != Line(j).P2()) && (Line(i).P2() != Line(j).P1()) && (Line(i).P2() != Line(j).P2()) ) 
				{
// THIS HAS TO BE CHANGED IN FINAL VERSION
// try only to cut with lines that are nearly normal to vector guide
					double ang_dev_b = TwoPointVectorAngleToReference(Line(j).P1(), Line(j).P2(), n_guide);
					if(abs(ang_dev_b) <= angular_variance_rad) 
					{
						double lambda = CutTwoLineSegments_Lambda(i,j);
						if(0.<=lambda && lambda<=1.)
						{
							double nu = CutTwoLineSegments_Lambda(j,i);
							int nearestpoint;
							if(nu<0.5)
								nearestpoint = Line(i).P1();
							else
								nearestpoint = Line(i).P2();
// test if angle to linesegment endpoints would be valid
							double ang_dev_p1 = TwoPointVectorAngleToReference(Line(j).P1(),nearestpoint, guide);
							double ang_dev_p2 = TwoPointVectorAngleToReference(Line(j).P2(),nearestpoint, guide);

							if( (abs(ang_dev_p1)<=abs(ang_dev_p2)) && (abs(ang_dev_p1)<=angular_variance_rad) )
							{
								budlines.Add(int2(nearestpoint,Line(j).P1())); // between two existing points
							}
							else if( (abs(ang_dev_p2)<=abs(ang_dev_p1)) && (abs(ang_dev_p2)<=angular_variance_rad) )
							{
								budlines.Add(int2(nearestpoint,Line(j).P2())); // between two existing points
							}
							else
							{
								Vector2D p1j = Point(Line(j).P1());
								Vector2D p2j = Point(Line(j).P2());
							  Vector2D newpoint = p1j*(1-lambda) + p2j*lambda;
								
								// make a new point 
								AddPointOnLine(newpoint,j);
								budlines.Add(int2(nearestpoint,Points().Length()));
							}
						}
					}
				}
			}
		}
	}
// erase double lines
	for(int i=1; i<=budlines.Length(); i++)
	{
		for(int j=i+1; j<=budlines.Length(); j++)
		{
			if( (budlines(i)(1) == budlines(j)(2) && budlines(i)(2) == budlines(j)(1)) || (budlines(i)(1) == budlines(j)(1) && budlines(i)(2) == budlines(j)(2)) )
			{
				budlines.Erase(j);
			}
		}
	}
// erase unnecessary lines ( connected through other paths )

	for(int i=1; i<=budlines.Length(); i++)
	{
		BudPolygon(budlines(i)(1), budlines(i)(2));
	}
	return 0;
}
#endif

int Crosssection::IsInnerLineOfPolygon(int polygonnr, int startingpointnr, int targetlinenr, Vector2D endpoint_coord)
{
  ClosedPolygon& the_polygon = Polygon(polygonnr);
  int nr_of_polygon_points = the_polygon.NPoints()-1;
	Vector2D startinpoint_coord = Point(startingpointnr);
	Vector2D budlinevector = endpoint_coord - startinpoint_coord;

// step I: determine if the new line starts as inside inside --> angle from previous borderline to next borderline is larger than angle from previous borderline to guide-vector
// not an inside line if angle to guide- >= angle next borderline
	int startingpoint_idx = the_polygon.PointNrs().Find(startingpointnr);

	int prevpoint_idx = startingpoint_idx-1;
	if(prevpoint_idx == 0) 
		prevpoint_idx = nr_of_polygon_points;
  int prevpointnr = the_polygon.PointNr(prevpoint_idx);

	int nextpoint_idx = startingpoint_idx+1;
	if(nextpoint_idx == nr_of_polygon_points +1) 
		nextpoint_idx = 1;
  int nextpointnr = the_polygon.PointNr(nextpoint_idx);

// TODO: check if the polygon is in clockwise or counterclockwise
	int flag_clockwise = 0;



	Vector2D minus_v_to =  Point(prevpointnr) - Point(startingpointnr);
	Vector2D v_away = Point(nextpointnr) - Point(startingpointnr);

	double polar_v_to = PolarAngle(minus_v_to);
	double polar_v_away = PolarAngle(v_away);
	double polar_v_budline = PolarAngle(budlinevector);

	double angle_border = polar_v_away - polar_v_to;
	while (angle_border<0.) angle_border += 2.*MY_PI;
	double angle_budline = polar_v_budline - polar_v_to;
	while (angle_budline<=0.) angle_budline += 2.*MY_PI;

	if(  flag_clockwise && angle_budline >= angle_border ) return 0; // <== DOES NOT START AS INSIDE LINE 
 	if( !flag_clockwise && angle_budline <= angle_border ) return 0; // <== DOES NOT START AS INSIDE LINE 
 
// step I: determine if no other line cuts the budline...
	for(int i=1; i <= NLines(); i++)
	{
		Vector2D lambda_mu = CutStraightLineWithLineSegment_LambdaMu(startingpointnr, budlinevector, targetlinenr);
		double lambda = lambda_mu.X();
		if( (0.+MESH_STD_TOL) < lambda && lambda < (1.-MESH_STD_TOL) )
		{
		  return 0;
		}
	}

  return 1;
}

int Crosssection::GetBudlineEndPointCoords(int startpointnr, int targetlinenr, Vector2D guide, Vector2D& rv_endpoint_coords, double angular_variance_rad)
{
	int endpoint = 0;
	Vector2D lambda_mu = CutStraightLineWithLineSegment_LambdaMu(startpointnr, guide, targetlinenr);
	double lambda = lambda_mu.X();
	double mu = lambda_mu.Y();
	double ang_dev_p1 = TwoPointVectorAngleToReference(startpointnr, Line(targetlinenr).P1(), guide);
	double ang_dev_p2 = TwoPointVectorAngleToReference(startpointnr, Line(targetlinenr).P2(), guide);

	if( abs(ang_dev_p1)<=angular_variance_rad && abs(ang_dev_p2)<=angular_variance_rad) // both points in other line are valid destinations, pick the nearer one
	{
		if(mu<0.5) 
			endpoint = Line(targetlinenr).P1();
		else
			endpoint = Line(targetlinenr).P2();
	}
	else if( abs(ang_dev_p1)<=angular_variance_rad ) // only p1 in angle range
	{
		endpoint = Line(targetlinenr).P1(); 
	}
	else if( abs(ang_dev_p2)<=angular_variance_rad ) // only p2 in angle range
	{
		endpoint = Line(targetlinenr).P2(); 
	}
	else 	// compute coordinates of the new point
	{
		Vector2D p1j = Point(Line(targetlinenr).P1());
		Vector2D p2j = Point(Line(targetlinenr).P2());
		rv_endpoint_coords = p1j*(1-mu) + p2j*mu;	
	}
	if(endpoint != 0)
		rv_endpoint_coords = Point(endpoint);
	return endpoint;
}

int Crosssection::LineSegmentExistsInPolygon(int startpoint, int endpoint, int polygonnumber)
{
	int existingline = 0;
	ClosedPolygon& the_polygon = Polygon(polygonnumber);
	for(int i=1; i<=the_polygon.NLines(); i++)
	{
		LineSegment& the_line = Line(the_polygon.LineNr(i));
		if( (the_line.P1() == startpoint && the_line.P2() == endpoint) ||
				(the_line.P1() == endpoint && the_line.P2() == startpoint) )
		{
			existingline = i;
		}
	}
  return existingline;
}

// returns the angle in rad of a vector specified by 2 points to the specified reference
double Crosssection::TwoPointVectorAngleToReference(int p1_idx, int p2_idx, Vector2D reference)
{	
	reference.Normalize();
	Vector2D p1 = Point(p1_idx);
	Vector2D p2 = Point(p2_idx);
	Vector2D v = p2-p1;
  double ang_dev = acos( abs(reference * v / v.Norm()) );

	return ang_dev;
}

// returns the angle in rad of the linesegment to the specified reference
double Crosssection::LineSegmentAngleToReference(int linesegmentnumber, Vector2D reference)
{
	return TwoPointVectorAngleToReference(Line(linesegmentnumber).P1(), Line(linesegmentnumber).P2(), reference);
}

// returns the parameter lambda of the intersection point on the second linesegment
Vector2D Crosssection::CutStraightLineWithLineSegment_LambdaMu(int pointnumber, Vector2D guide, int linesegmentnumber)
{
	Vector2D p1 = Point(pointnumber);
	Vector2D v1 = guide;
	Vector2D p2 = Point(Line(linesegmentnumber).P1());
	Vector2D v2 = Point(Line(linesegmentnumber).P2())-p2;

  Vector2D pi, lambda_mu;
	IntersectLines2D(p1,v1, p2,v2, pi, lambda_mu);
	return lambda_mu;
}

Vector2D Crosssection::CutTwoLineSegments_LambdaMu(int linesegmentnumber1, int linesegmentnumber2)
{
	int pt_nr = Line(linesegmentnumber1).P1();
	Vector2D p1 = Point(Line(linesegmentnumber1).P1());
	Vector2D v1 = Point(Line(linesegmentnumber1).P2())-p1;
	return CutStraightLineWithLineSegment_LambdaMu(pt_nr, v1, linesegmentnumber2);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// assign values for MBS::Draw
int Crosssection::AssingDrawPointsAndColorToPolygon(int polygonnumber)
{
	TArray<Vector2D> loccoords;
	ClosedPolygon& the_polygon = Polygon(polygonnumber);
	for(int i=1; i< the_polygon.NPoints(); i++)
	{
		loccoords.Add(points(the_polygon.PointNr(i)));
	}

	the_polygon.SetDrawPoints(loccoords);
	int color_nr = (the_polygon.Domain()-1) % (mycolors.N()-1) +2; // nonblack filling color
	the_polygon.SetDrawColor(&mycolors,color_nr);
	the_polygon.SetName(mystr("ClosedPolygon ")+mystr(polygonnumber));
	
	for(int i=1; i<=the_polygon.NLines(); i++)
	{
		LineSegment& the_line = Line(the_polygon.LineNr(i));
		Vector2D p1 = Point(the_line.P1());
		Vector2D p2 = Point(the_line.P2());
		Vector2D ps(0.,0.);
		if(the_line.IsSpline())
		{
			ps = Point(the_line.P3());
		}
		the_line.SetDrawPoints(p1,p2,ps);
		the_line.SetDrawColor(&mycolors); // boundary condition depending line color
		the_line.SetName(mystr("Polygon ")+mystr(polygonnumber)+mystr(" Edge ")+mystr(i));
	}
	
	return loccoords.Length();
}

// export the edges to create the 3D rotatory geometry in Netgen
// filling the arrays leftpoints, rightpoints, materials and refinements based on the information stored in InnerEdge, OuterEdge and Point of the polygons of the crosssection
// ATTENTION: the resulting coordinates are not in the xy-plane anymore! They are stored as [y,0,x] in leftpoints and rightpoints, this corresponds to the coordinates of a point on the rotor: [r,0,ax_pos]
//						leftpoints(i) and rightpoints(i) for i = 1,3,5... represent the outer edge of the geometry
//						leftpoints(i) and rightpoints(i) for i = 2,4,6... represent the inner edge of the geometry
int Crosssection::ExportEdgesForNetgen(TArray<Vector3D>& leftpoints, TArray<Vector3D>& rightpoints, TArray<int>& materials, TArray<int>& refinements, int draw_edges)
{
	leftpoints.Flush();
	rightpoints.Flush();
	materials.Flush();
	refinements.Flush();

	int p1,p2;
	int p1I,p2I;
	int maxI, maxO;
	TArray<double> xValues, nodeNr;
	TArray<int> flagIO;
	int count=1;
	int countIn=1;
	int countOut=1;
	int found=1;
	ClosedPolygon* the_polygon;
	double x,y;
	int next_edge;
	//int polyNr =4;
	//for(int i=polyNr; i<=polyNr; i++)
	for(int i=1; i<=Polygons().Length(); i++)
	{
		xValues.Flush();
		nodeNr.Flush();
		flagIO.Flush();

		the_polygon = &Polygon(i);
		maxI = the_polygon->InnerEdge().Length();
		maxO = the_polygon->OuterEdge().Length();

		for(int j=1; j<=maxI; j++)
		{
			p1 = the_polygon->PointNrInnerEdge(j);
			xValues.Add(Point(p1).X());
			nodeNr.Add(p1);
			flagIO.Add(0);	//0..InnerEdge
		}
		for(int j=1; j<=maxO; j++)
		{
			p1 = the_polygon->PointNrOuterEdge(j);
			xValues.Add(Point(p1).X());
			nodeNr.Add(p1);
			flagIO.Add(1);	//1..OuterEdge
		}
		// all x-coordinates are added to the TArrays
		QuicksortDouble(xValues,nodeNr,flagIO);

		//----------------------------------
		//----------------------------------

		if(draw_edges)
		{
			for(int j=1; j<=maxI-1; j++)
			{
				p1 = the_polygon->PointNrInnerEdge(j);
				p2 = the_polygon->PointNrInnerEdge(j+1);
				LineSegment line(GetMBS(),p1,p2,1,1,1,1);
				line.DrawStyle() = 16;
				line.SetDrawPoints(Point(p1),Point(p2));
				line.SetDrawParam(Vector3D(1.,0.05,0.02));
				line.SetCol(Vector3D(1,0,0));
				GetMBS()->Add(line); // Add GeomLine
				GetMBS()->UO(UO_LVL_dbg1) << "Polygons("<<i<<"): inner edge: line p1= "<< p1<<" = "<<Point(p1)<<", p2="<< p2 <<" = "<<Point(p2)<<"\n";
			}

			for(int j=1; j<=maxO-1; j++)
			{
				p1 = the_polygon->PointNrOuterEdge(j);
				p2 = the_polygon->PointNrOuterEdge(j+1);

				LineSegment line(GetMBS(),p1,p2,1,1,1,1);
				line.DrawStyle() = 16;
				line.SetDrawPoints(Point(p1),Point(p2));
				line.SetDrawParam(Vector3D(1.,0.05,0.02));
				line.SetCol(Vector3D(0,1,0));
				GetMBS()->Add(line); // Add GeomLine
				GetMBS()->UO(UO_LVL_dbg1) << "Polygons("<<i<<"): outer edge: line p1= "<< p1 <<" = "<<Point(p1)<<", p2="<< p2<<" = "<<Point(p2)<<"\n";	
			}
			// end of draw

			GetMBS()->UO(UO_LVL_dbg1) << "Polygons("<<i<<"): \n";
			for(int j=1;j<=flagIO.Length();j++)
			{
				GetMBS()->UO(UO_LVL_dbg1) << "j:" <<j<<", x= "<< xValues(j) <<", node="<< nodeNr(j) <<", flag = " << flagIO(j)<<"\n";	
			}
		} 
		//----------------------------------
		//----------------------------------

		// first element
		p1 = the_polygon->PointNrOuterEdge(1);
		p2 = the_polygon->PointNrOuterEdge(2);
		if(Point(p1).X()==Point(p2).X())		//starting with vertical edge
		{
			x=Point(p2).X();
			int k=2;
			while((Point(the_polygon->PointNrOuterEdge(k)).X()==x)&&(k<the_polygon->OuterEdge().Length())) {k++;}		// get the point with the highest y-coord
			leftpoints.Add(Vector3D(Point(the_polygon->PointNrOuterEdge(k-1)).Y(),0.,x));
			leftpoints.Add(Vector3D(Point(p1).Y(),0.,x));
			countIn=2;
			countOut=k-1;	//1
			next_edge=1;
		}
		else
		{
			// the beginning of the polygon is a triangle
			// this triangle is cut away
			GetMBS()->UO(UO_LVL_warn) << "WARNING: Crosssection::ExportEdgesForNetgen: Polygon("<<i<<") starts with a triangle, this is not implemented, so it is cut.\n";
			p1 = the_polygon->PointNrOuterEdge(1);
			p2 = the_polygon->PointNrOuterEdge(2);
			p1I=the_polygon->PointNrInnerEdge(1);
			p2I = the_polygon->PointNrInnerEdge(2);

			if(Point(p2).X()<Point(p2I).X())
			{
				x=Point(p1).X() + 0.5*(Point(p2).X()-Point(p1).X());
				y=Point(p1).Y() + 0.5*(Point(p2).Y()-Point(p1).Y());
				leftpoints.Add(Vector3D(y,0.,x));

				y=Point(p1I).Y()+(Point(p2I).Y()-Point(p1I).Y())*(x-Point(p1I).X())/(Point(p2I).X()-Point(p1I).X());
				leftpoints.Add(Vector3D(y,0.,x));
				countIn=1;
				countOut=2;
				next_edge=1;
			}
			else
			{
				x=Point(p1I).X()+0.5*(Point(p2I).X()-Point(p1I).X());

				y=Point(p1).Y()+(Point(p2).Y()-Point(p1).Y())*(x-Point(p1).X())/(Point(p2).X()-Point(p1).X());
				leftpoints.Add(Vector3D(y,0.,x));
				
				y=Point(p1I).Y() + 0.5*(Point(p2I).Y()-Point(p1I).Y());
				leftpoints.Add(Vector3D(y,0.,x));

				countIn=2;
				countOut=1;
				next_edge=0;
			}
		}
		
		while((countIn<=maxI)&&(countOut<=maxO))
		{
			// find next x-value for rightpoints
			x=leftpoints.Last().Z();
			for(int j=1;j<=xValues.Length();j++)
			{
				if(x < xValues(j))
				{
					next_edge = flagIO(j);
					count=j;
					j=xValues.Length()+1;			//to break for loop
				}
			}

			// set correct counter for next right point
			if(next_edge)
			{
				for(int j=1;j<=the_polygon->OuterEdge().Length();j++)
				{
					if(nodeNr(count)==the_polygon->PointNrOuterEdge(j))
					{
						countOut = j;
						j=the_polygon->OuterEdge().Length()+1;			//to break for loop
					}
				}
				for(int j=countOut-1;j>=1;j--)		//look for first node on outer edge that has the correct x-coord
				{
					if(Point(the_polygon->PointNrOuterEdge(j)).X()!=Point(the_polygon->PointNrOuterEdge(countOut)).X())
					{
						countOut=j+1;
						j=0;		// to break the for loop
					}
				}
			}
			else
			{
				for(int j=1;j<=the_polygon->InnerEdge().Length();j++)
				{
					if(nodeNr(count)==the_polygon->PointNrInnerEdge(j))
					{
						countIn = j;
						j=the_polygon->InnerEdge().Length()+1;			//to break for loop
					}
				}
				for(int j=countIn-1;j>=1;j--)		//look for first node on inner edge that has the correct x-coord
				{
					if(Point(the_polygon->PointNrInnerEdge(j)).X()!=Point(the_polygon->PointNrInnerEdge(countIn)).X())
					{
						countIn=j+1;
						j=0;		// to break the for loop
					}
				}
			}

			materials.Add(the_polygon->Domain());
			//materials.Add(the_polygon->Domain());
			materials.Add(0);
			refinements.Add(1); // dummy
			refinements.Add(1); // dummy

			//////// look for right end ////////////

			if(next_edge)	//outer edge comes first
			{
				p2 = the_polygon->PointNrOuterEdge(countOut);	// get the next node on outer surface
				x=Point(p2).X();
				rightpoints.Add(Vector3D(Point(p2).Y(),0.,x));		// outer edge

				for(int j=1; j<=maxI; j++)
				{
					p2I=the_polygon->PointNrInnerEdge(j);
					found = 1;
					if(Point(p2I).X()>=x)
					{
						found = j;
						j=maxI+1;	// to break for-loop
					}
				} // p2I and found correspond to point on inner edge

				if(Point(p2I).X()==x)
				{
					rightpoints.Add(Vector3D(Point(p2I).Y(),0.,x));	// inner edge
				}
				else
				{	// interpolation
					//p2I=the_polygon->PointNrInnerEdge(found);	//still valid
					p1I=the_polygon->PointNrInnerEdge(found-1);
					y=Point(p1I).Y()+(Point(p2I).Y()-Point(p1I).Y())*(x-Point(p1I).X())/(Point(p2I).X()-Point(p1I).X());
					rightpoints.Add(Vector3D(y,0.,x));	// inner edge
				}
			}
			else	// inner edge has smaller x-coord
			{
				p2I = the_polygon->PointNrInnerEdge(countIn);	// get the next node on inner surface
				x=Point(p2I).X();

				for(int j=1; j<=maxO; j++)
				{
					p2=the_polygon->PointNrOuterEdge(j);
					found = 1;
					if(Point(p2).X()>=x)
					{
						found = j;
						j=maxO+1;	// to break for-loop
					}
				} // p2O and found correspond to point on outer edge

				if(Point(p2).X()==x)
				{
					rightpoints.Add(Vector3D(Point(p2).Y(),0.,x));	// outer edge
				}
				else
				{	// interpolation
					//p2=the_polygon->PointNrOuterEdge(found);	//still valid
					p1=the_polygon->PointNrOuterEdge(found-1);
					y=Point(p1).Y()+(Point(p2).Y()-Point(p1).Y())*(x-Point(p1).X())/(Point(p2).X()-Point(p1).X());
					rightpoints.Add(Vector3D(y,0.,x));	// outer edge
				}
				rightpoints.Add(Vector3D(Point(p2I).Y(),0.,x));	// inner edge
			}// end of adding rightpoints

			// leftpoints of next element
			if(rightpoints.Last().Z()<xValues.Last())
			{
				//if(Point(nodeNr(count+1)).X()==x)	// same x-coord either on outer or inner edge again
				if((Point(the_polygon->PointNrOuterEdge(countOut+1)).X()==x)||((Point(the_polygon->PointNrInnerEdge(countIn+1)).X()==x)))
				{
					double yIn	= rightpoints.Last().X();		//X of rightpoints equal to Y of Vector2D
					double yOut = rightpoints(rightpoints.Length()-1).X(); //X of rightpoints equal to Y of Vector2D

					while(((Point(the_polygon->PointNrOuterEdge(countOut+1)).X()==x)||((Point(the_polygon->PointNrInnerEdge(countIn+1)).X()==x)))&&(count+1<flagIO.Length())) 
					{
						if((Point(the_polygon->PointNrOuterEdge(countOut+1)).X()==x)) //next outer edge has same x-coord
						{
							yOut =Point(the_polygon->PointNrOuterEdge(countOut+1)).Y();
							countOut++;
						}
						else	//next inner edge has same x-coord
						{
							yIn =Point(the_polygon->PointNrInnerEdge(countIn+1)).Y();
							countIn++;
						}
					}

					leftpoints.Add(Vector3D(yOut,0,x));
					leftpoints.Add(Vector3D(yIn,0,x));
				}
				else
				{	// actual rightpoints are next leftpoints
					leftpoints.Add(rightpoints(rightpoints.Length()-1));
					leftpoints.Add(rightpoints.Last());
				}
			}
			else
			{
				countIn = maxI+1;	// to break the while loop
			}
		}// while loop

		// clean away bad segments
		for(int j=1; j<=leftpoints.Length()-1; j=j+2)
		{
			if((leftpoints(j).Z()==rightpoints(j).Z())||(leftpoints(j+1).Z()==rightpoints(j+1).Z()))
			{
				leftpoints.Erase(j);	leftpoints.Erase(j);
				rightpoints.Erase(j);	rightpoints.Erase(j);
				materials.Erase(j);		materials.Erase(j);
				refinements.Erase(j);	refinements.Erase(j);
				GetMBS()->UO(UO_LVL_warn) << "WARNING: Crosssection::ExportEdgesForNetgen: Polygon("<<i<<"): ERASED because of bad edges!\n";
			}
		}
	}// for all polys
	return 1;
	// not yet implemented:
//					refinements.Add(the_line.Ref());
}

// simplification of the geometry: (cutted) cones are converted to (cutted) cylinders
// sets the radial component of the right points equal to those of the left points 
// can be called after ExportEdgesForNetgen, if there are problems e.g. with the meshing 
// the values in leftpoints and rightpoints are assumed to be coordinates of points on the rotor: [r,0,ax_pos]
int Crosssection::SimplifyConesToCylinders(TArray<Vector3D>& leftpoints, TArray<Vector3D>& rightpoints)
{
	for(int j=1; j<=leftpoints.Length()-1; j=j+2)
	{
		// outer edge
		if(leftpoints(j).X()!=rightpoints(j).X())
		{
			rightpoints(j).X()=leftpoints(j).X();
		}
		
		// inner edge
		if(leftpoints(j+1).X()!=rightpoints(j+1).X())
		{
			rightpoints(j+1).X()=leftpoints(j+1).X();
		}
	}
	return 1;
}

// determines if the linesegment is an "outer" edge of the polygon
int Crosssection::IsOuterEdge(int polygonnumber, int linesegmentnumber)
{
	ClosedPolygon& the_polygon = Polygon(polygonnumber);
  Vector2D p1 = Point(Line(the_polygon.LineNr(linesegmentnumber)).P1());
  Vector2D p2 = Point(Line(the_polygon.LineNr(linesegmentnumber)).P2());

	// * is not vertical
	if( abs((p2-p1).X()) <= MESH_STD_TOL ) return 0; // vertical line

// check if the value for y_bar of the line is the highest (from all non vertical..)
	double y_bar = (p1+p2).Y() * 0.5;
	TArray<double> axis_distance;
	for(int i=1; i<=the_polygon.NLines(); i++)
	{
		Vector2D lp1 = Point(Line(the_polygon.LineNr(i)).P1());
		Vector2D lp2 = Point(Line(the_polygon.LineNr(i)).P2());
		double dist = (lp1 + lp2).Y() *0.5;
		if( abs((lp2-lp1).X()) > MESH_STD_TOL )
			axis_distance.Add(dist);
	}

	QuicksortDouble(axis_distance,NaturalNumbers(the_polygon.NPoints()));

	if( abs(axis_distance.Last() - y_bar) <= MESH_STD_TOL )
//	if( abs(axis_distance.Last() - p1.Y())<=MESH_STD_TOL || abs(axis_distance.Last() - p2.Y())<=MESH_STD_TOL )
	{
		return 1;
	}
	return 0;
}

// determines if the linesegment is an "inner" edge of the polygon - line at y=0 (rotation axis) does not count as inner edge
int Crosssection::IsInnerEdge(int polygonnumber, int linesegmentnumber)
{
	ClosedPolygon& the_polygon = Polygon(polygonnumber);
	Vector2D p1 = Point(Line(the_polygon.LineNr(linesegmentnumber)).P1());
  Vector2D p2 = Point(Line(the_polygon.LineNr(linesegmentnumber)).P2());

	// * is not vertical
	if( abs((p2-p1).X()) <= MESH_STD_TOL ) return 0; // vertical line

	// check if the value for y_bar of the line is the highest (from all non vertical..)
	double y_bar = (p1+p2).Y() * 0.5;
	TArray<double> axis_distance;
	for(int i=1; i<=the_polygon.NLines(); i++)
	{
		Vector2D lp1 = Point(Line(the_polygon.LineNr(i)).P1());
		Vector2D lp2 = Point(Line(the_polygon.LineNr(i)).P2());
		double dist = (lp1 + lp2).Y() *0.5;
		if( abs((lp2-lp1).X()) > MESH_STD_TOL )
			axis_distance.Add(dist);
	}
	QuicksortDouble(axis_distance,NaturalNumbers(the_polygon.NPoints()));

	if( abs(axis_distance(1) - y_bar)  <=  MESH_STD_TOL )
//	if( abs(axis_distance(1) - p1.Y())<=MESH_STD_TOL || abs(axis_distance(1) - p2.Y())<=MESH_STD_TOL )
	{
		return 1;
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper functions

int Crosssection::CommonPointOfLineSegments(int linenr1, int linenr2)
{
	if( Line(linenr1).P1() == Line(linenr2).P1() ) return Line(linenr1).P1();
	if( Line(linenr1).P1() == Line(linenr2).P2() ) return Line(linenr1).P1();
	if( Line(linenr1).P2() == Line(linenr2).P1() ) return Line(linenr1).P2();
	if( Line(linenr1).P2() == Line(linenr2).P2() ) return Line(linenr1).P2();
  return 0;
}

int Crosssection::CommonPolygonOfLineSegments(int linenr1, int linenr2)
{
	for(int i=1; i<=NPolygons(); i++)
	{
		if( Polygon(i).LineNrs().Find(linenr1)>0 && Polygon(i).LineNrs().Find(linenr2)>0 )
			return i;
	}
  return 0;
}

int Crosssection::CommonPolygonOfPointAndLineSegment(int pointnr, int linenr)
{
	for(int i=1; i<=NPolygons(); i++)
	{
		if( Polygon(i).PointNrs().Find(pointnr)>0 && Polygon(i).LineNrs().Find(linenr)>0 )
			return i;
	}
  return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class SegmentedBlock: segmented hexahedral or cylinder (cutting planes)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodes, elements, cutting-planes, divisions (of segments), occupation (of segments)
// TODO: change access functions from int to THexDirections and TCylDirections enum 
// THINK ABOUT:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// + SegmentedHexBlock

// read a single hex block from edc !defined format!
void SegmentedHexBlock::ReadFromEDC(ElementDataContainer* edcblock)
{
	// disable popup warnings
	int edc_warn_lvl = edcblock->GetEDCWarningLevel();
	edcblock->SetEDCWarningLevel(0);

	// entries for the single block:   
	// geometry parameters: corner1, corner2
	edcblock->TreeGetVector3D("corner1", XMin(), YMin(), ZMin());
	edcblock->TreeGetVector3D("corner2", XMax(), YMax(), ZMax());
	// safety feature...
	if(XMin() > XMax()) 
	{ double swap = XMin(); XMin() = XMax(); XMax() = swap; } 
	if(YMin() > YMax()) 
	{ double swap = YMin(); YMin() = YMax(); YMax() = swap; } 
	if(ZMin() > ZMax()) 
	{ double swap = ZMin(); ZMin() = ZMax(); ZMax() = swap; } 
	// material, domain, color
	MatNr() = edcblock->TreeGetInt("material_number", MatNr());              
	Domain() = edcblock->TreeGetInt("domain_number", Domain());   
	edcblock->TreeGetVector3D("rgbcolor", Color().X(), Color().Y(), Color().Z());
	// mesh parameters:
	Vector3D dummy(XDivs_Block(), YDivs_Block(), ZDivs_Block());
	edcblock->TreeGetVector3D("divisions_xyz", dummy.X(), dummy.Y(), dummy.Z());
	BlockDivs().Set((int)dummy.X(), (int)dummy.Y(), (int)dummy.Z(),0);

	edcblock->SetEDCWarningLevel(edc_warn_lvl);
}

// creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void SegmentedHexBlock::GenerateIn(FEMesh* p_mesh, Vector3D meshsize)
{
	FEMesh_Generator& thegenerator = p_mesh->Generator();
	
	int xcl = XCut().Length(); // number of entries in cutting planes x
	int ycl = YCut().Length(); // number of entries in cutting planes y
	int zcl = ZCut().Length(); // number of entries in cutting planes z

	ResetBlockFaces();
	BlockElements().FlushArrays();
// generate subblocks
	for(int iz = 1; iz <= zcl-1; iz++)
	{
		for(int iy = 1; iy <= ycl-1; iy++)
		{
			for(int ix = 1; ix <=xcl-1; ix++)
			{
				if(Occupation(ix,iy,iz))
				{
					Vector3D pmin(XCut(ix), YCut(iy), ZCut(iz));
					Vector3D pmax(XCut(ix+1), YCut(iy+1), ZCut(iz+1));

					Box3D subbox(pmin, pmax);
					int3 divs(Division(ix,iy,iz));
					thegenerator.GenerateHexahedralBlock(subbox, divs, Domain(), MatNr(), Color());
				}
			}
		}
	}
	RefreshLists(p_mesh);
}

// recomputes the lists of elements, nodes, etc for the block
void SegmentedHexBlock::RefreshLists(FEMesh* p_mesh, int flag_faces) 
{
	BlockElements().FlushArrays();
  BlockNodes().FlushArrays();

	p_mesh->InitializeNodeSearchTree();
	p_mesh->GetElementsInBox(BlockElements(), GetBox());

	FEMesh_Set nodes_of_elements(TSetNodes);
	p_mesh->GetNodesOfElements(nodes_of_elements, BlockElements());
	
	p_mesh->GetNodesInBox(BlockNodes(), GetBox(), nodes_of_elements);

	if (flag_faces)
	{
		ResetBlockFaces();
		FEMesh_Set nodes(TSetNodes);
		FEMesh_Set faces(TSetFaces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(1.,0.,0.), XMin(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(1).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(1.,0.,0.), XMax(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(2).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,1.,0.), YMin(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(3).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,1.,0.), YMax(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(4).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,0.,1.), ZMin(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(5).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,0.,1.), ZMax(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(6).CopyFrom(faces);
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// + SegmentedCylBlock

// read a single cyl block from edc !defined format!
void SegmentedCylBlock::ReadFromEDC(ElementDataContainer* edcblock)
{
	// disable popup warnings
	int edc_warn_lvl = edcblock->GetEDCWarningLevel();
	edcblock->SetEDCWarningLevel(0);

	// entries for the single block:   
	// geometry parameters: inner radius, outer radius, bottom, top
	RI() = edcblock->TreeGetDouble("inner_radius", 0.);               // default: 0
	RO() = edcblock->TreeGetDouble("outer_radius", RI()+1.);          // default: thickness 1
	HB() = edcblock->TreeGetDouble("axis_pos_bottom", 0.);            // default: 0
	HT() = edcblock->TreeGetDouble("axis_pos_top", HB()+1.);          // default: height 1
	// material, domain
	MatNr() = edcblock->TreeGetInt("material_number", 1);             // default: material 1
	Domain() = edcblock->TreeGetInt("domain_number", 1);              // default: domain 1
	edcblock->TreeGetVector3D("rgbcolor", Color().X(), Color().Y(), Color().Z());
	// mesh parameters:
	RadDivs_Block() = edcblock->TreeGetDouble("radial_mesh_divisions", 1);  // default: 1 radial division
	TanDivs_Block() = edcblock->TreeGetDouble("angular_mesh_divisions", 1); // default: 1 angular segment for 45°
	AngRef() = edcblock->TreeGetDouble("outward_refinement", 0);      // default: no angular refinement in outward direction
	AxiDivs_Block() = edcblock->TreeGetDouble("axial_mesh_divisions", 1);   // default: 1 axial segment for entire height

	edcblock->SetLocalWarnTree(edc_warn_lvl);
}

// creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void SegmentedCylBlock::GenerateIn(FEMesh* p_mesh, Vector3D meshsize)
{
	FEMesh_Generator& thegenerator = p_mesh->Generator();
	int rcl = RadCut().Length(); // number of entries in cutting planes r
	int hcl = AxiCut().Length(); // number of entries in cutting planes h

	ResetBlockFaces();
	BlockElements().FlushArrays();
// generate subblocks
	for(int ih = 1; ih <= hcl-1; ih++)
	{
		for(int ir = 1; ir <=rcl-1; ir++)
		{
			if(Occupation(ir,ih))
			{
				double ri = RadCut(ir);
				double ro = RadCut(ir+1);
				double hb = AxiCut(ih);
				double ht = AxiCut(ih+1);
				int4 divs = Division(ir,ih);
				Vector3D p0(0.,0.,hb);
				Vector3D color(Color());
				
				if(ri == 0.0)
				{
					int angseg = TanDivs(ir,ih);
					int possible_inward_coarsening = 0;
					while ( angseg%2 == 0)
					{
						angseg /= 2;
						possible_inward_coarsening++;
					}
					int refinements = possible_inward_coarsening;
					
					int divisions_per_refinement;
					if(RadDivs(ir,ih)-2 >= (2*refinements))
						divisions_per_refinement = 2;           // "nice" double div
					else if (RadDivs(ir,ih)-2 >= refinements)
						divisions_per_refinement = 1;           // "flat" double div
					else
					{
						RadDivs(ir,ih)= 2+refinements;
						divisions_per_refinement = 1;           // requires additional radial segments (allowed since 
					}

					// compute sequence...
					int divisions_without_doublediv = RadDivs(ir,ih)-2 - refinements*divisions_per_refinement;
					double single_per_doublediv = (double)divisions_without_doublediv / (double)(1+refinements);

					double rad_inc = (ro-ri)/(double)RadDivs(ir,ih);
					thegenerator.GenerateCylinder(2*rad_inc, rad_inc, ht-hb, angseg, AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());

					double single = single_per_doublediv;               // number of single division rings to make (increases for each double div)
					int fact = 1;                                       // multiplication factor for starting angular segments 

					for(int i=3; i<=RadDivs(ir,ih); i++)
					{
						double r_i = ri +(i-1)*rad_inc;
						double r_o = ri +(i  )*rad_inc;
						double r_o_dd = ri +(i-1+divisions_per_refinement)*rad_inc;

						if(single < 1.)
						{
							thegenerator.GenerateHollowCylinderDoubleDiv(r_o_dd, r_i, ht-hb, angseg*fact, AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());
							if(divisions_per_refinement==2) i++;               // uses 2 radial divisions 
							fact *= 2;                                      // double the angular divisions from now
							single = single + single_per_doublediv;         // more single division rings
						}
						else       // IN CASE OF TOLERANCE PROBLEMS :    REMOVE ELSE; HAVE CONTINUE IN DOUBLEDIV BRANCH
						{
							thegenerator.GenerateHollowCylinder(r_o, r_i, ht-hb, angseg*fact, AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());
							single = single-1.;
						}
					}
				}
				else if (!Refine(ir,ih))   
				{
					double rad_inc = (ro-ri)/(double)RadDivs(ir,ih);
					for(int i=1; i<=RadDivs(ir,ih); i++)
					{
						double r_i = ri +(i-1)*rad_inc;
						double r_o = ri +(i  )*rad_inc;
						thegenerator.GenerateHollowCylinder(r_o, r_i, ht-hb, TanDivs(ir,ih), AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());
					}
				}
				else
				{
					int divisions_per_refinement;
					if(RadDivs(ir,ih) >= (2*Refine(ir,ih)))
						divisions_per_refinement = 2;           // "nice" double div
					else if (RadDivs(ir,ih) >= Refine(ir,ih))
						divisions_per_refinement = 1;           // "flat" double div
					else
					{
						RadDivs(ir,ih)= Refine(ir,ih);
						divisions_per_refinement = 1;           // requires additional radial segments (allowed since 
					}

					// compute sequence...
					int divisions_without_doublediv = RadDivs(ir,ih) - Refine(ir,ih)*divisions_per_refinement;
					double single_per_doublediv = (double)divisions_without_doublediv / (double)(1+Refine(ir,ih));
					
					double rad_inc = (ro-ri)/(double)RadDivs(ir,ih);
					double single = single_per_doublediv;               // number of single division rings to make (increases for each double div)
					int fact = 1;                                       // multiplication factor for starting angular segments 
					for(int i=1; i<=RadDivs(ir,ih); i++)
					{
						double r_i = ri +(i-1)*rad_inc;
						double r_o = ri +(i  )*rad_inc;
						double r_o_dd = ri +(i-1+divisions_per_refinement)*rad_inc;

						if(single < 1.)
						{
							thegenerator.GenerateHollowCylinderDoubleDiv(r_o_dd, r_i, ht-hb, TanDivs(ir,ih)*fact, AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());
							if(divisions_per_refinement) i++;               // uses 2 radial divisions 
							fact *= 2;                                      // double the angular divisions from now
							single = single + single_per_doublediv;         // more single division rings
						}
						else       // IN CASE OF TOLERANCE PROBLEMS :    REMOVE ELSE; HAVE CONTINUE IN DOUBLEDIV BRANCH
						{
							thegenerator.GenerateHollowCylinder(r_o, r_i, ht-hb, TanDivs(ir,ih)*fact, AxiDivs(ir,ih), p0, Domain(), MatNr(), color, IVector(), IVector());
							single = single-1.;
						}
					}
				}
				RefreshLists(p_mesh, Vector3D(0.,0.,0.), Vector3D(0.,0.,1)*(ht-hb));
			}
		}
	}
}

// recomputes the lists of elements, nodes, etc for the block
void SegmentedCylBlock::RefreshLists(FEMesh* p_mesh, Vector3D& base, Vector3D& axis, int flag_faces) 
{
	BlockElements().FlushArrays();
  BlockNodes().FlushArrays();

	p_mesh->GetNodesInCylinderShell(BlockNodes(), base, axis, RI(), RO(), HB(), HT());
	p_mesh->GetElementsInCylinderShell(BlockElements(), base, axis, RI(), RO(), HB(), HT());

	if (flag_faces)
	{
		ResetBlockFaces();
		FEMesh_Set nodes(TSetNodes);
		FEMesh_Set faces(TSetFaces);
		Vector3D a0 = axis; a0.Normalize();
		Vector3D m1 = base + HB()*a0;
		Vector3D m2 = base + HT()*a0;
		
		p_mesh->GetNodesOnCylinder(nodes, m1, m2, RI(), 0, BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(1).CopyFrom(faces);
		if(RO() > 1e-10)
		{
			p_mesh->GetNodesOnCylinder(nodes, m1, m2, RO(), 0, BlockNodes());
			p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
			BlockFace(2).CopyFrom(faces);
		}
		//p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,1.,0.), YMin(), Nodes());
		//p_mesh->FacesFromNodes(faces, nodes, Elements());
		//BlockFace(3).CopyFrom(faces);
		//p_mesh->GetNodesOnPlane(nodes, Vector3D(0.,1.,0.), YMax(), Nodes());
		//p_mesh->FacesFromNodes(faces, nodes, Elements());
		//BlockFace(4).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, axis, HB(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(5).CopyFrom(faces);
		p_mesh->GetNodesOnPlane(nodes, axis, HT(), BlockNodes());
		p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		BlockFace(6).CopyFrom(faces);
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// + SegmentedQuadBlock

// read a single quad block from edc !defined format!
void SegmentedQuadBlock::ReadFromEDC(ElementDataContainer* edcblock)
{
	// disable popup warnings
	int edc_warn_lvl = edcblock->GetEDCWarningLevel();
	edcblock->SetEDCWarningLevel(0);

	// entries for the single block:   
	// geometry parameters: corner1, corner2
	edcblock->TreeGetVector2D("corner1", XMin(), YMin());
	edcblock->TreeGetVector2D("corner2", XMax(), YMax());
	// safety feature...
	if(XMin() > XMax()) 
	{ double swap = XMin(); XMin() = XMax(); XMax() = swap; } 
	if(YMin() > YMax()) 
	{ double swap = YMin(); YMin() = YMax(); YMax() = swap; } 
	ZMin() = -1. * MESH_STD_TOL;
	ZMax() =  1. * MESH_STD_TOL;
	Thickness() = edcblock->TreeGetDouble("thickness", Thickness());
	// material, domain, color
	MatNr() = edcblock->TreeGetInt("material_number", MatNr());              
	Domain() = edcblock->TreeGetInt("domain_number", Domain());   
	edcblock->TreeGetVector3D("rgbcolor", Color().X(), Color().Y(), Color().Z());
	// mesh parameters:
	Vector2D dummy(XDivs_Block(), YDivs_Block());
	edcblock->TreeGetVector2D("divisions_xy", dummy.X(), dummy.Y());
	BlockDivs().Set((int)dummy.X(), (int)dummy.Y(), 1);

	edcblock->SetEDCWarningLevel(edc_warn_lvl);
}

// creates the SegmentedBlock in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void SegmentedQuadBlock::GenerateIn(FEMesh* p_mesh, Vector3D meshsize)
{
	FEMesh_Generator& thegenerator = p_mesh->Generator();
	
	int xcl = XCut().Length(); // number of entries in cutting planes x
	int ycl = YCut().Length(); // number of entries in cutting planes y

	ResetBlockFaces();
	BlockElements().FlushArrays();
// generate subblocks
	for(int iy = 1; iy <= ycl-1; iy++)
	{
		for(int ix = 1; ix <=xcl-1; ix++)
		{
			if(Occupation(ix,iy))
			{
				Vector2D pmin(XCut(ix), YCut(iy));
				Vector2D pmax(XCut(ix+1), YCut(iy+1));

				Box2D subbox(pmin, pmax);
				int2 divs(Division(ix,iy));
				thegenerator.GenerateQuadrilateralBlock(subbox, Thickness(), divs, Domain(), MatNr(), Color());
			}
		}
	}
	RefreshLists(p_mesh);
}

// recomputes the lists of elements, nodes, etc for the block
void SegmentedQuadBlock::RefreshLists(FEMesh* p_mesh, int flag_faces) 
{
	BlockElements().FlushArrays();
  BlockNodes().FlushArrays();

	p_mesh->InitializeNodeSearchTree();
	p_mesh->GetElementsInBox(BlockElements(), GetBox());

	FEMesh_Set nodes_of_elements(TSetNodes);
	p_mesh->GetNodesOfElements(nodes_of_elements, BlockElements());
	
	p_mesh->GetNodesInBox(BlockNodes(), GetBox(), nodes_of_elements);

	if (flag_faces)
	{
		ResetBlockFaces();
		FEMesh_Set nodes(TSetNodes);
		FEMesh_Set faces(TSetFaces);
// no faces from nodes 2d yet...
		//p_mesh->GetNodesOnPlane(nodes, Vector2D(1.,0.), XMin(), BlockNodes());
		//p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		//BlockFace(1).CopyFrom(faces);
		//p_mesh->GetNodesOnPlane(nodes, Vector2D(1.,0.), XMax(), BlockNodes());
		//p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		//BlockFace(2).CopyFrom(faces);
		//p_mesh->GetNodesOnPlane(nodes, Vector2D(0.,1.), YMin(), BlockNodes());
		//p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		//BlockFace(3).CopyFrom(faces);
		//p_mesh->GetNodesOnPlane(nodes, Vector2D(0.,1.), YMax(), BlockNodes());
		//p_mesh->FacesFromNodes(faces.Faces(), nodes.Nodes(), BlockElements().Elements());
		//BlockFace(4).CopyFrom(faces);
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  (base) class MeshedPart: several hexahedral or cylinder (SegmentedBlocks)with consistent mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contains: nodes, elements, blocks, cutting-planes, divisions (of segments), occupation (of segments),
// TODO: change access functions from int to THexDirections and TCylDirections enum 
// THINK ABOUT:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  MeshedHexPart

// read a sub-edc (part) containing several hexblocks !defined format!
void MeshedHexPart::ReadFromEDC(ElementDataContainer* edcpart)
{
	// disable popup warnings
	int edc_warn_lvl = edcpart->GetEDCWarningLevel();
	edcpart->SetEDCWarningLevel(0);

	// entries for the entire part:   
	// get all sub EDCs named "HexBlock###"
	int next_block_nr = 1;
	while(next_block_nr)
	{
		int idx = edcpart->Find(mystr("HexBlock")+mystr(next_block_nr));
		if(idx && (edcpart->Get(idx)).IsEDC())
		{
			ElementDataContainer* edcblock = (edcpart->Get(idx)).GetEDC();
			SegmentedHexBlock* newblock = new SegmentedHexBlock;
			newblock->ReadFromEDC(edcblock);
			Blocks().Add(newblock);
			next_block_nr++; // next iteration;
		}
		else
			next_block_nr = 0; // exit loop
	}
	edcpart->SetEDCWarningLevel(edc_warn_lvl);
}

// cuts the blocks in the list - prepares consistend meshing of the part
void MeshedHexPart::Cut()
{
	FindCuttingPlanes();
	ComputeOccupationArray();
	ComputeDivisionsArray();
	SetBlockData();
}

// computes all cutting planes (sorted), arrays include overall limits	
void MeshedHexPart::FindCuttingPlanes()
{
//size of part (surrounding box)
	Box3D overall = GetSurroundingBox();

// cutting planes due to the blocks
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedHexBlock& theblock = (SegmentedHexBlock&) Block(i); 
		XCut().Add(theblock.XMin());
		XCut().Add(theblock.XMax());
		YCut().Add(theblock.YMin());
		YCut().Add(theblock.YMax());
		ZCut().Add(theblock.ZMin());
		ZCut().Add(theblock.ZMax());
		for(int i=1; i<= theblock.XCut().Length(); i++)
			XCut().Add(theblock.XCut(i));
		for(int i=1; i<= theblock.YCut().Length(); i++)
			YCut().Add(theblock.YCut(i));
		for(int i=1; i<= theblock.ZCut().Length(); i++)
			ZCut().Add(theblock.ZCut(i));
	}

	RemoveRedundantEntries_ErrorTolerance(XCut(), 1e-10, 1);
	RemoveRedundantEntries_ErrorTolerance(YCut(), 1e-10, 1);
	RemoveRedundantEntries_ErrorTolerance(ZCut(), 1e-10, 1);

// remove all values out of surrounding box, only cutting planes 
	while (XCut(1) < overall.PMin().X()) XCut().Erase(1);
	while (XCut().Last() > overall.PMax().X()) XCut().Erase(XCut().Length());
	while (YCut(1) < overall.PMin().Y()) YCut().Erase(1);
	while (YCut().Last() > overall.PMax().Y()) YCut().Erase(YCut().Length());
	while (ZCut(1) < overall.PMin().Z()) ZCut().Erase(1);
	while (ZCut().Last() > overall.PMax().Z()) ZCut().Erase(ZCut().Length());
}

// computes which block occupies the volume-segment (the latter block always gets priority)
void MeshedHexPart::ComputeOccupationArray()
{
	int xcl = XCut().Length();
	int ycl = YCut().Length();
	int zcl = ZCut().Length();
//	Box3D overall = GetSurroundingBox();
	OccuArray().Flush();
	OccuArray().SetLen((xcl-1)*(ycl-1)*(zcl-1));
	OccuArray().SetAll(0);

	for(int i=1; i <= Blocks().Length(); i++) // all blocks
	{
		for(int iz = 1; iz <= zcl-1; iz++)
		{
			for(int iy = 1; iy <= ycl-1; iy++)
			{
				for(int ix = 1; ix <= xcl-1; ix++)
				{
					Vector3D center;
					center.X() = (XCut(ix) + XCut(ix+1)) * 0.5;
					center.Y() = (YCut(iy) + YCut(iy+1)) * 0.5;
					center.Z() = (ZCut(iz) + ZCut(iz+1)) * 0.5;

					if( GetBox(i).IsIn(center) ) // occupy volume section with block
					{
						int idx = 1 + (ix-1) + (iy-1)*(xcl-1) + (iz-1)*(xcl-1)*(ycl-1);
						Occupation( idx ) = i;
					}
				}
			}
		}
	}
}

// computes which block occupies the volume-segment (the latter block always gets priority)
void MeshedHexPart::ComputeDivisionsArray()
{
	int xcl = XCut().Length();
	int ycl = YCut().Length();
	int zcl = ZCut().Length();

	TArray<int> divs_x; divs_x.SetLen(xcl-1); divs_x.SetAll(0);
	TArray<int> divs_y; divs_y.SetLen(ycl-1); divs_y.SetAll(0);
	TArray<int> divs_z; divs_z.SetLen(zcl-1); divs_z.SetAll(0);

	for(int i=1; i <= Blocks().Length(); i++) // all blocks
	{
		SegmentedHexBlock& theblock = Block(i); 
		double xtotal = theblock.XMax()-theblock.XMin();
		double xfact = (double)theblock.XDivs_Block() / xtotal;
		double ytotal = theblock.YMax()-theblock.YMin();
		double yfact = (double)theblock.YDivs_Block() / ytotal;
		double ztotal = theblock.ZMax()-theblock.ZMin();
		double zfact = (double)theblock.ZDivs_Block() / ztotal;

		// compute indices range of this block in set
		int xfirst = XCut().Find(theblock.XMin());
		int xlast = XCut().Find(theblock.XMax());
		int yfirst = YCut().Find(theblock.YMin());
		int ylast = YCut().Find(theblock.YMax());
		int zfirst = ZCut().Find(theblock.ZMin());
		int zlast = ZCut().Find(theblock.ZMax());

		// USE SOME OTHER FUNCTION HERE FOR SMOOTHER BUT MORE ELEMENT INTENSIVE DIVISION
		// compute number of divisions x,y&z for each x,y&z segment of this block
		for(int j=xfirst; j<= xlast-1; j++)
		{
			double xthis = XCut(j+1) - XCut(j);
			int divs_divs = (int) ceil( xthis * xfact );      
			if (divs_x(j) < divs_divs) 
				divs_x(j) = divs_divs;
			if(MeshX() > 0.)		
			{
				int divs_mesh = (int) ceil( xthis / MeshX() );
				if (divs_x(j) < divs_mesh) 
					divs_x(j) = divs_mesh;
			}
		}
		for(int j=yfirst; j<= ylast-1; j++)
		{
			double ythis = YCut(j+1) - YCut(j);
			int divs_divs = (int) ceil( ythis * yfact );      
			if (divs_y(j) < divs_divs) 
				divs_y(j) = divs_divs;
			if(MeshY() > 0.)		
			{
				//EK 20120419 - bug fix 
				//int divs_mesh = (int) ceil( ythis / MeshX() );
				int divs_mesh = (int) ceil( ythis / MeshY() );
				if (divs_y(j) < divs_mesh) 
					divs_y(j) = divs_mesh;
			}
		}
		for(int j=zfirst; j<= zlast-1; j++)
		{
			double zthis = ZCut(j+1) - ZCut(j);
			int divs_divs = (int) ceil( zthis * zfact );      
			if (divs_z(j) < divs_divs) 
				divs_z(j) = divs_divs;
			if(MeshZ() > 0.)		
			{
				int divs_mesh = (int) ceil( zthis / MeshZ() );
				if (divs_z(j) < divs_mesh) 
					divs_z(j) = divs_mesh;
			}
		}
	}

// fill the DivsArray
	DivsArray().Flush();
	DivsArray().SetLen((xcl-1)*(ycl-1)*(zcl-1));
	DivsArray().SetAll(int4(1,1,1,1));
	for(int ix=1; ix<=divs_x.Length(); ix++)
	{
		for(int iy=1; iy<=divs_y.Length(); iy++)
		{
			for(int iz=1; iz<=divs_z.Length(); iz++)
			{
				int idx = 1 + (ix-1) + (iy-1)*(xcl-1) + (iz-1)*(ycl-1)*(xcl-1);
				int4 divs(divs_x(ix), divs_y(iy), divs_z(iz), 0);
				Division( idx ) = divs;
			}
		}
	}
}

// writes computed sectioning to the blocks
void MeshedHexPart::SetBlockData()
{
	int xcl = XCut().Length();
	int ycl = YCut().Length();
	int zcl = ZCut().Length();
//	Box3D overall = GetSurroundingBox();
	
	for(int i=1; i <= Blocks().Length(); i++)
	{
		SegmentedHexBlock& theblock = (SegmentedHexBlock&) Block(i); 

// compute indices range of this block in set
		int xfirst = XCut().Find(theblock.XMin());
		int xlast = XCut().Find(theblock.XMax());
		int yfirst = YCut().Find(theblock.YMin());
		int ylast = YCut().Find(theblock.YMax());
		int zfirst = ZCut().Find(theblock.ZMin());
		int zlast = ZCut().Find(theblock.ZMax());

		theblock.XCut().CopyFrom(XCut(),xfirst,xlast);
		theblock.YCut().CopyFrom(YCut(),yfirst,ylast);
		theblock.ZCut().CopyFrom(ZCut(),zfirst,zlast);
		theblock.SetAllOccupations(0);
// loop over all sections of the current block
		for(int iz = zfirst; iz <= zlast-1; iz++)
		{
			for(int iy = yfirst; iy <= ylast-1; iy++)
			{
				for(int ix = xfirst; ix <= xlast-1; ix++)
				{
					int idx =	1 + (ix-1) + (iy-1)*(xcl-1) + (iz-1)*(xcl-1)*(ycl-1)	;
					int idx_loc = theblock.IofXYZ(ix-xfirst+1, iy-yfirst+1, iz-zfirst+1); // dgb
					int occu = Occupation( idx ); 
					if(occu == i)
					{
						theblock.Occupation(ix-xfirst+1, iy-yfirst+1, iz-zfirst+1) = occu; // set occupation
					}
					int4 divs = Division( idx );
					theblock.Division(ix-xfirst+1, iy-yfirst+1, iz-zfirst+1) = divs;     // set divisions
				}
			}
		}
	}
}

// creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void MeshedHexPart::Generate()
{
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedHexBlock& theblock = (SegmentedHexBlock&)Block(i);
		theblock.GenerateIn(MeshPtr(), MeshWidth());
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();
	}
}

// recomputes the lists of elements, nodes, etc for the entire part
void MeshedHexPart::RefreshLists(int flag_faces)
{
	Box3D thebox = GetSurroundingBox();
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	if (flag_faces) ResetPartFaces(6);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedHexBlock& theblock = (SegmentedHexBlock&)Block(i);
		theblock.RefreshLists(MeshPtr(), flag_faces);
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();

		if (flag_faces)
		{
			if( thebox.PMin().X() == theblock.XMin() ) PartFace(1) += theblock.BlockFace(1);
			if( thebox.PMax().X() == theblock.XMax() ) PartFace(2) += theblock.BlockFace(2);
			if( thebox.PMin().Y() == theblock.YMin() ) PartFace(3) += theblock.BlockFace(3);
			if( thebox.PMax().Y() == theblock.YMax() ) PartFace(4) += theblock.BlockFace(4);
			if( thebox.PMin().Z() == theblock.ZMin() ) PartFace(5) += theblock.BlockFace(5);
			if( thebox.PMax().Z() == theblock.ZMax() ) PartFace(6) += theblock.BlockFace(6);
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  MeshedCylPart

// read a sub-edc (part) containing several cylblocks !defined format!
void MeshedCylPart::ReadFromEDC(ElementDataContainer* edcpart)
{
	// disable popup warnings
	int edc_warn_lvl = edcpart->GetEDCWarningLevel();
	edcpart->SetEDCWarningLevel(0);
	// entries for the entire part:   
	// axis-vector 
	Axis() = Vector3D(0., 0., 1.);    // default: z-axis
	edcpart->TreeGetVector3D("axis", Axis().X(), Axis().Y(), Axis().Z());
	if(Axis()==Vector3D(0.,0.,1.))
		Rotation() = Matrix3D( 1,0,0, 0,1,0, 0,0,1 );
	else
	{
	// GET THE ROTATION MATRIX !!!
	// R_n(a) = (
		Vector3D x = Vector3D(0.,0.,1.);									// z-axis (from generate cylinder)
		Vector3D n0 = x.Cross(Axis());										// rotation axis
		n0.Normalize(); 
		double a = acos((x*Axis()) / Axis().Norm());      // angle between z-axis and chosen axis
	  Rotation()(1,1) = cos(a) + n0.X()*n0.X()*(1-cos(a));
		Rotation()(1,2) = n0.X()*n0.Y()*(1-cos(a)) - n0.Z()*sin(a);
		Rotation()(1,3) = n0.X()*n0.Z()*(1-cos(a)) + n0.Y()*sin(a);
		Rotation()(2,1) = n0.Y()*n0.X()*(1-cos(a)) + n0.Z()*sin(a);
    Rotation()(2,2) = cos(a) + n0.Y()*n0.Y()*(1-cos(a));
		Rotation()(2,3) = n0.Y()*n0.Z()*(1-cos(a)) - n0.X()*sin(a);
		Rotation()(3,1) = n0.Z()*n0.X()*(1-cos(a)) - n0.Y()*sin(a);
		Rotation()(3,2) = n0.Z()*n0.Y()*(1-cos(a)) + n0.X()*sin(a);
    Rotation()(3,3) = cos(a) + n0.Z()*n0.Z()*(1-cos(a));
	}
	
	// basepoint
	edcpart->TreeGetVector3D("base", Base().X(), Base().Y(), Base().Z());
	Translation() = Base();

	// quadrants
	Quadrants() = edcpart->TreeGetInt("quadrants", 1); // default: create one quadrant
	if(!(Quadrants() & 7))
		MeshPtr()->GetMBS()->UO(UO_LVL_warn) << "FEMesh::ReadCylBlocksFromEDC: value QUADRANTS must be set to 1, 2 or 4 !!\n";

	// get all sub EDCs named "Cylinder###"
	int next_block_nr = 1;
	while(next_block_nr)
	{
		int idx = edcpart->Find(mystr("Cylinder")+mystr(next_block_nr));
		if(idx && (edcpart->Get(idx)).IsEDC())
		{
			ElementDataContainer* edcblock = (edcpart->Get(idx)).GetEDC();
			SegmentedCylBlock* newblock = new SegmentedCylBlock;
			newblock->ReadFromEDC(edcblock);
			Blocks().Add(newblock);
			next_block_nr++; // next iteration;
		}
		else
			next_block_nr = 0; // exit loop
	}
	edcpart->SetLocalWarnTree(edc_warn_lvl);
}

// cuts the blocks in the list - prepares consistend meshing of the part
void MeshedCylPart::Cut()
{
	FindCuttingPlanes();
	ComputeOccupationArray();
	ComputeDivisionsArray();
	SetBlockData();
}

// computes all cutting planes (sorted), arrays include overall limits	
void MeshedCylPart::FindCuttingPlanes()
{
	double tol = 1e-10;
//size of part (surrounding box)
	Box3D overall = GetSurroundingBox();  

// cutting planes due to the blocks
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&) Block(i); 
		RadCut().Add(theblock.RI());
		RadCut().Add(theblock.RO());
		AxiCut().Add(theblock.HB());
		AxiCut().Add(theblock.HT());
		for(int i=1; i<= theblock.RadCut().Length(); i++)
			RadCut().Add(theblock.RadCut(i));
		for(int i=1; i<= theblock.AxiCut().Length(); i++)
			AxiCut().Add(theblock.AxiCut(i));
	}
	RemoveRedundantEntries_ErrorTolerance(RadCut(), tol, 1);
	RemoveRedundantEntries_ErrorTolerance(AxiCut(), tol, 1);

// remove all values out of surrounding box, only cutting planes 
	while (RadCut(1) < overall.PMin().X()) RadCut().Erase(1);
	while (RadCut().Last() > overall.PMax().X()) RadCut().Erase(RadCut().Length());
	while (AxiCut(1) < overall.PMin().Z()) AxiCut().Erase(1);
	while (AxiCut().Last() > overall.PMax().Z()) AxiCut().Erase(AxiCut().Length());
}

// computes which block occupies the volume-segment (the latter block always gets priority)
void MeshedCylPart::ComputeOccupationArray()
{
	int rcl = RadCut().Length();
	int hcl = AxiCut().Length();
	OccuArray().Flush();
	OccuArray().SetLen((rcl-1)*(hcl-1));
	OccuArray().SetAll(0);

	for(int i=1; i <= Blocks().Length(); i++) // all blocks
	{
		for(int ih = 1; ih <= hcl-1; ih++)
		{
			for(int ir = 1; ir <= rcl-1; ir++)
			{
				Vector3D center;
				center.X() = (RadCut(ir) + RadCut(ir+1)) * 0.5;
				center.Y() = 0.5; // as default boxsize is 0..1
				center.Z() = (AxiCut(ih) + AxiCut(ih+1)) * 0.5;

				if( GetBox(i).IsIn(center) ) // occupy volume section with block
				{
					int idx = 1 + (ir-1) + (ih-1)*(rcl-1);
					Occupation( idx ) = i;
				}
			}
		}
	}
}

// computes the divisions for the blocks
void MeshedCylPart::ComputeDivisionsArray()
{
	int rcl = RadCut().Length();
	int hcl = AxiCut().Length();

	TArray<int> divs_rad; divs_rad.SetLen(rcl-1); divs_rad.SetAll(0);
	TArray<int> divs_tan; divs_tan.SetLen(rcl-1); divs_tan.SetAll(0);
	TArray<int> refine;   refine.SetLen(rcl-1);   divs_tan.SetAll(0);
	TArray<int> divs_axi; divs_axi.SetLen(hcl-1); divs_axi.SetAll(0);

	for(int i=1; i <= Blocks().Length(); i++) // all blocks
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&) Block(i); 
		double rtotal = theblock.RO()-theblock.RI();
		double rfact = (double)theblock.RadDivs_Block() / rtotal;
    double htotal = theblock.HT()-theblock.HB();
		double hfact = (double)theblock.AxiDivs_Block() / htotal;

// compute indices range of this block in set
		int rfirst = RadCut().Find(theblock.RI());
		int rlast = RadCut().Find(theblock.RO());
		int hfirst = AxiCut().Find(theblock.HB());
		int hlast = AxiCut().Find(theblock.HT());

// USE SOME OTHER FUNCTION HERE FOR SMOOTHER BUT MORE ELEMENT INTENSIVE DIVISION
// compute number of radial divisions (and tangential divisions) for each radial segment of this block
		for(int j=rfirst; j<= rlast-1; j++)
		{
			double rthis = RadCut(j+1) - RadCut(j);
			int divs_divs = (int) ceil( rthis * rfact );      
			if (divs_rad(j) < divs_divs) 
				divs_rad(j) = divs_divs;
			if(MeshRad() > 0.)		
			{
			int divs_mesh = (int) ceil( rthis / MeshRad() );
			if (divs_rad(j) < divs_mesh) 
				divs_rad(j) = divs_mesh;
			}
		}
// compute number of axial divisions for each axial segment of this block
		for(int j=hfirst; j<= hlast-1; j++)
		{
			double hthis = AxiCut(j+1) - AxiCut(j);
			int divs_divs = (int) ceil( hthis * hfact );      
			if (divs_axi(j) < divs_divs) 
				divs_axi(j) = divs_divs;
			if(MeshAxi() >0.) // meshsize 0:  no upper limit
			{
				int divs_mesh = (int) ceil( hthis / MeshAxi() );
				if (divs_axi(j) < divs_mesh) 
					divs_axi(j) = divs_mesh;
			}
		}
	}

// compute number of tangential divisions for each axial segment of this block
// STEP 1: BLOCKWISE fill the highest number of segments for each RadCut, fill in Array AngSeg()
  AngSeg().SetLen(RadCut().Length());
	AngSeg().SetAll(1);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&) Block(i); 
		int n_ri = RadCut().Find(theblock.RI());
		int n_ro = RadCut().Find(theblock.RO());
		
		int divs_divs_ri = theblock.TanDivs_Block();
		int divs_divs_ro = theblock.TanDivs_Block()*(int)(pow(2.,theblock.AngRef()));
		if(AngSeg(n_ri) < divs_divs_ri)
			AngSeg(n_ri) = divs_divs_ri;
		if(AngSeg(n_ro) < divs_divs_ro)
			AngSeg(n_ro) = divs_divs_ro;

		if(MeshTan() > 0.) // meshsize 0:  no upper limit
		{
			int divs_mesh_ri = (int) ceil( (RadCut(n_ri)*0.25*MY_PI) / MeshTan() );
			int divs_mesh_ro = (int) ceil( (RadCut(n_ro)*0.25*MY_PI) / MeshTan() );
			if(AngSeg(n_ri) < divs_mesh_ri)
				AngSeg(n_ri) = divs_mesh_ri;
			if(AngSeg(n_ro) < divs_mesh_ro)
				AngSeg(n_ro) = divs_mesh_ro;
		}
	}
// find the most segments for smallest nonzero radius (45°)	
	int segments_at_smallest_radius = AngSeg(1);
	if(RadCut(1) < 1e-10)
	{
		AngSeg(1) = AngSeg(2);
		segments_at_smallest_radius = AngSeg(2);
	}
	int segments = segments_at_smallest_radius;
// work outward and calculate segments for each RCut
	for(int i=2; i<= RadCut().Length(); i++)
	{
		if( AngSeg(i) > segments)                        // <-- use WHILE here to allow multiple refinements (LATER)
			segments *=2;
		AngSeg(i) = segments;
	}

// STEP 2: compute number of divisions and refineflag for each segment
// LATER: allow multiple angular refinements
	for(int i=1; i<= rcl-1; i++)
	{
// previous ring
		int prev_divs = AngSeg(1);
	  int prevflag = 0;
		if(i>1)
		{
			prev_divs = divs_tan(i-1);
			prevflag = refine(i-1);
		}
// this ring
		int d = AngSeg(i+1);
		if( d <= prev_divs) 
		{
			if(!prevflag) 
			{ divs_tan(i) = prev_divs; refine(i) = 0; }
			else 
			{	divs_tan(i) = 2*prev_divs; refine(i) = 0; }
		}
		else if ( d <= (2*prev_divs))
		{
			if(!prevflag)
			{	divs_tan(i) = prev_divs; refine(i) = 1; }
			else
			{	divs_tan(i) = 2*prev_divs; refine(i) = 0; }
		}
		else if ( d <= (4*prev_divs))
		{ 
			if(!prevflag)
			{	
				divs_tan(i) = prev_divs; refine(i) = 1; 
			  /*UO(UO_LVL_warn) << "angular mesh reinement limit met! - intruduce additional radial cuts !\n";*/ 
			}
			else
			{ divs_tan(i) = 2*prev_divs; refine(i) = 1; }
		}
		else
		{
			divs_tan(i) = 2*prev_divs; refine(i) = 1; 
			/*UO(UO_LVL_warn) << "angular mesh reinement limit met! - intruduce additional radial cuts !\n";*/
		}
	}

// fill the DivsArray
	DivsArray().Flush();
	DivsArray().SetLen((rcl-1)*(hcl-1));
	DivsArray().SetAll(int4(1,1,1,1));
	for(int ir=1; ir<=divs_rad.Length(); ir++)
	{
		for(int ih=1; ih<=divs_axi.Length(); ih++)
		{
			int idx = 1 + (ir-1) + (ih-1)*(rcl-1);
			int4 divs(divs_rad(ir), divs_tan(ir), divs_axi(ih), refine(ir));
      Division( idx ) = divs;
		}
	}
}

// writes computed sectioning to the blocks
void MeshedCylPart::SetBlockData()
{
	int rcl = RadCut().Length();
	int hcl = AxiCut().Length();
	
	for(int i=1; i <= Blocks().Length(); i++)
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&) Block(i); 

// compute indices range of this block in set
		int rfirst = RadCut().Find(theblock.RI());
		int rlast = RadCut().Find(theblock.RO());
		int hfirst = AxiCut().Find(theblock.HB());
		int hlast = AxiCut().Find(theblock.HT());

		theblock.RadCut().CopyFrom(RadCut(),rfirst,rlast);
		theblock.AxiCut().CopyFrom(AxiCut(),hfirst,hlast);
		theblock.SetAllOccupations(0);
// loop over all sections of the current block
		for(int ih = hfirst; ih <= hlast; ih++)
		{
			for(int ir = rfirst; ir <= rlast; ir++)
			{
				int idx =	1 + (ir-1) + (ih-1)*(rcl-1);
				int occu = Occupation( idx ); 
				if(occu == i)
				{
					theblock.Occupation(ir-rfirst+1, ih-hfirst+1) = occu; // set occupation
				}
				int4 divs = Division( idx );
				theblock.Division(ir-rfirst+1, ih-hfirst+1) = divs;         // set divisions
			}
		}
	}
}

// creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void MeshedCylPart::Generate()
{
// create 1st quadrant
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&)Block(i);
		theblock.GenerateIn(MeshPtr(), MeshWidth());
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();
	}
// quadrants - mirror operations
	if (Quadrants() == 2)
	{
		MeshPtr()->MirrorMeshAtPlane(Vector3D(1.,0.,0.), 0);
	}
	if(Quadrants() == 4)
	{
		MeshPtr()->MirrorMeshAtPlane(Vector3D(1.,0.,0.), 0);
		MeshPtr()->MirrorMeshAtPlane(Vector3D(0.,1.,0.), 0);
	}
// rotation and translation
	FEMesh3D& themesh3d = (FEMesh3D&) *MeshPtr();
	themesh3d.Transform(Translation(), Rotation()/*, Nodes()*/);
}

// recomputes the lists of elements, nodes, etc for the entire part
void MeshedCylPart::RefreshLists(int flag_faces)
{
	Box3D thebox = GetSurroundingBox();
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	if (flag_faces) ResetPartFaces(6);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedCylBlock& theblock = (SegmentedCylBlock&)Block(i);
		theblock.RefreshLists(MeshPtr(), Base(), Axis(), flag_faces);
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();

		if (flag_faces)
		{
			if( thebox.PMin().X() == theblock.RI() ) PartFace(1) += theblock.BlockFace(1);
			if( thebox.PMax().X() == theblock.RO() ) PartFace(2) += theblock.BlockFace(2);
			
			if( thebox.PMin().Z() == theblock.HB() ) PartFace(5) += theblock.BlockFace(5);
			if( thebox.PMax().Z() == theblock.HT() ) PartFace(6) += theblock.BlockFace(6);
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  MeshedQuadPart

// read a sub-edc (part) containing several hexblocks !defined format!
void MeshedQuadPart::ReadFromEDC(ElementDataContainer* edcpart)
{
	// disable popup warnings
	int edc_warn_lvl = edcpart->GetEDCWarningLevel();
	edcpart->SetEDCWarningLevel(0);

	// entries for the entire part:   
	// get all sub EDCs named "HexBlock###"
	int next_block_nr = 1;
	while(next_block_nr)
	{
		int idx = edcpart->Find(mystr("QuadBlock")+mystr(next_block_nr));
		if(idx && (edcpart->Get(idx)).IsEDC())
		{
			ElementDataContainer* edcblock = (edcpart->Get(idx)).GetEDC();
			SegmentedQuadBlock* newblock = new SegmentedQuadBlock;
			newblock->ReadFromEDC(edcblock);
			Blocks().Add(newblock);
			next_block_nr++; // next iteration;
		}
		else
			next_block_nr = 0; // exit loop
	}
	edcpart->SetEDCWarningLevel(edc_warn_lvl);
}

// cuts the blocks in the list - prepares consistend meshing of the part
void MeshedQuadPart::Cut()
{
	FindCuttingPlanes();
	ComputeOccupationArray();
	ComputeDivisionsArray();
	SetBlockData();
}

// computes all cutting planes (sorted), arrays include overall limits	
void MeshedQuadPart::FindCuttingPlanes()
{
	MeshedHexPart::FindCuttingPlanes();
}

// computes which block occupies the volume-segment (the latter block always gets priority)
void MeshedQuadPart::ComputeOccupationArray()
{
	MeshedHexPart::ComputeOccupationArray();
}

// computes which block occupies the volume-segment (the latter block always gets priority)
void MeshedQuadPart::ComputeDivisionsArray()
{
	MeshedHexPart::ComputeDivisionsArray();
}

// writes computed sectioning to the blocks
void MeshedQuadPart::SetBlockData()
{
	MeshedHexPart::SetBlockData();
}

// creates the Part in the defined Mesh - entries in Arrays Nodes, Elements, etc are remembered for LATEST call only
void MeshedQuadPart::Generate()
{
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedQuadBlock& theblock = (SegmentedQuadBlock&)Block(i);
		theblock.GenerateIn(MeshPtr(), MeshWidth());
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();
	}
}

// recomputes the lists of elements, nodes, etc for the entire part
void MeshedQuadPart::RefreshLists(int flag_faces)
{
	Box3D thebox = GetSurroundingBox();
	PartElements().SetFEMesh_Set(TSetElements);
	PartNodes().SetFEMesh_Set(TSetNodes);
	if (flag_faces) ResetPartFaces(6);
	for(int i=1; i<=Blocks().Length(); i++)
	{
		SegmentedQuadBlock& theblock = (SegmentedQuadBlock&)Block(i);
		theblock.RefreshLists(MeshPtr(), flag_faces);
		PartElements() += theblock.BlockElements();
		PartNodes() += theblock.BlockNodes();

// no faces from nodes 2d yet...
		//if (flag_faces)
		//{
		//	if( thebox.PMin().X() == theblock.XMin() ) PartFace(1) += theblock.BlockFace(1);
		//	if( thebox.PMax().X() == theblock.XMax() ) PartFace(2) += theblock.BlockFace(2);
		//	if( thebox.PMin().Y() == theblock.YMin() ) PartFace(3) += theblock.BlockFace(3);
		//	if( thebox.PMax().Y() == theblock.YMax() ) PartFace(4) += theblock.BlockFace(4);
		//}
	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ THE SCRAP YARD - these functions sould be eliminated a.s.a.p.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
