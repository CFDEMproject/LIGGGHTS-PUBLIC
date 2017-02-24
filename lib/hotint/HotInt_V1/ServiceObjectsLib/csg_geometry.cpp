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

#include "csg_geometry.h"
#include "windows.h" //for sleep

extern UserOutputInterface * global_uo;


// writes and executes bat-file to mesh geometry, defined in geo-file, in NETGEN
bool CreateMeshWithNETGEN(const char * geoFilename,const char * meshFilename,const char * NETGENPath,const char * WorkDir,const char * NETGENMeshFileType, const char* coarsity)
{
	// describtion:
	// 1) generates bat-file, that starts meshing in NETGEN
	// 2) executes this bat-file
	// 3) waits until meshing is finished
	//
	// it is assumed that geoFilename is in the directory WorkDir or an absolute path has to be provided
	// all other files will be stored in WorkDir directory
	//
	// examples for call:
	// CreateMeshWithNETGEN(mystr("geoTest.geo"),mystr("meshTest.txt"),mystr("\\Program Files (x86)\\Netgen-4.9.13_x64\\bin\\"),
	//						mystr("C:\\Users\\AK190974\\Documents\\cpp"),mystr("Neutral Format"), coarse);
	//
	// possible values of optional parameter CONST CHAR* COARSITY are: { "verycoarse" | "coarse" | "moderate" | "fine" | "veryfine" }
	//
	// DR 2012-02-24



	mystr batFilename("startNETGENmeshing.bat");
	batFilename = WorkDir + mystr("\\") + batFilename;

	// ############################################################
	//			write bat-File		
	ofstream f;
	f.open(batFilename);
	if(!f.is_open())
	{
		return false;
	}

	//$ PG 2012-3-5: new version
	f << "cd " << WorkDir << "\n";
	f << "echo 0 > wait.txt \n";
	f << "\"" << NETGENPath << "netgen.exe\" -geofile=" << geoFilename << " -meshfile=" << meshFilename << " -" << coarsity << " -meshfiletype=\""<< NETGENMeshFileType <<"\" -batchmode \n";
	f << "echo 1 > wait.txt \n";
	f << "exit\n";
	f.close();

	// ############################################################
	//			execute bat-File (generating mesh)
	system(mystr("START ")+batFilename);

	// ############################################################
	//			wait until meshing is finished
	int goon = 0;
	while (!goon)
	{
		{
			Sleep(100);
			ifstream wait(WorkDir + mystr("\\wait.txt"));
			int test;
			wait >> test;
			if (test == 1) goon = 1;
		}
	}

	return true;
}

// generates geo-file, for use in NETGEN, to define geometry of a rotor
bool ExportAsGeoFileCylinders(const char * filename, TArray<Vector3D> PointsLeft, TArray<Vector3D> PointsRight, TArray<int> Material, double maxh, const char * rotor_name, const char * title, const char * comment)
{
	//	writes geo-file for generating geometry of a rotor in NETGEN (see "csg" in NETGEN docu) 
	//	rotor axis is z-axis. plane to define cross-section of rotor is z-x-plane
	//	rotor can be composed of the following segments: cylinder and cutted cones
	//	each segment is defined by 2 lines: outer surface and inner surface
	//	each line i is defined by left point pil and right point pir of line
	//
	//	e.g. rotor with 2 segments:	
	//	2 lines, p1l-p1r and p3l-p3r = 4 points to define outer surface
	//	2 lines, p2l-p2r and p4l-p4r = 4 points to define outer surface
	//					
	//				x-axis					p3l   p3r
	//					|								\____/
	//					|								|####|
	//					|			 p1l	p1r	|####|	
	//					|			/________\|####| p4r			
	//					| 		|##############|/
	//					|			|#########|\p4l
	//					|	p2l	|#########|  
	//					|		 \|#########|/p2r
	//	z-axis	-.-.-.-.-.-.-.-.-.-.-.-.-.-.->
	//
	//	input:	filename				name of the geo-file, e.g. "rotor.geo"
	//					PointsLeft			coordinates of points defining the geometry (segments) of the rotor (x,0,z)
	//					PointsRight			coordinates of points defining the geometry (segments) of the rotor (x,0,z)
	//					Material				material number of segments, if material number == 0 than this segment is air (=cut away) 
	//						ATTENTION:		material number is just passed to NETGEN correctly if Material contains all(!) integers from 1 to m_max
	//					maxh						maximum size of mesh element (see docu of NETGEN for "more" info)
	//					rotor_name, title, comment	strings that will be written in the geo-file
	//
	//
	//	output:	geo-file "filename" to use in NETGEN 
	//
	//	ATTENTION: it is assumed, that:
	//	lines with odd numbers 2*i-1 define outer surface --> material number != 0
	//	lines with even numbers 2*i define corresponding inner surface
	//	z-coordinates of corresponing lines have to be equal
	//	if it is a full segment (no inner surface), than a dummy inner surface with x-coordinate == 0 has to be provided
	//
	//	DR 2012-04

	int n = PointsLeft.Length();
	if(n != PointsRight.Length())
	{
		(*global_uo) << "ERROR: ExportAsGeoFileCylinders: PointsLeft.Length() != PointsRight.Length()\n";
		return false;
	}
	if(n != Material.Length())
	{
		(*global_uo) << "ERROR: ExportAsGeoFileCylinders: PointsLeft.Length() != Material.Length()\n";
		return false;
	}
	if(n%2)
	{
		(*global_uo) << "ERROR: PointsLeft.Length()%2 != 0\n";
		return false;
	}
	
	ofstream f;
	f.open(filename);
	if(!f.is_open())
	{
		return false;
	}

	f << "#\n";
	f << "## " << title <<"\n";
	f << "#\n";
	f << "algebraic3d\n";
	f << "\n";
	f << "# " << comment <<"\n";

	Vector3D left,right; 
	
	for(int i = 1; i<= n; i++)
	{
		left = PointsLeft(i);
		right = PointsRight(i);

		if(left.X())
		{
			// create segments
			f << "\n";
			f << "solid ";

			if(left.X() == right.X())		// cylinder
			{
				f << "segment" << i <<" = cylinder (0, 0, " << left.Z() << "; 0, 0, " << right.Z() << "; " << left.X() << ")";
			}
			else												// cone
			{
				f << "segment" << i <<" = cone (0, 0, " << left.Z() << "; " << left.X() << "; 0, 0, " << right.Z() << "; " << right.X() << ")";
			}

			if(Material(i)!=0)					// add cutting planes
			{
				f << "\n		and plane (0, 0, " << left.Z()  << "; 0, 0, -1) \n";
				//f << "		and plane (0, 0, " << right.Z() << "; 0, 0, 1) \n";
				//f << "		and plane (0, 0, 0; 1, 0, 0) \n";
				//f << "		and plane (0, 0, 0; 0, 1, 0); \n";
				f << "		and plane (0, 0, " << right.Z() << "; 0, 0, 1); \n";
			}
			else												// if an hole is generated, no cutting planes are used
			{
				f << ";\n";
			}
		}
	}

	// a part consists of a segment with material and a cut away segment
	for(int i = 1; i<= (n/2); i++)
	{
		f << "\n";
		f << "solid part" << i << " = segment" << 2*i-1 ;

		// just add cut away segment if radius != 0
		if(PointsLeft(2*i).X()!=0)
		{
			f << " and not segment" << 2*i;
		}
		f <<  " -bc=" << i << "; \n";
	}

	//// compose rotor of all parts
	//f << "\n";
	//f << "solid rotor = ";
	//for(int i = 1; i<= (n/2); i++)
	//{
	//	f << "part" << i;
	//	if(i==(n/2))
	//	{
	//		f << ";";
	//	}
	//	else
	//	{
	//		f << " or ";
	//	}
	//}

	//f << "\n";
	//f << "tlo rotor -col=[1,0,0];";

	//// set all parts to tlo
	//f << "\n";
	//for(int i = 1; i<= (n/2); i++)
	//{
	//	double col = 1./i;
	//	f << "tlo part" << i << " -col=[" << col << ",0,0];\n";
	//}

	// compose parts of equal materials
	TArray<int> Mat_tmp;
	Mat_tmp.Flush();
	int mat_max=0;	// maximum material number

	for(int i =1; i<=n; i++)
	{
		if(Material(i))					// ignore Material(i)==0
		{
			if(!(Mat_tmp.Find(Material(i))))		// Material(i) is not yet part of Mat_tmp
			{
				Mat_tmp.Add(Material(i));
				if(Material(i)>mat_max) {mat_max = Material(i);}
			}
		}
	}
	// all the material numbers are stored in Mat_tmp

	int createQuarterModel = 0;
	int first = 1;
	for(int m=1; m<=Mat_tmp.Length(); m++)
	{
		f << "\n";
		if(createQuarterModel) 
		{		
			f << "solid partMaterial" << Mat_tmp(m) << " = (";
		}
		else 
		{		
			f << "solid partMaterial" << Mat_tmp(m) << " = ";
		}

		first = 1;
		for(int j=1; j<=n-1; j=j+2)
		{
			if(Material(j) == Mat_tmp(m))
			{
				if(!first)
				{
					f << " or";
				}
				f << " part" << (j+1)/2 ;
				first = 0;
			}
		}
		if(createQuarterModel) 
		{	
			f << ") and plane (0, 0, 0; -1, 0, 0) and plane (0, 0, 0; 0, -1, 0) -maxh=" << maxh << ";\n";
		}
		else
		{
			f << " -maxh=" << maxh << ";\n";
		}
	}

	// set partMaterial to tlo
	f << "\n";

	int forceOrderOfMaterials = 1;		// creates a warning if it is not possible to pass material numbers correctly to NETGEN
	int m_tmp;
	double col = 0;

	for(int m=1; m<=mat_max; m++)
	{
		m_tmp = Mat_tmp.Find(m);
		if((!m_tmp)&&forceOrderOfMaterials)
		{
			(*global_uo) << "WARNING: ExportAsGeoFileCylinders: some material number is missing, so material numbers in geo-file are not correct!!! \n";
			m--;
			forceOrderOfMaterials = 0;	// further processing without warnings
		}
		else
		{
			if(m_tmp)	// material exists in Mat_tmp
			{
				col = 1./m;
				f << "tlo partMaterial" << Mat_tmp(m_tmp) << " -col=[" << col << ",0,"<< 1-col <<"];\n";
			}
		}
	}

	return true;
}


bool ExportAsMeshedCylPartEDCtxtFile(const char * filename, TArray<Vector3D> PointsLeft, TArray<Vector3D> PointsRight, TArray<int> Material, const char * rotor_name, const char * title, const char * comment)
{
	//	writes txt-file for model file to generate mesh of a rotor
	//	rotor axis is z-axis. plane to define cross-section of rotor is z-x-plane
	//	rotor can be composed of the following segments: cylinder
	//	each segment is defined by 2 lines: outer surface and inner surface
	//	each line i is defined by left point pil and right point pir of line
	//
	//	e.g. rotor with 2 segments:	
	//	2 lines, p1l-p1r and p3l-p3r = 4 points to define outer surface
	//	2 lines, p2l-p2r and p4l-p4r = 4 points to define outer surface
	//					
	//				x-axis					p3l   p3r
	//					|								\____/
	//					|								|####|
	//					|			 p1l	p1r	|####|	
	//					|			/________\|####| p4r			
	//					| 		|##############|/
	//					|			|#########|\p4l
	//					|	p2l	|#########|  
	//					|		 \|#########|/p2r
	//	z-axis	-.-.-.-.-.-.-.-.-.-.-.-.-.-.->
	//
	//	input:	filename				name of the txt-file, e.g. "EDCMesh.txt"
	//					PointsLeft			coordinates of points defining the geometry (segments) of the rotor (x,0,z)
	//					PointsRight			coordinates of points defining the geometry (segments) of the rotor (x,0,z)
	//					Material				material number of segments, if material number == 0 than this segment is air (=cut away)
	//
	//	output:	txt-file "filename" to be read by HOTINT
	//
	//	ATTENTION: it is assumed, that:
	//	lines with odd numbers 2*i-1 define outer surface --> material number != 0
	//	lines with even numbers 2*i define corresponding inner surface
	//	z-coordinates of corresponing lines have to be equal
	//	if it is a full segment (no inner surface), than a dummy inner surface with x-coordinate == 0 has to be provided
	//	right points differ from left points only w.r.t. z-coordinate (only straight cylinders and no cones!), x- and y-coordinate of right points are ignored

	int n = PointsLeft.Length();
	if(n != PointsRight.Length())
	{
		//GetMBS()->UO(UO_LVL_err) << "ERROR: ExportAsGeoFileCylinders: PointsLeft.Length() != PointsRight.Length()\n";
		return false;
	}
	if(n != Material.Length())
	{
		//GetMBS()->UO(UO_LVL_err) << "ERROR: ExportAsGeoFileCylinders: PointsLeft.Length() != Material.Length()\n";
		return false;
	}
	if(n%2)
	{
		//GetMBS()->UO(UO_LVL_err) << "ERROR: PointsLeft.Length()%2 != 0\n";
		return false;
	}
	
	ofstream f;
	f.open(filename);
	if(!f.is_open())
	{
		return false;
	}

	f << "%%%%%%%%%%%%%%%%%%%" << title <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
	f << "%% " << comment <<"\n";
	f << "\n";
	f << rotor_name <<"\n";
	f << "{\n";
	f << " axis = [ 0., 0., 1. ]\n";
	f << " base = [ 0., 0., 0. ]\n";
	f << " quadrants = 1\n";

	for(int i = 1; i<= (n/2); i++)
	{
		int j = 2*i-1;
		f << "\n";
		f << "  Cylinder" << i << "\n";
		f << "  {\n";

		f << "    outer_radius         = " << PointsLeft(j).X()		<< "\n";
		f << "    inner_radius     = " << PointsLeft(j+1).X() << "\n";
		f << "    axis_pos_bottom  = "	<< PointsLeft(j).Z()		<< "\n";
		f << "    axis_pos_top     = "		<< PointsRight(j).Z()	<< "\n";

		f << "    material_number  = "		<< Material(j)	<< "\n";
		f << "    domain_number    = "		<< 1						<< "\n";
		f << "    rgbcolor         =  [ 0, 1, 0 ]	\n";

		f << "    radial_mesh_divisions  = 1 \n";
		f << "    angular_mesh_divisions = 1	\n";
		f << "    outward_refinement     = 0 \n";
		f << "    axial_mesh_divisions   = 1	\n";
		f << "  }\n";
	}
	f << "}\n";

	return true;
}
