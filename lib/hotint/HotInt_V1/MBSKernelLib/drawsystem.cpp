//#**************************************************************
//#
//# filename:             drawsystem.cpp
//# 
//# author:               Gerstmayr Johannes
//#
//# generated:						September 2010
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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

#include "mbs.h"
#include "element.h"
#include "graphicsConstants.h"
#include "node.h"
#include "sensors.h"

void MultiBodySystem::DrawSystem()
{
	if (GetProhibitRedraw()) return;

	GetRC()->PrintTextStruct(0,-1," HOTINT");

	TMStartTimer(2);

	GetRC()->SetBackgroundColor(1.,1.,1.);
	const Vector& x = GetDrawVector();

	GetRC()->SetTextColor(0.,0.,0.);

	if (GetCenterObject() != 0)
	{
		Vector3D bodypos = GetElement(GetCenterObject()).GetRefPosD();
		bodypos += GetCenterObjectOffset();
		GetRC()->SetCenterOffset((float)bodypos.X(),(float)bodypos.Y(),(float)bodypos.Z());
	}

	char str[100];
	sprintf(str,"    time=%4.4g s ", GetDrawTime());
	//sprintf(str,"    time=%3.3g days ", GetDrawTime()/(24*3600));

	//if (!GetIOption(116)) 
	{
		GetRC()->PrintTextStruct(3,-1,str);
	}

	if (GetIOption(106))
	{
		//Draw Origin ...
		SetColor(colgrey3);
		double s = GetDOption(103);
		Vector3D v1(-0.3*s, 0,0);
		Vector3D v2( s, 0,0);
		Vector3D v3( 0,-0.3*s,0);
		Vector3D v4( 0, s,0);
		Vector3D v5( 0,0,-0.3*s);
		Vector3D v6( 0,0, s);
		double d = GetDOption(114); //global line thickness
		MyDrawLine(v1,v2,d);
		MyDrawLine(v3,v4,d);
		MyDrawLine(v5,v6,d);

		GetRC()->PrintText3D((float)v2.X(), (float)v2.Y(), (float)v2.Z(), "X0");
		GetRC()->PrintText3D((float)v4.X(), (float)v4.Y(), (float)v4.Z(), "Y0");
		GetRC()->PrintText3D((float)v6.X(), (float)v6.Y(), (float)v6.Z(), "Z0");
	}

	//!AD 2012-12-14 new scene background use the entries in the options used for option Grid...
	if(GetIOption(126)) 
	{
		// UNDER CONSTRUCTION !
		// 6 options for size of RECTANGLES
		double refpos_x = GetDOption(107);
		double refpos_y = GetDOption(108);
		double refpos_z = GetDOption(109);
		double size_of_plane_x = GetDOption(112);
		double size_of_plane_y = GetDOption(113);
		double size_of_plane_z = GetDOption(142); 

		// 5 options for GRID DISCRETIZATION
		double size_of_grid_x = GetDOption(110);
		double size_of_grid_y = GetDOption(111);
		double size_of_grid_z = GetDOption(143); 
		int lines_crosses_dots = 1;         // HARDCODED: ALWAYS LINES ( code for dots stays in )
		int max_gridlines = 200;            // HARDCODED ON DEMAND 

		// 10 options for COLOR
		float color_p1r = (float) GetDOption(133); float color_p1g = (float) GetDOption(134); float color_p1b = (float) GetDOption(135); 
		float color_p2r = (float) GetDOption(136); float color_p2g = (float) GetDOption(137); float color_p2b = (float) GetDOption(138); 
		float color_p3r = (float) GetDOption(139); float color_p3g = (float) GetDOption(140); float color_p3b = (float) GetDOption(141); 
		float background_transparency_factor = GetDOption(144);

		// IOption 126 must be nonzero to activate
		double global_line_thickness = GetDOption(114); 
		double global_point_size = GetDOption(117);


		// account for maximum number of grid lines
		int nr_of_gridlines_x = (int) abs(size_of_plane_x/size_of_grid_x) +0.5; // excluding the 0-th line
		if (nr_of_gridlines_x > max_gridlines) { nr_of_gridlines_x = max_gridlines; size_of_grid_x = abs(size_of_plane_x / max_gridlines); }
		int nr_of_gridlines_y = (int) abs(size_of_plane_y/size_of_grid_y) +0.5; // excluding the 0-th line
		if (nr_of_gridlines_y > max_gridlines) { nr_of_gridlines_y = max_gridlines; size_of_grid_y = abs(size_of_plane_y / max_gridlines); }
		int nr_of_gridlines_z = (int) abs(size_of_plane_z/size_of_grid_z) +0.5; // excluding the 0-th line
		if (nr_of_gridlines_z > max_gridlines) { nr_of_gridlines_z = max_gridlines; size_of_grid_z = abs(size_of_plane_z / max_gridlines); }


		Vector3D p0(refpos_x, refpos_y, refpos_z); // intersection point of the planes
		int active_flag = GetIOption(126);
		// draw the planes
		if(active_flag & 1)
		{
			pCurrentRC->glColor4f(color_p1r, color_p1g, color_p1b, background_transparency_factor);
			TArray<Vector3D> ptsxy; ptsxy.Add(p0); ptsxy.Add(p0+Vector3D(size_of_plane_x,0.,0.)); ptsxy.Add(p0+Vector3D(size_of_plane_x,size_of_plane_y,0.)); ptsxy.Add(p0+Vector3D(0.,size_of_plane_y,0.));
			DrawQuad( ptsxy(1), ptsxy(2), ptsxy(3), ptsxy(4) );
			DrawPolygonOutline( ptsxy, global_line_thickness );
		}
		if(active_flag & 2)
		{
			pCurrentRC->glColor4f(color_p2r, color_p2g, color_p2b, background_transparency_factor);
			TArray<Vector3D> ptsxz; ptsxz.Add(p0); ptsxz.Add(p0+Vector3D(size_of_plane_x,0.,0.)); ptsxz.Add(p0+Vector3D(size_of_plane_x,0.,size_of_plane_z)); ptsxz.Add(p0+Vector3D(0.,0.,size_of_plane_z));
			DrawQuad( ptsxz(1), ptsxz(2), ptsxz(3), ptsxz(4) );
			DrawPolygonOutline( ptsxz, global_line_thickness );
		}
		if(active_flag & 4)
		{
			pCurrentRC->glColor4f(color_p3r, color_p3g, color_p3b, background_transparency_factor);
			TArray<Vector3D> ptsyz; ptsyz.Add(p0); ptsyz.Add(p0+Vector3D(0.,size_of_plane_y,0.)); ptsyz.Add(p0+Vector3D(0.,size_of_plane_y,size_of_plane_z)); ptsyz.Add(p0+Vector3D(0.,0.,size_of_plane_z));
			DrawQuad( ptsyz(1), ptsyz(2), ptsyz(3), ptsyz(4) );
			DrawPolygonOutline( ptsyz, global_line_thickness );
		}

		// account for negative size ( extends in negative axis direction )
		Vector3D ex(1.,0.,0.); if(size_of_plane_x<0) ex.X() = -1.;
		Vector3D ey(0.,1.,0.); if(size_of_plane_y<0) ey.Y() = -1.;
		Vector3D ez(0.,0.,1.); if(size_of_plane_z<0) ez.Z() = -1.;

		Vector3D old_color_line = colline;
		colline = Vector3D(0.8,0.8,0.8);
		if(lines_crosses_dots == 1) // GRID LINES
		{
			Vector3D p1,p2;
			for(int ix=0; ix<=nr_of_gridlines_x; ix++)
			{
				p1 = p0 + ex*size_of_grid_x*ix;
				p2 = p0 + ex*size_of_grid_x*ix + Vector3D(0.,size_of_plane_y,0.);
				if (active_flag & 1) 
					MyDrawLine(p1, p2, global_line_thickness); // this is a line on the XY plane
				p2 = p0 + ex*size_of_grid_x*ix + Vector3D(0.,0.,size_of_plane_z);
				if (active_flag & 2) 
					MyDrawLine(p1, p2, global_line_thickness); // this is a line on the XZ plane
			}
			for(int iy=0; iy<=nr_of_gridlines_y; iy++)
			{
				p1 = p0 + ey*size_of_grid_y*iy;
				p2 = p0 + ey*size_of_grid_y*iy + Vector3D(size_of_plane_x,0.,0.);
				if (active_flag & 1) 
					MyDrawLine(p1, p2, global_line_thickness); // this is a line on the XY plane
				p2 = p0 + ey*size_of_grid_y*iy + Vector3D(0.,0.,size_of_plane_z);
				if (active_flag & 4) 
					MyDrawLine(p1, p2, global_line_thickness); // this is a line on the YZ plane
			}
			for(int iz=0; iz<=nr_of_gridlines_z; iz++)
			{
				p1 = p0 + ez*size_of_grid_z*iz;
				p2 = p0 + ez*size_of_grid_z*iz + Vector3D(size_of_plane_x,0.,0.);
				if (active_flag & 2) 
					MyDrawLine(p1, p2, global_line_thickness); // thie is a line on the XZ plane
				p2 = p0 + ez*size_of_grid_z*iz + Vector3D(0.,size_of_plane_y,0.);
				if (active_flag & 4) 
					MyDrawLine(p1, p2, global_line_thickness); // this is a line on the YZ plane
			}
			colline = old_color_line;
		}
		/* else if(lines_crosses_dots == 2) { ; } // crosses not implemented yet */
		else // GRID DOTS
		{
			Vector3D p1;
			ChooseColor((float)colline[0],(float)colline[1],(float)colline[2]);
			GetRC()->ChoosePointSize((float)global_point_size);

			GetRC()->glBeginPoints();
			if (active_flag & 1)
			{
				for(int ix=0; ix<=nr_of_gridlines_x; ix++)
				{
					for(int iy=0; iy<=nr_of_gridlines_y; iy++)
					{
						p1 = p0 + ex*size_of_grid_x*ix + ey*size_of_grid_y*iy;
						GetRC()->glVertex((float)p1.X(), (float)p1.Y(), (float)p1.Z());				
					}
				}
			}
			if (active_flag & 2)
			{
				for(int ix=0; ix<=nr_of_gridlines_x; ix++)
				{
					for(int iz=0; iz<=nr_of_gridlines_z; iz++)
					{
						p1 = p0 + ex*size_of_grid_x*ix + ez*size_of_grid_z*iz;
						GetRC()->glVertex((float)p1.X(), (float)p1.Y(), (float)p1.Z());				
					}
				}
			}
			if (active_flag & 4)
			{
				for(int iy=0; iy<=nr_of_gridlines_y; iy++)
				{
					for(int iz=0; iz<=nr_of_gridlines_z; iz++)
					{
						p1 = p0 + ey*size_of_grid_y*iy + ez*size_of_grid_z*iz;
						GetRC()->glVertex((float)p1.X(), (float)p1.Y(), (float)p1.Z());				
					}
				}
			}
			GetRC()->glEnd();
		}
	}

	int old_transp_mode = GetTransparency();
	//+++++transparency is turned on here
	SetTransparency(0);

	//drawing preprocessing: for interpolation of stress, etc.

	for (int i=1; i<=nodes.Length(); i++) 
	{
		nodes(i)->SetDrawTemp(0);
		nodes(i)->SetDrawTempCnt(0);
	}

	for (int i=1; i<=elements.Length(); i++) 
	{
		if (GetElement(i).IsType(TCMSflag))
		{
			//UO() << "cms=" << i << ", nn=" << GetElement(i).NNodes() << "\n";
			for (int j=1; j <= GetElement(i).NCMSNodes(); j++)
			{
				GetElement(i).GetNode(j).SetDrawTemp(0);
				GetElement(i).GetNode(j).SetDrawTempCnt(0);
			}
		}
	}

	for (int i=1; i<=elements.Length(); i++) 
	{
		GetElementPtr(i)->DrawElementPreProc();
	}

	for (int i=1; i<=auxelements.Length(); i++) 
	{
		GetAuxElementPtr(i)->DrawElementPreProc();
	}

	//perform averaging of Finite element solution at nodes, using the number of values added to the node
	for (int i=1; i<=nodes.Length(); i++) 
	{
		double val = nodes(i)->GetDrawTempCnt();
		if (val == 0) val = 1;

		nodes(i)->GetDrawTemp() /= val;
	}

	for (int i=1; i<=elements.Length(); i++) 
	{
		if (GetElement(i).IsType(TCMSflag))
		{
			Element& el = GetElement(i);
			for (int j=1; j <= el.NCMSNodes(); j++)
			{
				double val = el.GetNode(j).GetDrawTempCnt();
				if (val == 0) val = 1;

				el.GetNode(j).GetDrawTemp() /= val;
			}
		}
	}

	//determine range of finite element contour colors
	if (GetActualPostProcessingFieldVariable() != NULL)
	{
		SetFlagComputeMinMaxFECol(1);
		if(GetIOption(102))		// auto adjust range
		{
			GetFEmincol() = 1e30;
			GetFEmaxcol() = -1e30;
		}
		DrawBodies();
		//global_uo << "aha, femincol = " << GetFEmincol() << ", max = " << GetFEmaxcol() << "\n";
		SetFlagComputeMinMaxFECol(0);
	}

	//compute loads, nodes, sensors, constraints (different order depending if transparent or not)
	if (GetIOption(142)) //draw loads
	{
		DrawLoads();
	}

	if (GetIOption(145)) //draw nodes
	{
		DrawNodes();
	}

	if (!GetIOption(130)) //sensors not transparent
	{
		DrawSensors();
	}
	if (!GetIOption(129)) //constraints not transparent
	{
		DrawConstraints();
	}
	if (!GetIOption(128)) //bodies not transparent
	{
		DrawBodies();
	}
	DrawElementAdd();
	//++++ up to here transparency is set to global value

	//++++ from here transparency is turned of
	SetTransparency(1);
	if (GetIOption(130)) //sensors transparent
	{
		DrawSensors();
	}
	if (GetIOption(129)) //constraints transparent
	{
		DrawConstraints();
	}
	if (GetIOption(128)) //bodies transparent
	{
		DrawBodies();
	}
	SetTransparency(0); //legend not transparent


	if (GetActualPostProcessingFieldVariable() != NULL && !GetIOption(141))
	{
		sprintf(str,"                       ");
		int toff = 3;
		int h2 = GetIOption(103);

		//*YV: new version of printing the title of the legend is based on the names of the variables
		bool flagMillimeter = GetIOption(143) == 1;		// whether "mm" or "m" should be printed
		char * dim_str = "";
		switch(GetActualPostProcessingFieldVariable()->GetDimensionality())
		{
		case FieldVariableDescriptor::FVD_length:										dim_str = flagMillimeter ? "mm" : "m"; break;
		case FieldVariableDescriptor::FVD_velocity:									dim_str = flagMillimeter ? "mm/s" : "m/s"; break;
		case FieldVariableDescriptor::FVD_acceleration:							dim_str = flagMillimeter ? "mm/s²" : "m/s²"; break;
		case FieldVariableDescriptor::FVD_force:										dim_str = "N"; break;
		case FieldVariableDescriptor::FVD_force_per_length:					dim_str = flagMillimeter ? "N/mm" : "N/m"; break;
		case FieldVariableDescriptor::FVD_force_per_length_square:	dim_str = flagMillimeter ? "N/mm²" : "N/m²"; break;
		case FieldVariableDescriptor::FVD_force_length:							dim_str = flagMillimeter ? "N*mm" : "N*m"; break;
		case FieldVariableDescriptor::FVD_1_per_length:							dim_str = flagMillimeter ? "1/mm" : "1/m"; break;
		case FieldVariableDescriptor::FVD_density:									dim_str = flagMillimeter ? "kg/mm³" : "kg/m³"; break;
		case FieldVariableDescriptor::FVD_pressure:									dim_str = flagMillimeter ? "N/mm²" : "N/m²"; break;
		}
		if(GetActualPostProcessingFieldVariable()->GetDimensionality() ==  FieldVariableDescriptor::FVD_none)
			sprintf(str, " %s max", GetActualPostProcessingFieldVariable()->GetTextualIdentifier().c_str());
		else
			sprintf(str, " %s max, %s", GetActualPostProcessingFieldVariable()->GetTextualIdentifier().c_str(), dim_str);
		pCurrentRC->PrintTextStruct(5+toff,-1,str);
		if(GetActualPostProcessingFieldVariable()->GetDimensionality() ==  FieldVariableDescriptor::FVD_none)
			sprintf(str, " %s min", GetActualPostProcessingFieldVariable()->GetTextualIdentifier().c_str());
		else
			sprintf(str, " %s min, %s", GetActualPostProcessingFieldVariable()->GetTextualIdentifier().c_str(), dim_str);
		pCurrentRC->PrintTextStruct(7+h2+toff,-1,str); //18
		//sprintf(str,"     %6.5gN/m²   ", TImincol);
		//pCurrentRC->PrintTextStruct(6+h2+toff,-1,str); //19
		DrawLegend(6+toff);

		for (double i=0.; i <= h2; i++)
		{
			double val = i/(double)h2;
			if (GetIOption(108) && val >= 0 && val <= 1.) 
			{
				//val = 1.-pow(1.-val,10.);
				val = 1.-pow(1.-val,5);
			}
			val = (GetFEmincol()-GetFEmaxcol())*val+GetFEmaxcol();
			//sprintf(str,"            %6.3e   ", val); 
			sprintf_s(str,"            %6.*e   ", GetIOption(160),val); //$ AD 2011-04-08: precision of legend adjustable

			pCurrentRC->PrintTextStruct(6+toff+(int)i,-1,str);
		}
	}
	SetTransparency(old_transp_mode);



	TMStopTimer(2);

}

//AD
void MultiBodySystem::DrawSystem2D()
{
	p2DDrawWindow->Reset();
	for (int i=1; i<=elements.Length(); i++) 
	{
		if (elements(i)->IsType(TConstraint) && elements(i)->DrawElementFlag())
		{
			elements(i)->DrawElement2D();
		}
	}
}

void MultiBodySystem::DrawElementAdd() 
{	
	int use_cutting_planes = UseCuttingPlanes();

	//draw only elements which do not belong to an element!
	for (int i=1; i <= drawelements.Length(); i++)
	{
		if(
			!drawelements(i)->GetElnum() &&
			(!use_cutting_planes || !GetOptions()->ViewingOptions()->CuttingPlaneCutGround() || CuttingPlanesAllow(drawelements(i)->GetRefPosD()))
			)
		{
			drawelements(i)->DrawYourself();
		}
	}	
}

void MultiBodySystem::DrawBodies()
{
	//drawbodiescnt++;
	//UO() << "cnt=" << drawbodiescnt << "\n";

	int use_cutting_planes = UseCuttingPlanes();
	int num_particles = 0.;

	for (int i=1; i<=elements.Length(); i++) 
	{
		if (elements(i)->DrawElementFlag())
		{
			if (!(elements(i)->IsType(TConstraint)))
			{
				if (!use_cutting_planes || elements(i)->GetAltShape() || !GetOptions()->ViewingOptions()->CuttingPlaneCutBodies() || CuttingPlanesAllow(elements(i)->GetRefPosD()))
				{
					// Draw body number if IOption(123) for all bodies
					if ((elements(i)->IsType(TBody) && GetIOption(123)) && !IsFlagComputeMinMaxFECol())
					{
						//does work for 2D-bodies ?
						Vector3D p = elements(i)->GetRefPosD();
						//global_uo << "p" << i << "=" << p << "\n";
						//Vector3D p = GetPosD(Vector3D(0));
						char text[40];
						sprintf(text,"element %d", i);

						GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
					}

					// 125 = show local frame
					if ((elements(i)->IsType(TBody) && GetIOption(125)) && !IsFlagComputeMinMaxFECol())
					{
						elements(i)->DrawElementLocalFrame();
					}

					// Draw element (also done in ComputeMinMaxFECol case for finite elements)
					if (!elements(i)->GetAltShape())
					{
						if (!(IsFlagComputeMinMaxFECol() && elements(i)->IsRigid()))
						{
							int no_surfaceupdate = GetIOption(119);
							if (use_cutting_planes && no_surfaceupdate==0 && elements(i)->IsFiniteElement() && elements(i)->Dim()==3) // probably it is Hex/Tet)
							{
								//$ YV 2012-12-11: treatment of outer faces in finite elements was placed to finite elements
								/*
								FiniteElement3D * fe = dynamic_cast<FiniteElement3D*>(elements(i));
								char copyofouterfaces = fe->Outer_Face();

								fe->SetOuterFacesCuttingPlane();
								elements(i)->DrawElement();

								fe->Outer_Face() = copyofouterfaces;
								*/
								elements(i)->DrawElement();
							}
							else
							{
								if (GetOptions()->PostProcOptions()->BodiesShowVelocityVector())
								{
									Element* el = elements(i);
									if (!GetOptions()->PostProcOptions()->BodiesVelocityVectorJustForParticles())
									{
										el->DrawElementVelocityVector();
									}
								}

								elements(i)->DrawElement();
							}
						}
					}
					// draw AltShape elements
					else if (!IsFlagComputeMinMaxFECol())
					{
						if (GetOptions()->PostProcOptions()->BodiesShowVelocityVector())
						{
							Element* el = elements(i);
							if (!GetOptions()->PostProcOptions()->BodiesVelocityVectorJustForParticles())
							{
								el->DrawElementVelocityVector();
							}
						}

						elements(i)->DrawElementAdd();
					}
				}
			}
		}
	}
	for (int i=1; i<=auxelements.Length(); i++) 
	{
		//$ MaSch 2013-10-28: moved drawing of SPH particles from elements to auxelements
		if (!auxelements(i)->IsType(TParticle) || GetIOption(155) <= 1 || ++num_particles % GetIOption(155) == 0)  //GetIOption(155) corresponds to PostProcOptions.Bodies.Particles.draw_every_nth
		{
			if (!use_cutting_planes || CuttingPlanesAllow(auxelements(i)->GetRefPosD()))  
			{
				// Draw auxelement body number if IOption(123) for all bodies
				if ((auxelements(i)->IsType(TBody) && GetIOption(123)) && !IsFlagComputeMinMaxFECol())
				{
					//does work for 2D-bodies ?
					Vector3D p = auxelements(i)->GetRefPosD();
					//global_uo << "p" << i << "=" << p << "\n";
					//Vector3D p = GetPosD(Vector3D(0));
					char text[40];
					sprintf(text,"aux-element %d", i);

					GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
				}
				// Draw auxelement (also done in ComputeMinMaxFECol case for finite elements)
				if (!auxelements(i)->GetAltShape())
				{
					if (!(IsFlagComputeMinMaxFECol() && auxelements(i)->IsRigid()))
					{
						if (GetOptions()->PostProcOptions()->BodiesShowVelocityVector() && auxelements(i)->IsType(TParticle))
						{
							auxelements(i)->DrawElementVelocityVector();
						}
						auxelements(i)->DrawElement();
					}
				}
				// draw AltShape auxelements
				else if (!IsFlagComputeMinMaxFECol())
				{
					if (GetOptions()->PostProcOptions()->BodiesShowVelocityVector() && auxelements(i)->IsType(TParticle))
					{
						auxelements(i)->DrawElementVelocityVector();
					}
					auxelements(i)->DrawElementAdd();
				}
			}
		}
	}
}

void MultiBodySystem::DrawConstraints()
{
	int use_cutting_planes = UseCuttingPlanes();

	if (GetIOption(121)) //do not draw constraints/joints
	{
		for (int i=1; i<=elements.Length(); i++) 
		{
			if (elements(i)->IsType(TConstraint) && elements(i)->DrawElementFlag())
			{
				if (!use_cutting_planes || CuttingPlanesAllow(elements(i)->GetRefPosD()))
				{
					if ((elements(i)->IsType(TConstraint) && GetIOption(124)))
					{
						//does work for 2D-bodies ?
						Vector3D p = elements(i)->GetRefPosD();
						//Vector3D p = GetPosD(Vector3D(0));
						char text[40];
						sprintf(text,"element %d", i);

						GetRC()->PrintText3D((float)p.X(),(float)p.Y(),(float)p.Z(),text);
					}

					if (!elements(i)->GetAltShape())
						elements(i)->DrawElement();

					elements(i)->DrawElementAdd();
				}
			}
		}
	}
}

void MultiBodySystem::DrawSensors()
{
	int use_cutting_planes = UseCuttingPlanes();
	Vector3D pos; 
	int np;

	if (GetIOption(127))
	{
		for (int i=1; i <= NSensors(); i++)
		{
			if(GetSensor(i).GetNumberOfDrawingPositions() > 0)
				if (!use_cutting_planes || CuttingPlanesAllow(GetSensor(i).GetDrawPosition(1)))
				{
					mystr label = "S" + mystr(i);
					GetSensor(i).Draw(label);
				}
		}
	}
}

void MultiBodySystem::DrawLoads()
{
	int use_cutting_planes = UseCuttingPlanes();

	for (int i=1; i<=elements.Length(); i++) 
	{
		if (!use_cutting_planes || CuttingPlanesAllow(elements(i)->GetRefPosD()))
		{
			if ((elements(i)->IsType(TBody)))
			{
				for (int j = 1; j <= elements(i)->NLoads(); j++)
				{
					elements(i)->GetLoad(j).DrawLoad(this,elements(i));
				}
			}
		}
	}
}

void MultiBodySystem::DrawNodes()
{
	int use_cutting_planes = UseCuttingPlanes();

	for (int i=1; i<=nodes.Length(); i++) 
	{
		if (!use_cutting_planes || CuttingPlanesAllow(nodes(i)->GetPosD())) //$ AD 2011-07-28: Cutplane for Nodes 
		{
			nodes(i)->DrawNode();
		}
	}
}


int MultiBodySystem::UseCuttingPlanes()
{
	int rv = 0;

	//$ PG 2012-3-6:[ Cutting plane may now be optionally handled by OpenGL (if GetOptions()->ViewingOptions()->CuttingPlaneCutWholeSceneByOpenGl()  (or GetIOption(157)) is true)
	if (!GetOptions()->ViewingOptions()->CuttingPlaneCutWholeSceneByOpenGl())
	{
		int n_cutting_planes = 2;
		ncut.SetLen(n_cutting_planes);
		pcut.SetLen(n_cutting_planes);

		use_cutting_plane(1) = GetIOption(149);
		ncut(1) = Vector3D(GetDOption(125),GetDOption(126),GetDOption(127));
		pcut(1) = ncut(1)*GetDOption(128);   // multiply by distance of plane from origin

		use_cutting_plane(2) = GetIOption(156);
		ncut(2) = Vector3D(GetDOption(129),GetDOption(130),GetDOption(131));
		pcut(2) = ncut(2)*GetDOption(132);

		for (int j=1; j<=n_cutting_planes; j++)   // multiply by distance of plane from origin
		{
			if (use_cutting_plane(j))
			{
				rv = 1;
			}
		}
	}
	//$ PG 2012-3-6:] Cutting plane may now be optionally handled by OpenGL

	return rv;
}

int MultiBodySystem::CuttingPlanesAllow(const Vector3D& pos)
{
	for (int j=1; j<=use_cutting_plane.Length(); j++)
	{
		if (use_cutting_plane(j) && (pos-pcut(j))*ncut(j) >= 0)   //$ DR 2011-05-20: GetRefPos changed to GetRefPosD, now control element is drawn correctly
		{
			return 0;     //cut_element = 1;
		}
	}

	return 1;
}

FieldVariableDescriptor * MultiBodySystem::GetActualPostProcessingFieldVariable()
{
	if(indexOfActualPostProcessingFieldVariable == 0)
		return NULL;
	return &(availableFieldVariables(indexOfActualPostProcessingFieldVariable));
}

void MultiBodySystem::SetIndexOfActualPostProcessingFieldVariable(int index)
{
	assert(index <= availableFieldVariables.Length());

	// piece of code, similar to that one from the old version of SetIOption() for the case "index == 102"
	if (indexOfActualPostProcessingFieldVariable != index)
	{
		indexOfActualPostProcessingFieldVariable = index;

		if (!GetIOption(100))
			GetFEmaxcol() = -1e15;
		else
			GetFEmaxcol() = GetDOption(100);

		if (!GetIOption(101))
			GetFEmincol() = 1e15;
		else
			GetFEmincol() = GetDOption(101);
	}

	// we set the textual identifier of the present variable
	if(index == 0)
		SetTOption(107, "");
	else
		SetTOption(107, GetActualPostProcessingFieldVariable()->GetTextualIdentifier());
}