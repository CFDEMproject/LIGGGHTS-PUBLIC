//#**************************************************************
//#
//# filename:             numnlsys_options.cpp
//#
//# author:               Gerstmayr Johannes, Peter Gruber
//#
//# date:									Feb 02 2012
//# description:          definition of [I|D|T]Option-access
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

#include "ioincludes.h"

#include <assert.h>
#include <memory.h>

#include <math.h>

#include "tarray.h"    
#include "mystring.h"  
#include "femath.h"  


//options:
// 1  ... 49 ... computational options
// 80 ... 99 ... file processing options
// 100...199 ... postprocessing options
// 200...249 ... OPENgl graphics options
// 250...299 ... Robot options (old)


void NumNLSys::DeleteOptions()
{
	delete hotint_options;
	for (int i=0; i <= 300; i++)
	{
		if (toptions[i]) delete[] toptions[i];
	}
}

const int min_value_legal_option = 80; //options < 100 will be illegal because they have moved to MBSSolSet()

void NumNLSys::InitializeOptions()
{
	//initial values for options
	for (int i = 0; i <= 300; i++)
	{
		toptions[i] = 0;
	}
	for (int i = min_value_legal_option; i <= 300; i++)
	{
		SetIOption(i,0);
		SetDOption(i,0);
	}
	InitializeOptionsAuto();
	/*
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//computation:
	//not provided in mbssolset:
	//SetIOption(  1,0); //reuse last eigenvectors
	//SetIOption( 10,0); //write mass and stiffness matrix to file ---> 	mbs->MBSSolSet().write_MK2file = 1;
	//SetIOption( 11,1); //Element-wise-Jacobian --> MBSSolSet().element_wise_jacobian
	//SetIOption( 20,1); // store mass matrix for FiniteElement3D --> mbs->MBSSolSet().store_FE_mass_matrix
	
	//SetIOption( 12,0); //Selected Model Function

	//SetTOption(  2,"---No Function---"); //Selected Model Function Name

	//SetDOption(  4,2.);			//computation_update_everyXsec 



	//from mbssolset:
	//SetIOption(  2,1); //writesolstep
	//SetIOption(  3,1); //writeresults
	//SetIOption(  4,0); //automatic step control (fulladaptive)
	//SetIOption(  5,2); //max stages
	//SetTOption(  1,"LobattoIIIA"); //tableauname; default integration method trapezoidal

	//
	//SetDOption(  1,0.001);	//max stepsize
	//SetDOption(  2,0.);			//start time
	//SetDOption(  3,10.);		//end time
	//SetDOption(  5,1e-8);		//nls_relativeaccuracy
	////SetDOption(  6,1e-6);		//nls_absoluteaccuracy: UNUSED
	//SetDOption(  7,1e-7);		//nls_numdiffepsi
	//SetDOption(  8,0.);			//	discontinuousaccuracy = 1E-4;			//accuracy for discontinuous problems (plasticity, contact, friction, ...)
	//SetDOption(  9,0.);			//	maxdiscontinuousit = 8;						//max. number of iterations for discont. problems


	//SetDOption( 10,1e-3);   //automatic stepsize: initstepsize
	//SetDOption( 11,1e-2);   //automatic stepsize: absaccuracy 
	//SetDOption( 12,1);      //automatic stepsize: relaccuracy #not used up to now!!!
	//SetDOption( 13,1e-4);   //automatic stepsize: minstepsize 
	//SetDOption( 14,0.01);   //store solution data manager every x step //# -2 == always, -1 == at max stepsize, 0 = never, x.x = at every time x.x


	////RL: parameter study
	////int:
	//SetIOption( 52,0);			//parameter study: m_check_parameter_variation, combo box flag: 0..."none", 1... "linear", 2..."quadratic" parameter variation
	////double:
	//SetDOption( 56,0.5);			//parameter study: m_parameter_variation_initial_step
	//SetDOption( 57,0.);			//parameter study: m_parameter_variation_start_value			
	//SetDOption( 58,1.);			//parameter study: m_parameter_variation_end_value		


	////static computation:
	////int:
	//SetIOption( 50,0);			//do static computation
	//SetIOption( 51,1);			//static computation: static_incloadinc			//if incloadinc successfull steps --> leads to load increment
	////double:
	//SetDOption( 50,1e-12);	//static computation: static_minloadinc     //minimal increment
	//SetDOption( 51,1.);			//static computation: static_maxloadinc			//maximum load increment, 0.1 for plasticity?
	//SetDOption( 52,1.);			//static computation: static_initloadinc		//initial load increment, 0.05 is good
	//SetDOption( 53,2.);			//static computation: static_loadinc_up			//increase load increment if success very often
	//SetDOption( 54,2.);			//static computation: static_loadinc_down		//decrease load increment if no success
	//SetDOption( 55,0.);			//static computation: spring-type regularisation parameter
	

	SetIOption(220,0);  // automatically close file after computation
	SetIOption(221,1);  // 1 .. show HOTINT window, 0 .. HOTINT-window minimized

	////+++++++++++++++++++++ Eigenvalue Computation
	//SetIOption(222, 3);    // number of eigenvalues to be computed for sparse iterative methods
	//SetIOption(223,1);  //maximal possible number of computed eigenvalues
	//SetIOption(224, 1000); // maximum number of iterations
	//SetIOption(225, 2);    // which solver to use (0..direct, 1..matlab, 2..lobpcg)
	//SetIOption(226, 0);    // number of zero eigenmodes
	//SetIOption(227, 0);    // uncheck number of zeromodes
	//SetIOption(228, 0);    // use preconditioning
	//SetDOption(225, 1e-6); // eigenvalue solver tolerance
	//SetDOption(226, 1);    // lambda for preconditioner (A + lambda M)^-1


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//file:
	////SetIOption(80, 0); //1 = always overwrite output files, 0 = append solution
	//SetIOption(81, 0); //1 = write solution file header, 0 = no header
	////SetTOption(80,"..\\..\\output\\");	//standard output directory **** save? use "..\\..\\output\\" for Release!
	//SetTOption(81,"sol.txt");						//standard output file
	//SetTOption(82,"solpar.txt");				//standard output file for parametric solution
	////SetTOption(83,"statein.txt");				//standard input file for state (initial conditions)
	////SetTOption(84,"stateout.txt");			//standard output file for state (final state of computation)
	//SetTOption(85,"testfile.mbs");			//default file name (for save, if save as has been applied, or if loaded already)
	//SetTOption(86,"");									//default directory for load/save as
	//SetTOption(87,"hotint.log");				//default file name for log file!
 // SetTOption(100,"");                 //solution file header comment

	SetTOption(101,"");									//recent file filename 1, old: 94
	SetTOption(102,"");									//recent file directory 1
	SetTOption(103,"");									//recent file filename 2
	SetTOption(104,"");									//recent file directory 2
	SetTOption(105,"");									//recent file filename 3
	SetTOption(106,"");									//recent file directory 3; old:99
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//postprocessing graphics:
	SetIOption(100,0);  //maxstress active
	SetIOption(101,0);  //minstress active
	SetIOption(102,0);  //stress/strain component, equal to drawmode, see GetDrawTypeComponents()
	//
	SetIOption(103,10); //tiling
	SetIOption(104,4); 	//redraw frequency: 0..off, 1..draw last frame, 2..100sec, 3..20sec, 4..2sec, 5..200ms, 6..50ms, 7..20ms, 8..every 10 frames, 9..every frame
	SetIOption(105,0); 	//grey mode (1=on, 0=off)
	SetIOption(106,0); 	//draw origin
	SetIOption(107,0); 	//invert isoplot colors
	SetIOption(108,0); 	//nonlinear isoplot color-bar
	SetIOption(109,1); 	//animation frames: show every i frame at animation
	SetIOption(110,1); 	//show mesh
	SetIOption(111,1); 	//show solution in mesh
	SetIOption(112,1); 	//show contact points
	SetIOption(113,1); 	//show contact body
	SetIOption(114,1); 	//show constraints --> show constraint faces ****
	SetIOption(115,1);  //start animation from beginning
	SetIOption(116,0);  //Show modes / Draw Chladni ISO-lines
	SetIOption(117,0);  //Draw Plate elements flat, only midplane (view from top only)
	SetIOption(118,0);  //Draw Stress/strain/etc. interpolated at nodes

	SetIOption(120,1);  //use degrees instead of radiant in edit dialogs for bodies and joints
	SetIOption(121,1);  //show joints
	SetIOption(122,1);  //draw 3D texts in front
	SetIOption(123,0);  //show element body numbers
	SetIOption(124,0);  //show constraint numbers
	SetIOption(125,0);  //show local body frame
	SetIOption(126,0);  //show grid (points e.g. every x,y=0.1 from -1 to +1)
	;										//1=YZ, 2=XZ, 3=XY
	SetIOption(127,0);  //show sensors

	SetIOption(128,1);  //bodies transparent
	SetIOption(129,1);  //constraints transparent
	SetIOption(130,1);  //sensors transparent
	SetIOption(131,1);  //background objects transparent ****
	SetIOption(132,1);  //draw bodies supersmooth
	SetIOption(133,1);  //draw bodies outline
	SetIOption(134,1);  //draw bodies faces
	SetIOption(135,1);  //draw constraints outline **** --> faces is IOption 114

	SetIOption(136,16); //axis tiling (for element face and outline, beams and plates)
	SetIOption(137,8);  //axis resolution: contour plot resolution along axis, beams and plates
	SetIOption(138,4);  //cross-section resolution: contour plot resolution at cross-section, beams and plates
	SetIOption(139,2);  //contour plot resolution for solid finite elements

	SetIOption(140,0);	//rotation input mode: 0=Euler angles, 1=Rotation X/Y/Z, 2=Euler parameters
	SetIOption(141,1);  //show startup banner

	SetIOption(142,0);	//draw loads
	SetIOption(143,0);  //Units in Legend: 0=SI(m,N, etc.); 1=mm, N, etc.
	SetIOption(144,1);	//draw control elements

	SetIOption(145,1);  //draw nodes
	SetIOption(146,1);  //draw surface elements only
	SetIOption(147,3);  //draw node resolution
	SetIOption(148,0);  //show node numbers

	SetIOption(149,0);  //use cutting plane

	SetIOption(150,0);  //animate deformation scaling
	SetIOption(151,0);  //scale deformation in rigid bodies


	//+++++++++++++++++++++++++++
	//double:
	SetDOption(100,0.); //maxstress
	SetDOption(101,0.); //minstress
	SetDOption(102,1.); //finite element line thickness (outline of 2D and 3D beam, plate)
	SetDOption(103,0.5);//origin size
	SetDOption(104,0.); //body local frame size ****

	SetDOption(105,1.); //Deformation scale factor
	SetDOption(106,1.); //shrinking factor

	SetDOption(107,-1.);//grid ref pos x
	SetDOption(108,-1.);//grid ref pos y
	SetDOption(109, 0.);//grid ref pos z
	SetDOption(110,0.1);//grid step1
	SetDOption(111,0.1);//grid step2
	SetDOption(112,2);	//grid size1
	SetDOption(113,2);	//grid size2
	SetDOption(114,1.); //global_line_thickness (coord system, etc.) ****
	SetDOption(115,1.); //rigid body (outline) line thickness ****
	SetDOption(116,1.); //constraint (outline) line thickness ****
	SetDOption(117,2.); //global_point_size (coord system, grid, etc.) ****
	SetDOption(118,0.2);//sensor cross size
	SetDOption(119,2.); //GeomElement (outline) line thickness ****

	SetDOption(120,0.1); //size of arrow for drawing of loads
	SetDOption(121,0.6); //R-value for drawing of loads
	SetDOption(122,0.6); //G-value for drawing of loads
	SetDOption(123,0. ); //B-value for drawing of loads

	SetDOption(124,0.001); //Draw node size

	SetDOption(125,1.); //Cutting plane normal-X
	SetDOption(126,0.); //Cutting plane normal-Y
	SetDOption(127,0.); //Cutting plane normal-Z
	SetDOption(128,0.); //Cutting plane distance

	//Graphics && OPENgl:
	SetIOption(200,1);  //rot axis for standard view angle_1 (rot axis 1, 2 or 3)
	SetIOption(201,2);  //rot axis for standard view angle_2 (rot axis 1, 2 or 3)
	SetIOption(202,3);  //rot axis for standard view angle_3 (rot axis 1, 2 or 3)
	SetIOption(203,1);  //rot axis for standard view angle2_1 (rot axis 1, 2 or 3)
	SetIOption(204,2);  //rot axis for standard view angle2_2 (rot axis 1, 2 or 3)
	SetIOption(205,3);  //rot axis for standard view angle2_3 (rot axis 1, 2 or 3)
	SetIOption(206,1);  //OpenGL lighting
	SetIOption(207,1);  //OpenGL SMOOTH ShadeModel smooth
	SetIOption(208,1);  //OpenGL enable light1
	SetIOption(209,1);  //OpenGL enable light2
	SetIOption(210,0);  //OpenGL light1 mode (0=standard, 1=use light position)
	SetIOption(211,0);  //OpenGL light2 mode (0=standard, 1=use light position)
	SetIOption(212,1);  //Immediate apply in openGL dialog
	//SetIOption(213,0);  //OpenGL 

	//double:
	SetDOption(200,0.);    //standard view angle_1
	SetDOption(201,0.);    //standard view angle_2
	SetDOption(202,0.);    //standard view angle_3
	SetDOption(203,0.);    //standard view angle2_1
	SetDOption(204,0.);    //standard view angle2_2
	SetDOption(205,0.);    //standard view angle2_3
	SetDOption(206,1.);    //global transparency for SetColor, 1=no translucency, 0=fully transparent
	SetDOption(207,0.25);  //light1 ambient parameter
	SetDOption(208,0.40);  //light1 diffuse parameter
	SetDOption(209,0.40);  //light1 specular parameter
	SetDOption(210,0.25);  //light2 ambient parameter
	SetDOption(211,0.40);  //light2 diffuse parameter
	SetDOption(212,0.00);  //light2 specular parameter
	SetDOption(213, 1.);   //light1 posx
	SetDOption(214, 1.);   //light1 posy
	SetDOption(215,-1.);   //light1 posz
	SetDOption(216, 0.);   //light2 posx
	SetDOption(217, 3.);   //light2 posy
	SetDOption(218, 2.);   //light2 posz
	SetDOption(219,60.);   //material shininess (0..128)
	SetDOption(220, 1.);   //material specular color intensity



	//+++++++++++++++++++++ to be made:
	//general tiling factor for GeomObjects, rigid bodies: cylinders, spheres, etc.
	//general tiling factor for Joints: cylinders, spheres


	//250..299:

	//int
	// draw options
	SetIOption(251, 1);  // draw_COG for Rigid3DMinCoord  (default: 1) 
	SetIOption(258, 12); // draw resolution for Rigid3DMinCoord
	//double
	SetDOption(259, 1);  // cog_factor  for Rigid3DMinCoord (default: 1)
	//-------------------------------
	// old:
	//SetIOption(250, 1);// draw_body  (default: 1)  (cmp. Rigid3DMinCoord)
	//SetIOption(252, 1);  // drawCOORD  (default: 1)
	//SetIOption(253, 1);  // drawJOINT  (default: 1)  

	//SetIOption(255, 1);  // show_text  (default: 1)
	//SetIOption(256, 0);  // startMCUManually (default: 0)
	//SetIOption(257, 2);  // controllerModus, describes, how the force/torque is computed w.r.t. the reference trajectory
	////double: 
	//SetDOption(250, 0.36);   // currentToTorque (Nm/Arms), KT (root mean square)  
	//SetDOption(251, 0.955);  // TP (Nm) peak torque                               
	//SetDOption(252,	0.318);  // TS(Nm) continuous stall torque                    
	//SetDOption(253, 3.6);    // IP(Arms) peak current                             
	//SetDOption(254,	1.0);    // IS (Arms) continuous stall current                
	//SetDOption(255,	0.00000233); //		driveInertia from data sheet
	// identified Factors
	//SetDOption(256,	4.5770e-4); //		drive1.currentFactor = MotorGain/currentToTorque; % 4.7222e-005
	//SetDOption(257, MY_PI/180.0);  //angleConversionFactor =  MY_PI/180.0;  // °  -> rad
	//SetDOption(258, 1.0e-3);  //lengthConversionFactor = 1.0e-3;       // mm -> m 		
	//
	//SetDOption(260,-1.0);  // drawing position controller offset x
	//SetDOption(261, 0.0);  // drawing position controller offset y
	//text:
	//SetTOption(250, "..\\..\\KB\\MotionControl"); //staRepPath
	//SetTOption(251, "02_Garbot\\B1WithDynamics.log"); //staRepFileName
	//SetTOption(252, "..\\..\\KB\\MotionControl"); //contrPath
	//SetTOption(253, "02_Garbot\\iosystem.cfg"); //contrFileName
	//---------------------------------------------------------------
*/
}

void NumNLSys::SetIOption(int index, int data)
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal ioption called, index=" << index << "\n";
#endif
	ioptions[index] = data;
}

const int& NumNLSys::GetIOption(int index) const
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal ioption called, index=" << index << "\n";
#endif
	return ioptions[index];	
}

int& NumNLSys::GetIOption(int index)
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal ioption called, index=" << index << "\n";
#endif
	return ioptions[index];	
}

void NumNLSys::SetDOption(int index, double data)
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal doption called, index=" << index << "\n";
#endif
	doptions[index] = data;
}

const double& NumNLSys::GetDOption(int index) const
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal doption called, index=" << index << "\n";
#endif
	return doptions[index];	
}

double& NumNLSys::GetDOption(int index)
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal doption called, index=" << index << "\n";
#endif
	return doptions[index];	
}

void NumNLSys::SetTOption(int index, const char* data)
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal toption called, index=" << index << "\n";
#endif
	if (toptions[index] != 0) delete[] toptions[index];

	int length = strlen(data);
	toptions[index] = new char[length + 1];
	strcpy(toptions[index], data);

	//toptions[index] = data;
}

const char* NumNLSys::GetTOption(int index) const
{
#ifdef _DEBUG
	if (index < min_value_legal_option) 
		UO() << "Warning: illegal toption called, index=" << index << "\n";
#endif
	return toptions[index];	
}


