//#***************************************************************************************
//#
//# filename:             flexible_web_contact_3D.cpp
//#
//# author:               Peter Gruber
//#
//# generated:            December 2012
//#
//# description:          Rigid obstacle in contact with a flexible web
//#
//# remarks:              An obstacle (Rigid3D) with initial (angular) velocity and under
//#                       gravity falls into a rectangular (nx x ny)-web made of flexible 
//#                       thin beam elements in 3D (see Fig. 1). At all of the web's nodes
//#                       rigid joints are connecting horizontal and vertical beams (which
//#                       themselves are consisting of beam elements and nodes).
//#                       Moreover, the web's corner nodes are spatially fixed, and nodes
//#                       on the web's edges (except for the corner nodes) are connected
//#                       to ground by spring-dampers.
//#
//#                       +-------+-------+-------+     -
//#                       |       |       |       |     |
//#                       |       |       |       |     |
//#                       |       |       |       |     |
//#                       +-------+-------+-------+     |
//#                       |       |       |       |     |
//#                       |       |   o   |       |     width
//#                       |       |       |       |     |
//#                       +-------+-------+-------+     |        
//#                       |       |       |       |     |        y
//#                       |       |       |       |     |        |
//#                       |       |       |       |     |        |
//#                       +-------+-------+-------+     -        +----x
//#                                                             /
//#                       |------- length --------|            z
//#
//#                       Fig 1.: web made out of connected beams of a certain length and
//#                       width, where 'o' denotes the center position. The refinement is
//#                       the same in x and y direction, i.e., one ends up with a
//#                       discretization of (nx x ny) elements and (nx+1) x (ny+1) nodes.
//#                       More details to be found in a GACM Report (to appear in 2013/14)
//#
//#                       Contact: -----------------------------------------
//#                       The German Association for Computational Mechanics
//#                       Institut für Kontinuumsmechanik
//#                       Leibniz Universität Hannover
//#                       Appelstraße 11, D-30167 Hannover, Germany
//#                       Phone: +49-511-762-3320, Fax: +49-511-762-5496
//#                       email: wriggers@ikm.uni-hannover.de
//#                       Internet: www.gacm.de
//#                       --------------------------------------------------
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


#include "model_tools.h"

// Model BEGIN
// *
// *
int GenerateContactWeb_InitModelData(MBS* mbs)
{
	return mbs->ReadModelData("..\\..\\ModelsLib\\open_source_examples\\flexible_web_contact_3D\\modeldata.txt");
}
// *
// *
int GenerateContactWeb(MBS* mbs)
{
	ElementDataContainer* edc = mbs->GetModelDataContainer();

	Vector3D center_position(0.);

	// geometrical data, see Fig. 1
	int nx = edc->TreeGetInt("ModelOptions.nx",10);  // number of beam elements in x direction
	int ny = edc->TreeGetInt("ModelOptions.ny",10);  // number of beam elements in y direction
	double ly = 0.02;   //thickness of ropes in local y-dir (local x-dir is in rope-direction)
	double lz = 0.02;   //thickness of ropes in local z-dir (local x-dir is in rope-direction)
	double length = 1.5;   //length of ropes which run in global x-dir
	double width = 1.5;    //length of ropes which run in global y-dir
	double dx = length/nx;
	double dy = width/ny;

	// material (hemp fiber)
	Vector3D beam_color(0.4,0.4,0.8);
	Beam3DProperties bp(mbs);
	double rho = 0.25*0.25*1.48e3;           // material density [kg/m^3]
	double obstacle_massdensity = 2700;      // Aluminium
	double E = 0.25*0.25*6.9e9;              // Young's modulus [Pa]
	double nu = 0.3;               // Poisson's ratio (needed just for calculation of torsional stiffness through shear modulus G) [0]
	double G = E/(2.+2.*nu);       // shear modulus [Pa]
	double Iy = pow(lz,3)*ly/12.;  // moment of inertia around local y-dir [m^4]
	double Iz = pow(ly,3)*lz/12.;  // moment of inertia around local z-dir [m^4]
	double Ix = Iy + Iz;           // moment of inertia around local x-dir (for a circular cross section) [m^4]
	double kx = 0.8436;            // correction factor between Ix for circular and \tilde Ix for quadratic cross section (see Schwab paper)
	double A = ly*lz;              // size of cross section [m^2]
	bp.SetBeam3DProperties(1, Vector(0.3*ly,0.3*lz), rho, E*A, E*Iy, E*Iz, G*Ix*kx, 0, 0, rho*A, rho*Ix, rho*Iy, rho*Iz);
	int matnr = mbs->AddMaterial(bp);

	// load
	MBSLoad gravity;
	gravity.SetGravity(-9.81, 3);  // gravitation in global z-dir

	// storage of node and element indices
	TArrayDynamic<IVector> node_idx_x;
	TArrayDynamic<IVector> node_idx_y;
	TArrayDynamic<IVector> element_idx_x;
	TArrayDynamic<IVector> element_idx_y;
	
	// generate ropes in global x-dir
	ANCFNodeS1rot1_3D initial_node_x(mbs);   // initial node of beam in x-direction
	ANCFNodeS1rot1_3D last_node_x(mbs);      // last node of beam in x-direction
	initial_node_x.SetANCFNodeS1rot1_3D(center_position + Vector3D(-.5*length, -.5*width, 0.), Vector3D(0.));
	last_node_x.SetANCFNodeS1rot1_3D(center_position + Vector3D(.5*length, -.5*width, 0.), Vector3D(0.));
	for (int i=0; i<=ny; i++)
	{
		node_idx_x.Add(GenerateNodesOnStraightLine(mbs, initial_node_x, last_node_x, nx+1));   // horizontal ropes (x-dir)

		element_idx_x.Add(GenerateANCFBeam3DTorsionBeam(mbs, node_idx_x.Last(), matnr, beam_color));
		ApplyLoadToElements(mbs, element_idx_x.Last(), gravity);
		
		initial_node_x.Pos().Y() += dy;
		last_node_x.Pos().Y() += dy;
	}

	// generate ropes in global y-dir
	ANCFNodeS1rot1_3D initial_node_y(mbs);   // initial node of beam in y-direction
	ANCFNodeS1rot1_3D last_node_y(mbs);      // last node of beam in y-direction
	initial_node_y.SetANCFNodeS1rot1_3D(center_position + Vector3D(-.5*length, -.5*width, 0.), Vector3D(0.,0.,MY_PI/2));
	last_node_y.SetANCFNodeS1rot1_3D(center_position + Vector3D(-.5*length, .5*width, 0.), Vector3D(0.,0.,MY_PI/2));
	for (int i=0; i<=nx; i++)
	{		
		node_idx_y.Add(GenerateNodesOnStraightLine(mbs, initial_node_y, last_node_y, ny+1));   // vertical ropes (x-dir)

		element_idx_y.Add(GenerateANCFBeam3DTorsionBeam(mbs, node_idx_y.Last(), matnr, beam_color));
		ApplyLoadToElements(mbs, element_idx_y.Last(), gravity);

		initial_node_y.Pos().X() += dx;
		last_node_y.Pos().X() += dx;
	}

	// rigid constraints at cross points
	double constraint_dim = -1;
	Vector3D constraint_col(0.4,0.7,0.3);
	for (int i=1; i<=element_idx_x.Length(); i++)
	{
		for (int j=1; j<=element_idx_y.Length(); j++)
		{
			Vector3D lp1(-dx/2,0,0); 
			Vector3D lp2(-dy/2,0,0);

			int k = i;
			int l = j;
			if (i == element_idx_x.Length())
			{
				k--;
				lp2.X() += dy;
			}
			if (j == element_idx_y.Length())
			{
				l--;
				lp1.X() += dx;
			}

			int en1 = element_idx_x(i)(l);
			int en2 = element_idx_y(j)(k);
			
			RigidJoint rj(mbs);
			rj.SetRigidJoint_LocalPos_to_LocalPos(en1, lp1, en2, lp2);
			rj.SetDrawDim(Vector3D(constraint_dim, 0, 0));
			mbs->AddElement(&rj);
		}
	}


	// constraints at outer corners of the web
	MathFunction mathfunc_minus, mathfunc_plus;
	mathfunc_plus.SetExpression("(1-cos(pi*t/SolverOptions.end_time))/2", "t", mbs);
	mathfunc_minus.SetExpression("(-1+cos(pi*t/SolverOptions.end_time))/2", "t", mbs);

	IVector pcc(0);
	double constraint_dim_save = constraint_dim;
	constraint_dim = lz;

	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1)(1), 1, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1)(1), 2, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1)(1), 3, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1).Last(), 1+6, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1).Last(), 2+6, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x(1).Last(), 3+6, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last()(1), 1, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last()(1), 2, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last()(1), 3, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last().Last(), 1+6, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last().Last(), 2+6, constraint_dim)));
	pcc.Add(mbs->AddElement(&PrescribedCoordConstraint(mbs, element_idx_x.Last().Last(), 3+6, constraint_dim)));

	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(1)))->GetMathFunction() = mathfunc_plus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(2)))->GetMathFunction() = mathfunc_plus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(4)))->GetMathFunction() = mathfunc_minus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(5)))->GetMathFunction() = mathfunc_plus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(7)))->GetMathFunction() = mathfunc_plus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(8)))->GetMathFunction() = mathfunc_minus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(10)))->GetMathFunction() = mathfunc_minus;
	static_cast<PrescribedCoordConstraint*>(mbs->GetElementPtr(pcc(11)))->GetMathFunction() = mathfunc_minus;

	// add spring-dampers at outer edges of the web
	double damping_parameter = 10;
	int elDOF = 14;
	for (int i=2; i<=element_idx_x.Length()-1; i++)
	{
		int elnr;
		SpringDamperActuator sda1(mbs);
		elnr = element_idx_x(i)(1);
		sda1.SetSpringDamperActuator_LocalPos_to_LocalPos(elnr, Vector3D(-dx/2,0,0), 0, mbs->GetElement(elnr).GetNode(1).Pos() + Vector3D(-length/20,0,0), 0, damping_parameter);
		mbs->AddElement(&sda1);

		SpringDamperActuator sda2(mbs);
		elnr = element_idx_x(i).Last();
		sda2.SetSpringDamperActuator_LocalPos_to_LocalPos(elnr, Vector3D(dx/2,0,0), 0, mbs->GetElement(elnr).GetNode(2).Pos() + Vector3D(length/20,0,0), 0, damping_parameter);
		mbs->AddElement(&sda2);
	}

	for (int i=2; i<=element_idx_y.Length()-1; i++)
	{
		int elnr;
		SpringDamperActuator sda1(mbs);
		elnr = element_idx_y(i)(1);
		sda1.SetSpringDamperActuator_LocalPos_to_LocalPos(elnr, Vector3D(-dy/2,0,0), 0, mbs->GetElement(elnr).GetNode(1).Pos() + Vector3D(0,-width/20,0), 0, damping_parameter);
		mbs->AddElement(&sda1);

		SpringDamperActuator sda2(mbs);
		elnr = element_idx_y(i).Last();
		sda2.SetSpringDamperActuator_LocalPos_to_LocalPos(elnr, Vector3D(dy/2,0,0), 0, mbs->GetElement(elnr).GetNode(2).Pos() + Vector3D(0,width/20,0), 0, damping_parameter);
		mbs->AddElement(&sda2);
	}

	constraint_dim = constraint_dim_save;

	
	

	// add contact	
	double contact_stiffness = 20e4;     // contact stiffness - the smaller this positive number is, the easier the numerical solver converges,
	                                     // determine the contact stiffness in particular dependency of the awaited maximum impulse!
	double radius = 0.4*lz;              // contact radius (used for penalty formulation - i.e., soft contact)
	int n_slave_points_per_beam = length/((double)nx*radius)*2; //50;
	mbs->UO() << "n_slave_points_per_beam = " << n_slave_points_per_beam << "\n";
	int slaveNODEmode = 0;		           // if slaveNODEmode==1 then use locnodenumbers, if slaveNODEmode==0 then use loccoords
	GeneralContact3D contact(mbs, slaveNODEmode, 2*radius, constraint_dim, constraint_col);
	contact.SetContactMode(3);           // 0 = with restitution, 3 = lin spring, 5=correct restitution
	contact.SetIsLagrange(0);            // Penalty instead of Lagrange formalism (Lagrange not implemented for contact radius > 0)
	contact.SetFriction(0, 0);           // no friction considered in this model
	contact.SetContactParams(0.5,1);     // restitution coefficient and hertzian contact parameter
	contact.SetContactMaxDist(2*radius); // if the nearest distance between the pair of a slave node and a master is larger than
	                                     // this value, then this pair is not considered for contact computation right away.
	contact.SetSearchTreeDim(40,40,12);  // (x,y,z)-subdivision of the bounding box for building search tree

	for (int k=1; k<=n_slave_points_per_beam; k++)
	{
		double xi = ((double)k - .5)/n_slave_points_per_beam - .5;   // distributed in [-1/2, 1/2]
		Vector3D ploc(xi,0,0);

		for (int i=1; i<=element_idx_x.Length(); i++)
		{
			for (int j=1; j<=element_idx_x(i).Length(); j++)
			{
				int elem_idx = element_idx_x(i)(j);
				contact.AddSlaveNode(elem_idx, ploc*mbs->GetElementPtr(elem_idx)->GetSize().X(), contact_stiffness, i, radius);
			}
		}
		for (int i=1; i<=element_idx_y.Length(); i++)
		{
			for (int j=1; j<=element_idx_y(i).Length(); j++)
			{
				int elem_idx = element_idx_y(i)(j);
				contact.AddSlaveNode(elem_idx, ploc*mbs->GetElementPtr(elem_idx)->GetSize().X(), contact_stiffness, element_idx_x.Length()+i, radius);
			}
		}
	}


	int idx_obstacle = 0;   // element index of obstacle (Rigid3D with alternative geometry)
	if(1) 
		// create obstacle and add trigs of obstacle to
		//   a) masterelements of contact, and
		//   b) the alternative geometry of a Rigid3D element
	{

		mystr workdir("..\\..\\ModelsLib\\open_source_examples\\flexible_web_contact_3D\\");
		mystr filename("obstacle");

		if(0) // create geometry
		{
			ofstream ofs((workdir+filename+mystr(".geo")).c_str());
			ofs.precision(4);

			// file header:
			ofs << "# file automatically generated geometry by HOTINT\n";
			ofs << "# obstacle\n";
			ofs << "#\n";
			ofs << "algebraic3d\n";

			double ol=1;
			double oh=.5;
			double ob=.5;
			double od=.05;
			double ox=.2;
			double oy=ox;
			double oz=2*ox;
			double or=ox/2;

			if (1)
				// create geo file (NETGEN standard geometry file type)
			{
				Vector3D box_center(od/2, ob/2, oh/2);
				Vector3D box_size(od, ob, oh);

				GEOCube mycube1;
				mycube1.SetOrthoCube(box_center, box_size);
				mycube1.Export(ofs, "mycube1");

				box_center.Set(ol/2, ob/2, od/2);
				box_size.Set(ol, ob, od);

				GEOCube mycube2;
				mycube2.SetOrthoCube(box_center, box_size);
				mycube2.Export(ofs, "mycube2");

				Vector3D p1(-od/2, ob/2, oh-oy);
				Vector3D p2(3*od/2, ob/2, oh-oy);

				GEOCyl mycyl1;
				mycyl1.SetCyl(p1, p2, or);
				mycyl1.Export(ofs, "mycyl1");

				p1.Set(ol-ox, ob/2, -oz);
				p2.Set(ol-ox, ob/2, oz+od);

				GEOCyl mycyl2;
				mycyl2.SetCyl(p1, p2, or);
				mycyl2.Export(ofs, "mycyl2");

				ofs << "solid myobstacle = ";
				ofs << "mycube1 and not mycyl1 or mycube2 or mycyl2;\n";
				ofs << "tlo myobstacle -col=[0.6,0.6,0.4];\n";
			}

			// close geo file
			ofs.close();
		}

		if (0)
			// create stl file from geo file with NETGEN
		{
			mystr netgen_path = "C:\\Program Files (x86)\\Netgen-4.9.10_Win32\\bin\\";
			CreateMeshWithNETGEN(filename+mystr(".geo"), filename+mystr(".stl"), netgen_path, workdir, "STL Format", "verycoarse");
		}

		if (1) 
			// a) create geommesh3D according to stl file,
			// b) translate/rotate/scale geommesh,
			// c) create rigid3D from mesh with corresponding density, initial position, initial rotation, initial velocity, and initial angular velocity;
			// d) create trigs from mesh;
			// e) add trigs to rigid3D as alternative geometry (used for drawing);
			// f) add trigs to contact as master elements;
		{
			GeomMesh3D mesh(mbs, 0, colblue);

			// read stl mesh for calculating center of mass (com) and mass moment of inertia (mmi) of rigid3dkardan
			//mesh.ReadSTLMesh(workdir+filename+mystr("_fine.stl"));
			mesh.ReadSTLMesh(workdir+filename+mystr(".stl"));

			// calculate volume, mass moment of inertia of mesh for given density, and then add rigid3d to mbs

			Vector3D obstacle_center_of_mass = mesh.ComputeCenterOfMass();
			mesh.Translate(-1.*obstacle_center_of_mass);
			mesh.Stretch(0.35);/*mesh.Stretch(0.5);*/

			Matrix3D obstacle_massmomentinertia = mesh.ComputeMassMomentOfInertia(obstacle_massdensity);
			double obstacle_volume = mesh.ComputeVolume();

			double phi = 0.125*MY_PI;
			double obstacle_initial_z_displacemenet = radius+3*ly+0.3;
			double turbation = 0.35*MY_PI/2;
			Vector3D obstacle_size(pow(obstacle_volume,1./3.));
			Vector obstacle_xg(6);       // 3*position, 3*velocity
			Vector obstacle_phi(6);      // 3*rotation, 3*angular vel
			obstacle_xg(1) = 0; //length/2.5*(cos(phi)*obstacle_center_of_mass(1)-sin(phi)*obstacle_center_of_mass(2));
			obstacle_xg(2) = 0; //width/2.5*(sin(phi)*obstacle_center_of_mass(1)+cos(phi)*obstacle_center_of_mass(2));
			obstacle_xg(3) = obstacle_center_of_mass(3) + obstacle_initial_z_displacemenet;
			obstacle_xg(6) = -2; //initial velocity
			obstacle_phi(1) = turbation;
			obstacle_phi(2) = MY_PI/2.+turbation;
			obstacle_phi(3) = 0*phi+turbation;
			obstacle_phi(4) = -turbation;
			obstacle_phi(5) = turbation/5;
			obstacle_phi(6) = -turbation/3;
			Vector3D obstacle_color(0.4, 0.4, 0.8);

			Rigid3D obstacle(mbs, obstacle_xg, obstacle_phi, obstacle_massdensity, obstacle_volume, obstacle_massmomentinertia, obstacle_size, obstacle_color);
			obstacle.AddLoad(gravity);
			obstacle.SetAltShape(1);
			idx_obstacle = mbs->AddElement(&obstacle);


			// get trigs from mesh, and add them 1) as GeomTrig3D to the rigid3d element, and 2) as masterelements to the constraint
			for(int i=1; i <= mesh.NTrigs(); i++)
			{
				Vector3D p1, p2, p3, n;
				mesh.GetTrig0(i, p1, p2, p3, n);
				GeomTrig3D trig(mbs, idx_obstacle, p1, p2, p3, Vector3D(0.6,0.6,0.6));
				mbs->GetElementPtr(idx_obstacle)->Add(trig);
				contact.AddMasterGeomElement(&trig, 0);
			}
		}
	}

	contact.FinishContactDefinition();
	int idx_contact = mbs->AddElement(&contact);


	// create sensors
	FieldVariableElementSensor sensorFV(mbs);

	// displacement sensors at center of web
	sensorFV.SetFVESPos3D(element_idx_y(1+nx/2)(1+ny/2), FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_x), Vector3D(-dy/2,0,0));
	sensorFV.SetSensorName(mystr("x-Pos web center"));
	mbs->AddSensor(&sensorFV);
	sensorFV.SetFVESPos3D(element_idx_y(1+nx/2)(1+ny/2), FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_y), Vector3D(-dy/2,0,0));
	sensorFV.SetSensorName(mystr("y-Pos web center"));
	mbs->AddSensor(&sensorFV);
	sensorFV.SetFVESPos3D(element_idx_y(1+nx/2)(1+ny/2), FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_z), Vector3D(-dy/2,0,0));
	sensorFV.SetSensorName(mystr("z-Pos web center"));
	mbs->AddSensor(&sensorFV);

	// displacement sensors at center of obstacle
	if (idx_obstacle)
	{
		sensorFV.SetFVESPos3D(idx_obstacle, FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_x), Vector3D(0,0,0));
		sensorFV.SetSensorName(mystr("x-Pos obstacle center of mass"));
		mbs->AddSensor(&sensorFV);
		sensorFV.SetFVESPos3D(idx_obstacle, FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_y), Vector3D(0,0,0));
		sensorFV.SetSensorName(mystr("y-Pos obstacle center of mass"));
		mbs->AddSensor(&sensorFV);
		sensorFV.SetFVESPos3D(idx_obstacle, FieldVariableDescriptor(FieldVariableDescriptor::FVT_position, FieldVariableDescriptor::FVCI_z), Vector3D(0,0,0));
		sensorFV.SetSensorName(mystr("z-Pos obstacle center of mass"));
		mbs->AddSensor(&sensorFV);
	}

	SingleElementDataSensor sensorSED(mbs);
	for (int i=1; i<=pcc.Length(); i++)
	{
		char dirchar = (i%3==1) ? 'x' : ((i%3==2) ? 'y' : 'z');
		sensorSED.SetSingleElementDataSensor(pcc(i), "Connector.constraint_force");
		sensorSED.SetSensorName(mystr("corner ")+mystr(i/3)+mystr(" force ")+mystr(dirchar));
		mbs->AddSensor(&sensorSED);
	}
	
	// assemble multibody system
	mbs->Assemble();

	// return error code (1 in case of no error)
	return 1;
}
// *
// *
// Model END




// Auto registration BEGIN
// *
// *
// Auto registration class
// (creating an instance of this class calls the default constructor,
//  and default constructor adds pointer to model function (GenerateContactWeb)
//  to global models list)
class models_ContactWeb_AutoRegistration
{
public:
	models_ContactWeb_AutoRegistration()
	{
		int (*pt2Function)(MBS* mbs);
		int (*pt2InitFunction)(MBS* mbs);

		ModelFunctionAutoRegistration(&GenerateContactWeb, "Flexible Web Contact 3D", "Contact between rigid obstacle and flexible web", 0, &GenerateContactWeb_InitModelData); 
	}
};
// *
// *
//this is the global variable, which is built at program start-up
//the constructor is called and the model functions are registered!
models_ContactWeb_AutoRegistration dummy_models_ContactWeb_AutoRegistration;
// *
// *
// Auto registration END


