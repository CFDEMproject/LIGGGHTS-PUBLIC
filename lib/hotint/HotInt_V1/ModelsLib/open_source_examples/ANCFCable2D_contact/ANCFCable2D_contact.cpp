//#**************************************************************
//#
//# filename:             ANCFCable2D_contact.cpp
//#
//# author:               Schörgenhumer Markus
//#
//# generated:						26.September 2013
//# description:          Test model of a flexible ANCF beam and a circular rigid body
//#												moving within a rigid frame in the gravity field, including
//#												mutual mechanical contact.
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


int Generate_InitModelData_ANCFCable2D_contact(MBS* mbs)
{
	return mbs->ReadModelData("../../ModelsLib/open_source_examples/ANCFCable2D_contact/ANCFCable2D_contact.txt");
}


int Generate_Model_ANCFCable2D_contact(MBS* mbs)
{
	ElementDataContainer* edc = mbs->GetModelDataContainer();

	int nel = edc->TreeGetInt("Geometry.n_fibers");
	double sx = edc->TreeGetDouble("Geometry.length");
	double sy = edc->TreeGetDouble("Geometry.width");
	int nx = edc->TreeGetInt("Geometry.nx");
	int ny = edc->TreeGetInt("Geometry.ny");

	double rho = edc->TreeGetDouble("Geometry.rho");
	double Em = edc->TreeGetDouble("Geometry.Em");
	double nu = edc->TreeGetDouble("Geometry.nu");

	double box_x = edc->TreeGetDouble("Geometry.box_x");
	double box_y = edc->TreeGetDouble("Geometry.box_y");
	int nbox_x = edc->TreeGetInt("Geometry.nres_x");
	int nbox_y = edc->TreeGetInt("Geometry.nres_y");

	Vector3D size(sx,sy,1.0);
	double cdim = sy/2;
	double wi = 1;	//width of GeomLine2D elements (in pts/pixel)

	//===============================2D fibers=========================================

	ANCFCable2D cable(mbs);
	Vector xc1(4);
	Vector xc2(4);
	double phi = -MY_PI/4.;
	xc1(1)=0.5*sx*cos(phi+MY_PI);
	xc1(2)=0.5*sx*sin(phi+MY_PI);
	xc1(3)=cos(phi);
	xc1(4)=sin(phi);
	xc2(1)=xc1(1)+sx*cos(phi);
	xc2(2)=xc1(2)+sx*sin(phi);
	xc2(3)=cos(phi);
	xc2(4)=sin(phi);

	//Material m1(mbs,rho,Em,nu);
	//int mat1 = mbs->AddMaterial(&m1);

	cable.SetANCFCable2D(xc1, xc2, rho, Em, size, Vector3D(0.,0.7,0.));
	//cable.SetANCFCable2D(xc1, xc2, vcenter, vcenter, n1, n2, rho, Em, size, Vector3D(0.,0.7,0.));
	int nr = mbs->AddElement(&cable);

	MBSLoad grav;
	grav.SetBodyLoad(-9.81*rho,2);
	mbs->GetElement(nr).AddLoad(grav);

	//MBSSensor force_x(mbs,TMBSSensor(TSElement+TSDOF),idx1,1); //measure force via Lagrange multiplier
	//force_x.SetSensorName(mystr("Node_")+mystr(i)+mystr("_force_x"));
	//mbs->AddSensor(&force_x);

	TArray<Vector2D> points;
	double dx = sx/nx;
	double dy = sy/ny;

	for(int i=0; i<nx; ++i)
		points.Add(Vector2D(-0.5*sx+i*dx,-0.5*sy));
	for(int i=0; i<ny; ++i)
		points.Add(Vector2D(0.5*sx,-0.5*sy+i*dy));
	for(int i=0; i<nx; ++i)
		points.Add(Vector2D(0.5*sx-i*dx,0.5*sy));
	for(int i=0; i<=ny; ++i)
		points.Add(Vector2D(-0.5*sx,0.5*sy-i*dy));

	//mbs->GetElement(nr).SetAltShape(1);

	////sensors for nodal positions and velocities
	//MBSSensor s1(mbs,TMBSSensor(TSElement+TSplanar+TSPos+TSX),nr,Vector3D(-0.5*size.X(),0.,0.));
	//s1.SetSensorName(mystr("Node_")+mystr(i)+mystr("_x"));
	//mbs->AddSensor(&s1);

	//sensors
	//field variables not available?
	//{
	//	FieldVariableElementSensor s1(mbs);
	//	s1.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_position,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x),Vector2D(0.));
	//	s1.SetSensorName(mystr("cable")+mystr("_x"));
	//	mbs->AddSensor(&s1);

	//	FieldVariableElementSensor s2(mbs);
	//	s2.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_position,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y),Vector2D(0.));
	//	s2.SetSensorName(mystr("cable")+mystr("_y"));
	//	mbs->AddSensor(&s2);

	//	FieldVariableElementSensor s3(mbs);
	//	s3.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_velocity,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x),Vector2D(0.));
	//	s3.SetSensorName(mystr("cable")+mystr("_vx"));
	//	mbs->AddSensor(&s3);

	//	FieldVariableElementSensor s4(mbs);
	//	s4.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_velocity,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y),Vector2D(0.));
	//	s4.SetSensorName(mystr("cable")+mystr("_vy"));
	//	mbs->AddSensor(&s4);
	//}

	//contact
	Vector3D contactcol(0.5,0,0.5);
	int slaveNODEmode = 0; //if NODEmode = 1, then use locnodenumbers, if NODEmode==0 then use loccoords
	double bordersize = 0.25*sy; //additional search radius for master and slave segments/nodes
	GeneralContact2D gc(mbs, slaveNODEmode, bordersize, Vector3D(0.0005,0,0), contactcol);
	gc.SetContactMode(0); //0 for Hertzian contact with restitution coefficient
	gc.SetIsLagrange(0);
	double friccoeff = 0.2;
	gc.SetFriction(1, edc->TreeGetDouble("Geometry.friction_coeff"));
	gc.SetContactParams(edc->TreeGetDouble("Geometry.restitution_coeff"),1); //coefficient of restitution, Hertzian contact parameter
	gc.SetContactMaxDist(0.5*sy); //max penetration; if exceeded, it is treated as if there where no contact
	gc.SetSearchTreeDim(20,20);
	double cstiff = edc->TreeGetDouble("Geometry.contact_stiffness");
	int bodyind = 1;

	for(int i=1; i<points.Length(); ++i)
	{
		gc.AddSlaveNode(nr, points(i), cstiff, bodyind);
	}
	
	if(edc->TreeGetInt("Geometry.mutual_contact"))
	{
		for(int i=1; i<points.Length(); ++i)
		{
			gc.AddMasterSegment(nr,points(i),points(i+1),bodyind); //be careful with orientation of master segments
		}
	}

	//===============================rigid body====================================
	{
		Vector x0i(6);
		x0i(1)=0.;
		x0i(2)=0.25*box_y;
		x0i(3)=0.;
		x0i(4)=0.;
		x0i(5)=0.;
		x0i(6)=0.;
		double r0=0.3*sx;
		Vector3D sizei(r0,r0,1.);
		Vector3D coli(1.,0.,0.);
		Rigid2D testbody(mbs,x0i,rho,sizei,coli);
		int nr = mbs->AddElement(&testbody);
		//MBSLoad load;
		//load.SetForceVector2D(Vector2D(1e-4,2e-4),Vector2D(0.));
		mbs->GetElement(nr).AddLoad(grav);

		TArray<Vector2D> points;
		int ni = 32;
		for(int j=0; j<=ni; ++j)
			points.Add(Vector2D( r0*cos(2*MY_PI/ni*j), r0*sin(2*MY_PI/ni*j) ));


		//for better visualization of rotation
		GeomLine2D c1(mbs,nr,Vector2D(0.,-0.5*r0),Vector2D(0.,0.5*r0),Vector3D(1.,0.,0.));
		GeomLine2D c2(mbs,nr,Vector2D(-0.5*r0,0.),Vector2D(0.5*r0,0.),Vector3D(1.,0.,0.));
		c1.SetDrawParam(Vector3D(2*wi, 10., 0.));
		c2.SetDrawParam(Vector3D(2*wi, 10., 0.));
		mbs->GetElement(nr).Add(c1);
		mbs->GetElement(nr).Add(c2);

		mbs->GetElement(nr).SetAltShape(1);
		
		//sensors
		//{
		//	FieldVariableElementSensor s1(mbs);
		//	s1.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_position,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x),Vector2D(0.));
		//	s1.SetSensorName(mystr("rigid")+mystr("_x"));
		//	mbs->AddSensor(&s1);

		//	FieldVariableElementSensor s2(mbs);
		//	s2.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_position,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y),Vector2D(0.));
		//	s2.SetSensorName(mystr("rigid")+mystr("_y"));
		//	mbs->AddSensor(&s2);

		//	FieldVariableElementSensor s3(mbs);
		//	s3.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_velocity,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x),Vector2D(0.));
		//	s3.SetSensorName(mystr("rigid")+mystr("_vx"));
		//	mbs->AddSensor(&s3);

		//	FieldVariableElementSensor s4(mbs);
		//	s4.SetFVESPos2D(nr,FieldVariableDescriptor(FieldVariableDescriptor::FieldVariableType::FVT_velocity,FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y),Vector2D(0.));
		//	s4.SetSensorName(mystr("rigid")+mystr("_vy"));
		//	mbs->AddSensor(&s4);
		//}

		////lock rotation
		//if(edc->TreeGetInt("Geometry.lock_rigid_body_rotation"))
		//{
		//	CoordConstraint cc1(mbs, nr, 3, cdim);
		//	mbs->AddElement(&cc1);
		//}

		//contact
		bodyind = 2;
		for(int i=1; i<=points.Length(); ++i)
		{
			gc.AddSlaveNode(nr, points(i), cstiff, bodyind);
		}

		if(edc->TreeGetInt("Geometry.mutual_contact"))
		{
			for(int i=1; i<points.Length(); ++i)
			{
				gc.AddMasterSegment(nr,points(i),points(i+1),bodyind); //be careful with orientation of master segments
				GeomLine2D c1(mbs,nr,points(i),points(i+1),Vector3D(1.,0.,0.));
				c1.SetDrawParam(Vector3D(2*wi, 10., 0.));
				mbs->GetElement(nr).Add(c1);
			}
		}

	}

	//===============================frame=========================================

	points.Flush();
	dx = box_x/nbox_x;
	dy = box_y/nbox_y;

	for(int i=0; i<nbox_x; ++i)
		points.Add(Vector2D(-0.5*box_x+i*dx,-0.5*box_y));
	for(int i=0; i<nbox_y; ++i)
		points.Add(Vector2D(0.5*box_x,-0.5*box_y+i*dy));
	for(int i=0; i<nbox_x; ++i)
		points.Add(Vector2D(0.5*box_x-i*dx,0.5*box_y));
	for(int i=0; i<=nbox_y; ++i)
		points.Add(Vector2D(-0.5*box_x,0.5*box_y-i*dy));

	//contact
	bodyind = 0; //must be 0 for mbs
	for(int i=1; i<points.Length(); ++i)
	{
		gc.AddMasterSegment(0,points(i+1),points(i), bodyind); //be careful with orientation of master segments
		GeomLine2D c1(mbs,0,points(i+1),points(i),Vector3D(0.,0.,0.));
		c1.SetDrawParam(Vector3D(2*wi, 10., 0.));
		mbs->Add(c1);
	}

	//finish contact
	gc.FinishContactDefinition();
	mbs->AddElement(&gc);

	mbs->Assemble();

	return 1;
};


class Models_AutoRegistration_ANCFCable2D_contact
{
public:
	Models_AutoRegistration_ANCFCable2D_contact()
	{
		ModelFunctionAutoRegistration(&Generate_Model_ANCFCable2D_contact, "ANCFCable2D_contact", "ANCFCable2D_contact", 0, &Generate_InitModelData_ANCFCable2D_contact); 
	}
};
Models_AutoRegistration_ANCFCable2D_contact dummy_ANCFCable2D_contact;