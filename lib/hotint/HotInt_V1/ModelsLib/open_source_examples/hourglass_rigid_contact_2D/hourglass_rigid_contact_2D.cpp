//#**************************************************************
//#
//# filename:             hourglass_rigid_contact_2D.cpp
//#
//# author:               Gerstmayr Johannes, Schörgenhumer Markus
//#
//# generated:						26.September 2013
//# description:          A set of randomly-sized rigid particles with two different
//#												densities and circular or random polygonal shapes is falling in an
//#												hourglass-shaped domain driven by gravity. Additionally, in the lower
//#												part of the domain, a rigid-body plate is rotating with a given 
//#												initial angular velocity. Using a contact formulation, all 
//#												interactions between the particles, the rotating body, and the rigid 
//#												walls are captured by the simulation.
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


int Generate_InitModelData_hourglass_rigid_contact_2D(MBS* mbs)
{
	return mbs->ReadModelData("../../ModelsLib/open_source_examples/hourglass_rigid_contact_2D/hourglass_rigid_contact_2D.txt");
}


int Generate_Model_hourglass_rigid_contact_2D(MBS* mbs)
{
	ElementDataContainer* edc = mbs->GetModelDataContainer();

	mbs->GetSolSet().withgraphics = 5;

	TArray<int> elnums(100);
	TArray<int> locelnums(100);
	TArray<int> elnums_top(100);
	TArray<int> elnums_bottom(100);
	TArray<int> elnums_right(100);
	TArray<int> elnums_left(100);

	TArray<Vector2D> boundarynode(100);
	TArray<Vector2D> contactnode(100);

	Vector2D xb[9];
	Vector2D vb[9];

	TArray<Vector2D> sdapos(100);
	TArray<int> sdaelnum(100);


	double f = 1.2;  //large:0.4
	double lz = 0.005*f; 
	double h = 0.05*f; //cube size
	double b = 0.05*f;
	double dist = 0.005*f+0*(0.005-0.004)*f;
	double disty = dist - 0*(0.005*1.4-0.001)*f;

	double hgh = 1.2; //hourglass height
	double hgw = 1.2; //hourglass width
	double hgw2 = 0.75*hgw; //hourglass width
	double hgdist = 0.2*f; //hourglass hole, 0.2 for quads ...


	int nfact = 2;
	int n2 = 8;//45/8; //large: 24 
	int n  = (int)(36./f)/1;	 //large: 44
	int iscontact = 1;
	int iscircle = 1;
	int isavalanche = 0;
	double circleradius = 0.5*h;


	//standard:
	double damp = 1;
	//double rho = 50000;
	double rho = 1000;
	double Em = 1e9*100;
	double nu = 0.2;
	double g = 9.81;

	double friccoeff = 0.2;
	double cstiff = 1e7;

	Vector3D colb(0.2,0.2,0.7);
	Vector3D colb2(0.7,0.2,0.2);
	Vector3D colc(0.3,0.8,0.3); 
	Vector3D contactcol(0.2,0.2,0.2); 
	Vector3D colground(0.5,0.5,0.1); 


	MBSLoad load1;


	int slavenodemode = 0; //no nodes ...
	GeneralContact2D gc(mbs, slavenodemode, 0.01, Vector3D(0.0025,0,0), contactcol);
	gc.SetContactMode(0);
	gc.SetIsLagrange(0);
	gc.SetFriction(1*0, friccoeff);
	gc.SetContactParams(0.5,1);
	gc.SetContactMaxDist(h*0.2);
	gc.SetSearchTreeDim((int)(250*0.2),(int)(400*0.2));


	srand(1); //initialize always with same number->always same model!!!

	int nbody = n*n2; //hourglass
	if (iscontact) 
	{
		int kk = 0;
		double cline = 1.*lz;
		double res = 4;
		TArray<Vector2D> p;

		for (int i = 1; i <= n; i++)
		{
			for (int j = 1; j <= n2; j++)
			{
				kk++;
				double offx = 0.1+0.2+b+b*0.5+0*0.5;
				double offy = -0.25          -0*2.3;

				if (isavalanche) offy += -0.75 + 0.5*(double)i*(b+dist);

				if (j % 2 == 0 && !isavalanche) 
				{
					offx += 0.5*(b+dist);
				}

				Vector xb1(6);
				//Vector3D size(b,h,lz);
				Vector3D size(b,h,1.0);
				xb1(1) = offx + (double)(i-1)*(b+dist);
				xb1(2) = offy + (double)(j-1)*(h+disty);

				double m = 0.125*h*0.5;
				Vector3D size2(size.X()-1*m,size.Y()-1*m,size.Z());
				Vector2D lp[18];
				Vector2D lp1[4];
				lp1[0] = Vector2D( 0.5*b, 0.5*h);
				lp1[1] = Vector2D(-0.5*b, 0.5*h);
				lp1[2] = Vector2D(-0.5*b,-0.5*h);
				lp1[3] = Vector2D( 0.5*b,-0.5*h);


				//additional points:
				lp[ 0] = Vector2D( 0.5 *b, 0.5 *h);
				lp[ 1] = Vector2D( 0.25*b, 0.5 *h);
				lp[ 2] = Vector2D( 0.  *b, 0.5 *h);
				lp[ 3] = Vector2D(-0.25*b, 0.5 *h);

				lp[ 4] = Vector2D(-0.5 *b, 0.5 *h);
				lp[ 5] = Vector2D(-0.5 *b, 0.25*h);
				lp[ 6] = Vector2D(-0.5 *b, 0.  *h);
				lp[ 7] = Vector2D(-0.5 *b,-0.25 *h);

				lp[ 8] = Vector2D(-0.5 *b,-0.5 *h);
				lp[ 9] = Vector2D(-0.25*b,-0.5 *h);
				lp[10] = Vector2D( 0.  *b,-0.5 *h);
				lp[11] = Vector2D( 0.25*b,-0.5 *h);

				lp[12] = Vector2D( 0.5 *b,-0.5 *h);
				lp[13] = Vector2D( 0.5 *b,-0.25*h);
				lp[14] = Vector2D( 0.5 *b, 0.  *h);
				lp[15] = Vector2D( 0.5 *b, 0.25*h);

				if (1)
				{
					//random shapes
					for (int jj = 1; jj <= 4; jj++)
					{
						double rnd = (double)rand()/(double)RAND_MAX;
						//lp[(jj-1)*4] *= 0.1+rnd*0.9;
						lp[(jj-1)*4] *= 0.5+rnd*0.5;
					}
					lp[16] = lp[0];

					for (int jj = 1; jj <= 4; jj++)
					{
						for (int jj2 = 1; jj2 <= 3; jj2++)
						{
							double f = jj2;
							lp[(jj-1)*4+jj2] = f/4.*lp[(jj-1)*4+4]+(1.-f/4.)*lp[(jj-1)*4];
						}
					}
					lp[17] = lp[1];
					for (int jj = 1; jj <= 4; jj++)
					{
						lp[(jj)*4] = 1./4.*(lp[(jj)*4-1]+2.*lp[(jj)*4]+lp[(jj)*4+1]);
					}
					lp[0] = lp[16];

				} 
				else
				{
					//quad shapes with rounded edges
					lp[ 0] += Vector2D(-m,-m);
					lp[ 4] += Vector2D( m,-m);
					lp[ 8] += Vector2D( m, m);
					lp[12] += Vector2D(-m, m);
				}

				Vector3D col = colb;

				double rho2 = rho;
				double circleradius2 = circleradius;
				if (iscircle)
				{
					if (1)
					{
						double rnd = (double)rand()/(double)RAND_MAX;
						//double f = 0.25+0.75*rnd; //original hourglass
						double f = 0.5+0.5*rnd;
						rho2 *= Sqr(f);
						circleradius2 *= f;
					}

					if ((i+j)%2 == 1 && !isavalanche) 
					{
						col = colb2;
						rho2 = 10.*rho2;
					}
				}
				else
				{
					if ((((i-1)/4)%2)==1) 
					{
						col = colb2;
					}
				}

				p.SetLen(16);
				for (int i=1; i <= 16; i++) p(i) = lp[i-1];

				load1.SetBodyLoad(-rho2*g,2);
				Rigid2D b1(mbs, xb1, rho2, size2, col);

				b1.AddLoad(load1);
				b1.SetMassDamping(damp);
				int nr = mbs->AddElement(&b1);

				mbs->GetElement(nr).SetAltShape(1);
				if (!iscircle)
				{
					//random polygons:
					GeomPolygon2D pg(mbs, nr, p, col);
					mbs->GetElement(nr).Add(pg);

					for (int ii=0; ii < 16; ii++)
					{
						gc.AddSlaveNode(nr, lp[ii], cstiff, kk);

						if (m != 0)
						{
							if (ii < 16)
							{
								GeomLine2D ge1(mbs,nr,lp[ii],lp[((ii+1)%16)],contactcol);
								ge1.SetDrawParam(Vector3D(cline,res,0));
								gc.AddMasterGeomElement(&ge1, kk);
							} 
						}
						else
						{
							if (ii < 4)
							{
								GeomLine2D ge1(mbs,nr,lp1[ii],lp1[((ii+1)%4)],contactcol);
								ge1.SetDrawParam(Vector3D(cline,res,0));
								gc.AddMasterGeomElement(&ge1, kk);
							}
						}
					}

				}
				else
				{
					//fixed circles:
					GeomCircle2D ci(mbs, nr, Vector2D(0.,0.), circleradius2, col);
					ci.SetDrawParam(Vector3D(cline,16,0));
					mbs->GetElement(nr).Add(ci);

					gc.AddSlaveNode(nr, Vector2D(0.,0.), cstiff, kk, circleradius2);
					gc.AddMasterGeomElement(&ci, kk);

					TArray<Vector2D> points;
					int cres = 6;
					for(int i=0; i<=cres; ++i)
					{
						points.Add(Vector2D(circleradius2*cos(2.0*MY_PI/cres*i),circleradius2*sin(2.0*MY_PI/cres*i)));
					}

				}

			}
		}


		cline = 1; // *= 4;
		if (!isavalanche)
		{
			//hourglass
			GeomLine2D ge1(mbs,0,Vector2D(hgw,-hgh),Vector2D(0,0),colground);
			ge1.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge1, 0);
			GeomLine2D ge2(mbs,0,Vector2D(hgw*2+hgdist,0),Vector2D(hgw+hgdist,-hgh),colground);
			ge2.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge2, 0);

			GeomLine2D ge3(mbs,0,Vector2D(hgw,-1.3*hgh),Vector2D(hgw,-hgh),colground);
			ge3.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge3, 0);
			GeomLine2D ge4(mbs,0,Vector2D(hgw+hgdist,-hgh),Vector2D(hgw+hgdist,-1.3*hgh),colground);
			ge4.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge4, 0);

			hgw2 *= 1.2;
			double hgh2 = 0.75*hgh;

			//top holder:
			GeomLine2D ge8(mbs,0,Vector2D(hgw-0.5*hgw2,-2*hgh2),Vector2D(hgw,-1.3*hgh),colground);
			ge8.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge8, 0);
			GeomLine2D ge9(mbs,0,Vector2D(hgw+hgdist,-1.3*hgh),Vector2D(hgw+0.5*hgw2+hgdist,-2*hgh2),colground);
			ge9.SetDrawParam(Vector3D(cline,res,0));
			gc.AddMasterGeomElement(&ge9, 0);

			if (1) //driven body:
			{
				Vector2D pc(hgw*1.15, -1.8*hgh);

				Vector xb1(6);
				double lbeam = 0.4;
				double hbeam = 0.025;
				//Vector3D size(lbeam,hbeam,lz);
				Vector3D size(lbeam,hbeam,1.0);
				Vector3D colrigid(0.4,0.4,0.4);
				xb1(1) = pc.X()+lbeam*0.5;
				xb1(2) = pc.Y();
				xb1(6) = 10; 

				Rigid2D b1(mbs, xb1, 1./166.*1e9, size, colrigid);

				int nr = mbs->AddElement(&b1);
				Pos2DConstraint pcc(mbs, nr, Vector2D(-0.5*lbeam,0), pc+Vector2D(-0.5*lbeam,0), 0.001,colc);
				mbs->AddElement(&pcc);


				GeomLine2D ge1(mbs,nr,Vector2D(0.5*lbeam,0.5*hbeam),Vector2D(-0.5*lbeam,0.5*hbeam),contactcol);
				ge1.SetDrawParam(Vector3D(1,1,0));
				gc.AddMasterGeomElement(&ge1, 0);

				GeomLine2D ge2(mbs,nr,Vector2D(-0.5*lbeam,0.5*hbeam),Vector2D(-0.5*lbeam,-0.5*hbeam),contactcol);
				ge2.SetDrawParam(Vector3D(1,1,0));
				gc.AddMasterGeomElement(&ge2, 0);

				GeomLine2D ge3(mbs,nr,Vector2D(-0.5*lbeam,-0.5*hbeam),Vector2D(0.5*lbeam,-0.5*hbeam),contactcol);
				ge3.SetDrawParam(Vector3D(1,1,0));
				gc.AddMasterGeomElement(&ge3, 0);

				GeomLine2D ge4(mbs,nr,Vector2D(0.5*lbeam,-0.5*hbeam),Vector2D(0.5*lbeam,0.5*hbeam),contactcol);
				ge4.SetDrawParam(Vector3D(1,1,0));
				gc.AddMasterGeomElement(&ge4, 0);
			}




			if (0)
			{
				//bottom cup:
				GeomLine2D ge5(mbs,0,Vector2D(hgw-0.5*hgw2,-3*hgh2),Vector2D(hgw-0.5*hgw2,-2*hgh2),colground);
				ge5.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge5, 0);
				GeomLine2D ge6(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-2*hgh2),Vector2D(hgw+0.5*hgw2+hgdist,-3*hgh2),colground);
				ge6.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge6, 0);
				GeomLine2D ge7(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-3*hgh2),Vector2D(hgw-0.5*hgw2,-3*hgh2),colground);
				ge7.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge7, 0);
			}
			else
			{
				//alternative bottom cup:
				double de = 0.2;
				double ww = 0.2;
				GeomLine2D ge5(mbs,0,Vector2D(hgw-0.5*hgw2,-(3-de)*hgh2),Vector2D(hgw-0.5*hgw2,-2*hgh2),colground);
				ge5.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge5, 0);
				GeomLine2D ge6(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-2*hgh2),Vector2D(hgw+0.5*hgw2+hgdist,-(3-de)*hgh2),colground);
				ge6.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge6, 0);
				GeomLine2D ge7(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist-ww*hgh2,-(3)*hgh2),Vector2D(hgw-0.5*hgw2+ww*hgh2,-(3)*hgh2),colground);
				ge7.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge7, 0);

				GeomLine2D ge8(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-(3-de)*hgh2),Vector2D(hgw+0.5*hgw2+hgdist-ww*hgh2,-(3)*hgh2),colground);
				ge8.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge8, 0);
				GeomLine2D ge9(mbs,0,Vector2D(hgw-0.5*hgw2+ww*hgh2,-(3)*hgh2),Vector2D(hgw-0.5*hgw2,-(3-de)*hgh2),colground);
				ge9.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge9, 0);
			}

			if (1)
			{
				//top cup:
				GeomLine2D ge5(mbs,0,Vector2D(0,0),Vector2D(0,hgh),colground);
				ge5.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge5, 0);
				GeomLine2D ge6(mbs,0,Vector2D(hgw*2+hgdist,hgh),Vector2D(hgw*2+hgdist,0),colground);
				ge6.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge6, 0);
				GeomLine2D ge7(mbs,0,Vector2D(0,hgh),Vector2D(hgw*2+hgdist,hgh),colground);
				ge7.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge7, 0);
			}

			if (0)
			{
				hgw2*=3;
				hgh*=1.25;
				//second bottom cup:
				GeomLine2D ge5(mbs,0,Vector2D(hgw-0.5*hgw2,-3*hgh),Vector2D(hgw-0.5*hgw2,-2.5*hgh),colground);
				ge5.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge5, 0);
				GeomLine2D ge6(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-2.5*hgh),Vector2D(hgw+0.5*hgw2+hgdist,-3*hgh),colground);
				ge6.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge6, 0);
				GeomLine2D ge7(mbs,0,Vector2D(hgw+0.5*hgw2+hgdist,-3*hgh),Vector2D(hgw-0.5*hgw2,-3*hgh),colground);
				ge7.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge7, 0);
			}
		}
		else
		{
			//avalanche
			cline *= 0.1;

			double s = 2*hgh;
			double w = 2*hgw;
			double ha = 0.2;
			double hb = 0.05;
			double ho = 0.2*s*0;

			Vector2D ps[50];
			int x = 0;
			ps[x++] = Vector2D(w*1.5,s*0.25-ho);
			ps[x++] = Vector2D(-w,-s-ho);
			ps[x++] = Vector2D(-w,-s+ha-ho);
			ps[x++] = Vector2D(-hb-w,-s+ha-ho);
			ps[x++] = Vector2D(-hb-w,-s-ho);
			ps[x++] = Vector2D(-hb-1.5*w,-s*1.25-ho);

			double offx = -hb-0.5*w;
			double offy = -s*0.25-ha*0;
			ha *= 0;
			ps[x++] = Vector2D(-w+offx,-s+ha-ho+offy);
			ps[x++] = Vector2D(-hb-w+offx,-s+ha-ho+offy);
			ps[x++] = Vector2D(-hb-w+offx,-s-ho+offy);
			ps[x++] = Vector2D(-hb-1.5*w+offx,-s*1.25-ho+offy);

			ps[x++] = Vector2D(-hb-1.5*w-0.5*w+offx,-s*1.25-ho+offy);
			ps[x++] = Vector2D(-hb-1.75*w-0.5*w+offx,-s*1.25+0.25*s-ho+offy);

			ps[x++] = Vector2D(-hb-2*w-0.5*w+offx,-s*1.5-ho+offy);
			ps[x++] = Vector2D(-hb-2*w-1*w+offx,-s*1.5-ho+offy);
			ps[x++] = Vector2D(-hb-2*w-1*w+offx,-s*1.5+s-ho+offy);

			int n=7+4+3;
			for (int ii = 0; ii < n; ii++)
			{
				GeomLine2D ge1(mbs,0,ps[ii],ps[ii+1],colground);
				ge1.SetDrawParam(Vector3D(cline,res,0));
				gc.AddMasterGeomElement(&ge1, 0);
			}
		}

		gc.FinishContactDefinition();
		mbs->AddElement(&gc);
	}

	mbs->Assemble();

	return 1;
};


class Models_AutoRegistration_hourglass_rigid_contact_2D
{
public:
	Models_AutoRegistration_hourglass_rigid_contact_2D()
	{
		ModelFunctionAutoRegistration(&Generate_Model_hourglass_rigid_contact_2D, "hourglass_rigid_contact_2D", "hourglass_rigid_contact_2D", 0, &Generate_InitModelData_hourglass_rigid_contact_2D); 
	}
};
Models_AutoRegistration_hourglass_rigid_contact_2D dummy_hourglass_rigid_contact_2D;