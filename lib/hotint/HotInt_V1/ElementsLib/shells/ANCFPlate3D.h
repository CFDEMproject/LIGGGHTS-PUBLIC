//#**************************************************************
//#
//# filename:             ANCFPlate3D.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						17.October 2004
//# description:          3D Element Library
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
 
#ifndef ANCFPlate3D__H
#define ANCFPlate3D__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFPlate3D ANCFPlate3D ANCFPlate3D ANCFPlate3D ANCFPlate3D ANCFPlate3D 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


const int ANCFPlateMaxIP=120; //7*4*4 or 4*5*5 for plate ??
//rigid cube
class ANCFPlate3D: public Body3D
{
public:
	//Body3D():Element() {mbs = NULL;};
	ANCFPlate3D(MBS* mbsi):Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T(), xg(), xgd(), e0() {};
	ANCFPlate3D(const ANCFPlate3D& e):Body3D(e.mbs),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T(), xg(), xgd(), e0() {CopyFrom(e);};

	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)
	//initial velocities!!!!!
	ANCFPlate3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& xc3, const Vector& xc4, 
		const Vector& vc1, const Vector& vc2, const Vector& vc3, const Vector& vc4, 
		double rhoi, double Emi, double nui, const Vector3D& si, const Vector3D& coli):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T(), xg(), xgd(), e0()
	{
		node[0]=0; node[1]=0; node[2]=0; node[3]=0;
		sos2=SOS();
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;

		int sadd = SOS()/NNodes()-4*Dim();
		x_init = (((xc1.Append(xc2)).Append(xc3)).Append(xc4)).Append(Vector(sadd*NNodes()));
		xg = x_init;
		Vector v_init = (((vc1.Append(vc2)).Append(vc3)).Append(vc4)).Append(Vector(sadd*NNodes()));

		nu = nui;
		Em = Emi;
		//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
		col = coli;
		concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

		BuildDSMatrices();
		Matrix Tinv = T;
		if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
		e0 = x_init;
		x_init = Tinv * x_init;
		v_init = Tinv * v_init;
		//x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
		x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!
	};

	//Element shares nodes with other elements, n1i, n2i, n3i, n4i are nodenumbers; element sets initial conditions for nodes
	ANCFPlate3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& xc3, const Vector& xc4, 
		int n1i, int n2i, int n3i, int n4i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T(), xg(), xgd(), e0()
	{
		node[0]=n1i; node[1]=n2i; node[2]=n3i; node[3]=n4i;
		sos2=(NS()-4*NNodes())*Dim();
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;

		int sadd = SOS()/NNodes()-4*Dim();
		x_init = (((xc1.Append(xc2)).Append(xc3)).Append(xc4)).Append(Vector(sadd*NNodes()));
		xg = x_init;

		nu = nui;
		Em = Emi;
		//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
		col = coli;
		concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

		BuildDSMatrices();
		Matrix Tinv = T;
		if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
		e0 = x_init;
		x_init = Tinv * x_init;
		x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
	};

	//Element shares nodes with other elements, n1i, n2i, n3i, n4i are nodenumbers; element sets initial conditions for nodes
	//for a flat plate element, nodal gradients are computed for a flat plate
	ANCFPlate3D(MBS* mbsi, const Vector3D& xc1, const Vector3D& xc2, const Vector3D& xc3, const Vector3D& xc4, 
		const Vector3D& vc1, const Vector3D& vc2, const Vector3D& vc3, const Vector3D& vc4,
		int n1i, int n2i, int n3i, int n4i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T(), xg(), xgd(), e0()
	{
		node[0]=n1i; node[1]=n2i; node[2]=n3i; node[3]=n4i;
		sos2=(NS()-4*NNodes())*Dim();
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;

		int sadd = SOS()/NNodes()-4*Dim();
		Vector3D xx[4];
		xx[0] = xc1;
		xx[1] = xc2;
		xx[2] = xc3;
		xx[3] = xc4;
		Vector nc[4]; //nodal coordinate, 12 components
		for (int i=0; i < 4; i++)
		{
			nc[i].Init(); nc[i].SetLen(12);
		}
		Vector3D nvx = xx[1]-xx[0]; //compute gradients for straight plate
		Vector3D nvy = xx[3]-xx[0];
		nvx.Normalize();
		nvy.Normalize();
		Vector3D nvz = nvx.Cross(nvy);
		for (int i=0; i < 4; i++)
		{
			nc[i]( 1) = xx[i](1);nc[i]( 2) = xx[i](2);nc[i]( 3) = xx[i](3);
			nc[i]( 4) = nvx(1);  nc[i]( 5) = nvx(2);  nc[i]( 6) = nvx(3); 
			nc[i]( 7) = nvy(1);  nc[i]( 8) = nvy(2);  nc[i]( 9) = nvy(3); 
			nc[i](10) = nvz(1);  nc[i](11) = nvz(2);  nc[i](12) = nvz(3); 
		}


		Vector nv1(12); //nodal velocities
		Vector nv2(12);
		Vector nv3(12);
		Vector nv4(12);
		nv1(1) = vc1(1); nv1(2) = vc1(2); nv1(3) = vc1(3); //set initial velocities, no initial gradient derivatives!
		nv2(1) = vc2(1); nv2(2) = vc2(2); nv2(3) = vc2(3);
		nv3(1) = vc3(1); nv3(2) = vc3(2); nv3(3) = vc3(3);
		nv4(1) = vc4(1); nv4(2) = vc4(2); nv4(3) = vc4(3);

		x_init = (((nc[0].Append(nc[1])).Append(nc[2])).Append(nc[3])).Append(Vector(sadd*NNodes()));
		xg = x_init;

		nu = nui;
		Em = Emi;
		rho = rhoi;
		col = coli;
		concentratedmass[0] = concentratedmass[1] = concentratedmass[2] = concentratedmass[3] = 0;

		BuildDSMatrices();
		Matrix Tinv = T;
		if (!Tinv.Invert2()) {UO() << "ERROR: T matrix not invertible!!!\n";}
		e0 = x_init;

		Vector v_init = (((nv1.Append(nv2)).Append(nv3)).Append(nv4)).Append(Vector(sadd*NNodes()));

		x_init = Tinv * x_init;
		v_init = Tinv * v_init;

		x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFPlate3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body3D::CopyFrom(e);
		const ANCFPlate3D& ce = (const ANCFPlate3D&)e;
		lx = ce.lx;
		ly = ce.ly;
		lz = ce.lz;
		Em = ce.Em;
		nu = ce.nu;
		//size = ce.size;
		xg = ce.xg;
		xgd = ce.xgd;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		temp = ce.temp;
		SV = ce.SV;
		DS = ce.DS;

		//integration points
		x1 = ce.x1; x2 = ce.x2; x3 = ce.x3; 
		w1 = ce.w1; w2 = ce.w2; w3 = ce.w3; 
		orderxy = ce.orderxy;
		orderz = ce.orderz;

		for (int i=0; i < ANCFPlateMaxIP; i++)
		{
			grad[i] = ce.grad[i];
			jacdet[i] = ce.jacdet[i];
		}
		T = ce.T;
		e0 = ce.e0;

		sos2 = ce.sos2;
		for (int i = 0; i < NNodes(); i++)
		{
			Ti[i] = ce.Ti[i];
			jac[i] = ce.jac[i];
			node[i] = ce.node[i];
			concentratedmass[i] = ce.concentratedmass[i];
		}

		rho = ce.rho; //DR 2013-02-04 deleted rho from class element
	}

	virtual void Initialize() 
	{
		Body3D::Initialize();
	}
	virtual void LinkToElements();
	virtual void BuildDSMatrices(); 

	virtual int SOS() const {return Dim()*NS();}; //size of K and M
	virtual int SOSowned() const {return sos2;}; //len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	virtual int IsRigid() const {return 0;}
	virtual void SetConcentratedMass(double cm, int i) {concentratedmass[i-1] = cm;}

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//change for other shape functions and dimension:
	virtual int Dim() const {return 3;} //default value
	virtual int NS() const;
	virtual int NNodes() const {return 4;}
	virtual int NodeSize() const {return 12;}
	virtual void GetS0(Vector& sf, const Vector3D& ploc) const;
	virtual void GetS0x(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sf) const;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//compute v = T*v
	virtual void ApplyT(Vector& v) const;

	virtual void ApplyTtp(Vector& v) const;

	virtual void ApplyTtp(Matrix& m) const;

	virtual void GetSMatrix(const Vector& sf, Matrix& sm) const
	{
		for (int j = 1; j <= Dim(); j++)
			for (int i = 1; i<=sf.Length(); i++)
				sm(j,i) = sf(i);
	}

	virtual void GetCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XG(i);
	}

	virtual void GetCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGP(i);
	}

	virtual void GetDrawCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGD(i);
	}

	virtual void GetDrawCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGPD(i);
	}

	virtual Vector3D GetPos(const Vector3D& p_loc) const;
	virtual Vector3D GetPos(const Vector3D& p_loc, const Vector& xg) const;

	virtual Vector3D GetVel(const Vector3D& p_loc) const;

	//from -1 to +1
	virtual Vector3D GetPosx(const Vector3D& p_loc, const Vector& xg) const;

	virtual void GetRot(Matrix3D& rot, const Vector& xg) const;

	virtual Vector3D GetPosD(const Vector3D& p_loc) const;

	//in reference element coordinates (-1..1)
	virtual Vector3D GetPos0D(const Vector3D& p_loc) const;

	virtual Vector3D GetVelD(const Vector3D& p_loc) const;

	virtual void SetComputeCoordinates()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}

	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);

	virtual void GetJacobi(Matrix3D& jac, const Vector3D& p, const Matrix& DS, const Vector& x0) const
	{
		int ns = NS();
		for (int j = 1; j <= Dim(); j++)
		{
			for (int i = 1; i <= Dim(); i++)
			{ 
				jac(i,j) = 0;
				for (int k=1; k <= ns; k++)
				{ 
					jac(i,j) += DS(j,k)*x0((k-1)*Dim()+i);
				}
			}
		}
		//global_uo << "jac=" << jac << "\n";
	}

	virtual void EvalF2(Vector& f, double t); 

	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetCoordinates(u);

		ApplyT(u);
		u -= e0;
		Vector3D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly,0.5*lz);
		Matrix3D G;
		Gradu(p0, u, G);
		return G+Matrix3D(1);
	}
	virtual Matrix3D GetRotMatrixP(const Vector3D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetCoordinatesP(u);

		ApplyT(u);
		Vector3D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly,0.5*lz);
		Matrix3D G;
		Gradu(p0, u, G);
		return G;
	}
	virtual Matrix3D GetRotMatrixD(const Vector3D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetDrawCoordinates(u);

		u = T*u;
		u -= e0;
		Vector3D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly,0.5*lz);
		Matrix3D G;
		GraduD(p0, u, G);
		return G+Matrix3D(1);
	}

	//ploc -1 ... +1
	virtual void Gradu(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const
	{
		static Matrix DS;
		GetDSMatrix0(ploc,DS);

		Matrix3D jac, jacinv;
		GetJacobi(jac,ploc,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		static Matrix grad;
		grad.SetSize(Dim(),NS());
		Mult(jacinv, DS, grad);

		gradu.SetAll(0);
		int dim = Dim();
		int l;
		for (int j = 1; j <= dim; j++) 
		{
			for (int i = 1; i <= NS(); i++)
			{
				l = (i-1)*dim+j;
				for (int k = 1; k <= dim; k++)
				{
					gradu(j,k) += grad(k,i)*u(l);
				}
			}
		}
	}

		//for body loads:
	//Computes f = d p_ref/d q * x
	virtual void ApplyDprefdq(Vector& f, const Vector3D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";

	}
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector3D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}
	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq)
	{
		GetH(dudq); 
	}
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
	{
		d.SetSize(NS()*Dim(),Dim());
		static Matrix d2;
		d2.SetSize(NS()*Dim(),Dim());
		double diffpar = 1e-8;
		GetdPosdqT(ploc+diffpar*vloc, d2);
		d = d2;
		GetdPosdqT(ploc,d2);
		d -= d2;
		d *= 1./diffpar;

	}

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d)
	{
		//p = S(p.x,p.y,p.z)*q; d(p)/dq 
		Vector3D p0=ploc;
		p0.Scale(0.5*lx,0.5*ly,0.5*lz);
		static Vector SV;
		GetS0(SV, p0);
		d.SetSize(NS()*Dim(),Dim());
		d.FillWithZeros();
		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*3+1,1)=SV(i);
			d((i-1)*3+2,2)=SV(i);
			d((i-1)*3+3,3)=SV(i);
		}
		ApplyTtp(d);
		//UO() << "Not yet implemented\n";
	}

	void GetdRotdqT(const Vector3D& ploc, Matrix& d)
	{
		//***********************************
		//vector rz for y-rotation
		//vector rz for x-rotation
		//vector rx/ry for z-rotation

		//damit gleich wie 
		//Compute Gradient ...


		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		ApplyT(xg);

		Vector3D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly,0.5*lz);

		static Matrix DS;
		GetDSMatrix0(p0,DS);

		Matrix3D jac, jacinv;
		GetJacobi(jac,p0,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		static Matrix grad;
		Matrix3D F;
		grad.SetSize(Dim(),NS());
		Mult(jacinv, DS, grad);

		F.SetAll(0);
		int dim = Dim();
		int l;
		for (int j = 1; j <= dim; j++) 
		{
			for (int i = 1; i <= NS(); i++)
			{
				l = (i-1)*dim+j;
				for (int k = 1; k <= dim; k++)
				{
					F(j,k) += grad(k,i)*xg(l);
				}
			}
		}

		static Vector Sx, Sy, Sz;
		Sx.SetLen(NS());
		Sy.SetLen(NS());
		Sz.SetLen(NS());
		for (int i = 1; i <= NS(); i++)
		{
			Sx(i) = grad(1,i);
			Sy(i) = grad(2,i);
			Sz(i) = grad(3,i);
		}


		Vector3D rx(F(1,1),F(2,1),F(3,1));
		Vector3D ry(F(1,2),F(2,2),F(3,2));
		Vector3D rz(F(1,3),F(2,3),F(3,3));

		double dyx = Sqr(ry.X())+Sqr(ry.Y());
		double dxy = Sqr(rx.Y())+Sqr(rx.X());

		double dzx = Sqr(rz.X())+Sqr(rz.Z());
		//double dxz = Sqr(rx.Z())+Sqr(rx.X());

		double dzy = Sqr(rz.Y())+Sqr(rz.Z());
		//double dyz = Sqr(ry.Y())+Sqr(ry.Z());

		for (int i = 1; i <= NS(); i++)
		{
			//moment x/y is applied to surface normal vector r,z
			//from delta gamma_x
			d((i-1)*3+1,1) = 0; //only delta r.X() terms!
			d((i-1)*3+2,1) =-(rz.Z()*Sz(i))/dzy; //only delta r.Y() terms!
			d((i-1)*3+3,1) = (rz.Y()*Sz(i))/dzy; //only delta r.Z() terms!
			
			//from delta gamma_y
			d((i-1)*3+1,2) = (rz.Z()*Sz(i))/dzx; //only delta r.X() terms!
			d((i-1)*3+2,2) = 0; //only delta r.Y() terms!
			d((i-1)*3+3,2) =-(rz.X()*Sz(i))/dzx; //only delta r.Z() terms!
			
			//from delta gamma_z
			d((i-1)*3+1,3) =-(ry.Y()*Sy(i))/dyx - (rx.Y()*Sx(i))/dyx; //only delta r.X() terms!
			d((i-1)*3+2,3) = (ry.X()*Sy(i))/dyx + (rx.X()*Sx(i))/dyx; //only delta r.Y() terms!
			d((i-1)*3+3,3) = 0; //only delta r.Z() terms!

		}
		//ApplyTtp(d);???

	}


	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) 
	{
		Vector3D p0(ploc.X()/(0.5*lx),0,0);

		SV.SetLen(NS());
		GetS0x(SV,p0);
		
		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		ApplyT(xg);
		dpdx = 0;
		for (int i=1; i <= 3; i++)
		{
			for (int j=1; j <= NS(); j++)
			{
				dpdx(i) += SV(j) * xg((j-1)*3+i);
			}
		}
	};

	//ploc -1 ... +1
	virtual void GraduD(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const
	{
		static Matrix DS;
		GetDSMatrix0(ploc,DS);

		Matrix3D jac, jacinv;
		GetJacobi(jac,ploc,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		static Matrix grad;
		grad.SetSize(Dim(),NS());
		Mult(jacinv, DS, grad);

		gradu.SetAll(0);
		int dim = Dim();
		int l;
		for (int j = 1; j <= dim; j++) 
		{
			for (int i = 1; i <= NS(); i++)
			{
				l = (i-1)*dim+j;
				for (int k = 1; k <= dim; k++)
				{
					gradu(j,k) += grad(k,i)*u(l);
				}
			}
		}
	}

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

	virtual void DrawElement();

protected:
	//mechanical:
	double lx, ly, lz, Em, nu;

	int node[4];
	int sos2;
	double concentratedmass[4];

	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix;
	Matrix Hmatrix;
	Matrix DS; //3*8
	Vector SV; //8
	Vector temp;

	//integration points
	Vector x1,x2,x3,w1,w2,w3;
	int orderxy, orderz;
	Matrix grad[ANCFPlateMaxIP];
	double jacdet[ANCFPlateMaxIP];

	Matrix Ti[4];
	Matrix T; //slope discontinuities
	Matrix3D jac[4]; //slope discontinuities
	Vector e0; //initial vector in e, not in p

	double rho; //DR 2013-02-04 deleted rho from class element
};




#endif

/*
		B.SetSize(ns*dim,9);
		B.SetAll(0);
		Vector piola1v;
		Matrix B;
*/