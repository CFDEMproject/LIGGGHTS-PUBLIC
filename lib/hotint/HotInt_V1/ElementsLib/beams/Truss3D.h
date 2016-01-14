//#**************************************************************
//#
//# filename:             Truss3D.h
//#
//# author:               Michael Stangl
//#
//# generated:			  13.01.2012
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
 
#ifndef TRUSS3D__H
#define TRUSS3D__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Truss3D Truss3D Truss3D Truss3D Truss3D Truss3D Truss3D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class Truss3D: public ANCFCable3D
{
public:
	//ANCFCable3D():Element() {mbs = NULL;};
	Truss3D(MBS* mbsi):ANCFCable3D(mbsi) {};
	Truss3D(const Truss3D& e):ANCFCable3D(e.mbs) {CopyFrom(e);};
	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)

	Truss3D(MBS* mbsi, const Vector3D& xc1, const Vector3D& xc2, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0):
	  ANCFCable3D(mbsi)
	{
		n1=0; n2=0; sos2=SOS();
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		Vector3D test1(xc1(1),xc1(2),xc1(3));
		Vector3D test2(xc2(1),xc2(2),xc2(3));

		x_init = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));
		xg = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));
		
		//x_init = xc1.Append(xc2);
		//xg = xc1.Append(xc2);

		Em = Emi;
		rho = rhoi;
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		BuildDSMatrices();
		e0 = x_init;

		x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
	};


	////Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	//Truss3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
	//	const Vector3D& si, const Vector3D& coli, int beammodel = 0):
	//  ANCFCable3D(mbsi,xc1,xc2,n1i,n2i,rhoi,Emi,si,coli,beammodel)
	//{
	//	x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
	//};

	////Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	//Truss3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
	//	int n1i, int n2i, double rhoi, double Emi,
	//	const Vector3D& si, const Vector3D& coli, int beammodel = 0):
	//  ANCFCable3D(mbsi,xc1,xc2,vc1,vc2,n1i,n2i,rhoi,Emi,si,coli,beammodel)
	//{
	//};

	  	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	Truss3D(MBS* mbsi, const Vector3D& xc1, const Vector3D& xc2, int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0):
	  ANCFCable3D(mbsi)
	{
		n1=n1i; n2=n2i; sos2=0;
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		x_init = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));
		xg = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));

		Em = Emi;
		rho = rhoi;
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		BuildDSMatrices();
		e0 = x_init;
		x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
	};

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	Truss3D(MBS* mbsi, const Vector3D& xc1, const Vector3D& xc2, const Vector3D& vc1, const Vector3D& vc2, 
		int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0):
	  ANCFCable3D(mbsi)
	{
		n1=n1i; n2=n2i; sos2=0;
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		x_init = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));
		xg = Vector(xc1.X(),xc1.Y(),xc1.Z()).Append(Vector(xc2.X(),xc2.Y(),xc2.Z()));

		//x_init = xc1.Append(xc2);
		//xg = xc1.Append(xc2);

		Vector v_init = Vector(vc1.X(),vc1.Y(),vc1.Z()).Append(Vector(vc2.X(),vc2.Y(),vc2.Z()));

		Em = Emi;
		rho = rhoi;
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		BuildDSMatrices();
		e0 = x_init;
		x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Truss3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		ANCFCable3D::CopyFrom(e);
		const Truss3D& ce = (const Truss3D&)e;
	}

	virtual void Initialize() 
	{
		ANCFCable3D::Initialize();
	}

	virtual int SOS() const {return 6;}; //size of K and M
//	virtual int SOSowned() const {return sos2;}; //len(u)

	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+6));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+6));}
	
	virtual int NS() const {return 2;}

	virtual void GetS0(Vector& sf, const double& ploc) const;
	virtual void GetS0x(Vector& sfx, const double& ploc) const;
	virtual void GetS0xx(Vector& sfx, const double& ploc) const;

	virtual void EvalM(Matrix& m, double t);

	virtual void EvalF2(Vector& f, double t); 

	virtual double GetKappa(const double& x, const Vector& xg) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const;

	virtual Vector3D GetRotMatv(const Vector3D& ploc, const Vector3D& v)
	{
		assert(0 && "Truss3D: GetRotMatv called!");return Vector3D(0.);
	}
	virtual Vector3D GetRotMatPv(const Vector3D& ploc, const Vector3D& v)
	{
		assert(0 && "Truss3D: GetRotMatPv called!");return Vector3D(0.);
	}
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
	{
		assert(0 && "Truss3D: GetdRotvdqT called!");return;
	}

	virtual void DrawElement();
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

protected:
	
};
#endif
