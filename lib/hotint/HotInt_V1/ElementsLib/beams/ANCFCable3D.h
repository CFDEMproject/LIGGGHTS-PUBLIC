//#**************************************************************
//#
//# filename:             ANCFCable3D.h
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
 
#ifndef ANCFCABLE3D__H
#define ANCFCABLE3D__H

#include "body3d.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFCable3D ANCFCable3D ANCFCable3D ANCFCable3D ANCFCable3D ANCFCable3D ANCFCable3D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


const int ANCFCableMaxIP=45;
//rigid cube
class ANCFCable3D: public Body3D    //$EDC$[beginclass,classname=ANCFCable3D,parentclassname=Body3D]
{
public:
	//Body3D():Element() {mbs = NULL;};
	ANCFCable3D(MBS* mbsi):Body3D(mbsi), massmatrix(), Hmatrix(), SV(), x1(), w1(),
			xg(), xgd(), e0() { SetElementName("ANCFCable3D"); };
	ANCFCable3D(const ANCFCable3D& e):Body3D(e.mbs),massmatrix(), Hmatrix(), SV(), x1(), w1(),
			xg(), xgd(), e0() {CopyFrom(e);};
	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)
	ANCFCable3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0, const Vector3D& directori=Vector3D(0,0,0)):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), x1(), w1(),
			xg(), xgd(), e0()
	{
		SetElementName("ANCFCable3D");
		n1=0; n2=0; sos2=SOS();
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		x_init = xc1.Append(xc2);
		xg = xc1.Append(xc2);

		Em = Emi;
		//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		BuildDSMatrices();
		e0 = x_init;
		x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!

		if (directori.X() != 0 || directori.Y() != 0 || directori.Z() != 0)
		{
			director.Set(directori.X(),directori.Y(),directori.Z());
			use_director = true;
		}
		else
		{
			use_director = false;
		}
	};

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFCable3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0, const Vector3D& directori=Vector3D(0,0,0)):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), x1(), w1(),
			xg(), xgd(), e0()
	{
		SetElementName("ANCFCable3D");
		n1=n1i; n2=n2i; sos2=0;
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		x_init = xc1.Append(xc2);
		xg = xc1.Append(xc2);

		Em = Emi;
		//rho = rhoi; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		//UO() << "Cable: rho=" << rho << ", E=" << Em << ", size=" << size 	<< ", n1=" << n1 << ", n2=" << n2 << "\n";

		BuildDSMatrices();
		e0 = x_init;
		x_init = x_init.Append(Vector(SOS())); //Velocity initial conditions can also be transformed by Tinv!!!
		
		if (directori.X() != 0 || directori.Y() != 0 || directori.Z() != 0)
		{
			director.Set(directori.X(),directori.Y(),directori.Z());
			use_director = true;
		}
		else
		{
			use_director = false;
		}
	};

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFCable3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
		int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0, const Vector3D& directori=Vector3D(0,0,0)):
	  Body3D(mbsi),massmatrix(), Hmatrix(), SV(), x1(), w1(),
			xg(), xgd(), e0()
	{
		SetElementName("ANCFCable3D");
		n1=n1i; n2=n2i; sos2=0;
		size = si;
		lx = si.X(); ly = si.Y(); lz = si.Z();

		mass = lx*ly*lz*rhoi;
		//lx=1;ly=1;lz=1;

		x_init = xc1.Append(xc2);
		xg = xc1.Append(xc2);
		Vector v_init = vc1.Append(vc2);

		Em = Emi;
		rho = rhoi;
		col = coli;
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		BuildDSMatrices();
		e0 = x_init;
		x_init = x_init.Append(v_init); //Velocity initial conditions can also be transformed by Tinv!!!

		if (directori.X() != 0 || directori.Y() != 0 || directori.Z() != 0)
		{
			director.Set(directori.X(),directori.Y(),directori.Z());
			use_director = true;
		}
		else
		{
			use_director = false;
		}
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFCable3D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body3D::CopyFrom(e);
		const ANCFCable3D& ce = (const ANCFCable3D&)e;
		lx = ce.lx;
		ly = ce.ly;
		lz = ce.lz;
		Em = ce.Em;
		xg = ce.xg;
		xgd = ce.xgd;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		SV = ce.SV;

		rho= ce.rho; //DR 2013-02-04 deleted rho from class element

		//integration points
		x1 = ce.x1; 
		w1 = ce.w1;
		orderx = ce.orderx;

		e0 = ce.e0;

		sos2 = ce.sos2;
		n1 = ce.n1; n2 = ce.n2;
		concentratedmass1 = ce.concentratedmass1;
		concentratedmass2 = ce.concentratedmass2;
		
		use_director = ce.use_director;
		if (use_director) 
		{
			director = ce.director;
		}
	}

	virtual void Initialize() 
	{
		Body3D::Initialize();
	}
	virtual void LinkToElements();
	virtual void BuildDSMatrices(); 

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const char* GetElementSpec() const { return "ANCFCable3D"; }

	virtual int SOS() const {return 12;}; //size of K and M
	virtual int SOSowned() const {return sos2;}; //len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	virtual int IsRigid() const {return 0;}

	virtual void SetConcentratedMass(double cm1, double cm2) {concentratedmass1=cm1; concentratedmass2=cm2;}

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+12));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+12));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+12));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+12));}

	virtual int NS() const {return 4;}
	virtual int NNodes() const {if (SOSowned() == 0) return 2; else return 0;};
	virtual const int& NodeNum(int i) const 
	{
		if (i == 1) return n1; else return n2;
	}
	virtual int& NodeNum(int i) 
	{
		if (i == 1) return n1; else return n2;
	}

	double& Rho() {return rho;} 
	const double& Rho() const {return rho;}


	virtual void GetS0(Vector& sf, const double& ploc) const;
	virtual void GetS0x(Vector& sfx, const double& ploc) const;
	virtual void GetS0xx(Vector& sfx, const double& ploc) const;

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

	virtual Vector3D GetPos(const double& p_loc) const;

	virtual Vector3D GetPosx0(const double& p_loc, const Vector& xg) const
	{
		double p0=p_loc;
		static Vector SV;
		GetS0x(SV, p0);
		Vector3D p(0.,0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*3+i);
			}
		}
		return p;
	}

	virtual Vector3D GetPosx0D(const double& p_loc) const
	{
		double p0=p_loc;
		static Vector SV;
		GetS0x(SV, p0);
		static Vector xgd;
		xgd.SetLen(SOS());
		GetDrawCoordinates(xgd);

		Vector3D p(0.,0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xgd((j-1)*3+i);
			}
		}
		return p;
	}

	virtual Vector3D GetPosxx0(const double& p_loc, const Vector& xg) const
	{
		double p0=p_loc;
		static Vector SV;
		GetS0xx(SV, p0);
		Vector3D p(0.,0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*3+i);
			}
		}
		return p;
	}

	virtual Vector3D GetVel(const double& p_loc) const;
	virtual Vector3D GetPosD(const double& p_loc) const;
	virtual Vector3D GetRefPosD() const {	return GetPosD(0); }     //PG & KN
	virtual Vector3D GetDOFPosD(int idof) const
	{
		//returns postion of i-th DOF
		if (idof < 7)
		{
			return GetPosD(-lx/2.);
		}
		return GetPosD(lx/2.);
	}

	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		Matrix3D rot;
		if (idof < 7)
			rot = Matrix3D(1.);//GetInitRotMatrix3D(-lx*.5);
		else
			rot = Matrix3D(1.);//GetInitRotMatrix3D(lx*.5);

		switch(idof)
		{
		case 1: case 7: return Vector3D(rot(1,1),rot(2,1),rot(3,1)); break;
		case 2: case 8: return Vector3D(rot(1,2),rot(2,2),rot(3,2)); break;
		case 3: case 9: return Vector3D(rot(1,3),rot(2,3),rot(3,3)); break;
		}
		return Vector3D(0.,0.,0.);
	}

	//in reference element coordinates (-1..1)
	virtual Vector3D GetPos0D(const double& p_loc) const;

	virtual Vector3D GetVelD(const double& p_loc) const;

	virtual Vector3D GetPos(const Vector3D& p_loc) const
	{
		if (use_director)
		{
			return GetPos(p_loc.X()) + GetRotMatrix(p_loc.X())*Vector3D(0.,p_loc.Y(),p_loc.Z());
		}
		return GetPos(p_loc.X());
	}
	virtual Vector3D GetPosD(const Vector3D& p_loc) const
	{
		if (use_director)
		{
			return GetPosD(p_loc.X()) + GetRotMatrixD(p_loc.X())*Vector3D(0.,p_loc.Y(),p_loc.Z());
		}
		return GetPosD(p_loc.X());
	}
	virtual Vector3D GetVel(const Vector3D& p_loc) const {return GetVel(p_loc.X());}
	virtual Vector3D GetVelD(const Vector3D& p_loc) const {return GetVelD(p_loc.X());}

	virtual void SetComputeCoordinates()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}

	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);

	virtual void EvalF2(Vector& f, double t); 

	virtual double GetKappa(const double& x, const Vector& xg) const;
	virtual double GetEpsAxial(const double& x, const Vector& xg) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const;
	virtual void GetDeltaEpsAxial(const double& x, const Vector& xg, Vector& depsaxial) const;

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
	virtual Vector3D GetRotMatv(const Vector3D& ploc, const Vector3D& v)
	{
		double diffpar = 1e-8;
		Vector3D vrot = (1./diffpar)*(GetPos(ploc + diffpar*v) - GetPos(ploc));
		return vrot;

		/*
		Vector3D vrot = GetPos(ploc + diffpar*v) - GetPos(ploc);
		vrot.Normalize();
		return vrot*v.Norm();*/
	}
	virtual Vector3D GetRotMatPv(const Vector3D& ploc, const Vector3D& v)
	{
		double diffpar = 1e-8;
		Vector3D vrot = (1./diffpar)*(GetVel(ploc + diffpar*v) - GetVel(ploc));
		return vrot;
		//vrot.Normalize();
		//return vrot*v.Norm();
	}
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
	{
		//Only vloc = (1,0,0) is possible!!!!!!!!!!!!!!!!!

		d.SetSize(NS()*Dim(),3);
		static Matrix d2;
		d2.SetSize(NS()*Dim(),3);
		double diffpar = 1e-8;
		GetdPosdqT(ploc+diffpar*vloc, d2);
		d = d2;
		//UO() << "d=" << d << "\n";
		GetdPosdqT(ploc,d2);
		d -= d2;
		d *= 1./diffpar;
	}

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d)
	{
		//p = S(p.x,p.y,p.z)*q; d(p)/dq 
		static Vector SV;
		GetS0(SV, ploc.X()/(0.5*lx));
		d.SetSize(NS()*Dim(),3);
		d.FillWithZeros();
		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*3+1,1)=SV(i);
			d((i-1)*3+2,2)=SV(i);
			d((i-1)*3+3,3)=SV(i);
		}
	}

	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) 
	{
		//optimierungspotenzial 500% !!!!!!!!!!!!!!!!!!!
		double p0 = ploc.X()/(0.5*lx);

		SV.SetLen(NS());
		GetS0x(SV,p0);
		
		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		dpdx = 0;
		for (int i=1; i <= 3; i++)
		{
			for (int j=1; j <= NS(); j++)
			{
				dpdx(i) += SV(j) * xg((j-1)*3+i);
			}
		}
	};

	virtual const Vector3D& GetDirector() const {return director;}
	virtual Vector3D& GetDirector() {return director;}
	virtual Matrix3D GetRotMatrix(const double p_loc) const
	{
		// use only, if director is defined (see constructor)
		assert (use_director);

		// p_loc in [-lx/2,lx/2]
		// p0 in [-1,1]
		double p0=p_loc/(0.5*lx);

		static Vector xg(SOS());
		for (int j=1; j<=SOS(); j++)
		{
			xg(j) = XG(j);
		}

		Vector3D rx = GetPosx0(p0, xg);
		rx.Normalize();

		Vector3D e20 = GetDirector();
		rx.GramSchmidt(e20); //d is projected into plane normal to rx (GramSchmidt() already returns a normalized Vector3D)

		Vector3D e30 = rx.Cross(e20); //length=1

		Matrix3D A(rx.X(),e20.X(),e30.X(),
			rx.Y(),e20.Y(),e30.Y(),
			rx.Z(),e20.Z(),e30.Z());
		
		return A;
	}
	virtual Matrix3D GetRotMatrixD(const double p_loc) const
	{		
		// use only, if director is defined (see constructor)
		assert (use_director);

		// p_loc in [-lx/2,lx/2]
		// p0 in [-1,1]
		double p0=p_loc/(0.5*lx);

		Vector3D rx = GetPosx0D(p0);
		rx.Normalize();

		Vector3D e20 = GetDirector();
		rx.GramSchmidt(e20); //d is projected into plane normal to rx (GramSchmidt() already returns a normalized Vector3D)

		Vector3D e30 = rx.Cross(e20); //length=1

		Matrix3D A(rx.X(),e20.X(),e30.X(),
			rx.Y(),e20.Y(),e30.Y(),
			rx.Z(),e20.Z(),e30.Z());
		
		return A;
	}

	virtual void DrawElement();

protected:
	//mechanical:
	double lx, ly, lz, Em;

	int n1, n2, sos2;

	double concentratedmass1, concentratedmass2;
	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix, Hmatrix;
	Vector SV; //4
	Vector temp;

	//integration points
	Vector x1,w1;
	int orderx;

	Vector e0; //initial vector in e, not in p

	Vector3D director;
	bool use_director;

	double rho; //DR 2013-02-04 deleted rho from class element

};   //$EDC$[endclass,ANCFCable3D]




#endif
