//#**************************************************************
//#
//# filename:             ANCFBeam2D.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						23.October 2007
//# description:          2D Element Library
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
 
#ifndef ANCFBEAM2D__H 
#define ANCFBEAM2D__H

#include "body2d.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeam2D ANCFBeam2D ANCFBeam2D ANCFBeam2D ANCFBeam2D ANCFBeam2D ANCFBeam2D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const int ANCFBeam2D_nshapes = 6;
const int ANCFBeam2DMaxIP = 45;//45;
//rigid cube
class ANCFBeam2D: public Body2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	ANCFBeam2D(MBS* mbsi):Body2D(mbsi), massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), w1(), w2(),
			T1(), T2(), T(), xg(), xgd(), e0() {};
	ANCFBeam2D(const ANCFBeam2D& e):Body2D(e.mbs),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), w1(), w2(),
			T1(), T2(), T(), xg(), xgd(), e0() {CopyFrom(e);};
	
	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
		int n1i, int n2i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0);

	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeam2D(*this);
		return ec;
	}
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const ANCFBeam2D& ce = (const ANCFBeam2D&)e;
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
		SV = ce.SV;
		DS = ce.DS;
		elasticforce_beam = ce.elasticforce_beam;

		//integration points
		x1 = ce.x1; x2 = ce.x2; 
		w1 = ce.w1; w2 = ce.w2; 
		orderx = ce.orderx;
		orderyz = ce.orderyz;

		for (int i=0; i < ANCFBeam2DMaxIP; i++)
		{
			grad[i] = ce.grad[i];
			jacdet[i] = ce.jacdet[i];
		}
		T1 = ce.T1;
		T2 = ce.T2;
		T = ce.T;
		e0 = ce.e0;
		jac1 = ce.jac1;
		jac2 = ce.jac2;

		sos2 = ce.sos2;
		n1 = ce.n1; n2 = ce.n2;
		concentratedmass1 = ce.concentratedmass1;
		concentratedmass2 = ce.concentratedmass2;

		/*
		beamEA = ce.beamEA;
		beamEI = ce.beamEI;
		beamGAks = ce.beamGAks;
		beamRhoA = ce.beamRhoA;
		beamRhoI = ce.beamRhoI;*/

		orderWl = ce.orderWl;
		orderWs = ce.orderWs;
		orderWsHR = ce.orderWsHR;
		orderWb = ce.orderWb;

		rho = ce.rho; //DR 2013-02-04 deleted rho from class element
	}

	virtual void Initialize() 
	{
		Body2D::Initialize();
	}
	virtual void LinkToElements();
	virtual void BuildDSMatrices(); 
	virtual void InitConstructor() //initialize while constructor is called
	{
		concentratedmass1 = 0;
		concentratedmass2 = 0;
	}

	virtual const char* GetElementSpec() const {return "ANCFBeam2D";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual int NS() const {return ANCFBeam2D_nshapes;}
	virtual int NodeSize() const {return ANCFBeam2D_nshapes;}

	virtual int SOS() const {return 12;}; //size of K and M
	virtual int SOSowned() const {return sos2;}; //len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	virtual int NNodes() const {if (SOSowned() == 0) return 2; else return 0;};
	virtual const int& NodeNum(int i) const 
	{
		if (i == 1) return n1; else return n2;
	}
	virtual int& NodeNum(int i) 
	{
		if (i == 1) return n1; else return n2;
	}

	virtual int IsRigid() const {return 0;}
	virtual void SetConcentratedMass(double cm1, double cm2) {concentratedmass1=cm1; concentratedmass2=cm2;}

	virtual void SetInitialCurvatures(); //set initial vectors for curvature, stretch, shear deformation for Schwab-Meijaard model
	virtual void SetRectangularBeamParameters(double& EA, double& EI, double& GAks, double& rhoA, double& rhoI); //set beam parameters for rectangular cross-section, Schwab-Meijaard model

	//beam Parameters including shear correction factor!
	virtual void SetBeamParameters(double EA, double EI, double GAks, double rhoA, double rhoI);

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}

	//compute v = T*v
	virtual void ApplyT(Vector& v) const;
	virtual void ApplyTD(Vector& v) const;

	virtual void ApplyTtp(Vector& v) const;

	virtual void ApplyTtp(Matrix& m) const;

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;

	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;

	virtual void GetS0x(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0xlin(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0xlin2(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0y(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0yx(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0xx(Vector& sfx, const Vector2D& ploc) const;

	virtual double GetKappa(const double& x, const Vector& xg) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const;

	virtual void GetRot(Matrix3D& rot, const Vector& xg) const;

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

	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc, const Vector& xg, int flagD) const;

	virtual Vector2D GetPos2D(const Vector2D& p_loc) const;

	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;

	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const;

	//in reference element coordinates (-1..1)
	virtual Vector2D GetPos2D0D(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2D0Dinit(const Vector2D& p_loc) const; //for eigenmodes

	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	virtual Vector2D GetPosx2D(const Vector2D& p_loc) const;
	virtual Vector2D GetPosxx2D(const Vector2D& p_loc) const;

	virtual void SetComputeCoordinates()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}

	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);

	virtual void GetJacobi(Matrix3D& jac, const Vector2D& p, const Matrix& DS, const Vector& x0) const
	{
		jac.SetSize(2,2);
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

	virtual void GetDMatrix(Matrix3D& D, double nu, double em) const;

	virtual void EvalF2(Vector& f, double t); 


	virtual Matrix3D GetRotMatrix(const Vector2D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetCoordinates(u);

		ApplyT(u);
		u -= e0;
		Vector2D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly);
		Matrix3D G;
		Gradu(p0, u, G);
		G(1,1) += 1;
		G(2,2) += 1;
		return G;
	}
	virtual Matrix3D GetRotMatrixP(const Vector2D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetCoordinatesP(u);

		ApplyT(u);
		Vector2D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly);
		Matrix3D G;
		Gradu(p0, u, G);
		return G;
	}
	virtual Matrix3D GetRotMatrixD(const Vector2D& ploc) const 
	{
		//Compute Gradient ...
		static Vector u;
		u.SetLen(SOS());
		GetDrawCoordinates(u);

		u = T*u;
		u -= e0;
		Vector2D p0(ploc);
		p0.Scale(0.5*lx,0.5*ly);
		Matrix3D G;
		GraduD(p0, u, G);

		G(1,1) += 1;
		G(2,2) += 1;
		return G;
	}
	//???? check with Cable2D, Rigid2D ?????????????????????
	//approximate angular velocity:
	virtual double GetAngularVel(const Vector2D& p_loc) const 
	{
		Matrix3D omega_skew = GetRotMatrixP(p_loc)*GetRotMatrix(p_loc).GetTp();

		return omega_skew(2,1);
	}


	//ploc -1 ... +1
	virtual void Gradu(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const
	{
		static Matrix DS;
		GetDSMatrix0(ploc,DS);

		Matrix3D jac, jacinv;
		GetJacobi(jac,ploc,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();
		int dim = Dim();

		static Matrix grad;
		grad.SetSize(dim,NS());
		Mult(jacinv, DS, grad);

		gradu.SetSize(2,2);
		gradu.SetAll(0);
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
	virtual void ApplyDprefdq(Vector& f, const Vector2D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";

	}
	//Computes f = d rot_ref/d q * x, rot bedeutet rotation um x, y, und z-Achse
	virtual void ApplyDrotrefdq(Vector& f, const Vector2D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}
	virtual void ApplyDrotdq(Vector& f, const Vector2D& x)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}
	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq)
	{
		GetH(dudq);
		//UO() << "Not yet implemented\n";
	}
	virtual void GetdRotvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d)
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

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d)
	{
		//p = S(p.x,p.y,p.z)*q; d(p)/dq 
		Vector2D p0=ploc;
		p0.Scale(0.5*lx,0.5*ly);
		static Vector SV;
		GetS0(SV, p0);
		d.SetSize(NS()*Dim(),Dim());
		d.FillWithZeros();
		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*Dim()+1,1)=SV(i);
			d((i-1)*Dim()+2,2)=SV(i);
		}
		ApplyTtp(d);
	}

	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& d)
	{
		Vector2D p0=ploc;
		p0.Scale(0.5*lx,0.5*ly);

		d.SetSize(NS()*Dim(),1);
		d.FillWithZeros(); //needed???
		
		static Vector Sy; Sy.SetLen(NS()); 

		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		//ApplyT(xg);

		Vector2D ry;

		GetS0y(Sy,p0);

		for (int i = 1; i <= NS(); i++)
		{
			ry.X() += Sy(i)*xg(Dim()*(i-1)+1);
			ry.Y() += Sy(i)*xg(Dim()*(i-1)+2);
		}

		double dyx = Sqr(ry.X())+Sqr(ry.Y());

		for (int i = 1; i <= NS(); i++)
		{
			//from delta gamma_z
			d((i-1)*Dim()+1,1) =-(ry.Y()*Sy(i))/dyx; //only delta r.X() terms!
			d((i-1)*Dim()+2,1) = (ry.X()*Sy(i))/dyx; //only delta r.Y() terms!
		}
		//ApplyTtp(d);
	}

	virtual void GetdPosdx(const Vector2D& ploc, Vector2D& dpdx) 
	{
		//at centerline!
		//optimierungspotenzial 500% !!!!!!!!!!!!!!!!!!!
		double p0 = ploc.X()/(0.5*lx);

		SV.SetLen(NS());
		GetS0x(SV,p0);
		
		static Vector xg;
		xg.SetLen(SOS());
		GetCoordinates(xg);
		ApplyT(xg);

		dpdx = 0;
		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= NS(); j++)
			{
				dpdx(i) += SV(j) * xg((j-1)*Dim()+i);
			}
		}

	};

	//ploc -1 ... +1
	virtual void GraduD(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const;

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);

	virtual void DrawElement();

protected:
	//mechanical:
	double lx, ly, lz, Em, nu;

	int n1, n2, sos2;
	double concentratedmass1, concentratedmass2;

	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix, Hmatrix;
	Matrix DS; //2*6
	Vector SV; //6
	Vector temp;
	int elasticforce_beam;

	//integration points
	Vector x1,x2,w1,w2;
	int orderx, orderyz;
	Matrix grad[ANCFBeam2DMaxIP]; //DS transformed, grad u = grad[i]*e
	double jacdet[ANCFBeam2DMaxIP];

	Matrix T1, T2, T; //slope discontinuities
	Matrix3D jac1, jac2; //slope discontinuities
	Vector e0; //initial vector in e, not in p

	Matrix K0; //stiffness matrix at initial configuration

	/*
	double beamEA; //for beam formulation
	double beamEI;
	double beamGAks;
	double beamRhoA;
	double beamRhoI;*/

	double factstiffWl; //reduce stiffness of thickness modes, standard==1

	int orderWl;   //5; 5 and 6 coincide with 8 up to 10 digits, 7=8
	int orderWs;   //2; more than 2 gives locking, 1 is rather inaccurate
	int orderWsHR; //4; for Hellinger Reissner, 4 is exact
	int orderWb;   //3; 3=4=6, 2 gives shear locking

	double rho; //DR 2013-02-04 deleted rho from class element

};





//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ANCFBeam2Dlin: public ANCFBeam2D
{
public:
	ANCFBeam2Dlin(const ANCFBeam2Dlin& e):ANCFBeam2D(e.mbs) {CopyFrom(e);};
	
	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFBeam2Dlin(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, 
		int n1i, int n2i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0);

	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeam2Dlin(*this);
		return ec;
	}

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual const char* GetElementSpec() const {return "ANCFBeam2Dlin";}

	virtual int NS() const {return 4;}
	virtual int NodeSize() const {return 4;}

	virtual int SOS() const {return 8;}; //size of K and M
	virtual int SOSowned() const {return sos2;}; //len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;
	virtual void GetS0x(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0y(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0yx(Vector& sfx, const Vector2D& ploc) const;
	virtual void GetS0xx(Vector& sfx, const Vector2D& ploc) const;

protected:


};




#endif

