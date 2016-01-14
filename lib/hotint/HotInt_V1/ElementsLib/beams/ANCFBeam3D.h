//#**************************************************************
//#
//# filename:             ANCFBeam3D.h
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
 
#ifndef ANCFBEAM3D__H 
#define ANCFBEAM3D__H

#include "body3d.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBEAM3D ANCFBEAM3D ANCFBEAM3D ANCFBEAM3D ANCFBEAM3D ANCFBEAM3D ANCFBEAM3D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


const int ANCFBeamMaxIP=45;//45;
//rigid cube
class ANCFBeam3D: public Body3D
{
public:
	//Body3D():Element() {mbs = NULL;};
	ANCFBeam3D(MBS* mbsi):Body3D(mbsi), massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T1(), T2(), T(), xg(), xgd(), e0() {};
	ANCFBeam3D(const ANCFBeam3D& e):Body3D(e.mbs),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), x3(), w1(), w2(), w3(),
			T1(), T2(), T(), xg(), xgd(), e0() {CopyFrom(e);};
	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)
	ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0);

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0);
	
	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFBeam3D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2, int n1i, int n2i, double rhoi, double Emi, double nui,
		const Vector3D& si, const Vector3D& coli, int beammodel = 0);

	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeam3D(*this);
		return ec;
	}
	virtual void CopyFrom(const Element& e)
	{
		Body3D::CopyFrom(e);
		const ANCFBeam3D& ce = (const ANCFBeam3D&)e;
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
		x1 = ce.x1; x2 = ce.x2; x3 = ce.x3; 
		w1 = ce.w1; w2 = ce.w2; w3 = ce.w3; 
		orderx = ce.orderx;
		orderyz = ce.orderyz;

		for (int i=0; i < ANCFBeamMaxIP; i++)
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

		beamEA = ce.beamEA;
		beamEIy = ce.beamEIy;
		beamEIz = ce.beamEIz;
		beamGJkx = ce.beamGJkx;
		beamGAky = ce.beamGAky;
		beamGAkz = ce.beamGAkz;

		orderWl = ce.orderWl;
		orderWs = ce.orderWs;
		orderWsHR = ce.orderWsHR;
		orderWt = ce.orderWt;
		orderWb = ce.orderWb;

		Wlepsx0 = ce.Wlepsx0 ;
		Wlepsy0 = ce.Wlepsy0 ;
		Wlepsz0 = ce.Wlepsz0 ;
		Wlgamyz0= ce.Wlgamyz0;
		WsHRk10 = ce.WsHRk10 ;
		WsHRk20 = ce.WsHRk20 ;
		Wtkapx0 = ce.Wtkapx0 ;
		Wbkapy0 = ce.Wbkapy0 ;
		Wbkapz0 = ce.Wbkapz0 ;

		factstiffWl = ce.factstiffWl;

		rho= ce.rho; //DR 2013-02-04 deleted rho from class element
	}

	virtual void Initialize() 
	{
		Body3D::Initialize();
	}
	virtual void LinkToElements();
	virtual void BuildDSMatrices(); 
	virtual void InitConstructor() //initialize while constructor is called
	{
		concentratedmass1 = 0;
		concentratedmass2 = 0;
	}

	virtual const char* GetElementSpec() const {return "ANCFBeam";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual int SOS() const {return 24;}; //size of K and M
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
	virtual Vector3D GetNodeLocPos(int i) const
	{
		if (i==1) return Vector3D(-0.5*lx,0.,0.);
		else if (i==2) return Vector3D(0.5*lx,0.,0.);
		
		return Vector3D(0.);
	}
	virtual int IsRigid() const {return 0;}
	virtual void SetConcentratedMass(double cm1, double cm2) {concentratedmass1=cm1; concentratedmass2=cm2;}

	virtual void SetInitialCurvatures(); //set initial vectors for curvature, stretch, shear deformation for Schwab-Meijaard model
	virtual void SetRectangularBeamParameters(); //set beam parameters for rectangular cross-section, Schwab-Meijaard model

	virtual void SetBeamParameters(double EA, double EIw, double EIv, double GJ, double GA);

	virtual void SetBeamParameters2(double beamEAi, double beamEIyi, double beamEIzi, double beamGJkxi, 
		double beamGAkyi, double beamGAkzi);

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+24));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+24));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+24));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+24));}

	virtual int NS() const {return 8;}

	//compute v = T*v
	virtual void ApplyT(Vector& v) const;
	virtual void ApplyTD(Vector& v) const;

	virtual void ApplyTtp(Vector& v) const;

	virtual void ApplyTtp(Matrix& m) const;

	virtual void GetS0(Vector& sf, const Vector3D& ploc) const;

	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sf) const;

	virtual void GetS0x(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetS0y(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetS0z(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetS0yx(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetS0zx(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetS0xx(Vector& sfx, const Vector3D& ploc) const;
	virtual void GetRot(Matrix3D& rot, const Vector& xg) const;

	virtual double GetTau(const Vector& xg) const;
	virtual void GetDeltaTau(const Vector& xg, Vector& deltatau) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const; //test only

	virtual void GetSMatrix(const Vector& sf, Matrix& sm) const
	{
		for (int j = 1; j <= 3; j++)
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

	virtual Vector3D GetVel(const Vector3D& p_loc) const;

	virtual Vector3D GetPosD(const Vector3D& p_loc) const;

	//in reference element coordinates (-1..1)
	virtual Vector3D GetPos0D(const Vector3D& p_loc, double def_scale = 1) const;

	virtual Vector3D GetVelD(const Vector3D& p_loc) const;

	virtual Vector3D GetPosx(const Vector3D& p_loc) const;

  virtual Vector3D GetDisplacement(const Vector3D& p_loc) const;
  virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const;

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
		for (int j = 1; j <= 3; j++)
		{
			for (int i = 1; i <= 3; i++)
			{ 
				jac(i,j) = 0;
				for (int k=1; k <= ns; k++)
				{ 
					jac(i,j) += DS(j,k)*x0((k-1)*3+i);
				}
			}
		}
		//global_uo << "jac=" << jac << "\n";
	}

	void GetDMatrix(Matrix& D, double nu, double em) const
	{
		//double nu = 1./2.*la/(la+mu);
		//double em = mu*(3.*la+2.*mu)/(la+mu);
		D(1,1)=em*(-1.+nu)/(1.+nu)/(-1.+2.*nu);D(1,2)=em*nu/(1.+nu)/(1-2*nu);D(1,3)=em*nu/(1.+nu)/(1-2*nu);D(1,4)=0;D(1,5)=0;D(1,6)=0;
		D(2,1)=em*nu/(1.+nu)/(1-2*nu);D(2,2)=em*(-1.+nu)/(1.+nu)/(-1.+2.*nu);D(2,3)=em*nu/(1.+nu)/(1-2*nu);D(2,4)=0;D(2,5)=0;D(2,6)=0;
		D(3,1)=em*nu/(1.+nu)/(1-2*nu);D(3,2)=em*nu/(1.+nu)/(1-2*nu);D(3,3)=em*(-1.+nu)/(1.+nu)/(-1.+2.*nu);D(3,4)=0;D(3,5)=0;D(3,6)=0;
		D(4,1)=0;D(4,2)=0;D(4,3)=0;D(4,4)=.5000000000*em/(1.+nu);D(4,5)=0;D(4,6)=0;
		D(5,1)=0;D(5,2)=0;D(5,3)=0;D(5,4)=0;D(5,5)=.5000000000*em/(1.+nu);D(5,6)=0;
		D(6,1)=0;D(6,2)=0;D(6,3)=0;D(6,4)=0;D(6,5)=0;D(6,6)=.5000000000*em/(1.+nu);
	}

	virtual void EvalF2(Vector& f, double t); 

	virtual int FastStiffnessMatrix() const; //1==Position derivatives done analytically, damping added with numerical differentiation
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components, m might be larger


	virtual Matrix3D GetRotMatrix(const Vector3D& ploc) const;
	virtual Matrix3D GetRotMatrixP(const Vector3D& ploc) const;
	virtual Matrix3D GetRotMatrixD(const Vector3D& ploc) const;

	virtual Matrix3D ComputeTangentFrame(double x, const Vector& xg) const;
	virtual Matrix3D ComputeTangentFrameP(double x, const Vector& xg, const Vector& xgp) const;
	virtual void GetTangentFramedRotvdqT(const Vector3D& vloc, double x, const Vector& xg, Matrix& d);

	virtual Matrix3D ComputeCrosssectionFrame(double x, const Vector& q) const;
	virtual Matrix3D ComputeCrosssectionFrameP(double x, const Vector& q, const Vector& qp) const;
	virtual void GetCrosssectionFramedRotvdqT(const Vector3D& vloc, double x, const Vector& q, Matrix& dAvdq);


	virtual Vector3D GetAngularVel(const Vector3D& p_loc) const 
	{
		Matrix3D omega_skew = GetRotMatrixP(p_loc)*GetRotMatrix(p_loc).GetTp();

		return Vector3D(omega_skew(2,3),omega_skew(1,3),-omega_skew(1,2));
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
		grad.SetSize(3,NS());
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
	virtual void ApplyDrotdq(Vector& f, const Vector3D& x)
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
	virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d);

	virtual void GetdPosdx(const Vector3D& ploc, Vector3D& dpdx) 
	{
		//optimierungspotenzial 500% !!!!!!!!!!!!!!!!!!!
		double p0 = ploc.X()/(0.5*lx);

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

	virtual void ComputeCorotationalFrame(Vector3D& pref, Matrix3D& Aref);
	virtual void ComputeCorotationalFrameD(Vector3D& pref, Matrix3D& Aref);
	virtual void ComputeCorotationalFrameDAvdq(Vector3D& v, Matrix& dAvdq);
	virtual void ComputeCorotationalFrameDATvdq(Vector3D& v, Matrix& dATvdq);

	//ploc -1 ... +1
	virtual void GraduD(const Vector3D& ploc, const Vector& u, Matrix3D& gradu) const;

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD);

	virtual void DrawElement();

protected:
	//mechanical:
	double lx, ly, lz, Em, nu;

	int n1, n2, sos2;
	double concentratedmass1, concentratedmass2;

	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix, Hmatrix;
	Matrix DS; //3*8
	Vector SV; //8
	Vector temp;
	int elasticforce_beam;

	//integration points
	Vector x1,x2,x3,w1,w2,w3;
	int orderx, orderyz;
	Matrix grad[ANCFBeamMaxIP]; //DS transformed, grad u = grad[i]*e
	double jacdet[ANCFBeamMaxIP];

	Matrix T1, T2, T; //slope discontinuities
	Matrix3D jac1, jac2; //slope discontinuities
	Vector e0; //initial vector in e, not in p

	Matrix K0; //stiffness matrix at initial configuration

	double beamEA; //for beam formulation
	double beamEIy;
	double beamEIz;
	double beamGJkx;
	double beamGAky;
	double beamGAkz;

	double factstiffWl; //reduce stiffness of thickness modes, standard==1

	int orderWl;   //5; 5 and 6 coincide with 8 up to 10 digits, 7=8
	int orderWs;   //2; more than 2 gives locking, 1 is rather inaccurate
	int orderWsHR; //4; for Hellinger Reissner, 4 is exact
	int orderWt;   //2; 1=2=4, 1 does not work in L-shape!!!
	int orderWb;   //3; 3=4=6, 2 gives shear locking


	Vector Wlepsx0;
	Vector Wlepsy0;
	Vector Wlepsz0;
	Vector Wlgamyz0;
	Vector WsHRk10;
	Vector WsHRk20;
	Vector Wtkapx0;
	Vector Wbkapy0;
	Vector Wbkapz0;

	double rho; //DR 2013-02-04 deleted rho from class element
};




#endif

