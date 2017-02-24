//#**************************************************************
//#
//# filename:             Beam2DFFRF.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						April 2006
//# description:          2D-plane stress plate element for floating frame of reference or absolute coordinates
//# remarks:						  
//#
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
 
#ifndef Beam2DFFR__H
#define Beam2DFFR__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


const int Beam2DMaxIP = 10; //maximum number of integration point
//rigid cube
class Beam2D: public Body2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Beam2D(MBS* mbsi):Body2D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), w1(),
			xg(), xgd() {};
	Beam2D(const Beam2D& e):Body2D(e.mbs),massmatrix(), Hmatrix(), SV(), DS(), x1(), w1(),
			xg(), xgd() {CopyFrom(e);};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Beam2D(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const Beam2D& ce = (const Beam2D&)e;

		beaml = ce.beaml;
		beamh = ce.beamh;
		beamt = ce.beamt;
		EI = ce.EI;
		rhoA = ce.rhoA;
		nfm = ce.nfm;

		//size = ce.size;
		xg = ce.xg;
		xgd = ce.xgd;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		temp = ce.temp;
		SV = ce.SV;
		DS = ce.DS;

		//integration points
		x1 = ce.x1;
		w1 = ce.w1;
		orderxy = ce.orderxy;
		orderxyM = ce.orderxyM;
		orderxyH = ce.orderxyH;

		//ffrfmode = ce.ffrfmode;

	}

	virtual void Initialize() 
	{
		Body2D::Initialize();
	}

	virtual void BuildDSMatrices(); 

	virtual int SOS() const {return NS();}; //size of K and M
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size
	virtual int FlexDOF() const {return NS();}
	virtual int NFM() const {return nfm;}
	virtual int FFRFmode() const {return ffrfmode;}; //model of FFRF eigenmodes

	virtual int IsRigid() const {return 0;}
	//virtual void SetConcentratedMass(double cm, int i) {concentratedmass[i-1] = cm;}

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	//virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	//virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg(iloc+SOS()));}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//change for other shape functions and dimension:
	virtual int Dim() const {return 2;} //must be dimension of element (2=2D, 3=3D)
	virtual int NS() const 
		{
			//if (FFRFmode()==3) return 4+NFM();
			return 2+NFM();
		}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//floating frame of reference formulation: for FFRF elements
	virtual double GetS0(double ploc, int shape) const;
	virtual void GetS0(Vector& sf, double ploc) const;
  double GetDS0(double ploc, int shape) const;
	virtual void GetDSMatrix0(Vector& sf, double ploc) const;
	virtual void GetDDSMatrix0(Vector& sf, double ploc) const;
	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components, m might be larger


	virtual Vector2D GetPos2D(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	virtual Vector2D GetPos2Drel(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2Drel(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DrelD(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DrelD(const Vector2D& p_loc) const;

	virtual Box3D GetElementBox() const
	{
		Box3D b;
		b.Add(ToP3D(GetPos2D(Vector2D(-beaml*0.5,-beamh*0.5))));
		b.Add(ToP3D(GetPos2D(Vector2D( beaml*0.5,-beamh*0.5))));
		b.Add(ToP3D(GetPos2D(Vector2D(-beaml*0.5, beamh*0.5))));
		b.Add(ToP3D(GetPos2D(Vector2D( beaml*0.5, beamh*0.5))));
		return b;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b;
		b.Add(ToP3D(GetPos2DD(Vector2D(-beaml*0.5,-beamh*0.5))));
		b.Add(ToP3D(GetPos2DD(Vector2D( beaml*0.5,-beamh*0.5))));
		b.Add(ToP3D(GetPos2DD(Vector2D(-beaml*0.5, beamh*0.5))));
		b.Add(ToP3D(GetPos2DD(Vector2D( beaml*0.5, beamh*0.5))));
		return b;
	}


  virtual void GetJacobi(double jac, double ploc) const
	{
		jac = beaml*0.5;
	}

	virtual void SetComputeCoordinates()
	{
		xg.SetLen(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}
	virtual void SetDrawCoordinates()
	{
		xgd.SetLen(SOS());
		for (int i = 1; i <= SOS(); i++)
			xgd(i) = XGD(i);
	}

	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);


	virtual void EvalF2(Vector& f, double t); 


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

	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq) //in fact it is DuDq Transposed
	{
		GetH(dudq);
	}

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d)
	{
		//fill in, f.Length is already set
		UO() << "Not yet implemented\n";
	}

	virtual void DrawElement();

protected:

	//mechanical:
	int nfm;  //Number Free Modes
	int ffrfmode; //model of FFRF eigenmode
	double EI, rhoA, beaml, beamh, beamt;

	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix;
	Matrix Hmatrix;
	Matrix DS; //
	Vector SV; //Shape vector
	Vector temp;

	//integration points
	Vector x1,w1;
	int orderxy, orderxyM, orderxyH;

};


/////////////////////////////////////
//2-node beam element with variable number of internal shape functions,  no axial deformation:
/////////////////////////////////////
class Beam2DFFRF: public Beam2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Beam2DFFRF(MBS* mbsi):Beam2D(mbsi) {};
	Beam2DFFRF(const Beam2DFFRF& e):Beam2D(e.mbs), K(),	SbarS(), Ibar12S(), I1S() {CopyFrom(e);};


	//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
	//initial velocities!!!!!
	Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, double L0, double H0, double T0, 
												double rhoi, double EmIi, const Vector3D& coli, int isCMSi = 0);

	Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, int ffrfmodei, double L0, double H0, double T0, 
												double rhoi, double EmIi, const Vector3D& coli, int isCMSi = 0);

	Beam2DFFRF(MBS* mbsi, int FFRFindex, int nfmi, int ffrfmodei, double L0, double H0, double T0, 
		double H20, double T20,	double rhoi, double EmIi, const Vector3D& coli, int isCMSi = 0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Beam2DFFRF(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Beam2D::CopyFrom(e);
		const Beam2DFFRF& ce = (const Beam2DFFRF&)e;


		FFRFind = ce.FFRFind;
		K = ce.K;

		SbarS = ce.SbarS;
		Ibar12S = ce.Ibar12S;
		I1S = ce.I1S;
		I1122S = ce.I1122S;

		ffrfmode = ce.ffrfmode;

		rho = ce.rho; //DR 2013-02-04 deleted rho from class element

		//include all data elements here!
	}

	virtual void LinkToElements();
	virtual void Initialize();

	virtual int FFRFDim() const {return 3;} //position X, Y, and angle
	virtual int SOS() const {return (1-IsCMS())*(FlexDOF()+FFRFDim());}; //size of K and M 
	virtual int SOSowned() const {return (1-IsCMS())*FlexDOF();}

	virtual int Dim() const {return 2;}

	virtual void GetI1(Vector& I1); //Shabana p. 209-211
	virtual double GetIkl(int k, int l);
	virtual void GetIbarkl(int k, int l, Vector& I1);
	virtual void GetSbar(Vector& Sbar);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for CMS Element
	virtual int IsCMS() const {return IsType(TCMS);}; //is part of component mode synthesis?

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}

	virtual const double& GetXact(int i) const 
	{
		if (IsCMS()) return ReferenceFrame().GetXactFull(i);
		else return mbs->GetXact(i);
	}
	virtual double& GetXact(int i)
	{
		if (IsCMS()) return ReferenceFrame().GetXactFull(i);
		else return mbs->GetXact(i);
	}
	virtual const double& GetDrawValue(int iloc) const 
	{
		if (IsCMS()) return ReferenceFrame().GetDrawValueFull(iloc);
		else return mbs->GetDrawValue(iloc);
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual const ReferenceFrame2D& ReferenceFrame() const {return (ReferenceFrame2D&)(GetMBS()->GetElement(FFRFind));}
	virtual ReferenceFrame2D& ReferenceFrame() {return (ReferenceFrame2D&)(GetMBS()->GetElement(FFRFind));}

	virtual double GetAngle2D() const {return XG(FlexDOF()+3);};
	virtual Vector2D GetRefPos2D() const {return Vector2D(XG(FlexDOF()+1),XG(FlexDOF()+2));};
	virtual Vector2D GetRefPos2DD() const {return Vector2D(XGD(FlexDOF()+1), XGD(FlexDOF()+2));};
	virtual Vector2D GetRefVel2D() const {return Vector2D(XGP(FlexDOF()+1),XGP(FlexDOF()+2));};

	virtual Vector2D GetPos2D(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	//insert all 9 entries of mass matrix
	virtual void EvalM(Matrix& m, double t);

	virtual void EvalMff(Matrix& m, double t);

	//insert quadratic velocity vector
	virtual void EvalF2(Vector& f, double t); 

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

	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq); //in fact it is DuDq Transposed

	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d);

	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& d);
	virtual double GetAngle2D(const Vector2D& ploc) const;
	virtual double GetAngle2DP(const Vector2D& ploc) const;
	virtual double GetAngle2DD(const Vector2D& ploc) const;

	virtual double GetW(double ploc);
	virtual double GetWx(double ploc);
	virtual double GetWxx(double ploc);
	virtual double GetWxxP(double ploc);
	virtual double GetPower();

	double testval;

private:
	int FFRFind; //index to Element of reference frame
	double beamhi, beamti;
	Matrix K;

	///+++++++++++++++++++++
	//store matrices FFRF:
	Vector SbarS;
	Vector Ibar12S;
	Vector I1S;
	double I1122S;

	double rho; //DR 2013-02-04 deleted rho from class element
};

#endif

