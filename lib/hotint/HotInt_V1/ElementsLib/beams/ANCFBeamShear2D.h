//#**************************************************************
//#
//# filename:             ANCFBeamShear2D.h
//#
//# author:               Astrid und Karin
//#
//# generated:						2009/2010
//# description:          2D ANCF beam element with shear deformation
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
 
#ifndef ANCFBeamShear2D__H
#define ANCFBeamShear2D__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFBeamShear2D ANCFBeamShear2D ANCFBeamShear2D ANCFBeamShear2D ANCFBeamShear2D 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



class ANCFBeamShear2D: public Body2D
{
public:
	//Body3D():Element() {mbs = NULL;};
	ANCFBeamShear2D(MBS* mbsi):Body2D(mbsi), massmatrix(), Hmatrix(), SV(), x1(), w1(),
		xg(), xgd(), q0() 
	{
		ks=1;		 //set to the default value 1, if there exists no value (means no shear correction)
		nnodes=2;//Default-Value (linear computation), in case nnodes is not set to a specific value
		RI=0;    //Default-Value: no reduced integration order
		BeamFormulation=0;
	}
	ANCFBeamShear2D(const ANCFBeamShear2D& e):Body2D(e.mbs),massmatrix(), Hmatrix(), SV(), x1(), w1(),
		xg(), xgd(), q0() 
	{
		CopyFrom(e);
	};

	//ANCFBeamShear2D(MBS* mbsi);
	void SetANCFBeamShear2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
		const Vector3D& si, const Vector3D& coli);
  //set-Fct. for xc3, n3i (for quadr. element):
	void SetANCFBeamShear2D(const Vector& xc1, const Vector& xc2, const Vector& xc3, int n1i, int n2i, int n3i, int materialnumi,
		const Vector3D& si, const Vector3D& coli);


	virtual Element* GetCopy()
	{
		Element* ec = new ANCFBeamShear2D(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)  //new variable -> copy it here
	{
		Body2D::CopyFrom(e);
		const ANCFBeamShear2D& ce = (const ANCFBeamShear2D&)e;
		lx = ce.lx;
		ly = ce.ly;
		lz = ce.lz;
		//Em = ce.Em;
		xg = ce.xg;
		xgd = ce.xgd;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		SV = ce.SV;

		//integration points
		x1 = ce.x1;
		w1 = ce.w1;
		order_Poisson = ce.order_Poisson;

		q0 = ce.q0;  //initial values

		nnodes=ce.nnodes;//nnodes=2 or =3
		sos2 = ce.sos2;
		n1 = ce.n1; n2 = ce.n2; n3=ce.n3;//add n3
		materialnum = ce.materialnum;
		ks = ce.ks;  //haben wir hinzugefügt: neue Variable -> im CopyFrom so eintragen (im .h oder .cpp-File)
								 //Konstruktor nullsetzen oder sinnvollen Wert, z.B. ks=0; (siehe weiter oben)
		RI=ce.RI;    //for Reduced Integration
		BeamFormulation=ce.BeamFormulation;//different formulations
	}

	virtual void Initialize()
	{
		Body2D::Initialize();
	}
	virtual const char* GetElementSpec() const {return "ANCFBeamShear2D";}

	virtual void LinkToElements();
	virtual void BuildDSMatrices();

	//nnodes =2 or =3
	//NS()=2*nnodes
	//SOS()=2*NS()

	virtual int SOS() const {return 2*NS();};//= 8 or 12  //size of K and M  
	virtual int SOSowned() const {return sos2;}; //len(u)
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size

	// one internal variable plasticstrain for perfect plasticity, two internal variables plasticstrain and effective plastic strain for hardening materials
	virtual int DataS() const {return 0;}  

	virtual int NNodes() const {if (SOSowned() == 0) return nnodes; else return 0;};//Get-Fkt, we need a Set-function
	virtual void SetNNodes(int nnodes_i)  //according Set-Fkt
	{
		nnodes = nnodes_i;
	}

	//n3 added in loop
	virtual const int& NodeNum(int i) const //this int stays const
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else if(i == 3) return n3;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear2D::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}
	virtual int& NodeNum(int i) //this int can be changed by user
	{
		if (i == 1) return n1;
		else if(i == 2) return n2;
		else if(i == 3) return n3;
		else 
		{
			mbs->UO() << "Error in ANCFBeamShear2D::NodeNum: node " << i << " does not exist!";
			return n1;
		}
	}

	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));}

	virtual void SetMaterialNum(int matnr)
	{
		materialnum = matnr;
	}


	void SetBeamParameters(double beamEIi, double beamEAi, double beamRhoAi)
	{
		GetMaterial().BeamEIy() = beamEIi;
		GetMaterial().BeamEA() = beamEAi;
		GetMaterial().BeamRhoA() = beamRhoAi;
	}
	int IsBeamParameters() const {if (GetBeamEIy() < 0) return 0; else return 1;}

	virtual int IsRigid() const {return 0;}

	virtual Box3D GetElementBox() const
	{
		return Box3D(ToP3D(Vector2D(XG(1),XG(2))),ToP3D(Vector2D(XG(5),XG(6))));
	}

	virtual Box3D GetElementBoxD() const
	{
		return Box3D(ToP3D(Vector2D(XGD(1),XGD(2))),ToP3D(Vector2D(XGD(5),XGD(6))));
	}

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}

	virtual int NS() const {return 2*nnodes;}//= 4 or 6

	//Shapefunctions
	virtual void GetShapes(Vector& sf, const Vector2D& ploc) const;
	virtual double GetSF(int i, const Vector2D& ploc) const;
	virtual void GetShapesx(Vector& sfx, const Vector2D& ploc) const;
	virtual double GetSFx(int i, const Vector2D& ploc) const;
	virtual void GetShapesxx(Vector& sfx, const Vector2D& ploc) const;
	virtual double GetSFxx(int i, const Vector2D& ploc) const;
	virtual void GetShapesy(Vector& sfy, const Vector2D& ploc) const;
	virtual double GetSFy(int i, const Vector2D& ploc) const;
	virtual void GetShapesxy(Vector& sfy, const Vector2D& ploc) const;
	virtual double GetSFxy(int i, const Vector2D& ploc) const;

	virtual void Gradu(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const;

	//Dimensions
	virtual double GetLx() const {return lx;}
	virtual double GetLy() const {return ly;}
	virtual double GetLz() const {return lz;}
	//for shear correction factor:
	virtual double Getks() const {return ks;}
	virtual void Setks(double ks_fact)  //Get-Set-Fct
	{
		ks = ks_fact;
	}
	//for reduced-integration-parameter:
	virtual double GetRI() const {return RI;}
	virtual void SetRI(int RI_fact)  //Get-Set-Fct
	{
		RI = RI_fact;
	}
	virtual void SetFormulation(int Form_fact) //Get-Set-Fct
	{
		BeamFormulation = Form_fact;
	}

	//Velocity, Displacement
	virtual Vector2D GetVel2D(const Vector2D& p_loc, int flagD=0) const;
	virtual Vector2D GetDisplacement2DD(const Vector2D& p_loc) const;

	//r,rx,rxx,rxy,ry (flagD for Visualization)
	virtual Vector2D GetPos2D(const Vector2D& p_loc, const Vector& xg) const;
	virtual Vector2D GetPos2D(const Vector2D& p_loc, int flagD) const;
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const { return GetPos2D(p_loc, 1); }
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const { return GetPos2D(p_loc, 0); }
	virtual Vector2D GetPosx2D(const Vector2D& p_loc, const Vector& xg) const;
	virtual Vector2D GetPosx2D(const Vector2D& p_loc, int flagD=0) const;
	virtual Vector2D GetPosy2D(const Vector2D& p_loc, const Vector& xg) const;
	virtual Vector2D GetPosy2D(const Vector2D& p_loc, int flagD=0) const;
	//virtual Vector2D GetPosxx2D(const Vector2D& p_loc, const Vector& xg) const;
	//virtual Vector2D GetPosxx2D(const Vector2D& p_loc, int flagD=0) const;
	virtual Vector2D GetPosxy2D(const Vector2D& p_loc, const Vector& xg) const;
	virtual Vector2D GetPosxy2D(const Vector2D& p_loc, int flagD=0) const;

	//r0 in reference element
	virtual Vector2D GetInitPos2D(const Vector2D& p_loc) const;

	//Plotscaling: l/2, h/2
	virtual Vector2D GetPos2D0D(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2D0D(const Vector2D& p_loc, double defscale) const;

	//t2(direction of the cross section), t1(perpendicular to the direction of cross section)
    virtual Vector2D GetT12D(const double& xbar) const;
	virtual Vector2D GetT22D(const double& xbar) const;

	//Gamma1, Gamma2 (and corresponding delta)
	virtual double GetGamma12D(const double& xbar) const;
	virtual double GetGamma22D(const double& xbar) const;
	virtual void GetDeltaGamma22D(const double& xbar, Vector& DeltaGamma2) const;
	virtual void GetDeltaGamma12D(const double& xbar, Vector& DeltaGamma1) const;

	//Angle: Theta' (and corresponding delta)
	virtual double GetThetax2D(const double& xbar) const;
	virtual void GetDeltaThetax2D(const double& xbar, Vector& DeltaThetax) const; //Vector.length=Number of DOFs

	//Eyy for thickness deformation (and corresponding delta)
	virtual double GetEyy2D(const double& xbar) const;
	virtual void GetDeltaEyy2D(const double& xbar, Vector& DeltaEyy) const;

//***************************************************************************************************************
//**********************************Änderungen bis hier**********************************************************
//***************************************************************************************************************

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


	virtual Vector2D GetRefPos2DD() const {return GetPos2D(Vector2D(0.,0.),1);};//

	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		if (idof <= 4) //left node: 1,...,4
		{
			return ToP3D(GetPos2D(Vector2D(-0.5*lx,0.),1));//1 für flagD
		}
		else if (idof <= 8) //right node: 5,...,8
		{
			return ToP3D(GetPos2D(Vector2D(0.5*lx,0.),1));
		}
		else if (idof <= 12) //middle node: 9,...,12
		{
			return ToP3D(GetPos2D(Vector2D(0.,0.),1));//middle node: (0,0)
		}
		else
		{
			return Vector3D(0, 0, 0);
		}
	}

	//Visualization: display the load-vector in the correct direction
	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		//Modulus-Calculus:
		if (idof>4) idof -= 4;//for 8
		if (idof>4) idof -= 4;//for 12

		if (idof == 1) return Vector3D(1.,0.,0.);
		else if (idof == 2) return Vector3D(0.,1.,0.);
		else if (idof == 3) return Vector3D(1.,0.,0.);
		else if (idof == 4) return ToP3D(Vector3D(0.,1.,0.));
		return Vector3D(0.,0.,0.);
	}


	virtual void SetComputeCoordinates()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);
	}
	virtual void SetComputeCoordinatesP()
	{
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGP(i);
	}

	virtual double GetBeamEIy() const {return GetMaterial().BeamEIy(); }
	virtual double GetBeamEA() const {return GetMaterial().BeamEA(); }
	virtual double GetRhoA() const {return GetMaterial().BeamRhoA(); }


	// Achtung
	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d);

	virtual void GetJacobi(Matrix3D& jac, const Vector2D& ploc, const Vector& xg0) const
	{
		jac.SetSize(2,2);
		jac.SetAll(0.);
		int ns = NS();
		for (int i = 1; i <= Dim(); i++)
		{ 
			for (int k=1; k <= ns; k++)
			{ 
				jac(i,1) += GetSFx(k, ploc)*xg0((k-1)*Dim()+i);
				jac(i,2) += GetSFy(k, ploc)*xg0((k-1)*Dim()+i);
			}
		}

		//mbs->UO() << "jac=" << jac << "\n";
	}

	virtual void GetH(Matrix& H);
	//Massmatrix
	virtual void EvalM(Matrix& m, double t);


	virtual void EvalF2(Vector& f, double t);


	virtual double PostNewtonStep(double t);		// changes plastic strains, returns yieldfunction(sigma)/Em
	virtual void PostprocessingStep();

	//for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq)
	{
		GetH(dudq);
	}

	// ploc is based [-lx/2, lx/2]
	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& d)
	{
		//{UO() << "called GetdAngle2DdqT";};

		d.SetSize(NS()*Dim(),1);
		d.FillWithZeros();
		
		Vector2D ry;
		ConstVector<12> Sy(NS()*Dim());
		ConstVector<12> xg(NS()*Dim());

		for (int i = 1; i <= SOS(); i++)
		{
			xg(i) = XG(i);
		}

		for (int i = 1; i <= NS(); i++)
		{
			Sy(i) = GetSFy(i,ploc);
		}
		ry=GetPosy2D(ploc,xg);

		double dyx = Sqr(ry.X())+Sqr(ry.Y());

		for (int i = 1; i <= NS(); i++)
		{
			// d theta/d q = -ry.y*Sy.x + ry.x*Sy.y / (ry.x^2+ry.y^2)
			d((i-1)*Dim()+1,1) =-(ry.Y()*Sy(i))/dyx;  // Sy.x
			d((i-1)*Dim()+2,1) = (ry.X()*Sy(i))/dyx;  // Sy.y
		}
	}
	//virtual void GetdAngle2Ddx(const Vector2D& ploc, double& dphidx);

	virtual double GetAngle2D(const Vector2D& ploc) const;
	//virtual double GetAngle2DP(const Vector2D& ploc) const;

	//virtual double GetCurvature2D(const Vector2D& p_loc) const;
	//virtual double GetCurvature2DP(const Vector2D& p_loc) const;
	//virtual double GetCurvature2DD(const Vector2D& p_loc) const {return 0;}


	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);
	virtual void DrawElement();

	//virtual double EA()
	//{if (beamEI != -1) return beamEA; else return Em*ly*lz; }


protected:
	//mechanical:
	double lx, ly, lz;
	double ks;  //shear correction factor
	int RI;
	int BeamFormulation;
	// into material:
	//double Em;
	//double beamEA, beamEI, beamRhoA;

	int n1, n2, n3, sos2;  //if nnodes=2, n3 is set to zero
  int nnodes;

	//temporary storage for acceleration:
	Vector xg, xgd; // [e]; e,x
	Matrix massmatrix, Hmatrix; //M = int(rho*((S)^T).S, dV,V); H = int(S,dV,V)
	Vector SV; 
	Vector temp;

	//integration points
	Vector x1,w1; //memory for numerical integration
	int order_Poisson;   //order of num. Int. in x-direction
	int order_Thickness;

	Vector q0; //initial vector for positions
};


#endif
