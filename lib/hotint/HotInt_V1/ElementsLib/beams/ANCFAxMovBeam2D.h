//#**************************************************************
//#
//# filename:             ANCFAxMovBeam2D.cpp
//#
//# author:               Astrid Sinwel, Gerstmayr Johannes
//#
//# generated:						August 4, 2009
//# description:          axially moving 2D ANCF - element
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
 
#ifndef ANCFAXMOVBEAM2D__H
#define ANCFAXMOVBEAM2D__H


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFAxMovBeam2D ANCFAxMovBeam2D ANCFAxMovBeam2D ANCFAxMovBeam2D ANCFAxMovBeam2D ANCFAxMovBeam2D ANCFAxMovBeam2D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// beam element including axial movement, see also similar element ANCFPipe2D
class ANCFAxMovBeam2D: public ANCFCable2D
{
public:

	ANCFAxMovBeam2D(MBS* mbsi):ANCFCable2D(mbsi),/* SV2(),*/
		velocity(0), useinitialcurvature(1), initcurv_endtime(0)
	{
		;
	};

	ANCFAxMovBeam2D(const ANCFAxMovBeam2D& e):ANCFCable2D(e.mbs)/*, SV2()*/

	{
		CopyFrom(e);
	};

	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)
	ANCFAxMovBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli):ANCFCable2D(mbsi), /*SV2(),*/
		velocity(0), useinitialcurvature(1), initcurv_endtime(0)
	{
		SetANCFCable2D(xc1, xc2, rhoi, Emi, si, coli);
	}

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFAxMovBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli):ANCFCable2D(mbsi),/* SV2(),*/
		velocity(0), useinitialcurvature(1), initcurv_endtime(0)
	{
		SetANCFCable2D(xc1, xc2, n1i, n2i, rhoi, Emi, si, coli);
	}

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFAxMovBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
		const Vector3D& si, const Vector3D& coli):ANCFCable2D(mbsi), /*SV2(),*/
		velocity(0), useinitialcurvature(1), initcurv_endtime(0)
	{
		SetANCFCable2D(xc1, xc2, n1i, n2i, materialnumi, si, coli);
	}

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFAxMovBeam2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2,
		int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli):ANCFCable2D(mbsi),/* SV2(),*/
		velocity(0), useinitialcurvature(1), initcurv_endtime(0)
	{
		SetANCFCable2D(xc1, xc2, vc1, vc2, n1i, n2i, rhoi, Emi, si, coli);
	}

	virtual ~ANCFAxMovBeam2D ()
	{
		delete velocity;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFAxMovBeam2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		ANCFCable2D::CopyFrom(e);
		const ANCFAxMovBeam2D& ce = (const ANCFAxMovBeam2D&)e;
		velocity = new MathFunction(*(ce.velocity));
		//SV2 = ce.SV2;
		useinitialcurvature = ce.useinitialcurvature;
		initcurv_endtime = ce.initcurv_endtime;
		order_eps = ce.order_eps;
		order_kappa = ce.order_kappa;
	}

	virtual void Initialize();

	virtual void EvalM(Matrix& m, double t);

	virtual void EvalF2(Vector& f, double t);

	virtual void StiffnessMatrix(Matrix& m) ;
	virtual int FastStiffnessMatrix() const;


	virtual void DrawElement();

	// axial velocity of the moving beam:
	virtual void SetAxialVelocity(const Vector& time_values, const Vector& vel_values)
	{
		delete velocity;
		velocity = new MathFunction();
		velocity->SetPiecewise(time_values, vel_values, 1);
	}
	virtual void SetAxialVelocity(double vel)
	{
		delete velocity;
		velocity = new MathFunction();
		velocity->SetConstant(vel);
	}

	virtual void SetAxialVelocity(const MathFunction& velfun)
	{
		delete velocity;
		velocity = new MathFunction(velfun);
	}

	virtual double AxialVelocity() const {return velocity->Evaluate(GetMBS()->GetTime());}
	virtual double AxialVelocity(double time) const {return velocity->Evaluate(time);}

	// initial curvature of elements is used,
	// curved element is free of stress
	virtual int& UseInitialCurvature() {return useinitialcurvature; }
	virtual const int& UseInitialCurvature() const {return useinitialcurvature; }

	virtual double& InitialCurvature_Endtime() {return initcurv_endtime;}
	virtual double GetEpsInit(const Vector2D& ploc, int flagD=0) const;


	virtual Vector2D GetGeneralizedVel(const Vector2D& ploc, Vector& xg, double time) const;
	virtual Vector2D GetGeneralizedVelD(const Vector2D& ploc, Vector& xg, double time) const;

	virtual void ComputeStressD(const Vector2D& ploc, int type, int c1, int c2, double& val)
	{
		ComputeStress(ploc, type, c1, c2, val, 1);
	}

	virtual void ComputeStress(const Vector2D& ploc, int type, int c1, int c2, double& val, bool flagD=0);

	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);

	virtual int GetNDof() const { return 8; }

	// get velocity vector v at local coordinate ploc
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	// weight for initial curvature: 0 means element is initially straight and stressless in straight config,
	//   1 means element is initially curved and stressless in curved config,
	//   initial curvature can be released over time
	virtual double GetInitLoadfact(int flagD=0) const
	{
		if (!UseInitialCurvature()) return 0;
		double loadfact = 1;
		if (initcurv_endtime > 0)
		{
			if (flagD) loadfact = (initcurv_endtime-GetMBS()->GetDrawTime()) / initcurv_endtime;
			else loadfact = (initcurv_endtime-GetMBS()->GetTime()) / initcurv_endtime;
			if (loadfact < 0) loadfact = 0;
		}
		return loadfact;
	}

protected:

	// axial velocity of the moving beam:
	MathFunction* velocity;
	// shape vector??
	//Vector SV2;
	// initial curvature of elements is used,
	// curved element is free of stress
	int useinitialcurvature;
	double initcurv_endtime;

	// these quantities are initialized Initialize()!
	int order_eps, order_kappa;
	// jacobi determinant at integration points for kappa and eps
	Vector det_intrule_kappa, det_intrule_eps;
	// initial values of kappa and eps
	Vector kappa0, eps0;
};



#endif
