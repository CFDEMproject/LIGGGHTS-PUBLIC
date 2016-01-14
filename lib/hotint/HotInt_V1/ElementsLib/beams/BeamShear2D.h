//#**************************************************************
//# filename:             BeamShear2D.h
//#
//# author:               Gruber, Nachbagauer
//#
//# generated:						
//# description:          
//# comments:
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
 

#pragma once
#include "FiniteElementGenericBeam2D.h"

const int BeamShear2DmaxDOF = 6;
const int BeamShear2DmaxIP = 2;

class BeamShear2D :	public FiniteElementGenericBeam2D
{
public:
	BeamShear2D(MBS* mbs) : FiniteElementGenericBeam2D(mbs) {};
	BeamShear2D(const BeamShear2D& e) : FiniteElementGenericBeam2D(e.mbs) { CopyFrom(e); }

	~BeamShear2D(void);

	virtual Element* GetCopy()
	{
		Element* e = new BeamShear2D(*this);
		return e;
	};

	void SetBeamShear2D(int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		int n1, 
		int n2, 
		int material_num,
		const Vector3D& size,
		const Vector3D& color);

	virtual const char* GetElementSpec() const { return "BeamShear2D"; }
	virtual TFiniteElementType GetElementType() const { return TFE_Beam2D; }
	virtual int GetActualInterpolationOrder() const { return 1; }
	virtual void DefineIntegrationRule(IntegrationRule& integrationRule);

	virtual void LinkToElements();

	virtual void EvalF2(Vector& f, double t);
	virtual void EvalM(Matrix& m, double t);

	/*virtual void GetH(Matrix& H);
	virtual void GetIntDuDq(Matrix& dudq) {	GetH(dudq);	}*/

	//virtual void GetdPosdqT(const Vector2D& p_loc, Matrix& dpdqi);

	virtual int NS() const { return 2; };
	virtual int DOFPerNode() const { return 3; }
	virtual Vector3D GetDOFDirD(int idof) const;
	virtual Vector3D GetDOFPosD(int idof) const;

	virtual Vector2D GetInplaneUnitVector2D(const double& p_loc) const { return GetDirector2(p_loc); }
	virtual Vector2D GetInplaneUnitVector2DD(const double& p_loc) const { return GetDirector2D(p_loc); }
	virtual Vector2D GetRefInplaneUnitVector2D(const double& p_loc) const { return GetRefDirector2(p_loc); }

	virtual Vector2D GetInplaneUnitVectorP2D(const double& p_loc) const { assert(0); return Vector2D(0.,0.); };

	virtual Vector2D GetDirector1(const double& p_loc) const;
	virtual Vector2D GetDirector1D(const double& p_loc) const;
	virtual Vector2D GetDirector2(const double& p_loc) const;
	virtual Vector2D GetDirector2D(const double& p_loc) const;
	virtual Vector2D GetRefDirector2(const double& p_loc) const;

	virtual double GetEpsAxial(const double& p_loc) const;
	virtual double GetEpsAxialD(const double& p_loc) const;
	virtual double GetTheta(const double& p_loc) const;
	virtual double GetThetaD(const double& p_loc) const;
	virtual double GetRefTheta(const double& p_loc) const;
	virtual double GetKappa(const double& p_loc) const;
	virtual double GetKappaD(const double& p_loc) const;
	virtual double GetGamma(const double& p_loc) const;
	virtual double GetGammaD(const double& p_loc) const;

	virtual void GetDeltaEpsAxial(const double& p_loc, Vector& delta_eps);
	virtual void GetDeltaKappa(const double& p_loc, Vector& delta_kappa);
	virtual void GetDeltaGamma(const double& p_loc, Vector& delta_gamma);

#pragma region shapefunction routines
	virtual double GetS0(const Vector2D& p_loc, int shape) const;
	virtual double GetS0x(const Vector2D& p_loc, int shape) const;
	virtual double GetS0xx(const Vector2D& p_loc, int shape) const;
#pragma endregion
};

