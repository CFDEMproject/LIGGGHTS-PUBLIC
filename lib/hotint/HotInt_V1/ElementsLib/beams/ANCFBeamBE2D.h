//#**************************************************************
//# filename:             ANCFBeamBE2D.h
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

const int ANCFBeamBE2DmaxDOF = 12;
const int ANCFBeamBE2DmaxIP = 10;

class ANCFBeamBE2D : public FiniteElementGenericBeam2D
{
public:
	ANCFBeamBE2D(MBS* mbs) : FiniteElementGenericBeam2D(mbs) {};
	ANCFBeamBE2D(const ANCFBeamBE2D& e) : FiniteElementGenericBeam2D(e.mbs) { CopyFrom(e); }
	~ANCFBeamBE2D(void);

	virtual Element* GetCopy()
	{
		Element* e = new ANCFBeamBE2D(*this);
		return e;
	};

	void SetANCFBeamBE2D(int bodyind,
		const Vector& xc1, 
		const Vector& xc2, 
		int n1, 
		int n2, 
		int material_num,
		const Vector3D& size,
		const Vector3D& color);

	virtual const char* GetElementSpec() const { return "ANCFBeamBE2D"; }
	virtual TFiniteElementType GetElementType() const { return TFE_ThinBeam2D; }
	virtual int GetActualInterpolationOrder() const { return 10; }

	virtual void DefineIntegrationRule(IntegrationRule& integrationRule);

	virtual void EvalF2(Vector& f, double t);
	virtual void EvalM(Matrix& m, double t);

	virtual void GetH(Matrix& H);
	virtual void GetIntDuDq(Matrix& dudq) {	GetH(dudq);	}

	virtual void GetdPosdqT(const Vector2D& p_loc, Matrix& dpdqi);

	virtual int NS() const { return 6; };
	virtual int DOFPerNode() const { return NS(); }
	virtual Vector3D GetDOFDirD(int idof) const;
	virtual Vector3D GetDOFPosD(int idof) const;

	virtual Vector2D GetInplaneUnitVector2D(const double& p_loc) const;
	virtual Vector2D GetInplaneUnitVector2DD(const double& p_loc) const;
	virtual Vector2D GetRefInplaneUnitVector2D(const double& p_loc) const;

	virtual Vector2D GetInplaneUnitVectorP2D(const double& p_loc) const;

	virtual double GetEpsAxial(const double& p_loc) const;
	virtual double GetEpsAxialD(const double& p_loc) const;
	virtual double GetKappa(const double& p_loc) const;
	virtual double GetKappaD(const double& p_loc) const;
	virtual void GetDeltaEpsAxial(const double& p_loc, Vector& delta_eps);
	virtual void GetDeltaKappa(const double& p_loc, Vector& delta_kappa);

#pragma region shapefunction routines
	virtual double GetS0(const Vector2D& p_loc, int shape) const;
	virtual double GetS0x(const Vector2D& p_loc, int shape) const;
	virtual double GetS0xx(const Vector2D& p_loc, int shape) const;
#pragma endregion
};
