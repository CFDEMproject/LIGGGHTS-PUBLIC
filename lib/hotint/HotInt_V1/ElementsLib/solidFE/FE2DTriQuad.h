//#**************************************************************
//#
//# filename:             FE2DTriQuad.h
//#
//# author:               Gerstmayr Johannes, YV
//#
//# generated:						October 2010
//# description:          2D - triangle and quadrilateral
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
 
#pragma once

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// QUAD  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// 4-node or 9-node quadrilateral
class Quadrilateral : public FiniteElement2D
{
public:
	Quadrilateral(MBS* mbsi) : FiniteElement2D(mbsi) {};
	Quadrilateral(const Quadrilateral& e) : FiniteElement2D(e.mbs) { CopyFrom(e); }

	void SetQuadrilateral(int bodyindi, const TArray<int>& nodelist, int material_num, double thickness, const Vector3D& coli);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Quadrilateral(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "Quadrilateral";}
	virtual TFiniteElementType GetElementType() const {return TFE_Quadrilateral;}

	virtual double GetS0(const Vector2D& ploc, int shape) const; 
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sm) const;
	virtual Vector2D GetNodeLocPos2D(int i) const;

protected:
	virtual int GetActualInterpolationOrder() const;

	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
};