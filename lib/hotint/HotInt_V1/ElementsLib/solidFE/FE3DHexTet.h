//#**************************************************************
//#
//# filename:             FE3DHexTet.h
//#
//# author:               Gerstmayr Johannes, YV
//#
//# generated:						October 2010
//# description:          3D-HexahedralGeneric finite element
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

#include "FiniteElement3D.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//HEXAHEDRAL  HEXAHEDRAL  HEXAHEDRAL  HEXAHEDRAL  HEXAHEDRAL  HEXAHEDRAL  HEXAHEDRAL  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<class FiniteElement3DBase>
class HexahedralGeneric : public FiniteElement3DBase
{
public:
	HexahedralGeneric(MBS* mbsi) : FiniteElement3DBase(mbsi) {};
	HexahedralGeneric(const HexahedralGeneric& e) : FiniteElement3DBase(e.mbs) { CopyFrom(e); }

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//depreciated functions without reference to material number:
	HexahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelememti=0);

	void SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelememti=0);

	void SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, const TArray<int>& nodelist,
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelememti=0);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	HexahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 int materialnumi, const Vector3D& coli, int CMSelememti=0);

	void SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc,
					 int material_num, const Vector3D& coli, int CMSelememti=0);

	void SetHexahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, const TArray<int>& nodelist,
					 int material_num, const Vector3D& coli, int CMSelememti=0);

	void SetHexahedral(int bodyindi, const TArray<int>& nodelist,
					 int material_num, const Vector3D& coli, int CMSelememti=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new HexahedralGeneric(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "HexahedralGeneric";}
	virtual TFiniteElementType GetElementType() const {return TFE_Hexahedral;}

	virtual double GetS0(const Vector3D& ploc, int shape) const; 
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const;
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const;

	virtual int NFaces() const {return 6;}
	virtual const FEFace & GetLocalFace(int i) const;
	//compute local face coordinates for drawing or for accessing face surface
	//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
	virtual void GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const;
	virtual Vector3D GetNodeLocPos(int i) const;

	
	//virtual int DataS() const  {return GetNIntegrationPoints(orderxyK)*6*IsInelasticMaterial();}

protected:
	virtual int GetActualInterpolationOrder() const;
	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
}; 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW TETRAHEDRAL  NEW TETRAHEDRAL  NEW TETRAHEDRAL  NEW TETRAHEDRAL  NEW TETRAHEDRAL
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<class FiniteElement3DBase>
class TetrahedralGeneric: public FiniteElement3DBase
{
public:
	TetrahedralGeneric(MBS* mbsi) : FiniteElement3DBase(mbsi) {};
	TetrahedralGeneric(const TetrahedralGeneric& e) : FiniteElement3DBase(e.mbs) {CopyFrom(e);}

	//nodal coordinates and initial velocities of points in 3D; 
	//order == 1 --> linear, order == 2 --> quadratic, etc.
	//CMSelement: 0, if no CMSelement, or otherwise the element number
	//depreciated:
	TetrahedralGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetTetrahedral(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetTetrahedral(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new TetrahedralGeneric(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "TetrahedralGeneric";}

	virtual TFiniteElementType GetElementType() const {return TFE_Tetrahedral;}

	//reference position (for estimates), -1..+1 based!!!
	virtual Vector3D GetRefPos() const {return GetPos(Vector3D(1./4.));};
	//reference position (for graphical issues), -1..+1 based!!!
	// virtual Vector3D GetRefPosD() const {return GetPosD(Vector3D(1./4.));};

	virtual double GetS0(const Vector3D& ploc, int shape) const; 
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const;
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const;

	virtual int NFaces() const {return 4;}
	virtual const FEFace & GetLocalFace(int i) const;

	//compute local face coordinates for drawing or for accessing face surface
	//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
	virtual void GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const;
	virtual Vector3D GetNodeLocPos(int i) const;

protected:
	virtual int GetActualInterpolationOrder() const;
	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
};

//$EK 2013-03-04: created class PrismGeneric
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW PRISM  NEW PRISM  NEW PRISM  NEW PRISM  NEW PRISM NEW PRISM NEW PRISM NEW PRISM
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<class FiniteElement3DBase>
class PrismGeneric: public FiniteElement3DBase
{
public:
	PrismGeneric(MBS* mbsi) : FiniteElement3DBase(mbsi) {};
	PrismGeneric(const PrismGeneric& e) : FiniteElement3DBase(e.mbs) {CopyFrom(e);}

	//nodal coordinates and initial velocities of points in 3D; 
	//order == 1 --> linear, order == 2 --> quadratic, etc.
	//CMSelement: 0, if no CMSelement, or otherwise the element number
	//depreciated:
	PrismGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetPrism(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetPrism(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PrismGeneric(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "PrismGeneric";}

	virtual TFiniteElementType GetElementType() const {return TFE_Prism;}

	//$EK ????????? what to use -> (1/3, 1/3, 1/2)
	//reference position (for estimates), -1..+1 based!!!
	virtual Vector3D GetRefPos() const {return GetPos(Vector3D(1./3., 1./3., 1./2.));};
	//reference position (for graphical issues), -1..+1 based!!!
	// virtual Vector3D GetRefPosD() const {return GetPosD(Vector3D(1./4.));};

	virtual double GetS0(const Vector3D& ploc, int shape) const; 
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const;
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const;

	virtual int NFaces() const {return 5;}
	virtual const FEFace & GetLocalFace(int i) const;

	//compute local face coordinates for drawing or for accessing face surface
	//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
	virtual void GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const;
	virtual Vector3D GetNodeLocPos(int i) const;

protected:
	virtual int GetActualInterpolationOrder() const;
	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//NEW PYRAMID  NEW PYRAMID  NEW PYRAMID  NEW PYRAMID  NEW PYRAMID NEW PYRAMID +++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<class FiniteElement3DBase>
class PyramidGeneric: public FiniteElement3DBase
{
public:
	PyramidGeneric(MBS* mbsi) : FiniteElement3DBase(mbsi) {};
	PyramidGeneric(const PyramidGeneric& e) : FiniteElement3DBase(e.mbs) {CopyFrom(e);}

	//nodal coordinates and initial velocities of points in 3D; 
	//order == 1 --> linear, order == 2 --> quadratic, etc.
	//CMSelement: 0, if no CMSelement, or otherwise the element number
	//depreciated:
	PyramidGeneric(MBS* mbsi, int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetPyramid(int bodyindi, const TArray<Vector3D>& xc, const TArray<Vector3D>& vc, 
					 double rhoi, double Emi, double nui, const Vector3D& coli, int CMSelementi=0);

	void SetPyramid(int bodyindi, const TArray<int>& nodelist, int material_num, const Vector3D& coli, int CMSelementi=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new PyramidGeneric(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "PyramidGeneric";}

	virtual TFiniteElementType GetElementType() const {return TFE_Pyramid; }

	
	//reference position (for estimates), -1..+1 based!!!
	// $EK center of mass!!!
	virtual Vector3D GetRefPos() const {return GetPos(Vector3D(3./8., 3./8., 1./4.));};
	//reference position (for graphical issues), -1..+1 based!!!
	// virtual Vector3D GetRefPosD() const {return GetPosD(Vector3D(1./4.));};

	virtual double GetS0(const Vector3D& ploc, int shape) const; 
	virtual void GetDSMatrix0(const Vector3D& ploc, Matrix& sm) const;
	//virtual double GetDS0(const Vector3D& ploc, int shape, int dxj) const;

	virtual int NFaces() const {return 5;}
	virtual const FEFace & GetLocalFace(int i) const;

	//compute local face coordinates for drawing or for accessing face surface
	//the vectors v1 and v2 are local vectors which point from refpos to the other nodes
	virtual void GetLocalFaceCoordinates(int face, Vector3D& refpos, Vector3D& v1, Vector3D& v2, Vector3D& v3) const;
	virtual Vector3D GetNodeLocPos(int i) const;

protected:
	virtual int GetActualInterpolationOrder() const;
	// implementation of IntegrationRuleProvider
	virtual void DefineIntegrationRule(IntegrationRule & integrationRule);
};

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// these are the actual classes to be used in models
typedef HexahedralGeneric<FiniteElement3D> Hexahedral;
typedef TetrahedralGeneric<FiniteElement3D> Tetrahedral;
typedef PrismGeneric<FiniteElement3D> Prism;
typedef PyramidGeneric<FiniteElement3D> Pyramid;

