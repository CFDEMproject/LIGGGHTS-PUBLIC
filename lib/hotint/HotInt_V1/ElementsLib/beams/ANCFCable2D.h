//#**************************************************************
//#
//# filename:             ANCFCable2D.h
//#
//#
//# author:               Gerstmayr Johannes & Rafael Ludwig
//#
//# generated:						24.April 2006
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
 
#ifndef ANCFCABLE2D__H
#define ANCFCABLE2D__H

#include "body2d.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ANCFCable2D ANCFCable2D ANCFCable2D ANCFCable2D ANCFCable2D ANCFCable2D ANCFCable2D
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//rigid cube
class ANCFCable2D: public Body2D //$EDC$[beginclass,classname=ANCFCable2D,parentclassname=Body2D,addelementtype=TAEBody+TAE2D+TAENotInRelease,addelementtypename=ANCFCable2D,texdescription="
//The ANCFCable2D is a beam based on the Absolute Nodal Coordinate Formulation (ANCF). 
//For details of this beam, see the literature J. Gerstmayr, H. Irschik. On the correct representation of bending and axial deformation in the absolute nodal coordinate formulation with an elastic line approach, Journal of Sound and Vibration, Vol. 318, pp. 461-487, 2008. DOI:10.1016/j.jsv.2008.04.019
//"]
{
public:

	enum {ANCFCable2D_MaxIP = 10};
	enum {ANCFCable2D_MaxDOF = 12};
	enum {ANCFCable2D_MaxNNodes = 2};
	enum {ANCFCable2D_MaxNShapes = 6};
	enum {ANCFCable2D_maxppy = 50}; // maximal number of plastic grid points in y direction
	// -------------------------------------------
	// CONSTRUCTOR, COPY-ROUTINES...

	// Constructor
	ANCFCable2D(MBS* mbsi):Body2D(mbsi), massmatrix(), Hmatrix(),
		//SV(), xg(), xgd(), 
		e0(), epsdamping(0), use_tangent_stiffness_matrix(0)
	{
		ElementDefaultConstructorInitialization();
	};

	// Copy-Constructor
	ANCFCable2D(const ANCFCable2D& e):Body2D(e.mbs),massmatrix(), Hmatrix(), 
		//SV(), xg(), xgd(), 
		e0() , epsdamping(0), use_tangent_stiffness_matrix(0)
	{
		CopyFrom(e);
	};

	// Different constructors with nodal/geometric information, initial values...

	////nodal coordinates of first and second point, no global nodes given
	//ANCFCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
	//	const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	//// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	//// a new material with rhoi and Emi is added
	//ANCFCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
	//	const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	//// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	//// material is given by number
	//ANCFCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
	//	const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	//// Element shares nodes n1 and n2 with other elements; element sets initial conditions and initial velocities for nodes
	//// a new material with rhoi and Emi is added
	//ANCFCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2,
	//	int n1i, int n2i, double rhoi, double Emi,const Vector3D& si, const Vector3D& coli, int kappanodalI=0);



	//nodal coordinates of first and second point, no global nodes given
	void SetANCFCable2D(const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	// a new material with rhoi and Emi is added
	void SetANCFCable2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions for nodes
	// material is given by number
	void SetANCFCable2D(const Vector& xc1, const Vector& xc2, int n1i, int n2i, int materialnumi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	// Element shares nodes n1 and n2 with other elements; element sets initial conditions and initial velocities for nodes
	// a new material with rhoi and Emi is added
	void SetANCFCable2D(const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2,
		int n1i, int n2i, double rhoi, double Emi,const Vector3D& si, const Vector3D& coli, int kappanodalI=0);

	//this function assigns default values to the element variables
	virtual void ElementDefaultConstructorInitialization()
	{
		plasticstrains_height = 0;
		plasticstrains_width = 0;
		kappanodal = 0;
		n1=0; n2=0; sos2=0;
		//lx = 0; ly = 0; lz = 0; //$JG2012-02: removed, because not needed! x,ly,lz is contained in size!
		size = Vector3D(0.,0.,0.);
		mass = 0;
		//rho = 0; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
		col = Vector3D(0.2,0.2,0.8);
		concentratedmass1 = 0;
		concentratedmass2 = 0;

		x_init = Vector(2*SOS()); //zero initialized
		elementname = GetElementSpec();
	}

	// returns copy of element
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFCable2D(*this);
		return ec;
	}

	// copy routine
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const ANCFCable2D& ce = (const ANCFCable2D&)e;
		//lx = ce.lx; //$JG2012-02: removed, because not needed! lx,ly,lz is contained in size!
		//ly = ce.ly;
		//lz = ce.lz;
		//xg = ce.xg;
		//xgd = ce.xgd;
		//massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		//SV = ce.SV;

		//orderx = ce.orderx;

		e0 = ce.e0;

		sos2 = ce.sos2;
		n1 = ce.n1; n2 = ce.n2;
		concentratedmass1 = ce.concentratedmass1;
		concentratedmass2 = ce.concentratedmass2;

		kappanodal = ce.kappanodal;

		plasticstrains_height = ce.plasticstrains_height;
		plasticstrains_width = ce.plasticstrains_width;

		epsdamping = ce.epsdamping;

		// arc-length method
		numConstraintElem = ce.numConstraintElem;
		gcload = ce.gcload;
		loaddof = ce.loaddof;

		use_tangent_stiffness_matrix = ce.use_tangent_stiffness_matrix;
	}

	// ------------------------------------------------
	// INITIALIZING

	// initializing routine: is called before time integration
	// internal data vector is initialized
	virtual void Initialize()
	{
		Body2D::Initialize();
		if (IsMaterialAvailable() && GetMaterial().IsInelasticMaterial())
		{
			Vector datainit(DataS()); // this vector is automatically initialized with zeros
			SetDataInit(datainit);    //initialize data variables with zero = initial inelastic strain == zero
		}
	}

	// build local-to-global map for element
	virtual void LinkToElements();
	// build D-Shape function matrices, NOT IMPLEMENTED!!!
	virtual void BuildDSMatrices();

	virtual const char* GetElementSpec() const {return "ANCFCable2D";}
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	// ----------------------------------------------
	// NUMBERS OF DEGREES OF FREEDOM,...

	// second order size, size of K and M
	virtual int SOS() const {return 8+4*kappanodal;};
	// number of elements own degrees of freedom
	virtual int SOSowned() const {return sos2;}; 
	 //size of first order explicit equations
	virtual int ES() const  {return 0;}; 
	//implicit (algebraic) size
	virtual int IS() const  {return 0;};  

	virtual int IsMaterialAvailable() const {return (GetMaterialNum() != 0 && GetMaterialNum() <= GetMBS()->NMaterials());} //check if material is available (e.g. if element is initialized)
	// Internal DATA VARIABLES:
	// * one internal variable plasticstrain for perfect plasticity
	//   OR
	// * two internal variables plasticstrain and effective plastic strain for hardening materials
	virtual int DataS() const {if (!IsMaterialAvailable()) {return 0;} else {return int(GetMaterial().IsHardeningMaterial()+1)*plasticstrains_height*plasticstrains_width; }}  

	// number of nodes of the element
	virtual int NNodes() const {if (SOSowned() == 0) return 2; else return 0;};

	// number of shape functions
	virtual int NS() const {return 4+2*kappanodal;}

	// global node number
	virtual const int& NodeNum(int i) const 
	{
		if (i == 1) return n1; else return n2;
	}
	virtual int& NodeNum(int i) 
	{
		if (i == 1) return n1; else return n2;
	}

	// element is not rigid body
	virtual int IsRigid() const {return 0;}

	// space dimension
	virtual int Dim() const {return 2;}


	// --------------------------------------------------------
	// Computational routines:

	// int du/dq dV
	virtual void GetH(Matrix& H);
	virtual void GetIntDuDq(Matrix& dudq)
	{
		GetH(dudq);
	}
	// mass matrix
	virtual void EvalM(Matrix& m, double t);
	// residual
	virtual void EvalF2(Vector& f, double t);

	// -------------------------------------------------------
	// Nonlinear iteration for plasticity
	virtual double PostNewtonStep(double t);		// changes plastic strains, returns yieldfunction(sigma)/Em
	virtual void PostprocessingStep();

	// ------------------------------------------------------
	// VISUALIZATION

	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);

	virtual void DrawElement();

		// element bounding box 
	virtual Box3D GetElementBox() const
	{
		return Box3D(ToP3D(Vector2D(XG(1),XG(2))),ToP3D(Vector2D(XG(5+2*kappanodal),XG(6+2*kappanodal))));
	}

	virtual Box3D GetElementBoxD() const
	{
		return Box3D(ToP3D(Vector2D(XGD(1),XGD(2))),ToP3D(Vector2D(XGD(5+2*kappanodal),XGD(6+2*kappanodal))));
	}


	// ------------------------------------------------------
	// GET-SET ROUTINES AND SIMILAR

	// Set beam material parameters
	void SetBeamParameters(double beamEIi, double beamEAi, double beamRhoAi)
	{
		if (IsMaterialAvailable())
		{
			GetMaterial().BeamEIy() = beamEIi;
			GetMaterial().BeamEA() = beamEAi;
			GetMaterial().BeamRhoA() = beamRhoAi;
		}
		else
		{
			GetMBS()->UO() << "ERROR: ANCFCable2D::SetBeamParameters: no material available\n";
		}
	}
	// check if beam material parameters are to be used
	int IsBeamParameters() const {if (GetBeamEIy() < 0) return 0; else return 1;}

	// Flag if kappa is nodal coordinate
	virtual int IsKappaNodal() const {return kappanodal;}

	// set concentrated masses in nodal points
	virtual void SetConcentratedMass(double cm1, double cm2) {concentratedmass1=cm1; concentratedmass2=cm2;}

	// Element dimensions
	virtual double GetLx() const {return size.X();}
	virtual double GetLy() const {return size.Y();}
	virtual double GetLz() const {return size.Z();}

	// routine for storing element coordinates into a vector
	virtual void GetCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XG(i);
	}

	// storing element velocity coordinates into a vector
	virtual void GetCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGP(i);
	}

	// storing element draw coordinates into a vector
	virtual void GetDrawCoordinates(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGD(i);
	}

	// storing element draw velocities into a vector
	virtual void GetDrawCoordinatesP(Vector& dc) const
	{
		for (int i = 1; i <= SOS(); i++)
			dc(i) = XGPD(i);
	}

	//// routines for setting local vector xg to element coordinates XG
	//virtual void SetComputeCoordinates()
	//{
	//	for (int i = 1; i <= SOS(); i++)
	//		xg(i) = XG(i);
	//}
	//// routines for setting local vector xg to element velocity coordinates XGP
	//virtual void SetComputeCoordinatesP()
	//{
	//	for (int i = 1; i <= SOS(); i++)
	//		xg(i) = XGP(i);
	//}




	// Shape function routines ------------------------------------
	virtual void GetS0(Vector& sf, const double& ploc) const;
	virtual void GetS0x(Vector& sfx, const double& ploc) const;
	virtual void GetS0xx(Vector& sfx, const double& ploc) const;
	// get Shape Function i
	virtual double GetS0(int i, const double& ploc) const;
	virtual double GetS0x(int i, const double& ploc) const;
	virtual double GetS0xx(int i, const double& ploc) const;


	// -----------------------------------------------
	// Routines for computing positions/velocities...

	// compute derivative of position vector with respect to generalized coordinates at local position ploc
	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d)
	{
		//p = S(p.x,p.y,p.z)*q; d(p)/dq
		ConstVector<ANCFCable2D_MaxNShapes> SV(NS());
		GetS0(SV, ploc.X()/(0.5*GetLx())); //$JG2012-02

		d.SetSize(NS()*Dim(),Dim());
		d.FillWithZeros();
		//d = S + ...
		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*Dim()+1,1)=SV(i);
			d((i-1)*Dim()+2,2)=SV(i);
		}
		double y = ploc.Y();
		if (y != 0)
		{
			GetS0x(SV, ploc.X()/(0.5*GetLx()));

			ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
			for (int i = 1; i <= SOS(); i++)
				xg(i) = XG(i);  //e=e1..e8

			Vector2D rx = GetPosx2D(ploc.X(), xg); //rx = S.e
			Vector2D n = Vector2D(-rx.Y(), rx.X()); //normal vector

			double rxnorm = rx.Norm();//sqrt(rxx^2+rxy^2)
			if(rxnorm == 0) rxnorm = 1;

			for (int i = 1; i <= NS(); i++)
			{
				//y/|n|*dn/dq
				d((i-1)*Dim()+1,2) +=  y*SV(i)/rxnorm;
				d((i-1)*Dim()+2,1) += -y*SV(i)/rxnorm;

				//y*n/|n|*(r_x1*S_x1 + r_x2*S_x2)
				d((i-1)*Dim()+1,1) +=  y*n.X()/rxnorm*(rx.X()*SV(i));
				d((i-1)*Dim()+1,2) +=  y*n.Y()/rxnorm*(rx.X()*SV(i));
				d((i-1)*Dim()+2,1) +=  y*n.X()/rxnorm*(rx.Y()*SV(i));
				d((i-1)*Dim()+2,2) +=  y*n.Y()/rxnorm*(rx.Y()*SV(i));
			}
		}
	}


	//gives derivative of d r / dx (r_x) at local position p_loc depending on coordinates xg
	virtual Vector2D GetPosx2D(const double& p_loc, const Vector& xg) const
	{
		double p0=p_loc;
		ConstVector<ANCFCable2D_MaxNShapes> SV(NS());
		GetS0x(SV, p0/(0.5*GetLx()));
		Vector2D p(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*Dim()+i);
			}
		}
		return p;
	}



	//gives second derivative of d2 r / dx2 (r_xx) at local position p_loc depending on coordinates xg
	virtual Vector2D GetPosxx2D(const double& p_loc, const Vector& xg) const
	{
		double p0=p_loc;
		ConstVector<ANCFCable2D_MaxNShapes> SV(NS());
		SV.SetLen(NS());
		GetS0xx(SV, p0/(GetLx()*0.5));
		Vector2D p(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*Dim()+i);
			}
		}
		return p;
	}


	//gives position r at local position p_loc depending on given coordinates xg, or on element local coordinates
	virtual Vector2D GetPos2D(const double& p_loc, const Vector& xg) const
	{
		double p0=p_loc/(0.5*GetLx());// x€[-GetLx()/2,GetLx()/2]=>p0€[-1,1]
		ConstVector<ANCFCable2D_MaxNShapes> SV(NS());
		GetS0(SV, p0);
		Vector2D p(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*Dim()+i);
			}
		}
		return p;
	};

	virtual Vector2D GetPos2D(const double& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		GetCoordinates(xg);
		double p0=p_loc/(0.5*GetLx());// x€[-GetLx()/2,GetLx()/2]=>p0€[-1,1]

		ConstVector<ANCFCable2D_MaxNShapes> SV(NS());
		GetS0(SV, p0);
		Vector2D p(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NS(); j++)
			{
				p(i) += SV(j)*xg((j-1)*Dim()+i);
			}
		}
		return p;
	};

	virtual Vector2D GetPos2DD(const double& p_loc) const;//

	// get velocity vector v at local coordinate ploc
	virtual Vector2D GetVel2D(const double& p_loc) const;

	// get displacement vector at local coordinate ploc
	virtual Vector2D GetDisplacement2D(const double& p_loc, bool flagD=0) const;//
	virtual Vector2D GetDisplacement2DD(const double& p_loc) const {return GetDisplacement2D(p_loc, 1);}

	// get reference position of element
	virtual Vector2D GetRefPos2DD() const {return GetPos2DD(Vector2D(0.,0.));};//

	// get position of idof degree of freedom
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		if (idof <= 4)
		{
			return ToP3D(GetPos2DD(Vector2D(-0.5*GetLx(),0.)));
		}
		else
		{
			return ToP3D(GetPos2DD(Vector2D(0.5*GetLx(),0.)));
		}
	}

	// get direction of of action of i-th DOF
	virtual Vector3D GetDOFDirD(int idof) const
	{
		if (idof > 4) idof -= 4;

		if (idof == 1) return ToP3D(Vector2D(1.,0.));
		else if (idof == 2) return ToP3D(Vector2D(0.,1.));
		else if (idof == 3) return Vector3D(0.,0.,0.);
		else if (idof == 4) return ToP3D(Vector3D(0.,0.,2.));
		return Vector3D(0.,0.,0.);
	}

	//in reference element coordinates (-1..1)
	virtual Vector2D GetPos2D0D(const double& p_loc) const;
	virtual Vector2D GetVel2DD(const double& p_loc) const;


	// get position at local point ploc away from axis
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);

		Vector2D rx=GetPosx2D(p_loc.X(),xg);
		Vector2D n = Vector2D(-rx(2),rx(1)); //Dachvektor
		n.Normalize();
		return GetPos2D(p_loc.X()) + n*p_loc.Y();
	};



	// get x-derivative of position at local point ploc away from axis
	virtual Vector2D GetPosx2D(const Vector2D& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);

		Vector2D rxx=GetPosxx2D(p_loc.X(),xg);
		Vector2D rx=GetPosx2D(p_loc.X(),xg);
		double rxnorm = rx.Norm();

		Vector2D nx = 1./rxnorm*Vector2D(-rxx(2),rxx(1)) + (rxx*rx)/Cub(rxnorm) * Vector2D(-rx(2), rx(1)); //Dachvektor
		return GetPosx2D(p_loc.X(), xg) + nx*p_loc.Y();
	};

	virtual Vector2D GetPosx2DD(const Vector2D& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGD(i);

		Vector2D rxx=GetPosxx2D(p_loc.X(),xg);
		Vector2D rx=GetPosx2D(p_loc.X(),xg);
		double rxnorm = rx.Norm();

		Vector2D nx = 1./rxnorm*Vector2D(-rxx(2),rxx(1)) + (rxx*rx)/Cub(rxnorm) * Vector2D(-rx(2), rx(1)); //Dachvektor
		return GetPosx2D(p_loc.X(), xg) + nx*p_loc.Y();
	};


	// get velocity at local point ploc away from axis
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XG(i);

		Vector2D rx=GetPosx2D(p_loc.X(),xg); //rx = S.e
		Vector2D n = Vector2D(-rx.Y(),rx.X()); //normal vector
		double rxnorm = rx.Norm();//sqrt(rxx^2+rxy^2)
		if (rxnorm == 0) return GetVel2D(p_loc.X());

		n.Normalize();

		ConstVector<ANCFCable2D_MaxDOF> xgp(SOS());
		for (int i = 1; i <= SOS(); i++)
			xgp(i) = XGP(i);

		Vector2D rxp=GetPosx2D(p_loc.X(),xgp); //rxp = S.e_dot
		Vector2D np = Vector2D(-rxp.Y(),rxp.X()); //normal vector
		double rxnormp = (rxp.X()*rx.X()+rxp.Y()*rx.Y())/rxnorm;

		Vector2D n0p = (1./Sqr(rxnorm))*(rxnorm*np-rxnormp*n);

		return GetVel2D(p_loc.X()) + p_loc.Y()*n0p;
	};

	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGD(i);

		Vector2D rx=GetPosx2D(p_loc.X(),xg); //rx = S.e
		Vector2D n = Vector2D(-rx.Y(),rx.X()); //normal vector
		double rxnorm = rx.Norm();//sqrt(rxx^2+rxy^2)
		if (rxnorm == 0) return GetVel2DD(p_loc.X());

		n.Normalize();

		ConstVector<ANCFCable2D_MaxDOF> xgp(SOS());
		for (int i = 1; i <= SOS(); i++)
			xgp(i) = XGPD(i);

		Vector2D rxp=GetPosx2D(p_loc.X(),xgp); //rxp = S.e_dot
		Vector2D np = Vector2D(-rxp.Y(),rxp.X()); //normal vector
		double rxnormp = (rxp.X()*rx.X()+rxp.Y()*rx.Y())/rxnorm;

		Vector2D n0p = (1./Sqr(rxnorm))*(rxnorm*np-rxnormp*n);

		return GetVel2DD(p_loc.X()) + p_loc.Y()*n0p;
	};

	// get position at local point ploc away from axis
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const 
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGD(i);

		Vector2D rx=GetPosx2D(p_loc.X(),xg);
		Vector2D n = Vector2D(-rx(2),rx(1)); //Dachvektor
		n.Normalize();
		return GetPos2DD(p_loc.X()) + n*p_loc.Y(); 
	};

	virtual Vector2D GetPos2D0D(const Vector2D& p_loc) const 
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGD(i);

		Vector2D rx=GetPosx2D(p_loc.X()*0.5*GetLx(),xg);
		Vector2D n = Vector2D(-rx(2),rx(1)); //Dachvektor
		n.Normalize();
		return GetPos2DD(p_loc.X()*0.5*GetLx()) + n*p_loc.Y()*0.5*GetLy(); 
	};

	// get position at local point ploc away from axis, when using deformation scale defscale, routine for drawing only
	virtual Vector2D GetPos2D0D(const Vector2D& p_loc, double defscale) const 
	{
		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		//compute position in reference configuration:
		Vector2D p0;


		Vector2D rx=GetPosx2D(p_loc.X()*0.5*GetLx(),x_init);

		Vector2D n = Vector2D(-rx(2),rx(1)); //Dachvektor
		n.Normalize();

		p0 = GetPos2D(p_loc.X()*0.5*GetLx(),x_init) + n*p_loc.Y()*0.5*GetLy(); 
		//++++++++++++++++++++++++++++++++++++++++++++++++++

		for (int i = 1; i <= SOS(); i++)
			xg(i) = XGD(i);

		rx=GetPosx2D(p_loc.X()*0.5*GetLx(),xg);

		n = Vector2D(-rx(2),rx(1)); //Dachvektor
		n.Normalize();

		return defscale*(GetPos2D(p_loc.X()*0.5*GetLx(),xg) + n*p_loc.Y()*0.5*GetLy() - p0) + p0; 
	};


	
	
	
	
	// Plasticity
	// PlasticStrainMatrix entries are equally distributed values of piecewise bilinear function on the element
	//  € [-GetLx()/2, GetLx()/2] x [-GetLy()/2, GetLy()/2]
	// Matrix entries specify function values, ordered as
	//   [ 1          2         3   ... width ]             [        y=GetLy()/2            ]
	//   [ w+1        w+2       w+3 ... 2w    ]   ----->    [x=-GetLx()/2             x=GetLx()/2]
	//
	//   [ (h-1)w+1  (h-1)w+1  (h-1)w+2  ...  hw]           [        y=-GetLy()/2           ]
	virtual double GetPlasticAxStrain(double x0, int flagD=0) const;    // integral mean value of plastic strains due to plastic cells, x0 € [-0.5GetLx(), 0.5GetLx()]
	virtual double GetPlasticAxStrainIP(int ip, int flagD=0) const;    // integral mean value of plastic strains due to plastic cells, ip is integration point number
	virtual double GetPlasticKappa(double x0, int flagD=0) const;				// first moment of plastic strains due to plastic cells
	virtual double GetPlasticStrain(Vector2D& ploc, int flagD=0) const;  // inter/extrapolated value of plastic strain at local coord. ploc due to plastic cells
	virtual void SetPlasticStrainMatrix( const Matrix& initstrains, int height, int width); // initialize plastic strain matrix
	virtual void SetPlasticStrainMatrix( int height, int width); // initialize plastic strain matrix with zero
	virtual void GetPlasticStrainMatrix(Matrix& mat, int flagD=0);
	virtual void GetPlasticStrainMatrix(const Matrix& mat, int flagD=0) const;
	virtual void GetPlasticStrainMatrixLastStep(const Matrix& mat) const;
#ifdef AND_INCLUDES
	// Transport scheme only available for AND
	virtual void GetPlasticStrainMatrixLastStepTransported(const Matrix& mat) const;
#endif // AND_INCLUDES

	// matrix containing internal variable xi for hardening, xi = int_0^t |plasticstrainP| ds
	virtual void GetInternalVariableMatrix( const Matrix& initvar, int flagD=0) const;
	virtual void GetInternalVariableMatrix( Matrix& initvar, int flagD=0);
	virtual void GetInternalVariableMatrixLastStep(const Matrix& mat) const;
#ifdef AND_INCLUDES
	// Transport scheme only available for AND
	virtual void GetInternalVariableMatrixLastStepTransported(const Matrix& mat) const;
#endif // AND_INCLUDES
	// inter/extrapolated value of plastic strain at local coord. ploc due to plastic cells
	virtual double GetInternalVariable(Vector2D& ploc, int flagD=0) const;



	// Energy
	virtual double GetKineticEnergy();
	virtual double GetPotentialEnergy();

	// useful functions: Compute Kappa, Eps_Axial, delta...
	// with Shape function vectors provided
	virtual double GetKappa(const double& x, const Vector& xg) const;
	virtual double GetKappaD(const double& x, const Vector& xg) const;
	virtual double GetEpsAxial(const double& x, const Vector& xg) const;
	virtual double GetEpsAxialD(const double& x, const Vector& xg) const;
	virtual void GetDeltaKappa(const double& x, const Vector& xg, const Vector& SVx, const Vector& SVxx, 
		Vector& dkappa, double& kappa) const;
	virtual void GetDeltaDeltaKappa(const double& x, const Vector& xg, const Vector& SVx, const Vector& SVxx, 
		double& kappa, Vector& dkappa, Matrix& ddkappa) const;
	virtual void GetDeltaEpsAxial(const double& x, const Vector& xg, const Vector& SVx, Vector& depsaxial) const;
	virtual void GetDeltaDeltaEpsAxial(const double& x, const Vector& xg, const Vector& SVx, 
		double& eps, Vector& deltaeps, Matrix& ddepsaxial) const;
	virtual double GetEps(const Vector2D& ploc) const;
	virtual double GetEpsD(const Vector2D& ploc) const;
	virtual double GetEpsInit(const Vector2D& ploc, int flagD=0) const;
	// shape function vectors not provided
	virtual void GetDeltaEpsAxial(const double& x, const Vector& xg, Vector& depsaxial) const
	{
		ConstVector<ANCFCable2D_MaxNShapes> SVx(NS());
		GetS0x(SVx,x);
		GetDeltaEpsAxial(x, xg, SVx, depsaxial);
	}
	virtual void GetDeltaDeltaEpsAxial(const double& x, const Vector& xg, 
		double& eps, Vector& deltaeps, Matrix& ddepsaxial) const
	{
		ConstVector<ANCFCable2D_MaxNShapes> SVx(NS());
		GetS0x(SVx,x);
		GetDeltaDeltaEpsAxial(x, xg, SVx, eps, deltaeps, ddepsaxial);
	}
	virtual void GetDeltaKappa(const double& x, const Vector& xg, Vector& dkappa, double& kappa) const
	{
		ConstVector<ANCFCable2D_MaxNShapes> SVx(NS()), SVxx(NS());
		GetS0x(SVx,x);
		GetS0xx(SVxx,x);
		GetDeltaKappa(x, xg, SVx, SVxx, dkappa, kappa);
	}
	virtual void GetDeltaDeltaKappa(const double& x, const Vector& xg,
		double& kappa, Vector& dkappa, Matrix& ddkappa) const
	{
		ConstVector<ANCFCable2D_MaxNShapes> SVx(NS()), SVxx(NS());
		GetS0x(SVx,x);
		GetS0xx(SVxx,x);
		GetDeltaDeltaKappa(x, xg, SVx, SVxx,kappa, dkappa, ddkappa);
	}

	// compute length of a curved element
	virtual double ComputeCurvedLength(const Vector& xg) const;

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

	virtual void GetIntDkappaDq2D(Vector& dudq);

	//z.B. für Gleiten entlang der x-Achse (wird derzeit nicht benötigt)
	virtual void GetdPosdx(const Vector2D& ploc, Vector2D& dpdx)
	{
		//optimierungspotenzial 500% !!!!!!!!!!!!!!!!!!!
		double p0 = ploc.X()/(0.5*GetLx());

		ConstVector<ANCFCable2D_MaxDOF> SV(NS());
		GetS0x(SV,p0);

		ConstVector<ANCFCable2D_MaxDOF> xg(SOS());
		GetCoordinates(xg);
		dpdx = 0;
		for (int i=1; i <= Dim(); i++)
		{
			for (int j=1; j <= NS(); j++)
			{
				dpdx(i) += SV(j) * xg((j-1)*Dim()+i);
			}
		}
	};
	virtual void GetdAngle2Ddx(const Vector2D& ploc, double& dphidx);
	virtual int GetKappaMode() const;

	virtual void GetdAngle2DdqT(const Vector2D& ploc, Matrix& d);
	virtual double GetAngle2D(const Vector2D& ploc) const;
	virtual double GetAngle2DP(const Vector2D& ploc) const;

	virtual double GetCurvature2D(const Vector2D& p_loc) const;
	virtual double GetCurvature2DP(const Vector2D& p_loc) const;
	virtual double GetCurvature2DD(const Vector2D& p_loc) const {return 0;}


	virtual double GetBeamEIy() const {return GetMaterial().BeamEIy(); }
	virtual double GetBeamEA() const {return GetMaterial().BeamEA(); }
	virtual double GetRhoA() const {return GetMaterial().BeamRhoA(); }


	virtual double& GetEpsDamping() {return epsdamping;}
	virtual const double& GetEpsDamping() const {return epsdamping;}

	// weight for initial curvature: 0 means element is initially straight and stressless in straight config,
	//   1 means element is initially curved and stressless in curved config
	//   this loadfact is not used everywhere! only tested in ANCFAxMovBeam2D!!
	virtual double GetInitLoadfact(int flagD=0) const {return 1;};

	virtual void SetArcLengthParameter(int elem, double load, int ldof);

	// AP: use tangent stiffness matrix for plastic computations
	// attention: experimental, does not work correctly
	virtual int& UseTangentStiffness() {return use_tangent_stiffness_matrix;}
	virtual const int& UseTangentStiffness() const {return use_tangent_stiffness_matrix;}
protected:
	//element dimensions
	//double lx, ly, lz;

	// node numbers
	int n1; //$EDC$[varaccess,EDCvarname="node_number1",EDCfolder="Geometry",tooltiptext="global number of node 1 (node must exist, add with AddNode(...))"] //DO NOT MODIFY THIS COMMENT!!!!
	int n2; //$EDC$[varaccess,EDCvarname="node_number2",EDCfolder="Geometry",tooltiptext="global number of node 2 (node must exist, add with AddNode(...))"] //DO NOT MODIFY THIS COMMENT!!!!

	// number of degrees of freedom added by element
	int sos2;

	int kappanodal; //add second derivative in node

	double concentratedmass1; //$EDC$[varaccess,EDCvarname="concentrated_mass_node1",EDCfolder="Physics",tooltiptext="a concentrated mass is added at position of node 1"] //DO NOT MODIFY THIS COMMENT!!!!
	double concentratedmass2; //$EDC$[varaccess,EDCvarname="concentrated_mass_node2",EDCfolder="Physics",tooltiptext="a concentrated mass is added at position of node 2"] //DO NOT MODIFY THIS COMMENT!!!!

	Matrix massmatrix, Hmatrix; //M = int(rho*((S)^T).S, dV,V); H = int(S,dV,V)

	Vector e0; //initial vector in e, not in p

	// Plasticity
	int plasticstrains_height;
	int plasticstrains_width;
	int use_tangent_stiffness_matrix;

	// Axial strain damping factor
	double epsdamping;

	// for arc-length method
	int numConstraintElem;
	double gcload;
	int loaddof;
  //EDC Vector3D size; //$EDC$[varaccess,EDCvarname="element_size",EDCfolder="Geometry",tooltiptext="size of element: length (lx), height (ly), out of plane width (lz)"] //DO NOT MODIFY THIS COMMENT!!!!
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node1_initial_position",EDCfolder="Initialization",vecstart=1,vecend=2,tooltiptext="initial values for position of node 1 [x1,y1]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node1_initial_slope",EDCfolder="Initialization",vecstart=3,vecend=4,tooltiptext="initial values for slope vector of node 1 [r1xx,r1xy]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node2_initial_position",EDCfolder="Initialization",vecstart=5,vecend=6,tooltiptext="initial values for position of node 2 [x2,y2]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node2_initial_slope",EDCfolder="Initialization",vecstart=7,vecend=8,tooltiptext="initial values for slope vector of node 2 [r2xx,r2xy]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node1_initial_velocity",EDCfolder="Initialization",vecstart=1+8,vecend=2+8,tooltiptext="initial values for velocity of node 1 [x1,y1]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node1_initial_slope_vel",EDCfolder="Initialization",vecstart=3+8,vecend=4+8,tooltiptext="initial values for slope velocity vector of node 1 [r1xx,r1xy]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node2_initial_velocity",EDCfolder="Initialization",vecstart=5+8,vecend=6+8,tooltiptext="initial values for velocity of node 2 [x2,y2]"]
	//EDC Vector x_init; //$EDC$[varaccess,EDCvarname="node2_initial_slope_vel",EDCfolder="Initialization",vecstart=7+8,vecend=8+8,tooltiptext="initial values for slope velocity vector of node 2 [r2xx,r2xy]"]
}; //$EDC$[endclass,ANCFCable2D]

const int npiezosensors = 1000;

//rigid cube
class ANCFPiezoCable2D: public ANCFCable2D
{
public:
	//Body3D():Element() {mbs = NULL;};
	ANCFPiezoCable2D(MBS* mbsi):ANCFCable2D(mbsi) {};
	ANCFPiezoCable2D(const ANCFPiezoCable2D& e):ANCFCable2D(e.mbs) {CopyFrom(e);};

	//nodal coordinates of first and second point (x1,x2, x1.Length()==12)
	ANCFPiezoCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI=0): ANCFCable2D(mbsi)
	{
		SetANCFCable2D(xc1, xc2, rhoi, Emi, si, coli, kappanodalI);
	};

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFPiezoCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, int n1i, int n2i, double rhoi, double Emi,
		const Vector3D& si, const Vector3D& coli, int kappanodalI=0): ANCFCable2D(mbsi)
	{
		SetANCFCable2D(xc1, xc2, n1i, n2i, rhoi, Emi, si, coli, kappanodalI); 
	}

	//Element shares nodes with other elements, n1 and n2 are nodenumbers; element sets initial conditions for nodes
	ANCFPiezoCable2D(MBS* mbsi, const Vector& xc1, const Vector& xc2, const Vector& vc1, const Vector& vc2,
		int n1i, int n2i, double rhoi, double Emi,const Vector3D& si, const Vector3D& coli, int kappanodalI=0): 
	ANCFCable2D(mbsi)
	{
		SetANCFCable2D(xc1, xc2, vc1, vc2, n1i, n2i, rhoi, Emi, si, coli, kappanodalI); 
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ANCFPiezoCable2D(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		ANCFCable2D::CopyFrom(e);
		const ANCFPiezoCable2D& ce = (const ANCFPiezoCable2D&)e;

		for (int i=0; i < npiezosensors; i++)
		{
			sensors[i] = ce.sensors[i];
		}
		for (int i=0; i < npiezosensors; i++)
		{
			par[i] = ce.par[i];
		}
	}

	virtual void SetSensor(int i, int sensornum) {sensors[i-1] = sensornum;}
	virtual void SetPar(int i, double val) {par[i-1] = val;}
	virtual double GetPar(int i) {return par[i-1];}
	virtual double EvaluatePiezoMoment(double x_loc);

	virtual void EvalF2(Vector& f, double t);

private:
	int sensors[npiezosensors];
	double par[npiezosensors];

};


#endif
