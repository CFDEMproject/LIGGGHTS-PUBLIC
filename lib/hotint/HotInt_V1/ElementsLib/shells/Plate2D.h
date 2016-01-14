//#**************************************************************
//#
//# filename:             Plate2D.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						April 2006
//# description:          2D-plane stress plate element for floating frame of reference or absolute coordinates
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
 
#ifndef Plate2D__H
#define Plate2D__H

#include "IntegrationRule.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  Plate2D  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


const int Plate2DMaxIP = 16; //4*4 for quadratic plate ...
const int Plate2DMaxnn = 9; //maximum number of nodes
//rigid cube
class Plate2D: public Body2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Plate2D(MBS* mbsi):Body2D(mbsi),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), w1(), w2(),
			xg(), xgd(), x_p0(), geom_nonlin(1)
	{
		elementname = GetElementSpec();
	};
	Plate2D(const Plate2D& e):Body2D(e.mbs),massmatrix(), Hmatrix(), SV(), DS(), x1(), x2(), w1(), w2(),
			xg(), xgd(), x_p0(), geom_nonlin(1) {CopyFrom(e);};

	void SetPlate2D(int bodyindi, const TArray<int>& nodelist, int material_num, double thickness, const Vector3D& coli, int CMSelememti=0);

	//To be overwritten in derived class:
	// not needed for an abstract class - YV
	/*
	virtual Element* GetCopy()
	{
		Element* ec = new Plate2D(*this);
		return ec;
	}*/

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Body2D::CopyFrom(e);
		const Plate2D& ce = (const Plate2D&)e;
		lz = ce.lz;
		// Em = ce.Em;
		// nu = ce.nu;
		//size = ce.size;
		xg = ce.xg;
		xgd = ce.xgd;
		massmatrix = ce.massmatrix;
		Hmatrix = ce.Hmatrix;
		temp = ce.temp;
		SV = ce.SV;
		DS = ce.DS;

		//integration points
		x1 = ce.x1; x2 = ce.x2; 
		w1 = ce.w1; w2 = ce.w2; 
		orderxy = ce.orderxy;
		orderxyM = ce.orderxyM;
		orderxyH = ce.orderxyH;

		for (int i=0; i < Plate2DMaxIP; i++)
		{
			grad[i] = ce.grad[i];
			jacdet[i] = ce.jacdet[i];
			//inel_strainIP[i] = ce.inel_strainIP[i];
		}

		sos2 = ce.sos2;
		for (int i = 1; i <= NNodes(); i++)
		{
			NodeNum(i) = ce.NodeNum(i);
			//concentratedmass[i-1] = ce.concentratedmass[i-1];
		}
		bodyind = ce.bodyind;
		x_p0 = ce.x_p0;

		geom_nonlin = ce.geom_nonlin;
	}

	virtual void Initialize() 
	{
		Body2D::Initialize();

		//compute mass of element:
		mass = 0;
		int dim = Dim(); 
		int ns = NS();

		Matrix3D jac;
		jac.SetSize(dim,dim);

		GetIntegrationRule(x1,w1,orderxyM);
		GetIntegrationRule(x2,w2,orderxyM);

		/*
		for (int i=0; i < Plate2DMaxIP; i++)
		{
			inel_strainIP[i] = Vector3D(0.,0.,0.);
		}
		*/

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));

				GetJacobi(jac,p);

				double jacdet = jac.Det();
				double rho = GetMBS()->GetMaterial(GetMaterialNum()).Density(); //$ DR 2013-02-04 deleted rho from class element, do not use it here!
				//mass += fabs (jacdet) * Rho() * w1(i1)*w2(i2) * lz; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
				mass += fabs (jacdet) * rho * w1(i1)*w2(i2) * lz;
			}
		}

	}

	virtual void LinkToElements();
	virtual void BuildDSMatrices(); 

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* GetElementSpec() const { return "Plate2D"; }
	virtual TFiniteElementType GetElementType() const = 0;


	virtual int SOS() const {return Dim()*NS();}; //size of K and M
	virtual int SOSowned() const {return 0;}; //number of internal degrees of freedom
	virtual int ES() const  {return 0;};  //size of first order explicit equations
	virtual int IS() const  {return 0;};  //implicit (algebraic) size
	//xxxxxxxxxxxxxxxxx

	virtual int Nip() const {return 1;} //default

	virtual int FlexDOF() const {return Dim()*NS();}

	virtual int IsRigid() const {return 0;}
	virtual int IsTrig() const {return 0;}
	//virtual void SetConcentratedMass(double cm, int i) {concentratedmass[i-1] = cm;}


// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+SOS()));}
	//virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg.Get(iloc+SOS()));}
	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg.Get(iloc+SOS()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+SOS()));}
//	//virtual const double& XGPD(int iloc) const {return mbs->GetDrawValue(ltg(iloc+SOS()));}
//	virtual const double& XGPD(int iloc) const {return GetDrawValue(ltg(iloc+SOS()));}

	virtual const bool& GeometricallyNonlin() const { return geom_nonlin; }
	virtual bool& GeometricallyNonlin() { return geom_nonlin; }


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//change for other shape functions and dimension:
	virtual int Dim() const {return 2;}
	virtual int NS() const {return NNodes();}

	virtual int NNodes() const {GetMBS()->UO() << "ERROR: Call to Plate2D::NNodes() illegal\n"; return 0;};
	virtual const int& NodeNum(int i) const {return node[i-1];}
	virtual int& NodeNum(int i) {return node[i-1];}

	virtual const Node& GetNode(int i) const {return GetMBS()->GetNode(NodeNum(i));}
	virtual Node& GetNode(int i) {return GetMBS()->GetNode(NodeNum(i));}

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const {GetMBS()->UO() << "ERROR: Call to Plate2D::GetS0 illegal\n";};
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const {GetMBS()->UO() << "ERROR: Call to Plate2D::GetDSMatrix0 illegal\n";};
	virtual double GetS0(const Vector2D& ploc, int i) const {return 0;};
	virtual double GetDS0(const Vector2D& ploc, int shape, int dxj) const {return 0;};
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//floating frame of reference formulation: for FFRF elements

	virtual void StiffnessMatrix(Matrix& m); //fill in sos x sos components, m might be larger
	virtual void GetI1(Vector& I1); //Shabana p. 209-211
	virtual double GetIkl(int k, int l);
	virtual void GetIbarkl(int k, int l, Vector& I1);
	virtual void GetSbar(Matrix& Sbar);
	virtual void GetSbarkl(int k, int l, Matrix& Sbar);

	//Position relative to frame, undeformed
	virtual Vector2D GetPos2Drel0(const Vector2D& p_loc) const;
	//relative positions (deformed):
	virtual Vector2D GetPos2Drel(const Vector2D& p_loc) const;
	virtual Vector2D GetVel2Drel(const Vector2D& p_loc) const;
	virtual Vector2D GetPos2DrelD(const Vector2D& p_loc, int use_magnification) const;
	virtual Vector2D GetVel2DrelD(const Vector2D& p_loc) const;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//-1..+1 based!!!
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetDisplacement2D(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetPos2DD(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetPos2DD(const Vector2D& p_loc, int use_magnification) const; //deformation scale factor
	//-1..+1 based!!!
	virtual Vector2D GetDisplacement2DD(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	virtual Vector2D GetRefPos2D() const;
	virtual Vector2D GetRefPos2DD() const; 


	virtual Box3D GetElementBox() const
	{
		Box3D b;
		for (int i = 1; i <= NS(); i++)
		{
			b.Add(ToP3D(GetPos2D(GetNodeLocPos2D(i))));
		}
		return b;
	}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b;
		for (int i = 1; i <= NS(); i++)
		{
			b.Add(ToP3D(GetPos2DD(GetNodeLocPos2D(i))));
		}
		return b;
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

	//for compute+draw:
	virtual void GetJacobi(Matrix3D& jac, const Vector2D& ploc) const
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
					jac(i,j) += GetDS0(ploc,k,j)*x_p0((k-1)*Dim()+i);
				}
			}
		}
	}

	virtual void JacobianF2(double t, Matrix& m, IVector& colref);

	virtual void EvalF2(Vector& f, double t); 

  virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	//ploc -1 ... +1
	virtual void Gradu(const Vector2D& ploc, const Vector& u, Matrix3D& gradu) const
	{
		static Matrix DS;
		GetDSMatrix0(ploc,DS);

		Matrix3D jac, jacinv;
		jac.SetSize(2,2); jacinv.SetSize(2,2);
		
		GetJacobi(jac,ploc);
		//GetJacobi(jac,ploc,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		static Matrix grad;
		grad.SetSize(Dim(),NS());
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

	virtual void StrainVector(int ind, const Vector& u, Vector3D& strain) const
	{
		int ns = NS();
		int dim = Dim();
		int ii = 0;

		if (GeometricallyNonlin() )
		{
			// compute F ... Deformation gradient
			Matrix3D F, FTF;
			F.SetSize(2,2); FTF.SetSize(2,2);
			F.SetAll(0);

			for (int i = 1; i <= ns; i++)
			{
				for (int j = 1; j <= dim; j++) 
				{
					// ii = (i-1)*dim+j;
					ii++;
					for (int k = 1; k <= dim; k++)
					{
						F(j,k) += grad[ind](k,i)*u(ii); //Verschiebungsgradient
					}
				}
			}
			
			F(1,1) += 1.;
			F(2,2) += 1.;

			// Green-Lagrange strain tensor
			FTF = 0.5* (F.GetTp() * F);
			strain(1) = FTF(1,1) - 0.5;
			strain(2) = FTF(2,2) - 0.5;
			strain(3) = 2.*FTF(1,2);
		}
		else // geometrically linear
		{
			Matrix3D G, GTG;  // displacement gradient
			G.SetSize(2,2); GTG.SetSize(2,2);
			G.SetAll(0);

			for (int i = 1; i <= ns; i++)
			{
				for (int j = 1; j <= dim; j++) 
				{
					// l = (i-1)*dim+j;
					ii++;
					for (int k = 1; k <= dim; k++)
					{
						G(j,k) += grad[ind](k,i)*u(ii); //Verschiebungsgradient
					}
				}
			}

			// linearized strain tensor
			strain(1) = G(1,1);
			strain(2) = G(2,2);
			strain(3) = ( G(1,2)+G(2,1) );
		}
	}

	virtual void FirstPiolaTensor(int ind, const Vector& u, Vector3D& sigma, Matrix3D& piola1) const
	{
		piola1(1,1) = sigma(1);
		piola1(2,2) = sigma(2);
		piola1(1,2) = sigma(3);
		piola1(2,1) = sigma(3);

		// in case of geometrically linear material: piola1 = stress tensor
		if (!GeometricallyNonlin() ) return;

		int ns = NS();
		int dim = Dim();
		int l=0;

		// compute F ... Deformation gradient
		Matrix3D G;
		G.SetSize(2,2);
		G.SetAll(0);

		for (int i = 1; i <= ns; i++)
		{
			for (int j = 1; j <= dim; j++) 
			{
				//l = (i-1)*dim+j;
				l++;
				for (int k = 1; k <= dim; k++)
				{
					G(j,k) += grad[ind](k,i)*u(l); //Verschiebungsgradient
				}
			}
		}

		piola1 += G * (piola1);

	}


	virtual void FirstPiolaTensor(int ind, const Vector& u, Matrix3D& piola1) const
	{
		int ns = NS();
		int dim = Dim();
		int ii = 0;

		Vector3D strain, sigma;
		Matrix3D F;

		//$!PG 2011-3-15:[
		//$!PG 2011-3-15: replace GetMaterialTensor2D(Matrix3D&) by ComputeElasticityMatrix(Matrix&)
		//Matrix3D C;
		//GetMaterial().GetMaterialTensor2D(C);
		Matrix C;
		GetMaterial().ComputeElasticityMatrix(C);
		//$!PG 2011-3-15:]

		if (GeometricallyNonlin() )
		{
			// compute F ... Deformation gradient
			Matrix3D FTF;
			F.SetSize(2,2); FTF.SetSize(2,2);
			F.SetAll(0);

			for (int i = 1; i <= ns; i++)
			{
				for (int j = 1; j <= dim; j++) 
				{
					// l = (i-1)*dim+j;
					ii++;
					for (int k = 1; k <= dim; k++)
					{
						F(j,k) += grad[ind](k,i)*u(ii); //Verschiebungsgradient
					}
				}
			}
			F(1,1) += 1.;
			F(2,2) += 1.;

			// Green-Lagrange strain tensor
			FTF = 0.5* (F.GetTp() * F);
			strain(1) = FTF(1,1) - 0.5;
			strain(2) = FTF(2,2) - 0.5;
			strain(3) = 2*FTF(1,2);

		}
		else // geometrically linear
		{
			Matrix3D G;  // displacement gradient
			G.SetSize(2,2);
			G.SetAll(0);

			for (int i = 1; i <= ns; i++)
			{
				for (int j = 1; j <= dim; j++) 
				{
					// l = (i-1)*dim+j;
					ii++;
					for (int k = 1; k <= dim; k++)
					{
						G(j,k) += grad[ind](k,i)*u(ii); //Verschiebungsgradient
					}
				}
			}

			// linearized strain tensor
			strain(1) = G(1,1);
			strain(2) = G(2,2);
			strain(3) = ( G(1,2)+G(2,1) );
		}

		sigma = C*strain;

		piola1(1,1) = sigma(1);
		piola1(2,2) = sigma(2);
		piola1(1,2) = sigma(3);
		piola1(2,1) = sigma(3);

		if(GeometricallyNonlin())
			piola1 = F*piola1;
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

	//only displacements, rotations makes no sense, even in rigid body
	//->only for volumeloads (gravity ...)
	virtual void GetIntDuDq(Matrix& dudq) //in fact it is DuDq Transposed
	{
		GetH(dudq);
	}
/*
	virtual void GetIntDuDqFCentrifugal(Matrix& dudq, const Vector3D& omega, const Vector3D& r0)
	{
		int dim = Dim();
		int ns = NS();

		Matrix H;	
		H.SetSize(ns*dim,dim);
		H.SetAll(0);
		Matrix3D jac;
		jac.SetSize(dim,dim);
	
		DS.SetSize(dim,ns);
		SV.SetLen(ns);

		GetIntegrationRule(x1,w1,orderxyH);
		GetIntegrationRule(x2,w2,orderxyH);

		for (int i1=1; i1<=x1.GetLen(); i1++)
		{
			for (int i2=1; i2<=x2.GetLen(); i2++)
			{
				Vector2D p(x1(i1),x2(i2));
				GetS0(SV,p);

				GetDSMatrix0(p,DS);
				GetJacobi(jac,p);
				//GetJacobi(jac,p,DS);
				double jacdet = jac.Det();
				double fact = fabs (jacdet) * w1(i1)*w2(i2)*lz;

				Vector3D r_rel = r0 - GetPos(p);
				Vector3D force = omega.Cross(omega.Cross(r_rel));

				for (int i=0; i<ns; i++)
				{
					double combinedfactor = fact*Rho();

					
					for (int j=1; j<=dim; j++)
					{
						H(i*dim+j,j)+=fact*SV(i+1);
					}
				}
			}
		}
		Hmatrix = H;
	}
	for (IntegrationPointsIterator ip(IntegrationRuleLoad); !ip.IsEnd(); ++ip)
	{
		GetJacobi(jac, ip.Point().MakeV2D());
		double jacdet = jac.Det();
		double fact = fabs (jacdet) * ip.Weight();
		Vector3D r_rel = r0 - GetPos(ip.Point());
		Vector3D force = omega.Cross(omega.Cross(r_rel));
		for (int i = 0; i < ns; i++)
		{
			double combinedfactor = fact*GetS0(ip.Point().MakeV2D(), i+1)*Rho();
			H(i*dim+1,1) += combinedfactor * force.X();
			H(i*dim+2,1) += combinedfactor * force.Y();
	//		H(i*dim+3,1) += combinedfactor * force.Z();
		}
	}
}
*/
	virtual void GetdPosdqT(const Vector2D& ploc, Matrix& d)
	{
		//p = S(p.x,p.y,p.z)*q; d(p)/dq 
		Vector2D p0=ploc;
		static Vector SV;
		GetS0(SV, p0);
		d.SetSize(NS()*Dim(),Dim());
		d.FillWithZeros();
		for (int i = 1; i <= NS(); i++)
		{
			d((i-1)*Dim()+1,1)=SV(i);
			d((i-1)*Dim()+2,2)=SV(i);
		}
	}

	//ploc -1 ... +1
	virtual void GraduD(const Vector2D& ploc, Matrix3D& gradu) const
	{
		//xgd is already set at beginning of DrawElement!!!!

		static Matrix DS;
		GetDSMatrix0(ploc,DS);

		Matrix3D jac, jacinv;
		GetJacobi(jac,ploc);
		//GetJacobi(jac,ploc,DS,e0);

		jac.GetInverse(jacinv);
		jacinv = jacinv.GetTp();

		static Matrix grad;
		grad.SetSize(Dim(),NS());
		Mult(jacinv, DS, grad);

		gradu.SetSize(2,2);
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
					gradu(j,k) += grad(k,i)*xgd(l);
				}
			}
		}
	}


	virtual void GetNodedPosdqT(int node, Matrix& dpdqi) 
	{
		dpdqi((node-1)*2+1,1) = 1; //dnx/dnx
		dpdqi((node-1)*2+2,2) = 1; //dny/dny
	};

	virtual void AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f) 
	{
		f((node-1)*2+1) += lambda(1);
		f((node-1)*2+2) += lambda(2);
	};   // f += dpdq*lambda

	virtual Vector2D GetNodeLocPos2D(int i) const
	{
		switch(i)
		{ //1..4 for linear element, 1..9 for quadratic element
		case 1: return Vector2D( 1., 1.); break;
		case 2: return Vector2D(-1., 1.); break;
		case 3: return Vector2D(-1.,-1.); break;
		case 4: return Vector2D( 1.,-1.); break;
		case 5: return Vector2D( 0., 1.); break;
		case 6: return Vector2D(-1., 0.); break;
		case 7: return Vector2D( 0.,-1.); break;
		case 8: return Vector2D( 1., 0.); break;
		case 9: return Vector2D( 0., 0.); break;
		default: return Vector2D(0,0);
		}
	}
	virtual Vector3D GetDOFPosD(int idof) const //returns postion of i-th DOF
	{
		int node = (idof-1)/2+1;
		return ToP3D(GetPos2DD(GetNodeLocPos2D(node)));
	}
	virtual Vector3D GetDOFDirD(int idof) const //returns direction of action of i-th DOF
	{
		int dir = (idof-1)%2;
		if (dir == 0) return Vector3D(1.,0.,0.);
		else if (dir == 1) return Vector3D(0.,1.,0.);
		else return Vector3D(0.,0.,1.);
	}

	virtual Vector2D GetNodePos2D(int i) const
	{
		return Vector2D(XG((i-1)*2+1)+x_p0((i-1)*2+1), XG((i-1)*2+2)+x_p0((i-1)*2+2));
	}

	virtual Vector2D GetNodePos2DD(int i) const 
	{
		return Vector2D(XGD((i-1)*2+1)+x_p0((i-1)*2+1), XGD((i-1)*2+2)+x_p0((i-1)*2+2));
	}

	virtual Vector2D GetNodeVel2D(int i) const
	{
		return Vector2D(XGP((i-1)*2+1), XGP((i-1)*2+2));
	}

	virtual int Is_Planestress() const {return !GetMaterial().IsPlaneStrain();}


	// variables, available for post-processing and sensing
	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables);
	// computation of the variables
	virtual double GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD);

	virtual void DrawElementPreProc();
	virtual void DrawElement();

protected:
	//mechanical: 
	// thickness
	double lz; 
	// Material properties --> into material
	//double Em, nu;

	bool geom_nonlin;

	int node[Plate2DMaxnn];
	int bodyind;
	int sos2;
	//double concentratedmass[Plate2DMaxnn];
	Vector x_p0; //nodal positions in undeformed configuration

	//temporary storage for acceleration:
	Vector xg, xgd;
	Matrix massmatrix;
	Matrix Hmatrix;
	Matrix DS; //
	Vector SV; //
	Vector temp;

	//integration points
	Vector x1,x2,w1,w2;
	int orderxy, orderxyM, orderxyH;
	Matrix grad[Plate2DMaxIP];
	double jacdet[Plate2DMaxIP]; 
};


class Plate2Dlin: public Plate2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Plate2Dlin(MBS* mbsi):Plate2D(mbsi) {};
	Plate2Dlin(const Plate2Dlin& e):Plate2D(e.mbs) {CopyFrom(e);};

	//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
	//initial velocities!!!!!
	Plate2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli, int* nodenums = 0);

	Plate2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		int matnr, double thickness, const Vector3D& coli, int* nodenums = 0);

	void SetPlate2Dlin(int bodyindi, const TArray<int>& nodelist, int material_num, double thickness, const Vector3D& coli, int CMSelememti=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Plate2Dlin(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Plate2D::CopyFrom(e);
		//const Plate2Dlin& ce = (const Plate2Dlin&)e;

		//include all data elements here!
	}

	virtual const char* GetElementSpec() const { return "Quad2Dlin"; }
	virtual TFiniteElementType GetElementType() const { return TFE_Quadrilateral; }

	virtual int NNodes() const {return 4;}
	virtual int Nip() const {return 2;} //integration points for one dimension

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;
	virtual double GetS0(const Vector2D& ploc, int i) const;
	virtual double GetDS0(const Vector2D& ploc, int shape, int dxj) const;
};


//nine-node plate element:
class Plate2Dquad: public Plate2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Plate2Dquad(MBS* mbsi):Plate2D(mbsi) {};
	Plate2Dquad(const Plate2Dquad& e):Plate2D(e.mbs) {CopyFrom(e);};

	//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
	//initial velocities!!!!!
	Plate2Dquad(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli);

	Plate2Dquad(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& xc5, const Vector2D& xc6, const Vector2D& xc7, const Vector2D& xc8, const Vector2D& xc9, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		const Vector2D& vc5, const Vector2D& vc6, const Vector2D& vc7, const Vector2D& vc8, const Vector2D& vc9, 
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli);

	void SetPlate2Dquad(int bodyindi, const TArray<int>& nodelist, int material_num, double thickness, const Vector3D& coli, int CMSelememti=0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Plate2Dquad(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Plate2D::CopyFrom(e);
		//const Plate2Dquad& ce = (const Plate2Dquad&)e;

		//include all data elements here!
	}

	virtual const char* GetElementSpec() const {return "Quad2Dquad";}
	virtual TFiniteElementType GetElementType() const { return TFE_Quadrilateral; }

	virtual int NNodes() const {return 9;}
	virtual int Nip() const {return 4;} //integration points for one dimension

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;
	virtual double GetS0(const Vector2D& ploc, int i) const;
	virtual double GetDS0(const Vector2D& ploc, int shape, int dxj) const;
};



class Trig2Dlin: public Plate2D
{
public:
	//Body2D():Element() {mbs = NULL;};
	Trig2Dlin(MBS* mbsi):Plate2D(mbsi) {};
	Trig2Dlin(const Trig2Dlin& e):Plate2D(e.mbs) {CopyFrom(e);};

	//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
	//initial velocities!!!!!
	Trig2Dlin(MBS* mbsi, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3,
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli);

  void Trig2Dlin::SetTrig2Dlin(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti = 0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Trig2Dlin(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Plate2D::CopyFrom(e);
		//const Trig2Dlin& ce = (const Trig2Dlin&)e;

		//include all data elements here!
	}
	virtual const char* GetElementSpec() const {return "Trig2Dlin";}
	virtual TFiniteElementType GetElementType() const { return TFE_Triangle; }

	virtual Vector2D GetRefPos2D() const 
	{
		return Vector2D(1./3.*(XG(1)+XG(3)+XG(5)+x_p0(1)+x_p0(3)+x_p0(5)), 1./3.*(XG(2)+XG(4)+XG(6)+x_p0(2)+x_p0(4)+x_p0(6)));
	};
	virtual Vector2D GetRefPos2DD() const
	{
		return Vector2D(1./3.*(XGD(1)+XGD(3)+XGD(5)+x_p0(1)+x_p0(3)+x_p0(5)), 1./3.*(XGD(2)+XGD(4)+XGD(6)+x_p0(2)+x_p0(4)+x_p0(6)));
	};

	virtual int IsTrig() const {return 1;}

	virtual int NNodes() const {return 3;}

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;
	virtual double GetS0(const Vector2D& ploc, int i) const;
	virtual double GetDS0(const Vector2D& ploc, int shape, int dxj) const;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//overwritten functions only for triangular elements:
	virtual void GetH(Matrix& H);

	virtual void EvalM(Matrix& m, double t);

	//for compute+draw:
	virtual void GetJacobi(Matrix3D& jac, const Vector2D& ploc) const
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
					jac(i,j) += GetDS0(ploc,k,j)*x_p0((k-1)*Dim()+i);
				}
			}
		}
	}

	virtual Vector2D GetNodeLocPos2D(int i) const
	{
		switch(i)
		{ //1..3 for linear element, 1..6 for quadratic element
		case 1: return Vector2D(1., 0.); break;
		case 2: return Vector2D(0., 1.); break;
		case 3: return Vector2D(0., 0.); break;
		case 4: return Vector2D(0.5, 0.5); break;
		case 5: return Vector2D(0., 0.5); break;
		case 6: return Vector2D(0.5, 0.); break;
		default: return Vector2D(0,0);
		}
	}

	virtual void EvalF2(Vector& f, double t); 

	virtual void DrawElement();


};

class Trig2Dquad: public Trig2Dlin
{
public:
	//Body2D():Element() {mbs = NULL;};
	Trig2Dquad(MBS* mbsi):Trig2Dlin(mbsi) {};
	Trig2Dquad(const Trig2Dquad& e):Trig2Dlin(e.mbs) {CopyFrom(e);};

	//nodal coordinates of points in 2D
	//if interpolate_points == 1, only 3 points are provided and the other are interpolated!
	Trig2Dquad(MBS* mbsi, int bodyindi, Vector2D* xc, Vector2D* vc, int interpolate_points,
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli);

	virtual void InitializeElement(int bodyindi, Vector2D* p, Vector2D* v,
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli);

	void Trig2Dquad::SetTrig2Dquad(int bodyindi, const TArray<int> &nodelist, int matnr, double thickness, const Vector3D &coli, int CMSelememti = 0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Trig2Dquad(*this);
		return ec;
	}

	virtual const char* GetElementSpec() const {return "Trig2Dquad";}
	virtual TFiniteElementType GetElementType() const { return TFE_Triangle; }

	//from triglin:
	//virtual void CopyFrom(const Element& e)
	//virtual int IsTrig() const {return 1;}

	virtual int NNodes() const {return 6;}

	virtual void GetS0(Vector& sf, const Vector2D& ploc) const;
	virtual void GetDSMatrix0(const Vector2D& ploc, Matrix& sf) const;
	virtual double GetS0(const Vector2D& ploc, int i) const;
	virtual double GetDS0(const Vector2D& ploc, int shape, int dxj) const;
};




//nine-node plate element:
class Plate2DquadFFRF: public Plate2Dquad
{
public:
	//Body2D():Element() {mbs = NULL;};
	Plate2DquadFFRF(MBS* mbsi):Plate2Dquad(mbsi) {};
	Plate2DquadFFRF(const Plate2DquadFFRF& e):Plate2Dquad(e.mbs), 
		Sbar_tilde(), Sbar_tildeSM(), K(),	SbarS(), SbarSM(), Ibar11S(), Ibar12S(), Ibar21S(), Ibar22S(), I1S() {CopyFrom(e);};


	//nodal coordinates of points in 2D; optional rotation matrix should be included to couple with 3D motion!
	//initial velocities!!!!!
	Plate2DquadFFRF(MBS* mbsi, int FFRFindex, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli, int isCMSi = 0);

	Plate2DquadFFRF(MBS* mbsi, int FFRFindex, int bodyindi, const Vector2D& xc1, const Vector2D& xc2, const Vector2D& xc3, const Vector2D& xc4, 
		const Vector2D& xc5, const Vector2D& xc6, const Vector2D& xc7, const Vector2D& xc8, const Vector2D& xc9, 
		const Vector2D& vc1, const Vector2D& vc2, const Vector2D& vc3, const Vector2D& vc4, 
		const Vector2D& vc5, const Vector2D& vc6, const Vector2D& vc7, const Vector2D& vc8, const Vector2D& vc9, 
		double rhoi, double Emi, double nui, double thickness, const Vector3D& coli, int isCMSi = 0);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Plate2DquadFFRF(*this);
		return ec;
	}

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Plate2Dquad::CopyFrom(e);
		const Plate2DquadFFRF& ce = (const Plate2DquadFFRF&)e;

		FFRFind = ce.FFRFind;
		K = ce.K;
		Sbar_tilde = ce.Sbar_tilde;
		SbarS = ce.SbarS;
		Sbar_tildeSM = ce.Sbar_tildeSM;
		SbarSM = ce.SbarSM;
		Ibar11S = ce.Ibar11S;
		Ibar12S = ce.Ibar12S;
		Ibar21S = ce.Ibar21S;
		Ibar22S = ce.Ibar22S;
		I1122S = ce.I1122S;

		//include all data elements here!
	}

	virtual const char* GetElementSpec() const {return "Quad2DquadFFRF";}
	virtual TFiniteElementType GetElementType() const { return TFE_Quadrilateral; }

	virtual void Initialize();

	virtual void LinkToElements();

	virtual int FFRFDim() const {return 3;} //position X, Y, and angle
	virtual int SOS() const {return (1-IsCMS())*(Dim()*NS()+FFRFDim());}; //size of K and M
	virtual int SOSowned() const {return 0;}; //size of unknowns added


	virtual void SetComputeCoordinates()
	{
		xg.SetLen(FlexDOF()+FFRFDim());
		for (int i = 1; i <= xg.Length(); i++)
			xg(i) = XG(i);
	}
	virtual void SetDrawCoordinates()
	{
		xgd.SetLen(FlexDOF()+FFRFDim());
		for (int i = 1; i <= xgd.Length(); i++)
			xgd(i) = XGD(i);
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//for CMS Element
	virtual int IsCMS() const {return IsType(TCMS);}; //is part of component mode synthesis?

	/*
	virtual const double& XG(int iloc) const {return GetXact(ltg(iloc));}
	virtual double& XG(int iloc) {return GetXact(ltg(iloc));}
	*/

// (AD) changed () to .Get()
	virtual const double& XGP(int iloc) const {return GetXact(ltg.Get(iloc+FlexDOF()+(1-IsCMS())*FFRFDim()));}
//	virtual const double& XGP(int iloc) const {return GetXact(ltg(iloc+FlexDOF()+(1-IsCMS())*FFRFDim()));}

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

	//virtual const CMSElement2D& CMSElement() const {return (CMSElement2D&)(GetMBS()->GetElement(FFRFind));}
	//virtual CMSElement2D& CMSElement() {return (CMSElement2D&)(GetMBS()->GetElement(FFRFind));}

	virtual double GetAngle2D() const {return XG(NS()*Dim()+3);};
	virtual Vector2D GetRefPos2D() const {return Vector2D(XG(NS()*Dim()+1),XG(NS()*Dim()+2));};
	virtual Vector2D GetRefPos2DD() const {return Vector2D(XGD(NS()*Dim()+1), XGD(NS()*Dim()+2));};
	virtual Vector2D GetRefVel2D() const {return Vector2D(XGP(NS()*Dim()+1),XGP(NS()*Dim()+2));};

	//-1..+1 based!!!
	virtual Vector2D GetPos2D(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetVel2D(const Vector2D& p_loc) const;
	//-1..+1 based!!!
	virtual Vector2D GetPos2DD(const Vector2D& p_loc, int use_magnification) const;
	//-1..+1 based!!!
	virtual Vector2D GetVel2DD(const Vector2D& p_loc) const;

	virtual Vector2D GetNodePos2D(int i) const;

	virtual Vector2D GetNodePos2DD(int i) const; 

	virtual Vector2D GetNodeVel2D(int i) const;

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

private:
	int FFRFind; //index to Element of reference frame
	Matrix K;

	///+++++++++++++++++++++
	//store matrices FFRF:
	Matrix Sbar_tilde;
	Matrix SbarS;
	SparseMatrix Sbar_tildeSM;
	SparseMatrix SbarSM;
	Vector Ibar11S,Ibar21S,Ibar12S,Ibar22S;
	Vector I1S;
	double I1122S;
	//++++++++++++++++++++++

};





#endif

