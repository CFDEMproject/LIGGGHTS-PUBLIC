//#**************************************************************
//#
//# filename:             GCMSRotorElement.cpp      
//#
//# author:               Pascal Ziegler, Alexander Humer
//#
//# generated:						2012
//# description:          Reference frame element for FFRF formulation
//#                       
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
 
//#include "MBS_includes_element.h"
#include "element.h"
#include "material.h"
#include "Node.h"
#include "rigid3D.h"
#include "Rigid3DKardan.h"
#include "femeshinterface.h"
#include "referenceframe3D.h"
#include "BaseCMSElement.h"
#include "GCMSRotorElement.h"
#include "lapack_routines.h"
#include "FiniteElementGeneric.h"

#include "femathhelperfunctions.h"


// set number of internal modes
// compute initial rigid body degrees of freedom,
// set size and color vectors
template <class RIGID>
void GCMSRotorElement<RIGID>::SetGCMSRotorElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D rotation_axis, double omega, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli)
{
	n_zeromodes = -1;
	type = (TMBSElement)(type|TCMSflag);

	nbmodes = 0;
	nimodes = nimodesi;
	col = coli;
	size = sizei;

	InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
	boundarydofelem.SetLen(0);
	boundarynode.SetLen(0);

	solverparameters.DefaultInitialize();

	EVfile = mystr("");

	// dofs_per_local_mode = 9;

	// set rotation axis
	if (rotation_axis.Norm() != 0.)
		rotation_axis.Normalize();
	else
		assert(0);
	this->rotation_axis = rotation_axis;

	// compute rigid dofs for p, v, phi, phip
	// other dofs will be appended to x_init later on
	x_init.SetLen(2*SOSRigid()+IS());
	x_init.SetAll(0.);
	Matrix3D A = ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z());


	//Matrix3D Ap = ComputeRotMatrixEulerP(phi.X(),phi.Y(), phi.Z(), phip.X(), phip.Y(), phip.Z()); //$ DR 201-03: old code

	// ====================================================
	//$ DR 2011-04-01:[ ComputeRotMatrixPfromOmegaVector (A,omega=phip) --> gehört in linalg3d
	//double b0, b1, b2, b3;
	////RotMatToQuaternions(ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z()),b0, b1, b2, b3);
	//RotMatToQuaternions(A,b0, b1, b2, b3);

	////initial values for time-wise derivative quaternions beta_t:
	//double beta0 = 2*b0; double beta1 = 2*b1; double beta2 = 2*b2; double beta3 = 2*b3; 
	////the betas are already multiplied with 2, compared to Shabana
	//Matrix G(4,4);
	//G(1,1) = -beta1; G(1,2) =  beta0; G(1,3) = -beta3; G(1,4) =  beta2;
	//G(2,1) = -beta2; G(2,2) =  beta3; G(2,3) =  beta0; G(2,4) = -beta1;
	//G(3,1) = -beta3; G(3,2) = -beta2; G(3,3) =  beta1; G(3,4) =  beta0;
	//G(4,1) =  beta0; G(4,2) =  beta1; G(4,3) =  beta2; G(4,4) =  beta3; //this is the time-derivative of the quaternion condition | |^2 = 1

	//Vector f(phip.X(),phip.Y(),phip.Z(),0.); //omega and zero
	//Vector betap(4); 
	//int rv = G.Solve(f, betap);
	//if (!rv) {GetMBS()->UO() << "GCMS:Initialization: could not determine initial Euler parameter velocities due to singularity!!!\n";}

	//Matrix3D Ap=ComputeRotMatrixPEulerParam(b0,b1,b2,b3,betap(1),betap(2),betap(3),betap(4));

	//$ DR 2011-04-01:]
	// ====================================================

	for (int i=1; i<=3; i++)
	{
		x_init(i) = p(i);
	}
	// rotation dofs
	// rotation starts (rotational angle is zero) from initial rotation given by A
	x_init(4) = 0.;	// cos phi - 1	(not to be confused with Vector phi!!!)
	x_init(5) = 0.;	// sin phi

	/*x_init(6) = A(3,1);
	x_init(7) = A(1,2);
	x_init(8) = A(2,2)-1.;
	x_init(9) = A(3,2);
	x_init(10) = A(1,3);
	x_init(11) = A(2,3);
	x_init(12) = A(3,3)-1.;*/
	// translation velocity dofs
	for (int i=1; i<=3; i++)
	{
		x_init(SOSRigid()+i) = v(i);
	}
	// compute rotational velocity dofs	(see paper)
	x_init(SOSRigid()+4) = 0.;
	x_init(SOSRigid()+5) = omega;
	
	/*x_init(SOSRigid()+6) = Ap(3,1);
	x_init(SOSRigid()+7) = Ap(1,2);
	x_init(SOSRigid()+8) = Ap(2,2);
	x_init(SOSRigid()+9) = Ap(3,2);
	x_init(SOSRigid()+10) = Ap(1,3);
	x_init(SOSRigid()+11) = Ap(2,3);
	x_init(SOSRigid()+12) = Ap(3,3);*/

	mbs->UO() << "x_init = " << x_init << "\n";
}

//must be called before Assemble!!!!!
//link elements, 
//ltg-elements, 
//compute nbmodes, 
//compute full M and K, 
//modal analysis, 
//transformation matrix and modal degrees of freedom
//modal matrices
template <class RIGID>
void GCMSRotorElement<RIGID>::DoModalAnalysis(const TArray<int2>& fixednodes)
{
	UO() << "GCMSRotorElement::DoModalAnalysis()\n";
	UO() << "n-FFRF-elements=" << NFFRFElements() << "\n";
	// -------------
	// link to elements

	UO() << "compute element LTG map \n";
	int sosfull = SOSFull();
	ComputeElementLTG();

	// -------------
	//compute necessary boundary nodes (DOF):
	UO() << "compute boundary modes\n";
	ComputeBoundaryDofIndexList(fixednodes, sosfull);

	// -------------
	//resort ltg according to (B) and (I) nodes:
	ResortLTG(sosfull);

	// -------------
	// compute number of zeromodes
	ComputeNZeromodes();

	// ---------------------
	// count total number of boundary modes
	nbmodes = 0;
	// number of classical boundary modes, which correspond to one dof in the full system
	int nb_modesclassic = 0;
	// number of additional static boundary modes, which correspond to a set of dofs
	int nb_modesstatic = staticmode_nodes.Length();
	// number of dofs in the full system corresponding to additional static modes
	int nb_dofsstatic = 0;
	// number of dofs in the full system which are fixed;
	int nb_dofsfixed = 0;
	// number of boundary nodes total, consisting of classic boundary nodes, additional static nodes, and fixed nodes
	int nb_dofstotal;

	for (int i=1; i <= boundarydof.Length(); i++)
	{
		if (boundarydof(i) == 1) nb_modesclassic++;
		if (boundarydof(i) == 2) nb_dofsstatic++;
		if (boundarydof(i) == 3) nb_dofsfixed++;
	}
	nbmodes = nb_modesclassic + nb_modesstatic;
	nb_dofstotal = nb_modesclassic+nb_dofsfixed+nb_dofsstatic;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//size of stiffness and mass-matrix: sosfull x sosfull
	int nb = NBModes();
	int nf = sosfull-nb_dofstotal; //number of internal (flexible) dof - unreduced
	int nr = NIModes();

	// preliminary number of flexible dofs, 
	// may be reduced due to linear dependency
	// now, SOS() is well-defined
	flexible_dofs = NModes()*GetNDofPerLocalMode();

	UO() << "sos_full=" << sosfull << "\n";
	UO() << "nb=" << nb << ", nf=" << nf << ", nr=" << nr << "\n";

	//+++++++++++++++++++++++
	// Compute classic Phi_CB matrix
	Matrix Phi_CB_classic;
	// check if eigenmodes can be read from file
	// if possible: Phi_CB_classic contains eigenmodes, read_eigenmodes_from_file = 1
	// if NOT possible: Phi_CB_classic is empty, read_eigenmodes_from_file = 0
	ReadEigenModesFile(Phi_CB_classic, read_eigenmodes_from_file);

	// if eigenmodes could not be read from file, they are computed
	if (read_eigenmodes_from_file == 0)
	{
		// Compute mass and stiffness matrix of unreduced system
		ComputeMassAndStiffnessMatrix(Msparse, Ksparse);

		BuildPhiCB_CMS(Phi_CB_classic, nb_modesclassic, nb_modesstatic, sosfull, nb_dofstotal, nb_dofsstatic);
	}


	//+++++++++++++++++++++++
	// Compute perpendicular Phi_CB matrix for GCMS method
	//+++++++++++++++++++++++

	Phi_CB.SetSize(sosfull, SOS());
	Phi_CB.SetAll(0.);

	Matrix3D A1, A2, A3;
	A1 = GetA1(); A2 = GetA2(); A3 = GetA3();

	// -------------
	// 12 rigid body modes: 3 translation and 9 rotation modes
	for( int node=1; node<=NCMSNodes(); node++)
	{
		Node& nodei = GetNode(node);
		Vector3D pos = nodei.Pos();
		Vector3D refpos = GetNode(RefNode1()).Pos();

		Vector3D v2 = A2*(pos - refpos);
		Vector3D v3 = A3*(pos - refpos);

		for (int i=1; i<=3; i++)
		{
			// translation: dofs 1,2,3 + rotations: 4,5
			Phi_CB(nodei.Get(i),i) = 1;
			Phi_CB(nodei.Get(i),4) = v2(i);
			Phi_CB(nodei.Get(i),5) = v3(i);
		}
	}

	// ------------
	// flexible modes: 3 modes for each local mode
	for (int locmode=1; locmode<=NModes(); locmode++)
	{
		for (int node=1; node<=NCMSNodes(); node++)
		{
			Node& nodei = GetNode(node);
			for (int i=1; i<=3; i++)
			{
				Phi_CB(nodei.Get(i),SOSRigid()+(locmode-1)*GetNDofPerLocalMode()+1) 
					= Phi_CB_classic(nodei.Get(1),locmode)*A1(i,1) + Phi_CB_classic(nodei.Get(2),locmode)*A1(i,2) + Phi_CB_classic(nodei.Get(3),locmode)*A1(i,3);
				Phi_CB(nodei.Get(i),SOSRigid()+(locmode-1)*GetNDofPerLocalMode()+2) 
					= Phi_CB_classic(nodei.Get(1),locmode)*A2(i,1) + Phi_CB_classic(nodei.Get(2),locmode)*A2(i,2) + Phi_CB_classic(nodei.Get(3),locmode)*A2(i,3);
				Phi_CB(nodei.Get(i),SOSRigid()+(locmode-1)*GetNDofPerLocalMode()+3) 
					= Phi_CB_classic(nodei.Get(1),locmode)*A3(i,1) + Phi_CB_classic(nodei.Get(2),locmode)*A3(i,2) + Phi_CB_classic(nodei.Get(3),locmode)*A3(i,3);
			}
		}
	}

	// ------------
	// find linear independent Phi_CB matrix
	// reset int flexible_dofs accordingly
	// this is not done in case the linear dependent flexible modes are constrained "UseFlexibleModesConstraints()"!!

	// old Phi_CB matrix containing linear dependent modes
	PhiCBallmodes = Phi_CB;
	project2indep.SetSize(SOSRigid(), SOSRigid());
	project2indep.SetAll(0.);
	for (int m = 1; m <= SOSRigid(); m++)
	{
		project2indep(m,m) = 1.;
	}

	if (1 && NModes() > 0 && !UseFlexibleModesConstraints())
	{
		// linear independent modes shall be computed:
		// 1) Matrix PhiTPhi contains the inner products of modes, where only each third mode is used 
		//    (consecutive three modes contain the same information in x, y, z components)
		// 2) Eigenvalues and modes of PhiTPhi are computed
		// 3) linear independent modes are linear combination of old Phi_CB according to Eigenvectors, if Eigenvalue > 1e-4
		//    otherwise, modes are dependent and not used
		// the first SOSRigid() modes are not included to get more efficient

		// 1) Compute PhiTPhi
		int n_dependent_modes = FlexDOF();
		int n_independent_modes = 0;
		int offset_rigid = SOSRigid();
		Matrix PhiTPhi(n_dependent_modes, n_dependent_modes);
		PhiTPhi.SetAll(0);

		for (int i=1; i<=n_dependent_modes; i++)
		{
			for (int j=1; j<=n_dependent_modes; j++)
			{
				for (int m=1; m<=sosfull; m++)
				{
					PhiTPhi(i,j) += PhiCBallmodes(m,i+offset_rigid)*PhiCBallmodes(m,j+offset_rigid);
				}
			}
		}

		// compute norms of mode vectors in PhiCB -> this are the sqrt of diagonal entries of PhiTPhi
		//  and scale PhiTPhi accordingly
		Vector norm_PhiCBallmodes(n_dependent_modes);
		for (int i=1; i<=n_dependent_modes; i++)
		{
			norm_PhiCBallmodes(i) = sqrt(PhiTPhi(i,i));
			// scale..
			for (int j=1; j<=n_dependent_modes; j++)
			{
				PhiTPhi(i,j) /= norm_PhiCBallmodes(i);
				PhiTPhi(j,i) /= norm_PhiCBallmodes(i);
			}
		}

		// compute eigenmodes of PhiTPhi
		Vector lami(n_dependent_modes);
		Vector work;
		work.SetLen(4*lami.Length());
		int info = LapackEVPSPD(lami.Length(), &(PhiTPhi(1,1)), &(lami(1)), &(work(1)), work.Length());
		//mbs->UO() << "aha, eigenvalues = " << lami << "\n";

		flexible_dofs = -SOSRigid()+offset_rigid;
		
		for (int i=lami.Length(); i>0; i--)
		{
			if (lami(i) >= 1e-5)
			{
				flexible_dofs++;
				n_independent_modes++;
				mbs->UO() << "* ";
			}
			mbs->UO() << "lambda(" << i << ") = " << lami(i) << "\n";
		}
		
		mbs->UO() << "aha, flexdof = " << flexible_dofs << "\n";

		if (n_independent_modes == n_dependent_modes) // do not project, keep original modes
		{
			Phi_CB = PhiCBallmodes;
		}
		else // project to linear independent modes, write matrix project2indep
		{
			Phi_CB.SetSize(sosfull,SOS());
			Phi_CB.SetAll(0);	
			for (int m=1; m<=sosfull; m++)
			{
				for (int i=1; i<=offset_rigid; i++)
				{
					Phi_CB(m,i) = PhiCBallmodes(m,i);
				}
				for (int i=1; i<=n_independent_modes; i++)
				{
					// use the i-th last mode (i.e. the i-th largest mode) with number i_mode
					int i_mode = n_dependent_modes+1-i;
					for (int j=1; j<=n_dependent_modes; j++)
					{
						Phi_CB(m,i+offset_rigid) += 1./norm_PhiCBallmodes(j)*PhiCBallmodes(m,j+offset_rigid)*PhiTPhi(i_mode,j);
					}
				}
			}
			// save projection from linear dependent to independent modes
			project2indep.SetSize(SOS()+n_dependent_modes-n_independent_modes, SOS());
			project2indep.SetAll(0.);
			for (int m = 1; m <= offset_rigid; m++)
			{
				project2indep(m,m) = 1.;
			}
			for (int m = 1; m <= n_dependent_modes; m++)
			{
				for (int i = 1; i <= n_independent_modes; i++)
				{
					// use the i-th last mode (i.e. the i-th largest mode) with number i_mode
					int i_mode = n_dependent_modes+1-i;
					project2indep(offset_rigid+m, offset_rigid+i) = 1./norm_PhiCBallmodes(m)*PhiTPhi(i_mode,m);
				}
			}
		}

		/*mbs->UO() << Phi_CB-PhiCBallmodes*project2indep << "\n";
		mbs->UO() << project2indep << "\n";*/
	} // end if (NModes() > 0)

	ofstream Phi_CB_out("phiCB.txt");
	Phi_CB_out << Phi_CB_classic << "\n";
	Phi_CB_out << "\n\n\n" << Phi_CB << "\n";
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// build K_r and M_r

	// case: sparse unreduced mass and stiffness were built
	if (Msparse.Getrows() > 0)
	{
		mbs->UO() << "build Mr...\n";
		double ts = -GetClockTime();
		Matrix mmid;
		Mult(Msparse,Phi_CB,mmid);
		MultTp(Phi_CB,mmid,Mr);

		mbs->UO() << "build Kr...\n";
		Mult(Ksparse,PhiCBallmodes,mmid);
		MultTp(PhiCBallmodes,mmid,Kr);

		ts += GetClockTime();
		mbs->UO() << "done in" << ts << "seconds\n";
	}
	// case: unreduced mass and stiffness are NOT available
	else
	{
		int sos_dependent = PhiCBallmodes.Getcols();
		Kr.SetSize(sos_dependent, sos_dependent);
		Kr.SetAll(0);		
		Mr.SetSize(SOS(), SOS());
		Mr.SetAll(0);
		ConstMatrix<FEmaxDOF*CMSmaxDOF> temp_PhiKM, temp_Phi, temp_PhiK;
		ConstMatrix<CMSmaxDOF*CMSmaxDOF> addKM;
		ConstMatrix<FEmaxDOF*FEmaxDOF> temp_K, temp_M;

		UO() << "Compute reduced mass and stiffness matrix in DoModalAnalysis:\n";
		for (int i=1; i<=NFFRFElements(); i++) 
		{
			if (NFFRFElements() > 500 && i%100 == 0) 
			{
				UO() << "\rElement " << i << " of " << NFFRFElements() << " ...";
			}

			Element& e = GetFFRFElement(i);
			int sos_el = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			temp_K.SetSize(sos_el, sos_el);
			temp_K.SetAll(0);
			temp_M.SetSize(sos_el, sos_el);
			temp_M.SetAll(0);
			// stiffness matrix
			// element stiffness matrix is negative semidef., multiply with -1 to get positive semidefinite matrix
			e.StiffnessMatrix(temp_K);
			temp_K *= -1;
			e.EvalMff(temp_M,0);

			// Kr = Phi_CB^T K_el Phi_CB

			// temp_Phi contains ltg-map of Phi_CB
			temp_Phi.SetSize(sos_el, SOS());
			temp_PhiK.SetSize(sos_el, sos_dependent);
			for (int j=1; j<=sos_el; j++)
			{
				for (int k=1; k<=SOS(); k++)
				{
					temp_Phi(j,k) = Phi_CB(ltg(j),k);
				}
				for (int k=1; k<=sos_dependent; k++)
				{
					temp_PhiK(j,k) = PhiCBallmodes(ltg(j),k);
				}
			}
			// ********** fast, Phi^T K Phi not exactly symmetric!
			Mult(temp_K, temp_PhiK, temp_PhiKM);
			MultTp(temp_PhiK,temp_PhiKM, addKM);
			////addKM.MakeSymmetric();
			Kr += addKM;
			Mult(temp_M, temp_Phi, temp_PhiKM);
			MultTp(temp_Phi,temp_PhiKM, addKM);
			Mr += addKM;

			// *********** very slow, but leads to really symmetric matrices
			//for (int ii=1; ii<=SOS(); ii++)
			//{
			//	for (int jj=1; jj<=ii; jj++)
			//	{
			//		for (int kk=1; kk<=sos_el; kk++)
			//		{
			//			for (int ll=1; ll<=sos_el; ll++)
			//			{
			//				Kr(ii,jj) += temp_Phi(kk,ii) * temp_Phi(ll,jj) * temp_K(kk,ll);
			//				Mr(ii,jj) += temp_Phi(kk,ii) * temp_Phi(ll,jj) * temp_M(kk,ll);
			//			}
			//		}
			//	}
			//}
			//for (int ii=1; ii<=SOS(); ii++)
			//{
			//	for (int jj=ii+1; jj<=SOS(); jj++)
			//	{
			//		Kr(ii,jj) = Kr(jj,ii);
			//		Mr(ii,jj) = Mr(jj,ii);
			//	}
			//}


		}

 
		// ------- Add additional masses at nodes to mass matrix
		for (int i=1; i<=massnodes.Length(); i++)
		{
			Node& node = GetNode(massnodes(i));
			for (int k=1; k<=node.SOS(); k++)
			{
				int nodedofk = node.Get(k);
				for (int j=1; j<=SOS(); j++)
				{
					for (int l=1; l<=SOS(); l++)
					{
						Mr(j,l) += massvalues(i) * Phi_CB(nodedofk,j) * Phi_CB(nodedofk,l);
					}
				}
			}
		}


	}


	//////int n_nonzero = 0;
	//////for (int i=1; i<=Kr.Getrows(); i++)
	//////{
	//////	for (int j=1; j<=Kr.Getcols(); j++)
	//////	{
	//////		if (fabs(Kr(i,j) ) > 1e-5)
	//////			n_nonzero ++;
	//////		else
	//////			Kr(i,j) = 0;
	//////	}
	//////}

	////////mbs->UO() << "Number of nonzeros = " << n_nonzero << "\n";
	//mbs->UO() << "Kr = " << Kr << "\n";
	Matrix test = Mr;
	int rv2 = test.Invert2();
	if (!rv2) {UO() << "ERROR: reduced CMS-Mass-matrix not invertable!!!\n";}

	test = Kr;
	rv2 = test.Invert2();
	if (!rv2) {UO() << "ERROR: reduced CMS-Stiffness-matrix not invertable!!!\n";}


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set new initial vector
	Vector xicopy = x_init; //length = 2*SOSRigid+1
	x_init.SetLen(2*SOS()+IS());
	x_init.SetAll(0.);

	for (int i=1; i<=SOSRigid(); i++)
	{
		// rigid body position:
		x_init(i) = xicopy(i);
		// rigid body velocity
		x_init(i+SOS()) = xicopy(SOSRigid()+i);
	}
	for (int i=1; i<=IS(); i++)
	{
		x_init(2*SOS()+i) = xicopy(2*SOSRigid()+i);
	}

	//ofstream xinitfile("x_init.txt");
	//xinitfile.precision(20);
	//xinitfile << "X_init = \n";
	//for (int i=1; i<=x_init.Length(); i++)
	//	xinitfile << x_init(i) << "\n";
	//UO() << "x_init=" << x_init << "\n";


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//intialize mass and volume
	ComputeMass();

	double totalvolume = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		totalvolume += e.GetMass()/e.Rho();
	}
	volume = totalvolume;

	UO() << "total mass=" << mass << "\n";
	UO() << "total volume=" << volume << "\n";
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetA1()
{
	Matrix3D m;
	m.SetSkew(rotation_axis);
	return Matrix3D(1.) + m*m;
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetA2()
{
	Matrix3D m;
	m.SetSkew(rotation_axis);
	return (-1)*m*m;
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetA3()
{
	Matrix3D m;
	m.SetSkew(rotation_axis);
	return m;
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetRotMatrix() const
{
	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgc(i) = XG(i);
	}
	return GetRotMatrix(xgc);
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetRotMatrix(const Vector& xgc) const
{
	Vector3D r1, r2;
	Vector3D r01, r02;
	r1 = GetNodePosRBModes(RefNode1(), xgc);
	r2 = GetNodePosRBModes(RefNode2(), xgc);
	//r3 = GetNodePosRBModes(RefNode3(), xgc);
	r01 = GetNode(RefNode1()).Pos();
	r02 = GetNode(RefNode2()).Pos();
	//r03 = GetNode(RefNode3()).Pos();
	Vector3D a, b, c, a0, b0, c0;
	a = rotation_axis;
	a.Normalize();
	b = (r2-r1);
	// GramSchmidt funktion  -->  b -= (a*b) * a;	b.Normalize();
	a.GramSchmidt(b);
	c = a.Cross(b);
	a0 = rotation_axis;
	a0.Normalize();
	b0 = (r02-r01);
	a0.GramSchmidt(b0);
	c0 = a0.Cross(b0);
	Matrix3D B(a(1), b(1), c(1),
						 a(2), b(2), c(2),
						 a(3), b(3), c(3));
	Matrix3D B0T(a0(1), a0(2), a0(3),
							b0(1), b0(2), b0(3),
							c0(1), c0(2), c0(3));
	return B*B0T;
}

template <class RIGID>
Matrix3D GCMSRotorElement<RIGID>::GetAbar(const Matrix3D& A) const
{
	Vector3D r01, r02;
	Vector3D a0, b0, c0;
	r01 = GetNode(RefNode1()).Pos();
	r02 = GetNode(RefNode2()).Pos();
	a0 = rotation_axis;
	a0.Normalize();
	b0 = (r02-r01);
	a0.GramSchmidt(b0);
	c0 = a0.Cross(b0);
	Matrix3D PTp(a0(1), a0(2), a0(3), b0(1), b0(2), b0(3), c0(1), c0(2), c0(3));
	Matrix3D Abar = PTp*A*PTp.GetTp();
	return Abar;
}

template <class RIGID>
void GCMSRotorElement<RIGID>::GetSinCosPhi(const Matrix3D& A, double* sinphi, double* cosphi) const
{
	Matrix3D Abar = GetAbar(A);
	*sinphi = Abar(3,2);
	*cosphi = Abar(2,2); 
	//return acos(Abar(2,2));
}

template <class RIGID>
void GCMSRotorElement<RIGID>::ApplyRotation(const Matrix3D& A, Vector& u) const
{
	Matrix3D Abar = GetAbar(A);

	// translational dofs
	int j = 1;
	Vector3D uloc2(u(j), u(j+1), u(j+2));
	j += 3;
	Vector3D uloc;
	uloc = A*uloc2;
	u(j-3) = uloc(1);
	u(j-2) = uloc(2);
	u(j-1) = uloc(3);

	// rotational dofs 
	Vector2D uloc_rot(u(4), u(5));
	u(4) = Abar(2,2) * uloc_rot(1) + Abar(2,3) * uloc_rot(2);
	u(5) = Abar(3,2) * uloc_rot(1) + Abar(3,3) * uloc_rot(2);

	// flexible dofs
	j = 6;
	while (j <= SOS())
	{
		Vector3D uloc2(u(j), u(j+1), u(j+2));
		j += 3;
		Vector3D uloc;
		uloc = Abar*uloc2;
		u(j-3) = uloc(1);
		u(j-2) = uloc(2);
		u(j-1) = uloc(3);
	}
}


// can be set to 2 if the stiffness matrix routine works fine 
template <class RIGID>
int GCMSRotorElement<RIGID>::FastStiffnessMatrix() const { if (Phi_CB.Getcols()==PhiCBallmodes.Getcols()) {return 2;} else {return 0; } }

template <class RIGID>
void GCMSRotorElement<RIGID>::ApplyRotationFromRight(const Matrix3D& A, Matrix& K) const
{
	Matrix3D AT = A.GetTp();
	Matrix3D Abar = GetAbar(A);
	Matrix3D AbarT = Abar.GetTp();

	for (int i = 1; i <= SOS(); i++)
	{
		// translational dofs
		int j = 1;
		Vector3D uloc2(K(i,j), K(i,j+1), K(i,j+2));
		j += 3;
		Vector3D uloc;
		uloc = AT*uloc2;
		K(i,j-3) = uloc(1);
		K(i,j-2) = uloc(2);
		K(i,j-1) = uloc(3);

		// rotational dofs 
		Vector2D uloc_rot(K(i,4), K(i,5));
		K(i,4) = AbarT(2,2) * uloc_rot(1) + AbarT(2,3) * uloc_rot(2);
		K(i,5) = AbarT(3,2) * uloc_rot(1) + AbarT(3,3) * uloc_rot(2);

		// flexible dofs
		j = 6;
		while (j <= SOS())
		{
			Vector3D uloc2(K(i,j), K(i,j+1), K(i,j+2));
			j += 3;
			Vector3D uloc;
			uloc = AbarT*uloc2;
			K(i,j-3) = uloc(1);
			K(i,j-2) = uloc(2);
			K(i,j-1) = uloc(3);	
		}
	}
}

template <class RIGID>
void GCMSRotorElement<RIGID>::ApplyRotationFromLeft(const Matrix3D& A, Matrix& K) const
{
	Matrix3D Abar = GetAbar(A);
	for (int i = 1; i <= SOS(); i++)
	{
		// translational dofs
		int j = 1;
		Vector3D uloc2(K(j,i), K(j+1,i), K(j+2,i));
		j += 3;
		Vector3D uloc;
		uloc = A*uloc2;
		K(j-3,i) = uloc(1);
		K(j-2,i) = uloc(2);
		K(j-1,i) = uloc(3);

		// rotational dofs 
		Vector2D uloc_rot(K(4,i), K(5,i));
		K(4,i) = Abar(2,2) * uloc_rot(1) + Abar(2,3) * uloc_rot(2);
		K(5,i) = Abar(3,2) * uloc_rot(1) + Abar(3,3) * uloc_rot(2);

		// flexible dofs
		j = 6;
		while (j <= SOS())
		{
			Vector3D uloc2(K(j,i), K(j+1,i), K(j+2,i));
			j += 3;
			Vector3D uloc;
			uloc = Abar*uloc2;
			K(j-3,i) = uloc(1);
			K(j-2,i) = uloc(2);
			K(j-1,i) = uloc(3);	
		}
	}
}

template <class RIGID>
void GCMSRotorElement<RIGID>::ApplyDADAkl(int k, int l, Vector& u) const
{
	Matrix3D A(0.);
	A(k,l) = 1.;

	// translational dofs
	int j = 1;
	// u(k) = u(l);
	double ul = u(l);
	for (int i=1; i<=3; i++)
	{
		if (i==k) u(i) = ul;
		else u(i) = 0;
	}
	j += 3;

	// rotational dofs 
	Vector2D uloc_rot(u(4), u(5));
	u(4) = A(2,2) * uloc_rot(1) + A(2,3) * uloc_rot(2);
	u(5) = A(3,2) * uloc_rot(1) + A(3,3) * uloc_rot(2);

	// flexible dofs
	j = 6;
	while (j <= SOS())
	{
		//Vector3D uloc2(u(j), u(j+1), u(j+2));
		//Vector3D uloc;
		//uloc = A*uloc2;
		//u(j-3) = uloc(1);
		//u(j-2) = uloc(2);
		//u(j-1) = uloc(3);
		// u(j-1+k) = u(j-1+l)
		double ul = u(j-1+l);
		for (int i=1; i<=3; i++)
		{
			if (i==k) u(j-1+i) = ul;
			else u(j-1+i) = 0;
		}
		j += 3;
	}
}

template <class RIGID>
void GCMSRotorElement<RIGID>::ComputeURigid(const Matrix3D& A, const Vector3D& uref1, Vector& urigid) const
{
	urigid.SetLen(SOS());
	urigid.SetAll(0.);
	// translation dofs: 1, 2, 3
	for (int i=1; i<=3; i++)
	{
		urigid(i) = uref1(i);
	}
	double sinphi = 0.;
	double cosphi = 0.;
	GetSinCosPhi(A, &sinphi, &cosphi);

	// rotation dofs
	urigid(4) = cosphi-1.;
	urigid(5) = sinphi;
}

// compute derivative of rotation matrix entries with respect to all reduced degrees of freedom
template <class RIGID>
void GCMSRotorElement<RIGID>::ComputeDADq(ConstVector<CMSmaxDOF> DAijDq[3][3]) const
{
	//mbs->UO() << "DADq not implemented yet!!!" << "\n";
	Matrix3D DAijDqk;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			DAijDq[i][j].SetLen(SOS());
			DAijDq[i][j].SetAll(0);
		}
	}
	for (int k=1; k<=SOS(); k++)
	{
		ComputeDADqk_NumericDiff(k, DAijDqk);
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				DAijDq[i][j](k) = DAijDqk(i+1,j+1);
			}
		}
	}

}


// compute derivative of rotation matrix entries with respect to reduced degree of freedom q_k
// using numerical differentiation
template <class RIGID>
void GCMSRotorElement<RIGID>::ComputeDADqk_NumericDiff(int k, Matrix3D& dAdq) const
{
	// otherwise - numerical differentiation
	// A depends only on rotation parameters dof 4 - 5
		if (k<4 || k>5)
		{
			dAdq.SetAll(0.);
			return;
		}
	double eps = 1e-10;
	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgc(i) = XG(i);
	}
	//Matrix3D dAdqtest;
	xgc(k) -= eps;
	Matrix3D A0 = GetRotMatrix(xgc);
	xgc(k) += 2*eps;
	Matrix3D Aeps = GetRotMatrix(xgc);
	xgc(k) -= eps;
	dAdq = Aeps - A0;
	dAdq *= 1./(2.*eps);
	return;

}

template <class RIGID>
void GCMSRotorElement<RIGID>::EvalF2(Vector& f, double t)
{
	Body3D::EvalF2(f,t);
	TMStartTimer(22);

	xg.SetLen(SOS());
	for (int i=1; i <= SOS(); i++) xg(i) = XG(i);

	ConstVector<CMSmaxDOF> temp; 
	ConstVector<CMSmaxDOF> qF;
	ConstVector<CMSmaxDOF> ATqF;
	ConstVector<CMSmaxDOF> urigid;
	ConstVector<CMSmaxDOF> KATqF;
	ConstVector<CMSmaxDOF> AKATqF;
	ConstVector<CMSmaxDOF> dAijdq;
	ConstVector<CMSmaxDOF> projectDADAklKATqF;
	
	// rigid body motion
	Matrix3D A = GetRotMatrix(xg);
	Matrix3D AT = A.GetTp();
	Vector3D uref1 = GetTranslation(xg);

	//+++++++++++++++++++++++++++
	// corotated stiffness matrix
	// f -= diag(A) K diag(A^T) (q-qR)

	TMStartTimer(24);
	// compute qR, qF = q-qR
	ComputeURigid(A,uref1,urigid);
	qF = xg;
	qF -= urigid;

	// compute diag(A^T)*project*(q-qR)
	// check if projection to independent modes necessary
	int do_project = 0;
	if (Phi_CB.Getcols() != PhiCBallmodes.Getcols()) {do_project = 1;}

	if (do_project) // do not use 
	{
		Mult(project2indep, qF, ATqF);
	}
	else
	{
		ATqF = qF;
	}
	ApplyRotation(AT,ATqF);

	// compute K diag(A^T)*project*(q-qR)
	Mult(Kr,ATqF,KATqF);

	// compute A K diag(A^T)*project*(q-qR)
	AKATqF = KATqF;
	ApplyRotation(A,AKATqF);
	// compute temp = project*diag(A) K diag(A^T)*project*(q-qR)
	if (do_project)
	{
		MultTp(project2indep, AKATqF, temp);
	}
	else
	{
		temp = AKATqF;
	}
	
	TMStopTimer(24);

	for (int i=1; i <= SOS(); i++) f(i) -= temp(i);


	//+++++++++++++++++++++++++++
	//compute f_NL:
	// f(k) -= dA_ij/dq_k [(q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))]

	if (1 && !GetMBS()->IsJacobianComputation())
	{
		Matrix3D Aij;
		double DfnlDAij;

		TMStartTimer(25);
		// compute DAij/dqk
		ConstVector<CMSmaxDOF> DAijDq[3][3];
		//Matrix3D DAijDqk;
		//for (int i=0; i<3; i++)
		//{
		//	for (int j=0; j<3; j++)
		//	{
		//		DAijDq[i][j].SetLen(SOS());
		//		DAijDq[i][j].SetAll(0);
		//	}
		//}
		//for (int k=1; k<=SOS(); k++)
		//{
		//	ComputeDADqk_NumericDiff(k, DAijDqk);
		//	for (int i=0; i<3; i++)
		//	{
		//		for (int j=0; j<3; j++)
		//		{
		//			DAijDq[i][j](k) = DAijDqk(i+1,j+1);
		//		}
		//	}
		//}
		ComputeDADq(DAijDq);
		TMStopTimer(25);
		TMStartTimer(26);
		
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= Dim(); j++)
			{
				/*Aij.SetAll(0);
				Aij(i,j) = 1;*/

				// compute temp = diag(dA/dA_ij) K diag(AT) (q-qR)
				temp = KATqF;
				ApplyDADAkl(i,j,temp);
				if (do_project)
				{
				MultTp(project2indep, temp, projectDADAklKATqF);
				}
				else
				{
					projectDADAklKATqF = temp;
				}
				// compute (q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))
				DfnlDAij = qF*projectDADAklKATqF;
				// subtract DfnlDAij * DAij/Dq = dA_ij/dq_k [(q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))]
				f.MultAdd(-DfnlDAij,DAijDq[i-1][j-1]);

			}
		}
		TMStopTimer(26);
	}
}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position / velocity of node in the GCMSRotor element
template <class RIGID>
Vector3D GCMSRotorElement<RIGID>::GetNodePos(int i, const Vector& xgc) const
{
	int sos = SOS();
	const Node& nodei = GetNode(i);
	Vector3D pglob = nodei.Pos();
	for (int j=1; j<=sos; j++)
	{
		for (int k=1; k<=Dim(); k++)
		{
			int nodedofk = nodei.Get(k);
			pglob(k) += Phi_CB(nodedofk, j) * (xgc(j));
		}
	}
	return pglob;
}

template <class RIGID>
void GCMSRotorElement<RIGID>::GetNodedPosdqT(int node, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();
	const Node& nodei = GetNode(node);
	
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=SOS(); j++)
		{
			dpdqi(j,k) += Phi_CB(nodedofk, j);
		}
	}
	return;
}

template <class RIGID>
Vector3D GCMSRotorElement<RIGID>::GetNodeVel(int i, const Vector& xgp) const
{
	int sos = SOS();
	const Node& nodei = GetNode(i);
	Vector3D vel(0.);
	for (int k=1; k<=Dim(); k++)
	{
		for (int j=1; j<=sos; j++)
		{
			int nodedofk = nodei.Get(k);
			vel(k) += Phi_CB(nodedofk, j) * xgp(j);
		}
	}	
	return vel;
}

template class GCMSRotorElement<Rigid3D>;
template class GCMSRotorElement<Rigid3DKardan>;