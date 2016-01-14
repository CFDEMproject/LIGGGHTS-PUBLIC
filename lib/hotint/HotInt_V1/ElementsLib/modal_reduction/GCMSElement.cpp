//#**************************************************************
//#
//# filename:             GCMSElement.cpp
//#
//# author:               Astrid Pechstein
//#
//# generated:						Okt 2010
//# description:          Element for generalized component mode synthesis with perpendicular modes
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
 

//#include "windows.h" //for shellexecute
#include "element.h"
#include "material.h"
#include "Node.h"
#include "rigid3D.h"
#include "femathhelperfunctions.h"
#include "rigid3Dkardan.h"
#include "femeshinterface.h"
#include "kinematicpairs.h"
#include "referenceFrame3D.h"
#include "BaseCMSElement.h"
#include "gcmselement.h"
//#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"
#include "lapack_routines.h"
#include "lapack_routines.h"
#include "finiteelement3d.h"
#include "myfile.h"

template <class RIGID>
void GCMSElement<RIGID>::SetUseRigidBodyConstraints(int use_rbc)
{
	use_rb_constraints = use_rbc;
	// in case of rigid body constraints, x_init vector has to be enlargened by IS() = 6
	if (use_rb_constraints)
	{
		Vector xicopy = x_init; 
		x_init.SetLen(2*SOS()+IS());
		x_init.SetAll(0.);

		for (int i=1; i<=xicopy.Length(); i++)
		{
			x_init(i) = xicopy(i);
		}
	}
	else // in case of no rigid body constraints, x_init vector is of size 2*SOS()
	{
		x_init.SetLen(2*SOS()+IS());
	}
}

template <class RIGID>
void GCMSElement<RIGID>::SetUseFlexibleModeConstraints(int use_fmc)
{
	use_mode_constraints = use_fmc;
	// in case of rigid body constraints, x_init vector has to be enlargened by IS()
	if (use_mode_constraints)
	{
		Vector xicopy = x_init; 
		x_init.SetLen(2*SOS()+IS());
		x_init.SetAll(0.);

		for (int i=1; i<=xicopy.Length(); i++)
		{
			x_init(i) = xicopy(i);
		}
	}
	else // in case of no rigid body constraints, x_init vector is of size 2*SOS()
	{
		x_init.SetLen(2*SOS()+IS());
	}
}


// use point reference point to define initial condition phip
// set number of internal modes
// compute initial rigid body degrees of freedom,
// set size and color vectors
template <class RIGID>
void GCMSElement<RIGID>::SetGCMSElement(const Vector3D& p, const Vector3D& v_refP, Vector3D phi, Vector3D phip_refP, const Vector3D& ref_node1TOref_pos, const int nimodesi,
		const Vector3D& sizei, const Vector3D& coli)
{
	Vector3D v = v_refP - phip_refP.Cross(ref_node1TOref_pos);
	SetGCMSElement(p, v, phi, phip_refP, nimodesi, sizei, coli);
}


// set number of internal modes
// compute initial rigid body degrees of freedom,
// set size and color vectors
template <class RIGID>
void GCMSElement<RIGID>::SetGCMSElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
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

	dofs_per_local_mode = 9;

	// compute rigid dofs for p, v, phi, phip
	// other dofs will be appended to x_init later on
	x_init.SetLen(2*SOSRigid()+IS());
	x_init.SetAll(0.);
	Matrix3D A = ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z());


	//Matrix3D Ap = ComputeRotMatrixEulerP(phi.X(),phi.Y(), phi.Z(), phip.X(), phip.Y(), phip.Z()); //$ DR 201-03: old code

	// ====================================================
	//$ DR 2011-04-01:[ ComputeRotMatrixPfromOmegaVector (A,omega=phip) --> gehört in linalg3d
	double b0, b1, b2, b3;
	//RotMatToQuaternions(ComputeRotMatrixEuler(phi.X(),phi.Y(),phi.Z()),b0, b1, b2, b3);
	RotMatToQuaternions(A,b0, b1, b2, b3);

	//initial values for time-wise derivative quaternions beta_t:
	double beta0 = 2*b0; double beta1 = 2*b1; double beta2 = 2*b2; double beta3 = 2*b3; 
	//the betas are already multiplied with 2, compared to Shabana
	Matrix G(4,4);
	G(1,1) = -beta1; G(1,2) =  beta0; G(1,3) = -beta3; G(1,4) =  beta2;
	G(2,1) = -beta2; G(2,2) =  beta3; G(2,3) =  beta0; G(2,4) = -beta1;
	G(3,1) = -beta3; G(3,2) = -beta2; G(3,3) =  beta1; G(3,4) =  beta0;
	G(4,1) =  beta0; G(4,2) =  beta1; G(4,3) =  beta2; G(4,4) =  beta3; //this is the time-derivative of the quaternion condition | |^2 = 1

	Vector f(phip.X(),phip.Y(),phip.Z(),0.); //omega and zero
	Vector betap(4); 
	int rv = G.Solve(f, betap);
	if (!rv) {GetMBS()->UO() << "GCMS:Initialization: could not determine initial Euler parameter velocities due to singularity!!!\n";}

	Matrix3D Ap=ComputeRotMatrixPEulerParam(b0,b1,b2,b3,betap(1),betap(2),betap(3),betap(4));

	//$ DR 2011-04-01:]
	// ====================================================

	for (int i=1; i<=3; i++)
	{
		x_init(i) = p(i);
	}
	// rotation dofs
	x_init(4) = A(1,1)-1.;
	x_init(5) = A(2,1);
	x_init(6) = A(3,1);
	x_init(7) = A(1,2);
	x_init(8) = A(2,2)-1.;
	x_init(9) = A(3,2);
	x_init(10) = A(1,3);
	x_init(11) = A(2,3);
	x_init(12) = A(3,3)-1.;
	// translation velocity dofs
	for (int i=1; i<=3; i++)
	{
		x_init(SOSRigid()+i) = v(i);
	}
	// compute rotational velocity dofs
	x_init(SOSRigid()+4) = Ap(1,1);
	x_init(SOSRigid()+5) = Ap(2,1);
	x_init(SOSRigid()+6) = Ap(3,1);
	x_init(SOSRigid()+7) = Ap(1,2);
	x_init(SOSRigid()+8) = Ap(2,2);
	x_init(SOSRigid()+9) = Ap(3,2);
	x_init(SOSRigid()+10) = Ap(1,3);
	x_init(SOSRigid()+11) = Ap(2,3);
	x_init(SOSRigid()+12) = Ap(3,3);

	GetMBS()->UO() << "x_init = " << x_init << "\n";
	//x_init.SetAll(0);
	//x_init(SOSRigid()+5)=-2*1000*MY_PI;
	//x_init(SOSRigid()+7)=2*1000*MY_PI;
	//x_init.SetRandom();

	// AH: testing purpose
	/*double Iphi1 = 3235.7967;
	double Iphi2 = 3235.7967;
	double Iphi3 = 4103.2571;
	Iphi(1,1) = Iphi1;
	Iphi(2,2) = Iphi2;
	Iphi(3,3) = Iphi3;*/
}


// do nothing, since constant matrices are precomputed in DoModalAnalysis
template <class RIGID>
void GCMSElement<RIGID>::Initialize()
{
	;
}



	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//must be called before Assemble!!!!!
//link elements, 
//ltg-elements, 
//compute nbmodes, 
//compute full M and K, 
//modal analysis, 
//transformation matrix and modal degrees of freedom
//modal matrices
template <class RIGID>
void GCMSElement<RIGID>::DoModalAnalysis(const TArray<int2>& fixednodes)
{
	UO() << "GCMSElement::DoModalAnalysis()\n";
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

	// -------------
	// 12 rigid body modes: 3 translation and 9 rotation modes
	for( int node=1; node<=NCMSNodes(); node++)
	{
		Node& nodei = GetNode(node);
		Vector3D pos = nodei.Pos();
		Vector3D refpos = GetNode(RefNode1()).Pos();
		for (int i=1; i<=3; i++)
		{
			// translation: dofs 1,2,3
			Phi_CB(nodei.Get(i),i) = 1;
			// rotation: x component dofs 4,5,6
			Phi_CB(nodei.Get(i),i+3) = pos.X()-refpos.X();
			// rotation: y component dofs 7,8,9
			Phi_CB(nodei.Get(i),i+6) = pos.Y()-refpos.Y();
			// rotation: z component dofs 10,11,12
			Phi_CB(nodei.Get(i),i+9) = pos.Z()-refpos.Z();
		}
	}

	// ------------
	// flexible modes: 9 modes for each local mode
	for (int locmode=1; locmode<=NModes(); locmode++)
	{
		for (int node=1; node<=NCMSNodes(); node++)
		{
			Node& nodei = GetNode(node);
			for (int i=1; i<=3; i++)
			{
				for (int j=1; j<=3; j++)
				{
					Phi_CB(nodei.Get(j),SOSRigid()+(locmode-1)*GetNDofPerLocalMode()+(i-1)*3+j) = Phi_CB_classic(nodei.Get(i),locmode);
				}
			}
		}
	}

	// ------------
	// find linear independent Phi_CB matrix
	// reset int flexible_dofs accordingly
	// this is not done in case the linear dependent flexible modes are constrained "UseFlexibleModesConstraints()"!!

	if (NModes() > 0 && !UseFlexibleModesConstraints())
	{
		// old Phi_CB matrix containing linear dependent modes
		Matrix PhiCBallmodes = Phi_CB;
		// linear independent modes shall be computed:
		// 1) Matrix PhiTPhi contains the inner products of modes, where only each third mode is used 
		//    (consecutive three modes contain the same information in x, y, z components)
		// 2) Eigenvalues and modes of PhiTPhi are computed
		// 3) linear independent modes are linear combination of old Phi_CB according to Eigenvectors, if Eigenvalue > 1e-4
		//    otherwise, modes are dependent and not used
		// the first SOSRigid() modes are not included to get more efficient

		// 1) Compute PhiTPhi
		int n_dependent_modes = FlexDOF()/3;
		int n_independent_modes; // in the end, FlexDOF is 3 times this number
		int offset_rigid = SOSRigid();
		Matrix PhiTPhi(n_dependent_modes, n_dependent_modes);
		PhiTPhi.SetAll(0);

		for (int i=1; i<=n_dependent_modes; i++)
		{
			for (int j=1; j<=n_dependent_modes; j++)
			{
				for (int m=1; m<=sosfull; m++)
				{
					PhiTPhi(i,j) += PhiCBallmodes(m,(i-1)*3+1+offset_rigid)*PhiCBallmodes(m,3*(j-1)+1+offset_rigid);
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
		//GetMBS()->UO() << "aha, eigenvalues = " << lami << "\n";

		flexible_dofs = -SOSRigid()+offset_rigid;
		n_independent_modes = 0;
		for (int i=lami.Length(); i>0; i--)
		{
			if (lami(i) >= 1e-7)
			{
				flexible_dofs += 3;
				n_independent_modes++;
				GetMBS()->UO() << "* ";
			}
			GetMBS()->UO() << "lambda(" << i << ") = " << lami(i) << "\n";
		}
		GetMBS()->UO() << "aha, flexdof = " << flexible_dofs << "\n";
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
					Phi_CB(m,3*(i-1)+1+offset_rigid) += 1./norm_PhiCBallmodes(j)*PhiCBallmodes(m,offset_rigid+3*(j-1)+1)*PhiTPhi(i_mode,j);
					Phi_CB(m,3*(i-1)+2+offset_rigid) += 1./norm_PhiCBallmodes(j)*PhiCBallmodes(m,offset_rigid+3*(j-1)+2)*PhiTPhi(i_mode,j);
					Phi_CB(m,3*(i-1)+3+offset_rigid) += 1./norm_PhiCBallmodes(j)*PhiCBallmodes(m,offset_rigid+3*(j-1)+3)*PhiTPhi(i_mode,j);
				}
			}
		}
	} // end if (NModes() > 0)

	ofstream Phi_CB_out("phiCB.txt");
	Phi_CB_out << Phi_CB_classic << "\n";
	Phi_CB_out << "\n\n\n" << Phi_CB << "\n";
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// build K_r and M_r

	// case: sparse unreduced mass and stiffness were built
	if (Msparse.Getrows() > 0)
	{
		GetMBS()->UO() << "build Mr...\n";
		double ts = -GetClockTime();
		Matrix mmid;
		Mult(Msparse,Phi_CB,mmid);
		MultTp(Phi_CB,mmid,Mr);

		GetMBS()->UO() << "build Kr...\n";
		Mult(Ksparse,Phi_CB,mmid);
		MultTp(Phi_CB,mmid,Kr);

		ts += GetClockTime();
		GetMBS()->UO() << "done in" << ts << "seconds\n";
	}
	// case: unreduced mass and stiffness are NOT available
	else
	{
		Kr.SetSize(SOS(), SOS());
		Kr.SetAll(0);		
		Mr.SetSize(SOS(), SOS());
		Mr.SetAll(0);
		ConstMatrix<FEmaxDOF*CMSmaxDOF> temp_PhiKM, temp_Phi;
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
			for (int j=1; j<=sos_el; j++)
			{
				for (int k=1; k<=SOS(); k++)
				{
					temp_Phi(j,k) = Phi_CB(ltg(j),k);
				}
			}
			// ********** fast, Phi^T K Phi not exactly symmetric!
			Mult(temp_K, temp_Phi, temp_PhiKM);
			MultTp(temp_Phi,temp_PhiKM, addKM);
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

	////////GetMBS()->UO() << "Number of nonzeros = " << n_nonzero << "\n";
	//GetMBS()->UO() << "Kr = " << Kr << "\n";
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

// storematrices not used -- would be very inefficient
// compute mass matrix of CMS element (constant!)
template <class RIGID>
void GCMSElement<RIGID>::EvalM(Matrix& m, double t)
{
	if (Mr.Getrows() == SOS() )
	{
		m = Mr;
	}
	else
	{
		Mr.SetSize(SOS(), SOS());
		Mr.SetAll(0);
		ConstMatrix<FEmaxDOF*CMSmaxDOF> temp_PhiM, temp_Phi;
		ConstMatrix<CMSmaxDOF*CMSmaxDOF> addM;
		ConstMatrix<FEmaxDOF*FEmaxDOF> temp_M;

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			if (NFFRFElements() > 500 && i%100 == 0) 
			{
				UO() << "\rElement " << i << " of " << NFFRFElements() << " ...";
			}

			Element& e = GetFFRFElement(i);
			int sos_el = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			temp_M.SetSize(sos_el, sos_el);
			temp_M.SetAll(0);
			e.EvalMff(temp_M,0);

			// Mr = Phi_CB^T M_el Phi_CB

			// temp_Phi contains ltg-map of Phi_CB
			temp_Phi.SetSize(sos_el, SOS());
			for (int j=1; j<=sos_el; j++)
			{
				for (int k=1; k<=SOS(); k++)
				{
					temp_Phi(j,k) = Phi_CB(ltg(j),k);
				}
			}
			Mult(temp_M, temp_Phi, temp_PhiM);
			MultTp(temp_Phi,temp_PhiM, addM);
			Mr += addM;
		}


		// ------- Add additional masses at nodes to mass matrix
		for (int i=1; i<=massnodes.Length(); i++)
		{
			Node& node = GetNode(massnodes(i));
			for (int k=1; k<=node.SOS(); k++)
			{
				for (int j=1; j<=SOS(); j++)
				{
					for (int l=1; l<=SOS(); l++)
					{
						Mr(j,l) += massvalues(i) * Phi_CB(node.Get(k),j) * Phi_CB(node.Get(k),l);
					}
				}
			}
		}
		m = Mr;
	}
}

// stiffness+quadratic velocity vector:
// Equations according to G.+Ambrosio, (34)(36)
// corotated stiffness matrix
// f -= diag(A) K diag(A^T) (q-qR)
// nonlinear terms
// f(m) -= dA_ij/dq_m [(q-qR)^T diag(dA/dA_ij) K diag(A) (q-qR^))]
template <class RIGID>
void GCMSElement<RIGID>::EvalF2(Vector& f, double t)
{
	Body3D::EvalF2(f,t);

	if (!mbs->IsJacobianComputation())
	{
		//GetMBS()->UO() << "F(" << t << ") = " << f << "\n";
	}
	TMStartTimer(22);

	xg.SetLen(SOS());
	for (int i=1; i <= SOS(); i++) xg(i) = XG(i);

	ConstVector<CMSmaxDOF> temp; 
	ConstVector<CMSmaxDOF> qF;
	ConstVector<CMSmaxDOF> ATqF;
	ConstVector<CMSmaxDOF> urigid;
	ConstVector<CMSmaxDOF> KATqF;
	ConstVector<CMSmaxDOF> dAijdq;

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

	// compute diag(A^T) (q-qR)
	ATqF = qF;
	ApplyRotation(AT,ATqF);

	// compute K diag(A^T) (q-qR)
	Mult(Kr,ATqF,KATqF);

	// compute temp = diag(A) K diag(A^T) (q-qR)
	temp = KATqF;
	ApplyRotation(A,temp);
	TMStopTimer(24);

	for (int i=1; i <= SOS(); i++) f(i) -= temp(i);


	//+++++++++++++++++++++++++++
	//compute f_NL:
	// f(k) -= dA_ij/dq_k [(q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))]

	if (/*1 || */!GetMBS()->IsJacobianComputation())
	{
		Matrix3D Aij;
		double DfnlDAij;

		TMStartTimer(25);
		// compute DAij/dqk
		ConstVector<CMSmaxDOF> DAijDq[3][3];
		Matrix3D DAijDqk;
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				DAijDq[i][j].SetLen(SOS());
				DAijDq[i][j].SetAll(0);
			}
		}
		//for (int k=1; k<=SOS(); k++)
		//{
		//	ComputeDADqk_NumericDiff(k, DAijDqk);
		//}
		ComputeDADq(DAijDq);
		TMStopTimer(25);
		TMStartTimer(26);

		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= Dim(); j++)
			{
				Aij.SetAll(0);
				Aij(i,j) = 1;

				// compute temp = diag(dA/dA_ij) K diag(AT) (q-qR)
				temp = KATqF;
				ApplyRotation(Aij,temp);
				// compute (q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))
				DfnlDAij = qF*temp;
				// subtract DfnlDAij * DAij/Dq = dA_ij/dq_k [(q-qR)^T diag(dA/dA_ij) K diag(AT) (q-qR))]
				f.MultAdd(-DfnlDAij,DAijDq[i-1][j-1]);

			}
		}
		TMStopTimer(26);
	}

	// number of constraints for rigid modes
	int offset_rigid = 0;
	if (UseRigidBodyConstraints())
	{ 
		offset_rigid = 6;
	// Constraint terms
	ConstVector<6> lambda(offset_rigid);
	for (int i=1; i<=offset_rigid; i++)
		lambda(i) = XG(SOS()*2+i);

	// 3 conditions:
	// Norm of i-th line of rotation matrix = 1, ||A_i|| = 1
	f(4) -= 2*lambda(1)*(xg(4)+1);
	f(5) -= 2*lambda(1)*(xg(5));
	f(6) -= 2*lambda(1)*(xg(6));
	f(7) -= 2*lambda(2)*(xg(7));
	f(8) -= 2*lambda(2)*(xg(8)+1);
	f(9) -= 2*lambda(2)*(xg(9));
	f(10) -= 2*lambda(3)*(xg(10));
	f(11) -= 2*lambda(3)*(xg(11));
	f(12) -= 2*lambda(3)*(xg(12)+1);
	// 3 conditions:
	// ith and jth line of rotation matrix are orthogonal, A_i . A_j = 0
	f(4) -= lambda(4)*xg(7);
	f(7) -= lambda(4)*(xg(4)+1);
	f(5) -= lambda(4)*(xg(8)+1);
	f(8) -= lambda(4)*xg(5);
	f(6) -= lambda(4)*xg(9);
	f(9) -= lambda(4)*xg(6);
	f(4) -= lambda(5)*xg(10);
	f(10) -= lambda(5)*(xg(4)+1);
	f(5) -= lambda(5)*xg(11);
	f(11) -= lambda(5)*xg(5);
	f(6) -= lambda(5)*(xg(12)+1);
	f(12) -= lambda(5)*xg(6);
	f(7) -= lambda(6)*xg(10);
	f(10) -= lambda(6)*xg(7);
	f(8) -= lambda(6)*xg(11);
	f(11) -= lambda(6)*(xg(8)+1);
	f(9) -= lambda(6)*(xg(12)+1);
	f(12) -= lambda(6)*xg(9);
	}

	// values of lambda concerning mode number mode
	ConstVector<9> lambda_mode(GetNDofPerLocalMode()-1);
	// values of degrees of freedom q concerning mode number mode
	ConstVector<9> q_mode(GetNDofPerLocalMode());
	// degrees of freedom concerning local mode number mode
	TArray<int> mode_index;
	mode_index.SetLen(GetNDofPerLocalMode());
	if (UseFlexibleModesConstraints())
	{
		for (int mode=1; mode<=NModes(); mode++)
		{
			// degrees of freedom belonging to local mode number mode
			for (int i=1; i<=GetNDofPerLocalMode()-1; i++)
			{
				mode_index(i) = SOSRigid()+(mode-1)*GetNDofPerLocalMode()+i;
				q_mode(i) = XG(mode_index(i));
				lambda_mode(i) = XG(SOS()*2+offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+i);
			}
			mode_index(GetNDofPerLocalMode()) = SOSRigid()+(mode)*GetNDofPerLocalMode();
			q_mode(GetNDofPerLocalMode()) = XG(mode_index(GetNDofPerLocalMode()));

			//if (fabs(xg(4)+1) > 0.2)
			//{
			// first constraint
			// q(1) = A_11 * qf_prime = (xg(4)+1)*qf_prime   AND  q(2) = A_21*qf_prime = xg(5)*qf_prime ---->
			// q(1)*xg(5) = q(2)*(xg(4)+1);
			// thus entries in f(mode_index(1)) due to q(1), f(mode_index(2)) due to q(2), f(5) due to xg(5), f(4) due to xg(4)
			f(mode_index(9)) -= lambda_mode(1)*xg(5);
			f(5) -= lambda_mode(1)*q_mode(9);
			f(mode_index(2)) += lambda_mode(1)*(xg(12)+1);
			f(12) += lambda_mode(1)*q_mode(2);

			// other constraints similar
			//2
			f(mode_index(9)) -= lambda_mode(2)*xg(6);
			f(6) -= lambda_mode(2)*q_mode(9);
			f(mode_index(3)) += lambda_mode(2)*(xg(12)+1);
			f(12) += lambda_mode(2)*q_mode(3);
			//3
			f(mode_index(9)) -= lambda_mode(3)*xg(7);
			f(7) -= lambda_mode(3)*q_mode(9);
			f(mode_index(4)) += lambda_mode(3)*(xg(12)+1);
			f(12) += lambda_mode(3)*q_mode(4);
			//4
			f(mode_index(9)) -= lambda_mode(4)*(xg(8)+1);
			f(8) -= lambda_mode(4)*q_mode(9);
			f(mode_index(5)) += lambda_mode(4)*(xg(12)+1);
			f(12) += lambda_mode(4)*q_mode(5);
			//5
			f(mode_index(9)) -= lambda_mode(5)*xg(9);
			f(9) -= lambda_mode(5)*q_mode(9);
			f(mode_index(6)) += lambda_mode(5)*(xg(12)+1);
			f(12) += lambda_mode(5)*q_mode(6);
			//6
			f(mode_index(9)) -= lambda_mode(6)*xg(10);
			f(10) -= lambda_mode(6)*q_mode(9);
			f(mode_index(7)) += lambda_mode(6)*(xg(12)+1);
			f(12) += lambda_mode(6)*q_mode(7);
			//7
			f(mode_index(9)) -= lambda_mode(7)*xg(11);
			f(11) -= lambda_mode(7)*q_mode(9);
			f(mode_index(8)) += lambda_mode(7)*(xg(12)+1);
			f(12) += lambda_mode(7)*q_mode(8);
			//8
			f(mode_index(9)) -= lambda_mode(8)*(xg(4)+1);
			f(4) -= lambda_mode(8)*q_mode(9);
			f(mode_index(1)) += lambda_mode(8)*(xg(12)+1);
			f(12) += lambda_mode(8)*q_mode(1);
			//}

		}
	}
	TMStopTimer(22);
}

template <class RIGID>
void GCMSElement<RIGID>::EvalG(Vector& g, double t)
{
	if (!UseRigidBodyConstraints() && !UseFlexibleModesConstraints())
	{ return; }

	int offset_rigid = 0;
	ConstVector<12> xg(SOSRigid()), xgp(SOSRigid());
	for (int i=1; i <= SOSRigid(); i++) 
	{
		xg(i) = XG(i);
		xgp(i) = XGP(i);
	}

	if (UseRigidBodyConstraints())
	{
		offset_rigid = 6;
		if (MaxIndex()==3)
		{ 
			// 3 conditions:
			// Norm of i-th line of rotation matrix = 1, ||A_i|| = 1
			g(1) = sqr(xg(4)+1)+sqr(xg(5))+sqr(xg(6));
			g(2) = sqr(xg(7))+sqr(xg(8)+1)+sqr(xg(9));
			g(3) = sqr(xg(10))+sqr(xg(11))+sqr(xg(12)+1);
			// 3 conditions:
			// ith and jth line of rotation matrix are orthogonal, A_i . A_j = 0
			g(4) = (xg(4)+1)*xg(7)+xg(5)*(xg(8)+1)+xg(6)*xg(9);
			g(5) = (xg(4)+1)*xg(10)+xg(5)*xg(11)+xg(6)*(xg(12)+1);
			g(6) = (xg(7))*xg(10)+(xg(8)+1)*xg(11)+xg(9)*(xg(12)+1);
		}
		else // Index <= 2
		{
			// 3 conditions:
			// Norm of i-th line of rotation matrix = 1, ||A_i|| = 1
			g(1) = 2*((xg(4)+1)*xgp(4)+xg(5)*xgp(5)+xg(6)*xgp(6));
			g(2) = 2*(xg(7)*xgp(7)+(xg(8)+1)*xgp(8)+xg(9)*xgp(9));
			g(3) = 2*(xg(10)*xgp(10)+xg(11)*xgp(11)+(xg(12)+1)*xgp(12));
			// 3 conditions:
			// ith and jth line of rotation matrix are orthogonal, A_i . A_j = 0
			g(4) = (xg(4)+1)*xgp(7)+(xgp(4))*xg(7)+xg(5)*xgp(8)+xgp(5)*(xg(8)+1)+xg(6)*xgp(9)+xgp(6)*xg(9);
			g(5) = (xg(4)+1)*xgp(10)+(xgp(4))*xg(10)+xg(5)*xgp(11)+xgp(5)*xg(11)+xg(6)*xgp(12)+xgp(6)*(xg(12)+1);
			g(6) = (xg(7))*xgp(10)+(xgp(7))*xg(10)+(xg(8)+1)*xgp(11)+xgp(8)*xg(11)+xg(9)*xgp(12)+xgp(9)*(xg(12)+1);
		}
	}

	// constraints for perpendicular mode shapes generated from one local mode
	ConstVector<9> q(GetNDofPerLocalMode()), qp(GetNDofPerLocalMode());
	if (UseFlexibleModesConstraints())
	{
		for (int mode=1; mode<=NModes(); mode++)
		{
			// degrees of freedom belonging to local mode number mode
			for (int i=1; i<=GetNDofPerLocalMode(); i++)
			{
				q(i) = XG(SOSRigid()+(mode-1)*GetNDofPerLocalMode()+i);
				qp(i) = XGP(SOSRigid()+(mode-1)*GetNDofPerLocalMode()+i);
			}
			// local degree of freedom from the FFRF
			// e.g. first constraint
			// q(1) = A_11 * qf_prime = (xg(4)+1)*qf_prime   AND  q(2) = A_21*qf_prime = xg(5)*qf_prime ---->
			// q(1)*xg(5) = q(2)*(xg(4)+1);
			// constraints -- index 3 problem
			if (MaxIndex()==3)
			{
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+1) = q(2)*(xg(12)+1.) - xg(5)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+2) = q(3)*(xg(12)+1.) - xg(6)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+3) = q(4)*(xg(12)+1.) - xg(7)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+4) = q(5)*(xg(12)+1.) - (xg(8)+1.)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+5) = q(6)*(xg(12)+1.) - xg(9)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+6) = q(7)*(xg(12)+1.) - xg(10)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+7) = q(8)*(xg(12)+1.) - xg(11)*q(9);
				g(offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1)+8) = q(1)*(xg(12)+1.) - (xg(4)+1.)*q(9);
			}
			else // constraints -- reduced index by time derivative
			{
				//if ( fabs(xg(4)+1.)>0.2 )
				//{
				int off = offset_rigid+(mode-1)*(GetNDofPerLocalMode()-1);
				g(off+1) = qp(2)*(xg(12)+1.)+q(2)*xgp(12) - xgp(5)*q(9)-xg(5)*qp(9);
				g(off+2) = qp(3)*(xg(12)+1.)+q(3)*xgp(12) - xgp(6)*q(9)-xg(6)*qp(9);
				g(off+3) = qp(4)*(xg(12)+1.)+q(4)*xgp(12) - xgp(7)*q(9)-xg(7)*qp(9);
				g(off+4) = qp(5)*(xg(12)+1.)+q(5)*xgp(12) - xgp(8)*q(9)-(xg(8)+1.)*qp(9);
				g(off+5) = qp(6)*(xg(12)+1.)+q(6)*xgp(12) - xgp(9)*q(9)-xg(9)*qp(9);
				g(off+6) = qp(7)*(xg(12)+1.)+q(7)*xgp(12) - xgp(10)*q(9)-xg(10)*qp(9);
				g(off+7) = qp(8)*(xg(12)+1.)+q(8)*xgp(12) - xgp(11)*q(9)-xg(11)*qp(9);
				g(off+8) = qp(1)*(xg(12)+1.)+q(1)*xgp(12) - xgp(4)*q(9)-(xg(4)+1.)*qp(9);
				//}

			}
		}
	}
}
	

// can be set to 2 if the stiffness matrix routine works fine 
template <class RIGID>
int GCMSElement<RIGID>::FastStiffnessMatrix() const {if (!UseRigidBodyConstraints()) {return 2;} else {return 0;} }
//int GCMSElement<RIGID>::FastStiffnessMatrix() const {return 0; }
// returns A Kr A^T //fill in sos x sos components only of stiffness matrix, m might be larger
// tests with other parts of stiffness matrix including differentiation of follower loads, dA/dq, ...
template <class RIGID>
void GCMSElement<RIGID>::StiffnessMatrix(Matrix& m)
{
	m.SetAll(0.);

	xg.SetLen(SOS());
	for (int i=1; i <= SOS(); i++) xg(i) = XG(i);

	Matrix3D A, AT;
	A = GetRotMatrix(xg);
	AT = A.GetTp();

	// Assemble reduced system matrix if not available
	if (Kr.Getrows() != SOS() )
	{
		Kr.SetSize(SOS(), SOS());
		Kr.SetAll(0);
		ConstMatrix<FEmaxDOF*CMSmaxDOF> temp_PhiK, temp_Phi;
		ConstMatrix<CMSmaxDOF*CMSmaxDOF> addK;
		ConstMatrix<FEmaxDOF*FEmaxDOF> temp_K;

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
			e.StiffnessMatrix(temp_K);

			// Kr = Phi_CB^T K_el Phi_CB

			// temp_Phi contains ltg-map of Phi_CB
			temp_Phi.SetSize(sos_el, SOS());
			for (int j=1; j<=sos_el; j++)
			{
				for (int k=1; k<=SOS(); k++)
				{
					temp_Phi(j,k) = Phi_CB(ltg(j),k);
				}
			}
			Mult(temp_K, temp_Phi, temp_PhiK);
			MultTp(temp_Phi,temp_PhiK, addK);
			// then Kr is positive semidefinite (element.StiffnessMatrix is negative semidef.)
			Kr -= addK;


		}
	}

	// Stiffness matrix = blockdiag(A) K blockdiag(A^T)
	m.AddSubmatrix(Kr,4,4,4,4,SOS()-3,SOS()-3,-1.);
	ApplyRotationFromRight(AT,m);
	ApplyRotationFromLeft(A,m);



	//Matrix AKrAT(Kr);
	//ApplyRotationFromRight(AT,AKrAT);
	//ApplyRotationFromLeft(A,AKrAT);
	////m.AddSubmatrix(AKrAT,1,1,1,1,SOS(),SOS(),-1.);



	////// dqr/dqi
	//ConstVector<CMSmaxDOF> DAijDq[3][3];
	//ComputeDADq(DAijDq);
	//for (int i=4; i<=12; i++)
	//{
	//	Matrix3D dAdq(0);
	//	for (int k=0; k<3; k++)
	//	{
	//		for (int j=0; j<3; j++)
	//		{
	//			dAdq(k+1,j+1) = DAijDq[k][j](i);
	//		}
	//		dAdq(k+1,k+1)+=1.;
	//	}
	//	ConstVector<CMSmaxDOF> AKrATqr;
	//	ConstVector<CMSmaxDOF> urigid;

	//	// dqr/dqi
	//	ComputeURigid(dAdq,Vector3D(0.),urigid);
	//	// AKrATqr = AKrAT * urigid;
	//	Mult(urigid, AKrAT, AKrATqr);
	//	for (int j=1; j<=SOS(); j++)
	//	{
	//		m(j,i) += AKrATqr(j);
	//	}

	//}

	//// ************  follower loads + constraints...
	//if (loads.Length() > 0 || constraints.Length() > 0)
	//{
	//static Vector load(SOS()), load2(SOS());
	//load.SetLen(SOS());
	//load2.SetLen(SOS());
	//load.SetAll(0);
	//load2.SetAll(0);

	//for (int il=1; il <= loads.Length(); il++)
	//{
	//	loads(il)->AddElementLoad(load,GetMBS()->GetTime());
	//}
	//for (int il=1; il <= constraints.Length(); il++)
	//{
	//	constraints(il)->AddElementCqTLambda(constraintindices(il),load);
	//}

	//double deps2=1e-10;
	//for (int i=1; i <= SOS(); i++)
	//{
	//	XG(i) += deps2;
	//	load2.SetAll(0);
	//	for (int il=1; il <= loads.Length(); il++)
	//	{
	//		loads(il)->AddElementLoad(load2,GetMBS()->GetTime());
	//	}
	//	for (int il=1; il <= constraints.Length(); il++)
	//	{
	//		constraints(il)->AddElementCqTLambda(constraintindices(il),load2);
	//	}
	//	load2 -= load;
	//	load2 *= 1./deps2;

	//	for (int j=1; j<=SOS(); j++)
	//		m(j,i) += load2(j);
	//	XG(i) -= deps2;
	//}
	//}

	//// ***********  nonlinear terms, should be very small!!
	//ConstVector<CMSmaxDOF> temp; 
	//ConstVector<CMSmaxDOF> qF;
	//ConstVector<CMSmaxDOF> ATqF;
	//ConstVector<CMSmaxDOF> urigid;
	//ConstVector<CMSmaxDOF> KATqF, KAijTF;
	//ConstVector<CMSmaxDOF> dAijdq;

	//// rigid body motion
	//Vector3D uref1 = GetTranslation(xg);


	////+++++++++++++++++++++++++++
	//// e_i^T diag(dA/dq_j) K diag(A^T) (q-qR)

	//// compute qR, qF = q-qR
	//ComputeURigid(A,uref1,urigid);
	//qF = xg;
	//qF -= urigid;

	//// compute diag(A^T) (q-qR)
	//ATqF = qF;
	//ApplyRotation(AT,ATqF);

	//// compute K diag(A^T) (q-qR)
	//Mult(Kr,ATqF,KATqF);

	//Matrix3D Aij, AijT;
	//double DfnlDAij;

	//// compute DAij/dqk
	//for (int i = 1; i <= Dim(); i++)
	//{
	//	for (int j = 1; j <= Dim(); j++)
	//	{
	//		Aij.SetAll(0);
	//		Aij(i,j) = 1;
	//		AijT.SetAll(0.);
	//		AijT(j,i) = 1.;

	//		// compute temp = diag(dA/dA_ij) K diag(AT) (q-qR)
	//		temp = KATqF;
	//		ApplyRotation(Aij,temp);
	//		// subtract  dA_ij/dq_k [diag(dA/dA_ij) K diag(AT) (q-qR))] e_l
	//		for (int k=1; k<=SOS(); k++)
	//		{
	//			for (int l=1; l<=SOS(); l++)
	//			{
	//				m(l,k) -= DAijDq[i-1][j-1](k)*temp(l);;
	//			}
	//		}
	//		// subtract  dA_ij/dq_k [diag(A) K diag(dA/dA_ij^T) (q-qR))] e_l
	//		// compute diag(A/dA_ij^T) (q-qR)
	//		ATqF = qF;
	//		ApplyRotation(AijT,ATqF);

	//		// compute K diag(A/dA_ij^T) (q-qR)
	//		Mult(Kr,ATqF,temp);
	//		// compute temp = diag(A) K diag(dAT/dA_ij) (q-qR)
	//		ApplyRotation(A,temp);
	//		for (int k=1; k<=SOS(); k++)
	//		{
	//			for (int l=1; l<=SOS(); l++)
	//			{
	//				m(l,k) -= DAijDq[j-1][i-1](k)*temp(l);;
	//			}
	//		}

	//	}
	//}


}


// apply rotation rot to all degrees of freedom (stored in vector u) --> diagonal matrix A
template <class RIGID>
void GCMSElement<RIGID>::ApplyRotation(const Matrix3D& rot, Vector& u) const
{
	if (SOS()%3) 
	{
		GetMBS()->UO() << "Error in GCMSElement::ApplyRotation: SOS = " << SOS() << " is not multiple of 3!\n";
		return;
	}
	int j=1;
	while (j<=SOS())
	{
		Vector3D uloc2(u(j), u(j+1), u(j+2));
		j += 3;
		Vector3D uloc;
		uloc = rot*uloc2;
		u(j-3) = uloc(1);
		u(j-2) = uloc(2);
		u(j-1) = uloc(3);
	}
}

// apply block-diagonal rotation rot to all degrees of freedom from left(stored in matrix K) --> block diagonal matrix A
// K = rot * K
template <class RIGID>
void GCMSElement<RIGID>::ApplyRotationFromLeft(const Matrix3D& rot, Matrix& K) const
{
	if (SOS()%3) 
	{
		GetMBS()->UO() << "Error in GCMSElement::ApplyRotation: SOS = " << SOS() << " is not multiple of 3!\n";
		return;
	}
	int nblocks = SOS()/3;
	for (int j=1; j<=SOS(); j++)
	{
		int k=1;
		while (k<=SOS())
		{
			Vector3D uloc(K(k,j), K(k+1,j), K(k+2,j));
			k += 3;
			uloc = rot*uloc;
			K(k-3,j) = uloc(1);
			K(k-2,j) = uloc(2);
			K(k-1,j) = uloc(3);
		}
	}
}
// apply block-diagonal rotation transposed from right to all degrees of freedom (stored in matrix K) --> block diagonal matrix A
// K = K * rot
template <class RIGID>
void GCMSElement<RIGID>::ApplyRotationFromRight(const Matrix3D& rot, Matrix& K) const
{
	if (SOS()%3) 
	{
		GetMBS()->UO() << "Error in GCMSElement::ApplyRotation: SOS = " << SOS() << " is not multiple of 3!\n";
		return;
	}
	Matrix3D rotT = rot.GetTp();
	int nblocks = SOS()/3;
	for (int j=1; j<=SOS(); j++)
	{
		int k=1;
		while (k<=SOS())
		{
			Vector3D uloc(K(j,k), K(j,k+1), K(j,k+2));
			k += 3;
			uloc = rotT*uloc;
			K(j,k-3) = uloc(1);
			K(j,k-2) = uloc(2);
			K(j,k-1) = uloc(3);
		}
	}
}

// compute degrees of freedom representing rigid body motion specified by A, uref1
template <class RIGID>
void GCMSElement<RIGID>::ComputeURigid(const Matrix3D& A, const Vector3D& uref1, Vector& urigid) const
{
	urigid.SetLen(SOS());
	urigid.SetAll(0.);
	// translation dofs: 1, 2, 3
	for (int i=1; i<=3; i++)
	{
		urigid(i) = uref1(i);
	}
	// rotation dofs
	urigid(4) = A(1,1)-1.;
	urigid(5) = A(2,1);
	urigid(6) = A(3,1);
	urigid(7) = A(1,2);
	urigid(8) = A(2,2)-1.;
	urigid(9) = A(3,2);
	urigid(10) = A(1,3);
	urigid(11) = A(2,3);
	urigid(12) = A(3,3)-1.;

}

// on Input: vec, dvecdq = d/dq(vec)
// on Output:
//		vec = vec / ||vec||
//		dvecdq = d/dq (vec/||vec||)
void DNormedVecdq(Vector3D& vec, Matrix& dvecdq)
{
	double normvec = vec.Norm();
	//dvecdq *= 1./normvec;
	for (int k=1; k<=dvecdq.Getrows(); k++)
	{
		double fac = ((vec(1)*dvecdq(k,1)+vec(2)*dvecdq(k,2)+vec(3)*dvecdq(k,3))/Cub(normvec));
		for (int i=1; i<=3; i++)
		{
			dvecdq(k,i) /= normvec;
			dvecdq(k,i) -= fac* vec(i);
		}
	}
	vec.Normalize();
}

// compute derivative of rotation matrix entries with respect to all reduced degrees of freedom
template <class RIGID>
void GCMSElement<RIGID>::ComputeDADq(ConstVector<CMSmaxDOF> dAijdq[3][3]) const
{
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			dAijdq[i][j].SetLen(SOS());
			dAijdq[i][j].SetAll(0.);		
		}
	}
	// in case of rigid body constraints, entries of A are degrees of freedom xg(4) .. xg(12)   (+1 on diagonal)
	if (UseRigidBodyConstraints())
	{
		for (int k=4; k<=12; k++)
		{
			// compute indices i,j of position A_ij of xg(k)
			int i = k%3;
			if (i==0) i=3;
			int j = (k-1)/3;
			dAijdq[i-1][j-1](k) = 1.;
		}
		return;
	}

	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgc(i) = XG(i);
	}
		
	Vector3D r1, r2, r3;
	Vector3D r01, r02, r03;
	r1 = GetNodePosRBModes(RefNode1(), xgc);
	r2 = GetNodePosRBModes(RefNode2(), xgc);
	r3 = GetNodePosRBModes(RefNode3(), xgc);
	r01 = GetNode(RefNode1()).Pos();
	r02 = GetNode(RefNode2()).Pos();
	r03 = GetNode(RefNode3()).Pos();

	Vector3D a, b, c, a0, b0, c0;
	// +++++ a0, b0, c0 for initial frame
	// +++++ a, b, c for current frame
	// +++++ dadq, dbdq, dcdq for deriveative of current frame with respect to dofs
	a0 = r02-r01;
	a0.Normalize();
	b0 = (r03-r01);
	a0.GramSchmidt(b0);
	c0 = a0.Cross(b0);
	Matrix3D B0T(a0(1), a0(2), a0(3),
		b0(1), b0(2), b0(3),
		c0(1), c0(2), c0(3));

	// dpos(refnode)/dq
	ConstMatrix<CMSmaxDOF*3> dr1dq, dr2dq, dr3dq;
	GetNodedPosdqTRBModes(RefNode1(), dr1dq);
	GetNodedPosdqTRBModes(RefNode2(), dr2dq);
	GetNodedPosdqTRBModes(RefNode3(), dr3dq);

	ConstMatrix<CMSmaxDOF*3> dadq(SOS(),3), dbdq(SOS(),3), dcdq(SOS(),3);
	a = r2-r1;
	dadq = dr2dq;
	dadq -= dr1dq;
	// compute a = a/||a||
	// compute dadq = d/dq(a/||a||)
	DNormedVecdq(a, dadq);

	b = (r3-r1);
	// diff GramSchmidt function
	// dbdq = (dr3dq-dr1dq) - ((dr3dq-dr1dq).a)a - ((r3-r1).dadq)a - ((r3-r1).a)dadq
	dbdq = dr3dq;
	dbdq -= dr1dq;
	for (int k=1; k<=dbdq.Getrows(); k++)
	{
		double faca=0;
		double facdadq=0;
		for (int i=1; i<=3; i++)
		{
			faca += a(i)*dbdq(k,i) + (r3(i)-r1(i))*dadq(k,i);
			facdadq += ((r3(i)-r1(i))*a(i));
		}
		for (int i=1; i<=3; i++)
			dbdq(k,i) -= faca*a(i) + facdadq*dadq(k,i);
	}
	b -= (a*b)*a;
	// b = b/||b||
	// dbdq = d/dq(b/||b||)
	DNormedVecdq(b, dbdq);

	c = a.Cross(b);
	//dcdq = a x dbdq + dadq x b
	for (int k=1; k<=dadq.Getrows(); k++)
	{
		dcdq(k,1) = a(2)*dbdq(k,3)-a(3)*dbdq(k,2) + dadq(k,2)*b(3)-dadq(k,3)*b(2);
		dcdq(k,2) = a(3)*dbdq(k,1)-a(1)*dbdq(k,3) + dadq(k,3)*b(1)-dadq(k,1)*b(3);
		dcdq(k,3) = a(1)*dbdq(k,2)-a(2)*dbdq(k,1) + dadq(k,1)*b(2)-dadq(k,2)*b(1);

		Matrix3D B(dadq(k,1), dbdq(k,1), dcdq(k,1),
			dadq(k,2), dbdq(k,2), dcdq(k,2),
			dadq(k,3), dbdq(k,3), dcdq(k,3));
		Matrix3D dAdq = B*B0T;
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				dAijdq[i][j](k) = dAdq(i+1,j+1);
			}
		}
	}

}
// compute derivative of rotation matrix entries with respect to reduced degree of freedom q_k
// using numerical differentiation
template <class RIGID>
void GCMSElement<RIGID>::ComputeDADqk_NumericDiff(int k, Matrix3D& dAdq)
{
	// in case of rigid body constraints, entries of A are degrees of freedom xg(4) .. xg(12)   (+1 on diagonal)
	if (UseRigidBodyConstraints())
	{
		dAdq.SetAll(0.);
		if (k>=4 && k<=12)
		{
			// compute indices i,j of position A_ij of xg(k)
			int i = k%3;
			if (i==0) i=3;
			int j = (k-1)/3;
			dAdq(i,j) = 1.;
		}
		return;
	}

	// otherwise - numerical differentiation
	// A depends only on rotation parameters dof 4 - 12
		if (k<4 || k>12)
		{
			dAdq.SetAll(0.);
			return;
		}
	double eps = 1e-10;
	Matrix3D dAdqtest;
	XG(k) -= eps;
	Matrix3D A0 = GetRotMatrix();
	XG(k) += 2*eps;
	Matrix3D Aeps = GetRotMatrix();
	XG(k) -= eps;
	dAdq = Aeps - A0;
	dAdq *= 1./(2.*eps);
	return;

}

// compute rotation matrix for actual set of reduced degrees of freedom
template <class RIGID>
Matrix3D GCMSElement<RIGID>::GetRotMatrix() const
{
	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgc(i) = XG(i);
	}
	return GetRotMatrix(xgc);
}
// compute rotation matrix for given set of reduced degrees of freedom xgc
template <class RIGID>
Matrix3D GCMSElement<RIGID>::GetRotMatrix(const Vector& xgc) const
{
	if (use_rb_constraints)
	{
		return Matrix3D(xgc(4)+1,xgc(7),xgc(10),
			              xgc(5),xgc(8)+1,xgc(11),
										xgc(6),xgc(9),xgc(12)+1);
	}
	Vector3D r1, r2, r3;
	Vector3D r01, r02, r03;
	r1 = GetNodePosRBModes(RefNode1(), xgc);
	r2 = GetNodePosRBModes(RefNode2(), xgc);
	r3 = GetNodePosRBModes(RefNode3(), xgc);
	r01 = GetNode(RefNode1()).Pos();
	r02 = GetNode(RefNode2()).Pos();
	r03 = GetNode(RefNode3()).Pos();
	Vector3D a, b, c, a0, b0, c0;
	a = r2-r1;
	a.Normalize();
	b = (r3-r1);
	// GramSchmidt funktion  -->  b -= (a*b) * a;	b.Normalize();
	a.GramSchmidt(b);
	c = a.Cross(b);
	a0 = r02-r01;
	a0.Normalize();
	b0 = (r03-r01);
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


// compute rotation matrix for actual set of reduced degrees of freedom
template <class RIGID>
Matrix3D GCMSElement<RIGID>::GetRotMatrixD() const
{
	ConstVector<CMSmaxDOF> xgd(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgd(i) = XGD(i);
	}
	return GetRotMatrix(xgd);
}

// compute rotation matrix for given set of reduced degrees of freedom xgc
template <class RIGID>
Matrix3D GCMSElement<RIGID>::GetRotMatrixP() const
{
	ConstVector<CMSmaxDOF> xgP(SOS());
	for (int i=1; i<=SOS(); i++)
	{
		xgP(i) = XGP(i);
	}

	Matrix3D AP(0.);
	ConstVector<CMSmaxDOF> DAijDq[3][3];
	ComputeDADq(DAijDq);

	// AP = sum_k DAij/dqk qkP
	GCMSElement<RIGID>& constthis = const_cast<GCMSElement<RIGID>&> (*this);
	for (int k=1; k<=SOSRigid(); k++)
	{
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				AP(i+1,j+1) += DAijDq[i][j](k)*xgP(k);
			}
		}
	}

	return AP;
}

// omega_tilde = AP * A^T
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetAngularVel(const Vector3D& p_loc) const 
{
	Matrix3D A = GetRotMatrix();
	Matrix3D AT = A.GetTp();
	Matrix3D AP = GetRotMatrixP();
	Matrix3D omegatilde = AP*AT;
	//GetMBS()->UO() << "omega tilde = " << omegatilde << "\n";

	Vector3D omega(omegatilde(3,2),omegatilde(1,3),omegatilde(2,1));
	return omega;
}

// omega_bar_tilde = A^T * AP
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetAngularVelLocal(const Vector3D& p_loc) const 
{
	Matrix3D A = GetRotMatrix();
	Matrix3D AT = A.GetTp();
	Matrix3D AP = GetRotMatrixP();
	Matrix3D omegatilde = AT*AP;

	Vector3D omegabar(-omegatilde(2,3),omegatilde(1,3),-omegatilde(1,2));
	return omegabar;
}

// compute translation vector for given set of reduced degrees of freedom xgc
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetTranslation(const Vector& xgc) const
{
	Vector3D r1;
	Vector3D r01;
	r1 = GetNodePos(RefNode1(), xgc);
	r01 = GetNode(RefNode1()).Pos();
	return (r1-r01);
}

// compute reference position = position of RefNode1
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetRefPosD() const
{
	return GetNodePosD(RefNode1());
}
// compute initial reference position = initial position of RefNode1
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetRefPosInit() const
{
	return GetNode(RefNode1()).Pos();
}
// compute reference position = position of RefNode1
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetRefPos() const
{
	return GetNodePos(RefNode1());
}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position / velocity of node in the CMS element
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodePos(int i, const Vector& xgc) const
{
	int sos = SOS();
	const Node& nodei = GetNode(i);
	Vector3D pglob = nodei.Pos();
	// AP: old, save version
	//for (int j=1; j<=sos; j++)
	//{
	//	for (int k=1; k<=Dim(); k++)
	//	{
	//		int nodedofk = nodei.Get(k);
	//		pglob(k) += Phi_CB(nodedofk, j) * (xgc(j));
	//	}
	//}
	// AP: end old, save version
	// AP: new version, uses that PHI_CB consists of multiples of 3-Id-matrix
	int sos3 = sos/3;
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=sos3; j++)
		{
			pglob(k) += Phi_CB(nodedofk, (j-1)*3+k) * (xgc((j-1)*3+k));
		}
	}
	return pglob;
}
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// position of node in the CMS element, only taking rigid body modes into account
template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodePosRBModes(int i, const Vector& xgc) const
{
	int sos_rigid = SOSRigid();
	const Node& nodei = GetNode(i);
	Vector3D pglob = nodei.Pos();
	for (int k=1; k<=Dim(); k++)
	{
		for (int j=1; j<=sos_rigid; j++)
		{
			int nodedofk = nodei.Get(k);
			pglob(k) += Phi_CB(nodedofk, j) * (xgc(j));
		}
	}
	return pglob;
}

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodePos(int i) const
{
	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int j=1; j<=SOS(); j++)
	{
		xgc(j) = XG(j);
	}
	return GetNodePos(i, xgc);
}


template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodeLocPos(int i) const
{
	GetMBS()->UO() << "local node position not available in GCMSElement!\n";
	return Vector3D(0.);
}

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodeVel(int i) const
{
	ConstVector<CMSmaxDOF> xgp(SOS());
	for (int j=1; j<=SOS(); j++)
	{
		xgp(j) = XGP(j);
	}
	return GetNodeVel(i, xgp);
}

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodeVel(int i, const Vector& xgp) const
{
	int sos = SOS();
	const Node& nodei = GetNode(i);
	Vector3D vel(0.);
	// AP: old, save version
	//for (int k=1; k<=Dim(); k++)
	//{
	//	for (int j=1; j<=sos; j++)
	//	{
	//		int nodedofk = nodei.Get(k);
	//		vel(k) += Phi_CB(nodedofk, j) * xgp(j);
	//	}
	//}
	// AP: old, save version
	// AP: new, fast version using that Phi_CB is multiple of 3-Id-blocks
	int sos3 = sos/3;
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=sos3; j++)
		{
			vel(k) += Phi_CB(nodedofk, 3*(j-1)+k) * xgp(3*(j-1)+k);
		}
	}
	return vel;
}

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetNodePosD(int i) const
{
	ConstVector<CMSmaxDOF> xgc(SOS());
	for (int j=1; j<=SOS(); j++)
	{
		xgc(j) = XGD(j);
	}
	return GetNodePos(i, xgc);
}


// d = d pos / dq(ploc)
// for rigid body dofs only, or not at all??
template <class RIGID>
void GCMSElement<RIGID>::GetdPosdqT(const Vector3D& ploc, Matrix& d)
{
	GetMBS()->UO() << "DPosDq not available in GCMSElement!\n";
	int sos = SOS();
	d.SetSize(sos,3);
	d.SetAll(0);
}

// d = d vel / dq(ploc,vloc)
// for rigid body dofs only, or not at all??
template <class RIGID>
void GCMSElement<RIGID>::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	GetMBS()->UO() << "DRotvDq not available in GCMSElement!\n";
	;
}

template <class RIGID>
void GCMSElement<RIGID>::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	GetMBS()->UO() << "DRotDq not available in GCMSElement!\n";
}

// differentiation of global position of Node node with respect to REDUCED degrees of freedom!
template <class RIGID>
void GCMSElement<RIGID>::GetNodedPosdqT(int node, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();
	int sos = SOS();
	const Node& nodei = GetNode(node);
	// AP: old, save version
	//for (int k=1; k<=Dim(); k++)
	//{
	//	int nodedofk = nodei.Get(k);
	//	for (int j=1; j<=sos; j++)
	//	{
	//		dpdqi(j,k) += Phi_CB(nodedofk, j);
	//	}
	//}
	// AP: end old, save version
	// AP: new fast version, using that Phi_CB contains multiples of 3-Id-matrix
	int sos3 = sos /3;
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=sos3; j++)
		{
			dpdqi((j-1)*3+k,k) += Phi_CB(nodedofk, (j-1)*3+k);
		}
	}
	// AP: end new fast version

	return;
}
	// differentiation of position of Node node due to the 12 RB-Modes
template <class RIGID>
void GCMSElement<RIGID>::GetNodedPosdqTRBModes(int node, Matrix& dpdqi) const
{
	dpdqi.SetSize(SOS(),Dim());
	dpdqi.FillWithZeros();
	int sosrigid = SOSRigid();
	const Node& nodei = GetNode(node);
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=sosrigid; j++)
		{
			dpdqi(j,k) += Phi_CB(nodedofk, j);
		}
	}
	return;
}

template <class RIGID>
void GCMSElement<RIGID>::AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f)   // f += dpdq*lambda
{
	int sos = SOS();
	const Node& nodei = GetNode(node);
	for (int k=1; k<=Dim(); k++)
	{
		int nodedofk = nodei.Get(k);
		for (int j=1; j<=sos; j++)
		{
			f(k) += Phi_CB(nodedofk, j)*lambda(k);
		}
	}
	return;
}

// int du/dq dV for the GCMS element, linear
template <class RIGID>
void GCMSElement<RIGID>::GetIntDuDq(Matrix& dudq) //in fact it is DuDq Transposed
{
	if (Hr.Getrows() == SOS())
	{
		dudq = Hr;
	}
	else
	{
		Hr.SetSize(SOS(), Dim());
		Hr.FillWithZeros();

		Matrix H_loc;
		for (int el=1; el<=NFFRFElements(); el++) 
		{
			Element& e = GetFFRFElement(el);
			int sos_loc = e.FlexDOF();
			const	TArray<int>& ltg = e.GetLTGArray();

			e.GetH(H_loc);
			// add Phi_CB^T * H_loc to global H matrix
			for (int j=1; j <= sos_loc; j++)
			{
				for (int k=1; k <= Dim(); k++)
				{
					for (int i=1; i<=SOS(); i++)
					{
						Hr(i,k) += Phi_CB(ltg(j),i)*H_loc(j,k);
					}
				}
			}
		}
		dudq = Hr;
	}
}

// int rho du/dq dV for the full cms element, nonlinear, depends on xg
template <class RIGID>
void GCMSElement<RIGID>::GetIntRhoDuDq(Matrix& rhodudq) //in fact it is DuDq Transposed
{
	if (intrhodudqr.Getrows() == SOS())
	{
		rhodudq = intrhodudqr;
	}
	else
	{
	rhodudq.SetSize(SOS(), Dim());
	rhodudq.FillWithZeros();

	Matrix rhodudq_loc;
	for (int el=1; el<=NFFRFElements(); el++) 
	{
		Element& e = GetFFRFElement(el);
		int sos_loc = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		e.GetIntRhoDuDq(rhodudq_loc);
		// add Phi_CB^T * H_loc to global H matrix
		for (int j=1; j <= sos_loc; j++)
		{
			for (int k=1; k <= Dim(); k++)
			{
				for (int i=1; i<=SOS(); i++)
				{
					rhodudq(i,k) += Phi_CB(ltg(j),i)*rhodudq_loc(j,k);
				}
			}
		}
	}

	// ------- Add additional masses at nodes
	for (int mn=1; mn<=massnodes.Length(); mn++)
	{
		Node& node = GetNode(massnodes(mn));
		for (int j=1; j<=node.SOS(); j++)
		{
			for (int i=1; i<=SOS(); i++)
			{
				rhodudq(i,j) += massvalues(i) * Phi_CB(node.Get(j),i);
			}
		}
	}

	intrhodudqr = rhodudq;
	}
}

template <class RIGID>
void GCMSElement<RIGID>::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	BaseCMSElement::GetElementData(edc);

	ElementData ed;

	ed.SetInt(refnode1, "RefNode1"); edc.Add(ed);
	ed.SetInt(refnode2, "RefNode1"); edc.Add(ed);
	ed.SetInt(refnode3, "RefNode1"); edc.Add(ed);
	ed.SetInt(flag_rotation_axis, "Flag_RotationAxis"); edc.Add(ed);
	ed.SetVector3D(rotation_axis(1), rotation_axis(2), rotation_axis(3), "RotationAxis"); edc.Add(ed);
}

template <class RIGID>
int GCMSElement<RIGID>::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = BaseCMSElement::SetElementData(edc);

	GetElemDataInt(mbs, edc, "RefNode1", refnode1, 1);
	GetElemDataInt(mbs, edc, "RefNode2", refnode2, 1);
	GetElemDataInt(mbs, edc, "RefNode3", refnode3, 1);
	GetElemDataInt(mbs, edc, "Flag_RotationAxis", flag_rotation_axis, 1);
	GetElemDataVector3D(mbs, edc, "RotationAxis", rotation_axis(1), rotation_axis(2), rotation_axis(3), 1);
	return rv;
}

template <class RIGID>
void GCMSElement<RIGID>::DrawElement()
{
	ReferenceFrame3D::DrawElement();

	////// HACK AP:	
	//// draw position of reference point 3 after full deformation and after rigid body transformation only
	//Vector3D refp(GetPosD(GetNode(RefNode3()).Pos()));
	//GetMBS()->DrawSphere(refp, 0.005);
	//Vector3D refp2(GetNodePosD(RefNode3()));
	//GetMBS()->DrawSphere(refp2, 0.005);

	//// HACK AP: draw also urigid, and reference frame for urigid (has to be the same as original)
	//ConstVector<CMSmaxDOF> urigid, xgdstore;
	//xgdstore.SetLen(SOS());
	//for (int i=1; i<=SOS(); i++)
	//{
	//	xgdstore(i) = XGD(i);
	//}
	//Matrix3D A = GetRotMatrix(xgdstore);
	//Vector3D uref1 = GetTranslation(xgdstore);
	//ComputeURigid(A,uref1,urigid);

	//for (int i=1; i<=SOS(); i++)
	//	mbs->GetDrawValue(ltg.Get(i)) = urigid(i);
	//ReferenceFrame3D::DrawElement();

	//for (int i=1; i<=NFFRFElements(); i++)
	//{
	//	Element& el = GetFFRFElement(i);
	//	el.DrawElement();
	//}

	//for (int i=1; i<=SOS(); i++)
	//	mbs->GetDrawValue(ltg.Get(i)) = xgdstore(i);

	////// HACK AP: draw modes
	//ConstVector<CMSmaxDOF> xgdstore;
	//xgdstore.SetLen(SOS());
	//for (int i=1; i<=SOS(); i++)
	//	xgdstore(i) = XGD(i);

	//for (int i=1; i<=SOS(); i++)
	//	mbs->GetDrawValue(ltg.Get(i)) = 0.0;
	//mbs->GetDrawValue(ltg.Get(14))=1.;
	//mbs->GetDrawValue(ltg.Get(16))=-1.;
	//mbs->GetDrawValue(ltg.Get(21))=1;
	//mbs->GetDrawValue(ltg.Get(4))=-1.;
	//mbs->GetDrawValue(ltg.Get(5))=1.;
	//mbs->GetDrawValue(ltg.Get(7))=-1;
	//mbs->GetDrawValue(ltg.Get(8))=-1;

	//for (int i=1; i<=NFFRFElements(); i++)
	//{
	//	Element& el = GetFFRFElement(i);
	//	el.DrawElement();
	//}
	//for (int i=1; i<=SOS(); i++)
	//	mbs->GetDrawValue(ltg.Get(i)) = 0.0;
	//mbs->GetDrawValue(ltg.Get(13))=1.;
	//mbs->GetDrawValue(ltg.Get(17))=1.;
	//mbs->GetDrawValue(ltg.Get(21))=1.;
	//for (int i=1; i<=NFFRFElements(); i++)
	//{
	//	Element& el = GetFFRFElement(i);
	//	el.DrawElement();
	//}

	//for (int i=1; i<=SOS(); i++)
	//	mbs->GetDrawValue(ltg.Get(i)) = xgdstore(i);
}


template <class RIGID>
Vector3D GCMSElement<RIGID>::GetPos(const Vector3D& p_loc) const
{
	return GetRotMatrix()*(p_loc-GetNode(RefNode1()).Pos())+GetRefPos();
};

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetVel(const Vector3D& p_loc) const
{
	GetMBS()->UO() << "Error in GCMSElement::GetVel: rigid body velocity not implemented for GCMSElement\n";
	return Vector3D(0.);
};

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetPosD(const Vector3D& p_loc) const
{
	return GetRotMatrixD()*(p_loc-GetNode(RefNode1()).Pos())+GetRefPosD();
};

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetVelD(const Vector3D& p_loc) const
{
	GetMBS()->UO() << "Error in GCMSElement::GetVelD: rigid body velocity not implemented for GCMSElement\n";
	return Vector3D(0.);
};

// Load Eigenmodes into Data Manager
// ATTENTION: these are unit vectors locally on the element
template <class RIGID>
void GCMSElement<RIGID>::EigenModesToDataManager()
{
	// number of Modes which is loaded into data manager = FlexDOF()
	// NModes is not correct for GCMSElement, correct for CMSElement
	ElementDataContainer edc;
	for (int j=1; j<=SOS(); j++)
	{
		// solution vector - is essentially a unit vector with 1 in the LTG(j)-th place
		Vector solvec;
		solvec.SetLen(GetMBS()->GetSystemSize());
		// only mode number j is used
		solvec.SetAll(0.);
		solvec(LTG(j)) = 1;

		ElementData ed; 
		ed.SetVector(solvec.GetVecPtr(),solvec.GetLen(),mystr("SV")+mystr(j));
		edc.Add(ed);

		// no data vector
		double *datavec = 0;
		ed.SetVector(datavec,0,mystr("DV")+mystr(j));
		edc.Add(ed);
	}
	mbs->UO().CallWCDriverFunction(20,SOS(),0,&edc);
}

template <class RIGID>
Vector3D GCMSElement<RIGID>::GetAngularMomentum(const Vector3D& p_ref) const
{
	Vector3D D(0.,0.,0.);
	for (int i = 1; i <= NFFRFElements(); i++)
	{
		D += GetFFRFElement(i).GetAngularMomentum(p_ref);
	}
	return D;
}

template <class RIGID>
double GCMSElement<RIGID>::GetKineticEnergy()
{
	Matrix m(SOS(),SOS());
	Vector qp(SOS());

	EvalM(m, 0.);

	for (int i = 1; i <= SOS(); i++)
	{
		qp(i) = XGP(i);
	}

	return 0.5*qp*(m*qp);
}

template <class RIGID>
double GCMSElement<RIGID>::GetPotentialEnergy()
{
	Vector q(SOS()), f(SOS());
	EvalF2(f, 0.);

	for (int i = 1; i <= SOS(); i++)
	{
		q(i) = -XG(i);
	}

	return f*q;
}

template class GCMSElement<Rigid3D>;
template class GCMSElement<Rigid3DKardan>;
