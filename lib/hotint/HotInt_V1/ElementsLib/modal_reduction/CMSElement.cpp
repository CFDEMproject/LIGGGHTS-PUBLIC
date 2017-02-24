//#**************************************************************
//#
//# filename:             CMSElement.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						11. July 2007
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
#include "cmselement.h"
//#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"
#include "lapack_routines.h"
#include "lapack_routines.h"
#include "finiteelement3d.h"
#include "myfile.h"


template<class RIGID>
void CMSElement<RIGID>::SetCMSElement(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
																			const Vector3D& sizei, const Vector3D& coli)
{
	n_zeromodes = -1;
	SetReferenceFrame(p, v, phi, phip, sizei, coli);
	type = (TMBSElement)(type|TCMSflag);

	nbmodes = 0;
	nimodes = 0;

	// compute initial conditions for nmodes=0, i.e. for 2*7+1 dofs
	//ComputeInitialConditions(p, v, phi, phip, x_init);

	// then set nimodes, x_init for system with flexible and rigid dofs will be computed later (in DoModalAnalysis)
	nimodes = nimodesi;
	//col = coli;
	//size = sizei;

	//InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
	boundarydofelem.SetLen(0);
	boundarynode.SetLen(0);

	solverparameters.DefaultInitialize();

	EVfile = mystr("");


	//////// HACK AP
	//double xinit_double[12] = {0.9771035473772125, 0.8119080368041125, 0.8149451665816472, 0.4345721941229590, 0.4196667118043127, 0.4654793137759562, 0.6468043502831272, 0.1137455299827921, 0.8447456427769386, 0.7828729837346822, -0.05939125082671293, -0.4017941539650536};
	//for (int i=1; i<=12; i++)
	//	x_init(i) = xinit_double[i-1];

};

template<>
void CMSElement<Rigid3DKardan>::SetCMSElementKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
																										const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs)
{
	n_zeromodes = -1;
	SetReferenceFrameKardan(p, v, phi, phip, sizei, coli, rs);

	type = (TMBSElement)(type|TCMSflag);

	nbmodes = 0;
	nimodes = 0;

	//// compute initial conditions for nmodes=0, i.e. for 2*7+1 dofs
	//ComputeInitialConditions(p, v, phi, phip, x_init);

	// then set nimodes, x_init for system with flexible and rigid dofs will be computed later (in DoModalAnalysis)
	nimodes = nimodesi;
	//col = coli;
	//size = sizei;

	//InitializeSearchtree(-0.6*size,0.6*size,10,10,10);
	boundarydofelem.SetLen(0);
	boundarynode.SetLen(0);

	solverparameters.DefaultInitialize();

	EVfile = mystr("");
}

template<>
void CMSElement<Rigid3D>::SetCMSElementKardan(const Vector3D& p, const Vector3D& v, Vector3D phi, Vector3D phip, const int nimodesi,
																							const Vector3D& sizei, const Vector3D& coli, Rigid3DKardan::RotationsSequence rs)
{
	GetMBS()->UO().InstantMessageText("SetCMSElementKardan called for CMSElement<Rigid3D>!\n");
};





//this function is called after assembly and AFTER filysis
template <class RIGID>
void CMSElement<RIGID>::Initialize() 
{
	//UO().InstantMessageText("begin CMSElement<RIGID>::Initialize");
	BaseCMSElement<RIGID>::Initialize();

	if(!animation_mode)
	{
		InitializeFFRFMatrices();
	}


	int sosfull = this->SOSFull();
}

template <class RIGID>
void CMSElement<RIGID>::InitializeFFRFMatrices()
{
	//if (NModes() != 0)
	{

		int dim = Dim();

		//GetMBS()->UO() << "I1..\n";
		GetI1(I1S);

		//GetMBS()->UO() << "Sbar..\n";
		GetSbar(SbarS);

		//GetMBS()->UO() << "Ikl..\n";
		for (int k=1; k <= dim; k++)
		{
			for (int l=1; l <= dim; l++)
			{
				IklS[k-1][l-1] = GetIkl(k, l);
				GetIbarkl(k, l, IbarklS[k-1][l-1]);
			}
		}


		//GetMBS()->UO() << "Sbarkl..\n";
		double ts = -GetClockTime();

		for (int k=1; k <= dim; k++)
		{
			for (int l=1; l <= dim; l++)
			{
				GetSbarkl(k, l, SbarklS[k-1][l-1]);
				//UO() << "SbarklS[" << k << "][" << l << "=" << SbarklS[k-1][l-1] << "\n";
				// use symmetry: Sbar_kl = Sbar_lk^T
				//if (k != l)
				//{
				//	SbarklS[l-1][k-1] = SbarklS[k-1][l-1];
				//	SbarklS[l-1][k-1].TpYs();
				//}
			}
		}
		ts += GetClockTime();
		UO() << "Sbarkl-time=" << ts << "\n";

		//++++++++++++++++++++++++++++++++++++++++++++++++++
		//compute reduced H matrix:
		//GetMBS()->UO() << "Hr..\n";

		ConstMatrix<CMSmaxDOF*3> add_H;
		ConstMatrix<FEmaxDOF*CMSmaxDOF> local_Phi;
		ConstMatrix<FEmaxDOF*3> local_H;
		Hr.SetSize(NModes(), Dim());
		Hr.FillWithZeros();

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();
			const	TArray<int>& ltg = e.GetLTGArray();

			e.GetH(local_H);
			// set element-local Phi_CB matrix
			GetLocalMatrix(Phi_CB, ltg, sos, local_Phi);

			// local_Phi^T * local_H
			MultTp(local_Phi, local_H, add_H);
			Hr += add_H;
		}

		GetMBS()->UO() << "I11 = " << IklS[2][2] +IklS[1][1]<< "\n";
		GetMBS()->UO() << "I22 = " << IklS[2][2] +IklS[0][0]<< "\n";
		GetMBS()->UO() << "I33 = " << IklS[1][1] +IklS[0][0]<< "\n";

		// AH
		Iphi(1,1) = IklS[2][2] + IklS[1][1];
		Iphi(2,2) = IklS[2][2] + IklS[0][0];
		Iphi(3,3) = IklS[1][1] + IklS[0][0];

		Vector3D v_cg(I1S[0], I1S[1], I1S[2]);
		v_cg /= GetMass();
		GetMBS()->UO() << "Mass = " << GetMass() << "\n";
		GetMBS()->UO() << "Center of gravity v_cg " << v_cg << "\n";;

		//mystr v1, v2, v3;
		//v1.SmartDouble2String(v_cg.X());
		//v2.SmartDouble2String(v_cg.Y());
		//v3.SmartDouble2String(v_cg.Z());
		//GetMBS()->UO() << "Center of gravity v_cg: [" << v1 << "," << v2 << "," << v3 << "]\n";;

		GetMBS()->UO() << "..init CMS finished \n";
	}

	ConstMatrix<CMSmaxDOF*CMSmaxDOF> m(SOS(),SOS());
	EvalM(m, 0.);
};


//link elements, ltg-elements, compute nbmodes, 
//  compute full matrices M and K, modal analysis, 
//  transformation matrix, 
//  reduced modal matrices Mr and Kr, x_init!
template <class RIGID>
void CMSElement<RIGID>::DoModalAnalysis(const TArray<int2>& fixednodes)
{

	UO() << "+++++++++++++++++\nDo Modal Analysis\n+++++++++++++++++\n";
	UO() << "n-FFRF-elements=" << NFFRFElements() << "\n";

	UO() << "compute reference numbers / init \n";

	int sosfull = SOSFull();
	//UO().InstantMessageText("ComputeElementLTG");
	ComputeElementLTG();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//compute necessary boundary nodes (DOF):
	UO() << "compute boundary modes\n";
	ComputeBoundaryDofIndexList(fixednodes, sosfull);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//resort ltg according to (B) and (I) nodes:
	//UO().InstantMessageText("ResortLTG");
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

	GetMBS()->UO() << "nbmodes classic " << nb_modesclassic << "\n";


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//size of stiffness and mass-matrix: sosfull x sosfull
	int nb = NBModes();
	int nf = sosfull-nb_dofstotal; //number of internal (flexible) dof - unreduced
	int nr = NIModes();

	UO() << "sos_full=" << sosfull << "\n";
	UO() << "nb=" << nb << ", nf=" << nf << ", nr=" << nr << "\n";



	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//sparse computation of Msparse and Ksparse, used for eigenvalues in Matlab and Hotint

	// check if eigenmodes can be read from file
	// if possible: Phi_CB contains eigenmodes, read_eigenmodes_from_file = 1
	// if NOT possible: Phi_CB is empty, read_eigenmodes_from_file = 0
	//UO().InstantMessageText("ReadEigenModesFile");
	ReadEigenModesFile(Phi_CB, read_eigenmodes_from_file);
	//UO().InstantMessageText("end ReadEigenModesFile");

	if(!animation_mode)
	{
		//+++++++++++++++++++++++
		// Case: Phi_CB is read from file
		//+++++++++++++++++++++++
		if (read_eigenmodes_from_file == 1)
		{
			//UO().InstantMessageText("read_eigenmodes_from_file == 1");
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// build K_r and M_r
			Kr.SetSize(NModes(), NModes());
			Kr.SetAll(0);		
			Mr.SetSize(NModes(), NModes());
			Mr.SetAll(0);

			ConstMatrix<FEmaxDOF*CMSmaxDOF> temp_PhiKM, temp_Phi;
			ConstMatrix<CMSmaxDOF*CMSmaxDOF> addKM;
			ConstMatrix<FEmaxDOF*FEmaxDOF> temp_K, temp_M;


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
				e.StiffnessMatrix(temp_K);
				temp_K *= -1;
				e.EvalMff(temp_M,0);

				// Kr = Phi_CB^T K_el Phi_CB

				// temp_Phi contains ltg-map of Phi_CB
				temp_Phi.SetSize(sos_el, NModes());


				for (int j=1; j<=sos_el; j++)
				{
					for (int k=1; k<=NModes(); k++)
					{
						temp_Phi(j,k) = Phi_CB(ltg(j),k);
					}
				}

				Mult(temp_K, temp_Phi, temp_PhiKM);
				MultTp(temp_Phi,temp_PhiKM, addKM);
				Kr += addKM;
				Mult(temp_M, temp_Phi, temp_PhiKM);
				MultTp(temp_Phi,temp_PhiKM, addKM);
				Mr += addKM;
			}
			// ------- Add additional masses at nodes to mass matrix
			for (int i=1; i<=massnodes.Length(); i++)
			{
				Node& node = GetNode(massnodes(i));
				for (int k=1; k<=node.SOS(); k++)
				{
					for (int j=1; j<=NModes(); j++)
					{
						for (int l=1; l<=NModes(); l++)
						{
							Mr(j,l) += massvalues(i) * Phi_CB(node.Get(k),j) * Phi_CB(node.Get(k),l);
						}
					}
				}
			}

		}

		//+++++++++++++++++++++++
		// Case: Eigenmodes are computed, if reading was not successful
		//+++++++++++++++++++++++
		else // if (read_eigenmodes_from_file == 0) // Eigenmodes are computed here!
		{
			//UO().InstantMessageText("read_eigenmodes_from_file == 0");
			// Compute mass and stiffness matrix of unreduced system
			ComputeMassAndStiffnessMatrix(Msparse, Ksparse);

			BuildPhiCB_CMS(Phi_CB, nb_modesclassic, nb_modesstatic, sosfull, nb_dofstotal, nb_dofsstatic);

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//build Kr and Mr:

			GetMBS()->UO() << "build Mr...\n";

			double ts = -GetClockTime();
			Matrix mmid;
			Mult(Msparse,Phi_CB,mmid);
			MultTp(Phi_CB,mmid,Mr);

			//ofstream massr(path3D+mystr("\\massR.txt"));
			//massr << Mr << "\n";

			GetMBS()->UO() << "build Kr...\n";
			Mult(Ksparse,Phi_CB,mmid);
			MultTp(Phi_CB,mmid,Kr);

			ts += GetClockTime();

			GetMBS()->UO() << "done in" << ts << "seconds\n";
		}

		UO() << "Mr=" << Mr << "\n";
		UO() << "Kr=" << Kr << "\n";
		UO() << "Mr=" << Mr << "\n";
		UO() << "max Phi=" << Phi_CB.MaxNorm() << "\n";
		UO() << "max Msparse=" << Msparse.MaxNorm() << "\n";

		Matrix test = Mr;
		int rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Mass-matrix not invertable!!!\n";}

		test = Kr;
		rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Stiffness-matrix not invertable!!!\n";}
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set new initial vector!!!
	Vector xicopy = x_init; //length = 2*SOSRigid+1
	int initlen = x_init.Length()/2; //only position
	x_init.SetLen(2*NModes()+xicopy.Length());

	//position:
	for (int i = 1; i <= NModes(); i++)
	{
		x_init(i) = 0;
	}
	for (int i = 1; i <= initlen; i++)
	{
		x_init(i+NModes()) = xicopy(i);
	}

	//velocity:
	for (int i = 1; i <= NModes(); i++)
	{
		x_init(i+NModes()+initlen) = 0;
	}
	for (int i = initlen+1; i <= xicopy.Length(); i++)
	{
		x_init(i+2*NModes()) = xicopy(i);
	}
	//UO() << "x_init=" << x_init << "\n";


	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//intialize a few other things:

	//InitializeFFRFMatrices();

	if(!animation_mode)
	{
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
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template <class RIGID>
void CMSElement<RIGID>::GetI1(Vector& I1) //Shabana p. 209-211
{
	//[3 x 1]
	int n = Dim();
	I1.SetLen(n);
	I1.SetAll(0);
	Vector tempv;

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();

		e.GetI1(tempv);
		for (int j=1; j <= n; j++)
		{
			I1(j) += tempv(j);
		}
	}
	// attached mass
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		for (int j=1; j <= n; j++)
		{
			I1(j) += massvalues(i)*node.Pos()(j);
		}
	}
}

template <class RIGID>
double CMSElement<RIGID>::GetIkl(int k, int l)
{
	//[1 x 1]
	double v = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		v += e.GetIkl(k,l);
	}
	// attached mass
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		v += massvalues(i)*node.Pos()(k)*node.Pos()(l);
	}

	return v;
}

template <class RIGID>
void CMSElement<RIGID>::GetIbarkl(int k, int l, Vector& I1)
{
	//[n x 1] --> [nr x 1]
	int sosfull = SOSFull();
	int nmodes = NModes();

	ConstVector<CMSmaxDOF> add_I1;
	ConstMatrix<FEmaxDOF*CMSmaxDOF> local_Phi;
	ConstVector<FEmaxDOF> local_I1;
	I1.SetLen(NModes());
	I1.SetAll(0);

	//UO() << "Ibarkl" << k << l << ":\n";

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		// set element-local Phi_CB matrix
		GetLocalMatrix(Phi_CB, ltg, sos, local_Phi);
		e.GetIbarkl(k,l,local_I1);
		MultTp(local_Phi, local_I1, add_I1);
		I1 += add_I1;
	}
	// attached mass
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		// mass * u_bar(k) * shape(l)
		// using that shape(l) = 1 for degree of freedom node.Get(l) and zero for all other dofs
		for (int n=1; n<=nmodes; n++)
		{
			I1(n) += massvalues(i) * node.Pos()(k) * Phi_CB(node.Get(l),n);
		}
	}
}

template <class RIGID>
void CMSElement<RIGID>::GetSbar(Matrix& Sbar)
{
	// Sbar * Phi_CB
	int nmodes = NModes();
	//[3 x n] --> [3 x nr]
	ConstMatrix<CMSmaxDOF*3> add_Sbar;
	ConstMatrix<FEmaxDOF*CMSmaxDOF> local_Phi;
	ConstMatrix<FEmaxDOF*3> local_Sbar;
	Sbar.SetSize(Dim(), NModes());
	Sbar.FillWithZeros();

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		e.GetSbar(local_Sbar);
		// set element-local Phi_CB matrix
		GetLocalMatrix(Phi_CB, ltg, sos, local_Phi);


		// local_Sbar * local_Phi
		Mult(local_Sbar, local_Phi, add_Sbar);
		Sbar += add_Sbar;
	}
	// attached mass
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		// mass * u_bar(k) * shape(l)
		// using that shape(l) = 1 for degree of freedom node.Get(l) and zero for all other dofs
		for (int n=1; n<=nmodes; n++)
		{
			for (int l=1; l<=Dim(); l++)
			{
				Sbar(l, n) += massvalues(i)*Phi_CB(node.Get(l),n);
			}
		}
	}

}

template <class RIGID>
void CMSElement<RIGID>::GetSbarkl(int k, int l, Matrix& Sbar)
{
	int sosfull = SOSFull();
	int nmodes = NModes();
	Sbar.SetSize(nmodes, nmodes);
	Sbar.SetAll(0);

	ConstMatrix<CMSmaxDOF*CMSmaxDOF> add_Sbar;
	ConstMatrix<FEmaxDOF*CMSmaxDOF> local_SbarPhi, local_Phi;
	ConstMatrix<FEmaxDOF*FEmaxDOF> local_Sbar;

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		// set element-local Phi_CB matrix
		GetLocalMatrix(Phi_CB, ltg, sos, local_Phi);

		e.GetSbarkl(k,l,local_Sbar); //this matrix is not symmetric!!!!

		//local_SbarPhi = local_Sbar*local_Phi;
		//add_Sbar = local_Phi.GetTp()*local_SbarPhi;
		Mult(local_Sbar, local_Phi, local_SbarPhi);
		MultTp(local_Phi, local_SbarPhi, add_Sbar);
		Sbar += add_Sbar;
	}
	// attached mass
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		// mass * shape(k) * shape(l)
		// using that shape(k) = 1 for dof node.get(k) and shape(l) = 1 for degree of freedom node.Get(l) and zero for all other dofs
		for (int m=1; m<=nmodes; m++)
		{
			for (int n=1; n<=nmodes; n++)
			{
				Sbar(m,n) += massvalues(i) * Phi_CB(node.Get(k),m) * Phi_CB(node.Get(l), n);
			}
		}
	}
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// FFRF: compute integrals int_V \rho ubar_k * ubar_l dV
template <class RIGID>
double CMSElement<RIGID>::GetIntRhoUkUl(int k, int l, const Vector& xgloc, Vector& temp, int use_mthetatheta_full)
{
	double uklbar = IklS[k-1][l-1];
	if (use_mthetatheta_full)
	{
		Mult(SbarklS[k-1][l-1],xgloc,temp);
		uklbar += IbarklS[k-1][l-1]*xgloc + IbarklS[l-1][k-1]*xgloc + temp*xgloc;
	}
	return uklbar;
}

// FFRF: compute integrals d/dt int_V \rho ubar_k * ubar_l dV
// attention: d/dt int_V \rho ubar_k * ubar_l dV = uklbarp(k,l) + uklbarp(l,k)
template <class RIGID>
double CMSElement<RIGID>::GetIntRhoUkUlP(int k, int l, const Vector& xgloc, const Vector& xgploc, Vector& temp)
{
	double uklbarp;
	Mult(SbarklS[k-1][l-1],xgloc,temp);
	uklbarp = IbarklS[l-1][k-1]*xgploc + temp*xgploc;
	return uklbarp;
}



template <class RIGID>
void CMSElement<RIGID>::EvalM(Matrix& m, double t) 
{

	TMStartTimer(23);

	Matrix3D A = GetRotMatrix();
	Matrix3D Gbar = GetGbar();
	Matrix3D GbarT = GetGbarT();
	Matrix3D G = GetG();

	ConstVector<CMSmaxDOF> temp;

	int off = NModes(); //offset where rigid body entries start!
	xg.SetLen(NModes());
	for (int i=1; i <= NModes(); i++) xg(i) = XG(i);

	int use_mrr=1;
	int use_mff=1;
	int use_mrf=1;
	int use_mthetatheta0=1;
	int use_mthetatheta_full=1;
	int use_mrtheta0=1;
	int use_mrtheta_full=1;
	int use_mftheta=1;

	// clear mass matrix
	m.SetAll(0.);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRR = Unit(3,3) * mass
	if (use_mrr)
	{
		m(1+off,1+off) = GetMass();
		m(2+off,2+off) = GetMass();
		m(3+off,3+off) = GetMass();
	}


	//UO() << "Mass=" << GetMass() << "\n";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//mRtheta = -A*S_tbartilde*Gbar, [3x4]
	//S_tbartilde is skew symmetric 3x3 matrix of S_tbar, which is integral over rho*(x + u_fbar)

	if (use_mrtheta0)
	{
		//compute S_tbar:
		Vector3D S_tbar;
		GetStbar(S_tbar, xg, use_mrtheta_full);

		Matrix3D S_tbartilde;
		S_tbartilde.SetSkew(S_tbar);

		Matrix3D mRtheta = -1*(A*(S_tbartilde*Gbar));

		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NRotParam(); j++)
			{
				m(i+off,3+j+off) = mRtheta(i,j);
				m(3+j+off,i+off) = mRtheta(i,j);
			}
		}

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_Rf: = A*Sbar //very important term

	if (use_mrf)
	{
		for (int j = 1; j <= NModes(); j++)
		{
			Vector3D v(SbarS(1,j),SbarS(2,j),SbarS(3,j)); //coincides with rho*GetH() ...
			v = A*v;
			m(j,1+off) = v(1);
			m(j,2+off) = v(2);
			m(j,3+off) = v(3);
			m(1+off,j) = v(1);
			m(2+off,j) = v(2);
			m(3+off,j) = v(3);
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_theta = Gbar^T I_bar_theta_theta Gbar
	if (use_mthetatheta0)
	{
		Matrix3D uklbar(0.); //contains integrals for I_thetatheta

		//compute only symmetric part:
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = i; j <= Dim(); j++)
			{
				uklbar(i,j) = GetIntRhoUkUl(i,j,xg, temp, use_mthetatheta_full);
			}
		}

		Matrix3D Ibartt;
		ComputeIbarThetaTheta(Ibartt, uklbar);

		Matrix3D mthetatheta;
		mthetatheta = GbarT * Ibartt * Gbar;

		for (int i = 1; i <= NRotParam(); i++) //4 Euler parameters
		{
			for (int j = 1; j <= NRotParam(); j++)
			{
				m(3+off+i,3+off+j) = mthetatheta(i,j); 
			}
		}
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_theta_f = - Gbar^T * I_theta_f

	if (use_mftheta)
	{
		//I_theta_f = int (rho * ubar_tilde * S) dV
		ConstVector<CMSmaxDOF> I_theta_f[3];
		ComputeIbarThetaF(I_theta_f[0], I_theta_f[1], I_theta_f[2], xg, temp);

		Matrix m0(NModes(), NRotParam());
		for (int i = 1; i <= NRotParam(); i++) //4 Euler parameters / 3 Kardan angles
		{
			for (int j = 1; j <= NModes(); j++)
			{
				double val = 0;
				for (int k = 1; k <= Dim(); k++)
				{
					val += GbarT(i,k)*I_theta_f[k-1](j);
				}
				m(j,off+3+i) = val; 
				m(off+3+i,j) = val; 
			}
		}
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//m_ff
	if (use_mff)
	{
		for (int i=1; i <= NModes(); i++)
		{
			for (int j=1; j <= NModes(); j++)
			{
				m(i,j) = Mr(i,j);
			}
		}
	}
	TMStopTimer(23);


}; 

int printdamp = 0;

template <class RIGID>
void CMSElement<RIGID>::EvalF2(Vector& f, double t) 
{
	Vector f0 = f;
	Body3D::EvalF2(f,t); //no follower forces!

	ConstVector<CMSmaxDOF> temp; 
	int n = NModes();
	int dim = Dim();
	xg.SetLen(n);
	for (int i=1; i <= n; i++) xg(i) = XG(i);

	// Constraint for Euler parameters
	AddEPCqTterms(f);

	TMStartTimer(22);



	//// for flexible coordinates
	//if (!(FastStiffnessMatrix() == 3 && GetMBS()->IsJacobianComputation()))
	//{
	Mult(Kr,xg,temp);
	for (int i=1; i <= n; i++) 
	{
		f(i) -= temp(i);
	}
	//}

	if (Dr.Length() == n)
	{
		for (int i=1; i <= n; i++) 
		{
			f(i) -= XGP(i) * Dr(i);
		}

		if (!printdamp)
		{
			printdamp = 1;
			GetMBS()->UO() << "damping=" << Dr << "\n";
		}
	}

	
	if (!GetMBS()->IsJacobianComputation()) //$!JG 2011-03-01: Auskommentiert um Dämpfungsmatrix zu ermitteln
	{
		AddQuadraticVelocityVector(f, t, temp);
	}

	TMStopTimer(22);
}; 




template<class RIGID>
void CMSElement<RIGID>::AddQuadraticVelocityVector(Vector& f, double t, Vector& temp)
{
	int nm = NModes();
	int dim = Dim();
	int nrotparam = NRotParam();

	int use_Qvr0 = 1;
	int use_Qvr_full = 1;
	int use_Qvtheta0 = 1; //here must be some error!!!
	int use_Qvtheta_full = 1;
	int use_Qvf = 1;

	TMStartTimer(17);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//quadratic velocity vector, Shabana p.229:

	ConstVector<CMSmaxDOF> Qvf;
	Qvf.SetLen(nm);
	Qvf.SetAll(0);
	ConstVector<CMSmaxDOF> xgp;
	xgp.SetLen(nm);
	for (int i = 1; i <= nm; i++) xgp(i) = XGP(i);

	ConstVector<4> Qvtheta;
	Qvtheta.SetLen(NRotParam());
	Qvtheta.SetAll(0);

	Vector3D QvR(0.);;

	Vector3D Rp = GetRefVel();

	ConstVector<4> beta(NRotParam()), betap(NRotParam());
	GetBeta(beta);
	GetBetaP(betap);
	Matrix3D A = GetRotMatrix();
	Matrix3D Gbar = GetGbar();
	Matrix3D GbarT = GetGbarT();
	Matrix3D Gbarp = GetGbarp();
	Matrix3D GbarpT = Gbarp.GetTp();
	Vector3D Gbarp_betap = Gbarp * betap;


	// d omega / d theta
	Matrix3D DomegaDtheta, DomegaDthetaT;
	DomegaDtheta = GetDOmegaDTheta(betap);
	DomegaDthetaT = DomegaDtheta.GetTp();

	Vector3D omegabar = Gbar*betap;

	Matrix3D omegabar_tilde, omegabar_tilde2; // \tilde \bar omega and (\tilde \bar omega)^2
	omegabar_tilde.SetSkew(omegabar);
	omegabar_tilde2 = omegabar_tilde * omegabar_tilde;

	Vector3D S_tbar;
	GetStbar(S_tbar, xg, use_Qvr_full);
	Matrix3D S_tbar_tilde;
	S_tbar_tilde.SetSkew(S_tbar);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Qv_R:
	// for general parametrization of rigid body
	// Qv_R = -A*[omegabar_tilde^2*S_t_bar+2*omega_bar_tilde*Sbar*qp_f] + A Stbartilde GbarP thetaP
	// for Euler Parameters (Rigid3D) the last term vanishes

	Vector3D Sbarqpf(0.,0.,0.);	// = Sbar*qp_f
	if (use_Qvr_full)
	{
		for (int i = 1; i <= dim; i++)
		{
			for (int j = 1; j <= nm; j++)
			{
				Sbarqpf(i) += SbarS(i,j)*xgp(j); 
			}
		}
	}

	if (use_Qvr0)
	{
		// term A Stbartilde GbarP thetaP (vanishes for Rigid3D)
		QvR = A*S_tbar_tilde*Gbarp_betap;
		// terms -A*[omegabar_tilde^2*S_t_bar+2*omega_bar_tilde*Sbar*qp_f]
		QvR -= A*omegabar_tilde*(omegabar_tilde*S_tbar + 2.*Sbarqpf);


		f(nm+1) += QvR(1);
		f(nm+2) += QvR(2);
		f(nm+3) += QvR(3);

	}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//             general                                                special case of Eulerparameters, then Gbarp = - domega/dtheta, Gbarp thetap = 0
	//Qv_theta: - (GbarpT-domega/dtheta) Ibar_theta_theta omegabar       -2(GbarpT) Ibar_theta_theta omegabar
	//          - (GbarpT-domega/dtheta)*Ibar_thetaf*qpf                 -2(GbarpT)*Ibar_thetaf*qpf
	//          - GbarT*Ibar_theta_thetaP*omegabar											 - GbarT*Ibar_theta_thetaP*omegabar
	//          - GbarT Ithetatheta Gbarp thetap												 0 = GbarT Ithetatheta 0

	// -- compute I_theta_f = int (rho * ubar_tilde * S) dV
	ConstVector<CMSmaxDOF> I_theta_f[3];
	if (use_Qvtheta_full || use_Qvf )
	{
		ComputeIbarThetaF(I_theta_f[0], I_theta_f[1], I_theta_f[2], xg, temp);
	}

	Matrix3D Ibarttp;
	Matrix3D Ibartt;
	if (use_Qvtheta0)
	{
		// -- compute Ibar_theta_theta, Ibar_theta_thetaP
		Matrix3D uklbar; //contains integrals for I_thetatheta
		for (int i = 1; i <= dim; i++)
		{
			for (int j = i; j <= dim; j++) //+++++++++ NOTE: only loop from i to dim !!! ++++++++++++
			{
				uklbar(i,j) = GetIntRhoUkUl(i,j,xg,temp,use_Qvtheta_full);
			}
		}
		//compute Ibar_theta_theta
		ComputeIbarThetaTheta(Ibartt, uklbar);

		if (use_Qvtheta_full)
		{
			Matrix3D uklbarp; //contains integrals for I_thetatheta_p, only one part!
			for (int i = 1; i <= Dim(); i++)
			{
				for (int j = 1; j <= Dim(); j++)
				{
					uklbarp(i,j) = GetIntRhoUkUlP(i,j,xg,xgp,temp);
				}
			}
			//compute Ibar_theta_theta_p
			ComputeIbarThetaThetaP(Ibarttp, uklbar, uklbarp);
		}

		ConstVector<4> tmp4(nrotparam);
		//// -- term 1:
		Matrix3D fact;
		Qvtheta.SetLen(3);
		Qvtheta.SetAll(0.);

		//// - (GbarpT-domega/dtheta) Ibar_theta_theta omegabar
		Vector3D Ibarttom = Ibartt * omegabar;
		fact = DomegaDthetaT - GbarpT;
		Mult(fact, Ibarttom, Qvtheta);
		Qvtheta += tmp4;


		// -- term 2, term 3:
		if (use_Qvtheta_full)
		{
			temp.SetLen(3);
			temp(1) = I_theta_f[0]*xgp;
			temp(2) = I_theta_f[1]*xgp;
			temp(3) = I_theta_f[2]*xgp;

			////////////  - (GbarpT-domega/dtheta)*Ibar_thetaf*qpf
			fact = GbarpT - DomegaDthetaT;
			Mult(fact, temp, tmp4);
			Qvtheta -=	tmp4;

			////////////  - GbarT*Ibar_theta_thetaP*omegabar
			Mult(GbarT, Ibarttp * (omegabar), tmp4);
			Qvtheta -=	tmp4;
		}

		//// these terms vanish in case of Rigid3D -----



		////////// -- term 5:   - GbarT Ithetatheta Gbarp thetap
		Vector3D term5; 
		Mult(GbarT, Ibartt*Gbarp_betap, tmp4);
		Qvtheta -= tmp4;

		// end vanishing terms -----------------------

		for (int i=1; i<=nrotparam; i++)
		{
			f(nm+dim+i) += Qvtheta(i);
		}

	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Qv_f
	// for general parametrization of rigid body
	// Qv_f: int rho*(S^T(-omegabartilde^2*ubar-2*omegabartilde*upbar_f) - Ibarthetaf^T Gbarp thetap
	//      = sum_jk [- Ibar_kj*omegabartilde2_jk - sum_l Sbar_kl*omegabartilde2_kl*qf_l - sum_l 2*omegabartilde_jk*S_kl*qpf_l] - Ibarthetaf^T Gbarp thetap
	// in case of Eulerparameters (Rigid3D), the last term vanishes

	if (use_Qvf)
	{ 
		for (int i = 1; i <= nm; i++)
		{
			for (int j = 1; j <= dim; j++)
			{
				// TERM int rho*(S^T(-omegabartilde^2*ubar+2*omegabartilde*upbar_f)
				for (int k = 1; k <= dim; k++)
				{
					// theta theta --  correct
					Qvf(i) -= IbarklS[k-1][j-1](i)*omegabar_tilde2(j,k);
					for (int l = 1; l <= nm; l++)
					{
						Qvf(i) -= SbarklS[k-1][j-1](i,l)*(2*omegabar_tilde(k,j)*xgp(l) + omegabar_tilde2(j,k)*xg(l)) ;
					}
				}
			}

			// term -2* Sbar^T omegabartilde A^T Rp !!!!!!!!!
			Vector3D help = omegabar_tilde*A.GetTp()*Rp;
			Qvf(i) -= 2*(SbarS(1,i)*help(1) +SbarS(2,i)*help(2) +SbarS(3,i)*help(3)); 


			// TERM  + Ibarthetaf^T Gbarp thetap   !!!!!!!!!!
			Vector3D IbarT_thetaf_i(I_theta_f[0](i), I_theta_f[1](i), I_theta_f[2](i));
			Qvf(i) -= IbarT_thetaf_i*Gbarp_betap;

		}

		for (int i=1; i <= nm; i++)
		{
			f(i) += Qvf(i);
		}

	}


	TMStopTimer(17);
}

template <class RIGID>
void CMSElement<RIGID>::ComputeIbarThetaTheta(Matrix3D& Ibartt, Matrix3D& uklbar)
{

	Ibartt(1,1) = uklbar(2,2) + uklbar(3,3);
	Ibartt(2,2) = uklbar(1,1) + uklbar(3,3);
	Ibartt(3,3) = uklbar(1,1) + uklbar(2,2);

	Ibartt(1,2) = -uklbar(1,2); //only upper triangle of uklbar computed!!!
	Ibartt(1,3) = -uklbar(1,3);
	Ibartt(2,3) = -uklbar(2,3);
	Ibartt(2,1) = Ibartt(1,2);
	Ibartt(3,1) = Ibartt(1,3);
	Ibartt(3,2) = Ibartt(2,3);
}

template <class RIGID>
void CMSElement<RIGID>::ComputeIbarThetaThetaP(Matrix3D& Ibarttp, Matrix3D& uklbar, Matrix3D& uklbarp)
{
	Ibarttp(1,1) = 2.*uklbarp(2,2) + 2.*uklbarp(3,3);
	Ibarttp(2,2) = 2.*uklbarp(1,1) + 2.*uklbarp(3,3);
	Ibarttp(3,3) = 2.*uklbarp(1,1) + 2.*uklbarp(2,2);

	Ibarttp(1,2) = -(uklbarp(2,1) + uklbarp(1,2));
	Ibarttp(1,3) = -(uklbarp(3,1) + uklbarp(1,3));
	Ibarttp(2,3) = -(uklbarp(3,2) + uklbarp(2,3));
	Ibarttp(2,1) = Ibarttp(1,2);
	Ibarttp(3,1) = Ibarttp(1,3);
	Ibarttp(3,2) = Ibarttp(2,3);

}

template <class RIGID>
void CMSElement<RIGID>::ComputeIbarThetaF(Vector& I_theta_f0, Vector& I_theta_f1, Vector& I_theta_f2, Vector& xg, Vector& temp)
{
	Mult(SbarklS[2][1],xg,temp);
	I_theta_f0  = IbarklS[2][1];
	I_theta_f0 += temp;
	Mult(SbarklS[1][2],xg,temp);
	I_theta_f0 -= IbarklS[1][2];
	I_theta_f0 -= temp;

	Mult(SbarklS[0][2],xg,temp);
	I_theta_f1  = IbarklS[0][2];
	I_theta_f1 += temp;
	Mult(SbarklS[2][0],xg,temp);
	I_theta_f1 -= IbarklS[2][0];
	I_theta_f1 -= temp;

	Mult(SbarklS[1][0],xg,temp);
	I_theta_f2  = IbarklS[1][0];
	I_theta_f2 += temp;
	Mult(SbarklS[0][1],xg,temp);
	I_theta_f2 -= IbarklS[0][1];
	I_theta_f2 -= temp;
}

template <class RIGID>
int CMSElement<RIGID>::FastStiffnessMatrix() const 
{
	return 0; //usually take 3
}

template <class RIGID>
void CMSElement<RIGID>::StiffnessMatrix(Matrix& m) //fill in sos x sos components, m might be larger
{
	//negative stiffness matrix!!!!
	for (int i=1; i <= NModes(); i++)
	{
		for (int j=1; j <= NModes(); j++)
		{
			m(i,j) = -Kr(i,j);
		}
	}
	//UO() << "Kr=" << Kr << "\n";

}

template <class RIGID>
void CMSElement<RIGID>::GetIntRhoDuDq(Matrix& rhodudq)
{
	if (intrhodudqr.Getrows() == SOS())
	{
		rhodudq = intrhodudqr;
	}
	else
	{
	Matrix3D A = GetRotMatrix();
	Matrix3D Gbar = GetGbar();


	rhodudq.FillWithZeros();

	xg.SetLen(NModes());
	for (int i=1; i <= NModes(); i++) xg(i) = XG(i);

	int off = NModes();

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// int rho d pos / d qR
	//same as mRR/rho
	rhodudq(1+off,1) = mass;
	rhodudq(2+off,2) = mass;
	rhodudq(3+off,3) = mass;


	//+++++++++++++++++++++++++++++++++++
	// int rho d pos / d qtheta
	//
	//mRtheta = -A*S_tbartilde*Gbar, [3x4]
	//S_tbartilde is skew symmetric 3x3 matrix of S_tbar, which is integral over rho*(x + u_fbar)
	//same as mRtheta = A_theta*[I1+Sbar*qf]
	//compute S_tbar = I1S + Sbar*qf:
	Vector3D S_tbar;
	GetStbar(S_tbar, xg, 1);

	Matrix3D S_tbartilde;
	S_tbartilde.SetSkew(S_tbar);
	Matrix3D ASbarGbar = (A*(S_tbartilde*Gbar));

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NRotParam(); j++)
		{
			rhodudq(3+j+off,i) = -ASbarGbar(i,j);
		}
	}

	//+++++++++++++++++++++++++++++++++++
	// int rho d pos / d qf
	//same as m_Rf = A*Sbar

	for (int j = 1; j <= NModes(); j++)
	{
		Vector3D v(SbarS(1,j),SbarS(2,j),SbarS(3,j)); 
		v = A*v;
		rhodudq(j,1) = v(1);
		rhodudq(j,2) = v(2);
		rhodudq(j,3) = v(3);
	}

	intrhodudqr = rhodudq; // save rhodudq for fast access
	}
}

template <class RIGID>
void CMSElement<RIGID>::GetIntDuDq(Matrix& dudq)
{
	Matrix3D A = GetRotMatrix();
	Matrix3D Gbar = GetGbar();


	dudq.FillWithZeros();

	int off = NModes();

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// int d pos / d qR
	//same as mRR/rho
	dudq(1+off,1) = volume;
	dudq(2+off,2) = volume;
	dudq(3+off,3) = volume;


	//+++++++++++++++++++++++++++++++++++
	// int d pos / d qtheta
	//
	//mRtheta = -A*S_tbartilde*Gbar, [3x4]
	//S_tbartilde is skew symmetric 3x3 matrix of S_tbar, which is integral over rho*(x + u_fbar)
	//same as mRtheta/rho = A_theta*[I1+Sbar*qf]/rho
	if (volume == 0) 
	{
		UO() << "Warning: volume of CMS element is zero!\n  ==>Body load can not be applied!!!\n";
		volume = 1;
	}
	double rho0 = mass / volume;

	//compute S_tbar/rho = I1S/rho0 + Hr*qf:
	Vector3D S_tbar;

	for (int i = 1; i <= Dim(); i++)
	{
		S_tbar(i) = I1S(i)/rho0;
		for (int j = 1; j <= NModes(); j++)
		{
			S_tbar(i) += Hr(j,i)*XG(j); // here, Hr is used instead of Sbar = Hr * rho0
		}
	}

	// Here, S_tbartilde is 1/rho0 S_tbartilde
	Matrix3D S_tbartilde;
	S_tbartilde.SetSkew(S_tbar);
	Matrix3D ASbarGbar = (A*(S_tbartilde*Gbar));

	for (int i = 1; i <= Dim(); i++)
	{
		for (int j = 1; j <= NRotParam(); j++)
		{
			dudq(3+j+off,i) = -ASbarGbar(i,j);
		}
	}

	//+++++++++++++++++++++++++++++++++++
	// int d pos / d qf
	//same as m_Rf/rho = A*Sbar/rho = A*Hr

	for (int j = 1; j <= NModes(); j++)
	{
		Vector3D v(Hr(j,1),Hr(j,2),Hr(j,3)); 
		v = A*v;
		dudq(j,1) = v(1);
		dudq(j,2) = v(2);
		dudq(j,3) = v(3);
	}


	/*

	if (!(GetMBS()->IsJacobianComputation()))
	{

	//+++++++++++++++++++++++++++++++++++
	//same as m_Rf = A*Sbar:
	Matrix3D A = GetRotMatrix2D();

	for (int j = 1; j <= NModes(); j++)
	{
	Vector2D v(SbarS(1,j),SbarS(2,j)); //stimmt mit rho*GetH() überein!!!
	v = A*v;
	dudq(j,1) = v(1)/rho;
	dudq(j,2) = v(2)/rho;
	}
	}
	*/
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// position / velocity of node in the CMS element
template <class RIGID>
Vector3D CMSElement<RIGID>::GetNodePos(int i) const
{
	int flexdof = FlexDOF();
	const Node& nodei = GetNode(i);
	Vector3D prel = nodei.Pos();
	for (int k=1; k<=Dim(); k++)
		for (int j=1; j<=flexdof; j++)
		{
			int nodedofk = nodei.Get(k);
			prel(k) += Phi_CB(nodedofk, j) * XG(j);
		}
		return GetRefPos() + GetRotMatrix()*prel;
}
template <class RIGID>
Vector3D CMSElement<RIGID>::GetNodeLocPos(int i) const 
{
	int nm = NModes();
	Vector3D prel = nodes(i)->Pos();
	for (int k=1; k<=Dim(); k++)
		for (int j=1; j<=nm; j++)
		{
			int nodedofk = GetNode(i).Get(k);
			prel(k) += Phi_CB(nodedofk, j) * XG(j);
		}
		return prel;
}

template <class RIGID>
Vector3D CMSElement<RIGID>::GetNodeVel(int i) const 
{
	int nm = NModes();
	Vector3D prel = nodes(i)->Pos();
	Vector3D vrel = Vector3D(0, 0, 0);
	for (int k=1; k<=Dim(); k++)
	{
		for (int j=1; j<=nm; j++)
		{
			int nodedofk = GetNode(i).Get(k);
			prel(k) += Phi_CB(nodedofk, j) * XG(j);
			vrel(k) += Phi_CB(nodedofk, j) * XGP(j);
		}
	}

	return GetRefVel() + GetRotMatrix()*vrel + GetRotMatrixP()*prel;
}

template <class RIGID>
Vector3D CMSElement<RIGID>::GetNodePosD(int i) const 
{
	int nm = NModes();
	Vector3D prel = nodes(i)->Pos();
	for (int k=1; k<=Dim(); k++)
		for (int j=1; j<=nm; j++)
		{
			int nodedofk = GetNode(i).Get(k);
			prel(k) += Phi_CB(nodedofk, j) * XGD(j);
		}
		return GetRefPosD() + GetRotMatrixD()*prel;
}

template <class RIGID>
void CMSElement<RIGID>::GetdPosdqT(const Vector3D& ploc, Matrix& d) 
{
	d.SetSize(SOS(),3);
	d.SetAll(0);
	int nm = NModes();

	d(1+nm,1)=1; d(1+nm,2)=0; d(1+nm,3)=0;
	d(2+nm,1)=0; d(2+nm,2)=1; d(2+nm,3)=0;
	d(3+nm,1)=0; d(3+nm,2)=0; d(3+nm,3)=1;
	Matrix3D G;
	GetH3T(ploc, G);

	d.SetSubmatrix(G,Dim()+1+nm,1);
};

template <class RIGID>
void CMSElement<RIGID>::GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d)
{
	int nm = NModes();
	int sosrigid = SOSRigid();
	d.SetSize(sosrigid+nm,3);
	for (int i=1; i<=3; i++)
		for (int j=1; j<=3; j++)
			d(i+nm,j) = 0;

	Matrix3D H;
	GetH3T(vloc,H);
	d.SetSubmatrix(H,Dim()+1+nm,1);
}

template <class RIGID>
void CMSElement<RIGID>::GetdRotdqT(const Vector3D& ploc, Matrix& d)
{
	int nm = NModes();
	int sosrigid = SOSRigid();
	d.SetSize(sosrigid+nm,3);
	for (int i=1; i<=3; i++)
		for (int j=1; j<=3; j++)
			d(i+nm,j) = 0;

	d.SetSubmatrix(GetGT(),Dim()+1+nm,1);
}

template <class RIGID>
void CMSElement<RIGID>::GetNodedPosdqT(int node, Matrix& d)
{
	int nm = NModes();
	// Rigid body dPos/dq
	// Pos = pref + A prel
	d.SetSize(SOS(),Dim());
	d.FillWithZeros();

	Vector3D ploc;
	ploc = GetNodeLocPos(node);//.Pos();

	//d_pref/dq:
	d(nm+1,1) = 1;
	d(nm+2,2) = 1;
	d(nm+3,3) = 1;

	//d_A/dq*prel:
	Matrix3D G;
	GetH3T(ploc, G);
	for (int i=1; i<=NRotParam(); i++)
		for (int j=1; j<=3; j++)
			d(nm+Dim()+i,j) = G(i,j);

	// Transpose(A dprel/dq)
	// node.dof(k) is kth global, unreduced degree of freedom for node
	// in this node, shape functions for node.dofs are 1, all other shape functions are zero
	// d ( j, i ) = Sum_k [ A(i,k) * Phi_CB( node.dof(k), j) ]
	Matrix3D A = GetRotMatrix();
	for (int i=1; i<=Dim(); i++)
		for (int j=1; j<=nm; j++)
			for (int k=1; k<=Dim(); k++)
			{
				int nodedofk = GetNode(node).Get(k);
				d(j, i) += A(i,k) * Phi_CB(nodedofk, j);
			}

			//int nm = NModes();
			//Vector3D prel = nodes(i)->Pos();
			//for (int k=1; k<=Dim(); k++)
			//	for (int j=1; j<=nm; j++)
			//	{
			//		int nodedofk = GetNode(i).Get(k);
			//		prel(k) += Phi_CB(nodedofk, j) * XGD(j);
			//	}
			//return GetRefPosD() + GetRotMatrixD()*prel;


			// COPIED FROM THE TWO-DIMENSIONAL ELEMENT
			// DPOSDQ THERE IS NOT REDUCED, BUT OF ORIGINAL FLOATING-FRAME SIZE!!!!!
			/*
			int off = NModes();
			dpdqi.SetSize(NModes()+3,2);
			dpdqi.FillWithZeros();

			//compute d/dq (GetRefPos2D() + GetRotMatrix2D()*prel);

			//compute relative node position:
			Vector2D prel = Vector2D(nodes(node)->Pos().X()+XG(nodes(node)->Get(1)),nodes(node)->Pos().Y()+XG(nodes(node)->Get(2)));

			//d/dq (GetRefPos2D()):
			//rigid body motion:
			dpdqi(NModes()+1,1) = 1;
			dpdqi(NModes()+1,2) = 0;
			dpdqi(NModes()+2,1) = 0;
			dpdqi(NModes()+2,2) = 1;

			// d/dq (GetRotMatrix2D())*prel:
			double phi = GetAngle2D();
			double sphi = sin(phi);
			double cphi = cos(phi);
			dpdqi(NModes()+3,1) = -prel.X()*sphi-prel.Y()*cphi;
			dpdqi(NModes()+3,2) =  prel.X()*cphi-prel.Y()*sphi;

			//GetRotMatrix2D() * d/dq prel 
			Matrix3D A = GetRotMatrix2D();

			Vector2D tmp = A*Vector2D(1,0);
			dpdqi(nodes(node)->LTG(1),1) = tmp.X();
			dpdqi(nodes(node)->LTG(1),2) = tmp.Y();

			tmp = A*Vector2D(0,1);
			dpdqi(nodes(node)->LTG(2),1) = tmp.X();
			dpdqi(nodes(node)->LTG(2),2) = tmp.Y();
			*/
}

template <class RIGID>
void CMSElement<RIGID>::AddNodedPosdqTLambda(int node, const Vector3D& lambda, Vector& f)   // f += dpdq*lambda
{
	int nm = NModes();

	Vector3D ploc;
	ploc = GetNode(node).Pos();

	// Rigid body dPos/dq
	// Pos = pref + A prel

	//d_pref/dq:
	f(nm+1) += lambda(1);
	f(nm+2) += lambda(2);
	f(nm+3) += lambda(3);

	//d_A/dq*prel:
	Matrix3D G;
	GetH3T(ploc, G);
	for (int i=1; i<=NRotParam(); i++)
		for (int j=1; j<=3; j++)
			f(nm+Dim()+i) += lambda(j)*G(i,j);

	// Transpose(A dprel/dq)
	// node.dof(k) is kth global, unreduced degree of freedom for node
	// in this node, shape functions for node.dofs are 1, all other shape functions are zero
	// d ( j, i ) = Sum_k [ A(i,k) * Phi_CB( node.dof(k), j) ]
	Matrix3D A = GetRotMatrix();
	for (int i=1; i<=Dim(); i++)
		for (int j=1; j<=nm; j++)
			for (int k=1; k<=Dim(); k++)
			{
				int nodedofk = GetNode(node).Get(k);
				f(j) += lambda(i) * A(i,k) * Phi_CB(nodedofk, j);
			}

			// COPIED FROM THE TWO-DIMENSIONAL ELEMENT
			// DPOSDQ THERE IS NOT REDUCED, BUT OF ORIGINAL FLOATING-FRAME SIZE!!!!!
			/*
			int nm = NModes();

			//compute d/dq (GetRefPos2D() + GetRotMatrix2D()*prel);

			//compute relative node position:
			Vector2D prel = Vector2D(nodes(node)->Pos().X()+XG(nodes(node)->Get(1)),nodes(node)->Pos().Y()+XG(nodes(node)->Get(2)));

			//d/dq (GetRefPos2D()):
			//rigid body motion:

			f(nm+1) += lambda(1);
			f(nm+2) += lambda(2);

			// d/dq (GetRotMatrix2D())*prel:
			double phi = GetAngle2D();
			double sphi = sin(phi);
			double cphi = cos(phi);
			f(nm+3) += lambda(1)*(-prel.X()*sphi-prel.Y()*cphi);
			f(nm+3) += lambda(2)*( prel.X()*cphi-prel.Y()*sphi);

			//GetRotMatrix2D() * d/dq prel 
			Matrix3D A = GetRotMatrix2D();

			Vector2D tmp = A*Vector2D(1,0);
			f(nodes(node)->LTG(1)) += lambda(1)*tmp.X();
			f(nodes(node)->LTG(1)) += lambda(2)*tmp.Y();

			tmp = A*Vector2D(0,1);
			f(nodes(node)->LTG(2)) += lambda(1)*tmp.X();
			f(nodes(node)->LTG(2)) += lambda(2)*tmp.Y();
			*/
}




template <class RIGID>
void CMSElement<RIGID>::GetBeta(Vector& beta) const
{
	//GetMBS()->UO() << "CMSElement<RIGID>::GetBeta(Vector& beta) is not tested yet!\n";
	for (int i=1; i<= NRotParam(); i++)
	{
		beta(i) = XG(GetIndBeta(i));
	}
}

template <class RIGID>
void CMSElement<RIGID>::GetBetaP(Vector& betap) const
{
	//GetMBS()->UO() << "CMSElement<RIGID>::GetBetaP(Vector& betap) is not tested yet!\n";
	for (int i=1; i<= NRotParam(); i++)
	{
		betap(i) = XGP(GetIndBeta(i));
	}
}

template <class RIGID>
void CMSElement<RIGID>::GetBetaD(Vector& beta) const
{
	//GetMBS()->UO() << "CMSElement<RIGID>::GetBetaD(Vector& beta) is not tested yet!\n";
	for (int i=1; i<= NRotParam(); i++)
	{
		beta(i) = XGD(GetIndBeta(i));
	}
}

template <class RIGID>
void CMSElement<RIGID>::GetBetaPD(Vector& betap) const
{		
	//GetMBS()->UO() << "CMSElement<RIGID>::GetBetaPD(Vector& betap) is not tested yet!\n";
	for (int i=1; i<= NRotParam(); i++)
	{
		betap(i) = XGPD(GetIndBeta(i));
	}
}	
//----------------------------------------------------


template <class RIGID>
void CMSElement<RIGID>::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	BaseCMSElement<RIGID>::GetElementData(edc);

	ElementData ed;

	//ed.SetInt(nbmodes, "NBModes"); edc.Add(ed);
	//ed.SetInt(n_zeromodes, "NZeroModes"); edc.Add(ed);
	//ed.SetInt(solverparameters.solvertype, "EigenvalueSolvertype"); ed.SetToolTipText("Eigenvalue Solver: 1 .. Matlab (iterative), 2 .. HOTINT (direct)"); edc.Add(ed);
	//ed.SetDouble(solverparameters.tolerance, "EVSolverTolerance"); ed.SetToolTipText("Eigenvalue Solver tolerance (only for iterative solver)"); edc.Add(ed);
	//ed.SetInt(solverparameters.maxiterations, "EVSolverMaxIterations"); ed.SetToolTipText ("Eigenvalue Solver maximum iterations (only for iterative solver)"); edc.Add(ed);
}

template <class RIGID>
int CMSElement<RIGID>::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = BaseCMSElement<RIGID>::SetElementData(edc);

	//GetElemDataInt(mbs, edc, "NIModes", nimodes, 1);
	//GetElemDataInt(mbs, edc, "NZeroModes", n_zeromodes, 1);
	//GetElemDataInt(mbs, edc, "EigenvalueSolvertype", solverparameters.solvertype, 1);
	//GetElemDataDouble(mbs, edc, "EVSolverTolerance", solverparameters.tolerance, 1);
	//GetElemDataInt(mbs, edc, "EVSolverMaxIterations", solverparameters.maxiterations, 1);
	return rv;
}


template <class RIGID>
void CMSElement<RIGID>::MatlabExportF2(Vector& f, double t)
{ 
	//not tested, especially the three subroutines may still contain bugs
	mystr filename(mystr("C:\\Dokumente und Einstellungen\\sinwel\\Eigene Dateien\\cpp\\test_Matlaboutput\\test_Matlaboutput\\CMSElement_")+mystr(GetOwnNum())+mystr(".h"));
	CMFile file_export(filename,TFMwrite);
	mystr buffer("");
	//header
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ automatically created by CMSElement::MatlabExportF2\n"));
	file_export.RWSoleStr(mystr("#include <math.h>"));

	//function headers
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("//+ function headers\n"));
	file_export.RWSoleStr(mystr("int MatrixTranspose(double* result, double* source, int scol, int srow);"));
	file_export.RWSoleStr(mystr("int MatrixAdd(double* result, double* source1, int s1col, int s1row, double* source2, int s2col, int s2row);"));
	file_export.RWSoleStr(mystr("int MatrixMultiply(double* result, double* source1, int s1col, int s1row, double* source2, int s2col, int s2row);"));

	//constants
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("// constants\n"));

	// number of modes
	file_export.RWSoleStr(mystr("const int nmodes = ")+mystr(NModes())+mystr(";"));
	// number of dofs for the rigid body
	file_export.RWSoleStr(mystr("const int nrigiddofs = ")+mystr(3+NRotParam())+mystr(";"));
	// dimension
	file_export.RWSoleStr(mystr("const int dim = 3;"));
	// number of rotation dofs
	file_export.RWSoleStr(mystr("const int dof_rot = ")+mystr(NRotParam())+mystr(";"));
	// number of dofs for the unreduced system
	file_export.RWSoleStr(mystr("const int ndof_full = ")+mystr(SOSFull())+mystr(";"));
	// volume
	file_export.RWSoleStr(mystr("double volume = ")+mystr(volume)+mystr(";"));
	// mass
	file_export.RWSoleStr(mystr("double mass = ")+mystr(mass)+mystr(";"));

	// Vector I1S
	buffer = I1S.MakeString("double I1");
	file_export.RWSoleStr(buffer);

	// Matrix SbarS
	buffer = SbarS.MakeString("double Sbar");
	file_export.RWSoleStr(buffer);

	// double[3][3] IklS	
	buffer = mystr("double Ikl[dim*dim] = {");
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			if (j) buffer += mystr(", ");
			buffer += mystr(IklS[i][j]);
		}
		if (i<2) buffer += mystr(",\n");
	}
	buffer += mystr("};");
	file_export.RWSoleStr(buffer);

	// Vector[3][3] IbarklS	
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			mystr string_ibarkls = mystr("double Ibarkl")+mystr(i)+mystr(j);
			buffer = IbarklS[i][j].MakeString(string_ibarkls);
			file_export.RWSoleStr(buffer);
		}
	}

	// Matrix[3][3] Sbarkl	
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			mystr string_sbarkls = mystr("double Sbarkl")+mystr(i)+mystr(j);
			buffer = SbarklS[i][j].MakeString(string_sbarkls);
			file_export.RWSoleStr(buffer);
		}
	}

	// Phi_CB
	buffer = Phi_CB.MakeString("double Phi_CB");
	file_export.RWSoleStr(buffer);
	// Kr
	buffer = Kr.MakeString("double Kr");
	file_export.RWSoleStr(buffer);
	// Mr
	buffer = Mr.MakeString("double Mr");
	file_export.RWSoleStr(buffer);
	// Hr
	buffer = Hr.MakeString("double Hr");
	file_export.RWSoleStr(buffer);

	//functions
	file_export.RWSoleStr(mystr("//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"));
	file_export.RWSoleStr(mystr("// functions\n"));

	file_export.RWSoleStr(mystr("int GetIndRefPos(int i) {return nmodes+i;}"));
	file_export.RWSoleStr(mystr("int GetIndBeta(int i) {return nmodes+i+3;}"));

	file_export.RWSoleStr(mystr("#include \"functions_Rigid3D.h\""));
	file_export.RWSoleStr(mystr("#include \"functions_CMSElement.h\""));

	file_export.RWSoleStr(mystr("int MatrixTranspose(double* result, double* source, int scol, int srow)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  int i,j;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<srow;i++)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    for(j=0;j<scol;j++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("       result[j*srow+i] = source[i*scol+j];"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}"));  
	file_export.RWSoleStr(mystr(""));  

	file_export.RWSoleStr(mystr("int MatrixAdd(double* result, double* source1, int s1col, int s1row, double* source2, int s2col, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  if (s1col != s2col ) return 0;"));  
	file_export.RWSoleStr(mystr("  if (s1row != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  int i;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<(s1col*s1row);i++) { result[i] = source1[i]+source2[i]; }")); 
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}"));  
	file_export.RWSoleStr(mystr(""));  

	file_export.RWSoleStr(mystr("int MatrixMultiply(double* result, double* source1, int s1col, int s1row, double* source2, int s2col, int s2row)"));
	file_export.RWSoleStr(mystr("{"));  
	file_export.RWSoleStr(mystr("  if (s1col != s2row ) return 0;"));  
	file_export.RWSoleStr(mystr("  int i,j,k;"));  
	file_export.RWSoleStr(mystr("  for(i=0;i<s1row;i++)"));  
	file_export.RWSoleStr(mystr("  {"));  
	file_export.RWSoleStr(mystr("    for(j=0;j<s2col;j++)"));  
	file_export.RWSoleStr(mystr("    {"));  
	file_export.RWSoleStr(mystr("      result[i*s2col+j] = 0.0;"));  
	file_export.RWSoleStr(mystr("      for(k=0;k<s1col;k++)"));  
	file_export.RWSoleStr(mystr("      {"));  
	file_export.RWSoleStr(mystr("        result[i*s2col+j] += source1[i*s1col+k]*source2[k*s2col+j];"));  
	file_export.RWSoleStr(mystr("      }"));  
	file_export.RWSoleStr(mystr("    }"));  
	file_export.RWSoleStr(mystr("  }"));  
	file_export.RWSoleStr(mystr("  return 1;"));  
	file_export.RWSoleStr(mystr("}"));  
	file_export.RWSoleStr(mystr(""));  

}

//$ AH 2012-02: compute angular momentum of the body with p_ref as referenc
template <class RIGID>
Vector3D CMSElement<RIGID>::GetAngularMomentum(const Vector3D& p_ref) const
{
	if (FlexDOF() == 0)
	{
		// AH: todo - satz von steiner
		Vector3D omegabar = GetAngularVelLocal();
		return GetRotMatrix() * (Iphi*omegabar);
	}
	else
	{
		Vector3D D(0.,0.,0.);
		for (int i = 1; i <= NFFRFElements(); i++)
		{
			D += GetFFRFElement(i).GetAngularMomentum(p_ref);
		}
		return D;
	}
}

template <class RIGID>
double CMSElement<RIGID>::GetKineticEnergy() 
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
double CMSElement<RIGID>::GetPotentialEnergy()
{
	Vector q(SOS()), f(SOS());
	EvalF2(f, 0.);

	for (int i = 1; i <= SOS()-6; i++)
	{
		q(i) = -XG(i);
	}

	return f*q;
}

template class CMSElement<Rigid3D>;
template class CMSElement<Rigid3DKardan>;
//template class NodalConstraintCMS<Rigid3D>;
//template class NodalConstraintCMS<Rigid3DKardan>;