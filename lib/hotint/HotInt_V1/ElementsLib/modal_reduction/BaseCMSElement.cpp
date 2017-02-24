//#**************************************************************
//#
//# filename:             BaseCMSElement.cpp
//#
//# author:               Astrid Pechstein
//#
//# generated:						Okt. 2010
//# description:          Base class for Component Mode Synthesis
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
#include "linalgeig.h"

mystr path3D = "..\\..\\CMSvectors";


template <class RIGID>
void BaseCMSElement<RIGID>::AddFEMesh(const FEMeshInterface & femesh)
{
	// add nodes from FEMeshInterface
	for (int n=1; n<=femesh.NPoints(); n++)
	{
//		Node node(3, 1, femesh.Point(n)); //$ AD 2011-02-24: capsuled Nodes in FEMesh, Point(i) function split in Get & Set, 2D & 3D
		Node node(3, 1, femesh.GetPoint3D(n));

		AddNode(&node);
	}
	mbs->UO() << "number of nodes " << femesh.NPoints() << "\n";

	for (int el=1; el<=femesh.NElements(); el++)
	{
		//JG, AD: use function from FEMesh to convert FEElements to MBSElements

		//get new element, force domain=1, force CMSElementNr=1
		//*YV{
		femesh.Settings().domain = 1;
		femesh.Settings().generate_ffrf_elements = true;
		femesh.Settings().cms_element_number = GetOwnNum();
		Element * element_ptr = femesh.GetPtrToMBSElement_new(el); //get new element, force domain=1, 
		//*YV}
		if (element_ptr)
		{
			int auxelnr = GetMBS()->AddAuxElement(element_ptr);
			AddFFRFElement(auxelnr);

			delete element_ptr; //delete element, has been generated with new in "FEMesh::GetPtrToMBSElement_new"!
		}
	}
}



//this function is called after assembly and AFTER filysis
template <class RIGID>
void BaseCMSElement<RIGID>::Initialize() 
{

	ReferenceFrame3D<RIGID>::Initialize();

}



template <class RIGID>
const double& BaseCMSElement<RIGID>::GetXactFull(int i) const 
{
	//this function usually should not be called, only for evaluation with WriteSol()!!!
	//and with constraints!!!

	//multiply xg with Phi_CB, just i-th row

	//static double v = 0;
	static Vector v;//(Phi_CB.Getrows());
	v.SetLen(Phi_CB.Getrows());
	int ii = i;
	int off = 0;
	if (i > Phi_CB.Getrows()) 
	{
		ii -= Phi_CB.Getrows();
		off = SOSRigid();
	}

	v(ii) = 0;
	//double off = 0;
	for (int j = 1; j <= Phi_CB.Getcols(); j++)
	{
		//AH: wtf?
		//if (j == 1) off = 1; else off = 0;

		double xg;
		if (i > Phi_CB.Getrows()) xg = XGP(j);
		else xg = XG(j);
		v(ii) += Phi_CB(ii,j)*xg; //(XG(j+off));
	}
	return v(ii);
}

template <class RIGID>
double& BaseCMSElement<RIGID>::GetXactFull(int i)
{
	//this function usually should not be called!!!

	//multiply xg with Phi_CB, just i-th row

	//static double v = 0;
	static Vector v;//(Phi_CB.Getrows());
	v.SetLen(Phi_CB.Getrows());
	int ii = i;
	if (i > Phi_CB.Getrows()) ii -= Phi_CB.Getrows();

	v(ii) = 0;
	for (int j = 1; j <= Phi_CB.Getcols(); j++)
	{
		v(ii) += Phi_CB(ii,j)*XG(j);
	}
	return v(ii);
}

template <class RIGID>
const double& BaseCMSElement<RIGID>::GetDrawValueFull(int i) const 
{
	//multiply xgd with Phi_CB, just i-th row
	//int mode = (GetMBS()->TestCnt()%NIModes()) + 1+NBModes();


	static Vector v;//(Phi_CB.Getrows());
	v.SetLen(Phi_CB.Getrows());
	int ii = i;
	if (i > Phi_CB.Getrows()) ii -= Phi_CB.Getrows();

	v(ii) = 0;

	for (int j = 1; j <= Phi_CB.Getcols(); j++)
	{
		if (GetMBS()->GetIOption(116))
		{
			if (j == GetMBS()->GetDOption(100)) v(ii) += Phi_CB(ii,j);
		}
		else
		{
			v(ii) += Phi_CB(ii,j)*XGD(j);
		}
	}
	return v(ii);

	/*	
	static double v = 0; //does not work properly!!!
	for (int j = 1; j <= Phi_CB.Getcols(); j++)
	{
	v += Phi_CB(i,j)*XGD(j);
	if (j == mode) v += Phi_CB(i,j)*1e-1;
	}
	return v;*/
}

// link elements and compute ltg map for element -> xgfull
template <class RIGID>
void BaseCMSElement<RIGID>::ComputeElementLTG()
{
	int sosfull = SOSFull();
	int sos1 = 1; //starting position for position DOF
	int sos2 = sosfull+1; //starting position for veloctiy DOF

	// Link nodes into FFRF System
	//add IDs to nodes --> then add IDs to elements

	//add references node_dof --> local dof
	for (int i=1; i<=nodes.Length(); i++) 
	{
		//displacements
		for (int j=1; j <= GetNode(i).SOS(); j++)
		{
			GetNode(i).AddLTG(sos1++);
		}
		//velocities
		for (int j=1; j <= GetNode(i).SOS(); j++)
		{
			GetNode(i).AddLTG(sos2++);
		}
	}

	// Link Elements into FFRF System
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		e.PreAssemble(); //$JG 2011-09-12

		e.Initialize();

		//only for internal dof .. not used
		for (int j=1; j <= e.SOSowned(); j++)
		{
			e.AddLTG(sos1++);
		}
		for (int j=1; j <= e.SOSowned(); j++)
		{
			e.AddLTG(sos2++);
		}

		//add local element  
		//set flag, that ltg list is build this time
		TMBSElement t = e.GetType();
		e.AddType(TCMSflag);
		e.LinkToElements();
		e.SetType(t);
	}

}

template <class RIGID>
void BaseCMSElement<RIGID>::ComputeBoundaryDofIndexList(const TArray<int2>& fixednodes, int sosfull)
{
	boundarydof.SetLen(sosfull);
	for (int i=1; i <= boundarydof.Length(); i++)
	{
		boundarydof(i) = 0;
	}

	// add all node-local fixed dofs to boundary dof list, fixed dofs are overruled by boundary dofs and staticmode dofs!
	for (int i=1; i <= fixednodes.Length(); i++)
	{
		boundarydof(GetNode(fixednodes(i).Get(1)).Get(fixednodes(i).Get(2))) = 3;
	}

	// add all dofs corresponding to additional static modes to the boundary dof list -- they are used additional to boundary dofs!
	for (int i=1; i <= staticmode_nodes.Length(); i++)
	{
		for (int noden=1; noden<=staticmode_nodes(i)->Length(); noden++)
		{
			boundarydof(GetNode((*staticmode_nodes(i))(noden)).Get(1)) = 2;
			boundarydof(GetNode((*staticmode_nodes(i))(noden)).Get(2)) = 2;
			boundarydof(GetNode((*staticmode_nodes(i))(noden)).Get(3)) = 2;
		}
	}

	// add all boundary nodes to boundary dof list
	for (int i=1; i <= boundarynode.Length(); i++)
	{
		boundarydof(GetNode(boundarynode(i).Get(1)).LTG(boundarynode(i).Get(2))) = 1; //fix the specified components of node
	}
	// add all element-local boundary dofs to boundary dof list
	for (int i=1; i <= boundarydofelem.Length(); i++)
	{
		boundarydof(GetFFRFElement(boundarydofelem(i).Get(1)).LTG(boundarydofelem(i).Get(2))) = 1;
	}

}

template <class RIGID>
void BaseCMSElement<RIGID>::ResortLTG(int sosfull)
{
	// resort ... sorted -> original
	IVector resort(sosfull*2);
	// resort2 ... original -> sorted
	IVector resort2(sosfull*2);

	resort.SetLen(0);
	for (int i=1; i <= sosfull; i++)
	{
		if (boundarydof(i) == 1) resort.Add(i); //add classical boundary nodes
	}
	for (int i=1; i <= sosfull; i++)
	{
		if (boundarydof(i) == 2) resort.Add(i); //add additional static boundary nodes
	}
	for (int i=1; i <= sosfull; i++)
	{
		if (boundarydof(i) == 3) resort.Add(i); //add fixed boundary nodes
	}
	for (int i=1; i <= sosfull; i++)
	{
		if (boundarydof(i) == 0) resort.Add(i); //add nodes which are not at boundary
	}
	//for velocities:
	for (int i=1; i<=sosfull; i++) resort.Add(resort(i)+sosfull);

	//UO() << "resort=" << resort << "\n";

	resort2.SetLen(sosfull);
	for (int i=1; i <= sosfull; i++)
	{
		resort2(resort(i)) = i;
	}

	//for velocities:
	for (int i=1; i<=sosfull; i++) resort2.Add(resort2(i)+sosfull);

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		for (int j=1; j <= e.LTGlength(); j++)
		{
			e.LTG(j) = resort2(e.LTG(j));
		}
	}

	//resort node DOF for constraints:
	for (int i=1; i<=nodes.Length(); i++) 
	{
		//resort node DOF for constraints
		for (int j=1; j <= GetNode(i).LTGLength(); j++)
		{
			GetNode(i).LTG(j) = resort2(GetNode(i).LTG(j));
		}
	}

}

template <class RIGID>
void BaseCMSElement<RIGID>::ComputeNZeromodes()
{
	//count number of fixed nodes (boundarydof is set)
	int n_fixednodes = 0;
	for (int i=1; i<=boundarydof.Length(); i++)
	{
		if (boundarydof(i)) {	n_fixednodes++;}
	}
	// number of zeromodes: in case of no fixednodes, n_zeromodes is 6 --> free-free modes are used
	// otherwise, n_zeromodes is 0, body has to be statically fixed to reference frame by fixednodes!!
	if (n_zeromodes < 0)
	{
		if (n_fixednodes == 0)
		{
			n_zeromodes = 6;
			mbs->UO() << "CMSElement: Computing free-free modes\n";
		}
		else if (n_fixednodes < 6)
		{
			n_zeromodes = 6 - n_fixednodes;
			mbs->UO() << "CMSElement: " << n_fixednodes << " degrees of freedom are fixed, " << n_zeromodes << " modes are eliminated!\n";
			mbs->UO() << "CMSElement: Attention: X = inv(K_II) K_IB cannot be computed correctly!\n";
		}
		else
		{
			n_zeromodes = 0;
			mbs->UO() << "CMSElement: Body is fixed w.r.t. reference frame\n";
		}
	}
	else
	{
		mbs->UO() << "CMSElement: Body is has to be fixed w.r.t. reference frame, otherwise inv(K_II) K_IB cannot be computed!\n";
	}
}

// compute/read eigenmodes, build classic Craig-Bampton matrix Phi_CB which contains modes columnwise
// degree of freedom management for Craig-Bampton matrix Phi_CB:
//
//             bd.class  bd.stat   modes
//            [ I_c = id      0        0 ]    boundary dofs for classical modes
// Phi_CB =   [     0   I_s = shapes   0 ]    boundary dofs for additional static modes
//          	[     0         0        0 ]    fixed dofs
//            [    X_c       X_s      Phi]    flexible dofs
//
// matrix X_c = - K_II^(-1) K_IB I_c =  - K_II^(-1) K_IB 
// matrix X_s = - K_II^(-1) K_IB I_s 
template <class RIGID>
void BaseCMSElement<RIGID>::BuildPhiCB_CMS(Matrix& Phi_CB, int nb_modesclassic, int nb_modesstatic, int sosfull, int nb_dofs, int nb_dofsstatic)
{
	int nb = NBModes();
	int nf = sosfull-nb_dofs; //number of internal (flexible) dof - unreduced
	int nr = NIModes();

	Matrix Phi_r(nr+n_zeromodes,nf); // = Phi_r
	Vector eigval(nr+n_zeromodes);


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// compute eigenvalues:
	// Matlab: sparse computation of nr smallest eigenvalues, 6 zero-eigenvalues are excluded!
	// Hotint: full computation of nr smallest eigenvalues, 6 zero-eigenvalues are excluded!
	//         Matrices K_II, M_II and K_IB are stored as full matrices!!!

	Matrix X; //K_II^(-1)*K_IB

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//sparse computation with MATLAB:
	if (nimodes)
	{
		if (solverparameters.solvertype==1)
		{
			Compute_Phi_X_Matlab(Msparse, Ksparse, nb_dofs, nf, nr, X, eigval, Phi_r);
		}
		else if (solverparameters.solvertype==0)
		{
			Compute_Phi_X(Msparse, Ksparse, nb_dofs, nf, nr, X, eigval, Phi_r);
		}
		else if (solverparameters.solvertype==2)
		{
			Compute_Phi_X_Sparse(Msparse, Ksparse, nb_dofs, nf, nr, X, eigval, Phi_r);
		}
		else
		{
			mbs->UO() << "CMSElement::solverparameters.solvertype Eigenvalue-Solver set to " << solverparameters.solvertype 
				<< ", only 1 (use iterative matlab) or 0 (use direct HOTINT solver) allowed!\n";
			return;
		}
	}
	GetMBS()->UO() << "Eigenfrequencies in Hz:\n";
	for (int ii=1; ii <= eigval.Length(); ii++)
	{
		GetMBS()->UO() << "mode " << ii << ": " << sqrt(fabs(eigval(ii)))/(2.*MY_PI) << "\n";
	}

	if (SOS() > CMSmaxDOF)
		mbs->UO().InstantMessageText("Error in CMSElement: SOS greater than maximum count CMSmaxDOF! Set CMSmaxDOF to greater value in file CMSElement.h");

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//build Phi_CB matrix:
	//
	// degree of freedom management for Craig-Bampton matrix Phi_CB:
	//
	//             bd.class  bd.stat   modes
	//            [ I_c = id      0        0 ]    boundary dofs for classical modes
	// Phi_CB =   [     0   I_s = shapes   0 ]    boundary dofs for additional static modes
	//          	[     0         0        0 ]    fixed dofs
	//            [    X_c       X_s      Phi]    flexible dofs
	//
	// matrix X_c = - K_II^(-1) K_IB I_c =  - K_II^(-1) K_IB 
	// matrix X_s = - K_II^(-1) K_IB I_s 

	Phi_CB.SetSize(sosfull, NModes());
	Phi_CB.SetAll(0);

	// Set Phi_CB_BB = I_c for the classical static modes
	for(int i=1; i <= nb_modesclassic; i++)
	{
		Phi_CB(i, i) = 1;
	}

	// AP: Here, insert other additional static modes for Phi_CB(nb_modesclassic + i, nb_modesclassic + j),  i = 1...nb_staticnodes,  j = 1...nb_staticmodes
	// Set Phi_CB_BB = I_s = shapes for the additional static modes
	for (int j=1; j<=nb_modesstatic; j++)
	{
		for (int ni=1; ni<=staticmode_nodes(j)->Length(); ni++)
		{
			int nodenum = (*staticmode_nodes(j))(ni);
			Node& node = GetNode(nodenum);
			switch (staticmode_type(j))
			{
			case TCMSStaticModeTranslation:
				{
					Vector3D direction = (*staticmode_loccoords(j))(1);
					Phi_CB(node.Get(1), j+nb_modesclassic) = direction(1);
					Phi_CB(node.Get(2), j+nb_modesclassic) = direction(2);
					Phi_CB(node.Get(3), j+nb_modesclassic) = direction(3);
					break;
				}
			case TCMSStaticModeRotation:
				{
					Vector3D direction = (*staticmode_loccoords(j))(1);
					Vector3D center = (*staticmode_loccoords(j))(2);
					Vector3D diff = node.Pos() - center;
					// direction of rotation = direction of axis times diff vector
					Vector3D rotdir = direction.Cross(diff);
					Phi_CB(node.Get(1), j+nb_modesclassic) = rotdir(1);
					Phi_CB(node.Get(2), j+nb_modesclassic) = rotdir(2);
					Phi_CB(node.Get(3), j+nb_modesclassic) = rotdir(3);
					break;
				}
			default:
				{
					mbs->UO() << "CMSElement::DoModalAnalysis: Static mode " << j << " has unknown type " << staticmode_type(j) << "\n";
				}
			}
		}
	}

	// Fill in Phi_CB_IB = -K_II^(-1) K_IB for the classical static modes j=1..nb_modesclassic
	for(int i=1; i <= nf; i++)
	{
		for(int j=1; j <= nb_modesclassic; j++)
		{
			Phi_CB(nb_dofs+i, j) = X(i,j);
		}
	}
	// Fill in Phi_CB_IB = -K_II^(-1) K_IB I_s = -K_II^(-1) K_IB Phi_CB_BB for the additional static modes j=1+nb_modesclassic .. nb
	for (int i=1; i<=nf; i++)
	{
		for (int j=1; j<=nb_modesstatic; j++)
		{
			double xij = 0;
			for (int k=1; k<=nb_dofsstatic; k++)
			{
				xij += X(i,k+nb_modesclassic) * Phi_CB(k+nb_modesclassic,j+nb_modesclassic);
			}
			Phi_CB(nb_dofs+i, j+nb_modesclassic) = xij;
		}
	}
	// fill dynamic modes obtained from eigenvalue problem into Phi_CB_II
	for(int i=1; i <= nf; i++)
	{
		for(int j=1; j <= nr; j++)
		{
			Phi_CB(nb_dofs+i, nb+j) = Phi_r(j+n_zeromodes,i);
		}
	}

	// write PHI_CB to specified file
	if (UseEigenmodeFile())
	{
		ofstream writeEVtofile(EVfile);
		writeEVtofile.precision(20);
		for (int j=1; j<=nbmodes; j++)
		{
			writeEVtofile << "0\n";
		}
		for (int j=1+nbmodes; j<=NModes(); j++)
		{
			writeEVtofile << sqrt(fabs(eigval(n_zeromodes+j-nbmodes)))/(2.*MY_PI) << "\n";
		}
		for (int i=1; i<=sosfull; i++)
		{
			for (int j=1; j<=NModes(); j++)
			{
				writeEVtofile << Phi_CB(i,j) << "\t";
			}
			writeEVtofile << "\n";
		}
		writeEVtofile.close();
	}

}

template <class RIGID>
void BaseCMSElement<RIGID>::ComputeMassAndStiffnessMatrix(SparseMatrix& Msparse, SparseMatrix& Ksparse)
{
	int sosfull = SOSFull();
	Msparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
	Msparse.FillWithZeros();
	Ksparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
	Ksparse.FillWithZeros();
	Matrix tempm;

	UO() << "assembling mass matrix and stiffness matrix ... \n";
	//UO() << "allocated mem in K=" << Ksparse.GetLAlloc()*8 << "bytes\n";



	for (int i=1; i<=NFFRFElements(); i++) 
	{
		if (NFFRFElements() > 500 && i%100 == 0) 
		{
			UO() << "\rElement " << i << " of " << NFFRFElements() << " ...";
		}

		Element& e = GetFFRFElement(i);

		int sos = e.FlexDOF();

		const	TArray<int>& ltg = e.GetLTGArray();

		tempm.SetSize(e.FlexDOF(), e.FlexDOF());
		tempm.SetAll(0);
		e.StiffnessMatrix(tempm);
		tempm *= -1;
		Ksparse.AddMatrix(ltg,ltg,sos,sos,tempm);

		e.EvalMff(tempm,0);
		Msparse.AddMatrix(ltg,ltg,sos,sos,tempm);
	}


	// ------- Add additional masses at nodes to mass matrix
	for (int i=1; i<=massnodes.Length(); i++)
	{
		Node& node = GetNode(massnodes(i));
		for (int k=1; k<=node.SOS(); k++)
			Msparse(node.Get(k), node.Get(k)) += massvalues(i);
	}
}
//$ DR 2011-05-27: changed initial value also for CMS to 0, before it was 1
// Normalization of the computed Eigenvectors, normalization mode = 0 (max(abs(v))=1) or 1 (v'*v = 1)
template <class RIGID>
void BaseCMSElement<RIGID>:: NormalizeEigenvectors(int normalization_mode, int nf, int nr, Matrix& eigmodes)
{
	if(normalization_mode == 0) // normalize such that max(abs(v)) = 1
	{
		for(int j=1; j <= nr+n_zeromodes; j++)
		{
			double norm = 0;
			for(int i=1; i <= nf; i++)
			{
				norm = Maximum(norm, fabs(eigmodes(j,i)));			
			}
			if (norm == 0) norm = 1;
			norm = 1./norm;
			for(int i=1; i <= nf; i++)
			{
				eigmodes(j,i) *= norm;
			}
		}
		UO(UO_LVL_ext) << "eigenmodes scaled such that max(abs(v)) = 1 \n";
	}

	else if(normalization_mode == 1) // normalize such that v'*v = 1
	{
		double sign = 1;
		int foundsign = 0;
		for(int j=1; j <= nr+n_zeromodes; j++)
		{
			double norm = 0;
			for(int i=1; i <= nf; i++)
			{
				if (!foundsign && eigmodes(j,i) != 0) 
				{
					sign = Sgn(eigmodes(j,i));
					foundsign = 1;
				}
				norm += Sqr(eigmodes(j,i));
			}
			if (norm == 0) norm = 1;
			norm = sign*1./sqrt(norm);
			for(int i=1; i <= nf; i++)
			{
				eigmodes(j,i) *= norm;
			}
			foundsign = 0;
			sign = 1;
		}
		UO(UO_LVL_ext) << "eigenmodes scaled such that v'*v = 1 \n";
	}
	return;
}



template <class RIGID>
void BaseCMSElement<RIGID>:: Compute_Phi_X(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, 
																					 Matrix& X, Vector& eigval, Matrix& eigmodes)
{
	Matrix K_II, M_II;
	Matrix K_IB;

	UO() << "CMSElement: Direct computation of eigenmodes ... \n";

	K_II.CopyFrom(Ksparse,nb+1,nb+1,nb+nf,nb+nf);
	M_II.CopyFrom(Msparse,nb+1,nb+1,nb+nf,nb+nf);
	K_IB.CopyFrom(Ksparse,nb+1,1,nb+nf,nb);

	// compute matrix X = - K_II^(-1) K_IB
	Vector help(nf), help2(nf); // reuse allocated data space
	X.SetSize(nf,nb);

	if (n_zeromodes==0)
	{
		// K_II has no kernel, inverse can be computed
		for (int i=1; i<=nb; i++)
		{
			// help is negative ith column of K_IB
			for (int j=1; j<=nf; j++)
				help(j) = -K_IB(j,i);
			// help is ith column of -K_II^-1 K_IB
			K_II.SolveLapack(help);
			for (int j=1; j<=nf; j++)
				X(j,i) = help(j);
		}
	}
	else
	{
		// K_II inverse cannot be computed!
		X.SetAll(0.);
	}




	if (nr==0) return;

	// -- -- compute generalized eigenvalues
	// solve gen. eigenvalue problem, eigenvectors are stored rowwise in all_eigmodes
	// eigenvalues are sorted by size, smallest EV first
	Vector work(4*nf);
	Matrix all_eigmodes;
	Vector all_eigvals(nf);
	all_eigmodes = K_II;

	int info = LapackGenEVPSPD(nf, &all_eigmodes(1,1), &M_II(1,1), &all_eigvals(1), &work(1), work.Length());

	// store eigenvectors in Matrix eigmodes
	for (int i = 1; i <= nr+n_zeromodes; i++)
	{
		eigval(i) = all_eigvals(i);
	}
	for(int i=1; i <= nf; i++)
	{
		for(int j=1; j <= nr+n_zeromodes; j++)
		{
			eigmodes(j,i) = all_eigmodes(j,i);
		}
	}


	NormalizeEigenvectors(this->GetMBS()->GetSolSet().eigsolv_normalization_mode, nf, nr, eigmodes);	//$ DR 2011-05-27

	// normalize eigenvectors
	//double sign = 1;
	//int foundsign = 0;
	//for(int j=1; j <= nr+n_zeromodes; j++)
	//{
	//	double norm = 0;
	//	for(int i=1; i <= nf; i++)
	//	{
	//		if (!foundsign && eigmodes(j,i) != 0) 
	//		{
	//			sign = Sgn(eigmodes(j,i));
	//			foundsign = 1;
	//		}
	//		norm += Sqr(eigmodes(j,i));
	//	}
	//	if (norm == 0) norm = 1;
	//	norm = sign*1./sqrt(norm);
	//	UO() << "norm = " << norm << "\n";
	//	for(int i=1; i <= nf; i++)
	//	{
	//		eigmodes(j,i) *= norm;
	//	}
	//	foundsign = 0;
	//	sign = 1;
	//}

	UO() << "Eigenvalues(Hz)=";
	for (int i=1; i <= eigval.Length(); i++)
	{
		UO() << sqrt(fabs(eigval(i)))/(2.*MY_PI) << "Hz, ";
	}
	UO() << "\n";

}

template <class RIGID>
void BaseCMSElement<RIGID>:: Compute_Phi_X_Sparse(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, 
																									Matrix& X, Vector& eigval, Matrix& eigmodes)
{
	X.SetSize(nf,nb);
	SparseMatrix *K_II, *M_II;
	if (nb > 0)
	{
		Matrix K_IB;

		UO() << "CMSElement: LOBPCG computation of eigenmodes ... \n";

		K_II = new SparseMatrix();
		M_II = new SparseMatrix();
		K_II->CopyFrom(Ksparse,nb+1,nb+1,nb+nf,nb+nf);
		M_II->CopyFrom(Msparse,nb+1,nb+1,nb+nf,nb+nf);
		K_IB.CopyFrom(Ksparse,nb+1,1,nb+nf,nb);

		SparseInverse *invK=0;
		if (n_zeromodes==0)
		{
			invK = new SparseInverse(*K_II);
			invK->Factorize();
		}

		// compute matrix X = - K_II^(-1) K_IB
		Vector help(nf), help2(nf); // reuse allocated data space

		if (n_zeromodes==0)
		{
			// K_II has no kernel, no zero eigenmodes -> inverse can be computed!
			for (int i=1; i<=nb; i++)
			{
				//// help is negative ith column of K_IB
				for (int j=1; j<=nf; j++)
					help(j) = -K_IB(j,i);
				// help is ith column of -K_II^-1 K_IB
				invK->Apply(help);
				for (int j=1; j<=nf; j++)
					X(j,i) = help(j);
			}
		}
		else
		{
			// inverse K_II cannot be computed!!
			// set X = 0!!
			X.SetAll(0.);
		}
		delete invK;
	}
	else
	{
		K_II = &Ksparse;
		M_II = &Msparse;
	}


	if (nr==0) return;




	SparseEigenvalueSolver EV(this->GetMBS());
	EV.Set(K_II, M_II, solverparameters.maxiterations, solverparameters.tolerance, 1); 
	if (solverparameters.use_precond)
		EV.SetPreconditioner(solverparameters.lambda_precond);

	TArray<int> unconstraineddofs(nf);
	for (int i=1; i<=nf; i++)
		unconstraineddofs(i) = i;
	EV.ComputeEigenModes(nr+n_zeromodes, eigval, eigmodes, unconstraineddofs);

	NormalizeEigenvectors(this->GetMBS()->GetSolSet().eigsolv_normalization_mode, nf, nr, eigmodes);	//$ DR 2011-05-27

	//// normalize eigenvectors
	//double sign = 1;
	//int foundsign = 0;
	//for(int j=1; j <= nr+n_zeromodes; j++)
	//{
	//	double norm = 0;
	//	for(int i=1; i <= nf; i++)
	//	{
	//		if (!foundsign && eigmodes(j,i) != 0) 
	//		{
	//			sign = Sgn(eigmodes(j,i));
	//			foundsign = 1;
	//		}
	//		norm += Sqr(eigmodes(j,i));
	//	}
	//	if (norm == 0) norm = 1;
	//	norm = sign*1./sqrt(norm);
	//	for(int i=1; i <= nf; i++)
	//	{
	//		eigmodes(j,i) *= norm;
	//	}
	//	foundsign = 0;
	//	sign = 1;
	//}

	UO() << "Eigenvalues(Hz)=";
	for (int i=1; i <= eigval.Length(); i++)
	{
		UO() << sqrt(fabs(eigval(i)))/(2.*MY_PI) << "Hz, ";
	}
	UO() << "\n";

	if (nb>0)
	{
		delete K_II;
		delete M_II;
	}

}

template <class RIGID>
void BaseCMSElement<RIGID>:: Compute_Phi_X_Matlab(SparseMatrix& Msparse, SparseMatrix& Ksparse, int nb, int nf, int nr, 
																									Matrix& X, Vector& eigval, Matrix& eigmodes)
{

	SparseMatrix K_II, M_II;
	Matrix K_IB;

	UO() << "CMSElement: Matlab computation of eigenmodes ... \n";

	K_II.CopyFrom(Ksparse,nb+1,nb+1,nb+nf,nb+nf);
	M_II.CopyFrom(Msparse,nb+1,nb+1,nb+nf,nb+nf);
	K_IB.CopyFrom(Ksparse,nb+1,1,nb+nf,nb);

	SparseInverse *invK=0;
	if (n_zeromodes==0)
	{
		invK = new SparseInverse(K_II);
		invK->Factorize();
	}

	// compute matrix X = - K_II^(-1) K_IB
	Vector help(nf); // reuse allocated data space
	X.SetSize(nf,nb);
	for (int i=1; i<=nb; i++)
	{
		// help is negative ith column of K_IB
		for (int j=1; j<=nf; j++)
			help(j) = -K_IB(j,i);
		// help is ith column of -K_II^-1 K_IB
		invK->Apply(help);
		for (int j=1; j<=nf; j++)
			X(j,i) = help(j);
	}


	if (nr==0) return;

	SparseEigenvalueSolver EV(this->GetMBS());
	EV.Set(&K_II, &M_II, solverparameters.maxiterations, solverparameters.tolerance, 1); 
	//Vector eval(nr+n_zeromodes);
	TArray<int> unconstraineddofs(nf);
	for (int i=1; i<=nf; i++)
		unconstraineddofs(i) = i;
	EV.ComputeEigenModesMatlab(path3D,nr+n_zeromodes, eigval, eigmodes, unconstraineddofs, 0);

	// store eigenvectors in Matrix eigmodes
	// plus write eigenvalues and modes to output file
	//if (UseEigenmodeFile())
	//{
	//	writeEVfile.precision(20);
	//	mbs->UO() << "Eigenmodes written to file \"" << EVfile << "\"\n";
	//	for (int i = 1; i <= nr+n_zeromodes; i++)
	//	{
	//		writeEVfile << eigval(i) << "\n";
	//	}
	//	for(int i=1; i <= nf; i++)
	//	{
	//		for(int j=1; j <= nr+n_zeromodes; j++)
	//		{
	//			writeEVfile << eigmodes(j,i) << "\t";
	//		}
	//		writeEVfile << "\n";
	//	}
	//}


	//}
	//else // eigenmode file
	//{

	//	for (int i = 1; i <= nr+n_zeromodes; i++)
	//	{
	//		infile >> eigval(i);
	//	}
	//	for(int i=1; i <= nf; i++)
	//	{
	//		for(int j=1; j <= nr+n_zeromodes; j++)
	//		{
	//			infile >> eigmodes(j,i);
	//		}
	//	}

	//}

	NormalizeEigenvectors(this->GetMBS()->GetSolSet().eigsolv_normalization_mode, nf, nr, eigmodes);	//$ DR 2011-05-27II

	//normalize eigenvectors, such that they are the same as internally computed eigenvectors:
	//double sign = 1;
	//int foundsign = 0;
	//for(int j=1; j <= nr+n_zeromodes; j++)
	//{
	//	double norm = 0;
	//	for(int i=1; i <= nf; i++)
	//	{
	//		if (!foundsign && eigmodes(j,i) != 0) 
	//		{
	//			sign = Sgn(eigmodes(j,i));
	//			foundsign = 1;
	//		}
	//		norm += Sqr(eigmodes(j,i));
	//	}
	//	if (norm == 0) norm = 1;
	//	norm = sign*1./sqrt(norm);
	//	UO() << "norm = " << norm << "\n";
	//	for(int i=1; i <= nf; i++)
	//	{
	//		eigmodes(j,i) *= norm;
	//	}
	//	foundsign = 0;
	//	sign = 1;
	//}
	UO() << "Eigenfrequencies(Hz)=";
	for (int i=1; i <= eigval.Length(); i++)
	{
		UO() << sqrt(fabs(eigval(i)))/(2.*MY_PI) << "Hz = " << sqrt(fabs(eigval(i))) << " rad/s\n";
	}
	UO() << "\n";

	delete invK;

}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class RIGID>
void BaseCMSElement<RIGID>::ComputeMass() 
{
	mass = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		//e.Initialize(); //already done
		mass += e.GetMass();
	}
	for (int i=1; i<=massnodes.Length(); i++)
	{
		mass += massvalues(i);
	}

}



template <class RIGID>
void BaseCMSElement<RIGID>::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ReferenceFrame3D::GetElementData(edc);

	ElementData ed;

	ed.SetInt(nbmodes, "NBModes"); edc.Add(ed);
	ed.SetInt(nimodes, "NIModes"); edc.Add(ed);		// DR 2011-05-11: added this line
	ed.SetInt(n_zeromodes, "NZeroModes"); edc.Add(ed);
	ed.SetInt(solverparameters.solvertype, "EigenvalueSolvertype"); ed.SetToolTipText("Eigenvalue Solver: 1 .. Matlab (iterative), 2 .. HOTINT (direct)"); edc.Add(ed);
	ed.SetDouble(solverparameters.tolerance, "EVSolverTolerance"); ed.SetToolTipText("Eigenvalue Solver tolerance (only for iterative solver)"); edc.Add(ed);
	ed.SetInt(solverparameters.maxiterations, "EVSolverMaxIterations"); ed.SetToolTipText ("Eigenvalue Solver maximum iterations (only for iterative solver)"); edc.Add(ed);
}

template <class RIGID>
int BaseCMSElement<RIGID>::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = ReferenceFrame3D::SetElementData(edc);

	GetElemDataInt(mbs, edc, "NBModes", nbmodes, 1);	// DR 2011-05-11: added this line
	GetElemDataInt(mbs, edc, "NIModes", nimodes, 1);
	GetElemDataInt(mbs, edc, "NZeroModes", n_zeromodes, 1);
	GetElemDataInt(mbs, edc, "EigenvalueSolvertype", solverparameters.solvertype, 1);
	GetElemDataDouble(mbs, edc, "EVSolverTolerance", solverparameters.tolerance, 1);
	GetElemDataInt(mbs, edc, "EVSolverMaxIterations", solverparameters.maxiterations, 1);
	return rv;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// NodalConstraintCMS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class RIGID>
void NodalConstraintCMS<RIGID>::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	NodalConstraint::GetElementData(edc);

	ElementData ed;
	ed.SetInt(cms_el, "CMSElement"); ed.SetToolTipText("Number of CMS element"); edc.Add(ed);
}

template <class RIGID>
int NodalConstraintCMS<RIGID>::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = NodalConstraint::SetElementData(edc);

	GetElemDataInt(mbs, edc, "CMSElement", cms_el, 1);
	return rv;
}

template <class RIGID>
void NodalConstraintCMS<RIGID>::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG

	Vector3D nodepos = GetCMSElement().GetNodePos(NodeNum(1));
	Vector3D nodevel = GetCMSElement().GetNodeVel(NodeNum(1));

	if (MaxIndex()==3)
	{
		if (!IsVelocityConstraint())
		{
			for (int j=1; j<=loccoords.Length(); j++)
				f(j) -= (nodepos(loccoords(j))-ground(j)->Evaluate(t));
		}
		else
		{ //velocity constraints:
			for (int j=1; j<=loccoords.Length(); j++)
				f(j) -= (nodevel(loccoords(j))-ground(j)->Evaluate(t));
		}
	}
	else
	{ //velocity constraints:
		if (!IsVelocityConstraint())
		{
			for (int j=1; j<=loccoords.Length(); j++)
				f(j) -= (nodevel(loccoords(j)));
		}
		else
		{ //velocity constraints:
			for (int j=1; j<=loccoords.Length(); j++)
				f(j) -= (nodevel(loccoords(j))-ground(j)->Evaluate(t));
		}
	}

};


template <class RIGID>
void NodalConstraintCMS<RIGID>::EvalF2(Vector& f, double t)
{
	TMStartTimer(27);


	ConstMatrix<CMSElement<RIGID>::CMSmaxDOF*3> dposdqT;
	GetCMSElement().GetNodedPosdqT(NodeNum(1), dposdqT);
	Vector3D nodepos = GetCMSElement().GetNodePos(NodeNum(1));
	Vector3D nodevel = GetCMSElement().GetNodeVel(NodeNum(1));

	//if (!mbs->IsJacobianComputation())
	//{
	//	mbs->UO() << "time " << t << ", pos = " << nodepos << ", vel = " << nodevel << "\n";
	//}

	if (!UsePenaltyFormulation())
	{
		// Term Cq^T lambda
		for (int i=1; i<=SOS(); i++)
			for (int j=1; j<=loccoords.Length(); j++)
			{
				f(i) -= dposdqT(i,loccoords(j))*XG(2*SOS()+j);
			}
	}
	else // penalty: spring+damping
	{		
		for (int i=1; i<=SOS(); i++)
			for (int j=1; j<=loccoords.Length(); j++)
			{
				f(i) -= GetPenaltyStiffness()*(nodepos(loccoords(j))-ground(j)->Evaluate(t))*dposdqT(i,loccoords(j));
				f(i) -= damping_coeff*nodevel(loccoords(j))*dposdqT(i,loccoords(j));
			}
	}
	TMStopTimer(27);
}

template <class RIGID>
void NodalConstraintCMS<RIGID>::DrawElement() 
{
	if (GetDrawSizeScalar() == 0) return;

	mbs->SetColor(GetCol());
	Vector3D dir(0.,0.,0.);
	for (int i=1; i<=loccoords.Length(); i++)
	{
		Vector3D offset(0.,0.,(i-1)*GetDrawSizeScalar());
		Vector3D p = GetCMSElement().GetNodePosD(NodeNum(1));
		dir.Set(0.,0.,0);
		if (loccoords(i) <=3)
		{
			dir(loccoords(i)) = 1.;
			mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
		}
		else
		{
			dir(3) = 1.;
			mbs->DrawCone(p + (loccoords(i)-2)*(-GetDrawSizeScalar()*0.75)*dir,p+(loccoords(i)-3)*(-GetDrawSizeScalar()*0.75)*dir,(GetDrawSizeScalar()*0.6),6,1);
		}
	}

};

template <class RIGID>
void NodalConstraintCMS<RIGID>::LinkToElements()
{
	TArray<int> storeltg(IS());
	for (int i=1; i <= IS(); i++)
	{
		storeltg.Add(LTG(i));
	}
	LTGreset();

	// add all SOS dofs from the CMS element
	const BaseCMSElement<RIGID>& cms = GetCMSElement();
	//Position(first SOS) and Velocity (second SOS):
	for (int i=1; i <= 2*cms.SOS(); i++)
	{
		AddLTG(cms.LTG(i));
	}
	// implicit dofs corresponding to constraint Lagrange parameters
	for (int i=1; i <= IS(); i++)
	{
		AddLTG(storeltg(i));
	}
}


	// ReadEigenModesFile reads eigenmodes from file columnwise into eigenmodemat
template <class RIGID>
void BaseCMSElement<RIGID>::ReadEigenModesFile(Matrix& eigenmodemat, int& successful)
{
	successful = 0;

	if (!GetMBS()->GetSolSet().eigsolv_reuse_last_eigvec) //$!DR 2012-01-27: bugfix, otherwise program crash in release mode
	{	
		mbs->UO() << "CMSElement: eigsolv_reuse_last_eigvec == 0 --> compute new eigenvectors\n";
		return;		
	}
	
	// check if eigenmodes can be read from file
	ifstream infile(EVfile);
	if (!infile.good())
	{
		mbs->UO() << "CMSElement: could not read eigenmode file \"" << EVfile << "\"\n";
		infile.close();
		UO().InstantMessageText("Could not read eigenmode file, see Info in Output-Window");
		return;
	}
	int sosfull = this->SOSFull();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// read Eigenvectors to matrix
	eigenmodemat.SetSize(sosfull, NModes());
	// Vector contains eigenvalues
	Vector eigval;
	eigval.SetLen(NModes());

	for (int j=1; j<=NModes(); j++)
	{
		infile >> eigval(j);
	}
	for (int i=1; i<=sosfull; i++)
	{
		for (int j=1; j<=NModes(); j++)
		{
			infile >> eigenmodemat(i,j);
		}
	}	

	infile.close();
	successful = 1;

}


// Load Eigenmodes into Data Manager
// ATTENTION: these are essentially unit vectors locally on the element
//            these flexible coordinates must be written into the correct local-to-global-positions!!!
//            do not forget about rigid body degrees of freedom of the CMSElement
//            do not forget that there are other elements as well (LTG necessary)!!
template <class RIGID>
void BaseCMSElement<RIGID>::EigenModesToDataManager()
{
	// number of Modes which is loaded into data manager = FlexDOF()
	// NModes is not correct for GCMSElement, correct for CMSElement
	ElementDataContainer edc;
	for (int j=1; j<=FlexDOF(); j++)
	{
		// solution vector - is essentially a unit vector with 1 in the LTG(j)-th place
		Vector solvec;
		solvec.SetLen(GetMBS()->GetSystemSize());
		// set flexible dofs: only mode number j is used
		solvec.SetAll(0.);
		solvec(LTG(j)) = 1;
		// set rigid body position dofs, as they are used in x_init
		// velocity dofs are all zero for this visualization, therefore not addressed here
		for (int i=1; i<=SOSRigid(); i++)
		{
			// ith rigid dof at position i+FlexDOF()
			solvec(LTG(i+FlexDOF())) = x_init(i+FlexDOF());
		}
		ElementData ed; 
		ed.SetVector(solvec.GetVecPtr(),solvec.GetLen(),mystr("SV")+mystr(j));
		edc.Add(ed);

		// no data vector
		double *datavec = 0;
		ed.SetVector(datavec,0,mystr("DV")+mystr(j));
		edc.Add(ed);
	}
	mbs->UO().CallWCDriverFunction(20,FlexDOF(),0,&edc);
}

template class BaseCMSElement<Rigid3D>;
template class BaseCMSElement<Rigid3DKardan>;
template class NodalConstraintCMS<Rigid3D>;
template class NodalConstraintCMS<Rigid3DKardan>;