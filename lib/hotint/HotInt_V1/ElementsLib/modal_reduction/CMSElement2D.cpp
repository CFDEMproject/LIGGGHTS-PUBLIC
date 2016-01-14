//#**************************************************************
//#
//# filename:             CMSElement2D.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						20. April  2006
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
 

#include "windows.h" //for shellexecute
#include "rigid2d.h"
#include "material.h"
#include "Node.h"
#include "femathhelperfunctions.h"
#include "femeshinterface.h"
#include "kinematicpairs.h"
#include "referenceFrame2D.h"
#include "cmselement2d.h"
//#include "graphicsconstants.h"
#include "elementdataaccess.h"
#include "solversettings_auto.h"
#include "lapack_routines.h"
#include "myfile.h"
#include "linalgeig.h"
#include "plate2d.h"

const int CMSmin_sparse_size = 72; //could be determined from maximum element size ...

mystr path = "C:\\Dokumente und Einstellungen\\gerstm\\Eigene Dateien\\cpp\\HotIntNew\\CMSvectors";
//mystr path = "C:\\hotint\\CMSvectors";

//this function is called after assembly and AFTER DoModalAnalysis
void CMSElement2D::Initialize() 
{
	ReferenceFrame2D::Initialize();

	if (NModes() != 0)
	{
		//UO() << "ERROR: Zero modes in component mode synthesis; Maybe the command DoModalAnalysis(int nimodes) has been forgotten?\n";

		Matrix tempm;

		//stored functions:
		GetSbarkl(1,2,Sbar_tilde);
		GetSbarkl(2,1,tempm);
		Sbar_tilde -= tempm; //Sbar12-Sbar21


		GetSbar(SbarS);
		GetIbarkl(1,1,Ibar11S);
		GetIbarkl(1,2,Ibar12S);
		GetIbarkl(2,1,Ibar21S); 
		GetIbarkl(2,2,Ibar22S);
		GetI1(I1S);

		I1122S = (GetIkl(1,1)+GetIkl(2,2));

		//Sbar_tildeSM.CopyFrom(Sbar_tilde);
		SbarSM.CopyFrom(SbarS);

		//++++++++++++++++++++++++++++++++++++++++++++++++++
		//compute reduced H matrix:
		Hr.SetSize(SOSFull(), Dim());
		Hr.FillWithZeros();

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();
			const	TArray<int>& ltg = e.GetLTGArray();

			e.GetH(tempm);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= Dim(); k++)
				{
					Hr(ltg(j),k) += tempm(j,k);
				}
			}
		}
		Hr = Phi_CB.GetTp()*Hr;
	}

};


//link elements, ltg-elements, compute nbmodes, compute M and K, modal analysis, 
//  transformation matrix, modal matrices, x_init!
void CMSElement2D::DoModalAnalysis(const TArray<int2>& fixednodes)
{
	UO() << "+++++++++++++++++\nDo Modal Analysis\n+++++++++++++++++\n";
	UO() << "n-FFRF-elements=" << NFFRFElements() << "\n";

	UO() << "compute reference numbers\n";

	//link elements + nodes!!!
	int sosfull = SOSFull();
	int sos1 = 1;
	int sos2 = sosfull+1;

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

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		
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

		//UO() << "ltg=" << e.GetLTGArray() << "\n";
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//compute necessary boundary nodes (DOF):
	UO() << "compute boundary modes\n";

	boundarydof.SetLen(sosfull);
	for (int i=1; i <= boundarydof.Length(); i++)
	{
		boundarydof(i) = 0;
	}

	for (int i=1; i <= boundarynode.Length(); i++)
	{
		boundarydof(GetNode(boundarynode(i)).LTG(1)) = 1; //fix all components of node
		boundarydof(GetNode(boundarynode(i)).LTG(2)) = 1;
	}

	for (int i=1; i <= boundarydofelem.Length(); i++)
	{
		//boundarydof.Add(GetFFRFElement(boundarydofelem(i).Get(1)).LTG(boundarydofelem(i).Get(2)));
		int eln = boundarydofelem(i).Get(1);
		int lc = boundarydofelem(i).Get(2);
		boundarydof(GetFFRFElement(eln).LTG(lc)) = 1;
	}

	nbmodes = 0;
	for (int i=1; i <= boundarydof.Length(); i++)
	{
		if (boundarydof(i) != 0) nbmodes++;
	}
	//UO() << "b-dof=" << boundarydof << "\n";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//resort ltg according to (B) and (I) nodes:

	IVector resort(sosfull*2);
	IVector resort2(sosfull*2);

	resort.SetLen(0);
	for (int i=1; i <= sosfull; i++)
	{
		if (boundarydof(i) != 0) resort.Add(i); //add boundary nodes
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
	//UO() << "resort2=" << resort2 << "\n";

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);

		for (int j=1; j <= e.LTGlength(); j++)
		{
			e.LTG(j) = resort2(e.LTG(j));
		}
		//UO() << "ltg_sort=" << e.GetLTGArray() << "\n";
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

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//size of stiffness and mass-matrix: sosfull x sosfull
	int nb = NBModes();
	int nf = sosfull-nb; //number of internal (flexible) dof - unreduced
	int nr = NIModes();

	UO() << "sos_full=" << sosfull << "\n";
	UO() << "nb=" << nb << ", nf=" << nf << ", nr=" << nr << "\n";


	Matrix Phi_r(nf,nr); // = Phi_r
	Vector eigval(nr);

	Matrix Phi_CBsave, Ksave, Msave;

	int usematlab = 1;

	if (usematlab)
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//sparse computation with MATLAB:
		Msparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
		Msparse.FillWithZeros();
		Ksparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
		Ksparse.FillWithZeros();
		Matrix tempm;

		UO() << "assembling mass matrix and stiffness matrix ... \n";
		//UO() << "allocated mem in K=" << Ksparse.GetLAlloc()*8 << "bytes\n";


		int issym = 1;

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			if (NFFRFElements() > 500 && i%100 == 0) 
			{
				UO() << "\rElement " << i << " of " << NFFRFElements() << " ...";
			}

			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			e.StiffnessMatrix(tempm);
			tempm.MakeSymmetric();
			if (!tempm.IsSymmetric(1e-6)) issym = 0;
			Ksparse.AddMatrix(ltg,ltg,sos,sos,tempm);

			e.EvalMff(tempm,0);
			tempm.MakeSymmetric();
			if (!tempm.IsSymmetric(1e-12)) issym = 0;
			Msparse.AddMatrix(ltg,ltg,sos,sos,tempm);
		}
		Msparse.EliminateZeroEntries(1e-15);
		Ksparse.EliminateZeroEntries(1e-15);
		if (NFFRFElements() > 500)
		{
			UO() << "\n";
		}
		//UO() << "allocated mem in K=" << Ksparse.GetLAlloc()*8 << "bytes\n";
		if (!issym) UO() << "Mass and stiffness matrix are not symmetric!\n";

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//extract K_II and M_II matrix

		SparseMatrix K_II, M_II;
		Matrix K_IB;
		Matrix X; //K_II^(-1)*K_IB

		UO() << "CMSElement: writing to matlab ... \n";

		K_II.CopyFrom(Ksparse,nb+1,nb+1,nb+nf,nb+nf);
		M_II.CopyFrom(Msparse,nb+1,nb+1,nb+nf,nb+nf);
		K_IB.CopyFrom(Ksparse,nb+1,1,nb+nf,nb);


		if (!GetMBS()->GetSolSet().eigsolv_reuse_last_eigvec) //reuse eigenvectors
		{
			{ //encapsulate!!!
				//ensure that stream is closed afterwards
				ofstream wait(path+mystr("\\wait.txt"));
				wait << "0\n";

				ofstream matout(path+mystr("\\matout.m"));
				ofstream M_IIdat(path+mystr("\\M_IId.dat"));
				ofstream K_IIdat(path+mystr("\\K_IId.dat"));
				ofstream K_IBdat(path+mystr("\\K_IBd.dat"));

				M_II.PrintToMatlabSparse(M_IIdat);
				K_II.PrintToMatlabSparse(K_IIdat);
				K_IB.PrintToMatlabSparse(K_IBdat);

				matout << "cd('" << path << "')\n";
				matout << "load M_IId.dat\n";
				matout << "M_II=spconvert(M_IId);\n";
				matout << ";\n";

				matout << "load K_IId.dat\n";
				matout << "K_II=spconvert(K_IId);\n";
				matout << ";\n";

				matout << "load K_IBd.dat\n";
				matout << "K_IB=spconvert(K_IBd);\n";
				matout << ";\n";

				matout << "neig=" << nr << "\n"; 
				matout << "comp_eig_CMS(K_II,M_II,K_IB," << nr << ");\n";
				matout << "exit\n"; 
				matout << "exit\n"; 
			}

			UO() << "compute eigenvectors in matlab ... \n";
			ShellExecute(NULL, "open", "matlab",	
				" /r matout", 
				//" -nosplash  /r matout", 
				path.c_str(), SW_HIDE);
				//path.c_str(), SW_SHOW);


			int goon = 0;
			while (!goon)
			{
				{
					Sleep(1000);
					ifstream wait(path+mystr("\\wait.txt"));
					int test;
					wait >> test;
					if (test == 1) goon = 1;
					//UO() << "not ready\n";
				}
			}
			UO() << "Eigenvectors computed in MATLAB!!!\n";
		}
		else
		{
			UO() << "*****************\nEigenvectors reused from last computation!\n*****************\n";
		}

		//GetMBS()->InstantMessageText("Compute eigenvectors with MATLAB\n AFTERWARDS press OK!");


		ifstream matin(path+mystr("\\matin.txt"));
		UO() << "read eigenvectors from matlab ...\n";
		for (int i=1; i <= nr; i++)
		{
			matin >> eigval(i);
		}

		int cnt = 0;
		for(int i=1; i <= nf; i++)
		{
			if (cnt > 20000 && nr*nf >= 50000) 
			{
				UO() << "\rread EV:" << (100*i)/nf << "%";
				cnt = 0;
			}
			cnt += nr;
			for(int j=1; j <= nr; j++)
			{
				matin >> Phi_r(i, j);
			}
		}
		UO() << "\r";

		//normalize eigenvectors, such that they are the same as internally computed eigenvectors:
		double sign = 1;
		int foundsign = 0;
		for(int j=1; j <= nr; j++)
		{
			double norm = 0;
			for(int i=1; i <= nf; i++)
			{
				if (!foundsign && Phi_r(i,j) != 0) 
				{
					sign = Sgn(Phi_r(i,j));
					foundsign = 1;
				}
				norm += Sqr(Phi_r(i, j));
			}
			if (norm == 0) norm = 1;
			norm = sign*1./sqrt(norm);
			for(int i=1; i <= nf; i++)
			{
				Phi_r(i, j) *= norm;
			}
			foundsign = 0;
			sign = 1;
		}

		X.SetSize(K_IB.Getrows(),K_IB.Getcols());
		int krows = X.Getrows();
		int kcols = X.Getcols();

		cnt = 0;
		for(int i=1; i <= krows; i++)
		{
			if (cnt > 50000 && krows*kcols >= 200000) 
			{
				UO() << "\rread X_IB:" << (100*i)/krows << "%";
				cnt = 0;
			}
			cnt += kcols;
			for(int j=1; j <= kcols; j++)
			{
				matin >> X(i, j);
			}
		}
		UO() << "\r";

		//UO() << "eigval1=" << eigval << "\n";
		//UO() << "Phi_r1=" << Phi_r << "\n";
		UO() << "Eigenvalues(Hz)=";
		for (int i=1; i <= eigval.Length(); i++)
		{
			UO() << sqrt(eigval(i))/(2.*MY_PI) << "Hz, ";
		}
		UO() << "\n";


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Phi_CB matrix:

		IVector fixed2(nb);
		for (int i=1; i <= nb; i++) fixed2(i) = 0;

		for (int i=1; i <= fixednodes.Length(); i++) 
		{
			fixed2(GetNode(fixednodes(i).Get(1)).Get(fixednodes(i).Get(2))) = 1;
		}

		IVector fixedresort(nb);
		IVector freeref;

		int nbnew = 0;
		for (int i=1; i <= nb; i++) 
		{
			if (fixed2(i))
			{
				fixedresort(i) = 0;
			}
			else
			{
				nbnew++;
				fixedresort(i) = nbnew; 
				freeref.Add(i);
			}
		}
		int cnt2 = nbnew;
		for (int i=1; i <= nb; i++) 
		{
			if (fixed2(i))
			{
				cnt2++;
				fixedresort(i) = cnt2; 
			}
		}

		nbmodes = nbnew;
		//GetMBS()->UO() << "nb_old=" << nb << "\n";
		//GetMBS()->UO() << "nb_new=" << nbnew << "\n";


		Phi_CB.SetSize(sosfull, NModes());
		Phi_CB.SetAll(0);

		for(int i=1; i <= nbnew; i++)
		{
			Phi_CB(i, i) = 1;
		}

		for(int i=1; i <= nf; i++)
		{
			for(int j=1; j <= nbnew; j++)
			{
				Phi_CB(nb+i, j) = X(i,freeref(j));
			}
		}

		for(int i=1; i <= nf; i++)
		{
			for(int j=1; j <= nr; j++)
			{
				Phi_CB(nb+i, nbnew+j) = Phi_r(i,j);
			}
		}

		//resort element DOF for fixed dofs
		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);

			for (int j=1; j <= e.LTGlength(); j++)
			{
				if (e.LTG(j) <= nb)
				{
				  e.LTG(j) = fixedresort(e.LTG(j));
				}
			}
		}

		//resort node DOF for fixed dofs:
		for (int i=1; i<=nodes.Length(); i++) 
		{
			//resort node DOF for constraints
			for (int j=1; j <= GetNode(i).LTGLength(); j++)
			{
				if (GetNode(i).LTG(j) <= nb)
				{
					//UO() << "Node" << i << "[" << j << "]=" << GetNode(i).LTG(j) << " --> " << fixedresort(GetNode(i).LTG(j)) << "\n";
					GetNode(i).LTG(j) = fixedresort(GetNode(i).LTG(j));
				} 
			}
		}
		//UO() << "Fixedresort=" << fixedresort << "\n";

		//UO() << "Phi_CB=" << Phi_CB << "\n";
		//UO() << "evec=" << evec << "\n";


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Kr and Mr:

		GetMBS()->UO() << "build Mr...\n";

		double ts = -GetClockTime();
		//Matrix Phi_CBtp = Phi_CB;//.GetTp();
		//Phi_CBtp.TpYs();

		//Mr = Phi_CBtp*(Msparse*Phi_CB);
		Matrix mmid;
		Mult(Msparse,Phi_CB,mmid);


		if (issym)
			MultSymTp(Phi_CB,mmid,Mr);
		else
			MultTp(Phi_CB,mmid,Mr);

		GetMBS()->UO() << "build Kr...\n";
		Mult(Ksparse,Phi_CB,mmid);

		if (issym)
			MultSymTp(Phi_CB,mmid,Kr);
		else
			MultTp(Phi_CB,mmid,Kr);
		//Kr = Phi_CBtp*(Ksparse*Phi_CB);

		ts += GetClockTime();

		GetMBS()->UO() << "done in" << ts << "seconds\n";

		UO() << "Mr=" << Mr << "\n";
		UO() << "Kr=" << Kr << "\n";

		Matrix test = Mr;
		int rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Mass-matrix not invertable!!!\n";}

		test = Kr;
		rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Stiffness-matrix not invertable!!!\n";}

	}
	else
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Direct computation with full matrices
		M.SetSize(sosfull,sosfull);
		M.FillWithZeros();
		K.SetSize(sosfull,sosfull);
		K.FillWithZeros();
		Matrix tempm;

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			e.StiffnessMatrix(tempm);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= sos; k++)
				{
// (AD) changed () to .Get()
					K(ltg.Get(j),ltg.Get(k)) += tempm(j,k);
//					K(ltg(j),ltg(k)) += tempm(j,k);
				}
			}
			e.EvalMff(tempm,0);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= sos; k++)
				{
// (AD) changed () to .Get()
					M(ltg.Get(j),ltg.Get(k)) += tempm(j,k);
//					M(ltg(j),ltg(k)) += tempm(j,k);
				}
			}
		}

		UO() << "K is symmetric=" << K.IsSymmetric(1e-8) << "\n";
		UO() << "M is symmetric=" << M.IsSymmetric(1e-15) << "\n";
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//extract K_II and M_II matrix

		int nb = NBModes();
		int nf = sosfull-nb; //number of internal (flexible) dof - unreduced
		int nr = NIModes();


		Matrix K_II, M_II, K_IB;
		Matrix X; //K_II^(-1)*K_IB

		K_II.SetSize(nf, nf);
		K_IB.SetSize(nf, nb);
		M_II.SetSize(nf, nf);

		K_II.SetAll(0);
		M_II.SetAll(0);
		K_IB.SetAll(0);

		K_II.AddSubmatrix(K,nb+1,nb+1,1,1,nf,nf,1.);
		M_II.AddSubmatrix(M,nb+1,nb+1,1,1,nf,nf,1.);
		K_IB.AddSubmatrix(K,nb+1,1,1,1,nf,nb,1.);

		//UO() << "K_II=" << K_II << "\n";
		//UO() << "M_II=" << M_II << "\n";

		Matrix M_IIinv = M_II;
		int rv = M_IIinv.Invert2();
		if (!rv) {UO() << "ERROR: CMS-Mass matrix not invertable!!!\n";}

		Matrix EVmat = M_IIinv*K_II;

		Eigenvalue EV(EVmat);
		Vector eval = EV.GetRealEigenvalues();
		Matrix evec = EV.GetV();

		//UO() << "eigenvalues = " << eval << "\n";
		//UO() << "eigenvectors = " << evec << "\n";

		//sort eigenvalues to get smallest nr eigenmodes:
		TArray<double> eval2;
		eval2.SetLen(eval.Length());
		for (int i=1; i <= eval.Length(); i++) {eval2(i) = eval(i);}

		TArray<int> eref;
		eref.SetLen(eval.Length());
		for (int i=1; i <= eval.Length(); i++) {eref(i) = i;}

		QuicksortDouble(eval2,eref);

		for (int i = 1; i <= nr; i++)
		{
			eigval(i) = eval2(i);
		}
		for(int i=1; i <= nf; i++)
		{
			for(int j=1; j <= nr; j++)
			{
				Phi_r(i, j) = evec(i,eref(j));
			}
		}
		double sign = 1;
		int foundsign = 0;
		for(int j=1; j <= nr; j++)
		{
			double norm = 0;
			for(int i=1; i <= nf; i++)
			{
				if (!foundsign && Phi_r(i,j) != 0) 
				{
					sign = Sgn(Phi_r(i,j));
					foundsign = 1;
				}
				norm += Sqr(Phi_r(i, j));
			}
			if (norm == 0) norm = 1;
			norm = sign*1./sqrt(norm);
			for(int i=1; i <= nf; i++)
			{
				Phi_r(i, j) *= norm;
			}
			foundsign = 0;
			sign = 1;
		}

		X = K_II;
		X.Invert2();
		X = X*K_IB;
		X *= -1;


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Phi_CB matrix:
		Phi_CB.SetSize(sosfull, NModes());
		Phi_CB.SetAll(0);

		for(int i=1; i <= nb; i++)
		{
			Phi_CB(i, i) = 1;
		}

		for(int i=1; i <= nf; i++)
		{
			for(int j=1; j <= nb; j++)
			{
				Phi_CB(nb+i, j) = X(i,j);
			}
		}

		for(int i=1; i <= nf; i++)
		{
			for(int j=1; j <= nr; j++)
			{
				Phi_CB(nb+i, nb+j) = Phi_r(i,j);
			}
		}

		//UO() << "Phi_CB=" << Phi_CB << "\n";
		//UO() << "evec=" << evec << "\n";


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Kr and Mr:
		Mr = Phi_CB.GetTp()*M*Phi_CB;
		Kr = Phi_CB.GetTp()*K*Phi_CB;


		Matrix test = Mr;
		int rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Mass-matrix not invertable!!!\n";}

		test = Kr;
		rv2 = test.Invert2();
		if (!rv2) {UO() << "ERROR: reduced CMS-Stiffness-matrix not invertable!!!\n";}


	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set new initial vector!!!
	Vector xicopy = x_init; //length = 6
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
	for (int i = 1; i <= initlen; i++)
	{
		x_init(i+2*NModes()+initlen) = xicopy(i+initlen);
	}
	//UO() << "x_init=" << x_init << "\n";

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//intialize a few other things:
	double totalmass = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		e.Initialize();
		totalmass += e.GetMass();
	}
	mass = totalmass;

	double totalvolume = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		e.Initialize();
		//totalvolume += e.GetMass()/e.GetRho();
		double v1 = e.GetVolume();
		totalvolume += e.GetVolume();

		//???rho of ffrf elements is zero!!!
	}
	volume = totalvolume;

}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CMSElement2D::GetI1(Vector& I1) //Shabana p. 209-211
{
	int n = Dim(); //maybe not dimension for 3D...
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
}

double CMSElement2D::GetIkl(int k, int l)
{
	double v = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		v += e.GetIkl(k,l);
	}
	return v;
}

void CMSElement2D::GetIbarkl(int k, int l, Vector& I1)
{
	I1.SetLen(SOSFull());
	I1.SetAll(0);
	Vector tempv;

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();

		const	TArray<int>& ltg = e.GetLTGArray();

		e.GetIbarkl(k,l,tempv);
		for (int j=1; j <= sos; j++)
		{
// (AD) changed () to .Get()
			I1(ltg.Get(j)) += tempv(j);
//			I1(ltg(j)) += tempv(j);
		}
	}
	I1 = Phi_CB.GetTp()*I1;
}

void CMSElement2D::GetSbar(Matrix& Sbar)
{
	Sbar.SetSize(Dim(), SOSFull());
	Sbar.FillWithZeros();
	Matrix tempm;

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		e.GetSbar(tempm);

		for (int j=1; j <= Dim(); j++)
		{
			for (int k=1; k <= sos; k++)
			{
// (AD) changed () to .Get()
				Sbar(j,ltg.Get(k)) += tempm(j,k);
//				Sbar(j,ltg(k)) += tempm(j,k);
			}
		}
	}

	Sbar = Sbar*Phi_CB;
}

void CMSElement2D::GetSbarkl(int k, int l, Matrix& Sbar)
{
	int sosfull = SOSFull();
	SparseMatrix SbarSM(sosfull, sosfull, CMSmin_sparse_size);
	SbarSM.FillWithZeros();
	Matrix tempm;

	//int issym = 1;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		int sos = e.FlexDOF();
		const	TArray<int>& ltg = e.GetLTGArray();

		e.GetSbarkl(k,l,tempm); //this matrix is not symmetric!!!!
		//if (!tempm.IsSymmetric(1e-8)) issym = 0;
		SbarSM.AddMatrix(ltg,ltg,sos,sos,tempm);
	}
	Sbar = Phi_CB.GetTp()*(SbarSM*Phi_CB);
}

/* void CMSElement2D::GetSbarkl(int k, int l, Matrix& Sbar)
{
Sbar.SetSize(SOSFull(), SOSFull());
Sbar.FillWithZeros();
Matrix tempm;

for (int i=1; i<=NFFRFElements(); i++) 
{
Element& e = GetFFRFElement(i);
int sos = e.FlexDOF();
const	TArray<int>& ltg = e.GetLTGArray();

e.GetSbarkl(k,l,tempm);
for (int j=1; j <= sos; j++)
{
for (int k=1; k <= sos; k++)
{
Sbar(ltg(j),ltg(k)) += tempm(j,k);
}
}
}
Sbar = Phi_CB.GetTp()*(Sbar*Phi_CB);
}*/

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int noFFRF = 0;
void CMSElement2D::EvalM(Matrix& m, double t) 
{
	TMStartTimer(23);
	//to be implemented
	//UO() << "RF EvalM\n";
	if (noFFRF)
	{
		m.SetAll(0);
		m.AddSubmatrix(Mr,1,1,1,1,NModes(),NModes(),1);;
		m(NModes()+1,NModes()+1) = 1;
		m(NModes()+2,NModes()+2) = 1;
		m(NModes()+3,NModes()+3) = 1;
	}
	else
	{
		static Vector temp;
		static Vector temp2;
		//static Matrix mtemp;

		int off = NModes(); //offset where rigid body entries start!
		xg.SetLen(NModes());
		for (int i=1; i <= NModes(); i++) xg(i) = XG(i);

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//mRR = Unit(2,2) * mass
		m(1+off,1+off) = GetMass();
		m(2+off,2+off) = GetMass();
		//UO() << "Mass=" << GetMass() << "\n";

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//mRtheta = A_theta*[I1+Sbar*qf]
		Vector2D I1;
		I1(1) = I1S(1);
		I1(2) = I1S(2);

		Vector2D pt(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NModes(); j++)
			{
				pt(i) += SbarS(i,j)*xg(j); //this term is important for jacobian
			}
		}

		I1 += pt;

		Matrix3D ADphi = GetRotMatrixDphi2D();
		I1 = ADphi*I1;
		m(1+off,3+off) = I1(1);
		m(2+off,3+off) = I1(2);
		m(3+off,1+off) = I1(1);
		m(3+off,2+off) = I1(2);


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//m_Rf: = A*Sbar //very important term

		Matrix3D A = GetRotMatrix2D();

		for (int j = 1; j <= NModes(); j++)
		{
			Vector2D v(SbarS(1,j),SbarS(2,j)); //stimmt mit rho*GetH() überein!!!
			v = A*v;
			m(j,1+off) = v(1);
			m(j,2+off) = v(2);
			m(1+off,j) = v(1);
			m(2+off,j) = v(2);
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//m_theta_theta_rr: I_11+I_22
		m(3+off,3+off) = I1122S; //(GetIkl(1,1)+GetIkl(2,2));

		//m_theta_theta_rf: 2*(Ibar_11+Ibar_22)*q_f
		double I11, I22;
		I11 = 0;
		for (int i=1; i <= NModes(); i++)
		{
			I11 += Ibar11S(i)*XG(i);
		}
		I22 = 0;
		for (int i=1; i <= NModes(); i++)
		{
			I22 += Ibar22S(i)*xg(i);
		}
		m(3+off,3+off) += 2.*(I11+I22);

		//m_theta_theta_ff: q_f^T*m_ff*q_f
		double mthth = 0;

		Mult(Mr,xg,temp);
		mthth = xg*temp;

		m(3+off,3+off) += mthth;

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//m_theta_f:
		temp.SetLen(NModes());
		temp = Ibar12S;
		temp -= Ibar21S; //-I_21+I_12

		for (int i=1; i <= NModes(); i++) // qf^T*(S12-S21)
		{
			for (int j=1; j <= NModes(); j++)
			{
				temp(i) += xg(j)*Sbar_tilde(j,i); //this term is important for jacobian!!!
			}
		}

		for (int i=1; i <= NModes(); i++) 
		{
			m(i,off+3) = temp(i);
			m(off+3,i) = temp(i);
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//m_ff
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

void CMSElement2D::EvalF2(Vector& f, double t) 
{
	Body2D::EvalF2(f,t);
	TMStartTimer(22);

	//TMStartTimer(23);
	//UO() << "RF EvalF2\n";
	//to be implemented
	xg.SetLen(NModes());
	for (int i=1; i <= NModes(); i++) xg(i) = XG(i);

	static Vector temp; 

	if (!(FastStiffnessMatrix() == 3 && GetMBS()->IsJacobianComputation()))
	{
		Mult(Kr,xg,temp);
		for (int i=1; i <= NModes(); i++) f(i) -= temp(i);
	}

	int quickjac = 0;

	//TMStopTimer(23);

	if (!noFFRF  && !GetMBS()->IsJacobianComputation())
	{

		TMStartTimer(17);
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//quadratic velocity vector, Shabana p.229:

		static Vector Qvf;
		Qvf.SetLen(NModes());
		Qvf.SetAll(0);

		double Qvtheta;
		Vector2D QvR(0.,0.);

		Matrix3D A = GetRotMatrix2D();
		Matrix3D Atheta = GetRotMatrixDphi2D();

		double theta = GetAngle2D();
		double thetap = GetAngle2DP();


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Qv_R:
		Vector2D I1;
		I1(1) = I1S(1);
		I1(2) = I1S(2);

		Vector2D pt(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NModes(); j++)
			{
				pt(i) += SbarS(i,j)*xg(j);
			}
		}

		I1 += pt;
		QvR = Sqr(thetap)*(A*I1);


		pt=Vector2D(0.,0.);
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NModes(); j++)
			{
				pt(i) += SbarS(i,j)*XGP(j);
			}
		}
		QvR += (-2.*thetap)*(Atheta*pt);


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Qv_theta: -2*theta_p* q_f_p*(m_ff*q_f+Ibar_0)
		temp.SetLen(NModes());
		temp.SetAll(0);

		double mthth = 0;

		if (!(quickjac && GetMBS()->IsJacobianComputation()))
		{
			Mult(Mr,xg,temp);
		}

		temp += Ibar11S;
		temp += Ibar22S; // Ibar_0 = Ibar_11+Ibar_22

		for (int i=1; i <= NModes(); i++)
		{
			mthth += XGP(i)*temp(i);
		}
		Qvtheta = -2.*thetap * mthth; //does not influence up very much



		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Qv_f: //temp is still m_ff*q_f+Ibar_0



		for (int i=1; i <= NModes(); i++)
		{
			Qvf(i) = Sqr(thetap)*temp(i);
		}

		for (int i=1; i <= NModes(); i++) xg(i) = XGP(i);

		if (!(quickjac && GetMBS()->IsJacobianComputation()))
		{
			Mult(Sbar_tilde,xg,temp);
			Qvf += (2.*thetap)*temp; 	//this term becomes important for high frequencies!
		}



		//fill in velocity vector terms (positive on right hand side):
		for (int i=1; i <= NModes(); i++)
		{
			f(i) += Qvf(i);
		}

		f(NModes()+1) += QvR(1);
		f(NModes()+2) += QvR(2);
		f(NModes()+3) += Qvtheta;

		TMStopTimer(17);
	}
	TMStopTimer(22);
}; 

int CMSElement2D::FastStiffnessMatrix() const 
{
	return 3; //usually take 3
}

void CMSElement2D::StiffnessMatrix(Matrix& m) //fill in sos x sos components, m might be larger
{
	//negative stiffness matrix!!!!
	for (int i=1; i <= NModes(); i++)
	{
		for (int j=1; j <= NModes(); j++)
		{
			m(i,j) = -Kr(i,j);
		}
	}

}


void CMSElement2D::GetIntDuDq(Matrix& dudq)
{
	dudq.FillWithZeros();

	if (!(GetMBS()->IsJacobianComputation()))
	{

		int off = NModes();

		//same as mRR/rho
		dudq(1+off,1) = volume;
		dudq(2+off,2) = volume;

		double rho = mass/volume;

		//+++++++++++++++++++++++++++++++++++
		//same as mRtheta/rho = A_theta*[I1+Sbar*qf]/rho
		Vector2D I1;
		I1(1) = I1S(1);
		I1(2) = I1S(2);

		Vector2D pt(0.,0.);

		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= NModes(); j++)
			{
				pt(i) += SbarS(i,j)*XG(j);
			}
		}

		I1 += pt;
		Matrix3D ADphi = GetRotMatrixDphi2D();
		I1 = ADphi*I1;
		dudq(3+off,1) = I1(1)/rho;
		dudq(3+off,2) = I1(2)/rho;

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
}

void CMSElement2D::GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi) 
{
	//Only for rigid body motion of frame
	int off = NModes();
	dpdqi.SetSize(NModes()+3,2);
	dpdqi.FillWithZeros();

	dpdqi(NModes()+1,1) = 1;
	dpdqi(NModes()+1,2) = 0;
	dpdqi(NModes()+2,1) = 0;
	dpdqi(NModes()+2,2) = 1;
	// dR/dphi*ploc:		
	double phi = GetAngle2D();
	double sphi = sin(phi);
	double cphi = cos(phi);
	dpdqi(NModes()+3,1) = -ploc.X()*sphi-ploc.Y()*cphi;
	dpdqi(NModes()+3,2) =  ploc.X()*cphi-ploc.Y()*sphi;

};

void CMSElement2D::GetNodedPosdqT(int node, Matrix& dpdqi)
{
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

}

void CMSElement2D::AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f)   // f += dpdq*lambda
{
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
}


const double& CMSElement2D::GetXactFull(int i) const 
{
	//this function usually should not be called, only for evaluation with WriteSol()!!!
	//and with constraints!!!

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

double& CMSElement2D::GetXactFull(int i)
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

const double& CMSElement2D::GetDrawValueFull(int i) const 
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
			if (j == GetMBS()->GetDOption(100)+0*NBModes()) v(ii)+=0.1*GetMBS()->GetDOption(105)*Phi_CB(ii,j);
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


















//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ACRSElement2D::Initialize() 
{
	ReferenceFrame2D::Initialize();
}


void ACRSElement2D::DoModalAnalysis()
{
	UO() << "+++++++++++++++++\nDo Modal Analysis\n+++++++++++++++++\n";
	UO() << "n-FFRF-elements=" << NFFRFElements() << "\n";

	//int sos1 = 1; //not used any more
	//int sos2 = sosfull+1;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//size of stiffness and mass-matrix: sosfull x sosfull
	int nb = NBModes();
	int nf = sosfull-nb; //number of internal (flexible) dof - unreduced
	int nr = NIModes();

	UO() << "sos_full=" << sosfull << "\n";
	UO() << "nb=" << nb << ", nf=" << nf << ", nr=" << nr << "\n";
 
	Matrix X; //K_II^(-1)*K_IB
	Matrix Phi_r(nf,nr); // = Phi_r
	Vector eigval(nr/2); 
	TArray<int> eref;
	Matrix evec(nf,nr/2); // = Phi_r

	Matrix Phi_CBsave, Ksave, Msave;
	Matrix M;
	Matrix K;

	if (UseSparseMK() || usematlab)
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//use matlab to compute eigenvectors!
		//extract K_II and M_II matrix

		int issym = 1; //faster operations ....

		SparseMatrix K_II, M_II;
		Matrix K_IB;

		UO() << "writing to matlab ... \n";

		K_II.CopyFrom(Ksparse,nb+1,nb+1,nb+nf,nb+nf);
		M_II.CopyFrom(Msparse,nb+1,nb+1,nb+nf,nb+nf);
		K_IB.CopyFrom(Ksparse,nb+1,1,nb+nf,nb);


		if (!GetMBS()->GetSolSet().eigsolv_reuse_last_eigvec) //reuse eigenvectors
		{
			{
				ofstream wait(path+mystr("\\wait.txt"));
				wait << "0\n";

				ofstream matout(path+mystr("\\matout.m"));
				ofstream M_IIdat(path+mystr("\\M_IId.dat"));
				ofstream K_IIdat(path+mystr("\\K_IId.dat"));
				ofstream K_IBdat(path+mystr("\\K_IBd.dat"));


				M_II.PrintToMatlabSparse(M_IIdat);
				K_II.PrintToMatlabSparse(K_IIdat);
				K_IB.PrintToMatlabSparse(K_IBdat);

				matout << "cd('" << path << "')\n";
				matout << "load M_IId.dat\n";
				matout << "M_II=spconvert(M_IId);\n";
				matout << ";\n";

				matout << "load K_IId.dat\n";
				matout << "K_II=spconvert(K_IId);\n";
				matout << ";\n";

				matout << "load K_IBd.dat\n";
				matout << "K_IB=spconvert(K_IBd);\n";
				matout << ";\n";

				matout << "neig=" << nr/2 << "\n"; 
				matout << "comp_eig_CMS(K_II,M_II,K_IB," << nr/2 << ");\n";

				matout << "exit\n"; 
				matout << "exit\n"; 
			}

			UO() << "compute eigenvectors in matlab ... \n";
			ShellExecute(NULL, "open", "matlab",	
				" /r matout", 
				//" -nodesktop /r matout", 
				//" -nosplash /r matout", 
				//path.c_str(), SW_HIDE);
				path.c_str(), SW_SHOW);


			int goon = 0;
			while (!goon)
			{
				{
					Sleep(1000);
					ifstream wait(path+mystr("\\wait.txt"));
					int test;
					wait >> test;
					if (test == 1) goon = 1;
					//UO() << "not ready\n";
				}
			}
			UO() << "Eigenvectors computed in MATLAB!!!\n";
		}
		else
		{
			UO() << "*****************\nEigenvectors reused from last computation!\n*****************\n";
		}

		//GetMBS()->InstantMessageText("Compute eigenvectors with MATLAB\n AFTERWARDS press OK!");


		ifstream matin(path+mystr("\\matin.txt"));

		UO() << "read eigenvectors from matlab ...\n";
		for (int i=1; i <= nr/2; i++)
		{
			matin >> eigval(i);
		}


		int cnt = 0;
		for(int i=1; i <= nf; i++)
		{
			if (cnt > 20000 && nr/2*nf >= 50000) 
			{
				UO() << "\rread EV:" << (100*i)/nf << "%";
				cnt = 0;
			}
			cnt += nr/2;
			for(int j=1; j <= nr/2; j++)
			{
				matin >> evec(i, j);
			}
		}
		UO() << "\r";

		//to be compatible with non-sparse version:
		eref.SetLen(nr/2);
		for (int i=1; i <= nr/2; i++) {eref(i) = i;}

		/*
		//normalize eigenvectors, such that they are the same as internally computed eigenvectors:
		double sign = 1;
		int foundsign = 0;
		for(int j=1; j <= nr/2; j++)
		{
			double norm = 0;
			for(int i=1; i <= nf; i++)
			{
				if (!foundsign && evec(i,j) != 0) 
				{
					sign = Sgn(evec(i,j));
					foundsign = 1;
				}
				norm += Sqr(Phi_r(i, j));
			}
			if (norm == 0) norm = 1;
			norm = sign*1./sqrt(norm);
			for(int i=1; i <= nf; i++)
			{
				evec(i, j) *= norm;
			}
			foundsign = 0;
			sign = 1;
		}

		//changed for ACRS:
		for(int i=1; i <= nf/2; i++)
		{
			for(int j=1; j <= nr/2; j++) 
			{
				Phi_r(2*i-1, 2*j-1) = evec(2*i-1,eref(j));
				Phi_r(2*i  , 2*j-1) = evec(2*i  ,eref(j));

				Phi_r(2*i-1, 2*j  ) =-evec(2*i  ,eref(j));
				Phi_r(2*i  , 2*j  ) = evec(2*i-1,eref(j));
				
			}
		}*/


		X.SetSize(K_IB.Getrows(),K_IB.Getcols());
		int krows = X.Getrows();
		int kcols = X.Getcols();

		cnt = 0;
		for(int i=1; i <= krows; i++)
		{
			if (cnt > 50000 && krows*kcols >= 200000) 
			{
				UO() << "\rread X_IB:" << (100*i)/krows << "%";
				cnt = 0;
			}
			cnt += kcols;
			for(int j=1; j <= kcols; j++)
			{
				matin >> X(i, j);
			}
		}
		UO() << "\r";

		//UO() << "eigval1=" << eigval << "\n";
		//UO() << "Phi_r1=" << Phi_r << "\n";
		UO() << "Eigenvalues(Hz)=";
		for (int i=1; i <= eigval.Length(); i++)
		{
			UO() << sqrt(eigval(i))/(2.*MY_PI) << "Hz, ";
		}
		UO() << "\n";


	}
	else
	{
		M = Mr;
		K = Kr;

		UO() << "K is symmetric=" << K.IsSymmetric(1e-6) << "\n";
		UO() << "M is symmetric=" << M.IsSymmetric(1e-10) << "\n";

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//extract K_II and M_II matrix

		int nb = NBModes();
		int nf = sosfull-nb; //number of internal (flexible) dof - unreduced
		int nr = NIModes();

		//changed for ACRS:
		if (((int)(nf/Dim()))*Dim() != nf) UO() << "ERROR: degrees of freedom not nodal, cannot use ACRS!!!\n";

		if (nr/2 > nf) 
		{
			GetMBS()->UO().InstantMessageText("ERROR: cannot compute more eigenvectors\nthan internal DOF!!!\n"); return;
		}

		Matrix K_II, M_II, K_IB;

		K_II.SetSize(nf, nf);
		K_IB.SetSize(nf, nb);
		M_II.SetSize(nf, nf);

		K_II.SetAll(0);
		M_II.SetAll(0);
		K_IB.SetAll(0);

		K_II.AddSubmatrix(K,nb+1,nb+1,1,1,nf,nf,1.);
		M_II.AddSubmatrix(M,nb+1,nb+1,1,1,nf,nf,1.);
		K_IB.AddSubmatrix(K,nb+1,1,1,1,nf,nb,1.);

		//UO() << "K_II=" << K_II << "\n";
		//UO() << "M_II=" << M_II << "\n";

		Matrix M_IIinv = M_II;
		int rv = M_IIinv.Invert2();
		if (!rv) {UO() << "ERROR: CMS-Mass matrix not invertable!!!\n";}

		Matrix EVmat = M_IIinv*K_II;

		Eigenvalue EV(EVmat);
		Vector eval = EV.GetRealEigenvalues();
		evec = EV.GetV();

		//UO() << "eigenvalues = " << eval << "\n";
		//UO() << "eigenvectors = " << evec << "\n";

		//sort eigenvalues to get smallest nr eigenmodes:
		TArray<double> eval2;
		eval2.SetLen(eval.Length());
		for (int i=1; i <= eval.Length(); i++) {eval2(i) = eval(i);}

		eref.SetLen(eval.Length());
		for (int i=1; i <= eval.Length(); i++) {eref(i) = i;}

		QuicksortDouble(eval2,eref);

		for (int i = 1; i <= nr/2; i++) 
		{
			eigval(i) = eval2(i);
		}

		X = K_II;
		X.Invert2();
		X = X*K_IB;
		X *= -1;

	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//build perpendicular eigenvectors:
	//changed for ACRS:
	double sign = 1;
	int foundsign = 0;
	for(int j=1; j <= nr/2; j++) 
	{
		double norm = 0;
		for(int i=1; i <= nf; i++)
		{
			if (!foundsign && evec(i,eref(j)) != 0) 
			{
				sign = Sgn(evec(i,eref(j)));
				foundsign = 1;
			}
			norm += Sqr(evec(i, eref(j)));
		}
		if (norm == 0) norm = 1;
		norm = sign*1./sqrt(norm);
		for(int i=1; i <= nf; i++)
		{
			evec(i, eref(j)) *= norm;
		}
		foundsign = 0;
		sign = 1;
	}

	//changed for ACRS:
	for(int i=1; i <= nf/2; i++)
	{
		for(int j=1; j <= nr/2; j++) 
		{
			Phi_r(2*i-1, 2*j-1) = evec(2*i-1,eref(j));
			Phi_r(2*i  , 2*j-1) = evec(2*i  ,eref(j));

			Phi_r(2*i-1, 2*j  ) =-evec(2*i  ,eref(j));
			Phi_r(2*i  , 2*j  ) = evec(2*i-1,eref(j));

		}
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//build Phi_CB matrix:

	int rigidmodes = 1;
	nbmodes += (nb-reduced1*2-reduced2*2)*rigidmodes; //2 rigid body modes ...

	Phi_CB.SetSize(sosfull, NModes());
	Phi_CB.SetAll(0);

	for(int i=1; i <= nb; i++)
	{
		Phi_CB(i, i) = 1;
	}

	for(int i=1; i <= nf; i++)
	{
		for(int j=1; j <= nb; j++)
		{
			Phi_CB(nb+i, j) = X(i,j);
		}
	}

	for(int i=1; i <= nf; i++)
	{
		for(int j=1; j <= nr; j++)
		{
			Phi_CB(nb+i, nb+j) = Phi_r(i,j);
		}
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//modify modes:
	if (1)
	{
		//add correct constant mode for first reference mode ?...
		//add correct linear mode for second reference mode ...
		int ref1x = GetNode(RefNode1()).LTG(1);
		int ref1y = GetNode(RefNode1()).LTG(2);
		int ref2x = GetNode(RefNode2()).LTG(1);
		int ref2y = GetNode(RefNode2()).LTG(2);
		double ref1px = GetNode(RefNode1()).Pos().X();
		double ref1py = GetNode(RefNode1()).Pos().Y();
		double ref2px = GetNode(RefNode2()).Pos().X();
		double ref2py = GetNode(RefNode2()).Pos().Y();

		double lnormx = fabs(GetNode(RefNode2()).Pos().X()-GetNode(RefNode1()).Pos().X());
		if (lnormx == 0.) 
			lnormx = 1;
		//GetMBS()->InstantMessageText("Reference nodes must have different x-coordinates!");
		double lnormy = fabs(GetNode(RefNode2()).Pos().Y()-GetNode(RefNode1()).Pos().Y());
		if (lnormy == 0.) 
			lnormy = 1;


		//UO() << "Phi_CB=" << Phi_CB << "\n";

		//change deformation modes of RefNode2() into 2 modes to model rotation: [X Y] and [-Y X]
		for(int i=1; i <= NNodes(); i++)
		{
			const Node& n = GetNode(i);
			/*
			//use simpler rigid modes, just need two modes for rigid body motion
			if (n.LTG(1) > 0 || n.LTG(1) == ref2x) 
			{
			Phi_CB(n.LTG(1), ref2x) = (n.Pos().X()-ref1px)/lnormx; 
			Phi_CB(n.LTG(1), ref2y) =-(n.Pos().Y()-ref1py)/lnormx; 
			}
			if (n.LTG(2) > 0 || n.LTG(2) == ref2y) 
			{
			Phi_CB(n.LTG(2), ref2x) = (n.Pos().Y()-ref1py)/lnormx; 
			Phi_CB(n.LTG(2), ref2y) = (n.Pos().X()-ref1px)/lnormx;
			}

			if (n.LTG(1) > 0 || n.LTG(1) == ref1x) 
			{
			Phi_CB(n.LTG(1), ref1x) = (-n.Pos().X()+ref2px)/lnormx; 
			Phi_CB(n.LTG(1), ref1y) =-(-n.Pos().Y()+ref2py)/lnormx; 
			}
			if (n.LTG(2) > 0 || n.LTG(2) == ref1y) 
			{
			Phi_CB(n.LTG(2), ref1x) = (-n.Pos().Y()+ref2py)/lnormx; 
			Phi_CB(n.LTG(2), ref1y) = (-n.Pos().X()+ref2px)/lnormx;
			}*/


			if (n.LTG(1) > 0 || n.LTG(1) == ref2x) 
			{
				Phi_CB(n.LTG(1), ref2x) = (n.Pos().X()-ref1px)/lnormx; 
				Phi_CB(n.LTG(1), ref2y) = 0; 
			}
			if (n.LTG(2) > 0 || n.LTG(2) == ref2y) 
			{
				Phi_CB(n.LTG(2), ref2x) = 0; 
				Phi_CB(n.LTG(2), ref2y) = (n.Pos().Y()-ref1py)/lnormy;
			}


			if (n.LTG(1) > 0 || n.LTG(1) == ref1x) 
			{
				Phi_CB(n.LTG(1), ref1x) = (-n.Pos().X()+ref2px)/(lnormx); 
				Phi_CB(n.LTG(1), ref1y) = 0; 
			}
			if (n.LTG(2) > 0 || n.LTG(2) == ref1y) 
			{
				Phi_CB(n.LTG(2), ref1x) = 0; 
				Phi_CB(n.LTG(2), ref1y) = (-n.Pos().Y()+ref2py)/(lnormy);
			}

		}

		//second correction for static modes:
		//recompute rotational deformation modes to exclude deformation in boundary nodes

		for(int j=1; j <= nb; j++)
		{
			if (j != ref1x && j != ref1y && j != ref2x && j != ref2y)
			{
				double factx = Phi_CB(j, ref2x);
				double facty = Phi_CB(j, ref2y);
				for(int i=1; i <= nf+nb; i++)
				{
					Phi_CB(i, ref2x) -=  factx*Phi_CB(i, j); 
					Phi_CB(i, ref2y) -=  facty*Phi_CB(i, j); 
				}
			}

			if (j != ref1x && j != ref1y && j != ref2x && j != ref2y)
			{
				double factx = Phi_CB(j, ref1x);
				double facty = Phi_CB(j, ref1y);
				for(int i=1; i <= nf+nb; i++)
				{
					Phi_CB(i, ref1x) -=  factx*Phi_CB(i, j); 
					Phi_CB(i, ref1y) -=  facty*Phi_CB(i, j); 
				}
			}

		}
	}

	/* //only test
	if (0)
	{
	//original shape of RefNode2 is corrected 
	int ref2x = GetNode(RefNode2()).LTG(1);
	int ref2y = GetNode(RefNode2()).LTG(2);

	int f = 0;
	for(int i=1; i <= (nf+nb*(1-f))/2; i++)
	{
	Phi_CB(2*i-1+f*nb, nb+nr+nb-1) -= Phi_CB(2*i-1+f*nb, ref2x); 
	Phi_CB(2*i  +f*nb, nb+nr+nb-1) -= Phi_CB(2*i-0+f*nb, ref2x); 
	Phi_CB(2*i-1+f*nb, nb+nr+nb  ) -=-Phi_CB(2*i-0+f*nb, ref2y); 
	Phi_CB(2*i  +f*nb, nb+nr+nb  ) -= Phi_CB(2*i-1+f*nb, ref2y); 
	}
	}*/

	//add rigid body mode:, does not work
	if (rigidmodes)
	{
		//add correct constant mode for first reference mode ?...
		//add correct linear mode for second reference mode ...
		int ref2x = GetNode(RefNode2()).LTG(1);
		int ref2y = GetNode(RefNode2()).LTG(2);
		double ref1x = GetNode(RefNode1()).Pos().X();
		double ref1y = GetNode(RefNode1()).Pos().Y();

		//change deformation modes of RefNode2() into 2 modes to model rotation: [X Y] and [-Y X]
		int ref3x = nb+nr+1;
		int ref3y = nb+nr+2;
		int ref4x = nb+nr+3;
		int ref4y = nb+nr+4;


		Matrix Phi2 = Phi_CB;

		//last node must be refnode2!!!
		IVector ndofresort(nb);
		int off = 0;
		int nodedofnum = 0;
		for(int j=1; j <= nb; j++)
		{
			int f = 0;
			ndofresort(j) = j*2-1-off;
			for(int i=1; i <= (nb*(1-f)+nf)/2; i++)
			{
				if (!((j == GetNode(RefNode1()).LTG(2) && reduced1) || (j == GetNode(RefNode2()).LTG(2) && reduced2)))
				{
					Phi_CB(2*i-1+f*nb, j*2-1-off) = Phi2(2*i-1+f*nb, j);
					Phi_CB(2*i-0+f*nb, j*2-1-off) = Phi2(2*i-0+f*nb, j);

					Phi_CB(2*i-1+f*nb, j*2  -off) =-Phi2(2*i-0+f*nb, j  );
					Phi_CB(2*i-0+f*nb, j*2  -off) = Phi2(2*i-1+f*nb, j  );
				}
			}
			if ((j == GetNode(RefNode1()).LTG(1) && reduced1) || (j == GetNode(RefNode2()).LTG(1) && reduced2)) off += 2;
		}
		for(int j=1; j <= nr; j++)
		{
			for(int i=1; i <= (nb+nf); i++)
			{
				Phi_CB(i, j+2*(nb-reduced1-reduced2)) = Phi2(i,j+nb);
			}
		}

		//resort node DOF for constraints:
		for (int i=1; i<=nodes.Length(); i++) 
		{
			//resort node DOF for constraints
			for (int j=1; j <= GetNode(i).LTGLength(); j++)
			{
				if (GetNode(i).LTG(j) <= nb)
				{
					GetNode(i).LTG(j) = ndofresort(GetNode(i).LTG(j));
				} 
			}
		}

	}

	//UO() << "Phi_CB=" << Phi_CB << "\n";



	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (UseSparseMK() || usematlab)
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Kr and Mr:

		GetMBS()->UO() << "build Mr...\n";

		double ts = -GetClockTime();

		int issym = 1;
		Matrix mmid;
		Mult(Msparse,Phi_CB,mmid);

		if (issym)
			MultSymTp(Phi_CB,mmid,Mr);
		else
			MultTp(Phi_CB,mmid,Mr);

		GetMBS()->UO() << "build Kr...\n";
		Mult(Ksparse,Phi_CB,mmid);

		if (issym)
			MultSymTp(Phi_CB,mmid,Kr);
		else
			MultTp(Phi_CB,mmid,Kr);

		ts += GetClockTime();

		GetMBS()->UO() << "done in" << ts << "seconds\n";

	}
	else
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build Kr and Mr:
		Mr = Phi_CB.GetTp()*M*Phi_CB;
		Kr = Phi_CB.GetTp()*K*Phi_CB;

	}

	//UO() << "Kr =" << Kr << "\n";
	//UO() << "Mr =" << Mr << "\n";

	/*
	ofstream otest("..\\..\\output\\test.m");
	otest << "phi=";
	Phi_CB.PrintToMatlab(otest);
	otest << ";\n";
	otest << "Kr=";
	Kr.PrintToMatlab(otest);
	otest << ";\n";
	otest << "Mr=";
	Mr.PrintToMatlab(otest);
	otest << ";\n";*/

	Matrix test = Mr;
	int rv2 = test.Invert2();
	if (!rv2) {UO() << "ERROR: reduced CMS-Mass-matrix not invertable!!!\n";}

	test = Kr;
	rv2 = test.Invert2();
	if (!rv2) {UO() << "ERROR: reduced CMS-Stiffness-matrix not invertable!!!\n";}

	//modal transformation for H-matrix:
	Hr = Phi_CB.GetTp()*Hr;
	//GetMBS()->InstantMessageText("ModalAnalysis 6"); //nr
	GetMBS()->UO() << "Modal Analysis done\n";
}




void ACRSElement2D::FinishAssembly() 
{
	//first command
	sosfull = SOSFull(); //very important for SOS() and SOSowned()!!!!!

	SetUseSparseMK(1/*GetMBS()->NLS_UseSparseSolver()*/); //muss übereinstimmen mit UseSparseSolver()
	if (TransformJacApply()) GetMBS()->SetTransformJacApply(1);

	usematlab = 1; //???compare usematlab --> error?

	int ss = GetMBS()->UseSparseSolver();

	newmode = 1*IsCMS();
	reduced1 = 1*IsCMS();
	reduced2 = 0*IsCMS();

	//UO() << "compute factorized matrices!\n";

	//GetMBS()->InstantMessageText("Finish Assembly 1");
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//give elements and nodes DOF-reference numbers

	int sos1 = 1;
	int sos2 = sosfull+1;

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

	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
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
		//set flag, that ltg list is built this time
		TMBSElement t = e.GetType();
		e.AddType(TCMSflag);
		e.LinkToElements();
		e.SetType(t);

		//UO() << "ltg=" << e.GetLTGArray() << "\n";
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//compute necessary boundary nodes (DOF):
	//only for CMS
	if (IsCMS())
	{
		UO() << "compute boundary modes\n";

		boundarydof.SetLen(sosfull);
		for (int i=1; i <= boundarydof.Length(); i++)
		{
			boundarydof(i) = 0;
		}

		for (int i=1; i <= boundarynode.Length(); i++)
		{
			boundarydof(GetNode(boundarynode(i)).LTG(1)) = 1; //fix all components of node
			boundarydof(GetNode(boundarynode(i)).LTG(2)) = 1;
		}

		nbmodes = 0;
		for (int i=1; i <= boundarydof.Length(); i++)
		{
			if (boundarydof(i) != 0) nbmodes++;
		}
		//UO() << "b-dof=" << boundarydof << "\n";

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//resort ltg according to (B) and (I) nodes:

		IVector resort(sosfull*2);
		IVector resort2(sosfull*2);

		resort.SetLen(0);
		for (int i=1; i <= sosfull; i++)
		{
			if (boundarydof(i) != 0) resort.Add(i); //add boundary nodes
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
		//UO() << "resort2=" << resort2 << "\n";

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);

			for (int j=1; j <= e.LTGlength(); j++)
			{
				e.LTG(j) = resort2(e.LTG(j));
			}
			//UO() << "ltg_sort=" << e.GetLTGArray() << "\n";
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

	if (!(UseSparseMK() || usematlab))
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build full system matrices: (in CMS they become reduced later on !!!
		Mr.SetSize(sosfull,sosfull);
		Mr.FillWithZeros();
		Kr.SetSize(sosfull,sosfull);
		Kr.FillWithZeros();
		Hr.SetSize(sosfull,Dim());
		Hr.FillWithZeros();
		Matrix tempm;

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			e.StiffnessMatrix(tempm);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= sos; k++)
				{
// (AD) changed () to .Get()
					Kr(ltg.Get(j),ltg.Get(k)) += tempm(j,k);
//					Kr(ltg(j),ltg(k)) += tempm(j,k);
				}
			}
			e.EvalMff(tempm,0);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= sos; k++)
				{
// (AD) changed () to .Get()
					Mr(ltg.Get(j),ltg.Get(k)) += tempm(j,k);
//					Mr(ltg(j),ltg(k)) += tempm(j,k);
				}
			}
			e.GetH(tempm);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= Dim(); k++)
				{
// (AD) changed () to .Get()
					Hr(ltg.Get(j),k) += tempm(j,k);
//					Hr(ltg(j),k) += tempm(j,k);
				}
			}

		//++++++++++++++++++++++++++++++++++++++++++++++++++
		//compute reduced H matrix:
		}
	}
	else
	{
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//build system matrices:
		usesparsepre = 1;
		if (usesparsepre)
		{
			Msparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
			Msparse.FillWithZeros();
			Ksparse.SetSize(sosfull,sosfull,CMSmin_sparse_size);
			Ksparse.FillWithZeros();
		}
		Hr.SetSize(sosfull,Dim());
		Hr.FillWithZeros();
		Matrix tempm;

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();

			const	TArray<int>& ltg = e.GetLTGArray();

			if (usesparsepre)
			{
				e.StiffnessMatrix(tempm);
				if (IsCMS()) tempm.MakeSymmetric();
				Ksparse.AddMatrix(ltg,ltg,sos,sos,tempm);

				e.EvalMff(tempm,0);
				if (IsCMS()) tempm.MakeSymmetric();
				Msparse.AddMatrix(ltg,ltg,sos,sos,tempm);
			}
			e.GetH(tempm);
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= Dim(); k++)
				{
// (AD) changed () to .Get()
					Hr(ltg.Get(j),k) += tempm(j,k);
//					Hr(ltg(j),k) += tempm(j,k);
				}
			}
		}
		//UO() << "Ksparse-alloc=" << Ksparse.GetLAlloc()*12./1.e6 << "MB\n";
		//UO() << "Msparse-alloc=" << Msparse.GetLAlloc()*12./1.e6 << "MB\n";
	}



	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set new initial vector!!!
	if (IsCMS())
	{
		DoModalAnalysis();
	}

	/*
	if (UseSparseMK())
	{
		Msparse.Destroy();
		Ksparse.Destroy();

		Msparse.CopyFrom(Mr);
		Ksparse.CopyFrom(Kr);
	}*/
	

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//set new initial vector!!!
	x_init.SetLen(2*SOS()); //sosfull or NModes

	//position:
	for (int i = 1; i <= SOS(); i++)
	{
		x_init(i) = 0;
	}

	//velocity:
	for (int i = 1; i <= SOS(); i++)
	{
		x_init(i+SOS()) = 0;
	}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//intialize a few other things:
	double totalmass = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		e.Initialize();
		totalmass += e.GetMass();
	}
	mass = totalmass;

	double totalvolume = 0;
	for (int i=1; i<=NFFRFElements(); i++) 
	{
		Element& e = GetFFRFElement(i);
		e.Initialize();
		//totalvolume += e.GetMass()/e.GetRho();
		totalvolume += e.GetVolume();
	}
	volume = totalvolume;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//initial rotation:

	Vector2D p1 = GetNodePosInit2D(RefNode1());
	Vector2D p2 = GetNodePosInit2D(RefNode2());
	double Lact = (p1-p2).Norm();
	double s = 1./Lact*(p2.Y()-p1.Y());
	double c = 1./Lact*(p2.X()-p1.X());

	initrot(1,1) = c;
	initrot(1,2) =-s;
	initrot(2,1) = s;
	initrot(2,2) = c;
	initrot.TpYs();
	initphi = atan2(s,c);

	GetMBS()->UO() << "Assembly done\n";
};

void ACRSElement2D::AddMSparse(SparseMatrix& m, double t)  //add sparse matrix into full system matrix
{
	if (t > 1.5*GetMBS()->GetStepSize() && TransformJacApply()) return;

	if (!IsCMS())
	{
		static IVector ltg2;
		static Matrix tempm;
		for (int i=1; i<=NFFRFElements(); i++) 
		{

			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();
			ltg2.SetLen(sos);

			const	TArray<int>& ltg = e.GetLTGArray();
			for (int j=1; j <= sos; j++)
// (AD) changed () to .Get()
				ltg2(j) = LTG(ltg.Get(j));
//				ltg2(j) = LTG(ltg(j));

			e.EvalMff(tempm,0);

			m.AddMatrix(ltg2,ltg2,sos,sos,tempm);
		}
	}
	else
	{
		m.AddMatrix(ltg,ltg,NModes(),NModes(),Mr);
	}
	
}

void ACRSElement2D::AddKSparse(SparseMatrix& m, double t)  //add sparse matrix into full system matrix
{
	if (t > 1.5*GetMBS()->GetStepSize() && TransformJacApply()) return;

	Matrix3D A = -1.*GetRotMatrix2D();
	Matrix3D AT = GetRotMatrix2D().GetTp();
	Matrix3D m22;
	m22.SetSize(2,2);
	static Matrix tempm;

	if (!IsCMS())
	{
		static IVector ltg2;

		//UO() << "A=" << A << "\n";

		for (int i=1; i<=NFFRFElements(); i++) 
		{
			Element& e = GetFFRFElement(i);
			int sos = e.FlexDOF();
			ltg2.SetLen(sos);
			const	TArray<int>& ltg = e.GetLTGArray();

			for (int j=1; j <= sos; j++)
// (AD) changed () to .Get()
				ltg2(j) = LTG(ltg.Get(j));
//				ltg2(j) = LTG(ltg(j));

			e.StiffnessMatrix(tempm);
			if (!TransformJacApply())
			{
				for (int i1=0; i1 < sos/2; i1++)
				{
					for (int j1=0; j1 < sos/2; j1++)
					{
						m22(1,1) = tempm(2*i1+1,2*j1+1);
						m22(1,2) = tempm(2*i1+1,2*j1+2);
						m22(2,1) = tempm(2*i1+2,2*j1+1);
						m22(2,2) = tempm(2*i1+2,2*j1+2);
						m22 = A*m22*AT;
						tempm(2*i1+1,2*j1+1) = m22(1,1);
						tempm(2*i1+1,2*j1+2) = m22(1,2);
						tempm(2*i1+2,2*j1+1) = m22(2,1);
						tempm(2*i1+2,2*j1+2) = m22(2,2);
					}
				}
			}
			else tempm*=-1;
			m.AddMatrix(ltg2,ltg2,sos,sos,tempm);
		}

	}
	else
	{
		tempm = Kr;
		int sos = NModes();
		if (!TransformJacApply())
		{
			for (int i1=0; i1 < sos/2; i1++)
			{
				for (int j1=0; j1 < sos/2; j1++)
				{
					m22(1,1) = tempm(2*i1+1,2*j1+1);
					m22(1,2) = tempm(2*i1+1,2*j1+2);
					m22(2,1) = tempm(2*i1+2,2*j1+1);
					m22(2,2) = tempm(2*i1+2,2*j1+2);
					m22 = A*m22*AT;
					tempm(2*i1+1,2*j1+1) = m22(1,1);
					tempm(2*i1+1,2*j1+2) = m22(1,2);
					tempm(2*i1+2,2*j1+1) = m22(2,1);
					tempm(2*i1+2,2*j1+2) = m22(2,2);
				}
			}
		}
		else tempm*=-1;

		m.AddMatrix(ltg,ltg,NModes(),NModes(),tempm);

	}

}

int ACRSElement2D::ElementBandwidth() const 
{
	//to set right initial size for Jacobian matrices
	if (IsCMS()) return SOS();
	else
	{
		int bw = 1;
		for (int i = 1; i <= NFFRFElements(); i++) 
		{
			bw = Maximum(bw, GetFFRFElement(i).FlexDOF());
		}

		bw *= 4; //4 elements could join in one node

		return bw;
	}
}

const double& ACRSElement2D::GetXactFull(int i) const 
{
	if (!IsCMS())
	{
		return XG(i);
	}
	else
	{
		//this function usually should not be called, only for evaluation with WriteSol()!!!
		//and with constraints!!!

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
}

double& ACRSElement2D::GetXactFull(int i)
{
	if (!IsCMS())
	{
		return XG(i);
	}
	else
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
}

const double& ACRSElement2D::GetDrawValueFull(int i) const 
{
	if (!IsCMS())
	{
		/*
		//rotated back:
		static Vector urigid;

		Matrix3D A=GetRotMatrix2DD();
		Matrix3D AT=GetRotMatrix2DD().GetTp();
		ComputeURigid(A,uref1,urigid);

		static Vector xgd;
		xgd.SetLen(SOS());
		for (int j = 1; j <= SOS(); j++)
		{
			xgd(j) = XGD(j)-urigid(j);
		}
		ApplyRotation(AT,xgd);

		return xgd(i);
		*/
		return XGD(i);
	}
	else
	{
		//multiply xgd with Phi_CB, just i-th row

		/*
		static Vector urigid;
		
		Matrix3D A;
		A.SetSize(2,2);
		double phi = 2.*MY_PI*0.25/4.*(GetMBS()->TestCnt()-3);
		A(1,1) = cos(phi);
		A(1,2) =-sin(phi);
		A(2,1) = sin(phi);
		A(2,2) = cos(phi);

		//GetMBS()->UO() << "rot=" << GetRotMatrix2DD() << ", phi=" << GetAngle2DD() << "\n";

		//A=GetRotMatrix2DD();
		Matrix3D AT=GetRotMatrix2DD().GetTp();

		Vector2D uref1 = GetNodePos2DD(RefNode1())-GetNodePosInit2D(RefNode1())+Vector2D(0,0.1);
		ComputeURigid(A,uref1,urigid);
		int mode = ((GetMBS()->TestCnt()-3)%(NModes())) + 1+0*NBModes();
		//int mode = 1+NBModes();*/

		static Vector v;//(Phi_CB.Getrows());
		v.SetLen(Phi_CB.Getrows());
		int ii = i;
		if (i > Phi_CB.Getrows()) ii -= Phi_CB.Getrows();

		v(ii) = 0;

		static Vector xgd;
		xgd.SetLen(Phi_CB.Getcols());
		for (int j = 1; j <= Phi_CB.Getcols(); j++)
		{
			xgd(j) = XGD(j);//+0*urigid(j);
		}
		//ApplyRotation(AT,xgd);

		for (int j = 1; j <= Phi_CB.Getcols(); j++)
		{
			v(ii) += Phi_CB(ii,j)*xgd(j);//(XGD(j)+0*urigid(j));
			//if (j == mode) v(ii)+=0.1*Phi_CB(ii,j);
		}
		return v(ii);
	}
}


void ACRSElement2D::GetdPosdqT(const Vector2D& ploc, Matrix& dpdqi) 
{
	//Only for rigid body motion of frame
	assert(0);
};

void ACRSElement2D::GetNodedPosdqT(int node, Matrix& dpdqi)
{
	dpdqi.SetSize(SOS(),2);
	dpdqi.FillWithZeros();

	//**

	if (newmode && IsCMS())
	{
		if ((node != RefNode1() || !reduced1) && (node != RefNode2() || !reduced2))
		{
			dpdqi(nodes(node)->LTG(1)  ,1) = 1; //dnx/dnx
			dpdqi(nodes(node)->LTG(1)+1,2) = 1; //dny/dny
			dpdqi(nodes(node)->LTG(2)  ,2) = 1; //dny/dny
			dpdqi(nodes(node)->LTG(2)+1,1) =-1; //dnx/dnx
		}
		else
		{
			dpdqi(nodes(node)->LTG(1)  ,1) = 1; //dnx/dnx
			dpdqi(nodes(node)->LTG(1)+1,2) = 1; //dny/dny
		}
	}
	else
	{
		dpdqi(nodes(node)->LTG(1),1) = 1; //dnx/dnx
		dpdqi(nodes(node)->LTG(2),2) = 1; //dny/dny
	}

}

void ACRSElement2D::AddNodedPosdqTLambda(int node, const Vector2D& lambda, Vector& f)   // f += dpdq*lambda
{
	//**
	//UO() << "l=" << lambda << ",node=" << node << "\n";
	if (newmode && IsCMS())
	{
		if ((node != RefNode1() || !reduced1) && (node != RefNode2() || !reduced2))
		{
			f(nodes(node)->LTG(1)  ) += lambda(1);
			f(nodes(node)->LTG(1)+1) += lambda(2);
			f(nodes(node)->LTG(2)  ) += lambda(2);
			f(nodes(node)->LTG(2)+1) +=-lambda(1);
		}
		else
		{
			f(nodes(node)->LTG(1)  ) += lambda(1);
			f(nodes(node)->LTG(1)+1) += lambda(2);
		}
	}
	else
	{
		f(nodes(node)->LTG(1)) += lambda(1);
		f(nodes(node)->LTG(2)) += lambda(2);

		//UO() << "n1=" << nodes(node)->LTG(1) << ", n2=" << nodes(node)->LTG(2) << "\n";
	}
}


void ACRSElement2D::EvalM(Matrix& m, double t) 
{
	TMStartTimer(23);

	//if (t < 1.*GetMBS()->GetStepSize()) GetMBS()->UO() << "***************\nERROR: Buggy Mass matrix used!!!!!!!!!!\n***************\n";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//
	int sos = SOS();
	for (int i=1; i <= sos; i++)
	{
		for (int j=1; j <= sos; j++)
		{
			m(i,j) = Mr(i,j);
		}
	}

	TMStopTimer(23);
}; 

void ACRSElement2D::EvalF2(Vector& f, double t) 
{
	Body2D::EvalF2(f,t);
	//UO() << "f2=" << f << "\n";
	
	if (GetMBS()->IsJacobianComputation()) return;

	//if (!UseSparseMK() && GetMBS()->IsJacobianComputation()) return;

	TMStartTimer(22);

	//TMStartTimer(23);
	//UO() << "RF EvalF2\n";
	//to be implemented
	xg.SetLen(SOS());
	for (int i=1; i <= SOS(); i++) xg(i) = XG(i);

	static Vector temp; 
	static Vector qF;
	static Vector ATqF;
	static Vector urigid;
	static Vector KATqF;
	static Vector dAijdq;

	Matrix3D A=GetRotMatrix2D();
	Matrix3D AT=GetRotMatrix2D().GetTp();


	//+++++++++++++++++++++++++++
	//add -A K A^T (q-q_R)

	//apply rigid body motion given by A to urigid
	Vector2D uref1 = GetNodeDisp2D(RefNode1());
	ComputeURigid(A,uref1,urigid);
	qF = xg;
	qF -= urigid;
	ATqF = qF;

	//static double lastprint = 0;
	//if (!GetMBS()->IsJacobianComputation() && t>lastprint) {UO() << "urigid=" << urigid << "\n";lastprint = t;}

	ApplyRotation(AT,ATqF);

	if (UseSparseMK() && !IsCMS())
	{
		if (usesparsepre)
		{
			Mult(Ksparse,ATqF,KATqF);
		}
		else
		{
			//if (IsCMS()) GetMBS()->InstantMessageText("ERROR: EvalF2 must have set usesparsepre!!!");

			KATqF.SetLen(ATqF.Length());
			KATqF.FillWithZeros();

			static Vector v;
			static Vector v2;
			static Matrix tempm;
			for (int i=1; i<=NFFRFElements(); i++) 
			{
				Element& e = GetFFRFElement(i);
				int sos = e.FlexDOF();

				if (sos)
				{
					const	TArray<int>& ltg = e.GetLTGArray();
					v.SetLen(sos);
// (AD) changed () to .Get()
					for (int j=1; j <= sos; j++) v(j) = ATqF(ltg.Get(j));
//					for (int j=1; j <= sos; j++) v(j) = ATqF(ltg(j));

					e.StiffnessMatrix(tempm);
					Mult(tempm,v,v2);

// (AD) changed () to .Get()
					for (int j=1; j <= sos; j++) KATqF(ltg.Get(j))+=v2(j);
//					for (int j=1; j <= sos; j++) KATqF(ltg(j))+=v2(j);
				}
			}
		}
	}
	else
	{
		//use this for CMS :
		Mult(Kr,ATqF,KATqF);
	}
	
	temp = KATqF;

	ApplyRotation(A,temp);

	for (int i=1; i <= SOS(); i++) f(i) -= temp(i);

	//TMStartTimer(17);
	//+++++++++++++++++++++++++++
	//compute f_NL:
	if (1)
	{
		Matrix3D Aij;
		Aij.SetSize(2,2);
		double DfnlDAij;
		for (int i = 1; i <= Dim(); i++)
		{
			for (int j = 1; j <= Dim(); j++)
			{
				Aij.SetAll(0);
				Aij(i,j) = 1;

				temp = KATqF;
				ApplyRotation(Aij,temp);
				DfnlDAij = qF*temp;

				ComputeDAijDq(i, j, dAijdq);
				temp = DfnlDAij * dAijdq;

				for (int i=1; i <= SOS(); i++) f(i) -= temp(i);
			}
		}
	}

	//TMStopTimer(17);
	TMStopTimer(22);
}; 


int ACRSElement2D::FastStiffnessMatrix() const 
{
	//only interesting if not sparse!!!
	return 2;
}

int ACRSElement2D::TransformJacApply() const
{
  return 1;
}

//compute Av=A^T*v in mode==0 and Av=A*v in mode==1
void ACRSElement2D::ApplyTransform(const Vector& v, Vector& Av, int mode)
{
	Matrix3D A;
	if (mode == 0) A = GetRotMatrix2D().GetTp();
	else A = GetRotMatrix2D();

	Vector2D r;
	int indx, indy;
	//UO() << "lv=" << v.Length() << ", Avl=" << Av.Length() << "\n";
	for (int i=0; i < SOS()/2; i++)
	{
		indx = LTG(2*i + 1);
		indy = LTG(2*i + 2);
		//UO() << "ix=" << indx << ", iy=" << indy << "\n";
		r.X() = v(indx);
		r.Y() = v(indy);
		r = A*r;
		Av(indx) = r.X();
		Av(indy) = r.Y();
	}

}

void ACRSElement2D::StiffnessMatrix(Matrix& m) //fill in sos x sos components, m might be larger
{
	//only used for non-sparse version!
	//negative stiffness matrix!!!!
	Matrix3D A = -1*GetRotMatrix2D();
	Matrix3D AT = GetRotMatrix2D().GetTp();
	Matrix3D m22;
	m22.SetSize(2,2);

	for (int i=0; i < SOS()/2; i++)
	{
		for (int j=0; j < SOS()/2; j++)
		{
			m22(1,1) = Kr(2*i+1,2*j+1);
			m22(1,2) = Kr(2*i+1,2*j+2);
			m22(2,1) = Kr(2*i+2,2*j+1);
			m22(2,2) = Kr(2*i+2,2*j+2);
			m22 = A*m22*AT;
			m(2*i+1,2*j+1) = m22(1,1);
			m(2*i+1,2*j+2) = m22(1,2);
			m(2*i+2,2*j+1) = m22(2,1);
			m(2*i+2,2*j+2) = m22(2,2);
		}
	}

}


//apply rotation to vector u
void ACRSElement2D::ApplyRotation(const Matrix3D& rot, Vector& u) const
{
	Vector2D r;
	for (int i=0; i < SOS()/2; i++)
	{
		r.X() = u(2*i+1);
		r.Y() = u(2*i+2);
		r = rot*r;
		u(2*i+1) = r.X();
		u(2*i+2) = r.Y();
	}
}

void ACRSElement2D::ComputeURigid(const Matrix3D& A, const Vector2D& uref1, Vector& urigid) const
{

	Vector2D u0, ur;
	Vector2D p1 = GetNodePosInit2D(RefNode1());
	urigid.SetLen(SOS());
	Matrix3D AI=A; //A-I
	AI(1,1) -= 1;
	AI(2,2) -= 1;

	double cosphi1 = AI(1,1);
	double cosphi = A(1,1);
	double sinphi = AI(2,1);

	if (IsCMS()) 
	{
		//only nbmodes affect rigid body motion!, no initial position available
		int nn = boundarynode.Length();
		urigid.FillWithZeros();
		for (int i=1; i <= nn; i++)
		{
			//**;
			//find rigid body motion for DOFs:

			if (newmode)
			{
				if ((boundarynode(i) != RefNode1() || !reduced1) && (boundarynode(i) != RefNode2() || !reduced2))
				{
					Vector2D ur2 = uref1;
					u0 = (GetNodePosInit2D(boundarynode(i))-p1);// + 0*uref1;
					urigid(GetNode(boundarynode(i)).LTG(1)  ) = u0.X()*cosphi1 + 1.*ur2.X(); //???rearrange
					urigid(GetNode(boundarynode(i)).LTG(1)+1) = u0.X()*sinphi  + 1.*ur2.Y();
					urigid(GetNode(boundarynode(i)).LTG(2)  ) = u0.Y()*cosphi1 + 0.*ur2.Y();
					urigid(GetNode(boundarynode(i)).LTG(2)+1) = u0.Y()*sinphi  - 0.*ur2.X();


				}
				else
				{
					u0 = AI * (GetNodePosInit2D(boundarynode(i))-p1) + uref1;
					urigid(GetNode(boundarynode(i)).LTG(1)  ) = u0.X();
					urigid(GetNode(boundarynode(i)).LTG(1)+1) = u0.Y();
				}
			}
			else
			{
				u0 = AI*(GetNodePosInit2D(boundarynode(i))-p1) + uref1;
				urigid(GetNode(boundarynode(i)).LTG(1)) = u0.X();
				urigid(GetNode(boundarynode(i)).LTG(2)) = u0.Y();
			}
		}

	} 
	else
	{
		int nn = NNodes();
		for (int i=1; i <= nn; i++)
		{
			u0 = AI*(GetNodePosInit2D(i)-p1) + uref1;
			urigid(GetNode(i).LTG(1)) = u0.X();
			urigid(GetNode(i).LTG(2)) = u0.Y();
		}
	}

}

void ACRSElement2D::ComputeDAijDq(int i, int j, Vector& v)
{
	v.SetLen(SOS());
	v.FillWithZeros();

	Vector2D p1 = GetNodePosInit2D(RefNode1());
	Vector2D p2 = GetNodePosInit2D(RefNode2());
	Vector2D u1 = GetNodeDisp2D(RefNode1());
	Vector2D u2 = GetNodeDisp2D(RefNode2());

	double p1x = p1.X();
	double p1y = p1.Y();
	double p2x = p2.X();
	double p2y = p2.Y();
	double u1x = u1.X();
	double u1y = u1.Y();
	double u2x = u2.X();
	double u2y = u2.Y();

	Vector2D dcsphi1, dcsphi2;

	double sign = 1; //for DsinphiDq
	if (j>i) sign = -1;
	if (i == j)
	{
		//DcosphiDq

		dcsphi1(1) = -(-2.0*u2y*p1y+2.0*p1y*u1y+2.0*p2y*u2y-2.0*p2y*p1y+p2y*p2y+u1y*u1y
			-2.0*u2y*u1y-2.0*p2y*u1y+u2y*u2y+p1y*p1y)/sqrt(pow(p2x*p2x+2.0*p2x*u2x-2.0*p2x*
			p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+u1x*u1x+p2y
			*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y*u1y+p1y*
			p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi1(2) = (-p2y-u2y+p1y+u1y)*(-p2x-u2x+p1x+u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x
			-2.0*p2x*p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+
			u1x*u1x+p2y*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y
			*u1y+p1y*p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi2(1) = (-2.0*u2y*p1y+2.0*p1y*u1y+2.0*p2y*u2y-2.0*p2y*p1y+p2y*p2y+u1y*u1y
			-2.0*u2y*u1y-2.0*p2y*u1y+u2y*u2y+p1y*p1y)/sqrt(pow(p2x*p2x+2.0*p2x*u2x-2.0*p2x*
			p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+u1x*u1x+p2y
			*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y*u1y+p1y*
			p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi2(2) = -(-p2y-u2y+p1y+u1y)*(-p2x-u2x+p1x+u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x
			-2.0*p2x*p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+
			u1x*u1x+p2y*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y
			*u1y+p1y*p1y+2.0*p1y*u1y+u1y*u1y,3.0));
	}
	else
	{
		//DsinphiDq

		dcsphi1(1) = (-p2y-u2y+p1y+u1y)*(-p2x-u2x+p1x+u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x
			-2.0*p2x*p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+
			u1x*u1x+p2y*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y
			*u1y+p1y*p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi1(2) = -(2.0*p2x*u2x+p1x*p1x+2.0*p1x*u1x+p2x*p2x-2.0*u2x*u1x-2.0*p2x*p1x
			-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x+u1x*u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x-2.0*p2x*
			p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+u1x*u1x+p2y
			*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y*u1y+p1y*
			p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi2(1) = -(-p2y-u2y+p1y+u1y)*(-p2x-u2x+p1x+u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x
			-2.0*p2x*p1x-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+
			u1x*u1x+p2y*p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y
			*u1y+p1y*p1y+2.0*p1y*u1y+u1y*u1y,3.0));
		dcsphi2(2) = (2.0*p2x*u2x+p1x*p1x+2.0*p1x*u1x+p2x*p2x-2.0*u2x*u1x-2.0*p2x*p1x-2.0
			*p2x*u1x+u2x*u2x-2.0*u2x*p1x+u1x*u1x)/sqrt(pow(p2x*p2x+2.0*p2x*u2x-2.0*p2x*p1x
			-2.0*p2x*u1x+u2x*u2x-2.0*u2x*p1x-2.0*u2x*u1x+p1x*p1x+2.0*p1x*u1x+u1x*u1x+p2y*
			p2y+2.0*p2y*u2y-2.0*p2y*p1y-2.0*p2y*u1y+u2y*u2y-2.0*u2y*p1y-2.0*u2y*u1y+p1y*p1y
			+2.0*p1y*u1y+u1y*u1y,3.0));
	}

	dcsphi1 *= sign;
	dcsphi2 *= sign;

	v(GetNode(RefNode1()).LTG(1)) = dcsphi1.X();
	v(GetNode(RefNode1()).LTG(2)) = dcsphi1.Y();
	v(GetNode(RefNode2()).LTG(1)) = dcsphi2.X();
	v(GetNode(RefNode2()).LTG(2)) = dcsphi2.Y();
}


void ACRSElement2D::DrawElement()
{

	mbs->SetColor(col);

	//UO() << "XGD(1)=" << XGD(1) << ", XGD(2)=" << XGD(2) << "\n";
/*
	Vector3D p1(GetPosD(Vector3D( 1.*size.X(), 1.*size.Y(),0)-GetNode(RefNode1()).Pos()));
	Vector3D p2(GetPosD(Vector3D(-0.*size.X(), 1.*size.Y(),0)-GetNode(RefNode1()).Pos()));
	Vector3D p3(GetPosD(Vector3D(-0.*size.X(),-0.*size.Y(),0)-GetNode(RefNode1()).Pos()));
	Vector3D p4(GetPosD(Vector3D( 1.*size.X(),-0.*size.Y(),0)-GetNode(RefNode1()).Pos()));

	//UO() << "p1=" << p1 << ",p3=" << p3 << "\n";

	double th = 2;
	mbs->MyDrawLine(p1,p2,th);
	mbs->MyDrawLine(p2,p3,th);
	mbs->MyDrawLine(p3,p4,th);
	mbs->MyDrawLine(p4,p1,th);
*/
};


void ACRSElement2D::PrintToAbaqus(ofstream& os)
{
	os.precision(17);

	IVector nodetrans(NNodes());
	IVector newnodes(0);
	IVector nodeused(NNodes());
	for (int i=1; i <= NNodes(); i++)
	{
		nodeused(i) = 0;
	}

	for (int i=1; i <= FFRFelements.Length(); i++)
	{
		Plate2D& pl = (Plate2D&)(GetMBS()->GetAuxElement(FFRFelements(i)));
		for (int j=1; j <= 8; j++)
		{
			nodeused(pl.NodeNum(j)) = 1;
		}
	}

	for (int i=1; i <= NNodes(); i++)
	{
		if (nodeused(i)) 
		{
			newnodes.Add(i);
			nodetrans(i) = newnodes.Length();
		}
		else
		{
			nodetrans(i) = 0;
		}
	}
	

	os << "*Heading\nSlidercrank\n**\n";
	os << "**PARTS\n";
	os << "**\n";
	os << "*Part, name=pleuel\n";
	os << "*End Part\n";
	os << "**\n";
	os << "**Material is E=5e8, nu=0.3, rho=5000\n";

	os << "**ASSEMBLY\n";
	os << "**\n";
	os << "*Assembly, name=Slidercrank\n";

	os << "*Instance, name=pleuel, part=pleuel\n";
	os << "*Node\n";
	for (int i=1; i <= newnodes.Length(); i++)
	{
		os << i << ", " << GetNode(newnodes(i)).Pos().X() << ", " << GetNode(newnodes(i)).Pos().Y() << "\n";
	}

	os << "*Element, type=CPS8,elset=el_pleuel\n";
	for (int i=1; i <= FFRFelements.Length(); i++)
	{
		Plate2D& pl = (Plate2D&)(GetMBS()->GetAuxElement(FFRFelements(i)));

		//8-noded element in Abaqus
		os << i << ", ";
		for (int j=1; j <= 8; j++)
		{
			os << nodetrans(pl.NodeNum(j));
			if (j != 8) os << ", ";
		}
		os << "\n";
	}
	os << "*Nset, nset=hinge1\n";
	for (int i=1; i <= boundarynode.Length(); i++)
	{
		os << nodetrans(boundarynode(i));
		if (i != boundarynode.Length()) 
		{
			if (i%16 == 0) os << "\n";
			else os << ", ";
		}
	}
	os << "\n";
	os << "*Solid Section, elset=el_pleuel,material=Material1\n";
	os << "1.,\n";
	os << "*End Instance\n";
	os << "*End Assembly\n";
	os << "*Material, name=Material1\n";
	os << "*Elastic\n";
	os << "500000000,0.3\n";

os << "** ----------------------------------------------------------------\n";
os << "** \n";
os << "** STEP: Step-1\n";
os << "** \n";
os << "*Step, name=Step-1\n";
os << "staticstep\n";
os << "*Static\n";
os << "1., 1., 1e-05, 1.\n";
os << "** \n";
os << "** \n";
os << "** LOADS\n";
os << "** \n";
os << "** Name: Load-1   Type: Body force\n";
os << "*Dload\n";
os << "PLEUEL.EL_PLEUEL, BY, -1000.\n";
os << "** \n";
os << "** OUTPUT REQUESTS\n";
os << "** \n";
os << "*Restart, write, frequency=1\n";
os << "** \n";
os << "** FIELD OUTPUT: F-Output-1\n";
os << "** \n";
os << "*Output, field, variable=PRESELECT\n";
os << "** \n";
os << "** HISTORY OUTPUT: H-Output-1\n";
os << "** \n";
os << "*Output, history, variable=PRESELECT\n";
os << "*El Print, freq=999999\n";
os << "*Node Print, freq=999999\n";
os << "*End Step\n";

}

