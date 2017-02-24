//#**************************************************************
//# filename:             computeeigenmodes.cpp
//#
//# author:								Astrid Pechstein and Karin Nachbagauer            
//#
//# generated:						
//# description:          
//# comments:	
//#							Computes generalized eigenvalues and eigenmodes of the Eigenvalue problem K*v = lambda * M * v
//#							Eigenmodes are loaded into Hotint Visualization
//#							three different options:
//#							 * direct computation: Use Lapack to compute all eigenvalues, full matrices are used
//#							 * sparse computation in Matlab:
//#							   - K and M are written in Matlab-(Sparse-)File
//#							   - gen. EVproblem is solved in Matlab (Arnoldi iteration)
//#							 * sparse computation by lobpcg iteration in Hotint
//#							   - gen. EVproblem is solved by Hotint implementation of lobpcg (Knyazev)
//#							 eigenvalues are loaded for Hotint visualization and can be stored using data manager
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


#include "mbs.h"


#include "windows.h"					//for shell execute
#include <direct.h>						//for getcwd
#include  <io.h>							//file operation _access
#include  <stdio.h>
#include  <stdlib.h>
#include "lapack_routines.h"

//$ YV 2012-12-14
//by now we explicitly reference particular elements:
// we need here to know, which degrees of freedom are constrained
#include "element.h"
//#include "constraint.h"
//#include "kinematicpairs.h"
//#include "specialconstraints.h"

int MultiBodySystem::ComputeEigenmodes()
{ 
	bComputeEigenmodes = 1; 
	bComputationIsInProgress = 1;
	bStopWhenPossible = 0;
	bPause = 0;

	SetProhibitRedraw(1);			//$ DR 2013-10-01 to avoid crash during computation (if animate scaling factor is active)

	SetComputationSolverOptions();

	int use_gyroscopic_terms = GetSolSet().eigsolv_gyro;

	//UO() << "old solution: =" << TIx0 << "\n";
	Vector sol = TIx0;

	//$ DR 2013-02-22 check if compute eigenmodes is called from a cpp-model or a hid-model
	mystr model_name = GetOptions()->GeneralOptions()->ModelFileInternalModelFunctionName();
	int nm = GetModelFunctionsIndex0();
	
	//$ SW 2013-11-5 [: in order to calculate the eigenmodes an initialization is necessary.
	// Doing an initialization at the original model would dismiss all the changes done since
	// the last calculation. Hence, if the model is from a HID file we can save the model
	// in a temporary hmc file and do the calculations there.
	// If the model is a cpp model we have to dismiss all changes anyway since we cannot
	// save cpp files.
	mystr hotint_input_data_file;
	if (nm == -1) //not an internal model
	{
		hotint_input_data_file = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.ModelFile.hotint_input_data_filename", "");

		//save current configuration
		ElementDataContainer edc;
		ElementData ed;

		mystr path = MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");
		mystr asv_model_name = "model_asv_tmp.hmc";
		ed.SetText(asv_model_name, "File_name"); 
		edc.Add(ed);
		ed.SetText(path, "Directory_name");
		edc.Add(ed);
		CallCompFunction(101, 1, 0, &edc); //save current configuration

		//set the model to be the temporary file and initialize this file
		GetMBS_EDC_Options()->TreeSetString("GeneralOptions.ModelFile.hotint_input_data_filename", path + asv_model_name);
		CallCompFunction(20,3); //Initialize(), each calculation of the eigensystem should start with XG=0 (otherwise Kmat wrong)
	}
	else
	{ // the current model is an internal model
		CallCompFunction(20,3); //Initialize(), each calculation of the eigensystem should start with XG=0 (otherwise Kmat wrong) 
	}
	//$ SW 2013-11-5:]


	//option = GetIOption(225);
	int option = GetSolSet().eigsolv_solvertype;
	UO(UO_LVL_ext) << "Chosen Option: " << option << "\n";
	

	// GSxs --> starting vector
	//TIx0 contains initial values
	mystatic Vector GSxs;
	//int ss = GetSystemSize();
	int sos= GetSecondOrderSize();//Number of 2nd order ODEs
	int es = GetFirstOrderSize(); //Number of 1st order ODEs
	int is = GetImplicitSize();   //Algebraic Equations (from boundary conditions) 

	//use actual solution for linearization (e.g. for pre-stressed problems):
	if (!GetSolSet().eigsolv_linearize_about_act_sol || sol.Length()!=TIx0.Length())
	{
		sol = TIx0;
	}
	else
	{
		UO(UO_LVL_ext) << "*********************************\n";
		UO(UO_LVL_ext) << "NOTE: Eigenvalues are computed for linearization around\nstored solution vector of last static/dynamic solution!\n";
		UO(UO_LVL_ext) << "*********************************\n";
	}

	GSxs.SetLen(sos+es+is);
  //GSxs is filled initial values
	for (int i=1; i <= sos; i++) GSxs(i) = sol(i);
	for (int i=1; i <= es; i++) GSxs(i+sos) = sol(i+2*sos);
	for (int i=1; i <= is; i++) GSxs(i+sos+es) = sol(i+2*sos+es);

/*
	mystatic Vector GSxs0; //with initial values
  //GSxs is filled with initial values
	GSxs0.SetLen(sos+es+is);
	for (int i=1; i <= sos; i++) GSxs0(i) = TIx0(i);
	for (int i=1; i <= es; i++) GSxs0(i+sos) = TIx0(i+2*sos);
	for (int i=1; i <= is; i++) GSxs0(i+sos+es) = TIx0(i+2*sos+es);
*/
	//UO() << "solution for linearization=" << sol << "\n";

	// Stiffness Matrix
	SparseMatrix Kmat; //sparse matrix build
	if (! use_gyroscopic_terms) // $ MSax 2013-06-25 : added because otherwise velocities whould be set to zero
	{
		StaticJacobian(Kmat,GSxs);//from starting vector GSxs the stiffness matrix Kmat is generated
		Kmat *= -1;
	}
	// Mass Matrix
	SparseMatrix Mmat;
	Mmat.SetSize((sos+es+is),(sos+es+is), MaxSparseBandwidth());
	if (! use_gyroscopic_terms) // $ MSax 2013-06-25 : added because otherwise velocities whould be set to zero
	{
		EvalM(GSxs, Mmat, 0);
	}
	else
	{
		EvalM(sol, Mmat, 0);
	}

	if (use_gyroscopic_terms  && (is != 0 || es != 0)) 
	{
		UO(UO_LVL_err) << "ERROR: Compute Eigenmodes: No constraints or first order equations allowed for damping or gyroscopic eigenvalue computation\n";
		use_gyroscopic_terms = 0;
	}
	if (use_gyroscopic_terms  && (sos > 1000)) 
	{
		UO(UO_LVL_err) << "ERROR: Compute Eigenmodes: system size must be smaller than 1000\n";
		use_gyroscopic_terms = 0;
	}

	//SparseMatrix Dmat;				//$ DR 2013-10-21 moved this matrix inside the if(use_gyro..) branch
	//Matrix Id2(sos*2,sos*2);	//$ DR 2013-10-21 removed this matrix, because not used at all
	//Matrix A(sos*2,sos*2);		//$ DR 2013-10-21 size of matrix is set now inside the if(use_gyro..) branch
	Matrix A;

	if (use_gyroscopic_terms) // MSax: use gyroscopic terms
	{
		SparseMatrix Dmat;
		A.SetSize(sos*2,sos*2);

		Kmat.SetSize((sos+es+is),(sos+es+is), MaxSparseBandwidth());
		Dmat.SetSize((sos+es+is),(sos+es+is), MaxSparseBandwidth());
		ComputeStiffnessAndDampingMatrix(Kmat, Dmat, sol);

		SparseMatrix GYmat;
		ComputeGyroscopicMatrix(GYmat); // GY..Gyroscopic matrix

		Matrix Mfull; 
		Mmat.GetMatrix(Mfull);

		Matrix Kfull; 
		Kmat.GetMatrix(Kfull);

		Matrix GYfull; 
		GYmat.GetMatrix(GYfull);

		Matrix Id1(sos,sos); 
		Id1.SetDiagMatrix(1.);

		// build  A = [0,E;-M^(-1)*K,-M^(-1)*G] matrix, see HID_USER docu
		A.SetSize(sos*2,sos*2);
		A.SetSubmatrix(Id1,1,1+sos);
		
		int rv=Mfull.Invert2();
		if (!rv) UO(UO_LVL_err) << "ERROR: ComputeEigenmodes: mass matrix not invertable!\n";

		Kfull= Mfull * Kfull;

		GYfull= Mfull * GYfull;

		A.SetSubmatrix(GYfull,1+sos,1+sos);
		A.SetSubmatrix(Kfull,1+sos,1);

	}
	int neig;
	int size;
	int numberzeromodes = 0;
	IVector unconstraineddofs(0);
	if (!use_gyroscopic_terms) // $ MSax 2013-07-04 : added
	{
		//$ DR 2013-02-22 error handling added
		if(Kmat.Getcols() * Mmat.Getcols()==0)
		{
			UO(UO_LVL_err) << "\n ERROR: The K- and/or M-matrix has size == 0! Does the system contain elements at all? \n Stopped computation of eigenmodes.\n";
			bComputationIsInProgress = 0;
			SetProhibitRedraw(0); //$ DR 2013-10-01 to avoid crash during computation (if animate scaling factor is active)

			//reset the filename to the original file name if changed
			if (nm == -1)
			{
				GetMBS_EDC_Options()->TreeSetString("GeneralOptions.ModelFile.hotint_input_data_filename", hotint_input_data_file);
			}
			return 0; 
		}

		// List of constrained dofs, which are eliminated from matrix
		TArray<int> CClist;
		CClist.SetLen(0);

		for(int i=1; i <= GetNElements(); i++)//pass all elements
		{
			Element& elem = GetElement(i);//set pointer elem to each element
			//$ YV 2013-01-03: constrained dofs are provided by the basis Element class, so that there is
			// no need to explicitly differentiate Element types here
			if(elem.IS() != 0)		// we consider only the constraints with Lagrange multipliers
				elem.GetGlobalConstrainedDOFs(CClist);
	#if 0
			//we look for CoordinateConstraints:
			CoordConstraint* ccelem = dynamic_cast<CoordConstraint*>(&elem);//we ask, if elem is of type CoordConstraint
			//->0, if elem is not a CoordConstraint, else ->Pointer to elem, if type is CoordConstraint
			if(ccelem!=0)//if pointer is part of ccelem, pointer is added to the list ConstrDOF
			{
				if (ccelem->IS() == 0) continue; // only coord-constraints with Lag.Mult. (IS==1) are of interest, cc with IS==0 are penalty, matrix row is kept
				int GlobConstrNr;
				GlobConstrNr=ccelem->GetGlobalConstrainedDOF();//from local number ccelem a global DOFnumber is generated
				//matout << GlobConstrNr << " ";//add global DOF-Number to the list
				CClist.Add(GlobConstrNr);
			}
			// $LA 2011-08-25:[
			//we look for NodalConstraints:
			NodalConstraint* ncelem = dynamic_cast<NodalConstraint*>(&elem);//we ask, if elem is of type NodalConstraint
			//->0, if elem is not a NodalConstraint, else ->Pointer to elem, if type is NodalConstraint
			if(ncelem!=0)//if pointer is part of ncelem, pointer is added to the list ConstrDOF
			{
				if (ncelem->IS() == 0) continue; // only nodal-constraints with Lag.Mult. (IS==1) are of interest, cc with IS==0 are penalty, matrix row is kept
				TArray<int> GlobConstrNrs;

				// LGA global constrainedDOF ausfindig machen!!!
				ncelem->GetGlobalConstrainedDOFs(GlobConstrNrs);//from local number ccelem a global DOFnumber is generated

				//matout << GlobConstrNr << " ";//add global DOF-Number to the list
				for(int j=1; j<= GlobConstrNrs.Length(); j++)
				{
					CClist.Add(GlobConstrNrs(j));
				}
			}
			// $LA 2011-08-25:]
	#endif
		}
		Quicksort(CClist);
		// Vector containing the degrees of freedom which are not constrained due to coord constraints
		for (int i=1; i<=sos; i++)
		{
			if (!CClist.Find(i))
				unconstraineddofs.Add(i);
		}
		//number of calculated eigenvectors
		if (option==0)
		{
			neig = unconstraineddofs.Length();
			// if (use_gyroscopic_terms) neig = 2*sos; $ MSax 2013-07-04 : removed
		}
		else
		{
			neig = GetSolSet().eigsolv_n_eigvals;
		}
		size = unconstraineddofs.Length();
	}
	else
	{
		neig = 2*sos; // $ MSax 2013-07-04 : added
	}

	UO(UO_LVL_ext) << "Number of eigenvalues computed =" << neig << "\n";

	int computedcorrectly = 0;
	// Matrix for eigenvectors, eigenvalues
	Matrix EigenvektorMatrix;
	//Vector eigval;	//$ RL 2011-02: mbs class variable eigval defined
	eigval.SetLen(neig);
	eigval.SetAll(0.);

	if (use_gyroscopic_terms) // $ MSax: added
	{
		EigenvektorMatrix.SetSize(neig,2*sos);
		EigenvektorMatrix.SetAll(0);

		Vector work(8*2*sos+16); // work vector for calculation (necessary for LAPACK)
		Vector wr(2*sos); // vector with real part of eigenvalues
		Vector wi(2*sos); // vector with imaginary parts of eigenvalues
		Matrix ev(2*sos,2*sos); // matrix with eigenvectors

		int lapackinfo = LapackGenEVP_CGEEV(2*sos, &A(1,1), &wr(1), &wi(1), &ev(1,1), &work(1), work.Length()); 
		if (lapackinfo) computedcorrectly=0;
		else computedcorrectly=1;

		TArray<double> evalues(0);
		evalues.SetLen(2*sos);
		evalues.SetAll(0);
		TArray<int> index(0);
		index.SetLen(2*sos);
		index.SetAll(0);

		for(int j=1; j<=2*sos; j++) 
		{
			evalues(j) = sqrt(sqr(wr(j))+sqr(wi(j))); // eigenvalues = (alphar + i alphai) / beta
			index(j) = j;
		}
		QuicksortDouble(evalues, index);	

		for(int j=1; j<=2*sos; j++)
		{
			eigval(j) = evalues(j);
			for(int k=1; k<=2*sos; k++)
			{
				EigenvektorMatrix(j,k) = ev(index(j),k);
			}
		}
	}
	else
	{
		EigenvektorMatrix.SetSize(neig,sos);
		EigenvektorMatrix.SetAll(0);
		
		//Solver Parameter (optional)
		int maxit = GetSolSet().eigsolv_maxiterations; //GetIOption(224);
		double tol = GetSolSet().eigsolv_accuracy; //GetDOption(225);
		numberzeromodes = GetSolSet().eigsolv_nzeromodes; //GetIOption(226);
		int use_preconditioning = GetSolSet().eigsolv_use_preconditioning; //GetIOption(228);
		double lambda_precond = GetSolSet().eigsolv_precond_lambda; //GetDOption(226);

		int init_random=1;
		SparseEigenvalueSolver evsolver(this);
		
		//evsolver.Set(&Kmat, &Mmat, maxit, tol, init_random);
		evsolver.Set(&Kmat, &Mmat, maxit, tol, init_random);

		// evsolver can be stopped by MBS
		evsolver.CanBeStopped() = 1;
		if (use_preconditioning)
			evsolver.SetPreconditioner(lambda_precond);
		//   MATLAB SPARSE SOLVER
		if (option == 1)
		{
			mystr pathMatlab="../../Matlab/ComputeEigensystem";
			computedcorrectly = evsolver.ComputeEigenModesMatlab(pathMatlab, neig, eigval, EigenvektorMatrix, unconstraineddofs,
															numberzeromodes);	
		}
		else if (option == 0) // Direct solver by LAPACK
		{
			Matrix Kfull, Mfull;
			neig = size; // compute all possible dofs
			Kfull.CopyFrom(Kmat,1,1,size,size, unconstraineddofs);
			Mfull.CopyFrom(Mmat,1,1,size,size, unconstraineddofs);
			if(UO().GetGlobalMessageLevel()==UO_LVL_dbg1)	//$ DR 2011-09-15
			{
				UO(UO_LVL_dbg1) << "K=" << Kfull << "\n";
				UO(UO_LVL_dbg1) << "M=" << Mfull << "\n";
			}
			Vector work(4*size);

			int lapackinfo = LapackGenEVPSPD(size, &Kfull(1,1), &Mfull(1,1), &eigval(1), &work(1), work.Length());
			if (lapackinfo) computedcorrectly=0;
			else computedcorrectly=1;

			// eigenvectors are stored rowwise in Matrix Kfull
			for (int j=1; j<=neig; j++)
			{
				for (int i=1; i<=size; i++)
				{
					EigenvektorMatrix(j,unconstraineddofs(i)) = Kfull(j,i);
				}
			}
		}
		else if (option == 2) // lobpcg - hack
		{
			eigval.SetLen(neig);
			computedcorrectly = evsolver.ComputeEigenModes(neig, eigval, EigenvektorMatrix, unconstraineddofs, numberzeromodes); 
		}
	}

	//normalize eigenvectors:
	//UO(UO_LVL_all) << "EigenvektorMatrix=" << EigenvektorMatrix << "\n";			//$ DR 2011-05-27

	if(GetSolSet().eigsolv_normalization_mode==0)		//$ DR 2011-05-27: normalize such that max(abs(v)) = 1
	{
		for(int j=1; j <= neig; j++)//number of eigenvectors
		{
			double maxval = 0;
			for(int i=1; i <= sos; i++)//length of eigenvectors
			{
				maxval = Maximum(maxval, fabs(EigenvektorMatrix(j, i)));
			}

			if (maxval!=0)
			{
				if (use_gyroscopic_terms) // $ MSax 2013-07-04 : added
				{
					for(int i=1; i <= sos*2; i++)//length of eigenvectors
					{
						EigenvektorMatrix(j, i) /= maxval;
					}
				}
				else
				{
					for(int i=1; i <= sos; i++)//length of eigenvectors
					{
						EigenvektorMatrix(j, i) /= maxval;
					}
				}
			}
		}
	}
	else if(GetSolSet().eigsolv_normalization_mode==1)	// normalize such that v'*v = 1
	{
		for(int j=1; j <= neig; j++)//number of eigenvectors
		{
			double norm = 0;
			for(int i=1; i <= sos; i++)//length of eigenvectors
			{
				norm += Sqr(EigenvektorMatrix(j, i));
			}
			if (norm!=0)
			{
				norm = 1./sqrt(norm);
				if (use_gyroscopic_terms) // $ MSax 2013-07-04 : added
				{
					for(int i=1; i <= sos*2; i++)//length of eigenvectors
					{
						EigenvektorMatrix(j, i) *= norm;
					}
				}
				else
				{
					for(int i=1; i <= sos; i++)//length of eigenvectors
					{
						EigenvektorMatrix(j, i) *= norm;
					}
				}
			}
		}
	}

	if (use_gyroscopic_terms) sos = sos*2; // $ MSax 2013-07-04 : added


	//UO(UO_LVL_all) << "EigenvektorMatrix (normalized)=" << EigenvektorMatrix << "\n"; //$ DR 2011-05-27
	EigenvektorMatrix *= GetSolSet().eigsolv_eigenmodes_scaling_factor; //scale eigenvectormatrix

	// eigenmode is eigenvector from linearized problem + actual state
	for(int i=1; i <= sos; i++)//length of eigenvectors
	{
		for(int j=1; j <= neig; j++)//number of eigenvectors
		{
			EigenvektorMatrix(j, i) += GetXact(i);
		}
	}

	if (!use_gyroscopic_terms) // $ MSax 2013-07-04: otherwise already in rad/seconds
	{
		// eigenfrequency in rad/s
		for(int j=1; j <= neig; j++)//number of eigenvectors
		{
			eigval(j) = sqrt(fabs(eigval(j)));
		}
	}


	//eigenvectors are sent to the DataManager (for automized animation without file-loading)
	int orig_storedata = GetSolSet().storedata;
	for(int i=1; i<=neig; i++)
	{
		SetTime(i);
		Vector& TIx0draw2 = GetSolVector();
		//Vector& TIx0draw2 = GetDrawVector();
		TIx0draw2.SetLen(sos);
		TIx0draw2.SetAll(0.);
		for(int j=1;j<=sos;j++)
		{
			TIx0draw2(j) = EigenvektorMatrix(i,j);
		}
		drawnow = 0;
		GetSolSet().storedata = -2;
		TIDrawAndStore();
	}

	GetSolSet().storedata = orig_storedata;

	//eigenvalues and eigenvectors are written to finished computationComputation Output:

	//Output: eigenvalues in a list, eigenvectors in matrixform (column-wise)
	//UO() << "Eigenfrequencies: " << eigval << "EigenvectorMatrix: " << EigenvektorMatrix << "\n";

	//Output: Eigenfrequencies
	if (computedcorrectly)
	{
		UO(UO_LVL_all) << "Eigenvalues computed up to prescribed tolerance\n";
	}
	else
	{
		UO(UO_LVL_err) << "Eigenvalue computation not successful!\n";
	}
	if (UO().GetGlobalMessageLevel() > UO_LVL_multsim)
	{
		UO(UO_LVL_all) << "\n" << "Eigenfrequencies:\n";
		for(int i=1; i<=neig; i++)
		{
			UO(UO_LVL_all) << i <<", "<< eigval(i)/(2.*MY_PI)<<" Hz = " << eigval(i) << "rad/s\n";
		}
	}
	else
	{
		UO(UO_LVL_multsim) << "\n" << "Eigenfrequencies (Hz):";
		for(int i=1; i<=neig; i++)
		{
			UO(UO_LVL_multsim) << eigval(i)/(2.*MY_PI);
			if (i < neig) UO(UO_LVL_multsim) << ", ";
		}
		UO(UO_LVL_multsim) << "\n" << "Eigenfrequencies (rad):";
		for(int i=1; i<=neig; i++)
		{
			UO(UO_LVL_multsim) << eigval(i);
			if (i < neig) UO(UO_LVL_multsim) << ", ";
		}
	}
	//Message for user, if claimed number of zeromodes is not reached:
	if(numberzeromodes!=0)
	{
		if(eigval(numberzeromodes)>0.1)
		{
			UO(UO_LVL_err) << "\n" << "Number of Zeromodes NOT obtained!" << "\n";
		}
	}

	// DR 2012-04-19: flag to change format of output added	
	int outp_format_flag = GetSolSet().eigval_outp_format_flag; // 1 .. eigenfreq., 2 .. eigenvec., 4 .. eigenfreq. in Hz (otherwise in rad/s)

	// output to file
	mystr outputfile = GetSolSet().sol_directory + mystr("eigenmodes_") + GetSolSet().sol_filename;
	ofstream output(outputfile);

	if(outp_format_flag & 1)
	{
		for(int i=1; i<=neig; i++)
		{
			if(outp_format_flag & 4) 
			{
				output << eigval(i)/2./MY_PI <<"\n";
			}
			else
			{
				output << eigval(i) <<"\n";				
			}
		}
	}

	if(outp_format_flag & 2)
	{
		for(int j=1; j<=sos; j++)
		{
			for(int i=1; i<=neig; i++)
			{
				output << EigenvektorMatrix(i,j) <<"\t";
			}
			output << "\n";
		}
	}

	output.close();

	bComputeEigenmodes = 0;
	bComputationIsInProgress = 0;

	SetProhibitRedraw(0); //$ DR 2013-10-01 to avoid crash during computation (if animate scaling factor is active)

	if (computedcorrectly)
	{
		uo.pUI->CallWCDriverFunction(7);
	}

	//reset the filename to the original file name if changed
	if (nm == -1)
	{
		GetMBS_EDC_Options()->TreeSetString("GeneralOptions.ModelFile.hotint_input_data_filename", hotint_input_data_file);
	}
	return 0;
}
