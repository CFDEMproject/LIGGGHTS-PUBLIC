//#**************************************************************
//# filename:             ComputeEigenmodes.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          ComputeEigenmodes-Dialog
//# comments:
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
 

#include "stdafx.h"
#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"
#include "ComputeEigenmodes.h"


extern UINT ThreadControlFunc(LPVOID pVOID);

IMPLEMENT_DYNAMIC(ComputeEigenmodes, CDialog)
ComputeEigenmodes::ComputeEigenmodes(CWnd* pParent /*=NULL*/)
	: CDialog(ComputeEigenmodes::IDD, pParent)
	, m_radio_compute_eigenmode_solver(0), pWCDI(0), pGLDrawWnd(0)
	, m_number_computed_eigenvalues(3)//Defaultvalue = compute the first 3 eigenvalues
	, m_number_maxiterations(100000)
	, m_convergence_tolerance(0.001)
	, m_number_zeromodes(0)
	, m_check_number_zeromodes(0)
	, m_check_usepreconditioner(0)
	, m_preconditioner_lambda(1)
{
	EnableAutomation();
}

ComputeEigenmodes::~ComputeEigenmodes()
{
}

void ComputeEigenmodes::OnFinalRelease()
{
	// Wenn der letzte Verweis auf ein Automatisierungsobjekt freigegeben wird,
	// wird OnFinalRelease aufgerufen. Die Basisklasse löscht das Objekt
	// automatisch. Fügen Sie zusätzlichen Bereinigungscode für Ihr Objekt
	// hinzu, bevor Sie die Basisklasse aufrufen.

	CDialog::OnFinalRelease();
}

//LoadData() 
void ComputeEigenmodes::LoadData()
{
	// So machen:
	pWCDI->MBS_EDC_TreeGetInt(m_number_computed_eigenvalues, "SolverOptions.Eigensolver.n_eigvals");
	pWCDI->MBS_EDC_TreeGetInt(m_number_maxiterations, "SolverOptions.Eigensolver.max_iterations");
	pWCDI->MBS_EDC_TreeGetDouble(m_convergence_tolerance, "SolverOptions.Eigensolver.accuracy");
	pWCDI->MBS_EDC_TreeGetInt(m_number_zeromodes, "SolverOptions.Eigensolver.n_zero_modes");
	pWCDI->MBS_EDC_TreeGetInt(m_check_number_zeromodes, "SolverOptions.Eigensolver.use_n_zero_modes");
	pWCDI->MBS_EDC_TreeGetInt(m_check_usepreconditioner, "SolverOptions.Eigensolver.use_preconditioning");
	pWCDI->MBS_EDC_TreeGetInt(m_radio_compute_eigenmode_solver, "SolverOptions.Eigensolver.solver_type");
	pWCDI->MBS_EDC_TreeGetDouble(m_preconditioner_lambda, "SolverOptions.Eigensolver.preconditioner_lambda");
	//m_number_computed_eigenvalues = pWCDI->GetIOption(222); //stores IOption(222)-Value in m_number_computed_eigenvalues, default=3
  //m_number_maxiterations = pWCDI->GetIOption(224); 
	//m_radio_compute_eigenmode_solver = pWCDI->GetIOption(225);
	//m_convergence_tolerance = pWCDI->GetDOption(225); 
	//m_preconditioner_lambda = pWCDI->GetDOption(226); 
	//m_number_zeromodes = pWCDI->GetIOption(226); 
	//m_check_number_zeromodes = pWCDI->GetIOption(227);
	//m_check_usepreconditioner = pWCDI->GetIOption(228);
}

//WriteData() 
void ComputeEigenmodes::WriteData()
{
	pWCDI->MBS_EDC_TreeSetInt(m_number_computed_eigenvalues, "SolverOptions.Eigensolver.n_eigvals");
	pWCDI->MBS_EDC_TreeSetInt(m_number_maxiterations, "SolverOptions.Eigensolver.max_iterations");
	pWCDI->MBS_EDC_TreeSetDouble(m_convergence_tolerance, "SolverOptions.Eigensolver.accuracy");
	pWCDI->MBS_EDC_TreeSetInt(m_number_zeromodes, "SolverOptions.Eigensolver.n_zero_modes");
	pWCDI->MBS_EDC_TreeSetInt(m_check_number_zeromodes, "SolverOptions.Eigensolver.use_n_zero_modes");
	pWCDI->MBS_EDC_TreeSetInt(m_check_usepreconditioner, "SolverOptions.Eigensolver.use_preconditioning");
	pWCDI->MBS_EDC_TreeSetInt(m_radio_compute_eigenmode_solver, "SolverOptions.Eigensolver.solver_type");
	pWCDI->MBS_EDC_TreeSetDouble(m_preconditioner_lambda, "SolverOptions.Eigensolver.preconditioner_lambda");
	//pWCDI->GetIOption(222) = m_number_computed_eigenvalues;//writes m_number_computed_eigenvalues in IOption(222)
	//pWCDI->GetIOption(224) = m_number_maxiterations;
	//pWCDI->GetDOption(225) = m_convergence_tolerance;
	//pWCDI->GetIOption(226) = m_number_zeromodes;
	//pWCDI->GetIOption(227) = m_check_number_zeromodes;
	//pWCDI->GetIOption(228) = m_check_usepreconditioner;
	//pWCDI->GetIOption(225) = m_radio_compute_eigenmode_solver;
	//pWCDI->GetDOption(226) = m_preconditioner_lambda;
}

void ComputeEigenmodes::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO_EIGENMODE_MATLAB_DIRECT, m_radio_compute_eigenmode_solver);
	//the value from dialog-box is stored in variable m_number_computed_eigenvalues:
	DDX_Text(pDX, IDC_NUMBEREV, m_number_computed_eigenvalues);
	//Solver Parameter:
	DDX_Text(pDX, IDC_EV_MAXIT, m_number_maxiterations);
	DDX_Text(pDX, IDC_EV_TOL, m_convergence_tolerance);
	DDX_Check(pDX, IDC_CHECK_NR_ZEROMODES, m_check_number_zeromodes);
	DDX_Text(pDX, IDC_EV_NR_ZEROMODES, m_number_zeromodes);
	DDX_Check(pDX, IDC_CHECK_EV_USE_PRECONDITIONER, m_check_usepreconditioner);
	DDX_Text(pDX, IDC_EDIT_EV_PRECONDITIONER_LAMBDA, m_preconditioner_lambda);
	
	DDV_MinMaxInt(pDX, m_number_computed_eigenvalues, 1, 2000000000);//range for the value
}

BEGIN_MESSAGE_MAP(ComputeEigenmodes, CDialog)
	ON_BN_CLICKED(IDOK, &ComputeEigenmodes::OnBnClickedOk)
	ON_BN_CLICKED(IDC_CHECK_NR_ZEROMODES, &ComputeEigenmodes::OnActivateZeromodes)
	ON_BN_CLICKED(IDC_CHECK_EV_USE_PRECONDITIONER, &ComputeEigenmodes::OnActivatePreconditioner)
	ON_BN_CLICKED(IDC_RADIO_EIGENMODE_MATLAB_DIRECT, &ComputeEigenmodes::OnIterativeComputation)
	ON_BN_CLICKED(IDC_RADIO_EIGENMODE_MATLAB_SPARSE, &ComputeEigenmodes::OnIterativeComputation)
	ON_BN_CLICKED(IDC_RADIO_EIGENMODE_LOBPCG_SPARSE, &ComputeEigenmodes::OnIterativeComputation)
	ON_MESSAGE(WM_CLOSE_EV_WINDOW, OnCloseWindow)
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(ComputeEigenmodes, CDialog)
END_DISPATCH_MAP()


// {B8FB3003-54C8-4ABA-858E-01CAA61057AB}
static const IID IID_IomputeEigenmodes =
{ 0xB8FB3003, 0x54C8, 0x4ABA, { 0x85, 0x8E, 0x1, 0xCA, 0xA6, 0x10, 0x57, 0xAB } };

BEGIN_INTERFACE_MAP(ComputeEigenmodes, CDialog)
	INTERFACE_PART(ComputeEigenmodes, IID_IomputeEigenmodes, Dispatch)
END_INTERFACE_MAP()

void ComputeEigenmodes::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	LoadData();

	// Zeromodes:
	CButton* check = (CButton*) GetDlgItem(IDC_CHECK_NR_ZEROMODES);
	CEdit* edit = (CEdit*) GetDlgItem(IDC_EV_NR_ZEROMODES);
	// set/unset checkbox
	check->SetCheck(m_check_number_zeromodes);//Check Box deaktivate
	// activate/deactivate edit box
	if (!m_check_number_zeromodes)
		edit->EnableWindow(0);//Edit Box deactivated
	else
		edit->EnableWindow(1);//Edit Box activated


	// solver parameter values
	SetDlgItemInt(IDC_NUMBEREV,m_number_computed_eigenvalues);
	SetDlgItemInt(IDC_EV_MAXIT,m_number_maxiterations);
	char buffer[50];
	sprintf(buffer, "%.0e",m_convergence_tolerance);
	SetDlgItemTextA(IDC_EV_TOL,buffer);//no SetDlgItemDouble available
	SetDlgItemInt(IDC_EV_NR_ZEROMODES,m_number_zeromodes);

	//// activate/deactivate solver parameter edit boxes
	EnableDisableSolverParameters();

	UpdateData(FALSE);

	ShowWindow(SW_SHOW);
}

UINT ThreadControlFuncEM(LPVOID pVOID)
{
	ComputeEigenmodes* pCE = (ComputeEigenmodes*)pVOID;
	pCE->GetWCDriver3DDlg()->GetWCDInterface()->CallCompFunction(320, 0);//send message to WCD interface to compute eigenmodes
	// set Text to Compute Eigenmodes in case that window is still open (i.e. in case that computation not converged)
	if(::IsWindow(pCE->m_hWnd))
		pCE->GetDlgItem(IDOK)->SetWindowText("COMPUTE EIGENMODES");
	return 0;
}

void ComputeEigenmodes::OnBnClickedOk()
{
	WCDDlg->GetWCDInterface()->SetPCFB(WCDDlg);
	UpdateData(TRUE);
	WriteData();//IOptions get the right values
	
	if(!pWCDI->IsComputationInProgress())
	{
		CWinThread * pWinThread  = AfxBeginThread(ThreadControlFuncEM,this);
		SetThreadPriority(pWinThread->m_hThread,THREAD_PRIORITY_BELOW_NORMAL);
		GetDlgItem(IDOK)->SetWindowText("Stop computation");
	}
	else
	{
		GetDlgItem(IDOK)->SetWindowText("Stopping..");
		pWCDI->StopWhenPossible();
		GetDlgItem(IDOK)->SetWindowText("COMPUTE EIGENMODES!");
	}
	
	OnOK();
}


//Specifying number of zeromodes or not:
void ComputeEigenmodes::OnActivateZeromodes()
{
	CButton* check = (CButton*) GetDlgItem(IDC_CHECK_NR_ZEROMODES);
	CEdit* edit = (CEdit*) GetDlgItem(IDC_EV_NR_ZEROMODES);
	if(check->GetCheck() == 1)//check
	{
		check->SetCheck(0);//no check
		edit->EnableWindow(0);//edit box deactivated -> no input possible
		SetDlgItemInt(IDC_EV_NR_ZEROMODES,0);//no entry -> default-value=0
	}
	else if(check->GetCheck() == 0)//no check
	{
		check->SetCheck(1);//check
		edit->EnableWindow(1);//edit box activated -> input data
	}
}

void ComputeEigenmodes::OnActivatePreconditioner()
{
	CButton* check = (CButton*) GetDlgItem(IDC_CHECK_EV_USE_PRECONDITIONER);
	CEdit* edit = (CEdit*) GetDlgItem(IDC_EDIT_EV_PRECONDITIONER_LAMBDA);
	if(check->GetCheck() == 1)//check
	{
		//check->SetCheck(0);//no check
		edit->EnableWindow(1);//edit box deactivated -> no input possible
	}
	else if(check->GetCheck() == 0)//no check
	{
		//check->SetCheck(1);//check
		edit->EnableWindow(0);//edit box activated -> input data
	}
}

//Specifying solver parameter options or not:
void ComputeEigenmodes::OnIterativeComputation()
{
	WCDDlg->GetWCDInterface()->SetPCFB(WCDDlg);
	UpdateData(TRUE);

	EnableDisableSolverParameters();

}

// enable or disable solver parameter edit boxes depending on chosen method
void ComputeEigenmodes::EnableDisableSolverParameters()
{
	int method = m_radio_compute_eigenmode_solver;
	CEdit* edit_tol = (CEdit*) GetDlgItem(IDC_EV_TOL);
	CEdit* edit_maxit = (CEdit*) GetDlgItem(IDC_EV_MAXIT);
	CEdit* edit_nev = (CEdit*) GetDlgItem(IDC_NUMBEREV);
	CEdit* edit_lambda = (CEdit*) GetDlgItem(IDC_EDIT_EV_PRECONDITIONER_LAMBDA);
	CButton* check_lambda = (CButton*) GetDlgItem(IDC_CHECK_EV_USE_PRECONDITIONER);

	if(method == 0) // direct method
	{
		edit_tol->EnableWindow(0);//edit box deactivated -> no input possible
		edit_maxit->EnableWindow(0);
		edit_nev->EnableWindow(0);
		edit_lambda->EnableWindow(0);
		check_lambda->EnableWindow(0);
	}
	else if (method == 1) // iterative method matlab
	{
		edit_tol->EnableWindow(1);//edit box activated -> input possible
		edit_maxit->EnableWindow(1);
		edit_nev->EnableWindow(1);
		edit_lambda->EnableWindow(0);
		check_lambda->EnableWindow(0);
	}
	else if (method == 2) // iterative method lobpcg
	{
		edit_tol->EnableWindow(1);//edit box activated -> input possible
		edit_maxit->EnableWindow(1);
		edit_nev->EnableWindow(1);
		check_lambda->EnableWindow(1);
		if (m_check_usepreconditioner)
			edit_lambda->EnableWindow(1);
		else
			edit_lambda->EnableWindow(0);
	}
}

void ComputeEigenmodes::OnOK()
{
	//DestroyWindow();
}

void ComputeEigenmodes::OnClose() 
{
	DestroyWindow();
}

void ComputeEigenmodes::OnCancel() 
{
	DestroyWindow();
}

LRESULT ComputeEigenmodes::OnCloseWindow(WPARAM, LPARAM)
{
	DestroyWindow();

	// Eat spurious WM_CLOSE_EV_WINDOW messages
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_CLOSE_EV_WINDOW, WM_CLOSE_EV_WINDOW, PM_REMOVE));

	return 0;
}