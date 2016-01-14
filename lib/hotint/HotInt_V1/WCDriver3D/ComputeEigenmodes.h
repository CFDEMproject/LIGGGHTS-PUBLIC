//#**************************************************************
//# filename:             ComputeEigenmodes.h
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
 

#pragma once

class CWCDriver3DDlg;

// ComputeEigenmodes-Dialogfeld

class ComputeEigenmodes : public CDialog
{
	DECLARE_DYNAMIC(ComputeEigenmodes)

public:
	ComputeEigenmodes(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~ComputeEigenmodes();

	virtual void OnFinalRelease();

// Dialogfelddaten
	enum { IDD = IDD_COMPUTEEIGENMODES };

	//functions not added by windows:
	void Create(CWnd * pParent);

	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	void SetWCDDlg(CWCDriver3DDlg* adial) {WCDDlg = adial;}
	void LoadData();  //Get data from WCDinterface
	void WriteData(); //Put data to WCDinterface

	CWCDriver3DDlg* GetWCDriver3DDlg() {return WCDDlg; }

private:
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;

	CWCDriver3DDlg* WCDDlg;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()
	DECLARE_DISPATCH_MAP()
	DECLARE_INTERFACE_MAP()

	void OnOK();
	void OnCancel();
	void OnClose();
	//void OnActivate();

public:
	int m_radio_compute_eigenmode_solver;
	int m_number_computed_eigenvalues;
	int m_number_maxiterations;
	double m_convergence_tolerance;
	int m_number_zeromodes;
	int m_check_number_zeromodes;
	int m_check_usepreconditioner;
	double m_preconditioner_lambda;
public:
	afx_msg void OnBnClickedOk();
  afx_msg void OnActivateZeromodes();
  afx_msg void OnActivatePreconditioner();
	afx_msg void OnIterativeComputation();
	afx_msg void EnableDisableSolverParameters();
	afx_msg LRESULT OnCloseWindow(WPARAM, LPARAM);
	//afx_msg void WriteValue();
};
