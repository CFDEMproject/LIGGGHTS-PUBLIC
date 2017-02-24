//#**************************************************************
//# filename:             DlgOneEditCtrl.h
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          
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

const int MAXCHAR_ONEEDITTEXT = 512;
// CDialogOneEditControl-Dialogfeld

class CDialogOneEditControl : public CDialog
{
	DECLARE_DYNAMIC(CDialogOneEditControl)

public:
	CDialogOneEditControl(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CDialogOneEditControl();

// Dialogfelddaten
	enum { IDD = IDD_DIALOGONEEDITCONTROL };

	CString statictext;
	CString dialogname;
	CString edittext;
	double* value;
	int* intvalue;
	int type; //1==int, 2==double, 3=text, 4=watch
	double minval, maxval;
	int idok;
	int sensornumber;
	int sensorwatchdialognumber;

	void InitDouble(const char* dialognameI, const char* statictextI, double* outputval, double min=-1e100, double max=1e100);
	
	void InitInt(const char* dialognameI, const char* statictextI, int* outputval, double min=-1e100, double max=1e100);

	void InitText(const char* dialognameI, const char* statictextI);

	void InitWatch(const char* dialognameI, int sensornumberI, int sensorwatchdialognumberI); //control to view certain numbers of the simulation

	int GetSensorNumber() const {return sensornumber;}
	void SetSensorText(const char* str);

	void Create(CWnd * pParent);
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }

protected:
	void OnOK();
	void OnCancel();
	void OnClose();

	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	virtual BOOL OnInitDialog();

private:
	WCDInterface* pWCDI; //only set for watch
	CGLDrawWnd* pGLDrawWnd; //only set for watch

};
