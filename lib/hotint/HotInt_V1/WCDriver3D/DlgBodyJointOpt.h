//#**************************************************************
//# filename:             DlgBodyJointOpt.h
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


// DialogBodyJointOptions-Dialogfeld

class DialogBodyJointOptions : public CDialog
{
	DECLARE_DYNAMIC(DialogBodyJointOptions)

public:
	DialogBodyJointOptions(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~DialogBodyJointOptions();

// Dialogfelddaten
	enum { IDD = IDD_DIALOGBODYJOINTOPTIONS };

//functions not added by windows:
	void Create(CWnd * pParent);
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	void LoadData();  //Get data from WCDinterface
	void WriteData(); //Put data to WCDinterface

private:
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	void OnOK();
	void OnCancel();
	void OnClose();

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedApply();
	afx_msg void OnBnClickedOk();
	BOOL m_check_show_joints;
	afx_msg void OnBnClickedCancel();
	BOOL m_show_body_numbers;
	BOOL m_check_usedegrees;
	BOOL m_check_show_constraint_numbers;
	BOOL m_check_showbodylocalframe;
	double m_bodylocalframesize;
	BOOL m_check_showsensors;
	double m_sensor_size;
	BOOL m_check_joints_transparent;
	BOOL m_check_sensors_transparent;
	BOOL m_check_bodies_transparent;
	BOOL m_check_bodies_supersmooth;
	BOOL m_check_show_body_outline;
	BOOL m_check_show_body_faces;
	int m_radio_eulerangles;
	BOOL m_check_showloads;
	double m_load_draw_size;
	afx_msg void OnBnClickedCheckShowJoints2();
	//BOOL m_check_show_control_ojects;
	afx_msg void OnEnChangeEditLoaddrawsize();
};
