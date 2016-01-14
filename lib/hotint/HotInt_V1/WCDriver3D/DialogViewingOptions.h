//#**************************************************************
//# filename:             DialogViewingOptions.h
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
#include "afxcmn.h"
#include "afxwin.h"


// DialogViewingOptions-Dialogfeld

class DialogViewingOptions : public CDialog
{
	DECLARE_DYNAMIC(DialogViewingOptions)

public:
	DialogViewingOptions(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~DialogViewingOptions();

// Dialogfelddaten
	enum { IDD = IDD_DIALOGVIEWINGOPTIONS };

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
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedApply();
	afx_msg void OnBnClickedOk();
	afx_msg void OnNMCustomdrawSliderRedrawfrequency(NMHDR *pNMHDR, LRESULT *pResult);
	int m_redrawfrequency;
	BOOL m_check_draworigin;
	//BOOL m_check_showcontactpoints;
	double m_standardview_angle1;
	double m_standardview_angle2;
	double m_standardview_angle3;
	int m_standardview_axis1;
	int m_standardview_axis2;
	double m_standardview_axis3;
	int m_WindowSizeX;
	int m_WindowSizeY;
	CSliderCtrl slider_redrawfrequency;
	afx_msg void OnBnClickedCheckDraworigin();
	afx_msg void OnNMReleasedcaptureSliderRedrawfrequency(NMHDR *pNMHDR, LRESULT *pResult);
	BOOL m_animate_beginning;
	int m_animationframes;
	CSliderCtrl slider_animationframes;
	afx_msg void OnNMCustomdrawSliderAnimationframes(NMHDR *pNMHDR, LRESULT *pResult);
	BOOL m_check_textstofront;
	double m_origin_size;
	double m_gridsize1;
	double m_gridstep1;
	CComboBox combo_gridtype;
	int m_combo_gridtype;
	double m_grid_posx;
	double m_grid_posy;
	double m_grid_posz;
	//BOOL m_check_show_startup_banner;
	BOOL m_check_use_cuttingplane;
	double m_cutplane_distance;
	double m_cutplane_nx;
	double m_cutplane_ny;
	double m_cutplane_nz;
};
