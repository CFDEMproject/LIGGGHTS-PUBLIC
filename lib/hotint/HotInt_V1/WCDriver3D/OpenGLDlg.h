//#**************************************************************
//# filename:             OpenGLDlg.h
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


// DialogOpenGLOptions-Dialogfeld

class DialogOpenGLOptions : public CDialog
{
	DECLARE_DYNAMIC(DialogOpenGLOptions)

public:
	DialogOpenGLOptions(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~DialogOpenGLOptions();

// Dialogfelddaten
	enum { IDD = IDD_OPENGL_OPTIONSDIALOG };

//functions not added by windows:
	void Create(CWnd * pParent);
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	void LoadData();  //Get data from WCDinterface
	void WriteData(); //Put data to WCDinterface
	void UpdateIfActivated(); //update data and redraw if activated by check

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
	BOOL m_check_immediate_apply;
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedApply();
	afx_msg void OnBnClickedCancel();
	BOOL m_check_smooth_shade_model;
protected:
	virtual BOOL OnNotify(WPARAM wParam, LPARAM lParam, LRESULT* pResult);
public:
	CSliderCtrl sliderctrl_transparency;
	afx_msg void OnNMCustomdrawSliderTransparency(NMHDR *pNMHDR, LRESULT *pResult);
	int m_transparency;
	BOOL m_check_enable_lighting;
	CSliderCtrl sliderctrl_shininess;
	int m_shininess;
	CSliderCtrl sliderctrl_speccolor;
	int m_speccolor;
	afx_msg void OnNMCustomdrawSliderShininess(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderSpeccolor(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl sliderctrl_ambient_light;
	int m_ambient_light;
	CSliderCtrl sliderctrl_diffuse_light;
	int m_diffuse_light;
	CSliderCtrl sliderctrl_specular_light;
	int m_specular_light;
	afx_msg void OnNMCustomdrawSliderAmbientLight(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderDiffuseLight(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderSpecularLight(NMHDR *pNMHDR, LRESULT *pResult);
	double m_lightposz;
	double m_lightposy;
	double m_lightposx;
	BOOL m_check_incllightpos;
	BOOL m_check_enable_light;
//	afx_msg void OnBnKillfocusCheckEnableLight();
//	afx_msg void OnBnKillfocusCheckIncllightpos();
	afx_msg void OnEnKillfocusEditLightposx();
	afx_msg void OnEnKillfocusEditLightposy();
	afx_msg void OnEnKillfocusEditLightposz();
	afx_msg void OnBnClickedCheckIncllightpos();
	afx_msg void OnBnClickedCheckEnableLight();
	afx_msg void OnBnClickedCheck1();
	BOOL m_check_enable_light2;
	BOOL m_check_incllightpos2;
	CSliderCtrl sliderctrl_ambient_light2;
	int m_ambient_light2;
	CSliderCtrl sliderctrl_diffuse_light2;
	int m_diffuse_light2;
	CSliderCtrl sliderctrl_specular_light2;
	int m_specular_light2;
	double m_lightposx2;
	double m_lightposy2;
	double m_lightposz2;
	afx_msg void OnBnClickedCheckEnableLight2();
	afx_msg void OnBnClickedCheckIncllightpos2();
	afx_msg void OnNMCustomdrawSliderAmbientLight2(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderDiffuseLight2(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderSpecularLight2(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnEnKillfocusEditLightposx2();
	afx_msg void OnEnKillfocusEditLightposy2();
	afx_msg void OnEnKillfocusEditLightposz2();
	afx_msg void OnBnClickedCheckSmoothshademodel();
	afx_msg void OnBnClickedCheckEnablelighting();
};
