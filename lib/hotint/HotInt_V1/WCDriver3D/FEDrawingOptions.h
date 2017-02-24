//#**************************************************************
//# filename:             FEDrawingOptions.h
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
 

#if !defined(GENERALOPTIONS__INCLUDED_)
#define GENERALOPTIONS__INCLUDED_


#pragma once

#include "gldrawwnd.h"
#include "afxwin.h"
#include "afxcmn.h"

// FEDrawingOptions-Dialogfeld

class FEDrawingOptions : public CDialog
{
	//DECLARE_DYNAMIC(FEDrawingOptions)

public:
	FEDrawingOptions();   // Standardkonstruktor
	//virtual ~FEDrawingOptions();

// Dialogfelddaten
	enum { IDD = IDD_DIALOG_GENERALOPTIONS };


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
	//should not be used by windows, data destroyed otherwise ...
	void OnOK();
	void OnCancel();
	void OnClose();

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnEnKillfocusMaxstress();
	double m_maxstress;
	afx_msg void OnEnKillfocusMinstress();
	double m_minstress;
	afx_msg void OnBnClickedCheckMinstress();
	afx_msg void OnBnClickedCheckMaxstress();
	BOOL m_check_maxstress;
	BOOL m_check_minstress;
	afx_msg void OnBnClickedButtonUpdateoptions();
	int m_colortiling;
//	afx_msg void OnNMReleasedcaptureSliderColortiling(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderColortiling(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl control_slider_colortiling;
	BOOL m_check_greymode;
	BOOL m_check_invertcolors;
	BOOL m_check_nonlinearscale;
	BOOL m_check_hidelegend;
	BOOL m_check_showmesh;
	BOOL m_check_showsolution;
	BOOL m_check_chladni_isolines;
	BOOL m_check_draw_flat_elements;
	double m_elem_line_thickness;
	double m_shrinking_factor;
	int m_axis_tiling;
	int m_axis_tiling_text;
	afx_msg void OnNMCustomdrawSliderAxistiling(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl slider_axis_tiling;
	CSliderCtrl slider_axis_resolution;
	int m_axis_resolution;
	int m_axis_resolution_text;
	afx_msg void OnNMCustomdrawSliderAxisresolution(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl slider_crosssection_res;
	int m_crosssection_res;
	int m_crosssection_res_text;
	afx_msg void OnNMCustomdrawSliderCrosssectionRes(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl slider_solidFE_resolution;
	int m_solidFE_resolution;
	int m_solidFE_resolution_text;
	afx_msg void OnNMCustomdrawSliderFeresolution(NMHDR *pNMHDR, LRESULT *pResult);
	CString m_colortilingtext;
	double m_deformation_scale_factor;
	BOOL m_check_plot_interpolated;
	int m_combo_units;
	CComboBox control_combo_units;
	BOOL m_check_shownodes;
	double m_node_draw_size;
	BOOL m_check_shownodenumbers;
	BOOL m_check_surface_elements_only;
public:
	BOOL m_check_animate_scaling_factor;
public:
	BOOL m_check_scale_rigid_body_disp;
	CComboBox m_comboFieldVariableTypes;
	CComboBox m_comboFieldVariableComponents;
	int m_current_field_variable_index;
	afx_msg void OnCbnSelchangeComboFieldVariableType();
	afx_msg void OnBnClickedCheckAutoAdjustRange();
	afx_msg void OnBnClickedButtonAdjustRange();
	BOOL m_check_auto_adjust_range;
};

#endif
