//#**************************************************************
//# filename:             GeneralOptions.cpp
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
 


#include "stdafx.h"
#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"
#include "FEdrawingOptions.h"
#include "fieldvariabledescriptor.h"


// GeneralOptions-Dialogfeld

//IMPLEMENT_DYNAMIC(FEDrawingOptions, CDialog)
FEDrawingOptions::FEDrawingOptions(/*CWnd* pParent*/)
	: CDialog(FEDrawingOptions::IDD, NULL)
	, m_maxstress(0), m_minstress(0)
	, m_check_maxstress(FALSE)
	, m_check_minstress(FALSE)
	, m_colortiling(0)
	, m_check_greymode(FALSE)
	, m_check_invertcolors(FALSE)
	, m_check_nonlinearscale(FALSE)
	, m_check_hidelegend(FALSE)
	, m_check_showmesh(TRUE)
	, m_check_showsolution(FALSE)
	, m_check_chladni_isolines(FALSE)
	, m_check_draw_flat_elements(FALSE)
	, m_elem_line_thickness(1)
	, m_shrinking_factor(0)
	, m_axis_tiling(0)
	, m_axis_tiling_text(0)
	, m_axis_resolution(0)
	, m_axis_resolution_text(0)
	, m_crosssection_res(0)
	, m_crosssection_res_text(0)
	, m_solidFE_resolution(0)
	, m_solidFE_resolution_text(0)
	, m_colortilingtext(_T(""))
	, m_deformation_scale_factor(0)
	, m_check_plot_interpolated(FALSE)
	, m_combo_units(0)
	, m_check_shownodes(FALSE)
	, m_node_draw_size(0)
	, m_check_shownodenumbers(FALSE)
	, m_check_surface_elements_only(FALSE)
	, m_check_animate_scaling_factor(FALSE)
	, m_check_scale_rigid_body_disp(FALSE)
	, m_check_auto_adjust_range(FALSE)
{
}


void FEDrawingOptions::Create(CWnd * pParent)
{
	LoadData();

	CDialog::Create(IDD,pParent);

	// we populate the combo box with the types for the post-processing
	const TArray<FieldVariableDescriptor> & variables = pWCDI->GetAvailableFieldVariables();
	m_comboFieldVariableTypes.AddString("body color");
	for(int i = 1; i <= variables.Length(); i++)
	{
		const FieldVariableDescriptor & fvd = variables(i);
		if(!fvd.IsNotForPlotting())
			if(m_comboFieldVariableTypes.FindString(0, fvd.GetTextualIdentifierWithoutComponents()) == CB_ERR)
				m_comboFieldVariableTypes.AddString(fvd.GetTextualIdentifierWithoutComponents());
	}

	control_combo_units.Clear();
	control_combo_units.AddString("N, m");
	control_combo_units.AddString("N, mm");

	control_slider_colortiling.SetRange(1,33); //needs to be done before LoadData!
	slider_axis_tiling.SetRange(1,32);
	slider_axis_resolution.SetRange(1,32);
	slider_crosssection_res.SetRange(1,8);
	slider_solidFE_resolution.SetRange(1,8);

	slider_axis_tiling.SetTicFreq(2); //draw less ticks than 32, otherwise it looks very dense
	slider_axis_resolution.SetTicFreq(2);

	//these commands help that the sliders are correctly refreshed when dialog is opened the second time
	slider_axis_tiling.SetPos(m_axis_tiling+100);
	slider_axis_resolution.SetPos(m_axis_resolution+100);
	slider_crosssection_res.SetPos(m_crosssection_res+100);
	slider_solidFE_resolution.SetPos(m_solidFE_resolution+100);
	control_slider_colortiling.SetPos(m_colortiling+100);

	if (pWCDI->GetIOption(103) != 33)
	{
		char str[20];
		sprintf(str, "%d", pWCDI->GetIOption(103));
		m_colortilingtext = CString(str);
	}
	else
		m_colortilingtext = CString("no color tiling");

	UpdateData(FALSE);

	CRect r;
	pParent->GetWindowRect(&r);
	CPoint p(r.right - 400,r.bottom - 0);
	GetWindowRect(&r);
//	ClientToScreen(&r);
	r.top = p.y - r.Height();
	r.left = p.x - r.Width();
	r.right = p.x;
	r.bottom = p.y;
	MoveWindow(r,FALSE);
	ShowWindow(SW_SHOW);

	RedrawWindow();

}

void FEDrawingOptions::LoadData()
{
	//+++++++++++++++++++++++++++++++++++++++++++++++
	//graphics/postprocessing:
	m_check_maxstress          = pWCDI->GetIOption(100);
	m_check_minstress          = pWCDI->GetIOption(101);
	m_maxstress                = pWCDI->GetDOption(100);
	m_minstress                = pWCDI->GetDOption(101);
	m_current_field_variable_index = pWCDI->GetIndexOfActualPostProcessingFieldVariable();
	m_colortiling              = pWCDI->GetIOption(103);
	m_check_greymode           = pWCDI->GetIOption(105);
	m_check_invertcolors       = pWCDI->GetIOption(107);
	m_check_nonlinearscale     = pWCDI->GetIOption(108);
	m_check_hidelegend         = pWCDI->GetIOption(141);
	m_check_auto_adjust_range	 = pWCDI->GetIOption(102);
												 
	m_check_showmesh           = pWCDI->GetIOption(110);
	m_check_showsolution       = pWCDI->GetIOption(111);
	m_check_chladni_isolines   = pWCDI->GetIOption(116);
	m_check_plot_interpolated  = pWCDI->GetIOption(118);

	m_check_draw_flat_elements = pWCDI->GetIOption(117);
	m_elem_line_thickness		   = pWCDI->GetDOption(102);
	m_shrinking_factor		     = pWCDI->GetDOption(106);
	m_deformation_scale_factor = pWCDI->GetDOption(105);

	//tiling:
	m_axis_tiling						= pWCDI->GetIOption(136)/2;
	m_axis_resolution				= pWCDI->GetIOption(137);
	m_crosssection_res			= pWCDI->GetIOption(138);
	m_solidFE_resolution		= pWCDI->GetIOption(139);

	m_combo_units						= pWCDI->GetIOption(143);
	
	m_check_shownodes       = pWCDI->GetIOption(145);
	m_node_draw_size				= pWCDI->GetDOption(124);
	m_check_shownodenumbers = pWCDI->GetIOption(148);
	m_check_surface_elements_only = pWCDI->GetIOption(146);

	m_check_animate_scaling_factor = pWCDI->GetIOption(150);
	m_check_scale_rigid_body_disp = pWCDI->GetIOption(151);
}

void FEDrawingOptions::WriteData()
{
	//graphics/postprocessing: 
	//integer/bool:
	pWCDI->SetIOption(100, (int)m_check_maxstress);
	pWCDI->SetIOption(101, (int)m_check_minstress);
	pWCDI->SetIOption(103, m_colortiling);
	pWCDI->SetIOption(105, (int)m_check_greymode);
	pWCDI->SetIOption(107, (int)m_check_invertcolors);
	pWCDI->SetIOption(108, (int)m_check_nonlinearscale);
	pWCDI->SetIOption(141, (int)m_check_hidelegend);
	pWCDI->SetIOption(102, (int)m_check_auto_adjust_range);

	pWCDI->SetIOption(110, m_check_showmesh       );
	pWCDI->SetIOption(111, m_check_showsolution   );
	pWCDI->SetIOption(116, m_check_chladni_isolines);
	pWCDI->SetIOption(117, m_check_draw_flat_elements);
	pWCDI->SetIOption(118, m_check_plot_interpolated);
	pWCDI->SetIndexOfActualPostProcessingFieldVariable(m_current_field_variable_index);

	// double values:
	pWCDI->SetDOption(100, m_maxstress);
	pWCDI->SetDOption(101, m_minstress);
	pWCDI->SetDOption(102, m_elem_line_thickness);
	pWCDI->SetDOption(106, m_shrinking_factor);
	pWCDI->SetDOption(105, m_deformation_scale_factor);

	//tiling:
	pWCDI->GetIOption(136)	= 2*m_axis_tiling					;
	pWCDI->GetIOption(137)	= m_axis_resolution			;
	pWCDI->GetIOption(138)	= m_crosssection_res		;
	pWCDI->GetIOption(139)	= m_solidFE_resolution	;

	pWCDI->GetIOption(143) = 	m_combo_units;
	pWCDI->GetIOption(145) =  m_check_shownodes;
	pWCDI->GetIOption(147) =  3;
	pWCDI->GetIOption(148) = m_check_shownodenumbers;
	pWCDI->GetDOption(124) =  m_node_draw_size;

	pWCDI->GetIOption(146) = m_check_surface_elements_only;

	pWCDI->GetIOption(150) = m_check_animate_scaling_factor;
	pWCDI->GetIOption(151) = m_check_scale_rigid_body_disp;

	pGLDrawWnd->SetAnimateScalingTimer(m_check_animate_scaling_factor);
}

void FEDrawingOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX, IDC_MAXSTRESS, m_maxstress);
	DDX_Text(pDX, IDC_MINSTRESS, m_minstress);
	DDX_Check(pDX, IDC_CHECK_MAXSTRESS, m_check_maxstress);
	DDX_Check(pDX, IDC_CHECK_MINSTRESS, m_check_minstress);
	DDX_Slider(pDX, IDC_SLIDER_COLORTILING, m_colortiling);
	DDX_Control(pDX, IDC_SLIDER_COLORTILING, control_slider_colortiling);
	DDX_Check(pDX, IDC_CHECK_GREYMODE, m_check_greymode);
	DDX_Check(pDX, IDC_CHECK_INVERTCOLORS, m_check_invertcolors);
	DDX_Check(pDX, IDC_CHECK_NONLINEARSCALE, m_check_nonlinearscale);
	DDX_Check(pDX, IDC_CHECK_HIDELEGEND, m_check_hidelegend);
	DDX_Check(pDX, IDC_CHECK_SHOWMESH, m_check_showmesh);
	DDX_Check(pDX, IDC_SHOWSOLUTION, m_check_showsolution);
	DDX_Check(pDX, IDC_CHECK_CHLADNI_ISOLINES, m_check_chladni_isolines);
	DDX_Check(pDX, IDC_DRAW_FLAT_ELEMENTS, m_check_draw_flat_elements);
	DDX_Text(pDX, IDC_ELEM_LINE_THICKNESS, m_elem_line_thickness);
	DDV_MinMaxDouble(pDX, m_elem_line_thickness, 0.01, 100);
	DDX_Text(pDX, IDC_EDIT_SHRINKING_FACTOR, m_shrinking_factor);
	DDX_Slider(pDX, IDC_SLIDER_AXISTILING, m_axis_tiling);
	DDX_Text(pDX, IDC_STATIC_AXISTILING, m_axis_tiling_text);
	DDX_Control(pDX, IDC_SLIDER_AXISTILING, slider_axis_tiling);
	DDX_Control(pDX, IDC_SLIDER_AXISRESOLUTION, slider_axis_resolution);
	DDX_Slider(pDX, IDC_SLIDER_AXISRESOLUTION, m_axis_resolution);
	DDX_Text(pDX, IDC_STATIC_AXISRESOLUTION, m_axis_resolution_text);
	DDX_Control(pDX, IDC_SLIDER_CROSSSECTION_RES, slider_crosssection_res);
	DDX_Slider(pDX, IDC_SLIDER_CROSSSECTION_RES, m_crosssection_res);
	DDX_Text(pDX, IDC_STATIC_CROSSSECTION_RES, m_crosssection_res_text);
	DDX_Control(pDX, IDC_SLIDER_FERESOLUTION, slider_solidFE_resolution);
	DDX_Slider(pDX, IDC_SLIDER_FERESOLUTION, m_solidFE_resolution);
	DDX_Text(pDX, IDC_STATIC_FERESOLUTION, m_solidFE_resolution_text);
	DDX_Text(pDX, IDC_STATIC_COLORTILINGTEXT, m_colortilingtext);

	DDX_Text(pDX, IDC_EDIT_DEFORMATION_SCALE_FACTOR, m_deformation_scale_factor);
	//DDV_MinMaxDouble(pDX, m_deformation_scale_factor, 0, 1e30);
	DDX_Check(pDX, IDC_PLOT_INTERPOLATED, m_check_plot_interpolated);
	DDX_CBIndex(pDX, IDC_COMBO_UNITS, m_combo_units);
	DDX_Control(pDX, IDC_COMBO_UNITS, control_combo_units);
	DDX_Check(pDX, IDC_CHECK_SHOWNODES, m_check_shownodes);
	DDX_Text(pDX, IDC_NODE_SIZE, m_node_draw_size);
	DDV_MinMaxDouble(pDX, m_node_draw_size, 0, 1e30);
	DDX_Check(pDX, IDC_CHECK_SHOWNODES_NUMBERS, m_check_shownodenumbers);
	DDX_Check(pDX, IDC_PLOT_SURFACEONLY, m_check_surface_elements_only);
	DDX_Check(pDX, IDC_CHECK_ANIMATE_SCALING_FACTOR, m_check_animate_scaling_factor);
	DDX_Check(pDX, IDC_CHECK_SCALE_RIGID_BODY_DISP, m_check_scale_rigid_body_disp);
	DDX_Control(pDX, IDC_COMBO_FIELD_VARIABLE_TYPE, m_comboFieldVariableTypes);
	DDX_Control(pDX, IDC_COMBO_FIELD_VARIABLE_COMPONENTS, m_comboFieldVariableComponents);

	DDX_Check(pDX, IDC_CHECK_AUTO_ADJUST_RANGE, m_check_auto_adjust_range);
	GetDlgItem(IDC_BUTTON_ADJUST_RANGE)->EnableWindow(!m_check_auto_adjust_range);

	// special treatment of the selector of the field variable:
	// the result is determined by the state of the two combo-boxes
	const TArray<FieldVariableDescriptor> & variables = pWCDI->GetAvailableFieldVariables();
	if(pDX->m_bSaveAndValidate)
	{
		// reading the data from the dialog
		if(m_comboFieldVariableTypes.GetCurSel() == 0)
			m_current_field_variable_index = 0;		// body color
		else
		{
			CString textType;
			m_comboFieldVariableTypes.GetLBText(m_comboFieldVariableTypes.GetCurSel(), textType);
			CString textComponents;
			if(m_comboFieldVariableComponents.GetCurSel() != CB_ERR)
				m_comboFieldVariableComponents.GetLBText(m_comboFieldVariableComponents.GetCurSel(), textComponents);
			int firstTypeMatchingIndex = 0;
			for(int i = 1; i <= variables.Length(); i++)
			{
				const FieldVariableDescriptor & fvd = variables(i);
				if(textType.Compare(fvd.GetTextualIdentifierWithoutComponents()) == 0)
				{
					// the type is the same - let's check the components
					if(firstTypeMatchingIndex == 0)
						firstTypeMatchingIndex = i;
					if(textComponents.Compare(fvd.GetTextualIdentifierComponentsOnly()) == 0)
					{
						// found a matching entry
						m_current_field_variable_index = i;
						firstTypeMatchingIndex = 0;
						break;
					}
				}
			}
			// did not find an exact matching entry - let's select the first entry with the matching type;
			// this behavior is suitable for the selection change in m_comboFieldVariableTypes
			if(firstTypeMatchingIndex != 0)
				m_current_field_variable_index = firstTypeMatchingIndex;
		}
	}
	else
	{
		// writing the data to the dialog
		m_comboFieldVariableComponents.ResetContent();
		m_comboFieldVariableComponents.EnableWindow(FALSE);
		if(m_current_field_variable_index == 0 || m_current_field_variable_index > variables.Length())
			m_comboFieldVariableTypes.SetCurSel(0);		// body color
		else
		{
			const FieldVariableDescriptor & fvd = variables(m_current_field_variable_index);
			CString textType = fvd.GetTextualIdentifierWithoutComponents();
			CString textComponents = fvd.GetTextualIdentifierComponentsOnly();
			m_comboFieldVariableTypes.SelectString(0, textType);
			// the combo with the components needs to be populated and set to the right position
			// we browse thru all variables with the same type name
			for(int i = 1; i <= variables.Length(); i++)
			{
				const FieldVariableDescriptor & fvd1 = variables(i);
				if(*fvd1.GetTextualIdentifierComponentsOnly() != '\0')		// there are some components
					if(textType.Compare(fvd1.GetTextualIdentifierWithoutComponents()) == 0)		// and the type is the same
					{
						m_comboFieldVariableComponents.EnableWindow(TRUE);
						m_comboFieldVariableComponents.AddString(fvd1.GetTextualIdentifierComponentsOnly());
					}
			}
			// and finally we select the corresponding string in the combo with components
			m_comboFieldVariableComponents.SelectString(-1, fvd.GetTextualIdentifierComponentsOnly());
		}
	}
}


BEGIN_MESSAGE_MAP(FEDrawingOptions, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_EN_KILLFOCUS(IDC_MAXSTRESS, OnEnKillfocusMaxstress)
	ON_EN_KILLFOCUS(IDC_MINSTRESS, OnEnKillfocusMinstress)
	ON_BN_CLICKED(IDC_CHECK_MINSTRESS, OnBnClickedCheckMinstress)
	ON_BN_CLICKED(IDC_CHECK_MAXSTRESS, OnBnClickedCheckMaxstress)
	ON_BN_CLICKED(IDC_BUTTON_UPDATEOPTIONS, OnBnClickedButtonUpdateoptions)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_COLORTILING, OnNMCustomdrawSliderColortiling)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_AXISTILING, OnNMCustomdrawSliderAxistiling)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_AXISRESOLUTION, OnNMCustomdrawSliderAxisresolution)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_CROSSSECTION_RES, OnNMCustomdrawSliderCrosssectionRes)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_FERESOLUTION, OnNMCustomdrawSliderFeresolution)
	ON_CBN_SELCHANGE(IDC_COMBO_FIELD_VARIABLE_TYPE, &FEDrawingOptions::OnCbnSelchangeComboFieldVariableType)
	ON_BN_CLICKED(IDC_CHECK_AUTO_ADJUST_RANGE, &FEDrawingOptions::OnBnClickedCheckAutoAdjustRange)
	ON_BN_CLICKED(IDC_BUTTON_ADJUST_RANGE, &FEDrawingOptions::OnBnClickedButtonAdjustRange)
END_MESSAGE_MAP()


// FEDrawingOptions-Meldungshandler

void FEDrawingOptions::OnBnClickedOk()
{
	OnBnClickedButtonUpdateoptions();
	
	OnOK();
}

void FEDrawingOptions::OnOK()
{
	DestroyWindow();
}

void FEDrawingOptions::OnClose() 
{
	DestroyWindow();
}

void FEDrawingOptions::OnCancel() 
{
	DestroyWindow();
}

void FEDrawingOptions::OnEnKillfocusMaxstress()
{

}

void FEDrawingOptions::OnEnKillfocusMinstress()
{
}

void FEDrawingOptions::OnBnClickedCheckMaxstress()
{
}

void FEDrawingOptions::OnBnClickedCheckMinstress()
{
}


void FEDrawingOptions::OnBnClickedButtonUpdateoptions()
{
	UpdateData(TRUE);
	WriteData();

	pGLDrawWnd->Redraw();
}

void FEDrawingOptions::OnNMCustomdrawSliderColortiling(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	if (m_colortiling != 33)
	{
		char str[20];
		sprintf(str, "%d", m_colortiling);
		m_colortilingtext = CString(str);
	}
	else
		m_colortilingtext = CString("no color tiling");

	UpdateData(FALSE);

	*pResult = 0;
}


void FEDrawingOptions::OnNMCustomdrawSliderAxistiling(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	m_axis_tiling_text = 2*m_axis_tiling;
	UpdateData(FALSE);

	*pResult = 0;

}


void FEDrawingOptions::OnNMCustomdrawSliderAxisresolution(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	m_axis_resolution_text = m_axis_resolution;
	UpdateData(FALSE);

	*pResult = 0;
}

void FEDrawingOptions::OnNMCustomdrawSliderCrosssectionRes(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	m_crosssection_res_text = m_crosssection_res;
	UpdateData(FALSE);

	*pResult = 0;
}

void FEDrawingOptions::OnNMCustomdrawSliderFeresolution(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	m_solidFE_resolution_text = m_solidFE_resolution;
	UpdateData(FALSE);

	*pResult = 0;
}


void FEDrawingOptions::OnCbnSelchangeComboFieldVariableType()
{
	UpdateData(TRUE);
	UpdateData(FALSE);
}

void FEDrawingOptions::OnBnClickedCheckAutoAdjustRange()
{
	UpdateData(TRUE);
	GetDlgItem(IDC_BUTTON_ADJUST_RANGE)->EnableWindow(!m_check_auto_adjust_range);
	m_check_auto_adjust_range	 = pWCDI->GetIOption(102);
}

void FEDrawingOptions::OnBnClickedButtonAdjustRange()
{
	pWCDI->SetIOption(102, 1);
	pGLDrawWnd->Redraw();
	pWCDI->SetIOption(102, 0);
}
