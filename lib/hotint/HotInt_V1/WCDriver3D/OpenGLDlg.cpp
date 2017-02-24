//#**************************************************************
//# filename:             OpenGLDlg.cpp
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
#include "OpenGLDlg.h"


// DialogOpenGLOptions-Dialogfeld

IMPLEMENT_DYNAMIC(DialogOpenGLOptions, CDialog)
DialogOpenGLOptions::DialogOpenGLOptions(CWnd* pParent /*=NULL*/)
	: CDialog(DialogOpenGLOptions::IDD, pParent)
	, m_check_immediate_apply(FALSE)
	, m_check_smooth_shade_model(FALSE)
	, m_transparency(0)
	, m_check_enable_lighting(FALSE)
	, m_shininess(0)
	, m_speccolor(0)
	, m_ambient_light(0)
	, m_diffuse_light(0)
	, m_specular_light(0)
	, m_lightposz(0)
	, m_lightposy(0)
	, m_lightposx(0)
	, m_check_incllightpos(FALSE)
	, m_check_enable_light(FALSE)
	, m_check_enable_light2(FALSE)
	, m_check_incllightpos2(FALSE)
	, m_ambient_light2(0)
	, m_diffuse_light2(0)
	, m_specular_light2(0)
	, m_lightposx2(0)
	, m_lightposy2(0)
	, m_lightposz2(0)
{
}

DialogOpenGLOptions::~DialogOpenGLOptions()
{
}

void DialogOpenGLOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK1, m_check_immediate_apply);
	DDX_Check(pDX, IDC_CHECK_SMOOTHSHADEMODEL, m_check_smooth_shade_model);

	DDX_Control(pDX, IDC_SLIDER_TRANSPARENCY, sliderctrl_transparency);
	DDX_Slider(pDX, IDC_SLIDER_TRANSPARENCY, m_transparency);
	DDX_Check(pDX, IDC_CHECK_ENABLELIGHTING, m_check_enable_lighting);
	DDX_Control(pDX, IDC_SLIDER_SHININESS, sliderctrl_shininess);
	DDX_Slider(pDX, IDC_SLIDER_SHININESS, m_shininess);
	DDX_Control(pDX, IDC_SLIDER_SPECCOLOR, sliderctrl_speccolor);
	DDX_Slider(pDX, IDC_SLIDER_SPECCOLOR, m_speccolor);
	DDX_Control(pDX, IDC_SLIDER_AMBIENT_LIGHT, sliderctrl_ambient_light);
	DDX_Slider(pDX, IDC_SLIDER_AMBIENT_LIGHT, m_ambient_light);
	DDX_Control(pDX, IDC_SLIDER_DIFFUSE_LIGHT, sliderctrl_diffuse_light);
	DDX_Slider(pDX, IDC_SLIDER_DIFFUSE_LIGHT, m_diffuse_light);
	DDX_Control(pDX, IDC_SLIDER_SPECULAR_LIGHT, sliderctrl_specular_light);
	DDX_Slider(pDX, IDC_SLIDER_SPECULAR_LIGHT, m_specular_light);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSZ, m_lightposz);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSY, m_lightposy);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSX, m_lightposx);
	DDX_Check(pDX, IDC_CHECK_INCLLIGHTPOS, m_check_incllightpos);
	DDX_Check(pDX, IDC_CHECK_ENABLE_LIGHT, m_check_enable_light);
	DDX_Check(pDX, IDC_CHECK_ENABLE_LIGHT2, m_check_enable_light2);
	DDX_Check(pDX, IDC_CHECK_INCLLIGHTPOS2, m_check_incllightpos2);
	DDX_Control(pDX, IDC_SLIDER_AMBIENT_LIGHT2, sliderctrl_ambient_light2);
	DDX_Slider(pDX, IDC_SLIDER_AMBIENT_LIGHT2, m_ambient_light2);
	DDX_Control(pDX, IDC_SLIDER_DIFFUSE_LIGHT2, sliderctrl_diffuse_light2);
	DDX_Slider(pDX, IDC_SLIDER_DIFFUSE_LIGHT2, m_diffuse_light2);
	DDX_Control(pDX, IDC_SLIDER_SPECULAR_LIGHT2, sliderctrl_specular_light2);
	DDX_Slider(pDX, IDC_SLIDER_SPECULAR_LIGHT2, m_specular_light2);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSX2, m_lightposx2);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSY2, m_lightposy2);
	DDX_Text(pDX, IDC_EDIT_LIGHTPOSZ2, m_lightposz2);
}


BEGIN_MESSAGE_MAP(DialogOpenGLOptions, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDAPPLY, OnBnClickedApply)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
	ON_WM_KILLFOCUS()
	ON_WM_LBUTTONUP()
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_TRANSPARENCY, OnNMCustomdrawSliderTransparency)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_SHININESS, OnNMCustomdrawSliderShininess)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_SPECCOLOR, OnNMCustomdrawSliderSpeccolor)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_AMBIENT_LIGHT, OnNMCustomdrawSliderAmbientLight)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_DIFFUSE_LIGHT, OnNMCustomdrawSliderDiffuseLight)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_SPECULAR_LIGHT, OnNMCustomdrawSliderSpecularLight)
//	ON_BN_KILLFOCUS(IDC_CHECK_ENABLE_LIGHT, OnBnKillfocusCheckEnableLight)
//	ON_BN_KILLFOCUS(IDC_CHECK_INCLLIGHTPOS, OnBnKillfocusCheckIncllightpos)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSX, OnEnKillfocusEditLightposx)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSY, OnEnKillfocusEditLightposy)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSZ, OnEnKillfocusEditLightposz)
	ON_BN_CLICKED(IDC_CHECK_INCLLIGHTPOS, OnBnClickedCheckIncllightpos)
	ON_BN_CLICKED(IDC_CHECK_ENABLE_LIGHT, OnBnClickedCheckEnableLight)
	ON_BN_CLICKED(IDC_CHECK1, OnBnClickedCheck1)
	ON_BN_CLICKED(IDC_CHECK_ENABLE_LIGHT2, OnBnClickedCheckEnableLight2)
	ON_BN_CLICKED(IDC_CHECK_INCLLIGHTPOS2, OnBnClickedCheckIncllightpos2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_AMBIENT_LIGHT2, OnNMCustomdrawSliderAmbientLight2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_DIFFUSE_LIGHT2, OnNMCustomdrawSliderDiffuseLight2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_SPECULAR_LIGHT2, OnNMCustomdrawSliderSpecularLight2)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSX2, OnEnKillfocusEditLightposx2)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSY2, OnEnKillfocusEditLightposy2)
	ON_EN_KILLFOCUS(IDC_EDIT_LIGHTPOSZ2, OnEnKillfocusEditLightposz2)
	ON_BN_CLICKED(IDC_CHECK_SMOOTHSHADEMODEL, OnBnClickedCheckSmoothshademodel)
	ON_BN_CLICKED(IDC_CHECK_ENABLELIGHTING, OnBnClickedCheckEnablelighting)
END_MESSAGE_MAP()


// DialogOpenGLOptions-Meldungshandler

void DialogOpenGLOptions::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	LoadData();

	char str[50];

	sliderctrl_transparency.SetRange(0,100);
	sprintf(str,"transparency: %d%%", 100-m_transparency);
	GetDlgItem(IDC_STATIC_TRANSPARENCY)->SetWindowText(str);

	sliderctrl_shininess.SetRange(0,100);
	sprintf(str,"shininess: %d%%", m_shininess);
	GetDlgItem(IDC_STATIC_SHININESS)->SetWindowText(str);

	sliderctrl_speccolor.SetRange(0,100);
	sprintf(str,"specular color intensity: %d%%", m_speccolor);
	GetDlgItem(IDC_STATIC_SPECCOLOR)->SetWindowText(str);

	//light1:
	sliderctrl_ambient_light.SetRange(0,100);
	sprintf(str,"ambient: %d%%", m_ambient_light);
	GetDlgItem(IDC_STATIC_AMBIENT_LIGHT)->SetWindowText(str);

	sliderctrl_diffuse_light.SetRange(0,100);
	sprintf(str,"diffuse: %d%%", m_diffuse_light);
	GetDlgItem(IDC_STATIC_DIFFUSE_LIGHT)->SetWindowText(str);

	sliderctrl_specular_light.SetRange(0,100);
	sprintf(str,"specular: %d%%", m_specular_light);
	GetDlgItem(IDC_STATIC_SPECULAR_LIGHT)->SetWindowText(str);

	//light2:
	sliderctrl_ambient_light2.SetRange(0,100);
	sprintf(str,"ambient: %d%%", m_ambient_light2);
	GetDlgItem(IDC_STATIC_AMBIENT_LIGHT2)->SetWindowText(str);

	sliderctrl_diffuse_light2.SetRange(0,100);
	sprintf(str,"diffuse: %d%%", m_diffuse_light2);
	GetDlgItem(IDC_STATIC_DIFFUSE_LIGHT2)->SetWindowText(str);

	sliderctrl_specular_light2.SetRange(0,100);
	sprintf(str,"specular: %d%%", m_specular_light2);
	GetDlgItem(IDC_STATIC_SPECULAR_LIGHT2)->SetWindowText(str);

	int initpos = 1000;//this initial position outside helps that the slider is redrawn correctly when opened the second time
	sliderctrl_transparency.SetPos(initpos);
	sliderctrl_shininess.SetPos(initpos);
	sliderctrl_speccolor.SetPos(initpos);
	sliderctrl_ambient_light.SetPos(initpos);
	sliderctrl_diffuse_light.SetPos(initpos);
	sliderctrl_specular_light.SetPos(initpos);
	sliderctrl_ambient_light2.SetPos(initpos);
	sliderctrl_diffuse_light2.SetPos(initpos);
	sliderctrl_specular_light2.SetPos(initpos);

	CRect r;
	pParent->GetWindowRect(&r);
	CPoint p(r.left,r.bottom - 5);
	GetWindowRect(&r);
	p.x += r.Width();

	r.top = p.y - r.Height();
	r.left = p.x - r.Width();
	r.right = p.x;
	r.bottom = p.y;
	MoveWindow(r,FALSE);
	ShowWindow(SW_SHOW);

	UpdateData(FALSE);

}

void DialogOpenGLOptions::LoadData()
{
	/*
	SetDOption(207,0.25);  //light1 ambient parameter
	SetDOption(208,0.40);  //light1 diffuse parameter
	SetDOption(209,0.40);  //light1 specular parameter
	SetDOption(210,0.25);  //light2 ambient parameter
	SetDOption(211,0.40);  //light2 diffuse parameter
	SetDOption(212,0.00);  //light2 specular parameter
	SetDOption(213, 1.);   //light1 posx
	SetDOption(214, 1.);   //light1 posy
	SetDOption(215,-3.);   //light1 posz
	SetDOption(216, 0.);   //light2 posx
	SetDOption(217, 3.);   //light2 posy
	SetDOption(218, 2.);   //light2 posz

	SetIOption(208,1);  //OpenGL enable light1
	SetIOption(209,1);  //OpenGL enable light2
	SetIOption(210,0);  //OpenGL light1 mode (0=standard, 1=use light position)
	SetIOption(211,0);  //OpenGL light2 mode (0=standard, 1=use light position)
*/
	m_check_immediate_apply = pWCDI->GetIOption(212);

	m_check_smooth_shade_model	= pWCDI->GetIOption(207);
	m_check_enable_lighting			= pWCDI->GetIOption(206);
	//double:
	m_transparency							= (int)(100.*pWCDI->GetDOption(206));
	m_speccolor									= (int)(100.*pWCDI->GetDOption(220));
	m_shininess									= (int)(100./128.*pWCDI->GetDOption(219));
	//light1:
	m_check_enable_light = pWCDI->GetIOption(208);
	m_check_enable_light2 = pWCDI->GetIOption(209);
	m_check_incllightpos = pWCDI->GetIOption(210);
	m_check_incllightpos2 = pWCDI->GetIOption(211);
	m_ambient_light	= (int)(100.*pWCDI->GetDOption(207));
	m_diffuse_light	= (int)(100.*pWCDI->GetDOption(208));
	m_specular_light= (int)(100.*pWCDI->GetDOption(209));
	m_ambient_light2= (int)(100.*pWCDI->GetDOption(210));
	m_diffuse_light2= (int)(100.*pWCDI->GetDOption(211));
	m_specular_light2=(int)(100.*pWCDI->GetDOption(212));
	m_lightposx = pWCDI->GetDOption(213);
	m_lightposy = pWCDI->GetDOption(214);
	m_lightposz = pWCDI->GetDOption(215);
	m_lightposx2 = pWCDI->GetDOption(216);
	m_lightposy2 = pWCDI->GetDOption(217);
	m_lightposz2 = pWCDI->GetDOption(218);

}

void DialogOpenGLOptions::WriteData()
{
	pWCDI->GetIOption(212) = m_check_immediate_apply;
	pWCDI->GetIOption(206) = m_check_enable_lighting;
	pWCDI->GetIOption(207) = m_check_smooth_shade_model;
	//double:
	pWCDI->GetDOption(206) = (double)m_transparency/100.;
	pWCDI->GetDOption(219) = 128./100.*(double)m_shininess;
	pWCDI->GetDOption(220) = (double)m_speccolor/100.;
	//light1:
	pWCDI->GetIOption(208) = m_check_enable_light;
	pWCDI->GetIOption(209) = m_check_enable_light2;
	pWCDI->GetIOption(210) = m_check_incllightpos;
	pWCDI->GetIOption(211) = m_check_incllightpos2;
	pWCDI->GetDOption(207) = 0.01*(double)m_ambient_light;
	pWCDI->GetDOption(208) = 0.01*(double)m_diffuse_light;
	pWCDI->GetDOption(209) = 0.01*(double)m_specular_light;
	pWCDI->GetDOption(210) = 0.01*(double)m_ambient_light2;
	pWCDI->GetDOption(211) = 0.01*(double)m_diffuse_light2;
	pWCDI->GetDOption(212) = 0.01*(double)m_specular_light2;
	pWCDI->GetDOption(213) = m_lightposx;
	pWCDI->GetDOption(214) = m_lightposy;
	pWCDI->GetDOption(215) = m_lightposz;
	pWCDI->GetDOption(216) = m_lightposx2;
	pWCDI->GetDOption(217) = m_lightposy2;
	pWCDI->GetDOption(218) = m_lightposz2;

	pWCDI->CallCompFunction(212,1); //1==SetOptions2EDCOptions
}



// DialogOpenGLOptions-Meldungshandler

void DialogOpenGLOptions::OnBnClickedCancel()
{
	OnCancel();
}

void DialogOpenGLOptions::UpdateIfActivated()
{
	UpdateData(TRUE);
	WriteData();

	if (m_check_immediate_apply)
	{
		pGLDrawWnd->Redraw();
	}
}

void DialogOpenGLOptions::OnBnClickedApply()
{
	UpdateData(TRUE);
	WriteData();
 
	pGLDrawWnd->Redraw();
}


void DialogOpenGLOptions::OnBnClickedOk()
{

	OnBnClickedApply();
	
	OnOK();
}

void DialogOpenGLOptions::OnOK()
{
	DestroyWindow();
}

void DialogOpenGLOptions::OnClose() 
{
	DestroyWindow();
}

void DialogOpenGLOptions::OnCancel() 
{
	DestroyWindow();
}
 

BOOL DialogOpenGLOptions::OnNotify(WPARAM wParam, LPARAM lParam, LRESULT* pResult)
{
	// TODO: Fügen Sie hier Ihren spezialisierten Code ein, und/oder rufen Sie die Basisklasse auf.
	//OnBnClickedApply(); //does not work???

	return CDialog::OnNotify(wParam, lParam, pResult);
}

void DialogOpenGLOptions::OnNMCustomdrawSliderTransparency(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"transparency: %d%%", 100-m_transparency);
	GetDlgItem(IDC_STATIC_TRANSPARENCY)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderShininess(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"shininess: %d%%", m_shininess);
	GetDlgItem(IDC_STATIC_SHININESS)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderSpeccolor(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"specular color intensity: %d%%", m_speccolor);
	GetDlgItem(IDC_STATIC_SPECCOLOR)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderAmbientLight(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"ambient: %d%%", m_ambient_light);
	GetDlgItem(IDC_STATIC_AMBIENT_LIGHT)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderDiffuseLight(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"diffuse: %d%%", m_diffuse_light);
	GetDlgItem(IDC_STATIC_DIFFUSE_LIGHT)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderSpecularLight(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"specular: %d%%", m_specular_light);
	GetDlgItem(IDC_STATIC_SPECULAR_LIGHT)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposx()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposy()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposz()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnBnClickedCheckIncllightpos()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnBnClickedCheckEnableLight()
{
	UpdateIfActivated();
}



void DialogOpenGLOptions::OnBnClickedCheck1() //automatic apply
{
	OnBnClickedApply();
}

void DialogOpenGLOptions::OnBnClickedCheckEnableLight2()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnBnClickedCheckIncllightpos2()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnNMCustomdrawSliderAmbientLight2(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"ambient: %d%%", m_ambient_light2);
	GetDlgItem(IDC_STATIC_AMBIENT_LIGHT2)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderDiffuseLight2(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"diffuse: %d%%", m_diffuse_light2);
	GetDlgItem(IDC_STATIC_DIFFUSE_LIGHT2)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnNMCustomdrawSliderSpecularLight2(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateIfActivated();
	char str[50];
	sprintf(str,"specular: %d%%", m_specular_light2);
	GetDlgItem(IDC_STATIC_SPECULAR_LIGHT2)->SetWindowText(str);

	*pResult = 0;
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposx2()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposy2()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnEnKillfocusEditLightposz2()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnBnClickedCheckSmoothshademodel()
{
	UpdateIfActivated();
}

void DialogOpenGLOptions::OnBnClickedCheckEnablelighting()
{
	UpdateIfActivated();
}
