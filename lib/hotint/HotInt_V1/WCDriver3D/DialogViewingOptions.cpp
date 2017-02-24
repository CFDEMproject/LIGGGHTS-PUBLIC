//#**************************************************************
//# filename:             DialogViewingOptions.cpp
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
#include "DialogViewingOptions.h"


// DialogViewingOptions-Dialogfeld

IMPLEMENT_DYNAMIC(DialogViewingOptions, CDialog)
DialogViewingOptions::DialogViewingOptions(CWnd* pParent /*=NULL*/)
	: CDialog(DialogViewingOptions::IDD, pParent)
	, m_redrawfrequency(0)
	, m_check_draworigin(FALSE)
	//, m_check_showcontactpoints(FALSE)
	, m_standardview_angle1(0)
	, m_standardview_angle2(0)
	, m_standardview_angle3(0)
	, m_standardview_axis1(0)
	, m_standardview_axis2(0)
	, m_standardview_axis3(0)
	, m_WindowSizeX(0)
	, m_WindowSizeY(0)
	, m_animate_beginning(FALSE)
	, m_animationframes(0)
	, m_check_textstofront(FALSE)
	, m_origin_size(0)
	, m_gridsize1(0)
	, m_gridstep1(0)
	, m_combo_gridtype(0)
	, m_grid_posx(0)
	, m_grid_posy(0)
	, m_grid_posz(0)
	//, m_check_show_startup_banner(FALSE)
	, m_check_use_cuttingplane(FALSE)
	, m_cutplane_distance(0)
	, m_cutplane_nx(0)
	, m_cutplane_ny(0)
	, m_cutplane_nz(0)
{
}

DialogViewingOptions::~DialogViewingOptions()
{
}

void DialogViewingOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Slider(pDX, IDC_SLIDER_REDRAWFREQUENCY, m_redrawfrequency);
	DDX_Check(pDX, IDC_CHECK_DRAWORIGIN, m_check_draworigin);
	//DDX_Check(pDX, IDC_CHECK_SHOWCONTACTPOINTS, m_check_showcontactpoints);
	DDX_Text(pDX, IDC_EDIT1, m_standardview_angle1);
	DDV_MinMaxDouble(pDX, m_standardview_angle1, -180, 360);
	DDX_Text(pDX, IDC_EDIT3, m_standardview_angle2);
	DDV_MinMaxDouble(pDX, m_standardview_angle2, -180, 360);
	DDX_Text(pDX, IDC_EDIT5, m_standardview_angle3);
	DDV_MinMaxDouble(pDX, m_standardview_angle3, -180, 360);
	DDX_Text(pDX, IDC_EDIT2, m_standardview_axis1);
	DDV_MinMaxInt(pDX, m_standardview_axis1, 1, 3);
	DDX_Text(pDX, IDC_EDIT4, m_standardview_axis2);
	DDV_MinMaxInt(pDX, m_standardview_axis2, 1, 3);
	DDX_Text(pDX, IDC_EDIT6, m_standardview_axis3);
	DDV_MinMaxDouble(pDX, m_standardview_axis3, 1, 3);
	DDX_Text(pDX, IDC_EDIT_WINDOWSIZEX, m_WindowSizeX);
	DDV_MinMaxInt(pDX, m_WindowSizeX, 20, 10000);
	DDX_Text(pDX, IDC_EDIT_WINDOWSIZEY, m_WindowSizeY);
	DDV_MinMaxInt(pDX, m_WindowSizeY, 20, 10000);
	DDX_Control(pDX, IDC_SLIDER_REDRAWFREQUENCY, slider_redrawfrequency);
	DDX_Check(pDX, IDC_ANIMATE_BEGINNING, m_animate_beginning);
	DDX_Slider(pDX, IDC_SLIDER_ANIMATIONFRAMES, m_animationframes);
	DDX_Control(pDX, IDC_SLIDER_ANIMATIONFRAMES, slider_animationframes);
	DDX_Check(pDX, IDC_CHECK_TEXTSTOFRONT, m_check_textstofront);
	DDX_Text(pDX, IDC_EDIT_ORIGINSIZE, m_origin_size);
	DDX_Text(pDX, IDC_EDIT_GRIDSIZE1, m_gridsize1);
	DDX_Text(pDX, IDC_EDIT_GRIDSTEP1, m_gridstep1);
	DDX_Control(pDX, IDC_COMBO_GRIDTYPE, combo_gridtype);
	DDX_CBIndex(pDX, IDC_COMBO_GRIDTYPE, m_combo_gridtype);
	DDX_Text(pDX, IDC_EDIT_GRIDPOSX, m_grid_posx);
	DDX_Text(pDX, IDC_EDIT_GRIDPOSY, m_grid_posy);
	DDX_Text(pDX, IDC_EDIT_GRIDPOSZ, m_grid_posz);
	//DDX_Check(pDX, IDC_STARTUP_BANNER, m_check_show_startup_banner);
	DDX_Check(pDX, IDC_USE_CUTTING_PLANE, m_check_use_cuttingplane);
	DDX_Text(pDX, IDC_EDIT_CUTPLANE_DISTANCE, m_cutplane_distance);
	DDX_Text(pDX, IDC_EDIT_CUTPLANE_X, m_cutplane_nx);
	DDX_Text(pDX, IDC_EDIT_CUTPLANE_Y, m_cutplane_ny);
	DDX_Text(pDX, IDC_EDIT_CUTPLANE_Z, m_cutplane_nz);
}


BEGIN_MESSAGE_MAP(DialogViewingOptions, CDialog)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
	ON_BN_CLICKED(IDCANCEL2, OnBnClickedApply)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_REDRAWFREQUENCY, OnNMCustomdrawSliderRedrawfrequency)
	ON_BN_CLICKED(IDC_CHECK_DRAWORIGIN, OnBnClickedCheckDraworigin)
	ON_NOTIFY(NM_RELEASEDCAPTURE, IDC_SLIDER_REDRAWFREQUENCY, OnNMReleasedcaptureSliderRedrawfrequency)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_ANIMATIONFRAMES, OnNMCustomdrawSliderAnimationframes)
END_MESSAGE_MAP()



void DialogViewingOptions::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	LoadData();

	combo_gridtype.Clear();
	combo_gridtype.AddString("no grid");
	combo_gridtype.AddString("grid YZ");
	combo_gridtype.AddString("grid XZ");
	combo_gridtype.AddString("grid XY");

	slider_redrawfrequency.SetRange(0,9);
	slider_animationframes.SetRange(1,13);

	slider_redrawfrequency.SetPos(100);
	slider_animationframes.SetPos(100);

	char str[20];
	int af = m_animationframes;
	if (m_animationframes == 11) af = 15;
	if (m_animationframes == 12) af = 20;
	if (m_animationframes == 13) af = 50;

	if (af==1) sprintf(str,"every frame");
	else sprintf(str,"every %d frames", af);
	GetDlgItem(IDC_STATIC_ANIMATIONFRAMES)->SetWindowText(str);

	UpdateData(FALSE);

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

}

void DialogViewingOptions::LoadData()
{

	//set actual WCDriver-dialog window size in options dialog
	CRect rd;
	pGLDrawWnd->GetWindowRect(&rd);
	m_WindowSizeX = rd.Width();
	m_WindowSizeY = rd.Height();

	//misc drawing options:
	m_redrawfrequency       = pWCDI->GetIOption(104);
	//m_check_showcontactpoints = pWCDI->GetIOption(112);
	m_check_draworigin      = pWCDI->GetIOption(106);

	m_animationframes = pWCDI->GetIOption(109);
	if (m_animationframes == 15) m_animationframes = 11;
	if (m_animationframes == 20) m_animationframes = 12;
	if (m_animationframes == 50) m_animationframes = 13;

	m_animate_beginning = pWCDI->GetIOption(115);
	m_check_textstofront = pWCDI->GetIOption(122);

	m_origin_size			= pWCDI->GetDOption(103);
	m_combo_gridtype	= pWCDI->GetIOption(126);
	m_gridsize1				= pWCDI->GetDOption(112);
	m_gridsize1				= pWCDI->GetDOption(113); //not implemented yet
	m_gridstep1				= pWCDI->GetDOption(110);
	m_gridstep1				= pWCDI->GetDOption(111); //not implemented yet
	m_grid_posx				= pWCDI->GetDOption(107);
	m_grid_posy				= pWCDI->GetDOption(108);
	m_grid_posz				= pWCDI->GetDOption(109);

	//standardview:
	m_standardview_angle1 = pWCDI->GetDOption(200);
	m_standardview_angle2 = pWCDI->GetDOption(201);
	m_standardview_angle3 = pWCDI->GetDOption(202);
	m_standardview_axis1 =  pWCDI->GetIOption(200);
	m_standardview_axis2 =  pWCDI->GetIOption(201);
	m_standardview_axis3 =  pWCDI->GetIOption(202);

	//m_check_show_startup_banner = pWCDI->GetIOption(141);

	//cutplane:
	m_check_use_cuttingplane = pWCDI->GetIOption(149);
	m_cutplane_distance = pWCDI->GetDOption(128);
	m_cutplane_nx = pWCDI->GetDOption(125);
	m_cutplane_ny = pWCDI->GetDOption(126);
	m_cutplane_nz = pWCDI->GetDOption(127);
}

void DialogViewingOptions::WriteData()
{

	pWCDI->SetIOption(104, m_redrawfrequency);
	pWCDI->SetIOption(106, (int)m_check_draworigin);
	//pWCDI->SetIOption(112, m_check_showcontactpoints);
	pWCDI->SetIOption(122, m_check_textstofront);

	int af = m_animationframes;
	if (m_animationframes == 11) af = 15;
	if (m_animationframes == 12) af = 20;
	if (m_animationframes == 13) af = 50;
	pWCDI->SetIOption(109, af);
	pWCDI->SetIOption(115, m_animate_beginning);

	pWCDI->GetDOption(103) = m_origin_size;
	pWCDI->GetIOption(126) = m_combo_gridtype;
	pWCDI->GetDOption(112) = m_gridsize1;
	pWCDI->GetDOption(113) = m_gridsize1; //not implemented yet
	pWCDI->GetDOption(110) = m_gridstep1;
	pWCDI->GetDOption(111) = m_gridstep1; //not implemented yet
	pWCDI->GetDOption(107) = m_grid_posx;
	pWCDI->GetDOption(108) = m_grid_posy;
	pWCDI->GetDOption(109) = m_grid_posz;

	//  standard view:
	pWCDI->SetDOption(200, m_standardview_angle1);
	pWCDI->SetDOption(201, m_standardview_angle2);
	pWCDI->SetDOption(202, m_standardview_angle3);
	pWCDI->SetIOption(200, m_standardview_axis1);
	pWCDI->SetIOption(201, m_standardview_axis2);
	pWCDI->SetIOption(202, (int)m_standardview_axis3);

	//pWCDI->GetIOption(141) = m_check_show_startup_banner;

	pWCDI->GetIOption(149) = m_check_use_cuttingplane;
	pWCDI->GetDOption(128) = m_cutplane_distance;

	pWCDI->GetDOption(125) = m_cutplane_nx;
	pWCDI->GetDOption(126) = m_cutplane_ny;
	pWCDI->GetDOption(127) = m_cutplane_nz;

	pWCDI->CallCompFunction(212,1); //1==SetOptions2EDCOptions
}



// DialogViewingOptions-Meldungshandler

void DialogViewingOptions::OnBnClickedCancel()
{
	OnCancel();
}

void DialogViewingOptions::OnBnClickedApply()
{
	UpdateData(TRUE);
	WriteData();

	//general operations on apply:
	CRect r, rd;
	GetParent()->GetWindowRect(&r);
	pGLDrawWnd->GetWindowRect(&rd);

	if (rd.Width() != m_WindowSizeX || rd.Height() != m_WindowSizeY)
	{
		int offx = r.Width() - rd.Width();
		int offy = r.Height() - rd.Height();

		CRect rect;
		GetWindowRect(&rect);
		GetParent()->MoveWindow(r.left,r.top,m_WindowSizeX+offx,m_WindowSizeY+offy);
	}

	pGLDrawWnd->Redraw();
}


void DialogViewingOptions::OnBnClickedOk()
{

	OnBnClickedApply();
	
	OnOK();
}

void DialogViewingOptions::OnOK()
{
	DestroyWindow();
}

void DialogViewingOptions::OnClose() 
{
	DestroyWindow();
}

void DialogViewingOptions::OnCancel() 
{
	DestroyWindow();
}
 

char* redrawFrequencyText2 [] = {"off", "draw last frame", "every 100sec", "every 20sec",
"every 2sec", "every 200ms", "every 50ms", "every 20ms", "every 10 frames", "every frame"};

void DialogViewingOptions::OnNMCustomdrawSliderRedrawfrequency(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	GetDlgItem(IDC_STATIC_REDRAWFREQUENCYTEXT)->SetWindowText(redrawFrequencyText2[m_redrawfrequency]);

	*pResult = 0;
}
 


void DialogViewingOptions::OnBnClickedCheckDraworigin()
{
	OnBnClickedApply(); //update data
}

void DialogViewingOptions::OnNMReleasedcaptureSliderRedrawfrequency(NMHDR *pNMHDR, LRESULT *pResult)
{
	OnBnClickedApply();
	*pResult = 0;
}

void DialogViewingOptions::OnNMCustomdrawSliderAnimationframes(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	UpdateData(TRUE);
	char str[20];
	int af = m_animationframes;
	if (m_animationframes == 11) af = 15;
	if (m_animationframes == 12) af = 20;
	if (m_animationframes == 13) af = 50;

	if (af==1) sprintf(str,"every frame");
	else sprintf(str,"every %d frames", af);
	GetDlgItem(IDC_STATIC_ANIMATIONFRAMES)->SetWindowText(str);

	*pResult = 0;
}
