//#**************************************************************
//# filename:             DlgBodyJointOpt.cpp
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
#include "DlgBodyJointOpt.h"


// DialogBodyJointOptions-Dialogfeld

IMPLEMENT_DYNAMIC(DialogBodyJointOptions, CDialog)
DialogBodyJointOptions::DialogBodyJointOptions(CWnd* pParent /*=NULL*/)
	: CDialog(DialogBodyJointOptions::IDD, pParent)
	, m_check_show_joints(FALSE)
	, m_show_body_numbers(FALSE)
	, m_check_usedegrees(FALSE)
	, m_check_show_constraint_numbers(FALSE)
	, m_check_showbodylocalframe(FALSE)
	, m_bodylocalframesize(0)
	, m_check_showsensors(FALSE)
	, m_sensor_size(0)
	, m_check_joints_transparent(FALSE)
	, m_check_sensors_transparent(FALSE)
	, m_check_bodies_transparent(FALSE)
	, m_check_bodies_supersmooth(FALSE)
	, m_check_show_body_outline(FALSE)
	, m_check_show_body_faces(FALSE)
	, m_radio_eulerangles(0)
	, m_check_showloads(FALSE)
	, m_load_draw_size(0)
	//, m_check_show_control_ojects(FALSE)
{
	EnableActiveAccessibility();
}

DialogBodyJointOptions::~DialogBodyJointOptions()
{
}

void DialogBodyJointOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK_SHOW_JOINTS, m_check_show_joints);
	DDX_Check(pDX, IDC_CHECK_SHOW_BODY_NUMBERS, m_show_body_numbers);
	DDX_Check(pDX, IDC_CHECK_SHOW_CONSTRAINT_NUMBERS, m_check_show_constraint_numbers);
	DDX_Check(pDX, IDC_CHECK_SHOW_BODY_LOCAL_FRAME, m_check_showbodylocalframe);
	DDX_Text(pDX, IDC_EDIT_BODYCOORDSIZE, m_bodylocalframesize);
	DDX_Check(pDX, IDC_CHECK_SHOW_SENSORS, m_check_showsensors);
	DDX_Text(pDX, IDC_EDIT_SENSORSIZE, m_sensor_size);
	DDX_Check(pDX, IDC_CHECK_JOINTS_TRANSPARENT, m_check_joints_transparent);
	DDX_Check(pDX, IDC_CHECK_SENSORS_TRANSPARENT2, m_check_sensors_transparent);
	DDX_Check(pDX, IDC_CHECK_SHOW_BODIES_TRANSPARENT, m_check_bodies_transparent);
	DDX_Check(pDX, IDC_CHECK_DRAW_BODIES_SUPERSMOOTH, m_check_bodies_supersmooth);
	DDX_Check(pDX, IDC_CHECK_SHOW_BODIES_OUTLINE, m_check_show_body_outline);
	DDX_Check(pDX, IDC_CHECK_SHOW_BODIES_FACES, m_check_show_body_faces);
	DDX_Check(pDX, IDC_CHECK_SHOW_LOADS, m_check_showloads);
	DDX_Text(pDX, IDC_EDIT_LOADDRAWSIZE, m_load_draw_size);
	DDV_MinMaxDouble(pDX, m_load_draw_size, 0, 1e32);
}


BEGIN_MESSAGE_MAP(DialogBodyJointOptions, CDialog)
	ON_BN_CLICKED(ID_APPLY, OnBnClickedApply)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
END_MESSAGE_MAP()


void DialogBodyJointOptions::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);

	LoadData();

	//initialize sliders

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

void DialogBodyJointOptions::LoadData()
{
	m_check_show_joints				= pWCDI->GetIOption(121);
	m_show_body_numbers				= pWCDI->GetIOption(123);
	m_check_show_constraint_numbers = pWCDI->GetIOption(124);

	m_check_showbodylocalframe= pWCDI->GetIOption(125);
	m_bodylocalframesize			= pWCDI->GetDOption(104);

	m_check_showsensors				= pWCDI->GetIOption(127);
	m_sensor_size							= pWCDI->GetDOption(118);

	m_check_bodies_transparent = pWCDI->GetIOption(128);
	m_check_joints_transparent = pWCDI->GetIOption(129);
	m_check_sensors_transparent= pWCDI->GetIOption(130);

	m_check_bodies_supersmooth = pWCDI->GetIOption(132);
	m_check_show_body_outline  = pWCDI->GetIOption(133);
	m_check_show_body_faces		 = pWCDI->GetIOption(134);

	m_check_showloads					 = pWCDI->GetIOption(142);
	m_load_draw_size					 = pWCDI->GetDOption(120);

	//m_check_show_control_ojects= pWCDI->GetIOption(144);
}

void DialogBodyJointOptions::WriteData()
{
	pWCDI->GetIOption(121) = m_check_show_joints;
	pWCDI->GetIOption(123) = m_show_body_numbers; 
	pWCDI->GetIOption(124) = m_check_show_constraint_numbers;

	pWCDI->GetIOption(125) = m_check_showbodylocalframe;
	pWCDI->GetDOption(104) = m_bodylocalframesize;

	pWCDI->GetIOption(127) = m_check_showsensors;
	pWCDI->GetDOption(118) = m_sensor_size;

	pWCDI->GetIOption(128) = m_check_bodies_transparent;
	pWCDI->GetIOption(129) = m_check_joints_transparent;
	pWCDI->GetIOption(130) = m_check_sensors_transparent;

	pWCDI->GetIOption(132) = m_check_bodies_supersmooth;
	pWCDI->GetIOption(133) = m_check_show_body_outline;
	pWCDI->GetIOption(134) = m_check_show_body_faces;

	pWCDI->GetIOption(142) = m_check_showloads;

	pWCDI->GetDOption(120) = m_load_draw_size;

	//pWCDI->GetIOption(144) = m_check_show_control_ojects;

	pWCDI->CallCompFunction(212,1); //1==SetOptions2EDCOptions

}



// DialogBodyJointOptions-Meldungshandler

void DialogBodyJointOptions::OnBnClickedCancel()
{
	OnCancel();
}

void DialogBodyJointOptions::OnBnClickedApply()
{
	UpdateData(TRUE);
	WriteData();

	//general operations on apply if exist

	pGLDrawWnd->Redraw();
}

void DialogBodyJointOptions::OnBnClickedOk()
{
	OnBnClickedApply();

	OnOK();
}

void DialogBodyJointOptions::OnOK()
{
	DestroyWindow();
}

void DialogBodyJointOptions::OnClose() 
{
	DestroyWindow();
}

void DialogBodyJointOptions::OnCancel() 
{
	DestroyWindow();
}
 

void DialogBodyJointOptions::OnBnClickedCheckShowJoints2()
{
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.
}


