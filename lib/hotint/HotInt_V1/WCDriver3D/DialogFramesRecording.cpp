//#**************************************************************
//# filename:             DialogFramesRecording.cpp
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
#include "DialogFramesRecording.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialogFramesRecording dialog


CDialogFramesRecording::CDialogFramesRecording(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogFramesRecording::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialogFramesRecording)
	m_bCheckRecordFrames = FALSE; //B163 // still required to be set to zero ( CGLDrawWnd::Init() called before DialogFrameRecording::LoadData() ) 
	//m_strPathToImageFiles = _T("D:\\tmp\\FramesRecording\\"); //T111&T113
	//m_nRecordEachFrameOf = 1; //I166
	m_nWindowSizeX = 0;
	m_nWindowSizeY = 0;
	m_nFrameCounter = 1;
	//m_bShowFrameNumbers = FALSE;//B165
	//m_bProcessImage = FALSE;//B164
	//}}AFX_DATA_INIT
}


void CDialogFramesRecording::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialogFramesRecording)
	DDX_Check(pDX, IDC_CHECK_RECORD_FRAMES, m_bCheckRecordFrames);                     //B163
	DDX_Text(pDX, IDC_EDIT_PATH_TO_IMAGE_FILES, m_strPathToSingleImageFiles);          //T120
	DDX_Text(pDX, IDC_EDIT_IMAGE_FILE_NAME, m_strSingleImageFileName);                 //T121
	DDX_Text(pDX, IDC_EDIT_PATH_TO_VIDEO_FILES, m_strPathToVideoImageFiles);           //T122
	DDX_Text(pDX, IDC_EDIT_VIDEO_FILE_NAME, m_strVideoImageFileName);                  //T123

	DDX_Text(pDX, IDC_EDIT_RECORD_EACH_FRAME_OF, m_nRecordEachFrameOf);                //I166          
	DDV_MinMaxInt(pDX, m_nRecordEachFrameOf, 1, 10000);
	DDX_Text(pDX, IDC_EDIT_WINDOW_SIZE_X, m_nWindowSizeX);
	DDX_Text(pDX, IDC_EDIT_WINDOW_SIZE_Y, m_nWindowSizeY);
	DDX_Text(pDX, IDC_EDIT_FRAME_COUNTER, m_nFrameCounter);
	DDV_MinMaxInt(pDX, m_nFrameCounter, 0, 100000);
	DDX_Check(pDX, IDC_CHECK_SHOW_FRAME_NUMBER, m_bShowFrameNumbers);                  //B165
//	DDX_Check(pDX, IDC_CHECK_PROCESS_IMAGE, m_bProcessImage);                          //B164
	//}}AFX_DATA_MAP
	DDX_Radio(pDX, IDC_RADIO_JPG, m_radio_fileformat);
}


BEGIN_MESSAGE_MAP(CDialogFramesRecording, CDialog)
	//{{AFX_MSG_MAP(CDialogFramesRecording)
	ON_BN_CLICKED(IDC_CHECK_RECORD_FRAMES, OnCheckRecordFrames)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDOK, &CDialogFramesRecording::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON_BROWSESINGLE, &CDialogFramesRecording::OnBnClickedButtonBrowsesingle)
	ON_BN_CLICKED(IDC_BUTTON_BROWSEVIDEO, &CDialogFramesRecording::OnBnClickedButtonBrowsevideo)
	ON_BN_CLICKED(IDC_RADIO_JPG, &CDialogFramesRecording::OnBnClickedRadioFormat)
	ON_BN_CLICKED(IDC_RADIO_PNG, &CDialogFramesRecording::OnBnClickedRadioFormat)
	ON_BN_CLICKED(IDC_RADIO_BMP, &CDialogFramesRecording::OnBnClickedRadioFormat)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialogFramesRecording message handlers

void CDialogFramesRecording::OnCheckRecordFrames() 
{
	bCurrentCheckRecordFrames = !bCurrentCheckRecordFrames;
	UpdateEnableControls();
}

BOOL CDialogFramesRecording::OnInitDialog() 
{
	CDialog::OnInitDialog();
	LoadData();

	bCurrentCheckRecordFrames = m_bCheckRecordFrames;
	UpdateEnableControls();
	
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CDialogFramesRecording::OnBnClickedOk()
{
	WriteData();
	OnOK();
}

void CDialogFramesRecording::UpdateEnableControls()
{
	// depending on CheckBox "Record Frames" the items are enabled or disabled
	GetDlgItem(IDC_EDIT_PATH_TO_IMAGE_FILES)->EnableWindow(!bCurrentCheckRecordFrames);
	GetDlgItem(IDC_EDIT_IMAGE_FILE_NAME)->EnableWindow(!bCurrentCheckRecordFrames);
	GetDlgItem(IDC_BUTTON_BROWSESINGLE)->EnableWindow(!bCurrentCheckRecordFrames);

	GetDlgItem(IDC_EDIT_PATH_TO_VIDEO_FILES)->EnableWindow(bCurrentCheckRecordFrames);
	GetDlgItem(IDC_EDIT_VIDEO_FILE_NAME)->EnableWindow(bCurrentCheckRecordFrames);
	GetDlgItem(IDC_BUTTON_BROWSEVIDEO)->EnableWindow(bCurrentCheckRecordFrames);

	GetDlgItem(IDC_EDIT_RECORD_EACH_FRAME_OF)->EnableWindow(bCurrentCheckRecordFrames);
	GetDlgItem(IDC_EDIT_FRAME_COUNTER)->EnableWindow(bCurrentCheckRecordFrames);
	GetDlgItem(IDC_CHECK_PROCESS_IMAGE)->EnableWindow(bCurrentCheckRecordFrames);
	GetDlgItem(IDC_CHECK_SHOW_FRAME_NUMBER)->EnableWindow(bCurrentCheckRecordFrames);
}

////will be erased:
//void CDialogFramesRecording::Serialize(CArchive & ar) ///AD removed: 2012-11-09
//{
//	if(ar.IsLoading())
//		ar >> m_bCheckRecordFrames >> m_bProcessImage >> m_bShowFrameNumbers >>
//			m_nRecordEachFrameOf >> m_strPathToImageFiles;
//	else
//		ar << m_bCheckRecordFrames << m_bProcessImage << m_bShowFrameNumbers <<
//			m_nRecordEachFrameOf << m_strPathToImageFiles;
//}

void CDialogFramesRecording::Configuration2EDC(ElementDataContainer& edc)
{
	// not in hotint options dialog
	// moved to SaveData()
	//edc.TreeSetBoolC("ViewingOptions.Animation.RecordSingleFrames.record",m_bCheckRecordFrames,"record frames");
	//edc.TreeSetBoolC("ViewingOptions.Animation.RecordSingleFrames.process_image",m_bProcessImage,"for each saved frame, call conversion program (Image Magick)");
	//edc.TreeSetBoolC("ViewingOptions.Animation.RecordSingleFrames.show_frame_numbers",m_bShowFrameNumbers,"show frame numbers in images");
	//edc.TreeSetIntC("ViewingOptions.Animation.RecordSingleFrames.record_every_x_frame",m_nRecordEachFrameOf,"record every x frames");
	//edc.TreeSetStringC("GeneralOptions.Paths.record_frames_path",m_strPathToImageFiles,"path to image files");
}

void CDialogFramesRecording::EDC2Configuration(const ElementDataContainer& edc)
{
	// not in hotint options dialog
	// moved to LoadData()
	//edc.TreeGetInt("ViewingOptions.Animation.RecordSingleFrames.record",m_bCheckRecordFrames);
	//edc.TreeGetInt("ViewingOptions.Animation.RecordSingleFrames.process_image",m_bProcessImage);
	//edc.TreeGetInt("ViewingOptions.Animation.RecordSingleFrames.show_frame_numbers",m_bShowFrameNumbers);
	//edc.TreeGetInt("ViewingOptions.Animation.RecordSingleFrames.record_every_x_frame",m_nRecordEachFrameOf);
	//m_strPathToImageFiles = edc.TreeGetString("GeneralOptions.Paths.record_frames_path");
}

void CDialogFramesRecording::LoadData()
{
// not in hotint options dialog
	//pWCDI->MBS_EDC_TreeGetInt(m_bCheckRecordFrames,"ViewingOptions.Animation.RecordSingleFrames.record");
	//pWCDI->MBS_EDC_TreeGetInt(m_bCheckIncludeOutputWindow,"ViewingOptions.Animation.RecordSingleFrames.include_output_window");
	//pWCDI->MBS_EDC_TreeGetInt(m_bProcessImage,"ViewingOptions.Animation.RecordSingleFrames.process_image");
	//pWCDI->MBS_EDC_TreeGetInt(m_bShowFrameNumbers,"ViewingOptions.Animation.RecordSingleFrames.show_frame_numbers");
	//pWCDI->MBS_EDC_TreeGetInt(m_nRecordEachFrameOf,"ViewingOptions.Animation.RecordSingleFrames.record_every_x_frame");
	m_bCheckRecordFrames = pWCDI->GetIOption(163);
	m_bCheckIncludeOutputWindow = pWCDI->GetIOption(168);
//	m_bProcessImage = pWCDI->GetIOption(164);
	m_bShowFrameNumbers = pWCDI->GetIOption(165);
	m_nRecordEachFrameOf = pWCDI->GetIOption(166);
	m_radio_fileformat = pWCDI->GetIOption(167);

	m_strPathToSingleImageFiles = pWCDI->GetTOption(120);
	m_strSingleImageFileName = pWCDI->GetTOption(121);
	m_strPathToVideoImageFiles = pWCDI->GetTOption(122);
	m_strVideoImageFileName = pWCDI->GetTOption(123);      
	
	//m_strPathToImageFiles = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.record_frames_path");
//  m_strPathToImageFiles = pWCDI->GetTOption(102);
	UpdateData(FALSE);
}

void CDialogFramesRecording::WriteData()
{
	UpdateData(TRUE);
	// not in hotint options dialog
	//pWCDI->MBS_EDC_TreeSetInt(m_bCheckRecordFrames,"ViewingOptions.Animation.RecordSingleFrames.record");
	//pWCDI->MBS_EDC_TreeSetInt(m_bCheckIncludeOutputWindow,"ViewingOptions.Animation.RecordSingleFrames.include_output_window");
	//pWCDI->MBS_EDC_TreeSetInt(m_bProcessImage,"ViewingOptions.Animation.RecordSingleFrames.process_image");
	//pWCDI->MBS_EDC_TreeSetInt(m_bShowFrameNumbers,"ViewingOptions.Animation.RecordSingleFrames.show_frame_numbers");
	//pWCDI->MBS_EDC_TreeSetInt(m_nRecordEachFrameOf,"ViewingOptions.Animation.RecordSingleFrames.record_every_x_frame");
	pWCDI->SetIOption(163,m_bCheckRecordFrames);
	pWCDI->SetIOption(168,m_bCheckIncludeOutputWindow);
//	pWCDI->SetIOption(164,m_bProcessImage);
	pWCDI->SetIOption(165,m_bShowFrameNumbers);
	pWCDI->SetIOption(166,m_nRecordEachFrameOf);
	pWCDI->SetIOption(167,m_radio_fileformat);

	pWCDI->SetTOption(120,m_strPathToSingleImageFiles);
	pWCDI->SetTOption(121,m_strSingleImageFileName);
	pWCDI->SetTOption(122,m_strPathToVideoImageFiles);
	pWCDI->SetTOption(123,m_strVideoImageFileName);

//	pWCDI->MBS_EDC_TreeSetString(m_strPathToImageFiles,"GeneralOptions.Paths.record_frames_path");
	//pWCDI->SetTOption(102,m_strPathToImageFiles);
}

void CDialogFramesRecording::OnBnClickedButtonBrowsesingle()
{
	UpdateData(TRUE);                       // in case the user modified the editbox manually
	OnBrowse(m_strPathToSingleImageFiles,"Select a Folder for Output (single image files)");
	UpdateData(FALSE);
}

void CDialogFramesRecording::OnBnClickedButtonBrowsevideo()
{
	UpdateData(TRUE);                       // in case the user modified the editbox manually
	OnBrowse(m_strPathToVideoImageFiles,"Select a Folder for Output (video files)");
	UpdateData(FALSE);
}

void CDialogFramesRecording::OnBrowse(CString& pathstring, const char* title)
{
	const CString dummy= (const CString) pathstring;
	char buffer[1024];
	sprintf(buffer,dummy);
  
	BROWSEINFO bi;
	bi.hwndOwner = m_hWnd; // Handle of the owner window
	bi.pidlRoot = NULL; // Desktop folder is used
	bi.lpszTitle = title;
	bi.pszDisplayName = buffer; // Buffer for selected folder name
//	bi.ulFlags = BIF_RETURNONLYFSDIRS; // Only returns file system directories
//	bi.ulFlags = BIF_RETURNONLYFSDIRS|BIF_NEWDIALOGSTYLE; // <--- this should create the NEW FOLDER button
	bi.ulFlags = BIF_EDITBOX | BIF_VALIDATE /*| BIF_BROWSEFORPRINTER*/ | BIF_RETURNONLYFSDIRS|0x0040; // <--- and this one works... :)
	bi.lpfn = NULL;
	bi.lParam = (LPARAM) buffer;

	LPITEMIDLIST pItemIDList = SHBrowseForFolder(&bi);
	if (pItemIDList)
	{
		if (SHGetPathFromIDList(pItemIDList,buffer)) 	
		{
		 pathstring = buffer;
		}
	}
}
void CDialogFramesRecording::OnBnClickedRadioFormat()
{
	UpdateData(TRUE);
	if(m_radio_fileformat == 0) //JPG
	{
		;
	}
	else if(m_radio_fileformat == 1) //PNG
	{
		;
	}
	else if(m_radio_fileformat == 2) //BMP
	{
		;
	}
	UpdateData(FALSE);
}
