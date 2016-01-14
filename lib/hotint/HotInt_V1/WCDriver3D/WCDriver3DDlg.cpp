//#**************************************************************
//# filename:             WCDriver3DDlg.cpp
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

#include <strsafe.h>

#include "ioincludes.h"
#include <time.h>
#include <sys/timeb.h>

#include "WCDriver3DDlg.h"
#include "WCDriver3D.h"
#include "DialogReadText.h"

#include "DlgOneEditCtrl.h"
#include "CustomEditDialog.h"
#include "CListDlg.h"
#include "TreeViewCustomEditDialog.h"

 


#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


// the name of the configuration file
#define CONFIG_FILE_NAME "hotint.cfg"


// width of the zone between the edit box and the image
#define NEUTRAL_ZONE_WIDTH 0
#define BUTTONS_OFFSET 4

// the windows message identifier
#define WM_UPDATE_TEXT (WM_USER + 1)
#define WM_UPDATE_PLOT (WM_USER + 11)
#define	WM_START_COMPUTATION (WM_USER + 2)

// the timer ID
#define OUTPUT_REFRESH_TIMER_ID (WM_USER + 101)
#define PLOTTOOL_REFRESH_TIMER_ID (WM_USER + 111)

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();
	void SetpWCDI(WCDInterface* pWCDIi);
	WCDInterface* pWCDI;
	// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL
 
	BOOL OnInitDialog();
		// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::SetpWCDI(WCDInterface* pWCDIi)
{
	pWCDI = pWCDIi;
}

BOOL CAboutDlg::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
	char* str = pWCDI->GetHotintVersion().GetString();

	CWnd * atext = GetDlgItem(IDC_ABOUTTEXT);
	CString cst = "HOTINT\n";
	cst += "A multibody system simulation software\n";
	cst += "Johannes Gerstmayr (main developer)\n\n";
	cst += "(HOTINT MBS) Version: ";
  cst += str;
  cst += "\nBuild Date: ";
	cst += __DATE__;
	//cst += " ";
	//cst += __TIME__;
	cst += "\n\n";
	cst += "Parts of HOTINT developed at: \nLCM, ACCM and Institute of Technical Mechanics (JKU Linz) \n(www.lcm.at, www.accm.co.at, tmech.mechatronik.uni-linz.ac.at)\n\n";
	cst += "Copyright (C) 2003-2013, Linz, Austria\n";
	cst += "www.hotint.org\n";
	cst += "email to: support@hotint.org";

	atext->SetWindowText(cst);	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
	// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

void ErrorExit(LPCTSTR lpszFunction) 
{ 
    // Retrieve the system error message for the last-error code

    LPVOID lpMsgBuf;
    LPVOID lpDisplayBuf;
    DWORD dw = GetLastError(); 

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL );

    // Display the error message and exit the process

    lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT, 
        (lstrlen((LPCTSTR)lpMsgBuf)+lstrlen((LPCTSTR)lpszFunction)+40)*sizeof(TCHAR)); 
    StringCchPrintf((LPTSTR)lpDisplayBuf, 
        LocalSize(lpDisplayBuf) / sizeof(TCHAR),
        TEXT("%s failed with error %d: %s"), 
        lpszFunction, dw, lpMsgBuf); 
    MessageBox(NULL, (LPCTSTR)lpDisplayBuf, TEXT("Error"), MB_OK); 

    LocalFree(lpMsgBuf);
    LocalFree(lpDisplayBuf);
    ExitProcess(dw); 
}


/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DDlg dialog

CWCDriver3DDlg::CWCDriver3DDlg(CWnd* pParent /*=NULL*/) :
CDialog(CWCDriver3DDlg::IDD, pParent),
bClosingDialog(false),
pWCDI(NULL)
{
	//{{AFX_DATA_INIT(CWCDriver3DDlg)
	m_ConsoleText = _T("");
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
	hCursCross = AfxGetApp()->LoadStandardCursor(IDC_CROSS);
	hCursDragHor = AfxGetApp()->LoadStandardCursor(IDC_SIZEWE);

	//char buf[2000];
	//GetCurrentDirectory(2000,buf);
	//SetDllDirectoryA(buf);

	hWorkingDll = AfxLoadLibrary(WORKING_DLL_NAME);
	if(!hWorkingDll)
	{
		//ofstream ofile("check_working_moduleDLL_path.tmp");
		//ofile.close();
		ErrorExit((LPCTSTR)(CString("HOTINT: Cannot load dynamic-link library: ") + WORKING_DLL_NAME + "!\n"));
	}

	pExpFn = (PExpFn)(GetProcAddress(hWorkingDll,WORKING_OBJECT_CREATION_FUNCTION_NAME));
	if(!pExpFn)
		ErrorExit((LPCTSTR)(CString("Wrong dynamic-link library ") + WORKING_DLL_NAME + "!\n"));

	asv_models.SetLength(0); //$ DR 2013-05-21
}

CWCDriver3DDlg::~CWCDriver3DDlg()
{
	//$ YV 2013-01-04: leads to a crash
	//DeleteCriticalSection(&uses_append_text);
	//if crashes: put ! in if and put the destroy window operations after FreeLibrary
	if(::IsWindow(DialogDataManager.m_hWnd))
		DialogDataManager.DestroyWindow();

	if(::IsWindow(DialogFEDrawingOptions.m_hWnd))
		DialogFEDrawingOptions.DestroyWindow();

	if(::IsWindow(dialogComputationOutput.m_hWnd))
		dialogComputationOutput.DestroyWindow();

	if(::IsWindow(dialogViewingOptions.m_hWnd))
		dialogViewingOptions.DestroyWindow();

	//if(::IsWindow(dialogComputationSettings.m_hWnd))		//$ DR 2013-10-16 removed
	//	dialogComputationSettings.DestroyWindow();

	if(::IsWindow(dialogBodyJointOptions.m_hWnd))
		dialogBodyJointOptions.DestroyWindow();

	if(::IsWindow(computeEigenmodes.m_hWnd))
		computeEigenmodes.DestroyWindow();

	if(::IsWindow(dialogProgressbar.m_hWnd))
		dialogProgressbar.DestroyWindow();

	if(::IsWindow(Dialog2DView.m_hWnd))            //!AD: maybe put that in WCDriverFunction too, nested views seem to be a problem here
		Dialog2DView.DestroyWindow();

	CallWCDriverFunction(9);											 // close all instances of PlotToolDlgs

	for (int i=0; i < MAX_SENSOR_WATCH; i++)
	{
		if (sensorwatchctrl[i])
		{
			if (::IsWindow(sensorwatchctrl[i]->m_hWnd))
				sensorwatchctrl[i]->DestroyWindow();
			delete sensorwatchctrl[i];
		}
	}

	if(pWCDI)
		delete pWCDI;

	FreeLibrary(hWorkingDll);


}

void CWCDriver3DDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	/*
	//{{AFX_DATA_MAP(CWCDriver3DDlg)
	DDX_Text(pDX, IDC_EDIT_CONSOLE, m_ConsoleText);
	//}}AFX_DATA_MAP
	*/
}

//ON_WM_INITMENUPOPUP() has been added such that menu check and change status works in the dialog!!!
BEGIN_MESSAGE_MAP(CWCDriver3DDlg, CDialog)
	//{{AFX_MSG_MAP(CWCDriver3DDlg)
	ON_WM_INITMENUPOPUP()
	ON_WM_SYSCOMMAND()
	ON_WM_QUERYDRAGICON()
	ON_WM_SIZE()
	ON_BN_CLICKED(IDC_BUTTON_GO, OnButtonGo)
	ON_BN_CLICKED(IDC_BUTTON_PAUSE, OnButtonPause)
	ON_WM_PAINT()
	ON_WM_LBUTTONDOWN()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONUP()
	ON_WM_CLOSE()
	ON_MESSAGE(WM_UPDATE_TEXT,OnUpdateText)
	ON_MESSAGE(WM_START_COMPUTATION, OnButtonGoMessage) //added by JG manually
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON_GENERALOPTIONS, OnBnClickedButtonGeneraloptions)
	ON_COMMAND(ID_RESULTS_DATAMANAGER, OnResultsDatamanager)
	ON_COMMAND(ID_RESULTS_VIEWOUTPUT, OnResultsViewoutput)
	ON_COMMAND(ID_COMPUTATION_RESET, OnComputationReset)
	ON_COMMAND(ID_COMPUTATION_PRINTTIMINGS, OnComputationPrinttimings)
	ON_COMMAND(ID_RESULTS_SHOWOUTPUTSTATICTEXT, OnResultsShowoutputstatictext)
	ON_COMMAND(ID_OPTIONS_SAVECONFIGURATION, OnOptionsSaveconfiguration)
	ON_COMMAND(ID_OPTIONS_SAVESOLVEROPTIONS, OnOptionsSaveSolverOptions)
	ON_COMMAND(ID_OPTIONS_SAVEHOTINTOPTIONS, OnOptionsSaveHotintOptions)
	ON_COMMAND(ID_RESULTS_ENABLEOUTPUT, OnResultsEnableoutput)
	ON_WM_MOVE()
	ON_COMMAND(ID_VIEW_FE_DRAWINGOPTIONS, OnViewFeDrawingoptions)
	ON_COMMAND(ID_VIEW_VIEWINGOPTIONS, OnViewViewingoptions)
	//	ON_WM_KEYUP()
	//ON_COMMAND(ID_COMPUTATION_COMPUTATIONOPTIONS, OnComputationComputationoptions)	//$ DR 2013-10-16 removed
  	ON_COMMAND(ID_COMPUTATION_COMPUTEEIGENMODES, OnComputationComputeEigenmodes)
	ON_COMMAND(ID_VIEW_RIGIDBODY, OnViewRigidbodyJointOptions)
	ON_COMMAND(ID_EDIT_EDITBODY, OnEditEditbody)
	//ON_COMMAND(ID_COMPUTATION_ASSEMBLESYSTEM, OnComputationAssemblesystem)
	ON_COMMAND(ID_VIEW_OPENGLOPTIONS, OnViewOpengloptions)
	//ON_COMMAND(ID_RIGIDBODIES_ADDMASS3D, OnRigidbodiesAddmass3d)
	//ON_COMMAND(ID_COMPUTATION_ASSIGNINITIALVECTOR, OnComputationAssigninitialvector)
	//ON_COMMAND(ID_RIGIDBODY_ADDRIGID3D, OnRigidbodyAddrigid3d)
	//ON_COMMAND(ID_3DJOINT_SPHERICALJOINT, On3djointSphericaljoint)
	//ON_COMMAND(ID_3DJOINT_REVOLUTEJOINT, On3djointRevolutejoint)
	ON_COMMAND(ID_FILE_EXIT, OnFileExit)
	ON_COMMAND(ID_FILE_NEWMBS, OnFileNewmbs)
	//ON_COMMAND(ID_ADDLOAD_GCLOAD, OnAddloadGcload)
	//ON_COMMAND(ID_ADDLOAD_BODYLOAD, OnAddloadBodyload)
	//ON_COMMAND(ID_ADDLOAD_FORCEVECTOR2D, OnAddloadForcevector2d)
	//ON_COMMAND(ID_ADDLOAD_MOMENT2D, OnAddloadMoment2d)
	//ON_COMMAND(ID_ADDLOAD_FORCEVECTOR3D, OnAddloadForcevector3d)
	//ON_COMMAND(ID_ADDLOAD_MOMENTVECTOR3D, OnAddloadMomentvector3d)
//	ON_COMMAND(ID_JOINTS_ADDJOINT, OnJointsAddjoint)
	ON_COMMAND(ID_FILE_SAVEMBS, OnFileSavembs)
	ON_COMMAND(ID_EDIT_EDITSENSOR, OnEditEditsensor)
	ON_COMMAND(ID_FILE_OPENMBS, OnFileOpenmbs)
	ON_COMMAND(ID_FILE_SAVEMBSAS, OnFileSavembsas)
	ON_COMMAND(ID_SPECIAL_ADDSENSOR, OnSpecialAddsensor)
	ON_COMMAND(ID_RESULTS_PLOTSENSOR, OnResultsPlotsensor)
	ON_COMMAND(ID_SPECIAL_ADDBODY, OnSpecialAddbody)
	ON_COMMAND(ID_RESULTS_SENSORWATCH, OnResultsSensorwatch)
	ON_COMMAND(ID_FILE_RECENTFILE1, OnFileRecentfile1)
	ON_BN_CLICKED(IDC_BUTTON_RECENT_FILE, OnBnClickedButtonRecentFile)
	ON_COMMAND(ID_SPECIAL_ADDGEOMELEMENT, OnSpecialAddgeomelement)
	ON_COMMAND(ID_EDIT_EDITGEOMELEMENT, OnEditEditgeomelement)
	//ON_COMMAND(ID_ADDCONNECTOR_JOINT, OnAddconnectorJoint2D)
	//ON_COMMAND(ID_ADDCONNECTOR_3DJOINT, OnAddConnectorJoint3D)
	//ON_COMMAND(ID_ADDCONNECTOR_SPECIAL, OnAddconnectorSpecial)
	ON_COMMAND(ID_SYSTEM_SHOWSYSTEMPROPERTIES, OnSystemShowsystemproperties)
	ON_COMMAND(ID_ADDOBJECT_ADDNODE, OnAddobjectAddnode)
	ON_COMMAND(ID_EDIT_EDITNODE, OnEditEditnode)
	ON_COMMAND(ID_SYSTEM_CHECKSYSTEM, OnSystemChecksystem)
	ON_COMMAND(ID__ABOUT, OnHelpAbout)
	ON_COMMAND(ID_VIEW_DEFAULTVIEW1, OnViewDefaultview1)
	ON_COMMAND(ID_VIEW_XY, OnViewXY)
	ON_COMMAND(ID_VIEW_XZ, OnViewXZ)
	ON_COMMAND(ID_VIEW_YZ, OnViewYZ)
	ON_WM_MOUSEWHEEL()
	ON_COMMAND(ID_COMPUTATION_LOADINITIALVECTOR, OnComputationLoadinitialvector)
	ON_COMMAND(ID_COMPUTATION_STORESOLUTIONVECTOR, OnComputationStoresolutionvector)
	//ON_COMMAND(ID_EDIT_CHANGEELEMENTNUMBERS, OnEditChangeelementnumbers)
	ON_COMMAND(ID__HELP, OnMenuHelp)
	ON_COMMAND(ID_VIEW_ROBOTOPTIONS, OnViewRobotoptions)
	ON_COMMAND(ID_EDIT_EDITMATERIAL, OnEditEditmaterial)
	ON_COMMAND(ID_VIEW_SAVECONFIGURATION, OnOptionsSaveconfiguration)
	ON_COMMAND(ID_SYSTEM_EDITMODELPARAMETERS, &CWCDriver3DDlg::OnSystemEditmodelparameters)
	ON_COMMAND(ID_SYSTEM_LOADMODELPARAMETER, &CWCDriver3DDlg::OnSystemLoadmodelparameter)
	ON_COMMAND(ID_SYSTEM_SAVEMODELDATA, &CWCDriver3DDlg::OnSystemSavemodeldata)
	ON_COMMAND(ID_VIEW_EDITALLOPTIONS, &CWCDriver3DDlg::OnViewEditalloptions)
	ON_WM_TIMER()
	ON_COMMAND(ID_RESULTS_PLOT2SENSORSXY, &CWCDriver3DDlg::OnResultsPlot2sensorsxy)
	//ON_COMMAND(ID_RESULTS_PLOTNSENSORS, &CWCDriver3DDlg::OnResultsPlotnsensors)
	ON_COMMAND(ID_VIEW_SHOWPROGRESSBAR, &CWCDriver3DDlg::OnViewShowprogressbar)
	ON_COMMAND(ID_RESULTS_PLOTTOOLDIALOG, &CWCDriver3DDlg::OnResultsPlottooldialog)
	ON_COMMAND(ID_FILE_SELECTMODEL, &CWCDriver3DDlg::OnFileSelectmodel)
	ON_COMMAND(ID_EDIT_EDITSOLVEROPTIONS, &CWCDriver3DDlg::OnEditEditsolveroptions)
	ON_COMMAND(ID_EDIT_EDITHOTINTOPTIONS, &CWCDriver3DDlg::OnEditEdithotintoptions)
	ON_COMMAND(ID_COMPUTATION_EDITCOMPUTATIONSTEPS, &CWCDriver3DDlg::OnEditComputationSteps)
	ON_COMMAND(ID_EDIT_EDITLOAD, &CWCDriver3DDlg::OnEditEditLoad)
	ON_COMMAND(ID_ADDOBJECT_ADDLOAD, &CWCDriver3DDlg::OnAddload)
	ON_COMMAND(ID_VIEW_SHOW_IO_BLOCKS_WINDOW, &CWCDriver3DDlg::OnViewShowContolWindow)
	ON_COMMAND(ID_ADDOBJECT_ADDMATERIAL, &CWCDriver3DDlg::OnAddobjectAddmaterial)
	ON_COMMAND(ID_ADDOBJECT_ADDBEAM3DPROPERTIES, &CWCDriver3DDlg::OnAddbeam3dproperties)
	ON_WM_KEYDOWN()
	ON_COMMAND(ID_ADDCONNECTOR_KINEMATICPAIRS, &CWCDriver3DDlg::OnAddconnectorKinematicpairs)
	ON_COMMAND(ID_ADDCONNECTOR_CONTROLELEMENT, &CWCDriver3DDlg::OnAddconnectorControlelement)
	ON_BN_CLICKED(IDC_BUTTON_SAVE_MBS, &CWCDriver3DDlg::OnBnClickedButtonSaveMbs)
	ON_COMMAND(ID_EDIT_DELETELOAD32963, &CWCDriver3DDlg::OnEditDeleteLoad)
	ON_COMMAND(ID_EDIT_DELETEGEOMELEMENT32964, &CWCDriver3DDlg::OnEditDeleteGeomelement)
	ON_COMMAND(ID_EDIT_DELETESENSOR32965, &CWCDriver3DDlg::OnEditDeleteSensor)
	ON_COMMAND(ID_EDIT_DELETENODE32966, &CWCDriver3DDlg::OnEditDeleteNode)
	ON_COMMAND(ID_EDIT_DELETEMATERIAL32967, &CWCDriver3DDlg::OnEditDeleteMaterial)
	ON_COMMAND(ID_EDIT_DELETEELEMENT32968, &CWCDriver3DDlg::OnEditDeleteElement)
	ON_COMMAND(ID_EDIT_UNDO32969, &CWCDriver3DDlg::OnEditUndo)
	ON_COMMAND(ID_SYSTEM_RUNMACRO, &CWCDriver3DDlg::OnSystemRunmacro)
	END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DDlg message handlers

BOOL CWCDriver3DDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	fileopened = 0;
	model_modified = 0;
	enableOutputText = 1;


	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	pWCDI = (*pExpFn)();
	pWCDI->SetUserInterface(this);
	GLDrawWnd.SetWCDI(pWCDI);

	GetClientRect(&PictureRect);
	//CRect EditRect;
	//GetDlgItem(IDC_EDIT_CONSOLE)->GetWindowRect(&EditRect);
	CRect GoButtonRect;
	GetDlgItem(IDC_BUTTON_GO)->GetWindowRect(&GoButtonRect);
	PictureRect.left = /*EditRect.Width() + */  NEUTRAL_ZONE_WIDTH;
	PictureRect.top = GoButtonRect.Height();

	//PositionButtons();

	DialogDataManager.SetWCDI(pWCDI);
	DialogDataManager.SetGLDrawWnd(&GLDrawWnd);

	DialogFEDrawingOptions.SetWCDI(pWCDI);
	DialogFEDrawingOptions.SetGLDrawWnd(&GLDrawWnd);

	dialogViewingOptions.SetWCDI(pWCDI);
	dialogViewingOptions.SetGLDrawWnd(&GLDrawWnd);

	//dialogComputationSettings.SetWCDI(pWCDI);		//$ DR 2013-10-16 removed
	//dialogComputationSettings.SetGLDrawWnd(&GLDrawWnd);

	dialogBodyJointOptions.SetWCDI(pWCDI);
	dialogBodyJointOptions.SetGLDrawWnd(&GLDrawWnd);

	dialogOpenGLOptions.SetWCDI(pWCDI);
	dialogOpenGLOptions.SetGLDrawWnd(&GLDrawWnd);

	computeEigenmodes.SetWCDI(pWCDI);
	computeEigenmodes.SetGLDrawWnd(&GLDrawWnd);
	computeEigenmodes.SetWCDDlg(this);

	dialogProgressbar.SetWCDI(pWCDI);
	dialogProgressbar.SetGLDrawWnd(&GLDrawWnd);

	Dialog2DView.SetWCDI(pWCDI);
	Dialog2DView.SetWCDDlg(this);
	//Dialog2DView.SetGLDrawWnd(&GLDrawWnd);

	//appPlotTool = AfxBeginThread((CRuntimeClass*)RUNTIME_CLASS(CPlotToolApp));
 // dialogPlotTool.SetWCDI(pWCDI);
	//dialogPlotTool.SetGLDrawWnd(&GLDrawWnd);

	//do this after setting user interface
	ReadConfigFile();
	pWCDI->InitializeAfterConfigLoaded();
	DisplaySelectedModelName();

	GLDrawWnd.Create(
		NULL,NULL,
		WS_CHILD | WS_BORDER | WS_VISIBLE,
		PictureRect,this,0
		);

  for (int i=0; i < MAX_SENSOR_WATCH; i++)
	{
		sensorwatchctrl[i] = 0;
	}

	if (strlen(pWCDI->GetTOption(101)) > 0)
	{
		CMenu* mmenu = GetMenu();
		mmenu->ModifyMenu(ID_FILE_RECENTFILE1, MF_BYCOMMAND, ID_FILE_RECENTFILE1, CString("&Recent: ")+CString(pWCDI->GetTOption(94)));
		mmenu->EnableMenuItem(ID_FILE_RECENTFILE1, MF_BYCOMMAND | MF_ENABLED); //MF_DISABLED | MF_GRAYED
	}
	else
	{
		CMenu* mmenu = GetMenu();
		mmenu->EnableMenuItem(ID_FILE_RECENTFILE1, MF_BYCOMMAND | MF_DISABLED | MF_GRAYED); //
	}


////#ifndef KB_INCLUDES
//	{
//	//robot options: --> deprecated ==> always disabled
//	CMenu* mmenu = GetMenu();
//	//mmenu->EnableMenuItem(ID_VIEW_ROBOTOPTIONS, MF_BYCOMMAND | MF_DISABLED | MF_GRAYED); //
//	mmenu->RemoveMenu(ID_VIEW_ROBOTOPTIONS, MF_BYCOMMAND);
//	}

//#endif


	RemoveExperimentalMenuItems();

	//start automatically:
	int autostart_flag;
	pWCDI->MBS_EDC_TreeGetInt(autostart_flag, "GeneralOptions.Application.start_computation_automatically");
	if (autostart_flag)
	{
		PostMessage(WM_START_COMPUTATION);
	}

	// minimize window if IOption(221) = "GeneralOptions.Application.show_hotint_window" is set to zero
	int showwindow_flag;
	pWCDI->MBS_EDC_TreeGetInt(showwindow_flag, "GeneralOptions.Application.show_hotint_window");
	if (!showwindow_flag)
		ShowWindow(SW_SHOWMINNOACTIVE);

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CWCDriver3DDlg::RemoveExperimentalMenuItems()
{
#ifdef __EXCLUDE_EXPERIMENTAL_MENU_ITEMS__
	CMenu* mmenu = GetMenu();

	//// ?
	mmenu->RemoveMenu(ID__HELP, MF_BYCOMMAND);

	//// Computation
	mmenu->RemoveMenu(ID_COMPUTATION_RESET, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_START, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_STOP, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_PAUSE, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_LOADINITIALVECTOR, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_STORESOLUTIONVECTOR, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_PRINTTIMINGS, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_COMPUTATION_COMPUTEEIGENMODES, MF_BYCOMMAND);

	//// System
	mmenu->RemoveMenu(5, MF_BYPOSITION);

	//// Add Object
	mmenu->RemoveMenu(4, MF_BYPOSITION);

	//// Edit
	mmenu->RemoveMenu(3, MF_BYPOSITION);

	//// View
	mmenu->RemoveMenu(ID_VIEW_SHOWPROGRESSBAR, MF_BYCOMMAND);

	//// File
	mmenu->RemoveMenu(ID_FILE_NEWMBS, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_FILE_OPENMBS, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_FILE_SAVEMBS, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_FILE_SAVEMBSAS, MF_BYCOMMAND);
	mmenu->RemoveMenu(ID_FILE_RECENTFILE1, MF_BYCOMMAND);

	//// finally, remove menu bars
	//mmenu->GetSubMenu(4)->RemoveMenu(2, MF_BYPOSITION);
	//mmenu->GetSubMenu(2)->RemoveMenu(4, MF_BYPOSITION);
	//mmenu->GetSubMenu(2)->RemoveMenu(3, MF_BYPOSITION);
	//mmenu->GetSubMenu(2)->RemoveMenu(2, MF_BYPOSITION);
	//mmenu->GetSubMenu(1)->RemoveMenu(0, MF_BYPOSITION);
	//mmenu->GetSubMenu(0)->RemoveMenu(4, MF_BYPOSITION);
	//mmenu->GetSubMenu(0)->RemoveMenu(4, MF_BYPOSITION);
	//mmenu->GetSubMenu(0)->RemoveMenu(2, MF_BYPOSITION);
	//mmenu->GetSubMenu(0)->RemoveMenu(1, MF_BYPOSITION);

	// remove edit-dialog which shows modelname
	GetDlgItem(IDC_EDIT_MODELNAME)->ModifyStyle(WS_VISIBLE,0,0);
#endif
}

//add this such that menu check, change text, etc. works in dialog with menu!!!
void CWCDriver3DDlg::OnInitMenuPopup(CMenu *pPopupMenu, UINT nIndex,BOOL bSysMenu)
{
	ASSERT(pPopupMenu != NULL);
	// Check the enabled state of various menu items.

	CCmdUI state;
	state.m_pMenu = pPopupMenu;
	ASSERT(state.m_pOther == NULL);
	ASSERT(state.m_pParentMenu == NULL);

	// Determine if menu is popup in top-level menu and set m_pOther to
	// it if so (m_pParentMenu == NULL indicates that it is secondary popup).
	HMENU hParentMenu;
	if (AfxGetThreadState()->m_hTrackingMenu == pPopupMenu->m_hMenu)
		state.m_pParentMenu = pPopupMenu;    // Parent == child for tracking popup.
	else if ((hParentMenu = ::GetMenu(m_hWnd)) != NULL)
	{
		CWnd* pParent = this;
		// Child windows don't have menus--need to go to the top!
		if (pParent != NULL &&
			(hParentMenu = ::GetMenu(pParent->m_hWnd)) != NULL)
		{
			int nIndexMax = ::GetMenuItemCount(hParentMenu);
			for (int nIndex = 0; nIndex < nIndexMax; nIndex++)
			{
				if (::GetSubMenu(hParentMenu, nIndex) == pPopupMenu->m_hMenu)
				{
					// When popup is found, m_pParentMenu is containing menu.
					state.m_pParentMenu = CMenu::FromHandle(hParentMenu);
					break;
				}
			}
		}
	}

	state.m_nIndexMax = pPopupMenu->GetMenuItemCount();
	for (state.m_nIndex = 0; state.m_nIndex < state.m_nIndexMax;
		state.m_nIndex++)
	{
		state.m_nID = pPopupMenu->GetMenuItemID(state.m_nIndex);
		if (state.m_nID == 0)
			continue; // Menu separator or invalid cmd - ignore it.

		ASSERT(state.m_pOther == NULL);
		ASSERT(state.m_pMenu != NULL);
		if (state.m_nID == (UINT)-1)
		{
			// Possibly a popup menu, route to first item of that popup.
			state.m_pSubMenu = pPopupMenu->GetSubMenu(state.m_nIndex);
			if (state.m_pSubMenu == NULL ||
				(state.m_nID = state.m_pSubMenu->GetMenuItemID(0)) == 0 ||
				state.m_nID == (UINT)-1)
			{
				continue;       // First item of popup can't be routed to.
			}
			state.DoUpdate(this, TRUE);   // Popups are never auto disabled.
		}
		else
		{
			// Normal menu item.
			// Auto enable/disable if frame window has m_bAutoMenuEnable
			// set and command is _not_ a system command.
			state.m_pSubMenu = NULL;
			state.DoUpdate(this, FALSE);
		}

		// Adjust for menu deletions and additions.
		UINT nCount = pPopupMenu->GetMenuItemCount();
		if (nCount < state.m_nIndexMax)
		{
			state.m_nIndex -= (state.m_nIndexMax - nCount);
			while (state.m_nIndex < nCount &&
				pPopupMenu->GetMenuItemID(state.m_nIndex) == state.m_nID)
			{
				state.m_nIndex++;
			}
		}
		state.m_nIndexMax = nCount;
	}
}


void CWCDriver3DDlg::PositionButtons()
{
	CRect Rect;
	CWnd * pB = GetDlgItem(IDC_BUTTON_GO);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += PictureRect.left - Rect.left;
	Rect.left = PictureRect.left;
	pB->MoveWindow(Rect,FALSE);

	int right;
	/*
	right = Rect.right;
	pB = GetDlgItem(IDC_BUTTON_DATA_MANAGER);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += right + BUTTONS_OFFSET - Rect.left;
	Rect.left = right + BUTTONS_OFFSET;
	pB->MoveWindow(Rect,FALSE);
	*/
	/*
	right = Rect.right;
	pB = GetDlgItem(IDC_BUTTON_READ_TEXT);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += right + BUTTONS_OFFSET - Rect.left;
	Rect.left = right + BUTTONS_OFFSET;
	pB->MoveWindow(Rect,FALSE);
	*/
	/*
	right = Rect.right;
	pB = GetDlgItem(IDC_BUTTON_SAVE_CONFIG);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += right + BUTTONS_OFFSET - Rect.left;
	Rect.left = right + BUTTONS_OFFSET;
	pB->MoveWindow(Rect,FALSE);
	*/


	right = Rect.right;
	pB = GetDlgItem(IDC_BUTTON_GENERALOPTIONS);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += right + BUTTONS_OFFSET - Rect.left;
	Rect.left = right + BUTTONS_OFFSET;
	pB->MoveWindow(Rect,FALSE);

	/*
	right = Rect.right;
	pB = GetDlgItem(IDC_CHECK_LIGHTING);
	pB->GetWindowRect(&Rect);
	ScreenToClient(Rect);
	Rect.right += right + BUTTONS_OFFSET - Rect.left;
	Rect.left = right + BUTTONS_OFFSET;
	pB->MoveWindow(Rect,FALSE);
	*/
}

LRESULT CWCDriver3DDlg::OnButtonGoMessage(WPARAM, LPARAM)
{
	OnButtonGo();
	return 0;
}


// OLD VERSION:
/*
LRESULT CWCDriver3DDlg::OnUpdateText(WPARAM, LPARAM)
{
	if(!DialogDataManager.RetrievingData())
	{
		UpdateData(FALSE);
		
		//CEdit * pEdit = (CEdit*)GetDlgItem(IDC_EDIT_CONSOLE);
		//pEdit->SetSel(100000000,100000000,FALSE);
		
		if(::IsWindow(dialogComputationOutput.m_hWnd))
		{
			dialogComputationOutput.UpdateText(m_ConsoleText);
		}
	}
	return 0;
}*/

double global_last_printtext_time = 0;
double GetWallClockTime()
{
	timeb tb;
	tb.time = 0;
	tb.millitm = 0;
	ftime(&tb);
	return (tb.time+(double)tb.millitm*0.001);	
}
// NEW VERSION (AD) calls Append, ReplaceLastLine or Update(full)
LRESULT CWCDriver3DDlg::OnUpdateText(WPARAM, LPARAM)
{
	if(!DialogDataManager.RetrievingData())
	{
		MSG msg;
		while(::PeekMessage(&msg, m_hWnd, WM_UPDATE_TEXT, WM_UPDATE_TEXT, PM_REMOVE)); // remove all pending update text messages

		UpdateData(FALSE);

		if(::IsWindow(dialogComputationOutput.m_hWnd))
		{
			if(update_text_action_flag == 1)	
			{
				dialogComputationOutput.AppendText(m_ConsoleText_append);
			}
			else if(update_text_action_flag == 2)
			{
				dialogComputationOutput.ReplaceLastLine(m_ConsoleText_append);
			}
			else if(update_text_action_flag == 0) 
			{
				dialogComputationOutput.UpdateText(m_ConsoleText);
			}
			else
			{
				int maxlength = m_ConsoleText.GetLength();

				if (pWCDI->GetIOption(152) != -1)
				{
					maxlength = min(pWCDI->GetIOption(152),maxlength);
					CString str = m_ConsoleText.Right(maxlength);
					dialogComputationOutput.UpdateText(str);
				}
				else
				{
					dialogComputationOutput.UpdateText(m_ConsoleText);
				}
			}
/*
			global_last_printtext_time = GetWallClockTime();	
			SetTextRefreshTimer(0);       // reset timer
// critical section - accessing variable with append_text
EnterCriticalSection(&uses_append_text);
			m_ConsoleText_buffer = m_ConsoleText_buffer.Mid(m_ConsoleText_append.GetLength());
			m_ConsoleText_append.Empty(); // delete append buffer after text was updated
LeaveCriticalSection(&uses_append_text);
*/
		}
	}
	return 0;
}

void CWCDriver3DDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
/*	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else*/
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

int hotint_startup_message_written = 0;
void CWCDriver3DDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}

#ifdef HOTINT_PUBLIC_DOMAIN
	if (!hotint_startup_message_written && pWCDI->GetIOption(141))
	{
		CString message = "";
		message += "Welcome to HOTINT V0.900 beta!\n";
		message += "This is a beta version and not all features have been yet fully tested!\n";
		message += "While the programmer acted with high caution and a lot of verification has been\n";
		message += "done already, it is recommended to verify the results!\n";
		message += "This program is distributed in the hope that it will be useful,\n"; 
		message += "but WITHOUT ANY WARRANTY; without even the implied warranty of \n";
		message += "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE!\n\n";

		message += "******   Bug reports and comments are welcome: jg@jku.at   ******";

		AfxMessageBox(message);
	}
	hotint_startup_message_written = 1;
#endif

}


// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CWCDriver3DDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CWCDriver3DDlg::OnSize(UINT nType, int cx, int cy) 
{
	if(cx != 0 && cy != 0)
	{
		int vertical_shift = PictureRect.bottom-cy;

		PictureRect.right = cx;
		PictureRect.bottom = cy;

		if(::IsWindow(GLDrawWnd.m_hWnd))
			GLDrawWnd.MoveWindow(PictureRect,FALSE);

		/*
		CRect r;

		CWnd * pEditConsole = GetDlgItem(IDC_EDIT_CONSOLE);
		if(pEditConsole)
		{
		pEditConsole->GetWindowRect(&r);
		ScreenToClient(&r);
		pEditConsole->MoveWindow(r.left,r.top,r.Width(),r.Height()-vertical_shift,FALSE);
		}
		*/

		RedrawWindow();
	}

	CDialog::OnSize(nType, cx, cy);
}


void CWCDriver3DDlg::OnMove(int x, int y)
{
	if(::IsWindow(dialogComputationOutput.m_hWnd))
	{
		CRect r, rw;
		GetWindowRect(&rw);
		dialogComputationOutput.GetWindowRect(&r);
		dialogComputationOutput.MoveWindow(rw.left-r.Width(), rw.top, r.Width(), r.Height(), TRUE); //no Repaint
	}

	__super::OnMove(x, y);

	// TODO: Fügen Sie hier Ihren Meldungsbehandlungscode ein.
}


void CWCDriver3DDlg::InstantMessageText(const char * pStr)
{
	AfxMessageBox(pStr);
}

void CWCDriver3DDlg::StatusText(const char * pStr)
{
	CString text = CString("HOTINT");
	if (fileopened)
	{
		if (strlen(pWCDI->GetTOption(85)) != 0)
		{
			text += " - ";
			text += pWCDI->GetTOption(85);
		}
		else
		{
			text += "Multibody simulation code";
		}
	}
	int len = strlen(pStr);
	if (len != 0)	text = CString(pStr) + CString(" - ") + text;

	char str[512];
	len = strlen(text);

	if (len > 511) text.SetAt(511, (char)0);
	strcpy(str, text);

	SendMessage(WM_SETTEXT,NULL,(LPARAM)str);
}

void CWCDriver3DDlg::DisplaySelectedModelName()
{
// version: display in Main Window 
	CString text = CString("HOTINT --- current Model: ") + CString(pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name"));
	this->SetWindowText(text);

// version: display in Iconbar - disabled as of 2013-09-30	
	//CString text = CString("Model: ") + CString(pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name"));
	//CEdit* eb = (CEdit*) GetDlgItem(IDC_EDIT_MODELNAME);
	//eb->SetWindowText(text);
}

void CWCDriver3DDlg::AddRecentFile(const char * filename) 	//$ DR 2012-12-06 add "filename" to list of recent files and resort this list
{

	TArrayDynamic<mystr> recent_files;
	for(int i=0; i<=4; i++)
	{
		recent_files.Add(pWCDI->GetTOption(113+i));
	}

	int pos = recent_files.Find(filename);
	if(pos==1)
	{
		return;		// is already at first position
	}
	else		// put the actual file at first position and add the others in old order
	{
		pWCDI->SetTOption(113,filename);
		int j;
		for(int i=1; i<=4; i++)
		{
			j=i;
			if(pos==i) {j++;}	// if already in list
			pWCDI->SetTOption(113+i,recent_files(j));
		}
		ElementDataContainer edc;
		pWCDI->CallCompFunction(210, 0, 0, &edc);		//$ DR 2012-12-17 update EDC of MBS
	}
}

void CWCDriver3DDlg::AddText(const char * pStr)
{
	if (!enableOutputText) return;

	update_text_action_flag = 3;
	const static char _EOL[] = {13,10}; //{13,13,10}
	const static CString EOL = _EOL;
	const static char _CARRET[] = {13,10}; //{13,13,10}
	const static CString CARRET = _CARRET;

	CString str(pStr);
	int nEOL = str.Find('\n',0);
	while(nEOL != -1)
	{
		m_ConsoleText += str.Left(nEOL) + EOL;
		str = str.Mid(nEOL+1);
		nEOL = str.Find('\n',0);
	}

	nEOL = str.Find('\r',0);
	if(nEOL != -1)
	{
		str = str.Mid(nEOL+1);
		int off = 1;
		int n = m_ConsoleText.Find((char)10, m_ConsoleText.GetLength()-off);
		while ( n == -1 && off < m_ConsoleText.GetLength())
		{
			off += 1;
			n = m_ConsoleText.Find((char)10, m_ConsoleText.GetLength()-off);
		}

		if (n != -1)
		{
			m_ConsoleText = m_ConsoleText.Left(n+1);
		}
	}
	m_ConsoleText += str;

	const double maxtime = 0;
	double acttime = 0;
	{
		timeb tb;
		tb.time = 0;
		tb.millitm = 0;
		ftime(&tb);
		//cout << "time=" << tb.time << ", millisec=" << tb.millitm << endl;
		acttime = tb.time+(double)tb.millitm*0.001;	
	}

	if (acttime - global_last_printtext_time > maxtime)
	{
		SendMessage(WM_UPDATE_TEXT);		// we can not just call UpdateData(FALSE) due to multithreading
	}
	else
	{
		SendMessage(WM_UPDATE_TEXT);		// we can not just call UpdateData(FALSE) due to multithreading
	}
}

// create directory (only, if it does not exist).
// returns 0 in case of error (wrong path), 1 in case it already existed, and 2 if created
int CWCDriver3DDlg::CreateMissingDirectory(const char * str)
{
	int rv = 2;

	if (!CreateDirectory(str, NULL))
	{
		DWORD err = GetLastError();
		if (err == ERROR_ALREADY_EXISTS)
		{
			rv = 1;
		}
		else   // err == ERROR_PATH_NOT_FOUND
		{
			rv = 0;
		}
	}
	return rv;
}

void CWCDriver3DDlg::SetPlotToolRedrawTimer(int flag_onoff, int delay_ms)
{
	if (flag_onoff)
	{
		int suc = SetTimer(PLOTTOOL_REFRESH_TIMER_ID, (UINT) delay_ms, NULL);
		if (!suc)
			::MessageBox(this->GetSafeHwnd(),"Unable to start timer","PLOTTOOL_REFRESH_TIMER_ID",MB_OK|MB_SYSTEMMODAL);
	}
	else
	{
		KillTimer(PLOTTOOL_REFRESH_TIMER_ID);
	}
}

void CWCDriver3DDlg::OnTimer(UINT_PTR nIDEvent)
{
	if (nIDEvent == PLOTTOOL_REFRESH_TIMER_ID)
	{
// check for pending PlotUpdateMessage
  	MSG msg;
		if (::PeekMessageA(&msg, m_hWnd, WM_UPDATE_PLOT, WM_UPDATE_PLOT, 0))
		{
			SendMessage(WM_UPDATE_PLOT);
		}
	}

	CWnd::OnTimer(nIDEvent);
}

void CWCDriver3DDlg::ResultsUpdated(int flag)
{
	if (flag==1 || flag==0)
		GLDrawWnd.SendMessage(WM_REDRAW);

	if((flag==2 || flag==0)&& GetWCDInterface()->StoreResultsIsOn())
	{
		DataStorage ds; double m_TimePoint;
		GetWCDInterface()->StoreResults(ds, m_TimePoint);
		DialogDataManager.AddEntry(ds, m_TimePoint);

		for (int i=0; i < MAX_SENSOR_WATCH; i++)
		{
			if (sensorwatchctrl[i]) 
			{
				if (::IsWindow(sensorwatchctrl[i]->m_hWnd))
				{
					ElementDataContainer edc;
					pWCDI->CallCompFunction(204, 0, sensorwatchctrl[i]->GetSensorNumber(), &edc);
					if (edc.Length() > 0) 
					{
						double v = edc.Get(1).GetDouble();
						char str[32];
						sprintf_s(str, "%.17g",v);
						sensorwatchctrl[i]->SetSensorText(str);
					}
				}
			}
		}
	}
}

void CWCDriver3DDlg::OnLButtonDown(UINT nFlags, CPoint point) 
{
	CDialog::OnLButtonDown(nFlags, point);
}


void CWCDriver3DDlg::OnMouseMove(UINT nFlags, CPoint point) 
{

	CDialog::OnMouseMove(nFlags, point);
}

void CWCDriver3DDlg::OnLButtonUp(UINT nFlags, CPoint point) 
{
	CDialog::OnLButtonUp(nFlags, point);
}


BOOL CWCDriver3DDlg::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	GLDrawWnd.OnMouseWheelGLDW(nFlags, zDelta, pt); //windows function does not pass the message ...

	return CDialog::OnMouseWheel(nFlags, zDelta, pt);
}


UINT ThreadControlFunc(LPVOID pVOID)
{
	CWCDriver3DDlg * pDlg = (CWCDriver3DDlg*)pVOID;
	pDlg->GetWCDInterface()->Go(pDlg);

	return 0;
}
/*
//parallel threads work well for parallelism!!!
int cnt;
UINT ThreadControlFunc2(LPVOID pVOID)
{
CWCDriver3DDlg * pDlg = (CWCDriver3DDlg*)pVOID;
//pDlg->GetWCDInterface()->Go(pDlg);

double x=1;
for (int i=1; i <= 100000; i++) 
{
for (int j=1; j <= 1; j++)
{
x+=sin((i+j)/100000.);
}
}
char str[200];
sprintf_s(str, "ready: %g",x);
//pDlg->AddText(str);
cnt--;
sprintf_s(str, "R%d %1.1f",cnt,(x+1.)/x);
pDlg->GetDlgItem(IDC_BUTTON_GO)->SetWindowText(str);

return 0;
}
*/

void CWCDriver3DDlg::OnButtonGo() 
{

	if(!pWCDI->IsComputationInProgress())
	{
		GetDlgItem(IDC_BUTTON_GO)->SetWindowText("Stop");
		//GetDlgItem(IDC_BUTTON_PAUSE)->SetForegroundWindow();
		CWinThread * pWinThread  = AfxBeginThread(ThreadControlFunc,this);
		SetThreadPriority(pWinThread->m_hThread,THREAD_PRIORITY_ABOVE_NORMAL);
		//SetThreadPriority(pWinThread->m_hThread,THREAD_PRIORITY_BELOW_NORMAL);
		/*
		cnt = 1000;
		for (int i=1; i <= 1000; i++)
		{
		CWinThread * pWinThread  = AfxBeginThread(ThreadControlFunc2,this);
		SetThreadPriority(pWinThread->m_hThread,THREAD_PRIORITY_BELOW_NORMAL);
		}*/
	}
	else
	{
		GetDlgItem(IDC_BUTTON_GO)->SetWindowText("Stopping..");
		pWCDI->StopWhenPossible();
	}
}

void CWCDriver3DDlg::FinishedComputation()
{
	GetDlgItem(IDC_BUTTON_GO)->SetWindowText("Restart");
}

void CWCDriver3DDlg::OnButtonPause() 
{
	if (pWCDI->IsPaused())
	{
		//set scroll position to last time point
		DialogDataManager.SetScrollPosToLastTimePoint();
		
		//resume computation
		GetDlgItem(IDC_BUTTON_GO)->EnableWindow(TRUE);
		GetDlgItem(IDC_BUTTON_PAUSE)->SetWindowText("Pause");	
		pWCDI->Resume();
	}
	else
	{
		if (!pWCDI->IsComputationInProgress())
			return;

		//pause computation
		pWCDI->Pause();
		GetDlgItem(IDC_BUTTON_PAUSE)->EnableWindow(FALSE);
		GetDlgItem(IDC_BUTTON_GO)->EnableWindow(FALSE);
		GetDlgItem(IDC_BUTTON_PAUSE)->SetWindowText("Pausing..");
	}
}

void CWCDriver3DDlg::PausedComputation()
{
	GetDlgItem(IDC_BUTTON_PAUSE)->EnableWindow(TRUE);
	GetDlgItem(IDC_BUTTON_PAUSE)->SetWindowText("Resume");
}

void CWCDriver3DDlg::SleepX(int xMilliseconds)
{
	Sleep(xMilliseconds);
}

void CWCDriver3DDlg::OnCancel()
{
#ifndef COMPILE_AND // AND: no questions whether to close or not
	if(pWCDI->IsComputationInProgress())
	{
		if(AfxMessageBox("The computation is in progress.\nExit program?",MB_YESNO | MB_ICONQUESTION) == IDNO)
			return;
	}
	else if (model_modified)
	{
		if(AfxMessageBox("The model has not been saved. All changes will be lost.\nExit program?",MB_YESNO | MB_ICONQUESTION) == IDNO)
			return;
	}
	else
	{
		if(!bClosingDialog)
			if(AfxMessageBox("Exit program?",MB_YESNO | MB_ICONQUESTION) == IDNO)
				return;
	}
#endif // COMPILE_AND
	CDialog::OnCancel();
}

void CWCDriver3DDlg::OnClose() 
{
	bClosingDialog = true;

	CDialog::OnClose();
}

const int min_options_index = 100; //old: 1


void CWCDriver3DDlg::ReadConfigFile()
{
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ElementDataContainer edc;

	//read EDC from config file:
	int read_success =	pWCDI->CallCompFunction(106,0,0,&edc);

	if(read_success)
	{
		if (edc.TreeFind("HOTINTConfiguration.HOTINTversion") != 0)
		{
			HotintVersionInfo config_version(edc.TreeGetDouble("HOTINTConfiguration.HOTINTversion",0.));

			if (config_version == pWCDI->GetHotintVersion())
			{		
				//write configuration parameters in WCDriver, GLDrawWnd
				EDC2Configuration(edc);
				return;
			}
			else
			{
				char str[512];
				sprintf_s(str, "The HOTINT Version %s of the configuration file\nis different from the actual HOTINT Version %s\nand is therefore ignored!", config_version.GetString(), pWCDI->GetHotintVersion().GetString());
#ifndef COMPILE_AND
				AfxMessageBox(str);
#endif	
			}
		}
	}
	
	// do this, when there is no config-file or wrong version number
	OnResultsViewoutput(); //$ DR 2013-02-26 added this line, to show output window
	OnResultsDatamanager(); //$ DR 2013-02-27 added this line, to show data manager
	return;

}

void CWCDriver3DDlg::Configuration2EDC(ElementDataContainer& edc)
{
	edc.TreeSetDoubleC("HOTINTConfiguration.HOTINTversion", pWCDI->GetHotintVersion().GetDoubleValue(), "Version with which this file has been generated");
	// first we write the positions of windows
	CRect rect;
	GetWindowRect(&rect);

	edc.TreeSetIntC("HOTINTConfiguration.Application_Window.rect_left",rect.left,"left coordinate of application window");
	edc.TreeSetIntC("HOTINTConfiguration.Application_Window.rect_top",rect.top,"top coordinate of application window");
	edc.TreeSetIntC("HOTINTConfiguration.Application_Window.rect_width",rect.Width(),"width of application window");
	edc.TreeSetIntC("HOTINTConfiguration.Application_Window.rect_height",rect.Height(),"left coordinate of application window");

	// now the parameters of the OpenGL scene and of the frames recording
	GLDrawWnd.Configuration2EDC(edc);
 

	// and finally the state of the control windows
	BOOL DataManagerWindowOpen = 0;
	if(::IsWindow(DialogDataManager.m_hWnd))
		DataManagerWindowOpen = 1;
	edc.TreeSetBoolC("HOTINTConfiguration.data_manager_open",DataManagerWindowOpen,"open data manager on startup");

	//ar << DialogDataManager.m_strDataRequest; //JG: not used anymore

	BOOL OutputDialogOpen = 0;
	int storedwidth = 200; //verify with outputdialog.cpp constructor!!!
	if(::IsWindow(dialogComputationOutput.m_hWnd))
	{
		dialogComputationOutput.GetWindowRect(&rect);
		storedwidth = rect.Width();
		OutputDialogOpen = 1;
	}
	edc.TreeSetBoolC("HOTINTConfiguration.Output_dialog.dialog_open",OutputDialogOpen,"open output dialog on startup");
	edc.TreeSetIntC("HOTINTConfiguration.Output_dialog.stored_width",storedwidth,"stored width of output dialog");
	edc.TreeSetBoolC("HOTINTConfiguration.Output_dialog.enable_output_text",enableOutputText,"enable output text in output dialog");

}

void CWCDriver3DDlg::EDC2Configuration(const ElementDataContainer& edc)
{
	int left, top, width, height;
	CRect rect;
	GetWindowRect(&rect);

	int left_default = 795;
	int top_default = 595;
	int width_default = 700;	//4000
	int height_default = 700;	//4000

	left = edc.TreeGetInt("HOTINTConfiguration.Application_Window.rect_left",left_default);
	top = edc.TreeGetInt("HOTINTConfiguration.Application_Window.rect_top",top_default);
	width = edc.TreeGetInt("HOTINTConfiguration.Application_Window.rect_width",width_default);
	height = edc.TreeGetInt("HOTINTConfiguration.Application_Window.rect_height",height_default);
	
	//$ DR 2013-03-27 removed setting to default screen
	////change coordinates, if outside of standard screen 800x600!!!!
	//if (left > left_default) left = left_default;
	//if (top > top_default) top = top_default;
	//if (width > width_default) width = width_default;
	//if (height > height_default) height = height_default;

	MoveWindow(left,top,width,height, TRUE);

	GLDrawWnd.EDC2Configuration(edc);

	BOOL DataManagerWindowOpen;
	DataManagerWindowOpen = edc.TreeGetBool("HOTINTConfiguration.data_manager_open",1); //$ DR 2013-03-18 changed default to 1
	DialogDataManager.m_strDataRequest = ""; //this will be erased in the future
	// and finally the state of the control windows
	if(DataManagerWindowOpen)
		OnResultsDatamanager();


	BOOL OutputDialogOpen;
	OutputDialogOpen = edc.TreeGetBool("HOTINTConfiguration.Output_dialog.dialog_open",1);	//$ DR 2013-02-26 changed default to 1
	dialogComputationOutput.storedwidth = edc.TreeGetInt("HOTINTConfiguration.Output_dialog.stored_width");
	if(OutputDialogOpen)
		OnResultsViewoutput();

	enableOutputText = edc.TreeGetBool("HOTINTConfiguration.Output_dialog.enable_output_text");
	enableOutputText = (enableOutputText+1)%2; //toggle 1.
	OnResultsEnableoutput();									 //toggle 2. ==> original
}



void CWCDriver3DDlg::OnOptionsSaveconfiguration()
{
	GLDrawWnd.SetProgramDirectoryAsCurrent();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ElementDataContainer edc;
	//generate EDC data and pass to MBS
	Configuration2EDC(edc);
	pWCDI->CallCompFunction(107,0,0,&edc);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	////old code, will be erased:
	//CFile file;
	//if(!file.Open(CONFIG_FILE_NAME,CFile::modeWrite | CFile::modeCreate))
	//{ 
	//	AfxMessageBox("Could not open the config file.", MB_ICONSTOP);
	//	return;
	//}

	//CArchive ar(&file,CArchive::store);

	//ar << HOTINT_version_info; //write the actual version of hotint!!!

	//// first we write the positions of windows
	//CRect rect;
	//GetWindowRect(&rect);
	//ar << (int)rect.left << (int)rect.top;
	//ar << (int)rect.Width() << (int)rect.Height();
	////GetDlgItem(IDC_EDIT_CONSOLE)->GetWindowRect(&rect);
	////ar << rect.Width();

	//// now the parameters of the OpenGL scene and of the frames recording
	//GLDrawWnd.SaveConfig(ar);

	////save options:
	////computation:
	//for (int i=001; i <= 300; i++)
	//{
	//	int val=0;
	//	if (i < min_options_index) ar << val;
	//	else ar << pWCDI->GetIOption(i); 
	//}
	//for (int i=001; i <= 300; i++)
	//{
	//	if (i!=2) //exclude start time
	//	{
	//		double val=0.;
	//		if (i < min_options_index) ar << val;
	//		else ar << pWCDI->GetDOption(i); 
	//	}
	//}
	////Text options
	//for (int i=1; i <= 300; i++)
	//{
	//	if (i < min_options_index || pWCDI->GetTOption(i) == 0)
	//		ar << CString(""); 
	//	else
	//		ar << CString(pWCDI->GetTOption(i));
	//} 
 //

	//// and finally the state of the control windows
	//BOOL DataManagerWindowOpen = 0;
	//if(::IsWindow(DialogDataManager.m_hWnd))
	//	DataManagerWindowOpen = 1;
	//ar << DataManagerWindowOpen;
	//ar << DialogDataManager.m_strDataRequest;

	//BOOL OutputDialogOpen = 0;
	//int storedwidth = 200; //verify with outputdialog.cpp constructor!!!
	//if(::IsWindow(dialogComputationOutput.m_hWnd))
	//{
	//	dialogComputationOutput.GetWindowRect(&rect);
	//	storedwidth = rect.Width();
	//	OutputDialogOpen = 1;
	//}
	//ar << storedwidth;
	//ar << OutputDialogOpen;

	//ar << enableOutputText;

}

void CWCDriver3DDlg::OnEditEditsolveroptions()
{
	CTreeViewCustomEditDialog ced; 
	ced.SetWCDI(pWCDI); 
	ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc_solveroptions;
	pWCDI->CallCompFunction(210,0,1,&edc_solveroptions); // call with value == 1 for solver options only

	ced.SetDialogName("Edit Solver Options");
	ced.SetElementDataContainer(&edc_solveroptions);
	ced.SetUseApply(1, 210, 1, 1);		// these are the codes of the functions for various actions available in the dialog
	ced.DoModal();
}

void CWCDriver3DDlg::OnOptionsSaveSolverOptions()
{
	GLDrawWnd.SetProgramDirectoryAsCurrent();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ElementDataContainer edc;
	//generate EDC data and pass to MBS
	Configuration2EDC(edc);
	pWCDI->CallCompFunction(107,0,1,&edc); // call with value == 1
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void CWCDriver3DDlg::OnEditEdithotintoptions()
{
	CTreeViewCustomEditDialog ced; 
	ced.SetWCDI(pWCDI); 
	ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc_hotintoptions;
	pWCDI->CallCompFunction(210,0,2,&edc_hotintoptions); // call with value == 2 for hotint options only

	ced.SetDialogName("Edit Hotint Options");
	ced.SetElementDataContainer(&edc_hotintoptions);
	ced.SetUseApply(1, 210, 1, 2);		// these are the codes of the functions for various actions available in the dialog
	ced.DoModal();
}

void CWCDriver3DDlg::OnOptionsSaveHotintOptions()
{
	GLDrawWnd.SetProgramDirectoryAsCurrent();

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	ElementDataContainer edc;
	//generate EDC data and pass to MBS
	Configuration2EDC(edc);
	pWCDI->CallCompFunction(107,0,2,&edc); // call with value == 2
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void CWCDriver3DDlg::OnFileExit()
{
	OnCancel();
}

void CWCDriver3DDlg::OnHelpAbout()
{
	CAboutDlg dlgAbout;
	dlgAbout.SetpWCDI(pWCDI);
	dlgAbout.DoModal();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CWCDriver3DDlg::CallWCDriverFunction(int action, int option, int value, ElementDataContainer* edc)
{
	//action=1: redraw
	//action=2: Edit loads of element 'value', load number 'option'

	int rv = 1;
	//insert actions performed by WCDriver on call by MBS-System

	switch(action)
	{
	case 1: //Redraw
		{
			GLDrawWnd.SendMessage(WM_REDRAW);
			Dialog2DView.ForceUpdate();
			//AD: Added update options and redraw all PlotTool Windows
			for (int i=1; i<= PlotToolArray.Length(); i++)
			{
				if(option==0) PlotToolArray(i)->Redraw();
				if(option==1) PlotToolArray(i)->UpdateDialogData();
			}
			break;
		}
	case 2: //Edit Dialog for load
		{
			rv = 0; //ok, nothing to delete

			ElementDataContainer edc;
			int nelem = value;
			int loadnum = option;
			pWCDI->GetElementData(nelem, 2, loadnum, edc); //type==2 --> element load
			CustomEditDialog ced; ced.SetWCDI(pWCDI); ced.SetGLDrawWnd(&GLDrawWnd);

			ced.SetDialogName("Edit load");

			ced.SetDeleteButton(1);

			ced.SetElementDataContainer(&edc);
			if (ced.DoModal() != IDCANCEL)
			{
				//write data to element;
				pWCDI->SetElementData(nelem, 2, loadnum, edc); //type2==load
			} 
			else
			{
				if (ced.GetDeleteFlag())
				{
					pWCDI->CallCompFunction(5, loadnum, nelem);
					//inform calling function that force is deleted (not needed)
					rv = ced.GetDeleteFlag();
				}
				else
				{
					rv = 0;
				}
			}

			break;
		}
	case 3: //Open File Requester for STL
		{
			static char BASED_CODE szFilter[] = "STL file (*.stl)|*.stl|All Files (*.*)|*.*||";

			CFileDialog fd(TRUE,"mbs",0,OFN_OVERWRITEPROMPT|OFN_ENABLESIZING|OFN_FILEMUSTEXIST ,szFilter); //true=load, false=save
			if(fd.DoModal() == IDCANCEL) 
			{
				rv = 0;
			}
			else
			{
				edc->Reset();
				ElementData ed;
				ed.SetText(fd.GetPathName(), "PathName");
				edc->Add(ed);
			}
			break;
		}
	case 4: //Set Program directory as current
		{
			GLDrawWnd.SetProgramDirectoryAsCurrent();
			break;
		}
	case 5: // Quit Program from MBS
		{
			CDialog::OnCancel();
			break;
		}
	case 6: //reset model
		{
			ModelChanged();
			break;
		}
	case 7: // close Eigenvalue computation window
		{
			if (::IsWindow(computeEigenmodes.m_hWnd)) computeEigenmodes.SendMessage(WM_CLOSE_EV_WINDOW);
			break;
		}
	case 8: // open sensorwatch for sensor "value" (PlotTool)
		{
			if (value>0)
			{
				CPlotToolDlg* activeplottool;
				int nr = GetPlotToolDlg(activeplottool,-1); // sensor always opens a new plottool ?
				CRect viewrect = ComputeDisplayRectForSensorWatch(nr);
				activeplottool->InitializeWithSensorTY(value, viewrect);
			}
			break;
		}
	case 9: // close all PlotTool instances
		{
			for (int i=1; i<= PlotToolArray.Length(); i++)
				PlotToolArray(i)->ExternClose();
			PlotToolArray.Flush();
			break;
		}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Progress Bar
//action=10: activate ProgressBar, options := (0)hide, (1)show, (2)set total ticks, (3)set actual ticks, (4)change status text, (5)pos/size
	case 10:
		{
			if (!::IsWindow(dialogProgressbar.m_hWnd) && !(option==1)) return 0; // no window and not "show" operation

			switch(option)
			{
			case 0: 	
//				if(::IsWindow(dialogProgressbar.m_hWnd))
					dialogProgressbar.ShowWindow(SW_HIDE);
				break;

			case 1: 
				if(!::IsWindow(dialogProgressbar.m_hWnd))
				{
					dialogProgressbar.Create(IDD_DIALOG_PROGRESSBAR,this);
					dialogProgressbar.Reset();
				}
				else
					dialogProgressbar.ShowWindow(SW_SHOW);
				break;

			case 2: 
//				if(::IsWindow(dialogProgressbar.m_hWnd))
					dialogProgressbar.SetTotalTicks(value);
				break;
			case 3:
//				if(::IsWindow(dialogProgressbar.m_hWnd))
					dialogProgressbar.SetActualTicks(value);
				break;
			case 4:
	//			if(::IsWindow(dialogProgressbar.m_hWnd))
				{
					mystr caption = edc->TreeGetString("Caption");
					mystr over = edc->TreeGetString("Over");
					mystr under = edc->TreeGetString("Under");
				
					if (caption.Length()) dialogProgressbar.SetCaptionText(caption);
					if (over.Length()) dialogProgressbar.SetTextOver(over);
					if (under.Length()) dialogProgressbar.SetTextUnder(under);
				}
				break;
//			case 5: break;
			default: assert(0); break;
			}
			break;
		}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Copy data (loaded from external file) to DataManager (e.g. Eigenmodes,...)
//action = 20, option:= number of vectors to store (e.g. number of eigenmodes), edc:= contains the vectors (entries "SV#" & "DV#")
	case 20:
		{
// (AD) under construction
			if(!::IsWindow(DialogDataManager.m_hWnd))
				DialogDataManager.Create(this);
		  
			DialogDataManager.RemoveAll(); // clear DialogDataManager


      double* solvec = NULL;
			double* datavec = NULL;
			int n_sol, n_data;
			
			for(int j=1; j<=option; j++)
			{
// entries for i'th solution vector
				edc->TreeGetVector(mystr("SV")+mystr(j), &solvec, n_sol);
// entries for i'th data vector
				edc->TreeGetVector(mystr("DV")+mystr(j), &datavec, n_data);

				DataStorage storage;
				WCDInterface::DataSaver& ds = storage; // must be cast to DataSaver to write data !!!
				ds.SetTime(j);
// solution vector
				ds << n_sol; 
				for (int i=0; i<n_sol; i++)
				{
					ds << solvec[i];
				}
// data vector
				ds << n_data; 
				for (int i=0; i<n_data; i++)
				{
					ds << datavec[i];
				}

				DialogDataManager.AddEntry(storage, j);
			}
			break;
		}
	default: ;
	}
	return rv;
}

void CWCDriver3DDlg::OnBnClickedButtonSaveMbs()
{
	OnFileSavembs();
}

void CWCDriver3DDlg::OnFileSavembs()
{
	mystr filename = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.hotint_input_data_filename");

	if (filename.Compare("")) 
	{
		OnFileSavembsas();
	}
	else
	{
		int pos = filename.Find(".hmc");
		int l = filename.Length();
		if(pos == l-4)
		{
			ElementDataContainer edc;
			ElementData ed;
			ed.SetText(filename, "File_name"); edc.Add(ed);
			ed.SetText("", "Directory_name"); edc.Add(ed);

			pWCDI->CallCompFunction(101, 1, 0, &edc); //no return value, should always work!

			SetModelModified(0);
			AddRecentFile(filename);
		}
		else
		{
			OnFileSavembsas();
		}
	}
}

void CWCDriver3DDlg::OnFileSavembsas()
{
	mystr path = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");
	static char BASED_CODE szFilter[] = "HMC file (*.hmc)|*.hmc|";
	CFileDialog fd(FALSE,"mbs",0,OFN_HIDEREADONLY|OFN_ENABLESIZING|OFN_OVERWRITEPROMPT,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) 
	{
		SetCurrentDirectory(path);	//$ DR 2012-12-06 reset path to old path
		return;
	}
	SetCurrentDirectory(path);	//$ DR 2012-12-06 reset path to old path

	if(fd.GetFileExt().Compare("hmc"))
	{
		AfxMessageBox("It is just possible to save *.hmc files!");
		return;
	}

	CString pathname = fd.GetPathName();

	if (strlen(fd.GetFileName()) != 0)
	{
	  //CMenu* mmenu = GetMenu();
		//mmenu->ModifyMenu(ID_FILE_RECENTFILE1, MF_BYCOMMAND, ID_FILE_RECENTFILE1, CString("&Recent: ")+CString(pWCDI->GetTOption(94)));
		//mmenu->EnableMenuItem(ID_FILE_RECENTFILE1, MF_BYCOMMAND | MF_ENABLED); //MF_DISABLED | MF_GRAYED

		ElementDataContainer edc;
		ElementData ed;
		ed.SetText(pathname, "File_name"); edc.Add(ed);
		ed.SetText("", "Directory_name"); edc.Add(ed);

		pWCDI->CallCompFunction(101, 1, 0, &edc); //no return value, should always work!

		pWCDI->SetTOption(110, pathname);
		pWCDI->MBS_EDC_TreeSetString(pathname,"GeneralOptions.ModelFile.hotint_input_data_filename");

		CString filename = fd.GetFileName();
		pWCDI->SetTOption(112,filename);
		pWCDI->MBS_EDC_TreeSetString(filename,"GeneralOptions.Paths.internal_model_function_name");

		SetModelModified(0);
		AddRecentFile(pathname);
	}
}

void CWCDriver3DDlg::AutoSavembs()
{
	ElementDataContainer edc;
	ElementData ed;
	mystr path = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");

	int activate = 0;
	pWCDI->MBS_EDC_TreeGetInt(activate,"GeneralOptions.Application.activate_autosave");
	if(activate)
	{
		if(asv_models.Length())
		{
			if(asv_models.Right(1).MakeInt() == MAX_AUTOSAVE_MODELS)
			{
				asv_models += mystr(1);
			}
			else
			{
				asv_models += mystr(asv_models.Right(1).MakeInt()+1);
			}

			if(asv_models.Length() > MAX_AUTOSAVE_MODELS)
			{
				asv_models.EraseChar(1);
			}
		}
		else
		{
			asv_models = "1";
		}
		mystr asv_model_name = "model_asv_" + asv_models.Right(1) + ".hmc";
		//ed.SetText("model_asv.hmc", "File_name");  edc.Add(ed);
		ed.SetText(asv_model_name, "File_name");  edc.Add(ed);
		ed.SetText(path, "Directory_name"); edc.Add(ed);

		pWCDI->CallCompFunction(101, 0, 0, &edc); //no return value, should always work!
	}
}

void CWCDriver3DDlg::OnFileOpenmbs()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, file can not be loaded!");
		return;
	}

	if(pWCDI->GetNElements() != 0 && !(AfxMessageBox("Discard all body information?",MB_YESNO | MB_ICONQUESTION) == IDYES))
	{
		return;
	}

	mystr path = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");
	//static char BASED_CODE szFilter[] = "HMC and HID files (*.hmc)|*.hmc|HID file (*.hid)|*.hid|All Files (*.*)|*.*||";
	static char BASED_CODE szFilter[] = "HMC and HID files (*.hmc;*.hid)|*.hmc; *.hid|All Files (*.*)|*.*||";
	CFileDialog fd(TRUE,"mbs",0,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ENABLESIZING|OFN_FILEMUSTEXIST ,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) 
	{
		SetCurrentDirectory(path);	//$ DR 2012-12-06 reset path to old path
		return;
	}
	SetCurrentDirectory(path);	//$ DR 2012-12-06 reset path to old path

	GLDrawWnd.prohibit_redraw = 1;

	CString pathname = fd.GetPathName();
	CString filename = fd.GetFileName();

	ElementDataContainer edc; ElementData ed;
	ed.SetText(pathname, "File_name"); edc.Add(ed);
	ed.SetText("", "Directory_name"); edc.Add(ed);


	pWCDI->SetTOption(110, pathname);
	pWCDI->MBS_EDC_TreeSetString(pathname,"GeneralOptions.ModelFile.hotint_input_data_filename");
	pWCDI->SetModelData_Initialized(0);
	pWCDI->CallCompFunction(1, 1, 0, &edc);				// Initialize MBS

	AddRecentFile(pathname);
	pWCDI->SetTOption(112,filename);
	pWCDI->MBS_EDC_TreeSetString(filename,"GeneralOptions.Paths.internal_model_function_name");

	pWCDI->CallCompFunction(210, 0, 0, &edc);		//$ DR 2012-12-06 update EDC of MBS
	ModelChanged();
}



void CWCDriver3DDlg::OnEditUndo()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, file can not be loaded!");
		return;
	}

	mystr pathname = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");

	if(asv_models.Length())
	{
		pathname += "model_asv_";
		pathname +=	asv_models.Right(1);
		pathname += ".hmc";
		asv_models.SetLength(asv_models.Length()-1);

		ElementDataContainer edc; ElementData ed;
		ed.SetText(pathname, "File_name"); edc.Add(ed);
		ed.SetText("", "Directory_name"); edc.Add(ed);

		pWCDI->SetTOption(110, pathname);
		pWCDI->MBS_EDC_TreeSetString(pathname,"GeneralOptions.ModelFile.hotint_input_data_filename");
		pWCDI->SetModelData_Initialized(0);
		pWCDI->CallCompFunction(1, 1, 0, &edc);				// Initialize MBS

		ModelChanged();
	}
	else
	{
		AfxMessageBox("No autosaved model available.");
		return;
	}
}


void CWCDriver3DDlg::OnFileRecentfile1()
{
	//return; //$ DR 2012-12-07 not implemented yet

	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, MBS can not be loaded!");
		return;
	}

	mystr filename = 	pWCDI->GetTOption(113);
	//mystr filename = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.recent_file1");

	if (filename.Find(":")>0)		//it is a file to load
	{
		ElementDataContainer edc; ElementData ed;
		ed.SetText(filename, "File_name"); edc.Add(ed);
		ed.SetText("", "Directory_name"); edc.Add(ed);

		pWCDI->SetTOption(110, filename);
		pWCDI->MBS_EDC_TreeSetString(filename,"GeneralOptions.ModelFile.hotint_input_data_filename");
		pWCDI->SetModelData_Initialized(0);
		pWCDI->CallCompFunction(1, 1, 0, &edc);				// Initialize MBS

		pWCDI->CallCompFunction(210, 0, 0, &edc);		// update EDC of MBS
		ModelChanged();
	}
	else	// internal model
	{
		ElementDataContainer edc;
		pWCDI->CallCompFunction(20,1,0,&edc); //edc contains model function names

		int found = 0; 
		mystr model = "";
		//found = edc.Find(filename); // not working

		for(int i=1; i<=edc.Length(); i++)
		{
			mystr model = edc.Get(i).GetText();
			if(model.Compare(filename))
			{
				found = i;
				i = edc.Length()+1; // to break for loop
			}
		}

		if (found)
		{
			pWCDI->MBS_EDC_TreeSetString(filename, "GeneralOptions.ModelFile.internal_model_function_name");
			pWCDI->MBS_EDC_TreeSetString("","GeneralOptions.ModelFile.hotint_input_data_filename");

			pWCDI->CallCompFunction(20,2); //tell MBS that model function has been changed!
			ModelChanged();
		}
		else
		{
			AfxMessageBox("The specified internal model is not available anymore!");
			return;
		}
	}
}



void CWCDriver3DDlg::OnBnClickedButtonRecentFile()
{
	OnFileRecentfile1();
}


void CWCDriver3DDlg::OnComputationLoadinitialvector()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running!\nInitial vector can not be loaded.");
		return;
	}

	static char BASED_CODE szFilter[] = "MBS file (*.sol)|*.sol|All Files (*.*)|*.*||";

	CFileDialog fd(TRUE,"mbs",0,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ENABLESIZING|OFN_FILEMUSTEXIST ,szFilter); //true=load, false=save

	if(fd.DoModal() == IDCANCEL) {return;}

	CString pathname = fd.GetPathName();
	int plen = fd.GetPathName().GetLength(); //dir+filename
	int flen = fd.GetFileName().GetLength(); //filename
	CString dirname = pathname.Left(plen-flen);

	ElementDataContainer edc; ElementData ed;
	ed.SetText(fd.GetFileName(), "File_name"); edc.Add(ed);
	ed.SetText(dirname, "Directory_name"); edc.Add(ed);

	int modified = GetModelModified();
	ModelChanged();
	SetModelModified(modified);

	if (!pWCDI->CallCompFunction(104, 0, 0, &edc))
	{
		AfxMessageBox("Error: Initial vector could not be loaded!");
		return;
	}
	pWCDI->CallCompFunction(1); //Set initial conditions

	GLDrawWnd.ContentsChanged(); //fit on screen and redraw
}

void CWCDriver3DDlg::OnComputationStoresolutionvector()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running!\nSolution can not be stored.");
		return;
	}

	static char BASED_CODE szFilter[] = "MBS file (*.sol)|*.sol|All Files (*.*)|*.*||";

	CFileDialog fd(FALSE,"mbs",0,OFN_HIDEREADONLY|OFN_ENABLESIZING|OFN_OVERWRITEPROMPT,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) {return;}


	CString pathname = fd.GetPathName();
	int plen = fd.GetPathName().GetLength(); //dir+filename
	int flen = fd.GetFileName().GetLength(); //filename
	CString dirname = pathname.Left(plen-flen);

	ElementDataContainer edc; ElementData ed;
	ed.SetText(fd.GetFileName(), "File_name"); edc.Add(ed);
	ed.SetText(dirname, "Directory_name"); edc.Add(ed);

	if (!pWCDI->CallCompFunction(103, 0, 0, &edc))
	{
		AfxMessageBox("Error: Initial vector could not be saved!");
		return;
	}

}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Computation:

void CWCDriver3DDlg::OnResultsDatamanager()
{
	if(!::IsWindow(DialogDataManager.m_hWnd))
	{
		DialogDataManager.Create(this);
	}
	else
	{
		DialogDataManager.CheckScreenLocation(); //check that data manager is not outside screen
	}
}

void CWCDriver3DDlg::OnViewShowContolWindow()
{
	if(!::IsWindow(Dialog2DView.m_hWnd))
	{
		Dialog2DView.Create(IDD_DIALOG_2D_VIEW,this);
		Dialog2DView.ShowWindow(SW_SHOW);
	}
	else
	{
		Dialog2DView.ShowWindow(SW_SHOW);
		Dialog2DView.CheckScreenLocation();
	}
}


void CWCDriver3DDlg::OnResultsViewoutput()
{
	if(!::IsWindow(dialogComputationOutput.m_hWnd))
	{
		dialogComputationOutput.Create(this);
		dialogComputationOutput.UpdateText(m_ConsoleText);
	}
}


void CWCDriver3DDlg::OnComputationReset()
{
	if(!pWCDI->IsComputationInProgress())
	{
		pWCDI->CallCompFunction(1); //Set initial conditions //$ DR 2013-10-11 not in ModelChanged anymore
		ModelChanged();
	}
	this->AddText("Computation reset\n");
}

void CWCDriver3DDlg::OnComputationPrinttimings()
{
	if(!pWCDI->IsComputationInProgress())
	{
		pWCDI->PrintTimingList();
	}
}

void CWCDriver3DDlg::OnResultsShowoutputstatictext()
{
	CDialogReadText drt;
	drt.m_Text = m_ConsoleText;
	drt.DoModal();
}

void CWCDriver3DDlg::OnResultsEnableoutput()
{
	enableOutputText = (enableOutputText+1)%2;

	if (enableOutputText) GetMenu()->CheckMenuItem(ID_RESULTS_ENABLEOUTPUT, MF_CHECKED);
	else GetMenu()->CheckMenuItem(ID_RESULTS_ENABLEOUTPUT, MF_UNCHECKED);

}


//void CWCDriver3DDlg::OnComputationAssemblesystem()
//{
//	if(!pWCDI->IsComputationInProgress())
//	{
//		pWCDI->CallCompFunction(3);
//		GLDrawWnd.ContentsChanged(); //fit on screen and redraw	}
//	}
//}

void CWCDriver3DDlg::OnSystemShowsystemproperties()
{
	pWCDI->CallCompFunction(301);
}


void CWCDriver3DDlg::OnSystemChecksystem()
{
	if (pWCDI->CallCompFunction(302))
	{
		AfxMessageBox("The system shows inconsistencies! See the output window and\ncheck the element, constraint, sensor and GeomElement numbers!");
	}
}


//void CWCDriver3DDlg::OnComputationAssigninitialvector()
//{
//	if(!pWCDI->IsComputationInProgress())
//	{
//		pWCDI->CallCompFunction(1);
//		GLDrawWnd.ContentsChanged(); //fit on screen and redraw
//	}
//}

void CWCDriver3DDlg::OnFileNewmbs()
{
	if(!pWCDI->IsComputationInProgress())
	{
		if(AfxMessageBox("Discard all body information?",MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->SetTOption(112, "new model");
			pWCDI->MBS_EDC_TreeSetString("new model", "GeneralOptions.ModelFile.internal_model_function_name");

			pWCDI->SetTOption(110, "");
			pWCDI->MBS_EDC_TreeSetString("", "GeneralOptions.ModelFile.hotint_input_data_filename");
			
			SetModelModified();
			pWCDI->CallCompFunction(20,2); //tell MBS that model function has been changed!

			pWCDI->SetTOption(112, "new model");
			pWCDI->MBS_EDC_TreeSetString("new model", "GeneralOptions.ModelFile.internal_model_function_name"); // not working yet

			ModelChanged();
			this->AddText("MBS System reset");
		}
	}
	else
	{
		AfxMessageBox("Computation is running, system can not be reset!");
	}
}

void CWCDriver3DDlg::OnResultsPlotsensor()
{
	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Select sensor", "Select sensor to plot:");
	list.DoModal();

	if (list.item_selected)
	{
		//$ RL 2011-6-9: [ define is replaced by hotint-option now.
		int use_plottool;
		pWCDI->MBS_EDC_TreeGetInt(use_plottool, "PlotToolOptions.activate");
		if(use_plottool)
		{			
			CPlotToolDlg* activeplottool;
			int nr = GetPlotToolDlg(activeplottool,-1); // sensor always opens a new plottool ?
			CRect viewrect = ComputeDisplayRectForSensorWatch(nr);
			activeplottool->InitializeWithSensorTY(list.item_selected, viewrect); //CRect(800,150,1100,350)); 
		}
		else
		{
			if (!pWCDI->CallCompFunction(121,0,list.item_selected)) //action=121: plot sensor 'list.item_selected'
			{
				AfxMessageBox("Warning: Sensor does not contain own solution file and therefore can not be plotted!");
			}
		}
		//$ RL 2011-6-9: ]
	}
}

void CWCDriver3DDlg::OnResultsPlot2sensorsxy()
{
	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	ElementDataContainer plotedc;

	list.SetUseEditItem(1);
	list.InitList("XY-Plot Select sensor X", "Select sensor to plot:");
	list.DoModal();

	if (list.item_selected)
	{
		int xsens_nr = list.item_selected;

		list.SetUseEditItem(1);
		list.InitList("XYPlot Select sensor Y", "Select sensor to plot:");
		list.DoModal();

		if (list.item_selected)
		{
			int ysens_nr = list.item_selected;

			int use_plottool;
			pWCDI->MBS_EDC_TreeGetInt(use_plottool, "PlotToolOptions.activate");
			if(use_plottool)
			{			
				CPlotToolDlg* activeplottool;
				int nr = GetPlotToolDlg(activeplottool,-1); // sensor always opens a new plottool ?
				CRect viewrect = ComputeDisplayRectForSensorWatch(nr);

				activeplottool->InitializeWithSensorXY(xsens_nr, ysens_nr, CRect(800,150,1100,350)); 
			}
			else
			{
				plotedc.TreeSetInt("xsensor",list.item_selected);
				plotedc.TreeSetInt("ysensor",list.item_selected);
				if (!pWCDI->CallCompFunction(122,0,0,&plotedc)) //plot two sensors 'edc.xsensor' vs 'edc.ysensor'  in matlab
				{
					AfxMessageBox("Warning: Sensor does not contain own solution file and therefore can not be plotted!");
				}
			}
		}
	}
}

////old code (used before PlotTool)
//void CWCDriver3DDlg::OnResultsPlotnsensors()
//{
//	CustomListDialog list;
//	list.SetPixelSize(200,400);
//  //list.SetUseEditItem(LBS_MULTIPLESEL);
//	ElementDataContainer edc;
//	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names
//
//	for (int i=1; i <= edc.Length(); i++)
//	{
//		char str[16];
//		sprintf_s(str,"%d-", i);
//		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
//	}
//
//	ElementDataContainer plotedc;
//
//	list.SetUseEditItem(1);
//	list.InitList("Time-Plot Select sensor Y1", "Select sensor to plot:");
//	
//	list.DoModal();
//
//	int count = 1;
//	while(list.item_selected)
//	{
//		plotedc.TreeSetInt("y" + mystr(count)+ "sensor",list.item_selected);
//		count++;
//		list.SetUseEditItem(1);
//		list.InitList("Time-Plot Select sensor Y" + mystr(count), "Select sensor to plot:");
//		list.DoModal();
//	}
//	if(!pWCDI->CallCompFunction(123,count-1,0,&plotedc)) //time plot of N sensors in matlab
//	{
//		AfxMessageBox("Warning: Sensor does not contain own solution file and therefore can not be plotted!");
//	}
//}

CRect CWCDriver3DDlg::ComputeDisplayRectForSensorWatch(int i)
{
	CRect desktop;
	GetDesktopWindow()->GetWindowRect(desktop);
	
	int totalheight = desktop.Height();
	int totalwidth = desktop.Width();

	int viewheight;
	int viewwidth;
	pWCDI->MBS_EDC_TreeGetInt(viewheight,"PlotToolOptions.Watches.initial_size_vertical");
	pWCDI->MBS_EDC_TreeGetInt(viewwidth,"PlotToolOptions.Watches.initial_size_horizontal");


	double factor_x = 1.1;  // hard coded factor for border,...
	double factor_y = 1.15; // hard coded factor for border,...

	int spots_x = (int) floor ((double)totalwidth / (viewwidth*factor_x));
	int spots_y = (int) floor ((double)totalheight / (viewheight*factor_y));

// position of window 1
	int p0x = totalwidth - (int) (viewwidth*(factor_x+1.)*0.5 +0.5);
	int p0y =  (int) (viewheight*(factor_y-1.) +0.5);

	
	int shift_y = (i-1)%spots_y;
	int shift_x = ((i-1)%(spots_x*spots_y)) / spots_y;
	int shift_z = (i-1) / (spots_x*spots_y);

	p0y = p0y + (int) (shift_y*viewheight*factor_y + 0.5);
	p0x = p0x - (int) (shift_x*viewwidth*factor_x + 0.5);

	p0x = p0x + (int) (shift_z*viewwidth*factor_x/5. +0.5);
	p0y = p0y + (int) (shift_y*viewheight*factor_y/5. + 0.5);

// place views such that 
	CRect view(p0x, p0y, p0x+viewwidth, p0y+viewheight);
	return view;
}

void CWCDriver3DDlg::GetSensorWatchDialog(CDialogOneEditControl*& watch, int& dialognum)
{
	for (int i=0; i < MAX_SENSOR_WATCH; i++)
	{
		if (!sensorwatchctrl[i]) 
		{
			sensorwatchctrl[i] = new CDialogOneEditControl();
			watch = sensorwatchctrl[i];
			dialognum = i;
			return;
		} 
		else if (!::IsWindow(sensorwatchctrl[i]->m_hWnd))
		{
			watch = sensorwatchctrl[i];
			dialognum = i;
			return;
		}
	}
	dialognum = 0;
	watch = 0;
	return;
}

void CWCDriver3DDlg::OnResultsSensorwatch()
{
	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[32];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);

	list.InitList("Select sensor", "Select sensor to watch:");

	list.DoModal();

	if (list.item_selected)
	{
		CDialogOneEditControl* cwatch;
		int dialognum;
		GetSensorWatchDialog(cwatch, dialognum);

		if (cwatch)
		{
			char str[32];
			sprintf_s(str, "%d: '", list.item_selected);
			CString name = CString("Sensor ")+CString(str)+CString(edc.Get(list.item_selected).GetDataName())+CString("'");
			cwatch->InitWatch(name, list.item_selected, dialognum);
			cwatch->SetWCDI(pWCDI);
			cwatch->Create(this);		
		}
		else
		{
			AfxMessageBox("Sorry, too many sensor watches opened!");		
		}
	}
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Options:

void CWCDriver3DDlg::OnBnClickedButtonGeneraloptions()
{
	OnViewViewingoptions();

	/*
	if(!::IsWindow(DialogFEDrawingOptions.m_hWnd))
	DialogFEDrawingOptions.Create(this);*/

}

void CWCDriver3DDlg::OnViewRigidbodyJointOptions()
{
	if(!::IsWindow(dialogBodyJointOptions.m_hWnd))
		dialogBodyJointOptions.Create(this);
}

void CWCDriver3DDlg::OnViewFeDrawingoptions()
{
	if(!::IsWindow(DialogFEDrawingOptions.m_hWnd))
		DialogFEDrawingOptions.Create(this);
}

void CWCDriver3DDlg::OnViewViewingoptions()
{
	if(!::IsWindow(dialogViewingOptions.m_hWnd))
		dialogViewingOptions.Create(this);
}

//void CWCDriver3DDlg::OnComputationComputationoptions()		//$ DR 2013-10-16 removed
//{
//	if(!::IsWindow(dialogComputationSettings.m_hWnd))
//		dialogComputationSettings.Create(this);
//}

void CWCDriver3DDlg::OnComputationComputeEigenmodes()
{
	if(!::IsWindow(computeEigenmodes.m_hWnd))
		computeEigenmodes.Create(this);
}

void CWCDriver3DDlg::OnViewOpengloptions()
{
	if(!::IsWindow(dialogOpenGLOptions.m_hWnd))
		dialogOpenGLOptions.Create(this);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CWCDriver3DDlg::ModelChanged()
{
	SetModelModified();
	pWCDI->CallCompFunction(105); //Clear stored initial vector

	//pWCDI->SetDOption2, 0.); 
	pWCDI->MBS_EDC_TreeSetDouble(0.,"SolverOptions.start_time"); //reset start time, otherwise nothing is computed the second time
	DialogDataManager.RemoveAll(); //remove stored output
	//pWCDI->CallCompFunction(3); //Assemble //$ 2013-02-04 not necessary?
	//pWCDI->CallCompFunction(1); //Set initial conditions //$ 2013-02-04 not necessary?
	//int rv = pWCDI->CallCompFunction(302); //Check Consistency //$ DR 2012-12-04 added
	//if(rv < 2)
	//{
	GLDrawWnd.prohibit_redraw = 0;
	GLDrawWnd.ContentsChanged(0); //fit on screen and redraw
	Dialog2DView.ForceUpdate(); //$ DR 2013-03-18 redraw 2D window
	//}
	StatusText(""); //set new status text according to file name, remove finished
	DisplaySelectedModelName();		//$ DR 2012-12-06 added
}


//Add/Edit elements:

//option&1 Edit/Add element
//option&4 Edit/Add load
//option&8 Edit/Add sensor
//option&16 Edit/Add GeomElement
//option&32 Edit/Add Node
//option&64 Edit/Add Material
//option&2 Add instead of edit
//option&256 edit Robot options menu (old)

int CWCDriver3DDlg::EditElementProperties(int nelem, int option)
{
	CTreeViewCustomEditDialog ced;  //$JG:2012-1-22: view element properties in tree view
	//CustomEditDialog ced; 
	ced.SetWCDI(pWCDI); ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc;
	if (option&1)
	{
		pWCDI->GetElementData(nelem, 1, 0, edc); //type=1->element, type=2->load, type=3->Sensor
		if (option&2)
			ced.SetDialogName("Add element");
		else
			ced.SetDialogName("Edit element");
	}
	else if (option&4)
	{
		pWCDI->GetElementData(nelem, 2, 0, edc); //type=1->element, type=2->load, type=3->Sensor
		if (option&2)
			ced.SetDialogName("Add Load");
		else
			ced.SetDialogName("Edit Load");
	}
	else if (option&8)
	{
		pWCDI->GetElementData(nelem, 3, 0, edc); //type=1->element, type=2->load, type=3->Sensor
		if (option&2)
			ced.SetDialogName("Add Sensor");
		else
			ced.SetDialogName("Edit Sensor");
	}
	else if (option&16)
	{
		pWCDI->GetElementData(nelem, 4, 0, edc); //type=1->element, type=2->load, type=3->Sensor, type=4->Sensor
		if (option&2)
			ced.SetDialogName("Add Geometric Element");
		else
			ced.SetDialogName("Edit Geometric Element");
	}
	else if (option&32)
	{
		pWCDI->GetElementData(nelem, 5, 0, edc); //type=1->element, type=2->load, type=3->Sensor, type=4->Sensor, type=5->Node
		if (option&2)
			ced.SetDialogName("Add Node");
		else
			ced.SetDialogName("Edit Node");
	}
	else if (option&64) //material
	{
		pWCDI->GetElementData(nelem, 6, 0, edc); //type=1->element, type=2->load, type=3->Sensor, type=4->Sensor, type=5->Node, type=6->Material
		if (option&2)
			ced.SetDialogName("Add Material");
		else
			ced.SetDialogName("Edit Material");
	}
	else if (option&256)
	{
		pWCDI->GetElementData(nelem, 50, 0, edc); //type=50->Robot options menu (old)
		ced.SetDialogName("Edit Robot options");
	}
	else
	{
		ced.SetDialogName("Edit");
	}

	ElementDataContainer edc_root; //in order to show root entries
	ElementData root;
	const char* rt = "root";
	root.SetEDC(&edc, rt);
	edc_root.Add(root);

	ElementData* edr = edc_root.TreeFind(rt);
	//const ElementDataContainer& edcnew = *edr->GetEDC();
	ElementDataContainer& edcnew = *edr->GetEDC();

	//ced.SetElementDataContainer(&edc);
	ced.SetElementDataContainer(&edc_root);

	int rv = 1;
	if (ced.DoModal() != IDCANCEL)
	{
		if (!(option&2))		// autosave before editing object. For saving before ADDING objects, AutoSave is called directly in these functions
		{
			AutoSavembs();	//$ DR 2013-03-19
		}
		SetModelModified();
		if (option&1)
		{
			//write data to element;
			rv = pWCDI->SetElementData(nelem, 1, 0, edcnew); //type 1==element data
		}
		else if (option&4)	//$ DR 2012-10
		{
			//write data to load;
			rv = pWCDI->SetElementData(nelem, 2, 0, edcnew); //type 2==load data
		} 
		else if (option&8)
		{
			//write data to sensor;
			rv = pWCDI->SetElementData(nelem, 3, 0, edcnew); //type 3==sensor data
		} 
		else if (option&16)
		{
			//write data to geomelement;
			rv = pWCDI->SetElementData(nelem, 4, 0, edcnew); //type 4==geomelement data
		} 
		else if (option&32)
		{
			//write data to node;
			rv = pWCDI->SetElementData(nelem, 5, 0, edcnew); //type 5==node data
		} 
		else if (option&64)
		{
			//write data to material;
			rv = pWCDI->SetElementData(nelem, 6, 0, edcnew); //type 6==material data
		} 
		else if (option&256)
		{
			//write data to geomelement;
			pWCDI->SetElementData(nelem, 50, 0, edcnew); //type 5==node data
		} 
		return rv;
	}
	else
	{
		//delete element:
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Add objects

void CWCDriver3DDlg::OnSpecialAddbody()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;

	int nElements = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCElement);
	for (int i=1; i <= nElements; i++)
	{
		{
			CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCElement,i);
			int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCElement,i);
			
			if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
			{
				if(!(typeflags & TAEconstraint) && !(typeflags & TAEinput_output))
				{
						list.AddString(name, i);
				}
			}
		}
	}

	//needs to be called at last:
	list.InitList("Add body", "Select body type to be added:");
	//add items to list:

	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nelem = pWCDI->GetObjectFactory()->AddObject(OFCElement, list.item_selected);

		if (nelem>0)
		{
			if (EditElementProperties(nelem,3))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(4,0,nelem);
			}
		}
		else
		{
			if (nelem == -1) AfxMessageBox("Error: Selected element type does not exist!!!\n");
		}
	}
}


void CWCDriver3DDlg::OnAddconnectorKinematicpairs()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;

	int nElements = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCElement);
	for (int i=1; i <= nElements; i++)
	{
		{
			CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCElement,i);
			int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCElement,i);
			
			if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
			{
				if(typeflags & TAEconstraint) 
				{
					list.AddString(name, i);
				}
			}
		}
	}

	//needs to be called at last:
	list.InitList("Add Connector", "Select connector type to be added:");
	//add items to list:

	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nelem = pWCDI->GetObjectFactory()->AddObject(OFCElement, list.item_selected);

		if (nelem>0)
		{
			if (EditElementProperties(nelem,3))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(4,0,nelem);
			}
		}
		else
		{
			if (nelem == -1) AfxMessageBox("Error: Selected connector type does not exist!!!\n");
		}
	}
}

void CWCDriver3DDlg::OnAddconnectorControlelement()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;
	int nElements = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCElement);
	for (int i=1; i <= nElements; i++)
	{
		{
			CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCElement,i);
			int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCElement,i);
			
			if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
			{
				if(typeflags & TAEinput_output)
				{
					list.AddString(name, i);
				}
			}
		}
	}

	//needs to be called at last:
	list.InitList("Add Connector", "Select connector type to be added:");
	//add items to list:

	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nelem = pWCDI->GetObjectFactory()->AddObject(OFCElement, list.item_selected);

		if (nelem>0)
		{
			if (EditElementProperties(nelem,3))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(4,0,nelem);
			}
		}
		else
		{
			if (nelem == -1) AfxMessageBox("Error: Selected connector type does not exist!!!\n");
		}
	}
}


void CWCDriver3DDlg::OnAddobjectAddnode()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	// create list for choosing the sensor-type
	CustomListDialog list;
	int nNodes = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCNode);
	for (int i=1; i <= nNodes; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCNode,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCNode,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add node", "Select node to be added:");

	//add items to list:
	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nnode = pWCDI->GetObjectFactory()->AddObject(OFCNode, list.item_selected);

		if (nnode>0)
		{
			if (EditElementProperties(nnode,32+2)) //32==Edit node, 2==Add
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(12,0,nnode);
			}
		}
		else
		{
			if (nnode == -1) AfxMessageBox("Error: Selected node type does not exist!!!\n");
		}
	}
}


void CWCDriver3DDlg::OnSpecialAddsensor()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	// create list for choosing the sensor-type
	CustomListDialog list;
	int nSensors = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCSensor);
	for (int i=1; i <= nSensors; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCSensor,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCSensor,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add sensor", "Select sensor to be added:");

	//add items to list:
	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		//int nsens = pWCDI->CallCompFunction(7,0,list.item_selected); //action=7: add sensor 'list.item_selected
		int nsens = pWCDI->GetObjectFactory()->AddObject(OFCSensor, list.item_selected);

		if (nsens>0)
		{
			if (EditElementProperties(nsens,8+2)) //8==Edit sensor, 2==Add
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(8,0,nsens);
			}
		}
		else
		{
			if (nsens == -1) AfxMessageBox("Error: Selected sensor type does not exist!!!\n");
		}
	}
}


void CWCDriver3DDlg::OnSpecialAddgeomelement()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;
	int nGeomElements = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCGeomElement);
	for (int i=1; i <= nGeomElements; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCGeomElement,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCGeomElement,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add Geometric Element", "Select GeomElement:");

	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		//int nelem = pWCDI->CallCompFunction(9,0,list.item_selected); //action=9: add geomelement 'list.item_selected'
		int nelem = pWCDI->GetObjectFactory()->AddObject(OFCGeomElement, list.item_selected);

		if (nelem>0)
		{
			if (EditElementProperties(nelem,16+2)) //16==Edit drawelement, 2==Add
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(10,0,nelem);
			}
		}
	}
}

//$ DR 2012-12 button added
void CWCDriver3DDlg::OnAddobjectAddmaterial()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;
	int nMaterials = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCMaterial);
	for (int i=1; i <= nMaterials; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCMaterial,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCMaterial,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add Material", "Select Material:");

	list.DoModal();
	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nmat = pWCDI->GetObjectFactory()->AddObject(OFCMaterial,list.item_selected);
		if (nmat>0)
		{
			if (EditElementProperties(nmat,66))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete material:
				pWCDI->CallCompFunction(15,0,nmat);
			}
		}
		else
		{
			if (nmat == -1) AfxMessageBox("Error: Selected material type does not exist!!!\n");
		}
	}
}

//$ DR 2013-01 button added
void CWCDriver3DDlg::OnAddbeam3dproperties()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	CustomListDialog list;
	int nBeamProp = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCBeamProperties);
	//int nMaterials = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCMaterial);
	for (int i=1; i <= nBeamProp; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCBeamProperties,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCBeamProperties,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add BeamProperties", "Select BeamProperties:");

	list.DoModal();
	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		int nmat = pWCDI->GetObjectFactory()->AddObject(OFCBeamProperties, list.item_selected);	

		if (nmat>0)
		{
			if (EditElementProperties(nmat,66))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete material:
				pWCDI->CallCompFunction(15,0,nmat);
			}
		}
		else
		{
			if (nmat == -1) AfxMessageBox("Error: Selected material type does not exist!!!\n");
		}
	}
}

//$ DR 2012-10 button added
void CWCDriver3DDlg::OnAddload()
{
	if(pWCDI->IsComputationInProgress()) {AfxMessageBox("Computation is running! Operation not possible."); return;}

	// create list for choosing the load-type
	CustomListDialog list;
	int nLoads = pWCDI->GetObjectFactory()->GetAvailableTypesCount(OFCLoad);
	for (int i=1; i <= nLoads; i++)
	{
		CString name = pWCDI->GetObjectFactory()->GetTypeName(OFCLoad,i);
		int typeflags = pWCDI->GetObjectFactory()->GetTypeFlags(OFCLoad,i);
		if( !(pWCDI->GetObjectFactory()->ExcludeExperimentalObjects()) || !(typeflags & TAENotInRelease) )
		{
			list.AddString(name, i);
		}
	}

	//needs to be called at last:
	list.InitList("Add load", "Select load type to be added:");

	//add items to list:
	list.DoModal();

	if (list.item_selected)
	{
		AutoSavembs();	//$ DR 2013-03-19
		//int nload = pWCDI->CallCompFunction(6,0,list.item_selected);	// add load of specified type 
		int nload = pWCDI->GetObjectFactory()->AddObject(OFCLoad, list.item_selected);

		if (nload>0)
		{
			if (EditElementProperties(nload,6))
			{
				//pWCDI->SetModelData_Initialized(0);						//$ DR 2013-02-04
				pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
				pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
				ModelChanged();																//$ DR 2013-02-04
			}
			else
			{
				//delete element:
				pWCDI->CallCompFunction(5,0,nload);
			}
		}
		else
		{
			if (nload == -1) AfxMessageBox("Error: Selected load type does not exist!!!\n");
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Edit existing objects

void CWCDriver3DDlg::OnEditEditbody()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}

	int nelem;

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(202,0,0,&edc); //edc contains element names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Element Properties", "Select element to edit:");

	list.DoModal();

	nelem = list.item_selected;

	if (nelem)
	{
		if (EditElementProperties(nelem,1))
		{
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();		
		}
	}

}

void CWCDriver3DDlg::OnEditEditnode()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	int nelem;

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(207,0,0,&edc); //edc contains node names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Node", "Select node to edit:");

	list.DoModal();

	nelem = list.item_selected;

	if (nelem)
	{
		if (EditElementProperties(nelem,32)) //edit node
		{
			ModelChanged();
		}
	}
	else
	{
		//AfxMessageBox("Invalid number selected!");
	}
}

void CWCDriver3DDlg::OnEditEditsensor()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	int nelem;

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Sensor", "Select sensor to edit:");

	list.DoModal();

	nelem = list.item_selected;

	if (nelem)
	{
		if (EditElementProperties(nelem,8)) //edit sensor
		{
			SetModelModified();
			//only redraw to show new sensors, nothing else needs to be done
			GLDrawWnd.Redraw();
		}
	}
	else
	{
		//AfxMessageBox("Invalid number selected!");
	}
}

void CWCDriver3DDlg::OnEditEditgeomelement()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	int nelem;

	CustomListDialog list;
	list.SetPixelSize(250,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(205,0,0,&edc); //edc contains GeomElement names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Geometric Element", "Select GeomElement to edit:");

	list.DoModal();

	nelem = list.item_selected;

	if (nelem)
	{
		if (EditElementProperties(nelem,16)) //edit sensor
		{
			SetModelModified();
			GLDrawWnd.ContentsChanged(0); //fit on screen and redraw
		}
	}
	else
	{
		//AfxMessageBox("Invalid number selected!");
	}
}

void CWCDriver3DDlg::OnEditEditmaterial()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	int nelem;

	CustomListDialog list;
	list.SetPixelSize(250,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(209,0,0,&edc); //edc contains material names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Materials", "Select Material to edit:");

	list.DoModal();

	nelem = list.item_selected;

	if (nelem)
	{
		if (EditElementProperties(nelem,64)) //edit material
		{
			SetModelModified();
			GLDrawWnd.ContentsChanged(0); //fit on screen and redraw
		}
	}
	else
	{
		//AfxMessageBox("Invalid number selected!");
	}
}


//$ DR 2012-10: loads moved from element to edc
void CWCDriver3DDlg::OnEditEditLoad()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}

	int nload;

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(214,0,0,&edc); //edc contains load names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Edit Load Properties", "Select load to edit:");

	list.DoModal();

	nload = list.item_selected;


	if (nload)
	{
		if (EditElementProperties(nload,4))
		{
			ModelChanged();
		}
	}
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Delete existing objects

void CWCDriver3DDlg::OnEditDeleteElement()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(202,0,0,&edc); //edc contains element names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Element", "DELETE element:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the element \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(4,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}


void CWCDriver3DDlg::OnEditDeleteNode()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}

	CustomListDialog list;
	list.SetPixelSize(220,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(207,0,0,&edc); //edc contains node names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Node", "Select node to delete:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the node \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(12,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}




void CWCDriver3DDlg::OnEditDeleteSensor()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	CustomListDialog list;
	list.SetPixelSize(220,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,0,&edc); //edc contains sensor names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Sensor", "Select sensor to delete:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the sensor \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(8,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}


void CWCDriver3DDlg::OnEditDeleteGeomelement()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	CustomListDialog list;
	list.SetPixelSize(250,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(205,0,0,&edc); //edc contains GeomElement names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Geometric Element", "Select GeomElement to delete:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the GeomElement \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(10,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}


void CWCDriver3DDlg::OnEditDeleteMaterial()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}
	CustomListDialog list;
	list.SetPixelSize(250,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(209,0,0,&edc); //edc contains material names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Material", "Select Material to delete:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the material \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(15,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}

void CWCDriver3DDlg::OnEditDeleteLoad()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running! Edit not possible.");
		return;
	}

	CustomListDialog list;
	list.SetPixelSize(200,400);

	ElementDataContainer edc;
	pWCDI->CallCompFunction(214,0,0,&edc); //edc contains load names

	for (int i=1; i <= edc.Length(); i++)
	{
		char str[16];
		sprintf_s(str,"%d-", i);
		list.AddString(CString(str)+CString(edc.Get(i).GetDataName()), i);
	}

	list.SetUseEditItem(1);
	list.InitList("Delete Load", "Select load to delete:");

	list.DoModal();
	int n = list.item_selected;
	if (n != 0)
	{
		if(AfxMessageBox(CString("Are you sure to delete the load \nand reset the computation?"),MB_YESNO | MB_ICONQUESTION) == IDYES)
		{
			pWCDI->CallCompFunction(5,0,n);
			pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
			pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
			ModelChanged();
		}
	}
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//other menu calls:

void CWCDriver3DDlg::OnViewDefaultview1()
{
	GLDrawWnd.OnButtonStandardViewXyz();
}

void CWCDriver3DDlg::OnViewXY()
{
	GLDrawWnd.OnButtonStandardViewXy();
}

void CWCDriver3DDlg::OnViewXZ()
{
	GLDrawWnd.OnButtonStandardViewXz();
}

void CWCDriver3DDlg::OnViewYZ()
{
	GLDrawWnd.OnButtonStandardViewYz();
}



void CWCDriver3DDlg::OnMenuHelp()
{
	//GLDrawWnd.SetProgramDirectoryAsCurrent();
	CString dir = GLDrawWnd.GetProgramDirectory();

	dir.Replace('\\','/');
	dir.Replace(" ","%20");
	
	//CString goal = CString("file:///") + dir + CString("/../documentation/hotint_docu_online.html");
	//AddText(goal);
	//ShellExecute(NULL, "open", "IEXPLORE.EXE", goal, NULL, SW_SHOW);

	//$ DR 2013-01-23 new docu:
	int l = dir.GetLength();
	//AddText(dir);
	dir.MakeReverse();
	//AddText(dir);
	int pos = dir.Find("/");	// position of last /
	dir.MakeReverse();
	CString goal = dir.Left(l-pos) + CString("documentation/");
	AddText(goal);
	//ShellExecute(NULL, "open", "D:/cpp/HotInt_V1/documentation/", NULL, NULL, SW_SHOWNORMAL);
	ShellExecute(NULL, "open", goal, NULL, NULL, SW_SHOWNORMAL);

}




void CWCDriver3DDlg::OnViewRobotoptions()
{
	//AfxMessageBox("robot options");

	// TODO: Fügen Sie hier Ihren Befehlsbehandlungscode ein.
	EditElementProperties(0, 256); //edit robot options

	GLDrawWnd.Redraw();
}


void CWCDriver3DDlg::OnSystemEditmodelparameters()
{

	CTreeViewCustomEditDialog ced; ced.SetWCDI(pWCDI); ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc;
	pWCDI->CallCompFunction(211, 0, 0, &edc); //get model data

	const char* options_text = "SolverOptions";
	//look for SolverOptions:
	int findso = edc.Find(options_text);
	ElementDataContainer* save_so = 0;
	if (findso)
	{
		save_so = edc.Get(findso).GetEDC()->GetCopy();
		edc.Delete(findso);
	}

	ElementDataContainer edc_root; //in order to show root entries
	ElementData root;
	const char* rt = "root";
	root.SetEDC(&edc, rt);
	edc_root.Add(root);

	//ced.SetDialogName("Edit Model Parameters");
	ced.SetDialogName("Show Global Variables");	//$ DR 2013-06-05

	ced.SetElementDataContainer(&edc_root);

	ElementData* edr = edc_root.TreeFind(rt);
	ElementDataContainer* edcnew = edr->GetEDC();

	if (ced.DoModal() != IDCANCEL)
	{
		if (save_so)
		{
			ElementData ed;
			ed.SetEDC(save_so, options_text);
			edcnew->Add(ed);
		}

		pWCDI->CallCompFunction(211, 1, 0, edcnew); //set model data
		//if SolverOptions have been there, replace them:

		pWCDI->GetUserInterface()->CallWCDriverFunction(6); //call CWCDriver3DDlg::ModelChanged()
		pWCDI->CallCompFunction(20,3); //Initialize()
		GLDrawWnd.SendMessage(WM_REDRAW);

	}
}

void CWCDriver3DDlg::OnSystemLoadmodelparameter()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, file can not be loaded!");
		return;
	}

	if(pWCDI->GetNElements() != 0 && !(AfxMessageBox("Discard all body information?",MB_YESNO | MB_ICONQUESTION) == IDYES))
	{
		return;
	}

	static char BASED_CODE szFilter[] = "ModelDataFile (*.txt)|*.txt|All Files (*.*)|*.*||";

	CFileDialog fd(TRUE,"mbs",0,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ENABLESIZING|OFN_FILEMUSTEXIST ,szFilter); //true=load, false=save
	//fd.m_pOFN->lpstrTitle = "Open MBS file";

	if(fd.DoModal() == IDCANCEL) {return;}

	if (fd.GetPathName().GetLength() != 0) //user did select file?
	{
		pWCDI->SetTOption(88, fd.GetPathName()); //store actual whole dirname + filename (=pathname)

		pWCDI->CallCompFunction(211, 2); //load model data

		pWCDI->GetUserInterface()->CallWCDriverFunction(6); //call CWCDriver3DDlg::ModelChanged()
		pWCDI->CallCompFunction(20,3); //Initialize()
		GLDrawWnd.SendMessage(WM_REDRAW);


	}


}

void CWCDriver3DDlg::OnSystemSavemodeldata()
{
	static char BASED_CODE szFilter[] = "ModelDataFile (*.txt)|*.txt|All Files (*.*)|*.*||";

	CFileDialog fd(FALSE,"mbs",0,OFN_HIDEREADONLY|OFN_ENABLESIZING|OFN_OVERWRITEPROMPT,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) {return;}


	if (fd.GetPathName().GetLength() != 0) //user did select file?
	{
		pWCDI->SetTOption(88, fd.GetPathName()); //store actual dirname + filename (=pathname)

		pWCDI->CallCompFunction(211, 3); //save model data
	}

}



//edit all model parameters within an EDC
void CWCDriver3DDlg::OnViewEditalloptions()
{
	CTreeViewCustomEditDialog ced; ced.SetWCDI(pWCDI); ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc;
	pWCDI->CallCompFunction(210, 0, 0, &edc); //get model data

	ced.SetDialogName("Edit All Options");
	ced.SetElementDataContainer(&edc);
	ced.SetUseApply(1, 210, 1, 0);		// these are the codes of the functions for various actions available in the dialog
	ced.DoModal();
}

void CWCDriver3DDlg::OnViewShowprogressbar()
{
	CallWCDriverFunction(10,1);
}

void CWCDriver3DDlg::OnResultsPlottooldialog()
{
	CPlotToolDlg* activeplottool;
	GetPlotToolDlg(activeplottool);
}

// select (autocreate) a PlotToolDialog
// chose number = -1 to force creation of a new dialog
// chose number = 0 to pick first visible
// chose number > 0 to pick a certain plottool dialog ( create a new one if it does not exist
int CWCDriver3DDlg::GetPlotToolDlg(CPlotToolDlg*& dialogPlotTool, int number)
{
	if( number <= PlotToolArray.Length() && number > 0) // pick existing
	{
		dialogPlotTool = PlotToolArray(number);
		dialogPlotTool->ShowWindow(SW_SHOW);
		return number;
	}

	if( number == 0 )
	{
		for (int i=1; i<= PlotToolArray.Length(); i++) // search for first visible
		{	
			if(::IsWindowVisible(PlotToolArray(i)->m_hWnd))
			{
				dialogPlotTool = PlotToolArray(i);
				dialogPlotTool->ShowWindow(SW_SHOW);
				return i;
			}
		}
	}

	// create new dialog
	return CreatePlotToolDlg(dialogPlotTool);
}

// create a new PlotToolDialog
int CWCDriver3DDlg::CreatePlotToolDlg(CPlotToolDlg*& dialogPlotTool)
{
	dialogPlotTool = new CPlotToolDlg(this);

	dialogPlotTool->SetWCDI(pWCDI);
	dialogPlotTool->Create(IDD_DIALOG_PLOTTOOL,this);
	dialogPlotTool->SetWindowTextA(mystr("PlotTool - Window #") + mystr(PlotToolArray.Length()+1));
	dialogPlotTool->ShowWindow(SW_SHOW);

// menu entry ?

	return PlotToolArray.Add(dialogPlotTool);
}

// remove the PlotToolDlg from the list
int CWCDriver3DDlg::DeletePlotToolDlg(CPlotToolDlg* dialogPlotTool)
{
	int nr = PlotToolArray.Find(dialogPlotTool);
	PlotToolArray.Erase(nr);

	for(int i=nr; i<=PlotToolArray.Length(); i++)
	{
		PlotToolArray(i)->SetWindowTextA(mystr("PlotTool - Window #") + mystr(PlotToolArray.Length()+1));  
	}
// menu entry ?

	return PlotToolArray.Length();
}


void CWCDriver3DDlg::OnFileSelectmodel()
{

	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, model can not be changed!");
		return;
	}

	CustomListDialog list;
	list.SetPixelSize(400,300);
	list.dialogname = "Select Model File";
	list.use_custom_check1=0;

	ElementDataContainer edc;
	pWCDI->CallCompFunction(20,1,0,&edc); //edc contains model function names

	for (int i=1; i <= edc.Length(); i++)
	{
		list.AddString(CString(edc.Get(i).GetText()), i);
	}
	if (list.DoModal() == IDOK)
	{
		int nelem = list.item_selected;
		if (CString(pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name")) != CString( edc.Get(nelem).GetText()))
		{
			pWCDI->MBS_EDC_TreeSetString(edc.Get(nelem).GetText(), "GeneralOptions.ModelFile.internal_model_function_name");
			pWCDI->MBS_EDC_TreeSetString("","GeneralOptions.ModelFile.hotint_input_data_filename");

			pWCDI->CallCompFunction(20,2); //tell MBS that model function has been changed!

			AddRecentFile(edc.Get(nelem).GetText());

			ModelChanged();

			//$ DR old code:
			//pWCDI->GetUserInterface()->CallWCDriverFunction(6); //call CWCDriver3DDlg::ModelChanged()
			//pWCDI->CallCompFunction(20,2); //tell MBS that model function has been changed!
			////pWCDI->CallCompFunction(1,1); //Initialize()
			//GLDrawWnd.ContentsChanged();
			////pGLDrawWnd->SendMessage(WM_REDRAW); // does not calculate a new value for Max
	
			//DisplaySelectedModelName();
		}
	}
}

void CWCDriver3DDlg::OnEditComputationSteps()
{
	CTreeViewCustomEditDialog ced; ced.SetWCDI(pWCDI); ced.SetGLDrawWnd(&GLDrawWnd);
	ElementDataContainer edc_computationsteps;

	ced.SetDialogName("NOT IMPLEMENTED YET - Edit Computation Steps - NOT IMPLEMENTED YET");
	ced.SetElementDataContainer(&edc_computationsteps);
	//ced.SetUseApply(1, 210, 1, 2);		// these are the codes of the functions for various actions available in the dialog
	ced.DoModal();


	// TODO: Fügen Sie hier Ihren Befehlsbehandlungscode ein.
}

void CWCDriver3DDlg::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	CDialog::OnKeyDown(nChar, nRepCnt, nFlags);
}

// DR 2013-06-12
void CWCDriver3DDlg::OnSystemRunmacro()
{
	if(pWCDI->IsComputationInProgress())
	{
		AfxMessageBox("Computation is running, macro can not be executed!");
		return;
	}

	mystr path = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");
	//static char BASED_CODE szFilter[] = "HMC and HID files (*.hmc)|*.hmc|HID file (*.hid)|*.hid|All Files (*.*)|*.*||";
	static char BASED_CODE szFilter[] = "HMC and HID files (*.hmc;*.hid)|*.hmc; *.hid|All Files (*.*)|*.*||";
	CFileDialog fd(TRUE,"mbs",0,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ENABLESIZING|OFN_FILEMUSTEXIST ,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) 
	{
		SetCurrentDirectory(path);	// reset path to old path
		return;
	}
	SetCurrentDirectory(path);	// reset path to old path

	GLDrawWnd.prohibit_redraw = 1;

	CString pathname = fd.GetPathName();
	CString filename = fd.GetFileName();

	ElementDataContainer edc; ElementData ed;
	ed.SetText(pathname, "File_name"); edc.Add(ed);
	//ed.SetText("", "Directory_name"); edc.Add(ed);
	ed.SetInt(1,"is_file"); edc.Add(ed);

	AutoSavembs();
	pWCDI->CallCompFunction(113,0,0,&edc);	// Add Model Data
	pWCDI->CallCompFunction(210, 0, 0, &edc);		//$ DR 2012-12-06 update EDC of MBS

	pWCDI->CallCompFunction(3);										// Assemble				//$ DR 2013-10-11 not in ModelChanged anymore 
	pWCDI->CallCompFunction(1);										// SetInitialConditions //$ DR 2013-10-11 not in ModelChanged anymore
	ModelChanged();
}
