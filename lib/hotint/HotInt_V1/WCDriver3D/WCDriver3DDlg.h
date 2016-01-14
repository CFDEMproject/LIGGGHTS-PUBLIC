//#**************************************************************
//# filename:             WCDriver3DDlg.h
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
 
#if !defined(AFX_WCDriver3DDLG_H__270A972D_93F5_46BD_8737_55C336E60107__INCLUDED_)
#define AFX_WCDriver3DDLG_H__270A972D_93F5_46BD_8737_55C336E60107__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CWCDriver3DDlg dialog

#include "mbs_interface.h"
#include "..\WorkingModule\WinCompDriverInterface.h"
#include "gldrawwnd.h"
#include "DialogDataManager.h"
#include "FEdrawingOptions.h"
#include "OutputDialog.h"
#include "DialogViewingOptions.h"
//#include "DialogComputationSettings.h" //$ DR 2013-10-16 removed
#include "DlgBodyJointOpt.h"
#include "OpenGLDlg.h"
#include "DlgOneEditCtrl.h"
#include "ComputeEigenmodes.h"
//#include "DialogComputationSettings.h" //$ DR 2013-10-16 removed
#include "DialogProgressBar.h"
#include "PlotToolDlg.h"
#include "DrawWindow2D.h"
#include "afxwin.h"

const int MAX_SENSOR_WATCH = 20;
const int MAX_AUTOSAVE_MODELS = 3;	//$ DR 2013-05-21
#define WM_CLOSE_EV_WINDOW (WM_USER + 3)

class CWCDriver3DDlg : public CDialog, public WCDInterface::UserInterface, public WCDInterface::ComputationFeedBack
{
	//variables:
	CRect PictureRect;
	int enableOutputText;

	typedef WCDInterface * (*PExpFn)();		// exported-from-dll function type
	HINSTANCE hWorkingDll;
	// pointer to a function exported from the dll
	PExpFn pExpFn;
	WCDInterface * pWCDI;

	CGLDrawWnd GLDrawWnd;

	HCURSOR hCursDragHor, hCursCross;

	void PositionButtons();

	CDialogDataManager DialogDataManager;
	FEDrawingOptions DialogFEDrawingOptions;
	COutputDialog dialogComputationOutput;
	DialogViewingOptions dialogViewingOptions;
	//DialogComputationSettings dialogComputationSettings;	//$ DR 2013-10-16 removed

	DialogBodyJointOptions dialogBodyJointOptions;
	DialogOpenGLOptions dialogOpenGLOptions;
	ComputeEigenmodes computeEigenmodes;
	CDialogProgressBar dialogProgressbar;
	
	CDrawWindow2DDlg Dialog2DView;
	TArray<CPlotToolDlg*> PlotToolArray;

	CDialogOneEditControl* sensorwatchctrl[MAX_SENSOR_WATCH];

	int fileopened;
	int model_modified;

	//TArray<int> asv_models; //$ DR 2013-05-21
	mystr asv_models;//$ DR 2013-05-21

	void ReadConfigFile();

// Construction
public:
	CWCDriver3DDlg(CWnd* pParent = NULL);	// standard constructor
	~CWCDriver3DDlg();

	WCDInterface * GetWCDInterface() { return pWCDI; }

	virtual void SetModelModified(int mod=1) {model_modified = mod;}
	virtual int GetModelModified() const {return model_modified;}

	CString m_ConsoleText_buffer;				// contains short string to append ( buffer for 
	CString m_ConsoleText_append;				// contains short string to append ( this is the stacked text
	int update_text_action_flag;				// choses routine to call in OutputDlg: 0 = full update, 1 = append, 2 = replace
	CRITICAL_SECTION uses_append_text;  // ensure that appendstring is not changed by AddText while processed by OnUpdateText

	virtual const CString& GetOutputDialogText() const {return m_ConsoleText;}

	//convert some internal settings (window position, open dialog windows, etc.) to EDC
	virtual void Configuration2EDC(ElementDataContainer& edc);
	//convert EDC to some internal settings (window position, open dialog windows, etc.)
	virtual void EDC2Configuration(const ElementDataContainer& edc);

// Dialog Data
	//{{AFX_DATA(CWCDriver3DDlg)
	enum { IDD = IDD_WCDRIVER3D_DIALOG };
	CString	m_ConsoleText;						// full output string
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWCDriver3DDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

public: //AD 2013-07-08: functions should be public so other Dialogs can use them (via pWCDI pointer)
// implementation of WCDInterface::UserInterface
	virtual void InstantMessageText(const char * );
	virtual void StatusText(const char * );
	virtual void AddText(const char * );
	virtual int CreateMissingDirectory(const char * );
	virtual void SleepX(int xMilliseconds);
	virtual void SetPlotToolRedrawTimer(int flag_onoff, int delay_ms = 1000);		//$ AD 2011-02 new timer: PlotToolRefreshTimer
//	virtual void OnTimer(UINT_PTR nIDEvent);
	//TODO
	virtual void AssureTextIsVisible() {}

// implementation of WCDInterface::ComputationFeedBack
	void ResultsUpdated(int flag);
	void FinishedComputation();
	void PausedComputation();
	bool RedrawNewResults() { return DialogDataManager.RedrawNewResults(); }

	virtual int CallWCDriverFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL);

	//add element; if option&1 -> element is added; returns 0 if dialog is canceled
	virtual int EditElementProperties(int nelem, int option = 0);

	virtual void GetSensorWatchDialog(CDialogOneEditControl*& watch, int& dialognum);
	
	// select (autocreate) a PlotToolDialog
	virtual int GetPlotToolDlg(CPlotToolDlg*& plottool, int number = 0);
	
	// create a new PlotToolDialog
	virtual int CreatePlotToolDlg(CPlotToolDlg*& plottool);

	// remove the PlotToolDlg from the list
	virtual int DeletePlotToolDlg(CPlotToolDlg* plottool);

protected:
  // display Rect for i-th PlotTool (Sensor Watch)
	CRect ComputeDisplayRectForSensorWatch(int i);

	//virtual void AddConstraint(int mode);
	virtual void RemoveExperimentalMenuItems();
	virtual void ModelChanged();
	virtual void DisplaySelectedModelName();
	virtual void AddRecentFile(const char * filename); 	//$ DR 2012-12-06 add "filename" to list of recent files and resort this list
	virtual void AutoSavembs();	//$ DR 2013-03-19

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CWCDriver3DDlg)
	virtual void OnInitMenuPopup(CMenu *pPopupMenu, UINT nIndex,BOOL bSysMenu);
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnButtonGo();
	afx_msg void OnButtonPause();
	afx_msg void OnPaint();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnClose();
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnTimer(UINT_PTR nIDEvent);

	//}}AFX_MSG
	afx_msg LRESULT OnUpdateText(WPARAM, LPARAM);
	afx_msg LRESULT OnButtonGoMessage(WPARAM, LPARAM);
  DECLARE_MESSAGE_MAP()


	bool bClosingDialog;
	void OnCancel();
public:
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnClickedButtonDrawres();
	afx_msg void OnBnClickedButtonDrawres2();
	afx_msg void OnBnClickedButtonTexres();
	afx_msg void OnBnClickedButtonTexres2();
	afx_msg void OnBnClickedButtonGeneraloptions();
	afx_msg void OnResultsDatamanager();
	afx_msg void OnResultsViewoutput();
	afx_msg void OnComputationReset();
	afx_msg void OnComputationPrinttimings();
	afx_msg void OnResultsShowoutputstatictext();
	afx_msg void OnOptionsSaveconfiguration();
	afx_msg void OnOptionsSaveSolverOptions();
	afx_msg void OnOptionsSaveHotintOptions();
	afx_msg void OnResultsEnableoutput();
	afx_msg void OnMove(int x, int y);
	afx_msg void OnViewFeDrawingoptions();
	afx_msg void OnViewViewingoptions();
//	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
	//afx_msg void OnComputationComputationoptions();	//$ DR 2013-10-16 removed
	afx_msg void OnComputationComputeEigenmodes();
	afx_msg void OnViewRigidbodyJointOptions();
	afx_msg void OnEditEditbody();
	//afx_msg void OnComputationAssemblesystem();
	afx_msg void OnViewOpengloptions();
	//afx_msg void OnComputationAssigninitialvector();
	afx_msg void OnFileExit();
	afx_msg void OnFileNewmbs();
	afx_msg void OnFileSavembs();
	afx_msg void OnEditEditsensor();
	afx_msg void OnFileOpenmbs();
	afx_msg void OnFileSavembsas();
	afx_msg void OnSpecialAddsensor();
	afx_msg void OnResultsPlotsensor();
	afx_msg void OnSpecialAddbody();
	afx_msg void OnResultsSensorwatch();
	afx_msg void OnFileRecentfile1();
	afx_msg void OnBnClickedButtonRecentFile();
	afx_msg void OnSpecialAddgeomelement();
	afx_msg void OnEditEditgeomelement();
	afx_msg void OnSystemShowsystemproperties();
	afx_msg void OnAddobjectAddnode();
	afx_msg void OnEditEditnode();
	afx_msg void OnSystemChecksystem();
	afx_msg void OnHelpAbout();
	afx_msg void OnViewDefaultview1();
	afx_msg void OnViewXY();
	afx_msg void OnViewXZ();
	afx_msg void OnViewYZ();
	afx_msg void OnComputationLoadinitialvector();
	afx_msg void OnComputationStoresolutionvector();
	afx_msg void OnMenuHelp();
	afx_msg void OnViewRobotoptions();
	afx_msg void OnEditEditmaterial();
public:
	afx_msg void OnSystemEditmodelparameters();
public:
	afx_msg void OnSystemLoadmodelparameter();
public:
	afx_msg void OnSystemSavemodeldata();
public:
	afx_msg void OnViewEditalloptions();
public:
	afx_msg void OnResultsPlot2sensorsxy();
public:
	afx_msg void OnResultsPlotnsensors();
public:
	afx_msg void OnViewShowprogressbar();
public:
	afx_msg void OnResultsPlottooldialog();
	afx_msg void OnBnClickedOptionsSaveconfiguration();
	afx_msg void OnFileSelectmodel();
	afx_msg void OnEditEditsolveroptions();
	afx_msg void OnEditEdithotintoptions();
	afx_msg void OnEditComputationSteps();
	afx_msg void OnEditEditLoad();
	afx_msg void OnAddload();
	afx_msg void OnViewShowContolWindow();
	afx_msg void OnAddobjectAddmaterial();
	afx_msg void OnAddbeam3dproperties();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnAddconnectorKinematicpairs();
	afx_msg void OnAddconnectorControlelement();
	afx_msg void OnBnClickedButtonSaveMbs();
	afx_msg void OnEditDeleteLoad();
	afx_msg void OnEditDeleteGeomelement();
	afx_msg void OnEditDeleteSensor();
	afx_msg void OnEditDeleteNode();
	afx_msg void OnEditDeleteMaterial();
	afx_msg void OnEditDeleteElement();
	afx_msg void OnEditUndo();
	afx_msg void OnSystemRunmacro();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_WCDriver3DDLG_H__270A972D_93F5_46BD_8737_55C336E60107__INCLUDED_)
