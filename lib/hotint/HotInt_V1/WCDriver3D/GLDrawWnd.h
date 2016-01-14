//#**************************************************************
//# filename:             GLDrawWnd.h
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
 


#if !defined(AFX_GLDRAWWND_H__C4A1F1E3_E308_11D4_BA16_000021038B9A__INCLUDED_)
#define AFX_GLDRAWWND_H__C4A1F1E3_E308_11D4_BA16_000021038B9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// GLDrawWnd.h : header file
//

#include <afxtempl.h>

#include "rendercontext.h"
#include "dialogframesrecording.h"

#include "extgl.h"
//#include "gl/gl.h"
#include "gl/glu.h"

#define WM_REDRAW (WM_USER + 1)

/////////////////////////////////////////////////////////////////////////////
// CGLDrawWnd window

class CGLDrawWnd : public CWnd, RenderContext_
{
	float saved_modelview_matrix[16];

	// flag: no rotation
	bool bNoRotationState;

	void MakeZoom(const CPoint & p1, const CPoint & p2);

	struct TextPortion
	{
		enum {Text3D, Text2D, TextStruct} type;
		float x,y,z;
		int nLineNo;
		int nXPos;
		CString text;
	};
	// in this list the texts are being stored
	CList<TextPortion,TextPortion&> Texts;
	float GetStructTextVertPos(int nLineNo);
	float GetStructTextHorPos(int nXPos,const CString & text);

	static CToolBar ToolBar;

	GLuint GL_SCENE_LIST;


	void Init();

	BOOL SetupPixelFormat();
	void StopOpenGL();
	void SetPerspective();

	CDC * m_pDC;
	HGLRC hrc;

	/*
	static GLsizei SelBufSize;
	void SelectRegion(const CPoint & p1, const CPoint & p2, bool bShift, bool bCtrl);
	CList<unsigned int,unsigned int&> Selection;
	*/

	// Viewing parameters
	float TranslationStepCoeff;		// how quick is the view pan
	float RotationStep;				// how quick is rotation
	float DetailLevelStepCoeff;		// how quick is zoom
	float PerspectiveStepCoeff;		// how quick is perspective
	GLdouble AspectRatio;			// currect aspect ratio

	enum 
	{
		ActionNone,
		ActionSelect,
		ActionMove,
		ActionRotate,
		ActionZoom,
		ActionPerspective
	} Action;
	CPoint CurrentPoint;	// here is the mouse
	CPoint StartingPoint;	// here it was as the operation started

	Vertex CenterPoint;		// here is the actual center of the construction
	Vertex CenterOffset;		// offset to the center of the construction (camera moving with object)
	float SCENE_OFFSET_COEFF;		// see SetPerspective()


	void DrawScene();

	int AxesPosition;		// 5 -- no axes
	bool bFittingTheView;	// is used in OnButtonFit(),
							// so that some details would not affect the scene

	bool bRotationTimerRunning;

	WCDInterface * pWCDI;

	int bShowLegend;

	CDialogFramesRecording DialogFramesRecording;
	void PerformVideoFramesRecord();


	// implementing the rest of the interface
	virtual void PrintText2D(float x, float y, const char * text);
	virtual void PrintText3D(float x, float y, float z, const char * text);
	virtual void PrintTextStruct(int nLineNo, int nXPos, const char * text);
	virtual int GetDetailLevel() { return 1; }
	virtual int	ShowLegend() { return bShowLegend; }
	virtual int GetConfigurationNumber() { return 1; }

	float BackgroundR,BackgroundG,BackgroundB;
	virtual void SetBackgroundColor(float r, float g, float b)
	{
		BackgroundR = r;
		BackgroundG = g;
		BackgroundB = b;
	}

	float ScrTextR,ScrTextG,ScrTextB;
	virtual void SetTextColor(float r, float g, float b)
	{
		ScrTextR = r;
		ScrTextG = g;
		ScrTextB = b;
	}
	virtual void SetCenterOffset(float x, float y, float z)
	{
		CenterOffset.x = x;
		CenterOffset.y = y;
		CenterOffset.z = z;
	}

	virtual void GetWindowSize(int& width, int& height)
	{
		CRect r;
		GetWindowRect(r);
		width = r.Width();
		height = r.Height();
	}

	virtual void SetSceneOffsetCoeff(float x);


	void PrintTextsFixed();

	// chooses the color for the axes to be painted
	void ChooseAxesColor();

	// these functions can be used in order to move the clipping plane, before
	// painting the lines or points that can appear in a surface of a solid
	virtual void BeginLinesOnSolid();
	virtual void EndLinesOnSolid();

	void CreateFEColorTexture(int ncols, int texnum, int grey); //create 1D texture for FE-color drawing

	CString ProgramDirectory;	// directory where the program is running

// Construction
public:
	CGLDrawWnd();

	void SetAnimateScalingTimer(int flag);
	bool prohibit_redraw;

	double old_scaling_factor;
	double virtual_animate_time;

	// the class should be not constructed without setting this pointer
	void SetWCDI(WCDInterface * p_wcdi);
	
	void Redraw();
	void ButtonFit();
	void ResetOpenGLParam(); //call if global light or material has changed
	void ContentsChanged(int forcefit=0); //Redraw and recompute scene size

	void SaveConfig(CArchive & ar);
	void LoadConfig(CArchive & ar);

	//convert some internal settings (window position, open dialog windows, etc.) to EDC
	virtual void Configuration2EDC(ElementDataContainer& edc);
	//convert EDC to some internal settings (window position, open dialog windows, etc.)
	virtual void EDC2Configuration(const ElementDataContainer& edc);

	void SetProgramDirectoryAsCurrent() { SetCurrentDirectory(ProgramDirectory); }
	const CString& GetProgramDirectory() { return ProgramDirectory; }

	BOOL & Lighting() { return pWCDI->GetIOption(206)/*bLighting*/; } //gl_lighting

	void OnMouseWheelGLDW(UINT nFlags, short zDelta, CPoint pt);

// Attributes
public:
	//HENHMETAFILE hEnhMetaFile;

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CGLDrawWnd)
	public:
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CGLDrawWnd();

	// Generated message map functions
//protected:
public:
	//{{AFX_MSG(CGLDrawWnd)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnDestroy();
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnButtonNoRotation();
	afx_msg void OnButtonStandardViewXy();
	afx_msg void OnButtonStandardViewXyz();
	afx_msg void OnButtonStandardViewXz();
	afx_msg void OnButtonStandardViewYz();
	afx_msg void OnButtonStandardViewRotate();
	afx_msg void OnButtonFit();
	afx_msg void OnButtonMoveAxes();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnButtonSaveImage();
	afx_msg void OnButtonFramesRecording();
	//}}AFX_MSG
	afx_msg LRESULT OnRedraw(WPARAM, LPARAM);
	BOOL ToolBarToolTipsSupport( UINT id, NMHDR * pTTTStruct, LRESULT * pResult );
	DECLARE_MESSAGE_MAP()
public:
//	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_GLDRAWWND_H__C4A1F1E3_E308_11D4_BA16_000021038B9A__INCLUDED_)
