//#**************************************************************
//# filename:             DrawWindow2D.h
//#
//# author:               Gerstmayr, Dorninger
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
#include "resource.h"
#include "afxcmn.h"
#include "afxwin.h"
#include "MyBaseView.h" // base class for View

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CDrawWindow2DView:                                       derived View class for the 2D - PlotWindow ( Control elements )   * 2012-12-12 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// derived class for CDrawWindow2DView
class CDrawWindow2DView : public CMyBase2DView
{
	DECLARE_DYNCREATE(CDrawWindow2DView)

public:
	CDrawWindow2DView(CWnd* pParent = NULL);       // Standardkonstruktor
	virtual ~CDrawWindow2DView() {};
public:	
	virtual void Reset();                // initialize with default values ( from EDC or hardcoded )
	virtual void ErrorMessage(mystr& msg);

// ACCESS // LINKED DIALOGS & OBJECTS
private:
	class CDrawWindow2DDlg* m_pParent;             // pointer to parent required, DrawElementList is stored in DLG
public:
	void SetParentDlg(class CDrawWindow2DDlg* parent) { m_pParent = parent; }
	class CDrawWindow2DDlg* GetParentDlg() { return m_pParent; }

protected:
	DECLARE_MESSAGE_MAP()


public:
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// CATCH MOUSE directly
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);          // zoom in and out 
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);                    // start drag/drop 
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);							        // end drag/drop // catch click
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);                      // value - tracking 
  afx_msg void OnNcMouseMove(UINT nHitTest, CPoint point);

// CATCH OTHER via Msg
	void OnKeyPressed(MSG* pMsg);																							// catch key - catching must be done by DLG

	
// DRAG - DROP:    
	void DragDropAction();																										// drag drop action - move and redraw


// THE DRAWING PROCESS
public:
//	virtual	BOOL RedrawWindow(LPCRECT lpRectUpdate = NULL,	CRgn* prgnUpdate = NULL,UINT flags = RDW_INVALIDATE | RDW_UPDATENOW | RDW_ERASE);
	virtual void OnDraw(CDC* pDC);                          // this routine is automatically called by the framework (redrawwindow etc.)
	virtual void DrawElements(CDC* pDC);										// draw all items to specified CDC 

	virtual BOOL IsVisible(CDC* pDC, ControlWindowContext_::DrawComponent& elem);
	virtual BOOL DrawLine(CDC* pDC, ControlWindowContext_::DrawComponent& elem);
	virtual BOOL DrawRectangle(CDC* pDC, ControlWindowContext_::DrawComponent& elem);
	virtual BOOL DrawEllipse(CDC* pDC, ControlWindowContext_::DrawComponent& elem);
	virtual BOOL DrawText(CDC* pDC, ControlWindowContext_::DrawComponent& elem);

// is redraw allowed/prohibited ?
	virtual int DoNotRedraw();
	virtual void ForbidRedraw();
	virtual void AllowRedraw();
// StatusBarText
	virtual void UpdateStatusBarInfo();
	virtual void SetStatusBarText_Main (mystr text);
	virtual void SetStatusBarText_X (mystr text);
	virtual void SetStatusBarText_Y (mystr text);

  Box2D initshownrange;
	double currentzoomfactor;
	virtual double ReComputeZoomFactor();

	virtual void SetInitialRange();
	virtual void MatchLogicalPoints();
	virtual void Rescale();
	CRect prevClientRect;

// REAL VALUES of COORDINATES
	Box2D fullrange;   // full range of all lines to plot (real values)                
	virtual void ComputeFullRange(int flag_equal = 0); // computes limits of datasets - real values
  virtual void SetShownRange(double xmin, double xmax, double ymin, double ymax); // sets range of plot

	// Catching 
	virtual int CatchDrawObject( CPoint capture);									// computes the object to be selected
	int selected_object_buffer;																		
	virtual int GetDistanceFromComponent(int i, CPoint capture);
	virtual int DistFromLine(CPoint p1, CPoint p2, CPoint pmouse); // replace this function with MinDispLP from linalg3d
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CDrawWindow2DDlg:                               this is the dialog that nests the Graph (CMyView object),                  * 2012-12-12 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// DrawWindow2DDlg-Dialogfeld

class CDrawWindow2DDlg : public MyBaseDialogForView, // provides the View::OnDraw, Zoom, etc.
	                       public ControlWindowContext_ // provides the list of Elements to be drawn
{
	enum { IDD = IDD_DIALOG_2D_VIEW };
	DECLARE_DYNAMIC(CDrawWindow2DDlg)

public:
	CDrawWindow2DDlg(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CDrawWindow2DDlg();
public:
	virtual BOOL OnInitDialog();
	virtual BOOL DestroyWindow();

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	DECLARE_MESSAGE_MAP()

// ACCESS // LINKED DIALOGS & OBJECTS
public:
	void SetView(class CDrawWindow2DView* myView) { m_pMyView = (CDrawWindow2DView*) myView; }
	class CDrawWindow2DView* GetView() { return (class CDrawWindow2DView*) m_pMyView; } // we know that the derived dialog also uses a derived view class

	WCDInterface* pWCDI;
	//CGLDrawWnd* pGLDrawWnd;
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	WCDInterface * GetWCDI() { return pWCDI; }
	//void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }


// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
public:
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);  
	afx_msg void OnMove(int x, int y);
	afx_msg void OnClose();
	afx_msg void OnPaint();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	
// catch keydown
	BOOL PreTranslateMessage(MSG* pMsg);
	virtual void ForceUpdate();

// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
	void PlaceElements(); 

	// CONTEXTMENU
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	afx_msg void OnCMSelectobject();
	afx_msg void OnCMZoomToFullRange();
	afx_msg void OnCMExportToFile();
	afx_msg void OnCMPrintGraph();
	afx_msg void OnCMShowStatusInfo();
	int m_flag_statusinfo;

	// functions to create a new construction node in the middle of the selected connection line
	virtual Vector2D ComputeSplitPointPosition(ControlWindowContext_::DrawComponent& elem); // new position in the middle
	virtual int ComputeSplitPointListIndex(ControlWindowContext_::DrawComponent& elem);      // list index to insert the new node

	// Diagnose functions - Functions that trigger redraw
	void SetCaller(mystr calleri) { GetView()->SetCaller(calleri); }
};