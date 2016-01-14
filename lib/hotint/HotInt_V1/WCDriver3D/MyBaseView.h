//#**************************************************************
//# filename:             MyBaseView.h
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
 

#pragma once
#include "resource.h"
#include "afxcmn.h"
#include "afxwin.h"
#ifndef MFEMATH__H
  #include "..\MBSKernelLib\femath.h"
#endif
//#include "PlotToolDlg.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyMemoryDC:       derived View class - base class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Memory device context

class CMyMemDC : public CDC
{
public: 
	CMyMemDC(CDC* pDC, const CRect* pRect = NULL) : CDC()
	{
		ASSERT(pDC!=NULL);

		m_pDC = pDC;
		m_BMPold = NULL;
		_isMemDC = !(pDC->IsPrinting());
		if(pRect==NULL) { pDC->GetClipBox(&m_DrawRectLogical); }
		else { m_DrawRectLogical = *pRect;	}

		if(_isMemDC)
		{
			CreateCompatibleDC(pDC);

			// create Bitmap
			CRect m_DrawRectPixels = m_DrawRectLogical;
			pDC->LPtoDP(&m_DrawRectPixels);							// convert logical -> pixel coordinates for the bitmap
			m_BMPoffscreen.CreateCompatibleBitmap(pDC, m_DrawRectPixels.Width(), m_DrawRectPixels.Height());
			m_BMPold = SelectObject(&m_BMPoffscreen);

			// properties of original DC
			SetMapMode(pDC->GetMapMode());
			SetWindowExt(pDC->GetWindowExt());
			SetViewportExt(pDC->GetViewportExt());
			SetWindowOrg(m_DrawRectLogical.left, m_DrawRectLogical.top);
		}
		else
		{
			m_bPrinting = pDC->m_bPrinting;
			m_hDC       = pDC->m_hDC;
			m_hAttribDC	= pDC->m_hAttribDC;
		}
		FillSolidRect(m_DrawRectLogical, pDC->GetBkColor());
	}

	~CMyMemDC()
	{
		if(_isMemDC)
		{
			//copy the Bitmap to the screen
			m_pDC->BitBlt(m_DrawRectLogical.left,    m_DrawRectLogical.top,
				            m_DrawRectLogical.Width(), m_DrawRectLogical.Height(),
										this, m_DrawRectLogical.left,    m_DrawRectLogical.top, SRCCOPY);
			SelectObject(m_BMPold);
		}
	}



private:
	CBitmap m_BMPoffscreen;                   // offscreen bitmap
	CBitmap* m_BMPold;                        // previous bitmap
	CDC* m_pDC;                               // valid pDC
	CRect m_DrawRectLogical;                  // drawing area
	int _isMemDC;                             // true if drawing, false if printing
};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyBase2DView:       derived View class - base class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// functionality in base Class:  
// * Zoom (both axes at once) & Move via drag drop
// * Draw (Lines, Rectangle, Arrow, CheckVisibility, grid?)
// * Text placement

class CMyBase2DView : public CView
{
	DECLARE_DYNCREATE(CMyBase2DView)

public:
	CMyBase2DView(CWnd* pParent = NULL);       // Standardkonstruktor
	virtual ~CMyBase2DView() {};
public:	
	virtual void Reset();                // initialize with default values ( from EDC or hardcoded )
	virtual void ErrorMessage(mystr& msg) {};

protected:
	DECLARE_MESSAGE_MAP()

public:
// COORDINATE SYSTEM CONVERSION   pixels <--> logical
	CPoint GetPixelsFromLogical(CPoint logical, mystr& caller = mystr());
	CPoint GetLogicalFromPixels(CPoint pixels, mystr& caller = mystr());
// COORDINATE SYSTEM CONVERSION   real <--> logical 
	virtual CPoint ValuesToLogical(double x, double y);     // map values to logical coordinates
	virtual double XofLogical(CPoint& logical); 
	virtual double YofLogical(CPoint& logical);


public:
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// CATCH MOUSE
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);          // zoom in and out 
	virtual void ReComputeShownRange_Zooming(UINT nFlags, short zDelta, CPoint pt); // Zooming centered around current mouse position
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);                    // start drag/drop 
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point); // center view       // end drag/drop 
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);                      // value - tracking 
	int m_mouseisinclient;

public:
// DRAG - DROP:    
	virtual void DragDropAction();																						// drag drop action - move and redraw
	CRect prevdragrect;
		CPoint dragdrop_start_device;
	CPoint dragdrop_final_device;
	CPoint mouse_pos;
	CPoint& DragDropStart() { return dragdrop_start_device; }
	CPoint& DragDropFinal() { return dragdrop_final_device; }
	CPoint& MousePosition() { return mouse_pos; }
	int skipdrawduringmove;

public:
// LOGICAL COORDINATES - A plot window with border
                    /* Zero-----+----------------------------------+-----+ */
                    /* |                                                 | */
  int pts_y_top;    /* 400: Title                                        |
                    /* |                    *TITLE*                      | */
                    /* +        +----------------------------------+     | */
                    /* |        |                           legend?|     | */
                    /* |        |                                  |     | */
 	int pts_y_plot;   /* 2000:Plot|       * THE PLOT *               |     | */
                    /* |        |                                  |     | */
                    /* |        |                                  |     | */
                    /* +        Origin-----------------------------+     | */
                    /* |         (0,0)      *axis*                       | */
  int pts_y_bottom; /* 300: Axis                                         | */
                    /* |                                                 | */
                    /* +--------+----------------------------------+-----+ */
                    /*  300:Axis                 3000: Plot          200   */
                  	int pts_x_left;     	int pts_x_plot;       int pts_x_right;

public:
	// coordinates of "zero" in respect to "origin"	
	virtual void GetWindowOrigin(CPoint& zero) { zero.SetPoint(-pts_x_left,pts_y_top+pts_y_plot);	}
	// coordinates of rectangle for "the plot" in respect to "origin"
	virtual void GetPlotRect(CRect& rect) { rect.left = 0; rect.bottom = 0;	rect.right = pts_x_plot; rect.top = pts_y_plot;	}
	virtual void GetXAxisRect(CRect& rect) { rect.left = 0;	rect.bottom = -pts_y_bottom; rect.right = pts_x_plot;	rect.top = 0;	}
	virtual void GetYAxisRect(CRect& rect) { rect.left = -pts_x_left; rect.bottom = 0; rect.right = 0; rect.top = pts_y_plot;	}
	virtual void GetTitleRect(CRect& rect) { rect.left = 0; rect.bottom = pts_y_plot;	rect.right = pts_x_plot; rect.top = pts_y_plot+pts_y_top;	}
	virtual void GetViewRect(CRect& viewrect)	{	viewrect.left = -pts_x_left; viewrect.bottom = -pts_y_bottom; viewrect.right = pts_x_plot+pts_x_right; viewrect.top = pts_y_plot+pts_y_top; }
	
public:
// THE DRAWING PROCESS
	virtual void OnDraw(CDC* pDC) {;}														// this routine is automatically called by the framework (redrawwindow etc.)
	int _is_executing_ondraw;

// Diagnose functions 
	mystr caller;																								// this string contains the caller function for all redraws - default "system", changed right before all RedrawWindow calls
	virtual void SetCaller(mystr calleri) { caller = calleri; }
	virtual mystr& GetCaller() { return caller; }

// TEXT PLACEMENT
	virtual void SelectDefaultFont(CDC* pDC, int size = 12, const char* name = "Arial"); // a default font ( Arial )
	virtual CFont* GetDefaultFont(CDC* pDC, int size = 12, const char* name = "Arial"); // a default font ( Arial )
	virtual void DrawTextAtPoint(CDC* pDC, const CString str, CPoint point, double angle_deg); 	// plot text at a given posiiton (this position is center of textbox)
	virtual void DrawTextInRect(CDC* pDC, const CString str, CRect rect, double angle_deg, TTextAllign allign = TTextAllign (HCenter+VCenter), UINT nOptions = 0);	// plot text at in given rectangle (text mayhave allignment)
	virtual void DrawTextAt(CDC* pDC, const CString str, CPoint center, CPoint vector_center2origin, double angle_deg, UINT nOptions = 0);	// this function actually writes the text
// VISIBILITY
	virtual int IsPointInRect(CPoint& point, CRect& rect); // computes if point is within rectangle (!including borders which CRect.PtInRect does not)
	virtual int IsPointOnRectBorder(CPoint& point, CRect& rect); // computes if point is exactly on the rectangles border

public:
// FUNCTIONS THAT PREPARE THE DRAWING DEVICE CONTENTS
	virtual void OnFilePrint() {	CView::OnFilePrint(); }
	virtual void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = NULL);
	virtual void PrepareDC_Screen(CDC* pDC, CPrintInfo* pInfo = NULL);              // client is actual CView
//	virtual void PrepareDC_Memory(CDC* pDC, CPrintInfo* pInfo = NULL);              // client is Rect with Pixels as output BMP - only in derived
	virtual void PrepareDC_Printer(CDC* pDC, CPrintInfo* pInfo = NULL);             // client is determined by Printer device
	virtual void PrepareDC_inRect(CDC* pDC, CPrintInfo* pInfo, CRect& clientRect);  // scales DC

public:
// FUNCTIONS FOR PRINTING
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

public:
// FILE EXPORT FUNCTIONS
	virtual void SaveAsBitmapGraphic(CString& path_and_filename_noext, int resX, int resY, int flag_jpg, int flag_png, int flag_bmp); // jpg, png, bmp
  virtual void SaveAsVectorGraphic(CString& path_and_filename_noext, int flag_emf); // emf

public:
// ranges
	int Get_XRange_logical() { return pts_x_plot; }
	int Get_YRange_logical() { return pts_y_plot; }
// range of values that is shown in viewport
	Box2D shownrange;         // must be double values
	double Get_XMin() { return shownrange.PMin().X(); } //returns lower x limit of view port
	double Get_XMax() { return shownrange.PMax().X(); } //returns upper x limit of view port
	double Get_YMin() { return shownrange.PMin().Y(); } //returns lower y limit of view port
	double Get_YMax() { return shownrange.PMax().Y(); } //returns upper y limit of view port
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class MyBaseDialogForView:       derived View class - base class for 2D Views                                                    * 29.11.2012 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MyBaseDialogForView-Dialogfeld
class MyBaseDialogForView : public CDialog
{
public:
	CWnd* CreateNewView(CCreateContext* pContext, CWnd *pParent, CRect& rect, int wID);

	enum { IDD = IDD_DIALOG_2D_VIEW };
	DECLARE_DYNAMIC(MyBaseDialogForView)

public:
	MyBaseDialogForView(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~MyBaseDialogForView();
public:
	virtual BOOL OnInitDialog();
	virtual BOOL DestroyWindow();

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	DECLARE_MESSAGE_MAP()

	// ACCESS // LINKED DIALOGS & OBJECTS
protected:
	class CMyBase2DView* m_pMyView; // must be derived from CView
	class CWCDriver3DDlg* m_wcdd;   // owning WCDriver3Dlg instance
public:
	void SetView(CMyBase2DView* myView) { m_pMyView = myView ; }
	void SetWCDDlg(CWCDriver3DDlg* wcdd) {m_wcdd = wcdd; }
	CMyBase2DView* GetView() { return m_pMyView; }
	CWCDriver3DDlg* GetWCDDlg() { return m_wcdd; }

// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
public:
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);        // defines the minimum size for the dialog - important to not encounter negative extensions for the nested view
	afx_msg void OnMove(int x, int y);
	afx_msg void OnClose();
	afx_msg void OnPaint();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	
	BOOL PreTranslateMessage(MSG* pMsg);
	
	int m_init_done;     // sometimes reqired during initializaiton when Statusbar is involved
	int slow_zoom;
	int _isclosing;

// actively prohibit redraw
private:
	int flag_no_redraw;
public:
	int DoNotRedraw() { return flag_no_redraw; }
	void ForbidRedraw() { flag_no_redraw = 1; }
	void AllowRedraw() { flag_no_redraw = 0; }

// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
	virtual void PlaceElements(); 

	virtual CRect PlaceView(CRect& viewrect); // places the View Rect in the Popup Dialog
	virtual void CheckScreenLocation();       // ensure that dialog is in visible range
	virtual void GetViewRect(CRect& viewrect) { GetClientRect(viewrect); if (m_flag_use_statusbar) { viewrect.bottom = viewrect.bottom - statusbar_height; } } 
	virtual void GetStatusRect(CRect& statusrect) { GetClientRect(statusrect); statusrect.top = statusrect.bottom - statusbar_height; }	

// status bar at the bottom of the graph with a slot for each line (shows last value)
	CStatusBarCtrl m_StatusBar;
	mystr statusbartextbuffer;
	//virtual void SetStatusBarText( mystr text, int slot) { m_StatusBar.SetText(text,slot,0); }
	//virtual void SetStatusBarText_Main (mystr text) {;}
	//virtual void SetStatusBarText_X (mystr text) {;}
	//virtual void SetStatusBarText_Y (mystr text) {;}
	int m_flag_use_statusbar;
	static const int statusbar_height = 25;
	virtual int& SetStatusBar(int flag_use) { return m_flag_use_statusbar; }

// name displayed in the window title bar
	CString def_name;
	void SetDefaultName(CString& defaultname) { def_name = defaultname; }
	CString& GetDefaultName() { return def_name; }
	void ReadDefaultName(CString& defaultname) { defaultname = def_name; }
	void SetDisplayName(CString& displayname) { this->SetWindowTextA(displayname); }
	void ReadDisplayname(CString& displayname) { this->GetWindowText(displayname); }

// FUNCTIONS FOR PRINTING
	void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo) { return; }
	void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo) { return; }
	void Print(); 

	virtual void SetCaller(mystr calleri) { GetView()->caller = calleri; }
};
