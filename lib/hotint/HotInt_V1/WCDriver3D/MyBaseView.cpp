//#**************************************************************
//# filename:             MyBaseView.cpp
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
#include "MyBaseView.h"
#include "savewindowbitmap.h"



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyBase2DView:       derived View class - base class for 2D Views                                                          * 29.11.2012 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLEMENT_DYNCREATE(CMyBase2DView, CView)

CMyBase2DView::CMyBase2DView(CWnd* parent)
{
	skipdrawduringmove = 0;
	_is_executing_ondraw = 0;
}

void CMyBase2DView::Reset()
{
	prevdragrect = CRect(0,0,0,0);
	caller = mystr("REDRAW triggered by SYSTEM\n");
	m_mouseisinclient = 0;
	skipdrawduringmove = 0;
}

BEGIN_MESSAGE_MAP(CMyBase2DView, CView)
	ON_WM_MOUSEWHEEL()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
END_MESSAGE_MAP()


// *****************************************************
// // COORDINATE SYSTEM CONVERSION   pixels <--> logical
// *****************************************************
CPoint CMyBase2DView::GetPixelsFromLogical(CPoint capture_LOGICAL, mystr& caller)
{
	CDC* pDC = this->GetDC(); 
	//AD: CATCH : sometimes the pDC is returned as null and the application crashes... inline void* CThreadSlotData::GetThreadValue(int nSlot)
	if(!pDC)                           // bailout in case the device context is not valid
	{
		ErrorMessage ( mystr("Illegal pDC - ") + caller );
		return 0;      
	}
// coordinates
	OnPrepareDC(pDC);                  // mapping function
  CPoint capture_PIXELS = capture_LOGICAL;
	pDC -> LPtoDP(&capture_PIXELS);    // mouse position in Logical Coordinates
	this->ReleaseDC(pDC);
	return capture_PIXELS;
}

CPoint CMyBase2DView::GetLogicalFromPixels(CPoint capture_PIXELS, mystr& caller)
{
	CDC* pDC = this->GetDC(); 
	//AD: CATCH : sometimes the pDC is returned as null and the application crashes... inline void* CThreadSlotData::GetThreadValue(int nSlot)
	if(!pDC)                           // bailout in case the device context is not valid
	{
		ErrorMessage ( mystr("Illegal pDC - ") + caller );
		return 0;      
	}
// coordinates
	OnPrepareDC(pDC);                   // mapping function
  CPoint capture_LOGICAL = capture_PIXELS;
	pDC -> DPtoLP(&capture_LOGICAL);    // mouse position in Pixel Coordinates
	this->ReleaseDC(pDC);
	return capture_LOGICAL;
}
// ************************************************
// COORDINATE SYSTEM CONVERSION   real <--> logical 
// ************************************************
CPoint CMyBase2DView::ValuesToLogical(double x, double y)
{
	CPoint xy_logical;
	double xrange = Get_XMax() - Get_XMin();
	double yrange = Get_YMax() - Get_YMin();
	xy_logical.x = (long) (((x - Get_XMin())/xrange) * Get_XRange_logical());
	xy_logical.y = (long) (((y - Get_YMin())/yrange) * Get_YRange_logical());
	return xy_logical;
}

double CMyBase2DView::XofLogical(CPoint& logical)
{
	return logical.x * (Get_XMax() - Get_XMin()) / (double) Get_XRange_logical() + Get_XMin();
}

double CMyBase2DView::YofLogical(CPoint& logical)
{
	return logical.y * (Get_YMax() - Get_YMin()) / (double) Get_YRange_logical() + Get_YMin();
}

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
BOOL CMyBase2DView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	// Zoom according to additional key pressed
	ReComputeShownRange_Zooming(nFlags, zDelta, pt);

	return 1; // otherwise zoom diverges
}

// Zooming centered around current mouse position
void CMyBase2DView::ReComputeShownRange_Zooming(UINT nFlags, short zDelta, CPoint pt)
{
	CPoint pixels = pt;
	ScreenToClient(&pixels);    // convert coordinates (screen) to client coordinates
	
	CRect cr;
	GetClientRect(cr);

	if (cr.PtInRect(pixels))
	{
		CPoint logical = GetLogicalFromPixels(pixels, mystr("PlotToolView::OnMouseWheel")); 

		double x = XofLogical(CPoint(logical.x,logical.y));
		double y = YofLogical(CPoint(logical.x,logical.y));

		// scaling factor
		double scale_tick = 1.2;
		int alt_down = GetKeyState(VK_MENU) & 0x80; // slow zoom when ALT key is down
		if(alt_down)
			scale_tick = 1.04;
		double scale_total = pow(scale_tick,(double) -zDelta/120); // one wheel tick is 120 units, is 20% zoom

		double xmin, xmax, ymin, ymax;
		if( (nFlags & MK_CONTROL) && !(nFlags & MK_SHIFT)) // only CTRL is down -> zoom Y
		{
			xmin = Get_XMin();
			xmax = Get_XMax();
			ymin = y + (Get_YMin() - y) * scale_total;
			ymax = y + (Get_YMax() - y) * scale_total;
		}
		else if( !(nFlags & MK_CONTROL) && (nFlags & MK_SHIFT)) // only SHIFT is down -> zoom X
		{
			xmin = x + (Get_XMin() - x) * scale_total;
			xmax = x + (Get_XMax() - x) * scale_total;
			ymin = Get_YMin();
			ymax = Get_YMax();
		}
		else // zoom both directions
		{
			xmin = x + (Get_XMin() - x) * scale_total;
			xmax = x + (Get_XMax() - x) * scale_total;
			ymin = y + (Get_YMin() - y) * scale_total;
			ymax = y + (Get_YMax() - y) * scale_total;
		}

		// Set the shown range
		shownrange = Box2D(Vector2D(xmin, ymin),Vector2D(xmax, ymax));
	}
	return; 
}

void CMyBase2DView::OnLButtonDown(UINT nFlags, CPoint point)
{
// should no longer be reached...

// standard drag drop - drag drop start
	DragDropStart() = point;
	CView::OnLButtonDown(nFlags, point);
}

void CMyBase2DView::OnLButtonUp(UINT nFlags, CPoint point)
{
// should no longer be reached...

	// standard drag drop - drag drop end
	DragDropFinal() = point;
	DragDropAction();
	SetCaller(mystr("REDRAW triggered by CMyBase2DView::OnLButtonUp\n"));
	RedrawWindow();
	CView::OnLButtonUp(nFlags, point);
}

void CMyBase2DView::OnMouseMove(UINT nFlags, CPoint point)
{
	if ((nFlags & MK_LBUTTON) && !(nFlags & MK_SHIFT)) // left mouse button is down, shift key is up --> dragdrop --> animate
	{
		DragDropFinal() = point; 
		// check if distance for redraw is reached
		double moved_x = DragDropFinal().x - DragDropStart().x;
		double moved_y = DragDropFinal().y - DragDropStart().y;
		double dragdropdistance = sqrt(moved_x*moved_x + moved_y*moved_y);

		if (dragdropdistance > 10.)
		{
			DragDropAction();
			SetCaller(mystr("REDRAW triggered by CMyBase2DView::OnMouseMove\n"));
			RedrawWindow();
			DragDropStart() = point;
		}
	}

// MORE IN DERIVED CLASS

	CView::OnMouseMove(nFlags, point);
}

// drag drop action - shift and redraw
void CMyBase2DView::DragDropAction()
{
// MOVE ENTIRE "SHEET"
	CPoint logical_start = GetLogicalFromPixels(DragDropStart(),mystr("DragDropAction - StartPoint"));
	CPoint logical_final = GetLogicalFromPixels(DragDropFinal(),mystr("DragDropAction - FinalPoint"));

	double x = XofLogical(logical_final) - XofLogical(logical_start);
	double y = YofLogical(logical_final) - YofLogical(logical_start);

	// rescale
	double xmin = Get_XMin() - x;
	double xmax = Get_XMax() - x;
	double ymin = Get_YMin() - y;
	double ymax = Get_YMax() - y;

	// ranges
	shownrange = Box2D(Vector2D(xmin, ymin),Vector2D(xmax, ymax));
	// MORE IN DERIVED CLASS
}

// **************
// TEXT PLACEMENT
// **************
void CMyBase2DView::SelectDefaultFont(CDC* pDC, int size, const char* name)
{
	CFont* defaultfont = GetDefaultFont(pDC, size, name);
	pDC->SelectObject(defaultfont);
}

CFont* CMyBase2DView::GetDefaultFont(CDC* pDC, int size, const char* name)
{
	// default font
	LOGFONT logfont;
	logfont.lfHeight = -14 * size;
	//logfont.lfHeight = -160; // ~ size 12
	logfont.lfWidth = 0;
	logfont.lfEscapement = 0;
	logfont.lfOrientation = 0;
	logfont.lfWeight = 400;
	logfont.lfItalic = 0;
	logfont.lfUnderline = 0;
	logfont.lfStrikeOut = 0;
	logfont.lfCharSet = 0;
	logfont.lfOutPrecision = 3;
	logfont.lfClipPrecision = 2;
	logfont.lfQuality = 1;
	logfont.lfPitchAndFamily = 34;
	strcpy(logfont.lfFaceName,name);
	//strcpy(logfont.lfFaceName,"Arial");
	
	CFont* defaultfont = new CFont;
	defaultfont->CreateFontA(logfont.lfHeight,logfont.lfWidth,logfont.lfOrientation,logfont.lfOrientation,logfont.lfWeight,logfont.lfItalic,logfont.lfUnderline,
			logfont.lfStrikeOut,logfont.lfCharSet,logfont.lfOutPrecision,logfont.lfClipPrecision,logfont.lfQuality,logfont.lfPitchAndFamily,logfont.lfFaceName);
  return defaultfont;
}

void CMyBase2DView::DrawTextAtPoint(CDC* pDC, const CString str, CPoint point, double angle_deg)
{
	DrawTextInRect(pDC,str,CRect(point,point),angle_deg);
}

void CMyBase2DView::DrawTextInRect(CDC* pDC, const CString str, CRect rect, double angle_deg, TTextAllign allign, UINT nOptions)
{
	// textbox in non rotated system
	CSize textsize = pDC->GetTextExtent(str);
	CSize rectsize(abs(rect.Width()),abs(rect.Height())); // height=top-bottom would be negative

	double angle_rad = angle_deg * 3.14159 / 180.0;
	// diagonal vectors of the rotated textbox
	CPoint vector_textbox_diag_ur2ll, vector_textbox_diag_ul2lr; 
	vector_textbox_diag_ur2ll.x = long (-textsize.cx * cos(angle_rad) + textsize.cy * sin(angle_rad));
	vector_textbox_diag_ur2ll.y = long (-textsize.cx * sin(angle_rad) - textsize.cy * cos(angle_rad));
	vector_textbox_diag_ul2lr.x = long (-textsize.cx * cos(angle_rad) - textsize.cy * sin(angle_rad));
	vector_textbox_diag_ul2lr.y = long (-textsize.cx * sin(angle_rad) + textsize.cy * cos(angle_rad));

	// calculate nonrotated rectangle that surrounds rotated textbox
	CSize used,unused;
	used.cx = abs(vector_textbox_diag_ur2ll.x);
	if (abs(vector_textbox_diag_ul2lr.x)>used.cx) used.cx=abs(vector_textbox_diag_ul2lr.x);
	used.cy = abs(vector_textbox_diag_ur2ll.y);
	if (abs(vector_textbox_diag_ul2lr.y)>used.cy) used.cy=abs(vector_textbox_diag_ul2lr.y);
	unused = rectsize-used;

	if( allign & ScaleDown2Fit ) // automatically scale font down to fit rectangle
	{
		if(unused.cx < 0 || unused.cy < 0)
		{
			double scaling_factor_x = 1.;
			if(unused.cx < 0)
			{
				scaling_factor_x = (double) (rectsize.cx) / (double) (used.cx);
			}

			double scaling_factor_y = 1.;
			if(unused.cy <0)
			{
				scaling_factor_y = (double) (rectsize.cy) / (double) (used.cy);
			}

			double scaling_factor = Minimum(scaling_factor_x,scaling_factor_y);

			LOGFONT logfont;
			// font must be scaled down to fit target rectangle
			pDC->GetCurrentFont()->GetLogFont(&logfont);
			logfont.lfHeight = (long) (logfont.lfHeight*scaling_factor);
			logfont.lfWeight = (long) (logfont.lfWeight*scaling_factor);
			logfont.lfWidth  = (long) (logfont.lfWidth *scaling_factor);
			CFont* scaledFont = new CFont();
			scaledFont->CreateFontA(logfont.lfHeight,logfont.lfWidth,logfont.lfOrientation,logfont.lfOrientation,logfont.lfWeight,logfont.lfItalic,logfont.lfUnderline,
				logfont.lfStrikeOut,logfont.lfCharSet,logfont.lfOutPrecision,logfont.lfClipPrecision,logfont.lfQuality,logfont.lfPitchAndFamily,logfont.lfFaceName);
			CFont* prevFont = pDC->SelectObject(scaledFont);

			// Catch if Font is too small 
			if (abs(logfont.lfHeight) >10 )
			{
			// second try - this time the fontsize will be correct 
				DrawTextInRect(pDC, str, rect, angle_deg, allign, nOptions);
			}
			else
			{
			// crashes occurred for lfHeight of -7
			// the font is too small anyway, skip printing text
			}
			pDC->SelectObject(prevFont);
			return;
		}
	}

	// compute center point of textbox - apply textallign
	CPoint center;
	int nr_h_pos = 3; // HRight-HLeft+1;
	int nr_v_pos = 3; // (VTop-VBottom) / (HRight-HLeft+1) + 1;
	int nr_h_center = 1; 
	int nr_v_center = 1;
	int h_allign = allign % (int) (nr_h_pos);  // 0.1.2 left.center.right
	int v_allign = allign / (int) (nr_h_pos);  // 0.1.2 top.center.bottom

	// here we assume that 
	center.x = rect.CenterPoint().x +												// center point of the placement rectangle  // "+" left < right       						 
			 (int) ( ((double)unused.cx/double(nr_h_pos-1)) *		// x-distance for a single degree of X-allignment ( number of horizontal intervals)
						   ((double)(h_allign-nr_h_center)) + 0.5 );	// number of shifts towards either direction from the center							
	
	center.y = rect.CenterPoint().y -												// center point of the placement rectangle  // "-" top < bottom						 
			 (int) ( ((double)unused.cy/double(nr_v_pos-1)) *		// y-distance for a single degree of Y-allignment ( number of vertical intervals)
						   ((double)(v_allign-nr_v_center)) + 0.5 );	// number of shifts towards either direction from the center			

	////center.x = rect.left + (int) (used.cx/2.0 + unused.cx*(h_allign-1)/(VCenter-1.) + 0.5);
	////center.y = rect.top  - (int) (used.cy/2.0 + unused.cy*(v_allign-1)/(VCenter-1.) + 0.5);

	// vector from center of text to origin of first letter
	CPoint vector_center2lowerleft;
	vector_center2lowerleft.x = (int) (vector_textbox_diag_ur2ll.x * 0.5);
	vector_center2lowerleft.y = (int) (vector_textbox_diag_ur2ll.y * 0.5);
	DrawTextAt(pDC,str,center,vector_center2lowerleft,angle_deg,nOptions);

	if(0) // hack - rectangle for text
	{
		CPen p(PS_DOT,5,RGB(0,0,0));
		CPen* prev = pDC->SelectObject(&p);
		pDC->MoveTo(rect.right,rect.top);
		pDC->LineTo(rect.left,rect.top);
		pDC->LineTo(rect.left,rect.bottom);
		pDC->LineTo(rect.right,rect.bottom);
		pDC->LineTo(rect.right,rect.top);
		pDC->SelectObject(prev);
	}
}

void CMyBase2DView::DrawTextAt(CDC* pDC, const CString str, CPoint textcenter, CPoint vector_center2origin, double angle_deg, UINT nOptions)
{
	CFont* actFont = pDC->GetCurrentFont();
	if (angle_deg!=0.0)
	{
		LOGFONT logfont;

		// font for rotated text
		pDC->GetCurrentFont()->GetLogFont(&logfont);
		logfont.lfOrientation = (long) (-angle_deg * 10.0 +0.5);
		CFont rotFont;
		rotFont.CreateFontA(logfont.lfHeight,logfont.lfWidth,logfont.lfOrientation,logfont.lfOrientation,logfont.lfWeight,logfont.lfItalic,logfont.lfUnderline,
			logfont.lfStrikeOut,logfont.lfCharSet,logfont.lfOutPrecision,logfont.lfClipPrecision,logfont.lfQuality,logfont.lfPitchAndFamily,logfont.lfFaceName);
		pDC->SelectObject(&rotFont);
	}

	CPoint textorigin;
	textorigin.x = textcenter.x + (int)(vector_center2origin.x + 0.5);
	textorigin.y = textcenter.y + (int)(vector_center2origin.y + 0.5);

	pDC->SetTextAlign(TA_BASELINE);
	pDC->SetBkMode(TRANSPARENT);
	pDC->ExtTextOut(textorigin.x, textorigin.y,
		nOptions, NULL, str, NULL);

	pDC->SelectObject(actFont); // back to original font
}

// **********
// VISIBILITY
// ***********
// computes if point is within rectangle (!including borders which CRect.PtInRect does not)
int CMyBase2DView::IsPointInRect(CPoint& point, CRect& rect)
{
	if(point.x < rect.left) return 0;
	if(point.x > rect.right) return 0;
	if(point.y < rect.bottom) return 0;
	if(point.y > rect.top) return 0;
	return 1;
}

// computes if point is exactly on the rectangles border
int CMyBase2DView::IsPointOnRectBorder(CPoint& point, CRect& rect)
{
	if( ((point.x == rect.left) || (point.x == rect.right)) && ((point.y >= rect.top) && (point.y <= rect.bottom)) ) return 1;
	if( ((point.y == rect.top) || (point.y == rect.bottom)) && ((point.x >= rect.left) && (point.x <= rect.right)) ) return 1;
	return 0;
}

// **************************************************
// FUNCTIONS THAT PREPARE THE DRAWING DEVICE CONTENTS
// **************************************************
void CMyBase2DView::OnPrepareDC(CDC* pDC, CPrintInfo* pInfo)
{
	if(!pDC->IsPrinting())
	{
		PrepareDC_Screen(pDC, pInfo);
	}
	CView::OnPrepareDC(pDC, pInfo);
}

// client is actual CView
void CMyBase2DView::PrepareDC_Screen(CDC* pDC, CPrintInfo* pInfo)
{
	CRect clientRect;
	GetClientRect(clientRect);
	PrepareDC_inRect(pDC, pInfo, clientRect);
}

// client is determined by Printer device
void CMyBase2DView::PrepareDC_Printer(CDC* pDC, CPrintInfo* pInfo)
{
	CRect clientRect = pInfo->m_rectDraw;
	PrepareDC_inRect(pDC, pInfo, clientRect);
}

// scales DC
void CMyBase2DView::PrepareDC_inRect(CDC* pDC, CPrintInfo* pInfo, CRect& clientRect)
{
	CRect L_ViewRect;
	CPoint L_Origin;

	this->GetViewRect(L_ViewRect);
	this->GetWindowOrigin(L_Origin);
	
	pDC->SetMapMode(MM_ISOTROPIC);
	pDC->SetWindowExt(L_ViewRect.Width(),L_ViewRect.Height());
	pDC->SetWindowOrg(L_Origin.x,L_Origin.y);
	
	pDC->SetViewportExt(clientRect.right-clientRect.left, clientRect.bottom-clientRect.top);
	pDC->SetViewportOrg(0,0);
}

// **********************
// FUNCTIONS FOR PRINTING
// **********************
void CMyBase2DView::OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo)
{
	CView::OnBeginPrinting(pDC, pInfo);
}

void CMyBase2DView::OnPrint(CDC* pDC, CPrintInfo* pInfo)
{
	PrepareDC_Printer(pDC, pInfo);
	CView::OnPrint(pDC, pInfo);
}

void CMyBase2DView::OnEndPrinting(CDC* pDC, CPrintInfo* pInfo)
{
	CView::OnEndPrinting(pDC, pInfo);
}

// *********************
// FILE EXPORT FUNCTIONS
// *********************
#include <gdiplus.h>
using namespace Gdiplus;
#include <atlimage.h>

void CMyBase2DView::SaveAsBitmapGraphic(CString& fn, int resX, int resY, int flag_jpg, int flag_png, int flag_bmp)
{
	// Coordinate Transformation
	CRect L_ViewRect;	this->GetViewRect(L_ViewRect);
	CPoint L_Origin;	this->GetWindowOrigin(L_Origin);
	CRect rect(0,0,resX, resY);

// Memory DC where the image is drawn
	CDC MemDC;
	MemDC.CreateCompatibleDC(NULL);
// set MapMode for device context
	PrepareDC_inRect(&MemDC,NULL,rect);
	MemDC.m_hAttribDC = MemDC.m_hDC;
// make CBitmap to draw to 
	CBitmap* pBitmap = new CBitmap();
	pBitmap->CreateBitmap(resX, resY, 1, MemDC.GetDeviceCaps(BITSPIXEL), NULL); // !! Device coordinates to create bitmap PIXELS
// initialize BMP for drawing
	MemDC.SelectObject(pBitmap);

// white backgroud color for entire bitmap region ( include the region that is not covered due to MM_ISOMETRIC
	MemDC.DPtoLP(rect);
	MemDC.PatBlt(rect.left, rect.top, rect.Width(), rect.Height(), WHITENESS);

	// draw into DC bitmap
	OnDraw(&MemDC);
	
	HBITMAP hBmp = (HBITMAP)pBitmap->GetSafeHandle();
	if(flag_jpg)
	{
	//	CImageSave(hBmp,fn+".jpg",ImageFormatJPEG);
		GDIPlusSave(hBmp,fn+".jpg",ImageFormatJPEG);
	}
	if(flag_png)
	{
	//	CImageSave(hBmp,fn+".png",ImageFormatPNG);
		GDIPlusSave(hBmp,fn+".png",ImageFormatPNG);
	}
	if(flag_bmp)
	{
	//	CImageSave(hBmp,fn+".bmp",ImageFormatBMP);
		GDIPlusSave(hBmp,fn+".bmp",ImageFormatBMP);
	}
}

void CMyBase2DView::SaveAsVectorGraphic(CString& fn, int flag_emf)
{
	// standard paint routine
	CPaintDC DC(this);
	CPaintDC* pDC = &DC;
	OnPrepareDC(pDC);

	// create metafile
	CMetaFileDC* pEmfDC = new CMetaFileDC;
	pEmfDC->CreateEnhanced(pDC, fn+".emf", NULL, NULL); // y-axis is swapped, chose MapMode
	pEmfDC->SetMapMode(MM_ISOTROPIC); 
	pEmfDC->m_hAttribDC = pEmfDC->m_hDC;

	OnDraw(pEmfDC);

//close metafile
	HENHMETAFILE hWMF = ((CMetaFileDC*)pEmfDC)->CloseEnhanced();
	::DeleteEnhMetaFile(hWMF);
	Invalidate();
}

#define ID_VIEW_ON_POPUP 1234
#define ID_MENU_PARENT_ON_POPUP 1235
#define ID_STATUSBAR_ON_POPUP 1236

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyBaseNestingDialog:       derived View class - base class for 2D Views                                                   * 29.11.2012 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// TODO: PlotToolDialogPopupChild could be derived
// TODO: derive 2D - IO - Plot
IMPLEMENT_DYNAMIC(MyBaseDialogForView, CDialog)

MyBaseDialogForView::MyBaseDialogForView(CWnd* pParent /*=NULL*/)
	: CDialog(MyBaseDialogForView::IDD, pParent)
{
}

MyBaseDialogForView::~MyBaseDialogForView()
{
}

BOOL MyBaseDialogForView::OnInitDialog()
{
	CDialog::OnInitDialog();

	// create the View
	CCreateContext cc;
	cc.m_pNewViewClass = RUNTIME_CLASS(CMyBase2DView);
	cc.m_pCurrentDoc = NULL;
	cc.m_pNewDocTemplate = NULL;
	cc.m_pLastView = NULL;
	cc.m_pCurrentFrame = NULL;

	m_pMyView = (CMyBase2DView*)CreateNewView(&cc, this, CRect(0, 0, 0, 0), ID_VIEW_ON_POPUP);
	if (m_pMyView == NULL)
		EndDialog(IDCANCEL);

	CRect rect;
	this->GetWindowRect(&rect);
	rect.top = rect.bottom- 25;
	int dbg = this->m_StatusBar.Create(WS_CHILD | WS_BORDER | WS_VISIBLE, rect, this, ID_STATUSBAR_ON_POPUP);

	GetView()->Reset();
	m_init_done = 1;
	slow_zoom = 0;
	_isclosing = 0;
	AllowRedraw();
	statusbartextbuffer = mystr("");

	SetCaller(mystr("REDRAW triggered by MyBaseDialogForView::OnInitDialog\n"));
	this->UpdateWindow();
	this->SetDisplayName(CString("2DBaseView"));

	return TRUE;  // return TRUE unless you set the focus to a control
}

BOOL MyBaseDialogForView::DestroyWindow()
{
	if (m_pMyView) m_pMyView->DestroyWindow();
	return CDialog::DestroyWindow();
}

void MyBaseDialogForView::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(MyBaseDialogForView, CDialog)
	ON_WM_PAINT()
	ON_WM_CLOSE()
	ON_WM_SIZE()
	ON_WM_GETMINMAXINFO()
	ON_WM_MOVE()
	ON_WM_SYSCOMMAND()
	ON_WM_MOUSEWHEEL()
END_MESSAGE_MAP()

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
void MyBaseDialogForView::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);
	if (m_init_done) PlaceElements();
}

void MyBaseDialogForView::OnGetMinMaxInfo(MINMAXINFO* lpMMI)
{
  // set the minimum tracking width
  // and the minimum tracking height of the window
  lpMMI->ptMinTrackSize.x = 40;
  lpMMI->ptMinTrackSize.y = 40;
}


void MyBaseDialogForView::OnMove(int x, int y)
{
	CDialog::OnMove(x, y);
	if (m_init_done) PlaceElements();
}

void MyBaseDialogForView::OnClose()
{	
	CDialog::DestroyWindow();
	CDialog::OnClose();
}

void MyBaseDialogForView::OnPaint()
{
	SetCaller(mystr("REDRAW triggered by MyBaseDialogForView::OnPaint\n"));
	GetView()->UpdateWindow();
}

void MyBaseDialogForView::OnSysCommand(UINT nID, LPARAM lParam) // catches selection in Menu
{
	CDialog::OnSysCommand(nID, lParam);
}

BOOL MyBaseDialogForView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	return GetView()->OnMouseWheel(nFlags, zDelta, pt);
}

BOOL MyBaseDialogForView::PreTranslateMessage(MSG* pMsg)
{
	return CDialog::PreTranslateMessage(pMsg);
}

// ***********************************
// placement of view within the dialog
// ***********************************

void MyBaseDialogForView::PlaceElements()
{
	// embedded view 
	CRect clientrect;
	this->GetClientRect(&clientrect);
	CRect viewrect;
	this->GetViewRect(viewrect);
	GetView()->MoveWindow(&viewrect);

// statusbar
	if (m_flag_use_statusbar)
	{
		CRect statusrect;
		this->GetStatusRect(statusrect);
		m_StatusBar.MoveWindow(&statusrect);
	}
}

CRect MyBaseDialogForView::PlaceView(CRect& viewrect)
{
	CRect oldviewrect;
	this->GetViewRect(oldviewrect);
	ClientToScreen(&oldviewrect);

	CRect oldwindowrect;
	this->GetWindowRect(oldwindowrect);

	CRect newwindowrect;
	newwindowrect.left   = viewrect.left   + oldwindowrect.left   - oldviewrect.left;
	newwindowrect.right  = viewrect.right  + oldwindowrect.right  - oldviewrect.right;
	newwindowrect.top    = viewrect.top    + oldwindowrect.top    - oldviewrect.top;
	newwindowrect.bottom = viewrect.bottom + oldwindowrect.bottom - oldviewrect.bottom;

	this->MoveWindow(newwindowrect.left, newwindowrect.top, newwindowrect.Width(), newwindowrect.Height(),1);
	SetCaller(mystr("REDRAW triggered by MyBaseDialogForView::PlaceView\n"));
	this->UpdateWindow();

	return CRect(oldviewrect);
}


void MyBaseDialogForView::CheckScreenLocation()
{
	LONG sx = (LONG)::GetSystemMetrics(SM_CXFULLSCREEN); 
	LONG sy = (LONG)::GetSystemMetrics(SM_CYFULLSCREEN); 

	CRect r;
	GetWindowRect(&r);

	CPoint p(r.right - 1,r.bottom - 1);

	if (p.x > sx-1) p.x = sx-1;
	if (p.y > sy-1) p.y = sy-1;

	r.top = p.y - r.Height();
	r.left = p.x - r.Width();
	r.right = p.x;
	r.bottom = p.y;
	MoveWindow(r);
	//ShowWindow(SW_SHOW);
}

// **************************
// *** FUNCTIONS FOR PRINTING
// **************************
void MyBaseDialogForView::OnPrint(CDC* pDC, CPrintInfo* pInfo)
{
	if (GetView())
	{
		GetView()->OnPrint(pDC, pInfo);
	}
}
// Print Graph Button
// prints only content of embedded (active) CView
void MyBaseDialogForView::Print()
{
	CDC dc;
	CPrintDialog printDlg(FALSE);
	if (printDlg.DoModal() == IDCANCEL)     // Get printer settings from user
		return;
	dc.Attach(printDlg.GetPrinterDC());     // Get and attach a printer DC
	dc.m_bPrinting = TRUE;

	// DOCINFO
	CString strTitle;                       // Get the application title
	strTitle.LoadString(AFX_IDS_APP_TITLE);
	DOCINFO di;                             // Initialise print document details
	::ZeroMemory (&di, sizeof (DOCINFO));
	di.cbSize = sizeof (DOCINFO);
	di.lpszDocName = strTitle;
	BOOL bPrintingOK = dc.StartDoc(&di);    // Begin a new print job

	// Get the printing extents and store in the m_rectDraw field of a 
	// CPrintInfo object
	CPrintInfo Info;
	Info.m_rectDraw.SetRect(0,0, 
		dc.GetDeviceCaps(HORZRES), 
		dc.GetDeviceCaps(VERTRES));
	OnBeginPrinting(&dc, &Info);            // Call your "Init printing" function

	dc.StartPage();                     // begin new page

	OnPrint(&dc, &Info);                // Call your "Print page" function

	bPrintingOK = (dc.EndPage() > 0);   // end page

	OnEndPrinting(&dc, &Info);            // Call your "Clean up" function

	if (bPrintingOK)
		dc.EndDoc();                        // end a print job
	else
		dc.AbortDoc();                      // abort job.

	dc.DeleteDC();                        // delete the printer DC
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ generic view creation function (not in CFrameWnd)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CWnd* MyBaseDialogForView::CreateNewView(CCreateContext* pContext, CWnd *pParent, CRect& rect, int wID)
{
	CWnd* pWnd = NULL;

	if (pContext != NULL)
	{
		if (pContext->m_pNewViewClass != NULL)
		{
			pWnd = (CWnd*)pContext->m_pNewViewClass->CreateObject();

			if (pWnd == NULL)
			{
				TRACE1("Error: Dynamic create of view %Fs failed\n", pContext->m_pNewViewClass->m_lpszClassName);
				return NULL;
			}
			ASSERT(pWnd->IsKindOf(RUNTIME_CLASS(CWnd)));

			if (!pWnd->Create(NULL, NULL, AFX_WS_DEFAULT_VIEW, rect, pParent, wID, pContext))
			{
				TRACE0("Error: couldn't create view \n");
				return NULL;
			}
			// send initial notification message
			pWnd->SendMessage(WM_INITIALUPDATE);
		}
	}
	return pWnd;
} 