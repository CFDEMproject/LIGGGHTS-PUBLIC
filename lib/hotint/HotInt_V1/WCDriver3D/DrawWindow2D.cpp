//#**************************************************************
//# filename:             BEBeamGetDKappa1.h
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
#include "DrawWindow2D.h"
#include "savewindowbitmap.h"



const int LOGICAL_PER_PIXEL = 10;
const double PIXELS_PER_REAL_UNIT = 8.;
const int DEFAULT_FONT_SIZE = 18;   // should be about PIXELS_PER_REAL_UNIT * 3

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CDrawWindow2DView:                                       derived View class for the 2D - PlotWindow ( Control elements )   * 2012-12-12 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// derived class for CDrawWindow2DView
IMPLEMENT_DYNCREATE(CDrawWindow2DView, MyBaseDialogForView)

CDrawWindow2DView::CDrawWindow2DView(CWnd* parent)
{
	m_pParent = (CDrawWindow2DDlg*) parent;
	Reset();
}

void CDrawWindow2DView::ErrorMessage(mystr& msg)
{
	this->GetParentDlg()->GetWCDI()->GetUserInterface()->AddText( mystr("Illegal pDC - ") + caller );
}

void CDrawWindow2DView::Reset()
{
  pts_y_top = 50;
	pts_y_plot = 2000;
	pts_y_bottom = 50;
	pts_x_left = 50;
	pts_x_plot = 2000;
	pts_x_right = 50;
}

BEGIN_MESSAGE_MAP(CDrawWindow2DView, CMyBase2DView)
	ON_WM_MOUSEWHEEL()
	ON_WM_MOUSEMOVE()
	ON_WM_NCMOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
END_MESSAGE_MAP()

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
BOOL CDrawWindow2DView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	// class specific 
	// no effect for Shift or Ctrl key, always zoom for both directions
  UINT nFlags_specific = nFlags &!MK_SHIFT &!MK_CONTROL;
	
	BOOL rv = CMyBase2DView::OnMouseWheel(nFlags_specific,zDelta,pt); // this already sets Shownrange in CMyBase2DView::ReComputeShownRange_Zooming

	ReComputeZoomFactor(); // compute the new currentzoomfactor
	
	SetCaller(mystr("REDRAW triggered by CDrawWindow2DView::OnMouseWheel\n"));
	
	return RedrawWindow();
}

void CDrawWindow2DView::OnLButtonDown(UINT nFlags, CPoint point)
{
	DragDropStart() = point;
	CView::OnLButtonDown(nFlags, point);
}

void CDrawWindow2DView::OnLButtonUp(UINT nFlags, CPoint point)
{
	AllowRedraw();
	DragDropFinal() = point;

	// catch Object at mouseposition
	int nr_nearest = this->CatchDrawObject(point);
	if(nr_nearest>0)  
	{
		// determine if the nearest element is a IOGUIResponse
		ControlWindowContext_::DrawComponent& elem = m_pParent->GetElement(nr_nearest);

		int elnr = elem.mbs_elnr;
		int is_response_element = GetParentDlg()->pWCDI->CallCompFunction(30, elnr, 0); // 0 to determine elementtype
		if (is_response_element)
		{ 
			this->GetParentDlg()->pWCDI->CallCompFunction(30, elnr, 1, NULL);   // decode "1"  to "single left button (mouse)click"
			GetParentDlg()->ForceUpdate();
		}
	}
	else
	{
	// standard drag drop - drag drop end
		DragDropAction();
		SetCaller(mystr("REDRAW triggered by CMyBase2DView::OnLButtonUp\n"));
		RedrawWindow();
	}
	CView::OnLButtonUp(nFlags, point);
}

void CDrawWindow2DView::OnMouseMove(UINT nFlags, CPoint point)
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

			int dbg_donotredraw = DoNotRedraw();
			if(!DoNotRedraw())                        // CATCH EXCEPTION - wincore.cpp l:248
			{
				SetCaller(mystr("REDRAW triggered by CDrawWindow2DView::OnMouseMove-DragDrop\n"));
				RedrawWindow();
			}
			else
			{
				int caught = 1;
			}
			DragDropStart() = point;
		}
	}

	CPoint capture_LOGICAL = GetLogicalFromPixels(point);

	double x = XofLogical(capture_LOGICAL);
	double y = YofLogical(capture_LOGICAL);		
	SetStatusBarText_X(mystr(" X: ") + mystr(x,4)); // to 4 digits
	SetStatusBarText_Y(mystr(" Y: ") + mystr(y,4)); // to 4 digits

	int nearest_obj = CatchDrawObject(point);
	UpdateStatusBarInfo();
	
	return CView::OnMouseMove(nFlags, point);
}

void CDrawWindow2DView::UpdateStatusBarInfo()
{
	if (selected_object_buffer > 0)
	{
		const ControlWindowContext_::DrawComponent& elem = GetParentDlg()->GetElement(selected_object_buffer);

		if(elem.IsSelectable())
		{
			mystr description = mystr("Selectable Element: ");
			int i = elem.sub_elnr / 65536; // input port number
			//int j = elem.sub_elnr % 65536; // list index nrunningumber
			switch (elem.sub_type)
			{
				case TElementFrame:				description += mystr("Frame of MBS-Element # ") + mystr(elem.sub_elnr) + mystr("."); break;
				case TElementName:				description += mystr("Label of MBS-Element # ") + mystr(elem.sub_elnr) + mystr("."); break;
				case TInputNode:					description += mystr("InputNode #") + mystr(elem.sub_elnr) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				case TOutputNode:					description += mystr("OutputNode #") + mystr(elem.sub_elnr) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				case TConstructionNode:		description += mystr("ConstructionNode for input # ") + mystr(i) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				case TConnectionLine:			description += mystr("ConnectionLine for input # ") + mystr(i) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				case TTextObject:					description += mystr("TextObject #") + mystr(elem.sub_elnr) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				case TSymbol:							description += mystr("ElementSymbol #") + mystr(elem.sub_elnr) + mystr(" (of MBS-Element # ") + mystr(elem.mbs_elnr) + mystr(")"); break;
				default:									description += mystr("NONE");
			}
			SetStatusBarText_Main(description);
		}
		else 
			SetStatusBarText_Main(mystr("")); // no selectable object
	}
	else
	{
		SetStatusBarText_Main(mystr("no selectable element"));
	}
}

void CDrawWindow2DView::OnNcMouseMove(UINT nHitTest, CPoint point)
{
	// this is reached when mouse moves OFF the client rect
	SetStatusBarText_X(mystr(""));
	SetStatusBarText_Y(mystr(""));
	CMyBase2DView::OnNcMouseMove(nHitTest, point);
}

void CDrawWindow2DView::DragDropAction()
{
	if (selected_object_buffer>0)
	{
		const ControlWindowContext_::DrawComponent& elem = GetParentDlg()->GetElement(selected_object_buffer);
		ControlWindowContext_::DrawComponent& elem_nc = GetParentDlg() -> GetElement(selected_object_buffer); 

		if (elem.IsSelectableMBSElement() || elem.IsSelectableConstructionNode() )
		{
// buffered element is MOVABLE independently - move only this element

			CPoint logical_start = GetLogicalFromPixels(DragDropStart(),mystr("DragDropAction - StartPoint"));
			CPoint logical_final = GetLogicalFromPixels(DragDropFinal(),mystr("DragDropAction - FinalPoint"));

			double x_start = XofLogical(logical_start);
			double y_start = YofLogical(logical_start);
			double x_final = XofLogical(logical_final);
			double y_final = YofLogical(logical_final);

			double x_shift = x_final-x_start;
			double y_shift = y_final-y_start;
			int remember_selected_object_buffer = selected_object_buffer;

			if (elem.IsSelectableMBSElement()) // entire Element is movable
			{
				int mbs_elem_nr = elem.mbs_elnr;
				m_pParent -> GetWCDDlg() -> GetWCDInterface() -> MoveElement(mbs_elem_nr, x_shift, y_shift, 0.);       // apply change of position
				m_pParent -> ForceUpdate();                                                                            // redraw wich updated position
				this->selected_object_buffer = remember_selected_object_buffer;                                        // keep "selected object", has been updated during redraw
			}

			if (elem.IsSelectableConstructionNode()) // construction node is movable
			{
				int mbs_elem_nr = elem.mbs_elnr;
				int list_idx = elem.sub_elnr % 65536;                 // last 16 bits of subelemnr code the array position to insert
				m_pParent -> GetWCDDlg() -> GetWCDInterface() -> MoveConNode2D(mbs_elem_nr, list_idx, x_shift, y_shift);
				
				//BOOL LeftMouseButtonDown = ( (GetKeyState(VK_LBUTTON) & 0x80) != 0 );
				//if (LeftMouseButtonDown)
				//{
				//	elem_nc.Center() += Vector2D(x_shift,y_shift);
				//	m_pParent->RedrawWindow();
				//}
				//else
				//{
						m_pParent -> ForceUpdate();                                                                            // redraw wich updated position
						this->selected_object_buffer = remember_selected_object_buffer;                                        // keep "selected object", has been updated during redraw
				//}
			}
		//	return CMyBase2DView::DragDropAction();	// under construction: replace this with moveobject2d ....

		}
		else
		{
// buffered element is not movable independently - move entire sheet
			return CMyBase2DView::DragDropAction();	
		}
	}
	else
	{
// no element selected - move entire sheet
		return CMyBase2DView::DragDropAction();
	}		
}

// *******************
// THE DRAWING PROCESS
// *******************
void CDrawWindow2DView::OnDraw(CDC* pDC)
{
	// TRY TO CATCH MANY REDRAW MESSAGES
		// delay redraw ??
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_DRAWITEM, WM_DRAWITEM, PM_REMOVE))
	{
		int dbg=1;
	}

// is redraw currently allowed ?
// Drawing could be prevented by: ForceUpdate, OnDraw
	if(DoNotRedraw()) 
	{
		return;
	}
	ForbidRedraw();

 	_is_executing_ondraw = 0;

// reset caller string for any calls by system
	SetCaller(mystr("REDRAW triggered by SYSTEM\n"));

	CRect plotrect;
	this->GetPlotRect(plotrect);

	if(shownrange.Empty())
	{
		ComputeFullRange(1);
		shownrange = fullrange;
	}

	if(0)
	{
		DrawElements(pDC);
	}
	else
	{
		CMyMemDC MDC(pDC);
		DrawElements(&MDC);
	}
	
	_is_executing_ondraw = 0;
	AllowRedraw();
}

void CDrawWindow2DView::DrawElements(CDC *pDC)
{
	this->SelectDefaultFont(pDC, 20);
	
	if(0) // HACK draw "+" at origin
	{
		pDC->MoveTo( ValuesToLogical(-.1,0.));
		pDC->LineTo( ValuesToLogical( .1,0.));
		pDC->MoveTo( ValuesToLogical(0.,-.1));
		pDC->LineTo( ValuesToLogical(0., .1));
	}

// DRAW the elements
	int n = GetParentDlg()->NElements();
	for (int i=1; i<=n; i++)
	{
		_is_executing_ondraw = n;
		ControlWindowContext_::DrawComponent& elem = m_pParent->GetElement(i);
		if (elem.IsLine()) { DrawLine(pDC,elem); }
		else if (elem.IsRectangle())	{	DrawRectangle(pDC,elem); }
		else if (elem.IsEllipse()) { DrawEllipse(pDC,elem); }
		else if (elem.IsText()) { DrawText(pDC,elem); }
	}
}

BOOL CDrawWindow2DView::IsVisible(CDC* pDC, ControlWindowContext_::DrawComponent& elem)
{
	// TODO: VISIBILITY CHECK
	return true;
}

// Element based on a Line
BOOL CDrawWindow2DView::DrawLine(CDC* pDC, ControlWindowContext_::DrawComponent& elem)
{
// draws a line from Point1 (elem::center) to Point2 (elem::size)
	// Pen for Line (color,thickness)
	COLORREF col = (COLORREF) ((int)(elem.col(1)*255+0.5) + ((int)(elem.col(2)*255+0.5)<<8) + ((int)(elem.col(3)*255+0.5)<<16));
	int thickness = ( elem.IsSymbol() ? LOGICAL_PER_PIXEL*50 : 1);								// thin line for connection lines, thicker for symblos
	CPen pen (PS_SOLID, thickness, col);
	CPen* prevPen = pDC->SelectObject(&pen);
	
	//draw
	CPoint startpoint = ValuesToLogical(elem.P1().X(), elem.P1().Y());
	CPoint endpoint = ValuesToLogical(elem.P2().X(), elem.P2().Y());
	pDC->MoveTo(startpoint);
	pDC->LineTo(endpoint);

	// previous pen
	pDC->SelectObject(prevPen);
	return true;
}

// Element based on a Rectangular ( Frame of all elements )
BOOL CDrawWindow2DView::DrawRectangle(CDC* pDC, ControlWindowContext_::DrawComponent& elem)
{
// draw a rectangle 
	// Pen for Border (color,thickness)
	COLORREF col = (COLORREF) ((int)(elem.col(1)*255+0.5) + ((int)(elem.col(2)*255+0.5)<<8) + ((int)(elem.col(3)*255+0.5)<<16));
	int thickness = ( elem.IsSymbol() ? LOGICAL_PER_PIXEL*50 : 1);								// thin line for connection lines, thicker for symblos
	CPen pen(PS_SOLID, thickness, col);
	CPen* prevPen = pDC->SelectObject(&pen);
	
	// color for area
	COLORREF col2;
	if ( (elem.col2(1) < 0. || elem.col2(2) < 0. || elem.col2(3) < 0.) ) col2 = pDC->GetBkColor();
	else col2 = (COLORREF) ((int)(elem.col2(1)*255+0.5) + ((int)(elem.col2(2)*255+0.5)<<8) + ((int)(elem.col2(3)*255+0.5)<<16)); 
	CBrush brush(col2);
	CBrush* prevBrush = pDC->SelectObject(&brush);

	// draw
	CPoint p1 = ValuesToLogical(elem.GetXMin(), elem.GetYMin());
	CPoint p2 = ValuesToLogical(elem.GetXMax(), elem.GetYMax());
	CRect rect(p1,p2);
	pDC->Rectangle(rect);

	// previous pen and brush
	pDC->SelectObject(prevPen);
	pDC->SelectObject(prevBrush);
	return true;
}

// Element based on an Ellipse ( e.g. Nodes )
BOOL CDrawWindow2DView::DrawEllipse(CDC* pDC, ControlWindowContext_::DrawComponent& elem)
{ 
// draw an ellipse
	// Pen for Border (color,thickness)
	COLORREF col = (COLORREF) ((int)(elem.col(1)*255+0.5) + ((int)(elem.col(2)*255+0.5)<<8) + ((int)(elem.col(3)*255+0.5)<<16));
	int thickness = ( elem.IsSymbol() ? LOGICAL_PER_PIXEL*50 : 1);								// thin line for connection lines, thicker for symblos
	CPen pen(PS_SOLID, thickness, col);
	CPen* prevPen = pDC->SelectObject(&pen);
	
	// color for area
	COLORREF col2;
	if ( (elem.col2(1) < 0. || elem.col2(2) < 0. || elem.col2(3) < 0.) ) col2 = pDC->GetBkColor();
	else col2 = (COLORREF) ((int)(elem.col2(1)*255+0.5) + ((int)(elem.col2(2)*255+0.5)<<8) + ((int)(elem.col2(3)*255+0.5)<<16)); 
	CBrush brush(col2);
	CBrush* prevBrush = pDC->SelectObject(&brush);

	//draw
	CPoint p1 = ValuesToLogical(elem.GetXMin(), elem.GetYMin());
	CPoint p2 = ValuesToLogical(elem.GetXMax(), elem.GetYMax());
	CRect rect(p1,p2);
	pDC->Ellipse(rect);

	// previous pen and brush
	pDC->SelectObject(prevPen);
	pDC->SelectObject(prevBrush);
	return true;
}

// Textbox
BOOL CDrawWindow2DView::DrawText(CDC* pDC, ControlWindowContext_::DrawComponent& elem)
{
	CPoint p1 = ValuesToLogical(elem.GetXMin(), elem.GetYMin());
	CPoint p2 = ValuesToLogical(elem.GetXMax(), elem.GetYMax());
	CRect rect(p1,p2);

	// CATCH ir rect is too smal
	if(abs(rect.Height()) < 15)
	{
			return false;
	}
	
	int fontsize= (int) (DEFAULT_FONT_SIZE * LOGICAL_PER_PIXEL * currentzoomfactor + 0.5);

	if (elem.sub_type == TElementName) // element name beneath the Element
	{
		fontsize = (int) (DEFAULT_FONT_SIZE * LOGICAL_PER_PIXEL * currentzoomfactor + 0.5);
	}
	if (elem.sub_type == TSymbol) // text in center of element
	{
		fontsize = (int) (DEFAULT_FONT_SIZE * LOGICAL_PER_PIXEL * currentzoomfactor * 1.3 + 0.5);
	}
  if (elem.sub_type == TTextObject)
	{
		fontsize = (int) (DEFAULT_FONT_SIZE * LOGICAL_PER_PIXEL * currentzoomfactor + 0.5);
	}

// check if font size is too small
	if (fontsize == 0) fontsize = 1;

	CFont* font = GetDefaultFont(pDC,fontsize);
	CFont* prevFont = pDC->SelectObject(font);

	// color for Text
	COLORREF col = (COLORREF) ((int)(elem.col(1)*255+0.5) + ((int)(elem.col(2)*255+0.5)<<8) + ((int)(elem.col(3)*255+0.5)<<16));
	COLORREF prevcolor = pDC->SetTextColor(col);

// Additional Debug output: Fontsize and textextent 
	CSize textsize = pDC->GetTextExtent(CString(elem.text));
	this->SetStatusBarText_Main(mystr("BASE:") + mystr(DEFAULT_FONT_SIZE * LOGICAL_PER_PIXEL) + mystr(" FS:") + mystr(fontsize) + mystr(" ZF:") + mystr(currentzoomfactor));
	this->SetStatusBarText_X(mystr(textsize.cx));
	this->SetStatusBarText_Y(mystr(textsize.cy));
// Additional Debug output:

	this->DrawTextInRect(pDC, CString(elem.text), rect, 0., elem.TextAllign());
//	this->DrawTextAt(pDC, CString(elem.text), rect.CenterPoint(), CPoint( (int)(-textsize.cx*0.5), (int)(-textsize.cy*0.5)), 0., elem.TextAllign());

	pDC->SelectObject(prevFont);
	pDC->SetTextColor(prevcolor);
	return true;
}

int CDrawWindow2DView::DoNotRedraw()
{ 
	return GetParentDlg()->DoNotRedraw(); 
}
void CDrawWindow2DView::ForbidRedraw() 
{ 
	GetParentDlg()->ForbidRedraw(); 
}
void CDrawWindow2DView::AllowRedraw() 
{ 
	GetParentDlg()->AllowRedraw(); 
}

void CDrawWindow2DView::SetStatusBarText_Main (mystr text) // slot1
{ 
	if (!GetParentDlg()->m_flag_statusinfo) text = "";
	(GetParentDlg()->m_StatusBar).SetText(text,1,0); 
} 
void CDrawWindow2DView::SetStatusBarText_X (mystr text) { (GetParentDlg()->m_StatusBar).SetText(text,2,0); } // slot2
void CDrawWindow2DView::SetStatusBarText_Y (mystr text) { (GetParentDlg()->m_StatusBar).SetText(text,3,0); } // slot3


double CDrawWindow2DView::ReComputeZoomFactor() 
{
	currentzoomfactor = initshownrange.SizeX() / shownrange.SizeX();
	if (currentzoomfactor<0 || currentzoomfactor>1e8)
		return currentzoomfactor;
	return currentzoomfactor;
}

void CDrawWindow2DView::SetInitialRange()
{
	// real coordinate boundaries of scene
	ComputeFullRange();

	// available pixels
	CRect curr; GetClientRect(curr);
	long currHeight = abs(curr.Height());
	long currWidth = abs(curr.Width());

	// initial displayed range - chosen such that default IOElement has agreeable size...
	double xmin = fullrange.Center().X() - (double)currWidth / PIXELS_PER_REAL_UNIT;
	double xmax = fullrange.Center().X() + (double)currWidth / PIXELS_PER_REAL_UNIT;
	double ymin = fullrange.Center().Y() - (double)currHeight / PIXELS_PER_REAL_UNIT;
	double ymax = fullrange.Center().Y() + (double)currHeight / PIXELS_PER_REAL_UNIT;
	SetShownRange(xmin, xmax, ymin, ymax);
	initshownrange = Box2D(Vector2D(xmin, ymin), Vector2D(xmax, ymax));
	currentzoomfactor = 1.;
}

// adjust the available logical points to the available pixels... 
void CDrawWindow2DView::MatchLogicalPoints()
{
	// available pixels
	CRect curr; GetClientRect(curr);
	long currHeight = abs(curr.Height());
	long currWidth = abs(curr.Width());

	// set the range for logical points
	int BM = 20; // base multiplicator -> 20 means that origin is at 5%
	this->pts_x_left =		0;																		// currWidth * LOGICAL_PER_PIXEL * 1;
	this->pts_x_plot =		currWidth * LOGICAL_PER_PIXEL * BM;		// currWidth * LOGICAL_PER_PIXEL * (BM-2);
	this->pts_x_right =		0;																		// currWidth * LOGICAL_PER_PIXEL * 1;
	this->pts_y_top =			0;																		// currHeight * LOGICAL_PER_PIXEL * 1;
	this->pts_y_plot =		currHeight * LOGICAL_PER_PIXEL * BM;	// currHeight * LOGICAL_PER_PIXEL * (BM-2);
	this->pts_y_bottom =	0;																		// currHeight * LOGICAL_PER_PIXEL * 1;
}

void CDrawWindow2DView::Rescale()
{
	CRect prev = prevClientRect;
	long prevHeight = abs(prev.Height());
	long prevWidth = abs(prev.Width());
	
	CRect curr; GetClientRect(curr);
	long currHeight = abs(curr.Height());
	long currWidth = abs(curr.Width());

	// scale the shown range to fit 
	// f.t.t.b. assume right/lower boundary is moved ( xmax, ymax )
	// later could compare WindowRect instead of ClientRect to find out which border was moved 
	double stretchfactor = (double)abs(curr.Width()) / (double)abs(prev.Width());

	Box2D prevshownrange = shownrange;
 	double near_x = shownrange.PMin().X();      // assumed to be same
	double far_x = shownrange.PMin().X() + shownrange.SizeX()*(double)abs(curr.Width()) / (double)abs(prev.Width());
	double near_y = shownrange.PMax().Y() - shownrange.SizeY()*(double)abs(curr.Height()) / (double)abs(prev.Height());
	double far_y = shownrange.PMax().Y();
	Box2D newshownrange(Vector2D(near_x,near_y),Vector2D(far_x,far_y));            
	SetShownRange(near_x, far_x, near_y, far_y);

	Box2D previnitshownrange = initshownrange;
 	near_x = initshownrange.PMin().X();      // assumed to be same
	far_x = initshownrange.PMin().X() + initshownrange.SizeX()*(double)abs(curr.Width()) / (double)abs(prev.Width());
	near_y = initshownrange.PMax().Y() - initshownrange.SizeY()*(double)abs(curr.Height()) / (double)abs(prev.Height());
	far_y = initshownrange.PMax().Y();
	Box2D newinitshownrange(Vector2D(near_x,near_y),Vector2D(far_x,far_y));            
	initshownrange = newinitshownrange;

	MatchLogicalPoints();
}

// computes limits of datasets - real values
void CDrawWindow2DView::ComputeFullRange(int flag_equal)
{
	double xmin,xmax,ymin,ymax;
	int n = GetParentDlg()->NElements();
	if (n==0) 
	{ 
		//fullrange.Clear(); 
		fullrange=Box2D(Vector2D(-1,-1),Vector2D(1,1));
	}
	else 
	{
		ControlWindowContext_::DrawComponent& elem1 = GetParentDlg()->GetElement(1);
		xmin = elem1.GetXMin();
		xmax = elem1.GetXMax();
		ymin = elem1.GetYMin();
		ymax = elem1.GetYMax();
		for (int i=2; i<=n; i++)
		{
			ControlWindowContext_::DrawComponent& elem = GetParentDlg()->GetElement(i);
			if (elem.GetXMin()<xmin) xmin = elem.GetXMin();
			if (xmax<elem.GetXMax()) xmax = elem.GetXMax();
			if (elem.GetYMin()<ymin) ymin = elem.GetYMin();
			if (ymax<elem.GetYMax()) ymax = elem.GetYMax();
		}
    fullrange = Box2D( Vector2D(xmin,ymin), Vector2D(xmax,ymax) );

		// make x and y axis intervals same size 
		if(flag_equal)
		{
			double ylog_by_xlog = (double) Get_YRange_logical() / (double) Get_XRange_logical();
			double yran_by_xran = fullrange.SizeY() / fullrange.SizeX();
			{
				if(ylog_by_xlog < yran_by_xran)
				{
					// y-axis is fits, adjust x-axis
					double middle = fullrange.Center().X();
					double newrange = fullrange.SizeX() * yran_by_xran / ylog_by_xlog;
					double newmin = middle - newrange * 0.5;
					double newmax = middle + newrange * 0.5;

					fullrange = Box2D(Vector2D(newmin,fullrange.PMin().Y()), Vector2D(newmax,fullrange.PMax().Y()));
				}
				else
				{
					// x-axis is fits, adjust y-axis
					double middle = fullrange.Center().Y();
					double newrange = fullrange.SizeY() / yran_by_xran * ylog_by_xlog;
					double newmin = middle - newrange * 0.5; 
					double newmax = middle + newrange * 0.5;

					fullrange = Box2D(Vector2D(fullrange.PMin().X(),newmin), Vector2D(fullrange.PMax().X(),newmax));
				}
			}
		}
	}
}

void CDrawWindow2DView::SetShownRange(double xmin, double xmax, double ymin, double ymax)
{
	shownrange = Box2D( Vector2D(xmin,ymin), Vector2D(xmax,ymax) );
}

// new code to catch the objects
int CDrawWindow2DView::CatchDrawObject( CPoint capture_PIXELS )
{
	//CPoint capture_LOGICAL = GetLogicalFromPixels(capture_PIXELS, mystr("CatchDrawObject"));
	//BOOL LeftMouseButtonDown = ( (GetKeyState(VK_LBUTTON) & 0x80) != 0 );

	selected_object_buffer = 0;
  int priority_buffer = 0;
	int priority = 0;                  // option for multiple valid objects, pick the one with highest priority
	int distance_buffer = 0xFFFFul;           // option for multiple valid objects, pick the nearest ( priority is checked first )
  int distance = 0xFFFFul;

// define a cutoff distance for each element type - distances that compute to more than this value cannot give focus
	const int cutoff_response = 0;   
	const int cutoff_nodes = 30;
	const int cutoff_lines = 20;
	const int cutoff_rects = 0;

// define a priority order - rects over nodes over lines...
  const int priority_response = 5;
	const int priority_nodes = 3;
	const int priority_lines = 1;
	const int priority_rects = 5;

	int n = GetParentDlg()->NElements();
	if(n==0)
	{
		return 0;                        // empty list
	}
	else
	{
//AD: only one loop, priority order as in the array ... 4 
		for(int i = 1; i<=n; i++)
		{
			ControlWindowContext_::DrawComponent& elem = GetParentDlg()->GetElement(i);
			int distance =  GetDistanceFromComponent(i,capture_PIXELS);

// only consider elements within the defined cutoff distance
			if (elem.IsSymbol())
			{
				priority = priority_response;
				if (distance>cutoff_response) 
					distance = 0xFFFFul;
			}
			else if (elem.IsEllipse())
			{
				priority = priority_nodes;
				if (distance>cutoff_nodes)
					distance = 0xFFFFul;
			}
			else if (elem.IsLine())
			{
				priority = priority_lines;
				if (distance>cutoff_lines)
					distance = 0xFFFFul;
			}
			else if (elem.IsRectangle())
			{
				priority = priority_response;
				if (distance>cutoff_rects)
					distance = 0xFFFFul;
			}
			else
			{
				priority = 0;
				distance = 0xFFFFul;
			}
// new selected object when Priority is same and distance is smaller OR priority is higher and distance is a valid number
			if( (priority==priority_buffer && distance<distance_buffer) || 
				  (priority>priority_buffer && distance<0xFFFFul) )
			{
				selected_object_buffer = i;
				priority_buffer = priority;
				distance_buffer = distance;
			}
		}
	return selected_object_buffer;
	}
}

int CDrawWindow2DView::GetDistanceFromComponent(int i, CPoint capture_PIXELS)
{
	//CPoint capture_LOGICAL = GetLogicalFromPixels(capture_PIXELS, mystr("CatchDrawObject"));

// checks distance IN PIXELS to the object ( line / rectangle / TODO: other shapes)
	ControlWindowContext_::DrawComponent& elem = GetParentDlg()->GetElement(i);
	int dist = 0xFFFFul;

// use the elements subtype to decide which shape must be checked
	if (elem.IsLine())
	{
		// catch line
		CPoint logical;
		logical = ValuesToLogical(elem.center.X(), elem.center.Y()); 
		CPoint p1 = GetPixelsFromLogical(logical);
		logical = ValuesToLogical(elem.size.X(), elem.size.Y());      		
		CPoint p2 = GetPixelsFromLogical(logical);
		dist = DistFromLine(p1,p2,capture_PIXELS); 
	}
	else if (elem.IsRectangle() || elem.IsEllipse() ) 
	{
		// catch rectangular
		CPoint logical;
		logical = ValuesToLogical(elem.GetXMin(), elem.GetYMin());				
		CPoint p1 = GetPixelsFromLogical(logical);		
		logical = ValuesToLogical(elem.GetXMax(), elem.GetYMin());				
		CPoint p2 = GetPixelsFromLogical(logical);
		logical = ValuesToLogical(elem.GetXMax(), elem.GetYMax());				
		CPoint p3 = GetPixelsFromLogical(logical);
		logical = ValuesToLogical(elem.GetXMin(), elem.GetYMax());				
		CPoint p4 = GetPixelsFromLogical(logical);

		if( ( Minimum(p1.x,p3.x) <= capture_PIXELS.x) && ( capture_PIXELS.x <= Maximum(p1.x,p3.x) ) &&
        ( Minimum(p1.y,p3.y) <= capture_PIXELS.y) && ( capture_PIXELS.y <= Maximum(p1.y,p3.y) ) )
		{
		// is IN the rectangle - some suitable constant value that allows to select nested elements above the main rectangle ( ?>0 )
			dist =  0; 
		}
		else
		{
			dist = Minimum( Minimum(DistFromLine(p1,p2,capture_PIXELS), DistFromLine(p2,p3,capture_PIXELS)),
				              Minimum(DistFromLine(p3,p4,capture_PIXELS), DistFromLine(p4,p1,capture_PIXELS)) );
		}
	}
	else
	{
		;
	}
	return dist;
}

int CDrawWindow2DView::DistFromLine(CPoint p1, CPoint p2, CPoint pmouse)
{
		Vector2D vp1(p1.x,p1.y);
		Vector2D vp2(p2.x,p2.y);
		Vector2D vpm(pmouse.x, pmouse.y);
	//	double dist2 = MinDistLP(vp1,vp2,vpm);
		double dist2;
		Vector2D v = vp2-vp1;
		Vector2D vlp = vpm-vp1;
		double num = v*vlp;
		double den = v*v;
		if(num<=0.) dist2 = (vpm-vp1).Norm();
		else if(num>=den) dist2 = (vpm-vp2).Norm();
		else if(den>0.) dist2 = sqrt(vlp*vlp - num * num /den);
		else dist2 = vlp.Norm();

		return int(dist2+0.5);
}

#define ID_VIEW_ON_POPUP 1234
#define ID_MENU_PARENT_ON_POPUP 1235
#define ID_STATUSBAR_ON_POPUP 1236

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CDrawWindow2DDlg:                               this is the dialog that nests the Graph (CMyView object),                  * 2012-12-12 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// DrawWindow2DDlg-Dialogfeld
IMPLEMENT_DYNAMIC(CDrawWindow2DDlg, MyBaseDialogForView)

CDrawWindow2DDlg::CDrawWindow2DDlg(CWnd* pParent /*=NULL*/)
	: MyBaseDialogForView(pParent)
{
	m_init_done = 0;
}

CDrawWindow2DDlg::~CDrawWindow2DDlg()
{
}

BOOL CDrawWindow2DDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// create the View
	CCreateContext cc;
	cc.m_pNewViewClass = RUNTIME_CLASS(CDrawWindow2DView);
	cc.m_pCurrentDoc = NULL;
	cc.m_pNewDocTemplate = NULL;
	cc.m_pLastView = NULL;
	cc.m_pCurrentFrame = NULL;

	m_pMyView = (CDrawWindow2DView*)CreateNewView(&cc, this, CRect(0, 0, 0, 0), ID_VIEW_ON_POPUP);
	if (m_pMyView == NULL)
		EndDialog(IDCANCEL);

	CRect rect;
	this->GetWindowRect(&rect);
	rect.top = rect.bottom - 25;
	int dbg = this->m_StatusBar.Create(WS_CHILD | WS_BORDER | WS_VISIBLE, rect, this, ID_STATUSBAR_ON_POPUP);

	// place the View in the client rectangle
	this->GetViewRect(GetView()->prevClientRect);              // compute the current client rectangle of the nesting dialog
	GetView()->prevClientRect.bottom -= 4;
	GetView()->prevClientRect.right -=4;

	PlaceElements();                                           
	GetView()->SetParentDlg(this);
	m_init_done = 1;

	// initial scene
	pWCDI->RenderControlWindow(this);													// Data is transfered to IO Blocks window
	GetView()->SetInitialRange();                             // initial range ( real coordinates ) is computed
	GetView()->MatchLogicalPoints();                          // logical points for the view are matched to the extend (pixels)
	// 

	SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::OnInitDialog\n"));
	this->UpdateWindow();
	this->SetDisplayName(CString("IO Blocks"));
	AllowRedraw();

	return TRUE;  // return TRUE unless you set the focus to a control
}

BOOL CDrawWindow2DDlg::DestroyWindow()
{
	if (GetView())
		GetView()->DestroyWindow();
	return CDialog::DestroyWindow();
}

void CDrawWindow2DDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CDrawWindow2DDlg, MyBaseDialogForView)
	ON_WM_PAINT()
	ON_WM_CLOSE()
	ON_WM_SIZE()
	ON_WM_GETMINMAXINFO()
	ON_WM_MOVE()
	ON_WM_SYSCOMMAND()
	ON_WM_MOUSEWHEEL()
	// CONTEXTMENU	
	ON_WM_CONTEXTMENU()
	ON_COMMAND(ID_CM_ZOOMTOFULLRANGE,		&CDrawWindow2DDlg::OnCMZoomToFullRange)
	ON_COMMAND(ID_CM_EXPORTTOFILE,			&CDrawWindow2DDlg::OnCMExportToFile)
	ON_COMMAND(ID_CM_PRINTGRAPH,				&CDrawWindow2DDlg::OnCMPrintGraph)
	ON_COMMAND(ID_CM_INFOBLOCKVIEW,			&CDrawWindow2DDlg::OnCMShowStatusInfo)
	ON_COMMAND(ID_CM_SELECTOBJECT,			&CDrawWindow2DDlg::OnCMSelectobject)
END_MESSAGE_MAP()

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
void CDrawWindow2DDlg::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);

	if (m_init_done)
	{
		PlaceElements();                              // all rescaling (zoomfactors, shoenrange, etc.) is done 
		// placeelements calls Redraw()
		CRect prevClientRect;	
		GetViewRect(prevClientRect);                  // get the current size ( before resize ) IN PIXELS   
		prevClientRect.bottom -=4;
		prevClientRect.right -=4;
		GetView()->prevClientRect = prevClientRect;   // remember previous client
	}
}

void CDrawWindow2DDlg::OnGetMinMaxInfo(MINMAXINFO* lpMMI)
{
  // set the minimum tracking width
  // and the minimum tracking height of the window
  lpMMI->ptMinTrackSize.x = 150;
  lpMMI->ptMinTrackSize.y = 100;
}


void CDrawWindow2DDlg::OnMove(int x, int y)
{
	CDialog::OnMove(x, y);
	//if (m_init_done) PlaceElements();
}

void CDrawWindow2DDlg::OnClose()
{	
	this->ShowWindow(SW_HIDE);
}

void CDrawWindow2DDlg::OnPaint()
{
	SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::OnPaint\n"));
	GetView()->UpdateWindow();
}

void CDrawWindow2DDlg::OnSysCommand(UINT nID, LPARAM lParam) // catches selection in Menu
{
	CDialog::OnSysCommand(nID, lParam);
}

BOOL CDrawWindow2DDlg::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	return GetView()->OnMouseWheel(nFlags, zDelta, pt);
}

BOOL CDrawWindow2DDlg::PreTranslateMessage(MSG* pMsg)
{
	if(pMsg->message==WM_KEYDOWN)
	{
		int key = pMsg->wParam;                                   // get the key-value
		int nr_of_changes = pWCDI->CallCompFunction(31, key);			// call the MBS-function
		if (nr_of_changes>0) ForceUpdate();                       // update if an MBS-element was changed 
	}
	return MyBaseDialogForView::PreTranslateMessage(pMsg);
}

// force a full update ( call DrawSystem2D )
void CDrawWindow2DDlg::ForceUpdate()
{
	// LOCK CALL OF DRAWING ROUTINE WHILE ARRAYS ARE UPDATED ( ACCESS TO ARRAYS NOT SAFE )
	ForbidRedraw();                               // manually block drawing of the dialog while the system is updated
	pWCDI->RenderControlWindow(this);
	// RELASE LOCK WHEN UPDATE IS DONE !
	AllowRedraw();
	if (m_init_done)
	{
		SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::ForceUpdate\n"));
		RedrawWindow();
	}
}

// ***************************************************************************************
// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
// ***************************************************************************************
void CDrawWindow2DDlg::PlaceElements()
{
  // embed view on dialog
	CRect viewrect;
	this->GetViewRect(viewrect);              // compute the current view rectangle 
	GetView()->MoveWindow(&viewrect);         // does move window trigger a call of OnDraw here ? --> NO
	
	// statusbar
	CRect statusrect;
	this->GetStatusRect(statusrect);
	m_StatusBar.MoveWindow(&statusrect);

// two 70 point slots for coordinates at right end with 30 points additional border
	int totalwidth = statusrect.Width();
	int m_widths[4] = {0, totalwidth-170, totalwidth-100, totalwidth-30};
	m_StatusBar.SetMinHeight(27);
	m_StatusBar.SetParts(4, m_widths);

	if(m_init_done)
		GetView()->Rescale();                 // function that rescales the logical extents and shownrange such that 
	

	RedrawWindow(0,0,RDW_FRAME|RDW_INVALIDATE);

//	if(!DoNotRedraw())                        // CATCH EXCEPTION - wincore.cpp l:248
//	{
//// this line can be reached when OnDraw is in progress...
//		//SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::PlaceElement\n"));
//		//this->RedrawWindow();                   // redraw Dialog
//		SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::PlaceElement II\n"));
//		GetView()->RedrawWindow();              // redraw View            
//	}
//	else
//	{
//		int caught = 1;
//	}
}

// ***************************************************************************************
// context menu -
// ***************************************************************************************
void CDrawWindow2DDlg::OnContextMenu(CWnd* pWnd, CPoint point)
{
	CMenu contextmenu;

	contextmenu.CreatePopupMenu();

// when an element can be selected
	int nr_nearest_drawelem = this->GetView()->selected_object_buffer;
	int nr_nearest_mbselem = 0;
	mystr str;
// identify the corresponding element number ( and type )
	if(nr_nearest_drawelem>0)
	{
		ControlWindowContext_::DrawComponent& elem = GetElement(nr_nearest_drawelem);
		if(elem.IsSelectableMBSElement())
		{
// specific code for a selected element frame --> display "Edit Element X"
			str = mystr("Edit Element ") + mystr(elem.mbs_elnr);
			contextmenu.AppendMenu(MF_STRING, ID_CM_SELECTOBJECT,	str.c_str());
		}
		else if(elem.IsSelectableConnectionLine())
		{
			str = mystr("Insert Construction Node") + mystr(elem.sub_elnr);
			contextmenu.AppendMenu(MF_STRING, ID_CM_SELECTOBJECT,	str.c_str());
		}
		else if(elem.IsSelectableConstructionNode())
		{
			str = mystr("Delete Construction Node") + mystr(elem.sub_elnr);
			contextmenu.AppendMenu(MF_STRING, ID_CM_SELECTOBJECT,	str.c_str());
		}
		contextmenu.AppendMenu(MF_SEPARATOR);
	}

	contextmenu.AppendMenu(MF_STRING, ID_CM_ZOOMTOFULLRANGE, "Default Range");
	contextmenu.AppendMenu(MF_STRING, ID_CM_EXPORTTOFILE,	"Save Screen");
	contextmenu.AppendMenu(MF_STRING, ID_CM_PRINTGRAPH,	"Print Screen");
  
	int nFlags = (MF_STRING | ( m_flag_statusinfo ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_INFOBLOCKVIEW,				"Show Status Bar Information"); //(V) or (-)

	contextmenu.TrackPopupMenu(TPM_LEFTALIGN| TPM_RIGHTBUTTON,point.x,point.y,this,NULL);
}

void CDrawWindow2DDlg::OnCMSelectobject()
{
	int nr_nearest_drawelem = this->GetView()->selected_object_buffer;
	int nr_nearest_mbselem = 0;
// identify the corresponding element number ( and type )
	if(nr_nearest_drawelem>0)
	{
		ControlWindowContext_::DrawComponent& elem = GetElement(nr_nearest_drawelem);
		if(elem.IsSelectable())
		{
			if(elem.IsSelectableMBSElement())
			{
// specific code for a selected element frame --> open edit element dialog
				int mbs_elem_nr = elem.mbs_elnr;

				this->ShowWindow(SW_HIDE); //$AD 2013-07-01: WORKAROUND make sure that EditElementProperties Dialog is visible   
				GetWCDDlg()->EditElementProperties(mbs_elem_nr, /*?*/1/*?*/);				
				this->ForceUpdate();	//$ AD 2013-07-01: refresh 2D window after Prpoerties are changed 
				this->ShowWindow(SW_SHOW); //$AD 2013-07-01: WORKAROUND make sure that EditElementProperties Dialog is visible  
			}
			else if(elem.IsSelectableConnectionLine())
			{
	// specific code for a selected connection line --> new node in the middle of the line
				int mbs_elem_nr = elem.mbs_elnr;
				int list_idx = elem.sub_elnr % 65536;                 // last 16 bits of subelemnr code the array position to insert
				int inputnr = elem.sub_elnr / 65536;                  // first 16 bits of subelemnr code the associated input number
				Vector2D pos = 0.5 *( elem.P1()+elem.P2() );					// position in the middle of the line

				GetWCDDlg()->GetWCDInterface()->InsertIOElemConNode(mbs_elem_nr, list_idx, inputnr, pos);
				this->ForceUpdate();			
			}
			else if(elem.IsSelectableConstructionNode())
			{
	// specific code for a selected construction node --> delete from both lists
				int mbs_elem_nr = elem.mbs_elnr;
				int list_idx = elem.sub_elnr % 65536;   
				int inputnr = elem.sub_elnr / 65536; 
				
				GetWCDDlg()->GetWCDInterface()->DeleteIOElemConNode(mbs_elem_nr,list_idx);
				this->ForceUpdate();			
			}
		}
	}
}

void CDrawWindow2DDlg::OnCMZoomToFullRange()
{
	GetView()->SetInitialRange();                             // initial range ( real coordinates ) is computed
	GetView()->MatchLogicalPoints();                          // logical points for the view are matched to the extend (pixels)
	SetCaller(mystr("REDRAW triggered by CDrawWindow2DDlg::OnCMZoomToFullRange\n"));
	GetView()->RedrawWindow();
}

void CDrawWindow2DDlg::OnCMExportToFile() 
{ 
	// determine filename - TODO File Save Dialog
	mystr filename;
	CFileDialog fd(FALSE, "jpg", "C:\\temp\\IOBlocks.jpg", OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, "JPG file (*.jpg)|*.jpg|PNG file (*.png)|*.png|BMP file (*.bmp)|*.bmp|EMF file (*.emf)|*.emf||");
	if(fd.DoModal() == IDOK) 
	{
		CString fp = fd.GetFileName();
		CString fn = fd.GetPathName();
		CString fe = fd.GetFileExt();
		filename = mystr(fn);

		int extension = 0;
		if (fe == CString("jpg")) extension = 1;
		if (fe == CString("png")) extension = 2;
		if (fe == CString("bmp")) extension = 4;
		if (fe == CString("emf")) extension = 8;

		CString path_and_filename_noext = fn.Left(fn.GetLength()-4);
		if (extension&(4+2+1))
		{
			CRect rect; GetView()->GetClientRect(rect);
			GetView()->SaveAsBitmapGraphic(path_and_filename_noext, abs(rect.Width() ), abs(rect.Height()), extension&1, extension&2, extension&4);
		}
		if (extension&8)
		{
			GetView()->SaveAsVectorGraphic(path_and_filename_noext, extension&8);
		}
	}
}

void CDrawWindow2DDlg::OnCMPrintGraph() 
{ 
	Print(); 
}

void CDrawWindow2DDlg::OnCMShowStatusInfo() 
{ 
	if(m_flag_statusinfo)
		m_flag_statusinfo = 0;  // change from show to dont show
	else
		m_flag_statusinfo = 1;  // change from dont show to show
	UpdateData(FALSE);        // keep the obsolete Checkbox in correct state ...
}

Vector2D CDrawWindow2DDlg::ComputeSplitPointPosition(ControlWindowContext_::DrawComponent& elem)
{
	Vector2D newconnodepos = (elem.P1()+elem.P2())*0.5;
	return newconnodepos;
}

int CDrawWindow2DDlg::ComputeSplitPointListIndex(ControlWindowContext_::DrawComponent& elem)
{
	int newconnodelistindex = elem.sub_elnr;
	if ( newconnodelistindex < 0 ) newconnodelistindex = (-elem.sub_elnr)+1;
	return newconnodelistindex;
}
