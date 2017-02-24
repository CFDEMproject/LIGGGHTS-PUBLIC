//#**************************************************************
//# filename:             PlotToolDlg.cpp
//#
//# author:               Dorninger
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
#include "PlotToolDlg.h"
#include "myfile.h"
#include "savewindowbitmap.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ owner drawn ComboBoxes:       ColorPicker, LineThicknessPicker, LineStylePicker
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CBColorPicker::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) // based on routine from documentation...
{
	ASSERT(lpDrawItemStruct->CtlType == ODT_COMBOBOX);
  CDC dc;

  dc.Attach(lpDrawItemStruct->hDC);

	// Save these value to restore them when done drawing.
		COLORREF crOldTextColor = dc.GetTextColor();
		COLORREF crOldBkColor = dc.GetBkColor();

  // If this item is selected, set the background color 
  // and the text color to appropriate values. Erase
  // the rect by filling it with the background color.
  if ((lpDrawItemStruct->itemAction & ODA_SELECT) &&
      (lpDrawItemStruct->itemState  & ODS_SELECTED))
  {
     dc.SetTextColor(::GetSysColor(COLOR_HIGHLIGHTTEXT));
     dc.SetBkColor(::GetSysColor(COLOR_HIGHLIGHT));
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, ::GetSysColor(COLOR_HIGHLIGHT));
  }
  else
  {
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, crOldBkColor);
  }

	 // OWN:
	int nr_color = lpDrawItemStruct->itemID; // nr_color is index in combobox 
	if(nr_color >= 0)
	{
		MyColorInfo& colorinfo = palette.Item(nr_color+1);

		CRect entirerect(lpDrawItemStruct->rcItem);
		int shrink = 1;
		
		CRect coloredboxrect(entirerect.left + shrink, entirerect.top + shrink, entirerect.left + shrink + (entirerect.Height() - 2* shrink), entirerect.bottom - shrink);
		dc.FillSolidRect(coloredboxrect, colorinfo.Color_COLREF());

		dc.SetTextColor(colorinfo.Color_COLREF());
    dc.SetBkColor(crOldBkColor);

		CRect textrect(entirerect.left + shrink + entirerect.Height(), entirerect.top, entirerect.right, entirerect.bottom);
		dc.DrawText(colorinfo.Name(), colorinfo.Name().Length(), textrect, DT_SINGLELINE);
	}
	else
    dc.SetBkColor(crOldBkColor);

  dc.Detach();
}

void CBColorPicker::MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct) // based on routine from documentation...
{
  ASSERT(lpMeasureItemStruct->CtlType == ODT_COMBOBOX);

	if (lpMeasureItemStruct->itemID != (UINT) -1)
	{
		LPCTSTR lpszText = (LPCTSTR) lpMeasureItemStruct->itemData;
		ASSERT(lpszText != NULL);
		CSize   sz;
		CDC*    pDC = GetDC();
		sz = pDC->GetTextExtent(lpszText);
		ReleaseDC(pDC);

		// OWN: 
		lpMeasureItemStruct->itemHeight = 2*sz.cy;
	}
}

int CBColorPicker::CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct)
{
	return 0;
}

void CBColorPicker::DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct)
{
	return ;
}


void CBLineStylePicker::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) // based on routine from documentation...
{
	ASSERT(lpDrawItemStruct->CtlType == ODT_COMBOBOX);
  CDC dc;

  dc.Attach(lpDrawItemStruct->hDC);

	// Save these value to restore them when done drawing.
		COLORREF crOldTextColor = dc.GetTextColor();
		COLORREF crOldBkColor = dc.GetBkColor();

  // If this item is selected, set the background color 
  // and the text color to appropriate values. Erase
  // the rect by filling it with the background color.
  if ((lpDrawItemStruct->itemAction & ODA_SELECT) &&
      (lpDrawItemStruct->itemState  & ODS_SELECTED))
  {
     dc.SetTextColor(::GetSysColor(COLOR_HIGHLIGHTTEXT));
     dc.SetBkColor(::GetSysColor(COLOR_HIGHLIGHT));
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, ::GetSysColor(COLOR_HIGHLIGHT));
  }
  else
  {
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, crOldBkColor);
  }

	 // OWN:
	int penstyle = (int) lpDrawItemStruct->itemData;
	if(penstyle >= 0)
	{
		CRect entirerect(lpDrawItemStruct->rcItem);
		int y_mid = (int) (entirerect.top+entirerect.bottom)/2;
		int shrink = 1;
	
		CPen pen(penstyle, 1, activecolor);
		dc.SelectObject(&pen);
		dc.MoveTo( entirerect.left+shrink, y_mid);
		dc.LineTo( entirerect.right-shrink, y_mid);
	}
  dc.SetBkColor(crOldBkColor);

  dc.Detach();
}

void CBLineStylePicker::MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct) // based on routine from documentation...
{
  ASSERT(lpMeasureItemStruct->CtlType == ODT_COMBOBOX);

	if (lpMeasureItemStruct->itemID != (UINT) -1)
	{
		LPCTSTR lpszText = (LPCTSTR) lpMeasureItemStruct->itemData;
		ASSERT(lpszText != NULL);
		CSize   sz;
		CDC*    pDC = GetDC();
		sz = pDC->GetTextExtent(lpszText);
		ReleaseDC(pDC);

		// OWN: 
		lpMeasureItemStruct->itemHeight = 1*sz.cy;
	}
}

int CBLineStylePicker::CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct)
{
	return 0;
}

void CBLineStylePicker::DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct)
{
	return ;
}

void CBPointStylePicker::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) // based on routine from documentation...
{
	ASSERT(lpDrawItemStruct->CtlType == ODT_COMBOBOX);
  CDC dc;

  dc.Attach(lpDrawItemStruct->hDC);

	// Save these value to restore them when done drawing.
		COLORREF crOldTextColor = dc.GetTextColor();
		COLORREF crOldBkColor = dc.GetBkColor();

  // If this item is selected, set the background color 
  // and the text color to appropriate values. Erase
  // the rect by filling it with the background color.
  if ((lpDrawItemStruct->itemAction & ODA_SELECT) &&
      (lpDrawItemStruct->itemState  & ODS_SELECTED))
  {
     dc.SetTextColor(::GetSysColor(COLOR_HIGHLIGHTTEXT));
     dc.SetBkColor(::GetSysColor(COLOR_HIGHLIGHT));
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, ::GetSysColor(COLOR_HIGHLIGHT));
  }
  else
  {
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, crOldBkColor);
  }

	 // OWN:
	int pointstyle = (int) lpDrawItemStruct->itemData;
	if(pointstyle >= 0)
	{
		CRect entirerect(lpDrawItemStruct->rcItem);
		int y_mid = (int) (entirerect.top+entirerect.bottom)/2;
		int x_mid = (int) (entirerect.left+entirerect.right)/2;
		int shrink = 1;
		int extent = (int) ( (entirerect.top-entirerect.bottom)/2.5 +0.5 );
// baseline	- MISSES LINK TO LINE STYLE PICKER
		CPen pen(activelinestyle, 1, activecolor);
		dc.SelectObject(&pen);
		dc.MoveTo( entirerect.left+shrink, y_mid);
		dc.LineTo( entirerect.right-shrink, y_mid);
		CBrush colorbrush(activecolor);
		CBrush whitebrush(crOldBkColor);

    switch (pointstyle)
		{
		case PTS_NONE: 
			break;
		case PTS_X: 
			dc.MoveTo( x_mid-extent, y_mid+extent);
			dc.LineTo( x_mid+extent, y_mid-extent);
			dc.MoveTo( x_mid-extent, y_mid-extent);
			dc.LineTo( x_mid+extent, y_mid+extent);
			break;
		case PTS_CROSS: 
			dc.MoveTo( x_mid+extent, y_mid       );
			dc.LineTo( x_mid-extent, y_mid       );
			dc.MoveTo( x_mid,        y_mid-extent);
			dc.LineTo( x_mid,        y_mid+extent);
			break;
		case PTS_DOT: 
			dc.SelectObject(&colorbrush);
      dc.Ellipse( x_mid-extent, y_mid+extent, x_mid+extent, y_mid-extent);
			break;
		case PTS_CIRCLE: 
			dc.SelectObject(&whitebrush);
      dc.Ellipse( x_mid-extent, y_mid+extent, x_mid+extent, y_mid-extent);
			break;
		default:
			break;
		}
	}
  dc.SetBkColor(crOldBkColor);

  dc.Detach();
}

void CBPointStylePicker::MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct) // based on routine from documentation...
{
  ASSERT(lpMeasureItemStruct->CtlType == ODT_COMBOBOX);

	if (lpMeasureItemStruct->itemID != (UINT) -1)
	{
		LPCTSTR lpszText = (LPCTSTR) lpMeasureItemStruct->itemData;
		ASSERT(lpszText != NULL);
		CSize   sz;
		CDC*    pDC = GetDC();
		sz = pDC->GetTextExtent(lpszText);
		ReleaseDC(pDC);

		// OWN: 
		lpMeasureItemStruct->itemHeight = 1*sz.cy;
	}
}

int CBPointStylePicker::CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct)
{
	return 0;
}

void CBPointStylePicker::DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct)
{
	return ;
}

void CBLineThicknessPicker::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) // based on routine from documentation...
{
	ASSERT(lpDrawItemStruct->CtlType == ODT_COMBOBOX);
  CDC dc;

  dc.Attach(lpDrawItemStruct->hDC);

	// Save these value to restore them when done drawing.
		COLORREF crOldTextColor = dc.GetTextColor();
		COLORREF crOldBkColor = dc.GetBkColor();

  // If this item is selected, set the background color 
  // and the text color to appropriate values. Erase
  // the rect by filling it with the background color.
  if ((lpDrawItemStruct->itemAction & ODA_SELECT) &&
      (lpDrawItemStruct->itemState  & ODS_SELECTED))
  {
     dc.SetTextColor(::GetSysColor(COLOR_HIGHLIGHTTEXT));
     dc.SetBkColor(::GetSysColor(COLOR_HIGHLIGHT));
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, ::GetSysColor(COLOR_HIGHLIGHT));
  }
  else
  {
     dc.FillSolidRect(&lpDrawItemStruct->rcItem, crOldBkColor);
  }

	 // OWN:
	int penthickness = (TPenWidth) lpDrawItemStruct->itemData;
	if(penthickness >= 0)
	{
		CRect entirerect(lpDrawItemStruct->rcItem);
		int y_mid = (int) (entirerect.top+entirerect.bottom)/2;
		int shrink = 1;
	
		CPen pen(PS_SOLID, (int) (sqrt((double)penthickness)+0.5), activecolor);
		dc.SelectObject(&pen);
		dc.MoveTo( entirerect.left+shrink, y_mid);
		dc.LineTo( entirerect.right-shrink, y_mid);
	}
  dc.SetBkColor(crOldBkColor);

  dc.Detach();
}

void CBLineThicknessPicker::MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct) // based on routine from documentation...
{
  ASSERT(lpMeasureItemStruct->CtlType == ODT_COMBOBOX);

	if (lpMeasureItemStruct->itemID != (UINT) -1)
	{
		LPCTSTR lpszText = (LPCTSTR) lpMeasureItemStruct->itemData;
		ASSERT(lpszText != NULL);
		CSize   sz;
		CDC*    pDC = GetDC();
		sz = pDC->GetTextExtent(lpszText);
		ReleaseDC(pDC);

		// OWN: 
		lpMeasureItemStruct->itemHeight = 1*sz.cy;
	}
}

int CBLineThicknessPicker::CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct)
{
	return 0;
}

void CBLineThicknessPicker::DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct)
{
	return ;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CPlotToolPopupChild:       this is the dialog that nests the Graph (CMyView object), 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLEMENT_DYNAMIC(CPlotToolPopupChild, MyBaseDialogForView)

CPlotToolPopupChild::CPlotToolPopupChild(CWnd* pParent /*=NULL*/)
: MyBaseDialogForView(pParent)
{
	m_init_done = 0;
}

CPlotToolPopupChild::~CPlotToolPopupChild()
{
}

BOOL CPlotToolPopupChild::OnInitDialog()
{
	CDialog::OnInitDialog();

	// create the View
	CCreateContext cc;
	cc.m_pNewViewClass = RUNTIME_CLASS(CMyPlotToolView);
	cc.m_pCurrentDoc = NULL;
	cc.m_pNewDocTemplate = NULL;
	cc.m_pLastView = NULL;
	cc.m_pCurrentFrame = NULL;

	SetView( (CMyPlotToolView*)CreateNewView(&cc, this, CRect(0, 0, 0, 0), ID_VIEW_ON_POPUP) );
	if (m_pMyView == NULL)
		EndDialog(IDCANCEL);
	GetView()->SetPlotToolDlg(this->m_pParent);

	CRect rect;
  this->GetWindowRect(&rect);
	rect.top = rect.bottom- 25;
	int dbg = this->m_StatusBar.Create(WS_CHILD | WS_BORDER | WS_VISIBLE, rect, this, IDC_PLOTTOOLPOPUP_STATUSBAR);

	// place the View in the client rectangle
	PlaceElements();
	GetView()->Reset();
	m_init_done = 1;
	slow_zoom = 0;
	_isclosing = 0;

// AD 2013-01-17: this is now in the context menu
	////// add item so SystemMenu
	////CMenu* sysmenu = this->GetSystemMenu(FALSE);
	////sysmenu->AppendMenu(MF_SEPARATOR);
	////sysmenu->AppendMenu(MF_STRING,ID_MENU_PARENT_ON_POPUP,"Show Parent");

	AllowRedraw();
	SetCaller(mystr("REDRAW triggered by CPlotToolPopupChild::OnInitDialog\n"));
	this->UpdateWindow(); // this is the first redraw
  this->SetDisplayName(CString("HOTINT Plot Tool"));

	return TRUE;  // return TRUE unless you set the focus to a control
	// AUSNAHME: OCX-Eigenschaftenseite muss FALSE zurückgeben.
}

BOOL CPlotToolPopupChild::DestroyWindow()
{
	// TODO: Fügen Sie hier Ihren spezialisierten Code ein, und/oder rufen Sie die Basisklasse auf.
	if (GetView())
		GetView()->DestroyWindow();
	GetParent()->m_hasPopupChild = 0;
	return CDialog::DestroyWindow();

	// PostNcDestroy() delete the pointer
}

//void CPlotToolPopupChild::DoDataExchange(CDataExchange* pDX)
//{
//	CDialog::DoDataExchange(pDX);
//}

BEGIN_MESSAGE_MAP(CPlotToolPopupChild, MyBaseDialogForView /*CDialog*/)
	ON_WM_PAINT()
	ON_WM_CLOSE()
	ON_WM_SIZE()
	ON_WM_GETMINMAXINFO()
	ON_WM_MOVE()
	ON_WM_LBUTTONUP()
	ON_WM_SYSCOMMAND()
	ON_WM_MOUSEWHEEL()
// CONTEXTMENU	
	ON_WM_CONTEXTMENU()
	ON_COMMAND(ID_CM_ZOOMTOFULLRANGE,		&CPlotToolPopupChild::OnCMZoomToFullRange)
	ON_COMMAND(ID_CM_SHOWLEGEND,				&CPlotToolPopupChild::OnCMShowLegend)
	ON_COMMAND(ID_CM_INFOPLOTTOOL,			&CPlotToolPopupChild::OnCMShowStatusInfo)
	ON_COMMAND(ID_CM_FASTDRAW,					&CPlotToolPopupChild::OnCMFastDraw)
	ON_COMMAND(ID_CM_AUTOUPDATE,				&CPlotToolPopupChild::OnCMAutoUpdate)
	ON_COMMAND(ID_CM_AUTORESCALE,				&CPlotToolPopupChild::OnCMAutoReScale)
	ON_COMMAND(ID_CM_DRAWSPARSE,				&CPlotToolPopupChild::OnCMDrawSparse)
	ON_COMMAND(ID_CM_MARKSPARSE,				&CPlotToolPopupChild::OnCMMarkSparse)
	ON_COMMAND(ID_CM_MARKCURRENT,				&CPlotToolPopupChild::OnCMMarkCurrent)
	ON_COMMAND(ID_CM_DRAWUPTOCURRENT,		&CPlotToolPopupChild::OnCMDrawUpToCurrent)
	ON_COMMAND(ID_CM_SHOWDIALOG,				&CPlotToolPopupChild::OnCMShowDialog)
	ON_COMMAND(ID_CM_HIDEDIALOG,				&CPlotToolPopupChild::OnCMHideDialog)
	ON_COMMAND(ID_CM_AXISEQUAL,					&CPlotToolPopupChild::OnCMAxisEqual)
	ON_COMMAND(ID_CM_EXPORTTOFILE,			&CPlotToolPopupChild::OnCMExportToFile)
	ON_COMMAND(ID_CM_PRINTGRAPH,				&CPlotToolPopupChild::OnCMPrintGraph)
//END CONTEXTMENU

END_MESSAGE_MAP()

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
void CPlotToolPopupChild::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);

	if (m_init_done) PlaceElements();
}

void CPlotToolPopupChild::OnGetMinMaxInfo(MINMAXINFO* lpMMI)
{
  // set the minimum tracking width
  // and the minimum tracking height of the window
  lpMMI->ptMinTrackSize.x = 100;
  lpMMI->ptMinTrackSize.y = 100;
}

void CPlotToolPopupChild::OnMove(int x, int y)
{
	if (m_init_done) PlaceElements();
	CDialog::OnMove(x, y);
	this->RedrawWindow(0,0,RDW_INVALIDATE|RDW_UPDATENOW|RDW_ERASE|RDW_FRAME);   // RDW_FRAME redraws the NonClient area too
}

void CPlotToolPopupChild::OnLButtonUp(UINT nFlags, CPoint point)
{
	this->GetView()->RedrawWindow();
}

void CPlotToolPopupChild::OnClose()
{	
	if(_isclosing) // dialog is already closing - prevent cascade
		return;
	_isclosing = 1;
// AD: why is the view illegal for an extern close ? but legal when either plottool dialog is closed directly ?
	CMyPlotToolView* view = GetView();    
// this if does saveguard
	if(::IsWindow(view->m_hWnd))          
		view->WriteData();

	if(::IsWindowVisible(GetParent()->m_hWnd)) // PlotToolDlg is visible -> HIDE View
	{
		this->ShowWindow(SW_HIDE);
	}
	else // PlotToolDlg is not visible -> DESTROY PlotTool, then DESTROY View
	{
		GetParent()->OnClose();
		CDialog::DestroyWindow();
		CDialog::OnClose();
	}
}

void CPlotToolPopupChild::OnPaint()
{
	SetCaller(mystr("REDRAW triggered by CPlotToolPopupChild::OnPaint\n"));
	GetView()->UpdateWindow();
}

BOOL CPlotToolPopupChild::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	return GetView()->OnMouseWheel(nFlags, zDelta, pt);
}

BOOL CPlotToolPopupChild::OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult)
{
	if (message == WM_SIZING)
	{
		this->GetView()->skipdrawduringmove = 1;
	}
	else if (message == WM_MOVING)
	{
		this->GetView()->skipdrawduringmove = 1;
	}
	else if(message == WM_EXITSIZEMOVE) // catch end of resize or move for additional redraw
	{
		this->GetView()->skipdrawduringmove = 0;
		this->GetView()->RedrawWindow();
	}
	else if (message == WM_MOUSELEAVE)
	{
// not used
	}
	else if (message == WM_NCMOUSELEAVE)
	{
// not used
	}

	return MyBaseDialogForView::OnWndMsg(message, wParam, lParam, pResult);
}

// ***************************************************************************************
// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
// ***************************************************************************************
void CPlotToolPopupChild::PlaceElements()
{
	if(m_init_done)
	{
		// allign with PlotTool MainDialog (controls, Lists, etc... )
		CRect parentrect,childrect;
		GetParent()->GetWindowRect(&parentrect);
		this->GetWindowRect(&childrect);
		int vshift = childrect.top-parentrect.top;
		int hshift = childrect.left-parentrect.right;
		parentrect.top += vshift; 
		parentrect.bottom += vshift;
		parentrect.left += hshift;
		parentrect.right += hshift;

		GetParent()->MoveWindow(parentrect);
	}

// embedded view on dialog
	CRect clientrect;
	this->GetClientRect(&clientrect);
	CRect viewrect;
	this->GetViewRect(viewrect);
	GetView()->MoveWindow(&viewrect);

// statusbar
	CRect statusrect;
	this->GetStatusRect(statusrect);
	m_StatusBar.MoveWindow(&statusrect);

// two 100 point slots for coordinates at right end with 30 points additional border
	int totalwidth = statusrect.Width();
	int m_widths[4] = {0, totalwidth-230,totalwidth-130,totalwidth-30-1};
	m_StatusBar.SetMinHeight(25);
	m_StatusBar.SetParts(4, m_widths);

	RedrawWindow(0,0,RDW_FRAME|RDW_INVALIDATE);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CGraphExportOptionsDlg:       this is the dialog that nests the Export Options, 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLEMENT_DYNAMIC(CGraphExportOptionsDlg, CDialog)

// "global" function and buffer for the starting path in the BrowseInfoDialog
char pathbuffer[1024];
int CALLBACK BrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData)
{
 	// If the BFFM_INITIALIZED message is received
	// set the path to the start path.
	switch (uMsg)
	{
		case BFFM_INITIALIZED:
		{

			if (NULL != lpData)
			{
				::SendMessage(hwnd, BFFM_SETSELECTION, TRUE, (LPARAM) pathbuffer);
			}
		}
	} 
	return 0;
}


CGraphExportOptionsDlg::CGraphExportOptionsDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CGraphExportOptionsDlg::IDD, pParent)
{
}

CGraphExportOptionsDlg::~CGraphExportOptionsDlg()
{
}

// GraphExportOptionsDlg-Meldungshandler
BOOL CGraphExportOptionsDlg::OnInitDialog()
{
	CDialog::OnInitDialog();
	Reset();
	LoadData();
	return true;
}

void CGraphExportOptionsDlg::Reset()
{
	m_export_path.Empty();   
	m_export_filename.Empty();
  m_pixels_horizontal = 0;
	m_pixels_vertical = 0;
	m_flag_export_jpg = 0;
	m_flag_export_png = 0;
	m_flag_export_bmp = 0;
	m_flag_export_emf = 0;

	//m_snapnow.SetBitmap(LoadBitmap(AfxGetInstanceHandle(), MAKEINTRESOURCE(IDB_SAVE_32)));
}

void CGraphExportOptionsDlg::LoadData()
{
	m_export_path = GetParent()->m_export_path;   
	m_export_filename = GetParent()->m_export_filename;
  m_pixels_horizontal = GetParent()->m_pixels_horizontal;
	m_pixels_vertical = GetParent()->m_pixels_vertical;
	m_flag_export_jpg = GetParent()->m_flag_export_jpg;
	m_flag_export_png = GetParent()->m_flag_export_png;
	m_flag_export_bmp = GetParent()->m_flag_export_bmp;
	m_flag_export_emf = GetParent()->m_flag_export_emf;
	UpdateData(FALSE);
	sprintf(pathbuffer, m_export_path);
}

void CGraphExportOptionsDlg::WriteData()
{
	GetParent()->m_export_path = m_export_path;   
	GetParent()->m_export_filename = m_export_filename;
  GetParent()->m_pixels_horizontal = m_pixels_horizontal;
	GetParent()->m_pixels_vertical = m_pixels_vertical;
	GetParent()->m_flag_export_jpg = m_flag_export_jpg;
	GetParent()->m_flag_export_png = m_flag_export_png;
	GetParent()->m_flag_export_bmp = m_flag_export_bmp;
	GetParent()->m_flag_export_emf = m_flag_export_emf;
	UpdateData(TRUE);
}


void CGraphExportOptionsDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	// *** Group: Export
	DDX_Text(pDX, IDC_EDIT_EXPORT_FILEPATH, m_export_path);
	DDX_Text(pDX, IDC_EDIT_EXPORT_FILENAME, m_export_filename);
	DDX_Text(pDX, IDC_EDIT_RESX, m_pixels_horizontal);
	DDX_Text(pDX, IDC_EDIT_RESY, m_pixels_vertical);
	DDX_Check(pDX, IDC_CHECK_EXPORT_JPG_SIZE, m_flag_export_jpg);
	DDX_Check(pDX, IDC_CHECK_EXPORT_PNG_SIZE, m_flag_export_png);
	DDX_Check(pDX, IDC_CHECK_EXPORT_BMP_SIZE, m_flag_export_bmp);
	DDX_Check(pDX, IDC_CHECK_EXPORT_EMF_SIZE, m_flag_export_emf);

	CString resolution = CString(mystr(" (") + mystr(m_pixels_horizontal) + mystr("x") + mystr(m_pixels_vertical) + mystr(")"));
	CString dummy;
	dummy = CString("JPEG") + resolution;
	DDX_Text(pDX, IDC_CHECK_EXPORT_JPG_SIZE, dummy);
	dummy = CString("PNG") + resolution;
	DDX_Text(pDX, IDC_CHECK_EXPORT_PNG_SIZE, dummy);
	dummy = CString("BMP") + resolution;
	DDX_Text(pDX, IDC_CHECK_EXPORT_BMP_SIZE, dummy);

	DDX_Control(pDX, IDC_BUTTON_SNAPSHOTNOW, m_snapnow);
}

BEGIN_MESSAGE_MAP(CGraphExportOptionsDlg, CDialog)
  // Group: Export
	ON_EN_KILLFOCUS(IDC_EDIT_RESX, &CGraphExportOptionsDlg::OnEnKillfocusEditResx)
	ON_EN_KILLFOCUS(IDC_EDIT_RESY, &CGraphExportOptionsDlg::OnEnKillfocusEditResy)
	ON_BN_CLICKED(IDC_BUTTON_FIGURESIZEASSCREEN, &CGraphExportOptionsDlg::OnBnClickedButtonFiguresizeasscreen)
	ON_BN_CLICKED(IDC_BUTTON_PATHBROWSE, &CGraphExportOptionsDlg::OnBnClickedButtonPathbrowse)
//	ON_BN_CLICKED(IDOK, &CGraphExportOptionsDlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON_SNAPSHOTNOW, &CGraphExportOptionsDlg::OnBnClickedButtonSnapshotnow)
END_MESSAGE_MAP()

void CGraphExportOptionsDlg::OnEnKillfocusEditResx()
{
	UpdateData(TRUE);    // read the edit boxes
	UpdateData(FALSE);   // write the labels
}

void CGraphExportOptionsDlg::OnEnKillfocusEditResy()
{
	UpdateData(TRUE);    // read the edit boxes
	UpdateData(FALSE);   // write the labels
}

void CGraphExportOptionsDlg::OnBnClickedButtonFiguresizeasscreen()
{
	CRect rect;	
	GetParent()->GetPopupChild()->GetView()->GetWindowRect(rect);
	m_pixels_horizontal = rect.Width();
	m_pixels_vertical = rect.Height();
	UpdateData(FALSE);
}

void CGraphExportOptionsDlg::OnBnClickedButtonPathbrowse()
{
	UpdateData(TRUE);
	const CString dummy= (const CString) m_export_path;
	char buffer[1024];
	sprintf(buffer,dummy);
  
	CString path(m_export_path);
	BROWSEINFO bi;
	bi.hwndOwner = m_hWnd; // Handle of the owner window
	bi.pidlRoot = NULL; // Desktop folder is used
	bi.lpszTitle = "Select output path ";
	bi.pszDisplayName = buffer; // Buffer for selected folder name
//	bi.ulFlags = BIF_RETURNONLYFSDIRS; // Only returns file system directories
//	bi.ulFlags = BIF_RETURNONLYFSDIRS|BIF_NEWDIALOGSTYLE; // <--- this should create the NEW FOLDER button
	bi.ulFlags = BIF_EDITBOX | BIF_VALIDATE /*| BIF_BROWSEFORPRINTER*/ | BIF_RETURNONLYFSDIRS|0x0040; // <--- and this one works... :)
	bi.lpfn = BrowseCallbackProc;
	bi.lParam = (LPARAM) buffer;

	LPITEMIDLIST pItemIDList = SHBrowseForFolder(&bi);
	if (pItemIDList)
	{
		if (SHGetPathFromIDList(pItemIDList,buffer)) 	
		{
			m_export_path = buffer;
			sprintf(pathbuffer,buffer);
		}
	}
	UpdateData(FALSE);
}

void CGraphExportOptionsDlg::OnBnClickedOk()
{
	WriteData();
	OnOK();
}

void CGraphExportOptionsDlg::OnBnClickedButtonSnapshotnow()
{
	UpdateData(TRUE);
	WriteData();
	GetParent()->ExportToFile();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CPlotToolDlg:       this is the dialog that nests the controls for the graph, 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLEMENT_DYNAMIC(CPlotToolDlg, CDialog)

CPlotToolDlg::CPlotToolDlg(CWnd* pParent /*=NULL*/)
: CDialog(CPlotToolDlg::IDD, pParent)
, m_radio_source(FALSE)
{
	m_init_done = 0;
	m_hasPopupChild = 0;
	SetWCDDlg( (CWCDriver3DDlg*) pParent);
	solfile = NULL;
	extfile = NULL;
}

CPlotToolDlg::~CPlotToolDlg()
{
}

// CPlotToolDlg-Meldungshandler
BOOL CPlotToolDlg::OnInitDialog()
{
	CDialog::OnInitDialog();
	Reset();
	LoadData();
	StartStopAutoUpdate(); // start or dont start refresh timer

	AssignOtherFonts();

	CCreateContext cc;
	cc.m_pNewViewClass = RUNTIME_CLASS(CMyPlotToolView);
	cc.m_pCurrentDoc = NULL;
	cc.m_pNewDocTemplate = NULL;
	cc.m_pCurrentFrame = NULL;

#ifdef preview
	SetPreview( (CMyView*)CreateNewView(&cc, this, CRect(0, 0, 0, 0), 1234) );
	if (m_pMyView == NULL)
		GetPreview(IDCANCEL);
#endif

	SetPopupChild( new CPlotToolPopupChild() );
	GetPopupChild()->SetParent(this); // before Create (::Create uses GetDoc() )
	if (!::IsWindow(GetPopupChild()->m_hWnd) )
	{
		GetPopupChild()->Create(IDD_DIALOG_2D_VIEW,this);
	}
	// initial size from EDC
	CRect childrect;
	CRect clientrect;  // <-- this shall match the values from EDC
	GetPopupChild()->GetWindowRect(&childrect);
	GetPopupChild()->GetClientRect(&clientrect);

	int border_x = childrect.Width() - clientrect.Width();
	int border_y = childrect.Height() - clientrect.Height();

	childrect.right = childrect.left + pWCDI->GetIOption(260) + border_x; // initial horizontal size 
	childrect.bottom = childrect.top + pWCDI->GetIOption(261) + border_y; // initial vertical size
	GetPopupChild()->MoveWindow(childrect);

	GetPopupChild()->ShowWindow(SW_SHOW);
	m_hasPopupChild = 1;

	// place the View in the view rectangle
	PlaceElements();
	m_init_done = 1;
	_isclosing = 0;
	_externclose = 0;
	_isloading = 0;

	UpdateData(FALSE);   // write e.g. bmp dimensions to dialog
	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnInitDialog\n"));
	this->UpdateWindow();
	return TRUE;  // return TRUE unless you set the focus to a control
	// AUSNAHME: OCX-Eigenschaftenseite muss FALSE zurückgeben.
}

BOOL CPlotToolDlg::DestroyWindow()
{
	// TODO: Fügen Sie hier Ihren spezialisierten Code ein, und/oder rufen Sie die Basisklasse auf.
#ifdef preview
	if (GetPreview()) m_pMyView->DestroyWindow();
#endif
	if (m_hasPopupChild)
		GetPopupChild()->DestroyWindow();

	return CDialog::DestroyWindow();
}

// *********************************
// COMMUNICATION WITH SETTINGs (EDC)
// *********************************
void CPlotToolDlg::Reset()
{
	// default font
	LOGFONT defaultfont;
	defaultfont.lfHeight = -160; // ~ size 12
	defaultfont.lfWidth = 0;
	defaultfont.lfEscapement = 0;
	defaultfont.lfOrientation = 0;
	defaultfont.lfWeight = 400;
	defaultfont.lfItalic = 0;
	defaultfont.lfUnderline = 0;
	defaultfont.lfStrikeOut = 0;
	defaultfont.lfCharSet = 0;
	defaultfont.lfOutPrecision = 3;
	defaultfont.lfClipPrecision = 2;
	defaultfont.lfQuality = 1;
	defaultfont.lfPitchAndFamily = 34;
	strcpy(defaultfont.lfFaceName,"Arial");
	m_axisfont = defaultfont;

	// Data Sources
	m_radio_source = 0;   
	
	SetSolFileNameFromEDC();
	if(solfile!=NULL) delete solfile;
	solfile_linesread = LoadSolutionToMatrixFile();
	
	filename_extern.Empty();
	
	if(m_list_columns.m_hWnd)
		m_list_columns.ResetContent();

	// Auxiliary Lines
//	m_auxline_function.Empty();

	// GeneralOptions
	m_flag_use_autoupdate = 0;
	m_update_interval = 1.0;
	m_flag_autorescale = 0;
	m_flag_show_legend = 0;
	m_flag_statusinfo = 0;
	m_flag_fastdraw = 0;

	// Data Sets to be drawn
	if(m_list_lines.m_hWnd)
	{
		m_list_lines.DeleteAllItems();
		m_list_lines.DeleteColumn(0);
		m_list_lines.InsertColumn(0,"Line Name",LVCFMT_CENTER,200,0);
	}
	lines.Flush();

	// properties of selected line
	m_edit_name_of_selection.Empty();
	// populate color picker
	line_color.Clear();
	palette.Reset();
	palette.AddDefaultColors(8);
	line_color.InitializePalette(palette);
	line_color_index = -1; //line_color.FindString(-1,"black");
	// populate thickness picker
	line_thickness.Clear(); // same sequencd as in TPenWidth
	line_thickness.SetItemData( line_thickness.AddString("very fine"), PT_VERYFINE );
	line_thickness.SetItemData( line_thickness.AddString("fine"), PT_FINE );
	line_thickness.SetItemData( line_thickness.AddString("normal"), PT_NORMAL );
	line_thickness.SetItemData( line_thickness.AddString("thick"), PT_THICK );
	line_thickness.SetItemData( line_thickness.AddString("very thick"), PT_VERYTHICK );
	line_thickness_index = line_thickness.FindString(-1,"fine");
	// populate line style picker
	line_style.Initialize();
	line_style_index = line_style.FindString(-1,"solid");
	//polpulate point style picker
	point_style.Initialize();
	point_style_index = line_style.FindString(-1,"none");
	point_size = 30;

	// Title and Axis
	m_title = CString(""); //("Title");

	m_x_axis_title = CString("");//("X-axis");
	x_min = 0.0;
	x_max = 0.0;
	flag_use_x_majorticks = 0;
	x_major_interval = 0.0;
	flag_use_x_minorticks = 0;
	x_minor_interval = 0.0;

	m_y_axis_title = CString("");//("Y-axis");
	y_min = 0.0;
	y_max = 0.0;
	flag_use_y_majorticks = 0;
	y_major_interval = 0.0;
	flag_use_y_minorticks = 0;
	y_minor_interval = 0.0;

	UpdateData(FALSE);
	OnBnClickedRadioSource();
	m_selected_line = 0;
	WriteToDlgPropertiesofLine(0);

	StartStopAutoUpdate();


	////CTabCtrl* TC_TabCtrl = (CTabCtrl*) GetDlgItem(IDC_TAB_PLOTTOOL);
	////TC_TabCtrl->DeleteAllItems();
	////TC_TabCtrl->InsertItem(0,"Labels");
	////TC_TabCtrl->InsertItem(1,"General");
	////TC_TabCtrl->InsertItem(2,"Export");
	////TC_TabCtrl->InsertItem(3,"Buttons");
	//////TC_TabCtrl->InsertItem(4,"Aux .Line");
	////TC_TabCtrl->SetCurSel(1);
	////	
	////DrawActiveTabElements();



	// disable all non-implemented dialog items (this list should become empty with time...)

	//GetDlgItem(IDC_CHECK_XAXIS_MAJOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_MAJOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_CHECK_XAXIS_MINOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_MINOR)->EnableWindow(FALSE);

	//GetDlgItem(IDC_CHECK_YAXIS_MAJOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_MAJOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_CHECK_YAXIS_MINOR)->EnableWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_MINOR)->EnableWindow(FALSE);

	/// HACK AD:: all dialog items not needed  // "to be removed" are made invisible for the time being !!!

	//GetDlgItem(IDC_BUTTON_FONT_XAXIS)->ShowWindow(FALSE);
	//GetDlgItem(IDC_STATIC_LIMITX)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_LOWER)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_UPPER)->ShowWindow(FALSE);
	//GetDlgItem(IDC_CHECK_XAXIS_MAJOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_MAJOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_CHECK_XAXIS_MINOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_XAXIS_MINOR)->ShowWindow(FALSE);

	//GetDlgItem(IDC_BUTTON_FONT_YAXIS)->ShowWindow(FALSE);
	//GetDlgItem(IDC_STATIC_LIMITY)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_LOWER)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_UPPER)->ShowWindow(FALSE);
	//GetDlgItem(IDC_CHECK_YAXIS_MAJOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_MAJOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_CHECK_YAXIS_MINOR)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_YAXIS_MINOR)->ShowWindow(FALSE);

	// Elements for aux line (parsed function)
	//GetDlgItem(IDC_STATIC_AUXLINE)->ShowWindow(FALSE);
	//GetDlgItem(IDC_EDIT_AUXLINE)->ShowWindow(FALSE);
	//GetDlgItem(IDC_BUTTON_ADDAUXLINE)->ShowWindow(FALSE);
}

// Load Settings From EDC
void CPlotToolDlg::LoadData()
{
	// general options
	pWCDI->MBS_EDC_TreeGetInt(m_flag_use_autoupdate,"PlotToolOptions.auto_redraw");
	pWCDI->MBS_EDC_TreeGetDouble(m_update_interval,"PlotToolOptions.auto_redraw_interval");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_autorescale,"PlotToolOptions.auto_rescale");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_show_legend,"PlotToolOptions.Legend.show");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_statusinfo,"PlotToolOptions.status_bar_info");
	
// data point display 
	pWCDI->MBS_EDC_TreeGetInt(m_flag_draw_every_nth,"PlotToolOptions.DataPoints.flag_draw_every_nth");
	pWCDI->MBS_EDC_TreeGetInt(m_draw_every_nth,"PlotToolOptions.DataPoints.draw_every_nth");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_mark_every_nth,"PlotToolOptions.DataPoints.flag_mark_every_nth");
	pWCDI->MBS_EDC_TreeGetInt(m_mark_every_nth,"PlotToolOptions.DataPoints.mark_every_nth");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_vertical_marker,"PlotToolOptions.DataPoints.vertical_marker");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_upto_current,"PlotToolOptions.DataPoints.draw_only_to_time");

	// relative size of the fonts
	pWCDI->MBS_EDC_TreeGetDouble(m_titlefont_scalingfactor,"PlotToolOptions.title_size_factor");
	pWCDI->MBS_EDC_TreeGetDouble(m_ticksfont_scalingfactor,"PlotToolOptions.ticks_size_factor");

  	// axis 
	pWCDI->MBS_EDC_TreeGetInt(m_border_thickness,"PlotToolOptions.line_thickness_border");
	pWCDI->MBS_EDC_TreeGetInt(m_axis_at_origin,"PlotToolOptions.Axis.draw_at_origin");
	pWCDI->MBS_EDC_TreeGetDouble(m_axis_overdraw,"PlotToolOptions.Axis.overdraw");

	pWCDI->MBS_EDC_TreeGetInt(m_axis_label_major,"PlotToolOptions.Axis.label_major");
	pWCDI->MBS_EDC_TreeGetInt(m_axis_label_minor,"PlotToolOptions.Axis.label_minor");
	pWCDI->MBS_EDC_TreeGetDouble(m_axis_ticks_overdraw,"PlotToolOptions.Axis.overdraw");
	pWCDI->MBS_EDC_TreeGetInt(m_digits_x_label,"PlotToolOptions.Axis.digits_x_labels");
	pWCDI->MBS_EDC_TreeGetInt(m_digits_y_label,"PlotToolOptions.Axis.digits_y_labels");

	pWCDI->MBS_EDC_TreeGetInt(m_axis_minor_ticks_x,"PlotToolOptions.Axis.minor_ticks_x");
	pWCDI->MBS_EDC_TreeGetInt(m_axis_minor_ticks_y,"PlotToolOptions.Axis.minor_ticks_y");

	// export
	m_export_path = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.plottool_image_path");
	m_export_filename = pWCDI->MBS_EDC_TreeGetString("PlotToolOptions.SavePicture.filename");

	pWCDI->MBS_EDC_TreeGetInt(m_pixels_horizontal,"PlotToolOptions.SavePicture.size_horizontal");
	pWCDI->MBS_EDC_TreeGetInt(m_pixels_vertical,"PlotToolOptions.SavePicture.size_vertical");
	
	pWCDI->MBS_EDC_TreeGetInt(m_flag_export_jpg,"PlotToolOptions.SavePicture.store_jpg");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_export_png,"PlotToolOptions.SavePicture.store_png");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_export_bmp,"PlotToolOptions.SavePicture.store_bmp");
	pWCDI->MBS_EDC_TreeGetInt(m_flag_export_emf,"PlotToolOptions.SavePicture.store_emf");
	
	UpdateData(FALSE);
}

// write changable settings back to EDC
void CPlotToolDlg::WriteData()
{
	if(!_externclose)
		UpdateData(TRUE);
	// general options
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_use_autoupdate,"PlotToolOptions.auto_redraw");
	//pWCDI->MBS_EDC_TreeSetDouble(m_update_interval,"PlotToolOptions.auto_redraw_interval");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_autorescale,"PlotToolOptions.auto_rescale");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_show_legend,"PlotToolOptions.Legend.show");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_statusinfo,"PlotToolOptions.status_bar_info");
	pWCDI->SetIOption(255,m_flag_use_autoupdate);
	pWCDI->SetDOption(255,m_update_interval);
	pWCDI->SetIOption(257,m_flag_autorescale);
	pWCDI->SetIOption(275,m_flag_show_legend);
	pWCDI->SetIOption(259,m_flag_statusinfo);
	
	// datapoints
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_draw_every_nth,"PlotToolOptions.DataPoints.flag_draw_every_nth");
	//pWCDI->MBS_EDC_TreeSetInt(m_draw_every_nth,"PlotToolOptions.DataPoints.draw_every_nth");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_mark_every_nth,"PlotToolOptions.DataPoints.flag_mark_every_nth");
	//pWCDI->MBS_EDC_TreeSetInt(m_mark_every_nth,"PlotToolOptions.DataPoints.mark_every_nth");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_vertical_marker,"PlotToolOptions.DataPoints.vertical_marker");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_upto_current,"PlotToolOptions.DataPoints.draw_only_to_time");
	pWCDI->SetIOption(290,m_flag_draw_every_nth);
	pWCDI->SetIOption(291,m_draw_every_nth);
	pWCDI->SetIOption(292,m_flag_mark_every_nth);
	pWCDI->SetIOption(293,m_mark_every_nth);
	pWCDI->SetIOption(294,m_flag_vertical_marker);
	pWCDI->SetIOption(295,m_flag_upto_current);

		// export
	//pWCDI->MBS_EDC_TreeSetString(m_export_path,"GeneralOptions.Paths.plottool_image_path");
	//pWCDI->MBS_EDC_TreeSetString(m_export_filename,"PlotToolOptions.SavePicture.filename");
	pWCDI->SetTOption(124,m_export_path);
	pWCDI->SetTOption(125,m_export_filename);

	//pWCDI->MBS_EDC_TreeSetInt(m_pixels_horizontal,"PlotToolOptions.SavePicture.size_horizontal");
	//pWCDI->MBS_EDC_TreeSetInt(m_pixels_vertical,"PlotToolOptions.SavePicture.size_vertical");
	pWCDI->SetIOption(280,m_pixels_horizontal);
	pWCDI->SetIOption(281,m_pixels_vertical);

	//pWCDI->MBS_EDC_TreeSetInt(m_flag_export_jpg,"PlotToolOptions.SavePicture.store_jpg");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_export_png,"PlotToolOptions.SavePicture.store_png");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_export_bmp,"PlotToolOptions.SavePicture.store_bmp");
	//pWCDI->MBS_EDC_TreeSetInt(m_flag_export_emf,"PlotToolOptions.SavePicture.store_emf");
	pWCDI->SetIOption(285,m_flag_export_jpg);
	pWCDI->SetIOption(286,m_flag_export_png);
	pWCDI->SetIOption(287,m_flag_export_bmp);
	pWCDI->SetIOption(288,m_flag_export_emf);

}

// Update from MBS_EDC in case the edc is changed while PlotToolDialog is open
void CPlotToolDlg::UpdateDialogData()
{
  LoadData();
	if(GetPopupChild()->GetView()!=NULL)
	{
		GetPopupChild()->GetView()->LoadData();
	}
	UpdateData(FALSE);
	ScaleToFull(); // automatic rescale of the graph !
}

// *********************************
// COMMUNICATION WITH SETTINGs (EDC)
// *********************************
void CPlotToolDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	// *** Group: Data Sources
	DDX_Radio(pDX, IDC_RADIO_SENSOR, m_radio_source);
	DDX_Text(pDX, IDC_EDIT_FILENAME, filename_extern);
	DDX_Text(pDX, IDC_LABEL_SOURCE, m_label_source);
	DDX_Control(pDX, IDC_LIST_COLUMNS, m_list_columns);

	// *** ListBox: Available Line
	// no DDX

	// *** Group: Auxillary Lines
//	DDX_Text(pDX, IDC_EDIT_AUXLINE, m_auxline_function);

	// *** Group: Title and Axis
	DDX_Text(pDX, IDC_EDIT_LABEL_TITLE, m_title);
	DDX_Text(pDX, IDC_EDIT_LABEL_XAXIS, m_x_axis_title);
	DDX_Text(pDX, IDC_EDIT_XAXIS_LOWER, x_min);
	DDX_Text(pDX, IDC_EDIT_XAXIS_UPPER, x_max);
	//DDX_Text(pDX, IDC_EDIT_XAXIS_MAJOR, x_major_interval);

	DDX_Text(pDX, IDC_EDIT_LABEL_YAXIS, m_y_axis_title);
	DDX_Text(pDX, IDC_EDIT_YAXIS_LOWER, y_min);
	DDX_Text(pDX, IDC_EDIT_YAXIS_UPPER, y_max);
	//DDX_Text(pDX, IDC_EDIT_YAXIS_MAJOR, y_major_interval);

	// *** ListBox: Drawn Lines
	DDX_Control(pDX, IDC_LIST_LINES, m_list_lines);
	DDX_Text(pDX, IDC_LABEL_LINEPROP, m_label_lineprop);
	DDX_Text(pDX, IDC_EDIT_LINE_NAME, m_edit_name_of_selection);
	DDX_Control(pDX, IDC_COMBO_LINE_COLOR, line_color);
	DDX_CBIndex(pDX, IDC_COMBO_LINE_COLOR, line_color_index);
	DDX_Control(pDX, IDC_COMBO_LINE_THICKNESS, line_thickness);
	DDX_CBIndex(pDX, IDC_COMBO_LINE_THICKNESS, line_thickness_index);
	DDX_Control(pDX, IDC_COMBO_LINE_STYLELINE, line_style);
	DDX_CBIndex(pDX, IDC_COMBO_LINE_STYLELINE, line_style_index);
	DDX_Control(pDX, IDC_COMBO_LINE_STYLEPOINT, point_style);
	DDX_CBIndex(pDX, IDC_COMBO_LINE_STYLEPOINT, point_style_index);
	DDX_Text(pDX, IDC_EDIT_POINTSIZE, point_size);

	// *** Group: General Options
	//DDX_Check(pDX, IDC_CHECK_UPDATE, m_flag_use_autoupdate);
	//DDX_Text(pDX, IDC_EDIT_UPDATE, m_update_interval);
	//DDX_Check(pDX, IDC_CHECK_AUTORESCALE, m_flag_autorescale);
	//DDX_Check(pDX, IDC_CHECK_LEGEND, m_flag_show_legend);
	//DDX_Check(pDX, IDC_CHECK_DRAW_EVERY_NTH, m_flag_draw_every_nth);
	//DDX_Text(pDX, IDC_EDIT_DRAW_EVERY_NTH, m_draw_every_nth);
	//DDX_Check(pDX, IDC_CHECK_VERT_AT_CURRENT, m_flag_vertical_marker);
	//DDX_Check(pDX, IDC_CHECK_DRAW_UPTO_CURRENT, m_flag_upto_current);
	//DDX_Check(pDX, IDC_CHECK_DRAW_FAST, m_flag_fastdraw);

	// *** Group: Export
	// *** NOW IN OWN DIALOG
	//DDX_Text(pDX, IDC_EDIT_EXPORT_FILEPATH, m_export_path);
	//DDX_Text(pDX, IDC_EDIT_EXPORT_FILENAME, m_export_filename);
	//DDX_Text(pDX, IDC_EDIT_RESX, m_pixels_horizontal);
	//DDX_Text(pDX, IDC_EDIT_RESY, m_pixels_vertical);
	//DDX_Check(pDX, IDC_CHECK_EXPORT_JPG_SIZE, m_flag_export_jpg);
	//DDX_Check(pDX, IDC_CHECK_EXPORT_PNG_SIZE, m_flag_export_png);
	//DDX_Check(pDX, IDC_CHECK_EXPORT_BMP_SIZE, m_flag_export_bmp);
	//DDX_Check(pDX, IDC_CHECK_EXPORT_EMF_SIZE, m_flag_export_emf);

	//CString resolution = CString(mystr(" (") + mystr(m_pixels_horizontal) + mystr("x") + mystr(m_pixels_vertical) + mystr(")"));
	//CString dummy;
	//dummy = CString("JPEG") + resolution;
	//DDX_Text(pDX, IDC_CHECK_EXPORT_JPG_SIZE, dummy);
	//dummy = CString("PNG") + resolution;
	//DDX_Text(pDX, IDC_CHECK_EXPORT_PNG_SIZE, dummy);
	//dummy = CString("BMP") + resolution;
	//DDX_Text(pDX, IDC_CHECK_EXPORT_BMP_SIZE, dummy);

	// redraw
//////////	if (this->m_init_done) 
//////////	{
//////////		SetCaller(mystr("REDRAW triggered by CPlotToolDlg::DoDataExchange\n"));
//////////		this->UpdateWindow(); // dialog
//////////		if (!m_do_not_redraw)
//////////		{
//////////			SetCaller(mystr("REDRAW triggered by CPlotToolDlg::DoDataExchange\n"));
//////////#ifdef preview
//////////			if(GetPreview()) GetPreview()->RedrawWindow(); // nested view
//////////#endif
//////////			if(GetPopupChild()!=NULL) GetPopupChild()->GetView()->RedrawWindow(); // popup view
//////////		}
//////////	}
	m_do_not_redraw = 0;
	//DDX_Control(pDX, IDC_TAB_PLOTTOOL, m_tabcontrol);
}

BEGIN_MESSAGE_MAP(CPlotToolDlg, CDialog)
	ON_WM_PAINT()
	ON_WM_CLOSE()
	ON_WM_TIMER()
	ON_MESSAGE(WM_UPDATE_PLOT,OnUpdateMsg)
	ON_WM_SIZE()
	ON_WM_MOVE()
	ON_WM_LBUTTONUP()

	// Group: Data Sources
	ON_BN_CLICKED(IDC_RADIO_SENSOR, &CPlotToolDlg::OnBnClickedRadioSource)
	ON_BN_CLICKED(IDC_RADIO_FILE, &CPlotToolDlg::OnBnClickedRadioSource)
	ON_BN_CLICKED(IDC_BUTTON_FILEBROWSE, &CPlotToolDlg::OnBnClickedButtonFilebrowse)
	ON_BN_CLICKED(IDC_BUTTON_LINEADD, &CPlotToolDlg::OnBnClickedButtonLineAdd)
	ON_LBN_DBLCLK(IDC_LIST_COLUMNS, &CPlotToolDlg::OnLbnDblclkListColumns)
	ON_BN_CLICKED(IDC_BUTTON_LINEADDXY, &CPlotToolDlg::OnBnClickedButtonLineAddXY)
	ON_LBN_SELCHANGE(IDC_LIST_COLUMNS, &CPlotToolDlg::OnLbnSelchangeListColumns)

	// Buttons...
	ON_BN_CLICKED(IDC_BUTTON_QUIT, &CPlotToolDlg::OnClose)

	// Group: Title and Axis
	ON_BN_CLICKED(IDC_BUTTON_FONT, &CPlotToolDlg::OnBnClickedButtonFont)
	ON_EN_KILLFOCUS(IDC_EDIT_LABEL_TITLE, &CPlotToolDlg::OnEnKillfocusEditTitle)
	ON_EN_KILLFOCUS(IDC_EDIT_LABEL_XAXIS, &CPlotToolDlg::OnEnKillfocusEditXAxis)
	ON_EN_KILLFOCUS(IDC_EDIT_LABEL_YAXIS, &CPlotToolDlg::OnEnKillfocusEditYAxis)

	// Group: Data Sets
	ON_NOTIFY(LVN_ITEMCHANGED, IDC_LIST_LINES, &CPlotToolDlg::OnLvnItemchangedListLines)
	ON_NOTIFY(LVN_KEYDOWN, IDC_LIST_LINES, &CPlotToolDlg::OnLvnKeydownListLines)

	// Group: Properties of Line
	ON_EN_KILLFOCUS(IDC_EDIT_LINE_NAME, &CPlotToolDlg::OnEnKillfocusEditLineName)
	ON_CBN_SELENDOK(IDC_COMBO_LINE_THICKNESS, &CPlotToolDlg::OnCbnSelendokComboLineThickness)
	ON_CBN_SELENDOK(IDC_COMBO_LINE_COLOR, &CPlotToolDlg::OnCbnSelendokComboLineColor)
	ON_CBN_SELENDOK(IDC_COMBO_LINE_STYLE, &CPlotToolDlg::OnCbnSelendokComboLineStyle)
	ON_CBN_SELENDOK(IDC_COMBO_LINE_STYLEPOINT, &CPlotToolDlg::OnCbnSelendokComboPointStyle)
	ON_EN_KILLFOCUS(IDC_EDIT_POINTSIZE, &CPlotToolDlg::OnEnKillfocusEditPointsize)

  // Group: Export
	ON_COMMAND(ID_LAYOUTS_SAVE, &CPlotToolDlg::OnLayoutsSave)
	ON_COMMAND(ID_LAYOUTS_LOAD, &CPlotToolDlg::OnLayoutsLoad)
	//ON_EN_KILLFOCUS(IDC_EDIT_RESX, &CPlotToolDlg::OnEnKillfocusEditResx)
	//ON_EN_KILLFOCUS(IDC_EDIT_RESY, &CPlotToolDlg::OnEnKillfocusEditResy)
	//ON_BN_CLICKED(IDC_BUTTON_FIGURESIZEASSCREEN, &CPlotToolDlg::OnBnClickedButtonFiguresizeasscreen)
	//ON_BN_CLICKED(IDC_BUTTON_PATHBROWSE, &CPlotToolDlg::OnBnClickedButtonPathbrowse)

	// CONTEXTMENU
	ON_WM_CONTEXTMENU()
	ON_COMMAND(ID_CM_UPDATEDATA, &CPlotToolDlg::OnCmUpdateData)
	ON_COMMAND(ID_CM_SHOWGRAPH, &CPlotToolDlg::OnCmShowGraph)
	ON_COMMAND(ID_CM_SAVESETTINGS, &CPlotToolDlg::OnCmSaveSettings)
	ON_COMMAND(ID_CM_MARKCURRENT, &CPlotToolDlg::OnCmMarkcurrent)
	
	//END CONTEXTMENU

// DISABLED CONTROLS
	ON_EN_KILLFOCUS(IDC_EDIT_XAXIS_LOWER, &CPlotToolDlg::OnEnKillfocusEditXaxisLower)
	ON_EN_KILLFOCUS(IDC_EDIT_XAXIS_UPPER, &CPlotToolDlg::OnEnKillfocusEditXaxisUpper)
	ON_EN_KILLFOCUS(IDC_EDIT_YAXIS_LOWER, &CPlotToolDlg::OnEnKillfocusEditYaxisLower)
	ON_EN_KILLFOCUS(IDC_EDIT_YAXIS_UPPER, &CPlotToolDlg::OnEnKillfocusEditYaxisUpper)
	//ON_BN_CLICKED(IDC_BUTTON_ADDAUXLINE, &CPlotToolDlg::OnBnClickedButtonAddAuxLine)
//TAB
	//ON_NOTIFY(TCN_SELCHANGE, IDC_TAB_PLOTTOOL, &CPlotToolDlg::OnTcnSelchangeTabPlottool)
// ELEMENTS MOVED TO CONTEXTMENU OF POPUPCHILD
	//ON_BN_CLICKED(IDC_AXISEQUAL, &CPlotToolDlg::OnBnClickedButtonAxisEqual)
	//ON_BN_CLICKED(IDC_CHECK_AUTORESCALE, &CPlotToolDlg::OnBnClickedCheckAutorescale)
	//ON_BN_CLICKED(IDC_CHECK_UPDATE, &CPlotToolDlg::StartStopAutoUpdate)
	//ON_EN_KILLFOCUS(IDC_EDIT_UPDATE, &CPlotToolDlg::OnEnKillfocusEditUpdate)
	//ON_BN_CLICKED(IDC_CHECK_DRAW_EVERY_NTH, &CPlotToolDlg::OnBnClickedCheckDrawEveryNth)
	//ON_EN_KILLFOCUS(IDC_EDIT_DRAW_EVERY_NTH, &CPlotToolDlg::OnEnKillfocusEditDrawEveryNth)
	//ON_BN_CLICKED(IDC_CHECK_LEGEND, &CPlotToolDlg::OnBnClickedCheckLegend)
	//ON_BN_CLICKED(IDC_CHECK_VERT_AT_CURRENT, &CPlotToolDlg::OnBnClickedCheckVertAtCurrent)
	//ON_BN_CLICKED(IDC_CHECK_DRAW_UPTO_CURRENT, &CPlotToolDlg::OnBnClickedCheckDrawUptoCurrent)
	//ON_BN_CLICKED(IDC_BUTTON_EXPORT, &CPlotToolDlg::OnBnClickedButtonExport)
	//ON_BN_CLICKED(IDC_CHECK_DRAW_FAST, &CPlotToolDlg::OnBnClickedCheckDrawFast)
	//ON_BN_CLICKED(IDPRINT, &CPlotToolDlg::OnBnClickedButtonPrint)
	//ON_BN_CLICKED(IDC_BUTTON_SHOW, &CPlotToolDlg::OnBnClickedButtonShow)
	//ON_BN_CLICKED(IDC_BUTTON_REDRAW, &CPlotToolDlg::OnBnClickedButtonRedraw)
	//ON_BN_CLICKED(IDC_UPDATE, &CPlotToolDlg::OnBnClickedUpdate)
	//ON_BN_CLICKED(IDC_BUTTON_SCALE_FULL, &CPlotToolDlg::OnBnClickedButtonScaleFull)

// ELEMENTS MOVED TO CONTEXT MENU
	//ON_BN_CLICKED(IDC_READOPTIONS, &CPlotToolDlg::OnBnClickedReadoptions)
	//ON_BN_CLICKED(IDC_BUTTON_HIDE, &CPlotToolDlg::OnBnClickedButtonHide)

	ON_BN_CLICKED(IDC_BUTTON_EXPORTOPTIONS, &CPlotToolDlg::OnBnClickedButtonExportOptions)
//	ON_BN_CLICKED(IDOK, &CPlotToolDlg::OnBnClickedOk)

END_MESSAGE_MAP()

// ************************************
// *** TIMER FUNCTIONS (redraw, update)
// ************************************
void CPlotToolDlg::StartRefreshTimer()
{
	int delay_ms = (int) (m_update_interval*1000.0); 
	int timer_id = (int) SetTimer(PLOT_REFRESH_TIMER_ID, delay_ms, 0);
	if (!timer_id)
		::MessageBox(this->GetSafeHwnd(),"Unable to start timer","PLOT_REFRESH_TIMER_ID",MB_OK|MB_SYSTEMMODAL);
}

void CPlotToolDlg::KillRefreshTimer()
{
	KillTimer(PLOT_REFRESH_TIMER_ID);
}

void CPlotToolDlg::OnTimer(UINT_PTR nIDEvent)
{
	if (nIDEvent == PLOT_REFRESH_TIMER_ID)
	{
		MSG msg;
		while(::PeekMessage(&msg, m_hWnd, PLOT_REFRESH_TIMER_ID, PLOT_REFRESH_TIMER_ID, PM_REMOVE));
	}
	if(m_flag_use_autoupdate)
	{
		this->Redraw();
	}
	CDialog::OnTimer(nIDEvent);
}

LRESULT CPlotToolDlg::OnUpdateMsg(WPARAM mode, LPARAM str)
{
	for(int i=1; i<=lines.Length(); i++)
	{
		if (lines(i).IsSensor()) lines(i).ReadLine();
	}
	Redraw();
	return 0;
}

// ***************************************************
// *** WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***************************************************
void CPlotToolDlg::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);
	if(m_init_done) PlaceElements();
}

void CPlotToolDlg::OnMove(int x, int y)
{
	if(m_init_done) PlaceElements();
	CDialog::OnMove(x, y);
	this->RedrawWindow(0,0,RDW_INVALIDATE|RDW_UPDATENOW|RDW_ERASE|RDW_FRAME);   // RDW_FRAME redraws the NonClient area too
}

void CPlotToolDlg::OnLButtonUp(UINT nFlags, CPoint point)
{
	this->GetPopupChild()->GetView()->RedrawWindow();
}

BOOL CPlotToolDlg::OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult)
{
	if (message == WM_SIZING)
	{
		this->GetPopupChild()->GetView()->skipdrawduringmove = 1;
	}
	if (message == WM_MOVING)
	{
		this->GetPopupChild()->GetView()->skipdrawduringmove = 1;
	}
	if(message == WM_EXITSIZEMOVE) // catch end of resize or move for additional redraw
	{
		this->GetPopupChild()->GetView()->skipdrawduringmove = 0;
		this->GetPopupChild()->GetView()->RedrawWindow();
	}
	return CDialog::OnWndMsg(message, wParam, lParam, pResult);
}


void CPlotToolDlg::OnClose()
{
	if(_isclosing) // dialog is already closing - prevent cascade
		return;
	_isclosing = 1;
	WriteData();

	// close popup
	if (m_hasPopupChild)
		GetPopupChild()->OnClose();
	m_init_done = 0;
	GetPopupChild()->OnClose();

  // remove from list in WCDriver
	if(!_externclose)
		this->GetWCDDlg() -> DeletePlotToolDlg(this);

	CDialog::DestroyWindow(); // make sure that OnInitDialog is called on next open 
	CDialog::OnClose();
}

void CPlotToolDlg::ExternClose()
{
	_externclose = 1;
	OnClose();
	_externclose = 0;
}

void CPlotToolDlg::OnPaint()
{
	// initialize view
	CDialog::OnPaint();
	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnPaint\n"));
	if(GetPopupChild()) GetPopupChild()->UpdateWindow();
}

// keep Dialog & View alligned
void CPlotToolDlg::PlaceElements()
{
#ifdef preview
	// embedded view 
	CRect rect;
	GetDlgItem(IDC_GRAPHVIEWWND)->GetWindowRect(&rect);
	ScreenToClient(rect);
	GetPreview()->MoveWindow(&rect);
#endif
	if(m_hasPopupChild)
	{
		// allign with child
		CRect parentrect,childrect;
		this->GetWindowRect(&parentrect);
		m_pChild->GetWindowRect(&childrect);

		int vshift = childrect.top-parentrect.top;
		int hshift = childrect.left-parentrect.right;

		childrect.top -= vshift; 
		childrect.bottom -= vshift;
		childrect.left -= hshift;
		childrect.right -= hshift;

		SetCaller(mystr("REDRAW triggered by CPlotToolDlg::PlaceElements - MoveWindow\n"));
	//	m_pChild->MoveWindow(childrect,FALSE); //AD: 2013-07-04 added Flag such that NO redraw is done for move
		m_pChild->MoveWindow(childrect,TRUE); 
	}
}

// **************************
// *** Menu Bar
// **************************
void CPlotToolDlg::OnLayoutsSave()
{
// determine filename - TODO File Save Dialog
	mystr filename;
	CFileDialog fd(FALSE, "plf", "C:\\temp\\layoutEDC.txt", OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, "TXT file (*.txt)|*.txt|PlotTool Layout file (*.plf)|*.plf|All Files (*.*)|*.*||");
	if(fd.DoModal() == IDOK) 
	{
		CString fn = fd.GetPathName();
		filename = mystr(fn);
	}
	SaveLayoutToFile(filename);
}

void CPlotToolDlg::SaveLayoutToFile(mystr& filename)
{
	ElementDataContainer main_edc, settings_edc, lines_edc, palette_edc;
	ElementData ed;
  
	UpdateData(TRUE);
// fill the EDC for properties
	ed.SetBool(m_flag_use_autoupdate,"auto_redraw"); settings_edc.Add(ed);
	ed.SetDouble(m_update_interval,"auto_redraw_interval"); settings_edc.Add(ed);
	ed.SetBool(m_flag_autorescale,"auto_rescale"); settings_edc.Add(ed);
	ed.SetBool(m_flag_draw_every_nth,"flag_draw_every_nth"); settings_edc.Add(ed);
	ed.SetInt(m_draw_every_nth,"draw_every_nth"); settings_edc.Add(ed);
	ed.SetBool(m_flag_mark_every_nth,"flag_mark_every_nth"); settings_edc.Add(ed);
	ed.SetInt(m_mark_every_nth,"mark_every_nth"); settings_edc.Add(ed);
	ed.SetInt(m_flag_statusinfo,"status_bar_info"); settings_edc.Add(ed);

	ed.SetDouble(GetPopupChild()->GetView()->pts_x_left,"distance_left_logical"); settings_edc.TreeAdd("Positioning",ed);

	ed.SetBool(m_axis_at_origin,"draw_at_origin"); settings_edc.TreeAdd("Axis",ed);
	ed.SetBool(m_axis_label_major,"label_major"); settings_edc.TreeAdd("Axis",ed);
	ed.SetDouble(m_axis_overdraw,"overdraw"); settings_edc.TreeAdd("Axis",ed);
	ed.SetDouble(m_axis_ticks_overdraw,"ticksize"); settings_edc.TreeAdd("Axis",ed);
	ed.SetInt(m_digits_x_label,"digits_x_labels"); settings_edc.TreeAdd("Axis",ed);
	ed.SetInt(m_digits_y_label,"digits_y_labels"); settings_edc.TreeAdd("Axis",ed);

	ed.SetBool(m_flag_show_legend,"show"); settings_edc.TreeAdd("Legend",ed);
	ed.SetDouble(pWCDI->GetDOption(275),"left"); settings_edc.TreeAdd("Legend",ed);
	ed.SetDouble(pWCDI->GetDOption(276),"right"); settings_edc.TreeAdd("Legend",ed);
	ed.SetDouble(pWCDI->GetDOption(277),"top"); settings_edc.TreeAdd("Legend",ed);
	ed.SetDouble(pWCDI->GetDOption(278),"bottom"); settings_edc.TreeAdd("Legend",ed);

	ed.SetText(m_export_filename,"filename"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetText(m_export_path,"path"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetInt(m_pixels_horizontal,"size_horizontal"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetInt(m_pixels_vertical,"size_vertical"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetBool(m_flag_export_jpg,"store_jpg"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetBool(m_flag_export_png,"store_png"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetBool(m_flag_export_bmp,"store_bmp"); settings_edc.TreeAdd("SavePictue",ed);
	ed.SetBool(m_flag_export_emf,"store_emf"); settings_edc.TreeAdd("SavePictue",ed);

	ed.SetText(m_title,"title"); settings_edc.TreeAdd("Captions_Labels",ed);
	ed.SetText(m_x_axis_title,"xaxis"); settings_edc.TreeAdd("Captions_Labels",ed);
	ed.SetText(m_y_axis_title,"yaxis"); settings_edc.TreeAdd("Captions_Labels",ed);

	ed.SetDouble(GetPopupChild()->GetView()->shownrange.PMin().X(),"xmin"); settings_edc.TreeAdd("Range",ed);
	ed.SetDouble(GetPopupChild()->GetView()->shownrange.PMax().X(),"xmax"); settings_edc.TreeAdd("Range",ed);
	ed.SetDouble(GetPopupChild()->GetView()->shownrange.PMin().Y(),"ymin"); settings_edc.TreeAdd("Range",ed);
	ed.SetDouble(GetPopupChild()->GetView()->shownrange.PMax().Y(),"ymax"); settings_edc.TreeAdd("Range",ed);

	ed.SetEDC(&settings_edc,"PlotToolOptions");
	main_edc.Add(ed);

// fill the EDC for axis settings & labels

// fill the EDC for lines
	for(int i=1; i<=this->lines.Length(); i++)
	{
		ElementDataContainer new_line_edc;
		lines(i).WriteToEDC(&new_line_edc);
		mystr containername = mystr("Line")+mystr(i);
		ed.SetEDC(&new_line_edc,containername);
		lines_edc.Add(ed);
	}
  ed.SetEDC(&lines_edc,"Data_Sets");
	main_edc.Add(ed);

// add the EDC containing the Palette
	palette.WriteToEDC(&palette_edc);
	ed.SetEDC(&palette_edc,"Color_Palette");
	main_edc.Add(ed);

// save the curent model name for a possible warning
	ed.SetText(pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name"),"modelname");
	main_edc.Add(ed);

// save EDC to File as last entry
	ed.SetText(filename,"File_name");
	main_edc.Add(ed);
	pWCDI->CallCompFunction(111,0,0,&main_edc);
}

void CPlotToolDlg::OnLayoutsLoad()
{
// determine filename 
	mystr filename("C:\\temp\\layoutEDC.txt");
	CFileDialog fd(TRUE, "txt", 0, OFN_FILEMUSTEXIST, "TXT file (*.txt)|*.txt|PlotTool Layout file (*.plf)|*.plf|All Files (*.*)|*.*||");
	if(fd.DoModal() == IDOK) 
	{
		CString fn = fd.GetPathName();
		filename = mystr(fn);
	}
	_isloading = 1;
	LoadLayoutFromFile(filename);
	_isloading = 0;
}

void CPlotToolDlg::LoadLayoutFromFile(mystr& filename)
{
//	Reset();
// Load File, get EDC
	ElementDataContainer main_edc;
	ElementData ed;
	ed.SetText(filename,"File_name"); main_edc.Add(ed);
	pWCDI->CallCompFunction(110,0,0,&main_edc);                           // "Load Text-File to EDC"
  
	int found = 0;
	found = main_edc.Find("modelname");
	if(found)
	{
		mystr mn_file = main_edc.TreeGetString("modelname");
		mystr mn_sele = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name");
		if (!mn_file.Compare(mn_sele))
		{
			pWCDI->GetUserInterface()->AddText("WARNING: Model name in layout file does not match name of selected model !! \n");
		}
	}
	
// parse sub-EDC properties
	found = main_edc.Find("PlotToolOptions");
  if(found && main_edc.Get(found).IsEDC())
	{
		ElementDataContainer* sub_edc = main_edc.Get(found).GetEDC();

// parse options in main EDC
		found = sub_edc->Find("auto_redraw");
		if(found) m_flag_use_autoupdate = sub_edc->TreeGetInt("auto_redraw");
		found = sub_edc->Find("auto_redraw_interval");
		if(found) m_update_interval = sub_edc->TreeGetDouble("auto_redraw_interval");
		found = sub_edc->Find("auto_rescale");
		if(found) m_flag_autorescale = sub_edc->TreeGetInt("auto_rescale");
		found = sub_edc->Find("flag_draw_every_nth");
		if(found) m_flag_draw_every_nth = sub_edc->TreeGetInt("flag_draw_every_nth");
		found = sub_edc->Find("draw_every_nth");
		if(found) m_draw_every_nth = sub_edc->TreeGetInt("draw_every_nth");
		found = sub_edc->Find("flag_mark_every_nth");
		if(found) m_flag_mark_every_nth = sub_edc->TreeGetInt("flag_mark_every_nth");
		found = sub_edc->Find("mark_every_nth");
		if(found) m_mark_every_nth = sub_edc->TreeGetInt("mark_every_nth");
		found = sub_edc->Find("status_bar_info");
		if(found) m_flag_statusinfo = sub_edc->TreeGetInt("status_bar_info");

// parse sub-EDC Positioning e.g. View
// THINK ABOUT: same names as in Options EDC ???
		found = sub_edc->Find("Positioning");
		if(found)
		{
			ElementDataContainer* positioning_edc = sub_edc->Get(found).GetEDC();

			found = positioning_edc->Find("distance_left_logical");
			if(found) 
				GetPopupChild()->GetView()->pts_x_left = positioning_edc->TreeGetDouble("distance_left_logical");
		}

		found = sub_edc->Find("Axis");
		if(found)
		{
			ElementDataContainer* axis_edc = sub_edc->Get(found).GetEDC();

			found = axis_edc->Find("draw_at_origin");
			if(found)	m_axis_at_origin = axis_edc->TreeGetBool("draw_at_origin");
			found = axis_edc->Find("label_major");
			if(found) m_axis_label_major = axis_edc->TreeGetBool("label_major");
			found = axis_edc->Find("overdraw");
			if(found) m_axis_overdraw = axis_edc->TreeGetDouble("overdraw");
			found = axis_edc->Find("ticksize");
			if(found) m_axis_ticks_overdraw = axis_edc->TreeGetDouble("ticksize");
			found = axis_edc->Find("digits_x_labels");
			if(found) m_digits_x_label = axis_edc->TreeGetInt("digits_x_labels");
			found = axis_edc->Find("digits_y_labels");
			if(found) m_digits_y_label = axis_edc->TreeGetInt("digits_y_labels");
		}

// parse sub-EDC Legend
		found = sub_edc->Find("Legend");
		if(found)
		{
			ElementDataContainer* legend_edc = sub_edc->Get(found).GetEDC();

			found = legend_edc->Find("show");
			if(found) m_flag_show_legend = legend_edc->TreeGetInt("show");
			found = legend_edc->Find("left");
			if(found) pWCDI->SetDOption(275,legend_edc->TreeGetInt("left"));
			//if(found) legend.left = (int) (pts_x_plot * legend_edc->TreeGetInt("left") /100. +0.5);			
			found = legend_edc->Find("right");
			if(found) pWCDI->SetDOption(276,legend_edc->TreeGetInt("right"));
			found = legend_edc->Find("top");
			if(found) pWCDI->SetDOption(277,legend_edc->TreeGetInt("top"));
			found = legend_edc->Find("bottom");
			if(found) pWCDI->SetDOption(278,legend_edc->TreeGetInt("bottom"));
		}

		found = sub_edc->Find("SavePicture");
		if(found)
		{
			ElementDataContainer* save_edc = sub_edc->Get(found).GetEDC();
			
			found = save_edc->Find("filename"); 
			if(found) m_export_filename = CString(save_edc->TreeGetString("filename"));
			found = save_edc->Find("path"); 
			if(found) m_export_path = CString(save_edc->TreeGetString("path"));
			found = save_edc->Find("size_horizontal"); 
			if(found) m_pixels_horizontal = save_edc->TreeGetInt("size_horizontal");
			found = save_edc->Find("size_vertical"); 
			if(found) m_pixels_vertical = save_edc->TreeGetInt("size_vertical");

			found = save_edc->Find("store_jpg"); 
			if(found) m_flag_export_jpg = save_edc->TreeGetBool("store_jpg");
			found = save_edc->Find("store_png"); 
			if(found) m_flag_export_png = save_edc->TreeGetBool("store_png");
			found = save_edc->Find("store_bmp"); 
			if(found) m_flag_export_bmp = save_edc->TreeGetBool("store_bmp");
			found = save_edc->Find("store_emf"); 
			if(found) m_flag_export_emf = save_edc->TreeGetBool("store_emf");
		}

// parse sub-EDC Range
		found = sub_edc->Find("Range");
		if(found)
		{
			ElementDataContainer* range_edc = sub_edc->Get(found).GetEDC();
			
			double xmin=0.; double xmax=0.; double ymin=0.; double ymax=0.;
			found = range_edc->Find("xmin");
			if(found) xmin = range_edc->TreeGetDouble("xmin"); 
			found = range_edc->Find("xmax");
			if(found) xmax = range_edc->TreeGetDouble("xmax"); 
			found = range_edc->Find("ymin");
			if(found) ymin = range_edc->TreeGetDouble("ymin"); 
			found = range_edc->Find("ymax");
			if(found) ymax = range_edc->TreeGetDouble("ymax"); 

			GetPopupChild()->GetView()->SetShownRange(xmin,xmax,ymin,ymax);
		}

// parse sub-EDC Caption_Labels
		found = sub_edc->Find("Captions_Labels");
		if(found)
		{
			ElementDataContainer* captions_edc = sub_edc->Get(found).GetEDC();

			found = captions_edc->Find("title");
			if(found) m_title = captions_edc->TreeGetString("title"); 
			found = captions_edc->Find("xaxis");
			if(found) m_x_axis_title = captions_edc->TreeGetString("xaxis"); 
			found = captions_edc->Find("yaxis");
			if(found) m_y_axis_title = captions_edc->TreeGetString("yaxis"); 
		}
	}
// Finished with parsing the properties
	UpdateData(FALSE);

	// remember the syves range in the layout
	double xmin = GetPopupChild()->GetView()->Get_XMin(); 
	double xmax = GetPopupChild()->GetView()->Get_XMax(); 
	double ymin = GetPopupChild()->GetView()->Get_YMin(); 
	double ymax = GetPopupChild()->GetView()->Get_YMax(); 

// parse sub-EDC Lines
	found = main_edc.Find("Data_Sets");
	if(found && main_edc.Get(found).IsEDC())
	{
		ElementDataContainer* sub_edc = main_edc.Get(found).GetEDC();
		lines.Flush();
		for(int i=1; i<= sub_edc->Length(); i++)
		{
			if(sub_edc->Get(i).IsEDC())
			{
				ElementDataContainer* line_edc = sub_edc->Get(i).GetEDC();
				Line_Data_Prop newline;
				newline.ReadFromEDC(line_edc);
				AddLine(newline.x_colnr, newline.y_colnr, newline.IsXY()); // use AddLine-Routine to correctly link data (drawback - no properties are applied)
				lines.Last().ReadFromEDC(line_edc); // read a 2nd time to apply all properties in a oneliner
			}
			// automatic rescale is done by AddLine functiopm
			if( m_flag_autorescale) 
				ScaleToFull(); // automatic rescale of the graph !
		}
	}

	if(!m_flag_autorescale) GetPopupChild()->GetView()->SetShownRange(xmin,xmax,ymin,ymax); // overrule values from lines of not on automatic

	// write the latest valid range into the editboxes
	UpdateDrawRange(GetPopupChild()->GetView()->Get_XMin(), GetPopupChild()->GetView()->Get_XMax(), 
									GetPopupChild()->GetView()->Get_YMin(), GetPopupChild()->GetView()->Get_YMax());
	UpdateData(FALSE);

	// parse sub-EDC Palette 
  found = main_edc.Find("Color_Palette");
  if(found && main_edc.Get(found).IsEDC())
	{
		ElementDataContainer* sub_edc = main_edc.Get(found).GetEDC();
    palette.ReadFromEDC(sub_edc);
		line_color.InitializePalette(palette);
	}

// Update the content of the HotintOptions EDC
	// WriteData();
// Write specific values that were just read in 

	RedrawWindow();
	GetPopupChild()->GetView()->RedrawWindow();
}


	// ***********************
	// *** Group: Data Sources
	// ***********************
void CPlotToolDlg::OnBnClickedRadioSource()
{
	UpdateData(TRUE);
	if (m_radio_source==0) // sensors
	{
		GetDlgItem(IDC_BUTTON_FILEBROWSE)->EnableWindow(FALSE);
		GetDlgItem(IDC_EDIT_FILENAME)->EnableWindow(FALSE);
		WriteSensorNamesToColList();
	}
	else // external file
	{
		GetDlgItem(IDC_BUTTON_FILEBROWSE)->EnableWindow(TRUE);
		GetDlgItem(IDC_EDIT_FILENAME)->EnableWindow(TRUE);
		WriteExternalColumnNamesToColList();
	}
	UpdateData(FALSE);
}

void CPlotToolDlg::OnBnClickedButtonFilebrowse()
{
	CFileDialog fd(TRUE, "dat", 0, OFN_FILEMUSTEXIST, "TXT file (*.txt)|*.txt|DAT file (*.dat)|*.dat|All Files (*.*)|*.*||");
	if(fd.DoModal() == IDOK) 
	{
		CString fn = fd.GetPathName();
		CString fn2 = fd.GetFileName();
		if(filename_extern != fn) // different file
		{
			filename_extern = fd.GetPathName().GetString();
			UpdateData(FALSE);
			LoadExternalToMatrixFile();	
			WriteExternalColumnNamesToColList();
		}
		else
		{
			filename_extern = fd.GetPathName().GetString();
			UpdateData(FALSE);
			WriteExternalColumnNamesToColList();
		}
	}
}

void CPlotToolDlg::OnBnClickedButtonLineAddXY()
{
	UpdateData(TRUE);
	int selected_lines = m_list_columns.GetSelCount();

	if(selected_lines!=2)
	{
		MessageBoxA("Select exactly two datasets for a x-y-line","ERROR",MB_OK);
		return;
	}

	int x_colnr = this->prev_selected_column+1;
	int y_colnr = this->act_selected_column+1;

	AddLine(x_colnr,y_colnr,1);
	ScaleToFull();
	UpdateData(FALSE);
	RedrawWindow();
	GetPopupChild()->GetView()->RedrawWindow();
}


void CPlotToolDlg::AddLine(int x_colnr, int y_colnr, int flag_isxy)
{
	Line_Data_Prop the_line(this);
	the_line.x_colnr = x_colnr;
	the_line.y_colnr = y_colnr;
	the_line.color = this->GetDefaultColor(lines.Length()+1);
	SetLineType(the_line, flag_isxy);
	SetLineName(the_line, x_colnr, y_colnr);
	LinkLineData(the_line);
	// Add to array and ListControl
	this->lines.Add(the_line);
	int nr_lines = m_list_lines.GetItemCount();
	m_list_lines.InsertItem(nr_lines,the_line.name);
	// ScaleToFull(); // automatic rescale of the graph ! //AD: removed 2013-07-04
}

void CPlotToolDlg::SetLineType(Line_Data_Prop& the_line, int flag_isxy)
{
	if(the_line.IsAuxLine()) return;
// determine # columns
	if(flag_isxy)
		the_line.SetType(TLT_xy);
	else
		the_line.SetType(TLT_ty);

// determine source
	if(m_radio_source == 1) // extern
	{
		the_line.filename = mystr((LPCTSTR)filename_extern);
		the_line.SetType( the_line.GetType() | TLT_external );
	}
	else // sensor
	{
		the_line.filename = mystr((LPCTSTR)filename_solution);
		the_line.SetType( the_line.GetType() | TLT_sens );
	}
}

void CPlotToolDlg::SetLineName(Line_Data_Prop& the_line, int x_colnr, int y_colnr)
{
	CString buffer_x;
	m_list_columns.GetText(x_colnr-1, buffer_x);
	CString buffer_y;
	m_list_columns.GetText(y_colnr-1, buffer_y);

	if(the_line.GetType() & TLT_xy)
	{
		the_line.name = mystr("X/Y- ") + (LPCTSTR)buffer_x + ", " + (LPCTSTR)buffer_y;
	}
	else
	{
		the_line.name = mystr("T/Y- ") + (LPCTSTR)buffer_y;
	}
}

void CPlotToolDlg::LinkLineData(Line_Data_Prop& the_line)
{
// for sensor with internal arrays: link to the array
	if(the_line.IsSensor())
	{
		if(the_line.IsTY())
		{
// +++ T/Y Plots +++
			if(the_line.y_colnr!=1) 
			{
				the_line.SetXDataPtr( pWCDI->GetSensorTimesArrayPtr(the_line.y_colnr-1) );     // x-values: sensor time array
				the_line.SetYDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.y_colnr-1) );    // y-values: sensor values array
				the_line.SetTDataPtr( pWCDI->GetSensorTimesArrayPtr(the_line.y_colnr-1) );     // t-values: sensot time array
			}
			else // this does not really make sense since this means drawing T over T but nevertheless
			{
				the_line.SetXDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
				the_line.SetYDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
				the_line.SetTDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
			}
		}
		else
		{
// +++ X/Y Plots +++
			if(the_line.y_colnr!=1) 
			{
				the_line.SetXDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.x_colnr-1) );		// x-values: sensor X values array
				the_line.SetYDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.y_colnr-1) );		// y-values: sensor Y values array
				the_line.SetTDataPtr( pWCDI->GetSensorTimesArrayPtr(the_line.x_colnr-1) );		// t-values: sensor X times array
			}
			else // one of the lines is "Time"
			{
			  if(the_line.x_colnr == 1) // x is time... why didn't chose a Y over T plot in the first place ?
				{
					the_line.SetXDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
					the_line.SetTDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
				}
				else
				{
					the_line.SetXDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.x_colnr-1) );
					the_line.SetTDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.x_colnr-1) );
				}
			 
				if(the_line.y_colnr == 1) // use time as y-values... 
				{
					the_line.SetYDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
					the_line.SetTDataPtr( pWCDI->GetSensorTimesArrayPtr(1) );
				}
				else
				{
					the_line.SetYDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.y_colnr-1) );
					the_line.SetTDataPtr( pWCDI->GetSensorValuesArrayPtr(the_line.y_colnr-1) );
				}
			}
		}
	}

// for data from solution file: link to solution file <-- this will be changed soon
	else if(the_line.IsSolFile())
	{
		the_line.SetXDataPtr(solfile->ColumnPtr(the_line.x_colnr));
		the_line.SetYDataPtr(solfile->ColumnPtr(the_line.y_colnr));
		the_line.SetTDataPtr(solfile->ColumnPtr(1));
	}

// for external file: COPY the data into the line to be able to use multiple files
// make sure to DELETE the arrays when the line is deleted
	else if(the_line.IsExtFile())
	{
		the_line.filename = extfile->GetNameNoPath();
		TArray<double>* xdata = (TArray<double>*) new TArray<double>;
		xdata->CopyFrom(*extfile->ColumnPtr(the_line.x_colnr));
		the_line.SetXDataPtr(xdata);
		TArray<double>* ydata = (TArray<double>*) new TArray<double>;
		ydata->CopyFrom(*extfile->ColumnPtr(the_line.y_colnr));
		the_line.SetYDataPtr(ydata);
		if(the_line.x_colnr==1) the_line.SetTDataPtr(xdata);
		else if (the_line.y_colnr==1) the_line.SetTDataPtr(ydata);
		else
		{
			TArray<double>* tdata = (TArray<double>*) new TArray<double>;
			tdata->CopyFrom(*extfile->ColumnPtr(1));
			the_line.SetTDataPtr(tdata);	
		}
	}
	the_line.UpdateMinMax();
}

int CPlotToolDlg::LoadExternalToMatrixFile()
{
	if (extfile!=NULL) delete extfile;
	extfile = new CMatrixFile(mystr(filename_extern), TFMread);
	if(extfile->IsGood())
		return extfile->ReadAllColumns();
	else return 0;
}

int CPlotToolDlg::LoadSolutionToMatrixFile()
{
	if (solfile!=NULL) delete solfile;
	solfile = new CMatrixFile(mystr(filename_solution), TFMread);
	if(solfile->IsGood())
		return solfile->ReadAllColumns();
	else return 0;
}

int CPlotToolDlg::UpdateSolutionToMatrixFile()
{
//	return solfile->ReadAllColumns();
	return solfile->AppendRead();
}

	// ****************************
	// *** ListBox: Available Lines
	// ****************************
// Adds T/Y plot(s) to the list of drawn lines
void CPlotToolDlg::OnBnClickedButtonLineAdd()
{
	UpdateData(TRUE);
	int selected_lines = m_list_columns.GetSelCount();

	// find column "Time"
	CString buffer;
	int t_colnr = 1;
	for(int i=1; i<=m_list_columns.GetCount(); i++)
	{
		m_list_columns.GetText(i-1,buffer);
		if(buffer == CString("Time"))
			t_colnr = i;
	}

	// get array of selected items
	TArray<int> selected;
	selected.SetLen(selected_lines);
	m_list_columns.GetSelItems(selected_lines, selected.GetDataPtr());

	for(int i=1; i<=selected_lines; i++)
	{
		int y_colnr = selected(i) + 1;  // array is zero based, 
		AddLine(t_colnr,y_colnr,0);
	}
	ScaleToFull();
	UpdateData(FALSE);
	RedrawWindow();
	GetPopupChild()->GetView()->UpdateWindow();
}	

void CPlotToolDlg::OnLbnDblclkListColumns()
{
	OnBnClickedButtonLineAdd();
}

void CPlotToolDlg::OnLbnSelchangeListColumns()
{
	prev_selected_column = act_selected_column;
	act_selected_column = m_list_columns.GetCurSel();
}

void CPlotToolDlg::WriteSensorNamesToColList()
{
	m_list_columns.ResetContent();
	mystr fn = mystr((LPCTSTR)filename_solution);
	mystr str;
	TArrayDynamic<mystr> colnames;

#ifdef sensors_from_solution_file
	if(LoadFile_LoadFile2String(fn,str))
	{
		int comments = LoadFile_CountCommentLines(str);
		int columns = LoadFile_CountColumns(str,comments+1);
		this->LoadFile_GetColumnNames(str,colnames,comments);
	}
	m_label_source = "Sensors (solution file)";

#else // directly from sensors
	ElementDataContainer edc;
	pWCDI->CallCompFunction(203,0,1,&edc); //edc contains sensor names - action 1 gives sensors writeresult property

	if(edc.Length() == 0)
	{
		m_label_source = "Sensors - No Sensors available in Model";
	}
	else
	{
		m_label_source = "Sensors (directly)";
		m_list_columns.InsertString(0,"Time"); // column "Time" is not in sensorlist
		for (int i=1; i <= edc.Length(); i++)
		{
// obsolete code - sensors are now read directly from 
			if( edc.Get(i).GetInt() & 1) // add only if the sensor is written to solution file
			{
				colnames.Add(mystr(edc.Get(i).GetDataName()));
			}
		}
	}
#endif
	for(int i=1; i<=colnames.Length(); i++)
		m_list_columns.InsertString(m_list_columns.GetCount(),colnames(i));
}

void CPlotToolDlg::WriteFileColumnNamesToColList()
{
	m_list_columns.ResetContent();
	mystr fn = mystr((LPCTSTR)filename_extern);
	mystr str;
	TArrayDynamic<mystr> colnames;

	CMatrixFile* the_file;
	if (m_radio_source==0) // sensors
    the_file = solfile;
	else
		the_file = extfile;

	for(int i=1; i<= the_file->NColumns(); i++)
	{
		colnames(i) = the_file->ColumnName(i);
	}
}


void CPlotToolDlg::WriteExternalColumnNamesToColList()
{
	m_list_columns.ResetContent();

	if (extfile!=NULL/* && extfile->IsGood()*/)
	{
		m_label_source = mystr(mystr(extfile->NColumns()) + " From File:" + extfile->GetNameNoPath());
		for(int i=1; i<=extfile->NColumns(); i++)
			m_list_columns.InsertString(m_list_columns.GetCount(),extfile->ColumnName(i));
	}
	UpdateData(FALSE);
}

//bugged...
//void CPlotToolDlg::WriteExternalColumnNamesToColList()
//{
//	m_list_columns.ResetContent();
//	mystr fn = mystr((LPCTSTR)filename_extern);
//	mystr str;
//
//	m_label_source = mystr(mystr(extfile->NColumns()) + " From File:" + extfile->name);
//	for(int i=1; i<=extfile->NColumns(); i++)
//		m_list_columns.InsertString(m_list_columns.GetCount(),extfile->ColumnName(i));
//}

////void CPlotToolDlg::WriteExternalColumnNamesToColList()
////{
////	m_list_columns.ResetContent();
////	mystr fn = mystr((LPCTSTR)filename_extern);
////	mystr str;
////	TArrayDynamic<mystr> colnames;
////
////	if((fn.Length() != 0) && LoadFile_LoadFile2String(fn, str))
////	{
////		int comments = LoadFile_CountCommentLines(str);
////		int columns = LoadFile_CountColumns(str,comments+1);
////		this->LoadFile_GetColumnNames(str,colnames,comments);
////	}
////
////	m_label_source = mystr(mystr(colnames.Length()) + " From File:" + filename_extern);
////	for(int i=1; i<=colnames.Length(); i++)
////		m_list_columns.InsertString(m_list_columns.GetCount(),colnames(i));
////}


	// *************************
	// *** Group: Title and Axis
	// *************************
void CPlotToolDlg::OnBnClickedButtonFont()
{
	SelectFontDialog(m_axisfont);    // font dialog
	AssignOtherFonts();
}

// dialog to select a font
void CPlotToolDlg::SelectFontDialog(LOGFONT& logfont)
{
	logfont.lfHeight /= 10;

	CHOOSEFONT cf;
	ZeroMemory( &cf, sizeof( CHOOSEFONT ));
	cf.lStructSize = sizeof( CHOOSEFONT );
	cf.hwndOwner = this->GetSafeHwnd();
	cf.lpLogFont = &logfont;
	cf.Flags = CF_SCREENFONTS | CF_EFFECTS | CF_TTONLY | CF_USESTYLE;

	int result = ChooseFont( &cf );

	logfont.lfHeight *= 10;
}

// Assign all fonts from chosen axisfont (scaled)
void CPlotToolDlg::AssignOtherFonts()
{
	int height = m_axisfont.lfHeight;

	m_titlefont = m_axisfont;
	m_titlefont.lfHeight = (int) (height * m_titlefont_scalingfactor + 0.5);

	m_ticksfont = m_axisfont;
	m_ticksfont.lfHeight = (int) (height * m_ticksfont_scalingfactor + 0.5);
}

// apply changes made in the textboxes to view 
void CPlotToolDlg::EditBoxSetsRange()
{
	UpdateData(TRUE);
	if(this->m_pChild)
	{
		SetCaller(mystr("REDRAW triggered by CPlotToolDlg::EditBoxSetsRange\n"));
		m_pChild->GetView()->SetShownRange(x_min, x_max, y_min, y_max);
		m_pChild->GetView()->RedrawWindow();
	}
}

// write actual limits of graph to textbox
void CPlotToolDlg::UpdateDrawRange(double xmin, double xmax, double ymin, double ymax)
{
	x_min = xmin;
	x_max = xmax;
	y_min = ymin;
	y_max = ymax;
	//if(this->IsWindowVisible()) // otherwise crashes in AssertValid()
	//	UpdateData(FALSE); 
}

// write artual ticksize to textbox
void CPlotToolDlg::UpdateMajorIntervals(double xmajor, double ymajor)
{
	x_major_interval = xmajor;
	y_major_interval = ymajor;
	//	UpdateData(TRUE);
}

	// ************************
	// *** ListBox: Drawn Lines
	// ************************
void CPlotToolDlg::OnLvnItemchangedListLines(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMLISTVIEW pNMLV = reinterpret_cast<LPNMLISTVIEW>(pNMHDR);

	POSITION pos = m_list_lines.GetFirstSelectedItemPosition();
	if(!pos)
	{
		WriteToDlgPropertiesofLine(0);
	}
	else
	{
		WriteToDlgPropertiesofLine(m_list_lines.GetNextSelectedItem(pos)+1);
	}
	UpdateData(FALSE);

	*pResult = 0;
}

void CPlotToolDlg::OnLvnKeydownListLines(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMLVKEYDOWN pLVKeyDow = reinterpret_cast<LPNMLVKEYDOWN>(pNMHDR);
	if (pLVKeyDow->wVKey == VK_DELETE) 
	{
		TArray<int> list_of_lines_to_delete(0);

		POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
		while(pos)
		{
			list_of_lines_to_delete.Add((int)pos);
			m_list_lines.GetNextSelectedItem(pos);
		}
		// as the list of selected line is sorted in ascending order...

		while(list_of_lines_to_delete.Length())
		{
			int n = list_of_lines_to_delete.Last();
			m_list_lines.DeleteItem(n-1);																				// ListCtrl is 0-based
			lines.Erase(n);																											// TArray is 1-based
			list_of_lines_to_delete.Erase(list_of_lines_to_delete.Length());
		}

		Redraw();
	}
	*pResult = 0;
}

void CPlotToolDlg::WriteToDlgPropertiesofLine(int selected_line)
{
	m_selected_line = selected_line;
	if(selected_line == 0) // no line selected
	{
		m_edit_name_of_selection = "";
		m_label_lineprop = mystr("no line selected");
		line_color_index = -1;
		line_thickness_index = -1;
		line_style_index = -1;
		point_style_index = -1;

		GetDlgItem(IDC_EDIT_LINE_NAME)->EnableWindow(FALSE);
		GetDlgItem(IDC_COMBO_LINE_COLOR)->EnableWindow(FALSE);
		GetDlgItem(IDC_COMBO_LINE_THICKNESS)->EnableWindow(FALSE);
		GetDlgItem(IDC_COMBO_LINE_STYLELINE)->EnableWindow(FALSE);
		GetDlgItem(IDC_COMBO_LINE_STYLEPOINT)->EnableWindow(FALSE);
		GetDlgItem(IDC_EDIT_POINTSIZE)->EnableWindow(FALSE);
	}
	else
	{
		Line_Data_Prop& the_line = lines(selected_line);
		
		m_edit_name_of_selection = the_line.name;
		
		COLORREF linecolor = the_line.color;
		line_color_index = palette.Find((int&)linecolor) -1 ;
				
		line_thickness_index = the_line.width -1; 
		line_thickness.activecolor = linecolor;

		line_style_index = the_line.penstyle;               // no correction since defined by WIN as 0-based
		line_style.activecolor = linecolor;

		point_style_index = the_line.pointstyle -1;
		point_style.activecolor = linecolor;

		point_size = lines(selected_line).pointsize;

		m_label_lineprop = mystr("Properties of Line ") + mystr(selected_line) + mystr(": ") + mystr(lines(selected_line).name.Left(16) + " ...");
		GetDlgItem(IDC_EDIT_LINE_NAME)->EnableWindow(TRUE);
		GetDlgItem(IDC_COMBO_LINE_COLOR)->EnableWindow(TRUE);
		GetDlgItem(IDC_COMBO_LINE_THICKNESS)->EnableWindow(TRUE);
		GetDlgItem(IDC_COMBO_LINE_STYLELINE)->EnableWindow(TRUE);
		GetDlgItem(IDC_COMBO_LINE_STYLEPOINT)->EnableWindow(TRUE);
		GetDlgItem(IDC_EDIT_POINTSIZE)->EnableWindow(TRUE);
	}
}

void CPlotToolDlg::OnEnKillfocusEditLineName()
{
	UpdateData(TRUE);
	lines(m_selected_line).name = mystr((LPCTSTR)m_edit_name_of_selection);
	m_list_lines.SetItemText(m_selected_line-1,0,lines(m_selected_line).name);
	WriteToDlgPropertiesofLine(m_selected_line);
	UpdateData(FALSE);
	return;
}

void CPlotToolDlg::OnCbnSelendokComboLineColor()
{
	UpdateData(TRUE);
	int selection = line_color_index;
	COLORREF colref = (COLORREF) line_color.GetItemData(selection);
	
	//this->lines(m_selected_line).color = colref;
	POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
	while(pos)
	{
		this->lines((int)pos).color = colref;
		m_list_lines.GetNextSelectedItem(pos);
	}

// change colors in other pickers
	line_thickness.activecolor = colref;
	line_thickness.RedrawWindow();
	line_style.activecolor = colref;
	line_style.RedrawWindow();
	point_style.activecolor = colref;
	point_style.RedrawWindow();
	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnCbnSelendokComboLineColor\n"));
	Redraw();
}

void CPlotToolDlg::OnCbnSelendokComboLineThickness()
{
	UpdateData(TRUE);
	int selection = line_thickness_index;
	int thickness = (int) line_thickness.GetItemData(selection);
	
	//this->lines(m_selected_line).width = thickness;
	POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
	while(pos)
	{
		this->lines((int)pos).width = thickness;
		m_list_lines.GetNextSelectedItem(pos);
	}

	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnCbnSelendokComboLineThickness\n"));
	Redraw();
}

void CPlotToolDlg::OnCbnSelendokComboLineStyle()
{
	UpdateData(TRUE);
	int selection = line_style_index;
	int penstyle = (int) line_style.GetItemData(selection);
	
	//this->lines(m_selected_line).penstyle = penstyle;
	POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
	while(pos)
	{
		this->lines((int)pos).penstyle = penstyle;
		m_list_lines.GetNextSelectedItem(pos);
	}

	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnCbnSelendokComboLineStyle\n"));
	Redraw();
}

void CPlotToolDlg::OnCbnSelendokComboPointStyle()
{
	UpdateData(TRUE);
	int selection = point_style_index;
	int pointstyle = (int) point_style.GetItemData(selection);
	
	//this->lines(m_selected_line).pointstyle = pointstyle;
	POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
	while(pos)
	{
		this->lines((int)pos).pointstyle = pointstyle;
		m_list_lines.GetNextSelectedItem(pos);
	}

	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnCbnSelendokComboPointStyle\n"));
	Redraw();
}

void CPlotToolDlg::OnEnKillfocusEditPointsize()
{
	UpdateData(TRUE);

	//this->lines(m_selected_line).pointsize = point_size;
	POSITION pos = this->m_list_lines.GetFirstSelectedItemPosition(); 
	while(pos)
	{
		this->lines((int)pos).pointsize = point_size;
		m_list_lines.GetNextSelectedItem(pos);
	}

	SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnEnKillfocusEditPointsize\n"));
	Redraw();
}

// **************************
// *** Group: General Options
// **************************
void CPlotToolDlg::SetUpdateRate(double interval_s)
{
	m_flag_use_autoupdate = 1;
	m_update_interval = interval_s;
	UpdateData(FALSE);
	StartStopAutoUpdate();
}

void CPlotToolDlg::StartStopAutoUpdate()
{
	UpdateData(TRUE);
	if (m_flag_use_autoupdate) 
	{
		//GetDlgItem(IDC_EDIT_UPDATE)->EnableWindow(TRUE);
		StartRefreshTimer();
	}
	else 
	{
		//GetDlgItem(IDC_EDIT_UPDATE)->EnableWindow(FALSE);
		KillRefreshTimer();
	}
}

void CPlotToolDlg::OnEnKillfocusEditUpdate()
{
	UpdateData(TRUE);
	StartStopAutoUpdate();
}


	// *****************
	// *** Group: Export
	// *****************
void CPlotToolDlg::ExportToFile()
{
	UpdateData(TRUE);

	if (GetPopupChild()->GetView())
	{
		CString path_and_filename_noext = m_export_path + mystr("\\") + m_export_filename;
		if(m_flag_export_jpg || m_flag_export_png || m_flag_export_bmp)
		{
			GetPopupChild()->GetView()->SaveAsBitmapGraphic(path_and_filename_noext, m_pixels_horizontal, m_pixels_vertical, m_flag_export_jpg, m_flag_export_png, m_flag_export_bmp);
		}
		if(m_flag_export_emf)
		{
			GetPopupChild()->GetView()->SaveAsVectorGraphic(path_and_filename_noext, m_flag_export_emf);
		}
	}
	else
		this->MessageBox("No View to export","ERROR:",IDOK);
}

	// ***************************** 
	// *** Group: Additional Buttons 
	// *****************************  
void CPlotToolDlg::ShowGraph()
{
	if(m_hasPopupChild)
	{
		this->m_pChild->ShowWindow(SW_SHOW);
	}
	PlaceElements();
}

void CPlotToolDlg::ScaleToFull()
{
	if (GetPopupChild()->GetView() != NULL)
	{
		GetPopupChild()->GetView()->ZoomToFullRange(0,1);
	}
}

void CPlotToolDlg::Redraw()
{
	if (GetPopupChild()->GetView() != NULL)
	{
		if(m_flag_autorescale)
		{
			GetPopupChild()->GetView()->ZoomToFullRange(0,1);
		}
		else
		{
			GetPopupChild()->GetView()->RedrawWindow(); // popup view can become illegal
		}
	}
#ifdef preview
	if (GetPreview())
	{
		GetPreview()->RedrawWindow(); // nested view is always valid
	}
#endif
}

void CPlotToolDlg::UpdateLineData()
{
	//LoadData();
	for(int i=1; i<=lines.Length(); i++)
	{
		if (lines(i).IsSensor()) lines(i).ReadLine();
	}
	Redraw();
}

  // **********************
	// *** calls from MenuBar
  // **********************
// initializes the dialog with a single t/y sensor already entered in the lines list
void CPlotToolDlg::InitializeWithSensorTY(int sensor_number, CRect viewrect)
{
	ShowWindow(SW_HIDE);
	// data line
	m_radio_source = 0;  // select solution
	WriteSensorNamesToColList();

// select only chosen line in dialog list
	for(int i=1; i<= m_list_columns.GetCount(); i++)
		m_list_columns.SetSel(i-1,0);
	m_list_columns.SetSel(sensor_number,1); // -1 for [0]-based,  +1 for "time" data

// add line and scale axis
	OnBnClickedButtonLineAdd();
	ScaleToFull();

// enter name in window, set as filename
	CString buffer;
	m_list_columns.GetText(sensor_number,buffer);

	mystr fn = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.record_frames_path");
	m_export_filename = CString(fn)+buffer;
	
	//buffer = "sensor direct - " + buffer;
	buffer = "HOTINT Plot Tool";
	GetPopupChild()->SetDefaultName(buffer);
	GetPopupChild()->SetDisplayName(buffer);

// place & set update
	if(viewrect != CRect(0,0,0,0))
		GetPopupChild()->PlaceView(viewrect);
	SetUpdateRate(m_update_interval); 
	StartStopAutoUpdate();
}

// initializes the dialog with a single x/y sensor already entered in the lines list
void CPlotToolDlg::InitializeWithSensorXY(int sensor_number_x, int sensor_number_y, CRect viewrect)
{
	ShowWindow(SW_HIDE);
// data line
	m_radio_source = 0;  // select solution
	WriteSensorNamesToColList();

	for(int i=1; i<= m_list_columns.GetCount(); i++)
		m_list_columns.SetSel(i-1,0);
	m_list_columns.SetSel(sensor_number_x,1); // -1 for [0]-based,  +1 for "time" data
	OnLbnSelchangeListColumns();              // to correctly remember the last 2 selected in order
	m_list_columns.SetSel(sensor_number_y,1); // -1 for [0]-based,  +1 for "time" data
	OnLbnSelchangeListColumns();

	OnBnClickedButtonLineAddXY();
	ScaleToFull();

	CString buffer,bufferx,buffery;
	m_list_columns.GetText(sensor_number_x,bufferx);
	m_list_columns.GetText(sensor_number_y,buffery);

	mystr fn = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.record_frames_path");
	m_export_filename = CString(fn) + buffery + "(" + bufferx + ")";

	//buffer = "sensor direct - " + buffery + "(" + bufferx + ")";
	buffer = "HOTINT Plot Tool";
	GetPopupChild()->SetDefaultName(buffer);
	GetPopupChild()->SetDisplayName(buffer);

	if(viewrect != CRect(0,0,0,0))
		GetPopupChild()->PlaceView(viewrect);
	SetUpdateRate(m_update_interval);
	StartStopAutoUpdate();
}

// initializes N dialogs with a single t/y sensor each
void CPlotToolDlg::InitializeWithNSensors(TArray<int>& sensor_nrs, TArray<CRect>& viewrects)
{
 
}


#ifdef zombies
	// ***************************************************************** 
	// *** GRAVEYARD - collect outdated code - disabled control elements
	// *****************************************************************
	// **************************
	// *** Group: Auxillary Lines
	// **************************
void CPlotToolDlg::OnBnClickedButtonAddAuxLine()
{
	UpdateData(TRUE);

	//// split into function name and 
	//CString left = m_auxline_function.Left(2);
	//CString parsestring = m_auxline_function;
	//parsestring.Delete(0,2);

	//Line_Data_Prop line(this);
	////line.InitAsMathFunction();

	//mystr expression;
	//mystr variable;
	////CParser& parser = (CParser&) pWCDI->MBSParser();
	//CMBSParser& parser = pWCDI->MBSParser();

	//if (left.CompareNoCase(CString("y=")) == 0) // !!! returns 0 for identical strings
	//{
	//	//line.p_mf->SetExpression((LPCTSTR)parsestring, "x", parser);
	//	//line.p_pf->SetParsedFunction1D(&parser, (LPCTSTR)parsestring, "x");
	//}
	//else if (left.CompareNoCase(CString("x=")) == 0) // !!! returns 0 for identical strings
	//{
	//	//line.p_mf->SetExpression((LPCTSTR)parsestring, "y", parser);
	//	//line.p_pf->SetParsedFunction1D(&parser, (LPCTSTR)parsestring, "y");
	//}
	//else
	//{
	//	;
	//}

	//this->lines.Add(line);
	//int nr_lines = m_list_lines.GetItemCount();
	//m_list_lines.InsertItem(nr_lines,line.name);

	ScaleToFull(); // automatic rescale of the graph !
}

void CPlotToolDlg::AddAuxLine(int flag_isxy)
{
  Line_Data_Prop the_line(this);
	the_line.SetType(flag_isxy);
	the_line.SetType( the_line.GetType() | TLT_aux );
}

	// ***************
	// *** TAB CONTROL 
	// ***************
void CPlotToolDlg::OnTcnSelchangeTabPlottool(NMHDR *pNMHDR, LRESULT *pResult)
{
	DrawActiveTabElements();

	*pResult = 0;
}

void CPlotToolDlg::DrawActiveTabElements()
{
	CTabCtrl* TC_TabCtrl = (CTabCtrl*) GetDlgItem(IDC_TAB_PLOTTOOL);
	int sel = TC_TabCtrl->GetCurSel();
	int tabnumber = 0;

	// sel goes from 0 .. 4

  // elements that are active on Tab "Title" (, Axis, Font)
	tabnumber = 0;
	((CStatic*) GetDlgItem(IDC_STATIC_LABEL_TITLE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_LABEL_TITLE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_LABEL_XAXIS)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_LABEL_XAXIS)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_LABEL_YAXIS)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_LABEL_YAXIS)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);

	((CButton*) GetDlgItem(IDC_BUTTON_FONT)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);

  // elements that are active on Tab "General Options"
  tabnumber = 1;
	((CButton*) GetDlgItem(IDC_CHECK_UPDATE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_UPDATE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_CHECK_UPDATE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_AUTORESCALE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_LEGEND)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_DRAW_EVERY_NTH)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_DRAW_EVERY_NTH)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_CHECK_DRAW_EVERY_NTH)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_VERT_AT_CURRENT)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_DRAW_UPTO_CURRENT)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_DRAW_FAST)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);

  // elements that are active on Tab "Export Options"
  tabnumber = 2;
	((CStatic*) GetDlgItem(IDC_STATIC_EXPORT_FILEPATH)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_EXPORT_FILEPATH)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_BUTTON_PATHBROWSE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_EXPORT_FILENAME)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_EXPORT_FILENAME)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_EXPORT_RESOLUTION)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_RESX)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CStatic*) GetDlgItem(IDC_STATIC_EXPORT_X)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CEdit*) GetDlgItem(IDC_EDIT_RESY)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_BUTTON_FIGURESIZEASSCREEN)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_EXPORT_JPG_SIZE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_EXPORT_PNG_SIZE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_EXPORT_BMP_SIZE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_CHECK_EXPORT_EMF_SIZE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);


	// elements that are active on Tab "Buttons"
  tabnumber = 3;
	((CButton*) GetDlgItem(IDC_BUTTON_SHOW)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_AXISEQUAL)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_BUTTON_HIDE)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_BUTTON_REDRAW)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDC_READOPTIONS)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);
	((CButton*) GetDlgItem(IDPRINT)) -> ShowWindow(sel==tabnumber ? SW_SHOW : SW_HIDE);

}


void CPlotToolDlg::OnBnClickedButtonAxisEqual()
{
	if (GetPopupChild()->GetView() != NULL)
	{
		GetPopupChild()->GetView()->ZoomToFullRange(1,1);
	}
}

void CPlotToolDlg::OnBnClickedButtonExport()
{
	ExportToFile();
}

// hide the Main dialog
void CPlotToolDlg::OnBnClickedButtonHide()
{
	this->ShowWindow(SW_HIDE);
}



// loads entire file to a string
int CPlotToolDlg::LoadFile_LoadFile2String(mystr& filename, mystr& str)
{
	CMFile file(filename, (TFileMode)TFMread);

	if (!file.IsGood())
	{
		
		MessageBox((mystr("Could not open File")+filename).c_str(),"ERROR",MB_OK);
		return 0;
	}
	str = "";
	file.RWF(str);

	return str.Length();
}

// count number of comment lines in the file
int CPlotToolDlg::LoadFile_CountCommentLines(mystr& filestring)
{
	int pos=0;
	int linecount=0;
	mystr line;

	while (pos!=-1)
	{
		filestring.GetUntil(pos,'\n',line,1);
		if(line.Left(1).Compare(mystr("%")))
			linecount++;
		else 
			pos=-1;
	}
	return linecount;
}

// counts number of columns in first data line
int CPlotToolDlg::LoadFile_CountColumns(mystr& filestring, int linenr)
{
	if(!linenr)
		linenr = LoadFile_CountCommentLines(filestring) + 1;

	int pos=0;
	mystr line;
	for(int i=1; i<=linenr; i++)
		filestring.GetUntil(pos,'\n',line,1);

	int pos_line=0;
	int numbercount=0;
	mystr number;
	while (pos_line!=-1)
	{
		number = line.GetWord(pos_line,1);
		if(pos_line!=-1)
			numbercount++;
	}
	return numbercount;
}

// extracts column names from last commented line
int CPlotToolDlg::LoadFile_GetColumnNames(mystr& filestring, TArrayDynamic<mystr>& columnnames, int linenr)
{
	columnnames.Flush();

	if(!linenr) 
		linenr = LoadFile_CountCommentLines(filestring);
	if(!linenr) 
	{
		// no line containing column names
		int cols = LoadFile_CountColumns(filestring,linenr);
		for(int i=1; i<= cols; i++)
			columnnames.Add(mystr("Column_")+mystr(i));
	}
	else
	{
		int pos=0;
		mystr line;
		for(int i=1; i<=linenr; i++)
			filestring.GetUntil(pos,'\n',line,1);

		line = line.Right(line.Length()-1); // remove "%"

		int pos_line=0;
		mystr word;
		while (pos_line!=-1)
		{
			word = line.GetWord(pos_line,1);
			if(pos_line!=-1)
				columnnames.Add(word);
		}
	}
	return columnnames.Length();
}

// reads a column to array
int CPlotToolDlg::LoadFile_ReadColumnToArray(int colnr, TArray<double>* column, mystr& filestring, int& append_at_line)
{
	if ((column->Length()+1) > append_at_line)
	{
		// data content of the column is too long - maybe computation restarted?
		append_at_line = 0;
		column->Flush();
	}

	if ((column->Length()+1) < append_at_line) // 
	{
		// data content of the columne is too short - ??? - missed a reload ?
		append_at_line = 0;
		column->Flush();
	}

	int rv;
	int pos=0;
	mystr line;
	int line_count = 0; // count data lines

	while (pos!=-1)
	{
		rv = filestring.GetUntil(pos,'\n',line,1);

		if(!line.Left(1).Compare(mystr("%")) && pos!=-1)
		{
			line_count++;
			if(line_count < append_at_line) continue; // skip all lines already read
			int pos_line = 0;
			if(line.Length() == 1) pos_line=-1; // skip last line.... 
			while (pos_line!=-1)
			{
				mystr word;
				for (int i = 1; i<= colnr; i++)
				{
					word = line.GetWord(pos_line,1);
					if (pos_line == -1)
					{
						MessageBox("column not found","PARSE ERROR",MB_OK);
						return 0;
					}
				}

				double nr = word.MakeDouble();
				column->Add(nr);
				pos_line = -1; //continue; 
			}			
		}
	}
	return column->Length();
}

// identify as Hotint - solution file
int CPlotToolDlg::LoadFile_IsSolFile(mystr& filestring)
{
	mystr line;
	int pos=0;
	int	rv = filestring.GetUntil(pos,'\n',line,1);
	if (line.Compare("%HOTINT_Solution_File"))
		return 1;
	return 0;
}

// identify as Hotint - solution parameter file
int CPlotToolDlg::LoadFile_IsSolParFile(mystr& filestring)
{
	mystr line;
	int pos=0;
	int	rv = filestring.GetUntil(pos,'\n',line,1);
	if (line.Compare("%HOTINT_Parameter_Solution_File"))
		return 1;
	return 0;
}
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyPlotToolView:       derived View class                                                                                  * 29.11.2012 +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLEMENT_DYNCREATE(CMyPlotToolView, CView)

CMyPlotToolView::CMyPlotToolView(CWnd* parent)
{
}

void CMyPlotToolView::ErrorMessage(mystr& msg)
{
	this->GetPlotToolDlg()->GetWCDI()->GetUserInterface()->AddText( mystr("Illegal pDC - ") + caller );
}

void CMyPlotToolView::Reset()
{
	LoadData();
	prevdragrect = CRect(0,0,0,0);
}

BEGIN_MESSAGE_MAP(CMyPlotToolView, CView)
	ON_WM_MOUSEWHEEL()
	ON_WM_MOUSEMOVE()
	ON_WM_NCMOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
END_MESSAGE_MAP()

// *********************************
// COMMUNICATION WITH SETTINGs (EDC)
// *********************************
void CMyPlotToolView::LoadData()
{
	WCDInterface* pWCDI = GetPlotToolDlg()->GetWCDI();

	// logical size of plot
	pts_x_plot = pWCDI->GetIOption(262);
	pts_y_plot = pWCDI->GetIOption(263);

	m_global_linethickness_factor = (int)pWCDI->GetDOption(258); 
	// logical size of lables
	pts_y_top = (int) (pts_y_plot*pWCDI->GetDOption(261)/100. + 0.5);
	pts_x_left = (int) (pts_x_plot*pWCDI->GetDOption(260)/100. + 0.5);
	pts_y_bottom = (int) (pts_y_plot*pWCDI->GetDOption(262)/100. + 0.5);
	pts_x_right = (int) (pts_y_plot*pWCDI->GetDOption(263)/100. + 0.5);

	// grid lines
	grid_x_major_linetype = pWCDI->GetIOption(270);
	grid_x_minor_linetype = pWCDI->GetIOption(271);
	grid_y_major_linetype = pWCDI->GetIOption(272);
	grid_y_minor_linetype = pWCDI->GetIOption(273);

	// legend
	legend.left = (int) (pts_x_plot*pWCDI->GetDOption(275)/100. + 0.5);
	legend.right = (int) (pts_x_plot*pWCDI->GetDOption(276)/100. + 0.5);
	legend.top = (int) (pts_y_plot*pWCDI->GetDOption(277)/100. + 0.5);
	legend.bottom = (int) (pts_y_plot*pWCDI->GetDOption(278)/100. + 0.5);
}

// write changeable settings back to EDC
void CMyPlotToolView::WriteData()
{
}

// ***********************************************
// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// ***********************************************
BOOL CMyPlotToolView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)    // Additional funcitonality: zoom only in one diection
{
	// class specific 
	// uses default zooming (with SHIFT and CTRL activated)
	UINT nFlags_specific = nFlags;

	BOOL rv = CMyBase2DView::OnMouseWheel(nFlags_specific,zDelta,pt); // this already sets Shownrange in CMyBase2DView::ReComputeShownRange_Zooming

	GetPlotToolDlg()->UpdateDrawRange(Get_XMin(), Get_XMax(), Get_YMin(), Get_YMax());
	
	SetCaller(mystr("REDRAW triggered by CMyPlotToolView::OnMouseWheel\n"));
		
	return RedrawWindow();
}

void CMyPlotToolView::OnLButtonDown(UINT nFlags, CPoint point)
{
	DragDropStart() = point;
	CView::OnLButtonDown(nFlags, point);
}

void CMyPlotToolView::OnLButtonUp(UINT nFlags, CPoint point)
{
	AllowRedraw();                                            // No Redraw of lines during move
	DragDropFinal() = point;

	// determine if rectangle zoom is used
	if(nFlags & MK_SHIFT)
	{
		RectZoomAction();
	}
	else
	{
	// standard drag drop - drag drop end
		DragDropAction();
		SetCaller(mystr("REDRAW triggered by CMyPlotToolView::OnLButtonUp\n"));
		RedrawWindow();
	}
	CView::OnLButtonUp(nFlags, point);
}

void CMyPlotToolView::OnMouseMove(UINT nFlags, CPoint point)
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
			SetCaller(mystr("REDRAW triggered by CMyPlotToolView::OnMouseMove-DragDrop\n"));
			RedrawWindow();
			DragDropStart() = point;
		}
	}

	if ((nFlags & MK_LBUTTON) && (nFlags & MK_SHIFT)) // left mouse button is down, shift key is down --> rectzoom --> animate
	{
		//// check if distance for redraw is reached
		//double moved_x = DragDropFinal().x - DragDropStart().x;
		//double moved_y = DragDropFinal().y - DragDropStart().y;
		//double dragdropdistance = sqrt(moved_x*moved_x + moved_y*moved_y);
		//
		//if (dragdropdistance > 10.)
		{
				DragDropFinal() = point;
				DrawZoomRect();
		}
	}

// status line entry for mouse position
	CRect rect;
	this->GetClientRect(&rect);
	
	if( rect.PtInRect(point) )
	{
		CPoint mousepos = GetLogicalFromPixels(point, mystr("CMyPlotToolView::OnMouseMove"));

		MousePosition() = mousepos;
		double x = XofLogical(point);
		double y = YofLogical(point);
	
		this->SetStatusBarText_X(mystr("X: ")+mystr(x)); 
		this->SetStatusBarText_Y(mystr("Y: ")+mystr(y)); 
		UpdateStatusBarInfo();

	}
	else
	{
		// ASSERT: mouse Pos is IN Client when MouseMove is triggered  thus
		// THIS CODE IS NEVER REACHED ... 
		// place this code in OnNcMouseMove
	}
	CView::OnMouseMove(nFlags, point);
}

void CMyPlotToolView::UpdateStatusBarInfo()
{
	// Status Bar Info is set in DrawLines()....
}

void CMyPlotToolView::OnNcMouseMove(UINT nHitTest, CPoint point)
{
	// this is reached when mouse moves OFF the client rect
	this->SetStatusBarText_X(mystr("")); 
	this->SetStatusBarText_Y(mystr("")); 
	CMyBase2DView::OnNcMouseMove(nHitTest, point);
	// DOES NOT TRIGGER WHEN MOVEMENT IS TOO FAST 
	// other approaches:
	// WM_MOUSELEAVE - in ::OnWndMsg
	// TrackMouseEvent - in ::OnMouseMove
}

BOOL CMyPlotToolView::OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult)
{
	if (message == WM_SIZING)
	{
		skipdrawduringmove = 1;
	}
	if (message == WM_MOVING)
	{
		skipdrawduringmove = 1;
	}
	if (message == WM_EXITSIZEMOVE)
	{
		skipdrawduringmove = 0;
		RedrawWindow();
	}
	if (message == WM_MOUSELEAVE)
	{
		this->SetStatusBarText_X(mystr("")); 
		this->SetStatusBarText_Y(mystr("")); 
		return CMyBase2DView::OnWndMsg(message, wParam, lParam, pResult);
	}
	return CMyBase2DView::OnWndMsg(message, wParam, lParam, pResult);
}


void CMyPlotToolView::DragDropAction()
{
	CMyBase2DView::DragDropAction(); 
	GetPlotToolDlg()->UpdateDrawRange(Get_XMin(), Get_XMax(), Get_YMin(), Get_YMax());
}

void CMyPlotToolView::RectZoomAction()
{
	CPoint logical_start = GetLogicalFromPixels(DragDropStart(),mystr("RectZoomAction - StartPoint"));
	CPoint logical_final = GetLogicalFromPixels(DragDropFinal(),mystr("RectZoomAction - FinalPoint"));
	CRect rect(logical_start,logical_final); // zoom rectangle in logical coordinates - just for min/max
	
	double xmin = XofLogical(rect.TopLeft());
	double xmax = XofLogical(rect.BottomRight());
	double ymin = YofLogical(rect.TopLeft());
	double ymax = YofLogical(rect.BottomRight());

	// ranges
	shownrange = Box2D(Vector2D(xmin, ymin),Vector2D(xmax, ymax));
	GetPlotToolDlg()->UpdateDrawRange(Get_XMin(), Get_XMax(), Get_YMin(), Get_YMax());

	prevdragrect = CRect(0,0,0,0); 

	SetCaller(mystr("REDRAW triggered by CMyPlotToolView::RectZoomAction\n"));
	RedrawWindow();
}

void CMyPlotToolView::DrawZoomRect()
{
	CPoint logical_start = GetLogicalFromPixels(DragDropStart(),mystr("RectZoomAction - StartPoint"));
	CPoint logical_final = GetLogicalFromPixels(DragDropFinal(),mystr("RectZoomAction - FinalPoint"));
	CRect dragrect(logical_start,logical_final); // zoom rectangle in logical coordinates 
	
	CDC* pDC = this->GetDC(); 
	if(!pDC)                           // bailout in case the device context is not valid
	{
		this->GetPlotToolDlg()->GetWCDI()->GetUserInterface()->AddText("Illegal pDC - CMyPlotToolView::DrawZoomRect");
		return;
	}
	OnPrepareDC(pDC);                  // mapping function

	// draw surrounding Rect
	CPen pen_zoomrect (PS_SOLID, 1 ,RGB(0,0,0));
	pDC->SelectObject(&pen_zoomrect);
	pDC->MoveTo(dragrect.left,dragrect.bottom);
	pDC->LineTo(dragrect.right,dragrect.bottom);
	pDC->LineTo(dragrect.right,dragrect.top);
	pDC->LineTo(dragrect.left,dragrect.top);
	pDC->LineTo(dragrect.left,dragrect.bottom);
	
	ReleaseDC(pDC);
}


// *******************
// THE DRAWING PROCESS
// *******************
void CMyPlotToolView::OnDraw(CDC* pDC)
{
// TRY TO CATCH MANY REDRAW MESSAGES
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_DRAWITEM, WM_DRAWITEM, PM_REMOVE));

// is redraw currently legal ?
	if(DoNotRedraw()) 
		return;
	ForbidRedraw();

// reset caller string for any calls by system
	SetCaller(mystr("REDRAW triggered by SYSTEM\n"));

	if(0)     //(GetPlotToolDlg()->m_flag_fastdraw == 0) 
	{
		DrawElements(pDC); // direct draw
	}
	else // 
	{
		CMyMemDC MDC(pDC);
		DrawElements(&MDC); // fast draw



////// offscreen draw
////		CRect L_ViewRect;	
////		this->GetViewRect(L_ViewRect);
////		CRect the_offscreenrect(0,0,abs(L_ViewRect.Width()),abs(L_ViewRect.Height()));
////
////// Memory DC where the image is drawn
////		CDC memDC;
////		memDC.CreateCompatibleDC(NULL);
////// set MapMode for device context
////		PrepareDC_inRect(&memDC, NULL, the_offscreenrect);
////		memDC.m_hAttribDC = memDC.m_hDC;
////		memDC.SetMapMode(pDC->GetMapMode());
////	
////// CBitmap in Memory DC to draw into
////		CBitmap* p_mem_Bitmap; 
////		p_mem_Bitmap = new CBitmap();
////// create bitmap in logical coordinates
////		p_mem_Bitmap->CreateBitmap(the_offscreenrect.Width(), the_offscreenrect.Height(), 1, memDC.GetDeviceCaps(BITSPIXEL), NULL);
////// initialize BMP for drawing
////		memDC.SelectObject(p_mem_Bitmap);
////// white backgroud color for entire bitmap region ( include the region that is not covered due to MM_ISOMETRIC )
////		memDC.PatBlt(L_ViewRect.left, L_ViewRect.bottom, L_ViewRect.Width(), -L_ViewRect.Height(), WHITENESS);
////// draw into DC bitmap
////		DrawElements(&memDC);
////
////// copy to real dc
////		pDC->BitBlt(L_ViewRect.left, L_ViewRect.bottom, abs(L_ViewRect.Width()), abs(L_ViewRect.Height()), &memDC, L_ViewRect.left, L_ViewRect.bottom, SRCCOPY);
////
////		delete p_mem_Bitmap;
	}
	AllowRedraw();
}

void CMyPlotToolView::DrawElements(CDC* pDC)
{
	if(shownrange.Empty())
	{
		shownrange = Box2D(Vector2D(0.,0.),Vector2D(1.,1.));
//		return;
	}

//// Get Mouse Button Status:  
//  BOOL LeftMouseButtonDown = ( (GetKeyState(VK_LBUTTON) & 0x80) != 0 );

	CPlotToolDlg* p_PTD = GetPlotToolDlg();
// BEGIN OF DRAWING
	// draw axis
	this->DrawAxis(pDC);

	// compute the ticks and labels
	this->DrawTicksAndLabels(pDC);

	// plot title
	this->DrawTitle(pDC);

	// axis labels
	this->DrawAxisLabels(pDC);

// Do not redraw the lines when window is moved/resized - variable skipdrawduringmove is changed in the function OnWndMsg
	if(!skipdrawduringmove)
	{
		// draw lines
		this->DrawLines(pDC);
	}

	// draw legend	
	if (p_PTD->m_flag_show_legend)
		this->DrawLegend(pDC);

// END OF DRAWING
}

void CMyPlotToolView::DrawAxis(CDC* pDC)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();
	CFont font_axislabel;
	font_axislabel.CreateFontIndirectA(&(p_PTD->m_axisfont));
	CRect plotrect;
	this->GetPlotRect(plotrect);
	CRect xaxisrect;
	this->GetXAxisRect(xaxisrect);
	CRect yaxisrect;
	this->GetYAxisRect(yaxisrect);

	// logical coordinates of point (0/0)
	CPoint origin = ValuesToLogical(0,0);
	// where to draw axes
	BOOL flag_x_axis_is_in_plotrect = p_PTD->m_axis_at_origin && (origin.y <= plotrect.top) && (origin.y >= plotrect.bottom); // bottom otherwise
	BOOL flag_y_axis_is_in_plotrect = p_PTD->m_axis_at_origin && (origin.x <= plotrect.right) && (origin.x >= plotrect.left); // left otherwise

	// draw Rect surrounding the graph
	CPen pen_border (PS_SOLID,p_PTD->m_border_thickness,RGB(0,0,0));
	pDC->SelectObject(&pen_border);
	pDC->MoveTo(plotrect.left,plotrect.bottom);
	pDC->LineTo(plotrect.right,plotrect.bottom);
	pDC->LineTo(plotrect.right,plotrect.top);
	pDC->LineTo(plotrect.left,plotrect.top);
	pDC->LineTo(plotrect.left,plotrect.bottom);

	// draw axis
	int overdraw ;
	CPen pen_axis(PS_SOLID,(int)(p_PTD->m_border_thickness*1.5),RGB(0,0,0));
	pDC->SelectObject(&pen_axis);

	overdraw = (int) ((plotrect.right-plotrect.left) * p_PTD->m_axis_overdraw / 100. + 0.5);
	CPoint p1(plotrect.left - overdraw, origin.y * flag_x_axis_is_in_plotrect);   // axis on lower border
	CPoint p2(plotrect.right + overdraw, origin.y * flag_x_axis_is_in_plotrect);
	pDC->MoveTo(p1);
	pDC->LineTo(p2);

	overdraw = (int) ((plotrect.top -plotrect.bottom) * p_PTD->m_axis_overdraw / 100. + 0.5);
	CPoint p3(origin.x * flag_y_axis_is_in_plotrect, plotrect.bottom - overdraw);  // axis on left border
	CPoint p4(origin.x * flag_y_axis_is_in_plotrect, plotrect.top + overdraw);
	pDC->MoveTo(p3);
	pDC->LineTo(p4);
}

void CMyPlotToolView::DrawTicksAndLabels(CDC* pDC)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();
	char buffer[32];

	// define all pens and fonts	
	CPen pen_axis(PS_SOLID,(int)(p_PTD->m_border_thickness*1.5),RGB(0,0,0));
	pDC->SelectObject(&pen_axis);

	CFont font_ticklabel;
	font_ticklabel.CreateFontIndirectA(&(p_PTD->m_ticksfont));
	pDC->SelectObject(&font_ticklabel);

	CRect plotrect;
	GetPlotRect(plotrect);
	// logical coordinates of point (0/0)
	CPoint origin = ValuesToLogical(0,0);
	// where to draw axes
	BOOL flag_x_axis_is_in_plotrect = p_PTD->m_axis_at_origin && (origin.y <= plotrect.top) && (origin.y >= plotrect.bottom); // bottom otherwise
	BOOL flag_y_axis_is_in_plotrect = p_PTD->m_axis_at_origin && (origin.x <= plotrect.right) && (origin.x >= plotrect.left); // left otherwise

// X AXIS
	CRect xaxisrect;
	GetXAxisRect(xaxisrect);
	int x_digits = GetPlotToolDlg()->m_digits_x_label;

// compute ticks, check if strings fit
	double x_tick_scaling_factor = 1.;
	TArray<double> xticks_major;
	CRect default_xaxis_label_rect;

	ComputeMajorTicksArray(xticks_major, Get_XMin(), Get_XMax(), 1, x_tick_scaling_factor);

// adjust the number of ticks such that the labels with a constant fontsize do not overlap
	if(p_PTD->m_axis_label_major)
	{
		int flag_x_ticks_scaled = 0;
		while(!flag_x_ticks_scaled)
		{
			// available space for single entry 
			double tick_abs = xticks_major(2) - xticks_major(1);
			long tick_logical = ValuesToLogical(xticks_major(2),0.).x - ValuesToLogical(xticks_major(1),0.).x;
			double filling_factor = 0.95;              // use up to 95% of available space --> min distance of 5% between texts
			long avaiable_logical = (long) (tick_logical * filling_factor);

			long text_max_x = 0;
			// test all label strings
			for(int i=1; i <= xticks_major.Length(); i++)
			{
				// textbox in non rotated system
				sprintf(buffer,"%.*g",x_digits,xticks_major(i));
				CString str(buffer);
				CSize textsize = pDC->GetTextExtent(str);
				if( textsize.cx > text_max_x)
				{
					text_max_x = textsize.cx;
				}
			}
			if(text_max_x > avaiable_logical)
			{
				// rescaling factor
				x_tick_scaling_factor = x_tick_scaling_factor * ( (double) avaiable_logical / (double) text_max_x );
				ComputeMajorTicksArray(xticks_major, Get_XMin(), Get_XMax(), 1, x_tick_scaling_factor);
			}
			else 
			{
				flag_x_ticks_scaled = 1;
				// compute default rectangle
				CPoint shift(0,-10); // hardcoded shift

				default_xaxis_label_rect.left   = (long) (-0.5 * filling_factor * tick_logical) + shift.x;
				default_xaxis_label_rect.right  = (long) ( 0.5 * filling_factor * tick_logical) + shift.x;

				default_xaxis_label_rect.top    = xaxisrect.top    + shift.y;
				default_xaxis_label_rect.bottom = xaxisrect.bottom + shift.y;
			}
		}
	}
	// scaling done - draw now

	// draw ticks
	long ticksize = (long) ((plotrect.right-plotrect.left) * p_PTD->m_axis_ticks_overdraw / 100. + 0.5);  // drawsize of the tick
	for(int i=1; i <= xticks_major.Length(); i++)
	{
		CPoint point_at_axis = ValuesToLogical(xticks_major(i), 0.);
		point_at_axis.y *= flag_x_axis_is_in_plotrect; 
		pDC->MoveTo( point_at_axis.x, point_at_axis.y + ticksize / 2 );
		pDC->LineTo( point_at_axis.x, point_at_axis.y - ticksize / 2 );
	}

	// write tick labels
	if(p_PTD->m_axis_label_major)
	{
		for(int i=1; i <= xticks_major.Length(); i++)
		{
			sprintf(buffer,"%.*g",x_digits,xticks_major(i)); // number of significant digits x_digits

			CRect labelrect;
			CPoint point_at_axis = ValuesToLogical(xticks_major(i), 0.);
			labelrect.left   = default_xaxis_label_rect.left   + point_at_axis.x;
			labelrect.right  = default_xaxis_label_rect.right  + point_at_axis.x;
			labelrect.top    = default_xaxis_label_rect.top;
			labelrect.bottom = default_xaxis_label_rect.bottom;
			
			this->DrawTextInRect(pDC, CString(buffer), labelrect, 0., TTextAllign(VTop+HCenter));
		}
	}

	// draw gridlines
	if(grid_x_major_linetype > 0)
	{
		DrawVertLines(pDC, xticks_major, grid_x_major_linetype);
	}

// Y AXIS
	CRect yaxisrect;
	GetYAxisRect(yaxisrect);
	int y_digits = GetPlotToolDlg()->m_digits_y_label;

	// compute ticks, check if strings fit
	double y_tick_scaling_factor = 1.;
	TArray<double> yticks_major;
	CRect default_yaxis_label_rect;

	ComputeMajorTicksArray(yticks_major, Get_YMin(), Get_YMax(), 1, y_tick_scaling_factor);

// adjust the number of ticks such that the labels with a constant fontsize do not overlap
// T.B.D. properly
			// available space for single entry 
			double tick_abs = yticks_major(2) - yticks_major(1);
			long tick_logical = ValuesToLogical(0.,yticks_major(2)).y - ValuesToLogical(0.,yticks_major(1)).y;
			double filling_factor = 0.95;              // use up to 95% of available space --> min distance of 5% between texts
			long avaiable_logical = (long) (tick_logical * filling_factor);

	CPoint shift(-40,20); // hardcoded shift

	default_yaxis_label_rect.left   = yaxisrect.left   + shift.x; 
	default_yaxis_label_rect.right  = yaxisrect.right  + shift.x;

	default_yaxis_label_rect.top    = (long) ( 0.5 * filling_factor * tick_logical) + shift.y;
	default_yaxis_label_rect.bottom = (long) (-0.5 * filling_factor * tick_logical) + shift.y;

	// draw ticks
	ticksize = (long) (abs(plotrect.top-plotrect.bottom) * p_PTD->m_axis_ticks_overdraw / 100. + 0.5);
	for(int i=1; i <= yticks_major.Length(); i++)
	{
		CPoint point_at_axis = ValuesToLogical(0., yticks_major(i));
		point_at_axis.x *= flag_y_axis_is_in_plotrect;
		pDC->MoveTo( point_at_axis.x + ticksize / 2, point_at_axis.y );
		pDC->LineTo( point_at_axis.x - ticksize / 2, point_at_axis.y );
	}

	// write tick labels
	if(p_PTD->m_axis_label_major)
	{
		for(int i=1; i <= yticks_major.Length(); i++)
		{
			sprintf(buffer,"%.*g",y_digits,yticks_major(i)); // number of significant digits y_digits

			CRect labelrect;
			CPoint point_at_axis = ValuesToLogical(0., yticks_major(i));
			labelrect.left   = default_yaxis_label_rect.left;
			labelrect.right  = default_yaxis_label_rect.right;
			labelrect.top    = default_yaxis_label_rect.top    + point_at_axis.y;
			labelrect.bottom = default_yaxis_label_rect.bottom + point_at_axis.y;
			
			this->DrawTextInRect(pDC, CString(buffer), labelrect, 0., TTextAllign(VCenter+HRight));
		}
	}

	// draw gridlines
	if(grid_y_major_linetype > 0)
	{
		DrawHoriLines(pDC, yticks_major, grid_y_major_linetype);
	}
}

void CMyPlotToolView::DrawVertLines(CDC* pDC, TArray<double> x, int line_type)
{
	CRect plotrect;
	this->GetPlotRect(plotrect);

	CPen pen_solid(PS_SOLID, 1, RGB(32,32,32));
	CPen pen_dash(PS_DASH, 1, RGB(32,32,32));
	CPen pen_dot(PS_DOT, 1, RGB(32,32,32));

	if(line_type == 1) 
		pDC->SelectObject(&pen_solid);
	else if(line_type == 2) 
		pDC->SelectObject(&pen_dash);
	else if(line_type == 3)
		pDC->SelectObject(&pen_solid);

	CPoint point_at_axis;
	for (int i=1; i <= x.Length(); i++)
	{
		point_at_axis = ValuesToLogical(x(i), 0.);
		pDC->MoveTo( point_at_axis.x, plotrect.bottom );
		pDC->LineTo( point_at_axis.x, plotrect.top );	
	}
}

void CMyPlotToolView::DrawHoriLines(CDC* pDC, TArray<double> y, int line_type)
{
	CRect plotrect;
	this->GetPlotRect(plotrect);

	CPen pen_solid(PS_SOLID, 1, RGB(32,32,32));
	CPen pen_dash(PS_DASH, 1, RGB(32,32,32));
	CPen pen_dot(PS_DOT, 1, RGB(32,32,32));

	if(line_type == 1) 
		pDC->SelectObject(&pen_solid);
	else if(line_type == 2) 
		pDC->SelectObject(&pen_dash);
	else if(line_type == 3)
		pDC->SelectObject(&pen_solid);

	CPoint point_at_axis;
	for (int i=1; i <= y.Length(); i++)
	{
		point_at_axis = ValuesToLogical(0.,y(i));
		pDC->MoveTo( plotrect.left, point_at_axis.y);
		pDC->LineTo( plotrect.right, point_at_axis.y);
	}
}

void CMyPlotToolView::DrawTitle(CDC* pDC)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();

	CFont font_title;
	font_title.CreateFontIndirectA(&(p_PTD->m_titlefont));
	CRect titlerect;
	this->GetTitleRect(titlerect);

	// plot title
	pDC->SelectObject(&font_title);
	DrawTextInRect(pDC,p_PTD->m_title,titlerect,0.0,TTextAllign(HCenter+VCenter));
}

void CMyPlotToolView::DrawAxisLabels(CDC* pDC)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();

	CFont font_x_axislabel;
	font_x_axislabel.CreateFontIndirectA(&(p_PTD->m_axisfont));
	CFont font_y_axislabel;
	font_y_axislabel.CreateFontIndirectA(&(p_PTD->m_axisfont));

	// axis labels
	CRect xaxisrect;
	this->GetXAxisRect(xaxisrect);
	pDC->SelectObject(&font_x_axislabel);
	DrawTextInRect(pDC,p_PTD->m_x_axis_title,xaxisrect,0.0,TTextAllign(HCenter+VCenter));

	CRect yaxisrect;
	this->GetYAxisRect(yaxisrect);
	pDC->SelectObject(&font_y_axislabel);
	DrawTextInRect(pDC,p_PTD->m_y_axis_title,yaxisrect,90.0,TTextAllign(HLeft+VCenter));
}

void CMyPlotToolView::DrawLines(CDC* pDC)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();
	CPen pen_plot (PS_SOLID,GetPenWidth(PT_FINE),RGB(255,0,0));
	CRect plotrect;
	this->GetPlotRect(plotrect);                                  // rectangle available for drawing

	int nr_lines = p_PTD->lines.Length();                         // number of lines to be drawn
	double tdraw = p_PTD->GetWCDI()->GetActualDrawTime();               // "current" time from datamanager (same as OpenGL-Window)

	mystr buffer("   ");
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//draw lines of each data set
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int last_point_to_be_drawn = 0;                                 // cutoff if only draw up to "current"
	for (int i=1; i<=nr_lines; i++)
	{
		CPoint linestartpoint; 
		CPoint lineendpoint;
		CPoint currentdatapoint;
		CPoint previousdatapoint;
		CRect markerrect;

// line properties from array
		COLORREF linecolor = p_PTD->lines(i).color;
		int linethickness = (int)(GetPenWidth(TPenWidth(p_PTD->lines(i).width))*m_global_linethickness_factor);
		int linestyle = p_PTD->lines(i).penstyle;
		int markerstyle = p_PTD->lines(i).pointstyle;
		int markersize = p_PTD->lines(i).pointsize;
		markerrect = CRect(-markersize, +markersize, +markersize, -markersize);

// define pens
		CPen pen_line (linestyle, linethickness, linecolor);	
		CPen pen_point (PS_SOLID, linethickness, linecolor);	
		pDC->SelectObject(&pen_line);

// concerning dataset
		TArray<double>& x_val = *(p_PTD->lines(i).x_data);
		TArray<double>& y_val = *(p_PTD->lines(i).y_data);
		int is_xy = p_PTD->lines(i).IsXY();
		int length_of_dataset = x_val.Length();
		if (length_of_dataset == 0) continue;
		if (y_val.Length()<length_of_dataset) length_of_dataset = y_val.Length();

// identify the position of the "current" time -
		int data_point_number_of_current = 0; // 0 for "not computed yet"
		last_point_to_be_drawn = length_of_dataset;  
		// does the entry number of the "current time" (from DataManager) have to be computed ?
		if (p_PTD->m_flag_upto_current || p_PTD->m_flag_vertical_marker)
		{
//!AD: THIS IS OPEN FOR DISCUSSION ( once or for every line - could become important if the time discretization of the data is not the same for all datalines )
//!AD: HERE WE ASSUME THAT TIME AXIS OF FIRST LINE COUNTS FOR ALL 			
			if(data_point_number_of_current<1) // enter only once
			{
				TArray<double>& t_val = *(p_PTD->lines(1).t_data);
			// assume sorted array for times, find last point to be shown
				double acc = 1.e-10;
				if ( t_val.Length()>1 ) acc = (t_val(1)-t_val(2))*1e-6; // accuracy for comparison of doubles
				if (t_val.Last()>tdraw) 
				{
					while ( (t_val(data_point_number_of_current) <= tdraw+acc) && (data_point_number_of_current<=t_val.Length()) )
						data_point_number_of_current++;
				}
				else data_point_number_of_current = last_point_to_be_drawn;
				if (p_PTD->m_flag_upto_current && data_point_number_of_current < last_point_to_be_drawn ) last_point_to_be_drawn = data_point_number_of_current;
			}
		}

// ****************
// DRAWING OF LINES
// ****************
// new procedure as of 03-07-2013
// 

// Specific Marker for the current time
		if(GetPlotToolDlg()->m_flag_vertical_marker)
		{
			currentdatapoint = ValuesToLogical(x_val(data_point_number_of_current),y_val(data_point_number_of_current));
			DrawMarkerCurrent(pDC, currentdatapoint, GetPenWidth(PT_THICK));
		}

// first point
		linestartpoint = ValuesToLogical(x_val(1),y_val(1));
		pDC->MoveTo( linestartpoint );
		pDC->SelectObject(&pen_point);
		this->DrawDataPoint( pDC, markerstyle, linestartpoint, markerrect, linecolor);          // always mark first data point
		pDC->SelectObject(&pen_line);

		int counter_skipdraw = p_PTD->m_draw_every_nth;
		int counter_skipmark = p_PTD->m_mark_every_nth;		
		
		previousdatapoint = linestartpoint;
		
// all the lines plus respective end points
		for (int j=2; j<=last_point_to_be_drawn; j++)
		{
			counter_skipdraw--;
			counter_skipmark--;

			currentdatapoint = ValuesToLogical(x_val(j),y_val(j));   // position of the next datapoint

//EXCEPTION 1:
// detect a Y(t) line going back to origin - can happen when using parameter variation
			if (!is_xy && (currentdatapoint.x < previousdatapoint.x) ) 
			{
// always finish previous line ( this line ends at the "previousdatapoint" )
				pDC->LineTo(previousdatapoint);
// always mark last point of previous line
				pDC->SelectObject(&pen_point);
				this->DrawDataPoint( pDC, markerstyle, previousdatapoint, markerrect, linecolor);          
				pDC->SelectObject(&pen_line);

// start a new line at the "currentdatapoint"
				pDC->MoveTo( currentdatapoint );
				linestartpoint = currentdatapoint;
// always mark first point of new line
				pDC->SelectObject(&pen_point);
				this->DrawDataPoint( pDC, markerstyle, currentdatapoint, markerrect, linecolor);          
				pDC->SelectObject(&pen_line);

// let the skipping counters start over
		    counter_skipdraw = p_PTD->m_draw_every_nth;
		    counter_skipmark = p_PTD->m_mark_every_nth;		
			}

			else
			{
// DATA POINT IS TO BE MARKED (in principle)
				if(!p_PTD->m_flag_mark_every_nth || !counter_skipmark) // mark a datapoint when skipping is off or counter is at 0
				{
					counter_skipmark = p_PTD->m_mark_every_nth;

					if(IsPointInRect(currentdatapoint,plotrect))
					{
						pDC->SelectObject(&pen_point);
						this->DrawDataPoint( pDC, markerstyle, currentdatapoint, markerrect, linecolor);          
						pDC->SelectObject(&pen_line);
					}
				}

// LINE IS TO BE DRAWN (in principle)
				if(!p_PTD->m_flag_draw_every_nth || !counter_skipdraw)  // draw a line when skipping is off or counter is at 0
				{
					counter_skipdraw = p_PTD->m_draw_every_nth;
					lineendpoint = currentdatapoint;
				
					if (IsPointInRect(linestartpoint,plotrect) && IsPointInRect(lineendpoint,plotrect))
					{
								pDC->LineTo( lineendpoint ); // entire line is visible - faster draw
					}
					else
					{
// EXCEPTION 2: line is partially out of range ( goes through bounding rectangle or is entirely invisible )
						DrawLineVisible(pDC, plotrect, linestartpoint, lineendpoint); // line only partly visible - slower draw routine
					}
					linestartpoint = lineendpoint;
				}
			}
		} // end loop(datapoints)
				
// entries in Status line
		if(i==1) 
			buffer = mystr(", y1= ") + mystr(y_val(last_point_to_be_drawn));
		else
			buffer += mystr(", y") + mystr(i) + mystr("= ") + mystr(y_val(last_point_to_be_drawn),4);
	} // end loop(datasets)

// x-value
	if(nr_lines)
	{
		if(p_PTD->lines(1).x_data->Length() >= last_point_to_be_drawn && p_PTD->lines(1).x_data->Length() > 0 && last_point_to_be_drawn > 0)
		{
			buffer = mystr("X= ") + mystr(p_PTD->lines(1).x_data->Get(last_point_to_be_drawn),4) + buffer;
		}
		else
		{
			buffer = mystr("X= undef") + buffer;
		}
	}
	this->SetStatusBarText_Main(buffer);
}

void CMyPlotToolView::DrawDataPoint(CDC* pDC, int pointstyle, CPoint& centerpoint, CRect& markerrect, COLORREF color)
{
	COLORREF crOldBkColor = pDC->GetBkColor();
	CBrush whitebrush(crOldBkColor);
 	CBrush colorbrush(color);
	
	pDC->SelectStockObject(HOLLOW_BRUSH);
	CPoint prevpoint = pDC->MoveTo(centerpoint);

	switch (pointstyle)
	{
	case PTS_NONE: 
		break;
	case PTS_X: 
		pDC->MoveTo(centerpoint.x + markerrect.left , centerpoint.y + markerrect.bottom);
		pDC->LineTo(centerpoint.x + markerrect.right, centerpoint.y + markerrect.top   );
		pDC->MoveTo(centerpoint.x + markerrect.left,  centerpoint.y + markerrect.top   );
		pDC->LineTo(centerpoint.x + markerrect.right, centerpoint.y + markerrect.bottom);
		break;
	case PTS_CROSS:
		pDC->MoveTo(centerpoint.x + markerrect.left , centerpoint.y                    );
		pDC->LineTo(centerpoint.x + markerrect.right, centerpoint.y                    );
		pDC->MoveTo(centerpoint.x,                    centerpoint.y + markerrect.top   );
		pDC->LineTo(centerpoint.x,                    centerpoint.y + markerrect.bottom);
		break;
	case PTS_DOT: 
 		pDC->SelectObject(&colorbrush);
	  pDC->Ellipse(centerpoint.x + markerrect.left, centerpoint.y + markerrect.top, centerpoint.x + markerrect.right, centerpoint.y + markerrect.bottom);
//		pDC->SelectObject(&whitebrush);
		pDC->SelectStockObject(HOLLOW_BRUSH);
		break;
	case PTS_CIRCLE: 
    pDC->Ellipse(centerpoint.x + markerrect.left, centerpoint.y + markerrect.top, centerpoint.x + markerrect.right, centerpoint.y + markerrect.bottom);
  	break;
	default:
		break;
	}
	pDC->MoveTo(prevpoint);
}

void CMyPlotToolView::DrawMarkerCurrent(CDC* pDC, CPoint datapoint, int thickness)
{
	////CRect plotrect;
	////this->GetPlotRect(plotrect);
	CPen pen_solid(PS_SOLID, thickness, RGB(32,0,0));
	CPen* prevpen = pDC->SelectObject(&pen_solid);

	////CPoint point_at_axis = ValuesToLogical(t, 0.);
	////pDC->MoveTo( point_at_axis.x, plotrect.bottom-100 );
	////pDC->LineTo( point_at_axis.x, plotrect.top+100 );	

	//COLORREF crOldBkColor = pDC->GetBkColor();
	//CBrush whitebrush(crOldBkColor);
 //	CBrush colorbrush(color);
	
	//pDC->SelectStockObject(HOLLOW_BRUSH);

		pDC->MoveTo(datapoint.x, datapoint.y + thickness*5);
		pDC->LineTo(datapoint.x, datapoint.y );
 
    pDC->LineTo(datapoint.x - thickness*1, datapoint.y + thickness*1);
		pDC->MoveTo(datapoint.x + thickness*1, datapoint.y + thickness*1);
		pDC->LineTo(datapoint.x , datapoint.y );

		pDC->SelectObject(prevpen);
}


// draws section of the line that is IN the plot-rectangle when startpoint or endpoint are NOT
void CMyPlotToolView::DrawLineVisible(CDC* pDC, CRect& plotrect, CPoint& start, CPoint& end, double oversize_percentage)
{
	CPoint startpoint,endpoint;
	int flag_draw=0;

	CPoint delta = end - start;

	if( (delta.x == 0) && (delta.y == 0) ) // never draw a line that has identical start & end points
	{
		flag_draw = 0;
	}

	else if( IsPointInRect(start,plotrect) && IsPointInRect(end,plotrect) ) // always draw line when both points are in rectangle
	{
		flag_draw=1;
		startpoint = start;
		endpoint = end;
	}

	else
	{
		// calculate 4 possible points at border
		CPoint points[4];
		// initialize intersection points outside of border
		for(int i=0; i<4; i++)
		{
			points[i] = CPoint(plotrect.left-1, plotrect.bottom-1); // make sure point is outside rect
		}

		// intersections with borders

		////if(delta.y == 0) // horizontal line
		////{
		////	points[0] = CPoint(plotrect.left, start.y); // point at x=0
		////	points[1] = CPoint(plotrect.right, start.y); // point at x=max
		////}
		////else if(delta.x == 0) // vertical line
		////{
		////	points[2] = CPoint(start.x, plotrect.bottom); // point at y=0
		////	points[3] = CPoint(start.x, plotrect.top); // point at y=max
		////}
		////else // find all true intersections of the line to draw with the border( lambda factor between zero and 1 )
		{
			double l_x0 = (plotrect.left - start.x)/(double)delta.x;
			if( (0. <= l_x0) && (l_x0 <= 1.) )
				//points[0] = CPoint(start.x+(int)(l_x0*delta.x+0.5), start.y+(int)(l_x0*delta.y+0.5)); // point at x=0
				points[0] = CPoint(0, start.y+(int)(l_x0*delta.y+0.5)); // point at x=0


			double l_xm = (plotrect.right - start.x)/(double)delta.x;
			if( (0. <= l_xm) && (l_xm <= 1.) )
				//points[1] = CPoint(start.x+(int)(l_xm*delta.x+0.5), start.y+(int)(l_xm*delta.y+0.5)); // point at x=max
				points[1] = CPoint(plotrect.right, start.y+(int)(l_xm*delta.y+0.5)); // point at x=max

			double l_y0 = (plotrect.bottom - start.y)/(double)delta.y;
			if( (0. <= l_y0) && (l_y0 <= 1.) )
				//points[2] = CPoint(start.x+(int)(l_y0*delta.x+0.5), start.y+(int)(l_y0*delta.y+0.5)); // point at y=0
				points[2] = CPoint(start.x+(int)(l_y0*delta.x+0.5), 0); // point at y=0

			double l_ym = (plotrect.top - start.y)/(double)delta.y;
			if( (0. <= l_ym) && (l_ym <= 1.) )
				//points[3] = CPoint(start.x+(int)(l_ym*delta.x+0.5), start.y+(int)(l_ym*delta.y+0.5)); // point at y=max
				points[3] = CPoint(start.x+(int)(l_ym*delta.x+0.5), plotrect.top); // point at y=max
		}

		// decide which of the 4 points to use
		int isin[4] = {0,0,0,0};
		for (int i=0; i<4; i++)
		{
			isin[i] = IsPointInRect(points[i],plotrect);
		}
		int arein_x = isin[0]+isin[1];
		int arein_y = isin[2]+isin[3];
		int arein = arein_x+arein_y;

		int i = 0;
		switch(arein)
		{
		case 0: // line is not visible - can only happen if: a) both points are out 
			break; 

		case 1: // can happen if: a) one point is within rect or at rect border, line goes outside (startpoint or endpoint)
			// strategy: set start AND end to intersection point, then replace with point that is inside rect
			for(i=0; i<=3; i++)
			{
				if(isin[i])
				{
					startpoint = points[i];
					endpoint = points[i];
				}
			}
			if( IsPointInRect(start,plotrect))
			{
				startpoint = start;
			}
			if( IsPointInRect(end,plotrect))
			{
				endpoint = end;
			}
			flag_draw = 1;
			break; 

		case 2: // can happen if: a) line has visible segment with both points outside rect or on border of rect
			//                b) one point is inside and line goes through corner
			// strategy: set start and end point form two different intersection points, then replace with point inside rect
			while(!isin[i] && i<4) i++;
			startpoint = points[i];
			i++;
			while(!isin[i] && i<4) i++;
			endpoint = points[i];
			if( IsPointInRect(start,plotrect))
			{
				startpoint = start;
			}
			if( IsPointInRect(end,plotrect))
			{
				endpoint = end;
			}
			flag_draw = 1;
			break; 

		case 3: // can happen if: a) two points coincide (1 corner) 
			// strategy: pick points on opposite sides of the rectangle, check delta-vector such that start and end point are placed in correct order
			// (then replace with point inside rect)
			if(arein_x ==2)
			{
				if (delta.x >= 0) // line goes s-->e
				{
					startpoint = points[0];
					endpoint = points[1];
				}
				else // line goes e<--s
				{
					startpoint = points[1];
					endpoint = points[0];
				}
			}
			else
			{
				if (delta.y >= 0) // line goes ^
				{
					startpoint = points[2];
					endpoint = points[3];
				}
				else
				{
					startpoint = points[3];
					endpoint = points[2];
				}
			}
			flag_draw = 1;
			break;

		case 4: // can happen if: a) two points coincide (2 corners) 
			// strategy: pick points on left and right borders, check delta-vector such that start and end point are placed in correct order
			// (then replace with point inside rect)
			if (delta.x >= 0) // line goes s-->e
			{
				startpoint = points[0];
				endpoint = points[1];
			}
			else // line goes e<--s
			{
				startpoint = points[1];
				endpoint = points[0];
			}
			flag_draw = 1;
			break; 

		default:
			break;
		} // switch
	} // else

	if(flag_draw)
	{
		pDC->MoveTo( startpoint );
		pDC->LineTo( endpoint );
	}
	pDC->MoveTo( end ); // starting point for next line
}

void CMyPlotToolView::DrawLegend(CDC* pDC)
{
	// +++ todo 
	// legend placement, Text,...
	CPlotToolDlg* p_PTD = GetPlotToolDlg();
	CPen pen_axis (PS_SOLID,GetPenWidth(PT_NORMAL),RGB(0,0,0));

	CRect legendrect;
	this->GetLegendRect(legendrect);

	// draw legend	
	pDC->FillSolidRect(legendrect, 0xFFFFFF);
	pDC->SelectObject(&pen_axis);
	pDC->MoveTo(legendrect.left,legendrect.top);
	pDC->LineTo(legendrect.right,legendrect.top);
	pDC->LineTo(legendrect.right,legendrect.bottom);
	pDC->LineTo(legendrect.left,legendrect.bottom);
	pDC->LineTo(legendrect.left,legendrect.top);

	long nr_lines = p_PTD->lines.Length();      // number of lines in plot

	if(nr_lines < 1)  
		return;

	long nr_slots = nr_lines;										// number of slots in legened
	long y_slot_height = (long) (legendrect.Height()/(nr_slots)); // height for a single slot in the legend
	double filling_factor = 0.95;              // use up to 95% of available space --> min distance of 5% between texts
	long available_logical = (long) abs(y_slot_height * filling_factor);

// adjust font size to available y - space
	double scaling_factor = 1.;

	long text_max_y = 0;
	for(int i=1; i<=nr_lines; i++)
	{
		CSize textsize = pDC->GetTextExtent(CString(p_PTD->lines(i).name));
		if(textsize.cy > text_max_y)
		{
			text_max_y = textsize.cy;
		}
	}
	if(text_max_y > available_logical)
	{
		scaling_factor = (double) available_logical / (double) text_max_y ;

		LOGFONT logfont;
		// font must be scaled down to fit target rectangle
		pDC->GetCurrentFont()->GetLogFont(&logfont);
		logfont.lfHeight = (long) (logfont.lfHeight*scaling_factor);
		logfont.lfWeight = (long) (logfont.lfWeight*scaling_factor);
		logfont.lfWidth  = (long) (logfont.lfWidth *scaling_factor);
		CFont scaledFont;
		scaledFont.CreateFontA(logfont.lfHeight,logfont.lfWidth,logfont.lfOrientation,logfont.lfOrientation,logfont.lfWeight,logfont.lfItalic,logfont.lfUnderline,
			logfont.lfStrikeOut,logfont.lfCharSet,logfont.lfOutPrecision,logfont.lfClipPrecision,logfont.lfQuality,logfont.lfPitchAndFamily,logfont.lfFaceName);
		pDC->SelectObject(&scaledFont);
	}

	for(int i=1; i<=nr_lines; i++)
	{
		long x_separation = (long) (0.8*legendrect.left + 0.2*legendrect.right); // left 20% reserved for line
		long y_slottop    = (long) (legendrect.top + (i-1) * y_slot_height);
		long y_slotcenter = (long) (legendrect.top + (i-0.5) * y_slot_height);
		long y_slotbottom = (long) (legendrect.top + i * y_slot_height);

		// line properties from array
		COLORREF linecolor = p_PTD->lines(i).color;
		CPen pen_line (p_PTD->lines(i).penstyle,GetPenWidth(TPenWidth(p_PTD->lines(i).width)),linecolor);	
		pDC->SelectObject(&pen_line);
		pDC->MoveTo(legendrect.left,y_slotcenter);
		pDC->LineTo(x_separation,y_slotcenter);

		CPoint middle( (legendrect.left+x_separation)/2, y_slotcenter);
		int markerstyle = p_PTD->lines(i).pointstyle;
		int markersize = p_PTD->lines(i).pointsize;
		CRect markerrect = CRect(-markersize, +markersize, +markersize, -markersize);
		DrawDataPoint( pDC, markerstyle, middle, markerrect, linecolor);          // always mark first

		CRect textrect(x_separation, y_slottop, legendrect.right, y_slotbottom); // for this label

		this->DrawTextInRect(pDC, CString(p_PTD->lines(i).name),textrect,0);
	}
}

int CMyPlotToolView::DoNotRedraw() 
{ 
	return GetPlotToolDlg()->GetPopupChild()->DoNotRedraw(); 
}
void CMyPlotToolView::ForbidRedraw() 
{ 
	GetPlotToolDlg()->GetPopupChild()->ForbidRedraw(); 
}
void CMyPlotToolView::AllowRedraw() 
{ 
	GetPlotToolDlg()->GetPopupChild()->AllowRedraw(); 
}

void CMyPlotToolView::SetStatusBarText_Main(mystr text) 
{ 
	m_pPTD->GetPopupChild()->statusbartextbuffer = text;
	if (!m_pPTD->m_flag_statusinfo) 
		(m_pPTD->GetPopupChild()->m_StatusBar).SetText("",1,0);
	else
		(m_pPTD->GetPopupChild()->m_StatusBar).SetText(text,1,0); 
}

void CMyPlotToolView::SetStatusBarText_X (mystr text) 
{ 
	(m_pPTD->GetPopupChild()->m_StatusBar).SetText(text,2,0); 
}

void CMyPlotToolView::SetStatusBarText_Y (mystr text)
{ 
	(m_pPTD->GetPopupChild()->m_StatusBar).SetText(text,3,0); 
}


// **************************************************
// FUNCTIONS THAT PREPARE THE DRAWING DEVICE CONTENTS
// **************************************************
// client is Rect with Pixels as output BMP
void CMyPlotToolView::PrepareDC_Memory(CDC* pDC, CPrintInfo* pInfo)
{
	CRect clientRect(0,0,GetPlotToolDlg()->m_pixels_horizontal, GetPlotToolDlg()->m_pixels_vertical);
	PrepareDC_inRect(pDC, pInfo, clientRect);
}


// ************************************************
// COORDINATE SYSTEM CONVERSION   real <--> logical 
// ************************************************

// computes limits of datasets
void CMyPlotToolView::ComputeFullRange(int flag_equal)
{
	CPlotToolDlg* p_PTD = GetPlotToolDlg();
	p_PTD->RecomputeAllLineRanges();
	CRect plotrect;
	this->GetPlotRect(plotrect);

	// defaults for "no line to plot"
	double x_min = 0.;
	double x_max = 1.;
	double y_min = 0.;
	double y_max = 1.;

	int i = 1;
	while (i <= p_PTD->lines.Length())
	{
		// start with limits of dataset 1 (first nonaux line)
		if(! p_PTD->lines(i).IsAuxLine() )
		{ 
			x_min = p_PTD->lines(i).x_min;
			x_max = p_PTD->lines(i).x_max;
			y_min = p_PTD->lines(i).y_min;
			y_max = p_PTD->lines(i).y_max;
			i++;
			break;
		}
		else
		{
			i++;
		}
	}
	// other lines
	for(int j=i; j<=p_PTD->lines.Length(); j++)
	{
		if(! p_PTD->lines(j).IsAuxLine() )
		{
			if(x_min > p_PTD->lines(j).x_min) x_min = p_PTD->lines(j).x_min;
			if(x_max < p_PTD->lines(j).x_max) x_max = p_PTD->lines(j).x_max;
			if(y_min > p_PTD->lines(j).y_min) y_min = p_PTD->lines(j).y_min;
			if(y_max < p_PTD->lines(j).y_max) y_max = p_PTD->lines(j).y_max;
		}
	}
	fullrange = Box2D(Vector2D(x_min,y_min),Vector2D(x_max,y_max));

	// make x and y axis intervals same size 
	if(flag_equal)
	{
		double ylog_by_xlog = (double) Get_YRange_logical() / (double) Get_XRange_logical();
		double yran_by_xran = fullrange.SizeY() / fullrange.SizeX();

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

void CMyPlotToolView::ViewSetsRange(Box2D& range)
{
	shownrange = range;
	GetPlotToolDlg()->UpdateDrawRange(Get_XMin(), Get_XMax(), Get_YMin(), Get_YMax());
}

void CMyPlotToolView::SetShownRange(double xmin, double xmax, double ymin, double ymax)
{
	shownrange = Box2D( Vector2D(xmin,ymin), Vector2D(xmax,ymax) );
}

// sets the limits such that all lines are shown 
void CMyPlotToolView::ZoomToFullRange(int flag_equal, int flag_redraw)
{
	ComputeFullRange(flag_equal);
	ViewSetsRange(fullrange);
	if (flag_redraw)
	{
		SetCaller(mystr("REDRAW triggered by CMyPlotToolView::ZoomToFullRange\n"));
		RedrawWindow();
	}
}

// computes the major ticks for the axis, flag choses the the scaled distance is ((max-min)*scalingfactor)
void CMyPlotToolView::ComputeMajorTicksArray(TArray<double>& ticks, double min, double max, int flag, double scalingfactor)
{
	double tol = 1E-5;
	ticks.Flush();
	if(max<=min)
	{
	  ticks(1) = 0.; ticks(2) = 1.; return;
	}

	if (flag == 1) // linear scaling
	{
		double tick = GetStandardTickLength((max-min)/scalingfactor); // automatically compute tick interval

		double x = (ceil(min/tick))*tick;  // find first x to draw a tick

		double x_floor = x - tick; // tick before interval
		if( abs( (x_floor - min) / tick ) < tol ) // add tick if it is close enough to border
		{
			ticks.Add(x_floor);     
		}

		while (x <= max) // ticks within the interval
		{
			if ( abs((x/tick)) < 1e-5 )  x=0.; // make sure that tick is at exactly 0.0
			ticks.Add(x);
			x += tick;
		}

		double x_ceil = x + tick; // tick after interval
		if( abs( (x_ceil - max) / tick ) < tol ) // add tick if it is close enough to border
		{
			ticks.Add(x_ceil);     
		}
	}
}

double CMyPlotToolView::GetStandardTickLength(double interval)
{
	if (interval == 0.) return 1;

	int exponent = (int) floor (log (fabs(interval)) / log(10.));
	double scalingfactor = pow(10., -exponent);
	double scaledval = interval * scalingfactor;  // gives 1<= scaledval <10

	double newval = 0.2;                // interval 1<= scaledval <= 2
	if (scaledval > 2.) newval = 0.5;   // interval 2<  scaledval <= 5
	if (scaledval > 5.) newval = 1.;    // interval 

	return newval / scalingfactor;
}

// NEW FUNCTIONS HERE - WHEN THEY ARE WORKING PUT THEM IN OTHR PLACE

// ***************************************************************************************
// context menu on the View side - 
// ***************************************************************************************
#define popups_as_resource_not
void CPlotToolPopupChild::OnContextMenu(CWnd* pWnd, CPoint point)
{
	CMenu contextmenu;
#ifdef popups_as_resource
	contextmenu.LoadMenu(IDR_CONTEXT_VIEW);
#else
	contextmenu.CreatePopupMenu();
	UINT nFlags;
	
	contextmenu.AppendMenu(MF_STRING, ID_CM_ZOOMTOFULLRANGE,	"Zoom to Full Range");
	
  nFlags = (MF_STRING | ( GetParent()->m_flag_show_legend ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_SHOWLEGEND,				"Show Legend"); //(V) or (-)
 //	nFlags = (MF_STRING | ( GetParent()->m_flag_fastdraw ? MF_CHECKED : MF_UNCHECKED ));
	//contextmenu.AppendMenu(nFlags,		ID_CM_FASTDRAW,					"Draw Fast");
  nFlags = (MF_STRING | ( GetParent()->m_flag_statusinfo ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_INFOPLOTTOOL,				"Show Status Bar Information"); //(V) or (-)

	nFlags = (MF_STRING | ( GetParent()->m_flag_use_autoupdate ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_AUTOUPDATE,				"Automatic Update");
	nFlags = (MF_STRING | ( GetParent()->m_flag_autorescale ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_AUTORESCALE,			"Automatic Rescale");
 	nFlags = (MF_STRING | ( GetParent()->m_flag_draw_every_nth ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_DRAWSPARSE,			"Sparse Data Drawn");
 	nFlags = (MF_STRING | ( GetParent()->m_flag_mark_every_nth ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_MARKSPARSE,			"Sparse Data Markers");
 	nFlags = (MF_STRING | ( GetParent()->m_flag_vertical_marker ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags,		ID_CM_MARKCURRENT,			"Mark Current Time");
	nFlags = (MF_STRING | ( GetParent()->m_flag_upto_current ? MF_CHECKED : MF_UNCHECKED ));
	contextmenu.AppendMenu(nFlags, ID_CM_DRAWUPTOCURRENT,	"Draw only up to current Time");
	
	if(! GetParent()->IsWindowVisible())
	{
		contextmenu.AppendMenu(MF_STRING, ID_CM_SHOWDIALOG,				"Show Dialog");
	}
	else
	{
		contextmenu.AppendMenu(MF_STRING, ID_CM_HIDEDIALOG,				"Hide Dialog");
	}
	contextmenu.AppendMenu(MF_STRING, ID_CM_AXISEQUAL,				"Axis Equal");
	contextmenu.AppendMenu(MF_STRING, ID_CM_EXPORTTOFILE,			"Export Graph to File");
	contextmenu.AppendMenu(MF_STRING, ID_CM_PRINTGRAPH,				"Print Graph");
	
	//CBitmap bmp1; bmp1.LoadBitmap(IDB_SAVE_32);
	//contextmenu.AppendMenu(MF_BITMAP, ID_CM_EXPORTTOFILE, &bmp1);
	//CBitmap bmp2; bmp2.LoadBitmap(IDB_PRINT_32);
	//contextmenu.AppendMenu(MF_BITMAP, ID_CM_PRINTGRAPH, &bmp2);

	//contextmenu.AppendMenu(MF_SEPARATOR,0,(char*)0);
#endif

	contextmenu.TrackPopupMenu(TPM_LEFTALIGN,point.x,point.y,this,NULL);
}

void CPlotToolPopupChild::OnCMZoomToFullRange() 
{	
	if (GetView() != NULL)
	{
		GetView()->ZoomToFullRange(0,1);
	}
}

void CPlotToolPopupChild::OnCMShowLegend() 
{ 
	if(GetParent()->m_flag_show_legend)
	{
		GetParent()->m_flag_show_legend = 0; // change from show to dont show
		GetView()->RedrawWindow();
	}
	else
	{
		GetParent()->m_flag_show_legend = 1; // change from dont show to show
		GetView()->RedrawWindow();
	}
	GetParent()->UpdateData(FALSE);        // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMShowStatusInfo() 
{ 
	if(GetParent()->m_flag_statusinfo)
	{
		GetParent()->m_flag_statusinfo = 0; // change from show to dont show
		m_StatusBar.SetText("",1,0);
	}
	else
	{
		GetParent()->m_flag_statusinfo = 1; // change from dont show to show
		m_StatusBar.SetText(statusbartextbuffer,1,0);
	}
	GetParent()->UpdateData(FALSE);        // keep the obsolete Checkbox in correct state ...
	m_StatusBar.UpdateWindow();
}


void CPlotToolPopupChild::OnCMFastDraw() 
{ 
	if(GetParent()->m_flag_fastdraw)
	{
		GetParent()->m_flag_fastdraw = 0; // change from show to dont show
	}
	else 
	{
		GetParent()->m_flag_fastdraw = 1; // change from dont show to show
	}
	GetParent()->UpdateData(FALSE);     // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMAutoUpdate() 
{ 
	if(GetParent()->m_flag_use_autoupdate)
		GetParent()->m_flag_use_autoupdate = 0; // change from show to dont show
	else
		GetParent()->m_flag_use_autoupdate = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);           // keep the obsolete Checkbox in correct state ...
	GetParent()->StartStopAutoUpdate();       // start or stop the update timer
}

void CPlotToolPopupChild::OnCMAutoReScale() 
{ 
	if(GetParent()->m_flag_autorescale)
		GetParent()->m_flag_autorescale = 0; // change from show to dont show
	else
		GetParent()->m_flag_autorescale = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);        // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMDrawSparse() 
{ 
	if(GetParent()->m_flag_draw_every_nth)
		GetParent()->m_flag_draw_every_nth = 0; // change from show to dont show
	else
		GetParent()->m_flag_draw_every_nth = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);           // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMMarkSparse() 
{ 
	if(GetParent()->m_flag_mark_every_nth)
		GetParent()->m_flag_mark_every_nth = 0; // change from show to dont show
	else
		GetParent()->m_flag_mark_every_nth = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);           // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMMarkCurrent() 
{ 
	if(GetParent()->m_flag_vertical_marker)
		GetParent()->m_flag_vertical_marker = 0; // change from show to dont show
	else
		GetParent()->m_flag_vertical_marker = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);            // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMDrawUpToCurrent() 
{ 
	if(GetParent()->m_flag_upto_current)
		GetParent()->m_flag_upto_current = 0; // change from show to dont show
	else
		GetParent()->m_flag_upto_current = 1; // change from dont show to show
	GetParent()->UpdateData(FALSE);         // keep the obsolete Checkbox in correct state ...
}

void CPlotToolPopupChild::OnCMShowDialog() 
{ 
	GetParent()->ShowWindow(SW_SHOW);
}
void CPlotToolPopupChild::OnCMHideDialog() 
{ 
	GetParent()->ShowWindow(SW_HIDE);
}

void CPlotToolPopupChild::OnCMAxisEqual() 
{ 
	if (GetView() != NULL)
	{
		GetView()->ZoomToFullRange(1,1);
	}
}

void CPlotToolPopupChild::OnCMExportToFile() 
{ 
	GetParent()->ExportToFile();
}

void CPlotToolPopupChild::OnCMPrintGraph() 
{ 
	/*GetParent()->*/Print(); 
}

void CPlotToolPopupChild::SetCaller(mystr calleri)
{
	GetView()->SetCaller(calleri); 
}


void CPlotToolDlg::OnContextMenu(CWnd* pWnd, CPoint point)
{
	CMenu contextmenu;
#ifdef popups_as_resource
	contextmenu.LoadMenu(IDR_CONTEXT_PLOTTOOLDIALOG);
#else
	contextmenu.CreatePopupMenu();
	contextmenu.AppendMenu(MF_STRING, ID_CM_UPDATEDATA, "Update Data");
	contextmenu.AppendMenu(MF_STRING, ID_CM_SHOWGRAPH, "Show Graph");
	contextmenu.AppendMenu(MF_STRING, ID_CM_SAVESETTINGS, "Save Settings");
	contextmenu.AppendMenu(MF_SEPARATOR,0,(char*)0);
#endif

	contextmenu.TrackPopupMenu(TPM_LEFTALIGN| TPM_RIGHTBUTTON,point.x,point.y,this,NULL);
}

void CPlotToolDlg::OnCmUpdateData() 
{ 
	for(int i=1; i<=lines.Length(); i++)
	{
		if (lines(i).IsSensor()) lines(i).ReadLine();
	}
	Redraw();
}

void CPlotToolDlg::OnCmShowGraph()
{
	ShowGraph();
}

void CPlotToolDlg::OnCmSaveSettings()
{
	UpdateData(TRUE);
	WriteData();
}

void CPlotToolDlg::OnCmMarkcurrent()
{
	// TODO: Fügen Sie hier Ihren Befehlsbehandlungscode ein.
}

void CPlotToolDlg::OnBnClickedButtonExportOptions()
{
	CGraphExportOptionsDlg options;
	options.SetParent(this);
	if(options.DoModal() == IDCANCEL)
		return;
}

void CPlotToolDlg::SetCaller(mystr calleri)
{
	m_pChild->GetView()->SetCaller(calleri); 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


