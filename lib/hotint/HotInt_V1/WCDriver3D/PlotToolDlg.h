//#**************************************************************
//# filename:             PlotToolDlg.h
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
 

#pragma once
#include "resource.h"
#include "afxcmn.h"
#include "afxwin.h"
#include "mathfunc.h"
#include "MyBaseView.h" // base class for View
#include "PlotToolDlg_aux.h" // MyColorList, default 

// !!! under construction !!!
#define preview_not                          // OBSOLETE, small preview window in PlotTool, Idea was that a single PlotTool Dialog nests several Graphs
                                             // as a consequence almost all variables for draw are also in the main dialog, not the graph dialog
#define sensors_from_solution_file_NOT       // TO BE REMOVED (when sensors can load from file...)

// TIMER
#define PLOT_REFRESH_TIMER_ID (WM_USER + 111)
#define WM_UPDATE_PLOT (WM_USER + 11)

// Dialog Elements that are vreated manually
#define ID_VIEW_ON_POPUP 1234
#define ID_MENU_PARENT_ON_POPUP 1235
#define IDC_PLOTTOOLPOPUP_STATUSBAR 1236

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ owner drawn ComboBoxes:       ColorPicker, LineThicknessPicker, LineStylePicker
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CBColorPicker : public CComboBox // item data == COLORREF
{
public: // lifecycle 
	CBColorPicker():CComboBox()	{ }
	//
public: // lifecycle II
	int InitializePalette(MyColorList& palettei, int cutoff = -1)
	{
	  palette.CopyFrom(palettei);
		return UpdateItems(cutoff);
	}
	int UpdateItems(int cutoff = -1)
	{
		ResetContent();//Clear();
		for(int i=1; i<=palette.N(); i++)
		{
			if ( cutoff > 0 && i > cutoff ) return cutoff;
			int nr = AddString(palette.Name(i));
			SetItemData(nr, palette.GetColRef(i) );
		}
		int count = GetCount();
		return count;
	}
	// OVERRIDES FOR OWNERDRAW
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct);               // Owner Draw: rectangle with picked
	virtual void MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct);
	virtual int CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct);
	virtual void DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct);
private:
	MyColorList palette;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Penstyles already defined  (wingdi.h)
class CBLineStylePicker : public CComboBox // item data == PS_STYLE
{
public: // lifecycle
	CBLineStylePicker():CComboBox() { activecolor = 0; }
	//
public: // lifecycle II
	void Initialize() // Add default line styles
	{
		Clear();
		SetItemData(AddString("solid"), PS_SOLID );
		SetItemData(AddString("dash"), PS_DASH );
		SetItemData(AddString("dot"), PS_DOT );
		SetItemData(AddString("dashdot"), PS_DASHDOT );
		SetItemData(AddString("dashdotdot"), PS_DASHDOTDOT );
		SetItemData(AddString("none"), PS_NULL );
	}
	// OVERRIDES FOR OWNERDRAW
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct);               // Owner Draw: rectangle with picked
	virtual void MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct);
	virtual int CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct);
	virtual void DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct);
public: // pass selected color into CB
	COLORREF activecolor;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef enum 
{
	PTS_NONE = 1,           // -----
	PTS_X = 2,              // --x--
	PTS_CROSS = 3,          // --+--
	PTS_DOT = 4,            // --@--   filled
	PTS_CIRCLE = 5          // --o--   empty
} TMarkerStyle;

class CBPointStylePicker : public CComboBox // itemdata == TPointStyle
{
public: // lifecycle
	CBPointStylePicker():CComboBox() { activecolor = 0; activelinestyle = 0; }
	//
public: // lifecycle II
	void Initialize() // Add default line styles
	{
		ResetContent();//Clear();
		SetItemData(AddString("none"), PTS_NONE );
		SetItemData(AddString("X"), PTS_X );
		SetItemData(AddString("+"), PTS_CROSS );
		SetItemData(AddString("@"), PTS_DOT );
		SetItemData(AddString("o"), PTS_CIRCLE );
	}
	// OVERRIDES FOR OWNERDRAW
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct);               // Owner Draw: rectangle with picked
	virtual void MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct);
	virtual int CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct);
	virtual void DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct);
public: // pass selected color into CB
	COLORREF activecolor;
	int activelinestyle;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef enum 
{
	PT_VERYFINE = 1,
	PT_FINE = 2,
	PT_NORMAL = 3,
	PT_THICK = 4,
	PT_VERYTHICK = 5
} TPenWidth;

class CBLineThicknessPicker : public CComboBox // item data == TPenWidth
{
public: // lifecycle
	CBLineThicknessPicker():CComboBox() { activecolor = 0; }
	//
	//
public: // lifecycle II
	//
	void Initialize() // Add default line styles
	{
		ResetContent();//Clear();
		SetItemData(AddString("very fine"), PT_VERYFINE );
		SetItemData(AddString("fine"), PT_FINE );
		SetItemData(AddString("normal"), PT_NORMAL );
		SetItemData(AddString("thick"), PT_THICK );
		SetItemData(AddString("very thick"), PT_VERYTHICK );
	}
	// OVERRIDES FOR OWNERDRAW
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct);               // Owner Draw: rectangle with picked
	virtual void MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct);
	virtual int CompareItem(LPCOMPAREITEMSTRUCT lpCompareItemStruct);
	virtual void DeleteItem(LPDELETEITEMSTRUCT lpDeleteItemStruct);
public: // draw in picked color
	COLORREF activecolor;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+ class Line_Data_Prop:   contains properties and data for a single line of the plot,
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
enum PlotTool_LineType 
{ 
	TLT_ty = 0,                   // plot 1 data column over time 
	TLT_xy = 1,                   // plot 2 data columns versus

	TLT_sens = 2,                 // data source is sensors own buffer
	TLT_solution = 4,             // data source is solution file
	TLT_external = 8,             // data source is external file
	TLT_aux = 16                  // use a parsed math function
};
// datastruct for line
class Line_Data_Prop 
{
public: // lifecycle
	Line_Data_Prop() : x_data(0), y_data(0)
	{
		Reset();
	}
	Line_Data_Prop(const Line_Data_Prop& other) : x_data(0), y_data(0)
	{
		Reset();
		CopyFrom(other);
	}
	Line_Data_Prop& operator = (const Line_Data_Prop& other)
	{
		if (this==&other) {return *this;}
		CopyFrom(other);
		return *this;
	}
	~Line_Data_Prop()
	{
	}

public: // lifecycle II
	Line_Data_Prop(class CPlotToolDlg* dlg) : x_data(0), y_data(0)
	{
		Reset();
		dialog = dlg;
	}
	void ResetData()
	{
		x_data = 0;
		y_data = 0;
		t_data = 0;
		x_min = 0.; x_max = 1.;
		y_min = 0.; y_max = 1.;
		lines_processed = 0;
	}
	void ResetProp()
	{
		name = "";
		color = RGB(0,0,0);
		width = PT_FINE;
		penstyle = PS_SOLID;
		pointstyle = PTS_NONE;
		pointsize = 40;
		filename = "";
		x_colnr = 0; y_colnr = 0;
		linetype = PlotTool_LineType(0);
	}
	void Reset()
	{
		dialog = NULL;
		ResetProp();
		ResetData();
	}
	void CopyFrom(const Line_Data_Prop& other)
	{
		dialog = other.dialog;
		name = other.name;

		color = other.color;
		width = other.width;
		penstyle = other.penstyle;
		pointstyle = other.pointstyle;
		pointsize = other.pointsize;

		x_data = other.x_data;
		y_data = other.y_data;
		t_data = other.t_data;
		x_min = other.x_min;
		x_max = other.x_max;
		y_min = other.y_min;
		y_max = other.y_max;
		lines_processed = other.lines_processed;

		filename = other.filename;
		x_colnr = other.x_colnr;
		y_colnr = other.y_colnr;
		linetype = other.linetype;
	}

	int ReadFromEDC(ElementDataContainer* edc)
	{
		// disable popup warnings
		int edc_warn_lvl = edc->GetEDCWarningLevel();
		edc->SetEDCWarningLevel(0);

		ResetProp();
		int found;
		
		found = edc->Find("name"); if(found && edc->Get(found).IsText()) name = edc->TreeGetString("name");
		
		found = edc->Find("color");	if(found)	color = edc->TreeGetInt("color");
		found = edc->Find("width");	if(found)	width = edc->TreeGetInt("width");
		found = edc->Find("penstyle"); if(found) penstyle = edc->TreeGetInt("penstyle");
		found = edc->Find("pointstyle"); if(found) pointstyle = edc->TreeGetInt("pointstyle");
		found = edc->Find("pointsize"); if(found)	pointsize = edc->TreeGetInt("pointsize");
		
		found = edc->Find("source_file_name"); if(found && edc->Get(found).IsText()) filename = edc->TreeGetString("source_file_name");
		found = edc->Find("x_column"); if(found) x_colnr = edc->TreeGetInt("x_column"); 
		found = edc->Find("y_column"); if(found) y_colnr = edc->TreeGetInt("y_column");
		found = edc->Find("linetype"); if(found) linetype = (PlotTool_LineType) edc->TreeGetInt("linetype");

		edc->SetEDCWarningLevel(edc_warn_lvl);
		return 1;
	}
	int WriteToEDC(ElementDataContainer* edc)
	{
		ElementData ed;

		ed.SetText(name,"name"); edc->Add(ed);
		
		ed.SetInt(color, "color"); ed.SetToolTipText("in COLORREF"); edc->Add(ed);
		ed.SetInt(width,"width"); ed.SetToolTipText("{1=VERY FINE .. 5=VERY THICK}"); edc->Add(ed);
		ed.SetInt(penstyle,"penstyle"); ed.SetToolTipText("{0=solid .. 5=noline}"); edc->Add(ed);
		ed.SetInt(pointstyle,"pointstyle"); ed.SetToolTipText("{1=NONE .. 5=CIRCLE}"); edc->Add(ed);
		ed.SetInt(pointsize,"pointsize"); ed.SetToolTipText("in logical units"); edc->Add(ed);

		ed.SetText(filename,"source_file_name"); edc->Add(ed);
		ed.SetInt(x_colnr,"x_column"); edc->Add(ed);
		ed.SetInt(y_colnr,"y_column"); edc->Add(ed);
		ed.SetInt(linetype,"linetype"); ed.SetToolTipText("binary flags"); edc->Add(ed);
	  return 1;
	}

public: // access
	PlotTool_LineType GetType() const { return linetype; }
	void SetType(int i) { linetype = PlotTool_LineType(i); }
	void SetType(PlotTool_LineType i) { linetype = i; }
	int IsTY() { return (!IsXY()); }
	int IsXY() { return (linetype & TLT_xy); }
	int IsSensor() { return (linetype & TLT_sens); }
	int IsSolFile() { return (linetype & TLT_solution); }
	int IsExtFile() { return (linetype & TLT_external); }
	int IsAuxLine() { return (linetype & TLT_aux); }

public: // link data
	void SetXDataPtr(TArray<double>* source_x) { x_data = source_x; }
	void SetYDataPtr(TArray<double>* source_y) { y_data = source_y; }
	void SetTDataPtr(TArray<double>* source_t) { t_data = source_t; }
	
public: // functionality
	int UpdateMinMax(int append_at_line=1)	// calculate value range
	{
		if(x_data->Length() == 0) // empty line
		{
			x_min = 0.0; x_max = 1.0;
			y_min = 0.0; y_max = 1.0;
			return 1;
		}
		else if(x_data->Length() == 1) // single point
		{
			x_min = x_data->Get(1)-0.5; x_max = x_data->Get(1)+0.5;
			y_min = y_data->Get(1)-0.5; y_max = y_data->Get(1)+0.5;
			return 1;
		}
		else if ( (x_data->Length() >= 2) && (lines_processed < 2)) // for longer datasets: reset min/max if lines are not processed ( can force update ) 
		{
			x_min = x_data->Get(1); x_max = x_data->Get(1);
			y_min = y_data->Get(1); y_max = y_data->Get(1);
		}

		int datalength = x_data->Length();
		int counter = 0;
		for(int i=lines_processed+1; i<=datalength; i++)
		{
			if(x_data->Get(i) < x_min) x_min = x_data->Get(i);
			if(x_data->Get(i) > x_max) x_max = x_data->Get(i);
			if(y_data->Get(i) < y_min) y_min = y_data->Get(i);
			if(y_data->Get(i) > y_max) y_max = y_data->Get(i);
			counter++;
		}
		lines_processed = datalength;
		return counter; // newly read lines
	}
	// reads x and y columns of data line from a file, may append
	int ReadLine(int flag_append = 1)
  {
    if(IsSensor()) return UpdateMinMax(lines_processed);
		else return 0;
	}
	
  // release the data for external file data ( which is COPIED here )
	void ReleaseData()
	{
		if(IsExtFile())
		{
			if(x_data) delete x_data;
			if(y_data) delete y_data;
			//if(t_data) delete t_data;
		}
	}

public: // variables
	class CPlotToolDlg* dialog;								// pointer to dlg ( use some member fctns )

	// dataset
	TArray<double>* x_data;	  					// array for x-values - holds the x-values for the datapoint (is time data in t/y-plot)
	TArray<double>* y_data; 			  		// array for y-values - holds the y-values for the datapoint
	TArray<double>* t_data;             // array for t-values - required for the option to mark the current time or draw only up to the current time. NOT used for drawing

	// displayed/graphical properties
	mystr name;													// displayed name of the line
	COLORREF color;											// color of the line
	int width;                          // width of the line ( picker-index, not logical points - this is done by CView )
	int penstyle;												// style of the line
	int pointstyle;                     // style for data point marker
	int pointsize;                      // size of the datapoint marker

	// processeed
	double x_min, x_max, y_min, y_max;	// min and max values
	int lines_processed;                // min & max checked up to this line

	// data properties
	PlotTool_LineType linetype;        // type of the line - binary flags, see enum

	mystr filename;											// name of the file where data is stored
	int x_colnr, y_colnr;								// column numbers within that file
}; 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CPlotToolPopupChild:       this is the dialog that nests the Graph (CMyView object), 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PlotToolPopupChild-Dialogfeld
class CPlotToolPopupChild : public MyBaseDialogForView
{
	enum { IDD = IDD_DIALOG_2D_VIEW };
	DECLARE_DYNAMIC(CPlotToolPopupChild)

public:
	CPlotToolPopupChild(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CPlotToolPopupChild();
public:
	virtual BOOL OnInitDialog();
	virtual BOOL DestroyWindow();

protected:
	//virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	DECLARE_MESSAGE_MAP()

// ACCESS // LINKED DIALOGS & OBJECTS
protected:
	class CPlotToolDlg* m_pParent;
public:
	void SetView(class CMyPlotToolView* myView) { m_pMyView = (CMyBase2DView*) myView; }
	void SetParent(CPlotToolDlg* parent) { m_pParent = parent; }
	class CMyPlotToolView* GetView() { return (class CMyPlotToolView*) m_pMyView; } // we know that the derived dialog also uses a derived view class
	class CPlotToolDlg* GetParent() { return m_pParent; }

// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
public:
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);  
	afx_msg void OnMove(int x, int y);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnClose();
	afx_msg void OnPaint();
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
protected:
// catch other framework messages
	virtual BOOL OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult); 

// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
	void PlaceElements(); 

// CONTEXT MENU
public:
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	afx_msg void OnCMZoomToFullRange();
	afx_msg void OnCMShowLegend();
	afx_msg void OnCMShowStatusInfo();
	afx_msg void OnCMFastDraw();
	afx_msg void OnCMAutoUpdate();
	afx_msg void OnCMAutoReScale();
	afx_msg void OnCMDrawSparse();
	afx_msg void OnCMMarkSparse();
	afx_msg void OnCMMarkCurrent();
	afx_msg void OnCMDrawUpToCurrent();
	afx_msg void OnCMShowDialog();
	afx_msg void OnCMHideDialog();
	afx_msg void OnCMAxisEqual();
	afx_msg void OnCMExportToFile();
	afx_msg void OnCMPrintGraph();
//END CONTEXTMENU


// Diagnose functions - Functions that trigger redraw
	virtual void SetCaller(mystr calleri);
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CGraphExportOptionsDlg:       this is the dialog that nests the Export Options, 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CGraphExportOptionsDlg-Dialogfeld

class CGraphExportOptionsDlg : public CDialog
{
	enum { IDD = IDD_DIALOG_EXPORTOPTIONS };
	DECLARE_DYNAMIC(CGraphExportOptionsDlg)

public:
	CGraphExportOptionsDlg(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CGraphExportOptionsDlg();
public:
	virtual BOOL OnInitDialog();
	virtual void Reset();
	void LoadData();																			 //Get data from parent dlg
	void WriteData();																			 //Put data to parent dlg

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung

	DECLARE_MESSAGE_MAP()

	// ACCESS // LINKED DIALOGS & OBJECTS
protected:
	class CPlotToolDlg* m_pParent;
public:
	void SetParent(CPlotToolDlg* parent) { m_pParent = parent; }
	class CPlotToolDlg* GetParent() { return m_pParent; }


	// *****************
	// *** Group: Export
	// *****************
	CString m_export_path;                                  // DDX: String containing the destination for plottool images
	CString m_export_filename;                              // DDX: String containing the file name for plottool images

	int m_pixels_horizontal;
	int m_pixels_vertical;

	int m_flag_export_jpg;
	int m_flag_export_png;
	int m_flag_export_bmp;
	int m_flag_export_emf;

	afx_msg void OnEnKillfocusEditResx();
	afx_msg void OnEnKillfocusEditResy();
	afx_msg void OnBnClickedButtonFiguresizeasscreen();
	afx_msg void OnBnClickedButtonPathbrowse();
	//int CALLBACK BrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData);

	afx_msg void OnBnClickedOk();
	CButton m_snapnow;
	afx_msg void OnBnClickedButtonSnapshotnow();
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CPlotToolDlg:       this is the dialog that nests the controls for the graph, 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CPlotToolDlg : public CDialog 
{
	enum { IDD = IDD_DIALOG_PLOTTOOL };
	DECLARE_DYNAMIC(CPlotToolDlg)

public:
	CPlotToolDlg(CWnd* pParent = NULL);   // Standardkonstruktor
	virtual ~CPlotToolDlg();
public:
	virtual BOOL OnInitDialog();
	virtual BOOL DestroyWindow();

	//COMMUNICATION WITH SETTINGs (EDC)
	virtual void Reset();
	void LoadData();                     //Get data from WCDinterface
	void WriteData();                    //Put data to WCDinterface
  void UpdateDialogData();             //Update from MBS_EDC in case the edc is changed while PlotToolDialog is open

	// Dialogfelddaten
protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV-Unterstützung
	DECLARE_MESSAGE_MAP()

	// ACCESS // LINKED DIALOGS & OBJECTS
private:  
	WCDInterface* pWCDI;
	CGLDrawWnd* pGLDrawWnd;
	CWCDriver3DDlg* pWCDDlg;
	class CPlotToolPopupChild* m_pChild;
public:
	void SetWCDI(WCDInterface * pWCDI_) { pWCDI = pWCDI_; }
	void SetGLDrawWnd(CGLDrawWnd * pGLDrawWnd_) { pGLDrawWnd = pGLDrawWnd_; }
	void SetWCDDlg(CWCDriver3DDlg* pWCDDlg_) { pWCDDlg = pWCDDlg_;}
	void SetPopupChild(CPlotToolPopupChild* pChild) { m_pChild = pChild; }
	WCDInterface* GetWCDI() { return pWCDI; }
	CGLDrawWnd* GetGLDrawWnd() { return pGLDrawWnd; }
	CWCDriver3DDlg* GetWCDDlg() { return pWCDDlg; }
	CPlotToolPopupChild* GetPopupChild() { return m_pChild; }
#ifdef preview
	class CMyView* m_pMyView; // must be derived from CView
	void SetPreview(CMyView* preview) { m_pMyView = preview; }
	CMyView* GetPreview() { return m_pMyView; }
#endif

// TIMER FUNCTIONS (redraw, update)
public:
	void StartRefreshTimer();
	void KillRefreshTimer();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg LRESULT OnUpdateMsg(WPARAM mode, LPARAM str);

// WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnMove(int x, int y);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);                       
	afx_msg void OnClose();
	void ExternClose();
	afx_msg void OnPaint();
	int m_init_done;
	int m_do_not_redraw;
	int m_hasPopupChild;
	int _isclosing;
	int _externclose;

// placement of the elements ( CMyView object & Dialog ) so the two dialogs stick together
	void PlaceElements(); 

// 
	void Redraw(); // redraw the entire view window
	void ShowGraph();
	void ScaleToFull(); // scales the Graph to show all data
  void UpdateLineData(); // updates the data content - manually e.g. during computaion
	void ExportToFile();


// FUNCTIONS FOR PRINTING
	void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo) {return;}
	void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo) {return;}
	void Print();

// AUXILARY FUNCTIONS
	mystr SetSolFileNameFromEDC()
	{
		filename_solution = pWCDI->MBS_EDC_TreeGetString("GeneralOptions.Paths.sensor_output_path");
		filename_solution += pWCDI->MBS_EDC_TreeGetString("SolverOptions.Solution.SolutionFile.output_filename");
		return (LPCTSTR) filename_solution;
	}

public: // Variables for the Dialog Items (DDX) , Functions for the Buttons, associated variables and functions 
	// ***********************
	// *** Menu Bar
	// ***********************
	afx_msg void OnLayoutsSave();
	void SaveLayoutToFile(mystr& filename);
	afx_msg void OnLayoutsLoad();
	void LoadLayoutFromFile(mystr& filename);
	int _isloading;
	int IsLoading() { return _isloading; }

	// ***********************
	// *** Group: Data Sources
	// ***********************
	int m_radio_source;                                     // DDX: radio button: source = from solution file OR from external file
	CString filename_extern;                                // DDX: name of the file for external
	CString filename_solution;															// name of solution file 
	CString m_label_source;																	// DDX: label over listbox displays source
	CListBox m_list_columns;																// DDX: ListBox displaying list of Sensors OR list of columns in selected file

	afx_msg void OnBnClickedRadioSource();									// radio buttons to pick source - use "old" external file
	afx_msg void OnBnClickedButtonFilebrowse();             // Button: "..." openfile dialog, update filename, fill listbox - file can change
	afx_msg void OnBnClickedButtonLineAddXY();							// adds a X/Y line to drawn lines

	void AddLine(int x_colnr, int y_colnr, int flag_isxy);  // 
	void SetLineType(class Line_Data_Prop& the_line, int flag_isxy);
	void SetLineName(class Line_Data_Prop& the_line, int x_colnr, int y_colnr);
	void LinkLineData(class Line_Data_Prop& the_line);

	class CMatrixFile* extfile;                             // holds the matrix the external file
	int LoadExternalToMatrixFile();                         // reads the external file

	class CMatrixFile* solfile;                             // holds the matrix the solution file 
	int LoadSolutionToMatrixFile();                         // reads the solution file
	int UpdateSolutionToMatrixFile();                       // for incremental read
	int solfile_linesread;                                  // for incremental read

	// ****************************
	// *** ListBox: Available Lines
	// ****************************
	afx_msg void OnBnClickedButtonLineAdd();                // adds selected line to drawn lines
	afx_msg void OnLbnDblclkListColumns();									// adds selected line to drawn lines
	afx_msg void OnLbnSelchangeListColumns();								// keep track of previously selected lines ( for x-y plot )
	int act_selected_column;                                // selection marker (ListBox)
	int prev_selected_column;                               // selection marker (ListBox)
	void WriteSensorNamesToColList();												// writes "Time" and all sensors to listbox
	void WriteFileColumnNamesToColList();			              // writes all columns from solution file or external file to listbox
	void WriteExternalColumnNamesToColList();								// writes all columns from external file to listbox
	
	// *************************
	// *** Group: Title and Axis
	// *************************
	CString m_title;                                        // DDX: string for graph title
	LOGFONT m_titlefont;                                    // font for title string
	LOGFONT m_axisfont;                                     // (dependent) font for axis text
	LOGFONT m_ticksfont;                                    // (dependent) font for axis labels (numbers)

	double m_titlefont_scalingfactor;                       // factor for relative fontsize
	double m_ticksfont_scalingfactor;                       // factor for relative fontsize

	int m_border_thickness;                                 // line thickness for the border
	int m_axis_at_origin;                                   // flag if axis move to within the rectangle if the origin is within the rectancle
  double m_axis_overdraw;                                 // axis are longer then rect

	int m_axis_label_major;                                 
  int m_axis_label_minor;
	double m_axis_ticks_overdraw;
	int m_digits_x_label;
	int m_digits_y_label;

	int m_axis_minor_ticks_x;
	int m_axis_minor_ticks_y;

	CString m_x_axis_title;                                 // DDX: string for the x-axis title
	double x_min,x_max;                                     // DDX: range for the x-axis 
	int flag_use_x_majorticks;
	double x_major_interval;
	int flag_use_x_minorticks;
	double x_minor_interval;
	
	CString m_y_axis_title;                                 // DDX: string for the y-axis title
	double y_min,y_max;                                     // DDX: range for the y-axis 
	int flag_use_y_majorticks;
	double y_major_interval;
	int flag_use_y_minorticks;
	double y_minor_interval;

	afx_msg void OnBnClickedButtonFont();                   // font selection dialog
	void SelectFontDialog(LOGFONT& logfont);                // creates the Font selection dialog
	void AssignOtherFonts();                                // Assign all fonts from chosen axisfont (scaled)

	afx_msg void OnEnKillfocusEditTitle() { UpdateData(TRUE); SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnEnKillfocusEditTitle\n")); Redraw(); }
	afx_msg void OnEnKillfocusEditXAxis() { UpdateData(TRUE); SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnEnKillfocusEditXAxis\n")); Redraw(); } 
	afx_msg void OnEnKillfocusEditYAxis() { UpdateData(TRUE); SetCaller(mystr("REDRAW triggered by CPlotToolDlg::OnEnKillfocusEditYAxis\n")); Redraw(); }
	afx_msg void OnEnKillfocusEditXaxisLower() { EditBoxSetsRange(); }
	afx_msg void OnEnKillfocusEditXaxisUpper() { EditBoxSetsRange(); }
	afx_msg void OnEnKillfocusEditYaxisLower() { EditBoxSetsRange(); }
	afx_msg void OnEnKillfocusEditYaxisUpper() { EditBoxSetsRange(); }
	virtual void EditBoxSetsRange();                       // apply changes made in the textboxes to view 
	void UpdateDrawRange(double xmin, double xmax, double ymin, double ymax); // write actual limits of graph to textbox 
	void UpdateMajorIntervals(double xmajor, double ymajor); // write artual ticksize to textbox

	// ************************
	// *** ListBox: Drawn Lines
	// ************************
	CListCtrl m_list_lines;                                 // DDX: ListBox displaying list of lines to be drawn
	afx_msg void OnLvnItemchangedListLines(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnLvnKeydownListLines(NMHDR *pNMHDR, LRESULT *pResult); // press DEL to delete a line
	int m_selected_line;

	void WriteToDlgPropertiesofLine(int selected_line);     // apply properties of the selected line to control elements
	afx_msg void OnEnKillfocusEditLineName();
	CString m_label_lineprop;
	CString m_edit_name_of_selection;                       // DDX: EditBox name of the selected line - also name of line in Legend

	afx_msg void OnCbnSelendokComboLineColor();
	afx_msg void OnCbnSelendokComboLineThickness();
	afx_msg void OnCbnSelendokComboLineStyle();
	afx_msg void OnCbnSelendokComboPointStyle();
	afx_msg void OnEnKillfocusEditPointsize();
  CBColorPicker line_color;                               // color picker dropdown - derived from CComboBox
	MyColorList palette;
	int line_color_index;
	MyColorList mycolors;
	COLORREF GetDefaultColor(int i) { return (COLORREF) mycolors.GetColRef(i); }
	CBLineThicknessPicker line_thickness;                   // line thickness picker dropdown - derived from CComboBox
	int line_thickness_index;
	CBLineStylePicker line_style;                           // line style picker dropdown - derived from CComboBox
	int line_style_index;
	CBPointStylePicker point_style;                         // point style picker dropdown - derived from CComboBox
	int point_style_index;
	int point_size;

	TArrayDynamic<Line_Data_Prop> lines;
	void RecomputeAllLineRanges()
	{
		for(int i=1; i<=lines.Length(); i++)
			lines(i).UpdateMinMax();
	}

	// **************************
	// *** Group: General Options
	// **************************
	int m_flag_use_autoupdate;                              // DDX: regular redraw ON/OFF
	void SetUpdateRate(double interval_s);
	double m_update_interval;                               // DDX: redraw interval in seconds
	int m_flag_autorescale;                                 // DDX: automatic scale to full size ON/OFF
	int m_flag_show_legend;                                 // DDX: Legend ON/OFF
	int m_flag_statusinfo;																	// DDX: Status Bar Info ON/OFF
	int m_flag_draw_every_nth;															// DDX: dont draw all points
	int m_draw_every_nth;                                   // DDX: draw every nth point
	int m_flag_mark_every_nth;															// DDX: dont mark all points
	int m_mark_every_nth;                                   // DDX: mark every nth point
	int m_flag_vertical_marker;															// DDX: draw a vertical marker
	int m_flag_upto_current;                                // DDX: draw only up to current
	int m_flag_fastdraw;																		// DDX: use fast draw routine - lower quality for unknown reasons...
	
	afx_msg void StartStopAutoUpdate(); 
	afx_msg void OnEnKillfocusEditUpdate(); 

	// *****************
	// *** Group: Export
	// *****************
	CString m_export_path;                                  // DDX: String containing the destination for plottool images
	CString m_export_filename;                              // DDX: String containing the file name for plottool images

	int m_pixels_horizontal;
	int m_pixels_vertical;

	int m_flag_export_jpg;
	int m_flag_export_png;
	int m_flag_export_bmp;
	int m_flag_export_emf;


  // **********************
	// *** calls from MenuBar
  // **********************
	void InitializeWithSensorTY(int sensor_number, CRect viewrect = CRect(0,0,0,0));                        // initializes the dialog with a single T/Y sensor already entered in the lines list
	void InitializeWithSensorXY(int sensor_number_x, int sensor_number_y, CRect viewrect = CRect(0,0,0,0)); // initializes the dialog with a single X/Y sensor already entered in the lines list
	void InitializeWithNSensors(TArray<int>& sensor_nrs, TArray<CRect>& viewrects);                         // initializes the dialog with a single x/y sensor already entered in the lines list

#ifdef zombies
	// ***************************************************************** 
	// *** GRAVEYARD - collect outdated code - disabled control elements
	// *****************************************************************


	// **************************
	// *** Group: Auxillary Lines
	// **************************
	afx_msg void OnBnClickedButtonAddAuxLine();             
  CString m_auxline_function;                             // DDX: string for aux function
	void AddAuxLine(int flag_isxy);                         // DOES NOT WORK YET !!! Adds an auxillary line to the List box ( Mathfunction )


	// ***************
	// *** TAB CONTROL 
	// ***************
	CTabCtrl m_tabcontrol;
	afx_msg void OnTcnSelchangeTabPlottool(NMHDR *pNMHDR, LRESULT *pResult);
	void DrawActiveTabElements();






	// **************************
	// *** Group: General Options
	// **************************
	afx_msg void OnBnClickedCheckAutorescale() { UpdateData(TRUE); }
	afx_msg void OnBnClickedCheckLegend() { UpdateData(TRUE); }
	afx_msg void OnBnClickedCheckDrawEveryNth() { UpdateData(TRUE); }
	afx_msg void OnEnKillfocusEditDrawEveryNth() { UpdateData(TRUE); }
	afx_msg void OnBnClickedCheckVertAtCurrent() { UpdateData(TRUE); }
	afx_msg void OnBnClickedCheckDrawUptoCurrent() { UpdateData(TRUE); }
	afx_msg void OnBnClickedCheckDrawFast() { UpdateData(TRUE); }

	afx_msg void OnBnClickedButtonExport() { ExportToFile(); }

	// *****************************      1  2
	// *** Group: Additional Buttons      3  4
	// *****************************      5  6  
	afx_msg void OnBnClickedButtonShow() { ShowGraph(); }    // show the Child dialog
	afx_msg void OnBnClickedButtonHide() { ShowWindow(SW_HIDE); }     // hide the Main dialog
	afx_msg void OnBnClickedButtonScaleFull() { ScaleToFull(); }
	afx_msg void OnBnClickedButtonRedraw() { Redraw(); }
	afx_msg void OnBnClickedButtonAxisEqual();
	afx_msg void OnBnClickedButtonPrint();     // call print routine
	afx_msg void OnBnClickedUpdate();          // manual trigger to read new data from file
	afx_msg void OnBnClickedReadoptions() { UpdateDialogData(); };     // extra button to grab changes in "edit hotint options"




// functions to read from file, filename is datamember of dialog
	int LoadFile_LoadFile2String(mystr& filename, mystr& str);	// loads entire file to a string
	int LoadFile_CountCommentLines(mystr& filestring);	// count number of comment lines in the file
	int LoadFile_CountColumns(mystr& filestring, int linenr = 0);	// counts number of columns in first data line
	int LoadFile_GetColumnNames(mystr& filestring, TArrayDynamic<mystr>& numbercount, int linenr = 0);	// extracts column names from last commented line
	int LoadFile_ReadColumnToArray(int colnr, TArray<double>* column, mystr& filestring, int& append_at_line);	// reads a column to array
	int LoadFile_IsSolFile(mystr& filestring); // identify as solution file
	int LoadFile_IsSolParFile(mystr& filestring); // identify as solution parameter file
	void NewSensorValues(double* vector, int len); 
#endif
// CONTEXTMENU
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	afx_msg void OnCmUpdateData();
	afx_msg void OnCmShowGraph();
	afx_msg void OnCmSaveSettings();
	afx_msg void OnCmMarkcurrent();
// END CONTEXTMENU
	afx_msg void OnBnClickedButtonExportOptions();
	//afx_msg void OnBnClickedOk();

	// Diagnose functions - Functions that trigger redraw
	virtual void SetCaller(mystr calleri);

protected:
	virtual BOOL OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult);
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class CMyPlotToolView:       derived View class 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// derived class for PlotToolView
class CMyPlotToolView : public CMyBase2DView
{
	DECLARE_DYNCREATE(CMyPlotToolView)

public:
	CMyPlotToolView(CWnd* pParent = NULL);       // Standardkonstruktor
	virtual ~CMyPlotToolView() {};
public:	
	virtual void Reset();                // initialize with default values ( from EDC or hardcoded )
	virtual void ErrorMessage(mystr& msg);

// ACCESS // LINKED DIALOGS & OBJECTS
private:
	CPlotToolDlg* m_pPTD;
public:
	void SetPlotToolDlg(CPlotToolDlg* pPTD) { m_pPTD = pPTD; } // must be set by 
	CPlotToolDlg* GetPlotToolDlg() { return m_pPTD; } // returns pointer to the Dialog containing all data

protected:
	DECLARE_MESSAGE_MAP()

public:
//COMMUNICATION WITH SETTINGs (EDC)
	virtual void LoadData();						 // load View-specific setting from options EDC
	virtual void WriteData();


public:
//WINDOW PROPERTIES & ADDITIONAL DISPLAY ELEMENTS
// CATCH MOUSE directly
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);          // zoom in and out 
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);                    // start drag/drop 
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point); // center view       // end drag/drop 
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);                      // value - tracking 
	afx_msg void OnNcMouseMove(UINT nHitTest, CPoint point);                  // ...

protected:
// CATCH OTHER via Msg
	virtual BOOL OnWndMsg(UINT message, WPARAM wParam, LPARAM lParam, LRESULT* pResult);

public:
	// DRAG - DROP:    
	void DragDropAction();																										// drag drop action - move and redraw
	void RectZoomAction();																										// zoom with a rectangle 
	void DrawZoomRect();



// THE DRAWING PROCESS
public:
	virtual void OnDraw(CDC* pDC);                          // this routine is automatically called by the framework (redrawwindow etc.)
private:
	virtual void DrawElements(CDC* pDC);                    // draws the elements directly on the screen 
// axis and grid
	virtual void DrawAxis(CDC* pDC);                        // axis and surrounding rectangle
	virtual void DrawTicksAndLabels(CDC* pDC);              // all labels and calls gridline
	virtual void DrawVertLines(CDC* pDC, TArray<double>x, int line_type );
	virtual void DrawHoriLines(CDC* pDC, TArray<double>y, int line_type );
	int grid_x_major_linetype, grid_x_minor_linetype;
	int grid_y_major_linetype, grid_y_minor_linetype;
// text
	virtual void DrawTitle(CDC* pDC);                       // graph title string
	virtual void DrawAxisLabels(CDC* pDC);                  // labels for x&y axis
// data
	virtual void DrawLines(CDC* pDC);                       // draws the data lines
  virtual void DrawDataPoint(CDC* pDC, int pointstyle, CPoint& centerpoint, CRect& markerrect, COLORREF color);
	virtual void DrawMarkerCurrent(CDC* pDC, CPoint datapoint, int thickness = 10);
	virtual void DrawLineVisible(CDC* pDC, CRect& plotrect, CPoint& start, CPoint& end, double oversize_percentage = 1.0); // draws section of the line that is IN the plot-rectangle when startpoint or endpoint are NOT
// legend
	virtual void DrawLegend(CDC* pDC);											// !!! UNDER CONSTRUCTION - NOT FINISHED !!!
	void GetLegendRect(CRect& legendrect)	{	legendrect = legend; }  
	CRect legend;

// is redraw allowed/prohibited ?
	virtual int DoNotRedraw(); 
	virtual void ForbidRedraw(); 
	virtual void AllowRedraw(); 

// StatusBarText
	virtual void UpdateStatusBarInfo();
	virtual void SetStatusBarText_Main(mystr text); 
	virtual void SetStatusBarText_X (mystr text);
	virtual void SetStatusBarText_Y (mystr text);

public:
// FUNCTIONS THAT PREPARE THE DRAWING DEVICE CONTENTS
	virtual void PrepareDC_Memory(CDC* pDC, CPrintInfo* pInfo = NULL);              // client is Rect with Pixels as output BMP - only in derived

public:
// REAL VALUES of COORDINATES
// full range of all lines to plot (real values)
	Box2D fullrange;                  
	double Get_XMin_glob() { return fullrange.PMin().X(); } //returns smallest x number in all lines that are plotted
	double Get_XMax_glob() { return fullrange.PMax().X(); } //returns largest x number in all lines that are plotted
	double Get_YMin_glob() { return fullrange.PMin().Y(); } //returns smallest y number in all lines that are plotted
	double Get_YMax_glob() { return fullrange.PMax().Y(); } //returns largest y number in all lines that are plotted

	virtual void ComputeFullRange(int flag_equal = 0); // compute limits of plot and the scaling factor
	virtual void ViewSetsRange(Box2D& range); // set range of plot - changes come from view
  virtual void SetShownRange(double xmin, double xmax, double ymin, double ymax); // sets range of plot
	virtual void ZoomToFullRange(int flag_equal = 0, int flag_redraw = 1); 	// sets the limits such that all lines are shown / optional equal size for x&y ticks / optional (do not) redraw

	// ticks
	virtual void ComputeMajorTicksArray(TArray<double>& ticks, double min, double max, int flag = 1, double scalingfactor = 1.); // computes the major ticks for the axis, the scaled distance is ((max-min)*scalingfactor)
	double GetStandardTickLength(double interval); // suggest suitable ticks for an interval (for linear scale)
	// TO BE DONE HERE:
	//double GetLogTickLength(double min, double max) {return 0.}; // suggest suitable ticks for an interval (for logarithmic scale)
typedef enum 
{
	PT_VERYFINE = 1,
	PT_FINE = 2,
	PT_NORMAL = 3,
	PT_THICK = 4,
	PT_VERYTHICK = 5
} TPenWidth;

public:
	// AUXILARY FUNCTIONS: LINE THICKNESSES, FONT SIZES, 
	int m_global_linethickness_factor;                      // scaling factor for all plotted lines
	int GetPenWidth(TPenWidth tag)                          // thickness in PTS for a selection in DropDownMenu
	{ 
		switch(tag)
		{
		case PT_VERYFINE: return 1;
		case PT_FINE: return 5;
		case PT_NORMAL: return 10;
		case PT_THICK: return 20;
		case PT_VERYTHICK: return 40;
		default: return 10;
		}
	}
	int GetFontSize(int fontsize) {return (int) (0.5 + fontsize * 6.0);} // font size in PTS
};


#pragma once


