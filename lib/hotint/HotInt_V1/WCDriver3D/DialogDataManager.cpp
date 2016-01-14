//#**************************************************************
//# filename:             DialogDataManager.cpp
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
#include "DialogDataManager.h"

//#include <iostream.h>
//#include <fstream.h>
//#include  <io.h>
//#include  <stdio.h>
//#include  <stdlib.h>

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <io.h>


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define WM_UPDATE (WM_USER + 2)
#define ANIMATE_TIMER_ID 1


//++++++++++++++++++++++++++++++++++++
//deleted elements in data manager:
//'Edit field' for submitting data: IDC_EDIT_DATA_REQUEST_COMMAND
//				DDX_Text(pDX, IDC_EDIT_DATA_REQUEST_COMMAND, m_strDataRequest);
//'Submit button': IDC_BUTTON_SUBMIT_PRINT_DATA_COMMAND
//                 ON_BN_CLICKED(IDC_BUTTON_SUBMIT_PRINT_DATA_COMMAND, OnButtonSubmitPrintDataCommand)
//
//Button 'retrieve data': IDC_BUTTON_RETRIEVE_DATA
//                        ON_BN_CLICKED(IDC_BUTTON_RETRIEVE_DATA, OnButtonRetrieveData)
//


/////////////////////////////////////////////////////////////////////////////
// CDialogDataManager dialog


CDialogDataManager::CDialogDataManager() : CDialog(CDialogDataManager::IDD, NULL),
bAnimating(false),
bRetrievingData(false)
, m_animation_delay(25)
{
	//{{AFX_DATA_INIT(CDialogDataManager)
	m_TimePointNumber = 0;
	m_strDataRequest = _T("");
	//}}AFX_DATA_INIT
}

void CDialogDataManager::Create(CWnd * pParent)
{
	CDialog::Create(IDD,pParent);
	CRect r;
	pParent->GetWindowRect(&r);
	CPoint p(r.right - 1,r.bottom - 1);

	LONG sx = (LONG)::GetSystemMetrics(SM_CXFULLSCREEN); 
	LONG sy = (LONG)::GetSystemMetrics(SM_CYFULLSCREEN); 

	if (p.x > sx-1) p.x = sx-1;
	if (p.y > sy-1) p.y = sy-1;

	GetWindowRect(&r);
	//	ClientToScreen(&r);
	r.top = p.y - r.Height();
	r.left = p.x - r.Width();
	r.right = p.x;
	r.bottom = p.y;
	MoveWindow(r,FALSE);
	ShowWindow(SW_SHOW);
}

void CDialogDataManager::CheckScreenLocation()
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


void CDialogDataManager::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	int MinTimePointNumber = 0;
	if(TimePointArray.GetSize())
		MinTimePointNumber = 1;
	int MaxTimePointNumber = TimePointArray.GetSize();

	if(::IsWindow(m_ScrollTimePoint.m_hWnd))
		m_ScrollTimePoint.SetScrollRange(MinTimePointNumber,MaxTimePointNumber);

	//{{AFX_DATA_MAP(CDialogDataManager)
	DDX_Control(pDX, IDC_SCROLLBAR_TIME_POINT, m_ScrollTimePoint);
	DDX_Text(pDX, IDC_EDIT_TIME_POINT_NUMBER, m_TimePointNumber);
	DDV_MinMaxInt(pDX, m_TimePointNumber, MinTimePointNumber, MaxTimePointNumber);
	//}}AFX_DATA_MAP

	if(!pDX->m_bSaveAndValidate)
		m_ScrollTimePoint.SetScrollPos(m_TimePointNumber);

	DDX_Text(pDX, IDC_EDIT_TIME_POINTS_SAVED, MaxTimePointNumber);

	long MemUsed = DataStorage::GetMemoryUsed()/1024;
	DDX_Text(pDX, IDC_EDIT_MEMORY_USED, MemUsed);

	CString CurrentTime;
	if(m_TimePointNumber)
	{
		CurrentTime.Format("%g", pWCDI->GetActualDrawTime());
	}
	DDX_Text(pDX, IDC_EDIT_CURRENT_TIME, CurrentTime);
	DDX_Text(pDX, IDC_EDIT_ANIMATION_DELAY, m_animation_delay);
	DDV_MinMaxDouble(pDX, m_animation_delay, 1., 1e8);
}


BEGIN_MESSAGE_MAP(CDialogDataManager, CDialog)
	//{{AFX_MSG_MAP(CDialogDataManager)
	ON_BN_CLICKED(IDC_BUTTON_ANIMATE, OnButtonAnimate)
	ON_BN_CLICKED(IDC_BUTTON_LOAD_FILE, OnButtonLoadFile)
	ON_BN_CLICKED(IDC_BUTTON_SAVE_FILE, OnButtonSaveFile)
	ON_EN_KILLFOCUS(IDC_EDIT_TIME_POINT_NUMBER, OnKillfocusEditTimePointNumber)
	ON_WM_CLOSE()
	ON_WM_HSCROLL()
	ON_WM_TIMER()
	ON_BN_CLICKED(IDC_BUTTON_SAVE_FILE_SPECIAL, OnButtonSaveFileSpecial)
	//}}AFX_MSG_MAP
	ON_MESSAGE(WM_UPDATE,OnUpdate)
	ON_EN_KILLFOCUS(IDC_EDIT_ANIMATION_DELAY, OnEnKillfocusEditAnimationDelay)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialogDataManager message handlers

BOOL CDialogDataManager::OnInitDialog() 
{
	CDialog::OnInitDialog();

	int MinTimePointNumber = 0;
	if(TimePointArray.GetSize())
		MinTimePointNumber = 1;
	int MaxTimePointNumber = TimePointArray.GetSize();

	m_ScrollTimePoint.SetScrollRange(MinTimePointNumber,MaxTimePointNumber);
	m_ScrollTimePoint.SetScrollPos(m_TimePointNumber);

	bAnimating = false;
	KillTimer(ANIMATE_TIMER_ID);

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void CDialogDataManager::OnButtonAnimate() 
{
	if(!TimePointArray.GetSize())
		return;
	if(bAnimating)
	{
		KillTimer(ANIMATE_TIMER_ID);
		bAnimating = false;
		return;
	}

	bAnimating = true;

	if (pWCDI->GetIOption(115))	m_TimePointNumber = 1;
	// AD: note - start with either first frame if the flag is specified  
	// b	115	ViewingOptions.Animation.animate_from_beginning	1	"1|(0) ... (Don't) start animation from beginning."

	SetTimer(ANIMATE_TIMER_ID,(int)m_animation_delay,NULL);
}

CString CDialogDataManager::GetHotintDataVersionString() const
{
	char str[256];
	sprintf_s(str, "HOTINTDataVersion %s",pWCDI->GetHotintVersion().GetString());
	//AfxMessageBox(str);

	return CString(str);
}

void CDialogDataManager::OnButtonSaveFile() 
{
	CDialogSaveSpecial dss;
	dss.m_FirstDataUnit = 1;
	dss.m_LastDataUnit = TimePointArray.GetSize();
	dss.m_SaveEachOf = 1;
	SaveFile(dss);
}

void CDialogDataManager::OnKillfocusEditTimePointNumber() 
{
	if(bAnimating)
		return;

	int OldTimePointNumber = m_TimePointNumber;
	UpdateData(TRUE);
	if(OldTimePointNumber != m_TimePointNumber)
		ActualizeTimePoint();
}

void CDialogDataManager::OnOK()
{
	UpdateData(TRUE);
	KillTimer(ANIMATE_TIMER_ID);
	DestroyWindow();
}

void CDialogDataManager::OnClose() 
{
	KillTimer(ANIMATE_TIMER_ID);
	DestroyWindow();
}

void CDialogDataManager::OnCancel() 
{
	KillTimer(ANIMATE_TIMER_ID);
	DestroyWindow();
}

void CDialogDataManager::AddEntry(DataStorage & ds, double m_TimePoint)
{
	bool bMoveToTheNewPosition = m_TimePointNumber == TimePointArray.GetSize();
	TimePointArray.Add(m_TimePoint);
	DataStorageArray.Add(ds);
	if(bMoveToTheNewPosition)
	{
		m_TimePointNumber = TimePointArray.GetSize();
	}
	if(::IsWindow(m_hWnd))
	{
		SendMessage(WM_UPDATE);
	}
}

void CDialogDataManager::RemoveAll()
{
	pWCDI->RemoveResults();  //PG new
	TimePointArray.RemoveAll();   //PG new
	DataStorageArray.RemoveAll();  //PG old

	m_TimePointNumber = TimePointArray.GetSize();
	if(::IsWindow(m_hWnd))
	{
		if(m_TimePointNumber)
			ActualizeTimePoint();
		else
			UpdateData(FALSE);
	}

}

LRESULT CDialogDataManager::OnUpdate(WPARAM, LPARAM)
{
	// Eat spurious WM_UPDATE messages
	MSG msg;
	while(::PeekMessage(&msg, m_hWnd, WM_UPDATE, WM_UPDATE, PM_REMOVE));

	UpdateData(FALSE);

	return 0;
}

void CDialogDataManager::ActualizeTimePoint()
{
	if(m_TimePointNumber > TimePointArray.GetSize() || (m_TimePointNumber < 1 && TimePointArray.GetSize() > 0))
	{
		//AfxMessageBox("Wrong time point number!",MB_OK | MB_ICONSTOP);	//$ DR 2013-05-24 removed because HOTINT hanged itself
		return;
	}

	
	int solset_sol_data_to_files;
	pWCDI->MBS_EDC_TreeGetInt(solset_sol_data_to_files, "SolverOptions.Solution.store_data_to_files");

	if (solset_sol_data_to_files)
	{
		//$ PG 2012-4-24: work-around... as long as DataStorage is part of arguments in LoadResults
		DataStorage ds;
		pWCDI->LoadResults((WCDInterface::DataLoader&) ds, m_TimePointNumber);
		
	}
	else
	{
		pWCDI->LoadResults((WCDInterface::DataLoader&)DataStorageArray[m_TimePointNumber - 1], m_TimePointNumber);
	}
	
	//!AD: Added update options and redraw all PlotTool Windows
	//pGLDrawWnd->SendMessage(WM_REDRAW);   //! AD 2013-02-08: BUGFIX no more double snapshots - INCLUDED IN WCDRIVERFUNCTION
	pWCDI->GetUserInterface()->CallWCDriverFunction(1);    // update all PlotTools

	UpdateData(FALSE);
}

void CDialogDataManager::SetScrollPosToLastTimePoint()
{
	m_TimePointNumber = TimePointArray.GetSize();
	m_ScrollTimePoint.SetScrollPos(m_TimePointNumber);
	ActualizeTimePoint();
}

void CDialogDataManager::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	if(bAnimating)
		return;

	int OldTimePointNumber = m_TimePointNumber;
	switch(nSBCode)
	{
	case 0:
		m_TimePointNumber--;
		break;
	case 1:
		m_TimePointNumber++;
		break;
	case 2:
		m_TimePointNumber -= TimePointArray.GetSize()/10 + 1;
		break;
	case 3:
		m_TimePointNumber += TimePointArray.GetSize()/10 + 1;
		break;
	case 5:
		m_TimePointNumber = nPos;
		break;
	}

	if(m_TimePointNumber < 1)
		m_TimePointNumber = 1;
	if(m_TimePointNumber > TimePointArray.GetSize())
		m_TimePointNumber = TimePointArray.GetSize();


	if(OldTimePointNumber != m_TimePointNumber)
		ActualizeTimePoint();

	CDialog::OnHScroll(nSBCode, nPos, pScrollBar);
}

void CDialogDataManager::OnTimer(UINT_PTR nIDEvent) 
{

	ActualizeTimePoint(); // this creates the file during redraw

	//AD: moved these lines down such that first frame is animated too ( increase AFTER first frame is drawn )
	m_TimePointNumber += pWCDI->GetIOption(109); 
	
	if(m_TimePointNumber >= TimePointArray.GetSize())
	{
		ActualizeTimePoint(); // always include last frame
		bAnimating = false;
		KillTimer(ANIMATE_TIMER_ID);
		return;
	}

	CDialog::OnTimer(nIDEvent);
}

void CDialogDataManager::OnButtonSubmitPrintDataCommand() 
{
	UpdateData(TRUE);
	pWCDI->PrintData(m_strDataRequest);
}

void CDialogDataManager::OnButtonRetrieveData() 
{
	SendMessage(WM_SETTEXT,NULL,(LPARAM)"Retrieving data - please wait");
	bRetrievingData = true;
	for(m_TimePointNumber = 1; m_TimePointNumber <= TimePointArray.GetSize(); m_TimePointNumber++)
		pWCDI->LoadResults((WCDInterface::DataLoader&)DataStorageArray[m_TimePointNumber - 1], m_TimePointNumber);
	bRetrievingData = false;
	SendMessage(WM_SETTEXT,NULL,(LPARAM)"Data operation is ready");

	m_TimePointNumber = TimePointArray.GetSize();
	ActualizeTimePoint();
}

void CDialogDataManager::OnButtonSaveFileSpecial() 
{
	CDialogSaveSpecial dss(TimePointArray.GetSize(),this);
	if(dss.DoModal() == IDCANCEL)
		return;
	SaveFile(dss);
}

void CDialogDataManager::SaveFile(CDialogSaveSpecial & info)
{
	// for the data file format specification see OnButtonLoadFile()

	static char BASED_CODE szFilter[] = "TXT file (*.txt)|*.txt|SOL file (*.sol)|*.sol|DAT file (*.dat)|*.dat|All Files (*.*)|*.*||";

	CFileDialog fd(FALSE,"txt",0,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ENABLESIZING,szFilter); //true=load, false=save
	if(fd.DoModal() == IDCANCEL) 
		return;

	//CFileDialog fd(FALSE);
	//fd.m_ofn.lpstrDefExt = "dat";
	//if(fd.DoModal() == IDCANCEL)
	//	return;

	SendMessage(WM_SETTEXT,NULL,(LPARAM)"Saving data - please wait");

	CString pathname = fd.GetPathName();

	int textmode = 0;
	if (pathname.GetLength() > 4 && ((pathname.Right(4).MakeLower() == CString(".txt")) || (pathname.Right(4).MakeLower() == CString(".sol"))))
	{
		textmode = 1;
		//AfxMessageBox("text-file!!!");
	}

	if (!textmode)
	{
		//.dat - archive mode
		CFile file;
		file.Open(pathname,CFile::modeWrite | CFile::modeCreate);

		CArchive ar(&file,CArchive::store);
		ar << CString(GetHotintDataVersionString()) <<
			(long)min(info.m_LastDataUnit,DataStorageArray.GetSize())/info.m_SaveEachOf;	//$ DR 2011-11-21: added typecast (long)
		for(	int i = max(info.m_FirstDataUnit - 1,0);
			i < min(info.m_LastDataUnit,DataStorageArray.GetSize());
			i += info.m_SaveEachOf)
		{
			DataStorageArray[i].Serialize(ar);
		}
	}
	else
	{
		//.txt mode
		//.dat - archive mode
		int opt = ios::out;

		ofstream fout;
		fout.open(pathname, opt); //with read sharing allowed!

		int rv;
		rv = _access(pathname, 0);
		if (_access(pathname, 0) == -1) 
		{
			AfxMessageBox("Error: the file could not be created.\nCheck if you have all permissions to write the file and\ncheck the filename for correctness!");
			SendMessage(WM_SETTEXT,NULL,(LPARAM)"Data NOT saved");
			return; 
		}

		fout.precision(17);
		int mbs_checksum1 = 0; //equal to total size of DOF loaded
		int mbs_checksum2 = 0; //checksum2 = 1*bodydof1+13*bodydof2+13^2*bodydof3+... - for later versions

		if (DataStorageArray.GetSize() != 0)
		{
			WCDInterface::DataLoader& data = (WCDInterface::DataLoader&)DataStorageArray[0];
			data.GetTime(); //this resets the current read position
			data >> mbs_checksum1;
		}

		int savedsteps = min(info.m_LastDataUnit-info.m_FirstDataUnit+1,DataStorageArray.GetSize())/info.m_SaveEachOf;

		fout << GetHotintDataVersionString() << "\n"; //version number
		fout << mbs_checksum1 << " " << mbs_checksum2 << "\n"; //checksum1=nDOF;
		fout << savedsteps << "\n"; //number of available steps

		int start = max(info.m_FirstDataUnit - 1,0);
		int end = min(info.m_LastDataUnit,DataStorageArray.GetSize());

		for(	int i = start;
			i < end;
			i += info.m_SaveEachOf)
		{
			WCDInterface::DataLoader& data = (WCDInterface::DataLoader&)DataStorageArray[i];
			int len;
			double x;

			//write data: time, len, data1, data2, .... \n
			fout << data.GetTime() << " "; //this resets the current read position

			data >> len;
			fout << len << " ";

			for (int j=1; j <= len; j++)
			{
				data >> x;
				fout << x << " ";
			}
			fout << "\n";
		}

	}
	SendMessage(WM_SETTEXT,NULL,(LPARAM)"Data saved");
}

void CDialogDataManager::OnButtonLoadFile() 
{
	int solset_sol_data_to_files;
	pWCDI->MBS_EDC_TreeGetInt(solset_sol_data_to_files, "SolverOptions.Solution.store_data_to_files");

	if (solset_sol_data_to_files)
	{
		int nTimePoints = pWCDI->ReadSolDataInfo();

		//$ PG 2012-4-20: TODO: still TimePointArray.GetSize() is used many times in the code to get the upper limit for a valid time point. formerly, this was done via DataStorageArray.GetSize(), which might be empty (if solset.sol_data_to_files == true). in future, TimePointArray (since this always consists of the entries [1,...,TimePointArray.GetSize()]) should be replaced by a single integer value, say int maxTimePoint.
		TimePointArray.RemoveAll();
		for(int i=1; i<=nTimePoints; i++)
		{
			TimePointArray.Add(i);
		}

		//$ PG 2012-4-24: fake - needed for the moment as long as LoadResults uses DataStorageArray as parameter, see in ActualizeTimePoint() the call pWCDI->LoadResults((WCDInterface::DataLoader&)DataStorageArray[m_TimePointNumber - 1], m_TimePointNumber);
		//DataStorageArray.RemoveAll();
		//DataStorage ds;
		//for(int i=1; i<=nTimePoints; i++)
		//{
		//	DataStorageArray.Add(ds);
		//}
		//$ PG 2012-4-20:]

		m_TimePointNumber = TimePointArray.GetSize();
		ActualizeTimePoint();
	}
	else
	{
		// data file format:
		// row1: text identifying the version of the data format: HOTINTDataVersionX.YYY
		// row2: totalDOF checksum2
		// row3: total number of records
		// row(3+i): record i: time size_of_record data1 data2 ...

		static char BASED_CODE szFilter[] = "TXT file (*.txt)|*.txt|DAT file (*.dat)|*.dat|All Files (*.*)|*.*||";

		CFileDialog fd(TRUE,"txt",0,OFN_ENABLESIZING,szFilter); //true=load, false=save

		//fd.m_ofn.lpstrDefExt = "dat";
		if(fd.DoModal() == IDCANCEL)
			return;


		CString pathname = fd.GetPathName();

		int textmode = 0;
		if (pathname.GetLength() > 4 && ((pathname.Right(4).MakeLower() == CString(".txt")) || (pathname.Right(4).MakeLower() == CString(".sol"))))
		{
			textmode = 1;
		}

		if (!textmode)
		{
			//archive mode
			CFile file;
			file.Open(fd.GetPathName(),CFile::modeRead);

			SendMessage(WM_SETTEXT,NULL,(LPARAM)"Loading data - please wait");

			CArchive ar(&file,CArchive::load);
			CString DataFormatVersionIdentifier;
			ar >> DataFormatVersionIdentifier;

			if(DataFormatVersionIdentifier != GetHotintDataVersionString())
			{
				CString msg = "Error: data format version mismatch.\nfile version: ";
				msg += DataFormatVersionIdentifier;
				msg += "\npresent version: ";
				msg += CString(GetHotintDataVersionString());
				AfxMessageBox(msg,MB_OK | MB_ICONSTOP);
				return;
			}

			int nTimePoints;
			ar >> nTimePoints;
			DataStorageArray.RemoveAll();
			for(int i = 0; i < nTimePoints; i++)
			{
				DataStorage ds;
				ds.Serialize(ar);
				DataStorageArray.Add(ds);
			}

			//$ PG 2012-4-20: TODO: still TimePointArray.GetSize() is used many times in the code to get the upper limit for a valid time point. formerly, this was done via DataStorageArray.GetSize(), which might be empty (if solset.sol_data_to_files == true). in future, TimePointArray (since this always consists of the entries [1,...,TimePointArray.GetSize()]) should be replaced by a single integer value, say int maxTimePoint.
			TimePointArray.RemoveAll();
			for(int i=1; i<=nTimePoints; i++)
			{
				TimePointArray.Add(i);
			}
			//$ PG 2012-4-20:]

			if(m_TimePointNumber = DataStorageArray.GetSize())
				ActualizeTimePoint();
			else
				UpdateData(FALSE);
		}
		else
		{
			//text mode:

#ifdef my_new_stdiostream
			int opt = ios::in;
			ifstream storein(pathname, opt);
#else
			int opt = ios::in|ios::nocreate;
			ifstream storein(pathname, opt, filebuf::sh_read);
#endif

			if (!(storein.good()&&(!storein.fail())))
			{
				AfxMessageBox("ERROR: could not open stored solution!");
			}
			else
			{
				//read header:
				CString header = "";

				int endit = 0; char ch;
				while(!endit)
				{
					if (storein.eof()) endit = 1;
					else
					{
						storein.get(ch);
						if (ch == (char)10 || ch == (char)13 || ch == (char)12 || ch == (char)11 || ch == '\n' || ch == '\t') //$ DR 2013-12-11 removed "|| ch == ' ' " from this list, because heades now has this format "HOTINTDataVersion 1.2.34"
							endit = 1;
						else {header += ch;}
					}
				}


				if (header.GetLength() > 17 && header.Left(17) == CString("HOTINTDataVersion"))
				{
					SendMessage(WM_SETTEXT,NULL,(LPARAM)"Loading data - please wait");

					CString version_str = header.Right(header.GetLength()-17);
					//AfxMessageBox(CString("Version number = ")+version);

					//HotintVersionInfo data_version(atof(version_str));
					mystr version_mystr(version_str);
					HotintVersionInfo data_version(version_mystr);

					if (data_version > pWCDI->GetHotintVersion())
					//if (data_version.GetDoubleValue() > pWCDI->GetHotintVersion().GetDoubleValue())	//$ DR 2013-12-11
					{
						CString msg = "Error: data format version mismatch.\nfile version: ";
						msg += "Hotint Data Version: ";
						msg += data_version.GetString();
						msg += "\npresent Hotint version: ";
						msg += CString(pWCDI->GetHotintVersion().GetString());
						AfxMessageBox(msg,MB_OK | MB_ICONSTOP);
					}

					int ndof, checksum2;
					int nTimePoints;
					double dnum;

					storein >> dnum; ndof = (int)dnum;
					storein >> dnum; checksum2 = (int)dnum; //actually not used!
					storein >> dnum; nTimePoints = (int)dnum;

					int nactdof = pWCDI->CallCompFunction(208); //size of solution vector (initial conditions) = SystemSize

					if (ndof > nactdof) 
					{
						AfxMessageBox("ERROR: Stored solution does have different size than actual model!\n Load appropriate model first!");
					}
					else
					{	//$ DR 2013-12-11 changed {}, such that this branch is only executed if model size is ok
						double time;
						int len;
						DataStorageArray.RemoveAll();
						int firstwarn=1;

						for(int i = 0; i < nTimePoints; i++)
						{
							DataStorage ds;

							storein >> time;
							storein >> dnum; len  = (int)dnum;;

							if (len < nactdof && firstwarn) 
							{
								AfxMessageBox("WARNING: Stored solution does have different size than actual model!\n HOTINT will fill up the solution vector with zeros.");
								firstwarn = 0;
							}

							//ds.Serialize(ar);
							//ds.AllocateNewData(sizeof(double)*(len+2+1)); //one in reserve ...
							ds.AllocateNewData(sizeof(double)*(nactdof+2+1)); //$ DR 2013-12-12 changed len to nactdof, to have enough space for zeros
							WCDInterface::DataSaver& data = (WCDInterface::DataSaver&)ds;

							data.SetTime(time);
							data << len;

							double d;
							for (int j = 1; j <= len; j++)
							{
								storein >> d;
								data << d;
							}
							
							for (int j = 1; j <= nactdof-len; j++)
							{
								data << (int)0;		//$ DR 2013-12-12 fill up with zeros, if stored solution vector is to short.
							}

							data << (int)0;		//$!DR 2011-11-22 otherwise TimeInt::LoadResults tries to load additional data
							DataStorageArray.Add(ds);

						}

						//$ PG 2012-4-20: TODO: still TimePointArray.GetSize() is used many times in the code to get the upper limit for a valid time point. formerly, this was done via DataStorageArray.GetSize(), which might be empty (if solset.sol_data_to_files == true). in future, TimePointArray (since this always consists of the entries [1,...,TimePointArray.GetSize()]) should be replaced by a single integer value, say int maxTimePoint.
						TimePointArray.RemoveAll();
						for(int i=1; i<=nTimePoints; i++)
						{
							TimePointArray.Add(i);
						}
						//$ PG 2012-4-20:]

						if(m_TimePointNumber = DataStorageArray.GetSize())
							ActualizeTimePoint();
						else
							UpdateData(FALSE);
					} // end of correct size
				}				
				else
				{
					AfxMessageBox("ERROR: No HOTINT header found in data file!");
				}
			}
		}

		SendMessage(WM_SETTEXT,NULL,(LPARAM)"Data loaded");
	}
}



bool CDialogDataManager::RedrawNewResults()
{
	if(m_TimePointNumber == TimePointArray.GetSize())
		return true;
	return false;
}
void CDialogDataManager::OnEnKillfocusEditAnimationDelay()
{
	UpdateData(TRUE);

}
