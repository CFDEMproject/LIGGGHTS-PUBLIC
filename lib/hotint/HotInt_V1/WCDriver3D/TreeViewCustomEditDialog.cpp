//#**************************************************************
//# filename:             TreeViewCustomEditDialog.cpp
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
#include "CustomEditDialog.h"
#include "mystring.h"
#include "TreeViewCustomEditDialog.h"

#define IDC_DELETE_BUTTON (CEid_offset1-1)
#define IDC_BACKGROUND (CEid_offset1-2)
const int TreeViewSize = 180;

IMPLEMENT_DYNAMIC(CStaticBackground, CStatic)
BEGIN_MESSAGE_MAP(CStaticBackground, CStatic)
	ON_NOTIFY_EX( TTN_NEEDTEXTA, 0, OnToolTipNotify )
	ON_NOTIFY_EX( TTN_NEEDTEXTW, 0, OnToolTipNotify )
END_MESSAGE_MAP()

BOOL CStaticBackground::OnToolTipNotify(UINT id, NMHDR *pNMHDR, LRESULT *pResult)
{
	CTreeViewCustomEditDialog * pTreeViewCustomEditDialog = dynamic_cast<CTreeViewCustomEditDialog*>(GetParent());

	// need to handle both ANSI and UNICODE versions of the message
	TOOLTIPTEXTA* pTTTA = (TOOLTIPTEXTA*)pNMHDR;
	TOOLTIPTEXTW* pTTTW = (TOOLTIPTEXTW*)pNMHDR;
	char * strTipText = NULL;
	UINT nID = pNMHDR->idFrom;
	if (pNMHDR->code == TTN_NEEDTEXTA && (pTTTA->uFlags & TTF_IDISHWND) ||
		pNMHDR->code == TTN_NEEDTEXTW && (pTTTW->uFlags & TTF_IDISHWND))
	{
		// idFrom is actually the HWND of the tool
		nID = ::GetDlgCtrlID((HWND)nID);
	}

	if (nID != 0) // will be zero on a separator
	{
		if (pNMHDR->code == TTN_NEEDTEXTA)
			pTTTA->lpszText = pTreeViewCustomEditDialog->FindToolTipText(nID);
		else
			pTTTW->lpszText = (LPWSTR)pTreeViewCustomEditDialog->FindToolTipText(nID);
	}

	*pResult = 0;

	return TRUE;    // message was handled
}


// CTreeViewCustomEditDialog

IMPLEMENT_DYNAMIC(CTreeViewCustomEditDialog, CDialog)

CTreeViewCustomEditDialog::CTreeViewCustomEditDialog(CWnd* pParent /*=NULL*/)
: CDialog(CTreeViewCustomEditDialog::IDD, pParent)
{
	dialogname = "Edit properties";
	add_deletebutton = 0;
	deleteflag = 0;
	useapply = 0;
	callcompfn_action = 0;
	callcompfn_option = 0;
	callcompfn_value = 0;
	nVerticalScroll = 0;
	flagDialogInitialized = false;
}

CTreeViewCustomEditDialog::~CTreeViewCustomEditDialog()
{
	for (int i=1; i <= treeItemsData.Length(); i++)
		if(treeItemsData(i).path != 0) delete treeItemsData(i).path;
	Delete();
}

void CTreeViewCustomEditDialog::Delete()
{
	for (int i=1; i <= cedit.Length(); i++)
		if (cedit(i) != 0) delete cedit(i);
	for (int i=1; i <= cstatic.Length(); i++)
		if (cstatic(i) != 0) delete cstatic(i);
	for (int i=1; i <= cbutton.Length(); i++)
		if (cbutton(i) != 0) delete cbutton(i);
	for (int i=1; i <= createdControlsData.Length(); i++)
	{
		if(createdControlsData(i).path != 0) delete createdControlsData(i).path;
		if(createdControlsData(i).toolTipText != 0) delete createdControlsData(i).toolTipText;
	}
	cedit.SetLen(0);
	cstatic.SetLen(0);
	cbutton.SetLen(0);
	createdControlsData.SetLen(0);
	groupbutton.SetLen(0);
}

void CTreeViewCustomEditDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TREECUSTOMEDIT, m_treectrl_control);
}


BEGIN_MESSAGE_MAP(CTreeViewCustomEditDialog, CDialog)
	ON_BN_CLICKED(IDOK, &CTreeViewCustomEditDialog::OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, &CTreeViewCustomEditDialog::OnBnClickedCancel)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREECUSTOMEDIT, &CTreeViewCustomEditDialog::OnTvnSelchangedTreecustomedit)
	ON_BN_CLICKED(IDC_APPLY, &CTreeViewCustomEditDialog::OnBnClickedApply)
	ON_WM_VSCROLL()
	ON_WM_SIZE()
	ON_WM_CONTEXTMENU()
	ON_COMMAND(ID_CM_CLIPBOARD, &CTreeViewCustomEditDialog::OnCopyToClipboard)

END_MESSAGE_MAP()

 

// CTreeViewCustomEditDialog-Meldungshandler

void CTreeViewCustomEditDialog::Create(CWnd * pParent)
{
	//not called when DoModal
	CDialog::Create(IDD,pParent);

	//position and size window:
	CRect r;
	pParent->GetWindowRect(&r);
	CPoint p(r.left,r.bottom - 5);
	GetWindowRect(&r);
	p.x += r.Width();

	r.top = p.y - r.Height();
	r.left = p.x - r.Width();
	r.right = p.x;
	r.bottom = p.y;
	MoveWindow(r,FALSE);
	ShowWindow(SW_SHOW);

	RedrawWindow();

	UpdateData(False);

}

void CTreeViewCustomEditDialog::OnBnClickedOk()
{
	if (useapply)
	{
		OnBnClickedApply();
	}
	else
	{
		WriteData();
	}
	OnOK();
}

void CTreeViewCustomEditDialog::OnBnClickedCancel()
{
	OnCancel();
}

void CTreeViewCustomEditDialog::WriteData()
{
	ElementDataContainer* edc = GetCurrentEDC();
	if (!edc) return;

	//fill in element data:
	int ceID = CEid_offset2;
	//char str[256];
	CString cstr;

	ElementDataContainer edc_tmp;									//$ DR 2013-06-05
	pWCDI->CallCompFunction(211, 0, 0, &edc_tmp); //get model data
	ElementDataContainer* edc_glob_var = edc_tmp.GetCopy();
	const ElementData* ed_glob_var = NULL;

	TArray<double> doublelist;

	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);

		if (ed.IsBool())
		{
			if (ed.IsGroup())
			{
				ed.SetBool(cbutton(groupbutton(ed.GetGroup()))->IsDlgButtonChecked(CEid_offset1+i));
			}
			else
			{
				ed.SetBool(wndStaticBackground.IsDlgButtonChecked(CEid_offset1+i));
			}
		}
		else if (ed.IsCompAction() || ed.IsWCDriverAction())
		{
			//do nothing!
		}
		else
		{
			int len = 1;
			if (ed.IsVector()) len = ed.GetVectorLen();

			if (ed.IsDouble() || ed.IsVector() || ed.IsMatrix() || ed.IsInt() || ed.IsText())
			{
				int l = len;
				int editmode = 0;
				if (l == 0 || l >= 4 || ed.IsMatrix() || ed.IsVariableLength()) 
				{
					l = 1;
					editmode = 1;
				}

				for (int j=1; j <= l; j++)
				{
					wndStaticBackground.GetDlgItemText(ceID, cstr);
					if (!ed.IsLocked())
					{
						if (ed.IsText()) 
						{
							ed.SetText((LPCTSTR)cstr);
						}
						else
						{
							//$ DR 2013-06-05 check for global variables
							int found = 0;
							ElementData* ed_tmp = edc_glob_var->TreeFind(cstr);
							double val_tmp=0;
							if(ed_tmp)
							{
								if(!((ed.IsInt() || ed.IsDouble() || ed.IsVectorXYZ()) && (ed_tmp->IsInt() || ed_tmp->IsDouble())))
								{
									mystr warning = mystr("WARNING: variable '") +mystr(ed.GetDataName())+mystr("' has been found, but is of wrong type and therefore ignored. \n") ;
									pWCDI->GetUserInterface()->AddText(warning);
									found = 0; // ignore this element data
								}
								else
								{
									if(ed_tmp->IsInt()) val_tmp = ed_tmp->GetInt();
									else val_tmp = ed_tmp->GetDouble();

									found = 1;
								}
							}
							//else { no warning if variable not found, because called too often}

							if (ed.IsInt()) 
							{
								if(found)	ed.SetInt(val_tmp);
								else		ed.SetInt(atoi((LPCTSTR)cstr));
							}
							else if (ed.IsDouble()) 
							{
								if(found)	ed.SetDouble(val_tmp);
								else		ed.SetDouble(atof((LPCTSTR)cstr));
							}
							else //Vector or matrix
							{
								if (!editmode)
								{
									if(found)	ed.SetVectorVal(j, val_tmp);
									else		ed.SetVectorVal(j, atof((LPCTSTR)cstr));
								}
								else if (ed.IsVector())
								{
									//+++++++++++++++++++++++++++++++++++++++++++++
									mystr mstr = (LPCTSTR)cstr;
									doublelist.SetLen(0);

									//$ DR 2013-06-05:[ added global variables for complete vector
									ed_tmp = edc_glob_var->TreeFind(mstr);
									found = 0;
									if(ed_tmp)
									{
										if(ed_tmp->IsVector()) found = 1;
									}
									if(found)
									{
										ed.SetVectorLen(ed_tmp->GetVectorLen());
										for(int jj=1; jj<=ed_tmp->GetVectorLen(); jj++)
										{
											ed.SetVectorVal(jj,ed_tmp->GetVectorVal(jj));
										}
									}		//$ DR 2013-06-05:] added global variables for complete vector
									else
									{
										int end = 0;
										int endpos;

										while (!end)
										{
											endpos = mstr.Find(',');
											if (endpos != -1 && endpos != 0 && endpos <= mstr.Length()-1)
											{
												mystr ss = mstr.SubString(0, endpos-1);
												//$ DR 2013-06-05:[ added global variables for single components
												ed_tmp = edc_glob_var->TreeFind(ss.c_str());
												found = 0;
												if(ed_tmp)
												{
													if(ed_tmp->IsInt()) 
													{
														val_tmp = ed_tmp->GetInt();
														found = 1;
													}
													else if(ed_tmp->IsDouble()) 
													{
														val_tmp = ed_tmp->GetDouble();
														found = 1;
													}
												}
												if(found)	doublelist.Add(val_tmp);
												else 			doublelist.Add(atof(ss.c_str()));
												//$ DR 2013-06-05:] added global variables
												mstr = mstr.SubString(endpos+1, mstr.Length()-1);
											}
											else 
											{
												mstr.EraseSpaces();
												if (mstr.Length() != 0) //last value
												{			
													//$ DR 2013-06-05:[ added global variables for single components
													ed_tmp = edc_glob_var->TreeFind(mstr.c_str());
													found = 0;
													if(ed_tmp)
													{
														if(ed_tmp->IsInt()) 
														{
															val_tmp = ed_tmp->GetInt();
															found = 1;
														}
														else if(ed_tmp->IsDouble()) 
														{
															val_tmp = ed_tmp->GetDouble();
															found = 1;
														}
													}
													if(found)	doublelist.Add(val_tmp);
													else 			doublelist.Add(atof(mstr.c_str())); 
													//$ DR 2013-06-05:] added global variables
												}
												end = 1;
											}
										}

										ed.SetVectorLen(doublelist.Length());
										for (int jj=1; jj <= doublelist.Length(); jj++)
										{
											ed.SetVectorVal(jj, doublelist(jj));
										}
									}
								}
								else if (ed.IsMatrix())
								{
									//$ DR 2013-06-05:[ added global variables for complete matrix
									ed_tmp = edc_glob_var->TreeFind(cstr);
									found = 0;
									if(ed_tmp)
									{
										if(ed_tmp->IsMatrix()) found = 1;
									}
									if(found)
									{
										ed.SetMatrixSize(ed_tmp->GetMatrixRows(),ed_tmp->GetMatrixCols());
										for(int jj=1; jj<=ed_tmp->GetMatrixRows(); jj++)
										{
											for(int kk=1; kk<=ed_tmp->GetMatrixCols(); kk++)
											{
												ed.SetMatrixVal(jj,kk,ed_tmp->GetMatrixVal(jj,kk));
											}
										}
									}		//$ DR 2013-06-05:] added global variables for complete matrix
									else
									{

										//+++++++++++++++++++++++++++++++++++++++++++++
										mystr mstr = (LPCTSTR)cstr+mystr((char)10)+mystr((char)13);
										doublelist.SetLen(0);

										//AfxMessageBox(mstr.c_str());
										int end = 0;
										int endpos, endpos2;
										int maxcol = 0;
										int error = 0;
										int nrows = 0;
										char str[258];
										int limit = 256;
										int startpos = 0;

										while (!end)
										{
											int emptystring = 1;
											int ncol = 0;
											int endcol = 0;

											while (!endcol)
											{
												endpos = mstr.Find(startpos,',');
												endpos2 = mstr.Find(startpos,(char)10);

												if (endpos2 < endpos || endpos == -1) 
												{
													endpos = endpos2;
													endcol = 1;
												}

												if (endpos != -1 && endpos != 0 && startpos < endpos && endpos+1 <= mstr.Length()-1 && mstr.Length() != 0)
												{
													int len = endpos-startpos;
													mstr.CopySubStringNoSpaces(str, startpos, endpos-1, limit);

													if (strlen(str) != 0)
													{
														//$ DR 2013-06-05:[ added global variables for single components
														ed_tmp = edc_glob_var->TreeFind(str);
														found = 0;
														if(ed_tmp)
														{
															if(ed_tmp->IsInt()) 
															{
																val_tmp = ed_tmp->GetInt();
																found = 1;
															}
															else if(ed_tmp->IsDouble()) 
															{
																val_tmp = ed_tmp->GetDouble();
																found = 1;
															}
														}
														if(found)	doublelist.Add(val_tmp);
														else 			doublelist.Add(atof(str)); 
														//$ DR 2013-06-05:] added global variables
														//doublelist.Add(atof(str));
														ncol++;

														startpos = endpos+1;
														emptystring = 0;
													}
													else
													{
														startpos = endpos+1;
														doublelist.Add(0);
													}
												}
												else 
												{
													end = 1;
													endcol = 1;
												}
											}
											if (!end && !emptystring) nrows++;
											if (ncol != maxcol) 
											{
												if (maxcol != 0 && ncol != 0)
												{
													error = 1;
												}

												if (ncol > maxcol) maxcol = ncol;
											}
										}
										if (error) AfxMessageBox("Error in Edit-field: number of\ncolumns must be equal in each row!");


										ed.SetMatrixSize(nrows, maxcol);
										for (int j1=1; j1 <= nrows; j1++)
										{
											for (int j2=1; j2 <= maxcol; j2++)
											{
												ed.SetMatrixVal(j1, j2, 0);
												if ((j1-1)*maxcol+j2 <= doublelist.Length())
												{
													ed.SetMatrixVal(j1, j2, doublelist((j1-1)*maxcol+j2));
												}
											}
										}									
									}	// end of reading each matrix component
								}// end of matrix
							}	// end of vector or matrix
						} // end of "not text"
					} // end of "not locked"
					ceID++;
				}
			}
		}
	}
}



BOOL CTreeViewCustomEditDialog::OnInitDialog()
{
	CDialog::OnInitDialog();

	flagDialogInitialized = true;

	//set dialog name
	SetWindowText(dialogname);

	wndStaticBackground.Create("", WS_CHILD | WS_VISIBLE /*| WS_BORDER*/, CRect(1,1,5,5), this, IDC_BACKGROUND);
	wndStaticBackground.EnableToolTips();

	//add delete button?
	if (add_deletebutton)
	{
		MyCButton* cb = new MyCButton();
		cbutton.Add(cb);

		cb->Create("Delete", BS_CENTER|BS_VCENTER|WS_CHILD|WS_VISIBLE, 
			CRect(1,1,5,5), this, IDC_DELETE_BUTTON);
		cb->SetCustomEditDialog(this, 0); //tell button where to call if button pressed, 0=delete button!
		cb->SetDlgCtrlID(IDC_DELETE_BUTTON);
		cb->SetFont(GetFont());
	}

	currentTreeItemData.edc = NULL;

	TreeItemData tidStart;
	tidStart.hti = TVI_ROOT;
	tidStart.path = "";
	tidStart.edc = GetRootEDC();
	InsertTreeItem(tidStart);

	HTREEITEM ht = m_treectrl_control.GetFirstVisibleItem();
	m_treectrl_control.SelectItem(ht);
	m_treectrl_control.SetFocus();
	m_treectrl_control.Expand(ht,TVE_EXPAND);

	//expand one level in the hierarchy of items:
	HTREEITEM hTreeItem = m_treectrl_control.GetChildItem(ht);
	if(hTreeItem) 
	{
		do 
		{
			m_treectrl_control.Expand(hTreeItem,TVE_EXPAND);	
		} while( (hTreeItem = m_treectrl_control.GetNextSiblingItem(hTreeItem)) != NULL );
	} // end if
	m_treectrl_control.EnsureVisible(ht);

	PositionControls();

	return TRUE;  // return TRUE unless you set the focus to a control
}

void CTreeViewCustomEditDialog::InsertTreeItem(TreeItemData tid)
{
	// here we come with the parent item
	for (int i=1; i <= tid.edc->Length(); i++)
	{
		if (tid.edc->Get(i).IsEDC())
		{
			TreeItemData tidChild;
			char * name = tid.edc->Get(i).GetDataName();
			tidChild.hti = m_treectrl_control.InsertItem(name, tid.hti);
			tidChild.edc = tid.edc->Get(i).GetEDC();
			tidChild.path = new char[strlen(tid.path) + strlen(name) + 2];
			int shift = 0;
			if(tid.hti != TVI_ROOT)
			{
				// for the lower items the path of the parent needs to be added at first
				strcpy(tidChild.path, tid.path);
				tidChild.path[strlen(tid.path)] = '.';
				shift = strlen(tid.path) + 1;
			}
			strcpy(tidChild.path + shift, name);
			treeItemsData.Add(tidChild);
			InsertTreeItem(tidChild);
		}
	}
}


void CTreeViewCustomEditDialog::UpdateDialogItems()
{
	ClearDialogItems();

	ElementDataContainer* edc = GetCurrentEDC();
	if (!edc) return;

	//set size of properties window depending on number of edit fields:
	int ned = edc->Length();
	if (ned < 1) ned = 1;

	TArray <int> groupnum;
	TArray <int> grouplength;

	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);
		if (ed.IsGroup() && ed.GetGroup() > groupnum.Length())
		{
			//groupnum.SetLen(ed.GetGroup());
			groupnum(ed.GetGroup()) = 0;
		}
	}
	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);
		if (ed.IsGroup())
		{
			groupnum(ed.GetGroup())++;
			if (groupnum(ed.GetGroup()) > 1) ned--;
		}
	}
	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);
		if (ed.IsMatrix())
		{
			ned += 2; //add two extra lines for matrix!
		}
	}

	nVerticalSpaceForTheControls = ned * (CEedit_height + CEedit_vspace);

	grouplength = groupnum;

	//fill in element data:
	CFont* pfont = GetDlgItem(IDOK)->GetFont();

	int ceID = CEid_offset2;
	char str[256];
	char str2[256];

	for (int i=1; i <= groupnum.Length(); i++) {groupnum(i) = 0;}
	int yoff_reduce = 0;

	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);
		int px1 = 0;
		int px2 = CEstatictext_width + CEedit_hspace;
		int py = (i-1-yoff_reduce)*(CEedit_height + CEedit_vspace) - nVerticalScroll;
		// now we form the path string for the item
		char path[1000];
		strcpy(path,currentTreeItemData.path);
		path[strlen(currentTreeItemData.path)] = '.';
		strcpy(path + strlen(currentTreeItemData.path) + 1, ed.GetDataName());

		if (ed.IsBool() || ed.IsCompAction() || ed.IsWCDriverAction())
		{
			//add event handle for button pressed!!!!

			//cstatic.Add(0); //do not delete this static entry

			AddCreatedControlData(ed.GetToolTipText(), path, CEid_offset1+i);

			if (ed.IsBool())
			{
				CButton* cb = new CButton();
				cbutton.Add(cb);
				//checkbox:
				if (!ed.IsGroup())
				{
					cb->Create(ed.GetDataName(), BS_AUTOCHECKBOX|BS_LEFT|WS_CHILD|WS_VISIBLE, 
						CRect(px1,py,px1+CEstatictext_width,py+CEedit_height), &wndStaticBackground, CEid_offset1+i);
					cb->SetCheck(ed.GetBool());
				} 
				else
				{
					int ng = ed.GetGroup();

					int groupwidth = CEbutton_width + CEedit_hspace;
					//if (groupwidth*grouplength(ng) > max_window_width) groupwidth = max_window_width / grouplength(ng);

					CButton* group = 0;
					if (groupbutton.Length() >= ng) group = cbutton(groupbutton(ng));
					else
					{ //create new group for radio buttons:
						group = new CButton();
						groupbutton(ng) = cbutton.Add(group);
						group->Create("", BS_GROUPBOX|WS_CHILD|WS_VISIBLE, 
							CRect(px1,py-7,px1+CEstatictext_width + 3*(CEedit_hspace + CEedit_width),py+2+CEedit_height), &wndStaticBackground, CEid_offset_group+ng);
					}

					int x = groupnum(ng)*(groupwidth);
					cb->Create(ed.GetDataName(), BS_AUTORADIOBUTTON|BS_LEFT|WS_CHILD|WS_VISIBLE, 
						CRect(px1+x,10,px1+x+groupwidth-CEedit_hspace,9+CEedit_height-2), group, CEid_offset1+i);
					cb->SetCheck(ed.GetBool());

					groupnum(ng)++;
					if (groupnum(ng) > 1) {yoff_reduce++;}
				}
				cb->SetFont(pfont);
			}
			else // is action
			{
				MyCButton* cb = new MyCButton();
				cbutton.Add(cb);
				//button + Action!
				int x = 0, y = 0;
				if (ed.IsGroup())
				{
					int ng = ed.GetGroup();
					x = groupnum(ng)*(CEbutton_width+CEedit_hspace);
					groupnum(ng)++;
					if (groupnum(ng) > 1) {yoff_reduce++; y = -(CEedit_height + CEedit_vspace);}
				}
				cb->Create(ed.GetDataName(), BS_CENTER|BS_VCENTER|WS_CHILD|WS_VISIBLE, 
					CRect(px1+x,py+y,px1+x+CEbutton_width,py+y+CEedit_height), &wndStaticBackground, CEid_offset1+i);
				cb->SetCustomEditDialog(this, i); //tell button where to call if button pressed!
				cb->SetFont(pfont);


			}

		}
		else
		{
			//cbutton.Add(0); //do not delete

			CStatic* cs = new CStatic();
			cstatic.Add(cs);

			cs->Create(ed.GetDataName(), SS_SUNKEN|WS_CHILD|WS_VISIBLE, 
				CRect(px1,py,px1+CEstatictext_width,py+CEedit_height), &wndStaticBackground, CEid_offset1+i);
			cs->SetFont(pfont);

			int len;

			len = 1;
			if (ed.IsVector()) len = ed.GetVectorLen();

			if (ed.IsDouble() || ed.IsVector() || ed.IsMatrix() || ed.IsInt() || ed.IsText())
			{
				int ew = CEedit_width;
				int editmode = 0; //vector is: "val1, val2, val3, ..." or matrix is multiline comma-separated
				double l = len;

				if (ed.IsVector() && (l > 3 || ed.IsVariableLength()))
				{
					if (l == 0 || l >= 4 || ed.IsVariableLength())
					{
						ew = 3 * CEedit_width + 2* CEedit_hspace;
						editmode = 1;
					}
					else
						ew = (int)(3.*(double)ew/(double)l);
				}
				else if (ed.IsVector() && l == 0 && !ed.IsVariableLength())
				{
					AfxMessageBox("HOTINT System Error: ElementData with vector length=0 not allowed. Use ed.SetVariableLength() in generation of elementdata.");
				}
				if (ed.IsMatrix())
				{
					ew = 3 * CEedit_width + 2* CEedit_hspace;
					editmode = 1;
				}

				if (editmode) l = 1;

				for (int j=1; j <= l; j++)
				{
					MyEdit* ce = new MyEdit();
					cedit.Add(ce);

					AddCreatedControlData(ed.GetToolTipText(), path, ceID);

					if (ed.IsInt() || ed.IsDouble() || ed.IsVector())
					{
						double mi, ma;
						bool oneLimit;
						ed.GetMinMaxVal(mi, ma, oneLimit);
						ce->SetMinMax(mi, ma,oneLimit);
					}

					//possibly WS_THICKFRAME|
					int xfact = 1;
					if (ed.IsText()) xfact = 3;

					DWORD dwstyle = ES_AUTOHSCROLL|WS_TABSTOP|WS_CHILD|WS_VISIBLE|WS_BORDER;

					if (ed.IsLocked()) dwstyle = dwstyle|ES_READONLY;
					if (ed.IsInt()) dwstyle = dwstyle|ES_NUMBER; //only digits can be entered!!!!
					if (ed.IsMatrix()) dwstyle = dwstyle|ES_LEFT|ES_MULTILINE|ES_AUTOHSCROLL|ES_AUTOVSCROLL|ES_WANTRETURN|WS_VSCROLL; //only digits can be entered!!!!
					else dwstyle = dwstyle|ES_RIGHT;

					int height = CEedit_height;
					if (ed.IsMatrix()) height += 2*(CEedit_height + CEedit_vspace);
					ce->Create(dwstyle, CRect(px2,py,px2+ew*xfact+CEedit_hspace*(xfact-1),py+height), &wndStaticBackground, ceID);
					ce->SetFont(pfont);

					if (ed.IsText()) 
					{
						CString strtmp = ed.GetText();
						ce->SetWindowText(strtmp);
					}
					else
					{
						if (ed.IsInt()) 
						{
							sprintf(str, "%d", ed.GetInt());
							ce->SetWindowText(str);
						}
						else if (!editmode)
						{
							sprintf(str, "%.16g", ed.GetVectorVal(j));
							sprintf(str2, "%.14g", ed.GetVectorVal(j));

							if (strlen(str) > strlen(str2)+4)
								sprintf(str, "%.14g", ed.GetVectorVal(j));

							ce->SetWindowText(str);
						}
						else if (editmode && ed.IsVector())
						{
							//fill in comma-separated vector:
							mystr mstr = "";

							for (int jj=1; jj <= len; jj++)
							{
								sprintf(str, "%.16g", ed.GetVectorVal(jj));
								sprintf(str2, "%.14g", ed.GetVectorVal(jj));

								if (strlen(str) > strlen(str2)+4)
									sprintf(str, "%.14g", ed.GetVectorVal(jj));

								mstr += str;
								if (jj != len) mstr += mystr(", ");
							}
							ce->SetWindowText(mstr.c_str());
						}
						else if (editmode && ed.IsMatrix())
						{
							yoff_reduce -= 2; //3 lines for a matrix
							//fill in comma-separated matrix in multi-line edit dialog!
							CString mstr = "";
							//mystr mstr = "";
							mystr endline = mystr((char)13) + mystr((char)10);

							int rows = ed.GetMatrixRows();
							int cols = ed.GetMatrixCols();
							for (int j1=1; j1 <= rows; j1++)
							{
								for (int j2=1; j2 <= cols; j2++)
								{
									double val = ed.GetMatrixVal(j1, j2);

									sprintf(str, "%.16g", val);
									sprintf(str2, "%.14g", val);

									if (strlen(str) > strlen(str2)+6)
									{
										if (j2 != cols)
											sprintf(str, "%.14g, ", val);
										else
											sprintf(str, "%.14g", val);
									}
									else
									{
										if (j2 != cols)
											sprintf(str, "%.16g, ", val);
										else
											sprintf(str, "%.16g", val);
									}
									mstr += str;
								}
								if (j1 != rows) mstr += endline.c_str();
							}
							ce->SetWindowText(mstr);
						}
					}
					ceID++;
					px2 +=  ew + CEedit_hspace;
				}
			}
		}
	}
	RedrawWindow();
}

char * CTreeViewCustomEditDialog::FindToolTipText(int itemID)
{
	for (int i=1; i <= createdControlsData.Length(); i++)
	{
		if (createdControlsData(i).id == itemID)
		{
			GetDlgItem(IDC_EDIT_PATH)->SetWindowText(createdControlsData(i).path);
			return createdControlsData(i).toolTipText;
		}
	}
	return NULL;
}

void CTreeViewCustomEditDialog::AddCreatedControlData(const char* toolTipText, const char * path, int itemID)
{
	CreatedControlData ccd;
	ccd.id = itemID;
	ccd.path = new char[strlen(path) + 1];
	strcpy(ccd.path,path);
	if(toolTipText != NULL)
	{
		ccd.toolTipText = new char[strlen(toolTipText) + 1];
		strcpy(ccd.toolTipText,toolTipText);
	}
	else
		ccd.toolTipText = NULL;
	createdControlsData.Add(ccd);
}

void CTreeViewCustomEditDialog::ClearDialogItems()
{
	for (int i=1; i <= cedit.Length(); i++)
	{
		if (cedit(i) != 0) cedit(i)->DestroyWindow();
	}
	for (int i=1; i <= cstatic.Length(); i++)
	{
		if (cstatic(i) != 0) cstatic(i)->DestroyWindow();
	}
	for (int i=1; i <= cbutton.Length(); i++)
	{
		if (cbutton(i) != 0) cbutton(i)->DestroyWindow();
	}

	RedrawWindow();

	Delete();
}

void CTreeViewCustomEditDialog::ButtonClicked(int edcnum)
{
	//char str[100];
	//sprintf(str, "Button clicked: %d", edcnum);
	//AfxMessageBox(str);
	ElementDataContainer* edc = GetCurrentEDC();
	if (!edc) return;

	int action;
	int elnum;
	int loadnum;

	if (edcnum != 0)
	{
		const ElementData& ed = edc->Get(edcnum);
		ed.GetWCDriverAction(action, loadnum, elnum);

		if (action == 2) //only for load
		{
			int flag = pWCDI->GetUserInterface()->CallWCDriverFunction(action, loadnum, elnum);
			if (flag) //delete load
			{
				edc->Delete(edcnum);
				UpdateDialogItems();
				//alternative: GetDlgItem(CEid_offset1 + edcnum)->DestroyWindow();
			}
		}
	}
	else
	{
		deleteflag = 1;//CEid_offset1 + edcnum;
		OnCancel();
		//ClearDialogItems();
		//CreateDialogItems();
	}
}



void CTreeViewCustomEditDialog::OnTvnSelchangedTreecustomedit(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	// TODO: Fügen Sie hier Ihren Kontrollbehandlungscode für die Benachrichtigung ein.

	//save changes to data!
	WriteData();

	nVerticalScroll = 0;

	// search for the selected tree item in the list
	HTREEITEM ht = m_treectrl_control.GetSelectedItem();
	currentTreeItemData.edc = NULL;
	for(int i = 1; i <= treeItemsData.Length(); i++)
		if(treeItemsData(i).hti == ht)
			currentTreeItemData = treeItemsData(i);

	UpdateDialogItems();

	UpdateScrollBarInfo();

	*pResult = 0;
}

void CTreeViewCustomEditDialog::OnBnClickedApply()
{
	if (useapply)
	{
		WriteData();
		pWCDI->CallCompFunction(callcompfn_action, callcompfn_option, callcompfn_value, edc); //set model data
		pGLDrawWnd->SendMessage(WM_REDRAW);
		//!AD: Added update options and redraw all PlotTool Windows
		pWCDI->GetUserInterface()->CallWCDriverFunction(1,1);    // update all PlotTools - update options from hotint options and redraw
	}
	else
	{
		OnBnClickedOk();
	}
}

void CTreeViewCustomEditDialog::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	int nDelta;
	int nMaxPos = nVerticalSpaceForTheControls - nCurrentBackgroundHeight;

	switch (nSBCode)
	{
	case SB_LINEDOWN:
		if (nVerticalScroll >= nMaxPos)
			return;
		nDelta = min(nMaxPos / 20, nMaxPos - nVerticalScroll);
		break;
	case SB_LINEUP:
		if (nVerticalScroll <= 0)
			return;
		nDelta = -min(nMaxPos / 20, nVerticalScroll);
		break;
  case SB_PAGEDOWN:
		if (nVerticalScroll >= nMaxPos)
			return;
		nDelta = min(nMaxPos / 5, nMaxPos - nVerticalScroll);
		break;
	case SB_PAGEUP:
		if (nVerticalScroll <= 0)
			return;
		nDelta = -min(nMaxPos / 5, nVerticalScroll);
		break;
	case SB_THUMBPOSITION:
	case SB_THUMBTRACK:
		nDelta = (int)nPos - nVerticalScroll;
		break;
	default:
		return;
	}
	nVerticalScroll += nDelta;
	SetScrollPos(SB_VERT,nVerticalScroll,TRUE);

	WriteData();
	SetRedraw(FALSE);
	UpdateDialogItems();
	SetRedraw(TRUE);
	RedrawWindow();
}

void CTreeViewCustomEditDialog::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);

	if(flagDialogInitialized)
		PositionControls();
}

void CTreeViewCustomEditDialog::PositionControls()
{
	CRect rectClient;
	GetClientRect(rectClient);

	// first we place the buttons in the bottom-right corner
	CRect rectButton;
	GetDlgItem(IDCANCEL)->GetWindowRect(rectButton);
	rectButton.MoveToXY(rectClient.BottomRight() - rectButton.Size() - CSize(CEedit_xoff, CEedit_yoff));
	GetDlgItem(IDCANCEL)->MoveWindow(rectButton);
	rectButton.MoveToX(rectButton.left - rectButton.Width() - CEedit_xoff);
	GetDlgItem(IDC_APPLY)->MoveWindow(rectButton);
	rectButton.MoveToX(rectButton.left - rectButton.Width() - CEedit_xoff);
	GetDlgItem(IDOK)->MoveWindow(rectButton);
	if(add_deletebutton)
	{
		rectButton.MoveToX(rectButton.left - rectButton.Width() - CEedit_xoff);
		GetDlgItem(IDC_DELETE_BUTTON)->MoveWindow(rectButton);
	}

	// the edit box for the path
	GetDlgItem(IDC_EDIT_PATH)->MoveWindow(CEedit_xoff, rectButton.top, rectButton.left - 2 * CEedit_xoff, rectButton.Height());

	nCurrentBackgroundHeight = rectClient.Height() - rectButton.Height() - 3 * CEedit_yoff;

	// now we place the tree control
	GetDlgItem(IDC_TREECUSTOMEDIT)->MoveWindow(CEedit_xoff, CEedit_yoff, TreeViewSize, nCurrentBackgroundHeight, FALSE);
	// and the background
	wndStaticBackground.MoveWindow(2 * CEedit_xoff + TreeViewSize, CEedit_yoff, rectClient.right - 3 * CEedit_xoff - TreeViewSize, nCurrentBackgroundHeight, false);

	UpdateScrollBarInfo();

	RedrawWindow();
}

void CTreeViewCustomEditDialog::UpdateScrollBarInfo()
{
	SCROLLINFO si;
	si.cbSize = sizeof(SCROLLINFO);
	si.fMask = SIF_ALL; // SIF_ALL = SIF_PAGE | SIF_RANGE | SIF_POS;
	si.nMin = 0;
	si.nMax = nVerticalSpaceForTheControls - nCurrentBackgroundHeight;
	si.nPage = si.nMax / 3 + 1;
	si.nPos = 0;
  SetScrollInfo(SB_VERT, &si, TRUE); 
}

void CTreeViewCustomEditDialog::OnContextMenu(CWnd* pWnd, CPoint point)
{
	CMenu contextmenu;
	contextmenu.CreatePopupMenu();
	contextmenu.AppendMenu(MF_STRING, ID_CM_CLIPBOARD, "Copy To Clipboard");
	contextmenu.TrackPopupMenu(TPM_LEFTALIGN| TPM_RIGHTBUTTON,point.x,point.y,this,NULL);
}

void CTreeViewCustomEditDialog::OnCopyToClipboard()
{
// open Clipboard
	if (!OpenClipboard())
	{
    AfxMessageBox( "Cannot open the Clipboard" );
    return;
  }
	if (!EmptyClipboard())
	{
    AfxMessageBox( "Cannot empty the Clipboard" );
    return;
	}

// Get Text
	//THINK ABOUT: include the current tree 
	ElementDataContainer * edc = GetCurrentEDC();
	pWCDI->CallCompFunction(112,0,0,edc); // this function changes the pionter edc !!!
	int entrynr = edc->Find("ClipBoardText");
	mystr text;
	if (entrynr > 0)
	{
		text = edc->TreeGetString("ClipBoardText");
		edc->Delete(entrynr);
	}

// write into global buffer - with lock
	int bytesize = (text.Length()+2)*sizeof(char);
	HGLOBAL data = GlobalAlloc(GMEM_MOVEABLE|GMEM_DDESHARE|GMEM_ZEROINIT,bytesize);
	char* buf = (char*) GlobalLock(data);  
	memcpy(buf,text.c_str(),bytesize);
	GlobalUnlock(data);

// Set Handle to the c_str
	HANDLE cbdata = SetClipboardData(CF_TEXT,data);
	
// close th Clipboard
	CloseClipboard();
}