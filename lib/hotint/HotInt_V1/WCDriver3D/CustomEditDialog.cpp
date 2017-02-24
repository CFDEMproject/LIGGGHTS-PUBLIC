//#**************************************************************
//# filename:             CustomEditDialog.cpp
//#
//# author:               Gerstmayr, Vetyukov
//#
//# generated:						
//# description:          ComputeEigenmodes-Dialog
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
#include "TreeViewCustomEditDialog.h"
#include "mystring.h"

// CustomEditDialog-Dialogfeld

IMPLEMENT_DYNAMIC(CustomEditDialog, CDialog)
CustomEditDialog::CustomEditDialog(CWnd* pParent /*=NULL*/)
: CDialog(CustomEditDialog::IDD, pParent), 	cedit(), cstatic(), cbutton(), ceditID(), cstaticID(), cbuttonID(), 
groupbutton(), tooltiptext(), tooltiptextID()

{
	dialogname = "Edit properties";
	add_deletebutton = 0;
	deleteflag = 0;
}

CustomEditDialog::~CustomEditDialog()
{
	Delete();
}

void CustomEditDialog::Delete()
{
	for (int i=1; i <= cedit.Length(); i++)
	{
		if (cedit(i) != 0) delete cedit(i);
	}
	for (int i=1; i <= cstatic.Length(); i++)
	{
		if (cstatic(i) != 0) delete cstatic(i);
	}
	for (int i=1; i <= cbutton.Length(); i++)
	{
		if (cbutton(i) != 0) delete cbutton(i);
	}
	for (int i=1; i <= tooltiptext.Length(); i++)
	{
		if (tooltiptext(i) != 0) delete tooltiptext(i);
	}
	cedit.SetLen(0);
	cstatic.SetLen(0);
	cbutton.SetLen(0);
	tooltiptext.SetLen(0);
	tooltiptextID.SetLen(0);

	ceditID.SetLen(0);
	cstaticID.SetLen(0);
	cbuttonID.SetLen(0);
	groupbutton.SetLen(0);
}


void CustomEditDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CustomEditDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
	ON_NOTIFY_EX( TTN_NEEDTEXT, 0, OnToolTipNotify )
END_MESSAGE_MAP()


// CustomEditDialog-Meldungshandler
void CustomEditDialog::Create(CWnd * pParent)
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

}

void CustomEditDialog::OnBnClickedOk()
{
	WriteData();
	OnOK();
}

void CustomEditDialog::OnBnClickedCancel()
{
	OnCancel();
}


void CustomEditDialog::WriteData()
{
	//fill in element data:

	int ceID = CEid_offset2;
	//char str[256];
	CString cstr;

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
				ed.SetBool(IsDlgButtonChecked(CEid_offset1+i));
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
					GetDlgItemText(ceID, cstr);

					if (!ed.IsLocked())
					{
						if (ed.IsText()) 
						{
							ed.SetText((LPCTSTR)cstr);
						}
						else
						{
							if (ed.IsInt()) 
							{
								ed.SetInt(atoi((LPCTSTR)cstr));
							}
							else if (ed.IsDouble()) 
							{
								ed.SetDouble(atof((LPCTSTR)cstr));
							}
							else //Vector
							{
								if (!editmode)
								{
									ed.SetVectorVal(j, atof((LPCTSTR)cstr));
								}
								else if (ed.IsVector())
								{
									//+++++++++++++++++++++++++++++++++++++++++++++
									mystr mstr = (LPCTSTR)cstr;
									doublelist.SetLen(0);

									//AfxMessageBox(mstr.c_str());
									int end = 0;
									int endpos;

									while (!end)
									{
										endpos = mstr.Find(',');
										if (endpos != -1 && endpos != 0 && endpos <= mstr.Length()-1)
										{
											mystr ss = mstr.SubString(0, endpos-1);
											doublelist.Add(atof(ss.c_str()));
											//AfxMessageBox(ss.c_str());
											mstr = mstr.SubString(endpos+1, mstr.Length()-1);
										}
										else 
										{
											mstr.EraseSpaces();
											if (mstr.Length() != 0) doublelist.Add(atof(mstr.c_str())); //last value
											end = 1;
										}
									}

									ed.SetVectorLen(doublelist.Length());
									for (int jj=1; jj <= doublelist.Length(); jj++)
									{
										ed.SetVectorVal(jj, doublelist(jj));
									}
								}
								else if (ed.IsMatrix())
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
													doublelist.Add(atof(str));
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
								}
							}
						}
					}
					ceID++;
				}
			}
		}
	}
}


BOOL CustomEditDialog::OnInitDialog()
{
	CDialog::OnInitDialog();

	CreateDialogItems();

	return TRUE;  // return TRUE unless you set the focus to a control
}

void CustomEditDialog::CreateDialogItems()
{
	//set dialog name
	SetWindowText(dialogname);

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
	grouplength = groupnum;

	int xoff = CEedit_xoff;
	int yoff = CEedit_yoff;

	int hwindow = ned*(CEedit_height + CEedit_vspace) + 2*yoff + CEbutton_offy;
	int max_window_width = CEstatictext_width + 3*(CEedit_hspace + CEedit_width);
	int wwindow = max_window_width + 2*xoff;

	//position and size window:
	CRect rw, cwr;

	GetWindowRect(&rw);
	GetClientRect(&cwr);
	hwindow += rw.Height() - cwr.Height();
	wwindow += rw.Width() - cwr.Width();

	if (hwindow < rw.Height()) hwindow = rw.Height();

	CPoint p((rw.left+rw.right)/2,(rw.top+rw.bottom)/2);

	rw.left = p.x - (long)(0.5*wwindow);
	rw.top = p.y - (long)(0.5*hwindow);
	rw.right = p.x + (long)(0.5*wwindow);
	rw.bottom = p.y + (long)(0.5*hwindow);

	int offx=0, offy=0;
	if (rw.top < 1) offy = 1-rw.top;
	rw.top += offy;
	rw.bottom += offy;
	if (rw.left < 1) offx = 1-rw.left;
	rw.left += offx;
	rw.right += offx;

	MoveWindow(rw,FALSE);
	ShowWindow(SW_SHOW);
	GetClientRect(cwr); //update new client window size

	//fill in element data:
	CFont* pfont = GetDlgItem(IDOK)->GetFont();

	//reposition IDOK and IDCANCEL

	CRect r;
	GetDlgItem(IDCANCEL)->GetWindowRect(r);

	int w = r.Width();
	int h = r.Height();

	r.left = CEedit_xoff;
	r.right = CEedit_xoff+w;
	r.bottom = cwr.bottom-CEedit_yoff;
	r.top = cwr.bottom-CEedit_yoff-h;

	GetDlgItem(IDCANCEL)->MoveWindow(r,TRUE);
	GetDlgItem(IDCANCEL)->ShowWindow(SW_SHOW);


	GetDlgItem(IDOK)->GetWindowRect(r);
	w = r.Width();
	h = r.Height();

	r.left = cwr.right-CEedit_xoff-w;
	r.right = cwr.right-CEedit_xoff;
	r.bottom = cwr.bottom-CEedit_yoff;
	r.top = cwr.bottom-CEedit_yoff-h;


	//char str3[256];
	//sprintf(str3, "cwr: x1=%d, y1=%d, x1=%d, y1=%d; wr: x1=%d, y1=%d, x1=%d, y1=%d", cwr.left, cwr.top, cwr.right, cwr.bottom, rw.left, rw.top, rw.right, rw.bottom);
	//SetWindowText(str3);

	GetDlgItem(IDOK)->MoveWindow(r,TRUE);
	GetDlgItem(IDOK)->ShowWindow(SW_SHOW);

	//add delete button?
	if (add_deletebutton)
	{
		MyCButton* cb = new MyCButton();
		cbutton.Add(cb);

		int x = (cwr.left+cwr.right)/2;

		cb->Create("Delete", BS_CENTER|BS_VCENTER|WS_CHILD|WS_VISIBLE, 
			CRect(x-w/2,cwr.bottom-CEedit_yoff-h,x+w/2,cwr.bottom-CEedit_yoff), this, CEid_offset1-1);
		cb->SetCustomEditDialog(this, 0); //tell button where to call if button pressed, 0=delete button!
		cb->SetFont(pfont);
	}

	EnableToolTips();

	int ceID = CEid_offset2;
	char str[256];
	char str2[256];

	for (int i=1; i <= groupnum.Length(); i++) {groupnum(i) = 0;}
	int yoff_reduce = 0;

	for (int i=1; i <= edc->Length(); i++)
	{
		ElementData& ed = edc->Get(i);
		int px1 = xoff;
		int px2 = xoff + CEstatictext_width + CEedit_hspace;
		int py = yoff + (i-1-yoff_reduce)*(CEedit_height + CEedit_vspace);

		if (ed.IsBool() || ed.IsCompAction() || ed.IsWCDriverAction())
		{
			//add event handle for button pressed!!!!

			//cstatic.Add(0); //do not delete this static entry

			if (ed.HasToolTip()) //add tooltip text to list including the itemID
			{
				AddToolTipText(ed.GetToolTipText(), CEid_offset1+i);
			}


			if (ed.IsBool())
			{
				CButton* cb = new CButton();
				cbutton.Add(cb);
				//checkbox:
				if (!ed.IsGroup())
				{
					cb->Create(ed.GetDataName(), BS_AUTOCHECKBOX|BS_LEFT|WS_CHILD|WS_VISIBLE, 
						CRect(px1,py,px1+CEstatictext_width,py+CEedit_height), this, CEid_offset1+i);
					cb->SetCheck(ed.GetBool());
				} 
				else
				{
					int ng = ed.GetGroup();

					int groupwidth = CEbutton_width+CEedit_hspace;
					if (groupwidth*grouplength(ng) > max_window_width) groupwidth = max_window_width / grouplength(ng);

					CButton* group = 0;
					if (groupbutton.Length() >= ng) group = cbutton(groupbutton(ng));
					else
					{ //create new group for radio buttons:
						group = new CButton();
						groupbutton(ng) = cbutton.Add(group);
						group->Create("", BS_GROUPBOX|WS_CHILD|WS_VISIBLE, 
							CRect(px1,py-7,px1+CEstatictext_width + 3*(CEedit_hspace + CEedit_width),py+2+CEedit_height), this, CEid_offset_group+ng);
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
					CRect(px1+x,py+y,px1+x+CEbutton_width,py+y+CEedit_height), this, CEid_offset1+i);
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
				CRect(px1,py,px1+CEstatictext_width,py+CEedit_height), this, CEid_offset1+i);
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

					if (ed.HasToolTip())//add tooltip text to list including the itemID
					{
						AddToolTipText(ed.GetToolTipText(), ceID);
					}

					if (ed.IsInt())
					{
						double mi, ma;
						bool oneL;
						ed.GetMinMaxVal(mi, ma, oneL);
						ce->SetMinMax(mi, ma,oneL);
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
					ce->Create(dwstyle, CRect(px2,py,px2+ew*xfact+CEedit_hspace*(xfact-1),py+height), this, ceID);
					ce->SetFont(pfont);

					if (ed.IsText()) 
					{
						SetDlgItemText(ceID, ed.GetText());
					}
					else
					{
						if (ed.IsInt()) 
						{
							sprintf(str, "%d", ed.GetInt());
							SetDlgItemText(ceID, str);
						}
						else if (!editmode)
						{
							sprintf(str, "%.16g", ed.GetVectorVal(j));
							sprintf(str2, "%.14g", ed.GetVectorVal(j));

							if (strlen(str) > strlen(str2)+4)
								sprintf(str, "%.14g", ed.GetVectorVal(j));

							SetDlgItemText(ceID, str);
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
							SetDlgItemText(ceID, mstr.c_str());
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
							SetDlgItemText(ceID, mstr);
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

//Notification handler
BOOL CustomEditDialog::OnToolTipNotify(UINT id, NMHDR *pNMHDR, LRESULT *pResult)
{
	// need to handle both ANSI and UNICODE versions of the message
	TOOLTIPTEXTA* pTTTA = (TOOLTIPTEXTA*)pNMHDR;
	TOOLTIPTEXTW* pTTTW = (TOOLTIPTEXTW*)pNMHDR;
	CString strTipText;
	UINT nID = pNMHDR->idFrom;
	if (pNMHDR->code == TTN_NEEDTEXTA && (pTTTA->uFlags & TTF_IDISHWND) ||
		pNMHDR->code == TTN_NEEDTEXTW && (pTTTW->uFlags & TTF_IDISHWND))
	{
		// idFrom is actually the HWND of the tool
		nID = ::GetDlgCtrlID((HWND)nID);
	}

	if (nID != 0) // will be zero on a separator
	{
		int n = FindToolTipID(nID);
		if (n)
		{
			strTipText.Format("%s", tooltiptext(n));
		}
	}

	if (pNMHDR->code == TTN_NEEDTEXTA)
	  //lstrcpyn(pTTTA->szText, strTipText, sizeof(pTTTA->szText));         //$ RL 2011-11-17: old: only 80 characters shown of ToolTipText
		lstrcpyn(pTTTA->lpszText, strTipText, strTipText.GetAllocLength()+50);//$ RL 2011-11-17: Set full ToolTipText
	else
		_mbstowcsz(pTTTW->szText, strTipText, sizeof(pTTTW->szText)); // sizeof(pTTTW->szText) ...  max. length of shown tool tip text
	*pResult = 0;

	return TRUE;    // message was handled
}

void CustomEditDialog::AddToolTipText(const char* str, int itemID)
{
	int length = strlen(str);
	char* str2 = new char[length + 1];
	strcpy(str2, str);

	tooltiptext.Add(str2);
	tooltiptextID.Add(itemID);
}


void CustomEditDialog::ClearDialogItems()
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



void CustomEditDialog::ButtonClicked(int edcnum)
{
	//char str[100];
	//sprintf(str, "Button clicked: %d", edcnum);
	//AfxMessageBox(str);

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
				ClearDialogItems();
				CreateDialogItems();
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


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MyEdit

IMPLEMENT_DYNAMIC(MyEdit, CEdit)
MyEdit::MyEdit()
{
	minval = 1;
	maxval = 0;
	oneLimit = 0;
}

MyEdit::~MyEdit()
{
}


BEGIN_MESSAGE_MAP(MyEdit, CEdit)
	ON_WM_KILLFOCUS()
END_MESSAGE_MAP()

void MyEdit::OnKillFocus(CWnd* pNewWnd)
{
	CEdit::OnKillFocus(pNewWnd);

	////assume integer value:
	//CString str;
	//GetWindowText(str);
	//int val = atoi((LPCTSTR)str);
	//char s[256];

	//assume double value:
	CString str;
	GetWindowText(str);
	double val = atof((LPCTSTR)str);
	char s[256];


	if(!oneLimit)
	{
		if (minval <= maxval)				// borders active for minimum and maximum
		{
			if (val < minval || val > maxval)
			{
				if (val < minval) sprintf(s, "%f", minval); 
				else sprintf(s, "%f", maxval);
				SetWindowText(s);

				sprintf(s, "The entered value %g is not in the Range %g...%g", val, minval, maxval);
				AfxMessageBox(s);
			}
		}
	}
	else		// just minimum OR maximum value is used
	{
		if (minval <= maxval)				// just use maximum value
		{
			if (val > maxval)
			{
				sprintf(s, "The entered value %g is greater than maximum value %g", val, maxval);
				AfxMessageBox(s);
			}
		}
		else												// just use minimum value
		{
			if (val < minval)
			{
				sprintf(s, "The entered value %g is smaller than minimum value %g", val, minval);
				AfxMessageBox(s);
			}
		}

	}
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MyCButton

IMPLEMENT_DYNAMIC(MyCButton, CButton)
MyCButton::MyCButton()
{
	ced = 0;
	tvced = 0;
}

MyCButton::~MyCButton()
{
}


BEGIN_MESSAGE_MAP(MyCButton, CButton)
	ON_CONTROL_REFLECT(BN_CLICKED, OnBnClicked)
END_MESSAGE_MAP()



// MyCButton-Meldungshandler
void MyCButton::OnBnClicked()
{
	if (ced) ced->ButtonClicked(edcnum);
	if (tvced) tvced->ButtonClicked(edcnum);
}

void MyCButton::SetCustomEditDialog(CustomEditDialog* cedI, int edcnumI)
{
	ced = cedI;
	edcnum = edcnumI;
}

void MyCButton::SetCustomEditDialog(CTreeViewCustomEditDialog* cedI, int edcnumI)
{
	tvced = cedI;
	edcnum = edcnumI;
}


