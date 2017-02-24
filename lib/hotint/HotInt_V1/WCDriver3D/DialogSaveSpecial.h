#if !defined(AFX_DIALOGSAVESPECIAL_H__F0D3CC57_A4A0_4E6A_8D90_B04CE8642EDE__INCLUDED_)
#define AFX_DIALOGSAVESPECIAL_H__F0D3CC57_A4A0_4E6A_8D90_B04CE8642EDE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialogSaveSpecial.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialogSaveSpecial dialog

class CDialogSaveSpecial : public CDialog
{
	int NumberOfDataUnits;

// Construction
public:
	CDialogSaveSpecial(int NumberOfDataUnits_ = 1,CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialogSaveSpecial)
	enum { IDD = IDD_DIALOG_SAVE_SPECIAL_PARAMS };
	int		m_FirstDataUnit;
	int		m_LastDataUnit;
	int		m_SaveEachOf;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialogSaveSpecial)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialogSaveSpecial)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALOGSAVESPECIAL_H__F0D3CC57_A4A0_4E6A_8D90_B04CE8642EDE__INCLUDED_)
