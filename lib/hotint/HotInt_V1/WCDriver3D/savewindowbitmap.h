//#**************************************************************
//# filename:             savewindowbitmap.h
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
 

#ifndef __SAVEWINDOWBITMAP_H_INCLUDED__
#define __SAVEWINDOWBITMAP_H_INCLUDED__

//AD&JG:
// New versions for the export of snapshots and graphs
// taking advantage of CImage class (ATL)
#include <atlimage.h>
#include <Gdiplus.h>

// Functions to write Bitmap - Object to File
int GetEncoderClsid(const WCHAR* format, CLSID* pClsid);
void CImageSave(HBITMAP& hBmp, const CString& FileName, const GUID& ImageFormat);
void GDIPlusSave(HBITMAP& hBmp, const CString& FileName, const GUID& ImageFormat);

// capture the screen region of the active window ( with any overlapping windows )
BOOL CaptureWindow(const CString& filename, const GUID& ImageFormat, CWnd* main , CWnd* aux = NULL); // identify the screen region for the screenshot
BOOL DoCapture(const POINT& coords, const SIZE& areaSize, const CString& FileName, const GUID& ImageFormat); // capture screenshot and save file
bool SaveWindowBitmap(CWnd* pWnd, const CString& FileName, const GUID& ImageFormat); // just calls the above

// keep the old code here, just in case
//PBITMAPINFO CreateBitmapInfoStruct(HBITMAP hBmp);
//void CreateBMPFile(const CString & File, PBITMAPINFO pbi, HBITMAP hBMP, HDC hDC);
//bool SaveWindowBitmap(CWnd * pWnd,const CString & FileName); // <-- OLD VERSION


#endif // __SAVEWINDOWBITMAP_H_INCLUDED__