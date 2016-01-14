//#**************************************************************
//# filename:             savewindowbitmap.cpp
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

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include "savewindowbitmap.h"
bool no_error;

void errhandler(char * p)
{
	AfxMessageBox(p,MB_ICONEXCLAMATION);
	no_error = false;
}


// Functions to write Bitmap - Object to File
#include <windows.h>
#include <stdio.h>
using namespace Gdiplus;

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
   UINT  num = 0;          // number of image encoders
   UINT  size = 0;         // size of the image encoder array in bytes

   ImageCodecInfo* pImageCodecInfo = NULL;

   GetImageEncodersSize(&num, &size);
   if(size == 0)
      return -1;  // Failure

   pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
   if(pImageCodecInfo == NULL)
      return -1;  // Failure

   GetImageEncoders(num, size, pImageCodecInfo);

   for(UINT j = 0; j < num; ++j)
   {
      if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return j;  // Success
      }    
   }
   free(pImageCodecInfo);
   return -1;  // Failure
}

void CImageSave(HBITMAP& hBmp, const CString& FileName, const GUID& ImageFormat)
{
	CImage image;
	image.Attach(hBmp);
	image.Save(FileName,ImageFormat);
}

void GDIPlusSave(HBITMAP& hBmp, const CString& FileName, const GUID& ImageFormat)
{
	// Initialize GDI+.
	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	CLSID clsid;
	if( ImageFormat == ImageFormatJPEG) GetEncoderClsid(L"image/jpeg", &clsid);
	if( ImageFormat == ImageFormatPNG) GetEncoderClsid(L"image/png", &clsid);
	if( ImageFormat == ImageFormatBMP) GetEncoderClsid(L"image/bmp", &clsid);
	
	Gdiplus::EncoderParameters encoderParameters;
  encoderParameters.Count = 1;
  encoderParameters.Parameter[0].Guid = Gdiplus::EncoderQuality;
  encoderParameters.Parameter[0].Type = Gdiplus::EncoderParameterValueTypeLong;
  encoderParameters.Parameter[0].NumberOfValues = 1;

  ULONG quality = 100;
  encoderParameters.Parameter[0].Value = &quality;

	size_t len = FileName.GetLength()+1;
	size_t convlen = 0;
	WCHAR fnw[1024];
	mbstowcs_s(&convlen, fnw, len, FileName, _TRUNCATE);

	Gdiplus::Bitmap bitmap(hBmp, NULL);
	bitmap.Save(fnw, &clsid, &encoderParameters);
}

// capture the screen region of the active window ( with any overlapping windows )

// identify the screen region for the screenshot
BOOL CaptureWindow(const CString& filename, const GUID& ImageFormat, CWnd* main, CWnd* aux)
{
    HWND hWnd = NULL;
// screen region of the foreground window ( original code ) 
		//  hWnd = ::GetForegroundWindow();   
		//  if(!hWnd) return FALSE;
		//  CRect rect;
		//  GetWindowRect(hWnd, &rect);
// compute Screen region to capture manually
// can not use foreground window, since for recording video the foreground window may be data manager 
		hWnd = main->GetSafeHwnd();
		if(!hWnd) return false;
		CRect rect;
		GetWindowRect(hWnd, &rect);
		if(aux)
		{
			hWnd = aux->GetSafeHwnd();
			if(!hWnd) return false;
			CRect rect_aux;
	 		GetWindowRect(hWnd, &rect_aux);
			rect.UnionRect(rect,rect_aux);
 		}

    rect.NormalizeRect();
    return DoCapture(CPoint(rect.left, rect.top), CSize(rect.Width(), rect.Height()), filename, ImageFormat);
}

// capture screenshot and save file
BOOL DoCapture(const POINT& coords, const SIZE& areaSize, const CString& FileName, const GUID& ImageFormat)
{
    CDC dc;
    HDC hdc = GetDC(NULL);  // <-- We use this instead of GetWindowDC. 
                            // This is the only thing I had to change other than 
                            // getting the window coordinates in CaptureWindow()
    dc.Attach(hdc);

    // Create a memory DC into which the bitmap will be captured
    CDC memDC;
    memDC.CreateCompatibleDC(&dc);

    // If there is already a bitmap, delete it as we are going to replace it
    CBitmap bmp;
    bmp.DeleteObject();

    ICONINFO info;
    GetIconInfo((HICON)::GetCursor(), &info);   

    CURSORINFO cursor;
    cursor.cbSize = sizeof(CURSORINFO);
    GetCursorInfo(&cursor);

    bmp.CreateCompatibleBitmap(&dc, areaSize.cx, areaSize.cy);
    CBitmap * oldbm = memDC.SelectObject(&bmp);

    // Before we copy the image in, we blank the bitmap to
    // the background fill color
    memDC.FillSolidRect(&CRect(0,0,areaSize.cx, areaSize.cy), RGB(255,255,255));

    // Copy the window image from the window DC into the memory DC
    memDC.BitBlt(0, 0, areaSize.cx, areaSize.cy, &dc, coords.x, coords.y, SRCCOPY|CAPTUREBLT);
    
		memDC.SelectObject(oldbm);  

		HBITMAP hBmp = (HBITMAP) bmp.GetSafeHandle();
		CImageSave(hBmp, FileName, ImageFormat);
		GDIPlusSave(hBmp, FileName, ImageFormat);

		// AD: could be useful...
		////// Optionally copy the image to the clipboard.
    ////if(programSettings.bWantClipboard)
    ////{
        ////if(OpenClipboard(NULL))
        ////{
        ////    EmptyClipboard();
        ////    SetClipboardData(CF_BITMAP, bmp);
        ////    CloseClipboard();
        ////}
    ////}

		BOOL success = true;

		DeleteObject(bmp.Detach());
    DeleteDC(dc.Detach());
    DeleteDC(memDC.Detach());
    return success;
}



//AD: new version using ATL and CImage
//create a BMP from the DC and pass it on to a CImage object
//save the CImage in a different routine
bool SaveWindowBitmap(CWnd* pWnd, const CString& FileName, const GUID& ImageFormat)  // saves the Main Window content
{
  return CaptureWindow(FileName, ImageFormat, pWnd);

 //// no_error = true;

	////HDC hdcScreen = pWnd->GetWindowDC()->m_hDC; 
	////HDC hdcCompatible = CreateCompatibleDC(hdcScreen); 

	////CRect r;
	////pWnd->GetWindowRect(&r);
	////HBITMAP hbmScreen = CreateCompatibleBitmap(hdcScreen,r.Width(),r.Height());

	////if (hbmScreen == 0) errhandler("hbmScreen"); 
	////if (!SelectObject(hdcCompatible, hbmScreen)) errhandler("Compatible Bitmap Selection"); 
	////if (!BitBlt(hdcCompatible, 0, 0, r.Width(), r.Height(), hdcScreen, 0, 0, SRCCOPY)) errhandler("Screen to Compat Blt Failed"); 
	////
	////CImageSave(hbmScreen,FileName,ImageFormat);
	////GDIPlusSave(hbmScreen,FileName,ImageFormat);
	////
	////DeleteObject(hbmScreen);
	////DeleteDC(hdcCompatible);
	////DeleteDC(hdcScreen);
	////
	////return no_error;
}




// The following example code defines a function that uses a BITMAPINFO structure and
//  allocates memory for and initializes members within a BITMAPINFOHEADER structure.
// Note that the BITMAPINFO structure cannot be used with either a BITMAPV4HEADER or a BITMAPV5HEADER structure.
//PBITMAPINFO CreateBitmapInfoStruct(HBITMAP hBmp)
//{ 
//    BITMAP bmp; 
//    PBITMAPINFO pbmi; 
//    WORD    cClrBits; 
//	
//    // Retrieve the bitmap color format, width, and height. 
//    if (!GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bmp)) 
//        errhandler("GetObject"); 
//	
//    // Convert the color format to a count of bits. 
//    cClrBits = (WORD)(bmp.bmPlanes * bmp.bmBitsPixel); 
//    if (cClrBits == 1) 
//        cClrBits = 1; 
//    else if (cClrBits <= 4) 
//        cClrBits = 4; 
//    else if (cClrBits <= 8) 
//        cClrBits = 8; 
//    else if (cClrBits <= 16) 
//        cClrBits = 16; 
//    else if (cClrBits <= 24) 
//        cClrBits = 24; 
//    else cClrBits = 32; 
//	
//    // Allocate memory for the BITMAPINFO structure. (This structure 
//    // contains a BITMAPINFOHEADER structure and an array of RGBQUAD 
//    // data structures.) 
//	
//	if (cClrBits != 24) 
//		pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
//		sizeof(BITMAPINFOHEADER) + 
//		sizeof(RGBQUAD) * (1<< cClrBits)); 
//	
//	// There is no RGBQUAD array for the 24-bit-per-pixel format. 
//	
//	else 
//		pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
//		sizeof(BITMAPINFOHEADER)); 
//	
//    // Initialize the fields in the BITMAPINFO structure. 
//	
//    pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER); 
//    pbmi->bmiHeader.biWidth = bmp.bmWidth; 
//    pbmi->bmiHeader.biHeight = bmp.bmHeight; 
//    pbmi->bmiHeader.biPlanes = bmp.bmPlanes; 
//    pbmi->bmiHeader.biBitCount = bmp.bmBitsPixel; 
//    if (cClrBits < 24) 
//        pbmi->bmiHeader.biClrUsed = (1<<cClrBits); 
//	
//    // If the bitmap is not compressed, set the BI_RGB flag. 
//    pbmi->bmiHeader.biCompression = BI_RGB; 
//	
//    // Compute the number of bytes in the array of color 
//    // indices and store the result in biSizeImage. 
//    // For Windows NT, the width must be DWORD aligned unless 
//    // the bitmap is RLE compressed. This example shows this. 
//    // For Windows 95/98/Me, the width must be WORD aligned unless the 
//    // bitmap is RLE compressed.
//    pbmi->bmiHeader.biSizeImage = ((pbmi->bmiHeader.biWidth * cClrBits +31) & ~31) /8
//		* pbmi->bmiHeader.biHeight; 
//    // Set biClrImportant to 0, indicating that all of the 
//    // device colors are important. 
//	pbmi->bmiHeader.biClrImportant = 0; 
//	return pbmi; 
//}


// The following example code defines a function that initializes the remaining structures,
//  retrieves the array of palette indices, opens the file, copies the data, and closes the file. 
//void CreateBMPFile(const CString & File, PBITMAPINFO pbi, 
//				   HBITMAP hBMP, HDC hDC) 
//{ 
//	HANDLE hf;                 // file handle 
//    BITMAPFILEHEADER hdr;       // bitmap file-header 
//    PBITMAPINFOHEADER pbih;     // bitmap info-header 
//    LPBYTE lpBits;              // memory pointer 
//    DWORD dwTotal;              // total count of bytes 
//    DWORD cb;                   // incremental count of bytes 
//    BYTE *hp;                   // byte pointer 
//    DWORD dwTmp; 
//	
//    pbih = (PBITMAPINFOHEADER) pbi; 
//    lpBits = (LPBYTE) GlobalAlloc(GMEM_FIXED, pbih->biSizeImage);
//	
//    if (!lpBits) 
//		errhandler("GlobalAlloc"); 
//	
//    // Retrieve the color table (RGBQUAD array) and the bits 
//    // (array of palette indices) from the DIB. 
//    if (!GetDIBits(hDC, hBMP, 0, (WORD) pbih->biHeight, lpBits, pbi, 
//        DIB_RGB_COLORS)) 
//    {
//        errhandler("GetDIBits"); 
//    }
//	
//    // Create the .BMP file. 
//    hf = CreateFile(File, 
//		GENERIC_READ | GENERIC_WRITE, 
//		(DWORD) 0, 
//		NULL, 
//		CREATE_ALWAYS, 
//		FILE_ATTRIBUTE_NORMAL, 
//		(HANDLE) NULL); 
//    if (hf == INVALID_HANDLE_VALUE) 
//        errhandler("CreateFile"); 
//    hdr.bfType = 0x4d42;        // 0x42 = "B" 0x4d = "M" 
//    // Compute the size of the entire file. 
//    hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) + 
//		pbih->biSize + pbih->biClrUsed 
//		* sizeof(RGBQUAD) + pbih->biSizeImage); 
//    hdr.bfReserved1 = 0; 
//    hdr.bfReserved2 = 0; 
//	
//    // Compute the offset to the array of color indices. 
//    hdr.bfOffBits = (DWORD) sizeof(BITMAPFILEHEADER) + 
//		pbih->biSize + pbih->biClrUsed 
//		* sizeof (RGBQUAD); 
//	
//    // Copy the BITMAPFILEHEADER into the .BMP file. 
//    if (!WriteFile(hf, (LPVOID) &hdr, sizeof(BITMAPFILEHEADER), 
//        (LPDWORD) &dwTmp,  NULL)) 
//    {
//		errhandler("WriteFile"); 
//    }
//	
//    // Copy the BITMAPINFOHEADER and RGBQUAD array into the file. 
//    if (!WriteFile(hf, (LPVOID) pbih, sizeof(BITMAPINFOHEADER) 
//		+ pbih->biClrUsed * sizeof (RGBQUAD), 
//		(LPDWORD) &dwTmp, NULL)) 
//        errhandler("WriteFile"); 
//	
//    // Copy the array of color indices into the .BMP file. 
//    dwTotal = cb = pbih->biSizeImage; 
//    hp = lpBits; 
//    if (!WriteFile(hf, (LPSTR) hp, (int) cb, (LPDWORD) &dwTmp,NULL)) 
//		errhandler("WriteFile"); 
//	
//    // Close the .BMP file. 
//	if (!CloseHandle(hf)) 
//		errhandler("CloseHandle"); 
//	
//    // Free memory. 
//    GlobalFree((HGLOBAL)lpBits);
//}
//
////bool SaveWindowBitmap(CWnd * pWnd,const CString & FileName)
////{
////	no_error = true;
////	
////	// Create a normal DC and a memory DC for the entire screen. The 
////	// normal DC provides a "snapshot" of the screen contents. The 
////	// memory DC keeps a copy of this "snapshot" in the associated 
////	// bitmap. 
////	
////	HDC hdcScreen = pWnd->GetWindowDC()->m_hDC; 
////	HDC hdcCompatible = CreateCompatibleDC(hdcScreen); 
////	
////	// Create a compatible bitmap for hdcScreen.
////
////	CRect r;
////	pWnd->GetWindowRect(&r);
////	
////	HBITMAP hbmScreen = CreateCompatibleBitmap(hdcScreen,r.Width(),r.Height());
////	
////	if (hbmScreen == 0) 
////		errhandler("hbmScreen"); 
////	
////	// Select the bitmaps into the compatible DC. 
////	
////	if (!SelectObject(hdcCompatible, hbmScreen)) 
////		errhandler("Compatible Bitmap Selection"); 
////	
////	//Copy color data for the entire display into a 
////	//bitmap that is selected into a compatible DC.
////	
////	if (!BitBlt(hdcCompatible, 
////		0,0, 
////		r.Width(), r.Height(), 
////		hdcScreen, 
////		0,0, 
////		SRCCOPY)) 
////		
////        errhandler("Screen to Compat Blt Failed"); 
////
////	// now the bitmap is ready, prepare for saving
////	PBITMAPINFO pbmi = CreateBitmapInfoStruct(hbmScreen);
////
////	// and finally
////	CreateBMPFile(FileName,pbmi,hbmScreen,hdcScreen);
////
////	DeleteObject(hbmScreen);
////	DeleteDC(hdcCompatible);
////	DeleteDC(hdcScreen);
////	LocalFree(pbmi);
////	
////	return no_error;
////}