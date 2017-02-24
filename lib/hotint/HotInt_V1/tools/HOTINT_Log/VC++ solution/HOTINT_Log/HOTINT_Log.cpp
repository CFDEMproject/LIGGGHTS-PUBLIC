//#**************************************************************
//#
//# filename:             HOTINT_Log.cpp
//#
//# author:               Erwin Karer
//#
//# generated:						
//# description:          
//#                       
//# remarks:						  
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

// HOTINT_Log.cpp: Hauptprojektdatei.

#include "stdafx.h"
#include "Form1.h"
//#include "..\..\..\MBSKernelLib\tarray.h"
//#include "..\..\..\MBSKernelLib\mystring.cpp"
//#include "..\..\..\MBSKernelLib\myfile.cpp"


using namespace HOTINT_Log;


[STAThreadAttribute]
int main(array<System::String ^> ^args)
{
	// Aktivieren visueller Effekte von Windows XP, bevor Steuerelemente erstellt werden
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false); 

	// Hauptfenster erstellen und ausführen
	//Form1* form = gcnew Form1();
	Application::Run(gcnew Form1());


	return 0;
}


