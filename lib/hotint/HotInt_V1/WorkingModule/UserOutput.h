//#**************************************************************
//# filename:             UserOutput.h
//#
//# author:               Yury Vetyukov
//#
//# generated:						2003
//# description:          
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
 
#ifndef __USER_OUTPUT_H__
#define __USER_OUTPUT_H__


#include "ioincludes.h"

#include "WinCompDriverInterface.h"


#include <assert.h>
#include <memory.h>
#include <math.h>

#include "tarray.h"    //+i2
#include "mystring.h"
#include "femath.h"

#include "useroutputinterface.h"

//#define ENDL "\n"
//const char* ENDL="\n";



class UserOutput : public UserOutputInterface
{
	char buf[100];
	char e_buf[100];
	char f_buf[100];

	int localmessagelevel;              // level of message, set by UO(UO_MSGLVL), valid until next call UO()
	int* globalmessagelevel;						// global message level
	//PG:
	int* globalfilemessagelevel;				// global log file message level
	ofstream* logfile;                  // hotint logfile
	double* critical_logfile_size;      // if the file size is passing a multiple of this factor, then a warning is displayed
	int counter_critical_logfile_size;  // counts, how often a multiple of *critical_logfile_size has been exceeded

	int output_precision;								// output precision for the next output
	int* oprec_double;								  // output precision of a double
	int* oprec_vector;								  // output precision of a double number in vector (s)
	int* oprec_matrix;								  // output precision of a double number in matrix (s)
 
	int* maxerrors;											// maximum number of displayed errors, linked to Solver_Options 
  int* maxwarnings;										// maximum number of displayed warningss, linked to Solver_Options 
// counters
	int errorcount;
	int warningcount;
// defaults, to be linked to options or set on access	
	int default_localmessagelevel;
	int default_globalmessagelevel;

	int default_oprec;
	int default_oprec_double;								 
	int default_oprec_vector;								 
	int default_oprec_matrix;								 

	int default_maxerror;
	int default_maxwarnings;

// used to save the current localmessagelevel, in order to reset to this level later
	int localmessagelevel_save;


public:

	WCDInterface::UserInterface * pUI;
	
	UserOutput() : pUI(NULL) 
	{
// reset counters
		errorcount = 0;
		warningcount = 0;
// defaults - valid until linked to options
		default_localmessagelevel = UO_LVL_all;
		default_globalmessagelevel = UO_LVL_all;
		default_oprec = 8;
		default_oprec_double = 8;								 
		default_oprec_vector = 8;								 
		default_oprec_matrix = 8;								 
		default_maxerror = 100;
		default_maxwarnings = 100;

		localmessagelevel = default_localmessagelevel;
		globalmessagelevel = &default_globalmessagelevel;
		globalfilemessagelevel = &default_globalmessagelevel;

		output_precision = default_oprec;
		oprec_double = &default_oprec_double;								  
		oprec_vector = &default_oprec_vector;								  
		oprec_matrix = &default_oprec_matrix;								  
		maxerrors = &default_maxerror;
		maxwarnings = &default_maxwarnings;

		logfile = NULL;
		critical_logfile_size = NULL;
		counter_critical_logfile_size = 0;
	}

// sets message level at access & until next access, increases counters
  void SetLocalMessageLevel(UO_MSGLVL message_level)
	{ 
		localmessagelevel=message_level;
		if (localmessagelevel==UO_LVL_err) CountError();
		if (localmessagelevel==UO_LVL_warn) CountWarning();
	}
// sets message level at access & until next access, increases counters
	void SetLocalMessageLevel(int message_level) { SetLocalMessageLevel( (UO_MSGLVL) message_level); }
	int GetLocalMessageLevel() { return localmessagelevel; }

// save the current localmessagelevel
	void SaveLocalMessageLevel()
	{
		localmessagelevel_save = localmessagelevel;
	}
//reset to the saved localmessagelevel
	void ResetLocalMessageLevel()
	{ 
		SetLocalMessageLevel(localmessagelevel_save);
	}

	int GetGlobalMessageLevel() { return *globalmessagelevel; }
	int GetGlobalFileMessageLevel() { return *globalfilemessagelevel; }

	void SetOutputPrec(int output_prec) { output_precision = output_prec; }

// links UserOutput to Solver_Options, to be called after UserOutput object is created
	void HookToSolverSettings(int* globalmessagelevel, int* globalfilemessagelevel, ofstream* logfile, double* critical_logfile_size,
		int* maxerrors, int* maxwarnings,	int* oprec_double, int* oprec_vector, int* oprec_matrix) 
	{ 
		UserOutput::globalmessagelevel = globalmessagelevel;
		UserOutput::globalfilemessagelevel = globalfilemessagelevel;
		UserOutput::logfile = logfile;
		UserOutput::critical_logfile_size = critical_logfile_size;
		UserOutput::maxerrors = maxerrors;
		UserOutput::maxwarnings = maxwarnings;
		UserOutput::oprec_double = oprec_double;
		UserOutput::oprec_vector = oprec_vector;
		UserOutput::oprec_matrix = oprec_matrix; 
	}

	int& MaxErrors() {return *maxerrors;}
	int& MaxWarnings() {return *maxwarnings;}

	int PrintMsg() // decides whether to print the message or not
	{ 
		if (pUI==NULL) return 0;																							// do not print in invalid UserInterface
		if (*globalmessagelevel < localmessagelevel) return 0;								// do not print when local message level is larger then global
		if ((localmessagelevel == UO_LVL_err) && ( errorcount > (*maxerrors) )) return 0;       // print no more errors
		if ((localmessagelevel == UO_LVL_warn) && ( warningcount > (*maxwarnings) )) return 0;  // print no more warnings
		return 1;
	}
//PG:
	int PrintMsgToLogFile() // decides whether to print the message to log file or not
	{ 
		if (logfile==NULL) return 0;																							// do not print in invalid UserInterface
		if (*globalfilemessagelevel < localmessagelevel && *globalmessagelevel < localmessagelevel) return 0;								// do not print when local message level is larger then global
		if ((localmessagelevel == UO_LVL_err) && ( errorcount > (*maxerrors) )) return 0;       // print no more errors
		if ((localmessagelevel == UO_LVL_warn) && ( warningcount > (*maxwarnings) )) return 0;  // print no more warnings
		return 1;
	}

	void CountError() // counts errors, displays message when maxerrors is reached
	{
		errorcount++;
		if (errorcount == (*maxerrors)+1) 
		{ 
			pUI->AddText("Maximum number of displayed errors reached, no more Errormessages will be displayed!\n");
		}
	}

	void CountWarning() // counts warnings, displays message when maxwarnings is reached
	{
		warningcount++;
		if (warningcount == (*maxwarnings)+1) 
		{ 
			pUI->AddText("Maximum number of displayed warnings reached, no more Warningmessages will be displayed!\n");
		}
	}

	virtual void InstantMessageText(const char* pStr)
	{
		pUI->InstantMessageText(pStr);
	}
/*
	void AddLevelTag() // On Access start line with a level tag for some
	{
		if (!PrintMsg()) return;
		switch (localmessagelevel)
		{
		case UO_LVL_0: break;  // no Output
		case UO_LVL_err: AddText("Error: "); break;  // necessary output (Errors, start/end simulation)
		case UO_LVL_warn: AddText("Warning: "); break;// almost necessary output (Warnings)
		case UO_LVL_ext: break;  // extended output (useful information)
		case UO_LVL_all: break;  // complete information
		case UO_LVL_dbg1: break; // debug level 1
		case UO_LVL_dbg2: break;  // debug level 2
		default: break;
		}
	}
*/
	//PG:
	//void AddText(const char * pStr)
	//{		
	//	if (!PrintMsg()) return;
	//	pUI->AddText(pStr);
	//}
	//void AddText(double x)
	//{
	//	if (!PrintMsg()) return;
	//	sprintf_s(buf,"%.*g",output_precision,x);
	//	pUI->AddText(buf);
	//}
	//void AddText(int x)
	//{
	//	if (!PrintMsg()) return;
	//	sprintf_s(buf,"%d",x);
	//	pUI->AddText(buf);
	//}
	void AddText(const char * pStr)
	{		
		if (PrintMsg())	pUI->AddText(pStr);
		if (PrintMsgToLogFile()) 
		{
			(*logfile) << pStr;
			logfile->flush();
			
			// print warning, if log file is exceeding critical size
			if (*globalfilemessagelevel >= UO_LVL_warn || *globalmessagelevel >= UO_LVL_warn)
			{
				if (logfile->rdbuf()->pubseekoff(0,ios::end,ios::out) > (*critical_logfile_size)*(counter_critical_logfile_size + 1)*1048576.)   //critical_logfile_size is given in megabytes, while pubseekoff returns size in bytes
				{	
					counter_critical_logfile_size++;
					int memo = localmessagelevel;
					localmessagelevel = UO_LVL_warn;
					AddText(mystr("WARNING: size of log file exceeds critical size of ") + mystr((*critical_logfile_size)*counter_critical_logfile_size) + mystr(" MB.\n"));
					localmessagelevel = memo;
				}
			}
		}
		return;
	}
	void AddText(double x)
	{
		sprintf_s(buf,"%.*g",output_precision,x);
		AddText(buf);
	}
	void AddText(int x)
	{
		sprintf_s(buf,"%d",x);
		AddText(buf);
	}

	virtual UserOutputInterface & operator <<(const char * pStr) 
	{
		AddText(pStr);
		return *this;
	}
	virtual UserOutputInterface & operator <<(int x) 
	{
		sprintf_s(buf,"%d",x);
		AddText(buf);
		return *this;
	}
	virtual UserOutputInterface & operator <<(double x)
	{
		if (output_precision==-1) output_precision = *oprec_double;
		sprintf_s(buf,"%.*g",output_precision,x);
		AddText(buf);
		return *this;
	}
	virtual UserOutputInterface & operator <<(const Vector3D & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		str += mystr("[");
		str += mystr(v.X(),output_precision);
		str += mystr(", ");
		str += mystr(v.Y(),output_precision);
		str += mystr(", ");
		str += mystr(v.Z(),output_precision);
		str += mystr("]");
		AddText(str.c_str());
		return *this;
	}
	virtual UserOutputInterface & operator <<(const Box3D & b)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		str += mystr("[");
		str += mystr(b.PMin().X(),output_precision);
		str += mystr(", ");
		str += mystr(b.PMin().Y(),output_precision);
		str += mystr(", ");
		str += mystr(b.PMin().Z(),output_precision);
		str += mystr("] - [");
		str += mystr(b.PMax().X(),output_precision);
		str += mystr(", ");
		str += mystr(b.PMax().Y(),output_precision);
		str += mystr(", ");
		str += mystr(b.PMax().Z(),output_precision);
		str += mystr("]");
		AddText(str.c_str());
		return *this;
	}
	virtual UserOutputInterface & operator <<(const Vector2D & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		str += mystr("[");
		str += mystr(v.X(),output_precision);
		str += mystr(", ");
		str += mystr(v.Y(),output_precision);
		str += mystr("]");
		AddText(str.c_str());
		return *this;
	}
	virtual UserOutputInterface & operator <<(const Vector & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		if (v.GetLen()<=70)
		{
			str+=mystr("[");
			for (int i=1; i<v.GetLen(); i++)
			{
				str+=mystr(v(i),output_precision);
				str+=mystr(", ");
			}
			if (v.GetLen()!=0)
			{
				str+=mystr(v(v.GetLen()),output_precision);
			}
			str+=mystr("]");
		} 
		else
		{

			str+=mystr("[");
			for (int i=1; i<=v.GetLen(); i++)
			{
				str+=mystr(v(i),output_precision);
				if (i!=v.GetLen()) {str+=mystr(", ");}
				if (i%16==0 && i!=v.GetLen()) {str+=mystr("\n  col("); str+=mystr(i); str+=mystr( "): ");}
			}
			str+=mystr("]");

		}
		AddText(str.c_str());
		return *this;
	}

	virtual UserOutputInterface & operator <<(const IVector & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		if (v.Length()<=70)
		{
			str+=mystr("[");
			for (int i=1; i<v.Length(); i++)
			{
				str+=mystr(v(i));
				str+=mystr(", ");
			}
			if (v.Length()!=0)
			{
				str+=mystr(v(v.Length()));
			}
			str+=mystr("]");
		}
		else
		{
			str+=mystr("[");
			for (int i=1; i<=v.Length(); i++)
			{
				str+=mystr(v(i));
				if (i!=v.Length()) {str+=mystr(", ");}
				if (i%16==0 && i!=v.Length()) {str+=mystr("\n  col("); str+=mystr(i); str+=mystr( "): ");}
			}
			str+=mystr("]");

		}
		AddText(str.c_str());
		return *this;
	}

	virtual UserOutputInterface & operator <<(const TArray<double> & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		if (v.Length()<=70)
		{
			str+=mystr("[");
			for (int i=1; i<v.Length(); i++)
			{
				str+=mystr(v(i),output_precision);
				str+=mystr(", ");
			}
			if (v.Length()!=0)
			{
				str+=mystr(v(v.Length()),output_precision);
			}
			str+=mystr("]");
		} 
		else
		{

			str+=mystr("[");
			for (int i=1; i<=v.Length(); i++)
			{
				str+=mystr(v(i),output_precision);
				if (i!=v.Length()) {str+=mystr(", ");}
				if (i%16==0 && i!=v.Length()) {str+=mystr("\n  col("); str+=mystr(i); str+=mystr( "): ");}
			}
			str+=mystr("]");

		}
		AddText(str.c_str());
		return *this;
	}


	virtual UserOutputInterface & operator <<(const SparseVector & v)
	{
		if (output_precision==-1) output_precision = *oprec_vector;
		mystr str;
		if (v.GetLen()<=70)
		{
			str+=mystr("[");
			for (int i=1; i<v.GetLen(); i++)
			{
				str+=mystr(v(i),output_precision);
				str+=mystr(", ");
			}
			if (v.GetLen()!=0)
			{
				str+=mystr(v(v.GetLen()),output_precision);
			}
			str+=mystr("]");
		} else
		{

			str+=mystr("[");
			for (int i=1; i<=v.GetLen(); i++)
			{
				str+=mystr(v(i),output_precision);
				if (i!=v.GetLen()) {str+=mystr(", ");}
				if (i%16==0 && i!=v.GetLen()) {str+=mystr("\n  col("); str+=mystr(i); str+=mystr( "): ");}
			}
			str+=mystr("]");

		}
		AddText(str.c_str());
		return *this;
	}

	virtual UserOutputInterface & operator <<(const Matrix3D & m)
	{
		if (output_precision==-1) output_precision = *oprec_matrix;
		mystr str;
		for (int i=0; i < m.Getrows(); i++)
		{
			str+=mystr("[");
			for (int j=0; j < m.Getcols()-1; j++)
			{
				str+=mystr(m.Get0(i,j),output_precision);
				str+=", ";
			}
			str+=mystr(m.Get0(i,m.Getcols()-1),output_precision);
			str+="]\n";
		}
		AddText(str.c_str());
		return *this;
	}

	virtual UserOutputInterface & operator <<(const Matrix & m)
	{
		if (output_precision==-1) output_precision = *oprec_matrix;
		mystr str;
		double max = m.MaxNorm();
		int i;

		if (max==0) {max=1;}
		max=(int)log10(max);
		max=pow(10,max);

		str+=mystr(max,output_precision); 
		str+=mystr(" *\n");

		for (i=1; i<=m.Getrows(); i++)
		{
			str+=mystr("[");
			for (int j=1; j<=m.Getcols(); j++)
			{
				char str1[32];
				sprintf_s(str1,"% 1.*f",output_precision,m(i,j)/max);
				str+=mystr(str1);

				if (j!=m.Getcols()) {str+=mystr(",");}
			}
			str+=mystr("]\n");
		}
		AddText(str.c_str());
		return *this;
	}

	virtual UserOutputInterface & operator <<(const SparseMatrix & m)
	{
		if (output_precision==-1) output_precision = *oprec_matrix;
		mystr str;
		double max = m.MaxNorm();
		int i;

		if (max==0) {max=1;}
		max=(int)log10(max);
		max=pow(10,max);

		str+=mystr(max); str+=mystr(" *\n");

		for (i=1; i<=m.Getrows(); i++)
		{
			str+=mystr("[");
			for (int j=1; j<=m.Getcols(); j++)
			{
				char str1[32];
				sprintf_s(str1,"% 1.*f",output_precision,m(i,j)/max);
				str+=mystr(str1);
	
				if (j!=m.Getcols()) {str+=mystr(",");}
			}
			str+=mystr("]\n");
		}
		AddText(str.c_str());
		return *this;
	}

	virtual int CallWCDriverFunction(int action, int option = 0, int value = 0, ElementDataContainer* edc = NULL)
	{
		return pUI->CallWCDriverFunction(action, option, value, edc);
	}
};


#endif	// __USER_OUTPUT_H__