//#***************************************************************************************
//# filename:     hotint_version_info.h
//# authors:      Peter Gruber
//# generated:    March 2013
//# description:  Data structure for storing and comparing the versions of Hotint
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
 


#ifndef __HOTINTVERSIONINFO__
#define __HOTINTVERSIONINFO__

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ PG 2013-3-7: old versioning system
//const double HOTINT_version_info = 1.001; //X.Y ... Major/minor revision
//const int HOTINT_log_number = 423;        //Z ... log-number

//$ PG 2013-3-7: new versioning system

// The following class handles version numbers of Hotint.
class HotintVersionInfo
{
	int major_revision;
	int minor_revision;
	int log_number;

public:
	HotintVersionInfo(int major_revision, int minor_revision, int log_number)
	{
		this->major_revision = major_revision;
		this->minor_revision = minor_revision;
		this->log_number = log_number;
	}
	
	HotintVersionInfo(const HotintVersionInfo& rhs)
	{
		major_revision = rhs.major_revision;
		minor_revision = rhs.minor_revision;
		log_number = rhs.log_number;
	}
	
	// In the edc and data files (etc.) double values are stored, which refer to
	// the version of Hotint. Such a double consists of the major revision number
	// berfore and the minor revision number past the comma
	// (Version 1.4.231 ->> 1.004).
	// The following constructor is used to convert a double of such type into an
	// object of HotintVersionInfo.
	HotintVersionInfo(double x)
	{
		major_revision = (int)x;
		minor_revision = (int)((x - major_revision + 1e-15)*1e3);
		log_number = 0;
	}

	HotintVersionInfo(mystr str)	//$ DR 2013-03-20
	{
		major_revision = 0; 		minor_revision = 0; 		log_number = 0;
		mystr substr;
		int pos = 0;
		char dot = '.';
		if(str.GetUntil(pos,dot,substr))
		{
			major_revision = substr.MakeInt();
			//pos += substr.Length();
			pos++;		//$ DR 2013-12-11, offset for "dot"
			if(str.GetUntil(pos,dot,substr))
			{
				minor_revision = substr.MakeInt();
				//pos += substr.Length();
				pos++;	//$ DR 2013-12-11, offset for "dot"
				str.GetUntil(pos,dot,substr);
				log_number = substr.MakeInt();
			}
		}
	}


	int GetMajorRevision() const {return major_revision;}
	int GetMinorRevision() const {return minor_revision;}
	int GetLogNumber() const {return log_number;}
	
	// The following method returns the version of Hotint in form of a double value.
	// See also the constructor HotintVersionInfo(double x).
	double GetDoubleValue() const
	{
		return (double)major_revision + ((double)minor_revision)*1e-3;
	}

	// Return the version string (e.g., "1.3.254")
	char* GetString() const
	{
		char* buffer = new char[32];
		
		if (log_number == 0)
		{
			sprintf(buffer, "%i.%i", major_revision, minor_revision);
		}
		else
		{
			sprintf(buffer, "%i.%i.%i", major_revision, minor_revision, log_number);
		}

		return buffer;
	}
};

static bool operator> (const HotintVersionInfo& lhs, const HotintVersionInfo& rhs)
{
	// log_number is not taken into account
	if (lhs.GetMajorRevision() > rhs.GetMajorRevision()) 
		return true;
	if (lhs.GetMajorRevision() < rhs.GetMajorRevision()) 
		return false;
	
	// so we know lhs.GetMajorRevision() == rhs.GetMajorRevision()
	if (lhs.GetMinorRevision() > rhs.GetMinorRevision())
		return true;
	
	return false;
}

static bool operator== (const HotintVersionInfo& lhs, const HotintVersionInfo& rhs)
{
	// log_number is not taken into account
	if (lhs.GetMajorRevision() == rhs.GetMajorRevision() && lhs.GetMinorRevision() == rhs.GetMinorRevision())
		return true;
	return false;
}

static bool operator!= (const HotintVersionInfo& lhs, const HotintVersionInfo& rhs)
{
	// log_number is not taken into account
	return !(lhs==rhs);
}

#include "hotint_version.h"
#endif  // __HOTINTVERSIONINFO__