//#**************************************************************
//#
//# filename:             mystring.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            20.05.97
//# description:          Variable size string class
//# remarks:						  This class has been developed, because at that time there was
//#												no powerful string class available for UNIX-machines.
//#												However, CString of Microsoft VS.NET can do the same!
//#												Note that a string with n characters has 0..n-1 characters 
//#												that there is a '0C' at position n
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
//#**************************************************************


#ifndef MYSTRING__H
#define MYSTRING__H

#include "ioincludes.h"

#include <string.h>
#include <math.h>

//for parsing:
int IsWordDelimiter(char d);
int IsSpace(char d);
int IsEOL(char d);
int IsNum(char d);
int IsUpperCaseLetter(char d);
int IsLowerCaseLetter(char d);
int IsLetter(char d);

class mystr
{
public:
	mystr();
	mystr(const char *);
	mystr(char);
	mystr(const mystr &);
	mystr(int);
	mystr(long);
	mystr(double);
	mystr(double,int); // explicit precision
	virtual ~mystr();
	mystr Left(int);
	mystr Right(int);
	void CopySubStringNoSpaces(char* deststr, int pos1, int pos2, int limit);
	mystr SubString(int pos1, int pos2) const;
	mystr GetWord(int& pos, int incl_delimiter=0); //returns word and moves p1 to position after delimiter or -1 for end of string
	mystr GetIdentifier(int& pos, int incl_delimiter=0); //returns identifier (word) and moves p1 to position after delimiter or -1 for end of string ==> in addition to Word, all separators excl. '.' are used (brackets, etc.)
	char PosGet(int& pos) const; //returns char and increments position 'pos' (or sets -1 if end of string)
	char PosPeek(int& pos) const; //returns char at position 'pos'
	int GetUntil(int& pos, char until, mystr& s, int incl_delimiter=0) const; //$JG2012-01-27; //pos starts at 0; reads string 's' between 'pos' and char 'until', increments pos to after until, sets -1 if end of string, returns 0 if 'until' not found
	int GetUntilEOL(int& pos, char until, mystr& s, int incl_delimiter=0) const; //pos starts at 0; reads string 's' between 'pos' and EOL or char 'until' (comment), increments pos to after until, sets -1 if end of string, returns 0 if 'until' not found
	int GoUntil(int& pos, char until) const; //reads string between 'pos' and char 'until', increments pos to after until, sets -1 if end of string, returns 0 if 'until' not found
	char ReadLeadingSpaces(int& pos) const; //$JG2012-01-27 added const//0-based pos!!! reads over leading spaces and moves pos to position after spaces, returns char at current position, pos=-1 for end of string
	char ReadLeadingSpacesAndCountLines(int& pos, int& linecnt) const; //$JG2012-01-27 added const//0-based pos!!! reads over leading spaces and moves pos to position after spaces, returns char at current position, pos=-1 for end of string; count lines in linecnt (only increase, no reset)
	mystr GetStringInBrackets(int& startpos, char ob, char cb)  const;	 //$JG2012-01-27 added const//reads string in brackets, counting opening and closing brackets, -1 if no success
	
//	void FromVector(const class Vector& vec); // writes data of vector in mystr-list: "{ , , }"
//	void FromMatrix(const class Matrix3D& mat); // writes data of matrix in mystr-list: "{ , , }"

	void EraseSpaces();
	void EraseSpacesHeadTail();
// returns position of first and last non-space letter (1-based)	
	void IdentifySpacesHeadTail(int& first, int& last); 
	void EraseChar(char c);
	void EraseChar(int pos); // removes a single character (at given position) from string
	int Replace(const mystr& searchstr, const mystr& replacestr); //find string searchstr, replace with string replacestr; this is done repeatedly; the return value counts the number of replacements
	mystr& InsertAt(int, const mystr &);
	mystr& WriteAt(int, const mystr &);
	void SmartDouble2String(double x, int nDigits = 0);
	int Length() const;
	void SetLength(int len);
	// $EK - 2013-01-17: added in order to resize (including memory reallocation) the string
	void ReSize(int len);
	int Find(const char) const; //$JG2012-01-27 added const
	int Find(int startpos, const char) const; //$JG2012-01-27 added const
	int Find(int startpos, const mystr &s) const; //$JG2012-01-27
	int Find(const char *) const; //$JG2012-01-27 added const
	int Find(const mystr &) const; //$JG2012-01-27 added const
	int FindEOL(int startpos) const; //$JG2012-01-27 added const
	mystr& operator = (const mystr &);
	friend mystr operator + (const mystr &, const mystr &);
	void operator += (const mystr &);
	char* c_str();
	const char* c_str() const;

	int IsValidNumber(); //$JG2012-01-27 added const // checks if string can contains a number
	int CountLines(); //$AD:2013-09-24 counts the lines in a multiline string

	mystr operator () (int, int) const; //$JG2012-01-27 added const
	int MakeInt() const; //$JG2012-01-27 added const
	double MakeDouble() const; //$JG2012-01-27 added const
	long MakeLong() const; //$JG2012-01-27 added const
	operator char *() const;
	char& operator [] (int);
	const char& operator [] (int) const;

// Compare function - can be case sensitive or case insensitive
// returns 1 for match
	int CStrCompare(const char* other) const //$JG2012-01-27 added const 
	{
	  return strcmp(str, other) == 0;
  }
	int Compare(const mystr& other, int casesensitive = 0) const //$JG2012-01-27 added const 
	{
	  if (casesensitive) return strcmp(str, other.str) == 0;
	  else return _stricmp(str, other.str) == 0;
  }
// Compare function - can be case sensitive or case insensitive
// returns 1 for match
//	int Compare(const char other, int casesensitive = 0) 
//	{
//	  if (casesensitive) return strcmp(str, &other) == 0;
//	  else return _stricmp(str, &other) == 0;
//  }

	friend int operator == (const mystr &, const mystr &);
	friend int operator < (const mystr &, const mystr &);
	friend int operator <= (const mystr &, const mystr &);
	friend int operator > (const mystr &, const mystr &);
	friend int operator >= (const mystr &, const mystr &);
	friend int operator != (const mystr &, const mystr &);
	friend ostream& operator << (ostream& os, const mystr& s);
	friend istream& operator >> (istream& is, mystr& s);
	static void SetToErrHandler(void (*)());

	int Trim_And_Get_Indices(int& comp1, int& comp2); // identifies and parses '[..,..]', shortens the string !! returns number of found components
	int Trim_And_Get_Indices(mystr& comp1, mystr& comp2); // identifies and parses '[..,..]', shortens the string !! returns number of found components

private:
	mystr(int, int);
	char *str;
	int length;
	static void(*ErrHandler)();
};

inline int IsWordDelimiter(char d)
{
	return (d == (char)10 || d == (char)13 || d == (char)12 || d == (char)11 || d == '\n'
		|| d == ' ' || d == ',' || d == ';' || d == '\t' 
		|| d == ':' || d == '{' || d == '}' || d == '[' || d == ']' || d == '=');
}

inline int IsWordDelimiterIdentifier(char d)
{
	return (d == (char)10 || d == (char)13 || d == (char)12 || d == (char)11 || d == '\n'
		|| d == ' ' || d == ',' || d == ';' || d == '\t' 
		|| d == ':' || d == '{' || d == '}' || d == '(' || d == ')' || d == '[' || d == ']' || d == '=');
}

inline int IsWordDelimiterBrackets(char d) // same as above but excluded Brackets (AD)
{
	return (d == (char)10 || d == (char)13 || d == (char)12 || d == (char)11 || d == '\n'
		|| d == ' ' || d == ',' || d == ';' || d == '\t' 
		|| d == ':' || d == '=');
}

inline int IsSpace(char d)
{
	return (d == (char)10 || d == (char)13 || d == (char)12 || d == (char)11 || d == '\n'
		|| d == ' ' || d == '\t');
}

inline int IsEOL(char d)
{
	return (d == (char)10 || d == (char)13 || d == (char)12 || d == '\n');
}

inline int IsNum(char d)
{
	if (d >= '0' && d <= '9') return 1;
	return 0;
}

inline int IsPoint(char d)
{
	if (d == '.') return 1;	
	return 0;
}

inline int IsSign(char d)
{
	if (d == '+' || d == '-' ) return 1;	
	return 0;
}

inline int IsScientificE(char d)
{
	if (d == 'E' || d == 'e') return 1;
	return 0;
}

inline int IsBlank(char d)
{
	if (d == ' ') return 1;
	return 0;
}

int IsCorrectNumberBetween(mystr num, int& first, int& last); 

inline int IsUpperCaseLetter(char d)
{
	if (d >= 'A' && d <= 'Z') return 1;
	return 0;
}
inline int IsLowerCaseLetter(char d) 
{
	if (d >= 'a' && d <= 'z') return 1;
	return 0;
}
inline int IsLetter(char d)
{
	return IsUpperCaseLetter(d) || IsLowerCaseLetter(d);
}

inline mystr::~mystr()
{
	delete[] str;
}

inline int mystr::Length() const
{
	return length;
}

inline void mystr::SetLength(int len)
{
	length = len;
}

inline int mystr::Find(const char c) const
{
	char *pos = strchr(str, int(c));
	return pos ? int(pos - str) : -1;
}

//startpos is 0-based!
inline int mystr::Find(int startpos, const char c) const
{
	char *pos = strchr(&str[startpos], int(c));
	return pos ? int(pos - str) : -1;
}

inline int mystr::Find(int startpos, const mystr &s) const //$JG2012-01-27 added const
{
	char *pos = strstr(&str[startpos], s.str);
	return pos ? int(pos - str) : -1;
}

inline int mystr::Find(const mystr &s) const //$JG2012-01-27 added const
{
	char *pos = strstr(str, s.str);
	return pos ? int(pos - str) : -1;
}

inline int mystr::Find(const char *s)  const //$JG2012-01-27 added const
{
	char *pos = strstr(str, s);
	return pos ? int(pos - str) : -1;
}

inline int mystr::FindEOL(int startpos)  const //$JG2012-01-27 added const
{
	for (int i=startpos; i < length; i++)
	{
		if (IsEOL(str[i])) return i;
	}
	return -1;
}

inline int mystr::MakeInt()  const //$JG2012-01-27 added const
{
	return atoi(str);
}

inline double mystr::MakeDouble()  const //$JG2012-01-27 added const
{
	return atof(str);
}

inline long mystr::MakeLong() const //$JG2012-01-27 added const
{
	return atol(str);
}

inline mystr::operator char *() const
{
	return str;
}

inline char* mystr::c_str()
{
	return str;
}

inline const char* mystr::c_str() const
{
	return str;
}

inline int operator == (const mystr &s1, const mystr& s2)
{
	return strcmp(s1.str, s2.str) == 0;
}

inline int operator < (const mystr &s1, const mystr& s2)
{
	return strcmp(s1.str, s2.str) < 0;
}

inline int operator <= (const mystr &s1, const mystr& s2)
{
	return strcmp(s1.str, s2.str) <= 0;
}

inline int operator > (const mystr &s1, const mystr& s2)
{
	return strcmp(s1.str, s2.str) > 0;
}

inline int operator >= (const mystr &s1, const mystr& s2)
{
	return strcmp(s1.str, s2.str) >= 0;
}

inline int operator != (const mystr &s1, const mystr& s2)
{
	return !(s1 == s2);
}

inline ostream& operator << (ostream& os, const mystr& s)
{
	return os << s.str;
}

inline void mystr::SetToErrHandler(void (*Handler)())
{
	ErrHandler = Handler;
};

#endif


