//#**************************************************************
//#
//# filename:             mystring.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            20.05.97
//# description:          Variable size string class
//# remarks:						  This class has been developed, because at that time there was
//#												no powerful string class available for UNIX-machines.
//#												However, CString of Microsoft VS.NET can do the same!
//#												Note that a string with n characters has 0..n-1 characters 
//#												that there is a '0C' at position n+1 (length=n)
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

//#include "stdafx.h" //for comment out
//#include "../WorkingModule/stdafx.h"
//#include <stdafx.h> //$!AD 2011-04-01 changed line such that tool:data_reader can compile... in case of error: add "../WorkingModule" to include path

#include "ioincludes.h"
#include <string.h>
#include <math.h>

#include "mystring.h"


extern void TIMBSWarningHandle(const char* warn, int use_instant_message_text);


void DefaultStringErrHandler()
{
	cerr << "ERROR: Invalid access to String-element!!!\n" << flush;
	//TIMBSWarningHandle(mystr("ERROR: Invalid access to String-element!!!\n"),1);
}

void (*mystr::ErrHandler)() = DefaultStringErrHandler;

mystr::mystr()
{
	length = 0;
	str = new char[1];
	str[0] = (char)0;
}

mystr::mystr(const char *s)
{
	length = strlen(s);
	str = new char[length + 1];
	strcpy(str, s);
}

mystr::mystr(char s)
{
	length = 1;
	str = new char[2];
	str[0] = s;
	str[1] = (char)0;
}

mystr::mystr(const mystr& s)
{
	length = s.length;
	str = new char[length + 1];
	strcpy(str, s.str);
}

mystr::mystr(int i)
{
	char buffer[32];
	sprintf_s(buffer, "%d", i);
	length = strlen(buffer);
	str = new char[length + 1];
	strcpy(str, buffer);
}

mystr::mystr(long l)
{
	char buffer[32];
	sprintf_s(buffer, "%ld", l);
	length = strlen(buffer);
	str = new char[length + 1];
	strcpy(str, buffer);
}

mystr::mystr(double d)
{
	char buffer[32];
	if (fabs(d) < 1E-100) {d = 0;}
	sprintf_s(buffer, "%g", d);
	length = strlen(buffer);
	str = new char[length + 1];
	strcpy(str, buffer);
}

mystr::mystr(double d, int prec) // with prec == significant digits
{
	char buffer[32];
	if (fabs(d) < 1E-100) {d = 0;}
	sprintf_s(buffer, "%.*g", prec+1,  d); //
	length = strlen(buffer);
	str = new char[length + 1];
	strcpy(str, buffer);
}

mystr::mystr(int n, int)
{
	length = n;
	str = new char[n + 1];
	str[n] = 0;
}

mystr mystr::Left(int r)
{
	if(r > length)
	{
		mystr::ErrHandler();
		mystr s;
		return s;
	}
	else
	{
		mystr tmp(r, 0);
		strncpy(tmp.str, str, r);
		return tmp;
	}
}

mystr mystr::Right(int l)
{
	if(l > length)
	{
		mystr::ErrHandler();
		mystr s;
		return s;
	}
	else
	{
		mystr tmp(l, 0);
		strncpy(tmp.str, str + length - l, l);
		return tmp;
	}
}

//inclusive the char at pos1 and pos2 (if pos1=pos2 one char is taken)
mystr mystr::SubString(int pos1, int pos2) const
{
	if (pos1<Length() && pos2 <Length() && pos1 <= pos2)
	{
		int l = pos2-pos1+1;
		mystr tmp(l, 0);

		if (l!=0)
			strncpy(tmp.str, str + pos1, l);

		return tmp;
	}
	else
	{
		mystr::ErrHandler();
		mystr s;
		return s;
	}
}

void mystr::CopySubStringNoSpaces(char* deststr, int pos1, int pos2, int limit)
{
	int j = 0;
	if (pos1<Length() && pos2 <Length() && pos1 <= pos2)
	{
		int l = pos2-pos1+1;
		for (int i=pos1; i <= pos2; i++)
		{
			if (j < limit-1 && !IsSpace(str[i]))
			{
				deststr[j++] = str[i];
			}
		}
		deststr[j] = (char)0;
	}
	else
	{
		mystr::ErrHandler();
	}
}

//0-based pos!!! reads over leading spaces and moves pos to position after spaces, returns char at current position, pos=-1 for end of string
char mystr::ReadLeadingSpaces(int& pos) const
{
	int endit = 0;
	char ch=(char)0;

	while(!endit) //read leading spaces, line feeds and carriage returns!
	{
		if (pos >= Length()) 
		{
			pos = -1;
			return (char)0;
		}
		ch = str[pos];
		if (!IsSpace(ch)) {endit = 1;}
		else {pos++;}
	}
	return ch;
}

//0-based pos!!! reads over leading spaces and moves pos to position after spaces, returns char at current position, pos=-1 for end of string
//count lines in linecnt
char mystr::ReadLeadingSpacesAndCountLines(int& pos, int& linecnt) const
{
	int endit = 0;
	char ch=(char)0;

	while(!endit) //read leading spaces, line feeds and carriage returns!
	{
		if (pos >= Length()) 
		{
			pos = -1;
			return (char)0;
		}
		ch = str[pos];
		if (!IsSpace(ch)) {endit = 1;}
		else {pos++;}
		if (IsEOL(ch)) {linecnt++;}
	}
	return ch;
}

//0-based pos!!! returns word and moves pos to position after delimiter or -1 for end of string
mystr mystr::GetWord(int& pos, int incl_delimiter)
{
	mystr s;
	char sbuf[102];
	int endit = 0;
	char ch=(char)0;

	//ch = ReadLeadingSpaces(pos); pos++
	while(!endit) //read leading spaces, line feeds and carriage returns!
	{
		if (pos >= Length()) 
		{
			pos = -1;
			return s;
		}
		ch = str[pos];
		if (!IsSpace(ch)) {endit = 1;}
		pos++;
	}

	sbuf[0] = ch;

	endit = 0;
	int bufpos = 1;
	while(!endit) //read word
	{
		if (pos >= Length()) {pos = -1; endit = 1;}
		else
		{
			ch = str[pos];

			if (IsWordDelimiter(ch)) 
			{
				endit = 1;
				if (incl_delimiter) 
				{
					sbuf[bufpos] = ch;
					bufpos++;
					pos++;
				}
			}
			else
			{
				if (bufpos >= 100) //100 !!!!!!!!!!!!
				{
					sbuf[bufpos] = (char)0;
					bufpos = 0;
					s += sbuf;
				}
				sbuf[bufpos] = ch;
				bufpos++;
				pos++;
			}
		}
	}
	sbuf[bufpos] = (char)0;
	s += sbuf;
	return s;
}

//0-based pos!!! returns word and moves pos to position after delimiter or -1 for end of string
mystr mystr::GetIdentifier(int& pos, int incl_delimiter)
{
	mystr s;
	char sbuf[102];
	int endit = 0;
	char ch=(char)0;

	//ch = ReadLeadingSpaces(pos); pos++
	while(!endit) //read leading spaces, line feeds and carriage returns!
	{
		if (pos >= Length()) 
		{
			pos = -1;
			return s;
		}
		ch = str[pos];
		if (!IsSpace(ch)) {endit = 1;}
		pos++;
	}

	sbuf[0] = ch;

	endit = 0;
	int bufpos = 1;
	while(!endit) //read word
	{
		if (pos >= Length()) {pos = -1; endit = 1;}
		else
		{
			ch = str[pos];

			if (IsWordDelimiterIdentifier(ch)) 
			{
				endit = 1;
				if (incl_delimiter) 
				{
					sbuf[bufpos] = ch;
					bufpos++;
					pos++;
				}
			}
			else
			{
				if (bufpos >= 100) //100 !!!!!!!!!!!!
				{
					sbuf[bufpos] = (char)0;
					bufpos = 0;
					s += sbuf;
				}
				sbuf[bufpos] = ch;
				bufpos++;
				pos++;
			}
		}
	}
	sbuf[bufpos] = (char)0;
	s += sbuf;
	return s;
}

//reads string in brackets, counting opening and closing brackets, -1 if no success
//works for infinitely many levels (stack limit)
//startpos must point at opening bracket
//pos finally points at position after cb
//also works for quotation marks "text", but only one level ...
mystr mystr::GetStringInBrackets(int& startpos, char ob, char cb) const
{
	int pos = startpos;
	int bracketcnt = 1;

	if (str[pos] != ob) 
	{
		startpos = -1;
		return mystr();
	}
	pos++;

	while (pos < Length() && bracketcnt != 0)
	{
		if (str[pos] == cb) bracketcnt--;
		else if (str[pos] == ob) bracketcnt++;
		pos++;
	}
	if (bracketcnt != 0) 
	{
		startpos = -1;
		return mystr();
	}

	int pos1 = startpos;
	startpos = pos;
	return SubString(pos1+1,pos-2);
}

char mystr::PosGet(int& pos) const //returns char and increments position 'pos' (or sets -1 if end of string)
{
	char c = str[pos++];
	if (pos >= Length()) pos = -1;
	return c;
}
char mystr::PosPeek(int& pos) const //returns char at position 'pos'
{
	return str[pos];
}
int mystr::GetUntil(int& pos, char until, mystr& s, int incl_delimiter) const //returns string between 'pos' and char 'until', increments pos to after until, sets -1 if end of string
{
	int rv = 1;
	s = "";
	int endit = 0;
	char sbuf[1002];
	char ch=(char)0;

	int bufpos = 0;
	while(!endit) 
	{
		if (pos >= Length()) 
		{
			pos = -1; 
			endit = 1;
			rv = 0;
		}
		else
		{
			ch = str[pos];

			if (ch == until) 
			{
				endit = 1;
				if (incl_delimiter) 
				{
					sbuf[bufpos] = ch;
					bufpos++;
					pos++;
				}
			}
			else
			{
				if (bufpos >= 1000) //1000 !!!!!!!!!!!!
				{
					sbuf[bufpos] = (char)0;
					bufpos = 0;
					s += sbuf;
				}
				sbuf[bufpos] = ch;
				bufpos++;
				pos++;
			}
		}
	}
	sbuf[bufpos] = (char)0;
	s += sbuf;
	return rv;
}

int mystr::GetUntilEOL(int& pos, char until, mystr& s, int incl_delimiter) const //returns string between 'pos' and char 'until', increments pos to after until, sets -1 if end of string; if end of file return 0, if until found return 2, else return 1
{
	int rv = 1;
	s = "";
	int endit = 0;
	char sbuf[1002];
	char ch=(char)0;

	int bufpos = 0;
	while(!endit) 
	{
		if (pos >= Length()) 
		{
			pos = -1; 
			endit = 1;
			rv = 0;
		}
		else
		{
			ch = str[pos];

			if ((ch == until && until != (char)0) || IsEOL(ch)) 
			{
				endit = 1;
				if (incl_delimiter) 
				{
					sbuf[bufpos] = ch;
					bufpos++;
					pos++;
				}
				if ((ch == until && until != (char)0)) rv = 2;
			}
			else
			{
				if (bufpos >= 1000) //1000 !!!!!!!!!!!!
				{
					sbuf[bufpos] = (char)0;
					bufpos = 0;
					s += sbuf;
				}
				sbuf[bufpos] = ch;
				bufpos++;
				pos++;
			}
		}
	}
	sbuf[bufpos] = (char)0;
	s += sbuf;
	return rv;
}

int mystr::GoUntil(int& pos, char until) const //moves 'pos' till char 'until', increments pos to after until, sets -1 if end of string
{
	while (1)
	{
		if (str[pos] == until) 
		{
			pos++;
			if (pos>=Length()) pos = -1;
			return 1;
		}
		pos++;
		if (pos>=Length())
		{
			pos = -1;
			return 0;
		}
	}
	return 0;
}

void mystr::EraseSpaces()
{
	int i = 0;
	int pos = 0;
	for (i = 0; i < length; i++)
	{
		str[pos] = str[i];
		//if (str[i] != ' ' && str[i] != (char)10 && str[i] != '\n') //old
		if (!IsSpace(str[i])) pos++;
	}
	str[pos] = 0;
	length = pos;
}

void mystr::EraseSpacesHeadTail()
{
// tail
	int pos = length-1;
	while (IsSpace(str[pos])) 
		pos--;
	length = pos+1;
	str[length] = 0;

//head
	pos = 0;
	while (IsSpace(str[pos]))
		pos++;
	for(int i=0; i <= length-pos; i++)
	  str[i] = str[i+pos];
	length = length - pos;

}

// returns position of first and last non-space letter (1-based)
void mystr::IdentifySpacesHeadTail(int& first, int& last)
{
// tail
	last = length;
	while (IsSpace(str[last-1])) 
		last--;
//head
	first = 1;
	while (IsSpace(str[first-1]))
		first++;
}

void mystr::EraseChar(char c)
{
	int i = 0;
	int pos = 0;
	for (i = 0; i < length; i++)
	{
		str[pos] = str[i];
		if (str[i] != c)
			pos++;
	}
	str[pos] = 0;
	length = pos;
}

void mystr::EraseChar(int pos) // 1 based ! 
{
	int i = 0;
	for (i = pos; i < length; i++)
	{
		str[i-1] = str[i];
	}
	str[length-1] = 0;
	length--;
}
int mystr::Replace(const mystr& searchstr, const mystr& replacestr) //find string searchstr, replace with string replacestr; this is done repeatedly; the return value counts the number of replacements
{
	int slen = searchstr.Length();
	int rlen = replacestr.Length();
	int pos = 0;
	int cnt = 0;
	while (pos < Length() && pos != -1)
	{
		pos = Find(pos, searchstr);
		if (pos != -1)
		{
			cnt++;
			for (int i=1; i<=slen; i++)
			{
				if (pos < Length()) EraseChar(pos+1);
			}
			InsertAt(pos, replacestr);
			pos += rlen;
			if (pos >= Length()) pos = -1;
		}
	}
	return cnt;
}

mystr& mystr::InsertAt(int pos, const mystr& s)
{
	if(pos > length)
	{
		mystr::ErrHandler();
		return *this;
	}
	int newLength = length + s.length;
	char *tmp = new char[newLength + 1];
	strncpy(tmp, str, pos);
	strcpy(tmp + pos, s.str);
	strcpy(tmp + pos + s.length, str + pos);
	delete[] str;
	length = newLength;
	str = tmp;
	return *this;
}

mystr &mystr::WriteAt(int pos, const mystr& s)
{
	if(pos > length)
	{
		mystr::ErrHandler();
		return *this;
	}
	int n = length - pos;
	if(s.length < n)
		n = s.length;
	strncpy(str + pos, s.str, n);
	return *this;
}

void mystr::SmartDouble2String(double x, int nDigits)
{
	delete[] str;

	char buffer[32];
	char buffer2[32];

	if (fabs(x) < 1E-100) {x = 0;}
	sprintf_s(buffer, "%.14g", x);
	int slen1 = strlen(buffer);
	sprintf_s(buffer2, "%.16g", x);
	int slen2 = strlen(buffer2);

	if(nDigits) //$ DR 2013-03-11: added in order to control the significant digits
	{
		mystr format = mystr("%.")+mystr(nDigits)+("g");
		sprintf_s(buffer,format,x);
		length = strlen(buffer);
		str = new char[length + 1];
		strcpy(str, buffer);
	}
	else	// original code
	{
		if (slen1 + 4 < slen2)
		{
			length = slen1;
			str = new char[length + 1];
			strcpy(str, buffer);
		}
		else
		{
			length = slen2;
			str = new char[length + 1];
			strcpy(str, buffer2);
		}
	}
}

// $EK - 2013-01-17: added in order to resize (including memory reallocation) the string
void mystr::ReSize(int len)
{
	delete[] str;
	length = len;
	str = new char[len + 1];
	str[len] = (char)0;
}

mystr& mystr::operator = (const mystr& s)
{
	if (&s != this)
	{
		delete[] str;
		length = s.length;
		str = new char[length + 1];
		strcpy(str, s.str);
	}
	return *this;
}

mystr operator + (const mystr& s1, const mystr& s2)
{
	mystr tmp(s1.length + s2.length, 0);
	if (s1.length != 0) strcpy(tmp.str, s1.str);
	if (s2.length != 0) strcpy(tmp.str + s1.length, s2.str);
	return tmp;
}

void mystr::operator += (const mystr& s)
{
	char *tmp = new char[length + s.length + 1];
	if (length != 0) strcpy(tmp, str);
	if (s.length != 0) strcpy(tmp + length, s.str);
	length += s.length;
	delete[] str;
	str = tmp;
}

char& mystr::operator [] (int n)
{
	static char dummy;
	if(n < length)
		return str[n];
	else
	{
		mystr::ErrHandler();
		return dummy;
	}
}

const char& mystr::operator [] (int n) const
{
	static char dummy;
	if(n < length)
		return str[n];
	else
	{
		mystr::ErrHandler();
		return dummy;
	}
}

mystr mystr::operator () (int l, int r) const
{
	if((l > r) || (r > length))
	{
		mystr::ErrHandler();
		mystr s;
		return s;
	}
	else
	{
		int n = r - l + 1;
		mystr tmp(n, 0);
		strncpy(tmp.str, str + 1, n);
		return tmp;
	}
}

istream& operator >> (istream& is, mystr& s)
{
	const int buflen = 1000;
	char buffer[buflen+1];

	int end = 0;
	s = "";
	mystr str;

	while (!end)
	{
		is.get(buffer, buflen);
		str = mystr(buffer);
		s += str;
		if (is.peek() == EOF) {end = 1;}
	}

	return is;
}
// checks if String can contains a number
// removes spaces around first 'E' or 'e' found in string 
// format of valid number:  [sign][digits][comma][digits][scientific: [ ] E|e [ ][sign]digits ]
// optional blanks before and after E in scientific notation
// return values:  =1: is a valid number; =0: not a valid number 

int mystr::IsValidNumber()
{
	this->EraseSpacesHeadTail();
// remove all spaces around 'E' and 'e'
	int pos = Find('E');
	if(pos==-1) pos = Find('e');
	
	if(pos!=-1) // found an 'E' or an 'e'
	{
		while (str[pos-1] == ' ') 
		{
			EraseChar(pos); // this function is 1 based
			pos--;
		}
		while (str[pos+1] == ' ')
		{
			EraseChar(pos+2);
		}
	}

// positions at 0 for not found
	int pos_dot = 0;      // one dot allowed
	int pos_sign1 = 0;    // leading sign
	int pos_sign2 = 0;    // sign at floating point exponent
	int pos_digit = 0;
	int pos_digitexp = 0;
	int begin_number = 0;
	int begin_exponent = 0;

	for(int i=1; i<=length; i++)
	{
		if(! begin_number) // begin of number
		{
			if (IsSign(str[i-1]))  { begin_number = i; pos_sign1 = i; }
			else if (IsNum(str[i-1]))  { begin_number = i; pos_digit = i; }
			else if (IsPoint(str[i-1])) { begin_number = i; pos_dot = i; }
			else return 0; // first non-space character is not a valid begin of a number
		}
		else if (!begin_exponent) // main part of number, before "E" or "e" of schientific notaion
		{
		 	if(IsNum(str[i-1])) { pos_digit = i; }
			else if(IsPoint(str[i-1])) // <-- point allowed only once
			{
				if (!pos_dot) pos_dot = i;
				else return 0; // decimal point for a 2nd time
			}
			else if (IsScientificE(str[i-1])) begin_exponent = i;
			else return 0; // any other character is not valid
		}
		else // if (begin_exponent)  // (optional) scientific notation
		{
			if (IsSign(str[i-1])) // <-- sign in exponent allowed only once
			{
				if (!pos_sign2 && !pos_digitexp) pos_sign2 = i;
				else return 0; // sign in exponent for a 2nd time
			}
			else if (IsNum(str[i-1])) pos_digitexp = i;
			else return 0;
		}
	}
	if (! begin_number) return 0; // could be empty string
	if (begin_exponent && !pos_digitexp) return 0; // could end with an "E" or "E-" ...
	return 1;
}




// check if the string holds a correct number (one-based)
// returns 1 if correct number was found & passes position-limits of number in integers first and last
// returns 0 otherwise
int IsCorrectNumberBetween(mystr num, int& first, int& last)
{
	num.IdentifySpacesHeadTail(first,last);
// positions at 0 for not found
	int pos_dot = 0;      // one dot allowed
	int pos_sign1 = 0;    // leading sign
	int pos_sign2 = 0;    // sign at floating point exponent
	int pos_digit = 0;
	int pos_digitexp = 0;
	int begin_number = 0;
	int begin_exponent = 0;


	for(int i=first; i<=last; i++)
	{
		if(! begin_number) // begin of number
		{
			if (IsSign(num[i-1]))  { begin_number = i; pos_sign1 = i; }
			else if (IsNum(num[i-1]))  { begin_number = i; pos_digit = i; }
			else if (IsPoint(num[i-1])) { begin_number = i; pos_dot = i; }
			else return 0; // first non-space character is not a valid begin of a number
		}
		else if (!begin_exponent) // main part of number, before "E" or "e" of schientific notaion
		{
		 	if(IsNum(num[i-1])) { pos_digit = i; }
			else if(IsPoint(num[i-1])) // <-- point allowed only once
			{
				if (!pos_dot) pos_dot = i;
				else return 0; // decimal point for a 2nd time
			}
			else if (IsScientificE(num[i-1])) begin_exponent = i;
			else return 0; // any other character is not valid
		}
		else // if (begin_exponent)  // (optional) scientific notation
		{
			if (IsSign(num[i-1])) // <-- sign in exponent allowed only once
			{
				if (!pos_sign2 && !pos_digitexp) pos_sign2 = i;
				else return 0; // sign in exponent for a 2nd time
			}
			else if (IsNum(num[i-1])) pos_digitexp = i;
			else return 0;
		}
	}
	if (! begin_number) return 0; // could be empty string
	if (begin_exponent && !pos_digitexp) return 0; // could end with an "E" or "E-" ...
	return 1;
}

int mystr::Trim_And_Get_Indices(int& comp1, int& comp2) // identifies and parses '[..,..]', shortens the string !! returns number of found components
{
	comp1=0; comp2=0;
  mystr buff1(32,32);
  mystr buff2(32,32);

	int nr_of_components = Trim_And_Get_Indices(buff1,buff2);
	if (nr_of_components==1 || nr_of_components==2) comp1 = buff1.MakeInt();
	if (nr_of_components==2) comp2 = buff2.MakeInt();

	return nr_of_components;

 // int pos_of_open_bracket = Find('[');
	//if (pos_of_open_bracket == -1) return 0;   // no bracket found

	//int pos_of_close_bracket = Find(']');
 // if (pos_of_close_bracket == -1) return -1; // error

	//// identify first component
	//int pos=pos_of_open_bracket+1;
 // mystr buffer = GetWord(pos);
	//comp1 = buffer.MakeInt();
	//int nr_of_components = 1;

	//int pos_of_comma = Find(',');
	//if (pos_of_comma > pos_of_open_bracket && pos_of_comma < pos_of_close_bracket && pos_of_comma >= pos)
	//{
	//	// identify second component
	//	pos = pos_of_comma+1;
	//	buffer = GetWord(pos);
	//	comp2 = buffer.MakeInt();
	//	nr_of_components = 2;
	//}
 // 
	//// trim the string
	//SetLength(pos_of_open_bracket);
	//this->str[pos_of_open_bracket] = 0;

	//return nr_of_components;
}

//////// DOES NOT WORK WITH MORE THEN ONE BRACKET a[b[1,2]]
//////int mystr::Trim_And_Get_Indices(mystr& comp1, mystr& comp2) // identifies and parses '[..,..]', shortens the string !! returns number of found components
//////{
//////	comp1=0; comp2=0;
//////  int pos_of_open_bracket = Find('[');
//////	if (pos_of_open_bracket == -1) return 0;   // no bracket found
//////	
//////
//////	int pos_of_close_bracket = Find(']');
//////  if (pos_of_close_bracket == -1) return -1; // error
//////
//////	// identify first component
//////	int pos=pos_of_open_bracket+1;
//////  comp1 = GetWord(pos);
//////	int nr_of_components = 1;
//////
//////	int pos_of_comma = Find(',');
//////	if (pos_of_comma > pos_of_open_bracket && pos_of_comma < pos_of_close_bracket && pos_of_comma >= pos)
//////	{
//////		// identify second component
//////		pos = pos_of_comma+1;
//////		comp2 = GetWord(pos);
//////		nr_of_components = 2;
//////	}
//////  
//////	// trim the string
//////	SetLength(pos_of_open_bracket);
//////	this->str[pos_of_open_bracket] = 0;
//////
//////	return nr_of_components;
//////}

int mystr::Trim_And_Get_Indices(mystr& comp1, mystr& comp2)
{
	comp1=0; comp2=0;

  int nr_open_brac = 0;        // number of open brackets
	int pos_of_first_open_bracket = length;
	int comp1_start = 0;         // position of beginning of component 1 in full string
	int comp1_end = 0;           // position of end on component 1 in full string
	int comp2_start = 0;
	int comp2_end = 0;

	int pos=0;
	int nr_comp=0;
	while (pos<length)
	{
		if (str[pos]=='[') nr_open_brac++;
		if (str[pos]==']') nr_open_brac--;

		if (nr_open_brac==1 && str[pos]=='[') { comp1_start = pos+1; nr_comp = 1; pos_of_first_open_bracket = pos;}
		if (nr_open_brac==1 && str[pos]==',') { comp1_end = pos-1, nr_comp = 2; comp2_start = pos+1; }
		if (nr_open_brac==0 && str[pos]==']') { if(nr_comp==1) comp1_end = pos-1; if(nr_comp==2) comp2_end = pos-1; }

		pos++;
	}
	if (nr_comp==1 || nr_comp==2) comp1 = this->SubString(comp1_start,comp1_end);
	if (nr_comp==2) comp2 = this->SubString(comp2_start,comp2_end);

// trim the string
	SetLength(pos_of_first_open_bracket);
	this->str[pos_of_first_open_bracket] = 0;

	return nr_comp;
}

int mystr::CountLines()
{
	//count lines:
	mystr line;
	int pos = 0;
	int eof = 0;
	int total_lines = 0;
	while (!eof && pos != -1)
	{
		total_lines++;
		if (GetUntilEOL(pos, (char)0, line, 1)) eof = 1;
		//pos++;
	}
	return total_lines;
}