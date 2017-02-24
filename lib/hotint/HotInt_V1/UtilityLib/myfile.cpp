//#**************************************************************
//#
//# filename:             myfile.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            11.06.97
//# description:          Class for Read and Write operations
//# remarks:						  do not try to understand the sense of every function
//#											  it was designed for a symbolic algebra program ...
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

//#include "../WorkingModule/stdafx.h" //for comment out
//#include "../WorkingModule/stdafx.h" 
//#include <stdafx.h> //$!AD 2011-04-01 changed line such that tool:data_reader can compile... in case of error: add "../WorkingModule" to include path

#include "ioincludes.h"

#include <string.h>
#include <math.h>

#include "mystring.h"
#include "tarray.h"
#include "myfile.h"


CMFile::CMFile(const mystr& namei, TFileMode fmi, int binary_mode)
{
	SetCMFile(namei, fmi, binary_mode);
}

void CMFile::SetCMFile(const mystr& namei, TFileMode fmi, int binary_mode) //(RL)
{
	ifile = NULL;
	ofile = NULL;
	name = namei;
	fm = fmi;
	if (fm == TFMread)
	{
#ifdef my_new_stdiostream
		int opt = ios::in;
		if (binary_mode) opt = opt|ios::binary;

		ifile = new ifstream(name.c_str(), opt);
#else
		int opt = ios::in|ios::nocreate;
		ifile = new ifstream(name.c_str(), opt, filebuf::sh_read);
#endif
	} 
	else
	{
		int opt = ios::out;
#ifdef my_new_stdiostream
		if (binary_mode) opt = opt|ios::binary;
		ofile = new ofstream(name.c_str(), opt);
#else
		ofile = new ofstream(name.c_str(), opt, filebuf::sh_write); //filebuf::sh_read corrected to filebuf::sh_read ==> not tested!
#endif
	}
}
//read one mystr, separated by ' ', EOL, EOF
void CMFile::RW(mystr& str)
{
	if (fm == TFMread) 
	{
		(*ifile) >> str; 
	} else
	{
		(*ofile) << str;
	}
}
//read whole file as mystr
void CMFile::RWF(mystr& str)
{
	if (fm == TFMread)
	{
		str="";
		char ch=' ';  //dummy char, not EOF
		int cnt=0;
		const int buflen = 10000;
		char buf[buflen+1];

		while (!ifile->eof())
		{
			int i=0;
			while (!ifile->eof() && i < buflen)
			{

				cnt++;
				ifile->get(ch);
//!AD 2012-02-10: [bugfix: force the string to end with a '\n'  
				if(ifile->gcount()!=0)       
				{				
					buf[i] = ch;            // original line - if get(ch) was not successful the same character was written twice
				}
				else
				{
					buf[i] = '\n';
				}
//!AD bugfix]
				i++;
			} 
			buf[i] = (char)0;
			str += mystr(buf);  
		}

		//cout << str;

	} else
	{
		(*ofile) << str;
	}
}

int CMFile :: CRW(const mystr& str)
{
	if (fm == TFMread)
	{
		mystr str2="";
		char ch=' ';
		while (ch != (char)EOF && IsDelimiter(ch))
		{
			(*ifile) >> ch;
		}
		str2 += ch;
		for (int i = 0; i < str.Length()-1; i++)
		{
			if (ch != (char)EOF)
			{
				(*ifile) >> ch;
				str2 += ch;
			}
		}
		if (str2 == str) {return 1;}
		else {return 0;}

	} else
	{
		(*ofile) << str;
		return 1;
	}
}

int CMFile :: CRW(const mystr& str, mystr& strreceived)
{
	if (fm == TFMread)
	{
		strreceived="";
		char ch=' ';
		while (ch != (char)EOF && IsDelimiter(ch))
		{
			(*ifile) >> ch;
		}
		strreceived += ch;
		for (int i = 0; i < str.Length()-1; i++)
		{
			if (ch != (char)EOF)
			{
				(*ifile) >> ch;
				strreceived += ch;
			}
		}
		if (strreceived == str) {return 1;}
		else {return 0;}

	} else
	{
		(*ofile) << str;
		return 1;
	}
}

int CMFile :: ReadTwoColumnsFromFile(TArray<double>& data_col1, TArray<double>& data_col2, const int ncolumns, const int column1, const int column2, const int nrOfHeaderLines, const double offset1, const double offset2, const int offset1_start_index, const int offset2_start_index)
{
	if (ifile!=0 && !mystr(name).Length() || fm != (TFileMode)TFMread || !IsGood())
	{
		// problems with global_uo --> commented out
		//if(!mystr(name).Length())global_uo << "No Filename for Elementsfile specified!\n";
		//else if(fm != (TFileMode)TFMread)global_uo << "File mode is not in read modus!\n";
		//else global_uo << "File not good!\n";
		return -1;
	}
  int rv = 0; //OK
	data_col1.SetLen(0);
	data_col2.SetLen(0);
	double val = -1e300;
	int lineCount = 0;
	if(IsGood())
	{
		while (!EndOfFile())
		{
			if(lineCount++ < nrOfHeaderLines)
			{
				// cut off header
				mystr s_tmp;
				RWuntilEOL(s_tmp, 1); // read line with last char
			}
			else
			{
				// read 
				double val1, val2;
				for(int col=1; col<=ncolumns && !EndOfFile(); col++)
				{
					RWDouble(val); // dummy
					if(val == -1e300)
					{
						rv = 1;
						return rv; // error
					}
					if(col == column1)
					{ 
						double offs = 0.;
						if(data_col1.Length() + 1 >= offset1_start_index)
						{
							offs = offset1; // consider offset at 
						}
						val1 = val+offs; // attention: no duplicate x-values should appear in file!!!
					}
					else
					{
						if(col == column2)
						{
							double offs = 0.;
							if(data_col2.Length() + 1 >= offset2_start_index)
							{
								offs = offset2;
							}
							val2 = val+offs;
						}
					}
				}
				if (!EndOfFile())
				{
					data_col1.Add(val1);
					data_col2.Add(val2);
				}			
			}
		}
	}
	else
	{		
		rv = 1; // error
	}
	return rv;
}


void CMFile :: RWchars(int n, mystr& str)
{
	if (fm == TFMread)
	{
		char ch=' ';  //dummy char, not EOF
		int cnt=0;
		char* buf = new char[n+1];

		while (!ifile->eof() && cnt < n)
		{
			ifile->get(ch);
			*(buf+cnt) = ch;
			cnt++;
		} 
		*(buf+cnt) = (char)0;
		str = mystr(*buf);

		//cout << str;

	}
	else
	{
// new AD - make sure to write exactly n characters
		int cnt = 1;
		int len = str.Length();
		if (len > n) len = n;

		while (cnt<=len)
		{
			(*ofile) << str.PosPeek(cnt);
			cnt++;
		}
		while (cnt<=n)
		{
			(*ofile) << ' ';
			cnt++;
		}
		//(*ofile) << str; // old
	}
}

void CMFile :: RWchar(char& ch)
{
	if (fm == TFMread)
	{
		ch=' ';
		if (!EndOfFile()) ifile->get(ch);
	} 
	else
	{
		(*ofile) << ch;
	}
}

void CMFile :: RWuntil(char until, mystr& str, int withlastchar)
{
	if (fm == TFMread)
	{
		str = "";
		char ch;
		ifile->get(ch);
		if (withlastchar) {str += ch;}

		while (!ifile->eof() && ch != until)
		{
			if (withlastchar) 
			{
				ifile->get(ch);
				str += ch;
			}
			else
			{
				str += ch;
				ifile->get(ch);
			}
		} 
		//cout << str;

	} 
	else
	{
		(*ofile) << str << until;
	}
}

//find char 'until'; for every 'increaseuntil', one more 'until' must be found
void CMFile :: RWuntil(char increaseuntil, char until, mystr& str, int withlastchar)
{
	int inc = 0;
	if (fm == TFMread)
	{
		str = "";
		char ch;
		ifile->get(ch);
		if (withlastchar) {str += ch;}

		while (!ifile->eof() && IsGood() && (ch != until || inc > 0))
		{
			if (ch == until && inc > 0) inc--; 
			if (ch == increaseuntil) inc++;

			if (withlastchar) 
			{
				ifile->get(ch);
				str += ch;
			}
			else
			{
				str += ch;
				ifile->get(ch);
			}
		} 
		//cout << str;

	} else
	{
		(*ofile) << str << until;
	}
}


void CMFile :: RWuntilEOL(mystr& str, int withlastchar)
{
	if (fm == TFMread)
	{
		str = "";
		char ch;
		ifile->get(ch);
		if (withlastchar) {str += ch;}

		while (!ifile->eof() && ch != (char)10 && ch != (char)13 && ch != '\n')
		{
			if (withlastchar)
			{
				ifile->get(ch);
				str += ch;
			}
			else
			{
				str += ch;
				ifile->get(ch);
			}
		} 
		//cout << str;
		if (ch == (char)13) 
		{
			ifile->get(ch);  //read carriage return code 10C
		}

	} 
	else
	{
		(*ofile) << str << endl;
	}
}


void CMFile :: RWSoleStr(mystr& str)
{
	const int withlastchar=0;
	if (fm == TFMread)
	{
		str = "";
		int endit = 0;
		char sbuf[101];
		char ch=0;
		while(!endit)
		{
			if (ifile->eof()) return;
			ifile->get(ch);
			if (!IsWordDelimiter(ch)) {endit = 1;}
		}
		sbuf[0] = ch;
		int pos = 1;

		endit = 0;
		while(!endit) 
		{
			if (ifile->eof()) {sbuf[pos] = 0; str = sbuf; return;}
			ifile->get(ch);
			if (IsDelimiter(ch)) 
			{
				endit = 1;
			}
			else
			{
				if (pos >= 100) //100 !!!!!!!!!!!!
				{
					sbuf[pos] = (char)0;
					pos = 0;
					str += sbuf;
				}
				sbuf[pos] = ch;
				pos++;
			}
		}
		sbuf[pos] = (char)0;
		str+=sbuf;
	} 
	else
	{
		(*ofile) << str << endl;
	}
}

void CMFile :: RWSoleStrBrackets(mystr& str) // same as above but excluded brackets (AD)
{
	const int withlastchar=0;
	if (fm == TFMread)
	{
		str = "";
		int endit = 0;
		char sbuf[101];
		char ch=0;
		while(!endit)
		{
			if (ifile->eof()) return;
			ifile->get(ch);
			if (!IsWordDelimiterBrackets(ch)) {endit = 1;}
		}
		sbuf[0] = ch;
		int pos = 1;

		endit = 0;
		while(!endit) 
		{
			if (ifile->eof()) {sbuf[pos] = 0; str = sbuf; return;}
			ifile->get(ch);
			if (IsWordDelimiterBrackets(ch)) 
			{
				endit = 1;
			}
			else
			{
				if (pos >= 100) //100 !!!!!!!!!!!!
				{
					sbuf[pos] = (char)0;
					pos = 0;
					str += sbuf;
				}
				sbuf[pos] = ch;
				pos++;
			}
		}
		sbuf[pos] = (char)0;
		str+=sbuf;
	} 
	else
	{
		(*ofile) << str << endl;
	}
}

void CMFile :: RWInt(int& i)
{
	if (fm == TFMread)
	{
		(*ifile) >> i;
	}
	else
	{
		(*ofile) << i << " ";
	}
}

void CMFile :: RWDouble(double& d)
{
	if (fm == TFMread)
	{
		(*ifile) >> d;
		/*		mystr dstr;
		RWSoleStr(dstr);
		//cout << " str='" << dstr << "'" << endl;
		d = dstr.MakeDouble();
		//cout.precision(17);
		//cout << "   d='" << d << "'" << endl;
		//cout.precision(6);
		*/
	}
	else
	{
		(*ofile) << d << " ";
	}
}

void CMFile::RWbinaryDouble(double& val)
{
	const int n = sizeof(double);
	int cnt=0;
	char buf[n];
	
	if (fm == TFMread)
	{
		ifile->read(buf,n);
		memcpy(&val,&buf,n);
	}
	else
	{
		memcpy(&buf,&val,n);
		ofile->write(buf,n);
	}
}

void CMFile :: RWbinaryFloat(float& val)
{
	if (fm == TFMread)
	{
		const int n = sizeof(float);
		int cnt=0;
		char buf[n];

		while (/*!ifile->eof() && */cnt < n)
		{
			ifile->get(buf[cnt++]);
		}
		memcpy(&val,&buf,n);
	}
	else
	{
		const int n = sizeof(float);
		int cnt=0;
		char buf[n];

		memcpy(&buf,&val,n);

		while (cnt < n)
		{
			ofile->put(buf[cnt++]);
		}
	}
}

void CMFile :: RWbinaryInt(int& val)
{
	if (fm == TFMread)
	{
		const int n = sizeof(int);
		int cnt=0;
		char buf[n];

		while (/*!ifile->eof() && */cnt < n)
		{
			ifile->get(buf[cnt++]);
		}
		memcpy(&val,&buf,n);
	}
	else
	{
		const int n = sizeof(float);
		int cnt=0;
		char buf[n];

		memcpy(&buf,&val,n);

		while (cnt < n)
		{
			ofile->put(buf[cnt++]);
		}
//		assert(0);
	}
}

void CMFile :: RWbinaryByte(char& ch)
{
	if (fm == TFMread)
	{
		ch=(char)0;
		ifile->read(&ch,1);
	} 
	else
	{
		(*ofile) << ch;
	}
}


int CMFile :: EndOfFile()
{
	if (fm == TFMread)
	{
		return ifile->eof();
	}
	else
	{
		return 0;
	}
}

int CMFile :: Rewind()
{
	if (fm = TFMread)
	{
		ifile->clear(); // resets all error-flags, most important: EOF flag
		return SetStreamPosition(0);
	}
	else
	{	
		ofile->clear();
		return SetStreamPosition(0);
	}
}

int CMFile :: GetStreamPosition()
{
	if (fm = TFMread)
	{
		return ifile->tellg();	
	}
	else
	{		
		return ofile->tellp();
	}
}

int CMFile :: SetStreamPosition(int pos)
{
	if (fm = TFMread)
	{
		ifile->seekg(pos, std::ios::beg);
		return GetStreamPosition();
	}
	else
	{		
		ofile->seekp(pos, std::ios::beg);
		return GetStreamPosition();
	}
}

mystr CMFile :: GetNameNoPath()
{
	int foundpos = 0;
	int pos = 0; 
	while (foundpos > -1)
	{
		foundpos = name.Find(pos,'\\');
		if(foundpos > -1) pos = foundpos+1;
	}
	mystr buffer;
	buffer = name.SubString(pos,name.Length()-1);

  return buffer;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// extension of CMFile to read matrix from file ( solution files )

// read a single column from the file
int CMatrixFile::ReadSingleColumn(int colnr)
{
	ParseHeader();
  return ReadColumns(IntVec1(colnr));
}

// read all columns
int CMatrixFile::ReadAllColumns()
{
	ParseHeader();
	return ReadColumns(NaturalNumbers(nr_columns));
}

int CMatrixFile::ReadSpecifiedColumns(TArray<int>& columns)
{
	ParseHeader();
	return ReadColumns(columns);
}

// continues to read the file - file is
int CMatrixFile::AppendRead()
{
	int lines_already_read = LinesRead();

// ASSUME: this function is called when the header etc has already been read, columns specified
	this->Rewind();

#define seekpos_rather_then_skiplines
#ifdef seekpos_rather_then_skiplines
  SetStreamPosition(lastvalidpos);
#else
	for(int i=1; i<= nr_headerlines; i++)
	{
		this->RWuntilEOL(linebuffer,1);
	}
	for(int i=1; i<= lines_already_read-1; i++)
	{
		this->RWuntilEOL(linebuffer,1);
	}
#endif

#define fromstream
	while( !EndOfFile() )
	{
#ifdef fromstream		
		StreamToColumns();
#else
		RWuntilEOL(linebuffer);
		MystrToColumns(linebuffer);
#endif
	}
	return Column(specified_cols.Last()).Length();
}

// read all columns specified in TArray - !!! assume that the column numbers are sorted ascendingin_specified_cols_sorted !!!
int CMatrixFile::ReadColumns(TArray<int> &in_specified_cols_sorted)
{
	specified_cols = in_specified_cols_sorted;
	if(specified_cols.Length() < 1) return 0;

// variable linebuffer contains the first non-commented line (from 
  MystrToColumns(linebuffer);

#define fromstream
	while( !EndOfFile() )
	{
#ifdef fromstream		
		StreamToColumns();
#else
		RWuntilEOL(linebuffer);
		MystrToColumns(linebuffer);
#endif
	}
	return Column(specified_cols.Last()).Length();
}

int CMatrixFile::ParseHeader()
{
	nr_headerlines = ReadHeader();                   // count number of headerlines, save header in buffer
	nr_columns = CountColumns();                     // count number of columns from first dataline
	Reset(nr_columns);                               // set number of columns
	ReadColumnNames();                               // column names are read or created as default

	return 0;
}

// writes the values from a single line mystr to the arrays
int CMatrixFile::MystrToColumns(mystr& line)
{
	if(line.Length() == 0) 
		line = linebuffer;

	int pos_line = 0;
	mystr word;
 	int spec_col_count = 1;
	int lastcol = specified_cols.Last();

	for (int i = 1; i<= lastcol; i++)                    // no need to parse entire line, only up to highest column
	{
		word = line.GetWord(pos_line,1);
		if (pos_line == -1)
		{
			return 0;
		}
		if( i == specified_cols(spec_col_count) )      
		{
			double val = word.MakeDouble();                 // this number must be written to file
			Column(i).Add(val);            
			spec_col_count++;
		}
	}
	return spec_col_count-1;
}

// writes the values from a single line of the stream to the arrays
int CMatrixFile::StreamToColumns()
{
	register double val;
	int spec_col_count = 1;
	int lastcol = specified_cols.Last();

	for (int i=1; i<= nr_columns; i++)
	{
		this->RWDouble(val);
		
		if(EndOfFile()) 
			return -1;                                     // reached end of file

		if( i == specified_cols(spec_col_count) )      
		{
			Column(i).Add(val);            
			if (spec_col_count < specified_cols.Length()) // reached last column that must be written, prevent TArray<T>::Resize 
				spec_col_count++;
		}
	}
// full row of the matrix was read
	lastvalidpos = GetStreamPosition();
	return spec_col_count-1;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper functions

int CMatrixFile::ReadHeader()
{
	headerstring = mystr();
  int linecount = 0;

	RWuntilEOL(linebuffer);
	while(IsCommentLine(linebuffer))
	{
	  headerstring += linebuffer;
		linecount++;
		RWuntilEOL(linebuffer);
	}
	return linecount;
}

int CMatrixFile::IsCommentLine(mystr& line)
{
	if(line.Length() == 0) 
		line = linebuffer;

	return (line.Left(CommentTag().Length()).Compare(CommentTag()));
}

// counts number of columns in first data line
int CMatrixFile::CountColumns(mystr& line)
{
	if(line.Length() == 0) 
		line = linebuffer;

	int pos=0;
	int numbercount=0;
	mystr number;
	while (pos!=-1)
	{
		number = line.GetWord(pos,1);
		if(pos!=-1)
			numbercount++;
	}
	return numbercount;
}

//////////// extracts column names from last header line
int CMatrixFile::ReadColumnNames()
{
	columnnames.Flush();

	if(nr_headerlines == 0) 
	{
		// no line containing column names
		for(int i=1; i<= nr_columns; i++)
			columnnames.Add(mystr("Column_")+mystr(i));
	}
	else
	{
		int pos=0;
		mystr line;
		for(int i=1; i<=nr_headerlines; i++)       // skip to last line of header
			headerstring.GetUntil(pos,'\n',line,1);

		line = line.Right(line.Length()-1); // remove "%"

		int pos_line=0;
		mystr word;
		while (pos_line!=-1)
		{
			word = line.GetWord(pos_line,1);
			if(pos_line!=-1)
				columnnames.Add(word);
		}
	}
	return columnnames.Length();
}

// returns number of (data) lines that were completely read into the arrays
int CMatrixFile::LinesRead()
{
  int linesread = 0;
	for (int i=1; i<=columns.Length(); i++)         // loop over all columns
	{
		int collength = Column(i).Length();
		if( collength != 0)                           // non empty column
		{
			if (!linesread)
			{
				linesread = collength;
			}
			else
			{
				if(collength <= linesread)
				{
  				linesread = collength;
				}
			}
		}
	}
	return linesread;
}

