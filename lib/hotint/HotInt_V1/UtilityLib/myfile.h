//#**************************************************************
//#
//# filename:             myfile.h
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

#ifndef MYFILE__H
#define MYFILE__H

typedef enum {TFMread = 1, TFMwrite = 2} TFileMode;

class CMFile
{
public:
	//open file namei in read or write mode
	CMFile(const mystr& namei, TFileMode fmi, int binary_mode = 0);
	void SetCMFile(const mystr& namei, TFileMode fmi, int binary_mode = 0);//used in constructor
	virtual ~CMFile() 
	{
		if (ofile) delete ofile;
		if (ifile) delete ifile;
	}
	int IsGood() const 
	{
		if (fm == TFMread)
		{
			return ifile->good()&&(!ifile->fail());
		} else
		{
			return ofile->good()&&(!ofile->fail());
		}
	}
	//read one mystr, separated by ' ', EOL, EOF
	virtual void RW(mystr& str);
	//read whole file as mystr
	virtual void RWF(mystr& str);

	//read or write mystr, but also compare, if read mystr is right
	int CRW(const mystr& str);
	int CRW(const mystr& str, mystr& strreceived);
  
	//$ RL 2012-1-17: [
	// return value is zero after succeccfully reading two columns from a data file
	// parameter description:
	//   data_col[1|2]...column data is stored here after reading from file, ncolumns...number of columns, first column   , second column    , number of header lines    , offset[1|2]...additional data value, second offset, offset[1|2]_start_index: after these numbers of lines read, the offset[1|2] is activated
	int ReadTwoColumnsFromFile(TArray<double>& data_col1, TArray<double>& data_col2, const int ncolumns, const int column1, const int column2, const int nrOfHeaderLines, const double offset1=0., const double offset2=0.   , const int offset1_start_index=1, const int offset2_start_index=1);
	//$ RL 2012-1-17: ]
	//read n characters
	void RWchars(int n, mystr& str);
	void RWchar(char& ch);
	void RWSoleStr(mystr& str);
	void RWSoleStrBrackets(mystr& str);

	void RWuntil(char until, mystr& str, int withlastchar = 1);
	void RWuntil(char increaseuntil, char until, mystr& str, int withlastchar = 1);
	void RWuntilEOL(mystr& str, int withlastchar = 1);
	
	void RWInt(int& i);
	int GetRWInt() { int i; RWInt(i); return i; };
	void RWDouble(double& d);
	double GetRWDouble() { double d; RWDouble(d); return d; };
	void SetPrecision(int prec) {if (ofile != NULL) ofile->precision(prec);}

	int IsEOL(char ch)
	{
		return (ch == (char)10 || ch == (char)13 || ch == '\n');
	}

	int IsDelimiter(char d)
	{
		return (d == (char)10 || d == (char)13 || d == (char)12 || d == (char)11 || d == '\n'
			|| d == ' ' || d == ',' || d == ';' || d == '\t' //original from parser
			|| d == ':' || d == '{' || d == '}' || d == '='); //this line added for MBS-system
	}

	int EndOfFile();

	int Rewind();
	int GetStreamPosition();
	int SetStreamPosition(int pos);

	void RWbinaryFloat(float& val);
	void RWbinaryDouble(double& val);
	void RWbinaryInt(int& val);
	void RWbinaryByte(char& ch);

	mystr name;
	mystr& GetName() { return name; }
	mystr GetNameNoPath();
	TFileMode fm;

	ofstream* ofile;
	ifstream* ifile;
};

// 
// extension of CMFile to read matrix from file ( solution files )
class CMatrixFile : public CMFile
{
public: // lifecycle
	// no standard constructor in base class: //CMatrixFile() {}
	// no copy constructor in base class //CMatrixFile(const CMatrixFile& other) {}
	~CMatrixFile() {}

public: // lifecycle II
	CMatrixFile(const mystr& namei, TFileMode fmi, mystr commenttagi = mystr("%"), int binary_mode = 0):CMFile(namei,fmi,binary_mode)
	{
		CommentTag() = commenttagi;
		Reset(0);
		linebuffer = mystr();                                    
	}
	void Reset(int nr_columns)
	{
		for(int i=1; i<=columns.Length(); i++)
		{
			Column(i).Flush();
		}
		columns.SetLen(nr_columns);
		columns.SetAll(TArray<double>(0));
	}

public: // access
	TArray<double>& Column(int i) { return columns(i); }
	TArray<double>* ColumnPtr(int i) { return &columns(i); }
	mystr& ColumnName(int i) { return columnnames(i); }
	int NColumns() { return columns.Length(); }
	int AddToCol(double val, int colnr) { return Column(colnr).Add(val); }
	mystr& CommentTag() { return commenttag; }

public: // short function calls
	int ReadSingleColumn(int colnr);                             // reads a single column from the file
	int ReadAllColumns();                                        // reads all columns from the file
	int ReadSpecifiedColumns(TArray<int>& columns);              // reads specified columns from the file
  int ReadHeader();                                            // loads all commented lines to header only to string, also counts headerlines
	int AppendRead();                                            // continues to read, assumes that header is read

protected: // functionality
	int ReadColumns(TArray<int>& in_specified_cols_sorted);      // read all columns specified in TArray - !!! assume that the column numbers are sorted ascendingin_specified_cols_sorted !!!
	int ParseHeader();

protected: // line parse functions
	int MystrToColumns(mystr& line);                             // writes the values from a single line mystr to the arrays
	int StreamToColumns();                                       // writes the values from a single line of the stream to the arrays

protected: // helper functions
	int IsCommentLine(mystr& line = mystr());                    // determine if line is a comment
  int CountColumns(mystr& line = mystr());                     // counts number of columns in first data line
  int ReadColumnNames();                                       // extracts column names from last header line ( if such a line exists )
	int LinesRead();                                             // returns number of (data) lines that were completely read into the arrays

private: // variables
	TArrayDynamic<TArray<double>> columns;                       // holds double values, stored column-wise, if column not read the array is empty
	TArrayDynamic<mystr> columnnames;
	mystr commenttag;
	TArray<int> specified_cols;                                  // stores the specified columns (for append)

	mystr headerstring;                                          // buffer for header ( only commented lines, no data lines included )
  mystr linebuffer;                                            // all lines are read to this buffer, can be used by several subroutines 

	int nr_headerlines;
	int nr_columns;
 	int lastvalidpos;
};

#endif
