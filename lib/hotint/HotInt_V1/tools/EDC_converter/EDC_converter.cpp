//#**************************************************************
//#
//# filename:             EDC_converter.cpp 
//#
//# author:               Gerstmayr, Reischl
//#
//# generated:						
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

#include "stdafx.h"
#include "math.h"
#include "..\..\UtilityLib\tarray.h"
#include "..\..\UtilityLib\mystring.cpp"
#include "..\..\UtilityLib\myfile.cpp"

//int _tmain(int argc, _TCHAR* argv[]) //this does not correctly pass arguments!
//{
//	return 0;
//}

//replace comma with dot
void ReplaceCommaStr(mystr& str)
{
	for (int j=1; j<=str.Length(); j++)
	{
		if (str[j-1] ==',') str[j-1] = '.';
	}
}

//erase double spaces from string
void EraseDoubleSpaces(mystr& str)
{
	int i = 0;
	int inew = 0;
	int lastspace = 0;
	while (i < str.Length())
	{
		str[inew] = str[i];
		if (!(lastspace && (str[i] == ' '))) 
		{
			inew++;
		}

		if (str[i] == ' ') 
		{
			lastspace = 1;
		}
		else 
		{
			lastspace = 0;
		}
		i++;
	}
	if (inew <= str.Length()-1) 
	{
		str[inew] = (char)0;
		str.SetLength(inew);
	}
}


int GetCommaSeparatedStrings(mystr str, TArrayDynamic<mystr>& substr_list, int showdetails)
{
	int lineok = 1;
	int searchcommands=1;
	int oldpos = 0;

	while (searchcommands)
	{
		int newpos = str.Find(oldpos,',');
		int stringstart = str.Find(oldpos,'"');
		if (stringstart != -1 && newpos != -1 && stringstart < newpos)
		{
			//there can be only one string between two commas ...
			int stringend = str.Find(stringstart+1,'"');
			newpos = stringend;
			if (stringend+1 < str.Length()) 
			{
				newpos = str.Find(stringend+1,',');
			}
			else 
			{
				newpos = -1;
			}
		}

		if (newpos != -1 && newpos>oldpos)
		{
		}
		else
		{
			newpos = str.Length();
			searchcommands = 0;
		}
		mystr substr = str.SubString(oldpos,newpos-1);
		substr_list.Add(substr);

		oldpos = newpos+1;
		if (showdetails) cout << "  command = '" << substr << "'\n";
	}
	if (substr_list.Length() == 0) lineok = 0;

	return lineok;
}

//identify a parameter, e.g.:
//int x
//int& x
//const int& x
//Vector x
//TArray<int> x
//==>[virtual][static][mutable][const] [type_name] [&|*] [parameter_name]

struct ParameterIDstruct
{
	mystr parameter_name;
	mystr type_name;
	int is_static;
	int is_mutable;
	int is_virtual; //only for functions
	int is_const;
	int is_pointer;   //*
	int is_reference; //&
};

int IdentifyParameter(mystr str, ParameterIDstruct& parameter_id, int showdetails)
{
	parameter_id.parameter_name = "";
	parameter_id.type_name = "";
	parameter_id.is_static = 0;
	parameter_id.is_mutable = 0;
	parameter_id.is_virtual = 0; //only for functions
	parameter_id.is_const = 0;
	parameter_id.is_pointer = 0;   //*
	parameter_id.is_reference = 0; //&

	int pos = 0; //read line from beginning;

	mystr wordstr;
	str.ReadLeadingSpaces(pos);
	if (pos == -1) return 0;

	int redo = 1;
	while (redo)
	{
		redo = 0;
		wordstr = str.GetWord(pos);
		if (showdetails) cout << "  wordstr=" << wordstr << "\n";

		if(wordstr == mystr("virtual")) {parameter_id.is_virtual = 1;redo = 1;}
		if(wordstr == mystr("mutable")) {parameter_id.is_mutable = 1;redo = 1;}
		if(wordstr == mystr("const")) {parameter_id.is_const = 1;redo = 1;}
		if(wordstr == mystr("static")) {parameter_id.is_static = 1;redo = 1;}
	}
	parameter_id.type_name = str.GetWord(pos); //now type should be there
	if (showdetails) cout << "vartype='" << parameter_id.type_name << "'\n";

	if (pos == -1) return 0;

	int pos_var = pos; //position where variable starts
	parameter_id.parameter_name = str.GetWord(pos);
	if (showdetails) cout << "parameter name='" << parameter_id.type_name << "'\n";

	if (parameter_id.parameter_name == mystr("")) return 0;
	return 1;
}


int main(int argc, char* argv[])
{
	cout << "EDC converter\n";

	//convert arguments into dynamic array
	const int max_narg = 10000;
	int nargs = argc;
	if (nargs > max_narg) nargs = max_narg;

	TArrayDynamic<mystr> argstrs(max_narg+1);
	for (int i=1; i <= nargs-1; i++)
	{
		argstrs(i) = argv[i];
	}

	//++++++++++++++++++++++++++++++++++++++++++++
	//output help information
	if (argstrs.Length() == 0 || argstrs.Find(mystr("-help")))
	{
		cout << "\nusage:\n";
		cout << "EDC converter -filelist " << '"' << "name_of_file_with_cpp_and_h_files" << '"' 
			<< "\n    -destination " << '"' << "name_of_destination_file_with_EDC_commands.cpp" << '"' << " [options]\n";
		cout << "possible options:\n";
		cout << "-help: list possible argument list\n";
		cout << "-showdetails: print detailed output during parsing\n";
		//cout << "-average N: average over N lines (default=1)\n";
		cout << "\n";

		cout << "exiting ...\n";

		return 0;
	}

	//++++++++++++++++++++++++++++++++++++++++++++
	//search for filename of file, which contains filelist
	int src = argstrs.Find(mystr("-filelist"));
	//search for destination filename
	int dest = argstrs.Find(mystr("-destination"));
	if (!dest) 
	{
		cout << "ERROR: no destination file defined!\n"; 
		return 1;
	}
	if (!src) 
	{
		cout << "ERROR: no filename for filelist defined!\n"; 
		return 1;
	}

	int maxlines = 100000; //maximum number of lines that are parsed in order to avoid crashes
	int showdetails = (argstrs.Find(mystr("-showdetails")) != 0);

	mystr sourcefile = argstrs(src+1);
	mystr destfile = argstrs(dest+1);

	//mystr sourcefile = "test2.txt";
	//mystr destfile = "dest.txt";

	cout << "filelist='" << sourcefile << "'\n";
	cout << "destination='" << destfile << "'\n";

	//++++++++++++++++++++++++++++++++++++++++++++
	//read filelist containing .h or .cpp files for parsing for EDC-access functions
	CMFile filelist(sourcefile, TFMread, 0);

	//int n_files = filelist.GetRWInt();
	int n_files = 0;
	int comment = 1;
	mystr filelistpath = "";
	while (comment)
	{
		mystr readline;
		filelist.RWuntilEOL(readline, 0);
		if (readline[0] == '%')
		{
			if (showdetails) cout << "read comment: " << readline << "\n";
		}
		else
		{
			comment = 0;
			filelistpath = readline;
			filelist.RWuntilEOL(readline, 0);
			n_files = readline.MakeInt();
		}

	}
	cout << "read " << n_files << " files ...\n";

	// read all filenames
	TArrayDynamic<mystr> filenames;
	mystr filename;
	for (int i=1; i<=n_files; i++)
	{
		filelist.RWuntilEOL(filename, 0);
		filenames.Add(filename);
		//cout << "file name '" << filename << "'\n";
	}

	//++++++++++++++++++++++++++++++++++++++++++++
	//file for import in HOTINT:
	ofstream outfile(destfile.c_str());

	outfile << "//#**************************************************************\n";
	outfile << "//#\n";
	outfile << "//# filename:             EDC_converter.cpp \n";
	outfile << "//#\n";
	outfile << "//# author:               Gerstmayr, Reischl\n";
	outfile << "//#\n";
	outfile << "//# generated:						\n";
	outfile << "//# description:          \n";
	outfile << "//# remarks:						  \n";
	outfile << "//#\n";
	outfile << "//# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian\n";
	outfile << "//# Center of Competence in Mechatronics GmbH, Institute of Technical Mechanics at the \n";
	outfile << "//# Johannes Kepler Universitaet Linz, Austria. All rights reserved.\n";
	outfile << "//#\n";
	outfile << "//# This file is part of HotInt.\n";
	outfile << "//# HotInt is free software: you can redistribute it and/or modify it under the terms of \n";
	outfile << "//# the HOTINT license. See folder 'licenses' for more details.\n";
	outfile << "//#\n";
	outfile << "//# bug reports are welcome!!!\n";
	outfile << "//# WWW:		www.hotint.org\n";
	outfile << "//# email:	bug_reports@hotint.org or support@hotint.org\n";
	outfile << "//#***************************************************************************************\n\n\n";


	outfile << "//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfile << "//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfile << "//+  automatically generated file for EDC class data exchange, do not modify!   +\n";
	outfile << "//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfile << "//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
	outfile << "//DESCRIPTION OF $EDC$ usage:                                                    \n";
	outfile << "//* mark beginning of class: $EDC$[beginclass,classname='class_name',parentclassname='parent_class_name',addelementtypename='EDCelementname',texdescription=" << '"' << "..." << '"' << ",]     \n";
	outfile << "//* mark end of class: $EDC$[endclass,'classname']                               \n";
	outfile << "//* with 'addelementtypename', the class is made available for adding on file load or in the menu of HOTINT\n";
	outfile << "//* EDC folders start with upper case letter, variables with lower case letter!!!\n";
	outfile << "//                                                                               \n";
	outfile << "//variable access:                                                               \n";
	outfile << "//  + syntax: [//EDC|const|static|mutable] [typename] [C++ variable name] //$EDC$[varaccess, {option1, {option2}}]\n";
	outfile << "//             ==>//EDC means that this variable is commented out totally indicated by symbol '//EDC' ==> e.g. because it is defined in another class\n";
	outfile << "//  + this generates an automated (EDC) access to this variable                  \n";
	outfile << "//* possible options:                                                            \n";
	outfile << "//  + internal_varaccess: create internal access functions for use in sensors and/or modifiers \n";
	outfile << "//  + readonly: this variable can be only read (e.g. for debugging, but not modified \n";
	outfile << "//  + EDCvarname=" << '"' << "variable_name" << '"' << ": assign variable name, if EDC-variable name should be different from C++ variable name\n";
	outfile << "//  + EDCfolder=" << '"' << "folder_name" << '"' << ": assign subfolder name of EDC, e.g. " << '"' << "Graphics" << '"' << " or " << '"' << "Physics" << '"' << "\n";
	outfile << "//  + tooltiptext=" << '"' << "text" << '"' << "]: tooltip (help) text, which describes the variable; this could also be used to document the variable name in HOTINT\n";
	outfile << "//  + int_bool: treat integer value as boolean (with check box)\n";
	outfile << "//  + variable_length_vector: vector has variable length ==> user can choose how many values are set\n";
	outfile << "//  + func_par1 = ...: assign a number or evaluable expression which is used as function parameter for double or int access\n";
	outfile << "//  + condition = ...: evaluable C++ expression, e.g. 'IsConstraint()' or 'i < 4', which is used as condition if Get/SetElementData entry is used or not\n";
	outfile << "//  + maxval = ..., minval = ... : define maximum and minimum integer values which are allowed\n";
	outfile << "//  + vecstart = ..., vecend = ... : define starting and ending of a vector by means of values or C++ expressions, if only a sub-vector should be used\n";
	outfile << "//  + remove: remove variable from an EDC of the parent class in a derived class. The elementdata is not really removed but locked and given a new tooltiptext\n";
	outfile << "//  \n";
	outfile << "//  + e.g.: static int switch; //$EDC$[varaccess, EDCvarname=" << '"' << "switch_for_HOTINT" << '"' << ", EDCfolder=" << '"' << "Debug" << '"' << ", tooltiptext=" << '"' << "this switch is nonsense" << '"' << ", readonly]\n";
	outfile << "//  \n";
	outfile << "//function access:                                                               \n";
	outfile << "//  + syntax: [const|static|mutable|virtual] [typename] [C++ function name]() //$EDC[funcaccess, {option1, {option2}}]\n";
	outfile << "//  + same options as for variable\n";
	outfile << "//  \n";
	outfile << "//  \n";

	for (int i=1; i<=n_files; i++)
	{
		filename = filenames(i);
		outfile << "#include " << '"' <<"..\\" << filename << '"' << "\n";
	}

	mystr add_class_descriptions = ""; //additional code for the adding of class descriptions; this are descriptions for .tex files
	mystr add_parse_init_code = ""; //additional code for the EDC parser in MBS class, during initialization
	mystr add_parse_add_code = ""; //additional code for the EDC parser in MBS class, for adding/generating an element or object

	// variables for generating the bibliography
	TArrayDynamic<mystr> ref;
	ref.SetLen(0);
	TArrayDynamic<mystr> ref_label;
	ref_label.SetLen(0);
	mystr bibfile = destfile;
	int posElLib = bibfile.Find("ElementsLib");
	bibfile = bibfile.Left(posElLib);
	bibfile = bibfile +mystr("documentation\\EDCauto_documentation\\bibliography_auto.tex");
	cout << "bibfile =  " << bibfile << " \n";

	int elementtypecnt = 0; //counter over all elements==>for all files same counter

	//go through all files:
	for (int i=1; i<=n_files; i++)
	{
		//mystr filename;
		//filelist.RWuntilEOL(filename, 0);
		filename = filenames(i);
		//cout << "file name '" << filename << "'\n";
		filename = filelistpath + filename;

		cout << "------------------------------------\n";
		cout << "parse file '" << filename << "'\n";

		//open single file for parsing:
		CMFile file(filename, TFMread, 0);

		mystr getEDstr; //string for the function GetElementData
		mystr setEDstr; //string for the function SetElementData
		mystr readSingle; //string for the function ReadSingleElementDataAuto
		mystr writeSingle; //string for the function WriteSingleElementDataAuto
		mystr availSingle; //string for the function IsAvailableSingleElementDataAuto
		mystr setEDstrRemove; //string for the function SetElementData, for removing entries of the edc

		mystr texDescr; //description of the class in tex-file
		mystr texDescrName; 

		//run through file and find EDC-identifiers
		int endfile = 0;
		int linecnt = 0;
		const int idstr_len = 5; //length of identifier string $EDC$
		mystr line;
		//maxlines = 1000;
		mystr actualclass=""; //store actual class name
		mystr parentclass=""; //store parent name of actual class 
		while (!endfile && linecnt++ <= maxlines)
		{
			file.RWuntilEOL(line, 0);
			//if (showdetails) cout << line << "\n";


			//+++++++++++++++++++++++++++++++++++
			//begin parse one line
			int edcpos = line.Find(mystr("$EDC$"));
			if (edcpos != -1)
			{
				int linelength = line.Length();
				if (showdetails) 
				{
					cout << "\nparse EDCline: " << line << ", ";
					cout << "EDC-pos=" << edcpos << "\n";
				}
				int lineok = 1;
				//search for commands after $EDC$:
				mystr commands=""; //
				if (edcpos + idstr_len >= linelength)
				{
					cout << "ERROR in line " << linecnt << ": no identifiers found after $EDC$ !\n";
					lineok = 0;
				}
				else
				{
					if (line[edcpos+idstr_len] != '[')
					{
						cout << "ERROR in line " << linecnt << ": no '[' found after $EDC$ !\n";
						lineok = 0;
					}
					else
					{
						//read commands in brackets:
						int pos = edcpos+idstr_len;
						int oldpos = pos;
						commands = line.GetStringInBrackets(pos,'[',']');
						
						//if opening and closing brackets go over several lines, also try to read consecutive lines
						const int maxbracketlines = 1000;
						int blcnt = 0;
						while (pos == -1 && blcnt++ < maxbracketlines)
						{
							mystr newline;
							file.RWuntilEOL(newline, 0);
							//line = line + mystr("\\n\\\n") + newline; //add the symbol "\n" for line-break in C++ string, a "\" for multi-line string, and the "\n" for the linebreak in the C++ code; all \ need to be doubled
							line = line /*+ mystr("\\n")*/ + newline; //$ 2013-02-06 new concept: line breaks are interpreted correctly, but not shown in ElementEDCauto.cpp

							pos = oldpos;
							commands = line.GetStringInBrackets(pos,'[',']');

							if (file.EndOfFile() || !file.IsGood()) blcnt = maxbracketlines;
						}
						if (blcnt > 0 && blcnt < maxbracketlines) 
						{
							commands.Replace("//", ""); //replace C++ comments
						}


						if (pos == -1) cout << "ERROR: no opening / closing bracket found in line" << linecnt << "\n";
						//cout << "pos=" << pos << "\n";

						//if (showdetails) cout << "commands=" << commands << "\n";
					}
				}
				//retrieve command list:
				TArrayDynamic<mystr> commandlist(0);
				lineok = lineok * GetCommaSeparatedStrings(commands, commandlist, showdetails);
				//cout << "commands=" << commands << "\n";

				if (lineok)
				{

					//line must be structured (symbol ";" must not be used in any text at the moment (only for separation of options!!!):
					//class A: public B //$EDC$[beginclass;A;B]  OR              (do not use space between $EDC$ and [
					//double x; //$EDC$[varaccess,tooltip="this is a test tooltiptext, which is for the variable x";readonly;minval=0;maxval=1e4]

					//cout << "text after $EDC$:" << line[edcpos+5] << "\n";
					//cout << "linelength=" << linelength << "\n";

					//trace class scope:
					//int fpos;
					int c_beginclass = commandlist.Find(mystr("beginclass"));
					int c_varacc = commandlist.Find(mystr("varaccess"));
					int c_funcacc = commandlist.Find(mystr("funcaccess")); //only non-parameter functions are supported at the moment: e.g. "virtual int SOS() const {...} //$EDC$[funcaccess, ...]"
					int c_endclass = commandlist.Find(mystr("endclass"));

					mystr addelementtypename = "";
					mystr addelementtype = "";
					mystr shortdescription = "";		//$ DR 2013-02-06 renamed from texdescription
					mystr texdescriptionDOF = "";
					mystr texdescriptionNode = "";
					mystr texdescriptionGeometry = "";
					mystr texdescriptionLimitations = "";
					mystr texdescriptionEquations = "";
					mystr texdescriptionExample = "";
					mystr texdescriptionComment = "";
					TArrayDynamic<mystr> figure;
					figure.SetLen(0);
					TArrayDynamic<mystr> caption;
					caption.SetLen(0);
					TArrayDynamic<mystr> modus;
					modus.SetLen(0);
					int ref_count = 0;

					if (c_beginclass)
					{
						for (int i=2; i <= commandlist.Length(); i++)
						{
							mystr str = commandlist(i);
							int pos = 0;
							mystr comword = str.GetWord(pos); //either read single word, or assignment
							//cout << "word='" << comword << "', pos = " << pos << ", strlength=" << str.Length() << "\n";
							if (pos < str.Length() && pos != -1) //==> read assignment; if pos is not at end of string and pos is not -1 (end of string reached), then this must be an assignment
							{
								if (str[pos] != '=')
								{
									cout << "ERROR: expected '=' after command, but received '" << str[pos] << "' in line " << linecnt << "\n";
								}
								pos ++;
								//cout << "pos=" << pos << ", str.Length()-1=" << str.Length()-1 << "\n";
								mystr rhs = str.SubString(pos, str.Length()-1);

								if (comword == mystr("classname"))
								{
									//mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
									actualclass = rhs;
									if (showdetails) cout << "    classname = '" << rhs << "'\n";
								}
								if (comword == mystr("parentclassname"))
								{
									parentclass = rhs;
									if (showdetails) cout << "    parentclassname = '" << rhs << "'\n";
								}
								if (comword == mystr("addelementtypename"))
								{
									//mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
									addelementtypename = rhs;
									if (showdetails) cout << "    addelementtypename = '" << rhs << "'\n";
								}
								if (comword == mystr("addelementtype"))
								{
									//mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
									addelementtype = rhs;
									if (showdetails) cout << "    addelementtype = '" << rhs << "'\n";
								}
								if (comword == mystr("texdescription"))
								{
									shortdescription = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    shortdescription = '" << shortdescription << "'\n";
								}
								if (comword == mystr("texdescriptionDOF"))
								{
									texdescriptionDOF = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionDOF = '" << texdescriptionDOF << "'\n";
								}
								if (comword == mystr("texdescriptionNode"))
								{
									texdescriptionNode = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionNode = '" << texdescriptionNode << "'\n";
								}
								if (comword == mystr("texdescriptionGeometry"))
								{
									texdescriptionGeometry = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionGeometry = '" << texdescriptionGeometry << "'\n";
								}
								if (comword == mystr("texdescriptionEquations"))
								{
									texdescriptionEquations = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionEquations = '" << texdescriptionEquations << "'\n";
								}
								if (comword == mystr("texdescriptionLimitations"))
								{
									texdescriptionLimitations = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionLimitations = '" << texdescriptionLimitations << "'\n";
								}
								if (comword == mystr("texdescriptionComments"))
								{
									texdescriptionComment = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionComment = '" << texdescriptionComment << "'\n";
								}
								if (comword == mystr("example"))
								{
									texdescriptionExample = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    texdescriptionExample = '" << texdescriptionExample << "'\n";
								}
								if (comword == mystr("figure"))
								{
									mystr figurename = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    total figure string = '" << figurename << "'\n";
									TArrayDynamic<mystr> captionname;
									GetCommaSeparatedStrings(figurename,captionname,showdetails);
									if (showdetails) cout << "    figure = '" << captionname(1) << "'\n";
									if(captionname.Length()==1) {captionname.Add(mystr(addelementtypename));}
									mystr caption_tmp = captionname(2);
									for(int c=3; c<=captionname.Length() ;c++)
									{
										caption_tmp = caption_tmp + mystr(",") + captionname(c);	// if there are "," in the caption, rebuild the string
									}
									if (showdetails) cout << "    caption = '" << captionname(2) << "'\n";
									figure.Add(captionname(1));
									caption.Add(caption_tmp);
								}
								if (comword == mystr("reference"))
								{
									mystr reference = rhs.SubString(1, rhs.Length()-2);
									ref.Add(reference);
									ref_count++;
									if (showdetails){ cout << "    ref.Length() = " << ref.Length() << ", ref_count = " << ref_count << ", reference = '" << reference << "'\n";}
								}
								if (comword == mystr("modus"))
								{
									mystr modusname = rhs.SubString(1, rhs.Length()-2);
									if (showdetails) cout << "    modus = '" << modusname << "'\n";
									modus.Add(modusname);
								}
								//actualclass = commandlist(2);
								//cout << "begin class '" << actualclass << "'\n";
								//if (commandlist.Length() >= 3)
								//{
								//	parentclass = commandlist(3);
								//	cout << "    parent class '" << parentclass << "'\n";
								//}

							}
							else //flag without assignment
							{
							}
						}//end for loop of commands
						if (actualclass == mystr(""))
						{
							cout << "ERROR: no classname found in beginclass statement in line " << linecnt << "\n";
						}
						//fill in all element data for begin class arguments:
						getEDstr += mystr("void ") + actualclass + mystr("::GetElementDataAuto(ElementDataContainer& edc)\n{\n");
						if (parentclass != mystr(""))
						{
							getEDstr += mystr("  ") + parentclass + mystr("::GetElementDataAuto(edc);\n");
						}
						getEDstr += mystr("  ElementData ed;\n");

						//setEDstr += mystr("int ") + actualclass + mystr("::SetElementDataAuto(const ElementDataContainer& edc)\n{\n");
						setEDstr += mystr("int ") + actualclass + mystr("::SetElementDataAuto(ElementDataContainer& edc)\n{\n");
						if (parentclass != mystr(""))
						{
							setEDstr += mystr("  ") + parentclass + mystr("::SetElementDataAuto(edc);\n");
						}
			
						readSingle += mystr("int ") + actualclass + mystr("::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable \n{\n");
						if (parentclass != mystr(""))
						{
							readSingle += mystr("  ") + mystr("int rv = ") + parentclass + mystr("::ReadSingleElementDataAuto(RWdata);\n");
							readSingle += mystr("  ") + mystr("if (rv==1) return 1;\n");
						}

						writeSingle += mystr("int ") + actualclass + mystr("::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable\n{\n");
						if (parentclass != mystr(""))
						{
							writeSingle += mystr("  ") + mystr("int rv = ") + parentclass + mystr("::WriteSingleElementDataAuto(RWdata);\n");
							writeSingle += mystr("  ") + mystr("if (rv==1) return 1;\n");
						}

						availSingle += mystr("int ") + actualclass + mystr("::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables\n{\n");
// AD: 2012-11-23:  recursion into base class is now in "manual" part
						//if (parentclass != mystr(""))
						//{
						//  availSingle += mystr("  ") + mystr("int rv = ") + parentclass + mystr("::GetAvailableSpecialValuesAuto(available_variables);\n");
						//}
						//else
						//{
						//  availSingle += mystr("  available_variables.Flush(); \n"); // flush the arrays at global parent
						//}

						int objFactory = 0;		//$ DR 0..not in ObjectFactory, 1..Element
						// in future: 2..Sensor,3..Node,...
						if (addelementtype.Find("TAEBody")!=-1) objFactory = 1;
						else if (addelementtype.Find("TAEconstraint")!=-1) objFactory = 1;
						else if (addelementtype.Find("TAEinput_output")!=-1) objFactory = 1;

						if (objFactory && (addelementtypename != mystr("")))
						//if (addelementtypename != mystr(""))
						{
							// AddObjectInfo
							add_class_descriptions += mystr("  AddObjectInfo(OFCElement,") + mystr('"') + addelementtypename + mystr('"') + mystr(", ") + mystr(addelementtype) + mystr(", ") + mystr('"') + texdescriptionExample + mystr('"')+ mystr(", ") + mystr('"') + shortdescription + mystr('"') + mystr(");\n");
							// AddElement
							elementtypecnt++;
							add_parse_add_code += mystr("    case ") + mystr(elementtypecnt) + mystr(": \n    {\n      ") +
								actualclass + mystr(" elem(GetMBS());\n      rv = GetMBS()->AddElement(&elem);\n      break;\n    }\n");
						}

						// always create entries for tex-file
						texDescrName = addelementtypename;
						if(shortdescription != mystr(""))
						{
							shortdescription.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Short description}{") + shortdescription +mystr("}\n\n");
						}
						if(texdescriptionDOF != mystr(""))
						{
							texdescriptionDOF.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Degrees of freedom}{") + texdescriptionDOF +mystr("}\n\n");
						}
						if(texdescriptionNode != mystr(""))
						{
							texdescriptionNode.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Nodes}{") + texdescriptionNode +mystr("}\n\n");
						}
						if(texdescriptionGeometry != mystr(""))
						{
							texdescriptionGeometry.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Geometry}{") + texdescriptionGeometry +mystr("}\n\n");
						}
						if(texdescriptionEquations != mystr(""))
						{
							texdescriptionEquations.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Equations}{") + texdescriptionEquations +mystr("}\n\n");
						}
						if(texdescriptionLimitations != mystr(""))
						{
							texdescriptionLimitations.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Limitations}{") + texdescriptionLimitations +mystr("}\n\n");
						}
						if(modus.Length())
						{
							texDescr += mystr("\\texdescriptionsubsec{Description of the different modi}{\n");
							texDescr += mystr("		\\modustablebegin\n");
							for(int mi = 1; mi <= modus.Length(); mi++)
							{
								texDescr += mystr("				\\modusline") +mystr(modus(mi))  +mystr("\n");
							}
							texDescr += mystr("		\\modustableend\n}\n\n");
						}
						if(texdescriptionComment != mystr(""))
						{
							texdescriptionComment.Replace("\\n","\n");
							texDescr += mystr("\\texdescriptionsubsec{Additional notes}{") + texdescriptionComment +mystr("}\n\n");
						}
						if(texdescriptionExample != mystr(""))
						{
							texDescr += mystr("% \\texdescriptionsubsec{Example}{\n");
							//texDescr += mystr("% \\verbatimtabinput[2]{D:/cpp/HotInt\_V1/documentation/EDCauto\_documentation/docu\_examples/") + mystr(texdescriptionExample) + mystr("}\n");
							texDescr += mystr("% \\verbatimtabinput[2]{EDCauto\_documentation/docu\_examples/") + mystr(texdescriptionExample) + mystr("}\n");	//$ DR 2013-09-20 relative path
							texDescr += mystr("% }\n\n");
						}

						for(int fi = 1; fi <= figure.Length(); fi++)
						{
							//texDescr += mystr("\\autodocufigure{D:/cpp/HotInt\_V1/documentation/EDCauto\_documentation/figures/") +mystr(figure(fi)) +mystr("}{") + mystr(caption(fi)) +mystr("}{") +mystr(addelementtypename) +mystr("figure") + mystr(fi) +mystr("}\n");
							texDescr += mystr("\\autodocufigure{./figures/") +mystr(figure(fi)) +mystr("}{") + mystr(caption(fi)) +mystr("}{") +mystr(addelementtypename) +mystr("figure") + mystr(fi) +mystr("}\n"); //$ DR 2013-09-20 relative path
						}
						// add labels for all references of this object
						//if(ref_count) {cout << " number of references = " << ref_count <<"\n";}
						for(int ri = 1; ri <= ref_count; ri++)
						{
							mystr rl = mystr(addelementtypename) +mystr("reference") + mystr(ri);
							ref_label.Add(rl);
							//cout << "ri = " << ri << " reference = " << ref_label.Last() << "\n"; 
						}
					}
					//fpos = commandlist.Find(mystr("endclass"));
					else if (c_endclass)
					{
						if (commandlist.Length() >= 2)
						{
							mystr str = commandlist(2);
							if (actualclass != str)
							{
								cout << "ERROR: class names of beginclass='" << actualclass << "' and endclass='" << str << "' do not match in line " << linecnt << "\n";
							}
							else
							{
								cout << "end class '" << actualclass << "'\n";
								getEDstr += "}\n\n";
								setEDstr += "  return 1;\n}\n\n";
								readSingle += "  return 0;\n}\n\n";
								writeSingle += "  return 0;\n}\n\n";
								availSingle += "  return 0;\n}\n\n";

								if (showdetails) cout << "    setEDstrRemove = '" << setEDstrRemove << "'\n";
								if(!setEDstrRemove.Compare(""))
								{
									setEDstrRemove = mystr("\n  ElementData ed;\n") + setEDstrRemove + mystr("\n");
									int i = setEDstr.Find("{");
									setEDstr.InsertAt(i+1,setEDstrRemove);
								}

								//finally write to HOTINT include file, separately for each class:
								outfile << "\n\n";
								outfile << "//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
								outfile << "//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
								outfile << "//automatically generated EDC-access\n";
								outfile << "//actual class=" << actualclass << "\n\n";
								outfile << getEDstr << "\n";
								outfile << setEDstr << "\n";
								outfile << readSingle << "\n";
								outfile << writeSingle << "\n";
								outfile << availSingle << "\n";

								getEDstr = "";
								setEDstr = "";
								readSingle = "";
								writeSingle = "";
								availSingle = "";
								setEDstrRemove = "";

								if(texDescrName!=mystr(""))
								{
									// write tex-file, separately for each class:
									mystr texfilename = mystr("../../../documentation/EDCauto\_documentation/classDescriptions/") + texDescrName + mystr(".tex");
									ofstream texfile(texfilename.c_str());
									//cout << "Writing to file " << texfilename.c_str() << "\n";
									texfile << "% This is an autogenerated file for the documetation of HOTINT\n";
									texfile << "% This file is only a small part of the documentation of the class " << texDescrName << "\n\n";
									texfile << texDescr << "\n\n";
									texfile << "% --------------------------------------------------------------------------------\n";
									texfile << "% end of autogenerated file\n";
									texfile.close();
								}
								texDescr = "";
								texDescrName = "";
							}

							actualclass = "";
						}
						else
						{
							cout << "ERROR: no class name found after $EDC$[endclass, ...] in line " << linecnt << "\n";
						}
					} else if (c_varacc || c_funcacc)
					{
						int pos = 0; //read line from beginning;
						//find type and name of variable:
						mystr add_type=""; //if function is virtual or variable is mutable
						mystr vartype; //variable type in C++
						mystr varname; //variable name in C++
						mystr EDCvarname; //variable name in EDC; without folder
						mystr EDCfolder=""; //folder name for EDC, e.g. "Geometry" for "Geometry.size"

						mystr funcname; //function name in C++
						int func_const = 0; //1==function is const
						int func_npar = 0; //number of parameters in function
						ParameterIDstruct par_id; //parameter ID structure of first function parameter

						int internalvar = 0; //internal access on this variable for (special value) sensors and/or modifiers
						int readonly = 0; //only GetElementData functions are written, data is not written
						int int_bool = 0; //use integer value as bool (only 0 / 1 possible) ==> as a check box
						int vecvarlen = 0; //vector has variable length (for input fields)
						int remove = 0;	//$ DR 2012-08-22: remove variable from an EDC of the parent class in a derived class. The elementdata is not really removed but locked and given a new tooltiptext

						//define possible minimum and maximum values for int/double values
						int usemaxval = 0;
						int useminval = 0;
						double maxval = 0;
						double minval = 0;

						//use subvector:
						mystr vecstart = "";
						mystr vecend = "";
						mystr func_par1 = ""; //function argument e.g. for Get(1)
						mystr condition = "";
						//int vecstart = 0;
						//int vecend = 0;
						//double func_par1 = 0; //function argument e.g. for Get(1)

						mystr tooltiptext = "";

						line.ReadLeadingSpaces(pos);
						vartype = line.GetWord(pos);
						if (showdetails) cout << "vartype=" << vartype << "\n";
						if(vartype == mystr("virtual")) //virtual is not interesting for parser (at least at the moment)
						{
							vartype = line.GetWord(pos); //now type should be there
							if (showdetails) cout << "vartype0=" << vartype << "\n";
						}
						if(vartype == mystr("mutable") || vartype == mystr("static") || vartype == mystr("const") || vartype == mystr("//EDC"))
						{
							add_type = vartype;
							vartype = line.GetWord(pos); //now type should be there
							if (showdetails) cout << "vartype=" << vartype << "\n";
						}
						int pos_var = pos; //position where variable starts
						varname = line.GetWord(pos);
						EDCvarname = varname; //this is the default, but it can be overwritten by the command "EDCvarname=..."
						mystr func_parameters = "";

						if (c_funcacc) 
						{
							//readonly = 1; //function is always read only!

							int l = varname.Length();
							if (varname[l-1] == ')' && varname[l-2] == '(')
							{
								EDCvarname = varname.SubString(0,l-3);
								funcname = EDCvarname;
								if (showdetails) cout << "func='" << EDCvarname << "'\n";
							}
							else
							{
								//read function:
								int posb = varname.Find('(');
								funcname = varname.SubString(0,posb);
								EDCvarname = funcname;
								if (showdetails) cout << "func='" << funcname << "'\n";

								mystr func_parameters = line.GetStringInBrackets(pos_var,'(',')');
								if (showdetails) cout << "func_par='" << func_parameters << "'\n";

								TArrayDynamic<mystr> parameter_list(0);
								GetCommaSeparatedStrings(func_parameters, parameter_list, showdetails);

								if (parameter_list.Length() == 1)
								{
									IdentifyParameter(parameter_list(1), par_id, showdetails);
								}
								else if (parameter_list.Length() > 1)
								{
									cout << "ERROR: only functions one parameter/argument can be used for EDC access, e.g. 'virtual int Get()' or 'int Get(int i)' in line " << linecnt << "\n";
									lineok = 0;
								}
							}
						}

						if (showdetails) 
						{
							if (c_varacc) 
							{
								cout << "read variable with name='" << varname << "' and type='" << vartype << "'\n";
							}
							else 
							{
								if (add_type == mystr("")) 
								{
									cout << "read function with name='" << varname << "' and type='" << vartype << "'\n";
								}
								else 
								{
									cout << "read function with name='" << varname << "' and type='" << add_type << " " << vartype << "'\n";
								}
							}
						}
						if (vartype.Length() == 0)
						{
							lineok = 0;
							cout << "ERROR: variable type of variable '" << varname << "' not found in line " << linecnt << "\n";
						}
						if (varname.Length() == 0)
						{
							lineok = 0;
							cout << "ERROR: variable name not found in line " << linecnt << "\n";
						}

						if (lineok)
						{
							//read further commands and options for varaccess
							//varaccess must be first command
							for (int i=2; i <= commandlist.Length(); i++)
							{
								mystr str = commandlist(i);
								int pos = 0;
								mystr comword = str.GetWord(pos); //either read single word, or assignment
								//cout << "word='" << comword << "', pos = " << pos << ", strlength=" << str.Length() << "\n";
								if (pos < str.Length() && pos != -1) //==> read assignment; if pos is not at end of string and pos is not -1 (end of string reached), then this must be an assignment
								{
									if (str[pos] != '=')
									{
										cout << "ERROR: expected '=' after command, but received '" << str[pos] << "' in line " << linecnt << "\n";
									}
									pos ++;
									//cout << "pos=" << pos << ", str.Length()-1=" << str.Length()-1 << "\n";
									mystr rhs = str.SubString(pos, str.Length()-1);

									if (comword == mystr("EDCvarname"))
									{
										mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
										EDCvarname = rhsstr;
										if (showdetails) cout << "    EDCvarname = '" << EDCvarname << "'\n";
									}
									else if (comword == mystr("EDCfolder"))
									{
										//cout << "rhs.Length()-2=" << rhs.Length()-2 << "\n";
										if (rhs.Length()-2 >= 1)
										{
											mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
											EDCfolder = rhsstr;
											if (showdetails) cout << "    EDCfolder = '" << EDCfolder << "'\n";
										}
										//else the EDCfolder string is ""
									}
									else if (comword == mystr("tooltiptext"))
									{
										mystr rhsstr = rhs.SubString(1, rhs.Length()-2);
										tooltiptext = rhsstr;
										if (showdetails) cout << "    tooltiptext = '" << tooltiptext << "'\n";
									}
									else if (comword == mystr("maxval"))
									{
										maxval = rhs.MakeDouble();
										usemaxval = 1;
										if (showdetails) cout << "    maxval = '" << maxval << "'\n";
									}
									else if (comword == mystr("minval"))
									{
										minval = rhs.MakeDouble();
										useminval = 1;
										if (showdetails) cout << "    minval = '" << minval << "'\n";
									}
									else if (comword == mystr("vecstart"))
									{
										vecstart = rhs;//.MakeInt();
										if (showdetails) cout << "    vecstart = '" << vecstart << "'\n";
									}
									else if (comword == mystr("vecend"))
									{
										vecend = rhs;//.MakeInt();
										if (showdetails) cout << "    vecend = '" << vecend << "'\n";
									}
									else if (comword == mystr("func_par1")) //for single function parameter1
									{
										func_par1 = rhs;//.MakeInt();
										if (showdetails) cout << "    func_par1 = '" << func_par1 << "'\n";
									}
									else if (comword == mystr("condition")) //for single function parameter1
									{
										condition = rhs;
										if (showdetails) cout << "    condition = '" << condition << "'\n";
									}

								}
								else //evaluate option
								{
									if (comword == mystr("internal_varaccess")) 
									{
										internalvar = 1;
										if (showdetails) cout << "    internal\n";
									}
									if (comword == mystr("readonly")) 
									{
										readonly = 1;
										if (showdetails) cout << "    readonly\n";
									}
									if (comword == mystr("int_bool")) 
									{
										int_bool = 1;
										if (showdetails) cout << "    int_bool\n";
									}
									if (comword == mystr("variable_length_vector"))
									{
										vecvarlen = 1;
										if (showdetails) cout << "    variable_length_vector\n";
									}
									if (comword == mystr("remove")) 
									{
										remove = 1;
										readonly = 1;
										if (showdetails) cout << "    remove\n";
									}

								}
							}

							if (c_funcacc)
							{
								//add "()", if there are no arguments
								if (func_par1 == mystr("") && (vartype == mystr("int") || vartype == mystr("double") || vartype == mystr("mystr") || vartype == mystr("char*"))) 
								{
									varname = funcname + mystr("()");
									//cout << "func-varname1=" << varname << "\n";
								}
								else if (func_par1 != mystr("") && (vartype == mystr("int") || vartype == mystr("double"))) 
								{
									varname = funcname + mystr("(") + func_par1 + mystr(")");
									//cout << "func-varname2=" << varname << "\n";
								}
								else if (vartype == mystr("Vector") || vartype == mystr("IVector") || vartype == mystr("TArray<int>") || vartype == mystr("Vector3D") || vartype == mystr("Vector2D") &&
									vecstart != mystr(""))
								{
									varname = funcname;
									//cout << "func-varname3=" << varname << "\n";
									//this is automatically handled by Vector types
								}
								else
								{
									cout << "ERROR: funcaccess not supported for this types in line " << linecnt << "\n";
								}
							}

							//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
							//now write variable access into Set/GetElementData functions:

							int required_value = 1; //this flag can be user-defined lateron; it indicates, that a variable is mandatory
							mystr required_value_str = mystr(required_value);
							mystr EDCtotalname;

							if (EDCfolder == mystr(""))
							{
								EDCtotalname = EDCvarname;
							}
							else
							{
								EDCtotalname = EDCfolder + mystr(".") + EDCvarname;
							}


							mystr readonly_str = "";
							if (readonly) {readonly_str = "ed.SetLocked(1); ";}

							mystr tooltipstr;
							if (tooltiptext != mystr("")) {tooltipstr = mystr("ed.SetToolTipText(")+mystr('"') + tooltiptext + mystr('"') + mystr("); ");}

							if (vecvarlen) {tooltipstr = mystr("ed.SetVariableLength(); ")+tooltipstr; }

							//generate string which always ends the getEDstr command; allowing for adjustable options
							mystr treeadd_str_noendl = readonly_str + tooltipstr + mystr("edc.TreeAdd(")+mystr('"')+EDCfolder+mystr('"')+mystr(",ed);");
							mystr treeadd_str = treeadd_str_noendl + mystr("\n");

							mystr getEDstr_inc = "";
							mystr setEDstr_inc = "";
							mystr readsingle_inc = "";
							mystr writesingle_inc = "";
							mystr availsingle_inc = "";
							mystr setEDstrRemove_inc ="";

              // cascade for Read&WriteSingleElementData - distinguish only 0,1,2 indices
							if (internalvar)
							{
                availsingle_inc += mystr("  available_variables.Add(ReadWriteElementDataVariableType(") + mystr('"') + EDCtotalname + mystr('"') + mystr(", ");
	
								if ((vartype == mystr("int")) || (vartype == mystr("double")))
								{
									//common "head"
		//							mystr common = mystr("  if (RWdata.variable_name == mystr(") +  mystr('"') + EDCtotalname + mystr('"') + mystr(")) \n  {\n"); // too slow
									mystr common = mystr("  if (!strcmp(RWdata.variable_name.c_str(),") +  mystr('"') + EDCtotalname + mystr('"') + mystr(")) \n  {\n"); // faster
									// no index, thus no index.check...
									readsingle_inc += common; writesingle_inc += common; 
									availsingle_inc += mystr("0, 0, 0., ");
									
									readsingle_inc += mystr("    RWdata.value = ") + varname;
									writesingle_inc += mystr("    ") + varname + mystr(" = RWdata.value");

									//common "tail"
									common = mystr("; return 1; \n  }\n");
									readsingle_inc += common; writesingle_inc += common; 
								}
								else if ((vartype == mystr("Vector2D")) || (vartype == mystr("Vector3D")) || (vartype == mystr("Vector")) || (vartype == mystr("IVector")))
								{
									//common "head"
									mystr common = mystr("  if (!strcmp(RWdata.variable_name.c_str(),") +  mystr('"') + EDCtotalname + mystr('"') + mystr(")) \n  {\n");
									common += mystr("    if ((RWdata.comp1 < 1) "); 	
									if (vartype == mystr("Vector2D")) common += mystr("|| (RWdata.comp1 > 2)");
									if (vartype == mystr("Vector3D")) common += mystr("|| (RWdata.comp1 > 3)");
									if (vartype == mystr("Vector"))   common += mystr("|| (RWdata.comp1 > ") + varname + mystr(".Length())");
									if (vartype == mystr("IVector"))  common += mystr("|| (RWdata.comp1 > ") + varname + mystr(".Length())");
									common += mystr(") return -2; \n");
									readsingle_inc += common; writesingle_inc += common; 
								  availsingle_inc += varname + mystr(".Length(), 0, 0., ");

									readsingle_inc += mystr("    RWdata.value = ") + varname + mystr("(RWdata.comp1)");
									writesingle_inc += mystr("    ") + varname + mystr("(RWdata.comp1) = RWdata.value");
								
									//common "tail"							
									common = mystr("; return 1;\n  }\n");
									readsingle_inc += common; writesingle_inc += common; 

								}
								else if ((vartype == mystr("Matrix3D")) || (vartype == mystr("Matrix")))
								{
									//common "head"
									mystr common = mystr("  if (!strcmp(RWdata.variable_name.c_str(),") +  mystr('"') + EDCtotalname + mystr('"') + mystr(")) \n  {\n");
								 	common += mystr("    if ((RWdata.comp1 < 1) || (RWdata.comp2 < 1) "); 	
									if(vartype == mystr("Matrix3D")) common += mystr("|| (RWdata.comp1 > 3) || (RWdata.comp2 > 3)");
									if(vartype == mystr("Matrix"))   common += mystr("|| (RWdata.comp1 > ") + varname + mystr(".Getrows()) || (RWdata.comp2 > ") + varname + mystr(".Getcols())");
									common += mystr(") return -2; \n");
									readsingle_inc += common; writesingle_inc += common; 
								  availsingle_inc += varname + mystr(".Getrows(), ") + varname + mystr(".Getcols(), 0., ");

									readsingle_inc += mystr("    RWdata.value = ") + varname + mystr("(RWdata.comp1,RWdata.comp2)");
									writesingle_inc += mystr("    ") + varname +  mystr("(RWdata.comp1,RWdata.comp2) = RWdata.value");

									//common "tail"							
									common = mystr("; return 1; \n  }\n");
									readsingle_inc += common; writesingle_inc += common; 
								}
								else
								{
								  // n.i.y.   e.g. mystr, char*
								}
								availsingle_inc += mystr('"') + tooltiptext + mystr('"') + mystr(")); \n");
							}

							// cascade for Get&SetElementData
							if (vartype == mystr("int") && int_bool)
							{
								getEDstr_inc += mystr("  ed.SetBool(") + varname + mystr(",")+mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str;
								setEDstr_inc += mystr("  GetElemDataBool(GetMBS(), edc, ") + mystr('"') + EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  int dummy=0; ed.SetBool(dummy,")+mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("int"))
							{
								getEDstr_inc += mystr("  ed.SetInt(") + varname + mystr(",")+mystr('"')+EDCvarname+mystr('"');
								if(useminval)
								{
									if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
									else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
								}
								else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
								getEDstr_inc += mystr("); ") + treeadd_str;

								setEDstr_inc += mystr("  GetElemDataInt(GetMBS(), edc, ") + mystr('"') + EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  int dummy=0; ed.SetInt(dummy,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("double"))
							{
								getEDstr_inc += mystr("  ed.SetDouble(") + varname + mystr(",")+mystr('"')+EDCvarname+mystr('"');
								if(useminval)
								{
									if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
									else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
								}
								else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
								getEDstr_inc += mystr("); ") + treeadd_str;

								setEDstr_inc += mystr("  GetElemDataDouble(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  double dummy= 0.; ed.SetDouble(dummy,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("Vector2D"))
							{
								if (vecstart != mystr(""))
								{
									mystr vecstr = varname + mystr("(") + mystr(vecstart) + mystr("),")
										+ varname + mystr("(") + mystr(vecstart) + mystr("+1)");
									getEDstr_inc += mystr("  ed.SetVector2D(") + vecstr + mystr("),") + mystr('"')+EDCvarname+mystr('"');
									if(useminval)
									{
										if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
										else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
									}
									else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
									getEDstr_inc += mystr("); ") + treeadd_str;

									setEDstr_inc += mystr("  {Vector2D vv; ");
									setEDstr_inc += mystr(" GetElemDataVector2D(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",vv, " + required_value_str+ ");");
									setEDstr_inc += mystr(" vv.Get(") + vecstr + mystr(");}\n");
								}
								else
								{
									getEDstr_inc += mystr("  ed.SetVector2D(") + varname + mystr(".X(),") + varname + mystr(".Y(),") + mystr('"')+EDCvarname+mystr('"');
									if(useminval)
									{
										if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
										else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
									}
									else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
									getEDstr_inc += mystr("); ") + treeadd_str;
									setEDstr_inc += mystr("  GetElemDataVector2D(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								}
								if (remove) {setEDstrRemove_inc += mystr("  {  Vector2D dummy(0.); ed.SetVector2D(dummy.X(),dummy.Y(),")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("Vector3D"))
							{
								if (vecstart != mystr(""))
								{
									mystr vecstr = varname + mystr("(") + mystr(vecstart) + mystr("),")
										+ varname + mystr("(") + mystr(vecstart) + mystr("+1),")
										+ varname + mystr("(") + mystr(vecstart) + mystr("+2)");
									getEDstr_inc += mystr("  ed.SetVector3D(") + vecstr + mystr(",") + mystr('"')+EDCvarname+mystr('"');
									if(useminval)
									{
										if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
										else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
									}
									else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
									getEDstr_inc += mystr("); ") + treeadd_str;

									setEDstr_inc += mystr("  {Vector3D vv; ");
									setEDstr_inc += mystr(" GetElemDataVector3D(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",vv, " + required_value_str+ ");");
									setEDstr_inc += mystr(" vv.Get(") + vecstr + mystr(");}\n");
								}
								else
								{
									getEDstr_inc += mystr("  ed.SetVector3D(") + varname + mystr(".X(),") + varname + mystr(".Y(),") + varname + mystr(".Z(),") + mystr('"')+EDCvarname+mystr('"');
									if(useminval)
									{
										if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
										else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
									}
									else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
									getEDstr_inc += mystr("); ") + treeadd_str;
									setEDstr_inc += mystr("  GetElemDataVector3D(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								}
								if (remove) {setEDstrRemove_inc += mystr("  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("Vector"))
							{
								if (vecstart != mystr("") && vecend != mystr(""))
								{
									mystr vecstr = mystr("Vector vv((") + mystr(vecend) + mystr(")-(") + mystr(vecstart) + mystr(")+1);");
									getEDstr_inc += mystr("\n  {") + vecstr + mystr(" for (int i=") + 
										mystr(vecstart) + mystr("; i<=") + mystr(vecend) + mystr("; i++) {vv(i+1-(") + mystr(vecstart) + mystr("))=")+varname+mystr("(i);}");
									getEDstr_inc += mystr("\n  ed.SetVector(vv.GetVecPtr(),vv.Length(),") + mystr('"')+EDCvarname+mystr('"');
									if(useminval)
									{
										if(usemaxval)			{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval;}
										else							{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",1");}
									}
									else if(usemaxval)	{ getEDstr_inc += mystr(",") + minval + mystr(",") + maxval + mystr(",0,1");}
									getEDstr_inc += mystr("); ") + treeadd_str_noendl + mystr("}\n");

									setEDstr_inc += mystr("  {") + vecstr;
									setEDstr_inc += mystr("  GetElemDataVector(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",vv, " + required_value_str+ ");\n");
									setEDstr_inc += mystr("  for (int i=") + mystr(vecstart) + mystr("; i<=") + mystr(vecend) + mystr("; i++) {")+varname+mystr("(i)=vv(i+1-(") + mystr(vecstart) + mystr("));}\n}\n");
								}
								else
								{
									getEDstr_inc += mystr("  ed.SetVector(") + varname + mystr(".GetVecPtr(),") + varname + mystr(".Length(),") + mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str;
									setEDstr_inc += mystr("  GetElemDataVector(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								}
								if (remove) {setEDstrRemove_inc += mystr("  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("IVector") || vartype == mystr("TArray<int>"))
							{
								if (vecstart != mystr("") || vecend != mystr(""))
								{
									cout << "ERROR: vecstart or vecend not implemented for IVector. Problem in line " << linecnt << "\n";
								}
								getEDstr_inc += mystr("\n  {Vector vv(") + varname + mystr(".Length()); for(int i = 1; i <= ") + varname + mystr(".Length(); i++) {vv(i) = ") + varname + mystr("(i);}\n") ;
								getEDstr_inc += mystr("  ed.SetVector(vv.GetVecPtr(),vv.Length(),") + mystr('"')+EDCvarname+mystr('"')+mystr("); ed.SetValuesInt(); ") + treeadd_str + mystr("}\n");

								setEDstr_inc += mystr("  GetElemDataIVector(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("Matrix3D"))
							{
								getEDstr_inc += mystr("\n  Matrix ") + varname + mystr("_1(") + varname + mystr(");\n");
								getEDstr_inc += mystr("  ed.SetMatrix(") + varname + mystr("_1.GetMatPtr(),") + varname + mystr("_1.Getrows(),") + varname + mystr("_1.Getcols(),") + mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str + mystr("\n");

								setEDstr_inc += mystr("\n  Matrix ") + varname + mystr("_1;\n");
								setEDstr_inc += mystr("  GetElemDataMatrix(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr("_1, " + required_value_str+ ");\n");
								setEDstr_inc += mystr("  ") + varname + mystr(" = Matrix3D(") + varname + mystr("_1);\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("Matrix"))
							{
								getEDstr_inc += mystr("\n  Matrix ") + varname + mystr("_1(") + varname + mystr(");\n");
								getEDstr_inc += mystr("  ed.SetMatrix(") + varname + mystr("_1.GetMatPtr(),") + varname + mystr("_1.Getrows(),") + varname + mystr("_1.Getcols(),") + mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str;

								setEDstr_inc += mystr("  GetElemDataMatrix(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  Matrix dummy(1,1); ed.SetMatrix(dummy.GetMatPtr(),1,1,")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else if (vartype == mystr("mystr") || vartype == mystr("char*"))
							{
								mystr cstr = mystr(".c_str()");
								if (vartype == mystr("char*")) cstr = "";

								getEDstr_inc += mystr("  ed.SetText(") + varname + cstr+mystr(',') + mystr('"')+EDCvarname+mystr('"')+mystr("); ") + treeadd_str;

								setEDstr_inc += mystr("  GetElemDataText(GetMBS(), edc, ") + mystr('"')+EDCtotalname+mystr('"') + mystr(",") + varname + mystr(", " + required_value_str+ ");\n");
								if (remove) {setEDstrRemove_inc += mystr("  {  ed.SetText(")+ mystr('"')+mystr("dummy")+mystr('"')+mystr(",")+mystr('"')+EDCvarname+mystr('"') + mystr("); ") + treeadd_str_noendl + mystr("}\n");}
							}
							else
							{
								cout << "ERROR: did not understand variable type '" << vartype << "' in line " << linecnt << "\n";
							}

							if (condition != mystr(""))
							{
								getEDstr_inc = mystr("  if(") + condition + mystr(")\n  {") + getEDstr_inc + mystr("  }\n");
								setEDstr_inc = mystr("  if(") + condition + mystr(")\n  {") + setEDstr_inc + mystr("  }\n");
							}

							if (readonly) 
							{
								setEDstr_inc = "" ; //reset to string before adding readonly variable (no set function needed/possible (e.g. in case of function access)
								writesingle_inc = "";
							}

							if (remove) //$ DR 2012-12-11
							{
								//getEDstr_inc = mystr("  {ElementData* edp = 	edc.TreeFind(") +mystr('"') +EDCtotalname +mystr('"') + mystr("); edp->SetHidden(1); edp->SetToolTipText(")+mystr('"')+mystr("data from parent class, not available in this element!")+mystr('"')+mystr(");	}\n");
								getEDstr_inc = mystr("  edc.TreeDelete(")	+mystr('"') +EDCtotalname +mystr('"') + mystr("); \n");
							}

							getEDstr += getEDstr_inc;
							setEDstr += setEDstr_inc;
							readSingle += readsingle_inc;
							writeSingle += writesingle_inc;
							availSingle += availsingle_inc;
							setEDstrRemove += setEDstrRemove_inc;

						}

					}
					else
					{
						cout << "ERROR: no command (funcaccess/varaccess) found " << "' in line " << linecnt << "!\n";
					}


				}
				//cout << "\n";
			}


			//end parse one line
			//+++++++++++++++++++++++++++++++++++

			if (file.EndOfFile() || !file.IsGood()) endfile = 1;
		}

		if (linecnt >= maxlines)
		{
			cout << "++++++++++++++++++++++++++++\n";
			cout << "++++++++++++++++++++++++++++\n";
			cout << "ERROR: maximum number of lines per file (100000) reached, parsing interupted!!!\n";
			cout << "++++++++++++++++++++++++++++\n";
			cout << "++++++++++++++++++++++++++++\n";
		}


	} //end go through all files

	//C:\Dokumente und Einstellungen\gerstm\Eigene Dateien\cpp\HotIntCVS\tools\data_reader\release>
	//data_reader.exe -filelist ..\..\EDC_converter\testfiles\filelist.txt -destination ..\..\EDC_converter\testfiles\testEDC.cpp

	//outfile << "void CEDCParser::AddClassDescriptions_Auto()\n{\n";
	//outfile << add_class_descriptions << "\n";
	//outfile << "}\n\n";

	outfile << "void MBSObjectFactory::AddObjectInfos_Auto()\n{\n";
	outfile << add_class_descriptions << "\n";
	outfile << "}\n\n";

	outfile << "int MBSObjectFactory::AddElement(int element_type_id)\n{\n";
	outfile << "  int rv = -1;\n";
	outfile << "  switch(element_type_id)\n  {\n";
	outfile << add_parse_add_code << "\n";
	outfile << "  default: ; \n  }\n";
	outfile << "  return rv;\n";
	outfile << "}\n";
	outfile.close();

	// write bibliography
	ofstream outfilebib(bibfile.c_str());
	outfilebib << "%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfilebib << "%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfilebib << "%+    automatically generated file for auto documentation, do not modify!!!    +\n";
	outfilebib << "%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	outfilebib << "%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";

	for(int i=1; i<=ref.Length(); i++)
	{
		outfilebib << "\\bibitem\{" << ref_label(i) <<"\} \n";
		outfilebib << "			" << ref(i) << "\n \n";
	}

	outfilebib.close();

	return 1;
}


