//#**************************************************************
//#
//# filename:             Sim2Hotint.cpp
//#
//# authors:              Rafael Ludwig
//#												Michael Stangl
//#
//# generated:						January 2011
//# description:          Import functions for Simulink Models into HOTINT
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
//#**************************************************************

#include "element.h"
#include "body2d.h"
#include "body3d.h"
#include "node.h"
#include "..\\ElementsLib\\deprecated_elements\\sensor_old.cpp"
#include "constraint.h"
#include "control.h"
#include "sim2hotint.h"
#include "myfile.h"
#include "elementdataaccess.h"
#include "options_class_auto.h"

//typedef int (*ModelFunctionAutoRegistration_pt2FunctionList)(MBS * mbs);//[modelFunctionAutoRegistration_max_n_fct];

//int(*function_ptr)(MBS * mbs)


#include "Sim2Hotint.h"
// add subsystem from simulink with inports and outports
void AddSimulinkSubSystem(MBS * mbs,mystr fname,TArray<int>& elnums_inputs,TArray<int>& elnums_outputs, Vector2D drawoffset, double scaling_factor)
{
	Simulink2HotintConversion S2Hc(mbs);
	S2Hc.SetDrawOffset(drawoffset);
	S2Hc.SetScaleFactor(scaling_factor);
	
	if(S2Hc.LoadSimulink2Hotint(fname))
	{
		// file exists!
		for(int i=1; i<=S2Hc.GetNumberOfInports(); i++ )
		{
			elnums_inputs.Add(S2Hc.GetInportMBSElementNumber(i));
			mbs->UO() << "inport" << i << " --> " << "element_number = " << elnums_inputs.Last() << "\n";
		}
		
		for(int i=1; i<=S2Hc.GetNumberOfOutports(); i++ )
		{
			elnums_outputs.Add(S2Hc.GetOutportMBSElementNumber(i));
			mbs->UO() << "outport" << i << " --> " << "element_number = " << elnums_outputs.Last() << "\n";
		}
	}
}

int Simulink2HotintConversion::LoadSimulink2Hotint(const mystr& fname)
{

	ElementDataContainer edc;
	//ElementData ed;
	if(ReadSimulinkModel(fname,&edc))
	{		
		mbs->UO(UO_LVL_0) << "read " << mystr(fname) << "\n";		
		edc_default = (edc.TreeFind("Model.BlockParameterDefaults"))->GetEDC();
		//ParseSimulinkEDCModel_Default(edc_default);

		
		//for(int i=1; i<=block_parameter_default_EDC_list.Length(); i++)
		//{
		//	mystr str;
		//	EDC2String(*block_parameter_default_EDC_list.Get(i), str);
		//	mbs->UO() << str << "\n";
		//}
	
		// find "system"
		ElementDataContainer* edc_system = (edc.TreeFind("Model.System"))->GetEDC();
		
		ParseSimulinkEDCModel(edc_system);
		return 1;
		//ConvertBlock2Element(&edc);
	}
	else
	{
		mbs->UO(UO_LVL_warn) << "can't read " << mystr(fname) << "\n";
		return 0;
	}
}

//ALTERNATIVE VERSION FOR PARSING MODELS IN SIMULINK DATA FORMAT (mdl-files)
// new version for parsing a string into an ElementDataContainer inclusive ToolTipText
int Simulink2HotintConversion::String2TreeEDC_SIM(mystr& str, ElementDataContainer& edc)
{
	int rv = 1; //success

	mystr el = "\n"; //end of line
	char elc = 'n';
	char obc = '{';  //open brace, for sub structure
	char cbc = '}';  //closing brace
	mystr ob(obc);
	mystr cb(cbc);
	char osbc = '[';  //open square bracket, for vector or matrix
	char csbc = ']';  //closing square bracket
	mystr osb = osbc;  //open square bracket, for vector or matrix
	mystr csb = csbc;  //closing square bracket
	mystr ds = ": ";	//double stop - not used any more
	mystr eq = "=";		//equal sign
	mystr semicolon = ";"; //double stop
	char semicolonc = ';';
	char commac = ','; //comma
	mystr comma = commac; //comma
	char doublequote = '"'; //for strings or text

	char commentc = '%'; //comment
	mystr comment = commentc; //comment
	mystr matrixsep = el; //may be changed to semicolon ...

	int oo = 0;
	int endstr = 0;
	int pos = 0;
	//int posold;
	mystr dataname;
	char ch;
	ElementData ed;
	mystr dummycomment; //dummy comment string

	//int donotuse_sign = (option & 1) == 1;

	while (!endstr)
	{
		ch = str.ReadLeadingSpaces(pos);
		while (ch == commentc && pos != -1)
		{
			str.GetUntilEOL(pos, (char)0, dummycomment);
			ch = str.ReadLeadingSpaces(pos);
		}


		if (pos == -1) //regular end of structure!
		{
			endstr = 1;
		}
		else
		{
			dataname = str.GetWord(pos, 0); //name of data structure or data variable
		}

		if (pos == -1) //regular end of structure!
		{
			endstr = 1;
		}
		else
		{
			if (oo) mbs->UO() << "Dataname='" << dataname << "'\n";

			ch = str.ReadLeadingSpaces(pos);

			while (ch == commentc && pos != -1)
			{
				str.GetUntilEOL(pos, (char)0, dummycomment);
				ch = str.ReadLeadingSpaces(pos);
			}

			if (ch == ob)
			{
				//read sub item
				//pos--;
				mystr substr = str.GetStringInBrackets(pos, obc, cbc);
				if (pos != -1)
				{
					ElementDataContainer edc_sub;
					String2TreeEDC_SIM(/*mbs, */substr, edc_sub);
					ed.SetEDC(&edc_sub, dataname.c_str());
					edc.Add(ed);
				}
				else
				{
					mbs->UO() << "ERROR:String2TreeEDC_SIM: could not read item '" << dataname << "' !\n";
					return 0; //end of file reached, not valid!
				}
			}			
			else 
			{
				if (ch == eq)
				{
					//read data:
					pos++; //go after eq sign
					ch = str.ReadLeadingSpaces(pos);
				}
				if (ch == osb) //vector/matrix
				{
					mystr text = str.GetStringInBrackets(pos, osbc, csbc);					

					if(text.Find(el)==-1 && text.Find(semicolon)==-1)
					{
						// vector has no semicolon and is written in one single line
						TArray<double> ta_vec;
						mystr s_value;						
						int posvec = 0;
						int possval = 0;
						int comma_found = 0;
						while(posvec != -1)
						{
							comma_found = text.GetUntil(posvec, commac, s_value); // 'value,'
							s_value.ReadLeadingSpaces(possval);
							ta_vec.Add(s_value.MakeDouble());

							if(comma_found)posvec++;
						}									
						Vector vec(ta_vec);
						//if(vec.Length() == 3)ed.SetVector3D(ta_vec(1), ta_vec(2), ta_vec(3), dataname); //3D
						ed.SetVector(vec.GetVecPtr(), vec.Length(), dataname);

						//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
						// comment --> ToolTipText
						mystr data;
						int rv = str.GetUntilEOL(pos, commentc, data);
						if (rv == 2) //comment found
						{
							if (str.PosPeek(pos) != commentc) mbs->UO() << "problem with comment\n";

							// comment --> ToolTipText (double)
							str.GetUntilEOL(pos, (char)0, dummycomment);
							mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
							ed.SetToolTipText(comment);
						}
						//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

						edc.Add(ed);
					}
					else
					{
						//no vector --> matrix (has semicolon)
						// vector has no semicolon and is written in one single line
						TArray<double> ta_vec;
						mystr s_value;						
						int posvec = 0;
						int possval = 0;
						int found = 0; //',' or ';' found
						int nCols = 0;
						while(posvec != -1)
						{
							int posvec_old = posvec;
							mystr text_old = text;
							found = text.GetUntil(posvec, commac, s_value); // 'value,'
							if(s_value.Find(semicolon)!=-1 && !nCols)
							{
								nCols = ta_vec.Length()+1;
							}
							if(s_value.Find(semicolon)!=-1)
							{
								text = text_old;
								posvec = posvec_old;
								found = text.GetUntil(posvec, semicolonc, s_value); // 'value;'
							}
							s_value.ReadLeadingSpaces(possval);
							ta_vec.Add(s_value.MakeDouble());

							if(found)posvec++;
						}	
						Vector vec(ta_vec);
						ed.SetMatrix(vec.GetVecPtr(), (int)(vec.Length()/nCols), nCols, dataname);
					}

					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// comment --> ToolTipText
					mystr data;
					int rv = str.GetUntilEOL(pos, commentc, data);
					if (rv == 2) //comment found
					{
						if (str.PosPeek(pos) != commentc) mbs->UO() << "problem with comment\n";

						// comment --> ToolTipText (double)
						str.GetUntilEOL(pos, (char)0, dummycomment);
						mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
						ed.SetToolTipText(comment);
					}
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

					edc.Add(ed);
				}
				else if (ch == doublequote)
				{
					mystr text = str.GetStringInBrackets(pos, doublequote, doublequote);
					ed.SetText(text, dataname);
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// comment --> ToolTipText
					mystr data;
					int rv = str.GetUntilEOL(pos, commentc, data);
					if (rv == 2) //comment found
					{
						if (str.PosPeek(pos) != commentc) mbs->UO() << "problem with comment\n";

						// comment --> ToolTipText (double)
						str.GetUntilEOL(pos, (char)0, dummycomment);
						mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
						ed.SetToolTipText(comment);
					}
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					edc.Add(ed);
				}
				else if (!IsNum(ch))
				{
					mystr text = str.GetWord(pos, 0);
					//mystr text = str.GetStringInBrackets(pos, doublequote, doublequote);
					ed.SetText(text, dataname);
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// comment --> ToolTipText
					mystr data;
					int rv = str.GetUntilEOL(pos, commentc, data);
					if (rv == 2) //comment found
					{
						if (str.PosPeek(pos) != commentc) mbs->UO() << "problem with comment\n";

						// comment --> ToolTipText (double)
						str.GetUntilEOL(pos, (char)0, dummycomment);
						mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
						ed.SetToolTipText(comment);
					}
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					edc.Add(ed);
				}
				else //must be double or integer number
				{
					mystr data;
					int rv = str.GetUntilEOL(pos, commentc, data);

					data.EraseSpaces();
					if (oo) mbs->UO() << "single data='" << data << "'\n";

					ed.SetDouble(data.MakeDouble(), dataname); 

					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// comment --> ToolTipText (double)
					if (rv == 2) //comment found
					{
						if (str.PosPeek(pos) != commentc) mbs->UO() << "problem with comment\n";

						str.GetUntilEOL(pos, (char)0, dummycomment);
						mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
						ed.SetToolTipText(comment);
					}
					//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					edc.Add(ed);
				}
			}
		}
	}

	return rv;
}

int Simulink2HotintConversion::ReadSimulinkModel(const mystr& fname, ElementDataContainer* edc)
{
	CMFile file(fname, TFMread);
	if (file.IsGood()) 
	{
		//read parameters from file	
		//mbs->UO()	<< "Read parameters from " << fname << "\n";
		mystr str;
		file.RWF(str); //read file	
		//mbs->UO()	<< str;

		String2TreeEDC_SIM(str, *edc); // write file contents into edc (multiple entries possible!!!)

		if(!mbs->GetModelDataContainer())
		{
			mbs->SetModelDataContainer(*edc);
		}
		else
		{
			for(int i=1;i<=edc->Length();i++)
			{
				ElementDataContainer* edc_model = mbs->GetModelDataContainer();
				edc_model->Add(edc->Get(i));
			}
		}
		return 1; //file found
	}
	return 0; //file not found
}

//void Simulink2HotintConversion::ParseSimulinkEDCModel_Default(const ElementDataContainer* edc)
//{
//	for(int i=1; i<=edc->Length(); i++)
//	{
//		//type 1==int, 2==double, 4==vector, 8==text, 16==BOOL, 32==vector2D/3D,
//		;									 //    64==elementtype, 128==actionDialogFunctionWCdriver, 256=actionCompFunction
//		;									 //    512==checkgroup=idata2, 1024==matrix, 2048==vdata is integer, 4096==variable length vector, 
//		;									 //    8192==ElementDataContainer
//		//if(edc->Get(i).GetType() == 8192)
//		//{
//		//	ParseSimulinkEDCModel(edc.Get(i));
//		//}
//		mbs->UO() << edc->Get(i).GetDataName() << "\n";
//		if(edc->Get(i).IsEDC())
//		{
//			mystr name = edc->Get(i).GetDataName();
//			if(!strcmp(name.c_str(), "Block"))
//			{
//				mbs->UO() << "Default found!" << "\n";
//  			const ElementDataContainer* edc_block = (edc->Get(i)).GetEDC();
//				block_parameter_default_EDC_list.Add(edc_block);
//			}
//		}
//	}
//}

// edc ... system edc of mdl-file
void Simulink2HotintConversion::ParseSimulinkEDCModel(ElementDataContainer* edc)
{
	for(int i=1; i<=edc->Length(); i++)
	{
		mbs->UO() << edc->Get(i).GetDataName() << "\n";
		if(edc->Get(i).IsEDC())
		{
			mystr name = edc->Get(i).GetDataName();
			if(!strcmp(name.c_str(), "Block"))
			{
				mbs->UO() << "Block found!" << "\n";
				ElementDataContainer* edc_block = edc->Get(i).GetEDC()->GetCopy();
				edc_blocks.Add(edc_block);
				Element* elptr = ConvertBlock2Element(edc_blocks.Last()); // returns 0 in case of MBSSensor 				

				/*int inport_num = edc_blocks.Last()->TreeGetInt("MBS_inport_elnum");
				mbs->UO() << "inport_num = " << inport_num << "\n";*/
			}
			else if(!strcmp(name.c_str(), "Line"))
			{
				mbs->UO() << "Line found!" << "\n";	
				ElementDataContainer* edc_line = edc->Get(i).GetEDC()->GetCopy();
				edc_lines.Add(edc_line);
				ConnectBlocks(edc_lines.Last());
			}
		}
	}
}

int Simulink2HotintConversion::GetRotation(MBS * mbs, ElementDataContainer* edc_block, const char* colkey, int defaultrot) //returns the numbers of 90° rotations
{
	if(!edc_block->TreeFind(colkey))
		return defaultrot;

	const char* colstring = edc_block->TreeGetString(colkey, "not_found");
	if(!strcmp(colstring, "right"))	return 0;
	else if(!strcmp(colstring, "up")) return 1; 
	else if(!strcmp(colstring, "left")) return 2; 
	else if(!strcmp(colstring, "down")) return 3; 
	else mbs->UO().InstantMessageText(mystr("Warning in GetRotation: ") + colkey + mystr(" not defined!"));
	return defaultrot;	
}

Vector3D Simulink2HotintConversion::GetColorVector(MBS * mbs, ElementDataContainer* edc_block, const char* colkey, const Vector3D& coldefault)
{
	if(!edc_block->TreeFind(colkey))
		return coldefault; // return default background color
	
	// color is defined
	const char* colstring = edc_block->TreeGetString(colkey, "not_found"); 
	if(!strcmp(colstring, "black"))	return colblack;
	else if(!strcmp(colstring, "white")) return colwhite ; 
	else if(!strcmp(colstring, "red")) return colred ; 
	else if(!strcmp(colstring, "green" )) return colgreen ; 
	else if(!strcmp(colstring, "blue")) return colblue ; 
	else if(!strcmp(colstring, "cyan")) return colcyan ; 
	else if(!strcmp(colstring, "magenta")) return colmagenta ; 
	else if(!strcmp(colstring, "yellow")) return colyellow ; 
	else if(!strcmp(colstring, "gray")) return colgrey3 ; 
	else if(!strcmp(colstring, "lightBlue")) return colblue3 ; 
	else if(!strcmp(colstring, "orange")) return colorange ; 
	else if(!strcmp(colstring, "darkGreen")) return Vector3D(0, 0.5, 0.25) ; 
	else if(edc_block->TreeFind(colkey)->IsText())
	{
		mystr str(edc_block->TreeGetString(colkey));
		Vector colcustom = String2Vector(mbs, str);
		if(colcustom.Length()==3)
		{
			// custom color found
			return Vector3D(colcustom(1), colcustom(2), colcustom(3));
		}			
	}

	// should not happen!
	mbs->UO().InstantMessageText(mystr("Warning in GetColorVector: ") + colstring + mystr(" not defined!"));
	return coldefault;	
}
	
Element* Simulink2HotintConversion::ConvertBlock2Element(ElementDataContainer* edc_block)
{
	mystr name = edc_block->TreeGetString("Name","not found");
	Sim2Hotint_FNPTR pt2Function;
	int rv = 0;

	mystr blocktype = edc_block->TreeGetString("BlockType", "not found");
	for(int i = 1; i<=sim2hotint_fnptr_list.Length(); i++)
	{
		if(!strcmp(blocktype.c_str(), (*sim2hotint_identifier_list(i)).c_str()))
		{
			pt2Function = sim2hotint_fnptr_list(i);
			if(!strcmp(blocktype.c_str(), "Outport"))
			{
				nrOfOutports++;
			}
			else if(!strcmp(blocktype.c_str(), "Inport"))
			{
				nrOfInports++;
			}
			rv = 1; // function found!
		}
	}

	if(rv)
	{
		Element* el = pt2Function(mbs, edc_block, edc_default, GetScaleFactor(), GetDrawOffset());
		if(el) // el == 0 in case of scope (scope is no element!!!)
		{
			InputOutputElement* IO = ((InputOutputElement*)el);
			IO->SetDrawBackgroundColor(GetColorVector(mbs, edc_block, "BackgroundColor", colwhite));
			IO->SetDrawForegroundColor(GetColorVector(mbs, edc_block, "ForegroundColor", colblack));
			IO->SetRotation(1.0*GetRotation(mbs, edc_block, "Orientation", 0));
		}
		return el;
	}
	else
	{
		mbs->UO().InstantMessageText("ConvertBlock2Element: BlockType = " + blocktype + " not implemented!\n");
		return 0;
	}
}

// search "DstBlock", recursively search "DstBlock" if in "Branch" is other "Branch"
void Simulink2HotintConversion::parseBranch(ElementDataContainer* edc_branch, TArray<mystr*>& dstBlock, TArray<int>& dstPort, TArray<Vector2D>& dstPoints, TArray<int>& NdstPoints) //, TArray<int>& dstPointLevel)
{
	// check "DstBlock"
	
	if(edc_branch->TreeFind("DstBlock"))
	{
		//"DstBlock" found
		mystr* strPtr = new mystr(edc_branch->TreeGetString("DstBlock", "notfound"));
		dstBlock.Add(strPtr);
		mbs->UO() << "DstBlock = " << *dstBlock.Last() << " added!\n";

		// check "DstPort"
		dstPort.Add(edc_branch->TreeGetInt("DstPort", 1));
		mbs->UO() << "DstPort = " << dstPort.Last() << " added!\n";

		//----------------------------------------------
		//b: store nodes for drawing of connection lines
		NdstPoints.Add(0);
		//if(dstPointLevel.Length() == 0)dstPointLevel.Add(1); //work around
		//else dstPointLevel.Add(dstPointLevel.Last()+1);

		ElementData* ed = edc_branch->TreeFind("Points");
		if(ed)
		{		
			if(ed->IsMatrix())
			{
				double x,y;
				for(int row=1;row<=ed->GetMatrixRows();row++)
				{
					x = ed->GetMatrixVal(row,1);
					y = ed->GetMatrixVal(row,2);
					dstPoints.Add(Vector2D(x,-y));
					NdstPoints.Last()= NdstPoints.Last()-1; // -1 ... branch --> consider previous points if they exist
					
					//if(oo)mbs->UO() << "Points = " << Vector2D(x,y) << " added!\n";
				}			
			}
			else	if(ed->IsVector())
			{
				if(ed->GetVectorLen()>2)
				{
					int len = ed->GetVectorLen();
				}
				double x,y;
				edc_branch->TreeGetVector2D("Points", x, y);
				dstPoints.Add(Vector2D(x,-y));
				NdstPoints.Last()= NdstPoints.Last()-1; // -1 ... branch --> consider previous points if they exist
				//if(oo)mbs->UO() << "Points = " << Vector2D(x,y) << " added!\n";
			}
		}
		//e: store nodes for drawing of connection lines
		//----------------------------------------------		

	}

	//if(edc_branch->TreeFind("Points"))
	//{
	//	dstPoints.Add(new mystr(edc_branch->TreeGetString("Points")));
	//	mbs->UO() << "Points = " << *dstPoints.Last() << " added!\n";
	//}

	// parse sub - "Branch"
	for(int i = 1; i<=edc_branch->Length();i++)
	{
		ElementData ed = edc_branch->Get(i);
		if(!strcmp(ed.GetDataName(), "Branch" ))
		{
			if(ed.IsEDC())
			{
				parseBranch(ed.GetEDC(), dstBlock, dstPort, dstPoints, NdstPoints); //, dstPointLevel);
			}
		}
	}
}

void Simulink2HotintConversion::ConnectBlocks(ElementDataContainer* edc_line)
{
	mystr srcBlock = edc_line->TreeGetString("SrcBlock", "not found");
	int oo = 1; //test output flag
	if(oo)mbs->UO() << "srcBlock = " << srcBlock << "\n";

	TArray<mystr*> dstBlock;
	// dstPoints and dstPort have same length (each Port)
	TArray<Vector2D> dstPoints;     // "Points" ... connection nodes for drawing
	TArray<int> NdstPoints;         // number of dstPoints for the specific node, negative numbers if in sub-Branches
	//TArray<int> dstPointLevel;    // increasing level of sub-Branches (0=edc_line), not used yet!
	TArray<int> dstPort;            // port number of destination
	int srcPort;

	// check "DstBlock"

	//mystr str = edc_line->TreeGetString("DstBlock", "notFound");
	if(edc_line->TreeFind("DstBlock"))//**RL strcmp(str.c_str(), "notFound"))
	{
		//"DstBlock" found
		mystr* strPtr = new mystr(edc_line->TreeGetString("DstBlock"));
		dstBlock.Add(strPtr);
		if(oo)mbs->UO() << "DstBlock = " << *dstBlock.Last() << " added!\n";

		// check "DstPort"
		dstPort.Add(edc_line->TreeGetInt("DstPort", 1));
		if(oo)mbs->UO() << "DstPort = " << dstPort.Last() << " added!\n";
	
		//----------------------------------------------
		//b: store nodes for drawing of connection lines
		NdstPoints.Add(0);
		//dstPointLevel.Add(0); // 0=edc_line
		ElementData* ed = edc_line->TreeFind("Points");
		if(ed)
		{		
			if(ed->IsMatrix())
			{
				double x,y;
				for(int row=1;row<=ed->GetMatrixRows();row++)
				{
					x = ed->GetMatrixVal(row,1);
					y = ed->GetMatrixVal(row,2);
					dstPoints.Add(Vector2D(x,-y));
					NdstPoints.Last()= NdstPoints.Last()+1;
					if(oo)mbs->UO() << "Points = " << Vector2D(x,y) << " added!\n";
				}			
			}
			else	if(ed->IsVector())
			{
				if(ed->GetVectorLen()>2)
				{
					int len = ed->GetVectorLen();
				}
				double x,y;
				edc_line->TreeGetVector2D("Points", x, y);
				dstPoints.Add(Vector2D(x,-y));
				NdstPoints.Last()= NdstPoints.Last()+1;
				if(oo)mbs->UO() << "Points = " << Vector2D(x,y) << " added!\n";
			}
			//dstPoints.Add(new Vector2D(edc_line->TreeGetString("Points")));
		}
		//e: store nodes for drawing of connection lines
		//----------------------------------------------		
	}

	

	srcPort = edc_line->TreeGetInt("SrcPort", 1);
	if(oo)mbs->UO() << "SrcPort = " << srcPort << "!\n";	


	// parse "Branch"
	for(int i = 1; i<=edc_line->Length();i++)
	{
		ElementData ed = edc_line->Get(i);
		if(!strcmp(ed.GetDataName(), "Branch" ))
		{
			if(ed.IsEDC())
			{
				parseBranch(ed.GetEDC(), dstBlock, dstPort, dstPoints, NdstPoints); //, dstPointLevel);
				//mystr str = "dstBlock = " + *dstBlock.Last() + "\n";
				if(dstBlock.Length() == 0)	mbs->UO() << "+++++++++++ dstBlock  has Length 0!!!!\n";

			}
		}
	}

	int isrc = 0;
	int idst = 0;
	int scopeNr = 1;
	int pointNum=1;
	for(int j=1; j<=dstBlock.Length(); j++)
	{
		for(int i=1; i<=edc_blocks.Length(); i++)
		{
			if(!strcmp(edc_blocks(i)->TreeGetString("Name","not found"), srcBlock))
			{
				isrc = i;
			}
			if(!strcmp(edc_blocks(i)->TreeGetString("Name","not found"), (*dstBlock(j)).c_str()))
			{
				idst = i;
			}
		}
		if(!isrc || !idst)
		{
			mbs->UO() << "Warning in ConnectBlocks: Can't make connection!\n";
			return;
		}
		int dst_elnum = edc_blocks(idst)->TreeGetInt("MBS_elnum");
		int src_elnum = edc_blocks(isrc)->TreeGetInt("MBS_elnum");

		if(dst_elnum == -1)
		{
			mbs->UO() << "MBS: add scope " << src_elnum << " with " << dst_elnum << " and Input Port " <<  dstPort(j) << "\n";
			MBSSensor iosens(mbs, TMBSSensor(TSOutputSensor), src_elnum, srcPort); //measure output "srcPort" of  InputOutpueElement
			//iosens.SetSensorName(mystr("Scope")+mystr(scopeNr++));
			iosens.SetSensorName(edc_blocks(idst)->TreeGetString("Name",mystr(mystr("Scope")+mystr(scopeNr++)).c_str()));
			iosens.SetFactor(1.); 
			mbs->AddSensor(&iosens);
		}
		else if(dst_elnum && src_elnum)
		{
			// connect two input - output elements (from 2 "Blocks")
			mbs->UO() << "MBS: connect element " << src_elnum << " with " << dst_elnum << " and Input Port " <<  dstPort(j) << "\n";
			InputOutputElement& dst = (InputOutputElement&)mbs->GetElement(dst_elnum);
			dst.AddInput(src_elnum, IOInputTypeElement, srcPort, dstPort(j));
			
		  
 
			//---------------------------------------------------------
			// TODO: check previous levels and add input nodes for them
			//       attention: problem with "Points" entry without "DstBlock" entry not solved
			//                  j ... not same for Points and DstBlock!!!!! 
			//---------------------------------------------------------
			// drawing position of input nodes for j-th "branch"
			//
			//TArray<Vector2D> nodepos; 
			//nodepos.Add(((InputOutputElement&)mbs->GetElement(src_elnum)).GetOutputPosD(srcPort));
			//
			//int jj;
			//int level = 0; //starting level
			//for(jj=j;jj>=1;jj--) 
			//{
			//	if(dstPointLevel(jj) == level)
			//		break;  // jj ... line level ("source")
			//}

			//for(;jj<=j;jj++) // count up from line level =0 to actual branch level
			//{
			//	if(NdstPoints(jj) > 0 && dstPointLevel(jj) == level + 1)
			//	{
			//		if(dstPointLevel(jj) == dstPointLevel(j))
			//		{	
			//			jj = j; // don't consider parallel branch
			//		}
			//		int firstnode=1;
			//		for(int jjj=1;jjj<=jj;jjj++)
			//		{
			//			firstnode += NdstPoints(jjj);
			//		}

			//		nodepos(firstnode) = nodepos(firstnode) + scalefactor*dstPoints(firstnode);
			//		for(int node=firstnode+1;node<=firstnode+NdstPoints(jj);node++)
			//		{
			//			nodepos.Add(nodepos(node-1) + scalefactor*dstPoints(node));  // from source to destination
			//		}


			//		//	NdstPoints.Add(0);
			//		//	dstPointLevel.Add(0); // 0=edc_line
			//	}
			//}

			// code works for simple drawing of nodes
			TArray<Vector2D> nodepos;
			if(NdstPoints(j) > 0)
			{
				nodepos.Add(((InputOutputElement&)mbs->GetElement(src_elnum)).GetOutputPosD(srcPort));
				nodepos(1) = nodepos(1) + scalefactor*dstPoints(1);
				for(int node=2;node<=NdstPoints(j);node++)
				{
					nodepos.Add(nodepos(node-1) + scalefactor*dstPoints(node));  // from source to destination
				}
			}
			// add input nodes for drawing
			for(int node=nodepos.Length();node>0;node--)
			{
				dst.AddInputNode(dstPort(j), nodepos(node)); 
			}

			pointNum+=NdstPoints(j);
		}
	}

	// delete allocated mystr*
	for (int i=1; i <= dstBlock.Length(); i++)
	{
		if (dstBlock(i) != 0) delete (*dstBlock(i));
		dstBlock(i) = 0;
	}
	dstBlock.SetLen(0);
}

// returns element number of inport
int Simulink2HotintConversion::GetInportMBSElementNumber(int inport_num)
{
	ElementData ed;
	for(int i=1; i<=edc_blocks.Length(); i++)
	{
		ed = edc_blocks(i)->Get(1);
		int inport_elnum = edc_blocks(i)->TreeGetInt("MBS_elnum");
		if(!strcmp(edc_blocks(i)->TreeGetString("BlockType", "not found"), "Inport") && inport_elnum)
		{
			// Inport element found
			int block_port_num = 1;
			int index = edc_blocks(i)->Find("Port");
			if(index) 
			{
				ElementData* ed_port = edc_blocks(i)->GetPtr(index);
				if (ed_port->IsText()) 
				{
					block_port_num = mystr(ed_port->GetText()).MakeInt();
				}
			}
			if(block_port_num == inport_num)
			{
				return inport_elnum;
			}
		}
	}
	return 0;
}

// returns element number of inport
int Simulink2HotintConversion::GetOutportMBSElementNumber(int outport_num)
{
	ElementData ed;
	for(int i=1; i<=edc_blocks.Length(); i++)
	{
		ed = edc_blocks(i)->Get(1);
		int outport_elnum = edc_blocks(i)->TreeGetInt("MBS_elnum");
		if(!strcmp(edc_blocks(i)->TreeGetString("BlockType", "not found"), "Outport") && outport_elnum)
		{
			// Inport element found
			//mystr port = edc_blocks(i)->TreeGetString("Port", "1"); // if "Port" not found, Simulink uses "1" as default value
			mystr port("1");
			if(edc_blocks(i)->TreeFind("Port") && edc_blocks(i)->TreeFind("Port")->GetText())
			{
				port = edc_blocks(i)->TreeGetString("Port", "1"); // if "Port" not found, Simulink uses "1" as default value
			}
			if(port.MakeInt() == outport_num)
			{
				return outport_elnum;
			}
		}
	}
	return 0;
}


// global functions
// parser of string 2 double
//
// mbs ... multibody system
// str_exp ... contains expressions of variables
// edc_par .. edc, where values of variables of mdl - file are defined (now: expressions like "MDL.a = 3" and "MDL.b * c=4" are also possible)
// operations +-/*, sin,... evaluated by ExpressionToDouble
double String2Double(MBS * mbs, mystr& str_exp, ElementDataContainer* edc_par)
{
	//if(1)
	//{
		// use parser
		return mbs->ExpressionToDouble(str_exp);
	//}

	//ElementDataContainer* edc = edc_par;
	//if(!edc){	edc = mbs->GetModelDataContainer();} // MDL-edc is in model - edc
	//else{	edc = edc_par; } // not used yet, maybe used in case of more than one subsystem with different edc_par for each mdl-subsystem (e.g. use same mdl-file with different parameters)	

	//double val = 0.;
	//if(edc->TreeFind(mystr(mystr("MDL.") + str_exp).c_str()))
	//{
	//	return edc->TreeGetDouble(mystr(mystr("MDL.") + str_exp).c_str(),0.);
	//}
	//else
	//{
	//	{
	//		int first = 1;
	//		int last = str_exp.Length();
	//		if(!IsCorrectNumberBetween(str_exp, first, last))			
	//		{
	//			// string doesn't hold a correct number
	//			mbs->UO(UO_LVL_err).InstantMessageText("Error: Could not parse " + str_exp + " correctly. Please add expression MDL." + str_exp + " to get correct simulation results.");
	//		}
	//	}
	//	val = str_exp.MakeDouble();
	//}
	//return val;
}

Vector String2Vector(MBS * mbs, mystr& str_exp, char osbc, char csbc, char commac, ElementDataContainer* edc_par) //$ RL 2011-01 evaluate string str and make vector
{
	//mystr must be formatted like following examples:  [1 2 4] or [1,2,3]
	mystr el = "\n"; //end of line
	char elc = 'n';
	mystr osb = osbc;  //open square bracket, for vector
	mystr csb = csbc;  //closing square bracket
	mystr semicolon = ";"; //double stop
	mystr comma = commac; //comma
	int endstr = 0;
	int pos = 0;
	mystr dataname;
	char ch;

	mystr dummycomment; //dummy comment string
	mystr text = str_exp.GetStringInBrackets(pos, osbc, csbc);

	TArray<double> ta_vec(0);
	if(text.Find(el)==-1  && text.Find(semicolon)==-1)
	{
		// vector has no semicolon and is written in one single line
		mystr s_value;						
		int posvec = 0;
		int possval = 0;
		int comma_found = 0;
		while(posvec != -1)
		{
			comma_found = text.GetUntil(posvec, commac, s_value); // 'value,'
			s_value.ReadLeadingSpaces(possval);
			//same as in class Vector ta_vec.Add(s_value.MakeDouble());
			ta_vec.Add(String2Double(mbs, s_value));
			if(comma_found)posvec++;
		}											
	}			
	return Vector(ta_vec);

}

Element* ConvertSim2Hotint_Gain(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	double gain = 0.0;
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Gain"))
	{
		// get gain value
		ElementData* ed = edc->TreeFind("Gain");

		// search in default edc
		if(!ed)
		{

			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),"Gain"))
					{
						ed = edc_def_block->TreeFind("Gain");
					}
				}
			}
		}

		if(ed)
		{
			if(ed->IsText())
			{
				mystr val = ed->GetText(); 
				gain = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_Gain!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_Gain!");
			return 0;
		}

		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_Gain!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_Gain!");
		return 0;
	}


	LinearTransformation IO(mbs); 
	IO.SetGain(gain);
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	return mbs->GetElementPtr(nIO);

}

Element* ConvertSim2Hotint_Inport(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	double gain = 1.0;
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Inport"))
	{
		// get gain value
		ElementData* ed = edc->TreeFind("Inport");

		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning1 in ConvertSim2Hotint_Inport!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_Inport!");
		return 0;
	}


	LinearTransformation IO(mbs); 
	IO.SetGain(gain);
	IO.SetElementName(edc->TreeGetString("Name","not found"));    
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	//IO.AddInput(nsensPhi1, 2, 1);
	int nIO = mbs->AddElement(&IO);


	ElementData ed;
	//ed.SetInt(nIO, "MBS_inport_elnum"); edc->Add(ed);    
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    


	//input and outputs set to 0;
	//IOLinearTransformation io(mbs, ....);

	return mbs->GetElementPtr(nIO);

}

Element* ConvertSim2Hotint_Outport(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	double gain = 1.0;
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Outport"))
	{	
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning1 in ConvertSim2Hotint_Outport!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_Outport!");
		return 0;
	}


	LinearTransformation IO(mbs); 
	IO.SetGain(gain);//TODO: edc->TreeGetDouble();
	IO.SetElementName(edc->TreeGetString("Name","not found"));    
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	//IO.AddInput(nsensPhi1, 2, 1);
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    

	//input and outputs set to 0;
	//IOLinearTransformation io(mbs, ....);

	return mbs->GetElementPtr(nIO);

}

Element* ConvertSim2Hotint_Sum(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{


	// EDC must be Block-EDC
	TArray<int> input_signs;
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
	TArray<double> inputs;

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Sum"))
	{	
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)
		// get gain value
		ElementData* ed = edc->TreeFind("Inputs");
		if(!ed)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),"Sum"))
					{
						ed = edc_def_block->TreeFind("Inputs");
					}
				}
			}
		}
		if(ed)
		{
			if(ed->IsText())
			{
				mystr val = ed->GetText(); // |++--
				for(int i=0; i<val.Length(); i++)
				{
					char c = val.PosPeek(i);
					if(c == '+')
					{
						inputs.Add(1.0);
					}
					else if(c == '-')
					{
						inputs.Add(-1.0);
					}
				}
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_Sum!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_Sum!");
			return 0;
		}


		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_Sum!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_Sum!");
		return 0;
	}

	LinearTransformation IO(mbs);
	IO.SetAdder(Vector(inputs));
	IO.SetElementName(edc->TreeGetString("Name","not found"));    
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	//IO.AddInput(nsensPhi1, 2, 1);
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    

	//input and outputs set to 0;
	//IOLinearTransformation io(mbs, ....);

	return mbs->GetElementPtr(nIO);

}


// global functions
Element* ConvertSim2Hotint_Integrator(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	//double initialCondition = 0.0; // TODO
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
	double initialCondition = 0.0;
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Integrator"))
	{

		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_Integrator!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down

		if(edc->TreeFind("InitialCondition"))
		{
			mystr val = edc->TreeGetString("InitialCondition","0.0");
			initialCondition = String2Double(mbs, val);
		}
	}
	else
	{
		mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_Integrator!");
		return 0;
	}

	STransferFunction IO(mbs); 
	IO.SetSTransferFunction(Vector(1.,0.),	Vector(0.,1.), Vector(initialCondition));
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	//IO.AddInput(nsensPhi1, 2, 1);
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	//input and outputs set to 0;
	//IOLinearTransformation io(mbs, ....);

	return mbs->GetElementPtr(nIO);

}

Element* ConvertSim2Hotint_Constant(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	double value = 0.0;
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),"Constant"))
	{
		// get value
		ElementData* ed = edc->TreeFind("Value");

		// search in default edc
		if(!ed)
		{
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),"Constant"))
					{
						ed = edc_def_block->TreeFind("Value");
					}
				}
			}
		}

		if(ed)
		{
			if(ed->IsText())
			{
				mystr val = ed->GetText(); 
				value = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_Constant!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_Constant!");
			return 0;
		}

		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_Constant!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_Constant!");
		return 0;
	}


	LinearTransformation IO(mbs); 
	IO.SetConstant(value);//TODO: edc->TreeGetDouble();
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	//IO.AddInput(nsensPhi1, 2, 1);
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	//input and outputs set to 0;
	//IOLinearTransformation io(mbs, ....);

	return mbs->GetElementPtr(nIO);

}
		
Element* ConvertSim2Hotint_Saturate(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{	
	// EDC must be Block-EDC
	double value_up = 0.0;   // lower limit
	double value_low = 0.0;  // upper limit

	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
  mystr blockName("Saturate");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_up = edc->TreeFind("UpperLimit");
		ElementData* ed_low = edc->TreeFind("LowerLimit");
		if(!ed_up || !ed_low)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_up)ed_up = edc_def_block->TreeFind("UpperLimit");
						if(!ed_low)ed_low = edc->TreeFind("LowerLimit");
					}
				}
			}
		}

		if(ed_up)
		{
			if(ed_up->IsText())
			{
				mystr val = ed_up->GetText();				
				value_up = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		if(ed_low)
		{
			if(ed_low->IsText())
			{
				mystr val = ed_low->GetText(); 
				value_low = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}


	// old:
	//MathFunction mf;
	//Vector xData(-1.0e30, value_low, value_up, 1.0e30);
	//Vector yData(value_low, value_low, value_up, value_up);
	//mf.SetPiecewise(xData, yData, 1); // piecewise linear
	//IOMathFunction IO(mbs); 
	//IO.SetIOMathFunction(mf);
  // new:
	IOSaturate IO(mbs);
	IO.SetIOSaturate(value_low, value_up);
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	
	return mbs->GetElementPtr(nIO);
}



Element* ConvertSim2Hotint_TransferFcn(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	//double initialCondition = 0.0; // TODO
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
  mystr blockName("TransferFcn");
	Vector num(0.);
	Vector den(0.);

	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		//BlockType		      TransferFcn
		//	Numerator		      "[1]"
		//	Denominator	      "[1 2 1]"
		// get value
		ElementData* ed_num = edc->TreeFind("Numerator");
		ElementData* ed_den = edc->TreeFind("Denominator");

		// search in default edc
		if(!ed_num || !ed_den)
		{
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_num)ed_num = edc_def_block->TreeFind("Numerator");
						if(!ed_den)ed_den = edc_def_block->TreeFind("Denominator");
					}
				}
			}
		}
		
		if(ed_num)
		{
			if(ed_num->IsText())
			{
				mystr val = ed_num->GetText();
				num = String2Vector(mbs, val, '[', ']', ' '); // space ' ' is delimiter of mdl-file vector
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		if(ed_den)
		{
			if(ed_den->IsText())
			{
				mystr val = ed_den->GetText();
				den = String2Vector(mbs, val, '[', ']', ' '); //Vector(val, '[', ']', ' '); // space ' ' is delimiter of mdl-file vector
			}
			else
			{
				mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}
  
	// add zeros, if numerator has lower lenght than denominator
	Vector numCorr(den.Length());
	int delta_exp = den.Length()-num.Length();
	for(int i = 1; i<=den.Length(); i++)
	{
		if(i <= delta_exp)
		{
			numCorr(i) = 0.;
		}
		else
		{
			numCorr(i) = num(i-delta_exp);
		}
	}

	STransferFunction IO(mbs);	
  IO.SetSTransferFunctionWithDescentPolynoms(numCorr,	den);//TODO: add initial Condition
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);
	
	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	
	return mbs->GetElementPtr(nIO);

}


Element* ConvertSim2Hotint_DeadZone(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset) //$ RL 2011-01: DeadZone (e.g.: "clearance")
{
	// EDC must be Block-EDC
	double value_up = 0.0;   // lower limit
	double value_low = 0.0;  // upper limit

	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
  mystr blockName("DeadZone");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_up = edc->TreeFind("UpperValue");
		ElementData* ed_low = edc->TreeFind("LowerValue");
		if(!ed_up || !ed_low)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_up)ed_up = edc_def_block->TreeFind("UpperValue");
						if(!ed_low)ed_low = edc_def_block->TreeFind("LowerValue");
					}
				}
			}
		}

		if(ed_up)
		{
			if(ed_up->IsText())
			{
				mystr val = ed_up->GetText();				
				value_up = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		if(ed_low)
		{
			if(ed_low->IsText())
			{
				mystr val = ed_low->GetText(); 
				value_low = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}

	IODeadZone IO(mbs);
	IO.SetIODeadZone(value_low, value_up);
	IO.SetElementName(edc->TreeGetString("Name","not found"));
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	
	return mbs->GetElementPtr(nIO);

}


Element* ConvertSim2Hotint_Scope(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)//$ RL 2011-01: scope is implemented as MBSSensor
{
	//return zero-pointer, because MBSSensor is no element
	ElementData ed;
	ed.SetInt(-1, "MBS_elnum"); edc->Add(ed);  //-1 ... sensor
	return 0;
}

Element* ConvertSim2Hotint_TransportDelay(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)//$ RL 2011-01: continuous delays of mdl-files are approximated as time discrete delays (in simulink, a buffer with old time variables is used)
{
  // EDC must be Block-EDC
	double value = 1.0;   // time delay (default 1.0 s)
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
  mystr blockName("TransportDelay");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_delay = edc->TreeFind("DelayTime");
		if(!ed_delay)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_delay)ed_delay = edc_def_block->TreeFind("DelayTime");
					}
				}
			}
		}

		if(ed_delay)
		{
			if(ed_delay->IsText())
			{
				mystr val = ed_delay->GetText();				
				value = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
	
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)
		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}
	int zexp = 1;         // exponent of z (number of discrete delays for approximation of continous delay)
	double maxTime = mbs->GetOptions()->PreProcOptions()->SimulinkModelFormatParserApproxPeriodCont2disc();
	if(value > maxTime)
	{
		zexp = ceil(value/maxTime); // approximize continuous time delays with time discrete delays with maximal delay "approx_period_cont2disc"
	}
 	assert(zexp>0);       // exponent of denominator must be greater than zero
	Vector num(1+zexp); num(1+zexp) = 1; //  1        num   
	//                                   //------  =  ---
	Vector den(1+zexp); den(1) = 1.;     //z^zexp     den

	ZTransferFunction IO(mbs);
	IO.SetSampleTime(value/(double)zexp); // whole delay "value" is discretized with zexp delays
	IO.SetZTransferFunctionWithDescentPolynoms(num,den); // numerator and denominator are sorted in descent order ==> [..., z^1, z^0]
	IO.SetElementName(edc->TreeGetString("Name","not found") + mystr(" 1/(z^") + mystr(zexp) + mystr(")"));
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	
	return mbs->GetElementPtr(nIO);
}

//$ RL 2011-05:[
Element* ConvertSim2Hotint_Product(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// Multiplication:	      "Matrix(*)" not implemented yet (only 1 dimensional signals <=> Elementwise)
	// Inputs		      "/"   ... division not implemented yet
  //
	// Inputs         "***"   ... 1D-product implemented
	// Inputs		      "3" ==> ... 1D-product implemented
  
  // EDC must be Block-EDC
	// Product of input(i)^operator(i)
	TArray<double> exp; // exponents: 1 <=> multiply, -1 <=> division
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
  
	
	mystr blockName("Product");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_inputs = edc->TreeFind("Inputs");
		ElementData* ed_multtype = edc->TreeFind("Multiplication");
		if(!ed_inputs || !ed_multtype)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_inputs) ed_inputs = edc_def_block->TreeFind("Inputs");
						if(!ed_multtype) ed_multtype = edc_def_block->TreeFind("Multiplication");
					}
				}
			}
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++
		//TODO: implement Matrix operations of Block "Product"
		mystr elementwise_text("Element-wise(.*)");
		if(!elementwise_text.Compare(mystr(ed_multtype->GetText())))
		{
			mbs->UO().InstantMessageText("Warning1 in ConvertSim2Hotint_" + blockName + ": Multiplication of type " + edc->TreeGetString("Multiplication") + " not implemented yet (only elementwise product).");
			return 0;	
		}
		
		mystr name("");
		if(ed_inputs)
		{
			if(ed_inputs->IsText())
			{
				// in "Inputs", the exp "*" and/or "/" are stored
				mystr val = ed_inputs->GetText();
				int j = 0;
				if(val.PosPeek(j) == '*' || val.PosPeek(j) == '/')
				{
					for(int i = 0; i< val.Length(); i++)
					{
						if(val.PosPeek(i) == '*')
						{
							exp.Add(1.); // multiplication
							name = name+"*";
						}
						else if(val.PosPeek(i) == '/')
						{
							exp.Add(-1.); // division
							name = name+"/";
						}
						else
						{
							mbs->UO().InstantMessageText(mystr("Warning2 in ConvertSim2Hotint_") + blockName + mystr(": ") + mystr(val.PosPeek(i)) + mystr(" not expected in Block \"Product\"."));
//						return 0;	 // maybe not possible (in case of spaces or tabs between * and /)
						}
					}
				}
				else
				{
					// in "Inputs", the number of multipliers is stored
					int value = String2Double(mbs, val);
					exp.SetLen(value);
					exp.SetAll(1.);
					for(int i=1;i<=value;i++)
					{
						name = name + "*";
					}
				}
			}
			else
			{
				mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}			

			// get position vector (upper left corner x1,y1; lower right corner x2, y2)
			Vector pos = EDCTreeGetVector(*edc, "Position");
			if(pos.Length() != 4)
			{
				mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
			x1 = pos(1);
			y1 = -pos(2); //Simulink y-Values are top/down
			x2 = pos(3);
			y2 = -pos(4); //Simulink y-Values are top/down
		}
		else
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		


		IOProduct IO(mbs);
		IO.SetIOProduct(Vector(exp));		
		IO.SetElementName(edc->TreeGetString("Name","not found") + mystr("") + name);
		IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
		IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
		int nIO = mbs->AddElement(&IO);
		ElementData ed;
		ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
		return mbs->GetElementPtr(nIO);
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}

}

Element* ConvertSim2Hotint_Fcn(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
// EDC must be Block-EDC
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
	mystr expr; // expression 
  mystr blockName("Fcn");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_expr = edc->TreeFind("Expr");
		if(!ed_expr)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_expr)ed_expr = edc_def_block->TreeFind("Expr");
					}
				}
			}
		}

		
		if(ed_expr)
		{
			if(ed_expr->IsText())
			{
				expr = ed_expr->GetText();
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}

	if(expr.Find('[')!=-1 || expr.Find(']')!=-1)
	{
		//TODO: only expressions with one input "u" implemented yet. u[1] is not valid yet.
		mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + ": '[' or ']' found, please rename u[1] to u in mdl-file! Further inputs are not implemented yet.");		
		return 0;
	}
	
	MathFunction mf;
	mf.SetExpression(expr, "u",mbs);

	IOMathFunction IO(mbs);
	IO.SetIOMathFunction(mf);
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    	
	return mbs->GetElementPtr(nIO);
}

// compute friction-Mathfunction
void StaticDynamicFrictionMathFunc(const double staticFriction, const double viscousFriction, const double zeroZone, MathFunction& mf)
{
	double inf = 1e10; //oo
	double limit = max(100*zeroZone, 10000.0); // value 'limit' must have bigger magnitude than zeroZone and unequal zero;  values with bigger magnitude than 'limit' are extrapolated by MathFunction

	Vector omega(4);
	omega(1) = -limit; 
	omega(2) = -zeroZone;
	omega(3) = zeroZone;
	omega(4) = limit;
	
	if(zeroZone < 1/inf && staticFriction < 1/inf)
	{
		// just viscous friction
		omega.SetLen(2);
		omega = Vector(-inf, inf);
	}
	Vector M_frict(omega.GetLen());
	for(int k = 1; k<= M_frict.Length(); k++)
	{
		M_frict(k) = staticFriction*Sgn(omega(k)) + viscousFriction*omega(k);
	}
	mf.SetPiecewise(omega, M_frict, 1);
}

Element* ConvertSim2Hotint_Reference(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)
{
	// EDC must be Block-EDC
	double gain = 0.0; // viscous friction
	double offset = 0.0; // static friction
	double zeroZone = mbs->GetOptions()->PreProcOptions()->SimulinkModelFormatParserDefaultValuesStaticDynamicFricionZeroZone();
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

  mystr blockName("Reference");
	if(!strcmp(edc->TreeGetString("Blocktype","not found"),blockName))
	{
		// get value
		ElementData* ed_gain = edc->TreeFind("gain");
		ElementData* ed_offset = edc->TreeFind("offset");
		ElementData* ed_sourceType = edc->TreeFind("SourceType");
		if(!ed_gain || !ed_offset || !ed_sourceType)
		{
			// search in default edc
			for(int i=1;i<=edc_default->Length();i++)
			{
				if(edc_default->Get(i).IsEDC())
				{
					ElementDataContainer* edc_def_block = edc_default->Get(i).GetEDC();
					if(!strcmp(edc_def_block->TreeGetString("Blocktype","not found"),blockName))
					{
						if(!ed_gain)ed_gain = edc_def_block->TreeFind("gain");
						if(!ed_offset)ed_offset = edc->TreeFind("offset");
						if(!ed_sourceType)ed_sourceType = edc->TreeFind("SourceType");
					}
				}
			}
		}

		if(!mystr(ed_sourceType->GetText()).Compare(mystr("Coulombic and Viscous Friction")))
		{
			mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "! Unknown SourceType of Block \"Reference\"" + mystr(ed_sourceType->GetText()) + "!");
			return 0;
		}

		if(ed_gain)
		{
			if(ed_gain->IsText())
			{
				mystr val = ed_gain->GetText();				
				gain = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}

		if(ed_offset)
		{
			if(ed_offset->IsText())
			{
				mystr val = ed_offset->GetText(); 
				offset = String2Double(mbs, val);
			}
			else
			{
				mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");
				return 0;
			}
		}
		else
		{
			mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		// get position vector (upper left corner x1,y1; lower right corner x2, y2)

		Vector pos = EDCTreeGetVector(*edc, "Position");
		if(pos.Length() != 4)
		{
			mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
			return 0;
		}
		x1 = pos(1);
		y1 = -pos(2); //Simulink y-Values are top/down
		x2 = pos(3);
		y2 = -pos(4); //Simulink y-Values are top/down
	}
	else
	{
		mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}

	MathFunction mf;
	StaticDynamicFrictionMathFunc(offset, gain, zeroZone, mf);
	
	IOMathFunction IO(mbs);
	IO.SetIOMathFunction(mf);
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	IO.SetSwitchOnlyInPostNewton(1);
	int nIO = mbs->AddElement(&IO);

	ElementData ed;
	ed.SetInt(nIO, "MBS_elnum"); edc->Add(ed);    
	return mbs->GetElementPtr(nIO);
}
//$ RL 2011-05:]

//$ RL 2011-7-12:[ use new EDC-functions for creating element (for new parse function, copy this function and change code after comments marked with '//**')
Element* ConvertSim2Hotint_DiscretePulseGenerator(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)//$ RL 2011-01: scope is implemented as MBSSensor
{ 
	// ==============================
	// assumptions: VectorParams1D on
	//              use in case "sample based" also the flag "use simulation time" - otherwise "SampleTime" is not found!
	// ==============================

	//** define block name
	mystr blockName("DiscretePulseGenerator");
	// search default edc
	ElementDataContainer* edc_def_block;
	for(int i=1;i<=edc_default->Length();i++)
	{
		if(edc_default->Get(i).IsEDC() &&  !strcmp(edc_default->Get(i).GetEDC()->TreeGetString("Blocktype","not found"),blockName))
		{
			edc_def_block = edc_default->Get(i).GetEDC();
		}
	}

	// replace defaults with edc-data
	edc_def_block->TreeReplaceEDCDataWith(edc);
	
	//** read variables
	double amplitude = 0.;    // amplitude of rectangle pulse
	double toffs = 0.;        // offset time for pulse sequence
	double period = 1.;       // period of rectangle
	double pulseWidth = 0.;   // pulse width - output is set to amplitude in this time span
	int useExternalTime = 0;  // set to nonzero value if input should be used as time source; otherwise use simulation time

	ElementData* ed = edc_def_block->TreeFind("Amplitude");
	if(ed && ed->IsText()){amplitude = String2Double(mbs, mystr(ed->GetText()));}
	else {mbs->UO().InstantMessageText("Warning1 in ConvertSim2Hotint_" + blockName + "!");return 0;}

	ed = edc_def_block->TreeFind("PhaseDelay");
	if(ed && ed->IsText()){toffs = String2Double(mbs, mystr(ed->GetText()));}
	else {mbs->UO().InstantMessageText("Warning2 in ConvertSim2Hotint_" + blockName + "!");return 0;}

	ed = edc_def_block->TreeFind("Period");
	if(ed && ed->IsText()){period = String2Double(mbs, mystr(ed->GetText()));}
  else {mbs->UO().InstantMessageText("Warning3 in ConvertSim2Hotint_" + blockName + "!");return 0;}

	ed = edc_def_block->TreeFind("PulseWidth");
	if(ed && ed->IsText()){pulseWidth = String2Double(mbs, mystr(ed->GetText()));}
	else {mbs->UO().InstantMessageText("Warning4 in ConvertSim2Hotint_" + blockName + "!");return 0;}

	ed = edc_def_block->TreeFind("TimeSource");
	if(ed && ed->IsText())
	{
		if(mystr(ed->GetText()).Compare("Use simulation time")){useExternalTime = 0;} // use simulation time (no input)
		else {useExternalTime = 1;} // use input port for time		
	}
	else
	{
		mbs->UO().InstantMessageText("Warning5 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}
	
	ed = edc_def_block->TreeFind("PulseType");
	if(ed && ed->IsText())
	{
		if(mystr(ed->GetText()).Compare(mystr("Sample based")))
		{
			// use percentage instead of times
			if(!useExternalTime)
			{
				ed = edc_def_block->TreeFind("SampleTime");
				if(ed && ed->IsText())
				{
					double sampletime = String2Double(mbs, mystr(ed->GetText()));
								
					period = ((int)period) * sampletime;
					pulseWidth = ((int)pulseWidth) * sampletime;
					toffs	= ((int)toffs) * sampletime;
				}
				else {mbs->UO().InstantMessageText("Warning6 in ConvertSim2Hotint_" + blockName + "!");return 0;}

			}
			else
			{
				// not implemented yet (unknown SampleTime)!
				mbs->UO().InstantMessageText("Warning7 in ConvertSim2Hotint_" + blockName + "!");return 0;
			}
		}
		else if(mystr(ed->GetText()).Compare(mystr("Time based")))
		{
			// nothing changes (time based parameters)
			pulseWidth = pulseWidth/100. * period; // pulse width in percent of period  --> (s)
		}
		else
		{
			mbs->UO().InstantMessageText("Warning8 in ConvertSim2Hotint_" + blockName + "!");return 0;
		}		
	}
	else {mbs->UO().InstantMessageText("Warning9 in ConvertSim2Hotint_" + blockName + "!");return 0;}

	// drawing parameters
	Vector pos = EDCTreeGetVector(*edc, "Position");
	if(pos.Length() != 4)
	{
		mbs->UO().InstantMessageText("Warning10 in ConvertSim2Hotint_" + blockName + "!");
		return 0;
	}
	double x1 = 0.0;
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;

	x1 = pos(1);
	y1 = -pos(2); //Simulink y-Values are top/down
	x2 = pos(3);
	y2 = -pos(4); //Simulink y-Values are top/down

	//** create element
	IOPulseGenerator IO(mbs);
	IO.SetIOPulseGenerator(amplitude, toffs, period,	pulseWidth, useExternalTime);
	IO.SetElementName(edc->TreeGetString("Name","not found"));  
	IO.SetRefPos2D(drawoffset+scalefactor*Vector2D((x1+x2)/2.0,(y1+y2)/2.0));
	IO.SetDrawDim(scalefactor*Vector3D(fabs(x2-x1),fabs(y2-y1),0.));
	int nIO = mbs->AddElement(&IO);

	// store element number in edc
	ed->SetInt(nIO, "MBS_elnum"); edc->Add(*ed);
	return mbs->GetElementPtr(nIO);
}
//$ RL 2011-7-12:]