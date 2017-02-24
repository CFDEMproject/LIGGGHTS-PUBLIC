//#**************************************************************
//#
//# filename:             mbs_communication.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						March 2010
//# description:          Load/Save mbs data, WCDriver and Windows communication 
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
 
#include "windows.h" //for shell execute
#include <direct.h>  //for getcwd
#include  <io.h>     //file operation _access
#include  <stdio.h>
#include  <stdlib.h>

#include "mbs.h"
#include "parser.h"
#include "script_parser.h"
#include "myfile.h"
#include "element.h"
#include "sensors.h"
#include "node.h"
#include "material.h"
#include "elementdataaccess.h"


int MultiBodySystem::SaveSolutionVector(const mystr& filename)  //rv=1 --> OK
{
	int opt = ios::out;


	ofstream* storeout = new ofstream;
#ifdef my_new_stdiostream
	storeout->open(filename.c_str(), opt); //with read sharing allowed!
#else
	storeout->open(filename.c_str(), opt, filebuf::sh_read); //with read sharing allowed!
#endif

	if ((int)storeout->good() == 0 || storeout->fail())  
	{
		UO() << "ERROR: could not write solution!!!\n";
		UO() << "good=" << storeout->good() << ", fail=" << storeout->fail() << "\n";
		UO() << "file=" << filename.c_str() << "\n";

		//storeout->close();
		//storeout->delete;
		//return 0;
	}
	storeout->precision(19);
	int len = GetSolVector().Length();

	char str[256];
	sprintf(str, "HOTINTDataVersion%.3f", GetHotintVersion().GetDoubleValue());
	(*storeout) << str << "\n";
	(*storeout) << len << " 0" << "\n"; //checksums!
	(*storeout) << "1\n"; //only one step stored

	char strtime[64];
	char strtime2[64];
	sprintf(strtime, "%.14f", GetTime());
	sprintf(strtime2, "%.17f", GetTime());

	if (strlen(strtime) < strlen(strtime2)-4)
		(*storeout) << strtime << "\n";
	else
		(*storeout) << strtime2 << "\n";

	(*storeout) << len << "\n";
	for (int i = 1; i <= len; i++)
	{
		(*storeout) << GetSol(i) << "\n";
	}

	storeout->close();
	delete storeout;

	return 1;
}

int MultiBodySystem::LoadInitialVector(const mystr& filename)  //rv=1 --> OK
{
	CMFile storein(filename, TFMread);
	Vector x0(GetSystemSize());

	stored_initialconditions = GetSolVector();

	if (!storein.IsGood())
	{
		UO() << "ERROR: could not open stored solution!!!\n";
		return 0;
	}
	else
	{
		//storein->seekg(0,ios::beg);

		//read either version or length:
		mystr header;
		storein.RWSoleStr(header);

		if (header.Length() > 17 && header.Left(17) == mystr("HOTINTDataVersion"))
		{
			//UO() << "Load initial vector with header\n"; //--> e.g. stored in Data manager

			mystr version_str = header.Right(header.Length()-17);
			//UO() << "Version number ='" << version << "'\n";

			HotintVersionInfo data_version(atof(version_str));

			if (data_version > GetHotintVersion())
			{
				UO() << "Warning: This is a solution vector of another HOTINT version, the data may be corrupt!\n";
			}

			int ndof, checksum2, nsteps;
			storein.RWInt(ndof);
			storein.RWInt(checksum2);
			storein.RWInt(nsteps);

			if (ndof != x0.Length()) {UO() << "Warning: Stored solution does have different size!!!\n";}
			if (nsteps != 1) {UO() << "Warning: Stored solution contains " << nsteps << " data units, only the first data unit is loaded!\n";}
		}

		double time;
		int len;

		storein.RWDouble(time);
		storein.RWInt(len);

		//UO() << "time=" << time << ", len=" << len << "\n";

		if (len > x0.Length()) {UO() << "ERROR: Stored solution is incompatible!!!\n";}
		else
		{
			for (int i=1; i <= len; i++) //additionally check for storein->good() ....???
			{
				storein.RWDouble(x0(i));
			}
			stored_initialconditions = x0;
			//UO() << "WARNING: solution vector added!!!!!\n";

			//UO() << "Stored solution loaded!\n";
		}
	}

	return 1;
}


//type=1: element
//type=2: force	i				//$ DR 2012-10 old:force value of element i
//type=3: sensor i
//type=4: geomelement i
//type=5: node i
//type=6: material i
//type=50: KB options (old)
void MultiBodySystem::GetElementData(int i, int type, int value, ElementDataContainer& edc) 
{
	switch (type)
	{
	case 1: //element data
		{
			GetElement(i).GetElementData(edc);
			break;
		}
	case 2: //force value of element i, force value
		{
			//GetElement(i).GetLoad(value).GetElementData(edc); //$ DR 2012-10 old code
			GetLoad(i).GetElementData(edc);
			break;
		}
	case 3: //sensor i data
		{
			GetSensor(i).GetElementData(edc);
			break;
		}
	case 4: //geomelement i data
		{
			GetDrawElement(i)->GetElementData(edc);
			break;
		}
	case 5: //Node i data
		{
			GetNode(i).GetElementData(edc);
			break;
		}
	case 6: //Material i data
		{
			GetMaterial(i).GetElementData(edc);
			break;
		}
	default: ;
	}
}

int MultiBodySystem::SetElementData(int i, int type, int value, ElementDataContainer& edc) 
{
	int rv = 1;
	switch (type)
	{
	case 1: //element data
		{
			rv = GetElement(i).SetElementData(edc);
			break;
		}
	case 2: //force value of element i, force value
		{
			//rv = GetElement(i).GetLoad(value).SetElementData(edc); //$ DR 2012-10 old code
			rv = GetLoad(i).SetElementData(edc); 
			break;
		}
	case 3: //sensor i data
		{
			rv = GetSensor(i).SetElementData(edc);
			break;
		}
	case 4: //geomelement i data
		{
			rv = GetDrawElement(i)->SetElementData(edc);
			break;
		}
	case 5: //Node i data
		{
			rv = GetNode(i).SetElementData(edc);
			break;
		}
	case 6: //Node i data
		{
			rv = GetMaterial(i).SetElementData((ElementDataContainer &)edc);
			break;
		}
		//#ifdef KB_INCLUDES
		//	case 50: //edit KB menu
		//		{
		//			kbmenu.SetElementData(edc);
		//			break;
		//		}
		//#endif
	default: ;
	}
	return rv;
}

//$AD 2013-07-08: Manipulate Content of arrays in IOElements from 2D Draw Window
void MultiBodySystem::InsertIOElemConNode(int elemnr, int list_idx, int input_nr, Vector2D pos)
{
	GetElement(elemnr).InsertConNode(list_idx, input_nr, pos);
}
void MultiBodySystem::DeleteIOElemConNode(int elemnr, int list_idx)
{
	GetElement(elemnr).DeleteConNode(list_idx);
}
//$AD 2013-07-10: Change Position of a single element (MBS element, conNode, ...) 
void MultiBodySystem::MoveConNode2D(int elemnr, int list_idx, double delta_x, double delta_y) 
{
	GetElement(elemnr).MoveConNode2D(list_idx, delta_x, delta_y);
}
void MultiBodySystem::MoveElement(int elemnr, double delta_x, double delta_y, double delta_z)
{
	GetElement(elemnr).MoveElement(delta_x, delta_y, delta_z);
}

//action=1: option=0: Set initial conditions, option=1: Initialize
//action=2: Add Element: deprecated
//action=3: Assemble system
//action=4: Delete Element value
//action=5: Delete force value
//action=6: Add force: deprecated
//action=7: Add sensor:: deprecated
//action=8: Delete sensor value
//action=9: Add GeomElement: deprecated
//action=10:Delete GeomElement value
//action=11: Add Node: deprecated
//action=12: Delete Node value
//action=13: Swap element: deprecated
//action=14: Add Material: deprecated
//action=15: Delete Material value

//action=20: ModelFunction: option=1: get ModelFunctionList
//
//action=30: Call Element functions in response to action in IOBlockView Window... 

//action=101: save MBS system to file (file name in edc.Find("File_name")+directory, .GetText() (option=1: with Options, =0: only MBS)
//action=102: load MBS system from file (file name in edc.Find("File_name")+directory, .GetText()
//action=103: save actual solution vector (file name in edc.Find("File_name")+directory, .GetText()
//action=104: load initial solution vector (file name in edc.Find("File_name")+directory, .GetText()
//action=105: clear initial solution vector (setlen = 0 means initial solution is not used)
//action=106: LoadMBSEDCConfiguration
//action=107: SaveMBSEDCConfiguration; ; all:value==0, solver:value==1, hotint:value==2; 
//
//action=110: Load EDC from file (file name in edc.Find("File_name")+directory, .GetText()
//action=111: Save EDC to file (file name in edc.Find("File_name")+directory, .GetText()
//action=112: The input-edc gets written into a mystr, the mystr is returned in the same edc
//
//action=121: plot single sensor value in matlab
//action=122: plot two sensors (x-y-plot) in matlab
//action=123: plot n sensors (vs time) in matlab
//
//action=201: GetNSensors()
//action=202: Get Element Names in edc-data
//action=203: Get Sensor Names in edc-data
//            option 1: return writeresult flag as int element data
//action=204: Get Sensor value of sensor number 'value' in edc-data
//action=205: Get GeomElement Names in edc-data
//action=206: Get Number of nodes
//action=207: Get Name of nodes
//action=208: Get Number of total DOF
//action=209: Get Material Names in edc-data
//action=210: Get/Set EditAllOptions in edc-data Get:option==0, Set:option==1;    all:value==0, solver:value==1, hotint:value==2; 
//action=211: Get/Set Model ElementDataContainer
//action=212: SetOptions: Options2EDC/EDC2Options, etc.
//action=213: Get Sensors data array(reference), DataArray:option=0, TimesArray:option=1
//action=214: Get Load Names in edc-data	//$ DR 2012-10
//
//action=301: Show System properties
//action=302: Check System properties

//action=320: Compute Eigenmodes (option=1: sparse, =0: direct)

int MultiBodySystem::CallCompFunction(int action, int option, int value, ElementDataContainer* edc)
{
	Vector3D defaultbodycol(0.2,0.2,0.8);
	Vector3D defaultconstraintcol(0.2,0.8,0.2);
	Vector3D defaultbodypos(0.,0.,0.);
	Vector3D defaultbodyvel(0.,0.,0.);
	Vector3D defaultbodyangle(0.,0.,0.);
	Vector3D defaultbodyangularvel(0.,0.,0.);
	Vector3D defaultbodysize(1.,1.,1.);
	double defaultbodyrho = 1000;
	Vector3D defaultconstraintdim(0.1,0.1,16.);

	int defaultsensorvisible = 1;
	SensorSignalStorageMode defaultSignalStorageMode = (SensorSignalStorageMode)(SSM_OwnFile | SSM_CommonFile);
	int defaultsensorprecision = 17;
	Vector3D defaultsensordrawdim(0.001,6,0);


	//flexible bodies:
	double defaultbodyEm = 1e8;
	double defaultbodynu = 0.;

	int rv = 1;
	switch(action)
	{
	case 1: //Set initial conditions
		{
			if (option == 0)
				SetInitialConditions();
			else if (option == 1)
				Initialize(); //this resets the whole MBS model

			break;
		}
	case 2: //Add Element of type "value"
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;
		}
	case 3: //Assemble system
		{
			Assemble();
			break;
		}
	case 4: //Delete element "value"
		{
			if(value<=NE())
			{
				DeleteElement(value);
			}
			else
			{
				UO() << "Error: element " << value << " does not exist!!!\n";
			}
			break;
		}
	case 5: //Delete force "value" 
		{
			if(value<=NLoads())
			{
				DeleteLoad(value);
			}
			else
			{
				UO() << "Error: load " << value << " does not exist!!!\n";
			}
			break;
		}
	case 6: //Add Force type "value" 
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;
		}
	case 7: //Add Sensor type "value"
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;		}
	case 8: //Delete sensor "value"
		{
			if (value >= 1 && value <= this->NSensors())
			{
				DeleteSensor(value);
			}
			else
			{
				UO() << "Error: sensor " << value << " does not exist!!!\n";
			}
			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 9: //Add GeomElement SensorTypeNameList above!!!
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;		}
	case 10: //action: Delete GeomElement "value"
		{
			if (value >= 1 && value <= NDrawElements())
			{
				DeleteDrawElement(value);
			}
			else
			{
				UO() << "Error: geom element " << value << " does not exist!!!\n";
			}
			break;
		}
	case 11: //action: Add Node
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;
		}
	case 12: //action: Delete Node "value"
		{
			if (value >= 1 && value <= NNodes())
			{
				DeleteNode(value); //delte this node and reorder node references in elements!
			}
			else
			{
				UO() << "Error: Node " << value << " does not exist!!!\n";
			}
			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 13: //action: Swap element "value" with element "option"
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 14: //action: Add Material
		{
			UO(UO_LVL_warn)<<"Warning: CallCompFunction: a deprecated function is called!";
			break;		}
	case 15: //action: Delete Material "value"
		{
			if (value >= 1 && value <= NMaterials())
			{
				DeleteMaterial(value); //delete this material and reorder material numbers in elements!
			}
			else
			{
				UO() << "Error: Material " << value << " does not exist!!!\n";
			}
			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//action=20: ModelFunction: option=1: get ModelFunctionList
	case 20: //action: Get ModelFunctionList
		{
			if (option == 1) //get ModelFunctionList
			{
				ElementData ed;
				//$ YV 2013-01-04: model function index is now 1-based
				for (int i=1; i <= GetModelsLibrary()->GetModelsCount(); i++)
				{
					ed.SetText(GetModelsLibrary()->GetModelInterface(i)->GetMBSModelName());
					ed.SetToolTipText(GetModelsLibrary()->GetModelInterface(i)->GetMBSModelDescription());
					edc->Add(ed);
				}
				//GetIOption(12) = GetModelFunctionsIndex0();
			}
			else if (option == 2) //ModelFunctionChanged
			{
				if (GetModelDataContainer()) {GetModelDataContainer()->Reset();}		// DR hack
				ModelFunctionChanged();
			}
			else if (option == 3) //Initialize
			{
				Initialize();
			}
			else if (option == 4) //get modelfunctionindex
			{
				ElementData ed;
				ed.SetInt(GetModelFunctionsIndex0(),"model_function_index");
				edc->Add(ed);
			}
			break;
		}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 30: //action: Call Element functions in response to action in IOBlockView Window... 
		{
			int elnr = option;
			Element& elem = GetElement(elnr);

			if (value == 0) 
			{ 
				if (GetElement(elnr).GetElementSpecification().Compare(mystr("IOMouseResponseElement")) )
					rv = 1; 
				else 
					rv = 0; 
			}

			// allow modifications 1 -> add 1, and -1 -> subtract 1
			if (value == 1)	 
			{ 
				elem.Increase(); 
			}
			if (value == -1) 
			{ 
				elem.Decrease(); 
			}
			break;
		}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 31: // action: call MBS-wide response to a Keyboard input
		{
			int key = option;
			int nr_of_changes = RespondToKey(key);
			rv = nr_of_changes;
			break;
		}
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//load/save etc.
	case 101: //save to file
		{
			mystr filename, dirname;
			GetElemDataText(this, *edc, "File_name", filename);
			GetElemDataText(this, *edc, "Directory_name", dirname); //include backslash --> dirname+filename==filepath

			isloadsavemode = 1;
			SaveToFile(dirname+filename,option);
			isloadsavemode = 0;
			break;
		}
	case 102: //load from file
		{
			mystr filename, dirname;
			GetElemDataText(this, *edc, "File_name", filename);
			GetElemDataText(this, *edc, "Directory_name", dirname); //include backslash --> dirname+filename==filepath

			isloadsavemode = 1;
			rv = LoadFromFile(dirname+filename);
			isloadsavemode = 0;

			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 103: //save actual solution vector
		{
			mystr filename, dirname;
			GetElemDataText(this, *edc, "File_name", filename);
			GetElemDataText(this, *edc, "Directory_name", dirname); //include backslash --> dirname+filename==filepath

			SaveSolutionVector(dirname+filename);
			break;
		}
	case 104: //load initial solution vector
		{
			mystr filename, dirname;
			GetElemDataText(this, *edc, "File_name", filename);
			GetElemDataText(this, *edc, "Directory_name", dirname); //include backslash --> dirname+filename==filepath

			LoadInitialVector(dirname+filename);
			break;
		}
	case 105: //clear initial solution vector (setlen = 0 means initial solution is not used)
		{
			//UO() << "Stored initial vector cleared!\n";
			stored_initialconditions.SetLen(0);
			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 106: //Load_MBS_EDC_config_file
		{
			mystr filename = "hotint_cfg.txt";
			CMFile file(filename, TFMread);

			//set default options, if file can not be loaded!
			SetSolverDialogOptions();

			SetOptions2EDCOptions();//copy numbered i/d/t options to EDC Tree to 

			OpenLogFile(0);

			ElementDataContainer fileedc;

			if (file.IsGood()) 
			{
				UO() << "load configuration file 'hotint_cfg.txt'\n";

				mystr str;
				file.RWF(str); //read file

				// AD 2013-10-22 use new routine to parse
				EDCParser().ParseAndExecuteString(str,fileedc);
//				EDCParser().String2TreeEDC(/*this, */str, fileedc);

				//$ DR 2013-03-18 for release version, do not use stored model name
				if(!fileedc.TreeGetBool("GeneralOptions.Application.reload_last_model",0))
				{
					fileedc.TreeDelete("GeneralOptions.ModelFile.hotint_input_data_filename");	
					fileedc.TreeDelete("GeneralOptions.ModelFile.internal_model_function_name");	
					mystr tmp("");
					fileedc.TreeSetString("GeneralOptions.ModelFile.hotint_input_data_filename", tmp.c_str());
					tmp = mystr("new model");
					fileedc.TreeSetString("GeneralOptions.ModelFile.internal_model_function_name", tmp.c_str());
				}

				//replace default option parameters with loaded parameters
				GetMBS_EDC_Options()->TreeReplaceEDCDataWith(&fileedc);
			}
			else
			{
				UO() << "Warning: Could not read hotint configuration file 'hotint_cfg.txt'. Default values are used instead. Press 'Save Hotint Options' to create this file.\n";
				rv = 0; 
			}

			//at this point, use arguments of program start to add some functionality:
			mystr modeldata_argstr = "";
			mystr mbsedc_argstr = "";
			int narg = __argc;

			for (int i=0; i < narg; i++)
			{
				mystr argstr = __argv[i];
				UO(UO_LVL_ext) << "hotint.exe argument " << i << "='" << argstr << "'";

				//$ RL 2012-7-25:[ 
				if(i == 0) // i==0 --> store application path
				{
					// get application path, (path of hotint.exe)
					mystr app_path("");
					{
						int until;
						for(until=argstr.Length()-1;until>=0;until--)
						{
							if(argstr.PosPeek(until) == '\\')
								break; // until <=> last backslash
						}
						for(int i=0;i<=until;i++)
						{
							app_path += argstr.PosPeek(i); // path until last backslash
						}
						argstr = "GeneralOptions.Paths.application_path=\"" + app_path+ "\"";
					}
				}
				else if(i == 1) // i==1 --> from drag&drop or command line
				{
					mystr hid_file("");
					int pos = argstr.Find(mystr(".hid"));
					if(pos == -1)
					{
						pos = argstr.Find(mystr(".hmc"));		//$ DR 2012-12-10 this ending is also supported
					}
					if(pos == -1 && GetMBS_EDC_Options()->TreeGetBool("GeneralOptions.ModelFile.accept_txt_file_as_model_file",0))
					{
						pos = argstr.Find(mystr(".txt"));
					}
					if (pos != -1)
					{
						if(DoesFileExist(argstr) && argstr.Find(mystr("=")) == -1)
						{
							hid_file = argstr;
							argstr = "GeneralOptions.ModelFile.hotint_input_data_filename=\"" + argstr+ "\"";
							// get path of Hotint Input Data file
							mystr path("");
							{
								int until;
								for(until=hid_file.Length()-1;until>=0;until--)
								{
									if(hid_file.PosPeek(until) == '\\')
										break; // until <=> last backslash
								}
								for(int i=0;i<=until;i++)
								{
									path += hid_file.PosPeek(i); // path until last backslash
								}
								GetMBS_EDC_Options()->TreeSetString("GeneralOptions.Paths.hotint_input_data_path", path.c_str());
							}
						}
					}
				}	
				//$ RL 2012-7-25:] 

				int pos = argstr.Find('=');
				if (pos != -1)
				{				
					mystr name = argstr.Left(pos);
					//UO(UO_LVL_ext) << "name = '" << name << "'\n";
					name.EraseSpaces(); //name anyway should not contain spaces ...

					if (GetMBS_EDC_Options()->TreeFind(name))
					{
						UO(UO_LVL_ext) << " ... parsed to MBS_EDC_Option";
						mbsedc_argstr += argstr + mystr("\n");
					}
					else
					{
						UO(UO_LVL_ext) << " ... parsed to ModelData";
						modeldata_argstr += argstr + mystr("\n");
					}
				}
				UO() << "\n";
			}	// end of for loop narg

			if (mbsedc_argstr.Length() != 0)
			{
				UO() << "mbs-edc-argstring=" << mbsedc_argstr << "\n";
				// AD 2013-10-22 use new routine to parse
				EDCParser().ParseAndExecuteString(mbsedc_argstr,mbs_edc_args);
//				EDCParser().String2TreeEDC(/*this, */mbsedc_argstr, mbs_edc_args);
				GetMBS_EDC_Options()->TreeReplaceEDCDataWith(&mbs_edc_args);
			}
			if (modeldata_argstr.Length() != 0) //write information in to modeldata container ==> will not be overwritten later
			{
				// AD 2013-10-22 use new routine to parse
				EDCParser().ParseAndExecuteString(modeldata_argstr,modeldata_edc_args);
//				EDCParser().String2TreeEDC(/*this, */modeldata_argstr, modeldata_edc_args);

				if (GetModelDataContainer())
				{
					GetModelDataContainer()->TreeReplaceEDCDataWith(&modeldata_edc_args);
				}
				else
				{
					SetModelDataContainer(modeldata_edc_args);
				}
			}

			edc->CopyFrom(*GetMBS_EDC_Options()); //some options are written in configuration

			SetEDCOptions2Options(); //copy EDC Tree to numbered i/d/t options for fast access
			SetComputationSolverOptions(); //copy EDC Tree to MBSSolSet()

			OpenLogFile(0);

			break;
		}
	case 107: //Save_MBS_EDC_config_file
		{
//!AD 2012-06-26: revised
//!AD 2012-04-17: solver options and hotint options can be stored separately - use autogenerated functions, changed names to better distinguish various EDCs
			mystr filename = "hotint_cfg.txt";
			CMFile file_read(filename, TFMread);
	
// read from current configuration file;
			ElementDataContainer edc_configfile;
			if (file_read.IsGood()) 
			{
			  mystr str;
				file_read.RWF(str);
				// AD 2013-10-22 use new routine to parse
				EDCParser().ParseAndExecuteString(str, edc_configfile);
//				EDCParser().String2TreeEDC(/*this, */str, edc_configfile);
			}
			else
			{
				UO() << "Error: could not read hotint configuration file 'hotint_cfg.txt'\n";
			}
		//	delete &file_read;


////// actualize solver option part of the options edc //AD: revision
////			SetComputationSolverOptions(); // <-- wrong
////			SetOptions2EDCOptions();
			SetSolverDialogOptions();
			SetOptions2EDCOptions();

// pick the data that will be written to file
			ElementDataContainer edc_tosave;
			if(value == 1) // save solver options only
			{
				SolverOptions2EDC(&edc_tosave);
				UO() << "save configuration file 'hotint_cfg.txt': solver options only\n";
			}
			else if(value == 2) // save hotint options only
			{
			  Options2EDC(&edc_tosave);
				UO() << "save configuration file 'hotint_cfg.txt': hotint options only\n";
			}
			else // rest of old code...
			{
				edc_tosave = *GetMBS_EDC_Options();
				UO() << "save configuration file 'hotint_cfg.txt': all options\n";
			}
// replace in the content of the config file
			edc_configfile.TreeReplaceEDCDataWith(&edc_tosave);

			//$ AD 2011-09: do not write ComputationSteps to the config file
			RemoveVariableFromEDCTree(edc_configfile, mystr("SolverOptions.ComputationSteps"), 0);

// edc contains window placement etc...
			if( value != 1 )
			{
				edc_configfile.TreeReplaceEDCDataWith(edc);
			}

// write to configfile
			CMFile file_write(filename, TFMwrite);

			//UO() << "load arrow=" << edcnew.TreeGetDouble("PostProcOptions.Loads.arrow_size") << "\n";

			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//set some default values, which should not be loaded:
			edc_configfile.TreeSetDouble("SolverOptions.start_time", 0.); //do not load any starttime except 0.
			edc_configfile.TreeSetBool("ViewingOptions.Animation.animate_deformation", 0); //animate at startup
			edc_configfile.TreeSetBool("ViewingOptions.Animation.RecordSingleFrames.record", 0); //do not record at beginning
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


			mystr str;
			EDCParser().EDC2String(edc_configfile, str);

			if(UO().GetGlobalMessageLevel()==UO_LVL_dbg2)	//$ DR 2011-09-15
			{
				UO(UO_LVL_dbg2) << str << "\n";
			}
			file_write.RW(mystr("%HOTINT configuration file 'hotint_cfg.txt'\n"));
			file_write.RW(mystr("%handle this file with care, it can crash the startup of HOTINT\n"));
			file_write.RW(mystr("%in case that nothing works anymore ==> delete this file ==> defaults will be created ==> save configuration ==> new 'hotint_cfg.txt' will be created\n"));
			file_write.RW(mystr("%\n"));
			file_write.RW(str);

			break;
		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//Load and Save ANY EDC
//action=110: Load EDC from file (file name in edc.Find("File_name")+directory, .GetText()
	case 110:
		{
			mystr filename = edc->TreeGetString("File_name");
			mystr the_string;
			edc->Reset();

			File2Str(filename, the_string);
			EDCParser().ParseAndExecuteString(the_string, *edc);
			break;

			//the_file.RWF(the_string);
			//ElementDataContainer the_edc,testedc;
			////EDCParser().String2EDC(the_string, the_edc, testedc);

			//EDCParser().String2TreeEDC(the_string, the_edc);
			//edc->CopyFrom(the_edc);
			//break;
		}
//action=111: Save EDC to file (file name in edc.Find("File_name")+directory, .GetText()
	case 111:
		{
			mystr filename = edc->TreeGetString("File_name");
			edc->Delete(edc->Find("File_name"));
			CMFile the_file(filename, TFMwrite);
			mystr the_string;

			EDCParser().EDC2String(*edc, the_string);
			the_file.RW(the_string);
	    break;
		}
//action=112: The input-edc gets written into a mystr, the mystr is returned in the same edc
	case 112:
		{
			mystr clipboardtext;
			EDCParser().EDC2String(*edc, clipboardtext, mystr("  "));

			//edc = new ElementDataContainer;
			ElementData ed;
			ed.SetText(clipboardtext, "ClipBoardText");
			edc->Add(ed);
			break;
		}
	case 113: //$ DR 2013-06-12
		{
			mystr str = edc->TreeGetString("File_name");
			int fromfile = edc->TreeGetInt("is_file");
			AddModelData(str/*,fromfile*/);
			break;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//plot in matlab
	case 121: //plot sensor 'value' in matlab
		{
			if (GetSensor(value).GetSignalStorageMode() & SSM_OwnFile == 0) return 0;

			UO().CallWCDriverFunction(4); //Set Programm directory as current

			mystr dir = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path"); //GetTOption80);
			//mystr filename = mystr("S")+mystr(GetSensor(value).GetSensorNumber())+mystr("-")+GetSensor(value).GetTypeName()+mystr(".txt");
			mystr filename = mystr("S")+mystr(value)+mystr("-")+GetSensor(value).GetSensorName()+mystr(".txt");
			mystr dir2 = dir.Left(dir.Length()-1); //without the '/' sign

			ofstream plotfunc((dir+mystr("MBSplotfunc.m")).c_str());

			char buffer[_MAX_PATH*4];

			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr oldpath = buffer;

			_chdir(dir.c_str());
			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr path = buffer;

			char str[1024]; 
			//if (path[path.Length()-1] != '\\') path += '\\';

			sprintf(str, "f1=load('%s');\n",(filename).c_str());
			plotfunc << str;
			plotfunc << "h=plot(f1(:,1),f1(:,2),'k');\n";

			plotfunc << "title(' Y: " << GetSensor(value).GetSensorName() << "','Fontsize',16);\n"; 
			plotfunc << "xlabel(' Time [s] ','Fontsize',16);\n";
			//$ YV 2012-06: there is no type name for the sensors
			//plotfunc << "ylabel('" << GetSensor(value).GetTypeName() << "','Fontsize',16);\n";
			plotfunc << "set(h,'LineWidth',1.0);\n";
			//plotfunc << "legend('" << GetSensor(value).GetTypeName() << "');\n";
			plotfunc << "waitfor(h);\n"; //wait till window is closed
			plotfunc << "exit;\n";       //exit instance of matlab; too many instances will need enormous memory and may crash system 
			plotfunc.close();
			ShellExecute(NULL, "open", "matlab",	
				" -nosplash /r MBSplotfunc",
				path.c_str(), SW_HIDE);
			//path.c_str(), SW_SHOW);

			_chdir(oldpath.c_str());

			break;
		}
	case 122: //plot two sensors 'edc.xsensor' vs 'edc.ysensor'  in matlab
		{
			int xsensor = edc->TreeGetInt("xsensor");
			int ysensor = edc->TreeGetInt("ysensor");

			if (GetSensor(xsensor).GetSignalStorageMode() & SSM_OwnFile == 0) return 0;
			if (GetSensor(ysensor).GetSignalStorageMode() & SSM_OwnFile == 0) return 0;

			UO().CallWCDriverFunction(4); //Set Programm directory as current

			mystr dir = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path"); //GetTOption80);
			mystr dir2 = dir.Left(dir.Length()-1); //without the '/' sign
			mystr filename = mystr("S")+mystr(xsensor)+mystr("-")+GetSensor(xsensor).GetSensorName()+mystr(".txt");
			mystr filename2 = mystr("S")+mystr(ysensor)+mystr("-")+GetSensor(ysensor).GetSensorName()+mystr(".txt");

			ofstream plotfunc((dir+mystr("MBSplotfunc.m")).c_str());

			char buffer[_MAX_PATH*4];

			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr oldpath = buffer;

			_chdir(dir.c_str());
			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr path = buffer;

			char str[1024]; 

			sprintf(str, "f1=load('%s');\n",(filename).c_str());
			plotfunc << str;
			sprintf(str, "f2=load('%s');\n",(filename2).c_str());
			plotfunc << str;
			plotfunc << "h=plot(f1(:,2),f2(:,2),'k');\n";
			plotfunc << "title({[' X: " << GetSensor(xsensor).GetSensorName() << " '];[' Y: " << GetSensor(ysensor).GetSensorName() << " ']},'Fontsize',16);\n"; 
			//$ YV 2012-06: there is no type name for the sensors
			//plotfunc << "xlabel('" << GetSensor(xsensor).GetTypeName() << "','Fontsize',16);\n";
			//plotfunc << "ylabel('" << GetSensor(ysensor).GetTypeName() << "','Fontsize',16);\n";
			plotfunc << "xlabel('" << GetSensor(xsensor).GetSensorName() << "','Fontsize',16);\n";
			plotfunc << "ylabel('" << GetSensor(ysensor).GetSensorName() << "','Fontsize',16);\n";
			plotfunc << "set(h,'LineWidth',1.0);\n";
			plotfunc << "waitfor(h);\n"; //wait till window is closed
			plotfunc << "exit;\n";       //exit instance of matlab; too many instances will need enormous memory and may crash system 
			plotfunc.close();
			ShellExecute(NULL, "open", "matlab",	
				" -nosplash /r MBSplotfunc",
				path.c_str(), SW_HIDE);
			//path.c_str(), SW_SHOW);

			_chdir(oldpath.c_str());

			break;

		}
	case 123: //time plot of N sensors in matlab
		{
			if(option <=0) 
			{
				return 0; //no sensors selected
			}
			TArray<int> ysensor;
			ysensor.SetLen(option);

			for(int sensNr=1;sensNr<=option;sensNr++)
			{
				ysensor(sensNr) = edc->TreeGetInt("y"+mystr(sensNr)+"sensor");
				if (GetSensor(ysensor(sensNr)).GetSignalStorageMode() & SSM_OwnFile == 0) return 0;
			}
			UO().CallWCDriverFunction(4); //Set Programm directory as current

			mystr dir = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path"); //GetTOption80);
			mystr dir2 = dir.Left(dir.Length()-1); //without the '/' sign


			ofstream plotfunc((dir+mystr("MBSplotfunc.m")).c_str());

			char buffer[_MAX_PATH*4];

			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr oldpath = buffer;

			_chdir(dir.c_str());
			_getcwd( buffer, _MAX_PATH*4); //about 1000 characters ...
			mystr path = buffer;

			char str[1024]; 
			plotfunc << "scrsz = get(0,'ScreenSize');\n";
			plotfunc << "fig = figure;\n";
			plotfunc << "set(fig,'Position',[1 (scrsz(4)/2-60) scrsz(3)*2/3 scrsz(4)/2]);\n";
			for(int sensNr=1;sensNr<=option;sensNr++)
			{
				mystr filename = mystr("S")+mystr(ysensor(sensNr))+mystr("-")+GetSensor(ysensor(sensNr)).GetSensorName()+mystr(".txt");
				sprintf(str, mystr("f") + mystr(sensNr) + +mystr("=load('%s');\n"),(filename).c_str());				
				plotfunc << str;
				plotfunc << "h(" << sensNr << ")=plot(f" << sensNr << "(:,1),f" << sensNr << "(:,2),";
				if(sensNr==1){plotfunc << "'k'";}
				else if(sensNr==1){plotfunc << "'k'";}
				else if(sensNr==2){plotfunc << "'b'";}			
				else if(sensNr==3){plotfunc << "'g'";}		
				else if(sensNr==4){plotfunc << "'r'";}
				else if(sensNr==5){plotfunc << "'c'";}
				else if(sensNr==6){plotfunc << "'m'";}
				else if(sensNr==7){plotfunc << "'k--'";}
				else if(sensNr==8){plotfunc << "'b--'";}
				else {plotfunc << "'k'";}
				plotfunc << ",'linewidth',1);\n";
				plotfunc << "hold on;\n";
			}
			//plotfunc << "title({[' Y1: " << GetSensor(ysensor(1)).GetSensorName() << " ']";
			//for(int sensNr=2;sensNr<=option;sensNr++)
			//{
			//	plotfunc << ";[' Y" << sensNr << ": " << GetSensor(ysensor(sensNr)).GetSensorName() << " ']";
			//}
			//plotfunc << "},'Fontsize',16);\n";
			int sameSensorType = 1;
			for(int sensNr=2;sensNr<=option;sensNr++)
			{
				if(!GetSensor(ysensor(sensNr-1)).GetTypeName().Compare(GetSensor(ysensor(sensNr)).GetTypeName()))
				{				
					sameSensorType = 0;
					break;
				}
			}
			if(sameSensorType)
			{
				plotfunc << "ylabel('" << GetSensor(ysensor(1)).GetTypeName() << "','Fontsize',16);\n";// same type
			}
			else
			{				
				plotfunc << "ylabel({['" << GetSensor(ysensor(1)).GetTypeName() << "']";
				for(int sensNr=2;sensNr<=option;sensNr++)
				{
					plotfunc << ";['" << GetSensor(ysensor(sensNr)).GetTypeName() << "']";	// different types
				}
				plotfunc << "},'Fontsize',16);\n"; //$ MS 2011-4-5
			}

			plotfunc << "xlabel('time','Fontsize',16);\n";
			plotfunc << "legend('" << GetSensor(ysensor(1)).GetSensorName();
			for(int sensNr=2;sensNr<=option;sensNr++)
			{
				plotfunc	<< "','" << GetSensor(ysensor(sensNr)).GetSensorName();
			}
			plotfunc << "');\n";
			//plotfunc << "set(h,'LineWidth',1.0);\n";
			plotfunc << "waitfor(h(1));\n"; //wait till window is closed
			plotfunc << "exit;\n";       //exit instance of matlab; too many instances will need enormous memory and may crash system 
			plotfunc.close();
			ShellExecute(NULL, "open", "matlab",	
				" -nosplash /r MBSplotfunc",
				dir.c_str(), SW_HIDE);
			//path.c_str(), SW_HIDE);
			//path.c_str(), SW_SHOW);

			_chdir(oldpath.c_str());

			break;

		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//information exchange, etc.
	case 201: //NSensors
		{
			return NSensors();
			break;
		}
	case 202: //Element names
		{
			for (int i=1; i <= NE(); i++)
			{
				int flag = 0; //add information on constraint, rigid, flexible, actuator, spring?
				ElementData ed;
				ed.SetInt(flag, GetElement(i).GetElementName().c_str()); edc->Add(ed);
			}
			break;
		}
	case 203: //Sensor names
		{
			for (int i=1; i <= NSensors(); i++)
			{
				int flag = 0; //add information on type?
				if(value==1)
				{
					flag = GetSensor(i).GetSignalStorageMode() != SSM_None; // write result data
				}
				ElementData ed;
				//$ YV 2012-06: display names do not seem to be actually used
				/*
				if(option ==1) // return the displaynames (plottool caption)
				{	
					ed.SetInt(flag, GetSensor(i).GetDisplayName().c_str()); edc->Add(ed);
				}
				else //if(option==0) // return sensornames
				{
					ed.SetInt(flag, GetSensor(i).GetSensorName().c_str()); edc->Add(ed);
				}
				*/
				ed.SetInt(flag, GetSensor(i).GetSensorName().c_str()); edc->Add(ed);
			}
			break;
		}		
	case 204: //Sensor value 
		{
			ElementData ed;
			if(option==1) // return time & all system values
			{
				Vector v(1+NSensors());
				v[0] = GetTime();
				for(int i=1; i <= NSensors(); i++)
				{
					v[i] = GetSensor(i).GetLastValue(); 
				}
				ed.SetVector(v.GetVecPtr(),v.GetLen(),"values");
				edc->Add(ed);
			}
			else //if(option==0) // return single value from sensor #value
			{
				double v = GetSensor(value).GetLastValue();
				ed.SetDouble(v, "value"); edc->Add(ed);
			}
			break;
		}
	case 205: //GeomElement names //MBS: DrawElement == GeomElement
		{
			for (int i=1; i <= NDrawElements(); i++)
			{
				int flag = 0; //add information on type?
				ElementData ed;
				ed.SetInt(flag, GetDrawElement(i)->GetName().c_str()); edc->Add(ed);
			}
			break;
		}
	case 206: //NNodes
		{
			return NNodes();
			break;
		}
	case 207: //Nodes names
		{
			for (int i=1; i <= NNodes(); i++)
			{
				int flag = 0; //add information on type?
				ElementData ed;
				ed.SetInt(flag, (GetNode(i).GetObjectName()+mystr("-")+mystr(GetNode(i).SOS())+("DOF")).c_str()); edc->Add(ed);
			}
			break;
		}
	case 208: //N initial values
		{
			int ss = GetSystemSize();
			return ss;
			break;
		}
	case 209: //Material names
		{
			for (int i=1; i <= NMaterials(); i++)
			{
				int flag = 0; //add information on type?
				ElementData ed;
				//ed.SetInt(flag, (GetMaterial(i).GetMaterialName()+mystr("-")+mystr(GetMaterial(i).GetTypeSpec())).c_str()); edc->Add(ed);
				ed.SetInt(flag, (GetMaterial(i).GetMaterialName()).c_str()); edc->Add(ed);	//$ DR 2013-01-21 change for script language
			}
			break;
		}
	case 210: //EditAllOptions
		//!AD: revised 2012-07-26
		{
			if (option==0) //Get all mbs_edc_options and mbssolset into edc
			{
				edc->Reset();

				if(value == 1) // save solver options only
				{
					SolverOptions2EDC(edc);
				}
				else if(value == 2) // save hotint options only
				{
					Options2EDC(edc);
				}
				else if(value == 0) // entire edc - currently no dialog active for all options
				{
					SetOptions2EDCOptions();//copy numbered i/d/t options to EDC Tree options ==> this is because old dialogs modify the i/d/toptions //JG
					edc->CopyFrom(*GetMBS_EDC_Options());
				}

				////--------------------------------------------
				////RL:b
				//ElementDataContainer edcnew;
				////copy mbssolset to edc:
				//Options2EDC(&edcnew);
				////replace edcnew data in edc:
				//edc->TreeReplaceEDCDataWith(&edcnew);
				////RL:e
				////--------------------------------------------
			}
			else if (option==1) //Write modified edc-data to mbs_edc_options
			{
				GetMBS_EDC_Options()->TreeReplaceEDCDataWith(edc); // overwrite MBS-EDC first

				if(value == 1) // update solver options only
				{
					SetComputationSolverOptions();
				}
				else if(value == 2) // update hotint options only
				{
					SetEDCOptions2Options();
				}
				else if(value == 0) // update entire edc - currently no dialog active for all options
				{
					SetComputationSolverOptions();
					SetEDCOptions2Options();
				}

				//==>this must always be done, otherwise the EDC options are overwritten by other functions!!!!
//<<<<<<< mbs_communication.cpp
//
//				//$ JG: instantly change log values:
//				MBSSolSet().nls_log_rel_error_goal = edc->TreeGetInt("SolverOptions.Newton.Log.rel_error_goal");
//				MBSSolSet().nls_log_jacobian_update_info = edc->TreeGetInt("SolverOptions.Newton.Log.jacobian_update_info");
//				MBSSolSet().nls_log_contractivity = edc->TreeGetInt("SolverOptions.Newton.Log.contractivity");
//				MBSSolSet().nls_log_iteration_error = edc->TreeGetInt("SolverOptions.Newton.Log.iteration_error");
//				MBSSolSet().nls_log_solution_vector = edc->TreeGetInt("SolverOptions.Newton.Log.solution_vector");
//				MBSSolSet().nls_log_jacobian = edc->TreeGetInt("SolverOptions.Newton.Log.jacobian");
//
//				MBSSolSet().output_level = edc->TreeGetInt("SolverOptions.Log.output_level");
//				MBSSolSet().outprec_double = edc->TreeGetInt("SolverOptions.Log.output_precision_double");
//				MBSSolSet().outprec_vector = edc->TreeGetInt("SolverOptions.Log.output_precision_vector");
//				MBSSolSet().outprec_matrix = edc->TreeGetInt("SolverOptions.Log.output_precision_matrix");
//				MBSSolSet().max_errors = edc->TreeGetInt("SolverOptions.Log.max_error_messages");
//				MBSSolSet().max_warnings = edc->TreeGetInt("SolverOptions.Log.max_warning_messages");
//				MBSSolSet().log_detailed_comp_output = edc->TreeGetInt("SolverOptions.Log.detailed_computation_output");
//				MBSSolSet().comp_output_every_x_sec = edc->TreeGetDouble("SolverOptions.Log.computation_output_every_x_sec");
//				MBSSolSet().write_MK2file = edc->TreeGetInt("SolverOptions.Log.write_mass_and_stiffness_matrix");
//
//				MBSSolSet().withgraphics = edc->TreeGetDouble("SolverOptions.Misc.with_graphics");
//=======

				OpenLogFile(1);
			}

			break;
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 211: //Model Data
		{
			if (option==0 && GetModelDataContainer() != 0) //Get Data in EDC
			{
				edc->Reset();
				edc->CopyFrom(*GetModelDataContainer());
			}
			else if (option==1) //Set Data in EDC
			{
				if (!GetModelDataContainer()) //generate empty modeldatacontainer
				{
					ElementDataContainer edcnew;
					SetModelDataContainer(edcnew);
				}
				else
				{
					GetModelDataContainer()->Reset();
				}
				GetModelDataContainer()->CopyFrom(*edc);
			}
			else if (option==2) //Load Model Data
			{
			}
			else if (option==3) //Save Model Data
			{
				ElementDataContainer* medc = GetModelDataContainer();
				mystr filestr;
				EDCParser().EDC2String(*medc, filestr);

				{
					mystr filename = GetTOption( 88); //directory+filename
					int opt = ios::out;


					ofstream* storeout = new ofstream;
#ifdef my_new_stdiostream
					storeout->open(filename.c_str(), opt); //with read sharing allowed!
#else
					storeout->open(filename.c_str(), opt, filebuf::sh_read); //with read sharing allowed!
#endif

					if ((int)storeout->good() == 0 || storeout->fail())  
					{
						UO() << "ERROR: could not write model data!!!\n";
						UO() << "good=" << storeout->good() << ", fail=" << storeout->fail() << "\n";
						UO() << "file=" << filename.c_str() << "\n";
					}
					else
					{
						(*storeout) << filestr.c_str() << "\n";
					}

					storeout->close();
					delete storeout;

				}
			}

			break;
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 212: //SetOptions
		{
			if (option==1) 
			{
				SetOptions2EDCOptions();//copy numbered i/d/t options to EDC Tree options
			}
			else if (option==2) 
			{
				SetEDCOptions2Options();//copy EDC Tree options to numbered i/d/t options
			}
			break;
		}
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 214: //Load names
		{
			for (int i=1; i <= NLoads(); i++)
			{
				int flag = 0; //add information on type?
				ElementData ed;
				ed.SetInt(flag, GetLoad(i).LoadName().c_str()); edc->Add(ed);
			}
			break;
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 301: //Compute system information
		{
			int rs = GetResortSize();
			int sos = GetSecondOrderSize();
			//int sos_rs = GetSecondOrderSize_RS();
			int es = GetFirstOrderSize();
			int is = GetImplicitSize();
			int ss = GetSystemSize();

			int nbodies=0;
			int njoints=0;
			for (int i=1; i <= NE(); i++)
			{
				if (GetElement(i).IsType(TBody)) nbodies++;
				else njoints++;
			}

			UO() << "++++++++++++++++++++++++++++++++++\n";
			UO() << "System properties:\n";
			UO() << "2nd order ODE=" << sos << " (sos)\n";
			UO() << "Algebraic Equ=" << is << " (is)\n";
			UO() << "1st order ODE=" << es << " (es)\n";
			UO() << "number of initial values=" << ss << " (=2*sos+es+is)\n";
			UO() << "Resorted equ =" << rs << "\n";
			UO() << "N. of bodies=" << nbodies << "\n";
			UO() << "N. of joints=" << njoints << "\n";

			UO() << "++++++++++++++++++++++++++++++++++\n";
			break;
		}
	case 302: //check system; rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
		{
			rv = CheckSystemConsistency();
			break;
		}

	case 320: 
		{
			rv = ComputeEigenmodes();
			break;
		}
	default: ;
	}

	return rv;
}

//$ DR 2012-11-29 new code according to skript-language
void MultiBodySystem::SaveToFile(mystr filename,int save_options)
{
	mystr el = "\n"; //end of line
	mystr ob = "{";  //open bracket
	mystr cb = "}";  //closing bracket
	mystr ds = ": "; //double stop
	mystr eq = "= "; //equal sign

	ofstream of(filename.c_str());

	//++++++++++++++++
	//header:
	char str[32];
	/*sprintf(str, "%0.3f\n\n", HOTINT_version_info);*/
	
	_SYSTEMTIME time;
	GetSystemTime(&time);

	// do not use HOTINT_data_file_version in header
	mystr header = "";

	header += mystr("%======================================================================\n");
	header += mystr("%                      HOTINT model \n");
	header += mystr("%\n");
	header += mystr("% model name: ") + mystr(MBS_EDC_TreeGetString("GeneralOptions.ModelFile.internal_model_function_name")) + mystr("\n");
	header += mystr("%\n");
	header += mystr("% saved on: ")+mystr(time.wDay)+mystr(".")+mystr(time.wMonth)+mystr(".")+mystr(time.wYear)+mystr(" at ")+mystr(time.wHour)+mystr(":")+mystr(time.wMinute)+mystr("\n");
	header += mystr("%======================================================================\n");
	header +=mystr("HOTINT_data_file_version = \"")+mystr(GetHotintVersion().GetString())+mystr("\" % version of HOTINT, that has been used to generate this file \n\n");

	of << header;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//		Check if the mbs can be saved at all
	bool error = 0;

	if(error)
	{
		UO() << "Writing of output file aborted due to the error(s) above! \n";
		of << "% Writing of output file aborted, because nodes are not available for skript language.";
		of.close();
		return;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//materials:
	for (int i=1; i<=NMaterials(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		Material& m = GetMaterial(i);
		m.GetElementData(edc);

		otext += "Material" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket

		if(m.GetType() & TMatBeam)
		{
			otext += "AddBeamProperties(Material" +mystr(i)+ ")" + el +el;	
		}
		else
		{
			otext += "AddMaterial(Material" +mystr(i)+ ")" + el +el;	
		}

		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// nodes:
	for (int i=1; i<=NNodes(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		Node& m = GetNode(i);
		m.GetElementData(edc);

		otext += "Node" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket
		otext += "AddNode(Node" +mystr(i)+ ")" + el +el;	
		
		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//elements:
	for (int i=1; i<=NE(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		Element& e = GetElement(i);
		e.GetElementData(edc);

		otext += "Element" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket
		
		if (e.IsType(TConstraint) && !e.IsType(TController))
		{
			otext += "AddConnector(Element" +mystr(i)+ ")" + el +el;	// add the element
		}
		else
		{
			otext += "AddElement(Element" +mystr(i)+ ")" + el +el;		// add the element
		}
		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//loads:
	for (int i=1; i<=NLoads(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		MBSLoad& l = GetLoad(i);
		l.GetElementData(edc);

		otext += "Load" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket
		otext += "AddLoad(Load" +mystr(i)+ ")" + el +el;	// add the element
		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//sensors:
	for (int i=1; i<=NSensors(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		Sensor& s = GetSensor(i);
		s.GetElementData(edc);

		otext += "Sensor" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket
		otext += "AddSensor(Sensor" +mystr(i)+ ")" + el +el;	// add the element
		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// GeomElements:
	for (int i=1; i<=NDrawElements(); i++) 
	{
		mystr otext = "";
		ElementDataContainer edc;
		GeomElement *e = GetDrawElement(i);
		e->GetElementData(edc);

		otext += "GeomElement" + mystr(i);
		otext += el + ob + el;			// open bracket

		mystr str;
		EDCParser().EDC2String(edc, str, "  ");
		otext += str;								// element data
		otext += cb + el;						// close bracket
		
		otext += "AddGeomElement(GeomElement" +mystr(i)+ ")" + el +el;		// add the element

		of << otext;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//options:
	if(save_options)
	{
		ElementDataContainer edc;

		//SolverOptions2EDC(edc); // save solver options only
		//Options2EDC(edc); // save hotint options only

		// entire edc - currently no dialog active for all options
		SetOptions2EDCOptions();//copy numbered i/d/t options to EDC Tree options ==> this is because old dialogs modify the i/d/toptions //JG
		edc.CopyFrom(*GetMBS_EDC_Options());

		////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		////set some default values, which should not be loaded:
		//edc_configfile.TreeSetDouble("SolverOptions.start_time", 0.); //do not load any starttime except 0.
		//edc_configfile.TreeSetBool("ViewingOptions.Animation.animate_deformation", 0); //animate at startup
		//edc_configfile.TreeSetBool("ViewingOptions.Animation.RecordSingleFrames.record", 0); //do not record at beginning
		////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		edc.TreeDelete("GeneralOptions.ModelFile.hotint_input_data_filename");	
		edc.TreeDelete("GeneralOptions.ModelFile.internal_model_function_name");
		edc.TreeDelete("GeneralOptions.ModelFile.recent_file1");
		edc.TreeDelete("GeneralOptions.ModelFile.recent_file2");
		edc.TreeDelete("GeneralOptions.ModelFile.recent_file3");
		edc.TreeDelete("GeneralOptions.ModelFile.recent_file4");
		edc.TreeDelete("GeneralOptions.ModelFile.recent_file5");

		mystr opts;
		EDCParser().EDC2String(edc, opts, "");
		of << "%======================================================================\n";
		of << "%                      Solver Options \n";
		of << opts;	// write the string to file
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//the end 
	of.close();

	UO(UO_LVL_all) << "\nMBS system saved as '" << filename << "'\n";
}

int MultiBodySystem::LoadFromFile(mystr filename)
{
	Destroy();
	int rv = ReadModelData(filename);
	if(rv)
	{
		SetModelData_Initialized(1);
		//Assemble(); // this is done by ModelChanged() wich is called later anyway
	}
	else
	{
		UO() << "ERROR in LoadFromFile during reading '" << filename << "'\n";
	}
	return rv;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//these functions are called from WCDriver via WinCompInterface:
MBSObjectFactoryInterface * MultiBodySystem::GetObjectFactory()
{
	return EDCParser().GetObjectFactory();
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++               EDC_MBS_OPTIONS             ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ElementDataContainer* MultiBodySystem::GetMBS_EDC_Options() 
{
	return mbs_edc_options;
}
const ElementDataContainer* MultiBodySystem::GetMBS_EDC_Options() const
{
	return mbs_edc_options;
}
void MultiBodySystem::SetMBS_EDC_Options(const ElementDataContainer& edc)
{
	if (mbs_edc_options != 0) delete mbs_edc_options;

	mbs_edc_options = new ElementDataContainer;

	mbs_edc_options->CopyFrom(edc);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ElementDataContainer* MultiBodySystem::GetMBS_EDC_Variables() 
{
	return mbs_edc_global_variables;
}
const ElementDataContainer* MultiBodySystem::GetMBS_EDC_Variables() const
{
	return mbs_edc_global_variables;
}
void MultiBodySystem::SetMBS_EDC_Variables(const ElementDataContainer& edc)
{
	if (mbs_edc_global_variables != 0) delete mbs_edc_global_variables;

	mbs_edc_global_variables = new ElementDataContainer;

	mbs_edc_global_variables->CopyFrom(edc);
}

//get values in MBS_EDC:
void MultiBodySystem::MBS_EDC_TreeGetDouble(double& data, const char* name) const
{
	data = GetMBS_EDC_Options()->TreeGetDouble(name);
}
void MultiBodySystem::MBS_EDC_TreeGetInt(int& data, const char* name) const
{
	data = GetMBS_EDC_Options()->TreeGetInt(name); //(int == bool internally)
}
const char* MultiBodySystem::MBS_EDC_TreeGetString(const char* name) const
{
	return GetMBS_EDC_Options()->TreeGetString(name);
}

//set values in MBS_EDC:
void MultiBodySystem::MBS_EDC_TreeSetDouble(double data, const char* name)
{
	GetMBS_EDC_Options()->TreeSetDouble(name, data);
}
void MultiBodySystem::MBS_EDC_TreeSetInt(int data, const char* name)
{
	GetMBS_EDC_Options()->TreeSetInt(name, data); //(int == bool internally)
}
void MultiBodySystem::MBS_EDC_TreeSetString(const char* data, const char* name)
{
	GetMBS_EDC_Options()->TreeSetString(name, data);
}

const CMBSParser& MultiBodySystem::MBSParser() const {return *mbs_parser;}
CMBSParser& MultiBodySystem::MBSParser() {return *mbs_parser;}

const CEDCParser& MultiBodySystem::EDCParser() const {return *edc_parser;}
CEDCParser& MultiBodySystem::EDCParser() {return *edc_parser;}

int MultiBodySystem::File2Str(const char* filename,	mystr& str)
{
	mystr filen = filename;
	CMFile file(filename, TFMread);
	UO()	<< "--------------------------------\n";

	//read parameters from file
	if (file.IsGood()) 
	{
		UO() << "read " << filename << "\n";

		//read file
		file.RWF(str);	
	}
	else
	{
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "ERROR: failed to read model data from " << filename << "\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";

		return 0;  // error
	}
	return 1;
}

//this function reads a file and stores it in an edc
int MultiBodySystem::File2EDC(const char* filename, ElementDataContainer* edc_file)
{
	//$ DR 2013-01-14:[ added call for stl-files
	mystr filen = filename;
	if(filen.Find(".stl")>0)
	{
		return STLFile2EDC(filen, 0, edc_file);
	}
	//$ DR 2013-01-14:]

	CMFile file(filename, TFMread);
	UO()	<< "--------------------------------\n";

	//read parameters from file
	if (file.IsGood()) 
	{
		UO() << "read " << filename << "\n";

		//read file
		mystr str;
		file.RWF(str);		

		// overwrite edc with contents from file
//		EDCParser().String2TreeEDC(/*this, */str, *edc_file);

		EDCParser().ParseAndExecuteString(str, *edc_file);
	}
	else
	{
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "ERROR: failed to read model data from " << filename << "\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";

		return 0;  // error
	}
	return 1;
}

//DR 2013-01-14 this function reads a stl-file and writes the data in the edc
int MultiBodySystem::STLFile2EDC(char* file, int binary, ElementDataContainer* return_value)
{
	int initsize = 10000;
	TArray<Vector3D> newnormals(initsize);
	TArray<Vector3D> newpoints(3*initsize);
	double stl_tolerance=1e-4;								//tolerance for equal points when reading

	if (!binary)
	{
		ifstream ist(file);

		char buf[200];
		Vector3D pts[3];
		Vector3D normal;

		int cntface = 0;
		int vertex = 0;
		int badnormals = 0;

		while (ist.good())
		{
			ist >> buf;
			if (strcmp (buf, "facet") == 0)
			{
				cntface++;
			}

			normal = 0;

			if (strcmp (buf, "normal") == 0)
			{
				ist >> normal.X() 
					>> normal.Y()
					>> normal.Z();
				if (normal.Norm() > 0)
					normal /= normal.Norm();
			}


			if (strcmp (buf, "vertex") == 0)
			{
				ist >> pts[vertex].X()
					>> pts[vertex].Y()
					>> pts[vertex].Z();


				vertex++;
				if (vertex == 3)
				{
					if (normal.Norm() <= 1e-5)
					{
						Vector3D v1 = pts[1]-pts[0];
						Vector3D v2 = pts[2]-pts[0];
						normal = v1.Cross(v2);
						if (normal.Norm() != 0) normal *= 1./normal.Norm();		  
					}

					vertex = 0;

					if ( (Dist (pts[0], pts[1]) < 1e-16) ||
						(Dist (pts[0], pts[2]) < 1e-16) ||
						(Dist (pts[1], pts[2]) < 1e-16) )
					{
						;
					}
					else
					{
						newnormals.Add(normal);
						newpoints.Add(pts[0]);
						newpoints.Add(pts[1]);
						newpoints.Add(pts[2]);
					}
				}
			}
		}

	}	//end read STL ASCII format
	else
	{ 
		//read STL BINARY format

		int binary_mode = 1;
		CMFile infile(file, TFMread, binary_mode);

		Vector3D pts[3];
		Vector3D normal;

		int nofacets = 0;
		int vertex = 0;

		if (sizeof(int) != 4 || sizeof(float) != 4)
			UO() << "for stl-binary compatibility only use 32 bit compilation!!!";

		//specific settings for stl-binary format
		const int namelen = 80; //length of name of header in file
		const int nospaces = 2; //number of spaces after a triangle

		//read header: name
		char buf[namelen+1];

		mystr str;
		//read leading 80 characters
		infile.RWchars(namelen, str);
		//global_uo << "STL-header=" << str << "\n";

		//read 4-byte int number of facets
		infile.RWbinaryInt(nofacets);
		//global_uo << "No. facets=" << nofacets << "\n";

		float f;
		char spaces[nospaces];

		//read triangles:
		int cntface, j;
		for (cntface = 1; cntface <= nofacets; cntface++)
		{
			//read 3 floats for triangle normal
			infile.RWbinaryFloat(f); normal.X() = (double)f;
			infile.RWbinaryFloat(f); normal.Y() = (double)f;
			infile.RWbinaryFloat(f); normal.Z() = (double)f;

			//read 3 x 3 floats for (x,y,z) coordinates of triangle points
			for (j = 0; j < 3; j++)
			{
				infile.RWbinaryFloat(f); pts[j].X() = (double)f;
				infile.RWbinaryFloat(f); pts[j].Y() = (double)f;
				infile.RWbinaryFloat(f); pts[j].Z() = (double)f;	  
			} 

			//read two bytes for attributes (usually unused)
			infile.RWchar(spaces[0]); 
			infile.RWchar(spaces[1]);
			//global_uo << "sp1=" << (int)spaces[0] << ", sp1=" << (int)spaces[1] << "\n";

			newnormals.Add(normal);
			newpoints.Add(pts[0]);
			newpoints.Add(pts[1]);
			newpoints.Add(pts[2]);
		}	    

	}

	//now sort points and generate a consistent surface mesh (not necessary for graphical representation)
	Box3D box;
	for (int i=1; i <= newpoints.Length(); i++) box.Add(newpoints(i));

	IVector pnums;
	pnums.SetLen(newpoints.Length());
	IVector items;

	double boxeps = stl_tolerance+1e-12;
	double eps2 = Sqr(stl_tolerance);

	TArray<Vector3D> locpoints;

	SearchTree st(20,20,20, box);
	for (int i=1; i <= newpoints.Length(); i++)
	{
		Box3D b(newpoints(i), newpoints(i));
		b.Increase(boxeps);
		st.GetItemsInBox(b, items);

		int found = 0;
		for (int j=1; j <= items.Length(); j++)
		{
			if (Dist2(newpoints(i), locpoints(items(j))) <= eps2) 
			{
				found = items(j);
				break;
			}
		}
		if (found)
		{
			pnums(i) = found;
		}
		else
		{
			pnums(i) = locpoints.Add(newpoints(i));
			st.AddItem(Box3D(newpoints(i),newpoints(i)), pnums(i));
		}
	}

	Matrix m(newnormals.Length(), 3);
	Matrix n(newnormals.Length(), 3);
	for (int i=1; i <= newnormals.Length(); i++)
	{
		m(i,1) = pnums(3*i-2);
		m(i,2) = pnums(3*i-1);
		m(i,3) = pnums(3*i);
		n(i,1) = newnormals(i).X();
		n(i,2) = newnormals(i).Y();
		n(i,3) = newnormals(i).Z();
	}

	//ElementDataContainer edc;

	ElementData ed;
	ed.SetMatrix(m.GetMatPtr(),m.Getrows(),m.Getcols(),"triangles");
	return_value->Add(ed);

	ed.SetMatrix(n.GetMatPtr(),n.Getrows(),n.Getcols(),"normals");
	return_value->Add(ed);

	Matrix p(newnormals.Length(), 3);
	for (int i=1; i <= locpoints.Length(); i++)
	{
		p(i,1) = locpoints(i).X();
		p(i,2) = locpoints(i).Y();
		p(i,3) = locpoints(i).Z();
	}

	ed.SetMatrix(p.GetMatPtr(),p.Getrows(),p.Getcols(),"points");
	return_value->Add(ed);

	//return_value.SetEDC(edc.GetCopy(),"MeshData");
	return 1;
}

//this function computes mass, moment of inertia, volume and center of gravity based on the data about the geometry and the material
int MultiBodySystem::ComputeInertia(ElementDataContainer* data, ElementDataContainer* return_value)
{
	int material_number = 0;
	double rho = 0;

	// set rho according to material_number or directly entered density, if there is no material_number
	if(!data->Find("material_number"))	// if there is no material_number search for entry "density"
	{
		rho = data->TreeGetDouble("density",0);
	}
	else
	{
		material_number = data->TreeGetInt("material_number",0);		
		if(material_number <= NMaterials())
		{
			rho = GetMaterial(material_number).Density();
		}
		else
		{
			return 0;		// not a valid material_number
		}
	}

	if(!rho) 	{	return 0; }// density == 0
	// rho is now set

	// start the computation of the inertia values
	Matrix3D Iphi(0.);
	double mass, volume;
	Vector3D cog;

	if(data->Find("MeshData"))
	{
		//see Maple file "maple/integration_triangle_volume.mws" and paper 'computing moments of piecewise polynomial surfaces'
		double im11=0, im12=0, im13=0, im22=0, im23=0, im33=0;
		double vol=0;
		Vector3D center(0.,0.,0.);

		int dummy = 0;
		int Npoints, Ntrigs;

		double * p;
		double * t;
		data->TreeGetMatrix("MeshData.points",&p,Npoints,dummy);
		if(dummy!=3) {return 0;}	// error
		data->TreeGetMatrix("MeshData.triangles",&t,Ntrigs,dummy);
		if(dummy!=3) {return 0;}	// error

		Matrix pts(Npoints,3,p);
		Matrix trigs(Ntrigs,3,t);

		for (int i=1; i <= Ntrigs; i++)
		{
			//const Vector3D& p1 = locpoints(trigs(i).Get(1));
			//const Vector3D& p2 = locpoints(trigs(i).Get(2));
			//const Vector3D& p3 = locpoints(trigs(i).Get(3));

			const Vector3D& p1 = Vector3D(pts(trigs(i,1),1),pts(trigs(i,1),2),pts(trigs(i,1),3));
			const Vector3D& p2 = Vector3D(pts(trigs(i,2),1),pts(trigs(i,2),2),pts(trigs(i,2),3));
			const Vector3D& p3 = Vector3D(pts(trigs(i,3),1),pts(trigs(i,3),2),pts(trigs(i,3),3));

			double v1x = p2.X()-p1.X();
			double v1y = p2.Y()-p1.Y();
			double v1z = p2.Z()-p1.Z();

			double v2x = p3.X()-p1.X();
			double v2y = p3.Y()-p1.Y();
			double v2z = p3.Z()-p1.Z();
			double px = p1.X();
			double py = p1.Y();
			double pz = p1.Z();

			im11 += v1x*v1x*v1z*(v1x*v2y-v2x*v1y)/20.0+px*px*pz*(v1x*v2y-v2x*v1y)/2.0+(
				2.0*px*v1x*v1z+v1x*v1x*pz)*(v1x*v2y-v2x*v1y)/12.0+(2.0*v2x*v1x*v1z+v1x*v1x*v2z)
				*(v1x*v2y-v2x*v1y)/60.0+(2.0*px*v2x*v1z+2.0*px*v1x*v2z+2.0*v2x*v1x*pz)*(v1x*v2y
				-v2x*v1y)/24.0+(v2x*v2x*v1z+2.0*v2x*v1x*v2z)*(v1x*v2y-v2x*v1y)/60.0+(px*px*v1z+
				2.0*px*v1x*pz)*(v1x*v2y-v2x*v1y)/6.0+(px*px*v2z+2.0*px*v2x*pz)*(v1x*v2y-v2x*v1y
				)/6.0+(2.0*px*v2x*v2z+v2x*v2x*pz)*(v1x*v2y-v2x*v1y)/12.0+v2x*v2x*v2z*(v1x*v2y-
				v2x*v1y)/20.0;
			double MapleGenVar1 = ((px*v1y+v1x*py)*v1z+v1x*v1y*pz)*(v1x*v2y-v2x*v1y)/12.0+((
				v2x*v1y+v1x*v2y)*v1z+v1x*v1y*v2z)*(v1x*v2y-v2x*v1y)/60.0+((px*v2y+v2x*py)*v1z+(
				px*v1y+v1x*py)*v2z+(v2x*v1y+v1x*v2y)*pz)*(v1x*v2y-v2x*v1y)/24.0+(v2x*v2y*v1z+(
				v2x*v1y+v1x*v2y)*v2z)*(v1x*v2y-v2x*v1y)/60.0+(px*py*v2z+(px*v2y+v2x*py)*pz)*(
				v1x*v2y-v2x*v1y)/6.0;
			im12 += MapleGenVar1+((px*v2y+v2x*py)*v2z+v2x*v2y*pz)*(v1x*v2y-v2x*v1y)/12.0
				+px*py*pz*(v1x*v2y-v2x*v1y)/2.0+(px*py*v1z+(px*v1y+v1x*py)*pz)*(v1x*v2y-v2x*v1y
				)/6.0+v1x*v1y*v1z*(v1x*v2y-v2x*v1y)/20.0+v2x*v2y*v2z*(v1x*v2y-v2x*v1y)/20.0;
			im13 += px*pz*pz*(v1x*v2y-v2x*v1y)/4.0+v1x*v1z*v1z*(v1x*v2y-v2x*v1y)/40.0+
				v2x*v2z*v2z*(v1x*v2y-v2x*v1y)/40.0+(px*v1z*v1z+2.0*v1x*pz*v1z)*(v1x*v2y-v2x*v1y
				)/24.0+(v2x*v1z*v1z+2.0*v1x*v2z*v1z)*(v1x*v2y-v2x*v1y)/120.0+(2.0*(px*v2z+v2x*
				pz)*v1z+2.0*v1x*pz*v2z)*(v1x*v2y-v2x*v1y)/48.0+(2.0*v2x*v2z*v1z+v1x*v2z*v2z)*(
				v1x*v2y-v2x*v1y)/120.0+(2.0*px*pz*v2z+v2x*pz*pz)*(v1x*v2y-v2x*v1y)/12.0+(px*v2z
				*v2z+2.0*v2x*pz*v2z)*(v1x*v2y-v2x*v1y)/24.0+(2.0*px*pz*v1z+v1x*pz*pz)*(v1x*v2y-
				v2x*v1y)/12.0;
			im22 += py*py*pz*(v1x*v2y-v2x*v1y)/2.0+v2y*v2y*v2z*(v1x*v2y-v2x*v1y)/20.0+
				v1y*v1y*v1z*(v1x*v2y-v2x*v1y)/20.0+(2.0*py*v1y*v1z+v1y*v1y*pz)*(v1x*v2y-v2x*v1y
				)/12.0+(2.0*v2y*v1y*v1z+v1y*v1y*v2z)*(v1x*v2y-v2x*v1y)/60.0+(2.0*py*v2y*v1z+2.0
				*py*v1y*v2z+2.0*v2y*v1y*pz)*(v1x*v2y-v2x*v1y)/24.0+(v2y*v2y*v1z+2.0*v2y*v1y*v2z
				)*(v1x*v2y-v2x*v1y)/60.0+(py*py*v1z+2.0*py*v1y*pz)*(v1x*v2y-v2x*v1y)/6.0+(py*py
				*v2z+2.0*py*v2y*pz)*(v1x*v2y-v2x*v1y)/6.0+(2.0*py*v2y*v2z+v2y*v2y*pz)*(v1x*v2y-
				v2x*v1y)/12.0;
			im23 += v2y*v2z*v2z*(v1x*v2y-v2x*v1y)/40.0+py*pz*pz*(v1x*v2y-v2x*v1y)/4.0+
				v1y*v1z*v1z*(v1x*v2y-v2x*v1y)/40.0+(py*v1z*v1z+2.0*v1y*pz*v1z)*(v1x*v2y-v2x*v1y
				)/24.0+(v2y*v1z*v1z+2.0*v1y*v2z*v1z)*(v1x*v2y-v2x*v1y)/120.0+(2.0*(py*v2z+v2y*
				pz)*v1z+2.0*v1y*pz*v2z)*(v1x*v2y-v2x*v1y)/48.0+(2.0*v2y*v2z*v1z+v1y*v2z*v2z)*(
				v1x*v2y-v2x*v1y)/120.0+(2.0*py*pz*v2z+v2y*pz*pz)*(v1x*v2y-v2x*v1y)/12.0+(py*v2z
				*v2z+2.0*v2y*pz*v2z)*(v1x*v2y-v2x*v1y)/24.0+(2.0*py*pz*v1z+v1y*pz*pz)*(v1x*v2y-
				v2x*v1y)/12.0;
			im33 += v2z*v1z*v1z*(v1x*v2y-v2x*v1y)/60.0+v1z*v1z*v1z*(v1x*v2y-v2x*v1y)/
				60.0+pz*v2z*v1z*(v1x*v2y-v2x*v1y)/12.0+pz*pz*v1z*(v1x*v2y-v2x*v1y)/6.0+v2z*v2z*
				v1z*(v1x*v2y-v2x*v1y)/60.0+pz*pz*v2z*(v1x*v2y-v2x*v1y)/6.0+pz*v2z*v2z*(v1x*v2y-
				v2x*v1y)/12.0+pz*pz*pz*(v1x*v2y-v2x*v1y)/6.0+v2z*v2z*v2z*(v1x*v2y-v2x*v1y)/60.0
				+pz*v1z*v1z*(v1x*v2y-v2x*v1y)/12.0;

			vol += 1./6.*v1z*(v1x*v2y-v2x*v1y)+1./6.*v2z*(v1x*v2y-v2x*v1y)+1./2.*pz*(v1x*v2y-v2x*v1y);

			center.X() += v1x*v1z*(v1x*v2y-v2x*v1y)/12.0+(v2x*v1z+v1x*v2z)*(v1x*v2y-v2x*v1y)/24.0+v2x*v2z*(v1x*v2y-v2x*v1y)/12.0+(px*v1z+v1x*pz)*(v1x*v2y-v2x*v1y)/6.0+(px*v2z+v2x*pz)*(v1x*v2y-v2x*v1y)/6.0+px*pz*(v1x*v2y-v2x*v1y)/2.0;
			center.Y() += v1y*v1z*(v1x*v2y-v2x*v1y)/12.0+(v2y*v1z+v1y*v2z)*(v1x*v2y-v2x*v1y)/24.0+v2y*v2z*(v1x*v2y-v2x*v1y)/12.0+(py*v1z+v1y*pz)*(v1x*v2y-v2x*v1y)/6.0+(py*v2z+v2y*pz)*(v1x*v2y-v2x*v1y)/6.0+py*pz*(v1x*v2y-v2x*v1y)/2.0;
			center.Z() += v1z*v1z*(v1x*v2y-v2x*v1y)/24.0+v2z*v1z*(v1x*v2y-v2x*v1y)/24.0+v2z*v2z*(v1x*v2y-v2x*v1y)/24.0+pz*v1z*(v1x*v2y-v2x*v1y)/6.0+pz*v2z*(v1x*v2y-v2x*v1y)/6.0+pz*pz*(v1x*v2y-v2x*v1y)/4.0;
		}
		volume = vol;
		mass = volume*rho;
		Iphi = rho * Matrix3D(im22+im33,-im12,-im13, -im12, im11+im33, -im23, -im13, -im23, im11+im22);
		if (vol != 0) 
		{
			cog = center * (1./vol);
		}
	}
	else if(data->Find("Cube"))
	{
		double x,y,z;
		data->TreeGetVector3D("Cube.body_dimensions",x,y,z);
		
		volume = x*y*z;
		mass = volume*rho;
		Iphi(1,1)= 1./12.*mass*(Sqr(y)+Sqr(z));
		Iphi(2,2)= 1./12.*mass*(Sqr(x)+Sqr(z));
		Iphi(3,3)= 1./12.*mass*(Sqr(x)+Sqr(y));
		cog = Vector3D(0.);
	}
	else if(data->Find("Cylinder"))
	{

	}

	// set the new edc
	ElementData ed;
	ed.SetDouble(volume,"volume"); return_value->Add(ed);
	ed.SetDouble(mass,"mass"); return_value->Add(ed);
	ed.SetVector3D(cog.X(),cog.Y(),cog.Z(),"center_of_mass"); return_value->Add(ed);
	ed.SetMatrix(Iphi.GetMatPtr(),3,3,"moment_of_inertia"); return_value->Add(ed);

	return 1;

}

//this function reads a column of a file and returns it as vector, rv = 1.. success, 0..file is not good
int MultiBodySystem::LoadVectorFromFile(const char* filename, int col, ElementDataContainer* return_value)
{
	CMatrixFile file( filename, TFMread );
	if(file.IsGood())
	{
		file.ReadSingleColumn(col);
		TArray<double>& coldata = file.Column(col);
		Vector test;
		test.SetVector(coldata);
		// set the new edc
		ElementData ed;
		ed.SetVector(test.GetVecPtr(), test.Length(),"vector"); return_value->Add(ed);
		return 1;
	}
	else
	{
		return 0;
	}
}

void MultiBodySystem::AddReplaceModelDataEDC(/*const */ElementDataContainer& edc)
{
	//$ DR 2013-06-12:[ code that has been in AddModelData
	// verbosity level is needed for further opertation
	if (edc.TreeFind("LoggingOptions.output_level"))    // temporary solution - best would be to have an entry like "GeneralOptions.show_general_warnings" in edc_file     //$ PG 2012-2-3: substituted SolverOptions.Log.output_level by LoggingOptions.output_level
	{
		GetOptions()->LoggingOptions()->OutputLevel() = edc.TreeGetInt("LoggingOptions.output_level");       //$ PG 2012-2-3: substituted solset.output_level by GetOptions()->LoggingOptions()->OutputLevel() and SolverOptions.Log.output_level by LoggingOptions.output_level
	}
	//$ DR 2013-06-12:] code that has been in AddModelData

	bool show_warnings = ( GetOptions()->LoggingOptions()->OutputLevel() >= UO_LVL_dbg1 );    //$ PG 2012-2-3: substituted solset.output_level by GetOptions()->LoggingOptions()->OutputLevel()

	ElementData* edmodel;

	// loop over all default main branches
	for (int i=1; i <= edc.Length(); i++) 
	{
		// get i-th main branch
		if(!GetMBS_EDC_Options()->TreeFind(edc.GetPtr(i)->GetDataName()))
		{
			ElementDataContainer edci;edci.TreeAdd("",edc.Get(i));
			edc_modeldata->TreeReplaceEDCDataWith(&edci, show_warnings, mystr(edc.GetPtr(i)->GetDataName())+mystr("."));
		}
	}

	//$ DR 2013-06-12:[ code that has been in AddModelData
	// replace solveroptions & other settings
	ElementData* eddef;
	ElementData* edfile;

	// loop over all default main branches
	for (int i=1; i <= GetMBS_EDC_Options()->Length(); i++) 
	{
		// get i-th main branch
		//eddef = edc_default.GetPtr(i);		
		eddef = GetMBS_EDC_Options()->GetPtr(i);
		edfile = edc.TreeFind(eddef->GetDataName());
		// replace the default settings by the user's file settings
		if (edfile != 0 && eddef->IsEDC() && edfile->IsEDC())
		{
			if(GetMBS_EDC_Options()->TreeFind(edfile->GetDataName()))
			{				
				UO(UO_LVL_ext) << "replacing hotint option: " << edfile->GetDataName() << "\n";
				// use file parameters (when existing) instead of default variables
				eddef->GetEDC()->TreeReplaceEDCDataWith(edfile->GetEDC(), show_warnings, mystr(edfile->GetDataName())+mystr("."));
			}
		}
	}
	// overwrite default options by user file data 
	EDC2SolverOptions(GetMBS_EDC_Options());
	EDC2Options(GetMBS_EDC_Options());       // copy EDC to i/d/toptions

	OpenLogFile(0);

	UO() << "done.\n";
	UO() << "--------------------------------\n";
	//$ DR 2013-06-12:] code that has been in AddModelData
}


int MultiBodySystem::ReadModelData(mystr filename)
{
	ElementDataContainer edc_file;

	// read parameters from file and store it in ElementDataContainer of mbs (details: see forum)
	CMFile file(filename, TFMread);
	UO()	<< "--------------------------------\n";
	int rv = 1;

	//read parameters from file
	if (file.IsGood()) 
	{
		UO(UO_LVL_0) << "read " << filename << "\n";

		//read file
		mystr str;
		file.RWF(str);		

		//$ DR 2013-03-09:[
		// Get the variables of the parameter variation
		ElementDataContainer edc_tmp;
		ElementData ed_tmp;
		if(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.activate"))
		{
			mystr varname = GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.MBS_EDC_variable_name");
			double val = edc_modeldata->TreeGetDouble(varname);
			ed_tmp.SetDouble(val,varname);
			edc_tmp.Add(ed_tmp);
			if(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.Var2.activate"))
			{
				varname = GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.Var2.MBS_EDC_variable_name");
				val = edc_modeldata->TreeGetDouble(varname);
				ed_tmp.SetDouble(val,varname);
				edc_tmp.Add(ed_tmp);
			}
		}		
		//$ DR 2013-03-09:]
		//$ RE & MSax 2013-06-27:[
		else if(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.activate") && edc_modeldata)
		{
			int number_of_params = GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Parameters.number_of_params");
			for (int col=1; col <= number_of_params; col++)
			{
				mystr varname =	GetMBS_EDC_Options()->TreeGetString(mystr("SolverOptions.Optimization.Parameters.param_name")+mystr(col));
				double val = edc_modeldata->TreeGetDouble(varname);
				ed_tmp.SetDouble(val,varname);
				edc_tmp.Add(ed_tmp);
			}
		}
		//$ RE & MSax 2013-06-27:]
		


//$!AD 26-09-2013: parse line by line directly into modeldata EDC [
// NEW CODE:
		edc_modeldata = 0;
		if (edc_modeldata != 0)
		{ 
			delete edc_modeldata;
			edc_modeldata = 0;
		}
		edc_modeldata = new ElementDataContainer;

		// overwrite edc with contents from file
		//edc_modeldata = &edc_file; //JG, this is just temporary in order to get variables, change to MBSParser->readedc_global/local!!!
		ElementData ed;
		ed.SetText(filename, "last_file_name");edc_modeldata->Add(ed); // path is stored e.g. for getting relative paths in include command

		//$ DR 2013-03-09:[
		// set the values for the parameter variation and save them with flagging them as readonly
		for(int i=1; i<=edc_tmp.Length();i++)
		{
			edc_modeldata->Add(edc_tmp.Get(i));
			edc_modeldata->Last().SetLocked(1);
		}
		edc_tmp.Reset();
		//$ DR 2013-03-09:]
 
		this->EDCParser().GetMBSParser()->SetLocalEDC2(edc_modeldata);
		rv = this->EDCParser().ParseAndExecuteString(str, *edc_modeldata);
		edc_file = *(edc_modeldata->GetCopy());	//$ DR+AD 2013-10-01 bugfix, copy is necessary for solver options

// OLD CODE:
		////////
		//////////UO() << "WARNING:set read EDC==>modeldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		////////rv = EDCParser().String2TreeEDC(/*this, */str, edc_file);  
		//////////UO() << "delete read EDC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		//////////delete read EDC
		////////edc_modeldata = 0;

		////////SetModelDataContainer(edc_file); //store empty edc to avoid crashes (due to null-pointer exeptions) at start of Generate - function
		////////
//$!AD ]
	}
	else // FileNotGood
	{
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "ERROR: failed to read model data from " << filename << "\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		UO()	<< "*****************\n";
		SetModelDataContainer(edc_file); //store empty edc to avoid crashes (due to null-pointer exeptions) at start of Generate - function

		return 0;  //$JG2012-02-21 changed from 2 to 0 in case of error
	}

	// verbosity level is needed for further opertation
	if (edc_file.TreeFind("LoggingOptions.output_level"))    // temporary solution - best would be to have an entry like "GeneralOptions.show_general_warnings" in edc_file         //$ PG 2012-2-3: substituted SolverOptions.Log.output_level by LoggingOptions.output_level
	{
		GetOptions()->LoggingOptions()->OutputLevel() = edc_file.TreeGetInt("LoggingOptions.output_level");    //$ PG 2012-2-3: substituted solset.output_level by GetOptions()->LoggingOptions()->OutputLevel(), and SolverOptions.Log.output_level by LoggingOptions.output_level
	}
	bool show_warnings = ( GetOptions()->LoggingOptions()->OutputLevel() >= UO_LVL_warn );     //$ PG 2012-2-3: substituted solset.output_level by GetOptions()->LoggingOptions()->OutputLevel()

	// load default parameter settings
	ElementDataContainer* edc_default = GetMBS_EDC_Options();

	RemoveVariableFromEDCTree(*edc_default, mystr("SolverOptions.ComputationSteps"), 1);

	ElementData* eddef;
	ElementData* edfile;

	//$ AD: 2011-11-08: special treatment for Computation Step entries in EDC 
	// as these entries would produce a popupwindow every time when using TreeReplace
	// copy computation steps manually, assume that entry does NOT exist in default EDC
	ElementData* edcsfile = edc_file.TreeFind("SolverOptions.ComputationSteps");
	if(edcsfile != 0 && edcsfile->IsEDC())
	{
		ElementDataContainer* edc_solopt_default = edc_default->TreeFind("SolverOptions")->GetEDC();
		edc_solopt_default->Add(*edcsfile);
		UO(UO_LVL_ext) << "Added Computation Steps: (" << edcsfile->GetEDC()->Length() << ")\n";
	// one could eliminate entries from FileEDC now...
		RemoveVariableFromEDCTree(edc_file, mystr("SolverOptions.ComputationSteps"), 0);
	}


	// loop over all default main branches
	for (int i=1; i <= edc_default->Length(); i++) 
	{

		// get i-th main branch
		eddef = edc_default->GetPtr(i);		
		edfile = edc_file.TreeFind(eddef->GetDataName());

		// replace the default settings by the user's file settings
		if (edfile != 0 && eddef->IsEDC() && edfile->IsEDC())
		{
			UO(UO_LVL_ext) << "replacing defaults: " << edfile->GetDataName() << "\n";

			// use file parameters (when existing) instead of default variables
			eddef->GetEDC()->TreeReplaceEDCDataWith(edfile->GetEDC(), show_warnings, mystr(edfile->GetDataName())+mystr("."));
		}
	}

	// overwrite default options by user file data 
	EDC2SolverOptions(edc_default);
	EDC2Options(edc_default);       // copy EDC to i/d/toptions

	UO() << "done.\n";
	UO() << "--------------------------------\n";

	return rv; //1==success, 0==error
}

// use add model data after command read model data; assumption: model data container is not empty!!!
//$ DR 2013-07-25: assumption is not necessary anymore.
//$ AD 2013-10-22: uses ParseAndExecuteString to instantly write to edc_modeldata
int MultiBodySystem::AddModelData(mystr filename)
{
	UO(UO_LVL_0) << "read " << filename << "\n";

	ElementDataContainer edc_file;
	
	// store file name for include command with relative paths
	ElementData ed;
	ed.SetText(filename, "last_file_name");
	if(edc_modeldata != 0)
	{
		edc_modeldata->TreeSet("",ed);
	}
	else
	{
		SetModelDataContainer(edc_file); //store empty edc to avoid crashes (due to null-pointer exeptions) 
	}

//$ AD 2013-10-22: read the string and parse it directly into the modeldatacontainer
	mystr str;
	int rv = File2Str(filename, str);
	if(!rv)return rv; 
	
	EDCParser().ParseAndExecuteString(str, *edc_modeldata);

	// replace solveroptions & other settings
	bool show_warnings = GetOptions()->LoggingOptions()->OutputLevel() >= UO_LVL_dbg1;
	ElementData* eddef;
	ElementData* edfile;

	// loop over all default main branches
	for (int i=1; i <= GetMBS_EDC_Options()->Length(); i++) 
	{
		// get i-th main branch
		//eddef = edc_default.GetPtr(i);		
		eddef = GetMBS_EDC_Options()->GetPtr(i);
		edfile = edc_modeldata->TreeFind(eddef->GetDataName());
		// replace the default settings by the user's file settings
		if (edfile != 0 && eddef->IsEDC() && edfile->IsEDC())
		{
			if(GetMBS_EDC_Options()->TreeFind(edfile->GetDataName()))
			{				
				UO(UO_LVL_ext) << "replacing hotint option: " << edfile->GetDataName() << "\n";
				// use file parameters (when existing) instead of default variables
				eddef->GetEDC()->TreeReplaceEDCDataWith(edfile->GetEDC(), show_warnings, mystr(edfile->GetDataName())+mystr("."));
			}
		}
	}
	// overwrite default options by user file data 
	EDC2SolverOptions(GetMBS_EDC_Options());
	EDC2Options(GetMBS_EDC_Options());       // copy EDC to i/d/toptions

	OpenLogFile(0);

	UO() << "done.\n";
	UO() << "--------------------------------\n";

	return 1; //success
}

// remove a variable or branch from the edc by name - this is used to prevent ComputationSteps to be written to file and to be remembered from a previous model		
int MultiBodySystem::RemoveVariableFromEDCTree(ElementDataContainer& edc, mystr& varname, int warning)
{
	mystr treename, elementname;
	edc.SplitIntoTreeAndElementName(varname, treename, elementname);
	ElementData* ed_branch = edc.TreeFind(varname);
	
	if (ed_branch != 0) // variable exists
	{
		ElementDataContainer* edc_nesting;
		if(treename.Length())
		{
			edc_nesting = edc.TreeFind(treename)->GetEDC();
		}
		else
		{
			edc_nesting = &edc;
		}

		if (warning)
		{
			if( ed_branch -> IsEDC() ) UO(UO_LVL_ext) << "removing sub-edc " + varname + " from edc\n";
			else UO(UO_LVL_ext) << "removing entry " + varname + "from edc \n";
		}

		edc_nesting->Delete(edc_nesting->Find(elementname));
		return 1;
	}
	return 0;
}

int MultiBodySystem::AddFileSensor(const char* filename, const mystr sensor_name, const int ncolumns, const int colTime, const int colSignal, const int interpolation, const int nrOfHeaderLines)
{
	/*
	MathFunction mf;
	if(mf.SetPiecewiseFromFile(filename, ncolumns, colTime, colSignal, interpolation, nrOfHeaderLines)) // interpolation: 0...steps, 1...piecewise linear
	{
		InstantMessageText(mystr("Error1: File " + mystr(filename) + " not found!"));
		return 1; // error
	}
	Vector2D posDraw(0.,0.); // drawing position
	Vector3D drawDim = Vector3D(0.2,0.2,0.); // draw dimension of control element
	IOTime lt_time(this);
	lt_time.SetElementName(sensor_name+"-IOTime");
	lt_time.SetRefPos2D(Vector2D(0.,3.0)+posDraw);//define reference position for drawing
	lt_time.SetDrawDim(drawDim);     //define size of rectangle for drawing (B, H, dummy)
	int nLt_time = AddElement(&lt_time);

	IOMathFunction iom_ty(this);
	iom_ty.SetIOMathFunction(mf);
	iom_ty.SetRefPos2D(Vector2D(0.5,3.0)+posDraw);
	iom_ty.SetDrawDim(drawDim);
	iom_ty.SetElementName(sensor_name+"-IOMathFcn");						//Change name of element from default name				
	iom_ty.AddInput(nLt_time, IOInputTypeElement, 1);  // connection(nLt_time,  typ=io, output 1)
	int nIo_ty = AddElement(&iom_ty);

	MBSSensor sens_mes(this, TMBSSensor(TSOutputSensor), nIo_ty, 1); //measure output 1
	sens_mes.SetSensorName(sensor_name);
	AddSensor(&sens_mes);
	*/
	return 0; //success
}
int MultiBodySystem::AddFileSensor2(const char* filename, const mystr sensor_name, const int colTime, const int colSignal, const int interpolation, mystr comment)
{
	/*
	MBSSensor sens_mes(this, TMBSSensor(TSFile), colTime, colSignal); //measure output 1
	if(sens_mes.SetTSFile(filename, colTime, colSignal, interpolation, comment)) // set file name
	{
		return 1; // error, don't add sensor
	}
	sens_mes.SetSensorName(sensor_name);
	AddSensor(&sens_mes);	
	*/
	return 0; //success
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void MultiBodySystem::SetOptions2EDCOptions()
{
	ElementDataContainer* edc = GetMBS_EDC_Options();
	ElementDataContainer edcnew;
	Options2EDC(&edcnew);                //copy options to edc:
	edc->TreeReplaceEDCDataWith(&edcnew);//replace edcnew data in edc:
}

void MultiBodySystem::SetEDCOptions2Options()
{
	ElementDataContainer* edc = GetMBS_EDC_Options();
	EDC2Options(edc);
}

void MultiBodySystem::SetComputationSolverOptions()
{
	ElementDataContainer* edc = GetMBS_EDC_Options();
	EDC2SolverOptions(edc);
}

void MultiBodySystem::SetSolverDialogOptions()
{
	//update computation options with solver settings:

	ElementDataContainer* edc = GetMBS_EDC_Options();
	ElementDataContainer edcnew;

	//copy mbssolset to edc:
	SolverOptions2EDC(&edcnew);
	//replace edcnew data in edc:
	edc->TreeReplaceEDCDataWith(&edcnew);

}

