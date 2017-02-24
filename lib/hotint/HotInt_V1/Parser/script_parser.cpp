//#**************************************************************
//#
//# filename:             script_parser.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
//# description:          parser for creating objects in script language
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
#include "mbs_interface.h"
#include "element.h"
#include "options_class_auto.h"
#include "script_parser.h"
//$ YV 2012-12-14: the following include is to be removed; general object factory should be used in the future
//#include "sensorsSpecific.h"	//$ DR 2012-10: necessary to include for CEDCParser
#include "elementdataaccess.h"
#include "graphicsconstants.h"
#include "material.h"
#include "geomelements.h"
#include "femathhelperfunctions.h"
#include "sensors.h"
#include "parser.h"
#include "elementsandmodelslibraryinterface.h"
#include "node.h"



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++EvaluateTextualParameterEntry++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDC_Command::EvaluateTextualParameterEntry(MBS *mbs, CEDCParser *edcParser, ElementDataContainer *param_edc)
{
// convert any textual entry to scalar,vector,matrix
// assume that the parameter EDC has only a single entry
	ElementData* edp = param_edc->GetPtr(1);
	if(edp->IsText())
	{
		mystr entry = mystr(edp->GetDataName()) + mystr("= ") + mystr(edp->GetText());
		ElementDataContainer tmp;

		// evaluate the textual entry
		edcParser->ParseAndExecuteString(entry,tmp);

		param_edc->Delete(1);
		param_edc->TreeReplaceEDCDataWith(&tmp);

		return 1;
	}
	return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++ExecuteCommand ++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDC_ComAddElement::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option) 
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int element_type_index = edc->Find("element_type"); //find index in element EDC
	if (element_type_index)
	{
		if (!edc->Get(element_type_index).IsText()) {assert(0);}

		int listindex;
		mystr element_type_str = edc->Get(element_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCElement, element_type_str);
		if (type_id < 0 ) 
		{
			edcParser->EDCError(mystr("Element type '") +  element_type_str + mystr("' not found in command AddElement(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		} //does not exist in element type list

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCElement, type_id); //this are the typeflags, e.g. constraint, 2D, etc.

		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + element_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			if (IsAddElement() && (typeflags & TAEconstraint)) //constraints are in the same list, but should not be added with "AddElement"
			{
				edcParser->EDCError(mystr("Connector '")+element_type_str+mystr("' can not be added with command AddElement(...), use AddConnector(...) instead!"));
			}
			else if (!IsAddElement() && !(typeflags & TAEconstraint)) //constraints are in the same list, but should not be added with "AddElement"
			{
				edcParser->EDCError(mystr("Element '")+element_type_str+mystr("' can not be added with command AddConnector(...), use AddElement(...) instead!"));
			}

			//add new element of selected type to MBS
			//int elnum = edcParser->AddElement(type_id, option); //old: AddElement(element_type_index, option);
			int elnum = edcParser->GetObjectFactory()->AddObject(OFCElement, type_id);

			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetElement(elnum).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetElement(elnum).SetElementData(edc);

			return_value.SetInt(elnum, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'element_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDC_ComAddGeomElement::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int element_type_index = edc->Find("geom_element_type"); //find index in element EDC
	if (element_type_index)
	{
		if (!edc->Get(element_type_index).IsText()) {assert(0);}

		int listindex;
		mystr element_type_str = edc->Get(element_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCGeomElement, element_type_str);
		if (type_id < 0) 
		{
			edcParser->EDCError(mystr("Element type '") +  element_type_str + mystr("' not found in command AddGeomElement(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCGeomElement, type_id); //this are the typeflags, e.g. constraint, 2D, etc.

		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + element_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new element of selected type to MBS
			//int geom_elnum=edcParser->AddGeomElement(type_id,0); //Add GeomElement with default values
			int geom_elnum = edcParser->GetObjectFactory()->AddObject(OFCGeomElement, type_id); //Add GeomElement with default values
			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetDrawElement(geom_elnum)->GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetDrawElement(geom_elnum)->SetElementData(edc);

			return_value.SetInt(geom_elnum, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'geom_element_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//$ DR 2012-11

mystr CEDC_ComPrint::GetCommandTexParameterDescription() const 
{
	mystr descr;
	descr += "There are three possibilities to use the command. The parameter can either be: \\\\ \n";
	descr += mystr("\\begin{itemize} \n");
	descr+= mystr("  \\item a text, e.g. Print(\"Hello world\")  \n");
	descr+= mystr("  \\item an ElementDataContainer, e.g. Print(my\\_mass)  \n");
	descr+= mystr("  \\item an ElementData, e.g. Print(my\\_mass.density)  \n");
	descr += mystr("\\end{itemize} \n");
	descr += mystr("In the case of a text or an ElementData, only the text itself is printed. In the case of an ElementDataContainer, also the name of the ElementData is printed.");
	return descr;
}


int CEDC_ComPrint::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	mystr str, num;

	//symbols defined in CEDCParser:
	mystr el = "\n"; //end of line
	mystr ob = "{";  //open brace
	mystr cb = "}";  //closing brace
	mystr osb = "[";  //open square bracket
	mystr csb = "]";  //closing square bracket
	mystr ds = ": "; //double stop
	mystr semicolon = ";"; //semicolon
	char commentc = '%'; //comment
	mystr doublequote = mystr('"'); //for strings or text

	//symbols with additional whitespaces:
	mystr eq = "= "; //equal sign
	mystr comma = ", "; //comma
	mystr comment = mystr("  ")+mystr(commentc); //comment
	mystr matrixsep = el; //may be changed to semicolon ...


	ElementData ed = parameter_EDCs(1)->Get(1);
	if ((parameter_EDCs(1)->Length()==1) && (! ed.IsEDC()))		// if it is not an EDC, than do not write "text = my text" but only "my text"
	{
		if (ed.IsText()) 
		{
			str = parameter_EDCs(1)->Get(1).GetText();
			str.Replace("\\n","\n");
		}
		else if (ed.IsBool())
		{
			str += mystr(ed.GetBool());
		}
		else if (ed.IsInt())
		{
			str += mystr(ed.GetInt());
		}
		else if (ed.IsDouble())
		{
			num.SmartDouble2String(ed.GetDouble());
			str += num;
		}
		else if (ed.IsText())
		{
			str += mystr(doublequote) + mystr(ed.GetText()) + mystr(doublequote);
		}
		else if (ed.IsVector())
		{
			str += osb;
			for (int j = 1; j <= ed.GetVectorLen(); j++)
			{
				num.SmartDouble2String(ed.GetVectorVal(j));
				str += num;
				if (j < ed.GetVectorLen()) str += comma;
			}
			str += csb;
		}
		else if (ed.IsMatrix())
		{
			str += osb;
			for (int j1 = 1; j1 <= ed.GetMatrixRows(); j1++)
			{
				for (int j2 = 1; j2 <= ed.GetMatrixCols(); j2++)
				{
					num.SmartDouble2String(ed.GetMatrixVal(j1, j2));
					str += num;
					if (j2 < ed.GetMatrixCols()) str += comma;
				}
				if (j1 < ed.GetMatrixRows()) str += matrixsep + mystr("    ");
			}
			str += csb;
		}
	}
	else
	{
		ElementDataContainer edc;
		edc.CopyFrom(*parameter_EDCs(1));
		edcParser->EDC2String(edc, str,"");
	}
	//mbs->UO(UO_LVL_warn) << str;
	mbs->UO(UO_LVL_0) << str;

	return 1;
}

int CEDC_ComInclude::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	if (!parameter_EDCs(1)->Get(1).IsText()) 
	{
		edcParser->EDCError(mystr("Error happened during parsing command ") + CommandName());
		return 0;
	}

	mystr param_text = parameter_EDCs(1)->Get(1).GetText(); // filename is absolute or alternatively relative to hid-file

	
	mystr filename = param_text; // param_text contains filename without double quotes

	int rv = 1; // success
	edcParser->MakeFilenameAbsolute(filename);		//$ DR 2013-07-04
	//$ MSax 2013-07-18 :[ begin
	//int pos = 1;
	//int isRelativePath = filename.PosPeek(pos) != ':';
	//if(isRelativePath)
	//{		
	//	mystr rel_path(""); 
	//	{
	//		mystr last_file_name = mbs->GetModelDataContainer()->TreeGetString("last_file_name");
	//		{
	//			int until;
	//			for(until=last_file_name.Length()-1;until>=0;until--)
	//			{
	//				if(last_file_name.PosPeek(until) == '\\')
	//					break; // until <=> last backslash
	//			}
	//			for(int i=0;i<=until;i++)
	//			{
	//				rel_path += last_file_name.PosPeek(i); // path until last backslash
	//			}
	//		}
	//	}
	//	
	//	// remove ..\ 
	//	{
	//		int pos=rel_path.Length()-1;
	//		if(rel_path.PosPeek(pos)=='\\')
	//		{
	//			rel_path.EraseChar(pos+1);
	//		}
	//	}

	//	int replaced = filename.Replace(mystr("..\\"), mystr("")); //find string searchstr, replace with string replacestr; this is done repeatedly; the return value counts the number of replacements
	//
	//	//for(int pos=rel_path.Length();pos--;pos>0)	//$ DR 2013-01-10 does not make sense
	//	for(int pos=rel_path.Length();pos>0;pos--)
	//	{
	//		if(!replaced){break;}
	//		int prev = pos-1;
	//		if(rel_path.PosPeek(prev)=='\\')
	//		{
	//			replaced --;
	//		}
	//		rel_path.EraseChar(pos);
	//	}
	//	filename = rel_path + mystr("\\") + filename;
	//}
	//$ MSax 2013-07-18 :] end
	if (DoesFileExist(filename))
	{
		ElementDataContainer edc_file;
		if(mbs->GetModelDataContainer() == 0)
		{
			rv = mbs->ReadModelData(filename); // empty model data container --> create new one
		}
		else
		{
			rv = mbs->AddModelData(filename);  // model data container has contents --> add and replace
		}
	}
	else
	{
		edcParser->EDCError(mystr("Error in Include - command: Can not find 'filename' '") + filename + mystr("'!"));
		rv = 0;
	}
	return rv;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2012-10
int CEDC_ComAddLoad::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int load_type_index = edc->Find("load_type"); //find index in element EDC
	if (load_type_index)
	{
		if (!edc->Get(load_type_index).IsText()) {assert(0);}

		//int listindex;
		mystr load_type_str = edc->Get(load_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCLoad, load_type_str);

		if (type_id < 0) 
		{
			edcParser->EDCError(mystr("Load type '") +  load_type_str + mystr("' not found in command AddLoad(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCLoad, type_id); //this are the typeflags, e.g. constraint, 2D, etc.

		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + load_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new load of selected type to MBS
			int load_elnum = edcParser->GetObjectFactory()->AddObject(OFCLoad, type_id); //Add Load with default values

			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetLoad(load_elnum).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetLoad(load_elnum).SetElementData(edc);

			return_value.SetInt(load_elnum, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'load_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2012-10
int CEDC_ComAddSensor::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int sensor_type_index = edc->Find("sensor_type"); //find index in element EDC
	if (sensor_type_index)
	{
		if (!edc->Get(sensor_type_index).IsText()) {assert(0);}

		//int listindex;
		mystr sensor_type_str = edc->Get(sensor_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCSensor, sensor_type_str);

		if (type_id < 0) 
		{
			edcParser->EDCError(mystr("Sensor type '") +  sensor_type_str + mystr("' not found in command AddSensor(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCSensor, type_id); //this are the typeflags, e.g. constraint, 2D, etc.

		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + sensor_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new sensor of selected type to MBS
			int sensor_elnum = edcParser->GetObjectFactory()->AddObject(OFCSensor, type_id); //Add Load with default values
			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetSensor(sensor_elnum).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetSensor(sensor_elnum).SetElementData(edc);

			return_value.SetInt(sensor_elnum, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'sensor_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2012-10
int CEDC_ComLoadSTL::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	if (!parameter_EDCs(1)->Get(1).IsText()) 
	{
		edcParser->EDCError(mystr("Error happened during parsing command ") + CommandName());
		return 0;
	}

	mystr param_text = parameter_EDCs(1)->Get(1).GetText(); // filename is absolute or alternatively relative to hid-file
	
	mystr filename = param_text; // param_text contains filename without double quotes

	int rv = 1; // success
	rv = edcParser->MakeFilenameAbsolute(filename);		//$ DR 2013-07-04
	//int pos = 1;
	//int isRelativePath = filename.PosPeek(pos) != ':';
	//if(isRelativePath)
	//{		
	//	mystr rel_path(""); 
	//	{
	//		mystr last_file_name = mbs->GetModelDataContainer()->TreeGetString("last_file_name");
	//		{
	//			int until;
	//			for(until=last_file_name.Length()-1;until>=0;until--)
	//			{
	//				if(last_file_name.PosPeek(until) == '\\')
	//					break; // until <=> last backslash
	//			}
	//			for(int i=0;i<=until;i++)
	//			{
	//				rel_path += last_file_name.PosPeek(i); // path until last backslash
	//			}
	//		}
	//	}
	//	
	//	// remove ..\ 
	//	{
	//		int pos=rel_path.Length()-1;
	//		if(rel_path.PosPeek(pos)=='\\')
	//		{
	//			rel_path.EraseChar(pos+1);
	//		}
	//	}

	//	int replaced = filename.Replace(mystr("..\\"), mystr("")); //find string searchstr, replace with string replacestr; this is done repeatedly; the return value counts the number of replacements
	//
	//	//for(int pos=rel_path.Length();pos--;pos>0)	//$ DR 2013-01-10 does not make sense
	//	for(int pos=rel_path.Length();pos>0;pos--)
	//	{
	//		if(!replaced){break;}
	//		int prev = pos-1;
	//		if(rel_path.PosPeek(prev)=='\\')
	//		{
	//			replaced --;
	//		}
	//		rel_path.EraseChar(pos);
	//	}
	//	filename = rel_path + mystr("\\") + filename;
	//}

	ElementDataContainer edc;

	if(rv)
	{
		rv = mbs->File2EDC(filename, &edc);
		return_value.SetEDC(&edc,"MeshData");
	}
	return rv;

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDC_ComAddMaterial::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int continuum_type_index = edc->Find("material_type"); //find index in element EDC
	mystr type = edc->TreeGetString("material_type");
	if (continuum_type_index)
	{
		if (!edc->Get(continuum_type_index).IsText()) {assert(0);}

		int listindex;
		mystr mat_type_str = edc->Get(continuum_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCMaterial, mat_type_str);

		if (type_id < 0)
		{
			edcParser->EDCError(mystr("Material type '") +  mat_type_str + mystr("' not found in command AddMaterial(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCMaterial, type_id); //this are the typeflags, e.g. constraint, 2D, etc.

		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + mat_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new material to MBS
			int mat_num = edcParser->GetObjectFactory()->AddObject(OFCMaterial, type_id); //Add Material with default values
			//retrieve the default EDC for the new material
			ElementDataContainer edc;
			mbs->GetMaterial(mat_num).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated material
			int rv = mbs->GetMaterial(mat_num).SetElementData(edc);

			return_value.SetInt(mat_num, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'material_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2013-01
int CEDC_ComAddBeam3DProperties::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int beam_property_type_index = edc->Find("material_type"); //find index in element EDC
	if (beam_property_type_index)
	{
		if (!edc->Get(beam_property_type_index).IsText()) {assert(0);}

		//int listindex;
		mystr beam_property_type_str = edc->Get(beam_property_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCBeamProperties, beam_property_type_str);

		if (type_id < 0) 
		{
			edcParser->EDCError(mystr("BeamProperty type '") +  beam_property_type_str + mystr("' not found in command AddBeamProperties(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}

		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCBeamProperties, type_id); //this are the typeflags, e.g. constraint, 2D, etc.
		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + beam_property_type_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new beam property of selected type to MBS
			int beamProp_matnum = edcParser->GetObjectFactory()->AddObject(OFCBeamProperties, type_id); //Add beam property with default values
			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetMaterial(beamProp_matnum).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetMaterial(beamProp_matnum).SetElementData(edc);

			return_value.SetInt(beamProp_matnum, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'beam_poperty_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2013-01
int CEDC_ComAddNode::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv_error = 0;
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	int node_type_index = edc->Find("node_type"); //find index in element EDC
	if (node_type_index)
	{
		if (!edc->Get(node_type_index).IsText()) {assert(0);}

		//int listindex;
		mystr node_property_str = edc->Get(node_type_index).GetText();
		int type_id = edcParser->GetObjectFactory()->GetObjectTypeId(OFCNode, node_property_str);

		if (type_id < 0) 
		{
			edcParser->EDCError(mystr("Node type '") +  node_property_str + mystr("' not found in command AddNode(...)!"));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		int typeflags = edcParser->GetObjectFactory()->GetTypeFlags(OFCNode, type_id); //this are the typeflags, e.g. constraint, 2D, etc.
		if( (typeflags & TAENotInRelease) && (edcParser->GetObjectFactory()->ExcludeExperimentalObjects()) )
		{
			edcParser->EDCError(mystr("ERROR: ") + node_property_str + mystr(" is not implemented yet. Check Homepage for available updates."));
			return_value.SetInt(rv_error, ""); //return value has no name!
			return 0;
		}
		else
		{
			//add new node of selected type to MBS
			int node_num = edcParser->GetObjectFactory()->AddObject(OFCNode, type_id); //Add Node with default values
			//retrieve the default EDC for the new element
			ElementDataContainer edc;
			mbs->GetNode(node_num).GetElementData(edc);
			//now overwrite the default with the user-defined values
			edc.TreeReplaceEDCDataWith(parameter_EDCs(1),1);
			//now set the user-defined EDC in the newly generated element
			int rv = mbs->GetNode(node_num).SetElementData(edc);

			return_value.SetInt(node_num, ""); //return value has no name!

			return 1;
		}
	}
	else
	{
		edcParser->EDCError(mystr("did not find 'node_type' in command ") + CommandName());
		return_value.SetInt(rv_error, ""); //return value has no name!
		return 0;
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2013-01-30
int CEDC_ComComputeInertia::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* edc = parameter_EDCs(1); //this contains the element Data
	ElementDataContainer edc_rv;
	int rv = mbs->ComputeInertia(edc,&edc_rv);
	if(rv)
	{
		return_value.SetEDC(&edc_rv,"InertiaValues");
	}
	return rv;
}
mystr CEDC_ComComputeInertia::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameter of this command is an ElementDataContainer, with the following entries: \n");
	descr +=			mystr("\\begin{itemize} \n");
	descr +=			mystr("  \\item density or material\\_number (one of these 2 has to be set!) \n");
	descr +=			mystr("	 \\item One of the following options to define the geometry: \n");
	descr +=			mystr("		  \\begin{itemize} \n");
	descr +=			mystr("				\\item  MeshData.triangles and MeshData.points \\\\ \n");
	descr +=			mystr("							  both entries are Matrices with 3 columns \n");
	descr +=			mystr("				\\item  Cube.body\\_dimensions \n");
	descr +=			mystr("			\\end{itemize} \n");
	descr +=			mystr("\\end{itemize} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2013-02-19

mystr CEDC_ComTransformPoints::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are as follows \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("  \\item points: Matrix of the points: Each line represents a point p. The 3 columns are the x-, y- and z-coordinate \n");
	descr +=			mystr("	 \\item trans: Vector of translation, 3 dimensions! \n");
	descr +=			mystr("	 \\item rot: rotation matrix (3x3), can be used for scaling as well as rotation \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};

int CEDC_ComTransformPoints::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* edcpoints = parameter_EDCs(1);	//the points that should be transformed
	ElementDataContainer* edcvec = parameter_EDCs(2);			//translation
	ElementDataContainer* edcMat = parameter_EDCs(3);			//rotation matrix
	
	int nRows = 0;
	int nCols = 3;
	double * p;

	ElementData* edp = edcMat->GetPtr(1);
	p = edp->GetMatrix();
	nRows = edp->GetMatrixRows();
	Matrix rotation(nRows,nCols,p);
	
	edp = edcpoints->GetPtr(1);
	p = edp->GetMatrix();
	nRows = edp->GetMatrixRows();
	Matrix points(nRows,nCols,p);

	Matrix new_points;
	new_points.SetSize(nRows,nCols);

	Vector3D translation;
	edp = edcvec->GetPtr(1);
	edp->GetVector(translation.X(),translation.Y(),translation.Z());

	Vector3D help;
	for (int i=1; i<=nRows; i++)
	{
		help = Vector3D(points(i,1),points(i,2),points(i,3));
		help = translation + rotation*help;
		new_points(i,1)=help.X();
		new_points(i,2)=help.Y();
		new_points(i,3)=help.Z();
	}

	return_value.SetMatrix(new_points.GetMatPtr(),new_points.Getrows(),new_points.Getcols(),"points");

	int rv = 1;
	return rv;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ DR 2013-07-04
int CEDC_ComLoadVecFromFile::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* edc = parameter_EDCs(1); 
	mystr filename;
	if(edc->Get(1).IsText())
	{
		filename = edc->Get(1).GetText();
	}
	else
	{
		filename = edc->TreeGetString("parameter1");
	}
	edcParser->MakeFilenameAbsolute(filename);

	edc = parameter_EDCs(2);
	mystr tmp;
	int col = 0;
	if(edc->Get(1).IsDouble())
	{
		col = edc->Get(1).GetDouble();
	}
	else if(edc->Get(1).IsInt())
	{
		col = edc->Get(1).GetInt();
	}
	else
	{
		tmp = edc->TreeGetString("parameter2");
	}

	if(tmp.IsValidNumber())
	{
		col = tmp.MakeInt();
	}
	else
	{
		if(col<=0)
		{
			edcParser->EDCError("Error while executing command 'LoadVectorFromFile': 2nd parameter has to be an integer!\n");
		}
	}

	ElementDataContainer edc_rv;
	int rv = mbs->LoadVectorFromFile(filename, col, &edc_rv);
	if(rv)
	{
		double* vptr;
		int len;
		edc_rv.TreeGetVector("Vector",&vptr,len);
		return_value.SetVector(vptr,len,"vector");
	}
	else
	{
		edcParser->EDCError("Error while executing command 'LoadVectorFromFile': could not read file!\n");
	}
	return rv;
}

mystr CEDC_ComLoadVecFromFile::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item The name of the file as string\n");
	descr +=			mystr("  \\item An integer defining in which column of the file the vector is stored\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-07-11
int CEDC_ComSum::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* summand1 = parameter_EDCs(1);   // first summand
	ElementDataContainer* summand2 = parameter_EDCs(2);   // second summand
	ElementData *edp1, *edp2;

	edp1 = summand1->GetPtr(1);
	edp2 = summand2->GetPtr(1);

	int rv = 1;
// scalar+scalar
	if (edp1->IsDouble() && edp2->IsDouble())
	{
		double s1 = edp1->GetDouble();
		double s2 = edp2->GetDouble();
		return_value.SetDouble(s1+s2,"sum");
	}
// vector+vector 
	else if (edp1->IsVector() && edp2->IsVector())
	{
		int l1 = edp1->GetVectorLen();
		int l2 = edp2->GetVectorLen();
		if(l1==l2)
		{
			Vector v1,v2;
			v1.LinkWith(edp1->GetVector(),l1);
			v2.LinkWith(edp2->GetVector(),l2);

			//+++++++++++++++++++++++++++++++++++++++++++
			Vector res = v1+v2; //internal representation of command
			//+++++++++++++++++++++++++++++++++++++++++++

			return_value.SetVector(res.GetVecPtr(), res.GetLen(), "sum");
		}
		else
		{
			edcParser->EDCError("Error while executing command 'Add': adding two vectors of different length\n");
		}
	}
// matrix + matrix ( this includes Vector Transposed )
	else if (edp1->IsMatrix() && edp2->IsMatrix())
	{
		int r1 = edp1->GetMatrixRows();
		int r2 = edp2->GetMatrixRows();
		int c1 = edp1->GetMatrixCols();
		int c2 = edp2->GetMatrixCols();

		if(r1==r2 && c1==c2)
		{
			Matrix m1(r1, c1, edp1->GetMatrix());
			Matrix m2(r2, c2, edp2->GetMatrix());

			//+++++++++++++++++++++++++++++++++++++++++++
			Matrix res = m1+m2; //internal representation of command
			//+++++++++++++++++++++++++++++++++++++++++++

			return_value.SetMatrix(res.GetMatPtr(), res.Getrows(), res.Getcols(), "sum");
		}
		else
		{
			edcParser->EDCError("Error while executing command 'Add': adding two matrixes of different length\n");
		}
	}
	else
	{
		edcParser->EDCError("Error while executing command 'Add': cannot add different types of scalar/vector/matrix\n");
	}
	return rv;
}

mystr CEDC_ComSum::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ summand, either scalar, vector or matrix \n");
	descr +=			mystr("  \\item $2^{nd}$ summand, either scalar, vector or matrix \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-07-11
int CEDC_ComProduct::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* summand1 = parameter_EDCs(1);   // first summand
	ElementDataContainer* summand2 = parameter_EDCs(2);   // second summand
	ElementData *edp1, *edp2;

	// when parameters are defined inline the parameter EDC contains text.
	if(summand1->GetPtr(1)->IsText())
		EvaluateTextualParameterEntry(mbs, edcParser, summand1);
	if(summand2->GetPtr(1)->IsText())
		EvaluateTextualParameterEntry(mbs, edcParser, summand2);

	edp1 = summand1->GetPtr(1);
	edp2 = summand2->GetPtr(1);


	int rv = 1;
// involving a scalar: {s*s, s*v, s*m, v*s, m*s}
	if (edp1->IsDouble() || edp2->IsDouble())
	{
		double s;
		ElementData *any;
		if (edp1->IsDouble())
		{
			s = edp1->GetDouble();
			any = edp2;
		}
		else
		{
			s = edp2->GetDouble();
			any = edp1;
		}
		
		if(any->IsDouble()) //s*s
		{
			double r = any->GetDouble();
			return_value.SetDouble(s*r,"product");
		}
		else //s*v, s*m,
		{
			if(any->IsVector())
			{
				Vector V;
				V.LinkWith(any->GetVector(),any->GetVectorLen());

				//+++++++++++++++++++++++++++++++++++++++++++
				Vector res = V*s; //internal representation of command
				//+++++++++++++++++++++++++++++++++++++++++++

				return_value.SetVector(res.GetVecPtr(), res.GetLen(), "product");
			}
			else
			{
				Matrix M(any->GetMatrixRows(), any->GetMatrixCols(), any->GetMatrix());

				//+++++++++++++++++++++++++++++++++++++++++++
				Matrix res = M*s; //internal representation of command
				//+++++++++++++++++++++++++++++++++++++++++++

				return_value.SetMatrix(res.GetMatPtr(), res.Getrows(), res.Getcols(), "product");
			}
		}
	}
// involving only vector and matrix: {v*v, v*m, m*v, m*m}
	else
	{
		Vector V1,V2;
		Matrix M1,M2;

//left factor
		if(edp1->IsVector())
		{
			V1.LinkWith(edp1->GetVector(), edp1->GetVectorLen());
		}
		else if(edp1->IsMatrix())
		{
			M1 = Matrix(edp1->GetMatrixRows(),edp1->GetMatrixCols(),edp1->GetMatrix());
		}
		else
		{
			edcParser->EDCError("Error while executing command 'Product': first factor neither scalar nor vector nor matrix \n");
		}
// right factor
		if(edp2->IsVector())
		{
			V2.LinkWith(edp2->GetVector(), edp2->GetVectorLen());
		}
		else if(edp2->IsMatrix())
		{
			M2 = Matrix(edp2->GetMatrixRows(),edp2->GetMatrixCols(),edp2->GetMatrix());
		}
		else
		{
			edcParser->EDCError("Error while executing command 'Product': second factor neither scalar nor vector nor matrix \n");
		}

// compute the product
// case m*m = m --> matrix product
		if(edp1->IsMatrix()&&edp2->IsMatrix()) 
		{
			if(M1.Getcols() == M2.Getrows())
			{
				//+++++++++++++++++++++++++++++++++++++++++++
				//internal representation of command
				Matrix res;
				Mult(M1,M2,res);
				//+++++++++++++++++++++++++++++++++++++++++++
				return_value.SetMatrix(res.GetMatPtr(), res.Getrows(), res.Getcols(), "product");
			}
			else
			{
				edcParser->EDCError("Error while executing command 'Product': dimension check failed - number of columns in left factor must be equal to number of rows in right factor \n");
			}
		}
// case v*m = v --> matrix product
		else if(edp1->IsVector()&&edp2->IsMatrix())
		{
			if(V1.GetLen() == M2.Getrows())
			{
				//+++++++++++++++++++++++++++++++++++++++++++
				//internal representation of command
				Vector res;
				Mult(V1,M2,res);
				//+++++++++++++++++++++++++++++++++++++++++++
				return_value.SetVector(res.GetVecPtr(), res.GetLen(), "product");
			}
			else
			{
				edcParser->EDCError("Error while executing command 'Product': dimension check failed - number of columns in left factor must be equal to number of rows in right factor \n");
			}
		}

//! AD: !!! MAYBE TAKE THESE CASES OUT !!! more commands eg transpose(), scalarproduct()
// case v*v =s -> scalar product
		else if(edp1->IsVector()&&edp2->IsVector())
		{
			//$ AD: warning removed, now in docu 
			//edcParser->EDCError("WARNING: multiplying two (row-)vectors is autmatically assumed to compute the scalar product of the two vectors \n");
			if(V1.GetLen() == V2.GetLen())
			{
				double scalarproduct = 0;
				for(int i=1; i<=V1.GetLen(); i++)
				{
					scalarproduct += ( V1(i)*V2(i) );
				}
				return_value.SetDouble(scalarproduct,"product");	
			}
			else
			{
				edcParser->EDCError("Error while executing command 'Product': factors (vector) must have same length \n");
			}
		}
// case m*v --> m*vT=vT
		else if(edp1->IsMatrix()&&edp2->IsVector())
		{
			//$ AD: warning removed, now in docu 
			//edcParser->EDCError("WARNING: multiplying a matrix with a (row-)vector is autmatically assumed to compute the product with the transposed vector \n");
			if(M1.Getcols() == V2.GetLen())
			{
				//+++++++++++++++++++++++++++++++++++++++++++
				//internal representation of command
				Vector res;
				Mult(M1,V2,res);
				//+++++++++++++++++++++++++++++++++++++++++++
				return_value.SetMatrix(res.GetVecPtr(), res.GetLen(), 1, "product");
			}
			else
			{
				edcParser->EDCError("Error while executing command 'Product': dimension check failed - number of columns in left factor must be equal to number of rows in right factor \n");
			}
		}
		else
		{
			edcParser->EDCError("Error while executing command 'Product': dimension check failed - number of columns in left factor must be equal to number of rows in right factor \n");
		}
	}
	return rv;
}

mystr CEDC_ComProduct::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ factor, either scalar, vector or matrix \n");
	descr +=			mystr("  \\item $2^{nd}$ factor, either scalar, vector or matrix \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("  product of two vectors is always computed as scalar product \n");
	descr +=			mystr("  for vector times Matrix the vector is automatically transposed if required \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-07-11
int CEDC_ComTranspose::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* edc1 = parameter_EDCs(1);   // sole input
	ElementData *edp1;
	edp1 = edc1->GetPtr(1);

	int rv = 1;
// transpose a vector -> vT is a matrix
	if (edp1->IsVector())
	{
		Vector V;
		V.LinkWith(edp1->GetVector(), edp1->GetVectorLen());;
		//+++++++++++++++++++++++++++++++++++++++++++
		//internal representation of command
		Matrix res(V.GetLen(), 1, V.GetVecPtr());
		//+++++++++++++++++++++++++++++++++++++++++++
		return_value.SetMatrix(res.GetMatPtr(), res.Getrows(), res.Getcols(), "transposed");
	}
	else if (edp1->IsMatrix())
	{
		Matrix M(edp1->GetMatrixRows(),edp1->GetMatrixCols(),edp1->GetMatrix());
		//+++++++++++++++++++++++++++++++++++++++++++
		//internal representation of command
		Matrix res = M.GetTp();
		//+++++++++++++++++++++++++++++++++++++++++++
		if(res.Getrows() == 1)
		{
		// this is actually a vector ...
			return_value.SetVector(res.GetMatPtr(), res.Getcols(), "product");
		}
		else
		{
			return_value.SetMatrix(res.GetMatPtr(), res.Getrows(), res.Getcols(), "product");
		}
	}
	else
	{
		edcParser->EDCError("Error while executing command 'Transpose': only vectors and matrices can be transposed\n");
	}
	return rv;
}

mystr CEDC_ComTranspose::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item vector or matrix to be transposed\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-07-11
//$! AD 2013-07-15 commented out Class Scalar product for the time being, functionality is in ComProduct...

////int CEDC_ComScalarProduct::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
////{
////	ElementDataContainer* summand1 = parameter_EDCs(1);   // first summand
////	ElementDataContainer* summand2 = parameter_EDCs(2);   // second summand
////	ElementData *edp1, *edp2;
////
////	edp1 = summand1->GetPtr(1);
////	edp2 = summand2->GetPtr(1);
////
////	int rv = 1;
////	if(edp1->IsVector()&&edp2->IsVector())
////	{	
////		Vector V1,V2;
////		V1.LinkWith(edp1->GetVector(), edp1->GetVectorLen());
////		V2.LinkWith(edp2->GetVector(), edp2->GetVectorLen());
////
////		if(V1.GetLen() == V2.GetLen())
////		{
////			double scalarproduct = 0;
////			for(int i=1; i<=V1.GetLen(); i++)
////			{
////				scalarproduct += ( V1(i)*V2(i) );
////			}
////			return_value.SetDouble(scalarproduct,"product");	
////		}
////		else
////		{
////			edcParser->EDCError("Error while executing command 'ScalarProduct': vectors must have same length \n");
////		}
////	}
////	else
////	{
////		edcParser->EDCError("Error while executing command 'ScalarProduct': scalar product only defined for two vectors\n");
////	}
////	return rv;
////}
////
////mystr CEDC_ComScalarProduct::GetCommandTexParameterDescription() const 
////{
////	mystr descr = mystr("The parameters of this command are \n");
////	descr +=			mystr("\\begin{enumerate} \n");
////	descr +=			mystr("	 \\item 1^{st} vector \n");
////	descr +=			mystr("  \\item 2^{nd} vector \n");
////	descr +=			mystr("\\end{enumerate} \n");
////	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
////
////	return descr;
////};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-07-11
int CEDC_ComCrossProduct::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementDataContainer* summand1 = parameter_EDCs(1);   // first summand
	ElementDataContainer* summand2 = parameter_EDCs(2);   // second summand
	ElementData *edp1, *edp2;

	edp1 = summand1->GetPtr(1);
	edp2 = summand2->GetPtr(1);

	int rv = 1;
	if(edp1->IsVector()&&edp2->IsVector())
	{
		Vector V1,V2;
		V1.LinkWith(edp1->GetVector(), edp1->GetVectorLen());
		V2.LinkWith(edp2->GetVector(), edp2->GetVectorLen());

		if(V1.GetLen() == V2.GetLen() && (V1.GetLen() == 3 || V1.GetLen() == 2) )   // allow  3Dx3D and 2Dx2D
		{
			Vector3D v1,v2;
			if(V1.GetLen() == 2)
			{
				v1 = Vector3D(V1(1), V1(2), 0.);
				v2 = Vector3D(V2(1), V2(2), 0.);
			}
			else 
			{
				v1 = Vector3D(V1(1), V1(2), V1(3));
				v2 = Vector3D(V2(1), V2(2), V2(3));
			}
						
			//+++++++++++++++++++++++++++++++++++++++++++
			//internal representation of command
			Vector3D res = v1.Cross(v2);
			//+++++++++++++++++++++++++++++++++++++++++++
			if(V1.GetLen() == 3)
			{
				return_value.SetVector(res.GetVecPtr(), 3, "crossproduct");	// return 3D Vector  3Dx3D -> 3D
			}
			else if(V1.GetLen() == 2)
			{
				return_value.SetDouble(res(3), "crossproduct");	// return scalar 2Dx2D -> scalar
			}
		}
		else
		{
			edcParser->EDCError("Error while executing command 'CrossProduct': vectors must have same length (2D or 3D)\n");
		}
	}
	else
	{
		edcParser->EDCError("Error while executing command 'CrossProduct': scalar product only defined for two vectors\n");
	}
	return rv;
}

mystr CEDC_ComCrossProduct::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ vector (2D or 3D)\n");
	descr +=			mystr("  \\item $2^{nd}$ vector (2D or 3D)\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("for two 3D vectors the retuen value is also a 3D vector. For two 2D vectors the return value is a scalar.\n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex

	return descr;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-09-13
int CEDC_ComFor::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;

	ElementData* ED_loop_init = parameter_EDCs(1)->GetPtr(1);
  ElementData* ED_loop_cond = parameter_EDCs(2)->GetPtr(1);
	ElementData* ED_loop_incr = parameter_EDCs(3)->GetPtr(1);
  ElementData* ED_loop_code = parameter_EDCs(4)->GetPtr(1);

	int err = 0;
	int suc = 1;

// CHECK PARAMETERS:
// parameter 1: loop initialization
	mystr str_loop_init;
	if(ED_loop_init->IsText())
	{
		str_loop_init = ED_loop_init->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing FOR loop: first parameter is not a text (loop initialization)");
		return rv;
	}
// parameter 2: loop condition
	mystr str_loop_cond;
	if(ED_loop_cond->IsText())
	{
		str_loop_cond = ED_loop_cond->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing FOR loop: second parameter is not a text (loop condition)");
		return rv;
	}
// parameter 3: loop increment
	mystr str_loop_incr;
	if(ED_loop_incr->IsText())
	{
		str_loop_incr = ED_loop_incr->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing FOR loop: third parameter is not a text (loop increment)");
		return rv;
	}

// parameter 4: loop code
	mystr str_loop_code;
	if(ED_loop_code->IsText())
	{
		str_loop_code = ED_loop_code->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing FOR loop: loop code is to text (text in brackets required)");
		return rv;
	}

// the actual FOR Loop

// initialize local loop variable container
	ElementDataContainer& localEDC = *(edcParser->GetMBSParser()->GetLocalEDC());  // existing local EDC for the FOR command ("outside" the FOR command)

	// initialize the loop variable
	suc = edcParser->ParseAndExecuteString(str_loop_init, localEDC);
	if(!suc)
	{
		edcParser->EDCError("Error while executing FOR loop: first parameter could not be translated to a loop initialization statement (\"i=1\")");
		return rv;
	}


// loop condition 1st time
	int flag_loop_condition = (int) (edcParser->GetMBSParser()->ExpressionToDouble(str_loop_cond,err)+0.5);				// loop condition

	while(flag_loop_condition) 
	{
// execute loop code 
		suc = edcParser->ParseAndExecuteString(str_loop_code, localEDC); 
		if(!suc)
		{
			edcParser->EDCError("Error while executing FOR loop: failed to execute loop code)");
			return rv;
		}
		rv++;

// increment
		suc = edcParser->ParseAndExecuteString(str_loop_incr, localEDC);

// condition
		flag_loop_condition = (int) (edcParser->GetMBSParser()->ExpressionToDouble(str_loop_cond,err)+0.5);
	}

	return rv;
}

mystr CEDC_ComFor::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ define and initialize loop variable (\"i=1\")\n");
	descr +=			mystr("  \\item $2^{nd}$ loop condition (\"i<5\")\n");
	descr +=			mystr("  \\item $3^{rd}$ loop increment (\"i=i+1\")\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("the command must be followed by a container for the loop code \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int CEDC_ComIf::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;

	ElementData* ED_loop_cond = parameter_EDCs(1)->GetPtr(1);
	ElementData* ED_loop_code = parameter_EDCs(2)->GetPtr(1);

	int err = 0;
	int suc = 1;

// CHECK PARAMETERS:
// parameter 1: loop condition
	mystr str_loop_cond;
	if(ED_loop_cond->IsText())
	{
		str_loop_cond = ED_loop_cond->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing IF condition: second parameter is not a text (loop condition)");
		return rv;
	}
// parameter 2: loop code
	mystr str_loop_code;
	if(ED_loop_code->IsText())
	{
		str_loop_code = ED_loop_code->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing IF condition: loop code is to text (text in brackets required)");
		return rv;
	}

// the actual IF condition

// initialize local loop variable container
	ElementDataContainer& localEDC = *(edcParser->GetMBSParser()->GetLocalEDC());  // existing local EDC for the IF command ("outside" the IF command)

// loop condition 1st time
	int flag_loop_condition = (int) (edcParser->GetMBSParser()->ExpressionToDouble(str_loop_cond,err)+0.5);				// loop condition

	if(flag_loop_condition)
	{
// execute conditional code 
		suc = edcParser->ParseAndExecuteString(str_loop_code, localEDC); 
	}
	return flag_loop_condition;
	/*return rv;*/
}

mystr CEDC_ComIf::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ condition (\"i<10\")\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("the command must be followed by a container for the conditional code \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ AD 2013-10-09
int CEDC_ComMesh_GenerateNewMesh::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
// MESH PROPERTIES
// Default values for a MeshHandle
// * this is to be autogenerated at a later time
	ElementDataContainer edc_meshproperties;
	ElementData ed;
	// why cant bool be locked ?
	ed.SetBool(1,"isHandle");												ed.SetLocked(1);	ed.SetToolTipText("This Container represents a Handle. Commands may be executed for this object");	
	edc_meshproperties.Add(ed);
	ed.SetText("generates discretized Beams","GenerateBeam");			ed.SetLocked(1);	ed.SetToolTipText("This Command is registerd for this handle");
	edc_meshproperties.TreeAdd("Functions",ed);
	ed.SetText("generates discretized Plates","GeneratePlate");		ed.SetLocked(1);	ed.SetToolTipText("This Command is registerd for this handle");
	edc_meshproperties.TreeAdd("Functions",ed);
	ed.SetText("remove redundant nodes or create constraints","GlueMesh");	ed.SetLocked(1);	ed.SetToolTipText("This Command is registerd for this handle");
	edc_meshproperties.TreeAdd("Functions",ed);
	ed.SetText("returns nodes in the defined box","GetNodesInBox");	ed.SetLocked(1);	ed.SetToolTipText("This Command is registerd for this handle");
	edc_meshproperties.TreeAdd("Functions",ed);

// Overwrite defaults with input from container 1
	ElementDataContainer* edc_override = parameter_EDCs(1);
	edc_meshproperties.TreeReplaceEDCDataWith(edc_override);

// TODO for script compatible mesh class
	// add a new MeshElement to the MBS

// MESH OBJECT AND VARIABLES
// * this is to be replaced by members of the script compatible mesh class
	double dummy = 0;
	ed.SetVector(&dummy, 0, "list_of_nodes");			ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Mesh or empty");			
	edc_meshproperties.Add(ed);
	ed.SetVector(&dummy, 0, "list_of_elements");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Elements associated with the Mesh or empty");
	edc_meshproperties.Add(ed);
	ed.SetVector(&dummy, 0, "list_of_redundant_nodes");	ed.SetValuesInt();	ed.SetVariableLength(); ed.SetToolTipText("List of redundant nodes in the mesh");
	edc_meshproperties.Add(ed);

	return_value.SetEDC(&edc_meshproperties,"MeshProperties");

	return 1;
}


mystr CEDC_ComMesh_GenerateNewMesh::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter EDC to overwrite the default properties\n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("the return vaule of the command MUST be assigned to a new variable(handle) \n\n");
	descr +=			mystr("overwritable enties in the properties EDC are: \n");
	descr +=			mystr("\\begin{itemize} \n");	
	descr +=			mystr("	 \\item ... \n");
	descr +=			mystr("\\end{itemize} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDC_ComMesh_GenerateBeam::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
	ElementDataContainer edc_rv;

// Default values for a BeamParameters
// * this is to be autogenerated at a later time
	ElementDataContainer edc_beamparameters;
	ElementData ed;
	ed.SetVector3D(0.,0.,0.,"P1");							edc_beamparameters.Add(ed);
	ed.SetVector3D(1.,0.,0.,"P2");							edc_beamparameters.Add(ed);
	ed.SetInt(1,"matnr");												edc_beamparameters.Add(ed);
	ed.SetInt(1,"discretization");							edc_beamparameters.Add(ed);
	ed.SetText("LinearBeam3D","element_type");	edc_beamparameters.Add(ed);
	ed.SetText("Node3DRxyz","node_type");				edc_beamparameters.Add(ed);

// Overwrite defaults with input from container 1
	ElementDataContainer* edc_override = parameter_EDCs(1);
	edc_beamparameters.TreeReplaceEDCDataWith(edc_override);

// container 2 contains the name of the handle ( can not pass the EDC* itself because that is going to be deleted... )
// consistency already checked in ParseCommandParameters...
	mystr handlename = parameter_EDCs(2)->TreeGetString("handle");
	ElementDataContainer* edc_local = edcParser->GetMBSParser()->GetLocalEDC(); 
	int nr_handle = edc_local->Find(handlename);
	ElementDataContainer* edc_handle_mesh = edc_local->Get(nr_handle).GetEDC();


// Actual Generation of BeamElements - to be a function of the MeshElement
// read geometry parameters
	Vector3D p1,p2;
	edc_beamparameters.TreeGetVector3D("P1",p1.X(), p1.Y(), p1.Z());
	edc_beamparameters.TreeGetVector3D("P2",p2.X(), p2.Y(), p2.Z());
	int matnr = edc_beamparameters.TreeGetInt("matnr");
	int discr = edc_beamparameters.TreeGetInt("discretization");

// when more then one type of element is implemented, perform a type check 
	mystr elem_type = edc_beamparameters.TreeGetString("element_type");
	mystr node_type = edc_beamparameters.TreeGetString("node_type");

// get current lists for mesh object
	double* vdata;
	int vlen;
	edc_handle_mesh->TreeGetVector("list_of_nodes",&vdata,vlen);            
	TArray<double> mesh_nnrs; mesh_nnrs.CopyFrom(vdata,vlen); mesh_nnrs.SetLen(vlen);			// node numbers already associated with the Mesh
	edc_handle_mesh->TreeGetVector("list_of_elements",&vdata,vlen);
	TArray<double> mesh_enrs; mesh_enrs.CopyFrom(vdata,vlen); mesh_enrs.SetLen(vlen);			// element numbers already associated with the Mesh

// *******************
// GENERATION OF NODES

// use existing script command AddNode
	mystr commandname("AddNode");
	int nparam, hasretval, c_option;
	CEDC_Command* command_addnode = edcParser->GetCommand( edcParser->GetCommandType(commandname, nparam, hasretval, c_option) );

// prepare EDC to be able to use existing script commands 
	ElementDataContainer edc_nodeproperties;
	TArray<ElementDataContainer*> paramEDCforNode;
	paramEDCforNode.Add(&edc_nodeproperties);
	ed.SetText(node_type,"node_type");																		edc_nodeproperties.Add(ed);

//  todo: automatically compute orientation for all nodes - depending on node_type ...
	// Node3DRxyz 
	Vector3D p12 = p2-p1;
	double pyth_xy= pow((p12.X()*p12.X()+p12.Y()*p12.Y()),(.5)); 

	double roll = 0;
	double pitch = atan2(p12.Y(), p12.X());
	double yaw = atan2(p12.Z(), pyth_xy);
	ed.SetVector3D(roll,-yaw,pitch,"reference_rot_angles");											edc_nodeproperties.TreeAdd("Geometry",ed);


// override discretization if entry element_size is present
	int j = edc_beamparameters.Find("element_size");											
	if (j != 0 && edc_beamparameters.Get(j).IsDouble())
	{
		double elemsize = edc_beamparameters.TreeGetDouble("element_size");
		double total_length = p12.Norm();
		double required = total_length / elemsize;
		int discr_new = (int) ceil(required);
		discr = discr_new;
	}


	TArray<double> nnrs;                                                  // node numbers for the current call GenerateBeam
	for(int i=0; i<=discr; i++)
	{
		Vector3D pos = ( (discr-i)*p1 + i*p2 ) * (1./double(discr));
		ed.SetVector3D(pos.X(),pos.Y(),pos.Z(),"reference_position");				edc_nodeproperties.TreeAdd("Geometry",ed);

// call AddNode with prepared EDC, returnvalue is number of added node
		command_addnode->ExecuteCommand(mbs,edcParser,paramEDCforNode,ed,0);
		
		nnrs.Add(ed.GetInt());
	}
// return value:nodes from current generation function
	ed.SetVector( nnrs.GetDataPtr(), nnrs.Length(), "list_of_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Beam or empty");
	edc_rv.Add(ed);
// update list_of_nodes in Mesh
	mesh_nnrs.Merge(nnrs);
	ed.SetVector( mesh_nnrs.GetDataPtr(), mesh_nnrs.Length(), "list_of_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Mesh or empty");
	ElementDataContainer update_nnrs; update_nnrs.Add(ed);
	edc_handle_mesh->TreeReplaceEDCDataWith(&update_nnrs);

// ***************************
// GENERATION OF BEAM ELEMENTS

// use existing script command AddElement
	commandname = mystr("AddElement");
	CEDC_Command* command_addelement = edcParser->GetCommand( edcParser->GetCommandType(commandname, nparam, hasretval, c_option) );
	
	// prepare EDC to be able to use existing script commands 
	ElementDataContainer edc_elementproperties;
	ed.SetText(elem_type,"element_type");																	edc_elementproperties.Add(ed);
	ed.SetInt(matnr,"material_number");																		edc_elementproperties.TreeAdd("Physics",ed);
	TArray<ElementDataContainer*> paramEDCforElement;
	paramEDCforElement.Add(&edc_elementproperties);

	TArray<double> enrs;
	for(int i=1; i<=discr; i++)
	{
		ed.SetInt(nnrs(i),"node_1");																				edc_elementproperties.TreeAdd("Geometry",ed);
		ed.SetInt(nnrs(i+1),"node_2");																			edc_elementproperties.TreeAdd("Geometry",ed);

// call AddElement with prepared EDC, returnvalue is number of added element
		command_addelement->ExecuteCommand(mbs,edcParser,paramEDCforElement,ed,0);

		enrs.Add(ed.GetInt());
	}

// return value:nodes from current generation function
	ed.SetVector( enrs.GetDataPtr(), enrs.Length(), "list_of_elements");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Elements associated with the Beam or empty");
	edc_rv.Add(ed);
// update list_of_nodes in Mesh
	mesh_enrs.Merge(enrs);
	ed.SetVector( mesh_enrs.GetDataPtr(), mesh_enrs.Length(), "list_of_elements");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Elements associated with the Mesh or empty");
	ElementDataContainer update_enrs; update_enrs.Add(ed);
	edc_handle_mesh->TreeReplaceEDCDataWith(&update_enrs);


	return_value.SetEDC(edc_rv.GetCopy(), "RV_Beam");

	return 1;
}

mystr CEDC_ComMesh_GenerateBeam::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter EDC containing the beam properties \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("enties in the properties EDC are: \n");
	descr +=			mystr("\\begin{itemize} \n");	
	descr +=			mystr("	 \\item P1 - position of left outer node \n");
	descr +=			mystr("	 \\item P2 - position of right outer node \n");
	descr +=      mystr("  \\item matnr - number of the material to be used for the beam elements (Beam3DProperties) \n");
	descr +=      mystr("  \\item discretization - number of the beam elements \n");
	descr +=      mystr("  \\item elementsize - (optional) maximum element size ( if set has priority over discretization ) \n");
	descr +=			mystr("\\end{itemize} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};


int CEDC_ComMesh_GeneratePlate::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
	ElementDataContainer edc_rv;

	// Default values for a PlateParameters
// * this is to be autogenerated at a later time
	ElementDataContainer edc_plateparameters;
	ElementData ed;
	ed.SetVector3D(0.,0.,0.,"P1");				edc_plateparameters.Add(ed);
	ed.SetVector3D(1.,0.,0.,"P2");				edc_plateparameters.Add(ed);
	ed.SetVector3D(0.,1.,0.,"P3");				edc_plateparameters.Add(ed);
	ed.SetInt(1,"matnr");						edc_plateparameters.Add(ed);
	ed.SetVector2D(1.,1.,"discretization");		ed.SetValuesInt();	edc_plateparameters.Add(ed);
	
	ed.SetText("ANCFThinPlate3D","element_type");	edc_plateparameters.Add(ed);
	ed.SetText("Node3DS1S2","node_type");			edc_plateparameters.Add(ed);
	ed.SetBool(1, "nonlinear");		edc_plateparameters.Add(ed);

// Overwrite defaults with input from container 1
	ElementDataContainer* edc_override = parameter_EDCs(1);
	edc_plateparameters.TreeReplaceEDCDataWith(edc_override);

// container 2 contains the name of the handle ( can not pass the EDC* itself because that is going to be deleted... )
// consistency already checked in ParseCommandParameters...
	mystr handlename = parameter_EDCs(2)->TreeGetString("handle");
	ElementDataContainer* edc_local = edcParser->GetMBSParser()->GetLocalEDC(); 
	int nr_handle = edc_local->Find(handlename);
	ElementDataContainer* edc_handle_mesh = edc_local->Get(nr_handle).GetEDC();


// Actual Generation of PlateElements - to be a function of the MeshElement
// read geometry parameters 
	Vector3D p1,p2,p3;
	edc_plateparameters.TreeGetVector3D("P1",p1.X(), p1.Y(), p1.Z());
	edc_plateparameters.TreeGetVector3D("P2",p2.X(), p2.Y(), p2.Z());
	edc_plateparameters.TreeGetVector3D("P3",p3.X(), p3.Y(), p3.Z());
	int matnr = edc_plateparameters.TreeGetInt("matnr");
	double thickness = edc_plateparameters.TreeGetDouble("thickness");
	Vector2D discr;
	edc_plateparameters.TreeGetVector2D("discretization", discr.X(), discr.Y());
	int is_nonlinear = edc_plateparameters.TreeGetBool("nonlinear", 1);
	
	TArray<double> loads;
	int i = edc_plateparameters.Find("load"); 
	if (i != 0)
	{
		int load = edc_plateparameters.TreeGetInt("load");
		if (load != 0)
		{
			loads.Add(load);
		}
	}

// when more then one type of element is implemented, perform a type check 
	mystr elem_type = edc_plateparameters.TreeGetString("element_type");
	mystr node_type = edc_plateparameters.TreeGetString("node_type");

// get current lists for mesh object
	double* vdata;
	int vlen;
	edc_handle_mesh->TreeGetVector("list_of_nodes",&vdata,vlen);            
	TArray<double> mesh_nnrs; mesh_nnrs.CopyFrom(vdata,vlen); mesh_nnrs.SetLen(vlen);			// node numbers already associated with the Mesh
	edc_handle_mesh->TreeGetVector("list_of_elements",&vdata,vlen);
	TArray<double> mesh_enrs; mesh_enrs.CopyFrom(vdata,vlen); mesh_enrs.SetLen(vlen);			// element numbers already associated with the Mesh

// *******************

// GENERATION OF NODES

// use existing script command AddNode
	mystr commandname("AddNode");
	int nparam, hasretval, c_option;
	CEDC_Command* command_addnode = edcParser->GetCommand( edcParser->GetCommandType(commandname, nparam, hasretval, c_option) );


// prepare EDC to be able to use existing script commands 
	ElementDataContainer edc_nodeproperties;
	TArray<ElementDataContainer*> paramEDCforNode;
	paramEDCforNode.Add(&edc_nodeproperties);
	ed.SetText(node_type,"node_type");																		edc_nodeproperties.Add(ed);

//  todo: automatically compute orientation for all nodes - depending on node_type ...
	Vector3D p12 = p2-p1;						//  e1 = p12/|p12|
	Vector3D p13 = p3-p1;						//  e2 = p13/|p13|
	Vector3D n(p12.Cross(p13));			// normal to the plane od the plate

	Vector3D e1(p12); e1.Normalize();
	ed.SetVector3D(e1.X(), e1.Y(), e1.Z(), "reference_slope1");	edc_nodeproperties.TreeAdd("Geometry",ed);
	Vector3D e2(p13); e2.Normalize();
	ed.SetVector3D(e2.X(), e2.Y(), e2.Z(), "reference_slope2");	edc_nodeproperties.TreeAdd("Geometry",ed);

// override discretization if entry element_size is present
	int j = edc_plateparameters.Find("element_size");										
	if (j != 0 && edc_plateparameters.Get(j).IsVector())
	{
		Vector2D elemsize;
		edc_plateparameters.TreeGetVector2D("element_size",elemsize.X(), elemsize.Y());
		double total_length, required;
		total_length = p12.Norm();
		required = total_length / elemsize.X();
		discr.X() = ceil(required);
		total_length = p13.Norm();
		required = total_length / elemsize.Y();
		discr.Y() = ceil(required);
	}


	TArray<double> nnrs;                                                  // node numbers for the current call GenerateBeam
	for(int i=0; i<=(int)discr.Y(); i++)
	{
		for(int j=0; j<=(int)discr.X(); j++)
		{
			Vector3D pos = p1 + p12*((double)j/discr.X()) + p13*((double)i/discr.Y());
			ed.SetVector3D(pos.X(),pos.Y(),pos.Z(),"reference_position");				edc_nodeproperties.TreeAdd("Geometry",ed);

// call AddNode with prepared EDC, returnvalue is number of added node
			command_addnode->ExecuteCommand(mbs,edcParser,paramEDCforNode,ed,0);

			nnrs.Add(ed.GetInt());
		}
	}
// return value:nodes from current generation function
	ed.SetVector( nnrs.GetDataPtr(), nnrs.Length(), "list_of_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Plate or empty");
	edc_rv.Add(ed);
// update list_of_nodes in Mesh
	mesh_nnrs.Merge(nnrs);
	ed.SetVector( mesh_nnrs.GetDataPtr(), mesh_nnrs.Length(), "list_of_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Mesh or empty");
	ElementDataContainer update_nnrs; update_nnrs.Add(ed);
	edc_handle_mesh->TreeReplaceEDCDataWith(&update_nnrs);


// ***************************
// GENERATION OF PLATE ELEMENTS

	double size1 = p12.Norm() / discr.X();
	double size2 = p13.Norm() / discr.Y();

// use existing script command AddElement
	commandname = mystr("AddElement");
	CEDC_Command* command_addelement = edcParser->GetCommand( edcParser->GetCommandType(commandname, nparam, hasretval, c_option) );
	
	// prepare EDC to be able to use existing script commands 
	ElementDataContainer edc_elementproperties;
	ed.SetText(elem_type,"element_type");	edc_elementproperties.Add(ed);
	ed.SetInt(matnr,"material_number");		edc_elementproperties.TreeAdd("Physics",ed);
	ed.SetDouble(thickness,"thickness");	edc_elementproperties.TreeAdd("Geometry",ed);
	ed.SetDouble(size1,"size1");	edc_elementproperties.TreeAdd("Geometry",ed);
	ed.SetDouble(size2,"size2");	edc_elementproperties.TreeAdd("Geometry",ed);
	ed.SetBool(is_nonlinear, "is_geometrically_nonlinear");		edc_elementproperties.TreeAdd("Physics",ed);
	
	if (loads.Length() > 0)
	{
		ed.SetVector(loads.GetDataPtr(), loads.Length(), "loads");	edc_elementproperties.Add(ed);
	}

	TArray<ElementDataContainer*> paramEDCforElement;
	paramEDCforElement.Add(&edc_elementproperties);


	TArray<double> enrs;
	for(int i=1; i<=discr.Y(); i++)
	{
		for(int j=1; j<=discr.X(); j++)
		{
			ed.SetInt(nnrs(j   +(i-1)*(discr.X()+1)),"node_1");		edc_elementproperties.TreeAdd("Geometry",ed);
			ed.SetInt(nnrs(j+1 +(i-1)*(discr.X()+1)),"node_2");		edc_elementproperties.TreeAdd("Geometry",ed);
			ed.SetInt(nnrs(j+1 +(i  )*(discr.X()+1)),"node_3");		edc_elementproperties.TreeAdd("Geometry",ed);
			ed.SetInt(nnrs(j   +(i  )*(discr.X()+1)),"node_4");		edc_elementproperties.TreeAdd("Geometry",ed);

	// call AddElement with prepared EDC, returnvalue is number of added element
			command_addelement->ExecuteCommand(mbs,edcParser,paramEDCforElement,ed,0);

			enrs.Add(ed.GetInt());
		}
	}
// return value:elements from current generation function
	ed.SetVector( enrs.GetDataPtr(), enrs.Length(), "list_of_elements");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Elements associated with the Beam or empty");
	edc_rv.Add(ed);
// update list_of_elements in Mesh
	mesh_enrs.Merge(enrs);
	ed.SetVector( mesh_enrs.GetDataPtr(), mesh_enrs.Length(), "list_of_elements");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Elements associated with the Mesh or empty");
	ElementDataContainer update_enrs; update_enrs.Add(ed);
	edc_handle_mesh->TreeReplaceEDCDataWith(&update_enrs);

	// first approach, slow

	return rv;
}


mystr CEDC_ComMesh_GeneratePlate::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter EDC containing the plate properties \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("enties in the properties EDC are: \n");
	descr +=			mystr("\\begin{itemize} \n");	
	descr +=			mystr("	 \\item P1 - position of first node \n");
	descr +=			mystr("	 \\item P2 - position of outer node along first direction\n");
	descr +=			mystr("	 \\item P3 - position of outer node along second direction \n");
	descr +=      mystr("  \\item matnr - number of the material \n"); 
	descr +=      mystr("  \\item discretization - number of the plate elements both\n");
	descr +=      mystr("  \\item elementsize - (optional) maximum element size ( if set has priority over discretization ) \n");
	descr +=			mystr("\\end{itemize} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};


int CEDC_ComMesh_GetNodesInBox::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
	ElementDataContainer edc_rv;

// container 1 contains the corners of the box
	Vector3D p1, p2;
	ElementDataContainer* edc_boxparameters = parameter_EDCs(1);
	edc_boxparameters->TreeGetVector3D("P1",p1.X(), p1.Y(), p1.Z());
	edc_boxparameters->TreeGetVector3D("P2",p2.X(), p2.Y(), p2.Z());
	Box3D box(p1,p2);

// container 2 contains the name of the handle ( can not pass the EDC* itself because that is going to be deleted... )
// consistency already checked in ParseCommandParameters...
	mystr handlename = parameter_EDCs(2)->TreeGetString("handle");
	ElementDataContainer* edc_local = edcParser->GetMBSParser()->GetLocalEDC(); 
	int nr_handle = edc_local->Find(handlename);
	ElementDataContainer* edc_handle_mesh = edc_local->Get(nr_handle).GetEDC();


// get current lists for mesh object
	double* vdata;
	int vlen;
	edc_handle_mesh->TreeGetVector("list_of_nodes",&vdata,vlen);            
	TArray<double> mesh_nnrs; mesh_nnrs.CopyFrom(vdata,vlen); mesh_nnrs.SetLen(vlen);			// node numbers already associated with the Mesh
	
	//edc_handle_mesh->TreeGetVector("list_of_redundant_nodes",&vdata,vlen);
	//TArray<double> mesh_redundant_nodes; mesh_redundant_nodes.CopyFrom(vdata,vlen); mesh_redundant_nodes.SetLen(vlen);	// redundant nodes in mesh

	TArray<double> nnrs;
// loop over all nodes
	for (int i=1; i<= mesh_nnrs.Length(); i++)
	{
		Vector3D nodepos = mbs->GetNode(mesh_nnrs(i)).RefConfPos();
		if (box.IsIn(nodepos) /*&& !mesh_redundant_nodes.Find(mesh_nnrs(i))*/)
		{
			nnrs.Add(mesh_nnrs(i));
		}
	}

// return value:nodes from current generation function
	return_value.SetVector( nnrs.GetDataPtr(), nnrs.Length(), "list_of_nodes");		return_value.SetValuesInt();	return_value.SetVariableLength();	return_value.SetToolTipText("List of Nodes returned from GetNodesInBox");

	return rv;
}

mystr CEDC_ComMesh_GetNodesInBox::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter EDC containing the box (defined by two corners) \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("enties in the properties EDC are: \n");
	descr +=			mystr("\\begin{itemize} \n");	
	descr +=			mystr("	 \\item P1 - position of first corner \n");
	descr +=			mystr("	 \\item P2 - position of second corner\n");
	descr +=			mystr("\\end{itemize} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

int CEDC_ComMesh_GlueMesh::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	// default parameter values, to be overridden by user
	ElementDataContainer edc_glue_parameters;
	ElementData ed;

	Vector nodes;
	ed.SetVector(nodes.GetVecPtr(), 0, "nodes_to_constrain");	ed.SetValuesInt(); ed.SetVariableLength();			edc_glue_parameters.Add(ed);
	ed.SetBool(0, "is_ground");	ed.SetToolTipText("default 0: constrain node pairs - set to 1 for ground constraint");	edc_glue_parameters.Add(ed);
	ed.SetBool(0, "remove_redundant_only");		ed.SetToolTipText("default 0: connect pairs - set to 1 to only remove redundant nodes"); edc_glue_parameters.Add(ed);
	ed.SetVector3D(1., 1., 1., "constrained_directions");	ed.SetValuesInt();										edc_glue_parameters.Add(ed);
	ed.SetVector3D(1., 1., 1., "constrained_rotations");	ed.SetValuesInt();										edc_glue_parameters.Add(ed);

	edc_glue_parameters.TreeReplaceEDCDataWith(parameter_EDCs(1));

	int len; double* p_vec;
	edc_glue_parameters.TreeGetVector("nodes_to_constrain", &p_vec, len);	// ist das so gut????
	nodes.LinkWith(p_vec, len);

	int is_ground = edc_glue_parameters.TreeGetBool("is_ground", 0);
	int redundant_only = edc_glue_parameters.TreeGetBool("remove_redundant_only", 0);


	Vector constrained_dirs, constrained_rots;
	edc_glue_parameters.TreeGetVector("constrained_directions", &p_vec, len);
	constrained_dirs.LinkWith(p_vec, len);
	edc_glue_parameters.TreeGetVector("constrained_rotations", &p_vec, len);
	constrained_rots.LinkWith(p_vec, len);

	// HACK: first of all, build a node2elements list
	mystr handlename = parameter_EDCs(2)->TreeGetString("handle");
	ElementDataContainer* edc_local = edcParser->GetMBSParser()->GetLocalEDC(); 
	int nr_handle = edc_local->Find(handlename);
	ElementDataContainer* edc_handle_mesh = edc_local->Get(nr_handle).GetEDC();

	double* elems_data, * nodes_data, * redundant_nodes_data; 
	int elems_len, nodes_len, redundant_nodes_len;
	
	edc_handle_mesh->TreeGetVector("list_of_nodes",&nodes_data,nodes_len);            
	TArray<double> mesh_nodes; mesh_nodes.CopyFrom(nodes_data,nodes_len); mesh_nodes.SetLen(nodes_len);			// node numbers already associated with the Mesh
	
	edc_handle_mesh->TreeGetVector("list_of_elements",&elems_data,elems_len);
	TArray<double> mesh_elems; mesh_elems.CopyFrom(elems_data,elems_len); mesh_elems.SetLen(elems_len);			// element numbers already associated with the Mesh
	
	edc_handle_mesh->TreeGetVector("list_of_redundant_nodes",&redundant_nodes_data,redundant_nodes_len);            
	TArray<double> mesh_redundant_nodes; mesh_redundant_nodes.CopyFrom(redundant_nodes_data,redundant_nodes_len); mesh_redundant_nodes.SetLen(redundant_nodes_len);			// node numbers already associated with the Mesh

	TArray<TArray<int>> nodes2elems;
	nodes2elems.SetLen(mesh_nodes.Length());
	for (int i = 1; i <= mesh_elems.Length(); i++)
	{
		Element& e = mbs->GetElement(mesh_elems(i));
		for (int j = 1; j <= e.NNodes(); j++)
		{
			nodes2elems(e.GetNode(j).NodeNum()).Add(mesh_elems(i));
		}
	}

	ElementDataContainer edc_elementproperties;		
	ed.SetText("GenericBodyJoint","element_type");		edc_elementproperties.Add(ed);

	int ed_num = edc_glue_parameters.Find("constrained_directions");
	if (ed_num)
	{
		edc_elementproperties.TreeAdd("Physics.Lagrange", edc_glue_parameters.Get(ed_num));
	}
	ed_num = edc_glue_parameters.Find("constrained_rotations");
	if (ed_num)
	{
		edc_elementproperties.TreeAdd("Physics.Lagrange", edc_glue_parameters.Get(ed_num));
	}
	/*ed.SetText("UniversalJoint", "element_type");	edc_elementproperties.Add(ed);
	ed.SetVector3D(0., 0., 1., "axis");				edc_elementproperties.TreeAdd("Position1", ed);
	ed.SetVector3D(0., 0., 1., "axis");				edc_elementproperties.TreeAdd("Position2", ed);*/

	mystr commandname = mystr("AddConnector");
	int nparam, hasretval, c_option;
	CEDC_Command* command_addrigidjoint = edcParser->GetCommand( edcParser->GetCommandType(commandname, nparam, hasretval, c_option) );

	TArray<double> new_redundant_nodes;		// nodes to remove from mbs
	// pairs
	if (!is_ground)
	{
		TArray<int2> node_pairs;
		
		for (int i = 1; i <= nodes.Length(); i++)
		{
			for (int j = i+1; j <= nodes.Length(); j++)
			{
				if ( nodes(j)!=0 && nodes(i)!=0 ) // skip redundant nodes ( marked in a previous iteration )
				{
					Node& n1 = mbs->GetNode((int)nodes(i));
					Node& n2 = mbs->GetNode((int)nodes(j));

					Vector3D dist = n1.RefConfPos() - n2.RefConfPos();

					if (dist.Norm() < 1e-10)
					{
						// check if nodes are identical; if so, replace n2 by n1
						// TODO: == operator for nodes!
						if (typeid(n1) == typeid(ANCFNodeS1S2_3D) && typeid(n2) == typeid(ANCFNodeS1S2_3D) &&
							(dynamic_cast<ANCFNodeS1S2_3D&>(n1)).GetRefSlope1() == (dynamic_cast<ANCFNodeS1S2_3D&>(n2)).GetRefSlope1() &&
							(dynamic_cast<ANCFNodeS1S2_3D&>(n1)).GetRefSlope2() == (dynamic_cast<ANCFNodeS1S2_3D&>(n2)).GetRefSlope2())
						{
							// replace in elements
							for (int k = 1; k <= nodes2elems(nodes(j)).Length(); k++)
							{
								Element& e = mbs->GetElement(nodes2elems(nodes(j))(k));
								for (int l = 1; l <= e.NNodes(); l++)
								{
									if (e.NodeNum(l) == n2.NodeNum())
									{
										e.NodeNum(l) = n1.NodeNum();
									}
								}
							}
							// remove redundant entries from pair list
							for (int m = 1; m <= node_pairs.Length(); m++)
							{
								if( node_pairs(m)(1) == nodes(j) || node_pairs(m)(2) == nodes(j))  
								{
									node_pairs.Erase(m);
								}
							}
							new_redundant_nodes.Add(nodes(j));
							nodes(j) = 0;             // ignore this nodelist entry - has beed identified as redundant
						}
						// if not, apply RigidBodyJoint later on
						else
						{
							if (!redundant_only)
							{
								int2 pair(nodes(i),nodes(j));
								node_pairs.Add(pair);
							}
						}
					}
				}
			}
		}

		for (int i = 1; i <= node_pairs.Length(); i++)
		{
			Node& n1 = mbs->GetNode(node_pairs(i)(1));
			Node& n2 = mbs->GetNode(node_pairs(i)(2));

			Element& e1 = mbs->GetElement(nodes2elems(n1.NodeNum())(1));
			Element& e2 = mbs->GetElement(nodes2elems(n2.NodeNum())(1));

			//mbs->UO() << n1.NodeNum() << " "  << n2.NodeNum() << "\n";
			//mbs->UO() << nodes2elems(n1.NodeNum())(1) << " "  << nodes2elems(n2.NodeNum())(1) << "\n";
			
			// find the local position of a node within an element
			Vector3D p_loc_n1, p_loc_n2;
			for (int j = 1; j <= e1.NNodes(); j++)
			{
				if (e1.NodeNum(j) == n1.NodeNum())
				{
					p_loc_n1 = e1.GetNodeLocPos(j);
				}
				if (e2.NodeNum(j) == n2.NodeNum())
				{
					p_loc_n2 = e2.GetNodeLocPos(j);
				}
			}
			//mbs->UO() << p_loc_n1 << " "  << p_loc_n2 << "\n";
			
			// prepare EDC to be able to use existing script commands 
			ed.SetInt(nodes2elems(n1.NodeNum())(1), "element_number");	edc_elementproperties.TreeAdd("Position1", ed);
			ed.SetInt(nodes2elems(n2.NodeNum())(1), "element_number");	edc_elementproperties.TreeAdd("Position2", ed);
			ed.SetVector3D(p_loc_n1.X(), p_loc_n1.Y(), p_loc_n1.Z(), "position");	edc_elementproperties.TreeAdd("Position1", ed);
			ed.SetVector3D(p_loc_n2.X(), p_loc_n2.Y(), p_loc_n2.Z(), "position");	edc_elementproperties.TreeAdd("Position2", ed);

			TArray<ElementDataContainer*> paramEDCforElement;
			paramEDCforElement.Add(&edc_elementproperties);

			command_addrigidjoint->ExecuteCommand(mbs,edcParser,paramEDCforElement,ed,0);
		}
	}
	else		// is_ground
	{
		// again, find node pairs in order to avoid redundant constraints
		// assume, no redundant nodes exist eny longer 
		TArray<int2> node_pairs;
		for (int i = 1; i <= nodes.Length(); i++)
		{
			for (int j = i+1; j <= nodes.Length(); j++)
			{
				Node& n1 = mbs->GetNode((int)nodes(i));
				Node& n2 = mbs->GetNode((int)nodes(j));

				Vector3D dist = n1.RefConfPos() - n2.RefConfPos();

				if (dist.Norm() < 1e-10)
				{
					int2 pair(nodes(i),nodes(j));
					node_pairs.Add(pair);
				}
			}
		}

		for (int i = 1; i <= nodes.Length(); i++)
		{
			Node& n = mbs->GetNode(nodes(i));
			Vector3D refpos = n.RefConfPos();

			Element& e = mbs->GetElement(nodes2elems(n.NodeNum())(1));
			int is_in_pair = 0;
			for (int j = 1; j <= node_pairs.Length(); j++)
			{
				if (node_pairs(j)(2) == nodes(i))
				{
					is_in_pair = 1; break;
				}
			}

			if (!is_in_pair)
			{
				// find the local position of a node within an element
				Vector3D p_loc;
				for (int j = 1; j <= e.NNodes(); j++)
				{
					if (e.NodeNum(j) == n.NodeNum())
					{
						p_loc = e.GetNodeLocPos(j); break;
					}
				}
				
				// prepare EDC to be able to use existing script commands 
				ed.SetInt(nodes2elems(n.NodeNum())(1), "element_number");		edc_elementproperties.TreeAdd("Position1", ed);
				ed.SetInt(0, "element_number");									edc_elementproperties.TreeAdd("Position2", ed);
				ed.SetVector3D(p_loc.X(), p_loc.Y(), p_loc.Z(), "position");	edc_elementproperties.TreeAdd("Position1", ed);
				ed.SetVector3D(refpos.X(), refpos.Y(), refpos.Z(), "position");	edc_elementproperties.TreeAdd("Position2", ed);

				TArray<ElementDataContainer*> paramEDCforElement;
				paramEDCforElement.Add(&edc_elementproperties);

				command_addrigidjoint->ExecuteCommand(mbs,edcParser,paramEDCforElement,ed,0);
			}
		}
	}

	TArray<int> permutation; permutation.SetLen(new_redundant_nodes.Length());
	QuicksortDouble(new_redundant_nodes, permutation);

	for (int i = new_redundant_nodes.Length(); i >= 1; i--)
	{
		int node_to_delete = new_redundant_nodes(i);
		mbs->DeleteNode(node_to_delete);
		mesh_nodes.Erase(mesh_nodes.Find(new_redundant_nodes(i)));
		
		for (int j = 1; j <= mesh_nodes.Length(); j++)
		{
			if (mesh_nodes(j) > node_to_delete)
			{
				mesh_nodes(j)--;
			}
		}
	}

	ElementDataContainer update_nodes;
	mesh_redundant_nodes.Merge(new_redundant_nodes);
	ed.SetVector( mesh_redundant_nodes.GetDataPtr(), mesh_redundant_nodes.Length(), "list_of_redundant_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of redundant nodes in mesh");
	update_nodes.Add(ed);
	ed.SetVector( mesh_nodes.GetDataPtr(), mesh_nodes.Length(), "list_of_nodes");	ed.SetValuesInt();	ed.SetVariableLength();	ed.SetToolTipText("List of Nodes associated with the Plate or empty");
	update_nodes.Add(ed);
	edc_handle_mesh->TreeReplaceEDCDataWith(&update_nodes);

	return 1;
}


mystr CEDC_Check_Entry_Exists::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter: string with the name and tree of the entry to check \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

int CEDC_Check_Entry_Exists::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	ElementData* ED_treename = parameter_EDCs(1)->GetPtr(1);
  mystr str_treename;
	if(ED_treename->IsText())
	{
		str_treename = ED_treename->GetText();
	}

	ElementDataContainer* localEDC = edcParser->GetMBSParser()->GetLocalEDC2();  // existing local EDC for the CHECK_ENTRY_EXISTS command

	ElementData* ed = localEDC->TreeFind(str_treename);
	if(ed==NULL)
	{
	// Entry not found
		return_value.SetInt(0,"rv");
	}
	else if(ed!=NULL && !ed->IsEDC())
	{
	// Entry found, not an EDC
		return_value.SetInt(1,"rv");
	}
	else if(ed!=NULL && ed->IsEDC())
	{
	// Entry found, also an EDC
		return_value.SetInt(2,"rv");
	}
	return 1;
}


mystr CEDC_Compare::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter: string A \n");
	descr +=			mystr("	 \\item $2^{nd}$ parameter: string B \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

int CEDC_Compare::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
	ElementData* ED_strA = parameter_EDCs(1)->GetPtr(1);
  mystr strA;
	if(ED_strA->IsText())
	{
		strA = ED_strA->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing COMPARE operation: first parameter is not a string variable");
		return rv;
	}
	ElementData* ED_strB = parameter_EDCs(2)->GetPtr(1);
  mystr strB;
	if(ED_strB->IsText())
	{
		strB = ED_strB->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing COMPARE operation: second parameter is not a string variable");
		return rv;
	}

	int result = strcmp(strA.c_str(), strB.c_str());
	
	return_value.SetInt(result,"rv");

	return 1;
}


mystr CEDC_StrCat::GetCommandTexParameterDescription() const 
{
	mystr descr = mystr("The parameters of this command are \n");
	descr +=			mystr("\\begin{enumerate} \n");
	descr +=			mystr("	 \\item $1^{st}$ parameter: string A \n");
	descr +=			mystr("	 \\item $2^{nd}$ parameter: string B \n");
	descr +=			mystr("\\end{enumerate} \n");
	descr +=			mystr("-");	//DR hack, to avoid bugs in latex
	return descr;
};

int CEDC_StrCat::ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int rv = 0;
	ElementData* ED_strA = parameter_EDCs(1)->GetPtr(1);
  mystr strA;
	if(ED_strA->IsText())
	{
		strA = ED_strA->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing COMPARE operation: first parameter is not a string variable");
		return rv;
	}
	ElementData* ED_strB = parameter_EDCs(2)->GetPtr(1);
  mystr strB;
	if(ED_strB->IsText())
	{
		strB = ED_strB->GetText();
	}
	else
	{
		edcParser->EDCError("Error while executing COMPARE operation: second parameter is not a string variable");
		return rv;
	}

	mystr result = strA + strB;

	return_value.SetText(result,"rv");

	return 1;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//void CEDCParser::EDC2TexTable(const ElementDataContainer& edc, mystr& str, const mystr& indent, const mystr& treename) //$ DR 2013-01-11 moved to models_auto_documentation

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int CEDCParser::GetParameterEDC(mystr text, TArray<ElementDataContainer*>& parameter_EDCs)
{
	int parPos = 0;
	int parCount = 1;
	mystr EDCtext;

	text.GetUntil(parPos, ',', EDCtext, 0);
	while(1)		// handle all parameters
	{
		bool success = false;
		//int future: this should be either an EDC name, or a structure defining an EDC {...}
		//at the moment, only an EDC-name is possible
		EDCtext.EraseSpacesHeadTail();

		//the EDC must be either global or local
		ElementDataContainer* edc = mbsParser->GetLocalEDC(); //local
		ElementDataContainer* edc2 = mbsParser->GetLocalEDC2(); //global
		ElementDataContainer* edcMBS = mbs->GetModelDataContainer();	//$ DR 2013-07-25 added for the case of multiple input files

		//$ RL 2012-7-19:[ new code with ExecuteCommand planned to be possible only with arguments
		ElementData* ed = 0; 
		if(edc)
		{
			ed = edc->TreeFind(EDCtext); // search only if edc exists!!
		}
		ElementData* ed2 = 0;
		if(edc2)
		{
			ed2 = edc2->TreeFind(EDCtext); // search only if edc exists!!
		}
		ElementData* edMBS = 0;
		if(edcMBS)
		{
			edMBS = edcMBS->TreeFind(EDCtext); //$ DR 2013-07-25 added for the case of multiple input files
		}

		if (!ed && ed2) {ed = ed2;}
		else if (!ed && !ed2 && edMBS) {ed = edMBS;}	//$ DR 2013-07-25 added for the case of multiple input files
		//else if (!ed && !ed2) 
		else if (!ed && !ed2 && !edMBS)
		{
			// no parameter edc found -> add parameter string to parameter_EDCs
			int pos=0;
			//if (text.PosPeek(pos) == doublequote)
			if (EDCtext.PosPeek(pos) == doublequote)
			{
				//pos = text.Length()-1;
				pos = EDCtext.Length()-1;
				//if(text.PosPeek(pos) == doublequote)
				if(EDCtext.PosPeek(pos) == doublequote)
				{
					// text is within doublequotes
					pos=0;
					//text = text.GetStringInBrackets(pos, doublequote, doublequote); // remove doublequotes
					EDCtext = EDCtext.GetStringInBrackets(pos, doublequote, doublequote); // remove doublequotes				
				}
			}
			ElementData ed3;
			mystr parName = "parameter"+mystr(parCount);
			//ed3.SetText(text,parName); // store parameter text
			ed3.SetText(EDCtext,parName); // store parameter text
			ElementDataContainer edc_text;
			edc_text.Add(ed3);
			parameter_EDCs.Add(edc_text.GetCopy()); //parameter_EDCs are deleted outside this function 
			//return 1;
			success = true;
		}

		if(!success)
		{
			//$ DR 2012-11-29: [ also use element data, e.g. for CEDC_ComPrint
			if (!ed->IsEDC() && (ed->IsBool() || ed->IsInt() || ed->IsDouble() || ed->IsMatrix() || ed->IsText() || ed->IsVector()))
			{
				ElementDataContainer edc_ed;
				edc_ed.Add(*ed->GetCopy());
				parameter_EDCs.Add(edc_ed.GetCopy()); //parameter_EDCs are deleted outside this function 
				//return 1;
				success = true;
			}
			//$ DR 2012-11-29: ]
		}

		if(!success)
		{
			//always use ed, ed2 overrides ed:
			//now the parameter should indicate an EDC
			if (ed->IsEDC())
			{
				parameter_EDCs.Add((ed->GetEDC())->GetCopy());  //parameter_EDCs are deleted outside this function 
			}
			else
			{
				//failed
				EDCError(mystr("ERROR: Expected a parameter set (ElementDataContainer) for function parameters"));
				return 0;
			}
		}

		if(parPos == -1) return 1;	// last parameter
		parPos++;
		text.GetUntil(parPos, ',', EDCtext, 0);
		parCount++;
	}

	return 1;
}	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//this function is already adapted for TreeEDCs!
//$ DR 2013-01-31 empty subEDCs are not written at all, see log 384
void CEDCParser::EDC2String(const ElementDataContainer& edc, mystr& str, const mystr& indent)
{
//symbols already defined in CEDCParser:
//	mystr el = "\n"; //end of line
//	mystr ob = "{";  //open brace
//	mystr cb = "}";  //closing brace
//	mystr osb = "[";  //open square bracket
//	mystr csb = "]";  //closing square bracket
//	mystr ds = ": "; //double stop
//	mystr semicolon = ";"; //semicolon
//	char commentc = '%'; //comment
//	mystr doublequote = mystr('"'); //for strings or text

//symbols with additional whitespaces:
	mystr eq = "= "; //equal sign
	mystr comma = ", "; //comma
	mystr comment = mystr("  ")+mystr(commentc); //comment

	mystr matrixsep = el; //may be changed to semicolon ...
	str = "";
	mystr num;

	for (int i=1; i <= edc.Length(); i++)
	{
		const ElementData& ed = edc.Get(i);

		if (!ed.IsLocked() && !ed.IsElementInfo() && !ed.IsCompAction() && !ed.IsWCDriverAction())
		{
			//str += indent + ed.GetDataName();
			if (ed.IsEDC())
			{ //hierarchically
				mystr substr;
				const ElementDataContainer* edcptr = ed.GetEDC();
				if (edcptr)
				{
					EDC2String(*edcptr, substr, indent+mystr("  "));
				}
				if(substr.Length())
				{
					str += indent + ed.GetDataName();
					str += el + indent + ob + el;
					str += substr;
					str += indent + cb + el;
				}
			}
			else
			{ //just this list
				str += indent + ed.GetDataName();//
				str += eq;
				if (ed.IsBool())
				{
					str += mystr(ed.GetBool());
				}
				else if (ed.IsInt())
				{
					str += mystr(ed.GetInt());
				}
				else if (ed.IsDouble())
				{
					num.SmartDouble2String(ed.GetDouble());
					str += num;
				}
				else if (ed.IsText())
				{
					str += mystr(doublequote) + mystr(ed.GetText()) + mystr(doublequote);
				}
				else if (ed.IsVector())
				{
					str += osb;
					for (int j = 1; j <= ed.GetVectorLen(); j++)
					{
						num.SmartDouble2String(ed.GetVectorVal(j));
						str += num;
						if (j < ed.GetVectorLen()) str += comma;
					}
					str += csb;
				}
				else if (ed.IsMatrix())
				{
					str += osb;
					for (int j1 = 1; j1 <= ed.GetMatrixRows(); j1++)
					{
						for (int j2 = 1; j2 <= ed.GetMatrixCols(); j2++)
						{
							num.SmartDouble2String(ed.GetMatrixVal(j1, j2));
							str += num;
							if (j2 < ed.GetMatrixCols()) str += comma;
						}
						if (j1 < ed.GetMatrixRows()) str += matrixsep + indent + mystr("    ");
					}
					str += csb;
				}
				if (ed.HasToolTip())
				{
					str += comment + ed.GetToolTipText();
				}
				str += el;
			}
		}
	}
}

//old version for saving MBS data, do not use!
//$JG2012-01-27: mbs is available in CEDCParser
int CEDCParser::String2EDC(/*MBS* mbs, */mystr& str, ElementDataContainer& edc, ElementDataContainer& testedc)
{
	char osb = '[';  //open square bracket
	char csb = ']';  //closing square bracket
	char ds = ':'; //double stop
	char eq = '='; //equal sign
	char comma = ','; //comma sign
	char semicolon = ';'; //semicolon sign

	int oo = 0;
	if (oo) mbs->UO() << "String2EDC\n";

	mystr dataname, data, data2, data3;
	int endstr = 0;
	int pos = 0;
	int posold;
	TArray<double> datalist;


	while (!endstr)
	{
		posold = pos;
		dataname = str.GetWord(pos, 0);
		if (pos != -1 && str[pos] == ds)
		{
			str = str.Right(str.Length()-posold); return 1;
		}
		if (pos != -1)
		{
			if (oo) mbs->UO() << "Dataname='" << dataname << "'\n";

			if (!str.GoUntil(pos, eq)) {return 0;}; //write error

			int x = testedc.Find(dataname.c_str());
			if (x)
			{
				ElementData& ed = testedc.Get(x);
				if (ed.IsVector())
				{
					//read [ ... ]
					if (!str.GoUntil(pos, osb)) {LoadError("Vector '[...]' expected!"); return 0;}

					if (!str.GetUntil(pos, csb, data, 0)) {LoadError(mystr("Vector '[...]' expected, but received '")+data+mystr("'")); return 0;};

					int len = ed.GetVectorLen();
					datalist.SetLen(0);

					int pos2 = 0;
					while (pos2 != -1 && (datalist.Length() <= len || ed.IsVariableLength()))
					{
						data.GetUntil(pos2, comma, data2, 0);

						datalist.Add(data2.MakeDouble());
						if (pos2 != -1) data.GoUntil(pos2, comma);
					}

					if (datalist.Length() != len && !ed.IsVariableLength()) {LoadError(mystr(ed.GetDataName())+mystr(": Vector of length ")+mystr(len)+mystr(" expected, but received length=")+mystr(datalist.Length())); return 0;};

					len = datalist.Length();
					Vector v(len);
					for (int i=1; i <= len; i++) v(i) = datalist(i);

					ed.SetVector(v.GetVecPtr(), len, dataname.c_str());
					edc.Add(ed);
					if (oo) mbs->UO() << v << "\n";

					str.GoUntil(pos, csb); //read final square bracket


				}
				else if (ed.IsMatrix())
				{
					//read [ ... ]
					if (!str.GoUntil(pos, osb)) {LoadError("Matrix '[...]' expected!"); return 0;}

					if (!str.GetUntil(pos, csb, data, 0)) {LoadError(mystr("Matrix '[...]' expected, but received '")+data+mystr("'")); return 0;};

					int rows = ed.GetMatrixRows();
					int cols = ed.GetMatrixCols();


					TArray<double> doublelist;

					//++++++++++++++++++++++++
					//same as in CustomEditDialog!
					char cstr[258];
					int limit = 256;
					int startpos = 0;

					int end = 0;
					int endpos, endpos2;
					int maxcol = 0;
					int error = 0;
					int nrows = 0;

					data = data+mystr("\n\n");

					while (!end)
					{
						int emptystring = 1;
						int ncol = 0;
						int endcol = 0;

						while (!endcol)
						{
							endpos = data.Find(startpos,',');
							endpos2 = data.FindEOL(startpos);

							if (endpos2 < endpos || endpos == -1) 
							{
								endpos = endpos2;
								endcol = 1;
							}

							if (endpos != -1 && endpos != 0 && startpos < endpos && endpos+1 <= data.Length()-1 && data.Length() != 0)
							{
								int len = endpos-startpos;
								data.CopySubStringNoSpaces(cstr, startpos, endpos-1, limit);

								if (strlen(cstr) != 0)
								{
									doublelist.Add(atof(cstr));
									ncol++;

									startpos = endpos+1;
									emptystring = 0;
								}
								else
								{
									startpos = endpos+1;
									doublelist.Add(0);
								}
							}
							else 
							{
								end = 1;
								endcol = 1;
							}
						}
						if (!end && !emptystring) nrows++;
						if (ncol != maxcol) 
						{
							if (maxcol != 0 && ncol != 0)
							{
								error = 1;
							}

							if (ncol > maxcol) maxcol = ncol;
						}
					}



					//++++++++++++++++++++++++


					/*
					int pos2 = 0;
					int foundcols = 0;
					TArray<double> doublelist;

					while (pos2 != -1)
					{
					data.GetUntil(pos2, semicolon, data3, 0);

					int pos3 = 0;
					while (pos3 != -1)
					{
					data3.GetUntil(pos3, comma, data2, 0);
					doublelist.Add(data2.MakeDouble());
					if (pos3 != -1) data3.GoUntil(pos3, comma);
					}
					if (!foundcols) foundcols = doublelist.Length(); //take first number of columns for matrix

					data.GoUntil(pos2, semicolon);
					}
					*/

					str.GoUntil(pos, csb); //read final square bracket

					if (cols == 0 || ed.IsVariableLength()) cols = maxcol;
					else if (cols != maxcol || error) {LoadError(mystr("Matrix should contain ")+mystr(rows)+mystr("rows and ")+mystr(cols)+mystr("columns !!!")); return 0;}

					Matrix m;
					if (cols != 0) 
					{
						if (rows == 0 || ed.IsVariableLength()) rows = doublelist.Length()/cols;
						else if (doublelist.Length()/cols != rows) {LoadError(mystr("Matrix should contain ")+mystr(rows)+mystr("rows and ")+mystr(cols)+mystr("columns !!!")); return 0;}

						m.SetSize(rows, cols);
						for (int i1=1; i1 <= rows; i1++)
						{
							for (int i2=1; i2 <= cols; i2++)
							{
								int ind = (i1-1)*cols+i2;
								if (ind <= doublelist.Length())
								{
									m(i1, i2) = doublelist(ind);
								}
							}
						}
					}
					else rows = 0;

					ed.SetMatrix(m.GetMatPtr(), rows, cols, dataname.c_str());
					edc.Add(ed);

					if (oo) mbs->UO() << m << "\n";

				}
				else
				{
					//data = str.GetWord(pos, 0);
					str.GetUntil(pos, '\n', data);
					data.EraseSpaces();
					if (oo) mbs->UO() << "single data='" << data << "'\n";
					if (ed.IsBool())
					{
						ed.SetBool(data.MakeInt()); edc.Add(ed);
					} 
					else if (ed.IsInt())
					{
						ed.SetInt(data.MakeInt()); edc.Add(ed);
					} 
					else if (ed.IsDouble())
					{
						ed.SetDouble(data.MakeDouble()); edc.Add(ed);
						//if (oo) mbs->UO() << "double data='" << data.MakeDouble() << ", " << ed.GetDouble() << "'\n";
					} 
					else if (ed.IsText())
					{
						ed.SetText(data.c_str()); edc.Add(ed);
					} 
					else
					{
						//error should not happen
					}
				}
			}
			else
			{
				LoadError(mystr("Specifier '")+dataname+mystr("' is unknown"));
				str.GetUntil(pos, '\n', data); //one line unknown keywords are no problem!!!
				//return 0;
			}
		}
		else
		{
			str = str.Right(str.Length()-pos); return 1;
		}


		if (pos == -1) endstr = 1;

	}
	return 1;
}

void CEDCParser::Initialize() 
{
	//here, things can be done during initialization of the CEDCparser
	el = "\n"; //end of line
	elc = '\n';
	obc = '{';  //open brace, for sub structure
	cbc = '}';  //closing brace
	ob = mystr(obc);
	cb = mystr(cbc);
	ofc = '(';  //open function (round '(') bracket
	cfc = ')';  //close function (round ')') bracket
	of = mystr(ofc);
	cf = mystr(cfc);
	osbc = '[';  //open square bracket, for vector or matrix
	csbc = ']';  //closing square bracket
	osb = osbc;  //open square bracket, for vector or matrix
	csb = csbc;  //closing square bracket
	ds = ": ";	//double stop - not used any more
	eq = "=";		//equal sign
	semicolon = ";"; //double stop
	semicolonc = ';'; //double stop
	commac = ','; //comma
	comma = commac; //comma
	doublequote = '"'; //for strings or text
	dotc = '.'; //dot
	dot = dotc; //dot
	space = ' '; //space
	spacec = ' ';
	tabc = '\t'; //tab
	tab = '\t';
	bool_yes = "yes";
	bool_no = "no";
	commentc = '%'; //comment
	comment = commentc; //comment
	matrixsep = el; //may be changed to semicolon ...

	//intialize switches and others:
	oo = 0;
	line_cnt = 0;
	total_lines = 0;

	version_check_performed = 0;
	
	InitializeCommands();
}

//read spaces, then comment sign, parse comment string into Tooltiptext of ed
//return 1 if success, 0 if fail
int CEDCParser::ReadComment(mystr& str, int& pos, ElementData& ed, const mystr& dataname)
{
	mystr data;
	int rv = str.GetUntilEOL(pos, commentc, data);
	if (rv == 2) //comment found
	{
		//mbs->UO() << "comment='" << data << "'+";
		if (str.PosPeek(pos) != commentc) 
		{
			mystr name = dataname;
			mbs->UO(UO_LVL_warn) << "Problem with comment in entry " << name <<". Check Textfile!\n";
			rv = 0;
		}

		// comment --> ToolTipText (double)
		mystr dummycomment;
		str.GetUntilEOL(pos, (char)0, dummycomment);
		mystr comment = dummycomment.SubString(1, dummycomment.Length()-1);
		ed.SetToolTipText(comment);
		//mbs->UO() << "'" << comment << "'\n";
	}

	return rv;
}

int CEDCParser::ReadBoolIntDoubleExpression(mystr& text, int& pos, ElementData& ed, const mystr& elemname)
{
	int rv = 1;
	text.EraseSpaces();
	if (PrintOutput()) mbs->UO() << "single text='" << text << "'\n";
	if (text[text.Length()-1] == ';') //accept ";" at end of statement
	{
		text.EraseChar(text.Length());
	}

	if(text.IsValidNumber())
	{
		ed.SetDouble(text.MakeDouble(), elemname); 
	}
	else // "expression"
	{
		if (PrintOutput()) mbs->UO() << "text='" << text << "' is no valid number\n";
		//bool variable
		if (bool_no.Compare(text))
		{
			ed.SetBool(0, elemname);
		}
		else if (bool_yes.Compare(text))
		{
			ed.SetBool(1, elemname);
		}
		else 
		{
			// DR 2012-08-02
			// AD 2013-10-16: collected calls to parse command
//AD: new procedure for commands
			mystr word = text; // "left hand side expression", not entire "str" and "pos" as for right hand side call
			int pos_tmp=0;

			mystr dataname = word.GetIdentifier(pos_tmp, 0); //name of data structure or data variable
			int com_id = IdentifyAndExecuteCommand(dataname, word, pos_tmp, ed, 1/*right hand side call*/);
			ed.SetDataName(elemname);

			if (com_id == -1)
			{ 
				rv = 0; // Parsing error occurred
			}
			else if (com_id > 0)
			{
				//rv = 1; // command executed
			}
			else //if (com_id == 0)
// NOT A COMMAND
			{
//!AD: 2012-12-11 assign entire EDCs: [
				// check if text is an EDC
// this would be already available in 
// CParseObj* CMBSParser::GetVariable(mystr& str)
// ElementDataContainer* CMBSParser::GetElementDataContainers(int i)
				int found_as_edc = 0;
				int found_as_entry = 0;
				////while (i <= NElementDataContainers() && ! found)
				////{
				////	ElementDataContainer* edc = GetElementDataContainers(i);
				
				//ElementDataContainer* edc_existing = this->GetMBS()->GetModelDataContainer();  // GLOBAL
				//if (edc_existing && edc_existing->TreeFind(text))
				//{
				//	ElementData* ed_other = edc_existing->TreeFind(text);
				//	if (ed_other->IsEDC())
				//	{
				//		found_as_edc = 1;
				//		ElementDataContainer* edc_new = new ElementDataContainer();
				//		edc_new->TreeReplaceEDCDataWith(ed_other->GetEDC());  // use TreeReplaceEDCDataWith to add the entries not clear the old ones
				//		ed.SetEDC(edc_new,elemname);
				//	}
				//}

				if (!(found_as_edc || found_as_entry))	
				{
			// check if RightHandSide exists as entry in LOCAL EDC first
					ElementDataContainer* localEDC = this->GetMBSParser()->GetLocalEDC();  // local context of the existing EDC
					if(localEDC && localEDC->TreeFind(text))
					{
						ElementData* ed_other = localEDC->TreeFind(text);
						if(ed_other->IsEDC())
						{
							found_as_edc = 1;
							ElementDataContainer* edc_new = new ElementDataContainer();
							edc_new->TreeReplaceEDCDataWith(ed_other->GetEDC());  // use TreeReplaceEDCDataWith to add the entries not clear the old ones
							ed.SetEDC(edc_new,elemname);
						}
						else
						{
							found_as_entry = 1;
							ed.CopyFrom(*ed_other);
							ed.SetDataName(elemname);
						}
					}
				}
				if (!(found_as_edc || found_as_entry))
				{
			// check if RightHandSide exists as entry in GLOBAL EDC first
					ElementDataContainer* globalEDC = this->GetMBSParser()->GetLocalEDC2();  // ROOT of the existing EDC 
					if(globalEDC && globalEDC->TreeFind(text))
					{
						ElementData* ed_other = globalEDC->TreeFind(text);
						if(ed_other->IsEDC())
						{
							found_as_edc = 1;
							ElementDataContainer* edc_new = new ElementDataContainer();
							edc_new->TreeReplaceEDCDataWith(ed_other->GetEDC());  // use TreeReplaceEDCDataWith to add the entries not clear the old ones
							ed.SetEDC(edc_new,elemname);
						}
						else
						{
							found_as_entry = 1;
							ed.CopyFrom(*ed_other);
							ed.SetDataName(elemname);
						}
					}
				}

				////	i++;
				////}
//!AD: ]
				if (!(found_as_edc || found_as_entry))
				{
					int errorflag;
					//////JG 2011-03-24
					//////Expression parser, using MBS_EDCs (edc_modeldata, mbs_edc_options, mbs_edc_global_variables) as variables
					////double value = ExpressionToDouble(mbs, text, errorflag, line_cnt);
					////	if (errorflag)
					////{
					////	ed.SetDouble(0, elemname); 
					////	rv = 0;
					////}
					////else
					////{
					////	ed.SetDouble(value, elemname); 
					////	rv = 1;
					////}	

//AD 2012-12-10
					//Expression parser extended to Vector and Matrix variables
					double d_value = 0;
					Vector v_value(0.);
					Matrix m_value(0.);
					int2 rv_rowscols = mbsParser->EvaluateExpression(text, errorflag, d_value, v_value, m_value, line_cnt);

					if (rv_rowscols == int2(1,1)) // scalar
					{
						if (errorflag) { ed.SetDouble(0, elemname); rv = 0; }
						else { ed.SetDouble(d_value, elemname); rv = 1; }			
					}
					else if (rv_rowscols(1)==1 || rv_rowscols(2)==1) // vector
					{
						if (errorflag) { ed.SetVector(v_value.GetVecPtr(), 1, elemname); rv = 0; }
						else { ed.SetVector(v_value.GetVecPtr(), v_value.Length(), elemname); rv = 1;}
					}
					else // matrix
					{
						if (errorflag) { ed.SetMatrix(m_value.GetMatPtr(), 1, 1, elemname); rv = 0; }
						else { ed.SetMatrix(m_value.GetMatPtr(), m_value.Getrows(), m_value.Getcols(), elemname); rv = 1; }
					}
				}
			}
		}
	}
	return rv;
}

//read a matrix or vector into ElementData
int CEDCParser::ReadMatrixVector(const mystr& text, ElementData& ed, const mystr& elemname, const mystr& dataname)
{
	int rv = 1;

	// col seperators : " "     ","
	// row seperators : "\n"    ";"
	TArray<double> ta_vec;
	mystr string_value;       // string containing the next value
	mystr string_entirerow;   // string containing the processed row 
	int pos_value = 0;        // position of next value to read
	int pos_in_string = 0;    // position to which the string is read
	int pos_in_row = 0;				// position to which the substring (row) is read
	int found_rowsep = 0;			// position of next row separator
	int found_colsep = 0;			// position of next column separator
	int rows=0; int cols=0; 	// number of rows, columns of the matrix 

	// check which separation chars are used
	text.ReadLeadingSpacesAndCountLines(pos_in_string, line_cnt);  // !skip spaces at the beginning of text to prevent pos_space==0 (for the case that ' ' is column separator)
	int pos_comma = text.Find(pos_in_string,commac);
	int pos_space = text.Find(pos_in_string,spacec);
	int pos_tab = text.Find(pos_in_string,tabc);
	int pos_semi = text.Find(pos_in_string,semicolonc);
	int pos_newline = text.Find(pos_in_string,elc);

	char colsep, rowsep;
	if ((pos_comma==-1)&&(pos_space==-1)) { colsep = commac; } // separation will not be found since this is a column vector
	if ((pos_comma > 0)&&(pos_space==-1)) { colsep = commac; } // separation with ','
	if ((pos_comma==-1)&&((pos_space > 0) || (pos_tab > 0))) { colsep = spacec; } // separation with ' '
	if ((pos_comma > 0)&&((pos_space > 0) || (pos_tab > 0))) { colsep = commac; } //mbs->UO(UO_LVL_warn) << "String2TreeEDC found column separators \",\" and \" \" in entry " << dataname <<". Check Textfile!\n"; }

	if ((pos_semi==-1)&&(pos_newline==-1)) { rowsep = semicolonc; } // separation will not be found since this is a row vector
	if ((pos_semi > 0)&&(pos_newline==-1)) { rowsep = semicolonc; } // separation with ';'
	if ((pos_semi==-1)&&(pos_newline > 0)) { rowsep = elc; }        // separation with '\n'
	if ((pos_semi > 0)&&(pos_newline > 0)) { rowsep = semicolonc; } // mbs->UO(UO_LVL_warn) << "String2TreeEDC found row separators \";\" and \"\\n\" in entry " << dataname <<". Check Textfile!\n"; }      

	while(pos_in_string != -1) // get -1 at end of text
	{

		rows++;
		found_rowsep = text.GetUntil(pos_in_string, rowsep, string_entirerow);
		string_entirerow.Replace(tab, space);

		int ccols = 0;	
		pos_in_row = 0;
		string_entirerow.EraseSpacesHeadTail();
		string_entirerow.ReadLeadingSpacesAndCountLines(pos_in_row, line_cnt); // !skip spaces at the beginning of a new row (for the case that ' ' is column separator)

		while(pos_in_row != -1 && rv) // get -1 at end of row
		{
			ccols++;
			string_entirerow.ReadLeadingSpaces(pos_in_row);
			found_colsep = string_entirerow.GetUntil(pos_in_row, colsep, string_value);

			if(string_value.IsValidNumber())
			{
				ta_vec.Add(string_value.MakeDouble());
			}
			else
			{				

				int errorflag;
				double value = mbsParser->ExpressionToDouble(string_value, errorflag, line_cnt);

				if (errorflag)
				{
					ta_vec.Add(0);
					rv = 0;
				}
				else
				{
					ta_vec.Add(value);
					rv = 1;
				}
			}

			if(found_colsep) pos_in_row++; // next value
			else // check consistency
			{
				if(rows == 1) 
				{
					cols = ccols;
				}
				else 
				{
					mystr name = dataname;
					if (cols!=ccols) 
					{
						mbs->UO(UO_LVL_warn) << "String2TreeEDC found mismatch in number of columns in entry " << name <<". Check Textfile!\n"; 
						rv = 0;
					}
				} 
			}
		}
		if(found_rowsep) {pos_in_string++;} // next row

	}
	Vector vec(ta_vec);
	if (rows <= 1) //if rows == 0, then assume vector
	{
		ed.SetVector(vec.GetVecPtr(), vec.Length(), elemname.c_str());  // vector
	}
	else 
	{
		ed.SetMatrix(vec.GetVecPtr(), rows, cols, elemname.c_str()); // matrix
	}

	return rv;
}


//####################################################################
//b: STANDARD VERSION of ParseAndExecuteString (f.k.a. String2TreeEDC)
// new version for parsing a string into an ElementDataContainer inclusive ToolTipText
// 6.9.10 added parameter treestr: name of tree structure for recursive calls (parent history...)
// 11.5.11 RL, AD: String2TreeEDC calls following recursive function "String2TreeEDC_recursive"
//                 Please use only "String2TreeEDC"!!!
// **.10.13 AD: functionality now only in ParseAndExecuteString
//              !! removed all calls to String2TreeEDC, removed implementation !!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//definition of EDC-structure:
// '%' | [identifier												 //comment or identifier appears
// if identifier is not a function:
//   ['{' | '=' | '(' | '%']                 //after identifier, either a structure opens '{' or a variable is assigned '=', a function is called '(' or a comment appears '%'
//   '{' ==> ['{ recursive_EDC }']
//   '=' ==> ['"' | '[...]' | [+-.1234567890] | mathobject/variable/function                        //read assignement: either string starts '"', vector/matrix '[...]', number, sin(2*pi*x)
//   ['%']																	 //comment might appear
//   ']'                                     //end of identifier structure
// if identifier is a function/command
//   '(function_parameters)'                 //read function parameters in paranthesis ()
//   '%'																		 //comment might appear
//   '{' ==> ['{ recursive_EDC }']           //if function supports a statement block (e.g. if/for/while/etc), then recursively read edc and execute function e.g. for-loop
// %
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// AD 2013-09-24: these functions parses a string ans applies it ENTRY BY ENTRY to the target element data container
// this function should replace the function StringToTreeEDC for the Commands "include", "for", etc...
// code is basically copied from function StringToTreeEDC
int CEDCParser::ParseAndExecuteString(mystr str, ElementDataContainer& edc)
{
// statistics on the string
	line_cnt = 0; //line counter in EDC
	total_lines = str.CountLines();
	mbs->UO(UO_LVL_dbg1) << "EDC total lines=" << total_lines << "\n";

// remember previous local EDC
	ElementDataContainer* prevlocaledc = mbsParser->GetLocalEDC();

// Set local variable context for the parser
	mbsParser->SetLocalEDC(&edc); //parsed (root-)edc is stored in local edc 

// cal the recursive part - the actual parsing
	int rv = ParseAndExecuteString_recursive(str, edc); 
	
// set local EDC back to previous 
	mbsParser->SetLocalEDC(prevlocaledc);

	return rv;
}

int CEDCParser::ParseAndExecuteString_recursive(mystr str, ElementDataContainer& edc)
{
	int rv = 1; //success
	int endstr = 0;
	int pos = 0;
	mystr dataname; // full name, may contain (sub-)tree
	mystr treename; // tree only
	mystr elemname; // varname only
	char ch;
	ElementData ed, return_value;
	mystr dummycomment; //dummy comment string
	
	
	int is_first_command = 1;
	TArray<ElementDataContainer*> parameter_EDCs;

// main loop
	while (!endstr && rv != 0)
	{
		ch = str.ReadLeadingSpacesAndCountLines(pos, line_cnt);
		
// skip lines that begin with comments  (e.g. file header, ...)
		while (ch == commentc && pos != -1)                               
		{
			str.GetUntilEOL(pos, (char)0, dummycomment);
			ch = str.ReadLeadingSpacesAndCountLines(pos, line_cnt);
		}

		int command_found = 0;
		if (pos == -1) //regular end of structure!
		{
			endstr = 1;
		}

//$AD 10-2013: new procedure for commands
		else
		{
			dataname = str.GetIdentifier(pos, 0); //name of data structure or data variable
			int com_id = IdentifyAndExecuteCommand(dataname, str, pos, ed, 0/*left hand side call*/);

			if (com_id == -1)
			{ 
				rv = 0; // Parsing error occurred
			}
			else if (com_id > 0)
			{
				command_found = 1;
				//rv = 1; // command executed
			}
			else //if (com_id == 0)
// NOT A COMMAND		
			{
				// split into tree name and element name for assignment without full tree structure
				edc.SplitIntoTreeAndElementName(dataname,treename,elemname); 
			}
		}
//$AD: end new procedure for commands

		if (pos == -1) //regular end of structure!
		{
			endstr = 1;
		}
		else if (command_found)
		{
			//do nothing, command already processed
		}

// BEGIN PARSING BEYOND DATANAME
		else
		{
			if (PrintOutput()) mbs->UO() << "Dataname='" << dataname << "'\n";

			ch = str.ReadLeadingSpacesAndCountLines(pos, line_cnt);

			while (ch == commentc && pos != -1)
			{
				str.GetUntilEOL(pos, (char)0, dummycomment);
				ch = str.ReadLeadingSpacesAndCountLines(pos, line_cnt);
			}

// FOUND "=":  ASSIGNMENT of various types
			if (ch == eq)
			{
				//read data:
				pos++; //go after eq sign
				ch = str.ReadLeadingSpacesAndCountLines(pos, line_cnt);
				
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				if (ch == osb) //vector/matrix
				{
					mystr text = str.GetStringInBrackets(pos, osbc, csbc);					
					ReadMatrixVector(text, ed, elemname, dataname);
					ReadComment(str, pos, ed, dataname);
					edc.TreeSet(treename,ed);//$JG2012-01-27: old if/else cases replaced //$ RL 2011-6-28: 
				}
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//string expected
				else if (ch == doublequote)
				{
					mystr text = str.GetStringInBrackets(pos, doublequote, doublequote);
					ed.SetText(text, elemname);
					ReadComment(str, pos, ed, dataname);
					edc.TreeSet(treename,ed);//$JG2012-01-27: old if/else cases replaced //$ RL 2011-6-28: 
				}
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//int, double, bool or regular expression
				else //must be double or integer number or "expression"
				{
					mystr text;
					int rv2 = str.GetUntilEOL(pos, commentc, text);
					if (rv2 == 2) {pos--;} //go back one character in order to find comment character in ReadComment(...) function hereafter
// EVALUATE RIGHT HAND SIDE OF THE ASSIGNMENT ( could include COMMANDS )
					rv = ReadBoolIntDoubleExpression(text, pos, ed, elemname);

					ReadComment(str, pos, ed, dataname);
// ADD OR REPLACE ENTRYs					
					edc.TreeSet(treename,ed); 
				}
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LINE IS PARSED - APPLY TO THE MODEL EDC
			}
			else if (ch == ob)
// FOUND "{": SUB TREE - RECURSION 
			{
				//read sub item
				mystr substr = str.GetStringInBrackets(pos, obc, cbc);
				if (pos != -1)
				{
					// check if tree already exists
					int i = edc.Find(elemname);
					if ((i != 0)&&(edc.Get(i).IsEDC())) // already exists as sub-edc
					{
						mbsParser->SetLocalEDC(edc.Get(i).GetEDC()); //$ RL 2011-5-11: //$ AD 2011-5-11: set local edc.
						rv = ParseAndExecuteString_recursive(substr, *(edc.Get(i).GetEDC()));
						mbsParser->SetLocalEDC(&edc); //$ RL 2011-5-11: //$ AD 2011-5-11: set local context edc back to value before recursion
					}
					else
					{
						ElementDataContainer* edc_sub = new ElementDataContainer(); //this is deleted when new EDC is set or EDC is deleted (old entry is always deleted)
						ed.SetBool(0,"dummy");
						int n = edc.Add(ed);
						edc.Get(n).SetEDCno_copy(edc_sub, elemname);
						mbsParser->SetLocalEDC(edc_sub); //$ RL 2011-5-11: //$ AD 2011-5-11: store sub-edc in local parser edc.
						rv = ParseAndExecuteString_recursive(substr, *edc_sub);
						mbsParser->SetLocalEDC(&edc); //$ RL 2011-5-11: //$ AD 2011-5-11: set local context edc back to value before recursion
					}
				}
				else
				{
					mbs->UO(UO_LVL_err) << "ParseAndExecuteString_recursive: could not read item '" << dataname << "' !\n";
					return 0; //end of file reached, not valid!
				}
			}
// ??? ADD LEFT SIDE LIST ACCESS RIGHT HERE ???
// DATANAME WAS FOLLOWED BY INVALID CHARACTER
			else
			{
				mbs->UO(UO_LVL_warn) << "did not expect '" << mystr(ch) << "' at line " << line_cnt+1 << " !\n";
				pos++;
				return 0; //end of file reached, not valid!
			}


		}
	} // end main loop
	return rv;
}

//e: STANDARD VERSION of ParseAndExecuteString
//############################################


// AD 2013-10-16: single call for commands - collect code from String2TreeEDC and ReadBoolIntDoubleExpression - returns "com_id" (including 0 for no command) or "-1" on error 
int CEDCParser::IdentifyAndExecuteCommand(mystr& dataname, mystr& str, int& pos, ElementData& return_value, int flag_RHS)
{
	int n_params, option, has_return_value;
 	int com_id = IdentifyCommand(dataname, n_params, has_return_value, option, flag_RHS);

	if (com_id > 0) // command found
	{
		if(!version_check_performed)
		{
			ElementDataContainer* edc_local = mbsParser->GetLocalEDC();
			HotintDataFileVersionCheck(version_check_performed, *edc_local);
		}

		TArray<ElementDataContainer*> parameter_EDCs;
		int parse_success = ParseCommandParameter(parameter_EDCs, dataname, n_params, str, pos, option);

		if (parse_success)
		{
// COMMAND IS IMMEDIATELY EXECUTED
			int rv = ExecuteCommand(com_id, parameter_EDCs, return_value, option);
			if(rv = 0)
			{
				com_id = -1;    // error occurred during execute command
			}
		}
		else
		{
			com_id = -1;      // error occurred during parsing of the command parameters
		}

		for(int np=1;np<=parameter_EDCs.Length();np++)
		{
			delete parameter_EDCs(np); // delete edc's
		}
	}
	return com_id;
}

// AD 2013-10-16: identifies command type and checks restrictions - returns "com_id" (including 0 for no command) or "-1" on error
int CEDCParser::IdentifyCommand(mystr& dataname, int& n_params, int& has_return_value, int& option, int flag_RHS )
{
	mystr handlename;				// treename
	mystr commandname;			// elementname
	ElementDataContainer* edc_local = mbsParser->GetLocalEDC();
	edc_local->SplitIntoTreeAndElementName(dataname, handlename, commandname);

// find registered command
  int com_id = GetCommandType(commandname, n_params, has_return_value, option);

	if(com_id > 0) 
// perform consistency checks
	{
// check 1: some commands are left hand side only (IF, FOR, ..)
		if ((option & TCO_LeftHandSideOnly) && flag_RHS)
		{
			EDCError(mystr("Right hand side call of Function '")+commandname+mystr("' is not allowed!"));
			return -1;
		}

// check 2: commands requiring a handle
		else if (option & TCO_RequiresHandle)
		{
	// call of a handle command is "HANDLE.COMMAND" so HANDLE must be registered in the current EDC
			if (handlename.Length())
			{
				int i = edc_local	->Find(handlename);
				if ((i != 0) && edc_local->Get(i).IsEDC())
				{
	// check if the handle is indeed a handle-object
					ElementDataContainer* edc_handle = edc_local->Get(i).GetEDC();
					int j = edc_handle->Find("isHandle");
					if (j==0)
					{
						EDCError(mystr("Handle '")+handlename+mystr("' for command '")+commandname+mystr("' is not a valid handle '"));
						return -1;
					}
					else 
					{
						j = edc_handle->Find("Functions");
						if ((j!=0) && edc_handle->Get(j).IsEDC())
						{
	// check if the command is registered to the handle
							ElementDataContainer* edc_registered = edc_handle->Get(j).GetEDC();
							int k = edc_registered->Find(commandname);
							if (k!=0)
							{
								;
							}
							else
							{
								EDCError(mystr("Function '")+commandname+mystr("' is not registered for Handle '")+handlename);
								return -1;
							}
						}
						else
						{
							EDCError(mystr("Handle '")+handlename+mystr("' has no allowed functions'"));
							return -1;
						}
					}
				}
			}
		}
// check 2b: when a handle is present and 
		else //(if !(option & TCO_RequiresHandle))
		{
			if (handlename.Length())
			{
				int i = edc_local	->Find(handlename);
				if ((i != 0) && edc_local->Get(i).IsEDC())
				{
	// check if the handle is indeed a handle-object
					ElementDataContainer* edc_handle = edc_local->Get(i).GetEDC();
					int j = edc_handle->Find("isHandle");
					if (j !=0)
					{
						EDCError(mystr("Function '")+commandname+mystr("' does not require a handle but is associated with one ('")+handlename+mystr("')"));
						return -1;
					}
				}
			}
		}
	}
	return com_id;
}

// AD 2013-10-16: parses parameters into EDCs - returns "1" on success and "0" on error
int CEDCParser::ParseCommandParameter(TArray<ElementDataContainer*>& parameter_EDCs, mystr& dataname, int n_params, mystr& str, int& pos, int option)
{
	int rv = 1;
	parameter_EDCs.SetLen(0);

//read parameter text
	mystr param_text = str.GetStringInBrackets(pos, ofc, cfc);
	if (param_text == mystr("") && n_params != 0)
	{
		EDCError(mystr("Expected ")+mystr(n_params)+mystr(" function parameters but received zero parameters in function '")+dataname+mystr("'"));
		rv = 0;
	}
	else 
	{	
		int prv = GetParameterEDC(param_text, parameter_EDCs); // fill parameter_EDCs with NEW edc's
		if (!prv)
		{
			EDCError(mystr("Function parameters '")+param_text+("' of function '")+dataname+mystr("' are not valid or incomplete"));
			rv = 0;
		}
		if (parameter_EDCs.Length() != n_params)
		{
			EDCError(mystr("Received ")+mystr(parameter_EDCs.Length())+(" function parameters, but function")+dataname+(" should have ")+mystr(n_params)+mystr(" parameters"));
			rv = 0;
		}
		else
		{
// additional parameter EDC contains (code-)container		
			if (option & TCO_FollowedByContainer)
			{
				while (str[pos] != obc) pos++;
				mystr codeincontainer = str.GetStringInBrackets(pos,obc,cbc);
				ElementData code;
				code.SetText(codeincontainer,"code");
				ElementDataContainer* codecontainer = new ElementDataContainer;
				codecontainer->Add(code);
				parameter_EDCs.Add(codecontainer);
			}
// additional parameter EDC contains the handle
			if (option & TCO_RequiresHandle)
			{
				ElementDataContainer* local = mbsParser->GetLocalEDC();
				mystr handlename, commandname;
				local->SplitIntoTreeAndElementName(dataname, handlename, commandname);

				//validity check only here
				int j = local->Find(handlename);
				if(j != 0 && local->Get(j).IsEDC())
				{
					//ElementDataContainer* edc_handle = local->Get(j).GetEDC();
					ElementData handle;
					handle.SetText(handlename,"handle");
					ElementDataContainer* handlecontainer = new ElementDataContainer;
					handlecontainer->Add(handle);
					parameter_EDCs.Add(handlecontainer);
				}
				else
				{
					EDCError(mystr("Handle '")+handlename+mystr("' could not be linked to Command '")+commandname+mystr("'. call only in local context!"));
					rv = 0;
				}
			}
		}
	}
	return rv;
}

//execute a command with type_id, a list of EDCs for command parameters and an option
int CEDCParser::ExecuteCommand(int command_type_id, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	return commands(command_type_id)->ExecuteCommand(GetMBS(), this, parameter_EDCs, return_value, option);
}

// AD 2013-10-16: perform a check for the datafile version (once)	
int CEDCParser::HotintDataFileVersionCheck(int& version_check_performed, ElementDataContainer& edc)
{
	version_check_performed = 1;
	mystr tmp_str = edc.TreeGetString("HOTINT_data_file_version");
	if(tmp_str.Length())
	{
		HotintVersionInfo ver_file(tmp_str);		// version in the input file
		HotintVersionInfo ver = mbs->GetHotintVersion();	
		int ver_log = ver.GetLogNumber();
		int ver_file_log = ver_file.GetLogNumber();
		if((ver > ver_file) || ((ver == ver_file) && (ver_log > ver_file_log)))
		{
			mbs->UO(UO_LVL_warn) << "WARNING: The input file was created with an old version of HOTINT. This may be the reason for errors. The actual version of HOTINT is '" << mbs->GetHotintVersion().GetString() << "'. To avoid this warning, check the functionality of the input file and correct the HOTINT_data_file_version in the input file afterwards.\n";
		}
	}
	else
	{
		mbs->UO(UO_LVL_warn) << "WARNING: The input file does not contain the variable HOTINT_data_file_version which indicates, that it has probably been created with an old version of HOTINT. This may be the reason for errors. The actual version of HOTINT is '" << mbs->GetHotintVersion().GetString() << "'. To avoid this warning, check the functionality of the input file and set the variable HOTINT_data_file_version in the input file right at the beginning.\n";
	}
	return 1;
}



// DR 2013-01-14: this function evaluates a command of the parser, using the parameter_EDC as additional input data
// AD 2013-10-16: remark: can be called from .cpp model for a) testing new script function, b) use existing script routines
int CEDCParser::ExecuteParserCommand(mystr & command, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	int pos = 0;
	mystr elemname = "element_name";
	ElementDataContainer *edc;
	ElementData ed;

	if(parameter_EDCs.Length())
	{
		edc = parameter_EDCs(1);
		ed.SetEDC(edc, "parameter_EDC");
		mbsParser->SetLocalEDC(edc);
	}
	else
	{
		ed.SetBool(1,"empty");
	}

	int rv = ReadBoolIntDoubleExpression(command, pos, ed, elemname);

	return rv;
}

//$ DR 2013-07-04:
// checks if filename has relative path, and if yes, converts it into an absolute one
// relative means relative to last_file_name and this is the file in which e.g. the command "Include" is used
int CEDCParser::MakeFilenameAbsolute(mystr& filename)
{
	int rv = 0; //1=success
	int pos = 1;
	int isRelativePath = filename.PosPeek(pos) != ':';
	if(isRelativePath)
	{		
		mystr rel_path(""); 
		mystr last_file_name = mbs->GetModelDataContainer()->TreeGetString("last_file_name");

		int until;
		for(until=last_file_name.Length()-1;until>=0;until--)
		{
			if(last_file_name.PosPeek(until) == '\\')
				break; // until <=> last backslash
		}
		for(int i=0;i<=until;i++)
		{
			rel_path += last_file_name.PosPeek(i); // path until last backslash
		}
		
		// remove ..\ 
		int pos=rel_path.Length()-1;
		if(rel_path.PosPeek(pos)=='\\')
		{
			rel_path.EraseChar(pos+1);
		}

		int replaced = filename.Replace(mystr("..\\"), mystr("")); //find string searchstr, replace with string replacestr; this is done repeatedly; the return value counts the number of replacements
	
		for(int pos=rel_path.Length();pos>0;pos--)
		{
			if(!replaced){break;}
			int prev = pos-1;
			if(rel_path.PosPeek(prev)=='\\')
			{
				replaced --;
			}
			rel_path.EraseChar(pos);
		}
		filename = rel_path + mystr("\\") + filename;
		rv = MakeFilenameAbsolute(filename);	// recursive check
	}
	else
	{
		rv = 1;	// no change --> success
	}
	if(rv<1)
	{
		//edcParser->EDCError(mystr("Error happened during CEDCParser::MakeFilenameAbsolute ");
	}
	return rv;
}