//#**************************************************************
//#
//# filename:             Sim2Hotint.h
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

#ifndef SIM2HOTINT__H
#define SIM2HOTINT__H

#include "graphicsConstants.h"

extern double String2Double(MBS * mbs, mystr& str_exp, ElementDataContainer* edc_par=0); //$ RL 2011-01: evaluate string str and make double; //$ RL 2011-6-3: parser with standard operations like +-*/ available, todo: subedc's
extern Vector String2Vector(MBS * mbs, mystr& str_exp, char osbc = '[', char csbc = ']', char commac = ',', ElementDataContainer* edc_par = 0); //$ RL 2011-01: evaluate string str and make vector
extern Element* ConvertSim2Hotint_Gain(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Inport(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Outport(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Sum(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Integrator(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Constant(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_Saturate(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);
extern Element* ConvertSim2Hotint_TransferFcn(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset); //$ RL 2011-01: transfer function added.
extern Element* ConvertSim2Hotint_DeadZone(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$ RL 2011-01: DeadZone implemented (e.g.: "clearance").
extern Element* ConvertSim2Hotint_Scope(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$ RL 2011-01: Scope is translated into MBSSensor.
extern Element* ConvertSim2Hotint_TransportDelay(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$!RL 2011-01: continuous delay is approximized as time discrete delay.

extern Element* ConvertSim2Hotint_Product(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$!RL 2011-05: time continuous block for multiplication and division.
extern Element* ConvertSim2Hotint_Fcn(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$!RL 2011-05: time continuous block for 1-dimensional function (input variable = u).
extern Element* ConvertSim2Hotint_Reference(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$!RL 2011-05: time continuous block.
extern Element* ConvertSim2Hotint_DiscretePulseGenerator(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset);//$!RL 2011-07: time continuous rectangular pulse generator.



typedef Element* (*Sim2Hotint_FNPTR)(MBS * mbs, ElementDataContainer* edc, ElementDataContainer* edc_default, double scalefactor, Vector2D drawoffset)  ; 
// add subsystem from simulink with inports and outports
extern void AddSimulinkSubSystem(MBS * mbs,mystr fname,TArray<int>& elnums_inputs,TArray<int>& elnums_outputs, Vector2D drawoffset=Vector2D(0.0), double scaling_factor=0.005);

class LineConnection
{

};

class Simulink2HotintConversion
{

public:
	Simulink2HotintConversion(MBS * mbsi)
	{
		mbs = mbsi;
		Initialize();
		scalefactor = 0.005;
		drawoffset = Vector2D(0.,0.);
		nrOfInports = 0;
		nrOfOutports = 0;
		edc_default = 0;
		//variable_edc_name = "";
		//edc_par = 0;
		//fname_par = "";
	}

	~Simulink2HotintConversion()
	{		
		for(int i=1;i<=sim2hotint_identifier_list.Length();i++)
		{
			delete sim2hotint_identifier_list(i);
		}
	
		for(int i=1;i<=edc_blocks.Length();i++)
		{
			delete edc_blocks(i);
		}

		for(int i=1;i<=edc_lines.Length();i++)
		{
			delete edc_lines(i);
		}
	}

	virtual void Initialize()
	{
		mystr* str;

		str = new mystr("Gain");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Gain); 

		str = new mystr("Inport");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Inport);


		str = new mystr("Outport");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Outport); 

		str = new mystr("Sum");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Sum); 

		str = new mystr("Integrator");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Integrator); 	

		str = new mystr("Constant");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Constant);

		str = new mystr("Saturate");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Saturate);


		str = new mystr("TransferFcn");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_TransferFcn);	
	
		str = new mystr("DeadZone");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_DeadZone);	

		str = new mystr("Scope");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Scope);	

		str = new mystr("TransportDelay");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_TransportDelay);	

		str = new mystr("Product");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Product);	
		
		str = new mystr("Fcn");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Fcn);

		str = new mystr("Reference");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_Reference);	

		str = new mystr("DiscretePulseGenerator");
		sim2hotint_identifier_list.Add(str);
		sim2hotint_fnptr_list.Add(*ConvertSim2Hotint_DiscretePulseGenerator);	
	}

	virtual int LoadSimulink2Hotint(const mystr& fname); //main function for creating the 
	//virtual void SetParameterFile(const mystr& fname_txtI)
	//{
	//	fname_par = fname_txtI;
	//}
	virtual int ReadSimulinkModel(const mystr& fname,ElementDataContainer* edc); //this function reads the model data from mdl-file into an edc

	//virtual void ParseSimulinkEDCModel_Default(const ElementDataContainer* edc); 
	virtual void ParseSimulinkEDCModel(ElementDataContainer* edc);  //this function parses the edc and creates connected InputOutputElements from edc-information
	virtual Vector3D GetColorVector(MBS * mbs, ElementDataContainer* edc_block, const char* colkey, const Vector3D& coldefault=colblack); //returns the color vector of an element
	virtual int GetRotation(MBS * mbs, ElementDataContainer* edc_block, const char* colkey , int defaultrot=0); //returns the numbers of 90° rotations
	virtual Element* ConvertBlock2Element(ElementDataContainer* edc_block); //create InputOutputElements element from information of edc_block
	virtual void parseBranch(ElementDataContainer* edc_branch, TArray<mystr*>& dstBlock, TArray<int>& dstPort, TArray<Vector2D>& dstPoints, TArray<int>& NdstPoints); //, TArray<int>& dstPointLevel); // parse keyword "Branch" in mdl-file and get destination(s) and destination line(s)
	virtual void ConnectBlocks(ElementDataContainer* edc_line); // connect two InputOutputElements due to the information in edc_line
	virtual double GetScaleFactor() {return scalefactor;} // get drawing scale factor
	virtual void SetScaleFactor(double scalefactori) {scalefactor=scalefactori;} // set scalefactor
	virtual Vector2D GetDrawOffset() {return drawoffset;} // get drawing offset
	virtual void SetDrawOffset(Vector2D& drawoffseti) {drawoffset=drawoffseti;} // set drawing offset
	virtual int GetInportMBSElementNumber(int inport_num);   //compute element number of input in mbs from input number
	virtual int GetOutportMBSElementNumber(int outport_num); //compute element number of output in mbs from output number
	virtual int GetNumberOfInports(){return nrOfInports;};   //returns number of inputs
	virtual int GetNumberOfOutports(){return nrOfOutports;}; //returns number of outputs

	
	//virtual void SetEDCName(mystr& variable_edc_namei){variable_edc_name = variable_edc_namei;}; //set data container name, where variables of model are stored. If no name is set, the variables are searched in the model-edc.

private:

	MBS * mbs;
	//ElementDataContainer* edc_par; // name of parameter edc, where variables of mdl-file are defined.
	//mystr fname_par; // name of parameter file with variables of mdl-file in container "MDL" e.g. MDL.a = 2

	TArray<mystr*> sim2hotint_identifier_list;      // names of input/output objects in mdl-file
	TArray<Sim2Hotint_FNPTR> sim2hotint_fnptr_list; // list of function pointers for creating control elements from block-edc

	ElementDataContainer* edc_default;              // edc to default values from mdl-file
	TArray<ElementDataContainer*> edc_blocks;       // blocks-information of input/output objects from mdl-file is stored in this edc
	TArray<ElementDataContainer*> edc_lines;        // line-information (connection of input/output objects) is stored in this edc
	//TArray<Element*>
	//TArray<LineConnection*>
	double scalefactor; // scalefactor for drawing 
	Vector2D drawoffset;// offset for drawing

	int nrOfInports; // number of inports/inputs
	int nrOfOutports;// number of outports/outputs

	//mystr variable_edc_name; // name of edc, where variables of model are stored

	int String2TreeEDC_SIM(mystr& str, ElementDataContainer& edc);

};
#endif