//#**************************************************************
//#
//# filename:             script_parser.h
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


#ifndef SCRIPT_PARSER__H
#define SCRIPT_PARSER__H

class CEDCParser;
class CMBSParser;
struct MBSObjectFactoryInterface;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Parsing input files:

//typedef enum {TAEflexible=1, TAEconstraint=2, TAE2D = 4, TAEspecial_connector=8} TAddElementType; //$ DR 2012-07: old code

////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////parser for ElementDataContainer:
////$JG2012-01-27: created
////string-->EDC
////EDC-->string

//these are default values for the creation of bodies, sensors, etc.
class Default_ParserObjectValues
{
public:
	Default_ParserObjectValues()
	{
		defaultbodycol = Vector3D(0.2,0.2,0.8);
		defaultconstraintcol = Vector3D(0.2,0.8,0.2);
		defaultbodypos = Vector3D(0.,0.,0.);
		defaultbodyvel = Vector3D(0.,0.,0.);
		defaultbodyangle = Vector3D(0.,0.,0.);
		defaultbodyangularvel = Vector3D(0.,0.,0.);
		defaultbodysize = Vector3D(1.,1.,1.);
		defaultbodyrho = 1000;
		defaultconstraintdim = Vector3D(0.1,0.1,16.);

		defaultsensorvisible = 1;
		defaultsensorwriteresults = 3;
		defaultsensorprecision = 17;
		defaultsensordrawdim = Vector3D(0.001,6,0);

		defaultgeomtile = 16;
		defaultgeomcol = Vector3D(0.5,0.5,0.5);
    defaultgeomdim = Vector3D(1,1,1); 
    defaulttransparency=-1;

		//flexible bodies:
		defaultbodyEm = 1e8;
		defaultbodynu = 0.;
	}
public:
	Vector3D defaultbodycol;
	Vector3D defaultconstraintcol;
	Vector3D defaultbodypos;
	Vector3D defaultbodyvel;
	Vector3D defaultbodyangle;
	Vector3D defaultbodyangularvel;
	Vector3D defaultbodysize;
	double defaultbodyrho;
	Vector3D defaultconstraintdim;

	int defaultsensorvisible;
	int defaultsensorwriteresults;
	int defaultsensorprecision;
	Vector3D defaultsensordrawdim;

	int defaultgeomtile;
	Vector3D defaultgeomcol;
	Vector3D defaultgeomdim; 
  double defaulttransparency;

	//flexible bodies:
	double defaultbodyEm;
	double defaultbodynu;
};

//typedef enum {TCT_AddElement=1, TCT_AddConstraint=2, TCT_AddSensor=3, TCT_AddLoad=4, TCT_AddGeomElement=5, } TEDCParser_CommandType;

typedef enum { TCO_None = 0, TCO_FollowedByContainer = 1, TCO_LeftHandSideOnly = 2, TCO_RequiresHandle = 4} TEDCParser_CommandOption;

//class for EDC-based commands, which can be parsed
class CEDC_Command
{
public:
	//initialization
	CEDC_Command(int id) 
	{
		command_id = id;
	}

	//unit command id, which is automatically assigned
	virtual int GetCommandID() const {return command_id;}

	//name of command in parser
	virtual const char* CommandName() const {assert(0); return "";}

	//Tex description for the docu
	virtual mystr GetCommandTexDescription() const {return "";};
	virtual mystr GetCommandTexParameterDescription() const {return "";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "";};
	virtual mystr GetCommandExampleFileName() const {return "";};

	//number of fixed parameters, which are necessary for the execution of the command
	virtual int NParameters() const {return 0;}

	virtual int HasReturnValue() const {return 0;} //0: no return value, 1: has a return value (ElementData)
	virtual int GetOptions() const {return 0;} //some options, which are not yet specified

	//execute command, 1==success, 0==failed
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option) {assert(0); return 0;}

	//convert textual (inline, not previously defined as variables) parameter to numerical value
	//NOTE: inline definition of Vector is NOT working: 
	//in function GetParameterEDC spaces get contracted, "," are interpreted as seperators for arguments
	virtual int EvaluateTextualParameterEntry(MBS* mbs, CEDCParser * edcParser, ElementDataContainer* param_edc); 

private:
	int command_id;
};

class CEDC_ComAddElement: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddElement(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddElement";}
	
	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "Adds an element to the system. See the description of the elements above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the element. ATTENTION: the entry element\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the element in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddElement.txt";};

	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return 0;}

	virtual int IsAddElement() const {return 1;} //this only distinguishes the AddElement from the AddConstraint command

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComAddConstraint: public CEDC_ComAddElement
{
public:
	//initialization
	CEDC_ComAddConstraint(int id): CEDC_ComAddElement(id) {}

	virtual const char* CommandName() const {return "AddConnector";}
	virtual int IsAddElement() const {return 0;} //this only distinguishes the AddElement from the AddConstraint command
	virtual mystr GetCommandTexDescription() const {return "Adds a connector to the system. See the description of the connectors above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the connector. ATTENTION: the entry element\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the connector in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddConnector.txt";};
};
//$ RL 2012-7-19:[  
class CEDC_ComAddGeomElement: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddGeomElement(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddGeomElement";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	
	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command adds an geometric element.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the geometric element. ATTENTION: the entry geom\\_element\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the geometric element in the MBS.";}
	virtual mystr GetCommandExampleFileName() const {return "AddGeomElement.txt";};


	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComInclude: public CEDC_Command
{
public:
	//initialization
	CEDC_ComInclude(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "Include";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 0;}


	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command includes a file.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is the absolut or relative filename. If a relative filename is used, than the path is relative to the last file! Be carefull, if you use this command more than one time in a file.";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "There is no return value defined yet.";}
	virtual mystr GetCommandExampleFileName() const {return "Include.txt";};


	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ RL 2012-7-19:] 

class CEDC_ComPrint: public CEDC_Command
{
public:
	//initialization
	CEDC_ComPrint(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "Print";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 0;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "Prints a text to the output window";};
	virtual mystr GetCommandTexParameterDescription() const;
	virtual mystr GetCommandTexReturnValueDescription() const {return "There is no return value for this command";};
	virtual mystr GetCommandExampleFileName() const {return "Print.txt";};


	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComAddLoad: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddLoad(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddLoad";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
	virtual mystr GetCommandTexDescription() const {return "Adds a load to the system. See the description of the loads above in order to get the available options. You have to adjust the value 'loads' in the element to assign the load to the element.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the load. ATTENTION: the entry load\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the load in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddLoad.txt";};
private:
};

//$ DR 2012-10
class CEDC_ComAddSensor: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddSensor(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddSensor";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
	virtual mystr GetCommandTexDescription() const {return "Adds a sensor to the system. See the description of the sensors above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the sensor. ATTENTION: the entry sensor\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the sensor in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddSensor.txt";};
private:
};

//$ DR 2012-12-11:[
class CEDC_ComLoadSTL: public CEDC_Command
{
public:
	//initialization
	CEDC_ComLoadSTL(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "ReadSTLFile";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}


	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command reads a stl-mesh from a file and stores the data in an ElementDataContainer.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is the absolut or relative filename.";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is an ElementDataContainer with 2 entries: triangles and points.";}
	virtual mystr GetCommandExampleFileName() const {return "ReadSTLFile.txt";};


	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ DR 2012-12-11:] 

class CEDC_ComAddMaterial: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddMaterial(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddMaterial";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
	virtual mystr GetCommandTexDescription() const {return "Adds a material to the system. See the description of the materials above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the material. ATTENTION: the entry material\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the material in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddMaterial.txt";};
private:
};

//$ DR 2013-01
class CEDC_ComAddBeam3DProperties: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddBeam3DProperties(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddBeamProperties";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
	virtual mystr GetCommandTexDescription() const {return "Adds a BeamProperty to the system. See the description of the BeamProperties above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the BeamProperties. ATTENTION: the entry material\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the node in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddBeamProperties.txt";};
private:
};

//$ DR 2013-01
class CEDC_ComAddNode: public CEDC_Command
{
public:
	//initialization
	CEDC_ComAddNode(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "AddNode";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
	virtual mystr GetCommandTexDescription() const {return "Adds a node to the system. See the description of the nodes above in order to get the available options.";};
	virtual mystr GetCommandTexParameterDescription() const {return "The parameter of this command is an ElementDataContainer with the data of the node. ATTENTION: the entry node\\_type must exist!";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value of this command is the number of the node in the MBS.";};
	virtual mystr GetCommandExampleFileName() const {return "AddNode.txt";};
private:
};


//$ DR 2013-01-30:[
class CEDC_ComComputeInertia: public CEDC_Command
{
public:
	//initialization
	CEDC_ComComputeInertia(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "ComputeInertia";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command computes the mass, moment of inertia, volume and center of mass based on the information about the geometry and the material of a body";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is an ElementDataContainer with 4 entries: volume, mass, moment\\_of\\_inertia and center\\_of\\_mass";}
	virtual mystr GetCommandExampleFileName() const {return "ComputeInertia.txt";};


	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ DR 2013-01-30:] 

//$ DR 2013-02-19:[
class CEDC_ComTransformPoints: public CEDC_Command
{
public:
	//initialization
	CEDC_ComTransformPoints(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "TransformPoints";}
	virtual int NParameters() const {return 3;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "With this command, the geometry described by the points can be transformed. It is possible to apply rotation and/or translation and/or scaling. The new point pN is computed according to the formula pN = trans + rot*p.";};
	virtual mystr GetCommandTexParameterDescription() const;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a Matrix containing the transformed points pN.";}
	virtual mystr GetCommandExampleFileName() const {return "TransformPoints.txt";};

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ DR 2013-02-19:] 

//$ DR 2013-07-04:[
class CEDC_ComLoadVecFromFile: public CEDC_Command
{
public:
	//initialization
	CEDC_ComLoadVecFromFile(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "LoadVectorFromFile";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command reads a vector from a file and returns this vector.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the vector.";}
	virtual mystr GetCommandExampleFileName() const {return "LoadVectorFromFile.txt";};

	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ DR 2013-07-04:] 

//$ AD 2013-07-11:[
class CEDC_ComSum: public CEDC_Command
{
public:
	//initialization
	CEDC_ComSum(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "Sum";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command adds two components of the same type ( scalar, vector or matrix ).";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the sum of the two inputs.";}
	virtual mystr GetCommandExampleFileName() const {return "Add_Mult.txt";};
	
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComProduct: public CEDC_Command
{
public:
	//initialization
	CEDC_ComProduct(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "Product";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command multiplies two components of the type ( scalar, vector or matrix ) when the operation is defined.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the product of the two inputs.";}
	virtual mystr GetCommandExampleFileName() const {return "Add_Mult.txt";};
	
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComTranspose: public CEDC_Command
{
public:
	//initialization
	CEDC_ComTranspose(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "Transpose";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command transposes a matrix or vector.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a matrix or a vector.";}
	virtual mystr GetCommandExampleFileName() const {return "Add_Mult.txt";};
	
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

//$! AD 2013-07-15 commented out Class Scalar product for the time being, functionality is in ComProduct...
////class CEDC_ComScalarProduct: public CEDC_Command
////{
////public:
////	//initialization
////	CEDC_ComScalarProduct(int id): CEDC_Command(id) {}
////
////	virtual const char* CommandName() const {return "ScalarProduct";}
////	virtual int NParameters() const {return 2;}
////	virtual int HasReturnValue() const {return 1;}
////
////	//description for the docu
////	virtual mystr GetCommandTexDescription() const {return "This command computes the scalar product of two vectors.";};
////	virtual mystr GetCommandTexParameterDescription() const ;
////	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the scalar product.";}
////	virtual mystr GetCommandExampleFileName() const {return "Add_Mult.txt";};
////	
////	//execute command, 1==success, 0==failed
////	//it is already checked, that the number of parameters is equal to the desired number of parameters
////	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
////private:
////};		

class CEDC_ComCrossProduct: public CEDC_Command
{
public:
	//initialization
	CEDC_ComCrossProduct(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "CrossProduct";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command computes the cross product of two vectors.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the scalar cross product.";}
	virtual mystr GetCommandExampleFileName() const {return "Add_Mult.txt";};
	
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ AD 2013-07-11:] 

//$ AD 2013-09-13:[ loop and condition
class CEDC_ComFor: public CEDC_Command
{
public:
	//initialization
	CEDC_ComFor(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "for";}
	virtual int NParameters() const {return 3;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_FollowedByContainer + TCO_LeftHandSideOnly;} // the command is followed by a code block in a { } container

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command starts a FOR loop for the subsequent block.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the number of iterations.";}
	virtual mystr GetCommandExampleFileName() const {return "Loop_Cond.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComIf: public CEDC_Command
{
public:
	//initialization
	CEDC_ComIf(int id): CEDC_Command(id) {}

	virtual const char* CommandName() const {return "if";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_FollowedByContainer + TCO_LeftHandSideOnly;} // the command is followed by a code block in a { } container

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command evaluates an IF condition for the subsequent block.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is the 1 for true and 0 for false.";}
	virtual mystr GetCommandExampleFileName() const {return "Loop_Cond.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ AD 2013-09-13:] 

//$ AD 2013-10-09:[ functions to generate specific handles
class CEDC_ComMesh_GenerateNewMesh: public CEDC_Command
{
public:
	//initialization
	CEDC_ComMesh_GenerateNewMesh(int id): CEDC_Command(id) {}
	
	virtual const char* CommandName() const {return "GenerateNewMesh";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command generates a Handle to a Mesh Object for further operations.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a special EDC (Handle) that must be assigned to a new variable.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ AD 2013-10-09:]

//$ AD 2013-10-09:[ Mesh commands
class CEDC_ComMesh_GenerateBeam: public CEDC_Command
{
public:
	//initialization
	CEDC_ComMesh_GenerateBeam(int id): CEDC_Command(id) {}
	
	virtual const char* CommandName() const {return "GenerateBeam";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_RequiresHandle;} 

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command generates a Beam within a Mesh Object.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a list of element numbers for the newly generated beam elements.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComMesh_GeneratePlate: public CEDC_Command
{
public:
	//initialization
	CEDC_ComMesh_GeneratePlate(int id): CEDC_Command(id) {}
	
	virtual const char* CommandName() const {return "GeneratePlate";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_RequiresHandle;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command generates a Plate within a Mesh Object.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a list of element numbers for the newly generated plate elements.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComMesh_GetNodesInBox: public CEDC_Command
{
public:
	//initialization
	CEDC_ComMesh_GetNodesInBox(int id): CEDC_Command(id) {}
	
	virtual const char* CommandName() const {return "GetNodesInBox";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_RequiresHandle;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command returns a list of nodes (registered to the mesh) in a given box.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "The return value is a list of node numbers.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

class CEDC_ComMesh_GlueMesh: public CEDC_Command
{
	public:
	//initialization
	CEDC_ComMesh_GlueMesh(int id): CEDC_Command(id) {}
	
	virtual const char* CommandName() const {return "GlueMesh";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_RequiresHandle;}

	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command glues mesh together.";};
	virtual mystr GetCommandTexParameterDescription() const {return "TODO.";};
	virtual mystr GetCommandTexReturnValueDescription() const {return "TODO.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};
//$ AD 2013-10-09:]

//$ AD 2013-11-27: command to determine if a specific 
class CEDC_Check_Entry_Exists: public CEDC_Command
{
	public:
	//initialization
	CEDC_Check_Entry_Exists(int id): CEDC_Command(id) {}
	virtual const char* CommandName() const {return "DoesEntryExist";}
	virtual int NParameters() const {return 1;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_None;}
	
	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command glues mesh together.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "returns 0 if the entry does not exist, returns 1 if the entry exists.";}
	virtual mystr GetCommandExampleFileName() const {return "Mesh.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

//$ AD 2013-12-13: command compare strings - could be used for other comparisons as well
class CEDC_Compare: public CEDC_Command
{
	public:
	//initialization
	CEDC_Compare(int id): CEDC_Command(id) {}
	virtual const char* CommandName() const {return "Compare";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_None;}
	
	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command compares two strings.";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "returns 0 if both strings are identical, returns >0 or <0 otherwise indicating which string has higher value .";}
	virtual mystr GetCommandExampleFileName() const {return "strings.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};

//$ AD 2013-12-13: command concatinate strings 
class CEDC_StrCat: public CEDC_Command
{
	public:
	//initialization
	CEDC_StrCat(int id): CEDC_Command(id) {}
	virtual const char* CommandName() const {return "StrCat";}
	virtual int NParameters() const {return 2;}
	virtual int HasReturnValue() const {return 1;}
	virtual int GetOptions() const {return TCO_None;}
	
	//description for the docu
	virtual mystr GetCommandTexDescription() const {return "This command joins two strings together";};
	virtual mystr GetCommandTexParameterDescription() const ;
	virtual mystr GetCommandTexReturnValueDescription() const {return "returns a single string - strA+strB.";}
	virtual mystr GetCommandExampleFileName() const {return "strings.txt";};
	//execute command, 1==success, 0==failed
	//it is already checked, that the number of parameters is equal to the desired number of parameters
	virtual int ExecuteCommand(MBS * mbs, CEDCParser * edcParser, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);
private:
};


class CEDCParser 
{
public:
	//initialize parser
	//CEDCParser(MBS * mbs);
	CEDCParser(): commands(), class_names(), class_descriptions()
	{
 		mbs = NULL;
		Initialize();
	};

	//$ YV 2012-12-30: expression parser is also stored as a member here
	virtual void SetMBS(MBS * mbs, CMBSParser * mbsParser)
	{
		this->mbs = mbs;
		this->mbsParser = mbsParser;
	}
	virtual const MBS * GetMBS() const {return mbs;}
	virtual MBS* GetMBS() {return mbs;}

	virtual void Initialize();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//auxiliary functions:
	virtual int ReadComment(mystr& str, int& pos, ElementData& ed, const mystr& dataname); //read comment at end of statement into ElementData
	virtual int ReadBoolIntDoubleExpression(mystr& text, int& pos, ElementData& ed, const mystr& elemname); //read expression containing bool (yes/no), int/double or a mathematical expression
	virtual int ReadMatrixVector(const mystr& text, ElementData& ed, const mystr& elemname, const mystr& dataname); //read structure in [] containing matrix or vector
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	virtual void EDCError(const mystr& str) {mbs->UO(UO_LVL_err) << "EDC Parser error: '" << str << "' at line " << line_cnt << "\n";}
	virtual void EDCWarning(const mystr& str) {mbs->UO(UO_LVL_err) << "EDC Parser warning: '" << str << "' at line " << line_cnt << "\n";}

	virtual int MakeFilenameAbsolute(mystr& filename);	// checks if filename has relative path, and if yes, converts it into an absolute one. returns 1 if success
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//commands, set up all possible commands:
	virtual void InitializeCommands() 
	{
		//add all EDC-commands manually here!
		int comcnt = 0;

		for(int i=1; i<=commands.Length(); i++) {delete commands(i);}	//$ DR 2012-10: clear array, because this function is called several times
		commands.Flush();

		CEDC_ComAddElement* cae = new CEDC_ComAddElement(++comcnt);
		commands.Add(cae);

		CEDC_ComAddGeomElement* cage = new CEDC_ComAddGeomElement(++comcnt);
		commands.Add(cage);

		CEDC_ComAddConstraint* cac = new CEDC_ComAddConstraint(++comcnt);		//$ DR 2012-08
		commands.Add(cac);

		CEDC_ComAddLoad* cal = new CEDC_ComAddLoad(++comcnt);		//$ DR 2012-10
		commands.Add(cal);

		CEDC_ComAddSensor* cas = new CEDC_ComAddSensor(++comcnt);		//$ DR 2012-10
		commands.Add(cas);

		CEDC_ComAddMaterial* cam = new CEDC_ComAddMaterial(++comcnt);		//$ DR+PG 2012-12
		commands.Add(cam);

		CEDC_ComAddBeam3DProperties* cabp = new CEDC_ComAddBeam3DProperties(++comcnt);		//$ DR+PG 2012-12
		commands.Add(cabp);

		CEDC_ComAddNode* can = new CEDC_ComAddNode(++comcnt);			//$ DR 2013-01
		commands.Add(can);

		CEDC_ComInclude* caincl = new CEDC_ComInclude(++comcnt);
		commands.Add(caincl);

		CEDC_ComPrint* cpr = new CEDC_ComPrint(++comcnt);		//$ DR 2012-11
		commands.Add(cpr);

		CEDC_ComLoadSTL* clstl = new CEDC_ComLoadSTL(++comcnt);		//$ DR 2012-12
		commands.Add(clstl);

		CEDC_ComLoadVecFromFile* clvff = new CEDC_ComLoadVecFromFile(++comcnt);			//$ DR 2013-07
		commands.Add(clvff);

		CEDC_ComTransformPoints* ctp = new CEDC_ComTransformPoints(++comcnt);		//$ DR 2013-02-19
		commands.Add(ctp);
		
		CEDC_ComComputeInertia* cci = new CEDC_ComComputeInertia(++comcnt);			//$ DR 2013-01
		commands.Add(cci);

		CEDC_ComSum* cs = new CEDC_ComSum(++comcnt);			//$ AD 2013-07
		commands.Add(cs);

		CEDC_ComProduct* cp = new CEDC_ComProduct(++comcnt);    //$ AD 2013-07
		commands.Add(cp);

		CEDC_ComTranspose* ct = new CEDC_ComTranspose(++comcnt);    //$ AD 2013-07
		commands.Add(ct);

		//CEDC_ComScalarProduct* csp = new CEDC_ComScalarProduct(++comcnt);    //$ AD 2013-07
		//commands.Add(csp);
		
		CEDC_ComCrossProduct* ccp = new CEDC_ComCrossProduct(++comcnt);    //$ AD 2013-07
		commands.Add(ccp);

		CEDC_ComFor* cfor = new CEDC_ComFor(++comcnt);    //$ AD 2013-09
		commands.Add(cfor);
		
		CEDC_ComIf* cif = new CEDC_ComIf(++comcnt);    //$ AD 2013-09
		commands.Add(cif);

		CEDC_ComMesh_GenerateNewMesh* cmgmesh = new CEDC_ComMesh_GenerateNewMesh(++comcnt);    //$ AD 2013-09
		commands.Add(cmgmesh);

		CEDC_ComMesh_GenerateBeam* cmgbeam = new CEDC_ComMesh_GenerateBeam(++comcnt);    //$ AD 2013-09
		commands.Add(cmgbeam);

		CEDC_ComMesh_GeneratePlate* cmgplate = new CEDC_ComMesh_GeneratePlate(++comcnt);    //$ AD 2013-10
		commands.Add(cmgplate);

		CEDC_ComMesh_GetNodesInBox* cmgetnodes = new CEDC_ComMesh_GetNodesInBox(++comcnt);    //$ AD 2013-10
		commands.Add(cmgetnodes);

		CEDC_ComMesh_GlueMesh* cmgluemesh = new CEDC_ComMesh_GlueMesh(++comcnt);    //$ AD 2013-10
		commands.Add(cmgluemesh);

		CEDC_Check_Entry_Exists* cmentryexists = new CEDC_Check_Entry_Exists(++comcnt);    //$ AD 2013-11
		commands.Add(cmentryexists);

		CEDC_Compare* cmcompare = new CEDC_Compare(++comcnt);    //$ AD 2013-12
		commands.Add(cmcompare);

		CEDC_StrCat* cmstrcat = new CEDC_StrCat(++comcnt);    //$ AD 2013-12
		commands.Add(cmstrcat);

	};

	//get the command type_id from the command_name; also provide the number of parameters and optional flags for the command
	//return 0, if failed
	virtual int GetCommandType(const mystr& command_name, int& n_parameters, int& has_return_value, int& option)
	{
		n_parameters = 0;
		has_return_value = 0;
		option = 0;
		for (int i=1; i <= commands.Length(); i++)
		{
			if (mystr(commands(i)->CommandName()) == command_name) //==>STringCMP would be more efficient
			{
				n_parameters = commands(i)->NParameters();
				has_return_value = commands(i)->HasReturnValue();
				option = commands(i)->GetOptions();
				return i;
			}
		}
		return 0;
	}

	virtual mystr GetCommandTypeName(int i) const {return mystr(commands(i)->CommandName());}
	virtual int NCommands() const {return commands.Length();}
	virtual CEDC_Command* GetCommand(int i) const {return commands(i);}

// AD 2013-10-16: single call for commands - collect code from String2TreeEDC and ReadBoolIntDoubleExpression - returns "com_id" (including 0 for no command) or "-1" on error 
	virtual int IdentifyAndExecuteCommand(mystr& dataname, mystr& str, int& pos, ElementData& return_value, int flag_RHS = 0);

// AD 2013-10-16: identifies command type and checks restrictions - returns "com_id" (including 0 for no command) or "-1" on error
	virtual int IdentifyCommand(mystr& dataname, int& n_params, int& has_return_value, int& option, int flag_RHS );

// AD 2013-10-16: parses parameters into EDCs - returns "1" on success and "0" on error
	virtual int ParseCommandParameter(TArray<ElementDataContainer*>& parameter_EDCs, mystr& dataname, int n_params, mystr& str, int& pos, int option);

//execute a command with type_id, a list of EDCs for command parameters and an option 
	virtual int ExecuteCommand(int command_type_id, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);

// AD 2013-10-16: perform a check for the datafile version (once)	
	virtual int HotintDataFileVersionCheck(int& version_check_performed, ElementDataContainer& edc_version);




	//get TArray<ElementDataContainer*> of comma-separated parameter sets for a command
	virtual int GetParameterEDC(mystr text, TArray<ElementDataContainer*>& parameter_EDCs);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//main EDC parsing functions:
	//transform edc into string, works already for TreeEDCs!
	virtual void EDC2String(const ElementDataContainer& edc, mystr& str, const mystr& indent = "");
	//virtual void EDC2TexTable(const ElementDataContainer& edc, mystr& str, const mystr& indent = "", const mystr& treename = ""); //$ DR 2013-01-11 moved to models_auto_documentation

	//transforms string into edc until a keyword is read, returns remaining string:
	virtual int String2EDC(/*MBS * mbs, */mystr& str, ElementDataContainer& edc, ElementDataContainer& testedc); 

	////////transforms a string automatically into a tree edc; only accepts: double, string (""), vector, matrix
	//////virtual int String2TreeEDC(/*MBS * mbs, */mystr str, ElementDataContainer& edc); 
	////////recursive version of string to TreeEDC conversion
	//////virtual int String2TreeEDC_recursive(mystr str, ElementDataContainer& edc);

  // AD 2013-09-24: these functions parses a string ans applies it ENTRY BY ENTRY to the target element data container
	virtual int ParseAndExecuteString(mystr str, ElementDataContainer& edc); 
	virtual int ParseAndExecuteString_recursive(mystr str, ElementDataContainer& edc);
	//virtual int AddReplaceSingleEntry(ElementData& entryED, ElementDataContainer& targetEDC);

	virtual Default_ParserObjectValues& DPOV() {return dpov;}
	virtual const Default_ParserObjectValues& DPOV() const {return dpov;}
	
	// DR 2013-01-14: this function evaluates a command of the parser, using the parameter_EDC as additional input data 
	// AD 2013-10-16: remark: can be called from .cpp model for a) testing new script function, b) use existing script routines
	virtual int ExecuteParserCommand(mystr & command, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option);

	// +++++++++++++ the classs descriptions are not used at the moment +++++++++++++++++++++++	DR 2013-01-11
	//virtual void AddClassDescriptions_Auto();
	//add a class name with a description to the parser
	virtual void AddClassDescription(const mystr& name, const mystr& description)
	{
		class_names.Add(name);
		class_descriptions.Add(description);
	}
	//get a description for a certain class name
	virtual mystr GetClassDescription(const mystr& name)
	{
		for (int i=1; i <= class_names.Length(); i++)
		{
			if (name == mystr(class_names(i))) return mystr(class_descriptions(i));
		}
		return "";
	}
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 	virtual bool PrintOutput() {
 		bool val = oo;
 		if (mbs)
 		{
 			if (mbs->GetOptions()->LoggingOptions()->EDCParserGeneralInformation() == 1)
 			{
 				val = true;
 			}
 			else
 			{
 				val = false;
 			}
 		}

 		return val;
 	}

public:
	MBSObjectFactoryInterface * GetObjectFactory() { return pObjectFactory; }
	void SetObjectFactory(MBSObjectFactoryInterface * pObjectFactory) { this->pObjectFactory = pObjectFactory; }

public: 
	CMBSParser* GetMBSParser() { return mbsParser; }

private:
	MBS * mbs;
	CMBSParser * mbsParser;
	MBSObjectFactoryInterface * pObjectFactory;		// a pointer to the instance in the models and elements module

	Default_ParserObjectValues dpov; //default parser object values

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//commands:
	TArray<CEDC_Command*> commands;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//$ DR 2013-01-11 not used at the moment:
	MyStrList class_names, class_descriptions; //this are two lists, which provide a detailed .tex based description for a single class

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//definitions of characters for parsing:
	mystr el; //end of line
	char elc;
	char obc;  //open brace, for sub structure
	char cbc;  //closing brace
	mystr ob;
	mystr cb;
	char ofc;  //open function (round '(') bracket
	char cfc;  //close function (round ')') bracket
	mystr of;
	mystr cf;
	char osbc;  //open square bracket, for vector or matrix
	char csbc;  //closing square bracket
	mystr osb;  //open square bracket, for vector or matrix
	mystr csb;  //closing square bracket
	mystr ds;	//double stop - not used any more
	mystr eq;		//equal sign
	mystr semicolon; //double stop
	char semicolonc; //double stop
	char commac; //comma
	mystr comma; //comma
	char doublequote; //for strings or text
	char dotc; //dot
	mystr dot; //dot
	mystr space; //space
	char spacec;
	char tabc;
	mystr tab;
	mystr bool_yes;
	mystr bool_no;

	char commentc; //comment
	mystr comment; //comment
	mystr matrixsep; //may be changed to semicolon ...

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//other variables and switches:
	int line_cnt; //line counter in EDC
	int total_lines; //total number of lines in string
	int oo; //turn on / off detailed output

	int version_check_performed;    // flag to perform a version check when first command is encountered

	void LoadError(const mystr& str) {mbs->UO() << "Error in load file: " << str.c_str() << "\n";}
};
#endif
