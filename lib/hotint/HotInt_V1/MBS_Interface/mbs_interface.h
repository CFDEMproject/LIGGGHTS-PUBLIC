//#***************************************************************************************
//# filename:				mbs_interface.h
//# authors:				Johannes Gerstmayr, Yuri Vetyukov
//# generated:    
//# rewritten:    
//# description:  
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

#pragma once

// this is an interface for the MultiBodySystem class;
// it contains functions, which may be needed by elements and model generating functions
// at the stages of model creation, computation and drawing;
// it is placed into the class hierarchy at the level of NumNLSys,
// as several functions, which belong to the interface,
// are defined there or in TimeInt
// in future, deeper decoupling is foreseen, so that we
// separate some of the parts of the common interface into base
// classes

#include "preprocessor_includes.h"    // Hotint preprocessor definitions
#include "mystring.h"
#include "tarray.h"
#include "..\WorkingModule\hotint_version_info.h"

// we need linear algebra everywhere - and prior to including the headers we need these definitions
class Matrix;
class MatrixXD;
class SparseMatrix;
class Vector3D;
class Vector2D;
class Matrix2D;
class Matrix3D;
class Vector;
class int3;
typedef TArray<int> IVector;
typedef TArray<double> DVector;
struct SuperMatrix;
#include "mathaux.h"
#include "linalg3d.h"
#include "linalg.h"

#include "UserOutputInterface.h"
#include "NumSolverInterface.h"

class Element;
class Sensor;
class Node;
class GeomElement;
class MBSLoad;
class Material;

class IntegrationRulesLibrary;
//class NumNLSolver;
struct FieldVariableDescriptor;
class ElementDataContainer;
struct RenderContext;
struct ControlWindowContext;
class HOTINTOptions;
class SolverSettings;
class TimeIntLog;

class SearchTree;

class CMBSParser;

//$ MaSch 2013-08-26
class TCPIPHotInt;

// old definitions
//for faster operations, parallel threads impossible!!!
#define mystatic static
//#define mystatic  
//#define gencnt

struct MBSSolutionAccessInterface
{
	// solution vector
	virtual const double& GetXact(int i) const = 0;
	virtual double& GetXact(int i) = 0;
	virtual double GetEigenValue(int index) const = 0;
	virtual const double& GetDrawValue(int i) const  = 0;
	virtual double& GetDrawValue(int i) = 0;
	virtual const Vector & GetLastSolVector() const = 0;
	virtual const Vector& GetSolVector() const = 0;
	virtual const Vector& GetXact() const = 0;		// probably this is the same as the function below
	virtual Vector& GetSolVector() = 0;
	virtual void SetActState(const Vector& x) = 0;
	virtual void SetActState(double* ptr, int len) = 0;
	virtual const Vector& GetLastNLItSolVector() const = 0;

	// extra data
	virtual const double& GetDataAct(int i) const = 0;
	virtual double& GetDataAct(int i) = 0;
	virtual const double& GetDataDraw(int i) const = 0;
	virtual double& GetDataDraw(int i) = 0;
	virtual Vector& GetLastDataVector() = 0;
	virtual const Vector& GetLastDataVector() const = 0;
	virtual Vector& GetLastNLItDataVector() = 0;
	virtual const Vector& GetLastNLItDataVector() const = 0;
};

struct MBSObjectsAccessInterface
{
	// elements
	virtual int AddElement(Element* e) = 0;
	virtual void DeleteElement(int i) = 0;
	virtual Element& GetElement(int i) = 0;
	virtual const Element& GetElement(int i) const = 0;
	virtual Element* GetElementPtr(int i) = 0;
	virtual const Element* GetElementPtr(int i) const = 0;
	virtual int NE() const = 0;
	virtual int GetNElements() const = 0;		// the same as above

	// constraints
	virtual int NConstraints() const = 0;	// DR 2013-01-15

	// auxiliary elements
	virtual Element& GetAuxElement(int i) = 0;
	virtual const Element& GetAuxElement(int i) const = 0;
	virtual Element* GetAuxElementPtr(int i) = 0;
	virtual const Element* GetAuxElementPtr(int i) const = 0;
	virtual int AddAuxElement(Element* e) = 0;
	virtual int NAuxE() const = 0;

	// sensors
	virtual int AddSensor(Sensor * s) = 0;
	virtual Sensor & GetSensor(int i) = 0;
	virtual const Sensor & GetSensor(int i) const = 0;
	virtual int NSensors() const = 0;
	virtual void DeleteSensor(int i) = 0;
	virtual Sensor* GetSensorPtr(int i) = 0;							//$ DR 2013-01-11 added
	virtual const Sensor* GetSensorPtr(int i) const = 0;  //$ DR 2013-01-11 added

	// loads
	virtual int AddLoad(const MBSLoad& li) = 0;
	virtual MBSLoad& GetLoad(int i) = 0;
	virtual const MBSLoad& GetLoad(int i) const = 0;
	virtual void DeleteLoad(int i) = 0;
	virtual int NLoads() const = 0;
	virtual MBSLoad* GetLoadPtr(int i) = 0;
	virtual const MBSLoad* GetLoadPtr(int i) const = 0;

	// nodes
	virtual void DeleteNode(int i) = 0;
	virtual int AddNode(Node* n) = 0;
	virtual int AddNode(Node* n, SearchTree& tree) = 0;
	virtual Node& GetNode(int i) = 0;
	virtual const Node& GetNode(int i) const = 0;
	virtual int NNodes() const = 0;
	virtual int AddBodyNode(Node* n) = 0;
	virtual int AddBodyNode(Node* n, SearchTree& tree) = 0;
	virtual Node* GetNodePtr(int i) = 0;
	virtual const Node* GetNodePtr(int i) const = 0;


	// materials
	virtual int AddMaterial(const Material& m) = 0;
	virtual Material& GetMaterial(int i) = 0;
	virtual const Material& GetMaterial(int i) const = 0;
	virtual int NMaterials() const = 0;
	virtual int DeleteMaterial(int i) = 0;
	virtual Material* GetMaterialPtr(int i) = 0;
	virtual const Material* GetMaterialPtr(int i) const = 0;


	// geometric (draw) elements
	virtual int Add(const GeomElement& de) = 0;
	virtual int Add(const GeomElement& de, int elnum) = 0;
	virtual GeomElement* GetDrawElement(int i) const = 0;
	virtual GeomElement* GetDrawElement(int i) = 0;
	virtual int NDrawElements() const = 0;
	virtual void DeleteDrawElement(int i) = 0;
};

struct MBS3DDrawingInterface
{
	virtual RenderContext * GetRC() = 0;
	//virtual ControlWindowContext * GetCWC() = 0;
	virtual void MyDrawLineH(const Vector3D& p1, const Vector3D&p2, const Vector3D& vy2, 
		double t, double h, int drawouterface=1) const = 0;
	virtual void DrawColorQuads(const TArray<Vector3D>& p, const TArray<double>& v, int n1, int n2, 
		int colormode, int drawlines=0, int vres=1) = 0;
	virtual void DrawQuad(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4) const = 0;
	virtual void DrawTrig(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3) const = 0;
	virtual void DrawHex(const Vector3D& p1,const Vector3D& p2,const Vector3D& p3,const Vector3D& p4,
		const Vector3D& p5,const Vector3D& p6,const Vector3D& p7,const Vector3D& p8, int drawouterfaces=1) const = 0;	
	virtual void DrawCube(const Vector3D& p0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3) const = 0;
	virtual void MyDrawRectangle(const Vector3D& p1, const Vector3D& p2, const Vector3D& p3, const Vector3D& p4, double thickness, const Vector3D* colline, const Vector3D* colfill=0) const = 0;
	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double thickness) const = 0;
	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double thickness, const Vector3D& col) const = 0;
	virtual void MyDrawLine(const Vector3D& p1, const Vector3D&p2, double t, double h) const = 0;
	virtual void MyDrawCircleXY(const Vector3D p, double r, const Vector3D& col, int res=12, double thickness=1) const = 0;
	virtual void MyDrawArrow(const Vector3D& p1, const Vector3D&p2, const Vector3D& col, 
		double linethickness = -1, double headsize = -1, int resolution=8) const = 0;
	virtual void DrawArrow(const Vector3D& p1, const Vector3D&p2, double linethickness = -1, double headsize = -1, int resolution = 8) const = 0;
	virtual void DrawColorArrow(double v, const Vector3D& p1, const Vector3D&p2, double linethickness = -1, double headsize = -1, int resolution = 8) = 0;	
	virtual void DrawColorZyl(double v, const Vector3D& pz1,const Vector3D& pz2, double rz, int tile) = 0;
	virtual void DrawZyl(const Vector3D& pz1, const Vector3D& pz2, double rz, int tile=8) const = 0;
	virtual void DrawCone(const Vector3D& pz1, const Vector3D& pz2, double rz, int tile=8, int drawconelines = 0) const = 0;
	virtual void DrawZyl(const Vector3D& pz1, const Vector3D& pz2, const Vector3D& pz1dir,const Vector3D& pz2dir, double rz, int leftend = 1, int rightend = 1, int tile=8) const = 0;
	virtual void DrawSphere(const Vector3D& p1, double r, int tile=8, double fill=1) const = 0;
	virtual void DrawColorSphere(double value, const Vector3D& p, double r, int tile=8, double fill=1) = 0;
	virtual void DrawPolygon(const TArray<Vector3D>& p, int drawlines=0, double linewidth=1) const = 0;
	virtual void DrawPolygonOutline(const TArray<Vector3D>& p, double linewidth=1) const = 0;

	virtual void SetColor(const Vector3D& col) = 0;
	virtual void SetLineColor(const Vector3D& col) = 0;
	virtual const Vector3D& GetLineColor() = 0;
	virtual void SetTransparency(int transp_on) = 0;
	virtual int GetTransparency() const = 0;
	virtual void ChooseColor(float R, float G, float B) const = 0;
	virtual void SetDrawlines(int i) = 0;
	virtual void SetDrawlinesH(int i) = 0;

	// for geom objects
	virtual double GetObjData(int row, int col) const = 0;
	virtual double GetObjDataStepSize() const = 0;
	virtual Vector3D & ColLine() = 0;

	virtual int UseCuttingPlanes() = 0;
	virtual int CuttingPlanesAllow(const Vector3D& pos) = 0;

	//$ PG 2013-6-12: those functions were added in order to get some old models (in model.cpp) working again:
	virtual void SetCenterObject(int co, const Vector3D& offset) = 0;
	virtual int GetCenterObject() const = 0;
	virtual double GetMagnifyYZ() const = 0;
	virtual void SetMagnifyYZ(double amagyz) = 0;
	virtual void SetShowelem(const int __se) = 0;  // depreciated
	virtual void SetShowelem2(const int __se) = 0;  // depreciated
	virtual void SetShowelem3(const int __se) = 0;  // depreciated
};


//// AD: [ f.t.t.b. enums for 2D window here
	// identify the type and number of the element on the MBS-side
	typedef enum { TNoMBSElem = 0, TIOBlock = 1, TSensor = 2 } TMBSElementType;                          

	// identify the type and number of the element on the Draw-side
	typedef enum { TNoSubType = 0, 
		     TElementFrame = 1, TElementName = 2,
				 TInputNode = 3, TOutputNode = 4,
				 TConstructionNode = 5,  TConnectionLine = 6,
				 TTextObject = 7,
				 TSymbol = 8} 	TDrawElementType;   

	typedef enum { TNoGeoType = 0, TGLine = 1, TGRectangle = 2, TGEllipse = 3, TGText = 4} TGeomElementType;

	typedef enum {HLeft=0, HCenter=1, HRight=2,               // horizontal allingment
								VTop=0,  VCenter=3, VBottom=6,              // vertical allignment
							  ScaleDown2Fit = 128													// additional  properties : "scale to fit" scale down fontsize to fit in rectangle - use a binary value here
							} TTextAllign;

	//typedef enum {HLeft=1, HCenter=5, HRight=9,
	//              VTop=10, VCenter=50, VBottom=90} TTextAllign;

	//// AD: ]


struct MBS2DDrawingInterface
{
	virtual ControlWindowContext * GetCWC() = 0;

	virtual void AddDrawComponent_Line(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D p1, Vector2D p2, Vector3D col) = 0;
	virtual void AddDrawComponent_Rect(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background) = 0;
	virtual void AddDrawComponent_Ellipse(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_border, Vector3D col_background) = 0;
	virtual void AddDrawComponent_Text(int mbs_elnr, TMBSElementType mbs_type, int sub_elnr, TDrawElementType sub_type, Vector2D center, Vector2D size, Vector3D col_text, mystr& text, TTextAllign positioning) = 0;
};

struct MBSOptionsInterface
{
	virtual void SetIOption(int index, int data) = 0;
	virtual const int& GetIOption(int index) const = 0;
	virtual int& GetIOption(int index) = 0;
	virtual void SetDOption(int index, double data) = 0;
	virtual const double& GetDOption(int index) const = 0;
	virtual double& GetDOption(int index) = 0;
	virtual void SetTOption(int index, const char* data) = 0;
	virtual const char* GetTOption(int index) const = 0;

	virtual HOTINTOptions* GetOptions() = 0;
	virtual const HOTINTOptions* GetOptions() const = 0;

	// these two functions were required by the parser
	virtual ElementDataContainer* GetMBS_EDC_Options() = 0;
	virtual ElementDataContainer* GetMBS_EDC_Variables() = 0;
	// $EK 2013-01-11 added two members to interface
	virtual void SolverOptions2EDC(ElementDataContainer* edc) = 0;
	virtual void EDC2SolverOptions(const ElementDataContainer* edc) = 0;
};

// is often needed alone
struct MBSUOInterface
{
	virtual UserOutputInterface & UO(int message_level = UO_LVL_all, int output_prec = -1) const = 0;
};

struct MBSModelDataInterface
{
	virtual ElementDataContainer* GetModelDataContainer() = 0;
	// $EK 2013-01-11 added GetModelDataContainer_args
	virtual ElementDataContainer* GetModelDataContainer_args() = 0;
	virtual void SetModelDataContainer(const ElementDataContainer& edc) = 0;
	virtual int ReadModelData(mystr filename) = 0;
	virtual int AddModelData(mystr filename) = 0;
	virtual void AddReplaceModelDataEDC(ElementDataContainer& edc) = 0; // $AD: 2013-09-13 added for FOR command of script language ( no external file is read here )
};

// functionality related to parsing expressions
struct MBSParserInterface
{
	virtual double EvaluateParsedFunction1D(const mystr & parsedFunctionExpression, const mystr & parsedFunctionVariables, double t) = 0;
	// this function, required by Sim2Hotint, invokes expression parser to evaluate an expression
	virtual double ExpressionToDouble(mystr & expression) = 0;
	//this function reads a file and stores it in an edc - required by femesh and GeomMesh3D
	virtual int File2EDC(const char* filename, ElementDataContainer* edc_file) = 0;
	// this function evaluates a command of the parser, using the parameter_EDC as additional input data
	virtual int ExecuteParserCommand(mystr & command, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option) = 0;
	//this function computes mass, moment of inertia, volume and center of gravity based on the data about the geometry and the material
	virtual int ComputeInertia(ElementDataContainer* data, ElementDataContainer* return_value) = 0; //$ DR 2013-01-30
	//this function reads a column of a file and stores the vector in an edc
	virtual int LoadVectorFromFile(const char* filename, int col, ElementDataContainer* return_value) = 0;	//$ DR 2013-07-04 
	//this function returns a pointer to the parser (needed in MathFunction)
	virtual CMBSParser& MBSParser() = 0;
};

typedef enum
{
	TSimulationNotStarted        =  1,
	TSimulationRunning           =  TSimulationNotStarted << 1,
	TSimulationEndedRegularly    =  TSimulationNotStarted << 2,
	TSimulationStoppedDueError   =  TSimulationNotStarted << 3,
	TSimulationStoppedByUser     =  TSimulationNotStarted << 4,
	TSimulationStoppedByElement  =  TSimulationNotStarted << 5,
	TSimulationProcessFinished   =  TSimulationNotStarted << 6				// $ MaSch 2013-08-08: this flag is set in MultiBodySystem::PerformComputation() after the whole simulation process
																																		//(e.g. one single computation, or all computations of a parameter variation) has finished
} TSimulationStatus;

class SimulationStatus
{
protected:
	unsigned int status;

public:
	SimulationStatus(): status(unsigned int(0)) {};
	SimulationStatus(const TSimulationStatus & s): status(s) {};
	void SetStatusFlag(const TSimulationStatus & s) {status = s;}										//reset status to s
	void AddSetStatusFlag(const TSimulationStatus & s) {status |= s;}								//(additionally) set flag s
	void UnsetStatusFlag(const TSimulationStatus & s) {status &= (~s);}							//delete flag s
	unsigned int GetStatusFlag(const TSimulationStatus & s) {return(status & s);}		//returns 0 if flag s is not set, otherwise s (!=0) 
	void ResetStatusFlags() {status &= unsigned int(0);}
};


struct MBS :
	public MBSSolutionAccessInterface,
	public MBSObjectsAccessInterface,
	public MBS3DDrawingInterface,
	public MBS2DDrawingInterface,                 
	public MBSOptionsInterface,
	public MBSUOInterface,
	public MBSModelDataInterface,
	public MBSParserInterface
{
	// specific functions
	virtual double GetTime() const = 0;
	virtual TimeIntLog GetTimeIntLog() const = 0;
	virtual double GetDrawTime() const = 0;
	virtual int GetCSNumber(double globaltime) = 0;
  virtual double GetCSTime(double globaltime) = 0;
	virtual double LoadFact() = 0;

	virtual void Assemble() = 0;
	virtual void BuildLTGLists() = 0;
	virtual int MaxIndex() const = 0;
	virtual void SetMaxIndex(int i) = 0;
	virtual int UseDependencies() const = 0;
	virtual void ClearSystem() = 0;

	virtual const double& GetVelocityAndAcceleration(int i) const = 0; // $ MSax 2013-07-16 : changed from GetAcceleration to GetVelocityAndAcceleration

	virtual const SolverSettings& GetSolSet() const = 0;
	virtual SolverSettings& GetSolSet() = 0;

	// parameter variation access functions
	virtual const double& GetVarparameter() const = 0;
	virtual double& GetVarparameter() = 0;
	virtual const double& GetVarparameter2() const = 0;
	virtual double& GetVarparameter2() = 0;

	virtual void StopByElement() = 0;
	virtual SimulationStatus & GetSimulationStatus() = 0;
	//virtual const NumNLSolver& NumSolver() const = 0;
	virtual NumSolverInterface& NumSolver() = 0;
	virtual void EDCError(const mystr& str) = 0;
	virtual double GetMagnifyYZ() const = 0;
	virtual FieldVariableDescriptor * GetActualPostProcessingFieldVariable() = 0;
	virtual int IsLoadSaveMode() const = 0;
	virtual int IsJacobianComputation() const = 0;
	virtual double DiscontinuousAccuracy() const = 0;
	virtual double& DiscontinuousAccuracy() = 0;
	virtual int StopCalculation() const = 0;
	virtual int GetSystemSize() const = 0;	// for BaseCMSElement
	virtual void SetTransformJacApply(int i) = 0;	// for CMSElement2D
	virtual int UseSparseSolver() const = 0;		// ACRSElement2D
	virtual double GetStepSize() const = 0;			// ACRSElement2D
	virtual int GetSecondOrderSize() const = 0;		// specialconstraints

	virtual void SetUseDependencies(int flag) = 0;

	virtual void SetReducedBandsize(int reducedbandsizeI) = 0;

	virtual double CharacteristicLength() const = 0;

	// functions below appeared due to ancf...2d
	virtual void UpdateFEMinMaxCol(double val) = 0;
	virtual const double& GetTImincol() const = 0;
	virtual double& GetFEmincol() = 0;
	virtual const double& GetTImaxcol() const = 0;
	virtual double& GetFEmaxcol() = 0;
	virtual Vector3D FEColor(double val) const = 0;
	virtual int GetDrawResolution() const = 0;

	// contact2D/3D, control
	virtual int DoStaticComputation() const = 0;
	virtual double GetStepRecommendation() const = 0;
	virtual void SetStepRecommendation(double s) = 0;
	virtual double GetStepSizeNew() = 0;
	virtual double GetStepEndTime() const = 0;
	virtual void ForceJacobianRecomputation() = 0;

	// this function returns the current hotint version, as defined in "WorkingModule\hotint_version.h"
	virtual const HotintVersionInfo& GetHotintVersion() const = 0;

	//$ MaSch 2013-08-26
	//TCP/IP
	virtual TCPIPHotInt* & GetServerSocket() = 0;

};

//Use the following function in order to register Model functions automatically in the MBS system
void ModelFunctionAutoRegistration(int(*function_ptr)(MBS* mbs), const char* functionName, const char* description, int option,
																	 int(*function_ptr_init_modeldata)(MBS* mbs) = 0);