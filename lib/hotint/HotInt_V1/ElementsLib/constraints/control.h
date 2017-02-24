//#**************************************************************
//#
//# filename:             control.h           
//#
//# author:               Gerstmayr Johannes, Rafael Ludwig
//#
//# generated:						9 October 2008
//# description:          Elements for control, transmission elements, etc.
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
 
#ifndef CONTROL__H
#define CONTROL__H

#include "TcpIpRoutines.h"

const int IOInputTypeElement = 1;
const int IOInputTypeSensor = 2;

//$ RE 2012-11-15: autogeneration functions and variables added

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: InputOutputElement (Constraint)
// short description: Superclass for continuous input-output elements
// available formulations: Continuous state-space elements
// type: SISO (single input-single output); MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description:	Superclass providing basic functions for CONTINUOUS input-output elements ==> only subclasses are useful for simulation (e.g. for controller circuits)
// class variables:
//   -inputs: sensor number i: if (input_type==1==IOInputTypeElement) --> global number of InputOutputElement; if (input_type==2==IOInputTypeSensor) --> global sensor number
//   -input_types: 1==InputOutputElement, 2==Sensor
//   -input_localnum: output number of sensor: if (input_type==1) --> output number of InputOutputElement
//   -n_output: number of ouputs
//   -n_state: number of state variables (integrators)
//   -ref_pos: reference position of InputOutputSymbol
//   -draw_dim: drawing parameters
//   -rotation: rotation, 1==90°, 2==180°, 3==270°, 4=360°
//   -colbackground: background color of control objects
//   -colforeground: foreground color of control objects
//   -input_nodes: nodal positions of connections of inputs
//   -input_nodes_num: input number for according input_nodes
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class InputOutputElement: public Constraint  //$EDC$[beginclass,classname=InputOutputElement,parentclassname=Constraint,texdescription="
//The InputOutputElement is the superclass for continuous input-output elements (e.g.: gain, transfer functions, etc.). 
//Specialized subclasses of the InputOutputElement can be used in the model. Subclasses of these port-block elements can be connected by 
//their inputs and outputs. The values of the outputs can be measured by sensors.
//"]
{
public:

	InputOutputElement(MBS* mbsi):Constraint(mbsi), inputs(0), input_types(0), input_localnum(0)
	{	
		mbs = mbsi; //is set in ::Contraint::Element(mbsi)
		InitConstructor();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new InputOutputElement(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor();

	virtual const char* GetElementSpec() const {return "InputOutputElement";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
	virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables);
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata);

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual void EvalF(Vector& f, double t) 
	{
		//to be overwritten in specific class!

	};
	virtual void EvalG(Vector& f, double t) {};

	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f) {};

	virtual void ComputationFinished() {}; //function is called when computation is finished (e.g. in order to free memory, write results, close connections, etc.)

	virtual int IS() const {return 0;};
	virtual int ES() const {return n_state;};
	virtual int Dim() const {return 2;}  //drawing is 2D

	//get a reference position of the body in 3d
	virtual Vector2D GetRefPos2DD() const 
	{
		return Vector2D(ref_pos.X(), ref_pos.Y());
	}

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const 
	{
		return Vector3D(GetRefPos2DD().X(), GetRefPos2DD().Y(), 0.);
	}

	virtual void SetDrawBackgroundColor(const Vector3D& coli){colbackground = coli;}
	virtual const Vector3D& GetDrawBackgroundColor(){return colbackground;}

	virtual void SetDrawForegroundColor(const Vector3D& coli){colforeground = coli;}
	virtual const Vector3D& GetDrawForegroundColor(){return colforeground;}

	virtual void AddInput(int ioelem_num, int input_type, int localnumI=1)
	{
		inputs.Add(ioelem_num);
		input_types.Add(input_type);
		input_localnum.Add(localnumI);

		if (input_type == IOInputTypeElement) //input is other input-output element
		{
			AddElement(ioelem_num); //add Dependency to this element!
		}
		else if (input_type == IOInputTypeSensor) //input is sensor
		{
			AddSensor(ioelem_num);
		}
	}

	// attention: all input_portnum must have correct connection to input element
	virtual void AddInput(int ioelem_num, int input_type, int localnumI, int input_portnum)
	{
		inputs(input_portnum) = ioelem_num;
		input_types(input_portnum) = input_type;
		input_localnum(input_portnum) = localnumI;

		if (input_type == IOInputTypeElement) //input is other input-output element
		{
			AddElement(ioelem_num); //add Dependency to this element!
		}
		else if (input_type == IOInputTypeSensor) //input is sensor
		{
			AddSensor(ioelem_num);
		}
	}

	virtual int GetInputNum(int i) const {return inputs(i);}
	virtual int GetInputLocalNum(int i) const {return input_localnum(i);}
	virtual int GetInputType(int i) const {return input_types(i);}
	virtual int GetNInputs() const {return inputs.Length();} //$EDC$[funcaccess,readonly,EDCvarname="number_of_inputs",EDCfolder="IOBlock",tooltiptext="number of inputs"]
	virtual int GetNOutputs() const {return n_output;}       //$EDC$[funcaccess,readonly,EDCvarname="number_of_outputs",EDCfolder="IOBlock",tooltiptext="number of outputs"]
	virtual int GetNStates() const {return n_state;}         //$EDC$[funcaccess,readonly,EDCvarname="number_of_states",EDCfolder="IOBlock",tooltiptext="number of states"]

	virtual void SetNOutputs(int i) {n_output = i;}

	virtual void SetOutputName(int output_nr, mystr& name);
	virtual char* GetOutputName(int output_nr);

	virtual void SetNStates(int i) {n_state = i;}

	virtual void SetRefPos2D(Vector2D rp) {ref_pos = rp;}
	virtual void SetDrawDim(const Vector3D& drawdim) {draw_dim = drawdim;}
	virtual void SetRotation(double rot) {rotation = rot;}; //1=90°, 2=180°, ...
	virtual double GetRotation() const {return rotation;}

	virtual double GetInput(double t, int i=1) const
	{
		if(input_types.Length() < i || inputs.Length() < i)
		{
			mbs->UO(UO_LVL_err).InstantMessageText(mystr("Error: input ") + mystr(i) + mystr(" of element with number ") + mystr(GetOwnNum()) + mystr(" not found. Input is set to zero!!!\n"));
			return 0.;
		}
		//if(inputs(i) < 1 || inputs(i) > mbs->NE())
		if(inputs(i) < 1 || (inputs(i) > mbs->NE() && input_types(i) == IOInputTypeElement) || (inputs(i) > mbs->NSensors() && input_types(i) == IOInputTypeSensor)) // $ MSax 2013-07-29: bugfix
		{	
			mbs->UO(UO_LVL_err).InstantMessageText(mystr("Error: element with number ") + mystr(inputs(i)) + mystr(" not found --> input of element with number ") + mystr(GetOwnNum()) + mystr(" is set to zero!!!\n"));
			return 0.;
		}

		if (input_types(i) == IOInputTypeElement) //InputOutputElement
		{
			const InputOutputElement& ioe = (const InputOutputElement&)GetMBS()->GetElement(inputs(i));
			//const InputOutputElement* ioe = (InputOutputElement*)GetMBS()->GetElementPtr(inputs(i));
			return ioe.GetOutput(t, input_localnum(i));
		}
		else if (input_types(i) == IOInputTypeSensor)
		{
			//$ YV 2012-06: the sensors may produce just one scalar value
			const Sensor & s = GetMBS()->GetSensor(inputs(i));
			//$ YV 2012-06: here we use an explicit conversion to avoid difficulties with the ``const'' modifier
			return ((Sensor &)s).GetCurrentValueWithSensorProcessing(t);
		}
		else 
		{
			//should not happen:
			mbs->UO().InstantMessageText("Unknown input_type --> input of element " + elementname + " was set to zero!");
			return 0;
		}
	}

	virtual double GetOutput(double t, int i=1) const {return 0;}  //to be overwritten in specific class!		
	

	virtual void GetDirectFeedThroughElements(TArray<int>& elnums) const; //return all elements which are depending on this element by direct feed through
	virtual int IsDirectFeedThrough() const {return 0;}

	virtual void SetInitialValues(const Vector& xi) {x_init = xi;}

	virtual Box3D GetElementBox() const
	{
		if (GetMBS()->GetIOption(144))
		{
			// orientation
			double phi = rotation*MY_PI*0.5;
			Matrix3D rot = RotMatrix3(phi);

			Vector2D rp = GetRefPos2D();		// reference position
			Vector3D rp3 = ToP3D(rp);				// reference position 3D

			// rectangle
			double b = draw_dim.X();
			double h = draw_dim.Y();
			Vector3D p1(-0.75*b,-0.75*h,0.);
			Vector3D p2( 0.75*b,-0.75*h,0.);
			Vector3D p3( 0.75*b, 0.75*h,0.);
			Vector3D p4(-0.75*b, 0.75*h,0.);

			// rotate and translate rectangle
			p1 = rp3+rot*p1;
			p2 = rp3+rot*p2;	
			p3 = rp3+rot*p3;
			p4 = rp3+rot*p4;

			// box
			double dimZ = sqrt(b*h);
			Box3D box;
			box.Add(Vector3D( p1.X(), p1.Y(),-dimZ));
			box.Add(Vector3D( p2.X(), p2.Y(),-dimZ));
			box.Add(Vector3D( p3.X(), p3.Y(),-dimZ));
			box.Add(Vector3D( p4.X(), p4.Y(),-dimZ));
			box.Add(Vector3D( p1.X(), p1.Y(), dimZ));
			box.Add(Vector3D( p2.X(), p2.Y(), dimZ));
			box.Add(Vector3D( p3.X(), p3.Y(), dimZ));
			box.Add(Vector3D( p4.X(), p4.Y(), dimZ));
			box.Increase(dimZ);
			return box;	
		}
		else
		{
			Box3D b;
			return b;
		}
	}

	virtual Box3D GetElementBoxD() const
	{
		if (GetMBS()->GetIOption(144))
		{
			// orientation
			double phi = rotation*MY_PI*0.5;
			Matrix3D rot = RotMatrix3(phi);

			Vector2D rp = GetRefPos2DD();		// reference position
			Vector3D rp3 = ToP3D(rp);				// reference position 3D

			// rectangle
			double b = draw_dim.X();
			double h = draw_dim.Y();
			Vector3D p1(-0.75*b,-0.75*h,0.);
			Vector3D p2( 0.75*b,-0.75*h,0.);
			Vector3D p3( 0.75*b, 0.75*h,0.);
			Vector3D p4(-0.75*b, 0.75*h,0.);

			// rotate and translate rectangle
			p1 = rp3+rot*p1;
			p2 = rp3+rot*p2;	
			p3 = rp3+rot*p3;
			p4 = rp3+rot*p4;

			// box
			double dimZ = sqrt(b*h);
			Box3D box;
			box.Add(Vector3D( p1.X(), p1.Y(),-dimZ));
			box.Add(Vector3D( p2.X(), p2.Y(),-dimZ));
			box.Add(Vector3D( p3.X(), p3.Y(),-dimZ));
			box.Add(Vector3D( p4.X(), p4.Y(),-dimZ));
			box.Add(Vector3D( p1.X(), p1.Y(), dimZ));
			box.Add(Vector3D( p2.X(), p2.Y(), dimZ));
			box.Add(Vector3D( p3.X(), p3.Y(), dimZ));
			box.Add(Vector3D( p4.X(), p4.Y(), dimZ));
			box.Increase(dimZ);
			return box;
		}
		else
		{
			Box3D b;
			return b;
		}
	}

	virtual void DrawElement();
	virtual void DrawElement2D();
	virtual void DrawBlockSymbol();

	virtual const char* SymbolText() const;

	virtual void AddInputNode(int input_num, const Vector2D& node) // Adds a construction Node to an input ( at "outer" position )
	{
		input_nodes.Add(node);
		input_nodes_num.Add(input_num);
	}

	virtual Vector2D ConNodePos(int i) { return input_nodes(i); }
// list operation: insert a construction node into the list ( both lists, at position i )
	virtual void InsertConNode(int list_idx, int input_nr, Vector2D& node) 
	{
// insert at array position defined by list_idx
		input_nodes.Insert(list_idx, node);
		input_nodes_num.Insert(list_idx, input_nr); 
	}
// list operation: delete the construction node from both arays
	virtual void DeleteConNode(int list_idx)  // 
	{
		input_nodes.Erase(list_idx);
		input_nodes_num.Erase(list_idx);
	}
	virtual void MoveConNode2D(int list_idx, double delta_x, double delta_y)
	{
		Vector2D& p = input_nodes.Elem(list_idx);
		p.X() += delta_x;
		p.Y() += delta_y;

	}
	virtual void MoveElement(double delta_x, double delta_y, double delta_z)
	{
		//Vector2D& refpos = ref_pos // IOElement has only 2D refpos
		ref_pos.X() += delta_x;
		ref_pos.Y() += delta_y;
//		refpos.Z() += delta_z;

	}


	virtual Vector2D GetInputPosD(int i) const; //return absolute position of input #i
	virtual Vector3D GetInputPos3DD(int i) const {assert(0 && "Only used for 3D-Elements!"); return Vector3D(0.);}; //return absolute position of input #i

	virtual Vector2D GetOutputPosD(int i) const; //return absolute position of input #i
	virtual Vector3D GetOutputPos3DD(int i) const {assert(0 && "Only used for 3D-Elements!"); return Vector3D(0.);}; //return absolute position of input #i

	virtual Vector3D ToP3D(const Vector2D& p) const 
	{
		return Vector3D(p.X(),p.Y(),0.);
	}

	// functions to catch the input - implemented in derived class
	virtual int RespondToKey(int key) {return false;}

protected:

	//system inputs	
	TArray<int> inputs;			//$EDC$[varaccess,EDCvarname="input_element_numbers",EDCfolder="IOBlock",variable_length_vector,tooltiptext="vector of element(s) or sensor number(s) connected to input, only valid element numbers permitted!"] // sensor number i: if (input_type==1==IOInputTypeElement) --> global number of IOutputElement
	//if (input_type==2==IOInputTypeSensor) --> global sensor number
	TArray<int> input_types;	//$EDC$[varaccess,EDCvarname="input_element_types",EDCfolder="IOBlock",variable_length_vector,tooltiptext="vector with types of connected inputs; 1=IOElement, 2=Sensor"] // 1==Input-OutputElement, 2==Sensor
	TArray<int> input_localnum;	//$EDC$[varaccess,EDCvarname="input_local_number",EDCfolder="IOBlock",variable_length_vector,tooltiptext="vector with i-th number of output of previous IOelement connected to this element"] // output number of sensor: if (input_type==1) --> output number of IOutputElement

	int n_output;    //number of system outputs, use dedicated function [S|G]etNOutputs()
	int n_state;     // number of state variables, use dedicated function GetNInputs()

	//drawing properties:
	Vector2D ref_pos;       //$EDC$[varaccess,EDCvarname="position",EDCfolder="Graphics",tooltiptext="reference drawing position"]
	Vector3D draw_dim;      //$EDC$[varaccess,EDCvarname="draw_size",EDCfolder="Graphics",tooltiptext="draw size"]
	double rotation;	    //$EDC$[varaccess,EDCvarname="rotation",EDCfolder="Graphics",tooltiptext="rotation: 1==90°, 2==180°, 3==270°, 4=360°"]
	Vector3D colbackground; //$EDC$[varaccess,EDCvarname="background_color",EDCfolder="Graphics",tooltiptext="background color; -1=transparent"]
	Vector3D colforeground; //$EDC$[varaccess,EDCvarname="foreground_color",EDCfolder="Graphics",tooltiptext="foreground color"]

	// further drawing information:
	TArray<Vector2D> input_nodes; //!EDC$[varaccess,readonly_varaccess,EDCfolder="Graphics",EDCvarname="input_nodes",variable_length_vector,tooltiptext="position of drawing nodes of connection line to corresponding to the input with number \"input_node_number\""] //nodal positions of connections of inputs
	TArray<int> input_nodes_num;  //$EDC$[varaccess,EDCfolder="Graphics",EDCvarname="input_nodes_num",variable_length_vector,tooltiptext="number of input of drawing position \"input_nodes\""]

	// list containing names of outputs
	MyStrList output_names;

	//// remove entries from EDC
	//EDC	int use_penalty_formulation;		//$EDC$[varaccess,remove,EDCvarname="use_penalty_formulation",EDCfolder="Physics"]
	//EDC	double spring_stiffness;   		    //$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC	int use_local_coordinate_system;	//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC   Vector3D col;                       //$EDC$[varaccess,remove,EDCvarname="RGB_color",EDCfolder="Graphics"]

};//$EDC$[endclass,InputOutputElement]

const double discrete_time_tol = 1e-10;
const double discrete_time_tol_inv = 1e10;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: InputOutputElementDiscrete (InputOutputElement)
// short description: Superclass for discontinuous input-output elements
// available formulations: Discontinuous state-space elements
// type: SISO (single input-single output); MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description:	Superclass providing basic functions for DISCONTINUOUS input-output elements ==> only subclasses are useful in simulation e.g. for time discrete controller circuits
// discrete system with constant sample time
// yk = f(u, y, k)
// u = [u_(k-n), u_(k-(n-1)), ..., u_(k-n+m) ]
// y = [y_(k-n), y_(k-(n-1)), ..., y_(k-(n-1)), y_(k-1)]
// m < n
// class variables:
//   -deltaT: discrete time step
//   -toff: offset to sample time (Tk = k*dT + off)
//   -n_disc_states: number of time discrete variables
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class InputOutputElementDiscrete: public InputOutputElement//$EDC$[beginclass,classname=InputOutputElementDiscrete,parentclassname=InputOutputElement]
{
public:

	InputOutputElementDiscrete(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetSampleTime(double val){		deltaT = val;	};
	virtual double GetSampleTime(){return	deltaT;	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new InputOutputElementDiscrete(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const InputOutputElementDiscrete& ce = (const InputOutputElementDiscrete&)e;

		//k = ce.k;
		deltaT = ce.deltaT;
		toff = ce.toff;
		n_disc_states = ce.n_disc_states;
	}

	virtual int GetNDiscreteStates() const {return n_disc_states;}
	virtual void SetNDiscreteStates(int n_disc_statesI) {n_disc_states = n_disc_statesI;}

	virtual void InitConstructor()
	{
		SetNStates(0);
		SetNOutputs(1);
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		toff = 0;
		deltaT = 0;
		//k = 0;

		SetNDiscreteStates(0);
		Vector datainit(DataS());
		datainit.SetAll(0.);
		SetDataInit(datainit);
	}
	//###########################
	// DISCRETE ELEMENT FUNCTIONS
	virtual int DataS() const {return 1+GetNDiscreteStates();} // size of data vector, which will is automatically restored if solver-step not valid

	virtual int GetK() const 
	{
		return (int)XData(1);
	}
	virtual void SetK(int k) {XData(1) = (double)k;}

	virtual const double& XDiscrete(int i) const {return XData(1+i);}
	virtual double& XDiscrete(int i) {return XData(1+i);}

	virtual double Roundval(double x) const;

	virtual double GetDiscreteTime() const;
	virtual double GetNextDiscreteTime() const;

	virtual void StartTimeStep(); // update index "K" here and make shorter step size if it is too big
	virtual void EndTimeStep();

	virtual int isDiscreteEvent(double discrete_time) const;

	// to be overwritten in children elements
	virtual void UpdateDiscreteElementState()	{assert(0);}; //*** xk+1 = f(xk,uk), define function in child class
	virtual void UpdateDiscreteElementInput(double t) {assert(0);}; //*** read actual uk , define function in child class
	virtual double GetOutput(double t, int i=1) const 	{	assert(0); return 0; }  //*** define function in child class
	//###########################

	virtual const char* GetElementSpec() const {return "IODiscrete";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual int IsDirectFeedThrough() const {	assert(0); return 0;	}   //*** define function in child class
	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();

protected:
	double deltaT, toff; 
	int n_disc_states; // number of discrete variables (without k)
};//$EDC$[endclass,InputOutputElementDiscrete]

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: ZTransferFunction (InputOutputElementDiscrete)
// short description: Discontinuous Transfer function
// available formulations: Discontinuous state-space elements
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description:	 Discontinuous Transfer function in z-Space
// Realization of Z-Transfer Function (Z-transform):
// theorical background:
// y(z) = G(z)*u(z); G(z) = (num(1)+num(2)*z + ... + num(n)*z^(n+1)) / (den(1) + den(2)*z + ... + den(n)*z^(n+1)); 
// Z-Space
// y_k * z^n  o---------0 y_(k+n) 
// yk = f(u, y, k)
// u = [u_(k-n), u_(k-(n-1)), ..., u_(k-n+m) ]^t
// y = [y_(k-n), y_(k-(n-1)), ..., y_(k-(n-1)), y_(k-1)]^t
// class variables:
//   -num: numerator of transfer function	
//   -den: denominator of transfer function	
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ZTransferFunction: public InputOutputElementDiscrete//$EDC$[beginclass,classname=ZTransferFunction,parentclassname=InputOutputElementDiscrete,addelementtype=TAEinput_output,addelementtypename=IODiscreteTransferFunction,
//texdescription="Discontinuous transfer function in z-space. It is a SISO (single input-single output) control element. Inital state is zero.",
//figure="IODiscreteTransferFunction,IODiscreteTransferFunction",
//texdescriptionEquations="
//$y(z)=\mathbf{G}(z)u(z)$ \\ \\
//$\mathbf{G}(z)=\frac{\mathbf{num}}{\mathbf{den}}$ \\ \\
//user input: \\ 
//$\mathbf{num}(z) = num_1 + num_2z + num_3z^2 +...+ num_{n+1}z^n$ \\ 
//$\mathbf{den}(z) = den_1 + den_2z + den_3z^2 +...+ den_{n+1}z^n$ \\ 
//Theoretical background: Realization of z-transfer function as time discrete state space model \\ 
//\begin{eqnarray}
//	\left[\begin{array}{c}
//		z_{k+1,1} \medskip \\
//		. \medskip \\
//		. \medskip \\
//		. \medskip \\
//		z_{k+1,n}
//	\end{array} \right]
//	=\left[
//	\begin{array}{ccccc}
//		0 & . & . & 0 & -den_1 \medskip \\
//		1 & 0 & . & 0 & -den_2 \medskip \\
//		. & . & . & . & . \medskip \\
//		0 & 0 & . & 1 & -den_n
//	\end{array} \right] .
//	\left[\begin{array}{c}
//		z_{k,1} \medskip \\
//		. \medskip \\
//		. \medskip \\
//		. \medskip \\
//		z_{k,n}
//	\end{array} \right] +
//	\left[\begin{array}{c}
//		num_1 - num_{n+1}den_1 \medskip \\
//		. \medskip \\
//		. \medskip \\
//		num_n - num_{n+1}den_n
//	\end{array} \right] u_k
//\end{eqnarray}
//\begin{equation}
//  y_k = z_{k,n}+num_{n+1}u_z;
//\end{equation}",
//example="ZTransferFunction.txt"]



//\begin{eqnarray}
//	\mathbf{B} &=&\left[
//	\begin{array}{c}
//		num(1) - num(n+1)*den(1) \medskip \\
//		. \medskip \\
//		. \medskip \\
//		num(n) - num(n+1)*den(n)
//	\end{array}
//	\right] \quad
//	\end{array}

//$y(z)=\mathbf{G}(z)u(z); \mathbf{G}(z) = \frac{num(1)+num(2)z+...+num(n)z^{n+1}}{den(1)+den(2)+...+den(n)z^{n+1}}$ \\ \\
//z-space \\ 
//$y_k*z^n$ $\circ---\circ$ $y_{(k+n)}$ \\ 
//$y_k = f\left(u,y,k\right)$ \\ 
//$u = \left[u_{(k-n)}, u_{(k-(n-1))},...,u_{(k-n+m)}\right]^T$ \\ 
//$y = \left[y_{(k-n)}, y_{(k-(n-1))},...,y_{(k-(n-1))},y_{(k-1)}\right]^T$
{
public:

	ZTransferFunction(MBS* mbsi):InputOutputElementDiscrete(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetZTransferFunction(const Vector& numI, const Vector& denI, const Vector& initVector)
	{
		//-------------------------------------------------------------------
		// G(z) = p_num(z)/p_den(z)
		// p_num(z) = num(1) + num(2) * z + num(3) * z^2 + ... + num(n) * z^m
		// p_den(z) = den(1) + den(2) * z + den(3) * z^2 + ... + den(n) * z^n
		// m<=n
		//-------------------------------------------------------------------
		//check if sizes are consistent!
		assert(num.Length() == den.Length()); 

		//inputs must be added separately	
		num = numI;
		den = denI;
		double fact = den(den.Length());
		assert(fact != 0);

		num *= (1./fact);
		den *= (1./fact);

		SetNOutputs(1);
		SetNStates(0);                	
		int dimIODiscrete = DataS(); // dimension of root element
		SetNDiscreteStates(initVector.Length());

		Vector init;
		init.SetLen(0); // new init vector inclusive base element variables
		init = init.Append(GetDataInit());
		init = init.Append(initVector); 
		init = init.Append(0.); // u0=0.0

		SetDataInit(init);		

		//UO() << "init=" << init << "\n";
	}

	virtual void SetZTransferFunctionWithDescentPolynoms(const Vector& numI, const Vector& denI)
	{
		//-------------------------------------------------------------------
		// G(z) = p_num(z)/p_den(z)
		// p_num(z) = num(m) + num(m-1) * z + num(m-2) * z^2 + ... + num(1) * z^m
		// p_den(z) = den(n) + den(n-1) * z + den(n-2) * z^2 + ... + den(1) * z^n
		// m<=n
		//-------------------------------------------------------------------
		//check if sizes are consistent!
		assert(numI.Length() == denI.Length());

		// revert numerator and denominator
		int len = numI.Length();
		Vector tmp1(len), tmp2(len);
		for(int i = 1; i <= len;i++)
		{				
			tmp1(len - i + 1) = numI(i);
			tmp2(len - i + 1) = denI(i);
		}
		SetZTransferFunction(tmp1, tmp2);
	}

	virtual void SetZTransferFunction(const Vector& numI, const Vector& denI)
	{
		//check if sizes are consistent!
		assert(numI.Length() == denI.Length());

		Vector initVector(denI.Length()-1);
		initVector.SetAll(0);

		SetZTransferFunction(numI, denI, initVector);

	}

	virtual int DataS() const {return 2+GetNDiscreteStates();}    // dim(k)=1, dim(xk)=n, dim(uk) = 1 --> 2+n
	virtual void SetUk(double valI)	{		XData(DataS()) = valI;	} //set input u
	virtual double GetUk() const {		return XData(DataS());	}//get input u

	// set uk
	virtual void UpdateDiscreteElementInput(double t)
	{
		SetUk(GetInput(t,1));
	};

	//function is called every discrete event
	// xk+1 = f(xk,uk)
	virtual void UpdateDiscreteElementState()
	{
		//x(k+1) = A*x(k)+b*u(k);
		//A=[0 ..... 0 -den1    ]
		//  [1 0 ... 0 -den2    ]
		//  [0 0 ... 1 -den(n)]
		//bb(i) = num(i) - num(n+1)*den(i)
		int n = GetNDiscreteStates();
		double u = GetUk();

		xk1_temp.SetLen(n);
		if(n>0)	{	xk1_temp(1) = -den(1)*XDiscrete(n) + (num(1)-num(n+1)*den(1))*u;}
		for(int i = 2; i<=n; i++ )
		{
			xk1_temp(i) = - den(i)*XDiscrete(n) + (num(i)-num(n+1)*den(i))*u ; //A(i,i-1) * x(i-1) + b(i)*u;
			xk1_temp(i) += XDiscrete(i-1);   // ones in A-matrix
		}

		for(int i = 1; i<=n; i++ )
		{
			XDiscrete(i) = xk1_temp(i);
		}
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ZTransferFunction(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElementDiscrete::CopyFrom(e);
		const ZTransferFunction& ce = (const ZTransferFunction&)e;

		num = ce.num;
		den = ce.den;
		//init_vec = ce.init_vec; // $ MSax 2013-03-01: added
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();		

		num = Vector(1); num(1) = 1;// $MSax 2013-02-28: added
		den = Vector(1); den(1) = 1;
		//init_vec = Vector(1);
	}

	virtual const char* GetElementSpec() const {return "IODiscreteTransferFunction";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const 
	{
		//to be overwritten in specific class!
		if (i != 1) 
		{
			assert(0);
			return 0;
		}
		int n = GetNDiscreteStates();
		double u = 0;
		if (num(n+1) != 0.)
		{
			u = GetUk();
		}

		double y;
		if (n==0)
		{
			y = num(n+1)*u;
		}
		else 
		{
			y = XDiscrete(n);
			y += num(n+1)*u;
		}
		return y;
	}

	virtual int IsDirectFeedThrough() const 
	{
		//assert(num.Length() == den.Length());
		return num(num.Length()) != 0;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();

protected:
	Vector num; //!EDC$[varaccess,EDCvarname="num",EDCfolder="IOBlock",tooltiptext="Coefficients of numerator polynomial of z-function, b0+b1*z+b2*z^2+...+bn*z^n"]
	Vector den; //!EDC$[varaccess,EDCvarname="den",EDCfolder="IOBlock",tooltiptext="Coefficients of denominator polynomial of z-function, a0+a1*z+a2*z^2+...+am*z^m"] //numerator and denominator of transfer function	
	Vector xk1_temp; //temporary variable
	//Vector init_vec; // $ MSax 2013-03-01: added
};//$EDC$[endclass,ZTransferFunction]


//--------------------------------------------------------------------------------------------------------------
// Realization of (pseudo)random value generator

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: Random Source (InputOutputElementDiscrete)
// short description: Discontinuous random source
// available formulations: Discontinuous state-space elements; internal random generator and Linear Feedback Shift Register (LFSR)
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Discontinuous random source using alternatively an internal C++ based pseudo random generator or a Linear Feedback Shift Register
// class variables:
//   -amplitude, offset, seed: maximal amplitude, additional offset of random signal, seed € [0.,1.]... initialization of random generator
//	 -TRandomType method: method of random value generation; method = Trand --> use built-in command "rand"; method = TLFSR --> Linear Feedback Shift Register
//   -rand_max: maximal value for random value generator RAND_MAX, alternative: 2^N-1
//   -isConstAmplitude: set to 1 if output value should be only +amplitude or -amplitude
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef enum {Trand = 0, TLFSR = 1} TRandomType; // Trand ... use built-in command "rand", TLFSR ... Linear Feedback Shift Register

class RandomSource: public InputOutputElementDiscrete//$EDC$[beginclass,classname=RandomSource,parentclassname=InputOutputElementDiscrete,addelementtype=TAEinput_output,addelementtypename=IORandomSource,
//texdescription="Discontinuous random source using alternatively an internal C++ based pseudo random generator or a linear feedback shift register. It has no input and one output.",
//figure="IORandomSource,IORandomSource",example="addRandomSource.txt",
//modus="{method $0$}{IOBlock.method must be set to $0$. The built-in random generator is used.}",
//modus="{method $1$}{IOBlock.method must be set to $1$. Generate a pseudo random binary signal by using Linear Feedback Shift Register.}"]

{
public:

	// constructor
	RandomSource(MBS* mbsi):InputOutputElementDiscrete(mbsi)
	{	
		InitConstructor();
	};

	// set function with number of bits of random signal
	virtual void SetRandomSource(double amplitudeI, double offsetI = 0., TRandomType methodI = (TRandomType)Trand, double seedI = 0, double init_val = 0, int Nbits = 15);
	// InitRandomSource is used in SetRandomSource
	virtual void InitRandomSource(double amplitudeI, double offsetI = 0., TRandomType methodI = (TRandomType)Trand, double seedI = 0, double init_val = 0, int rand_maxi = RAND_MAX);

	// random integer number € [0, RAND_MAX]; method TLFSR returns numbers unequal zero
	virtual int GetRandomInt();

	void SetRandMax(int val){rand_max = val;} // set max. value

	int GetRandMax(){return rand_max;	}       // get max. value

	void SetConstantAmplitude(int val = 1){isConstAmplitude = val;}
	// size of data vector
	virtual int DataS() const {return 1+GetNDiscreteStates();}    // dim(k)=1, dim(xk)=1 --> 2

	// set uk
	virtual void UpdateDiscreteElementInput(double t)	{	}; // no input

	//function is called every discrete event
	// xk+1 = f(xk,uk)
	virtual void UpdateDiscreteElementState();

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new RandomSource(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElementDiscrete::CopyFrom(e);
		const RandomSource& ce = (const RandomSource&)e;

		amplitude = ce.amplitude;
		offset = ce.offset;
		seed = 0;
		method = ce.method;
		rand_max = ce.rand_max;
		isConstAmplitude = ce.isConstAmplitude;
		init_val = ce.init_val;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();		
		amplitude = 0.;
		offset = 0.;
		seed = 0;
		method = Trand;
		rand_max = 32767; //RAND_MAX<=>0x7fff
		isConstAmplitude = 0;
		init_val = 0;
	}

	virtual const char* GetElementSpec() const {return "IORandomSource";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 0;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const 
	{
		//to be overwritten in specific class!
		if(isConstAmplitude)
		{
			return amplitude*Sgn(XDiscrete(1));
		}
		return XDiscrete(1); //yk = xk
	}

	virtual int IsDirectFeedThrough() const 
	{
		return 0;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();

protected:

	double amplitude; //!EDC$[varaccess,EDCvarname="amplitude",EDCfolder="IOBlock",tooltiptext="maximal amplitude of random signal"] 
	double offset; //!EDC$[varaccess,EDCvarname="offset",EDCfolder="IOBlock",tooltiptext="additional offset of random signal"] 
	double seed; //!EDC$[varaccess,EDCvarname="seed",EDCfolder="IOBlock",tooltiptext="seed € [0.,1.]... initialization of random generator"]  // maximal amplitude, additional offset of random signal, seed € [0.,1.]... initialization of random generator
	TRandomType method;  //!EDC$[varaccess,EDCvarname="method",EDCfolder="IOBlock",tooltiptext="method of random value generation: Trand = 0 ... use built-in command "rand", TLFSR = 1... Linear Feedback Shift Register"]  // method of random value generation
	int rand_max; //!EDC$[varaccess,EDCvarname="rand_max",EDCfolder="IOBlock",tooltiptext="maximum value that can be returned by the rand function."]                   // standard: RAND_MAX, alternative: 2^N-1
	int isConstAmplitude;           // set to 1 if output value should be only +amplitude or -amplitude
	double init_val; //$ MSax 2013-02-28: added
	//old: int random_value;             // value between [0, RAND_MAX], RAND_MAX <=> 0x7fff <=> binary: 0111 1111 1111 1111 ==> stored in XDiscrete(2)

	//EDC TArray<int> inputs;			//$EDC$[varaccess,remove,EDCvarname="input_element_numbers",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed
	//EDC TArray<int> input_types;	//$EDC$[varaccess,remove,EDCvarname="input_element_types",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed
	//EDC TArray<int> input_localnum;	//$EDC$[varaccess,remove,EDCvarname="input_local_number",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed

};//$EDC$[endclass,RandomSource]

//--------------------------------------------------------------------------------------------------------------
// DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS DISCRETE SYSTEMS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: LinearTransformation (InputOutputElementDiscrete)
// short description: Continuous linear transformation 
// available formulations: Continuous state-space elements
// type: SISO (single input-single output), MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description: Continuous linear transformation y = A*u+b
//y = A*u+b; //a..output vector, A...transformation matrix, u..input vector, b=offset vector
//--> linear transformation y = A*u
//--> simple gain y_1 = A_11*u_1
//--> constant y_1 = b_1
// class variables:
//   -A_coeff: transformation matrix A
//   -b_coeff: offset vector
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class LinearTransformation: public InputOutputElement//$EDC$[beginclass,classname=LinearTransformation,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOLinearTransformation,
//texdescription="Continuous linear transformation. The transfer function type is SISO (single input-single output) or MIMO (multi input-multi output).",
//figure="IOLinearTransformation,IOLinearTransformation",
//modus="{linear transformation $y=A\,u$}{Set $\mathbf{b}$ to zero.}",
//modus="{gain $y_1=A_{1,1}u_1$}{Set A as scalar value and b is zero.}",
//modus="{constant $y_1=b_1$}{Set A to zero and b to the constant value.}",
//texdescriptionEquations="
//\begin{equation}
//  \mathbf{y} = \mathbf{A}\mathbf{u}+\mathbf{b}; \\ 
//\end{equation}
//Matrix $\mathbf{A}$ and vector $\mathbf{b}$ are user defined.",
//example="LinearTransformation.txt"]
{
public:

	LinearTransformation(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetLinearTransformation(const Matrix& A, const Vector& b);

	virtual void SetConstant(double val);

	virtual void SetGain(double val);

	virtual void SetAdder(const Vector& signs);

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new LinearTransformation(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor();

	virtual const char* GetElementSpec() const {return "IOLinearTransformation";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const;

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();
	virtual int IsDirectFeedThrough() const 
	{
		if (A_coeff.Getrows() == 0 || A_coeff.Getcols() == 0)
		{
			return 0;
		}
		if (A_coeff.Norm2() == 0)
		{
			return 0;
		}
		return 1;
	}

protected:
	Matrix A_coeff; //$EDC$[varaccess,EDCvarname="A_matrix",EDCfolder="IOBlock",variable_length_vector,tooltiptext="transformation matrix A: y=A.u+b"]
	Vector b_coeff; //$EDC$[varaccess,EDCvarname="b_vector",EDCfolder="IOBlock",variable_length_vector,tooltiptext="offset vector b: y=A.u+b"]
};//$EDC$[endclass,LinearTransformation]


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: Quantizer (InputOutputElement)
// short description: Quantizer
// available formulations: Continuous state-space elements
// type: SISO(single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Continuous linear transformation y = A*u+b
//y = A*u+b; //a..output vector, A...transformation matrix, u..input vector, b=offset vector
//--> linear transformation y = A*u
//--> simple gain y_1 = A_11*u_1
//--> constant y_1 = b_1
// class variables:
//   -A_coeff: transformation matrix A
//   -b_coeff: offset vector
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Quantizer: public InputOutputElement//$EDC$[beginclass,classname=Quantizer,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOQuantizer,
//texdescription="A quantizer block passes its input signal through a stair-step function so that many neighboring points on the input axis are mapped to one point on the output axis. The effect is to quantize a smooth signal into a stair-step output. It is a SISO (single input-single output) control element.",
//figure="IOQuantizer,IOQuantizer",
//texdescriptionEquations="
//\begin{equation}
//y(u) =
//\begin{cases}
// r\,floor\left(\frac{u}{r}+0.5r\right) , & \text{if } r != 0 \\
// u, & \text{if } r = 0
//\end{cases}
//\end{equation}
//The user defined rounding value is $r$.",
//example="Quantizer.txt"]
{
public:

	Quantizer(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetQuantizer(double roundI)
	{
		roundval = roundI;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new Quantizer(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const Quantizer& ce = (const Quantizer&)e;

		roundval = ce.roundval;
	}

	virtual void InitConstructor()
	{
		SetNOutputs(1);
		SetNStates(0);
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		roundval = 0.1; //$ MSax 2013-02-28: added
	}

	virtual const char* GetElementSpec() const {return "IOQuantizer";}
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const
	{
		if (GetNInputs() == 1)
		{
			if (roundval != 0)
				return roundval * floor(GetInput(t, 1)/roundval + 0.5*roundval);
			else
				return GetInput(t, 1);
		}
		else
		{
			assert(0 && "Quantizer has no input!!!");
		}
		assert(0 && "return value added for compiler");		
		return 0;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();
	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}

protected:
	double roundval;	
};//$EDC$[endclass,Quantizer]


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: STransferFunction (InputOutputElement)
// short description: STransferFunction
// available formulations: Continuous state-space elements, linear transfer function
// type: SISO(single input-single output), linear transfer function
// development status: complete, auto-generation functions not complete
// long description: Continuous linear transformation y = A*u+b
// Realization of S-Transfer Function (Laplace Transform):
//    G(s) = (num(1)+num(2)*s + ... + num(n)*s^(n+1)) / (den(1) + den(2)*s + ... + den(n)*s^(n+1)); 
// class variables:
// -num, den: numerator and denominator of transfer function in ascent order of Laplace variable "s"; must have same lenght (fill numerator with zero's if it has lower order than denominator)
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class STransferFunction: public InputOutputElement//$EDC$[beginclass,classname=STransferFunction,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOContinuousTransferFunction,
//texdescription="The STransferFunction is a linear transfer function for continuous state-space elements. It is a SISO (single input-single output) type.",
//figure="IOContinuousTransferFunction,IOContinuousTransferFunction",
//texdescriptionEquations="
//$y(s)=\mathbf{G}(s)u(s)$ \\ \\
//$\mathbf{G}(s)=\frac{\mathbf{num(s)}}{\mathbf{den(s)}}$ \\ \\
//user input: \\ 
//$\mathbf{num}(s) = num_1 + num_2s + num_3s^2 +...+ num_{n+1}s^n$ \\ 
//$\mathbf{den}(s) = den_1 + den_2s + den_3s^2 +...+ den_{n+1}s^n$",
//example="STransferFunction.txt"]
{
public:

	STransferFunction(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetSTransferFunctionWithDescentPolynoms(const Vector& numI, const Vector& denI)
	{
		//-------------------------------------------------------------------
		// G(s) = p_num(s)/p_den(s)
		// p_num(s) = num(m) + num(m-1) * s + num(m-2) * s^2 + ... + num(1) * s^m
		// p_den(s) = den(n) + den(n-1) * s + den(n-2) * s^2 + ... + den(1) * s^n
		// m<=n
		//-------------------------------------------------------------------
		//check if sizes are consistent!
		if(numI.Length() != denI.Length())
		{
			mbs->UO(UO_LVL_err).InstantMessageText("Error#1 while creating S-Transfer function! Numerator and denominator vectors must have same size.\n"); 
			return;
		}

		// revert numerator and denominator
		int len = numI.Length();
		Vector tmp1(len), tmp2(len);
		for(int i = 1; i <= len;i++)
		{				
			tmp1(len - i + 1) = numI(i);
			tmp2(len - i + 1) = denI(i);
		}
		SetSTransferFunction(tmp1, tmp2);

	}
	virtual void SetSTransferFunctionWithDescentPolynoms(const Vector& numI, const Vector& denI, const Vector& initVector)
	{
		//-------------------------------------------------------------------
		// G(s) = p_num(s)/p_den(s)
		// p_num(s) = num(m) + num(m-1) * s + num(m-2) * s^2 + ... + num(1) * s^m
		// p_den(s) = den(n) + den(n-1) * s + den(n-2) * s^2 + ... + den(1) * s^n
		// m<=n
		//-------------------------------------------------------------------
		//check if sizes are consistent!
		if(numI.Length() != denI.Length())
		{
			mbs->UO(UO_LVL_err).InstantMessageText("Error#2 during creating S-Transfer function! Numerator and denominator vectors must have same size.\n"); 
			SetSTransferFunction(Vector(1.0), Vector(1.0), Vector(0.));
			return;
		}

		// revert numerator and denominator
		int len = numI.Length();
		Vector tmp1(len), tmp2(len);
		for(int i = 1; i <= len;i++)
		{				
			tmp1(len - i + 1) = numI(i);
			tmp2(len - i + 1) = denI(i);
		}
		SetSTransferFunction(tmp1, tmp2, initVector);
	}
	virtual void SetSTransferFunction(const Vector& numI, const Vector& denI, const Vector& initVector)
	{
		//check if sizes are consistent!
		if(numI.Length() != denI.Length())
		{
			mbs->UO(UO_LVL_err).InstantMessageText("Error#3 during creating S-Transfer function! Numerator and denominator vectors must have same size.\n"); 
			SetSTransferFunction(Vector(1.0), Vector(1.0), Vector(0.));
			return;
		}

		//inputs must be added separately
		num = numI;
		den = denI;
		
		double fact = den(den.Length());
		if(fact == 0)
		{
			mbs->UO(UO_LVL_err).InstantMessageText("Error#4 while creating S-Transfer function! Denominator coefficient corresponding to highest exponent of Laplace variable is zero.\n"); 
			fact = 1.;
		}

		num *= (1./fact);
		den *= (1./fact);

		SetNStates(den.Length()-1);
		if (GetNStates() == 0) UO(UO_LVL_err).InstantMessageText("Error#5 in S-Transfer Function: is fully algebraic!\n");
		SetInitialValues(initVector);

		//check if sizes are consistent!
		if(GetNStates() != initVector.Length()){UO(UO_LVL_err).InstantMessageText("Error5 in S-Transfer Function: Problems during inizialisation!\n");}
	}

	virtual void SetSTransferFunction(const Vector& numI, const Vector& denI)
	{
		Vector initVector(numI.Length()-1);
		initVector.SetAll(0);

		SetSTransferFunction(numI, denI, initVector);
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new STransferFunction(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const STransferFunction& ce = (const STransferFunction&)e;

		num = ce.num;
		den = ce.den;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		SetNOutputs(1);		
		elementname = GetElementSpec();
		SetSTransferFunction(Vector(1.,0.,0.,0.), Vector(0.,0.,0.,1.0)); // set to valid value
	}

	virtual const char* GetElementSpec() const {return "IOContinuousTransferFunction";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual void EvalF(Vector& f, double t);


	virtual double GetOutput(double t, int i=1) const;

	virtual int IsDirectFeedThrough() const 
	{
		int n = ES();
		if (num.Length() != 0)
		{
			if (num(n+1) != 0) return 1;
		}
		return 0;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	Vector num; //$EDC$[varaccess,variable_length_vector,EDCvarname="numerator",EDCfolder="IOBlock",tooltiptext="ascending numerator coefficients n of transfer-function. TF = num/den with num = n(1)*1+n(2)*s+n(3)*s*s+... Will be normalized automatically!"] 
	Vector den; //$EDC$[varaccess,variable_length_vector,EDCvarname="denominator",EDCfolder="IOBlock",tooltiptext="ascending denominator coeffs d of transfer-function. TF = num/den with den = d(1)*1+d(2)*s+d(3)*s*s+... Will be normalized automatically!"]
};//$EDC$[endclass,STransferFunction]


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: LinearODE (InputOutputElement)
// short description: Linear ODE
// available formulations: Continuous state-space elements, linear ODE
// type: SISO(single input-single output), linear differential equations
// development status: complete, auto-generation functions not complete
// long description: Linear ordinary differential equation
// Realization of Linear ODE
// x_dot = A*x + B*u
// y = C*x + D*u
// class variables:
// -A_coeff, B_coeff, C_coeff, D_coeff: matrices of linear differential equation system
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class LinearODE: public InputOutputElement//$EDC$[beginclass,classname=LinearODE,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOLinearODE,
//texdescription="The LinearODE Element represents a linear ordinary differential equation of SISO (single input-single output) or MIMO (multi input-multi output) type.",
//figure="IOLinearODE,IOLinearODE",
//texdescriptionEquations="
//$\mathbf{\dot{x}}=\mathbf{A}\,\mathbf{x}+\mathbf{B}\,\mathbf{u}$ \\
//$\mathbf{y}=\mathbf{C}\,\mathbf{x}+\mathbf{D}\,\mathbf{u}$ \\ \\
//Matrices $\mathbf{A}$, $\mathbf{B}$, $\mathbf{C}$ and $\mathbf{D}$ are user defined.",
//example="LinearODE.txt"]
{
public:

	LinearODE(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetLinearODE(const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& D, const Vector& initVector)
	{
		//set linear ODE
		int ni = B.Getcols(); //number of inputs
		int ns = A.Getcols(); //number of states
		int no = C.Getrows(); //number of outputs

		assert(A.Getrows() == ns);
		assert(B.Getrows() == ns);
		assert(C.Getcols() == ns);
		assert(D.Getrows() == no);
		assert(D.Getcols() == ni);

		assert(initVector.Length() == ns);

		A_coeff = A;
		B_coeff = B;
		C_coeff = C;
		D_coeff = D;

		SetNOutputs(no);
		SetNStates(ns);
		SetInitialValues(initVector);

		/*UO() << "Linear ODE:\n";
		UO() << "A=" << A_coeff << "\n";
		UO() << "B=" << B_coeff << "\n";
		UO() << "C=" << C_coeff << "\n";
		UO() << "D=" << D_coeff << "\n";
		UO() << "x_init=" << x_init << "\n";*/
	}

	virtual void SetFirstOrderODE(double a, double b, double x_0) //x_dot = a*x + b*u, y=x
	{
		Matrix A(a);
		Matrix B(b);
		Matrix C(1.);
		Matrix D(0.);

		Vector initVector(x_0);

		SetLinearODE(A, B, C, D, initVector);
	}

	virtual void SetSecondOrderODE(double inertia_fact, double damp_fact, double stiff_fact, double input_fact, double x_0, double v_0) //m*x_ddot + d*x_dot + k*x = b*u, y(1) = x, y(2) = x_dot (=v)
	{
		// this defines a linear oscillator, with one input (input gain b), two outputs (x, v)
		//[x] [   0        1   ] [x] [  0 ]
		//[ ]=[                ]*[ ]+[    ] * [u]
		//[v] [-m^-1*k  -m^-1*d] [v] [m^-1]
		//
		//y = C * x + D * u, C=diag(1, 1), D=zero matrix

		assert(inertia_fact != 0);
		double minv = 1./inertia_fact;

		Matrix A(0., 1., -minv*stiff_fact, -minv*damp_fact);

		Matrix B(2,1);
		B(1,1) = 0.;
		B(2,1) = minv*input_fact;

		Matrix C(2,2);
		C.SetDiagMatrix(1.);

		Matrix D(2,1);
		D(1,1) = 0.; D(2,1) = 0.;

		Vector initVector(x_0, v_0);

		SetLinearODE(A, B, C, D, initVector);
	}


	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new LinearODE(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const LinearODE& ce = (const LinearODE&)e;

		A_coeff = ce.A_coeff;
		B_coeff = ce.B_coeff;
		C_coeff = ce.C_coeff;
		D_coeff = ce.D_coeff;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();

		A_coeff = Matrix(1,1);  //$ MSax 2013-02-28 added
		B_coeff = Matrix(1,1);
		C_coeff = Matrix(1,1);
		D_coeff = Matrix(1,1);
	}

	virtual const char* GetElementSpec() const {return "IOLinearODE";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual void EvalF(Vector& f, double t) 
	{
		for (int i=1; i <= A_coeff.Getrows(); i++)
		{
			double val = 0;
			for (int j=1; j <= A_coeff.Getcols(); j++)
			{
				val += A_coeff(i,j)*XG(j);
			}
			for (int j=1; j <= B_coeff.Getcols(); j++)
			{
				val += B_coeff(i,j)*GetInput(t, j);
			}
			f(i) += val;
		}
	};

	virtual double GetOutput(double t, int i=1) const 
	{
		//to be overwritten in specific class!
		if (i <= 0 || i > C_coeff.Getrows()) 
		{
			assert(0);
			return 0;
		}
		//return 0;

		//y(i) = Sum_j {C(i,j)*x(j) + D(i,j)*u(j)}
		double val = 0;
		for (int j=1; j <= C_coeff.Getcols(); j++)
		{
			val += C_coeff(i,j)*XG(j);
		}
		for (int j=1; j <= D_coeff.Getcols(); j++)
		{
			if (D_coeff(i,j) != 0.) 
				val += D_coeff(i,j)*GetInput(t, j);
		}
		return val;
	}
	virtual int IsDirectFeedThrough() const 
	{
		if (D_coeff.Getrows() == 0 || D_coeff.Getcols() == 0)
		{
			return 0;
		}
		if (D_coeff.Norm2() == 0)
		{
			return 0;
		}
		return 1;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	Matrix A_coeff, B_coeff, C_coeff, D_coeff;
};//$EDC$[endclass,LinearODE]


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOMathFunction (InputOutputElement)
// short description: Continuous MathFunction and look-up tables
// available formulations: mathematical functions and look-up tables
// type: continuous element, SISO
// development status: auto-generation functions not complete
// long description: Continuous MathFunction and look-up tables (x-variable: Input, y-variable: Output)
// output = Mathfunction(input); //1 input, output is computed via MathFunction
// class variables:
// - mathfunc: MathFunction
// - elementspec: specification of element; variable used for name
// - pieceWiseSwitchOnlyInPostNewton:  flag for activation / deactivate switching of index only in post newton
// - pieceWiseIndex:  index for inter-(extrapolation) of values
// - pieceWiseIteration: number of iterations (if max. number of iterations is reached, the postnewton error should be zero and this number is reduced)
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IOMathFunction: public InputOutputElement //$EDC$[beginclass,classname=IOMathFunction,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOMathFunction,
//texdescription="A IOMathFunction contains a mathematical expression or a lookup table with different modes for piecewise interpolation. The output is result of the evalutation of the MathFunction as a function of input. It is a SISO (single input-single output) control element.",
//figure="IOMathFunction,IOMathFunction",example="MathFunction.txt",
//modus="{parsed function}{IOBlock.MathFunction.piecewise\_mode must be set to $-1$. In IOBlock.MathFunction.parsed\_function one specifies a string representing parsed function, e.g. '$A*sin(u)$' with funtion parameter $u$ defined in IOBlock.MathFunction.parsed\_function\_parameter.}",
//modus="{piecewise mode - constant}{IOBlock.MathFunction.piecewise\_mode must be set to $0$. The vectors IOBlock.MathFunction.piecewise\_points and IOBlock.MathFunction.piecewise\_values are used. The output value is piecewise constant with jumps at the supporting points.}",
//modus="{piecewise mode - linear}{IOBlock.MathFunction.piecewise\_mode must be set to $1$. The vectors IOBlock.MathFunction.piecewise\_points and IOBlock.MathFunction.piecewise\_values are used. The output value is piecewise linear between the supporting points.}",
//modus="{piecewise mode - quadratic}{IOBlock.MathFunction.piecewise\_mode must be set to $2$ and in addition to the other piecwise modes the vector IOBlock.MathFunction.piecewise\_diff\_values is needed. The output is a quadratic interpolation between the supporting points.}"]
{
public:

	IOMathFunction(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	virtual void SetIOMathFunction(const MathFunction& mathfuncI);

	virtual void SetSwitchOnlyInPostNewton(int val = 1){pieceWiseSwitchOnlyInPostNewton = val;} //$ RL 2012-7-25: 1... change interval after newton computation and evaluate new jacobi matrix in case of new interval, activate this flag only if the IOMathFunction is no function of time (TItime)

	//To be overwritten in derived class:
	virtual Element* GetCopy();

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor();

	virtual const char* GetElementSpec() const {return "IOMathFunction";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const;

	virtual int DataS() const {return 0;} //Data size for non-state variables (length of XData-vector)

	virtual double PostNewtonStep(double t);

	virtual void IOMathFunction::PostprocessingStep();

	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	MathFunction mathfunc;
	mystr elementspec;                   //EDC$[varaccess,readonly,EDCfolder="",EDCvarname="elementspec",tooltiptext="element specification"]
	int pieceWiseSwitchOnlyInPostNewton; // flag for activation / deactivate switching of index only in post newton
	int pieceWiseIndex;                  // index for inter-(extrapolation) of values
	int pieceWiseIteration;              // number of iterations (if max. number of iterations is reached, the postnewton error should be zero and this number is reduced)
}; //$EDC$[endclass,IOMathFunction]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOSaturate (InputOutputElement)
// short description: Saturation
// available formulations: Continuous state-space elements
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Continuous Saturation element
// class variables:
// -lowerLimit: lower limit; values below are saturated
// -upperLimit: upper limit; values above are saturated
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IOSaturate: public InputOutputElement//$EDC$[beginclass,classname=IOSaturate,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOSaturate,
//texdescription="Continuous saturation element for upper and lower limits. It is a SISO (single input-single output) control element.",
//figure="IOSaturate,IOSaturate",
//texdescriptionEquations="
//\begin{equation}
//y(u) =
//\begin{cases}
// ul, & \text{if } u>ul \\
// u, & \text{if } ll\leq u \leq ul \\
// ll, & \text{if } u<ll
//\end{cases}
//\end{equation}
//In the defined equation $ul$ is the upper limit and $ll$ is the lower limit.",
//example="Saturate.txt"]
{
public:

	IOSaturate(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	// lower/-upperLimit... lower/upper limits of saturation; max_val.. 
	virtual void SetIOSaturate(double lowerLimiti, double upperLimiti)
	{
		lowerLimit = lowerLimiti;
		upperLimit = upperLimiti;

		SetNOutputs(1);
		SetNStates(0);
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOSaturate(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOSaturate& ce = (const IOSaturate&)e;
		lowerLimit = ce.lowerLimit;
		upperLimit = ce.upperLimit;

	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		lowerLimit = 0;	//$ MSax 2013-02-28 added
		upperLimit = 0.1;	//$ MSax 2013-02-28 added
	}

	virtual double GetOutput(double t, int i=1) const
	{
		double u = GetInput(t, i);
		if(u < lowerLimit)return lowerLimit;
		if(u > upperLimit)return upperLimit;
		return u;
	}

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual const char* GetElementSpec() const {return "IOSaturate";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual const char* IOSaturate::SymbolText() const;
	virtual void DrawBlockSymbol();

protected: 
	double lowerLimit;
	double upperLimit;
};//$EDC$[endclass,IOSaturate]
//$ RL 2011-11:] IOSaturate

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IODeadZone (InputOutputElement)
// short description: Deadzone
// available formulations: Continuous state-space elements
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Continuous Deadzone element
//    The outputs between upper and lower limit is the value zero.
//    This leads to an offset of the input signal by the corresponding lower or upper limit.
// class variables:
// -start_deadzone: lower limit of deadzone
// -end_deadzone: upper limit of deadzone
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IODeadZone: public InputOutputElement//$EDC$[beginclass,classname=IODeadZone,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IODeadZone,
//texdescription="Continuous dead-zone element. The outputs between upper and lower limit is zero. This leads to an offset of the input signal by the corresponding lower or upper limit. It is a SISO (single input-single output) control element.",
//figure="IODeadZone,IODeadZone",
//texdescriptionEquations="
//\begin{equation}
//y(u) =
//\begin{cases}
// u-sd, & \text{if } u<sd \\
// 0, & \text{if } u \geq sd \text{ and } u \leq ed \\
// u-ed, & \text{if } u>ed
//\end{cases}
//\end{equation}
//In the defined equation $sd$ is the start dead-zone value, $ed$ is the end dead-zone value.",
//example="DeadZone.txt"]
{

public:

	IODeadZone(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
		SetNStates(0);
		SetNOutputs(1);
	};

	void SetIODeadZone(double start_deadzoneI, double end_deadzoneI)
	{
		assert(start_deadzoneI < end_deadzoneI);
		start_deadzone = start_deadzoneI;
		end_deadzone = end_deadzoneI;		
	}
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IODeadZone(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IODeadZone& ce = (const IODeadZone&)e;
		start_deadzone = ce.start_deadzone;
		end_deadzone = ce.end_deadzone;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		start_deadzone = 0.;
		end_deadzone = 0.;
	}

	virtual const char* GetElementSpec() const {return "IODeadZone";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	//
	virtual double GetOutput(double t, int i=1) const 
	{
		double u = GetInput(t, i);
		if(u < start_deadzone)
		{
			return u - start_deadzone;
		}
		else if(u > end_deadzone)
		{
			return u - end_deadzone;
		}
		else 
		{
			return 0.; //u >= start_deadzone && u =< end_deadzone
		}
	}

	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}


	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	double start_deadzone; // lower limit of deadzone
	double end_deadzone;  // upper limit of deadzone
};//$EDC$[endclass,IODeadZone]
//$ RL 2011-01:] IODeadzone added

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOProduct (InputOutputElement)
// short description: Product
// available formulations: Continuous state-space elements
// type: SISO (single input-single output), MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description: Continuous product (or division) of one or more inputs
//    y=u1^exp1*u2^exp2*...*un^expn+offset, y...output, u...input
//    exp = -1 ... multiply inverse of input
//    exp =  1 ... multiply input
// class variables:
// -exp: exponent of inputs, see formula in long description
// -offset: offset, see formula in long description
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IOProduct: public InputOutputElement//$EDC$[beginclass,classname=IOProduct,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOProduct,
//texdescription="Continuous product (or division) of one or more inputs. A dedicated exponent for every input and a offset can be applied.",
//figure="IOProduct,IOProduct",
//texdescriptionEquations="
//\begin{equation}
//y(\mathbf{u}) = u_1^{exp_1}u_2^{exp_2}...u_n^{exp_n}+offset
//\end{equation}
//All exponents are stored in a vector. For a simple multiplication with a input the dedicated exponent is set to 1, for a division the exponent is set to -1. The offset is a scalar value.",
//example="Product.txt"]
{
public:
	IOProduct(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
		SetNStates(0);
		SetNOutputs(1);
	};
	// exp (x) ... 1 = multiplicator, -1 ... division
	void SetIOProduct(Vector& expi, double offseti = 0.)
	{
		exp = expi;
		if(exp.Length() <= 0)
		{
			exp = Vector(1.,1.);
			mbs->UO(UO_LVL_err).InstantMessageText("Error in SetIOProduct: exponent vector is empty.");
		}
		offset = offseti;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOProduct(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOProduct& ce = (const IOProduct&)e;
		exp = ce.exp;
		offset = ce.offset;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();

		exp = Vector(1);  //$ MSax 2013-02-28
		offset = 0;  //$ MSax 2013-02-28
	}

	virtual const char* GetElementSpec() const {return "IOProduct";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1; /*MSax 2013-02-28 changed from 1 to -1*/}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual double GetOutput(double t, int i=1) const;

	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	Vector exp;    // exponents of inputs
	double offset; // output offset
};//$EDC$[endclass,IOProduct]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOTime (InputOutputElement)
// short description: Time source
// available formulations: Continuous state-space elements
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Continuous time source
//   This element simply outputs the time.
// class variables:
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IOTime: public InputOutputElement //$EDC$[beginclass,classname=IOTime,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOTime,
//texdescription="Continuous time source. This element simply outputs the time.",
//figure="IOTime,IOTime",example="addTime.txt"]
{
public:

	IOTime(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOTime(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOTime& ce = (const IOTime&)e;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		SetNStates(0);
		SetNOutputs(1);
	}

	virtual const char* GetElementSpec() const {return "IOTime";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 0;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------


	virtual double GetOutput(double t, int i=1) const 
	{
		return t;
	}
	virtual int IsDirectFeedThrough() const 
	{
		return 0;
	}


	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();
protected:

	//EDC TArray<int> inputs;			//$EDC$[varaccess,remove,EDCvarname="input_element_numbers",EDCfolder="IOBlock"] // $ MSax 2013-08-28: removed
	//EDC TArray<int> input_types;	//$EDC$[varaccess,remove,EDCvarname="input_element_types",EDCfolder="IOBlock"] // $ MSax 2013-08-28: removed
	//EDC TArray<int> input_localnum;	//$EDC$[varaccess,remove,EDCvarname="input_local_number",EDCfolder="IOBlock"] // $ MSax 2013-08-28: removed

}; //$EDC$[endclass,IOTime]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOPulseGenerator (InputOutputElement)
// short description: Pulse Generator
// available formulations: Continuous state-space elements
// type: SISO (single input-single output)
// development status: complete, auto-generation functions not complete
// long description: Continuous Pulse Generator
//  This element outputs repeating sequence or rectangular pulses after a certain delay
//  parameters: amplitude                pulse width (s)
//                                         ________    ________    __ ....
//                                         |       |   |       |   |
//                            t=0__________|       |___|       |___|
//                              phase delay (s)
//                                           period (s)
//                                         |<--------->|
// class variables:
// -amplitude: amplitude of rectangle pulse
// -toffs: offset time for pulse sequence
// -period: period of rectangle
// -pulseWidth: pulse width - output is set to amplitude in this time span
// -useExternalTime: set to nonzero value if input should be used as time source; otherwise use simulation time
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class IOPulseGenerator: public InputOutputElement//$EDC$[beginclass,classname=IOPulseGenerator,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOPulseGenerator,
//texdescription="Continuous pulse generator. This element outputs repeating sequence or rectangular pulses after a certain delay. It has no input and one output.",
//figure="IOPulseGenerator,IOPulseGenerator",
//texdescriptionEquations="
//\begin{equation}
//\Delta t = t-t_{offset}
//\end{equation}
//\begin{equation}
//t_{rest} = \Delta t \text{ mod } p
//\end{equation}
//\begin{equation}
//y(t) =
//\begin{cases}
// a, & \text{if } \Delta t \geq 0 \text{ and } t_{rest} < pw \\
// 0, & \text{else }
//\end{cases}
//\end{equation}
//User defined variables are pulse amplitude $a$, time offset $t_{offset}$, signal period $p$ and pulse width $pw$.",example="addPulseGenerator.txt"
//]

{
public:

	IOPulseGenerator(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
		SetNStates(0);
		SetNOutputs(1);
	};

	void SetIOPulseGenerator(double amplitudeI, double toffsI, double periodI,	double pulseWidthI, int useExternalTimeI = 0)
	{
		amplitude = amplitudeI;
		toffs = toffsI;
		period = periodI;
		pulseWidth = pulseWidthI;
		useExternalTime = useExternalTimeI;
		assert(period > 0. && pulseWidth >= 0.);
	}	

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOPulseGenerator(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOPulseGenerator& ce = (const IOPulseGenerator&)e;
		amplitude = ce.amplitude;
		toffs = ce.toffs;
		period = ce.period;
		pulseWidth = ce.pulseWidth;
		useExternalTime = ce.useExternalTime;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();

		amplitude = 1; // $ MSax 2013-02-28: added
		toffs = 0;
		period = 1;
		pulseWidth = 0.5;
		useExternalTime = 0;
	}

	virtual const char* GetElementSpec() const {return "IOPulseGenerator";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 0;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual double GetOutput(double t, int i=1) const 
	{
		double time;
		if(useExternalTime)time = GetInput(t, i);
		else time = t;

		double tt = time-toffs; // shifted time (without offset)
		double trest = fmod(tt,period);

		if(tt >= 0. && trest < pulseWidth)
		{
			return amplitude;
		}
		else
		{
			return 0.;
		}
	}
	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

protected:
	double amplitude;    // amplitude of rectangle pulse
	double toffs;        // offset time for pulse sequence
	double period;       // period of rectangle
	double pulseWidth;   // pulse width - output is set to amplitude in this time span
	int useExternalTime; // set to nonzero value if input should be used as time source; otherwise use simulation time

	//EDC TArray<int> inputs;			//$EDC$[varaccess,remove,EDCvarname="input_element_numbers",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed
	//EDC TArray<int> input_types;	//$EDC$[varaccess,remove,EDCvarname="input_element_types",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed
	//EDC TArray<int> input_localnum;	//$EDC$[varaccess,remove,EDCvarname="input_local_number",EDCfolder="IOBlock"] // $ MSax 2013-04-23: removed

};//$EDC$[endclass,IOPulseGenerator]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: ControllerInterfaceData (Element)
// short description: Superclass for specific controller circuits
// available formulations: Discontinuous state-space elements
// type: SISO (single input-single output), MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description: Superclass for discontinuous controller interface
//  Specific time discrete input-output elements can be implemented as subclass of ControllerInterfaceData (e.g. from external code)
//  The function x(k+1)=f(xk,uk) is called every discrete event from the class "ControllerInterface". 
//  Example: This Interface contains pointers to global update and initialize functions, sample time and number of element outputs.
//   
//   The following functions have to be implemented in the subclass:
//     virtual const double& Uk(int num=1) const;
//     virtual double& Uk(int num=1);
//     virtual void SetUk(double valI, int num=1);
//     virtual double GetOutput(double t, int i=1) const;
//
// class variables:
//	-void (*Reg_Init)(): pointer on initialize function
//	-void (*Reg_Update)(): pointer on update function
//	-deltaT: Reg_Update is called every deltaT seconds
//	-nOut: number of outputs
//	-nInp: number of inputs
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ControllerInterfaceData: public Element//!EDC$[beginclass,classname=ControllerInterfaceData,parentclassname=Element]
{
public:
	ControllerInterfaceData(MBS* mbsi) : Element(mbsi)
	{
		//problem with debug modus if following two lines are no comment
		//Reg_Init = NULL;			 
		//Reg_Update = NULL; 
		deltaT = 0.0;
		nOut = 0;
		nInp = 0;
	}

	ControllerInterfaceData(const ControllerInterfaceData& e) : Element(e.mbs)	{CopyFrom(e);}

	//  pointer on global functions can be added with &function, e.g. Controller from simulink
	virtual void SetControllerInterfaceData(void (*Reg_InitI)(), void (*Reg_UpdateI)(), double deltaTI, int nOutI, int nInpI)
	{
		Reg_Init = Reg_InitI;
		Reg_Update = Reg_UpdateI;
		deltaT = deltaTI;
		nOut = nOutI;
		nInp = nInpI;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ControllerInterfaceData(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e)
	{
		Element::CopyFrom(e);
		const ControllerInterfaceData& ce = (const ControllerInterfaceData&)e;
		Reg_Init = ce.Reg_Init;
		Reg_Update = ce.Reg_Update;
		deltaT = ce.deltaT;		
		nOut = ce.nOut;
		nInp = ce.nInp;
	}  

	virtual void Controller_Init()  {(*Reg_Init)();}
	virtual void Controller_Update(){(*Reg_Update)();}

	virtual int GetNOutputs() const{return nOut;}
	virtual int GetNInputs() const{return nInp;}

	virtual void SetSampleTime(double val){		deltaT = val;	};
	virtual double GetSampleTime() const{return deltaT;};

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	//virtual void GetElementDataAuto(ElementDataContainer& edc);
	//virtual int SetElementDataAuto(ElementDataContainer& edc);	
	//virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	//virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	//virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	//-----------------------------------------------------
	//b: EXAMPLE INTERFACE DATA (IN CHILD CLASS)
	////implement input and outputs in child class; example: 
	////datah is child of class ControllerInterfaceData with specific data storage of inputs and outputs
	////function is called every discrete event
	//// xk+1 = f(xk,uk)
	// virtual const double& Uk(int num=1) const;
	// virtual double& Uk(int num=1);
	// virtual void SetUk(double valI, int num=1);
	// virtual double GetOutput(double t, int i=1) const;
	//e: EXAMPLE INTERFACE DATA (IN CHILD CLASS)
	//-----------------------------------------------------

private:
	void (*Reg_Init)();	  // pointer on initialize function
	void (*Reg_Update)();	// pointer on update function
	double deltaT;				// Reg_Update is called every deltaT seconds
	int nOut;							// number of outputs
	int nInp;             // number of inputs
};//!EDC$[endclass,ControllerInterfaceData]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: ControllerInterface (InputOutputElementDiscrete)
// short description: Interface for time discrete input output-elements
// available formulations: Discontinuous state-space elements
// type: SISO (single input-single output), MIMO (multi input-multi output)
// development status: complete, auto-generation functions not complete
// long description: Interface for time discrete input output-elements (e.g. controller circuits) using special ControllerInterfaceData
//  This interface element executes specific time discrete input-output elements in form of a subclass of ControllerInterfaceData (e.g. from external code).
//  The function x(k+1)=f(xk,uk) is called every discrete event from this interface.
// class variables:
//	-ControllerInterfaceData data: data of special subclass derived from "ControllerInterfaceData" with user defined input-outputs.
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ControllerInterface: public InputOutputElementDiscrete//!EDC$[beginclass,classname=ControllerInterface,parentclassname=InputOutputElementDiscrete]
{
public: //$ RL 2010-04

	ControllerInterface(ControllerInterfaceData& dataI): InputOutputElementDiscrete(dataI.GetMBS()), data(dataI.GetMBS())
	{	
		InitConstructor();
		data.CopyFrom(dataI);
	};

	ControllerInterface(MBS* mbsi) : InputOutputElementDiscrete(mbsi), data(mbsi)	{	};

	ControllerInterface(const ControllerInterface& e) : InputOutputElementDiscrete(e.mbs), data(e.mbs)	{CopyFrom(e);}

	virtual void SetControllerInterface(const ControllerInterfaceData& dataI)		 //old: const Vector& initVector
	{
		InitConstructor();
		data.CopyFrom(dataI); // copy data container

		//SetNOutputs(data.nOut);
		SetNStates(0);                	
		int dimIODiscrete = DataS(); // dimension of root element
		SetNDiscreteStates(0);       //old: SetNDiscreteStates(initVector.Length())
		InputOutputElementDiscrete::SetSampleTime(data.GetSampleTime());		 //sample time already stored
		SetNOutputs(data.GetNOutputs());

		Vector init;
		init.SetLen(0);                    // new init vector inclusive base element variables
		init = init.Append(GetDataInit());

		Vector initVector;		
		initVector.SetLen(0); 
		init = init.Append(initVector);	   // initial vector xk(t=0) (could be inserted here if necessary)

		Vector initU;
		initU.SetLen(data.GetNInputs());
		init = init.Append(initU); // u0=0.0

		SetDataInit(init);
		data.Controller_Init(); 		
	}

	virtual void SetSampleTime(double val){		deltaT = val;	};
	virtual double GetSampleTime(){return	data.GetSampleTime();	};

	virtual int DataS() const {return 1+GetNDiscreteStates()+ data.GetNInputs();}    // dim(k)=1, dim(xk)=n, dim(uk) = 1 --> 2+n

	virtual void SetUk(double valI, int num = 1)	{		XData(1+GetNDiscreteStates()+num) = valI;	}		//set num'th input u	
	virtual double GetUk(int num = 1) const {		return XData(1+GetNDiscreteStates()+num);	}   //get num'th input u	

	// set uk-vector
	virtual void UpdateDiscreteElementInput(double t)
	{		
		for(int i=1; i<=data.GetNInputs(); i++)
		{
			SetUk(GetInput(t, i), i); 
		}
	};

	//function is called every discrete event
	// xk+1 = f(xk,uk)
	virtual void UpdateDiscreteElementState()	 	{	GetMBS()->UO().InstantMessageText("ControllerInterface_UpdateDiscreteElementState_Error: Define function in child class.\n");assert(0);};//	data.Controller_Update();	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new ControllerInterface(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElementDiscrete::CopyFrom(e);
		const ControllerInterface& ce = (const ControllerInterface&)e;
		data.CopyFrom(ce.data);		
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();		
	}

	virtual const char* GetElementSpec() const {return "ControllerInterface";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	//virtual void GetElementDataAuto(ElementDataContainer& edc);
	//virtual int SetElementDataAuto(ElementDataContainer& edc);	
	//virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	//virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	//virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual double GetOutput(double t, int i=1) const { GetMBS()->UO().InstantMessageText("ControllerInterface_GetOutput_Error: Define function in child class.\n");assert(0);	return 0;	}

	virtual int IsDirectFeedThrough() const 	{		return 0;	}

	virtual void DrawElement() 	{		InputOutputElement::DrawElement();	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

	// --------------- EXAMPLE INTERFACE (IN CHILD CLASS)---------------------
	////implement input and outputs in child class; example: 
	////datah is child of class ControllerInterfaceData with specific data storage of inputs and outputs
	////function is called every discrete event
	//// xk+1 = f(xk,uk)
	//virtual void UpdateDiscreteElementState()	 	{	datah.Controller_Update();	}

	//virtual double GetOutput(double t, int i=1) const
	//{		
	//	return datah.GetOutput(double t, i);
	//}
	//
	//virtual void SetUk(double valI, int num = 1)	//set num'th input u
	//{		
	//	ControllerInterface::SetUk(valI, num);	// HOTINT storage
	//	datah.SetUk(valI, num);			            // controller specific storage
	//}
	// ------------------------------------------------------------------------

protected:
	ControllerInterfaceData data; // contains data for interface	
};//!EDC$[endclass,ControllerInterface]

// INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// IOTimeWindow
class IOTimeWindow: public InputOutputElement//$EDC$[beginclass,classname=IOTimeWindow,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOTimeWindow,
//texdescription="This element helps to capture a special time window. It has two inputs and one output.",figure="IOTimeWindow,IOTimeWindow",
//modus="{$t_{end} > t_{start}$}{Output is determined with inequation (a).}",
//modus="{$t_{end} \leq t_{start}$}{Output is determined with inequation (b).}",
//texdescriptionEquations="
//\begin{equation}
//\left(a\right) \qquad y(\mathbf{u}) =
//\begin{cases}
// u_2, & \text{if } t_{start} \leq u_1 \leq t_{end} \\
// 0, & \text{else }
//\end{cases}
//\end{equation}
//\begin{equation}
//\left(b\right) \qquad y(\mathbf{u}) =
//\begin{cases}
// u_2, & \text{if } t_{start} \leq u_1 \\
// 0, & \text{else }
//\end{cases}
//\end{equation}",example="TimeWindow.txt"]

{//$ DR 2011-12:[ IOTimeWindow added

public:

	IOTimeWindow(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	// y = u(2) if t_start <= u(1) <= t_end, y = 0 else
	virtual void SetIOTimeWindow(const double t_starti, const double t_endi)
	{
		t_start = t_starti;
		delta_t = t_endi - t_starti;
		SetNOutputs(1);
		SetNStates(0);
		reached_end = 0;
	}

	// y = u(2) if t_start <= u(1), y = 0 else
	virtual void SetIOTimeWindow(const double t_starti)
	{
		t_start = t_starti;
		delta_t = -1;
		SetNOutputs(1);
		SetNStates(0);
		reached_end = 0;
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOTimeWindow(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOTimeWindow& ce = (const IOTimeWindow&)e;

		t_start = ce.t_start;
		delta_t = ce.delta_t;
		t_end = ce.t_end; // $MSax 2013-03-01: added
		reached_end = ce.reached_end;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		t_start = 0; // $MSax 2013-03-01: added
		t_end = 0; // $MSax 2013-03-01: added
	}

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual const char* GetElementSpec() const {return "IOTimeWindow";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 2;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	virtual double GetOutput(double t, int i=1) const;

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();
	virtual int IsDirectFeedThrough() const 
	{
		return 1;
	}

protected:
	double t_start, delta_t, t_end; //$ MSax 2013-03-01: added t_end for script language
	mutable bool reached_end;
};//$EDC$[endclass,IOTimeWindow]
//$ DR 2011-12:] IOTimeWindow added

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ RL 2012-2-1: //$ MS 2012-2-1:[ this element stops the computation, if input is unequal zero
class IOStopComputation: public InputOutputElement//$EDC$[beginclass,classname=IOStopComputation,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOStopComputation,
//texdescription="This element stops the computation, if input is unequal zero. It has one input and no output.",figure="IOStopComputation,IOStopComputation",example="StopComputation.txt"]
{
public:

	IOStopComputation(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOStopComputation(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOStopComputation& ce = (const IOStopComputation&)e;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		SetNStates(0);
		SetNOutputs(0);
	}

	virtual const char* GetElementSpec() const {return "IOStopComputation";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual double GetOutput(double t, int i=1) const {	return 0.;} //dummy, not used yet

	virtual int IsDirectFeedThrough() const 
	{
		return 0;
	}


	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();

	void StartTimeStep()	
	{
		if (GetInput(mbs->GetTime(),1)>0)
		{
			//GetMBS()->TIFinished();// computation stops, if input is positive
			GetMBS()->StopByElement();// computation stops, if input is positive
		}
	}

protected:
};//$EDC$[endclass,IOStopComputation]
//$ RL 2012-2-1: //$ MS 2012-2-1:]
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: Modifier (InputOutputElement)
// short description: Modfier is used to change element data variables during static or dynamic simulation
// available formulations:
// type: continuous element
// development status: missing check for closed loops of evaluation (if modifier points to same object data as input)
// long description: 
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$ MSax 2013-1:[ added element data modifier
class IOElementDataModifier: public InputOutputElement //$EDC$[beginclass,classname=IOElementDataModifier,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOElementDataModifier,
//texdescription="This element can be used to modify data of a constraint or element. It has one input and no output.",figure="IOElementDataModifier",example="ElementDataModifier.txt"]
{
public:

	IOElementDataModifier(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};

	void SetIOElementDataModifier(int element_numberI, const char* variable_nameI)
	{
		AddElement(element_numberI); //add element to constraints

		successfully_converted = RWdata.GetVariableNameAndComponents(variable_nameI); //convert variable name into name and components in RWdata structure
		RWdata.RWaccess = TRWElementDataWrite; //JG 2013-01-11: we only want to read values in sensors

		if (!successfully_converted) 
		{
			UO() << "ERROR: IOElementDataModifier: variable name='" << variable_nameI << "' is not legal. Modifier is inactive.\n";
		}
		element_number = element_numberI;
		variable_name = mystr(variable_nameI);
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy();

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor();

	virtual const char* GetElementSpec() const {return "IOElementDataModifier";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};
	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol();
	virtual void ModifyAction() //call this function to execute modification
	{
		if (successfully_converted) 
		{
			double time = GetMBS()->GetStepEndTime(); //modification with e.g. a mathfunction should be made with evaluated function at end of time step

			RWdata.value = GetInput(time);
			GetElem(1).WriteSingleElementData(RWdata);
		}
	}

	virtual void StartTimeStep() //modification at the very beginning of the time step (note that modifier elements should be sorted to the beginning of element list)
	{
		if (successfully_converted) ModifyAction();
	}
	virtual void PrecomputeEvalFunctions() //modification at the very beginning of every iteration of static of dynamic simulation
	{
		if (successfully_converted && !modify_at_start_time_step_only) ModifyAction();
	}

protected:
	ReadWriteElementDataVariableType RWdata; //this is the structure for the write access
	int successfully_converted; //store, if elementdata variable name has been successfully converted

	int element_number;
	mystr variable_name; 

	int modify_at_start_time_step_only; //$EDC$[varaccess,int_bool,EDCvarname="start_of_timestep_only",EDCfolder="IOBlock",tooltiptext="modify element data at start time step only."]

}; //$EDC$[endclass,IOElementDataModifier]
//$ MSax 2013-1:]



class IODisplay: public InputOutputElement //$EDC$[beginclass,classname=IODisplay,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IODisplay,
//texdescription="This element can be used to display any (single) numberical value fed into the (single) input.",
//figure="IODisplay,IODisplay",example="Display.txt"]
{
public:
	IODisplay(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
		Vector datainit(DataS());		// $MSax 2013-04-02: added
		datainit.SetAll(0.);				// $MSax 2013-04-02: added
		SetDataInit(datainit);			// $MSax 2013-04-02: added
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IODisplay(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IODisplay& ce = (const IODisplay&)e;
		ndigits = ce.ndigits;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		SetNStates(0);
		SetNOutputs(0);

		draw_dim = Vector3D(3*draw_dim.X(),draw_dim.Y(),0.);
		SetNDigits(3);
	}

	virtual int DataS() const { return 1; } // for initialization of XData vector ( remember the current_value for all times ) $MSax 2013-04-02: added

	virtual const char* GetElementSpec() const {return "IODisplay";}
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual double GetOutput(double t, int i=1) const {return 0.;}

  virtual void EndTimeStep() // $ MSax 2013-04-02: added
	{
			XData(1) = GetInput(GetMBS()->GetStepEndTime());
	}

	virtual int IsDirectFeedThrough() const 
	{
		return 0;
	}

	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();

	virtual void SetNDigits(int n) { ndigits = n; }  
protected:
	int ndigits;                                        //$EDC$[varaccess,EDCvarname="number_of_digits",EDCfolder="IOBlock",tooltiptext="number of digits"]

}; //$EDC$[endclass,IODisplay]


class IOResponseElement: public InputOutputElement //$EDC$[beginclass,classname=IOResponseElement,parentclassname=InputOutputElement,addelementtype=TAEinput_output+TAENotInRelease,addelementtypename=IOResponseElement,
//texdescription="This element allows the IOBlock to respond to a key or mouse input. ( e.g. pressing '+' or '-' in the IOBlockView Window ) ",
//figure="IOResponseElement,IOResponseElement"]
{
public:
	typedef enum { TMouseResponse = 1, TKeyResponse = 2 } TGuiResponseInputMode; 
	typedef enum { TToggle = 1, TSlider = 2 } TGuiResponseRangeMode;

public:
	IOResponseElement(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOResponseElement(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOResponseElement& ce = (const IOResponseElement&)e;
		current_value = ce.current_value;
		lower_bound = ce.lower_bound;
		upper_bound = ce.upper_bound;
		increment = ce.increment;
		input_mode = ce.input_mode;
		range_mode = ce.range_mode;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		SetNStates(0);
		SetNOutputs(1);
// default mode: mouse response - Off/On
		current_value = 0.;
		lower_bound = 0.;
		upper_bound = 1.;
		increment = 1.;
		range_mode = TToggle;

		InitXData();
	}

	virtual const char* GetElementSpec() const {return "IOResponseElement";}
	virtual int GetExpectedNumberOfInputs(){return 0;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);			//set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	// functinos to store the current value in XData
	virtual int DataS() const { return 1; } // for initialization of XData vector ( remember the current_value for all times )
	virtual void InitXData() 
	{ 
		Vector datainit(DataS());
		datainit(1) = current_value;
		SetDataInit(datainit);
	}
	virtual void WriteToXData()
	{ 
		XData(1) = current_value; 
	}
	virtual void Initialize() 
	{ 
		if (lower_bound == 0. && upper_bound == 1. && increment == 1.) 
			range_mode = TToggle;
		else range_mode = TSlider;

		WriteToXData();							// set XData after System is assembled
	} 
	virtual void ApplyChange()
	{
		SimulationStatus status = GetMBS()->GetSimulationStatus(); //$ MaSch 2013-08-19
	
		//if ( status == TSimulationNotStarted || status == TSimulationEndedRegularly)
		if ( status.GetStatusFlag(TSimulationNotStarted) || status.GetStatusFlag(TSimulationEndedRegularly) ) //$ MaSch 2013-08-19
			WriteToXData();   // Apply data imediately ONLY if the simulation has not been started ... 
	}
	virtual void StartTimeStep() 
	{ 
		WriteToXData();            // when the simulation is running, apply changes at the beginning of a step ( no iterations etc.. )
	} 

	virtual double GetOutput(double t, int i=1) const {	return XData(1); /*return current_value;*/ }
	virtual int IsDirectFeedThrough() const {	return 0;	}

	virtual const char* SymbolText() const;
  virtual void DrawBlockSymbol();

	// Get Funcitons
	virtual TGuiResponseInputMode GetInputMode() { return input_mode; }																						// response on mouse or key
	virtual double GetCurrentValue() { return current_value; }																						// get the current value
	virtual int IsToggle() { if (range_mode == TToggle) return true; else return false; }																																									// Toggle or Slide

	// Set Functions without update ( can be used before Assemble )
	virtual void SetInputMode(TGuiResponseInputMode mode) { input_mode = mode; }
	virtual void SetIncrement(double inc) { increment = inc; } 
	virtual int SetCurrentVaule(double val) 
	{ 
		if (lower_bound > val || val > upper_bound) { return false; } // conflict range <-> current_value
		else { current_value = val;  return true; }
	}
	virtual int SetRange(double lower, double upper)
	{
		if ( lower > current_value || current_value > upper ) { return false; } // conflict range <-> current_value
		else { lower_bound = lower; upper_bound = upper; return true; }
	}

	// Set Functions with update ( use only after Assemble - these are the resoponse to the GUI )
	virtual void Increase() { if (!IsToggle()) SetCurrentVaule(GetCurrentValue()+1); else FlipState(); ApplyChange(); }
	virtual void Decrease() { if (!IsToggle()) SetCurrentVaule(GetCurrentValue()-1); else FlipState(); ApplyChange(); }
	virtual void FlipState() { if (GetCurrentValue() == 0) SetCurrentVaule(1); else SetCurrentVaule(0); ApplyChange(); }

	// functions to catch the input - implemented in derived class
	virtual int KeyPressed(int key) { return false; }


protected:
	double current_value;
	double lower_bound, upper_bound, increment;
	TGuiResponseInputMode input_mode;
	TGuiResponseRangeMode range_mode;
}; //$EDC$[endclass,IOResponseElement]

class IOKeyResponseElement: public IOResponseElement //$EDC$[beginclass,classname=IOKeyResponseElement,parentclassname=IOResponseElement,addelementtype=TAEinput_output+TAENotInRelease,addelementtypename=IOKeyResponseElement,
//texdescription="This element allows the IOBlock to respond to a key input. ( e.g. pressing '+' or '-' ) ",
//figure="IOKeyResponseElement,IOKeyResponseElement"]
{
public:
	IOKeyResponseElement(MBS* mbsi):IOResponseElement(mbsi)
	{	
		InitConstructor();
	};

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOKeyResponseElement(*this);
		return ec;
	}
		//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		IOResponseElement::CopyFrom(e);
		const IOKeyResponseElement& ce = (const IOKeyResponseElement&)e;
		inc_key = ce.inc_key;
		dec_key = ce.dec_key;
	}

// copied from winuser.h
// chould be <included>
#define VK_OEM_PLUS       0xBB   // '+' any country
#define VK_OEM_MINUS      0xBD   // '-' any country

	virtual void InitConstructor()
	{
		IOResponseElement::InitConstructor();
		input_mode = IOResponseElement::TKeyResponse;
		inc_key = VK_OEM_PLUS;
		dec_key = VK_OEM_MINUS;
	}
	virtual const char* GetElementSpec() const {return "IOKeyResponseElement";}

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	// set the keys to respond to
	virtual void SetIncKey(int key) { inc_key = key; }
	virtual void SetDecKey(int key) { dec_key = key; }

	// trigger response to a key pressed
	virtual int RespondToKey(int key) 
	{
		if (key == inc_key) 
		{ 
			Increase(); 
			return 1;
		}
		if (key == dec_key) 
		{ 
			Decrease(); 
			return 1; 
		}
		return 0;
	}


private:
	int inc_key;								// key that triggers increase
	int dec_key;								// key that triggers decrease
}; //$EDC$[endclass,IOKeyResponseElement]

class IOMouseResponseElement: public IOResponseElement //$EDC$[beginclass,classname=IOMouseResponseElement,parentclassname=IOResponseElement,addelementtype=TAEinput_output+TAENotInRelease,addelementtypename=IOMouseResponseElement,
//texdescription="This element allows the IOBlock to respond to a mouse input. ( e.g. clicking the IOElement in the IOBlocksView) ",
//figure="IOMouseResponseElement,IOMouseResponseElement"]
{
public:
	IOMouseResponseElement(MBS* mbsi):IOResponseElement(mbsi)
	{	
		InitConstructor();
	};
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOMouseResponseElement(*this);
		return ec;
	}
		//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		IOResponseElement::CopyFrom(e);
		const IOMouseResponseElement& ce = (const IOMouseResponseElement&)e;
	}

	virtual void InitConstructor()
	{
		IOResponseElement::InitConstructor();
		input_mode = IOResponseElement::TMouseResponse;
	}
	virtual const char* GetElementSpec() const {return "IOMouseResponseElement";}

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

}; //$EDC$[endclass,IOMouseResponseElement]


class IOMinMax: public InputOutputElement//$EDC$[beginclass,classname=IOMinMax,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOMinMax,
//texdescription="This block returns the minimum, maximum or average value of the input. Up to a specific point of time, this functionality is switched off and the output y is equal to the input u. This block can be used to postprocess sensor values.",
//modus="{1 = minimum}{$y = u$ for $t \leq t_0$ \newline $y = \min_{t \geq t_0}(u)$ for $t > t_0 $ \newline with $t_0 = $ IOBlock.start\_time}",
//modus="{2 = maximum}{$y = u$ for $t \leq t_0$ \newline $y = \max_{t \geq t_0}(u)$ for $t > t_0 $}",
//modus="{3 = average}{$y = u$ for $t \leq t_0$ \newline $y = \frac{1}{N} \sum_{t_i \geq t_0} u_i $ for $t_i > t_0 $}",
//modus="{4 = minimum(abs)}{$y = u$ for $t \leq t_0$ \newline $y = \min_{t \geq t_0}(\left|u\right|)$ for $t > t_0 $}",
//modus="{5 = maximum(abs)}{$y = u$ for $t \leq t_0$ \newline $y = \min_{t \geq t_0}(\left|u\right|)$ for $t > t_0 $}",
//modus="{6 = average(abs)}{$y = u$ for $t \leq t_0$ \newline $y = \frac{1}{N} \sum_{t_i \geq t_0} \left|u_i\right| $ for $t_i > t_0 $}",
//example="IOMinMax.txt",
//figure="IOMinMax,IOMinMax"]
{
public:

	IOMinMax(MBS* mbsi):InputOutputElement(mbsi)
	{	
		InitConstructor();
		SetNStates(0);
		SetNOutputs(1);
	};

	void SetIOMinMax(int mode_i, double start_time_i = 0.)
	{
		start_time = start_time_i;
		mode = mode_i;		
	}
	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new IOMinMax(*this);
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		InputOutputElement::CopyFrom(e);
		const IOMinMax& ce = (const IOMinMax&)e;
		start_time = ce.start_time;
		mode = ce.mode;
		first_time_step = ce.first_time_step;
		current_value = ce.current_value;
	}

	virtual void InitConstructor()
	{
		InputOutputElement::InitConstructor();
		elementname = GetElementSpec();
		start_time = 0.;
		mode = 1;
		InitXData();
		first_time_step = 1;
	}

	// functions to store the min/max value in XData
	virtual int DataS() const { return 2; } // for initialization of XData vector ( remember the min/max value for all times )
	virtual void InitXData() 
	{ 
		Vector datainit(DataS());
		datainit(1) = 0.;	// just dummy value
		datainit(2) = 0.;	// start counter with 0
		SetDataInit(datainit);
	}
	virtual void WriteToXData()
	{ 
		XData(1) = current_value; 
		XData(2) = XData(2)+1;			// counter for average
	}

	virtual const char* GetElementSpec() const {return "IOMinMax";}
	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return 1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 
	virtual void GetElementData(ElementDataContainer& edc) 		//fill in all element data
	{
		InputOutputElement::GetElementData(edc);
	}
	virtual int SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
	{
		return InputOutputElement::SetElementData(edc);
	}

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual void StartTimeStep();
	virtual double GetOutput(double t, int i=1) const 
	{
		if(first_time_step) 
			return GetInput(t,1);
		else 	
		{
			if(mode == 3 || mode == 6)	// average
			{
				return current_value/XData(2);
			}
			else
				return current_value; 
		}
	}
	virtual int IsDirectFeedThrough() const {	return 1;	}

	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol() {	return InputOutputElement::DrawBlockSymbol(); }

protected:
	int mode;						//$EDC$[varaccess,EDCvarname="mode", minval=1,maxval=6,EDCfolder="IOBlock",tooltiptext="1..min, 2..max, 3..avg, 4..min(abs), 5..max(abs), 6..avg(abs)"]
	double start_time;  //$EDC$[varaccess,EDCvarname="start_time", EDCfolder="IOBlock",tooltiptext="Up to this point of time, the output is equal to the input. Afterwards the output is computed according to the mode."]
	int first_time_step;
	double current_value;
};//$EDC$[endclass,IOMinMax]


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: IOTCPIPBlock
// short description: Communication block for TCP/IP
// available formulations:
// type:
// development status: not finished
// long description: 
//**end(ued)**
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

enum CommunicationFlag {Fneutral = 0x00000000, Freset = 0x00000001, Ferror = 0x00000002, Fclose = 0x00000003};

//$ MSax 2013-5:[ added
class IOTCPIPBlock: public InputOutputElement //$EDC$[beginclass,classname=IOTCPIPBlock,parentclassname=InputOutputElement,addelementtype=TAEinput_output,addelementtypename=IOTCPIPBlock,
//texdescription="This I/O element is a communication block based on TCP/IP which allows HOTINT to connect to other programs or tools, opening up a large range of possible applications
// including external control, user-defined ``add-ons'', or even co-simulation.//
// Based on the specified IP (v4) address and port number the IOTCPIPBlock sets up a server socket and waits for a
// connection request from a client. Hence, HOTINT here plays the server role, and the external program is the client application. Data exchange is performed at a stage before every time step in HOTINT, following below protocol:\\ 
// The outgoing data, i.e.~the data sent from HOTINT to the client, is an array of 8-byte double precision numbers corresponding to the current values of the inputs of the I/O element, and additionally, an 8-byte sequence
// which is appended to that array and includes communication control flags (see the \textit{Communication flags} - section below for more details). Hence, the total amount of outgoing data is (number of inputs\,+\,1) times 8 bytes (double precision numbers).
// After the client has received and processed that data, it sends back a data package to HOTINT -- the incoming data for the I/O element --  which again consists of an array of double precision numbers, this time with the length (number of outputs\,+\,1).
// The first (number of output) double precision values determine the outputs of the I/O element, and the last 8 bytes again are used for the transfer of communication flags.\\
// HOTINT now begins the computation of one time step, where the transmitted data from the client is accessible via the outputs of the IOTCPIPBlock. \vspace*{12pt} \\
// \textit{Important notes} \vspace*{6pt} \\
// -- The waiting procedure for the client connection request, as well as the send and receive operations all are so-called ``blocking calls''. This means that HOTINT will wait for those operations to
// finish, and during that time, not respond to any user input. Therefore, a reasonable timeout (default is 10 seconds) should be specified for the IOTCPIPBlock to allow TCP/IP connection or transmission error handling. \vspace*{6pt}\\
// -- You will probably have to adjust your firewall settings and set appropriate permissions for HOTINT and the client application. \vspace*{6pt}\\
// --  Depending on the implementation of the client, it might be neccessary to start the server, i.e., HOTINT, first. \vspace*{6pt} \\
// -- Since HOTINT is running on Microsoft Windows, the memory byte order, also called ``endianness'', is ``Little Endian'', which means that the least significant bytes/digits are stored ``first'' in memory, i.e., on the smallest memory address.
// Therefore, any data sent from or received by the IOTCPIPBlock has or must have that byte order, respectively.
// You probably have to take that into account on the client side, especially if the client is running on a different platform and/or architecture on another computer. \vspace*{12pt}\\
// \textit{Communication flags} \vspace*{6pt} \\
// Currently, the following 4-byte flags are implemented: \vspace*{6pt}\\
// (1) Neutral flag: \texttt{0x00000000} (integer value: \texttt{0}). This flag signals that the application is running (properly) and no further action is required. \\
// (2) Reset flag: \texttt{0x00000001} (integer value: \texttt{1}). This flag is sent from HOTINT to the client in the first step of the computation. This can be used, for instance, to
// reset the client application.\\
// (3) Error flag: \texttt{0x00000002} (integer value: \texttt{2}). Indicates that an error has occurred. If HOTINT receives the error flag, an error message is issued, the connection is closed and the program execution terminated. \\
// (4) Close flag: \texttt{0x00000003} (integer value: \texttt{3}). This flag is sent from HOTINT to the client to indicate that the computation has finished and the connection will be closed, which is the case when the computation has actually finished, or the ``Stop''-button has been hit. \\
// (5) Any other value: Treated as error flag (3). \vspace*{6pt}\\  
// One of these flags is stored in and read from the last 8 bytes of the exchanged data -- corresponding to one additional double precision number -- in either direction in every time step. Currently, for simplicity, the flag is just casted explicitly from an integer to
// a double precision number which then can be transmitted and casted back to an integer exactly. Of course, this procedure must be followed on both the server and the client side.
//",
//texdescriptionComments="For the use of this element an active network adapter is required. If you only want to communicate locally on your computer and do not have an active network adapter,
// you can use a so-called ``loopback adapter'' which emulates an active real adapter in a real network and can be configured and used as such. The following steps summarize how 
//  a ``loopback adapter'' can be installed on Microsoft Windows: \vspace*{6pt} \\
//(1) Open the device manager\\
//(2) Select the network category and choose ``Action $\rightarrow$ Add legacy hardware'' via the menu\\
//(3) Choose the option for manual installation and select the category ``network adapters'' from the list\\
//(4) In the next dialog select ``Microsoft'' as vendor and ``Microsoft loopback adapter'' as hardware component\\
//(5) Proceed and finish the installation\\ \\
//\textbf{Example: Communication with Simulink/Matlab} \\
//This example demonstrates how to realize a connection between HOTINT and Matlab/Simulink. The purpose of the TCP/IP block is to use other powerful tools for some computations. For example it is possible to do the control law calculations for the actuation of the multibody system in Simulink (as alternative to the IOBlocks in HOTINT). It's also very simple to do a parameter variation, see the advanced example in the folder ``examples/balancing\_cart\_TCPIP''. \\
//In ``examples/TCPIP'' a very simple communication example is included: From the HOTINT side four different double values (simulation time t multiplied by the gain factors one to four) are transmitted to Simulink. Simulink summates the first and second respectively the third and fourth value and sends the two double values back to HOTINT. The values are captured by sensors, stored in the solution file and can be visualized in the plot tool. \\ 
//Comment: For testing purposes you can also use the executable ``TCPIP\_client.exe'' which has the same functionality as the Simulink example. To use this client executable create a ``IP.txt'' file in the same folder. The first four lines represent the IP address of the HOTINT computer, the fifth line is the port number. \\ \\
//To start this example, following things have to be done: \\
//(1) Start Matlab/Simulink and open the file \textsf{communication.mdl} in the folder ``examples/TCPIP'', see figure \ref{structure_Simulink}.\\
//Comment: If the ``Instrument Control Toolbox'' is not installed the TCP/IP communication blocks appear red and indicate an unresolved reference to a library block (bad link).  The figure shows the basic structure that should not be changed. The output of the ``TCP/IP Receive'' block is a vector $y_{rec}=[t,x_1,...,x_n,f]^{T}$ with HOTINT time $t$, data variables $x_1$ to $x_n$ and the handling flag $f$. The ``Selector'' block outputs the last element of the vector (flag $f$) for the flag handling. You have to adapt the these two blocks if you want to change the number of received variables. There is no need to change the ``flag handling in'' block.\\
//\begin{figure}[H]
//	\begin{center}
//		\includegraphics[width=0.9\textwidth]{D:/cpp/HotInt_V1/documentation/EDCauto_documentation/figures/TCPIP.png}
//		\caption{TCP/IP communication with Matlab/Simulink (do not change this structure)}
//		\label{structure_Simulink}
//	\end{center}
//\end{figure}
//\begin{figure}[H]
//	\begin{center}
//		\includegraphics[width=0.9\textwidth]{D:/cpp/HotInt_V1/documentation/EDCauto_documentation/figures/TCPIP_SS_computations.png}
//		\caption{Subsystem computations}
//		\label{SS_comp}
//	\end{center}
//\end{figure}
//(2) Make sure that the ``Current Folder'' is the folder which include the \textsf{communication.mdl} file.\\
//(3) Double click the ``TCP/IP Receive'' block and select the ``Remote address'' (i.e., the IP address) of the computer HOTINT is running on and select a ``Port''. Repeat this point for the ``TCP/IP Send'' block. \\
//Comment: If HOTINT and Simulink is running on the same computer you can also choose localhost (``127.0.0.1''). \\
//(4) Set the ``Sample Time'' of every block (TCP/IP Receive, Constants,...) and choose fixed step size in the ``Solver Options''. \\
//(5) Open the subsystem ``computations'', see figure \ref{SS_comp}. This subsystem contains all computations $\mathbf{y=f\left(u\right)}$ with input $\mathbf{u}$ and output $\mathbf{y}$. \\
//Comment: Change this subsystem to your needs. \\
//(6) Open the subsystem ``flag handling out'', see figure \ref{fho}. In default no handling flags are transmitted to HOTINT. \\
//Comment: Change this subsystem to your needs. \\
//(7) Save the mdl file. 
//\begin{figure}[H]
//	\begin{center}
//		\includegraphics[width=0.9\textwidth]{D:/cpp/HotInt_V1/documentation/EDCauto_documentation/figures/TCPIP_flagHandlingIn.png}
//		\caption{TCP/IP subsystem ``flag handling in''}
//	\end{center}
//\end{figure}
//\begin{figure}[H]
//	\begin{center}
//		\includegraphics[width=0.9\textwidth]{D:/cpp/HotInt_V1/documentation/EDCauto_documentation/figures/TCPIP_flagHandlingOut.png}
//		\caption{TCP/IP subsystem ``flag handling out''}
//		\label{fho}
//	\end{center}
//\end{figure}
//(8) Open the \textsf{TCPIP.hid} HOTINT file and type in the same ``ip\_address'' and ``port\_number'' as for the Matlab/Simulink side. \\
//(9) Make sure that ``max\_step\_size'' and ``min\_step\_size'' in the subtree ``SolverOptions.Timeint'' are set to the same value as the fixed ``Sample Time'' in Simulink. \\
//Comment: This is very important especially for the case of time dependent blocks like integrators in Simulink. \\
//(10) Save the file. \\
//(11) Load the \textsf{communication.hid} file in HOTINT. \\
//(12) Click the ``Start simulation'' button in Simulink. \\
//(13) Click the ``Start!'' button in HOTINT. \\ \\
//Comment: The points 11-13 have to be executed within the timeout limits. You can change the latter in the TCP/IP blocks for both HOTINT and Simulink. During these steps connection errors might occur due to firewall restrictions; you will probably have to set the corresponding permissions in your firewall(s). \\
//It is also recommended to choose the Simulink ``Simulation stop time'' higher as the ``end\_time'' in HOTINT. The reason is that HOTINT sends a stop flag after the last simulation step and in Simulink this flag is used to execute a ``Stop'' block which ends the communication and simulation.
//",example="TCPIP.txt"]
{
public:

	IOTCPIPBlock(MBS* mbsi):InputOutputElement(mbsi),serversocket(mbsi->GetServerSocket())
	{	
		InitConstructor();
	};

	//void SetIOTCPIPBlock(....) // TODO
	//{

	//};

	//To be overwritten in derived class:
	virtual Element* GetCopy();

	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e);

	virtual void InitConstructor();

	virtual const char* GetElementSpec() const {return "IOTCPIPBlock";}
	//virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual int GetExpectedNumberOfInputs(){return -1;}; // return the needed expected inputs, if not correct, the check of consistency reports error; -1 ... arbitrary number of inputs expected 

	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer

	//-----------------------------------------------------
	//b: auto-generated functions
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);	
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); // Write access to a single element variable
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); // returns all available directly (ReadWrite-) accessable variables
	//e: auto-generated functions
	//-----------------------------------------------------

	virtual void DrawElement() 
	{
		InputOutputElement::DrawElement();
	};
	virtual const char* SymbolText() const;
	virtual void DrawBlockSymbol() {	return InputOutputElement::DrawBlockSymbol(); }

	virtual double GetOutput(double t, int i=1) const
	{
		return XData(i);
	}

	void SetOutput(double val, int i=1) //write-access to XData
	{
		XData(i)=val;
	}

	virtual void StartTimeStep() // communication -> receive and send data
	{
		//double t_end = GetMBS()->GetStepEndTime(); //t_end = t+deltaT
		double t = GetMBS()->GetTime();

		//int i = 1;
		//double u = GetInput(t, i);
		// TODO

		Communicate(t);

		//double u1 = GetInput(t, 1);
		//double u2 = GetInput(t, 2);
		//for (int i=1; i<=n_output; i++)
		//{
		//	if (i<3) XData(i) = u1+i;
		//	else XData(i) = u2+i;
		//}
	}

	virtual void Initialize(); // initialization TCPIP connection

	virtual int DataS() const { return n_output; } // for initialization of XData vector ( remember the current_value for all times )

	virtual void ComputationFinished();

  virtual ~IOTCPIPBlock();

	//set a CommunicationFlag on byte position pos (starting from 0) in the array addtcpdata (i.e., at the address addtcpdata+i)
	//if convertdouble is set, the flag is converted explicitly to double before
	void SetCommunicationFlag(int pos, CommunicationFlag flag, int convertdouble=0);

	//get a CommunicationFlag on byte position pos in addtcpdata
	//if convertdouble is set, it is assumed that the flag has been converted to a an 8-byte double-precision number before
	CommunicationFlag GetCommunicationFlag(int pos, int convertdouble=0);

protected:
	//variables
	int port_number; //$EDC$[varaccess,EDCvarname="port_number", EDCfolder="IOBlock",tooltiptext="Port number, e.g. '50000'."]
	mystr ip_address; //$EDC$[varaccess,EDCvarname="ip_address", EDCfolder="IOBlock",tooltiptext="IP address, e.g. '127.0.0.1' (localhost). Do not neglect the dots between the numbers."]
	//EDC int n_output;   //$EDC$[varaccess,EDCvarname="received_data_size", EDCfolder="IOBlock",tooltiptext="Number of received data values (outputs). This number has to be consistent with the transmitted data values of the other communication side (the additional double for the communication flags is not corresponding to this number)."]

	int timeout; //$EDC$[varaccess,EDCvarname="timeout", EDCfolder="IOBlock",tooltiptext="TCP/IP timeout in milliseconds; default is 10000."]

	TCPIPHotInt* & serversocket;
	double* tcpdata;  // auxiliary array for TCP/IP data transfer (data must be contiguous in memory)
	char* addtcpdata; // auxiliary array for additional TCP/IP data 
	
	int use_default_protocol; //enabled by default; if so, the following flags and parameters are set to their defaults and the default protocol is used (for more details see element description)	
	int use_additional_communication;  // flag to enable/disable transfer of additional data (e.g. protocol information); 1 by default
	int separate_additional_communication; // if this flag is set, first, the additional information is transferred, and then the actual data separately; otherwise, the additional info is just appended to the data; 0 by default
	int additional_communication_length; // number of data units of additional data (same for incoming and outgoing data); zero, if use_additional_communication is not set; 1 by default
	int additional_communication_unit_size; // size of data units of additional data; possible values: 1,2,4,8 bytes; default value is 8; if separate_additional_communication is not set, the default is used

	int isinitialized; //internally used for this element; status of the TCPIP socket via serversocket.GetStatus()
	void Communicate(double t);
	void OutgoingDataCommunication(double t);
	void IncomingDataCommunication();
	int ReactToAdditionalCommunication(double t); // for support of any communication scheme / protocol

	//close connection, memory clean-up
	void CloseAndCleanUp();

}; //$EDC$[endclass,IOTCPIPBlock]
//$ MSax 2013-5:]

#endif
