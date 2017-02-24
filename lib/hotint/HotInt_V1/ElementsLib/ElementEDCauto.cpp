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


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+  automatically generated file for EDC class data exchange, do not modify!   +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//DESCRIPTION OF $EDC$ usage:                                                    
//* mark beginning of class: $EDC$[beginclass,classname='class_name',parentclassname='parent_class_name',addelementtypename='EDCelementname',texdescription="...",]     
//* mark end of class: $EDC$[endclass,'classname']                               
//* with 'addelementtypename', the class is made available for adding on file load or in the menu of HOTINT
//* EDC folders start with upper case letter, variables with lower case letter!!!
//                                                                               
//variable access:                                                               
//  + syntax: [//EDC|const|static|mutable] [typename] [C++ variable name] //$EDC$[varaccess, {option1, {option2}}]
//             ==>//EDC means that this variable is commented out totally indicated by symbol '//EDC' ==> e.g. because it is defined in another class
//  + this generates an automated (EDC) access to this variable                  
//* possible options:                                                            
//  + internal_varaccess: create internal access functions for use in sensors and/or modifiers 
//  + readonly: this variable can be only read (e.g. for debugging, but not modified 
//  + EDCvarname="variable_name": assign variable name, if EDC-variable name should be different from C++ variable name
//  + EDCfolder="folder_name": assign subfolder name of EDC, e.g. "Graphics" or "Physics"
//  + tooltiptext="text"]: tooltip (help) text, which describes the variable; this could also be used to document the variable name in HOTINT
//  + int_bool: treat integer value as boolean (with check box)
//  + variable_length_vector: vector has variable length ==> user can choose how many values are set
//  + func_par1 = ...: assign a number or evaluable expression which is used as function parameter for double or int access
//  + condition = ...: evaluable C++ expression, e.g. 'IsConstraint()' or 'i < 4', which is used as condition if Get/SetElementData entry is used or not
//  + maxval = ..., minval = ... : define maximum and minimum integer values which are allowed
//  + vecstart = ..., vecend = ... : define starting and ending of a vector by means of values or C++ expressions, if only a sub-vector should be used
//  + remove: remove variable from an EDC of the parent class in a derived class. The elementdata is not really removed but locked and given a new tooltiptext
//  
//  + e.g.: static int switch; //$EDC$[varaccess, EDCvarname="switch_for_HOTINT", EDCfolder="Debug", tooltiptext="this switch is nonsense", readonly]
//  
//function access:                                                               
//  + syntax: [const|static|mutable|virtual] [typename] [C++ function name]() //$EDC[funcaccess, {option1, {option2}}]
//  + same options as for variable
//  
//  
#include "..\MBS_Interface\element.h"
#include "..\UtilityLib\elementdataaccess.h"
#include "..\MBS_Interface\node.h"
#include "..\MBS_Interface\material.h"
#include "..\MBSElementsAndModels\MBSObjectFactory.h"
#include "..\MBS_Interface\mbs_interface.h"
#include "..\ElementsLib\rigid\mass1D.h"
#include "..\ElementsLib\body2D.h "
#include "..\ElementsLib\rigid\rigid2D.h"
#include "..\ElementsLib\body3D.h"
#include "..\ElementsLib\rigid\rigid3D.h"
#include "..\ElementsLib\rigid\Rigid3DKardan.h"
#include "..\ElementsLib\beams\Beam3D.h"
#include "..\ElementsLib\beams\ANCFCable2D.h"
#include "..\ElementsLib\beams\ANCFCable3D.h"
#include "..\ElementsLib\beams\ANCFBeamShear3D.h"
#include "..\ElementsLib\beams\ANCFBeam3DTorsion.h"
#include "..\ElementsLib\beams\FiniteElementGenericBeam2D.h"
#include "..\ElementsLib\beams\ANCFBeamShearFE2D.h"
#include "..\ElementsLib\shells\ANCFSimpleThinPlate3D.h"
#include "..\ElementsLib\constraints\constraint.h"
#include "..\ElementsLib\constraints\KinematicPairs.h"
#include "..\ElementsLib\constraints\SpecialConstraints.h"
#include "..\ElementsLib\constraints\RigidBodyJoints.h"
#include "..\ElementsLib\constraints\AverageConstraint.h"
#include "..\MBS_Interface\sensors.h"
#include "..\ServiceObjectsLib\sensorsSpecific.h"
#include "..\ElementsLib\constraints\control.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Element

void Element::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  ed.SetText(GetElementSpecification().c_str(),"element_type"); ed.SetToolTipText("specification of element type. Once the element is added to the mbs, you MUST NOT change this type anymore!"); edc.TreeAdd("",ed);
  ed.SetInt(IsRigid(),"is_rigid"); ed.SetLocked(1); ed.SetToolTipText("flag if body is a rigid body"); edc.TreeAdd("Info",ed);
  ed.SetInt(IsFiniteElement(),"is_finite_element"); ed.SetLocked(1); ed.SetToolTipText("flag if body is a finite element"); edc.TreeAdd("Info",ed);
  ed.SetInt(IS(),"n_algebraic_equations"); ed.SetLocked(1); ed.SetToolTipText("number of algebraic equations"); edc.TreeAdd("Info",ed);
  ed.SetInt(SOS(),"n_second_order_ODEs"); ed.SetLocked(1); ed.SetToolTipText("number of second order ODEs in which the element contributes with some terms to the mass matrix and RHS"); edc.TreeAdd("Info",ed);
  ed.SetInt(SOSowned(),"n_second_order_ODE_unknowns"); ed.SetLocked(1); ed.SetToolTipText("number of second order ODE unknowns, which are owned (caused) by one element (not borrowed e.g. from another element or node"); edc.TreeAdd("Info",ed);
  ed.SetInt(ES(),"n_first_order_ODEs"); ed.SetLocked(1); ed.SetToolTipText("number of first order ODEs"); edc.TreeAdd("Info",ed);
  ed.SetInt(SS(),"element_equation_size"); ed.SetLocked(1); ed.SetToolTipText("element system size=2*SOS()+ES()+IS()"); edc.TreeAdd("Info",ed);
  ed.SetInt(DataS(),"data_size"); ed.SetLocked(1); ed.SetToolTipText("number of data variables (e.g. plastic strains), for which there are no separate equations"); edc.TreeAdd("Info",ed);
  ed.SetInt(Dim(),"element_dimension"); ed.SetLocked(1); ed.SetToolTipText("dimensionality of element (2D/3D)"); edc.TreeAdd("Info",ed);
  ed.SetInt(NNodes(),"number_of_nodes"); ed.SetLocked(1); ed.SetToolTipText("number of nodes in finite elements"); edc.TreeAdd("Info",ed);
  ed.SetInt(NLoads(),"number_of_loads"); ed.SetLocked(1); edc.TreeAdd("Info",ed);
  ed.SetText(elementname.c_str(),"name"); ed.SetToolTipText("name of the element"); edc.TreeAdd("",ed);
  ed.SetInt(elnum,"element_number"); ed.SetLocked(1); ed.SetToolTipText("number of the element in the mbs"); edc.TreeAdd("",ed);
  if(!IsType(TConstraint))
  {
  {Vector vv(loads.Length()); for(int i = 1; i <= loads.Length(); i++) {vv(i) = loads(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"loads"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("Set loads attached to this element: 'nr_load1, nr_load2, ...' or empty"); edc.TreeAdd("",ed);
}
  }

  {Vector vv(sensors.Length()); for(int i = 1; i <= sensors.Length(); i++) {vv(i) = sensors(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"sensors"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetVariableLength(); ed.SetToolTipText("attached sensors"); edc.TreeAdd("Info",ed);
}
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"RGB_color"); ed.SetToolTipText("[red, green, blue] color of element, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  if(!IsType(TConstraint))
  {
  {Vector vv(drawelements.Length()); for(int i = 1; i <= drawelements.Length(); i++) {vv(i) = drawelements(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"geom_elements"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("Set Geometric elements to represent body 'geomelem1, geomelem2, ...' or empty"); edc.TreeAdd("Graphics",ed);
}
  }
  if(!IsType(TConstraint))
  {  ed.SetBool(altshape,"use_alternative_shape"); ed.SetToolTipText("Graphical representation of element with geom-objects that are attached to the element"); edc.TreeAdd("Graphics",ed);
  }
  ed.SetBool(draw_element,"show_element"); ed.SetToolTipText("Flag to draw element"); edc.TreeAdd("Graphics",ed);
}


int Element::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GetElemDataText(GetMBS(), edc, "element_type",GetElementSpecification(), 1);
  GetElemDataText(GetMBS(), edc, "name",elementname, 1);
  if(!IsType(TConstraint))
  {  GetElemDataIVector(GetMBS(), edc, "loads",loads, 1);
  }
  GetElemDataVector3D(GetMBS(), edc, "Graphics.RGB_color",col, 1);
  if(!IsType(TConstraint))
  {  GetElemDataIVector(GetMBS(), edc, "Graphics.geom_elements",drawelements, 1);
  }
  if(!IsType(TConstraint))
  {  GetElemDataBool(GetMBS(), edc, "Graphics.use_alternative_shape",altshape, 1);
  }
  GetElemDataBool(GetMBS(), edc, "Graphics.show_element",draw_element, 1);
  return 1;
}


int Element::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  return 0;
}


int Element::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  return 0;
}


int Element::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=BaseNode

void BaseNode::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  ed.SetText(GetTypeName().c_str(),"node_type"); ed.SetToolTipText("specification of node type. Once the node is added to the mbs, you MUST NOT change this type anymore!"); edc.TreeAdd("",ed);
  ed.SetText(nodename.c_str(),"name"); ed.SetToolTipText("Node identifier."); edc.TreeAdd("",ed);
  ed.SetInt(nodenum,"node_number"); ed.SetToolTipText("Node Number."); edc.TreeAdd("",ed);
  ed.SetVector3D(pos.X(),pos.Y(),pos.Z(),"reference_position"); ed.SetToolTipText("Position (2D/3D) in reference configuration."); edc.TreeAdd("Geometry",ed);

  {Vector vv((2*SOS())-(1)+1); for (int i=1; i<=2*SOS(); i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_initial_values"); ed.SetToolTipText("initial values for all degrees of freedom of node"); edc.TreeAdd("Initialization",ed);}
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"RGB_color"); ed.SetToolTipText("[red, green, blue] color of element, range = 0..1,"); edc.TreeAdd("Graphics",ed);
  ed.SetInt(visible,"visible"); ed.SetToolTipText("Visibility of node."); edc.TreeAdd("Graphics",ed);
}


int BaseNode::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GetElemDataText(GetMBS(), edc, "node_type",GetTypeName(), 1);
  GetElemDataText(GetMBS(), edc, "name",nodename, 1);
  GetElemDataInt(GetMBS(), edc, "node_number",nodenum, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_position",pos, 1);
  {Vector vv((2*SOS())-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node_initial_values",vv, 1);
  for (int i=1; i<=2*SOS(); i++) {x_init(i)=vv(i+1-(1));}
}
  GetElemDataVector3D(GetMBS(), edc, "Graphics.RGB_color",col, 1);
  GetElemDataInt(GetMBS(), edc, "Graphics.visible",visible, 1);
  return 1;
}


int BaseNode::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  return 0;
}


int BaseNode::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  return 0;
}


int BaseNode::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Node

void Node::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseNode::GetElementDataAuto(edc);
  ElementData ed;
}


int Node::SetElementDataAuto(ElementDataContainer& edc)
{
  BaseNode::SetElementDataAuto(edc);
  return 1;
}


int Node::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseNode::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseNode::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=GenericNode

void GenericNode::GetElementDataAuto(ElementDataContainer& edc)
{
  Node::GetElementDataAuto(edc);
  ElementData ed;
}


int GenericNode::SetElementDataAuto(ElementDataContainer& edc)
{
  Node::SetElementDataAuto(edc);
  return 1;
}


int GenericNode::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int GenericNode::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int GenericNode::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Node2D

void Node2D::GetElementDataAuto(ElementDataContainer& edc)
{
  GenericNode::GetElementDataAuto(edc);
  ElementData ed;
}


int Node2D::SetElementDataAuto(ElementDataContainer& edc)
{
  GenericNode::SetElementDataAuto(edc);
  return 1;
}


int Node2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = GenericNode::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = GenericNode::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Node3D

void Node3D::GetElementDataAuto(ElementDataContainer& edc)
{
  GenericNode::GetElementDataAuto(edc);
  ElementData ed;
}


int Node3D::SetElementDataAuto(ElementDataContainer& edc)
{
  GenericNode::SetElementDataAuto(edc);
  return 1;
}


int Node3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = GenericNode::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = GenericNode::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Node3DRxyz

void Node3DRxyz::GetElementDataAuto(ElementDataContainer& edc)
{
  Node3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(ref_angles.X(),ref_angles.Y(),ref_angles.Z(),"reference_rot_angles"); ed.SetToolTipText("Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int Node3DRxyz::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Node3D::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_rot_angles",ref_angles, 1);
  return 1;
}


int Node3DRxyz::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3DRxyz::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3DRxyz::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Node3DR123

void Node3DR123::GetElementDataAuto(ElementDataContainer& edc)
{
  Node3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(ref_angles.X(),ref_angles.Y(),ref_angles.Z(),"reference_rot_angles"); ed.SetToolTipText("Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int Node3DR123::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Node3D::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_rot_angles",ref_angles, 1);
  return 1;
}


int Node3DR123::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3DR123::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Node3DR123::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFNodeS1rot1_3D

void ANCFNodeS1rot1_3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Node3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(ref_angles.X(),ref_angles.Y(),ref_angles.Z(),"reference_rot_angles"); ed.SetToolTipText("Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int ANCFNodeS1rot1_3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Node3D::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_rot_angles",ref_angles, 1);
  return 1;
}


int ANCFNodeS1rot1_3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS1rot1_3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS1rot1_3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFNodeS2S3_3D

void ANCFNodeS2S3_3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Node3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(ref_slope2.X(),ref_slope2.Y(),ref_slope2.Z(),"reference_slope2"); ed.SetToolTipText("slope 2 of node in reference configuration."); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(ref_slope3.X(),ref_slope3.Y(),ref_slope3.Z(),"reference_slope3"); ed.SetToolTipText("slope 3 of node in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int ANCFNodeS2S3_3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Node3D::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_slope2",ref_slope2, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_slope3",ref_slope3, 1);
  return 1;
}


int ANCFNodeS2S3_3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS2S3_3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS2S3_3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFNodeS1S2_3D

void ANCFNodeS1S2_3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Node3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(ref_slope1.X(),ref_slope1.Y(),ref_slope1.Z(),"reference_slope1"); ed.SetToolTipText("slope 1 of node in reference configuration."); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(ref_slope2.X(),ref_slope2.Y(),ref_slope2.Z(),"reference_slope2"); ed.SetToolTipText("slope 2 of node in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int ANCFNodeS1S2_3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Node3D::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_slope1",ref_slope1, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_slope2",ref_slope2, 1);
  return 1;
}


int ANCFNodeS1S2_3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Node3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS1S2_3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Node3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS1S2_3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFNodeS2_2D

void ANCFNodeS2_2D::GetElementDataAuto(ElementDataContainer& edc)
{
  GenericNode::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector2D(ref_slope2.X(),ref_slope2.Y(),"reference_slope2"); ed.SetToolTipText("Slope of cross section (direction eta) in reference configuration."); edc.TreeAdd("Geometry",ed);
}


int ANCFNodeS2_2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GenericNode::SetElementDataAuto(edc);
  GetElemDataVector2D(GetMBS(), edc, "Geometry.reference_slope2",ref_slope2, 1);
  return 1;
}


int ANCFNodeS2_2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = GenericNode::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS2_2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = GenericNode::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFNodeS2_2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Material

void Material::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  ed.SetText(GetTypeName().c_str(),"material_type"); ed.SetToolTipText("specification of material type. Once the material is added to the mbs, you MUST NOT change this type anymore!"); edc.TreeAdd("",ed);
  ed.SetVector3D(materialcolor.X(),materialcolor.Y(),materialcolor.Z(),"color"); ed.SetToolTipText("material color (as yet used with FEMesh, only)"); edc.TreeAdd("Graphics",ed);
  ed.SetText(materialname.c_str(),"name"); ed.SetToolTipText("name of the material"); edc.TreeAdd("",ed);
  ed.SetDouble(density,"density"); ed.SetToolTipText("density (rho) for gravitational force"); edc.TreeAdd("Solid",ed);
  ed.SetDouble(youngs_modulus,"youngs_modulus"); ed.SetToolTipText("Youngs modulus"); edc.TreeAdd("Solid",ed);
  ed.SetDouble(poisson_ratio,"poisson_ratio"); ed.SetToolTipText("Poisson ratio"); edc.TreeAdd("Solid",ed);
  ed.SetBool(plane,"plane"); ed.SetToolTipText("true: 2D, false: 3D"); edc.TreeAdd("Solid",ed);
  ed.SetBool(plane_stress,"plane_stress"); ed.SetToolTipText("for 2D-Elements only; 1: plane stress, 0: plane strain"); edc.TreeAdd("Solid",ed);
  ed.SetDouble(yieldstress,"yield_stress"); ed.SetToolTipText("Yield Stress s_y, e.g., |dev s| <= s_y"); edc.TreeAdd("Inelasticity",ed);
  ed.SetDouble(tangentmodule,"tangent_module"); ed.SetToolTipText("Modulus of hardening H"); edc.TreeAdd("Inelasticity",ed);
  ed.SetText(inelasticity_type_str.c_str(),"inelasticity_type"); ed.SetToolTipText("linear_elastic, elasto_plastic (= Prandtl Reuss plasticity + isotropic hardening), nonlinear_elastic_Simo_Hughes (see Simo and Hughes, Computational Inelasticity 1998: S=lambda/2*(J*J-1)/C + mu*(1-1/C))"); edc.TreeAdd("Inelasticity",ed);
  ed.SetText(inelasticity_solution_method_str.c_str(),"inelasticity_solution_method"); ed.SetToolTipText("fixed_point, return_mapping, consistent_tangent_stiffness (see Simo and Hughes, Computational Inelasticity 1998)"); edc.TreeAdd("Inelasticity",ed);
}


int Material::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GetElemDataText(GetMBS(), edc, "material_type",GetTypeName(), 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color",materialcolor, 1);
  GetElemDataText(GetMBS(), edc, "name",materialname, 1);
  GetElemDataDouble(GetMBS(), edc, "Solid.density",density, 1);
  GetElemDataDouble(GetMBS(), edc, "Solid.youngs_modulus",youngs_modulus, 1);
  GetElemDataDouble(GetMBS(), edc, "Solid.poisson_ratio",poisson_ratio, 1);
  GetElemDataBool(GetMBS(), edc, "Solid.plane",plane, 1);
  GetElemDataBool(GetMBS(), edc, "Solid.plane_stress",plane_stress, 1);
  GetElemDataDouble(GetMBS(), edc, "Inelasticity.yield_stress",yieldstress, 1);
  GetElemDataDouble(GetMBS(), edc, "Inelasticity.tangent_module",tangentmodule, 1);
  GetElemDataText(GetMBS(), edc, "Inelasticity.inelasticity_type",inelasticity_type_str, 1);
  GetElemDataText(GetMBS(), edc, "Inelasticity.inelasticity_solution_method",inelasticity_solution_method_str, 1);
  return 1;
}


int Material::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  return 0;
}


int Material::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  return 0;
}


int Material::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Beam3DProperties

void Beam3DProperties::GetElementDataAuto(ElementDataContainer& edc)
{
  Material::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(beamCrossSectionType,"cross_section_type"); ed.SetToolTipText("1: rectangular, 2: circular, 3: polygonal"); edc.TreeAdd("",ed);
  ed.SetVector(beamCrossSectionSize.GetVecPtr(),beamCrossSectionSize.Length(),"cross_section_size"); ed.SetVariableLength(); ed.SetToolTipText("vector length of cross_section_size depends on cross_section_type: length 1 for circular cross section, length 2 for rectangular cross section (y and z extension), and length 2*n for polygonal cross section (p1y,p1z,p2y,p2z,...,pny,pnz)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamEA,"EA"); ed.SetToolTipText("youngs modulus * area"); edc.TreeAdd("",ed);
  ed.SetDouble(beamEIy,"EIy"); ed.SetToolTipText("bending stiffness w.r.t. y-axis (2D-beam)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamEIz,"EIz"); ed.SetToolTipText("bending stiffness w.r.t. z-axis"); edc.TreeAdd("",ed);
  ed.SetDouble(beamGAky,"GAky"); ed.SetToolTipText("shear stiffness including shear correction factor ky (2D-beam)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamGAkz,"GAkz"); ed.SetToolTipText("shear stiffness including shear correction factor kz"); edc.TreeAdd("",ed);
  ed.SetDouble(beamGJkx,"GJkx"); ed.SetToolTipText("torsional stiffness including shear correction factor kx"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoA,"RhoA"); ed.SetToolTipText("density * area"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoIx,"RhoIx"); ed.SetToolTipText("density * second area of moment w.r.t. x-axis"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoIy,"RhoIy"); ed.SetToolTipText("density * second area of moment w.r.t. y-axis (2D-beam)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoIz,"RhoIz"); ed.SetToolTipText("density * second area of moment w.r.t. z-axis"); edc.TreeAdd("",ed);
  edc.TreeDelete("Solid.density"); 
  edc.TreeDelete("Solid.youngs_modulus"); 
  edc.TreeDelete("Solid.poisson_ratio"); 
  edc.TreeDelete("Solid.plane_stress"); 
  edc.TreeDelete("Inelasticity.yield_stress"); 
  edc.TreeDelete("Inelasticity.tangent_module"); 
  edc.TreeDelete("Solid.plane"); 
  ed.SetDouble(density,"density"); ed.SetToolTipText("density (rho) for gravitational force"); edc.TreeAdd("",ed);
}


int Beam3DProperties::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"density"); ed.SetLocked(1); ed.SetToolTipText("density (rho) for gravitational force"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"youngs_modulus"); ed.SetLocked(1); ed.SetToolTipText("Youngs modulus"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"poisson_ratio"); ed.SetLocked(1); ed.SetToolTipText("Poisson ratio"); edc.TreeAdd("Solid",ed);}
  {  int dummy=0; ed.SetBool(dummy,"plane_stress"); ed.SetLocked(1); ed.SetToolTipText("for 2D-Elements only; 1: plane stress, 0: plane strain"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"yield_stress"); ed.SetLocked(1); ed.SetToolTipText("Yield Stress s_y, e.g., |dev s| <= s_y"); edc.TreeAdd("Inelasticity",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"tangent_module"); ed.SetLocked(1); ed.SetToolTipText("Modulus of hardening H"); edc.TreeAdd("Inelasticity",ed);}
  {  int dummy=0; ed.SetInt(dummy,"plane"); ed.SetLocked(1); edc.TreeAdd("Solid",ed);}


  Material::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "cross_section_type",beamCrossSectionType, 1);
  GetElemDataVector(GetMBS(), edc, "cross_section_size",beamCrossSectionSize, 1);
  GetElemDataDouble(GetMBS(), edc, "EA",beamEA, 1);
  GetElemDataDouble(GetMBS(), edc, "EIy",beamEIy, 1);
  GetElemDataDouble(GetMBS(), edc, "EIz",beamEIz, 1);
  GetElemDataDouble(GetMBS(), edc, "GAky",beamGAky, 1);
  GetElemDataDouble(GetMBS(), edc, "GAkz",beamGAkz, 1);
  GetElemDataDouble(GetMBS(), edc, "GJkx",beamGJkx, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoA",beamRhoA, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoIx",beamRhoIx, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoIy",beamRhoIy, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoIz",beamRhoIz, 1);
  GetElemDataDouble(GetMBS(), edc, "density",density, 1);
  return 1;
}


int Beam3DProperties::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Material::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam3DProperties::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Material::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam3DProperties::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Beam2DProperties

void Beam2DProperties::GetElementDataAuto(ElementDataContainer& edc)
{
  Material::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(beamCrossSectionType,"cross_section_type"); ed.SetToolTipText("1: rectangular, 2: circular, 3: polygonal"); edc.TreeAdd("",ed);
  ed.SetVector(beamCrossSectionSize.GetVecPtr(),beamCrossSectionSize.Length(),"cross_section_size"); ed.SetVariableLength(); ed.SetToolTipText("vector length of cross_section_size depends on cross_section_type: length 1 for circular cross section, length 2 for rectangular cross section (y and z extension), and length 2*n for polygonal cross section (p1y,p1z,p2y,p2z,...,pny,pnz)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamEA,"EA"); ed.SetToolTipText("youngs modulus * area"); edc.TreeAdd("",ed);
  ed.SetDouble(beamEIy,"EIy"); ed.SetToolTipText("bending stiffness w.r.t. y-axis (2D-beam)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamGAky,"GAky"); ed.SetToolTipText("shear stiffness including shear correction factor ky (2D-beam)"); edc.TreeAdd("",ed);
  ed.SetDouble(beamGJkx,"GJkx"); ed.SetToolTipText("torsional stiffness including shear correction factor kx"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoA,"RhoA"); ed.SetToolTipText("density * area"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoIx,"RhoIx"); ed.SetToolTipText("density * second area of moment w.r.t. x-axis"); edc.TreeAdd("",ed);
  ed.SetDouble(beamRhoIy,"RhoIy"); ed.SetToolTipText("density * second area of moment w.r.t. y-axis (2D-beam)"); edc.TreeAdd("",ed);
  edc.TreeDelete("Solid.density"); 
  edc.TreeDelete("Solid.youngs_modulus"); 
  edc.TreeDelete("Solid.poisson_ratio"); 
  edc.TreeDelete("Solid.plane_stress"); 
  edc.TreeDelete("Inelasticity.yield_stress"); 
  edc.TreeDelete("Inelasticity.tangent_module"); 
  edc.TreeDelete("Inelasticity.inelasticity_type"); 
  edc.TreeDelete("Solid.plane"); 
}


int Beam2DProperties::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"density"); ed.SetLocked(1); ed.SetToolTipText("density (rho) for gravitational force"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"youngs_modulus"); ed.SetLocked(1); ed.SetToolTipText("Youngs modulus"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"poisson_ratio"); ed.SetLocked(1); ed.SetToolTipText("Poisson ratio"); edc.TreeAdd("Solid",ed);}
  {  int dummy=0; ed.SetBool(dummy,"plane_stress"); ed.SetLocked(1); ed.SetToolTipText("for 2D-Elements only; 1: plane stress, 0: plane strain"); edc.TreeAdd("Solid",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"yield_stress"); ed.SetLocked(1); ed.SetToolTipText("Yield Stress s_y, e.g., |dev s| <= s_y"); edc.TreeAdd("Inelasticity",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"tangent_module"); ed.SetLocked(1); ed.SetToolTipText("Modulus of hardening H"); edc.TreeAdd("Inelasticity",ed);}
  {  int dummy=0; ed.SetInt(dummy,"inelasticity_type"); ed.SetLocked(1); ed.SetToolTipText("0 = not specified, 1 = Prandtl Reuss plasticity + isotropic hardening"); edc.TreeAdd("Inelasticity",ed);}
  {  int dummy=0; ed.SetInt(dummy,"plane"); ed.SetLocked(1); edc.TreeAdd("Solid",ed);}


  Material::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "cross_section_type",beamCrossSectionType, 1);
  GetElemDataVector(GetMBS(), edc, "cross_section_size",beamCrossSectionSize, 1);
  GetElemDataDouble(GetMBS(), edc, "EA",beamEA, 1);
  GetElemDataDouble(GetMBS(), edc, "EIy",beamEIy, 1);
  GetElemDataDouble(GetMBS(), edc, "GAky",beamGAky, 1);
  GetElemDataDouble(GetMBS(), edc, "GJkx",beamGJkx, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoA",beamRhoA, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoIx",beamRhoIx, 1);
  GetElemDataDouble(GetMBS(), edc, "RhoIy",beamRhoIy, 1);
  return 1;
}


int Beam2DProperties::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Material::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam2DProperties::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Material::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam2DProperties::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Mass1D

void Mass1D::GetElementDataAuto(ElementDataContainer& edc)
{
  Element::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(drawres,"drawing_tiling"); ed.SetToolTipText("tiling of circle/sphere to represent Mass1D;	the drawing_tiling should be set small in order to improve efficiency,	but large for nice graphical represenations"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(radius,"radius"); ed.SetToolTipText("drawing radius of mass"); edc.TreeAdd("Graphics",ed);

  {Vector vv((1)-(1)+1); for (int i=1; i<=1; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_position"); ed.SetToolTipText("initial values for position [x]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((2)-(2)+1); for (int i=2; i<=2; i++) {vv(i+1-(2))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_velocity"); ed.SetToolTipText("initial values for velocity [v]"); edc.TreeAdd("Initialization",ed);}
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("total mass of point mass"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(pref3D.X(),pref3D.Y(),pref3D.Z(),"reference_position"); ed.SetToolTipText("Reference point for transformation of 1D objects to 3D; p = [X, Y, Z]"); edc.TreeAdd("Graphics",ed);

  Matrix rotref3D_1(rotref3D);
  ed.SetMatrix(rotref3D_1.GetMatPtr(),rotref3D_1.Getrows(),rotref3D_1.Getcols(),"rotation_matrix"); ed.SetToolTipText("Rotation matrix for transformation of 1D objects to 3D"); edc.TreeAdd("Graphics",ed);

}


int Mass1D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Element::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Graphics.drawing_tiling",drawres, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.radius",radius, 1);
  {Vector vv((1)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_position",vv, 1);
  for (int i=1; i<=1; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((2)-(2)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_velocity",vv, 1);
  for (int i=2; i<=2; i++) {x_init(i)=vv(i+1-(2));}
}
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.reference_position",pref3D, 1);

  Matrix rotref3D_1;
  GetElemDataMatrix(GetMBS(), edc, "Graphics.rotation_matrix",rotref3D_1, 1);
  rotref3D = Matrix3D(rotref3D_1);
  return 1;
}


int Mass1D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Element::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass1D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Element::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass1D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Rotor1D

void Rotor1D::GetElementDataAuto(ElementDataContainer& edc)
{
  Mass1D::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Initialization.initial_position"); 
  edc.TreeDelete("Initialization.initial_velocity"); 

  {Vector vv((1)-(1)+1); for (int i=1; i<=1; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_rotation"); ed.SetToolTipText("initial value for rotation"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((2)-(2)+1); for (int i=2; i<=2; i++) {vv(i+1-(2))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_angular_velocity"); ed.SetToolTipText("initial value for angular velocity"); edc.TreeAdd("Initialization",ed);}
  edc.TreeDelete("Physics.mass"); 
  ed.SetDouble(mass,"moment_of_inertia"); ed.SetToolTipText("mass moment of inertia in kg*m*m"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Graphics.drawing_tiling"); 
  edc.TreeDelete("Graphics.radius"); 
  ed.SetDouble(radius,"radius"); ed.SetToolTipText("radius of rotor"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(length,"length"); ed.SetToolTipText("length of rotor"); edc.TreeAdd("Graphics",ed);
}


int Rotor1D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"initial_position"); ed.SetLocked(1); edc.TreeAdd("Initialization",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"initial_velocity"); ed.SetLocked(1); edc.TreeAdd("Initialization",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"mass"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"drawing_tiling"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"radius"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  Mass1D::SetElementDataAuto(edc);
  {Vector vv((1)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_rotation",vv, 1);
  for (int i=1; i<=1; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((2)-(2)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_angular_velocity",vv, 1);
  for (int i=2; i<=2; i++) {x_init(i)=vv(i+1-(2));}
}
  GetElemDataDouble(GetMBS(), edc, "Physics.moment_of_inertia",mass, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.radius",radius, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.length",length, 1);
  return 1;
}


int Rotor1D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Mass1D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rotor1D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Mass1D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rotor1D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Body2D

void Body2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Element::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(pref3D.X(),pref3D.Y(),pref3D.Z(),"reference_position"); ed.SetToolTipText("Reference point for transformation of planar objects to 3D; p = [X, Y, Z]"); edc.TreeAdd("Geometry",ed);

  Matrix rotref3D_1(rotref3D);
  ed.SetMatrix(rotref3D_1.GetMatPtr(),rotref3D_1.Getrows(),rotref3D_1.Getcols(),"rotation_matrix"); ed.SetToolTipText("Rotation matrix for transformation of planar objects to 3D"); edc.TreeAdd("Geometry",ed);

  ed.SetVector3D(size.X(),size.Y(),size.Z(),"general_size"); ed.SetToolTipText("Dimensions of a regular cube [L_x, L_y, L_z]"); edc.TreeAdd("Geometry",ed);
}


int Body2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Element::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.reference_position",pref3D, 1);

  Matrix rotref3D_1;
  GetElemDataMatrix(GetMBS(), edc, "Geometry.rotation_matrix",rotref3D_1, 1);
  rotref3D = Matrix3D(rotref3D_1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.general_size",size, 1);
  return 1;
}


int Body2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Element::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Body2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Element::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Body2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Mass2D

void Mass2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body2D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(drawres,"drawing_tiling"); ed.SetToolTipText("tiling of circle/sphere to represent Mass2D;	the drawing_tiling should be set small in order to improve efficiency,	but large for nice graphical representations"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(size.X(),"radius"); ed.SetToolTipText("drawing radius of mass"); edc.TreeAdd("Graphics",ed);

  {Vector vv((2)-(1)+1); for (int i=1; i<=2; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_position"); ed.SetToolTipText("initial values for position [x,y]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((4)-(3)+1); for (int i=3; i<=4; i++) {vv(i+1-(3))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_velocity"); ed.SetToolTipText("initial values for velocity [vx,vy]"); edc.TreeAdd("Initialization",ed);}
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("total mass of point mass"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Geometry.general_size"); 
}


int Mass2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"general_size"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  Body2D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Graphics.drawing_tiling",drawres, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.radius",size.X(), 1);
  {Vector vv((2)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_position",vv, 1);
  for (int i=1; i<=2; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((4)-(3)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_velocity",vv, 1);
  for (int i=3; i<=4; i++) {x_init(i)=vv(i+1-(3));}
}
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  return 1;
}


int Mass2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body2D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body2D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Rigid2D

void Rigid2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body2D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(Iphi,"moment_of_inertia"); ed.SetToolTipText("[I_ZZ]"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("mass of the body in kg"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(size.X(),size.Y(),size.Z(),"body_dimensions"); ed.SetToolTipText("Dimensions of a regular cube [L_x, L_y, (L_z)]"); edc.TreeAdd("Graphics",ed);

  {Vector vv((2)-(1)+1); for (int i=1; i<=2; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_position"); ed.SetToolTipText("[X, Y]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((5)-(4)+1); for (int i=4; i<=5; i++) {vv(i+1-(4))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_velocity"); ed.SetToolTipText("[vX, vY]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((3)-(3)+1); for (int i=3; i<=3; i++) {vv(i+1-(3))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_rotation"); ed.SetToolTipText("rotation in rad"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((6)-(6)+1); for (int i=6; i<=6; i++) {vv(i+1-(6))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_angular_velocity"); ed.SetToolTipText("Angular velocity in rad/s"); edc.TreeAdd("Initialization",ed);}
  edc.TreeDelete("Geometry.general_size"); 
}


int Rigid2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"general_size"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  Body2D::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Physics.moment_of_inertia",Iphi, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.body_dimensions",size, 1);
  {Vector vv((2)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_position",vv, 1);
  for (int i=1; i<=2; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((5)-(4)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_velocity",vv, 1);
  for (int i=4; i<=5; i++) {x_init(i)=vv(i+1-(4));}
}
  {Vector vv((3)-(3)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_rotation",vv, 1);
  for (int i=3; i<=3; i++) {x_init(i)=vv(i+1-(3));}
}
  {Vector vv((6)-(6)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_angular_velocity",vv, 1);
  for (int i=6; i<=6; i++) {x_init(i)=vv(i+1-(6));}
}
  return 1;
}


int Rigid2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body2D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body2D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Body3D

void Body3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Element::GetElementDataAuto(edc);
  ElementData ed;
}


int Body3D::SetElementDataAuto(ElementDataContainer& edc)
{
  Element::SetElementDataAuto(edc);
  return 1;
}


int Body3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Element::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Body3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Element::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Body3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Mass3D

void Mass3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(drawres,"drawing_tiling"); ed.SetToolTipText("tiling of circle/sphere to represent Sphere"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(size.X(),"radius"); ed.SetToolTipText("drawing radius of mass"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(x_init(1),x_init(1+1),x_init(1+2),"initial_position"); ed.SetToolTipText("coordinates for initial position of mass [X Y Z]"); edc.TreeAdd("Initialization",ed);
  ed.SetVector3D(x_init(4),x_init(4+1),x_init(4+2),"initial_velocity"); ed.SetToolTipText("coordinates for initial velocity of mass [X Y Z]"); edc.TreeAdd("Initialization",ed);
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("total mass of point mass"); edc.TreeAdd("Physics",ed);
}


int Mass3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Body3D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Graphics.drawing_tiling",drawres, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.radius",size.X(), 1);
  {Vector3D vv;  GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_position",vv, 1); vv.Get(x_init(1),x_init(1+1),x_init(1+2));}
  {Vector3D vv;  GetElemDataVector3D(GetMBS(), edc, "Initialization.initial_velocity",vv, 1); vv.Get(x_init(4),x_init(4+1),x_init(4+2));}
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  return 1;
}


int Mass3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Mass3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=NodalMass3D

void NodalMass3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Mass3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(nodenum,"node_number"); ed.SetToolTipText("node number to which the mass refers"); edc.TreeAdd("",ed);
}


int NodalMass3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Mass3D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "node_number",nodenum, 1);
  return 1;
}


int NodalMass3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Mass3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int NodalMass3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Mass3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int NodalMass3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=NodalDiskMass3D

void NodalDiskMass3D::GetElementDataAuto(ElementDataContainer& edc)
{
  NodalMass3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetBool(full_mass_matrix,"full_mass_matrix"); ed.SetToolTipText("set to 1 if influence of tilted mass should be considered in the mass matrix"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(Ip.X(),Ip.Y(),Ip.Z(),"moment_of_inertia"); ed.SetToolTipText("moments of inertia of the disk"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Graphics.radius"); 
  ed.SetDouble(size.X(),"thickness"); ed.SetToolTipText("drawing thickness of disk mass"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(size.Y(),"radius"); ed.SetToolTipText("drawing radius of disk mass"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.mass"); 
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("total mass of disk"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Initialization.initial_position"); 
  edc.TreeDelete("Initialization.initial_velocity"); 
}


int NodalDiskMass3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"radius"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"mass"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"initial_position"); ed.SetLocked(1); edc.TreeAdd("Initialization",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"initial_velocity"); ed.SetLocked(1); edc.TreeAdd("Initialization",ed);}


  NodalMass3D::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "Physics.full_mass_matrix",full_mass_matrix, 1);
  GetElemDataVector3D(GetMBS(), edc, "Physics.moment_of_inertia",Ip, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.thickness",size.X(), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.radius",size.Y(), 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  return 1;
}


int NodalDiskMass3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = NodalMass3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int NodalDiskMass3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = NodalMass3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int NodalDiskMass3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Rigid3D

void Rigid3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;

  Matrix Iphi_1(Iphi);
  ed.SetMatrix(Iphi_1.GetMatPtr(),Iphi_1.Getrows(),Iphi_1.Getcols(),"moment_of_inertia"); ed.SetToolTipText("[I_XX,I_XY,I_XZ; ...]"); edc.TreeAdd("Physics",ed);

  ed.SetDouble(volume,"volume"); ed.SetToolTipText("volume of the body in m*m*m"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(size.X(),size.Y(),size.Z(),"body_dimensions"); ed.SetToolTipText("Dimensions of a regular cube [L_x, L_y, L_z] in m"); edc.TreeAdd("Graphics",ed);

  {Vector vv((3)-(1)+1); for (int i=1; i<=3; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_position"); ed.SetToolTipText("[X, Y, Z]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((NRotParam()+6)-(NRotParam()+4)+1); for (int i=NRotParam()+4; i<=NRotParam()+6; i++) {vv(i+1-(NRotParam()+4))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"initial_velocity"); ed.SetToolTipText("[X, Y, Z]"); edc.TreeAdd("Initialization",ed);}
  ed.SetDouble(mass,"mass"); ed.SetToolTipText("mass of the body in kg"); edc.TreeAdd("Physics",ed);
}


int Rigid3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Body3D::SetElementDataAuto(edc);

  Matrix Iphi_1;
  GetElemDataMatrix(GetMBS(), edc, "Physics.moment_of_inertia",Iphi_1, 1);
  Iphi = Matrix3D(Iphi_1);
  GetElemDataDouble(GetMBS(), edc, "Physics.volume",volume, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.body_dimensions",size, 1);
  {Vector vv((3)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_position",vv, 1);
  for (int i=1; i<=3; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((NRotParam()+6)-(NRotParam()+4)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.initial_velocity",vv, 1);
  for (int i=NRotParam()+4; i<=NRotParam()+6; i++) {x_init(i)=vv(i+1-(NRotParam()+4));}
}
  GetElemDataDouble(GetMBS(), edc, "Physics.mass",mass, 1);
  return 1;
}


int Rigid3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Rigid3DKardan

void Rigid3DKardan::GetElementDataAuto(ElementDataContainer& edc)
{
  Rigid3D::GetElementDataAuto(edc);
  ElementData ed;
}


int Rigid3DKardan::SetElementDataAuto(ElementDataContainer& edc)
{
  Rigid3D::SetElementDataAuto(edc);
  return 1;
}


int Rigid3DKardan::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Rigid3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid3DKardan::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Rigid3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rigid3DKardan::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Beam3D

void Beam3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(n1,"node_1",0,0,1); ed.SetToolTipText("number of Node 1"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(n2,"node_2",0,0,1); ed.SetToolTipText("number of Node 2"); edc.TreeAdd("Geometry",ed);
  ed.SetBool(axialdeformation,"axial_deformation"); ed.SetToolTipText("include effect of axial deformation"); edc.TreeAdd("Physics",ed);
  ed.SetInt(materialnum,"material_number"); ed.SetToolTipText("material number which contains the main material properties of the beam"); edc.TreeAdd("Physics",ed);
}


int Beam3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Body3D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_1",n1, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_2",n2, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.axial_deformation",axialdeformation, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.material_number",materialnum, 1);
  return 1;
}


int Beam3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Beam3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RotorBeamXAxis

void RotorBeamXAxis::GetElementDataAuto(ElementDataContainer& edc)
{
  Beam3D::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Geometry.body_dimensions"); 
}


int RotorBeamXAxis::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"body_dimensions"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  Beam3D::SetElementDataAuto(edc);
  return 1;
}


int RotorBeamXAxis::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Beam3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RotorBeamXAxis::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Beam3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RotorBeamXAxis::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFCable2D

void ANCFCable2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body2D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(n1,"node_number1"); ed.SetToolTipText("global number of node 1 (node must exist, add with AddNode(...))"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(n2,"node_number2"); ed.SetToolTipText("global number of node 2 (node must exist, add with AddNode(...))"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(concentratedmass1,"concentrated_mass_node1"); ed.SetToolTipText("a concentrated mass is added at position of node 1"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(concentratedmass2,"concentrated_mass_node2"); ed.SetToolTipText("a concentrated mass is added at position of node 2"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(size.X(),size.Y(),size.Z(),"element_size"); ed.SetToolTipText("size of element: length (lx), height (ly), out of plane width (lz)"); edc.TreeAdd("Geometry",ed);

  {Vector vv((2)-(1)+1); for (int i=1; i<=2; i++) {vv(i+1-(1))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_position"); ed.SetToolTipText("initial values for position of node 1 [x1,y1]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((4)-(3)+1); for (int i=3; i<=4; i++) {vv(i+1-(3))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_slope"); ed.SetToolTipText("initial values for slope vector of node 1 [r1xx,r1xy]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((6)-(5)+1); for (int i=5; i<=6; i++) {vv(i+1-(5))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_position"); ed.SetToolTipText("initial values for position of node 2 [x2,y2]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((8)-(7)+1); for (int i=7; i<=8; i++) {vv(i+1-(7))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_slope"); ed.SetToolTipText("initial values for slope vector of node 2 [r2xx,r2xy]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((2+8)-(1+8)+1); for (int i=1+8; i<=2+8; i++) {vv(i+1-(1+8))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_velocity"); ed.SetToolTipText("initial values for velocity of node 1 [x1,y1]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((4+8)-(3+8)+1); for (int i=3+8; i<=4+8; i++) {vv(i+1-(3+8))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_slope_vel"); ed.SetToolTipText("initial values for slope velocity vector of node 1 [r1xx,r1xy]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((6+8)-(5+8)+1); for (int i=5+8; i<=6+8; i++) {vv(i+1-(5+8))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_velocity"); ed.SetToolTipText("initial values for velocity of node 2 [x2,y2]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((8+8)-(7+8)+1); for (int i=7+8; i<=8+8; i++) {vv(i+1-(7+8))=x_init(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_slope_vel"); ed.SetToolTipText("initial values for slope velocity vector of node 2 [r2xx,r2xy]"); edc.TreeAdd("Initialization",ed);}
}


int ANCFCable2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Body2D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number1",n1, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number2",n2, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.concentrated_mass_node1",concentratedmass1, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.concentrated_mass_node2",concentratedmass2, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.element_size",size, 1);
  {Vector vv((2)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_initial_position",vv, 1);
  for (int i=1; i<=2; i++) {x_init(i)=vv(i+1-(1));}
}
  {Vector vv((4)-(3)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_initial_slope",vv, 1);
  for (int i=3; i<=4; i++) {x_init(i)=vv(i+1-(3));}
}
  {Vector vv((6)-(5)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_initial_position",vv, 1);
  for (int i=5; i<=6; i++) {x_init(i)=vv(i+1-(5));}
}
  {Vector vv((8)-(7)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_initial_slope",vv, 1);
  for (int i=7; i<=8; i++) {x_init(i)=vv(i+1-(7));}
}
  {Vector vv((2+8)-(1+8)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_initial_velocity",vv, 1);
  for (int i=1+8; i<=2+8; i++) {x_init(i)=vv(i+1-(1+8));}
}
  {Vector vv((4+8)-(3+8)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_initial_slope_vel",vv, 1);
  for (int i=3+8; i<=4+8; i++) {x_init(i)=vv(i+1-(3+8));}
}
  {Vector vv((6+8)-(5+8)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_initial_velocity",vv, 1);
  for (int i=5+8; i<=6+8; i++) {x_init(i)=vv(i+1-(5+8));}
}
  {Vector vv((8+8)-(7+8)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_initial_slope_vel",vv, 1);
  for (int i=7+8; i<=8+8; i++) {x_init(i)=vv(i+1-(7+8));}
}
  return 1;
}


int ANCFCable2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body2D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFCable2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body2D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFCable2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFCable3D

void ANCFCable3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;
}


int ANCFCable3D::SetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::SetElementDataAuto(edc);
  return 1;
}


int ANCFCable3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFCable3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFCable3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShear3DGeneric

void ANCFBeamShear3DGeneric::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetBool(is_straight_beam_in_reference_configuration,"straight_beam"); ed.SetToolTipText("is straight beam in reference configuration"); edc.TreeAdd("ShearBeam",ed);
  ed.SetVector3D(size.X(),size.Y(),size.Z(),"body_dimensions",0,0,1); ed.SetToolTipText("dimensions of the beam. [L_x (length), L_y (width), L_z (height)]"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(materialnum,"material_number"); ed.SetToolTipText("material number which contains the main material properties of the beam"); edc.TreeAdd("Physics",ed);
  ed.SetBool(perform_reduced_integration,"reduced_integration"); ed.SetToolTipText("reduced integration in cont. mech. formulation (CMF)"); edc.TreeAdd("ShearBeam",ed);
  ed.SetInt(BeamFormulation,"beamformulation"); ed.SetToolTipText("2 = CMF, 4 = SMF"); edc.TreeAdd("ShearBeam",ed);
  ed.SetBool(calclinear,"calc_linear"); ed.SetToolTipText("linearized strain computation in cont. mech. formulation (CMF)"); edc.TreeAdd("ShearBeam",ed);
  ed.SetInt(n1,"node_number1"); ed.SetToolTipText("global number of node 1 (left), node must already exist"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(n2,"node_number2"); ed.SetToolTipText("global number of node 2 (right), node must already exist"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(nnodes,"nnodes"); ed.SetLocked(1); ed.SetToolTipText("number of nodes"); edc.TreeAdd("ShearBeam",ed);
  ed.SetInt(order_axial,"integration_order_axial"); ed.SetToolTipText("axial integration order"); edc.TreeAdd("ShearBeam",ed);
  ed.SetInt(order_crosssectional,"integration_order_cross_section"); ed.SetToolTipText("cross section integration order, taks effect only in cont. mech. formulation (CMF)"); edc.TreeAdd("ShearBeam",ed);

  {Vector vv((3)-(1)+1); for (int i=1; i<=3; i++) {vv(i+1-(1))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_reference_position"); ed.SetToolTipText("position of node 1 (left) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((6)-(4)+1); for (int i=4; i<=6; i++) {vv(i+1-(4))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_reference_slope_2"); ed.SetToolTipText("slope vector 2 of node 1 (left) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((9)-(7)+1); for (int i=7; i<=9; i++) {vv(i+1-(7))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_reference_slope_3"); ed.SetToolTipText("slope vector 3 of node 1 (left) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((12)-(10)+1); for (int i=10; i<=12; i++) {vv(i+1-(10))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_reference_position"); ed.SetToolTipText("position of node 2 (right) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((15)-(13)+1); for (int i=13; i<=15; i++) {vv(i+1-(13))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_reference_slope_2"); ed.SetToolTipText("slope vector 2 of node 2 (right) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((18)-(16)+1); for (int i=16; i<=18; i++) {vv(i+1-(16))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_reference_slope_3"); ed.SetToolTipText("slope vector 3 of node 2 (right) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((3)-(1)+1); for (int i=1; i<=3; i++) {vv(i+1-(1))=penalty(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"penalty_kappa"); ed.SetToolTipText("penalty term for kappa [kappa1,kappa2,kappa3]"); edc.TreeAdd("ShearBeam",ed);}

  {Vector vv((6)-(4)+1); for (int i=4; i<=6; i++) {vv(i+1-(4))=penalty(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"penalty_gamma"); ed.SetToolTipText("penalty term for gamma [gamma1,gamma2,gamma3]"); edc.TreeAdd("ShearBeam",ed);}

  {Vector vv((3)-(1)+1); for (int i=1; i<=3; i++) {vv(i+1-(1))=penalty(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"penalty_E"); ed.SetToolTipText("penalty term for green lagrange strains (E) [Eyy,Ezz,Eyz]"); edc.TreeAdd("ShearBeam",ed);}
  edc.TreeDelete("Physics.density"); 
  edc.TreeDelete("Physics.mass"); 
  edc.TreeDelete("Physics.mass_prop_damping"); 
}


int ANCFBeamShear3DGeneric::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"density"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"mass"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"mass_prop_damping"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}


  Body3D::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "ShearBeam.straight_beam",is_straight_beam_in_reference_configuration, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.body_dimensions",size, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.material_number",materialnum, 1);
  GetElemDataBool(GetMBS(), edc, "ShearBeam.reduced_integration",perform_reduced_integration, 1);
  GetElemDataInt(GetMBS(), edc, "ShearBeam.beamformulation",BeamFormulation, 1);
  GetElemDataBool(GetMBS(), edc, "ShearBeam.calc_linear",calclinear, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number1",n1, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number2",n2, 1);
  GetElemDataInt(GetMBS(), edc, "ShearBeam.integration_order_axial",order_axial, 1);
  GetElemDataInt(GetMBS(), edc, "ShearBeam.integration_order_cross_section",order_crosssectional, 1);
  {Vector vv((3)-(1)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_reference_position",vv, 1);
  for (int i=1; i<=3; i++) {q0(i)=vv(i+1-(1));}
}
  {Vector vv((6)-(4)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_reference_slope_2",vv, 1);
  for (int i=4; i<=6; i++) {q0(i)=vv(i+1-(4));}
}
  {Vector vv((9)-(7)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node1_reference_slope_3",vv, 1);
  for (int i=7; i<=9; i++) {q0(i)=vv(i+1-(7));}
}
  {Vector vv((12)-(10)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_reference_position",vv, 1);
  for (int i=10; i<=12; i++) {q0(i)=vv(i+1-(10));}
}
  {Vector vv((15)-(13)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_reference_slope_2",vv, 1);
  for (int i=13; i<=15; i++) {q0(i)=vv(i+1-(13));}
}
  {Vector vv((18)-(16)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node2_reference_slope_3",vv, 1);
  for (int i=16; i<=18; i++) {q0(i)=vv(i+1-(16));}
}
  {Vector vv((3)-(1)+1);  GetElemDataVector(GetMBS(), edc, "ShearBeam.penalty_kappa",vv, 1);
  for (int i=1; i<=3; i++) {penalty(i)=vv(i+1-(1));}
}
  {Vector vv((6)-(4)+1);  GetElemDataVector(GetMBS(), edc, "ShearBeam.penalty_gamma",vv, 1);
  for (int i=4; i<=6; i++) {penalty(i)=vv(i+1-(4));}
}
  {Vector vv((3)-(1)+1);  GetElemDataVector(GetMBS(), edc, "ShearBeam.penalty_E",vv, 1);
  for (int i=1; i<=3; i++) {penalty(i)=vv(i+1-(1));}
}
  return 1;
}


int ANCFBeamShear3DGeneric::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DGeneric::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DGeneric::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShear3DLinear

void ANCFBeamShear3DLinear::GetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShear3DGeneric::GetElementDataAuto(edc);
  ElementData ed;
}


int ANCFBeamShear3DLinear::SetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShear3DGeneric::SetElementDataAuto(edc);
  return 1;
}


int ANCFBeamShear3DLinear::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = ANCFBeamShear3DGeneric::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DLinear::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = ANCFBeamShear3DGeneric::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DLinear::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShear3DQuadratic

void ANCFBeamShear3DQuadratic::GetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShear3DGeneric::GetElementDataAuto(edc);
  ElementData ed;

  {Vector vv((21)-(19)+1); for (int i=19; i<=21; i++) {vv(i+1-(19))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node3_reference_position"); ed.SetToolTipText("position of node 3 (middle) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((24)-(22)+1); for (int i=22; i<=24; i++) {vv(i+1-(22))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node3_reference_slope_2"); ed.SetToolTipText("slope vector 2 of node 3 (middle) in reference configuration."); edc.TreeAdd("Initialization",ed);}

  {Vector vv((27)-(25)+1); for (int i=25; i<=27; i++) {vv(i+1-(25))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node3_reference_slope_3"); ed.SetToolTipText("slope vector 3 of node 3 (middle) in reference configuration."); edc.TreeAdd("Initialization",ed);}
  ed.SetInt(n3,"node_number3"); ed.SetToolTipText("global number of node 3 (middle), node must already exist"); edc.TreeAdd("Geometry",ed);
}


int ANCFBeamShear3DQuadratic::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  ANCFBeamShear3DGeneric::SetElementDataAuto(edc);
  {Vector vv((21)-(19)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node3_reference_position",vv, 1);
  for (int i=19; i<=21; i++) {q0(i)=vv(i+1-(19));}
}
  {Vector vv((24)-(22)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node3_reference_slope_2",vv, 1);
  for (int i=22; i<=24; i++) {q0(i)=vv(i+1-(22));}
}
  {Vector vv((27)-(25)+1);  GetElemDataVector(GetMBS(), edc, "Initialization.node3_reference_slope_3",vv, 1);
  for (int i=25; i<=27; i++) {q0(i)=vv(i+1-(25));}
}
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number3",n3, 1);
  return 1;
}


int ANCFBeamShear3DQuadratic::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = ANCFBeamShear3DGeneric::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DQuadratic::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = ANCFBeamShear3DGeneric::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShear3DQuadratic::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeam3DTorsion

void ANCFBeam3DTorsion::GetElementDataAuto(ElementDataContainer& edc)
{
  Body3D::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(n1,"node_number1"); ed.SetToolTipText("global number of node 1 (left), node must already exist"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(n2,"node_number2"); ed.SetToolTipText("global number of node 2 (right), node must already exist"); edc.TreeAdd("Geometry",ed);
  ed.SetBool(do_update_directors,"update_directors"); ed.SetToolTipText("update directors during calculation"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(kinematic_computation_mode,"kinematic_computation_mode"); ed.SetToolTipText("0 .. exact kinematic terms + 5th order gaussian integration (slow), 1 .. exact terms + 1st order lobatto integration (fast), 2 .. constant mass matrix approximation (fastest)"); edc.TreeAdd("Computation",ed);
  ed.SetInt(intorder_mass,"mass"); ed.SetToolTipText("integration order for mass terms"); edc.TreeAdd("Computation.IntegrationOrder",ed);
  ed.SetInt(intorder_axial_strain,"axial_strain"); ed.SetToolTipText("integration order for work of axial strain"); edc.TreeAdd("Computation.IntegrationOrder",ed);
  ed.SetInt(intorder_curvature,"curvature"); ed.SetToolTipText("integration order for work of curvature"); edc.TreeAdd("Computation.IntegrationOrder",ed);
  ed.SetInt(materialnum,"material_number"); ed.SetToolTipText("material number which contains the main material properties of the beam"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(size.X(),"beam_length"); ed.SetLocked(1); ed.SetToolTipText("length of the beam, calculated by numerical integration of |r'|"); edc.TreeAdd("Info",ed);
}


int ANCFBeam3DTorsion::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Body3D::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number1",n1, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_number2",n2, 1);
  GetElemDataBool(GetMBS(), edc, "Geometry.update_directors",do_update_directors, 1);
  GetElemDataInt(GetMBS(), edc, "Computation.kinematic_computation_mode",kinematic_computation_mode, 1);
  GetElemDataInt(GetMBS(), edc, "Computation.IntegrationOrder.mass",intorder_mass, 1);
  GetElemDataInt(GetMBS(), edc, "Computation.IntegrationOrder.axial_strain",intorder_axial_strain, 1);
  GetElemDataInt(GetMBS(), edc, "Computation.IntegrationOrder.curvature",intorder_curvature, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.material_number",materialnum, 1);
  return 1;
}


int ANCFBeam3DTorsion::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Body3D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeam3DTorsion::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Body3D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeam3DTorsion::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=FiniteElementGenericBeam2D

void FiniteElementGenericBeam2D::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  ed.SetText(GetElementSpecification().c_str(),"element_type"); ed.SetToolTipText("specification of element type. Once the element is added to the mbs, you MUST NOT change this type anymore!"); edc.TreeAdd("",ed);
  ed.SetInt(plasticip_x,"number_plastic_intpoints_x"); ed.SetToolTipText("Number of integration points for plastic variable in x direction"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(plasticip_y,"number_plastic_intpoints_y"); ed.SetToolTipText("Number of integration points for plastic variable in y direction"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(elnum,"element_number"); ed.SetLocked(1); ed.SetToolTipText("number of the element in the mbs"); edc.TreeAdd("",ed);
  ed.SetText(elementname.c_str(),"name"); ed.SetToolTipText("name of the element"); edc.TreeAdd("",ed);

  {Vector vv(nodes.Length()); for(int i = 1; i <= nodes.Length(); i++) {vv(i) = nodes(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_numbers"); ed.SetValuesInt(); ed.SetToolTipText("Global node numbers"); edc.TreeAdd("Geometry",ed);
}
  ed.SetInt(materialnum,"material_number"); ed.SetToolTipText("material number which contains the main material properties of the body or element"); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"RGB_color"); ed.SetToolTipText("[red, green, blue] color of element, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  if(!IsType(TConstraint))
  {
  {Vector vv(loads.Length()); for(int i = 1; i <= loads.Length(); i++) {vv(i) = loads(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"loads"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("Set loads attached to this element: 'nr_load1, nr_load2, ...' or empty"); edc.TreeAdd("",ed);
}
  }

  {Vector vv(sensors.Length()); for(int i = 1; i <= sensors.Length(); i++) {vv(i) = sensors(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"sensors"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetVariableLength(); ed.SetToolTipText("attached sensors"); edc.TreeAdd("Info",ed);
}
}


int FiniteElementGenericBeam2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GetElemDataText(GetMBS(), edc, "element_type",GetElementSpecification(), 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.number_plastic_intpoints_x",plasticip_x, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.number_plastic_intpoints_y",plasticip_y, 1);
  GetElemDataText(GetMBS(), edc, "name",elementname, 1);
  GetElemDataIVector(GetMBS(), edc, "Geometry.node_numbers",nodes, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.material_number",materialnum, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.RGB_color",col, 1);
  if(!IsType(TConstraint))
  {  GetElemDataIVector(GetMBS(), edc, "loads",loads, 1);
  }
  return 1;
}


int FiniteElementGenericBeam2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  return 0;
}


int FiniteElementGenericBeam2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  return 0;
}


int FiniteElementGenericBeam2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShearFE2DGeneric

void ANCFBeamShearFE2DGeneric::GetElementDataAuto(ElementDataContainer& edc)
{
  FiniteElementGenericBeam2D::GetElementDataAuto(edc);
  ElementData ed;

  {Vector vv((2)-(1)+1); for (int i=1; i<=2; i++) {vv(i+1-(1))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_position"); ed.SetLocked(1); ed.SetToolTipText("initial values for position of node 1 (left). r1 = [x1,y1]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((4)-(3)+1); for (int i=3; i<=4; i++) {vv(i+1-(3))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node1_initial_slope_eta"); ed.SetLocked(1); ed.SetToolTipText("initial values for slope vector of node 1 (left). [d(r1)/d(eta)]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((6)-(5)+1); for (int i=5; i<=6; i++) {vv(i+1-(5))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_position"); ed.SetLocked(1); ed.SetToolTipText("initial values for position of node 2 (right). r2 = [x2,y2]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((8)-(7)+1); for (int i=7; i<=8; i++) {vv(i+1-(7))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node2_initial_slope_eta"); ed.SetLocked(1); ed.SetToolTipText("initial values for slope vector of node 2 (right). [d(r2)/d(eta)]"); edc.TreeAdd("Initialization",ed);}
  ed.SetDouble(size(1),"lx"); ed.SetToolTipText("length of undeformed beam element l_x"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(size(2),"ly"); ed.SetToolTipText("length of undeformed beam element l_y"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(size(3),"lz"); ed.SetToolTipText("length of undeformed beam element l_z"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(use_reduced_integration,"use_reduced_integration"); ed.SetToolTipText("use reduced integration along the beam axis"); edc.TreeAdd("Initialization",ed);
  ed.SetInt(use_reduced_integration_poisson,"use_reduced_integration_poisson"); ed.SetToolTipText("use reduced integration across beam thickness to avoid Poisson locking"); edc.TreeAdd("Initialization",ed);
  ed.SetInt(use_contmech_formulation,"use_continuum_mechanics"); ed.SetToolTipText("use continuum mechanics formulation [1] or structural mechanics formulation [0]"); edc.TreeAdd("Initialization",ed);
  ed.SetDouble(shear_correction_factor,"shear_correction_factor"); ed.SetToolTipText("shear correction factor"); edc.TreeAdd("Initialization",ed);
}


int ANCFBeamShearFE2DGeneric::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  FiniteElementGenericBeam2D::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Geometry.lx",size(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.ly",size(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.lz",size(3), 1);
  GetElemDataInt(GetMBS(), edc, "Initialization.use_reduced_integration",use_reduced_integration, 1);
  GetElemDataInt(GetMBS(), edc, "Initialization.use_reduced_integration_poisson",use_reduced_integration_poisson, 1);
  GetElemDataInt(GetMBS(), edc, "Initialization.use_continuum_mechanics",use_contmech_formulation, 1);
  GetElemDataDouble(GetMBS(), edc, "Initialization.shear_correction_factor",shear_correction_factor, 1);
  return 1;
}


int ANCFBeamShearFE2DGeneric::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = FiniteElementGenericBeam2D::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DGeneric::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = FiniteElementGenericBeam2D::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DGeneric::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShearFE2DLinear

void ANCFBeamShearFE2DLinear::GetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShearFE2DGeneric::GetElementDataAuto(edc);
  ElementData ed;
}


int ANCFBeamShearFE2DLinear::SetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShearFE2DGeneric::SetElementDataAuto(edc);
  return 1;
}


int ANCFBeamShearFE2DLinear::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = ANCFBeamShearFE2DGeneric::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DLinear::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = ANCFBeamShearFE2DGeneric::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DLinear::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFBeamShearFE2DQuadratic

void ANCFBeamShearFE2DQuadratic::GetElementDataAuto(ElementDataContainer& edc)
{
  ANCFBeamShearFE2DGeneric::GetElementDataAuto(edc);
  ElementData ed;

  {Vector vv((10)-(9)+1); for (int i=9; i<=10; i++) {vv(i+1-(9))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node3_initial_position"); ed.SetLocked(1); ed.SetToolTipText("initial values for position of node 3 (middle). r3 = [x3,y3]"); edc.TreeAdd("Initialization",ed);}

  {Vector vv((12)-(11)+1); for (int i=11; i<=12; i++) {vv(i+1-(11))=q0(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node3_initial_slope_eta"); ed.SetLocked(1); ed.SetToolTipText("initial values for slope vector of node 3 (right). [d(r3)/d(eta)]"); edc.TreeAdd("Initialization",ed);}
}


int ANCFBeamShearFE2DQuadratic::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  ANCFBeamShearFE2DGeneric::SetElementDataAuto(edc);
  return 1;
}


int ANCFBeamShearFE2DQuadratic::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = ANCFBeamShearFE2DGeneric::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DQuadratic::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = ANCFBeamShearFE2DGeneric::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFBeamShearFE2DQuadratic::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ANCFSimpleThinPlate3D

void ANCFSimpleThinPlate3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Element::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(size1,"size1"); ed.SetToolTipText("plate dimension in x-direction"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(size2,"size2"); ed.SetToolTipText("plate dimension in y-direction"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(thickness,"thickness"); ed.SetToolTipText("plate thickness"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(nodes(1),"node_1"); ed.SetToolTipText("node1 of the element"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(nodes(2),"node_2"); ed.SetToolTipText("node2 of the element"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(nodes(3),"node_3"); ed.SetToolTipText("node3 of the element"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(nodes(4),"node_4"); ed.SetToolTipText("node4 of the element"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(materialnum,"material_number"); ed.SetToolTipText("material number containing the material properties of the plate"); edc.TreeAdd("Physics",ed);
  ed.SetBool((int&)geometricNonlinearityStatus,"is_geometrically_nonlinear"); ed.SetToolTipText("0 ... geometrically linear (small deformation), 1 ... geometrically nonlinear"); edc.TreeAdd("Physics",ed);
}


int ANCFSimpleThinPlate3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Element::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Geometry.size1",size1, 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.size2",size2, 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.thickness",thickness, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_1",nodes(1), 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_2",nodes(2), 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_3",nodes(3), 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.node_4",nodes(4), 1);
  GetElemDataInt(GetMBS(), edc, "Physics.material_number",materialnum, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.is_geometrically_nonlinear",(int&)geometricNonlinearityStatus, 1);
  return 1;
}


int ANCFSimpleThinPlate3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Element::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFSimpleThinPlate3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Element::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ANCFSimpleThinPlate3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Constraint

void Constraint::GetElementDataAuto(ElementDataContainer& edc)
{
  Element::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetBool(use_local_coordinate_system,"use_local_coordinate_system"); ed.SetToolTipText("0=use global coordinates, 1=use local coordinate system of Body 1"); edc.TreeAdd("Geometry",ed);
  ed.SetBool(use_penalty_formulation,"use_penalty_formulation"); ed.SetToolTipText("0 = use lagrange multipliers (index 3 DAE, exact), 1 = use penalty formulation (no additional equation added, approximate constraint)"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(spring_stiffness,"spring_stiffness"); ed.SetToolTipText("general or penalty stiffness parameter"); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Graphics.show_element"); 
  ed.SetBool(draw_element,"show_connector"); ed.SetToolTipText("Flag to draw connector"); edc.TreeAdd("Graphics",ed);
}


int Constraint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"show_element"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  Element::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "Geometry.use_local_coordinate_system",use_local_coordinate_system, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.use_penalty_formulation",use_penalty_formulation, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.spring_stiffness",spring_stiffness, 1);
  GetElemDataBool(GetMBS(), edc, "Graphics.show_connector",draw_element, 1);
  return 1;
}


int Constraint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Element::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Constraint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Element::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Constraint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=BasePointJoint

void BasePointJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetBool(stiffness_in_joint_local_frame,"use_joint_local_frame"); ed.SetToolTipText("Use a special joint local frame"); edc.TreeAdd("Geometry",ed);
  ed.SetInt(elements(1),"element_number",1,0,1); ed.SetToolTipText("Number of constrained element"); edc.TreeAdd("Position1",ed);
  ed.SetInt(elements(2),"element_number",0,0,1); ed.SetToolTipText("Number of constrained element"); edc.TreeAdd("Position2",ed);
  ed.SetVector3D(loccoords(1).X(),loccoords(1).Y(),loccoords(1).Z(),"position"); ed.SetToolTipText("local position. Only used if node_number == 0!"); edc.TreeAdd("Position1",ed);
  ed.SetVector3D(loccoords(2).X(),loccoords(2).Y(),loccoords(2).Z(),"position"); ed.SetToolTipText("local or global (if element_number == 0) position. Only used if node_number == 0!"); edc.TreeAdd("Position2",ed);
  ed.SetInt(nodes(1),"node_number",1,0,1); ed.SetToolTipText("local or global (if element_number == 0) node number."); edc.TreeAdd("Position1",ed);
  ed.SetInt(nodes(2),"node_number",1,0,1); ed.SetToolTipText("local or global (if element_number == 0) node number."); edc.TreeAdd("Position2",ed);

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetToolTipText("[x,y,z]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetVector3D(spring_stiffness3.X(),spring_stiffness3.Y(),spring_stiffness3.Z(),"spring_stiffness_vector"); ed.SetToolTipText("penalty stiffness parameter [kx,ky,kz]. Just used if scalar spring_stiffness == 0, otherwise kx=ky=kz=spring_stiffness"); edc.TreeAdd("Physics.Penalty",ed);

  Matrix JA0i_1(JA0i);
  ed.SetMatrix(JA0i_1.GetMatPtr(),JA0i_1.Getrows(),JA0i_1.Getcols(),"joint_local_frame"); ed.SetToolTipText("Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body 1. Just used if use_joint_local_frame == 1"); edc.TreeAdd("Geometry",ed);

  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  ed.SetBool(use_local_coordinate_system,"use_local_coordinate_system"); ed.SetToolTipText("0=use global coordinates, 1=use local coordinate system of Body 1"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient for viscous damping (F = d*v), applied in all constrained directions"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_local_frame_size,"draw_size_joint_local_frame"); ed.SetToolTipText("drawing dimensions of joint local frame. If set to -1, than global_draw_scalar_size is used. If set to 0, than no joint local frame is drawn."); edc.TreeAdd("Graphics",ed);
}


int BasePointJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "Geometry.use_joint_local_frame",stiffness_in_joint_local_frame, 1);
  GetElemDataInt(GetMBS(), edc, "Position1.element_number",elements(1), 1);
  GetElemDataInt(GetMBS(), edc, "Position2.element_number",elements(2), 1);
  GetElemDataVector3D(GetMBS(), edc, "Position1.position",loccoords(1), 1);
  GetElemDataVector3D(GetMBS(), edc, "Position2.position",loccoords(2), 1);
  GetElemDataInt(GetMBS(), edc, "Position1.node_number",nodes(1), 1);
  GetElemDataInt(GetMBS(), edc, "Position2.node_number",nodes(2), 1);
  GetElemDataIVector(GetMBS(), edc, "Physics.Lagrange.constrained_directions",dir, 1);
  GetElemDataVector3D(GetMBS(), edc, "Physics.Penalty.spring_stiffness_vector",spring_stiffness3, 1);

  Matrix JA0i_1;
  GetElemDataMatrix(GetMBS(), edc, "Geometry.joint_local_frame",JA0i_1, 1);
  JA0i = Matrix3D(JA0i_1);
  GetElemDataBool(GetMBS(), edc, "Geometry.use_local_coordinate_system",use_local_coordinate_system, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size_joint_local_frame",draw_local_frame_size, 1);
  return 1;
}


int BasePointJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int BasePointJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int BasePointJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=CoordConstraint

void CoordConstraint::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(coord_offset,"coord_offset"); ed.SetToolTipText("coordinate offset d, see documentation section equation"); edc.TreeAdd("",ed);
  ed.SetDouble(coord_gain_factor,"coord_gain_factor"); ed.SetToolTipText("coordinate gain factor k, see documentation section equation"); edc.TreeAdd("",ed);
  ed.SetBool(relative_to_inital_values,"relative_to_inital_values"); ed.SetToolTipText("flag == 1: full equation is used, see documentation; flag == 0: the init state values qi0 and qj0 are neglected."); edc.TreeAdd("",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient Dp for viscous damping"); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(spring_stiffness,"spring_stiffness"); ed.SetToolTipText("general or penalty stiffness parameter Sp"); edc.TreeAdd("Physics.Penalty",ed);
}


int CoordConstraint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "coord_offset",coord_offset, 1);
  GetElemDataDouble(GetMBS(), edc, "coord_gain_factor",coord_gain_factor, 1);
  GetElemDataBool(GetMBS(), edc, "relative_to_inital_values",relative_to_inital_values, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.spring_stiffness",spring_stiffness, 1);
  return 1;
}


int CoordConstraint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int CoordConstraint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int CoordConstraint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=VelCoordConstraint

void VelCoordConstraint::GetElementDataAuto(ElementDataContainer& edc)
{
  CoordConstraint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Penalty.damping"); 
  edc.TreeDelete("relative_to_inital_values"); 
  ed.SetBool(relative_to_inital_values,"relative_to_inital_values"); ed.SetToolTipText("flag == 1: full equation is used, see documentation; flag == 0: the init state derivatives d(qi0)/dt and d(qj0)/dt are neglected."); edc.TreeAdd("",ed);
}


int VelCoordConstraint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"relative_to_inital_values"); ed.SetLocked(1); edc.TreeAdd("",ed);}


  CoordConstraint::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "relative_to_inital_values",relative_to_inital_values, 1);
  return 1;
}


int VelCoordConstraint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = CoordConstraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int VelCoordConstraint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = CoordConstraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int VelCoordConstraint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RollingJoint3D

void RollingJoint3D::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(rolling_plane_point.X(),rolling_plane_point.Y(),rolling_plane_point.Z(),"rolling_plane_point"); ed.SetToolTipText("any point at plane on that wheel is rolling"); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(rolling_plane_normal.X(),rolling_plane_normal.Y(),rolling_plane_normal.Z(),"rolling_plane_normal"); ed.SetToolTipText("normal vector of plane on that wheel is rolling"); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(wheel_local_center_point.X(),wheel_local_center_point.Y(),wheel_local_center_point.Z(),"wheel_local_center_point"); ed.SetToolTipText("wheel center point in local coordinates of wheel"); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(wheel_local_axis.X(),wheel_local_axis.Y(),wheel_local_axis.Z(),"wheel_local_axis"); ed.SetToolTipText("wheel axis in local coordinates of wheel"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(wheel_radius,"wheel_radius"); ed.SetToolTipText("radius of the idealized wheel"); edc.TreeAdd("Geometry",ed);
  ed.SetBool(use_penalty_inplane_dir,"use_penalty_inplane"); ed.SetToolTipText("1 = enabled, 0 = disabled"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetBool(use_penalty_contact,"use_penalty_contact"); ed.SetToolTipText("1 = enabled, 0 = disabled"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(penalty_stiffness_inplane,"penalty_stiffness_inplane",0,0,1); ed.SetToolTipText("penalty stiffness in rolling plane"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(penalty_stiffness_contact,"penalty_stiffness_contact",0,0,1); ed.SetToolTipText("penalty stiffness in contact direction (contact wheel-ground)"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(contact_damping,"contact_damping",0,0,1); ed.SetToolTipText("linear contact damping coefficient"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetBool(use_contact_condition,"use_contact_condition"); ed.SetToolTipText("consider contact condition of wheel"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetBool(use_friction,"use_friction"); ed.SetToolTipText("friction force inplane against sliding velocity (in case if not sticking)"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(friction_coeff_inplane,"friction_coeff_inplane"); ed.SetToolTipText("friction coefficient in rolling forward/lateral direction (must be combined)"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(rolling_friction,"rolling_friction"); ed.SetToolTipText("velocity proportional friction"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(draw_dim(1),"draw_dim"); ed.SetToolTipText("size of sphere (not used for computation)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"force_dim"); ed.SetToolTipText("force scaling factor in [m/N] (not used for computation)"); edc.TreeAdd("Graphics",ed);
  ed.SetInt(elements(1),"element_number",1,0,1); ed.SetToolTipText("number of constrained element"); edc.TreeAdd("Element",ed);
  ed.SetDouble(data_init(1),"sticking"); ed.SetToolTipText("sticking condition at initial time step (1 if sticking, 0 if not sticking)"); edc.TreeAdd("Initialization",ed);
  ed.SetDouble(data_init(2),"contact"); ed.SetToolTipText("contact condition at initial time step (1 if in contact, 0 if not in contact)"); edc.TreeAdd("Initialization",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
}


int RollingJoint3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.rolling_plane_point",rolling_plane_point, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.rolling_plane_normal",rolling_plane_normal, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.wheel_local_center_point",wheel_local_center_point, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.wheel_local_axis",wheel_local_axis, 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.wheel_radius",wheel_radius, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.Penalty.use_penalty_inplane",use_penalty_inplane_dir, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.Penalty.use_penalty_contact",use_penalty_contact, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.penalty_stiffness_inplane",penalty_stiffness_inplane, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.penalty_stiffness_contact",penalty_stiffness_contact, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.contact_damping",contact_damping, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.Penalty.use_contact_condition",use_contact_condition, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.Penalty.use_friction",use_friction, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.friction_coeff_inplane",friction_coeff_inplane, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.rolling_friction",rolling_friction, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_dim",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.force_dim",draw_dim(2), 1);
  GetElemDataInt(GetMBS(), edc, "Element.element_number",elements(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Initialization.sticking",data_init(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Initialization.contact",data_init(2), 1);
  return 1;
}


int RollingJoint3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RollingJoint3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RollingJoint3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=SlidingPointJoint

void SlidingPointJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(elemind,"elemind",1,0,1); ed.SetToolTipText("Index of the initial sliding body."); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(loccoords(1).X(),loccoords(1).Y(),loccoords(1).Z(),"position_1"); ed.SetToolTipText("Vector from the center of body number 1 (en1) to the sliding point in the local body 1 coordinate system."); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(loccoords(2).X(),loccoords(2).Y(),loccoords(2).Z(),"position_2"); ed.SetToolTipText("Vector from the center of the first body of en2 array to the sliding point in the local body 2 coordinate system."); edc.TreeAdd("Geometry",ed);

  {Vector vv(elements.Length()); for(int i = 1; i <= elements.Length(); i++) {vv(i) = elements(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("Element numbers: [en1,en2_1,en2_2,...,en2_n]."); edc.TreeAdd("Geometry",ed);
}
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
}


int SlidingPointJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Geometry.elemind",elemind, 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.position_1",loccoords(1), 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.position_2",loccoords(2), 1);
  GetElemDataIVector(GetMBS(), edc, "Geometry.element_numbers",elements, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  return 1;
}


int SlidingPointJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SlidingPointJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SlidingPointJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=SlidingPrismaticJoint

void SlidingPrismaticJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector3D(loccoords(1).X(),loccoords(1).Y(),loccoords(1).Z(),"position_1"); ed.SetToolTipText("Vector from the center of body number 1 (en1) to the sliding point in the local body 1 coordinate system."); edc.TreeAdd("Geometry",ed);
  ed.SetVector3D(loccoords(2).X(),loccoords(2).Y(),loccoords(2).Z(),"position_2"); ed.SetToolTipText("Vector from the center of the first body of en2 array to the sliding point in the local body 2 coordinate system."); edc.TreeAdd("Geometry",ed);

  {Vector vv(elements.Length()); for(int i = 1; i <= elements.Length(); i++) {vv(i) = elements(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("Element numbers: [en1,en2_1,en2_2,...,en2_n]."); edc.TreeAdd("Geometry",ed);
}
  ed.SetInt(elemind,"elemind",1,0,1); ed.SetToolTipText("Index of the initial sliding body."); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(k1,"k1",0,0,1); ed.SetToolTipText("Stiffness for rotation about global x - axis."); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(k2,"k2",0,0,1); ed.SetToolTipText("Stiffness for rotation about global y - axis."); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(k3,"k3",0,0,1); ed.SetToolTipText("Stiffness for rotation about global z - axis."); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(d1,"d1",0,0,1); ed.SetToolTipText("Damping of rotation about global x - axis."); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(d2,"d2",0,0,1); ed.SetToolTipText("Damping of rotation about global x - axis."); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(d3,"d3",0,0,1); ed.SetToolTipText("Damping of rotation about global x - axis."); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
}


int SlidingPrismaticJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.position_1",loccoords(1), 1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.position_2",loccoords(2), 1);
  GetElemDataIVector(GetMBS(), edc, "Geometry.element_numbers",elements, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.elemind",elemind, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.k1",k1, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.k2",k2, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.k3",k3, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.d1",d1, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.d2",d2, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.d3",d3, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  return 1;
}


int SlidingPrismaticJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SlidingPrismaticJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SlidingPrismaticJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=MultiNodalSphericalJoint

void MultiNodalSphericalJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(draw_dim(1),"draw_size_center_sphere"); ed.SetToolTipText("drawing dimensions of the sphere indicating the average position of the constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"draw_size_circle_spheres"); ed.SetToolTipText("drawing dimensions of the spheres indicating the positions of the multiple nodes."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(3),"draw_size_line_thickness"); ed.SetToolTipText("thickness of the lines from the outer spheres to the center sphere"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(loccoords(2).X(),loccoords(2).Y(),loccoords(2).Z(),"ground"); ed.SetToolTipText("global (average) position, just used if node_numbers.Length()==0"); edc.TreeAdd("Position2",ed);

  {Vector vv(nodes1.Length()); for(int i = 1; i <= nodes1.Length(); i++) {vv(i) = nodes1(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_numbers"); ed.SetValuesInt(); ed.SetToolTipText("global node numbers of kinematic pair 1"); edc.TreeAdd("Position1",ed);
}

  {Vector vv(nodes2.Length()); for(int i = 1; i <= nodes2.Length(); i++) {vv(i) = nodes2(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_numbers"); ed.SetValuesInt(); ed.SetToolTipText("global node numbers of kinematic pair 2"); edc.TreeAdd("Position2",ed);
}
  edc.TreeDelete("Graphics.draw_size"); 
  edc.TreeDelete("Position1.node_number"); 
  edc.TreeDelete("Position2.node_number"); 
  edc.TreeDelete("Position1.element_number"); 
  edc.TreeDelete("Position2.element_number"); 
  edc.TreeDelete("Position1.position"); 
  edc.TreeDelete("Position2.position"); 
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
}


int MultiNodalSphericalJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size_center_sphere",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size_circle_spheres",draw_dim(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size_line_thickness",draw_dim(3), 1);
  GetElemDataVector3D(GetMBS(), edc, "Position2.ground",loccoords(2), 1);
  GetElemDataIVector(GetMBS(), edc, "Position1.node_numbers",nodes1, 1);
  GetElemDataIVector(GetMBS(), edc, "Position2.node_numbers",nodes2, 1);
  return 1;
}


int MultiNodalSphericalJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int MultiNodalSphericalJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int MultiNodalSphericalJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=FrictionConstraintUniDir

void FrictionConstraintUniDir::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(Fn0,"Fn"); ed.SetToolTipText("constant normal force Fn"); edc.TreeAdd("Physics",ed);
  ed.SetBool(constantNormalF,"use_constant_Fn"); ed.SetToolTipText("if this flag is set, normal force is assumed to be constant"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(velocity_tolerance,"velocity_tolerance"); ed.SetToolTipText("if velocity is below this value, than velocity is assumed to be zero"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fr_coeff_st,"fr_coeff_st"); ed.SetToolTipText("friction coefficient 탎, used when elements are sticking (F = 탎*Fn)"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fr_coeff_kin,"fr_coeff_kin"); ed.SetToolTipText("friction coefficient 탃, used when elements are sliding (F = 탃*Fn)"); edc.TreeAdd("Physics",ed);
}


int FrictionConstraintUniDir::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Physics.Fn",Fn0, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.use_constant_Fn",constantNormalF, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.velocity_tolerance",velocity_tolerance, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.fr_coeff_st",fr_coeff_st, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.fr_coeff_kin",fr_coeff_kin, 1);
  return 1;
}


int FrictionConstraintUniDir::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FrictionConstraintUniDir::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FrictionConstraintUniDir::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Rope3D

void Rope3D::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(initial_length,"rope_length"); ed.SetLocked(1); ed.SetToolTipText("initial length l0 of rope (computed automatically)"); edc.TreeAdd("Geometry",ed);
  edc.TreeDelete("Graphics.draw_size_joint_local_frame"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Geometry.use_joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Position1.element_number"); 
  edc.TreeDelete("Position2.element_number"); 
  edc.TreeDelete("Position1.position"); 
  edc.TreeDelete("Position2.position"); 
  edc.TreeDelete("Position1.node_number"); 
  edc.TreeDelete("Position2.node_number"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 

  {Vector vv(elements.Length()); for(int i = 1; i <= elements.Length(); i++) {vv(i) = elements(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("element numbers of the suspension points"); edc.TreeAdd("Geometry",ed);
}
}


int Rope3D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataIVector(GetMBS(), edc, "Geometry.element_numbers",elements, 1);
  return 1;
}


int Rope3D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rope3D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Rope3D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=FrictionConstraint

void FrictionConstraint::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(Fn0,"normal_force"); ed.SetToolTipText("constant normal force Fn"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(velocity_tolerance,"velocity_tolerance",0,0,1); ed.SetToolTipText("If velocity is below this value, sticking starts, or if 'keep sliding' is active, the transition region is used."); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fr_coeff_st,"fr_coeff_st"); ed.SetToolTipText("static friction coefficient, used to determine the threshold when sliding starts."); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fr_coeff_kin,"fr_coeff_kin"); ed.SetToolTipText("kinematic friction coefficient, used to calculate the constant force during sliding phase."); edc.TreeAdd("Physics",ed);
  ed.SetBool(keep_sliding,"keep_sliding"); ed.SetToolTipText("The constraint will never go to modus 'stick'."); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(spring_stiffness,"spring_stiffness"); ed.SetToolTipText("spring stiffness c, only used during sticking phase!"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient d for viscous damping, only used during sticking phase!"); edc.TreeAdd("Physics.Penalty",ed);
}


int FrictionConstraint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Physics.normal_force",Fn0, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.velocity_tolerance",velocity_tolerance, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.fr_coeff_st",fr_coeff_st, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.fr_coeff_kin",fr_coeff_kin, 1);
  GetElemDataBool(GetMBS(), edc, "Physics.keep_sliding",keep_sliding, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.spring_stiffness",spring_stiffness, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  return 1;
}


int FrictionConstraint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FrictionConstraint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FrictionConstraint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Contact1D

void Contact1D::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(lcoord1,"local_coordinate"); ed.SetToolTipText("Local coordinate of element 1 to be constrained"); edc.TreeAdd("Coordinate1",ed);
  ed.SetInt(lcoord2,"local_coordinate"); ed.SetToolTipText("Local coordinate of element 2 to be constrained (not used if ground constraint)"); edc.TreeAdd("Coordinate2",ed);
  ed.SetDouble(lpos1,"position"); ed.SetToolTipText("Local position at which contact occurs"); edc.TreeAdd("Coordinate1",ed);
  ed.SetDouble(lpos2,"position"); ed.SetToolTipText("Local (or global if ground) position at which contact occurs"); edc.TreeAdd("Coordinate2",ed);
  ed.SetDouble(contact_direction,"direction"); ed.SetToolTipText("Direction of the contact: +1 if the first body is on top, or else -1"); edc.TreeAdd("Physics",ed);
  ed.SetInt(mode,"mode"); ed.SetToolTipText("mode of computation"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("Drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(spring_stiffness,"spring_stiffness"); ed.SetToolTipText("spring stiffness c"); edc.TreeAdd("Physics.Mode1",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient d for viscous damping"); edc.TreeAdd("Physics.Mode1",ed);
}


int Contact1D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Coordinate1.local_coordinate",lcoord1, 1);
  GetElemDataInt(GetMBS(), edc, "Coordinate2.local_coordinate",lcoord2, 1);
  GetElemDataDouble(GetMBS(), edc, "Coordinate1.position",lpos1, 1);
  GetElemDataDouble(GetMBS(), edc, "Coordinate2.position",lpos2, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.direction",contact_direction, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.mode",mode, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Mode1.spring_stiffness",spring_stiffness, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Mode1.damping",damping_coeff, 1);
  return 1;
}


int Contact1D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Contact1D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Contact1D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=BaseBodyJoint

void BaseBodyJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;

  {Vector vv(dirRot.Length()); for(int i = 1; i <= dirRot.Length(); i++) {vv(i) = dirRot(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_rotations"); ed.SetValuesInt(); ed.SetToolTipText("[angle about x axis,angle about y axis,angle about z axis]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"); edc.TreeAdd("Physics.Lagrange",ed);
}

  Matrix penaltyStiffness_1(penaltyStiffness);
  ed.SetMatrix(penaltyStiffness_1.GetMatPtr(),penaltyStiffness_1.Getrows(),penaltyStiffness_1.Getcols(),"stiffness_matrix"); ed.SetToolTipText("3x3 matrix with stiffness parameters"); edc.TreeAdd("Physics.Penalty",ed);


  Matrix penaltyDamping_1(penaltyDamping);
  ed.SetMatrix(penaltyDamping_1.GetMatPtr(),penaltyDamping_1.Getrows(),penaltyDamping_1.Getcols(),"damping_matrix"); ed.SetToolTipText("3x3 matrix with damping parameters"); edc.TreeAdd("Physics.Penalty",ed);


  Matrix penaltyStiffnessRot_1(penaltyStiffnessRot);
  ed.SetMatrix(penaltyStiffnessRot_1.GetMatPtr(),penaltyStiffnessRot_1.Getrows(),penaltyStiffnessRot_1.Getcols(),"stiffness_matrix_rotation"); ed.SetToolTipText("3x3 matrix with stiffness parameters for rotation"); edc.TreeAdd("Physics.Penalty",ed);


  Matrix penaltyDampingRot_1(penaltyDampingRot);
  ed.SetMatrix(penaltyDampingRot_1.GetMatPtr(),penaltyDampingRot_1.Getrows(),penaltyDampingRot_1.Getcols(),"damping_matrix_rotation"); ed.SetToolTipText("3x3 matrix with damping parameters for rotation"); edc.TreeAdd("Physics.Penalty",ed);

  edc.TreeDelete("Geometry.joint_local_frame"); 

  Matrix JA0i_1(JA0i);
  ed.SetMatrix(JA0i_1.GetMatPtr(),JA0i_1.Getrows(),JA0i_1.Getcols(),"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);

  ed.SetVector3D(bryantAngles.X(),bryantAngles.Y(),bryantAngles.Z(),"joint_local_frame_in_bryant_angles"); ed.SetToolTipText("Prerotate joint coordinate system w.r.t. local coordinate system of body 1 [phi x, phi y, phi z]. Rot. sequence: JA0i=A(phi z)A(phi y)A(phi x)"); edc.TreeAdd("Geometry",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Geometry.use_joint_local_frame"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Physics.Penalty.damping"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  edc.TreeDelete("Position1.node_number"); 
  edc.TreeDelete("Position2.node_number"); 
  edc.TreeDelete("Graphics.draw_size"); 
  ed.SetDouble(draw_cone_size,"cone_size"); ed.SetToolTipText("cone size for standard joint drawing"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"color_body1"); ed.SetToolTipText("[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_ext.X(),col_ext.Y(),col_ext.Z(),"color_body2"); ed.SetToolTipText("[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
}


int BaseBodyJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataIVector(GetMBS(), edc, "Physics.Lagrange.constrained_rotations",dirRot, 1);

  Matrix penaltyStiffness_1;
  GetElemDataMatrix(GetMBS(), edc, "Physics.Penalty.stiffness_matrix",penaltyStiffness_1, 1);
  penaltyStiffness = Matrix3D(penaltyStiffness_1);

  Matrix penaltyDamping_1;
  GetElemDataMatrix(GetMBS(), edc, "Physics.Penalty.damping_matrix",penaltyDamping_1, 1);
  penaltyDamping = Matrix3D(penaltyDamping_1);

  Matrix penaltyStiffnessRot_1;
  GetElemDataMatrix(GetMBS(), edc, "Physics.Penalty.stiffness_matrix_rotation",penaltyStiffnessRot_1, 1);
  penaltyStiffnessRot = Matrix3D(penaltyStiffnessRot_1);

  Matrix penaltyDampingRot_1;
  GetElemDataMatrix(GetMBS(), edc, "Physics.Penalty.damping_matrix_rotation",penaltyDampingRot_1, 1);
  penaltyDampingRot = Matrix3D(penaltyDampingRot_1);
  GetElemDataVector3D(GetMBS(), edc, "Geometry.joint_local_frame_in_bryant_angles",bryantAngles, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.cone_size",draw_cone_size, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body1",col, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body2",col_ext, 1);
  return 1;
}


int BaseBodyJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int BaseBodyJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int BaseBodyJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RevoluteJoint

void RevoluteJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained directions cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}

  {Vector vv(dirRot.Length()); for(int i = 1; i <= dirRot.Length(); i++) {vv(i) = dirRot(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_rotations"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained rotations cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetVector3D(rotAxis.X(),rotAxis.Y(),rotAxis.Z(),"rotation_axis"); ed.SetToolTipText("local rotation axis w.r.t body 1 coordinate system"); edc.TreeAdd("Physics",ed);
  ed.SetBool(standard_joint_drawing,"standard_joint_drawing"); ed.SetToolTipText("flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_diameter,"diameter"); ed.SetToolTipText("diameter of the revolute joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_axis_length,"axis_length"); ed.SetToolTipText("axis length of the revolute joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(spring_stiffness,"stiffness",0,0,1); ed.SetToolTipText("stiffness parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
}


int RevoluteJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Physics.rotation_axis",rotAxis, 1);
  GetElemDataBool(GetMBS(), edc, "Graphics.standard_joint_drawing",standard_joint_drawing, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.diameter",draw_diameter, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.axis_length",draw_axis_length, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.stiffness",spring_stiffness, 1);
  return 1;
}


int RevoluteJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RevoluteJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RevoluteJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=PrismaticJoint

void PrismaticJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained directions cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}

  {Vector vv(dirRot.Length()); for(int i = 1; i <= dirRot.Length(); i++) {vv(i) = dirRot(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_rotations"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained rotations cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetVector3D(slidingDirection.X(),slidingDirection.Y(),slidingDirection.Z(),"sliding_direction"); ed.SetToolTipText("local sliding direction w.r.t body 1 coordinate system"); edc.TreeAdd("Physics",ed);
  ed.SetBool(standard_joint_drawing,"standard_joint_drawing"); ed.SetToolTipText("flag for drawing mode; 1 == draw constraint nicely; 0 == show constrained directions and rotations;"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_length,"rail_length"); ed.SetToolTipText("length of the prismatic joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(draw_cube_size.X(),draw_cube_size.Y(),draw_cube_size.Z(),"joint_cube_size"); ed.SetToolTipText("cube dimension of prismatic joint (for drawing); [lx (in sl. dir.),ly (normal to sl. dir.),lz (normal to sl. dir.)]"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(spring_stiffness,"stiffness",0,0,1); ed.SetToolTipText("stiffness parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
}


int PrismaticJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Physics.sliding_direction",slidingDirection, 1);
  GetElemDataBool(GetMBS(), edc, "Graphics.standard_joint_drawing",standard_joint_drawing, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.rail_length",draw_length, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.joint_cube_size",draw_cube_size, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.stiffness",spring_stiffness, 1);
  return 1;
}


int PrismaticJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int PrismaticJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int PrismaticJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=UniversalJoint

void UniversalJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 
  edc.TreeDelete("Graphics.color_body1"); 
  edc.TreeDelete("Graphics.color_body2"); 
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"color_body1"); ed.SetToolTipText("[red, green, blue] color of the hinge connected to the first body, range = 0..1"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_ext.X(),col_ext.Y(),col_ext.Z(),"color_body2"); ed.SetToolTipText("[red, green, blue] color of the hinge connected to the first body, range = 0..1"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_t.X(),col_t.Y(),col_t.Z(),"color_cross"); ed.SetToolTipText("[red, green, blue] color of the cross shaft"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Position1.position"); 
  edc.TreeDelete("Position2.position"); 
  ed.SetVector3D(loccoords(1).X(),loccoords(1).Y(),loccoords(1).Z(),"position"); ed.SetToolTipText("local position"); edc.TreeAdd("Position1",ed);
  ed.SetVector3D(loccoords(2).X(),loccoords(2).Y(),loccoords(2).Z(),"position"); ed.SetToolTipText("local or global (if element_number == 0) position"); edc.TreeAdd("Position2",ed);
  ed.SetVector3D(axis_1.X(),axis_1.Y(),axis_1.Z(),"axis"); ed.SetToolTipText("the axis of the cross connected to body 1 in local coordinates"); edc.TreeAdd("Position1",ed);
  ed.SetVector3D(axis_2.X(),axis_2.Y(),axis_2.Z(),"axis"); ed.SetToolTipText("the axis of the cross connected to body 2 in local coordinates"); edc.TreeAdd("Position2",ed);
  ed.SetDouble(draw_length,"draw_length"); ed.SetToolTipText("length of the universal joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_width,"draw_width"); ed.SetToolTipText("width of the universal joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(draw_direction_1.X(),draw_direction_1.Y(),draw_direction_1.Z(),"draw_direction_1"); ed.SetToolTipText("direction from body 1 to joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(draw_direction_2.X(),draw_direction_2.Y(),draw_direction_2.Z(),"draw_direction_2"); ed.SetToolTipText("direction from body 2 to joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Graphics.draw_size_joint_local_frame"); 
  edc.TreeDelete("Physics.Penalty.damping"); 
  edc.TreeDelete("Physics.Penalty.stiffness"); 
  edc.TreeDelete("Graphics.cone_size"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
}


int UniversalJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"color_body1"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"color_body2"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"cone_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body1",col, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body2",col_ext, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_cross",col_t, 1);
  GetElemDataVector3D(GetMBS(), edc, "Position1.position",loccoords(1), 1);
  GetElemDataVector3D(GetMBS(), edc, "Position2.position",loccoords(2), 1);
  GetElemDataVector3D(GetMBS(), edc, "Position1.axis",axis_1, 1);
  GetElemDataVector3D(GetMBS(), edc, "Position2.axis",axis_2, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_length",draw_length, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_width",draw_width, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.draw_direction_1",draw_direction_1, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.draw_direction_2",draw_direction_2, 1);
  return 1;
}


int UniversalJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int UniversalJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int UniversalJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RigidJoint

void RigidJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained directions cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}

  {Vector vv(dirRot.Length()); for(int i = 1; i <= dirRot.Length(); i++) {vv(i) = dirRot(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_rotations"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained rotations cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetBool(standard_joint_drawing,"standard_joint_drawing"); ed.SetToolTipText("flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dimension,"cube_length"); ed.SetToolTipText("rigid joint dimension (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(spring_stiffness,"stiffness",0,0,1); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
}


int RigidJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "Graphics.standard_joint_drawing",standard_joint_drawing, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.cube_length",draw_dimension, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.stiffness",spring_stiffness, 1);
  return 1;
}


int RigidJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RigidJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RigidJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=CylindricalJoint

void CylindricalJoint::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained directions cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}

  {Vector vv(dirRot.Length()); for(int i = 1; i <= dirRot.Length(); i++) {vv(i) = dirRot(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_rotations"); ed.SetValuesInt(); ed.SetLocked(1); ed.SetToolTipText("constrained rotations cannot be changed"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetVector3D(rotSlideAxis.X(),rotSlideAxis.Y(),rotSlideAxis.Z(),"rotation_sliding_axis"); ed.SetToolTipText("local rotation/sliding axis w.r.t body 1 coordinate system"); edc.TreeAdd("Physics",ed);
  ed.SetBool(standard_joint_drawing,"standard_joint_drawing"); ed.SetToolTipText("flag for drawing mode; 1 == draw constraint element; 0 == show constrained directions and rotations;"); edc.TreeAdd("Graphics",ed);
  ed.SetVector2D(draw_cylinder_size.X(),draw_cylinder_size.Y(),"joint_cylinder_size"); ed.SetToolTipText("cylinder dimension of cylindrical joint (for drawing); [lx (cyl. length, in sl. dir.),d (cylinder diameter)]"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_axis_length,"axis_length"); ed.SetToolTipText("axis length of the revolute joint (for drawing)"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetDouble(spring_stiffness,"stiffness",0,0,1); ed.SetToolTipText("stiffness parameter used for translation and rotation"); edc.TreeAdd("Physics.Penalty",ed);
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
}


int CylindricalJoint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataVector3D(GetMBS(), edc, "Physics.rotation_sliding_axis",rotSlideAxis, 1);
  GetElemDataBool(GetMBS(), edc, "Graphics.standard_joint_drawing",standard_joint_drawing, 1);
  GetElemDataVector2D(GetMBS(), edc, "Graphics.joint_cylinder_size",draw_cylinder_size, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.axis_length",draw_axis_length, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.stiffness",spring_stiffness, 1);
  return 1;
}


int CylindricalJoint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int CylindricalJoint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int CylindricalJoint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=SpringDamperActuator

void SpringDamperActuator::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(l0,"spring_length"); ed.SetToolTipText("length of the spring in the initial configuration"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fa,"actor_force"); ed.SetToolTipText("constant force acting on the spring"); edc.TreeAdd("Physics",ed);
  ed.SetInt(forcemode,"forcemode"); ed.SetToolTipText("defines how the spring and damper force is computed: 0..constant coefficient, 1..MathFunction (stiffness and damping), 2..IOElementDataModifier, 3..MathFunction (spring force and damping force)"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"color_body1"); ed.SetToolTipText("[red, green, blue] first color of constraint (spring), range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_ext.X(),col_ext.Y(),col_ext.Z(),"color_body2"); ed.SetToolTipText("[red, green, blue] second color of constraint (damper), range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(spring_stiffness,"spring_stiffness",0,0,1); ed.SetToolTipText("stiffness coefficient of the linear spring. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  edc.TreeDelete("Physics.Penalty.damping"); 
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient for viscous damping. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Geometry.use_joint_local_frame"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Graphics.draw_size"); 
  edc.TreeDelete("Graphics.draw_size_joint_local_frame"); 
  ed.SetDouble(draw_dim(1),"spring_diameter",0,0,1); ed.SetToolTipText("spring diameter used for drawing only."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"spring_coils",0,0,1); ed.SetToolTipText("spring coils used for drawing. If set to 0, then a cylinder with the value 'spring_diameter' as diameter is shown instead of the coils."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(spring_res,"spring_resolution",1,10); ed.SetToolTipText("spring resolution used for drawing (very coarse = 1, very smooth = 10)."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(3),"damper_diameter",0,0,1); ed.SetToolTipText("damper diameter used for drawing only. If set to 0, then the damper is not shown. It's recommended to choose the value smaller then the spring diameter."); edc.TreeAdd("Graphics",ed);
}


int SpringDamperActuator::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Physics.spring_length",l0, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.actor_force",fa, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.forcemode",forcemode, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body1",col, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body2",col_ext, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.spring_stiffness",spring_stiffness, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_diameter",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_coils",draw_dim(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_resolution",spring_res, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.damper_diameter",draw_dim(3), 1);
  return 1;
}


int SpringDamperActuator::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SpringDamperActuator::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SpringDamperActuator::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RigidLink

void RigidLink::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(distancemode,"distancemode"); ed.SetToolTipText("defines the distance: 0..constant distance, 1..MathFunction, 2..IOElementDataModifier"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(distance,"link_length"); ed.SetToolTipText("constant distance is used, when distancemode = 0"); edc.TreeAdd("Physics.Constant",ed);
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"color_body1"); ed.SetToolTipText("[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_ext.X(),col_ext.Y(),col_ext.Z(),"color_body2"); ed.SetToolTipText("[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Geometry.use_joint_local_frame"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Graphics.draw_size"); 
  edc.TreeDelete("Graphics.draw_size_joint_local_frame"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  edc.TreeDelete("Physics.Penalty.damping"); 
  ed.SetDouble(draw_dim(1),"cylinder1_diameter",0,0,1); ed.SetToolTipText("cylinder one diameter (drawing only)."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"cylinder2_diameter",0,0,1); ed.SetToolTipText("cylinder two diameter (drawing only). Only used if distance not constant = distancemode 1 or 2."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(3),"cylinder1_length",0,0,1); ed.SetToolTipText("cylinder one length (drawing only). Only used if distance not constant = distancemode 1 or 2."); edc.TreeAdd("Graphics",ed);
}


int RigidLink::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Physics.distancemode",distancemode, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Constant.link_length",distance, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body1",col, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body2",col_ext, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.cylinder1_diameter",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.cylinder2_diameter",draw_dim(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.cylinder1_length",draw_dim(3), 1);
  return 1;
}


int RigidLink::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RigidLink::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RigidLink::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RotatorySpringDamperActuator

void RotatorySpringDamperActuator::GetElementDataAuto(ElementDataContainer& edc)
{
  BaseBodyJoint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetDouble(phi0,"spring_angle_offset"); ed.SetToolTipText("spring angle offset is used if constant_spring_angle_offset is enabled. A positive offset equates a positve angle about the rotation axis."); edc.TreeAdd("Physics",ed);
  ed.SetDouble(ma,"actuator_torque"); ed.SetToolTipText("constant torque of an actuator. A positive torque is acting about the rotation axis in a positive sense."); edc.TreeAdd("Physics",ed);
  ed.SetVector3D(glob_rot_axis.X(),glob_rot_axis.Y(),glob_rot_axis.Z(),"rotation_axis"); ed.SetToolTipText("local axis of rotation w.r.t. body 1 coordinate system in inital configuration"); edc.TreeAdd("Physics",ed);
  ed.SetInt(forcemode,"forcemode"); ed.SetToolTipText("defines how the spring and damper moment is computed: 0..constant coefficient, 1..MathFunction, 2..IOElementDataModifier"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(spring_stiffness,"spring_stiffness",0,0,1); ed.SetToolTipText("stiffness parameter of the rotatory spring. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient for viscous damping. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Geometry.use_joint_local_frame"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Graphics.draw_size"); 
  edc.TreeDelete("Graphics.draw_size_joint_local_frame"); 
  edc.TreeDelete("Physics.Lagrange.constrained_rotations"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix"); 
  edc.TreeDelete("Physics.Penalty.stiffness_matrix_rotation"); 
  edc.TreeDelete("Physics.Penalty.damping_matrix_rotation"); 
  edc.TreeDelete("Geometry.joint_local_frame_in_bryant_angles"); 
  edc.TreeDelete("Graphics.cone_size"); 
  ed.SetDouble(draw_dim(1),"spring_size",0,0,1); ed.SetToolTipText("radius of torsional spring. This parameter is used for drawing only."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"windings",0,0,1); ed.SetToolTipText("number of windings of torsional spring. This parameter is used for drawing only."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(spring_res,"spring_resolution",1,10); ed.SetToolTipText("spring resolution used for drawing (very coarse = 1, very smooth = 10)."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(3),"axis_radius",0,0,1); ed.SetToolTipText("radius of torsional spring axis (cylinder). This parameter is used for drawing only."); edc.TreeAdd("Graphics",ed);
}


int RotatorySpringDamperActuator::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"draw_size_joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_rotations"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"stiffness_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"damping_matrix_rotation"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"joint_local_frame_in_bryant_angles"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"cone_size"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  BaseBodyJoint::SetElementDataAuto(edc);
  GetElemDataDouble(GetMBS(), edc, "Physics.spring_angle_offset",phi0, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.actuator_torque",ma, 1);
  GetElemDataVector3D(GetMBS(), edc, "Physics.rotation_axis",glob_rot_axis, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.forcemode",forcemode, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.spring_stiffness",spring_stiffness, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_size",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.windings",draw_dim(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_resolution",spring_res, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.axis_radius",draw_dim(3), 1);
  return 1;
}


int RotatorySpringDamperActuator::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BaseBodyJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RotatorySpringDamperActuator::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BaseBodyJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RotatorySpringDamperActuator::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=SpringDamperActuator2D

void SpringDamperActuator2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(elements(1),"element_number",1,0,1); ed.SetToolTipText("Number of constrained element"); edc.TreeAdd("Position1",ed);
  ed.SetInt(elements(2),"element_number",0,0,1); ed.SetToolTipText("Number of constrained element (0 if ground joint)"); edc.TreeAdd("Position2",ed);
  ed.SetVector2D(loccoords(1).X(),loccoords(1).Y(),"position"); ed.SetToolTipText("local position 1"); edc.TreeAdd("Position1",ed);
  ed.SetVector2D(loccoords(2).X(),loccoords(2).Y(),"position"); ed.SetToolTipText("local or global position 2"); edc.TreeAdd("Position2",ed);
  edc.TreeDelete("Graphics.RGB_color"); 
  ed.SetVector3D(col.X(),col.Y(),col.Z(),"color_body1"); ed.SetToolTipText("[red, green, blue] first color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(col_ext.X(),col_ext.Y(),col_ext.Z(),"color_body2"); ed.SetToolTipText("[red, green, blue] second color of constraint, range = 0..1, use default color:[-1,-1,-1]"); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  ed.SetDouble(l0,"spring_length"); ed.SetToolTipText("length of the spring in the initial configuration"); edc.TreeAdd("Physics",ed);
  ed.SetDouble(fa,"actor_force"); ed.SetToolTipText("constant force acting on the spring"); edc.TreeAdd("Physics",ed);
  ed.SetInt(forcemode,"forcemode"); ed.SetToolTipText("defines how the spring and damper force is computed: 0..constant coefficient, 1..MathFunction, 2..IOElementDataModifier"); edc.TreeAdd("Physics",ed);
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  ed.SetDouble(spring_stiffness,"spring_stiffness",0,0,1); ed.SetToolTipText("stiffness coefficient of the linear spring. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient for viscous damping. Only used if forcemode is 0."); edc.TreeAdd("Physics.Linear",ed);
  ed.SetDouble(draw_dim(1),"spring_diameter",0,0,1); ed.SetToolTipText("spring diameter used for drawing only."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(2),"spring_coils",0,0,1); ed.SetToolTipText("spring coils used for drawing. If set to 0, then a cylinder with the value 'spring_diameter' as diameter is shown instead of the coils."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(spring_res,"spring_resolution",1,10); ed.SetToolTipText("spring resolution used for drawing (very coarse = 1, very smooth = 10)."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(3),"damper_diameter",0,0,1); ed.SetToolTipText("damper diameter used for drawing only. If set to 0, then the damper is not shown. It's recommended to choose the value smaller then the spring diameter."); edc.TreeAdd("Graphics",ed);
}


int SpringDamperActuator2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Position1.element_number",elements(1), 1);
  GetElemDataInt(GetMBS(), edc, "Position2.element_number",elements(2), 1);
  GetElemDataVector2D(GetMBS(), edc, "Position1.position",loccoords(1), 1);
  GetElemDataVector2D(GetMBS(), edc, "Position2.position",loccoords(2), 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body1",col, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.color_body2",col_ext, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.spring_length",l0, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.actor_force",fa, 1);
  GetElemDataInt(GetMBS(), edc, "Physics.forcemode",forcemode, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.spring_stiffness",spring_stiffness, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Linear.damping",damping_coeff, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_diameter",draw_dim(1), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_coils",draw_dim(2), 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.spring_resolution",spring_res, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.damper_diameter",draw_dim(3), 1);
  return 1;
}


int SpringDamperActuator2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SpringDamperActuator2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SpringDamperActuator2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=PointJoint2D

void PointJoint2D::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(elements(1),"element_number",1,0,1); ed.SetToolTipText("Number of constrained element"); edc.TreeAdd("Position1",ed);
  ed.SetInt(elements(2),"element_number",0,0,1); ed.SetToolTipText("Number of constrained element (0 if ground joint)"); edc.TreeAdd("Position2",ed);
  ed.SetVector2D(loccoords(1).X(),loccoords(1).Y(),"position"); ed.SetToolTipText("local position 1"); edc.TreeAdd("Position1",ed);
  ed.SetVector2D(loccoords(2).X(),loccoords(2).Y(),"position"); ed.SetToolTipText("local or global position 2"); edc.TreeAdd("Position2",ed);
  ed.SetVector2D(spring_stiffness2.X(),spring_stiffness2.Y(),"spring_stiffness_vector"); ed.SetToolTipText("penalty stiffness parameter [kx,ky]. Just used if scalar spring_stiffness == 0, otherwise kx=ky=spring_stiffness"); edc.TreeAdd("Physics.Penalty",ed);

  {Vector vv(dir.Length()); for(int i = 1; i <= dir.Length(); i++) {vv(i) = dir(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"constrained_directions"); ed.SetValuesInt(); ed.SetToolTipText("[x,y]...(1 = constrained, 0 = free), can be defined as local or global directions (see Geometry)"); edc.TreeAdd("Physics.Lagrange",ed);
}
  ed.SetDouble(damping_coeff,"damping",0,0,1); ed.SetToolTipText("damping coefficient for viscous damping (F = d*v), applied in all constrained directions"); edc.TreeAdd("Physics.Penalty",ed);
  ed.SetBool(stiffness_in_joint_local_frame,"use_joint_local_frame"); ed.SetToolTipText("Use a special joint local frame"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(phi_z,"joint_local_frame"); ed.SetToolTipText("Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body 1 with angle phi_z about the z axis. Just used if use_joint_local_frame == 1"); edc.TreeAdd("Geometry",ed);
  ed.SetDouble(draw_local_frame_size,"draw_size_joint_local_frame"); ed.SetToolTipText("drawing dimensions of joint local frame. If set to -1, than global_draw_scalar_size is used. If set to 0, than no joint local frame is drawn."); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(draw_dim(1),"draw_size"); ed.SetToolTipText("drawing dimensions of constraint. If set to -1, than global_draw_scalar_size is used."); edc.TreeAdd("Graphics",ed);
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  ed.SetBool(use_local_coordinate_system,"use_local_coordinate_system"); ed.SetToolTipText("0=use global coordinates, 1=use local coordinate system of Body 1"); edc.TreeAdd("Geometry",ed);
}


int PointJoint2D::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "Position1.element_number",elements(1), 1);
  GetElemDataInt(GetMBS(), edc, "Position2.element_number",elements(2), 1);
  GetElemDataVector2D(GetMBS(), edc, "Position1.position",loccoords(1), 1);
  GetElemDataVector2D(GetMBS(), edc, "Position2.position",loccoords(2), 1);
  GetElemDataVector2D(GetMBS(), edc, "Physics.Penalty.spring_stiffness_vector",spring_stiffness2, 1);
  GetElemDataIVector(GetMBS(), edc, "Physics.Lagrange.constrained_directions",dir, 1);
  GetElemDataDouble(GetMBS(), edc, "Physics.Penalty.damping",damping_coeff, 1);
  GetElemDataBool(GetMBS(), edc, "Geometry.use_joint_local_frame",stiffness_in_joint_local_frame, 1);
  GetElemDataDouble(GetMBS(), edc, "Geometry.joint_local_frame",phi_z, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size_joint_local_frame",draw_local_frame_size, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.draw_size",draw_dim(1), 1);
  GetElemDataBool(GetMBS(), edc, "Geometry.use_local_coordinate_system",use_local_coordinate_system, 1);
  return 1;
}


int PointJoint2D::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int PointJoint2D::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int PointJoint2D::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=AverageConstraint

void AverageConstraint::GetElementDataAuto(ElementDataContainer& edc)
{
  BasePointJoint::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Position1.element_number"); 
  edc.TreeDelete("Position2.element_number"); 
  edc.TreeDelete("Position1.position"); 
  edc.TreeDelete("Position2.position"); 
  edc.TreeDelete("Position1.node_number"); 
  edc.TreeDelete("Position2.node_number"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness_vector"); 
  edc.TreeDelete("Physics.Penalty.damping"); 
  edc.TreeDelete("Physics.Lagrange.constrained_directions"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Geometry.joint_local_frame"); 

  {Vector vv(nodes.Length()); for(int i = 1; i <= nodes.Length(); i++) {vv(i) = nodes(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("(local) node numbers of the kinematic pair 1"); edc.TreeAdd("Position1",ed);
}

  {Vector vv(nodes2.Length()); for(int i = 1; i <= nodes2.Length(); i++) {vv(i) = nodes2(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"node_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("(local) node numbers of the kinematic pair 2"); edc.TreeAdd("Position2",ed);
}
  ed.SetInt(use_local_coordinate_system,"use_local_coordinate_system"); ed.SetToolTipText("0 ...use global coordinates, i=1,2,..n ...use local coordinate system of body i, evaluated at local position (0,0,0)"); edc.TreeAdd("Geometry",ed);

  Matrix JA0i_1(JA0i);
  ed.SetMatrix(JA0i_1.GetMatPtr(),JA0i_1.Getrows(),JA0i_1.Getcols(),"joint_local_frame"); ed.SetToolTipText("Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body i. Just used if use_joint_local_frame == 1"); edc.TreeAdd("Geometry",ed);


  {Vector vv(elements.Length()); for(int i = 1; i <= elements.Length(); i++) {vv(i) = elements(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("element numbers of the kinematic pair 1"); edc.TreeAdd("Position1",ed);
}

  {Vector vv(elements2.Length()); for(int i = 1; i <= elements2.Length(); i++) {vv(i) = elements2(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("element numbers of the kinematic pair 2"); edc.TreeAdd("Position2",ed);
}
  ed.SetVector(weights1.GetVecPtr(),weights1.Length(),"weights"); ed.SetVariableLength(); ed.SetToolTipText("weights of the points of kinematic pair 1"); edc.TreeAdd("Position1",ed);
  ed.SetVector(weights2.GetVecPtr(),weights2.Length(),"weights"); ed.SetVariableLength(); ed.SetToolTipText("weights of the points of kinematic pair 2"); edc.TreeAdd("Position2",ed);
  ed.SetBool(elements_share_dofs,"elements_share_dofs"); ed.SetToolTipText("check, if not just one, but a set of elements are constrained, and some of those elements share the same dofs"); edc.TreeAdd("Geometry",ed);
  ed.SetVector(stiffnesses.GetVecPtr(),stiffnesses.Length(),"weights"); ed.SetVariableLength(); ed.SetToolTipText("each moment can be weighted with a stiffness, e.g. c1*fx + c2*fy"); edc.TreeAdd("Physics.Penalty",ed);
}


int AverageConstraint::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"element_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"position"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position1",ed);}
  {  int dummy=0; ed.SetInt(dummy,"node_number"); ed.SetLocked(1); edc.TreeAdd("Position2",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"spring_stiffness_vector"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"damping"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"constrained_directions"); ed.SetLocked(1); edc.TreeAdd("Physics.Lagrange",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Matrix dummy(3,3); ed.SetMatrix(dummy.GetMatPtr(),3,3,"joint_local_frame"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}


  BasePointJoint::SetElementDataAuto(edc);
  GetElemDataIVector(GetMBS(), edc, "Position1.node_numbers",nodes, 1);
  GetElemDataIVector(GetMBS(), edc, "Position2.node_numbers",nodes2, 1);
  GetElemDataInt(GetMBS(), edc, "Geometry.use_local_coordinate_system",use_local_coordinate_system, 1);

  Matrix JA0i_1;
  GetElemDataMatrix(GetMBS(), edc, "Geometry.joint_local_frame",JA0i_1, 1);
  JA0i = Matrix3D(JA0i_1);
  GetElemDataIVector(GetMBS(), edc, "Position1.element_numbers",elements, 1);
  GetElemDataIVector(GetMBS(), edc, "Position2.element_numbers",elements2, 1);
  GetElemDataVector(GetMBS(), edc, "Position1.weights",weights1, 1);
  GetElemDataVector(GetMBS(), edc, "Position2.weights",weights2, 1);
  GetElemDataBool(GetMBS(), edc, "Geometry.elements_share_dofs",elements_share_dofs, 1);
  GetElemDataVector(GetMBS(), edc, "Physics.Penalty.weights",stiffnesses, 1);
  return 1;
}


int AverageConstraint::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = BasePointJoint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int AverageConstraint::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = BasePointJoint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int AverageConstraint::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Sensor

void Sensor::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  ed.SetInt(sensnum,"sensor_number"); ed.SetLocked(1); ed.SetToolTipText("number of the sensor in the mbs"); edc.TreeAdd("",ed);
  ed.SetText(name.c_str(),"name"); ed.SetToolTipText("name of the sensor for the output files and for the plot tool"); edc.TreeAdd("",ed);
  ed.SetText(GetTypeName().c_str(),"sensor_type"); ed.SetToolTipText("specification of sensor type. Once the sensor is added to the mbs, you MUST NOT change this type anymore!"); edc.TreeAdd("",ed);
}


int Sensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  GetElemDataText(GetMBS(), edc, "name",name, 1);
  GetElemDataText(GetMBS(), edc, "sensor_type",GetTypeName(), 1);
  return 1;
}


int Sensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  return 0;
}


int Sensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  return 0;
}


int Sensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ElementSensor

void ElementSensor::GetElementDataAuto(ElementDataContainer& edc)
{
  Sensor::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(elementNumber,"element_number"); ed.SetToolTipText("number of the element, to which the sensor is applied"); edc.TreeAdd("",ed);
}


int ElementSensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Sensor::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "element_number",elementNumber, 1);
  return 1;
}


int ElementSensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Sensor::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ElementSensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Sensor::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ElementSensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=FieldVariableElementSensor

void FieldVariableElementSensor::GetElementDataAuto(ElementDataContainer& edc)
{
  ElementSensor::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(localNodeNumber,"node_number"); ed.SetToolTipText("local node number. If > 0, then the position of this node is used."); edc.TreeAdd("",ed);
  ed.SetVector3D(localPosition.X(),localPosition.Y(),localPosition.Z(),"local_position"); ed.SetToolTipText("local position at which the field variable is evaluated."); edc.TreeAdd("",ed);
}


int FieldVariableElementSensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  ElementSensor::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "node_number",localNodeNumber, 1);
  GetElemDataVector3D(GetMBS(), edc, "local_position",localPosition, 1);
  return 1;
}


int FieldVariableElementSensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = ElementSensor::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FieldVariableElementSensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = ElementSensor::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int FieldVariableElementSensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=MultipleSensor

void MultipleSensor::GetElementDataAuto(ElementDataContainer& edc)
{
  Sensor::GetElementDataAuto(edc);
  ElementData ed;

  {Vector vv(sensorNrs.Length()); for(int i = 1; i <= sensorNrs.Length(); i++) {vv(i) = sensorNrs(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"sensor_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("number of the sensors, that are used for computation"); edc.TreeAdd("",ed);
}
  ed.SetVector(weights.GetVecPtr(),weights.Length(),"weights"); ed.SetVariableLength(); ed.SetToolTipText("weights for e.g. a weighted sum. This vector must have the same length as sensor_numbers or must be empty!"); edc.TreeAdd("",ed);
}


int MultipleSensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Sensor::SetElementDataAuto(edc);
  GetElemDataIVector(GetMBS(), edc, "sensor_numbers",sensorNrs, 1);
  GetElemDataVector(GetMBS(), edc, "weights",weights, 1);
  return 1;
}


int MultipleSensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Sensor::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int MultipleSensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Sensor::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int MultipleSensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=LoadSensor

void LoadSensor::GetElementDataAuto(ElementDataContainer& edc)
{
  Sensor::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(loadNumber,"load_number"); ed.SetToolTipText("number of the load, to which the sensor is applied"); edc.TreeAdd("",ed);
}


int LoadSensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Sensor::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "load_number",loadNumber, 1);
  return 1;
}


int LoadSensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Sensor::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LoadSensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Sensor::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LoadSensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=SystemSensor

void SystemSensor::GetElementDataAuto(ElementDataContainer& edc)
{
  Sensor::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetText(object_str.c_str(),"object"); ed.SetToolTipText("Object tracked by systemsensor. Is either 'DOF' (global degree of freedom), 'EV' (global eigenvalue), 'jacobians', 'newton_iterations', 'discontinuous_iterations', 'rhs_evaluations', or 'rhs_evaluations_jacobian'"); edc.TreeAdd("",ed);
  ed.SetInt(global_index,"global_index"); ed.SetToolTipText("Number of the global index. Has to be set if (and only if) object is 'DOF' or 'EV'."); edc.TreeAdd("",ed);
}


int SystemSensor::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  Sensor::SetElementDataAuto(edc);
  GetElemDataText(GetMBS(), edc, "object",object_str, 1);
  GetElemDataInt(GetMBS(), edc, "global_index",global_index, 1);
  return 1;
}


int SystemSensor::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Sensor::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SystemSensor::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Sensor::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int SystemSensor::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=InputOutputElement

void InputOutputElement::GetElementDataAuto(ElementDataContainer& edc)
{
  Constraint::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(GetNInputs(),"number_of_inputs"); ed.SetLocked(1); ed.SetToolTipText("number of inputs"); edc.TreeAdd("IOBlock",ed);
  ed.SetInt(GetNOutputs(),"number_of_outputs"); ed.SetLocked(1); ed.SetToolTipText("number of outputs"); edc.TreeAdd("IOBlock",ed);
  ed.SetInt(GetNStates(),"number_of_states"); ed.SetLocked(1); ed.SetToolTipText("number of states"); edc.TreeAdd("IOBlock",ed);

  {Vector vv(inputs.Length()); for(int i = 1; i <= inputs.Length(); i++) {vv(i) = inputs(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"input_element_numbers"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("vector of element(s) or sensor number(s) connected to input, only valid element numbers permitted!"); edc.TreeAdd("IOBlock",ed);
}

  {Vector vv(input_types.Length()); for(int i = 1; i <= input_types.Length(); i++) {vv(i) = input_types(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"input_element_types"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("vector with types of connected inputs; 1=IOElement, 2=Sensor"); edc.TreeAdd("IOBlock",ed);
}

  {Vector vv(input_localnum.Length()); for(int i = 1; i <= input_localnum.Length(); i++) {vv(i) = input_localnum(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"input_local_number"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("vector with i-th number of output of previous IOelement connected to this element"); edc.TreeAdd("IOBlock",ed);
}
  ed.SetVector2D(ref_pos.X(),ref_pos.Y(),"position"); ed.SetToolTipText("reference drawing position"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(draw_dim.X(),draw_dim.Y(),draw_dim.Z(),"draw_size"); ed.SetToolTipText("draw size"); edc.TreeAdd("Graphics",ed);
  ed.SetDouble(rotation,"rotation"); ed.SetToolTipText("rotation: 1==90�, 2==180�, 3==270�, 4=360�"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(colbackground.X(),colbackground.Y(),colbackground.Z(),"background_color"); ed.SetToolTipText("background color; -1=transparent"); edc.TreeAdd("Graphics",ed);
  ed.SetVector3D(colforeground.X(),colforeground.Y(),colforeground.Z(),"foreground_color"); ed.SetToolTipText("foreground color"); edc.TreeAdd("Graphics",ed);

  {Vector vv(input_nodes_num.Length()); for(int i = 1; i <= input_nodes_num.Length(); i++) {vv(i) = input_nodes_num(i);}
  ed.SetVector(vv.GetVecPtr(),vv.Length(),"input_nodes_num"); ed.SetValuesInt(); ed.SetVariableLength(); ed.SetToolTipText("number of input of drawing position \"input_nodes\""); edc.TreeAdd("Graphics",ed);
}
  edc.TreeDelete("Physics.use_penalty_formulation"); 
  edc.TreeDelete("Physics.Penalty.spring_stiffness"); 
  edc.TreeDelete("Geometry.use_local_coordinate_system"); 
  edc.TreeDelete("Graphics.RGB_color"); 
}


int InputOutputElement::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  int dummy=0; ed.SetInt(dummy,"use_penalty_formulation"); ed.SetLocked(1); edc.TreeAdd("Physics",ed);}
  {  double dummy= 0.; ed.SetDouble(dummy,"spring_stiffness"); ed.SetLocked(1); edc.TreeAdd("Physics.Penalty",ed);}
  {  int dummy=0; ed.SetInt(dummy,"use_local_coordinate_system"); ed.SetLocked(1); edc.TreeAdd("Geometry",ed);}
  {  Vector3D dummy(0.); ed.SetVector3D(dummy.X(),dummy.Y(),dummy.Z(),"RGB_color"); ed.SetLocked(1); edc.TreeAdd("Graphics",ed);}


  Constraint::SetElementDataAuto(edc);
  GetElemDataIVector(GetMBS(), edc, "IOBlock.input_element_numbers",inputs, 1);
  GetElemDataIVector(GetMBS(), edc, "IOBlock.input_element_types",input_types, 1);
  GetElemDataIVector(GetMBS(), edc, "IOBlock.input_local_number",input_localnum, 1);
  GetElemDataVector2D(GetMBS(), edc, "Graphics.position",ref_pos, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.draw_size",draw_dim, 1);
  GetElemDataDouble(GetMBS(), edc, "Graphics.rotation",rotation, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.background_color",colbackground, 1);
  GetElemDataVector3D(GetMBS(), edc, "Graphics.foreground_color",colforeground, 1);
  GetElemDataIVector(GetMBS(), edc, "Graphics.input_nodes_num",input_nodes_num, 1);
  return 1;
}


int InputOutputElement::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = Constraint::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int InputOutputElement::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = Constraint::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int InputOutputElement::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=InputOutputElementDiscrete

void InputOutputElementDiscrete::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int InputOutputElementDiscrete::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int InputOutputElementDiscrete::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int InputOutputElementDiscrete::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int InputOutputElementDiscrete::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=ZTransferFunction

void ZTransferFunction::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElementDiscrete::GetElementDataAuto(edc);
  ElementData ed;
}


int ZTransferFunction::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElementDiscrete::SetElementDataAuto(edc);
  return 1;
}


int ZTransferFunction::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElementDiscrete::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ZTransferFunction::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElementDiscrete::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int ZTransferFunction::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=RandomSource

void RandomSource::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElementDiscrete::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("IOBlock.input_element_numbers"); 
  edc.TreeDelete("IOBlock.input_element_types"); 
  edc.TreeDelete("IOBlock.input_local_number"); 
}


int RandomSource::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_numbers"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_types"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_local_number"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}


  InputOutputElementDiscrete::SetElementDataAuto(edc);
  return 1;
}


int RandomSource::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElementDiscrete::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RandomSource::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElementDiscrete::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int RandomSource::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=LinearTransformation

void LinearTransformation::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;

  Matrix A_coeff_1(A_coeff);
  ed.SetMatrix(A_coeff_1.GetMatPtr(),A_coeff_1.Getrows(),A_coeff_1.Getcols(),"A_matrix"); ed.SetVariableLength(); ed.SetToolTipText("transformation matrix A: y=A.u+b"); edc.TreeAdd("IOBlock",ed);
  ed.SetVector(b_coeff.GetVecPtr(),b_coeff.Length(),"b_vector"); ed.SetVariableLength(); ed.SetToolTipText("offset vector b: y=A.u+b"); edc.TreeAdd("IOBlock",ed);
}


int LinearTransformation::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataMatrix(GetMBS(), edc, "IOBlock.A_matrix",A_coeff, 1);
  GetElemDataVector(GetMBS(), edc, "IOBlock.b_vector",b_coeff, 1);
  return 1;
}


int LinearTransformation::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LinearTransformation::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LinearTransformation::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=Quantizer

void Quantizer::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int Quantizer::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int Quantizer::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Quantizer::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int Quantizer::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=STransferFunction

void STransferFunction::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetVector(num.GetVecPtr(),num.Length(),"numerator"); ed.SetVariableLength(); ed.SetToolTipText("ascending numerator coefficients n of transfer-function. TF = num/den with num = n(1)*1+n(2)*s+n(3)*s*s+... Will be normalized automatically!"); edc.TreeAdd("IOBlock",ed);
  ed.SetVector(den.GetVecPtr(),den.Length(),"denominator"); ed.SetVariableLength(); ed.SetToolTipText("ascending denominator coeffs d of transfer-function. TF = num/den with den = d(1)*1+d(2)*s+d(3)*s*s+... Will be normalized automatically!"); edc.TreeAdd("IOBlock",ed);
}


int STransferFunction::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataVector(GetMBS(), edc, "IOBlock.numerator",num, 1);
  GetElemDataVector(GetMBS(), edc, "IOBlock.denominator",den, 1);
  return 1;
}


int STransferFunction::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int STransferFunction::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int STransferFunction::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=LinearODE

void LinearODE::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int LinearODE::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int LinearODE::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LinearODE::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int LinearODE::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOMathFunction

void IOMathFunction::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOMathFunction::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOMathFunction::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMathFunction::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMathFunction::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOSaturate

void IOSaturate::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOSaturate::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOSaturate::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOSaturate::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOSaturate::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IODeadZone

void IODeadZone::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IODeadZone::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IODeadZone::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IODeadZone::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IODeadZone::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOProduct

void IOProduct::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOProduct::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOProduct::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOProduct::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOProduct::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOTime

void IOTime::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("IOBlock.input_element_numbers"); 
  edc.TreeDelete("IOBlock.input_element_types"); 
  edc.TreeDelete("IOBlock.input_local_number"); 
}


int IOTime::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_numbers"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_types"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_local_number"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}


  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOTime::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTime::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTime::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOPulseGenerator

void IOPulseGenerator::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  edc.TreeDelete("IOBlock.input_element_numbers"); 
  edc.TreeDelete("IOBlock.input_element_types"); 
  edc.TreeDelete("IOBlock.input_local_number"); 
}


int IOPulseGenerator::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_numbers"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_element_types"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}
  {  Vector dummy(1); ed.SetVector(dummy.GetVecPtr(),1,"input_local_number"); ed.SetLocked(1); edc.TreeAdd("IOBlock",ed);}


  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOPulseGenerator::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOPulseGenerator::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOPulseGenerator::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOTimeWindow

void IOTimeWindow::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOTimeWindow::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOTimeWindow::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTimeWindow::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTimeWindow::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOStopComputation

void IOStopComputation::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOStopComputation::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOStopComputation::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOStopComputation::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOStopComputation::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOElementDataModifier

void IOElementDataModifier::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetBool(modify_at_start_time_step_only,"start_of_timestep_only"); ed.SetToolTipText("modify element data at start time step only."); edc.TreeAdd("IOBlock",ed);
}


int IOElementDataModifier::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataBool(GetMBS(), edc, "IOBlock.start_of_timestep_only",modify_at_start_time_step_only, 1);
  return 1;
}


int IOElementDataModifier::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOElementDataModifier::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOElementDataModifier::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IODisplay

void IODisplay::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(ndigits,"number_of_digits"); ed.SetToolTipText("number of digits"); edc.TreeAdd("IOBlock",ed);
}


int IODisplay::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "IOBlock.number_of_digits",ndigits, 1);
  return 1;
}


int IODisplay::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IODisplay::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IODisplay::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOResponseElement

void IOResponseElement::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOResponseElement::SetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::SetElementDataAuto(edc);
  return 1;
}


int IOResponseElement::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOResponseElement::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOResponseElement::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOKeyResponseElement

void IOKeyResponseElement::GetElementDataAuto(ElementDataContainer& edc)
{
  IOResponseElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOKeyResponseElement::SetElementDataAuto(ElementDataContainer& edc)
{
  IOResponseElement::SetElementDataAuto(edc);
  return 1;
}


int IOKeyResponseElement::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = IOResponseElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOKeyResponseElement::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = IOResponseElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOKeyResponseElement::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOMouseResponseElement

void IOMouseResponseElement::GetElementDataAuto(ElementDataContainer& edc)
{
  IOResponseElement::GetElementDataAuto(edc);
  ElementData ed;
}


int IOMouseResponseElement::SetElementDataAuto(ElementDataContainer& edc)
{
  IOResponseElement::SetElementDataAuto(edc);
  return 1;
}


int IOMouseResponseElement::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = IOResponseElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMouseResponseElement::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = IOResponseElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMouseResponseElement::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOMinMax

void IOMinMax::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(mode,"mode",1,6); ed.SetToolTipText("1..min, 2..max, 3..avg, 4..min(abs), 5..max(abs), 6..avg(abs)"); edc.TreeAdd("IOBlock",ed);
  ed.SetDouble(start_time,"start_time"); ed.SetToolTipText("Up to this point of time, the output is equal to the input. Afterwards the output is computed according to the mode."); edc.TreeAdd("IOBlock",ed);
}


int IOMinMax::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "IOBlock.mode",mode, 1);
  GetElemDataDouble(GetMBS(), edc, "IOBlock.start_time",start_time, 1);
  return 1;
}


int IOMinMax::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMinMax::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOMinMax::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//automatically generated EDC-access
//actual class=IOTCPIPBlock

void IOTCPIPBlock::GetElementDataAuto(ElementDataContainer& edc)
{
  InputOutputElement::GetElementDataAuto(edc);
  ElementData ed;
  ed.SetInt(port_number,"port_number"); ed.SetToolTipText("Port number, e.g. '50000'."); edc.TreeAdd("IOBlock",ed);
  ed.SetText(ip_address.c_str(),"ip_address"); ed.SetToolTipText("IP address, e.g. '127.0.0.1' (localhost). Do not neglect the dots between the numbers."); edc.TreeAdd("IOBlock",ed);
  ed.SetInt(n_output,"received_data_size"); ed.SetToolTipText("Number of received data values (outputs). This number has to be consistent with the transmitted data values of the other communication side (the additional double for the communication flags is not corresponding to this number)."); edc.TreeAdd("IOBlock",ed);
  ed.SetInt(timeout,"timeout"); ed.SetToolTipText("TCP/IP timeout in milliseconds; default is 10000."); edc.TreeAdd("IOBlock",ed);
}


int IOTCPIPBlock::SetElementDataAuto(ElementDataContainer& edc)
{
  ElementData ed;


  InputOutputElement::SetElementDataAuto(edc);
  GetElemDataInt(GetMBS(), edc, "IOBlock.port_number",port_number, 1);
  GetElemDataText(GetMBS(), edc, "IOBlock.ip_address",ip_address, 1);
  GetElemDataInt(GetMBS(), edc, "IOBlock.received_data_size",n_output, 1);
  GetElemDataInt(GetMBS(), edc, "IOBlock.timeout",timeout, 1);
  return 1;
}


int IOTCPIPBlock::ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata) // Read access to a single element variable 
{
  int rv = InputOutputElement::ReadSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTCPIPBlock::WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata) // Write access to a single element variable
{
  int rv = InputOutputElement::WriteSingleElementDataAuto(RWdata);
  if (rv==1) return 1;
  return 0;
}


int IOTCPIPBlock::GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables) // returns all available directly (ReadWrite-) accessable variables
{
  return 0;
}


void MBSObjectFactory::AddObjectInfos_Auto()
{
  AddObjectInfo(OFCElement,"Mass1D", TAEBody+TAE1D, "Mass1D.txt", "A point mass in one dimensions with 1 position coordinate. The computation of the dynamics of the point mass is extremely simple. The Mass1D can be used for a lot of applications which can be represented by the same type of equations. If you interpret the 'mass' to be 'moment of inertia' and the 'position' to be 'angle', then you can realize a 1D rotatory element as well.");
  AddObjectInfo(OFCElement,"Rotor1D", TAEBody+TAE1D, "Rotor1D.txt", "A rotor with 1 degree of freedom (the rotation). Mathematically implemented like Mass1D but different geometric representation.");
  AddObjectInfo(OFCElement,"Mass2D", TAEBody+TAE2D, "mass2D.txt", "A point mass in two dimensions with 2 position coordinates. The computation of the dynamics of the point mass is extremely simple, thus the Mass2D can be used for many body simulations (e.g. particles).");
  AddObjectInfo(OFCElement,"Rigid2D", TAEBody+TAE2D, "Rigid2D.txt", "A rigid body in 2D.");
  AddObjectInfo(OFCElement,"Mass3D", TAEBody, "AddElement.txt", "A point mass in three dimensions with 3 position coordinates. The computation of the dynamics of the point mass is extremely simple, thus the Mass3D can be used for many body simulations (e.g. particles).");
  AddObjectInfo(OFCElement,"NodalDiskMass3D", TAEBody, "NodalDiskMass3D.txt", "This is a disk mass for the purpose of rotordynamics applications and should be used together with the RotorBeamXAxis element.");
  AddObjectInfo(OFCElement,"Rigid3D", TAEBody, "Rigid3D.txt", "A rigid body in 3D.");
  AddObjectInfo(OFCElement,"Rigid3DKardan", TAEBody, "Rigid3DKardan.txt", "A rigid body in 3D, implemented with bryant angles.");
  AddObjectInfo(OFCElement,"LinearBeam3D", TAEBody, "LinearBeam3D.txt", "The Beam3D element is a three dimensional elastic beam element which is aligned along the local x axis. It provides a decoupled calculation of bending in y and z direction, axial deformation in x direction and torsion about the x axis. Shear deformation is not considered. The decoupled calculation is a simplification of the real, nonlinear problem, but for small deformations the results coincidence highly with the exact solution.");
  AddObjectInfo(OFCElement,"RotorBeamXAxis", TAEBody, "RotorBeamXAxis.txt", "The RotorBeamXAxis element is a three dimensional elastic rotor beam element. It has exact the same characteristics and properties as the LinearBeam3D element except two differences. The first difference is that for a rotor element it is necessary to enable big rotation about the rotor axis instead of the small rotation of the LinearBeam3D. The second difference is that all element DOF are stored w.r.t. local beam coordinate system.");
  AddObjectInfo(OFCElement,"ANCFCable2D", TAEBody+TAE2D+TAENotInRelease, "", "The ANCFCable2D is a beam based on the Absolute Nodal Coordinate Formulation (ANCF). For details of this beam, see the literature J. Gerstmayr, H. Irschik. On the correct representation of bending and axial deformation in the absolute nodal coordinate formulation with an elastic line approach, Journal of Sound and Vibration, Vol. 318, pp. 461-487, 2008. DOI:10.1016/j.jsv.2008.04.019");
  AddObjectInfo(OFCElement,"ANCFBeamShear3DLinear", TAEBody, "ANCFBeamShear3DLinear.txt", "ANCFBeamShear3DLinear is an ANCF beam finite element for multibody dynamics systems which is capable of large deformations and can be used for static as well as dynamic investigations. The beam finite element can reproduce axial, bending, shear and torsional deformation. A linear interpolation for the geometry and the displacement along the beam axis is chosen.\\The definition of the beam finite element is based on the absolute nodal coordinate formulation (ANCF), which uses slope vectors for the parameterization of the orientation of the cross section instead of rotational parameters. Two different formulations for the elastic forces of the beam elements are presented:\\(1) A structural mechanics based formulation of the elastic forces based on Reissner's nonlinear rod theory including generalized strain measures. A term accounting for thickness and cross section deformation is included and shear locking is prevented.\\(2) A continuum mechanics based formulation of the elastic forces for a St.Venant Kirchhoff material which avoids the Poisson and shear locking phenomenon.");
  AddObjectInfo(OFCElement,"ANCFBeamShear3DQuadratic", TAEBody, "ANCFBeamShear3DQuadratic.txt", "ANCFBeamShear3DQuadratic is an ANCF beam finite element for multibody dynamics systems which is capable of large deformations and can be used for static as well as dynamic investigations. The beam finite element can reproduce axial, bending, shear and torsional deformation. A quadratic interpolation for the geometry and the displacement along the beam axis is chosen.\\The definition of the beam finite element is based on the absolute nodal coordinate formulation (ANCF), which uses slope vectors for the parameterization of the orientation of the cross section instead of rotational parameters. Two different formulations for the elastic forces of the beam elements are presented:\\(1) A structural mechanics based formulation of the elastic forces based on Reissner's nonlinear rod theory including generalized strain measures. A term accounting for thickness and cross section deformation is included and shear locking is prevented.\\(2) A continuum mechanics based formulation of the elastic forces for a St.Venant Kirchhoff material which avoids the Poisson and shear locking phenomenon.");
  AddObjectInfo(OFCElement,"ANCFBeam3DTorsion", TAEBody, "ANCFBeam3DTorsion.txt", "ANCFBeam3DTorsion is a Bernoulli-Euler beam finite element in ANCF (Absolute Nodal Coordinate Formulation) capable of large axial, bending, and torsional deformations.");
  AddObjectInfo(OFCElement,"ANCFBeamShearFE2DLinear", TAEBody+TAENotInRelease, "", "");
  AddObjectInfo(OFCElement,"ANCFBeamShearFE2DQuadratic", TAEBody+TAENotInRelease, "", "");
  AddObjectInfo(OFCElement,"ANCFThinPlate3D", TAEBody+TAENotInRelease, "", "");
  AddObjectInfo(OFCElement,"PointJoint", TAEconstraint, "PointJointShort.txt", "The PointJoint constrains two elements at a local position or node each. If only one element is specified (second element 0), a ground PointJoint is realized.");
  AddObjectInfo(OFCElement,"CoordinateConstraint", TAEconstraint+TAEspecial_connector, "CoordinateConstraint.txt", "The CoordinateConstraint constrains two elements by constraining a single coordinate of each element, e.g. the x-displacement of two different elements. If the second element number is zero, a groundjoint can be realized. The CoordinateConstraint uses the lagrange multiplier formulation by default, which means that there is no constraint violation at all. For static problems, the lagrange multiplier constraint formulation is applied directly, by adding the kinematical conditions to the nonlinear system equations. In dynamic (time dependent) simulations, the constraint is solved on the position (displacement) levelwith index 3 solvers and on the velocity level with index 2 solvers. Alternatively, the penalty formulation can be used, which means that a certain (very high) spring stiffness is used instead of lagrange multipliers. Thus, no additional equation is added, however, the systemequations may become unsolvable stiff (ill conditioned) in case of static problems; for dynamical problems, the very high stiffness might lead to high-frequency oscillations, inaccurate solutions or no convergence.");
  AddObjectInfo(OFCElement,"VelocityCoordinateConstraint", TAEconstraint+TAEspecial_connector, "VelocityCoordinateConstraint.txt", "Similar to CoordinateConstraint. Lagrangian constraint implemented for index 3 and index 2 solvers. A penalty formulation is also implemented.");
  AddObjectInfo(OFCElement,"RollingJoint3D", TAEconstraint+TAENotInRelease, "", "A point at fixed relative position (in global coordinates) with respect to the center of the body must have zero velocity (or nearly zero). The attached body must be a rigid body!");
  AddObjectInfo(OFCElement,"SlidingPointJoint", TAEconstraint, "SlidingPointJoint.txt", "This joint enables sliding of a fixed point of a body i along the x - axis of another body j. Both body i and body j can be flexible or rigid. Body j can contain more than one elements. No rotations are constrained at all. Only a Lagrangian formulation is implemented, the penalty formulation is not implemented yet. A MaxIndex 2 and 3 formulation exists.");
  AddObjectInfo(OFCElement,"SlidingPrismaticJoint", TAEconstraint, "SlidingPrismaticJoint.txt", "This joint enables sliding of a fixed point of a body i along the x - axis of another body j. Both body i and body j can be flexible or rigid. Body j can contain more than one element. The difference to the SlidingPointJoint is that the relative rotation between the bodies is also constrained. A Lagrangian formulation is used for both stiff and springy constrained rotation. For the position constraint only a stiff formulation exists. A penalty formulation is not implemented yet. There is a MaxIndex 2 and 3 formulation implemented.");
  AddObjectInfo(OFCElement,"Rope3D", TAEconstraint, "Rope3D.txt", "Elastic rope that is always under tension and can be fixed to multiple bodies and ground. There are 2 different kinds of suspensions points. Suspension points fixed on the ground are defined with the element number 0 and the global position. Suspension points on bodies are defined with the element number and the corresponding local position.");
  AddObjectInfo(OFCElement,"FrictionConstraint", TAEconstraint+TAEspecial_connector, "FrictionConstraint.txt", "The FrictionConstraint is acting on an arbitraty coordinate, including rotations. It can be used to connect two elements to each other or one element to ground. Up to a specified threshold of the force, the constraint is sticking, which is realized by a spring-damper formulation. Above this threshold, a constant friction force is applied during the sliding phase. Alternatively sticking can be switched off and a coulomb friction force, with a transition region for very small velocities, can be applied.");
  AddObjectInfo(OFCElement,"Contact1D", TAEconstraint+TAEspecial_connector, "Contact1D.txt", "Contact1D realizes a contact formulation between two elements or one element and ground. Only one coordinate (direction) is considered per element.");
  AddObjectInfo(OFCElement,"GenericBodyJoint", TAEconstraint, "GenericBodyJointShort.txt", "The GenericBodyJoint constrains two elements at a local position each. If only one element is specified (second element 0), a ground GenericBodyJoint is realized. A penalty and lagrange formulation is available.");
  AddObjectInfo(OFCElement,"RevoluteJoint", TAEconstraint, "RevoluteJointShort.txt", "The RevoluteJoint constrains all relative degrees of freedom between two bodies except the rotation about a local rotation axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.");
  AddObjectInfo(OFCElement,"PrismaticJoint", TAEconstraint, "PrismaticJointShort.txt", "The PrismaticJoint constrains all relative degrees of freedom between two bodies except the translation along a local sliding axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.");
  AddObjectInfo(OFCElement,"UniversalJoint", TAEconstraint, "UniversalJoint.txt", "The UniversalJoint constains the local position of two elements and keeps two axes, one on each body, perpendicular to each other.");
  AddObjectInfo(OFCElement,"RigidJoint", TAEconstraint, "RigidJointShort.txt", "The RigidJoint constrains the position and relative angles of an element at a specified local position. If only one element is specified, a ground joint is realized. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.");
  AddObjectInfo(OFCElement,"CylindricalJoint", TAEconstraint, "CylindricalJointShort.txt", "The CylindricalJoint constrains like the RevoluteJoint, but allows additionally translation along the rotational axis. A penalty formulation exists, which replaces the exact lagrange constraint by a approximation with joint stiffness and damping.");
  AddObjectInfo(OFCElement,"SpringDamperActuator", TAEconstraint, "SpringDamperActuator.txt", "The Spring-Damper-Actuator connects two points with a spring, a damper and a actor element, in which actuator force fa remains constant. The resultant force is applied in the connection line of these points. There are different modes available, how the spring and damper force is calculated. It is also possible to change the neutral spring length. This joint is realized in Penalty formulation only.");
  AddObjectInfo(OFCElement,"RigidLink", TAEconstraint, "RigidLink.txt", "A rigid link is a rigid constraint element that provides a stiff connection between nodes or positions in the model. In standard mode the distance between the connected points remains constant. In extended mode it is possible to change the distance as a function of time or input.  There is only a Lagrange formulation implemented.");
  AddObjectInfo(OFCElement,"RotatorySpringDamperActuator", TAEconstraint, "RotationalSpringDamperActuator.txt", "The RotatorySpringDamperActuator connects two elements with rotatory spring, damper and a constant actuator moment ma. Positive rotation around rotation axis according to right hand rule. There are different modes available, how the spring and damper moment is calculated. It is also possible to change the neutral spring angle. This joint is realized in Penalty formulation only.");
  AddObjectInfo(OFCElement,"SpringDamperActuator2D", TAEconstraint, "SpringDamperActuator2D.txt", "The SpringDamperActuator2D is a simplified version of the SpringDamperActuator for 2D elements. Nodes are not supported in the 2D version. Apart from that the constraint has the same functionality as the 3D version. See  the SpringDamperActuator documentation for more information.");
  AddObjectInfo(OFCElement,"PointJoint2D", TAEconstraint, "PointJoint2D.txt", "The PointJoint2D is a simplified version of the PointJoint for 2D elements. It constrains two elements at at local element positions. If only one element is specified (second element 0), a ground PointJoint is realized. It provides both Lagrangian and penalty formulation.");
  AddObjectInfo(OFCElement,"AverageConstraint", TAEconstraint+TAENotInRelease, "", "Bla The AverageConstraint is a Connector which acts on two sets A and B of weighted body points. Each set corresponds two a seperate body of the MBS, set B may also correspond to the global coordinate system (ground). Body points are wheter specified as pairs of element number and local position, or as element number and local node number. If an element number is zero, then a node number is interpreted as global node number and a local position is interpreted as global (ground) position.");
  AddObjectInfo(OFCElement,"IODiscreteTransferFunction", TAEinput_output, "ZTransferFunction.txt", "Discontinuous transfer function in z-space. It is a SISO (single input-single output) control element. Inital state is zero.");
  AddObjectInfo(OFCElement,"IORandomSource", TAEinput_output, "addRandomSource.txt", "Discontinuous random source using alternatively an internal C++ based pseudo random generator or a linear feedback shift register. It has no input and one output.");
  AddObjectInfo(OFCElement,"IOLinearTransformation", TAEinput_output, "LinearTransformation.txt", "Continuous linear transformation. The transfer function type is SISO (single input-single output) or MIMO (multi input-multi output).");
  AddObjectInfo(OFCElement,"IOQuantizer", TAEinput_output, "Quantizer.txt", "A quantizer block passes its input signal through a stair-step function so that many neighboring points on the input axis are mapped to one point on the output axis. The effect is to quantize a smooth signal into a stair-step output. It is a SISO (single input-single output) control element.");
  AddObjectInfo(OFCElement,"IOContinuousTransferFunction", TAEinput_output, "STransferFunction.txt", "The STransferFunction is a linear transfer function for continuous state-space elements. It is a SISO (single input-single output) type.");
  AddObjectInfo(OFCElement,"IOLinearODE", TAEinput_output, "LinearODE.txt", "The LinearODE Element represents a linear ordinary differential equation of SISO (single input-single output) or MIMO (multi input-multi output) type.");
  AddObjectInfo(OFCElement,"IOMathFunction", TAEinput_output, "MathFunction.txt", "A IOMathFunction contains a mathematical expression or a lookup table with different modes for piecewise interpolation. The output is result of the evalutation of the MathFunction as a function of input. It is a SISO (single input-single output) control element.");
  AddObjectInfo(OFCElement,"IOSaturate", TAEinput_output, "Saturate.txt", "Continuous saturation element for upper and lower limits. It is a SISO (single input-single output) control element.");
  AddObjectInfo(OFCElement,"IODeadZone", TAEinput_output, "DeadZone.txt", "Continuous dead-zone element. The outputs between upper and lower limit is zero. This leads to an offset of the input signal by the corresponding lower or upper limit. It is a SISO (single input-single output) control element.");
  AddObjectInfo(OFCElement,"IOProduct", TAEinput_output, "Product.txt", "Continuous product (or division) of one or more inputs. A dedicated exponent for every input and a offset can be applied.");
  AddObjectInfo(OFCElement,"IOTime", TAEinput_output, "addTime.txt", "Continuous time source. This element simply outputs the time.");
  AddObjectInfo(OFCElement,"IOPulseGenerator", TAEinput_output, "addPulseGenerator.txt", "Continuous pulse generator. This element outputs repeating sequence or rectangular pulses after a certain delay. It has no input and one output.");
  AddObjectInfo(OFCElement,"IOTimeWindow", TAEinput_output, "TimeWindow.txt", "This element helps to capture a special time window. It has two inputs and one output.");
  AddObjectInfo(OFCElement,"IOStopComputation", TAEinput_output, "StopComputation.txt", "This element stops the computation, if input is unequal zero. It has one input and no output.");
  AddObjectInfo(OFCElement,"IOElementDataModifier", TAEinput_output, "ElementDataModifier.txt", "This element can be used to modify data of a constraint or element. It has one input and no output.");
  AddObjectInfo(OFCElement,"IODisplay", TAEinput_output, "Display.txt", "This element can be used to display any (single) numberical value fed into the (single) input.");
  AddObjectInfo(OFCElement,"IOResponseElement", TAEinput_output+TAENotInRelease, "", "This element allows the IOBlock to respond to a key or mouse input. ( e.g. pressing '+' or '-' in the IOBlockView Window ) ");
  AddObjectInfo(OFCElement,"IOKeyResponseElement", TAEinput_output+TAENotInRelease, "", "This element allows the IOBlock to respond to a key input. ( e.g. pressing '+' or '-' ) ");
  AddObjectInfo(OFCElement,"IOMouseResponseElement", TAEinput_output+TAENotInRelease, "", "This element allows the IOBlock to respond to a mouse input. ( e.g. clicking the IOElement in the IOBlocksView) ");
  AddObjectInfo(OFCElement,"IOMinMax", TAEinput_output, "IOMinMax.txt", "This block returns the minimum, maximum or average value of the input. Up to a specific point of time, this functionality is switched off and the output y is equal to the input u. This block can be used to postprocess sensor values.");
  AddObjectInfo(OFCElement,"IOTCPIPBlock", TAEinput_output, "TCPIP.txt", "This I/O element is a communication block based on TCP/IP which allows HOTINT to connect to other programs or tools, opening up a large range of possible applications including external control, user-defined ``add-ons'', or even co-simulation. Based on the specified IP (v4) address and port number the IOTCPIPBlock sets up a server socket and waits for a connection request from a client. Hence, HOTINT here plays the server role, and the external program is the client application. Data exchange is performed at a stage before every time step in HOTINT, following below protocol:\\  The outgoing data, i.e.~the data sent from HOTINT to the client, is an array of 8-byte double precision numbers corresponding to the current values of the inputs of the I/O element, and additionally, an 8-byte sequence which is appended to that array and includes communication control flags (see the \textit{Communication flags} - section below for more details). Hence, the total amount of outgoing data is (number of inputs\,+\,1) times 8 bytes (double precision numbers). After the client has received and processed that data, it sends back a data package to HOTINT -- the incoming data for the I/O element --  which again consists of an array of double precision numbers, this time with the length (number of outputs\,+\,1). The first (number of output) double precision values determine the outputs of the I/O element, and the last 8 bytes again are used for the transfer of communication flags.\\ HOTINT now begins the computation of one time step, where the transmitted data from the client is accessible via the outputs of the IOTCPIPBlock. \vspace*{12pt} \\ \textit{Important notes} \vspace*{6pt} \\ -- The waiting procedure for the client connection request, as well as the send and receive operations all are so-called ``blocking calls''. This means that HOTINT will wait for those operations to finish, and during that time, not respond to any user input. Therefore, a reasonable timeout (default is 10 seconds) should be specified for the IOTCPIPBlock to allow TCP/IP connection or transmission error handling. \vspace*{6pt}\\ -- You will probably have to adjust your firewall settings and set appropriate permissions for HOTINT and the client application. \vspace*{6pt}\\ --  Depending on the implementation of the client, it might be neccessary to start the server, i.e., HOTINT, first. \vspace*{6pt} \\ -- Since HOTINT is running on Microsoft Windows, the memory byte order, also called ``endianness'', is ``Little Endian'', which means that the least significant bytes/digits are stored ``first'' in memory, i.e., on the smallest memory address. Therefore, any data sent from or received by the IOTCPIPBlock has or must have that byte order, respectively. You probably have to take that into account on the client side, especially if the client is running on a different platform and/or architecture on another computer. \vspace*{12pt}\\ \textit{Communication flags} \vspace*{6pt} \\ Currently, the following 4-byte flags are implemented: \vspace*{6pt}\\ (1) Neutral flag: \texttt{0x00000000} (integer value: \texttt{0}). This flag signals that the application is running (properly) and no further action is required. \\ (2) Reset flag: \texttt{0x00000001} (integer value: \texttt{1}). This flag is sent from HOTINT to the client in the first step of the computation. This can be used, for instance, to reset the client application.\\ (3) Error flag: \texttt{0x00000002} (integer value: \texttt{2}). Indicates that an error has occurred. If HOTINT receives the error flag, an error message is issued, the connection is closed and the program execution terminated. \\ (4) Close flag: \texttt{0x00000003} (integer value: \texttt{3}). This flag is sent from HOTINT to the client to indicate that the computation has finished and the connection will be closed, which is the case when the computation has actually finished, or the ``Stop''-button has been hit. \\ (5) Any other value: Treated as error flag (3). \vspace*{6pt}\\   One of these flags is stored in and read from the last 8 bytes of the exchanged data -- corresponding to one additional double precision number -- in either direction in every time step. Currently, for simplicity, the flag is just casted explicitly from an integer to a double precision number which then can be transmitted and casted back to an integer exactly. Of course, this procedure must be followed on both the server and the client side.");

}

int MBSObjectFactory::AddElement(int element_type_id)
{
  int rv = -1;
  switch(element_type_id)
  {
    case 1: 
    {
      Mass1D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 2: 
    {
      Rotor1D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 3: 
    {
      Mass2D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 4: 
    {
      Rigid2D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 5: 
    {
      Mass3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 6: 
    {
      NodalDiskMass3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 7: 
    {
      Rigid3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 8: 
    {
      Rigid3DKardan elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 9: 
    {
      Beam3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 10: 
    {
      RotorBeamXAxis elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 11: 
    {
      ANCFCable2D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 12: 
    {
      ANCFBeamShear3DLinear elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 13: 
    {
      ANCFBeamShear3DQuadratic elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 14: 
    {
      ANCFBeam3DTorsion elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 15: 
    {
      ANCFBeamShearFE2DLinear elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 16: 
    {
      ANCFBeamShearFE2DQuadratic elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 17: 
    {
      ANCFSimpleThinPlate3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 18: 
    {
      BasePointJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 19: 
    {
      CoordConstraint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 20: 
    {
      VelCoordConstraint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 21: 
    {
      RollingJoint3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 22: 
    {
      SlidingPointJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 23: 
    {
      SlidingPrismaticJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 24: 
    {
      Rope3D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 25: 
    {
      FrictionConstraint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 26: 
    {
      Contact1D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 27: 
    {
      BaseBodyJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 28: 
    {
      RevoluteJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 29: 
    {
      PrismaticJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 30: 
    {
      UniversalJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 31: 
    {
      RigidJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 32: 
    {
      CylindricalJoint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 33: 
    {
      SpringDamperActuator elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 34: 
    {
      RigidLink elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 35: 
    {
      RotatorySpringDamperActuator elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 36: 
    {
      SpringDamperActuator2D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 37: 
    {
      PointJoint2D elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 38: 
    {
      AverageConstraint elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 39: 
    {
      ZTransferFunction elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 40: 
    {
      RandomSource elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 41: 
    {
      LinearTransformation elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 42: 
    {
      Quantizer elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 43: 
    {
      STransferFunction elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 44: 
    {
      LinearODE elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 45: 
    {
      IOMathFunction elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 46: 
    {
      IOSaturate elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 47: 
    {
      IODeadZone elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 48: 
    {
      IOProduct elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 49: 
    {
      IOTime elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 50: 
    {
      IOPulseGenerator elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 51: 
    {
      IOTimeWindow elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 52: 
    {
      IOStopComputation elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 53: 
    {
      IOElementDataModifier elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 54: 
    {
      IODisplay elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 55: 
    {
      IOResponseElement elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 56: 
    {
      IOKeyResponseElement elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 57: 
    {
      IOMouseResponseElement elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 58: 
    {
      IOMinMax elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }
    case 59: 
    {
      IOTCPIPBlock elem(GetMBS());
      rv = GetMBS()->AddElement(&elem);
      break;
    }

  default: ; 
  }
  return rv;
}
