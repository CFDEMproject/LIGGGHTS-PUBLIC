//#**************************************************************
//#
//# filename:             element.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
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
 
#include "element.h"
#include "constraint.h"
#include "node.h"
#include "elementdataaccess.h"
#include "material.h"
#include "femathhelperfunctions.h"
#include "sensors.h"
#include "rendercontext.h"
#include "options_class_auto.h"
#include "geomelements.h"
#include "solversettings_auto.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Element   Element   Element   Element   Element   Element   
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Element::~Element() 
{
	////$ DR 2012-10:[ loads moved from element to mbs, old code
	//for (int i=1; i <= loads.Length(); i++)
	//{
	//	delete loads(i);
	//}
	////$ DR 2012-10:] loads moved from element to mbs, old code


	/*
	for (int i=1; i <= drawelements.Length(); i++)
	{
	delete drawelements(i); //drawelements are now integer values!
	}*/

};

Element* Element::GetCopy()
{
	Element* ec = new Element();
	ec->CopyFrom(*this);
	return ec;
}

void Element::CopyFrom(const Element& e)
{
	mbs = e.mbs;
	x_init = e.x_init;
	data_init = e.data_init;
	//rho = e.rho; //$ DR 2013-02-04 deleted rho from class element, do not use it here!
	materialnum = e.materialnum;
	col = e.col;
	mass = e.mass;
	damping_m = e.damping_m;
	altshape = e.altshape;
	type = e.type;
	elementname = e.elementname;
	draw_element = e.draw_element;


	//int arrays can be simply copied
	ltg = e.ltg;
	ltgdata = e.ltgdata;
	constraintindices = e.constraintindices;

	//pointer arrays must be treated carefully: delete and copy!!!

	////$ DR 2012-10:[ loads moved from element to mbs, old code
	//for (int i=1; i <= loads.Length(); i++)
	//{
	//	delete loads(i);
	//}
	////$ DR 2012-10:] loads moved from element to mbs, old code

	for (int i=1; i <= constraints.Length(); i++)
	{
		delete constraints(i);
	}
	for (int i=1; i <= constraints_nodouble.Length(); i++)
	{
		delete constraints_nodouble(i);
	}
	loads.Flush();
	constraints.Flush();
	constraints_nodouble.Flush();

	////$ DR 2012-10:[ loads moved from element to mbs, old code
	//for (int i=1; i <= e.loads.Length(); i++)
	//{
	//	MBSLoad* l = e.loads(i)->GetCopy();
	//	loads.Add(l);
	//}
	////$ DR 2012-10:] loads moved from element to mbs, old code
	loads = e.loads;	//$ DR 2012-10: loads moved from element to mbs

	//doesnt really work, because 1 constraint belongs to several elements!!!
	for (int i=1; i <= e.constraints.Length(); i++)
	{
		Constraint* c = (Constraint*)(e.constraints(i)->GetCopy());
		constraints.Add(c);
	}
	for (int i=1; i <= e.constraints_nodouble.Length(); i++)
	{
		Constraint* c = (Constraint*)(e.constraints_nodouble(i)->GetCopy());
		constraints_nodouble.Add(c);
	}

	drawelements = e.drawelements;

	sensors = e.sensors;
	elements = e.elements;

	/*
	drawelements.Flush();
	for (int i=1; i <= e.drawelements.Length(); i++)
	{
	GeomElement* de = e.drawelements(i)->GetCopy();
	drawelements.Add(de);
	}
	*/

}


int Element::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = 0;

	//check node numbers
	if (PerformNodeCheck())
	{
		for (int i = 1; i <= NNodes(); i++)
		{
			if (NodeNum(i) < 1 || (NodeNum(i) > GetMBS()->NNodes() && !GetNode(i).IsAuxNode())) 
			{
				errorstr += mystr("    Element ") + mystr(GetOwnNum()) + mystr(": Node ") + 
					mystr(i) + mystr(" has an invalid node number ") + mystr(NodeNum(i)) + mystr("!\n");
				rv = Maximum(rv,1); 
			}
		}
	}

	//check constraint numbers
	for (int i=1; i <= NE(); i++)
	{
		/*
		if (elements(i) <= 0 || elements(i) > GetMBS()->NE())
		{
		errorstr += mystr("    Element ") + mystr(GetOwnNum()) + mystr(" (possibly a constraint) has invalid element number:") + mystr(elements(i)) + mystr("!\n");
		rv = 1;
		}
		//does not work for contact
		if (!GetMBS()->GetElement(elements(i)).IsType(TBody) && !GetMBS()->GetElement(elements(i)).IsType(TController))
		{
		errorstr += mystr("    Element ") + mystr(GetOwnNum()) + mystr(" (possibly a constraint) refers to an invalid element:") + mystr(elements(i)) + mystr(" which is not a valid body!\n");
		rv = 1;
		}*/
	}

	//$ DR 2012-10: loads moved from element to mbs
	int rv_l=0;
	for (int i=1; i <= NLoads(); i++)
	{
		if(loads(i)<=mbs->NLoads())
		{
			rv_l= mbs->GetLoadPtr(loads(i))->CheckConsistency(*this, errorstr);
			rv = Maximum(rv,rv_l);
		}
		else
		{
			errorstr += mystr("    The load number ") + mystr(loads(i)) +mystr(" does not exist in the mbs!\n");
			rv = Maximum(rv,2); 
		}
	}

	if(GetAltShape())
	{
		for (int i=1; i <= NGeomElements(); i++)
		{
			//$MaSch 12-07-2012: bug-fix in consistency checking of draw elements
			//if(drawelements(i)<=GetMBS()->NDrawElements()) 
			if(drawelements(i)>GetMBS()->NDrawElements()) 
			{
				errorstr += mystr("    Element ") + mystr(GetOwnNum()) + mystr(": Graphics.geom_elements contains the invalid number ") + 
					mystr(drawelements(i)) + mystr("!\n");
				rv = Maximum(rv,2); 				
			}
		}
	}


	//other checks need to be performed in derived class
	return rv;
}

void Element::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	//ed.SetInt(GetOwnNum(), "Element_number"); ed.SetLocked(1); edc.Add(ed);

	//mystr eleminfo = GetElementSpec();
	////if (SOS()) eleminfo = eleminfo + mystr(",  ") + mystr(SOS()) + mystr(" ODE(2nd order)");
	////if (ES()) eleminfo = eleminfo + mystr(", ") + mystr(ES()) + mystr(" ODE");
	////if (IS()) eleminfo = eleminfo + mystr(", ") + mystr(IS()) + mystr(" AE");
 // ed.SetText(eleminfo.c_str(), "element_specification"); ed.SetLocked(1); 
	//edc.TreeAdd("Info",ed);
	
	//ed.SetText(GetElementName(), "Element_name"); edc.Add(ed);

	//constraints:
	if (NC())
	{
		mystr cinfo;
		for (int i=1; i<= NC(); i++)
		{
			cinfo = cinfo + mystr(GetConstraint(i)->GetOwnNum());
			if (i != NC()) cinfo = cinfo + mystr(", ");
		}
		ed.SetText(cinfo.c_str(), "constraint_elements"); ed.SetToolTipText("These constraints are attached to the element"); ed.SetLocked(1); edc.TreeAdd("Links",ed);

	}

	//ed.SetVector3D(col.X(), col.Y(), col.Z(), "RGB_color"); ed.SetToolTipText("[red, green, blue], range = 0..1"); edc.Add(ed);

	//if (!IsType(TConstraint))
	//{
		//ed.SetBool(altshape, "Use_alternative_shape"); ed.SetToolTipText("Only draw Geom-Objects that are attached to the element"); edc.Add(ed);
		//SetElemDataIVector(ed, drawelements, "geom_elements"); ed.SetVariableLength(); ed.SetToolTipText("Set Geometric elements to represent body 'elem1, elem2, ...' or empty"); edc.TreeAdd("Graphics", ed);

		//if (GetMaterialNum() == 0)
		//{
		//	ed.SetDouble(rho, "density"); edc.TreeAdd("Physics",ed);
		//}
		//else
		//{
		//	ed.SetInt(GetMaterialNum(), "material_number"); edc.TreeAdd("Physics",ed);
		//}
	//}
	//ed.SetBool(draw_element, "Draw_element"); ed.SetToolTipText("Flag to draw element"); edc.Add(ed);

	//$ DR 2012-10:[ moved loads from element to mbs
	//action buttons for loads:
	//for (int i = 1; i <= NLoads(); i++)
	//{
	//	ed.SetWCDriverAction(2, i, GetOwnNum(), GetLoad(i).LoadName()); 
	//	ed.SetToolTipText("Open load dialog"); 
	//	ed.SetGroup(1);
	//	edc.Add(ed);
	//}
	//$ DR 2012-10:] moved loads from element to mbs


	GetElementDataAuto(edc);
}

int Element::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = 1;
	int pos;

	mystr type_old = GetElementSpec();

	if (!IsType(TConstraint))
	{
		//pos = edc.Find("Use_alternative_shape");
		//if (pos)
		//{
		//	const ElementData& ed = edc.Get(pos);
		//	if (ed.IsBool())
		//		altshape = ed.GetBool();
		//}
		//GetElemDataIVector(GetMBS(), edc, "Geom_elements", drawelements, 0);		//$ DR 2012-12-05	removed, because already in auto function

		//GetElemDataInt(GetMBS(), edc, "Material_number", GetMaterialNum(), 0);	//$ DR 2012-12-05	removed, because already in auto function


		////old: should be finally erased!
		////$ DR 2013-02-04 deleted rho from class element, do not use it here!
		//pos = edc.Find("Density");
		//if (pos)
		//{
		//	const ElementData& ed = edc.Get(pos);
		//	if (ed.IsDouble())
		//		rho = ed.GetDouble();
		//}

	}
	//GetElemDataBool(GetMBS(), edc, "Draw_element", draw_element);

	SetElementDataAuto(edc);

	mystr type_new = edc.TreeGetString("element_type");
	if(!type_new.Compare(type_old))
	{
		UO().InstantMessageText("ERROR: You MUST NOT change the type of the element!");
		return 0;
	}

	return rv;
}

//$ DR 2012-10
int Element::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// no parent class, flush the array here ( and only here )
	available_variables.Flush();

	// Automatic entries for this class 
	Element::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	//$ DR 2013-02-21: added if-conditions
	if(SS())		{	available_variables.Add(ReadWriteElementDataVariableType("Internal.DOF",SS(),0,0.,mystr("degrees of freedom (or generalized unknowns) of the element. range: 1-") + mystr(SS()) ));}
	if(SOS())		{
		available_variables.Add(ReadWriteElementDataVariableType("Internal.second_order_variable",SOS(),0,0.,mystr("second order variables of the element. range: 1-") + mystr(SOS()) )) ;
		available_variables.Add(ReadWriteElementDataVariableType("Internal.second_order_variable_velocity",SOS(),0,0.,mystr("velocities of second order variables of the element. range: 1-") + mystr(SOS()) )) ;
	}
	if(ES())		{available_variables.Add(ReadWriteElementDataVariableType("Internal.first_order_variable",ES(),0,0.,mystr("first order variables of the element. range: 1-") + mystr(ES()) )) ;}
	if(IS())		{available_variables.Add(ReadWriteElementDataVariableType("Internal.algebraic_variable",IS(),0,0.,mystr("algebraic variables of the element. range: 1-") + mystr(IS()) )) ;}
	if(DataS())	{available_variables.Add(ReadWriteElementDataVariableType("Internal.data_variable",DataS(),0,0.,mystr("data varibales of the element which are no degrees of freedom (e.g. inelastic strain, contact state, friction state, etc.). range: 1-") + mystr(DataS()) )) ;}

	// Manual WRITE entries for this class
	//available_variables.Add(ReadWriteElementDataVariableType("Internal.mass",0,0,0.,mystr("total mass of rigid body and point mass"), TRWElementDataReadWrite)); //test only!!!!

	return 0;
}


//retrieve value of element variable/function named RWdata.variable_name with specified components (vector/matrix); return 1, if variable found and read; 0 if not found, -1 if no readaccess, -2 if range-fault
int Element::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
// call base class routine ( not required in Element )	
	//int rv = ParentClass::ReadSingleElementData(RWdata);
	//if (rv == 1) return 1;

// manual things to read
	// in case of Element: access the entries in XG and Xdata, with the following strings
	//DOF (<SOS), DOF_velocity (<SOS), second_order_variable (<SOS), second_order_variable_velocity (<SOS), first_order_variable (<ES), algebraic_variable (<IS), Lagrange_multiplier (<IS), XData
	//hint: XG={SOS(),SOS(),ES(),IS()} to introduce the offsets
  
	if(strcmp(RWdata.variable_name.c_str(), "Internal.DOF") == 0)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= SS()) //range check
		{
			RWdata.value = XG(RWdata.comp1); return 1; 
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Internal.second_order_variable") /*RWdata.variable_name == mystr("Internal.second_order_variable")*/)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= SOS()) //range check
		{
			RWdata.value = XG(RWdata.comp1); return 1; 
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Internal.second_order_variable_velocity") /*RWdata.variable_name == mystr("Internal.second_order_variable_velocity")*/)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= SOS()) //range check
		{
			RWdata.value = XG(RWdata.comp1+SOS()); return 1;
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Internal.first_order_variable") /*RWdata.variable_name == mystr("Internal.first_order_variable")*/)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= ES()) //range check
		{
			RWdata.value = XG(RWdata.comp1+2*SOS()); return 1;
		}
		else return -2; 
	}

	if (RWdata.variable_name.CStrCompare("Internal.algebraic_variable") /*RWdata.variable_name ==  mystr("Internal.algebraic_variable")*/) /*|| (RWdata.variable_name == "Internal.Lagrange_multiplier")*/ 
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= IS()) //range check
		{
			RWdata.value = XG(RWdata.comp1+2*SOS()+ES()); return 1;
		}
		else return -2; 
	}

	if (RWdata.variable_name.CStrCompare("Internal.XData") /*RWdata.variable_name ==  mystr("Internal.XData")*/)
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= this->DataS() ) //range check -  ATTENTION DadaS() is a ardcoded value and not necessarily the actual length of the Vector ltgdata
		{
			RWdata.value = XData(RWdata.comp1); return 1;
		}
		else return -2; 

	}
	return ReadSingleElementDataAuto(RWdata);
}


//write value of element variable/function named RWdata.variable_name with specified components (vector/matrix); return 1, if variable found and written; 0 if not found, -1 if no writeaccess, -2 if range-fault
int Element::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 
{
	// call base class routine ( not required in Element )	
	//int rv = ParentClass::WriteSingleElementData(RWdata);
	//if (rv == 1) return 1;
	//manual things to write

	//if(RWdata.variable_name.CStrCompare("Internal.mass") /*RWdata.variable_name == mystr("Internal.mass")*/)
	//{
	//	mass = RWdata.value;
	//}


	return WriteSingleElementDataAuto(RWdata);
}

//$ DR 2013-02
// add available fieldvariables depending on the available kinematicsAccessFunctions
void Element::GetAvailableFieldVariables(TArray<FieldVariableDescriptor> & variables)
{
	// get the max_component
	FieldVariableDescriptor::FieldVariableComponentIndex max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_none;
	if (Dim() == 3) max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_z;
	else if(Dim()==2) max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_y;
	else if(Dim()==1) max_component = FieldVariableDescriptor::FieldVariableComponentIndex::FVCI_x;

	// position, displacement and velocity
	if (GetKinematicsAccessFunctions()&TKAF_position) 
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_position, max_component);
	}
	if (GetKinematicsAccessFunctions()&TKAF_displacement) 
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_displacement, max_component);
	}
	if ((GetKinematicsAccessFunctions()&TKAF_velocity)&&(!GetMBS()->GetSolSet().dostaticcomputation)) 
	{
		FieldVariableDescriptor::AddTypeIntoArray(variables, FieldVariableDescriptor::FVT_velocity, max_component);
	}
}

double Element::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector3D & local_position, bool flagD)
{
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_position: return fvd.GetComponent(flagD ? GetPosD(local_position) : GetPos(local_position));
	case FieldVariableDescriptor::FVT_displacement: return fvd.GetComponent(flagD ? GetDisplacementD(local_position) : GetDisplacement(local_position));
	case FieldVariableDescriptor::FVT_velocity: return fvd.GetComponent(flagD ? GetVelD(local_position) : GetVel(local_position));
	default: assert(0); return 0;
	}
}
double Element::GetFieldVariableValue(const FieldVariableDescriptor & fvd, const Vector2D & local_position, bool flagD)
{
	switch(fvd.VariableType())
	{
	case FieldVariableDescriptor::FVT_position: return fvd.GetComponent(flagD ? GetPos2DD(local_position) : GetPos2D(local_position));
	case FieldVariableDescriptor::FVT_displacement: return fvd.GetComponent(flagD ? GetDisplacement2DD(local_position) : GetDisplacement2D(local_position));
	case FieldVariableDescriptor::FVT_velocity: return fvd.GetComponent(flagD ? GetVel2DD(local_position) : GetVel2D(local_position));
	default: assert(0); return 0;
	}
}

double Element::GetConstraintDrift(double t) const //$ JG 2011-02
{
	if (IS() == 0) return 0; //no constraints

	int oldindex = mbs->MaxIndex();
	mbs->SetMaxIndex(3);
	
	Vector f(IS()); //slow version!!!
	int elnum = GetOwnNum(); //element number
	mbs->GetElement(elnum).EvalG(f, t);
	mbs->SetMaxIndex(oldindex);

	return f.MaxNorm();
}


void Element::SetMaterialNum(int mnum) 
{
	materialnum = mnum; 
}

 
void Element::AddMaterial(const Material& m) 
{
	materialnum = GetMBS()->AddMaterial(m); 
}

const double& Element::Rho() const {return GetMaterial().Density();}
const double& Element::Em() const {return GetMaterial().YoungsModulus();}
const double& Element::Nu() const {return GetMaterial().PoissonRatio();}

double& Element::Rho() {return GetMaterial().Density();}
double& Element::Em()  {return GetMaterial().YoungsModulus();}
double& Element::Nu()  {return GetMaterial().PoissonRatio();}


void Element::AddElement(int en) //add depending element
{
	elements.Add(en);
}

//$ DR 2012-10:[ loads moved from element to mbs
// old code:
//void Element::AddLoad(const MBSLoad& li)
//{
//	MBSLoad* l = new MBSLoad(li);
//	loads.Add(l);
//
//	int elem = li.HasElementDependence();
//	if (elem)
//	{
//		TArray<int> elnums;
//
//		elnums.Add(elem);
//		GetMBS()->GetElement(elem).GetDirectFeedThroughElements(elnums);
//		//UO() << "Elem" << GetOwnNum() << ": Add load dependence elements: " << elnums << "\n";
//
//		for (int i=1; i<=elnums.Length(); i++)
//		{
//			elements.Add(elnums(i));
//		}
//	}
//}

// new code:

void Element::AddLoad(const MBSLoad& li)
{
	// the following code is now in mbs.h/cpp:
	//MBSLoad* l = new MBSLoad(li);
	//loads.Add(l);
	int loadNr = GetMBS()->AddLoad(li);	// add load in mbs
	loads.Add(loadNr);									// add entry to load list of element

	int elem = li.HasElementDependence();
	if (elem)
	{
		TArray<int> elnums;

		elnums.Add(elem);
		GetMBS()->GetElement(elem).GetDirectFeedThroughElements(elnums);

		for (int i=1; i<=elnums.Length(); i++)
		{
			elements.Add(elnums(i));
		}
	}
}

// deletes the i-th load of the element
void Element::DeleteLoad(int i) 
{
	// mbs->DeleteLoad(loads(i));	// do not delete load from mbs, it may be used by another element!

	for (int j=i; j < NLoads(); j++)	// adjust load list of element
	{
		loads(j) = loads(j+1);
	}
	loads.SetLen(loads.Length()-1);
}

//$ DR 2012-10:] loads moved from element to mbs

void Element::JacobianF2(double t, Matrix& m, IVector& colref)
{
	//Default element Jacobian with divided differences

	//store x-Vector!!!
	//set Matrix size m: element->sos x element->ss
	//do not call this function, if sos == 0!!!
  
	int sos = SOS();

	static Vector locjacf20; //for local jacobians
	static Vector locjacf21; //for local jacobians 
	locjacf20.SetLen(sos);
	locjacf21.SetLen(sos);


	double numdiffepsi = GetMBS()->NumSolver().NumDiffepsi();
	double eps;

	//only non-symmetric mode!!!
	locjacf20.SetAll(0);
	EvalF2(locjacf20,t);


	colref.SetLen(0);

	int usesparseK = UseSparseK();
	if (!usesparseK)
	{
		const	TArray<int>& ltg = GetLTGArray();
		for (int i = 1; i <= SS(); i++)
		{
// (AD) changed () to .Get()
			colref.Add(ltg.Get(i));
//			colref.Add(ltg(i));
		}
	}

	//use constraints_nodouble:
	//special for contact where one element is linked several times to the same contact constraint 
	for (int i=1; i <= constraints_nodouble.Length(); i++)
	{
		const Constraint& c = *(constraints_nodouble(i));
		for (int j=1; j <= c.SS(); j++)
		{
			colref.Add(c.LTG(j));
		}
	}

	m.SetSize(sos,colref.Length());

	int begcnt = 1;
	if (!usesparseK)
	{
		if (FastStiffnessMatrix() == 1)
		{
			//fill in d EvalF2/dq for position DOF, remaining matrix filled afterwards
			begcnt = sos+1;
			StiffnessMatrix(m);
		} 
		else if (FastStiffnessMatrix() == 2)
		{
			//fill in d EvalF2/dq for position and velocity DOF, remaining matrix filled afterwards
			begcnt = 2*sos+1;
			StiffnessMatrix(m);
		} 
		else if (FastStiffnessMatrix() == 3)
		{
			//fill in part of d EvalF2/dq for position and velocity DOF, remaining matrix Added afterwards
			m.FillWithZeros();
			StiffnessMatrix(m);
		}
	} 

	//Matrix mold = m;

	double storex;
	for (int i = begcnt; i <= colref.Length(); i++)
	{
		eps = numdiffepsi*Maximum(1e-2,fabs(XGG(colref(i))));

		if (GetMBS()->NumSolver().SymmetricJacobian())
		{
			storex = XGG(colref(i));
			XGG(colref(i)) += eps;

			locjacf21.SetAll(0);
			EvalF2(locjacf21,t);

			XGG(colref(i)) -= 2*eps;
			locjacf20.SetAll(0);
			EvalF2(locjacf20,t);
			XGG(colref(i)) = storex;
		}
		else
		{
			storex = XGG(colref(i));
			XGG(colref(i)) += 2*eps;
			locjacf21.SetAll(0);
			EvalF2(locjacf21,t);
			XGG(colref(i)) = storex;
		}

		//UO() << "f1=" << locjacf21 << "\n";

		if (FastStiffnessMatrix() == 3)
		{
			for (int j=1; j<=sos;j++)
			{
				m(j,i) += (locjacf21(j)-locjacf20(j))/(2.*eps);
			}
		}
		else
		{
			for (int j=1; j<=sos;j++)
			{
				m(j,i) = (locjacf21(j)-locjacf20(j))/(2.*eps);
			}
		}
	}

	/*
	if (sos == 24 && GetMBS()->GetTime() > 1.99)
	{
	Matrix K(24,24);
	StiffnessMatrix(K);
	Matrix Kdiff(24,24);
	Kdiff.CopyFrom(m,1,1,24,24);

	UO() << "diff=" << K-Kdiff << "\n";
	//UO() << "K=    " << K << "\n";
	UO() << "diff.Norm()=" << (Kdiff-K).Norm2() << "\n";
	UO() << "Kdiff.Norm()=" << Kdiff.Norm2() << "\n";
	}*/

}




void Element::JacobianG(double t, Matrix& m, IVector& colref)
{
	//Default element Jacobian with divided differences

	//store x-Vector!!!
	//set Matrix size m: element->is x element->ss

	//do not call this function, if is == 0!!!
	int is = IS();

	static Vector locjacf20; //for local jacobians
	static Vector locjacf21; //for local jacobians 
	locjacf20.SetLen(is);
	locjacf21.SetLen(is);


	double numdiffepsi = GetMBS()->NumSolver().NumDiffepsi();
	double eps;

	//only non-symmetric mode!!!
	locjacf20.SetAll(0);
	EvalG(locjacf20,t);
	double storex;

	const	TArray<int>& ltg = GetLTGArray();
	colref.SetLen(0);

	//constraint equations depend on element coordinates!!!
	if (IsType(TConstraint))
	{
		Constraint* c = (Constraint*)this;
		/*
		for (int i=1; i <= c->NE(); i++)
		{
		const Element& e = c->GetElem(i);
		*/

		for (int i=1; i <= c->NE_nodouble(); i++)
		{
			const Element& e = c->GetElem_nodouble(i);

			for (int j=1; j <= e.SS(); j++)
			{
				colref.Add(e.LTG(j));
				//if (c->ElementsShareDOFs())
				//{
				//	
				//	colref.AddIfNotExists(e.LTG(j));
				//}
				//else
				//{
				//	colref.Add(e.LTG(j));    
				//}
			}
		}
		if (c->ElementsShareDOFs())
		{
			RemoveRedundantEntries(colref, 1);  //$ PG 2013-6-13: um multi-point-constraints in mode=2 zu ermoeglichen:
		}
	}

	//practically add lagrange multipliers:
	for (int i = 1; i <= SS(); i++)
	{
// (AD) changed () to .Get()
		colref.Add(ltg.Get(i));
//		colref.Add(ltg(i));
	}

	//practically unrealistic that a constraint has constraints:
	for (int i=1; i <= constraints.Length(); i++)
	{
		const Constraint& c = *(constraints(i));
		for (int j=1; j <= c.SS(); j++)
		{
			colref.Add(c.LTG(j));
		}
	}

	//UO() << "colref2=" << colref << "\n";

	//UO() << "dependencies=" << colref << "\n";
	m.SetSize(is,colref.Length());

	for (int i = 1; i <= colref.Length(); i++)
	{
		eps = numdiffepsi*Maximum(1e-2,fabs(XGG(colref(i))));

		storex = XGG(colref(i));
		XGG(colref(i)) += 2*eps;
		locjacf21.SetAll(0);
		EvalG(locjacf21,t);
		XGG(colref(i)) = storex;
		//UO() << "f1=" << locjacf21 << "\n";

		for (int j=1; j<=is;j++)
		{
			m(j,i)= (locjacf21(j)-locjacf20(j))/(2.*eps);
		}
	}
	//UO() << "elem-jac=" << m << "\n";

}

void Element::EvalF2(Vector& f, double t) 
{
	// add -C_q^T \lambda
	//UO() << "load-f1=" << f << "\n";

	TMStartTimer(16);//12%
	for (int i=1; i <= loads.Length(); i++)
	{
		//loads(i)->AddElementLoad(f,t);	//$ DR 2012-10:] loads moved from element to mbs, old code
		mbs->GetLoad(loads(i)).AddElementLoad(f,t,this);
	}
	//UO() << "load-f2=" << f << "\n";

	for (int i=1; i <= constraints.Length(); i++)
	{
		constraints(i)->AddElementCqTLambda(t, constraintindices(i),f);
	}
	//UO() << "load-f3=" << f << "\n";
	TMStopTimer(16);
}; 

int slowminvf2warned = 0;
void Element::EvalMinvF2(Vector& f, double t) 
{
	EvalF2(f, t);

	//slow version, should be avoided by optimized elements!
	if (!slowminvf2warned) 
	{
		UO() << "Warning: Very slow version of MinvF2 called, method might be inefficient!\n";
		slowminvf2warned = 1;
	}

	Matrix m(SOS(),SOS());
	EvalM(m, t);
	m.Invert();
	f = m*f;
}; 


double Element::GetPotentialEnergy()
{
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//energy from loads:
	//works at least for constant non-follower loads
	//penalty constraints need to be added

	static Vector f;
	static Vector xg;
	f.SetLen(SOS());
	f.SetAll(0);

	xg.SetLen(SOS());
	for (int i=1; i <= SOS(); i++)
	{
		xg(i) = XG(i);
	}

	double t = GetMBS()->GetTime();

	for (int i=1; i <= loads.Length(); i++)
	{
		//loads(i)->AddElementLoad(f,t);		//$ DR 2012-10: loads moved from element to mbs, old code
		GetMBS()->GetLoad(loads(i)).AddElementLoad(f,t,this);		//$ DR 2012-10: loads moved from element to mbs
	}

	//constraints do not add work, except if penalty approach ...

	//forces add negative energy, because external work ...
	double Epot = -(xg*f);

	return Epot;


}


void Element::BuildDependencies() 
{
	int mss = mbs->GetSystemSize();
	dependencies.SetLen(mss);
	for (int i=1; i <= mss; i++)
	{
		dependencies(i) = 0;
	}
	for (int i=1; i <= SS(); i++)
	{
		dependencies(LTG(i))=1;
	}
	for (int i=1; i <= constraints.Length(); i++)
	{
		Constraint* c = constraints(i);
		for (int j=1; j <= c->SS(); j++)
		{
			dependencies(c->LTG(j)) = 1;
		}
	}
	for (int i=1; i <= elements.Length(); i++)
	{
		if (elements(i))
		{
			const Element& e = GetElem(i);

			for (int j=1; j <=e.SS(); j++)
			{
				dependencies(e.LTG(j))=1;
			}
		}
	}
	for (int i=1; i <= sensors.Length(); i++)
	{
		if (sensors(i) > 0 && sensors(i) <= NSensors())
		{
			//$ YV 2012-06: changed here from MBSSensor to Sensor
			Sensor & s = GetMBS()->GetSensor(i);
			for (int j=1; j<=s.GetNumberOfRelatedElements(); j++)
			{
				UO() << "Elem" << GetOwnNum() << ": Add elem dependence " << s.GetRelatedElementNumber(j) << " for sensor " << sensors(i) << "\n";
				const Element& e = mbs->GetElement(s.GetRelatedElementNumber(j));

				for (int j = 1; j <= e.SS(); j++)
				{
					dependencies(e.LTG(j)) = 1;
				}
			}
		}

	}

}


void Element::DrawElement() 
{
};

void Element::DrawElement2D() 
{
};


void Element::DrawElementLocalFrame() 
{
	double s = GetMBS()->GetDOption(104);

	GetMBS()->ChooseColor(0.3f,0.3f,0.3f);

	Vector3D v1(-0.2*s, 0,0);
	Vector3D v2( s, 0,0);
	Vector3D v3( 0,-0.2*s,0);
	Vector3D v4( 0, s,0);
	Vector3D v5( 0,0,-0.2*s);
	Vector3D v6( 0,0, s);
	v1 = GetPosD(v1);
	v2 = GetPosD(v2);
	v3 = GetPosD(v3);
	v4 = GetPosD(v4);
	v5 = GetPosD(v5);
	v6 = GetPosD(v6);
	double d = GetMBS()->GetDOption(114);
	GetMBS()->MyDrawLine(v1,v2,d);
	GetMBS()->MyDrawLine(v3,v4,d);
	GetMBS()->MyDrawLine(v5,v6,d);

	char str[20];
	sprintf(str, "X%d", GetOwnNum());
	GetMBS()->GetRC()->PrintText3D((float)v2.X(), (float)v2.Y(), (float)v2.Z(), str);
	sprintf(str, "Y%d", GetOwnNum());
	GetMBS()->GetRC()->PrintText3D((float)v4.X(), (float)v4.Y(), (float)v4.Z(), str);
	sprintf(str, "Z%d", GetOwnNum());
	GetMBS()->GetRC()->PrintText3D((float)v6.X(), (float)v6.Y(), (float)v6.Z(), str);
}


void Element::DrawElementVelocityVector() 
{
	double epsilon = 1e-12;
	int modus = GetMBS()->GetOptions()->PostProcOptions()->BodiesVelocityVectorScalingMode(); //1=constant (a), 2=linear (ax), 3=exponential (a(1-exp(-x/b)))
	double a = GetMBS()->GetOptions()->PostProcOptions()->BodiesVelocityVectorScalingA(); //scaling factor
	double b = GetMBS()->GetOptions()->PostProcOptions()->BodiesVelocityVectorScalingB(); //knee factor (for mode 3)
	double thickness = GetMBS()->GetOptions()->PostProcOptions()->BodiesVelocityVectorScalingThickness()*0.0025; //thickness for vectors
	double r = 1./5.; //= min(headsize/length); for scaling headsize appropriately
	if (fabs(a) < epsilon || (modus == 3 && fabs(b) < epsilon)) 
	{
		return;
	}

	GetMBS()->ChooseColor(0.3f,0.3f,0.3f);

	Vector3D p = GetPosD(Vector3D(0.));
	Vector3D vel = GetVelD(Vector3D(0.));
	double absvel = vel.Norm();
	if (absvel < epsilon)
	{
		return;
	}

	double scaling_factor = a;
	if (modus == 1)
	{
		scaling_factor *= 1./absvel;
	}
	else if (modus == 3)
	{
		scaling_factor *= (1.-exp(-absvel/b))/absvel;
	}
	vel *= scaling_factor;
	
	double linethickness = thickness/3.;
	double headsize = thickness;

	// rescaling of headsize and linethickness, if length of vector is smaller than headsize/VelocityVectorMinThicknessRatio
	double length_over_headsize = absvel*scaling_factor*r/headsize;
	if (length_over_headsize < 1)
	{
		linethickness *= length_over_headsize;
		headsize *= length_over_headsize;
	}

	//adjust color of velocity vector - only if field-variable velocity is plotted
	FieldVariableDescriptor* pFvd = GetMBS()->GetActualPostProcessingFieldVariable();
	TArray<FieldVariableDescriptor> elementVariables;
	GetAvailableFieldVariables(elementVariables);
	if (pFvd && pFvd->VariableType() == FieldVariableDescriptor::FVT_velocity && FieldVariableDescriptor::FindTypeInArray(elementVariables, pFvd->VariableType(), pFvd->ComponentIndex1()))
	{
		GetMBS()->DrawColorArrow(GetFieldVariableValue(*pFvd, p, 1), p, p+vel, linethickness, headsize);
	}
	else
	{
		Vector3D col(0.3,0.3,0.3);
		GetMBS()->MyDrawArrow(p, p+vel, col, linethickness, headsize);
	}	
}

Box3D Element::GetBoundingBox() const
{
	if (GetAltShape())
	{
		Box3D b;
		int NgeomElementsInMBS = GetMBS()->NDrawElements();
		for (int i=1; i <= drawelements.Length(); i++)
		{
			if(drawelements(i) <= NgeomElementsInMBS)	//$ DR 2012-12-05 added "if" for bugfix
			{
				GetMBS()->GetDrawElement(drawelements(i))->SetElnum(GetOwnNum());
				GetMBS()->GetDrawElement(drawelements(i))->AddBoundingBox(b);
			}
		}
		return b;
	}
	else
	{
		return GetElementBox();
	}
}

Box3D Element::GetBoundingBoxD() const
{
	if (GetAltShape())
	{
		Box3D b;
		int NgeomElementsInMBS = GetMBS()->NDrawElements();
		for (int i=1; i <= drawelements.Length(); i++)
		{
			if(drawelements(i) <= NgeomElementsInMBS)	//$ DR 2012-12-05 added "if" for bugfix
			{
				GetMBS()->GetDrawElement(drawelements(i))->SetElnum(GetOwnNum());
				GetMBS()->GetDrawElement(drawelements(i))->AddBoundingBoxD(b);
			}
		}
		return b;
	}
	else
	{
		return GetElementBoxD();
	}
}

Box2D Element::GetBoundingBox2D() const
{
	if (GetAltShape())
	{
		Box2D b;
		int NgeomElementsInMBS = GetMBS()->NDrawElements();
		for (int i=1; i <= drawelements.Length(); i++)
		{
			if(drawelements(i) <= NgeomElementsInMBS)	//$ DR 2012-12-05 added "if" for bugfix
			{
				GetMBS()->GetDrawElement(drawelements(i))->SetElnum(GetOwnNum());
				GetMBS()->GetDrawElement(drawelements(i))->AddBoundingBox2D(b);
			}
		}
		return b;
	}
	else
	{
		return GetElementBox2D();
	}
}

Box2D Element::GetBoundingBox2DD() const
{
	if (GetAltShape())
	{
		Box2D b;
		int NgeomElementsInMBS = GetMBS()->NDrawElements();
		for (int i=1; i <= drawelements.Length(); i++)
		{
			if(drawelements(i) <= NgeomElementsInMBS)	//$ DR 2012-12-05 added "if" for bugfix
			{
				GetMBS()->GetDrawElement(drawelements(i))->SetElnum(GetOwnNum());
				GetMBS()->GetDrawElement(drawelements(i))->AddBoundingBox2D(b);
			}
		}
		return b;
	}
	else
	{
		return GetElementBox2DD();
	}
}

void Element::DrawElementAdd() 
{
	int use_cutting_planes = mbs->UseCuttingPlanes();

	if (GetMBS()->GetIOption(134))
	{
		int NgeomElementsInMBS = GetMBS()->NDrawElements();
		for (int i=1; i <= drawelements.Length(); i++)
		{
			if(drawelements(i) <= NgeomElementsInMBS)	//$ DR 2012-12-05 added "if" for bugfix
			{
				GetMBS()->GetDrawElement(drawelements(i))->SetElnum(GetOwnNum());

				if (!use_cutting_planes || !GetMBS()->GetOptions()->ViewingOptions()->CuttingPlaneCutBodiesAltshapes() || mbs->CuttingPlanesAllow(GetMBS()->GetDrawElement(drawelements(i))->GetRefPosD()))
				{
					GetMBS()->GetDrawElement(drawelements(i))->DrawYourself();
				}
			}
		}
	}
}

