//#**************************************************************
//#
//# filename:             AverageConstraint.cpp
//#
//# author:               Gerstmayr Johannes, Reischl Daniel, Peter Gruber
//#
//# generated:						17.April 2013
//# description:          Avearaging Constraint, which computes a weighted average of positions or higher moment_settings
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
 





// TODO:

//
// Velocity, max_index = 2, ComputeMomentumP, etc.  --> DR + PG
// Damping (?) --> DR


//Fragen JG:

//virtual void GetdRotdqT(const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdRotdqT\n";};
////return the derivative d(Rot*vloc)/dq
//virtual void GetdRotvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d) {UO() << "ERROR:called Body3D::GetdRotvdqT\n";};

#include "element.h"
#include "constraint.h"
#include "kinematicpairs.h"
#include "specialconstraints.h"
#include "geomelements.h"
#include "elementdataaccess.h"
#include "graphicsconstants.h"
#include "rendercontext.h"
#include "material.h"
#include "finiteelement3d.h"
#include "averageconstraint.h"
#include "femathhelperfunctions.h"

int AverageConstraint::NE_nodouble() const
{
	return elements_nodouble.Length();
}

void AverageConstraint::SetElement_nodouble()
{
	//$ PG 2013-4-24: fill elements_nodouble with all element numbers except 0 (ground)
	elements_nodouble.Flush();
	sumWeights1 = 0;
	sumWeights2 = 0;

	for(int i=1; i<=elements.Length(); i++)
	{
		if(elements(i))
		{
			elements_nodouble.AddIfNotExists(elements(i));
		}
		sumWeights1 += weights1(i);
	}
	for(int i=1; i<=elements2.Length(); i++)
	{
		if(elements2(i))
		{
			elements_nodouble.AddIfNotExists(elements2(i));
		}
		sumWeights2 += weights2(i);
	}
}

int AverageConstraint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	if (IS() > MAX_IS)
	{
		errorstr = mystr("ERROR: AverageConstraint: Maximum allowed implicit size is ")+mystr(MAX_IS)+mystr("!\n");
		rv = 1;
	}
	
	if(sumWeights1 <= 0)
	{
		errorstr = mystr("ERROR: AverageConstraint: The sum of the weights of kinematic pair 1 has to be greater than zero!\n");
		rv = 1;
	}
	if(sumWeights2 <= 0)
	{
		errorstr = mystr("ERROR: AverageConstraint: The sum of the weights of kinematic pair 2 has to be greater than zero!\n");
		rv = 1;
	}
	for(int i=1; i<=moment_settings.Length(); i++)
	{
		if(GetMomentOrder(i)<0)
		{
			errorstr = mystr("ERROR: AverageConstraint: The order p of the polynomial ") + mystr(i) + mystr(" has to be an integer with p >= 0 !\n");
			rv = 1;
		}
	}

	// can not draw errors ==> rv = 2
	if(!(elements.Length() == weights1.Length() && ((elements.Length() == loccoords.Length())||(elements.Length() == nodes.Length()))))
	{					
		errorstr = mystr("ERROR: AverageConstraint: The number of elements of kinematic pair 1 has to be equal to the number of weights and coordinates (or nodes)!\n");
		rv = 2;
	}
	if(!(elements2.Length() == weights2.Length() && ((elements2.Length() == loccoords2.Length())||(elements2.Length() == nodes2.Length()))))
	{					
		errorstr = mystr("ERROR: AverageConstraint: The number of elements of kinematic pair 2 has to be equal to the number of weights and coordinates (or nodes)!\n");
		rv = 2;
	}
	if(GetAverageConstraintType() == ACT_None)
	{
		errorstr = mystr("ERROR: AverageConstraint: You must not mix formulations with local nodes and local coordinates.\n");
		rv = 2;
	}

	return rv;
}

// for constraints only, these are the necessary (!) access functions!	
void AverageConstraint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));
	if (is_nodal_constraint)
	{
		kaf = (int)(TKinematicsAccessFunctions(kaf+TKAF_node_position+TKAF_D_node_pos_D_q+TKAF_node_ref_conf_position));	
	}
	else
	{
		kaf = (int)(TKinematicsAccessFunctions(kaf+TKAF_position+TKAF_D_pos_D_q+TKAF_ref_conf_position));
		if (UseLocalCoordinateSystem())
		{
			kaf = (int)(TKinematicsAccessFunctions(kaf+TKAF_D_rot_v_D_q));
		}
	}
	
	KinAccFunc.Add(kaf);
}

void AverageConstraint::ElementDefaultConstructorInitialization()
{
	SetPenaltyFormulation(0);
	SetJointLocalFrame(Matrix3D(1.));
	stiffness_in_joint_local_frame = 0;
	SetUseLocalCoordinateSystem(0);
	is_nodal_constraint = 0;
	elements_share_dofs = 0;

	velocity_constraint = 0;
	elementname = GetElementSpec();

	Vector datainit(DataS());
	datainit.SetAll(0.);
	SetDataInit(datainit);

	// dummy values
	// nodes, elements and loccoords are always of length 2 because this is the concept of BasePointJoint
	nodes.Set2(0,0);
	elements.Set2(1,2);

	scaling_factor = 1.;
}

void AverageConstraint::Initialize() 
{
	//do not use BasePointJoint::Initialize because then the wrong SetPos2ToGlobalCoord is called
	Constraint::Initialize();

	if(GetAverageConstraintType() == ACT_Nodes)   // if there exists a '0' in IVector nodes, then this constraint is not a is_nodal_constraint
	{
		is_nodal_constraint = 1;
	}

	SetMomentumPolynomial();

	// find if there are dofs which are shared by more than one element
	if (!ElementsShareDOFs())
	{
		TArray<int> dofs;
		TArray<int> flags_redundant;
		for (int i=1; i<=NE_nodouble(); i++)
		{
			//Element* el_ptr = &GetElem_nodouble(i);
			dofs.Append(GetElem_nodouble(i).GetLTGArray());
			/*for (int j=1; j<=el_ptr->SS(); j++)
			{
				dofs.Add(el_ptr->LTG(j));
			}*/
		}
		FindRedundantEntries_QSort(dofs, flags_redundant, 0);
		if (flags_redundant.Find(1))
		{
			SetElementsShareDOFS();
		}
	}
}

// returns the type of constraint (coord based or node based) depending on the data in the TArrays
AverageConstraintType AverageConstraint::GetAverageConstraintType() const
{
	TArray<int> all_nodes(nodes); all_nodes.Append(nodes2);

	if (all_nodes.Length() > 0)
	{
		bool has_positive_entries = false;
		bool all_entries_are_positive = true;

		for (int i=1; i<=all_nodes.Length(); i++)
		{
			if (all_nodes(i)>0)
			{
				has_positive_entries = true;
			}
			if (all_nodes(i)<=0)
			{
				all_entries_are_positive = false;
			}
		}

		if (has_positive_entries)
		{
			if (all_entries_are_positive)
			{
				return ACT_Nodes;
			}
			else
			{
				return ACT_None;
			}
		}
	}

	if (loccoords.Length() == elements.Length() && loccoords2.Length() == elements2.Length())
	{
		return ACT_Loccoords;
	}

	return ACT_None;
}

void AverageConstraint::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const AverageConstraint& ce = (const AverageConstraint&)e;

	nodes2			= ce.nodes2;
	loccoords2	= ce.loccoords2;
	elements2		= ce.elements2;
	weights1		= ce.weights1;
	weights2		= ce.weights2;
	sumWeights1	= ce.sumWeights1;
	sumWeights2	= ce.sumWeights2;
	polynomial_1 = ce.polynomial_1;
	polynomial_2 = ce.polynomial_2;

	moment_settings			= ce.moment_settings;
	stiffnesses	= ce.stiffnesses;
	elements_nodouble	= ce.elements_nodouble;
	is_nodal_constraint = ce.is_nodal_constraint;
	elements_share_dofs = ce.elements_share_dofs;

	scaling_factor = ce.scaling_factor;
}

void AverageConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	BasePointJoint::GetElementData(edc);
	ElementData ed;
	
	Matrix m;
	m.SetSize(loccoords.Length(),3);
	for(int i=1; i<=loccoords.Length();i++)
	{
		Vector v;
		v.Set3D(loccoords(i).X(),loccoords(i).Y(),loccoords(i).Z());
		m.SetRowVec(v,i);
	}
	int cols = 3;
	int rows = loccoords.Length();
	ed.SetMatrix(m.GetMatPtr(),rows,cols,"positions");
	ed.SetToolTipText("(local) positions of the points of kinematic pair 1");
	edc.TreeAdd("Position1",ed);

	m.SetSize(loccoords2.Length(),3);
	for(int i=1; i<=loccoords2.Length();i++)
	{
		Vector v;
		v.Set3D(loccoords2(i).X(),loccoords2(i).Y(),loccoords2(i).Z());
		m.SetRowVec(v,i);
	}
	cols = 3;
	rows = loccoords2.Length();
	ed.SetMatrix(m.GetMatPtr(),rows,cols,"positions");
	ed.SetToolTipText("(local) positions of the points of kinematic pair 2");
	edc.TreeAdd("Position2",ed);

	m.SetSize(moment_settings.Length(),3);
	for(int i=1; i<=moment_settings.Length();i++)
	{
		Vector v;
		v.Set3D(moment_settings(i).Get(1),moment_settings(i).Get(2),moment_settings(i).Get(3));
		m.SetRowVec(v,i);
	}
	cols = 3;
	rows = moment_settings.Length();
	ed.SetMatrix(m.GetMatPtr(),rows,cols,"moment_settings");
	ed.SetToolTipText("moment_settings(i) = [constrained coordinate, moment order, moment coordinate], e.g. rz*pow(y,4) = [3,4,2]; rx*pow(z,2) = [1,2,3]; ry*pow(s,3) = [2,3,0] (with 's': arc length); ry = [2,0,0]");
	edc.TreeAdd("Physics",ed);
}

int AverageConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = BasePointJoint::SetElementData(edc);

	int rows,cols;
	double *mp;

	if(edc.TreeGetMatrix("Position1.positions",&mp,rows,cols))
	{
		Matrix m(rows,cols,mp);
		loccoords.SetLen(rows);
		for(int i=1; i<=rows; i++)
		{
			loccoords(i) = Vector3D(m(i,1),m(i,2),m(i,3));
		}
	}

	if(edc.TreeGetMatrix("Position2.positions",&mp,rows,cols))
	{
		Matrix m(rows,cols,mp);
		loccoords2.SetLen(rows);
		for(int i=1; i<=rows; i++)
		{
			loccoords2(i) = Vector3D(m(i,1),m(i,2),m(i,3));
		}
	}

	if(edc.TreeGetMatrix("Physics.moment_settings",&mp,rows,cols))
	{
		Matrix m(rows,cols,mp);
		moment_settings.SetLen(rows);
		for(int i=1; i<=rows; i++)
		{
			moment_settings(i) = int3(m(i,1),m(i,2),m(i,3));
		}
	}

	return rv;
}

int AverageConstraint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	//BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 

	// Manual READ entries for this class

	// Manual WRITE entries for this class

	return 0;
}

int AverageConstraint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	//// call base class routine 	
	////int rv = BasePointJoint::ReadSingleElementData(RWdata);
	////if (rv == 1) return 1;

	//// manual things to read  
	//if(RWdata.variable_name.CStrCompare("Internal.actor_force") )
	//{
	//	RWdata.value = ComputeRopeForce(mbs->GetTime()); 
	//	return 1; 
	//}
	//if(RWdata.variable_name.CStrCompare("Internal.rope_length") )
	//{
	//	RWdata.value = GetLengthOfRope();
	//	return 1; 
	//}
	//if(RWdata.variable_name.CStrCompare("Internal.coiled_length") )
	//{
	//	RWdata.value = coiled_length;
	//	return 1; 
	//}

	//return ReadSingleElementDataAuto(RWdata);
	return 1;
}

int AverageConstraint::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	//// call base class routine ( not required in Element )	
	//int rv = BasePointJoint::WriteSingleElementData(RWdata);
	//if (rv == 1) return 1;
	////manual things to write
	//if(RWdata.variable_name.CStrCompare("Internal.coiled_length"))
	//{
	//	coiled_length = RWdata.value;
	//	return 1; 
	//}

	//return WriteSingleElementDataAuto(RWdata);
	return 1;
}

void AverageConstraint::SetAverageConstraint_LocalNodes_to_LocalNodes(TArray<int> element_numbers1, TArray<int> node_numbers1, Vector weights1_i, TArray<int> element_numbers2, TArray<int> node_numbers2, Vector weights2_i)
{

}
void AverageConstraint::SetAverageConstraint_LocalCoords_to_LocalCoords(TArray<int> element_numbers1, TArray<Vector3D> loccoords1_i, Vector weights1_i, TArray<int> element_numbers2, TArray<Vector3D> loccoords2_i, Vector weights2_i)
{
	elements = element_numbers1;
	loccoords = loccoords1_i;
	weights1 = weights1_i;

	elements2 = element_numbers2;
	loccoords2 = loccoords2_i;
	weights2 = weights2_i;

	SetElement_nodouble();
}
void AverageConstraint::SetAverageConstraint_LocalNodes_to_GlobalPos(TArray<int> element_numbers1, TArray<int> node_numbers1, Vector weights1_i, TArray<Vector3D> ground, Vector weights2_i)
{
	elements = element_numbers1;
	nodes = node_numbers1;
	weights1 = weights1_i;

	elements2.SetLen(weights2_i.Length());
	elements2.SetAll(0);
	loccoords2 = ground;
	weights2 = weights2_i;
	
	SetElement_nodouble();
}
void AverageConstraint::SetAverageConstraint_LocalCoords_to_GlobalPos(TArray<int> element_numbers1, TArray<Vector3D> loccoords1, Vector weights1_i, TArray<Vector3D> ground, Vector weights2_i)
{
	elements = element_numbers1;
	loccoords = loccoords1;
	weights1 = weights1_i;

	elements2.SetLen(weights2_i.Length());
	elements2.SetAll(0);
	loccoords2 = ground;
	weights2 = weights2_i;

	SetElement_nodouble();
}

//$ PG 2013-4-24: [ just for testing with existing model
void AverageConstraint::SetEdgeConstraint_Displacement(const IVector& element_numbers1, const TArray<Vector3D>& loccoords1, const Vector& weights, const TArray<int3>& moments_i, const Matrix& prescribed_displacements, const Vector& stiffnesses)
{
	// is input valid?
	int len=element_numbers1.Length();
	assert(len > 0);
	assert(loccoords1.Length() == len);
	assert(prescribed_displacements.Getrows() == 0 || prescribed_displacements.Getrows() == len);
	
	// define ground points based on element_numbers1, loccoords1, constrained_direction, and prescribed_displacements
	TArray<Vector3D> ground;
	ground.SetLen(len);
	for(int i=1; i<=ground.Length(); i++)
	{
		ground(i) = GetMBS()->GetElement(element_numbers1(i)).GetRefConfPos(loccoords1(i));
	}
	if (prescribed_displacements.Getrows() == len)
	{
		for(int i=1; i<=ground.Length(); i++)
		{
			ground(i) += Vector3D(prescribed_displacements(i,1), prescribed_displacements(i,2), 0.);
		}
	}

	// call set function (same weights for material and ground points)
	SetAverageConstraint_LocalCoords_to_GlobalPos(element_numbers1, loccoords1, weights, ground, weights);
	SetMoments(moments_i, stiffnesses); //Lagrange (Penalty would require stiffnesses for each of the moment_settings here)
}
//$ PG 2013-4-24: ]

void AverageConstraint::EvalF2(Vector& f, double t)
{
	// f = [f_el1, f_el2,... f_el1_vel, f_el2_vel]
	// eli is the i-th element in the list of elements_no_double

	if(!UsePenaltyFormulation())
		return;	// no penalty method --> Lagrange multiplier --> no EvalF2

	ConstVector<FEmaxDOF> f_el;
	int offset = 0;
	int sos;
	
	Vector momentum_times_stiffness(stiffnesses);
	for(int m=1; m<=moment_settings.Length();m++)
	{
		momentum_times_stiffness(m)*=ComputeMomentum(t,m);
	}

	// compute the f-vector f_el for each element
	// add f_el to the correct position in the f-vector of the constraint
	for (int locel=1; locel <= NE_nodouble(); locel++)
	{
		sos = GetElem_nodouble(locel).SOS();
		f_el.SetLen(sos);
		f_el.SetAll(0.);
		ComputeMomentumDerivative(t,locel,momentum_times_stiffness,f_el);	// compute f_el
		
		for(int j=1; j<=sos; j++)
		{
			f(offset+j) -= f_el(j);
		}
		offset+=sos;
	}

}

// needed for sensor
//Vector3D AverageConstraint::ComputeForce(double t, int moment_index) const
// returns force per unit length(area/volume/...  depends on definition of weights by user) at specified position and time
Vector3D AverageConstraint::ComputeForce(double t, int pos_nr, int kinpair) const
{
	Vector3D force(0.);

	if (UsePenaltyFormulation())
	{ 
		for (int moment_idx=1; moment_idx<=moment_settings.Length(); moment_idx++)
		{
			// f += c*M*f(pos)*v
			force += stiffnesses(moment_idx)*ComputeMomentum(t,moment_idx)*GetMomentumPolynomial(kinpair, pos_nr, moment_idx)*GetFrameVector(GetConstrainedDirection(moment_idx));
		}
	}
	else   // Lagrange
	{
		for (int moment_idx=1; moment_idx<=moment_settings.Length(); moment_idx++)
		{
			// f += lambda*f(pos)*v
			force += (XG(moment_idx)*GetMomentumPolynomial(kinpair, pos_nr, moment_idx))*GetFrameVector(GetConstrainedDirection(moment_idx));
		}
		force *= 1./scaling_factor;
	}

	return force;
}

//Vector3D AverageConstraint::ComputeForce(double t, int kinpair) const
//{
//	const TArray<int>& elements_kinpair = (kinpair==1) ? elements : elements2 ;
//	const Vector& weights = (kinpair==1) ? weights1 : weights2 ;
//	
//	Vector3D force(0.);
//	for(int i=1; i<=elements_kinpair.Length(); i++)
//	{
//		force += ComputeForce(t, kinpair, i)*weights(i);
//	}
//
//	return force;
//}

// rotation to the correct coordinate system (local frame)
// subtraction of the prescribed positions
// evaluation of the polynomial momenti
// returns "C", see definition of "C" in section "equation" of pdf user-docu
double AverageConstraint::ComputeMomentum(double t, int moment_index) const
{
	double val = 0;
	double component = 0;
	Vector3D pos;
	
	Vector3D dv = GetFrameVector(GetConstrainedDirection(moment_index));

	for(int i=1; i<=elements.Length(); i++)
	{
		const Body3D& b = GetBody3D(1,i);
		if(is_nodal_constraint)
		{
			pos = b.GetNodePos(nodes(i));
		}
		else								
		{
			pos = b.GetPos(loccoords(i)); 
		}
		
		component = pos*dv;

		// computation of the momenti of higher orders
		component *= GetMomentumPolynomial(1, i, moment_index);
		component *= weights1(i);
		val += component;
	}

	for (int i=1; i<=elements2.Length(); i++)
	{
		if (elements2(i))
		{
			const Body3D& b = GetBody3D(2,i);
			if (is_nodal_constraint)
			{
				pos = b.GetNodePos(nodes2(i));
			}
			else
			{
				pos = b.GetPos(loccoords2(i));
			}
		}
		else	// ground joint --> relative position
		{
			pos = loccoords2(i);
		}
		
		component = pos*dv;
		
		// computation of the momenti of higher orders
		component *= GetMomentumPolynomial(2, i, moment_index);
		component *= weights2(i);
		val -= component;
	}
	return val;
}

void AverageConstraint::ComputeMomentumDerivative(double t, int locelemind, const Vector& multiplier, Vector& f)
{	
	// this function computes the summation over all moment derivatives \nabla_q M_k(q)
	// the parameter "multiplier" is a Vector of lagrange multipliers or stiffness parameters (for penalty - one for each moment -- see member TArray<int2> moment_settings),
	

	// In case the the global directions shall be constrained, only the GetPos()-terms in ComputeMomentum are dependent on element-DOFs, i.e. the derivative
	// basically involves the computations of Body3D::GetdPosdqT(..) for each of the elements addressed by the array elements_nodouble.
	// Corresponding to these elements the local residual vector f is provided as an input-output argument.

	//
	// In case the local (and not the global) directions of the body(ies) shall be constrained (UseLocalCoordinateSystem() == 1), then another term  (~ Body3D::GetdRotvdqT, since the rotation depends on q of the first kinpair)
	// has to be added corresponding to the derivative of the constrained direction (see GetFrameVector(..)) w.r.t. the element dofs.


	int global_element_number = elements_nodouble(locelemind);

	ConstMatrix<FEmaxDOF*3> drotvdqT(GetMBS()->GetElementPtr(global_element_number)->SOS(), GetMBS()->GetElementPtr(global_element_number)->Dim(), 1);
	ConstMatrix<FEmaxDOF*3> dposdqT(GetMBS()->GetElementPtr(global_element_number)->SOS(), GetMBS()->GetElementPtr(global_element_number)->Dim(), 1);
	ConstVector<FEmaxDOF> temp_vector(GetMBS()->GetElementPtr(global_element_number)->SOS(), 1);

	for (int k=1; k<=NKinPairs(); k++)   //k= 1 or 2 of kinematic pair
	{
		TArray<int>& elements_k = (k==1) ? elements : elements2 ;
		TArray<Vector3D>& loccoords_k = (k==1) ? loccoords : loccoords2 ;
		TArray<int>& nodes_k = (k==1) ? nodes : nodes2 ;
		Vector& weights_k = (k==1) ? weights1 : weights2 ;
		int sign_k = (k==1) ? 1 : -1 ;

		for (int i=1; i<=elements_k.Length(); i++)    //loop over positions
		{
			if (elements_k.Get(i) == global_element_number)
			{
				Body3D& b = GetBody3D(k,i);
				
				if (is_nodal_constraint)
				{
					b.GetNodedPosdqT(nodes_k(i), dposdqT);
				}
				else								
				{
					b.GetdPosdqT(loccoords_k(i), dposdqT);
				}
				
				for (int m=1; m<=moment_settings.Length(); m++)    //loop over (direction, polynomial-order) - couples
				{
					if (fabs(multiplier(m))>1e-14)   // non-zero (numerically)
					{
						Mult(dposdqT, GetFrameVector(GetConstrainedDirection(m)), temp_vector);
						f += (temp_vector)*(sign_k*GetMomentumPolynomial(k, i, m)*weights_k(i)*multiplier(m));

						// if constrained direction d is dependent on element-dofs, then the product rule causes a second term
						if (UseLocalCoordinateSystem())
						{
							Vector3D pos;			

							// determine constant part e_m of rotation (according to AverageConstraint::GetFrameVector), and apply it to GetdRotvdqT as argument 'vloc'
							Vector3D e_m(0.);
							e_m(GetConstrainedDirection(m)) = 1.;

							if(stiffness_in_joint_local_frame)
							{
								e_m = e_m*JA0i;
							}
														
							if(is_nodal_constraint)
							{
								GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint::ComputeMomentumDerivative not yet implemented for case (UseLocalCoordinateSystem() && is_nodal_constraint)!\n";
							}
							else								
							{
								b.GetdRotvdqT(e_m, loccoords_k(i), drotvdqT);
								pos = b.GetPos(loccoords_k(i));
							}

							Mult(drotvdqT, pos, temp_vector);							
							f += (temp_vector)*(sign_k*GetMomentumPolynomial(k, i, m)*weights_k(i)*multiplier(m));
						}
					}
				}
			}
		}
	}
}

// get x^p (lever_direction=1), y^p (lever_direction=2), z^p (lever_direction=3), or \xi^p (lever_direction=0)
double AverageConstraint::GetMomentumPolynomial(const int kinpair, const int pos_nr, const int moment_index) const
{
	const Matrix& polynomial = (kinpair==1) ? polynomial_1 : polynomial_2;
	return polynomial(moment_index, pos_nr);
}

// computation of x^p (lever_direction=1), y^p (lever_direction=2), z^p (lever_direction=3), or \xi^p (lever_direction=0)
// PG will improve this matrices, such that the numerical solver is happy (orthogonalize moments --> integral of higher moments = 0)
void AverageConstraint::SetMomentumPolynomial()
{
	for (int kinpair=1; kinpair <= 2; kinpair++)
	{
		Matrix& polynomial = (kinpair==1) ? polynomial_1 : polynomial_2;
		Vector& weights = (kinpair==1) ? weights1 : weights2;
		polynomial.SetSize(moment_settings.Length(), NElements(kinpair));

		// precomputation of center and deviation of pointcloud (\max_{i,j} |(p_i-p_j) . d| where d is e1, e2, and e3 of frame (see GetFrameVector(i))
		Vector3D pos = GetRefConfPos(kinpair, 1);
		Vector3D pointcloud_center(pos);
		Vector3D pointcloud_size_max;
		for (int i=1; i<=3; i++)
		{
			pointcloud_size_max(i) = pos*GetFrameVector(i);
		}
		Vector3D pointcloud_size_min(pointcloud_size_max);
		for (int pos_nr=2; pos_nr <= NElements(kinpair); pos_nr++)
		{
			pos = GetRefConfPos(kinpair, pos_nr);
			pointcloud_center += pos;
			for (int i=1; i<=3; i++)
			{
				pointcloud_size_max(i) = max(pointcloud_size_max(i), pos*GetFrameVector(i));
				pointcloud_size_min(i) = min(pointcloud_size_min(i), pos*GetFrameVector(i));
			}
		}
		pointcloud_center *= 1./NElements(kinpair);  // arithmetic midpoint
		Vector3D pointcloud_size = pointcloud_size_max - pointcloud_size_min;

		// computation of polynomial values
		for (int pos_nr=1; pos_nr <= NElements(kinpair); pos_nr++)
		{
			for (int moment_index=1; moment_index <= moment_settings.Length(); moment_index++)
			{
				int moment_order = GetMomentOrder(moment_index);
				int lever_direction = GetLeverDirection(moment_index);
				if (moment_order == 0)
				{
					polynomial(moment_index, pos_nr) = 1.;
				}
				else
				{
					double xi;   // position in lever direction (if lever_direction != 0) or arc length along curve (if lever_direction == 0)
					if (lever_direction == 0)   // assuming trapecoid rule for back-calculation of xi from weights1/2
					{
						// first compute xi out of the information of weights and pos_nr
						// additional information: xi e [-1,1]

						if (pos_nr == 1)
						{
							xi = -1;
						}
						else if (pos_nr == weights.Length())
						{
							xi = 1;
						}
						else
						{
							xi = weights.Sum(pos_nr)-weights(pos_nr)/2;	// at this moment: xi e [0,sumWeights1]
							xi = 2*xi/sumWeights1 - 1;	// xi e ]-1,1[
						}
					}		
					else // lever direction is set to 1, 2, or 3
					{
						if (pointcloud_size(lever_direction))
						{
						Vector3D pos = GetRefConfPos(kinpair, pos_nr);
						xi = ((pos - pointcloud_center) * GetFrameVector(lever_direction)) * 2. / pointcloud_size(lever_direction);
						}
						else
						{
							GetMBS()->UO(UO_LVL_err) << "ERROR: either just one point specified, or point cloud is plane in lever direction " << lever_direction << "!\n";
							assert(0);
						}
					}

					// then compute the polynomial
					polynomial(moment_index, pos_nr) = pow(xi, moment_order);
				}
			}
		}

		GetMBS()->UO(UO_LVL_dbg1) << "kinpair = " << kinpair << ", polynomial = " << polynomial << "\n";
	}
}

Vector3D AverageConstraint::GetRefConfPos(int kinpair, int pos_nr) const
{
	// return position of a material point (pos_nr) of body 1/2 (kinpair) in reference configuration
	const IVector& element_array = (kinpair == 1) ? elements : elements2;
	const TArray<Vector3D>& loccoords_array = (kinpair == 1) ? loccoords : loccoords2;
	if (element_array(pos_nr) == 0)  // ground
	{
		return loccoords_array(pos_nr);
	}

	return GetBody3D(kinpair, pos_nr).GetRefConfPos(loccoords_array(pos_nr));
}

Vector3D AverageConstraint::GetFrameVector(int i) const
{
	Matrix3D A = GetRotMati();
	return Vector3D(A(1,i),A(2,i),A(3,i));
}

void AverageConstraint::AddElementCqTLambda(double t, int locelemind, Vector&f)
{
// we calculate \sum_m lambda_m (\nabla M_m(q)) with respect to element-DOFs:
	ConstVector<MAX_IS> lambda(IS(),1);
	for (int m=1; m<=IS(); m++) 
	{
		lambda(m) = XG(m);
	}

	ComputeMomentumDerivative(t, locelemind, lambda, f);
	f *= scaling_factor;
	//mbs->UO() << "locelem " << locelemind << ": " << f << "\n";
	//$ PG 2013-4-26: suggest to rename ComputeMomentumDerivative --> ComputeSumOfMomentumDerivative
}


void AverageConstraint::EvalG(Vector &f, double t)
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG

	if (MaxIndex() == 3)
	{
		for(int m=1; m <= moment_settings.Length(); m++)
		{
			f(m) = ComputeMomentum(t, m);
		}
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: in AverageConstraint: MaxIndex() != 3 not implemented yet!\n";
	}
	f *= scaling_factor;
}

// old
//void AverageConstraint::LinkToElementsPenalty()
//{
//	LTGreset();
//
//	int nElems;
//	int sos;
//	// add all SOS dofs from the elements
//	//Position(first SOS) 
//
//	for (int k=1; k <= NKinPairs(); k++)
//	{
//		if(k==1) nElems = elements.Length();
//		else nElems = elements2.Length();
//
//		for(int i=1; i<=nElems; i++)
//		{
//			sos = GetElement(k,i).SOS();
//			for (int j=1; j <= sos; j++)
//			{
//				AddLTG(GetElement(k,i).LTG(j));
//			}
//		}
//	}
//
//	//and Velocity (second SOS):
//	for (int k=1; k <= NKinPairs(); k++)
//	{
//		if(k==1) nElems = elements.Length();
//		else nElems = elements2.Length();
//
//		for(int i=1; i<=nElems; i++)
//		{
//			sos = GetElement(k,i).SOS();
//			for (int j=1; j <= sos; j++)
//			{
//				AddLTG(GetElement(k,i).LTG(j+sos));
//			}
//		}
//	}
//};


void AverageConstraint::LinkToElementsPenalty()
// f = [f_el1, f_el2,... f_el1_vel, f_el2_vel]
// eli is the i-th element in the list of elements_no_double
{
	LTGreset();
	int sos;

	// add all SOS dofs from the elements
	//Position(first SOS) 
	for(int i=1; i<=NE_nodouble(); i++)
	{
		sos = GetElem_nodouble(i).SOS();
		for (int j=1; j <= sos; j++)
		{
			AddLTG(GetElem_nodouble(i).LTG(j));
		}
	}

	//and Velocity (second SOS):
	for(int i=1; i<=NE_nodouble(); i++)
	{
		sos = GetElem_nodouble(i).SOS();
		for (int j=1; j <= sos; j++)
		{
			AddLTG(GetElem_nodouble(i).LTG(j+sos));
		}
	}
};

void AverageConstraint::LinkToElementsLagrange()
{
	if (IsLocalNodeConstraint())
	{
		//$ PG 2013-4-24: [
		for (int i = 1; i <= NE_nodouble(); i++)
		{
			GetElem_nodouble(i).AddConstraint(this, i);
		}
		// particularly needed for evaluation of AverageConstraint::AddElementCqTLambda(..., int locelemindex, ...)
		// which is called in Element::EvalF2 via
		// for (int i=1; i <= constraints.Length(); i++)
		// {
		//   constraints(i)->AddElementCqTLambda(t, constraintindices(i),f);
		// }
		// where TArray<Constraint*> constraints and TArray<int> constraintindices are both assembled through calling
		// Element::AddConstraint(Constraint* c, int index), which we do right here!
		//$ PG 2013-4-24: ]
	}
}

void AverageConstraint::DrawElement() 
{
	Constraint::DrawElement();

	Vector3D pos;
	int res = 8;
	mbs->SetColor(Vector3D(0.,0.,0.8));
	for(int i=1; i<=elements.Length(); i++)
	{
		
		if(is_nodal_constraint){ pos = mbs->GetElementPtr(elements(i))->GetNodePosD(nodes(i)); }
		else								{ pos = mbs->GetElementPtr(elements(i))->GetPosD(loccoords(i)); }
		mbs->DrawSphere(pos,GetDrawSizeScalar(),res);
	}

	mbs->SetColor(Vector3D(0.8,0.1,0.1)); 
	res = 7;
	for(int i=1; i<=elements2.Length(); i++)
	{
		if(elements2(i)) 
		{
			if(is_nodal_constraint){ pos = mbs->GetElementPtr(elements2(i))->GetNodePosD(nodes2(i)); }
			else								{ pos = mbs->GetElementPtr(elements2(i))->GetPosD(loccoords2(i)); }
		}
		else
		{
			pos = loccoords2(i);	// ground
		}
		mbs->DrawSphere(pos,GetDrawSizeScalar(),res);
	}

	if(draw_local_frame_size > 0 || draw_local_frame_size == -1)
	{
		double s = draw_local_frame_size;
		if(s == -1) {s = GetMBS()->GetDOption(104); }	// use default value

		GetMBS()->ChooseColor(0.3f,0.3f,0.3f);

		// size of the frame
		Vector3D v1(-0.2*s, 0,0);
		Vector3D v2( s, 0,0);
		Vector3D v3( 0,-0.2*s,0);
		Vector3D v4( 0, s,0);
		Vector3D v5( 0,0,-0.2*s);
		Vector3D v6( 0,0, s);

		Matrix3D A = GetRotMatiD();
		Vector3D p;
		if(is_nodal_constraint){ p = mbs->GetElementPtr(elements(1))->GetNodePosD(nodes(1)); }
		else								{ p = mbs->GetElementPtr(elements(1))->GetPosD(loccoords(1)); }

		v1 = p + A*v1;
		v2 = p + A*v2;
		v3 = p + A*v3;
		v4 = p + A*v4;
		v5 = p + A*v5;
		v6 = p + A*v6;

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
};