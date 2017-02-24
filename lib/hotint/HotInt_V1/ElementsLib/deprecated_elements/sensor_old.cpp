//#**************************************************************
//#
//# filename:             sensor_old.cpp
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


#ifndef __SENSOR_OLD_CPP__
#define __SENSOR_OLD_CPP__


#include "element.h"
#include "body2d.h"
#include "body3d.h"
#include "node.h"
#include "sensor_old.h"
#include "constraint.h"
#include "control.h"
#include "fft_utilities.h"


//$ YV 2013-01-02: moved here from optimization.cpp
int FFToptimization_averaging=0; //0=no averaging, 1, 2, n.... is average over -n ... +n values (0=no averaging, 1=2 values, 2=5values)


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//MBSSensor:

//$ YV 2012-06: the sensors may produce just one scalar value
double MBSSensor::GetCurrentValue(double time)
{
	if(type & TSFile)
	{
		if(tsfileValues)
		{
			actvalue = tsfileValues->Evaluate(time); // attention: this is only at start or end of solver time step correct due to use of double t = GetMBS()->GetTime();
		}
	}

	if (type & TSDOF)
	{
		if (type & TSVel) //$JG 18-11-2011
		{
			actvalue = GetMBS()->GetElement(elementlist(1)).XGP(nodelist(1));
		}
		else
		{
			actvalue = GetMBS()->GetElement(elementlist(1)).XG(nodelist(1));
		}
	}

	if (type & TSXData)		//$ AH 27-11-2011
	{
		actvalue = GetMBS()->GetElement(elementlist(1)).XData(nodelist(1));
	}

	if (type & TSConstraintDrift)	//$ JG 2011-02
	{
		//node number not used!!!
		actvalue = GetMBS()->GetElement(elementlist(1)).GetConstraintDrift(time);
	}

	if (type & TSAxialExtension)
	{
		if(type & TSElement)
		{
			Body2D *e2D;
			Body3D *e3D;
			if (type & TSplanar)
			{
				e2D = (Body2D*)GetMBS()->GetElementPtr(elementlist(1));
				actvalue = e2D->GetFieldVariableValue(FieldVariableDescriptor::FVT_beam_axial_extension, poslist2D(1), false);
			}
			else
			{
				e3D = (Body3D*)GetMBS()->GetElementPtr(elementlist(1));
				actvalue = e3D->GetFieldVariableValue(FieldVariableDescriptor::FVT_beam_axial_extension, poslist(1), false);
			}
		}
		else
		{
			GetMBS()->UO() << "MBSSensor::TSAxialExtension: not implemented yet for nodes!\n";
		}
	}

	if (type & TSForce)
	{
		if ((GetMBS()->GetElement(elementlist(1)).GetType() & TController) ||
			(GetMBS()->GetElement(elementlist(1)).GetType() & TConstraint))  // DR 2011-11-17: changed "==" to "&"
		{
			const Constraint& c = (const Constraint&)GetMBS()->GetElement(elementlist(1));
			int dir = 0;
			if (type & TSX) dir = 1;
			if (type & TSY) dir = 2;
			if (type & TSZ) dir = 3;

			actvalue = c.GetActorForce(time, dir);
		}
	}

	if (type & TSSpecialSensorValue)		// DR 2012-01-12: arbitrary sensor value working for all elements
	{
		if(type & TSElement)
		{
			assert(elementlist.Length() == 1);        // single element 
			assert(nodelist.Length() == 1);           // single special value

			Element* el = (Element*) GetMBS()->GetElementPtr(elementlist(1)); 
			actvalue = el->GetSpecialSensorValue(nodelist(1), time);
		}
	}

	if (type & TSOutputSensor)
	{
		const InputOutputElement& ioe = (const InputOutputElement&)GetMBS()->GetElement(elementlist(1));
		actvalue = ioe.GetOutput(time, nodelist(1));
	}

	if (type & TSAccel)
	{
		if(type & TSElement)//$ RL 2011-3-15:[ measurement of components of global / local acceleration vector; note: TSAccelDOF renamed to TSAccel.
		{
			assert(elementlist.Length() == 1);
			const Element& e = GetMBS()->GetElement(elementlist(1));
			Matrix3D AT(1.); //unit matrix
			if (type & TSLocalAxis)
			{
				AT = e.GetRotMatrix(poslist(1)).GetTp(); // ploc = AT*pglob
			}
			if (type & TSX)
			{
				actvalue = (AT*e.GetAcceleration(poslist(1))).X();
			}
			else if (type & TSY)
			{
				actvalue = (AT*e.GetAcceleration(poslist(1))).Y();
			}
			else if (type & TSZ)
			{
				actvalue = (AT*e.GetAcceleration(poslist(1))).Z();
			}
			else
			{
				//actvalue = (AT*e.GetAcceleration(poslist(1))).Norm();
				actvalue = (e.GetAcceleration(poslist(1))).Norm(); //$ RL 2011-3-16:  rotation matrix doesn't change length of vector! --> faster computation
			}
		}
		else if(type & TSDOF) //$ RL 2011-3-15:] measurement of components of global / local acceleration vector; note: TSAccelDOF renamed to TSAccel.
		{
			const Element& e = GetMBS()->GetElement(elementlist(1));
			int acc_i = nodelist(1);
			//int k = e.LTG(acc_i + e.SOS()); // $ MSax 2013-07-16 : removed
			//actvalue = GetMBS()->GetAcceleration(k); // $ MSax 2013-07-16 : removed
			actvalue = e.XGPP(acc_i); // $ MSax 2013-07-16 : added
		}
		else//$ RL 2011-3-15:[ 
		{
			assert(0 && "sensor type not known (use TSDOF or TSElement)");
		}//$ RL 2011-3-15:] 

	}

	if (type & (TSStress | TSStrain | TSInelStrain))
	{
		//if ((type & TSAuxElem))
		//{
		//}
		//else
		//{ 
		if (type&TSElement || type&TSAuxElem)
		{
			Body2D* e2D;
			Body3D* e3D;
			if (!(type & TSplanar) ) 
			{
				//3D:
				if (type&TSElement && (elementlist(1) < 0 || elementlist(1) > GetMBS()->NE() || 
					GetMBS()->GetElementPtr(elementlist(1))->Dim() != 3 || 
					!GetMBS()->GetElementPtr(elementlist(1))->IsType(TBody)))
				{
					GetMBS()->UO() << "Error in Sensor: sensor element number is not a valid 3D body!!!\n";
					return actvalue;
				}
				if (type&TSAuxElem && (elementlist(1) < 0 || elementlist(1) > GetMBS()->NAuxE() || 
					GetMBS()->GetAuxElementPtr(elementlist(1))->Dim() != 3 || 
					!GetMBS()->GetAuxElementPtr(elementlist(1))->IsType(TBody)))
				{
					GetMBS()->UO() << "Error in Sensor: sensor auxiliary element number is not a valid 3D body!!!\n";
					return actvalue;
				}
				if (type&TSAuxElem) e3D = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(1));
				else e3D = (Body3D*)GetMBS()->GetElementPtr(elementlist(1));

			}
			else
			{
				//2D:
				if (type&TSElement && (elementlist(1) < 0 || elementlist(1) > GetMBS()->NE() || 
					GetMBS()->GetElementPtr(elementlist(1))->Dim() != 2 || 
					!GetMBS()->GetElementPtr(elementlist(1))->IsType(TBody)) )
				{
					GetMBS()->UO() << "Error in Sensor: sensor element number is not a valid 2D body!!!\n";
					return actvalue;
				}
				if (type&TSAuxElem && (elementlist(1) < 0 || elementlist(1) > GetMBS()->NAuxE() || 
					GetMBS()->GetAuxElementPtr(elementlist(1))->Dim() != 2 || 
					!GetMBS()->GetAuxElementPtr(elementlist(1))->IsType(TBody)) )
				{
					GetMBS()->UO() << "Error in Sensor: sensor element number is not a valid 2D body!!!\n";
					return actvalue;
				}
				if (type&TSAuxElem) e2D = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(1));
				else e2D = (Body2D*)GetMBS()->GetElementPtr(elementlist(1));
			}


			if (type & TSElement || type&TSAuxElem)
			{
				//				GetMBS()->UO() << "in TSElement!!!!!!!!!!!" << "\n";
				FieldVariableDescriptor::FieldVariableType variable_type = FieldVariableDescriptor::FVT_problem_specific;
				if (type & TSStress)
				{
					variable_type = FieldVariableDescriptor::FVT_stress;
					if(tensor_comp == 0)
						variable_type = FieldVariableDescriptor::FVT_stress_mises;
				}
				else if (type & TSStrain)
					variable_type = FieldVariableDescriptor::FVT_total_strain;
				else if (type & TSInelStrain)
					variable_type = FieldVariableDescriptor::FVT_inelastic_strain;

				if(variable_type != FieldVariableDescriptor::FVT_problem_specific)
				{
					// some type was evaluated
					FieldVariableDescriptor::FieldVariableComponentIndex component_index_1;
					FieldVariableDescriptor::FieldVariableComponentIndex component_index_2;
					switch(tensor_comp)
					{
					case 0: component_index_1 = FieldVariableDescriptor::FVCI_none; component_index_2 = FieldVariableDescriptor::FVCI_none; break;
					case 1: component_index_1 = FieldVariableDescriptor::FVCI_x; component_index_2 = FieldVariableDescriptor::FVCI_x; break;
					case 2: component_index_1 = FieldVariableDescriptor::FVCI_y; component_index_2 = FieldVariableDescriptor::FVCI_y; break;
					case 3: component_index_1 = FieldVariableDescriptor::FVCI_z; component_index_2 = FieldVariableDescriptor::FVCI_z; break;
					case 4: component_index_1 = FieldVariableDescriptor::FVCI_y; component_index_2 = FieldVariableDescriptor::FVCI_z; break;
					case 5: component_index_1 = FieldVariableDescriptor::FVCI_x; component_index_2 = FieldVariableDescriptor::FVCI_z; break;
					case 6: component_index_1 = FieldVariableDescriptor::FVCI_x; component_index_2 = FieldVariableDescriptor::FVCI_y; break;
					default: assert(0);
					}

					FieldVariableDescriptor fvd(variable_type, component_index_1, component_index_2);
					// ASTRID: replaced e23D->GetPos(poslist(1)) by poslist(1) below:
					//         then, local positions are used in element-routine GetFieldVariableValue
					//         otherwise, global positions were used, which did not fit with element routine
					if( !(type & TSplanar) )
					{
						actvalue = e3D->GetFieldVariableValue(fvd, /*e3D->GetPos(*/poslist(1)/*)*/, false);
					}
					else
					{
						actvalue = e2D->GetFieldVariableValue(fvd, /*e2D->GetPos2D(*/poslist2D(1)/*)*/, false);
					}

					//GetMBS()->UO() << "tensorcomp = " << tensor_comp << ", comp = " << comp << "\n";
					//GetMBS()->UO() << "actvalue = " << actvalue << "\n";
				}
			}
			else if (type & TSNode)
			{
				//missing
				GetMBS()->UO() << "MBSSensor::TensorComponents: not implemented yet for nodes!\n";
			}
		}
	}

	if (!(type & TSDOF)) //$JG 18-11-2011: in order to allow (TSDOF+TSvel)
	{
		//pure nodal sensor: measure position, displacement, velocity
		if ((type & TSNode) && !(type & TSElement) && elementlist.Length() == 0)
		{
			int dir = 0;
			if (type & TSX) dir = 1;
			if (type & TSY) dir = 2;
			if (type & TSZ) dir = 3;

			if (dir)
			{
				if (type & TSPos)
				{
					actvalue = GetMBS()->GetNode(nodelist(1)).GetPos()(dir);
				}
				else if (type & TSVel)
				{
					actvalue = GetMBS()->GetNode(nodelist(1)).GetVel()(dir);
				}
				else if (type & TSDisplacement)
				{
					actvalue = GetMBS()->GetNode(nodelist(1)).GetDisplacement()(dir);
				}
			}
			else 
				actvalue = 0;
		}
		else if (type & (TSPos | TSVel | TSAngle | TSKardanAngle | TSCurvature | TSDeflection | TSDist))
		{
			if (!(type & TSplanar))
			{ //3D:
				Body3D* e;
				if (!(type & TSAuxElem) ) 
				{
					if (elementlist(1) < 0 || elementlist(1) > GetMBS()->NE() || 
						GetMBS()->GetElementPtr(elementlist(1))->Dim() != 3 || 
						!GetMBS()->GetElementPtr(elementlist(1))->IsType(TBody))
					{
						GetMBS()->UO() << "Error in Sensor: sensor element number is not a valid 3D body!!!\n";
						return actvalue;
					}
					e = (Body3D*)GetMBS()->GetElementPtr(elementlist(1));
				}
				else e = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(1));

				if ((type & TSElement || type&TSAuxElem) && !(type & TSNode))
				{
					Matrix3D A(1.); //unit matrix
					if (type & TSLocalAxis)
					{
						if (!(type & TSOtherLocalAxisElement))
						{
							A = e->GetRotMatrix(poslist(1)).GetTp();
						}
						else if (elementlist.Length() >= 2)
						{

							Body3D* e2;
							if (type&TSElement) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
							if (type&TSAuxElem) e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));
							A = e2->GetRotMatrix(poslist(2)).GetTp();
						}
					}
					if (type & TSPos)
					{
						if (type & TSX)
							actvalue = (A*e->GetPos(poslist(1))).X();
						else if (type & TSY)
							actvalue = (A*e->GetPos(poslist(1))).Y();
						else if (type & TSZ)
							actvalue = (A*e->GetPos(poslist(1))).Z();
						else //no xyz
							actvalue = (A*e->GetPos(poslist(1))).Norm();
					}
					else if ((type & TSVel) && !(type & TSAngle))
					{
						if (type & TSX)
							actvalue = (A*e->GetVel(poslist(1))).X();
						else if (type & TSY)
							actvalue = (A*e->GetVel(poslist(1))).Y();
						else if (type & TSZ)
							actvalue = (A*e->GetVel(poslist(1))).Z();
						else //no xyz
							actvalue = (A*e->GetVel(poslist(1))).Norm();
					}
					else if ((type & TSAngle) && !(type & TSVel))
					{
						//missing
						//GetMBS()->UO() << "MBSSensor::Angle: not implemented yet!\n";
						//get rotation with respect to a certain axis:

						/*if (strcmp(e->GetElementSpec(),"Rigid3D") != 0)
						{
						GetMBS()->UO() << "Error in Angle sensor: sensor only implemented for Rigid3D!!!\n";
						return actvalues;
						}*/

						//depreciated: axis!!! --> axis of rotation
						int axis = 0; //x-axis
						if (type & TSX) axis = 1;
						else if (type & TSY) axis = 2;
						else if (type & TSZ) axis = 3;

						if (axis)
						{
							//depreciated: axis of rotation must be specified!, TSX, TSY, TSZ not valid!!!
							//only in Plate3D example with static deformation (with Arend)
							Matrix3D A = e->GetRotMatrix(poslist(1));
							actvalue = GetRotation(axis, A);
						}
						else
						{
							//double oldangle = actvalue;   
							double oldangle = (actvalue - offset)/factor;

							//TSAngle sensor with 4 vectors: position in Body1, rotation axis (global/local), body-fixed orientation and global orientation vector
							if (elementlist.Length() == 1)
							{
								Matrix3D A1 = e->GetRotMatrix(poslist(1));
								Vector3D rot = poslist(2);
								if (type & TSLocalAxis) rot = A1*rot;

								Vector3D vec1 = A1*poslist(3); //body-fixed orientation (reference) vector
								Vector3D vec2 = poslist(4); //global orientation (reference) vector

								ProjectInPlane(Vector3D(0.,0.,0.), rot, vec1);
								ProjectInPlane(Vector3D(0.,0.,0.), rot, vec2);
								vec1.Normalize();
								vec2.Normalize();

								actvalue = NormalizedVectorAngle(vec1,vec2);

								if (rot*(vec1.Cross(vec2)) < 0) actvalue = -actvalue;
							}
							//TSAngle sensor with 5 vectors: position in Body1, position in Body2, rotation axis (global/local in Body1), body1-fixed orientation vector, body2-fixed orientation vector
							else if (elementlist.Length() == 2)
							{
								Body3D* e2;
								if (type&TSElement) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
								else if (type&TSAuxElem) e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));

								Matrix3D A1 = e->GetRotMatrix(poslist(1));
								Matrix3D A2 = e2->GetRotMatrix(poslist(2));
								Vector3D rot = poslist(3);
								if (type & TSLocalAxis)	rot = A1*rot;

								Vector3D vec1 = A1*poslist(4); //body1-fixed orientation (reference) vector
								Vector3D vec2 = A2*poslist(5); //body2-fixed orientation (reference) vector

								ProjectInPlane(Vector3D(0.,0.,0.), rot, vec1);
								ProjectInPlane(Vector3D(0.,0.,0.), rot, vec2);
								vec1.Normalize();
								vec2.Normalize();

								actvalue = NormalizedVectorAngle(vec1,vec2);

								if (rot*(vec1.Cross(vec2)) < 0) actvalue = -actvalue;


							}

							if (fabs(actvalue - oldangle) > MY_PI)
							{
								int fact = (int) (actvalue - oldangle)/(2*MY_PI);
								actvalue=actvalue - (fact + Sgn<double>(actvalue - oldangle))*2*MY_PI;
							}
						}
					}
					else if ((type & TSAngle) && (type & TSVel))
					{
						/*
						if (strcmp(e->GetElementSpec(),"Rigid3D") != 0)
						{
						GetMBS()->UO() << "Error in Angle sensor: sensor only implemented for Rigid3D!!!\n";
						return actvalues;
						}*/

						if (elementlist.Length() == 1)
						{
							if (type & TSLocalAxis)
							{
								actvalue = e->GetAngularVel(poslist(1)) * (e->GetRotMatrix()* poslist(2)); //project angular vel. into co-rotated axis of rotation
							}
							else
							{
								actvalue = e->GetAngularVel(poslist(1)) * (poslist(2)); //project angular vel. into axis of rotation
							}
						}
						else if (elementlist.Length() == 2) 
						{
							Body3D* e2;
							if (type&TSElement) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
							else if (type&TSAuxElem) e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));
							Vector3D rot = poslist(3);
							if (type & TSLocalAxis)	rot = e->GetRotMatrix(poslist(1))*rot;						
							actvalue = (e->GetAngularVel(poslist(1)) - e2->GetAngularVel(poslist(2))) * rot;  //project relative angular vel. into axis of rotation
						}
					}
					else if (type & TSKardanAngle)
					{
						int axis = 3; //z-axis is standard
						if (type & TSX) axis = 1;
						else if (type & TSY) axis = 2;
						else if (type & TSZ) axis = 3;

						Matrix3D A = e->GetRotMatrix(poslist(1));
						// compute quaternions
						double b0, b1, b2, b3;
						RotMatToQuaternions(A, b0, b1, b2, b3);
						// compute kardan angles
						Vector3D phi;
						QuaternionsToKardanAngles(b0, b1, b2, b3, phi);
						actvalue = phi(axis);		
					}
					else if (type & TSDist)
					{
						Body3D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						actvalue = (e->GetPos(poslist(1)) - e2->GetPos(poslist(2))).Norm();
					}
					else if (type & TSDeflection)
					{
						Body3D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						Body3D* e3;
						if (!(type & TSAuxElem) ) e3 = (Body3D*)GetMBS()->GetElementPtr(elementlist(3));
						else e3 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(3));

						actvalue = Deflection(e->GetPos(poslist(1)), e2->GetPos(poslist(2)), e3->GetPos(poslist(3)), poslist(4));
					}
				}
				else if (type & TSNode)
				{
					if (type & TSPos)
					{
						if (type & TSX)
							actvalue = e->GetNodePos(nodelist(1)).X();
						else if (type & TSY)
							actvalue = e->GetNodePos(nodelist(1)).Y();
						else if (type & TSZ)
							actvalue = e->GetNodePos(nodelist(1)).Z();
					}
					else if (type & TSVel)
					{
						if (type & TSX)
							actvalue = e->GetNodeVel(nodelist(1)).X();
						else if (type & TSY)
							actvalue = e->GetNodeVel(nodelist(1)).Y();
						else if (type & TSZ)
							actvalue = e->GetNodeVel(nodelist(1)).Z();
					}
					else if (type & TSAngle)
					{
						actvalue = 0;
						//missing
						GetMBS()->UO() << "MBSSensor::Angle: not implemented yet for nodes!\n";
					}
					else if (type & TSDist)
					{
						Body3D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));

						actvalue = (e->GetNodePos(nodelist(1)) - e2->GetNodePos(nodelist(2))).Norm();
					}
					else if (type & TSDeflection)
					{
						Body3D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body3D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						Body3D* e3;
						if (!(type & TSAuxElem) ) e3 = (Body3D*)GetMBS()->GetElementPtr(elementlist(3));
						else e3 = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(3));

						actvalue = Deflection(e->GetNodePos(nodelist(1)), e2->GetNodePos(nodelist(2)), 
							e3->GetNodePos(nodelist(3)), poslist(1));
					}
					//$ PG 2011-4-11: moved TSKardanAngle from TSNode (here) to TSElement
					//else if (type & TSKardanAngle)
					//{
					//	int axis = 3; //z-axis is standard
					//	if (type & TSX) axis = 1;
					//	else if (type & TSY) axis = 2;
					//	else if (type & TSZ) axis = 3;

					//	Matrix3D A = e->GetRotMatrix(poslist(1));
					//	// compute quaternions
					//	double b0, b1, b2, b3;
					//	RotMatToQuaternions(A, b0, b1, b2, b3);
					//	// compute kardan angles
					//	Vector3D phi;
					//	QuaternionsToKardanAngles(b0, b1, b2, b3, phi);
					//	actvalue = phi(axis);		
					//}
				}
			}
			else
			{ //2D:
				Body2D* e;
				if (!(type & TSAuxElem) ) 
				{
					if (elementlist(1) < 0 || elementlist(1) > GetMBS()->NE() || 
						GetMBS()->GetElementPtr(elementlist(1))->Dim() != 2 || 
						!GetMBS()->GetElementPtr(elementlist(1))->IsType(TBody))
					{
						GetMBS()->UO() << "Error in Sensor: sensor element number is not a valid 3D rigid body!!!\n";
						return actvalue;
					}
					e = (Body2D*)GetMBS()->GetElementPtr(elementlist(1));
				}
				else e = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(1));

				if (type & TSElement)
				{
					if (type & TSPos)
					{
						if (type & TSX)
							actvalue = e->GetPos2D(poslist2D(1)).X();
						else if (type & TSY)
							actvalue = e->GetPos2D(poslist2D(1)).Y();
						else if (type & TSZ) // Component Z --> compute absolute value of position vector
							actvalue = e->GetPos2D(poslist2D(1)).Norm();
					}
					else if ((type & TSVel) && !(type & TSAngle) && !(type & TSCurvature))
					{
						if (type & TSX)
							actvalue = e->GetVel2D(poslist2D(1)).X();
						else if (type & TSY)
							actvalue = e->GetVel2D(poslist2D(1)).Y();
						else if (type & TSZ) // Component Z --> compute absolute value of velocity vector
							actvalue = e->GetVel2D(poslist2D(1)).Norm();
					}
					else if (type & TSAngle && !(type & TSVel))
					{
						actvalue = e->GetAngle2D(poslist2D(1));
					}
					else if ((type & TSAngle) && (type & TSVel))
					{
						actvalue = e->GetAngle2DP(poslist2D(1));
					}
					else if (type & TSCurvature && !(type & TSVel))
					{
						actvalue = e->GetCurvature2D(poslist2D(1));
					}
					else if (type & TSCurvature && (type & TSVel))
					{
						actvalue = e->GetCurvature2DP(poslist2D(1));
					}
					else if (type & TSDist)
					{
						Body2D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body2D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(2));

						actvalue = (e->GetPos2D(poslist2D(1)) - e2->GetPos2D(poslist2D(2))).Norm();
					}
					else if (type & TSDeflection)
					{
						Body2D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body2D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						Body2D* e3;
						if (!(type & TSAuxElem) ) e3 = (Body2D*)GetMBS()->GetElementPtr(elementlist(3));
						else e3 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(3));

						actvalue = Deflection(e->ToP3D(e->GetPos2D(poslist2D(1))), e2->ToP3D(e2->GetPos2D(poslist2D(2))), 
							e3->ToP3D(e3->GetPos2D(poslist2D(3))), Vector3D(0,0,1));
					}
				}
				else if (type & TSNode)
				{
					if (type & TSPos)
					{
						if (type & TSX)
							actvalue = e->GetNodePos2D(nodelist(1)).X();
						else if (type & TSY)
							actvalue = e->GetNodePos2D(nodelist(1)).Y();
					}
					else if (type & TSVel)
					{
						if (type & TSX)
							actvalue = e->GetNodeVel2D(nodelist(1)).X();
						else if (type & TSY)
							actvalue = e->GetNodeVel2D(nodelist(1)).Y();
					}
					else if (type & TSAngle)
					{
						//missing
						GetMBS()->UO() << "MBSSensor::Angle: not implemented yet for nodes!\n";
					}
					else if (type & TSDist)
					{
						Body2D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body2D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						actvalue = (e->GetNodePos2D(nodelist(1)) - e2->GetNodePos2D(nodelist(2))).Norm();
					}
					else if (type & TSDeflection)
					{
						Body2D* e2;
						if (!(type & TSAuxElem) ) e2 = (Body2D*)GetMBS()->GetElementPtr(elementlist(2));
						else e2 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(2));
						Body2D* e3;
						if (!(type & TSAuxElem) ) e3 = (Body2D*)GetMBS()->GetElementPtr(elementlist(3));
						else e3 = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(3));

						actvalue = Deflection(e->ToP3D(e->GetNodePos2D(nodelist(1))), e2->ToP3D(e2->GetNodePos2D(nodelist(2))), 
							e3->ToP3D(e3->GetNodePos2D(nodelist(3))), Vector3D(0,0,1));
					}
				}
			}
		}
		else if (type & TSEigenValue)
		{
			assert(nodelist.Length() > 0);
			actvalue = GetMBS()->GetEigenValue(nodelist(1));
		}
	}
	//$ RL 2011-7-14: FFT code
	else if(type & TSFFT)
	{
		// sensor value does not change (fft is post processed after sensor-files are written, result is in sensor file with postfix "-fft")		
	}

	//$ AD: just testing sensors on Load-factor // parameter
	if (type & TSLoad)
	{
		assert(elementlist.Length() == 1);        // single element 
		assert(nodelist.Length() == 1);           // single load

		assert(elementlist(1) <= GetMBS()->NE()); // element exists
		Element* el = (Element*) GetMBS()->GetElementPtr(elementlist(1));
	
		assert(nodelist(1) <= el->NLoads());      // load is defined on element
		MBSLoad& ld = el->GetLoad(nodelist(1));

		actvalue = ld.Evaluate(time);			// load factor

		//$ DR 2012-10: the MBSSensor for Tpointload is currently not available. There will be a new Sensor soon!"; 
		//$ DR 2012-03-02: extension for pointloads
		//if((!el->IsRigid())&&(ld.GetLoadType() == Tpointload))
		//{
			//Vector3D force,pos;
			//ld.GetForceVector3D(force, pos);
			//force *= actvalue;
			//if (type & TSX)				{ actvalue = force.X();}
			//else if(type & TSY)		{ actvalue = force.Y();}
			//else if(type & TSZ)		{ actvalue = force.Z();}
		//}
	}

	//typedef enum {TSElement = 1, TSNode = 2, TSPos = 4, TSVel = 8, TSAngle = 16,
	//	TSX = 32, TSY = 64, TSZ = 128, TSDist = 256, TSDeflection = 512, TSDOF = 1024, TSplanar = 2048

	actvalue = factor * actvalue + offset;
	
	if(HasReferenceValues())	
	{
		double val= referenceValues->Evaluate(time);
		actvalue -= val;
	}
	if(GetReferenceSensorNumber()) //$ RL 2011-02:[ 
	{
		assert(referenceSensorNumber < GetMBS()->NSensors());
		//$ YV 2012-06: the sensors may produce just one scalar value
		actvalue -= GetMBS()->GetSensor(referenceSensorNumber).GetCurrentValueWithSensorProcessing(GetMBS()->GetTime());
	}
	
	//$ RL 2011-02:]
	
	return actvalue;

}





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int MBSSensor::GetNumberOfDrawingPositions()
{
	if(type & TSDeflection)
		return 3;
	if(type & TSDist)
		return 2;
	return 1;
}

Vector3D MBSSensor::GetDrawPosition(int i)
{
	if (type & TSDOF)
	{
		return GetMBS()->GetElement(elementlist(1)).GetRefPosD();
	}
	else if(type & TSAccel)
	{
		return GetMBS()->GetElement(elementlist(1)).GetPosD(poslist(1));
	}

	if (type & TSOutputSensor)
	{
		const InputOutputElement& ioe = (const InputOutputElement&)GetMBS()->GetElement(elementlist(1));
		return ioe.ToP3D(ioe.GetOutputPosD(nodelist(1)));
	}

	if(type & TSForce)		//$ DR 2011-12-20
	{
		return GetMBS()->GetElement(elementlist(1)).GetRefPosD();
	}

	if ((type & TSNode) && !(type & TSElement) && elementlist.Length() == 0) //FE node sensor
	{
		return GetMBS()->GetNode(nodelist(1)).GetPosD();
	}
	else //other sensors
		if (type & (TSPos | TSVel | TSDeflection | TSDist | TSStress | TSStrain | TSInelStrain | TSAxialExtension)) //TSStress | TSStrain | TSInelStrain not tested!
		{
			if (!(type & TSplanar))		// 3D
			{
				Body3D* e;
				if (!(type & TSAuxElem) )
					e = (Body3D*)GetMBS()->GetElementPtr(elementlist(i));
				else
					e = (Body3D*)GetMBS()->GetAuxElementPtr(elementlist(i));

				if (type & TSElement)
					return e->GetPosD(poslist(i));
				else if (type & TSNode)
					return e->GetNodePosD(nodelist(i));
			}
			else		// 2D
			{
				Body2D* e;
				if (!(type & TSAuxElem) )
					e = (Body2D*)GetMBS()->GetElementPtr(elementlist(i));
				else
					e = (Body2D*)GetMBS()->GetAuxElementPtr(elementlist(i));

				if (type & TSElement)
					return e->ToP3D(e->GetPos2DD(poslist2D(i)));
				else if (type & TSNode)
					return e->ToP3D(e->GetNodePos2DD(nodelist(i)));
			}
		}
	return 0;
}

mystr MBSSensor::GetTypeName()
{
// replaced "_" with "-" (for axis labels in matlab)
	mystr name = "";

	if (type&TSAuxElem) name = "Aux";

	if(type&TSEigenValue)
	{
		name = mystr("EigenFrequency")+mystr("-DOF")+mystr(nodelist(1));
	}
	else if (type&TSDOF) 
	{
		name += mystr("Element")+mystr(elementlist(1))+mystr("-DOF")+mystr(nodelist(1));
	}
	else if (type&TSOutputSensor) 
	{
		name += mystr("Element")+mystr(elementlist(1))+mystr("-Output")+mystr(nodelist(1));
	}
	else if (type&TSStress || type&TSStrain || type&TSInelStrain) 
	{
		if (type&TSStress) name += mystr("Element")+mystr(elementlist(1))+mystr("-Stress");
		if (type&TSStrain) name += mystr("Element")+mystr(elementlist(1))+mystr("-Strain");
		if (type&TSInelStrain) name += mystr("Element")+mystr(elementlist(1))+mystr("-InelStrain");

		if (tensor_comp == TSMises) name += mystr("-Mises");
		else if (tensor_comp == TSCompXX) name += mystr("-XX");
		else if (tensor_comp == TSCompYY) name += mystr("-YY");
		else if (tensor_comp == TSCompZZ) name += mystr("-ZZ");
		else if (tensor_comp == TSCompYZ) name += mystr("-YZ");
		else if (tensor_comp == TSCompXZ) name += mystr("-XZ");
		else if (tensor_comp == TSCompXY) name += mystr("-XY");
	}
	else if (type&TSAxialExtension)
	{
		name += mystr("Element")+mystr(elementlist(1))+mystr("-AxialExtension");
	}
	else if (type&TSForce) 
	{
		name += mystr("Element")+mystr(elementlist(1))+mystr("-Force");//+mystr(nodelist(1));
	}
	else if (type&TSSpecialSensorValue) 
	{
		name += mystr("Element")+mystr(elementlist(1))+mystr("-SpecialSensorValue-")+mystr(nodelist(1));
	}
	else if (type&TSAccel) 
	{
		if(type&TSElement)
		{
			name += mystr("Element")+mystr(elementlist(1))+mystr("-");
			if(type&TSLocalAxis)
			{
				name += mystr("Loc-");
			}
		}
		else if(type&TSDOF)
		{
			name += mystr("Element")+mystr(elementlist(1))+mystr("-AccelDOF")+mystr(nodelist(1));
		}
		else
		{
			assert(0); // type not known (use TSDOF or TSElement)
		}
	}
/*	else if (type&TSMultSensor) 
	{
		name += mystr("MultipleSensor");
	}*/
	else
	{
		if (!((type & TSNode) && !(type & TSElement) && elementlist.Length() == 0))
		{
			name += mystr("Element")+mystr(elementlist(1))+mystr("-");
		}
	}

	if (type&TSNode)
	{
		name += "Node" + mystr(nodelist(1)) + mystr("-");
	}

	if (type&TSDist) 
	{
		name += mystr(elementlist(2)) + mystr("-Dist");
	}
	if (type&TSDeflection) 
	{
		name += mystr(elementlist(2)) + mystr("-") + mystr(elementlist(3)) + mystr("-Deflection");
	}

	if (!((type&TSDist) || (type&TSDeflection) || (type&TSAxialExtension)))
	{
/* (AD): old
		if (type&TSPos) name += "p";
		if (type&TSDisplacement) name += "u";
		if (type&TSVel) name += "v";
		if (type&TSAngle) name += "angle";
		if (type&TSCurvature) name += "curvature";
		if (type&TSX) name += "x";
		if (type&TSY) name += "y";
		if (type&TSZ) name += "z";
		if (type&TSForce) name += "force";
*/
// (AD): new
		if (type&TSplanar) name += "2D-";
	
		if (type&TSX) name += "x-";
		if (type&TSY) name += "y-";
		if (type&TSZ) name += "z-";
		if (type&TSCurvature) name += "curvature";

		if ((type&TSAngle) && !(type&TSVel)) name += "angle";
		if ((type&TSAngle) && (type&TSVel)) name += "angular velocity";

		if (type&TSPos) name += "position";
		if (type&TSDisplacement) name += "displacement";
		if ((type&TSVel) && !(type&TSAngle)) name += "velocity";
		if (type&TSForce) name += "force";
		if ((type&TSAccel) && (type&TSElement)) name += "acceleration";//$ RL 2011-3-15: 
	}

	return name;
}

//$ RL 2011-7-25:[ 

// reference data vector from MathFunction are evaluated at given time-points "times" and stored in the vector "refvalues"
void MBSSensor::GetReferenceValues(const Vector& times, Vector& refvalues) const
{		
	refvalues.SetLen(times.Length());
	if(hasReferenceValues)
	{
		for(int i=1;i<=times.Length();i++)
		{			
			refvalues(i) = referenceValues->Evaluate(times.Get(i));
		}
	}
	else
	{
		assert(0); // should not happen!!!
		refvalues.SetAll(0.);
	}	
}
//$ RL 2011-7-25:]
//$ RL 2011-8-18:[ 
//$ RL 2012-1-12: moved to linalg.h, use GetExtremeValues instead
//Vector MBSSensor::GetNExtremeValuesFFT(const Vector& vec, TArray<int>& indices, int minmax) const
//{
//	Vector v(fft_NLocalExtremeValues);
//	indices.SetLen(0);
//
//	double vecoldold = vec(1)*minmax;
//	double vecold = vec(2)*minmax;
//
//	for(int i=3;i<=vec.Length();i++)
//	{
//		if(vecoldold < vecold && vecold > vec(i)*minmax)
//		{
//			// vec(i-1)*factor is already a maximum 
//			indices.Add(i-1);
//			v(indices.Length()) = vecold;
//		}
//		if(indices.Length() == v.Length())
//		{
//			break; // N extreme values found
//		}
//		vecoldold = vecold;
//		vecold = vec(i)*minmax;
//	}
//	return v;
//}

//$ RL 2011-7-14: FFT code
// get one global (minimum) maximum or N local (minima ) maxima with respect to flag 'minmax' = 1 (-1)
Vector MBSSensor::GetExtremeValuesFFT(const Vector& vec, TArray<int>& indices, int minmax) const
{
	if(GetFlag_SensorComputationFFT() & TSCFFTlocalExtremeVals)
	{
		//evaluate N extrema
		TArray<double> extr;
		TArray<int> ind;

		TArray<int> dummy1;
		TArray<double> dummy2;
		if(minmax > 0)
		{
			// max
			GetExtremeValues(vec,ind,extr,dummy1,dummy2,0);  //old:  return GetNExtremeValuesFFT(vec, indices, minmax);
		}
		else
		{
			// min(v) = max(-v)
			Vector vec2 = -1.0*vec;
			GetExtremeValues(vec2,ind,extr,dummy1,dummy2,0); //old: return GetNExtremeValuesFFT(vec, indices, minmax);
		}
		if(fft_NLocalExtremeValues>extr.Length())
		{
			//if fft_NLocalExtremeValues is set to big value, all min/max values are used.
			indices = ind;
			return extr;
		}
		else
		{
			indices.CopyFrom(ind, fft_NLocalExtremeValues);
			TArray<double> tmp;
			tmp.CopyFrom(extr, fft_NLocalExtremeValues);
			return Vector(tmp);
		}
	}
	else
	{
		// evaluate 1 global extremum
		Vector v(1);
		v(1) = vec(1);
		indices(1) = 1;
		for(int i=2; i<=vec.Length();i++)
		{
			double vold=v(1);
			v(1) = Maximum(minmax*v(1), vec(i));
			if(vold!=v(1))
			{
				indices(1) = i; // new index
			}
		}
		return v;
	}
}
//$ RL 2011-8-18:] 

bool MBSSensor::IsConsistent(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	//check element numbers
	int rv = 0;
	for (int i=1; i <= GetNumberOfRelatedElements(); i++)
	{
		if (type&TSElement && (GetRelatedElementNumber(i) <= 0 || GetRelatedElementNumber(i) > GetMBS()->NE()))
		{
			errorstr += "    Sensor has invalid element numbers!\n";
			rv = 1;
		}
		else if (type&TSAuxElem && (GetRelatedElementNumber(i) <= 0 || GetRelatedElementNumber(i) > GetMBS()->NAuxE()))
		{
			errorstr += "    Sensor has invalid element numbers!\n";
			rv = 1;
		}
	}

	if (!rv)
	{
		if ((type&TSDOF)) 
		{
			if (nodelist(1) <= 0 || nodelist(1) > ( GetMBS()->GetElement(GetRelatedElementNumber(1)).SOS()*2 + GetMBS()->GetElement(GetRelatedElementNumber(1)).ES() + GetMBS()->GetElement(GetRelatedElementNumber(1)).IS() ))
			{
				errorstr += "    DOF Sensor has invalid DOF reference numbers!\n";
				rv = 1;
			}
		}
		if (type&TSAccel)//$ RL 2011-3-15:[ TSAccelDOF --> TSAccel
		{
			if ((type&TSDOF))
			{
				if (nodelist(1) <= 0 || nodelist(1) > GetMBS()->GetElement(GetRelatedElementNumber(1)).SS())
				{
					errorstr += "    DOF Sensor has invalid DOF reference numbers!\n";
					rv = 1;
				}
			}
			else if(!type&TSElement)
			{
				errorstr += "    Acceleration Sensor has invalid type!\n";//$ RL 2011-3-15: use TSDOF or TSElement!
				rv = 1;
			}
		}//$ RL 2011-3-15:] TSAccelDOF --> TSAccel
		if (type&TSElement) //sensor elements must be bodies with position ...
		{
			for (int i=1; i <= GetNumberOfRelatedElements(); i++)
			{
				Element* e = GetMBS()->GetElementPtr(GetRelatedElementNumber(i));

				if ((((e->IsType(TConstraint) && !(type&TSForce || type&TSDOF))) && !(e->IsType(TController))) && !(type&TSSpecialSensorValue))	// DR 2012-01-12 TSSpecialSensorValue added
				{
					errorstr += mystr("    Invalid sensor element ")+mystr(GetRelatedElementNumber(i))+mystr("! Only body position/velocity/angle can be measured!\n");
					rv = 1;
				}
				if (!(type&TSDOF) && !(type&TSOutputSensor) && !(type&TSXData)) 
				//$ RL 2011-3-15: 'type&' added. old://if (!TSDOF && !TSAccelDOF && !TSOutputSensor)//$ RL 2011-3-15: commented out -> makes no sense!
				//$ AH 27-12-2011: added TSXData
				{ //check if 2D/3D consistent
					if (((type&TSplanar) && e->Dim() != 2) || (!(type&TSplanar) && e->Dim() != 3))
					{
						errorstr += mystr("    Invalid sensor element ")+mystr(GetRelatedElementNumber(i))+mystr("! Body and sensor dimension are different (2D <==> 3D)!\n");
						rv = 1;
					}
				}
			}
		}
		else if (type&TSAuxElem) //sensor elements must be bodies with position ...
		{
			for (int i=1; i <= GetNumberOfRelatedElements(); i++)
			{
				Element* e = GetMBS()->GetAuxElementPtr(GetRelatedElementNumber(i));

				if (((e->IsType(TConstraint) && !(type&TSForce || type&TSDOF))) && !(e->IsType(TController)))
				{
					errorstr += mystr("    Invalid sensor aux element ")+mystr(GetRelatedElementNumber(i))+mystr("! Only body position/velocity/angle can be measured!\n");
					rv = 1;
				}
				if (!(type&TSDOF) && !(type&TSOutputSensor)) //$ RL 2011-3-15: 'type&' added.
				{ //check if 2D/3D consistent
					if (((type&TSplanar) && e->Dim() != 2) || (!(type&TSplanar) && e->Dim() != 3))
					{
						errorstr += mystr("    Invalid sensor aux element ")+mystr(GetRelatedElementNumber(i))+mystr("! Body and sensor dimension are different (2D <==> 3D)!\n");
						rv = 1;
					}
				}
			}
		}
	}

	//$ YV 2012-06: the sense of this flag has changed
	return !rv;
}


void MBSSensor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	mystr sname = "";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//general Sensor data:
	ed.SetInt(0,"Sensor"); edc.Add(ed); //dummy field

//	ed.SetInt(GetSensorNumber(), "Sensor_number"); ed.SetLocked(1); edc.Add(ed);		//$ YV 2012-06: sensor number has been removed from edc, as it is not used in SetElementData anyway
	ed.SetText(GetSensorName(), "Sensor_name"); edc.Add(ed);
	
	//ed.SetText(GetDisplayName(), "Display_name"); ed.SetToolTipText("Used for PlotTool caption"); edc.Add(ed);
	/*
	ed.SetBool(flag_initialize_plottool_at_assemble, "Initialize_PlotTool"); ed.SetToolTipText("If checked, a PlotTool instance for the Sensor will be opened"); edc.Add(ed);

	mystr enstr = "Element_number";
	mystr str1 = "";

	ed.SetBool((type&TSAuxElem) > 0, "Use_aux_elements"); edc.Add(ed);

	if ((type&TSDist) || (type&TSDeflection)) str1 = "1";

	ed.SetBool(visible, "Visible"); ed.SetToolTipText("If unchecked, the sensor will never be drawn"); edc.Add(ed);
	ed.SetVector2D(draw_dim.X(), draw_dim.Y(), "Draw_parameters"); ed.SetToolTipText("[sphere size, sphere resolution]"); edc.Add(ed);

	//depreciated: ed.SetInt(writeresults, "Write_results"); ed.SetToolTipText("1=write results to general file, 2=write results to single file, 3=write both"); edc.Add(ed);
	int writegeneralfile = writeresults&1;
	int writeownsolfile = (writeresults&2) > 1;
	ed.SetBool(writegeneralfile, "Write_results_general"); ed.SetToolTipText("Write results to general solution file"); edc.Add(ed);
	ed.SetBool(writeownsolfile, "Write_results_own_file"); ed.SetToolTipText("Write results to separate solution file"); edc.Add(ed);

	ed.SetInt(precision, "Output_precision",1,17); ed.SetToolTipText("The number of digits used in the output file"); edc.Add(ed);
	ed.SetDouble(factor, "Output_factor"); ed.SetToolTipText("Output_sensor_value = offset + factor * sensor_value"); edc.Add(ed);
	ed.SetDouble(offset, "Output_offset"); ed.SetToolTipText("Output_sensor_value = offset + factor * sensor_value"); edc.Add(ed);


	if (elementlist.Length() != 0) ed.SetInt(elementlist(1), (enstr+str1).c_str(), 1, GetMBS()->NE()); edc.Add(ed);
	if (type&TSDist) 
	{ed.SetInt(elementlist(2), (enstr+mystr("2")).c_str()); edc.Add(ed);}
	if (type&TSDeflection) 
	{ed.SetInt(elementlist(2), (enstr+mystr("2")).c_str()); edc.Add(ed);
	ed.SetInt(elementlist(3), (enstr+mystr("3")).c_str()); edc.Add(ed);}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//specific Sensor data:

	if (type&TSDOF) 
	{
		sname = "Local_DOF_Sensor";
		if(type&TSAccel)sname="Accel_DOF_Sensor";
		ed.SetInt(nodelist(1), "Local_DOF", 1, MYMAXINT); edc.Add(ed);
	}
	//$ RL 2011-3-15:[ not used any more
	//else if ((type&TSAccel) && (type&TSDOF)) 
	//{
	//	sname = "Local_DOF_Acceleration_Sensor";
	//	ed.SetInt(nodelist(1), "Local_DOF", 1, MYMAXINT); edc.Add(ed);
	//}
	//$ RL 2011-3-15:] not used any more

	if (type&TSNode)
	{
		mystr str = "Local_node_number";
		ed.SetInt(nodelist(1), (str+str1).c_str()); edc.Add(ed);
	}
	//else if (!((type&TSDOF) || (type&TSAccelDOF)))$ RL 2011-3-15: TSAccelDOF not used any more
	else if(!(type&TSDOF))
	{
		if (!(type&TSplanar))
		{
			SetElemDataVector3D(edc, poslist(1), (mystr("Local_pos")+str1)); edc.Get(edc.Length()).SetToolTipText("point1 [X, Y, Z]");
		}
		else
		{
			SetElemDataVector2D(edc, poslist2D(1), (mystr("Local_pos2D")+str1)); edc.Get(edc.Length()).SetToolTipText("point1 [X, Y]");
		}
	}

	//Outputsensor
	if (type&TSOutputSensor) 
	{
		sname = "IOElement_Sensor";
		GetMBS()->UO() << "ERROR: IOElement_Sensor not implemented for saving!\n";
	}
	//Stress/strain sensor
	else if (type&TSStress || type&TSStrain || type&TSInelStrain) 
	{
		if (type&TSStress) 
		{
			sname = "Stress_Sensor";
		}
		if (type&TSStrain) 
		{
			sname = "Strain_Sensor";
		}
		if (type&TSInelStrain) 
		{
			sname = "InelStrain_Sensor";
		}
		int c1 = 0;
		int c2 = 0;
		switch(tensor_comp)
		{
		case 0: c1=0; c2=0; break;
		case 1: c1=1; c2=1; break;
		case 2: c1=2; c2=2; break;
		case 3: c1=3; c2=3; break;
		case 4: c1=2; c2=3; break;
		case 5: c1=1; c2=3; break;
		case 6: c1=1; c2=2; break;
		default: ;
		}
		Vector2D comps(c1,c2);
		SetElemDataVector2D(edc, comps, "Tensor_Components");
	}
	else if (type&TSAxialExtension) 
	{
		sname = "AxialExtension_Sensor";
		GetMBS()->UO() << "ERROR: Force_Sensor not implemented for saving!\n";
	}
	else if (type&TSForce) 
	{
		sname = "Force_Sensor";
		GetMBS()->UO() << "ERROR: Force_Sensor not implemented for saving!\n";
	}
	else if (type&TSSpecialSensorValue) 
	{
		sname = "Special-Value_Sensor";
		GetMBS()->UO() << "ERROR: Special-Value_Sensor not implemented for saving!\n";
	}
	//position/velocity/acceleration: +++++++++++++++++++++++++++++++++++++++++++++++
	else if ((type&TSPos || type&TSVel || type&TSAccel) && !(type&TSAngle))
	{
		if(type&TSPos)sname = "Position_Sensor";
		if(type&TSVel)sname = "Velocity_Sensor";
		if(type&TSAccel)sname = "Acceleration_Sensor";

		ed.SetBoolGroup((type&TSPos) > 0, 1, "Position"); edc.Add(ed); //group 1
		ed.SetBoolGroup((type&TSVel) > 0, 1, "Velocity"); edc.Add(ed); //group 1
		ed.SetBoolGroup((type&TSAccel) > 0, 1, "Acceleration"); edc.Add(ed); //group 1

		int max = 3;
		if (type&TSplanar) max = 2;

		ed.SetBoolGroup((type&TSX) > 0, 2, "x-component"); edc.Add(ed); //group 2
		ed.SetBoolGroup((type&TSY) > 0, 2, "y-component"); edc.Add(ed); //group 2
		if (max==3) 
		{ 
			ed.SetBoolGroup((type&TSZ) > 0, 2, "z-component"); edc.Add(ed);  //group 2
		}
	}
	//angle/angular velocity: +++++++++++++++++++++++++++++++++++++++++++++++
	else if (type&TSAngle)
	{
		sname = "Angle_Sensor";

		ed.SetBoolGroup((type&TSVel) == 0, 1, "Angle_of_rotation"); edc.Add(ed); //group 1
		ed.SetBoolGroup((type&TSVel) > 0, 1, "Angular_velocity"); edc.Add(ed); //group 1

		if (!(type&TSplanar))
		{
			SetElemDataVector3D(edc, poslist(2), "Rotation_axis"); edc.Get(edc.Length()).SetToolTipText("axis of rotation [X, Y, Z]");
			ed.SetBool((type&TSLocalAxis) > 0, "Body_fixed_axis"); ed.SetToolTipText("Rotation axis is given in local body-fixed coordinates"); edc.Add(ed);
		}
	}
	//distance between two points: +++++++++++++++++++++++++++++++++++++++++++++++
	else if (type&TSDist) 
	{
		if (type&TSplanar) GetMBS()->UO() << "ERROR: Distance sensor ElementData not implemented for 2D!!!\n";

		sname = "Distance_Sensor";

		if (type&TSNode)
		{
			ed.SetInt(nodelist(2), "Local_node_number2"); edc.Add(ed);
		}

		if (type&TSNode)
		{
			ed.SetInt(nodelist(2), "Local_node_number2"); edc.Add(ed);
		}
		else
		{
			SetElemDataVector3D(edc, poslist(2), "Local_pos2"); edc.Get(edc.Length()).SetToolTipText("point2 [X, Y, Z]");
		}
	}
	else if (type&TSDeflection) 
	{
		if (type&TSplanar) GetMBS()->UO() << "ERROR: Deflection sensor ElementData not implemented for 2D!!!\n";

		sname = "Deflection_Sensor";

		if (type&TSNode)
		{
			ed.SetInt(nodelist(2), "Local_node_number2"); edc.Add(ed);
			ed.SetInt(nodelist(3), "Local_node_number3"); edc.Add(ed);
			SetElemDataVector3D(edc, poslist(1), "Normal_vector"); edc.Get(edc.Length()).SetToolTipText("normal to the deflection that defines the sign of the deflection [X, Y, Z]"); // D.R. added this line, 13.01.2011 
		}
		else
		{
			SetElemDataVector3D(edc, poslist(2), "Local_pos2"); edc.Get(edc.Length()).SetToolTipText("point2 [X, Y, Z]");
			SetElemDataVector3D(edc, poslist(3), "Local_pos3"); edc.Get(edc.Length()).SetToolTipText("midpoint [X, Y, Z]");
			SetElemDataVector3D(edc, poslist(4), "Normal_vector"); edc.Get(edc.Length()).SetToolTipText("normal to the deflection that defines the sign of the deflection [X, Y, Z]");
		}
	}
	ed.SetText(sname, "Sensor_type"); ed.SetLocked(1); edc.Get(1) = ed;
	*/
}

int MBSSensor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer //$ DR 2012-10 return value changed from bool to int
{
	int rv = 1;

	/*

	mystr oldsensorname = GetSensorName();
	mystr oldtypename = GetTypeName();

	//initial data:
	writeresults = 1;
	precision = 17;
	factor = 1;
	offset = 0;
	draw_dim = 0;
	visible = 0;
	sensorname = "";


	TMBSSensor oldtype = type;
	type = (TMBSSensor)0;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//general Sensor data:

	//mystr ttype;
	//GetElemDataText(mbs, edc, "Sensor_type", ttype, 1); //is not stored!

	int auxflag, flag;
	GetElemDataBool(mbs, edc, "Use_aux_elements", auxflag, 0);
	if (auxflag) type = (TMBSSensor)(type+TSAuxElem); 
	if (oldtype&TSplanar) type = (TMBSSensor)(type+TSplanar); //not changeable!
	int nodebased = 0;

	GetElemDataBool(mbs, edc, "Initialize_PlotTool", flag_initialize_plottool_at_assemble, 0);
	GetElemDataBool(mbs, edc, "Visible", visible,0);

	Vector2D vec2;
	GetElemDataVector2D(mbs, edc, "Draw_parameters", vec2,0); draw_dim.X() = vec2(1); draw_dim.Y() = vec2(2);

	int res1 = 0; 
	int res2 = 0;
	GetElemDataBool(mbs, edc, "Write_results_general", res1, 0);
	GetElemDataBool(mbs, edc, "Write_results_own_file", res2, 0);
	writeresults = res1 + 2*res2;

	GetElemDataInt(mbs, edc, "Output_precision", precision,0);
	GetElemDataDouble(mbs, edc, "Output_factor", factor,0);
	GetElemDataDouble(mbs, edc, "Output_offset", offset,0);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//specific Sensor data:

	//if (oldtype&TSAccelDOF)//$ RL 2011-3-15:[ not needed any more (no code duplication!).
	//{
	//	type = (TMBSSensor)(oldtype); //not changeable!
	//	GetElemDataInt(mbs, edc, "Element_number", elementlist(1), 1);
	//	GetElemDataInt(mbs, edc, "Local_DOF", nodelist(1), 1);

	//	//check DOF:
	//	Element* e;
	//	if (auxflag) e = GetMBS()->GetAuxElementPtr(elementlist(1));
	//	else e = GetMBS()->GetElementPtr(elementlist(1));

	//	if (e)
	//	{
	//		if (nodelist(1) < 1 || nodelist(1) > e->SS())
	//		{
	//			nodelist(1) = 1;
	//			char str[256];
	//			sprintf(str, "The local DOF number %d is invalid for element number %d and needs to be corrected!", nodelist(1), elementlist(1));
	//			GetMBS()->UO().InstantMessageText(str);
	//			rv = 0;
	//		}
	//	}
	//}
	//else if (oldtype&TSDOF)//$ RL 2011-3-15:] not needed any more (no code duplication!).
	if (oldtype&TSDOF)
	{
		type = (TMBSSensor)(type+TSDOF);
		if(oldtype&TSAccel)
		{
			type = (TMBSSensor)(type+TSAccel);
		}
		GetElemDataInt(mbs, edc, "Element_number", elementlist(1), 1);
		GetElemDataInt(mbs, edc, "Local_DOF", nodelist(1), 1);

		//check DOF:
		Element* e;
		if (auxflag) e = GetMBS()->GetAuxElementPtr(elementlist(1));
		else e = GetMBS()->GetElementPtr(elementlist(1));

		if (e)
		{
			if (nodelist(1) < 1 || nodelist(1) > e->SS())
			{
				nodelist(1) = 1;
				char str[256];
				sprintf(str, "The local DOF number %d is invalid for element number %d and needs to be corrected!", nodelist(1), elementlist(1));
				GetMBS()->UO().InstantMessageText(str);
				rv = 0;
			}
		}
	}
	//position/velocity: ++++++++++++++++++++++++++++++++++++++++++++
	else if (((oldtype&TSPos) || (oldtype&TSVel) || (oldtype&TSAccel)) && !(oldtype&TSAngle))
	{
		GetElemDataBool(mbs, edc, "Position", flag);
		if (flag) type = (TMBSSensor)(type+TSPos);
		GetElemDataBool(mbs, edc, "Velocity", flag);
		if (flag) type = (TMBSSensor)(type+TSVel);
		GetElemDataBool(mbs, edc, "Acceleration", flag);
		if (flag) type = (TMBSSensor)(type+TSAccel);

		if (GetElemDataInt(mbs, edc, "Element_number", elementlist(1), 1))
		{
			if (GetMBS()->GetElement(elementlist(1)).Dim() == 2 && !(type&TSplanar)
				|| GetMBS()->GetElement(elementlist(1)).Dim() == 3 && (type&TSplanar))
				GetMBS()->EDCError("Planar sensors only work for planar elements, spatial sensors for spatial elements!");

		}

		if (type&TSplanar) GetMBS()->EDCError("Sensors for planar elements do not yet work properly!!!!!");


		int elementtype = 0;
		if (type&TSplanar)
		{
			Vector2D v;
			if (GetElemDataVector2D(mbs, edc, "Local_pos2D", v, 0))
			{
				poslist2D(1) = v;
				type = (TMBSSensor)(type+TSElement);
				elementtype = 1;
			}
		}
		else
		{
			Vector3D v;
			if (GetElemDataVector3D(mbs, edc, "Local_pos", v, 0))
			{
				poslist(1) = v;
				type = (TMBSSensor)(type+TSElement);
				elementtype = 1;
			}
		}

		if (!elementtype)
		{
			GetElemDataInt(mbs, edc, "Local_node_number", nodelist(1), 1);
			type = (TMBSSensor)(type+TSNode);
		}

		GetElemDataBool(mbs, edc, "x-component", flag);
		if (flag) type = (TMBSSensor)(type+TSX);
		GetElemDataBool(mbs, edc, "y-component", flag);
		if (flag) type = (TMBSSensor)(type+TSY);
		GetElemDataBool(mbs, edc, "z-component", flag);
		if (flag) type = (TMBSSensor)(type+TSZ);


	}
	//angle/angular velocity: ++++++++++++++++++++++++++++++++++++++++++++
	else if (oldtype&TSAngle)
	{
		flag = 0;
		GetElemDataBool(mbs, edc, "Angle_of_rotation", flag);
		if (flag) type = (TMBSSensor)(type+TSAngle);

		flag = 0;
		GetElemDataBool(mbs, edc, "Angular_velocity", flag);
		if (flag) type = (TMBSSensor)(type+TSAngle+TSVel);

		if (GetElemDataInt(mbs, edc, "Element_number", elementlist(1), 1))
		{
			if (!(type&TSplanar))
			{
				//3D--> need rotation axis!!!
				GetElemDataVector3D(mbs, edc, "Rotation_axis", poslist(2), 1);
				flag = 0;
				GetElemDataBool(mbs, edc, "Body_fixed_axis", flag, 0);
				if (flag) type = (TMBSSensor)(type+TSLocalAxis);
			}
		}

		if (type&TSplanar)
		{
			Vector2D v;
			if (GetElemDataVector2D(mbs, edc, "Local_pos2D", v, 0))
			{
				poslist2D(1) = v;
			}
		}
		else
		{
			Vector3D v;
			if (GetElemDataVector3D(mbs, edc, "Local_pos", v, 0))
			{
				poslist(1) = v;
			}
		}

		type = (TMBSSensor)(type|TSElement); //nodes not possible (up to now)
	}
	//Distance: ++++++++++++++++++++++++++++++++++++++++++++
	else if (oldtype&TSDist)
	{
		if (type&TSplanar) GetMBS()->UO() << "ERROR: Distance sensor ElementData not implemented for 2D!!!\n";
		type = (TMBSSensor)(type+TSDist);
		elementlist.SetLen(2);
		poslist.SetLen(2);
		GetElemDataInt(mbs, edc, "Element_number1", elementlist(1), 1);
		GetElemDataInt(mbs, edc, "Element_number2", elementlist(2), 1);

		if (GetElemDataVector3D(mbs, edc, "Local_pos1", poslist(1), 0))
		{
			type = (TMBSSensor)(type+TSElement);
			GetElemDataVector3D(mbs, edc, "Local_pos2", poslist(2), 1);
		}
		else
		{
			GetElemDataInt(mbs, edc, "Local_node_number1", nodelist(1), 1);
			GetElemDataInt(mbs, edc, "Local_node_number2", nodelist(2), 1);
			type = (TMBSSensor)(type+TSNode);
		}
	}
	//Deflection: ++++++++++++++++++++++++++++++++++++++++++++
	else if (oldtype&TSDeflection)
	{
		if (type&TSplanar) GetMBS()->UO() << "ERROR: Deflection sensor ElementData not implemented for 2D!!!\n";
		type = (TMBSSensor)(type+TSDeflection);
		elementlist.SetLen(3);
		GetElemDataInt(mbs, edc, "Element_number1", elementlist(1), 1);
		GetElemDataInt(mbs, edc, "Element_number2", elementlist(2), 1);
		GetElemDataInt(mbs, edc, "Element_number3", elementlist(3), 1);

		if (GetElemDataVector3D(mbs, edc, "Local_pos1", poslist(1), 0))
		{
			type = (TMBSSensor)(type+TSElement);
			GetElemDataVector3D(mbs, edc, "Local_pos2", poslist(2), 1);
			GetElemDataVector3D(mbs, edc, "Local_pos3", poslist(3), 1);
			GetElemDataVector3D(mbs, edc, "Normal_vector", poslist(4), 1);
		}
		else
		{
			GetElemDataInt(mbs, edc, "Local_node_number1", nodelist(1), 1);
			GetElemDataInt(mbs, edc, "Local_node_number2", nodelist(2), 1);
			GetElemDataInt(mbs, edc, "Local_node_number3", nodelist(3), 1);
			GetElemDataVector3D(mbs, edc, "Normal_vector", poslist(1), 1);
			type = (TMBSSensor)(type+TSNode);
		}
	}
	if (oldtype&TSOutputSensor) 
	{
		type = (TMBSSensor)(type+TSOutputSensor);
		GetMBS()->UO() << "ERROR: IOElement_Sensor not implemented for reading!\n";
	}
	//Stress/strain sensor
	else if (oldtype&TSStress || oldtype&TSStrain || oldtype&TSInelStrain) 
	{
		if (oldtype&TSStress) 
		{
			type = (TMBSSensor)(type+TSStress);
		}
		if (type&TSStrain) 
		{
			type = (TMBSSensor)(type+TSStrain);
		}
		if (type&TSInelStrain) 
		{
			type = (TMBSSensor)(type+TSInelStrain);
		}

		Vector2D comps(0.,0.);
		GetElemDataVector2D(GetMBS(),edc, "Tensor_Components", comps);

		if (comps.X()==0 && comps.Y()==0) {tensor_comp = (TMBSTensorcomponent)0; }
		if (comps.X()==1 && comps.Y()==1) {tensor_comp = TSCompXX; }
		if (comps.X()==2 && comps.Y()==2) {tensor_comp = TSCompYY; }
		if (comps.X()==3 && comps.Y()==3) {tensor_comp = TSCompZZ; }
		if (comps.X()==2 && comps.Y()==3) {tensor_comp = TSCompYZ; }
		if (comps.X()==1 && comps.Y()==3) {tensor_comp = TSCompXZ; }
		if (comps.X()==1 && comps.Y()==2) {tensor_comp = TSCompXY; }

	}
	else if (type&TSAxialExtension) 
	{
		type = (TMBSSensor)(type+TSAxialExtension);
		GetMBS()->UO() << "ERROR: AxialExtension not implemented for reading!\n";
	}
	else if (type&TSForce) 
	{
		type = (TMBSSensor)(type+TSForce);
		GetMBS()->UO() << "ERROR: Force_Sensor not implemented for reading!\n";
	}

	//write sensor name at the end:
	GetElemDataText(mbs, edc, "Sensor_name", sensorname, 0);
	GetElemDataText(mbs, edc, "Display_name", displayname, 0);
	*/

	/*
	GetMBS()->UO() << "oldsname=" << oldsensorname << "\n";
	GetMBS()->UO() << "newsname=" << sensorname << "\n";
	GetMBS()->UO() << "newname==oldname" << (int)(sensorname == oldsensorname) << "\n";
	GetMBS()->UO() << "oldname==typename" << (int)(oldtypename == oldsensorname) << "\n";
	GetMBS()->UO() << "typename=" << GetTypeName() << "\n";
	*/

/*
	if ((sensorname == oldsensorname) && (oldtypename == oldsensorname))
	{
		//set new sensorname according to new elements!!!
		sensorname = GetTypeName();
	}
	*/

	return rv;
}

Vector MBSSensor::GetSensorComputationValue() const
{
	Vector v(1); //one value, initialized with 0

	if(GetFlag_SensorComputationFFT() && !FFT_amplitudes().Length())
	{
		v(1) = 0.;
		GetMBS()->UO().InstantMessageText("Warning: Empty FFT data. Sensor cost function value is set to zero.");
		return v;
	}

	switch (FlagSensorComputation())	
	{
	case TSCmin: 
		{
			if(GetFlag_SensorComputationFFT())
			{
				// frequency domain
				// use amplitude or phase for evaluation of computation value		
				// evaluate 1 or N exreme values				
				TArray<int> ind;
				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					v = GetExtremeValuesFFT(FFT_amplitudes(), ind, -1); // -1.0...search minimum
				}
				else if(GetFlag_SensorComputationFFT() & TSCFFTphase)
				{
					v = GetExtremeValuesFFT(FFT_phase(), ind, -1); // -1.0...search minimum
				}
				// return freqencies instead of extreme value
				if(GetFlag_SensorComputationFFT() & TSCFFTfrequency)
				{
					v = Vector(ind.Length());
					for(int i=1;i<=ind.Length();i++)
					{
						v(i) = FFT_frequencies()(ind(i));	
					}
				}
			}
			else
			{
				// time domain
				v(1) = GetMinVal(); 
			}
			break;
		}
	case TSCmax: 
		{
			if(GetFlag_SensorComputationFFT()) 
			{
				// frequency domain
				// use amplitude or phase for evaluation of computation value
				// evaluate 1 or N exreme values
				TArray<int> ind;
				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					v = GetExtremeValuesFFT(FFT_amplitudes(), ind, 1); //1.0...search maximum
				}
				else if(GetFlag_SensorComputationFFT() & TSCFFTphase)
				{
					v = GetExtremeValuesFFT(FFT_phase(), ind, 1); //1.0...search maximum
				}


				// return freqencies instead of extreme value
				if(GetFlag_SensorComputationFFT() & TSCFFTfrequency)
				{
					v = Vector(ind.Length());
					for(int i=1;i<=ind.Length();i++)
					{
						v(i) = FFT_frequencies()(ind(i));	
					}
				}
			}
			else
			{
				// time domain
				v(1) = GetMaxVal(); 
			}
			break;
		}
	case TSCminmax: 
		{
			if(GetFlag_SensorComputationFFT())
			{
				// frequency domain
				// evaluate 2 or 2*N exreme values: 
				TArray<int> indmin, indmax, ind;
				Vector vmin, vmax;
				// use amplitude or phase for evaluation of computation value
				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					vmin = GetExtremeValuesFFT(FFT_amplitudes(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_amplitudes(), indmax, 1); // 1.0...search maximum
				}
				else if(GetFlag_SensorComputationFFT() & TSCFFTphase)
				{
					vmin = GetExtremeValuesFFT(FFT_phase(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_phase(), indmax, 1); // 1.0...search maximum
				}
				// v= <min1, max1>...global maxima or v = <min1, min2,...,minN, max1, max2,...,maxN>
				v = Vector(indmin.Length()+indmax.Length());
				ind.SetLen(0);
				// set minimum values
				for(int i=1;i<=indmin.Length();i++)
				{
					v(i) = vmin(i);
					ind.Add(indmin(i));
				}
				// set maximum values
				for(int i=1;i<=indmax.Length();i++)
				{
					v(i+vmin.Length()) = vmax(i);
					ind.Add(indmax(i));
				}

				// return freqencies instead of extreme values
				if(GetFlag_SensorComputationFFT() & TSCFFTfrequency)
				{
					v = Vector(ind.Length()); 
					for(int i=1;i<=ind.Length();i++)
					{
						v(i) = FFT_frequencies()(ind(i));	//v = <fmin1, fmin2,...,fminN, fmax1, fmax2,...,fmaxN>
					}
				}
			}
			else
			{
				// time domain
				v = Vector(2); v(1) = GetMinVal(); v(2) = GetMaxVal(); 
			}
			break;
		}
	case TSCamplitude: 
		{
			if(GetFlag_SensorComputationFFT())
			{
				// use amplitude or phase for evaluation of computation value
				// evaluate 2 or 2*N exreme values: 
				TArray<int> indmin, indmax, ind;
				Vector vmin, vmax;	

				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					vmin = GetExtremeValuesFFT(FFT_amplitudes(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_amplitudes(), indmax, 1); // 1.0...search maximum
				}
				else if(GetFlag_SensorComputationFFT() & TSCFFTphase)
				{
					vmin = GetExtremeValuesFFT(FFT_phase(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_phase(), indmax, 1); // 1.0...search maximum
				}
				 
				// v= <min1, max1>...global maxima or v = <min1, min2,...,minN, max1, max2,...,maxN>
				
				int imax = Minimum(indmax.Length(),indmin.Length());//if not all extreme values found, the TSCamplitude is computed only for the values where both, maximum and minimum, are found.
				
				assert(imax == fft_NLocalExtremeValues); //TODO: insert warning instead of "assert" if necessary (not all extreme values found)
				
				v = Vector(imax);  
				ind.SetLen(0);
				// set minimum values
				for(int i=1;i<=imax;i++)
				{
					v(i) = vmax(i)-vmin(i);
				}
				
				// return freqencies instead of extreme values
				if(GetFlag_SensorComputationFFT() & TSCFFTfrequency)
				{
					v = Vector(imax); 
					for(int i=1;i<=imax;i++)
					{
						v(i) = FFT_frequencies()(indmax(i))-FFT_frequencies()(indmin(i));	//v = <fmax1-fmin1, fmax2-fmin2,...,fmaxN-fminN>
					}
				}		
			}
			else
			{
				// time domain
				v(1) = GetMaxVal() - GetMinVal();
			}
			break;
		}
	case TSCmaxabs: 
		{
			if(GetFlag_SensorComputationFFT())
			{
				// frequency domain
				// use amplitude or phase for evaluation of computation value
				TArray<int> indmin, indmax, ind;
				Vector vmin, vmax;
				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					vmin = GetExtremeValuesFFT(FFT_amplitudes(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_amplitudes(), indmax, 1); // 1.0...search maximum
				}
				else if(GetFlag_SensorComputationFFT() & TSCFFTphase)
				{
					vmin = GetExtremeValuesFFT(FFT_phase(), indmin,-1); //-1.0...search minimum
					vmax = GetExtremeValuesFFT(FFT_phase(), indmax, 1); // 1.0...search maximum
				}

				v = Vector(1);  
				
				// set minimum values
						
				int cvindex = 0; // cvindex > 0 ==> computation value is from vmax, cvindex < 0 ==> computation value is from vmin
				double vold;
				for(int i=1;i<=vmin.Length();i++)
				{
					vold = v(1);
					v(1) = Maximum(fabs(v(1)), fabs(vmin(i)));
					if(v(1)>vold)
					{
						cvindex = -i;
					}
				}

				for(int i=1;i<=vmax.Length();i++)
				{
					vold = v(1);
					v(1) = Maximum(fabs(v(1)), fabs(vmax(i)));
					if(v(1)>vold)
					{
						cvindex = i;
					}
				}

				// return freqencies instead of extreme values
				if(GetFlag_SensorComputationFFT() & TSCFFTfrequency)
				{ 
					if(cvindex < 0)
					{
						v(1) = FFT_frequencies()(indmin(-cvindex));
					}
					else if(cvindex > 0)
					{
						v(1) = FFT_frequencies()(indmax(cvindex));
					}
					else
					{
						v(1) = 0.; // no minimum or maximum found
					}
				}
			}
			else
			{
				// time domain
				v(1) = Maximum(fabs(GetMaxVal()), fabs(GetMinVal())); 
			}
			break;
		}
	case TSCmeanL2:
		{
			//returns average error e
			//return value: e = sqrt[1/N *( e1^2+e2^2+e3^3+...eN^2)]
			//with ei= xi-xrefi at timestep i; xi is value and xrefi is reference value, N is total number of timesteps
			if(GetFlag_SensorComputationFFT())
			{
				assert(!(GetFlag_SensorComputationFFT() & TSCFFTfrequency) && !(GetFlag_SensorComputationFFT() & TSCFFTlocalExtremeVals)&& !(GetFlag_SensorComputationFFT() & TSCFFTphase)); // makes no sense for L2-norm
				// frequency domain				
				//TArray<double> frequencies, amplitudes, phase;
				//LoadSensorDataFFT(frequencies, amplitudes, phase);

				//frequencies(1) = 0
				//frequencies(2) = 1/(N*fft_constant_sample_time)
				if(GetFlag_SensorComputationFFT() & TSCFFTamplitude)
				{
					//analyze amplitude spectrum for evaluation of computation value
					// only test case: int useLog = 0; // 0...mean quadratic error (default), 1...mean logarithmic error
					double sum = 0;
					for(int i=2;i<=FFT_amplitudes().Length();i++)
					{
						if(GetStartFrequencySensorComputation() <= FFT_frequencies()(i) && FFT_frequencies()(i) <= GetEndFrequencySensorComputation())
						{
							// only values of frequencies f€[fmin, fmax] are considered in the computation value
							//if(useLog)
							//{
							//	sum += log10(fabs(FFT_amplitudes().Get(i)));
							//}
							//else
							//{
								sum += Sqr(FFT_amplitudes().Get(i));
							//}
							if(FFT_frequencies()(i) > GetEndFrequencySensorComputation())
							{
								break;
							}
							//mbs->UO() << "sum" << mystr(i) << "=" << mystr(sum) << "\t";
							//mbs->UO() << "amplitude" << mystr(i) << "=" << amplitudes.Get(i) << "\t";
							//mbs->UO() << "log10(amplitudes(" << mystr(i) << ")=" << log10(amplitudes.Get(i)) << "\n";
						}
					}
					int add_constant_component = 0; // (0)/1 ...(Don't) consider 0 Hz
					if(GetStartFrequencySensorComputation() <= 0 && 0 <= GetEndFrequencySensorComputation())
					{	
						add_constant_component = 1;	
					}

					//if(useLog)
					//{
					//	v(1) = (add_constant_component*log10(fabs(FFT_amplitudes()(1)))+2.0*sum)/(double)n_sensor_computations;
					//}
					//else
					//{					
						//v(1) = 0.5*sqrt(Sqr(add_constant_component*FFT_amplitudes()(1))+2.0*sum); //L2-norm
						v(1) = 0.5*sqrt(Sqr(2.*add_constant_component*FFT_amplitudes()(1))+2.0*sum); //L2-norm

						//Differenzphase Mittelwert des Phasenfehlers würde ev. sinn machen, sqrt(1/N*sum(sqr(phase))
					//}
				}
			}
			else
			{
				// time domain
				v(1) = sqrt(GetRefVal() / (double)n_sensor_computations); //L2norm
			}
			break;
		}
	default: ;
	}
	// weight of sensor computation values (standard: 1.0)
	for(int i=1;i<=v.Length();i++)
	{
		v(i) = v(i)*computationValueWeight.Get(i);
	}
	return v;
}

void MBSSensor::DoSensorComputations()
{
	int fmm = FlagSensorComputation();

	/*
	if(GetWriteResults() & 4)
	{
		StoreSensorTimeValueData();
	}
	*/

	if(IsSensorComputation())
	{

		if (FlagSensorComputation() >= 1 && FlagSensorComputation() <= 5)
		{							 
			// just two values are compared
			SetMaxVal(Maximum(GetMaxVal(), GetLastValue()));
			SetMinVal(Minimum(GetMinVal(), GetLastValue()));
		}
		else if (FlagSensorComputation() == 6) //L2 norm
		{
			double val = GetRefVal();
			SetRefVal(val + Sqr(GetLastValue()));
		}
		IncrementNSensorComputation();

	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ RL 2011-7-14:[ FFT-sensor files are created
void MBSSensor::PostComputationOperationsFFT()
{
	//++++++++++++++++++++++++++++++++++++++++
	// options for fast fourier transformation
	int fft_use_window = 0; // 0... no window
	// 1... hamming window: 0.54-0.46*cos(2*pi*n/N) //hamming 0<=n<N
	//++++++++++++++++++++++++++++++++++++++++

	//do evaluations of fast fourier transformation after computation

	if(IsSensorComputationFFT()) //Sensor stores the times and values in Tarrays now
	{
		/*
		TArray<double> times;           // store time points for fft-computation
		TArray<double> values;          // store time points for fft-computation
		LoadSensorData(times, values); // load data from file, possible difference of reference data and measurement data
		*/
		if(signalHistoryTimes.Length() < 2 || signalHistoryValues.Length() < 2)
		{
			GetMBS()->UO(UO_LVL_err).InstantMessageText("Error in PostComputationOperationsFFT: No recorded sensor data found for fft-computation.");
			return; // should not happen!!!
		}

		double sampleTime = GetConstantSampleTimeFFT();
		if(sampleTime<=0)
		{	
			GetMBS()->UO(UO_LVL_warn) << "Warning: Invalid sample Time for fast fourier transformation found \n Sample time 0.001 s is used since now.";
			sampleTime = 0.001;
			SetConstantSampleTimeFFT(sampleTime);
		}
		//sample sensor output --> get constant time step	
		Vector vtime(signalHistoryTimes);
		Vector samples(signalHistoryValues);

		constTimeStep(vtime, samples, sampleTime, GetStartTimeSensorComputation(), GetEndTimeSensorComputation()); // content of samples has constant time step
		if(samples.Length()==0 || vtime.Length() == 0)
		{
			GetMBS()->UO(UO_LVL_warn).InstantMessageText("Warning: FFT computation aborted due to empty time or sample vector.");
			return;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//b: get reference values, which are already substracted from the output of the sensor values "samples"
		int flag_addrefsamples=0; // set to value unequal zero to compute amplitude spectrum of reference and measurement data seperately
		flag_addrefsamples = HasReferenceValues() + GetReferenceSensorNumber(); //unequal zero if reference math function or reference sensor exists

		// define OptionsEDC-flag, if necessary!!!


		flag_addrefsamples = flag_addrefsamples*GetFlag_Diff_FFTs_MesRef(); //switch on/off


		Vector refsamples;
		refsamples.SetLen(0);
		if(flag_addrefsamples)
		{				
			if(GetReferenceSensorNumber())
			{
				// load data of reference sensor from file
				int refSensNum = GetReferenceSensorNumber();
				Sensor & refSensor = GetMBS()->GetSensor(refSensNum);

				// the reference sensor has the values in memory (the flag SSM_InternalArray is set in MBSSensor::SetReferenceSensor())
				Vector vtime(signalHistoryTimes);
				refsamples = Vector(signalHistoryValues);
				constTimeStep(vtime, refsamples, sampleTime,GetStartTimeSensorComputation(), GetEndTimeSensorComputation());       // content of samples has constant time step
			}
			else
			{
				// evaluate math function for reference values
				GetReferenceValues(vtime, refsamples); // refsamples has already constant time step				
			}		
			//e: get reference values, which are already substracted from the output of the sensor values "samples"
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//b: make seperate ffts for reference samples and measurement data (use this, if phase shifts should not be encounted!)
			if(refsamples.Length())
			{					
				samples += refsamples;  // measurement - reference = difference ==> measurement = difference + reference
			}
		}
		//e: make seperate ffts for reference samples and measurement data
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//b: compute fft

		//-----------windowing-----------
		int oo = 0; // print fft-results in output window
		if(fft_use_window)
		{
			int N = samples.Length();
			for(int n=0;n<N;n++)
			{
				samples(n+1)*=0.54-0.46*cos(2*MY_PI*n/N);
			}
		}
		//----------windowing------------
		//Vector frequencies, amplitudes, phase;
		makefft(sampleTime, samples, FFT_frequencies(), FFT_amplitudes(), FFT_phase()); // samples ... "measurement" data(flag_addrefsamples=1 or no reference data available) or difference data (flag_addrefsamples=0 and reference data available)

		if(oo)
		{
			GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");
			GetMBS()->UO(UO_LVL_0) << mystr("%FFT of sensor values\n");
			GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");
			for(int j=1; j<=FFT_frequencies().Length();j++)
			{
				GetMBS()->UO(UO_LVL_0) << FFT_frequencies().Get(j) << mystr(" ") << FFT_amplitudes().Get(j) << mystr(" ") << FFT_phase().Get(j) << mystr("\n");
			}
		}

		if(flag_addrefsamples)
		{
			//-----------windowing-----------
			if(fft_use_window)
			{
				int N = refsamples.Length();
				for(int n=0;n<N;n++)
				{
					refsamples(n+1)*=0.54-0.46*cos(2*MY_PI*n/N);
				}
			}
			//----------windowing------------
			Vector amp2, phase2;
			makefft(sampleTime, refsamples, FFT_frequencies(), amp2, phase2); // refsamples ... reference data

			if(oo)
			{
				GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");	
				GetMBS()->UO(UO_LVL_0) << mystr("%FFT of reference values\n");
				GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");	
				for(int j=1; j<=FFT_frequencies().Length();j++)
				{
					GetMBS()->UO(UO_LVL_0) << FFT_frequencies().Get(j) << mystr(" ") << amp2(j) << mystr(" ") << phase2(j) << mystr("\n");
				}
			}

			FFT_amplitudes() -= amp2;  // difference of amplitude spectrums (measurement - reference)
			FFT_phase() -= phase2;     // difference of phase spectrums (measurement - reference)


			if(oo)
			{
				GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");	
				GetMBS()->UO(UO_LVL_0) << mystr("%FFT of sensor values minus FFT of reference values\n");
				GetMBS()->UO(UO_LVL_0) << mystr("%-------------------------------------------------------------------------\n");	
				for(int j=1; j<=FFT_frequencies().Length();j++)
				{
					GetMBS()->UO(UO_LVL_0) << FFT_frequencies().Get(j) << mystr(" ") << fabs(FFT_amplitudes().Get(j)) << mystr(" ") << FFT_phase().Get(j) << mystr("\n");
				}
			}
		}

		// make sure to have only positive amplitudes
		for(int j=1; j<=FFT_amplitudes().Length();j++)
		{
			FFT_amplitudes()(j) = fabs(FFT_amplitudes()(j));
		}

		// do smoothing here
		Vector tmp(FFT_amplitudes().Length());
		SmoothenDoubleArray(FFT_amplitudes(), tmp, FFToptimization_averaging);
		FFT_amplitudes()=Vector(tmp);
		SmoothenDoubleArray(FFT_phase(), tmp, FFToptimization_averaging);
		FFT_phase()=Vector(tmp);

		// write into file
		if((GetSignalStorageMode() & SSM_OwnFile) && GetOutputFileFFT())
		{
			//store fft data
			ofstream* ssol = GetOutputFileFFT();
			ssol->precision(GetPrecision());

			for (int j = 1; j <= FFT_frequencies().Length(); j++)
			{
				(*ssol) << FFT_frequencies().Get(j) << " " << FFT_amplitudes().Get(j)	<< " " <<  FFT_phase().Get(j) << "\n" << flush;
			}			
			(*ssol) << flush;
		}
	}

}
//$ RL 2011-7-14:] FFT-sensor files are created
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//$ YV 2012-06: MBSSensor::LoadSensorData() was used only in FFT for the actual sensor and for the reference one. But both sensors are explicitly switched into SSM_InternalArray mode. Therefore we don’t need the function and can take the values from the internal arrays.
/*
int MBSSensor::LoadSensorData(TArray<double>& timesi, TArray<double>& valuesi)
{
	//+++++++++++++++++++++++++++++++++++++++
	// get data from stored values
	if(times.Length() > 0 && values.Length() == times.Length())
	{
		timesi.CopyFrom(times);
		valuesi.CopyFrom(values);
		return 0; // ok
	}

	//+++++++++++++++++++++++++++++++++++++++
	// get values from sensor files
	timesi.SetLen(0);
	valuesi.SetLen(0);

	mystr dir = GetMBS()->GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path");
  mystr pathname = (dir+mystr("S")+mystr(GetSensorNumber())+mystr("-")+GetSensorName()+mystr(".txt"));
	CMatrixFile cf(pathname, (TFileMode)TFMread);
	if(!cf.IsGood()) 
	{
		GetMBS()->UO(UO_LVL_err).InstantMessageText("Warning: Problems during load sensor data from file!");
		return 1; // error
	}	
	cf.ReadAllColumns(); //only 2 columns in file
	//cf.ReadSpecifiedColumns(IntVec2(1,2)); //alternative: read column 1 and 2
	timesi.CopyFrom(cf.Column(1));
	valuesi.CopyFrom(cf.Column(2));
	return 0; // ok
}
*/

//// load sensor data from fft - file
//// frequencies, amplitudes, phase
//void MBSSensor::LoadSensorDataFFT(Vector& frequencies, Vector& amplitudes, Vector& phase) const
//{
//	//+++++++++++++++++++++++++++++++++++++++
//	// get data from stored values
//	if(fft_frequencies.Length() > 0 && fft_amplitudes.Length() == fft_phase.Length() && fft_frequencies.Length() == fft_phase.Length())
//	{
//		frequencies = fft_frequencies;
//		amplitudes = fft_amplitudes;
//		phase = fft_phase;
//		return;
//	}
//	
//	//+++++++++++++++++++++++++++++++++++++++
//	// get values from sensor-fft files
//	frequencies.SetLen(0);
//	amplitudes.SetLen(0);
//	phase.SetLen(0);
//
//	mystr dir = mbs->GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path");
//  mystr pathname = (dir+mystr("S")+mystr(GetSensorNumber())+mystr("-")+GetSensorName()+mystr("-fft.txt"));
//	CMatrixFile cf(pathname, (TFileMode)TFMread);
//	if(!cf.IsGood()) 
//	{
//		mbs->UO(UO_LVL_err).InstantMessageText("ERROR: could not load sensor-fft data from file!");
//		return;
//	}
//	cf.ReadAllColumns(); 
//	//cf.ReadSpecifiedColumns(IntVec2(1,2));
//	frequencies.CopyFrom(cf.Column(1));
//	amplitudes.CopyFrom(cf.Column(2));
//	phase.CopyFrom(cf.Column(3));
//}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MBSMultipleSensor::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	/*
	ElementData ed;
	mystr sname = "";

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//general Sensor data:
	ed.SetInt(0,"Sensor"); edc.Add(ed); //dummy field, lateron filled with name

	ed.SetInt(GetSensorNumber(), "Sensor_number"); ed.SetLocked(1); edc.Add(ed);
	ed.SetText(GetSensorName(), "Sensor_name"); edc.Add(ed);

	//change: make a list, add a list of int/double editable in edit field, comma-separated
	SetElemDataIVector(edc, sensors, "Sub_sensors"); 
	edc.Last().SetToolTipText("Comma-separated list of sensor numbers for the multiple sensor [sens1, sens2, ...]");
	edc.Last().SetVariableLength();
	//ed.SetInt(sensors(1), "Sub_Sensor1", 1, GetMBS()->NSensors()); edc.Add(ed);
	//ed.SetInt(sensors(2), "Sub_Sensor2", 1, GetMBS()->NSensors()); edc.Add(ed);

	int writegeneralfile = writeresults&1;
	int writeownsolfile = (writeresults&2) > 1;
	ed.SetBool(writegeneralfile, "Write_results_general"); ed.SetToolTipText("Write results to general solution file"); edc.Add(ed);
	ed.SetBool(writeownsolfile, "Write_results_own_file"); ed.SetToolTipText("Write results to separate solution file"); edc.Add(ed);

	ed.SetInt(precision, "Output_precision",1,17); ed.SetToolTipText("The number of digits used in the output file"); edc.Add(ed);
	ed.SetDouble(factor, "Output_factor"); ed.SetToolTipText("Output_sensor_value = offset + factor * sensor_value"); edc.Add(ed);
	ed.SetDouble(offset, "Output_offset"); ed.SetToolTipText("Output_sensor_value = offset + factor * sensor_value"); edc.Add(ed);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//specific Sensor data:

	ed.SetBoolGroup((multtype&TSMAverage) > 0, 1, "Average"); ed.SetToolTipText("result = 1 / nsensors *(sensor1 + sensor2 + sensor3 + ...)"); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMSum) > 0, 1, "Sum"); ed.SetToolTipText("result = sensor1 + sensor2 + sensor3 + ..."); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMDifference) > 0, 1, "Difference"); ed.SetToolTipText("result = sensor1 - (sensor2 + sensor3 + ...)"); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMMult) > 0, 1, "Mult"); ed.SetToolTipText("result = sensor1 * sensor2 * sensor3 * ..."); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMDiv) > 0, 1, "Div"); ed.SetToolTipText("result = sensor1/(sensor2 * sensor3 * ...)"); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMNorm) > 0, 1, "Norm"); ed.SetToolTipText("result = sqrt(sensor1^2 + sensor2^2 + sensor3^2 + ...)"); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMMax) > 0, 1, "Maximum"); ed.SetToolTipText("result = max(sensor1, sensor2, sensor3, ...)"); edc.Add(ed);
	ed.SetBoolGroup((multtype&TSMMin) > 0, 1, "Minimum"); ed.SetToolTipText("result = min(sensor1, sensor2, sensor3, ...)"); edc.Add(ed);

	ed.SetText("Multiple_Sensor", "Sensor_type"); ed.SetLocked(1); edc.Get(1) = ed;
	*/
}

int MBSMultipleSensor::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer //$ DR 2012-10 return value changed from bool to int
{
	int rv = 1;

	/*
	//initial data:
	writeresults = 1;
	precision = 17;
	factor = 1;
	offset = 0;
	draw_dim = 0;
	visible = 0;
	sensorname = "";

	//TMBSSensor oldtype = type;
	type = TSMultSensor;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//general Sensor data:
	GetElemDataText(mbs, edc, "Sensor_name", sensorname, 0);


	int flag;

	int res1 = 0; 
	int res2 = 0;
	GetElemDataBool(mbs, edc, "Write_results_general", res1, 0);
	GetElemDataBool(mbs, edc, "Write_results_own_file", res2, 0);
	writeresults = res1 + 2*res2;

	GetElemDataInt(mbs, edc, "Output_precision", precision,0);
	GetElemDataDouble(mbs, edc, "Output_factor", factor,0); 
	GetElemDataDouble(mbs, edc, "Output_offset", offset,0);

	//change: make a list, add a list of int/double editable in edit field, comma-separated
	GetElemDataIVector(GetMBS(), edc, "Sub_sensors", sensors, 1); 
	if (sensors.Length() == 0) 
	{
		GetMBS()->EDCError("The number of Sub_sensors is zero");
		rv = 0;
	}
	if (GetMBS()->NSensors() == 1)
	{
		GetMBS()->EDCError("There are no valid Sub_sensors available in the MBS system, the Sub_sensors are cleared");
		sensors.SetLen(0);
		rv = 0;
	}

	//GetElemDataInt(mbs, edc, "Sub_Sensor1", sensors(1), 1);
	//GetElemDataInt(mbs, edc, "Sub_Sensor2", sensors(2), 1);

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//specific Sensor data:

	GetElemDataBool(mbs, edc, "Average", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMAverage);
	GetElemDataBool(mbs, edc, "Sum", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMSum);
	GetElemDataBool(mbs, edc, "Difference", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMDifference);
	GetElemDataBool(mbs, edc, "Mult", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMMult);
	GetElemDataBool(mbs, edc, "Div", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMDiv);
	GetElemDataBool(mbs, edc, "Norm", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMNorm);
  GetElemDataBool(mbs, edc, "Maximum", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMMax);
  GetElemDataBool(mbs, edc, "Minimum", flag, 0);
	if (flag) multtype = (TMBSMultSensor)(type+TSMMin);
	*/

	return rv;
}





double MBSMultipleSensor::GetCurrentValue(double time)
{
	double actvalue = 0; //for case of 0 sensors

	//compute average value of all sensors:
	if (multtype & TSMAverage)
	{
		double av = 0;
		for (int i = 1; i <= sensors.Length(); i++)
		{
			av += GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
		}
		if (sensors.Length() != 0) av /= (double)sensors.Length();

		actvalue = av;
	}
	else if (multtype & TSMSum)
	{
		double av = 0;
		for (int i = 1; i <= sensors.Length(); i++)
		{
			//GetMBS()->UO() << "s" << i << "=" << GetMBS()->GetSensor(sensors(i)).GetValue() << "\n";
			av += GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
		}
		//GetMBS()->UO() << "sum=" << av << "\n";

		actvalue = av;
	}
	else if (multtype & TSMDifference)
	{
		double av = 0;
		for (int i = 1; i <= sensors.Length(); i++)
		{
			if (i == 1) av = GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
			else av -= GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
		}

		actvalue = av;
	}

	else if (multtype & TSMMult)
	{
		double av = 1;
		for (int i = 1; i <= sensors.Length(); i++)
		{
			av *= GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
		}

		actvalue = av;
	}
	else if (multtype & TSMDiv)
	{
		double av = 1;
		for (int i = 1; i <= sensors.Length(); i++)
		{
			if (i == 1) av = GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
			else av /= GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
		}

		actvalue = av;
	}
	else if (multtype & TSMMax)
	{
		double max = GetMBS()->GetSensor(sensors(1)).GetCurrentValueWithSensorProcessing(time);
		for (int i = 2; i <= sensors.Length(); i++)
		{
			double val = GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
			if (max < val) max = val;
		}
		actvalue = max;
	}
	else if (multtype & TSMMin)
	{
		double min = GetMBS()->GetSensor(sensors(1)).GetCurrentValueWithSensorProcessing(time);
		for (int i = 2; i <= sensors.Length(); i++)
		{
			double val = GetMBS()->GetSensor(sensors(i)).GetCurrentValueWithSensorProcessing(time);
			if (min > val) min = val;
		}
		actvalue = min;
	}
	//transformation of return values:
	actvalue = factor * actvalue + offset;

	return actvalue;
}

#endif // __SENSOR_OLD_CPP__