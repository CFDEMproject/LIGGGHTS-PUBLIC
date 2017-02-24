//#**************************************************************
//#
//# filename:             RigidBodyJoints.cpp
//#
//# author:               Saxinger Martin
//#
//# generated:						19. December 2012
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
 
#include "RigidBodyJoints.h"
#include "rendercontext.h"
#include "elementdataaccess.h"
#include "sensors.h"
#include "control.h"

#include "femathhelperfunctions.h"

void BaseBodyJoint::ElementDefaultConstructorInitialization()
{
	BasePointJoint::ElementDefaultConstructorInitialization();

	stiffness_in_joint_local_frame = 0;
	SetUseLocalCoordinateSystem(1);

	IVector dirRot(3);
	dirRot.Set3(1,1,1);
	SetConstrainedDirectionsRot(dirRot);
	JA0i = JA0j = Matrix3D(1);
	penaltyStiffness = Matrix3D(0);
	penaltyDamping = Matrix3D(0);
	penaltyStiffnessRot = Matrix3D(0);
	penaltyDampingRot = Matrix3D(0);

	standard_joint_drawing = 0;
	draw_cone_size = -1;
}

void BaseBodyJoint::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const BaseBodyJoint& ce = (const BaseBodyJoint&)e;

	JA0i = ce.JA0i;
	JA0j = ce.JA0j;
	dirRot = ce.dirRot;
	penaltyStiffness = ce.penaltyStiffness;
	penaltyDamping = ce.penaltyDamping;
	penaltyStiffnessRot = ce.penaltyStiffnessRot;
	penaltyDampingRot = ce.penaltyDampingRot;
	drotTvdq = ce.drotTvdq;
	standard_joint_drawing = ce.standard_joint_drawing;
	draw_cone_size = ce.draw_cone_size;
}

void BaseBodyJoint::SetBaseBodyJoint_LocalPos_to_LocalPos(int elem1, Vector3D lc1, int elem2, Vector3D lc2, IVector constrDir, IVector constrDirRot, Vector3D bryantAnglesI, Matrix3D stiffness, Matrix3D stiffnessRot, Matrix3D damping, Matrix3D dampingRot)
{
	bryantAngles = bryantAnglesI;
	SetPos1ToLocalCoord(elem1, lc1);
	SetPos2ToLocalCoord(elem2, lc2);
	if(!IsNullMatrix(stiffness) || !IsNullMatrix(stiffnessRot)) 
	{
		SetPenaltyFormulation(1);
		SetStiffnessMatrix(stiffness); // sets flag for penalty formulation
		SetStiffnessMatrixRot(stiffnessRot); // sets flag for penalty formulation
		SetDampingMatrix(damping);
		SetDampingMatrixRot(dampingRot);
	} else // no stiffness matrices defined
	{  
		SetConstrainedDirections(constrDir);
		SetConstrainedDirectionsRot(constrDirRot);
	}
}

int BaseBodyJoint::IsNullMatrix(Matrix3D mat) const // return value is 1 if null matrix
{
	double sum = 0;
	for(int i=1; i<=3; i++) 
	{
		for(int j=1; j<=3; j++) 
		{
			sum = sum + abs(mat(i,j));
		}
	}
	if (sum >= 1e-10) return 0;
	return 1;
}

void BaseBodyJoint::Initialize() 
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	}

	if(bryantAngles.Norm()>0)
	{
		JA0i = ComputeRotMatrixWithKardanAngles(bryantAngles(1), bryantAngles(2), bryantAngles(3));
	}
	else
	{
		JA0i = Matrix3D(1.);
	}

	Matrix3D A0i = GetRotMatBodyi();
	Matrix3D A0j = GetRotMatBodyj();
	JA0j = A0j.GetTp()*A0i*JA0i;

	if(UsePenaltyFormulation())
	{
		freeRotPenalty = PenaltyRotModus();
	}
}

int BaseBodyJoint::PenaltyRotModus() const // returns the rotation axis: 0 = rigid constraint, 1 = rot about x-axis, 2 = rot about y-axis, 3 = rot about z-axis, 4 = free rotations
{
	// looking for zero rows and colums
	// based on this the rotation axis can be found
	if(IsNullMatrix(GetStiffnessMatrixRot()) && IsNullMatrix(GetDampingMatrixRot()))
	{
		return 4; // free rotation about all axes
	}

	Matrix3D mat(0,1,1);
	Matrix3D stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
	Matrix3D dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

	if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about x-axis
	{
		return 1;
	}

	mat = Matrix3D(1,0,1);
	stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
	dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

	if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about y-axis
	{
		return 2;
	}

	mat = Matrix3D(1,1,0);
	stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
	dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

	if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about z-axis
	{
		return 3;
	}

	return 0;
}

int BaseBodyJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BasePointJoint::CheckConsistency(errorstr);
	if (rv) return rv;

	// check if invalid dirRot vector for lagrange constraint
	int cnt = 0;

	if(!UsePenaltyFormulation())
	{
		for(int i=1; i<=3; i++)
		{
			if (!dirRot(i)) cnt++;
		}
		if (cnt==2)
		{
				errorstr = mystr("ERROR: BaseBodyJoint: invalid lagrange constraint for rotation!\n");
				rv = 1;
		}
	}
	else
	{
		Matrix3D mat(0,1,1);
		Matrix3D stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
		Matrix3D dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

		if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about x-axis
		{
			cnt++;
		}

		mat = Matrix3D(1,0,1);
		stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
		dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

		if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about y-axis
		{
			cnt++;
		}

		mat = Matrix3D(1,1,0);
		stiffRot = GetStiffnessMatrixRot()-mat*GetStiffnessMatrixRot()*mat;
		dampRot = GetDampingMatrixRot()-mat*GetDampingMatrixRot()*mat;

		if(IsNullMatrix(stiffRot) && IsNullMatrix(dampRot)) // rot about z-axis
		{
			cnt++;
		}

		if (cnt == 2)
		{
				errorstr = mystr("ERROR: BaseBodyJoint: invalid penalty constraint for rotation!\n");
				rv = 1;
		}
	}
	return rv;
}

void BaseBodyJoint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	// penalty formulation
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_angular_velocity+TKAF_D_pos_D_q+TKAF_D_rot_D_q));	
	KinAccFunc.Add(kaf);

	// lagrange formulation
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_pos_D_q+TKAF_D_rot_v_D_q));	
	KinAccFunc.Add(kaf);
}

int BaseBodyJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	// angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	// Manual WRITE/READWRITE entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.joint_bryant_angle",3,0,0.,mystr("prescribe the angles of the joint coordinate system (for actuation, penalty formulation ONLY!)"), TRWElementDataWrite));

	return 0;
}

int BaseBodyJoint::GetSignOfLagrangeParam(int dir)
{
	int sign = 1;
	if (!dirRot(1) && dirRot(2) && dirRot(3)) // rot. about x axis
	{
		if (dir==3) {sign = -1;};
	} else if (dirRot(1) && !dirRot(2) && dirRot(3)) // rot. about y axis
	{
		if (dir==1) {sign = -1;};
	} else if (dirRot(1) && dirRot(2) && !dirRot(3)) // rot. about z axis
	{
		if (dir==2) {sign = -1;};
	} else if (dirRot(1) && dirRot(2) && dirRot(3)) // constrain all rotations
	{
		if (dir == 1 || dir==3) {sign = -1;};
	}
	return sign;
}

int BaseBodyJoint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = BasePointJoint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.constraint_force_global"))
	{
		if (RWdata.comp1 >= 0 && RWdata.comp1 <= 3) //range check
		{
			if(UsePenaltyFormulation())
			{
				if (RWdata.comp1 == 0)
				{
					RWdata.value = ComputeForce(mbs->GetTime()).Norm();
				}
				else
				{
					RWdata.value = ComputeForce(mbs->GetTime())(RWdata.comp1); // compute force returns the force in global coordinate system
				}
				return 1;
			}
			else // return lagrange parameters
			{
				if (RWdata.comp1 == 0)
				{
					double sum = 0;
					int next = 0;
					for(int i=1; i<=3; i++) 
					{
						if(dir(i))
						{
							sum += Sqr(XG(++next));
						}
					}
					RWdata.value = sqrt(sum);
				}
				else
				{
					if (dir(1)!=0 && dir(2)!=0 && dir(3)!=0) // lagrange parameters are already the global force components
					{
						RWdata.value = XG(RWdata.comp1);
					}
					else
					{
						int next = 0;
						Vector3D loc_force(0.);
						Vector3D glob_force;
						for(int i=1; i<=3; i++) 
						{
							if(dir(i)) {loc_force(i) = XG(++next);};
						}
						glob_force = GetRotMati()*loc_force;
						RWdata.value = glob_force(RWdata.comp1);
					}
				}
				return 1; 
			}
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Connector.constraint_force_local"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			if(UsePenaltyFormulation())
			{
				Vector3D loc_force = GetRotMati().GetTp()*ComputeForce(mbs->GetTime());
				RWdata.value = loc_force(RWdata.comp1); // compute force returns the force in global coordinate system
				return 1;
			}
			else // return lagrange parameters
			{
				if (dir(1)!=0 && dir(2)!=0 && dir(3)!=0) // lagrange parameters are the global force components
				{
					Vector3D loc_force = GetRotMati().GetTp()*Vector3D(XG(1),XG(2),XG(3));
					RWdata.value = loc_force(RWdata.comp1);
				}
				else
				{
					if (dir(RWdata.comp1) == 0)
					{
						RWdata.value = 0;
					}
					else
					{
						int next = 0;
						for(int i=1; i<=RWdata.comp1; i++) 
						{
							if(dir(i)) {++next;};
						}
						RWdata.value = XG(next);
					}
				}
				return 1; 
			}
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Connector.constraint_moment_global"))
	{
		if (RWdata.comp1 >= 0 && RWdata.comp1 <= 3) //range check
		{
			if(UsePenaltyFormulation())
			{
				if (RWdata.comp1 == 0)
				{
					RWdata.value = ComputeMoment(mbs->GetTime()).Norm();
				}
				else
				{
					Vector3D tmp = /*GetRotMati()*/ComputeMoment(mbs->GetTime());
					RWdata.value = tmp(RWdata.comp1); 
				}
				return 1;
			}
			else // return lagrange parameters
			{
				// calulate the offset for the lagrange parameters
				int offset = 0;
				for(int i=1; i<=3; i++) 
				{
					if(dir(i)) {offset += 1;};
				}

				if (RWdata.comp1 == 0)
				{
					double sum = 0;
					int next = 0;
					for(int i=1; i<=3; i++) 
					{
						if(dirRot(i))
						{
							++next;
							sum += Sqr(XG(next+offset));
						}
					}
					RWdata.value = sqrt(sum);
				}
				else
				{
					int next = 0;
					Vector3D loc_moment(0.);
					Vector3D glob_moment;
					for(int i=1; i<=3; i++) 
					{
						if(dirRot(i)) 
						{
							++next;
							loc_moment(i) = GetSignOfLagrangeParam(i)*XG(next+offset);
						};
					}
					glob_moment = GetRotMati()*loc_moment;
					RWdata.value = glob_moment(RWdata.comp1);
				}
				return 1; 
			}
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Connector.constraint_moment_local"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			if(UsePenaltyFormulation())
			{
				Vector3D loc_moment = GetRotMati().GetTp()*ComputeMoment(mbs->GetTime());
				RWdata.value = loc_moment(RWdata.comp1); 
				return 1;
			}
			else // return lagrange parameters
			{
				if (dirRot(RWdata.comp1) == 0)
				{
					RWdata.value = 0;
				}
				else
				{
					// calulate the offset for the lagrange parameters
					int offset = 0;
					for(int i=1; i<=3; i++) 
					{
						if(dir(i)) {offset += 1;};
					}
					int next = 0;
					Vector3D loc_moment(0.);
					for(int i=1; i<=3; i++) 
					{
						if(dirRot(i)) 
						{
							++next;
							loc_moment(i) = GetSignOfLagrangeParam(i)*XG(next+offset);
						};
					}
					RWdata.value = loc_moment(RWdata.comp1);
				}

				return 1; 
			}
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Connector.constraint_angle"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			Vector3D phi;

			if(UsePenaltyFormulation())
			{
				Matrix3D A = GetRotMatj().GetTp()*GetRotMati();
				if(!freeRotPenalty) // all rotations constrained
				{
					if(1) // linearized angles
					{
						phi = Vector3D(-A(2,3),A(1,3),-A(1,2));
						phi = -1*phi;
						RWdata.value = phi(RWdata.comp1);
					}
					else // exact bryant angles
					{
						RotMatToKardanAngles(A, phi);
						phi = -1*phi;
						RWdata.value = phi(RWdata.comp1);
					}
					return 1;
				}
				else
				{
					if(freeRotPenalty == 1)
					{
						Matrix3D Mat_transform = Matrix3D(0,0,1,-1,0,0,0,-1,0); // use z-axis as rotation axis
						Matrix3D Mati = GetRotMati()*Mat_transform;
						Matrix3D Matj = GetRotMatj()*Mat_transform;
						Vector3D phi_;
						RotMatToKardanAngles(Matj.GetTp()*Mati, phi_); //in cases of small angles about x or z, the functions ComputeRotMatrixWithKardanAngles and RotMatToKardanAngles lead to the same angles in the end
						phi.Set(phi_(3),-1*phi_(1),-1*phi_(2));
					} 
					else if(freeRotPenalty == 2)
					{
						Matrix3D Mat_transform = Matrix3D(1,0,0,0,0,1,0,-1,0); // use z-axis as rotation axis
						Matrix3D Mati = GetRotMati()*Mat_transform;
						Matrix3D Matj = GetRotMatj()*Mat_transform;
						Vector3D phi_;
						RotMatToKardanAngles(Matj.GetTp()*Mati, phi_); //in cases of small angles about x or z, the functions ComputeRotMatrixWithKardanAngles and RotMatToKardanAngles lead to the same angles in the end
						phi.Set(phi_(1),phi_(3),-1*phi_(2));
					} 
					else
					{
						RotMatToKardanAngles(A, phi); // MSax 2013-02-12, warning: wrong angles when rotation about y-axis
					}
					phi = -1*phi;
					RWdata.value = phi(RWdata.comp1);
					return 1;
				}	 
			}
			else // lagrange formulation
			{
				double angle = 0;

				if (dirRot(1)==0 && dirRot(2) != 0 && dirRot(3) != 0 && RWdata.comp1 == 1) // rot. about x-axis
				{
						Vector3D n1 = GetRotMati()*Vector3D(0,1,0);
						Vector3D n2 = GetRotMatj()*Vector3D(0,1,0);
						Vector3D rot1 = GetRotMati()*Vector3D(1,0,0);
						Vector3D rot2 = GetRotMatj()*Vector3D(1,0,0);
						angle = NormalizedVectorAngle(n1,n2); //phi should always be positive
						if ((0.5*(rot1+rot2))*(n1.Cross(n2)) < 0) angle = -angle; 
				}

				if (dirRot(1) != 0 && dirRot(2) == 0 && dirRot(3) != 0 && RWdata.comp1 == 2) // rot. about y-axis
				{
						Vector3D n1 = GetRotMati()*Vector3D(1,0,0);
						Vector3D n2 = GetRotMatj()*Vector3D(1,0,0);
						Vector3D rot1 = GetRotMati()*Vector3D(0,1,0);
						Vector3D rot2 = GetRotMatj()*Vector3D(0,1,0);
						angle = NormalizedVectorAngle(n1,n2); //phi should always be positive
						if ((0.5*(rot1+rot2))*(n1.Cross(n2)) < 0) angle = -angle;  
				}

				if (dirRot(1) != 0 && dirRot(2) != 0 && dirRot(3) == 0 && RWdata.comp1 == 3) // rot. about z-axis
				{
						Vector3D n1 = GetRotMati()*Vector3D(0,1,0);
						Vector3D n2 = GetRotMatj()*Vector3D(0,1,0);
						Vector3D rot1 = GetRotMati()*Vector3D(0,0,1);
						Vector3D rot2 = GetRotMatj()*Vector3D(0,0,1);
						angle = NormalizedVectorAngle(n1,n2); //phi should always be positive
						if ((0.5*(rot1+rot2))*(n1.Cross(n2)) < 0) angle = -angle; 
				}

				if (dirRot(1)==0 && dirRot(2)==0 && dirRot(3)==0) // free rotation
				{
					RotMatToKardanAngles(GetRotMatj().GetTp()*GetRotMati(), phi);
					phi = -1*phi;
					angle = phi(RWdata.comp1);
				}
				RWdata.value = angle;
				return 1; 
			}
		}
		else return -2; 
	}

	if(RWdata.variable_name.CStrCompare("Connector.constraint_displacement"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			if(UsePenaltyFormulation())
			{
				RWdata.value = (GetRotMati().GetTp()*(GetPos2()-GetPos1()))(RWdata.comp1);
				return 1; 
			}
			else // lagrange formulation
			{
				double displ = 0;
				if (dir(RWdata.comp1) == 0) // displacement in direction RWdata.comp1
				{
					displ = (GetRotMati().GetTp()*(GetPos2()-GetPos1()))(RWdata.comp1);
				}
				RWdata.value = displ;
				return 1; 
			}
		}
		else return -2; 
	}
	return ReadSingleElementDataAuto(RWdata);
}

int BaseBodyJoint::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = Constraint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write
	if(RWdata.variable_name.CStrCompare("Connector.joint_bryant_angle"))
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			bryantAngles(RWdata.comp1) = RWdata.value;
			JA0i = ComputeRotMatrixWithKardanAngles(bryantAngles(1), bryantAngles(2), bryantAngles(3));
			return 1; 
		}
		else return -2; 
	}
	return WriteSingleElementDataAuto(RWdata);
}

int BaseBodyJoint::GetNumberOfConstrainedCoords() const // used for lagrange method
{
	int counter = 0;
	for(int i=1; i<=3; i++)
	{
		if(dir(i)!=0) counter++;
	}
	for(int i=1; i<=3; i++)
	{
		if(dirRot(i)!=0) counter++;
	}
	return counter;
}	

void BaseBodyJoint::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG
	Vector3D delTrans;
	Vector3D delRot;
	Vector3D lv1 = Vector3D(1,0,0);
	Vector3D lv2 = Vector3D(0,1,0);
	Vector3D lv3 = Vector3D(0,0,1);

	Matrix3D roti = GetRotMati();
	Matrix3D rotj = GetRotMatj();
	// g... global coordinate system; l... local coordinate system
	Vector3D gvi1(roti(1,1), roti(2,1), roti(3,1));
	Vector3D gvi2(roti(1,2), roti(2,2), roti(3,2));
	Vector3D gvi3(roti(1,3), roti(2,3), roti(3,3));
	Vector3D gvj1(rotj(1,1), rotj(2,1), rotj(3,1));
	Vector3D gvj2(rotj(1,2), rotj(2,2), rotj(3,2));
	Vector3D gvj3(rotj(1,3), rotj(2,3), rotj(3,3));

	if (MaxIndex()==3)
	{
		if (!IsVelocityConstraint())
		{
			if (dir(1) && dir(2) && dir(3)) // all translations are constrained (no need to transform into local joint i coordinate system)
			{
				delTrans = GetPos1()-GetPos2();
			}
			else
			{
				delTrans = GetRotMati().GetTp()*(GetPos1()-GetPos2());
			}

			if (!dirRot(1) && dirRot(2) && dirRot(3)) // rot. about x axis
			{
				delRot(2) = gvj1 * gvi3; // Myi
				delRot(3) = gvj1 * gvi2; // -Mzi
			} 
			else if (dirRot(1) && !dirRot(2) && dirRot(3)) // rot. about y axis
			{
				delRot(1) = gvj2 * gvi3; // -Mxi
				delRot(3) = gvj2 * gvi1; // Mzi
			} 
			else if (dirRot(1) && dirRot(2) && !dirRot(3)) // rot. about z axis
			{
				delRot(1) = gvj3 * gvi2; // Mxi
				delRot(2) = gvj3 * gvi1; // -Myi
			} 
			else if (dirRot(1) && dirRot(2) && dirRot(3)) // constrain all rotations
			{
				delRot(1) = gvj2 * gvi3; // -Mxi
				delRot(2) = gvj1 * gvi3; // Myi
				delRot(3) = gvj1 * gvi2; // -Mzi
			} 
			else
			{
				delRot.Set(0.,0.,0.);
			}
		}
		else
		{ 
			if(displ.Length()) // if ground is not constant and described by matfunction
			{
				GetMBS()->UO() << "ERROR: BaseBodyJoint: you are using velocity constraints AND displacement";
			}
		}
	}
	else if(MaxIndex() <3 || IsVelocityConstraint())
	{ //velocity constraints:
		// translation
		if (dir(1) && dir(2) && dir(3)) // all translations are constrained (no need to transform into local joint i coordinate system)
		{
			delTrans = GetVel1()-GetVel2();
		}
		else
		{
			delTrans = GetRotMati().GetTp()*(GetVel1()-GetVel2())+GetRotMatiP().GetTp()*(GetPos1()-GetPos2()); // time derivative of delTrans = GetRotMati().GetTp()*(GetPos1()-GetPos2());
		}

		Matrix3D rotip = GetRotMatiP();
		Matrix3D rotjp = GetRotMatjP();
		// g... global coordinate system; l... local coordinate system
		Vector3D gvi1p(rotip(1,1), rotip(2,1), rotip(3,1));
		Vector3D gvi2p(rotip(1,2), rotip(2,2), rotip(3,2));
		Vector3D gvi3p(rotip(1,3), rotip(2,3), rotip(3,3));
		Vector3D gvj1p(rotjp(1,1), rotjp(2,1), rotjp(3,1));
		Vector3D gvj2p(rotjp(1,2), rotjp(2,2), rotjp(3,2));
		Vector3D gvj3p(rotjp(1,3), rotjp(2,3), rotjp(3,3));

		// the following equations are the time derivatives of the indes 3 formulations
		if (!dirRot(1) && dirRot(2) && dirRot(3)) // rot. about x axis
		{
			delRot(2) = gvj1p * gvi3 + gvj1 * gvi3p; // Myi
			delRot(3) = gvj1p * gvi2 + gvj1 * gvi2p; // Mzi
		} 
		else if (dirRot(1) && !dirRot(2) && dirRot(3)) // rot. about y axis
		{
			delRot(1) = gvj2p * gvi3 + gvj2 * gvi3p; // Mxi
			delRot(3) = gvj2p * gvi1 + gvj2 * gvi1p; // Mzi
		} 
		else if (dirRot(1) && dirRot(2) && !dirRot(3)) // rot. about z axis
		{
			delRot(1) = gvj3p * gvi2 + gvj3 * gvi2p; // Mxi
			delRot(2) = gvj3p * gvi1 + gvj3 * gvi1p; // Myi
		} 
		else if (dirRot(1) && dirRot(2) && dirRot(3)) // constrain all rotations
		{
			delRot(1) = gvj2p * gvi3 + gvj2 * gvi3p; // Mxi
			delRot(2) = gvj1p * gvi3 + gvj1 * gvi3p; // Myi
			delRot(3) = gvj1p * gvi2 + gvj1 * gvi2p; // Mzi
		} 
		else
		{
			delRot.Set(0.,0.,0.);
		}
	}

	int next=0;
	for(int i=1; i<=3; i++) 
	{
		if(dir(i))
		{
			f(++next)=delTrans(i);
		}
	}

	for(int i=1; i<=3; i++) 
	{
		if(dirRot(i))
		{
			f(++next)=delRot(i);
		}
	}
};

//$ SW 2013-10-10: trying to avoid as much as possible matrix and vector allocations
void BaseBodyJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	Matrix &hTrans = tmpmat1;
	Matrix &hRot = tmpmat2;
	Vector &hvec = tmpvec1;

	hTrans.SetSize(f.Length(),3);
	hRot.SetSize(f.Length(),3);

	if(IsLocalNodeConstraint())
	{

		if (UsePenaltyFormulation()) {return;}

		double sign = 1;
		if (locelemind == 2) {sign = -1;}


		if (nodes(locelemind)!=0)						// position defined by local node number (also for CMS-Element)
		{
			GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
			return;
		}
		else
		{
			// translation
			GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq); // dpdq is dpdqT

			if (dir(1) && dir(2) && dir(3)) // all directions are constrained
			{
				hTrans = dpdq;
			}
			else
			{
				if (locelemind==1)
				{
					GetdRotTvdqT(GetPos1()-GetPos2(), loccoords(locelemind), drotTvdq, locelemind);
					Matrix3D AiB = GetRotMatBodyi();
					//Matrix JA0i_ = JA0i;
					//hTrans = (drotTvdq+dpdq*AiB_)*JA0i_;
					hTrans = dpdq;
					hTrans.ApplySqrMat(AiB);
					hTrans += drotTvdq;
					hTrans.ApplySqrMat(JA0i);
				}
				else
				{
					Matrix3D Ai = GetRotMati();
					Matrix Ai_ = Ai;
					hTrans = dpdq*Ai_;
				}
			}

			// rotation
			Vector3D lv1 = Vector3D(1,0,0);
			Vector3D lv2 = Vector3D(0,1,0);
			Vector3D lv3 = Vector3D(0,0,1);

			Matrix3D roti = GetRotMati();
			Matrix3D rotj = GetRotMatj();
			// g... global coordinate system; l... local coordinate system
			Vector3D gvi1(roti(1,1), roti(2,1), roti(3,1));
			Vector3D gvi2(roti(1,2), roti(2,2), roti(3,2));
			Vector3D gvi3(roti(1,3), roti(2,3), roti(3,3));
			Vector3D gvj1(rotj(1,1), rotj(2,1), rotj(3,1));
			Vector3D gvj2(rotj(1,2), rotj(2,2), rotj(3,2));
			Vector3D gvj3(rotj(1,3), rotj(2,3), rotj(3,3));
			

			if (locelemind==1)
			{
				if (!dirRot(1) && dirRot(2) && dirRot(3)) // rot. about x axis
				{
					//delRot(2) = gvj1 * gvi3; // Myi
					//delRot(3) = gvj1 * gvi2; // -Mzi
					GetBody3D(1).GetdRotvdqT(JA0i*lv3,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj1,hvec);
					hRot.SetColVec(hvec,2);

					GetBody3D(1).GetdRotvdqT(JA0i*lv2,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj1,hvec);
					hRot.SetColVec(hvec,3);

				} 
				else if (dirRot(1) && !dirRot(2) && dirRot(3)) // rot. about y axis
				{
					//delRot(1) = gvj2 * gvi3; // -Mxi
					//delRot(3) = gvj2 * gvi1; // Mzi
					GetBody3D(1).GetdRotvdqT(JA0i*lv3,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj2,hvec);
					hRot.SetColVec(hvec,1);

					GetBody3D(1).GetdRotvdqT(JA0i*lv1,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj2,hvec);
					hRot.SetColVec(hvec,3);

				} 
				else if (dirRot(1) && dirRot(2) && !dirRot(3)) // rot. about z axis
				{
					//delRot(1) = gvj3 * gvi2; // Mxi
					//delRot(2) = gvj3 * gvi1; // -Myi
					GetBody3D(1).GetdRotvdqT(JA0i*lv2,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj3,hvec);
					hRot.SetColVec(hvec,1);

					GetBody3D(1).GetdRotvdqT(JA0i*lv1,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj3,hvec);
					hRot.SetColVec(hvec,2);

				} 
				else if (dirRot(1) && dirRot(2) && dirRot(3)) // constrain all rotations
				{
					//delRot(1) = gvj2 * gvi3; // -Mxi
					//delRot(2) = gvj1 * gvi3; // Myi
					//delRot(3) = gvj1 * gvi2; // -Mzi
					GetBody3D(1).GetdRotvdqT(JA0i*lv3,loccoords(1),dpdq);  // dpdq is used as drotvdqt
					Mult(dpdq,gvj2,hvec);
					hRot.SetColVec(hvec,1);
					Mult(dpdq,gvj1,hvec);
					hRot.SetColVec(hvec,2);

					GetBody3D(1).GetdRotvdqT(JA0i*lv2,loccoords(1),dpdq);
					Mult(dpdq,gvj1,hvec);
					hRot.SetColVec(hvec,3);
				}

			} 
			else
			{
				if (!dirRot(1) && dirRot(2) && dirRot(3)) // rot. about x axis
				{
					//delRot(2) = gvj1 * gvi3; // Myi
					//delRot(3) = gvj1 * gvi2; // -Mzi
					GetBody3D(2).GetdRotvdqT(JA0j*lv1,loccoords(2),dpdq);
					Mult(dpdq,gvi3,hvec);
					hRot.SetColVec(hvec,2);

					Mult(dpdq,gvi2,hvec);
					hRot.SetColVec(hvec,3);

				} 
				else if (dirRot(1) && !dirRot(2) && dirRot(3)) // rot. about y axis
				{
					//delRot(1) = gvj2 * gvi3; // -Mxi
					//delRot(3) = gvj2 * gvi1; // Mzi
					GetBody3D(2).GetdRotvdqT(JA0j*lv2,loccoords(2),dpdq);
					Mult(dpdq,gvi3,hvec);
					hRot.SetColVec(hvec,1);

					Mult(dpdq,gvi1,hvec);
					hRot.SetColVec(hvec,3);

				} 
				else if (dirRot(1) && dirRot(2) && !dirRot(3)) // rot. about z axis
				{
					//delRot(1) = gvj3 * gvi2; // Mxi
					//delRot(2) = gvj3 * gvi1; // -Myi
					GetBody3D(2).GetdRotvdqT(JA0j*lv3,loccoords(2),dpdq);
					Mult(dpdq,gvi2,hvec);
					hRot.SetColVec(hvec,1);

					Mult(dpdq,gvi1,hvec);
					hRot.SetColVec(hvec,2);

				} 
				else if (dirRot(1) && dirRot(2) && dirRot(3)) // constrain all rotations
				{
					//delRot(1) = gvj2 * gvi3; // -Mxi
					//delRot(2) = gvj1 * gvi3; // Myi
					//delRot(3) = gvj1 * gvi2; // -Mzi
					GetBody3D(2).GetdRotvdqT(JA0j*lv2,loccoords(2),dpdq);
					Mult(dpdq,gvi3,hvec);
					hRot.SetColVec(hvec,1);

					GetBody3D(2).GetdRotvdqT(JA0j*lv1,loccoords(2),dpdq);
					Mult(dpdq,gvi3,hvec);
					hRot.SetColVec(hvec,2);
					Mult(dpdq,gvi2,hvec);
					hRot.SetColVec(hvec,3);
				}
			}
		}

		for (int i=1; i <= f.Length(); i++)
		{
			int next = 0;
			for(int c=1; c<=3; c++)
			{
				if(dir(c))
				{
					f(i) -= sign*(hTrans(i,c)*XG(++next));
				}
			}
			for(int c=1; c<=3; c++)
			{
				if(dirRot(c))
				{
					f(i) -= hRot(i,c)*XG(++next);
				}
			}
		}
	}
	else		// global nodes
	{
			// nothing at all?
	}
};

void BaseBodyJoint::EvalF2(Vector& f, double t)
{
	//f = [f1 f2], where f1 is the residual vector of constraint element1 and f2 of constraint element 2
	//Matrix dpdq;
	double sign = 1.;
	int offset = 0;

	if(!UsePenaltyFormulation()) return;	// no penalty method --> Lagrange multiplier --> no EvalF2
	
	Vector3D force = ComputeForce(t);
	Vector3D moment = ComputeMoment(t);

	//UO(UO_LVL_dbg1) << "force: " << force << "\n";

	for (int i=1; i <= NKinPairs(); i++)
	{
		if (i==2) sign = -1.;
		if(elements(i)!=0)
		{
			//======== get Matrix dpdq 
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)	// CMS
			{
				//GetBody3D(i).GetNodedPosdqT(nodes(i), dpdq);
				GetMBS()->UO() << "ERROR: BaseBodyJoint: not implemented yet";
			}
			else
			{
				if (nodes(i)!=0)						// position defined by node number
				{
					//GetBody3D(i).GetNodedPosdqT(nodes(i), dpdq);
					GetMBS()->UO() << "ERROR: BaseBodyJoint: not implemented yet";
				}
				else
				{
					if (force.Norm2() != 0)
					{
						GetBody3D(i).GetdPosdqT(loccoords(i),dpdq);
					}
					if (moment.Norm2() != 0)
					{
						GetBody3D(i).GetdRotdqT(loccoords(i),drotTvdq); // drotTvdq is here drotdqT
					}
				}
			}

			if (i==2) offset = GetElem(1).SOS();			// has to be changed for mixed global/local node constraints
			for (int j=1; j<=GetElem(i).SOS(); j++)
			{
				if (force.Norm2() != 0)
				{
					f(j+offset) -= sign*(dpdq(j,1)*force.X() + dpdq(j,2)*force.Y() + dpdq(j,3)*force.Z());
				}
				if (moment.Norm2() != 0)
				{
					f(j+offset) -= sign*(drotTvdq(j,1)*moment.X() + drotTvdq(j,2)*moment.Y() + drotTvdq(j,3)*moment.Z());
				}
				//UO(UO_LVL_dbg1) << "i = " << i << ", elem = "  << elements(i) << "ltg(" << j+offset << ") = " << ltg(j+offset) << "\n";
			}
		}
		else
		{
			//if (i==2) offset = 3;
			//f(1+offset) -= sign*force.X();
			//f(2+offset) -= sign*force.Y();
			//f(3+offset) -= sign*force.Z();
			GetMBS()->UO() << "ERROR: BaseBodyJoint: not implemented yet";
		}
	}
};

Vector3D BaseBodyJoint::ComputeForce(double t) const
{
	Vector3D u; // displacement
	Vector3D v; // velocity
	Vector3D f; // resulting force
	Matrix3D rot; // reference rotation
	Matrix3D k_mat;
	Matrix3D d_mat;
	int use_damping = 0;

	if (UsePenaltyFormulation())
	{
		if (!IsNullMatrix(GetDampingMatrix())) {use_damping = 1;};

		u = GetPos1() - GetPos2();
		k_mat = GetStiffnessMatrix();

		if(use_damping) 
		{	
			v = GetVel1() - GetVel2();
			d_mat = GetDampingMatrix();
		}

		// F = A*K*A'*u
		// transform stiffness and damping from local to global coordinate system
		rot = GetRotMati();
		rot.Transpose();
		k_mat = k_mat*rot;
		if(use_damping) {d_mat = d_mat*rot;};
		rot = GetRotMati();
		k_mat = rot*k_mat;
		if(use_damping) {d_mat = rot*d_mat;};

		f = k_mat*u;	
		if(use_damping) {f += d_mat*v;};

		return f;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: ComputeForce just implemented for PenaltyFormulation";
		return Vector3D(0.0);
	}
};

Vector3D BaseBodyJoint::ComputeMoment(double t) const
{
	Vector3D phi; // angle
	Vector3D omega; // angular velocity
	Vector3D m; // resulting moment
	Matrix3D rot; // reference rotation
	Matrix3D k_mat;
	Matrix3D d_mat;
	Matrix3D A, Ap;
	int use_damping = 0;

	if (UsePenaltyFormulation())
	{
		if (!IsNullMatrix(GetDampingMatrixRot())) {use_damping = 1;};

		if(!freeRotPenalty) // all rotations are constrained
		{
			A = GetRotMatj().GetTp()*GetRotMati();
			//compute local linearized angle and angular velocity
			phi = Vector3D(-A(2,3),A(1,3),-A(1,2)); // local angles
		}
		else // constrain about some axis
		{
			if(freeRotPenalty == 1) // rotation about x-axis
			{
				// looking for angles about y and z- axes;
				Matrix3D Mat_transform = Matrix3D(0,0,1,-1,0,0,0,-1,0); // use z-axis as rotation axis to avoid problems with singularity of the rot matrix
				Matrix3D Mati = GetRotMati()*Mat_transform;
				Matrix3D Matj = GetRotMatj()*Mat_transform;
				Vector3D phi_;
				RotMatToKardanAngles(Matj.GetTp()*Mati, phi_);  // phi_ contains the angles about the new joint coordinate system
				phi.Set(phi_(3),-1*phi_(1),-1*phi_(2));  // transform angles to original joint coordinate system
				phi.X() = 0;
			} 
			else if(freeRotPenalty == 2) // rotation about y-axis
			{
				Matrix3D Mat_transform = Matrix3D(1,0,0,0,0,1,0,-1,0); // use z-axis as rotation axis to avoid problems with singularity of the rot matrix
				Matrix3D Mati = GetRotMati()*Mat_transform;
				Matrix3D Matj = GetRotMatj()*Mat_transform;
				Vector3D phi_;
				RotMatToKardanAngles(Matj.GetTp()*Mati, phi_);
				phi.Set(phi_(1),phi_(3),-1*phi_(2));
				phi.Y() = 0;
			} 
			else // freeRotPenalty == 3 // rotation about z-axis
			{
				RotMatToKardanAngles(GetRotMatj().GetTp()*GetRotMati(), phi);
				phi.Z() = 0;
			}
		}

		k_mat = GetStiffnessMatrixRot();

		if(use_damping) 
		{	
			// small linearized velocities
			Ap = GetRotMatjP().GetTp()*GetRotMati()+GetRotMatj().GetTp()*GetRotMatiP();
			Vector3D omega_loc = Vector3D(-Ap(2,3),Ap(1,3),-Ap(1,2));

			if(!freeRotPenalty)
			{
				omega = GetRotMati().GetTp()*GetRotMatiP()*phi+omega_loc;
			}
			else // rotation axis used
			{
				Vector3D angVel1, angVel2; 
				if (elements(1) != 0)
				{
					angVel1 = GetBody3D(1).GetAngularVel(loccoords(1));
				}
				else
				{
					angVel1 = Vector3D(0.);
				}
				if (elements(2) != 0)
				{
					angVel2 = GetBody3D(2).GetAngularVel(loccoords(2));
				}
				else
				{
					angVel2 = Vector3D(0.);
				}

				omega = GetRotMati().GetTp()*(angVel1 - angVel2);

				// set omega about rotation axis to zero because no damping about this rotation axis
				if(freeRotPenalty == 1)
				{
					omega.X() = 0;
				} 
				else if(freeRotPenalty == 2)
				{
					omega.Y() = 0;
				} 
				else // freeRotPenalty == 3
				{
					omega.Z() = 0;
				}
			}
			d_mat = GetDampingMatrixRot();
		}

		m = k_mat*phi;	
		if(use_damping) {m += d_mat*omega;};

		// return global moment vector
		return GetRotMati()*m;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: ComputeMoment just implemented for PenaltyFormulation";
		return Vector3D(0.0);
	}
};

Matrix3D BaseBodyJoint::GetRotMati() const
{
	return GetRotMatBodyi()*JA0i;
}

Matrix3D BaseBodyJoint::GetRotMatiD() const
{
	return GetRotMatrixBodyD(1)*JA0i; // 1 --> Body i
}

Matrix3D BaseBodyJoint::GetRotMatj() const
{
	return GetRotMatBodyj()*JA0j;
}

Matrix3D BaseBodyJoint::GetRotMatjD() const
{
	return GetRotMatrixBodyD(2)*JA0j;
}

Matrix3D BaseBodyJoint::GetRotMatiP() const
{
	return GetRotMatBodyiP()*JA0i;
}

Matrix3D BaseBodyJoint::GetRotMatjP() const
{
	return GetRotMatBodyjP()*JA0j;
}

Matrix3D BaseBodyJoint::GetRotMatrixBody(int i) const
{
	Matrix3D rot;

	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			//pos = (GetMBS()->GetNode(nodes(i))).GetPos();
			GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
			return Matrix3D(0.0);
		}
		else											
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS) // CMS
			{
				//pos = GetBody3D(i).GetNodePos(nodes(i));
				GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
				return Matrix3D(0.0);
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					//pos = GetBody3D(i).GetNodePos(nodes(i));
					GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
					return Matrix3D(0.0);
				}
				else
				{
					rot = GetBody3D(i).GetRotMatrix(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: no element number defined for body 1";
		return Matrix3D(0.0);
	}

	else if (i==2)								// pos2 is ground	
	{
		rot = Matrix3D(1.);				

	}
	return rot;
}

Matrix3D BaseBodyJoint::GetRotMatrixBodyD(int i) const
{
	Matrix3D rot;

	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			//pos = (GetMBS()->GetNode(nodes(i))).GetPos();
			GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
			return Matrix3D(0.0);
		}
		else											
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS) // CMS
			{
				//pos = GetBody3D(i).GetNodePos(nodes(i));
				GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
				return Matrix3D(0.0);
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					//pos = GetBody3D(i).GetNodePos(nodes(i));
					GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
					return Matrix3D(0.0);
				}
				else
				{
					rot = GetBody3D(i).GetRotMatrixD(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: no element number defined for body 1";
		return Matrix3D(0.0);
	}
	else if (i==2)								// pos2 is ground	
	{
		rot = Matrix3D(1.);				
	}
	return rot;
}

Matrix3D BaseBodyJoint::GetRotMatrixBodyP(int i) const
{
	Matrix3D rot;

	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			//pos = (GetMBS()->GetNode(nodes(i))).GetPos();
			GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
			return Matrix3D(0.0);
		}
		else											
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS) // CMS
			{
				//pos = GetBody3D(i).GetNodePos(nodes(i));
				GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
				return Matrix3D(0.0);
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					//pos = GetBody3D(i).GetNodePos(nodes(i));
					GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: not implemented yet";
					return Matrix3D(0.0);
				}
				else
				{
					rot = GetBody3D(i).GetRotMatrixP(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BaseBodyJoint: no element number defined for body 1";
		return Matrix3D(0.0);
	}
	else if (i==2)								// pos2 is ground	
	{
		rot = Matrix3D(0.);				
	}

	return rot;
}

void BaseBodyJoint::DrawElement() 
{
	Constraint::DrawElement();
	//double t_draw = Constraint::GetGlobalTime();

	//if(mbs->GetSimulationStatus() != TSimulationRunning && t_draw == 0.) return;

	Vector3D p;
	Matrix3D rot_mat;
	Vector3D k;
	Vector3D krot;
	int flag = 1;

		if(UsePenaltyFormulation()) 
		{
			k = Vector3D(penaltyStiffness(1,1),penaltyStiffness(2,2),penaltyStiffness(3,3));
			krot = Vector3D(penaltyStiffnessRot(1,1),penaltyStiffnessRot(2,2),penaltyStiffnessRot(3,3));
		}
		else 
		{
			k = Vector3D(dir(1),dir(2),dir(3));		
			krot = Vector3D(dirRot(1),dirRot(2),dirRot(3));
		}

	if(standard_joint_drawing == 0)
	{
		for (int i=1; i <= NKinPairs(); i++)
		{
			if ( i== 2) flag = -1;
			p = GetDrawPosition(i);
			if (i==1)
			{
				rot_mat = GetRotMatiD();
			}
			else
			{
				rot_mat = GetRotMatjD();
			}
			StandardJointDrawing(rot_mat, p, k, krot, flag, GetDrawSizeCone());
		}
	}

	if(draw_local_frame_size > 0 || draw_local_frame_size == -1)
	{
		double s = draw_local_frame_size;
		if(s == -1) 
		{
			s = GetMBS()->GetDOption(104);	// use default value
		}	

		mbs->ChooseColor(0.3f,0.3f,0.3f);

		// size of the frame
		Vector3D v1(-0.2*s, 0,0);
		Vector3D v2( s, 0,0);
		Vector3D v3( 0,-0.2*s,0);
		Vector3D v4( 0, s,0);
		Vector3D v5( 0,0,-0.2*s);
		Vector3D v6( 0,0, s);

		Vector3D v1i,v2i,v3i,v4i,v5i,v6i,v1j,v2j,v3j,v4j,v5j,v6j;
		Matrix3D Ai = GetRotMatiD();
		Vector3D pi = GetDrawPosition(1);

		v1i = pi + Ai*v1;
		v2i = pi + Ai*v2;
		v3i = pi + Ai*v3;
		v4i = pi + Ai*v4;
		v5i = pi + Ai*v5;
		v6i = pi + Ai*v6;

		Matrix3D Aj = GetRotMatjD();
		Vector3D pj = GetDrawPosition(2);

		v1j = pj + Aj*v1;
		v2j = pj + Aj*v2;
		v3j = pj + Aj*v3;
		v4j = pj + Aj*v4;
		v5j = pj + Aj*v5;
		v6j = pj + Aj*v6;

		double d = GetMBS()->GetDOption(114);
		GetMBS()->MyDrawLine(v1i,v2i,d);
		GetMBS()->MyDrawLine(v3i,v4i,d);
		GetMBS()->MyDrawLine(v5i,v6i,d);
		GetMBS()->MyDrawLine(v1j,v2j,d);
		GetMBS()->MyDrawLine(v3j,v4j,d);
		GetMBS()->MyDrawLine(v5j,v6j,d);

		char str[20];
		sprintf(str, "Xi%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v2i.X(), (float)v2i.Y(), (float)v2i.Z(), str);
		sprintf(str, "Yi%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v4i.X(), (float)v4i.Y(), (float)v4i.Z(), str);
		sprintf(str, "Zi%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v6i.X(), (float)v6i.Y(), (float)v6i.Z(), str);
		sprintf(str, "Xj%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v2j.X(), (float)v2j.Y(), (float)v2j.Z(), str);
		sprintf(str, "Yj%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v4j.X(), (float)v4j.Y(), (float)v4j.Z(), str);
		sprintf(str, "Zj%d", GetOwnNum());
		GetMBS()->GetRC()->PrintText3D((float)v6j.X(), (float)v6j.Y(), (float)v6j.Z(), str);
	}
};

double BaseBodyJoint::GetDrawSizeCone()
{
	if(draw_cone_size<=0)
		return GetMBS()->GetDOption(171);			
	else
		return draw_cone_size;
}

//***************************************************************************************************************************************************
// RevoluteJoint
//***************************************************************************************************************************************************
//$ SW 2013-08-29: changed class name from RevoluteJoint_ to RevoluteJoint

void RevoluteJoint::ElementDefaultConstructorInitialization()
{
	BaseBodyJoint::ElementDefaultConstructorInitialization();

	rotAxis =  Vector3D(1.,0.,0.);
	IVector dirRot(3);
	dirRot.Set3(0,1,1);
	SetConstrainedDirectionsRot(dirRot);
	draw_diameter = 0.001;
	draw_axis_length = 0.002;

	damping_coeff = 1e2;
	spring_stiffness = 1e6;

	standard_joint_drawing = 1;
}

void RevoluteJoint::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const RevoluteJoint& ce = (const RevoluteJoint&)e;
	
	rotAxis = ce.rotAxis;
	draw_diameter = ce.draw_diameter;
	draw_axis_length = ce.draw_axis_length;
}

int RevoluteJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BaseBodyJoint::CheckConsistency(errorstr);
	if (rv) return rv;
	return rv;
}

void RevoluteJoint::Initialize() 
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	} 
	Vector3D axis = rotAxis;
	axis.Normalize();
	Vector3D n1, n2;
	axis.SetNormalBasis(n1,n2);
	JA0i.Set(axis,n1,n2);

	Matrix3D A0i = GetRotMatBodyi();
	Matrix3D A0j = GetRotMatBodyj();

	JA0j = A0j.GetTp()*A0i*JA0i;

	if(UsePenaltyFormulation())
	{
		penaltyStiffness = Matrix3D(spring_stiffness);
		penaltyDamping = Matrix3D(damping_coeff);
		penaltyStiffnessRot = Matrix3D(0.,spring_stiffness,spring_stiffness);
		penaltyDampingRot = Matrix3D(0.,damping_coeff,damping_coeff);

		freeRotPenalty = 1; // x-axis is rotation axis
	}
}

void RevoluteJoint::DrawElement() 
{
	BaseBodyJoint::DrawElement();
	
	if(standard_joint_drawing == 1) // StandardJointDrawing
	{
		int res = GetDrawSizeResolution();

		if (NKinPairs()==1)
		{
			Vector3D rot;
			rot = GetRotMatiD()*Vector3D(1,0,0);
			rot.Normalize();
			rot *= 0.5*GetDrawSizeAxisLength();
			Vector3D p = GetDrawPosition(1);
			Vector3D p_ground = GetDrawPosition(2);
			mbs->SetColor(GetCol());
			mbs->DrawZyl(p+rot,p-rot,0.5*GetDrawSizeScalar(),res);
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(p_ground+1.5*rot,p_ground-1.5*rot,0.07*GetDrawSizeScalar(),res);
		} else
		{
			Vector3D p1 = GetDrawPosition(1);
			Vector3D p2 = GetDrawPosition(2);

			Vector3D rot1, rot2;
			rot1 = GetRotMatiD()*Vector3D(1,0,0);
			rot1.Normalize();
			rot1 *= 0.5*GetDrawSizeAxisLength();
			rot2 = GetRotMatjD()*Vector3D(1,0,0);
			rot2.Normalize();
			rot2 *= 0.5*GetDrawSizeAxisLength();

			mbs->SetColor(GetCol());
			mbs->DrawZyl(p1+rot1,p1+rot1*0.1,0.5*GetDrawSizeScalar(),res);
			mbs->DrawZyl(p1-rot1,p1-rot1*0.1,0.5*GetDrawSizeScalar(),res);
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(p2+1.2*rot2,p2-1.2*rot2,0.4*GetDrawSizeScalar(),res);
		}
	}
};


double RevoluteJoint::GetDrawSizeAxisLength()
{
	if(draw_axis_length<=0)
		return GetMBS()->GetDOption(172);
	else
		return draw_axis_length;
}

double RevoluteJoint::GetDrawSizeScalar()
{
	if(draw_diameter<=0)
		return GetMBS()->GetDOption(171);			
	else
		return draw_diameter;
}

int RevoluteJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	// angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	return 0;
}

//***************************************************************************************************************************************************
// PrismaticJoint
//***************************************************************************************************************************************************
//$ SW 2013-08-29: changed class name from PrismaticJoint_ to PrismaticJoint

void PrismaticJoint::ElementDefaultConstructorInitialization()
{
	BaseBodyJoint::ElementDefaultConstructorInitialization();

	slidingDirection =  Vector3D(1.,0.,0.);
	
	IVector dir(3);
	dir.Set3(0,1,1);
	SetConstrainedDirections(dir);

	damping_coeff = 1e2;
	spring_stiffness = 1e6;

	standard_joint_drawing = 1;

	draw_length = 0.02;
	draw_cube_size = Vector3D(0.005);
}

void PrismaticJoint::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const PrismaticJoint& ce = (const PrismaticJoint&)e;
	
	slidingDirection = ce.slidingDirection;
	draw_length = ce.draw_length;
	draw_cube_size = ce.draw_cube_size;
}

int PrismaticJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BaseBodyJoint::CheckConsistency(errorstr);
	if (rv) return rv;
	return rv;
}

void PrismaticJoint::Initialize() 
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	} 
	Vector3D sliding_dir = slidingDirection;
	sliding_dir.Normalize();
	Vector3D n1, n2;
	sliding_dir.SetNormalBasis(n1,n2);
	JA0i.Set(sliding_dir,n1,n2);

	Matrix3D A0i = GetRotMatBodyi();
	Matrix3D A0j = GetRotMatBodyj();

	JA0j = A0j.GetTp()*A0i*JA0i;

	if(UsePenaltyFormulation())
	{
		penaltyStiffness = Matrix3D(0.,spring_stiffness,spring_stiffness);
		penaltyDamping = Matrix3D(0.,damping_coeff,damping_coeff);
		penaltyStiffnessRot = Matrix3D(spring_stiffness);
		penaltyDampingRot = Matrix3D(damping_coeff);

		freeRotPenalty = 0; // no rotation
	}
}

void PrismaticJoint::DrawElement() 
{
	BaseBodyJoint::DrawElement();
	
	if(standard_joint_drawing == 1) // StandardJointDrawing
	{
		Vector3D p1 = GetDrawPosition(1);
		Vector3D p2 = GetDrawPosition(2);

		Vector3D v1 = GetRotMatiD()*Vector3D(1,0,0);
		Vector3D v2 = GetRotMatiD()*Vector3D(0,1,0);
		Vector3D v3 = GetRotMatiD()*Vector3D(0,0,1);

		Vector3D draw_size = GetDrawSize();

		v1 *= 1.2*draw_size(1);
		v2 *= 1.2*draw_size(2);
		v3 *= 1.2*draw_size(3);

		p1 -= 0.5*(v1+v2+v3);

		mbs->SetColor(GetCol());
		mbs->DrawCube(p1, v1,v2,v3);		

		v1 = GetRotMatjD()*Vector3D(1,0,0);
		v2 = GetRotMatjD()*Vector3D(0,1,0);
		v3 = GetRotMatjD()*Vector3D(0,0,1);

		v1 *= GetDrawSizeAxisLength();
		v2 *= draw_size(2);
		v3 *= draw_size(3);

		p2 -= 0.5*(v1+v2+v3);

		mbs->SetColor(GetColExt());
		mbs->DrawCube(p2, v1,v2,v3);	
	}
};

double PrismaticJoint::GetDrawSizeAxisLength()
{
	if(draw_length<=0)
		return GetMBS()->GetDOption(172);
	else
		return draw_length;
}

Vector3D PrismaticJoint::GetDrawSize()
{
	if(draw_cube_size(1) <= 0 || draw_cube_size(2) <= 0 || draw_cube_size(3) <= 0)
	{
		double val = GetMBS()->GetDOption(171);			
		return Vector3D(val);
	}
	else
		return draw_cube_size;
}

int PrismaticJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	// angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	return 0;
}


//***************************************************************************************************************************************************
// UniversalJoint
//***************************************************************************************************************************************************
//$ SW 2013-08-29: added class

void UniversalJoint::ElementDefaultConstructorInitialization()
{
	
	BaseBodyJoint::ElementDefaultConstructorInitialization();
	
	//default initial rotation is about the x-axis
	axis_1 = Vector3D(0.,1.,0.);
	axis_2 = Vector3D(0.,0.,1.);
	col_t = Vector3D(0.2,0.2,0.2);

	draw_direction_1 = Vector3D(1.,0.,0.);
	draw_direction_2 = Vector3D(-1.,0.,0.);

	draw_length = 0.05;
	draw_width = 0.002;

	standard_joint_drawing = 1;
}

void UniversalJoint::SetUniversalJoint(int elem1, Vector3D lc1, Vector3D axis_1, int elem2, Vector3D lc2, Vector3D axis_2)
{
	this->ElementDefaultConstructorInitialization();
	
	SetPos1ToLocalCoord(elem1, lc1);
	SetPos2ToLocalCoord(elem2, lc2);

	this->axis_1 = axis_1;
	this->axis_2 = axis_2;
}

void UniversalJoint::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const UniversalJoint& ce = (const UniversalJoint&)e;

	axis_1 = ce.axis_1;
	axis_2 = ce.axis_2;
	draw_length = ce.draw_length;
	draw_width = ce.draw_width;
}

int UniversalJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BaseBodyJoint::CheckConsistency(errorstr);
	if (rv) return rv;

	if (axis_1.Norm()==0)
	{
		GetMBS()->UO(UO_LVL_err) << "Axis 1 of the cross in the universal joint is (0, 0, 0)";
		rv = 2;
	}
	if (axis_2.Norm()==0)
	{
		GetMBS()->UO(UO_LVL_err) << "Axis 2 of the cross in the universal joint is (0, 0, 0)";
		rv = 2;
	}
	if (rv) return rv;

	Mult(GetRotMatBodyi(),axis_1,tmp_Vec2);
	Vector3D axis_1_global(tmp_Vec2(1),tmp_Vec2(2),tmp_Vec2(3));

	Mult(GetRotMatBodyj(),axis_2,tmp_Vec2);
	Vector3D axis_2_global(tmp_Vec2(1),tmp_Vec2(2),tmp_Vec2(3));

	if (axis_1_global*axis_2_global>0.05)
	{
		//the axes have to be (nearly) orthogonal
		GetMBS()->UO(UO_LVL_err) << "The two axis of the cross in the universal joint must be perpendicular.";
		rv = 1;
	}
	return rv;
}

void UniversalJoint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair){
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	// lagrange formulation
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_rotation_matrix+TKAF_rotation_matrix_P+TKAF_D_pos_D_q+TKAF_D_rot_v_D_q));	
	KinAccFunc.Add(kaf);
}

void UniversalJoint::Initialize()
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	}

	axis_1.Normalize();

	//calculate the global coordinates of the axis
	Mult(GetRotMatBodyi(),axis_1,tmp_Vec1);
	Vector3D axis_1_global(tmp_Vec1(1),tmp_Vec1(2),tmp_Vec1(3));

	Mult(GetRotMatBodyj(),axis_2,tmp_Vec1);
	Vector3D axis_2_global(tmp_Vec1(1),tmp_Vec1(2),tmp_Vec1(3));

	axis_1_global.GramSchmidt(axis_2_global);
	if (axis_2_global.Norm()==0)
	{
		assert(false);
		GetMBS()->UO(UO_LVL_err) << "The two axis of the cross in the universal joint should be perpendicular.";
		axis_2_global(1) = axis_1_global(2);
		axis_2_global(2) = -axis_1_global(1);
		axis_2_global(3) = 0.;
	}

	//Transform the axis back to local coordinates 

	Mult(GetRotMatBodyj().GetTp(),axis_2_global,tmp_Vec1);
	axis_2(1) = tmp_Vec1(1);
	axis_2(2) = tmp_Vec1(2);
	axis_2(3) = tmp_Vec1(3);

	if(UsePenaltyFormulation())
	{
		//NOT IMPLEMENTED YET
		assert (false);
		GetMBS()->UO(UO_LVL_err) << "Penalty formulation for UniversalJoint is not implemented yet";
		return;
	}
}

void UniversalJoint::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) {return;}  // penalty method --> no Lagrange multiplier --> no EvalG
	if (MaxIndex() == 3) // we use positions
	{
		Vector3D delTrans;

		delTrans = GetPos1()-GetPos2();

		for(int i=1; i<=3; i++) 
		{
			f(i)=delTrans(i);
		}

		Mult(GetRotMatBodyi(),axis_1,tmp_Vec1);
		Vector3D axis_1_global(tmp_Vec1(1),tmp_Vec1(2),tmp_Vec1(3));

		Mult(GetRotMatBodyj(),axis_2,tmp_Vec1);
		Vector3D axis_2_global(tmp_Vec1(1),tmp_Vec1(2),tmp_Vec1(3));

		//axis_1_global*axis_2_global
		f(4)=0;
		for (int i=1; i<=3; i++)
		{
			f(4) -= axis_1_global(i)*axis_2_global(i);
		}
	}
	else if (MaxIndex() == 2) //we use velocities
	{
		Vector3D delVel = GetBody3D(1).GetVel(loccoords(1))-GetBody3D(2).GetVel(loccoords(2));
		for(int i=1; i<=3; i++) 
		{
			f(i)=delVel(i);
		}

		//lrot1^T*Ap(1)^T*A(2)*lrot2 + lrot1^T*A(1)^T*Ap(2)*lrot2 = 0
		Matrix3D A1=GetBody3D(1).GetRotMatrix(loccoords(1)); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Matrix3D A2=GetBody3D(2).GetRotMatrix(loccoords(2)); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Matrix3D A1p=GetBody3D(1).GetRotMatrixP(loccoords(1)); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Matrix3D A2p=GetBody3D(2).GetRotMatrixP(loccoords(2)); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);

		Mult(A1,axis_1,tmp_v1);
		Mult(A1p,axis_1,tmp_vp1);
		Mult(A2,axis_2,tmp_v2);
		Mult(A2p,axis_2,tmp_vp2);

		f(4) = tmp_v1*tmp_vp2+tmp_vp1*tmp_v2;
	}
}

//calculates (dC/dq)T*lambda for the element with the local id locelemid
void UniversalJoint::AddElementCqTLambda(double t, int locelemind, Vector& f)
{
	if (UsePenaltyFormulation())
	{
		assert (false);
		GetMBS()->UO(UO_LVL_err) << "Penalty formulation for UniversalJoint is not implemented yet";
		return;
	}

	if (IsLocalNodeConstraint())
	{
		// the sign of element 2 is -1 since the force
		// acts in the opposite direction then
		double sign = 1;
		if (locelemind==2)
		{
			sign = -1;
		}
		//position defined by local node number
		if (nodes(locelemind) != 0)
		{
			Vector3D lambda = Vector3D(XG(1)*sign, XG(2)*sign, XG(3)*sign);
			GetBody3D(locelemind).AddNodedPosdqTLambda(nodes(locelemind), lambda, f);


			//GetBody3D(locelemind).GetdRotdqT(loccoords(locelemind),dRotdqT);
			if (locelemind == 1)
			{
				//get d(Rot*rot_axis_1)dqT
				GetBody3D(1).GetdRotvdqT(axis_1,loccoords(1),tmp_dRotvdqT);
				tmp_Vec1.SetLen(3);
				tmp_Vec2.SetLen(f.Length());
				Mult(GetRotMatBodyj(),axis_2, tmp_Vec1);
				Mult(tmp_dRotvdqT,tmp_Vec1,tmp_Vec2);

				for (int i=1;i <= f.Length();i++)
				{
					f(i) += tmp_Vec2(i)*XG(4);
				}
			}
			else
			{
				//the same again but with negative sign and rotation axis 2
				GetBody3D(2).GetdRotvdqT(axis_2,loccoords(2),tmp_dRotvdqT);
				tmp_Vec1.SetLen(3);
				tmp_Vec2.SetLen(f.Length());
				Mult(GetRotMatBodyi(),axis_1, tmp_Vec1);
				Mult(tmp_dRotvdqT,tmp_Vec1,tmp_Vec2);

				for (int i=1;i <= f.Length();i++)
				{
					f(i) += tmp_Vec2(i)*XG(4);
				}
			}
		}
		else
		{
			// the positions are not given by nodes
			// get dpdq
			GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
			Vector3D lambda = Vector3D(XG(1)*sign, XG(2)*sign, XG(3)*sign);

			//f += dpdq*lambda

			Mult(dpdq,lambda,tmp_Vec1);
			assert(tmp_Vec1.Length()==f.Length());
			for (int i=1; i <= f.Length(); i++)
			{
				f(i) += tmp_Vec1(i);
			}

			//get dRotdqt

			if (locelemind == 1)
			{
				//get d(Rot*rot_axis_1)dqT
				GetBody3D(1).GetdRotvdqT(axis_1,loccoords(1),tmp_dRotvdqT);
				tmp_Vec1.SetLen(3);
				tmp_Vec2.SetLen(f.Length());
				Mult(GetRotMatBodyj(),axis_2, tmp_Vec1);
				Mult(tmp_dRotvdqT,tmp_Vec1,tmp_Vec2);

				for (int i=1;i <= f.Length();i++)
				{
					f(i) += tmp_Vec2(i)*XG(4);

				}
			}
			else
			{
				//the same again but with negative sign and rotation axis 2
				GetBody3D(2).GetdRotvdqT(axis_2,loccoords(2),tmp_dRotvdqT);
				tmp_Vec1.SetLen(3);
				tmp_Vec2.SetLen(f.Length());
				Mult(GetRotMatBodyi(),axis_1, tmp_Vec1);
				Mult(tmp_dRotvdqT,tmp_Vec1,tmp_Vec2);

				for (int i=1;i <= f.Length();i++)
				{
					f(i) += tmp_Vec2(i)*XG(4);
				}
			}
		}
	}
	else
	{
		//not implemented yet
		assert(false);
	}
}


void UniversalJoint::EvalF2(Vector &f, double t)
{
	if (!UsePenaltyFormulation())
	{
		return;
	}

	GetMBS()->UO(UO_LVL_err) << "Penalty formulation for UniversalJoint is not implemented yet";
	assert(false);
}

void UniversalJoint::DrawElement() 
{

	BaseBodyJoint::DrawElement();
	
	if(standard_joint_drawing == 1) // StandardJointDrawing
	{
		Vector3D p1 = GetDrawPosition(1);
		Matrix3D roti = GetRotMatrixBodyD(1);

		Vector3D v1 = roti*draw_direction_1; //unit vector in drawing direction (direction of the length)
		Vector3D a1 = roti*axis_1; //unit vector in axis direction
		a1.Normalize();
		Vector3D t1 = v1.Cross(a1); // unit vector normal to v1 and a1


		Vector3D p2 = GetDrawPosition(2);
		Matrix3D rotj = GetRotMatrixBodyD(2);

		Vector3D v2 = rotj*draw_direction_2; //unit vector in drawing direction (direction of the length)
		Vector3D a2 = rotj*axis_2; //unit vector in axis direction
		a2.Normalize();
		Vector3D t2 = v2.Cross(a2); // unit vector normal to v2 and a2

		if (v1*a1>0.1 || v2*a2>0.1)
		{
			mbs->UO(UO_LVL_err) << "Problem drawing UniversalJoint: the drawing direction and the axis should be perpendicular.";
			return;
		}

		mbs->SetColor(col);
		mbs->DrawCube(p1-v1*draw_length-a1*(draw_width*0.5)-t1*(draw_width*0.5),v1*(draw_length*0.5),a1*draw_width,t1*draw_width);
		mbs->DrawCube(p1-v1*(draw_length*0.5)-a1*(draw_width*0.5)-t1*(draw_width*0.25),v1*(draw_length*0.5+draw_width*0.2),a1*(draw_width*0.2),t1*(draw_width*0.5));
		mbs->DrawCube(p1-v1*(draw_length*0.5)+a1*(draw_width*0.5)-t1*(draw_width*0.25),v1*(draw_length*0.5+draw_width*0.2),a1*(draw_width*-0.2),t1*(draw_width*0.5));
		mbs->SetColor(col_t);
		mbs->DrawZyl(p1-a1*(draw_width*0.49), p1+a1*(draw_width*0.49), draw_width*0.2);

		mbs->SetColor(col_ext);
		mbs->DrawCube(p2-v2*draw_length-a2*(draw_width*0.5)-t2*(draw_width*0.5),v2*(draw_length*0.5),a2*draw_width,t2*draw_width);
		mbs->DrawCube(p2-v2*(draw_length*0.5)-a2*(draw_width*0.5)-t2*(draw_width*0.25),v2*(draw_length*0.5+draw_width*0.2),a2*(draw_width*0.2),t2*(draw_width*0.5));
		mbs->DrawCube(p2-v2*(draw_length*0.5)+a2*(draw_width*0.5)-t2*(draw_width*0.25),v2*(draw_length*0.5+draw_width*0.2),a2*(draw_width*-0.2),t2*(draw_width*0.5));
		mbs->SetColor(col_t);
		mbs->DrawZyl(p2-a2*(draw_width*0.49), p2+a2*(draw_width*0.49), draw_width*0.2);
	}
}

double UniversalJoint::GetDrawSizeAxisLength()
{
	//REMOVE IF METHOD IS PROPERLY IMPLEMENTED
//	assert(false);

	if(draw_length<=0)
		return GetMBS()->GetDOption(172);
	else
		return draw_length;
}

Vector3D UniversalJoint::GetDrawSize()
{
	//REMOVE IF METHOD IS PROPERLY IMPLEMENTED
	//assert(false);
	return Vector3D(2*draw_length,draw_width,draw_width);
}

int UniversalJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Constraint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BasePointJoint::GetAvailableSpecialValuesAuto(available_variables);


	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// displacement
	//available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	return 0;
}

int UniversalJoint::IS() const
{
	if (UsePenaltyFormulation())
	{
		return 0;
	}
	else
	{
		return 4;
	}
}


//***************************************************************************************************************************************************
// RigidJoint
//***************************************************************************************************************************************************
//$ SW 2013-08-29: changed class name from RigidJoint_ to RigidJoint

void RigidJoint::ElementDefaultConstructorInitialization()
{
	BaseBodyJoint::ElementDefaultConstructorInitialization();

	damping_coeff = 1e2;
	spring_stiffness = 1e6;

	standard_joint_drawing = 1;
	draw_dimension = 0.01;
}

void RigidJoint::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const RigidJoint& ce = (const RigidJoint&)e;
	draw_dimension = ce.draw_dimension;
}

int RigidJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BaseBodyJoint::CheckConsistency(errorstr);
	if (rv) return rv;
	return rv;
}

void RigidJoint::Initialize() 
{
	BaseBodyJoint::Initialize();
	if(UsePenaltyFormulation())
	{
		penaltyStiffness = Matrix3D(spring_stiffness);
		penaltyDamping = Matrix3D(damping_coeff);
		penaltyStiffnessRot = penaltyStiffness;
		penaltyDampingRot = penaltyDamping;

		freeRotPenalty = 0; // no rotation
	}
}

void RigidJoint::DrawElement() 
{
	BaseBodyJoint::DrawElement();
	
	mbs->SetColor(GetCol());
	int res = (int)GetDrawSizeResolution();
	
	if(standard_joint_drawing == 1) // StandardJointDrawing
	{
		Vector3D pi;
		Vector3D pj;

		/*Vector3D roti1, roti2, roti3;
		Vector3D rotj1, rotj2, rotj3;

		Matrix3D rotmati = GetRotMatiD();

		roti1 = GetRotMatiD()*Vector3D(1,0,0);
		roti1 *= GetDrawSizeScalar();
		roti2 = GetRotMatiD()*Vector3D(0,1,0);
		roti2 *= GetDrawSizeScalar();
		roti3 = GetRotMatiD()*Vector3D(0,0,1);
		roti3 *= GetDrawSizeScalar();*/

		//$ AH 2013-12: faster
		Matrix3D rotmati = GetRotMatiD();
		Vector3D roti1(rotmati(1,1), rotmati(2,1), rotmati(3,1));
		Vector3D roti2(rotmati(1,2), rotmati(2,2), rotmati(3,2));
		Vector3D roti3(rotmati(1,3), rotmati(2,3), rotmati(3,3));
		roti1 *= GetDrawSizeScalar(); 
		roti2 *= GetDrawSizeScalar();
		roti3 *= GetDrawSizeScalar();

		/*rotj1 = GetRotMatjD()*Vector3D(1,0,0);
		rotj2 = GetRotMatjD()*Vector3D(0,1,0);
		rotj3 = GetRotMatjD()*Vector3D(0,0,1);*/

		Matrix3D rotmatj = GetRotMatjD();
		Vector3D rotj1(rotmatj(1,1), rotmatj(2,1), rotmatj(3,1));
		Vector3D rotj2(rotmatj(1,2), rotmatj(2,2), rotmatj(3,2));
		Vector3D rotj3(rotmatj(1,3), rotmatj(2,3), rotmatj(3,3));

		roti1 *= 1.05;
		rotj1 *= 1.00*GetDrawSizeScalar();
		rotj2 *= 0.95*GetDrawSizeScalar();
		rotj3 *= 1.05*GetDrawSizeScalar();

		pi = GetDrawPosition(1);
		pj = GetDrawPosition(2);

		pi -= 0.5*(roti1+roti2+roti3);
		pj -= 0.5*(rotj1+rotj2+rotj3);

		mbs->SetColor(GetCol());
		mbs->DrawCube(pi, roti1, roti2, roti3);
		mbs->SetColor(GetColExt());
		mbs->DrawCube(pj, rotj1, rotj2, rotj3);
	}
};


double RigidJoint::GetDrawSizeScalar()
{
	if(draw_dimension<=0)
		return GetMBS()->GetDOption(171);			
	else
		return draw_dimension;
}

int RigidJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	// angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	return 0;
}

//***************************************************************************************************************************************************
// CylindricalJoint
//***************************************************************************************************************************************************
//$ SW 2013-08-29: changed class name from CylindricalJoint_ to CylindricalJoint

void CylindricalJoint::ElementDefaultConstructorInitialization()
{
	BaseBodyJoint::ElementDefaultConstructorInitialization();

	rotSlideAxis =  Vector3D(1.,0.,0.);

	IVector dir_dirRot(3);
	dir_dirRot.Set3(0,1,1);

	SetConstrainedDirectionsRot(dir_dirRot);
	SetConstrainedDirections(dir_dirRot);
	draw_cylinder_size = Vector2D(0.01);
	draw_axis_length = 0.02;

	damping_coeff = 1e2;
	spring_stiffness = 1e6;

	standard_joint_drawing = 1;
}

void CylindricalJoint::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const CylindricalJoint& ce = (const CylindricalJoint&)e;
	
	rotSlideAxis = ce.rotSlideAxis;
	draw_cylinder_size = ce.draw_cylinder_size;
	draw_axis_length = ce.draw_axis_length;
}

int CylindricalJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = BaseBodyJoint::CheckConsistency(errorstr);
	if (rv) return rv;
	return rv;
}

void CylindricalJoint::Initialize() 
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	} 
	Vector3D rot_sliding_axis = rotSlideAxis;
	rot_sliding_axis.Normalize();
	Vector3D n1, n2;
	rot_sliding_axis.SetNormalBasis(n1,n2);
	JA0i.Set(rot_sliding_axis,n1,n2);

	Matrix3D A0i = GetRotMatBodyi();
	Matrix3D A0j = GetRotMatBodyj();

	JA0j = A0j.GetTp()*A0i*JA0i;

	if(UsePenaltyFormulation())
	{
		penaltyStiffness = Matrix3D(0.,spring_stiffness,spring_stiffness);
		penaltyDamping = Matrix3D(0.,damping_coeff,damping_coeff);
		penaltyStiffnessRot = penaltyStiffness;
		penaltyDampingRot = penaltyDamping;

		freeRotPenalty = 1; // x-axis is rotation axis
	}
}

void CylindricalJoint::DrawElement() 
{
	BaseBodyJoint::DrawElement();
	
	mbs->SetColor(GetCol());
	int res = (int)GetDrawSizeResolution();
	Vector2D draw_size = GetDrawSize();
	
	if(standard_joint_drawing == 1) // StandardJointDrawing
	{
		if (NKinPairs()==1)
		{
			Vector3D rot, rot2;
			rot = GetRotMatiD()*Vector3D(1,0,0);
			rot.Normalize();
			rot *= 0.5*GetDrawSizeAxisLength();	
			rot2 = GetRotMatjD()*Vector3D(1,0,0);
			rot2.Normalize();
			rot2 *= draw_size(1);	

			Vector3D p = GetDrawPosition(1);
			Vector3D p_ground = GetDrawPosition(2);
			mbs->SetColor(GetCol());
			mbs->DrawZyl(p+rot2,p-rot2,0.5*draw_size(2),res);
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(p_ground+1.5*rot,p_ground-1.5*rot,0.35*draw_size(2),res); 
		} else
		{
			Vector3D p1 = GetDrawPosition(1);
			Vector3D p2 = GetDrawPosition(2);

			Vector3D rot1, rot2;
			rot1 = GetRotMatiD()*Vector3D(1,0,0);
			rot1.Normalize();
			rot1 *= 0.5*draw_size(1);
			rot2 = GetRotMatjD()*Vector3D(1,0,0);
			rot2.Normalize();
			rot2 *= 0.5*GetDrawSizeAxisLength();

			mbs->SetColor(GetCol());
			mbs->DrawZyl(p1-rot1,p1+rot1,0.5*draw_size(2),res);
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(p2+1.5*rot2,p2-1.5*rot2,0.35*draw_size(2),res);
		}
	}
};

double CylindricalJoint::GetDrawSizeAxisLength()
{
	if(draw_axis_length<=0)
		return GetMBS()->GetDOption(172);
	else
		return draw_axis_length;
}

Vector2D CylindricalJoint::GetDrawSize()
{
	if(draw_cylinder_size(1) <= 0 || draw_cylinder_size(2) <= 0)
	{
		double val = GetMBS()->GetDOption(171);
		return Vector2D(val);
	}
	else
		return draw_cylinder_size;
}

int CylindricalJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	BasePointJoint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BaseBodyJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	// angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	return 0;
}

//***************************************************************************************************************************************************
// SpringDamperActuator
//***************************************************************************************************************************************************

void SpringDamperActuator::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const SpringDamperActuator& ce = (const SpringDamperActuator&)e;

	l0 = ce.l0;
	fa = ce.fa;

	forcemode = ce.forcemode;
	mathfunc_k = ce.mathfunc_k;
	mathfunc_d = ce.mathfunc_d;
	mathfunc_fk = ce.mathfunc_fk;
	mathfunc_fd = ce.mathfunc_fd;
	force_k = ce.force_k;
	force_d = ce.force_d;

	spring_res = ce.spring_res;
}

void SpringDamperActuator::ElementDefaultConstructorInitialization()
{
	BasePointJoint::ElementDefaultConstructorInitialization();

	SetDampingCoeff(1);
	SetPenaltyStiffness(1e2);
	SetPenaltyFormulation(1);
	l0 = 0;
	fa = 0;
	forcemode = 0;
	draw_dim = Vector3D(0.01,10,0.01);
	spring_res = 5;

	mathfunc_k = MathFunction();
	mathfunc_d = MathFunction();
	mathfunc_fk = MathFunction();
	mathfunc_fd = MathFunction();

	force_k = 0;
	force_d = 0;

	auto_comp_ground = 0;
}

void SpringDamperActuator::Initialize()
{
	//if(auto_comp_ground)
	//{
	//	SetPos2ToGlobalCoord(GetPos1());
	//}	
}

void SpringDamperActuator::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	// penalty + elemNr + loc node
	kaf = (int)(TKinematicsAccessFunctions(TKAF_node_position+TKAF_node_velocity+TKAF_D_node_pos_D_q));	
	KinAccFunc.Add(kaf);
	
	// penalty + elemNr + locCoord
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_D_pos_D_q));	
	KinAccFunc.Add(kaf);
}

void SpringDamperActuator::SetStiffnessDamping()
{
	Matrix data;
	int pos = 0;
	if (pos) mathfunc_k.SetData(TMFpiecewiselinear, data);
	if (pos) mathfunc_d.SetData(TMFpiecewiselinear, data);
}

void SpringDamperActuator::SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force)
{
	forcemode = 1;
	mathfunc_k.SetPiecewise(xpos, force, 1);
}
void SpringDamperActuator::SetPiecewiseLinearDamper(const Vector& vel, const Vector& force)
{
	forcemode = 1;
	mathfunc_d.SetPiecewise(vel, force, 1);
}

void SpringDamperActuator::SetPiecewiseLinearSpringForce(const Vector& xpos, const Vector& force)
{
	forcemode = 3;
	mathfunc_fk.SetPiecewise(xpos, force, 1);
}
void SpringDamperActuator::SetPiecewiseLinearDamperForce(const Vector& vel, const Vector& force)
{
	forcemode = 3;
	mathfunc_fd.SetPiecewise(vel, force, 1);
}

void SpringDamperActuator::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	
	Vector3D pa = GetDrawPosition(1);
	Vector3D pb = GetDrawPosition(2);

	if (NKinPairs()==1) Swap(pa,pb);

	if (draw_dim.Y() == 0) // if coils == 0 ==> draw cylinder with spring diameter
	{
		mbs->DrawZyl(pa, pb, 0.5*draw_dim.X(),12);
	}
	else
	{
		//draw spring:
		double phi,phi2,xx;
		Vector3D vf;
		double steps = 200*0.5*spring_res; //0.2 ->> very coarse
		double rots = draw_dim.Y();
		double r = draw_dim.X()/2.;
		double rz = draw_dim.X()/10.;
		int ztile = 6;

		vf = pa-pb;
		double lenvf = vf.Norm();
		Vector3D dir = vf;
		vf = ComputeForceDirectionD();
		vf.Normalize();
		vf *= lenvf;
		if (vf*dir < 0) vf *= -1;

		Vector3D vfn = vf;
		Vector3D vz,vy;
		vfn.SetNormalBasis(vz,vy);
		Vector3D pz1, pz2;

		if(GetPenaltyStiffness() || (forcemode > 0)) //$ DR 2013-05-03 linear spring damper with stiffness = 0, do not draw spring
		{
			double off = steps*0.1;
			mbs->SetColor(GetCol());
			for (xx = off; xx <= steps-off; xx++)
			{
				phi = xx*2*MY_PI/steps*rots;
				phi2 = (xx+1.2)*2*MY_PI/steps*rots;

				Vector3D pr1=cos(phi)*r*vz+sin(phi)*r*vy;
				Vector3D pr2=cos(phi2)*r*vz+sin(phi2)*r*vy;

				pz1=pb+xx/steps*vf+pr1;
				pz2=pb+(xx+1)/steps*vf+pr2;
				mbs->DrawZyl(pz1,pz2,rz,ztile);
				if (xx==off) 
				{
					mbs->DrawZyl(pb,pz1,rz,ztile);
					mbs->DrawSphere(pz1,1.05*rz,12);
				}
			}
			mbs->DrawZyl(pa,pz2,rz,ztile);
			mbs->DrawSphere(pz2,1.05*rz,12);
		}

		if (draw_dim.Z() != 0 && GetDampingCoeff() != 0)
		{
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(pa, pa-0.5*l0*vfn, 0.5*draw_dim.Z(),12);
			mbs->DrawZyl(pa-0.5*l0*vfn, pb, 0.5*draw_dim.Z()*0.6,12);

			mbs->DrawSphere(pa,0.65*draw_dim.Z(),12);
			mbs->DrawSphere(pb,0.65*draw_dim.Z(),12);
		}
	}
};

Vector3D SpringDamperActuator::ComputeForce(double t)  const
{
	TMStartTimer(28);
	double l = 0;
	double xp = 0;
	Vector3D dir = ComputeForceDirection();
	double f = 0;

	if (UsePenaltyFormulation())
	{
		Vector3D v12 = GetPos1() - GetPos2();
		l = v12*dir;
		//l = v12.Norm(); // MSax: correct but not efficient
		if(UseDamping()) {	
			Vector3D v12p = GetVel1() - GetVel2();
			xp = v12p*dir;
			//xp = (v12*v12p)/l; // MSax: correct but not efficient
		}

		if (forcemode == 0) // linear stiffness and damper
		{
			double x = (l-l0);
			f = GetPenaltyStiffness()*x;	
			if(UseDamping()) { f += GetDampingCoeff() * xp; }
		}
		else if (forcemode == 1) // nonlinear stiffness and damping by MathFunction
		{
			double x = (l-l0);
			if (mathfunc_k.GetFuncMode() == TMFpiecewiselinear)
			{ //mathfunc_k ... stiffness depending on 'x'
				f = mathfunc_k.Evaluate(x)*x;
			}
			else
			{
				GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: The stiffness function must be piecewise linear";
			}
			if (mathfunc_d.GetFuncMode() == TMFpiecewiselinear)
			{ //mathfunc_k ... damping depending on 'xp'
				f += mathfunc_d.Evaluate(xp)*xp;
			}
			else
			{
				GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: The damping function must be piecewise linear";
			}
		}
		//$ SW 2013-11-20: Added forcemode 3
		else if (forcemode ==3) // nonlinear spring and damping force by MathFunction
		{
			double x = (l-l0);
			if (mathfunc_fk.GetFuncMode() == TMFpiecewiselinear)
			{ //mathfunc_fk ... spring force depending on 'x'
				f = mathfunc_fk.Evaluate(x);
			}
			else
			{
				GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: The spring force function must be piecewise linear";
			}
			if (mathfunc_fd.GetFuncMode() == TMFpiecewiselinear)
			{ //mathfunc_fd ... damping force depending on 'xp'
				f += mathfunc_fd.Evaluate(xp);
			} 
			else
			{
				GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: The damping force function must be piecewise linear";
			}
		}
		else
		{ // nonlinear stiffness and damping
			f =	force_k + force_d;
		}

		f += fa; // add constant actor force
		TMStopTimer(28);
		return f*dir;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: ComputeForce just implemented for PenaltyFormulation";
		TMStopTimer(28);
		return Vector3D(0.0);
	}
	TMStopTimer(28);
}

Vector3D SpringDamperActuator::ComputeForceDirection() const   //return global force direction, normalized
{
	Vector3D v12;
	v12 = GetPos1() - GetPos2();
	v12.Normalize();
	
	return v12;
}

Vector3D SpringDamperActuator::ComputeForceDirectionD()   //return global force direction, normalized
{
	Vector3D v12;
	v12 = GetDrawPosition(1)-GetDrawPosition(2);
	v12.Normalize();

	return v12;
}

void SpringDamperActuator::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	GetElementDataAuto(edc);

	// MathFunction
	ElementDataContainer edc_mf_k;
	mathfunc_k.GetElementData(edc_mf_k);
	ed.SetEDC(&edc_mf_k,"MathFunction_k"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the stiffness as function of l-l0."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	

	ElementDataContainer edc_mf_d;
	mathfunc_d.GetElementData(edc_mf_d);
	ed.SetEDC(&edc_mf_d,"MathFunction_d"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the damping as function of v."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);

	ElementDataContainer edc_mf_fk;
	mathfunc_fk.GetElementData(edc_mf_fk);
	ed.SetEDC(&edc_mf_fk,"MathFunction_fk"); ed.SetToolTipText("MathFunction is used, if forcemode = 3. Gives the spring force as function of l-l0."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	

	ElementDataContainer edc_mf_fd;
	mathfunc_fd.GetElementData(edc_mf_fd);
	ed.SetEDC(&edc_mf_fd,"MathFunction_fd"); ed.SetToolTipText("MathFunction is used, if forcemode = 3. Gives the damping force as function of v."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);
}

int SpringDamperActuator::SetElementData(ElementDataContainer& edc) 		//update element data
{
	int rv = 1;
	ElementData* mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_k"));
	if (mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* mathFuncEdc = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *mathFuncEdc;
		mathfunc_k.SetElementData(mbs,edc_mf);
	}

	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_d"));
	if (mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* mathFuncEdc = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *mathFuncEdc;
		mathfunc_d.SetElementData(mbs,edc_mf);
	}

	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_fk"));
	if (mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* mathFuncEdc = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *mathFuncEdc;
		mathfunc_fk.SetElementData(mbs,edc_mf);
	}

	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_fd"));
	if (mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* mathFuncEdc = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *mathFuncEdc;
		mathfunc_fd.SetElementData(mbs,edc_mf);
	}

	rv = SetElementDataAuto(edc);
	return rv;
}

int SpringDamperActuator::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// Automatic entries for this class 
	BasePointJoint::GetAvailableSpecialValues(available_variables);
	SpringDamperActuator::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_acting_force",0,0,0.,mystr("internal resultant force of connector"), TRWElementDataRead));

	available_variables.Add(ReadWriteElementDataVariableType("Connector.spring_length_offset",0,0,0.,mystr("prescribe the neutral spring length"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.spring_force",0,0,0.,mystr("prescribe the stiffness force"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.damper_force",0,0,0.,mystr("prescribe the damping force"), TRWElementDataWrite));

	return 0;
}

int SpringDamperActuator::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = BasePointJoint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.constraint_acting_force"))
	{
		RWdata.value = ComputeForce(mbs->GetTime())*ComputeForceDirection();
		return 1;
	}

	return ReadSingleElementDataAuto(RWdata);
}

int SpringDamperActuator::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = BasePointJoint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write
	if(RWdata.variable_name.CStrCompare("Connector.spring_length_offset"))
	{
		l0 = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.spring_force"))
	{
		force_k = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.damper_force"))
	{
		force_d = RWdata.value;
		return 1; 
	}

	return WriteSingleElementDataAuto(RWdata);
}

//***************************************************************************************************************************************************
// RigidLink
//***************************************************************************************************************************************************

void RigidLink::CopyFrom(const Element& e)
{
	BasePointJoint::CopyFrom(e);
	const RigidLink& ce = (const RigidLink&)e;

	distancemode = ce.distancemode;

	distance = ce.distance;

	mathfunc_l = ce.mathfunc_l;
	mathfunc_v = ce.mathfunc_v;

	IO_dist = ce.IO_dist;
	IO_vel = ce.IO_vel;

}

void RigidLink::ElementDefaultConstructorInitialization()
{
	BasePointJoint::ElementDefaultConstructorInitialization();

	SetPenaltyFormulation(0);
	distance = 0;
	distancemode = 0;
	draw_dim = Vector3D(0.01,0,0);

	mathfunc_l = MathFunction();
	mathfunc_v = MathFunction();

	IO_dist = 0;
	IO_vel = 0;

	elementname = GetElementSpec();
}

void RigidLink::Initialize()
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	}	
}

int RigidLink::SOS() const // MSax: SOS() from class BasePointJoint cannot be used, see *) because vector dir is not used
{
	if(IsLocalNodeConstraint()) 
	{
		return 0;
	}
	else 
	{
		return 1; // *) return value would be 0 because vector dir is not used
	}
}

void RigidLink::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	// lagrange + elemNr + loc node
	kaf = (int)(TKinematicsAccessFunctions(TKAF_node_position+TKAF_node_velocity+TKAF_D_node_pos_D_q));	
	KinAccFunc.Add(kaf);
	
	// lagrange + elemNr + locCoord
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_D_pos_D_q));	
	KinAccFunc.Add(kaf);
}

void RigidLink::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;

	GetElementDataAuto(edc);

	// MathFunction
	ElementDataContainer edc_mf_l;
	mathfunc_l.GetElementData(edc_mf_l);
	ed.SetEDC(&edc_mf_l,"MathFunction_l"); ed.SetToolTipText("MathFunction is used, if distancemode = 1. Gives the distance at the time t."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	

	ElementDataContainer edc_mf_v;
	mathfunc_v.GetElementData(edc_mf_v);
	ed.SetEDC(&edc_mf_v,"MathFunction_v"); ed.SetToolTipText("MathFunction is used, if distancemode = 1. Gives the derivative of the distance with respect to time at the time t."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	
}

int RigidLink::SetElementData(ElementDataContainer& edc) 		//update element data
{
	int rv = 1;

	ElementData* mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_l"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf_v = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf_l = *edcp_mf_v;
		mathfunc_l.SetElementData(mbs,edc_mf_l);
	}

	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_v"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf_v = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf_v = *edcp_mf_v;
		mathfunc_v.SetElementData(mbs,edc_mf_v);
	}

	SetElementDataAuto(edc);
	return rv;
}

void RigidLink::DrawElement() 
{
	Constraint::DrawElement();
	
	Vector3D pa = GetDrawPosition(1);
	Vector3D pb = GetDrawPosition(2);

	//if (NKinPairs()==1) Swap(pa,pb); // cylinder near ground point should be bigger than the other

	mbs->SetColor(GetCol());
	if (draw_dim.Y() == 0 || draw_dim.Z() == 0 || distancemode == 0) // draw one cylinder
	{
		mbs->DrawZyl(pa, pb, 0.5*draw_dim.X(),12);
		mbs->DrawSphere(pa,0.65*draw_dim.X(),12);
		mbs->DrawSphere(pb,0.65*draw_dim.X(),12);
	}
	else
	{
		Vector3D dirD = ComputeForceDirectionD();
		mbs->DrawZyl(pa, pa-draw_dim.Z()*dirD, 0.5*draw_dim.X(),12);
		mbs->DrawSphere(pa,0.65*draw_dim.X(),12);

		mbs->SetColor(GetColExt());
		mbs->DrawZyl(pa, pb, 0.5*draw_dim.Y(),12);
		mbs->DrawSphere(pb,0.65*draw_dim.Y(),12);
	}
};

Vector3D RigidLink::ComputeForceDirection() const   //return global force direction, normalized
{
	Vector3D v12;
	v12 = GetPos1() - GetPos2();
	v12.Normalize();
	
	return v12;
}

Vector3D RigidLink::ComputeForceDirectionD()   //return global force direction, normalized
{
	Vector3D v12;
	v12 = GetDrawPosition(1)-GetDrawPosition(2);
	v12.Normalize();

	return v12;
}

void RigidLink::EvalG(Vector& f, double t)  
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG

	// get distance at time t
	double offset;
	if (MaxIndex()==3 && !IsVelocityConstraint()) // offset is the distance between the points
	{
		if (distancemode == 0) // constant distance
		{
			offset = distance;
		}
		else if (distancemode == 1)
		{
			if (mathfunc_l.GetFuncMode() == TMFpiecewiselinear) {offset = mathfunc_l.Evaluate(t);} // distance at time t
		} else
		{
			offset = IO_dist;
		}
	}
	else // offset is the velocity between the points
	{
		if (distancemode == 0) // constant distance
		{
			offset = 0;
		}
		else if (distancemode == 1)
		{
			if (mathfunc_v.GetFuncMode() == TMFpiecewiselinear) {offset = mathfunc_v.Evaluate(t);} // distance at time t
		} else
		{
			offset = IO_vel;
		}
	}

	double l = (GetPos1()-GetPos2()).Norm();

	if (MaxIndex()==3)
	{
		if (!IsVelocityConstraint())
		{
			f(1) = l-offset; // position constraint
		}
		else
		{ //velocity constraints:
			f(1) = ((GetPos1()-GetPos2())*(GetVel1()-GetVel2()))/l-offset;
			if(displ.Length()) 
			{
				GetMBS()->UO() << "ERROR: RigidLink: you are using velocity constraints AND displacement";
			}
		}
	}
	else	// MaxIndex < 3
	{ //velocity constraints:
		f(1) = ((GetPos1()-GetPos2())*(GetVel1()-GetVel2()))/l-offset;
	}
};

void RigidLink::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	Vector hTrans;

	if(IsLocalNodeConstraint())
	{
		if (UsePenaltyFormulation()) return;

		double sign = 1;
		if (locelemind == 2) sign = -1;

		if (nodes(locelemind)!=0)						// position defined by local node number (also for CMS-Element)
		{
			GetBody3D(locelemind).GetNodedPosdqT(nodes(locelemind), dpdq);
		}
		else
		{
			GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
		}

		Mult(dpdq,ComputeForceDirection(),hTrans);

		for (int i=1; i <= f.Length(); i++)
		{
			f(i) -= sign*(hTrans(i)*XG(1));
		}
	}
	else		// global nodes
	{
			// nothing at all?
	}
};


int RigidLink::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// Automatic entries for this class 
	BasePointJoint::GetAvailableSpecialValues(available_variables);
	RigidLink::GetAvailableSpecialValuesAuto(available_variables);

	// Manual WRITE entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.link_length",0,0,0.,mystr("distance between the connected points (l0)"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.link_velocity",0,0,0.,mystr("derivative of the distance with respect to time (v)"), TRWElementDataWrite));

	return 0;
}

int RigidLink::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = BasePointJoint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write

	if(RWdata.variable_name.CStrCompare("Connector.link_length"))
	{
		IO_dist = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.link_velocity"))
	{
		IO_vel = RWdata.value;
		return 1; 
	}

	return WriteSingleElementDataAuto(RWdata);
}

//***************************************************************************************************************************************************
// RotatorySpringDamperActuator
//***************************************************************************************************************************************************

void RotatorySpringDamperActuator::CopyFrom(const Element& e)
{
	BaseBodyJoint::CopyFrom(e);
	const RotatorySpringDamperActuator& ce = (const RotatorySpringDamperActuator&)e;

	phi0 = ce.phi0;
	ma = ce.ma;
	glob_rot_axis = ce.glob_rot_axis;

	forcemode = ce.forcemode;

	mathfunc_k = ce.mathfunc_k;
	mathfunc_d = ce.mathfunc_d;

	moment_k = ce.moment_k;
	moment_d = ce.moment_d;

	spring_res = ce.spring_res;
}

void RotatorySpringDamperActuator::ElementDefaultConstructorInitialization()
{
	BaseBodyJoint::ElementDefaultConstructorInitialization();

	SetDampingCoeff(1);
	SetPenaltyStiffness(1e2);
	SetPenaltyFormulation(1);
	phi0 = 0;
	ma = 0;
	forcemode = 0;
	draw_dim = Vector3D(0.01,10,0.006); //drawdim: X=spring drawsize, Y=windings, Z=axis radius (cylinder)
	spring_res = 5;

	mathfunc_k = MathFunction();
	mathfunc_d = MathFunction();

	moment_k = 0;
	moment_d = 0;
}

void RotatorySpringDamperActuator::Initialize() 
{
	if(auto_comp_ground)
	{
		SetPos2ToGlobalCoord(GetPos1());
	}
	Vector3D axis = glob_rot_axis;
	axis.Normalize();
	Vector3D n1, n2;
	axis.SetNormalBasis(n1,n2);

	JA0i.Set(axis,n1,n2);

	Matrix3D A0i = GetRotMatBodyi();
	Matrix3D A0j = GetRotMatBodyj();

	JA0j = A0j.GetTp()*A0i*JA0i;
}

Vector3D RotatorySpringDamperActuator::ComputeMoment(double t)  const
{
	double phi = 0;
	double xp = 0;
	double m = 0;

	if (UsePenaltyFormulation())
	{
		Vector3D n1 = GetRotMati()*Vector3D(0.,1.,0.);
		Vector3D n2 = GetRotMatj()*Vector3D(0.,1.,0.);
		Vector3D rot1 = GetRotMati()*Vector3D(1.,0.,0.);
		Vector3D rot2 = GetRotMatj()*Vector3D(1.,0.,0.);
		phi = NormalizedVectorAngle(n1,n2); //phi should always be positive
		if (phi < 0) phi += 2.*MY_PI; //should not happen!
		if ((0.5*(rot1+rot2))*(n1.Cross(n2)) < 0) // defines the rotation direction
		{
			phi = -phi;
		}

		if(UseDamping() && forcemode != 2) { // 0..constant coefficient, 1..MathFunction
			Vector3D angvel1;
			Vector3D angvel2;
			if (elements(1) != 0)
			{
				angvel1 = GetBody3D(1).GetAngularVel(loccoords(1));
			}
			else
			{
				angvel1 = Vector3D(0.);
			}
			if (elements(2) != 0)
			{
				angvel2 = GetBody3D(2).GetAngularVel(loccoords(2));
			}
			else
			{
				angvel2 = Vector3D(0.);
			}
			xp = (angvel1-angvel2)*rot1;
		}

		if (forcemode == 0) // linear stiffness and damper
		{
			double x = (phi-phi0);
			m = -1*GetPenaltyStiffness()*x;	
			if(UseDamping()) { m += GetDampingCoeff() * xp; }
		}
		else if (forcemode == 1) // nonlinear stiffness and damping by MathFunction
		{
			double x = (phi-phi0);
			if (mathfunc_k.GetFuncMode() == TMFpiecewiselinear) {m = -1*mathfunc_k.Evaluate(x)*x;} //mathfunc_k ... stiffness depending on 'x'
			if (mathfunc_d.GetFuncMode() == TMFpiecewiselinear) {m += mathfunc_d.Evaluate(xp)*xp;} //mathfunc_k ... damping depending on 'xp'
		} else
		{ // nonlinear stiffness and damping
			m = moment_k + moment_d;
		}

		m += ma; // add constant actor force

		return m*rot1;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: RotatorySpringDamperActuator: ComputeMoment just implemented for PenaltyFormulation";
		return Vector3D(0.0);
	}
}

Vector3D RotatorySpringDamperActuator::ComputeForce(double t) const
{
	if (UsePenaltyFormulation())
	{
		return Vector3D(0.);
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: RotatorySpringDamperActuator: ComputeForce just implemented for PenaltyFormulation";
		return Vector3D(0.0);
	}
};

void RotatorySpringDamperActuator::DrawElement() 
{
	Constraint::DrawElement();

	double r = draw_dim.X(); //radius of torsional spring
	double rev = draw_dim.Y(); //number of revolutions for torsional spring
	double t = draw_dim.Z(); //radius of axis cylinder

	mbs->SetColor(GetCol());

	Vector3D p1;
	Vector3D p2, rot, n1, t1, t2;

	Matrix3D R1=GetRotMatiD();
	Matrix3D R2=GetRotMatjD();

	p1 = GetDrawPosition(1);

	rot = R1*Vector3D(1,0,0);
	n1 = R1*Vector3D(0,1,0);
	t1 = rot.Cross(n1);
	t2 = rot.Cross(R2*Vector3D(0,1,0));

	p2 = GetDrawPosition(2);

	double dphi = NormalizedVectorAngle(t2,n1);

	double phioff = 0;
	dphi=-MY_PI*0.5+dphi; //only 90°+90°
	if (NKinPairs() < 2) dphi += MY_PI; // is ground joint
	else 
	{
		phioff = dphi;
		dphi = -dphi + MY_PI;
	}

	if (r != 0) 
	{
		int mode = 2; //mode=1: zig-zag, mode=2: revolutions
		double ntile = 16*mode*spring_res; //tiling
		if (mode == 2) dphi += rev*MY_PI; //6 revolutions
		for (double i = 1; i <= ntile; i++)
		{
			double phia = (i-1.)/ntile*dphi+phioff;
			double phib = (i)/ntile*dphi+phioff;
			double ra = r;
			double rb = r;

			if (mode == 1)
			{
				if ((i-1) > 1 && (i-1) < ntile-1 && ((int)i-1)% 2 == 0) ra = r*0.85;
				if (i > 1 && i < ntile-1 && (int)i % 2 == 0) rb = r*0.85;
			}
			else
			{
				ra *= (1.-0.33*phia/dphi);
				rb *= (1.-0.33*phib/dphi);
			}
			Vector3D va = p1 + ra*cos(phia)*n1 + ra*sin(phia)*t1;
			Vector3D vb = p1 + rb*cos(phib)*n1 + rb*sin(phib)*t1;

			mbs->MyDrawLine(va,vb,2);
		}
	}

	mbs->DrawZyl(p1-(0.5*r)*rot,p1+(0.5*r)*rot,t,8*(int)spring_res);

	//double rs = 0.65*draw_dim.Z();//original
	double rs = 0.02*draw_dim.X();
	mbs->DrawSphere(p1,rs,12);
	mbs->SetColor(GetColExt());
	mbs->DrawSphere(p2,rs,12);
	
	/*
	//draw n1, t1, t2:
	mbs->SetColor(Vector3D(0,0,1)); // blue
	mbs->DrawZyl(p1-(0.5*r)*n1,p1+(0.5*r)*n1,t,8);
	mbs->SetColor(Vector3D(1,0,0)); // red
	mbs->DrawZyl(p1-(0.5*r)*t1,p1+(0.5*r)*t1,t,8);

	mbs->SetColor(Vector3D(0.8,0.8,0.2)); // braun
	mbs->DrawZyl(p1-(0.5*r)*t2,p1+(0.5*r)*t2,t,8);
	*/
};

void RotatorySpringDamperActuator::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));

	// not implemented yet for nodes
	
	// penalty + elemNr + locCoord
	kaf = (int)(TKinematicsAccessFunctions(TKAF_rotation_matrix+TKAF_angular_velocity+TKAF_D_rot_D_q));	
	KinAccFunc.Add(kaf);
}

void RotatorySpringDamperActuator::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;

	GetElementDataAuto(edc);

	// MathFunction
	ElementDataContainer edc_mf_k;
	mathfunc_k.GetElementData(edc_mf_k);
	ed.SetEDC(&edc_mf_k,"MathFunction_k"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the stiffness as function of phi-phi0."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	

	ElementDataContainer edc_mf_d;
	mathfunc_d.GetElementData(edc_mf_d);
	ed.SetEDC(&edc_mf_d,"MathFunction_d"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the damping as function of omega."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	
}

int RotatorySpringDamperActuator::SetElementData(ElementDataContainer& edc) 		//update element data
{
	int rv = 1;

	ElementData* mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_k"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *edcp_mf;
		mathfunc_k.SetElementData(mbs,edc_mf);
	}

	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_d"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *edcp_mf;
		mathfunc_d.SetElementData(mbs,edc_mf);
	}

	SetElementDataAuto(edc);
	return rv;
}

int RotatorySpringDamperActuator::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// Automatic entries for this class 
	BasePointJoint::GetAvailableSpecialValues(available_variables);
	RotatorySpringDamperActuator::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_acting_moment",0,0,0.,mystr("internal moment of connector"), TRWElementDataRead));

	available_variables.Add(ReadWriteElementDataVariableType("Connector.angle_offset",0,0,0.,mystr("prescribe the angle offset"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.spring_moment",0,0,0.,mystr("prescribe the stiffness moment"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.damper_moment",0,0,0.,mystr("prescribe the damping moment"), TRWElementDataWrite));

	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_global",3,0,0.,mystr("internal global force of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_global",3,0,0.,mystr("internal global moment of connector"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_force_local",3,0,0.,mystr("internal local force of connector (joint coordinate system JAi)"), TRWElementDataRead));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_moment_local",3,0,0.,mystr("internal local moment of connector (joint coordinate system JAi)"), TRWElementDataRead));

	//angle
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_angle",3,0,0.,mystr("bryant angles between the joint coordinate systems JAi and JAj. All constrained components are zero."), TRWElementDataRead));

	// displacement
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_displacement",3,0,0.,mystr("displacement between the joint coordinate systems JAi and JAj expressed in coordinate system JAi"), TRWElementDataRead));

	return 0;
}

int RotatorySpringDamperActuator::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = BaseBodyJoint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.constraint_acting_moment"))
	{
		RWdata.value = ComputeMoment(mbs->GetTime())*(GetRotMati()*Vector3D(1.,0.,0.));
		return 1;
	}

	return ReadSingleElementDataAuto(RWdata);
}

int RotatorySpringDamperActuator::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = BaseBodyJoint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write
	if(RWdata.variable_name.CStrCompare("Connector.angle_offset"))
	{
		phi0 = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.spring_moment"))
	{
		moment_k = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.damper_moment"))
	{
		moment_d = RWdata.value;
		return 1; 
	}
	return WriteSingleElementDataAuto(RWdata);
}

int RotatorySpringDamperActuator::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	return BasePointJoint::CheckConsistency(errorstr);
}


void SpringDamperActuator2D::ElementDefaultConstructorInitialization()
{
	Constraint::ElementDefaultConstructorInitialization();

	use_penalty_formulation = 1;

	SetPenaltyStiffness(1e2);

	damping_coeff = 1;

	l0 = 0;
	fa = 0;
	forcemode = 0;
	draw_dim = Vector3D(0.01,10,0.01);
	spring_res = 5;

	mathfunc_k = MathFunction();
	mathfunc_d = MathFunction();

	force_k = 0;
	force_d = 0;

	loccoords.Set2(0,0);
	elements.Set2(1,0);

	elementname = GetElementSpec();
}

void SpringDamperActuator2D::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const SpringDamperActuator2D& ce = (const SpringDamperActuator2D&)e;

	loccoords = ce.loccoords;
	damping_coeff = ce.damping_coeff;

	l0 = ce.l0;
	fa = ce.fa;

	forcemode = ce.forcemode;
	mathfunc_k = ce.mathfunc_k;
	mathfunc_d = ce.mathfunc_d;

	force_k = ce.force_k;
	force_d = ce.force_d;

	spring_res = ce.spring_res;
}

int SpringDamperActuator2D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);

	if(elements(1)>mbs->NE() || elements(2)>mbs->NE())
	{
		UO(UO_LVL_err) << "An element defined in the joint with element number " << GetOwnNum() << " does not exist!\n";
		rv = 1;
	}
	if(elements(1)==GetOwnNum() || elements(2)==GetOwnNum())
	{
		UO(UO_LVL_err) << "An element defined in the joint with element number " << GetOwnNum() << " has the same number as the element number of the joint!\n";
		rv = 1;
	}

	if (rv) return rv;
}

void SpringDamperActuator2D::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair) // MSax 2013-08-05: added
{
	KinAccFunc.SetLen(0);

	int kaf = (int)(TKinematicsAccessFunctions(TKAF_position_2D+TKAF_velocity_2D+TKAF_D_pos_D_q_2D));	
	KinAccFunc.Add(kaf);
}

int SpringDamperActuator2D::SOS() const 
	{
		int nsos = 0;

		for (int i=1; i <= elements.Length(); i++)
		{
			if(elements(i)!=0)
			{
				if(elements(i)<=mbs->NE() && elements(i)>0 && elements(i) != GetOwnNum())
				{
					nsos += GetElem(i).SOS();
				}
				else
				{
					return 0; // in this case the elements are not correctly defined yet
				}
			}
		}
		return nsos;
	};  // explicit size, number of constrained dofs

void SpringDamperActuator2D::EvalF2(Vector& f, double t)
{
	double sign = 1.;
	int offset = 0;

	if(!UsePenaltyFormulation()) return;	// no penalty method --> Lagrange multiplier --> no EvalF2
	
	Vector2D force = ComputeForce(t);

	for (int i=1; i <= elements.Length(); i++)
	{
		if (i==2) sign = -1.;
		if(elements(i)!=0)
		{
			GetBody2D(i).GetdPosdqT(loccoords(i),dpdq);

			if (i==2) offset = GetElem(1).SOS();
			for (int j=1; j<=GetElem(i).SOS(); j++)
			{
				f(j+offset) -= sign*(dpdq(j,1)*force.X() + dpdq(j,2)*force.Y());
			}
		}
	}
};

Vector2D SpringDamperActuator2D::ComputeForce(double t) const
{
	double l = 0;
	double xp = 0;
	Vector2D dir = ComputeForceDirection();
	double f = 0;

	if (UsePenaltyFormulation())
	{
		Vector2D v12 = GetPosition(1) - GetPosition(2);
		l = v12*dir;
		if(UseDamping()) {	
			Vector2D v12p = GetVelocity(1) - GetVelocity(2);
			xp = v12p*dir;
			//xp = (v12*v12p)/l; // MSax: correct but not efficient
		}

		if (forcemode == 0) // linear stiffness and damper
		{
			double x = (l-l0);
			f = GetPenaltyStiffness()*x;	
			if(UseDamping()) { f += GetDampingCoeff() * xp; }
		}
		else if (forcemode == 1) // nonlinear stiffness and damping by MathFunction
		{
			double x = (l-l0);
			if (mathfunc_k.GetFuncMode() == TMFpiecewiselinear) {f = mathfunc_k.Evaluate(x)*x;} //mathfunc_k ... stiffness depending on 'x'
			if (mathfunc_d.GetFuncMode() == TMFpiecewiselinear) {f += mathfunc_d.Evaluate(xp)*xp;} //mathfunc_k ... damping depending on 'xp'
		} else
		{ // nonlinear stiffness and damping
			f =	force_k + force_d;
		}

		f += fa; // add constant actor force

		return f*dir;
	}
	else
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: SpringDamperActuator: ComputeForce just implemented for PenaltyFormulation";
		return Vector2D(0.0);
	}
}

void SpringDamperActuator2D::SetStiffnessDamping()
{
	Matrix data;
	int pos = 0;
	if (pos) mathfunc_k.SetData(TMFpiecewiselinear, data);
	if (pos) mathfunc_d.SetData(TMFpiecewiselinear, data);
}

void SpringDamperActuator2D::SetPiecewiseLinearSpring(const Vector& xpos, const Vector& force)
{
	forcemode = 1;
	mathfunc_k.SetPiecewise(xpos, force, 1);
}
void SpringDamperActuator2D::SetPiecewiseLinearDamper(const Vector& vel, const Vector& force)
{
	forcemode = 1;
	mathfunc_d.SetPiecewise(vel, force, 1);
}

void SpringDamperActuator2D::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	
	Vector3D pa = GetDrawPosition(1);
	Vector3D pb = GetDrawPosition(2);

	if (elements.Length()==1) Swap(pa,pb);

	if (draw_dim.Y() == 0) // if coils == 0 ==> draw cylinder with spring diameter
	{
		mbs->DrawZyl(pa, pb, 0.5*draw_dim.X(),12);
	}
	else
	{
		//draw spring:
		double phi,phi2,xx;
		Vector3D vf;
		double steps = 200*0.5*spring_res; //0.2 ->> very coarse
		double rots = draw_dim.Y();
		double r = draw_dim.X()/2.;
		double rz = draw_dim.X()/10.;
		int ztile = 6;

		vf = pa-pb;
		double lenvf = vf.Norm();
		Vector3D dir = vf;
		vf = ComputeForceDirectionD();
		vf.Normalize();
		vf *= lenvf;
		if (vf*dir < 0) vf *= -1;

		Vector3D vfn = vf;
		Vector3D vz,vy;
		vfn.SetNormalBasis(vz,vy);
		Vector3D pz1, pz2;

		if(GetPenaltyStiffness() || (forcemode > 0)) //$ DR 2013-05-03 linear spring damper with stiffness = 0, do not draw spring
		{
			double off = steps*0.1;
			mbs->SetColor(GetCol());
			for (xx = off; xx <= steps-off; xx++)
			{
				phi = xx*2*MY_PI/steps*rots;
				phi2 = (xx+1.2)*2*MY_PI/steps*rots;

				Vector3D pr1=cos(phi)*r*vz+sin(phi)*r*vy;
				Vector3D pr2=cos(phi2)*r*vz+sin(phi2)*r*vy;

				pz1=pb+xx/steps*vf+pr1;
				pz2=pb+(xx+1)/steps*vf+pr2;
				mbs->DrawZyl(pz1,pz2,rz,ztile);
				if (xx==off) 
				{
					mbs->DrawZyl(pb,pz1,rz,ztile);
					mbs->DrawSphere(pz1,1.05*rz,12);
				}
			}
			mbs->DrawZyl(pa,pz2,rz,ztile);
			mbs->DrawSphere(pz2,1.05*rz,12);
		}

		if (draw_dim.Z() != 0 && GetDampingCoeff() != 0)
		{
			mbs->SetColor(GetColExt());
			mbs->DrawZyl(pa, pa-0.5*l0*vfn, 0.5*draw_dim.Z(),12);
			mbs->DrawZyl(pa-0.5*l0*vfn, pb, 0.5*draw_dim.Z()*0.6,12);

			mbs->DrawSphere(pa,0.65*draw_dim.Z(),12);
			mbs->DrawSphere(pb,0.65*draw_dim.Z(),12);
		}
	}
};

Vector2D SpringDamperActuator2D::ComputeForceDirection() const   //return global force direction, normalized
{
	Vector2D v12;
	v12 = GetPosition(1) - GetPosition(2);
	v12.Normalize();
	
	return v12;
}

Vector3D SpringDamperActuator2D::ComputeForceDirectionD()   //return global force direction, normalized
{
	Vector3D v12;
	v12 = GetDrawPosition(1)-GetDrawPosition(2);
	v12.Normalize();

	return v12;
}

void SpringDamperActuator2D::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	ElementData ed;
	GetElementDataAuto(edc);

	// MathFunction
	ElementDataContainer edc_mf_k;
	mathfunc_k.GetElementData(edc_mf_k);
	ed.SetEDC(&edc_mf_k,"MathFunction_k"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the stiffness as function of l-l0."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	

	ElementDataContainer edc_mf_d;
	mathfunc_d.GetElementData(edc_mf_d);
	ed.SetEDC(&edc_mf_d,"MathFunction_d"); ed.SetToolTipText("MathFunction is used, if forcemode = 1. Gives the damping as function of v."); edc.TreeAdd("Physics.MathFunction",ed);// edc.Add(ed);	
}

int SpringDamperActuator2D::SetElementData(ElementDataContainer& edc) 		//update element data
{
	int rv = 1;

	ElementData* mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_k"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *edcp_mf;
		mathfunc_k.SetElementData(mbs,edc_mf);
	}
	mathFuncEd = edc.TreeFind(mystr("Physics.MathFunction.MathFunction_d"));
	if(mathFuncEd != NULL && mathFuncEd->IsEDC())
	{
		ElementDataContainer* edcp_mf = mathFuncEd->GetEDC();
		ElementDataContainer edc_mf = *edcp_mf;
		mathfunc_d.SetElementData(mbs,edc_mf);
	}

	SetElementDataAuto(edc);
	return rv;
}

int SpringDamperActuator2D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// Automatic entries for this class 
	Constraint::GetAvailableSpecialValues(available_variables);
	SpringDamperActuator2D::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Connector.constraint_acting_force",0,0,0.,mystr("internal resultant force of connector"), TRWElementDataRead));

	available_variables.Add(ReadWriteElementDataVariableType("Connector.spring_length_offset",0,0,0.,mystr("prescribe the neutral spring length"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.spring_force",0,0,0.,mystr("prescribe the stiffness force"), TRWElementDataWrite));
	available_variables.Add(ReadWriteElementDataVariableType("Connector.damper_force",0,0,0.,mystr("prescribe the damping force"), TRWElementDataWrite));

	return 0;
}

int SpringDamperActuator2D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Connector.constraint_acting_force"))
	{
		RWdata.value = ComputeForce(mbs->GetTime())*ComputeForceDirection();
		return 1;
	}

	return ReadSingleElementDataAuto(RWdata);
}

int SpringDamperActuator2D::WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine ( not required in Element )	
	int rv = Constraint::WriteSingleElementData(RWdata);
	if (rv == 1) return 1;
	//manual things to write
	if(RWdata.variable_name.CStrCompare("Connector.spring_length_offset"))
	{
		l0 = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.spring_force"))
	{
		force_k = RWdata.value;
		return 1; 
	}
	if(RWdata.variable_name.CStrCompare("Connector.damper_force"))
	{
		force_d = RWdata.value;
		return 1; 
	}

	return WriteSingleElementDataAuto(RWdata);
}


Vector2D SpringDamperActuator2D::GetPosition(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return loccoords(2);
	}
	else
	{
		return GetBody2D(i).GetPos2D(loccoords(i));
	}
}

Vector2D SpringDamperActuator2D::GetVelocity(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return Vector2D(0.);
	}
	else
	{
		return GetBody2D(i).GetVel2D(loccoords(i));
	}
}

Vector3D SpringDamperActuator2D::GetDrawPosition(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return Vector3D(loccoords(2).X(),loccoords(2).Y(),0.);
	}
	else
	{
		return GetBody2D(i).ToP3D(GetBody2D(i).GetPos2DD(loccoords(i)));
		//return GetBody2D(i).GetPosD(loccoords(i));
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PointJoint2D

void PointJoint2D::ElementDefaultConstructorInitialization()
{
	Constraint::ElementDefaultConstructorInitialization();

	loccoords.Set2(0,0);
	elements.Set2(1,0);

	use_penalty_formulation = 0;

	SetPenaltyStiffness(1e8);
	damping_coeff = 1;

	spring_stiffness2 = Vector2D(0.);
	
	dir.Set2(1,1);

	draw_local_frame_size = 0.01;

	stiffness_in_joint_local_frame = 0;
	phi_z = 0;

	x_init = Vector(SS());
}

void PointJoint2D::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const PointJoint2D& ce = (const PointJoint2D&)e;

	loccoords = ce.loccoords;
	dir = ce.dir;
	draw_local_frame_size = ce.draw_local_frame_size;
	damping_coeff = ce.damping_coeff;
	spring_stiffness2 = ce.spring_stiffness2;

	stiffness_in_joint_local_frame = ce.stiffness_in_joint_local_frame;
	phi_z = ce.phi_z;
}

int PointJoint2D::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);

	// todo: check for same initial positions

	if (rv) return rv;
}

void PointJoint2D::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);

	int kaf = (int)(TKinematicsAccessFunctions(TKAF_position_2D+TKAF_velocity_2D+TKAF_D_pos_D_q_2D));	
	KinAccFunc.Add(kaf);
}

int PointJoint2D::SOS() const 
{
	if (!UsePenaltyFormulation()) 
	{
		return 0;
	}
	else
	{
		int nsos = 0;
		for (int i=1; i <= elements.Length(); i++)
		{
			if(elements(i)!=0)
			{
				if(elements(i)<=mbs->NE() && elements(i)>0 && elements(i) != GetOwnNum())
				{
					nsos += GetElem(i).SOS();
				}
				else
				{
					return 0; // in this case the elements are not correctly defined yet
				}
			}
		}
		return nsos;
	}
};  // explicit size, number of constrained dofs

Vector2D PointJoint2D::GetPenaltyStiffness2(double t) const 
{
	if(NSteps()==0)
		return spring_stiffness2;
	else
		return spring_stiffness2 * Constraint::GetCStepsFact(t);
}

void PointJoint2D::EvalF2(Vector& f, double t)
{
	double sign = 1.;
	int offset = 0;

	if(!UsePenaltyFormulation()) return;	// no penalty method --> Lagrange multiplier --> no EvalF2
	
	Vector2D force = ComputeForce(t);

	for (int i=1; i <= elements.Length(); i++)
	{
		if (i==2) sign = -1.;
		if(elements(i)!=0)
		{
			GetBody2D(i).GetdPosdqT(loccoords(i),dpdq);

			if (i==2) offset = GetElem(1).SOS();
			for (int j=1; j<=GetElem(i).SOS(); j++)
			{
				f(j+offset) -= sign*(dpdq(j,1)*force.X() + dpdq(j,2)*force.Y());
			}
		}
	}
};

Vector2D PointJoint2D::ComputeForce(double t) const
{
	Vector2D u; // displacement
	Vector2D v; // velocity
	Vector2D f; // resultant force
	Vector2D k = this->GetPenaltyStiffness2(t); // spring stiffness
	Matrix3D rot; // reference rotation
	Matrix3D k_mat;

	if (UsePenaltyFormulation())
	{
		u = GetPosition(1) - GetPosition(2);
		if(UseDamping()) {	v = GetVelocity(1) - GetVelocity(2);}

		k_mat.SetSize(2,2); 
		k_mat(1,1) = k.X();
		k_mat(2,2) = k.Y();

		if(UseLocalCoordinateSystem()||stiffness_in_joint_local_frame)
		{
			rot = GetRotMati();
			rot.Transpose();
			k_mat = k_mat*rot;
			rot = GetRotMati();
			k_mat = rot*k_mat;
		}

		f = k_mat*u;	

		if(UseDamping()) { f += GetDampingCoeff() * v; }

		return f;
	}
	else
	{
		Vector2D force(XG(1)*dir(1), XG(2)*dir(2));
		return force;
	}
}

void PointJoint2D::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG
	Vector2D f2;

	if (!UseLocalCoordinateSystem()) // all directions constrained
	{
		if (MaxIndex()==3)
		{
			f2 = GetPosition(1)-GetPosition(2);
		}
		else	// MaxIndex < 3
		{ //velocity constraints:
			f2 = GetVelocity(1)-GetVelocity(2);
		}
	}
	else // local coordinate system is used
	{
		if (MaxIndex()==3)
		{
			f2 = GetRotMati().GetTp()*(GetPosition(1)-GetPosition(2));
		}
		else 	// MaxIndex < 3
		{ //velocity constraints:
			f2 = GetRotMatiP().GetTp()*(GetPosition(1)-GetPosition(2))+GetRotMati().GetTp()*(GetVelocity(1)-GetVelocity(2));
		}
	}

	int next=0;
	for(int i=1; i<=2; i++) 
	{
		if(dir(i))
		{
			f(++next)=f2(i);
		}
	}
};

//$ SW 2013-10-10: trying to avoid matrix and vector allocations
void PointJoint2D::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	if (UsePenaltyFormulation()) return;

	double sign = 1;
	if (locelemind == 2) sign = -1;
	//Matrix hTrans;

	if (!UseLocalCoordinateSystem()) // and all directions are constrained (else EvalG writes an error)
	{
		GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	}
	else //UseLocalCoordinateSystem()
	{
		GetBody2D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);

		if (dir(1) && dir(2)) // all directions constraint
		{
			Matrix3D Ai = GetRotMati();
			//$ SW 2013-10-9: avoid new allocation of matrices to increase the performance
			//Matrix Ai_;
			//Ai_.SetSize(2,2);
			//Ai_(1,1) = Ai(1,1); Ai_(1,2) = Ai(1,2); Ai_(2,1) = Ai(2,1); Ai_(2,2) = Ai(2,2);
			Matrix2D Ai_;
			Ai_(1,1) = Ai(1,1); Ai_(1,2) = Ai(1,2); Ai_(2,1) = Ai(2,1); Ai_(2,2) = Ai(2,2);

			//hTrans = dpdq*Ai_;  => dpdq *=Ai_
			dpdq.ApplySqrMat(Ai_);
		}
		else // not all directions constraint
		{	
			if (locelemind==1)
			{
				GetdRotTvdqT(GetPosition(1)-GetPosition(2), loccoords(locelemind), drotTvdq, locelemind);
				//Matrix AiB_;
				//AiB_= GetBody2D(locelemind).GetRotMatrix2D();
				Matrix3D AiB_ = GetBody2D(locelemind).GetRotMatrix2D();
				//Matrix JA0i_;
				//Matrix3D rotPhi;
				//rotPhi.SetSize(2,2);
				Matrix2D rotPhi;
				if (stiffness_in_joint_local_frame)
				{
					double cosphi = cos(phi_z); 
					double sinphi = sin(phi_z); 
					rotPhi(1,1) = cosphi;
					rotPhi(1,2) =-sinphi;
					rotPhi(2,1) = sinphi;
					rotPhi(2,2) = cosphi;
					//JA0i_ = rotPhi;
				}
				else
				{
					rotPhi(1,1) = 1; rotPhi(1,2) = 0; rotPhi(2,1) = 0; rotPhi(2,2) = 1;
					//JA0i_ = rotPhi; // not in joint local frame
				}
				//JA0i_.SetSize(2,2);
				AiB_.SetSize(2,2);

				//$ SW 2013-10-9: this would allocate three new matrices during computation
				//hTrans = (drotTvdq+dpdq*AiB_)*JA0i_;
				dpdq.ApplySqrMat(AiB_);
				dpdq += drotTvdq;
				dpdq.ApplySqrMat(rotPhi);
			}
			else
			{
				Matrix3D Ai = GetRotMati();
				//Matrix Ai_;
				//Ai_.SetSize(2,2);
				//Ai_(1,1) = Ai(1,1); Ai_(1,2) = Ai(1,2); Ai_(2,1) = Ai(2,1); Ai_(2,2) = Ai(2,2);
				//hTrans = dpdq*Ai_;
				Matrix2D Ai_(
					Ai(1,1),Ai(1,2),
					Ai(2,1),Ai(2,2));
				dpdq.ApplySqrMat(Ai_);
			}
		}
	}

	for (int i=1; i <= f.Length(); i++)
	{
		int next = 0;
		for(int c=1; c<=2; c++)
		{
			if(dir(c)==1.)
			{
				//f(i) -= sign*(hTrans(i,c)*XG(++next));
				f(i) -= sign*(dpdq(i,c)*XG(++next));
			}
		}
	}
};

void PointJoint2D::GetdRotTvdqT(const Vector2D& vloc, const Vector2D& ploc, Matrix& d, int bodyindex)
{
	// MSax: computes the dRotTvdqT matrix with dRotvdqT element methods. See BasePointJoint for more information
	//$ SW 2013-10-10: uses a class member instead of allocating a new matrix
	Matrix &coeffMat = tmpmat;
	GetBody2D(bodyindex).GetdRotvdqT(Vector2D(1.,0.),ploc,coeffMat); //coefficient matrix 1
	
	d.SetSize(GetBody2D(bodyindex).SOS(),Dim());
	d.SetAll(0);

	for (int i=1; i<=SOS(); i++)
	{
		for (int j=1; j<=2; j++)
		{
			d(i,1) += coeffMat(i,j)*vloc(j);
		}
	}

	GetBody2D(bodyindex).GetdRotvdqT(Vector2D(0.,1.),ploc,coeffMat); //coefficient matrix 2

	for (int i=1; i<=SOS(); i++)
	{
		for (int j=1; j<=2; j++)
		{
			d(i,2) += coeffMat(i,j)*vloc(j);
		}
	}
};

void PointJoint2D::DrawElement() 
{
	Constraint::DrawElement();
	double t_draw = Constraint::GetGlobalTime();

	if (GetDrawSizeScalar() == 0) return;

	mbs->SetColor(GetCol());

	TArray<Vector3D> constr_dir;
	Vector3D tmp_dir;
	Vector3D p;
	Vector3D k;
	int flag = 1;

	for (int i=1; i <= NE(); i++)
	{
		p = GetDrawPosition(i);
		constr_dir.Flush();

		if ( i== 2) flag = -1;
		if(!mbs->GetSimulationStatus().GetStatusFlag(TSimulationRunning) && t_draw == 0.) k = Vector3D(1.0); //$ AD 2011-09-15: draw all constraints before computation starts (even if they are off at t=0); //$ MaSch 2013-08-19
		else if(UsePenaltyFormulation())
		{
			Vector2D k2D = GetPenaltyStiffness2(t_draw);
			k = Vector3D(k2D.X(),k2D.Y(),1.);
		}
		else k = Vector3D(dir(1),dir(2),1.);		

		if (UseLocalCoordinateSystem())
		{
			tmp_dir.Set(0.0,0.0,0.0);
			Vector3D tmp_dir2;
			Matrix3D k_mat(k(1),k(2),k(3));
			k_mat = GetRotMatiD() * k_mat;

			for(int j=1; j<= 3; j++)
			{
				tmp_dir.Set(0.0,0.0,0.0);
				tmp_dir(j) = flag;
				tmp_dir2 = k_mat * tmp_dir;

				if(!(tmp_dir2 == Vector3D(0.0)))
				{
					tmp_dir2.Normalize();
					constr_dir.Add(tmp_dir2);
				}
			}
		}
		else
		{
			for(int j=1; j<= 3; j++)
			{
				tmp_dir.Set(0.0,0.0,0.0);
				if(k(j)) 
				{
					tmp_dir(j) = flag;
					constr_dir.Add(tmp_dir);
				}
			}
		}

		if(constr_dir.Length()==3)
		{
			mbs->DrawSphere(p, GetDrawSizeScalar());
		}
		else	// just some directions are constrained
		{
			for (int j=1; j<=constr_dir.Length(); j++)
			{
				mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*constr_dir(j),p,(GetDrawSizeScalar()*0.6),6,1);
			}
		}
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
		Vector3D p = GetDrawPosition(1);

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

int PointJoint2D::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// Automatic entries for this class 
	Constraint::GetAvailableSpecialValues(available_variables);
	PointJoint2D::GetAvailableSpecialValuesAuto(available_variables);

	// Manual READ entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.acting_force",0,0,0.,mystr("internal resultant force of connector"), TRWElementDataRead));

	return 0;
}

int PointJoint2D::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	//manual things to read  
	if(RWdata.variable_name.CStrCompare("Internal.acting_force"))
	{
		RWdata.value = ComputeForce(mbs->GetTime()).Norm();
		return 1;
	}

	return ReadSingleElementDataAuto(RWdata);
}


Matrix3D PointJoint2D::GetRotMati() const
{
	Matrix3D A(1);

	if(UseLocalCoordinateSystem())
	{
		A = GetBody2D(UseLocalCoordinateSystem()).GetRotMatrix2D();  // $ MSax 2013-08-08: in some joints it is possible to set UseLocalCoordinateSystem to an integer > 1, then this body is used!
	}

	// additional (constant) rotation 
	if(stiffness_in_joint_local_frame)
	{
		Matrix3D rotPhi;
    rotPhi.SetSize(2,2);
		double cosphi = cos(phi_z); 
		double sinphi = sin(phi_z); 
		rotPhi(1,1) = cosphi;
		rotPhi(1,2) =-sinphi;
		rotPhi(2,1) = sinphi;
		rotPhi(2,2) = cosphi;
		A = A*rotPhi;
	}
	return A;
}

Matrix3D PointJoint2D::GetRotMatiP() const
{
	Matrix3D A(0);

	if(UseLocalCoordinateSystem())
	{
		A = GetBody2D(UseLocalCoordinateSystem()).GetRotMatrix2DP();
	}

	if(stiffness_in_joint_local_frame)
	{
		Matrix3D rotPhi;
    rotPhi.SetSize(2,2);
		double cosphi = cos(phi_z); 
		double sinphi = sin(phi_z); 
		rotPhi(1,1) = cosphi;
		rotPhi(1,2) =-sinphi;
		rotPhi(2,1) = sinphi;
		rotPhi(2,2) = cosphi;
		A = A*rotPhi;
	}
	return A;
}

Matrix3D PointJoint2D::GetRotMatiD() const
{
	Matrix3D A(1);

	if(UseLocalCoordinateSystem())
	{
		A = GetBody2D(UseLocalCoordinateSystem()).GetRotMatrixD();
	}

	if(stiffness_in_joint_local_frame)
	{
		A = A*ComputeRotMatrixWithKardanAngles(0,0,phi_z);
	}
	return A;
}

Vector2D PointJoint2D::GetPosition(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return loccoords(2);
	}
	else
	{
		return GetBody2D(i).GetPos2D(loccoords(i));
	}
}

Vector2D PointJoint2D::GetVelocity(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return Vector2D(0.);
	}
	else
	{
		return GetBody2D(i).GetVel2D(loccoords(i));
	}
}

Vector3D PointJoint2D::GetDrawPosition(int i) const
{
	if (!elements(2) && i==2) // ground joint
	{
		return Vector3D(loccoords(2).X(),loccoords(2).Y(),0.);
	}
	else
	{
		return GetBody2D(i).ToP3D(GetBody2D(i).GetPos2DD(loccoords(i)));
	}
}
