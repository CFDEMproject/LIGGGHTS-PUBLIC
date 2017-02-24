//#**************************************************************
//#
//# filename:             KinematicPairs.cpp
//#
//# author:               Gerstmayr Johannes, Reischl Daniel
//#
//# generated:						20.April 2011
//# description:          Kinematic pairs can be divided in lower (surface contact) and higher (point or line contact) kinematic pairs. 
//#												In this file the lower kinematic pairs and the rigid joints (all d.o.f. are constrained) are provided.
//#												Other implemented constraints are in SpecialConstraints.h and some old kinematic pairs can be found in
//#												KinematicPairsDeprecated.h.
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
 
//# Rigid:				all d.o.f. are constrained, e.g. welding

//# lower kinematic pairs:

//# Spherical:		all translatory d.o.f. are constrained
//# Rotary:				all translatory and 2 rotatory d.o.f. are constrained, only rotation about 1 axis possible
//# Cylindrical:	like Rotary, but additionally translation along rotational-axis possible
//# Translatory:	all rotary and 2 translatory d.o.f are constrained
//# Plane:				not implemented
//# Skrew-type:		not implemented

#include "element.h"
#include "constraint.h"
#include "kinematicpairs.h"
#include "elementdataaccess.h"
#include "graphicsconstants.h"
#include "rendercontext.h"

#include "femathhelperfunctions.h"

//$ SW 2013-08-29: Moved CylindricalJoint, PrismaticJoint, RevoluteJoint and RigidJoint to KinematicPairsDeprecated.cpp
// and renamed them into CylindricalJointOLD, PrismaticJointOLD, RevoluteJointOLD and RigidJointOLD respectively.
// The classes CylindricalJoint_, PrismaticJoint_, RevoluteJoint_ and RigidJoint_ have been renamed to CylindricalJoint,
// PrismaticJoint, RevoluteJoint and RigidJoint (these classes can be found in RigidBodyJoints.h)


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical	Spherical	Spherical	Spherical	Spherical	Spherical	Spherical	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: NodalConstraint	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void NodalConstraint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	//ed.SetInt((int)draw_dim.Z(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 16"); edc.Add(ed);

	int i;
	for (i=1; i <= NConstrainedNodes(); i++)
	{
		ed.SetInt(NodeNum(i), mystr("NodeNum")+mystr(i)); ed.SetToolTipText("Number of node"); edc.Add(ed);
	}

	for (i=1; i <= loccoords.Length(); i++)
	{
		ed.SetInt(loccoords(i), mystr("Node_Coordinate")+mystr(i)); ed.SetToolTipText("Local coordinate of node to be constrained"); edc.Add(ed);
	}

	if (UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if(GetPenalty())
	{
		ed.SetDouble(GetPenaltyStiffness(), "Penalty_stiffness"); ed.SetToolTipText("Penalty stiffness of coordinate constraint"); edc.Add(ed);
		ed.SetDouble(damping_coeff, "Penalty_damping"); ed.SetToolTipText("Penalty damping of coordinate constraint"); edc.Add(ed);
	}
}

int NodalConstraint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);
	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0);
	SetDrawSizeScalar(tmp);
	//int dd;
	//if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) draw_dim.Z() = dd;

	int i;
	for (i=1; i <= NConstrainedNodes(); i++)
	{
		GetElemDataInt(mbs, edc, mystr("NodeNum")+mystr(i), NodeNum(i), 1);
	}

	for (i=1; i <= loccoords.Length(); i++)
	{
		GetElemDataInt(mbs, edc, mystr("Node_Coordinate")+mystr(i), loccoords(i), 1);
		if (loccoords(i) < 1 )
		{
			loccoords(i) = 1;
			GetMBS()->EDCError("Illegal element coordinate in NodalConstraint");
			rv = 0;
		}
	}

	if (UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if(GetPenalty())
	{
		GetElemDataDouble(mbs, edc, "Penalty_stiffness", tmp);
		SetPenaltyStiffness(tmp);
	  GetElemDataDouble(mbs, edc, "Penalty_damping", damping_coeff);
	}

	return rv;
}

void NodalConstraint::EvalG(Vector& f, double t) 
{
	if (UsePenaltyFormulation()) return; //(RL) 	if (UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if(penalty)

//<<<<<<< KinematicPairs.cpp
	if((NSteps()!=0)&&(Constraint::GetCStepsFact(t) == 0.))
//=======
//	if(NSteps() != 0 && Constraint::GetCStepsFact(t) == 0.)
//>>>>>>> 1.13
	{
		f(1) = XG(1);
		return; //$ AD 2011-09-09: hard turnoff - disable when constraint-factor is 0.
	}

	if (IsLocalNodeConstraint())
	{
		f(1) = 0;
		//only for 3D nodes!
		if (Dim() != 3) GetMBS()->UO().InstantMessageText("ERROR: NodalConstraint::EvalG only implemented for 3D Nodes");
		
		for (int i=1; i <= NConstrainedNodes(); i++)
		{
			double flag = 1;
			if (i==2) flag = -1;

			Vector3D nodepos;
			Vector3D nodevel;

			nodepos = GetElem(i).GetNodePos(nodes(i));
			nodevel = GetElem(i).GetNodeVel(nodes(i));

			//UO() << "i=" << i << ": ********\n";
			//UO() << "evalg::nodepos=" << nodepos << "\n";
			//UO() << "evalg::nodevel=" << nodevel << "\n";

			if (MaxIndex()==3)
			{
				if (!IsVelocityConstraint())
				{
					f(1) += -flag * nodepos(loccoords(i));
				}
				else
				{ //velocity constraints:
					f(1) += -flag * nodevel(loccoords(i));
				}
			}
			else
			{ //velocity constraints:
				if (!IsVelocityConstraint())
				{
					f(1) += -flag * nodevel(loccoords(i));
				}
				else
				{ //velocity constraints:
					f(1) += -flag * nodevel(loccoords(i));
				}
			}
		}
	}
	else
	{
		//only for 1 node!
		if (NConstrainedNodes() != 1) GetMBS()->UO().InstantMessageText("ERROR: NodalConstraint::EvalG only implemented for ground nodes");

		if (MaxIndex()==3)
		{
			if (!IsVelocityConstraint())
			{
				for (int i=1; i<=loccoords.Length(); i++)
					f(i) = -(XG(i)-ground(i)->Evaluate(t));
			}
			else
			{ //velocity constraints:
				for (int i=1; i<=loccoords.Length(); i++)
					f(i) = -(XGP(i)-ground(i)->Evaluate(t));
			}
		}
		else
		{ //velocity constraints:
			if (!IsVelocityConstraint())
			{
				for (int i=1; i<=loccoords.Length(); i++)
					f(i) = -XGP(i);
			}
			else
			{ //velocity constraints:
				for (int i=1; i<=loccoords.Length(); i++)
					f(i) = -(XGP(i)-ground(i)->Evaluate(t));
			}
		}
	}

};


void NodalConstraint::EvalF2(Vector& f, double t)
{
	//if (!mbs->IsJacobianComputation())
	//{
	//	Vector3D nodepos(XG(1),XG(2),XG(3));
	//	Vector3D nodevel(XGP(1), XGP(2), XGP(3));
	//	global_uo << "time " << t << ", pos = " << nodepos << ", vel = " << nodevel << "\n";
	//}

	if (IsLocalNodeConstraint())
	{
		if (Dim() != 3) GetMBS()->UO().InstantMessageText("ERROR: NodalConstraint::EvalF2 only implemented for 3D Nodes");

		Vector3D pos(0.);
		Vector3D vel(0.);
		for (int i=1; i <= NConstrainedNodes(); i++)
		{
			double flag = 1;
			if (i==2) {	flag = -1; }
			Vector3D nodepos = GetBody3D(i).GetNodePos(NodeNum(i));
			Vector3D nodevel = GetBody3D(i).GetNodeVel(NodeNum(i));
			pos += flag*nodepos;
			vel += flag*nodevel;
		}

		for (int i=1; i <= NConstrainedNodes(); i++)
		{
			int sosoff = 0;
			GetBody3D(i).GetNodedPosdqT(NodeNum(i), dpdq);

			double flag = 1;
			if (i==2) 
			{
				flag = -1;
				sosoff = GetBody3D(1).SOS();
			}

			if (!UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if (!GetPenalty())
			{
				double lambda = XG(2*SOS()+1);
				//UO() << "lambda=" << lambda << "\n";

				for (int j=1; j<=dpdq.Getrows(); j++)
				{
					f(j + sosoff) -= flag * dpdq(j,loccoords(i))*lambda;
				}
			}
			else
			{
				for (int j=1; j<=dpdq.Getrows(); j++)
				{
					f(j + sosoff) -= flag * GetPenaltyStiffness(t)*pos(loccoords(i))*dpdq(j,loccoords(i));
					f(j + sosoff) -= flag * damping_coeff*vel(loccoords(i))*dpdq(j,loccoords(i));

					if(loccoords.Length()>3){ //$ DR 2011-03-09:[ if you are using IVector to define more than 1 direction
						f(j + sosoff) -= flag * GetPenaltyStiffness(t)*pos(loccoords(i+2))*dpdq(j,loccoords(i+2));
						f(j + sosoff) -= flag * damping_coeff*vel(loccoords(i+2))*dpdq(j,loccoords(i+2));
						if(loccoords.Length()>5){
							f(j + sosoff) -= flag * GetPenaltyStiffness(t)*pos(loccoords(i+4))*dpdq(j,loccoords(i+4));
							f(j + sosoff) -= flag * damping_coeff*vel(loccoords(i+4))*dpdq(j,loccoords(i+4));
						}
					} //$ DR 2011-03-09:]
				}
			}
		}
	}
	else
	{
		if (!UsePenaltyFormulation()) //$ DR 2011-04-20: old code: if (!GetPenalty())
		{
			// Term Cq^T lambda
			for (int i=1; i<=loccoords.Length(); i++)
				f(i) -= XG(2*SOS()+i);
		}
		else  // penalty: spring+damping
		{		
			for (int i=1; i<=loccoords.Length(); i++)
			{
				f(i) -= GetPenaltyStiffness(t)*(XG(i)-ground(i)->Evaluate(t));
				f(i) -= damping_coeff*XGP(i);
			}
		}
	}
}

//$ AD 2011-09-08: intoduce computation steos for constraints:   NodalConstraint && penalty formulation
double NodalConstraint::GetPenaltyStiffness(double t)
{
	if(t < 0.)
		return Constraint::GetPenaltyStiffness();                // case no t passed for evaluation (NodalConstraintCMS)

  if(NSteps()==0)
		return Constraint::GetPenaltyStiffness();
	else
		return Constraint::GetPenaltyStiffness() * Constraint::GetCStepsFact(t);
}

void NodalConstraint::DrawElement() 
{
	Constraint::DrawElement();

	if (GetDrawSizeScalar() == 0) return;

	mbs->SetColor(GetCol());

	Vector3D dir(0.,0.,0.);

	if (!IsLocalNodeConstraint())
	{
		for (int i=1; i<=loccoords.Length(); i++)
		{	
			Vector3D offset(0.,0.,(i-1)*GetDrawSizeScalar());
			Vector3D p = GetNode(1).GetPosD();
			dir.Set(0.,0.,0);
			if (loccoords(i) <=3)
			{
				dir(loccoords(i)) = 1.;
				mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
			}
			else
			{
				dir(3) = 1.;
				mbs->DrawCone(p + (loccoords(i)-2)*(-GetDrawSizeScalar()*0.75)*dir,p+(loccoords(i)-3)*(-GetDrawSizeScalar()*0.75)*dir,(GetDrawSizeScalar()*0.6),6,1);
			}
		}
	}
	else
	{
		// DR 2011-05-11: draw cones correctly if more than 1 direction is constrained
		for (int j=1; j<=loccoords.Length(); j++)
		{	
			int i= ((j+1)%2)+1;					// i = 1 or 2
			double flag = 1;
			if (i==2) flag = -1;
			Vector3D p = GetBody3D(i).GetNodePosD(NodeNum(i));

			dir.Set(0.,0.,0);
			if (loccoords(j) <=3)
			{
				dir(loccoords(j)) = flag;
				mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
			}
		}
		// DR 2011-05-11: old code:
		//for (int i=1; i<=nodes.Length(); i++)
		//{	
		//	double flag = 1;
		//	if (i == 2) flag = -1;

		//	Vector3D p = GetBody3D(i).GetNodePosD(NodeNum(i));

		//	dir.Set(0.,0.,0);
		//	if (loccoords(i) <=3)
		//	{
		//		dir(loccoords(i)) = flag;
		//		mbs->DrawCone(p + (-GetDrawSizeScalar()*0.75)*dir,p,(GetDrawSizeScalar()*0.6),6,1);
		//	}
		//}
	}

};

void NodalConstraint::LinkToElements()
{
	if (!IsLocalNodeConstraint())
	{
		TArray<int> storeltg(IS());
		for (int i=1; i <= IS(); i++)
		{
			storeltg.Add(LTG(i));
		}
		LTGreset();

		// Node dofs
		const Node& node = GetNode(1);
		//Position:
		for (int i=1; i <= loccoords.Length(); i++)
		{
			AddLTG(node.Get(loccoords(i)));
		}
		for (int i=1; i <= loccoords.Length(); i++)
		{
			AddLTG(node.Get(loccoords(i)+node.SOS()));
		}

		// remaining implicit size dofs
		for (int i=1; i <= IS(); i++)
		{
			AddLTG(storeltg(i));
		}
	}
	else
	{
		//store lagrange multiplicator
		TArray<int> storeltg(IS());
		for (int i=1; i <= IS(); i++)
		{
			storeltg.Add(LTG(i));
		}
		LTGreset();

		// add all SOS dofs from the elements
		//Position(first SOS) 
		for (int k=1; k <= NE(); k++)
		{
			for (int i=1; i <= GetBody3D(k).SOS(); i++)
			{
				AddLTG(GetBody3D(k).LTG(i));
			}
		}
		//and Velocity (second SOS):
		for (int k=1; k <= NE(); k++)
		{
			for (int i=1; i <= GetBody3D(k).SOS(); i++)
			{
				AddLTG(GetElem(k).LTG(i+GetBody3D(k).SOS()));
			}
		}
		// implicit dofs corresponding to constraint Lagrange parameters
		for (int i=1; i <= IS(); i++)
		{
			AddLTG(storeltg(i));
		}
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: BasePointJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**begin(ued)**
// name: BasePointJoint (Constraint)
// short description: 
// available formulations: Lagrange multiplier, Penalty Formulation
// type of constraint equation: nonlinear constraint equations
// index formulation available: 
// ground constraint: yes
// element-element constraint: yes
// development status: in progress
// long description:	should replace spherical joint and nodal constraint in future
//								
// class variables:
//  
//**end(ued)**
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//$ DR 2012-07: CheckConsistency added
int BasePointJoint::CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	int rv = Constraint::CheckConsistency(errorstr);
	if (rv) return rv;

	// mixing local and global definition of positions
	int check_localGlobal = 0;
	if(elements(1)) {check_localGlobal = 1;} // local formulation
	for(int i=2; i<=NKinPairs(); i++)
	{
		if( (elements(i) && check_localGlobal) || (!elements(i) && !check_localGlobal) ) {}
		else
		{
			errorstr = mystr("ERROR: BasePointJoint: you are mixing local and global formulation for defining the positions!\n");
			rv = 1;
		}
	}

	// global node to ground for Lagrange multiplier
	if(!(UsePenaltyFormulation()+elements(1)+elements(2)))	// this is just 1 if it is a ground joint with global node and Lagrange multiplier
	{
			errorstr = mystr("ERROR: SphericalJoint: global node to ground not implemented yet for Lagrange multiplier!\n");
			rv = 1;
	}

	// using elements that do not exist
	if(elements(1)>mbs->NE() || elements(2)>mbs->NE())
	{
		UO(UO_LVL_err) << "An element defined in the joint with element number " << GetOwnNum() << " does not exist!\n";
		rv = 1;
	}

	TKinematicsAccessFunctions kaf_elem;
	for(int i=1; i<=NKinPairs(); i++)
	{
		if(elements(i))
		{
			Body3D& b = GetBody3D(i);
			kaf_elem = b.GetKinematicsAccessFunctions();
			if(!IsSuitableElement(kaf_elem))
			{
				UO(UO_LVL_err) << "Element " << elements(i) << ", which is used in the joint " << GetOwnNum() << ", is not suitable for this constraint!\n";
				rv = 1;
			}
		}
	}

	return rv;
}

// for constraints only, these are the necessary (!) access functions!	//$ DR 2013-02-13
void BasePointJoint::GetNecessaryKinematicAccessFunctions(TArray<int> &KinAccFunc, int numberOfKinematicPair)
{
	KinAccFunc.SetLen(0);
	int kaf = (int)(TKinematicsAccessFunctions(TKAF_none));
	// ---------- there are different modi available:
	// penalty + global node
	// to do!

	// penalty + elemNr + loc node
	// lagrange + elemNr + loc node
	kaf = (int)(TKinematicsAccessFunctions(TKAF_node_position+TKAF_node_velocity+TKAF_D_node_pos_D_q));	
	KinAccFunc.Add(kaf);
	
	// penalty + elemNr + locCoord
	// lagrange + elemNr + locCoord
	kaf = (int)(TKinematicsAccessFunctions(TKAF_position+TKAF_velocity+TKAF_D_pos_D_q));	
	KinAccFunc.Add(kaf);
}

//	TKinematicsAccessFunctions kaf = TKinematicsAccessFunctions(TKAF_none);
//	
//	if(elements(numberOfKinematicPair))	// otherwise global nodes
//		kaf = TKinematicsAccessFunctions(kaf+TKAF_position+TKAF_velocity);	
//	if((numberOfKinematicPair == 1) && UseLocalCoordinateSystem()) 
//		kaf = TKinematicsAccessFunctions(kaf+TKAF_rotation_matrix);
//	if(elements(numberOfKinematicPair) && UsePenaltyFormulation())
//		kaf = TKinematicsAccessFunctions(kaf+TKAF_D_pos_D_q);
//	
//	return kaf;
//}

void BasePointJoint::ElementDefaultConstructorInitialization()
{
	IVector dir(3);
	dir.Set3(1,1,1);

	SetConstrainedDirections(dir);
	SetPenaltyFormulation(0);
	SetUseLocalCoordinateSystem(0);
	loccoords.Set2(0,0);
	nodes.Set2(0,0);
	velocity_constraint = 0;
	SetDampingCoeff(0.0);
	elementname = GetElementSpec();
	auto_comp_ground = 0;
	stiffness_in_joint_local_frame = 0;
	draw_local_frame_size = 0;	// do not draw local frame
	//elements.Set2(0,0);	// this would mean elem1 = node1 = 0 --> Pos1 would be ground
	elements.Set2(1,0);
	elementname = GetElementSpec();
}

void BasePointJoint::CopyFrom(const Element& e)
{
	Constraint::CopyFrom(e);
	const BasePointJoint& ce = (const BasePointJoint&)e;

	nodes = ce.nodes;
	loccoords = ce.loccoords;
	velocity_constraint = ce.velocity_constraint;
	spring_stiffness3 = ce.spring_stiffness3;
	dir = ce.dir;
	dpdq = ce.dpdq;
	damping_coeff = ce.damping_coeff;
	JA0i = ce.JA0i;
	auto_comp_ground = ce.auto_comp_ground;
	stiffness_in_joint_local_frame = ce.stiffness_in_joint_local_frame;
	draw_local_frame_size = ce.draw_local_frame_size;

	displ.SetLen(ce.displ.Length());
	displ.SetAll(NULL);
	for(int i=1; i<=ce.displ.Length(); i++)
	{
		MathFunction* mathfun= new MathFunction;
		mathfun=(ce.displ(i));
		displ(i)= mathfun;
	}
}

Vector3D BasePointJoint::GetDrawPosition(int i) 
{
	Vector3D pos;
	Body3D body(GetMBS());


	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			pos = (GetMBS()->GetNode(nodes(i))).GetPosD();
		}
		else
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)	// CMS-Element
			{
				pos = GetBody3D(i).GetNodePosD(nodes(i));
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					pos = GetBody3D(i).GetNodePosD(nodes(i));
				}
				else
				{
					pos = GetBody3D(i).GetPosD(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: GetDrawPosition: no element number defined for body 1";
		return Vector3D(0.0);
	}
	else if (i==2)								// pos2 is ground	
	{
		//return loccoords(2);				// if SetPos2ToGlobalCoord is used, then loccords(2) contains the global(!) coords of ground
		
		Vector3D pos = loccoords(i);					// if SetPos2ToGlobalCoord is used, then loccords(2) contains the global(!) coords of ground
		if(displ.Length() != 0)
		{
			double t_glob = Constraint::GetGlobalTime();
			for(int i=1; i<=displ.Length(); i++)
			{
				if(displ(i) != NULL)
				{
					pos(i) = pos(i)+ (displ(i)->Evaluate(t_glob));
				}
			}
		}

		return pos;				

	}

	return pos;
}



Vector3D BasePointJoint::GetPosition(int i) const
{
	Vector3D pos;
	//Body3D body(GetMBS());

	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			pos = (GetMBS()->GetNode(nodes(i))).GetPos();
		}
		else											
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS) // CMS
			{
				pos = GetBody3D(i).GetNodePos(nodes(i));
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					pos = GetBody3D(i).GetNodePos(nodes(i));
				}
				else
				{
					pos = GetBody3D(i).GetPos(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: no element number defined for body 1";
		return Vector3D(0.0);
	}

	else if (i==2)								// pos2 is ground	
	{
		pos = loccoords(i);					// if SetPos2ToGlobalCoord is used, then loccords(2) contains the global(!) coords of ground
		if(displ.Length() != 0)
		{
			double t_glob = Constraint::GetGlobalTime();
			for(int i=1; i<=displ.Length(); i++)
			{
				if(displ(i) != NULL)
				{
					pos(i) = pos(i)+ (displ(i)->Evaluate(t_glob));
				}
			}
		}

		return pos;				

	}

	return pos;
}


Vector3D BasePointJoint::GetVelocity(int i) const
{
	Vector3D vel;
	//Body3D body(GetMBS());
	if (i <= NKinPairs())					// not ground
	{
		if(elements(i)==0)				// global node number
		{
			vel = (GetMBS()->GetNode(nodes(i))).GetVel();
		}
		else
		{
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)	//CMS
			{
				vel = GetBody3D(i).GetNodeVel(nodes(i));
			}
			else
			{
				if (nodes(i)!=0)						// position defined by local node number
				{
					//vel = GetNode(i).GetVel(); //$ DR 2011-10-31: bugfix
					vel = GetBody3D(i).GetNodeVel(nodes(i));
				}
				else
				{
					vel = GetBody3D(i).GetVel(loccoords(i));
				}
			}
		}
	}
	else if (i==1)								// pos1 is ground	--> not possible
	{
    GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: no element number defined for body 1";
		return Vector3D(0.0);
	}

	else if (i==2)								// pos2 is ground	
	{
		return Vector3D(0.0);
	}

	return vel;
}

void BasePointJoint::SetPosToLocalNode(int i, int elem_nr, int loc_node_nr)
{
	// if i > elements.length resize is done automatically by TArray
	nodes(i) = loc_node_nr;
	elements(i) = elem_nr;
	//if(NKinPairs() < i-1)
	//{ 
	//	GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: you have to specify pos1 before pos2";
	//	return;
	//}
	//else if(NKinPairs() == (i-1))
	//{
	//	elements.Add(elem_nr);
	//}
	//else
	//{
	//	elements.Erase(i);
	//	elements.Insert(i,elem_nr);
	//}
}

void BasePointJoint::SetPosToGlobalNode(int i, int glob_node_nr)
{
	// for global nodes the element number '0' is inserted
	// if i > elements.length resize is done automatically by TArray
	nodes(i) = glob_node_nr;
	elements(i) = 0;

	//if(NKinPairs() < i-1)
	//{ 
	//	GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: you have to specify pos1 before pos2";
	//	return;
	//}
	//else if(NKinPairs() == (i-1))
	//{
	//	elements.Add(0);
	//}
	//else
	//{
	//	elements.Erase(i);
	//	elements.Insert(i,0);
	//}
}

void BasePointJoint::SetPosToLocalCoord(int i, int elem_nr, Vector3D loccoordi)
{
	// if i > elements.length resize is done automatically by TArray
	loccoords(i) = loccoordi;
	elements(i) = elem_nr;
	nodes(i) = 0;
	
	//if(NKinPairs() < i-1)
	//{ 
	//	GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: you have to specify pos1 before pos2";
	//	return;
	//}
	//else if(NKinPairs() == (i-1))
	//{
	//	elements.Add(elem_nr);
	//}
	//else
	//{
	//	elements.Erase(i);
	//	elements.Insert(i,elem_nr);
	//}
}

void BasePointJoint::SetDisplacement(TArray<MathFunction*>& disp)
{
	displ.SetLen(disp.Length());
	for(int i=1; i<=displ.Length(); i++)
	{
		if (disp(i) == NULL)
		{
			displ(i) = NULL;
		}
		else
		{
			MathFunction* mathfun= new MathFunction;
			mathfun=(disp(i));
			displ(i)= mathfun;
		}
	}
}

void BasePointJoint::SetPos2ToGlobalCoord(Vector3D ground)
{
	loccoords(2) = ground;
	elements(2) = 0;
	nodes(2) = 0;
}

// is called if ground = Vector3D(0.)
// the ground will be computed during Initialization
void BasePointJoint::AutoComputeGround()
{
	// SetPos2ToGlobalCoord(GetPos1());		// --> not possible here, because mbs is not yet assembled
	SetPos2ToGlobalCoord(Vector3D(0.));
	auto_comp_ground = 1;
}

void BasePointJoint::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);		//$ DR 2013-05-07 bugfix changed "Element::" to "Constraint::"

	//$JG2012-01
	//simple decision if point is nodal / position based (default)
	//simple decision if point is local (e.g.nodal/pos) (default) or global (nodal/pos)
	//always element numbers in EDC
	//always node numbers in EDC
	//use scalar spring stiffness / damping OR
	//  3+3 Mathfunctions stiffness / damping: Kxx, Kyy, Kzz, Kyz, Kxz, Kxy; Dxx, ...
	//reusable functions for other constraints!
}

int BasePointJoint::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	return Constraint::SetElementData(edc); //$ DR 2013-05-07 bugfix changed "Element::" to "Constraint::"
}

int BasePointJoint::GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables)
{
	// call base class routines
	Constraint::GetAvailableSpecialValues(available_variables);

	// Automatic entries for this class 
	BasePointJoint::GetAvailableSpecialValuesAuto(available_variables);

	// Manual entries for this class
	available_variables.Add(ReadWriteElementDataVariableType("Internal.actor_force",3,0,0.,mystr("force applied to the kinematic pairs due to the spring. range: 1-3 corresponds to force in x-y-z direction"))) ;
	return 0;
}

int BasePointJoint::ReadSingleElementData(ReadWriteElementDataVariableType& RWdata) 		
{
	// call base class routine 	
	int rv = Constraint::ReadSingleElementData(RWdata);
	if (rv == 1) return 1;

	// manual things to read  
	if(RWdata.variable_name.CStrCompare("Internal.actor_force") /*RWdata.variable_name == mystr("Internal.actor_force")*/ )
	{
		if (RWdata.comp1 > 0 && RWdata.comp1 <= 3) //range check
		{
			RWdata.value = GetActorForce(mbs->GetTime(),RWdata.comp1); 
			return 1; 
		}
		else return -2; 
	}

	return ReadSingleElementDataAuto(RWdata);
}

void BasePointJoint::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG
	//PG test: mbs->UO() << "dim of f = " << f.Length() << "\n";
	Vector3D f3;

	//Matrix3D k_mat;
	if(!UseLocalCoordinateSystem() && (dir(1)==0 || dir(2)==0 || dir(3)==0))
	{
		GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: You can only use global coordinate system for three constrained directions in lagrange formulation\n";	
		return;
	}
	else if (!UseLocalCoordinateSystem()) // all directions constrained
	{
		if (MaxIndex()==3)
		{
//<<<<<<< KinematicPairs.cpp
//			if (/*(mbs->DoStaticComputation())||*/(MaxIndex ()==3))
//=======
			if (!IsVelocityConstraint())
//>>>>>>> 1.30
			{
				f3 = GetPos1()-GetPos2();
			}
			else
			{ //velocity constraints:
				f3 = GetVel1()-GetVel2();
				if(displ.Length()) 
				{
					GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: you are using velocity constraints AND displacement";
				}
			}
		}
		else	// MaxIndex < 3
		{ //velocity constraints:
			f3 = GetVel1()-GetVel2();
		}
	}
	else // local coordinate system is used
	{
		if (MaxIndex()==3)
		{
			if (!IsVelocityConstraint())
			{
				f3 = GetRotMati().GetTp()*(GetPos1()-GetPos2());
			}
			else
			{ //velocity constraints:
				f3 = GetRotMatiP().GetTp()*(GetPos1()-GetPos2())+GetRotMati().GetTp()*(GetVel1()-GetVel2());
				if(displ.Length()) 
				{
					GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: you are using velocity constraints AND displacement";
				}
			}
		}
		else	// MaxIndex < 3
		{ //velocity constraints:
			f3 = GetRotMatiP().GetTp()*(GetPos1()-GetPos2())+GetRotMati().GetTp()*(GetVel1()-GetVel2());
		}
	}

	int next=0;
	for(int i=1; i<=3; i++) 
	{
		if(dir(i))
		{
			f(++next)=f3(i);
		}
	}


//	if(UsePenaltyFormulation()) return;  // penalty method --> no Lagrange multiplier --> no EvalG
//	//PG test: mbs->UO() << "dim of f = " << f.Length() << "\n";
//	Vector3D f3;
//
//			Matrix3D k_mat;
//	if(UseLocalCoordinateSystem()||stiffness_in_joint_local_frame)
//	{
//		GetMBS()->UO().InstantMessageText("ERROR: BasePointJoint: EvalG not implemented yet for LocalCoordinateSystem or joint_local_frame");	
//		////$!DR 2012-05-11 test[
//		//Matrix3D rot; // reference rotation
//
//		//Vector3D u = GetPos1()-GetPos2();
//		////global_uo << "u: " << u <<"\n";
//		//k_mat = Matrix3D(dir.X(),dir.Y(),dir.Z());
//
//		//// 0 != A*K*A'*u
//		//rot = GetRotMati();
//		////global_uo << "rot: " << rot <<"\n";
//		//rot.Transpose();
//		//k_mat = k_mat*rot;
//
//		//rot = GetRotMati();
//		//k_mat = rot*k_mat;
//
//		////global_uo << "k_mat: " << k_mat <<"\n";
//
//		//f3 = k_mat*u;	
//		////$!DR 2012-05-11 test]
//	}
//	else
//	{
//		if (MaxIndex()==3)
//		{
////<<<<<<< KinematicPairs.cpp
////			if (/*(mbs->DoStaticComputation())||*/(MaxIndex ()==3))
////=======
//			if (!IsVelocityConstraint())
////>>>>>>> 1.30
//			{
//				f3 = GetPos1()-GetPos2();
//			}
//			else
//			{ //velocity constraints:
//				f3 = GetVel1()-GetVel2();
//				if(displ.Length()) 
//				{
//					GetMBS()->UO() << "ERROR: BasePointJoint: you are using velocity constraints AND displacement";
//				}
//			}
//		}
//		else	// MaxIndex < 3
//		{ //velocity constraints:
//			f3 = GetVel1()-GetVel2();
//		}
//	}
//
//	//int next=0;
//	//for(int i=1; i<=3; i++) 
//	//{
//	//	if(dir(i)==1)
//	//	{
//	//		f(++next)=f3(i);
//	//	}
//	//}
//
//	Vector3D use_this_line = dir;		// if there is no rotation
//
//	///// for rotation still missing
//	//if(UseLocalCoordinateSystem()||stiffness_in_joint_local_frame)
//	//{
//	//	//eig(k_mat)
//	//}
//
//	int next=0;
//	for(int i=1; i<=3; i++) 
//	{
//		if(use_this_line(i))
//		{
//			f(++next)=f3(i);
//		}
//	}

};


void BasePointJoint::EvalF2(Vector& f, double t)
{
	TMStartTimer(29);
	//f = [f1 f2], where f1 is the residual vector of constraint element1 and f2 of constraint element 2
	//Matrix dpdq;
	double sign = 1.;
	int offset = 0;
	int NotAllDirectionsConstrained = 0;

	//if (spring_stiffness3.X()==0 || spring_stiffness3.Y()==0 || spring_stiffness3.Z()==0)  // $ MSax 2013-02-01: added
	//{
	//	NotAllDirectionsConstrained = 1;
	//}

	if(!UsePenaltyFormulation()) return;	// no penalty method --> Lagrange multiplier --> no EvalF2
	
	Vector3D force = ComputeForce(t);
	//Vector3D moment;
	//if (NotAllDirectionsConstrained)  // $ MSax 2013-02-01: added
	//{
	//	moment = ComputeMoment(force);
	//}
	//UO(UO_LVL_dbg1) << "force: " << force << "\n";

	//if(IsLocalNodeConstraint())
	//{
	for (int i=1; i <= NKinPairs(); i++)
	{
		if (i==2) sign = -1.;
		if(elements(i)!=0)
		{
			//======== get Matrix dpdq 
			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)	// CMS
			{
				GetBody3D(i).GetNodedPosdqT(nodes(i), dpdq);
			}
			else
			{
				if (nodes(i)!=0)						// position defined by node number
				{
					GetBody3D(i).GetNodedPosdqT(nodes(i), dpdq);
				}
				else
				{
					GetBody3D(i).GetdPosdqT(loccoords(i),dpdq);
				}
			}

			//======== get Matrix drotdq , $ MSax 2013-02-01: added 
			//if (i==1 && NotAllDirectionsConstrained)
			//{
			//	if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)	// CMS
			//	{
			//		GetMBS()->UO() << "ERROR: BasePointJoint: not implemented yet for unconstrained directions";
			//	}
			//	else
			//	{
			//		if (nodes(i)!=0)						// position defined by node number
			//		{
			//			GetMBS()->UO() << "ERROR: BasePointJoint: not implemented yet for unconstrained directions";
			//		}
			//		else
			//		{
			//			GetBody3D(i).GetdRotdqT(loccoords(i),drotTvdq); // drotTvdq is here drotdqT
			//		}
			//	}
			//}



			//UO(UO_LVL_dbg1) << "i = " << i << ", elem = " << elements(i) << ", dpdq: " << dpdq << "\n";

			if (i==2) offset = GetElem(1).SOS();			// has to be changed for mixed global/local node constraints
			for (int j=1; j<=GetElem(i).SOS(); j++)
			{
				f(j+offset) -= sign*(dpdq(j,1)*force.X() + dpdq(j,2)*force.Y() + dpdq(j,3)*force.Z());
				//if (i==1 && NotAllDirectionsConstrained)
				//{
				//	f(j+offset) -= sign*(drotTvdq(j,1)*moment.X() + drotTvdq(j,2)*moment.Y() + drotTvdq(j,3)*moment.Z());
				//}
				//UO(UO_LVL_dbg1) << "i = " << i << ", elem = "  << elements(i) << "ltg(" << j+offset << ") = " << ltg(j+offset) << "\n";
			}
		}
		else
		{
			//if (NotAllDirectionsConstrained)
			//{
			//	GetMBS()->UO() << "ERROR: BasePointJoint: not implemented yet for unconstrained directions";
			//}
			if (i==2) offset = 3;
			f(1+offset) -= sign*force.X();
			f(2+offset) -= sign*force.Y();
			f(3+offset) -= sign*force.Z();
		}
	}
	//}
	TMStopTimer(29);
};

Vector3D BasePointJoint::ComputeForce(double t) const
{
	Vector3D u; // displacement
	Vector3D v; // velocity
	Vector3D f; // resulting force
	Vector3D k = this->GetPenaltyStiffness3(t); // spring stiffness
	Matrix3D rot; // reference rotation
	Matrix3D k_mat;

	if (UsePenaltyFormulation())
	{
		u = GetPos1() - GetPos2();
		if(UseDamping()) {	v = GetVel1() - GetVel2();}

		//if(UseLocalCoordinateSystem())
		//{
		//	if((mbs->GetElementPtr(elements(1)))->GetType() >= TCMS)
		//	{
		//		GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: ComputeForce: UseLocalCoordinateSystem not implemented for CMS-Element yet";
		//		return Vector3D(0.0);
		//	}
		//	else
		//	{
		//		rot = GetRotMati();
		//		rot.Transpose();
		//		u = rot * u;
		//		//u = (GetRotMati().Transpose()) * u;
		//	}
		//}

		//f.X() = u.X() * k.X();
		//f.Y() = u.Y() * k.Y();
		//f.Z() = u.Z() * k.Z();
		//if(UseDamping()) { f += GetDampingCoeff() * v; }

		//if(UseLocalCoordinateSystem())
		//{
		//	f = GetRotMati()*f;
		//}

		k_mat = Matrix3D(k.X(),k.Y(),k.Z());

		// F = A*K*A'*u
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
		//$ PG 2013-4-20: case of Lagrangian constraint formalism, i.e., if UsePenaltyFormulation()==false
		Vector3D force(XG(1)*dir(1), XG(2)*dir(2), XG(3)*dir(3));
		return force;

		//GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: ComputeForce just implemented for PenaltyFormulation";
		//return Vector3D(0.0);
	}

};

//Vector3D BasePointJoint::ComputeMoment(Vector3D& force) const
//{
//	Vector3D u; // displacement
//	//Vector3D f;
//
//	if (UsePenaltyFormulation())
//	{
//		//f = ComputeForce(t);
//		u = GetPos2() - GetPos1();
//
//		//mbs->MyDrawLine(Vector3D(1,0,0),Vector3D(1,0,0)+u.Cross(f)*0.01, 1);
//		return u.Cross(force);
//	}
//	else
//	{
//		GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: ComputeForce just implemented for PenaltyFormulation";
//		return Vector3D(0.0);
//	}
//};


	double BasePointJoint::GetActorForce(double computation_time, int dir) const //get actor force, with optional direction (x,y,z, etc.)
	{
		if (dir == 0) return ComputeForce(computation_time).Norm();
		if (dir > 3) assert(0 && "BasePointJoint::GetActorForce");

		return ComputeForce(computation_time)(dir);
	}
	void BasePointJoint::LinkToElementsPenalty()
	{
		if (IS()!=0 && SOS()!=0)
		{
			UO(UO_LVL_err) << "Error: BasePointJoint::LinkToElementsPenalty() is not possible for mixed penalty/Lagrange elements\n";
		}
		LTGreset();

		// add all SOS dofs from the elements
		//Position(first SOS) 
		for (int k=1; k <= NKinPairs(); k++)
		{
			if(elements(k)==0)						// global node number
			{
				const Node& node = GetMBS()->GetNode(nodes(k));
				for (int i=1; i <= node.Dim(); i++)
				{
					AddLTG(node.Get(i));
				}
			}
			else													// local node number or local coordinate
			{
				for (int i=1; i <= GetElem(k).SOS(); i++)
				{
					AddLTG(GetElem(k).LTG(i));
				}
			}
		}
		//and Velocity (second SOS):
		for (int k=1; k <= NKinPairs(); k++)
		{
			if(elements(k)==0)						// global node number
			{
				const Node& node = GetMBS()->GetNode(nodes(k));
				for (int i=1; i <= node.Dim(); i++)
				{
					AddLTG(node.Get(i+node.SOS()));
				}
			}
			else													// local node number or local coordinate
			{
				for (int i=1; i <= GetElem(k).SOS(); i++)
				{
					AddLTG(GetElem(k).LTG(i+GetElem(k).SOS()));
				}
			}
		}
	}

	void BasePointJoint::LinkToElementsLagrange()
	{
		if(IsLocalNodeConstraint()){
			for (int i = 1; i <= NE(); i++)	//$ DR 2012-05-02 NE() instead of elements.Length()
			{
				GetElem(i).AddConstraint(this,i);
			}
		}
		else
		{
			//UO(UO_LVL_err) << "Error: BasePointJoint::LinkToElementsLagrange() is not yet implemented for global nodes\n";

			TArray<int> storeltg(IS());
			for (int i=1; i <= IS(); i++)
			{
				storeltg.Add(LTG(i));
			}
			LTGreset();

			//// Node dofs
			const Node& node = GetMBS()->GetNode(nodes(1));

			//Position:
			for (int i=1; i <= dir.Length(); i++)
			{
				if(dir(i)) 
				{
					AddLTG(node.Get(i));
				}
			}

			//Velocity
			for (int i=1; i <= dir.Length(); i++)
			{
				if(dir(i)) 
				{
					AddLTG(node.Get(i+node.SOS()));
				}
			}

			// remaining implicit size dofs
			for (int i=1; i <= IS(); i++)
			{
				AddLTG(storeltg(i));
			}
		}
	}


void BasePointJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	Matrix hTrans;

	if(IsLocalNodeConstraint())
	{

		if (UsePenaltyFormulation()) return;

		double sign = 1;
		if (locelemind == 2) sign = -1;


		if (!UseLocalCoordinateSystem()) // and all directions are constrained (else EvalG writes an error)
		{
			if (nodes(locelemind)!=0)						// position defined by local node number (also for CMS-Element)
			{
				GetBody3D(locelemind).GetNodedPosdqT(nodes(locelemind), dpdq);
			}
			else
			{
				GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
			}
			hTrans = dpdq;
		}
		else //UseLocalCoordinateSystem()
		{
			if (nodes(locelemind)!=0)						// position defined by local node number (also for CMS-Element)
			{
				GetBody3D(locelemind).GetNodedPosdqT(nodes(locelemind), dpdq);
			}
			else
			{
				GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
			}
			// $ MSax 2013-08-08 : [ // sollte auch funktionieren, wenn alle Richtungen gesperrt sind und das Element keine Rot Matrix hat!!!! ==> korrigiert
			if (dir(1) && dir(2) && dir(3)) // all directions constraint
			{
				Matrix Ai_ = GetRotMati();
				hTrans = dpdq*Ai_;
			}
			else // not all directions constraint
			{	
				if (locelemind==1)
				{ // sollte auch funktionieren, wenn alle Richtungen gesperrt sind und das Element keine Rot Matrix hat!!!!
					GetdRotTvdqT(GetPos1()-GetPos2(), loccoords(locelemind), drotTvdq, locelemind);
					Matrix AiB_= GetBody3D(locelemind).GetRotMatrix();
					Matrix JA0i_;

					if (stiffness_in_joint_local_frame)
					{
						JA0i_= JA0i;
					}
					else
					{
						JA0i_ = Matrix3D(1); // not in joint local frame
					}
					hTrans = (drotTvdq+dpdq*AiB_)*JA0i_;
				}
				else
				{
					Matrix3D Ai = GetRotMati();
					Matrix Ai_ = Ai;
					hTrans = dpdq*Ai_;
				}
			}
		}
		// $ MSax 2013-08-08 : ] // sollte auch funktionieren, wenn alle Richtungen gesperrt sind und das Element keine Rot Matrix hat!!!! ==> korrigiert

		//Vector3D pos = GetBody3D(locelemind).GetPos(loccoords(locelemind));
		//Vector3D vel = GetBody3D(locelemind).GetVel(loccoords(locelemind));

		//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
		//C = p_ref+R*p_loc

		for (int i=1; i <= f.Length(); i++)
		{
			//$!DR 2012-03-01: change if just some directions are constrained
			//f(i) -= sign*(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
			int next = 0;
			for(int c=1; c<=3; c++)
			{
				if(dir(c)==1.)
				{
					f(i) -= sign*(hTrans(i,c)*XG(++next));
				}
			}
		}
	}
	else		// global nodes
	{
			// nothing at all?
	}

	//if(IsLocalNodeConstraint())
	//{

	//	if (UsePenaltyFormulation()) return;

	//	double sign = 1;
	//	if (locelemind == 2) sign = -1;


	//	if (nodes(locelemind)!=0)						// position defined by local node number (also for CMS-Element)
	//	{
	//		GetBody3D(locelemind).GetNodedPosdqT(nodes(locelemind), dpdq);
	//	}
	//	else
	//	{
	//		GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	//	}

	//	//Vector3D pos = GetBody3D(locelemind).GetPos(loccoords(locelemind));
	//	//Vector3D vel = GetBody3D(locelemind).GetVel(loccoords(locelemind));

	//	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//	//C = p_ref+R*p_loc

	//	for (int i=1; i <= f.Length(); i++)
	//	{
	//		//$!DR 2012-03-01: change if just some directions are constrained
	//		//f(i) -= sign*(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
	//		int next = 0;
	//		for(int c=1; c<=3; c++)
	//		{
	//			if(dir(c)==1.)
	//			{
	//				f(i) -= sign*(dpdq(i,c)*XG(++next));
	//			}
	//		}
	//	}
	//}
	//else		// global nodes
	//{
	//		// nothing at all?
	//}
};

void BasePointJoint::GetdRotTvdqT(const Vector3D& vloc, const Vector3D& ploc, Matrix& d, int bodyindex)
{
	// MSax: GetdRotTvdqT is derived from GetdRotvdqT, this method is only used in BaseBodyJoint; for more information see scans "GetdRotTvdqT"
	Matrix coeffMat1, coeffMat2, coeffMat3;
	GetBody3D(bodyindex).GetdRotvdqT(Vector3D(1.,0.,0.),ploc,coeffMat1); //coefficient matrix 1
	GetBody3D(bodyindex).GetdRotvdqT(Vector3D(0.,1.,0.),ploc,coeffMat2); //coefficient matrix 2
	GetBody3D(bodyindex).GetdRotvdqT(Vector3D(0.,0.,1.),ploc,coeffMat3); //coefficient matrix 3

	
	d.SetSize(GetBody3D(bodyindex).SOS(),Dim());
	d.SetAll(0);

	for (int i=1; i<=SOS(); i++)
	{
		for (int j=1; j<=3; j++)
		{
			d(i,1) += coeffMat1(i,j)*vloc(j);
		}
	}

	for (int i=1; i<=SOS(); i++)
	{
		for (int j=1; j<=3; j++)
		{
			d(i,2) += coeffMat2(i,j)*vloc(j);
		}
	}

	for (int i=1; i<=SOS(); i++)
	{
		for (int j=1; j<=3; j++)
		{
			d(i,3) += coeffMat3(i,j)*vloc(j);
		}
	}
};

int BasePointJoint::SOS() const 
	{
		if (!UsePenaltyFormulation()) 
		{
			//return 0;
			//return GetNumberOfConstrainedCoords();		// should this be 0?
			if(IsLocalNodeConstraint()) 
			{
				//$!DR 2012-03-01: change if just some directions are constrained
				//return (3 - GetNumberOfConstrainedCoords()) * elements.Length();
				return 0;
			}
			else 
			{
				return GetNumberOfConstrainedCoords();
			}
		}
		else
		{
			int nsos = 0;
			for (int i=1; i <= NKinPairs(); i++)
			{
				if(elements(i)==0) nsos += 3;		// global node number

 				else 
				{
					if(elements(i)<=mbs->NE() && elements(i) != GetOwnNum()) 
					{
						nsos += GetElem(i).SOS(); // $ MSax 2013-08-07 : bugfix --> check if element exist
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


void BasePointJoint::DrawElement() 
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

	for (int i=1; i <= NKinPairs(); i++)
	{
		p = GetDrawPosition(i);
		constr_dir.Flush();

		if ( i== 2) flag = -1;
		if(!mbs->GetSimulationStatus().GetStatusFlag(TSimulationRunning) && t_draw == 0.) k = Vector3D(dir(1),dir(2),dir(3)); //$ AD 2011-09-15: draw all constraints before computation starts (even if they are off at t=0); //$ MaSch 2013-08-19 //$ AH 2013-10-25 
		else if(UsePenaltyFormulation()) k = GetPenaltyStiffness3(t_draw);
		else k = Vector3D(dir(1),dir(2),dir(3));		

		if (UseLocalCoordinateSystem())
		{
			//UO(UO_LVL_err) << "ERROR: BasePointJoint: DrawElement just implemented for global definition of stiffness";
			tmp_dir.Set(0.0,0.0,0.0);

			if((mbs->GetElementPtr(elements(i)))->GetType() >= TCMS)
			{
				//GetMBS()->UO(UO_LVL_err) << "ERROR: BasePointJoint: DrawElement: UseLocalCoordinateSystem not implemented for CMS-Element yet";
				//return;
				mbs->DrawSphere(p, GetDrawSizeScalar());
			}
			else
			{
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

		if(constr_dir.Length()==3)	// all directions are constrained, draw spherical joint //$ DR 2012-12-17 change according to JG
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

Vector3D BasePointJoint::GetPenaltyStiffness3(double t) const 
{
	if(NSteps()==0)
		return spring_stiffness3;
	else
		return spring_stiffness3 * Constraint::GetCStepsFact(t);
}

Matrix3D BasePointJoint::GetRotMati() const
{
	Matrix3D A(1);

	// frame is corotated with body
	//$ DR 2013-05-17: in some joints it is possible to set UseLocalCoordinateSystem to an integer > 1, than this body is used!
	if(UseLocalCoordinateSystem())
	{
		//A = GetBody3D(1).GetRotMatrix(loccoords(1));	//$ DR 2013-05-17 old code
		A = GetBody3D(UseLocalCoordinateSystem()).GetRotMatrix(Vector3D(0));	//$ DR 2013-05-17 old code
	}

	// additional (constant) rotation 
	if(stiffness_in_joint_local_frame)
	{
		A = A*JA0i;
	}
	return A;
}

Matrix3D BasePointJoint::GetRotMatiP() const
{
	Matrix3D A(0);

	if(UseLocalCoordinateSystem())
	{
		A = GetBody3D(1).GetRotMatrixP(loccoords(1));
	}

	if(stiffness_in_joint_local_frame)
	{
		A = A*JA0i;
	}
	return A;
}

Matrix3D BasePointJoint::GetRotMatiD() const
{
	Matrix3D A(1);

	if(UseLocalCoordinateSystem())
	{
		A = GetBody3D(1).GetRotMatrixD(loccoords(1));
	}

	if(stiffness_in_joint_local_frame)
	{
		A = A*JA0i;
	}
	return A;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Spherical: Spherical Joint 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//#define use_sphericaljoint_new	//$!DR 2011-07-06: removed old code
//#ifdef use_sphericaljoint_old		//$!DR 2011-07-06: removed old code

void SphericalJointDEPRECATED::GetElementData(ElementDataContainer& edc) 		//fill in all element data
{
	Constraint::GetElementData(edc);

	ElementData ed;
	ed.SetDouble(GetDrawSizeScalar(), "Draw_size"); ed.SetToolTipText("General drawing size of constraint"); edc.Add(ed);
	ed.SetInt(GetDrawSizeResolution(), "Draw_resolution", 0, 1000); ed.SetToolTipText("Number of quads to approximate round geometries, suggested = 3"); edc.Add(ed);

	if (elements.Length()==1)
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, p_global, "Global_joint_pos"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	} else
	{
		SetElemDataVector3D(edc, loccoords(1), "Local_joint_pos_body1"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
		SetElemDataVector3D(edc, loccoords(2), "Local_joint_pos_body2"); edc.Get(edc.Length()).SetToolTipText("[X, Y, Z]");
	}
	if(UsePenaltyFormulation())
	{
		ed.SetDouble(GetPenaltyStiffness(), "Penalty_stiffness"); ed.SetToolTipText("Penalty stiffness of coordinate constraint"); edc.Add(ed);
	}
}

int SphericalJointDEPRECATED::SetElementData(ElementDataContainer& edc) //set element data according to ElementDataContainer
{
	int rv = Constraint::SetElementData(edc);

	double tmp;
	GetElemDataDouble(mbs, edc, "Draw_size", tmp, 0);
	SetDrawSizeScalar(tmp);
	int dd;
	if (GetElemDataInt(mbs, edc, "Draw_resolution", dd, 0)) SetDrawSizeResolution(dd);

	if (elements.Length()==1)
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Global_joint_pos", p_global);
	} else
	{
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body1", loccoords(1));
		GetElemDataVector3D(GetMBS(), edc, "Local_joint_pos_body2", loccoords(2));
	}

	if(UsePenaltyFormulation())
	{
		GetElemDataDouble(mbs, edc, "Penalty_stiffness", tmp); 		
		SetPenaltyStiffness(tmp);
	}
	return rv;
}


void SphericalJointDEPRECATED::EvalG(Vector& f, double t) 
{
	if(UsePenaltyFormulation()) return;

	if(Constraint::GetCStepsFact(t) == 0.)
	{
		f(1) = XG(1);
		f(2) = XG(2);
		f(3) = XG(3);
		return; //$ AD 2011-09-09: hard turnoff - disable when constraint-factor is 0.
	}

	if (elements.Length()<1 || elements.Length()>2)
	{
		mbs->UO() << "ERROR: SphericalJointDEPRECATED::EvalG, number of elements != 1 or 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		if (elements.Length()==1)
		{
			Vector3D v = GetBody3D(1).GetPos(loccoords(1))-p_global;
			f = v;
		}
		else
			if (elements.Length()==2) 
			{
				Vector3D v = GetBody3D(1).GetPos(loccoords(1))-GetBody3D(2).GetPos(loccoords(2));
				//UO() << "p1=" << GetBody3D(1).GetPos(loccoords(1)) << "\n";
				//UO() << "p2=" << GetBody3D(2).GetPos(loccoords(2)) << "\n";
				f = v;
			}
	}
	else if (MaxIndex()<=2)
	{
		if (elements.Length()==1)
		{
			Vector3D v = GetBody3D(1).GetVel(loccoords(1));
			//UO() << "vel=" << v << "\n";
			f = v;
		}
		else
			if (elements.Length()==2) 
			{
				Vector3D v = GetBody3D(1).GetVel(loccoords(1))-GetBody3D(2).GetVel(loccoords(2));
				f = v;
			}
	}
};

void SphericalJointDEPRECATED::EvalF2(Vector& f, double t)
{
	//f = [f1 f2], where f1 is the residual vector of constraint element1 and f2 of constraint element 2

	if(!UsePenaltyFormulation()) return; 

//	Vector3D force = ComputeForce(t);
	Vector3D force = ComputeForce(t); //$ AD 2011-09-08 conputation steps
	for (int i=1; i <= NE(); i++)
	{
		double sign = 1.;
		int offset = 0;
		if (i==2) 
		{
			sign = -1.;
		  offset = GetElem(1).SOS();
		}
		GetBody3D(i).GetdPosdqT(loccoords(i),dpdq);

		for (int j=1; j<=GetElem(i).SOS(); j++)
		{
			f(j+offset) -= sign*(dpdq(j,1)*force.X() + dpdq(j,2)*force.Y() + dpdq(j,3)*force.Z());
		}
	}
}

//<<<<<<< KinematicPairs.cpp
//Vector3D SphericalJointDEPRECATED::ComputeForce(t) const
//=======
//Vector3D SphericalJoint::ComputeForce(double t) const
//>>>>>>> 1.8
Vector3D SphericalJointDEPRECATED::ComputeForce(double t) const		//$ DR 2011-09-12: this is the solution of the conflict
{
	Vector3D u; // displacement
	Vector3D f; // resulting force
	Vector3D k = this->GetPenaltyStiffness3(t); // spring stiffness
	Matrix3D rot; // reference rotation

	if (UsePenaltyFormulation())
	{
		if (IsGroundJoint())
		{
			Vector3D p1 = GetBody3D(1).GetPos(loccoords(1));
			u =  p1 - p_global;
		}
		else
		{
			Vector3D p1 = GetBody3D(1).GetPos(loccoords(1));
			Vector3D p2 = GetBody3D(2).GetPos(loccoords(2));
			u = p1 - p2;
		}
		if(this->UseLocalCoordinateSystem())
		{
			rot = GetBody3D(1).GetRotMatrix();
			rot.Transpose();
			u = rot * u;
		}
		f.X() = u.X() * k.X();
		f.Y() = u.Y() * k.Y();
		f.Z() = u.Z() * k.Z();

		if(this->UseLocalCoordinateSystem())
		{
			rot = GetBody3D(1).GetRotMatrix();
			f = rot * f;
		}

		return f;
	}
	else
	{
		return this->GetElem(1).XG(1);
	}
	//if (UsePenaltyFormulation())
	//{
	//	Vector3D u1(0.,0.,0.);
	//	for(int i=1; i<=Dim(); i++)
	//		u1(i) = GetElem(1).XG(loccoords(i)) - GetElem(1).GetXInit()(loccoords(i)));
	//	
	//	Vector3D u2(0.,0.,0.);
	//	if (elements.Length()==2) 
	//	{
	//		for(int i=1; i<=Dim(); i++)
	//			u2(i) = GetElem(2).XG(loccoords(i)) - GetElem(2).GetXInit()(loccoords(i)));
	//	}
	//	return GetPenaltyStiffness() * (u1-u2);
	//}
	//else
	//{
	//	Vector3D rv;
	//	for(int i=1; i<=Dim(); i++)
	//			rv(i) = GetElem(1).XG(i);
	//	return rv; //Lagrange parameter
	//}
}


void SphericalJointDEPRECATED::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	if (UsePenaltyFormulation()) return;

	double sign = 1;
	if (locelemind == 2) sign = -1;

	GetBody3D(locelemind).GetdPosdqT(loccoords(locelemind),dpdq);
	Vector3D pos = GetBody3D(locelemind).GetPos(loccoords(locelemind));
	Vector3D vel = GetBody3D(locelemind).GetVel(loccoords(locelemind));

	//-C_q^T \lambda = [dC/dx; dC/dy; dC/dphi]^T [\lambda1;\lambda2] 
	//C = p_ref+R*p_loc

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= sign*(dpdq(i,1)*XG(1)+dpdq(i,2)*XG(2)+dpdq(i,3)*XG(3));
	}

};
Vector3D SphericalJointDEPRECATED::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void SphericalJointDEPRECATED::DrawElement() 
{
	Constraint::DrawElement();

	int res = GetDrawSizeResolution();
	if (res < 2) res = 3;

	if (GetDrawSizeScalar() != 0)
	{
		mbs->SetColor(GetCol());
		mbs->DrawSphere(GetBody3D(1).GetPosD(loccoords(1)),0.5*GetDrawSizeScalar(),res);
		if (elements.Length()==1)
		{
			mbs->SetColor(Vector3D(0.8,0.1,0.1));
			mbs->DrawSphere(p_global,0.53*GetDrawSizeScalar(),5);
		}
		else
		{
			mbs->SetColor(Vector3D(0.8,0.1,0.1));
			mbs->DrawSphere(GetBody3D(2).GetPosD(loccoords(2)),0.53*GetDrawSizeScalar(),5);
		}
	}
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical	Cylindrical			
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Cylindrical: CylindricalPointJoint
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CylindricalPointJoint::Initialize() 
{
	//find orthogonal vectors:
	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	Vector3D vi1,vi2,vi3;
	Vector3D gp1 = GetBody3D(1).GetPos(lp1);
	Vector3D gp2 = GetBody3D(2).GetPos(lp2);
	vi3 = gp1 - gp2;

	vi3.SetNormalBasis(vi1,vi2);

	Matrix3D RTi=GetBody3D(1).GetRotMatrix(lp1).GetTp();

	vi1 = RTi*vi1;
	vi2 = RTi*vi2;

	loccoords(3) = vi1;
	loccoords(4) = vi2;

};

void CylindricalPointJoint::EvalG(Vector& f, double t) 
{
	if (elements.Length() != 2)
	{
		mbs->UO() << "ERROR: CylindricalPointJoint::EvalG, number of elements != 2\n"; return;
	}

	if (MaxIndex()==3)
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lvi1 = loccoords(3);
		const Vector3D& lvi2 = loccoords(4);

		Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);

		Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A1*lvi1;
		Vector3D gvi2 = A1*lvi2;

		f(1) = gvi1 * rpij;
		f(2) = gvi2 * rpij;
	}
	else if (MaxIndex()<=2)
	{
		const Vector3D& lp1 = loccoords(1);
		const Vector3D& lp2 = loccoords(2);
		const Vector3D& lvi1 = loccoords(3);
		const Vector3D& lvi2 = loccoords(4);

		Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);
		Vector3D rpijp = GetBody3D(1).GetVel(lp1)-GetBody3D(2).GetVel(lp2);

		Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Matrix3D A1p= GetBody3D(1).GetRotMatrixP(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A1*lvi1;
		Vector3D gvi2 = A1*lvi2;
		Vector3D gvi1p = A1p*lvi1;
		Vector3D gvi2p = A1p*lvi2;

		f(1) = gvi1 * rpijp + gvi1p * rpij;
		f(2) = gvi2 * rpijp + gvi2p * rpij;
	}
};

void CylindricalPointJoint::AddElementCqTLambda(double t, int locelemind, Vector& f) 
{
	hmat.SetSize(f.Length(),2);


	const Vector3D& lp1 = loccoords(1);
	const Vector3D& lp2 = loccoords(2);
	const Vector3D& lvi1 = loccoords(3);
	const Vector3D& lvi2 = loccoords(4);

	if (locelemind==1)
	{
		Vector3D rpij = GetBody3D(1).GetPos(lp1)-GetBody3D(2).GetPos(lp2);

		Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A1*lvi1;
		Vector3D gvi2 = A1*lvi2;


		GetBody3D(1).GetdRotvdqT(lvi1,lp1,dpdq); //Hi1
		Mult(dpdq,rpij,hvec);
		hmat.SetColVec(hvec,1);

		GetBody3D(1).GetdRotvdqT(lvi2,lp1,dpdq); //Hi2
		Mult(dpdq,rpij,hvec);
		hmat.SetColVec(hvec,2);

		GetBody3D(1).GetdPosdqT(lp1,dpdq); //Hip
		Mult(dpdq,gvi1,hvec);
		hmat.AddColVec(1,hvec);
		Mult(dpdq,gvi2,hvec);
		hmat.AddColVec(2,hvec);
		//UO() << "hmat1=" << hmat << "\n";
	}
	else
	{
		Matrix3D A1 = GetBody3D(1).GetRotMatrix(lp1); //for deformable bodies: needs to be changed to Body3D().GetRotVec(lpos,locvec);
		Vector3D gvi1 = A1*lvi1;
		Vector3D gvi2 = A1*lvi2;

		GetBody3D(2).GetdPosdqT(lp2,dpdq); //Hip
		Mult(dpdq,gvi1,hvec); hvec *= -1;
		hmat.SetColVec(hvec,1);
		Mult(dpdq,gvi2,hvec); hvec *= -1;
		hmat.SetColVec(hvec,2);
		//UO() << "hmat2=" << hmat << "\n";
	}
	//UO() << "hmat=" << hmat;
	hvec.SetLen(2);
	hvec(1) = XG(1);
	hvec(2) = XG(2);

	for (int i=1; i <= f.Length(); i++)
	{
		f(i) -= hmat(i,1)*XG(1)+hmat(i,2)*XG(2);
	}
};

Vector3D CylindricalPointJoint::GetRefPosD()	const 
{
	return GetBody3D(1).GetPosD(loccoords(1));
}

void CylindricalPointJoint::DrawElement() 
{
	Constraint::DrawElement();

	mbs->SetColor(GetCol());
	//...
	if (GetDrawSizeScalar() == 0) return;

	Vector3D rot1, rot2;
	rot1 = GetBody3D(1).GetRotMatrixD(loccoords(1))*(loccoords(3).Cross(loccoords(4)));
	rot1.Normalize();
	rot1 *= 0.5*GetDrawSizeAxisLength();
	rot2 = loccoords(6);
	rot2.Normalize();
	rot2 *= 0.5*GetDrawSizeAxisLength();
	Vector3D p = GetBody3D(1).GetPosD(loccoords(1));
	mbs->DrawZyl(p+rot1,p-rot1,GetDrawSizeScalar(),20);
	mbs->SetColor(colgrey2);
	mbs->DrawZyl(p+1.5*rot2,p-1.5*rot2,0.14*GetDrawSizeScalar(),8);
};