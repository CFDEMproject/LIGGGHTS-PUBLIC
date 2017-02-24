//#**************************************************************
//#
//# filename:             AverageConstraint.h
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
 

#ifndef __AVERAGE_CONSTRAINT__
#define __AVERAGE_CONSTRAINT__

#include "KinematicPairs.h"

typedef enum {ACT_None = 0, ACT_Loccoords=1, ACT_Nodes=2} AverageConstraintType;

class AverageConstraint: public BasePointJoint//$EDC$[beginclass,classname=AverageConstraint,parentclassname=BasePointJoint,addelementtype=TAEconstraint+TAENotInRelease,addelementtypename=AverageConstraint,
//texdescription="
//Bla The AverageConstraint is a Connector which acts on two sets A and B of weighted body points. 
//Each set corresponds two a seperate body of the MBS, set B may also correspond to the global coordinate system (ground). 
//Body points are wheter specified as pairs of element number and local position, or as element number and local node number. 
//If an element number is zero, then a node number is interpreted as global node number and a local position is interpreted as global (ground) position.
//",
//figure="AverageConstraint_1,Geometry center $\mathbf r^A =\sum_{n=1}^{N^A} w^A_n \mathbf r^A_n$ and local frame $\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right)$ of the Average-Constraint.",
//texdescriptionDOF="
//If the constraint follows the Lagrange formalism, then one scalar (Lagrange) multiplier for each \emph{moment setting} is used as DOF. Otherwise (penalty formalism), there are no DOFs.
//",
//texdescriptionGeometry="
//The geometry of the constraint is defined by the center (weighted average of global positions) $\mathbf r^A =\sum_i w^A_i \mathbf r^A_i$ of the weighted point-set $({\mathbf r}_i^A, w_i^A)_{i=1}^{N^A}$ on body $A$, and a local frame $\left(\mathbf e_1^A, \mathbf e_2^A, \mathbf e_3^A\right)$.
//This local frame is defined by the user by setting a frame $\left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)$ and the integer \code{use\_local\_coordinate\_system} \developer{(see the method \code{SetUseLocalCoordinateSystem(int i)})}. If the user sets \code{use\_local\_coordinate\_system}$=0$, then the local frame of the constraint equals the user defined frame in the global coordinate system, i.e., $$\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right) = \left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)\,$$ but if this integer is positive, say \code{use\_local\_coordinate\_system}$=i>0$, then the local frame of the constraint equals the user defined frame corotated by the rotation matrix of element $i$ (in the global element index list of the multibody system):
//$$\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right) = \mathbf{A}_i \left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)\,.$$
//If the element number $i$ is a flexible body, then the rotation matrix is taken from its local position $(0,0,0)$. Note, that the rotation matrix $\mathbf{A}_i$ depends on the generalized coordinates $\mathbf q_i$ of element number $i$.\\
//The user further determines one or more moment settings \developer{(see the method \code{SetMoments(\ldots)})}, each of which consists of 
//\begin{itemize} 
//  \item an integer $d \in \{1,2,3\}$ defining the constrained direction $\mathbf e_d$ (of the local frame),\\ 
//  \item an integer $p \in \{0,1,2,\ldots\}$ addressing the polynomial deggree of the moment,\\ 
//  \item an integer $l \in \{1,2,3\}$ specifying the lever direction $\mathbf e_l$ (again of the local frame) for calculation of the integrated moment, and\\ 
//  \item potentially a scalar penalty coefficient $k>0$, if the constraint follows the penalty instead of the Lagrange formalism.
//\end{itemize} 
//Those parameters are collected in a \code{TArray<int3> moment\_settings} of length $M>0$, s.t. for $m \in \{0,\ldots,M\}$ we collect \code{moment\_settings(m)} $= (d_m,p_m,l_m)$, and the \code{Vector penalty\_coefficients} is either of length $0$ (Lagrange) or length $M$ with \code{penalty\_coefficients(m)} $=k_m$. 
//As an example one may define \code{moment\_settings} as $$\begin{pmatrix}\left(1,0,0\right)\\\left(2,0,0\right)\\\left(3,0,0\right)\\\left(3,1,2\right)\end{pmatrix}\,,$$ which results in constraining all components of the center position difference $\mathbf r^{AB}$, plus the first moment for components in $\mathbf e_3$-direction, and a lever in $\mathbf e_2$-direction.
//",
//texdescriptionEquations="
//Let us first assume,  
//\begin{itemize}
//  \item that for a fixed $m$ the degree of moment $p_m$ equals zero, \\
//  \item and that the local coordinate system of the constraint is not co-rotated with the local frame of the weighted body point set A, i.e., $\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right) = \left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)$.
//\end{itemize}
//Then, constraining a certain component of the center position difference $\mathbf r^{AB} = \mathbf r^B - \mathbf r^A$ between the both sets of weighted body points may be acchieved by the following constraint equation:
//\begin{equation*}
//  C_m(\mathbf q_A, \mathbf q_B) = 
//  \sum_{n=1}^{N^A} \mathbf r_n^A(\mathbf q_A)^T \mathbf{\bar e}_{d_m} w_n^A -
//  \sum_{n=1}^{N^B} \mathbf r_n^B(\mathbf q_B)^T \mathbf{\bar e}_{d_m} w_n^B\,,
//\end{equation*} 
//where the two vectors $\mathbf q_A$ and $\mathbf q_B$ denote the DOFs corresponding to body A or B of the multibody system. 
//Apart from the position, also the moment of a certain degree $p_m$ w.r.t. a lever direction $\mathbf e_{l_m}$ may be constrained: 
//Let us assume,  
//\begin{itemize}
//  \item that for a fixed $m$ the degree of moment $p_m$ is a positive integer, \\
//  \item but still the local coordinate system of the constraint is not co-rotated with the local frame of the weighted body point set A, i.e., $\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right) = \left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)$.
//\end{itemize}
//Then the corresponding constraint equations read
//\begin{equation*}
//  C_m(\mathbf q_A, \mathbf q_B) = 
//  \sum_{n=1}^{N^A} \mathbf r_n^A(\mathbf q_A)^T \mathbf{\bar e}_{d_m} \, w_n^A \left(y_{m,n}^A\right)^{p_{m}} -
//  \sum_{n=1}^{N^B} \mathbf r_n^B(\mathbf q_B)^T \mathbf{\bar e}_{d_m} \, w_n^B \left(y_{m,n}^B\right)^{p_{m}}
//\end{equation*}
//with $y_{m,n}^{A/B}$ denoting the distance between the position of a body point and the weighted center position (both in reference configuration) projected into lever direction, e.g., for set A:
//\begin{equation*}
//  y_{m,n}^{A} = \left(\mathbf r_{0,n}^A-\mathbf r_0^A\right)^T \mathbf{\bar e}_{l_m}\,,\quad \text{where}\quad \mathbf r_0^A = \sum_{n=1}^{N^A} w_n^A \mathbf r_{0,n}^A\,.
//\end{equation*}
//Finally, the definition of $C_m$ in the most complicated case
//\begin{itemize}
//  \item for a fixed $m$ the degree of moment $p_m$ is a positive integer, \\
//  \item and the local coordinate system of the constraint is co-rotated with element number $i$, i.e.,
//	$$\left(\mathbf e_1, \mathbf e_2, \mathbf e_3\right) = \mathbf{A}_i \left(\mathbf{\bar e}_1, \mathbf{\bar e}_2, \mathbf{\bar e}_3\right)\,,$$
//\end{itemize}
//includes the rotation matrix $\mathbf{A}_i$ of element $i$ (at its local position $(0,0,0)$, and depending on the generalized coordinates $\mathbf q_i$):
//\begin{equation*}
//  C_m(\mathbf q_A, \mathbf q_B) = 
//  \sum_{n=1}^{N^A} \mathbf r_n^A(\mathbf q_A)^T \mathbf{A}_i(\mathbf q_i)\,\mathbf{\bar e}_{d_m} \, w_n^A \left(y_{m,n}^A\right)^{p_{m}} -
//  \sum_{n=1}^{N^B} \mathbf r_n^B(\mathbf q_B)^T \mathbf{A}_i(\mathbf q_i)\,\mathbf{\bar e}_{d_m} \, w_n^B \left(y_{m,n}^B\right)^{p_{m}}\,.
//\end{equation*}
//with $y_{m,n}^{A/B}$ now defined
//\begin{equation*}
//  y_{m,n}^{A} = \left(\mathbf r_{0,n}^A-\mathbf r_0^A\right)^T \mathbf{A}_{0,i} \mathbf{\bar e}_{l_m}\,.
//\end{equation*}
// and $\mathbf{A}_{0,i}$ denoting the rotation matrix of element number $i$ in reference configuration at its local position $(0,0,0)$.
//",
//modus="{element to ground}{Position2.element\_numbers must only contain zeros.}",
//modus="{element to element}{Position2.element\_numbers must not contain any 0}",
//modus="{Lagrange}{If Physics.use\_penalty\_formulation = 0, then no stiffness (=Physics.Penalty.weights) is used.}",
//modus="{Penalty}{If Physics.use\_penalty\_formulation = 1, then the stiffnesses (=Physics.Penalty.weights) are used in order to calculate the forces, that are applied to the kinematic pairs.}"]

//------------------------------- old docu, not true anymore, but valuable --------------------------
//figure="AverageConstraint_2,Computation of the rotation matrix",
// //\developer{See Fig.~\ref{AverageConstraintfigure2}
//\begin{equation*}
//  \mathbf r_0^A = \sum w_i \mathbf r_{i,0}^A 
//\end{equation*}  //\begin{equation*}
//  \mathbf r^A = \sum w_i \mathbf r_{i}^A 
//\end{equation*}  \begin{equation*}
//  \mathbf r_{i,0}^A = \mathbf r_0^A + \mathbf A_0 \mathbf x_{i,0}
//\end{equation*}  \begin{equation*}
//  \mathbf r_{i}^A = \mathbf r^A + \mathbf A \mathbf x_{i}
//\end{equation*}  
//for small deflections:
//\begin{equation*}
//  \mathbf x_{i,0} = \mathbf x_{i}
//\end{equation*}  \begin{equation*}
//  \mathbf A_0^{-1}(\mathbf r_{i,0}^A - \mathbf r_0^A) = \mathbf A^{-1}(\mathbf r_{i}^A - \mathbf r^A)
//\end{equation*}  \begin{equation*}
//  \mathbf A \mathbf A_0^{-1}(\mathbf r_{i,0}^A - \mathbf r_0^A) = \mathbf r_{i}^A - \mathbf r^A
//\end{equation*}  \begin{equation*}
//  \mathbf A \mathbf b_i = \mathbf c_i 
//\end{equation*} 
//with constant vector $\mathbf b_i$ and easy to compute vector $\mathbf c_i$ and $i=1,2,3,...$
//
//Question 1: explicit - implicit\\
//explicit: per tk in StartTimeStep: A(tk) = f(tk-1) and therefore constant for time step tk\\
//implicit: A(tk= = f(r(tk)) and therefore in ComputeMomentum additional f \\
//
//Question 2: (How) is it possible to determine the rotation Matrix if the points are a) in a plane and b) on a line\\
//
//Question 3: Should question 2 should be solved in AverageConstraint, linalg or ReferenceFrame?\\
//
//concerning Q2:\\
//a) coord-system is well defined with 3 points (r1, r2, r3) in the plane, except of the sign of the vector perpendicular to the plane. Solution: define an artificial Point outside the plane (P  = r2+n)\\
//
//b) transformation of a vector in reference config (the line) to a vector in the actual config (the line now). it should by possible to get the transformation matrix.}
//------------------------------- end of old docu------------------------------------------------
{
public:
	AverageConstraint(MBS* mbsi):BasePointJoint(mbsi)
	{	
		ElementDefaultConstructorInitialization();
	};

	AverageConstraint(const AverageConstraint& ct): BasePointJoint(ct.mbs)
	{
		ElementDefaultConstructorInitialization();
		CopyFrom(ct);
	};

	~AverageConstraint()
	{
	};

	virtual void ElementDefaultConstructorInitialization();
	virtual void Initialize();

	virtual int CheckConsistency(mystr& errorstr); //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	virtual void GetNecessaryKinematicAccessFunctions(TArray<int> &KAF, int numberOfKinematicPair);

	virtual Element* GetCopy()
	{
		Element* ec = new AverageConstraint(*this);
		return ec;
	}

	virtual void CopyFrom(const Element& e);
	virtual const char* GetElementSpec() const {return "AverageConstraint";}

	virtual int IS() const // ok
	{
		if(UsePenaltyFormulation()) return 0;
		else
		{
			return moment_settings.Length();	
		}
	}

	virtual int SOS() const // ok
	{
		if(UsePenaltyFormulation())
		{
			int sos=0;

			for(int i=1; i<=NE_nodouble(); i++)
			{
				sos += GetElem_nodouble(i).SOS();
			}

			return sos;
		}
		else return 0; 
	}

	virtual int ES() const { return 0; }
	virtual int Dim() const 
	{
		//if(mbs->NE())
		//	return GetElement(1,1).Dim();
		//else
			return 3;
	}

	virtual void SetElement_nodouble();
	virtual int NE_nodouble() const;
	
	virtual int NKinPairs() const 
	{
		if(elements2(1)) return 2;
		else return 1;
	}

	// k = kinematic pair, i=internal number
	virtual const Element& GetElement(int k, int i) const 
	{
		if(k==1) return GetElem(i);
		else return mbs->GetElement(elements2(i));
	}
	virtual Element& GetElement(int k, int i)
	{
		if(k==1) return GetElem(i);
		else return mbs->GetElement(elements2(i));
	}
	virtual const Body3D& GetBody3D(int k, int i) const
	{
		return (Body3D&) GetElement(k,i);
	}
	virtual Body3D& GetBody3D(int k, int i)
	{
		return (Body3D&) GetElement(k,i);
	}
	virtual const Element& GetElem_nodouble(int i) const {return mbs->GetElement(elements_nodouble(i));}
	virtual Element& GetElem_nodouble(int i) {return mbs->GetElement(elements_nodouble(i));}
	virtual int NElements(int k) const 
	{
		if (k==1)
		{
			return elements.Length();
		}
		return elements2.Length();
	}

	virtual Vector3D GetRefConfPos(int kinpair, int pos_nr) const;

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data according to ElementDataContainer
  virtual int GetAvailableSpecialValues(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //for special value sensor
	virtual int ReadSingleElementData(ReadWriteElementDataVariableType& RWdata); 		//for special value sensor
	virtual int WriteSingleElementData(const ReadWriteElementDataVariableType& RWdata); // for IOElementDataModifier

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
  virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void DrawElement();
	//virtual void StartTimeStep();  // compute center of material points on each body 1/2 --> center1 / center2 (at start of each time step)

	//virtual Vector3D GetPosition(int i) const;
	//virtual Vector3D GetDrawPosition(int i) const;

	virtual void EvalF2(Vector& f, double t);	
	//virtual Vector3D ComputeForce(double t, int moment_idx) const;
	//virtual double ComputeForceComponent(double t, int kinpair, int pos_nr, int moment_idx) const;
	//virtual Vector3D ComputeForce(double t, int kinpair, int pos_nr) const;
	//virtual Vector3D ComputeForce(double t, int kinpair=1) const;
	virtual Vector3D ComputeForce(double t, int pos_nr, int kinpair=1) const;

	virtual double ComputeMomentum(double t, int moment_index) const;
	virtual void ComputeMomentumDerivative(double t, int locelemind, const Vector& multiplier, Vector& f);
	
	virtual double GetMomentumPolynomial(const int kinpair, const int pos_nr, const int moment_index) const;

	virtual void AddElementCqTLambda(double t, int locelemind, Vector&f);
	virtual void EvalG(Vector &f, double t);

	virtual void LinkToElementsPenalty();
	virtual void LinkToElementsLagrange();

	//$ PG 2013-6-14:[ ElementsShareDOFs(): 
	// for multi-point constraints, which constrain a set of elements
	// to ground or to another set of elements (e.g. AverageConstraint),
	// elements of the same set might share the same global dofs. 
	// if this is the case, this function returns 1, else zero. 
	// this information is used in the assembling of the global system 
	// see also void Element::JacobianG(double t, Matrix& m, IVector& colref)
	// called by void MultiBodySystem::LocalJacobianG(SparseMatrix& m, Vector& x).
	virtual int ElementsShareDOFs() {	return elements_share_dofs; }
	virtual void SetElementsShareDOFS() { elements_share_dofs = 1; }

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const { return GetPosition(1);}

	// do not use these functions!
	void SetPos1ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint: SetPos1ToLocalNode not implemented.";}
	void SetPos2ToLocalNode(int elem_nr, int loc_node_nr)		{GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint: SetPos2ToLocalNode not implemented.";}
	void SetPos1ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint: SetPos1ToLocalCoord not implemented.";}
	void SetPos2ToLocalCoord(int elem_nr, Vector3D loccoord){GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint: SetPos2ToLocalCoord not implemented.";}
	void SetPos2ToGlobalCoord(Vector3D ground)						  {GetMBS()->UO(UO_LVL_err) << "ERROR: AverageConstraint: SetPos2ToGlobalCoord not implemented.";}
	void SetPositions(TArray<int> element_numbers,  TArray<Vector3D> coordinates);
	
	void SetAverageConstraint_LocalNodes_to_LocalNodes(TArray<int> element_numbers1, TArray<int> node_numbers1, Vector weights1, TArray<int> element_numbers2, TArray<int> node_numbers2, Vector weights2);
	void SetAverageConstraint_LocalCoords_to_LocalCoords(TArray<int> element_numbers1, TArray<Vector3D> loccoords1, Vector weights1, TArray<int> element_numbers2, TArray<Vector3D> loccoords2, Vector weights2);
	void SetAverageConstraint_LocalNodes_to_GlobalPos(TArray<int> element_numbers1, TArray<int> node_numbers1, Vector weights1, TArray<Vector3D> ground, Vector weights2);
	void SetAverageConstraint_LocalCoords_to_GlobalPos(TArray<int> element_numbers1, TArray<Vector3D> loccoords1, Vector weights1, TArray<Vector3D> ground, Vector weights2);
	//$ PG 2013-4-24: [ just for testing with existing model
	void SetEdgeConstraint_Displacement(const IVector& element_numbers1, const TArray<Vector3D>& loccoords1, const Vector& weights, const TArray<int3>& moments_i, const Matrix& prescribed_displacements = Matrix(0,0), const Vector& stiffnesses = Vector(0));
	void SetScalingFactor(double scaling_factor_i) {scaling_factor = scaling_factor_i; }
	//$ PG 2013-4-24: ]
	void SetMoments(const TArray<int3>& moments_i, const Vector& stiffnesses_i = Vector(0))
	{
		moment_settings = moments_i;						// moment_settings(i) = [constrained coordinate, moment order, lever direction]
		if (stiffnesses_i.Length() && stiffnesses_i*stiffnesses_i > 1e-14)
		{
			SetPenaltyFormulation(1);
			stiffnesses = stiffnesses_i;		// just for penalty formulation: if moment_settings.Length()>0, each moment can be weighted with a stiffness, e.g. c1*fx + c2*fy
		}
	}	

	int GetConstrainedDirection(int i) const {return moment_settings(i).Get(1);}
	int GetLeverDirection(int i) const {return moment_settings(i).Get(3);}
	int GetMomentOrder(int i) const {return moment_settings(i).Get(2);}
	Vector3D GetFrameVector(int i) const;

protected:
	virtual void SetMomentumPolynomial();  // should only be called once (in Initialize())
	virtual AverageConstraintType GetAverageConstraintType() const;

	// remove entries
	//EDC Vector3D spring_stiffness3;						//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]	
	//EDC int elements(1);											//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position1"]
	//EDC int elements(2);											//$EDC$[varaccess,remove,EDCvarname="element_number",EDCfolder="Position2"]
	//EDC Vector3D loccoords(1);								//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position1"]
	//EDC Vector3D loccoords(2);								//$EDC$[varaccess,remove,EDCvarname="position",EDCfolder="Position2"]
	//EDC int nodes(1);													//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position1"]
	//EDC int nodes(2);													//$EDC$[varaccess,remove,EDCvarname="node_number",EDCfolder="Position2"]
	//EDC	double spring_stiffness;							//$EDC$[varaccess,remove,EDCvarname="spring_stiffness",EDCfolder="Physics.Penalty"]
	//EDC Vector3D spring_stiffness3;						//$EDC$[varaccess,remove,EDCvarname="spring_stiffness_vector",EDCfolder="Physics.Penalty"]	
	//EDC	double damping_coeff;									//$EDC$[varaccess,remove,EDCvarname="damping",EDCfolder="Physics.Penalty"]
	//EDC TArray<int> dir;											//$EDC$[varaccess,remove,EDCvarname="constrained_directions",EDCfolder="Physics.Lagrange"]
	//EDC int use_local_coordinate_system;			//$EDC$[varaccess,remove,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry"]
	//EDC Matrix3D JA0i;												//$EDC$[varaccess,remove,EDCvarname="joint_local_frame",EDCfolder="Geometry"]

	// Geometry ---------------------------:
	//EDC TArray<int> nodes;										//$EDC$[varaccess,EDCvarname="node_numbers",EDCfolder="Position1",variable_length_vector,tooltiptext="(local) node numbers of the kinematic pair 1"]
	TArray<int>	nodes2;													//$EDC$[varaccess,EDCvarname="node_numbers",EDCfolder="Position2",variable_length_vector,tooltiptext="(local) node numbers of the kinematic pair 2"]
	//EDC int use_local_coordinate_system;			//$EDC$[varaccess,EDCvarname="use_local_coordinate_system",EDCfolder="Geometry",tooltiptext="0 ...use global coordinates, i=1,2,..n ...use local coordinate system of body i, evaluated at local position (0,0,0)"]
  //EDC Matrix3D JA0i;												//$EDC$[varaccess,EDCvarname="joint_local_frame",EDCfolder="Geometry",tooltiptext="Prerotate stiffness vector w.r.t. global coordinate system or local coordinate system of body i. Just used if use_joint_local_frame == 1"]	

	//TArray<Vector3D> loccoords;	already exists in BasePointJoint
	TArray<Vector3D> loccoords2;		// no auto access possible, see GetElementData
	//EDC TArray<int> elements;									//$EDC$[varaccess,EDCvarname="element_numbers",EDCfolder="Position1",variable_length_vector,tooltiptext="element numbers of the kinematic pair 1"]
	TArray<int> elements2;											//$EDC$[varaccess,EDCvarname="element_numbers",EDCfolder="Position2",variable_length_vector,tooltiptext="element numbers of the kinematic pair 2"]
	Vector weights1;														//$EDC$[varaccess,EDCvarname="weights",EDCfolder="Position1",variable_length_vector,tooltiptext="weights of the points of kinematic pair 1"]
	Vector weights2;														//$EDC$[varaccess,EDCvarname="weights",EDCfolder="Position2",variable_length_vector,tooltiptext="weights of the points of kinematic pair 2"]
	double sumWeights1;
	double sumWeights2;
	Matrix polynomial_1;
	Matrix polynomial_2;

	int elements_share_dofs;                    //$EDC$[varaccess,EDCvarname="elements_share_dofs",EDCfolder="Geometry",int_bool,tooltiptext="check, if not just one, but a set of elements are constrained, and some of those elements share the same dofs"]

	// Physics ---------------------------:
	TArray<int3> moment_settings;		// moment_settings(i) = [constrained coordinate, moment order, lever direction], e.g. r_z*y^4 = [3,4,2]; r_x*z^2 = [1,2,3]; r_y*s^3 = [2,3,0] (with 's': arc_length); r_y = [2,0,0]; // no auto access possible, see GetElementData
	Vector stiffnesses;			//$EDC$[varaccess,EDCvarname="weights",EDCfolder="Physics.Penalty",variable_length_vector,tooltiptext="each moment can be weighted with a stiffness, e.g. c1*fx + c2*fy"]

	// internal:
	TArray<int> elements_nodouble;
	
	double scaling_factor;

	int	is_nodal_constraint;		// 0.. local coordinates, 1..local nodes 
	enum {MAX_IS = 100};
};//$EDC$[endclass,AverageConstraint]


#endif  //__AVERAGE_CONSTRAINT__