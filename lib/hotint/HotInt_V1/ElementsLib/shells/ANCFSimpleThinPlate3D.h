#ifndef ANCFSimpleThinPlate3D__H
#define ANCFSimpleThinPlate3D__H

#include "FiniteElement3D.h" 
#include "femathhelperfunctions.h"
#include "Material.h"
#include "node.h"
#include "IntegrationRule.h"

//$EDC$[beginclass,classname=ANCFSimpleThinPlate3D,addelementtype=TAEBody+TAENotInRelease,addelementtypename=ANCFThinPlate3D,parentclassname=Element]
class ANCFSimpleThinPlate3D: public FiniteElementGeneric<Body3D>
{
public:
	ANCFSimpleThinPlate3D(MBS* mbsi) : FiniteElementGeneric<Body3D>(mbsi) 
	{
		ElementDefaultConstructorInitialization();
	}

	ANCFSimpleThinPlate3D(const ANCFSimpleThinPlate3D& e) : FiniteElementGeneric<Body3D>(e.mbs)
	{
		CopyFrom(e);
	}

	virtual void ElementDefaultConstructorInitialization();

	virtual Element* GetCopy()
	{
		Element* e = new ANCFSimpleThinPlate3D(*this);
		return e;
	}
	virtual void CopyFrom(const Element& e)
	{
		FiniteElementGeneric<Body3D>::CopyFrom(e);
		const ANCFSimpleThinPlate3D & ce = (const ANCFSimpleThinPlate3D&)e;
		this->size1 = ce.size1;
		this->size2 = ce.size2;
		this->thickness = ce.thickness;
		this->q_ref = ce.q_ref;
		this->jac = ce.jac;
	}

	virtual void SetANCFSimpleThinPlate3D(int body, const TArray<int>& nodes, int material, 
		double size1, double size2, double thickness, const Vector3D& color);

	virtual const char* GetElementSpec() const { return "ANCFThinPlate3D"; }
	virtual TFiniteElementType GetElementType() const { return TFE_ThinPlate; }
	
	virtual int Dim() const { return 3; }
	virtual int DOFPerNode() const { return 9; }
	virtual int NS() const { return 12; }
	virtual int DataS() const { return 0; }

	virtual void DefineIntegrationRule(IntegrationRule& integrationRule);
	virtual int GetActualInterpolationOrder() const	{ return 4;	}

	virtual double GetThickness() const { return this->thickness; }

	virtual void EvalF2(Vector& f, double t);
	virtual void EvalF2GeomLin(Vector& f, double t);
	virtual void EvalF2GeomNonlin(Vector& f, double t);

	virtual void EvalM(Matrix& m, double t);
	virtual void EvalMff(Matrix& m, double t) { EvalM(m, t); }

	virtual void StiffnessMatrix(Matrix& m);
	virtual void StiffnessMatrixLin(Matrix& m);
	virtual void StiffnessMatrixNonlin(Matrix& m);
	virtual int FastStiffnessMatrix() const { return 1*2; }

	virtual double GetJacInvDS(const Vector3D& p_loc, Matrix& jacinvDS) const;
	virtual Matrix3D GetElementJacobian(const Vector2D& p_loc) const;

	virtual Vector3D GetMidPlaneStrains(const Vector2D& p_loc, TComputeDrawInitFlag flag) const;
	//virtual Vector3D GetMidPlaneStrains(const Vector2D& p_loc, const Matrix& grad, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetCurvatures(const Vector2D& p_loc, TComputeDrawInitFlag flag) const;
	//virtual Vector3D GetCurvatures(const Vector2D& p_loc, const Matrix& grad, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetStrains(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetStresses(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;

	virtual void GetDeltaEps(const Vector2D& p_loc, Matrix& delta_eps, TComputeDrawInitFlag flag) const;
	virtual void GetDeltaKappa(const Vector2D& p_loc, Matrix& delta_kappa, TComputeDrawInitFlag flag) const;

	//virtual void GetDeltaEps(const Vector2D& p_loc, const Matrix& grad, Matrix& delta_eps, TComputeDrawInitFlag flag) const;
	//virtual void GetDeltaKappa(const Vector2D& p_loc, const Matrix& grad, Matrix& delta_kappa, TComputeDrawInitFlag flag) const;

	virtual Vector3D GetPos_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetDisplacement_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetVel_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetRefPos(const Vector3D& p_loc) const { return GetPos_dc(p_loc, TCD_reference_configuration); }
	virtual Vector3D GetRefPosD() const { return GetRefPos(Vector3D(-0.5,-1.,5.)); };

	// redirect to new function
	virtual Vector3D GetPos(const Vector3D& p_loc) const { return GetPos_dc(p_loc, TCD_compute); }
	//virtual Vector3D GetPos(const Vector2D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetPosD(const Vector3D& p_loc) const { return GetPos_dc(p_loc, TCD_draw_magnified); }
	// ref, init, magnified, etc.
	virtual Vector3D GetDisplacementD(const Vector3D& p_loc) const { return GetDisplacement_dc(p_loc, TCD_draw); }

	virtual Vector3D GetNormal(const Vector2D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetNormal(const Vector3D& p_loc, TComputeDrawInitFlag flag) const { return GetNormal((const Vector2D&)p_loc, flag); }

	virtual Vector3D GetDNormalDqi(const Vector2D& p_loc, int dof, TComputeDrawInitFlag flag) const;

	virtual Vector3D GetDPosDAlpha(const Vector2D& p_loc, int alpha, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetDPosPDAlpha(const Vector2D& p_loc, int alpha, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetDDPosDAlphaDBeta(const Vector2D& p_loc, int alpha_beta, TComputeDrawInitFlag flag) const;

	virtual Vector3D GetDPosDAlpha(const Vector2D& p_loc, int alpha, Matrix& grad, TComputeDrawInitFlag flag) const;
	virtual Vector3D GetDDPosDAlphaDBeta(const Vector2D& p_loc, int alpha_beta, Matrix& grad, TComputeDrawInitFlag flag) const;

	virtual Vector3D GetPosx(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDPosDAlpha(p_loc, 1, flag); }
	virtual Vector3D GetPosy(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDPosDAlpha(p_loc, 2, flag); }  
	virtual Vector3D GetPosxP(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDPosPDAlpha(p_loc, 1, flag); }
	virtual Vector3D GetPosyP(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDPosPDAlpha(p_loc, 2, flag); }  

	virtual Vector3D GetPosxx(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDDPosDAlphaDBeta(p_loc, 1, flag); }
	virtual Vector3D GetPosxy(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDDPosDAlphaDBeta(p_loc, 2, flag); }
	virtual Vector3D GetPosyy(const Vector2D& p_loc, TComputeDrawInitFlag flag) const { return GetDDPosDAlphaDBeta(p_loc, 3, flag); }

	virtual Vector3D GetNodeLocPos(int local_node_number) const;

	virtual Matrix3D GetRotMatrix_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Matrix3D GetRotMatrix(const Vector3D& p_loc) const { return GetRotMatrix_dc(p_loc, TCD_compute); }
	virtual Matrix3D GetRotMatrixD(const Vector3D& p_loc) const { return GetRotMatrix_dc(p_loc, TCD_draw); }
	virtual Matrix3D GetRotMatrixP_dc(const Vector3D& p_loc, TComputeDrawInitFlag flag) const;
	virtual Matrix3D GetRotMatrixP(const Vector3D& p_loc) { return GetRotMatrix_dc(p_loc, TCD_compute); }
	
	virtual void GetdPosdqT(const Vector3D& p_loc, Matrix& m);
	virtual void GetdRotdqT(const Vector3D& p_loc, Matrix& m);

	virtual void GetdRotvdqT(const Vector3D& v_loc, const Vector3D& p_loc, Matrix& m);

	// shape functions as they are, accessible per index nShapeFunction = 1 .. 12
	// local coordinates q1, q2 vary from -1 to +1
	double GetS0(const Vector2D& p_loc, int nsf) const;
	double GetDS0(const Vector2D& p_loc, int nsf, int alpha) const;
	double GetDDS0(const Vector2D& p_loc, int nsf, int alpha) const;

	virtual void GetAvailableFieldVariables(TArray<FieldVariableDescriptor>& variables);
	virtual double GetFieldVariableValue_dc(const FieldVariableDescriptor& fvd, const Vector3D& p_loc, TComputeDrawInitFlag flag);
	virtual double GetFieldVariableValue(const FieldVariableDescriptor& fvd, const Vector3D& p_loc, bool flagD)
	{
		if (flagD)
		{
			return GetFieldVariableValue_dc(fvd, p_loc, TCD_draw);
		}
		else
		{
			return GetFieldVariableValue_dc(fvd, p_loc, TCD_compute);
		}
	}

	virtual TKinematicsAccessFunctions GetKinematicsAccessFunctions(int mode = 1) const
	{
		return TKinematicsAccessFunctions(
			TKAF_position + TKAF_displacement + TKAF_velocity + 
			TKAF_rotation_matrix + TKAF_rotation_matrix_P + TKAF_D_rot_v_D_q +
			TKAF_D_pos_D_q + TKAF_D_rot_D_q
			);
	}

	virtual Vector3D GetDOFPosD(int idof) const;
	virtual Vector3D GetDOFDirD(int idof) const;

	virtual void AddSurfacePressure(Vector& f, double pressure, int dir);
	virtual Vector3D GetSurfaceNormalD(int dir);
	

	virtual Box3D GetElementBox() const
	{
		return Box3D(GetPos_dc(Vector3D(-1.,-1.,-1.), TCD_reference_configuration), GetPos_dc(Vector3D(1.,1.,1.),TCD_reference_configuration));
	}

	virtual Box3D GetElementBoxD() const
	{
		return Box3D(GetPos_dc(Vector3D(-1.,-1.,-1.), TCD_reference_configuration), GetPos_dc(Vector3D(1.,1.,1.),TCD_reference_configuration));
	}

	virtual void DrawElement();

	// keep here for now, since q_ref still needs to be implemented in element
	virtual const double XG_dc(int iloc, TComputeDrawInitFlag flag) const 
	{
		switch (flag)
		{
			case TCD_compute:
				return q_ref(iloc) + GetXact(ltg.Get(iloc));
			case TCD_draw:
				return q_ref(iloc) + GetDrawValue(ltg.Get(iloc));
			case TCD_draw_magnified:
			{
				double scaling = mbs->GetDOption(105);
				return q_ref(iloc) + scaling*GetDrawValue(ltg.Get(iloc));
			}
			case TCD_reference_configuration:
				return q_ref(iloc);
			case TCD_initial_values:
				return q_ref(iloc) + GetXInit(iloc);
			case TCD_cached:
				return xg_cached(iloc);			
			default: 
			{
				assert("XG_dc called with illegal flag" && 0);
				return 0.;
			}
		}
	}

	virtual const double XGP_dc(int iloc, TComputeDrawInitFlag flag) const 
	{
		return XG_dc(iloc+SOS(), flag);
	}

	virtual void SetXG_cached(TComputeDrawInitFlag flag = TCD_compute)
	{
		for (int i = 1; i <= 2*SOS(); i++)
		{
			xg_cached.SetLen(2*SOS());
			xg_cached(i) = XG_dc(i, flag);
		}
	}

	virtual void GetElementData(ElementDataContainer& edc)
	{
		Element::GetElementData(edc);
	}
	virtual int SetElementData(ElementDataContainer& edc)
	{
		int rv = Element::SetElementData(edc);
		this->SetANCFSimpleThinPlate3D(bodyind, nodes, materialnum, size1, size2, thickness, col);
		return rv;
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc);
	virtual int SetElementDataAuto(ElementDataContainer& edc);
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata);
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata);
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables);

protected:
	double size1;		//$EDC$[varaccess,EDCvarname="size1",EDCfolder="Geometry",tooltiptext="plate dimension in x-direction"]
	double size2;		//$EDC$[varaccess,EDCvarname="size2",EDCfolder="Geometry",tooltiptext="plate dimension in y-direction"]
	double thickness;	//$EDC$[varaccess,EDCvarname="thickness",EDCfolder="Geometry",tooltiptext="plate thickness"]

	ConstVector<5> jac;
	ConstVector<36*2> q_ref;
	ConstVector<36*2> xg_cached;

	//EDC int nodes(1) //$EDC$[varaccess, EDCvarname="node_1", EDCfolder="Geometry",tooltiptext="node1 of the element"]
	//EDC int nodes(2) //$EDC$[varaccess, EDCvarname="node_2", EDCfolder="Geometry",tooltiptext="node2 of the element"]
	//EDC int nodes(3) //$EDC$[varaccess, EDCvarname="node_3", EDCfolder="Geometry",tooltiptext="node3 of the element"]
	//EDC int nodes(4) //$EDC$[varaccess, EDCvarname="node_4", EDCfolder="Geometry",tooltiptext="node4 of the element"]
	//EDC int materialnum; //$EDC$[varaccess,EDCvarname="material_number",EDCfolder="Physics",tooltiptext="material number containing the material properties of the plate"] 
	//EDC int (int&)geometricNonlinearityStatus; //$EDC$[varaccess,EDCvarname="is_geometrically_nonlinear",EDCfolder="Physics",int_bool,tooltiptext="0 ... geometrically linear (small deformation), 1 ... geometrically nonlinear"] 
};
//$EDC$[endclass,ANCFSimpleThinPlate3D]

#endif