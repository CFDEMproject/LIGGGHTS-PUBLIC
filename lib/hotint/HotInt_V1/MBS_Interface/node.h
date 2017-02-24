//#**************************************************************
//#
//# filename:             node.h
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						July 2004
//# description:          Driver and model for timeintegration
//#                       Model of a rigid arm and a hydraulic zylinder
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
//#**************************************************************


#ifndef NODE__H
#define NODE__H

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                   NODE
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Vector3D node_default_col(0.4,0.4,0.1);

//$EDC$[beginclass,classname=BaseNode]
class BaseNode
{
public:
	virtual Node* GetCopy() = 0;
	virtual void CopyFrom(const Node& n) = 0;
	virtual const mystr& GetObjectName() const = 0;
	virtual mystr& GetObjectName() = 0;
	virtual void SetObjectName(const char* name) = 0;

	virtual mystr GetTypeName() {return "Node";};	//$EDC$[funcaccess,EDCvarname="node_type",tooltiptext="specification of node type. Once the node is added to the mbs, you MUST NOT change this type anymore!"]

	virtual const MBS* GetMBS() const = 0;
	virtual MBS* GetMBS() = 0;
	virtual void SetMBS(MBS* mbsI) = 0;
	virtual UserOutputInterface & UO(int message_level = UO_LVL_all, int output_prec = -1) = 0;

	virtual const int& NodeNum() const = 0;
	virtual int& NodeNum() = 0;

	virtual const char* GetObjectSpec() const = 0;

	virtual void GetElementData(ElementDataContainer& edc) = 0;
	virtual int SetElementData(ElementDataContainer& edc) = 0;

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void Set(int i, int val) = 0;
	virtual void AddLTG(int val) = 0;
	virtual int Get(int i) const = 0;
	virtual int Length() const = 0;
	virtual int LTGLength() const = 0;
	virtual void LTGreset() = 0;

	virtual const int& LTG(int val) const = 0;
	virtual int& LTG(int val) = 0;

	virtual const double& XG(int iloc) const = 0;
	virtual double& XG(int iloc) = 0;
	virtual const double& XGP(int iloc) const = 0;
	virtual double& XGP(int iloc) = 0;
	virtual const double& XGD(int iloc) const = 0;
	virtual const double& XGPD(int iloc) const = 0;

	virtual const int& SOS() const = 0;
	virtual int& SOS() = 0;

	virtual const int& ES() const = 0;
	virtual int& ES() = 0;
	virtual const int& IS() const = 0;
	virtual int& IS() = 0;

	virtual int GetBodyInd() const = 0;
	virtual void SetBodyInd(int bodyind) = 0;
	virtual const Vector3D& Pos() const = 0;
	virtual Vector3D& Pos() = 0;
	virtual Vector& X_Init() = 0;
	virtual const Vector& X_Init() const = 0;
	virtual int SetX_Init(Vector3D pos, Vector3D vel) = 0;
	virtual int SetX_Init(Vector2D pos, Vector2D vel) = 0;
	virtual int SetX_Init(const Vector& v) = 0;

	virtual double GetDrawTemp() const = 0;
	virtual double& GetDrawTemp() = 0;
	virtual void SetDrawTemp(double draw_tempI) = 0;

	virtual int GetDrawTempCnt() const = 0;
	virtual int& GetDrawTempCnt() = 0;
	virtual void SetDrawTempCnt(int draw_temp_cntI) = 0;

	virtual Vector3D RefConfPos() const = 0;
	virtual Vector2D RefConfPos2D() const = 0;

	virtual void ResetNodeToElementList() = 0;
	virtual int AddElementNumber(int elnum) = 0;
	virtual int NElements() const = 0;
	virtual int ElemNum(int i) const = 0;

	virtual Vector3D GetDisplacement() const = 0;
	virtual Vector2D GetDisplacement2D() const = 0;
	virtual Vector3D GetDisplacementD() const = 0;
	virtual Vector2D GetDisplacement2DD() const = 0;
	virtual Vector3D GetVel() const = 0;
	virtual Vector2D GetVel2D() const = 0;
	virtual Vector3D GetVelD() const = 0;
	virtual Vector2D GetVel2DD() const = 0;
	virtual Vector3D GetPos() const = 0;
	virtual Vector2D GetPos2D() const = 0;
	virtual Vector3D GetPosD() const = 0;
	virtual Vector2D GetPos2DD() const = 0;
	virtual Vector3D ToP3D(const Vector2D& p) const = 0;

	virtual int Dim() const = 0;
	virtual const Vector3D& GetCol() const = 0;
	virtual void DrawNode() = 0;
	virtual int IsVisible() const = 0;
	virtual int IsAuxNode() const = 0;
	virtual void SetAuxNode(int flag=1) = 0;

protected:
	BaseNode() : ltg(0), x_init(), elements(), sos(0), es(0), is(0) {}
	BaseNode(int ndof) : ltg(2*ndof), x_init(), elements(), sos(0), es(0), is(0) {}

	mystr nodename;			//$EDC$[varaccess,EDCvarname="name",EDCfolder="",tooltiptext="Node identifier."]
	int nodenum;			  //$EDC$[varaccess,EDCvarname="node_number",EDCfolder="",tooltiptext="Node Number."]
	int body;				    //should be the domain number in the future //EDC$[varaccess,EDCvarname="body_index",EDCfolder="",tooltiptext="Body index, where node belongs to."]

	int isauxnode;			//if node belongs to CMS element, it is not really contain DOFs!

	TArray<int> ltg;		  //local to global DOF list
	TArray<int> elements;	//node to element list
	int sos;            //EDC$[varaccess,EDCvarname="sos",EDCfolder="",tooltiptext="Number of second order degrees of freedom."]
	int es;
	int is;
	

	Vector3D pos;			  //$EDC$[varaccess,EDCvarname="reference_position",EDCfolder="Geometry",tooltiptext="Position (2D/3D) in reference configuration."]

	MBS* mbs;
	Vector x_init;			//$EDC$[varaccess,EDCvarname="node_initial_values",EDCfolder="Initialization",vecstart=1,vecend=2*SOS(),tooltiptext="initial values for all degrees of freedom of node"]
	//initial conditions (displacements)
							//*YV: x_init may also contain reference configuration,
							// when more that 3 degrees of freedom are assiciated with a node (see ANCFThinPlate3D)

	//drawing:
	Vector3D col;        //$EDC$[varaccess,EDCvarname="RGB_color",EDCfolder="Graphics",tooltiptext="[red, green, blue] color of element, range = 0..1,"]
	int visible;        //$EDC$[varaccess,EDCvarname="visible",EDCfolder="Graphics",tooltiptext="Visibility of node."]
	double draw_temp;		  //temporary storage for drawing (stress, etc.)
	int draw_temp_cnt;		//temporary counter for drawing, number of elements per node
};
//$EDC$[endclass,BaseNode]

//$EDC$[beginclass,classname=Node,parentclassname=BaseNode]
class Node : public BaseNode
{
public:
	Node() : BaseNode() //: ltg(0), x_init(), elements()
	{
		InitConstructor();
	};
	Node(MBS* mbs) : BaseNode()
	{
		this->mbs = mbs;
		InitConstructor();
	}
	void InitConstructor()
	{
		x_init = Vector(2*SOS());
		ltg.SetLen(0);
		elements.SetLen(0);
		sos=0; 
		pos = Vector3D(0,0,0); 
		body = 0; 
		nodename = mystr(GetObjectSpec()); 
		nodenum = 0;
		mbs = 0;
		col = node_default_col;
		visible = 1;
		isauxnode = 0;
	}
	Node(int ndof) : BaseNode(ndof) //: ltg(2*ndof), x_init(), elements()
	{
		InitConstructor();
		sos = ndof;
	};
	Node(int ndof, int bodyind, const Vector3D& p) : BaseNode(ndof) //ltg(2*ndof), x_init(), elements()
	{
		InitConstructor();
		sos = ndof; 
		body = bodyind; 
		pos = p;
	};
	Node(const Node& n) : BaseNode(n.ltg.Length()) //ltg(n.ltg.Length()), x_init(), elements()
	{
		CopyFrom(n);
	};
	//$ AH 2013-01-20: not needed any longer
	/*Node& operator=(const Node& n) 
	{
		if (this == &n) {return *this;}
		CopyFrom(n);
		return *this;
	}*/
	virtual Node* GetCopy()
	{
		Node* nc = new Node();
		nc->CopyFrom(*this);
		return nc;
	}
	virtual void CopyFrom(const Node& n)
	{
		ltg.SetLen(n.Length());
		for (int i=1; i <= n.Length(); i++)
			ltg(i) = n.ltg(i);
		sos = n.sos; es = n.es; is = n.is;

		body = n.body;
		pos = n.pos;
		nodename = n.nodename;
		nodenum = n.nodenum;

		mbs = n.mbs;
		x_init = n.x_init;
		elements = n.elements;

		col = n.col;
		visible = n.visible;

		isauxnode = n.isauxnode;
	}

	virtual ~Node() {ltg.Flush();};

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const mystr& GetObjectName() const {return nodename;}
	virtual mystr& GetObjectName() {return nodename;}
	virtual void SetObjectName(const char* name) {nodename = name;}

	virtual const MBS* GetMBS() const {return mbs;}
	virtual MBS* GetMBS() {return mbs;}
	virtual void SetMBS(MBS* mbsI) {mbs = mbsI;}
	//virtual UserOutput& UO() {return mbs->uout;};
	//virtual UserOutput& UO(int message_level = UO_LVL_all) { mbs->uout.SetLocalMessageLevel(message_level); return mbs->uout;}; //(AD)
	//$ YV 2012-11-28
	//virtual UserOutput& UO(int message_level = UO_LVL_all, int output_prec = -1) { mbs->uout.SetLocalMessageLevel(message_level); mbs->uout.SetOutputPrec(output_prec); return mbs->uout;}; //$ AD 2011-02 output_prec
	virtual UserOutputInterface & UO(int message_level = UO_LVL_all, int output_prec = -1) { return mbs->UO(message_level,output_prec); }

	virtual const int& NodeNum() const {return nodenum;}
	virtual int& NodeNum() {return nodenum;}

	virtual const char* GetObjectSpec() const {return "Node";}	//identifier of object

	virtual void GetElementData(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!


	virtual void Set(int i, int val) {ltg.Elem(i) = val;} //index
	virtual void AddLTG(int val) {ltg.Add(val);}
// (AD) changed () to .Get()
	virtual int Get(int i) const {return ltg.Get(i);} //index
//	virtual int Get(int i) const {return ltg(i);} //index
	virtual int Length() const {return ltg.Length();}
	virtual int LTGLength() const {return ltg.Length();} //=Length()
	virtual void LTGreset() {ltg.SetLen(0);}

// (AD) changed () to .Get() & Elem()
	virtual const int& LTG(int val) const {return ltg.Get(val);} //=Get(int i)
	virtual int& LTG(int val) {return ltg.Elem(val);} 
//	virtual const int& LTG(int val) const {return ltg(val);} //=Get(int i)
//	virtual int& LTG(int val) {return ltg(val);} 

// (AD) changed () to .Get()
	virtual const double& XG(int iloc) const 
	{
		return GetMBS()->GetXact(ltg.Get(iloc));
	}
	virtual double& XG(int iloc) 
	{
		return GetMBS()->GetXact(ltg.Get(iloc));
	}
	virtual const double& XGP(int iloc) const 
	{
		return GetMBS()->GetXact(ltg.Get(iloc+SOS()));
	}
	virtual double& XGP(int iloc) 
	{
		return GetMBS()->GetXact(ltg.Get(iloc+SOS()));
	}
	virtual const double& XGD(int iloc) const 
	{
		return GetMBS()->GetDrawValue(ltg.Get(iloc));
	}

	virtual const double& XGPD(int iloc) const 
	{
		return GetMBS()->GetDrawValue(ltg.Get(iloc+SOS()));
	}

	virtual const int& SOS() const {return sos;}
	virtual int& SOS() {return sos;}

	//$ AH 2013-01-20: additionally needed by MBS
	virtual const int& ES() const {return es;}
	virtual int& ES() {return es;}
	virtual const int& IS() const {return is;}
	virtual int& IS() {return is;}

	virtual int GetBodyInd() const {return body;}
	virtual void SetBodyInd(int bodyind) {body = bodyind;};
	virtual const Vector3D& Pos() const {return pos;}
	virtual Vector3D& Pos() {return pos;}
	virtual Vector& X_Init() {return x_init;}
	virtual const Vector& X_Init() const {return x_init;}
	virtual int SetX_Init(Vector3D pos, Vector3D vel)
	{
		Vector vec;
		vec.SetLen(6);
		vec(1) = pos.X();
		vec(2) = pos.Y();
		vec(3) = pos.Z();
		vec(4) = vel.X();
		vec(5) = vel.Y();
		vec(6) = vel.Z();
		return SetX_Init(vec);
	}
	virtual int SetX_Init(Vector2D pos, Vector2D vel)
	{
		Vector vec;
		vec.SetLen(4);
		vec(1) = pos.X();
		vec(2) = pos.Y();
		vec(3) = vel.X();
		vec(4) = vel.Y();
		return SetX_Init(vec);
	}
	virtual int SetX_Init(const Vector& v)
	{
		x_init = v;
		return 0;
	}

	virtual double GetDrawTemp() const {return draw_temp;}
	virtual double& GetDrawTemp() {return draw_temp;}
	virtual void SetDrawTemp(double draw_tempI) {draw_temp = draw_tempI;}

	virtual int GetDrawTempCnt() const {return draw_temp_cnt;}
	virtual int& GetDrawTempCnt() {return draw_temp_cnt;}
	virtual void SetDrawTempCnt(int draw_temp_cntI) {draw_temp_cnt = draw_temp_cntI;}

	virtual Vector3D RefConfPos() const {if (Dim()==2) return ToP3D(RefConfPos2D()); else return pos;}
	virtual Vector2D RefConfPos2D() const {return Vector2D(pos.X(),pos.Y());}

	virtual void ResetNodeToElementList() {elements.SetLen(0);}
	virtual int AddElementNumber(int elnum) {return elements.Add(elnum);}
	virtual int NElements() const {return elements.Length();}
	virtual int ElemNum(int i) const {return elements(i);}

	virtual Vector3D GetDisplacement() const
	{
		if (Dim()==2) return ToP3D(GetDisplacement2D()); else return Vector3D(XG(1),XG(2),XG(3));
	}

	virtual Vector2D GetDisplacement2D() const
	{
		return Vector2D(XG(1),XG(2));
	}

	virtual Vector3D GetDisplacementD() const
	{
		if (Dim()==2)
		{
			return ToP3D(GetDisplacement2DD());
		}

		double fact = GetMBS()->GetDOption(105); //deformation scaling
		return fact * Vector3D(XGD(1),XGD(2),XGD(3)); //$ AD fix node displacement scale bug
	}

	virtual Vector2D GetDisplacement2DD() const
	{
		double fact = GetMBS()->GetDOption(105); //deformation scaling
		return fact * Vector2D(XGD(1),XGD(2)); //$ AD fix node displacement scale bug
	}

	virtual Vector3D GetVel() const
	{
		if (Dim()==2) return ToP3D(GetVel2D()); else return Vector3D(XGP(1),XGP(2),XGP(3));
	}

	virtual Vector2D GetVel2D() const 
	{
		return Vector2D(XGP(1),XGP(2));
	}

	virtual Vector3D GetVelD() const 
	{
		if (Dim()==2) 
		{
			return ToP3D(GetVel2DD());
		}
		
		return Vector3D(XGPD(1),XGPD(2),XGPD(3));
	}

	virtual Vector2D GetVel2DD() const
	{
		return Vector2D(XGPD(1),XGPD(2));
	}

	virtual Vector3D GetPos() const {return RefConfPos() + GetDisplacement();}
	virtual Vector2D GetPos2D() const {return RefConfPos2D() + GetDisplacement2D();}
	virtual Vector3D GetPosD() const 
	{


		return RefConfPos() + GetDisplacementD();
	}
	virtual Vector2D GetPos2DD() const {return RefConfPos2D() + GetDisplacement2DD();}


	virtual Vector3D ToP3D(const Vector2D& p) const 
	{
		return Vector3D(p.X(),p.Y(),0.);
	}

	virtual int Dim() const;


	virtual const Vector3D& GetCol() const {return col;}

	virtual void DrawNode();

	virtual int IsVisible() const {return visible != 0;}

	virtual int IsAuxNode() const {return isauxnode;}
	virtual void SetAuxNode(int flag=1) {isauxnode = flag;}


/*
private:
	mystr nodename;				//node identifier
	int nodenum;					//node number
	int isauxnode;				//if node belongs to CMS element, it is not really contain DOFs!

	double draw_temp;			//temporary storage for drawing (stress, etc.)
	int draw_temp_cnt;		//temporary counter for drawing, number of elements per node

	TArray<int> ltg; //local to global DOF list
	int sos;							//number of second order degrees of freedom
	int body;							//body where node belongs to

	Vector3D pos;					//Reference configuration position (2D/3D)

	MBS* mbs;
	Vector x_init;				//initial conditions (displacements)
												//*YV: x_init may also contain reference configuration,
												// when more that 3 degrees of freedom are assiciated with a node (see ANCFThinPlate3D)

	TArray<int> elements; //node to element list

	//drawing:
	Vector3D col;
	int visible;					//flag, if visible ==0 ==> do not draw, else draw node
*/
};
//$EDC$[endclass,Node]

typedef enum { NDT_position = 1, NDT_displacement = 2, NDT_angle = 4, NDT_slope1 = 8, NDT_slope2 = 16, NDT_slope3 = 32, 
	NDT_temperature = 64, NDT_electric_potential = 128, NDT_velocity = 256, 
	NDT_comp_1 = 512, NDT_comp_2 = 1024, NDT_comp_3 = 2048, 
	NDT_sos = 4096, NDT_es = 4096*2, NDT_is = 4096*2*2} NodeDOFType;

//$EDC$[beginclass,classname=GenericNode,parentclassname=Node]
class GenericNode : public Node
{
public:
	GenericNode() {}
	GenericNode(MBS* mbs) {}
	GenericNode(const GenericNode& n) //: ltg(), x_init(), elements(), doftypes()
	{
		CopyFrom(n);
	}
	virtual GenericNode* GetCopy() 
	{
		GenericNode* n = new GenericNode(*this);
		return n;
	}

	void InitConstructor() 
	{
		x_init = Vector(2*SOS());
	}

	virtual void CopyFrom(const GenericNode& n)
	{
		body = n.body;
		//ref_pos = n.ref_pos;
		nodename = n.nodename;
		nodenum = n.nodenum;
		mbs = n.mbs;
		sos = n.sos; es = n.es; is = n.is;

		pos = n.pos;
		ltg = n.ltg;
		
		elements = n.elements;
		doftypes = n.doftypes;

		x_init = n.x_init;

		col = n.col;
		visible = n.visible;
		isauxnode = n.isauxnode;
	}
	virtual ~GenericNode() { ltg.Flush(); doftypes.Flush(); elements.Flush(); };

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const char* GetObjectSpec() const { return "GenericNode"; }	//identifier of object
	virtual const int& NodeNum() const { return nodenum; }
	virtual int& NodeNum() { return nodenum; }

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual void Set(int i, int val) { ltg.Elem(i) = val; }	//index
	virtual int Get(int i) const { return ltg.Get(i); }		//index
	virtual void AddLTG(int val) { ltg.Add(val); }
	virtual int LTGLength() const { return ltg.Length(); }	
	virtual void LTGreset() { ltg.SetLen(0); }
	virtual const int& LTG(int val) const { return ltg.Get(val); } //=Get(int i)
	virtual int& LTG(int val) { return ltg.Elem(val); }
	virtual int Length() const { return ltg.Length(); }

	virtual const int& SOS() const { return sos; }
	virtual int& SOS() { return sos; }
	virtual const int& ES() const {return es;}
	virtual int& ES() {return es;}
	virtual const int& IS() const {return is;}
	virtual int& IS() {return is;}

	virtual int Dim() const;
	virtual void SetInitValues(const Vector& x_init) 
	{
		if (x_init.Length() == 2*sos+es+is)
			this->x_init = x_init;
		else
			mbs->UO() << "inconsistent inital values of Node " << this->NodeNum() << "\n";
	}

	virtual const double& XG(int iloc) const
	{ 
		return mbs->GetXact(ltg.Get(iloc)); 
	}
	virtual double& XG(int iloc) 
	{ 
		return mbs->GetXact(ltg.Get(iloc));
	}
	virtual const double& XGP(int iloc) const 
	{ 
		return mbs->GetXact(ltg.Get(iloc+SOS())); 
	}
	virtual double& XGP(int iloc) 
	{ 
		return mbs->GetXact(ltg.Get(iloc+SOS())); 
	}
	virtual const double& XGD(int iloc) const 
	{
		return mbs->GetDrawValue(ltg.Get(iloc)); 
	}
	virtual const double& XGPD(int iloc) const 
	{ 
		return mbs->GetDrawValue(ltg.Get(iloc+SOS())); 
	}

	virtual Vector3D GetRefPos() const { return pos; }
	virtual Vector2D GetRefPos2D() const { return Vector2D(pos(1), pos(2)); }

	virtual Vector3D GetPos() const { return GetRefPos() + GetDisplacement(); }
	virtual Vector2D GetPos2D() const { return GetRefPos2D() + GetDisplacement2D(); }
	virtual Vector3D GetPosD() const { return GetRefPos() + GetDisplacementD(); }
	virtual Vector2D GetPos2DD() const { return GetRefPos2D() + GetDisplacement2DD(); }

	virtual Vector3D GetDisplacement() const;
	virtual Vector2D GetDisplacement2D() const;
	virtual Vector3D GetDisplacementD() const; 
	virtual Vector2D GetDisplacement2DD() const; 

	virtual Vector3D GetVel() const; 
	virtual Vector2D GetVel2D() const;
	virtual Vector3D GetVelD() const;
	virtual Vector2D GetVel2DD() const;

	virtual int AddElementNumber(int elnum) { return elements.Add(elnum); }
	virtual int NElements() const { return elements.Length(); }
	virtual int ElemNum(int i) const { return elements(i); }
	virtual void ResetNodeToElementList() { elements.SetLen(0); }

	virtual int GetBodyInd() const { return body; }
	virtual void SetBodyInd(int bodyind) { body = bodyind; };

	// draw routines
	virtual double GetDrawTemp() const { return draw_temp; }
	virtual double& GetDrawTemp() { return draw_temp; }
	virtual void SetDrawTemp(double draw_temp) { this->draw_temp = draw_temp; }

	virtual int GetDrawTempCnt() const { return draw_temp_cnt; }
	virtual int& GetDrawTempCnt() { return draw_temp_cnt; }
	virtual void SetDrawTempCnt(int draw_temp_cnt) { this->draw_temp_cnt = draw_temp_cnt; }

	virtual void DrawNode();

	virtual int AddNodeDOF(NodeDOFType ndt);

protected:
	/*
	MBS* mbs;
	int nodenum;
	TArray<int> ltg;		//local to global DOF list
	TArray<int> elements;
	TArray<NodeDOFType> doftypes;	

	int sos, es, is;
	int body;	// check again, if this is necessary

	Vector3D ref_pos;
	Vector x_init;

	Vector3D color;
	double draw_temp;		//temporary storage for drawing (stress, etc.)
	int draw_temp_cnt;		//temporary counter for drawing, number of elements per node
	*/
	TArray<NodeDOFType> doftypes;
	//int es, is;
	//Vector3D ref_pos;  //this exists already in the base class as "pos" //EDC$[varaccess,EDCvarname="Pos",EDCfolder="Geometry",tooltiptext="Position of node in reference configuration."]
};
//$EDC$[endclass,GenericNode]

//$EDC$[beginclass,classname=Node2D,parentclassname=GenericNode]
class Node2D : public GenericNode
{
public: 
	Node2D(MBS* mbs) : GenericNode(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
		InitConstructor();
	}

	Node2D(const Node2D& n) : GenericNode(n.mbs) { CopyFrom(*this); }

	virtual GenericNode* GetCopy() 
	{
		GenericNode* n = new Node2D(*this);
		return n;
	}
	~Node2D(void) {};
	virtual void CopyFrom(const Node2D& n)
	{
		GenericNode::CopyFrom(n);
	}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual int Dim() const { return 2; }

	virtual Vector3D GetDisplacement() const { return Vector3D(XG(1), XG(2), 0.); } //JG2013
	virtual Vector2D GetDisplacement2D() const { return Vector2D(XG(1), XG(2)); }
	virtual Vector3D GetDisplacementD() const 
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector3D(0.);
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector3D(XGD(1), XGD(2), 0.); //JG2013
	}
	virtual Vector2D GetDisplacement2DD() const
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector2D(0.);
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector2D(XGD(1), XGD(2)); 
	} 

	virtual Vector3D GetVel() const	
	{ 
		return Vector3D(XGP(1), XGP(2), 0.); 
	}
	virtual Vector2D GetVel2D() const 
	{
		return Vector2D(XGP(1), XGP(2)); 
	}
	virtual Vector3D GetVelD() const 
	{
		return Vector3D(XGPD(1), XGPD(2), 0.); 
	}
	virtual Vector2D GetVel2DD() const 
	{ 
		return Vector2D(XGPD(1), XGPD(2)); 
	}
};
//$EDC$[endclass,Node2D]

//$EDC$[beginclass,classname=Node3D,parentclassname=GenericNode,addelementtypename=Node3D,
//texdescription="Node3D is the basic finite element node in 3D. It owns a reference position in 3D, and 3 degrees of freedom resembling the displacement in 3D.",
//texdescriptionGeometry="The geometry of the node is defined by its current position $\mathbf{r}$ measured in the global frame of the multibody system, which is the sum of the user defined reference position $\mathbf{r}_0$ and the displacement vector $\mathbf{u}$, which is composed of the nodal degrees of freedom.",
//texdescriptionDOF="This node provides three degrees of freedom, all of which are components of the displacement vector $\mathbf{u}=(q_1,q_2,q_3)^T$ measured in the global frame of the multibody system.",
//example="Node3D.txt"]
class Node3D : public GenericNode
{
public: 
	Node3D() : GenericNode() {}
	Node3D(MBS* mbs) : GenericNode(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_3 + NDT_sos));
		InitConstructor();
	}


	Node3D(const Node3D& n) : GenericNode(n.mbs) { CopyFrom(*this); }
	virtual Node3D* GetCopy() 
	{
		Node3D* n = new Node3D();
		n->CopyFrom(*this);
		return n;
	}
	~Node3D(void) {};
	virtual void CopyFrom(const Node3D& n)
	{
		GenericNode::CopyFrom(n);
	}
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual mystr GetTypeName() {return "Node3D";};	
	virtual mystr GetObjectSpec() {return "Node3D";}
	virtual int Dim() const { return 3; }

	virtual Vector3D GetDisplacement() const { return Vector3D(XG(1), XG(2), XG(3)); }
	virtual Vector2D GetDisplacement2D() const { return Vector2D(XG(1), XG(2)); }
	virtual Vector3D GetDisplacementD() const 
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector3D(0.);
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector3D(XGD(1), XGD(2), XGD(3)); 
	}
	virtual Vector2D GetDisplacement2DD() const
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector2D();
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector2D(XGD(1), XGD(2)); 
	} 

	virtual Vector3D GetVel() const	
	{
		return Vector3D(XGP(1), XGP(2), XGP(3)); 
	}
	virtual Vector2D GetVel2D() const
	{
		return Vector2D(XGP(1), XGP(2)); 
	}
	virtual Vector3D GetVelD() const 
	{
		return Vector3D(XGPD(1), XGPD(2), XGPD(3)); 
	}
	virtual Vector2D GetVel2DD() const
	{ 
		return Vector2D(XGPD(1), XGPD(2)); 
	}
};
//$EDC$[endclass,Node3D]

//$EDC$[beginclass,classname=Node3DRxyz,parentclassname=Node3D,addelementtypename=Node3DRxyz,
//texdescription="Node3DRxyz is a finite element node in 3D. It has a 3D reference position, a reference orientation described by bryant angles and 6 degrees of freedom.",
//texdescriptionDOF="The first 3 degrees of freedom are used to describe the displacement $(q_1,q_2,q_3)^T = \mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the last 3 are used for the description of linearized (small) angles $(\phi_x,\phi_y,\phi_z)^T$. All degrees of freedom are w.r.t. the global coordinate system.",
//texdescriptionGeometry="The reference position of the node is defined by the user via $\mathtt{Geometry.reference\_position}$ and the reference orientation via $\mathtt{Geometry.reference\_rot\_angles}$. The current position is evaluated by adding the displacement (the first three degrees of freedom $\left(q_1,q_2,q_3\right)^T$) to the reference position of the node."]
class Node3DRxyz : public Node3D
{
public: 
	Node3DRxyz() : Node3D() {}
	Node3DRxyz(MBS* mbs) : Node3D(mbs)
	{
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_3 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_3 + NDT_sos));
		InitConstructor();
	}
	void InitConstructor() 
	{
		GenericNode::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_angles = Vector3D(0.);
	}

	Node3DRxyz(const Node3DRxyz& n) : Node3D(n.mbs) { CopyFrom(*this); }

	virtual Node3DRxyz* GetCopy() 
	{
		Node3DRxyz* n = new Node3DRxyz();
		n->CopyFrom(*this);
		return n;
	}
	~Node3DRxyz(void) {};
	virtual void CopyFrom(const Node3DRxyz& n)
	{
		GenericNode::CopyFrom(n);
		ref_angles = n.ref_angles;
	}

	virtual void Set(Vector3D posI, Vector3D ref_anglesI)
	{
		pos = posI;
		ref_angles = ref_anglesI;
	}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const Vector3D& GetRefAngles() const { return ref_angles; }

	virtual Matrix3D GetLocalFrame() const
	{ 
		return ComputeRotMatrixWithKardanAngles(ref_angles.X(), ref_angles.Y(), ref_angles.Z());
	}

	virtual double GetLocalFrame(int i, int j) const
	{
		//$ PG 2013-1-21: quick solution for now... has to be done more efficiently
		return GetLocalFrame()(i,j); 
	}

	virtual mystr GetTypeName() {return "Node3DRxyz";};	
	virtual mystr GetObjectSpec() {return "Node3DRxyz";}
	virtual int Dim() const { return 6; }

protected:
	Vector3D ref_angles;     //$EDC$[varaccess,EDCvarname="reference_rot_angles",EDCfolder="Geometry",tooltiptext="Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."]
};
//$EDC$[endclass,Node3DRxyz]

//$EDC$[beginclass,classname=Node3DR123,parentclassname=Node3D,addelementtypename=Node3DR123,
//texdescription="Node3DR123 is a finite element node in 3D. It has a 3D reference position, a reference orientation described by bryant angles and 6 degrees of freedom.",
//texdescriptionDOF="The first 3 degrees of freedom are used to describe the displacement $(q_1,q_2,q_3)^T = \mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the last 3 are used for the description of linearized (small) angles $(\phi_x,\phi_y,\phi_z)^T$. All degrees of freedom are w.r.t. the reference coordinate system of the node.",
//texdescriptionGeometry="The reference position of the node is defined by the user via $\mathtt{Geometry.reference\_position}$ and the orientation via $\mathtt{Geometry.reference\_rot\_angles}$. The current position is evaluated by adding the displacement (the first three degrees of freedom $\left(q_1,q_2,q_3\right)^T$ transformed into the global coordinate system) to the reference position of the node."]
class Node3DR123 : public Node3D
{
public: 
	Node3DR123() : Node3D() {}
	Node3DR123(MBS* mbs) : Node3D(mbs)
	{
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
		//AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_3 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_3 + NDT_sos));
		InitConstructor();
	}
	void InitConstructor() 
	{
		GenericNode::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_angles = Vector3D(0.);
	}

	Node3DR123(const Node3DR123& n) : Node3D(n.mbs) { CopyFrom(*this); }

	virtual Node3DR123* GetCopy() 
	{
		Node3DR123* n = new Node3DR123();
		n->CopyFrom(*this);
		return n;
	}
	~Node3DR123(void) {};
	virtual void CopyFrom(const Node3DR123& n)
	{
		GenericNode::CopyFrom(n);
		ref_angles = n.ref_angles;
	}

	virtual void Set(Vector3D posI, Vector3D ref_anglesI)
	{
		pos = posI;
		ref_angles = ref_anglesI;
	}

	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual const Vector3D& GetRefAngles() const { return ref_angles; }

	virtual Matrix3D GetLocalFrame() const
	{ 
		return ComputeRotMatrixWithKardanAngles(ref_angles.X(), ref_angles.Y(), ref_angles.Z());
	}

	virtual double GetLocalFrame(int i, int j) const
	{
		//$ PG 2013-1-21: quick solution for now... has to be done more efficiently
		return GetLocalFrame()(i,j); 
	}

	virtual mystr GetTypeName() {return "Node3DR123";};	
	virtual mystr GetObjectSpec() {return "Node3DR123";}
	virtual int Dim() const { return 6; }

	virtual Vector3D GetDisplacement() const { return GetLocalFrame()*Vector3D(XG(1), XG(2), XG(3)); }
	virtual Vector2D GetDisplacement2D() const
	{
		Vector3D tmp = GetLocalFrame()*Vector3D(XG(1), XG(2), XG(3));
		return Vector2D(tmp(1), tmp(2)); 
	}
	virtual Vector3D GetDisplacementD() const 
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector3D(0.);
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * GetLocalFrame()*Vector3D(XGD(1), XGD(2), XGD(3)); 
	}
	virtual Vector2D GetDisplacement2DD() const
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector2D();
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		Vector3D tmp = GetLocalFrame()*Vector3D(XGD(1), XGD(2), XGD(3));
		return fact * Vector2D(tmp(1), tmp(2)); 
	} 

	virtual Vector3D GetVel() const	
	{
		return GetLocalFrame()*Vector3D(XGP(1), XGP(2), XGP(3)); 
	}
	virtual Vector2D GetVel2D() const
	{
		Vector3D tmp = GetLocalFrame()*Vector3D(XGP(1), XGP(2), XGP(3));
		return Vector2D(tmp(1), tmp(2)); 
	}
	virtual Vector3D GetVelD() const 
	{
		return GetLocalFrame()*Vector3D(XGPD(1), XGPD(2), XGPD(3)); 
	}
	virtual Vector2D GetVel2DD() const
	{ 
		Vector3D tmp = GetLocalFrame()*Vector3D(XGPD(1), XGPD(2), XGPD(3));
		return Vector2D(tmp(1), tmp(2)); 
	}

protected:
	Vector3D ref_angles;     //$EDC$[varaccess,EDCvarname="reference_rot_angles",EDCfolder="Geometry",tooltiptext="Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."]
};
//$EDC$[endclass,Node3DR123]


//!EDC$[beginclass,classname=NodeBeam3D,parentclassname=GenericNode,addelementtypename=NodeLinearBeam3DTest,
//texdescription="NodeLinearBeam3DTest is a temporary node and it should be used with LinearBeam3D element. It provides $6$ degrees of freedom. $3$ DOF are for translation, $3$ DOF for rotation.",
//example="NodeLinearBeam3D_test.txt"]
//class NodeBeam3D : public GenericNode
//{
//public:
//	NodeBeam3D():GenericNode() {}
//
//	NodeBeam3D(MBS* mbs) : GenericNode(mbs)
//	{
//		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
//		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_3 + NDT_sos));
//		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_3 + NDT_sos));
//		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_2 + NDT_sos));
//		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_1 + NDT_sos));
//		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
//		InitConstructor();
//	}
//
//	void InitConstructor() 
//	{
//		GenericNode::InitConstructor();
//		nodename = mystr(GetObjectSpec()); 
//		ref_angles = Vector3D(0.);
//	}
//
//	NodeBeam3D(const NodeBeam3D& n) : GenericNode(n.mbs) { CopyFrom(*this); }
//	virtual NodeBeam3D* GetCopy() 
//	{
//		NodeBeam3D* n = new NodeBeam3D();
//		n->CopyFrom(*this);
//		return n;
//	}
//	~NodeBeam3D() {};
//
//	virtual void CopyFrom(const NodeBeam3D& n)
//	{
//		GenericNode::CopyFrom(n);
//		ref_angles = n.ref_angles;
//	}
//
//	virtual int Dim() const { return 3; }
//
//	void SetNodeBeam3D(const Vector3D& ref_pos, const Vector3D& ref_kardan_angles, int bodyind = 0) //ltg(2*ndof), x_init(), elements()
//	{
//		//sos = ndof; 
//		body = bodyind; 
//		pos = ref_pos;
//		ref_angles = ref_kardan_angles;
//	};
//	
//	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
//	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
//	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
//	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
//	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter
//
//	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
//	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
//
//	virtual const Vector3D& GetRefAngles() const { return ref_angles; }
//
//	virtual Matrix3D GetLocalFrame() const
//	{ 
//		return ComputeRotMatrixWithKardanAngles(ref_angles.X(), ref_angles.Y(), ref_angles.Z());
//	}
//
//	virtual double GetLocalFrame(int i, int j) const
//	{
//		//$ PG 2013-1-21: quick solution for now... has to be done more efficiently
//		return GetLocalFrame()(i,j); 
//	}
//
//	virtual mystr GetTypeName() {return "NodeLinearBeam3DTest";}	
//	virtual mystr GetObjectSpec() {return "NodeLinearBeam3DTest";}
//
//	virtual Vector3D GetDisplacement() const { return Vector3D(XG(6), XG(1), XG(3)); }
//	virtual Vector2D GetDisplacement2D() const { return Vector2D(XG(6), XG(1)); }
//	virtual Vector3D GetDisplacementD() const 
//	{
//		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
//		{
//			return Vector3D(0.);
//		}
//
//		double fact = mbs->GetDOption(105);		//deformation scaling
//		return fact * Vector3D(XGD(6), XGD(1), XGD(3)); 
//	}
//	virtual Vector2D GetDisplacement2DD() const
//	{
//		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
//		{
//			return Vector2D();
//		}
//
//		double fact = mbs->GetDOption(105);		//deformation scaling
//		return fact * Vector2D(XGD(6), XGD(1)); 
//	} 
//
//	virtual Vector3D GetVel() const	
//	{
//		return Vector3D(XGP(6), XGP(1), XGP(3)); 
//	}
//	virtual Vector2D GetVel2D() const
//	{
//		return Vector2D(XGP(6), XGP(1)); 
//	}
//	virtual Vector3D GetVelD() const 
//	{
//		return Vector3D(XGPD(6), XGPD(1), XGPD(3)); 
//	}
//	virtual Vector2D GetVel2DD() const
//	{ 
//		return Vector2D(XGPD(6), XGPD(1)); 
//	}
//
//
//protected:
//	Vector3D ref_angles;     //!EDC$[varaccess,EDCvarname="reference_rot_angles",EDCfolder="Geometry",tooltiptext="Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."]
//};
//!EDC$[endclass,NodeBeam3D]


//$EDC$[beginclass,classname=ANCFNodeS1rot1_3D,parentclassname=Node3D,addelementtypename=Node3DS1rot1,
//texdescription="Node3DS1rot1 is a finite element node for elements in 3D, and provides $7$ degrees of freedom.",
//texdescriptionGeometry="The reference geometry of the node is defined by the user via (a) $\mathtt{Geometry.reference\_position}$ and (b) the rotation $\mathtt{Geometry.reference\_rot\_angles}$. The rotation is prescribed by the user in form of kardan angles (initially, local $(S_1,S_2,S_3)$ and global frame $(x,y,z)$ are identical, then rotate local frame around $S1$, then $S2$ and finally $S3$). The current position is evaluated by adding displacement (the first three degrees of freedom) to the reference position of the node (degrees of freedom: $\left(q_1,q_2,q_3\right)^T$), and the current rotation of the node is obtained by adding the change of the first axis of the local frame (DOFs: $\left(q_4,q_5,q_6\right)^T$) to the first axis of the local frame in reference configuration of the node, and finally rotating the two other axes around the first axis of the local frame by the amount of the $7$th degree of freedom $q_7 = \theta$.",
//texdescriptionDOF="This node provides $7$ degrees of freedom: the first $3$ degrees of freedom are the displacement $(q_1,q_2,q_3)^T = \mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the next $3$ DOFs denote the change of the \emph{first slope}, which is the partial derivative of the position $\left(q_4,q_5,q_6\right)^T = \mathbf{r}_{,\xi}-\mathbf{r}_{0,\xi}$ with $\xi$ denoting the first of the three coordinates $(\xi,\eta,\zeta)$ of the reference element, and the $7$th degree of freedom is the local rotation $q_7 = \theta$ of the node around its current direction $S1$.",
//example="Node3DS1rot1.txt"]
class ANCFNodeS1rot1_3D : public Node3D
{
public:
	ANCFNodeS1rot1_3D():Node3D() {}

	ANCFNodeS1rot1_3D(MBS* mbs) : Node3D(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_3 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_angle + NDT_comp_1 + NDT_sos));
		InitConstructor();
	}

	void InitConstructor() 
	{
		Node3D::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_angles = Vector3D(0.);
	}

	ANCFNodeS1rot1_3D(const ANCFNodeS1rot1_3D& n) : Node3D(n.mbs) { CopyFrom(*this); }
	virtual ANCFNodeS1rot1_3D* GetCopy() 
	{
		ANCFNodeS1rot1_3D* n = new ANCFNodeS1rot1_3D();
		n->CopyFrom(*this);
		return n;
	}
	~ANCFNodeS1rot1_3D() {};

	virtual void CopyFrom(const ANCFNodeS1rot1_3D& n)
	{
		Node3D::CopyFrom(n);
		ref_angles = n.ref_angles;
	}


	void SetANCFNodeS1rot1_3D(const Vector3D& ref_pos, const Vector3D& ref_kardan_angles, int bodyind = 0) //ltg(2*ndof), x_init(), elements()
	{
		//sos = ndof; 
		body = bodyind; 
		pos = ref_pos;
		ref_angles = ref_kardan_angles;
	};

	
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual const Vector3D& GetRefAngles() const { return ref_angles; }

	virtual Matrix3D GetLocalFrame() const
	{ 
		return ComputeRotMatrixWithKardanAngles(ref_angles.X(), ref_angles.Y(), ref_angles.Z());
	}

	virtual double GetLocalFrame(int i, int j) const
	{
		//$ PG 2013-1-21: quick solution for now... has to be done more efficiently
		return GetLocalFrame()(i,j); 
	}

	virtual mystr GetTypeName() {return "Node3DS1rot1";}	
	virtual mystr GetObjectSpec() {return "Node3DS1rot1";}

protected:
	Vector3D ref_angles;     //$EDC$[varaccess,EDCvarname="reference_rot_angles",EDCfolder="Geometry",tooltiptext="Kardan rotation angles (X,Y,Z) in rad in global frame of node in reference configuration."]
};
//$EDC$[endclass,ANCFNodeS1rot1_3D]






//$EDC$[beginclass,classname=ANCFNodeS2S3_3D,parentclassname=Node3D,addelementtypename=Node3DS2S3,
//texdescription="Node3DS2S3 is a finite element node for elements in 3D, and provides $9$ degrees of freedom.",
//texdescriptionGeometry="The reference geometry of the node is defined by the user via (a) $\mathtt{Geometry.reference\_position}$ and (b) the slopes $\mathtt{Geometry.ref\_slope2}$ and $\mathtt{Geometry.ref\_slope3}$. The current position is evaluated by adding the displacement (the first three degrees of freedom $\left(q_1,q_2,q_3\right)^T$) to the reference position of the node, and further the current slopes of the node are obtained by adding the change of the second and third slopes (DOFs: $\left(q_4,q_5,q_6\right)^T$ and $\left(q_7,q_8,q_9\right)^T$ ) to the second and third slopes in reference configuration of the node.",
//texdescriptionDOF="This node provides $9$ degrees of freedom: the first $3$ degrees of freedom are the displacement $(q_1,q_2,q_3)^T = \mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the next $3$ DOFs denote the change of the \emph{second slope}, which are the partial derivatives of the position $\left(q_4,q_5,q_6\right)^T = \mathbf{r}_{,\eta}-\mathbf{r}_{0,\eta}$ and $\left(q_7,q_8,q_9\right)^T = \mathbf{r}_{,\zeta}-\mathbf{r}_{0,zeta}$, where $\eta$ and $\zeta$ denote the second and third of the three coordinates $(\xi,\eta,\zeta)$ of the reference element."]
class ANCFNodeS2S3_3D : public Node3D
{
public:
	ANCFNodeS2S3_3D():Node3D() {}

	ANCFNodeS2S3_3D(MBS* mbs) : Node3D(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_3 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope3 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope3 + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope3 + NDT_comp_3 + NDT_sos));

		InitConstructor();
	}

	void InitConstructor() 
	{
		Node3D::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_slope2 = Vector3D(0,1,0);
		ref_slope3 = Vector3D(0,0,1);	
	}

	ANCFNodeS2S3_3D(const ANCFNodeS2S3_3D& n) : Node3D(n.mbs) { CopyFrom(*this); }
	virtual ANCFNodeS2S3_3D* GetCopy() 
	{
		ANCFNodeS2S3_3D* n = new ANCFNodeS2S3_3D();
		n->CopyFrom(*this);
		return n;
	}
	~ANCFNodeS2S3_3D() {};

	virtual void CopyFrom(const ANCFNodeS2S3_3D& n)
	{
		Node3D::CopyFrom(n);
		ref_slope2 = n.ref_slope2;
		ref_slope3 = n.ref_slope3;
	}


	void SetANCFNodeS2S3_3D(const Vector3D& ref_pos, const Vector3D& ref_slope2, const Vector3D& ref_slope3, int bodyind = 0) //ltg(2*ndof), x_init(), elements()
	{
		//sos = ndof; 
		body = bodyind; 
		pos = ref_pos;
		this->ref_slope2 = ref_slope2;
		this->ref_slope3 = ref_slope3;
	};

	
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual const Vector3D& GetRefSlope2() const { return ref_slope2; }
	virtual const Vector3D& GetRefSlope3() const { return ref_slope3; }

	//virtual Matrix3D GetLocalFrame() const
	//{ 
	//	return ComputeRotMatrixWithKardanAngles(ref_angles.X(), ref_angles.Y(), ref_angles.Z());
	//}

	//virtual double GetLocalFrame(int i, int j) const
	//{
	//	//$ PG 2013-1-21: quick solution for now... has to be done more efficiently
	//	return GetLocalFrame()(i,j); 
	//}

	virtual mystr GetTypeName() {return "Node3DS2S3";}	
	virtual mystr GetObjectSpec() {return "Node3DS2S3";}

protected:
	Vector3D ref_slope2;     //$EDC$[varaccess,EDCvarname="reference_slope2",EDCfolder="Geometry",tooltiptext="slope 2 of node in reference configuration."]
	Vector3D ref_slope3;     //$EDC$[varaccess,EDCvarname="reference_slope3",EDCfolder="Geometry",tooltiptext="slope 3 of node in reference configuration."]
};
//$EDC$[endclass,ANCFNodeS2S3_3D]


//$EDC$[beginclass,classname=ANCFNodeS1S2_3D,parentclassname=Node3D,addelementtypename=Node3DS1S2,
//texdescription="Node3DS1S2 is a finite element node for elements in 3D, and provides $9$ degrees of freedom.",
//texdescriptionGeometry="The reference geometry of the node is defined by the user via (a) $\mathtt{Geometry.reference\_position}$ and (b) the slopes $\mathtt{Geometry.ref\_slope1}$ and $\mathtt{Geometry.ref\_slope2}$. The current position is evaluated by adding the displacement (the first three degrees of freedom $\left(q_1,q_2,q_3\right)^T$) to the reference position of the node, and further the current slopes of the node are obtained by adding the change of the first and second slopes (DOFs: $\left(q_4,q_5,q_6\right)^T$ and $\left(q_7,q_8,q_9\right)^T$ ) to the first and second slopes in reference configuration of the node.",
//texdescriptionDOF="This node provides $9$ degrees of freedom: the first $3$ degrees of freedom are the displacement $(q_1,q_2,q_3)^T = \mathbf{u} = \mathbf{r}-\mathbf{r}_0$, the next $3$ DOFs denote the change of the \emph{first slope}, which are the partial derivatives of the position $\left(q_4,q_5,q_6\right)^T = \mathbf{r}_{,\xi}-\mathbf{r}_{0,\xi}$ and $\left(q_7,q_8,q_9\right)^T = \mathbf{r}_{,\eta}-\mathbf{r}_{0,eta}$, where $\xi$ and $\eta$ denote the first and second of the three coordinates $(\xi,\eta,\zeta)$ of the reference element."]
class ANCFNodeS1S2_3D : public Node3D
{
public:
	ANCFNodeS1S2_3D():Node3D() {}

	ANCFNodeS1S2_3D(MBS* mbs) : Node3D(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope1 + NDT_comp_3 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_3 + NDT_sos));

		InitConstructor();
	}

	void InitConstructor() 
	{
		Node3D::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_slope1 = Vector3D(1,0,0);
		ref_slope2 = Vector3D(0,1,0);	
	}

	ANCFNodeS1S2_3D(const ANCFNodeS1S2_3D& n) : Node3D(n.mbs) { CopyFrom(*this); }

	virtual ANCFNodeS1S2_3D* GetCopy() 
	{
		ANCFNodeS1S2_3D* n = new ANCFNodeS1S2_3D();
		n->CopyFrom(*this);
		return n;
	}
	~ANCFNodeS1S2_3D() {};

	virtual void CopyFrom(const ANCFNodeS1S2_3D& n)
	{
		Node3D::CopyFrom(n);
		ref_slope1 = n.ref_slope1;
		ref_slope2 = n.ref_slope2;
	}


	void SetANCFNodeS1S2_3D(const Vector3D& ref_pos, const Vector3D& ref_slope1, const Vector3D& ref_slope2, int bodyind = 0) //ltg(2*ndof), x_init(), elements()
	{
		//sos = ndof; 
		body = bodyind; 
		pos = ref_pos;
		this->ref_slope1 = ref_slope1;
		this->ref_slope2 = ref_slope2;
	};

	
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual const Vector3D& GetRefSlope1() const { return ref_slope1; }
	virtual const Vector3D& GetRefSlope2() const { return ref_slope2; }

	virtual mystr GetTypeName() {return "Node3DS1S2";}	
	virtual mystr GetObjectSpec() {return "Node3DS1S2";}

protected:
	Vector3D ref_slope1;     //$EDC$[varaccess,EDCvarname="reference_slope1",EDCfolder="Geometry",tooltiptext="slope 1 of node in reference configuration."]
	Vector3D ref_slope2;     //$EDC$[varaccess,EDCvarname="reference_slope2",EDCfolder="Geometry",tooltiptext="slope 2 of node in reference configuration."]
};
//$EDC$[endclass,ANCFNodeS1S2_3D]



//$EDC$[beginclass,classname=ANCFNodeS2_2D,parentclassname=GenericNode,addelementtypename=Node2DS2]
class ANCFNodeS2_2D : public GenericNode
{
public:
	ANCFNodeS2_2D():GenericNode() {}

	ANCFNodeS2_2D(MBS* mbs) : GenericNode(mbs)
	{
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_displacement + NDT_comp_2 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_1 + NDT_sos));
		AddNodeDOF(NodeDOFType(NDT_slope2 + NDT_comp_2 + NDT_sos));
		InitConstructor();
	}

	void InitConstructor() 
	{
		GenericNode::InitConstructor();
		nodename = mystr(GetObjectSpec()); 
		ref_slope2 = Vector2D(0.);
	}

	ANCFNodeS2_2D(const ANCFNodeS2_2D& n) : GenericNode(n.mbs) { CopyFrom(*this); }
	virtual ANCFNodeS2_2D* GetCopy() 
	{
		ANCFNodeS2_2D* n = new ANCFNodeS2_2D();
		n->CopyFrom(*this);
		return n;
	}
	~ANCFNodeS2_2D() {};

	virtual void CopyFrom(const ANCFNodeS2_2D& n)
	{
		GenericNode::CopyFrom(n);
		ref_slope2 = n.ref_slope2;
	}

	virtual int Dim() const { return 3; }

	void SetANCFNodeS2_2D(const Vector3D& ref_pos, const Vector2D& ref_slope_eta, int bodyind = 0) //ltg(2*ndof), x_init(), elements()
	{
		body = bodyind; 
		pos = ref_pos;
		ref_slope2 = ref_slope_eta;
	};
	
	virtual void GetElementDataAuto(ElementDataContainer& edc); 		//fill in all element data
	virtual int SetElementDataAuto(ElementDataContainer& edc); //set element data acc. to ElementDataContainer, return 0 if failed or values invalid!
	virtual int ReadSingleElementDataAuto(ReadWriteElementDataVariableType& RWdata); 		//automatically generated function from EDC_converter
	virtual int WriteSingleElementDataAuto(const ReadWriteElementDataVariableType& RWdata); //automatically generated function from EDC_converter
	virtual int GetAvailableSpecialValuesAuto(TArrayDynamic<ReadWriteElementDataVariableType>& available_variables); //automatically generated function from EDC_converter

	virtual void GetElementData(ElementDataContainer& edc); 				//fill in all element data
	virtual int SetElementData(ElementDataContainer& edc);	//set element data acc. to ElementDataContainer, return 0 if failed or values invalid!

	virtual const Vector2D& GetRefSlope2() const { return ref_slope2; }


	virtual mystr GetTypeName() {return "Node2DS2";}	
	virtual mystr GetObjectSpec() {return "Node2DS2";}

	virtual Vector3D GetDisplacement() const { return Vector3D(XG(1), XG(2), 0); }
	virtual Vector2D GetDisplacement2D() const { return Vector2D(XG(1), XG(2)); }
	virtual Vector3D GetDisplacementD() const 
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector3D(0.);
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector3D(XGD(1), XGD(2), 0); 
	}
	virtual Vector2D GetDisplacement2DD() const
	{
		if (ltg.Length() == 0) //$ DR+PG 2013-01-22 XGD() not available as long as ltg has not been assembled.... is this solution ok?
		{
			return Vector2D();
		}

		double fact = mbs->GetDOption(105);		//deformation scaling
		return fact * Vector2D(XGD(1), XGD(2)); 
	} 

	virtual Vector3D GetVel() const	
	{
		return Vector3D(XGP(1), XGP(2),0); 
	}
	virtual Vector2D GetVel2D() const
	{
		return Vector2D(XGP(1), XGP(2)); 
	}
	virtual Vector3D GetVelD() const 
	{
		return Vector3D(XGPD(1), XGPD(2),0); 
	}
	virtual Vector2D GetVel2DD() const
	{ 
		return Vector2D(XGPD(1), XGPD(2)); 
	}


protected:
	Vector2D ref_slope2;     //$EDC$[varaccess,EDCvarname="reference_slope2",EDCfolder="Geometry",tooltiptext="Slope of cross section (direction eta) in reference configuration."]
};
//$EDC$[endclass,ANCFNodeS2_2D]

#endif

