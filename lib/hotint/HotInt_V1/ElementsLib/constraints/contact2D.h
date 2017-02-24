//#**************************************************************
//#
//# filename:             contact2D.h
//#
//# author:               Gerstmayr Johannes, Peter Gruber
//#
//# generated:						17. October 2006
//# description:          contact formulations
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
 
#ifndef CONTACT2D__H
#define CONTACT2D__H

#include "constraint.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GeneralContact 2D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GeneralContact2D: public Constraint
{
public:
	//Pos2DConstraint():loccoords(), dpdq() {mbs = NULL; };
	GeneralContact2D(MBS* mbsi):Constraint(mbsi),
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), mathfunctions(), master_adhesion_coefficient()
	{
	};
	GeneralContact2D(const GeneralContact2D& gc):Constraint(gc.mbs),
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), mathfunctions(), master_adhesion_coefficient()
	{CopyFrom(gc);};

	GeneralContact2D(MBS* mbsi, int slaveNODEmodei, double bordersizei,
		const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), 
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), mathfunctions(), master_adhesion_coefficient()
	{	

		slaveNODEmode = slaveNODEmodei;
		contactmode = 3;
		nlstepcnt = 0;

		bordersize = bordersizei; //increased search radius
		searchtreeix = 20;
		searchtreeiy = 20;
		contactmaxdist = 0.1;

		isfriction = 0;
		frictioncoeff = 0; //slip
		frictioncoeff_stick = 0; //stick
		restitution_coeff = 0.95; //coeff of restitution; 0.95 ususally
		hertzian_contact_param = 1; //coeff for Hertzian contact; 1 works good

		islagrange = 1;

		GetCol() = coli;
		draw_dim = cdim;
	};

	virtual ~GeneralContact2D() 
	{
		for (int i=1; i <= masterelements.Length(); i++)
		{
			delete masterelements(i);
		}
		for (int i=1; i <= mathfunctions.Length(); i++)
		{
			delete mathfunctions(i);
		}
		for (int i = 1; i <= pointradius.Length(); i++)
		{
			if (pointradius(i) != 0)
			{
				delete multimastercontact(i);
			}
		}
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new GeneralContact2D(*this);
		//Element* ec = new GeneralContact2D(mbs);
		//CopyFrom(*this); //funktioniert ohne nicht ... ?
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const GeneralContact2D& ce = (const GeneralContact2D&)e;

		elements_nodouble = ce.elements_nodouble;

		//slave:
		loccoords = ce.loccoords;
		locnodenums = ce.locnodenums;

		slaveNODEmode = ce.slaveNODEmode;
		slaveelements = ce.slaveelements;
		slavebodyindex = ce.slavebodyindex;
		pointradius = ce.pointradius;

		multimastercontact.SetLen(ce.multimastercontact.Length());
		for (int i=1; i <= ce.multimastercontact.Length(); i++)
		{
			if (ce.pointradius(i) != 0)
			{
				multimastercontact(i) = new IVector(*ce.multimastercontact(i));
			}
		}

		//master:
		masterlocind = ce.masterlocind;
		masterbodyindex = ce.masterbodyindex;

		masterelements.SetLen(ce.masterelements.Length());
		for (int i=1; i <= ce.masterelements.Length(); i++)
		{
			masterelements(i) = ce.masterelements(i)->GetCopy();
		}
		tempitemlist = ce.tempitemlist;

		masterisactive = ce.masterisactive;
		master_adhesion_coefficient = ce.master_adhesion_coefficient;
		mathfunctions.SetLen(ce.mathfunctions.Length());
		for (int i=1; i <= ce.mathfunctions.Length(); i++)
		{
			mathfunctions(i) = ce.mathfunctions(i)->GetCopy();
		}


		mastersearchtree = ce.mastersearchtree;
		slavesearchtree = ce.slavesearchtree;
		bordersize = ce.bordersize;
		searchtreeix = ce.searchtreeix;
		searchtreeiy = ce.searchtreeiy;

		c_contact = ce.c_contact;
		gapp_init = ce.gapp_init;
		gapp_init_min = ce.gapp_init_min;
		iscontact = ce.iscontact;

		contactmode = ce.contactmode;
		nlstepcnt = ce.nlstepcnt;
		contactmaxdist = ce.contactmaxdist;
		islagrange = ce.islagrange;

		restitution_coeff = ce.restitution_coeff;
		hertzian_contact_param = ce.hertzian_contact_param;

		isfriction = ce.isfriction;
		frictioncoeff = ce.frictioncoeff;
		frictioncoeff_stick = ce.frictioncoeff_stick;
		isstick = ce.isstick;
		islaststick = ce.islaststick;
		laststickp = ce.laststickp;
		slipdir = ce.slipdir;


		if (ce.NMC())
		{
			mastertoslave = new IVector[ce.NMC()]();
			for (int i = 0; i < ce.NMC(); i++)
			{
				mastertoslave[i].Init();
				mastertoslave[i] = ce.mastertoslave[i];
			}
		}

		dpdq = ce.dpdq;


	}

	virtual void LinkToElements() 
	{
		for (int i = 1; i <= elements.Length(); i++)
		{
			if (elements(i) != 0)	GetElem(i).AddConstraint(this,i);
		}
	}

	//elements not doubled: only for certain elements!!!
	virtual int NE_nodouble() const {return elements_nodouble.Length();}
	virtual const Element& GetElem_nodouble(int i) const {return mbs->GetElement(elements_nodouble(i));}
	virtual Element& GetElem_nodouble(int i) {return mbs->GetElement(elements_nodouble(i));}

	virtual const char* GetElementSpec() const {return "Contact2D";}
	virtual int CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	{ return Element::CheckConsistency(errorstr);}

	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void EvalG(Vector& f, double t);

	//To be replaced in derived class
	virtual void AddElementCqTLambda(double t, int locelemind, Vector& f);

	//implicit (algebraic) size: normal contacts+tangential contacts
	virtual int IS() const 
	{
		return slaveelements.Length()*islagrange;
	};

	virtual int Dim() const {return GetElem(1).Dim();}  //has 3D position???

	//get a reference position of the body in 3d
	virtual Vector3D GetRefPosD() const;

	virtual void DrawElement();

	//+++++++++++++++++++++++++++++++++++++++++++++++++++
	//general contact specific functions:

	virtual void AddElementCoord(int en, const Vector2D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void AddSlaveNode(int elnum, const Vector2D& loccoord, double cstiff, int bodyindex = 1, double pointradiusi = 0)
	{
		slaveelements.Add(elnum);
		loccoords.Add(loccoord);
		pointradius.Add(pointradiusi);
		c_contact.Add(cstiff);
		slavebodyindex.Add(bodyindex);
	}

	virtual void AddSlaveNode(int elnum, int locnodenum, double cstiff, int bodyindex = 1, double pointradiusi = 0)
	{
		slaveelements.Add(elnum);
		locnodenums.Add(locnodenum);
		pointradius.Add(pointradiusi);
		c_contact.Add(cstiff);
		slavebodyindex.Add(bodyindex);
	}

	virtual int AddMathFunction(MathFunction* mf)
	{
		mathfunctions.Add(mf->GetCopy());
		return mathfunctions.Length();
	}

	virtual void AddMasterGeomElement(GeomElement* ge, int bodyindex, int isactivefunction = 0, double adhesion_coeff=0)
	{
		GeomElement* gen = ge->GetCopy();
		masterelements.Add(gen);
		masterlocind.Add(1);
		masterisactive.Add(isactivefunction);
		masterbodyindex.Add(bodyindex);
		master_adhesion_coefficient.Add(adhesion_coeff);
	}
	
	virtual void AddMasterSegment(int elnum, const Vector2D& locsegpos1i, const Vector2D& locsegpos2i, int bodyindex = 0)
	{
		GeomLine2D gl(GetMBS(), elnum, locsegpos1i, locsegpos2i, Vector3D(0.8, 0.2, 0.2));
		AddMasterGeomElement(&gl, bodyindex);
	}

	virtual void AddMasterSegment(int elnum, int locsegnum1i, int locsegnum2i, int bodyindex = 0)
	{
		GeomLine2D gl(GetMBS(), elnum, locsegnum1i, locsegnum2i, Vector3D(0.8, 0.2, 0.2));
		AddMasterGeomElement(&gl, bodyindex);
	}

	// AP 2013-01-11: removed source code to cpp file, otherwise FinishContactDefinition cannot be called in model file
	virtual void FinishContactDefinition();

	virtual void Initialize() 
	{

		BuildSearchTrees();

	}
	
	//number of slave contact points
	virtual int NSC() const 
	{
		return slaveelements.Length();
	}
	//number of master contact segments
	virtual int NMC() const 
	{
		return masterelements.Length();
	}

	virtual void BuildSearchTrees();

	virtual void SetContactMode(int contactmodeI) {contactmode = contactmodeI;}
	virtual void SetIsLagrange(int islagrangeI) {islagrange = islagrangeI;}
	virtual void SetFriction(int isfrictionI, double frictioncoeffI, double frictioncoeff_stickI = -1) 
	{
		isfriction = isfrictionI; 
		frictioncoeff = frictioncoeffI;
		if (frictioncoeff_stickI == -1) frictioncoeff_stick = frictioncoeff;
		else frictioncoeff_stick = frictioncoeff_stickI;
	}
	virtual void SetContactParams(double restitution_coeffi, double hertzian_contact_parami)
	{
		restitution_coeff = restitution_coeffi;
		hertzian_contact_param = hertzian_contact_parami;
	}
	virtual void SetContactMaxDist(double d) {contactmaxdist = d;}
	virtual void SetSearchTreeDim(int ix, int iy) {searchtreeix = ix; searchtreeiy = iy;}

	virtual const Body2D& GetSlaveBody(int i) const;
	virtual Body2D& GetSlaveBody(int i) ;
	virtual const Body2D& GetMasterBody(int i) const;
	virtual Body2D& GetMasterBody(int i) ;

	virtual Vector2D GetSlaveNodePos(int i) const;
	virtual Vector2D GetSlaveNodePosD(int i) const;
	virtual Vector2D GetSlaveNodeVel(int i) const;

	//approximate velocity at segment point pglob due linear interpolation
	virtual Vector2D GetMasterSegVel(const Vector2D& pglob, int j, int locind) const;
	virtual Vector2D GetMasterLocPos(const Vector2D& pp, int j, int locind) const;
	
	//get nearest possible master segment
	virtual void GetNearestMasterSegment(int slavenum, const Vector2D& p, int& masternum, int& locind, double& dist, Vector2D& pp);
	
	//get all master segments with negative gap, including circleradius of slave node
	virtual void GetNearestMasterSegments(int slavenum, const Vector2D& p, double circlerad, TArray<int>& a_masternum, 
		TArray<int>& a_locind, TArray<double>& a_dist, TArray<Vector2D>& a_pp);

	virtual double GetContactForce(int i, double gap, double gapp); //penalty contact force for slave node i
	virtual void GetTangentStickForce(int i, const Vector2D& t, double tvel, const Vector2D& pp, const Vector2D& lastpp, Vector2D& tforce); //get tangent force for stick

	virtual GeomElement* GeomElem(int i) const {return masterelements(i);}
	virtual int GetMasterElnum(int i) const 
	{
		return GeomElem(i)->GetElnum();
	}

	virtual MathFunction* GetMasterIsActiveMathFunction(int i) const
	{
		if (i > mathfunctions.Length()) 
		{
			GetMBS()->UO() << "warning: invalid Mathfunctions!\n";
			return mathfunctions.Last();
		}
		return mathfunctions(i);
	}
	virtual int GetMasterIsActiveCondition(double t, int masterelem) const;

	// maximum number of post-newton steps, afterwards error=0 is returned
	// to avoid non-converging iterations
	virtual int MaxNLIt() { return 3; }

protected:
	TArray<int> elements_nodouble; //every element is added only once, for automatic jacobian!!!

	//Slave Elements:
	TArray<int> slaveelements;
	TArray<Vector2D> loccoords; //Contact node positions
	TArray<int> locnodenums;    //Contact node numbers
	int slaveNODEmode;					//if NODEmode = 1, then use locnodenumbers, if NODEmode==0 then use loccoords

	TArray<double> c_contact;				//contact stiffness for every slave node
	TArray<double> gapp_init;				//initial contact velocity, slave node
	TArray<double> gapp_init_min;		//minimum initial contact vel, slave node
	TArray<int> iscontact;					//for each slave node: if no penetration in last nonlinit -> iscontact = 0, else iscontact = mastersegment number
	TArray<int> slavebodyindex;			//body number for slave element
	TArray<double> pointradius;			//for slave node, radius around each point, for circle contact
	TArray<IVector*> multimastercontact; //master contact elements, if pointradius != 0, contains master element numbers if contact

	//temporary for GetNearestMasterSegments for multiple contact!
	TArray<int> tempmasternum;			//will not be copied!!!
	TArray<double> tempdist;
	TArray<int> templocind;
	TArray<Vector2D> temp_pp;
	IVector tempitemlist;						//masteritemlist, for sorting

	//new Master Elements:
	TArray<GeomElement*> masterelements;  //geometric element (line, circle, polygon, etc.)
	TArray<int> masterlocind;							  //identifies master local node / segment etc.
	IVector* mastertoslave;								//mastertoslave[i-1] list contains all slave elements which are in contact with master element i
	TArray<int> masterbodyindex;					//body number for each master segment/element -->redundant with GeomElement???
	IVector masterisactive;								//if 0 --> master is always active, otherwise use Condition of function i
	TArray<MathFunction*> mathfunctions;	//for masterisactive
	TArray<double> master_adhesion_coefficient; //use this in order to add glue to master surface

	//SearchTree
	SearchTree2D mastersearchtree; //for optimized node fill
	SearchTree2D slavesearchtree;  //for optimized node fill, not so far needed in GeneralContact2D, but in SPHContact2D!
	TArray<int> itemlist;      //item list for GetItem command

	double bordersize; //additional search radius for master and slave segments/nodes
	int searchtreeix;  //number of hash entries for searchtree in x-direction
	int searchtreeiy;  //number of hash entries for searchtree in y-direction

	int contactmode; //0= contact, 1=no loss of contact (always contact stiffness), 2=no loss of contact, 5=improved rest.coeff, 4=lagrange mult., no stiffness
	int nlstepcnt;
	double contactmaxdist; //maximum distance of penetration, after that it is treated as no contact
	int islagrange;		//islagrange = 1: contact formulation with lagrange multipliers 
	double restitution_coeff; //coeff of restitution; 0.95 ususally
	double hertzian_contact_param; //coeff for Hertzian contact; 1 works good

	//FRICTION:
	int	isfriction;		//friction on/off
	double frictioncoeff;	//friction parameter, slip
	double frictioncoeff_stick; //friction parameter stick
	TArray<int> slipdir;         //for each slave node: slipdir = 1: in direction master segment node 1 to node 2, else: slipdir = -1
	TArray<int> isstick;         //for each slave node: stick = 1 / slip = 0
	TArray<int> islaststick;     //for each slave node: last successfull step was stick = 1 / slip = 0
	TArray<Vector2D> laststickp; //for each slave node: last (global) stick position

	Matrix dpdq; //temporary element, no need to copy???
};


#endif
