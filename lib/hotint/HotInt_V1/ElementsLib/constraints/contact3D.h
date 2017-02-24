//#**************************************************************
//#
//# filename:             contact3D.h
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
 
#ifndef CONTACT3D__H
#define CONTACT3D__H

#include "constraint.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GeneralContact 3D
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class GeneralContact3D: public Constraint
{
public:
	GeneralContact3D(MBS* mbsi):Constraint(mbsi),
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		multimastercontact_last(), multimastercontact_actlast(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), 
		mathfunctions(), constantmastersearchtree()
	{
	};
	GeneralContact3D(const GeneralContact3D& gc):Constraint(gc.mbs),
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		multimastercontact_last(), multimastercontact_actlast(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), 
		mathfunctions(), constantmastersearchtree()
	{CopyFrom(gc);};

	GeneralContact3D(MBS* mbsi, int slaveNODEmodei, double bordersizei,
		const Vector3D& cdim, const Vector3D& coli):Constraint(mbsi), 
		loccoords(), dpdq(), locnodenums(), masterelements(), masterlocind(), slaveelements(),
		mastersearchtree(), slavesearchtree(), c_contact(), gapp_init(), gapp_init_min(), iscontact(), itemlist(),
		slavebodyindex(), masterbodyindex(), isstick(), islaststick(), laststickp(), slipdir(), pointradius(), multimastercontact(),
		multimastercontact_last(), multimastercontact_actlast(),
		tempmasternum(), tempdist(), templocind(), temp_pp(), tempitemlist(), masterisactive(), 
		mathfunctions(), constantmastersearchtree()
	{	

		slaveNODEmode = slaveNODEmodei;
		contactmode = 3;
		nlstepcnt = 0;

		bordersize = bordersizei; //increased search radius
		searchtreeix = 20;
		searchtreeiy = 20;
		searchtreeiz = 20;
		mastersearchtreeinitialized = 0;

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

	virtual ~GeneralContact3D() 
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
				delete multimastercontact_last(i);
				delete multimastercontact_actlast(i);
			}
		}
	}

	//To be overwritten in derived class:
	virtual Element* GetCopy()
	{
		Element* ec = new GeneralContact3D(*this);
		//Element* ec = new GeneralContact3D(mbs);
		//CopyFrom(*this); //funktioniert ohne nicht ... ?
		return ec;
	}
	//To be overwritten in derived class:
	virtual void CopyFrom(const Element& e)
	{
		Constraint::CopyFrom(e);
		const GeneralContact3D& ce = (const GeneralContact3D&)e;

		elements_nodouble = ce.elements_nodouble;

		//slave:
		loccoords = ce.loccoords;
		locnodenums = ce.locnodenums;

		slaveNODEmode = ce.slaveNODEmode;
		slaveelements = ce.slaveelements;
		slavebodyindex = ce.slavebodyindex;
		pointradius = ce.pointradius;

		multimastercontact.SetLen(ce.multimastercontact.Length());
		multimastercontact_last.SetLen(ce.multimastercontact_last.Length());
		multimastercontact_actlast.SetLen(ce.multimastercontact_actlast.Length());
		for (int i=1; i <= ce.multimastercontact.Length(); i++)
		{
			if (ce.pointradius(i) != 0)
				multimastercontact(i) = new IVector(*ce.multimastercontact(i));
		}
		for (int i=1; i <= ce.multimastercontact_last.Length(); i++)
		{
			if (ce.pointradius(i) != 0)
				multimastercontact_last(i) = new IVector(*ce.multimastercontact_last(i));
		}
		for (int i=1; i <= ce.multimastercontact_actlast.Length(); i++)
		{
			if (ce.pointradius(i) != 0)
				multimastercontact_actlast(i) = new IVector(*ce.multimastercontact_actlast(i));
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
		mathfunctions.SetLen(ce.mathfunctions.Length());
		for (int i=1; i <= ce.mathfunctions.Length(); i++)
		{
			mathfunctions(i) = ce.mathfunctions(i)->GetCopy();
		}


		mastersearchtree = ce.mastersearchtree;
		constantmastersearchtree = ce.constantmastersearchtree;
		slavesearchtree = ce.slavesearchtree;
		bordersize = ce.bordersize;
		searchtreeix = ce.searchtreeix;
		searchtreeiy = ce.searchtreeiy;
		searchtreeiz = ce.searchtreeiz;
		mastersearchtreebox = ce.mastersearchtreebox;
		mastersearchtreeinitialized = ce.mastersearchtreeinitialized;			

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

			mastertoslave_last = new IVector[ce.NMC()]();
			for (int i = 0; i < ce.NMC(); i++)
			{
				mastertoslave_last[i].Init();
				mastertoslave_last[i] = ce.mastertoslave_last[i];
			}

			mastertoslave_actlast = new IVector[ce.NMC()]();
			for (int i = 0; i < ce.NMC(); i++)
			{
				mastertoslave_actlast[i].Init();
				mastertoslave_actlast[i] = ce.mastertoslave_actlast[i];
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

	virtual const char* GetElementSpec() const {return "Contact";}
	virtual int CheckConsistency(mystr& errorstr) //rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
	{ return Element::CheckConsistency(errorstr);}

	virtual Box3D GetElementBoxD() const
	{
		Box3D b  = Constraint::GetElementBoxD();

		for (int i = 1; i <= NMC(); i++)
		{
			b.Add(GeomElem(i)->GetBoundingBoxD());
		}
		return b;
	}

	//elements not doubled: only for certain elements!!
	virtual int NE_nodouble() const {return elements_nodouble.Length();}
	virtual const Element& GetElem_nodouble(int i) const {return mbs->GetElement(elements_nodouble(i));}
	virtual Element& GetElem_nodouble(int i) {return mbs->GetElement(elements_nodouble(i));}


	virtual double PostNewtonStep(double t);
	virtual void PostprocessingStep();

	virtual void EvalG(Vector& f, double t);

	//compute gap, gapp, normal, tangents, etc. some computations that are done in AddElementCqTLambda
	virtual void ComputeGap(int i, int mnum, double& gap, double& gapp, Vector3D& p, Vector3D& pp, 
		Vector3D& ploc, Vector3D& t1, Vector3D& t2, Vector3D& n, double& tvel, double& forcefact, double& forceadd, int nlstep);

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

	virtual void AddElementCoord(int en, const Vector3D& loccoord)
	{
		AddElement(en);
		loccoords.Add(loccoord);
	}

	virtual void AddSlaveNode(int elnum, const Vector3D& loccoord, double cstiff, int bodyindex = 1, double pointradiusi = 0)
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

	virtual void AddMasterGeomElement(GeomElement* ge, int bodyindex, int isactivefunction = 0)
	{
		GeomElement* gen = ge->GetCopy();
		masterelements.Add(gen);
		masterlocind.Add(1);
		masterisactive.Add(isactivefunction);
		masterbodyindex.Add(bodyindex);
	}

	virtual void FinishContactDefinition()
	{
		x_init = Vector(IS()); //number of normal contact lagrange multipliers == number of slave nodes

		int nnc = NSC(); //number of normal contacts
		iscontact.SetLen(nnc);
		gapp_init.SetLen(nnc);
		gapp_init_min.SetLen(nnc);

		isstick.SetLen(nnc);
		islaststick.SetLen(nnc);
		laststickp.SetLen(nnc);
		slipdir.SetLen(nnc);
		multimastercontact.SetLen(nnc); //only elements set where pointradius != 0
		multimastercontact_last.SetLen(nnc); //only elements set where pointradius != 0
		multimastercontact_actlast.SetLen(nnc); //only elements set where pointradius != 0

		for (int i=1; i <= NSC(); i++)
		{
			iscontact(i) = 0;
			isstick(i) = 0;
			islaststick(i) = 0;
			laststickp(i) = Vector3D(0.,0.,0.);
			slipdir(i) = 1;

			gapp_init_min(i) = 1e-6; //1e-2 for hourglass
			gapp_init(i) = gapp_init_min(i);

			if (pointradius(i) != 0)
			{
				multimastercontact(i) = new IVector(2);
				multimastercontact_last(i) = new IVector(2);
				multimastercontact_actlast(i) = new IVector(2);
			}
		}
		elements.SetLen(0);
		for (int i = 1; i <= NSC(); i++)
		{	
			elements.Add(slaveelements(i));
		}
		for (int i = 1; i <= NMC(); i++)
		{	
			elements.Add(masterelements(i)->GetElnum());
		}

		//for element jacobian:
		elements_nodouble.SetLen(0);
		for (int i = 1; i <= NSC(); i++)
		{	
			if (!Find(slaveelements(i), elements_nodouble))
			{
				elements_nodouble.Add(slaveelements(i));
			}
		}

		for (int i = 1; i <= NMC(); i++)
		{	
			if (!Find(masterelements(i)->GetElnum(), elements_nodouble))
			{
				elements_nodouble.Add(masterelements(i)->GetElnum());
			}
		}

		mastertoslave = new IVector[NMC()]();
		for (int i = 0; i < NMC(); i++)
		{
			mastertoslave[i].Init();
		}

		mastertoslave_last = new IVector[NMC()]();
		for (int i = 0; i < NMC(); i++)
		{
			mastertoslave_last[i].Init();
		}

		mastertoslave_actlast = new IVector[NMC()]();
		for (int i = 0; i < NMC(); i++)
		{
			mastertoslave_actlast[i].Init();
		}
		
		tempitemlist.SetLen(NMC());
		for (int i = 1; i <= NMC(); i++)
			tempitemlist(i) = 0;
	}

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
		if (contactmode == 5)
			restitution_coeff = ConvertRestitution(restitution_coeffi);
		else
			restitution_coeff = restitution_coeffi;

		hertzian_contact_param = hertzian_contact_parami;
	}
	virtual double ConvertRestitution(double restcoeff);

	virtual void SetContactMaxDist(double d) {contactmaxdist = d;}
	virtual void SetSearchTreeDim(int ix, int iy, int iz) {searchtreeix = ix; searchtreeiy = iy; searchtreeiz = iz;}

	virtual const Body3D& GetSlaveBody(int i) const;
	virtual Body3D& GetSlaveBody(int i) ;
	virtual const Body3D& GetMasterBody(int i) const;
	virtual Body3D& GetMasterBody(int i) ;

	virtual Vector3D GetSlaveNodePos(int i) const;
	virtual Vector3D GetSlaveNodePosD(int i) const;
	virtual Vector3D GetSlaveNodeVel(int i) const;

	//approximate velocity at segment point pglob due linear interpolation
	virtual Vector3D GetMasterSegVel(const Vector3D& pglob, int j, int locind) const;
	virtual Vector3D GetMasterLocPos(const Vector3D& pp, int j, int locind) const;
	
	//get nearest possible master segment
	virtual void GetNearestMasterSegment(int slavenum, const Vector3D& p, int& masternum, int& locind, 
		double& dist, Vector3D& pp);
	
	//get all master segments with negative gap, including circleradius of slave node
	virtual void GetNearestMasterSegments(int slavenum, const Vector3D& p, double circlerad, TArray<int>& a_masternum, 
		TArray<int>& a_locind, TArray<double>& a_dist, TArray<Vector3D>& a_pp);

	virtual double GetContactForce(int i, double gap, double gapp); //penalty contact force for slave node i
	virtual void GetTangentStickForce(int i, const Vector3D& t, double tvel, const Vector3D& pp, const Vector3D& lastpp, Vector3D& tforce); //get tangent force for stick

	virtual GeomElement* GeomElem(int i) const {return masterelements(i);}
	virtual int GetMasterElnum(int i) const 
	{
		return GeomElem(i)->GetElnum();
	}

protected:
	TArray<int> elements_nodouble; //every element is added only once, for automatic jacobian!!!

	//Slave Elements:
	TArray<int> slaveelements;
	TArray<Vector3D> loccoords; //Contact node positions
	TArray<int> locnodenums;    //Contact node numbers
	int slaveNODEmode;					//if NODEmode = 1, then use locnodenumbers, if NODEmode==0 then use loccoords

	TArray<double> c_contact;				//contact stiffness for every slave node
	TArray<double> gapp_init;				//initial contact velocity, slave node
	TArray<double> gapp_init_min;		//minimum initial contact vel, slave node
	TArray<int> iscontact;					//for each slave node: if no penetration in last nonlinit -> iscontact = 0, else iscontact = mastersegment number
	TArray<int> slavebodyindex;			//body number for slave element
	TArray<double> pointradius;			//for slave node, radius around each point, for circle contact
	TArray<IVector*> multimastercontact; //master contact elements, if pointradius != 0, contains master element numbers if contact
	TArray<IVector*> multimastercontact_last; //master contact elements of last time step
	TArray<IVector*> multimastercontact_actlast; //master contact elements which where active in last step or are active in current step

	//temporary for GetNearestMasterSegments for multiple contact!
	TArray<int> tempmasternum;			//will not be copied!!!
	TArray<double> tempdist;
	TArray<int> templocind;
	TArray<Vector3D> temp_pp;
	IVector tempitemlist;						//masteritemlist, for sorting

	//new Master Elements:
	TArray<GeomElement*> masterelements;  //geometric element (line, circle, polygon, etc.)
	TArray<int> masterlocind;							//identifies master local node / segment etc.
	IVector* mastertoslave;								//mastertoslave[i-1] list contains all slave elements which are in contact with master element i
	IVector* mastertoslave_last;				  //mastertoslave[i-1] list from end of last time step
	IVector* mastertoslave_actlast;				//mastertoslave[i-1] list from, compilation of contacts active at beginning or end of actual time step
	TArray<int> masterbodyindex;					//body number for each master segment/element -->redundant with GeomElement???
	IVector masterisactive;								//if 0 --> master is always active, otherwise use Condition of function i
	TArray<MathFunction*> mathfunctions;	//for masterisactive

	//SearchTree
	SearchTree mastersearchtree;					//for optimized node fill
	SearchTree constantmastersearchtree;	//for elements that do not move
	SearchTree slavesearchtree;						//for optimized node fill, NOT NEEDED!!!
	TArray<int> itemlist;									//item list for GetItem command
	Box3D mastersearchtreebox;						//initial box for mastersearchtree
	int mastersearchtreeinitialized;			//for first initialization

	double bordersize; //additional search radius for master and slave segments/nodes
	int searchtreeix;  //number of hash entries for searchtree in x-direction
	int searchtreeiy;  //number of hash entries for searchtree in y-direction
	int searchtreeiz;  //number of hash entries for searchtree in z-direction

	int contactmode; //0= contact, 1=no loss of contact (always contact stiffness), 2=no loss of contact, lagrange mult., no stiffness
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
	TArray<Vector3D> laststickp; //for each slave node: last (global) stick position

	Matrix dpdq; //temporary element, no need to copy???
};


#endif
