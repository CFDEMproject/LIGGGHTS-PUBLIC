//#**************************************************************
//#
//# filename:             mbs.cpp
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
//#***************************************************************************************
 

#include "mbs.h"
//#include "windows.h"					//for shell execute
#include <direct.h>						//for getcwd
#include  <io.h>							//file operation _access
#include  <stdio.h>
#include  <stdlib.h>

#include "parser.h"
#include "script_parser.h"
#include "element.h"
#include "sensors.h"
#include "node.h"
#include "material.h"
#include "parser.h"

#define MBS_use_parallelization_NOT

#ifdef MBS_use_parallelization
	#include <omp.h> //for parallelization
	const int OMP_MAX_N_THREADS=12;
#endif

#include "lapack_routines.h"
#include "options_mbs_auto.h"

//the following file contains the conversion MBSSolSet <==> EDC
#include "mbssolset_convert_auto.h"

extern MultiBodySystem * TIMBS;

void TIMBSWarningHandle(const char* warn, int use_instant_message_text)
{
//$ AD 2011-11: flag may have negative values, then no output
	if(use_instant_message_text < 0) return;

	if(!use_instant_message_text)
		TIMBS->UO(UO_LVL_warn) << "WARNING:" << warn << "\n";
	else
		TIMBS->UO(UO_LVL_warn).InstantMessageText(mystr("WARNING:")+mystr(warn));

}

//AP 2013-01-15: constructor TimeInt() is called explicitely (here *global_uo is set, thus it is essential)
MultiBodySystem::MultiBodySystem(): TimeInt(), elements(), nodes(), auxelements(), sensors(), indexOfActualPostProcessingFieldVariable(0), eigval()
{
	maxindex = 2;
	magnify_yz = 1;

	//++++++++++ parameter variation/optimization
	varparameter = 1;
	varparameter2 = 0;
	ptr2PerformComputation_Function = 0;
	//++++++++++

	resortsize_add = 0;
	//evalfcnt = 0;
	//evalf_jaccnt = 0;
	//writesolcnt = 0;
	//testcnt = 0;

	centerobject = 0;
	centerobject_offset = Vector3D(0.,0.,0.);

	transformJacApply = 0;
	prohibit_redraw = 0;

	elementbandwidth = 16;
	use_dependencies = 0;

	isfirstcomputation = 1;
	isinconsistent = 0;
	isloadsavemode = 0;

	resort_active = 1;

	edc_modeldata = 0; //empty!
	modeldata_initialized = 0;

	mbs_edc_options = 0;
	mbs_edc_global_variables = 0;

	ElementDataContainer edc_init;
	SetMBS_EDC_Options(edc_init);

	ElementDataContainer edc_variables;
	SetMBS_EDC_Variables(edc_variables);

	mbs_parser = new CMBSParser();
	MBSParser().SetMBS(this);
	MBSParser().Init();

	edc_parser = new CEDCParser();
	edc_parser->SetMBS(this, mbs_parser);
	edc_parser->Initialize();

	stored_initialconditions.SetLen(0);
	isSolutionFileHeaderWritten = 0;
	isSolParFileHeaderWritten = 0;

	//$ MS + RL 2012-2-29: coreComputationFinished=0; //$ MS 2011-3-16: 
	simulationStatus.SetStatusFlag(TSimulationNotStarted); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//this call adds all class .tex descriptions for online or offline access
	//EDCParser().AddClassDescriptions_Auto();

	serversocket = NULL;
}

MultiBodySystem::~MultiBodySystem()
{
	//delete everyhing newed!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Destroy();

	if (mbs_edc_options) delete mbs_edc_options;
	if (edc_modeldata) delete edc_modeldata;
	if (mbs_edc_global_variables) delete mbs_edc_global_variables;
	if (mbs_parser) delete mbs_parser;
}


void MultiBodySystem::Destroy()
{
	for (int i=1; i <= elements.Length(); i++)
	{
		delete elements(i);
	}
	for (int i=1; i <= auxelements.Length(); i++)
	{
		delete auxelements(i);
	}
	for (int i=1; i <= nodes.Length(); i++)
	{
		delete nodes(i);
	}
	for (int i=1; i <= sensors.Length(); i++)
	{
			delete sensors(i);
	}
	for (int i=1; i <= loads.Length(); i++)	//$ DR 21012-10
	{
		delete loads(i);
	}
	elements.Flush();
	nodes.Flush();
	auxelements.Flush();
	sensors.Flush();
	loads.Flush(); //$ DR 21012-10

	for (int i=1; i <= drawelements.Length(); i++)
	{
		delete drawelements(i);
	}
	drawelements.Flush();

	for (int i=1; i <= materials.Length(); i++)
	{
		delete materials(i);
	}
	materials.Flush();

	NumSolver().DestroyOldJacs();

	mbs_parser->GetVariableList()->Flush();
}

void MultiBodySystem::InitFirst() 
{
	TimeInt::InitFirst();

	//init global variables
	InitGlobalEDCVariables();
}

//add global variables
void MultiBodySystem::InitGlobalEDCVariables()
{
	ElementDataContainer* edc = GetMBS_EDC_Variables();
	ElementData ed;

	ed.SetDouble(MY_PI, "pi"); edc->Add(ed);
}

int MultiBodySystem::CheckConsistencyLTG()
{
	//check LTG (due to delete or internal error)

	MBSsos = GetSecondOrderSize();
	MBSes = GetFirstOrderSize();
	MBSis = GetImplicitSize();
	MBSss = GetSystemSize();

	int l2problem = 0; //level 2 problem: can not compute and not draw

	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		const	TArray<int>& ltg = e->GetLTGArray();
		if (ltg.Length() != e->SS()) 
		{
			l2problem = 1;
			UO(UO_LVL_err) << "System inconsistency: Element " << i << " has an inconsistency in the local to global DOF reference numbers!\n";
		}

		for (int j=1; j<=ltg.Length(); j++)
		{
			if (ltg.Get(j) <= 0 || ltg.Get(j) > MBSss)
			{
				l2problem = 1;
				UO(UO_LVL_err) << "System inconsistency: Element " << i << " has invalid local to global DOF reference numbers!\n";
			}
		}
	}

	if(l2problem) return 2;
	else return 0;
}

int MultiBodySystem::CheckConsistencyElements()
{
	int l1problem = 0; //level 1 problem: can not compute
	int l2problem = 0; //level 2 problem: can not compute and not draw
	mystr errorstr = "";

	for (int i=1; i<=elements.Length(); i++) 
	{
		int rv = elements(i)->CheckConsistency(errorstr);
		if (rv == 1) l1problem = 1;
		else if (rv == 2) l2problem = 1;

		if (rv)
		{
			UO(UO_LVL_err) << "System inconsistency in element " << i << ":\n";
			UO(UO_LVL_err) << errorstr;
			errorstr = "";
			if (rv == 1) l1problem = 1;
			else if (rv == 2) l2problem = 1;
		}
	}

	if (l2problem) return 2;
	else if (l1problem) return 1;

	return 0;
}

int MultiBodySystem::CheckConsistencySensors()
{
	int l1problem = 0; //level 1 problem: can not compute
	int l2problem = 0; //level 2 problem: can not compute and not draw
	mystr errorstr = "";

	for (int i=1; i <= NSensors(); i++)
	{
		//$ YV 2012-06: for sensors true is consistent and false is not
		int rv = sensors(i)->IsConsistent(errorstr) ? 0 : 1;
		if (rv == 1) l1problem = 1;
		//else if (rv == 2) l2problem = 1;

		if (rv)
		{
			UO(UO_LVL_err) << "System inconsistency in sensor " << i << ":\n";
			UO(UO_LVL_err) << errorstr;
			errorstr = "";
			if (rv == 1) l1problem = 1;
			else if (rv == 2) l2problem = 1;
		}
	}

	if (l2problem) return 2;
	else if (l1problem) return 1;

	return 0;
}

int MultiBodySystem::CheckConsistencyNodes()
{
	int l1problem = 0; //level 1 problem: can not compute
	int l2problem = 0; //level 2 problem: can not compute and not draw
	mystr errorstr = "";

	//check if all nodes are used:
	TArray<int> nodescheck;
	nodescheck.SetLen(NNodes());
	for (int i=1; i <= nodescheck.Length(); i++)
	{
		nodescheck(i) = 0;
	}
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		if (e->PerformNodeCheck())
		{
			for (int j=1; j <= e->NNodes(); j++)
			{
				//nodenumber is already verified in Element->CheckConsistency, her it is checked if node is used ...
				if (e->NodeNum(j) > 0 && e->NodeNum(j) <= nodescheck.Length())
				{
					nodescheck(e->NodeNum(j)) = 1;
				}
			}
		}
	}
	for (int i=1; i <= nodescheck.Length(); i++)
	{
		if (nodescheck(i) == 0)
		{
			UO(UO_LVL_err) << "System inconsistency in node number " << i << ":\n    Node is not used. You need to delete the node, otherwise the computation is not possible!!!\n";
			l1problem = 1;
		} 
	}

	if (l2problem) return 2;
	else if (l1problem) return 1;

	return 0;
}

int MultiBodySystem::CheckConsistencyMaterials()
{
	int l1problem = 0; //level 1 problem: can not compute
	int l2problem = 0; //level 2 problem: can not compute and not draw
	mystr errorstr = "";

	for (int i=1; i<=NMaterials(); i++) 
	{
		int rv = materials(i)->CheckConsistency(errorstr);
		if (rv == 1) l1problem = 1;
		else if (rv == 2) l2problem = 1;

		if (rv)
		{
			UO(UO_LVL_err) << "System inconsistency in material " << i << ":\n";
			UO(UO_LVL_err) << errorstr;
			errorstr = "";
			if (rv == 1) l1problem = 1;
			else if (rv == 2) l2problem = 1;
		}
	}

	if (l2problem) return 2;
	else if (l1problem) return 1;

	return 0;
}

int MultiBodySystem::CheckSystemConsistency() //check system; rv==0 --> OK, rv==1 --> can not compute, rv==2 --> can not draw and not compute
{
	//int l1problem = 0; //level 1 problem: can not compute
	//int l2problem = 0; //level 2 problem: can not compute and not draw

	//++++++++++++++++++++++++++++++++++++++++++++++++++++
	//first check if assemble and initialize is possible
	//check general inconsistencies
	//check LTG (due to delete or internal error)

	//MBSsos = GetSecondOrderSize();
	//MBSes = GetFirstOrderSize();
	//MBSis = GetImplicitSize();
	//MBSss = GetSystemSize();
//
//	for (int i=1; i<=elements.Length(); i++) 
//	{
//		Element* e = elements(i);
//		const	TArray<int>& ltg = e->GetLTGArray();
//		if (ltg.Length() != e->SS()) 
//		{
//			l2problem = 1;
//			UO(UO_LVL_err) << "System inconsistency: Element " << i << " has an inconsistency in the local to global DOF reference numbers!\n";
//		}
//
//		for (int j=1; j<=ltg.Length(); j++)
//		{
//// (AD) changed () to .Get()
//			if (ltg.Get(j) <= 0 || ltg.Get(j) > MBSss)
////			if (ltg(j) <= 0 || ltg(j) > MBSss)
//			{
//				l2problem = 1;
//				UO(UO_LVL_err) << "System inconsistency: Element " << i << " has invalid local to global DOF reference numbers!\n";
//			}
//		}
//	}
//
//	if (l2problem) return 2;

	int rv = CheckConsistencyLTG();
	if(rv==2) return 2;

	//now call element consistency check:
	//mystr errorstr = "";
	//for (int i=1; i<=elements.Length(); i++) 
	//{
	//	int rv = elements(i)->CheckConsistency(errorstr);
	//	if (rv == 1) l1problem = 1;
	//	else if (rv == 2) l2problem = 1;

	//	if (rv)
	//	{
	//		UO(UO_LVL_err) << "System inconsistency in element " << i << ":\n";
	//		UO(UO_LVL_err) << errorstr;
	//		errorstr = "";
	//		if (rv == 1) l1problem = 1;
	//		else if (rv == 2) l2problem = 1;
	//	}
	//}

	rv = Maximum(rv,CheckConsistencyElements());

	//for (int i=1; i <= NSensors(); i++)
	//{
	//	//$ YV 2012-06: for sensors true is consistent and false is not
	//	int rv = sensors(i)->IsConsistent(errorstr) ? 0 : 1;
	//	if (rv == 1) l1problem = 1;
	//	//else if (rv == 2) l2problem = 1;

	//	if (rv)
	//	{
	//		UO(UO_LVL_err) << "System inconsistency in sensor " << i << ":\n";
	//		UO(UO_LVL_err) << errorstr;
	//		errorstr = "";
	//		if (rv == 1) l1problem = 1;
	//		else if (rv == 2) l2problem = 1;
	//	}
	//}
	
	rv = Maximum(rv,CheckConsistencySensors());

	//check if all nodes are used:
	//TArray<int> nodescheck;
	//nodescheck.SetLen(NNodes());
	//for (int i=1; i <= nodescheck.Length(); i++)
	//{
	//	nodescheck(i) = 0;
	//}
	//for (int i=1; i<=elements.Length(); i++) 
	//{
	//	Element* e = elements(i);
	//	if (e->PerformNodeCheck())
	//	{
	//		for (int j=1; j <= e->NNodes(); j++)
	//		{
	//			//nodenumber is already verified in Element->CheckConsistency, her it is checked if node is used ...
	//			if (e->NodeNum(j) > 0 && e->NodeNum(j) <= nodescheck.Length())
	//			{
	//				nodescheck(e->NodeNum(j)) = 1;
	//			}
	//		}
	//	}
	//}
	//for (int i=1; i <= nodescheck.Length(); i++)
	//{
	//	if (nodescheck(i) == 0)
	//	{
	//		UO(UO_LVL_err) << "System inconsistency in node number " << i << ":\n    Node is not used. You need to delete the node, otherwise the computation is not possible!!!\n";
	//		l1problem = 1;
	//	} 
	//}
	rv = Maximum(rv,CheckConsistencyNodes());
	
	//$ PG 2012-12-19: added consistency check for materials
	rv = Maximum(rv,CheckConsistencyMaterials());

	if (rv==2) return 2;
	else if (rv==1) return 1;

	UO() << "System is consistent!\n";
	return 0;

	//if (l2problem) return 2;
	//else if (l1problem) return 1;

	//UO() << "System is consistent!\n";
	//return 0;

}

void MultiBodySystem::PrecomputeEvalFunctions()   //PG: precomputation for EvalF(), EvalF2(), EvalG(), EvalM(), and also for CqTLambda() in case of constraints
{
	//adjusted for IOElementDataModifier by JG 2013-17-01 and MSax 

	for (int i=1; i<=elements.Length(); i++) 
	{
		//first initialize all modifier elements, because these elements modify the data of other elements!
		if (elements(i)->GetType() & TIOElementDataModifier) elements(i)->PrecomputeEvalFunctions();
	}
	for (int i=1; i<=elements.Length(); i++) 
	{
		if (!(elements(i)->GetType() & TIOElementDataModifier)) elements(i)->PrecomputeEvalFunctions();
		//elements(i)->StartTimeStep();
	}


	//for (int i=1; i<=elements.Length(); i++)
	//{
	//	elements(i)->PrecomputeEvalFunctions();
	//}
}


void MultiBodySystem::EvalM(const Vector& x, Matrix& m, double t)
{		
	TMStartTimer(6);
	SetActState(x);

	//old: set all 0
	//m.SetAll(0);
	//new: clear only element entries:
	int maxbw = 0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int sos = e->SOS();
		if (sos)
		{
			const	TArray<int>& ltg = e->GetLTGArray();
			for (int j=1; j <= sos; j++)
			{
				for (int k=1; k <= sos; k++)
				{
// (AD) changed () to .Get()
					maxbw = Maximum(maxbw,ltg.Get(j)-ltg.Get(k)+1);
					m(ltg.Get(j),ltg.Get(k)) = 0;
//					maxbw = Maximum(maxbw,ltg(j)-ltg(k)+1);
//					m(ltg(j),ltg(k)) = 0;
				}
			}
		}
	}
	SetBWM(maxbw);
	//UO() << "max Bandwidth=" << maxbw << "\n";
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int sos = e->SOS();
		if (sos)
		{
			tempm.SetSize(sos,sos);
			tempm.FillWithZeros();
			e->EvalM(tempm,t);
			//uo << "tempm=" << tempm;

			m.AddMatrix(e->GetLTGArray(), e->GetLTGArray(), tempm);

			/*
			const	TArray<int>& ltg = e->GetLTGArray();
			for (int j=1; j <= sos; j++)
			{
			for (int k=1; k <= sos; k++)
			{
			//m(e->LTG(j),e->LTG(k)) += tempm(j,k);
			m(ltg(j),ltg(k)) += tempm(j,k);
			}
			}*/
		}
	}
	TMStopTimer(6);
	//uo << "m=" << m << "\n";
}

void MultiBodySystem::EvalM(const Vector& x, SparseMatrix& m, double t)
{
	TMStartTimer(6);
	SetActState(x);

	m.FillWithZeros();

	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);

		if (e->UseSparseM())
		{
			e->AddMSparse(m,t);
		}
		else
		{
			int sos = e->SOS();
			if (sos)
			{
				tempm.SetSize(sos,sos);
				tempm.FillWithZeros();
				e->EvalM(tempm,t);
				const	TArray<int>& ltg = e->GetLTGArray();
				m.AddMatrix(ltg,ltg,sos,sos,tempm);
			}
		}
	}
	TMStopTimer(6);
}

void MultiBodySystem::EvalF2(const Vector& x, Vector& f, double t)
{
	TMStartTimer(5);
	if (NLS_GetJacCol() == 0) {log.evalfcnt++;}
	if (NLS_GetJacCol() != 0) {log.evalf_jaccnt++;}

	int jc = GetJacCol();

	f.SetAll(0);
	SetActState(x); 
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int sos = e->SOS();
		if (sos)
		{
			if(e->IsDependent(jc))
			{
				tempf.SetLen(sos);
				tempf.SetAll(0);
				e->EvalF2(tempf,t);
				for (int j=1; j <= sos; j++)
				{
					f(e->LTG(j)) += tempf(j);
				}
			}
		}
	}
	//uo << "f2=" << f << "\n";
	TMStopTimer(5);
}

#ifndef MBS_use_parallelization
//for rapid Explicit integration (e.g. particles)
void MultiBodySystem::EvalMinvF2(const Vector& x, Vector& f, double t)
{
	TMStartTimer(5);
	if (NLS_GetJacCol() == 0) {log.evalfcnt++;}
	if (NLS_GetJacCol() != 0) {log.evalf_jaccnt++;}

	f.SetAll(0);
	SetActState(x); 
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int sos = e->SOS();
		if (sos)
		{
			tempf.SetLen(sos);
			tempf.SetAll(0);
			e->EvalMinvF2(tempf,t);
			for (int j=1; j <= sos; j++)
			{
				f(e->LTG(j)) += tempf(j);
			}
		}
	}
	//uo << "f2=" << f << "\n";
	TMStopTimer(5);
}
#endif

#ifdef MBS_use_parallelization
//for rapid Explicit integration (e.g. particles)
void MultiBodySystem::EvalMinvF2(const Vector& x, Vector& f, double t)
{
	TMStartTimer(5);
	if (NLS_GetJacCol() == 0) {evalfcnt++;}
	if (NLS_GetJacCol() != 0) {evalf_jaccnt++;}

	f.SetAll(0);
	SetActState(x); 

	//omp_set_num_threads(MBSSolSet().parallel_number_of_threads);
	omp_set_num_threads(2);

	int i,j;

	Vector tempf[OMP_MAX_N_THREADS];
	int nt;
	int ne = elements.Length();
	//#pragma omp parallel for private(nt,j,sos,e)
/*
	double a[OMP_MAX_N_THREADS];
	a[0]=0;a[1]=0;a[2]=0;a[3]=0;
	double b=0;
	int n = 20000;


//#pragma omp parallel shared(tempf, elements, f, t, ne, a, b) private(i,j,nt)
#pragma omp parallel for shared(a) private(nt,j)
	for (i=1; i<=n; i++) 
	{
	  nt = omp_get_thread_num();
		for (j=1; j<= 100; j++)
		{
			a[nt] += sin((double)j*(double)i);
		}
	}
*/


#pragma omp parallel shared(tempf, f, t, ne) private(nt, j)
	{
#pragma omp for
		for (i=1; i<=ne; i++) 
		{
			nt = omp_get_thread_num();
			if (elements(i)->SOS())
			{
				tempf[nt].SetLen(elements(i)->SOS());
				tempf[nt].SetAll(0);
				elements(i)->EvalMinvF2(tempf[nt],t);
				for (j=1; j <= elements(i)->SOS(); j++)
				{
					f(elements(i)->LTG(j)) += tempf[nt](j);
				}
			}
		}
	} //End of parallel region 

	/*
	for (i=1; i<=ne; i++) 
	{
		int nt = 0;//omp_get_thread_num();
		if (elements(i)->SOS())
		{
			tempf[nt].SetLen(elements(i)->SOS());
			tempf[nt].SetAll(0);
			elements(i)->EvalMinvF2(tempf[nt],t);
			//tempf[nt] *= a[nt];
			for (j=1; j <= elements(i)->SOS(); j++)
			{
				f(elements(i)->LTG(j)) += tempf[nt](j);
			}
		}
	}
*/
	//uo << "f2=" << f << "\n";
	TMStopTimer(5);
}
#endif

void MultiBodySystem::EvalF2(const Vector& x, SparseVector& f, double t, IVector& rowind, IVector& clearind, IVector& elems)
{
	//rowind has f.Length() is initialized with zeros and returned with zeros!
	//clearind is arbitrary and contains the indices of rowind to be set 0
	//rowind and clearind is used to fastly access the sparsevector f

	//elems is used to see if the second vector for numerical differentiation needs to be computed

	TMStartTimer(5);
	TMStartTimer(9);

	clearind.SetLen(0);
	f.FillWithZeros();
	elems.SetLen(0);

	if (NLS_GetJacCol() == 0) {log.evalfcnt++;}
	if (NLS_GetJacCol() != 0) {log.evalf_jaccnt++;}
	int ltg;

	int jc = GetJacCol();


	for (int i=1; i<=elements.Length(); i++) 
	{
		const Element& e = *elements(i);
		int sos = e.SOS();
		if (sos)
		{
			if(e.IsDependent(jc))
			{
				elems.Add(i);
				for (int j=1; j <= sos; j++)
				{
					ltg = e.LTG(j);
					if (!rowind(ltg))
					{
						clearind.Add(ltg);
						f.AddEntry(ltg,0);
						rowind(ltg) = f.NEntries();
					}
				}
			}
		}
	}
	TMStopTimer(9);

	SetActState(x); 
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int sos = e->SOS();
		if (sos)
		{
			if(e->IsDependent(jc))
			{
				tempf.SetLen(sos);
				tempf.SetAll(0);
				e->EvalF2(tempf,t);
				for (int j=1; j <= sos; j++)
				{
					f.Entry(rowind(e->LTG(j))) += tempf(j);
					//f(e->LTG(j)) += tempf(j);
				}
			}
		}
	}

	for (int i=1; i <= clearind.Length(); i++) rowind(clearind(i)) = 0;

	TMStopTimer(5);
}

#define SpecialAdd Add //AddIfNotExists
//compute entries of each row in sparse matrix
void MultiBodySystem::ComputeSparseMatrixUsage(IVector& usage_per_dof) 
{
	int sos= GetSecondOrderSize();//Number of 2nd order ODEs
	int es = GetFirstOrderSize(); //Number of 1st order ODEs
	int is = GetImplicitSize();   //Algebraic Equations (from boundary conditions) 
	int matsize = sos+es+is;
	int init_col_size=4; //initial size of columns

	usage_per_dof.SetLen(matsize);
	usage_per_dof.SetAll(0);

	TArray<IVector*> usage_per_dof_matrix(sos+es+is);
	for (int i=1; i<=matsize; i++)
	{
		usage_per_dof_matrix(i) = new IVector(init_col_size);
	}

	//this function is called before computation of sparse stiffness or mass matrix
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);

		int elem_sos = e->SOS();
		int elem_es = e->ES();
		int elem_is = e->IS();
		int elem_ss = elem_sos+elem_es+elem_is;

		const	TArray<int>& ltgrows = e->GetLTGArray();

		//UO(UO_LVL_all) << "elem" << i << ": n=" << ltgrows.Length() << ", ltg=" << ltgrows << "\n"; 

		//classical stiffness terms:
		for (int j=1; j <= elem_sos; j++)
		{
			int ind = ltgrows(j);
			usage_per_dof(ind) += elem_ss;

			IVector* updm = usage_per_dof_matrix(ind);
			for (int k=1; k <= elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}
			for (int k=elem_sos*2+1; k <= elem_ss+elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}
		}
		//stiffness terms from first order equations:
		for (int j=1; j <= elem_es; j++)
		{
			int ind = ltgrows(elem_sos*2+j)-sos;
			usage_per_dof(ind) += elem_ss;

			IVector* updm = usage_per_dof_matrix(ind);
			for (int k=1; k <= elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}
			for (int k=elem_sos*2+1; k <= elem_ss+elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}
		}
		//constraint jacobian matrices:
		for (int j=1; j <= elem_is; j++)
		{
			int ind = ltgrows(elem_sos*2+elem_es+j)-sos;
			usage_per_dof(ind) += elem_ss;

			IVector* updm = usage_per_dof_matrix(ind);
			for (int k=1; k <= elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}
			for (int k=elem_sos*2+1; k <= elem_ss+elem_sos; k++)
			{
				updm->SpecialAdd(ltgrows(k));
			}

			for (int k=1; k <= e->NE(); k++)
			{
				//UO(UO_LVL_all) << "  elemnum=" << e->GetElnum(k) << "\n";
				Element* ce = GetElementPtr(e->GetElnum(k));
				int celem_sos = ce->SOS();
				int celem_es = ce->ES();
				int celem_is = ce->IS();
				int celem_ss = celem_sos+celem_es+celem_is;

				int ind = ltgrows(elem_sos*2+elem_es+j)-sos;
				usage_per_dof(ind) += celem_ss;

				const	TArray<int>& cltgrows = ce->GetLTGArray();
				IVector* updm = usage_per_dof_matrix(ind);
				for (int k=1; k <= celem_sos; k++)
				{
					updm->SpecialAdd(cltgrows(k));
				}
				for (int k=celem_sos*2+1; k <= celem_ss+celem_sos; k++)
				{
					updm->SpecialAdd(cltgrows(k));
				}
			}
		}
		if (e->IsType(TConstraint))
		{
			//UO(UO_LVL_all) << "elem" << i << ", NE=" << e->NE() << "\n";
			for (int k=1; k <= e->NE(); k++)
			{
				//UO(UO_LVL_all) << "  elemnum=" << e->GetElnum(k) << "\n";
				Element* ce = GetElementPtr(e->GetElnum(k));
				int celem_sos = ce->SOS();
				int celem_es = ce->ES();
				int celem_is = ce->IS();
				int celem_ss = celem_sos+celem_es+celem_is;

				const	TArray<int>& cltgrows = ce->GetLTGArray();

				for (int m=1; m <= celem_sos; m++)
				{
					int ind = cltgrows(m);
					usage_per_dof(ind) += elem_is;

					IVector* updm = usage_per_dof_matrix(ind);
					for (int k=elem_sos*2+elem_es+1; k <= elem_ss+elem_sos; k++)
					{
						updm->SpecialAdd(ltgrows(k));
					}
				}
				for (int m=1; m <= celem_es; m++)
				{
					int ind = cltgrows(celem_sos*2+m)-sos;
					usage_per_dof(ind) += elem_is;

					IVector* updm = usage_per_dof_matrix(ind);
					for (int k=elem_sos*2+elem_es+1; k <= elem_ss+elem_sos; k++)
					{
						updm->SpecialAdd(ltgrows(k));
					}
				}
				for (int m=1; m <= celem_is; m++)
				{
					int ind = cltgrows(celem_sos*2+celem_es+m)-sos;
					usage_per_dof(ind) += elem_is;

					IVector* updm = usage_per_dof_matrix(ind);
					for (int k=elem_sos*2+elem_es+1; k <= elem_ss+elem_sos; k++)
					{
						updm->SpecialAdd(ltgrows(k));
					}
				}
			}
		}

	}

	int totalsize = 0;
	int totalsize2 = 0;
	for (int i=1; i <= usage_per_dof.Length(); i++)
	{
		totalsize += usage_per_dof(i);
		RemoveRedundantEntries(*usage_per_dof_matrix(i), 1);
		totalsize2 += usage_per_dof_matrix(i)->Length();

		usage_per_dof(i) = usage_per_dof_matrix(i)->Length();
	}

	//UO(UO_LVL_all) << "ComputeSparseMatrixUsage=" << totalsize << "\n";

	//UO(UO_LVL_all) << "ComputeSparseMatrixUsage_Exact=" << totalsize2 << "\n";

	//delete allocated arrays:
	for (int i=1; i<=matsize; i++)
	{
		delete usage_per_dof_matrix(i);
	}
}

void MultiBodySystem::ComputeGyroscopicMatrix(SparseMatrix& gy) // $ MSax 2013-07-25 : added
{
	//int is = GetImplicitSize();
	//int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();
	//int ss = GetSystemSize()

	SparseMatrix locgy;
	Matrix fulllocgy;
	gy.SetSize(sos,sos,MaxSparseBandwidth()*2);
	gy.FillWithZeros();
		
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);

		int sos = e->SOS();
		if (sos)
		{
			e->GyroscopicMatrix(locgy); 
			if (locgy.CountEntries() != 0) // fill in to global gyroscopic matrix (only if not empty)
			{
				const	TArray<int>& ltgrows = e->GetLTGArray();
				locgy.GetMatrix(fulllocgy);
				gy.AddMatrix(ltgrows,ltgrows,sos,sos,fulllocgy);
			}
		}
	}
}

void MultiBodySystem::LocalJacobianF2(SparseMatrix& m, Vector& x)
{
	int mode = 2; //1=old mode, 2=element-wise mode
	if ((!numsol.ModifiedNewtonActive() || !GetSolSet().nls_modifiednewton) && !(GetSolSet().nls_usesparsesolver)) 
	{
		mode = 1;
	}

	if (mode == 1) UO() << "*****************\nwarning: switched MBS:LocalJacobianF2 to old mode\n******************\n";
  
	double t = GetTime() + 0.5*GetStepSize();
	//UO() << "localjacf2\n";

	if (mode == 1)
	{
		TimeInt::LocalJacobianF2(m,x);
	}
	else
	{
		TMStartTimer(5); // otherwise EvalF2 not right counted
		SetActState(x);

		m.FillWithZeros();

		for (int i=1; i<=elements.Length(); i++) 
		{
			Element* e = elements(i);

			int sos = e->SOS();
			if (sos)
			{
				if (e->UseSparseK())
				{
					e->AddKSparse(m,t); //follower forces are missing?!!!

					//e->AddDSparse(m,t); 

					//only constraint part:
					e->JacobianF2(t,tempm,colref);

					const	TArray<int>& ltgrows = e->GetLTGArray();
					m.AddMatrix(ltgrows,colref,sos,colref.Length(),tempm);
				}
				else
				{

					e->JacobianF2(t,tempm,colref);
					const	TArray<int>& ltgrows = e->GetLTGArray();
					m.AddMatrix(ltgrows,colref,sos,colref.Length(),tempm);
				}
			}
		}
		TMStopTimer(5);

		//int n = m.CountEntries();
		//UO() << "K-alloc=" << m.GetLAlloc()*8./1.e6 << "MB, entries=" << (double)n/1e6 << "\n";
		
	}

}

void MultiBodySystem::StaticJacobianF2(SparseMatrix& m, Vector& x)
{
	double t = GetTime() + 0.5*GetStepSize();
	//UO() << "localjacf2\n";

	TMStartTimer(5); // otherwise EvalF2 not right counted
	SetActState(x);

	int fullsos = GetSecondOrderSize();

	m.FillWithZeros();
	//TArray<int> ltg_static; //this is the ltg-array, but velocity part is deleted
	TArray<int> colref_static; //this is the colref-array, but velocity part is deleted
	Matrix tempmnew;

	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);

		int sos = e->SOS();
		int es = e->ES();
		int is = e->IS();
		if (sos)
		{
			//the following function should be extended by e->UseSparseK(), see LocalJacobianF2(sparse...
			if (e->UseSparseK())
			{
				assert(0);
			}

			colref_static.SetLen(0);

			e->JacobianF2(t,tempm,colref);
			const	TArray<int>& ltgrows = e->GetLTGArray();
			for (int j=1; j <= colref.Length(); j++)
			{
				if (colref(j) <= fullsos) 
				{
					colref_static.Add(colref(j)); //first part of colref should be same as e->ltg
				}
				else if (colref(j) > 2*fullsos) 
				{
					colref_static.Add(colref(j)-fullsos); //first part of colref should be same as e->ltg
				}
			}
			//for (int j=2*sos+1; j <= colref.Length(); j++)
			//{
			//	colref_static.Add(colref(j)-fullsos);
			//}
			//ltg_static not needed!
			tempmnew.SetSize(tempm.Getrows(), tempm.Getcols()-sos);
			tempmnew.SetAll(0.);
			tempmnew.AddSubmatrix(tempm,1,1,1,1,tempm.Getrows(), sos, 1.);
			tempmnew.AddSubmatrix(tempm,1,1+2*sos,1,1+sos,tempm.Getrows(), tempm.Getcols()-2*sos, 1.);

			m.AddMatrix(ltgrows,colref_static,sos,colref_static.Length(),tempmnew);
		}
	}
	TMStopTimer(5);
}

void MultiBodySystem::LocalJacobianF2(Matrix& m, Vector& x)
{
	int mode = 2; //1=old mode, 2=element-wise mode
	if (!UseJacobianElementWise()) mode = 1;

	double t = GetTime() + 0.5*GetStepSize();

	//UO() << "localjacf2\n";

	if (mode == 1)
	{
		TimeInt::LocalJacobianF2(m,x);
	}
	else
	{
		SetActState(x);

		m.FillWithZeros();

		TMStartTimer(5); // otherwise EvalF2 not right counted
		for (int i=1; i<=elements.Length(); i++) 
		{
			Element* e = elements(i);
			int sos = e->SOS();
			if (sos)
			{
				e->JacobianF2(t,tempm,colref);
				//UO() << "elem=" << i << ": JacF2=" << tempm << "\n";
				const	TArray<int>& ltgrows = e->GetLTGArray();

				for (int j=1; j <= sos; j++)
				{
					for (int k=1; k <= colref.Length(); k++)
					{
// (AD) changed () to .Get()
						m(ltgrows.Get(j),colref.Get(k)) += tempm(j,k);
//						m(ltgrows(j),colref(k)) += tempm(j,k);
					}
				}
			}
		}
		TMStopTimer(5);
	}
}




void MultiBodySystem::LocalJacobianG(SparseMatrix& m, Vector& x)
{
	int mode = 2; //1=old mode, 2=element-wise mode
	double t = GetTime() + 0.5*GetStepSize();
	if (mode == 1)
	{
		TimeInt::LocalJacobianG(m,x);
	}
	else
	{
		static IVector rowref(100);
		int mbsoff = GetSecondOrderSize()*2+GetFirstOrderSize();

		SetActState(x);

		m.FillWithZeros();

		TMStartTimer(8); // otherwise EvalF2 not right counted
		for (int i=1; i<=elements.Length(); i++) 
		{
			rowref.SetLen(0);
			Element* e = elements(i);
			int is = e->IS();
			if (is)
			{
				e->JacobianG(t,tempm,colref);
				const	TArray<int>& ltgrows = e->GetLTGArray();
				for (int j=1; j <= is; j++)
				{
					rowref.Add(ltgrows.Get(j+2*e->SOS()+e->ES())-mbsoff);
				}
				m.AddMatrix(rowref,colref,is,colref.Length(),tempm);
			}
		}
		TMStopTimer(8);
	}

}




void MultiBodySystem::EvalF(const Vector& x, Vector& f, double t)
{
	TMStartTimer(7);
	f.SetAll(0);
	int off = -2*MBSsos;
	SetActState(x);
	int jc = GetJacCol();
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int es = e->ES();
		if (es)
		{
			if (e->IsDependent(jc))
			{
				int off2=2*e->SOS();
				tempf.SetLen(es);
				tempf.SetAll(0);
				e->EvalF(tempf,t);
				for (int j=1; j <= es; j++)
				{
					f(e->LTG(j+off2)+off) = tempf(j);
				}
			}
		}
	}
	TMStopTimer(7);
}

void MultiBodySystem::EvalG(const Vector& x, Vector& f, double t)
{
	TMStartTimer(8);
	int jc = GetJacCol();

	f.SetAll(0);
	int off = -2*MBSsos-MBSes;
	SetActState(x);
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int is = e->IS();
		if (is)
		{
			if (e->IsDependent(jc))
			{
				int off2=2*e->SOS()+e->ES();
				tempf.SetLen(is);
				tempf.SetAll(0);
				e->EvalG(tempf,t);
				for (int j=1; j <= is; j++)
				{
					f(e->LTG(j+off2)+off) = tempf(j);
				}
			}
		}
	}
	//uo << "g=" << f << "\n";

	//f*=-2./(GetStepSize()); //2./StepSize() originates from KVersion time integration with acceleration instead of position, sign originates because J=-M-tau^2/4*K-CqT L

	TMStopTimer(8);
}

void MultiBodySystem::EvalG(const Vector& x, SparseVector& f, double t, IVector& rowind, IVector& clearind)
{
	//rowind has f.Length() is initialized with zeros and returned with zeros!
	//clearind is arbitrary and contains the indices of rowind to be set 0
	TMStartTimer(8);
	TMStartTimer(9);

	clearind.SetLen(0);
	f.FillWithZeros();

	int off = -2*MBSsos-MBSes;

	//if (NLS_GetJacCol() == 0) {evalfcnt++;}
	//if (NLS_GetJacCol() != 0) {evalf_jaccnt++;}

	int ltg, off2;

	int jc = GetJacCol();

	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int is = e->IS();
		if (is)
		{
			if (e->IsDependent(jc))
			{
				off2=2*e->SOS()+e->ES();
				tempf.SetLen(is);
				tempf.SetAll(0);
				e->EvalG(tempf,t);
				for (int j=1; j <= is; j++)
				{
					ltg = e->LTG(j+off2)+off;
					clearind.Add(ltg);
					f.AddEntry(ltg,0);
					rowind(ltg) = f.NEntries();
				}
			}
		}
	}
	TMStopTimer(9);

	SetActState(x); 
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		int is = e->IS();
		if (is)
		{
			if(e->IsDependent(jc))
			{
				off2=2*e->SOS()+e->ES();
				tempf.SetLen(is);
				tempf.SetAll(0);
				e->EvalG(tempf,t);
				for (int j=1; j <= is; j++)
				{
					f.Entry(rowind(e->LTG(j+off2)+off)) += tempf(j);
				}
			}
		}
	}

	for (int i=1; i <= clearind.Length(); i++) rowind(clearind(i)) = 0;

	//for (int i=1; i <= f.NEntries(); i++) f.Entry(i) *= -1;

	TMStopTimer(8);
}

double MultiBodySystem::GetKineticEnergy()
{
	double en=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		en += elements(i)->GetKineticEnergy();
	}
	return en;
}

double MultiBodySystem::GetPotentialEnergy()
{
	double en=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		en += elements(i)->GetPotentialEnergy();
	}
	return en;
}



int MultiBodySystem::GetImplicitSize() const
{
	int s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->IS();
	}
	return s;
}

int MultiBodySystem::GetSecondOrderSize() const
{
	int s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->SOSowned();
	}
	for (int i=1; i<=nodes.Length(); i++) 
	{
		s+=nodes(i)->SOS();
	}
	return s;
}

int MultiBodySystem::GetSecondOrderSize_RS() const
{
	int s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->SOSowned_RS();
	}
	for (int i=1; i<=nodes.Length(); i++) 
	{
		s+=nodes(i)->SOS();
	}
	return s+resortsize_add;
}

int MultiBodySystem::GetFirstOrderSize() const
{
	int s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->ES();
	}
	return s;
}

int MultiBodySystem::GetResortSize() const
{
	int s=0;
	if (!ResortActive()) return 0;

	for (int i=1; i<=elements.Length(); i++) 
	{
		if (ResortMode2())
		{
			s += abs((elements(i)->SOSowned_RS()-elements(i)->SOSowned()));
		}
		else
		{
			s += (elements(i)->SOSowned_RS()-elements(i)->SOSowned());
		}
	}
	return s+resortsize_add;
}

int MultiBodySystem::GetDataSize() const //size of data variables for each time step of system
{
	int s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->DataS();
	}
	//$ MaSch 2013-10-28: added support of SPH particle data via auxiliary elements
	for (int i=1; i<=auxelements.Length(); i++) 
	{
		if(auxelements(i)->IsType(TParticle))
		{
			s+=auxelements(i)->DataS();
		}
	}
	return s;
}

double MultiBodySystem::PostNewtonStep(double t) 
{
	if (SolverPrintsDetailedOutput())
	{
		SolverUO() << "-->PostNewtonStep " << log.TInonlinit << ":";
	}

	TMStartTimer(20);
	double s=0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		s+=elements(i)->PostNewtonStep(t);
	}
	TMStopTimer(20);

	if (SolverPrintsDetailedOutput())
	{
		mystr str(" err=");
		str += mystr(s);
		str += mystr(", goal=");
		str += mystr(DiscontinuousAccuracy());
		str += mystr(", t=");
		str += mystr(t);
		str += mystr("\n");
		SolverUO() << str;

		if (GetOptions()->LoggingOptions()->SolverPostNewtonIterationDataVector())
		{
			SolverUO() << "Data vector = " << GetDataVector() << "\n";
		}
	}

	return s;
};

void MultiBodySystem::StartTimeStep() 
{
	//adjusted for IOElementDataModifier by JG 2013-17-01 and MSax 

	for (int i=1; i<=elements.Length(); i++) 
	{
		//first initialize all modifier elements, because these elements modify the data of other elements!
		if (elements(i)->GetType() & TIOElementDataModifier) elements(i)->StartTimeStep();
	}
	for (int i=1; i<=elements.Length(); i++) 
	{
		if (!(elements(i)->GetType() & TIOElementDataModifier)) elements(i)->StartTimeStep();
		//elements(i)->StartTimeStep();
	}
};

void MultiBodySystem::EndTimeStep() 
{
	for (int i=1; i<=elements.Length(); i++) 
	{
		elements(i)->EndTimeStep();
	}
};

void MultiBodySystem::ComputationFinished() 
{
	for (int i=1; i<=elements.Length(); i++) 
	{
		elements(i)->ComputationFinished();
	}
};

void MultiBodySystem::PostprocessingStep() 
{
	SetActState(GetSolVector());
	for (int i=1; i<=elements.Length(); i++) 
	{
		elements(i)->PostprocessingStep();
	}
};

double MultiBodySystem::GetError()
{
	SetActState(GetSolVector());
	double err=0;
	for (int i=1; i <= elements.Length(); i++) 
	{
		err+=elements(i)->GetError();
	}
	int numc = MBSss-MBSis;
	if (numc!=0) err = sqrt(err)/(double)numc;
	else err = sqrt(err);

	return err;
}

ElementDataContainer* MultiBodySystem::GetComputationStepEDC(int computation_step_number)
{
	ElementDataContainer* edc = GetMBS_EDC_Options();
	ElementDataContainer* CSEDC;
	ElementData* ed = edc->TreeFind("SolverOptions.ComputationSteps");
	if(ed && ed->IsEDC() )
	{
		CSEDC = ed->GetEDC();
	}
	else 
		return NULL; // no EDC ComputationSteps

	mystr step_tag = mystr("Step") + mystr(computation_step_number);
  ed = CSEDC->TreeFind(step_tag);
	if(ed && ed->IsEDC())
	{
    return ed->GetEDC();
	}
	else
		return NULL; // no EDC Step#
}

int MultiBodySystem::ParseStepEndTimes()
{
	CSEndTimes.Flush();

// declaration & first iteration 
	ElementDataContainer* stepEDC = GetComputationStepEDC(1);
	double stependtime = 0.;
	double laststependtime = 0.;

	while (stepEDC)
	{
		stependtime = stepEDC->TreeGetDouble("step_end_time");
		if(stependtime > laststependtime) // thus all step_end_times > 0.
		{
			CSEndTimes.Add(stependtime);

// next iteration
			laststependtime = stependtime;
			stepEDC = GetComputationStepEDC(CSEndTimes.Length()+1);
		}
		else
		{
			UO().InstantMessageText(mystr("step_end_time of Step") + mystr(CSEndTimes.Length()+1) + mystr("is NOT higher then in previous step!\n parsing stopped, all following Computation steps are not accounted for"));
			break;
		}
	}
// in case no computation step is defined ... 
	if(CSEndTimes.Length() == 0)
	{
		CSEndTimes.Add(solset.endtime);
	}

	return CSEndTimes.Length();
}

int MultiBodySystem::DoubleCheckWithLoadSettings()
{
	int nr_cs = CSEndTimes.Length();
	int nr_ls = GetMaxLoadSteps();

	while (CSEndTimes.Length() < nr_ls)
	{
		AddDummyComputationStep();
		nr_cs = CSEndTimes.Length();
	}
	return nr_cs;
}

// maximal number of loadsteps from all loads - there MUST be at least that many computation steps...
int MultiBodySystem::GetMaxLoadSteps() const
{
	int max = 0;

	for(int i=1; i <= GetNElements(); i++)
	{
		const Element& el=GetElement(i);
		
		for (int j=1; j <= el.NLoads(); j++)
		{
			const MBSLoad& load = el.GetLoad(j);
			if(max < load.NLoadSteps())
				max=load.NLoadSteps();
		}
	}
	return max;
}

// additional computationsteps if loadsteps supercede computation steps, interval always 1 sec*
int MultiBodySystem::AddDummyComputationStep()
{
  int nr_cs = CSEndTimes.Length();               // existing computation step entries
	double endtime = 1;                            
	if (nr_cs) endtime = CSEndTimes.Last() + 1.;  // end time

  CSEndTimes.Add(endtime);

// entry in EDC - entries not temporary !!!
	ElementDataContainer* edc = GetMBS_EDC_Options();
	mystr tree = mystr("SolverOptions.ComputationSteps.Step") + mystr(nr_cs+1);
	ElementData ed;
	ed.SetDouble(endtime,"step_end_time");
	edc->TreeAdd(tree,ed);
  
	return CSEndTimes.Length();
}


int MultiBodySystem::ApplyComputationStepSettings(int stepnumber)
{
	UO(UO_LVL_ext) << mystr("*** Computation Step ") + mystr(stepnumber) + mystr("@ t= ") + mystr(GetTime()) + mystr(" ***\n");
	// get step EDC
	ElementDataContainer* stepEDC = GetComputationStepEDC(stepnumber);
  
	if(stepEDC != NULL)
	{
		ElementDataContainer edc;
// solset -> EDC
		SolverOptions2EDC(&edc);	
// EDC + changes
		ElementData* ed_solveroptions = edc.TreeFind("SolverOptions");     // variant a) pick subedc "SolverOptions" to write to     rather then insert a level "SolverOptions" the step EDC ...
		if(ed_solveroptions && ed_solveroptions->IsEDC() )
		{
		//$ AD: 2011-11-08: add step_end_time manually to prevent popups
			ElementData* dummy = stepEDC->TreeFind("step_end_time");
			ed_solveroptions->GetEDC()->Add(*dummy);
			ed_solveroptions->GetEDC()->TreeReplaceEDCDataWith(stepEDC,1);
		}
		else
			return -1; // rv "-1" for "no SolverOptions found, no changes applied"
// EDC -> solset
		EDC2SolverOptions(&edc);

		return stepnumber;
	}
	else
		return 0; // rv "0" for "no EDC found, no changes applied"
}


// add sensors for all "TController"-elements
void MultiBodySystem::AddControllerSensors()
{
	/*
	for(int i=1; i <= GetNElements(); i++)
	{
		const Element& el=GetElement(i);
			
		if(el.IsType(TController))
		{
			UO() << "elem" << i << ", no=" << ((InputOutputElement&)el).GetNOutputs() << "\n";
			for(int j=1; j <= ((InputOutputElement&)el).GetNOutputs(); j++)
			{
				MBSSensor s(this, TMBSSensor(TSOutputSensor), i, j);					
				s.SetSensorName(mystr("ControllerElement") + mystr(i) + mystr("_output") + mystr(j)); 
				AddSensor(&s);		
				UO() << "Add controller sensor (elementNum = " << i << ", outputNum = " << j << "\n";
			}
		}
	}	
	*/
}

//$ YV 2011-06-07:[ placed this piece of code into a separate function
void MultiBodySystem::BuildLTGLists()
{
	//reset nodes and elements LTG, for multiple call of assemble!
	for (int i=1; i<=nodes.Length(); i++) 
	{
		GetNode(i).LTGreset();
		GetNode(i).ResetNodeToElementList();
	}

	int warn1 = 0;
	for (int i=1; i<=elements.Length(); i++) 
	{
		//build node_to_element lists in nodes
		for (int j=1; j<=GetElement(i).NNodes(); j++)
		{
			//if (!GetElement(i).GetNode(j).IsAuxNode())
			GetElement(i).GetNode(j).AddElementNumber(i);
		}
		//initialize LTG lists
		/*if (GetElement(i).IsType(TCMSflag) && !warn1) 
		{
			warn1 = 1;
			UO() << "Assemble: check if still works with CMS elements!!!\n";
		}*/
		GetElement(i).LTGreset();
		GetElement(i).RemoveConstraints();
		GetElement(i).LTGdataReset();
	}

	for (int i=1; i<=auxelements.Length(); i++) 
	{
		for (int j=1; j<=GetAuxElement(i).NNodes(); j++)
		{
			//if (!GetElement(i).GetNode(j).IsAuxNode())
			GetAuxElement(i).GetNode(j).AddElementNumber(i);
		}
	}

	//These numbers are offsets!!!!!
	int sos1 = 1;
	int sos2 = GetSecondOrderSize()+1;
	int es = GetSecondOrderSize()*2+1;
	int is = GetSecondOrderSize()*2+GetFirstOrderSize()+1;
	int datas = 1;

	//DOF-references
	for (int i=1; i<=nodes.Length(); i++) 
	{
		for (int j=1; j <= GetNode(i).SOS(); j++)
		{
			GetNode(i).AddLTG(sos1++);
		}
		for (int j=1; j <= GetNode(i).SOS(); j++)
		{
			GetNode(i).AddLTG(sos2++);
		}
	}

	for (int i=1; i<=elements.Length(); i++) 
	{
		//do not clear ltg-lists, because elements might have already built up lists (CMS-element)??
		Element* e = elements(i);
		for (int j=1; j <= e->SOSowned(); j++)
		{
			e->AddLTG(sos1++);
		}
		for (int j=1; j <= e->SOSowned(); j++)
		{
			e->AddLTG(sos2++);
		}
		for (int j=1; j <= e->ES(); j++)
		{
			e->AddLTG(es++);
		}
		for (int j=1; j <= e->IS(); j++)
		{
			e->AddLTG(is++);
		}
		for (int j=1; j <= e->DataS(); j++)
		{
			e->AddLTGdata(datas++);
		}
		e->LinkToElements();

		e->LinkLoads();
		//uo << "ltg" << i << "=" << e->LTG(1) << ", " << e->LTG(2) << "\n";
	}

	//$ MaSch 2013-10-28: added support of SPH particle data via auxiliary elements
	for (int i=1; i<=auxelements.Length(); i++) 
	{
		Element* e = auxelements(i);
		if(e->IsType(TParticle))
		{
			for (int j=1; j <= e->DataS(); j++)
			{
				e->AddLTGdata(datas++);
			}
		}
	}


}
//$ YV 2011-06-07:]

//create LTG lists:
void MultiBodySystem::Assemble()
{
	int info = 0;
	if (info) UO() << "start assemble\n";

	for (int i=1; i<=elements.Length() + auxelements.Length(); i++)
	{
		Element * e;
		if(i <= elements.Length())
			e = &GetElement(i);
		else
			e = &GetAuxElement(i - elements.Length());
		e->PreAssemble();
	}
	
	BuildLTGLists();
	
	if (UseDependencies())
	{
		for (int i=1; i<=elements.Length(); i++) 
		{
			Element* e = elements(i);
			e->BuildDependencies();
		}
	}
	else
	{
	//	uo << "**************\nWARNING: dependencies not built!!!\n**************\n";
	}

	LabelIOElementOutputs();   // 
  FilterIOElements();        // create a shortlist of all IO elements in the MBS

	//$ YV 2011-06-07:[ the following line is commented out, as resorting is not needed any more
	//ComputeResortVector();
	
	//set outer faces for drawing
	//$ YV 2012-12-13: setting outer faces for drawing has been moved to FiniteElement3D::PreAssemble()

	CollectAvailableFieldVariables();
	if (info)
		UO() << "finished assembling\n";
}

//$ AD: new 2013-02-01 see Log 388
void MultiBodySystem::LabelIOElementOutputs()
{
// output to Load
	for (int i=1; i<=NLoads(); i++)
	{
		if (GetLoad(i).HasIOElement())
		{
			int elemnr = -1;
			int outputnr = -1;
			GetLoad(i).GetIOElementNr(elemnr, outputnr);
			if (elemnr > 0 && outputnr > 0)
			{
				if (elemnr <= NE())
				{
					if (GetElement(elemnr).IsType(TMBSElement(TController+TConstraint)))
					{
						//InputOutputElement* ioelem = (InputOutputElement*)GetElementPtr(elemnr); // does not agree with "interface idea" ...
						Element* ioelem = GetElementPtr(elemnr); 

						if (outputnr <= ioelem->GetNOutputs())
						{
							mystr str = mystr("L") + mystr(i) + mystr(": ") + GetLoad(i).LoadName();
//							mystr str = mystr("L") + mystr(i) + mystr(": ") + GetLoad(i).GetLoadTypeString();
							ioelem->SetOutputName(outputnr, str);
						}
						else
							UO() << mystr("ERROR: Load ") + mystr(i) + mystr(": referenced output does not exist on element. \n");
					}
					else
						UO() << mystr("ERROR: Load ") + mystr(i) + mystr(": referenced element is not an IOElement. \n");
				}
				else
					UO() << mystr("ERROR: Load ") + mystr(i) + mystr(" has an illegal entry for IOElement number. \n");
			}
			else
				UO() << mystr("ERROR: Load ") + mystr(i) + mystr(" has an illegal entry for IOElement number and/or output number. \n");
		}
	}
// Sensor
}

int MultiBodySystem::FilterIOElements()
{
	ioelements.Flush();
	for (int i=1; i<=NE(); i++)
	{
		if (GetElement(i).IsType(TMBSElement(TController+TConstraint)))
		{
			ioelements.Add(i);
		}
	}
	return ioelements.Length();
}

int MultiBodySystem::RespondToKey(int key)
{
	int nr_changes = 0;
	for (int i=1; i<=ioelements.Length(); i++)
	{
		nr_changes += GetElement(ioelements(i)).RespondToKey(key); 
	}
	return nr_changes;
}

//void MultiBodySystem::RespondToMouse(double x, double y)
//{
//
//}



//Vector for resorting the band mass- and stiffness matrix such that it is regular
//-> for quaternion elements, the constraints are resorted to the SOS part!!!
void MultiBodySystem::ComputeResortVector()
{
	//These numbers are offsets!!!!!
	int rs = GetResortSize();
	int sos = GetSecondOrderSize();
	int sos_rs = GetSecondOrderSize_RS();
	int es = GetFirstOrderSize();
	int is = GetImplicitSize();
	int ne = elements.Length();

	//the resorted vector has only components sos+es+is, not 2*sos; 
	//therefore, the sorting subtracts sos for es and is components!!!!

	//UO() << "Re-sort size=" << rs << "\n";

	resortsize_add = 0;

	IVector resortconstraint;
	resortconstraint.SetLen(ne);
	for (int i = 1; i <= ne; i++)
		resortconstraint(i) = 0;

	resortvector.SetLen(0);

	int actdofind = 1;
	int actnode = 1;

	for (int i=1; i<=ne; i++) 
	{
		Element* e = elements(i);

		//UO() << "El" << i << ", SOS2=" << e->SOSowned() << ", SOS2_RS=" << e->SOSowned_RS() << "\n";

		//DOF-references
		for (int k=actnode; k<=nodes.Length(); k++) 
		{
			actnode = k;
			int nsos = GetNode(k).SOS();
			for (int j=1; j <= nsos; j++)
			{
				resortvector.Add(GetNode(k).Get(j));
			}
			if ((nsos && GetNode(k).Get(nsos) > actdofind) && !(i == ne)) break;
		}
		actnode++;


		if (ResortMode2() && (e->SOSowned() > e->SOSowned_RS()))
		{
			//fill in part which goes to sos2 part	
			for (int j=1; j <= e->SOSowned_RS(); j++)
			{
				resortvector.Add(e->LTG(j));
			}
			//other part is filled together with constraints
		} 
		else //old mode, can only put constraints into sos-part
		{
			for (int j=1; j <= e->SOSowned(); j++)
			{
				resortvector.Add(e->LTG(j));
			}
			for (int j=e->SOSowned()+1; j <= e->SOSowned_RS(); j++)
			{
				resortvector.Add(e->LTG(j+e->SOSowned()+e->ES())-sos);
			}
		}
		if (resortvector.Length() != 0) actdofind = resortvector.Last();

		if (!TransformJacApply()) //sort constraints into multibody system -> better band structure
		{
			for (int k=1; k <= Minimum(e->NC(),2); k++) //cccc
			{
				//UO() << "constraint " << k << "\n";
				//$ YV 2012-12-13: we can be sure that Constraint is derived from Element
				Element * c = (Element*)(e->GetConstraint(k));
				int cind = c->GetOwnNum();

				int dist = 0;
				const int maxdist = 2; //cccc
				if (c->NE() == 2) 
				{
					dist = abs(c->GetElnum(1)-c->GetElnum(2));
				}
				else if (c->NE() == 1)
				{
					dist = 1;
				}
				else 
				{
					dist = maxdist+1; //do not resort
				}

				//UO() << "c-NE=" << c->NE() << "\n";

				if (c->ES() > 0 || c->SOS() > 0) 
				{
					dist = maxdist+1; // do not resort if special actuators, etc.
				}
				//UO() << "constraint" << cind << ": dist=" << dist << "\n";

				if (!resortconstraint(cind) && dist <= maxdist)
				{
					resortconstraint(cind) = 1;
					for (int j=1; j <= c->ES(); j++)
					{
						resortvector.Add(c->LTG(c->SOSowned_RS()*2+j)-sos);
						resortsize_add++;
					}

					for (int j=1; j <= c->IS_RS(); j++)
					{
						resortvector.Add(c->LTG(c->SOSowned_RS()*2+c->ES()+j)-sos);
						resortsize_add++;
					}
				}
			}
		}
	}
	for (int i=1; i<=ne; i++) 
	{
		Element* e = elements(i);
		if (!resortconstraint(i))
			for (int j=1; j <= e->ES(); j++)
			{
				resortvector.Add(e->LTG(e->SOSowned_RS()*2+j)-sos);
			}
	}

	/* //not necessary, everything is done with IS_RS()
	if (ResortMode2() && 0)
	{
	//add part which goes to constraint part:
	for (int i=1; i<=ne; i++) 
	{
	Element* e = elements(i);
	if (!resortconstraint(i) && (e->SOSowned() > e->SOSowned_RS()))
	{
	for (int j=e->SOSowned_RS()+1; j <= e->SOSowned(); j++)
	{
	resortvector.Add(e->LTG(j));
	}
	}
	}
	}*/

	for (int i=1; i<=ne; i++) 
	{
		Element* e = elements(i);
		if (!resortconstraint(i))
			for (int j=1; j <= e->IS_RS(); j++)
			{
				resortvector.Add(e->LTG(e->SOSowned()+e->SOSowned_RS()+e->ES()+j)-sos); //changed e->SOSowned_RS()*2 to e->SOSowned_RS()+e->SOSowned() .. not tested
			}
	}

	if (!ResortActive())
	{
		//UO() << "Resort deactivated!!!\n";
		for (int i=1; i <= sos+es+is; i++)
			resortvector(i) = i;
	}

	/*
	UO() << "Re-sort size new=" << GetResortSize() << "\n";

	UO() << "resortvector=[";
	for (int i=1; i <= resortvector.Length(); i++)
		UO() << resortvector(i) << ",";
	UO() << "]\n";

	UO() << "resortconstraint=[";
	for (int i=1; i <= resortconstraint.Length(); i++)
	UO() << resortconstraint(i) << ",";
	UO() << "]\n";
	
	//UO().InstantMessageText("Resortvector");
*/

}

void MultiBodySystem::ComputeMaxSparseBandwidth()
{
	elementbandwidth = 4; //minimum
	for (int i=1; i<=elements.Length(); i++) 
	{
		Element* e = elements(i);
		elementbandwidth = Maximum(elementbandwidth,e->ElementBandwidth());
	} 
	//UO() << "ComputeMaxSparseBandwidth=" << elementbandwidth << "\n";
}


int MultiBodySystem::TransformJacApply() const
{
	return transformJacApply;
}

void MultiBodySystem::ApplyTransform(const Vector& v, Vector& Av, int mode) //compute Av=A^T*v in mode==0 and Av=A*v in mode==1
{
	//should only be called by NumNLSys if transformJacApply is set
	Av = v;
	for (int i=1; i<=elements.Length(); i++)
	{
		if (elements(i)->TransformJacApply())
			elements(i)->ApplyTransform(v, Av, mode);
	}
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int MultiBodySystem::AddElement(Element* e)
{
	Element* ec = e->GetCopy();
	//ec->GetCol()=ec->GetMaterial().GetMaterialColor(); 			// element gets color from material:
	elements.Add(ec);
	elements.Last()->SetOwnNum(elements.Length());

	return elements.Length();
}

//$ DR 2013-05-23 removed deprecated function AddBody
//int MultiBodySystem::AddBody(Element* e)
//{
//	//search bodies:
//	int found = 0;
//	int i = 1; 
//	while (i <= NE() && !found)
//	{
//		if (!GetElement(i).IsType(TBody)) found = i;
//		i++;
//	}
//
//	if (found)
//	{
//		InsertElement(found, e);
//		return found; //must return element number!
//	}
//	else
//	{
//		Element* ec = e->GetCopy();
//		elements.Add(ec);
//		elements.Last()->SetOwnNum(elements.Length());
//		return elements.Length();
//	}
//
//}

//$ DR 2013-05-23 removed deprecated function InsertElement
//void MultiBodySystem::InsertElement(int i, Element* e)
//{
//	//resort element list:
//	elements.SetLen(NE()+1);
//
//	for (int j=NE(); j > i; j--)
//	{
//		elements(j) = elements(j-1);
//	}
//
//	Element* ec = e->GetCopy();
//	elements(i) = ec;
//
//	//set new element numbers:
//	for (int j=i; j <= NE(); j++)
//	{
//		elements(j)->SetOwnNum(j);
//	}
//
//	//update element numbers of constraints:
//	for (int j=1; j <= NE(); j++)
//	{
//		if (elements(j)->IsType(TConstraint))
//		{
//			//$ YV 2012-12-13: Constraint is just an Element here
//			Element* c = elements(j);
//			for (int k=1; k <= c->NE(); k++)
//			{
//				if (c->GetElnum(k) >= i)
//				{
//					c->SetElnum(k, c->GetElnum(k)+1);
//				}
//			}
//		}
//	}
//	for (int j=NSensors(); j >= 1; j--)
//	{
//		for (int k=1; k <= sensors(j)->GetNumberOfRelatedElements(); k++)
//		{
//			if (sensors(j)->GetRelatedElementNumber(k) >= i)
//			{
//				mystr tname = sensors(j)->GetTypeName();
//				sensors(j)->GetRelatedElementNumber(k)++;
//				//set new sensor names with new element numbers ...
//				if (sensors(j)->GetSensorName() == tname)
//				{
//					sensors(j)->SetSensorName(sensors(j)->GetTypeName());
//				}
//			} 
//		}
//	}
//}


int MultiBodySystem::AddAuxElement(Element* e)
{
	Element* ec = e->GetCopy();
	auxelements.Add(ec);
	//UO() << "add aux-elem=" << auxelements.Length() << "\n";

	return auxelements.Length();
}

int MultiBodySystem::AddNode(Node* n)
{
	Node* nc = n->GetCopy();
	nc->SetMBS(this);
	int j =  nodes.Add(nc);
	nc->NodeNum() = j;
	return j;
}

int MultiBodySystem::AddNode(Node* n, SearchTree& tree)
{
	//only works, if searchtree contains all existing nodes!

	double tol = 1e-8*CharacteristicLength();

	IVector items;

	Box3D box(n->Pos(),n->Pos());
	box.Increase(tol);

	tree.GetItemsInBox(box,items);
	for (int i = 1; i <= items.Length(); i++)
	{
		if (Dist(n->Pos(),nodes(items(i))->Pos()) <= tol) return items(i);
	}
	Node* nc = n->GetCopy();
	nc->SetMBS(this);

	int index = nodes.Add(nc);
	nc->NodeNum() = index;
	tree.AddItem(Box3D(n->Pos(),n->Pos()), index);

	return index;
}

int MultiBodySystem::AddBodyNode(Node* n)
{
	for (int i = 1; i <=nodes.Length(); i++)
	{
		if (Dist(n->Pos(),nodes(i)->Pos()) <= 1e-8*CharacteristicLength() && n->GetBodyInd() == nodes(i)->GetBodyInd())
		{
			return i;
		}
	}
	Node* nc = n->GetCopy();
	nc->SetMBS(this);
	int j =  nodes.Add(nc);
	nc->NodeNum() = j;
	return j;
}

int MultiBodySystem::AddBodyNode(Node* n, SearchTree& tree)
{
//$ AD 2011-02: AddBodyNode with Searchtree 

	double tol = 1e-8*CharacteristicLength();

	IVector items;

	Box3D box(n->Pos(),n->Pos());
	box.Increase(tol);

	tree.GetItemsInBox(box,items);
	for (int i = 1; i <= items.Length(); i++)
	{
		if ( (Dist(n->Pos(),nodes(items(i))->Pos()) <= tol) && (n->GetBodyInd() == nodes(items(i))->GetBodyInd()) ) return items(i);
	}
	Node* nc = n->GetCopy();
	nc->SetMBS(this);

	int index = nodes.Add(nc);
	nc->NodeNum() = index;
	tree.AddItem(Box3D(n->Pos(),n->Pos()), index);

	return index;
}



int MultiBodySystem::AddSensor(Sensor* s)
{
	Sensor* sc = s->GetCopy();
	sensors.Add(sc);
	sensors.Last()->SetOwnNum(sensors.Length());

	return sensors.Length();
}

//$ DR 2012-10: loads moved from element to mbs
int MultiBodySystem::AddLoad(const MBSLoad& li)
{
	//$ YV 2013-01-03: here we use GetCopy instead
	//MBSLoad* l = new MBSLoad(li);
	MBSLoad* l = li.GetCopy();
	l->SetOwnNum(loads.Length()+1);
	return loads.Add(l);

	// the following code is still in element.h/cpp:

	//int elem = li.HasElementDependence();
	//if (elem)
	//{
	//	TArray<int> elnums;

	//	elnums.Add(elem);
	//	GetElement(elem).GetDirectFeedThroughElements(elnums);

	//	for (int i=1; i<=elnums.Length(); i++)
	//	{
	//		elements.Add(elnums(i));
	//	}
	//}
}

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

void MultiBodySystem::AddRigidBodyDOFSensor(int bodynum, const Vector3D& pos_rel, mystr general_str)
{
	/*
	MBSSensor sensAPosX(this, (TMBSSensor)(TSElement+TSPos+TSX), bodynum, pos_rel);
	sensAPosX.SetSensorName(general_str + mystr("_pos_x"));
	AddSensor(&sensAPosX);
	MBSSensor sensAPosY(this, (TMBSSensor)(TSElement+TSPos+TSY), bodynum, pos_rel);
	sensAPosY.SetSensorName(general_str + mystr("_pos_y"));
	AddSensor(&sensAPosY);
	MBSSensor sensAPosZ(this, (TMBSSensor)(TSElement+TSPos+TSZ), bodynum, pos_rel);
	sensAPosZ.SetSensorName(general_str + mystr("_pos_z"));
	AddSensor(&sensAPosZ);

	MBSSensor sensAngX(this, (TMBSSensor)(TSElement+TSAngle), bodynum, pos_rel, Vector3D(1.,0.,0.), Vector3D(0.,0.,1.), Vector3D(0.,0.,1.));
	sensAngX.SetSensorName(general_str + mystr("_glob_ang_x"));			
	AddSensor(&sensAngX);
	MBSSensor sensAngY(this, (TMBSSensor)(TSElement+TSAngle), bodynum, pos_rel, Vector3D(0.,1.,0.), Vector3D(0.,0.,1.), Vector3D(0.,0.,1.));
	sensAngY.SetSensorName(general_str + mystr("_glob_ang_y"));			
	AddSensor(&sensAngY);
	MBSSensor sensAngZ(this, (TMBSSensor)(TSElement+TSAngle), bodynum, pos_rel, Vector3D(0.,0.,1.), Vector3D(1.,0.,0.), Vector3D(1.,0.,0.));
	sensAngZ.SetSensorName(general_str + mystr("_glob_ang_z"));			
	AddSensor(&sensAngZ);
	*/
}

// old version (AD)
//void MultiBodySystem::AddRigid2DBodyDOFSensor(int bodynum, mystr general_str, int usepos, int usevel, int displacement) // old version
//{
//	AddRigid2DBodyDOFSensor(bodynum,Vector2D(0.0,0.0),general_str,usepos,usevel,displacement);
//}

// new version (AD), please change code accordingly (insert Vector2D& pos_rel)
void MultiBodySystem::AddRigid2DBodyDOFSensor(int bodynum, Vector2D& pos_rel, mystr general_str, int usepos, int usevel, int displacement)
{
	/*
	if(usepos)
	{
		MBSSensor sensAPosX(this, TMBSSensor(TSElement + TSPos + TSX + TSplanar), bodynum, pos_rel); 
		sensAPosX.SetSensorName(general_str + mystr("_pos_x"));
		if (displacement)
		{
			sensAPosX.SetOffset(-GetElement(bodynum).GetXInit().Get(1)); // quick solution, replace with GetPosInit()...
			sensAPosX.SetSensorName(general_str + mystr("_disp_x"));
		}
		AddSensor(&sensAPosX);
	
		MBSSensor sensAPosY(this, TMBSSensor(TSElement + TSPos + TSY + TSplanar), bodynum, pos_rel); 
		sensAPosY.SetSensorName(general_str + mystr("_pos_y"));
		if (displacement) 
		{
			sensAPosY.SetOffset(-GetElement(bodynum).GetXInit().Get(2)); // quick solution, replace with GetPosInit()...
			sensAPosY.SetSensorName(general_str + mystr("_disp_y"));
		}
		AddSensor(&sensAPosY);

		MBSSensor sensAngX(this, TMBSSensor(TSElement + TSAngle + TSplanar), bodynum, pos_rel); 
		sensAngX.SetSensorName(general_str + mystr("_ang_phi"));
		AddSensor(&sensAngX);
	}
	if(usevel)
	{
		MBSSensor sensAVelX(this, TMBSSensor(TSElement + TSVel + TSX + TSplanar), bodynum, pos_rel); 
		sensAVelX.SetSensorName(general_str + mystr("_vel_x"));
		AddSensor(&sensAVelX);

		MBSSensor sensAVelY(this, TMBSSensor(TSElement + TSVel + TSY + TSplanar), bodynum, pos_rel); 
		sensAVelY.SetSensorName(general_str + mystr("_vel_y"));
		AddSensor(&sensAVelY);

		MBSSensor sensAngVelX(this, TMBSSensor(TSElement + TSVel + TSAngle + TSplanar), bodynum, pos_rel); 
		sensAngVelX.SetSensorName(general_str + mystr("_ang_vel_phi"));
		AddSensor(&sensAngVelX);
	}
	*/
}

int MultiBodySystem::AddMaterial(const Material& m)
{
	Material* mc = m.GetCopy();

	materials.Add(mc);
	materials.Last()->SetMaterialNumber(materials.Length());

	return materials.Length();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MultiBodySystem::DeleteElement(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if element is used by a Geomelement, not necessary, just set to 0
	// check if element is used by a sensor
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedElements(); k++)
		{
			if(GetSensor(j).GetRelatedElementNumber(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete element " << i << ", because it is still used by sensor " << j << "!\n";
				return;
			}
		}
	}
	// check if element is used by a load, not necessary, will be adjusted to non-dependent load

	// check if element is used by an other element
	for (int j = 1; j <= NE(); j++)
	{
		if(j != i) // maybe some IO-Elements have themselves as inputs somewhen, then they should also be deleted
		{
			for (int k = 1; k<=GetElement(j).NE(); k++)
			{
				if(GetElement(j).GetElnum(k) == i)
				{
					UO(UO_LVL_err) << "ERROR: Can not delete element " << i << ", because it is still used by element " << j << "!\n";
					return;
				}
			}
		}
	}

	///////////----------- DELETE ---------------------------//////////////////
	// delete the element in the MBS
	delete elements(i);	// just deletes the element, does not change the length of the TArray
	elements.Erase(i);		// decreases the length of the TArray
	UO(UO_LVL_warn) << "Element "<< i << "deleted\n";

	///////////----------- ADJUST ---------------------------//////////////////
	//set new element numbers:
	for (int j=i; j <= NE(); j++)
	{
		elements(j)->SetOwnNum(j);
	}
	
	// adjust for GeomElements
	for (int j = 1; j <= NDrawElements(); j++)
	{
		int num = GetDrawElement(j)->GetElnum();
		if (num == i)
		{
			GetDrawElement(j)->SetElnum(0);
			UO(UO_LVL_warn) << "WARNING: Deleted element " << i << " had been used by GeomElement " << j << ". Check the properties of the GeomElement!\n";
		}
		else if (num > i)
		{
			GetDrawElement(j)->SetElnum(num-1);
		}
	}

	// adjust for sensors
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedElements(); k++)
		{
			if(GetSensor(j).GetRelatedElementNumber(k) > i)
			{
				GetSensor(j).GetRelatedElementNumber(k)--;
			}
		}
	}

	// adjust for loads
	for (int j = 1; j <= NLoads(); j++)
	{
		if(GetLoad(j).HasElementDependence())
		{
			int num, output;
			GetLoad(j).GetIOElementNr(num, output);
			if (num == i)
			{
				GetLoad(j).SetLoadFunc(0);	// no dependence on element i anymore
				UO(UO_LVL_warn) << "WARNING: Deleted element " << i << " had been used by load " << j << ". Check the properties of the load!\n";
			}
			else if (num > i)
			{
				GetLoad(j).SetIOElement(num-1, output);	// adjust number
			}
		}
	}

	// adjust for elements
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = 1; k<=GetElement(j).NE(); k++)
		{
			int num = GetElement(j).GetElnum(k);
			if (num > i)
			{
				GetElement(j).SetElnum(k,num-1);
			}
		}
	}
}

void MultiBodySystem::DeleteSensor(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if sensor is used by an element
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = 1; k<=GetElement(j).NSensors(); k++)
		{
			if(GetElement(j).GetSensorNum(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete sensor " << i << ", because it is still used by element " << j << "!\n If this is a control element, just remove the sensor from the IOBlock inputs, than you can delete the sensor.\n";
				return;
			}
		}
	}
	// check if sensor is used by another sensor
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedSensors(); k++)
		{
			if(GetSensor(j).GetRelatedSensorNumber(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete sensor " << i << ", because it is still used by sensor " << j << "!\n";
				return;
			}
		}
	}
	///////////----------- DELETE ---------------------------//////////////////
	// delete the sensor in the MBS
	delete sensors(i);	// just deletes the sensor, does not change the length of the TArray
	sensors.Erase(i);		// decreases the length of the TArray
	UO(UO_LVL_warn) << "Sensor "<< i << "deleted\n";

	///////////----------- ADJUST ---------------------------//////////////////
	//set new sensor numbers:
	for (int j=i; j <= NSensors(); j++)
	{
		sensors(j)->SetOwnNum(j);
	}

	// adjust for sensors
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedSensors(); k++)
		{
			if(GetSensor(j).GetRelatedSensorNumber(k) > i)
			{
				GetSensor(j).GetRelatedSensorNumber(k)--;
			}
		}
	}

	// adjust for elements
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = 1; k<=GetElement(j).NSensors(); k++)
		{
			int num = GetElement(j).GetSensorNum(k);
			if (num > i)
			{
				GetElement(j).SetSensorNum(k,num-1);
			}
		}
	}


}

//$ DR 2012-10: loads moved from element to mbs
void MultiBodySystem::DeleteLoad(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if load is used by a sensor
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedLoads(); k++)
		{
			if(GetSensor(j).GetRelatedLoadNumber(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete load " << i << ", because it is still used by sensor " << j << "!\n";
				return;
			}
		}
	}

	// not necessary to check for elements. if element uses load, than load is removed from this element!

	///////////----------- DELETE ---------------------------//////////////////
	// delete the load in the MBS
	delete loads(i);	// just deletes the MBSLoad, does not change the length of the TArray
	loads.Erase(i);		// decreases the length of the TArray
	UO(UO_LVL_warn) << "Load "<< i << "deleted\n";

	// correct the entries in the elements
	Element* elP;
	TArray<int> el_loads;
	int elHadLoad;

	for (int j=1; j < NE(); j++)
	{
		elP=GetElementPtr(j);
		el_loads.Flush();
		elHadLoad = 0;
		el_loads = elP->GetLoadNrs();
		for(int l=1; l<=el_loads.Length(); l++)
		{
			// find and delete the load in the element
			if(el_loads(l)==i)
			{
				el_loads.Erase(l);
				elHadLoad = 1;
			}
		}

		///////////----------- ADJUST ---------------------------//////////////////
		// correct the numbers of the loads in the load list of the element
		for(int l=1; l<=el_loads.Length(); l++)
		{
			if(el_loads(l)>=i)
			{
				el_loads(l)=el_loads(l)-1;
			}
		}

		elP->SetLoadNrs(el_loads);
	}
	
	//set new load numbers:
	for (int j=i; j <= NLoads(); j++)
	{
		loads(j)->SetOwnNum(j);
	}

	// adjust for sensors
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedLoads(); k++)
		{
			if(GetSensor(j).GetRelatedLoadNumber(k) > i)
			{
				GetSensor(j).GetRelatedLoadNumber(k)--;
			}
		}
	}
}

void MultiBodySystem::DeleteDrawElement(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if geom-element is used by an element --> not necessary, just set to 0

	///////////----------- DELETE ---------------------------//////////////////
	delete drawelements(i);			// just deletes the Element, does not change the length of the TArray
	drawelements.Erase(i);
	UO(UO_LVL_warn) << "GeomElement "<< i << "deleted\n";
	
	///////////----------- ADJUST ---------------------------//////////////////
	// only elements use GeomElements
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = GetElement(j).NGeomElements(); k >= 1; k--)
		{
			if (GetElement(j).GetGeomElementNum(k) == i)
			{
				GetElement(j).GetGeomElementNum(k) = 0;
				GetElement(j).SetAltShape(0);
			}
			else if (GetElement(j).GetGeomElementNum(k) > i)
			{
				GetElement(j).GetGeomElementNum(k)--;
			}
		}
	}
	// adjust of own number not necessary, GeomElements do not know their number!
}

int MultiBodySystem::DeleteMaterial(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if material is used by an element
	for (int j = 1; j <= NE(); j++)
	{
		if(GetElement(j).GetMaterialNum() == i)
		{
			UO(UO_LVL_err) << "ERROR: Can not delete material " << i << ", because it is still used by element " << j << "!\n";
			return 0;
		}
	}

	///////////----------- DELETE ---------------------------//////////////////
	delete materials(i);			// just deletes the material, does not change the length of the TArray
	materials.Erase(i);
	UO(UO_LVL_warn) << "Material "<< i << "deleted\n";

	///////////----------- ADJUST ---------------------------//////////////////
	for (int j=i; j <= NMaterials(); j++)
	{
		materials(j)->SetMaterialNumber(j);
	}

	// adjust for elements
	for (int j = 1; j <= NE(); j++)
	{
		if (GetElement(j).GetMaterialNum() > i)
		{
			GetElement(j).GetMaterialNum()--;
		}
	}
	return 1;
}

void MultiBodySystem::DeleteNode(int i)
{
	///////////----------- CHECK ---------------------------//////////////////
	// check if node is used by an element
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = 1; k<=GetElement(j).NNodes(); k++)
		{
			if(GetElement(j).NodeNum(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete node " << i << ", because it is still used by element " << j << "!\n";
				return;
			}
		}
	}

	// check if node is used by a sensor
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedNodes(); k++)
		{
			if(GetSensor(j).GetRelatedNodeNumber(k) == i)
			{
				UO(UO_LVL_err) << "ERROR: Can not delete node " << i << ", because it is still used by sensor " << j << "!\n";
				return;
			}
		}
	}

	///////////----------- DELETE ---------------------------//////////////////
	delete nodes(i);			// just deletes the Element, does not change the length of the TArray
	nodes.Erase(i);
	UO(UO_LVL_warn) << "Node "<< i << "deleted\n";

	///////////----------- ADJUST ---------------------------//////////////////
	for (int j=i; j <= NNodes(); j++)
	{
		nodes(j)->NodeNum() = j;
	}

	// adjust for elements
	for (int j = 1; j <= NE(); j++)
	{
		for (int k = 1; k <= GetElement(j).NNodes(); k++)
		{
			if (GetElement(j).NodeNum(k) > i)
			{
				GetElement(j).NodeNum(k)--;
			}
		}
	}

	// adjust for sensors
	for (int j = 1; j <= NSensors(); j++)
	{
		for (int k = 1; k<=GetSensor(j).GetNumberOfRelatedNodes(); k++)
		{
			if(GetSensor(j).GetRelatedNodeNumber(k) > i)
			{
				GetSensor(j).GetRelatedNodeNumber(k)--;
			}
		}
	}

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// get the i-th sensors (as registered in mbs) own stored values from internal array
TArray<double>* MultiBodySystem::GetSensorValuesArrayPtr(int i)
{
	Sensor& the_sensor = GetSensor(i);
	TArray<double>* pValues = the_sensor.GetSignalHistoryValuesPtr(); 
	return pValues;
}

// get the i-th sensors (as registered in mbs) own stored times from internal array
TArray<double>* MultiBodySystem::GetSensorTimesArrayPtr(int i)
{
	Sensor& the_sensor = GetSensor(i);
	TArray<double>* pTimes = the_sensor.GetSignalHistoryTimesPtr(); 
	return pTimes;
}

// get the i-th sensors (as registered in mbs) name
mystr MultiBodySystem::GetSensorName(int i)
{ 
	return sensors(i)->GetSensorName(); 
}

// get the number of the (first) sensor that matches the sensorname
int MultiBodySystem::GetSensorNumber(mystr& name)                                                                   
{
	for(int i=1; i<=NSensors(); i++)
	{
		if(name.Compare(sensors(i)->GetSensorName()))
			return i;
	}
	return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




void MultiBodySystem::SetGlobalInitConditions(Vector& x0)
{
	for (int i=1; i<=elements.Length(); i++)
	{
		elements(i)->SetGlobalInitConditions(x0);
	}
}

void MultiBodySystem::SetGlobalInitData(Vector& xdata)
{
	for (int i=1; i<=elements.Length(); i++)
	{
		elements(i)->SetGlobalInitData(xdata);
	}
	//$ MaSch 2013-10-28: added support of SPH particle data via auxiliary elements
	for (int i=1; i<=auxelements.Length(); i++)
	{
		if(auxelements(i)->IsType(TParticle))
		{
			auxelements(i)->SetGlobalInitData(xdata);
		}
	}
}

int MultiBodySystem::NConstraints() const
{
	int nc = 0;
	for (int i=1; i <= NE(); i++) if (elements(i)->IsType(TConstraint)) nc++;
	return nc;
}

Box3D MultiBodySystem::GetBoundingBoxD() const
{
	Box3D b;
	for (int i=1; i<=elements.Length(); i++)
	{
		b.Add(elements(i)->GetBoundingBoxD());
	}

	//global_uo << "box1=" << b.PMin() << "," << b.PMax() << "\n";
	//bounding box of drawelements that belong to an element is managed in element function!!!
	for (int i=1; i<=drawelements.Length(); i++)
	{
		if (!drawelements(i)->GetElnum())
			b.Add(drawelements(i)->GetBoundingBoxD());
	}
	//global_uo << "box2=" << b.PMin() << "," << b.PMax() << "\n";

// AD 2013-01-18:  adapted for new "grid" == scene background
	if (GetIOption(126)) //draw grid:
	{
		double refpos_x = GetDOption(107);
		double refpos_y = GetDOption(108);
		double refpos_z = GetDOption(109);
    double size_of_plane_x = GetDOption(112);
    double size_of_plane_y = GetDOption(113);
    double size_of_plane_z = GetDOption(142); 

		Vector3D p1(refpos_x, refpos_y, refpos_z);
		Vector3D sz(size_of_plane_x, size_of_plane_y, size_of_plane_z);
		Box3D bg(p1, p1+sz); 
		b.Add(bg);
	}

	b.Increase(0.001); //original: for better graphics; $JG2013-06-24
	b.InflateFactor(1.5); //increase bounding box, such that it fits better into the sceen $JG2013-06-24

	return b;
}



void MultiBodySystem::CollectAvailableFieldVariables()
{
	availableFieldVariables.Flush();

	for(int i = 1; i <= elements.Length() + auxelements.Length(); i++)
	{
		Element * element;
		if(i <= elements.Length())
			element = elements(i);
		else
			element = auxelements(i - elements.Length());
		TArray<FieldVariableDescriptor> elementVariables;
		element->GetAvailableFieldVariables(elementVariables);
		FieldVariableDescriptor::MergeArrays(availableFieldVariables, elementVariables);
	}

	// and now we update the index of the actually selected variable for post-processing
	// based on its textual identifier stored in the cfg file
	bool found = false;
	for(int i = 1; i <= availableFieldVariables.Length(); i++)
		if(strcmp(GetTOption(107), availableFieldVariables(i).GetTextualIdentifier()) == 0)
		{
			SetIndexOfActualPostProcessingFieldVariable(i);
			found = true;
			break;
		}
	if(!found)
		SetIndexOfActualPostProcessingFieldVariable(0);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Procedure for registration of model generation functions in the MBS-System
//$ YV 2013-01-02: the old mechanism is now replaced by the new one: the models are now created in the client dll;
// instead of the functions, which are commented out below, we simply access the models library from the dll
#if 0

const int modelFunctionAutoRegistration_max_n_fct = 256; //number of functions
int modelFunctionAutoRegistration_n_fct = 0; //number of functions
char* modelFunctionAutoRegistration_fname[modelFunctionAutoRegistration_max_n_fct];
char* modelFunctionAutoRegistration_description[modelFunctionAutoRegistration_max_n_fct];
int modelFunctionAutoRegistration_option[modelFunctionAutoRegistration_max_n_fct];
typedef int (*ModelFunctionAutoRegistration_pt2FunctionList)(MBS* mbs);//[modelFunctionAutoRegistration_max_n_fct];
ModelFunctionAutoRegistration_pt2FunctionList modelFunctionAutoRegistration_pt2FunctionList[modelFunctionAutoRegistration_max_n_fct];
ModelFunctionAutoRegistration_pt2FunctionList modelFunctionAutoRegistration_pt2Fcn_InitModelData[modelFunctionAutoRegistration_max_n_fct];

int MultiBodySystem::GetNumberOfModelFunctions() const 					     //return number of available model functions
{
	return modelFunctionAutoRegistration_n_fct;
}

char* MultiBodySystem::GetModelFunctionName(int index0) const        //return model name string
{
	return modelFunctionAutoRegistration_fname[index0];
}

char* MultiBodySystem::GetModelFunctionDescription(int index0) const //return model description string
{
	return modelFunctionAutoRegistration_description[index0];
}

int MultiBodySystem::GetModelFunctionOption(int index0) const				 //return option of model function
{
	return modelFunctionAutoRegistration_option[index0];
}

int MultiBodySystem::CallModelFunction(int index0)				     //Call model function with index "index0"
{
	return (modelFunctionAutoRegistration_pt2FunctionList[index0])(this);
}

int MultiBodySystem::CallModelFunction_InitModelData(int index0)				     //Call model function with index "index0"
{
	return (modelFunctionAutoRegistration_pt2Fcn_InitModelData[index0])(this);
}

int MultiBodySystem::ModelFunctionHasInitModelData(int index0)				     //check if model function has initialization procedure
{
	if (modelFunctionAutoRegistration_pt2Fcn_InitModelData[index0] != 0)
	{
		return 1;
	}
	else 
	{
		return 0;
	}
}

void ModelFunctionAutoRegistration(int(*function_ptr)(MBS* mbs), const char* functionName, const char* description, int option,
																	 int(*function_ptr_init_modeldata)(MBS* mbs))
{
	if (modelFunctionAutoRegistration_n_fct < modelFunctionAutoRegistration_max_n_fct)
	{
		int index = modelFunctionAutoRegistration_n_fct;

		//modelFunctionAutoRegistration_pt2Function/*[index]*//* = function_ptr;
		modelFunctionAutoRegistration_pt2FunctionList[index] = function_ptr;

		//function for initilization of model data:
		modelFunctionAutoRegistration_pt2Fcn_InitModelData[index] = function_ptr_init_modeldata;

		int length = strlen(functionName);
		modelFunctionAutoRegistration_fname[index] = new char[length + 1];
		strcpy(modelFunctionAutoRegistration_fname[index], functionName);

		length = strlen(description);
		modelFunctionAutoRegistration_description[index] = new char[length + 1];
		strcpy(modelFunctionAutoRegistration_description[index], description);

		modelFunctionAutoRegistration_option[index] = option;

		//increase index after all initialization!
		modelFunctionAutoRegistration_n_fct++;
	}
}

#endif

void MultiBodySystem::ModelFunctionChanged()
{
	SetModelData_Initialized(0);
	Initialize();
}


int MultiBodySystem::GetModelFunctionsIndex0()
{
	int ind0 = -1;
	int nmodel = GetModelsLibrary()->GetModelsCount();

	HOTINTOptions hotint_options_nonconstant_copy(*GetOptions());
	mystr model_str = hotint_options_nonconstant_copy.GeneralOptions()->ModelFileInternalModelFunctionName();
	
	//$ YV 2013-01-04: model indices are 1-based
	for (int i=1; i <= nmodel; i++)
	{
		//UO() << "Model " << i << "= '" << GetModelFunctionName(i) << "', description='" << GetModelFunctionDescription(i) << "'\n";
		if (model_str == GetModelsLibrary()->GetModelInterface(i)->GetMBSModelName())
		{
			ind0 = i;
		}
	}
	return ind0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





double MultiBodySystem::EvaluateParsedFunction1D(const mystr & parsedFunctionExpression, const mystr & parsedFunctionVariables, double t)
{
	CParsedFunction pf;
	//$ YV 2012-13-12: changed MBSParser() to &MBSParser()
	pf.SetParsedFunction1D(&MBSParser(), parsedFunctionExpression, parsedFunctionVariables);
	return pf.Evaluate(t);
}

double MultiBodySystem::ExpressionToDouble(mystr & expression)
{
	int error;
	return MBSParser().ExpressionToDouble(expression, error);
}

//$ DR 2013-01-14: added as new interface to script parser
int MultiBodySystem::ExecuteParserCommand(mystr & command, TArray<ElementDataContainer*>& parameter_EDCs, ElementData& return_value, int option)
{
	return EDCParser().ExecuteParserCommand(command, parameter_EDCs, return_value, option);
}
