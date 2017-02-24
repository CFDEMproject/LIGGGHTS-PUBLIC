//#**************************************************************
//#
//# filename:             timeint.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						09.06.99
//# description:          Class for implicit and explicit Runge Kutta Time integration, 
//#												variable stepsize and arbitrary order
//# remarks:						  The file tableaus.txt must be provided in Project/Release or Debug !!!
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
 
//
//
// TODO:
// Add option: continue anyway, even if not converged
//
//
//
//


//#include "stdafx.h"
#include "ioincludes.h"

#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <math.h>

#include "mystring.h"  

#include "tarray.h"    //i1
#include "myfile.h"    //i1

//#include "femath.h"    //+i2
#include "timeint.h"
#include "options_auto.h"

//$ YV 2013-01-03: this cannot be done in this way any longer
//#include  "..\ElementsLib\contact_globalvariables.h"    //access these (extern global) variables for reasons of log-output
//double global_time; //global time

#ifdef gencnt
int generated_vectors;
int generated_matrices;
int * genvec = &generated_vectors;
int * genmat = &generated_matrices; 
#endif

UserOutputInterface * global_uo;

TimeInt::TimeInt():stage_tableaus(), tableaus()
{
	//AP 2013-01-15: set global_uo pointer already in constructor, 
	//      it is needed in CreateWCDObject() in WorkingModule.cpp, 
	//      line 96: pModelsLibrary->SetGlobalVariablesPointers(global_uo, &TMtspent, &TMtstart);
	global_uo = &UO();

	TIInit(); 
};

void TimeInt::InitFirst() 
{
	//AD: Initialize the paused flag to zero - otherwise pressing the pause button may result in crashes
	bPause = 0;  

	//drawing:
	transparency_on = 1;

	reduceimage = 1;
	TIdrawtexres = 4;
	texImage = NULL;
	texstorenn = 0;
	SetTexStoreResolution(3); //32

	InitializeOptions();

	SetJacobianComputationFlag(0);

	//AP 2013-01-15: global_uo pointer is set already in constructor TimeInt()
	//global_uo = &UO();

  // hook UserOutput object to solver settings
	//$ YV 2012-12-14: UO() and the actual one are the same, no need to hook
	/*UO().HookToSolverSettings(& GetOptions()->LoggingOptions()->OutputLevel(), & GetOptions()->LoggingOptions()->FileOutputLevel(), &logout, & GetOptions()->LoggingOptions()->CriticalLogFileSize(), & GetOptions()->LoggingOptions()->MaxErrorMessages(), & GetOptions()->LoggingOptions()->MaxWarningMessages(),
		& GetOptions()->LoggingOptions()->OutputPrecisionDouble(), & GetOptions()->LoggingOptions()->OutputPrecisionVector(), & GetOptions()->LoggingOptions()->OutputPrecisionMatrix()); //$ PG 2012-2-3: substituted solset.output_level by GetOptions()->LoggingOptions()->OutputLevel()*/

	//called from Working module base class after first initialization
	//do only once initialization:
	numsol = NumNLSolver(this, &solset); //solset: set a reference to solset for variables, which show immediate effect of modification during simulation (e.g. logging of newton iterations)
	//numsol.SetReferenceToSolverOptions(solset); 

	TIkv = 1;
	drawlines = 0;
	drawlinesh = 0;

	loadfact = 1;

	NLS_SetJacCol(0);

	SetTime(0);
	TIstep = 0;
	act_tableau = 1;

	//$ RL 2012-5-29:[ get correct path of tableaus.txt
	mystr tableaus_file("tableaus.txt");
	if(__argc>0)
	{
		mystr strarg(__argv[0]); // first argument is path of hotint.exe, where also tableaus.txt is located
		int until;
		for(until=strarg.Length()-1;until>=0;until--)
		{
			if(strarg.PosPeek(until) == '\\')
				break; // until <=> last backslash
		}
		mystr path("");
		for(int i=0;i<=until;i++)
		{
			path += strarg.PosPeek(i); // path until last backslash
		}
		tableaus_file = path + "tableaus.txt";
	}
	LoadTableaus(tableaus_file);
	//$ RL 2012-5-29:] get correct path of tableaus.txt

	TMResetTimer();
}

void TimeInt::TIInit()
{
	TIsteprecommendation = 1e100;

	TIfinished = 0;
	TIdrawtime = 0;

	TImincol = 1e30;
	TImaxcol = -1e30;
	flag_compute_minmax_FEcol = 0;

	TInoincs = 0;
	TIwritesolevery = 0.;

	ada_err_sum = 0;

	acterror = 0;
	acterror_damp = 0;
	TIm_initialized = 0;
	bandwidthm = 0;
	reducedbandsize = 0;

	drawnow = 1;
	TInumber_of_solution_data_steps = 0; //PG old: .SetLen(0);

	lastloadresults = -1e100;
	laststoredata = -1e100;

	tidraw_offx = 0.;
	tidraw_offy = 0.;
	tidraw_offz = 0.;

	TIissubstep = 0;
	oldenergy = 0;

	TIcomptime = 0;
	TItimeperstep = 0;
	colline = Vector3D(0.1,0.1,0.1);

	log.Init();
	logout_name = "hotint.log";
}

void TimeInt::SetStartVector(const Vector& x0) 
{
	TIx0 = x0; 
	TIx0draw = x0;

	TISaveState(TIlastnonlinit_state);

	if (TIkv) 
	{
		int sos = GetSecondOrderSize();
		int es = GetFirstOrderSize();
		int is = GetImplicitSize();

		TIk0.SetLen(2*sos+es+is);
		TIk0.SetAll(0);
		int m_invertable = 0;

		if (sos != 0 && m_invertable)
		{
			evalf.SetLen(sos);
			Matrix m(sos,sos);
			EvalF2(x0, evalf, TItime); //compute F-Ku-Gv
			EvalM(x0, m, TItime);
			Vector xh(sos);
			if (!m.Invert2()) {uo << "K-version, Init: Mass matrix not invertible!!!" << ENDL;}
			else {xh = m*evalf;}

			//second order start vector
			TIk0.Copy(TIx0,1+sos,1,sos);
			TIk0.Copy(xh,1,1+sos,sos);
		}
		//first order start vector
		if (es != 0)
		{
			evalf.SetLen(es);
			EvalF(x0, evalf, TItime); //compute First order equ.
			TIk0.Copy(evalf,1,1+2*sos,es);
		}

		//implicit start vector
		TIk0.Copy(TIx0,1+2*sos+es,1+2*sos+es,is);

	}
}

int TimeInt::NonlinSubStep()
{
	TMStartTimer(21);
	int rv = 1;
	int end = 0;
	double corerr;
	int it = 0;
	TISaveState(TIlastnonlinit_state);


	TIdiscontstep = 0;

	//global_time = TItime;

	int olddiscont;
	while (!end && rv)
	{
		olddiscont = TIdiscontstep;

		TIRestoreState(TIlastnonlinit_state);

		it++;
		log.TInonlinit = it;
		//Vector f(GetSystemSize());

		rv = GeneralIStep();
 		corerr = PostNewtonStep(TItime);


		if (corerr == -1) 
		{
			rv = 0; end = 1;
			reduce_order = 2;
		}
		else if (corerr < DiscontinuousAccuracy()) 
		{
			end = 1;
		}
		else if (TIdiscontstep != olddiscont) 
		{
			reduce_order = 2;
			numsol.DestroyOldJac(TIstep);
		}

		
		if (GetStepRecommendation() < TIstep && TIstep > TIminstep)
		{
			end = 2;
			rv = 0;
			//uo << "Step Recommendation:" << GetStepRecommendation() << "\n";
		}

		if (it >= MaxDiscontinuousIt()) 
		{
			end = 1;/* uo << "ERROR to many iterations in nonlinear-corrector iteration!!!" << ENDL;*/
			if(!solset.ignore_maxdiscontinuousit) rv = 0;
		}
	}

	if (end != 2) 
	{
		PostprocessingStep();
	}
	TInonls = it;

	TMStopTimer(21);

	return rv;
}

int flag_applyjac = 0;
int minstepwarned = 0;
//ofstream errout("..\\..\\output\\comperr.txt");
int TimeInt::FullStep()
{
	TISaveState(TIlaststep_state); //save state of end of last step (or initial step)
	SetStepRecommendation(1e100);

	int mode = 3; //1 is single step, 2 is adaptive substep if failure, 3 is adaptive stepreduction
	StartTimeStep();
	TIstep = TIstepnew;
	// AP: TItemp_state is used to save state after StartTimeStep(), for AND this means after transport
	TISaveState(TItemp_state);

	if (GetStepRecommendation() < TIstepnew)
	{
		TIstep = GetStepRecommendation();
	}

	if (mode == 1)
	{
		log.TInewtonit = 0;
		int rv = NonlinSubStep();
		FinishStep();

		return rv;
	}
	else if (mode == 3)
	{

		double savestep = TIstep;
		log.TInewtonit = 0;
		double TItimeold = TItime;
		int rv = 0;
		double it = 0;
		while (it < 15 && !rv)
		{	
			it++;
			if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: sub-iteration = ")+mystr(it)+mystr(", time increment = ")+mystr(TIstep)+mystr(", time = ")+mystr(TItime)+mystr("\n");
			SetStepRecommendation(1e100);
			//if (TIstep != savestep) numsol.DestroyOldJac(TIstep);
			if (TItime+TIstep > solset.endtime)
			{
				TIstep = solset.endtime - TItime;
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: set time increment = ")+mystr(TIstep)+mystr(", otherwise end time would be exceed.\n");
			}

			rv = NonlinSubStep(); 

			if (!rv && it < 15) 
			{

				TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables

				if (GetStepRecommendation() < TIstep) 
				{
					TIstep = GetStepRecommendation();
					if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: recommended time increment: ")+mystr(TIstep)+mystr("\n");
	
					if (it > 3)
					{
						TIstep *= 0.5;
						if (SolverPrintsDetailedOutput()) SolverUO() << "step: time increment multiplied by 0.5\n";
					}
					else if (it > 7)
					{
						TIstep *= 0.1;
						if (SolverPrintsDetailedOutput()) SolverUO() << "step: time increment multiplied by 0.1\n";
					}
				}
				else
				{
					TIstep *= 0.5;
					if (SolverPrintsDetailedOutput()) SolverUO() << "step: time increment multiplied by 0.5\n";
					
					if (TIstep <= TIminstep)
					{
						it = 14;
						if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: set sub-iteration = 14 to perform final sub-iteration (max-sub-iterations=15), since time increment ")+mystr(TIstep)+mystr(" is less equal TIminstep=")+mystr(TIminstep)+mystr("\n");
					}
				}
					//?only limit minimum step size in case that step is reduced due to bad convergence
				if (TIstep < TIminstep) 
				{

					TIstep = TIminstep; 
					if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: set time increment to TIminstep=")+mystr(TIminstep)+mystr(", since it's been already smaller than that.\n");
				}	
				//}
			}
		}
		//if (it > 1) uo << "stepsize = " << TIstep << "\n";

		if (!rv)
		{
			if (!minstepwarned) uo << "ERROR: Step at minimal stepsize, step not fully converged. Further warning suppressed!" << ENDL;
			minstepwarned = 1;
			//rv = 1;
		}

		FinishStep();

		TIstep = savestep;

		return rv;
	}
	else
	{
		double savestep = TIstep;
		log.TInewtonit = 0;
		double TItimeold = TItime;
		int rv = 0;

		double it = 1;
		double it2 = 1;
		while ((TIstep > TIminstep || TIstep == TIstepnew) && rv == 0)
		{
			SetStepRecommendation(1);

			TItime = TItimeold;
			TIstep = savestep/it;
			if (it > 1) numsol.DestroyOldJac(TIstep);
			if (TIstep > TIminstep && it < 256) numsol.MaxFullNewtonSteps() = 0;
			else numsol.MaxFullNewtonSteps() = NLS_MaxFullNewtonSteps();

			for (int i=1; i <= it2; i++)
			{
				rv += abs(NonlinSubStep()-1); 
				TItime += TIstep;
			}
			if (rv == 0) rv = 1;
			else rv = 0;

			if (!rv) 
			{
				TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables
				if (mode == 2) 
				{
					it *= 2.;
					it2 *= 2;
				}
				else
				{
					if (TIstepnew/GetStepRecommendation() > it) 
					{
						it = TIstepnew/GetStepRecommendation();
						//uo << "reducestep=" << TIstepnew/it << ", it=" << it << "\n";
					}
					else
						it *= 2.;

					if (GetStepRecommendation() < TIminstep) it = TIstepnew/TIminstep;
				}
			}
		}
		//if (it > 1) uo << "stepsize = " << TIstep << "\n";

		TIstep = savestep;
		TItime = TItimeold;	

		if (!rv) uo << "ERROR: Step at minimal stepsize, computation terminated" << ENDL;
		FinishStep();

		return rv;
	}


}

void TimeInt::TISaveState(TIstepsave& stepsave)
{
	stepsave.x0 = TIx0;
	stepsave.k0 = TIk0;
	stepsave.time = TItime;
	stepsave.data = TIdata;	
	/*
	TIx0save = TIx0;
	TIk0save = TIk0;
	TItimesave = TItime;
	*/
	SaveState(); //for multibody system, if internal states exist (should not happen)
}

void TimeInt::TIRestoreState(const TIstepsave& stepsave)
{
	TIx0 = stepsave.x0;
	TIk0 = stepsave.k0;
	TItime = stepsave.time;
	TIdata = stepsave.data;
	/*
	TIx0 = TIx0save;
	TIk0 = TIk0save;
	TItime = TItimesave;
	*/
	RestoreState(); //for multibody system, if internal states exist (should not happen)
}

int TimeInt::FullAdaptiveStep()
{
	int comp_ref=0; //no exact reference-value computation

	TIstep = TIstepnew;
	int success = 1;
	TInoincs=1;
	int end = 0;
	int firstrun = 1;
	log.TInewtonit = 0;
	TIactlocerr = AbsAccuracy();

	double value1;
	double value2;
	double value3;

	double acterror_damp_save = acterror_damp;
	oldenergy = GetEnergy();

	int rv; //0=success, 1=failure

	//parameters for adaptive strategy
	//double smin=tableaus[act_tableau]->ODE_order;
	//if (smin < 2) {smin = 2;}
	double smin=2;              //4
	double smax=smin*smin;      //smin^2
	double tsmaxinc = 5;        //5
	double tsmininc = 2;        //2
	double tsmindec = 2;        //2
	double err_damp_fact = 0.1; //0.1

	//safety parameter
	double sts=1;

	TISaveState(TIlaststep_state); //save state of end of last step (or initial step)

	TIrejectedsteps--;
	//try until successful timestep
	while (!end)
	{
		TIrejectedsteps++;
		//restore old state:
		TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables
		acterror_damp = acterror_damp_save;

		int largestep_discont=0;

		//compute error in this timestep:
		TIstepnew = TIstep; //suggestion for next timestep
		double prevstep = TIstep;
		comperr=0;

		double disttoraster = ((int(TItime/TImaxstep))*TImaxstep)+TImaxstep-TItime;
		if (fabs(disttoraster) <= 0.1*TIminstep)
		{disttoraster += TImaxstep;}

		if (TIstep >= (1-0.1*TIminstep)*disttoraster)
		{
			TIstep = disttoraster;
			prevstep = TIstep;
		}
		else if (TIstep > 0.8*disttoraster)
		{
			TIstep = 0.5*disttoraster;
			prevstep = TIstep;
		}

		TIissubstep = 0;
		rv = abs(NonlinSubStep()-1);  //0 = ok, 1 = failure

		if (rv==0) 
		{
			//1. Vergleichswert ok
			value1=GetError();
			if (TIdiscontstep) largestep_discont = 1;
			TIstep *= 0.5;
		}
		else
		{
			//not converged
			if (TIstep == TIminstep)
			{
				uo << "ERROR: Step at minimal stepsize, computation terminated" << ENDL;
				return 0;
			}

			log.changestep++;
			reduce_order += 1;

			TIstep *= 0.5;
			TIstep = Maximum(TIstep, TIminstep);
			rv=0;
			continue;
		}

		//restore old state:
		TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables


		TIissubstep = 1;

		rv = abs(NonlinSubStep()-1); TItime+=TIstep;
		rv += abs(NonlinSubStep()-1);

		if (rv!=0) 
		{
			log.changestep++;
			reduce_order += 1;

			if (prevstep == TIminstep)
			{
				uo << "ERROR: Step at minimal stepsize, computation terminated" << ENDL;
				return 0;
			}

			TIstep *= 0.5;
			TIstep = Maximum(TIstep, TIminstep);
			continue;

		}

		//store/compare error:
		value2=GetError();

		TIissubstep = 0;

		if (comp_ref)
		{

			//restore old state:
			TIstep = prevstep;
			TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables

			//Sehr genauen wert berechnen für Vergleichszwecke:
			double fineness = 10;
			TIstep /= fineness;
			rv = 0;
			for (int ii = 1; ii <= fineness ; ii++)
			{
				rv += abs(NonlinSubStep()-1); 
				TItime += TIstep;
			}
			value3=GetError();

			if (rv==1) 
			{
				uo << "ERROR: FullAdaptive Step-Failure at 0.1 base-step" << ENDL;
				uo << "Exit" << ENDL;
				return 0;
			}
			comperr=sts*fabs(value2-value3);
		}

		TItime = TIlaststep_state.time; //old starttime


		TIstep = prevstep;
		
		//compute Errors

		double D1=fabs(value1-value2);

		double rowsum = 1./(pow(2.,tableaus[act_tableau]->ODE_order)-1.);
		double ada_error=sts*D1*rowsum;

		acterror = ada_error;
		if (acterror < acterror_damp) {acterror_damp = acterror;}
		else {acterror_damp = (acterror_damp*(1.-err_damp_fact)+acterror*err_damp_fact);}

		ada_error = acterror_damp; 

		//uo << "val1=" << value1  << ",val2=" << value2 << ", step=" << TIstep << ", err=" << acterror << ", adaerr=" << ada_error << ", crit=" << AbsAccuracy() << "\n";


		//double CalcErr;
		double TIstepold = TIstep;

		//suggest new (larger) timestep-size
		if (ada_error <= AbsAccuracy()/smax && TIstep<TImaxstep)
		{
			end=1; //timestep completed

			if (ada_error == 0.) {ada_error = 1e-16;}
			if (tableaus[act_tableau]->nstages == 0.) 
			{(*global_uo) << "Internal Error: tableaus[act_tableau]->nstages is Zero!" << "\n";}

			//cout << "Aktuelle Zeitschrittweite: " << TIstep << endl;
			TIstepnew=TIstep*pow(AbsAccuracy()/(ada_error*smin),1./(double)tableaus[act_tableau]->ODE_order);
			log.changestep++;

			if ((1.9*TIstep < TIstepnew) && (TIstepnew < 3*TIstep))
			{
				TIstepnew = 2.*TIstep;
			}

			if (TIstepnew>TImaxstep)
			{
				TIstepnew=TImaxstep;
			}

			if (TIstepnew > tsmaxinc*TIstep) 
			{	
				TIstepnew = tsmaxinc*TIstep;
			}

			if (TIstepnew <= tsmininc*TIstep) 
			{
				TIstepnew = TIstep; 
				log.changestep--;
			}
		}
		else if(ada_error > AbsAccuracy())
		{
			if (TIstep==TIminstep)
			{
				uo << "ERROR: Step at minimal stepsize, computation terminated\n";
				uo << "**** Consider changing error tolerance and Newton accuracy ****\n";
				return 0;
			}
			//suggest smaller timestep (> minstepsize, otherwise success=0)
			log.changestep++;

			{
				TIstepold = TIstep;
				//cout << "Aktuelle Zeitschrittweite: " << TIstep << endl;
				double od = tableaus[act_tableau]->ODE_order;

				if (largestep_discont && od > 2) 
				{
					od = 2; 
				}

				TIstep=TIstepold*pow(AbsAccuracy()/(smin*ada_error),1./od);

				if (largestep_discont && TIstep < 0.25*TIstepold) 
				{
					TIstep = 0.25*TIstepold;
				}

				if (TIstep*tsmindec > TIstepold)
				{
					TIstep = TIstepold/tsmindec;
				}
			}
			if (TIstep < TIminstep) {TIstep = TIminstep;}

			end=0;
			TInoincs++;

		}
		else
		{
			end=1;

		}
	}

	//Needs average (min/max) variables of last 20 to 50 steps:
	//  newtonits, nonlinits, isdiscontinuous, ...
	//  also steps/Jacobianupdate --> stage selection lower if frequently ???
	steps_since_orderchange++;
	min_newtonits = Minimum(min_newtonits, GetNewtonIt());
	newjacobians_since_orderchange = log.jaccount - last_jaccount;

	if (solset.variable_order)
	{
		const int minsteps = 30;
		int old_stage = act_stage;
		if (newjacobians_since_orderchange == 0) newjacobians_since_orderchange = 1;

		double steps_per_jac = (double)steps_since_orderchange/(double)newjacobians_since_orderchange;

		if (act_stage > 1)
		{
			if (reduce_order || (steps_since_orderchange > minsteps && (steps_per_jac <= 2 || min_newtonits > 11)))
			{
				act_stage--;
				if (reduce_order > 1) act_stage = 1;
				else TIstepnew *= 0.5; 
				act_tableau = stage_tableaus(act_stage);
			}
		}
		if (act_stage < stage_tableaus.Length())
		{
			if (!reduce_order && steps_per_jac >= 4 &&
				min_newtonits <= 6 && (steps_since_orderchange > minsteps))
			{
				act_stage++;
				int old_tab = act_tableau;
				act_tableau = stage_tableaus(act_stage);

				TIstepnew *= 4;
				//pow(2,tableaus(act_tableau)->ODE_order)/tableaus(act_tableau)->ODE_order;
			}
		}

		if (old_stage != act_stage)
		{numsol.DestroyOldJacs();}

		if (old_stage != act_stage || reduce_order || steps_since_orderchange > 1.5*minsteps)
		{
			steps_since_orderchange = 0;
			min_newtonits = NLS_MaxNewtonSteps();
			newjacobians_since_orderchange = 0;
			last_jaccount = log.jaccount;
			reduce_order = 0;
		}
		if (old_stage != act_stage)
		{
			uo << "change ODE-order to = " << tableaus(act_tableau)->ODE_order << ENDL;
			/*		uo << "steps_per_jac=" << steps_per_jac <<
			", min_newtonits = " << min_newtonits << 
			", newjacs = " << newjacobians_since_orderchange << 
			", steps_since=" << steps_since_orderchange <<
			", reduce_order=" << reduce_order << ENDL;
			*/
		}
	}

	FinishStep();

	//return success; //1= succeeded, 0=failed
	return 1;

}

void TimeInt::FinishStep()
{
	// new!
	TItime += TIstep;
	
	EndTimeStep();

	//if (solset.writeresults) 
	//{
		//does not work for IO discrete elements
		//if (!FullAdaptive() && TIwritesolevery != 0)
		//{
		//	//interpolated output:
		//	TISaveState(TItemp_state);
		//	double t0 = TItime;
		//	IRK_Tableau* t = tableaus[act_tableau];
		//	int ns = t->nstages;
		//	const Vector& ct = t->c;

		//	while (TItime+TIstep >= (TIlastwritesol+TIwritesolevery-0.1*TIminstep))
		//	{
		//		TIlastwritesol+=TIwritesolevery;

		//		double tau = (TIlastwritesol-t0)/TIstep;
		//		//compute interpolant for states, start and end point
		//		mystatic Vector c;
		//		//mystatic Vector f;
		//		mystatic Vector xh;
		//		int offs = 0;
		//		int nip = ns;
		//		if (ct(1) != 0) {nip++; offs=1; }
		//		if (ct(ns) != 1) {nip++;}
		//		c.SetLen(nip);
		//		//f.SetLen(nip);
		//		c(1) = 0;
		//		c(nip) = 1;
		//		for (int i=2; i <= nip; i++)
		//		{
		//			if (!(i == nip && ct(ns) != 1)) c(i) = ct(i-offs);
		//		}
		//		//uo << "c=" << c << ENDL;

		//		TIx0.SetAll(0);
		//		//xi[i] contain stage solutions
		//		//TIx0start: solution at t_rel=0
		//		//TIx0save: solution at t_rel=1
		//		for (int i=1; i <= nip; i++)
		//		{
		//			double Li = LagrangeWeight(i,tau,c);
		//			//uo << "L" << i << "=" << Li << ENDL;
		//			if (i == 1) TIx0 += Li*TIx0start;
		//			else if (i == nip && ct(ns) != 1) TIx0 += Li*TItemp_state.x0;
		//			else
		//			{
		//				TIx0 += Li*xi[i-1];
		//			}
		//		}

		//		WriteSol();
		//	}
		//	TIRestoreState(TItemp_state);
		//	TItime += TIstep; //possibly this is too late ...
		//}
		//else 
		WriteSol();
	//}
	TIDrawAndStore();
}

void TimeInt::TIDrawAndStore()
{
	//redraw frequency: 0..off, 1..draw last frame, 2..100sec, 3..20sec, 4..2sec, 5..200ms, 6..50ms, 7..20ms, 8..every 10 frames, 9..every frame
	if (((solset.withgraphics >= 1 && solset.withgraphics <= 8) && (drawnow != 0)) || solset.withgraphics == 9)
	{
		if (GetCompTime() > 2 + lastloadresults || TIfinished!=0)
		{
			TIx0draw = TIx0;
			TIdrawtime = TItime;
			// Set XData vector
			TIdrawdata = TIdata;
		}
		drawnow=0;
		pCFB->ResultsUpdated(1); //redraw only, no results updated
	}

	pCFB->ResultsUpdated(2); //update results only, no redraw
			
	//$ PG 2012-6-18:[
	// user may pause calculation by pause-button (see CWCDriver3DDlg::OnButtonPause)
	// the pause functionality is primarily needed for save operations with the data manager during a simulation.
	// it is safest to pause the calculation immediately (and only) after the draw-configuration is copied from the computation-configuration (TIx0draw = TIx0, etc.).
	if (PauseCalculation())
	{
		pCFB->PausedComputation();
		int memo = bComputationIsInProgress;
		bComputationIsInProgress = 0;
		while (PauseCalculation()) {/*computation thread keeps stuck here*/ uo.pUI->SleepX(100); }
		bComputationIsInProgress = memo;
	}
	//$ PG 2012-6-18:]
}

int TimeInt::StoreResultsIsOn()
{
	if (solset.storedata == 0) return 0;
	else if (solset.storedata == -2) return 1;
	else if (solset.storedata == -1 && (GetTime()+0.1*TIminstep >= laststoredata+TImaxstep)) {laststoredata = GetTime(); return 1;}
	else if (solset.storedata > 0 && (GetTime()+0.1*TIminstep) >= laststoredata+solset.storedata) {laststoredata = GetTime(); return 1;}
	return 0;
}

int TimeInt::RemoveResults()
{
	return 0;
}

void TimeInt::StoreResults(DataSaver & storage, double &m_TimePoint)
{
	const Vector& v = GetSolVector();
	const Vector& d = GetDataVector();

	if (solset.sol_data_to_files)
	{
		TArray<double> data_container(3+v.Length()+d.Length());
		int idx = 1;
		
		m_TimePoint = GetTime();
		int vlen = (double) v.Length();
		int dlen = (double) d.Length();
		
		data_container(idx++) = m_TimePoint;
		data_container(idx++) = vlen;
		data_container(idx++) = dlen;
		for (int i=1; i<=vlen; i++)
		{
			data_container(idx++) = v(i);
		}		
		for (int i=1; i<=dlen; i++)
		{
			data_container(idx++) = d(i);
		}
		
		// store iteration number in array (which is written to file info.txt at the end of the simulation)
		//TInumber_of_solution_data_steps.Add(TIit);
		TInumber_of_solution_data_steps++;

		// store solution and data of this iteration in own file
		CMFile file(mystr(solset.sol_directory)+mystr("solution_data\\")+mystr(TInumber_of_solution_data_steps)+mystr(".dat"), TFMwrite, 1);
		for (int i=1; i<=data_container.Length(); i++)
		{
			file.RWbinaryDouble(data_container(i));
		}
	}
	else
	{
		storage.SetTime(GetTime());

		storage << v.GetLen();
		for (int i=1; i<=v.GetLen(); i++)
		{
			storage << v(i);
		}

		storage << d.Length();
		for (int i=1; i<=d.GetLen(); i++)
		{
			storage << d(i);
		}
	}
}

void TimeInt::LoadResults(DataLoader & loader, int m_TimePointNumber)
{
	lastloadresults = GetCompTime();

	if (solset.sol_data_to_files)
	{
		CMFile file(mystr(solset.sol_directory)+mystr("solution_data\\")+mystr(m_TimePointNumber)+mystr(".dat"), TFMread, 1);
		
		//int drawtimestep; file.RWbinaryDouble(drawtimestep);

		double drawtime; file.RWbinaryDouble(drawtime);
		SetDrawTime(drawtime);

		double len_double; file.RWbinaryDouble(len_double);
		int vlen = (int) len_double;
		TIx0draw.SetLen(vlen);
		TIx0draw.SetAll(0.);
		file.RWbinaryDouble(len_double);
		int dlen = (int) len_double;
		TIdrawdata.SetLen(dlen);
		TIdrawdata.SetAll(0.);
		for (int i=1; i<=vlen; i++)
		{
			file.RWbinaryDouble(TIx0draw(i));
		}
		TIdrawdata.SetAll(0.);
		for (int i=1; i<=dlen; i++)
		{
			file.RWbinaryDouble(TIdrawdata(i));
		}
	}
	else
	{
		SetDrawTime(loader.GetTime());

		int len;
		loader >> len;

		TIx0draw.SetLen(len);
		TIx0draw.SetAll(0.);

		for (int i=1; i<=len; i++)
		{
			loader >> TIx0draw(i);
		}

		loader >> len;

		TIdrawdata.SetLen(len);
		TIdrawdata.SetAll(0.);

		for (int i=1; i<=len; i++)
		{
			loader >> TIdrawdata(i);
		}
	}
}

int TimeInt::FirstStep()
{
	TIlastwritesol = 0;

	if (solset.writeresults)
	{
		double step = TIstep; //...necessary for output file
		TIstep = 0; //==>GetStepEndTime = TItime = 0 //...necessary for output file
		WriteSol();
		TIstep = step;//...necessary for output file
	}

	//Compute initial conditions for algebraic part
	return 1;
}

int TimeInt::GeneralIStep()
{
	//macht einen Schritt

	int ss = GetSystemSize();
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();
	int sos_rs = GetSecondOrderSize_RS();
	//int rs = GetResortSize();
	//if (GetTime() < 0.9*GetStepSize()) uo << "*****************\nWARNING: changed bandsize-6\n******************\n";
	numsol.SetBandSize(sos_rs - reducedbandsize); //for banded-solver, which part is banded!

	//UO() << "sos_rs=" << sos_rs << ", sos=" << sos << "\n";

	int rv = 1;

	//als Start für beide Lösungen den letzten Zustand verwenden 
	int i,j;
	IRK_Tableau* t = tableaus[act_tableau];
	int ns = t->nstages;

	//Kversion of Runge Kutta schemes:
	if (TIkv && t->implicit)
	{
		//TI-KVersion!!!!!!!!!!!
		int ns = t->nstages;
		int c0f = 0; //for lobattoIIIC und RadauIA methods!!!!
		if (t->c(1) == 0)
		{
			c0f = 1;
		}

		mystatic Vector GSxs;
		mystatic Vector ke[global_maxstages]; //explicit equations
		mystatic Vector ku[global_maxstages]; //sos equations, u
		mystatic Vector kv[global_maxstages]; //sos equations, v
		mystatic Vector evalf;
		mystatic Vector xh2;
		mystatic Vector xh3;
		mystatic Vector iv[global_maxstages];

		mystatic Vector g1e; //first order
		mystatic Vector g1u; //second order u
		mystatic Vector g1v; //second order v

		TIu0.SetLen(sos);
		TIu0e.SetLen(es);
		TIku0e.SetLen(es);
		TIv0.SetLen(sos);
		TIi0.SetLen(is);

		TIu0.Copy(TIx0,1,1,sos);        //needed for NLF
		TIv0.Copy(TIx0,1+sos,1,sos);    //needed for NLF
		TIu0e.Copy(TIx0,1+2*sos,1,es);  //needed for NLF
		TIi0.Copy(TIx0,1+es+2*sos,1,is);//needed for NLF

		TIkv0.SetLen(sos); 
		TIku0 = TIv0;
		TIkv0.Copy(TIk0,1+sos,1,sos);
		TIku0e.Copy(TIk0,1+2*sos,1,es);

		GSxs.SetLen((ns-c0f)*(sos+es+is));
		for (i = 1; i <= (ns-c0f); i++)
		{
			GSxs.Copy(TIkv0,1,1+(i-1)*sos,sos);
			GSxs.Copy(TIku0e,1,1+(ns-c0f)*sos+(i-1)*es,es);
			GSxs.Copy(TIi0,1,1+(ns-c0f)*(sos+es)+(i-1)*is,is);
		}

		impl_equ_flag = 0;
		//mystatic Vector xss = GSxs;
		if (!flag_applyjac)
		{
			numsol.NLSolveInfo() = TIstep;
			rv = numsol.NLSolve(GSxs); //call Newton
		} else
		{
			mystatic Vector F;
			mystatic Vector dx;
			F.SetLen(GSxs.GetLen());
			dx.SetLen(GSxs.GetLen());
			NLF(GSxs,F);
			rv = numsol.ApplyJac(F, dx); //apply only Jacobian
			GSxs -= dx;
			if (!rv) uo << "ApplyJac=0\n";
		}

		log.jaccount += numsol.GetJacCount();
		log.TInewtonitsum += numsol.GetNewtonIts();
		log.TInewtonit = Maximum(numsol.GetNewtonIts(), log.TInewtonit);
		TIcondnum = numsol.GetJacCondnum();


		if (c0f)
		{
			ku[1] = TIku0;
			kv[1] = TIkv0;
			ke[1] = TIku0e;
		}

		for (i = 1; i <= (ns-c0f); i++)
		{
			kv[i+c0f].SetLen(sos);
			kv[i+c0f].Copy(GSxs,1+sos*(i-1),1,sos);
		}

		for (i = 1; i <= (ns-c0f); i++)
		{
			ke[i+c0f].SetLen(es);
			ke[i+c0f].Copy(GSxs,1+sos*(ns-c0f)+(i-1)*es,1,es);
		}

		//Compute ku[i]
		for (i = 1; i <= (ns-c0f); i++)
		{
			ku[i+c0f] = TIv0;
			for (j = 1; j <= ns; j++)
			{
				xh = kv[j];
				xh *= TIstep * t->A(i+c0f,j);
				ku[i+c0f] += xh;
			}
		}

		//Evaluation step with coefficients b_i
		//with decomposition gu/gv

		g1u = TIu0;
		g1v = TIv0;
		g1e = TIu0e;
		xh.SetLen(sos);
		xh2.SetLen(es);

		for (i = 1; i <= ns; i++)
		{
			xh = ku[i];
			xh *= TIstep * t->b(i);
			g1u += xh;

			xh = kv[i];
			xh *= TIstep * t->b(i);
			g1v += xh;

			xh2 = ke[i];
			xh2 *= TIstep * t->b(i);
			g1e += xh2;
		}

		//compute new TIx0:
		TIx0.Copy(g1u,1,1,sos);
		TIx0.Copy(g1v,1,1+sos,sos);
		TIx0.Copy(g1e,1,1+2*sos,es);

		//generate new TIk0 for next timestep:
		TIk0.Copy(ku[ns],1,1,sos);
		TIk0.Copy(kv[ns],1,1+sos,sos);
		TIk0.Copy(ke[ns],1,1+2*sos,es);
		TIk0.Copy(GSxs,(ns-c0f)*(sos+es)+1+is*((ns-c0f)-1),1+2*sos+es,is); //iv[ns]

		//compute final value of implicit solution, extrapolate
		//maybe there is some error, esp for c1=0 ???
		if (t->Ainvertable && t->c(ns) != 1 && is!=0)
		{
			double k0 = 1. - t->bAinv.Sum();
			xh3.SetLen(is);

			if (c0f) {k0+=t->bAinv(1);}
			TIi0 *= k0;
			for (i = 1; i <= (ns-c0f); i++)
			{
				xh3.Copy(GSxs,1+(ns-c0f)*(sos+es)+is*(i-1),1,is);
				xh3 *= t->bAinv(i+c0f);
				TIi0 += xh3;
			}
			TIx0.Copy(TIi0,1,1+2*sos+es,is); //iv[ns]
		}
		else
		{
			TIx0.Copy(GSxs,(ns-c0f)*(sos+es)+1+is*((ns-c0f)-1),1+2*sos+es,is); //iv[ns]
		}
	}
	//pure first order equations (with constraints):
	else if (GetSecondOrderSize() == 0 && t->implicit)
	{

		mystatic Vector gv[global_maxstages];
		mystatic Vector iv[global_maxstages];
		mystatic Vector evalf;
		mystatic Vector GSxs;
		mystatic Vector gb;
		mystatic Vector gx;

		TIu0e.SetLen(es);
		TIu0e.Copy(TIx0,1,1,es);

		TIi0.SetLen(is);
		TIi0.Copy(TIx0,es+1,1,is);

		GSxs.SetLen(ns*(es+is));
		//find start-vector for newton
		for (i = 0; i < ns; i++)
		{
			GSxs.Copy(TIu0e,1,1+es*i,es);
		}
		for (i = 0; i < ns; i++)
		{
			GSxs.Copy(TIi0,1,1+es*ns+is*i,is);
		}

		impl_equ_flag = 0;

		rv = numsol.NLSolve(GSxs); //call Newton

		log.jaccount += numsol.GetJacCount();
		log.TInewtonitsum += numsol.GetNewtonIts();
		log.TInewtonit = Maximum(numsol.GetNewtonIts(), log.TInewtonit);
		TIcondnum = numsol.GetJacCondnum();

		for (i = 1; i <= ns; i++)
		{
			gv[i].SetLen(es);
			gv[i].Copy(GSxs,1+es*(i-1),1,es);
			//gv[i] = GSxs.SubVector(1+es*(i-1), i*es);
		}
		for (i = 1; i <= ns; i++)
		{
			iv[i].SetLen(is);
			iv[i].Copy(GSxs,1+es*ns+is*(i-1),1,is);
			//iv[i] = GSxs.SubVector(ns*es+1+is*(i-1), ns*es+is*i);
		}

		gb.SetLen(es);
		gb = TIu0e;

		evalf.SetLen(es);
		gx.SetLen(es+is);

		for (i = 1; i <= ns; i++)
		{
			gx.Copy(gv[i],1,1,es);
			gx.Copy(iv[i],1,es+1,is);
			EvalF(gx,evalf,TItime+TIstep*t->c(i));
			evalf.Mult(TIstep*t->b(i));
			gb += evalf;
		}      

		TIx0.Copy(gb,1,1,es);
		TIx0.Copy(iv[ns],1,es+1,is);
		EvalF(TIx0, evalf,TItime+TIstep); //write new state into system variables

	}
	else if (!t->implicit && is==0)
	{
		//explicit Runge Kutta!!!

		//very inefficient mass matrix!!!
		mystatic Matrix mass[global_maxstages];

		evalf.SetLen(sos);
		TIu0.SetLen(sos);
		TIu0.Copy(TIx0,1,1,sos);
		TIv0.SetLen(sos);
		TIv0.Copy(TIx0,1+sos,1,sos);
		TIu0e.SetLen(es);
		TIu0e.Copy(TIx0,1+2*sos,1,es);

		//build gu, gv, gue
		gu[1] = TIu0;
		gv[1] = TIv0;
		gue[1] = TIu0e;

		for (i=1; i <= ns; i++)
		{
			//evaluate right hand sides for first and second order equations
			kv[i].SetLen(sos);
			kue[i].SetLen(es);
			xi[i].SetLen(2*sos+es);

			xi[i].Copy(gu[i],1,1,sos);
			xi[i].Copy(gv[i],1,1+sos,sos);
			xi[i].Copy(gue[i],1,1+2*sos,es);

			//compute RK-f(x,t) for gu,gue and gv 
			ku[i] = gv[i];
			evalfe.SetLen(es);
			if (es != 0)
			{
				EvalF(xi[i],kue[i],TItime+TIstep*t->c(i));
			}

			if (sos != 0)
			{
				mass[i].SetSize(sos,sos);
				mass[i].FillWithZeros();
				evalf.SetLen(sos);
				EvalM(xi[i], mass[i], TItime+TIstep*t->c(i));
				EvalF2(xi[i], evalf, TItime+TIstep*t->c(i));

				if (!mass[i].Solve(evalf,kv[i])) {uo << "Error in GeneralIStep: Mass matrix is singular!" << ENDL; return 0;}
			}

			//compute next stage
			if (i < ns)
			{
				gu[i+1] = TIu0;
				gv[i+1] = TIv0;
				gue[i+1] = TIu0e;
				for (j = 1; j <= i; j++)
				{
					gu[i+1] += TIstep*t->A(i+1,j)*ku[j];
					gv[i+1] += TIstep*t->A(i+1,j)*kv[j];
					gue[i+1]+= TIstep*t->A(i+1,j)*kue[j];
				}
			}
		}

		//compute evaluation step
		for (i = 1; i <= ns; i++)
		{
			TIu0  += TIstep*t->b(i)*ku[i];
			TIv0  += TIstep*t->b(i)*kv[i];
			TIu0e += TIstep*t->b(i)*kue[i];
		}

		//compute new TIx0:
		TIx0.Copy(TIu0,1,1,sos);
		TIx0.Copy(TIv0,1,1+sos,sos);
		TIx0.Copy(TIu0e,1,1+2*sos,es);
	}
	else
	{
		uo << "Mixed first-second order system not yet implemented!!!" << ENDL;
	}


	//TInewtonitsum += TInewtonit;
	if (rv != 1) 
	{
		//uo << "NL-Solve (Newton) failed" << ENDL; 
		rv=0; 
	}
	return rv;
}


//Compute (df/dx)
void TimeInt::LocalJacobianF(Matrix& m, Vector& x)
{
	int i,j;
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int ss = GetSystemSize();
	if (es == 0) return;

	locjacf0.SetLen(es);
	locjacf1.SetLen(es);
	locjacf2.SetLen(es);

	double numdiffepsi = numsol.NumDiffepsi();
	double eps;
	NLS_SetJacCol(0);

	double t = TItime+0.5*TIstep;
	if (!numsol.SymmetricJacobian())
		EvalF(x,locjacf0,t);

	for (i = 1; i <= ss; i++)
	{
		NLS_SetJacCol(i);
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (numsol.SymmetricJacobian())
		{
			x(i) += eps;
			EvalF(x,locjacf1,t);

			x(i) -= 2*eps;
			EvalF(x,locjacf2,t);
			x(i) += eps;
		}
		else
		{
			x(i) += 2*eps;
			EvalF(x,locjacf1,t);
			x(i) -= 2*eps;
			//EvalF(x,locjacf2,t);
			locjacf2 = locjacf0;
		}

		for (j=1; j<=es;j++)
		{
			m(j,i)=0.5/eps*(locjacf1(j)-locjacf2(j));
		}
	}
	NLS_SetJacCol(0);
}

//Compute (df2/dx)
void TimeInt::LocalJacobianF2(Matrix& m, Vector& x)
{
	int i,j;
	int sos= GetSecondOrderSize();
	if (sos == 0) return;
	int ss = GetSystemSize();
	int is = GetImplicitSize();

	locjacf20.SetLen(sos);
	locjacf21.SetLen(sos);
	locjacf22.SetLen(sos);


	double numdiffepsi = numsol.NumDiffepsi();
	double eps;

	double t = TItime+0.5*TIstep;
	NLS_SetJacCol(0);

	if (!numsol.SymmetricJacobian())
		EvalF2(x,locjacf20,t);
	double storex;

	//uo << "x0=" << x << "\n";

	for (i = 1; i <= ss; i++)
	{
		NLS_SetJacCol(i);
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (numsol.SymmetricJacobian())
		{
			storex = x(i);
			x(i) += eps;
			EvalF2(x,locjacf21,t);

			x(i) -= 2*eps;
			EvalF2(x,locjacf22,t);
			x(i) = storex;
		}
		else
		{
			x(i) += 2*eps;
			EvalF2(x,locjacf21,t);

			x(i) -= 2*eps;
			locjacf22 = locjacf20;
		}

		for (j=1; j<=sos;j++)
		{
			m(j,i)=0.5/eps*(locjacf21(j)-locjacf22(j));
		}
		//uo << "jac1=" << 0.5/eps*(locjacf21-locjacf22) << "\n";
		//uo << "m=" << m << "\n";
	}
	NLS_SetJacCol(0);
}

//Compute (df2/dx)
void TimeInt::LocalJacobianF2(SparseMatrix& m, Vector& x)
{
	//uo << "m=" << m.Getrows() << ", " << m.Getcols() << "\n";

	int i,j;
	int sos= GetSecondOrderSize();
	if (sos == 0) return;
	int ss = GetSystemSize();
	int is = GetImplicitSize();

	slocjacf20.SetLen(sos,32);
	slocjacf21.SetLen(sos,32);
	slocjacf22.SetLen(sos,32);

	locjacf20.SetLen(sos); //for un-symmetric computation
	locjacf21.SetLen(sos);
	locjacf22.SetLen(sos);


	temprowind.SetLen(sos);
	for (i = 1; i <= temprowind.Length(); i++) temprowind(i) = 0;

	double numdiffepsi = numsol.NumDiffepsi();
	double eps;

	double t = TItime+0.5*TIstep;
	static TArray<int> elems(100);
	static TArray<int> oldelems(100);
	elems.SetLen(0);
	oldelems.SetLen(0);

	NLS_SetJacCol(0);
	m.FillWithZeros(); //no costs, but necessary, because not everything filled in!!!

	if (!numsol.SymmetricJacobian())
		EvalF2(x,locjacf20,t);

	double storex;
	int nprint = 6000;
	if (ss > nprint) uo << "JacF2(" << ss <<")";

	//int testcnt = 0;


	for (i = 1; i <= ss; i++)
	{
		if (ss > nprint && i%1000==0) uo << ".";
		//if (i%1000 == 0) uo << "i=" << i << "\n";
		NLS_SetJacCol(i);
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		//eps2 = eps;

		if (i==sos+1 || i==ss-is+1) oldelems.SetLen(0);

		if (numsol.SymmetricJacobian())
		{
			if (0)
			{
				storex = x(i); 
				x(i) += 2*eps;
				EvalF2(x,slocjacf21,t,temprowind, tempclearind, elems);
				if (!IsEqual(elems, oldelems))	//if equal then use old slocjacf22!!!
				{
					oldelems = elems;
					x(i) -= 2*eps;
					EvalF2(x,slocjacf22,t,temprowind, tempclearind, elems);
				}
				x(i) = storex;
			}
			else
			{
				storex = x(i); 
				x(i) += eps;
				EvalF2(x,slocjacf21,t,temprowind, tempclearind, elems);
				x(i) -= 2*eps;
				EvalF2(x,slocjacf22,t,temprowind, tempclearind, elems);
				x(i) = storex;
			}
		}
		else
		{
			x(i) += 2*eps;
			EvalF2(x,locjacf21,t);

			x(i) -= 2*eps;
		}

		if (numsol.SymmetricJacobian())
		{
			for (j=1; j<=slocjacf21.NEntries();j++)
			{
				slocjacf21.Entry(j) = (slocjacf21.Entry(j)-slocjacf22.Entry(j))/(2*eps);
			}
			m.SetColVector(slocjacf21,i);
		}
		else
		{
			for (j=1; j <= sos;j++)
			{
				locjacf21(j)=(locjacf21(j)-locjacf20(j))/eps;
			}
			m.SetColVector(locjacf21,1,sos,i);
		}
	}
	if (ss > nprint) uo << "\n";

	//uo << "m=" << m << "\n";
	//uo << "jac-evalf=" << testcnt << "\n";
	NLS_SetJacCol(0);
}

void TimeInt::LocalJacobianG(Matrix& m, Vector& x)
{
	int is = GetImplicitSize();
	if (is == 0) return;

	int i,j;
	int es = GetFirstOrderSize();
	int ss = GetSystemSize();

	locjacg0.SetLen(is);
	locjacg1.SetLen(is);
	locjacg2.SetLen(is);

	double numdiffepsi = numsol.NumDiffepsi();
	double eps;

	double t = TItime+0.5*TIstep;
	NLS_SetJacCol(0);
	if (!numsol.SymmetricJacobian())
		EvalG(x,locjacg0,t);

	for (i = 1; i <= ss; i++)
	{
		NLS_SetJacCol(i);
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (numsol.SymmetricJacobian())
		{
			x(i) += eps;
			EvalG(x,locjacg1,t);

			x(i) -= 2*eps;
			EvalG(x,locjacg2,t);
			x(i) += eps;
		}
		else
		{
			x(i) += 2*eps;
			EvalG(x,locjacg1,t);
			x(i) -= 2*eps;
			locjacg2 = locjacg0;
		}

		for (j=1; j<=is;j++)
		{
			m(j,i)=0.5/eps*(locjacg1(j)-locjacg2(j));
		}
	}
	NLS_SetJacCol(0);

}

void TimeInt::LocalJacobianG(SparseMatrix& m, Vector& x)
{
	int i,j;
	int sos= GetSecondOrderSize();
	int ss = GetSystemSize();
	int is = GetImplicitSize();
	if (is == 0) return;

	slocjacg0.SetLen(is,32);
	slocjacg1.SetLen(is,32);
	slocjacg2.SetLen(is,32);

	locjacg0.SetLen(is); //for un-symertric computation
	locjacg1.SetLen(is);
	locjacg2.SetLen(is);

	temprowind.SetLen(is);
	for (i = 1; i <= temprowind.Length(); i++) temprowind(i) = 0;

	double numdiffepsi = numsol.NumDiffepsi();
	double eps;

	double t = TItime+0.5*TIstep;
	NLS_SetJacCol(0);
	m.FillWithZeros(); //no costs, but necessary, because not everything filled in!!!

	if (!numsol.SymmetricJacobian())
		EvalG(x,locjacg0,t);

	for (i = 1; i <= ss; i++)
	{
		NLS_SetJacCol(i);
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (numsol.SymmetricJacobian())
		{
			x(i) += eps;
			EvalG(x,slocjacg1,t,temprowind, tempclearind);

			x(i) -= 2*eps;
			EvalG(x,slocjacg2,t,temprowind, tempclearind);
			x(i) += eps;
		}
		else
		{
			x(i) += 2*eps;
			EvalG(x,locjacg1,t);
			x(i) -= 2*eps;
		}

		if (numsol.SymmetricJacobian())
		{
			for (j=1; j<=slocjacg1.NEntries();j++)
			{
				slocjacg1.Entry(j) = 0.5/eps*(slocjacg1.Entry(j)-slocjacg2.Entry(j));
			}
			m.SetColVector(slocjacg1,i);
		}
		else
		{
			for (j=1; j <= is;j++)
			{
				locjacg1(j)=0.5/eps*(locjacg1(j)-locjacg0(j));
			}
			m.SetColVector(locjacg1,1,is,i);
		}
	}
	NLS_SetJacCol(0);
}



void TimeInt::LocalJacobianM(Matrix& m, Vector& x)
{
	NLS_SetJacCol(0);
	int i,j;

	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int ss = GetSystemSize();
	int sos= GetSecondOrderSize();
	if (sos == 0) return;

	locjacf20.SetLen(sos);
	locjacf21.SetLen(sos);
	locjacf22.SetLen(sos);

	mystatic Vector TIkv0(sos,sos);
	mystatic Matrix locm;
	locm.SetSize(sos,sos);
	locm = 0.;

	TIkv0.SetLen(sos);
	TIkv0.Copy(TIk0,1+sos,1,sos);


	double numdiffepsi = numsol.NumDiffepsi();
	double eps;

	double t = TItime+0.5*TIstep;
	NLS_SetJacCol(0);

	if (!numsol.SymmetricJacobian())
	{
		EvalM(x,locm,t);
		locjacf20 = locm*TIkv0;
	}

	for (i = 1; i <= sos; i++)
	{
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (numsol.SymmetricJacobian())
		{
			x(i) += eps;
			EvalM(x,locm,t);
			locjacf21 = locm*TIkv0;

			x(i) -= 2*eps;
			EvalM(x,locm,t);
			locjacf22 = locm*TIkv0;

			x(i) += eps;
		}
		else
		{
			x(i) += 2*eps;
			EvalM(x,locm,t);
			locjacf21 = locm*TIkv0;

			x(i) -= 2*eps;
			locjacf22 = locjacf20;
		}

		for (j=1; j<=sos;j++)
		{
			m(j,i)=0.5/eps*(locjacf21(j)-locjacf22(j));
		}
	}
	NLS_SetJacCol(0);
}

//Compute Full Jacobian
void TimeInt::FullJacobian(Matrix& m, Vector& x)
{
}

void TimeInt::ComputeStiffnessAndDampingMatrix(SparseMatrix& k, SparseMatrix& d, Vector& x)
{
	TMStartTimer(19);
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();
	int ss = GetSystemSize();
	
	// $ MSax 2013-07-25 : [ removed because not used
	//mystatic Vector f1;
	//mystatic Vector f2;
	//f1.SetLen(x.GetLen());
	//f2.SetLen(x.GetLen());
	// $ MSax 2013-07-25 : ] removed because not used

	NLS_SetJacCol(0);

	int m_initsize = MaxSparseBandwidth()*2;//2; //initial sparse matrix size
	mlocgs.SetSize(is,ss,m_initsize);

	xx = x; // positions AND velocities !!

	k.SetSize((sos+es+is),(sos+es+is), m_initsize);
	k.FillWithZeros();
	d.SetSize((sos+es+is),(sos+es+is), m_initsize);
	d.FillWithZeros();

	mlocf2s.SetSize(sos,ss,m_initsize);
	if (sos) LocalJacobianF2(mlocf2s,xx);

	//Matrix K(sos,ss);
	//LocalJacobianF2(K,xx); //!!!!!!!!!!!
	//uo << "KD=" << K << "\n";


	k.AddSubmatrix(mlocf2s,1,1,1,1,sos,sos,1.);
	d.AddSubmatrix(mlocf2s,1,1+sos,1,1,sos,sos,1.);
	
	////Stiffness part
	//m.AddSubmatrix(mlocf2s,1,1+sos,1,1,sos,sos,1);
	//m.AddSubmatrix(mlocf2s,1,1+2*sos+es,1,1+sos+es,sos,is,1);

	//if (is) LocalJacobianG(mlocgs,xx);
	//if (es) {UO() << "Sparse static solver not implemented for first order ODEs!\n"; return; }

	////Constraint part
	//m.AddSubmatrix(mlocgs, 1,1,sos+es+1,1,is,sos,1);
	//m.AddSubmatrix(mlocgs, 1,1+2*sos+es,sos+es+1,1+sos+es,is,is,1);

	TMStopTimer(19);
}

void TimeInt::ComputeGyroscopicMatrix(SparseMatrix& gy)//, Vector& x) // $ MSax 2013-07-25 : added
{
}

void TimeInt::StaticJacobian(SparseMatrix& m, Vector& x)
{
	int newmode = solset.static_experimental_sparse_jacobian; //new mode using less memory (using mbs->StaticJacobianF2(...))

	TMStartTimer(19);
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();
	int ss = GetSystemSize();
	
	mystatic Vector f1;
	mystatic Vector f2;
	f1.SetLen(x.GetLen());
	f2.SetLen(x.GetLen());

	IVector usage_per_row;
	ComputeSparseMatrixUsage(usage_per_row);

	NLS_SetJacCol(0);

	int m_initsize = MaxSparseBandwidth()*2;//2; //initial sparse matrix size
	mlocgs.SetSize(is,ss,m_initsize);

	xx.SetLen(ss);
	//xx = TIx0;
	for (int i=1; i <= sos; i++)
	{
		xx(i) = x(i);
		xx(i+sos) = 0.; //velocities set to zero!
	}
	for (int i=1; i <= es; i++)
		xx(i+2*sos) = x(i+sos);
	for (int i=1; i <= is; i++)
		xx(i+2*sos+es) = x(i+sos+es);

	//uo << "JacobianF2" << ", initsize=" << m_initsize << "\n";
	if (newmode)
	{
		m.SetSizePerColumn((sos+es+is),(sos+es+is), usage_per_row);
		m.FillWithZeros();

		//fill in stiffness matrix directly:
		if (sos) StaticJacobianF2(m,xx);
	}
	else
	{
		m.SetSize((sos+es+is),(sos+es+is), m_initsize);
		m.FillWithZeros();

		mlocf2s.SetSize(sos,ss,m_initsize);
		if (sos) LocalJacobianF2(mlocf2s,xx);
		//Stiffness part
		m.AddSubmatrix(mlocf2s,1,1,1,1,sos,sos,1);
		m.AddSubmatrix(mlocf2s,1,1+2*sos+es,1,1+sos+es,sos,is,1);
	}

	if (is) LocalJacobianG(mlocgs,xx);
	if (es) {UO() << "Sparse static solver not implemented for first order ODEs!\n"; return; }

	//Constraint part
	m.AddSubmatrix(mlocgs, 1,1,sos+es+1,1,is,sos,1);
	m.AddSubmatrix(mlocgs, 1,1+2*sos+es,sos+es+1,1+sos+es,is,is,1);

	//UO() << "sparse jac=" << m << "\n";

	//double msize = m.GetLAlloc();
	//double mused = m.CountEntries();
	//UO(UO_LVL_all) << "sparse matrix entries=" << mused << "\n";
	//UO(UO_LVL_all) << "sparse matrix allocated=" << msize << "\n";

	//UO(UO_LVL_all) << "sparse matrix size=" << msize*8/1.e6 << "MB\n";
	//UO(UO_LVL_all) << "sparse matrix used=" << mused*8/1.e6 << "MB\n";
	//UO(UO_LVL_all) << "factor=" << msize/mused << "\n";

	//for (int i=1; i <= m.Getrows(); i++)
	//{
	//	if (m.RowLen(i) > usage_per_row(i)) 
	//	{
	//		UO(UO_LVL_all) << "rowlen(" << i << ")=" << m.RowLen(i) << ", est.=" << usage_per_row(i);
	//		UO(UO_LVL_all) << " WARNING!";
	//		UO(UO_LVL_all) << "\n";
	//	}

	//}

	//uo << "ready\n";


	//UO(UO_LVL_all) << "mdiff=" << (mf-mf2).MaxNorm() << "\n";
	//UO(UO_LVL_all) << "mnew=" << mf2 << "\n";
	//UO(UO_LVL_all) << "inverse returns=" << mf.Invert2() << "\n";

	//UO().InstantMessageText("wait");

	//UO() << "jac_F2=" << mlocf2s << "\n";
	//UO() << "sparse jac=" << mf << "\n";

	TMStopTimer(19);
}


int fullnewtonwarnedsparse = 0;
//Compute approximate Jacobian (if modified newton) or call base function
void TimeInt::Jacobian(SparseMatrix& m, Vector& x)
{
	//UO() << "jac\n";

	SetJacobianComputationFlag(1);
	//return;
	NLS_SetJacCol(0);
	IRK_Tableau* t = tableaus[act_tableau];
	int ns = t->nstages;

	if (DoStaticComputation())
	{
		//uo << "jacobian...";
		StaticJacobian(m,x);
		//uo << "computed!\n";
	}
	else
	{
		if (!numsol.ModifiedNewtonActive() && !fullnewtonwarnedsparse)
		{
			uo << "ERROR: Full Newton called for Sparse Matrices!!!\n";
			fullnewtonwarnedsparse = 1;

			//SetJacobianComputationFlag(0); //do this before return!!!
			//return;
		}


		TMStartTimer(19);

		int m_initsize = MaxSparseBandwidth(); //initial sparse matrix size
		//UO() << "maxbandwidth=" << m_initsize << "\n";
		int i;
		int j;

		int c0f = 0; //for lobattoIIIC und RadauIA methods!!!!
		if (t->c(1) == 0)
		{
			c0f = 1;
		}

		int is = GetImplicitSize();
		int es = GetFirstOrderSize();
		int sos = GetSecondOrderSize();
		int ss = GetSystemSize();
		double tau = TIstep;
		double tau2 = TIstep*TIstep;

		xx.SetLen(ss);
		xx = TIx0;

		mlocms.SetSize(sos,sos,m_initsize); //usually 32
		mlocf2s.SetSize(sos,ss,m_initsize);
		mlocgs.SetSize(is,ss,m_initsize);

		int locjacm_on = 0; //compute dM/TIkv ... normally not necessary!!!

		mlocf.SetSize(es,ss);
		if (locjacm_on) mlocm2.SetSize(sos,sos);

		diag.SetSize(es,es);
		diag.FillWithZeros();
		for (i=1; i<=es; i++)
		{
			diag(i,i) = 1.;
		}

		if (es) 
		{
			LocalJacobianF(mlocf,xx);
		}

		if (sos) 
		{
			LocalJacobianF2(mlocf2s,xx);
			EvalM(xx,mlocms,TItime+0.5*TIstep);

			if (locjacm_on) LocalJacobianM(mlocm2,xx);
		}

		if (is)
		{
			LocalJacobianG(mlocgs,xx);
		}

		TMStopTimer(19);
		TMStartTimer(4);

		m.SetSize((ns-c0f)*(sos+es+is),(ns-c0f)*(sos+es+is), m_initsize*(ns-c0f));
		m.FillWithZeros();


		for(i=1; i<= (ns-c0f); i++)
		{
			for(j=1; j<= (ns-c0f); j++)
			{
				int offri = sos*(i-1)+1;
				int offcj = sos*(j-1)+1;
				int offrif = sos*(ns-c0f)+es*(i-1)+1;
				int offcjf = sos*(ns-c0f)+es*(j-1)+1;
				int offrig = (sos+es)*(ns-c0f)+is*(i-1)+1;
				int offcjg = (sos+es)*(ns-c0f)+is*(j-1)+1;

				//Part A: sos
				if (i==j) m.AddSubmatrix(mlocms,1,1,offri,offcj,sos,sos,-1);
				if (locjacm_on && i==j) m.AddSubmatrix(mlocm2,1,1,offri,offcj,sos,sos,-1*tau2*t->A2(i+c0f,j+c0f));

				m.AddSubmatrix(mlocf2s,1,1,         offri,offcj,sos,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2s,1,1+sos,     offri,offcj,sos,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2s,1,1+sos*2,   offri,offcjf,sos,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocf2s,1,1+sos*2+es,offri,offcjg,sos,is,1);

				//Part B: es
				m.AddSubmatrix(mlocf,1,1,         offrif,offcj,es,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf,1,1+sos,     offrif,offcj,es,sos,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(diag,1,1,offrif,offcjf,es,es,-1);
				m.AddSubmatrix(mlocf,1,1+sos*2,   offrif,offcjf,es,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocf,1,1+sos*2+es,offrif,offcjg,es,is,1);

				//Part C: is
				m.AddSubmatrix(mlocgs,1,1,         offrig,offcj,is,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocgs,1,1+sos,     offrig,offcj,is,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocgs,1,1+sos*2,   offrig,offcjf,is,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocgs,1,1+sos*2+es,offrig,offcjg,is,is,1);

			}
		}

		/*
		uo << "sos=" << sos << ",es=" << es << ",is=" << is << ",ns=" << ns << ",c0f=" << c0f << ",m=" << m.Getrows() << "x" << m.Getcols() << "\n";
		uo << "m-alloc=" << m.GetLAlloc() << "\n";
		uo << "mlocms-alloc=" << mlocms.GetLAlloc() << "\n";
		uo << "mlocf2s-alloc=" << mlocf2s.GetLAlloc() << "\n";
		uo << "mlocgs-alloc=" << mlocgs.GetLAlloc() << "\n";
		*/
		//uo << "K=" << mlocf2s << "\n";
		//uo << "M=" << mlocms << "\n";
	}

	//if (!system_matrices_written && solset.write_MK2file)  //$ PG 2012-2-3: modified (see next line)
	if (!system_matrices_written && GetOptions()->LoggingOptions()->WriteMassAndStiffnessMatrix())
	{
		system_matrices_written= 1;

		mystr dir = solset.sol_directory;
		mystr fileload = "loadMK.m";
		mystr fileK = "Kdat.dat";
		mystr fileM = "Mdat.dat";

		ofstream fout((dir+fileload).c_str());
		fout << "load Kdat.dat\nK = spconvert(Kdat);\n";
		ofstream Kdat((dir+fileK).c_str());
		ofstream Mdat((dir+fileM).c_str());

		if (!DoStaticComputation())
		{
			int dim = mlocms.Getcols();
			SparseMatrix KS;
			KS.CopyFrom(mlocf2s,1,1,dim,dim);
			KS.PrintToMatlabSparse(Kdat);

			mlocms.PrintToMatlabSparse(Mdat);//mlocms=massenmatrix, schreibe in Mdat-File im Matlab-Format

			fout << "load Mdat.dat\nM = spconvert(Mdat);\n";
			//uo << "M=" << mlocms.GetMatrix() << "\n";
			//uo << "K=" << KS.GetMatrix() << "\n";
		}
		else
		{
			m.PrintToMatlabSparse(Kdat);
			//uo << "K=" << m << "\n";
		}
		fout << flush;

		UO() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		UO() << "M and K written to file '" << (dir+fileM) << "'\nand '" << (dir+fileK) << "'!\n";
		UO() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

		//PrintMatrix01(m.GetMatrix());
	}

	SetJacobianComputationFlag(0);
}

//Compute approximate Jacobian (if modified newton) or call base function
int fulljacobianwarned = 0;
void TimeInt::Jacobian(Matrix& m, Vector& x)
{
	SetJacobianComputationFlag(1);

	NLS_SetJacCol(0);
	IRK_Tableau* t = tableaus[act_tableau];
	int ns = t->nstages;

	if (DoStaticComputation())
	{
		TMStartTimer(3);
		int i,j;
		mystatic Vector f1;
		mystatic Vector f2;
		f1.SetLen(x.GetLen());
		f2.SetLen(x.GetLen());

		m.SetSize(x.GetLen(),x.GetLen());
		
		int usefullnewton = 1; //0 is faster

		if (usefullnewton)
		{
			double numdiffepsi = this->NLS_NumDiffepsi();
			NLS_SetJacCol(0);
			NLS_SetJacFullNewton(1);

			FullNewtonCnt()++;

			if (!NLS_SymmetricJacobian())
			{
				NLF(x,f2);
			}
			double eps;
			double xstore;
			//UO() << "len=" << m.Getcols() << "\n";
			for (i = 1; i <= m.Getcols(); i++)
			{
				eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
				if (NLS_SymmetricJacobian())
				{
					int jc = i;
					if (i > GetSecondOrderSize()) jc += GetSecondOrderSize();
					NLS_SetJacCol(jc);

					xstore = x(i);
					x(i) += eps;
					NLF(x,f1);

					x(i) -= 2*eps;
					NLF(x,f2);

					//$ PG 2013-4-24: [ debugging output
					//int gsize=f2.Length()-GetSecondOrderSize();
					//Vector differenc(gsize);
					//for(int j=1; j<=gsize; j++) differenc(j) = f2(GetSecondOrderSize()+j) - f1(GetSecondOrderSize()+j);
					//UO(UO_LVL_dbg1) << "i=" <<i<<", diff = " << differenc << "\n";
					//$ PG 2013-4-24: ]

					x(i) = xstore;
				}
				else
				{
					xstore = x(i);
					x(i) += 2*eps;
					NLF(x,f1);
					x(i) = xstore;
				}

				for (j=1; j<=f1.GetLen(); j++)
				{
					m(j,i)=0.5/eps*(f1(j)-f2(j));
				}
			}

			NLS_SetJacFullNewton(0);
			NLS_SetJacCol(0);
		}
		else
		{
			int is = GetImplicitSize();
			int es = GetFirstOrderSize();
			int sos = GetSecondOrderSize();
			int ss = GetSystemSize();

			mlocf.SetSize(es,ss);
			mlocf2.SetSize(sos,ss);
			mlocg.SetSize(is,ss);

			xx.SetLen(ss);
			for (int i=1; i <= sos; i++)
			{
				xx(i) = x(i);
				xx(i+sos) = 0.; //velocities set to zero!
			}
			for (int i=1; i <= es; i++)
				xx(i+2*sos) = x(i+sos);
			for (int i=1; i <= is; i++)
				xx(i+2*sos+es) = x(i+sos+es);

			if (sos) 
			{
				LocalJacobianF2(mlocf2,xx);
			}
			if (es) 
			{
				LocalJacobianF(mlocf,xx);
			}
			if (is)
			{
				LocalJacobianG(mlocg,xx);
			}

			int j1;
			for (int j=1; j <= sos+es+is; j++)
			{
				j1 = j;
				if (j > sos) j1+=sos;
				for (int i = 1; i <= sos; i++)
					m(i,j) = mlocf2(i,j1);

				for (int i = 1; i <= es; i++)
					m(i+sos,j) = mlocf(i,j1);

				for (int i = 1; i <= is; i++)
					m(i+sos+es,j) = mlocg(i,j1);
			}

		}
		TMStopTimer(3);

		//UO() << "jac_static=" << m << "\n";
		//uo << "Jac=\n"; PrintMatrix01(m);
		
		//if (!system_matrices_written && solset.write_MK2file)  //$ PG 2012-2-3: modified (see next line)
		if (!system_matrices_written && GetOptions()->LoggingOptions()->WriteMassAndStiffnessMatrix())
		{
			system_matrices_written= 1;

			mystr dir = solset.sol_directory;
			mystr fileK = "Kdat.dat";
			mystr fileload = "loadK.m";

			ofstream fout((dir+fileload).c_str());
			ofstream Kdat((dir+fileK).c_str());

			SparseMatrix KS(m);
			KS.PrintToMatlabSparse(Kdat);
			fout << "load Kdat.dat\nK = spconvert(Kdat);\n";

			uo << "jacobian written to '" << (dir+fileK) << "\n";
		}

		//uo << "Jac=\n"; PrintMatrix01(m);

		SetJacobianComputationFlag(0);
		return;
	}

	if (!numsol.ModifiedNewtonActive())
	{
		if (0)
		{
			NumNLSys::NLS_SetTIstages(ns);
			NumNLSys::Jacobian(m,x);
			SetJacobianComputationFlag(0);
			return;
		}
		else
		{
			if (!fulljacobianwarned) 
			{
				fulljacobianwarned=1;
				uo << "Warning: Full Jacobian turned off!!!\n";
			}
		}
	}

	//compute kversionjacobian +++
	if (KVersion())
	{
		TMStartTimer(19);
		int i;
		int j;

		int c0f = 0; //for lobattoIIIC und RadauIA methods!!!!
		if (t->c(1) == 0)
		{
			c0f = 1;
		}

		int is = GetImplicitSize();
		int es = GetFirstOrderSize();
		int sos = GetSecondOrderSize();
		int ss = GetSystemSize();
		double tau = TIstep;
		double tau2 = TIstep*TIstep;

		xx.SetLen(ss);
		xx = TIx0;

		mlocf.SetSize(es,ss);
		mlocf2.SetSize(sos,ss);

		int minit = 0;
		if (mlocm.Getcols() != sos) minit = 1;
		mlocm.SetSize(sos,sos);
		if (minit) mlocm = 0.;

		mlocg.SetSize(is,ss);
		mlocm2.SetSize(sos,sos);
		int locjacm_on = 0; //compute dM/TIkv ... normally not necessary!!!

		diag.SetSize(es,es);
		diag.FillWithZeros();
		for (i=1; i<=es; i++)
		{
			diag(i,i) = 1.;
		}

		if (es) 
		{
			LocalJacobianF(mlocf,xx);
		}
		if (sos) 
		{
			LocalJacobianF2(mlocf2,xx);
			EvalM(xx,mlocm,TItime+0.5*TIstep);

			if (locjacm_on) LocalJacobianM(mlocm2,xx);
		}
		if (is)
		{
			LocalJacobianG(mlocg,xx);
		}
		TMStopTimer(19);
		TMStartTimer(4);

		m.SetSize((ns-c0f)*(sos+es+is),(ns-c0f)*(sos+es+is));
		m.FillWithZeros(1, 1, (ns-c0f)*(sos+es+is), (ns-c0f)*(sos+es+is));


		for(i=1; i<= (ns-c0f); i++)
		{
			for(j=1; j<= (ns-c0f); j++)
			{
				int offri = sos*(i-1)+1;
				int offcj = sos*(j-1)+1;
				int offrif = sos*(ns-c0f)+es*(i-1)+1;
				int offcjf = sos*(ns-c0f)+es*(j-1)+1;
				int offrig = (sos+es)*(ns-c0f)+is*(i-1)+1;
				int offcjg = (sos+es)*(ns-c0f)+is*(j-1)+1;

				//Part A: sos
				if (i==j) m.AddSubmatrix(mlocm,1,1,offri,offcj,sos,sos,-1);
				if (locjacm_on && i==j) m.AddSubmatrix(mlocm2,1,1,offri,offcj,sos,sos,-1*tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1,         offri,offcj,sos,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1+sos,     offri,offcj,sos,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1+sos*2,   offri,offcjf,sos,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocf2,1,1+sos*2+es,offri,offcjg,sos,is,1);

				//Part B: es
				m.AddSubmatrix(mlocf,1,1,         offrif,offcj,es,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf,1,1+sos,     offrif,offcj,es,sos,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(diag,1,1,offrif,offcjf,es,es,-1);
				m.AddSubmatrix(mlocf,1,1+sos*2,   offrif,offcjf,es,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocf,1,1+sos*2+es,offrif,offcjg,es,is,1);


				//Part C: is, original
				m.AddSubmatrix(mlocg,1,1,         offrig,offcj,is,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocg,1,1+sos,     offrig,offcj,is,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocg,1,1+sos*2,   offrig,offcjf,is,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocg,1,1+sos*2+es,offrig,offcjg,is,is,1);
/*
				//Part A: sos
				if (i==j) m.AddSubmatrix(mlocm,1,1,offri,offcj,sos,sos,-1);
				if (locjacm_on && i==j) m.AddSubmatrix(mlocm2,1,1,offri,offcj,sos,sos,-1*tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1,         offri,offcj,sos,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1+sos,     offri,offcj,sos,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocf2,1,1+sos*2,   offri,offcjf,sos,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocf2,1,1+sos*2+es,offri,offcjg,sos,is,1);

				//Part C: is
				m.AddSubmatrix(mlocg,1,1,         offrig,offcj,is,sos,tau2*t->A2(i+c0f,j+c0f));
				m.AddSubmatrix(mlocg,1,1+sos,     offrig,offcj,is,sos,tau*t->A(i+c0f,j+c0f));
				m.AddSubmatrix(mlocg,1,1+sos*2,   offrig,offcjf,is,es,tau*t->A(i+c0f,j+c0f));
				if (i==j) m.AddSubmatrix(mlocg,1,1+sos*2+es,offrig,offcjg,is,is,1);
*/
			}
		}

		TMStopTimer(4);


		//if (!system_matrices_written && solset.write_MK2file /*&& RelApproxi(GetTime(), 0.320)*/) //$ PG 2012-2-3: modified (see next line)
		if (!system_matrices_written && GetOptions()->LoggingOptions()->WriteMassAndStiffnessMatrix() /*&& RelApproxi(GetTime(), 0.320)*/)
		{
			system_matrices_written= 1;

			mystr dir = solset.sol_directory;
			mystr fileload = "loadMK.m";
			mystr fileK = "Kdat.dat";
			mystr fileM = "Mdat.dat";
			mystr fileD = "Ddat.dat";

			ofstream fout((dir+fileload).c_str());
			ofstream Kdat((dir+fileK).c_str());
			ofstream Mdat((dir+fileM).c_str());
			ofstream Ddat((dir+fileD).c_str());

			int dim = mlocm.Getcols();
			Matrix KK; KK.CopyFrom(mlocf2,1,1,dim,dim);

			//int dim = m.Getcols(); //full jacobian
			//Matrix KK; KK.CopyFrom(m,1,1,dim,dim);


			SparseMatrix KS(KK);
			KS.PrintToMatlabSparse(Kdat);
			fout << "load Kdat.dat\nK = spconvert(Kdat);\n";

			Matrix DD; DD.CopyFrom(mlocf2,1,1+dim,dim,dim*2);
			SparseMatrix DS(DD);
			DS.PrintToMatlabSparse(Ddat);
			fout << "load Ddat.dat\nD = spconvert(Ddat);\n";

			SparseMatrix MS(mlocm);
			MS.PrintToMatlabSparse(Mdat);
			fout << "load Mdat.dat\nM = spconvert(Mdat);\n";
			fout << "% use command FULL to show complete matrices\n";

			fout << flush;

			UO() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			UO() << "M and K written to file '" << (dir+fileM) << "'\nand '" << (dir+fileK) << "'!\n";
			UO() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		}


		/*
		uo << "Jac=" << m << "\n";
		uo << "K=" << mlocf2 << "\n";
		uo << "M=" << mlocm << "\n";
		PrintMatrix01(m); */
		/*
		Matrix mi = m;
		mi.Invert();
		uo << "minv=" << mi << "\n";
		uo << "m*minv=" << (m*mi) << "\n";
		*/
		/*
		uo << "TIx0=" << TIx0 << "\n";
		//PrintMatrix01(m); 
		uo << m << "\n";
		*/

		//if (RelApproxi(TItime, 0.0) || RelApproxi(TItime, 0.320))
		{
			//uo << "Jac=" << m << "\n";
			//uo << "K=" << mlocf2 << "\n";
			//uo << "M=" << mlocm << "\n";
		}

		/*
		Matrix minv = m;
		uo << "Jac-invertible=" << minv.Invert2() << "\n";
		uo << "M^-1=" << minv << "\n";*/




	}
	else if (GetSecondOrderSize()==0)
	{
		IRK_Tableau* t = tableaus[act_tableau];

		int i;
		int j;

		int is = GetImplicitSize();
		int es = GetFirstOrderSize();
		int ss = GetSystemSize();


		mlocf.SetSize(es,ss);
		mlocg.SetSize(is,ss);
		mlocf2.SetSize(es,ss);

		diag.SetSize(es,es);
		diag.FillWithZeros();
		for (i=1; i<=es; i++)
		{
			diag(i,i) = -1.;
		}

		LocalJacobianF(mlocf,x);
		LocalJacobianG(mlocg,x);

		m.SetSize(ns*(es+is),ns*(es+is));

		for(i=1; i<= ns; i++)
		{
			for(j=1; j<= ns; j++)
			{
				mlocf2 = mlocf;
				mlocf2 *= (TIstep * t->A(i,j));

				if (i==j)
				{
					mlocf2.AddMatrix(1,1,diag);
					m.InsertMatrix(es*ns+1+(i-1)*is,1+(j-1)*es,mlocg,1,1,is,es);
					m.InsertMatrix(es*ns+1+(i-1)*is,1+es*ns+(j-1)*is,mlocg,1,es+1,is,is);
				}
				else
				{
					m.FillWithZeros(es*ns+1+(i-1)*is,1+(j-1)*es,is,es);
					m.FillWithZeros(es*ns+1+(i-1)*is,1+es*ns+(j-1)*is,is,is);
				}
				m.InsertMatrix(1+(i-1)*es,1+(j-1)*es,mlocf2,1,1,es,es);
				m.InsertMatrix(1+(i-1)*es,1+es*ns+(j-1)*is,mlocf2,1,1+es,es,is);
			}
		}
	}
	else
	{
		NumNLSys::Jacobian(m,x);
		SetJacobianComputationFlag(0);
		return;
	}

	//uo << "Jac=\n"; PrintMatrix01(m);

	SetJacobianComputationFlag(0);
}

void TimeInt::NLF(const Vector& x, Vector& f)
{
	TMStartTimer(13);
	int ss = GetSystemSize();
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();

	if (DoStaticComputation())
	{
		//$$$
		evalf.SetLen(sos);
		evalfe.SetLen(es);
		evalg.SetLen(is);
		xh.SetLen(2*sos+es+is);

		//use vector xh - help vector
		//copy components and set velocities 0
		xh.Copy(x,1,1,sos);
		xh.Copy(x,1+sos,1+2*sos,es);
		xh.Copy(x,1+sos+es,1+2*sos+es,is);
		for (int i = 1+sos; i <= 2*sos; i++) xh(i) = 0;
		//uo << "xh=" << xh << "\n";
		EvalF2(xh, evalf,   TItime);
		EvalF(xh,evalfe, TItime);
		EvalG(xh,evalg,  TItime);
		//uo << "evalf=" << evalf << "\n";

		//does not work, spring-type regularisation:
		if (solset.static_spring_type_reg_param != 0.)
		{ 
			for (int i=1; i <= sos; i++) evalf(i) -= LoadFact()*solset.static_spring_type_reg_param*x(i);
		}
		f.SetLen(sos+es+is);
		f.Copy(evalf, 1,1,sos);
		f.Copy(evalfe,1,1+sos,es);
		f.Copy(evalg, 1,1+sos+es,is);

		//uo << "NLF=" << f << "\n";
		TMStopTimer(13);
		return;
	}


	IRK_Tableau* t = tableaus[act_tableau];
	int ns = t->nstages;

	if (!impl_equ_flag)
	{

		int i, j;

		//Kversion of Runge Kutta schemes:
		if (TIkv)
		{
			evalf.SetLen(sos);
			evalfe.SetLen(es);
			evalg.SetLen(is);
			int ns = t->nstages;
			int c0f = 0; //for lobattoIIIC und RadauIA methods!!!!
			if (t->c(1) == 0)
			{
				c0f = 1;
			}

			if (c0f)
			{
				ku[1] = TIku0;
				kv[1] = TIkv0;
				kue[1] = TIku0e;
			}

			//construct kv[i]
			for (i = 1; i <= (ns-c0f); i++)
			{
				kv[i+c0f].SetLen(sos);
				kv[i+c0f].Copy(x,1+sos*(i-1),1,sos);
			}
			//construct kue[i]
			for (i = 1; i <= (ns-c0f); i++)
			{
				kue[i+c0f].SetLen(es);
				kue[i+c0f].Copy(x,1+sos*(ns-c0f)+es*(i-1),1,es);
			}

			xh.SetLen(sos);
			//compute ku[i]
			for (i = 1; i <= (ns-c0f); i++)
			{
				ku[i+c0f] = TIv0;
				for (j = 1; j <= ns; j++)
				{
					/*
					xh = kv[j];
					xh *= TIstep * t->A(i+c0f,j);
					ku[i+c0f] += xh;*/
					ku[i+c0f].MultAdd(TIstep * t->A(i+c0f,j), kv[j]);
				}
			}
			//TIu0, TIv0, TIu0e passed from GeneralIStep

			//++++++++++++++++++++++++++++++++++++++++++++++++
			//construct xi[i]:

			//Compute gu[i]
			for (i = 1; i <= ns; i++)
			{
				gu[i] = TIu0;
				for (j = 1; j <= ns; j++)
				{
					/*xh = ku[j];
					xh *= TIstep * t->A(i,j);
					gu[i] += xh;*/
					gu[i].MultAdd(TIstep * t->A(i,j),ku[j]);
				}
			}

			//Compute gv[i]
			for (i = 1; i <= ns; i++)
			{
				gv[i] = TIv0;
				for (j = 1; j <= ns; j++)
				{
					/*xh = kv[j];
					xh *= TIstep * t->A(i,j);
					gv[i] += xh;*/
					gv[i].MultAdd(TIstep * t->A(i,j),kv[j]);
				}
			}

			//Compute gue[i]
			if (es!=0)
			{
				xh2.SetLen(es);
				for (i = 1; i <= ns; i++)
				{
					gue[i] = TIu0e;
					for (j = 1; j <= ns; j++)
					{
						/*xh2 = kue[j];
						xh2 *= TIstep * t->A(i,j);
						gue[i] += xh2;*/
						gue[i].MultAdd(TIstep * t->A(i,j),kue[j]);
					}
				}
			}
			//construct iv[i]
			for (i = 1; i <= (ns-c0f); i++)
			{
				iv[i].SetLen(is);
				iv[i].Copy(x,(ns-c0f)*(sos+es)+1+is*(i-1),1,is);
			}

			for (i = 1; i <= (ns-c0f); i++)
			{
				xi[i].SetLen(2*sos+es+is);
				xi[i].Copy(gu[i+c0f],1,1,sos);
				xi[i].Copy(gv[i+c0f],1,1+sos,sos);
				xi[i].Copy(gue[i+c0f],1,1+2*sos,es);
				xi[i].Copy(iv[i],1,1+2*sos+es,is);
			}
			//++++++++++++++++++++++++++++++++++++++++++++++++

			f.SetLen((sos+is+es)*(ns-c0f));
			if (sos != 0)
			{
				//Compute -M_i K_i + F(u_i,v_i)
				evalf.SetLen(sos);
				xh.SetLen(sos);
				for (i = 1; i <= (ns-c0f); i++)
				{
					PrecomputeEvalFunctions(); //$ PG 2012-1-15: precomputation for EvalF(), EvalF2(), EvalG(), EvalM(), and also for CqTLambda() in case of constraints //BUG: at the moment only available for EvalF2() and EvalM() //Has to be fixed via reorganizing loops in NLF(..)
					if (!UseSparseMass())
					{
						TIm.SetSize(sos,sos);
						if (!TIm_initialized)
						{
							TIm.FillWithZeros();
							TIm_initialized = 1;
							EvalM(xi[i], TIm, TItime+TIstep*t->c(i+c0f));

							if (TransformJacApply()) 
							{
								//uo << "***********************\nERROR: Buggy Mass matrix used!!!!!!!!!!\n***********************\n";
								uo << "***********************\nWARNING: Modified Mass matrix used!!!!!!!!!!\n***********************\n";
							}
						}
						if (!TransformJacApply() && !solset.assume_constant_mass_matrix) 
						{
							EvalM(xi[i], TIm, TItime+TIstep*t->c(i+c0f));
						}
					}
					else
					{
						//use sparse mass matrix:
						if (!TIm_initialized || !TransformJacApply())
						{
							TIm_initialized = 1;
							TIm_sparse.SetSize(sos,sos,MaxSparseBandwidth());
							TIm_sparse.FillWithZeros();
							EvalM(xi[i], TIm_sparse, TItime+TIstep*t->c(i+c0f));
							if (TransformJacApply()) 
								uo << "***********************\nWARNING: Modified Mass matrix used!!!!!!!!!!\n***********************\n";
							//uo << "***************\nERROR: Buggy Mass matrix used!!!!!!!!!!\n***************\n";
						}
						else if (!TransformJacApply() && !solset.assume_constant_mass_matrix)
						{
							//sometimes not needed!!!
							TIm_sparse.FillWithZeros();
							EvalM(xi[i], TIm_sparse, TItime+TIstep*t->c(i+c0f));
						}
					}
					EvalF2(xi[i], evalf, TItime+TIstep*t->c(i+c0f));

					//for (int k=1; k<=20; k++)
					if (!UseSparseMass())
					{ 
						if (bandwidthm)
						{
							MultBW(TIm,kv[i+c0f],xh,bandwidthm);
						}
						else
						{
							Mult(TIm,kv[i+c0f],xh);
						}
					}
					else
					{
						Mult(TIm_sparse,kv[i+c0f],xh);
					}
					evalf -= xh; //evalf = -M_i K_i + F(u_i,v_i)
					f.Copy(evalf,1,1+(i-1)*sos,sos);
				}
			}
			if (es != 0)
			{
				evalfe.SetLen(es);
				for (i = 1; i <= (ns-c0f); i++)
				{
					EvalF(xi[i],evalfe,TItime+TIstep*t->c(i+c0f));
					evalfe-=kue[i+c0f];
					f.Copy(evalfe,1,1+(ns-c0f)*sos+es*(i-1),es);
				}
			}
			if (is != 0)
			{
				for (i = 1; i <= (ns-c0f); i++)
				{
					EvalG(xi[i],evalg,TItime+TIstep*t->c(i+c0f));
					f.Copy(evalg,1,1+(ns-c0f)*(sos+es)+is*(i-1),is);
				}
			}
		}
		//pure first order equations (with constraints):
		else if (GetSecondOrderSize() == 0)
		{

			evalg.SetLen(is);
			evalf.SetLen(es);

			for (i = 1; i <= ns; i++)
			{
				gv[i].SetLen(es);
				gv[i].Copy(x,1+es*(i-1),1,es);

				iv[i].SetLen(is);
				iv[i].Copy(x,1+ns*es+is*(i-1),1,is);

				xv[i].SetLen(es+is);
				xv[i].Copy(x,1+es*(i-1),1,es);
				xv[i].Copy(x,1+ns*es+is*(i-1),1+es,is);
			}

			TIu0.SetLen(es);
			TIu0.Copy(TIx0,1,1,es);

			f.SetLen(ns*(es+is));

			for (j = 1; j <= ns; j++)
			{
				evalfs[j].SetLen(es);
				EvalF(xv[j], evalfs[j], TItime+TIstep*t->c(j));
			}

			for (i = 1; i <= ns; i++)
			{
				fg = TIu0;
				fg -= gv[i];
				for (j = 1; j <= ns; j++)
				{
					//EvalF(xv[j], evalf, TItime+TIstep*t->c[j]);
					evalf = evalfs[j];
					evalf.Mult(TIstep * t->A(i,j));
					fg += evalf;
				}
				f.Copy(fg,1,1+(i-1)*es,es);
			}

			if (is != 0)
			{
				for (i = 1; i <= ns; i++)
				{
					EvalG(xv[i],evalg,TItime+TIstep*t->c(i));
					f.Copy(evalg,1,1+ns*es+(i-1)*is,is);
				}
			}
		}
	}
	else
	{
		EvalG((TIx0.SubVector(1,es)).Append(x),f,TItime+TIstep);
	}
	TMStopTimer(13);
}

//convert seconds to format days hours:min:secs
mystr GetTString(double s)
{
	int days = (int)(s/(3600.*24.));
	int hours = (int)((s-days*(3600*24))/3600.);
	int mins = (int)((s-days*(3600*24)-hours*3600)/60.);
	int secs = (int)((s-days*(3600*24)-hours*3600-mins*60));
	int millis = (int)((s-(double)((int)s))*10.);

	mystr ts="";
	if (days) {ts+=mystr(days)+mystr(" days, ");};
	if (hours || days) 
	{
		ts+=mystr(hours)+mystr(":");
		if (mins < 10) {ts+="0";}
	};
	ts+=mystr(mins)+mystr(":");
	if (secs < 10) {ts+="0";}
	ts+=mystr(secs)+mystr(".");
	ts+=mystr(millis);

	return ts;
}

double TimeInt::GetClockTime() const
{
	timeb tb;
	tb.time = 0;
	tb.millitm = 0;
	ftime(&tb);
	//cout << "time=" << tb.time << ", millisec=" << tb.millitm << endl;
	return tb.time+(double)tb.millitm*0.001;
}

void TimeInt::LoadTableaus(const mystr& tableauname) 
{


	//----------------------------------- read tableau
	CMFile tin(tableauname, TFMread);

	int ntab, i, j, k;

	mystr method;
	IRK_Tableau* tableau;
	tableaus.SetLen(0);

#ifdef COMPILE_AND
#include "../ANDlib/LoadTableaus_And.h"
#else

	tin.RWInt(ntab); // #tableaus
	if(ntab < 1 || ntab > 10000)
	{
		uo.pUI->InstantMessageText("Could not load the tables from \"" + tableauname + "\"!! Check if working directory of WCDriver is set.");
		return;
	}

	i = 1;
	while(i <= ntab)
	{

		tin.RWuntilEOL(method,0);

		if (method == mystr("#Runge-Kutta"))
		{
			tableau = new IRK_Tableau();

			tableaus.Add(tableau);

			i++;

			tin.RWuntilEOL(tableau->name,0);
			tin.RWuntilEOL(tableau->info,0);
			tin.RWInt(tableau->ODE_order);
			//uo << "load tableau '" << tableau->name << "', order=" << tableau->ODE_order << ENDL;

			tin.RWInt(tableau->DAE_order_y[1]);
			tin.RWInt(tableau->DAE_order_y[2]);
			tin.RWInt(tableau->DAE_order_z[1]);
			tin.RWInt(tableau->DAE_order_z[2]);
			tin.RWInt(tableau->nstages);

			int ns=tableau->nstages;
			tableau->b = Vector(ns);
			for (j = 1; j <= ns; j++)
			{
				tin.RWDouble(tableau->b(j));
			}

			tableau->c = Vector(ns);
			for (j = 1; j <= ns; j++)
			{
				tin.RWDouble(tableau->c(j));	      
			}

			tableau->A = Matrix(ns, ns);
			for (j = 1; j <= ns; j++)
			{
				for (k = 1; k <= ns; k++)
				{
					tin.RWDouble(tableau->A(j,k));
				}
			}

			//second order equations:
			tableau->A2 = tableau->A*tableau->A;
			//for methods which do not fulfill c1=1:
			tableau->Ainv = tableau->A;
			tableau->Ainvertable = tableau->Ainv.Invert();
			if (tableau->Ainvertable)
			{
				tableau->bAinv = tableau->b*tableau->Ainv;
			}

			//check if tableau is implicit
			tableau->implicit = 0;
			if (tableau->c(1) != 0) {tableau->implicit = 1;}
			for (j = 1; j <= ns; j++)
			{
				for (k = j; k <= ns; k++)
				{
					if (tableau->A(j,k) != 0) {tableau->implicit = 1;} 
				}
			}
			/*			
			if (!tableau->implicit) 
			{
			uo << tableau->info << "\n" 
			<< "A=" << tableau->A << "\n" 
			<< "b=" << tableau->b << "\n" 
			<< "c=" << tableau->c << "\n"; 
			}*/
		}
	}
#endif
	uo << "loaded " << ntab << " Runge Kutta Tableaus\n";
};


int TimeInt::GetTableau(mystr tabname, int stages)
{

	int i=1;
	while (i<=tableaus.Length() && (tableaus[i]->name != tabname || tableaus[i]->nstages != stages))
	{
		i++;
	}

	if (i>tableaus.Length())
	{
		uo << "ERROR: Could not find an according tableau '" << tabname.c_str() 
			<< "'\n   to the number of stages" << stages << " !" << ENDL;
		return 0;
	}
	else
	{
		return i;
	}
}

int TimeInt::LoadTableauList(mystr tabname, TArray<int>& list, int minstage, int maxstage)
{
	int i;
	int found = 0;
	for (i = 1; i <= tableaus.Length(); i++)
	{
		if (tableaus(i)->name == tabname &&
			tableaus(i)->nstages >= minstage && tableaus(i)->nstages <= maxstage) 
		{
			list.Add(i);
			found = 1;
		} 
	}
	return found;
}

double TimeInt::LagrangeWeight(int l, double tau, const Vector& c) const
{
	double w=1;
	for (int i=1; i <= c.Length(); i++)
	{
		if (i!=l)
		{
			w *= (tau-c(i))/(c(l)-c(i));
		}
	}
	return w;
}

void TimeInt::PrintTimingList()
{
	TMPrintTimer();
}


void TimeInt::ResetComputation()
{
	TMResetTimer();
	Initialize();
	TIInit();
	//uo.pUI->CallWCDriverFunction(1); //redraw 
}



//print step information and draw results each step
void TimeInt::StepPrintAndDraw(int rv)
{
	//calc time:
	double lt = GetClockTime();
	TIcomptime = TIV.start_clock_time + lt;

	//do graphics if it's time for it!
	solset.withgraphics = GetIOption(104);
	//redraw frequency: 0..off, 1..draw last frame, 2..100sec, 3..20sec, 4..2sec, 5..200ms, 6..50ms, 7..20ms, 8..every 10 frames, 9..every frame
	if (solset.withgraphics == 1 && TItime > solset.endtime-0.1*TIminstep) {TIV.lastdraw = TItime; drawnow++;}
	if (solset.withgraphics == 2 && (TIcomptime - TIV.lastdraw) > 100.) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 3 && (TIcomptime - TIV.lastdraw) >  20.) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 4 && (TIcomptime - TIV.lastdraw) >   2.) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 5 && (TIcomptime - TIV.lastdraw) > 0.2) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 6 && (TIcomptime - TIV.lastdraw) > 0.05) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 7 && (TIcomptime - TIV.lastdraw) > 0.02) {TIV.lastdraw = TIcomptime; drawnow++;}
	if (solset.withgraphics == 8 && (TIit - TIV.lastdraw) >= 10) {TIV.lastdraw = TIit; drawnow++;}

	TIV.timetogo2 = TIV.timetogo;
	TIV.timetogo = TIcomptime / TItime * (solset.endtime-TItime);
	TIV.timetogo = 0.2*TIV.timetogo+0.8*TIV.timetogo2;

	TIV.maxnewtonit = Maximum(GetNewtonIt(),TIV.maxnewtonit);

	if (TItime-solset.endtime >= -0.1*TIminstep)
	{
		drawnow++;
		TIDrawAndStore();
	}
	if ((TIcomptime - TIV.last_showstatustext) >= 1 || //update text every GetDOption(  4) seconds
		((TItime-solset.endtime)>=-0.1*TIminstep) ||
		(rv == 0 || StopCalculation() || TIfinished!=0)) 
	{
		char progress_str[64];
		double time_total = solset.endtime - solset.starttime;
		double time_left = solset.endtime - TItime;

		double progress = 100;
		if (time_total != 0) progress = 100.*(time_total-time_left)/time_total;

		char pc = '%';
		sprintf(progress_str, "%2.1f%c", progress, pc);
		uo.pUI->StatusText(progress_str);
	
#define test_progressbar_not  // hack (AD) test ProgressBar
#ifdef test_progressbar
		this->UO().pUI->CallWCDriverFunction(10,2,100);                // 100 ticks
		this->UO().pUI->CallWCDriverFunction(10,3,(int)progress+0.5);  // actual position to %value
    
		ElementDataContainer edc;																			 // text fields
		ElementData ed;
		
		ed.SetText("Progress bar");
		ed.SetDataName("Caption");
		edc.Add(ed);

		ed.SetText("Subroutine: StepPrintAndDraw");
		ed.SetDataName("Over");
		edc.Add(ed);

		ed.SetText(mystr(progress) + mystr("% progress (estimate)")); 
		ed.SetDataName("Under");
		edc.Add(ed);
		
		this->UO().pUI->CallWCDriverFunction(10,4,0,&edc);
#endif

		TIV.last_showstatustext = TIcomptime;
	}

	//print results, if it's time for it:
	if (SolverPrintsDetailedOutput() ||
		(TIcomptime - TIV.lastprintres) >= GetOptions()->LoggingOptions()->ComputationOutputEveryXSec() || //update text every GetOptions()->LoggingOptions()->ComputationOutputEveryXSec() seconds
		((TItime-solset.endtime)>=-0.1*TIminstep) ||
		(rv == 0 || StopCalculation() || TIfinished!=0)) 
	{
		//char progress_str[64];
		double time_total = solset.endtime - solset.starttime;
		double time_left = solset.endtime - TItime;

//		#ifdef gencnt
//		uo << "number of new calls for Vector/Matrix: n_gen_vec=" << *genvec << ", n_gen_mat=" << *genmat << "\n";
//		#endif

		//double progress = 100;
		//if (time_total != 0) progress = 100.*(time_total-time_left)/time_total;

		//char pc = '%';
		//sprintf(progress_str, "%2.1f%c", progress, pc);
		//uo.pUI->StatusText(progress_str);

		TIV.lastprintres = TIcomptime;

		if (TIV.outputlevel >=3)
		{
			//uo << "Generated vectors = " << GenVecCnt() << ", matrices = " << GenMatCnt() << ENDL;

			SolverUO() << 
				("step=")+mystr(TIit)+
				(", t=")+mystr(TItime)+
				(", tinc=")+mystr(TIstep)+
				//mystr(", x1=")+mystr(GetCharacteristicSol())+
				(", Nits=")+mystr(TIV.maxnewtonit)+
				(", Dits=")+mystr(log.TInonlinit)+
				(", Jacs=")+mystr(GetJacCount())+
				(", real_time=")+mystr(GetTString(TIcomptime))+
				(", real_time_togo=")+mystr(GetTString(TIV.timetogo))+
				("\n");
			uo.ResetLocalMessageLevel();
		}
		//$ YV 2013-01-03: the code below was anyway inactive - this needs to be done differently in the future
#if 0
		if (0 && contact_cnt_mastersearch!=0) //write contact information
		{
			uo << "n_gap=" << contact_cnt_gap << ", n_ms=" << contact_cnt_mastersearch;
      uo << ", itemlist=" << (double)contact_cnt_itemlist/(double)contact_cnt_mastersearch;
      uo << ", foundmaster=" << (double)contact_cnt_foundmaster/(double)contact_cnt_mastersearch;
      uo << ", list2=" << (double)contact_cnt_list2/(double)contact_cnt_mastersearch;
      uo << ", list3=" << (double)contact_cnt_list3/(double)contact_cnt_mastersearch;
      uo << ", list4=" << (double)contact_cnt_list4/(double)contact_cnt_mastersearch;
			uo << "\n";
		}
#endif
		TIV.maxnewtonit = 0;
	}
}

int TimeInt::Integrate() //JG const SolverSettings& solversettings)
{

	SetSolver(&numsol); //in order to be able to use solver option right away

	//solset = solversettings; //JG, changed because timeint and MBS now have same solversettings

	SetComputationSolverOptions(); //copy general options to MBSSolSet(), because MBSSolSet() is faster!
	log.Init(); //initialize counters
#ifdef gencnt
	(*genvec) = 0;
	(*genmat) = 0;
#endif


	//+++++++++
	//==>put into log!!!!
	//SetNewtonIt(0);
	//SetNewtonItSum(0);
	//jaccount = 0;
	//changestep=0;
  //++++++++++


	TItime = solset.starttime;
	TISaveState(TIlaststep_state); //save state of end of last step (or initial step) ==> last step available


	if (lastloadresults > TItime) lastloadresults = -1e100;
	if (laststoredata > TItime)	laststoredata = -1e100;

	//uo << "storedata=" << solset.storedata << "\n";

	if (solset.maxstages >= global_maxstages)
	{
		uo.pUI->InstantMessageText("Maximum number of stages (20) is succeeded!!!");
	}

	//convert to old parameters:
	double remainigtime = solset.endtime-TItime;

	double steps = floor(remainigtime/solset.maxstepsize+0.1*solset.minstepsize);
	TIbasestep = remainigtime/(double)steps;
	TIstep = TIbasestep;

	if (solset.initstepsize<=solset.maxstepsize && solset.initstepsize> solset.minstepsize && FullAdaptive())
	{
		TIstep = solset.initstepsize;
	}

	system_matrices_written = 0; //flag, system matrices are written at first time step

	//numsol.DestroyOldJacs();
	TIm_initialized = 0;
	bandwidthm = 0;

	TImaxstep = solset.maxstepsize;
	TIminstep = solset.minstepsize;
	TIstepnew = TIstep;

	TInoincs = 1;
	TInonls = 1;

	int rv = 1;

	TIfinished = 0;

	TIit = 0;
	TIrejectedsteps = 0;
	pCFB->ResultsUpdated(0); //initial values are also stored

	//initialize timer variables:
	double lt = GetClockTime(); //temporary variable, store actual clock
	TIV.timetogo = 0;
	TIV.lastdraw = 0;
	TIV.lastprintres = 0;
	TIV.outputlevel = 3; //2
	TIV.maxnewtonit=0;

	//time(&lt);
	TIV.start_clock_time=-lt;

	steps_since_orderchange = 0;
	min_newtonits = NLS_MaxNewtonSteps();
	newjacobians_since_orderchange;
	last_jaccount = 0;
	reduce_order = 0;
	
	numsol.ModifiedNewton() = NLS_ModifiedNewton();
	numsol.AbsoluteAccuracy() = NLS_AbsoluteAccuracy();
	numsol.RelativeAccuracy() = NLS_RelativeAccuracy()/sqrt(TIstepnew);
	numsol.NumDiffepsi() = NLS_NumDiffepsi();
	numsol.MaxModNewtonSteps() = NLS_MaxModNewtonSteps();
	numsol.MaxRestartNewtonSteps() = NLS_MaxRestartNewtonSteps();
	numsol.MaxFullNewtonSteps() = NLS_MaxFullNewtonSteps();
	numsol.TrustRegion() = NLS_TrustRegion();
	numsol.TrustRegionDiv() = NLS_TrustRegionDiv();
	numsol.SymmetricJacobian() = NLS_SymmetricJacobian();
	numsol.UseSparseSolver() = NLS_UseSparseSolver();
	numsol.SolveUndeterminedSystem() = NLS_SolveUndeterminedSystem();
	numsol.EstimatedConditionNumber() = NLS_EstimatedConditionNumber();

	uo << "Use sparsesolver=" << UseSparseSolver() << "\n";

	if (!UseSparseMass())
	{
		mlocm.SetSize(GetSecondOrderSize(),GetSecondOrderSize());
		mlocm.FillWithZeros();
 	}

	if (DoStaticComputation())
	{
		//write initial state:
		FirstStep();
		rv = StaticComputation();

		drawnow = 1;
		TIDrawAndStore();

		//write solution:
		FirstStep();
	}
	else if (!DoImplicitIntegration())
	{
		//write initial state:
		FirstStep();
		rv = ExplicitIntegration();

		drawnow = 1;
		TIDrawAndStore();

		//write solution:
		FirstStep();
	}
	else
	{
		drawnow = 0;
		if (solset.variable_order)
		{
			stage_tableaus.Flush();
			int success = LoadTableauList(solset.tableauname, stage_tableaus, solset.minstages, solset.maxstages);

			if (!success) {uo << "Could not find tableaus '" << solset.tableauname << "'!!!" << ENDL; return 0;}

			act_stage = 1;
			act_tableau = stage_tableaus(act_stage);
		}
		else
		{
			act_tableau=GetTableau(solset.tableauname, solset.maxstages);
		}

		if (!act_tableau) return 0;

		if (TIV.outputlevel > 2)
			uo << "Integration: " << tableaus[act_tableau]->name << ", " 
			<< tableaus[act_tableau]->info << ", ODE-order=" << tableaus[act_tableau]->ODE_order << ENDL;

		//write solution:
		FirstStep();

		while (((solset.endtime-TItime)>0.1*TIminstep) && rv != 0 && !StopCalculation() && TIfinished==0)
		{
			if (SolverPrintsDetailedOutput())
			{
				mystr str("**************** step ");
				str+=mystr(TIit+1);
				str+=mystr(" (implicit");
				if (FullAdaptive()) str+=mystr(", adaptive");
				str+=mystr(") begin\n");
				SolverUO() << str;
			}

			if(!FullAdaptive() && RelApproxi(solset.endtime,TItime+TIstep))
			{
				TIstep = solset.endtime - TItime;
			}


			rv = IntegrateStep();
			TIit++;
			//uo << "Veccnt=" << GenVecCnt() << "\n";
			//uo << "Matcnt=" << GenMatCnt() << "\n";

			StepPrintAndDraw(rv);

			if (SolverPrintsDetailedOutput())
			{
				mystr str("**************** step ");
				str+=mystr(TIit);
				str+=mystr(" (implicit");
				if (FullAdaptive()) str+=mystr(", adaptive");
				str+=mystr(") end\n");
				SolverUO() << str;
			}
		}
	}


	//set new start time = actual end time (or breaking point)
	MBS_EDC_TreeSetDouble(TItime, "SolverOptions.start_time"); //set endtime as start time--> if start clicked several times
	//==> also reset initial conditions??


	if (!rv && NumSolver().GetErrorMsg().Length()!=0) 
	{uo << "step=" << TIit-1 << ", ERROR in Newton method: " << NumSolver().GetErrorMsg().c_str() << ENDL;}
	else if (!rv && !DoStaticComputation()) {uo << "ERROR: Time integration not successful!!!" << ENDL;}
	else if (!rv && DoStaticComputation())  {uo << "ERROR: Static computation not successful!!!" << ENDL;}

	if (TIfinished == -1) {uo << "Time integration stopped with error." << ENDL; rv = 0; }

	return rv;
}


int TimeInt::IntegrateStep()
{
	TIx0start = TIx0;
	TIk0start = TIk0;
	int rv;

	if(!FullAdaptive())
	{
		rv = FullStep();
	}
	else
	{
		rv = FullAdaptiveStep();
	}

	return rv;
}

int TimeInt::StaticComputation()
{
	int rv = 1; //return value
	int oldmodnewton = NLS_ModifiedNewton(); // returns solver setting nls_modifiednewton
	if (!numsol.UseSparseSolver()) NLS_ModifiedNewton() = 1; // activates modified newton method if sparse solver is not used

	//--------------------------------------------------------------
	// read system size, set initial values
	mystatic Vector GSxs;
	int ss  = GetSystemSize();
	int sos = GetSecondOrderSize();
	int es  = GetFirstOrderSize();
	int is  = GetImplicitSize();

	GSxs.SetLen(sos+es+is);

	//TIx0 contains initial values
	for (int i=1; i <= sos; i++) GSxs(i) = TIx0(i);
	for (int i=1; i <= es;  i++) GSxs(i+sos) = TIx0(i+2*sos);
	for (int i=1; i <= is;  i++) GSxs(i+sos+es) = TIx0(i+2*sos+es);

	//--------------------------------------------------------------
	//tell solver that it is a static step, time-step -1 will not happen
	numsol.NLSolveInfo() = -1; 

	//--------------------------------------------------------------
	// declare and/or initiate counters and incrememts
	loadfact = 0.;
	int successcnt = 0;
	int failcnt = 0;
	int newtonits = 0;
	int jaccnt = 0;
	int nlit = 0;
	int lastnewtonits = 0;
	int lastjaccnt = 0;

	double lt = GetClockTime();
	double start_clock_time=-lt;
	double lastdraw = 0;
	double lastprintres = 0;
	drawnow = 0;
	TIit = 0;

#define computationsteps //$ AD 2011-08: ComputationSteps replace LoadSteps in the long run...

#ifndef computationsteps // must be set AFTER ApplyComputationStepSettings for correct values
	solset.static_initloadinc = Minimum(solset.static_initloadinc,solset.static_maxloadinc);
	double loadinc = solset.static_initloadinc;
#endif

#ifdef computationsteps
	ParseStepEndTimes();																	// extract all computation step end times from EDC
	DoubleCheckWithLoadSettings();
  double cslength = GetTMaxCompSteps();                 // last step end time
	solset.endtime=max(cslength,solset.endtime);          // could be forced to 
	ComputationStepNumber() = 1;                          // begin computation in first step
	ApplyComputationStepSettings(ComputationStepNumber());
// get correct load inc for new ComputationStep, DO NOT CHANGE value of solset.static_initloadinc any more
	double loadinc = Minimum(solset.static_initloadinc,solset.static_maxloadinc);
#else
	//--------------------------------------------------------------
	//vary load steps from starttime to endtime(or maximal load steps):
	TItime = solset.starttime; // set initial artificial time to starttime (default = 0)
	double lslength = GetTMaxLoadStepsLength();
	solset.endtime=max(lslength, solset.endtime); // endtime default=1; endtime otherwise maximal number of loadsteps in example
	//$ YV 2011-06:	according to the change in GetTMaxLoadStepsLength(), added the following test to avoid situations when endtime is nowhere specified

#endif
	if(solset.endtime <= 0)
		solset.endtime = 1;
	UO(UO_LVL_err) << "static computation, artificial time interval [" << TItime << ", " << solset.endtime << "]\n";
	logout.flush();

	//--------------------------------------------------------------
	// while
	//		endtime not reached
	//		return value is 1
	//		calculation is not stopped
	//		time integration is not finished
	// do 
	//		Load iteration (load corresponds to artificial time)
	while (TItime < solset.endtime-1e-12 && rv == 1  && !StopCalculation() && TIfinished==0)
	{
		if (SolverPrintsDetailedOutput())
		{
			SolverUO() << mystr("**************** step ")+mystr(TIit+1)+mystr(" (static, adaptive) begin\n");
		}

		//+++++++++++++++++++++++++++++
		// Time of last converged step, in case this step does not converge
		double TItime_last_converged_step = TItime;

		//+++++++++++++++++++++++++++++
		// add load increment to TItime

#ifdef computationsteps
		double thisstependtime = GetEndTimeCompStep(ComputationStepNumber()); 
		// if next timevalue is higher or marginally smaller then step end time, use step end time
		if( TItime + loadinc + 1e-12 - thisstependtime > 0.)
			TItime = thisstependtime;
		else
			TItime += loadinc;
#else
		// use next integer of TItime if the difference is already smaller than 1e-12
		if (floor(TItime + loadinc + 1e-12) != floor(TItime))
			loadinc = 1.-(TItime-floor(TItime));
		// utilize all integers TItime passes and do not miss them out due to a higher load increment
		if ( floor(TItime + loadinc) != floor(TItime))
			TItime=floor(TItime + loadinc);
		//else if ( floor(TItime + loadinc + 1e-12) != floor(TItime)) 
		//	TItime=floor(TItime + loadinc);
		else
			TItime += loadinc;
#endif

		// hit endtime, do not overshoot
		if (TItime > solset.endtime) 
		{
			TItime = solset.endtime;
		}

// obsolete - loadfactor raises from 0 to 1 in each computation step
		//////////// load factor is only used for relaxation of error tolerance goal; loadfact has to be <= 1
		//////////loadfact = TItime;
		//////////if (loadfact > 1) loadfact = 1;

		// compute step time - set variable for later use (by elements, loads, ...) in Evaluate() functions
		ComputationStepTime() = GetCSTime(TItime);
		loadfact = ComputationStepTime();
	
		//+++++++++++++++++++++++++++++
		// reset counters
		lastnewtonits = 0;
		lastjaccnt = 0;

		//nonlinear iteration (e.g. for plastic components):
		int nonlinfin = 0;
		nlit = 0;
		int maxlni = 0;

		//+++++++++++++++++++++++++++++
		//save start condition at iteration
		TISaveState(TIlaststep_state);


		//+++++++++++++++++++++++++++++
		StartTimeStep();

		// AND: AP: TItemp_state is used to store plastic strains at beginning of nonlinear iteration,
		//      after transport which occurs in StartTimerStep
		// otherwise: TItemp_state is not used so far
		TISaveState(TItemp_state);
		TISaveState(TIlastnonlinit_state);

		double init_error = -1;

		//+++++++++++++++++++++++++++++
		// while
		//		nonlinear interation is not finished
		//		return value is 1
		//		calculation is not stopped
		// do
		//		Nonlinear iteration
		while (!nonlinfin && rv == 1 && !StopCalculation())
		{
			//save start condition at nonlinear iteration (e.g. for error estimation)
			//TISaveState(TIlastnonlinit_state); 

			int sos_rs = GetSecondOrderSize_RS(); //size of resorted second order equations
			numsol.SetBandSize(sos_rs - reducedbandsize); //for banded-solver, which part is banded!

			nlit++;
			log.TInonlinit = nlit;

			// set initial values
			for (int i=1; i <= sos; i++) GSxs(i) = TIx0(i);
			for (int i=1; i <= es;  i++) GSxs(i+sos) = TIx0(i+2*sos);
			for (int i=1; i <= is;  i++) GSxs(i+sos+es) = TIx0(i+2*sos+es);

			numsol.NLSolveInfo() = 1;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// call Newton, NLF and Jacobian modified for flag DoStaticComputation
			// return value is 1 if Newton method was successful

			rv = numsol.NLSolve(GSxs);

			// set interation counter for number of newton iterations (only maximal number is retained)
			int lni = numsol.GetNewtonIts();
			if(lni>maxlni)	maxlni = lni;

			lastnewtonits += lni;
			lastjaccnt += numsol.GetJacCount();
			if(UO().GetGlobalMessageLevel()==UO_LVL_dbg2)	//$ DR 2011-09-15
			{
				UO(UO_LVL_dbg2) << "rv=" << rv << "\n";
			}
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// if Newton was successful
			if (rv == 1)
			{
				//.............................
				// set local vector for nonlinear iteration step
				xh.SetLen(2*sos+es+is);
				for (int i=1; i <= sos; i++) xh(i) = GSxs(i);
				for (int i=1; i <= sos; i++) xh(i+sos) = TIx0(i+sos);
				for (int i=1; i <= es; i++) xh(i+2*sos) = GSxs(i+sos);
				for (int i=1; i <= is; i++) xh(i+2*sos+es) = GSxs(i+sos+es);

				SetStartVector(xh);


				//.............................
				// nonlinear step for inelastic strain
				//
				// inelastic update:
				// nlerror = ||yield||
				// inel_strain += tau * yield * direction_of_yield
				
				double nlerr = PostNewtonStep(TItime); 

				//.............................
				if ( init_error < 0 ) init_error = nlerr;
				// if ( nlit%1 == 0) UO() << "step " << nlit << ", nonlinerror = " << nlerr << ", discontacc " << DiscontinuousAccuracy() <<  "\n";

				//.............................
				// relax accuracy in iterations with small loadfactor
				double relaxfact = 1;
				if (solset.static_use_tol_relax_factor)
					relaxfact = Sqr(1./loadfact);
				double max_relaxfactor = solset.static_maxtol_relax_factor;
				if (relaxfact > max_relaxfactor) relaxfact = max_relaxfactor;

				//.............................
				double tolgoal = DiscontinuousAccuracy();

				if (loadfact <= 1-1e-8) // same accuracy in loadfact as in TItime <= (endtime-1e-8)!
					tolgoal *= relaxfact; 

				//.............................
				if ( nlerr <= tolgoal) 
				{
					// nonlinear iteration was successful, finished
					nonlinfin = 1;

					// postprocessing for nonlinear step
					PostprocessingStep();
					if (SolverPrintsDetailedOutput())
					{
						if (GetOptions()->LoggingOptions()->SolverStepSolutionVectorIncrement() && GetOptions()->LoggingOptions()->SolverStepSolutionVector())
						{
							SolverUO() << "step: solution vector = " << TIx0 << "\n";
						}
						if (GetOptions()->LoggingOptions()->SolverStepSolutionVectorIncrement())
						{
							SolverUO() << "step: increment vector = " << TIx0 - GSxs << "\n";
						}
					}
				}
				else
				{
					// print reached error_goal in percentage; for large systems each 10th iteration or for very large systems every iteration
					if ((ss > 2000 && (nlit % 10 == 9)) || (ss > 10000))
					{
						SolverUO() << "\rDit=" << nlit << ", discont. error reached " << tolgoal/nlerr*100. << "% of goal\n";
					}
				}

				//******** to draw in every iteration 
				//lastdraw = TIcomptime; drawnow++;
				//TIDrawAndStore();
				
				//******** to enable pause and resume by click
				//uo.InstantMessageText("PostNewtonStep");
				//if(UO().GetGlobalMessageLevel()==UO_LVL_dbg1)	//$ DR 2011-09-15
				//{
				//	UO(UO_LVL_dbg1) << "  disc. step " << nlit << ": newton_its=" << lni 
				//		<< /*", disc. err="*/", nlerr=" << nlerr 
				//		<< /*", tol. goal="*/", tol=" << tolgoal << "\n";
				//}
			}
			else
			{
				//UO(UO_LVL_dbg1) << "  newton did not converge\n";   //PG no more necessary, switch on LoggingOptions.Solver.general_information
			}

			if (nlit > solset.maxdiscontinuousit)
			{
				rv = 0;
			}
		}
		//+++++++++++++++++++++++++++++
		// increase counters
		TIit++;
		if ((ss > 2000 && (nlit >= 9)) || ss > 10000) UO() << "\r";

		newtonits += lastnewtonits;
		jaccnt += lastjaccnt;
		if (nlit != 0)
		{
			lastnewtonits /= (int)sqrt((double)nlit);
			lastjaccnt /= (int)sqrt((double)nlit);
		}

		//+++++++++++++++++++++++++++++
		// adapting load increment size
		if (!rv || lastnewtonits > 20*5 || lastjaccnt > 15*5 || nlit > 20*5) // not converged, too many iterations needed
		{
#ifdef COMPILE_AND
			// AP: TIRestoreState only in case of rv=0 (i.e. in case of no convergence of time step
			if (!rv)
#endif
			TIRestoreState(TIlaststep_state);

			successcnt = -1;
			failcnt++;
			if (SolverPrintsDetailedOutput())
			{
				mystr str("");
				if (!rv)
				{
					if (nlit > solset.maxdiscontinuousit)
					{
						str+=mystr("  * Discontinuous.max_iterations have been reached\n");
					}
					else
					{
						str+=mystr("  * Newton did not converge\n");
					}
				}
				if (lastnewtonits > 15*5) str+=mystr("  * more than 75 Jacobians during this step\n");
				if (lastnewtonits > 20*5) str+=mystr("  * more than 100 Newton iterations during this step\n");
				if (nlit > 20*5) str+=mystr("  * more than 100 Discontinuous iterations during this step\n");
				SolverUO() << mystr("step: failcnt++, since\n")+str;
			}

			if (loadinc > solset.static_minloadinc)
			{
				if (!rv) TItime = TItime_last_converged_step; // adapted artificial time to last converged step		
				
				loadinc /= solset.static_loadinc_down;				// reduce load increment
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: load increment is multiplied by ")+mystr(1/solset.static_loadinc_down)+mystr("\n");

				if (failcnt>=2)
				{
					loadinc *= 0.25;
					if (SolverPrintsDetailedOutput()) SolverUO() << "step: load increment is moreover multiplied by 0.25, since failcnt > 2\n";
				}
				rv = 1;
			}
			else if (0) //0=save mode; try to overcome servere instabilities!!!!  buggy
			{
				//just keep on going ...
				rv = 1;
				successcnt = 0;
			}

			// hit minimal load increment, do not "under"shoot
			if (loadinc < solset.static_minloadinc)
			{
				loadinc = solset.static_minloadinc;
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: load increment is set to minimum value of ")+mystr(solset.static_minloadinc)+mystr(", since it's been already smaller than that.\n");
			}
		}
#ifndef COMPILE_AND
		else if (lastjaccnt <= 7 && nlit <= 4) // fast converged
#else
		else if (lastjaccnt <= 7 && nlit <= 4*3) // AP
#endif
		{
			failcnt = 0;
			successcnt ++;
			if (SolverPrintsDetailedOutput())	SolverUO() << mystr("step: failcnt=0 and successcnt++, since lastjaccnt=")+mystr(lastjaccnt)+mystr("<7, and nlit=")+mystr(nlit)+mystr("<4.\n");

			// increase load increment after solset.static_incloadinc successfull steps
			if (successcnt >= solset.static_incloadinc) 
			{
				successcnt = 0;
				loadinc *= solset.static_loadinc_up;
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: load increment is multiplied by ")+mystr(solset.static_loadinc_up)+mystr("\n");
			}

			// hit maximal load increment, do not overshoot	
			if (loadinc > solset.static_maxloadinc)
			{
				loadinc = solset.static_maxloadinc;
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step: load increment is set to max value of ")+mystr(solset.static_maxloadinc)+mystr(", since it's been already bigger than that.\n");
			}
		}
		else
		{			
			if (SolverPrintsDetailedOutput() && failcnt != 0) SolverUO() << mystr("step: failcnt set to 0 again.\n");
			failcnt = 0;
		}

		//+++++++++++++++++++++++++++++
		TIx0draw = TIx0;



		//+++++++++++++++++++++++++++++
		lt = GetClockTime();
		TIcomptime = start_clock_time + lt;

		//+++++++++++++++++++++++++++++
		solset.withgraphics = GetIOption(104);
		//redraw frequency: 0..off, 1..draw last frame, 2..100sec, 3..20sec, 4..2sec, 5..200ms, 6..50ms, 7..20ms, 8..every 10 frames, 9..every frame
		if (solset.withgraphics == 1 && TItime > solset.endtime-0.1*TIminstep) {lastdraw = TItime; drawnow++;}
		if (solset.withgraphics == 2 && (TIcomptime - lastdraw) > 100.) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 3 && (TIcomptime - lastdraw) >  20.) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 4 && (TIcomptime - lastdraw) >   2.) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 5 && (TIcomptime - lastdraw) > 0.2) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 6 && (TIcomptime - lastdraw) > 0.05) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 7 && (TIcomptime - lastdraw) > 0.02) {lastdraw = TIcomptime; drawnow++;}
		if (solset.withgraphics == 8 && (TIit - lastdraw) >= 10) {lastdraw = TIit; drawnow++;}


		//+++++++++++++++++++++++++++++
		int endflag = 0;
		if ((TItime >= solset.endtime-1e-8) || rv == 0 || StopCalculation() || TIfinished!=0 )
		{
			drawnow = 1;
			endflag = 1;
#ifdef COMPILE_AND
			// AP: for AND, return 1 for successful computation, -1 else
			if (rv) TIfinished=1;
			else TIfinished=-1;
#endif
		}

		//+++++++++++++++++++++++++++++
		if (SolverPrintsDetailedOutput() || (TIcomptime - lastprintres) > GetOptions()->LoggingOptions()->ComputationOutputEveryXSec() || endflag) //update text every GetDOption(  4) seconds
		{
			lastprintres = TIcomptime;

			TIV.timetogo2 = TIV.timetogo;
			TIV.timetogo = TIcomptime / TItime * (solset.endtime-TItime);
			TIV.timetogo = 0.2*TIV.timetogo+0.8*TIV.timetogo2;

			SolverUO() <<
				("step=")+mystr(TIit)+
				(", f=")+mystr(loadfact)+
				(", finc=")+mystr(loadinc)+
				//mystr(", x1=")+mystr(GetCharacteristicSol())+
				(", Nits=")+mystr(maxlni)+
				(", Dits=")+mystr(nlit)+
				(", Jacs=")+mystr(jaccnt)+
				(", real_time=")+mystr(GetTString(TIcomptime))+
				(", real_time_togo=")+mystr(GetTString(TIV.timetogo))+
				("\n");
		}

		if(UO().GetGlobalMessageLevel()==UO_LVL_dbg2 && SolverPrintsDetailedOutput())	//$ DR 2011-09-15
		{
			UO(UO_LVL_dbg2) << "step lastdraw = " << lastdraw << ", TIcomptime = " << TIcomptime << ", drawnow = " << drawnow << "\n";
		}
		TIDrawAndStore();
		
		//+++++++++++++++++++++++++++++
		if (!failcnt) WriteSol();
		
		//+++++++++++++++++++++++++++++
		EndTimeStep();
#ifdef computationsteps
// additional code for change of computation step here
		if(IsComputationStepFinished(TItime))
		{
			ComputationStepFinished();

			ComputationStepNumber()++;                              // increase ComputationStepNumber by one since next TItime is in next Step 
			ApplyComputationStepSettings(ComputationStepNumber());
			// get correct load inc for new ComputationStep
			loadinc = Minimum(solset.static_initloadinc,solset.static_maxloadinc);
		}
#endif

		if (SolverPrintsDetailedOutput())
		{
			SolverUO() << mystr("**************** step ")+mystr(TIit)+mystr(" (static, adaptive) end\n");
		}
	}

	//--------------------------------------------------------------
	if (rv == 1 && TItime == solset.endtime) 
	{
		SolverUO() << "static computation completed successfully!\n";
	}

	SetNewtonItSum(newtonits);
	SetJacCount(jaccnt);

	SolverUO() << newtonits << " Newton iterations and " << jaccnt << " Jacobians computed\n";

	//--------------------------------------------------------------
	xh.SetLen(2*sos+es+is);
	for (int i=1; i <= sos; i++) xh(i) = GSxs(i);
	for (int i=1; i <= sos; i++) xh(i+sos) = TIx0(i+sos);
	for (int i=1; i <= es; i++) xh(i+2*sos) = GSxs(i+sos);
	for (int i=1; i <= is; i++) xh(i+2*sos+es) = GSxs(i+sos+es);

	SetStartVector(xh);
	
	if (SolverPrintsDetailedOutput())
	{
		if (GetOptions()->LoggingOptions()->SolverStepSolutionVector())
		{
			SolverUO() << "solution vector = " << TIx0 << "\n";
		}
	}
	
	//--------------------------------------------------------------
	loadfact = 1.;
	TItime = solset.endtime;
	NLS_ModifiedNewton() = oldmodnewton;

	//--------------------------------------------------------------
	logout.flush();
	return rv;
}

int TimeInt::ExplicitIntegration()
{
	int rv = 1;

	if (GetImplicitSize() != 0) 
	{
		uo << "ERROR: explicit integration failed: algebraic equations cannot be solved with explicit solver!" << ENDL;
	}


	drawnow = 0;

	if (TIV.outputlevel > 2)
		uo << "explicit integration: mid-point rule, order = 1" << ENDL;

	//write solution:
	FirstStep();

	while (((solset.endtime-TItime)>0.1*TIminstep) && rv != 0 && !StopCalculation() && TIfinished==0)
	{
		if (SolverPrintsDetailedOutput())
		{
			SolverUO() << mystr("**************** step ")+mystr(TIit+1)+mystr(" (explicit, adaptive) begin\n");
		}

		if(RelApproxi(solset.endtime,TItime+TIstep))
		{
			TIstep = solset.endtime - TItime;
		}

		//uo << "step = " << TIit << ", time=" << TItime << ", step=" << TIstep << "\n";
		rv = AdaptiveExplicitStep();
		TIit++;

		StepPrintAndDraw(rv);

		if (SolverPrintsDetailedOutput())
		{
			SolverUO() << mystr("**************** step ")+mystr(TIit)+mystr(" (explicit, adaptive) end\n");
		}
	}

	if (TIV.outputlevel >= 2)
	{
		uo << "computational time = " << TIcomptime << " seconds" << ENDL;
		uo << "rejected steps = " << TIrejectedsteps << ENDL;
	}
	if (TIit != 0 && TIcomptime!=0) 
	{
		if (TIV.outputlevel >= 3)
		{
			uo << "timesteps per second = " << (double)TIit/TIcomptime << " steps/s" << ENDL;
		}
	}
	return rv;
}

int TimeInt::AdaptiveExplicitStep()
{
	TISaveState(TIlaststep_state); //save state of end of last step (or initial step)

	StartTimeStep();
	
	TISaveState(TItemp_state); // AP: state after StartTimeStep(), for AND this is after transport
	TIstep = TIstepnew;
	double savestep = TIstep;
	log.TInewtonit = 0;
	double TItimeold = TItime;
	int rv = 0;

	double it = 0;
	while (it < 15 && !rv)
	{
		it++;
		if (SolverPrintsDetailedOutput()) SolverUO() << mystr("step it=")+mystr(it)+mystr(", step increment=")+mystr(TIstep)+mystr(", time=")+mystr(TItime)+mystr("\n");

		SetStepRecommendation(1e100);

		if (TItime+TIstep > solset.endtime) TIstep = solset.endtime - TItime;

		rv = ExplicitStep(); 

		if (!rv && it < 15) 
		{
			TIRestoreState(TIlaststep_state);   //set to old TIx0, reset all internal nonlinear variables
			if (GetStepRecommendation() < TIstep) 
			{
				TIstep = GetStepRecommendation();
				if (it > 3) 
				{
					TIstep *= 0.5;
					if (SolverPrintsDetailedOutput()) SolverUO() << "step increment multiplied by 0.5\n";
				}
				else if (it > 7)
				{
					TIstep *= 0.1;
					if (SolverPrintsDetailedOutput()) SolverUO() << "step increment multiplied by 0.1\n";
				}
			}
			else
			{
				TIstep *= 0.5;
				if (SolverPrintsDetailedOutput()) SolverUO() << "step increment multiplied by 0.5\n";
				
				if (TIstep <= TIminstep)
				{
					it = 14;
					if (SolverPrintsDetailedOutput()) SolverUO() << mystr("perform last step, since TIstep=")+mystr(TIstep)+mystr(" is less equal TIminstep=")+mystr(TIminstep)+mystr("\n");
				}
			}

			if (TIstep < TIminstep) 
			{
				TIstep = TIminstep; 
				if (SolverPrintsDetailedOutput()) SolverUO() << mystr("TIstep set to TIminstep=")+mystr(TIminstep)+mystr(", since it's been already smaller than that.\n");
			}

		}
	}

	if (!rv)
	{
		if (!minstepwarned) uo << "ERROR: Step at minimal stepsize, step not fully converged. Further warning suppressed!" << ENDL;
		minstepwarned = 1;
		//rv = 1;
	}

	FinishStep();
	TIstep = savestep;
	EndTimeStep();

	return rv;
}

int TimeInt::ExplicitStep()
{
	TMStartTimer(21);
	int rv = 1;
	double corerr;
	TInonls = 0;

	TIdiscontstep = 0;
	TISaveState(TIlastnonlinit_state);

	rv = GeneralEStep();
	corerr = PostNewtonStep(TItime);

	if (corerr == -1) 
	{
		rv = 0;
	}
		
	if (GetStepRecommendation() < TIstep && TIstep > TIminstep)
	{
		rv = 0;
	}

	if (rv != 0) 
	{
		PostprocessingStep();
	}

	TMStopTimer(21);

	return rv;
}

//only for test: mass matrix
Matrix expl_mass;

int TimeInt::GeneralEStep()
{
	//macht einen Schritt

	int ss = GetSystemSize();
	int is = GetImplicitSize();
	int es = GetFirstOrderSize();
	int sos = GetSecondOrderSize();
	int sos_rs = GetSecondOrderSize_RS();

	int rv = 1;

	int i,j;
	IRK_Tableau* t = tableaus[act_tableau];
	int ns = t->nstages;

	int oldMinv = 0;

	//explicit Runge Kutta!!!

	evalf.SetLen(sos);
	TIu0.SetLen(sos);
	TIu0.Copy(TIx0,1,1,sos);
	TIv0.SetLen(sos);
	TIv0.Copy(TIx0,1+sos,1,sos);
	TIu0e.SetLen(es);
	TIu0e.Copy(TIx0,1+2*sos,1,es);

	//build gu, gv, gue
	gu[1] = TIu0;
	gv[1] = TIv0;
	gue[1] = TIu0e;

	for (i=1; i <= ns; i++)
	{
		PrecomputeEvalFunctions(); //$ PG 2012-1-15: precomputation for EvalF(), EvalF2(), EvalG(), EvalM(), and also for CqTLambda() in case of constraints

		//evaluate right hand sides for first and second order equations
		kv[i].SetLen(sos);
		kue[i].SetLen(es);
		xi[i].SetLen(2*sos+es);

		xi[i].Copy(gu[i],1,1,sos);
		xi[i].Copy(gv[i],1,1+sos,sos);
		xi[i].Copy(gue[i],1,1+2*sos,es);

		//compute RK-f(x,t) for gu,gue and gv 
		ku[i] = gv[i];
		evalfe.SetLen(es);
		if (es != 0)
		{
			EvalF(xi[i],kue[i],TItime+TIstep*t->c(i));
		}

		if (sos != 0)
		{
			if (oldMinv)
			{
				expl_mass.SetSize(sos,sos);
				expl_mass.FillWithZeros();
				evalf.SetLen(sos);
				EvalM(xi[i], expl_mass, TItime+TIstep*t->c(i));
				EvalF2(xi[i], evalf, TItime+TIstep*t->c(i));

				TMStartTimer(11); //solve
				int lrv = expl_mass.Solve(evalf,kv[i]);
				TMStopTimer(11); //solve
				if (!lrv) {uo << "Error in GeneralEStep: Mass matrix is singular!" << ENDL; return 0;}
			}
			else
			{
				EvalMinvF2(xi[i], kv[i], TItime+TIstep*t->c(i));
			}
		}

		//compute next stage
		if (i < ns)
		{
			gu[i+1] = TIu0;
			gv[i+1] = TIv0;
			gue[i+1] = TIu0e;
			for (j = 1; j <= i; j++)
			{
				gu[i+1] += TIstep*t->A(i+1,j)*ku[j];
				gv[i+1] += TIstep*t->A(i+1,j)*kv[j];
				gue[i+1]+= TIstep*t->A(i+1,j)*kue[j];
			}
		}
	}

	//compute evaluation step
	for (i = 1; i <= ns; i++)
	{
		TIu0  += TIstep*t->b(i)*ku[i];
		TIv0  += TIstep*t->b(i)*kv[i];
		TIu0e += TIstep*t->b(i)*kue[i];
	}

	//compute new TIx0:
	TIx0.Copy(TIu0,1,1,sos);
	TIx0.Copy(TIv0,1,1+sos,sos);
	TIx0.Copy(TIu0e,1,1+2*sos,es);

	return rv;
}

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //# Computation Steps with changing SolverOptions...
	//Time mapping

// returns number of the Computation Step for given global time // returns -1 it globaltime < 0 or globaltime > last step_end_time
int TimeInt::GetCSNumber(double globaltime)
{
//rv :   -1         if t < 0
//       n+1        if t > last_endtime
//       i [1..n}   "else" - if in 0 <= t <= last_endtime

	if(globaltime < 0.)
		return -1; // error - time negative !

// in first ComputationStep
	if(globaltime <= CSEndTimes(1)) 
		return 1;   

// after last computation step endtime 
	if(globaltime > CSEndTimes.Last()) // includes 0.0
		return CSEndTimes.Length()+1;

// most likely step for actual time
	int selected = ComputationStepNumber(); 
  if(selected < 1)
		selected = 1; 

// in selected ComputationStep
	if(CSEndTimes(selected-1) < globaltime && globaltime <= CSEndTimes(selected)) // includes last_endtime
		return selected;    

// in i-th ComputationStep (search array)
	for(int i=2; i <= CSEndTimes.Length(); i++)
	{
		if(CSEndTimes(i-1) < globaltime && globaltime <= CSEndTimes(i)) // includes last_endtime
		return i;     
	}

// global time is larger then last ComputationStep endtime...
	return -1; 
}

// returns relative time in corresponding step   rv(0..1] // returns 0. only for globaltime == 0 or globaltime > last step_end_time
double TimeInt::GetCSTime(double globaltime)
{
//rv:    -1            if t < 0
//       0             if t == 0
//       1.            if t > last_endtime
//       t' (0,1.]     "else" t is percentage within interval

	if(globaltime < 0.)
		return -1.; // error - time negative !

	if(globaltime == 0.)
		return 0.; // only way to obtain "0.00000000000000000000000"

	int stepnumber = GetCSNumber(globaltime);
	if(stepnumber > CSEndTimes.Length()) // stepnumber "-1" can not be rv ! (catched 	by "if(globaltime < 0.)")
		return 1.;  // global time is larger then last ComputationStep endtime...

// time within Computation Step intervals
	double intervalstart = 0.0;  // valid for first interval
	if (stepnumber > 1)
		intervalstart = CSEndTimes(stepnumber-1); // valid for all other intervals
	double intervalend = CSEndTimes(stepnumber);

	if(globaltime >= intervalend - 1e-12)
		return 1.; // get a good "1.0000000000000000000000000000" at end of interval

	double relativetime = (globaltime - intervalstart) / (intervalend - intervalstart);
	return relativetime; // (0 .. 1)
}

// returns 1 if globaltime is any step_end_time
int TimeInt::IsComputationStepFinished(double globaltime)
{
	if( GetCSTime(globaltime) == 1. && globaltime <= CSEndTimes.Last() ) 
		return 1.;

  return 0;
}





