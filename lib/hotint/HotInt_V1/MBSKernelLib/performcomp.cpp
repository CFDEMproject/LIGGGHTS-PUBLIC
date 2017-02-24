//#**************************************************************
//#
//# filename:             performcomp.cpp
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
 
#include "MBS.h"
#include "myfile.h"
#include "sensitivity.h"
#include "optimization.h"
#include "element.h"
#include "sensors.h"

////#include "..\SuperLU\slu_ddefs.h"
//#include "..\SuperLU\SuperLUmain.h"
//
////ofstream sol("..\\..\\output\\sol.txt"); //Ö
//ofstream solpar2("..\\output\\par_rest.txt");

void MultiBodySystem::PerformComputation()
{
	//check if parameter variation or compute eigenmodes:
	if(ptr2PerformComputation_Function != 0)
	{
		UO(UO_LVL_ext) << "\n****************************\nDo user-defined computation!\n****************************\n\n";
		CallPtr2PerformComputation_Function(this);
	}
	else if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.activate"))
	{
		UO(UO_LVL_ext) << "\n****************************\nDo parameter variation\n****************************\n\n";
		PerformParameterVariation();
	}
	else if(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.activate"))
	{
		UO(UO_LVL_ext) << "\n*****************************************\nDo parameter optimization\n*****************************************\n\n";
		Optimization(this).PerformOptimization();		//$ YV 2012-11-28
	}
	else if(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Sensitivity.activate"))
	{
		UO(UO_LVL_ext) << "\n*****************************************\nDo sensitivity analysis\n*****************************************\n\n";		
		Sensitivity(this).PerformSensitivityAnalysis();
	}
	else 
	{
		//normal mode: do a single computation
		PerformSingleStaticDynamicComputation();
	}

	simulationStatus.AddSetStatusFlag(TSimulationProcessFinished); //$ MaSch 2013-08-19
	ComputationFinished();

	//if option is set==>close application after computation!
	if (GetMBS_EDC_Options()->TreeGetInt("GeneralOptions.Application.close_application_when_finished")) // GetIOption(220))
	{
		UO().CallWCDriverFunction(5); //quit
	}
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MultiBodySystem::PerformSingleStaticDynamicComputation()
{
	int rv = PreComputationOperations();
	if (!rv) return;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	rv = ComputeSystem(); //Main solution procedure !!! (time integration/static solver)
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if (!rv) simulationStatus.SetStatusFlag(TSimulationStoppedDueError); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19 
	
	//OpenFiles() is done in following procedure:
	PostComputationOperations(rv);

	WriteSolutionVector();

	DoFinalSensorEvaluations();

	CloseFiles();

	pCFB->FinishedComputation();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int MultiBodySystem::ComputeSystem()
{
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Eigensolver.do_eigenmode_computation"))
	{
		UO(UO_LVL_ext) << "\n****************************\nDo eigenvalue/eigenmode computation\n****************************\n\n";
		ComputeEigenmodes();
		return 1;
	}
	else
	{
		return Integrate();
	}
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//$!PG 2011-2-22:[  new local function in order to handle more cases of parameter variation (also varstep<0 and varstepfact<1) 	
bool MultiBodySystem::DoParameterVariationComputation(bool const isgeometric, const double par, const double varstep, const double varstepfact, const double varend) const
{
	double some_space = 1e-8;
	if (isgeometric)
	{
		if ((varstepfact >= 1) && (par <= varend*(1.+some_space)))
			return true;
		else if ((varstepfact < 1) && (par >= varend*(1.-some_space)))
			return true;
	}
	else
	{
		if ((varstep >= 0) && (par <= varend*(1.+some_space)))
			return true;
		else if ((varstep < 0) && (par >= varend*(1.-some_space)))
			return true;
	}

	return false;
}
//$!PG 2011-2-22:]

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MultiBodySystem::PerformParameterVariation()
{
	//###################################################################################################################
	//
	//     do parametric variation
	//
	//###################################################################################################################
	

	
	
	isSolParFileHeaderWritten = 0;

	double varstart = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.start_value"); //
	double varend = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.end_value"); //=11;

	double varstep = 0; //!!!!!!!!!!!!!!!!!!
	double varstepfact = 1; //1 = incremental

	//$!PG 2011-2-22:[  new passage in order to handle more cases of parameter variation (also varstep<0 and varstepfact<1) 
	bool isgeometric = false;
	//$!PG 2011-2-22:]

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.geometric") == 0) 
	{ //arithmetic
		varstep = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.arithmetic_step");
	}
	else //geometric 
	{
		//$!PG 2011-2-22:[  new line in order to handle more cases of parameter variation (also varstep<0 and varstepfact<1) 
		isgeometric = true;
		//$!PG 2011-2-22:]

		varstep = 0;
		varstepfact = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.geometric_step");
	}

	//###################################################################################################################
	//second variation

	double varstart2 = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.Var2.start_value"); //
	double varend2 = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.Var2.end_value"); //=11;

	double varstep2 = 0; //!!!!!!!!!!!!!!!!!!
	double varstepfact2 = 1; //1 = incremental

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.Var2.geometric") == 0) 
	{ //arithmetic
		varstep2 = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.Var2.arithmetic_step");
	}
	else //geometric 
	{
		varstep2 = 0;
		varstepfact2 = GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.ParameterVariation.Var2.geometric_step");
	}

	varparameter2 = varstart2;
	int var2_activate = GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.Var2.activate");

	if (!var2_activate)
	{
		varstart2 = 0;
		varparameter2 = 0;
		varend2 = 0;
		varstep2 = 1;
		varstepfact2 = 1;
	}

	while (DoParameterVariationComputation(isgeometric, varparameter2, varstep2, varstepfact2, varend2) && !StopCalculation())
	{
		if (var2_activate)
		{
			mystr varname2 = GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.Var2.MBS_EDC_variable_name");
			// RL 2011-02:[ TODO: function of class PerformComputation should be used --> no duplicate code!
			if (varname2 != mystr(""))
			{				
				if (GetMBS_EDC_Options()->TreeFind(varname2))
				{
					GetMBS_EDC_Options()->TreeFind(varname2)->SetLocked(0);	//$ DR 2013-04-09: parameter variation for script language
					GetMBS_EDC_Options()->TreeSetDouble(varname2, varparameter2);
					GetMBS_EDC_Options()->TreeFind(varname2)->SetLocked(1); //$ DR 2013-04-09: parameter variation for script language
				}
				else if (GetModelDataContainer() != 0 && GetModelDataContainer()->TreeFind(varname2))
				{
					GetModelDataContainer()->TreeFind(varname2)->SetLocked(0);	//$ DR 2013-04-09: parameter variation for script language
					GetModelDataContainer()->TreeSetDouble(varname2, varparameter2);		
					GetModelDataContainer()->TreeFind(varname2)->SetLocked(1); //$ DR 2013-04-09: parameter variation for script language
				} 
				else
				{
					UO() << "*********************\nWarning: can not access variable '" << varname2 << "' for parameter variation!\n*********************\n";
				}
			}
			// RL 2011-02:] TODO: function of class PerformComputation should be used --> no duplicate code!

			UO(UO_LVL_all) << "\n\n\n\n********************\n********************\n********************\n********************\n";
			UO(UO_LVL_all) << "\n********************\nStart simulation with varparameter 2\n" << varname2 << " = " << varparameter2 
				<< "\n********************\n********************\n********************\n********************\n";
		}

		double par = varstart;


		//while (par <= varend*(1.+1e-8) && !StopCalculation())
		//$!PG 2011-2-22: why only for (par <= varend)? parameter variation should also handle the case (varstepfact < 1)
		//$!PG 2011-2-22: Instead I implemented MultiBodySystem::DoParameterVariationComputation(...)
		while (DoParameterVariationComputation(isgeometric, par, varstep, varstepfact, varend) && !StopCalculation())
		{
			varparameter = par;	//par and varparameter is effectively the same

			GetMBS_EDC_Options()->TreeSetDouble("SolverOptions.start_time", 0.);

			mystr varname = GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.MBS_EDC_variable_name");
			if (varname != mystr(""))
			{
				if (GetMBS_EDC_Options()->TreeFind(varname))
				{
					GetMBS_EDC_Options()->TreeFind(varname)->SetLocked(0); //$ DR 2013-04-09: parameter variation for script language
					GetMBS_EDC_Options()->TreeSetDouble(varname, varparameter);
					GetMBS_EDC_Options()->TreeFind(varname)->SetLocked(1); //$ DR 2013-04-09: parameter variation for script language
				}
				else if (GetModelDataContainer() != 0 && GetModelDataContainer()->TreeFind(varname))
				{
					GetModelDataContainer()->TreeFind(varname)->SetLocked(0); //$ DR 2013-04-09: parameter variation for script language
					GetModelDataContainer()->TreeSetDouble(varname, varparameter);
					GetModelDataContainer()->TreeFind(varname)->SetLocked(1); //$ DR 2013-04-09: parameter variation for script language
				} 
				else
				{
					UO() << "*********************\nWarning: can not access variable '" << varname << "' for parameter variation!\n*********************\n";
				}
			}


			UO(UO_LVL_all) << "\n********************\nStart simulation with varparameter\n" << varname << " = " << varparameter << "\n********************\n";

			//********************************************************************
			//********************************************************************
			//repeatedly call computation:
			RepeatedlyPerformComputation();

			//********************************************************************
			//next iteration:
			par += varstep;
			par *= varstepfact;
		}
		//********************************************************************
		//next iteration, varparameter2:
		varparameter2 += varstep2;
		varparameter2 *= varstepfact2;
	}

	pCFB->FinishedComputation();

	if (GetMBS_EDC_Options()->TreeGetInt("GeneralOptions.Application.close_application_when_finished")) // GetIOption(220))
	{
		UO().CallWCDriverFunction(5); //quit
	}

}

double MultiBodySystem::EvalSensorCostFunctionVal() //$ RL 2011-02: evaluate cost function (e.g. for parameter optimization)
{
	double cost_function_val = 0.;
	if (GetTIFinished() == -1)//$ MS 2011-9-8: if time-integration failed  //$ YV -> MS: 2011?
	{
		cost_function_val += 1.e10;
	}
	else
	{
		for (int i=1; i <= NSensors(); i++)
		{
			if(GetSensor(i).HasSensorProcessingEvaluationData())
			{				
				TArray<double> signalProcessingEvaluationData = GetSensor(i).GetSignalProcessingEvaluationData();
				for(int i = 1; i <= signalProcessingEvaluationData.Length(); i++)
					cost_function_val += signalProcessingEvaluationData(i);

				if(GetSolSet().paropt_sensors != 0)
				{
					UO(UO_LVL_warn) << "Depriciated data from sensor-processing is used as cost function.\n Change option to \"SolverOptions.Optimization.sensors=0\".\n";
				}
			}
			// for script language, the end value of the sensor is used
			else if(GetSolSet().paropt_sensors == i)
			{
				cost_function_val += GetSensor(i).GetLastValue();
			}
		}
	}
	return cost_function_val;
}

int MultiBodySystem::RepeatedlyPerformComputation(int flags)
{

	//Model must be always initialized! ==> always isfirst
	int isfirst = isfirstcomputation;
	Initialize();
	isfirstcomputation = isfirst; //because is overwritten in PreComputationOperations / Initialize

	//OpenFiles() is done in following procedure:
	int rv = PreComputationOperations(1+2); //1..always set initial conditions; 2..always append to parameter solution file
	if(!rv)
	{
		UO(UO_LVL_err).InstantMessageText("Error in pre-computation operations, simulation stopped!\n");
		return rv; //$ RL 2011-8-5: if system is inconsistent, computation must not be started,
	}

	isfirstcomputation = 0; //append to files after this point!
	rv = ComputeSystem();
	if (!rv) simulationStatus.SetStatusFlag(TSimulationStoppedDueError); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19 
	PostComputationOperations(rv);
	if (!(flags&1)) //decide wheter to write solpar file directly
	{
		DoFinalSensorEvaluations();
	}
	CloseFiles();

	return rv;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int MultiBodySystem::PreComputationOperations(int flags)
{
	if (isfirstcomputation)
	{
		//isfirstcomputation = 0;

		UO(UO_LVL_ext) << "solution output file=" << GetSolSet().sol_filename << "\n";
		UO(UO_LVL_ext) << "output directory = " << GetSolSet().sol_directory << "\n";
		
		//$ PG 2012-4-19:[ for saving solution data to files
		if (solset.sol_data_to_files)
		{
			mystr sol_data_dir = solset.sol_directory+mystr("solution_data\\");
			mystr sol_data_dir_name("solution data directory"); 
			int rv = UOfull().pUI->CreateMissingDirectory(sol_data_dir.c_str());
			switch (rv)
			{
			case 0:				// path to directory is wrong
				UO(UO_LVL_err) << mystr("ERROR: ") + sol_data_dir_name + mystr(" cannot be created, since path ") + solset.sol_directory + mystr(" is wrong!\n");
				return 0;
			case 1:				// directory existed
				UO(UO_LVL_ext) << sol_data_dir_name + mystr(" = ") + sol_data_dir + mystr("\n");
				break;
			case 2:				// directory created
				UO(UO_LVL_ext) << sol_data_dir_name + mystr(" = ") + sol_data_dir + mystr(" (created) \n");
				break;
			default:
				assert(0);
				break;
			}
		}
		//$ PG 2012-4-19:] for saving solution data to files 


		//UO(UO_LVL_ext) << "log file = " << GetSolSet().log_filename << "\n";  //$ PG 2012-2-3: modified (see next line)
		UO(UO_LVL_ext) << "log file = " << GetOptions()->LoggingOptions()->DefaultLogFilename() << "\n";

		int appendflag = 0;
		if ((flags&2)) appendflag += 2; //always append to parameter file (e.g. in parameter variation)
		OpenFiles(appendflag);

		SetInitialConditions();

	}
	else
	{
		//the following option is used for continuation of a computation (PAUSE): ==> NEEDED???
		if (GetMBS_EDC_Options()->TreeGetBool("SolverOptions.Timeint.reset_after_simulation") || (flags&1) != 0)
		{
			SetInitialConditions(); //this is set in Initialize(), for the case that initial con
		}

		//decide whether to append solution to existing file, or to replace it:
		int appendflag = 1;
		if (GetSolSet().sol_replace_files) appendflag = 0;
		if ((flags&2)) appendflag += 2; //always append to parameter file (e.g. in parameter variation)

		OpenFiles(appendflag);
		//PG: OpenLogFile(1);
	}

	//$ DR 2012-11-08 moved from "if (isfirstcomputation)-branch" to here. Consistency should be checked always when there are changes!
	int modified = UOfull().pUI->GetModelModified();
	if(isfirstcomputation || modified)
	{
		isinconsistent = CheckSystemConsistency();
		if (isinconsistent) 
		{
			UO(UO_LVL_err).InstantMessageText("Error: Due to an inconsistent system, all output files are closed!\n"); //$ RL 2012-5-22: user should now, what happens. 
			CloseFiles();
			return 0;
		}
	}
	isfirstcomputation = 0;

	SetMaxIndex(GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Timeint.max_index"));
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.do_static_computation")) SetMaxIndex(3); //dostaticcomputation

	simulationStatus.SetStatusFlag(TSimulationRunning); //$ MaSch 2013-11-20

	TMResetTimer(); //reset timers for evaluation of function computational time
	TMStartTimer(1); //timer for total computational time, only solution/integration!!!
	return 1;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void MultiBodySystem::PostComputationOperations(int integrate_rv)
{

	sol << flush;
	solpar << flush;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //$ RL 2011-7-14: FFT code
	//$ RL 2011-7-14:[ FFT-sensor files are created
	for(int i=1; i<=sensors.Length();i++)
	{
		if(GetSensor(i).GetOwnOutputFile())
		{
			(*GetSensor(i).GetOwnOutputFile()) << flush;
		}

		// now we perform post-computation evaluation of the sensor signal

		if(GetSensor(i).NeedsFileForPostComputationProcessing())
		{
			// we need to treat the file with the results here
			ofstream post_computation_output_file;
			mystr dir = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path");
			post_computation_output_file.open((dir+mystr("S")+mystr(i)+mystr("-")+GetSensor(i).GetSensorName()+mystr("-post_comp.txt")).c_str(), ios::out);
			GetSensor(i).ApplyPostComputationSensorProcessing(&post_computation_output_file);
			post_computation_output_file.close();
		}
		else
			GetSensor(i).ApplyPostComputationSensorProcessing(NULL);		// postcomputation processing is performed without an output file stream
	}
	//$ RL 2011-7-14:] FFT-sensor files are created
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	

	TMStopTimer(1);

	//$ PG 2012-4-19: write info-file if solution data has been saved to files instead of memory
	WriteSolDataInfo();

	if (!integrate_rv)
	{
		UO(UO_LVL_err) << "*** time integration not finished! ***" << "\n";
		simulationStatus.SetStatusFlag(TSimulationStoppedDueError); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19 
	}
	else 
	{
		if (!StopCalculation())
		{
			UO(UO_LVL_err) << "computation finished\n";
			simulationStatus.SetStatusFlag(TSimulationEndedRegularly); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19
		}
		else
		{
			UO(UO_LVL_err) << "computation stopped by user\n";
			simulationStatus.SetStatusFlag(TSimulationStoppedByUser); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19
		}
	}

	UO(UO_LVL_err) << "---------------------------------------------------------------\nCOMPUTATION SUMMARY:\n";
	UO(UO_LVL_err) << "computational time = " << TIcomptime << " seconds" << ENDL;
	UO(UO_LVL_err) << "rejected steps = " << TIrejectedsteps << ENDL;
	if (TIit != 0 && TIcomptime!=0) UO(UO_LVL_all) << "timesteps per second = " << (double)TIit/TIcomptime << " steps/s" << ENDL;

	UO(UO_LVL_err) << "number of timesteps = "  << GetTIit()<< "\n";
	UO(UO_LVL_err) << "number of jacobians = " << GetJacCount() << "\n";
	//UO(UO_LVL_ext) << FullNewtonCnt() << " Full Jacobians computed." << "\n";
	UO(UO_LVL_dbg1) << "changes in the timestep = " << GetStepChanges() << "\n";
	UO(UO_LVL_err) << "total newton iterations = " << GetNewtonItSum() << "\n";
	UO(UO_LVL_err) << "full evaluations of EvalF2 = " << log.evalfcnt << "\n";
	UO(UO_LVL_err) << "partial Jacobi evaluations of EvalF2 = " << log.evalf_jaccnt << "\n";
#ifdef gencnt
	UO(UO_LVL_err) << "allocated matrices = " << GenMatCnt() << "\n";
	UO(UO_LVL_err) << "allocated vectors = " << GenVecCnt() << "\n";
#endif
	UO(UO_LVL_err) << "---------------------------------------------------------------\n";

	//if (StopCalculation()) return; //this can not be earlier - otherwise the timer is not accurate and the output is not done on manual stop

	ComputationFinished();
}

void MultiBodySystem::DoFinalSensorEvaluations() //do final sensor evaluations (min/max/diff, etc.) for sensors after computation
{
	if (StopCalculation()) return;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//output sensor data at final step:

	if(DoFinalEigenValueComputation())
	{
		CloseFiles();
		double simulation_time = GetTime(); //$ PG 2013-11-27: memo needed, since time is set to length of eigensystem (-> for vizualizing Eigenmodes with DataManager)
		// sensor of type TSEigenValue exist
		ComputeEigenmodes();

		//$ PG 2013-11-26: perform recalculation of all system sensors after this final eigenvalue computation (at least for all ev-systemsensors this is required)
		for (int i=1; i <= NSensors(); i++)
		{
			Sensor& s = GetSensor(i);
			if(s.IsSensorOnEigensystem())
			{
				s.Evaluate(simulation_time);
			}
		}

		int append = 1;
		OpenFiles(append);
		//PG: OpenLogFile(1);
	}

	int is_sensor_computation = 0;
	for (int i=1; i <= NSensors(); i++)
	{
		// we need to know if the sensor has post evaluation data;
		if(GetSensor(i).HasSensorProcessingEvaluationData())
			is_sensor_computation = 1;
	}
	
	solpar.precision(16);

	if(WriteSolParFileHeader())
	{
		solpar << GetSolParFileHeaderString();
		isSolParFileHeaderWritten = 1;
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.activate")) {solpar << varparameter << " ";}
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.Var2.activate")) {solpar << varparameter2 << " ";}

	//write size of system (only second order equations):
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_second_order_size")) //$ RL 2011-02: OPTIONAL default: no
	{
		solpar << GetSecondOrderSize() << " "; 
	}

                                            
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_CPU_time")) //$ RL 2011-02: OPTIONAL default: no
	{
		//write CPU-time:
		solpar << TMGetTimer(1) << " "; //computational time (CPU-costs)
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Eigensolver.do_eigenmode_computation")) // MSax 2013-07-04 : added eigenvalues, for campbell diagram
	{
		for(int i=1; i <= eigval.Length(); i++)
		{
			solpar << eigval(i)/(2.*MY_PI) << " ";
		}
	}

	//write final sensor values:
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_final_sensor_values")) //$ RL 2011-02: OPTIONAL default: yes
	{
		for (int i=1; i <= NSensors(); i++)
		{
			if(GetSensor(i).GetSignalStorageMode() != SSM_None)
			{
				//$ YV 2012-06: the sensors may produce just one scalar value
				solpar << GetSensor(i).GetLastValue() << " ";
			}
		}
	}

	//write sensor evaluation values:
	if (is_sensor_computation)
	{
		UO(UO_LVL_ext) << "***************************\nSensor computations:\n";
		for (int i=1; i <= NSensors(); i++)
		{
			if(GetSensor(i).HasSensorProcessingEvaluationData())
			{
				UO(UO_LVL_ext) << "sensor " << i << " : [ ";
				TArray<double> & signalProcessingEvaluationData = GetSensor(i).GetSignalProcessingEvaluationData();
				for (int j = 1; j <= signalProcessingEvaluationData.Length(); j++)
				{
					double v = signalProcessingEvaluationData(i);
					solpar << v << " ";
					UO(UO_LVL_ext) << v << " ";
				}
				UO(UO_LVL_ext) << "]\n";
			}
		}
		UO(UO_LVL_ext) << "***************************\n";
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_cost_function")) //$ RL 2011-02: OPTIONAL default: yes
	{
		// write cost function:
		double costfunc = EvalSensorCostFunctionVal();		
		solpar << costfunc << " ";
	}

	solpar << "\n" << flush;

}

void MultiBodySystem::WriteSolutionVector() //write solution vector into file if option is set
{
	if (StopCalculation()) return;
	
	//store solution
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.store_solution_state"))
	{
		mystr path = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path");

		int strlen = path.Length();
		if (strlen != 0)
		{

			if (path[strlen-1] != '\\' && path[strlen-1] != '/') path += mystr('\\');
			mystr file = path+GetMBS_EDC_Options()->TreeGetString("SolverOptions.Solution.store_solution_state_name");

			//OLD: mystr file = path+mystr("store_out.txt");
			ofstream storeout(file.c_str()); 

			UO(UO_LVL_ext) << "final solution vector stored in file '" << file << "'\n";

			Vector x = GetSolVector();

			if ((int)storeout.good() == 0 || storeout.fail())
			{
				UO(UO_LVL_err) << "ERROR: could not save solution, check filename!!!\n";
			}
			else
			{

				storeout.precision(19);
				(storeout) << x.Length() << "\n";
				for (int i=1; i <= x.Length(); i++)
				{
					(storeout) << x(i) << " ";
				}
				(storeout) << endl;
				(storeout) << GetTime() << endl;
			}
			storeout.close();
		}
	}

}

