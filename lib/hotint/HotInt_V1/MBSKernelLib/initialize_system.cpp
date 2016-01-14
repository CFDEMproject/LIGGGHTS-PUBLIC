//#**************************************************************
//#
//# filename:             initialize_system.cpp
//# 
//# author:               Gerstmayr Johannes
//#
//# generated:						September 2010
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
#include "element.h"
#include "sensors.h"


void MultiBodySystem::Initialize() //is called upon full reset of the computation (erase data in data manager, assemble system, new initial conditions, etc.
{
	//PG: OpenLogFile(0);
	SetProhibitRedraw(1);
	Destroy(); //if called several times!!!
	simulationStatus.SetStatusFlag(TSimulationNotStarted); //$ MS + RL 2012-2-29: //$ MaSch 2013-08-19
	MBS_EDC_TreeGetDouble(GetSolSet().withgraphics,"ViewingOptions.Misc.redraw_frequency"); //should be erased, this exists twice in mbssolset and in i/d/toptions
	
	SetComputationSolverOptions(); //copy general options to MBSSolSet(); this is done because old models use MBSSolSetOptions!!
	SetEDCOptions2Options(); //copy numbered i/d/t options to EDC Tree options

//$ AH 2011-11-28: reset pointer to perform computation function
	ptr2PerformComputation_Function = 0;

	isfirstcomputation = 1;
	isSolutionFileHeaderWritten = 0;
	isSolParFileHeaderWritten = 0;

	colref.Init();

	//without that you can not write from an element to uout!!!!
	uout.pUI = uo.pUI;

  // hook UserOutput object to solver settings
	UOfull().HookToSolverSettings(& GetOptions()->LoggingOptions()->OutputLevel(), & GetOptions()->LoggingOptions()->FileOutputLevel(), &logout, & GetOptions()->LoggingOptions()->CriticalLogFileSize(), & GetOptions()->LoggingOptions()->MaxErrorMessages(), & GetOptions()->LoggingOptions()->MaxWarningMessages(),
		& GetOptions()->LoggingOptions()->OutputPrecisionDouble(), & GetOptions()->LoggingOptions()->OutputPrecisionVector(), & GetOptions()->LoggingOptions()->OutputPrecisionMatrix()); //$ PG 2012-2-3: substituted MBSSolSet().output_level by GetOptions()->LoggingOptions()->OutputLevel()
	uout.HookToSolverSettings(& GetOptions()->LoggingOptions()->OutputLevel(), & GetOptions()->LoggingOptions()->FileOutputLevel(), &logout, & GetOptions()->LoggingOptions()->CriticalLogFileSize(), & GetOptions()->LoggingOptions()->MaxErrorMessages(), & GetOptions()->LoggingOptions()->MaxWarningMessages(),
		& GetOptions()->LoggingOptions()->OutputPrecisionDouble(), & GetOptions()->LoggingOptions()->OutputPrecisionVector(), & GetOptions()->LoggingOptions()->OutputPrecisionMatrix()); //$ PG 2012-2-3: substituted MBSSolSet().output_level by GetOptions()->LoggingOptions()->OutputLevel()
	
	int modelfile_good = 1; //if error occurs, cancel initialization

	
	//$ RL 2012-7-31:[ this is necessary for drag&drop and doule click action of Hotint Input Data files (.hid)
	//                 ReadModelData is called here, because in the selected model name should be defined in this hid-file 
	//                 in order to select the model from the file
	int model_from_file = 0;	//$ DR 2013-02-21 
	mystr hotint_input_data_file = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.ModelFile.hotint_input_data_filename","");
	if(hotint_input_data_file.Length() > 0 /*&& !IsModelData_Initialized()*/)	//$ DR removed "&& !IsModelData_Initialized("
	{
		model_from_file = 1;
		modelfile_good = ReadModelData(hotint_input_data_file);		//$ DR 2013-02-21 added "modelfile_good =" 
		Assemble();
		SetModelData_Initialized(1);
		// AddRecentFile(filename); //DR not working yet
	}
	//$ RL 2012-7-31:]
	
	if (!model_from_file)
	{
		//Choose Model (or empty model):
		int nmodel = GetModelsLibrary()->GetModelsCount();
		UO() << "number of models=" << nmodel << "\n";

		UO() << "selected model name=" << GetOptions()->GeneralOptions()->ModelFileInternalModelFunctionName() << "\n";
		int nm = GetModelFunctionsIndex0();
		//GetIOption12 = nm;

		UO() << "selected model number=" << nm << "\n";


		//if ((nm != -1) && (!initialized_from_file))
		if (nm != -1)
		{
			if (!IsModelData_Initialized())
			{
				ElementDataContainer old_modeldata_edc;

				if (GetModelsLibrary()->GetModelInterface(nm)->HasMBSModelInitData()) 
				{
					modelfile_good = GetModelsLibrary()->GetModelInterface(nm)->InitializeMBSModelData(this);
					if (modelfile_good)
					{
						SetModelData_Initialized(1);
					}
					else
					{
						UO().InstantMessageText("WARNING: your InitModelData function returned 0, which indicates an error. The Initialization is terminated! see HOTINT change: 2012-02-21: JG/PG ==> set return value to 1");
					}			
				}
				else
				{
					if (GetSolSet().default_model_data_file.Length() != 0)
					{
						UO(UO_LVL_ext) << "Load default model data file ='" << GetSolSet().default_model_data_file << "\n";
						ReadModelData(GetSolSet().default_model_data_file);
					}
				}

				//write old model data == arguments from startup
				if (modeldata_edc_args.Length())
				{
					GetModelDataContainer()->TreeReplaceEDCDataWith(&modeldata_edc_args);
					modeldata_edc_args.Reset(); //reset ==> not used any more at next call
				}
				if (mbs_edc_args.Length())
				{
					//at this point of initialization, the i/d/t options and the MBSSolSet() contain actual configuration&option data
					SetOptions2EDCOptions(); //copy numbered i/d/t options to EDC Tree options
					SetSolverDialogOptions(); //copy MBSSolSet() to general options; this is done because old models use MBSSolSetOptions!!
					GetMBS_EDC_Options()->TreeReplaceEDCDataWith(&mbs_edc_args);
					SetEDCOptions2Options(); //copy EDC Tree to numbered i/d/t options for fast access
					SetComputationSolverOptions(); //copy EDC Tree to MBSSolSet()

					mbs_edc_args.Reset(); //reset ==> not used any more at next call
				}


			}

			if (modelfile_good) //$JG2012-02-21
			{
				if (!bComputeEigenmodes)		//$ AD 2011-09-14: do not reassign direct Sensor watches when calculating Eigenmodes...
				{
					// close all instances of plottool
					uo.pUI->CallWCDriverFunction(9); // close all instances of PlotToolDlgs
				} //$JG2012-02-21

				//create model
				modelfile_good = GetModelsLibrary()->GetModelInterface(nm)->CreateMBSModel(this);
				if (!modelfile_good)
				{
					UO().InstantMessageText("WARNING: your model function returned 0, which indicates an error. The Initialization is terminated! see HOTINT change: 2012-02-21: JG/PG ==> set return value to 1");
				}
				// create all sensorwatches //$JG2012-02-21
				if (!bComputeEigenmodes)
				{
					for (int i=1; i<=sensors.Length(); i++)
					{
						if(sensors(i)->GetOpenSensorWatchPlotToolAtStartUpFlag() )
						{
							uo.pUI->CallWCDriverFunction(8,0,i);
						}
					}
				}
			}
		}
		else
		{
			UO() << "No model selected!\n";
		}
	}


	if (modelfile_good) //$JG2012-02-21
	{
		SetInitialConditions();
	}
	else
	{
		Destroy();
	}

	SetSolverDialogOptions(); //copy MBSSolSet() to general options; this is done because old models use MBSSolSetOptions!!

	MBS_EDC_TreeSetDouble(GetSolSet().withgraphics,"ViewingOptions.Misc.redraw_frequency"); //should be erased, this exists twice in mbssolset and in i/d/toptions

	//UO() << "done\n";
	SetProhibitRedraw(0);

	//UO() << "finish initialization\n";

	//GetMBS_EDC_Options()->TreeGetString
	//mystr path = MBS_EDC_TreeGetString("GeneralOptions.Paths.application_path");	//$ DR 2013-05-27
	//UO(UO_LVL_multsim) << "Application path: " << path <<" \n"; 

}

void MultiBodySystem::SetInitialConditions()
{
	MBSsos = GetSecondOrderSize();
	MBSes = GetFirstOrderSize();
	MBSis = GetImplicitSize();
	MBSss = GetSystemSize();


	// Vector initialized with zero!:
	Vector x0(GetSystemSize());
	SetGlobalInitConditions(x0);
	//uo << "Global InitVector=" << x0 << "\n";

	Vector xdata(GetDataSize());
	SetGlobalInitData(xdata);

	SetTime(0.);								//reset TItime
	//simulationStatus = TSimulationRunning; //$ MS + RL 2012-2-29:
	GetMBS_EDC_Options()->TreeSetDouble("SolverOptions.start_time",0.);//reset start time, otherwise old end time is used as new start time

	SetDataVector(xdata); 
	SetStartVector(x0); 


	SetActState(GetSolVector());

	//after initial conditions are set
	//be careful with constraints! --> the original initial conditions (stored in the element) need to be set first!
	//-->also important for finite elements!
	for (int i=1; i<=elements.Length(); i++) 
	{
		elements(i)->Initialize();
	} 

	//load solution
	if (0) //(0 && !MBSSolSet().dostaticcomputation)
	{
		ifstream* storein = new ifstream;

#ifdef my_new_stdiostream
		int opt = ios::in;
		storein->open("..\\..\\output\\store_in13.txt", opt); 
#else
		int opt = ios::in|ios::nocreate;
		storein->open("..\\..\\output\\store_in13.txt", opt, filebuf::sh_read); 
#endif

		if (storein->fail() || !storein->good())
		{
			UO() << "ERROR: could not open stored solution!!!\n";
		}
		else
		{
			storein->seekg(0,ios::beg);

			int len = 0;
			(*storein) >> len;
			UO() << "len=" << len << "\n";
			UO() << "x0-len=" << x0.Length() << "\n";

			//if (len != x0.Length()) {UO() << "No stored solution available!!!\n";}
			if (len > x0.Length()) {UO() << "No stored solution available!!!\n";}
			else
			{
				for (int i=1; i <= len /*156 , 84 len*/; i++) //pipec 4:20, 8:36, 16:68; 32:132, 64:260   pipe2:84, pipe4: 156
				{
					(*storein) >> x0(i);
				}
				double timei;
				(*storein) >> timei;
				//SetTime(timei);

				UO() << "**********************\nStored solution loaded!!!\n**********************\n";
			}
		}
		storein->close();
		delete storein;
	}

	if (stored_initialconditions.Length() != 0 && stored_initialconditions.Length() <= x0.Length())
	{
		for (int i=1; i <= stored_initialconditions.Length(); i++)
		{
			x0(i) = stored_initialconditions(i);
		}
		UO() << "Stored solution loaded!!!\n";
	}

	//nore change state to a loaded initial vector:
	SetStartVector(x0); 
	SetActState(GetSolVector());

	ComputeMaxSparseBandwidth();
}

void MultiBodySystem::InitializeAfterConfigLoaded()
{
	TimeInt::InitializeAfterConfigLoaded();

	
	Initialize(); //initialize models


	
	//*JG, 29 Aug, 2010; this should not be necessary anymore:
	//if (NE() != 0)								//if empty model is loaded then do not override stored computational settings!!!
	//	SetSolverDialogOptions();		//copy general options to MBSsolset options

}


