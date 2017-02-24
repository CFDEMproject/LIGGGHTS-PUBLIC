//#**************************************************************
//#
//# filename:             writesolution.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:						September 2010
//# description:          
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
 
#include "mbs.h"
#include "myfile.h"
#include "element.h"
#include "sensors.h"

//#include "windows.h"					//for shell execute
#include <direct.h>						//for getcwd
#include  <io.h>							//file operation _access
#include  <stdio.h>
#include  <stdlib.h>
#include <windows.h>


void MultiBodySystem::WriteSolDataInfo()
{	//$ PG 2012-4-19: write info-file if solution data has been saved to files instead of memory	
	if (solset.sol_data_to_files)
	{
		mystr sol_data_dir = solset.sol_directory+mystr("solution_data\\");
		CMFile info_file(sol_data_dir+mystr("info.txt"), TFMwrite);
		mystr header(GetSolutionDataInfoFileHeaderString());

		info_file.RW(header);
		info_file.RW(mystr("\n"));

		int n_steps = GetNumberOfSolutionDataSteps();
		info_file.RWInt(n_steps);
		info_file.RW(mystr("\n"));
		
		//write solver options edc
		
		// ...
	}
}

int MultiBodySystem::ReadSolDataInfo()
{	//$ PG 2012-4-19: read info-file (if solution data is saved to files instead of memory), and return number of stored states	
	if (solset.sol_data_to_files)
	{
		mystr sol_data_dir = solset.sol_directory+mystr("solution_data\\");
		mystr info_file_name = sol_data_dir+mystr("info.txt");
		CMatrixFile info_file(info_file_name, TFMread);
		UO(UO_LVL_ext) << mystr("Read solution data info from ")+info_file_name+mystr("\n");

		info_file.ReadHeader();
		
		int n_steps=0;
		info_file.RWInt(n_steps);
		SetNumberOfSolutionDataSteps(n_steps);
		UO(UO_LVL_ext) << mystr("Number or solution data steps set to ")+mystr(n_steps)+mystr(".\n");
		
		//read solver options edc
		
		// ...

		return  n_steps;
	}
	return 0;
}

void MultiBodySystem::WriteSol()
{
	log.writesolcnt++;

	//do sensor computations in every step!
	SetActState(GetSolVector());
	for (int i=1; i <= NSensors(); i++)
	{
		Sensor & s = GetSensor(i);
		s.Evaluate(GetTime());
	}

	if (solset.writeresults == 0 || GetSolSet().writesolstep == 0) return;//$ RL 2011-10-20: avoid crashes due to modulo zero operations

	if ((log.writesolcnt-1)%GetSolSet().writesolstep != 0 /*&& !DoStaticComputation()*/) return;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// write header
	if(WriteSolutionFileHeader())
	{
		// write header
		sol << GetSolutionFileHeaderString().c_str();
		isSolutionFileHeaderWritten = 1;
	}

	//const Vector& x = GetSolVector();

	int sos = GetSecondOrderSize();
	sol.precision(16);
	char strtime[100];
	sprintf(strtime,"% 12.9f ", GetTime());
	sol << strtime;


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//output sensor data:

	for (int i=1; i <= NSensors(); i++)
	{
		if (GetSensor(i).GetSignalStorageMode() & SSM_CommonFile)
		{
			sol.precision(GetSensor(i).GetPrecision());
			sol << GetSensor(i).GetLastValue() << " ";
		}
		if (GetSensor(i).GetSignalStorageMode() & SSM_OwnFile)
		{
			ofstream* ssol = GetSensor(i).GetOwnOutputFile();
			ssol->precision(GetSensor(i).GetPrecision());
			(*ssol) << strtime;
			(*ssol) << GetSensor(i).GetLastValue() << " ";
			(*ssol) << "\n";
			if (GetSolSet().sol_immediately_write) (*ssol) << flush;
		}
	}
	sol.precision(16); //changed in sensor writing

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (0) //write kinetic, potential and total energy
	{
		double EK;
		EK = GetKineticEnergy();
		double EP;
		EP = GetPotentialEnergy();

		sol << EK << " " << EP << " " << EK+EP << " ";
	}

	sol << "\n";
	if (GetSolSet().sol_immediately_write) sol << flush;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mystr MultiBodySystem::GetSolutionDataInfoFileHeaderString() // returns the header string for the solution-data info file
{
	mystr str("%HOTINT_Solution-Data_Info-File\n");
	if(GetModelFunctionsIndex0()!=-1)
	{
		str = str + "%" + GetModelsLibrary()->GetModelInterface(GetModelFunctionsIndex0())->GetMBSModelName() + "\n";	
	}
	else
	{
		str = str + "%\n"; //do not write model name because it does not exist!     OLD: "Please select \"ModelFunction\" and press \"Save configuration\" button. Then restart HOTINT.\n";
	}
  time_t t;	time(&t);  struct tm *ts;  ts = localtime(&t); // get current date and time
	str = str + "%" + asctime(ts);
	str = str + "%Comment: " + GetSolSet().sol_file_header_comment + "\n";
	str = str + "%\n"; // unused yet

	//mystr colnumbers("%1    ");
	//mystr   colnames("%Time ");
	//for(int i=1; i <= NSensors(); i++)
	//{
	//	if(GetSensor(i).GetWriteResults() & 1)
	//	{
	//		colnames = colnames + GetSensor(i).GetSensorName() + " ";
	//		colnumbers = colnumbers + mystr(i+1);
	//		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
	//			colnumbers = colnumbers + " ";
	//	}
	//}
	//str = str + colnumbers + mystr("\n") + colnames + mystr("\n");
	return str;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mystr MultiBodySystem::GetSolutionFileHeaderString() // returns the header string
{
	mystr str("%HOTINT_Solution_File\n");
	if(GetModelFunctionsIndex0()!=-1)
	{
		str = str + "%" + GetModelsLibrary()->GetModelInterface(GetModelFunctionsIndex0())->GetMBSModelName() + "\n";	
	}
	else
	{
		str = str + "%\n"; //do not write model name because it does not exist!     OLD: "Please select \"ModelFunction\" and press \"Save configuration\" button. Then restart HOTINT.\n";
	}

  time_t t;	time(&t);  struct tm *ts;  ts = localtime(&t); // get current date and time
	str = str + "%" + asctime(ts);
	str = str + "%Comment: " + GetSolSet().sol_file_header_comment + "\n";
	str = str + "%\n"; // unused yet

	mystr colnumbers("%1    ");
	mystr   colnames("%Time ");
	for(int i=1; i <= NSensors(); i++)
	{
		if(GetSensor(i).GetSignalStorageMode() & SSM_CommonFile)
		{
			colnames = colnames + GetSensor(i).GetSensorName() + " ";
			colnumbers = colnumbers + mystr(i+1);
			for(int j=colnumbers.Length(); j<colnames.Length(); j++)
				colnumbers = colnumbers + " ";
		}
	}
	str = str + colnumbers + mystr("\n") + colnames + mystr("\n");
	return str;
};


mystr MultiBodySystem::GetSolParFileHeaderString() // returns the header string for SolParFile
{
	mystr str("%HOTINT_Parameter_Solution_File\n");
	if(GetModelFunctionsIndex0()!=-1)
	{
		str = str + "%" + GetModelsLibrary()->GetModelInterface(GetModelFunctionsIndex0())->GetMBSModelName() + "\n";
	}
	else
	{
		str = str + "%\n"; //do not write model name because it does not exist!     OLD: "Please select \"ModelFunction\" and press \"Save configuration\" button. Then restart HOTINT.\n";
	}
  time_t t;	time(&t);  struct tm *ts;  ts = localtime(&t); // get current date and time
	str = str + "%" + asctime(ts);
	str = str + "%Comment: " + GetSolSet().sol_file_header_comment + "\n";
	str = str + "%\n"; // unused yet

	mystr colnumbers = "%";
	mystr colnames = "%";
	int nr_vp = 0;

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.activate"))
	{
		colnames = colnames + GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.MBS_EDC_variable_name") + mystr(" ");
		colnumbers = colnumbers + mystr(nr_vp+1) + " ";
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
			colnumbers = colnumbers + " ";
		nr_vp++;
	}
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.ParameterVariation.Var2.activate"))
	{
		colnames = colnames + GetMBS_EDC_Options()->TreeGetString("SolverOptions.ParameterVariation.Var2.MBS_EDC_variable_name") + mystr(" ");
		colnumbers = colnumbers + mystr(nr_vp+1) + " ";
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
			colnumbers = colnumbers + " ";
		nr_vp++;
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_second_order_size")) 
	{
		colnames   = colnames   + "second_order_size  ";
		colnumbers = colnumbers + mystr(nr_vp+1) + " ";
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
		colnumbers = colnumbers + " ";
		nr_vp++;
	}
                                       
	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_CPU_time")) 
	{
		colnames   = colnames   + "computation_time  ";
		colnumbers = colnumbers + mystr(nr_vp+1) + " ";
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
		colnumbers = colnumbers + " ";
		nr_vp++;
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Eigensolver.do_eigenmode_computation")) // MSax 2013-07-04 : added eigenvalues, for campbell diagram
	{
		for(int i=1; i <= eigval.Length(); i++)
		{
			colnames = colnames + "eigval" + mystr(i) + " ";
			colnumbers = colnumbers + mystr(nr_vp+1) + " ";
			for(int j=colnumbers.Length(); j<colnames.Length(); j++)
				colnumbers = colnumbers + " ";
			nr_vp++;
		}
	}

	for(int i=1; i <= NSensors(); i++)
	{
		if(GetSensor(i).GetSignalStorageMode() & SSM_CommonFile)
		{
			colnames = colnames + GetSensor(i).GetSensorName() + " ";
			colnumbers = colnumbers + mystr(nr_vp+1) + " ";
			for(int j=colnumbers.Length(); j<colnames.Length(); j++)
				colnumbers = colnumbers + " ";
			nr_vp++;
		}
	}

	if (GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.ParameterFile.write_cost_function")) //$ RL 2011-3-22: OPTIONAL default: yes //$ MS 2011-3-22: 
	{
		colnames   = colnames   + "costfunction_value  ";
		colnumbers = colnumbers + mystr(nr_vp+1) + " ";
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
		colnumbers = colnumbers + " ";
		nr_vp++;
	}
	str = str + colnumbers + mystr("\n") + colnames + mystr("\n");
	return str;
};

bool MultiBodySystem::LogFileNameChanged() { return strcmp(logout_name,GetOptions()->LoggingOptions()->DefaultLogFilename()); }
void MultiBodySystem::OpenLogFile(int mode)
{
	// if log file not open: open it
	// else, if log file name has changed, then if
	// mode == 0: rename log file
	// mode == 1: close old log file and open new log file
	if (!logout.is_open())
	{
		logout_name = GetOptions()->LoggingOptions()->DefaultLogFilename();
		logout.open(logout_name.c_str(), ios::out);
	}
	else
	{
		if (LogFileNameChanged())
		{
			assert (mode == 0 || mode == 1);
			//get new log file name, and save old one
			mystr old_logout_name(logout_name);
			logout_name = GetOptions()->LoggingOptions()->DefaultLogFilename();
			
			//print information
			if (mode == 1)
			{
				logout << "\n<logfile continued in " << logout_name << ">\n";
			}

			//close old log file
			logout.close();

			//if file with new file name already exists, then backup it
			mystr new_name = mystr("backup.") + mystr((int)time(NULL)) + mystr('.') + mystr(logout_name);
			bool file_existed = ( rename(logout_name , new_name.c_str()) == 0 );

			//optionally, rename old log file with new name
			if (mode == 0)
			{
				rename(old_logout_name , logout_name);
			}

			//open new log file
			logout.open(logout_name.c_str(), ios::out); //log file, with read sharing allowed!
			
			//print information
			if (mode == 0)
			{
				if (file_existed)
				{
					logout << "old " << logout_name << " was renamed to " << new_name << "\n";
				}
			}
			else
			{
				logout << "<logfile continued from " << old_logout_name << ">";
				if (file_existed)
				{
					logout << " -- old " << logout_name << " was renamed to " << new_name;
				}
				logout << "\n\n";
			}
		}
	}
}




void MultiBodySystem::OpenFiles(int flag)
{
	mystr dir = GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path"); //old: GetTOption80);
	mystr soln = GetMBS_EDC_Options()->TreeGetString("SolverOptions.Solution.SolutionFile.output_filename"); //GetTOption81);	// solution file
	mystr solparn = GetMBS_EDC_Options()->TreeGetString("SolverOptions.Solution.ParameterFile.parameter_variation_filename"); //GetTOption82);

	UO().CallWCDriverFunction(4); //Set Programm directory as current

	UO() << "open files with flag=" << flag << "\n";

	int opt = ios::out;
	int opt2 = ios::out;

	if (flag&1) //append to existing file
	{
		opt = ios::out | ios::app;
		opt2 = ios::out | ios::app;
	}
	if (flag&2) //append to existing file
	{
		opt2 = ios::out | ios::app;
	}

#ifndef COMPILE_AND
	int rv;
	rv = _access(dir.c_str(), 0);
	if (rv == -1 && dir.Length() != 0)
	{
		bool directoryCreated = CreateDirectoryA(dir.c_str(), 0);
		if (!directoryCreated)
		{
			UO() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			UO() << "ERROR: the output directory '" << dir << "' does not exist, and it could not be created. No output will be written!!!\n";
			UO() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		}
	}
#endif

#ifdef my_new_stdiostream
	sol.open((dir+soln).c_str(), opt); //with read sharing allowed!
	solpar.open((dir+solparn).c_str(), opt2); //with read sharing allowed!
#else
	sol.open((dir+soln).c_str(), opt, filebuf::sh_read); //with read sharing allowed!
	solpar.open((dir+solparn).c_str(), opt2, filebuf::sh_read); //with read sharing allowed!
#endif

	//sensor files:
	for (int i=1; i <= NSensors(); i++)
	{
		if(GetSensor(i).GetSignalStorageMode() & SSM_OwnFile)
		{
			ofstream* file = new ofstream();
#ifdef my_new_stdiostream
			//file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetTypeName()+mystr(".txt")).c_str(), opt);
			file->open((dir+mystr("S")+mystr(i)+mystr("-")+GetSensor(i).GetSensorName()+mystr(".txt")).c_str(), opt);
#else
			//file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetTypeName()+mystr(".txt")).c_str(), opt, filebuf::sh_read);
			file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetSensorName()+mystr(".txt")).c_str(), opt, filebuf::sh_read);
#endif
			GetSensor(i).SetOwnOutputFile(file);

			//$ RL 2011-7-14:[ create file for results of fft-computation
			//$ YV 2012-06: files for sensor post-computation evaluation results are not opened in advance any longer
			/*
			if(GetSensor(i).IsSensorComputationFFT())
			{
				ofstream* fft_file = new ofstream();
  #ifdef my_new_stdiostream
				//fft_file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetTypeName()+mystr("-fft.txt")).c_str(), opt);
				fft_file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetSensorName()+mystr("-fft.txt")).c_str(), opt);
  #else		
				//fft_file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetTypeName()+mystr("-fft.txt")).c_str(), opt, filebuf::sh_read);
				fft_file->open((dir+mystr("S")+mystr(GetSensor(i).GetSensorNumber())+mystr("-")+GetSensor(i).GetSensorName()+mystr("-fft.txt")).c_str(), opt, filebuf::sh_read);
  #endif
				GetSensor(i).SetOutputFileFFT(fft_file);				
			}
			*/
			//$ RL 2011-7-14:]
		}
	}
}

void MultiBodySystem::CloseLogFile()
{
	logout.close();
}

void MultiBodySystem::CloseFiles()
{
	sol.close();
	solpar.close();
	logout.close();

	for (int i=1; i <= NSensors(); i++)
	{
		if (GetSensor(i).GetSignalStorageMode() & SSM_OwnFile)
		{
			if(GetSensor(i).GetOwnOutputFile())
			{
				GetSensor(i).GetOwnOutputFile()->close();
				delete GetSensor(i).GetOwnOutputFile();
				GetSensor(i).SetOwnOutputFile(NULL);
			}
			//$ RL 2011-7-14:[
			//$ YV 2012-06: files for sensor post-computation evaluation results are closed immediately after the data is written
			/*
			if(GetSensor(i).GetOutputFileFFT())
			{
				GetSensor(i).GetOutputFileFFT()->close();
				delete GetSensor(i).GetOutputFileFFT();
				GetSensor(i).SetOutputFileFFT(0);
			}
			*/
			//$ RL 2011-7-14:]
		}
	}
}


