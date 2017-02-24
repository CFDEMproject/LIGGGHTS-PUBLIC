//#**************************************************************
//#
//# filename:             sensitivity.cpp
//#
//# author:               Rafael Ludwig
//#
//# generated:						February 2011
//# description:          Sensitivity analysis of sensor values from certain model parameters.
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
 
#include "sensitivity.h"
#include "element.h"
#include "sensors.h"

mystr Sensitivity::GetHeaderString() const  //returns header for output file
{
	mystr str("%HOTINT_Sensitivity_Analysis_Solution_File\n");
	if(mbs->GetModelFunctionsIndex0()!=-1)
	{
		str = str + "%" + mbs->GetModelsLibrary()->GetModelInterface(mbs->GetModelFunctionsIndex0())->GetMBSModelName() + "\n";
	}
  time_t t;	time(&t);  struct tm *ts;  ts = localtime(&t); // get current date and time
	str = str + "%" + asctime(ts);
	str = str + "%Comment: " + mbs->GetSolSet().sol_file_header_comment + "\n";
	str = str + "%Compute sensitivity of sensor values with respect to parameters\n";
	
	mystr colnumbers, colnames;
	// sensors
	str = str + "%Sensors: (rows of sensitivity matrix)\n";
	colnames = "%";
	colnumbers = "%";
	for(int i=1; i<=sensor_names.Length();i++)
	{
		colnames = colnames + *sensor_names(i) + mystr(" ");
		colnumbers = colnumbers + mystr(i);
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
		{
			colnumbers = colnumbers + " ";
		}
	}
	str = str + colnumbers + mystr("\n");
	str = str + colnames + mystr("\n");	
	
	
	// parameters
	str = str + "%Parameters: (columns of sensitivity matrix)\n";
	colnumbers = "%";
	colnames = "%";
	for(int i=1; i<=param_names.Length();i++)
	{
		colnames = colnames + *param_names(i) + mystr(" ");
		colnumbers = colnumbers + mystr(i);
		for(int j=colnumbers.Length(); j<colnames.Length(); j++)
		{
			colnumbers = colnumbers + " ";
		}
	}
	str = str + colnumbers + mystr("\n");
	str = str + colnames + mystr("\n");
	return str;
}



void Sensitivity::Init()
{
	mbs->OpenFiles(0); //open files already here, do not append
	ElementDataContainer* edc = mbs->GetMBS_EDC_Options();

	mystr method = edc->TreeGetString("SolverOptions.Sensitivity.method");

	// choose type of analysis
	if(method.Compare(mystr("Forward")))
	{
		type_diff = 1;  //  1 ... df/dx =~ (f(x+dx)-f(x))/dx        forward differentiation
	}
	else if(method.Compare(mystr("Backward")))
	{
		type_diff = -1; // -1 ... df/dx =~ (f(x)-f(x-dx))/dx        backward differentiation
	}
	else
	{
		type_diff = 0;  //  0 ... df/dx =~ (f(x+dx/2)-f(x-dx/2))/dx central differentiation		
	}
	
	// get parameter names
	if(edc->TreeGetInt("SolverOptions.Sensitivity.use_optimization_parameters"))
	{
		number_of_params = edc->TreeGetInt("SolverOptions.Optimization.Parameters.number_of_params");
		for (int col=1; col <= number_of_params; col++)
		{
			mystr str =	edc->TreeGetString(mystr("SolverOptions.Optimization.Parameters.param_name")+mystr(col));
			param_names(col) = new mystr(str);
		}	
	}
	else
	{
		number_of_params = edc->TreeGetInt("SolverOptions.Sensitivity.Parameters.number_of_params");
		for (int col=1; col <= number_of_params; col++)
		{
			mystr str =	edc->TreeGetString(mystr("SolverOptions.Sensitivity.Parameters.param_name")+mystr(col));
			param_names(col) = new mystr(str);
		}
	}

	abs_diff_val = edc->TreeGetDouble("SolverOptions.Sensitivity.num_diff_parameter_absolute");
	rel_diff_val = edc->TreeGetDouble("SolverOptions.Sensitivity.num_diff_parameter_relative");

	NCompSensors = 0; // sensor with computation value (case use_final_sensor_values == 0)

	for(int i = 1; i<=mbs->NSensors();i++)
	{
		mystr str =	"";		
		if(mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Sensitivity.use_final_sensor_values",0))
		{
			str = mbs->GetSensor(i).GetSensorName();
			sensor_names(i) = new mystr(str); // store names for header
		}
		else
		{
			if(mbs->GetSensor(i).HasSensorProcessingEvaluationData())
			{
				str = mbs->GetSensor(i).GetSensorName();
				NCompSensors++;
				sensor_names(NCompSensors) = new mystr(str); // store names for header
			}
		}	
	}	

	mbs->SolParFile() << GetHeaderString(); // write header into output file (use sensor_names and parameter_names)

}

// get function values from sensors
void Sensitivity::GetFunctionValues(Vector& f)
{
	if(mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Sensitivity.use_final_sensor_values",0))
	{
		f.SetLen(mbs->NSensors());
		for(int i = 1; i<=mbs->NSensors(); i++)
		{
			//fill f with all sensor end values
			//$ YV 2012-06: the sensors may produce just one scalar value
			f(i) = mbs->GetSensor(i).GetLastValue(); // use final sensor values
		}
	}
	else
	{
		f.SetLen(NCompSensors);
		int i = 1; // index of f-vector
		for(int isens = 1; isens<=mbs->NSensors();isens++)
		{
			if(mbs->GetSensor(i).HasSensorProcessingEvaluationData())			// the data is produced by sensor processors in MBS::PostComputationOperations()
			{
				f(i) = ((Vector)mbs->GetSensor(isens).GetSignalProcessingEvaluationData()).GetNorm(); // use computation sensor values
				i++;
			}
		}
	}
}


// main part of Sensitivity
int Sensitivity::PerformSensitivityAnalysis()
{
	Init(); // initialize variables
	// fill vector x with model parameters
	Vector x = 1.0 * GetModelParameters();
	// dx=K*x + D
	//Vector factors(1.0, 2.0, 3.0, 4.0, 5.0, 6.0); //test of parameter sensitivity w.r.t. other parameter set: 
	//Vector factors(1.0); // sensitivity w.r.t. nominal parameters
	int indexf_max = 0; // number of sensors of sensitivity analysis

	//test[for(int nFactor = 1; nFactor <= factors.Length(); nFactor++)
	//{
	//test]	x = factors(nFactor) * x; // multiply nominal parameters with factor (for sensitivity analysis with modified parameter set - useful, if nominal parameters "x" are in (local) minima)

	dx = rel_diff_val * x  + abs_diff_val * Vector(number_of_params, 0, 1);
	// do sensitivity computation
	//testmbs->SolParFile() <<  "\n%Evaluate sensitivity for nominal parameters multiplied with factor " << factors(nFactor) << "\n";

	for(int jPar = 1; jPar <= number_of_params; jPar++)
	{
		//++++++++++++++++++++++++++++++++++++
		//b: difference quotient variables
		// dfi/dxj =~ (fi(x2j)-fi(x1j))/dxj
		double dxj = dx(jPar);
		Vector x1j = x;
		Vector x2j = x;
		if(type_diff == 1)  //  1 ... df/dx =~ (f(x+dx)-f(x))/dx,       forward differentiation
		{
			//f2j = f(x+dx), f1j = f(x) 
			x2j(jPar) += dxj;
			x1j(jPar) += 0.;
		}
		else if(type_diff == -1) // -1 ... df/dx =~ (f(x)-f(x-dx))/dx         backward differentiation
		{
			//f2j = f(x), f1j = f(x-dx)
			x2j(jPar) += 0.;
			x1j(jPar) -= dxj;
		}
		else //type_diff == 0 ... df/dx =~ (f(x+dx/2)-f(x-dx/2))/dx,  central differentiation		
		{
			//f2j = f(x+dx/2), f1j = f(x-dx/2)
			x2j(jPar) += 0.5 * dxj;
			x1j(jPar) -= 0.5 * dxj;				
		}
		//e: difference quotient variables
		//++++++++++++++++++++++++++++++++++++

		//++++++++++++++++++++++++++++++++++++
		//b: compute f2j
		SetParameters(&x2j);		

		mbs->CloseFiles(); //because it is opened in following function
		//*****
		//repeatedly call computation:
		mbs->UO(UO_LVL_multsim) << mystr("\n****************************");
		mbs->UO(UO_LVL_multsim) << mystr("\nStart Simulation f(x2): \n\n");
		mbs->RepeatedlyPerformComputation(1); //do not write into solpar file
		mbs->OpenFiles(1); //append		
		//*****
		Vector f2i;
		GetFunctionValues(f2i); // get function values from sensors
		//e: compute f2i
		//++++++++++++++++++++++++++++++++++++

		//++++++++++++++++++++++++++++++++++++
		//b: compute f1i
		SetParameters(&x1j);		
		mbs->CloseFiles(); //because it is opened in following function
		//*****
		//repeatedly call computation:
		mbs->UO(UO_LVL_multsim) << mystr("\nStart Simulation f(x1): \n\n");
		mbs->RepeatedlyPerformComputation(1); //do not write into solpar file
		mbs->OpenFiles(1); //append
		//*****
		Vector f1i;
		GetFunctionValues(f1i); // get function values from sensors
		assert(f2i.Length() == f1i.Length());		
		//e: compute f2j
		//++++++++++++++++++++++++++++++++++++
		if(dxj == 0.)
		{
			mbs->UO(UO_LVL_warn) << "Warning: dx" << jPar << " is zero. Check, Sensitivity.num_diff_parameter_absolute and Sensitivity.num_diff_parameter_relative.\n";
			dxj = 1.;
		}

		if(f1i.Length() == 0)
		{
			mbs->UO(UO_LVL_err).InstantMessageText("error: no matching sensitivity values from sensor found. Sensitivity analysis stopped.\n");
			return 1; // error
		}
		if (mbs->StopCalculation()) break;
		indexf_max = f2i.Length(); // length of vectors f2i and f1i
		if(jPar == 1)//test && nFactor == 1)
		{
			//test S.SetSize(indexf_max, number_of_params*factors.Length()); // initialize size
			S.SetSize(indexf_max, number_of_params); // initialize size
			S.SetAll(0.);
			//test Smax.SetSize(indexf_max, number_of_params); // initialize size
			//test Smax.SetAll(0.);
		}
		int indexf = 1; // index of vectors f2i and f1i
		for(int i = 1; i <= mbs->NSensors(); i++)
		{
			if(mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Sensitivity.use_final_sensor_values",0) || mbs->GetSensor(i).HasSensorProcessingEvaluationData())
			{
				// compute sensitivity matrix df_i/dparameter_j
				double sensitivity = fabs((f2i(indexf) - f1i(indexf))/dxj);
				S(indexf, jPar) = sensitivity; //case: only one analysis --> Smax, nFactor, factors can be deleted in this case (Smax --> S)!!

				/*test[		S(indexf, jPar + (nFactor-1)*number_of_params) = sensitivity; // multiple sensitivity analysis possible (for tests)
				if(Smax(indexf, jPar) < sensitivity)
				{
				mbs->UO(UO_LVL_0) << "Smax(" << indexf << ", " << jPar << ") = " << sensitivity << "\n";
				Smax(indexf, jPar) = sensitivity;
				}
				test]*/
				mbs->UO(UO_LVL_0) << "sensor: " + mbs->GetSensor(i).GetSensorName() + mystr("\nparameter :") + (*param_names(jPar)) + mystr("\nsensitivity: (") + mystr(f2i(indexf)) + mystr("-") + mystr(f1i(indexf)) + mystr(")/") + mystr(dxj) + mystr("=") + mystr(sensitivity) + mystr("\n");
				mbs->UO(UO_LVL_0) << "\nS = df" << indexf << "/dp" << jPar << " = (" << mystr(f2i(indexf)) + mystr("-") + mystr(f1i(indexf)) + mystr(")/") + mystr(dxj) + mystr("=")  << sensitivity << "\n\n";
				//old: --> output is transposed!!! mbs->SolParFile() <<  mystr(mystr("d_sens") + mystr(indexf) + mystr("/d_par") + mystr(jPar) + mystr("=") + mystr(sensitivity) + mystr("\t"));
				indexf++;
			}
		}
		if (mbs->StopCalculation()) break;
	} // end for
	mbs->SolParFile() <<  "%Sensitivity matrix: d_sensorValue(row)/d_parameter(column)\n";
	mbs->UO(UO_LVL_0) << "%Sensitivity matrix: d_sensorValue(row)/d_parameter(column):\n";

	if(S.GetMatPtr())
	{
		for(int row = 1; row <= sensor_names.Length(); row++)
		{
			for(int col = 1; col <= param_names.Length(); col++)
			{
				mbs->UO(UO_LVL_0) << S(row, col) << "\t";
				mbs->SolParFile() << S(row, col) << "\t";
			}
			mbs->SolParFile() << "\n";
			mbs->UO(UO_LVL_0) << "\n";
		}
	}

  //old:[	if(factors.Length() > 1)
	//{
	//	mbs->SolParFile() << "\n%==========================\n";
	//	mbs->SolParFile() << "%maximal sensitivity values: \n";
	//	for(int jPar = 1; jPar <= number_of_params; jPar++)
	//	{
	//		for(int indexf = 1; indexf <= indexf_max; indexf++)
	//		{
	//				mbs->SolParFile() <<  mystr(mystr("d_sens") + mystr(indexf) + mystr("/d_par") + mystr(jPar) + mystr("=") + mystr(Smax(indexf, jPar)) + mystr("\t"));
	//		}
	//		mbs->SolParFile() <<  "\n";
	//	}
	//}
	//mbs->SolParFile() << "%==========================\n";
	//old:]
	mbs->SolParFile() <<  flush;

	mbs->Get_pCFB()->FinishedComputation();

	return 0;
}
