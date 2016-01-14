//#**************************************************************
//#
//# filename:             performcomputation.h
//#
//# author:               Rafael Ludwig
//#
//# generated:						February 2011
//# description:          class PerformComputation created
//#
//# remarks:              functions used in different computation types (e.g.: Optimization) are in this class to avoid duplicate code
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
 
#ifndef PERFORMCOMPUTATION__H
#define PERFORMCOMPUTATION__H
#include "MBS.h"
// RL 2011-02: use this class for functions used in multiple perform computations types --> no duplicate code!
class PerformComputation
{
public:
	PerformComputation(MultiBodySystem* mbsi)
	{
		mbs = mbsi;
		param_names(0);   // names of parameters for sensitivity analysis
	}
	// destructor
	~PerformComputation()
	{
		for (int i=1; i <= param_names.Length(); i++)
		{
			delete param_names(i);
			param_names(i) = 0;
		}
		//param_elementnumbers.SetLen(0);
	}
	// get model parameter vector
	virtual const Vector GetModelParameters() const
	{
		int number_of_params = param_names.Length();
		Vector modelParams(number_of_params);
		mbs->UO(UO_LVL_multsim) << "\n********************\n";
		for (int col=1; col <= number_of_params && !mbs->StopCalculation(); col++)
		{
			mystr varname = *param_names(col);
			if (varname != mystr(""))
			{
				if(!mbs->GetModelDataContainer())
				{
					mbs->UO(UO_LVL_err) << "\nError: no model parameters found!\n*********************\n";
					modelParams(col) = 0.;
				}
				else
				{
					ElementData* elfound = mbs->GetModelDataContainer()->TreeFind(varname);
					if(!elfound)
					{
						mbs->UO(UO_LVL_err) << "\nError: can not access variable '" << varname << "' for parameter variation!\n*********************\n";
						modelParams(col) = 0.;
					}
					else
					{
						modelParams(col) = elfound->GetDouble();
					}
				}
				mbs->UO(UO_LVL_multsim) << "    " << varname << " = " << modelParams(col) << "\n";
			}
		}
		return modelParams;
	}
	
	// replace mbs-parameter or model parameter with parameter val
	virtual void SetParameter(mystr varname, double varparameter)// , int element_number=0)
	{
		if (varname != mystr(""))
		{
			if (mbs->GetMBS_EDC_Options()->TreeFind(varname))
			{
				mbs->GetMBS_EDC_Options()->TreeFind(varname)->SetLocked(0); //$ RE, MSax, DR 2013-06-27: set parameters for script language
				mbs->GetMBS_EDC_Options()->TreeSetDouble(varname, varparameter);
				mbs->GetMBS_EDC_Options()->TreeFind(varname)->SetLocked(1); //$ RE, MSax, DR 2013-06-27: set parameters for script language
			}
			else if (mbs->GetModelDataContainer() != 0 && mbs->GetModelDataContainer()->TreeFind(varname))
			{
				mbs->GetModelDataContainer()->TreeFind(varname)->SetLocked(0); //$ RE, MSax, DR 2013-06-27: set parameters for script language
				mbs->GetModelDataContainer()->TreeSetDouble(varname, varparameter);
				mbs->GetModelDataContainer()->TreeFind(varname)->SetLocked(1); //$ RE, MSax, DR 2013-06-27: set parameters for script language
			} 
			else
			{
				mbs->UO(UO_LVL_err) << "*********************\nError: can not access variable '" << varname << "' for parameter variation!\n*********************\n";
			}
		}
	}
	
	// replace mbs-parameter or model parameters with parameters from vector paramset
	virtual void SetParameters(const Vector* paramset)
	{

		mbs->UO(UO_LVL_multsim) << "\n********************\n";
		for (int col=1; col <= paramset->Length() && !mbs->StopCalculation(); col++)
		{
			mystr varname = *param_names(col);
			SetParameter(varname, (*paramset)(col));//, param_elementnumbers(col));
			mbs->UO(UO_LVL_multsim) << "    " << varname << " = " << (*paramset)(col) << "\n";
		}
	}

  virtual TArray<mystr*>& GetModelParameterNames()
	{
		return param_names;
	}
public:
	MultiBodySystem* mbs;
	TArray<mystr*> param_names; // names of parameters
	//TArray<int> param_elementnumbers;
};
#endif