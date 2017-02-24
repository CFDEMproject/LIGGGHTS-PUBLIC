//#**************************************************************
//#
//# filename:             Optimization.h
//#
//# author:               Rafael Ludwig, Johannes Gerstmayr  
//#
//# generated:						January 2011
//# description:          OPTIMIZATION 
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
 
#ifndef OPTIMIZATION__H
#define OPTIMIZATION__H
#include "mbs.h"
#include "performcomputation.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//b: NewtonOptimization
//$ RL 2011-11-3:[ parent class for newton minima search: f(x)->min (NewtonFunction f(x) can be defined in child class)
// 
// Formulas from "Skriptum zur Vorlesung Numerik und Optimierung für Mechatroniker, Walter Zulehner, Institut für Mathematik, Johannes Kepler Universität Linz, Wintersemester 1994/95, p.86
// dx = (x-x(k))
// f(x) =~ f(x(k))+g(x(k))dx+0.5*dx^t*H*dx = qk(x) --> min.
// grad(qk) = g(x(k))+H*dx, g(x) = [df(x)/dxi]|i
// grad(grad(qk))=H(x(k))=[d^f(x)/(dxi*dxj)]|i,j  ... assumption: H>0
// iteration process
// x(k+1) = x(k) - H(x(k))^(-1)*g(x(k))
// function names of class:
// f ... NewtonFunction
// g ... NewtonFunctionGradient
// H ... NewtonFunctionHesseMatrix


class NewtonOptimization
{
public:
	NewtonOptimization(MultiBodySystem* mbsi)
	{
		abs_diff_val = 0.;
	  rel_diff_val = 0.;
		mbs = mbsi;
		NewtonFunctionCalls = 0;
	}

	virtual void SetNewtonOptimization(double abs_diff_valI, double rel_diff_valI)
	{
		assert(abs_diff_valI > 0 || rel_diff_valI > 0);
		abs_diff_val = abs_diff_valI;
	  rel_diff_val = rel_diff_valI;
	}
	//evaluation of function value f(x)
	virtual double NewtonFunction(const Vector& x); //R^n==>R, 2x differentiable, define this in child class
	
	//dx=K*x + D
	double EvalDx(const Vector& x, const int i);

	// for saving time, some function evaluations are done here in order to use the values multiple times.
	void PrecomputeNewtonFunctionValues(const Vector& x, double& functionValue, TArray<double>& functionValues);

	//g=df/dx, vs...storage of vectors
	Vector NewtonFunctionGradient(const Vector& x, const double& functionValue, const TArray<double>& functionValues);

	// symmetric Hesse matrix
	Matrix NewtonFunctionHesseMatrix(const Vector& x, const double& functionValue, const TArray<double>& functionValues);
	
	//performs one newton iteration step depending on step factor alpha €(0,1), return value is the vector of next iteration x(k+1)
	virtual Vector PerformNewtonIteration(const Vector& xk, const double alpha, double& functionValue);
	
	// do some iterations to reduce functional f(x): optimal value is f(x)=0, x is returned, the residual after the newton iterations is evaluated
	virtual Vector PerformNewtonIterations(const Vector& x, double& newtonFunctionResidual, const int number_of_iterations = 1);

	virtual double NewtonErrorValue(){return 1e100;} // this dummy value is used in case of error (e.g. inversion of jacobian not successful)
	// check the performance 
	void IncreaseNewtonFunctionCalls(){NewtonFunctionCalls++;};
	virtual int GetNewtonFunctionCalls(){return NewtonFunctionCalls;};

private:
	double abs_diff_val; // absolute value "D" for computation of dx=K*x + D (==>f'(x) = df/dx). 
	double rel_diff_val; // relative tolerance for step size	
	MultiBodySystem * mbs;
  int NewtonFunctionCalls; //number of newton function calls (e.g. for checking the performance of an algorithm)
};

class Optimization: public PerformComputation, public NewtonOptimization
{
public:
	//constructor
	Optimization(MultiBodySystem* mbsi):PerformComputation(mbsi), NewtonOptimization(mbsi)
	{
		//mbs = mbsi;
		method="";
		number_of_params=0;
		initial_population_size=0;
		surviving_population_size=0;
		number_of_children=0;
		number_of_generations=0;
		range_reduction_factor=0.;
		n_comp=0; //number of computations (e.g.: children or initial_population_size)
		param_minvals(0);
		param_maxvals(0);
		actparameters(100);
		actcost_function_values(0);
		parent_parameters(100);
		parent_parameters.SetAll(0);
		parent_cost_function_values(0);
		isRestart = 0;
		iparfile = 0;
		edc = 0;
	};
	//destructor
	~Optimization()
	{
		for(int i=1; i<=actparameters.Length(); i++)
		{
			//delete actparameters(i);
			actparameters(i) = 0;
		}				
	}
	//initialize Optimization
	void Init();

	// initialize parameters of first generation or mutate surviving parameters (ga)
	void InitializeFirstGenerationOrMutateSurvivors(const int gen, const int child, const int np,const Matrix& rand_vals, Vector* paramset);

	//evaluate optimization cost function based on sensors
	void EvalCostFunction(double& cost_function_val);
  //return distance to nearest parameter of parameter "candidate"
	double GetDistanceToNearestParameter(Vector* candidate);
  // radius of (hyper-)sphere in normed parameter space; in the inner of this radius, only one parameter survives (helpful in case of more than one (local) minima)
	double GetNearestAllowedDistance(const int gen);
  // select the surviving parameters
	void ChooseSurvivingPopulation(const int gen);
	// read a line from parameter file and stores values in vals, returns 1 if done
	int ReadOptValuesFromSolParFile(MBS* mbs, TArray<double>& vals, int& gen);
	// read initial population from file
	void DoRestart(TArray<Vector*>& allparameters, int& gen);
	// compute first generation or one further generation based on parent parametes, return value is 0 if successful
	int ComputeOneGeneration(int& gen, Vector& nominal_parameters, TArray<Vector*>& allparameters, int run_with_nominal_parameters);
	// select optimization type in this function
	int PerformOptimization();
  // main part of genetic optimization
	int PerformGeneticOptimization();
	//***************************************************
	// main part of newton optimization
	int PerformNewtonOptimization();

	virtual double NewtonFunction(const Vector& x); // from class NewtonOptimization
	double PenaltyFunctionValue(const Vector& x) const; // additional penalty value for NewtonFunction is computed, in case of components of x are not in the correct limits
private:
	ElementDataContainer* edc;  // solver options edc
	mystr method; // e.g. "Genetic", "Newton"
	int number_of_params; // number of parameters, which are optimized
	int initial_population_size; // initial parameter population size
	int surviving_population_size; // number of surviving parameters each generation
	int number_of_children; // number of child parameters of the surviving generation
	int number_of_generations; //max. number of parameter generations
	double range_reduction_factor;	 // reduction factor, which reduces the distance between parent parameters and child parameters each generation
	int n_comp; //number of computations (children or initial_population_size)
	TArray<double> param_minvals; // lower limit of parameter intervals
	TArray<double> param_maxvals; // upper limit of parameter intervals
	//TArray<mystr*> param_names;   // names of optimized parameters (same as edc-name)
	//store values of actual generation: 
	TArray<Vector*> actparameters; // actual generation parameters
	TArray<double> actcost_function_values; // actual cost function values
	TArray<Vector*> parent_parameters; // parend (surviving) parameters
	TArray<double> parent_cost_function_values; // cost function values due to sensor computation (sensor value is error, which is minimized)
	int isRestart; // set to 1 if already existing parameters (in parameter file) should be used as initial population
	CMFile* iparfile; //parameterfile for reading
};
#endif
