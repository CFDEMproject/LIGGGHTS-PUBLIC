//#**************************************************************
//#
//# filename:             optimization.cpp
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
 
#include "mbs.h"
#include "myfile.h"
#include "optimization.h"

//$ YV 2013-01-02: the variable below belongs to the elements and models dll, and is now defined in sensor.cpp
//int FFToptimization_averaging=0; //0=no averaging, 1, 2, n.... is average over -n ... +n values (0=no averaging, 1=2 values, 2=5values)
int opt_simulation_number = 0; // this number is increased every optimization; in order to be able to distinguish the solutions of sensors
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//b: NewtonOptimization
double NewtonOptimization::NewtonFunction(const Vector& x) //f(x): x€R^n==>R, 2x differentiable
{
	IncreaseNewtonFunctionCalls();
	double x1 = x(1);
	double x2 = x(2);
	return Sqr((x1-0.15)*(x2-0.2))+To4(x1-0.15);
	//return To4(sqrt(fabs(x1))-0.25)+Sqr((sqrt(fabs(x1))-0.25)*(sqrt(fabs(x2))-0.09));
	//return To4(sqrt(fabs(x1))-0.25)+Sqr((sqrt(fabs(x1))-0.25)*(sqrt(fabs(x2))-0.09)); //[0.0625, 0.0081]; [0.06270733203465699, 0.0084608300288260049], f12x1opt=7.043346855914265e-013 after 13 iterations
	//return Sqr(x1)+Sqr(x2); // [0,0];   x0 = [-0.0001002997743224654, -0.0001002997743224654], f0xopt=2.012008945827498e-008 after 1 iterations  
}; // ... 

//dx=K*x + D
double NewtonOptimization::EvalDx(const Vector& x, const int i)
{
	return rel_diff_val*x(i)  + abs_diff_val;
}

//pre-compute function values, which are multiple used and store them; functionValue ... f(x), functionValues(j) ... f(x1,x2,...,xi+dxi, ... xn)
void NewtonOptimization::PrecomputeNewtonFunctionValues(const Vector& x, double& functionValue, TArray<double>& functionValues)
{
	functionValue = NewtonFunction(x); //f(x)
	
	int len = x.Length();
	for(int i = 1; i <= len; i++)
	{
		Vector xi = x;
		xi(i) += EvalDx(x,i); 
		//++++++++++++++++++++++++++++++++++++
		double fi = NewtonFunction(xi); //f(x1,x2,...,xi+dxi, ... xn)
		if(i == 1)
		{
			functionValues.SetLen(0);
			if(mbs->UO(UO_LVL_dbg2).PrintMsg())
				mbs->UO(UO_LVL_dbg2, 16) << "f(x) = " << functionValue << "\n";			
		}
		functionValues.Add(fi);
		if(mbs->UO(UO_LVL_dbg2).PrintMsg())
			mbs->UO(UO_LVL_dbg2, 16) << "f(x" << i << ") = " << fi << "\n";
	}
}

//g=df/dx, vs...storage of vectors
Vector NewtonOptimization::NewtonFunctionGradient(const Vector& x, const double& functionValue, const TArray<double>& functionValues)
{
	//if(0)
	//{
	//	return 2.*x(1); //test
	//}
	//else
	//{
		int len = x.Length();
		Vector gvec(len);

		int type_diff = 1; //1 ... forward difference quotient
		for(int i = 1; i <= len; i++)
		{
			//++++++++++++++++++++++++++++++++++++
			//b: difference quotient variables
			double dxi =  EvalDx(x, i);
			double f2 = functionValues(i);
			double f1 = functionValue;
			gvec(i) = (f2-f1)/dxi;
		}
		return gvec;
	//}
}

// symmetric Hesse matrix
Matrix NewtonOptimization::NewtonFunctionHesseMatrix(const Vector& x, const double& functionValue, const TArray<double>& functionValues)
{
	// grad(grad(f))=~grad(grad(qk)):=H(x(k))
	// H(x(k))=[d^f(x)/(dxi*dxj)]|i,j  ... assumption: H>0
	//if(0)
	//{
	//	return Matrix(2.); //test
	//}
	//else
	//{
		// create symmetric Hesse matrix
		//
		// H(i,j)=1/(dxi*dxj)*(f(x1,x2,...,xi+dxi,...,xj+dxj,...,xn)-f(x1,x2,...,xi+dxi,...,xn)-f(x1,x2,...xj+dxj,...,xn)+f(x1,x2,...,xn))
		//       =1/(dxi*dxj)* (fij - fi - fj + f )
		// case: i==j
		// H(i,i)=1/(dxi*dxi)*(f(x1,x2,...,xi+2*dxi,..,xn)-2*f(x1,x2,...xi+dxi,...,xn)+f(x1,x2,...,xn))
		// add stored vectors to Hesse matrix
		int len = x.Length();
		Matrix Hmat(len,len);
		for(int i = 1; i <= len; i++)
		{
			for(int j = i; j <= len; j++)
			{
				//    f  - fi - fj
				Hmat(i,j)=functionValue-functionValues(i)-functionValues(j);
				//if(mbs->UO(UO_lvl_dbg2).PrintMsg())
				//{
				//	mbs->UO(UO_lvl_dbg2, 16) << "f - f" << i << "-f" << j << ")=" << Hmat(i,j) << "\n";
				//}
			}
		}

		//compute derivatives f(x1,x2,...,xi+dxi,...,xj+dxj,...,xn)
		for(int i = 1; i <= len; i++)
		{
			for(int j = i; j <= len; j++)
			{
				Vector xij = x;	
				double dxi = EvalDx(x, i);
				double dxj = EvalDx(x, j);
				xij(i) += dxi;
				xij(j) += dxj;
				double fval = NewtonFunction(xij);


				Hmat(i,j) = Hmat(i,j)+fval;
				//if(mbs->UO(UO_lvl_dbg2).PrintMsg())
				//{
				//	mbs->UO(UO_lvl_dbg2,16) << "f" << i << j << "=" << fval << "\n";
				//	mbs->UO(UO_lvl_dbg2,16) << "sum f" << i << j << "=" << Hmat(i,j) << "\n";
				//}				
				Hmat(i,j) = Hmat(i,j)/(dxi*dxj);
				if(mbs->UO(UO_LVL_dbg2).PrintMsg())
				{
					mbs->UO(UO_LVL_dbg2,16) << "H" << i << j << "=" << Hmat(i,j) << "\n";
				}				
				

				if(i!=j)
				{
					Hmat(j,i) = Hmat(i,j); // symmetric Hesse matrix
				}
			}
		}
		return Hmat;
	//}
}

//performs one newton iteration step depending on step factor alpha €(0,1), return value is the vector of next iteration x(k+1)
Vector NewtonOptimization::PerformNewtonIteration(const Vector& xk, const double alpha, double& functionValue)
{
	TArray<double> functionValues;
	PrecomputeNewtonFunctionValues(xk, functionValue, functionValues);
	// compute gradient and Hesse matrix by use of precomputed values
	Vector gk = NewtonFunctionGradient(xk,functionValue, functionValues);
	Matrix invHk = NewtonFunctionHesseMatrix(xk, functionValue, functionValues);
	if(mbs->UO(UO_LVL_dbg1).PrintMsg())
	{
		mbs->UO(UO_LVL_dbg1, 16) << "    Hk=" << invHk << "\n";
		mbs->UO(UO_LVL_dbg1, 16) << "    gk=" << gk << "\n";
	}

	if(mbs->UO(UO_LVL_sim).PrintMsg())
	{
		Vector ev;
		invHk.Eigenvalues(ev);
		double max_ev = 0.;
		double min_ev = 1e100;
		for(int i=1;i<=ev.Length();i++)
		{
			max_ev = max(max_ev, ev(i));
			min_ev = min(min_ev, ev(i));
		}
		if(max_ev > 0. && min_ev < 1e100)
		{
			mbs->UO(UO_LVL_multsim) << "H = " << invHk << "\n";
			double cond_H = abs(max_ev)/abs(min_ev);
			mbs->UO(UO_LVL_multsim) << "cond(H) = " << cond_H << " (if this number high, small (e.g. numerical) errors have big influence on the iteration step size.)\n";
			mbs->UO(UO_LVL_multsim) << "det(H) = " << invHk.Det() << " (determinant must not vanish ==> otherwise solution does not exist.)\n";
		}
	}


	
	// inverse exists
	Vector dx;
	if(!invHk.Invert())
	{
		// should not happen
		mbs->UO(UO_LVL_multsim) << "WARNING: Jacobian inversion NOT successful!\n";// , coefficients of parameter vector of dx set to " << NewtonErrorValue() << "!";
		dx.SetLen(xk.Length());
		dx.SetAll(NewtonErrorValue());
		return dx;
	}
	else
	{
		Mult(invHk, gk, dx);
		if(mbs->UO(UO_LVL_dbg1).PrintMsg())
		{
			mbs->UO(UO_LVL_dbg1, 16) << "    dx=" << dx << "\n";
		}
		return xk - alpha*dx; // x(k) --> x(k+1): x(k+1) = x(k) - H(x(k))^(-1)*g(x(k))
	}
}

// do some iterations to reduce functional f(x): optimal value is f(x)=0, x is returned
Vector NewtonOptimization::PerformNewtonIterations(const Vector& x, double& newtonFunctionResidual, const int number_of_iterations)
{
	int Nit = 0;
	int kmax = 1;
	Vector xopt = x; //xopt...optimal x-value with lowest f(xopt)
	Vector xk = x;   //xk ... x-vector of iteration k
	double fopt = 1e100;
	for(int k=0;k<=number_of_iterations;k++)
	{
		// storage of evaluated function values during one iteration (needed if f(x) is very time expensive)	     
		double alpha = 1.0; // step size factor, Newton formula: alpha == 1 (see script "zulehner")
		int maxit = 5; // maximum of iterations with decreased step size (for search of f(x(k+1)<f(xk))
		int stop=1;    // stop iterations of loop with iterator 'k', if xopt == x (no better solution than starting value found)
		for(int it = 1; it<=maxit; it++)
		{				
			Nit++;
			double functionValue; // function value of starting value f(xk)
			Vector xact = PerformNewtonIteration(xk, alpha, functionValue);
			if(k==0&& it==1) 
			{
				// store start values of newton: xk, f(xk)
				fopt = functionValue;
				xopt=xk;
				if(mbs->UO(UO_LVL_multsim).PrintMsg()){mbs->UO(UO_LVL_multsim, 16) << "x0" << " = " << xk << ", f0" << "=" << fopt << "\n";}
			}
			// check, if new function value is smaller than other function values
			double f_act=NewtonErrorValue();
			if(xact.Get(1) != NewtonErrorValue())
			{
								// successful computation of NewtonIteration-Step
				 f_act = NewtonFunction(xact); // evaluate new function value for comparison
			}
			//if(mbs->UO(UO_LVL_multsim).PrintMsg()){	mbs->UO(UO_LVL_multsim, 16) << "k = " << k <<", it=" << it << ",xact = " << xact << ", f_act=" << f_act << ", alpha = " << alpha << "\n";	}
				
			if(f_act < fopt)
			{
				fopt = f_act;
				xopt = xact;
				if(mbs->GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.Optimization.Newton.absolute_accuracy") < fopt)
				{
					stop = 0;
				}
				break;
			}
			if (mbs->StopCalculation()) break;
			//mbs->UO(UO_LVL_dbg1, 16) << "  alpha = " << alpha << ", f_act=" << fnext << "\n";			
			alpha *= 0.25	;	
		}// alpha - loop
	
		xk=xopt;
		if(mbs->UO(UO_LVL_multsim).PrintMsg()){mbs->UO(UO_LVL_multsim, 16) << "x" << k+1 << " = " << xk << ", f"<< k+1 << "=" << fopt << "\n";}
		if(stop)
		{
			if(mbs->UO(UO_LVL_0).PrintMsg()){mbs->UO(UO_LVL_0, 16) << "++++++ NEWTONS method finished after " << Nit << " iterations! +++++++\n\n";}
			break;
		}
		if (mbs->StopCalculation()) break;
	}
	newtonFunctionResidual = fopt;
	return xopt;
}
//e: NewtonOptimization
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




//initialize Optimization
void Optimization::Init()
{
	edc = PerformComputation::mbs->GetMBS_EDC_Options();
	//init optimization data
	method = edc->TreeGetString("SolverOptions.Optimization.method");
	number_of_params = edc->TreeGetInt("SolverOptions.Optimization.Parameters.number_of_params");
	initial_population_size = edc->TreeGetInt("SolverOptions.Optimization." + method + ".initial_population_size");
	surviving_population_size = edc->TreeGetInt("SolverOptions.Optimization." + method + ".surviving_population_size");
	number_of_children = edc->TreeGetInt("SolverOptions.Optimization.Genetic.number_of_children");
	number_of_generations = edc->TreeGetInt("SolverOptions.Optimization.Genetic.number_of_generations");
	range_reduction_factor = edc->TreeGetDouble("SolverOptions.Optimization.Genetic.range_reduction_factor");
	if (surviving_population_size > initial_population_size) initial_population_size = surviving_population_size;

	// init random
	srand((int)(edc->TreeGetDouble("SolverOptions.Optimization.Genetic.randomizer_initialization")*(double)RAND_MAX));
	// init solution file
	isRestart = edc->TreeGetInt("SolverOptions.Optimization.restart",0);

	actparameters(100);

	parent_parameters(100);
	if(isRestart)
	{
		PerformComputation::mbs->OpenFiles(1); //open files already here, append to existing file
		mystr dir = PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetString("GeneralOptions.Paths.sensor_output_path"); //old: GetTOption80);
		mystr solparn = PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetString("SolverOptions.Solution.ParameterFile.parameter_variation_filename"); //GetTOption82);
		mystr file = dir + solparn;
		iparfile = new CMFile(file.c_str(), TFMread);
		if(!iparfile->IsGood())
		{
			PerformComputation::mbs->UO(UO_LVL_multsim) << mystr(mystr("Optimization: Could not open \"") + file + mystr("\". New parameter file is created by optimization."));
			PerformComputation::mbs->CloseFiles();
			isRestart = 0; //abort reading from file and create new parameter file
			delete iparfile;
		}
		else
		{
			mystr line;
			mystr word;
			int nrOfHeaderLines = 3;
			for(int i=1;i<=nrOfHeaderLines;i++)
			{
				iparfile->RWuntilEOL(line); // header
			}
			int pos = 0;
			word = line.GetWord(pos,0); // first word of this line must be "%generation" --> otherwise, the optimization will be started in the normal modus
			if(!word.Compare(mystr("%generation")))
			{
				PerformComputation::mbs->UO(UO_LVL_multsim) << "Optimization: restart was not possible. New parameter file is created by optimization.";
				PerformComputation::mbs->CloseFiles();
				isRestart = 0; //abort reading from file and create new parameter file
				delete iparfile;
			}
		}
	}

	if(!isRestart)
	{
		//attention: change also nrOfHeaderLines, if header lines are inserted of removed
    int flag_append=!PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.always_replace_files");
		PerformComputation::mbs->OpenFiles(0); //open files already here, do not append to any file
		PerformComputation::mbs->SolParFile() << "%Optimization results file\n";
		PerformComputation::mbs->SolParFile() << "%\n";
		PerformComputation::mbs->SolParFile() << "%generation child n ";
	}

	for (int col=1; col <= number_of_params; col++)
	{
		mystr str =	edc->TreeGetString(mystr("SolverOptions.Optimization.Parameters.param_name")+mystr(col));
		param_names(col) = new mystr(str);
		//param_elementnumbers(col) = edc->TreeGetInt(mystr("SolverOptions.Optimization.Parameters.element_number")+mystr(col));
		param_minvals(col) = edc->TreeGetDouble(mystr("SolverOptions.Optimization.Parameters.param_minval")+mystr(col));
		param_maxvals(col) = edc->TreeGetDouble(mystr("SolverOptions.Optimization.Parameters.param_maxval")+mystr(col));
		if(!isRestart) PerformComputation::mbs->SolParFile() << str << " ";

		if(param_minvals(col) >= param_maxvals(col))
		{
			PerformComputation::mbs->UO().InstantMessageText(mystr("WARNING: SolverOptions.Optimization.Parameters.param_minval")+mystr(col)+mystr(" bigger or equal to ")+mystr("SolverOptions.Optimization.Parameters.param_maxval")+mystr(col));
		}
	}
	
	if(!isRestart) PerformComputation::mbs->SolParFile() << "cost_function_value\n";// thread\n";
}


// initialize parameters of first generation or mutate surviving parameters (ga)
void Optimization::InitializeFirstGenerationOrMutateSurvivors(const int gen, const int child, const int np,const Matrix& rand_vals, Vector* paramset)//, TArray<mystr*>& str_file)
{
	paramset->SetLen(number_of_params);
	for (int col=1; col <= number_of_params && !PerformComputation::mbs->StopCalculation(); col++)
	{
		// old: randval = (double)rand()/(double)RAND_MAX; //value between 0 and 1;  ... bad for multiple processing
		double randval = rand_vals(np,col)/(double)RAND_MAX;// value between 0 and 1;
		double range = (param_maxvals(col) - param_minvals(col));
		double start = param_minvals(col);

		if (gen == 1)
		{
			//in first generation, generate random values over total range
			(*paramset)(col) = randval * range + start;
		}
		else
		{
			//in further generations, mutate recent populations
			randval -= 0.5; //now goes from -0.5 .. 0.5
			int use_gauss_dist = 1; //use gaussian distribution
			if (use_gauss_dist)
			{
				double val = 2.*fabs(randval); // 0..1
				double f = 1e10; //maximum number for case of val == 0
				if (val > 0 && val <= 1)
				{
					//x=+-sqrt(-log(y));  //Gauss function
					f = sqrt(-log(val));  // val=]0,1] -->f(val)=]oo,0]
				}
				else if (val > 1)
				{
					f = 1;
				}
				//variance = sqrt(.5)
				//f(x)=exp(-x^2), probabilities: f € [-2,2] ==> 99 percent; f € [-1,1] ==> 84 percent; f € [-0.5,0.5] ==> 52 percent
				randval = 0.5*Sgn(randval)*f;  // randval € [-1,1] ==> 99 percent; randval € [-0.5,0.5] ==> 84 percent; randval € [-0.25,0.25] ==> 52 percent
			}

			(*paramset)(col) = parent_parameters(child)->Get(col) + randval *	range	*	pow(range_reduction_factor, (double)gen-1.);
			if ((*paramset)(col) < start) (*paramset)(col) = start;
			if ((*paramset)(col) > start + range) (*paramset)(col) = start + range;
		}
	}
}

//evaluate optimization cost function based on sensors
void Optimization::EvalCostFunction(double& cost_function_val)
{
	cost_function_val = PerformComputation::mbs->EvalSensorCostFunctionVal();
  //$ RL 2011-02: moved to MBS
	//for (int i=1; i <= mbs->NSensors(); i++)
	//{
	//	if (mbs->GetSensor(i).GetWriteResults() && mbs->GetSensor(i).FlagSensorComputation() != 0)
	//	{
	//		Vector v = mbs->GetSensor(i).GetSensorComputationValue();

	//		cost_function_val += v.GetNorm(); //L2-norm of all values (usually fabs of single value)
	//	}
	//}
	
	actcost_function_values.Add(cost_function_val);
}
//return distance to nearest parameter of parameter "candidate"
double Optimization::GetDistanceToNearestParameter(Vector* candidate)
{
	double min_dist = 1e100; //initialize to oo
	// compute distance to nearest parent
	for (int iParent=1; iParent<=parent_parameters.Length(); iParent++)
	{
		Vector* parent = parent_parameters(iParent);
		// compute normed distance to one parent
		double dist = 0.;
		for (int col=1; col<=number_of_params; col++)
		{
			double range = param_maxvals(col) - param_minvals(col);
			double delta = candidate->Get(col) - parent->Get(col);
			if(range!=0.){dist += Sqr(delta/range);} // parameter range=0 would give infinite distance ==> these parameters omitted for calculation of distance
		}
		dist = sqrt(dist); // sqrt(sum(delta²(col)/range²(col),col=col...nParameters))
		min_dist = min(min_dist,dist);
	}
	return min_dist;
}
// radius of (hyper-)sphere in normed parameter space; in the inner of this radius, only one parameter survives (helpful in case of more than one (local) minima)
double Optimization::GetNearestAllowedDistance(int gen)
{
	//                   dimension-times e.g. 2D ==> sqrt(2)
	double val = edc->TreeGetDouble("SolverOptions.Optimization.Genetic.min_allowed_distance_factor", 0.0); // >=0, 1 = normed length of side of parameter-hyper-cube

	val *= pow(range_reduction_factor, (double)gen-1.); // depends on generation

	return val;
}
// select the surviving parameters
void Optimization::ChooseSurvivingPopulation(const int gen)
{
	TArray<int> surviving_ind(actcost_function_values.Length());
	for (int i=1; i<=actcost_function_values.Length(); i++)
	{
		surviving_ind(i) = i;
	}

	//sort the cost function values ==> get best ones!
	TArray<double> sortedcostfn = actcost_function_values;
	QuicksortDouble(sortedcostfn, surviving_ind);
	parent_parameters.SetAll(0);
	parent_parameters.SetLen(0);
	parent_cost_function_values.SetAll(0);
	parent_cost_function_values.SetLen(0);

	//use the first "surviving_population_size" values as the best values
	
	PerformComputation::mbs->UO(UO_LVL_multsim) << "----------\nSelection:\n";
	int nsurvivors = 0;
	for (int i=1; i<=surviving_ind.Length(); i++)
	{
		Vector* candidate = actparameters(surviving_ind(i));
		//mbs->UO() << "Candidate: " << mystr(i) << mystr(" ") << candidate->Get(1) << ", Par2 = "  <<  candidate->Get(2) << mystr("\n");
		if(i <= surviving_ind.Length() && // if stop button is pressed during simulation
			parent_parameters.Length() < surviving_population_size) // not all parents found
		{
			Vector* candidate = actparameters(surviving_ind(i));
			if(GetDistanceToNearestParameter(candidate) > GetNearestAllowedDistance(gen))
			{// distance of parameter is not too close to already chosen parents (with better costfunction))
				parent_parameters.Add(candidate);
				parent_cost_function_values.Add(sortedcostfn(i));
				nsurvivors++;				
			}
			else
			{
				if(nsurvivors<surviving_population_size)
				{
					PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("Parameter set with cost-function ranking ")+mystr(i) + mystr(" is too close to other survivors. Reduce the minimal allowed distance, if it should be considered as parent parameter of the next generation.\n");				
				}
			}
		}
	}
	
	if(PerformComputation::mbs->UO(UO_LVL_multsim).PrintMsg())
	{			
		mystr str("The chosen surviving parameters (parent parameters) from generation ");
		str = str + mystr(gen) + mystr(" are:\n");
		PerformComputation::mbs->UO(UO_LVL_multsim) << str << "\n";
		str = mystr("");
		for(int nParent=1;nParent<=surviving_population_size&&nParent<=parent_parameters.Length();nParent++)
		{
			if(parent_parameters(nParent))
			{
				for(int col=1;col<=parent_parameters(nParent)->Length();col++)
				{
					str = str + mystr(parent_parameters(nParent)->Get(col)) + mystr("; ");
				}
				str = str + mystr(" cost function=") + mystr(parent_cost_function_values(nParent)) + mystr("\n");
			}
		}
		PerformComputation::mbs->UO(UO_LVL_multsim) << str << "\n";			
		PerformComputation::mbs->UO(UO_LVL_multsim) << "End of Selection\n----------\n";
	}
}

// read a line from parameter file and stores values in vals, returns 1 if done
int Optimization::ReadOptValuesFromSolParFile(MBS* mbs, TArray<double>& vals, int& gen)
{
	vals.SetLen(0);
	if(iparfile->IsGood())
	{
		mystr line;
		iparfile->RWuntilEOL(line);
		if(line.Length() == 0)
		{
			return 1; //done
		}
		mystr word;
		int pos = 0;
		word = line.GetWord(pos,0); // gen
		if(pos!=-1)gen = word.MakeInt();
		word = line.GetWord(pos,0); // child
		word = line.GetWord(pos,0); // np
		word = line.GetWord(pos,0); // first parameter
		while(pos!=-1)
		{
			//parameters and costfunctionvalue
			vals.Add(word.MakeDouble());
			word = line.GetWord(pos,0);
		}
		if(vals.Length() != number_of_params + 1)
		{
			return 1; //done
		}
		return 0; //success
	}

	return 1; //done

}

void Optimization::DoRestart(TArray<Vector*>& allparameters, int& gen)
{
	n_comp=1; // number of computed parameter
	int done = 0; // end of file
	while(!done) // "candidates" for next generation
	{
		TArray<double> vals;
		done = ReadOptValuesFromSolParFile(PerformComputation::mbs, vals, gen);
		if(!done)
		{
			Vector* paramset = new Vector(number_of_params);
			mystr str_file("");
			for(int col=1;col<=number_of_params;col++)
			{
				double val = vals.Get(col);
				(*paramset)(col) = val; // copy values
				str_file = str_file + val + mystr(" ");
			}
			// read
			allparameters.Add(paramset);
			actparameters.Add(paramset);
			actcost_function_values.Add(vals.Last());
			if(PerformComputation::mbs->UO(UO_LVL_multsim).PrintMsg())
			{
				PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("initial population from file: generation = " + mystr(gen) + "\n  parameter set: " + str_file + " cost function=" + mystr(actcost_function_values.Last()) + mystr("\n\n"));// possibly problematic during multi-processor computation ==> directly output in file
			}
			n_comp++; //counter of paramers
		}
	} // "candidates" for next generation
	delete iparfile;
}

// main part of genetic optimization
int Optimization::PerformOptimization()
{	
	Init(); // initialize variables
	PerformComputation::mbs->UO(UO_LVL_ext) << "\n*****************************************\nmethod = " + method +   "\n*****************************************\n\n";
	//if(method.Compare("Newton"))
	//{
	//	return PerformNewtonOptimization(mbs);
	//}
	//else if(method.Compare("Genetic"))
	if(method.Compare("Genetic"))
	{
		return PerformGeneticOptimization();
	}
	else
	{
		//mbs->UO().InstantMessageText("ERROR: No optimization method chosen!\nmethods=[\"Newton\" | \"Genetic\"]");
		PerformComputation::mbs->UO().InstantMessageText("ERROR: Please choose optimization method \"Genetic\".");
		return 0;
	}
}


int Optimization::ComputeOneGeneration(int& gen, Vector& nominal_parameters, TArray<Vector*>& allparameters, int run_with_nominal_parameters)
{
	// computation of generation bigger than 1 and initialisation in case of no restart with initial population from parameter file
	actparameters.SetLen(0);
	actcost_function_values.SetLen(0);

	int n_child = surviving_population_size;
	if (gen == 1) //first population has special size
	{
		n_child = 1;
	}
	else
	{
		n_child = parent_parameters.Length();

		// consider parent in actual population (if best parameter set is already found, it survives)
		actparameters.SetLen(n_child);
		actcost_function_values.SetLen(n_child);
		for(int child=1; child <= n_child; child++)
		{
			//(RL) parents are also part of the generation
			actparameters.Add(parent_parameters(child));
			actcost_function_values.Add(parent_cost_function_values(child));
		}
	}

	for (int child=1; child <= n_child; child++)
	{
		n_comp = number_of_children; //number of computations
		if (gen == 1) //first population has special size
		{
			n_comp = initial_population_size;
		}

		if(PerformComputation::mbs->GetModelDataContainer()->TreeFind("Optimization.run_once"))
		{
			PerformComputation::mbs->UO(UO_LVL_err).InstantMessageText("Error in Optimization: Old variable name \"Optimization.run_once\" is deprecated! Use new variable \"SolverOptions.Optimization.run_with_nominal_parameters\"\n");
			return 1;
		}

		if(run_with_nominal_parameters)
		{
			//only single computation
			number_of_generations = 1;
			n_child = 1;
			n_comp = 1;
			PerformComputation::mbs->UO(UO_LVL_multsim) << "Run with nominal parameters is activated  (test of nominal optimization cost function value)!\n";
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// compute random variables before (parallel) computation
		Matrix rand_vals(n_comp, number_of_params); // for parallelization
		for (int np=1; np <= n_comp; np++)
		{
			for (int col=1; col <= number_of_params; col++)
			{
				rand_vals(np,col) = (double)rand();
			}
		}

		for (int np=1; np <= n_comp; np++) // "candidates" for next generation
		{		
			Vector* paramset = new Vector(number_of_params);
			allparameters.Add(paramset);
			//+++++++++++++++++++++++++++++++++++++
			//set random values:
			InitializeFirstGenerationOrMutateSurvivors(gen, child, np, rand_vals, paramset);

			// use nominal parameters when following flag is activated
			if(run_with_nominal_parameters)
			{
				PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("\nUse nominal parameters: \n\n");
				paramset = &nominal_parameters; // use nominal parameters for test purpose
			}
			//+++++++++++++++++++++++++++++++++++++
			actparameters.Add(paramset);
		} // "candidates" for next generation

		//omp_set_num_threads(1); // single thread (parallel mbs - computation not possible yet)
		//#pragma omp parallel for
		// compute costfunctions of "candidates" for next generation
		for (int np=1; np <= n_comp && !PerformComputation::mbs->StopCalculation(); np++)
		{
			Vector* paramset = actparameters(actparameters.Length() - n_comp + np);//actparameters((child-1)*n_child + np);
			//+++++++++++++++++++++++++++++++++++++
			// write "candidates" for next generation in model edc
			SetParameters(paramset); // set different parameters only in case of run_with_nominal_parameters = 0 !!!

			//+++++++++++++++++++++++++++++++++++++
			PerformComputation::mbs->CloseFiles(); //because it is opened in following function
			//********************************************************************
			
			//repeatedly call computation:
			PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("\nStart Simulation: \n\n");
			opt_simulation_number++;
			if(!PerformComputation::mbs->RepeatedlyPerformComputation(1)) //do not write into solpar file
			{
				PerformComputation::mbs->UO(UO_LVL_err) << "Error: Optimization stopped. Check consistency of system!\n";
				return 1;
			}


			//compute optimization cost function based on sensors:
			double cost_function_val = 0.;
			EvalCostFunction(cost_function_val);
			//********************************************************************
			int flag_append=!PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Solution.always_replace_files");
			if(!flag_append) flag_append += 2; // append only parameter file
			if(run_with_nominal_parameters)
			{
				flag_append =1+2; // don't overwrite files
			}
			PerformComputation::mbs->OpenFiles(flag_append);
			//********************************************************************

			mystr str_file("");
			for(int col=1;col<=paramset->Length();col++)
			{
				str_file = str_file + paramset->Get(col) + mystr(" ");
			}

			if(PerformComputation::mbs->UO(UO_LVL_multsim).PrintMsg())
			{
				PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("Simulation result: generation = " + mystr(gen) + mystr(", child=") + mystr(child) +	mystr(", n=") + mystr(np) + "\n  parameter set: " + str_file + " cost function=" + mystr(cost_function_val) + mystr("\n\n"));// possibly problematic during multi-processor computation ==> directly output in file
			}
			PerformComputation::mbs->SolParFile() <<  mystr(mystr(gen) + mystr(" ") + mystr(child) +mystr(" ") + mystr(np) + mystr(" ")+ str_file +mystr(cost_function_val) + mystr("\n")) << flush;// + mystr(" ") + mystr(omp_get_thread_num())		
		}// "candidates" for next generation (multi-threaded processing)
		if (PerformComputation::mbs->StopCalculation()) break;
	}//"childs" <==> furtile population
	return 0;
}

int Optimization::PerformGeneticOptimization()
{
  if(!PerformComputation::mbs->GetModelDataContainer())
	{
		PerformComputation::mbs->UO(UO_LVL_multsim) << "No model parameters found.\nOptimization stopped.";
		PerformComputation::mbs->CloseFiles(); //because it is opened in following function

		return 0;
	}

  TArray<Vector*> allparameters; // for deleting all with new generated Vector's after computation
	int gen=1;
	Vector nominal_parameters = GetModelParameters();
	int run_with_nominal_parameters = PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.run_with_nominal_parameters");
	for (gen=1; gen <= number_of_generations; gen++)
	{
		//todo: if(!run_with_nominal_parameters)PerformComputation::mbs->UO().pUI->StatusText(mystr("[")+mystr(gen-1)+"|"+mystr(number_of_generations)+(mystr("]"))); // show finished generations in status text ... gen-1
		
		if(isRestart && !run_with_nominal_parameters)
		{
			// restart parameter optimization from file
			actparameters.SetLen(0);
			actcost_function_values.SetLen(0);
			int tmp;
			DoRestart(allparameters, tmp); // read generation number and initialize parameters from already computed file
			isRestart = 0; // restart done!
		}
		else
		{
			if(ComputeOneGeneration(gen, nominal_parameters, allparameters, run_with_nominal_parameters))
			{
				return 0; // error during computation
			}	
		}
		ChooseSurvivingPopulation(gen);
		if (PerformComputation::mbs->StopCalculation()) break;
	}

	if(parent_parameters.Length()>0)
	{
		mystr str_file("");
		Vector* paramset = parent_parameters(1);
		for(int col=1;col<=paramset->Length();col++)
		{
			str_file = str_file + paramset->Get(col) + mystr(" ");
		}
		PerformComputation::mbs->SolParFile() <<  mystr(mystr(gen+1) + mystr(" ") + mystr(1) +mystr(" ") + mystr(1) + mystr(" ")+ str_file +mystr(parent_cost_function_values(1)) + mystr("\n"));// << flush;// + mystr(" ") + mystr(omp_get_thread_num())
	
		// final simulation with optimal parameter set
		if( !run_with_nominal_parameters)
		{
			// write optimal parameter in model edc
			SetParameters(paramset);
			PerformComputation::mbs->CloseFiles(); //because it is opened in following function
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n******************************************\n");
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("Final simulation with optimized parameters \n\n");
			PerformComputation::mbs->UO(UO_LVL_0) <<  mystr( str_file +mystr(parent_cost_function_values(1)) + mystr("\n"));
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n******************************************\n");
			if(!PerformComputation::mbs->RepeatedlyPerformComputation(1)) //do not write into solpar file
			{
				PerformComputation::mbs->UO(UO_LVL_err).InstantMessageText("Error: Optimization stopped. Check consistency of system!\n");
				return 1;
			}
		}
	}
	for(int i=1; i <= allparameters.Length(); i++)
	{
		delete [] allparameters(i)->GetVecPtr();
		allparameters(i)->Init();
	}
	allparameters.SetLen(0);

	PerformComputation::mbs->Get_pCFB()->FinishedComputation();

	return 0;
}

// main part of newton optimization
int Optimization::PerformNewtonOptimization()
{
	SetNewtonOptimization(PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.Optimization.Newton.param_epsilon_abs") , PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetDouble("SolverOptions.Optimization.Newton.param_epsilon_rel"));
	
  if(!PerformComputation::mbs->GetModelDataContainer())
	{
		PerformComputation::mbs->UO(UO_LVL_multsim) << "No model parameters found.\nOptimization stopped.";
		PerformComputation::mbs->CloseFiles(); //because it is opened in following function

		return 0;
	}

  TArray<Vector*> allparameters; // for deleting all with new generated Vector's after computation
	int gen=1;  // compute one generation with randomly distributed parameters (shooting method)
	Vector nominal_parameters = GetModelParameters();
	int run_with_nominal_parameters = PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.run_with_nominal_parameters");
	if(!run_with_nominal_parameters)
	{
		 //Optimization.Newton.random_starting_values ... 1: Use random initialization values from Genetic optimization and select best values; 
		 //                                               0: Use model parameter values as starting values.
		run_with_nominal_parameters = !PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Newton.random_starting_values");
		if(!PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Newton.random_starting_values"))
		{
			PerformComputation::mbs->UO(UO_LVL_multsim) << "SolverOptions.Optimization.Newton.random_starting_values not activated --> nominal parameters are used as start values of Newton algorithm \n";
		}
	}
		
	int test = 0;
	if(test)
	{	
		// perform newton computations
		int max_newton_iterations = 30;
		for(int i =1; i<=1; i++)
		{
			double functionValue;
			//PerformNewtonIteration(Vector(0.4,0.2),1.0, functionValue);
			PerformNewtonIterations(Vector(0.8,0.6), functionValue, PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Newton.max_number_of_iterations"));
		}
		PerformComputation::mbs->Get_pCFB()->FinishedComputation();
		return 0;
	}	
	
	if(isRestart && !run_with_nominal_parameters)
	{
		// restart parameter optimization from file
		actparameters.SetLen(0);
		actcost_function_values.SetLen(0);
		DoRestart(allparameters, gen); // read generation number and initialize parameters from already computed file
		isRestart = 0; // restart done!
	}
	else
	{
		if(ComputeOneGeneration(gen, nominal_parameters, allparameters, run_with_nominal_parameters))
		{
			return 0; // error during computation
		}
	}
	ChooseSurvivingPopulation(gen); //chose surviving population as initial values for Newtons method
	
	
	double dummy_val = 1e100;
	Vector act_parameter(parent_parameters(1)->Length());
	Vector best_parameter(parent_parameters(1)->Length());
	best_parameter(1) = dummy_val;

	double best_cost_function = parent_cost_function_values(1); // best cost function value from initial population
	
	if(!run_with_nominal_parameters)
	{
		// perform newton iterations	
		for(int i =1; i<=parent_parameters.Length(); i++)
		{		
			double newton_cost_function;
			Vector act_parameter = PerformNewtonIterations(*parent_parameters(i), newton_cost_function, PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Newton.max_number_of_iterations"));
			if(newton_cost_function < best_cost_function)
			{
				best_cost_function = newton_cost_function;
				best_parameter = act_parameter;	
			}
		}
	}

	if(best_parameter(1) == dummy_val)
	{
		// no better parameter found by newton algorithm
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n****************************************************************************\n");
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\nNewton optimization done after ") + mystr(GetNewtonFunctionCalls()) + mystr(" simulations.\n");
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\nBest parameters from initial parameters: ") << *parent_parameters(1) << ", best cost function value = " << best_cost_function << "\n";
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n****************************************************************************\n");
	}
	else
	{
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n****************************************************************************\n");
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\nNewton optimization finished after ") + mystr(GetNewtonFunctionCalls()) + mystr(" simulations.\n");
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\nBest parameters: ") << best_parameter << ", best cost function value = " << best_cost_function << "\n";
		PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n****************************************************************************\n");
	}

	if(0 && parent_parameters.Length()>0)
	{
		mystr str_file("");
		Vector* paramset = parent_parameters(1);
		for(int col=1;col<=paramset->Length();col++)
		{
			str_file = str_file + paramset->Get(col) + mystr(" ");
		}
		PerformComputation::mbs->SolParFile() <<  mystr(mystr(gen+1) + mystr(" ") + mystr(1) +mystr(" ") + mystr(1) + mystr(" ")+ str_file +mystr(parent_cost_function_values(1)) + mystr("\n"));// << flush;// + mystr(" ") + mystr(omp_get_thread_num())
	
	
		if( !PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.run_with_nominal_parameters"))
		{
			// write optimal parameter in model edc
			SetParameters(paramset);
			PerformComputation::mbs->CloseFiles(); //because it is opened in following function
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n******************************************\n");
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("Final simulation with optimized parameters \n\n");
			PerformComputation::mbs->UO(UO_LVL_0) <<  mystr( str_file +mystr(parent_cost_function_values(1)) + mystr("\n"));
			PerformComputation::mbs->UO(UO_LVL_0) << mystr("\n******************************************\n");
			if(!PerformComputation::mbs->RepeatedlyPerformComputation(1)) //do not write into solpar file
			{
				PerformComputation::mbs->UO(UO_LVL_err) << "Optimization stopped!\n";
				return 1;
			}
		}
	}

	for(int i=1; i <= allparameters.Length(); i++)
	{
		delete [] allparameters(i)->GetVecPtr();
		allparameters(i)->Init();
	}
	allparameters.SetLen(0);

	PerformComputation::mbs->Get_pCFB()->FinishedComputation();

	return 0;
}

// additional penalty value for NewtonFunction is computed, in case of components of x are not in the correct limits
double Optimization::PenaltyFunctionValue(const Vector& x) const
{
	double penalty_val = 0.;    // value of penalty function

	double penalty_weight = 1.; // factor for weightning of penalty function

	int flag = PerformComputation::mbs->GetMBS_EDC_Options()->TreeGetInt("SolverOptions.Optimization.Newton.use_param_limits",0);

	if(flag)
	{
		// 0...no limit of optimized parameter values, 1...use param_[min|max]val, 2...assume all parameters positive, -1...assume all parameters negative.
		for(int i=1;i<=x.Length();i++)
		{
			double infinity = 1e300;
			double lowerlimit = -infinity; // lower parameter limit
			double upperlimit =  infinity; // upper parameter limit
			if(flag == 1)
			{
				lowerlimit = param_minvals(i);
				upperlimit = param_maxvals(i);
			}
			else if(flag == 2)
			{
				// all parameters positive
				lowerlimit = 0.;
			}
			else if(flag == -1)
			{
				// all parameters negative
				upperlimit = 0.;
			}
		
			double delta = 0.; //distance, how far the i-th compontent of parameter vector is over the limit)
			if(x.Get(i)>upperlimit)
			{
				delta = x.Get(i) - upperlimit;
			}
			else if(x.Get(i)<lowerlimit)
			{
				delta = lowerlimit - x.Get(i);
			}
			
			if(delta>0)
			{
				//required characteristics of penalty function g: g(delta=0) = 0, g'(delta=0)=0, ()'=d/d_delta
				penalty_val += exp(delta)-delta-1;
				//penalty_val += 0.5*Sqr(delta);   //quadratic function also possible, but exponential function chosen (faster increasing)
			}
		}
	}
	return penalty_weight*penalty_val;
}

double Optimization::NewtonFunction(const Vector& x) // from class NewtonOptimization
{
	IncreaseNewtonFunctionCalls();
	SetParameters(&x);

	PerformComputation::mbs->CloseFiles(); //because it is opened in following function
	//repeatedly call computation:
	PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("\n****************************");
	PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("\nEvaluate Newton Function: \n\n");
	PerformComputation::mbs->RepeatedlyPerformComputation(1); //do not write into solpar file
	PerformComputation::mbs->OpenFiles(1); //append		
	//*****
	double cost_function_val = PerformComputation::mbs->EvalSensorCostFunctionVal();

	mystr str_file;
	for(int col=1;col<=x.Length();col++) { str_file = str_file + x.Get(col) + mystr(" ");	 }
	if(PerformComputation::mbs->UO(UO_LVL_multsim).PrintMsg()) { PerformComputation::mbs->UO(UO_LVL_multsim) << mystr("Simulation result: ")+ mystr(" generation = ") + mystr(2) + mystr(", child=") + mystr(0) +	mystr(", n=") + mystr(0) +  mystr("\n  parameter set: ") + str_file + " cost function=" << mystr(cost_function_val) << mystr("\n\n");}// possibly problematic during multi-processor computation ==> directly output in file
	PerformComputation::mbs->SolParFile() <<  mystr(mystr(2) + mystr(" ") + mystr(0) +mystr(" ") + mystr(0) + mystr(" ")+ str_file +mystr(cost_function_val) + mystr("\n")) << flush;
	
	return cost_function_val;
}