//#**************************************************************
//#
//# filename:             femath.cpp
//#
//# author:               Gerstmayr Johannes
//#
//# generated:            15.05.97
//# description:          Classes for linear and nonlinear algebra which is
//#												thought to be used in finite element or similar applications
//#												There are 2D and 3D and arbitrary size Vectors (Vector2D, Vector3D, Vector),
//#												arbitrary size matrices (Matrix), and a nonlinear system (NumNLSys)
//#												and a nonlinear solver (NumNLSolver)
//# remarks:							Indizes run from 1 to n except in Vector3D/2D
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
 
#include "ioincludes.h"

#include <assert.h>
#include <memory.h>

#include <math.h>

#include "tarray.h"    
#include "mystring.h"  
#include "femath.h"  
#include "femathHelperFunctions.h"




//for output possibility within femath
#include "..\workingmodule\WorkingModuleBaseClass.h"
extern UserOutputInterface * global_uo;
//#include "..\SuperLU\slu_ddefs.h"

#ifdef gencnt 
extern int * genvec;
extern int * genmat; 
#endif

//extern double global_time;

//ofstream testout("..\\..\\output\\test.txt");

TArray<double> TMtspent(TMlen);
TArray<double> TMtstart(TMlen);
TArray<mystr> TMmsg(TMlen);

void TMResetTimer()
{
	TMmsg(1) = mystr("total time");
	TMmsg(2) = mystr("DrawSystem()");
	TMmsg(3) = mystr("Full Jacobian");
	TMmsg(4) = mystr("Partial Jacobian Additioning");
	TMmsg(5) = mystr("mbs->EvalF2");
	TMmsg(6) = mystr("mbs->EvalM");
	TMmsg(7) = mystr("mbs->EvalF");
	TMmsg(8) = mystr("mbs->EvalG"); 
	TMmsg(9) = mystr("Dependencies");
	TMmsg(10) = mystr("Newton Solve");
	TMmsg(11) = mystr("Solve");
	TMmsg(12) = mystr("Mult(Matrix,Vector,Vector&) general");
	TMmsg(13) = mystr("NLF");
	TMmsg(14) = mystr("Factorize (sparse+dense jacobi)");
	TMmsg(15) = mystr("Apply (sparse+dense jacobi)");
	TMmsg(16) = mystr("dCdq+AddElementLoad (Element::EvalF2)");
	TMmsg(17) = mystr("Quadratic VV");
	TMmsg(18) = mystr("Rigid3D-RotMat");
	TMmsg(19) = mystr("Partial Jacobian");
	TMmsg(20) = mystr("PostNewtonStep"); //$ LA 2011-2-1: Renameing of old NonlinStep
	TMmsg(21) = mystr("TimeStep");
	TMmsg(22) = mystr("EvalF2 CMS");
	TMmsg(23) = mystr("EvalM CMS");
	TMmsg(24) = mystr("Build contact search trees");
	TMmsg(25) = mystr("Contact Search");
	TMmsg(26) = mystr("Compute gap");
	TMmsg(27) = mystr("Nearest Segment");
	TMmsg(28) = mystr("test3");
	TMmsg(29) = mystr("test4");

	for(int i=1; i<= TMlen; i++)
	{
		TMtstart(i) = 0;
		TMtspent(i) = 0;
		//(*global_uo) << i << " " << TMmsg(i) << "\n";
	}
}

//ofstream timings("..\\..\\output\\timings.txt"); //Ö
//ofstream timings2("..\\..\\output\\timings2.txt"); //Ö


void TMPrintTimer()
{
#ifdef timeron
	global_uo->SaveLocalMessageLevel();
	global_uo->SetLocalMessageLevel(UO_LVL_all);

	double div = TMGetTimer(1);
	if (TMGetTimer(1) == 0) {div = 1;}
	for(int i=1; i<= TMlen; i++)
	{
		double tt = TMGetTimer(i);
		if (i>1) {tt=tt/div;}
		char tstr[30];
		if (i!=1)  
		{
			if (tt*100 < 10)
				sprintf_s(tstr,"  %06.4f",tt*100);
			else
				sprintf_s(tstr,"%06.4f",tt*100); 
			(*global_uo) << mystr("[")+mystr(i)+mystr("]")+mystr(" :\t")+mystr(tstr)+mystr("% :  ")+mystr(TMmsg(i))+mystr("\n");
		}
		else
		{
			sprintf_s(tstr,"%4.2f s",tt); 
			(*global_uo) << mystr("[")+mystr(i)+mystr("]")+mystr(" :\t")+mystr(tstr)+mystr(" :  ")+mystr(TMmsg(i))+mystr("\n");
		}
	}
	double acc = TMGetTimer(1);
	if (acc != 0)
	{
		(*global_uo) << "measurement acc.=" << 2./acc << "%\n";
	}
	(*global_uo) << "TimeStep should be " << (TMGetTimer(3)+TMGetTimer(4)+TMGetTimer(5)+
		TMGetTimer(6)+TMGetTimer(7)+TMGetTimer(8)+TMGetTimer(9)+TMGetTimer(11)+TMGetTimer(12)+TMGetTimer(20))/div*100 << "%.\n";

	double tM=TMGetTimer(6);
	double tK=TMGetTimer(5)-TMGetTimer(16);
	double tC=TMGetTimer(8)+TMGetTimer(16);
	double tS=TMGetTimer(14)+TMGetTimer(15);
	double sum = tM+tK+tC+tS;
	double sum2 = TMGetTimer(1);
	double tR = sum2-sum;

	if(sum2 != 0)
	{		
		(*global_uo) << "CPU_time=" << TMGetTimer(1) << "\n";
		(*global_uo) << "CPU%MKCSR=" << " " << tM/sum2*100 << "% " << tK/sum2*100 << "% " << tC/sum2*100 
		<< "% " << tS/sum2*100 << "% " << tR/sum2*100 << "%\n"; 
	}
/*
	timings << TMGetTimer(1) << " " << tM/sum*100 << "% " << tK/sum*100 << "% " << tC/sum*100 << "% " << tS/sum*100 << "%\n" << flush; 
	timings2 << TMGetTimer(1) << " " << tM/sum2*100 << "% " << tK/sum2*100 << "% " << tC/sum2*100 
		<< "% " << tS/sum2*100 << "% " << tR/sum2*100 << "%\n" << flush; 
*/
	global_uo->ResetLocalMessageLevel();
#endif
}
//$ YV 2013-01-02: this is the implementation of timer functions for the kernel, the arrays can be directly accessed
#ifdef timeron
inline void TMStartTimer(const int& i)
{
	TMtstart(i) = GetClockTime();
}

inline void TMStopTimer(const int& i)
{
	TMtspent(i) += GetClockTime() - TMtstart(i);
}
#endif
double TMGetTimer(int i)
{
	if (TMtspent(i) < 1e-6) return 0;
	else return TMtspent(i);
}

/*double Maximum(double a, double b, double c)
{
if (a>=b && a>=c) return a;
else if (b>=c) return b;
return c;
}
double Minimum(double a, double b, double c)
{
if (a<=b && a<=c) return a;
else if (b<=c) return b;
return c;
}*/

Matrix GetDiag(int n, double val)
{
	int i;
	Matrix m(n,n);
	for (i=1; i<=n; i++)
		m(i,i)=val;

	return m;
}

int IsEqual(const IVector& v1, const IVector& v2)
{
	if (v1.Length() != v2.Length()) return 0;

	for (int i=1; i <= v1.Length(); i++)
		if (v1(i) != v2(i)) return 0;

	return 1;
}

int Find(int i, const IVector& v)
{
	for (int j=1; j <= v.Length(); j++)
	{
		if (v(j) == i) return j;
	}
	return 0;
}

//interpolate several values by means of cosine functions
double InterpolateSmooth(double t, double t0, double t1, double A)
{
	double period = (t1-t0); 
	double tt = (t-t0)/period;
	//double omega = MY_PI/period; //only half period of cosine fcn
	//return A*(1.0-cos(omega*(t-t0)))/2.0; 
	return A*(10.*Cub(tt)-15.*To4(tt)+6*To5(tt)); 
}

//smooth interpolation over time t of function values f0,f1,f2,f3 at times t0,t1,t2,t3
double InterpolateFnt2(double t, double t0, double t1, double t2, double t3, double f0, double f1, double f2, double f3)
{
	if(t < t0)	{		return f0;	}
	else if(t < t1)
	{
		return InterpolateSmooth(t, t0, t1, f1-f0) + f0; 
	}
	else if(t < t2)
	{
		return InterpolateSmooth(t, t1, t2, f2-f1) + f1; 
	}
	else if(t < t3)
	{
		return InterpolateSmooth(t, t2, t3, f3-f2) + f2;
	}
	else return f3;
}

//smooth interpolation over time t of function values f0,f1,f2,f3 at times t0,t1,t2,t3
double InterpolateFntArray(double t, TArray<double> t_vec, TArray<double> f_vec)
{
	if (t_vec.Length() == 0 || t_vec.Length() != f_vec.Length()) return 0;

	if(t < t_vec(1))	{return f_vec(1);}
	for (int i=2; i <= t_vec.Length(); i++)
	{
		if (t < t_vec(i))
		{
			return InterpolateSmooth(t, t_vec(i-1), t_vec(i), f_vec(i)-f_vec(i-1)) + f_vec(i-1); 
		}
	}
	return f_vec.Last();
}

double Mises(const Matrix& m)
{
	if (m.Getcols()!=3) return 0;

	Matrix s = m-1./3.*m.Trace()*GetDiag(3);
	return sqrt(3./2.*(s*s).Trace());

}

#ifdef gencnt
int GenVecCnt() 
{
	return *genvec;
}
int GenMatCnt() 
{
	return *genmat;
}
#endif


void SparseJacMat::ApplyTransform(Vector& v, int mode, NumNLSys* nls)
{
	if (nls->TransformJacApply())
	{
		vtrans = v;
		nls->ApplyTransform(vtrans,v,mode);
	}
}

//********************** SPARSEJACMAT **********************************
void SparseJacMat::Apply(const Vector& R, Vector& d, NumNLSys* nls)
{

	TMStartTimer(15);

	if (nls->UseSuperLU())
	{

		d = R;
		//if (solvertype == PARDISOINVERSE)
		//{
		//	pardisoinv->Apply(d);
		//}
		//else
		//{
		//	superluinv->SuperLUApply(d);
		//}
		sparseinv->Apply(d);
	}
	else
	{
		//JG2013-10: not available open source
		assert(0 && "Error: Band LU not available!");

		//d.SetLen(J_vv.Getrows()+J_zz.Getrows());
		////mystatic Vector temp;
		////mystatic Vector temp2;
		////mystatic Vector Rsort;

		//Rsort.SetLen(R.Length());
		//if (resortsize)
		//{
		//	//(*global_uo) << "resort R\n";
		//	for (int i = 1; i <= R.Length(); i++)
		//		Rsort(i) = R(resortvector(i));
		//} 
		//else
		//{
		//	Rsort = R;
		//}


		//d.SetAll(0);

		//Vector R_v, R_z;
		//R_v.LinkWith(Rsort,1,J_vv.Getcols());
		//R_z.LinkWith(Rsort,1+J_vv.Getcols(),J_zz.Getcols());

		//Vector q_v; q_v.LinkWith(d,1,J_vv.Getcols());
		//Vector q_z; q_z.LinkWith(d,1+J_vv.Getcols(),J_zz.Getcols());

		////q_z = J_sch^(-1) * (R_z - (J_zv * (J_vv^(-1) * R_v)));

		//temp = R_v;
		//ApplyTransform(temp,0,nls);

		//if (useband)
		//{
		//	Band_LU_Bks0(J_vv_band,J_vv_band.Getrows(),lu,lu,J_vv_band_aux, J_vv_band_pivot, temp);
		//}
		//else
		//{
		//	J_vv_band.LUBCKSUB(LUindx, temp);
		//}

		//ApplyTransform(temp,1,nls);


		////temp2 = R_v; //compare old and new Bks
		////Band_LU_Bks(J_vv_band,J_vv_band.Getrows(),lu,lu,J_vv_band_aux, J_vv_band_pivot, temp2);
		////(*global_uo) << "diff=" << (temp-temp2).MaxNorm() << "\n";

		//Mult(J_zv,temp,temp2);
		//temp = R_z;
		//temp -= temp2;
		//Mult(J_sch,temp,q_z);

		////q_v = J_vv^(-1) * (R_v - J_vz * q_z);

		//Mult(J_vz,q_z,temp);
		//temp2 = R_v;
		//temp2 -= temp;

		//q_v = temp2;
		//ApplyTransform(q_v,0,nls);

		//if (useband)
		//{
		//	Band_LU_Bks0(J_vv_band,J_vv_band.Getrows(),lu,lu,J_vv_band_aux, J_vv_band_pivot, q_v);
		//}
		//else
		//{
		//	J_vv_band.LUBCKSUB(LUindx, q_v);
		//}

		//ApplyTransform(q_v,1,nls);

		//if (resortsize)
		//{
		//	temp.SetLen(d.Length());
		//	temp = d;
		//	for (int i = 1; i <= d.Length(); i++)
		//		d(resortvector(i)) = temp(i);
		//}
	}
	TMStopTimer(15);
}



int SparseJacMat::Factorize(NumNLSys* nls)
{
	//Daniel: comment out this function!!!
	assert(0 && "Error: SparseJacMat::Factorize(NumNLSys* nls) not available any more!");

	//double epszero = 1e-16; //optimize LUBCKSUB

	////(*global_uo) << "minv=" << minv_test << "\n";
	////mystatic Matrix temp1, temp2;
	////mystatic Vector tempcol;

	//int n = J_vv.Getcols();
	//int bw;
	//int n_zz = J_zz.Getcols();

	//if (LUcomputed) bw = LUcomputed_bw;
	//else bw = J_vv.GetBandWidth();

	////(*global_uo) << "Jvv=" << J_vv << "\n";

	//if (bw < (8*n)/10)
	//{
	//	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//	useband = 1;

	//	int col_band = bw*2-1;

	//	lu = bw-1;

	//	if (!(LUcomputed && nls->TransformJacApply()))
	//	{
	//		LUcomputed = 1;
	//		LUcomputed_bw = bw;
	//		TMStartTimer(27);

	//		J_vv_band_aux.SetSize(n, lu);

	//		J_vv.GetBandMatrix(J_vv_band, lu);
	//		//(*global_uo) << "Jzzband=" << J_vv_band << "\n";

	//		if (!n_vv_written) (*global_uo) << "LU-lalloc=" << (J_vv_band_aux.Getrows()*J_vv_band_aux.Getcols()+J_vv_band.Getrows()*J_vv_band.Getcols())*8./(1024.*1024.) << "MB\n";

	//		J_vv_band_pivot.SetLen(n);
	//		double signd;
	//		int rvb;
	//		rvb = Band_LU_Dec(J_vv_band,n,lu,lu,J_vv_band_aux, J_vv_band_pivot, signd);
	//		if (rvb == 0) {(*global_uo) << "banddec returned " << rvb << "\n"; return 0;}
	//	TMStopTimer(27);
	//	}

	//	//TMStartTimer(28);
	//	//J_sch = J_zz - J_zv * (J_vv^-1 * J_vz):
	//	mtemp1.SetSize(J_vv.Getrows(),J_vz.Getcols());
	//	tempcol.SetLen(J_vz.Getrows());
	//	for (int i = 1; i <= J_vz.Getcols(); i++)
	//	{
	//		J_vz.GetColVec(i,tempcol);
	//		if (!tempcol.IsZero(epszero))
	//		{
	//			ApplyTransform(tempcol,0,nls);
	//			Band_LU_Bks(J_vv_band,n,lu,lu,J_vv_band_aux, J_vv_band_pivot, tempcol);
	//			ApplyTransform(tempcol,1,nls);
	//		}
	//		mtemp1.SetColVec(tempcol,i);
	//	}
	//	//TMStopTimer(28);
	//}
	//else
	//{
	//	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//	//full version, J_vv_band is a full matrix and no banded matrix!!!
	//	useband = 0;

	//	TMStartTimer(27); //11.3% --> 14% with 500x500

	//	if (!(LUcomputed && nls->TransformJacApply()))
	//	{
	//		//(*global_uo) << "Jvv=" << J_vv << "\n";
	//		J_vv_band.CopyFrom(J_vv,1,1,n,n);

	//		//if (!n_vv_written) (*global_uo) << "J_vv=" << J_vv << "\n";
	//		//if (!n_vv_written) (*global_uo) << "J_vv_band=" << J_vv_band << "\n";

	//		LUcomputed = 1;
	//		LUcomputed_bw = bw;

	//		//double signd;
	//		int rvb = 0;
	//		rvb = J_vv_band.LU(LUindx, LUvv);
	//		if (rvb == 0) {(*global_uo) << "LU-decomposition returned " << rvb << "\n"; return 0;}
	//	}
	//	TMStopTimer(27);

	//	//J_sch = J_zz - J_zv * (J_vv^-1 * J_vz)
	//	//part: mtemp = (J_vv^-1 * J_vz):


	//	//TMStartTimer(28); //6% -> 3% with optimized LUBCKSUB
	//	mtemp1.SetSize(J_vv.Getrows(),J_vz.Getcols());
	//	tempcol.SetLen(J_vz.Getrows());
	//	//int cnt = 0;
	//	for (int i = 1; i <= J_vz.Getcols(); i++)
	//	{
	//		J_vz.GetColVec(i,tempcol);

	//		if (!tempcol.IsZero(epszero))
	//		{
	//			//cnt++;
	//			ApplyTransform(tempcol,0,nls);
	//			J_vv_band.LUBCKSUB(LUindx, tempcol);
	//			ApplyTransform(tempcol,1,nls);
	//		}
	//		mtemp1.SetColVec(tempcol,i);
	//	}
	//	//TMStopTimer(28);

	//	//(*global_uo) << "nonzero-cols=" << (double)cnt/((double)J_vz.Getcols()+1e-9)*100. << "% \n";

	//}

	//if (n_vv_written < 1) (*global_uo) << "n_vv=" << n << ", n_zz=" << J_zz.Getcols() << ", bw=" << bw << ", useband=" << useband << "\n";

	////J_sch = J_zz - J_zv * (J_vv^-1 * J_vz);
	////TMStartTimer(25); //dense: 12%, sparse version needs 0.5%
	//mtemp2.SetSize(J_zz.Getrows(),J_zz.Getcols()); //needed???

	//J_zvS.CopyFrom(J_zv);       //sparse version
	//Mult(J_zvS,mtemp1,mtemp2);

	////Mult(J_zv,mtemp1,mtemp2); //dense version

	//J_sch = J_zz;
	//J_sch -= mtemp2;
	////TMStopTimer(25);


	//TMStartTimer(26); //<1%
	//if (n_zz != 0)
	//{
	//	if (!J_sch.Invert2()) {(*global_uo) << "Schur-complement could not be inverted!\n"; return 0;}
	//}
	//TMStopTimer(26);
	//n_vv_written ++;

	return 1;
}




//********************** NUMNLSYS**********************************
UserOutput& NumNLSys::UOfull(int message_level, int output_prec) const
{
	uo.SetLocalMessageLevel(message_level); 
	uo.SetOutputPrec(output_prec); 
	return uo;
};

void NumNLSys::SetSolver(NumNLSolver* s) 
{
	solver = s;
}

int NumNLSys::NLS_GetJacCol() const 
{
	if (solver->SymmetricJacobian()) 
	{
		return jaccol; 
	}
	else {return 0;}
}

int NumNLSys::UseSparseSolver() const 
{
	//if (TransformJacApply()) return 1;
	//return 0;
	return solver->UseSparseSolver();
}

int NumNLSys::SolveUndeterminedSystem() const 
{
	return solver->SolveUndeterminedSystem();
}

int& NumNLSys::SolveUndeterminedSystem() 
{
	return solver->SolveUndeterminedSystem();
}

double NumNLSys::EstimatedConditionNumber() const 
{
	return solver->EstimatedConditionNumber();
}

double& NumNLSys::EstimatedConditionNumber() 
{
	return solver->EstimatedConditionNumber();
}

int NumNLSys::UseSparseJac() const 
{
	//return 0;
	return solver->UseSparseSolver();
}

int NumNLSys::UseSuperLU() const 
{
	return 1;
}

void NumNLSys::Jacobian(Matrix& m, Vector& x)
{
	TMStartTimer(3);
	int i,j;
	//mystatic Vector f1;
	//mystatic Vector f2;
	f1.SetLen(x.GetLen());
	f2.SetLen(x.GetLen());

	m.SetSize(x.GetLen(),x.GetLen());

	double numdiffepsi = solver->NumDiffepsi();
	NLS_SetJacCol(0);
	NLS_SetJacFullNewton(1);
	fullnewtoncnt++;

	if (!solver->SymmetricJacobian())
	{
		NLF(x,f2);
	}
	double eps;
	double xstore;
	(*global_uo) << "len=" << m.Getcols() << "\n";
	for (i = 1; i <= m.Getcols(); i++)
	{
		eps = numdiffepsi*Maximum(1e-2,fabs(x(i)));
		if (solver->SymmetricJacobian())
		{
			//if (TIstages==1) 
			//{
			NLS_SetJacCol(i);
			//}

			xstore = x(i);
			x(i) += eps;
			NLF(x,f1);
			//(*global_uo) << "f1=" << f1 << "\n";

			x(i) -= 2*eps;
			NLF(x,f2);

			//(*global_uo) << "f2=" << f2 << "\n";
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
	TMStopTimer(3);
}

int lastmatrix = 0;
Matrix lastm;

int NumNLSolver::Factorize(Matrix& minv, SparseMatrix& sminv, SparseJacMat& sparseminv)
{
	//Matrix mm = sminv.GetMatrix();
	//PrintMatrix01(mm);
	//global_uo.InstantMessageText("NumSolver::Factorize begin");
	
	global_uo->SaveLocalMessageLevel();
	global_uo->SetLocalMessageLevel(UO_LVL_dbg1);

	TMStartTimer(14);
	int rv = 0;
	if (nls->UseSparseSolver())
	{
		if (nls->UseSparseJac()) //use sminv
		{
			if (nls->UseSuperLU())
			{
				double facttime;
				delete sparseminv.sparseinv;
				sparseminv.sparseinv = new SparseInverse(sminv);
				sparseminv.sparseinv->SetSparseMatrix(sminv);
				rv = 1;
				if (sminv.Getcols()>10000 && global_uo->PrintMsg()) (*global_uo) << "Factorize Matrix with sparse solver .. ";
				facttime = -GetClockTime();
				rv = sparseminv.sparseinv->Factorize();
				facttime += GetClockTime();
				if (sminv.Getcols()>10000 && global_uo->PrintMsg()) (*global_uo) << "done, factorization time " << facttime << "s\n";
				sparseminv.sparseinv->IsFirstStep()=0;

			}
			else
			{
				//JG2013-10: functions not available open source
				assert(0 && "Error: Band Sparse Factorize not available any more!");
			//	//(*global_uo) << "factorize changed\n";
			//	sparseminv.resortsize = nls->GetResortSize();
			//	int rs = nls->GetResortVector().Length();
			//	sparseminv.resortvector.SetLen(rs);

			//	for (int i=1; i <= rs; i++)
			//	{
			//		sparseminv.resortvector(i) = nls->GetResortVector()(i);
			//	} 
			//	/*
			//	(*global_uo) << "J=" << sminv << "\n";
			//	(*global_uo) << "resortsize=" << sparseminv.resortsize << "\n";
			//	(*global_uo) << "rs=" << rs << "\n";
			//	(*global_uo) << "rsvector=" << sparseminv.resortvector << "\n";
			//	*/
			//	if (nls->GetResortSize())
			//	{
			//		sparseminv.J_vv.CopyFrom(sminv,1,1,bandsize,bandsize,nls->GetResortVector());
			//		sparseminv.J_vz.CopyFrom(sminv,1,1+bandsize,bandsize,sminv.Getcols(),nls->GetResortVector());
			//		sparseminv.J_zv.CopyFrom(sminv,1+bandsize,1,sminv.Getrows(),bandsize,nls->GetResortVector());
			//		sparseminv.J_zz.CopyFrom(sminv,1+bandsize,1+bandsize,
			//			sminv.Getrows(),sminv.Getcols(),nls->GetResortVector());
			//	}
			//	else
			//	{
			//		sparseminv.J_vv.CopyFrom(sminv,1,1,bandsize,bandsize);
			//		//global_uo.InstantMessageText("NumSolver::Factorize begin2");
			//		sparseminv.J_vz.CopyFrom(sminv,1,1+bandsize,bandsize,sminv.Getcols());
			//		sparseminv.J_zv.CopyFrom(sminv,1+bandsize,1,sminv.Getrows(),bandsize);
			//		sparseminv.J_zz.CopyFrom(sminv,1+bandsize,1+bandsize,
			//			sminv.Getrows(),sminv.Getcols());

			//		//(*global_uo) << "sminv=" << sminv << "\n";
			//		//(*global_uo) << "J_vv=" << sparseminv.J_vv << "\n";
			//		//(*global_uo) << "J_vz=" << sparseminv.J_vz << "\n";
			//		//(*global_uo) << "J_zv=" << sparseminv.J_zv << "\n";
			//		//(*global_uo) << "J_zz=" << sparseminv.J_zz << "\n";


			//	}

			//	rv = sparseminv.Factorize(nls);
			}
		}
		else
		{
			//JG2013-10: case cannot happen
				assert(0 && "Error: Band Sparse Factorize not available any more!");
			//// bandsize = minv.Getcols(); //for testing
			//sparseminv.resortsize = nls->GetResortSize();
			//if (nls->GetResortSize())
			//{
			//	//TMStartTimer(24); //1% time lost

			//	int rs = nls->GetResortVector().Length();
			//	sparseminv.resortvector.SetLen(rs);
			//	for (int i=1; i <= rs; i++)
			//	{
			//		sparseminv.resortvector(i) = nls->GetResortVector()(i);
			//	}


			//	//(*global_uo) << "sm-diff=" << (sm.GetMatrix()-minv).MaxNorm() << "\n";

			//	//original:
			//	//SparseMatrix sm(minv);
			//	//sparseminv.J_vv.CopyFrom(minv,1,1,bandsize,bandsize,nls->GetResortVector());
			//	
			//	sparseminv.mtemp1.CopyFrom(minv,1,1,bandsize,bandsize,nls->GetResortVector()); //1% time lost
			//	sparseminv.J_vv.CopyFrom(sparseminv.mtemp1);

			//	//sparseminv.J_vv.CopyFrom(minv,1,1,bandsize,bandsize,nls->GetResortVector());
			//	sparseminv.J_vz.CopyFrom(minv,1,1+bandsize,bandsize,minv.Getcols(),nls->GetResortVector());
			//	sparseminv.J_zv.CopyFrom(minv,1+bandsize,1,minv.Getrows(),bandsize,nls->GetResortVector());
			//	sparseminv.J_zz.CopyFrom(minv,1+bandsize,1+bandsize,minv.Getrows(),minv.Getcols(),nls->GetResortVector());

			//	//TMStopTimer(24);

			//}
			//else
			//{
			//	TMStartTimer(24);

			//	sparseminv.J_vv.CopyFrom(minv,1,1,bandsize,bandsize);
			//	sparseminv.J_vz.CopyFrom(minv,1,1+bandsize,bandsize,minv.Getcols());
			//	sparseminv.J_zv.CopyFrom(minv,1+bandsize,1,minv.Getrows(),bandsize);
			//	sparseminv.J_zz.CopyFrom(minv,1+bandsize,1+bandsize,
			//		minv.Getrows(),minv.Getcols());

			//	TMStopTimer(24);

			//}
			////sparseminv.minv_test = minv;
			////sparseminv.minv_test.Invert2();

			//rv = sparseminv.Factorize(nls);
		}
	}
	else
	{
		if (nls->UseSparseJac()) //use sminv
		{
			minv.CopyFrom(sminv,1,1,sminv.Getrows(),sminv.Getcols());
		}

		rv = minv.Invert2();

	}
	TMStopTimer(14);

	global_uo->ResetLocalMessageLevel();

	return rv;
}


//********************** APPLYJAC **********************************
int NumNLSolver::ApplyJac(const Vector& f, Vector& res)
{
	SaveJac* jac;
	int chj = ChooseJac();
	//(*global_uo) << "nlsinfo=" << nonlinsolveinfo
	//	<< ", nls1=" << sjac[0].nlsinfo << ", nls2=" << sjac[1].nlsinfo << "\n";
	if (chj != 0)
	{
		jac = &(sjac[chj-1]);
		Mult(jac->oldjacmat,f,res);
		return 1;
	}
	else return 0;
}


//********************** NUMNLSOLVER**********************************
int NumNLSolver :: NLSolve(Vector& x0)
{
	if (SolverPrintsDetailedOutput())
	{
		SolverUO() << "-->";
		if (ModifiedNewton())
		{
			SolverUO() << "Modified Newton:";
		}
		else
		{
			SolverUO() << "Full Newton:";
		}
	}

	TMStartTimer(10);
	//(*global_uo) << "NLSolve\n";
	//(*global_uo) << "genvec=" << genvec << ", ";
	//(*global_uo) << "genmat=" << genmat << "\n";
	//(*global_uo) << "nls-inf=" << this->nonlinsolveinfo << "\n";
	//(*global_uo) << "  xstart0=" << x0 << "\n";

	jaccount = 0;
	nls->SetSolver(this);
	//mystatic Vector x0start;
	x0start = x0;  
	nls->NLS_SetJacCol(0);

	//mystatic Vector f;
	f.SetLen(x0.GetLen());
	//mystatic Vector xd;
	xd.SetLen(x0.GetLen());

	const double MAXERR = 1e100;
	const double MAXERRINC = 1e10; //$JG 2012-11-07 changed from 1e10 to 1e25 because of problems with Cable2D element in first step
	int it = 0;

	nls->NLF(x0,f); 


	double errinit = f.MaxNorm();
	double err = errinit;
	if (errinit < absoluteaccuracy)
	{
		errinit = absoluteaccuracy;
	}
	double error_goal = relativeaccuracy*errinit;    // error_goal = relativeaccuracy * min(absoluteaccuracy, errinit)

	if (SolverPrintsDetailedOutput())
	{
		mystr str(" goal=");
		str += mystr(error_goal);
		str += mystr(" (err=");
		str += mystr(err);
		str += mystr(", rel_acc=");
		str += mystr(relativeaccuracy);
		str += mystr(", abs_acc=");
		str += mystr(absoluteaccuracy);
		str += mystr(", [goal=max(err,abs_acc)rel_acc, err=MaxNorm(residual)])\n");

		SolverUO() << str;

		if (nls->GetOptions()->LoggingOptions()->SolverNewtonIterationResidualVector())
		{
			SolverUO() << "residual vector = " << f << "\n";
		}
	}

	double lasterr = err;
	int jacsingular = 0;
	int stopnewton = 0;

	double maxcontractivity = 0;
	contractivity = 0.5;
	Matrix* minv; //full matrix
	SparseMatrix* sminv; //sparse matrix build
	SparseJacMat* sparseminv; //sparse matrix factorize

	//initialize for full newton:
	SaveJac* jac = &sjac[0];
	minv = &(jac->oldjacmat);
	sminv = &(jac->oldsjacmat);
	sparseminv = &(jac->sparsejacmat); 

	if (ModifiedNewton())
	{
		int chj = ChooseJac();
		//chj=0; //every timestep a Jacobian
		if (chj != 0)
		{
			jac = &(sjac[chj-1]);
		}
		else
		{
			if (sjac[0].oldjac == 0) {jac = &sjac[0]; }
			else if (sjac[1].oldjac == 0) {jac = &sjac[1];}
			else if ((fabs(sjac[0].nlsinfo) >= fabs(nonlinsolveinfo) || fabs(sjac[1].nlsinfo) >= fabs(nonlinsolveinfo))
				&& sjac[0].lastbuilt) {jac = &sjac[1];}
			else if ((fabs(sjac[0].nlsinfo) >= fabs(nonlinsolveinfo) || fabs(sjac[1].nlsinfo) >= fabs(nonlinsolveinfo))
				&& sjac[1].lastbuilt) {jac = &sjac[0];}
				//else if ((sjac[1].oldjacage < sjac[0].oldjacage && sjac[1].lastbuilt)) {jac = &sjac[0];}
				//else if ((sjac[0].oldjacage < sjac[1].oldjacage && sjac[0].lastbuilt)) {jac = &sjac[1];}
			else if ((fabs(sjac[1].nlsinfo - nonlinsolveinfo) > fabs(sjac[0].nlsinfo - nonlinsolveinfo))) {jac = &sjac[1]; }
			else {jac = &sjac[0];}

			jac->oldjac = 0;
		}

		minv = &(jac->oldjacmat);
		sminv = &(jac->oldsjacmat);
		sparseminv = &(jac->sparsejacmat);

		if (jac->oldjac == 0 || jac->oldjacage > jac->maxjacage || jac->oldjacsize != x0.GetLen())
		{

			// report on the reason to update jacobian
			if (SolverPrintsDetailedOutput())
			{
				mystr jacobian_update_info_sparsity = (nls->UseSparseJac() ? mystr("Sparse ") : mystr("Full "));
				mystr jacobian_update_info = jacobian_update_info_sparsity + mystr("Jacobian (#") + mystr(jaccount) + mystr(") is computed");
				if (jac->oldjac == 0)
				{
					jacobian_update_info += mystr(".\n");
				}
				else if (jac->oldjacage > jac->maxjacage)
				{
					jacobian_update_info += mystr(", since the age of the Jacobian was greater than ") + mystr(jac->maxjacage) + mystr(".\n");
				}
				else //jac->oldjacsize != x0.GetLen()
				{
					jacobian_update_info += mystr(", since the length of the solution vector did not fit the size of the Jacobian.\n");
				}
				SolverUO() << jacobian_update_info;
			}

			AssembleJacobian(minv, sminv, x0);

			/*
			nls->Jacobian(*sminv,x0); 
			nls->Jacobian(*minv,x0);
			//(*global_uo) << "Jac=" << *minv;
			//(*global_uo) << "Jac_s=" << sminv->GetMatrix();
			//(*global_uo) << "Jac-diff=" << *minv-sminv->GetMatrix();
			(*global_uo) << "N-Jac-diff=" << (*minv-sminv->GetMatrix()).MaxNorm() << "\n";
			*/

			if (!Factorize(*minv,*sminv,*sparseminv))
			{
				jacsingular = 1000;
			};

			jac->oldjac = 1;
			jac->nlsinfo = nonlinsolveinfo;
			jac->oldjacage = 0;
			jac->oldjacsize = x0.GetLen();
			sjac[0].lastbuilt = 0;
			sjac[1].lastbuilt = 0;
			jac->lastbuilt = 1;
		}
		else
		{
			jac->oldjacage++;
		}

	}

	int n_modstopped = 0;
	stopmnr = 0;
	int contit = 0;

	while ((it < MaxNewtonSteps() && err >= error_goal && jacsingular < 2 && !stopnewton) || it==0)
	{
		it++;
		contit++;
		int maxmodnewtonsteps1 = MaxModNewtonSteps();
		int maxmodnewtonsteps2 = MaxModNewtonSteps()+MaxRestartNewtonSteps();
	
		double low_contractivity_tolerance = 0.7;
		double high_contractivity_tolerance = 2.;
		if (((it == maxmodnewtonsteps1) || (contit>2 && contractivity > low_contractivity_tolerance)) && !stopmnr && ModifiedNewton())  //contractivity > 0.7 is optimum for many cases
		{
//	(*global_uo) << "maxmodnewtonsteps1=" << maxmodnewtonsteps1 << ",contr=" << contractivity << "\n";
			double memo_contractivity = contractivity;
			if (n_modstopped==0) 
			{
				x0=x0start;
				contractivity = 0.5;
			}
			if ((contractivity >= high_contractivity_tolerance && n_modstopped>0) || n_modstopped>2)
			{
				if (SolverPrintsDetailedOutput())
				{
					SolverUO() << mystr("Newton error was set to MAXERR=")+mystr(MAXERR)+mystr(", since ");
					if (contractivity >= high_contractivity_tolerance)
					{
						SolverUO() << mystr("contractivity=")+mystr(contractivity)+mystr(" exceeded highly critical value of ")+mystr(high_contractivity_tolerance)+mystr(".\n");
					}
					else
					{
						SolverUO() << mystr("Jacobian has been updated more than twice, but still ");
						if (contractivity > low_contractivity_tolerance)
						{
							SolverUO() << mystr("contractivity=")+mystr(memo_contractivity)+mystr(" exceeded critical value of ")+mystr(low_contractivity_tolerance)+mystr(".\n");
						}
						else
						{
							SolverUO() << mystr("iterations=")+mystr(it)+mystr(" exceeded critical value of ")+mystr(maxmodnewtonsteps1)+mystr(".\n");
						}
					}
				}
				err=MAXERR;
			}
			else
			{
				contit = 0;
				n_modstopped++;


//	(*global_uo) << "jac1\n";

				if (SolverPrintsDetailedOutput())
				{
					mystr jacobian_update_info_sparsity = (nls->UseSparseJac() ? mystr("Sparse ") : mystr("Full "));
					mystr jacobian_update_info = jacobian_update_info_sparsity + mystr("Jacobian (# ") + mystr(jaccount) + mystr(") is computed");
					if (it == maxmodnewtonsteps1)
					{
						jacobian_update_info += mystr(", since iterations=") + mystr(it) + mystr(" exceeded critical value of ") + mystr(maxmodnewtonsteps1) + mystr(".\n");
					}
					else
					{
						jacobian_update_info += mystr(", since contractivity ") + mystr(memo_contractivity) + mystr(" exceeded critical value of ") + mystr(low_contractivity_tolerance) + mystr(", and contit is greater than 2.\n");
					}
					SolverUO() << jacobian_update_info;
				}

				AssembleJacobian(minv, sminv, x0);

				nls->NLF(x0,f); 
				lasterr = f.MaxNorm();

				if (!Factorize(*minv,*sminv,*sparseminv))
				{
					jacsingular = 1000;
				};


				jac->oldjac = 1;
				jac->nlsinfo = nonlinsolveinfo;
				jac->oldjacage = 0;
				jac->oldjacsize = x0.GetLen();
				sjac[0].lastbuilt = 0;
				sjac[1].lastbuilt = 0;
				jac->lastbuilt = 1;
			}
		}
		if (ModifiedNewton() && (it > maxmodnewtonsteps2 || err == MAXERR) && !stopmnr) 
		{
			if (SolverPrintsDetailedOutput())
			{
				SolverUO() << mystr("Switch from Modified to Full Newton, since ");
				if (it > maxmodnewtonsteps2)
				{
					SolverUO() << mystr("iterations=")+mystr(it)+mystr(" exceed highly critical value of ")+mystr(maxmodnewtonsteps2)+mystr(".\n");
				}
				else
				{
					SolverUO() << mystr("Newton error=MAXERR=")+mystr(MAXERR)+mystr(".\n");
				}
			}

			stopmnr = 1; 
			x0 = x0start; 
			nls->NLF(x0,f); 
			jac->oldjac = 0; 
			jac->jacfailedcnt++;
		}
		//if (jac->jacfailedcnt > 20) {jac->maxjacage = jac->maxjacage/2+1; jac->jacfailedcnt = 0;}

		if (!ModifiedNewton() || stopmnr) 
		{
			AssembleJacobian(minv, sminv, x0);

			// factorize Jacobian & solve system
			contractivity = 0;
			fulljaccnt++;
			//mystatic IVector iv;
			int rv;

			if (nls->UseSparseSolver())
			{
				rv = Factorize(*minv,*sminv,*sparseminv);
				if (rv) sparseminv->Apply(f,xd,nls);
			}
			else
			{
				if (nls->SolveUndeterminedSystem())
				{
					int minv_rank = 0;
					double rcond = 1./nls->EstimatedConditionNumber();
					rv = minv->UndeterminedSystemSolve(f, xd, rcond, minv_rank);
					int minv_dim = min(minv->Getrows(), minv->Getcols());
					if (SolverPrintsDetailedOutput() && minv_rank < minv_dim)
					{
						SolverUO() << "WARNING: Estimated rank of Jacobian (" << minv_rank << ") smaller than dimension (" << minv_dim << ")\n";
					}
				}
				else
				{
					rv = minv->Solve(f,xd);
				}
			}
			
			
			if (0) // old branch
			{
				rv = minv->LU(iv);
				if (rv < 0) rv=0;
				xd=f;
				if (rv) minv->LUBCKSUB(iv,xd);
			}

			if (!rv) 
			{
				x0 = x0start; 
				xd.SetAll(0);
				nls->NLF(x0,f); 
				jacsingular+=1;
			};
		}
		else
		{
			//TMStartTimer(24);
			if (nls->UseSparseSolver())
			{
				sparseminv->Apply(f,xd,nls);
			}
			else
			{
				TMStartTimer(15);
				Mult(*minv,f,xd);
				TMStopTimer(15);
				//(*global_uo) << "xdf=" << xd << "\n";
			}

			//TMStopTimer(24);
		}

		double trust = 1;

		if (trustregion)
		{
			Vector x;
			double lasterr, acterr;
			trust = trustregiondiv;
			lasterr = 1e50;

			double ftrust = 0;

			while (trust <= 2)
			{
				x=x0-trust*xd; //Vector generation!!!!
				nls->NLF(x,f);

				acterr = f.MaxNorm();

				if (acterr < lasterr) {ftrust = trust; lasterr = acterr;}

				trust += trustregiondiv;
			}

			trust = ftrust;
			if (ftrust == 0) {trust = 1;}

		}

		//trust = 0.8;//0.8 work best, DASSL uses 0.75 - but this leads to many Jacobians sometimes!!!; 1-jac->oldjacage*0.025 ??;
		if (trust != 1) {xd.Mult(trust);}

		x0 -= xd;
		nls->NLF(x0,f);
		err = f.MaxNorm();

		if (SolverPrintsDetailedOutput()) 
		{
			mystr str("Nit=");
			str += mystr(it);
			str += mystr(": err=");
			str += mystr(err);
			if (ModifiedNewton() && !stopmnr)
			{
				str += mystr(", contractivity=");
				str += mystr(contractivity);
			}
			str += mystr("\n");
			SolverUO() << str;
		
			if (nls->GetOptions()->LoggingOptions()->SolverNewtonIterationSolutionVector())
			{
				SolverUO() << "solution vector = " << xd << "\n";
			}
			if (nls->GetOptions()->LoggingOptions()->SolverNewtonIterationResidualVector())
			{
				SolverUO() << "residual vector = " << f << "\n";
			}
			//EK & PG 2012-05-18 at the end of the computation no recomputation of jacobian necessary - see comment below
			//if (nls->GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiMatrix())
			//{
			//	if (nls->UseSparseJac())
			//	{
			//		//SparseMatrix m; nls->Jacobian(m,f);   //PG 2012-05-18: recalculation of sparse jacobian lead to different results for (SolverNewtonIterationJacobiMatrix()==0/1)
			//		Matrix m = sminv->GetMatrix();
			//		SolverUO() << "sparse jacobi matrix = " << m << "\n";
			//	}
			//	else
			//	{
			//		//Matrix m; nls->Jacobian(m,f);      //PG 2012-05-18: recalculation not necessary, however no difference behaviour (as in sparse case) was encountered
			//		Matrix m(*minv);
			//		SolverUO() << "jacobi matrix = " << m << "\n";
			//	}
			//}
		}

		if (lasterr != 0)
		{
			if (contractivity == 0)
			{
				contractivity = err/lasterr; 
			}
			else
			{
				contractivity = sqrt(contractivity*err/lasterr);
			}
		}
		else
		{
			contractivity = 0.1;
		}

		if (it>4) maxcontractivity=Maximum(maxcontractivity,contractivity);
		lasterr = err;

		//MaSch 2013-02-04: temporarily changed the MAXERRINC-test to avoid comparison to absolute numbers in the first Newton step (where x0start = 0)
		//due to problems with very small and thin (length ~ 1e-6 m) ANCFCable2D elements
		//(test case: finite bending deformation by constant load in x-direction with clamped end, r' in x-dir gets very large in the first step (scaling linearly with force/density),
		//such that x0.MaxNorm() > MAXERRINC*(x0start.MaxNorm()+1.) = MAXERRINC in the first Newton step
		//static computation without problems; 

		//if (!x0.IsValid(MAXERR) || IsNaN(err) || (err >= MAXERR) || (err <= -MAXERR) || x0.MaxNorm() > MAXERRINC*(x0start.MaxNorm()+1.)) //err is Nan or out of range
		if (!x0.IsValid(MAXERR) || IsNaN(err) || (err >= MAXERR) || (err <= -MAXERR) || (x0start.MaxNorm()==0. ? false : (x0.MaxNorm() > MAXERRINC*(x0start.MaxNorm())))  ) //err is Nan or out of range
		{
			if(SolverPrintsDetailedOutput())
			{
				mystr str("err was set to MAXERR=");
				str+=mystr(MAXERR);
				str+=mystr(", since\n");
				if (!x0.IsValid(MAXERR)) str+=mystr("  * max-norm of x0 exceeds critical value of MAXERR=")+mystr(MAXERR)+mystr("\n");
				if (IsNaN(err)) str+=mystr("  * max-norm of residual is nan\n");
				if (err >= MAXERR || err <= -MAXERR) str+=mystr("  * max-norm of residual exceeds critical value of MAXERR=")+mystr(MAXERR)+mystr("\n");
				if (x0.MaxNorm() > MAXERRINC*(x0start.MaxNorm()+1.)) str+=mystr("  * max-norm of x0 exceeds critical value of MAXERRINC*(x0start.MaxNorm()+1.)=")+mystr(MAXERRINC*(x0start.MaxNorm()+1.))+mystr("\n");
				SolverUO() << str;
			}
			err = MAXERR;

			if (!ModifiedNewton() || stopmnr) 
			{
				stopnewton = 1;

				if(SolverPrintsDetailedOutput())
				{
					mystr str("Full Newton method stopped!\n");
					SolverUO() << str;
				}
			}
			else
			{
			}
		}
	}

	newtonits = it;
	contractivity = maxcontractivity;
	TMStopTimer(10);

	if (!x0.IsValid(MAXERR)) //values of solution are out of range
	{
		error_msg += "Solution diverged! ";
		if (SolverPrintsDetailedOutput())
		{
			SolverUO() << "  error message=Solution diverged! " << "\n";
		}
		if (jacsingular != 0)
		{
			error_msg += "Jacobian is singular! ";
			if (SolverPrintsDetailedOutput())
			{
				SolverUO() << "  error message=Jacobian is singular! " << "\n";
			}
		}
		return 0;
	}

	if (err == MAXERR)  //err is NaN
	{
		error_msg += "Residual went to infinity! ";
		if (SolverPrintsDetailedOutput()) {SolverUO() << "  error message=Residual went to infinity! " << "\n";}
		if (jacsingular != 0)
		{
			error_msg += "Jacobian is singular! ";
			if (SolverPrintsDetailedOutput()) {SolverUO() << "  error message=Jacobian is singular! " << "\n";}
		}
		x0 = x0start;
		return 0;
	}

	if ((err < error_goal) && !jacsingular && !stopnewton) 
	{
		return 1;
	} 
	else 
	{
		if (jacsingular != 0)
		{
			error_msg += "Jacobian is singular! ";
			if (SolverPrintsDetailedOutput()) {SolverUO() << "  error message=Jacobian is singular! " << "\n";}
		}
		else if (err >= error_goal)
		{
			error_msg += "Error tolerance not reached (with Full Newton)!\n";
			error_msg += mystr("Error = ")+mystr(err)+mystr(", err-goal=") + mystr(error_goal);
			if (SolverPrintsDetailedOutput()) {SolverUO() << "  error message=Error tolerance not reached (with Full Newton)!" << "\n";}
		}
		else
		{
			error_msg += "Newton method not successful! ";
			if (SolverPrintsDetailedOutput()) {SolverUO() << "  error message=Newton method not successful! " << "\n";}
		}

		x0 = x0start;
		return 0;
	}

}

//$ PG 2013-11-6: assemble Jacobian (and print some information), set jaccount++
void NumNLSolver::AssembleJacobian(Matrix* minv, SparseMatrix* sminv, Vector& x0)
{
	if (nls->UseSparseJac())
	{
		nls->Jacobian(*sminv,x0);	
		jaccount++;
		if (GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiMatrix() && SolverPrintsDetailedOutput())
		{
			OutputJacobian(*sminv);
		}
		if (GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiCondition() && SolverPrintsDetailedOutput())
		{
			OutputConditionNumber(*sminv);
		}
	} 
	else
	{
		nls->Jacobian(*minv,x0);
		jaccount++;
		if (GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiMatrix() && SolverPrintsDetailedOutput())
		{
			OutputJacobian(*minv);
		}
		if (GetOptions()->LoggingOptions()->SolverNewtonIterationJacobiCondition() && SolverPrintsDetailedOutput())
		{
			OutputConditionNumber(*minv);
		}
	}
}
	
//$ PG 2013-10-17: calculate condition number of jacobi matrix, IS SLOW, should only be used for designing elements, constraints, for adjusting models, or for general debugging purpose
void NumNLSolver::OutputConditionNumber(const Matrix& M) const
{
	double rcond;
	int info = M.EstimateReciprocalConditionNumber(rcond);
	
	if (info < 0)
	{
		SolverUO() << "ERROR: Illegal value in matrix at position " << -info << " of the array!\n";
	}
	else if (info > 0)
	{
		SolverUO() << "ERROR: Zero division in LU factorization. M(" << info << "," << info << ") = 0, matrix is singular!\n";
	}
	else //if (info==0)
	{
		SolverUO() << "kappa = " << 1./rcond << "\n";
	}
}

//$ PG 2013-10-17: calculate condition number of jacobi matrix, IS SLOW, should only be used for designing elements, constraints, for adjusting models, or for general debugging purpose
void NumNLSolver::OutputConditionNumber(const SparseMatrix& M) const
{
	SolverUO() << "WARNING: Calculation of condition number for sparse matrices uses routine for full matrices!\n";
	OutputConditionNumber(M.GetMatrix());
}

//$ PG 2013-10-17: output of jacobi matrix, IS SLOW, should only be used for debugging reasons
void NumNLSolver::OutputJacobian(const Matrix& M) const
{
	ostringstream string;
	M.PrintToMatlab(string);    //extract string from stream via string.str().c_str()
	SolverUO() << "Jacobian (# " <<  jaccount << ") = " << string.str().c_str() << "\n";
}

//$ PG 2013-10-17: output of sparse jacobi matrix, IS SLOW, should only be used for debugging reasons
void NumNLSolver::OutputJacobian(const SparseMatrix& M) const
{
	SolverUO() << "WARNING: Output of sparse jacobi matrices uses routine for full matrices!\n";
	OutputJacobian(M.GetMatrix());
}




