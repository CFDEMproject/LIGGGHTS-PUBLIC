//#***************************************************************************************
//# filename:     femathhelperfunctions.cpp
//#
//# author:				Johannes Gerstmayr, Yuri Vetyukov
//# 
//# generated:      
//# description:  
//#                       
//# comments:      
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


#include "mbs_interface.h"
#include "myfile.h"
#include "femathhelperfunctions.h"

int DoesFileExist(mystr filename, UserOutputInterface * pUO, ofstream * logFile)
{
	CMFile file(filename, TFMread);
	//read parameters from file
	if (!file.IsGood()) 
	{
		if(pUO != NULL)
			pUO->InstantMessageText("WARNING: file " + filename + " does not exist \n");
		if(logFile != NULL)
			(*logFile) << "WARNING: " << filename << " does not exist \n";
		return 0; // file doesn't exist
	}
	else
	{
		return 1; // file exists
	}
}

void Quicksort(TArray<int>& x)
//Sorts an array x(1..x.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). x is replaced on output by its sorted rearrangement. Normally, the argument len should
//be set to the size of array x, but if len is smaller than this, then only the first len elements of a
//are sorted.
{
	int len = x.Length();
	int i,j,inc;
	int v;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= len);

	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc+1;i<=len;i++) 
		{ //Outer loop of straight insertion.
			v=x(i);
			j=i;
			while (x(j-inc) > v) 
			{ //Inner loop of straight insertion.
				x(j)=x(j-inc);
				j -= inc;
				if (j <= inc) break;
			}
			x(j)=v;
		}
	} while (inc > 1);
}

//template <typename T>
//void Quicksort(TArray<int>& a, TArray<T>& b)
////Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
////sort). a is replaced on output by its sorted rearrangement. b is sorted according to the order of a.
//{
//	int n = a.Length();
//	int i,j,inc;
//	int v;
//	double vb;
//	inc=1; //Determine the starting increment.
//	do 
//	{
//		inc *= 3;
//		inc++;
//	} while (inc <= n);
//	do 
//	{ //Loop over the partial sorts.
//		inc /= 3;
//		for (i=inc;i<n;i++) 
//		{ //Outer loop of straight insertion.
//			v=a.Get0(i);
//			vb=b.Get0(i);
//			j=i;
//			while (a.Get0(j-inc) > v) 
//			{ //Inner loop of straight insertion.
//				a.Elem0(j)=a.Get0(j-inc);
//				b.Elem0(j)=b.Get0(j-inc);
//				j -= inc;
//				if (j < inc) break;
//			}
//			a.Elem0(j)=v;
//			b.Elem0(j)=vb;
//		}
//	} while (inc > 1);
//
//}

void QuicksortDouble(TArray<double>& a, TArray<int>& b)
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b is sorted according to the order of a.
{
	int n = a.Length();
	int i,j,inc;
	double v;
	int vb;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= n);
	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc;i<n;i++) 
		{ //Outer loop of straight insertion.
			v=a.Get0(i);
			vb=b.Get0(i);
			j=i;
			while (a.Get0(j-inc) > v) 
			{ //Inner loop of straight insertion.
				a.Elem0(j)=a.Get0(j-inc);
				b.Elem0(j)=b.Get0(j-inc);
				j -= inc;
				if (j < inc) break;
			}
			a.Elem0(j)=v;
			b.Elem0(j)=vb;
		}
	} while (inc > 1);

}

void Quicksort(TArray<int>& a, TArray<int>& b, TArray<double>& c)
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b and c are sorted according to the order of a.
{
	int n = a.Length();
	int i,j,inc;
	int v;
	int vb;
	double vc;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= n);

	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc+1;i<=n;i++) 
		{ //Outer loop of straight insertion.
			v=a(i);
			vb=b(i);
			vc=c(i);
			j=i;
			while (a(j-inc) > v) 
			{ //Inner loop of straight insertion.
				a(j)=a(j-inc);
				b(j)=b(j-inc);
				c(j)=c(j-inc);
				j -= inc;
				if (j <= inc) break;
			}
			a(j)=v;
			b(j)=vb;
			c(j)=vc;
		}
	} while (inc > 1);
}

void QuicksortDouble(TArray<double>& a, TArray<double>& b, TArray<int>& c)
//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement. b and c are sorted according to the order of a.
{
	int n = a.Length();
	int i,j,inc;
	double v;
	double vb;
	int vc;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= n);

	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc+1;i<=n;i++) 
		{ //Outer loop of straight insertion.
			v=a(i);
			vb=b(i);
			vc=c(i);
			j=i;
			while (a(j-inc) > v) 
			{ //Inner loop of straight insertion.
				a(j)=a(j-inc);
				b(j)=b(j-inc);
				c(j)=c(j-inc);
				j -= inc;
				if (j <= inc) break;
			}
			a(j)=v;
			b(j)=vb;
			c(j)=vc;
		}
	} while (inc > 1);
}

// returns a list of redundant, optional ascending sort of the array (QSort sorts anyway...)
int FindRedundantEntries_QSort(TArray<int>& arr, TArray<int>& flags_redundant, int flag_sortascending)
{
	TArray<int> indices_sorted;
	indices_sorted.SetLen(arr.Length());
	for(int i=1; i <= arr.Length(); i++) indices_sorted(i) = i;

	flags_redundant.SetLen(arr.Length());
	flags_redundant.SetAll(0);

	if (flag_sortascending) // ascending sort of original array
	{
		Quicksort(arr,indices_sorted);
		for(int i=2; i <= arr.Length(); i++)
		{
			if(arr(i) == arr(i-1))
			{
				flags_redundant(i) = 1;
			}
		}
	}
	else // work on copy
	{
		TArray<int> copyofarr(arr);
		Quicksort(copyofarr,indices_sorted);
		for(int i=2; i <= arr.Length(); i++)
		{
			if(copyofarr(i) == copyofarr(i-1))
			{
				flags_redundant(indices_sorted(i)) = 1;
			}
		}
	}
	return 1; 
}

// removes redundant entries from a TArray<int>, optional ascending sort 
int RemoveRedundantEntries(TArray<int>& arr, int flag_sortascending)
{
	TArray<int> flags_redundant;    // 1 - redundant, 0 - unique && from !!Q-sorted!! sequence 

	FindRedundantEntries_QSort(arr,flags_redundant,flag_sortascending);
	arr.EraseMany(flags_redundant);

	return arr.Length();
}


// returns a list of redundant, optional ascending sort of the array (QSort sorts anyway...)
int FindRedundantEntries_QSort(TArray<double>& arr, TArray<int>& flags_redundant, int flag_sortascending)
{
	TArray<int> indices_sorted;
	indices_sorted.SetLen(arr.Length());
	for(int i=1; i <= arr.Length(); i++) indices_sorted(i) = i;

	flags_redundant.SetLen(arr.Length());
	flags_redundant.SetAll(0);

	if (flag_sortascending) // ascending sort of original array
	{
		QuicksortDouble(arr,indices_sorted);
		for(int i=2; i <= arr.Length(); i++)
		{
			if(arr(i) == arr(i-1))
			{
				flags_redundant(i) = 1;
			}
		}
	}
	else // work on copy
	{
		TArray<double> copyofarr(arr);
		QuicksortDouble(copyofarr,indices_sorted);
		for(int i=2; i <= arr.Length(); i++)
		{
			if(copyofarr(i) == copyofarr(i-1))
			{
				flags_redundant(indices_sorted(i)) = 1;
			}
		}
	}
	return 1; 
}

// removes redundant entries from a TArray<double>, optional ascending sort 
int RemoveRedundantEntries(TArray<double>& arr, int flag_sortascending)
{
	TArray<int> flags_redundant;    // 1 - redundant, 0 - unique && from !!Q-sorted!! sequence 

	FindRedundantEntries_QSort(arr,flags_redundant,flag_sortascending);
	arr.EraseMany(flags_redundant);

	return arr.Length();
}

// removes redundant entries (with error tolerance) from a TArray<double>, optional ascending sort 
int RemoveRedundantEntries_ErrorTolerance(TArray<double>& arr, double tol, int flag_sortascending)
{
	RemoveRedundantEntries(arr, flag_sortascending);

	for(int i=2; i<= arr.Length(); i++)
	{
		if( abs(arr(i-1) - arr(i)) < tol )
		{
			arr.Erase(i);
			i--;
		}
	}
	return arr.Length();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GetIntegrationRule(Vector& x, Vector& w, int order)
{
	if (order==0)
	{
		x.SetLen(0);
		w.SetLen(0);
		return;
	} 
	else if (order==1)
	{
		x.SetLen(1);
		w.SetLen(1);
		x(1) = 0;
		w(1) = 2.;
		return;
	} 
	else if (order<=3)
	{
		x.SetLen(2);
		w.SetLen(2);
		x(1) = -sqrt(1./3.);
		x(2) = sqrt(1./3.);
		w(1) = 1.;
		w(2) = 1.;
		return;
	} 
	else if (order<=5)
	{
		x.SetLen(3);
		w.SetLen(3);
		x(1) = -sqrt(3./5.);
		x(2) = 0.;
		x(3) = sqrt(3./5.);
		w(1) = 5./9.;
		w(2) = 8./9.;
		w(3) = 5./9.;
		return;
	}
	else if (order<=7)
	{
		x.SetLen(4);
		w.SetLen(4);
		x(1) = -sqrt(3./7.+sqrt(120.)/35.);
		x(2) = -sqrt(3./7.-sqrt(120.)/35.);
		x(3) = sqrt(3./7.-sqrt(120.)/35.);
		x(4) = sqrt(3./7.+sqrt(120.)/35.);
		w(1) = 1./2.-5./(3.*sqrt(120.));
		w(2) = 1./2.+5./(3.*sqrt(120.));
		w(3) = 1./2.+5./(3.*sqrt(120.));
		w(4) = 1./2.-5./(3.*sqrt(120.));		
		return;
	}
	else if (order<=9)
	{
		x.SetLen(5);
		w.SetLen(5);
		x(1) = -0.906179845938664;
		x(2) = -0.5384693101056831;
		x(3) = 0.;
		x(4) = 0.5384693101056831;
		x(5) = 0.906179845938664;
		w(1) = 0.23692688505618914;
		w(2) = 0.47862867049936636;
		w(3) = 0.5688888888888889;
		w(4) = 0.47862867049936636;
		w(5) = 0.23692688505618914;
		return;
	} 
	else if (order<=11)
	{
		x.SetLen(6);
		w.SetLen(6);
		x(1) = -0.932469514203152028;
		x(2) = -0.661209386466264514;
		x(3) = -0.238619186083196908;
		x(4) = 0.238619186083196910;
		x(5) = 0.661209386466264510;
		x(6) = 0.932469514203152030;
		w(1) = 0.171324492379170345;
		w(2) = 0.360761573048138608;
		w(3) = 0.467913934572691048;
		w(4) = 0.467913934572691048;
		w(5) = 0.360761573048138608;
		w(6) = 0.171324492379170345;
		return;
	}
	else if (order<=13)
	{
		x.SetLen(7);
		w.SetLen(7);
		x(1) = -0.949107912342758525;
		x(2) = -0.741531185599394440;
		x(3) = -0.405845151377397166;
		x(4) = 0.000000000000000000;
		x(5) = 0.405845151377397170;
		x(6) = 0.741531185599394440;
		x(7) = 0.949107912342758520;

		w(1) = 0.129484966168869693;
		w(2) = 0.279705391489276668;
		w(3) = 0.381830050505118944;
		w(4) = 0.417959183673469388;
		w(5) = 0.381830050505118944;
		w(6) = 0.279705391489276668;
		w(7) = 0.129484966168869693;

		return;
	}
	else if (order<=15)
	{
		x.SetLen(8);
		w.SetLen(8);
		x(1) = -0.960289856497536232;
		x(2) = -0.796666477413626740;
		x(3) = -0.525532409916328986;
		x(4) = -0.183434642495649804;
		x(5) = 0.183434642495649800;
		x(6) = 0.525532409916328990;
		x(7) = 0.796666477413626740;
		x(8) = 0.960289856497536230;

		w(1) = 0.101228536290376259;
		w(2) = 0.222381034453374470;
		w(3) = 0.313706645877887288;
		w(4) = 0.362683783378361982;
		w(5) = 0.362683783378361982;
		w(6) = 0.313706645877887288;
		w(7) = 0.222381034453374470;
		w(8) = 0.101228536290376259;
		return;
	}
	else if (order<=17)
	{
		x.SetLen(9);
		w.SetLen(9);

		x(1) = -0.9681602395076260898;
		x(2) = -0.8360311073266357943;
		x(3) = -0.6133714327005903973;
		x(4) = -0.3242534234038089290;
		x(5) = 0;
		x(6) = 0.32425342340380892904;
		x(7) = 0.61337143270059039731;
		x(8) = 0.83603110732663579430;
		x(9) = 0.96816023950762608984;

		w(1) = 0.0812743883615744120;
		w(2) = 0.1806481606948574041;
		w(3) = 0.2606106964029354623;
		w(4) = 0.3123470770400028401;
		w(5) = 0.3302393550012597632;
		w(6) = 0.3123470770400028401;
		w(7) = 0.2606106964029354623;
		w(8) = 0.1806481606948574041;
		w(9) = 0.081274388361574412 ;

		return;
	}
	else if (order<=19)
	{
		x.SetLen(10);
		w.SetLen(10);

		x(1) = -0.97390652851717174343;
		x(2) = -0.86506336668898464737;
		x(3) = -0.67940956829902443559;
		x(4) = -0.43339539412924721340;
		x(5) = -0.14887433898163116019;
		x(6) = 0.14887433898163116019;
		x(7) = 0.43339539412924721340;
		x(8) = 0.67940956829902443559;
		x(9) = 0.86506336668898464737;
		x(10) = 0.97390652851717174343;

		w(1) = 0.066671344308687693903;
		w(2) = 0.14945134915057972647;
		w(3) = 0.21908636251598201383;
		w(4) = 0.26926671930999629412;
		w(5) = 0.29552422471475303656;
		w(6) = 0.29552422471475303656;
		w(7) = 0.26926671930999629412;
		w(8) = 0.21908636251598201383;
		w(9) = 0.14945134915057972647 ;
		w(10) = 0.066671344308687693903;
		return;
	}

	else if (order <= 21)
	{
		x.SetLen(11);
		w.SetLen(11);

		x(1) = -0.97822865814605686197;
		x(2) = -0.88706259976809542778;
		x(3) = -0.73015200557404913440;
		x(4) = -0.51909612920681169612;
		x(5) = -0.26954315595234490388;
		x(6) = 0;
		x(7) = 0.26954315595234490388;
		x(8) = 0.51909612920681169612;
		x(9) = 0.73015200557404913440;
		x(10) = 0.88706259976809542778;
		x(11) = 0.97822865814605686197;

		w(1) = 0.055668567116175696197;
		w(2) = 0.12558036946490366836;
		w(3) = 0.18629021092773498380;
		w(4) = 0.23319376459199075979;
		w(5) = 0.26280454451024659601;
		w(6) = 0.27292508677790061622;
		w(7) = 0.26280454451024659601;
		w(8) = 0.23319376459199075979;
		w(9) = 0.18629021092773498380;
		w(10) = 0.12558036946490366836;
		w(11) = 0.055668567116175696197;
		return;
	}
	else if (order <= 23)
	{
		x.SetLen(12);
		w.SetLen(12);

		x(1) = -0.98156063424671913253;
		x(2) = -0.90411725637047490878;
		x(3) = -0.76990267419430469253;
		x(4) = -0.58731795428661737191;
		x(5) = -0.36783149899818012862;
		x(6) = -0.12523340851146880226;
		x(7) = 0.12523340851146880226;
		x(8) = 0.36783149899818012862;
		x(9) = 0.58731795428661737191;
		x(10) = 0.76990267419430469253;
		x(11) = 0.90411725637047490878;
		x(12) = 0.98156063424671913253;

		w(1) = 0.047175336386514117593;
		w(2) = 0.10693932599531794092;
		w(3) = 0.16007832854334633210;
		w(4) = 0.20316742672306598028;
		w(5) = 0.23349253653835491673;
		w(6) = 0.24914704581340282874;
		w(7) = 0.24914704581340282874;
		w(8) = 0.23349253653835491673;
		w(9) = 0.20316742672306598028; 
		w(10) = 0.16007832854334633210;
		w(11) = 0.10693932599531794092;
		w(12) = 0.047175336386514117593;
		return;
	}
	else if (order <= 25)
	{
		x.SetLen(13); w.SetLen(13);

		x(1) = -0.98418305471858813505;
		x(2) = -0.91759839922297792292;
		x(3) = -0.80157809073330987815;
		x(4) = -0.64234933944034011688;
		x(5) = -0.44849275103644670182;
		x(6) = -0.23045831595513477374;
		x(7) = 0;
		x(8) = 0.2304583159551347737;
		x(9) = 0.44849275103644670182;
		x(10) = 0.64234933944034011688;
		x(11) = 0.80157809073330987815;
		x(12) = 0.91759839922297792292;
		x(13) = 0.98418305471858813505;

		w(1) = 0.040484004765316182473;
		w(2) = 0.092121499837728437754;
		w(3) = 0.13887351021978752708;
		w(4) = 0.17814598076194601561;
		w(5) = 0.20781604753688867615;
		w(6) = 0.22628318026289734322;
		w(7) = 0.23255155323087389752;
		w(8) = 0.22628318026289734322;
		w(9) = 0.20781604753688867615;
		w(10) = 0.17814598076194601561;
		w(11) = 0.13887351021978752708;
		w(12) = 0.092121499837728437754;
		w(13) = 0.040484004765316182473;

		return;
	}
	else if (order <= 27)
	{
		x.SetLen(14); w.SetLen(14);
		w(1) = 0.035119460331750236570;
		x(2) = -0.92843488366357362906;
		w(2) = 0.080158087159759222606;
		x(3) = -0.82720131506976501967;
		w(3) = 0.12151857068790296312;
		x(4) = -0.68729290481168536786;
		w(4) = 0.15720316715819379616;
		x(5) = -0.51524863635815409957;
		w(5) = 0.18553839747793793302;
		x(6) = -0.31911236892788963360;
		w(6) = 0.20519846372129557643;
		x(7) = -0.10805494870734366764;
		w(7) = 0.21526385346315776714;
		x(8) = 0.10805494870734366764;
		w(8) = 0.21526385346315776714;
		x(9) = 0.31911236892788963360;
		w(9) = 0.20519846372129557643;
		x(10) = 0.51524863635815409957;
		w(10) = 0.18553839747793793302;
		x(11) = 0.68729290481168536786;
		w(11) = 0.15720316715819379616;
		x(12) = 0.82720131506976501967;
		w(12) = 0.12151857068790296312;
		x(13) = 0.92843488366357362906;
		w(13) = 0.080158087159759222606;
		x(14) = 0.98628380869681242515;
		w(14) = 0.035119460331750236570;
		return;
	}

	else if (order <= 29)
	{
		x.SetLen(15); w.SetLen(15);
		x(1) = -0.98799251802048537741;
		w(1) = 0.030753241996119028145;
		x(2) = -0.93727339240070595139;
		w(2) = 0.070366047488107749674;
		x(3) = -0.84820658341042731720;
		w(3) = 0.10715922046717152316;
		x(4) = -0.72441773136017006962;
		w(4) = 0.13957067792615432400;
		x(5) = -0.57097217260853883047;
		w(5) = 0.16626920581699403123;
		x(6) = -0.39415134707756349641;
		w(6) = 0.18616100001556210031;
		x(7) = -0.20119409399743459765;
		w(7) = 0.19843148532711168963;
		x(8) = 0;
		w(8) = 0.20257824192556128651;
		x(9) = 0.20119409399743459765;
		w(9) = 0.19843148532711168963;
		x(10) = 0.39415134707756349641;
		w(10) = 0.18616100001556210031;
		x(11) = 0.57097217260853883047;
		w(11) = 0.16626920581699403123;
		x(12) = 0.72441773136017006962;
		w(12) = 0.13957067792615432400;
		x(13) = 0.84820658341042731720;
		w(13) = 0.10715922046717152316;
		x(14) = 0.93727339240070595139;
		w(14) = 0.070366047488107749674;
		x(15) = 0.98799251802048537741;
		w(15) = 0.030753241996119028145;
		return;
	}
	assert(0);
	return;

}

void GetIntegrationRuleLobatto(Vector& x, Vector& w, int order)
{
	//maybe the order is actually one less, only for time integration the order is such high!
	if (order<=2)
	{
		x.SetLen(2);
		w.SetLen(2);
		x(1) = -1.000000000000000000;
		x(2) = 1.000000000000000000;
		w(1) = 1.000000000000000000;
		w(2) = 1.000000000000000000;
		return;
	} 
	else if (order<=4)
	{
		x.SetLen(3);
		w.SetLen(3);
		x(1) = -1.00000000000000000;
		x(2) = 0.00000000000000000;
		x(3) = 1.00000000000000000;
		w(1) = 0.33333333333333333;
		w(2) = 1.33333333333333333;
		w(3) = 0.33333333333333333;
		return;
	} 
	else if (order<=6)
	{
		x.SetLen(4);
		w.SetLen(4);
		x(1) = -1.00000000000000000;
		x(2) = -0.44721359549995794;
		x(3) = 0.44721359549995794;
		x(4) = 1.00000000000000000;
		w(1) = 0.16666666666666667;
		w(2) = 0.83333333333333333;
		w(3) = 0.83333333333333333;
		w(4) = 0.16666666666666667;

		return;
	} 
	else if (order<=8)
	{
		x.SetLen(5);
		w.SetLen(5);
		x(1) = -1.00000000000000000;
		x(2) = -0.65465367070797714;
		x(3) = 0.00000000000000000;
		x(4) = 0.65465367070797714;
		x(5) = 1.00000000000000000;
		w(1) = 0.10000000000000000;
		w(2) = 0.54444444444444444;
		w(3) = 0.71111111111111111;
		w(4) = 0.54444444444444444;
		w(5) = 0.10000000000000000;
		return;
	} 
	else if (order<=10)
	{
		x.SetLen(6);
		w.SetLen(6);
		x(1) = -1.00000000000000000;
		x(2) = -0.76505532392946469;
		x(3) = -0.28523151648064510;
		x(4) = 0.28523151648064510;
		x(5) = 0.76505532392946469;
		x(6) = 1.00000000000000000;
		w(1) = 0.06666666666666667;
		w(2) = 0.37847495629784698;
		w(3) = 0.55485837703548635;
		w(4) = 0.55485837703548635;
		w(5) = 0.37847495629784698;
		w(6) = 0.06666666666666667;
		return;
	} 
	else if (order<=12)
	{
		x.SetLen(7);
		w.SetLen(7);
		x(1) = -1.00000000000000000;
		x(2) = -0.83022389627856693;
		x(3) = -0.46884879347071421;
		x(4) = 0.00000000000000000;
		x(5) = 0.46884879347071421;
		x(6) = 0.83022389627856693;
		x(7) = 1.00000000000000000;
		w(1) = 0.04761904761904762;
		w(2) = 0.27682604736156595;
		w(3) = 0.43174538120986262;
		w(4) = 0.48761904761904762;
		w(5) = 0.43174538120986262;
		w(6) = 0.27682604736156595;
		w(7) = 0.04761904761904762;
		return;
	} 
	assert(0);
	return;
}

void GetIntegrationRuleTrig(Vector& x1, Vector& x2, Vector& w, int order)
{
	//int F(r,s) dr ds = sum(F(ri, si)*w(i)), factor 1/2 already included
	//according to Hughes, finite elements, 2000, p. 173
	if (order<=1) //?????
	{
		x1.SetLen(1);
		x2.SetLen(1);
		w.SetLen(1);
		x1(1) = 1./3.;
		x2(1) = 1./3.;
		w(1) = 0.5000000000000000000;
		return;
	} 
	else if (order<=2)
	{
		x1.SetLen(3);
		x2.SetLen(3);
		w.SetLen(3);
		x1(1) = 2./3.;
		x2(1) = 1./6.;
		x1(2) = 1./6.;
		x2(2) = 2./3.;
		x1(3) = 1./6.;
		x2(3) = 1./6.;
		w(1) = 1./6.;
		w(2) = 1./6.;
		w(3) = 1./6.;
		return;
	} 
	else if (order<=3)
	{
		x1.SetLen(4);
		x2.SetLen(4);
		w.SetLen(4);
		x1(1) = 1./3.;
		x2(1) = 1./3.;
		x1(2) = 0.6;
		x2(2) = 0.2;
		x1(3) = 0.2;
		x2(3) = 0.6;
		x1(4) = 0.2;
		x2(4) = 0.2;

		w(1) = -0.5*0.5625;
		w(2) =  0.5*0.52083333333333333;
		w(3) =  0.5*0.52083333333333333;
		w(4) =  0.5*0.52083333333333333;
		return;
	} 
	else if (order<=4)
	{
		x1.SetLen(6);
		x2.SetLen(6);
		w.SetLen(6);
		double a1 = 0.816847572980459;
		double a2 = 0.091576213509771;
		double b1 = 0.108103018168070;
		double b2 = 0.445948490915965;
		x1(1) = a1;
		x2(1) = a2;
		x1(2) = a2;
		x2(2) = a1;
		x1(3) = a2;
		x2(3) = a2;
		x1(4) = b1;
		x2(4) = b2;
		x1(5) = b2;
		x2(5) = b1;
		x1(6) = b2;
		x2(6) = b2;

		w(1) = 0.5*0.109951743655322;
		w(2) = 0.5*0.109951743655322;
		w(3) = 0.5*0.109951743655322;
		w(4) = 0.5*0.223381589678011;
		w(5) = 0.5*0.223381589678011;
		w(6) = 0.5*0.223381589678011;
		return;
	} 
	assert(0);
	return;
}







void GetIntegrationRuleTet(Vector& x1, Vector& x2, Vector& x3, Vector& w, int order)
{
	static double tet_order1_points1[] = 
	{ 0.25 };

	static double tet_order1_points2[] = 
	{ 0.25 };

	static double tet_order1_points3[] = 
	{ 0.25 };

	static double tet_order1_weights[] = 
	{
		1.0/6.0
	};

	static double tet_order2_points1[] = 
	{
		0.585410196624969,
			0.138196601125011,
			0.138196601125011,
			0.138196601125011
	};
	static double tet_order2_points2[] = 
	{
		0.138196601125011,
			0.585410196624969,
			0.138196601125011,
			0.138196601125011
	};
	static double tet_order2_points3[] = 
	{
		0.138196601125011,
			0.138196601125011,
			0.585410196624969,
			0.138196601125011
	};

	static double tet_order2_weights[] = 
	{ 1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0 };

	static double tet_order3_points1[] = 
	{
		0.25, 
			1./2.,
			1./6.,
			1./6.,
			1./6.
	};
	static double tet_order3_points2[] = 
	{
		0.25, 
			1./6.,
			1./2.,
			1./6.,
			1./6.
	};
	static double tet_order3_points3[] = 
	{
		0.25, 
			1./6.,
			1./6.,
			1./2.,
			1./6.
	};

	//$EK 2012-11-05 bug fix in quadrature rule; 5th factor 9./120. was missing
	//static double tet_order3_weights[] = 
	//{ -4./30., 9./120., 9./120., 9./120. }; changed to (according to Hughes' book; see below)
	static double tet_order3_weights[] = 
	{ -4./30., 9./120., 9./120., 9./120., 9./120.};

	static double tet_order5_points1[] = 
	{
		0.454496295874350,
			0.454496295874350,
			0.454496295874350,
			0.045503704125650,
			0.045503704125650,
			0.045503704125650,
			0.310885919263301,
			0.067342242210098,
			0.310885919263301,
			0.310885919263301,
			0.092735250310891,
			0.721794249067326,
			0.092735250310891,
			0.092735250310891
	};
	static double tet_order5_points2[] = 
	{
		0.454496295874350,
			0.045503704125650,
			0.045503704125650,
			0.454496295874350,
			0.454496295874350,
			0.045503704125650,
			0.310885919263301,
			0.310885919263301,
			0.067342242210098,
			0.310885919263301,
			0.092735250310891,
			0.092735250310891,
			0.721794249067326,
			0.092735250310891
	};
	static double tet_order5_points3[] = 
	{
		0.045503704125650,
			0.454496295874350,
			0.045503704125650,
			0.454496295874350,
			0.045503704125650,
			0.454496295874350,
			0.310885919263301,
			0.310885919263301,
			0.310885919263301,
			0.067342242210098,
			0.092735250310891,
			0.092735250310891,
			0.092735250310891,
			0.721794249067326
	};
	static double tet_order5_weights[] = 
	{
		0.007091003462847, 0.007091003462847, 0.007091003462847,
			0.007091003462847, 0.007091003462847, 0.007091003462847,
			0.018781320953003, 0.018781320953003, 0.018781320953003, 0.018781320953003,
			0.012248840519394, 0.012248840519394, 0.012248840519394, 0.012248840519394
	};


	//int F(r,s,t) dr ds dt = sum(F(ri, si, ti)*w(i)), factor 1/6 already included in weights
	//according to Hughes, finite elements, 2000, p. 174 and
	//Netgen, intrule.cc file
	if (order<=1)
	{
		x1.LinkWith(&tet_order1_points1[0],1);
		x2.LinkWith(&tet_order1_points2[0],1);
		x3.LinkWith(&tet_order1_points3[0],1);
		w.LinkWith( &tet_order1_weights[0],1);
		return;
	} 
	else if (order<=2)
	{
		x1.LinkWith(&tet_order2_points1[0],4);
		x2.LinkWith(&tet_order2_points2[0],4);
		x3.LinkWith(&tet_order2_points3[0],4);
		w.LinkWith( &tet_order2_weights[0],4);
		return;
	} 
	else if (order<=3)
	{
		x1.LinkWith(&tet_order3_points1[0],5);
		x2.LinkWith(&tet_order3_points2[0],5);
		x3.LinkWith(&tet_order3_points3[0],5);
		w.LinkWith( &tet_order3_weights[0],5);
		return;
	} 
	else if (order<=5)
	{
		x1.LinkWith(&tet_order5_points1[0],14);
		x2.LinkWith(&tet_order5_points2[0],14);
		x3.LinkWith(&tet_order5_points3[0],14);
		w.LinkWith( &tet_order5_weights[0],14);
		return;
	} 
	assert(0);
	return ;
}

double GetClockTime()
{
	timeb tb;
	tb.time = 0;
	tb.millitm = 0;
	ftime(&tb);
	return tb.time+(double)tb.millitm*0.001;
}